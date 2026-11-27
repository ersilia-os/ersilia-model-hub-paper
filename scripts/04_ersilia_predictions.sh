#!/usr/bin/env bash
#
# 04_ersilia_predictions.sh
#
# Fetch, serve and run each antimicrobial model over the EU OpenScreen and CoAdd
# compound libraries using the Ersilia CLI. This is the first live Ersilia-CLI
# step in the repo (all earlier predictions came from the Isaura precalc cache).
#
# Models: the 15 pathogen models in config/pathogens_of_interest.csv (eosid column)
# plus the CoAdd model eos3dys (COADD_MODEL_ID in src/default.py), minus EXCLUDE.
#
# Ordering is MODEL-MAJOR (changed 2026-11-09, was library-major): each model is
# fetched once, served once, run on BOTH libraries, then closed and DELETED (image,
# containers and EOS folder) before the next model is fetched. Only one model image
# (~5-11 GB) is on disk at any time, instead of all 16 (~133 GB).
#
# Versioning: each output file carries the model's major version in its name,
# {eosid}_{vMajor}.csv, the same convention as the Isaura files in
# data/processed/annotation_preds_ref_library. The version is the vX.Y.Z tag that Docker
# Hub's `latest` points to (checked before and after the pull) and must match the
# major version in the Airtable snapshot (src/airtable_utils.prediction_version), which
# is what steps 05-06 look up; a mismatch fails that model instead of writing a file
# the downstream steps would not find. The full release, the Docker image digest, the
# row count and the run date are appended to _versions.tsv.
#
# Skip-if-exists: a (model, library) pair is skipped when its versioned file exists;
# a model with both files present is not fetched at all. Predictions are written to
# _partial/ first and moved into place only when `ersilia run` succeeds, so an
# interrupted run never leaves a half-written file that would later be skipped.
#
# Inputs (built here, skip-if-exists):
#   output/04_ersilia_predictions/inputs/euopenscreen_smiles.csv   (single `smiles` column)
#   output/04_ersilia_predictions/inputs/coadd_smiles.csv          (single `smiles` column)
# Outputs:
#   output/04_ersilia_predictions/euopenscreen/{eosid}_{vMajor}.csv
#   output/04_ersilia_predictions/coadd/{eosid}_{vMajor}.csv
#   output/04_ersilia_predictions/_versions.tsv                    (append-only)
#   output/04_ersilia_predictions/_failures.log                    (append-only)
#
# Requirements:
#   `ersilia` installed in a conda environment (default name: `ersilia`; override with
#   ERSILIA_ENV) and a running Docker daemon. `ersilia` is intentionally NOT in
#   requirements.txt. The version helpers run on the ersilia env's python using only the
#   standard library (they read src/default.py and the Airtable snapshot directly).
#
# Usage (wrap in caffeinate so the Mac does not sleep during the multi-day run):
#   caffeinate -i bash scripts/04_ersilia_predictions.sh
#   SMOKE=1 bash scripts/04_ersilia_predictions.sh      # 1 model, ~1,500-compound subsets
#   MODELS="eos8lcw eos3dys" bash scripts/04_ersilia_predictions.sh
#   EXCLUDE="" bash scripts/04_ersilia_predictions.sh   # include eos43d6 again

set -uo pipefail   # NOT -e: a single model's failure must not abort the whole run

# --- Paths -------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

CONFIG_CSV="$REPO_ROOT/config/pathogens_of_interest.csv"
EUOS_MERGED_DIR="$REPO_ROOT/data/raw/euopenscreen_data/02_merged"
COADD_DIR="$REPO_ROOT/data/raw/coadd_data"

OUT_DIR="$REPO_ROOT/output/04_ersilia_predictions"
INPUTS_DIR="$OUT_DIR/inputs"
PARTIAL_DIR="$OUT_DIR/_partial"
FAILURES_LOG="$OUT_DIR/_failures.log"
VERSIONS_TSV="$OUT_DIR/_versions.tsv"
LIBRARIES=(euopenscreen coadd)   # EU OpenScreen first within each model

# --- Config ------------------------------------------------------------------
ERSILIA_ENV="${ERSILIA_ENV:-ersilia}"
SMOKE="${SMOKE:-0}"
# eos43d6 (M. tuberculosis) is "In maintenance" in Airtable; held back until it is checked
# by hand (user decision 2026-11-09). EXCLUDE="" puts it back.
EXCLUDE="${EXCLUDE-eos43d6}"

mkdir -p "$INPUTS_DIR" "$PARTIAL_DIR" "$OUT_DIR/euopenscreen" "$OUT_DIR/coadd"
[ -f "$VERSIONS_TSV" ] || printf 'date\teosid\tlibrary\trelease\thub_digest\tn_rows\tfile\n' > "$VERSIONS_TSV"

log() { echo "[$(date '+%F %T')] $*"; }

# --- Environment checks --------------------------------------------------------
if ! command -v conda >/dev/null 2>&1; then
    echo "ERROR: conda not found on PATH; cannot activate env '$ERSILIA_ENV'." >&2
    exit 1
fi
eval "$(conda shell.bash hook)"
conda activate "$ERSILIA_ENV" || { echo "ERROR: could not activate conda env '$ERSILIA_ENV'." >&2; exit 1; }
if ! command -v ersilia >/dev/null 2>&1; then
    echo "ERROR: 'ersilia' CLI not found in env '$ERSILIA_ENV'. Install it there or set ERSILIA_ENV." >&2
    exit 1
fi
if ! docker info >/dev/null 2>&1; then
    echo "ERROR: Docker daemon is not running. Start Docker Desktop first." >&2
    exit 1
fi
log "Using $(ersilia --version 2>/dev/null | tail -1) from $(command -v ersilia)  (env: $ERSILIA_ENV)"

# --- Version helpers -------------------------------------------------------------
# Expected major version from the Airtable snapshot. Mirrors
# src/airtable_utils.prediction_version (standard library only, so it runs in the ersilia env).
expected_version() {
    python - "$REPO_ROOT" "$1" <<'PY'
import csv, os, re, sys
root, model = sys.argv[1], sys.argv[2]
src = open(os.path.join(root, "src", "default.py")).read()
date = re.search(r'AIRTABLE_SNAPSHOT_DATE\s*=\s*"([^"]+)"', src).group(1)
path = os.path.join(root, "data", "raw", f"airtable_metadata_{date.replace('-', '')}.csv")
with open(path, encoding="utf-8-sig") as fh:
    for row in csv.DictReader(fh):
        if row["Identifier"] == model and row["Release"].strip():
            print(row["Release"].strip().split(".")[0])
            break
PY
}

# Release that Docker Hub's `latest` tag points to: the vX.Y.Z tag sharing its digest.
# Prints "<release> <digest>" (release empty if no version tag matches). The model's own
# packed card has no "Release" field — the image is built before the GitHub release adds it
# to metadata.yml, and the release step only re-tags that image (eos21dr v5.0.0, 2026-11-09).
hub_release() {
    python - "$1" <<'PY'
import json, re, sys, urllib.request
url = f"https://hub.docker.com/v2/repositories/ersiliaos/{sys.argv[1]}/tags?page_size=100"
tags = []
try:
    while url:
        with urllib.request.urlopen(url, timeout=30) as fh:
            page = json.load(fh)
        tags += page.get("results", [])
        url = page.get("next")
except OSError:
    print(" ")
    sys.exit(0)
digest = next((t.get("digest") for t in tags if t["name"] == "latest"), None)
versions = [t["name"] for t in tags
            if digest and t.get("digest") == digest and re.fullmatch(r"v\d+\.\d+\.\d+", t["name"])]
versions.sort(key=lambda v: tuple(int(x) for x in v[1:].split(".")))
print(versions[-1] if versions else "", digest or "")
PY
}

# --- Build consolidated SMILES inputs (skip-if-exists) -----------------------
euos_input="$INPUTS_DIR/euopenscreen_smiles.csv"
coadd_input="$INPUTS_DIR/coadd_smiles.csv"
if [ "$SMOKE" = "1" ]; then
    # Smoke inputs live under their own names so they never shadow the full libraries.
    euos_input="$INPUTS_DIR/smoke_euopenscreen_smiles.csv"
    coadd_input="$INPUTS_DIR/smoke_coadd_smiles.csv"
fi

# EU OpenScreen: already consolidated upstream to a single deduplicated `smiles`
# column. SMOKE mode takes its first 1,500 compounds (the former short.csv subset is no
# longer staged).
euos_source="$EUOS_MERGED_DIR/02_only_smiles.csv"
if [ -f "$euos_input" ]; then
    log "[skip] EU OpenScreen input exists ($euos_input)"
else
    [ -f "$euos_source" ] || { echo "ERROR: EU OpenScreen source not found: $euos_source" >&2; exit 1; }
    # already deduplicated upstream; preserve order, just drop any blank SMILES lines
    n_keep=999999999; [ "$SMOKE" = "1" ] && n_keep=1500
    { head -1 "$euos_source"; tail -n +2 "$euos_source" | grep -v '^[[:space:]]*$' | head -n "$n_keep"; } > "$euos_input"
    log "EU OpenScreen input -> $euos_input ($(( $(wc -l < "$euos_input") - 1 )) compounds)"
fi

# CoAdd: the full screening library from the canonical 00_smiles_info.csv (staged by
# 00_download_data.py). Predict on the standardized SMILES (`std_smiles`), matching
# ../new-modelling; dedup on the exact std_smiles string. Rows with an EMPTY std_smiles
# (199 exotic structures that failed upstream standardization) are dropped here so no
# blank SMILES is fed to Ersilia -> ~100,005 compounds. SMILES contain no commas, so the
# column is safe to `cut`.
coadd_info="$COADD_DIR/00_smiles_info.csv"
if [ -f "$coadd_input" ]; then
    log "[skip] CoAdd input exists ($coadd_input)"
else
    [ -f "$coadd_info" ] || { echo "ERROR: CoAdd SMILES list not found: $coadd_info (run 00_download_data.py)" >&2; exit 1; }
    std_col="$(head -1 "$coadd_info" | tr ',' '\n' | grep -nx 'std_smiles' | cut -d: -f1)"
    [ -n "$std_col" ] || { echo "ERROR: no 'std_smiles' column in $coadd_info" >&2; exit 1; }
    if [ "$SMOKE" = "1" ]; then
        { echo "smiles"; tail -n +2 "$coadd_info" | cut -d, -f"$std_col" | grep -v '^[[:space:]]*$' | sort -u | head -n 1500; } > "$coadd_input"
    else
        { echo "smiles"; tail -n +2 "$coadd_info" | cut -d, -f"$std_col" | grep -v '^[[:space:]]*$' | sort -u; } > "$coadd_input"
    fi
    log "CoAdd input -> $coadd_input ($(( $(wc -l < "$coadd_input") - 1 )) unique compounds)"
fi

input_for() { [ "$1" = "euopenscreen" ] && echo "$euos_input" || echo "$coadd_input"; }

# --- Model list --------------------------------------------------------------
[ -f "$CONFIG_CSV" ] || { echo "ERROR: config not found: $CONFIG_CSV" >&2; exit 1; }
# Read loop rather than `mapfile`: mapfile is a bash 4+ builtin and macOS ships
# bash 3.2, where it fails silently under `set -u` and leaves the list empty.
models=()
while IFS= read -r m || [ -n "$m" ]; do   # `|| [ -n ... ]` keeps a last line with no trailing newline
    [ -n "$m" ] && models+=("$m")
done < <(tail -n +2 "$CONFIG_CSV" | cut -d, -f3 | sed '/^$/d')

coadd_model="$(grep -oE 'COADD_MODEL_ID[[:space:]]*=[[:space:]]*"[^"]+"' "$REPO_ROOT/src/default.py" \
               | head -1 | grep -oE 'eos[0-9a-z]+' || true)"
models+=("${coadd_model:-eos3dys}")

# MODELS replaces the default list (order or subset) without editing the config.
if [ -n "${MODELS:-}" ]; then
    models=()
    for m in $MODELS; do models+=("$m"); done
    log "MODELS override active (${#models[@]}): ${models[*]}"
fi

if [ -n "$EXCLUDE" ]; then
    kept=()
    for m in "${models[@]}"; do
        case " $EXCLUDE " in *" $m "*) log "Excluded: $m";; *) kept+=("$m");; esac
    done
    models=("${kept[@]}")
fi

if [ "$SMOKE" = "1" ]; then
    models=("${models[0]}")
    log "SMOKE mode: single model (${models[0]}), ~1,500-compound subsets."
fi
log "Models to run (${#models[@]}): ${models[*]}"

# --- One model: fetch -> serve -> run every library -> close -> delete ----------
cleanup_model() {
    ersilia close >/dev/null 2>&1 || true
    log "  delete $1"
    ersilia delete "$1" >/dev/null 2>&1 || log "  [warn] could not delete $1; remove it by hand"
}

fail() {   # eosid, message
    log "[FAIL] $1: $2"
    printf '%s\t%s\t%s\n' "$(date '+%F %T')" "$1" "$2" >> "$FAILURES_LOG"
}

run_model() {
    local eosid="$1" lib out partial todo=() vexp release vmaj digest n
    vexp="$(expected_version "$eosid")"
    [ -n "$vexp" ] || { fail "$eosid" "no Release in the Airtable snapshot"; return 1; }

    for lib in "${LIBRARIES[@]}"; do
        out="$OUT_DIR/$lib/${eosid}_${vexp}.csv"
        if [ -f "$out" ] && [ "$SMOKE" != "1" ]; then
            log "  [skip] $lib/${eosid}_${vexp}.csv exists"
        else
            todo+=("$lib")
        fi
    done
    [ "${#todo[@]}" -gt 0 ] || return 0

    local before after
    before="$(hub_release "$eosid")"
    log "  fetch $eosid (expecting $vexp; Docker Hub latest = ${before% *})"
    ersilia fetch "$eosid" || { fail "$eosid" "fetch failed"; cleanup_model "$eosid"; return 1; }
    after="$(hub_release "$eosid")"
    release="${after% *}"
    digest="${after#* }"
    vmaj="${release%%.*}"
    log "  fetched release ${release:-unknown} (Docker Hub digest ${digest:-unknown})"
    # `latest` must not move during the pull, or the recorded release may not be what was pulled.
    if [ "$before" != "$after" ]; then
        fail "$eosid" "Docker Hub latest changed during the fetch ($before -> $after)"
        cleanup_model "$eosid"
        return 1
    fi
    if [ -z "$release" ] || [ "$vmaj" != "$vexp" ]; then
        fail "$eosid" "Docker Hub release '${release:-unknown}' does not match Airtable $vexp"
        cleanup_model "$eosid"
        return 1
    fi

    log "  serve $eosid"
    ersilia serve "$eosid" || { fail "$eosid" "serve failed"; cleanup_model "$eosid"; return 1; }

    local status=0
    for lib in "${todo[@]}"; do
        out="$OUT_DIR/$lib/${eosid}_${vmaj}.csv"
        [ "$SMOKE" = "1" ] && out="$PARTIAL_DIR/smoke_${lib}_${eosid}_${vmaj}.csv"
        partial="$PARTIAL_DIR/${lib}_${eosid}.csv"
        rm -f "$partial"
        log "  run   $lib / $eosid"
        if ersilia run -i "$(input_for "$lib")" -o "$partial" && [ -s "$partial" ]; then
            mv "$partial" "$out"
            n=$(( $(wc -l < "$out") - 1 ))
            [ "$SMOKE" = "1" ] || printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$(date '+%F %T')" "$eosid" "$lib" \
                "$release" "${digest:-unknown}" "$n" "$(basename "$out")" >> "$VERSIONS_TSV"
            log "  -> $out ($n rows)"
        else
            fail "$eosid" "$lib run failed"
            status=1
        fi
    done

    cleanup_model "$eosid"
    return $status
}

done_count=0
fail_count=0
for eosid in "${models[@]}"; do
    log "=== $eosid ==="
    if run_model "$eosid"; then
        done_count=$((done_count + 1))
    else
        fail_count=$((fail_count + 1))
    fi
done

log "Done. ${done_count} model(s) complete or skipped, ${fail_count} with failures."
[ "$fail_count" -gt 0 ] && log "See $FAILURES_LOG"
exit 0
