#!/usr/bin/env bash
#
# 04_ersilia_predictions.sh
#
# Fetch, serve and run each antimicrobial model over the EU OpenScreen and CoAdd
# compound libraries using the Ersilia CLI. This is the first live Ersilia-CLI
# step in the repo (all earlier predictions came from the Isaura precalc cache);
# it regenerates the EU OpenScreen predictions that were removed as stale.
#
# Models (16): the 15 pathogen models in config/pathogens_of_interest.csv
# (eosid column) plus the CoAdd model eos3dys (COADD_MODEL_ID in src/default.py).
# Each model predicts BOTH libraries -> 32 prediction runs.
#
# Ordering is LIBRARY-MAJOR: all EU OpenScreen predictions (every model) run FIRST,
# then all CoAdd predictions. EU OpenScreen has priority. Fetches are cached, so the
# CoAdd pass re-serves each model without re-downloading. Skip-if-exists leaves any
# already-computed {library}/{eosid}.csv untouched.
#
# Inputs (built here, skip-if-exists):
#   output/04_ersilia_predictions/inputs/euopenscreen_smiles.csv   (single `smiles` column)
#   output/04_ersilia_predictions/inputs/coadd_smiles.csv          (single `smiles` column)
# Outputs (one CSV per model, per library, bare eosid):
#   output/04_ersilia_predictions/euopenscreen/{eosid}.csv
#   output/04_ersilia_predictions/coadd/{eosid}.csv
#   output/04_ersilia_predictions/_failures.log                    (append-only)
#
# Requirements:
#   `ersilia` must be installed in a conda environment (default name: `ersilia`;
#   override with the ERSILIA_ENV variable). `ersilia` is intentionally NOT in
#   requirements.txt -- it is a separate CLI, not a pip dep of the analysis env.
#
# Usage:
#   bash scripts/04_ersilia_predictions.sh              # full run (overnight-scale)
#   SMOKE=1 bash scripts/04_ersilia_predictions.sh      # quick 1-model end-to-end check
#   ERSILIA_ENV=my_env bash scripts/04_ersilia_predictions.sh

set -uo pipefail   # NOT -e: a single model's failure must not abort the whole run

# --- Paths -------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

CONFIG_CSV="$REPO_ROOT/config/pathogens_of_interest.csv"
EUOS_MERGED_DIR="$REPO_ROOT/data/raw/euopenscreen_data/02_merged"
COADD_DIR="$REPO_ROOT/data/raw/coadd_data"

OUT_DIR="$REPO_ROOT/output/04_ersilia_predictions"
INPUTS_DIR="$OUT_DIR/inputs"
FAILURES_LOG="$OUT_DIR/_failures.log"

# --- Config ------------------------------------------------------------------
ERSILIA_ENV="${ERSILIA_ENV:-ersilia}"
SMOKE="${SMOKE:-0}"

mkdir -p "$INPUTS_DIR" "$OUT_DIR/euopenscreen" "$OUT_DIR/coadd"

log() { echo "[$(date '+%F %T')] $*"; }

# --- Activate the Ersilia conda env once -------------------------------------
# Activating once (rather than `conda run` per command) keeps the served-model
# session alive across the fetch/serve/run/close sequence.
# Fallback if activation is unavailable: prefix each ersilia call with
#   conda run --no-capture-output -n "$ERSILIA_ENV" ersilia ...
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
log "Using ersilia from: $(command -v ersilia)  (env: $ERSILIA_ENV)"

# --- Build consolidated SMILES inputs (skip-if-exists) -----------------------
euos_input="$INPUTS_DIR/euopenscreen_smiles.csv"
coadd_input="$INPUTS_DIR/coadd_smiles.csv"

# EU OpenScreen: already consolidated upstream to a single deduplicated `smiles`
# column. SMOKE mode uses the 1,500-row short.csv subset for a fast check.
if [ "$SMOKE" = "1" ]; then
    euos_source="$EUOS_MERGED_DIR/short.csv"
else
    euos_source="$EUOS_MERGED_DIR/02_only_smiles.csv"
fi
if [ -f "$euos_input" ]; then
    log "[skip] EU OpenScreen input exists ($euos_input)"
else
    [ -f "$euos_source" ] || { echo "ERROR: EU OpenScreen source not found: $euos_source" >&2; exit 1; }
    # already deduplicated upstream; preserve order, just drop any blank SMILES lines
    { head -1 "$euos_source"; tail -n +2 "$euos_source" | grep -v '^[[:space:]]*$'; } > "$euos_input"
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

# --- Model list --------------------------------------------------------------
[ -f "$CONFIG_CSV" ] || { echo "ERROR: config not found: $CONFIG_CSV" >&2; exit 1; }
# Read loop rather than `mapfile`: mapfile is a bash 4+ builtin and macOS ships
# bash 3.2, where it fails silently under `set -u` and leaves the list empty.
models=()
while IFS= read -r m || [ -n "$m" ]; do   # `|| [ -n ... ]` keeps a last line with no trailing newline
    [ -n "$m" ] && models+=("$m")
done < <(tail -n +2 "$CONFIG_CSV" | cut -d, -f3 | sed '/^$/d')

# Append the CoAdd model (read from src/default.py; fall back to eos3dys).
coadd_model="$(grep -oE 'COADD_MODEL_ID[[:space:]]*=[[:space:]]*"[^"]+"' "$REPO_ROOT/src/default.py" \
               | head -1 | grep -oE 'eos[0-9a-z]+' || true)"
coadd_model="${coadd_model:-eos3dys}"
models+=("$coadd_model")

# Optional MODELS override: a space-separated list of eosids that replaces the
# default (config order + CoAdd model). Lets you control ordering or run a subset
# without editing config/pathogens_of_interest.csv. Motivation: per-batch cost
# varies ~36x across these models (measured on the CoAdd library, batch size 100:
# eos5qya ~5 s, eos2e3s ~11 s, eos4an7 ~182 s -> ~50 h for one model), so a slow
# model early in the list blocks every faster one behind it.
#   MODELS="eos8lcw eos3dys" bash scripts/04_ersilia_predictions.sh
if [ -n "${MODELS:-}" ]; then
    models=()
    for m in $MODELS; do models+=("$m"); done
    log "MODELS override active (${#models[@]}): ${models[*]}"
fi

if [ "$SMOKE" = "1" ]; then
    models=("${models[0]}")
    log "SMOKE mode: single model (${models[0]}), EU OpenScreen subset."
fi
log "Models to run (${#models[@]}): ${models[*]}"

# --- One (library, model) prediction: skip / fetch / serve / run / close -----
run_one() {
    local eosid="$1" lib="$2"
    local out="$OUT_DIR/$lib/$eosid.csv"
    if [ -f "$out" ]; then
        log "  [skip] $lib / $eosid exists"
        return 0
    fi
    log "  fetch $eosid"; ersilia fetch "$eosid" || return 1
    log "  serve $eosid"; ersilia serve "$eosid" || { ersilia close >/dev/null 2>&1 || true; return 1; }
    log "  run   $lib / $eosid"
    # `ersilia run` is the current CLI; older builds use `ersilia api -i .. -o ..`.
    ersilia run -i "$INPUTS_DIR/${lib}_smiles.csv" -o "$out" || { ersilia close >/dev/null 2>&1 || true; return 1; }
    ersilia close >/dev/null 2>&1 || true
    return 0
}

# Library-major: all EU OpenScreen first (priority), then all CoAdd.
done_count=0
fail_count=0
for lib in euopenscreen coadd; do
    log "########## LIBRARY: $lib ##########"
    for eosid in "${models[@]}"; do
        log "=== $lib / $eosid ==="
        if run_one "$eosid" "$lib"; then
            done_count=$((done_count + 1))
        else
            log "[FAIL] $lib / $eosid"
            printf '%s\t%s\t%s\n' "$(date '+%F %T')" "$eosid" "$lib fetch/serve/run failed" >> "$FAILURES_LOG"
            ersilia close >/dev/null 2>&1 || true
            fail_count=$((fail_count + 1))
        fi
    done
done

log "Done. ${done_count} library-run(s) succeeded/skipped, ${fail_count} failed."
[ "$fail_count" -gt 0 ] && log "See failures in $FAILURES_LOG"
exit 0
