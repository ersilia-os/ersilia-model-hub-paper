"""Prediction-quality scoring of the ChEMBL/PubChem pathogen models with `eosquality`.

`eosquality` (https://github.com/ersilia-os/eosquality) quantifies whether a given Ersilia
run output is "trustworthy" — explicitly NOT the probability that a prediction is correct.
It works in two stages:

  fit  — once per Ersilia model, against that model's predictions on the ENTIRE canonical
         Ersilia reference library (1,355,109 molecules). Learns the distributional and
         neighbourhood shape of the model's own output space. The reference population is
         the model's own predictions, NOT ground truth.
  run  — scores an arbitrary query CSV of the SAME model's predictions against that fit.

Four scores are used here, each calibrated to (0, 1] against the reference's own
distribution (so the reference median is ~0.5 BY CONSTRUCTION), each with a pre-calibration
`*_raw` companion:

  typicality   per-column: how common is this predicted value in the reference?
  extremity    per-column: how far from the column centre does it sit?
  support      fingerprint space: how close is this molecule to the reference set?
  consistency  output space: how quiet is its fingerprint-neighbourhood?

The fifth score, `signal` (XGBoost/SHAP Gini attribution), is DELIBERATELY EXCLUDED
(user-directed, 2026-09-15). It is opt-in upstream, self-labelled provisional
(SIGNAL_FORMULA_VERSION = "gini_v2"), is not exported from eosquality's __init__.py, and
defaults to only 1,000 XGBoost training rows.

Two upstream bugs in eosquality 0.1.0 are worked around here, neither by patching it:

1. BROKEN DOWNLOAD URL. eosquality/library/identity.py:52 hardcodes the S3 prefix
   ".../eosvc-public/eosquality/indices/", but eosvc actually publishes under
   ".../eosvc-public/eosquality/data/indices/" — note the missing `data/` segment. Every
   file 403s at the baked-in URL. EOSQUALITY_BASE_URL below supplies the corrected prefix
   via the EOSQUALITY_REFERENCE_BASE_URL env var, which fixes BOTH the index folder and
   the sibling library CSV (eosquality's library_csv_url() derives the CSV path by swapping
   the trailing "indices/" for "libraries/").

2. STRICT RDKIT VERSION GATE. VectorIndex._check_rdkit_version (vectorindex.py:363-377)
   compares the running RDKit against the index's build version with a strict string `!=`
   and raises RuntimeError on ANY difference. The published index was built with RDKit
   2026.03.1. check_rdkit_version() below fails early with the exact remediation command
   rather than letting `fit` die deep in a subprocess.

WHY THE FULL 1.36M LIBRARY, with no subsampling: this is forced by eosquality, not chosen.
quality.py:204-209 requires len(reference) == len(index) exactly, and
_validate_reference_against_library_csv (quality.py:914-938) requires reference["input"] to
match the canonical library SMILES row-for-row. The Isaura reference-library predictions
already staged in this repo satisfy both as-is — verified row-aligned against
data/raw/compound_lists/reference_library_smiles.csv — so they are passed to `fit`
DIRECTLY, with no copying or reformatting. Note also that eosquality's fit input filename
is load-bearing (cli/fit.py:45 parses `eos<id>_v<n>` out of it); `eos5eya_v3.csv` already
parses correctly.

NOTE ON DATA DOWNLOAD (CLAUDE.md convention): the ~2.2 GB eosquality reference library is
fetched into the ~/.eosquality/ PACKAGE CACHE, not into this repo's data/ tree — the same
status as the ersilia models under ~/eos, which are likewise not routed through
00_download_data.py. It is fetched here (user-directed, 2026-09-15) rather than in
00_download_data.py. A NOTE: comment in Section 1 of 00_download_data.py points here.

Run this script in the `paper` conda env, which includes eosquality (it needs the
`eosquality` CLI on PATH). Once the merged summary covers all 15 pathogens it also draws the
per-pathogen panel (src/plots_eosquality.py); `--plot-only` redraws it without re-scoring.

Requires:
    data/processed/annotation_preds_ref_library/{eosid}_v*.csv
        — Isaura reference-library predictions for the pathogen model
          (staged by 00_download_data.py, Section 4).
    output/11_reference_library_projection/11_top1000_per_pathogen.csv
        — produced by 11_reference_library_projection.py; supplies the top-1000 query set.
    config/pathogens_of_interest.csv
    ~/.eosquality/  — fetched automatically by this script on first run.

Usage:
    conda activate paper
    python xx_eosquality.py                       # E. coli / eos5eya (default)
    python xx_eosquality.py --pathogen saureus
    python xx_eosquality.py --force               # refit / rescore (asks before deleting)
    python xx_eosquality.py --plot-only           # redraw the panel from the summary CSV

Outputs:
    output/xx_eosquality/artifacts_{eosid}_{version}/          fitted reference population
    output/xx_eosquality/queries/{queryset}_{eosid}_{version}.csv
    output/xx_eosquality/xx_scores_{queryset}_{eosid}_{version}.csv
    output/xx_eosquality/xx_summary_{eosid}_{version}.csv
"""

import argparse
import glob
import json
import os
import shutil
import subprocess
import sys
import time

import pandas as pd
import rdkit

root = os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.abspath(os.path.join(root, ".."))
sys.path.append(os.path.join(root, "..", "src"))
from chembl_models_analyses_common import merge_into_combined_summary  # noqa: E402
from default import RANDOM_SEED  # noqa: E402
from eval_abx_enrichment import pathogen_phylo_order  # noqa: E402
from plots_eosquality import save_eosquality_figure  # noqa: E402

REFLIB_DIR = os.path.join(repo_root, "data", "processed", "annotation_preds_ref_library")
PATHOGENS_CONFIG_PATH = os.path.join(repo_root, "config", "pathogens_of_interest.csv")
# Same two config inputs step 12 passes pathogen_phylo_order, so the panel's row order
# matches the step 10 / 12 matrices.
ENDPOINT_SELECTION_PATH = os.path.join(repo_root, "config", "08_endpoint_selection.csv")
TAXONOMY_PATH = os.path.join(repo_root, "config", "organism_taxonomy.csv")
TOP1000_PATH = os.path.join(
    repo_root, "output", "11_reference_library_projection", "11_top1000_per_pathogen.csv"
)

output_dir = os.path.join(repo_root, "output", "xx_eosquality")
QUERIES_DIR = os.path.join(output_dir, "queries")
# One row per (pathogen, query set, score), merged across runs so the 15 pathogens can be
# run one at a time without earlier results disappearing (see merge_combined_summary).
COMBINED_SUMMARY_PATH = os.path.join(output_dir, "xx_all_pathogens_summary.csv")
os.makedirs(output_dir, exist_ok=True)
os.makedirs(QUERIES_DIR, exist_ok=True)

# Corrected S3 prefix — see bug 1 in the module docstring. eosquality's own
# DEFAULT_REFERENCE_BASE_URL omits the `data/` segment and 403s on every file.
EOSQUALITY_BASE_URL = (
    "https://eosvc-public.s3.amazonaws.com/eosquality/data/indices/"
)

# The four calibrated scores to fit. `signal` is excluded on purpose — see the module
# docstring. eosquality's own DEFAULT_SCORES happens to be this same set.
SCORES = "typicality,support,consistency,extremity"

# RANDOM_QUERY_SIZE: size of the null-control query set. NOTE this is an IN-SAMPLE random
# subsample, not a held-out set — eosquality forces the fit to cover the entire library
# (see the module docstring), so nothing is left out to hold out. It still serves its
# purpose: because every calibrated score is a CDF lookup against the reference's own
# distribution, a random draw from that same reference MUST centre near 0.5. That is the
# check that the pipeline is behaving, before any number from the top-1000 set is read.
# 10,000 is comfortably enough to resolve a median to ~0.005 and is the same order as
# eosquality's own MIN_REFERENCE_SAMPLES.
RANDOM_QUERY_SIZE = 10_000

pathogens_df = pd.read_csv(PATHOGENS_CONFIG_PATH)
PATHOGEN_TO_EOSID = dict(zip(pathogens_df["code"], pathogens_df["eosid"]))
ALL_PATHOGENS = list(PATHOGEN_TO_EOSID.keys())

parser = argparse.ArgumentParser(
    description="Score an Ersilia pathogen model's predictions with eosquality. "
                "See the top of this file for the full methodology."
)
parser.add_argument(
    "--pathogen", default="ecoli", choices=ALL_PATHOGENS,
    help="Pathogen code from config/pathogens_of_interest.csv (default: ecoli).",
)
parser.add_argument(
    "--force", action="store_true",
    help="Recompute the fit and the scores even if they already exist. Prompts for "
         "confirmation before removing anything.",
)
parser.add_argument(
    "--plot-only", action="store_true",
    help="Skip fitting and scoring; only redraw the per-pathogen panel from the merged "
         "summary CSV.",
)
args = parser.parse_args()


def sh(cmd, env_extra=None, stream=True):
    """Run *cmd* (a list), streaming its output. Returns the CompletedProcess."""
    env = dict(os.environ)
    if env_extra:
        env.update(env_extra)
    print(f"  $ {' '.join(cmd)}", flush=True)
    return subprocess.run(cmd, env=env, check=False, text=True,
                          capture_output=not stream)


def confirm_removal(path):
    """Ask before deleting anything (CLAUDE.md: never delete without confirmation).

    Applies to files as well as folders — a superseded score CSV is still a result, and
    old analysis files may have scientific value.
    """
    if not sys.stdin.isatty():
        sys.exit(
            f"ERROR: --force wants to delete {path}\n"
            "       Refusing to do that non-interactively. Re-run in a terminal, or "
            "remove it yourself."
        )
    reply = input(f"  --force will DELETE {path}\n  Type 'yes' to confirm: ").strip()
    if reply != "yes":
        sys.exit("  Aborted; nothing was deleted.")
    shutil.rmtree(path) if os.path.isdir(path) else os.remove(path)
    print(f"  removed {path}")


def ensure_library():
    """Fetch the eosquality reference library into ~/.eosquality/ if not already valid.

    `eosquality download` is a no-op (and touches no network) when a complete, valid
    cached copy is already present, so this is safe to call on every run.
    """
    print("\n[1/5] eosquality reference library")
    t0 = time.time()
    res = sh(["eosquality", "download"],
             env_extra={"EOSQUALITY_REFERENCE_BASE_URL": EOSQUALITY_BASE_URL})
    if res.returncode != 0:
        sys.exit(
            f"ERROR: `eosquality download` failed (exit {res.returncode}).\n"
            f"       Base URL used: {EOSQUALITY_BASE_URL}\n"
            "       If this 403s, the published S3 prefix has moved again — see bug 1 "
            "in this file's docstring."
        )
    print(f"  done in {time.time() - t0:.0f}s")


def check_rdkit_version():
    """Fail early on eosquality's strict RDKit gate (bug 2 in the module docstring)."""
    from eosquality.library.identity import reference_library_path

    meta_path = os.path.join(reference_library_path(), "metadata.json")
    with open(meta_path) as f:
        stored = json.load(f).get("rdkit_version")
    current = rdkit.__version__
    print(f"  index built with RDKit {stored}; environment has RDKit {current}")
    if stored is not None and stored != current:
        sys.exit(
            f"ERROR: RDKit version mismatch — eosquality compares these with a strict "
            f"string `!=` (vectorindex.py:363-377) and will raise RuntimeError.\n"
            f"       Fix, inside the `paper` conda env:\n"
            f"           pip install rdkit=={stored.replace('.0', '.')}\n"
            f"       (the published index is pinned to RDKit {stored})"
        )


def resolve_reference_csv(eosid):
    """Locate the Isaura reference-library prediction CSV for *eosid*, newest version."""
    matches = sorted(glob.glob(os.path.join(REFLIB_DIR, f"{eosid}_v*.csv")))
    if not matches:
        sys.exit(
            f"ERROR: no reference-library predictions found at "
            f"{os.path.join(REFLIB_DIR, eosid + '_v*.csv')}\n"
            "       Run 00_download_data.py (Section 4 — Isaura) first."
        )
    path = matches[-1]
    version = os.path.basename(path).rsplit("_", 1)[1].split(".")[0]
    return path, version


def fit_reference(reference_csv, artifacts_dir):
    """Fit the reference population. Skipped if artifacts already exist."""
    print("\n[2/5] fit reference population")
    if os.path.exists(artifacts_dir):
        if args.force:
            confirm_removal(artifacts_dir)
        else:
            print(f"  artifacts already exist, skipping fit -> {artifacts_dir}")
            print("  (pass --force to refit)")
            return
    t0 = time.time()
    res = sh(["eosquality", "fit", "-i", reference_csv, "-o", artifacts_dir,
              "--scores", SCORES, "--verbose"],
             env_extra={"EOSQUALITY_REFERENCE_BASE_URL": EOSQUALITY_BASE_URL})
    if res.returncode != 0:
        sys.exit(f"ERROR: `eosquality fit` failed (exit {res.returncode}).")
    print(f"  fit done in {time.time() - t0:.0f}s")


def report_feature_selection(artifacts_dir):
    """Print which prediction columns survived eosquality's correlation-cluster reduction.

    eosquality's --max-features default of 10 is kept as-is (user-directed, 2026-09-15),
    and is the #1 open TODO in eosquality's own README. Any model with more than 10 numeric
    columns therefore has some silently dropped, so we name them explicitly here.
    """
    shared = os.path.join(artifacts_dir, "shared")
    with open(os.path.join(shared, "schema.json")) as f:
        all_cols = [c["name"] for c in json.load(f)["columns"]]
    with open(os.path.join(shared, "selected_columns.json")) as f:
        kept = json.load(f)["selected_columns"]
    dropped = [c for c in all_cols if c not in kept]
    print(f"  feature selection: {len(all_cols)} columns -> {len(kept)} kept")
    print(f"    kept    : {', '.join(kept)}")
    if dropped:
        print(f"    DROPPED : {', '.join(dropped)}")
        print("    (eosquality --max-features default = 10, correlation-cluster medoids)")
    return all_cols, kept, dropped


def build_queries(reference_csv, eosid, version, pathogen):
    """Build the two query CSVs. Reads the reference prediction CSV exactly once."""
    print("\n[3/5] build query sets")
    t0 = time.time()
    print(f"  reading {os.path.basename(reference_csv)} ...")
    df = pd.read_csv(reference_csv)
    print(f"  {len(df):,} rows x {len(df.columns)} columns in {time.time() - t0:.0f}s")

    n_nan = int(df.isna().sum().sum())
    if n_nan:
        print(f"  NOTE: {n_nan:,} NaN cells present in the reference predictions. "
              "Nothing is dropped; eosquality handles NaNs per-column.")

    queries = {}

    random_path = os.path.join(QUERIES_DIR, f"random{RANDOM_QUERY_SIZE}_{eosid}_{version}.csv")
    if os.path.exists(random_path) and not args.force:
        print(f"  exists, reusing -> {os.path.basename(random_path)}")
    else:
        sample = df.sample(n=RANDOM_QUERY_SIZE, random_state=RANDOM_SEED)
        sample.to_csv(random_path, index=False)
        print(f"  random control: {len(sample):,} rows (seed {RANDOM_SEED}) "
              f"-> {os.path.basename(random_path)}")
    queries[f"random{RANDOM_QUERY_SIZE}"] = random_path

    top_path = os.path.join(QUERIES_DIR, f"top1000_{eosid}_{version}.csv")
    if os.path.exists(top_path) and not args.force:
        print(f"  exists, reusing -> {os.path.basename(top_path)}")
    else:
        top = pd.read_csv(TOP1000_PATH)
        keys = top.loc[top["pathogen_code"] == pathogen, "key"]
        merged = df[df["key"].isin(set(keys))]
        print(f"  top-1000 hits: {len(keys):,} keys for '{pathogen}' -> "
              f"{len(merged):,} matched in the prediction file")
        if len(merged) != len(keys):
            print(f"    NOTE: {len(keys) - len(merged):,} key(s) did not match. "
                  "Nothing dropped silently — this count is the discrepancy.")
        merged.to_csv(top_path, index=False)
    queries["top1000"] = top_path

    return queries


def score_queries(queries, artifacts_dir, eosid, version):
    """Run each query set against the fitted artifacts."""
    print("\n[4/5] score query sets")
    outputs = {}
    for name, query_path in queries.items():
        out_path = os.path.join(output_dir, f"xx_scores_{name}_{eosid}_{version}.csv")
        outputs[name] = out_path
        if os.path.exists(out_path):
            if args.force:
                confirm_removal(out_path)
            else:
                print(f"  scores exist, skipping -> {os.path.basename(out_path)}")
                continue
        t0 = time.time()
        res = sh(["eosquality", "run", "-i", query_path, "-a", artifacts_dir,
                  "-o", out_path, "--verbose"],
                 env_extra={"EOSQUALITY_REFERENCE_BASE_URL": EOSQUALITY_BASE_URL})
        if res.returncode != 0:
            sys.exit(f"ERROR: `eosquality run` failed on '{name}' (exit {res.returncode}).")
        print(f"  {name} scored in {time.time() - t0:.0f}s -> {os.path.basename(out_path)}")
    return outputs


def merge_combined_summary(new_rows):
    """Merge *new_rows* into the cross-pathogen table, replacing only its own pathogen.

    Delegates the merge to src/chembl_models_analyses_common.py's
    merge_into_combined_summary, then re-sorts so the table order is stable across runs.
    """
    combined = merge_into_combined_summary(new_rows, COMBINED_SUMMARY_PATH)
    combined = combined.sort_values(["pathogen", "query_set", "score"])
    combined.to_csv(COMBINED_SUMMARY_PATH, index=False)
    return combined


def summarise(outputs, eosid, version, pathogen, artifacts_dir):
    """Write per-query-set descriptive statistics. Deliberately draws no conclusions."""
    print("\n[5/5] summary")
    with open(os.path.join(artifacts_dir, "manifest.json")) as f:
        manifest = json.load(f)
    rows = []
    for name, path in outputs.items():
        df = pd.read_csv(path)
        score_cols = [c for c in df.columns if c not in ("key", "input")]
        for col in score_cols:
            s = df[col]
            rows.append({
                "pathogen": pathogen, "eosid": eosid, "version": version,
                "n_features": manifest["n_features"],
                "n_features_selected": manifest["n_features_selected"],
                "query_set": name, "score": col, "n": int(s.notna().sum()),
                "n_nan": int(s.isna().sum()), "mean": s.mean(), "std": s.std(),
                "min": s.min(), "q25": s.quantile(0.25), "median": s.median(),
                "q75": s.quantile(0.75), "max": s.max(),
            })
    summary = pd.DataFrame(rows)
    summary_path = os.path.join(output_dir, f"xx_summary_{eosid}_{version}.csv")
    summary.to_csv(summary_path, index=False)

    combined = merge_combined_summary(summary)
    n_path = combined["pathogen"].nunique()
    print(f"  combined table now covers {n_path} pathogen(s) "
          f"-> {os.path.basename(COMBINED_SUMMARY_PATH)}")

    calibrated = summary[~summary["score"].str.endswith("_raw")]
    print(calibrated.to_string(
        index=False,
        columns=["query_set", "score", "n", "n_nan", "median", "mean", "min", "max"],
        float_format=lambda v: f"{v:.4f}",
    ))
    print(f"\n  -> {summary_path}")
    print("\n  Reminder: calibrated scores are CDF lookups against the reference's own "
          "distribution,\n  so the random control's medians are expected near 0.5 BY "
          "CONSTRUCTION. Interpreting\n  what the top-1000 numbers mean is a call for the "
          "user, not this script.")
    return summary_path


def plot_panel():
    """Draw the per-pathogen panel from the merged summary, once it covers all 15 pathogens."""
    print("\n[plot] per-pathogen panel")
    if not os.path.exists(COMBINED_SUMMARY_PATH):
        print(f"  [skip] {os.path.basename(COMBINED_SUMMARY_PATH)} not found — score a pathogen first")
        return
    summary = pd.read_csv(COMBINED_SUMMARY_PATH)
    missing = sorted(set(ALL_PATHOGENS) - set(summary["pathogen"]))
    if missing:
        print(f"  [skip] summary still lacks {len(missing)} pathogen(s): {', '.join(missing)}")
        return
    order = pathogen_phylo_order(ENDPOINT_SELECTION_PATH, PATHOGENS_CONFIG_PATH, TAXONOMY_PATH)
    save_eosquality_figure(summary, order, output_dir)


def main():
    if args.plot_only:
        plot_panel()
        return
    pathogen = args.pathogen
    eosid = PATHOGEN_TO_EOSID[pathogen]
    print(f"=== eosquality :: {pathogen} ({eosid}) ===")

    ensure_library()
    check_rdkit_version()

    reference_csv, version = resolve_reference_csv(eosid)
    print(f"  reference predictions: {os.path.basename(reference_csv)}")
    artifacts_dir = os.path.join(output_dir, f"artifacts_{eosid}_{version}")

    fit_reference(reference_csv, artifacts_dir)
    report_feature_selection(artifacts_dir)
    queries = build_queries(reference_csv, eosid, version, pathogen)
    outputs = score_queries(queries, artifacts_dir, eosid, version)
    summarise(outputs, eosid, version, pathogen, artifacts_dir)
    plot_panel()


if __name__ == "__main__":
    main()
