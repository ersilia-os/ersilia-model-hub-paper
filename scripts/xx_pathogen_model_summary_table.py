"""Per-model dataset/performance summary table, one row per retained ChEMBL/
PubChem sub-model across the 15 curated pathogen models (193 rows total).

Reproduces a table that had been hand-built (as
~/Desktop/ersilia-paper/SupplementaryTable1.csv) by joining and renaming
values already available from files this repo already stages locally — no
new download needed:

  config/pathogens_of_interest.csv                                    (git-committed)
      pathogen (full species name), code, eosid.
  data/raw/chembl_model_reports/10_reports/10_reports.csv             (staged by
      00_download_data.py; same file xx_chembl_models_drugbank.py/reflib.py
      use) — one row per retained sub-model: internal name (e.g. DR_0001),
      n_compounds (final, i.e. after augmentation), n_positives,
      n_added_negatives, n_added_decoys, auroc_mean. Row order (DR_ before
      SP_, size-descending within each) is the table's canonical order —
      preserved as-is, not re-sorted.
  data/processed/annotation_preds_ref_library/{eosid}_reports.csv        (staged by
      xx_chembl_models_reflib.py's fetch_internal_to_public_mapping())
      used only for its internal-name -> public model_name mapping (e.g.
      DR_0001 -> chembl_dose_response_0) — its own row order differs from
      10_reports.csv's and is not used.
  data/raw/chembl_model_reports/09c_scaffold_vs_random/09c_delta_auroc.csv (staged
      by 00_download_data.py, from chembl-antimicrobial-models' step 09c —
      AUROC delta between random (09) and scaffold-grouped (09b) CV splits)
      used for auroc_scaffold_mean, joined on (pathogen code, name). Every one
      of the 193 retained sub-models has a match (verified); the file's 3
      extra rows belong to discarded models not in this table.

n_compounds_original = n_compounds - n_added_negatives - n_added_decoys (same
formula scripts/03_chembl_models_performance.py:241 uses, for the same reason:
n_compounds already counts added negatives/decoys). original_ratio/final_ratio
are computed by division — neither source file carries a ratio column.

Verified row-for-row (193/193 rows, all 12 columns) against the hand-built
table before this script was finalized.

Usage:
    python xx_pathogen_model_summary_table.py

Outputs:
    output/xx_pathogen_model_summary_table/pathogen_model_summary_table.csv
"""

import os
import sys

import pandas as pd

root = os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.abspath(os.path.join(root, ".."))
sys.path.append(os.path.join(root, "..", "src"))
from chembl_models_analyses_common import require

PATHOGENS_CONFIG_PATH = os.path.join(repo_root, "config", "pathogens_of_interest.csv")
REPORTS_10_PATH = os.path.join(repo_root, "data", "raw", "chembl_model_reports", "10_reports", "10_reports.csv")
REFLIB_DIR = os.path.join(repo_root, "data", "processed", "annotation_preds_ref_library")
SCAFFOLD_PATH = os.path.join(repo_root, "data", "raw", "chembl_model_reports",
                              "09c_scaffold_vs_random", "09c_delta_auroc.csv")

output_dir = os.path.join(repo_root, "output", "xx_pathogen_model_summary_table")
os.makedirs(output_dir, exist_ok=True)

pathogens_df = pd.read_csv(PATHOGENS_CONFIG_PATH)
reports_10 = pd.read_csv(require(REPORTS_10_PATH, "Run 00_download_data.py first."), dtype={"name": str})
scaffold_df = pd.read_csv(require(SCAFFOLD_PATH, "Run 00_download_data.py first."), dtype={"name": str})
scaffold_lookup = {
    (row.pathogen, row.name): row.auroc_scaffold_mean for row in scaffold_df.itertuples()
}

rows = []
for _, prow in pathogens_df.iterrows():
    reports_path = require(
        os.path.join(REFLIB_DIR, f"{prow.eosid}_reports.csv"),
        "Run xx_chembl_models_reflib.py first to fetch each eos repo's model/checkpoints/reports.csv.",
    )
    name_to_public = pd.read_csv(reports_path).set_index("original_name")["model_name"].to_dict()

    pathogen_reports = reports_10[reports_10["pathogen"] == prow.code]
    for _, r in pathogen_reports.iterrows():
        dataset_name = r["name"]  # r.name is the row INDEX, not the "name" column — must subscript
        n_compounds_final = int(r.n_compounds)
        n_added_inactives = int(r.n_added_negatives)
        n_added_decoys = int(r.n_added_decoys)
        n_compounds_original = n_compounds_final - n_added_inactives - n_added_decoys

        rows.append({
            "pathogen": prow.pathogen,
            "eos_id": prow.eosid,
            "dataset_name": dataset_name,
            "model_name": name_to_public[dataset_name],
            "n_compounds_original": n_compounds_original,
            "n_positives": int(r.n_positives),
            "original_ratio": round(r.n_positives / n_compounds_original, 4),
            "n_added_inactives": n_added_inactives,
            "n_added_decoys": n_added_decoys,
            "n_compounds_final": n_compounds_final,
            "final_ratio": round(r.n_positives / n_compounds_final, 4),
            "auroc_mean": r.auroc_mean,
            "auroc_mean_scaffold": scaffold_lookup.get((prow.code, dataset_name)),
        })

table = pd.DataFrame(rows)
n_missing_scaffold = table["auroc_mean_scaffold"].isna().sum()
if n_missing_scaffold:
    print(f"WARNING: {n_missing_scaffold} rows have no scaffold-split AUROC match — left as NaN, not dropped.")

out_path = os.path.join(output_dir, "pathogen_model_summary_table.csv")
table.to_csv(out_path, index=False)
print(f"Saved: {out_path}  ({len(table)} rows across {table['pathogen'].nunique()} pathogens)")
