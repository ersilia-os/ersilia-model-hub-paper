"""Download and stage all input data needed for the analysis pipeline.

Run this script once before running any other analysis scripts.
Organised in four sections matching the data source type.

Section 1 – Repos / eosvc
    Data from companion git repos (expected cloned at the same level as this repo)
    or from this repo's eosvc-tracked storage. Add new sources to SECTION1_SOURCES.

Section 2 – GitHub
    Raw files fetched directly from public GitHub URLs.

Section 3 – Airtable
    Metadata downloaded from the Airtable share link.

Section 4 – Isaura
    Pre-computed model predictions fetched from the Isaura cache.

Usage
-----
    python 00_download_data.py            # copy Section 1 data from companion repos
    python 00_download_data.py --eosvc    # pull Section 1 data from eosvc instead

Outputs
-------
    data/raw/euopenscreen_data/02_merged/{code}.csv    (one per pathogen)
    data/raw/euopenscreen_data/02_binarised_assays/{assay_id}.csv
    data/raw/euopenscreen_data/primary_assays_manual.csv
    data/raw/euopenscreen_data/06_subset_data/exclusivity/{code}_{exclusive,nonexclusive}.csv
    data/raw/euopenscreen_data/06_subset_data/secondary/{code}_secondary.csv
    data/raw/coadd_data/00_smiles_info.csv                        (CoAdd screening library SMILES)
    data/raw/coadd_data/03_binarised_inhibition/{strain}.csv
    data/raw/coadd_data/05_binarised_mic/{strain}.csv
    data/raw/chembl_model_reports/{pathogen}/{name}.csv           (per-fold CV reports)
    data/raw/chembl_model_reports/{pathogen}/{name}_folds.json
    data/raw/chembl_model_reports/10_reports/10_reports.csv       (aggregated summary)
    data/raw/chembl_model_reports/10_reports/10_discarded_models.csv
    data/raw/compound_lists/reference_library_smiles.csv
    data/raw/compound_lists/drugbank_smiles.csv
    data/raw/airtable_metadata.csv
    data/processed/annotation_preds_ref_library/{model_id}_{version}.csv
"""

import argparse
import io
import json
import os
import re
import shutil
import subprocess
import sys

import pandas as pd
import requests

root = os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.abspath(os.path.join(root, ".."))
sys.path.append(os.path.join(root, "..", "src"))

from default import (
    AIRTABLE_BASE_ID,
    AIRTABLE_SHARE_URL,
    AIRTABLE_VIEW_ID,
    DRUGBANK_URL,
    REFERENCE_LIBRARY_URL,
)
from isaura_utils import download_from_isaura

parser = argparse.ArgumentParser(
    description="Download and stage all input data for the analysis pipeline."
)
parser.add_argument(
    "--eosvc",
    action="store_true",
    help="Pull Section 1 data from eosvc storage instead of companion repos.",
)
parser.add_argument(
    "--skip-isaura",
    action="store_true",
    help="Skip Section 4 (Isaura precalc predictions) and the validation summary — much faster.",
)
args = parser.parse_args()

raw_dir = os.path.join(repo_root, "data", "raw")
compound_lists_dir = os.path.join(raw_dir, "compound_lists")
annotation_dir = os.path.join(repo_root, "data", "processed", "annotation_preds_ref_library")
os.makedirs(compound_lists_dir, exist_ok=True)
os.makedirs(annotation_dir, exist_ok=True)


# =============================================================================
# Section 1 helpers
# =============================================================================

def _copy_files(src_path, dst_dir, include_fn=None, rename_fn=None, label_prefix=""):
    """Copy matching files from src_path into dst_dir with skip-if-exists."""
    os.makedirs(dst_dir, exist_ok=True)
    for fname in sorted(os.listdir(src_path)):
        if not os.path.isfile(os.path.join(src_path, fname)):
            continue
        if include_fn and not include_fn(fname):
            continue
        out_name = rename_fn(fname) if rename_fn else fname
        out_path = os.path.join(dst_dir, out_name)
        if os.path.exists(out_path):
            print(f"  Already exists: {label_prefix}{out_name}")
            continue
        shutil.copy2(os.path.join(src_path, fname), out_path)
        print(f"  -> {out_path}")


def copy_from_repo(repo_name, src_dir, dst_dir, include_fn=None, rename_fn=None, recursive=False):
    """Copy matching files from a sibling repo directory into a local directory.

    If recursive=True, iterates over immediate subdirectories of src_dir and
    copies their contents into matching subdirectories under dst_dir.
    """
    repo_path = os.path.join(repo_root, "..", repo_name)
    if not os.path.isdir(repo_path):
        print(
            f"  [SKIP] {repo_name} not found at {repo_path}\n"
            f"         Clone: git clone https://github.com/ersilia-os/{repo_name}\n"
            f"         Or use: python 00_download_data.py --eosvc"
        )
        return
    src_path = os.path.join(repo_path, src_dir)
    if not os.path.isdir(src_path):
        print(
            f"  [SKIP] {src_dir} not found in {repo_name} at {src_path}\n"
            f"         Source not published yet — keeping any existing cached copy."
        )
        return
    if recursive:
        for subdir in sorted(os.listdir(src_path)):
            if not os.path.isdir(os.path.join(src_path, subdir)):
                continue
            _copy_files(
                src_path=os.path.join(src_path, subdir),
                dst_dir=os.path.join(dst_dir, subdir),
                include_fn=include_fn,
                rename_fn=rename_fn,
                label_prefix=f"{subdir}/",
            )
    else:
        _copy_files(src_path, dst_dir, include_fn, rename_fn)


def download_from_eosvc(eosvc_path, dst_dir):
    """Pull a path from this repo's eosvc storage."""
    os.makedirs(dst_dir, exist_ok=True)
    subprocess.run(
        ["eosvc", "download", "--path", eosvc_path],
        cwd=repo_root,
        check=True,
    )


# =============================================================================
# Section 1 sources
# =============================================================================
# Each entry describes one dataset. To add a new source, append a dict with:
#   description : human-readable label printed during download
#   repo        : companion repo name (cloned at the same level as this repo)
#   src_dir     : path within the companion repo to copy from
#   dst_dir     : local destination directory (absolute path)
#   eosvc_path  : path passed to `eosvc download --path` when --eosvc is used
#   recursive   : (optional) if True, copies per-subdirectory instead of flat
#   include     : (optional) callable(filename) → bool, filters which files to copy
#   rename      : (optional) callable(filename) → str, renames files on copy
#
# NOTE: EU OpenScreen H5 prediction files are currently placed manually in
# config/eu-openscreen_preds_h5/. When published in a companion repo or
# eosvc, add an entry here.

# Whitelist of per-pathogen curation summary CSVs staged from chembl-antimicrobial-tasks
# (output/stage4/<pathogen>/). These are per-assay / per-pool SUMMARY tables, not the full
# per-molecule datasets — sufficient to rebuild the curation figures without copying ~1 GB.
_CURATION_SUMMARY_FILES = {
    "21_curation_summary.csv",
    "21_curation_stats.csv",
    "21_curation_categories.csv",
    "22_binarisation_summary.csv",
    "22_cutoff_sensitivity.csv",
    "23_pool_summary.csv",
    "23_chemspace_partition.csv",
    "23_first_pass.csv",
    "24_cv_summary.csv",
    "25_pool_summary.csv",
    "25_merge_log.csv",
    "26_cv_summary.csv",
}

SECTION1_SOURCES = [
    {
        "description": "EU OpenScreen antimicrobial tasks — merged",
        "repo": "eu-openscreen-antimicrobial-tasks",
        "src_dir": "data/processed/02_merged",
        "dst_dir": os.path.join(raw_dir, "euopenscreen_data", "02_merged"),
        "eosvc_path": "data/raw/euopenscreen_data/02_merged",
    },
    {
        "description": "EU OpenScreen antimicrobial tasks — binarised assays",
        "repo": "eu-openscreen-antimicrobial-tasks",
        "src_dir": "data/processed/02_binarised_assays",
        "dst_dir": os.path.join(raw_dir, "euopenscreen_data", "02_binarised_assays"),
        "eosvc_path": "data/raw/euopenscreen_data/02_binarised_assays",
    },
    {
        "description": "EU OpenScreen primary assays manual annotation",
        "repo": "eu-openscreen-antimicrobial-tasks",
        "src_dir": "data/config",
        "dst_dir": os.path.join(raw_dir, "euopenscreen_data"),
        "eosvc_path": "data/raw/euopenscreen_data",
        "include": lambda f: f == "primary_assays_manual.csv",
    },
    {
        "description": "EU OpenScreen subset data (exclusive/non-exclusive/secondary)",
        "repo": "eu-openscreen-antimicrobial-tasks",
        "src_dir": "output/06_subset_data",
        "dst_dir": os.path.join(raw_dir, "euopenscreen_data", "06_subset_data"),
        "eosvc_path": "data/raw/euopenscreen_data/06_subset_data",
        "recursive": True,
    },
    {
        "description": "ChEMBL antimicrobial model reports — per-fold CV",
        "repo": "chembl-antimicrobial-models",
        "src_dir": "output/09_reports",
        "dst_dir": os.path.join(raw_dir, "chembl_model_reports"),
        "eosvc_path": "data/raw/chembl_model_reports",
        "recursive": True,
    },
    {
        "description": "ChEMBL antimicrobial model reports — aggregated summary",
        "repo": "chembl-antimicrobial-models",
        "src_dir": "output/10_reports",
        "dst_dir": os.path.join(raw_dir, "chembl_model_reports", "10_reports"),
        "eosvc_path": "data/raw/chembl_model_reports/10_reports",
        "include": lambda f: f.endswith(".csv"),  # top-level CSVs only; skips plots/
    },
    # CoAdd binarised task data, staged into per-type subfolders because
    # 03_binarised_inhibition and 05_binarised_mic share 12 per-strain filenames
    # and would collide if flattened into one directory.
    {
        "description": "CoAdd binarised inhibition tasks",
        "repo": "coadd-binary-tasks",
        "src_dir": "data/processed/coadd/03_binarised_inhibition",
        "dst_dir": os.path.join(raw_dir, "coadd_data", "03_binarised_inhibition"),
        "eosvc_path": "data/raw/coadd_data/03_binarised_inhibition",
    },
    {
        "description": "CoAdd binarised MIC tasks",
        "repo": "coadd-binary-tasks",
        "src_dir": "data/processed/coadd/05_binarised_mic",
        "dst_dir": os.path.join(raw_dir, "coadd_data", "05_binarised_mic"),
        "eosvc_path": "data/raw/coadd_data/05_binarised_mic",
    },
    # Canonical CoAdd screening-library SMILES list (full ~100k compounds, columns
    # smiles,std_smiles,inchikey,mw). Feeds 04_ersilia_predictions.sh as the CoAdd
    # prediction library (predict on std_smiles). Distinct from the binarised task
    # folders above, which cover only the compounds carrying activity labels.
    {
        "description": "CoAdd screening library — canonical SMILES list",
        "repo": "coadd-binary-tasks",
        "src_dir": "data/processed/coadd",
        "dst_dir": os.path.join(raw_dir, "coadd_data"),
        "eosvc_path": "data/raw/coadd_data",
        "include": lambda f: f == "00_smiles_info.csv",
    },
    # ChEMBL data-curation summaries (chembl-antimicrobial-tasks, Stage 4). Only the small
    # per-pathogen and aggregate SUMMARY CSVs are copied — never the full cleaned/binarised
    # molecule datasets — so the curation figures (02_chembl_data_curation.py) are rebuilt
    # from summaries alone (~40 MB total). See _CURATION_SUMMARY_FILES for the whitelist.
    {
        "description": "ChEMBL curation — per-pathogen step 21-26 summaries",
        "repo": "chembl-antimicrobial-tasks",
        "src_dir": "output/stage4",
        "dst_dir": os.path.join(raw_dir, "chembl_curation"),
        "eosvc_path": "data/raw/chembl_curation",
        "recursive": True,
        "include": lambda f: f in _CURATION_SUMMARY_FILES,
    },
    {
        "description": "ChEMBL curation — step 27 aggregate tables",
        "repo": "chembl-antimicrobial-tasks",
        "src_dir": "output/stage4/general_plots",
        "dst_dir": os.path.join(raw_dir, "chembl_curation", "general"),
        "eosvc_path": "data/raw/chembl_curation/general",
        "include": lambda f: f in {
            "27_master_table.csv", "27_cutoff_sensitivity.csv",
            "27_final_data_overlap.csv", "27_chembl_space.json",
            "27_chembl_coverage.csv",
        },
    },
]


# =============================================================================
# Section 1 — Repos / eosvc
# =============================================================================

for source in SECTION1_SOURCES:
    print(f"\n{source['description']}...")
    if args.eosvc:
        download_from_eosvc(source["eosvc_path"], source["dst_dir"])

    else:
        copy_from_repo(
            repo_name=source["repo"],
            src_dir=source["src_dir"],
            dst_dir=source["dst_dir"],
            include_fn=source.get("include"),
            rename_fn=source.get("rename"),
            recursive=source.get("recursive", False),
        )


# =============================================================================
# Section 2 — GitHub
# =============================================================================
# Raw files fetched directly from public GitHub URLs.

# Reference library
ref_path = os.path.join(compound_lists_dir, "reference_library_smiles.csv")
if os.path.exists(ref_path):
    print("\nReference library already exists, skipping download.")
    df_ref = pd.read_csv(ref_path)
else:
    print("\nDownloading reference library...")
    resp = requests.get(REFERENCE_LIBRARY_URL, timeout=60)
    resp.raise_for_status()
    df_ref = pd.read_csv(io.StringIO(resp.text))
    df_ref = df_ref.rename(columns={"standardized_smiles": "input"})
    n_before = len(df_ref)
    df_ref = df_ref.drop_duplicates(subset=["input"])
    n_dropped = n_before - len(df_ref)
    if n_dropped:
        print(f"  Removed {n_dropped} duplicate SMILES from reference library.")
    df_ref.to_csv(ref_path, index=False)
    print(f"  -> {ref_path} ({len(df_ref)} compounds)")

# DrugBank
db_path = os.path.join(compound_lists_dir, "drugbank_smiles.csv")
if os.path.exists(db_path):
    print("DrugBank already exists, skipping download.")
else:
    print("Downloading DrugBank...")
    resp = requests.get(DRUGBANK_URL, timeout=60)
    resp.raise_for_status()
    df_db = pd.read_csv(io.StringIO(resp.text))
    df_db = df_db.rename(columns={"Smiles": "smiles"})
    n_before = len(df_db)
    df_db = df_db.drop_duplicates(subset=["smiles"])
    n_dropped = n_before - len(df_db)
    if n_dropped:
        print(f"  Removed {n_dropped} duplicate SMILES from DrugBank.")
    df_db.to_csv(db_path, index=False)
    print(f"  -> {db_path} ({len(df_db)} compounds)")


# =============================================================================
# Section 3 — Airtable
# =============================================================================
# Metadata downloaded from the Airtable share link.

meta_path = os.path.join(raw_dir, "airtable_metadata.csv")
if os.path.exists(meta_path):
    print("\nAirtable metadata already exists, skipping download.")
    df_meta = pd.read_csv(meta_path)
else:
    print("\nDownloading Airtable metadata...")
    session = requests.Session()
    session.headers.update({
        "User-Agent": (
            "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 "
            "Chrome/120.0.0.0 Safari/537.36"
        )
    })
    resp = session.get(AIRTABLE_SHARE_URL, timeout=30)
    resp.raise_for_status()
    m = re.search(r"initData\s*=\s*(\{)", resp.text)
    if not m:
        raise RuntimeError(
            "Could not parse Airtable share page — page structure may have changed."
        )
    init, _ = json.JSONDecoder().raw_decode(resp.text, m.start(1))
    csrf_token = init["csrfToken"]
    access_policy = init["accessPolicy"]

    csv_url = f"https://airtable.com/v0.3/view/{AIRTABLE_VIEW_ID}/downloadCsv"
    csv_resp = session.get(
        csv_url,
        params={"accessPolicy": access_policy},
        headers={
            "Accept": "*/*",
            "x-csrf-token": csrf_token,
            "x-airtable-application-id": AIRTABLE_BASE_ID,
            "x-requested-with": "XMLHttpRequest",
            "x-time-zone": "UTC",
            "x-user-locale": "en-US",
            "Referer": AIRTABLE_SHARE_URL,
        },
        timeout=60,
    )
    csv_resp.raise_for_status()
    df_meta = pd.read_csv(io.BytesIO(csv_resp.content), encoding="utf-8-sig")
    df_meta.to_csv(meta_path, index=False)
    print(f"  -> {meta_path} ({len(df_meta)} models)")

# =============================================================================
# Section 4 — Isaura
# =============================================================================
# Pre-computed model predictions fetched from the Isaura cache.

if args.skip_isaura:
    print("\n[--skip-isaura] Skipping Section 4 (Isaura predictions) and validation.")
    sys.exit(0)

print("\nDownloading annotation model predictions for reference library...")

annotation_models = df_meta[
    (df_meta["Task"] == "Annotation") & (df_meta["Status"] == "Ready")
][["Identifier", "Release"]].reset_index(drop=True)

skipped = []
for _, row in annotation_models.iterrows():
    model_id = row["Identifier"]
    version = row["Release"]
    if pd.isna(version) or not str(version).strip():
        print(f"  Skipping {model_id}: no release version in metadata.")
        skipped.append(model_id)
        continue
    version = str(version).strip()
    isaura_version = version.split(".")[0]
    output_csv = os.path.join(annotation_dir, f"{model_id}_{isaura_version}.csv")
    if os.path.exists(output_csv):
        print(f"  Already exists: {model_id} {isaura_version}")
        continue
    print(f"  Fetching {model_id} {isaura_version}...")
    try:
        download_from_isaura(
            model_id=model_id,
            model_version=isaura_version,
            input_csv=ref_path,
            output_path=annotation_dir,
        )
        print(f"    -> {output_csv}")
    except (Exception, SystemExit) as e:
        print(f"  WARNING: {model_id} {version} not available or failed: {e}")

if skipped:
    print(f"  Skipped {len(skipped)} models with no version: {skipped}")
print(f"  Done. Annotation predictions in {annotation_dir}")


# =============================================================================
# Validation
# =============================================================================

print("\n--- Validation summary ---")

n_ref = len(df_ref)
print(f"Reference library: {n_ref} compounds")

missing_files = []
incomplete = []

for _, row in annotation_models.iterrows():
    model_id = row["Identifier"]
    version = row["Release"]
    if pd.isna(version) or not str(version).strip():
        missing_files.append((model_id, "no version"))
        continue
    isaura_version = str(version).strip().split(".")[0]
    output_csv = os.path.join(annotation_dir, f"{model_id}_{isaura_version}.csv")
    if not os.path.exists(output_csv):
        missing_files.append((model_id, isaura_version))
        continue
    df_pred = pd.read_csv(output_csv)
    n_pred = df_pred["input"].nunique() if "input" in df_pred.columns else len(df_pred)
    if n_pred < n_ref:
        incomplete.append((model_id, isaura_version, n_pred, n_ref - n_pred))

if missing_files:
    print(f"\nMissing prediction files ({len(missing_files)}):")
    for model_id, ver in missing_files:
        print(f"  {model_id} ({ver})")
else:
    print(f"All {len(annotation_models)} annotation models have prediction files.")

if incomplete:
    print(f"\nIncomplete predictions ({len(incomplete)} models):")
    for model_id, ver, n_pred, n_missing in incomplete:
        print(f"  {model_id} {ver}: {n_pred}/{n_ref} compounds ({n_missing} missing)")
else:
    print("All available prediction files cover the full reference library.")
