"""Download and stage all input data needed for the analysis pipeline.

Run this script once before running any other analysis scripts.
Organised in four sections matching the data source type.

Section 1 – Companion repos (via their eosvc)
    Published data of companion repos, downloaded from their eosvc storage in the public
    `eosvc-public` S3 bucket. Sibling clones are never read or written. Add new sources to
    SECTION1_SOURCES.

Section 2 – GitHub
    Raw files fetched directly from public GitHub URLs.

Section 3 – Airtable
    Metadata downloaded from the Airtable share link into a dated, frozen file. The snapshot date
    lives in src/default.py (AIRTABLE_SNAPSHOT_DATE); bump it there to take a new snapshot.

Section 4 – Isaura
    Pre-computed model predictions fetched from the Isaura cache.

Usage
-----
    python 00_download_data.py                 # all four sections
    python 00_download_data.py --skip-isaura   # Sections 1-3 only (much faster)

Every section is skip-if-exists. To re-stage something after an upstream update, delete the local
file(s) and re-run. This repo's own eosvc storage is not used here: pull or push it with the
`eosvc` CLI directly.

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
    data/raw/lazyqsar_benchmark/all_results.csv
    data/raw/airtable_metadata_{YYYYMMDD}.csv   (frozen snapshot, date set in src/default.py)
    data/processed/annotation_preds_ref_library/{model_id}_{version}.csv
    data/processed/annotation_preds_ref_library/eos4djh_v1.csv
    data/processed/eos1klk_projection/eos1klk_v1.csv
"""

import argparse
import io
import json
import os
import re
import sys
import xml.etree.ElementTree as ET
from urllib.parse import quote

import pandas as pd
import requests

root = os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.abspath(os.path.join(root, ".."))
sys.path.append(os.path.join(root, "..", "src"))

from default import (
    AIRTABLE_BASE_ID,
    AIRTABLE_METADATA_FILE,
    AIRTABLE_SHARE_URL,
    AIRTABLE_SNAPSHOT_DATE,
    AIRTABLE_VIEW_ID,
    DRUGBANK_URL,
    EOSVC_PUBLIC_URL,
    PHYSCHEM_MODEL_ID,
    PROJECTION_MODEL_ID,
    PROJECTION_PREDS_SUBDIR,
    REFERENCE_LIBRARY_URL,
)
from airtable_utils import load_airtable_metadata
from isaura_utils import download_from_isaura

parser = argparse.ArgumentParser(
    description="Download and stage all input data for the analysis pipeline."
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
projection_dir = os.path.join(repo_root, "data", "processed", PROJECTION_PREDS_SUBDIR)
os.makedirs(compound_lists_dir, exist_ok=True)
os.makedirs(annotation_dir, exist_ok=True)
os.makedirs(projection_dir, exist_ok=True)


# =============================================================================
# Section 1 helpers
# =============================================================================
# Companion-repo data is read ONLY from each repo's published eosvc storage in the public
# `eosvc-public` S3 bucket, over plain HTTPS (no credentials, no sibling clone needed). Sibling
# working copies are never read and never written: what is staged here is always the published
# version. Objects live at {EOSVC_PUBLIC_URL}/{repo}/{path-within-repo}.

_S3_NS = {"s3": "http://s3.amazonaws.com/doc/2006-03-01/"}


def list_eosvc_files(repo, src_dir, recursive=False):
    """List published files under a companion repo folder in the public eosvc bucket.

    Parameters
    ----------
    repo : str
        Companion repo name, which is also its top-level prefix in the bucket.
    src_dir : str
        Folder within the companion repo, e.g. ``output/10_reports``.
    recursive : bool
        If False, return only files directly in ``src_dir``. If True, return only files one
        subfolder down (``{subdir}/{file}``), mirroring the per-pathogen folder layout.

    Returns
    -------
    list of str
        Paths relative to ``src_dir``, sorted.
    """
    prefix = f"{repo}/{src_dir.strip('/')}/"
    params = {"list-type": "2", "prefix": prefix}
    if not recursive:
        params["delimiter"] = "/"
    rel_paths = []
    while True:
        resp = requests.get(f"{EOSVC_PUBLIC_URL}/", params=params, timeout=60)
        resp.raise_for_status()
        tree = ET.fromstring(resp.content)
        for key in tree.iterfind("s3:Contents/s3:Key", _S3_NS):
            rel = key.text[len(prefix):]
            depth = rel.count("/")
            if rel and not rel.endswith("/") and depth == (1 if recursive else 0):
                rel_paths.append(rel)
        token = tree.find("s3:NextContinuationToken", _S3_NS)
        if token is None:
            break
        params["continuation-token"] = token.text
    return sorted(rel_paths)


def fetch_from_eosvc(repo, src_dir, dst_dir, include_fn=None, rename_fn=None, recursive=False):
    """Download a companion repo's published files into dst_dir, skipping files already present.

    To re-stage a file (e.g. after an upstream update), delete the local copy and re-run.
    Files are written to a temporary name first, so an interrupted download never leaves a
    truncated file that a later run would mistake for a complete one.
    """
    rel_paths = list_eosvc_files(repo, src_dir, recursive=recursive)
    if include_fn:
        rel_paths = [r for r in rel_paths if include_fn(os.path.basename(r))]
    base_url = f"{EOSVC_PUBLIC_URL}/{repo}/{src_dir.strip('/')}"
    print(f"  Source: {base_url}/")
    if not rel_paths:
        print("  [SKIP] No matching files published in eosvc — keeping any existing local copy.")
        return
    n_existing = 0
    for rel in rel_paths:
        subdir, fname = os.path.split(rel)
        out_name = rename_fn(fname) if rename_fn else fname
        out_path = os.path.join(dst_dir, subdir, out_name)
        if os.path.exists(out_path):
            n_existing += 1
            continue
        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        resp = requests.get(f"{base_url}/{quote(rel)}", stream=True, timeout=120)
        resp.raise_for_status()
        tmp_path = out_path + ".part"
        with open(tmp_path, "wb") as f:
            f.writelines(resp.iter_content(chunk_size=1 << 20))
        os.replace(tmp_path, out_path)
        print(f"  -> {os.path.relpath(out_path, repo_root)}")
    print(f"  {len(rel_paths) - n_existing} downloaded, {n_existing} already present "
          f"(of {len(rel_paths)} published).")


# =============================================================================
# Section 1 sources
# =============================================================================
# Each entry describes one dataset. To add a new source, append a dict with:
#   description : human-readable label printed during download
#   repo        : companion repo name (its prefix in the public eosvc bucket)
#   src_dir     : path within the companion repo to download from
#   dst_dir     : local destination directory (absolute path)
#   recursive   : (optional) if True, downloads per-subdirectory instead of flat
#   include     : (optional) callable(filename) → bool, filters which files to download
#   rename      : (optional) callable(filename) → str, renames files on download
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
    },
    # NOTE: these per-assay files are MEASUREMENT-level — one row per assay well, so a structure the
    # library registered as two separate compounds appears as two rows with the same SMILES. Upstream
    # applies its "Active prevails" (max bin) rule only when building `02_merged`
    # (`scripts/02_binarise_and_merge.py::merge_pathogen_rows`), never to these files, so
    # `load_euos_primary` re-applies the SAME rule on read to get one row per compound. If upstream
    # ever ships a compound-level per-assay output, that dedupe becomes a no-op (its log line will
    # report 0 collapses) and can be dropped.
    #
    # NOTE: two upstream metadata artefacts, neither of which affects the numbers we compute, but both
    # of which make `assays_annotated_manual.csv` untrustworthy as a COUNT source. Root cause for both:
    # `00_extract_assays.py` carries the ECBD `deprecated` field through to `00_assay_summary.csv` as
    # metadata but never FILTERS on it, and the hand-built annotation sheet inherited the extra rows.
    #
    #  1. EOS300078 (S. aureus primary) is listed TWICE — 981 actives and 379. Not a transcription
    #     error: the ECBD dump holds a SUPERSEDED record (deprecated 2025-02-10, intended_target
    #     EOS300077, 981 actives) beside the current one (intended_target EOS300170, 379). 379 is the
    #     live figure and the one we use, confirmed independently against the raw extract:
    #     `activity == "active"` is 379 AND `value >= 70` is 379, matching the assay's documented
    #     "active = >=70% inhibition" cutoff. 981 corresponds to no cutoff in the data (>=35% gives
    #     1020, >=40% gives 826) — it is a stale count from the superseded analysis.
    #  2. EOS300161 (E. coli, non-primary) is a fully deprecated assay — deprecated 2014-11-20, with
    #     no current record under that id. It was re-registered as EOS300159: identical compound set
    #     (3,866, 100% overlap) and identical activity labels, only the `value` numbers differ. Both
    #     ids are still extracted and annotated, so E. coli looks like 6 assays but has 5 live ones.
    #
    # Neither touches our analysis: we take labels from `primary_assays_manual.csv` (one correct row
    # per pathogen, saureus = 379) and EOS300161 is not a primary assay, so it reaches only `02_merged`
    # (InChIKey lookup, set-valued) and the unused `secondary/` subset. Report upstream; do not
    # work around it here.
    {
        "description": "EU OpenScreen antimicrobial tasks — binarised assays",
        "repo": "eu-openscreen-antimicrobial-tasks",
        "src_dir": "data/processed/02_binarised_assays",
        "dst_dir": os.path.join(raw_dir, "euopenscreen_data", "02_binarised_assays"),
    },
    {
        "description": "EU OpenScreen primary assays manual annotation",
        "repo": "eu-openscreen-antimicrobial-tasks",
        "src_dir": "data/config",
        "dst_dir": os.path.join(raw_dir, "euopenscreen_data"),
        "include": lambda f: f == "primary_assays_manual.csv",
    },
    # NOTE: only the `exclusivity/` subtree (112 KB) is read by any code — `_load_exclusivity_task` in
    # src/eval_euopenscreen.py. `secondary/` (2.7 MB) and `abx_similarity/` (50 MB) are staged but
    # currently UNUSED. Left in place deliberately rather than filtered out: they are cheap to keep
    # next to the subset they belong with, and dropping a staged input is a data decision. Narrow the
    # `include` here if the 52.7 MB ever matters.
    #
    # These subsets are COMPOUND-level (one row per structure), unlike the measurement-level per-assay
    # files above. That mismatch is why `load_euos_primary` collapses duplicate SMILES — without it the
    # own-assay active count and `exclusive + shared` disagreed by up to 1 per organism.
    {
        "description": "EU OpenScreen subset data (exclusive/non-exclusive/secondary)",
        "repo": "eu-openscreen-antimicrobial-tasks",
        "src_dir": "output/06_subset_data",
        "dst_dir": os.path.join(raw_dir, "euopenscreen_data", "06_subset_data"),
        "recursive": True,
    },
    {
        "description": "ChEMBL antimicrobial model reports — per-fold CV",
        "repo": "chembl-antimicrobial-models",
        "src_dir": "output/09_reports",
        "dst_dir": os.path.join(raw_dir, "chembl_model_reports"),
        "recursive": True,
    },
    {
        "description": "ChEMBL antimicrobial model reports — aggregated summary",
        "repo": "chembl-antimicrobial-models",
        "src_dir": "output/10_reports",
        "dst_dir": os.path.join(raw_dir, "chembl_model_reports", "10_reports"),
        "include": lambda f: f.endswith(".csv"),  # top-level CSVs only; skips plots/
    },
    {
        "description": "ChEMBL antimicrobial models — scaffold-split vs random-split AUROC",
        "repo": "chembl-antimicrobial-models",
        "src_dir": "output/09c_scaffold_vs_random",
        "dst_dir": os.path.join(raw_dir, "chembl_model_reports", "09c_scaffold_vs_random"),
        "include": lambda f: f == "09c_delta_auroc.csv",  # skip the PNG, only need the CSV
    },
    # CoAdd binarised task data, staged into per-type subfolders because
    # 03_binarised_inhibition and 05_binarised_mic share 12 per-strain filenames
    # and would collide if flattened into one directory.
    {
        "description": "CoAdd binarised inhibition tasks",
        "repo": "coadd-binary-tasks",
        "src_dir": "data/processed/coadd/03_binarised_inhibition",
        "dst_dir": os.path.join(raw_dir, "coadd_data", "03_binarised_inhibition"),
    },
    {
        "description": "CoAdd binarised MIC tasks",
        "repo": "coadd-binary-tasks",
        "src_dir": "data/processed/coadd/05_binarised_mic",
        "dst_dir": os.path.join(raw_dir, "coadd_data", "05_binarised_mic"),
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
        "recursive": True,
        "include": lambda f: f in _CURATION_SUMMARY_FILES,
    },
    {
        "description": "ChEMBL curation — step 27 aggregate tables",
        "repo": "chembl-antimicrobial-tasks",
        "src_dir": "output/stage4/general_plots",
        "dst_dir": os.path.join(raw_dir, "chembl_curation", "general"),
        "include": lambda f: f in {
            "27_master_table.csv", "27_cutoff_sensitivity.csv",
            "27_final_data_overlap.csv", "27_chembl_space.json",
            "27_chembl_coverage.csv",
        },
    },
    {
        "description": "LazyQSAR v3.4.2 benchmark results",
        "repo": "ersilia-ml-benchmark",
        "src_dir": "output/lazyqsar/v3.4.2",
        "dst_dir": os.path.join(raw_dir, "lazyqsar_benchmark"),
        "include": lambda f: f == "all_results.csv",
    },
]


# =============================================================================
# Section 1 — Companion repos (via their eosvc)
# =============================================================================

# NOTE: the eosquality reference library (~2.2 GB: Morgan FP index, 50-NN self-neighbour
# arrays, physchem/MACCS matrices, plus the library SMILES CSV) is NOT fetched here. It
# lands in the ~/.eosquality/ PACKAGE CACHE, not in this repo's data/ tree — the same
# status as the ersilia models under ~/eos, which are likewise not routed through this
# script. scripts/xx_eosquality.py fetches it on first run.
#
# Two things to know if that fetch ever needs to move here (user-directed, 2026-09-15):
#   - eosquality 0.1.0 has a BROKEN default URL: library/identity.py:52 points at
#     .../eosvc-public/eosquality/indices/, but eosvc publishes to
#     .../eosvc-public/eosquality/data/indices/ (missing `data/` segment; every file 403s).
#     Override it with EOSQUALITY_REFERENCE_BASE_URL — that single env var fixes both the
#     index folder and the sibling libraries/ CSV.
#   - the index's smiles.csv is byte-identical to the reference_library_smiles.csv fetched
#     in Section 2 below, apart from the `input` -> `smiles` header rename. The other
#     ~2.1 GB are eosquality-specific precalculations that exist nowhere in this repo.

for source in SECTION1_SOURCES:
    print(f"\n{source['description']}...")
    fetch_from_eosvc(
        repo=source["repo"],
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
# Metadata downloaded from the Airtable share link into a dated, frozen file named by
# AIRTABLE_METADATA_FILE. A new snapshot only happens when AIRTABLE_SNAPSHOT_DATE is bumped by hand
# in src/default.py; an existing dated file is never re-downloaded or overwritten.
#
# The file is saved exactly as exported (raw data is read-only). `keep_default_na=False` is
# load-bearing: `None` is a legitimate License value, and pandas' default NA set would blank it on
# the round-trip. The Status remap (AIRTABLE_STATUS_REMAP) is applied on load, not here.
snapshot_path = os.path.join(raw_dir, AIRTABLE_METADATA_FILE)
if os.path.exists(snapshot_path):
    print(f"\n{AIRTABLE_METADATA_FILE} already exists, skipping download.")
else:
    print(f"\nDownloading Airtable metadata (new snapshot {AIRTABLE_SNAPSHOT_DATE})...")
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
    df_raw = pd.read_csv(
        io.BytesIO(csv_resp.content), encoding="utf-8-sig", keep_default_na=False
    )
    df_raw.to_csv(snapshot_path, index=False)
    print(f"  -> {snapshot_path} ({len(df_raw)} models)")

df_meta = load_airtable_metadata(snapshot_path)
print(f"  Metadata {AIRTABLE_METADATA_FILE} (Airtable pull {AIRTABLE_SNAPSHOT_DATE}): "
      f"{len(df_meta)} models, {(df_meta['Status'] == 'Ready').sum()} Ready after remap")

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

# eos4djh (datamol-basic-descriptors, Task=Representation/Featurization) is not "Annotation", so
# the loop above never sees it — fetched here explicitly instead. Release is pinned to v1 (Airtable
# "v1.1.0"), following the same major-version-only convention as the annotation loop. Unlike the
# projector below it lands in annotation_dir, because step 08 reads every property model from that
# one folder; its Task is recorded in src/default.py, not in the path.
print(f"\nDownloading {PHYSCHEM_MODEL_ID} (basic physicochemical descriptors, reference library)...")
physchem_output_csv = os.path.join(annotation_dir, f"{PHYSCHEM_MODEL_ID}_v1.csv")
if os.path.exists(physchem_output_csv):
    print(f"  Already exists: {PHYSCHEM_MODEL_ID} v1")
else:
    try:
        download_from_isaura(
            model_id=PHYSCHEM_MODEL_ID,
            model_version="v1",
            input_csv=ref_path,
            output_path=annotation_dir,
        )
        print(f"    -> {physchem_output_csv}")
    except (Exception, SystemExit) as e:
        print(f"  WARNING: {PHYSCHEM_MODEL_ID} not available or failed: {e}")

# eos1klk (2D projector, Task=Representation/Projection) is not "Annotation", so the loop above
# never sees it — fetched here explicitly instead. Release is pinned to v1 (Airtable "v1.2.0"),
# following the same major-version-only convention as the annotation loop.
print(f"\nDownloading {PROJECTION_MODEL_ID} (2D projector, reference library)...")
projection_output_csv = os.path.join(projection_dir, f"{PROJECTION_MODEL_ID}_v1.csv")
if os.path.exists(projection_output_csv):
    print(f"  Already exists: {PROJECTION_MODEL_ID} v1")
else:
    try:
        download_from_isaura(
            model_id=PROJECTION_MODEL_ID,
            model_version="v1",
            input_csv=ref_path,
            output_path=projection_dir,
        )
        print(f"    -> {projection_output_csv}")
    except (Exception, SystemExit) as e:
        print(f"  WARNING: {PROJECTION_MODEL_ID} not available or failed: {e}")


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
