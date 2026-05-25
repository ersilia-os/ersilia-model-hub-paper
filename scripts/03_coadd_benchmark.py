"""Evaluate CoAdd model (eos3dys) endpoints against EU OpenScreen binary tasks.

For each CoAdd endpoint, computes AUROC against each pathogen's primary assay
(as defined in data/raw/euopenscreen_tasks/primary_assays_manual.csv). Compounds
present in the CoAdd training set for a given endpoint are removed before computing
AUROC (leakage check). Requires the coadd_training data downloaded via 00_download_data.py.

Requires
--------
    data/processed/02_euopenscreen_preds/eos3dys.csv
    data/raw/euopenscreen_tasks/primary_assays_manual.csv
    data/raw/euopenscreen_tasks/02_binarised_assays/{assay_eos_id}.csv
    data/raw/euopenscreen_tasks/02_merged/02_{code}.csv  (for InChIKey enrichment)
    data/raw/euopenscreen_tasks/06_subset_data/exclusivity/{code}_{exclusive,nonexclusive}.csv
    data/raw/euopenscreen_tasks/06_subset_data/secondary/{code}_secondary.csv
    data/raw/coadd_training/{strain}.csv

Outputs
-------
    output/03_coadd_benchmark/auroc_matrix.csv
    output/03_coadd_benchmark/auroc_heatmap.png
    output/03_coadd_benchmark/auroc_swarm_same_vs_different.png
    output/03_coadd_benchmark/auroc_matrix_{exclusive,nonexclusive,secondary}.csv
    output/03_coadd_benchmark/auroc_heatmap_{exclusive,nonexclusive,secondary}.png
    output/03_coadd_benchmark/auroc_heatmap_exclusive_vs_nonexclusive.png
    output/03_coadd_benchmark/auroc_heatmap_primary_vs_secondary.png
"""

import os
import sys

import numpy as np
import pandas as pd
import stylia
from sklearn.metrics import roc_auc_score

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))

from default import COADD_MODEL_ID
from plotting_utils import abbrev, plot_auroc_heatmap, plot_same_vs_diff_swarm

pred_csv = os.path.join(
    root, "..", "data", "processed", "02_euopenscreen_preds", f"{COADD_MODEL_ID}.csv"
)
tasks_dir = os.path.join(root, "..", "data", "raw", "euopenscreen_tasks")
binarised_dir = os.path.join(tasks_dir, "02_binarised_assays")
merged_dir = os.path.join(tasks_dir, "02_merged")
subset_dir = os.path.join(tasks_dir, "06_subset_data")
primary_csv = os.path.join(tasks_dir, "primary_assays_manual.csv")
coadd_training_dir = os.path.join(root, "..", "data", "raw", "coadd_training")
output_dir = os.path.join(root, "..", "output", "03_coadd_benchmark")

os.makedirs(output_dir, exist_ok=True)


# --- Step 1: Load CoAdd predictions ---

if not os.path.exists(pred_csv):
    print(
        f"Predictions not found at {pred_csv}. "
        "Run 02_euopenscreen_auroc.py first to convert the H5 file."
    )
    sys.exit(1)

pred_df = pd.read_csv(pred_csv)
endpoint_cols = [c for c in pred_df.columns if c not in ("key", "smiles")]
print(f"[{COADD_MODEL_ID}] {len(pred_df)} compounds, {len(endpoint_cols)} endpoints")


# --- Step 2: Build endpoint → training InChIKeys map (for leakage removal) ---

# Endpoint names follow the pattern {strain}_{metric}; strip the metric suffix
# to identify the matching CoAdd training file.
ENDPOINT_SUFFIXES = ["_inhib_50", "_mic_25", "_ic50"]


def _endpoint_to_strain(endpoint):
    for suffix in ENDPOINT_SUFFIXES:
        if endpoint.endswith(suffix):
            return endpoint[: -len(suffix)]
    return None


endpoint_train_keys = {}
coadd_training_available = os.path.isdir(coadd_training_dir)

for ep in endpoint_cols:
    if not coadd_training_available:
        endpoint_train_keys[ep] = set()
        continue
    strain = _endpoint_to_strain(ep)
    if strain is None:
        print(f"[leakage] {ep}: no known metric suffix — skipping leakage check")
        endpoint_train_keys[ep] = set()
        continue
    train_file = os.path.join(coadd_training_dir, f"{strain}.csv")
    if not os.path.exists(train_file):
        print(f"[leakage] {ep}: training file not found ({strain}.csv) — skipping leakage check")
        endpoint_train_keys[ep] = set()
        continue
    train_df = pd.read_csv(train_file, usecols=["inchikey"])
    endpoint_train_keys[ep] = set(train_df["inchikey"].dropna())
    print(f"[leakage] {ep}: {len(endpoint_train_keys[ep])} training compounds")

if not coadd_training_available:
    print(
        f"[WARN] CoAdd training data not found at {coadd_training_dir}. "
        "Run 00_download_data.py to enable leakage removal."
    )


# --- Step 3: Load EU OpenScreen tasks (primary assay per pathogen) ---

if not os.path.exists(primary_csv):
    print(
        "primary_assays_manual.csv not found. "
        "Run 00_download_data.py to fetch eu-openscreen-antimicrobial-tasks."
    )
    sys.exit(0)

primary_df = pd.read_csv(primary_csv)
task_files = {}
for _, mrow in primary_df.iterrows():
    tpath = os.path.join(binarised_dir, f"{mrow['assay_eos_id']}.csv")
    if os.path.exists(tpath):
        task_files[abbrev(mrow["pathogen"])] = {
            "path": tpath,
            "code": mrow["pathogen_code"],
        }

if not task_files:
    print(
        "No EU OpenScreen task files found. "
        "Run 00_download_data.py to fetch eu-openscreen-antimicrobial-tasks."
    )
    sys.exit(0)

print(f"\nEvaluating {len(endpoint_cols)} endpoints × {len(task_files)} tasks...")


# --- Step 4: Compute AUROC for each (endpoint, task) pair ---

matrix_records = []
for ep in endpoint_cols:
    train_keys = endpoint_train_keys.get(ep, set())
    ep_pred = pred_df[["smiles", ep]].rename(columns={ep: "score"})
    row = {"endpoint": ep}

    for tname, tinfo in task_files.items():
        task_df = pd.read_csv(tinfo["path"])  # smiles, bin
        mpath = os.path.join(merged_dir, f"02_{tinfo['code']}.csv")
        if os.path.exists(mpath):
            inchi_df = pd.read_csv(mpath, usecols=["smiles", "inchikey"])
            task_df = task_df.merge(inchi_df, on="smiles", how="left")
        merged = task_df.merge(ep_pred, on="smiles", how="inner")
        merged = merged[merged["bin"].isin([0, 1])]
        if train_keys and "inchikey" in merged.columns:
            n_before = len(merged)
            merged = merged[~merged["inchikey"].isin(train_keys)]
            n_removed = n_before - len(merged)
            if n_removed:
                print(f"  [{ep} / {tname}] removed {n_removed} training overlap compounds")
        if merged["bin"].nunique() < 2:
            row[tname] = float("nan")
            print(f"  [{ep} / {tname}] insufficient classes after deduplication — NaN")
        else:
            n_pos = int(merged["bin"].sum())
            n_neg = int((merged["bin"] == 0).sum())
            auc = roc_auc_score(merged["bin"], merged["score"])
            row[tname] = round(auc, 4)
            print(f"  [{ep} / {tname}] AUROC={auc:.4f} ({n_pos}+ {n_neg}-)")

    matrix_records.append(row)

matrix_df = pd.DataFrame(matrix_records).set_index("endpoint")
matrix_df.to_csv(os.path.join(output_dir, "auroc_matrix.csv"))
print("\n--- AUROC matrix (deduplicated) ---")
print(matrix_df.to_string())


# Endpoint prefix → abbreviated pathogen label (must match matrix columns)
_ENDPOINT_TO_PATHOGEN = {
    "abaumannii": "A. baumannii",
    "calbicans": "C. albicans",
    "ecoli": "E. coli",
    "efaecium": "E. faecium",
    "kpneumoniae": "K. pneumoniae",
    "paeruginosa": "P. aeruginosa",
    "saureus": "S. aureus",
}


def _build_highlight_cells(mdf):
    """Return (row_idx, col_idx) pairs where endpoint prefix matches task pathogen."""
    cols = list(mdf.columns)
    cells = []
    for row_idx, ep in enumerate(mdf.index):
        for prefix, pathogen in _ENDPOINT_TO_PATHOGEN.items():
            if ep.startswith(prefix + "_") or ep == prefix:
                if pathogen in cols:
                    cells.append((row_idx, cols.index(pathogen)))
                break
    return cells


# --- Step 5: Heatmap ---

# Format: print | Style: article
stylia.set_format("print")
stylia.set_style("article")

fig, axs = stylia.create_figure(1, 1, width=0.5, height=0.5)
plot_auroc_heatmap(
    axs.next(),
    matrix_df,
    title="CoAdd (eos3dys) AUROC per endpoint vs EU OpenScreen tasks",
    highlight_cells=_build_highlight_cells(matrix_df),
)
stylia.save_figure(os.path.join(output_dir, "auroc_heatmap.png"))


# --- Step 5b: Same vs different pathogen AUROC ---

same_aurocs = []
diff_aurocs = []
for ep, row in matrix_df.iterrows():
    matched_pathogen = None
    for prefix, pathogen in _ENDPOINT_TO_PATHOGEN.items():
        if ep.startswith(prefix + "_") or ep == prefix:
            matched_pathogen = pathogen
            break
    for col in matrix_df.columns:
        v = row[col]
        if np.isnan(v):
            continue
        if col == matched_pathogen:
            same_aurocs.append(v)
        else:
            diff_aurocs.append(v)

fig, axs = stylia.create_figure(1, 1)
plot_same_vs_diff_swarm(
    axs.next(),
    same_aurocs,
    diff_aurocs,
    title="CoAdd (eos3dys) — same vs different pathogen AUROC",
)
stylia.save_figure(os.path.join(output_dir, "auroc_swarm_same_vs_different.png"))


# --- Step 6: Subset analyses (exclusive, nonexclusive, secondary) ---

def _load_subset_task(primary_csv_path, merged_path, mode, subset_dir, code):
    """Load eval task DataFrame for a given analysis mode.

    Returns a DataFrame with smiles, bin, inchikey (where available), or None.
    exclusive/nonexclusive: actives from subset file + primary inactives.
    secondary: secondary assay file with its own bin labels.
    """
    if mode == "secondary":
        path = os.path.join(subset_dir, "secondary", f"{code}_secondary.csv")
        return pd.read_csv(path) if os.path.exists(path) else None
    act_path = os.path.join(subset_dir, "exclusivity", f"{code}_{mode}.csv")
    if not os.path.exists(act_path):
        return None
    actives = pd.read_csv(act_path)  # smiles, inchikey
    if actives.empty:
        return None
    actives = actives.assign(bin=1)
    primary = pd.read_csv(primary_csv_path)  # smiles, bin
    if os.path.exists(merged_path):
        inchi_df = pd.read_csv(merged_path, usecols=["smiles", "inchikey"])
        primary = primary.merge(inchi_df, on="smiles", how="left")
    inactives = primary[primary["bin"] == 0].copy()
    return pd.concat([actives, inactives], ignore_index=True)


_SUBSET_TITLES = {
    "exclusive": "CoAdd (eos3dys) AUROC — exclusive actives vs primary inactives",
    "nonexclusive": "CoAdd (eos3dys) AUROC — non-exclusive actives vs primary inactives",
    "secondary": "CoAdd (eos3dys) AUROC — secondary assay",
}

for mode in ("exclusive", "nonexclusive", "secondary"):
    subset_task_dfs = {}
    for _, mrow in primary_df.iterrows():
        task_df_s = _load_subset_task(
            primary_csv_path=os.path.join(binarised_dir, f"{mrow['assay_eos_id']}.csv"),
            merged_path=os.path.join(merged_dir, f"02_{mrow['pathogen_code']}.csv"),
            mode=mode,
            subset_dir=subset_dir,
            code=mrow["pathogen_code"],
        )
        if task_df_s is not None and not task_df_s.empty:
            subset_task_dfs[abbrev(mrow["pathogen"])] = task_df_s

    if not subset_task_dfs:
        print(f"\n[{mode}] No task data available — subset data may not be downloaded yet.")
        continue

    print(f"\nEvaluating {len(endpoint_cols)} endpoints × {len(subset_task_dfs)} {mode} tasks...")
    subset_matrix_records = []
    for ep in endpoint_cols:
        train_keys_s = endpoint_train_keys.get(ep, set())
        ep_pred_s = pred_df[["smiles", ep]].rename(columns={ep: "score"})
        row_s = {"endpoint": ep}

        for tname, task_df_s in subset_task_dfs.items():
            merged_s = task_df_s.merge(ep_pred_s, on="smiles", how="inner")
            merged_s = merged_s[merged_s["bin"].isin([0, 1])]
            if train_keys_s and "inchikey" in merged_s.columns:
                merged_s = merged_s[~merged_s["inchikey"].isin(train_keys_s)]
            if merged_s["bin"].nunique() < 2:
                row_s[tname] = float("nan")
            else:
                auc_s = roc_auc_score(merged_s["bin"], merged_s["score"])
                row_s[tname] = round(auc_s, 4)

        subset_matrix_records.append(row_s)

    subset_matrix_df = pd.DataFrame(subset_matrix_records).set_index("endpoint")
    subset_matrix_df.to_csv(os.path.join(output_dir, f"auroc_matrix_{mode}.csv"))
    print(f"[{mode}] matrix saved → auroc_matrix_{mode}.csv")

    fig, axs = stylia.create_figure(1, 1, width=0.5, height=0.5)
    plot_auroc_heatmap(
        axs.next(), subset_matrix_df,
        title=_SUBSET_TITLES[mode],
        highlight_cells=_build_highlight_cells(subset_matrix_df),
    )
    stylia.save_figure(os.path.join(output_dir, f"auroc_heatmap_{mode}.png"))


# --- Step 7: Paired comparison heatmaps ---

def _load_matrix(path):
    if not os.path.exists(path):
        return None
    return pd.read_csv(path, index_col=0)


_m_primary = _load_matrix(os.path.join(output_dir, "auroc_matrix.csv"))
_m_exclusive = _load_matrix(os.path.join(output_dir, "auroc_matrix_exclusive.csv"))
_m_nonexclusive = _load_matrix(os.path.join(output_dir, "auroc_matrix_nonexclusive.csv"))
_m_secondary = _load_matrix(os.path.join(output_dir, "auroc_matrix_secondary.csv"))

if _m_exclusive is not None and _m_nonexclusive is not None:
    fig, axs = stylia.create_figure(1, 2)
    plot_auroc_heatmap(
        axs.next(), _m_exclusive,
        title="CoAdd (eos3dys) — exclusive actives",
        highlight_cells=_build_highlight_cells(_m_exclusive),
    )
    plot_auroc_heatmap(
        axs.next(), _m_nonexclusive,
        title="CoAdd (eos3dys) — non-exclusive actives",
        highlight_cells=_build_highlight_cells(_m_nonexclusive),
    )
    stylia.save_figure(os.path.join(output_dir, "auroc_heatmap_exclusive_vs_nonexclusive.png"))
    print("\n[Step 7] Saved auroc_heatmap_exclusive_vs_nonexclusive.png")
else:
    print("\n[Step 7] exclusive or nonexclusive matrix not available — skipping paired heatmap.")

if _m_primary is not None and _m_secondary is not None:
    fig, axs = stylia.create_figure(1, 2)
    plot_auroc_heatmap(
        axs.next(), _m_primary,
        title="CoAdd (eos3dys) — primary assay",
        highlight_cells=_build_highlight_cells(_m_primary),
    )
    plot_auroc_heatmap(
        axs.next(), _m_secondary,
        title="CoAdd (eos3dys) — secondary assay",
        highlight_cells=_build_highlight_cells(_m_secondary),
    )
    stylia.save_figure(os.path.join(output_dir, "auroc_heatmap_primary_vs_secondary.png"))
    print("[Step 7] Saved auroc_heatmap_primary_vs_secondary.png")
else:
    print("\n[Step 7] primary or secondary matrix not available — skipping paired heatmap.")
