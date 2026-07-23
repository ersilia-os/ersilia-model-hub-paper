"""Convert EU OpenScreen H5 predictions to CSV and compute AUROCs against pathogen tasks.

H5 files in config/eu-openscreen_preds_h5/ are converted to flat CSVs in
data/processed/xx_euopenscreen_preds/. For each pathogen, only the primary assay
(as defined in data/raw/euopenscreen_data/primary_assays_manual.csv) is used.

Also produces a cross-pathogen heatmap: consensus_score AUROC of every pathogen model
against every available EU OpenScreen task, using only deduplicated compounds (training-
set molecules excluded per row pathogen).

Requires
--------
    config/eu-openscreen_preds_h5/*.h5
    config/pathogens_of_interest.csv
    data/raw/euopenscreen_data/primary_assays_manual.csv
    data/raw/euopenscreen_data/02_binarised_assays/{assay_eos_id}.csv
    data/raw/euopenscreen_data/02_merged/02_{code}.csv  (for InChIKey enrichment)
    data/raw/euopenscreen_data/06_subset_data/exclusivity/{code}_{exclusive,nonexclusive}.csv
    data/raw/euopenscreen_data/06_subset_data/secondary/{code}_secondary.csv

Outputs
-------
    data/processed/xx_euopenscreen_preds/{eosid}.csv        (one per model)
    output/xx_euopenscreen_preds/leakage_report.csv
    output/xx_euopenscreen_preds/auroc_scores.csv
    output/xx_euopenscreen_preds/auroc_dotplot.png
    output/xx_euopenscreen_preds/roc_curves.png
    output/xx_euopenscreen_preds/auroc_scores_deduplicated.csv
    output/xx_euopenscreen_preds/auroc_dotplot_deduplicated.png
    output/xx_euopenscreen_preds/roc_curves_deduplicated.png
    output/xx_euopenscreen_preds/auroc_heatmap_deduplicated.csv
    output/xx_euopenscreen_preds/auroc_heatmap_deduplicated.png
    output/xx_euopenscreen_preds/auroc_scores_{exclusive,nonexclusive,secondary}.csv
    output/xx_euopenscreen_preds/auroc_dotplot_{exclusive,nonexclusive,secondary}.png
"""

import os
import sys

import math

import h5py
import pandas as pd
import stylia
from sklearn.metrics import roc_auc_score

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))

from plotting_utils import abbrev, plot_auroc_dotplot, plot_auroc_heatmap, plot_roc_single

raw_h5_dir = os.path.join(root, "..", "config", "eu-openscreen_preds_h5")
processed_dir = os.path.join(root, "..", "data", "processed", "xx_euopenscreen_preds")
tasks_dir = os.path.join(root, "..", "data", "raw", "euopenscreen_data")
binarised_dir = os.path.join(tasks_dir, "02_binarised_assays")
merged_dir = os.path.join(tasks_dir, "02_merged")
subset_dir = os.path.join(tasks_dir, "06_subset_data")
models_dir = os.path.join(root, "..", "..", "chembl-antimicrobial-models")
output_dir = os.path.join(root, "..", "output", "xx_euopenscreen_preds")

os.makedirs(processed_dir, exist_ok=True)
os.makedirs(output_dir, exist_ok=True)

pathogens = pd.read_csv(
    os.path.join(root, "..", "config", "pathogens_of_interest.csv")
)

primary_df = pd.read_csv(os.path.join(tasks_dir, "primary_assays_manual.csv"))
code_to_assay = dict(zip(primary_df["pathogen_code"], primary_df["assay_eos_id"]))


# --- Step 1: Convert H5 to CSV ---

def h5_to_dataframe(h5_path):
    with h5py.File(h5_path, "r") as f:
        keys = [k.decode() for k in f["Keys"][:]]
        smiles = [s.decode() for s in f["Inputs"][:]]
        features = [ft.decode() for ft in f["Features"][:]]
        values = f["Values"][:]
    df = pd.DataFrame(values, columns=features)
    df.insert(0, "smiles", smiles)
    df.insert(0, "key", keys)
    return df


for h5_file in sorted(os.listdir(raw_h5_dir)):
    if not h5_file.endswith(".h5"):
        continue
    eosid = h5_file.replace("_v1.h5", "")
    out_csv = os.path.join(processed_dir, f"{eosid}.csv")
    if os.path.exists(out_csv):
        print(f"[SKIP] {eosid}: already converted")
        continue
    df = h5_to_dataframe(os.path.join(raw_h5_dir, h5_file))
    df.to_csv(out_csv, index=False)
    print(f"[OK] {eosid}: {len(df)} rows → {out_csv}")


# --- Step 2: Leakage check + AUROCs ---

records = []
roc_data = {}  # {pathogen: (y_true, y_pred_consensus)}
records_dedup = []
roc_data_dedup = {}
leakage_records = []
pathogen_train_keys = {}  # {code: set of InChIKeys in training set}

for _, row in pathogens.iterrows():
    eosid = row["eosid"]
    code = row["code"]
    pathogen = row["pathogen"]

    assay_id = code_to_assay.get(code)
    if assay_id is None:
        print(f"[SKIP] {pathogen}: no primary assay defined")
        continue
    task_csv = os.path.join(binarised_dir, f"{assay_id}.csv")
    pred_csv = os.path.join(processed_dir, f"{eosid}.csv")

    if not os.path.exists(task_csv):
        print(f"[SKIP] {pathogen}: task file not found at {task_csv}")
        continue
    if not os.path.exists(pred_csv):
        print(f"[SKIP] {pathogen}: predictions not found for {eosid}")
        continue

    task_df = pd.read_csv(task_csv)  # smiles, bin
    mpath = os.path.join(merged_dir, f"02_{code}.csv")
    if os.path.exists(mpath):
        inchi_df = pd.read_csv(mpath, usecols=["smiles", "inchikey"])
        task_df = task_df.merge(inchi_df, on="smiles", how="left")
    pred_df = pd.read_csv(pred_csv)
    feature_cols = [c for c in pred_df.columns if c not in ("key", "smiles")]

    merged = task_df.merge(pred_df, on="smiles", how="inner")
    if len(merged) == 0:
        print(f"[WARN] {pathogen}: no SMILES overlap between predictions and task data")
        continue

    # Drop inconclusive (-1) and undefined (NaN) labels; keep only 0/1
    merged = merged[merged["bin"].isin([0, 1])]
    if len(merged) == 0:
        print(f"[WARN] {pathogen}: no conclusive labels after filtering")
        continue

    n_pos = int(merged["bin"].sum())
    n_neg = int((merged["bin"] == 0).sum())
    print(f"[{pathogen}] {len(merged)} conclusive molecules — {n_pos} active, {n_neg} inactive")

    # Leakage check: compare evaluation InChIKeys against training data
    overlap_keys = set()
    train_csv = os.path.join(
        models_dir, "output", "17_quality_checks", code, "all_smiles_decoys.csv"
    )
    if os.path.exists(train_csv):
        train_df = pd.read_csv(
            train_csv, usecols=["inchikey", "n_active", "n_inactive", "n_decoy"]
        )
        eval_keys = set(task_df[task_df["bin"].isin([0, 1])]["inchikey"])
        train_keys = set(train_df["inchikey"])
        pathogen_train_keys[code] = train_keys
        overlap_keys = eval_keys & train_keys

        n_ov = len(overlap_keys)
        if n_ov > 0:
            ov_eval = task_df[task_df["inchikey"].isin(overlap_keys)]
            ov_train = train_df[train_df["inchikey"].isin(overlap_keys)]
            n_ov_active_eval = int((ov_eval["bin"] == 1).sum())
            n_ov_inactive_eval = int((ov_eval["bin"] == 0).sum())
            n_ov_active_train = int((ov_train["n_active"] > 0).sum())
            n_ov_inactive_train = int((ov_train["n_inactive"] > 0).sum())
            n_ov_decoy_train = int((ov_train["n_decoy"] > 0).sum())
            print(
                f"[WARN] {pathogen}: {n_ov} overlapping molecules "
                f"(eval active={n_ov_active_eval}, eval inactive={n_ov_inactive_eval} | "
                f"train active={n_ov_active_train}, train inactive={n_ov_inactive_train}, "
                f"train decoy={n_ov_decoy_train})"
            )
        else:
            n_ov_active_eval = n_ov_inactive_eval = 0
            n_ov_active_train = n_ov_inactive_train = n_ov_decoy_train = 0

        leakage_records.append({
            "pathogen": pathogen,
            "code": code,
            "n_train": len(train_df),
            "n_eval_conclusive": len(eval_keys),
            "n_overlap": n_ov,
            "n_overlap_active_eval": n_ov_active_eval,
            "n_overlap_inactive_eval": n_ov_inactive_eval,
            "n_overlap_active_train": n_ov_active_train,
            "n_overlap_inactive_train": n_ov_inactive_train,
            "n_overlap_decoy_train": n_ov_decoy_train,
        })
    else:
        print(f"[SKIP leakage] {pathogen}: training data not found at {train_csv}")

    for feat in feature_cols:
        auc = roc_auc_score(merged["bin"], merged[feat])
        records.append(
            {
                "pathogen": pathogen,
                "code": code,
                "eosid": eosid,
                "feature": feat,
                "auroc": round(auc, 4),
            }
        )

    if "consensus_score" in feature_cols:
        roc_data[pathogen] = (
            merged["bin"].values,
            merged["consensus_score"].values,
        )

    # Deduplicated: remove molecules seen in training
    merged_dedup = merged[~merged["inchikey"].isin(overlap_keys)]
    if merged_dedup["bin"].nunique() < 2:
        print(f"[SKIP dedup] {pathogen}: only one class remaining after deduplication")
    else:
        n_pos_d = int(merged_dedup["bin"].sum())
        n_neg_d = int((merged_dedup["bin"] == 0).sum())
        print(f"[{pathogen} dedup] {len(merged_dedup)} molecules — {n_pos_d} active, {n_neg_d} inactive")
        for feat in feature_cols:
            auc_d = roc_auc_score(merged_dedup["bin"], merged_dedup[feat])
            records_dedup.append({
                "pathogen": pathogen,
                "code": code,
                "eosid": eosid,
                "feature": feat,
                "auroc": round(auc_d, 4),
            })
        if "consensus_score" in feature_cols:
            roc_data_dedup[pathogen] = (
                merged_dedup["bin"].values,
                merged_dedup["consensus_score"].values,
            )

if leakage_records:
    leakage_df = pd.DataFrame(leakage_records)
    leakage_df.to_csv(os.path.join(output_dir, "leakage_report.csv"), index=False)
    print("\n--- Leakage report ---")
    print(leakage_df.to_string(index=False))

if not records:
    print(
        "No AUROC results — task data not yet available. "
        "Run again once eu-openscreen-antimicrobial-tasks is downloaded."
    )
    sys.exit(0)

auroc_df = pd.DataFrame(records)
auroc_df.to_csv(os.path.join(output_dir, "auroc_scores.csv"), index=False)
print(auroc_df.to_string(index=False))


# --- Step 3: Plots ---

# Format: print | Style: article
stylia.set_format("print")
stylia.set_style("article")

# 3a: dot plot — all assays, all pathogens
fig, axs = stylia.create_figure(1, 1)
plot_auroc_dotplot(axs.next(), auroc_df)
stylia.save_figure(os.path.join(output_dir, "auroc_dotplot.png"))

# 3b: ROC curves — consensus_score, one panel per pathogen
n = len(roc_data)
ncols = min(n, 4)
nrows = math.ceil(n / ncols)

fig, axs = stylia.create_figure(nrows, ncols)
for pathogen, (y_true, y_pred) in sorted(roc_data.items()):
    plot_roc_single(axs.next(), y_true, y_pred, title=abbrev(pathogen))
stylia.save_figure(os.path.join(output_dir, "roc_curves.png"))


# --- Step 4: Deduplicated outputs ---

if records_dedup:
    auroc_dedup_df = pd.DataFrame(records_dedup)
    auroc_dedup_df.to_csv(
        os.path.join(output_dir, "auroc_scores_deduplicated.csv"), index=False
    )

    fig, axs = stylia.create_figure(1, 1)
    plot_auroc_dotplot(axs.next(), auroc_dedup_df)
    stylia.save_figure(os.path.join(output_dir, "auroc_dotplot_deduplicated.png"))

    if roc_data_dedup:
        n_d = len(roc_data_dedup)
        ncols_d = min(n_d, 4)
        nrows_d = math.ceil(n_d / ncols_d)
        fig, axs = stylia.create_figure(nrows_d, ncols_d)
        for pathogen, (y_true, y_pred) in sorted(roc_data_dedup.items()):
            plot_roc_single(axs.next(), y_true, y_pred, title=abbrev(pathogen))
        stylia.save_figure(os.path.join(output_dir, "roc_curves_deduplicated.png"))


# --- Step 5: Cross-pathogen AUROC heatmap (deduplicated) ---

# Collect primary assay task files, keyed by abbreviated pathogen name
task_files = {}
for _, mrow in primary_df.iterrows():
    tpath = os.path.join(binarised_dir, f"{mrow['assay_eos_id']}.csv")
    if os.path.exists(tpath):
        task_files[abbrev(mrow["pathogen"])] = {
            "path": tpath,
            "code": mrow["pathogen_code"],
        }

heatmap_records = []
for _, prow in pathogens.iterrows():
    eosid = prow["eosid"]
    pcode = prow["code"]
    pathogen = prow["pathogen"]

    pred_csv = os.path.join(processed_dir, f"{eosid}.csv")
    train_keys = pathogen_train_keys.get(pcode, set())
    row_record = {"pathogen": pathogen}

    if not os.path.exists(pred_csv):
        print(f"[SKIP heatmap] {pathogen}: predictions not found for {eosid}")
        for tname in task_files:
            row_record[tname] = float("nan")
        heatmap_records.append(row_record)
        continue

    pred_df = pd.read_csv(pred_csv)
    if "consensus_score" not in pred_df.columns:
        print(f"[SKIP heatmap] {pathogen}: no consensus_score in {eosid}")
        for tname in task_files:
            row_record[tname] = float("nan")
        heatmap_records.append(row_record)
        continue
    pred_df = pred_df[["smiles", "consensus_score"]]

    for tname, tinfo in task_files.items():
        task_df_h = pd.read_csv(tinfo["path"])  # smiles, bin
        mpath_h = os.path.join(merged_dir, f"02_{tinfo['code']}.csv")
        if os.path.exists(mpath_h):
            inchi_h = pd.read_csv(mpath_h, usecols=["smiles", "inchikey"])
            task_df_h = task_df_h.merge(inchi_h, on="smiles", how="left")
        merged_h = task_df_h.merge(pred_df, on="smiles", how="inner")
        merged_h = merged_h[merged_h["bin"].isin([0, 1])]
        if train_keys and "inchikey" in merged_h.columns:
            merged_h = merged_h[~merged_h["inchikey"].isin(train_keys)]
        if merged_h["bin"].nunique() < 2:
            row_record[tname] = float("nan")
        else:
            row_record[tname] = round(
                roc_auc_score(merged_h["bin"], merged_h["consensus_score"]), 4
            )

    heatmap_records.append(row_record)

heatmap_df = pd.DataFrame(heatmap_records).set_index("pathogen")
heatmap_df.to_csv(os.path.join(output_dir, "auroc_heatmap_deduplicated.csv"))
print("\n--- AUROC heatmap (deduplicated) ---")
print(heatmap_df.to_string())

fig, axs = stylia.create_figure(1, 1, width=0.5, height=0.5)
plot_auroc_heatmap(
    axs.next(), heatmap_df,
    title="AUROC — consensus score vs EU OpenScreen tasks (deduplicated)"
)
stylia.save_figure(os.path.join(output_dir, "auroc_heatmap_deduplicated.png"))


# --- Step 6: Subset analyses (exclusive, nonexclusive, secondary) ---

def _load_subset_task(primary_csv, merged_path, mode, subset_dir, code):
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
    primary = pd.read_csv(primary_csv)  # smiles, bin
    if os.path.exists(merged_path):
        inchi_df = pd.read_csv(merged_path, usecols=["smiles", "inchikey"])
        primary = primary.merge(inchi_df, on="smiles", how="left")
    inactives = primary[primary["bin"] == 0].copy()
    return pd.concat([actives, inactives], ignore_index=True)


_SUBSET_PLOT_TITLES = {
    "exclusive": "AUROC — exclusive actives vs primary inactives",
    "nonexclusive": "AUROC — non-exclusive actives vs primary inactives",
    "secondary": "AUROC — secondary assay",
}

for mode in ("exclusive", "nonexclusive", "secondary"):
    subset_records = []

    for _, row in pathogens.iterrows():
        eosid = row["eosid"]
        code = row["code"]
        pathogen = row["pathogen"]

        assay_id = code_to_assay.get(code)
        if assay_id is None:
            continue

        task_df_s = _load_subset_task(
            primary_csv=os.path.join(binarised_dir, f"{assay_id}.csv"),
            merged_path=os.path.join(merged_dir, f"02_{code}.csv"),
            mode=mode,
            subset_dir=subset_dir,
            code=code,
        )
        if task_df_s is None:
            continue

        pred_csv_s = os.path.join(processed_dir, f"{eosid}.csv")
        if not os.path.exists(pred_csv_s):
            continue

        pred_df_s = pd.read_csv(pred_csv_s)
        feature_cols_s = [c for c in pred_df_s.columns if c not in ("key", "smiles")]

        merged_s = task_df_s.merge(pred_df_s, on="smiles", how="inner")
        merged_s = merged_s[merged_s["bin"].isin([0, 1])]

        train_keys_s = pathogen_train_keys.get(code, set())
        if train_keys_s and "inchikey" in merged_s.columns:
            merged_s = merged_s[~merged_s["inchikey"].isin(train_keys_s)]

        if merged_s["bin"].nunique() < 2:
            print(f"[SKIP {mode}] {pathogen}: fewer than 2 classes")
            continue

        n_pos_s = int(merged_s["bin"].sum())
        n_neg_s = int((merged_s["bin"] == 0).sum())
        print(f"[{pathogen} {mode}] {len(merged_s)} molecules — {n_pos_s} active, {n_neg_s} inactive")

        for feat in feature_cols_s:
            auc_s = roc_auc_score(merged_s["bin"], merged_s[feat])
            subset_records.append({
                "pathogen": pathogen,
                "code": code,
                "eosid": eosid,
                "feature": feat,
                "auroc": round(auc_s, 4),
            })

    if subset_records:
        df_s = pd.DataFrame(subset_records)
        df_s.to_csv(os.path.join(output_dir, f"auroc_scores_{mode}.csv"), index=False)
        fig, axs = stylia.create_figure(1, 1)
        plot_auroc_dotplot(axs.next(), df_s, title=_SUBSET_PLOT_TITLES[mode])
        stylia.save_figure(os.path.join(output_dir, f"auroc_dotplot_{mode}.png"))
        print(f"[{mode}] {len(df_s)} AUROC records → auroc_scores_{mode}.csv")
    else:
        print(f"[{mode}] No results — subset data may not be downloaded yet.")


# --- Step 7: Comparison heatmaps (consensus_score only) ---

def _consensus_series(csv_path):
    """Extract consensus_score AUROC per pathogen from an auroc_scores CSV."""
    if not os.path.exists(csv_path):
        return None
    df = pd.read_csv(csv_path)
    return df[df["feature"] == "consensus_score"].set_index("pathogen")["auroc"]


cs_primary = _consensus_series(os.path.join(output_dir, "auroc_scores.csv"))
cs_exclusive = _consensus_series(os.path.join(output_dir, "auroc_scores_exclusive.csv"))
cs_nonexclusive = _consensus_series(os.path.join(output_dir, "auroc_scores_nonexclusive.csv"))
cs_secondary = _consensus_series(os.path.join(output_dir, "auroc_scores_secondary.csv"))

# Sort pathogens by primary AUROC so both figures share the same y-axis order
if cs_primary is not None:
    pathogen_order = cs_primary.sort_values(ascending=False).index.tolist()
else:
    pathogen_order = None


def _make_panel_df(series, col_name, order):
    s = series.reindex(order) if order is not None else series
    return pd.DataFrame({col_name: s})


# Figure 1: exclusive vs non-exclusive
if cs_exclusive is not None and cs_nonexclusive is not None:
    order = pathogen_order or cs_exclusive.sort_values(ascending=False).index.tolist()
    fig, axs = stylia.create_figure(1, 2)
    plot_auroc_heatmap(
        axs.next(), _make_panel_df(cs_exclusive, "Exclusive", order),
        title="Exclusive actives"
    )
    plot_auroc_heatmap(
        axs.next(), _make_panel_df(cs_nonexclusive, "Non-exclusive", order),
        title="Non-exclusive actives"
    )
    stylia.save_figure(os.path.join(output_dir, "auroc_comparison_exclusivity.png"))
    print("Step 7: saved auroc_comparison_exclusivity.png")
else:
    print("Step 7: exclusive/nonexclusive scores not found — skipping exclusivity comparison.")

# Figure 2: primary vs secondary
if cs_primary is not None and cs_secondary is not None:
    order = pathogen_order or cs_primary.sort_values(ascending=False).index.tolist()
    fig, axs = stylia.create_figure(1, 2)
    plot_auroc_heatmap(
        axs.next(), _make_panel_df(cs_primary, "Primary", order),
        title="Primary assay"
    )
    plot_auroc_heatmap(
        axs.next(), _make_panel_df(cs_secondary, "Secondary", order),
        title="Secondary assay"
    )
    stylia.save_figure(os.path.join(output_dir, "auroc_comparison_primary_secondary.png"))
    print("Step 7: saved auroc_comparison_primary_secondary.png")
else:
    print("Step 7: primary/secondary scores not found — skipping primary vs secondary comparison.")
