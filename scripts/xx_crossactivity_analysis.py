"""Assess cross-pathogen discriminatory power of ChEMBL and CoAdd models.

Three analyses address whether pathogen-specific models distinguish pathogen-
specific activity or simply learn general antimicrobial features:

  1. Active compound overlap — Jaccard similarity between active sets of the 7
     EU OpenScreen tasks. High cross-pathogen overlap explains high cross-model
     AUROCs without implying model failure.

  2. Specificity index (ChEMBL) — for each of the 7 ChEMBL pathogen models with
     a matching EU OpenScreen task, computes: same-pathogen AUROC minus mean
     cross-pathogen AUROC. Positive values indicate the model is better on its
     own pathogen.

  3. Pan-active vs specific-active AUROC — for each task, splits its actives
     into compounds also active in ≥1 other task (pan-actives) and compounds
     unique to that task (specific-actives). AUROCs computed for both subsets
     using the matched ChEMBL consensus_score and the best-matching CoAdd
     endpoint. A drop toward 0.5 for specific-actives signals that models
     cannot distinguish truly pathogen-specific compounds.

Requires
--------
    data/raw/euopenscreen_data/{code}.csv
    config/pathogens_of_interest.csv
    data/processed/xx_euopenscreen_preds/{eosid}.csv   (ChEMBL predictions)
    data/processed/xx_euopenscreen_preds/eos3dys.csv   (CoAdd predictions)
    output/xx_euopenscreen_preds/auroc_heatmap_deduplicated.csv
    output/xx_coadd_benchmark/auroc_matrix.csv

Outputs
-------
    output/xx_crossactivity_analysis/active_overlap_counts.csv
    output/xx_crossactivity_analysis/active_overlap_jaccard.csv
    output/xx_crossactivity_analysis/active_overlap_heatmap.png
    output/xx_crossactivity_analysis/specificity_index.csv
    output/xx_crossactivity_analysis/specificity_index.png
    output/xx_crossactivity_analysis/panactive_auroc.csv
    output/xx_crossactivity_analysis/panactive_auroc.png
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
from plotting_utils import abbrev, plot_auroc_heatmap, plot_specificity_bars

tasks_dir = os.path.join(root, "..", "data", "raw", "euopenscreen_data")
preds_dir = os.path.join(root, "..", "data", "processed", "xx_euopenscreen_preds")
config_dir = os.path.join(root, "..", "config")
chembl_heatmap_csv = os.path.join(
    root, "..", "output", "xx_euopenscreen_preds", "auroc_heatmap_deduplicated.csv"
)
coadd_matrix_csv = os.path.join(
    root, "..", "output", "xx_coadd_benchmark", "auroc_matrix.csv"
)
output_dir = os.path.join(root, "..", "output", "xx_crossactivity_analysis")

os.makedirs(output_dir, exist_ok=True)

# Format: print | Style: article
stylia.set_format("print")
stylia.set_style("article")


# =============================================================================
# Load shared inputs
# =============================================================================

pathogens = pd.read_csv(os.path.join(config_dir, "pathogens_of_interest.csv"))
pathogen_to_code = dict(zip(pathogens["pathogen"], pathogens["code"]))
code_to_eosid = dict(zip(pathogens["code"], pathogens["eosid"]))

task_files = {
    fname[:-4]: os.path.join(tasks_dir, fname)
    for fname in sorted(os.listdir(tasks_dir))
    if fname.endswith(".csv") and os.path.isfile(os.path.join(tasks_dir, fname))
}
task_codes = sorted(task_files.keys())

if not task_files:
    print("No EU OpenScreen task files found. Run 00_download_data.py first.")
    sys.exit(0)

# Load all task dataframes once
task_dfs = {code: pd.read_csv(path) for code, path in task_files.items()}

# Active SMILES per task (bin==1)
active_smiles = {
    code: set(df[df["bin"] == 1]["smiles"]) for code, df in task_dfs.items()
}

print(f"Loaded {len(task_codes)} task files: {task_codes}")
for code, actives in active_smiles.items():
    print(f"  {code}: {len(actives)} actives")


# =============================================================================
# Step 1 — Active compound overlap matrix
# =============================================================================

print("\n--- Step 1: Active compound overlap ---")

counts = pd.DataFrame(index=task_codes, columns=task_codes, dtype=float)
jaccard = pd.DataFrame(index=task_codes, columns=task_codes, dtype=float)

for a in task_codes:
    for b in task_codes:
        inter = len(active_smiles[a] & active_smiles[b])
        union = len(active_smiles[a] | active_smiles[b])
        counts.loc[a, b] = inter
        jaccard.loc[a, b] = round(inter / union, 4) if union > 0 else 0.0

counts.to_csv(os.path.join(output_dir, "active_overlap_counts.csv"))
jaccard.to_csv(os.path.join(output_dir, "active_overlap_jaccard.csv"))

print("Active overlap (Jaccard):")
print(jaccard.to_string())

fig, axs = stylia.create_figure(1, 1, width=0.5, height=0.5)
plot_auroc_heatmap(
    axs.next(),
    jaccard.astype(float),
    title="Active compound overlap — Jaccard similarity between tasks",
)
stylia.save_figure(os.path.join(output_dir, "active_overlap_heatmap.png"))


# =============================================================================
# Step 2 — Specificity index (ChEMBL models)
# =============================================================================

print("\n--- Step 2: ChEMBL specificity index ---")

if not os.path.exists(chembl_heatmap_csv):
    print(f"[SKIP] ChEMBL heatmap not found at {chembl_heatmap_csv}")
    specificity_df = pd.DataFrame()
else:
    heatmap_df = pd.read_csv(chembl_heatmap_csv, index_col=0)
    task_col_set = set(heatmap_df.columns)

    records = []
    for pathogen in heatmap_df.index:
        code = pathogen_to_code.get(pathogen)
        row = heatmap_df.loc[pathogen]
        valid_cols = [c for c in heatmap_df.columns if not pd.isna(row[c])]

        if code and code in task_col_set:
            same = row[code] if not pd.isna(row[code]) else float("nan")
            cross_vals = [row[c] for c in valid_cols if c != code]
            mean_cross = float(np.mean(cross_vals)) if cross_vals else float("nan")
            spec_idx = round(float(same) - mean_cross, 4) if not np.isnan(same) and not np.isnan(mean_cross) else float("nan")
        else:
            same = float("nan")
            mean_cross = float(np.mean([row[c] for c in valid_cols])) if valid_cols else float("nan")
            spec_idx = float("nan")

        records.append({
            "pathogen": pathogen,
            "code": code,
            "same_pathogen_auroc": round(float(same), 4) if not np.isnan(same) else float("nan"),
            "mean_cross_auroc": round(mean_cross, 4) if not np.isnan(mean_cross) else float("nan"),
            "specificity_index": spec_idx,
        })

    specificity_df = pd.DataFrame(records)
    specificity_df.to_csv(os.path.join(output_dir, "specificity_index.csv"), index=False)
    print(specificity_df.to_string(index=False))

    fig, axs = stylia.create_figure(1, 1)
    plot_specificity_bars(
        axs.next(),
        specificity_df,
        title="ChEMBL model specificity — same vs cross-pathogen AUROC",
    )
    stylia.save_figure(os.path.join(output_dir, "specificity_index.png"))


# =============================================================================
# Step 3 — Pan-active vs specific-active AUROC
# =============================================================================

print("\n--- Step 3: Pan-active vs specific-active AUROC ---")

# Label each active SMILES as pan-active (in ≥2 tasks) or specific (in exactly 1)
all_active_smiles = {}
for code in task_codes:
    for smi in active_smiles[code]:
        all_active_smiles.setdefault(smi, set()).add(code)

pan_active_smiles = {smi for smi, tasks in all_active_smiles.items() if len(tasks) >= 2}
specific_smiles = {smi for smi, tasks in all_active_smiles.items() if len(tasks) == 1}

n_pan = len(pan_active_smiles)
n_spec = len(specific_smiles)
print(f"Pan-active compounds (active in ≥2 tasks): {n_pan}")
print(f"Specific-active compounds (active in exactly 1 task): {n_spec}")

# Load CoAdd predictions and select best-matching endpoint per task
coadd_pred = pd.read_csv(os.path.join(preds_dir, f"{COADD_MODEL_ID}.csv"))
endpoint_cols = [c for c in coadd_pred.columns if c not in ("key", "smiles")]

best_coadd_ep = {}
if os.path.exists(coadd_matrix_csv):
    coadd_matrix = pd.read_csv(coadd_matrix_csv, index_col=0)
    for tcode in task_codes:
        matching = [ep for ep in endpoint_cols if ep.startswith(tcode + "_")]
        if not matching:
            best_coadd_ep[tcode] = None
            print(f"  [CoAdd] No endpoint found for {tcode} — will be NaN")
            continue
        # Best = highest AUROC on that specific task column
        if tcode in coadd_matrix.columns:
            ep_aurocs = coadd_matrix.loc[matching, tcode]
            best_ep = ep_aurocs.dropna().idxmax() if ep_aurocs.dropna().size > 0 else matching[0]
        else:
            best_ep = matching[0]
        best_coadd_ep[tcode] = best_ep
        print(f"  [CoAdd] Best endpoint for {tcode}: {best_ep}")
else:
    print(f"[WARN] CoAdd matrix not found at {coadd_matrix_csv} — CoAdd results will be NaN")
    best_coadd_ep = {tcode: None for tcode in task_codes}

MIN_SPECIFIC_ACTIVES = 5


def _compute_auroc(positives_smiles, neg_smiles, pred_df, score_col):
    """Merge positives + negatives with predictions and compute AUROC."""
    pos_set = set(positives_smiles)
    neg_set = set(neg_smiles)
    subset = pred_df[pred_df["smiles"].isin(pos_set | neg_set)][["smiles", score_col]]
    subset = subset.drop_duplicates("smiles")
    subset = subset.copy()
    subset["label"] = subset["smiles"].map(lambda s: 1 if s in pos_set else 0)
    subset = subset.dropna(subset=[score_col])
    if subset["label"].nunique() < 2:
        return float("nan")
    return round(roc_auc_score(subset["label"], subset[score_col]), 4)


pan_records = []
for tcode in task_codes:
    eosid = code_to_eosid.get(tcode)
    task_df = task_dfs[tcode]
    inactives = task_df[task_df["bin"] == 0]["smiles"].tolist()

    # Pan-actives for this task: active here AND in ≥1 other task
    task_pan = active_smiles[tcode] & pan_active_smiles
    # Specific-actives for this task: active here AND in no other task
    task_spec = active_smiles[tcode] & specific_smiles

    n_pan_task = len(task_pan)
    n_spec_task = len(task_spec)
    print(f"\n{tcode}: {n_pan_task} pan-actives, {n_spec_task} specific-actives")

    row = {
        "task": tcode,
        "n_panactive": n_pan_task,
        "n_specific": n_spec_task,
        "auroc_panactive_chembl": float("nan"),
        "auroc_specific_chembl": float("nan"),
        "auroc_panactive_coadd": float("nan"),
        "auroc_specific_coadd": float("nan"),
    }

    # --- ChEMBL ---
    if eosid:
        chembl_csv = os.path.join(preds_dir, f"{eosid}.csv")
        if os.path.exists(chembl_csv):
            chembl_pred = pd.read_csv(chembl_csv)[["smiles", "consensus_score"]]
            if n_pan_task >= 1:
                row["auroc_panactive_chembl"] = _compute_auroc(
                    task_pan, inactives, chembl_pred, "consensus_score"
                )
            if n_spec_task >= MIN_SPECIFIC_ACTIVES:
                row["auroc_specific_chembl"] = _compute_auroc(
                    task_spec, inactives, chembl_pred, "consensus_score"
                )
            else:
                print(f"  [ChEMBL] {tcode}: only {n_spec_task} specific-actives (< {MIN_SPECIFIC_ACTIVES}) — AUROC unreliable, set to NaN")
        else:
            print(f"  [ChEMBL] predictions not found for {eosid}")

    # --- CoAdd ---
    best_ep = best_coadd_ep.get(tcode)
    if best_ep:
        coadd_col = coadd_pred[["smiles", best_ep]]
        if n_pan_task >= 1:
            row["auroc_panactive_coadd"] = _compute_auroc(
                task_pan, inactives, coadd_col, best_ep
            )
        if n_spec_task >= MIN_SPECIFIC_ACTIVES:
            row["auroc_specific_coadd"] = _compute_auroc(
                task_spec, inactives, coadd_col, best_ep
            )
        else:
            print(f"  [CoAdd]  {tcode}: only {n_spec_task} specific-actives (< {MIN_SPECIFIC_ACTIVES}) — AUROC unreliable, set to NaN")

    print(f"  ChEMBL — pan: {row['auroc_panactive_chembl']}, specific: {row['auroc_specific_chembl']}")
    print(f"  CoAdd  — pan: {row['auroc_panactive_coadd']}, specific: {row['auroc_specific_coadd']}")
    pan_records.append(row)

pan_df = pd.DataFrame(pan_records)
pan_df.to_csv(os.path.join(output_dir, "panactive_auroc.csv"), index=False)
print("\n--- Pan-active vs specific-active AUROC ---")
print(pan_df.to_string(index=False))


# =============================================================================
# Pan-active AUROC plot: paired scatter, one panel per model type
# =============================================================================

nc = stylia.NamedColors()

def _plot_panactive_panel(ax, df, pan_col, spec_col, title):
    x = range(len(df))
    x_pan = [i - 0.1 for i in x]
    x_spec = [i + 0.1 for i in x]

    pan_vals = df[pan_col].values
    spec_vals = df[spec_col].values

    # Lines connecting pan/specific pairs
    for i in range(len(df)):
        if not np.isnan(pan_vals[i]) and not np.isnan(spec_vals[i]):
            ax.plot([x_pan[i], x_spec[i]], [pan_vals[i], spec_vals[i]],
                    color=nc.silver, linewidth=0.8)

    ax.scatter(x_pan, pan_vals, color=nc.turquoise, zorder=3, label="Pan-active")
    ax.scatter(x_spec, spec_vals, color=nc.crimson, zorder=3, label="Specific-active")
    ax.axhline(0.5, color=nc.silver, linestyle="--")
    ax.set_xticks(list(x))
    ax.set_xticklabels(df["task"].tolist(), rotation=45, ha="right")
    ax.set_ylim(0, 1)
    ax.legend(loc="lower right")
    stylia.label(ax, xlabel=" ", ylabel="AUROC", title=title)


fig, axs = stylia.create_figure(1, 2)
_plot_panactive_panel(
    axs.next(), pan_df,
    "auroc_panactive_chembl", "auroc_specific_chembl",
    "ChEMBL — pan-active vs specific-active",
)
_plot_panactive_panel(
    axs.next(), pan_df,
    "auroc_panactive_coadd", "auroc_specific_coadd",
    "CoAdd — pan-active vs specific-active",
)
stylia.save_figure(os.path.join(output_dir, "panactive_auroc.png"))

print("\nDone. Outputs in", output_dir)
