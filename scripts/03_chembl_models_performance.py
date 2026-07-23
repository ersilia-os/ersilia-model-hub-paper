"""Visualise ChEMBL antimicrobial model performance per pathogen.

For each pathogen subdirectory in data/raw/chembl_model_reports/, loads all
per-model reports (5-fold cross-validation), computes mean ± std AUROC, and
produces two plots: one ROC curve per model in a subplot grid (all folds
concatenated, coloured by mean AUROC) and horizontal paired rank boxplots
(all folds pooled). Per-fold AUROC comes from each {name}.csv; the per-compound
probability and rank arrays come from the matching {name}_folds.json sidecar.

Requires
--------
    data/raw/chembl_model_reports/{pathogen}/{name}.csv
    data/raw/chembl_model_reports/{pathogen}/{name}_folds.json

Outputs
-------
    output/03_chembl_models_performance/{pathogen}_auroc_summary.csv
    output/03_chembl_models_performance/{pathogen}_roc_curves.png
    output/03_chembl_models_performance/{pathogen}_rank_boxplots.png
"""

import json
import math
import os
import sys

import numpy as np
import pandas as pd
import stylia

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))

from plotting_utils import plot_rank_boxplots, plot_roc_single

reports_dir = os.path.join(root, "..", "data", "raw", "chembl_model_reports")
output_dir = os.path.join(root, "..", "output", "03_chembl_models_performance")
os.makedirs(output_dir, exist_ok=True)

# Format: print | Style: article
stylia.set_format("print")
stylia.set_style("article")


def load_models(pathogen_dir):
    """Load all model reports in a pathogen directory.

    Per-fold summary metrics come from each {name}.csv; the per-compound
    probability and rank arrays come from the matching {name}_folds.json
    sidecar (keyed by fold, each with y_true / y_hat / y_rank).

    Returns a list of dicts, one per model:
        name, mean_auroc, std_auroc, y_true, y_pred,
        rank_actives, rank_inactives   (all folds pooled, out-of-fold)
    """
    models = []
    for fname in sorted(os.listdir(pathogen_dir)):
        if not fname.endswith(".csv"):
            continue
        fpath = os.path.join(pathogen_dir, fname)
        try:
            df = pd.read_csv(fpath)
        except Exception as e:
            print(f"  [SKIP] {fname}: {e}")
            continue

        required = {"model_name", "fold", "auroc"}
        if not required.issubset(df.columns):
            print(f"  [SKIP] {fname}: missing required columns")
            continue

        folds_path = os.path.join(pathogen_dir, f"{fname[:-4]}_folds.json")
        if not os.path.exists(folds_path):
            print(f"  [SKIP] {fname}: missing {os.path.basename(folds_path)}")
            continue
        with open(folds_path) as fh:
            folds = json.load(fh)

        model_name = df["model_name"].iloc[0]
        mean_auroc = df["auroc"].mean()
        std_auroc = df["auroc"].std()

        # Pool per-compound predictions across all folds; each compound is
        # held out (out-of-fold) exactly once.
        y_true, y_hat, y_rank = [], [], []
        for fd in folds.values():
            y_true.extend(fd["y_true"])
            y_hat.extend(fd["y_hat"])
            y_rank.extend(fd["y_rank"])
        y_true = np.asarray(y_true)
        y_hat = np.asarray(y_hat, dtype=float)
        y_rank = np.asarray(y_rank, dtype=float)

        rank_actives = y_rank[y_true == 1]
        rank_inactives = y_rank[y_true == 0]

        models.append({
            "name": model_name,
            "mean_auroc": mean_auroc,
            "std_auroc": std_auroc,
            "y_true": y_true.tolist(),
            "y_pred": y_hat.tolist(),
            "rank_actives": rank_actives,
            "rank_inactives": rank_inactives,
        })

    return models


# "10_reports" is the aggregated-summary subdir staged by 00_download_data.py,
# not a pathogen — exclude it from the per-pathogen iteration.
NON_PATHOGEN_DIRS = {"10_reports"}
pathogens = sorted(
    d for d in os.listdir(reports_dir)
    if os.path.isdir(os.path.join(reports_dir, d)) and d not in NON_PATHOGEN_DIRS
)

if not pathogens:
    print(f"No pathogen subdirectories found in {reports_dir}. Exiting.")
    sys.exit(0)

for pathogen in pathogens:
    pathogen_dir = os.path.join(reports_dir, pathogen)
    print(f"\n[{pathogen}] Loading models...")
    models = load_models(pathogen_dir)

    if not models:
        print(f"  No valid model CSVs found — skipping.")
        continue

    print(f"  {len(models)} models loaded.")

    # AUROC summary CSV
    summary_df = pd.DataFrame([
        {"model": m["name"], "mean_auroc": round(m["mean_auroc"], 4),
         "std_auroc": round(m["std_auroc"], 4)}
        for m in models
    ]).sort_values("mean_auroc", ascending=False)
    summary_path = os.path.join(output_dir, f"{pathogen}_auroc_summary.csv")
    summary_df.to_csv(summary_path, index=False)
    print(f"  -> {summary_path}")

    # ROC curves — one subplot per model, sorted descending by mean AUROC
    n_models = len(models)
    ncols = min(n_models, 4)
    nrows = math.ceil(n_models / ncols)
    cm = stylia.FadingColormap("cobalt")
    cm.fit(np.array([0.5, 1.0]))
    fig, axs = stylia.create_figure(nrows, ncols)
    for m in sorted(models, key=lambda x: x["mean_auroc"], reverse=True):
        color = cm.transform(np.array([m["mean_auroc"]]))[0]
        plot_roc_single(axs.next(), m["y_true"], m["y_pred"], title=m["name"], color=color)
    roc_path = os.path.join(output_dir, f"{pathogen}_roc_curves.png")
    stylia.save_figure(roc_path)
    print(f"  -> {roc_path}")

    # Rank boxplots — height scales with number of models
    bx_height = max(0.5, n_models * 0.07)
    fig, axs = stylia.create_figure(1, 1, height=bx_height)
    plot_rank_boxplots(axs.next(), models, title=pathogen)
    bx_path = os.path.join(output_dir, f"{pathogen}_rank_boxplots.png")
    stylia.save_figure(bx_path)
    print(f"  -> {bx_path}")

print("\nDone.")
