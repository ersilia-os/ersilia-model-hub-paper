"""Visualise ChEMBL antimicrobial model performance per pathogen.

For each pathogen subdirectory in data/raw/chembl_model_reports/, loads all
per-model reports (5-fold cross-validation), computes mean ± std AUROC, and
produces two plots: one ROC curve per model in a subplot grid (all folds
concatenated, coloured by mean AUROC) and horizontal paired rank boxplots
(all folds pooled). Per-fold AUROC comes from each {name}.csv; the per-compound
probability and rank arrays come from the matching {name}_folds.json sidecar.

The per-pathogen folders are the pipeline's **step 09** reports (196 trained models),
which is a superset of the 193 models retained by step 10: three models were trained
and reported but then discarded for mean AUROC < 0.7 (calbicans/588506,
hpylori/SP_catchall, pfalciparum/743093_merged2). All 196 are plotted; the summary
CSVs carry the step-10 verdict in `retained` / `discard_reason` so the distinction is
auditable. (The other three step-10 discards are `untrainable` and never produced a
step-09 report, so they cannot appear here at all.)

Requires
--------
    config/pathogens_of_interest.csv
    data/raw/chembl_model_reports/{pathogen}/{name}.csv
    data/raw/chembl_model_reports/{pathogen}/{name}_folds.json
    data/raw/chembl_model_reports/10_reports/10_reports.csv
    data/raw/chembl_model_reports/10_reports/10_discarded_models.csv

Outputs
-------
    output/03_chembl_models_performance/{pathogen}_auroc_summary.csv
    output/03_chembl_models_performance/individual_plots/png/{pathogen}_*.png
    output/03_chembl_models_performance/individual_plots/pdf/{pathogen}_*.pdf
    output/03_chembl_models_performance/individual_plots/figure_cells.json

The per-pathogen panels are intermediate results — 30 figures, too many for the paper.
They land in `individual_plots/`, leaving the top level of the output dir for the
condensed cross-pathogen figures.
"""

import json
import os
import sys

import numpy as np
import pandas as pd

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))

# Importing the plotting stack applies the publication presets (print/article).
from plots_chembl_performance import (
    save_activity_ratios_figure,
    save_consensus_figure,
    save_performance_figures,
    write_figure_cells,
)

config_dir = os.path.join(root, "..", "config")
reports_dir = os.path.join(root, "..", "data", "raw", "chembl_model_reports")
curation_dir = os.path.join(root, "..", "data", "raw", "chembl_curation")
output_dir = os.path.join(root, "..", "output", "03_chembl_models_performance")
# One panel per pathogen is an intermediate result, not a paper figure — they are
# kept apart from the condensed cross-pathogen figures that go in output_dir itself.
individual_dir = os.path.join(output_dir, "individual_plots")
os.makedirs(output_dir, exist_ok=True)
os.makedirs(individual_dir, exist_ok=True)


def load_models(pathogen_dir):
    """Load all model reports in a pathogen directory.

    Per-fold summary metrics come from each {name}.csv; the per-compound
    probability and rank arrays come from the matching {name}_folds.json
    sidecar (keyed by fold, each with y_true / y_hat / y_rank).

    Returns a list of dicts, one per model:
        name, stem, mean_auroc, std_auroc, y_true, y_pred,
        rank_actives, rank_inactives   (all folds pooled, out-of-fold)

    ``stem`` is the file stem, always a string, and is the join key against the
    step-10 tables. The in-file ``model_name`` column is *not* usable for that: 26
    models are named with digits only (e.g. 1242), so pandas types them as int64
    here and the join would silently miss them.
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
            "stem": fname[:-4],
            "mean_auroc": mean_auroc,
            "std_auroc": std_auroc,
            "y_true": y_true.tolist(),
            "y_pred": y_hat.tolist(),
            "rank_actives": rank_actives,
            "rank_inactives": rank_inactives,
        })

    return models


# Pathogen code -> full binomial, used for the panel titles ("calbicans" -> "C. albicans").
pathogen_names = (
    pd.read_csv(os.path.join(config_dir, "pathogens_of_interest.csv"))
    .set_index("code")["pathogen"]
    .to_dict()
)

# Step-10 verdict per model. `name` is forced to str on both sides because the
# purely-numeric model names (1242, 588506, ...) would otherwise be typed as int64
# in one table and str in the other, and the lookup would miss them.
step10_dir = os.path.join(reports_dir, "10_reports")
retained_keys = set(
    pd.read_csv(os.path.join(step10_dir, "10_reports.csv"), dtype={"name": str})
    .apply(lambda r: (r["pathogen"], r["name"]), axis=1)
)
discard_reasons = {
    (r["pathogen"], r["name"]): r["reason"]
    for _, r in pd.read_csv(
        os.path.join(step10_dir, "10_discarded_models.csv"), dtype={"name": str}
    ).iterrows()
}

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

# Unique cleaned molecules per pathogen (for the consensus-dot size); staged step-02 summary.
coverage_path = os.path.join(curation_dir, "general", "27_chembl_coverage.csv")
mol_counts = {}
if os.path.exists(coverage_path):
    cov = pd.read_csv(coverage_path)
    if "is_union" in cov.columns:
        cov = cov[cov["is_union"] == False]  # noqa: E712 — per-pathogen rows, not the union
    mol_counts = dict(zip(cov["pathogen"], cov["n_cleaned_inchikeys"]))
else:
    print(f"[WARN] {coverage_path} not found — consensus dots will use a default size")

footprints = {}
plotted_keys = set()
consensus_entries = []  # per-pathogen retained CV AUROCs for the condensed summary figure

for pathogen in pathogens:
    pathogen_dir = os.path.join(reports_dir, pathogen)
    print(f"\n[{pathogen}] Loading models...")
    models = load_models(pathogen_dir)

    if not models:
        print("  No valid model CSVs found — skipping.")
        continue

    print(f"  {len(models)} models loaded.")
    plotted_keys.update((pathogen, m["stem"]) for m in models)

    # AUROC summary CSV, carrying the step-10 keep/discard verdict
    summary_df = pd.DataFrame([
        {"model": m["name"], "mean_auroc": round(m["mean_auroc"], 4),
         "std_auroc": round(m["std_auroc"], 4),
         "retained": (pathogen, m["stem"]) in retained_keys,
         "discard_reason": discard_reasons.get((pathogen, m["stem"]), "")}
        for m in models
    ]).sort_values("mean_auroc", ascending=False)
    summary_path = os.path.join(output_dir, f"{pathogen}_auroc_summary.csv")
    summary_df.to_csv(summary_path, index=False)
    n_dropped = int((~summary_df["retained"]).sum())
    note = f" ({n_dropped} discarded at step 10)" if n_dropped else ""
    print(f"  -> {summary_path}{note}")

    # Retained-only CV AUROCs feed the condensed cross-pathogen consensus figure.
    consensus_entries.append({
        "code": pathogen,
        "pathogen": pathogen_names.get(pathogen, pathogen),
        "aurocs": summary_df.loc[summary_df["retained"], "mean_auroc"].tolist(),
        "n_molecules": mol_counts.get(pathogen),
    })

    # ROC grid + rank boxplots, both as PNG + vector PDF on the 3 cm cell grid
    figs = save_performance_figures(pathogen, models, individual_dir, pathogen_names)
    footprints.update(figs)
    for name, cells in figs.items():
        print(f"  -> {name} ({cells[0]}x{cells[1]} cells)")

manifest = write_figure_cells(footprints, individual_dir)
print(f"\n-> {manifest}")

# Condensed cross-pathogen figures at the top level of output_dir.
#
# Dataset size and balance per pathogen, from the step-10 report table. Added negatives are excluded
# from BOTH the size and the active fraction: n_compounds counts them, so subtracting them (and the
# decoys, which are zero for every model) reproduces the curation pipeline's own n_mol_after and
# ar_after. Only the 54 models that were given negatives are affected, but for those the change is
# large — mtuberculosis/DR_0012 goes from 2450 compounds at 0.50 active to 1411 at 0.87.
step10 = pd.read_csv(os.path.join(step10_dir, "10_reports.csv"), dtype={"name": str})
datasets = pd.DataFrame({
    "pathogen": step10["pathogen"],
    "dataset": step10["name"],
    "size": step10["n_compounds"] - step10["n_added_negatives"] - step10["n_added_decoys"],
    "n_active": step10["n_positives"],
})
datasets["n_inactive"] = datasets["size"] - datasets["n_active"]
datasets["active_fraction"] = datasets["n_active"] / datasets["size"]
datasets.sort_values(["pathogen", "size"], ascending=[True, False]).to_csv(
    os.path.join(output_dir, "dataset_sizes.csv"), index=False)

top_fp = {}
top_fp.update(save_activity_ratios_figure(datasets, output_dir, pathogen_names))
top_fp.update(save_consensus_figure(consensus_entries, output_dir))
if top_fp:
    write_figure_cells(top_fp, output_dir)
    for name, cells in top_fp.items():
        print(f"-> {name} ({cells[0]}x{cells[1]} cells)")

# Reconcile step 09 (plotted) against step 10 (retained + discarded). Anything in
# neither table means the staged 10_reports/ is out of sync with 09_reports/.
unaccounted = plotted_keys - retained_keys - set(discard_reasons)
missing_reports = retained_keys - plotted_keys
print(f"\nPlotted (step 09): {len(plotted_keys)} | retained (step 10): "
      f"{len(plotted_keys & retained_keys)} | discarded: "
      f"{len(plotted_keys & set(discard_reasons))}")
if unaccounted:
    print(f"  [WARN] {len(unaccounted)} plotted models in neither step-10 table: "
          f"{sorted(unaccounted)[:5]}")
if missing_reports:
    print(f"  [WARN] {len(missing_reports)} retained models have no step-09 report: "
          f"{sorted(missing_reports)[:5]}")
print("Done.")
