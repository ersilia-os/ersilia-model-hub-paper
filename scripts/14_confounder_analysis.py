"""Step 14 — property/resemblance columns as activity predictors, and the endpoint-confounder check.

Two parts in one script, both reading step 08's cached property matrices directly (never rebuilding
one, never touching a raw prediction file):

1. **Predictor performance** — treats every column of the physchem, abx (step 08) and cytotox
   (step 08) property blocks as a PREDICTOR, and every curated activity endpoint as a binary TARGET,
   giving one performance value per (predictor, target) pair — **101 x 307 = 31,007** as the
   selection now stands.

2. **Endpoint confounder check** (absorbed from the former standalone step 15, 2026-09-02) — for
   every endpoint in step 09's ``09_endpoint_quality.csv`` ranking, the strongest non-same-model
   physchem/abx/cytotox predictor of its own top-1000, cross-referenced against step 09's weakest
   endpoints. An endpoint whose top-1000 is well recovered by a single property column is, whatever
   its peer agreement, largely reproducing a property filter rather than a bioactivity ranking. Kept
   as its own output for the moment rather than merged into step 09's table, since that table is
   built before this step has run.

Physchem's own per-descriptor stats, distributions grid and UMAP panel (formerly a "part 1" here,
before that a standalone ``10_physchem_matrix.py``) were removed 2026-09-02 (user-directed): they
tested nothing about activity, so they didn't belong in a step whose whole point is the confound
check. Physchem's raw matrix is unaffected and still feeds part 1 below as one of the three
predictor families; only its standalone descriptive figures are gone. Recoverable from `git log` if
the per-descriptor stats/distributions are ever needed again.

Both parts are a descriptive association measure, NOT a trained model: nothing is fitted, nothing is
split, and no random seed is involved. Every value is a rank statistic over the full 1.35M-compound
library.

Binarization of the targets: the ``ACTIVITY_BINARIZE_TOP_N`` = 1000 highest-scoring compounds are the
positive class, all remaining 1,354,109 are negative — a user-directed RANK cutoff on a fixed count,
never a score threshold. All selected endpoints are ``direction == higher`` (asserted), so "top" is
unambiguous.

Metric per predictor, from its own value type resolved on the full column:
    continuous -> AUROC (Mann-Whitney rank-sum identity; exact, and ~100x faster at this scale)
    binary     -> balanced accuracy
Both share a 0.5 chance baseline, which is what lets them share a y-axis. AUROC is reported RAW and
may fall below 0.5, so an anti-correlated predictor stays visible as such.

    python 08_property_matrices.py         # physchem/abx/cytotox predictors
    python 07_score_matrices.py            # activity targets
    python 09_bioactivity_endpoints.py     # endpoint quality ranking (part 3)
    python 14_confounder_analysis.py

Outputs
-------
    output/14_confounder_analysis/14_predictor_performance.csv
    output/14_confounder_analysis/14_predictor_summary.csv
    output/14_confounder_analysis/14_curated_predictor_performance.csv
    output/14_confounder_analysis/14_endpoint_confounders.csv
    output/14_confounder_analysis/png|pdf/14_performance_{family}.*
    output/14_confounder_analysis/png|pdf/14_performance_curated_predictors.*
    output/14_confounder_analysis/figure_cells.json
"""

import os
import sys

import pandas as pd

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))

from eval_endpoint_quality import confounder_stats  # noqa: E402
from eval_predictor_performance import run_predictor_performance  # noqa: E402
from plots_predictor_performance import (  # noqa: E402
    save_curated_predictor_figure, save_predictor_performance_figures,
)

config_dir = os.path.join(root, "..", "config")
output_dir = os.path.join(root, "..", "output", "14_confounder_analysis")
os.makedirs(output_dir, exist_ok=True)

property_matrices_dir = os.path.join(root, "..", "output", "08_property_matrices")
jaccard_output_dir = os.path.join(root, "..", "output", "09_bioactivity_endpoints")

#: One matrix CSV per property family, in family order. The family of each column is read from its
#: own `{family}__` prefix, so this list only decides which files are opened, never how a predictor
#: is labelled.
physchem_csv = os.path.join(property_matrices_dir, "08_physchem_matrix_named.csv")
abx_csv = os.path.join(property_matrices_dir, "08_abx_matrix_named.csv")
cytotox_csv = os.path.join(property_matrices_dir, "08_cytotox_matrix_named.csv")
property_csvs = [physchem_csv, abx_csv, cytotox_csv]
parquet_path = os.path.join(
    root, "..", "output", "07_score_matrices", "07_score_matrix_full.parquet")
selection_csv = os.path.join(config_dir, "08_endpoint_selection.csv")
pathogens_csv = os.path.join(config_dir, "pathogens_of_interest.csv")
endpoint_quality_csv = os.path.join(jaccard_output_dir, "09_endpoint_quality.csv")

for path, step in [(physchem_csv, "08_property_matrices.py"),
                   (abx_csv, "08_property_matrices.py"),
                   (cytotox_csv, "08_property_matrices.py"),
                   (parquet_path, "07_score_matrices.py"),
                   (endpoint_quality_csv, "09_bioactivity_endpoints.py")]:
    if not os.path.exists(path):
        sys.exit(f"Missing {path}. Run `python {step}` first.")

# --------------------------------------------------------------------------- #
# 1. Predictor performance                                                     #
# --------------------------------------------------------------------------- #
perf, summary, curated = run_predictor_performance(
    property_csvs=property_csvs, parquet_path=parquet_path,
    selection_csv=selection_csv, pathogens_csv=pathogens_csv, output_dir=output_dir)
save_predictor_performance_figures(output_dir, perf)
save_curated_predictor_figure(output_dir, curated)

# --------------------------------------------------------------------------- #
# 2. Endpoint confounder check (former step 15)                                #
# --------------------------------------------------------------------------- #
# Scope and ranking come from step 9's endpoint-quality table (pathogens of interest with >5
# endpoints, weakest AUROC-upranking first); the confounder itself is computed fresh here, since it
# needs THIS step's property-predictor data, which does not exist at step-09 time.
endpoint_quality = pd.read_csv(endpoint_quality_csv)
conf = confounder_stats(perf, endpoint_quality).reset_index()
conf_path = os.path.join(output_dir, "14_endpoint_confounders.csv")
conf.to_csv(conf_path, index=False)
print(f"\n[endpoint-confounders] {len(conf)} of {len(endpoint_quality)} in-scope endpoints have a "
      f"non-same-model property predictor -> {os.path.basename(conf_path)}")

# Not merged into step 9's table (per the current pass, "for the moment") — just joined here for
# the printed report, so the weakest endpoints by AUROC-upranking can be read against their
# strongest confounder side by side.
report = endpoint_quality.merge(conf, on="endpoint", how="left")
cols = ["organism", "model_id", "column_name", "auroc_out_same_median", "auroc_out_specificity",
        "confounder_predictor", "confounder_family", "confounder_value", "confounder_abs_dev"]
print("\n[endpoint-confounders] weakest 15 endpoints (step 9's ranking) with their strongest "
      "confounder:")
print(report[cols].head(15).to_string(index=False, float_format=lambda v: f"{v:.4f}"))

print(f"\nDone -> {output_dir}")
