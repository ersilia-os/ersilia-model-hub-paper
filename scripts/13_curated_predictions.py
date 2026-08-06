"""Step 13 — do the property/resemblance columns predict pathogen activity?

Treats every column of the step-10 (physchem), step-11 (abx) and step-12 (cytotox) property blocks as a
PREDICTOR, and every curated activity endpoint as a binary TARGET, giving one performance value per
(predictor, target) pair — **101 x 300 = 30,300** as the selection now stands (was 101 x 260 = 26,260
before the M. tuberculosis endpoints were enabled). One figure per predictor family shows each
predictor's distribution of performance across all targets.

This is a descriptive association measure, NOT a trained model: nothing is fitted, nothing is split,
and no random seed is involved. Every value is a rank statistic over the full 1.35M-compound library.

Binarization of the targets: the ``ACTIVITY_BINARIZE_TOP_N`` = 1000 highest-scoring compounds are the
positive class, all remaining 1,354,109 are negative — a user-directed RANK cutoff on a fixed count,
never a score threshold. All 260 selected endpoints are ``direction == higher`` (asserted), so "top"
is unambiguous.

Metric per predictor, from its own value type resolved on the full column:
    continuous -> AUROC (Mann-Whitney rank-sum identity; exact, and ~100x faster at this scale)
    binary     -> balanced accuracy
Both share a 0.5 chance baseline, which is what lets them share a y-axis. AUROC is reported RAW and
may fall below 0.5, so an anti-correlated predictor stays visible as such.

    python 10_physchem_matrix.py           # physchem predictors
    python 11_abx_resemblance_matrix.py    # abx predictors
    python 12_cytotox_matrix.py            # cytotox predictors
    python 07_score_matrices.py                   # activity targets
    python 13_curated_predictions.py

A fourth figure turns the same machinery on the activity endpoints themselves: each endpoint's RAW
score as the predictor against every endpoint's binarized version as the target (300 x 300), grouped
on the x-axis by the predictor endpoint's organism and coloured by whether the target belongs to the
same organism. Self-pairs (an endpoint against its own binarization) are 1.0 by construction and are
excluded from that figure, retained in the CSV only as a correctness check.

Outputs
-------
    output/13_curated_predictions/13_predictor_performance.csv
    output/13_curated_predictions/13_predictor_summary.csv
    output/13_curated_predictions/13_activity_self_performance.csv
    output/13_curated_predictions/png|pdf/13_performance_{family}.*
    output/13_curated_predictions/png|pdf/13_performance_activity_by_organism.*
    output/13_curated_predictions/figure_cells.json
"""

import os
import sys

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))

from eval_predictor_performance import run_all  # noqa: E402
from plots_predictor_performance import (  # noqa: E402
    save_activity_self_figure, save_curated_predictor_figure, save_pathogen_subset_figure,
    save_predictor_performance_figures,
)

config_dir = os.path.join(root, "..", "config")
output_dir = os.path.join(root, "..", "output", "13_curated_predictions")
os.makedirs(output_dir, exist_ok=True)

#: One matrix CSV per property family, in family order. The family of each column is read from its
#: own `{family}__` prefix, so this list only decides which files are opened, never how a predictor
#: is labelled.
physchem_csv = os.path.join(
    root, "..", "output", "10_physchem_matrix", "10_physchem_matrix_named.csv")
abx_csv = os.path.join(
    root, "..", "output", "11_abx_resemblance_matrix", "11_abx_matrix_named.csv")
cytotox_csv = os.path.join(
    root, "..", "output", "12_cytotox_matrix", "12_cytotox_matrix_named.csv")
property_csvs = [physchem_csv, abx_csv, cytotox_csv]
parquet_path = os.path.join(
    root, "..", "output", "07_score_matrices", "07_score_matrix_full.parquet")
selection_csv = os.path.join(config_dir, "08_endpoint_selection.csv")
pathogens_csv = os.path.join(config_dir, "pathogens_of_interest.csv")

for path, step in [(physchem_csv, "10_physchem_matrix.py"),
                   (abx_csv, "11_abx_resemblance_matrix.py"),
                   (cytotox_csv, "12_cytotox_matrix.py"),
                   (parquet_path, "07_score_matrices.py")]:
    if not os.path.exists(path):
        sys.exit(f"Missing {path}. Run `python {step}` first.")

perf, summary, self_perf, subset, curated = run_all(
    property_csvs=property_csvs, parquet_path=parquet_path,
    selection_csv=selection_csv, pathogens_csv=pathogens_csv, output_dir=output_dir)
footprints = save_predictor_performance_figures(output_dir, perf)
footprints = save_activity_self_figure(output_dir, self_perf, footprints)
footprints = save_pathogen_subset_figure(output_dir, subset, footprints)
save_curated_predictor_figure(output_dir, curated, footprints)
print(f"\nDone -> {output_dir}")
