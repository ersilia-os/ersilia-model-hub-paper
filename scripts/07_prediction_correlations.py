"""Step 07 — Inter-model prediction correlations on the reference library (initial assessment).

The ~130 annotation models were all run on the SAME ~1.35M-compound reference library and staged by
``00_download_data.py`` as ``data/processed/annotation_preds_ref_library/{model_id}_{version}.csv``.
This step asks how much the models AGREE, at the output-column level (each output column is a node):

  - a library-wide **Spearman** correlation matrix across all ~1500 output columns, computed on a
    fixed-seed ``CORR_SAMPLE_N``-compound sample (raise it for a final run);
  - **top-N Jaccard overlap** (N = 100, 500) of the highest-scoring compounds, for probability-type
    columns only (a known "higher = more" direction);
  - auto-assigned focus groups (Target Organism; cytotoxicity via a regex over Airtable text) that
    the user REVIEWS before figures — motivating the "do cytotoxicity models / same-organism models
    correlate?" questions.

Two stages, gated on a manual review of the auto group assignments:

    python 07_prediction_correlations.py                 # build stage (default): matrix + groups
    #   -> review output/07_prediction_correlations/07_group_assignments.csv
    python 07_prediction_correlations.py --analyze        # correlation + overlap + figures
    python 07_prediction_correlations.py --all            # both, without stopping (skips review)

The build stage caches ``07_score_matrix.parquet``; re-runs reuse it. This is an EXPLORATORY step:
its outputs guide a narrower follow-up.

Outputs
-------
    output/07_prediction_correlations/07_score_matrix.parquet        (sampled key × node matrix)
    output/07_prediction_correlations/07_column_index.csv            (per-node value type + groups)
    output/07_prediction_correlations/07_group_assignments.csv       (auto groups — REVIEW THIS)
    output/07_prediction_correlations/07_spearman_corr.csv           (node × node Spearman)
    output/07_prediction_correlations/07_topn_overlap_N100.csv        (probability-node Jaccard)
    output/07_prediction_correlations/07_topn_overlap_N500.csv
    output/07_prediction_correlations/07_group_correlation_summary.csv
    output/07_prediction_correlations/png/<panel>.png
    output/07_prediction_correlations/pdf/<panel>.pdf
    output/07_prediction_correlations/figure_cells.json
"""

import argparse
import os
import sys

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))

from default import ANNOTATION_PREDS_SUBDIR  # noqa: E402
from eval_correlations import run_analyze, run_build  # noqa: E402
from plots_correlations import save_correlation_figures  # noqa: E402

parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
parser.add_argument("--analyze", action="store_true",
                    help="Run the analyze stage (needs a built matrix + reviewed groups).")
parser.add_argument("--all", action="store_true",
                    help="Run build then analyze without stopping for group review.")
args = parser.parse_args()

pred_dir = os.path.join(root, "..", "data", "processed", ANNOTATION_PREDS_SUBDIR)
meta_path = os.path.join(root, "..", "data", "raw", "airtable_metadata.csv")
output_dir = os.path.join(root, "..", "output", "07_prediction_correlations")
os.makedirs(output_dir, exist_ok=True)

if not args.analyze:
    run_build(pred_dir=pred_dir, meta_path=meta_path, output_dir=output_dir)

if args.analyze or args.all:
    run_analyze(output_dir=output_dir)
    save_correlation_figures(output_dir)
    print(f"\nDone → {output_dir}")
