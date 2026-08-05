"""Step 07 — the full-library score matrices: named, column-scaled, row-normalized.

The foundation of the correlation analysis. Every annotation model was run on the SAME
~1.35M-compound reference library and staged by ``00_download_data.py`` as
``data/processed/annotation_preds_ref_library/{model_id}_{version}.csv``, all aligned on ``key``, so
the predictions form a clean rectangular matrix with no compound alignment needed.

This script produces **five matrices**, each 1,355,109 rows x 260 columns:

1. **named** — the raw scores for every ``selected == Yes`` endpoint in
   ``config/08_endpoint_selection.csv``, columns renamed
   ``{pathogen_code}__{model_id}__{column_name}`` so pathogen and source model read straight off the
   column name.
2. **z-score** — ``(x - mean) / std`` per column. Standard ML preprocessing; sensitive to outliers,
   allows negative values.
3. **rank-percentile** — each value replaced by its percentile rank within its own column, bounded
   ``[0, 1]``. Robust to outliers and to columns on wildly different native scales; the same
   rank-based idea behind every metric downstream.
4. **z-score + L2 row-normalized** — each compound's 260-endpoint profile divided by its Euclidean
   norm, so profiles compare by shape rather than magnitude. L2 is the conventional norm for signed,
   mean-centred data and the one cosine similarity between profiles is built on.
5. **rank-percentile + L1 row-normalized** — each profile divided by the sum of its absolute values,
   so it sums to 1. Percentile values are non-negative, which makes this a clean compositional
   vector: each endpoint's relative share of that compound's total activity.

Steps 2/3 act on columns, steps 4/5 on rows; the two are independent and composed in that order. The
norms are deliberately *different* between 4 and 5 because the value ranges differ — L1 on signed
z-scores has no equally clean reading.

Reads the ~15 GB of raw prediction CSVs at most **once**, caching the result as
``07_score_matrix_full.parquet`` (every endpoint the selection CSV references, Yes AND No, so
re-enabling one later needs no re-extraction). Every subsequent run, and every downstream step,
re-slices that cache in under a second. Each output CSV is skipped if it already exists.

    python 07_score_matrices.py

Outputs
-------
    output/07_score_matrices/07_score_matrix_full.parquet   (cache: all referenced endpoints)
    output/07_score_matrices/07_score_matrix_named.csv
    output/07_score_matrices/07_score_matrix_named_zscore.csv
    output/07_score_matrices/07_score_matrix_named_rankpct.csv
    output/07_score_matrices/07_score_matrix_named_zscore_l2rownorm.csv
    output/07_score_matrices/07_score_matrix_named_rankpct_l1rownorm.csv
"""

import os
import sys
import time

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))

from default import ANNOTATION_PREDS_SUBDIR  # noqa: E402
from eval_correlations import (  # noqa: E402
    build_full_library_matrix,
    build_named_score_matrix,
    row_normalize,
    scale_matrix,
)

pred_dir = os.path.join(root, "..", "data", "processed", ANNOTATION_PREDS_SUBDIR)
output_dir = os.path.join(root, "..", "output", "07_score_matrices")
config_dir = os.path.join(root, "..", "config")
os.makedirs(output_dir, exist_ok=True)

endpoint_selection_path = os.path.join(config_dir, "08_endpoint_selection.csv")
pathogens_of_interest_path = os.path.join(config_dir, "pathogens_of_interest.csv")
full_matrix_cache_path = os.path.join(output_dir, "07_score_matrix_full.parquet")

#: label -> (output filename, transform applied to the named matrix)
VARIANTS = {
    "named": ("07_score_matrix_named.csv", lambda m: m),
    "zscore": ("07_score_matrix_named_zscore.csv", lambda m: scale_matrix(m, "zscore")),
    "rankpct": ("07_score_matrix_named_rankpct.csv", lambda m: scale_matrix(m, "rank_pct")),
    "zscore+l2": ("07_score_matrix_named_zscore_l2rownorm.csv",
                  lambda m: row_normalize(scale_matrix(m, "zscore"), "l2")),
    "rankpct+l1": ("07_score_matrix_named_rankpct_l1rownorm.csv",
                   lambda m: row_normalize(scale_matrix(m, "rank_pct"), "l1")),
}

if not os.path.exists(full_matrix_cache_path):
    print(f"[score-matrices] building full-library cache from {pred_dir} "
          "(reads ~15 GB of raw prediction files — this can take several minutes)...")
    t0 = time.time()
    full = build_full_library_matrix(pred_dir, endpoint_selection_path)
    full.to_parquet(full_matrix_cache_path)
    print(f"[score-matrices] wrote {full_matrix_cache_path} {full.shape} in "
          f"{time.time() - t0:.1f}s")
    del full

pending = {k: v for k, v in VARIANTS.items()
           if not os.path.exists(os.path.join(output_dir, v[0]))}
if not pending:
    print(f"[skip] all {len(VARIANTS)} matrices already exist in {output_dir}. "
          "Delete one to rebuild it.")
else:
    t0 = time.time()
    base = build_named_score_matrix(
        pred_dir=pred_dir, endpoint_selection_path=endpoint_selection_path,
        pathogens_of_interest_path=pathogens_of_interest_path,
        full_matrix_cache_path=full_matrix_cache_path)
    print(f"[score-matrices] named matrix {base.shape} ready in {time.time() - t0:.1f}s")

    for label, (fname, transform) in pending.items():
        out_path = os.path.join(output_dir, fname)
        t0 = time.time()
        matrix = transform(base)
        print(f"[score-matrices] {label} computed in {time.time() - t0:.1f}s {matrix.shape} "
              "— writing CSV...")
        t0 = time.time()
        matrix.to_csv(out_path)
        print(f"[score-matrices] wrote {fname} ({os.path.getsize(out_path) / 1e9:.2f} GB) in "
              f"{time.time() - t0:.1f}s")
        del matrix

print(f"\nDone → {output_dir}")
