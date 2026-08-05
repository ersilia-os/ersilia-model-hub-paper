"""Step 08 — per-pathogen top-1000 Jaccard: same pathogen vs. different pathogen.

Asks whether endpoints of the SAME pathogen pick out the same compounds more than endpoints of
DIFFERENT pathogens do. For every pair of the 260 selected endpoint columns, the Jaccard overlap of
their top-1000 highest-scoring compounds (out of 1,355,109), then aggregated per pathogen: turquoise =
every unordered pair of that pathogen's own columns, crimson = its columns against every other
pathogen's.

**Three matrix variants are computed, not five.** Top-N Jaccard depends only on each column's own
internal ranking, and both column scalings from step 07 are strictly increasing per column — so
z-scoring and rank-percentiling cannot change any column's top-1000 set. The unscaled matrix and both
scaled matrices therefore give the same result, computed once and labelled as covering all three. Only
the row-normalized matrices change rankings, because each row is divided by a different scalar. The
identity is **asserted at runtime**, not assumed — and the assertion earns its keep: baseline ==
rank-percentiled holds exactly, but baseline == z-scored comes back FALSE, because ``(x - mean) / std``
in float32 reorders near-tied values in one column of 260 (~0.001 on 138 of 67,600 Jaccard cells).

The 260x260 Jaccard matrices are written out for reuse by later analysis.

Reporting choices worth knowing:

  - **Minimum endpoints per pathogen** is a command-line argument (default 5). A pathogen below it is
    removed from the analysis ENTIRELY — it stops being a different-pathogen partner too, not merely
    loses its own box. Pass several values to compare thresholds side by side; each writes its own
    ``min<K>`` outputs rather than overwriting.
  - **Same-model pairs are INCLUDED** in the boxes — the literal "each column against all others". Two
    output columns of one model agreeing says little about pathogen specificity, so the summary CSV
    carries ``same_median_excl_same_model`` alongside ``same_median``. Read them together: at
    ``min5``, three of the eleven surviving pathogens have NO cross-model same-pathogen pair at all.
  - **Linear x-axis** — values bunch near zero, but exact-zero pairs render instead of being silently
    dropped by a log axis. Nothing is filtered.

    python 08_pathogen_jaccard.py          # min 5
    python 08_pathogen_jaccard.py 2 5      # both, for comparison

Outputs (per threshold K)
-------
    output/08_pathogen_jaccard/08_jaccard_top1000_<variant>_matrix.csv          (260x260, reused)
    output/08_pathogen_jaccard/08_pathogen_jaccard_top1000_min<K>_<variant>_summary.csv
    output/08_pathogen_jaccard/png/08_pathogen_jaccard_top1000_min<K>_<variant>.png
    output/08_pathogen_jaccard/pdf/08_pathogen_jaccard_top1000_min<K>_<variant>.pdf
"""

import os
import sys
import time

import numpy as np
import pandas as pd

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))

from default import ANNOTATION_PREDS_SUBDIR  # noqa: E402
from eval_correlations import (  # noqa: E402
    build_named_score_matrix,
    column_metric_pairs,
    multi_column_pathogen_nodes,
    parse_named_column,
    pathogen_metric_boxes,
    pathogen_metric_summary,
    row_normalize,
    scale_matrix,
    topn_jaccard_matrix,
)
from plots_matrix_analyses import pathogen_jaccard_figure  # noqa: E402

CUTOFF = 1000
#: Minimum endpoints a pathogen must have OF ITS OWN to enter the analysis. User-directed, not
#: fitted; overridable on the command line (see the module docstring).
DEFAULT_MIN_COLUMNS = (5,)
min_columns_values = [int(a) for a in sys.argv[1:]] or list(DEFAULT_MIN_COLUMNS)

pred_dir = os.path.join(root, "..", "data", "processed", ANNOTATION_PREDS_SUBDIR)
matrix_dir = os.path.join(root, "..", "output", "07_score_matrices")
output_dir = os.path.join(root, "..", "output", "08_pathogen_jaccard")
config_dir = os.path.join(root, "..", "config")
os.makedirs(output_dir, exist_ok=True)

endpoint_selection_path = os.path.join(config_dir, "08_endpoint_selection.csv")
pathogens_of_interest_path = os.path.join(config_dir, "pathogens_of_interest.csv")
full_matrix_cache_path = os.path.join(matrix_dir, "07_score_matrix_full.parquet")

# (slug, human label, how to derive the matrix from the base named matrix)
VARIANTS = [
    ("baseline", "unscaled scores (= rank-percentiled; = z-scored bar one column)", lambda m: m),
    ("zscore_l2rownorm", "z-scored columns, L2 row-normalized",
     lambda m: row_normalize(scale_matrix(m, "zscore"), "l2")),
    ("rankpct_l1rownorm", "rank-percentiled columns, L1 row-normalized",
     lambda m: row_normalize(scale_matrix(m, "rank_pct"), "l1")),
]

# ----------------------------------------------------------------------------- #
# 1. The 260x260 Jaccard matrices (cached — computing them is the expensive part)
# ----------------------------------------------------------------------------- #
jaccards = {}
missing = [v for v in VARIANTS
           if not os.path.exists(os.path.join(output_dir, f"08_jaccard_top{CUTOFF}_{v[0]}_matrix.csv"))]
base = None
if missing:
    t0 = time.time()
    base = build_named_score_matrix(
        pred_dir=pred_dir, endpoint_selection_path=endpoint_selection_path,
        pathogens_of_interest_path=pathogens_of_interest_path,
        full_matrix_cache_path=full_matrix_cache_path)
    print(f"[pathogen-jaccard] base named matrix {base.shape} ready in {time.time() - t0:.1f}s")

for slug, label, derive in VARIANTS:
    jac_path = os.path.join(output_dir, f"08_jaccard_top{CUTOFF}_{slug}_matrix.csv")
    if os.path.exists(jac_path):
        jaccards[slug] = pd.read_csv(jac_path, index_col=0)
        print(f"[pathogen-jaccard] {slug}: reusing cached {os.path.basename(jac_path)}")
        continue
    t0 = time.time()
    variant = derive(base)
    jac = topn_jaccard_matrix(variant, CUTOFF)
    del variant
    jac.to_csv(jac_path)
    jaccards[slug] = jac
    print(f"[pathogen-jaccard] {slug}: top-{CUTOFF} Jaccard {jac.shape} in "
          f"{time.time() - t0:.1f}s -> {os.path.basename(jac_path)}")

# The scaling-invariance claim is asserted, not assumed: if a future change to scale_matrix broke
# monotonicity, the "baseline covers all three" label would silently become a lie.
if base is not None:
    for name, method in (("z-scored", "zscore"), ("rank-percentiled", "rank_pct")):
        other = topn_jaccard_matrix(scale_matrix(base, method), CUTOFF)
        same = np.allclose(jaccards["baseline"].to_numpy(), other.to_numpy(), equal_nan=True)
        print(f"[pathogen-jaccard] check: baseline == {name} top-{CUTOFF} Jaccard -> {same}")
    del base

# ----------------------------------------------------------------------------- #
# 2. Per-pathogen aggregation and figure, per threshold
# ----------------------------------------------------------------------------- #
for min_columns in min_columns_values:
    print(f"\n===== minimum {min_columns} endpoints per pathogen =====")
    for slug, label, _ in VARIANTS:
        jac = jaccards[slug]
        nodes = multi_column_pathogen_nodes(jac, min_columns=min_columns)
        all_pathogens = pd.Series([parse_named_column(c)[0] for c in jac.columns]).value_counts()
        n_columns = pd.Series([parse_named_column(n)[0] for n in nodes]).value_counts()
        dropped = all_pathogens.drop(index=n_columns.index)

        pairs = column_metric_pairs(jac.loc[nodes, nodes])
        boxes = pathogen_metric_boxes(pairs)
        summary = pathogen_metric_summary(boxes, n_columns)
        stem = f"08_pathogen_jaccard_top{CUTOFF}_min{min_columns}_{slug}"
        summary.to_csv(os.path.join(output_dir, f"{stem}_summary.csv"), index=False)

        png_path, pdf_path = pathogen_jaccard_figure(
            boxes, summary, cutoff=CUTOFF,
            matrix_label=f"{label}  ({min_columns}+ endpoints per pathogen)",
            name=stem, output_dir=output_dir)

        print(f"[pathogen-jaccard] {slug}")
        print(f"    kept {len(nodes)} of {len(jac.columns)} columns, "
              f"{len(summary)} of {len(all_pathogens)} pathogens "
              f"(>= {min_columns} endpoints each)")
        print(f"    dropped {len(dropped)} pathogens / {int(dropped.sum())} columns")
        print(f"    wrote {os.path.basename(png_path)} + {os.path.basename(pdf_path)}")
        print(summary.to_string(index=False, float_format=lambda v: f"{v:.4f}"))
        print()

print(f"Done → {output_dir}")
