"""Step 08 — per-pathogen top-1000 Jaccard: same pathogen vs. different pathogen.

Asks whether endpoints of the SAME pathogen pick out the same compounds more than endpoints of
DIFFERENT pathogens do. For every pair of the selected endpoint columns, the Jaccard overlap of their
top-1000 highest-scoring compounds (out of 1,355,109), then aggregated per pathogen: turquoise =
every unordered pair of that pathogen's own columns, crimson = its columns against every other
pathogen's.

**Scope: the 15 curated pathogens of interest** (``config/pathogens_of_interest.csv``), replacing the
former ``min<K>``-endpoint thresholds (2026-08-07). Restricting the node set removes every other
pathogen ENTIRELY — the ~40 gut-microbiome and other organisms stop being different-pathogen partners
too — so each pathogen's crimson box is against the other 14 of interest only, **not** against all 57.
That is a deliberately narrower comparator than the old ``min5`` figure used, and a caption must say
so: "specific to this pathogen" here means "relative to the other priority pathogens".

**Two of the 15 have a single endpoint** (*Campylobacter*, *H. pylori*), so they have no
same-pathogen pair and no turquoise box — only a crimson one, and a NaN ``same_median``. They are
kept rather than dropped: for a pathogen on the priority list, having too little in the hub to assess
is itself the result, and hiding the row would hide it.

**Three matrix variants are computed, not five.** Top-N Jaccard depends only on each column's own
internal ranking, and both column scalings from step 07 are strictly increasing per column — so
z-scoring and rank-percentiling cannot change any column's top-1000 set. The unscaled matrix and both
scaled matrices therefore give the same result, computed once and labelled as covering all three. Only
the row-normalized matrices change rankings, because each row is divided by a different scalar. The
identity is **asserted at runtime**, not assumed — and the assertion earns its keep: baseline ==
rank-percentiled holds exactly, but baseline == z-scored comes back FALSE, because ``(x - mean) / std``
in float32 reorders near-tied values in one column of 300 (~0.001 on 156 of 90,000 Jaccard cells).

The full 300x300 Jaccard matrices are written out for reuse by later analysis — they are computed over
ALL pathogens, and only the per-pathogen aggregation is restricted to the 15.

Reporting choices worth knowing:

  - **Same-model pairs are INCLUDED** in the boxes — the literal "each column against all others". Two
    output columns of one model agreeing says little about pathogen specificity, so the summary CSV
    carries ``same_median_excl_same_model`` alongside ``same_median``. Read them together: several
    pathogens have NO cross-model same-pathogen pair at all, so their turquoise box is one model
    agreeing with itself.
  - **Linear x-axis** — values bunch near zero, but exact-zero pairs render instead of being silently
    dropped by a log axis. Nothing is filtered.

  This is a DIAGNOSTIC, not a paper panel: plain matplotlib rather than the 3 cm cell grid, because
  the y axis carries per-pathogen endpoint and pair counts that go illegible at page width.

    python 08_pathogen_jaccard.py

Outputs
-------
    output/08_pathogen_jaccard/08_jaccard_top1000_<variant>_matrix.csv          (300x300, reused)
    output/08_pathogen_jaccard/08_pathogen_jaccard_top1000_<variant>_summary.csv
    output/08_pathogen_jaccard/png/08_pathogen_jaccard_top1000_<variant>.png
    output/08_pathogen_jaccard/pdf/08_pathogen_jaccard_top1000_<variant>.pdf
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
    pathogens_of_interest_nodes,
    parse_named_column,
    pathogen_metric_boxes,
    pathogen_metric_summary,
    row_normalize,
    scale_matrix,
    topn_jaccard_matrix,
)
from plots_matrix_analyses import pathogen_jaccard_figure  # noqa: E402

CUTOFF = 1000

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
for slug, label, _ in VARIANTS:
    jac = jaccards[slug]
    nodes = pathogens_of_interest_nodes(jac, pathogens_of_interest_path, endpoint_selection_path)
    all_pathogens = pd.Series([parse_named_column(c)[0] for c in jac.columns]).value_counts()
    n_columns = pd.Series([parse_named_column(n)[0] for n in nodes]).value_counts()
    dropped = all_pathogens.drop(index=n_columns.index)

    pairs = column_metric_pairs(jac.loc[nodes, nodes])
    boxes = pathogen_metric_boxes(pairs)
    summary = pathogen_metric_summary(boxes, n_columns)
    stem = f"08_pathogen_jaccard_top{CUTOFF}_{slug}"
    summary.to_csv(os.path.join(output_dir, f"{stem}_summary.csv"), index=False)

    png_path, pdf_path = pathogen_jaccard_figure(
        boxes, summary, cutoff=CUTOFF,
        matrix_label=f"{label}  (15 pathogens of interest)",
        name=stem, output_dir=output_dir)

    print(f"[pathogen-jaccard] {slug}")
    print(f"    kept {len(nodes)} of {len(jac.columns)} columns, "
          f"{len(summary)} of {len(all_pathogens)} pathogens (the curated 15)")
    print(f"    excluded {len(dropped)} other pathogens / {int(dropped.sum())} columns — they are "
          "no longer different-pathogen partners either")
    # A priority pathogen with one endpoint has no same-pathogen pair. It is KEPT, with a diff box
    # and a NaN same_median, because the gap is the finding for a pathogen on the priority list.
    no_same = summary[summary["n_same_pairs"] == 0]
    if len(no_same):
        print(f"    {len(no_same)} pathogen(s) with a single endpoint, so no same-pathogen pair "
              f"and no same box: {', '.join(no_same['pathogen'])}")
    print(f"    wrote {os.path.basename(png_path)} + {os.path.basename(pdf_path)}")
    print(summary.to_string(index=False, float_format=lambda v: f"{v:.4f}"))
    print()

print(f"Done → {output_dir}")
