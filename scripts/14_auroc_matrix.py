"""Step 14 — one score per organism, then the aggregated AUROC matrix.

Each organism's activity endpoints are collapsed into a SINGLE per-compound score before any AUROC is
computed, so a pathogen's weight comes from the pathogen rather than from how many assays it happens
to have (P. falciparum contributes 13 endpoints to this set, H. pylori 1).

    rows (y)  15 organisms, each binarized at ACTIVITY_BINARIZE_TOP_N = 1000 on its aggregate score
    cols (x)  the same 15 aggregates, then cytotoxicity (6), abx resemblance (3), physchem (3) = 27

The merge: an organism's endpoint columns are scaled to percentile ranks within the full library
(`ORGANISM_MERGE_METHOD`) and then averaged (`ORGANISM_MERGE_AGG`). Scaling BEFORE aggregating is the
point — raw scores from different models sit on unrelated ranges, and averaging them directly would
weight whichever endpoint has the widest one.

**Five organisms have exactly one endpoint** (Campylobacter, Enterobacter, E. faecium, H. pylori,
S. pneumoniae), so nothing is merged and their score IS that endpoint's percentile rank; their row is
not the same kind of quantity as E. coli's 11-endpoint mean. A ChEMBL `consensus_score` is also itself
an aggregate over sub-models, so averaging it with individual assay endpoints gives it equal weight to
a single assay. Both follow from merging the endpoints as selected, and are recorded, not corrected.

Scores come from step 07's parquet CACHE, re-scaled here, not from 07_score_matrix_named_rankpct.csv:
that 6.8 GB CSV is stale (260 columns, predating two config changes) and lacks several of these
endpoints. Step 07's own mean-rank section rebuilds from the parquet for the same reason.

Cells carry their printed AUROC on a discrete spectral scale clipped at 0.5 (below-chance values all
render in the light-grey bottom bin; the unclipped numbers are in 14_auroc_matrix.csv).

    python 07_score_matrices.py    # writes the parquet cache
    python 10_physchem_matrix.py
    python 11_abx_resemblance_matrix.py
    python 12_cytotox_matrix.py
    python 14_auroc_matrix.py

Outputs
-------
    output/14_auroc_matrix/14_auroc_matrix.csv        the 15 x 27 matrix, wide form
    output/14_auroc_matrix/14_auroc_matrix_axes.csv   axis order + class / block / model / n_endpoints
    output/14_auroc_matrix/14_organism_scores.parquet the 1.35M x 15 aggregate score matrix
    output/14_auroc_matrix/png|pdf/14_auroc_matrix.*
    output/14_auroc_matrix/figure_cells.json
"""

import os
import sys

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))

from eval_auroc_matrix import (  # noqa: E402
    aggregated_matrix, axes_table, bioactivity_order, diagonal_check, organism_order,
    organism_scores, predictor_order,
)
from eval_auroc_matrix import overlap_matrix  # noqa: E402
from plots_auroc_matrix import (  # noqa: E402
    save_auroc_matrix_figure, save_overlap_matrix_figure,
)

config_dir = os.path.join(root, "..", "config")
output_dir = os.path.join(root, "..", "output", "14_auroc_matrix")
os.makedirs(output_dir, exist_ok=True)

selection_csv = os.path.join(config_dir, "08_endpoint_selection.csv")
pathogens_csv = os.path.join(config_dir, "pathogens_of_interest.csv")
parquet_path = os.path.join(
    root, "..", "output", "07_score_matrices", "07_score_matrix_full.parquet")
#: One matrix CSV per property family. Column families come from each column's own `{family}__`
#: prefix, so this list only decides which files are opened.
physchem_csv = os.path.join(
    root, "..", "output", "10_physchem_matrix", "10_physchem_matrix_named.csv")
abx_csv = os.path.join(
    root, "..", "output", "11_abx_resemblance_matrix", "11_abx_matrix_named.csv")
cytotox_csv = os.path.join(
    root, "..", "output", "12_cytotox_matrix", "12_cytotox_matrix_named.csv")
property_csvs = [physchem_csv, abx_csv, cytotox_csv]

for path, step in [(parquet_path, "07_score_matrices.py"),
                   (physchem_csv, "10_physchem_matrix.py"),
                   (abx_csv, "11_abx_resemblance_matrix.py"),
                   (cytotox_csv, "12_cytotox_matrix.py")]:
    if not os.path.exists(path):
        sys.exit(f"Missing {path}. Run `python {step}` first.")

# Only endpoints actually present in the cache; the config can select one staged after it was built.
cached = set(pq.ParquetFile(parquet_path).schema.names)
endpoints = bioactivity_order(selection_csv, pathogens_csv, available=cached)
organisms = organism_order(endpoints)
print(f"[auroc-matrix] {len(endpoints)} endpoints -> {len(organisms)} organisms")

scores = organism_scores(parquet_path, endpoints)

# A mean of percentile ranks must sit at ~0.5; a value far off means the scaling ran along the wrong
# axis, which would otherwise produce a plausible-looking but meaningless matrix.
means = scores.mean()
if not np.allclose(means, 0.5, atol=0.02):
    off = means[(means - 0.5).abs() > 0.02]
    sys.exit(f"Aggregate score column mean(s) far from 0.5 — scaling axis is wrong:\n{off}")
print(f"[organism-scores] column means {means.min():.4f}-{means.max():.4f} (must be ~0.5)")

scores.to_parquet(os.path.join(output_dir, "14_organism_scores.parquet"))

cols = predictor_order()
matrix, cols = aggregated_matrix(scores, organisms, cols, property_csvs)

diag = diagonal_check(matrix, organisms["organism"].tolist())
print(f"  diagonal (self-pairs): min {diag.min():.6f}, max {diag.max():.6f} "
      f"({len(diag)} cells — must all be 1.0)")
if not np.allclose(diag, 1.0, atol=1e-3):
    sys.exit("Diagonal is not 1.0 — rows and columns are misaligned; refusing to draw the matrix.")

matrix.to_csv(os.path.join(output_dir, "14_auroc_matrix.csv"))
axes_table(organisms, cols).to_csv(
    os.path.join(output_dir, "14_auroc_matrix_axes.csv"), index=False)
print(f"  -> 14_auroc_matrix.csv ({matrix.shape[0]} x {matrix.shape[1]})")

footprints = save_auroc_matrix_figure(output_dir, matrix, organisms, cols)

# --- Second view: how many of the row's actives are in the column's top 1000 ---
overlap = overlap_matrix(scores, organisms, cols, property_csvs)
overlap.to_csv(os.path.join(output_dir, "14_overlap_matrix.csv"))
print(f"  -> 14_overlap_matrix.csv ({overlap.shape[0]} x {overlap.shape[1]})")
save_overlap_matrix_figure(output_dir, overlap, organisms, cols, footprints)
print(f"\nDone -> {output_dir}")
