"""Step 10 — one score per organism, then the aggregated AUROC matrix.

Each organism's activity endpoints are collapsed into a SINGLE per-compound score before any AUROC is
computed, so a pathogen's weight comes from the pathogen rather than from how many assays it happens
to have (P. falciparum contributes 13 endpoints to this set, H. pylori 1).

    rows (y)  15 organisms, each binarized at ACTIVITY_BINARIZE_TOP_N = 1000 on its aggregate score
    cols (x)  the same 15 aggregates, then cytotoxicity (1) and abx resemblance (1) = 17

The merge: an organism's endpoint columns are scaled to percentile ranks within the full library
(`ORGANISM_MERGE_METHOD`) and then averaged (`ORGANISM_MERGE_AGG`). Scaling BEFORE aggregating is the
point — raw scores from different models sit on unrelated ranges, and averaging them directly would
weight whichever endpoint has the widest one.

**Five organisms have exactly one endpoint** (Campylobacter, Enterobacter, E. faecium, H. pylori,
S. pneumoniae), so nothing is merged and their score IS that endpoint's percentile rank; their row is
not the same kind of quantity as E. coli's 12-endpoint mean. A ChEMBL `consensus_score` is also itself
an aggregate over sub-models, so averaging it with individual assay endpoints gives it equal weight to
a single assay. Both follow from merging the endpoints as selected, and are recorded, not corrected.

Scores come from step 07's parquet CACHE, re-scaled here, not from 07_score_matrix_named_rankpct.csv:
that 6.8 GB CSV is stale (260 columns, predating two config changes) and lacks several of these
endpoints. Step 07's own mean-rank section rebuilds from the parquet for the same reason.

**Property predictors are two MERGED rank-sum columns, not 12 raw ones** (2026-09-02, user-directed —
see `eval_auroc_matrix.merged_predictor_scores` and `scripts/README.md`): cytotoxicity's 6 columns and
abx-resemblance's 3 are each percentile-ranked then summed per compound, collapsing to one column per
family. Physchem (mw/clogp/tpsa) was dropped rather than merged, so it is not a dependency here.

**Depends on nothing but steps 07 and 08.** The abx and cytotox matrices are read straight from step
08's cache, not from steps 12/13's own output — those two now only draw figures, on a background from
step 11 (reference-library projection), which is why this used to be a circular dependency before the
2026-09 reorg: the old numbering had this analysis last (step 14), needing steps that themselves
needed the projection step, which needed THIS analysis for one panel. Reading 08 directly instead of
12/13 breaks that; moving this step to slot 10 (right after 09, right before 11) means the projection
step's one panel that needs this analysis now works on the very first run. See `scripts/README.md`'s
"Pipeline reorganization" notes.

Cells carry their printed AUROC on a discrete spectral scale clipped at 0.5 (below-chance values all
render in the light-grey bottom bin; the unclipped numbers are in 10_auroc_matrix_phylo.csv).

    python 07_score_matrices.py    # writes the parquet cache
    python 08_property_matrices.py # abx/cytotox predictor matrices
    python 10_auroc_matrix.py

A SECOND view follows the AUROC matrix: how many of the row organism's actives are among the
column's own top N, at the standard `ACTIVITY_BINARIZE_TOP_N` = 1000 cutoff (the other two cutoffs
this view was compared at, 100 and 10000, were dropped 2026-09-01 — see `scripts/README.md`).

Rows are grouped by `ORGANISM_CLASS_ORDER`; within each class, both matrices are drawn ordered by real
NCBI-taxonomy lineage (`eval_auroc_matrix.phylogeny_organism_order`), with a matching dendrogram
diagnostic. The plain alphabetical-within-class baseline and the alternative hierarchical-clustering
(`hclust`) order were both tried as a comparison and dropped once it had done its job — `phylo` is
what was kept. See `scripts/README.md` for the full history.

Non-abx robustness check (absorbed from the former standalone `xx_non_abx_matrix.py`, 2026-09-02)
----------------------------------------------------------------------------------------------------
Repeats the whole analysis above — the exact same two figures, same `phylo` row order, same
`continuous_color` overlap scale — on a library filtered to purge antibiotic-like compounds
(`ABX_FILTERS`, 8 AND'd conditions from step 08's abx block, keeps 97.62% / 1,322,835 of 1,355,109).
If the organism-to-organism structure above were driven by known antibiotics converging on shared
scaffolds, the matrix would be reporting the library's composition rather than the models' behaviour.
Two figures were dropped folding this in — the 100/10000 overlap cutoffs and the discrete colour
scale the standalone script used — so this section draws literally the same two plots as the main
analysis, not a fuller variant; both were already retired from the MAIN analysis on 2026-09-01/02, so
this only brings the non-abx section into line with it, not a new decision. `ABX_FILTERS`'
thresholds, including the one graded value (`abx_score < 0.5`), are unchanged from the standalone
script — kept as-is for now. See `scripts/README.md`'s "Cutoffs reviewed and kept on 2026-08-11" for
why `eos19mt`'s 38 antibiotic-class flags and `eos6ojg`'s other 9 similarity-count columns were
deliberately left out (sulfonamides alone is 76,473 of 80,234 eos19mt hits — too promiscuous a motif
outside antibiotics to be a good "looks like a known antibiotic" proxy).

**Ranks are recomputed within the subset** (user-directed): the row mask is applied before scaling,
so each endpoint's percentile ranks are relative to the 1.32M surviving compounds, not inherited from
the full library — the counterfactual "what if the library had never contained these compounds",
not "which of the full library's top-N are antibiotic-like".

**The taxonomy-based row order is reused, not recomputed** — it depends only on the fixed 15-organism
list and `config/organism_taxonomy.csv`, neither of which changes under the mask, so
`row_orders["phylo"]` from the main analysis is applied here directly, and the dendrogram diagnostic
is not redrawn (it would be identical).

**Ends with a cell-by-cell comparison against the main (unfiltered) matrix**, reported rather than
asserted: with 97.6% of the library retained the honest expectation is "barely moved", and a large
shift is more likely a masking bug than a finding — but which it is, is a call for a human to make
against the numbers.

Outputs
-------
    output/10_auroc_matrix/10_organism_scores.parquet    the 1.35M x 15 aggregate score matrix
    output/10_auroc_matrix/10_merged_predictors.csv      the 2 merged rank-sum predictor columns
    output/10_auroc_matrix/10_row_order_comparison.csv   baseline/phylo position per organism
    output/10_auroc_matrix/10_auroc_matrix_phylo.csv
    output/10_auroc_matrix/10_auroc_matrix_axes_phylo.csv
    output/10_auroc_matrix/10_overlap_matrix_top1000_phylo.csv
    output/10_auroc_matrix/png|pdf/10_auroc_matrix_phylo.*
    output/10_auroc_matrix/png|pdf/10_overlap_matrix_top1000_phylo.*
    output/10_auroc_matrix/png|pdf/10_phylo_dendrogram.*    diagnostic only, not in figure_cells.json
    output/10_auroc_matrix/10_nonabx_filter_summary.csv       per-filter pass counts
    output/10_auroc_matrix/10_nonabx_merged_predictors.csv
    output/10_auroc_matrix/10_nonabx_auroc_matrix_phylo.csv
    output/10_auroc_matrix/10_nonabx_auroc_matrix_axes_phylo.csv
    output/10_auroc_matrix/10_nonabx_overlap_matrix_top1000_phylo.csv
    output/10_auroc_matrix/10_nonabx_auroc_delta_vs_full.csv  cell-by-cell diff from the main matrix
    output/10_auroc_matrix/png|pdf/10_nonabx_auroc_matrix_phylo.*
    output/10_auroc_matrix/png|pdf/10_nonabx_overlap_matrix_top1000_phylo.*
    output/10_auroc_matrix/figure_cells.json             4 entries (2 main + 2 non-abx)
"""

import os
import sys

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))

from eval_auroc_matrix import (  # noqa: E402
    aggregated_matrix, axes_table, bioactivity_order, diagonal_check, merged_predictor_scores,
    organism_order, organism_scores, phylogeny_class_linkages, phylogeny_organism_order,
    predictor_order, reorder_bioactivity_axes,
)
from eval_auroc_matrix import overlap_matrix  # noqa: E402
from default import ACTIVITY_BINARIZE_TOP_N  # noqa: E402
from eval_predictor_performance import _assert_key_alignment  # noqa: E402
from plots_auroc_matrix import (  # noqa: E402
    save_auroc_matrix_figure, save_dendrogram_figure, save_overlap_matrix_figure,
)

#: The non-abx filter, as (parquet column, operator, threshold). Thresholds are user-directed; the
#: 0.5 on `abx_score` is the only non-zero one and is the single judgement call in the set. Kept
#: unchanged from the former standalone `xx_non_abx_matrix.py`. `eos19mt`'s 38 antibiotic-class flags
#: and `eos6ojg`'s other 9 similarity-count columns were considered and deliberately left out — see
#: scripts/README.md's "Cutoffs reviewed and kept on 2026-08-11".
ABX_FILTERS = [
    ("eos11sm:abx_score", "lt", 0.5),
    ("eos2xeq:is_sim_known_ab", "eq", 0),
    ("eos2xeq:nitrofuran_motif", "eq", 0),
    ("eos2xeq:fluoroquinolone_motif", "eq", 0),
    ("eos2xeq:carbepenem_motif", "eq", 0),
    ("eos2xeq:betalactam_motif", "eq", 0),
    ("eos6ojg:num_sim_0_5_all", "eq", 0),
    ("eos6ojg:num_sim_0_5_subset", "eq", 0),
]

config_dir = os.path.join(root, "..", "config")
output_dir = os.path.join(root, "..", "output", "10_auroc_matrix")
os.makedirs(output_dir, exist_ok=True)

selection_csv = os.path.join(config_dir, "08_endpoint_selection.csv")
pathogens_csv = os.path.join(config_dir, "pathogens_of_interest.csv")
taxonomy_csv = os.path.join(config_dir, "organism_taxonomy.csv")
parquet_path = os.path.join(
    root, "..", "output", "07_score_matrices", "07_score_matrix_full.parquet")
#: Raw per-family matrices, from step 08's cache — inputs to the rank-sum merge below, not read
#: again after. Reading them from step 08 directly (rather than from steps 12/13, which only draw
#: figures now) means this step depends on nothing but 07 and 08 — no figure-drawing script has to
#: run first. See scripts/README.md for why that matters (it is what breaks a real circular
#: dependency this pipeline used to have).
property_matrices_dir = os.path.join(root, "..", "output", "08_property_matrices")
abx_csv = os.path.join(property_matrices_dir, "08_abx_matrix_named.csv")
cytotox_csv = os.path.join(property_matrices_dir, "08_cytotox_matrix_named.csv")
raw_property_csvs = [abx_csv, cytotox_csv]
#: The abx family's own parquet cache — the non-abx section reads the 8 ABX_FILTERS columns from
#: this rather than from `abx_csv`, since it needs positional row-masking, not a name-indexed read.
abx_parquet_path = os.path.join(property_matrices_dir, "08_abx_matrix_full.parquet")

for path, step in [(parquet_path, "07_score_matrices.py"),
                   (abx_csv, "08_property_matrices.py"),
                   (cytotox_csv, "08_property_matrices.py"),
                   (abx_parquet_path, "08_property_matrices.py")]:
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

scores.to_parquet(os.path.join(output_dir, "10_organism_scores.parquet"))

# Cytotoxicity's 6 raw columns and abx-resemblance's 3 each collapse to one rank-sum column, so the
# rest of the pipeline reads only this small derived file, never the raw per-family matrices again.
merged = merged_predictor_scores(raw_property_csvs)
merged.to_csv(os.path.join(output_dir, "10_merged_predictors.csv"), index=False)
print(f"  -> 10_merged_predictors.csv ({merged.shape[0]} x {merged.shape[1]})")
property_csvs = [os.path.join(output_dir, "10_merged_predictors.csv")]

cols = predictor_order()
matrix, cols, predictor_top_indices = aggregated_matrix(
    scores, organisms, cols, property_csvs, return_top_indices=True)

diag = diagonal_check(matrix, organisms["organism"].tolist())
print(f"  diagonal (self-pairs): min {diag.min():.6f}, max {diag.max():.6f} "
      f"({len(diag)} cells — must all be 1.0)")
if not np.allclose(diag, 1.0, atol=1e-3):
    sys.exit("Diagonal is not 1.0 before reordering — rows and columns are misaligned; refusing to "
             "proceed.")

# --- Row order: same class grouping as always, ordered by real NCBI-taxonomy lineage ---
# `row_orders["baseline"]` (alphabetical-within-class) is kept as a reference column in
# 10_row_order_comparison.csv only — no heatmap is drawn from it. An alternative hierarchical-
# clustering-on-AUROC-profile order was tried and dropped 2026-09-02 (see scripts/README.md).
class_of = organisms.set_index("organism")["organism_class"]
row_orders = {"baseline": organisms["organism"].tolist()}
phylo_linkages = phylogeny_class_linkages(organisms, taxonomy_csv)
row_orders["phylo"] = phylogeny_organism_order(organisms, taxonomy_csv)

# Diagnostic only (plain matplotlib, not figure_cells.json) — reads the SAME linkages used to build
# row_orders["phylo"] above, so it cannot show a different tree than the one the heatmap was
# reordered by.
save_dendrogram_figure(
    output_dir, phylo_linkages,
    "taxonomic rank distance (of 5: phylum, class, order, family, genus)",
    "10_phylo_dendrogram")

grouped = {}
for o in row_orders["phylo"]:
    grouped.setdefault(class_of[o], []).append(o)
print("  [phylo] order:")
for cls, members in grouped.items():
    print(f"    {cls}: {', '.join(members)}")

m2, org2 = reorder_bioactivity_axes(matrix, organisms, row_orders["phylo"])
diag2 = diagonal_check(m2, org2["organism"].tolist())
if not np.allclose(diag2, 1.0, atol=1e-3):
    sys.exit("[phylo] diagonal is not 1.0 after reordering — refusing to draw the matrix.")

m2.to_csv(os.path.join(output_dir, "10_auroc_matrix_phylo.csv"))
axes_table(org2, cols).to_csv(
    os.path.join(output_dir, "10_auroc_matrix_axes_phylo.csv"), index=False)
print(f"  -> 10_auroc_matrix_phylo.csv ({m2.shape[0]} x {m2.shape[1]})")
footprints = save_auroc_matrix_figure(output_dir, m2, org2, cols, name="10_auroc_matrix_phylo")

comparison = pd.DataFrame({"organism": row_orders["baseline"]})
comparison["organism_class"] = comparison["organism"].map(class_of)
for suffix, order in row_orders.items():
    position = {o: i for i, o in enumerate(order)}
    comparison[f"{suffix}_position"] = comparison["organism"].map(position)
comparison = comparison.sort_values("baseline_position").reset_index(drop=True)
comparison.to_csv(os.path.join(output_dir, "10_row_order_comparison.csv"), index=False)
print(f"  -> 10_row_order_comparison.csv ({len(comparison)} organisms)")

# --- Second view: how many of the row's actives are in the column's top N ---
# Single cutoff, ACTIVITY_BINARIZE_TOP_N = 1000 (overlap_matrix's own default) — the 100 and 10000
# cutoffs this view was compared at were dropped 2026-09-01 (see scripts/README.md).
overlap = overlap_matrix(scores, organisms, cols, property_csvs,
                         predictor_top_indices=predictor_top_indices)
top_n = ACTIVITY_BINARIZE_TOP_N
name = f"10_overlap_matrix_top{top_n}"

# Cheap invariants, both true by construction — a failure means the tops were mismatched. Checked on
# the un-reordered matrix (no output written for it — see the AUROC-matrix block above) before any
# permutation, then again after reordering.
n_bio = len(organisms)
diag = np.diag(overlap.to_numpy()[:, :n_bio])
assert (diag == top_n).all(), f"top-{top_n} overlap diagonal is not {top_n}: {diag}"
bio = overlap.to_numpy()[:, :n_bio]
assert (bio == bio.T).all(), f"top-{top_n} bioactivity block is not symmetric"

ov2, org2 = reorder_bioactivity_axes(overlap, organisms, row_orders["phylo"])
diag2 = np.diag(ov2.to_numpy()[:, :n_bio])
if not (diag2 == top_n).all():
    sys.exit(f"[phylo] top-{top_n} overlap diagonal is not {top_n} after reordering.")
bio2 = ov2.to_numpy()[:, :n_bio]
if not (bio2 == bio2.T).all():
    sys.exit(f"[phylo] top-{top_n} bioactivity block is not symmetric after reordering.")

ov2.to_csv(os.path.join(output_dir, f"{name}_phylo.csv"))
print(f"  -> {name}_phylo.csv ({ov2.shape[0]} x {ov2.shape[1]})")
# continuous_color=True (user-directed, 2026-09-02): a plain linear 0-1000 colour scale, not the
# discrete non-uniform bins (see plots_auroc_matrix.OverlapMatrixPlot) — those looked log-scaled
# since equal visual steps covered wildly unequal value ranges.
footprints = save_overlap_matrix_figure(output_dir, ov2, org2, cols, footprints,
                                        top_n=top_n, name=f"{name}_phylo", continuous_color=True)

# ----------------------------------------------------------------------------- #
# Non-abx robustness check (absorbed from the former standalone                  #
# xx_non_abx_matrix.py, 2026-09-02) — the exact same two plots, on a library     #
# filtered to purge antibiotic-like compounds.                                   #
# ----------------------------------------------------------------------------- #
# Step 10's own main computation above never checks this; on a filtered run it is the difference
# between a correct matrix and a plausible wrong one.
n_total = _assert_key_alignment(raw_property_csvs, parquet_path)

keys_07 = pd.read_parquet(parquet_path, columns=[]).index.to_series().reset_index(drop=True)
keys_abx = pd.read_parquet(abx_parquet_path, columns=[]).index.to_series().reset_index(drop=True)
if not keys_abx.equals(keys_07):
    sys.exit("Key order in 08_abx_matrix_full.parquet differs from the step-07 matrix — the mask "
             "would select the wrong compounds. Rebuild both from the same reference library.")
print(f"\n[non-abx] key order verified across 3 sources ({n_total:,} compounds)")

# --- The mask ---
available = set(pq.ParquetFile(abx_parquet_path).schema.names)
missing = [c for c, _, _ in ABX_FILTERS if c not in available]
if missing:
    sys.exit(f"Filter column(s) absent from {os.path.basename(abx_parquet_path)}: {missing}")

abx = pd.read_parquet(abx_parquet_path, columns=[c for c, _, _ in ABX_FILTERS])
mask = np.ones(len(abx), dtype=bool)
rows = []
for column, op, threshold in ABX_FILTERS:
    v = abx[column].to_numpy(dtype=float)
    # NaN fails BOTH comparisons, so an unscored compound is excluded rather than assumed clean.
    # That is the conservative reading and it is reported, not silently applied.
    keep = v < threshold if op == "lt" else v == threshold
    n_nan = int(np.isnan(v).sum())
    rows.append({"column": column, "operator": "<" if op == "lt" else "==",
                 "threshold": threshold, "n_pass": int(keep.sum()),
                 "pct_pass": round(100 * keep.sum() / len(v), 4), "n_nan": n_nan})
    mask &= keep
del abx

n_kept = int(mask.sum())
filter_summary = pd.DataFrame(rows)
filter_summary.loc[len(filter_summary)] = {
    "column": "ALL (AND)", "operator": "", "threshold": np.nan, "n_pass": n_kept,
    "pct_pass": round(100 * n_kept / len(mask), 4),
    "n_nan": int(filter_summary["n_nan"].max())}
filter_summary.to_csv(os.path.join(output_dir, "10_nonabx_filter_summary.csv"), index=False)

print(f"\n[non-abx] filter: {n_kept:,} of {len(mask):,} compounds retained "
      f"({100 * n_kept / len(mask):.2f}%), {len(mask) - n_kept:,} removed")
print(filter_summary.to_string(index=False))
n_nan_total = int(filter_summary["n_nan"].iloc[:-1].max())
if n_nan_total:
    print(f"[non-abx] NOTE: {n_nan_total} compound(s) have a NaN in a filter column and are "
          "EXCLUDED (NaN fails `== 0`). Conservative; flip the comparison if they should be kept.")

# --- The exact same two plots as above, recomputed on the masked subset ---
scores_nonabx = organism_scores(parquet_path, endpoints, row_mask=mask)
assert len(scores_nonabx) == n_kept, f"scores has {len(scores_nonabx)} rows, expected {n_kept}"

means_nonabx = scores_nonabx.mean()
if not np.allclose(means_nonabx, 0.5, atol=0.02):
    off = means_nonabx[(means_nonabx - 0.5).abs() > 0.02]
    sys.exit(f"[non-abx] Aggregate score column mean(s) far from 0.5 — scaling axis is wrong:\n{off}")
print(f"[non-abx] organism-score column means {means_nonabx.min():.4f}-{means_nonabx.max():.4f} "
      "(must be ~0.5, and this is meaningful only because the ranks were recomputed AFTER masking)")

# Ranks recomputed WITHIN the mask, same rule as every other column here.
merged_nonabx = merged_predictor_scores(raw_property_csvs, row_mask=mask)
merged_nonabx.to_csv(os.path.join(output_dir, "10_nonabx_merged_predictors.csv"), index=False)
print(f"  -> 10_nonabx_merged_predictors.csv ({merged_nonabx.shape[0]} x {merged_nonabx.shape[1]})")
property_csvs_nonabx = [os.path.join(output_dir, "10_nonabx_merged_predictors.csv")]

cols_nonabx = predictor_order()
matrix_nonabx, cols_nonabx, predictor_top_indices_nonabx = aggregated_matrix(
    scores_nonabx, organisms, cols_nonabx, property_csvs_nonabx, row_mask=None,
    return_top_indices=True)

diag_nonabx = diagonal_check(matrix_nonabx, organisms["organism"].tolist())
print(f"  diagonal (self-pairs): min {diag_nonabx.min():.6f}, max {diag_nonabx.max():.6f} "
      f"({len(diag_nonabx)} cells — must all be 1.0)")
if not np.allclose(diag_nonabx, 1.0, atol=1e-3):
    sys.exit("[non-abx] Diagonal is not 1.0 — rows and columns are misaligned; refusing to draw.")

m2_nonabx, org2_nonabx = reorder_bioactivity_axes(matrix_nonabx, organisms, row_orders["phylo"])
diag2_nonabx = diagonal_check(m2_nonabx, org2_nonabx["organism"].tolist())
if not np.allclose(diag2_nonabx, 1.0, atol=1e-3):
    sys.exit("[non-abx] [phylo] diagonal is not 1.0 after reordering — refusing to draw.")

m2_nonabx.to_csv(os.path.join(output_dir, "10_nonabx_auroc_matrix_phylo.csv"))
axes_table(org2_nonabx, cols_nonabx).to_csv(
    os.path.join(output_dir, "10_nonabx_auroc_matrix_axes_phylo.csv"), index=False)
print(f"  -> 10_nonabx_auroc_matrix_phylo.csv ({m2_nonabx.shape[0]} x {m2_nonabx.shape[1]})")
footprints = save_auroc_matrix_figure(output_dir, m2_nonabx, org2_nonabx, cols_nonabx,
                                      name="10_nonabx_auroc_matrix_phylo")

overlap_nonabx = overlap_matrix(scores_nonabx, organisms, cols_nonabx, property_csvs_nonabx,
                                predictor_top_indices=predictor_top_indices_nonabx)
n_bio = len(organisms)
diag_ov = np.diag(overlap_nonabx.to_numpy()[:, :n_bio])
assert (diag_ov == top_n).all(), f"[non-abx] top-{top_n} overlap diagonal is not {top_n}: {diag_ov}"
bio_ov = overlap_nonabx.to_numpy()[:, :n_bio]
assert (bio_ov == bio_ov.T).all(), f"[non-abx] top-{top_n} bioactivity block is not symmetric"

ov2_nonabx, org2_nonabx = reorder_bioactivity_axes(overlap_nonabx, organisms, row_orders["phylo"])
diag2_ov = np.diag(ov2_nonabx.to_numpy()[:, :n_bio])
if not (diag2_ov == top_n).all():
    sys.exit(f"[non-abx] [phylo] top-{top_n} overlap diagonal is not {top_n} after reordering.")
bio2_ov = ov2_nonabx.to_numpy()[:, :n_bio]
if not (bio2_ov == bio2_ov.T).all():
    sys.exit(f"[non-abx] [phylo] top-{top_n} bioactivity block is not symmetric after reordering.")

nonabx_overlap_name = f"10_nonabx_overlap_matrix_top{top_n}"
ov2_nonabx.to_csv(os.path.join(output_dir, f"{nonabx_overlap_name}_phylo.csv"))
print(f"  -> {nonabx_overlap_name}_phylo.csv ({ov2_nonabx.shape[0]} x {ov2_nonabx.shape[1]})")
footprints = save_overlap_matrix_figure(output_dir, ov2_nonabx, org2_nonabx, cols_nonabx, footprints,
                                        top_n=top_n, name=f"{nonabx_overlap_name}_phylo",
                                        continuous_color=True)

# --- How far did it actually move? Reported, not asserted — see the module docstring. ---
shared_cols = [c for c in m2_nonabx.columns if c in m2.columns]
delta = m2_nonabx[shared_cols] - m2.loc[m2_nonabx.index, shared_cols]
delta.to_csv(os.path.join(output_dir, "10_nonabx_auroc_delta_vs_full.csv"))
a = delta.to_numpy(dtype=float)
print(f"\n[non-abx] vs the main matrix over {len(shared_cols)} shared columns "
      f"({a.size} cells): max |delta| {np.nanmax(np.abs(a)):.4f}, "
      f"median |delta| {np.nanmedian(np.abs(a)):.4f}")
worst = delta.abs().stack().sort_values(ascending=False).head(5)
print("  largest shifts:")
for (r, c), v in worst.items():
    print(f"    {r:<28} {c:<38} {m2.loc[r, c]:.4f} -> {m2_nonabx.loc[r, c]:.4f}  ({v:+.4f})")
dropped_cols = [c for c in m2.columns if c not in m2_nonabx.columns]
print(f"  columns absent here (constant under the filter): {dropped_cols}")

print(f"\nDone -> {output_dir}")
