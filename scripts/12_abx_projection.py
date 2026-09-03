"""Step 12 — the antibiotic-resemblance endpoints on the library UMAP, and their class enrichment.

Reads the ``abx`` block step 08 already built and cached (``config/antibiotic_resemblance.csv``,
``selected == Yes`` rows: 55 endpoints across 4 models) — no raw prediction file is read here, and no
matrix is rebuilt. This step is the figure/statistics half of what used to be one script; the matrix-
building half moved to :mod:`08_property_matrices` so it can be shared without pulling in this step's
own dependency on step 11's chemical-space background (see that step's docstring and
``scripts/README.md`` for why that split matters — it is what breaks a real circular dependency this
pipeline used to have).

The figure: every endpoint's compounds on eos1klk's UMAP, silver full-library density behind,
crimson highlights on top — reusing step 11's ``11_umap_background.csv`` rather than recomputing it,
so the panels are directly comparable to the pathogen ones.

**The highlight rule is not a rank cutoff.** 54 of the 55 endpoints are binary flags or small
integer counts (only ``abx_score`` is continuous), so each panel shows every compound with a value
> 0, capped at ``PROJECTION_TOP_N``, never padded with zeros. See
:func:`eval_abx_matrix.endpoint_highlights` for why, and read every panel's ``n_shown/n_nonzero``
annotation: where the cap binds on a binary column the kept subset is arbitrary.

Endpoints with zero non-zero compounds library-wide are omitted from the FIGURE only (user-directed
— they have nothing to draw). They remain in the matrix, in the stats CSV, and are printed below.

Pathogen top-1000 x eos19mt antibiotic-class Fisher enrichment (replaces the former Jaccard/UMAP
overlap analysis, 2026-09-01)
--------------------------------------------------------------------------------------------------
For every pathogen and every one of eos19mt's 38 antibiotic structural classes, tests whether being
in the pathogen's top ``PROJECTION_TOP_N`` (by ``consensus_score``, from step 11) is associated with
belonging to that class (value > 0), via a one-sided Fisher exact test over the full reference
library. Scoped to eos19mt only — not the mixed multi-model set the earlier overlap version used —
because "belonging to specific antibiotic classes as defined by eos19mt" is the question asked.

Rows are eos19mt's 38 classes in the model's own (== the config file's) order; columns are the 15
pathogens in step 10's phylogenetic order, reusing :mod:`eval_auroc_matrix`'s ordering functions so
this matrix's columns can never disagree with the AUROC matrix's rows. See
:mod:`eval_abx_enrichment` for the contingency-table definition (missing values excluded per class,
degenerate constant-zero classes report ``NaN``, not a misleading computed value).

    python 08_property_matrices.py
    python 11_reference_library_projection.py
    python 12_abx_projection.py

Outputs
-------
    output/12_abx_projection/12_abx_endpoint_stats.csv
    output/12_abx_projection/12_abx_umap_highlights.csv
    output/12_abx_projection/12_abx_enrichment_long.csv        (570 rows: 15 pathogens x 38 classes)
    output/12_abx_projection/12_abx_enrichment_odds_ratio.csv  (38 x 15 wide matrix)
    output/12_abx_projection/12_abx_enrichment_pvalue.csv      (38 x 15 wide matrix)
    output/12_abx_projection/png|pdf/12_umap_abx_endpoints_max1000.*
    output/12_abx_projection/png|pdf/12_abx_enrichment_matrix.*
    output/12_abx_projection/figure_cells.json
"""

import os
import sys
import time

import pandas as pd

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))

from default import (  # noqa: E402
    PROJECTION_MODEL_ID, PROJECTION_PREDS_SUBDIR, PROJECTION_TOP_N,
)
from eval_abx_enrichment import run_all as run_abx_enrichment  # noqa: E402
from eval_abx_matrix import endpoint_highlights, load_umap  # noqa: E402
from eval_property_matrix import property_endpoint_stats  # noqa: E402
from plots_abx_enrichment import save_enrichment_figure  # noqa: E402
from plots_abx_projection import save_abx_projection_figure  # noqa: E402

config_dir = os.path.join(root, "..", "config")
output_dir = os.path.join(root, "..", "output", "12_abx_projection")
property_matrices_dir = os.path.join(root, "..", "output", "08_property_matrices")
projection_output_dir = os.path.join(root, "..", "output", "11_reference_library_projection")
os.makedirs(output_dir, exist_ok=True)

selection_path = os.path.join(config_dir, "antibiotic_resemblance.csv")
pathogens_csv = os.path.join(config_dir, "pathogens_of_interest.csv")
endpoint_selection_csv = os.path.join(config_dir, "08_endpoint_selection.csv")
taxonomy_csv = os.path.join(config_dir, "organism_taxonomy.csv")
projection_file = os.path.join(
    root, "..", "data", "processed", PROJECTION_PREDS_SUBDIR, f"{PROJECTION_MODEL_ID}_v1.csv")
pathogen_top_csv = os.path.join(
    projection_output_dir, f"11_top{PROJECTION_TOP_N}_per_pathogen.csv")
abx_named_path = os.path.join(property_matrices_dir, "08_abx_matrix_named.csv")
stats_path = os.path.join(output_dir, "12_abx_endpoint_stats.csv")
highlights_path = os.path.join(output_dir, "12_abx_umap_highlights.csv")
background_path = os.path.join(projection_output_dir, "11_umap_background.csv")

for path, step in [(abx_named_path, "08_property_matrices.py"),
                   (pathogen_top_csv, "11_reference_library_projection.py"),
                   (background_path, "11_reference_library_projection.py"),
                   (projection_file, "00_download_data.py")]:
    if not os.path.exists(path):
        sys.exit(f"Missing {path}. Run `python {step}` first.")

# --------------------------------------------------------------------------- #
# 1. The named matrix and its selection metadata, from step 08's cache         #
# --------------------------------------------------------------------------- #
t0 = time.time()
matrix = pd.read_csv(abx_named_path, index_col="key")
sel_all = pd.read_csv(selection_path)
sel = sel_all[sel_all["selected"] == "Yes"].copy()
print(f"[abx-matrix] loaded {os.path.basename(abx_named_path)} {matrix.shape} in "
      f"{time.time() - t0:.1f}s")

# --------------------------------------------------------------------------- #
# 2. Per-endpoint stats over the full library                                  #
# --------------------------------------------------------------------------- #
stats = property_endpoint_stats(matrix, list(matrix.columns), config_csv=selection_path)

# Reported, never acted on: a column that is constant over 1.35M compounds cannot be z-scored
# (std = 0) and has no top-N, so whoever chooses a transform after the join needs to see it here.
constant = stats[stats["n_nonzero"] == 0]
sparse = stats[(stats["n_nonzero"] > 0) & (stats["n_nonzero"] < PROJECTION_TOP_N)]
binary = stats[stats["n_unique"] <= 2]
print(f"\n[abx-matrix] {len(stats)} endpoints over {len(matrix):,} compounds")
print(f"    binary (<=2 distinct values): {len(binary)}")
print(f"    constant zero library-wide:   {len(constant)}"
      + (f" -> {', '.join(constant['column_name'])}" if len(constant) else ""))
print(f"    fewer than {PROJECTION_TOP_N} non-zero:      {len(sparse)}")
total_nan = int(stats["n_nan"].sum())
print(f"    missing values (kept, never imputed): {total_nan}")

# --------------------------------------------------------------------------- #
# 3. UMAP highlights                                                           #
# --------------------------------------------------------------------------- #
t0 = time.time()
proj = load_umap(projection_file)
print(f"\n[abx-matrix] loaded UMAP coordinates for {len(proj):,} compounds "
      f"in {time.time() - t0:.1f}s")

highlights, stats = endpoint_highlights(matrix, stats, proj, top_n=PROJECTION_TOP_N)
stats.to_csv(stats_path, index=False)
highlights.to_csv(highlights_path, index=False)
print(f"[abx-matrix] wrote {os.path.basename(stats_path)} and "
      f"{os.path.basename(highlights_path)} ({len(highlights):,} highlighted rows)")

capped = stats[stats["n_shown"] < stats["n_nonzero"]]
print(f"    {len(capped)} endpoint(s) hit the {PROJECTION_TOP_N} cap — their panels show an "
      "arbitrary (key-ordered) subset, not the full flagged set")

# --------------------------------------------------------------------------- #
# 4. Figure — one panel per endpoint with something to draw                    #
# --------------------------------------------------------------------------- #
drawable = stats[stats["n_shown"] > 0].reset_index(drop=True)
print(f"\n[abx-matrix] figure: {len(drawable)} of {len(stats)} endpoints "
      f"({len(stats) - len(drawable)} omitted for having no non-zero compound)")
save_abx_projection_figure(output_dir, background_path, highlights, drawable,
                           top_n=PROJECTION_TOP_N)

# --------------------------------------------------------------------------- #
# 5. Pathogen top-1000 x eos19mt antibiotic-class Fisher enrichment            #
# --------------------------------------------------------------------------- #
# Reuses the already-loaded named matrix (section 1) and sel (its 55-row endpoint metadata) —
# eos19mt's 38 classes are a subset of the same columns, no re-read of the raw prediction files.
print(f"\n[abx-matrix] pathogen top-1000 x eos19mt Fisher enrichment...")
classes, pathogens, long_df, odds_ratio, p_value = run_abx_enrichment(
    matrix=matrix, sel=sel, pathogen_top_csv=pathogen_top_csv,
    selection_csv=endpoint_selection_csv, pathogens_csv=pathogens_csv, taxonomy_csv=taxonomy_csv,
    output_dir=output_dir)

# Figure only (user-directed): a class with odds_ratio 0 or NaN against every pathogen has nothing
# to show — no enrichment signal and no colour but the neutral/NaN one. Dropped from the FIGURE
# only; all 38 classes stay in the CSVs above, same convention as section 4's "drawable" endpoints.
informative = ~(odds_ratio.fillna(0) == 0).all(axis=1)
omitted = odds_ratio.index[~informative].tolist()
print(f"\n[abx-matrix] figure: {informative.sum()} of {len(odds_ratio)} eos19mt classes have "
      f"signal ({len(omitted)} omitted, all-NaN/zero against every pathogen): {omitted}")
save_enrichment_figure(output_dir, odds_ratio.loc[informative], p_value.loc[informative], pathogens)

print(f"\nDone -> {output_dir}")
