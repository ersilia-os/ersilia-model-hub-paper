"""Step 11 — the antibiotic-resemblance score matrix, and its endpoints on the library UMAP.

Companion to step 07's pathogen matrices, built from ``config/antibiotic_resemblance.csv``
(``selected == Yes`` rows: 55 endpoints across 4 models). Same reference library, same ``key``
alignment, so this block can be joined column-wise to the pathogen and cytotoxicity blocks later —
column names are ``abx__{model_id}__{column_name}``, the same three-part shape the pathogen columns
use, with a constant group code in the pathogen slot.

**Un-normalized only.** No scaling and no row normalization are applied here, by request: the
transforms are a decision for after the blocks are joined, and the stats CSV this writes is what
that decision should be made from.

The figure: every endpoint's compounds on eos1klk's UMAP, silver full-library density behind,
crimson highlights on top — reusing step 09's ``09_umap_background.csv`` rather than recomputing it,
so the panels are directly comparable to the pathogen ones.

**The highlight rule is not a rank cutoff.** 54 of the 55 endpoints are binary flags or small
integer counts (only ``abx_score`` is continuous), so each panel shows every compound with a value
> 0, capped at ``PROJECTION_TOP_N``, never padded with zeros. See
:func:`eval_abx_matrix.endpoint_highlights` for why, and read every panel's ``n_shown/n_nonzero``
annotation: where the cap binds on a binary column the kept subset is arbitrary.

Endpoints with zero non-zero compounds library-wide are omitted from the FIGURE only (user-directed
— they have nothing to draw). They remain in the matrix, in the stats CSV, and are printed below.

Pathogen x antibiotic-resemblance overlap (merged in from the former step 14, 2026-08-06)
-----------------------------------------------------------------------------------------
Also overlays the two compound sets this pipeline has already selected, on the same UMAP layout and
background: one figure per pathogen (15), each a 3x3 grid with one panel per antibiotic-resemblance
endpoint (``default.ABX_OVERLAP_ENDPOINTS``, 9), showing the pathogen's hits, the endpoint's
antibiotic-like compounds and — drawn last and larger — their intersection. Merged here because the
abx side of that intersection is the highlights table this script writes a few lines earlier; as a
separate step it re-read it from disk.

**No threshold is chosen at that stage** — it is set arithmetic over two existing selections. But the
two sides are selected DIFFERENTLY and a panel cannot be read without knowing it:

*   Pathogen side: a RANK CUTOFF, the top ``PROJECTION_TOP_N`` by ``consensus_score`` from step 09.
    Every pathogen contributes exactly ``PROJECTION_TOP_N`` compounds.
*   Abx side: the value > 0 rule above, capped at ``PROJECTION_TOP_N``. Three of the nine
    (carbepenem_motif, ansamycins_rifamycins_macrolides, b_lactams_all) have fewer positives than the
    cap and are drawn exhaustively; the other six hit the cap and show an arbitrary subset of their
    positives, marked ``*`` on the panel.

    python 11_abx_resemblance_matrix.py

Outputs
-------
    output/11_abx_resemblance_matrix/11_abx_matrix_full.parquet   (cache: all referenced)
    output/11_abx_resemblance_matrix/11_abx_matrix_named.csv      (1,355,109 x 55)
    output/11_abx_resemblance_matrix/11_abx_endpoint_stats.csv
    output/11_abx_resemblance_matrix/11_abx_umap_highlights.csv
    output/11_abx_resemblance_matrix/11_overlap_points.csv
    output/11_abx_resemblance_matrix/11_overlap_counts.csv
    output/11_abx_resemblance_matrix/png|pdf/11_umap_abx_endpoints_max1000.*
    output/11_abx_resemblance_matrix/png|pdf/11_umap_abx_overlap_{pathogen}.*   (15)
    output/11_abx_resemblance_matrix/figure_cells.json
"""

import os
import sys
import time

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))

from default import (  # noqa: E402
    ANNOTATION_PREDS_SUBDIR, PROJECTION_MODEL_ID, PROJECTION_PREDS_SUBDIR, PROJECTION_TOP_N,
)
from eval_abx_matrix import (  # noqa: E402
    build_abx_named_matrix, endpoint_highlights, endpoint_stats, load_umap,
)
from eval_pathogen_abx_overlap import run_all as run_pathogen_overlap  # noqa: E402
from plots_abx_projection import save_abx_projection_figure  # noqa: E402
from plots_pathogen_abx_overlap import save_overlap_figures  # noqa: E402

pred_dir = os.path.join(root, "..", "data", "processed", ANNOTATION_PREDS_SUBDIR)
projection_file = os.path.join(
    root, "..", "data", "processed", PROJECTION_PREDS_SUBDIR, f"{PROJECTION_MODEL_ID}_v1.csv")
config_dir = os.path.join(root, "..", "config")
output_dir = os.path.join(root, "..", "output", "11_abx_resemblance_matrix")
projection_output_dir = os.path.join(root, "..", "output", "09_reference_library_projection")
os.makedirs(output_dir, exist_ok=True)

selection_path = os.path.join(config_dir, "antibiotic_resemblance.csv")
pathogens_csv = os.path.join(config_dir, "pathogens_of_interest.csv")
pathogen_top_csv = os.path.join(
    projection_output_dir, f"09_top{PROJECTION_TOP_N}_per_pathogen.csv")
cache_path = os.path.join(output_dir, "11_abx_matrix_full.parquet")
named_path = os.path.join(output_dir, "11_abx_matrix_named.csv")
stats_path = os.path.join(output_dir, "11_abx_endpoint_stats.csv")
highlights_path = os.path.join(output_dir, "11_abx_umap_highlights.csv")
background_path = os.path.join(projection_output_dir, "09_umap_background.csv")

if not os.path.exists(projection_file):
    sys.exit(f"Missing {projection_file}. Run `python 00_download_data.py` first "
             f"(it fetches {PROJECTION_MODEL_ID} in Section 4).")

# --------------------------------------------------------------------------- #
# 1. The named matrix (un-normalized)                                          #
# --------------------------------------------------------------------------- #
t0 = time.time()
matrix, sel = build_abx_named_matrix(pred_dir, selection_path, full_matrix_cache_path=cache_path)
print(f"[abx-matrix] named matrix {matrix.shape} ready in {time.time() - t0:.1f}s")

if os.path.exists(named_path):
    print(f"[skip] {os.path.basename(named_path)} already exists. Delete it to rebuild.")
else:
    t0 = time.time()
    matrix.to_csv(named_path)
    print(f"[abx-matrix] wrote {os.path.basename(named_path)} "
          f"({os.path.getsize(named_path) / 1e9:.2f} GB) in {time.time() - t0:.1f}s")

# --------------------------------------------------------------------------- #
# 2. Per-endpoint stats over the full library                                  #
# --------------------------------------------------------------------------- #
stats = endpoint_stats(matrix, sel)

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
# 5. Pathogen hits vs antibiotic-resemblance hits (former step 14)             #
# --------------------------------------------------------------------------- #
# Set arithmetic over two selections this pipeline has already made — no threshold is chosen here.
# The two sides are selected DIFFERENTLY and a panel cannot be read without knowing it: the pathogen
# side is a rank cutoff (exactly PROJECTION_TOP_N per pathogen, by consensus_score, from step 09),
# while the abx side is "every compound with a value > 0, capped at PROJECTION_TOP_N" (section 3
# above), so an abx endpoint can contribute far fewer than the cap, or an arbitrary subset of its
# positives where the cap binds.
if not os.path.exists(pathogen_top_csv):
    sys.exit(f"Missing {pathogen_top_csv}. Run `python 09_reference_library_projection.py` first.")

print(f"\n[abx-matrix] pathogen x abx overlap...")
pathogens, overlap_endpoints, counts = run_pathogen_overlap(
    pathogen_top_csv=pathogen_top_csv, abx_highlights_csv=highlights_path,
    pathogens_csv=pathogens_csv, abx_config_csv=selection_path, output_dir=output_dir)
save_overlap_figures(output_dir, background_path, pathogens, overlap_endpoints, counts)

print(f"\nDone -> {output_dir}")
