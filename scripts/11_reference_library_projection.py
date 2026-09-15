"""Step 11 — reference-library chemical-space projection, coloured by pathogen activity.

``eos1klk`` (2D projector, Task=Representation/Projection) computes four 2D layouts — PCA, UMAP,
t-SNE, TMAP — of the same ~1.35M-compound reference library staged by ``00_download_data.py``.
This script draws three figure families on those layouts, all over the same silver full-library
density background and all using a rank cutoff, never a score threshold:

1. **Per pathogen (PCA / t-SNE / TMAP).** Each of the 15 pathogens in
   ``config/pathogens_of_interest.csv`` highlighted in crimson at its ``PROJECTION_TOP_N``
   highest-scoring compounds by ``consensus_score``, from that pathogen's own ChEMBL-derived
   annotation model (predictions staged into
   ``data/processed/annotation_preds_ref_library/``).

2. **Per organism, UMAP only, matched to step 10's AUROC matrix.** Instead of this script's own
   alphabetical order and per-pathogen ``consensus_score`` ranking, the UMAP panel reads step 10's
   OWN top-N compounds and phylogeny-within-class row order directly from its output
   (``output/10_auroc_matrix/10_organism_scores.parquet`` and ``10_row_order_comparison.csv``), so
   this figure and the AUROC matrix show the exact same 15-organism order and the exact same
   highlighted compounds. **Optional, not fatal, if ``10_auroc_matrix.py`` has not run yet** — the
   UMAP panel then falls back to this script's own per-pathogen ``consensus_score`` ranking (like
   the PCA/t-SNE/TMAP panels) instead of the AUROC matrix's matched order, with a printed note. In
   the current numbering this is a non-issue in practice: step 10 runs immediately before this one,
   so the matched panel is available on the very first pass. See
   ``src/plots_projection.AurocMatchedUmapGridPlot`` for the figure and
   ``src/eval_projection.auroc_matched_top_n_per_organism`` for the selection.

3. **Per organism, from ``COADD_MODEL_ID`` (eos3dys), UMAP only.** The same view driven by one
   independent CoAdd-trained model instead of the 15 ChEMBL ones, over the 9 organisms it covers.
   eos3dys publishes no ``consensus_score`` — its 22 outputs are independent per-strain, per-assay
   endpoints — so here the top ``PROJECTION_TOP_N`` is ranked **per endpoint** by that endpoint's
   own score. The 20 bioactivity endpoints (``assay_type == "bioactivity"`` in
   ``config/08_endpoint_selection.csv``; its two ``Homo sapiens`` toxicity columns are excluded and
   belong to step 13) are grouped onto their organism, so an organism with several endpoints shows
   them in one panel as shades of one hue.

See ``src/eval_projection.py`` for the memory approach: prediction columns are read
``usecols``-restricted and in chunks, reduced immediately to top-N rows, and the SMILES ``input``
column is never read.

    python 10_auroc_matrix.py                    # optional, only upgrades the AUROC-matched UMAP panel
    python 11_reference_library_projection.py

Outputs
-------
    output/11_reference_library_projection/11_{method}_background.csv
    output/11_reference_library_projection/11_top{PROJECTION_TOP_N}_per_pathogen.csv
    output/11_reference_library_projection/11_coadd_top{PROJECTION_TOP_N}_per_endpoint.csv
    output/11_reference_library_projection/png/11_{method}_top{PROJECTION_TOP_N}_pathogens.png
    output/11_reference_library_projection/pdf/11_{method}_top{PROJECTION_TOP_N}_pathogens.pdf
    output/11_reference_library_projection/png/11_umap_coadd_top{PROJECTION_TOP_N}_pathogens.png
    output/11_reference_library_projection/pdf/11_umap_coadd_top{PROJECTION_TOP_N}_pathogens.pdf
    output/11_reference_library_projection/figure_cells.json
"""

import os
import sys

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))

from default import (COADD_MODEL_ID, PROJECTION_MODEL_ID,  # noqa: E402
                     PROJECTION_PREDS_SUBDIR)
from eval_projection import run_all, run_coadd, run_auroc_matched_umap  # noqa: E402
from plots_projection import (save_coadd_projection_figures,  # noqa: E402
                              save_projection_figures)

config_dir = os.path.join(root, "..", "config")
projection_file = os.path.join(
    root, "..", "data", "processed", PROJECTION_PREDS_SUBDIR, f"{PROJECTION_MODEL_ID}_v1.csv")
pred_dir = os.path.join(root, "..", "data", "processed", "annotation_preds_ref_library")
coadd_pred_file = os.path.join(pred_dir, f"{COADD_MODEL_ID}_v1.csv")
pathogens_csv = os.path.join(config_dir, "pathogens_of_interest.csv")
selection_csv = os.path.join(config_dir, "08_endpoint_selection.csv")
output_dir = os.path.join(root, "..", "output", "11_reference_library_projection")
os.makedirs(output_dir, exist_ok=True)

# The UMAP panel is redrawn to match step 10's AUROC matrix exactly (see the module docstring), so
# it reads step 10's OWN output rather than recomputing the merge here.
auroc_matrix_output_dir = os.path.join(root, "..", "output", "10_auroc_matrix")
auroc_matched_scores_parquet = os.path.join(auroc_matrix_output_dir, "10_organism_scores.parquet")
auroc_matched_row_order_csv = os.path.join(auroc_matrix_output_dir, "10_row_order_comparison.csv")

if not os.path.exists(projection_file):
    sys.exit(
        f"Missing {projection_file}. Run `python 00_download_data.py` first "
        f"(it fetches {PROJECTION_MODEL_ID} explicitly in Section 4)."
    )

if not os.path.exists(coadd_pred_file):
    sys.exit(
        f"Missing {coadd_pred_file}. Run `python 00_download_data.py` first — {COADD_MODEL_ID} is "
        "staged by the automatic Section 4 Isaura loop (Task=Annotation, Status=Ready)."
    )

# run_all returns the loaded ~1.35M-row projection table so the eos3dys family reuses it rather
# than reading the same file a second time.
proj = run_all(projection_file=projection_file, pred_dir=pred_dir, pathogens_csv=pathogens_csv,
               output_dir=output_dir)

# The AUROC-matched UMAP panel is the ONLY output of this script that needs step 10 — everything
# above (background grids, per-pathogen PCA/t-SNE/TMAP panels) does not. So this is checked here,
# just before the one call that needs it, rather than at the top of the script: a missing step 10
# skips this one panel rather than the whole run (in practice step 10 always runs first, so this
# is just a safety fallback, not the expected path).
auroc_matched_table = None
missing_auroc_matched = [p for p in (auroc_matched_scores_parquet, auroc_matched_row_order_csv)
                         if not os.path.exists(p)]
if missing_auroc_matched:
    print(f"\n[projection] step 10's output is missing ({', '.join(missing_auroc_matched)}) — the "
          "UMAP panel falls back to this script's own per-pathogen consensus_score ranking instead "
          "of the AUROC matrix's matched order. Run `python 10_auroc_matrix.py` then rerun this "
          "script for the matched panel; every other output above is unaffected either way.")
else:
    auroc_matched_table = run_auroc_matched_umap(proj, auroc_matched_scores_parquet,
                                                 auroc_matched_row_order_csv)
save_projection_figures(output_dir, pathogens_csv, auroc_matched_table=auroc_matched_table)

coadd_endpoints_df = run_coadd(proj, coadd_pred_file, selection_csv, output_dir)
save_coadd_projection_figures(output_dir, coadd_endpoints_df)
print(f"\nDone -> {output_dir}")
