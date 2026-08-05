"""Step 10 — reference-library chemical-space projection, coloured by pathogen activity.

``eos1klk`` (2D projector, Task=Representation/Projection) computes four 2D layouts — PCA, UMAP,
t-SNE, TMAP — of the same ~1.35M-compound reference library staged by ``00_download_data.py``.
This script draws each projection as a silver full-library density background and highlights, in
crimson, each of the 15 pathogens' (``config/pathogens_of_interest.csv``) ``PROJECTION_TOP_N``
highest-scoring compounds by ``consensus_score`` (predictions already staged by
``00_download_data.py`` into ``data/processed/annotation_preds_ref_library/``) — a rank cutoff,
never a score threshold — producing one small-multiples figure per method.

See ``src/eval_projection.py`` for the memory approach: only one pathogen's two score columns are
ever read at a time, reduced immediately to its top-N rows.

    python 10_reference_library_projection.py

Outputs
-------
    output/10_reference_library_projection/10_{method}_background.csv
    output/10_reference_library_projection/10_top{PROJECTION_TOP_N}_per_pathogen.csv
    output/10_reference_library_projection/png/10_{method}_top{PROJECTION_TOP_N}_pathogens.png
    output/10_reference_library_projection/pdf/10_{method}_top{PROJECTION_TOP_N}_pathogens.pdf
    output/10_reference_library_projection/figure_cells.json
"""

import os
import sys

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))

from default import PROJECTION_MODEL_ID, PROJECTION_PREDS_SUBDIR  # noqa: E402
from eval_projection import run_all  # noqa: E402
from plots_projection import save_projection_figures  # noqa: E402

config_dir = os.path.join(root, "..", "config")
projection_file = os.path.join(
    root, "..", "data", "processed", PROJECTION_PREDS_SUBDIR, f"{PROJECTION_MODEL_ID}_v1.csv")
pred_dir = os.path.join(root, "..", "data", "processed", "annotation_preds_ref_library")
pathogens_csv = os.path.join(config_dir, "pathogens_of_interest.csv")
output_dir = os.path.join(root, "..", "output", "10_reference_library_projection")
os.makedirs(output_dir, exist_ok=True)

if not os.path.exists(projection_file):
    sys.exit(
        f"Missing {projection_file}. Run `python 00_download_data.py` first "
        f"(it fetches {PROJECTION_MODEL_ID} explicitly in Section 4)."
    )

run_all(projection_file=projection_file, pred_dir=pred_dir, pathogens_csv=pathogens_csv,
        output_dir=output_dir)
save_projection_figures(output_dir, pathogens_csv)
print(f"\nDone -> {output_dir}")
