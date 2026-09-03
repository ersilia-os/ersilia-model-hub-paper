"""Step 13 — cytotoxicity/safety endpoint stats, and the toxicity projection.

Reads the ``cytotox`` block step 08 already built and cached (``config/cytotoxicity_models.csv``,
``selected == "Yes"`` rows: 24 endpoints across 4 models) — no raw prediction file is read here, and
no matrix is rebuilt. The matrix-building half of this step moved to :mod:`08_property_matrices` so
it can be shared with step 14 without pulling in this step's own dependency on step 11's chemical-
space background (see that step's docstring and ``scripts/README.md`` for why that split matters —
it is what breaks a real circular dependency this pipeline used to have).

The four source models:
    eos42ez  antibiotics-ai-cytotox         (HepG2 / HSkMC / IMR90 cytotoxicity)
    eos7m30  admet-ai-exact                 (Tox21 panels, AMES, DILI, ClinTox, carcinogenicity, LD50)
    eos3le9  hepg2-mmv                      (HepG2 IC50 at 5 uM and 10 uM)
    eos3dys  coadd-antimicrobial-activity   (cytotoxicity IC50, haemolytic activity)

**Every column here is a prediction**, unlike the physchem block, which is deterministic RDKit
arithmetic. The two must not be treated interchangeably.

Toxicity projection (merged in from the former step 13, 2026-08-06)
-------------------------------------------------------------------
Also draws the reference-library chemical-space projection coloured by predicted toxicity: a silver
full-library density background with, in crimson, each endpoint's ``PROJECTION_TOP_N`` most toxic
compounds, as one small-multiples figure with one panel per endpoint. The pathogen counterpart is
step 11, on the same ``eos1klk`` layout.

**The background grid is step 11's, reused not recomputed** (2026-08-07). All three property
families now sit on ONE full-library density grid — pathogens (step 11), abx (step 12) and toxicity
(here) — so their panels are directly comparable. This step now REQUIRES step 11 to have run.

**A rank cutoff, never a score threshold** — the top ``PROJECTION_TOP_N`` compounds per endpoint by
count, with no score value chosen or reviewed.

**"Most toxic" is the HIGHEST value for all 24 endpoints** — the predicted probability of the toxic
class for the 23 classification endpoints, and for the one regression endpoint ``ld50_zhu`` the
convention its own model uses. See ``src/eval_tox_projection.py``.

**Only the UMAP layout is drawn** (``default.TOX_PROJECTION_METHOD``): 24 endpoints x 4 methods would
be 96 panels, and UMAP is the layout step 11's pathogen figures are read on.

    python 08_property_matrices.py
    python 11_reference_library_projection.py
    python 13_toxicity_projection.py

Outputs
-------
    output/13_toxicity_projection/13_cytotox_endpoint_stats.csv
    output/13_toxicity_projection/13_top{PROJECTION_TOP_N}_per_endpoint.csv
    output/13_toxicity_projection/{png,pdf}/13_umap_top{PROJECTION_TOP_N}_toxicity.{png,pdf}
    output/13_toxicity_projection/figure_cells.json
"""

import os
import sys

import pandas as pd

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))

from default import (  # noqa: E402
    PROJECTION_MODEL_ID, PROJECTION_PREDS_SUBDIR, TOX_PREFIX, TOX_PROJECTION_METHOD,
)
from eval_property_matrix import property_endpoint_stats, report_missing  # noqa: E402
from eval_tox_projection import run_all as run_tox_projection  # noqa: E402
from plots_tox_projection import save_tox_projection_figures  # noqa: E402

projection_file = os.path.join(
    root, "..", "data", "processed", PROJECTION_PREDS_SUBDIR, f"{PROJECTION_MODEL_ID}_v1.csv")
config_dir = os.path.join(root, "..", "config")
output_dir = os.path.join(root, "..", "output", "13_toxicity_projection")
property_matrices_dir = os.path.join(root, "..", "output", "08_property_matrices")
#: Step 11's full-library density grid, reused rather than recomputed — see the docstring above.
projection_output_dir = os.path.join(root, "..", "output", "11_reference_library_projection")
background_path = os.path.join(
    projection_output_dir, f"11_{TOX_PROJECTION_METHOD}_background.csv")
os.makedirs(output_dir, exist_ok=True)

PREFIX = TOX_PREFIX
config_csv = os.path.join(config_dir, "cytotoxicity_models.csv")
named_path = os.path.join(property_matrices_dir, "08_cytotox_matrix_named.csv")
stats_path = os.path.join(output_dir, "13_cytotox_endpoint_stats.csv")

for path, step in [(named_path, "08_property_matrices.py"),
                   (projection_file, "00_download_data.py"),
                   (background_path, "11_reference_library_projection.py")]:
    if not os.path.exists(path):
        sys.exit(f"Missing {path}. Run `python {step}` first.")

# --------------------------------------------------------------------------- #
# 1. Per-endpoint stats over the full library, from step 08's cached matrix    #
# --------------------------------------------------------------------------- #
matrix = pd.read_csv(named_path)
endpoints = [c for c in matrix.columns if c.startswith(f"{PREFIX}__")]
print(f"[cytotox-matrix] loaded {os.path.basename(named_path)} "
      f"({len(matrix):,} x {len(endpoints)})")

stats = property_endpoint_stats(matrix, endpoints, config_csv=config_csv)
stats.to_csv(stats_path, index=False)
print(f"\n[cytotox-matrix] {len(stats)} endpoints over {len(matrix):,} compounds "
      f"-> {os.path.basename(stats_path)}")

constant = stats[stats["n_unique"] <= 1]
binary = stats[stats["n_unique"] <= 2]
print(f"    binary (<=2 distinct values):  {len(binary)}")
print(f"    constant library-wide:         {len(constant)}"
      + (f" -> {', '.join(constant['column_name'])}" if len(constant) else ""))
print(f"    missing values (kept, never imputed): {int(stats['n_nan'].sum())}")

report_missing(matrix, endpoints)
del matrix

# --------------------------------------------------------------------------- #
# 2. Toxicity projection (former step 13)                                      #
# --------------------------------------------------------------------------- #
# Reads back step 08's matrix CSV rather than an in-memory frame: the projection streams it in
# chunks (key + the 24 score columns only), which is what keeps a 1.35M x 24 top-N scan off the
# peak memory of this process.
print(f"\n[cytotox-matrix] toxicity projection over {len(endpoints)} endpoints...")
tox_endpoints = run_tox_projection(projection_file=projection_file, properties_csv=named_path,
                                   config_csv=config_csv, output_dir=output_dir,
                                   background_path=background_path)
save_tox_projection_figures(output_dir, tox_endpoints, background_path=background_path)

print(f"\nDone -> {output_dir}")
