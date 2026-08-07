"""Step 12 — the cytotoxicity and safety-endpoint matrix.

Builds the ``cytotox`` block of the reference library from ``config/cytotoxicity_models.csv``
(``selected == "Yes"`` rows: 24 endpoints across 4 models), one row per compound, columns named
``cytotox__{model_id}__{column_name}`` — the same three-part shape as the pathogen matrix (step 07,
prefix = pathogen code), the physchem block (step 10) and the abx block (step 11), so the families can
be joined column-wise on ``key`` later.

Split out of the former ``11_additional_properties.py`` (2026-08-06), which built the physchem and
cytotox blocks in one table. Each family now constructs its own matrix.

The four source models:
    eos42ez  antibiotics-ai-cytotox         (HepG2 / HSkMC / IMR90 cytotoxicity)
    eos7m30  admet-ai-exact                 (Tox21 panels, AMES, DILI, ClinTox, carcinogenicity, LD50)
    eos3le9  hepg2-mmv                      (HepG2 IC50 at 5 uM and 10 uM)
    eos3dys  coadd-antimicrobial-activity   (cytotoxicity IC50, haemolytic activity)

**Every column here is a prediction**, unlike the physchem block (step 10), which is deterministic
RDKit arithmetic. The two must not be treated interchangeably.

**`eos7m30`'s own 8 physicochemical columns stay `No`** in the config: ``eos4djh`` is the single
physchem source, so no two blocks can disagree about molecular weight or logP.

**Carrying the model ID matters here.** Two different models score HepG2 —
``cytotox__eos42ez__cytotoxicity_hepg2`` and ``cytotox__eos3le9__ic50_hepg2_72h_5um`` — so the column
name is what keeps them distinguishable.

**Endpoint selection is not made here** — it is read from the ``selected`` column of the config, which
is manually curated. No threshold or cutoff is applied to any score; raw model outputs are carried
through unchanged, and nothing is dropped or imputed.

Toxicity projection (merged in from the former step 13, 2026-08-06)
-------------------------------------------------------------------
Also draws the reference-library chemical-space projection coloured by predicted toxicity: a silver
full-library density background with, in crimson, each endpoint's ``PROJECTION_TOP_N`` most toxic
compounds, as one small-multiples figure with one panel per endpoint. The pathogen counterpart is
step 09, on the same ``eos1klk`` layout. Merged here because it consumes exactly this matrix and
nothing else; as a separate step it re-read the 720 MB CSV this script has just written.

**The background grid is step 09's, reused not recomputed** (2026-08-07). All three property
families now sit on ONE full-library density grid — pathogens (step 09), abx (step 11) and toxicity
(here) — so their panels are directly comparable. This step used to compute its own
``12_umap_background.csv``; it came out byte-identical to ``09_umap_background.csv``, so the
duplicate bought nothing but a second chance to drift. This step now REQUIRES step 09 to have run.

**A rank cutoff, never a score threshold** — the top ``PROJECTION_TOP_N`` compounds per endpoint by
count, with no score value chosen or reviewed.

**"Most toxic" is the HIGHEST value for all 24 endpoints** — the predicted probability of the toxic
class for the 23 classification endpoints, and for the one regression endpoint ``ld50_zhu`` the
convention its own model uses. See ``src/eval_tox_projection.py``.

**Only the UMAP layout is drawn** (``default.TOX_PROJECTION_METHOD``): 24 endpoints x 4 methods would
be 96 panels, and UMAP is the layout step 09's pathogen figures are read on.

    python 12_cytotox_matrix.py

Outputs
-------
    output/12_cytotox_matrix/12_cytotox_matrix_named.csv      (1,355,109 x 24 + key/input)
    output/12_cytotox_matrix/12_cytotox_endpoint_stats.csv
    output/12_cytotox_matrix/12_top{PROJECTION_TOP_N}_per_endpoint.csv
    output/12_cytotox_matrix/{png,pdf}/12_umap_top{PROJECTION_TOP_N}_toxicity.{png,pdf}
    output/12_cytotox_matrix/figure_cells.json
"""

import os
import sys
import time

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))

from default import (  # noqa: E402
    ANNOTATION_PREDS_SUBDIR, PROJECTION_MODEL_ID, PROJECTION_PREDS_SUBDIR,
    TOX_PROJECTION_METHOD,
)
from eval_property_matrix import (  # noqa: E402
    build_property_matrix, property_endpoint_stats, report_missing,
)
from eval_tox_projection import run_all as run_tox_projection  # noqa: E402
from plots_tox_projection import save_tox_projection_figures  # noqa: E402

pred_dir = os.path.join(root, "..", "data", "processed", ANNOTATION_PREDS_SUBDIR)
projection_file = os.path.join(
    root, "..", "data", "processed", PROJECTION_PREDS_SUBDIR, f"{PROJECTION_MODEL_ID}_v1.csv")
config_dir = os.path.join(root, "..", "config")
output_dir = os.path.join(root, "..", "output", "12_cytotox_matrix")
#: Step 09's full-library density grid, reused rather than recomputed — see section 3.
projection_output_dir = os.path.join(root, "..", "output", "09_reference_library_projection")
background_path = os.path.join(
    projection_output_dir, f"09_{TOX_PROJECTION_METHOD}_background.csv")
os.makedirs(output_dir, exist_ok=True)

PREFIX = "cytotox"
config_csv = os.path.join(config_dir, "cytotoxicity_models.csv")
named_path = os.path.join(output_dir, "12_cytotox_matrix_named.csv")
stats_path = os.path.join(output_dir, "12_cytotox_endpoint_stats.csv")

# --------------------------------------------------------------------------- #
# 1. The named matrix (un-normalized)                                          #
# --------------------------------------------------------------------------- #
t0 = time.time()
matrix, endpoints = build_property_matrix(pred_dir, config_csv, PREFIX)
print(f"[cytotox-matrix] matrix ready in {time.time() - t0:.1f}s")

if os.path.exists(named_path):
    print(f"[skip] {os.path.basename(named_path)} already exists. Delete it to rebuild.")
else:
    t0 = time.time()
    matrix.to_csv(named_path, index=False)
    print(f"[cytotox-matrix] wrote {os.path.basename(named_path)} "
          f"({os.path.getsize(named_path) / 1e9:.2f} GB) in {time.time() - t0:.1f}s")

# --------------------------------------------------------------------------- #
# 2. Per-endpoint stats over the full library                                  #
# --------------------------------------------------------------------------- #
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

# --------------------------------------------------------------------------- #
# 3. Toxicity projection (former step 13)                                      #
# --------------------------------------------------------------------------- #
if not os.path.exists(projection_file):
    sys.exit(
        f"Missing {projection_file}. Run `python 00_download_data.py` first "
        f"(it fetches {PROJECTION_MODEL_ID} explicitly in Section 4)."
    )

# Reads back the matrix CSV this script just wrote rather than the in-memory frame: the projection
# streams it in chunks (key + the 24 score columns only), which is what keeps a 1.35M x 24 top-N
# scan off the peak memory of this process.
#
# The background grid comes from step 09, not from here. All three property families now sit on ONE
# density grid — pathogens (09), abx (11) and toxicity (12) — so their panels are directly
# comparable. Recomputing it here produced a byte-identical file, so the duplicate bought nothing
# but a second chance to drift.
if not os.path.exists(background_path):
    sys.exit(f"Missing {background_path}. Run `python 09_reference_library_projection.py` first — "
             "step 12 reuses its background grid rather than recomputing it.")

print(f"\n[cytotox-matrix] toxicity projection over {len(endpoints)} endpoints...")
tox_endpoints = run_tox_projection(projection_file=projection_file, properties_csv=named_path,
                                   config_csv=config_csv, output_dir=output_dir,
                                   background_path=background_path)
save_tox_projection_figures(output_dir, tox_endpoints, background_path=background_path)

print(f"\nDone -> {output_dir}")
