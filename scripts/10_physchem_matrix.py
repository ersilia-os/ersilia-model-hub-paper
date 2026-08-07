"""Step 10 — the basic physicochemical descriptor matrix, and each descriptor's distribution.

Builds the ``physchem`` block of the reference library from ``config/physchem_models.csv``
(``selected == "Yes"`` rows: all 22 descriptors of ``eos4djh``), one row per compound, columns named
``physchem__{model_id}__{column_name}`` — the same three-part shape as the pathogen matrix (step 07,
prefix = pathogen code), the abx block (step 11) and the cytotox block (step 12), so the families can
be joined column-wise on ``key`` later.

Split out of the former ``11_additional_properties.py`` (2026-08-06), which built the physchem and
cytotox blocks in one table. Each family now constructs its own matrix and its own figures.

**These are calculations, not predictions.** ``eos4djh`` (datamol-basic-descriptors) wraps
deterministic RDKit-family arithmetic, so unlike every other property block there is no model error,
no training set and therefore no leakage dimension. A consumer must not treat these columns
interchangeably with the predicted ones, and a caption must not describe them as "predicted".

**All 22 descriptors are kept even though 5 are near-redundant** (signed off, not an oversight).
``n_rings``, ``n_aliphatic_rings``, ``n_aromatic_rings`` and ``n_saturated_rings`` are each **exactly**
the sum of their two carbocycle/heterocycle components on 100% of rows, and ``n_radical_electrons`` is
0 for all but ~1 compound in 5,000. Anything fitting a model or reading a correlation matrix off this
block should drop those five first — the independent set is the other 17.
(``n_saturated_rings`` is *not* a duplicate of ``n_aliphatic_rings``; they differ on 34% of rows.)

**Observed ranges are the library's filter, not the model's.** MW 198.9-699.5, cLogP -2.0-7.0,
``n_rings`` <= 6 — the Ersilia reference library is already drug-like-filtered, so these bounds
describe the input set. A caption must not read them as a property of ``eos4djh``.

**Upstream spelling kept:** ``n_aliphatic_heterocyles``, ``n_aromatic_heterocyles`` and
``n_saturated_heterocyles`` are misspelled in Datamol itself. The config matches the real CSV header
rather than correcting it, or the lookup would miss.

**Un-normalized only**, matching the abx and cytotox blocks: choosing a transform is a decision for
after the blocks are joined, and the endpoint stats written here are what that decision should be made
from. Nothing is dropped and no value is imputed.

    python 10_physchem_matrix.py

Outputs
-------
    output/10_physchem_matrix/10_physchem_matrix_named.csv     (1,355,109 x 22 + key/input)
    output/10_physchem_matrix/10_physchem_endpoint_stats.csv
    output/10_physchem_matrix/{png,pdf}/10_physchem_distributions.{png,pdf}
    output/10_physchem_matrix/figure_cells.json
"""

import os
import sys
import time

import pandas as pd

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))

from default import (  # noqa: E402
    ANNOTATION_PREDS_SUBDIR, PHYSCHEM_PREFIX, PHYSCHEM_PROJECTION_ENDPOINTS,
    PROJECTION_MODEL_ID, PROJECTION_PREDS_SUBDIR, PROJECTION_TOP_N, TOX_PROJECTION_METHOD,
)
from eval_property_matrix import (  # noqa: E402
    build_property_matrix, property_endpoint_stats, report_missing,
)
# Generic top-N helpers; they live in eval_tox_projection only because step 12 needed them first.
from eval_projection import load_projection  # noqa: E402
from eval_tox_projection import attach_coordinates, endpoint_top_n  # noqa: E402
from plots_property_matrix import (  # noqa: E402
    save_physchem_projection_figure, save_property_distribution_figure,
)

pred_dir = os.path.join(root, "..", "data", "processed", ANNOTATION_PREDS_SUBDIR)
projection_file = os.path.join(
    root, "..", "data", "processed", PROJECTION_PREDS_SUBDIR, f"{PROJECTION_MODEL_ID}_v1.csv")
config_dir = os.path.join(root, "..", "config")
output_dir = os.path.join(root, "..", "output", "10_physchem_matrix")
os.makedirs(output_dir, exist_ok=True)
#: Step 09's full-library density grid, shared with steps 11 and 12 rather than recomputed.
projection_output_dir = os.path.join(root, "..", "output", "09_reference_library_projection")
background_path = os.path.join(
    projection_output_dir, f"09_{TOX_PROJECTION_METHOD}_background.csv")
top_n_path = os.path.join(output_dir, f"10_top{PROJECTION_TOP_N}_per_descriptor.csv")

PREFIX = PHYSCHEM_PREFIX
config_csv = os.path.join(config_dir, "physchem_models.csv")
named_path = os.path.join(output_dir, "10_physchem_matrix_named.csv")
stats_path = os.path.join(output_dir, "10_physchem_endpoint_stats.csv")

# --------------------------------------------------------------------------- #
# 1. The named matrix (un-normalized)                                          #
# --------------------------------------------------------------------------- #
t0 = time.time()
matrix, endpoints = build_property_matrix(pred_dir, config_csv, PREFIX)
print(f"[physchem-matrix] matrix ready in {time.time() - t0:.1f}s")

if os.path.exists(named_path):
    print(f"[skip] {os.path.basename(named_path)} already exists. Delete it to rebuild.")
else:
    t0 = time.time()
    matrix.to_csv(named_path, index=False)
    print(f"[physchem-matrix] wrote {os.path.basename(named_path)} "
          f"({os.path.getsize(named_path) / 1e9:.2f} GB) in {time.time() - t0:.1f}s")

# --------------------------------------------------------------------------- #
# 2. Per-endpoint stats over the full library                                  #
# --------------------------------------------------------------------------- #
stats = property_endpoint_stats(matrix, endpoints, config_csv=config_csv)
stats.to_csv(stats_path, index=False)
print(f"\n[physchem-matrix] {len(stats)} endpoints over {len(matrix):,} compounds "
      f"-> {os.path.basename(stats_path)}")

# Reported, never acted on: a constant column cannot be z-scored (std = 0) and has no top-N, so
# whoever picks a transform after the join needs to see it here.
constant = stats[stats["n_unique"] <= 1]
near_constant = stats[(stats["n_unique"] > 1) & (stats["pct_nonzero"] < 1.0)]
print(f"    constant library-wide:        {len(constant)}"
      + (f" -> {', '.join(constant['column_name'])}" if len(constant) else ""))
print(f"    non-zero in <1% of compounds: {len(near_constant)}"
      + (f" -> {', '.join(near_constant['column_name'])}" if len(near_constant) else ""))
print(f"    missing values (kept, never imputed): {int(stats['n_nan'].sum())}")

report_missing(matrix, endpoints)

# --------------------------------------------------------------------------- #
# 3. Figure — one distribution panel per descriptor                            #
# --------------------------------------------------------------------------- #
print(f"\n[physchem-matrix] figure: {len(endpoints)} descriptor distributions")
save_property_distribution_figure(output_dir, matrix, endpoints,
                                 name="10_physchem_distributions")

# --------------------------------------------------------------------------- #
# 4. UMAP panels for three descriptors (user-directed)                         #
# --------------------------------------------------------------------------- #
# Only MW, TPSA and cLogP, on step 09's shared background — the same grid the abx (11) and toxicity
# (12) panels use, so all four families are directly comparable. NOTE these three answer a weaker
# question: the descriptors are continuous and unimodal, so "top 1000" is the extreme TAIL (the
# heaviest / most polar / most lipophilic molecules), not a selected set. See
# default.PHYSCHEM_PROJECTION_ENDPOINTS.
for path, step in [(projection_file, "00_download_data.py"),
                   (background_path, "09_reference_library_projection.py")]:
    if not os.path.exists(path):
        sys.exit(f"Missing {path}. Run `python {step}` first.")

proj_endpoints = pd.DataFrame([
    {"model_id": "eos4djh", "column_name": c, "endpoint": f"{PREFIX}__eos4djh__{c}"}
    for c in PHYSCHEM_PROJECTION_ENDPOINTS
])
missing_cols = [e for e in proj_endpoints["endpoint"] if e not in matrix.columns]
if missing_cols:
    sys.exit(f"PHYSCHEM_PROJECTION_ENDPOINTS names columns absent from the matrix: {missing_cols}")

print(f"\n[physchem-matrix] UMAP panels for {len(proj_endpoints)} descriptors: "
      f"{', '.join(PHYSCHEM_PROJECTION_ENDPOINTS)}")
proj = load_projection(projection_file)
tops = endpoint_top_n(named_path, proj_endpoints["endpoint"].tolist(), n=PROJECTION_TOP_N)
top_table = attach_coordinates(tops, proj, method=TOX_PROJECTION_METHOD)
top_table.to_csv(top_n_path, index=False)
print(f"  -> {os.path.basename(top_n_path)} ({len(top_table):,} rows)")
for r in proj_endpoints.itertuples():
    g = top_table[top_table["endpoint"] == r.endpoint]
    print(f"    {r.column_name}: top {len(g)} (value {g['score'].min():.4g}-{g['score'].max():.4g})")

save_physchem_projection_figure(output_dir, background_path, top_table, proj_endpoints,
                                top_n=PROJECTION_TOP_N, method=TOX_PROJECTION_METHOD)

print(f"\nDone -> {output_dir}")
