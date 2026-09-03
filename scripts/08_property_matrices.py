"""Step 08 — the physchem, abx and cytotox full-library score matrices, cached.

Step 10 and steps 12-14 downstream all need one or more of these three raw score matrices, and building them —
reading their raw prediction CSVs (physchem: 1 model; abx: 4 models; cytotox: 4 models) over the full
~1.35M-compound reference library — used to be duplicated across steps 10, 11 and 12, each doing its
own read on its own run. That is *why* those three steps could not run without step 09's chemical-
space background already existing (they built their matrix and drew a background-dependent figure in
one script), which in turn could not run without step 14, which needed 11 and 12's matrices — a real
circular dependency (see ``scripts/README.md``).

This step pulls the matrix-building — and nothing else — out of steps 10-12, the same way step 07
already isolates the bioactivity/pathogen matrix from everything downstream of it. Nothing here reads
a figure, a background grid, or another step's output; only ``data/processed/annotation_preds_ref_library/``
raw prediction files and the three per-family selection configs. Every step from 09 onward is now a
pure consumer of a cache — either this one or step 07's.

Three families, three configs, ``selected == "Yes"`` rows only:

    physchem   config/physchem_models.csv         (22 columns, eos4djh — deterministic RDKit arithmetic)
    abx        config/antibiotic_resemblance.csv  (55 columns across 4 models)
    cytotox    config/cytotoxicity_models.csv     (24 columns across 4 models)

**Un-normalized only**, matching steps 10-12's existing convention: no scaling or row normalization —
that is a decision for whoever joins the blocks downstream, not this step. Nothing is dropped or
imputed; missing values are counted in the stats CSVs.

    python 08_property_matrices.py

Outputs
-------
    output/08_property_matrices/08_physchem_matrix_named.csv     (1,355,109 x 22 + key/input)
    output/08_property_matrices/08_physchem_endpoint_stats.csv
    output/08_property_matrices/08_abx_matrix_full.parquet       (cache: all referenced abx columns)
    output/08_property_matrices/08_abx_matrix_named.csv          (1,355,109 x 55)
    output/08_property_matrices/08_abx_endpoint_stats.csv
    output/08_property_matrices/08_cytotox_matrix_named.csv      (1,355,109 x 24 + key/input)
    output/08_property_matrices/08_cytotox_endpoint_stats.csv
"""

import os
import sys
import time

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))

from default import ANNOTATION_PREDS_SUBDIR, PHYSCHEM_PREFIX, TOX_PREFIX  # noqa: E402
from eval_abx_matrix import build_abx_named_matrix  # noqa: E402
from eval_property_matrix import (  # noqa: E402
    build_property_matrix, property_endpoint_stats, report_missing,
)

pred_dir = os.path.join(root, "..", "data", "processed", ANNOTATION_PREDS_SUBDIR)
config_dir = os.path.join(root, "..", "config")
output_dir = os.path.join(root, "..", "output", "08_property_matrices")
os.makedirs(output_dir, exist_ok=True)

# --------------------------------------------------------------------------- #
# 1. Physchem                                                                  #
# --------------------------------------------------------------------------- #
physchem_config = os.path.join(config_dir, "physchem_models.csv")
physchem_named_path = os.path.join(output_dir, "08_physchem_matrix_named.csv")
physchem_stats_path = os.path.join(output_dir, "08_physchem_endpoint_stats.csv")

t0 = time.time()
physchem_matrix, physchem_endpoints = build_property_matrix(pred_dir, physchem_config, PHYSCHEM_PREFIX)
print(f"[property-matrices] physchem matrix ready in {time.time() - t0:.1f}s")

if os.path.exists(physchem_named_path):
    print(f"[skip] {os.path.basename(physchem_named_path)} already exists. Delete it to rebuild.")
else:
    t0 = time.time()
    physchem_matrix.to_csv(physchem_named_path, index=False)
    print(f"[property-matrices] wrote {os.path.basename(physchem_named_path)} "
          f"({os.path.getsize(physchem_named_path) / 1e9:.2f} GB) in {time.time() - t0:.1f}s")

physchem_stats = property_endpoint_stats(physchem_matrix, physchem_endpoints, config_csv=physchem_config)
physchem_stats.to_csv(physchem_stats_path, index=False)
print(f"[property-matrices] {len(physchem_stats)} physchem endpoints -> "
      f"{os.path.basename(physchem_stats_path)}")
report_missing(physchem_matrix, physchem_endpoints)
del physchem_matrix

# --------------------------------------------------------------------------- #
# 2. Antibiotic resemblance (abx)                                              #
# --------------------------------------------------------------------------- #
abx_selection = os.path.join(config_dir, "antibiotic_resemblance.csv")
abx_cache_path = os.path.join(output_dir, "08_abx_matrix_full.parquet")
abx_named_path = os.path.join(output_dir, "08_abx_matrix_named.csv")
abx_stats_path = os.path.join(output_dir, "08_abx_endpoint_stats.csv")

t0 = time.time()
abx_matrix, _ = build_abx_named_matrix(pred_dir, abx_selection, full_matrix_cache_path=abx_cache_path)
print(f"[property-matrices] abx matrix {abx_matrix.shape} ready in {time.time() - t0:.1f}s")

if os.path.exists(abx_named_path):
    print(f"[skip] {os.path.basename(abx_named_path)} already exists. Delete it to rebuild.")
else:
    t0 = time.time()
    abx_matrix.to_csv(abx_named_path)
    print(f"[property-matrices] wrote {os.path.basename(abx_named_path)} "
          f"({os.path.getsize(abx_named_path) / 1e9:.2f} GB) in {time.time() - t0:.1f}s")

abx_stats = property_endpoint_stats(abx_matrix, list(abx_matrix.columns), config_csv=abx_selection)
abx_stats.to_csv(abx_stats_path, index=False)
print(f"[property-matrices] {len(abx_stats)} abx endpoints -> {os.path.basename(abx_stats_path)}")
del abx_matrix

# --------------------------------------------------------------------------- #
# 3. Cytotoxicity                                                              #
# --------------------------------------------------------------------------- #
cytotox_config = os.path.join(config_dir, "cytotoxicity_models.csv")
cytotox_named_path = os.path.join(output_dir, "08_cytotox_matrix_named.csv")
cytotox_stats_path = os.path.join(output_dir, "08_cytotox_endpoint_stats.csv")

t0 = time.time()
cytotox_matrix, cytotox_endpoints = build_property_matrix(pred_dir, cytotox_config, TOX_PREFIX)
print(f"[property-matrices] cytotox matrix ready in {time.time() - t0:.1f}s")

if os.path.exists(cytotox_named_path):
    print(f"[skip] {os.path.basename(cytotox_named_path)} already exists. Delete it to rebuild.")
else:
    t0 = time.time()
    cytotox_matrix.to_csv(cytotox_named_path, index=False)
    print(f"[property-matrices] wrote {os.path.basename(cytotox_named_path)} "
          f"({os.path.getsize(cytotox_named_path) / 1e9:.2f} GB) in {time.time() - t0:.1f}s")

cytotox_stats = property_endpoint_stats(cytotox_matrix, cytotox_endpoints, config_csv=cytotox_config)
cytotox_stats.to_csv(cytotox_stats_path, index=False)
print(f"[property-matrices] {len(cytotox_stats)} cytotox endpoints -> "
      f"{os.path.basename(cytotox_stats_path)}")
report_missing(cytotox_matrix, cytotox_endpoints)
del cytotox_matrix

print(f"\nDone -> {output_dir}")
