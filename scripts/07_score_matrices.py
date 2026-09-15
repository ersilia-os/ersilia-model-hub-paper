"""Step 07 — the full-library score matrix cache.

The foundation of the correlation analysis. Every annotation model was run on the SAME
~1.35M-compound reference library and staged by ``00_download_data.py`` as
``data/processed/annotation_preds_ref_library/{model_id}_{version}.csv``, all aligned on ``key``, so
the predictions form a clean rectangular matrix with no compound alignment needed.

Reads the raw prediction CSVs and derives the **named matrix** — the raw scores for every
``selected == Yes`` endpoint in ``config/08_endpoint_selection.csv``, columns renamed
``{pathogen_code}__{model_id}__{column_name}`` so pathogen and source model read straight off the
column name. Downstream steps derive their own scaled/row-normalized variants from this in memory
(column z-score or rank-percentile, then row L2 or L1 normalization — see ``eval_correlations.
scale_matrix``/``row_normalize``) rather than reading a pre-materialised copy; no scaled variant is
written to disk by this script, since nothing in the repo reads one that isn't already in memory
(a prior CSV-export flag for this was removed 2026-09-01 after going unused for a month — see
``scripts/README.md``'s "Pipeline reorganization" note).

Reads the ~15 GB of raw prediction CSVs at most **once**, caching the result as
``07_score_matrix_full.parquet`` (every endpoint the selection CSV references, Yes AND No, so
re-enabling one later needs no re-extraction). Every subsequent run, and every downstream step,
re-slices that cache in under a second. Every output is skipped if it already exists, so deleting one
file rebuilds exactly that file.

Mean percentile rank per compound
---------------------------------
Also collapses the rank-percentile matrix along its columns to one number per compound: that
compound's **mean percentile rank** over every selected endpoint. A percentile rank is a compound's
standing within one endpoint over the whole library, so the mean answers "how highly does this
molecule rank, on average, across everything we predict". Merged in from the former standalone step
09 (2026-08-06) — it reads nothing this script does not already have in memory.

The request was phrased as a column-wise sum; the reported quantity is the **mean**, which is what
"average rank" means and what keeps the value on the interpretable [0, 1] percentile scale. The two
differ only by a constant factor — multiply by ``n_endpoints`` in the CSV for the sum.

**The centre of this distribution is fixed by construction, not a finding.** Each percentile-rank
column has mean ~= 0.5, so the grand mean of the row means is ~= 0.5 necessarily. Only the spread and
shape carry information: at 260 endpoints the observed SD was 0.133 against the ~0.018 that mutually
independent endpoints would give, i.e. the endpoints share a great deal of signal. What that means
substantively — shared chemistry, correlated training data, genuine broad-spectrum compounds — is not
something these numbers settle.

Rows are NOT dropped for missing values: 15 compounds have one unscored endpoint
(``pfalciparum__eos4zfy__maip_score``), so their mean is taken over one fewer endpoint. The
per-compound CSV carries ``n_endpoints`` per row, making those visible rather than silently averaged
differently.

    python 07_score_matrices.py

Outputs
-------
    output/07_score_matrices/07_score_matrix_full.parquet   (cache: all referenced endpoints)
    output/07_score_matrices/07_mean_rank_per_compound.csv
    output/07_score_matrices/07_mean_rank_quantiles.csv
    output/07_score_matrices/{png,pdf}/07_mean_rank_distribution.{png,pdf}
    output/07_score_matrices/figure_cells.json
"""

import json
import os
import sys
import time

import numpy as np
import pandas as pd

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))

from default import ANNOTATION_PREDS_SUBDIR  # noqa: E402
from eval_correlations import build_full_library_matrix, build_named_score_matrix, scale_matrix  # noqa: E402
from plots_matrix_analyses import MeanRankDistributionPlot  # noqa: E402

pred_dir = os.path.join(root, "..", "data", "processed", ANNOTATION_PREDS_SUBDIR)
output_dir = os.path.join(root, "..", "output", "07_score_matrices")
config_dir = os.path.join(root, "..", "config")
os.makedirs(output_dir, exist_ok=True)

endpoint_selection_path = os.path.join(config_dir, "08_endpoint_selection.csv")
pathogens_of_interest_path = os.path.join(config_dir, "pathogens_of_interest.csv")
full_matrix_cache_path = os.path.join(output_dir, "07_score_matrix_full.parquet")

if not os.path.exists(full_matrix_cache_path):
    print(f"[score-matrices] building full-library cache from {pred_dir} "
          "(reads ~15 GB of raw prediction files — this can take several minutes)...")
    t0 = time.time()
    full = build_full_library_matrix(pred_dir, endpoint_selection_path)
    full.to_parquet(full_matrix_cache_path)
    print(f"[score-matrices] wrote {full_matrix_cache_path} {full.shape} in "
          f"{time.time() - t0:.1f}s")
    del full

#: Mean-percentile-rank outputs (former step 09). Skip-if-exists; the per-compound CSV is the
#: marker, since it is the expensive one to write.
mean_rank_csv = os.path.join(output_dir, "07_mean_rank_per_compound.csv")
mean_rank_quantiles_csv = os.path.join(output_dir, "07_mean_rank_quantiles.csv")
QUANTILES = [0.001, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 0.999]

mean_rank_pending = not os.path.exists(mean_rank_csv)

# --------------------------------------------------------------------------- #
# Mean percentile rank per compound (former step 09)                          #
# --------------------------------------------------------------------------- #
if not mean_rank_pending:
    print(f"[skip] {os.path.basename(mean_rank_csv)} already exists. Delete it to rebuild.")
else:
    base = build_named_score_matrix(
        pred_dir=pred_dir, endpoint_selection_path=endpoint_selection_path,
        pathogens_of_interest_path=pathogens_of_interest_path,
        full_matrix_cache_path=full_matrix_cache_path)
    t0 = time.time()
    ranked = scale_matrix(base, "rank_pct")
    print(f"[mean-rank] rank_pct computed in {time.time() - t0:.1f}s")

    # NaN-skipping mean: a compound with an unscored endpoint is averaged over the endpoints it
    # HAS, never dropped and never treated as zero. n_endpoints records which.
    mean_rank = ranked.mean(axis=1)
    n_endpoints = ranked.notna().sum(axis=1)
    partial = int((n_endpoints < ranked.shape[1]).sum())
    print(f"[mean-rank] {len(mean_rank):,} compounds; {partial} averaged over fewer than "
          f"{ranked.shape[1]} endpoints (missing predictions, kept not dropped)")

    per_compound = pd.DataFrame({"mean_percentile_rank": mean_rank, "n_endpoints": n_endpoints})
    t0 = time.time()
    per_compound.to_csv(mean_rank_csv)
    print(f"[mean-rank] wrote {os.path.basename(mean_rank_csv)} "
          f"({os.path.getsize(mean_rank_csv) / 1e9:.2f} GB) in {time.time() - t0:.1f}s")

    vals = mean_rank.to_numpy(dtype=float)
    q = pd.DataFrame({"quantile": QUANTILES,
                      "mean_percentile_rank": np.nanquantile(vals, QUANTILES)})
    q.to_csv(mean_rank_quantiles_csv, index=False)
    print(f"[mean-rank] mean={np.nanmean(vals):.4f}  SD={np.nanstd(vals):.4f}  "
          f"min={np.nanmin(vals):.4f}  max={np.nanmax(vals):.4f}")
    print(q.to_string(index=False, float_format=lambda v: f"{v:.4f}"))

    plot = MeanRankDistributionPlot(
        vals, label_note=f"{ranked.shape[1]} endpoints, {len(vals):,} compounds")
    plot.save(output_dir)
    with open(os.path.join(output_dir, "figure_cells.json"), "w") as f:
        json.dump({plot.name: list(plot.cells)}, f, indent=2)
    print(f"[mean-rank] figure: {plot.name}")
    del ranked

print(f"\nDone → {output_dir}")
