"""Step 09 — distribution of each molecule's average percentile rank across all endpoints.

Takes the rank-percentile matrix (step 07's ``rank_pct`` scaling) and collapses it along its columns
to one number per compound: that compound's **mean percentile rank** over the 260 selected endpoints.
A percentile rank is a compound's standing within one endpoint over the whole library, so the mean is
"how highly does this molecule rank, on average, across everything we predict".

The request was phrased as a column-wise sum; the reported quantity is the **mean** (sum / 260), which
is what "average rank" means and what keeps the value on the interpretable [0, 1] percentile scale.
The two differ only by that constant factor — multiply by ``n_endpoints`` in the CSV for the sum.

**The centre of this distribution is fixed by construction, not a finding.** Each percentile-rank
column has mean ~= 0.5, so the grand mean of the row means is ~= 0.5 necessarily. Only the spread and
shape carry information: observed SD 0.133 against the ~0.018 that 260 mutually independent endpoints
would give, i.e. the endpoints share a great deal of signal. What that means substantively — shared
chemistry, correlated training data, genuine broad-spectrum compounds — is not something these numbers
settle.

Rows are NOT dropped for missing values: 15 compounds have one unscored endpoint
(``pfalciparum__eos4zfy__maip_score``), so their mean is taken over 259 endpoints instead of 260. The
per-compound CSV carries ``n_endpoints`` per row, making those visible rather than silently averaged
differently.

Rebuilds the matrix in memory from step 07's parquet cache (~40 s for the ranking) rather than reading
back the 6.8 GB rank-percentile CSV.

    python 09_mean_rank_distribution.py

Outputs
-------
    output/09_mean_rank_distribution/png/09_mean_rank_distribution.png
    output/09_mean_rank_distribution/pdf/09_mean_rank_distribution.pdf
    output/09_mean_rank_distribution/figure_cells.json
    output/09_mean_rank_distribution/09_mean_rank_per_compound.csv
    output/09_mean_rank_distribution/09_mean_rank_quantiles.csv
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
from eval_correlations import build_named_score_matrix, scale_matrix  # noqa: E402
from plots_matrix_analyses import MeanRankDistributionPlot  # noqa: E402

pred_dir = os.path.join(root, "..", "data", "processed", ANNOTATION_PREDS_SUBDIR)
matrix_dir = os.path.join(root, "..", "output", "07_score_matrices")
output_dir = os.path.join(root, "..", "output", "09_mean_rank_distribution")
config_dir = os.path.join(root, "..", "config")
os.makedirs(output_dir, exist_ok=True)

endpoint_selection_path = os.path.join(config_dir, "08_endpoint_selection.csv")
pathogens_of_interest_path = os.path.join(config_dir, "pathogens_of_interest.csv")
full_matrix_cache_path = os.path.join(matrix_dir, "07_score_matrix_full.parquet")

QUANTILES = [0.001, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 0.999]

t0 = time.time()
base = build_named_score_matrix(
    pred_dir=pred_dir, endpoint_selection_path=endpoint_selection_path,
    pathogens_of_interest_path=pathogens_of_interest_path,
    full_matrix_cache_path=full_matrix_cache_path)
print(f"[mean-rank] base named matrix {base.shape} ready in {time.time() - t0:.1f}s")

t0 = time.time()
ranked = scale_matrix(base, "rank_pct")
print(f"[mean-rank] rank_pct computed in {time.time() - t0:.1f}s")

mean_rank = ranked.mean(axis=1)          # skips NaN -> mean over the endpoints a compound has
n_endpoints = ranked.notna().sum(axis=1)
partial = int((n_endpoints < ranked.shape[1]).sum())
print(f"[mean-rank] {len(mean_rank):,} compounds; {partial} averaged over fewer than "
      f"{ranked.shape[1]} endpoints (missing predictions, kept not dropped)")

per_compound = pd.DataFrame({"mean_percentile_rank": mean_rank, "n_endpoints": n_endpoints})
out_csv = os.path.join(output_dir, "09_mean_rank_per_compound.csv")
t0 = time.time()
per_compound.to_csv(out_csv)
print(f"[mean-rank] wrote {out_csv} ({os.path.getsize(out_csv) / 1e9:.2f} GB) in "
      f"{time.time() - t0:.1f}s")

vals = mean_rank.to_numpy(dtype=float)
q = pd.DataFrame({"quantile": QUANTILES, "mean_percentile_rank": np.nanquantile(vals, QUANTILES)})
q.to_csv(os.path.join(output_dir, "09_mean_rank_quantiles.csv"), index=False)
print(f"[mean-rank] mean={np.nanmean(vals):.4f}  SD={np.nanstd(vals):.4f}  "
      f"min={np.nanmin(vals):.4f}  max={np.nanmax(vals):.4f}")
print(q.to_string(index=False, float_format=lambda v: f"{v:.4f}"))

plot = MeanRankDistributionPlot(
    vals, label_note=f"{ranked.shape[1]} endpoints, {len(vals):,} compounds")
plot.save(output_dir)
with open(os.path.join(output_dir, "figure_cells.json"), "w") as f:
    json.dump({plot.name: list(plot.cells)}, f, indent=2)
print(f"[mean-rank] figure: {plot.name}")

print(f"\nDone → {output_dir}")
