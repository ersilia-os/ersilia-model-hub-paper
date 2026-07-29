"""Step 06 figure — CoAdd own-strain AUROC per organism, both endpoints.

Reads the ``06_coadd_auroc.csv`` summary written by :func:`eval_coadd.run_all` and draws it with
the shared :class:`plots_euopenscreen.MetricByOrganismPlot` base, so the look matches every
other step. Colours come only from :mod:`plotting_colors`.
"""

import json
import os

import pandas as pd

from plotting_colors import CATEGORY_COLORS
from plots_euopenscreen import MetricByOrganismPlot

# CoAdd endpoints: mic_10 ~ dose-response (DR), inhib_50 ~ single-point (SP).
ENDPOINT_COLORS = {"mic_10": CATEGORY_COLORS["DR"], "inhib_50": CATEGORY_COLORS["SP"]}


class CoaddSharedAurocPlot(MetricByOrganismPlot):
    """CoAdd own-strain AUROC per shared organism, both endpoints (dedup)."""

    def __init__(self, coadd_df, ax=None, cells=(3, 3)):
        super().__init__(
            coadd_df, group_col="endpoint", group_colors=ENDPOINT_COLORS,
            group_order=["inhib_50", "mic_10"], metric="auroc", xlabel="AUROC",
            title="", name="coadd_shared_auroc",
            ref=0.5, xlim=(0, 1), prefer_best=True, cells=cells, ax=ax)


def save_coadd_figures(output_dir):
    """Build the CoAdd own-strain panel from ``06_coadd_auroc.csv`` in ``output_dir`` and record
    its footprint in ``figure_cells.json``."""
    path = os.path.join(output_dir, "06_coadd_auroc.csv")
    coadd = pd.read_csv(path) if os.path.exists(path) else pd.DataFrame()
    footprints = {}
    p = CoaddSharedAurocPlot(coadd, cells=(3, 3))
    if p.is_available:
        p.save(output_dir)
        footprints[p.name] = list(p.cells)
        print(f"  figure: {p.name}")
    else:
        print(f"  [skip figure] {p.name}: no data")
    with open(os.path.join(output_dir, "figure_cells.json"), "w") as f:
        json.dump(footprints, f, indent=2)
