"""Step 13 figures — reference-library projection coloured by predicted toxicity.

Reads ONLY the two small summary CSVs written by :func:`eval_tox_projection.run_all` — never the
step-11 property table or the eos1klk coordinate table. One :class:`GridPlot` small-multiples
figure: one panel per selected toxicity endpoint, silver background = full-library density,
crimson overlay = that endpoint's ``PROJECTION_TOP_N`` most toxic compounds (a rank cutoff, never
a score threshold), all on the same :data:`default.TOX_PROJECTION_METHOD` layout used by step 10's
pathogen figures, so the two are directly comparable.
"""

import json
import os

import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
import stylia as st

from default import PROJECTION_TOP_N, TOX_PROJECTION_METHOD
from plotting_base import GridPlot
from plotting_colors import hue
from plotting_utils import marker_legend, merge_figure_cells, sequential_cmap
from plots_projection import _grid_extent, _pivot

BACKGROUND_CMAP = sequential_cmap("silver")
POINT_COLOR = hue("crimson")
POINT_SIZE = 4

#: 24 endpoints over 6 columns → a 4x6-cell footprint, i.e. the full 180 mm print width at the
#: same ~30x30 mm panel size as step 10's pathogen grid.
GRID_COLS = 6


class ToxicityProjectionGridPlot(GridPlot):
    """One panel per toxicity endpoint, all sharing one projection layout and background.

    Every panel carries the same silver full-library density, so panels are comparable to each
    other and to step 10's pathogen panels at a glance. The overlay is the endpoint's top-N points
    only, never a continuous score, so no colour scale is needed — just the one marker key.
    """

    def __init__(self, method, background, top_n_table, endpoints, top_n=PROJECTION_TOP_N,
                 cols=GRID_COLS):
        n_bins = int(max(background["bin_i"].max(), background["bin_j"].max()) + 1) \
            if len(background) else 0
        self.bg_arr = _pivot(background, "n_compounds", n_bins) if n_bins else None
        self._bg_extent = _grid_extent(background) if n_bins else None
        bg_vmax = float(np.nanpercentile(self.bg_arr, 99)) if self.bg_arr is not None else 1.0
        self._bg_norm = mcolors.Normalize(vmin=0.0, vmax=bg_vmax or 1.0)

        items = []
        for r in endpoints.itertuples():
            g = top_n_table[top_n_table["endpoint"] == r.endpoint] \
                if len(top_n_table) else top_n_table
            items.append({
                "column_name": r.column_name, "model_id": r.model_id,
                "x": g[f"{method}_x"].to_numpy() if len(g) else np.array([]),
                "y": g[f"{method}_y"].to_numpy() if len(g) else np.array([]),
            })

        self.build_grid(items, cols=cols, name=f"12_{method}_top{top_n}_toxicity",
                        panel_fn=self._panel, edge_xlabel=f"{method.upper()} 1",
                        edge_ylabel=f"{method.upper()} 2")
        if self.is_available:
            marker_legend(self.ax, [{"label": f"top {top_n} most toxic",
                                     "color": POINT_COLOR, "markersize": 3}])

    def _panel(self, ax, item, color, xlabel, ylabel):
        if self.bg_arr is not None:
            ax.imshow(self.bg_arr.T, origin="lower", cmap=BACKGROUND_CMAP, norm=self._bg_norm,
                      extent=self._bg_extent, aspect="auto")
        ax.scatter(item["x"], item["y"], s=POINT_SIZE, color=POINT_COLOR, edgecolors="none",
                   zorder=3)
        ax.set_xticks([])
        ax.set_yticks([])
        # Endpoint name over its model ID: several endpoints measure overlapping biology from
        # different models (two HepG2 readouts, three cytotoxicity readouts), so the panel is only
        # unambiguous with the provenance shown.
        ax.text(0.03, 0.95, item["column_name"], transform=ax.transAxes,
                ha="left", va="top", fontsize=st.FONTSIZE_SMALL)
        ax.text(0.03, 0.82, item["model_id"], transform=ax.transAxes,
                ha="left", va="top", fontsize=st.FONTSIZE_SMALL, alpha=0.6)
        st.label(ax, xlabel=xlabel, ylabel=ylabel)


# --------------------------------------------------------------------------- #
# Entry point                                                                  #
# --------------------------------------------------------------------------- #
def _read(output_dir, fname, **kw):
    path = os.path.join(output_dir, fname)
    return pd.read_csv(path, **kw) if os.path.exists(path) else pd.DataFrame()


def save_tox_projection_figures(output_dir, endpoints, method=TOX_PROJECTION_METHOD,
                                top_n=PROJECTION_TOP_N):
    """Build the step-12 toxicity figure from the summary CSVs in ``output_dir`` and record its footprint
    in ``figure_cells.json``."""
    background = _read(output_dir, f"12_{method}_background.csv")
    top_n_table = _read(output_dir, f"12_top{top_n}_per_endpoint.csv")
    footprints = {}
    plot = ToxicityProjectionGridPlot(method, background, top_n_table, endpoints, top_n=top_n)
    if plot.is_available:
        plot.save(output_dir)
        footprints[plot.name] = list(plot.cells)
        print(f"  figure: {plot.name}")
    else:
        print(f"  [skip figure] 12_{method}_top{top_n}_toxicity: insufficient data")
    merge_figure_cells(output_dir, footprints)
