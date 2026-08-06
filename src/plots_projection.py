"""Step 09 figures — reference-library projection coloured by pathogen activity.

Reads ONLY the small per-method background CSVs and the single top-N-per-pathogen CSV written by
:func:`eval_projection.run_all` — never the raw prediction files or the eos1klk coordinate table.
One figure per projection method (:data:`default.PROJECTION_METHODS`): a :class:`GridPlot`
small-multiples grid, one panel per pathogen (``config/pathogens_of_interest.csv``), silver
background = full-library density, crimson overlay = that pathogen's ``PROJECTION_TOP_N``
highest-scoring compounds (a rank cutoff, never a score threshold) at that method's coordinates.
"""

import json
import os

import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
import stylia as st

from default import PROJECTION_METHODS, PROJECTION_TOP_N
from plotting_base import GridPlot
from plotting_colors import hue
from plotting_utils import abbrev, marker_legend, sequential_cmap

BACKGROUND_CMAP = sequential_cmap("silver")
POINT_COLOR = hue("crimson")
POINT_SIZE = 4

GRID_COLS = 4


def _pivot(df, value_col, n_bins):
    """Tidy ``bin_i, bin_j, value_col`` rows -> a ``(n_bins, n_bins)`` array, NaN where missing."""
    arr = np.full((n_bins, n_bins), np.nan)
    if len(df):
        arr[df["bin_i"].to_numpy(), df["bin_j"].to_numpy()] = df[value_col].to_numpy()
    return arr


def _grid_extent(df):
    """``(xmin, xmax, ymin, ymax)`` cell EDGES from a tidy grid's cell-CENTER columns.

    ``imshow`` needs this passed explicitly as ``extent=`` — left to its default it draws in
    pixel-index coordinates (``0..n_bins``), not the real projection coordinates the top-N scatter
    is plotted in, and the two would silently disagree (the scatter would collapse into whatever
    tiny corner ``[-1, 1]``-ish data occupies inside a ``[0, 60]``-ish image).
    """
    xs, ys = np.sort(df["x_center"].unique()), np.sort(df["y_center"].unique())
    dx = (xs[-1] - xs[0]) / (len(xs) - 1) if len(xs) > 1 else 1.0
    dy = (ys[-1] - ys[0]) / (len(ys) - 1) if len(ys) > 1 else 1.0
    return (xs[0] - dx / 2, xs[-1] + dx / 2, ys[0] - dy / 2, ys[-1] + dy / 2)


class PathogenProjectionGridPlot(GridPlot):
    """One projection method's small-multiples grid: one panel per pathogen.

    Every panel shares the same silver background (the full library's density in this method's
    projection), so the 15 panels — and the same panels across the other three methods' figures —
    are directly comparable at a glance. The overlay is that pathogen's top-N points only, never a
    continuous score, so no colour scale/colourbar is needed — just the one marker key.
    """

    def __init__(self, method, background, top_n_table, pathogens, top_n=PROJECTION_TOP_N,
                cols=GRID_COLS):
        n_bins = int(max(background["bin_i"].max(), background["bin_j"].max()) + 1) \
            if len(background) else 0
        self.bg_arr = _pivot(background, "n_compounds", n_bins) if n_bins else None
        self._bg_extent = _grid_extent(background) if n_bins else None
        bg_vmax = float(np.nanpercentile(self.bg_arr, 99)) if self.bg_arr is not None else 1.0
        self._bg_norm = mcolors.Normalize(vmin=0.0, vmax=bg_vmax or 1.0)

        items = []
        for _, row in pathogens.iterrows():
            code, name = row["code"], row["pathogen"]
            g = top_n_table[top_n_table["pathogen_code"] == code] \
                if len(top_n_table) else top_n_table
            items.append({
                "code": code, "pathogen": name,
                "x": g[f"{method}_x"].to_numpy() if len(g) else np.array([]),
                "y": g[f"{method}_y"].to_numpy() if len(g) else np.array([]),
            })

        self.build_grid(items, cols=cols, name=f"09_{method}_top{top_n}_pathogens",
                        panel_fn=self._panel, edge_xlabel=f"{method.upper()} 1",
                        edge_ylabel=f"{method.upper()} 2")
        if self.is_available:
            marker_legend(self.ax, [{"label": f"top {top_n} predicted",
                                     "color": POINT_COLOR, "markersize": 3}])

    def _panel(self, ax, item, color, xlabel, ylabel):
        if self.bg_arr is not None:
            ax.imshow(self.bg_arr.T, origin="lower", cmap=BACKGROUND_CMAP, norm=self._bg_norm,
                      extent=self._bg_extent, aspect="auto")
        ax.scatter(item["x"], item["y"], s=POINT_SIZE, color=POINT_COLOR, edgecolors="none",
                  zorder=3)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.text(0.03, 0.94, abbrev(item["pathogen"]), transform=ax.transAxes,
               ha="left", va="top", fontsize=st.FONTSIZE_SMALL)
        st.label(ax, xlabel=xlabel, ylabel=ylabel)


# --------------------------------------------------------------------------- #
# Entry point                                                                  #
# --------------------------------------------------------------------------- #
def _read(output_dir, fname, **kw):
    path = os.path.join(output_dir, fname)
    return pd.read_csv(path, **kw) if os.path.exists(path) else pd.DataFrame()


def save_projection_figures(output_dir, pathogens_csv, top_n=PROJECTION_TOP_N):
    """Build every step-09 figure (one per projection method) from the summary CSVs in
    ``output_dir`` and record their footprints in ``figure_cells.json``."""
    pathogens = pd.read_csv(pathogens_csv)
    top_n_table = _read(output_dir, f"09_top{top_n}_per_pathogen.csv")
    footprints = {}
    for method in PROJECTION_METHODS:
        background = _read(output_dir, f"09_{method}_background.csv")
        plot = PathogenProjectionGridPlot(method, background, top_n_table, pathogens, top_n=top_n)
        if plot.is_available:
            plot.save(output_dir)
            footprints[plot.name] = list(plot.cells)
            print(f"  figure: {plot.name}")
        else:
            print(f"  [skip figure] 09_{method}_top{top_n}_pathogens: insufficient data")
    with open(os.path.join(output_dir, "figure_cells.json"), "w") as f:
        json.dump(footprints, f, indent=2)
