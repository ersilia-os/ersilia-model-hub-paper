"""Step 14 figures — pathogen hits vs antibiotic-resemblance hits on the library UMAP.

Reads ONLY the small summary CSVs written by :mod:`eval_pathogen_abx_overlap` and step 10's
existing ``09_umap_background.csv`` — never a prediction file or a full matrix.

One figure PER PATHOGEN, each a small-multiples grid with one panel per abx endpoint. Every panel
carries the same silver full-library density behind three point layers:

    crimson   this pathogen's top-N by consensus_score, NOT in the abx set
    cobalt    this abx endpoint's highlighted compounds, NOT in the pathogen set
    lime      the intersection — in both

The intersection draws last and larger, because it is the subject of the figure and would
otherwise be buried under whichever single-membership layer happened to be drawn over it.
"""

import json
import os

import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
import stylia as st

from default import PROJECTION_TOP_N
from plotting_base import GridPlot
from plotting_colors import hue
from plotting_utils import abbrev, marker_legend, merge_figure_cells, sequential_cmap
from plots_abx_projection import _wrap_name
from plots_projection import _grid_extent, _pivot

BACKGROUND_CMAP = sequential_cmap("silver")

#: group -> (colour, marker size, z-order). The intersection is drawn last, largest and in a hue
#: that is neither of the two single-membership colours, so it can never read as a shade of one of
#: them. crimson keeps its step-10 meaning ("this pathogen's top-N").
GROUP_STYLE = {
    "pathogen": (hue("crimson"), 3, 3),
    "abx": (hue("cobalt"), 3, 4),
    "both": (hue("lime"), 8, 5),
}
DRAW_ORDER = ("pathogen", "abx", "both")

#: Panel annotations must clear every scatter layer (z 3-5) or the points bury them.
LABEL_Z = 6

#: Nine endpoints over three columns → a 3x3-cell footprint (90 mm square, half the print width),
#: keeping each panel at the same ~30 mm as steps 10 and 13.
GRID_COLS = 3


class PathogenAbxOverlapGridPlot(GridPlot):
    """One pathogen's grid: a panel per abx endpoint, three point groups each."""

    def __init__(self, pathogen_code, pathogen_name, background, points, counts, endpoints,
                 cols=GRID_COLS):
        n_bins = int(max(background["bin_i"].max(), background["bin_j"].max()) + 1) \
            if len(background) else 0
        self.bg_arr = _pivot(background, "n_compounds", n_bins) if n_bins else None
        self._bg_extent = _grid_extent(background) if n_bins else None
        bg_vmax = float(np.nanpercentile(self.bg_arr, 99)) if self.bg_arr is not None else 1.0
        self._bg_norm = mcolors.Normalize(vmin=0.0, vmax=bg_vmax or 1.0)

        items = []
        for e in endpoints.itertuples():
            g = points[points["endpoint"] == e.endpoint]
            c = counts[counts["endpoint"] == e.endpoint]
            n_both = int(c["n_both"].iloc[0]) if len(c) else 0
            capped = bool(c["abx_capped"].iloc[0]) if len(c) else False
            items.append({
                "column_name": e.column_name, "model_id": e.model_id,
                "n_both": n_both, "capped": capped,
                "groups": {k: g[g["group"] == k] for k in DRAW_ORDER},
            })

        self.build_grid(items, cols=cols, name=f"11_umap_abx_overlap_{pathogen_code}",
                        panel_fn=self._panel, edge_xlabel="UMAP 1", edge_ylabel="UMAP 2")
        if self.is_available:
            # Abbreviated genus, and the "top N" left to the caption: the full binomial plus the
            # count is wider than a 3 cm panel and the legend box spills outside the grid.
            marker_legend(self.ax, [
                {"label": abbrev(pathogen_name),
                 "color": GROUP_STYLE["pathogen"][0], "markersize": 3},
                {"label": "antibiotic-like", "color": GROUP_STYLE["abx"][0], "markersize": 3},
                {"label": "both", "color": GROUP_STYLE["both"][0], "markersize": 4},
            ])

    def _panel(self, ax, item, color, xlabel, ylabel):
        if self.bg_arr is not None:
            ax.imshow(self.bg_arr.T, origin="lower", cmap=BACKGROUND_CMAP, norm=self._bg_norm,
                      extent=self._bg_extent, aspect="auto")
        for group in DRAW_ORDER:
            g = item["groups"][group]
            if not len(g):
                continue
            c, size, z = GROUP_STYLE[group]
            ax.scatter(g["umap_x"], g["umap_y"], s=size, color=c, edgecolors="none", zorder=z)
        ax.set_xticks([])
        ax.set_yticks([])
        name = _wrap_name(item["column_name"])
        ax.text(0.03, 0.95, name, transform=ax.transAxes,
                ha="left", va="top", fontsize=st.FONTSIZE_SMALL, zorder=LABEL_Z)
        # The intersection size is the number the panel exists to convey, and at these point sizes
        # a handful of lime dots is easy to miss. A trailing * marks an endpoint whose abx set hit
        # step 12's cap, i.e. one where a small overlap may just reflect an arbitrary subset.
        # Sits below however many lines the wrapped name took, so the two never collide.
        ax.text(0.03, 0.95 - 0.13 * (name.count("\n") + 1),
                f"{item['model_id']}  n={item['n_both']}{'*' if item['capped'] else ''}",
                transform=ax.transAxes, ha="left", va="top",
                fontsize=st.FONTSIZE_SMALL, alpha=0.6, zorder=LABEL_Z)
        st.label(ax, xlabel=xlabel, ylabel=ylabel)


# --------------------------------------------------------------------------- #
# Entry point                                                                  #
# --------------------------------------------------------------------------- #
def save_overlap_figures(output_dir, background_csv, pathogens, endpoints, counts):
    """One figure per pathogen, plus a shared ``figure_cells.json`` of their footprints."""
    background = pd.read_csv(background_csv) if os.path.exists(background_csv) else pd.DataFrame()
    points = pd.read_csv(os.path.join(output_dir, "11_overlap_points.csv"))

    footprints = {}
    for p in pathogens.itertuples():
        plot = PathogenAbxOverlapGridPlot(
            p.code, p.pathogen, background,
            points[points["pathogen_code"] == p.code],
            counts[counts["pathogen_code"] == p.code], endpoints)
        if plot.is_available:
            plot.save(output_dir)
            footprints[plot.name] = list(plot.cells)
            print(f"  figure: {plot.name}")
        else:
            print(f"  [skip figure] {p.code}: insufficient data")
    merge_figure_cells(output_dir, footprints)
