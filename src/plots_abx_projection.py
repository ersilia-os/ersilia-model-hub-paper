"""Step 12 figures — antibiotic-resemblance endpoints on the reference-library UMAP.

Reads ONLY the small summary CSVs written by :mod:`eval_abx_matrix` (the per-endpoint highlight
table) and step 11's existing ``11_umap_background.csv`` — never the 1.35M x 55 matrix or the raw
prediction files, per the repo's "feed figures from summary CSVs" rule.

One small-multiples grid, one panel per endpoint: silver background = full-library density in UMAP
space, crimson overlay = that endpoint's highlighted compounds.

**The highlight rule is NOT step 11's rank cutoff**, because these endpoints are not continuous.
54 of the 55 selected columns are binary flags or small integer counts, so a plain "top 1000" would
pad most panels with arbitrarily chosen zero-valued compounds. Instead each panel shows every
compound with a NON-ZERO value, capped at :data:`default.PROJECTION_TOP_N` (highest value first).
Nothing is ever padded with zeros, so a panel with 3 hits draws 3 points. Each panel is annotated
``n_shown/n_nonzero`` so a capped panel is never mistaken for an exhaustive one.
"""

import os

import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
import stylia as st

from default import PROJECTION_TOP_N
from plotting_base import GridPlot
from plotting_colors import hue
from plotting_utils import marker_legend, merge_figure_cells
from plots_projection import _grid_extent, _pivot

BACKGROUND_CMAP = None  # set lazily; sequential_cmap import kept local to avoid a cycle
POINT_COLOR = hue("crimson")
POINT_SIZE = 3

#: Six columns is the one width that lands the grid on exactly stylia's 180 mm print width
#: (``CELLS_PER_WIDTH`` = 6), giving true 3 cm panels. Seven would render 210 mm wide — over the
#: Nature two-column width the whole cell grid exists to respect.
GRID_COLS = 6


#: Endpoint names longer than this are wrapped onto a second line at an underscore. Names like
#: ``aminopyrimidines_trimethoprim_der`` are wider than a 3 cm panel and would otherwise run across
#: the neighbouring one — nothing is abbreviated away, only wrapped.
NAME_WRAP = 20


def _wrap_name(name, width=NAME_WRAP):
    """Break a long ``snake_case`` endpoint name onto multiple lines at underscore boundaries."""
    if len(name) <= width:
        return name
    lines, current = [], ""
    for part in name.split("_"):
        candidate = f"{current}_{part}" if current else part
        if len(candidate) > width and current:
            lines.append(current)
            current = part
        else:
            current = candidate
    if current:
        lines.append(current)
    return "\n".join(lines)


class AbxProjectionGridPlot(GridPlot):
    """Small-multiples UMAP grid, one panel per antibiotic-resemblance endpoint.

    Every panel shares the same silver background (the full library's UMAP density), so panels are
    directly comparable at a glance. The overlay is a point set, never a continuous score, so no
    colour scale is needed — just the one marker key.

    Endpoints with ZERO non-zero compounds over the whole library are expected to have been dropped
    upstream by :func:`eval_abx_matrix.endpoint_highlights` (user-directed): they have nothing to
    draw, and 1000 tie-broken zeros would render as a fabricated cloud. They remain listed in the
    stats CSV and are printed by the script.
    """

    def __init__(self, background, highlights, stats, cols=GRID_COLS, top_n=PROJECTION_TOP_N):
        from plotting_utils import sequential_cmap

        bg_cmap = sequential_cmap("silver")
        n_bins = int(max(background["bin_i"].max(), background["bin_j"].max()) + 1) \
            if len(background) else 0
        self.bg_arr = _pivot(background, "n_compounds", n_bins) if n_bins else None
        self._bg_extent = _grid_extent(background) if n_bins else None
        self._bg_cmap = bg_cmap
        bg_vmax = float(np.nanpercentile(self.bg_arr, 99)) if self.bg_arr is not None else 1.0
        self._bg_norm = mcolors.Normalize(vmin=0.0, vmax=bg_vmax or 1.0)

        items = []
        for _, row in stats.iterrows():
            g = highlights[highlights["endpoint"] == row["endpoint"]] if len(highlights) \
                else highlights
            items.append({
                "endpoint": row["endpoint"],
                "column_name": row["column_name"],
                "n_shown": int(row["n_shown"]),
                "n_nonzero": int(row["n_nonzero"]),
                "x": g["umap_x"].to_numpy() if len(g) else np.array([]),
                "y": g["umap_y"].to_numpy() if len(g) else np.array([]),
            })

        self.build_grid(items, cols=cols, name=f"12_umap_abx_endpoints_max{top_n}",
                        panel_fn=self._panel, edge_xlabel="UMAP 1", edge_ylabel="UMAP 2")
        if self.is_available:
            # Into a trailing EMPTY cell, not over the first panel: every panel's lower-left
            # carries its n_shown/n_nonzero annotation, and a legend there hides the one number
            # that stops a capped panel being read as exhaustive. Falls back to the first panel
            # only when the grid happens to be exactly full.
            spare = self.fig.axes[len(items):]
            marker_legend(spare[0] if spare else self.ax,
                          [{"label": f"non-zero compounds (max {top_n})",
                            "color": POINT_COLOR, "markersize": 3}],
                          loc="center" if spare else "lower right")

    def _panel(self, ax, item, color, xlabel, ylabel):
        if self.bg_arr is not None:
            ax.imshow(self.bg_arr.T, origin="lower", cmap=self._bg_cmap, norm=self._bg_norm,
                      extent=self._bg_extent, aspect="auto")
        ax.scatter(item["x"], item["y"], s=POINT_SIZE, color=POINT_COLOR, edgecolors="none",
                   zorder=3)
        ax.set_xticks([])
        ax.set_yticks([])
        # The count is part of the panel's meaning, not decoration: "1000" alone cannot be told
        # apart from "all of them" without the total beside it.
        ax.text(0.03, 0.97, _wrap_name(item["column_name"]), transform=ax.transAxes,
                ha="left", va="top", fontsize=st.FONTSIZE_SMALL, linespacing=1.1)
        ax.text(0.03, 0.03, f"{item['n_shown']:,}/{item['n_nonzero']:,}", transform=ax.transAxes,
                ha="left", va="bottom", fontsize=st.FONTSIZE_SMALL)
        st.label(ax, xlabel=xlabel, ylabel=ylabel)


# --------------------------------------------------------------------------- #
# Entry point                                                                  #
# --------------------------------------------------------------------------- #
def save_abx_projection_figure(output_dir, background_path, highlights, stats,
                               top_n=PROJECTION_TOP_N):
    """Build the step-12 UMAP grid from the summary CSVs and record its footprint.

    ``background_path`` points at step 11's ``11_umap_background.csv`` — the same full-library
    density grid, reused rather than recomputed so the two steps' panels are directly comparable.
    """
    if not os.path.exists(background_path):
        raise FileNotFoundError(
            f"Missing {background_path}. Run `python 11_reference_library_projection.py` first — "
            "step 12 reuses its UMAP background grid rather than recomputing it.")
    background = pd.read_csv(background_path)

    plot = AbxProjectionGridPlot(background, highlights, stats, top_n=top_n)
    footprints = {}
    if plot.is_available:
        plot.save(output_dir)
        footprints[plot.name] = list(plot.cells)
        print(f"  figure: {plot.name}  ({plot.cells[0]} x {plot.cells[1]} cells)")
    else:
        print("  [skip figure] 11_umap_abx_endpoints: no endpoints with data")
    merge_figure_cells(output_dir, footprints)
    return plot
