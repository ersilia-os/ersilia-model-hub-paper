"""Step 11 figures — reference-library projection coloured by pathogen activity.

Reads ONLY the small per-method background CSVs and the top-N summary CSVs written by
:func:`eval_projection.run_all` / :func:`eval_projection.run_auroc_matched_umap` — never the raw
prediction files or the eos1klk coordinate table. One figure per projection method
(:data:`default.PROJECTION_METHODS`): a :class:`GridPlot` small-multiples grid, one panel per
pathogen (``config/pathogens_of_interest.csv``), silver background = full-library density, crimson
overlay = that pathogen's ``PROJECTION_TOP_N`` highest-scoring compounds (a rank cutoff, never a
score threshold) at that method's coordinates.

The UMAP panel is the exception: :class:`AurocMatchedUmapGridPlot` replaces it with a 3x5 grid
ordered and ranked to match step 10's AUROC matrix exactly (the AUROC matrix's own top-N per
organism, its own phylogeny-within-class row order) rather than this script's own alphabetical
order and per-pathogen ``consensus_score`` ranking — see that class's docstring. PCA/t-SNE/TMAP keep
:class:`PathogenProjectionGridPlot` unchanged.

Plus a second family on the same background — see :class:`CoaddPathogenProjectionGridPlot` — driven
by ``COADD_MODEL_ID`` (eos3dys) instead of the per-pathogen ChEMBL models, so the two can be read
side by side.
"""

import os

import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
import stylia as st
from scipy.stats import gaussian_kde

from default import COADD_MODEL_ID, COADD_PROJECTION_METHOD, PROJECTION_METHODS, PROJECTION_TOP_N
from plotting_base import GridPlot, MultiPanelPlot
from plotting_colors import hue, shades
from plotting_utils import abbrev, marker_legend, merge_figure_cells, sequential_cmap

BACKGROUND_CMAP = sequential_cmap("silver")
POINT_COLOR = hue("crimson")
POINT_SIZE = 4

GRID_COLS = 4

#: eos3dys's 9 organisms over three columns → a 3x3-cell footprint (90 mm square, half the print
#: width), keeping each panel at the same ~30 mm as the pathogen grid above and steps 12-14.
COADD_GRID_COLS = 3

#: Total footprint for the AUROC-matched UMAP grid: 0.5 of stylia's two-column print width
#: (180 mm) is 90 mm; a strictly-square 3x5 grid at 90 mm needs 18 mm panels (90/5), so height =
#: 3 x 18 = 54 mm. On the repo's 3 cm/6-cell reference grid that is cells=(1.8, 3.0) — a fractional
#: footprint, the sanctioned exception documented in docs/figure_conventions.md (also used by
#: task_subtask_waffle cells=(1.5,1.5) and pathogen_circles cells=(2.5,2.5)).
AUROC_MATCHED_UMAP_ROWS = 3
AUROC_MATCHED_UMAP_COLS = 5
AUROC_MATCHED_UMAP_CELLS = (1.8, 3.0)

#: Fixed axis range (user-directed). The real UMAP extent is slightly wider on both axes (roughly
#: -1.17 to 1.09, see 11_umap_background.csv), so a thin sliver of background and possibly a few
#: top-1000 points near the outer edge fall outside this fixed frame and are not drawn.
AUROC_MATCHED_UMAP_LIM = (-1, 1)

#: One shared density colormap across all 15 panels (user-directed, revised 2026-09-01 — the
#: original 4-hue OVERLAP_MATRIX_SPECTRUM read as busy/uneven at this scale): the house
#: "coolwarm"-equivalent, ``stylia.DivergingColormap("crimson_cobalt")`` — the same diverging pair
#: the AUROC matrix itself uses — REVERSED so low density sits at the cool cobalt end and high
#: density at the warm crimson end (crimson_cobalt's un-reversed direction runs crimson at 0). A
#: sequential reading of a nominally-diverging colormap is intentional: it borrows the pleasant,
#: legible two-hue gradient without implying a meaningful zero-centre the way the AUROC matrix's own
#: use of it does.
AUROC_MATCHED_UMAP_DENSITY_CMAP = st.DivergingColormap("crimson_cobalt").cmap.reversed()

#: Scatter marker size for this grid only (smaller than the shared POINT_SIZE used by
#: PathogenProjectionGridPlot — user-directed 2026-09-01): 18 mm panels read better with a finer
#: mark than the 30 mm panels those grids use.
AUROC_MATCHED_UMAP_POINT_SIZE = 2.5

#: Marker size per endpoint layer, ``(darkest, palest)``. An organism's endpoints are strongly
#: correlated readouts of the same biology — ``inhib_50`` and ``mic_25`` of one strain especially —
#: so their layers overlap heavily, and at one size whichever drew last would hide the rest.
#:
#: The pairing runs PALEST-largest-underneath, DARKEST-smallest-on-top, so an overlap reads as a
#: dark core inside a pale halo. The other way round (pale on top of dark) was tried first and is
#: worse twice over: pale-on-dark is the harder contrast to read, and it puts the palest layer —
#: already the weakest — at the smallest size as well.
#:
#: Both values sit BELOW the single-layer ``POINT_SIZE`` of 4: a panel here carries up to 4000
#: points rather than 1000, and at 30 mm the density, not the dot size, is what saturates it.
COADD_SIZE_RANGE = (1.2, 3.0)

#: Panel annotations must clear every scatter layer or the points bury them.
LABEL_Z = 6

#: Palest tint allowed for an endpoint shade — capped above the generic
#: :data:`plotting_colors.SHADE_LIGHTEN_FLOOR` (0.4) because here the shade also colours the panel's
#: TEXT key, and thin letter strokes need more contrast against white than a filled dot does.
COADD_SHADE_FLOOR = 0.5

#: Vertical step between the stacked text lines in a panel, in axes fraction. Four endpoints is the
#: most any organism has, so the key stays within the panel's top half.
LABEL_STEP = 0.09

#: The key sits ON the data, so each line gets the same semi-transparent white backing the house
#: legends use (``LEGEND_KW``) — without it the labels are unreadable wherever a point cloud reaches
#: the top-left corner, which in these panels is most of them.
LABEL_BBOX = dict(facecolor="white", alpha=0.7, edgecolor="none", pad=0.8)


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


def _background_layer(background):
    """``(density array, imshow extent, norm)`` for the shared silver full-library layer.

    The norm clips at the 99th percentile so that a handful of very dense cells cannot wash out the
    rest of the map. A missing/empty background gives ``(None, None, norm)``, which every panel
    treats as "draw the overlay only" rather than failing.
    """
    n_bins = int(max(background["bin_i"].max(), background["bin_j"].max()) + 1) \
        if len(background) else 0
    arr = _pivot(background, "n_compounds", n_bins) if n_bins else None
    extent = _grid_extent(background) if n_bins else None
    vmax = float(np.nanpercentile(arr, 99)) if arr is not None else 1.0
    return arr, extent, mcolors.Normalize(vmin=0.0, vmax=vmax or 1.0)


def _draw_background(ax, arr, extent, norm):
    """Draw the silver density layer. ``.T`` because ``_pivot`` is indexed ``[bin_i, bin_j]``."""
    if arr is not None:
        ax.imshow(arr.T, origin="lower", cmap=BACKGROUND_CMAP, norm=norm, extent=extent,
                  aspect="auto")


class PathogenProjectionGridPlot(GridPlot):
    """One projection method's small-multiples grid: one panel per pathogen.

    Every panel shares the same silver background (the full library's density in this method's
    projection), so the 15 panels — and the same panels across the other three methods' figures —
    are directly comparable at a glance. The overlay is that pathogen's top-N points only, never a
    continuous score, so no colour scale/colourbar is needed — just the one marker key.
    """

    def __init__(self, method, background, top_n_table, pathogens, top_n=PROJECTION_TOP_N,
                cols=GRID_COLS):
        self.bg_arr, self._bg_extent, self._bg_norm = _background_layer(background)

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

        self.build_grid(items, cols=cols, name=f"11_{method}_top{top_n}_pathogens",
                        panel_fn=self._panel, edge_xlabel=f"{method.upper()} 1",
                        edge_ylabel=f"{method.upper()} 2")
        if self.is_available:
            marker_legend(self.ax, [{"label": f"top {top_n} predicted",
                                     "color": POINT_COLOR, "markersize": 3}])

    def _panel(self, ax, item, color, xlabel, ylabel):
        _draw_background(ax, self.bg_arr, self._bg_extent, self._bg_norm)
        ax.scatter(item["x"], item["y"], s=POINT_SIZE, color=POINT_COLOR, edgecolors="none",
                  zorder=3)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.text(0.03, 0.94, abbrev(item["pathogen"]), transform=ax.transAxes,
               ha="left", va="top", fontsize=st.FONTSIZE_SMALL)
        st.label(ax, xlabel=xlabel, ylabel=ylabel)


class AurocMatchedUmapGridPlot(MultiPanelPlot):
    """UMAP top-1000-per-organism, ordered and coloured to match step 10's AUROC matrix exactly.

    Reads :func:`eval_projection.auroc_matched_top_n_per_organism`'s table: the AUROC matrix's OWN
    rank-pct-averaged per-organism top-N compounds (not this script's per-pathogen
    ``consensus_score`` top-N — see :class:`PathogenProjectionGridPlot`), in the AUROC matrix's own
    phylogeny-within-class row order. A 3x5 grid at exactly half stylia's two-column print width
    (``cells=AUROC_MATCHED_UMAP_CELLS``), with a density-shaded (not flat-colour) scatter, so it
    reads as the AUROC matrix's direct visual companion rather than as a fourth
    :class:`PathogenProjectionGridPlot` method.

    Axis chrome (ticks, frame, x/y labels) is dropped entirely — there are 15 identical panels and
    nothing method- or unit-specific to label, unlike :class:`PathogenProjectionGridPlot`'s edge
    ``"{METHOD} 1"/"{METHOD} 2"`` labels.
    """

    def __init__(self, background, table, method="umap", cells=AUROC_MATCHED_UMAP_CELLS,
                name="11_umap_top1000_pathogens"):
        self.name = name
        self.is_available = len(table) > 0
        if not self.is_available:
            self.cells = cells
            return

        self.bg_arr, self._bg_extent, self._bg_norm = _background_layer(background)
        axs = self._new_figure(AUROC_MATCHED_UMAP_ROWS, AUROC_MATCHED_UMAP_COLS, cells, name)

        organisms = (table[["organism", "phylo_position"]]
                    .drop_duplicates().sort_values("phylo_position"))
        for ax, row in zip(axs.axs_flat, organisms.itertuples()):
            g = table[table["organism"] == row.organism]
            self._panel(ax, row.organism,
                       g[f"{method}_x"].to_numpy(), g[f"{method}_y"].to_numpy())
        self.ax = axs[0]

    def _panel(self, ax, organism, x, y):
        _draw_background(ax, self.bg_arr, self._bg_extent, self._bg_norm)

        valid = np.isfinite(x) & np.isfinite(y)
        x, y = x[valid], y[valid]
        # gaussian_kde needs at least a few points to fit a covariance; a near-empty organism (data
        # gap, not expected in practice) falls back to a flat colour rather than raising.
        if len(x) >= 3:
            density = gaussian_kde(np.vstack([x, y]))(np.vstack([x, y]))
            # Draw sparsest points first, densest last: with a single scatter() call, array order IS
            # z-order, so without this a crowded, high-density cluster can be partially painted over
            # by whatever sparse points happen to sit later in the (arbitrary, top-N-merge) order.
            order = np.argsort(density)
            x, y, density = x[order], y[order], density[order]
            norm = mcolors.Normalize(vmin=float(density.min()), vmax=float(density.max()))
            ax.scatter(x, y, c=density, cmap=AUROC_MATCHED_UMAP_DENSITY_CMAP, norm=norm,
                      s=AUROC_MATCHED_UMAP_POINT_SIZE, edgecolors="none", zorder=3)
        else:
            ax.scatter(x, y, color=hue("cobalt"), s=AUROC_MATCHED_UMAP_POINT_SIZE, edgecolors="none",
                      zorder=3)

        ax.set_xlim(*AUROC_MATCHED_UMAP_LIM)
        ax.set_ylim(*AUROC_MATCHED_UMAP_LIM)
        ax.set_aspect("equal")
        ax.axis("off")

        # A real axis title (loc="left"), not in-panel text: it sits in the space matplotlib reserves
        # above the axes, so it can never land on top of a scatter point — the same mechanism
        # plotting_utils.roc_panel uses for its own per-cell identifier. axis("off") still leaves it
        # visible (it hides ticks/spines/x-y labels, not the title).
        ax.set_title(abbrev(organism), loc="left", fontsize=st.FONTSIZE_SMALL, pad=2)


class CoaddPathogenProjectionGridPlot(GridPlot):
    """:data:`default.COADD_MODEL_ID`'s top-N on the same layout: one panel per organism.

    The companion to :class:`PathogenProjectionGridPlot`, on a byte-identical background, so the two
    can be read side by side: same library, same projection, same rank cutoff, but an independent
    CoAdd-trained predictor instead of the per-pathogen ChEMBL models.

    eos3dys has no ``consensus_score``, so its top-N is per ENDPOINT. Where an organism has several
    endpoints (extra strains, or an inhibition and an MIC readout) they share one panel as shades of
    the same hue: they are variants of one organism's activity and should read as that organism at a
    glance, not as unrelated categories. The shade also colours the endpoint's name in the panel, so
    the annotation doubles as the key — a shared legend cannot work here because the endpoint names
    differ from panel to panel.
    """

    def __init__(self, method, background, top_n_table, endpoints, top_n=PROJECTION_TOP_N,
                 cols=COADD_GRID_COLS):
        self.bg_arr, self._bg_extent, self._bg_norm = _background_layer(background)

        items = []
        # unique() keeps config order of first appearance, which is what sets the shade ramp.
        for organism in endpoints["organism"].unique():
            eps = endpoints[endpoints["organism"] == organism]
            colors = shades("crimson", len(eps), floor=COADD_SHADE_FLOOR)
            # Ascending with the shade ramp: darkest endpoint smallest, palest largest. See
            # COADD_SIZE_RANGE.
            sizes = np.linspace(COADD_SIZE_RANGE[0], COADD_SIZE_RANGE[1], len(eps))
            layers = []
            for (_, e), color, size in zip(eps.iterrows(), colors, sizes):
                g = top_n_table[top_n_table["column_name"] == e["column_name"]] \
                    if len(top_n_table) else top_n_table
                layers.append({
                    "label": e["label"], "color": color, "size": float(size),
                    "x": g[f"{method}_x"].to_numpy() if len(g) else np.array([]),
                    "y": g[f"{method}_y"].to_numpy() if len(g) else np.array([]),
                })
            items.append({"organism": organism, "layers": layers})

        self.build_grid(items, cols=cols,
                        name=f"11_{method}_coadd_top{top_n}_pathogens",
                        panel_fn=self._panel, edge_xlabel=f"{method.upper()} 1",
                        edge_ylabel=f"{method.upper()} 2")

    def _panel(self, ax, item, color, xlabel, ylabel):
        _draw_background(ax, self.bg_arr, self._bg_extent, self._bg_norm)
        # Palest layer at the lowest z, darkest last on top — the reverse of the shade order, so the
        # z-index runs down as i runs up. See COADD_SIZE_RANGE.
        n = len(item["layers"])
        for i, layer in enumerate(item["layers"]):
            ax.scatter(layer["x"], layer["y"], s=layer["size"], color=layer["color"],
                       edgecolors="none", zorder=3 + (n - 1 - i))
        ax.set_xticks([])
        ax.set_yticks([])
        ax.text(0.03, 0.96, abbrev(item["organism"]), transform=ax.transAxes,
                ha="left", va="top", fontsize=st.FONTSIZE_SMALL, zorder=LABEL_Z, bbox=LABEL_BBOX)
        for i, layer in enumerate(item["layers"]):
            ax.text(0.03, 0.96 - LABEL_STEP * (i + 1), layer["label"], transform=ax.transAxes,
                    ha="left", va="top", fontsize=st.FONTSIZE_SMALL, color=layer["color"],
                    zorder=LABEL_Z, bbox=LABEL_BBOX)
        st.label(ax, xlabel=xlabel, ylabel=ylabel)


# --------------------------------------------------------------------------- #
# Entry point                                                                  #
# --------------------------------------------------------------------------- #
def _read(output_dir, fname, **kw):
    path = os.path.join(output_dir, fname)
    return pd.read_csv(path, **kw) if os.path.exists(path) else pd.DataFrame()


def save_projection_figures(output_dir, pathogens_csv, top_n=PROJECTION_TOP_N,
                            auroc_matched_table=None):
    """Build the per-pathogen step-11 figures (one per projection method) from the summary CSVs in
    ``output_dir`` and record their footprints in ``figure_cells.json``.

    ``auroc_matched_table`` is :func:`eval_projection.auroc_matched_top_n_per_organism`'s table. When
    given, the UMAP figure is built as :class:`AurocMatchedUmapGridPlot` (AUROC-matrix order/top-N/
    styling) instead of :class:`PathogenProjectionGridPlot` — PCA/t-SNE/TMAP always use the latter
    unchanged.
    """
    pathogens = pd.read_csv(pathogens_csv)
    top_n_table = _read(output_dir, f"11_top{top_n}_per_pathogen.csv")
    footprints = {}
    for method in PROJECTION_METHODS:
        background = _read(output_dir, f"11_{method}_background.csv")
        if method == "umap" and auroc_matched_table is not None:
            plot = AurocMatchedUmapGridPlot(background, auroc_matched_table, method=method)
        else:
            plot = PathogenProjectionGridPlot(method, background, top_n_table, pathogens, top_n=top_n)
        if plot.is_available:
            plot.save(output_dir)
            footprints[plot.name] = list(plot.cells)
            print(f"  figure: {plot.name}")
        else:
            print(f"  [skip figure] 11_{method}_top{top_n}_pathogens: insufficient data")
    # merge, not dump: step 11 writes a second figure family into this same output dir, and a
    # raw dump here would truncate the manifest to whichever family saved last.
    merge_figure_cells(output_dir, footprints)


def save_coadd_projection_figures(output_dir, endpoints, method=COADD_PROJECTION_METHOD,
                                  top_n=PROJECTION_TOP_N):
    """Build the :data:`default.COADD_MODEL_ID` figure and record its footprint.

    ``endpoints`` is :func:`eval_projection.coadd_endpoints`' frame, threaded through by the script
    so the config is parsed once. The background grid is step 11's own
    ``11_{method}_background.csv``, reused rather than recomputed.
    """
    background = _read(output_dir, f"11_{method}_background.csv")
    top_n_table = _read(output_dir, f"11_coadd_top{top_n}_per_endpoint.csv")
    footprints = {}
    plot = CoaddPathogenProjectionGridPlot(method, background, top_n_table, endpoints, top_n=top_n)
    if plot.is_available:
        plot.save(output_dir)
        footprints[plot.name] = list(plot.cells)
        print(f"  figure: {plot.name} ({COADD_MODEL_ID})")
    else:
        print(f"  [skip figure] 11_{method}_coadd_top{top_n}_pathogens: insufficient data")
    merge_figure_cells(output_dir, footprints)
