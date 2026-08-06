"""Step 16 figure — the AUROC matrix with peripheral annotation tracks.

One cell per (predictor, activity endpoint) pair, coloured by AUROC on a DISCRETE 0.1-wide scale,
with colour bands outside the grid carrying the groupings a reader needs to navigate 71 columns and
59 rows: block category, organism class, organism, and predictor model.

Reads only the matrix assembled by :mod:`eval_auroc_matrix`.
"""

import json
import os

import matplotlib.colors as mcolors
import numpy as np
import stylia as st
from matplotlib.cm import ScalarMappable
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1 import make_axes_locatable

from default import (AUROC_MATRIX_BINS, AUROC_MATRIX_CENTER, AUROC_MATRIX_CMAP,
                     OVERLAP_BLANK_DIAGONAL, OVERLAP_MATRIX_BINS, OVERLAP_MATRIX_HUE,
                     ORGANISM_CLASS_ORDER)
from plotting_base import BasePlot
from plotting_colors import INK, REFERENCE_LINE, distinct_colors, hue
from plotting_utils import abbrev, heatmap, sequential_cmap

#: Track band thickness, as a fraction of the main axes, and the gap between bands.
TRACK_SIZE = "2.2%"
TRACK_PAD = 0.035

#: A track cell for a category that does not apply to that row/column (e.g. organism class over the
#: physchem columns) is drawn in the neutral colour — never left to imply membership.
TRACK_NA = REFERENCE_LINE

BLOCK_HUES = {"bioactivity": "cobalt", "cytotoxicity": "crimson",
              "abx resemblance": "amber", "physchem": "turquoise"}

#: Extra pad (points) for the ROW labels. The two left-hand tracks are inserted into exactly the
#: space matplotlib reserves for y tick labels, so without this the bands are drawn over the ends of
#: the labels. Sized to clear both bands plus their pads (2 x 2.2% of the axes width + 2 x TRACK_PAD).
Y_LABEL_PAD = 14

#: Minimum span, in cells, before a track band is labelled in place. Narrower spans rely on the
#: per-row/column tick labels instead of shrinking the text below legibility.
MIN_SPAN_TO_LABEL = 4


def discrete_auroc_cmap(bins=AUROC_MATRIX_BINS, name=AUROC_MATRIX_CMAP,
                        center=AUROC_MATRIX_CENTER):
    """``(ListedColormap, BoundaryNorm)`` for the discrete AUROC scale.

    A DIVERGING ramp pinned to chance: cool below 0.5, near-white at it, warm above. This is the
    right family here because 0.5 is a real neutral with data on both sides — 44 cells rank actives
    below inactives — not merely because it is easy to read.

    The two arms are unequal (0.2-0.5 below, 0.5-1.0 above), so bin midpoints are mapped through a
    TWO-SLOPE normalisation into colormap space: each arm is stretched to its own half. Sampling the
    colormap linearly across the whole range instead would put the white point at 0.6 and quietly
    assert that 0.6 is chance.
    """
    from stylia import DivergingColormap

    base = DivergingColormap(name).cmap
    lo, hi = bins[0], bins[-1]
    positions = []
    for a, b in zip(bins[:-1], bins[1:]):
        mid = (a + b) / 2
        if mid <= center:
            positions.append(0.5 * (mid - lo) / (center - lo))
        else:
            positions.append(0.5 + 0.5 * (mid - center) / (hi - center))
    # Reversed: stylia's crimson_cobalt runs crimson -> cobalt, and the cool-low / warm-high
    # convention is the one a reader applies without being told.
    colors = [base(1.0 - p) for p in positions]
    return mcolors.ListedColormap(colors), mcolors.BoundaryNorm(bins, len(positions), clip=True)


def discrete_overlap_cmap(bins=OVERLAP_MATRIX_BINS, name=OVERLAP_MATRIX_HUE):
    """``(ListedColormap, BoundaryNorm)`` for the discrete Jaccard scale.

    SEQUENTIAL, not diverging: an overlap count's neutral is 0, at the end of the scale rather than
    in the middle, so there is nothing to diverge around. The bins are non-uniform (see
    :data:`default.OVERLAP_MATRIX_BINS`), and colours are sampled at EVEN spacing across the ramp
    rather than at the boundary values — the bins already encode the skew, and spacing the colours by
    value too would make the narrow low bins nearly identical.
    """
    n = len(bins) - 1
    ramp = sequential_cmap(name)
    colors = [ramp(v) for v in np.linspace(0.18, 1.0, n)]
    return mcolors.ListedColormap(colors), mcolors.BoundaryNorm(bins, n, clip=True)


def _blank_diagonal(matrix, n):
    """Copy of ``matrix`` with its leading ``n x n`` diagonal set to NaN.

    Those cells are an entity against ITSELF — 1.0 by construction for both metrics — so they are a
    property of the axes, not a measurement. Blanking keeps them from anchoring the eye along a
    perfect diagonal, and for Jaccard also stops them dominating a data-driven colour scale.
    The true values, diagonal included, stay in the CSV.
    """
    out = matrix.copy()
    for i in range(n):
        out.iloc[i, i] = np.nan
    return out


def _mark_diagonal(ax, n, color=None):
    """Dashed outline on each blanked diagonal cell, so its position stays visible without colour.

    Without the outline a blank diagonal is indistinguishable from missing data; the dashes say
    "deliberately not shown" rather than "no value".
    """
    color = INK if color is None else color
    for i in range(n):
        ax.add_patch(Rectangle((i - 0.5, i - 0.5), 1, 1, fill=False, edgecolor=color,
                               linestyle=(0, (1.5, 1.5)), linewidth=0.5, zorder=5))


def _organism_label(name):
    """Axis label for an organism.

    ``abbrev`` shortens the genus, which is right for a binomial but wrong for a genus-level entry:
    "Campylobacter spp" would become "C. spp", indistinguishable from *C. albicans* in the same
    figure. Genus-level names therefore keep their genus in full.
    """
    return name.split()[0] if name.endswith(" spp") else abbrev(name)


def _spans(values):
    """Contiguous runs in a label sequence, as ``(label, start, stop_exclusive)``."""
    out, start = [], 0
    for i in range(1, len(values) + 1):
        if i == len(values) or values[i] != values[start]:
            out.append((values[start], start, i))
            start = i
    return out


class _MatrixPlotBase(BasePlot):
    """Chrome shared by the two step-14 matrices: peripheral annotation tracks and the discrete
    colourbar. Both draw the same axes with the same groupings and differ only in the quantity in
    the cells, so the tracks live here rather than being duplicated per metric.
    """

    def _track(self, divider, side, values, colors, name, *, vertical=False):
        """One annotation band pinned to the grid edge, plus in-place labels for wide spans.

        ``sharex``/``sharey`` with the main axes is what guarantees the band stays cell-aligned:
        without it the band draws in its own data coordinates and a half-cell offset — which is
        exactly the kind of error that still looks like a plausible figure — becomes possible.
        """
        shared = {"sharey": self.ax} if vertical else {"sharex": self.ax}
        ax = divider.append_axes(side, size=TRACK_SIZE, pad=TRACK_PAD, **shared)
        rgba = np.array([mcolors.to_rgba(TRACK_NA if v is None or (isinstance(v, float)
                                                                   and np.isnan(v))
                                         else colors.get(v, TRACK_NA)) for v in values])
        n = len(values)
        if vertical:
            ax.imshow(rgba.reshape(-1, 1, 4), aspect="auto", interpolation="none",
                      extent=(-0.5, 0.5, n - 0.5, -0.5))
        else:
            ax.imshow(rgba.reshape(1, -1, 4), aspect="auto", interpolation="none",
                      extent=(-0.5, n - 0.5, 0.5, -0.5))
        # NOT set_xticks([]): these axes SHARE their x (or y) with the main grid, and a shared axis
        # shares its locator — clearing ticks here would strip the matrix's own row/column labels.
        # Hiding the axis is per-Axes and leaves the parent's labels intact.
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)
        ax.grid(False)
        for spine in ax.spines.values():
            spine.set_visible(False)

        # Bands wider than MIN_SPAN_TO_LABEL name themselves; narrower ones rely on the per-row and
        # per-column tick labels rather than shrinking text below legibility.
        if not vertical:
            for label, start, stop in _spans(values):
                if label is None or (stop - start) < MIN_SPAN_TO_LABEL:
                    continue
                text = abbrev(label) if name == "organism" else str(label)
                ax.text((start + stop - 1) / 2, 0, text, ha="center", va="center",
                        fontsize=st.FONTSIZE_SMALL, color=INK, clip_on=True)
        return ax

    def _summary_row(self, divider, values, cmap, norm, *, value_fmt, text_light_when,
                     col_labels, label="mean"):
        """A separate one-row band under the grid: the column-wise mean of the matrix.

        Kept OUTSIDE the main axes rather than appended as a 16th row, so it is never mistaken for
        another organism. It shares x with the grid, so it stays column-aligned.

        The diagonal is excluded from the mean — those cells are an entity against itself and are
        blanked in the grid for that reason, so letting them into the average would inflate every
        bioactivity column by a guaranteed maximum.

        This band also carries the column tick labels: appending an axes below the grid puts it
        exactly where the grid's own labels would go, so they are moved down rather than overlapped.
        """
        ax = divider.append_axes("bottom", size="7%", pad=0.06, sharex=self.ax)
        arr = np.asarray(values, dtype=float)
        ax.imshow(arr.reshape(1, -1), cmap=cmap, norm=norm, aspect="auto",
                  interpolation="none", extent=(-0.5, len(arr) - 0.5, 0.5, -0.5))
        for i, v in enumerate(arr):
            if not np.isfinite(v):
                continue
            ax.text(i, 0, value_fmt.format(v), ha="center", va="center",
                    fontsize=st.FONTSIZE_SMALL,
                    color="white" if text_light_when(v) else INK)
        ax.set_yticks([0])
        ax.set_yticklabels([label], fontsize=st.FONTSIZE_SMALL)
        ax.tick_params(axis="y", length=0, pad=Y_LABEL_PAD)
        ax.set_xticks(range(len(arr)))
        ax.set_xticklabels(col_labels, rotation=90, ha="center", fontsize=st.FONTSIZE_SMALL)
        ax.tick_params(axis="x", length=0)
        # stylia's style puts a grid on every new axes; here it strikes through the printed means.
        ax.grid(False)
        ax.set_axisbelow(False)
        # The grid's own column labels move here; leaving both would print them twice.
        self.ax.tick_params(labelbottom=False)
        return ax

    def _colorbar(self, cmap, norm, divider, bins=None, label="AUROC"):
        """Discrete key: one swatch per 0.1 bin, labelled at the boundaries.

        Appended to the grid's RIGHT edge through the same divider as the tracks, so it is pinned
        outside the matrix whatever the axes do. Fixed figure coordinates were tried and broke as
        soon as the swatch legend was dropped and the axes grew into the space. The bottom edge is
        not an option either: the column tick labels run long and a bottom bar lands on them.
        """
        bins = AUROC_MATRIX_BINS if bins is None else bins
        cax = divider.append_axes("right", size="1.2%", pad=0.12)
        # "uniform" so every bin gets equal height in the key: the Jaccard bins are non-uniform, and
        # proportional spacing would squeeze its four low bins into an unreadable sliver.
        cb = self.fig.colorbar(ScalarMappable(cmap=cmap, norm=norm), cax=cax,
                               orientation="vertical", ticks=bins, spacing="uniform")
        cb.set_label(label, fontsize=st.FONTSIZE_SMALL)
        cb.ax.tick_params(labelsize=st.FONTSIZE_SMALL)
        # The scale is clipped below its first boundary, so that tick is an inequality: everything
        # from the matrix minimum up to 0.5 sits in this bin.
        fmt = "{:.2f}" if max(bins) <= 1.0 and min(np.diff(bins)) < 0.05 else "{:.1f}"
        ticks = [f"\u2264{fmt.format(bins[0])}"] + [fmt.format(b) for b in bins[1:]]
        cb.ax.set_yticklabels(ticks)



class AurocMatrixPlot(_MatrixPlotBase):
    """The 59 x 71 AUROC matrix, with four annotation tracks on the top edge and two on the left."""

    #: Landscape footprint (3 rows x 6 cols = 180 x 90 mm). At 15 rows the labels still have room,
    #: and the column tick labels take their space below the grid rather than inside it.
    def __init__(self, matrix, rows, cols, ax=None, cells=(3, 6)):
        super().__init__(ax=ax, cells=cells)
        self.name = "14_auroc_matrix"
        if not matrix.size:
            self._unavailable()
            return

        cmap, norm = discrete_auroc_cmap()
        n_bio = len(rows)
        shown = _blank_diagonal(matrix, n_bio)

        # --- Category -> colour, built once so the top and left tracks always agree ---
        # Classes take the LIGHTENED pass of the categorical palette. Block category and organism
        # class share one legend, and at the base level `bioactivity` and `Gram-negative bacteria`
        # both land on cobalt — two different meanings in one key, in the same colour.
        class_colors = dict(zip(ORGANISM_CLASS_ORDER,
                                distinct_colors(len(ORGANISM_CLASS_ORDER), levels=(0.55, 0.3))))
        model_ids = list(dict.fromkeys(cols["model_id"]))
        model_colors = dict(zip(model_ids, distinct_colors(len(model_ids))))

        # --- Column-axis category sequences (bioactivity span, then the property blocks) ---
        col_block = ["bioactivity"] * n_bio + cols["block"].tolist()
        col_class = rows["organism_class"].tolist() + [None] * len(cols)
        col_model = [None] * n_bio + cols["model_id"].tolist()

        labels = [_organism_label(o) for o in rows["organism"]]
        col_labels = labels + cols["column_name"].tolist()

        # Values are printed: at 15 x 27 each cell is ~6 mm and carries its AUROC. On a DIVERGING
        # scale the saturated cells are at BOTH ends, so light text is keyed to distance from chance
        # rather than to magnitude — white in the middle of the scale would be invisible.
        heatmap(self.ax, shown, cmap=cmap, norm=norm, annotate=True, value_fmt="{:.2f}",
                text_light_when=lambda v: abs(v - AUROC_MATRIX_CENTER) >= 0.35,
                nan_color=TRACK_NA,
                x_rotation=90, row_labels=labels, col_labels=col_labels,
                annot_fontsize=st.FONTSIZE_SMALL)
        _mark_diagonal(self.ax, n_bio)
        self.ax.tick_params(labelsize=st.FONTSIZE_SMALL, length=0)
        # stylia's article style draws a grid; over a heatmap it lands on top of the cells and
        # strikes through the printed values. A filled matrix needs no gridlines.
        self.ax.grid(False)
        self.ax.set_axisbelow(False)
        self.ax.tick_params(axis="y", pad=Y_LABEL_PAD)

        # stylia's blank figure ships placeholder axis labels; a matrix labels every row and column
        # individually and needs neither.
        self.ax.set_xlabel("")
        self.ax.set_ylabel("")

        divider = make_axes_locatable(self.ax)
        # Appended innermost-first, so the LAST call ends up outermost. Order runs fine -> coarse
        # going outward: model, organism, class, block.
        # No organism track: every row and every bioactivity column now IS an organism and is
        # named on the axis, so a 15-hue band would only repeat the label.
        self._track(divider, "top", col_model, model_colors, "model")
        self._track(divider, "top", col_class, class_colors, "class")
        self._track(divider, "top", col_block, {b: hue(h) for b, h in BLOCK_HUES.items()}, "block")
        self._track(divider, "left", rows["organism_class"].tolist(), class_colors, "class",
                    vertical=True)
        self._summary_row(divider, shown.mean(axis=0, skipna=True), cmap, norm,
                          value_fmt="{:.2f}",
                          text_light_when=lambda v: abs(v - AUROC_MATRIX_CENTER) >= 0.35,
                          col_labels=col_labels)

        # No swatch legend (user-directed). The block track labels itself in place, and organism
        # class is inferable from the organism names on both axes — but note that the narrow class
        # bands (Fungi, Protozoa, Helminths, Mycobacteria) are now unkeyed colour only.
        self._colorbar(cmap, norm, divider)

    # ---------------------------------------------------------------- tracks
# --------------------------------------------------------------------------- #
# Entry point                                                                  #
# --------------------------------------------------------------------------- #
def save_auroc_matrix_figure(output_dir, matrix, rows, cols):
    """Build the step-14 figure and record its footprint in ``figure_cells.json``."""
    plot = AurocMatrixPlot(matrix, rows, cols)
    footprints = {}
    if plot.is_available:
        plot.save(output_dir)
        footprints[plot.name] = list(plot.cells)
        print(f"  figure: {plot.name}")
    else:
        print("  [skip figure] 14_auroc_matrix: empty matrix")
    with open(os.path.join(output_dir, "figure_cells.json"), "w") as f:
        json.dump(footprints, f, indent=2)
    return footprints


class OverlapMatrixPlot(_MatrixPlotBase):
    """Top-N Jaccard overlap between each organism's actives and every column's actives.

    Same axes and the same annotation tracks as :class:`AurocMatrixPlot`, different quantity: how
    many of the row's 1000 actives are among the column's own top 1000. The bioactivity block is
    symmetric here, which the AUROC block is not.

    The self-overlap diagonal is 1000/1000 by construction and is BLANKED, so the colour scale is set
    by real overlaps (max 724) rather than by a cell that cannot be anything else.
    """

    def __init__(self, matrix, rows, cols, ax=None, cells=(3, 6),
                 blank_diagonal=OVERLAP_BLANK_DIAGONAL):
        super().__init__(ax=ax, cells=cells)
        self.name = "14_overlap_matrix"
        if not matrix.size:
            self._unavailable()
            return

        n_bio = len(rows)
        shown = _blank_diagonal(matrix, n_bio) if blank_diagonal else matrix
        cmap, norm = discrete_overlap_cmap()

        class_colors = dict(zip(ORGANISM_CLASS_ORDER,
                                distinct_colors(len(ORGANISM_CLASS_ORDER), levels=(0.55, 0.3))))
        model_ids = list(dict.fromkeys(cols["model_id"]))
        model_colors = dict(zip(model_ids, distinct_colors(len(model_ids))))

        col_block = ["bioactivity"] * n_bio + cols["block"].tolist()
        col_class = rows["organism_class"].tolist() + [None] * len(cols)
        col_model = [None] * n_bio + cols["model_id"].tolist()

        labels = [_organism_label(o) for o in rows["organism"]]
        col_labels = labels + cols["column_name"].tolist()

        # Whole numbers: the cell is a count of compounds out of ACTIVITY_BINARIZE_TOP_N.
        heatmap(self.ax, shown, cmap=cmap, norm=norm, annotate=True, value_fmt="{:.0f}",
                text_light_when=lambda v: v >= OVERLAP_MATRIX_BINS[-3], nan_color=TRACK_NA,
                x_rotation=90, row_labels=labels, col_labels=col_labels,
                annot_fontsize=st.FONTSIZE_SMALL)
        if blank_diagonal:
            _mark_diagonal(self.ax, n_bio)
        self.ax.tick_params(labelsize=st.FONTSIZE_SMALL, length=0)
        self.ax.tick_params(axis="y", pad=Y_LABEL_PAD)
        self.ax.grid(False)
        self.ax.set_axisbelow(False)
        self.ax.set_xlabel("")
        self.ax.set_ylabel("")

        divider = make_axes_locatable(self.ax)
        self._track(divider, "top", col_model, model_colors, "model")
        self._track(divider, "top", col_class, class_colors, "class")
        self._track(divider, "top", col_block, {b: hue(h) for b, h in BLOCK_HUES.items()}, "block")
        self._track(divider, "left", rows["organism_class"].tolist(), class_colors, "class",
                    vertical=True)
        self._summary_row(divider, shown.mean(axis=0, skipna=True), cmap, norm,
                          value_fmt="{:.0f}",
                          text_light_when=lambda v: v >= OVERLAP_MATRIX_BINS[-3],
                          col_labels=col_labels)
        self._colorbar(cmap, norm, divider, bins=OVERLAP_MATRIX_BINS,
                       label="Shared actives (of 1000)")


def save_overlap_matrix_figure(output_dir, matrix, rows, cols, footprints=None):
    """The shared-actives figure, appended to ``figure_cells.json``."""
    plot = OverlapMatrixPlot(matrix, rows, cols)
    footprints = {} if footprints is None else footprints
    if plot.is_available:
        plot.save(output_dir)
        footprints[plot.name] = list(plot.cells)
        print(f"  figure: {plot.name}")
    else:
        print("  [skip figure] 14_overlap_matrix: empty matrix")
    with open(os.path.join(output_dir, "figure_cells.json"), "w") as f:
        json.dump(footprints, f, indent=2)
    return footprints
