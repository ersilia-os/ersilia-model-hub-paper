"""Step 10 figure — the AUROC matrix with peripheral annotation tracks.

One cell per (predictor, activity endpoint) pair, coloured by AUROC on a DISCRETE 0.1-wide scale,
with colour bands outside the grid carrying the groupings a reader needs to navigate 71 columns and
59 rows: block category, organism class, organism, and predictor model.

Reads only the matrix assembled by :mod:`eval_auroc_matrix`.
"""

import os

import matplotlib.colors as mcolors
import numpy as np
import stylia as st
from matplotlib.cm import ScalarMappable
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1 import make_axes_locatable

from default import (ACTIVITY_BINARIZE_TOP_N, AUROC_MATRIX_BINS, AUROC_MATRIX_CENTER,
                     AUROC_MATRIX_CMAP, ORGANISM_CLASS_ORDER, OVERLAP_BLANK_DIAGONAL,
                     OVERLAP_MATRIX_BINS, OVERLAP_MATRIX_SPECTRUM)
from plotting_base import BasePlot
from plotting_colors import INK, REFERENCE_LINE, distinct_colors
from plotting_utils import abbrev, heatmap, merge_figure_cells, spectrum_cmap

#: Track band thickness, as a fraction of the main axes, and the gap between bands.
TRACK_SIZE = "2.2%"
TRACK_PAD = 0.035

#: A track cell for a category that does not apply to that row/column (e.g. organism class over the
#: physchem columns) is drawn in the neutral colour — never left to imply membership.
TRACK_NA = REFERENCE_LINE

#: Display label per property family on the AUROC/overlap matrix column axis (user-directed,
#: 2026-09-02) — since the cytotox/abx merge each family is a single column named "rank_sum" (see
#: default.AUROC_MATRIX_BLOCKS), so the axis shows this instead of the bare family key or a
#: "family_rank_sum" concatenation. "abx sim" rather than "abx": the merged column folds in
#: abx_score alongside the two similarity counts, so "resemblance"/"sim" reads more accurately than
#: the bare family key would.
PREDICTOR_FAMILY_LABELS = {"cytotox": "cytotox", "abx": "abx sim"}

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


def overlap_bins(top_n, base=OVERLAP_MATRIX_BINS, base_n=ACTIVITY_BINARIZE_TOP_N):
    """:data:`default.OVERLAP_MATRIX_BINS` rescaled from ``base_n`` to ``top_n``.

    The published bins are hand-tuned to top-1000 (counts 0-724, median 3), so reusing them at
    another cutoff would waste the scale: ``BoundaryNorm(clip=True)`` would collapse a top-100 matrix
    into the lowest two bins and saturate a top-10000 one in the highest. Scaling by ``top_n /
    base_n`` keeps each bin at the same FRACTION of the cutoff, which is what makes the three panels
    comparable — a cell in the third bin means the same thing on all of them.

    Two boundaries are pinned rather than scaled. ``0`` is the floor, and ``1`` is the informative
    one: it separates cells sharing NO compound at all from cells sharing something, a distinction
    that exists at every cutoff and would be scaled away (at top-100 the second bin would land on 1
    and collide with it). After rounding, the remaining boundaries are forced strictly increasing —
    at top-100 the naive scale gives 10 -> 1, 25 -> 2.5, which round into each other.

    Returns a list of ints. ``overlap_bins(base_n)`` is ``base`` unchanged.
    """
    scale = top_n / base_n
    out = [int(base[0]), int(base[1])]          # 0 and 1, pinned
    for b in base[2:]:
        out.append(max(int(round(b * scale)), out[-1] + 1))
    if list(out) != sorted(set(out)):
        raise ValueError(f"overlap_bins({top_n}) is not strictly increasing: {out}")
    return out


def discrete_overlap_cmap(bins=OVERLAP_MATRIX_BINS, hues=OVERLAP_MATRIX_SPECTRUM):
    """``(ListedColormap, BoundaryNorm)`` for the discrete Jaccard scale.

    SEQUENTIAL, not diverging: an overlap count's neutral is 0, at the end of the scale rather than
    in the middle, so there is nothing to diverge around. The bins are non-uniform (see
    :data:`default.OVERLAP_MATRIX_BINS`), and colours are sampled at EVEN spacing across the ramp
    rather than at the boundary values — the bins already encode the skew, and spacing the colours by
    value too would make the narrow low bins nearly identical.

    **User-directed, 2026-08-31:** uses :func:`plotting_utils.spectrum_cmap` (white → several house
    hues, see :data:`default.OVERLAP_MATRIX_SPECTRUM`) rather than the plain 2-stop
    :func:`plotting_utils.sequential_cmap` — a single hue's white-to-base range left the higher bins
    of an 8-bin scale hard to tell apart. Sampled from **0.0**, not the usual legibility-margin offset
    other sequential scales in this repo use (see :func:`plotting_colors.count_shades`'s
    ``headroom``): the 0-shared-actives bin is meant to read as true, blank white here, not a pale
    tint — 0 is a common, meaningful value for this matrix (76.2% of bioactivity off-diagonal cells
    at top-100, 26.7% at top-1000, 0% at top-10000), not an edge case to keep visible.
    """
    n = len(bins) - 1
    ramp = spectrum_cmap(hues)
    colors = [ramp(v) for v in np.linspace(0.0, 1.0, n)]
    return mcolors.ListedColormap(colors), mcolors.BoundaryNorm(bins, n, clip=True)


def _dark_threshold(cmap, vmax, target_luminance=0.5, resolution=256):
    """The value at which ``cmap`` (sampled over ``[0, vmax]``) first drops below
    ``target_luminance`` — the point past which white cell-annotation text reads better than ink.

    Computed from the colormap itself rather than guessed as a fraction of ``vmax``, so it stays
    correct if the palette or the scale's maximum ever changes — relevant here because
    :func:`plotting_utils.spectrum_cmap` mixes several house hues of different intrinsic lightness,
    unlike a single-hue ramp where "darker" and "higher value" track each other exactly.
    """
    t = np.linspace(0.0, 1.0, resolution)
    rgba = np.array([cmap(v) for v in t])
    luminance = 0.2126 * rgba[:, 0] + 0.7152 * rgba[:, 1] + 0.0722 * rgba[:, 2]
    below = np.flatnonzero(luminance < target_luminance)
    frac = t[below[0]] if len(below) else 1.0
    return frac * vmax


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
    """Chrome shared by the two step-10 matrices: peripheral annotation tracks and the discrete
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

    def _colorbar(self, cmap, norm, divider, bins=None, label="AUROC", continuous=False):
        """Discrete key: one swatch per 0.1 bin, labelled at the boundaries.

        Appended to the grid's RIGHT edge through the same divider as the tracks, so it is pinned
        outside the matrix whatever the axes do. Fixed figure coordinates were tried and broke as
        soon as the swatch legend was dropped and the axes grew into the space. The bottom edge is
        not an option either: the column tick labels run long and a bottom bar lands on them.

        ``continuous=True`` (user-directed, 2026-09-02 \u2014 step 10's overlap matrix) skips the
        discrete-bin machinery entirely: a plain linear colourbar with matplotlib's own tick
        placement, matching a ``Normalize``d ``norm`` rather than a ``BoundaryNorm``. The default
        (``False``) keeps the original discrete, evenly-spaced-per-bin key used by the AUROC matrix
        and by every OTHER caller of this method \u2014 unchanged.
        """
        cax = divider.append_axes("right", size="1.2%", pad=0.12)
        if continuous:
            cb = self.fig.colorbar(ScalarMappable(cmap=cmap, norm=norm), cax=cax,
                                   orientation="vertical")
            cb.set_label(label, fontsize=st.FONTSIZE_SMALL)
            cb.ax.tick_params(labelsize=st.FONTSIZE_SMALL)
            return

        bins = AUROC_MATRIX_BINS if bins is None else bins
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
    """The 15 x 17 AUROC matrix, with one annotation track on the left."""

    #: Near-square footprint. Narrowed 2026-09-02 (user-directed) from the 27-column-era (3, 6)
    #: landscape shape once the property blocks collapsed to cytotox+abx sim (17 columns total) —
    #: that width made every cell a wide rectangle instead of a square. `aspect="equal"` on the
    #: heatmap call below enforces literally square cells regardless of any residual mismatch.
    def __init__(self, matrix, rows, cols, ax=None, cells=(3, 3.6), name=None):
        super().__init__(ax=ax, cells=cells)
        # Overridable so a variant run (e.g. step 10's own non-abx section, the same matrix over a
        # filtered library) does not write a figure claiming to be the main one.
        self.name = name or "10_auroc_matrix"
        if not matrix.size:
            self._unavailable()
            return

        cmap, norm = discrete_auroc_cmap()
        n_bio = len(rows)
        shown = _blank_diagonal(matrix, n_bio)

        # Classes take the LIGHTENED pass of the categorical palette, for the one surviving
        # (left, row) track — see the class docstring for why the top tracks are gone.
        class_colors = dict(zip(ORGANISM_CLASS_ORDER,
                                distinct_colors(len(ORGANISM_CLASS_ORDER), levels=(0.55, 0.3))))

        labels = [_organism_label(o) for o in rows["organism"]]
        # PREDICTOR_FAMILY_LABELS, not raw column_name: since the 2026-09-02 cytotox/abx merge every
        # property family is a single column literally named "rank_sum" (default.AUROC_MATRIX_BLOCKS)
        # — the family's own display label is what actually distinguishes the two on the axis.
        col_labels = labels + [PREDICTOR_FAMILY_LABELS.get(f, f) for f in cols["family"]]

        # Values are printed: at 15 x 27 each cell is ~6 mm and carries its AUROC. On a DIVERGING
        # scale the saturated cells are at BOTH ends, so light text is keyed to distance from chance
        # rather than to magnitude — white in the middle of the scale would be invisible.
        heatmap(self.ax, shown, cmap=cmap, norm=norm, annotate=True, value_fmt="{:.2f}",
                text_light_when=lambda v: abs(v - AUROC_MATRIX_CENTER) >= 0.35,
                nan_color=TRACK_NA,
                x_rotation=90, row_labels=labels, col_labels=col_labels,
                annot_fontsize=st.FONTSIZE_SMALL, aspect="equal")
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
        # The top model/class/block tracks were dropped 2026-09-02 (user-directed, "not so
        # informative now"): with only 2 property columns left post-merge, "model" was always
        # "merged", "class" over the property columns was always blank (TRACK_NA), and "block" only
        # ever separated 15 bioactivity columns from 2 property ones — a distinction the
        # PREDICTOR_FAMILY_LABELS column labels already make directly. Only the row-side class band
        # (below) still carries information a reader cannot already get from the axis labels.
        self._track(divider, "left", rows["organism_class"].tolist(), class_colors, "class",
                    vertical=True)
        self._summary_row(divider, shown.mean(axis=0, skipna=True), cmap, norm,
                          value_fmt="{:.2f}",
                          text_light_when=lambda v: abs(v - AUROC_MATRIX_CENTER) >= 0.35,
                          col_labels=col_labels)

        # No swatch legend (user-directed). Organism class is inferable from the organism names on
        # both axes — but note that the narrow class bands (Fungi, Protozoa, Helminths, Mycobacteria)
        # are unkeyed colour only.
        self._colorbar(cmap, norm, divider)
# --------------------------------------------------------------------------- #
# Entry point                                                                  #
# --------------------------------------------------------------------------- #
def save_auroc_matrix_figure(output_dir, matrix, rows, cols, name=None):
    """Build the AUROC matrix figure and merge its footprint into ``figure_cells.json``.

    Routes through :func:`plotting_utils.merge_figure_cells` rather than writing the manifest
    directly: step 10 now draws more than one AUROC-matrix figure into the same output dir (the
    baseline order plus row-order comparison variants), and writing with ``"w"`` here would let
    whichever call ran last truncate the manifest to its own single entry — the exact bug already
    fixed for step 12's abx grid.
    """
    plot = AurocMatrixPlot(matrix, rows, cols, name=name)
    footprints = {}
    if plot.is_available:
        plot.save(output_dir)
        footprints[plot.name] = list(plot.cells)
        print(f"  figure: {plot.name}")
    else:
        print(f"  [skip figure] {plot.name}: empty matrix")
    return merge_figure_cells(output_dir, footprints)


class OverlapMatrixPlot(_MatrixPlotBase):
    """Top-N Jaccard overlap between each organism's actives and every column's actives.

    Same axes and the same annotation tracks as :class:`AurocMatrixPlot`, different quantity: how
    many of the row's ``top_n`` actives are among the column's own top ``top_n``. The bioactivity
    block is symmetric here, which the AUROC block is not.

    The self-overlap diagonal is ``top_n``/``top_n`` by construction and is BLANKED, so the colour
    scale is set by real overlaps rather than by a cell that cannot be anything else.

    ``top_n`` drives the figure name, the colour bins (:func:`overlap_bins`) and the colourbar
    label, so a different cutoff produces an independent figure on a comparable scale. The class
    stays generic over ``top_n`` for that reason, even though the pipeline currently only ever
    draws it at the single ``ACTIVITY_BINARIZE_TOP_N`` cutoff (1000) — the 100 and 10000 cutoffs
    this view was once compared at were dropped 2026-09-01 (see `scripts/README.md`). **They were
    not interchangeable**: chance overlap for two random top-N sets is ``top_n^2 / n_total``, which
    at top-10000 exceeds the first non-zero bin, so a handful of shared actives reads as *below*
    chance there and as meaningful at top-100.
    """

    #: Near-square footprint (see AurocMatrixPlot's, same rationale — narrowed 2026-09-02 from the
    #: 27-column-era (3, 6) once the property blocks collapsed to cytotox+abx sim).
    def __init__(self, matrix, rows, cols, ax=None, cells=(3, 3.6),
                 blank_diagonal=OVERLAP_BLANK_DIAGONAL, top_n=ACTIVITY_BINARIZE_TOP_N,
                 name=None, continuous_color=False):
        super().__init__(ax=ax, cells=cells)
        # The cutoff is in the name so that a figure drawn at a different top_n cannot overwrite
        # this one's PNG/PDF or collide with it in figure_cells.json.
        self.name = name or f"10_overlap_matrix_top{top_n}"
        if not matrix.size:
            self._unavailable()
            return

        n_bio = len(rows)
        shown = _blank_diagonal(matrix, n_bio) if blank_diagonal else matrix
        # `continuous_color` (user-directed, 2026-09-02): the DISCRETE non-uniform bins below were
        # tuned to a heavily skewed distribution and drawn with `spacing="uniform"` in the colourbar
        # — equal visual height per bin regardless of its numeric width, which is what a LOG-scale
        # colourbar looks like even though nothing here is actually log-transformed. Step 10's
        # top-1000 figure — both the main one and its non-abx counterpart, since the latter was
        # folded into step 10 to draw the exact same plots (2026-09-02) — opts into a plain
        # continuous linear scale instead (same OVERLAP_MATRIX_SPECTRUM palette, just
        # `Normalize(0, top_n)` instead of `BoundaryNorm`). The discrete non-uniform-bin default
        # remains available for any future caller that still needs the original fix for this skew
        # (see OVERLAP_MATRIX_BINS).
        bins = overlap_bins(top_n)
        if continuous_color:
            cmap, norm = spectrum_cmap(OVERLAP_MATRIX_SPECTRUM), mcolors.Normalize(vmin=0, vmax=top_n)
            # Read the switch-to-white point off the ACTUAL colormap (luminance crossing 0.5) rather
            # than guessing a fraction of top_n — stays correct if the palette or top_n changes.
            dark_at = _dark_threshold(cmap, top_n)
            text_light_when = lambda v: v >= dark_at
        else:
            cmap, norm = discrete_overlap_cmap(bins=bins)
            text_light_when = lambda v: v >= bins[-3]

        class_colors = dict(zip(ORGANISM_CLASS_ORDER,
                                distinct_colors(len(ORGANISM_CLASS_ORDER), levels=(0.55, 0.3))))

        labels = [_organism_label(o) for o in rows["organism"]]
        # PREDICTOR_FAMILY_LABELS, not raw column_name: since the 2026-09-02 cytotox/abx merge every
        # property family is a single column literally named "rank_sum" — the family's own display
        # label is what actually distinguishes the two on the axis.
        col_labels = labels + [PREDICTOR_FAMILY_LABELS.get(f, f) for f in cols["family"]]

        # Whole numbers: the cell is a count of compounds out of ACTIVITY_BINARIZE_TOP_N.
        heatmap(self.ax, shown, cmap=cmap, norm=norm, annotate=True, value_fmt="{:.0f}",
                text_light_when=text_light_when, nan_color=TRACK_NA,
                x_rotation=90, row_labels=labels, col_labels=col_labels,
                annot_fontsize=st.FONTSIZE_SMALL, aspect="equal")
        if blank_diagonal:
            _mark_diagonal(self.ax, n_bio)
        self.ax.tick_params(labelsize=st.FONTSIZE_SMALL, length=0)
        self.ax.tick_params(axis="y", pad=Y_LABEL_PAD)
        self.ax.grid(False)
        self.ax.set_axisbelow(False)
        self.ax.set_xlabel("")
        self.ax.set_ylabel("")

        # The top model/class/block tracks were dropped 2026-09-02 (user-directed, "not so
        # informative now" — see AurocMatrixPlot). Only the row-side class band still carries
        # information a reader cannot already get from the axis labels.
        divider = make_axes_locatable(self.ax)
        self._track(divider, "left", rows["organism_class"].tolist(), class_colors, "class",
                    vertical=True)
        self._summary_row(divider, shown.mean(axis=0, skipna=True), cmap, norm,
                          value_fmt="{:.0f}",
                          text_light_when=text_light_when,
                          col_labels=col_labels)
        self._colorbar(cmap, norm, divider, bins=bins,
                       label=f"Shared actives (of {top_n})", continuous=continuous_color)


def save_overlap_matrix_figure(output_dir, matrix, rows, cols, footprints=None,
                               top_n=ACTIVITY_BINARIZE_TOP_N, name=None, continuous_color=False):
    """One shared-actives figure, merged into ``figure_cells.json``.

    Called once per cutoff. Routes through :func:`plotting_utils.merge_figure_cells` — like
    :func:`save_auroc_matrix_figure` above — rather than writing the manifest directly, so a call on
    its own (e.g. re-running just this figure) cannot truncate entries another figure already wrote.
    The ``footprints`` parameter is accepted for backward compatibility with callers that still
    thread it through by hand; it is otherwise unused.
    """
    plot = OverlapMatrixPlot(matrix, rows, cols, top_n=top_n, name=name,
                             continuous_color=continuous_color)
    new = {}
    if plot.is_available:
        plot.save(output_dir)
        new[plot.name] = list(plot.cells)
        print(f"  figure: {plot.name}")
    else:
        print(f"  [skip figure] {plot.name}: empty matrix")
    return merge_figure_cells(output_dir, new)


def save_dendrogram_figure(output_dir, linkages, xlabel, name):
    """One horizontal dendrogram per class with 2+ organisms, reading ``linkages`` from
    :func:`eval_auroc_matrix.phylogeny_class_linkages` — generic over what the merge height MEANS
    (passed in via ``xlabel``) so it isn't tied to one particular tree-building method. An earlier
    hierarchical-clustering-on-AUROC-profile alternative used this same function; it was dropped
    2026-09-02 in favour of the taxonomy tree only.

    A diagnostic, not a :class:`BasePlot` on the 3 cm cell grid — same departure as step 08's
    ``pathogen_jaccard_figure`` (:mod:`plots_matrix_analyses`): sized in inches by leaf count, and
    not entered into ``figure_cells.json``. The 4 singleton classes (Mycobacteria, Fungi, Protozoa,
    Helminths at 15 organisms) have no tree to draw and are simply absent from ``linkages``. Reads
    ``Z`` directly from the SAME function that ordered the corresponding heatmap's rows, so a
    dendrogram can never show a different tree than the one its heatmap was permuted by.
    """
    import matplotlib.pyplot as plt
    from scipy.cluster.hierarchy import dendrogram

    if not linkages:
        print(f"  [skip figure] {name}: no class has more than one organism to cluster")
        return None

    classes = list(linkages.keys())
    heights = [len(members) for _, members in linkages.values()]
    fig, axes = plt.subplots(len(classes), 1, figsize=(5, sum(h * 0.35 + 0.55 for h in heights)),
                             gridspec_kw={"height_ratios": heights}, squeeze=False)
    for ax, cls in zip(axes[:, 0], classes):
        z, members = linkages[cls]
        dendrogram(z, ax=ax, orientation="left", labels=[_organism_label(o) for o in members],
                  color_threshold=0, above_threshold_color=INK)
        # scipy draws leaf 0 at the BOTTOM for orientation="left" (leaf n-1 at the top) — the
        # opposite of `members`' top-to-bottom order and of the heatmap's row order. Without this,
        # the dendrogram reads bottom-to-top while its heatmap reads top-to-bottom, and the two
        # can't be visually cross-referenced by position.
        ax.invert_yaxis()
        ax.set_title(cls, fontsize=st.FONTSIZE_SMALL, loc="left")
        ax.tick_params(labelsize=st.FONTSIZE_SMALL, length=0)
        for spine in ("top", "right", "bottom"):
            ax.spines[spine].set_visible(False)

    axes[-1, 0].set_xlabel(xlabel, fontsize=st.FONTSIZE_SMALL)
    fig.tight_layout()

    png_dir = os.path.join(output_dir, "png")
    pdf_dir = os.path.join(output_dir, "pdf")
    os.makedirs(png_dir, exist_ok=True)
    os.makedirs(pdf_dir, exist_ok=True)
    png_path = os.path.join(png_dir, f"{name}.png")
    pdf_path = os.path.join(pdf_dir, f"{name}.pdf")
    fig.savefig(png_path, dpi=300)
    fig.savefig(pdf_path)
    plt.close(fig)
    print(f"  figure: {name} (plain matplotlib diagnostic, not in figure_cells.json)")
    return png_path, pdf_path
