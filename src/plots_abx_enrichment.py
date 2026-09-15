"""Step 12 figure — the eos19mt Fisher-enrichment matrix.

One cell per (eos19mt antibiotic class, pathogen) pair, from :mod:`eval_abx_enrichment`'s wide
odds-ratio/p-value matrices. Colour is a diverging scale on log2(odds ratio), centred on 0 (odds
ratio == 1, no association) — cool where a class is under-represented in the pathogen's top-N,
warm where it is over-represented, matching the AUROC matrix's cool-low/warm-high convention.
Some pairs have an odds ratio of exactly 0 or infinity (e.g. a class with a single non-zero
compound library-wide, which lands either entirely inside or entirely outside a pathogen's top-N)
— the colour scale is clipped at :data:`default.ABX_ENRICHMENT_LOG2OR_CAP` so these cells still
render at the scale's extreme rather than breaking the colormap, but the printed cell text is
always the TRUE odds ratio (``0``/``inf``), never the clipped colour value.

A trailing ``*``/``**``/``***`` marks the raw (unadjusted) per-cell p-value at
:data:`default.ABX_ENRICHMENT_SIG_THRESHOLDS` — a display convenience, not a multiple-testing
claim (the long-format CSV carries a Benjamini-Hochberg-adjusted column for that).
"""

import matplotlib.colors as mcolors
import numpy as np
import stylia as st
from matplotlib.patches import Circle
from mpl_toolkits.axes_grid1 import make_axes_locatable

from default import (ABX_ENRICHMENT_LOG2OR_CAP, ABX_ENRICHMENT_OR_SIZE_CAP,
                     ABX_ENRICHMENT_SIG_THRESHOLDS, ORGANISM_CLASS_ORDER)
from plotting_base import BasePlot
from plotting_colors import INK, REFERENCE_LINE, distinct_colors
from plotting_utils import diverging_cmap, heatmap, merge_figure_cells, nested_size_legend, sentence_case
from plots_auroc_matrix import _MatrixPlotBase, _organism_label

#: Near-full-page footprint: 38 rows need a lot of vertical room, and the cell text (odds ratios up
#: to 5 digits plus significance asterisks) needs the full print width's worth of column room too —
#: tuned visually against the rendered figure (an initial (7.2, 3) let neighbouring cells' text run
#: into each other), not derived from the row/column counts.
CELLS = (7.2, 6)

#: Bubble-grid variant shares these across its main panel and its key panel.
BUBBLE_EDGE_COLOR = INK
BUBBLE_EDGE_LINEWIDTH = 0.5
#: No in-cell text to fit here (colour + size + a separate key panel carry the value instead), so the
#: pitch only has to stay legible for bubbles/labels rather than for fitting digits — tuned visually
#: against the rendered figure, same as CELLS above, still a real reduction from the heatmap's
#: (7.2, 6).
BUBBLE_CELLS = (6.0, 3.2)


def _or_text(v):
    """The true odds ratio as display text — ``inf``/``0`` shown as such, never clipped."""
    if np.isnan(v):
        return ""
    if np.isinf(v):
        return "inf" if v > 0 else "0"
    if v >= 100:
        return f"{v:.0f}"
    if v >= 10:
        return f"{v:.1f}"
    return f"{v:.2f}"


def _sig_stars(p, thresholds=ABX_ENRICHMENT_SIG_THRESHOLDS):
    """``*``/``**``/``***`` for the raw p-value at each successive threshold, ``""`` if not
    significant or missing (degenerate class).
    """
    if not np.isfinite(p):
        return ""
    return "*" * sum(p < t for t in thresholds)


#: Tokens in a raw eos19mt class name (config/antibiotic_resemblance.csv's column_name, snake_case)
#: that are genuinely shortened forms rather than dictionary words — spelled out with a trailing dot
#: ("der" -> "derivative(s)", "cpds" -> "compounds") or as the acronym they actually are ("lc"/"mc"
#: -> long-/medium-chain). Everything else in these names (amino, acid, all, quat, ...) already reads
#: as a real word once de-snake-cased, so it is left alone rather than guessed at.
CLASS_LABEL_ABBREVIATIONS = {"der": "der.", "cpds": "cpds.", "lc": "LC", "mc": "MC"}


def _class_label(name):
    """An eos19mt class name as natural-language axis text: underscores become spaces, known
    abbreviated tokens are rewritten (:data:`CLASS_LABEL_ABBREVIATIONS`), and the result is
    sentence-cased (:func:`plotting_utils.sentence_case`) rather than left in snake_case.
    """
    words = [CLASS_LABEL_ABBREVIATIONS.get(w, w) for w in name.split("_")]
    return sentence_case(" ".join(words))


def _bubble_area(odds_ratio, *, min_area, max_area, cap):
    """Marker area (pt^2) for one cell's TRUE odds ratio, for :class:`AbxEnrichmentBubblePlot`.

    ``0`` (tested, no association) maps to ``min_area`` — small but visible, distinguishing it from
    a degenerate untested cell (NaN, no bubble drawn at all — the caller's responsibility to skip
    before calling this). Everything from just above 0 up to ``cap`` (including ``+inf``, the one
    all-in-or-out-of-top-N cell) is linear IN AREA onto ``[min_area, max_area]`` — area, not
    diameter, is what a reader's eye integrates as "size" — then clips above ``cap``, the same
    convention the colour scale's own ``clip(-cap, cap)`` uses.
    """
    if odds_ratio <= 0:
        return min_area
    frac = min(float(odds_ratio), cap) / cap
    return min_area + frac * (max_area - min_area)


def _sig_edge_legend(ax, x0, y, *, radius, spacing=0.3, color=None):
    """Two example circles — filled vs unfilled-with-a-coloured-edge — for the significance key.

    Mirrors :class:`AbxEnrichmentBubblePlot`'s own encoding: a significant cell is filled in its
    log2OR colour; a non-significant one is drawn hollow, with that same colour on the edge instead
    of a fill, so the colour information survives even where significance doesn't hold.

    ``x0``/``y``/``spacing`` are in the CALLER's data coordinates — both example circles must land
    within the caller's own axes limits, so pass a ``spacing`` that fits (e.g. on a 0-1 axis, not the
    literal ``1.0`` a caller might reach for as "one unit apart").
    """
    color = REFERENCE_LINE if color is None else color
    for dx, facecolor, edgecolor, label in (
        (0, color, BUBBLE_EDGE_COLOR, f"p < {ABX_ENRICHMENT_SIG_THRESHOLDS[0]}"),
        (spacing, "none", color, "n.s."),
    ):
        ax.add_patch(Circle((x0 + dx, y), radius, facecolor=facecolor, edgecolor=edgecolor,
                            linewidth=BUBBLE_EDGE_LINEWIDTH))
        ax.text(x0 + dx, y - radius * 1.6, label, ha="center", va="top", fontsize=st.FONTSIZE_SMALL)


class AbxEnrichmentMatrixPlot(_MatrixPlotBase):
    """The 38 (eos19mt class) x 15 (pathogen) Fisher-enrichment matrix, with one top track for
    each pathogen's organism class (:data:`default.ORGANISM_CLASS_ORDER`) — reusing
    :class:`plots_auroc_matrix._MatrixPlotBase`'s track/colorbar chrome, on the same page-budget
    3 cm cell grid as every other matrix figure in this repo.
    """

    def __init__(self, odds_ratio, p_value, pathogens, ax=None, cells=CELLS, name=None,
                cap=ABX_ENRICHMENT_LOG2OR_CAP):
        super().__init__(ax=ax, cells=cells)
        self.name = name or "12_abx_enrichment_matrix"
        if not odds_ratio.size:
            self._unavailable()
            return

        with np.errstate(divide="ignore", invalid="ignore"):
            log2or = np.log2(odds_ratio.astype(float))
        clipped = log2or.clip(-cap, cap)

        cmap = diverging_cmap(low="cobalt", mid="white", high="crimson").copy()
        cmap.set_bad(REFERENCE_LINE)
        norm = mcolors.Normalize(vmin=-cap, vmax=cap)

        col_labels = [_organism_label(p) for p in pathogens["pathogen"]]
        row_labels = [_class_label(c) for c in odds_ratio.index]
        heatmap(self.ax, clipped, cmap=cmap, norm=norm, annotate=False, nan_color=REFERENCE_LINE,
               x_rotation=90, row_labels=row_labels, col_labels=col_labels,
               annot_fontsize=st.FONTSIZE_SMALL, aspect="auto")

        clip_arr = clipped.to_numpy()
        or_arr = odds_ratio.to_numpy()
        p_arr = p_value.to_numpy()
        for i in range(clip_arr.shape[0]):
            for j in range(clip_arr.shape[1]):
                if np.isnan(clip_arr[i, j]):
                    continue  # degenerate class: no test, nothing to print — colour marks it NaN
                light = abs(clip_arr[i, j]) >= cap * 0.5
                text = _or_text(or_arr[i, j]) + _sig_stars(p_arr[i, j])
                self.ax.text(j, i, text, ha="center", va="center",
                            fontsize=st.FONTSIZE_SMALL, color="white" if light else INK)

        self.ax.tick_params(labelsize=st.FONTSIZE_SMALL, length=0)
        self.ax.grid(False)
        self.ax.set_axisbelow(False)
        self.ax.set_xlabel("")
        self.ax.set_ylabel("")

        divider = make_axes_locatable(self.ax)
        class_colors = dict(zip(ORGANISM_CLASS_ORDER,
                                distinct_colors(len(ORGANISM_CLASS_ORDER), levels=(0.55, 0.3))))
        self._track(divider, "top", pathogens["organism_class"].tolist(), class_colors,
                   "organism_class", vertical=False)
        self._colorbar(cmap, norm, divider, label="log2(odds ratio)", continuous=True)


def save_enrichment_figure(output_dir, odds_ratio, p_value, pathogens, name=None):
    """Build the Fisher-enrichment matrix figure and merge its footprint into
    ``figure_cells.json``.
    """
    plot = AbxEnrichmentMatrixPlot(odds_ratio, p_value, pathogens, name=name)
    footprints = {}
    if plot.is_available:
        plot.save(output_dir)
        footprints[plot.name] = list(plot.cells)
        print(f"  figure: {plot.name}")
    else:
        print(f"  [skip figure] {plot.name}: empty matrix")
    return merge_figure_cells(output_dir, footprints)


class AbxEnrichmentBubblePlot(_MatrixPlotBase):
    """A bubble-grid alternative to :class:`AbxEnrichmentMatrixPlot`, same 38 x 15 grid, but each
    cell is a circle rather than a filled cell and there is no organism-class top track:

    - **Colour** is the same clipped log2(odds ratio) diverging scale as the heatmap.
    - **Size** is proportional to the TRUE (uncapped) odds ratio, capped at
      :data:`default.ABX_ENRICHMENT_OR_SIZE_CAP` — a separate cap from the colour scale's, tuned to
      the raw odds-ratio distribution rather than its log2.
    - **Fill** marks significance: a significant cell (raw p < the first
      :data:`default.ABX_ENRICHMENT_SIG_THRESHOLDS` bound) is filled in its log2OR colour; a
      non-significant one is left UNFILLED, with that same colour drawn as the edge instead — so the
      colour information survives even where the association isn't significant.

    No in-cell text (colour + size + the companion :class:`AbxEnrichmentBubbleKeyPlot` carry the
    value instead), which is what lets this run far more compact than the heatmap.

    Bubbles are allowed to overlap where a large one sits next to a small one — rather than shrink
    every marker to rule that out, ``zorder`` is set INVERSELY to size, so a bigger bubble always
    sits behind a smaller one it overlaps rather than covering it.
    """

    #: Bubble diameter as a fraction of the (measured) cell pitch. ABOVE 1.0 (user-directed,
    #: 2026-09-03: "some of the bubbles bigger... ok with them occupying regions occupied by other
    #: bubbles") — unlike plots_euopenscreen.ActiveOverlapPiePlot's identically-named constant, which
    #: stays under 1.0 specifically to keep bubbles from ever touching, here the handful of cells at
    #: or near ABX_ENRICHMENT_OR_SIZE_CAP are meant to visibly spill into neighbouring cells; the
    #: per-cell zorder in the draw loop below is what keeps a smaller neighbour from being hidden by
    #: that overflow. Local to this class, not a shared constant, for the same reason
    #: ActiveOverlapPiePlot's is: a purely visual tuning value with no scientific meaning outside its
    #: own figure.
    MAX_DIAM_FRAC = 1.6
    #: Marker area (pt^2) for an odds ratio of exactly 0 — small but visible, so "tested, no
    #: association" reads differently from a degenerate untested cell (no bubble at all).
    MIN_AREA = 4.0

    def __init__(self, odds_ratio, p_value, pathogens, ax=None, cells=BUBBLE_CELLS, name=None,
                cap=ABX_ENRICHMENT_LOG2OR_CAP, size_cap=ABX_ENRICHMENT_OR_SIZE_CAP):
        super().__init__(ax=ax, cells=cells)
        self.name = name or "12_abx_enrichment_bubbles"
        if not odds_ratio.size:
            self._unavailable()
            return

        with np.errstate(divide="ignore", invalid="ignore"):
            log2or = np.log2(odds_ratio.astype(float))
        clipped = log2or.clip(-cap, cap)

        cmap = diverging_cmap(low="cobalt", mid="white", high="crimson")
        norm = mcolors.Normalize(vmin=-cap, vmax=cap)

        nrows, ncols = clipped.shape
        col_labels = [_organism_label(p) for p in pathogens["pathogen"]]
        row_labels = [_class_label(c) for c in odds_ratio.index]

        # Ticks/limits/aspect only — no heatmap()/imshow background, the bubbles ARE the data.
        self.ax.set_xlim(-0.5, ncols - 0.5)
        self.ax.set_ylim(nrows - 0.5, -0.5)  # row 0 on top, matching heatmap()'s imshow convention
        self.ax.set_xticks(range(ncols))
        self.ax.set_xticklabels(col_labels, rotation=90, ha="center")
        self.ax.set_yticks(range(nrows))
        self.ax.set_yticklabels(row_labels)
        self.ax.set_aspect("equal")
        self.ax.tick_params(labelsize=st.FONTSIZE_SMALL, length=0)
        self.ax.grid(False)
        self.ax.set_axisbelow(False)
        self.ax.set_xlabel("")
        self.ax.set_ylabel("")

        divider = make_axes_locatable(self.ax)
        self._colorbar(cmap, norm, divider, label="log2(odds ratio)", continuous=True)

        # Layout must be final (including the colourbar band just carved out of self.ax) before
        # bubble radii are fixed in points — same requirement as
        # plots_euopenscreen.ActiveOverlapPiePlot's own sizing. Unlike that class, no explicit
        # tight_layout() call here: the sibling AbxEnrichmentMatrixPlot uses the same divider-based
        # colourbar chrome with no tight_layout() and renders correctly, and an explicit call here
        # fought with it when this class also carried a top track (verified: it clipped the band off
        # the canvas) — a plain draw() is enough to force final axes positions.
        self.fig.canvas.draw()
        bb = self.ax.get_window_extent()
        pitch = min(bb.width * 72.0 / self.fig.dpi / ncols, bb.height * 72.0 / self.fig.dpi / nrows)
        max_area = (self.MAX_DIAM_FRAC * pitch) ** 2
        self.max_area = max_area  # exposed so the companion key panel can match it if ever needed

        clip_arr = clipped.to_numpy()
        or_arr = odds_ratio.to_numpy()
        p_arr = p_value.to_numpy()
        cells_to_draw = []
        for i in range(nrows):
            for j in range(ncols):
                if np.isnan(clip_arr[i, j]):
                    continue  # degenerate class: no test, no bubble at all
                area = _bubble_area(or_arr[i, j], min_area=self.MIN_AREA, max_area=max_area,
                                    cap=size_cap)
                significant = np.isfinite(p_arr[i, j]) and p_arr[i, j] < ABX_ENRICHMENT_SIG_THRESHOLDS[0]
                cells_to_draw.append((area, i, j, significant))

        # Largest first, ascending zorder: rank 0 (the largest bubble) gets the LOWEST zorder, so a
        # big bubble always sits BEHIND any smaller one it overlaps rather than covering it — the
        # zorder is what actually guarantees this (user-directed, 2026-09-03), the descending draw
        # order is just a belt-and-braces match to it.
        cells_to_draw.sort(key=lambda t: -t[0])
        pt_to_data = 1.0 / pitch
        for rank, (area, i, j, significant) in enumerate(cells_to_draw):
            radius = (np.sqrt(area) / 2.0) * pt_to_data  # matplotlib area convention, see
                                                          # plotting_utils.nested_size_legend
            color = cmap(norm(clip_arr[i, j]))
            # Significant: filled in the data colour. Not significant: left UNFILLED, with that same
            # colour on the edge instead — so the colour information survives either way.
            facecolor = color if significant else "none"
            edgecolor = BUBBLE_EDGE_COLOR if significant else color
            self.ax.add_patch(Circle(
                (j, i), radius=radius, facecolor=facecolor, edgecolor=edgecolor,
                linewidth=BUBBLE_EDGE_LINEWIDTH, zorder=3 + rank))


class AbxEnrichmentBubbleKeyPlot(BasePlot):
    """Standalone key for :class:`AbxEnrichmentBubblePlot`: bubble size (odds ratio) and edge style
    (significance) — same precedent as plots_euopenscreen.ActiveOverlapKeyPlot, a separate small
    panel rather than crowding the already-compact main grid.
    """

    #: Odds-ratio values shown in the size key, spanning 0 up to the size cap.
    SIZE_KEY_VALUES = (0, 25, 200, ABX_ENRICHMENT_OR_SIZE_CAP)

    def __init__(self, ax=None, cells=(1.0, 2.0), max_area=None,
                min_area=AbxEnrichmentBubblePlot.MIN_AREA, size_cap=ABX_ENRICHMENT_OR_SIZE_CAP):
        super().__init__(ax=ax, cells=cells)
        self.name = "12_abx_enrichment_bubbles_key"
        self.ax.set_xlim(0, 1)
        self.ax.set_ylim(0, 1)
        # Equal aspect: this panel's cells=(1.0, 2.0) footprint is not square, and without this the
        # manually-drawn significance circles in _sig_edge_legend (unlike nested_size_legend's own
        # per-axis point conversion) would render as ellipses.
        self.ax.set_aspect("equal")
        self.ax.set_axis_off()

        # A fixed reference area (independent of any one figure's live pitch measurement) so this
        # key renders identically regardless of what page it ends up on — same reasoning as
        # ActiveOverlapKeyPlot.WEDGE_AREA.
        max_area = 900.0 if max_area is None else max_area
        key_values = list(self.SIZE_KEY_VALUES)
        key_areas = [_bubble_area(v, min_area=min_area, max_area=max_area, cap=size_cap)
                    for v in key_values]
        nested_size_legend(self.ax, key_values, key_areas, x=0.3, y_base=0.08,
                           color=REFERENCE_LINE, title="odds ratio", label_fmt="{:.0f}")
        _sig_edge_legend(self.ax, 0.55, 0.75, radius=0.06, spacing=0.25)


def save_enrichment_bubble_figure(output_dir, odds_ratio, p_value, pathogens, name=None):
    """Build the bubble-grid enrichment figure plus its key panel, and merge both footprints into
    ``figure_cells.json`` alongside :func:`save_enrichment_figure`'s heatmap entry.
    """
    plot = AbxEnrichmentBubblePlot(odds_ratio, p_value, pathogens, name=name)
    key = AbxEnrichmentBubbleKeyPlot()
    footprints = {}
    if plot.is_available:
        plot.save(output_dir)
        key.save(output_dir)
        footprints[plot.name] = list(plot.cells)
        footprints[key.name] = list(key.cells)
        print(f"  figure: {plot.name}")
        print(f"  figure: {key.name}")
    else:
        print(f"  [skip figure] {plot.name}: empty matrix")
    return merge_figure_cells(output_dir, footprints)
