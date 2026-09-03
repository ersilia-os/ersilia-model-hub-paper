"""Concrete plot classes and the figure entry point for script 01 (model metadata).

Mirrors ``zairachem/report/plots.py`` (one class per plot type, all subclassing
:class:`plotting_base.BasePlot`). ``save_metadata_figures`` renders each panel as its own
Nature-sized figure and saves it as both PNG and vector PDF (see the function docstring).
"""

import json
import math
import os

import circlify
import numpy as np
import matplotlib
import matplotlib.patches as mpatches
import matplotlib.patheffects as patheffects
import stylia

from plotting_base import BasePlot, CELLS_PER_WIDTH, CELL_MM, pdf_page_mm
from plotting_colors import (TASK_COLORS, SUBTASK_COLORS, SOURCE_TYPE_COLORS, BAR_DEFAULT,
                             BIOAREA_GROUP_HATCH,
                             ARCH_COLORS, ARCH_DISPLAY, LICENSE_CLASS_COLORS, REFERENCE_LINE,
                             catchall_colors, hue)
from plotting_utils import (abbrev, box_with_jitter, hbar,
                            place_inside_labels, stacked_hbar)
from default import RANDOM_SEED, RUNTIME_BATCH, RUNTIME_COLUMN

# --------------------------------------------------------------------------------------
# Shared footprints (3 cm cells; see plotting_base for the grid). Defined here, above the
# classes, because they are used as default arguments — those are evaluated at import time.
# --------------------------------------------------------------------------------------
# A quarter of the 183 mm page width, square: 183/4 = 45.75 mm = 1.525 cells. Two sit side by side in
# half a page row. At this size a 10-bar panel gets a ~3.4 mm row pitch against 2.1 mm tick labels —
# 1.3 mm of clearance, so it compresses without going illegible.
QUARTER_SQUARE = (1.525, 1.525)

# Smallest square this figure uses: 45 mm (1.5 cells), for the technical box panels — four fit across
# a page row. Going below this is a legibility decision, not a sizing one: the axis furniture (y
# label, log tick labels, rotated category labels) is fixed in points and does not shrink, so the
# chart area is what gives way.
SMALL_SQUARE = (1.5, 1.5)

# A quarter of stylia's actual canvas: create_figure(width=0.25) over SIZE = 7.09 in = 180.09 mm is
# 45.02 mm = 1.5 cells. NOTE this is 0.75 mm narrower than QUARTER_SQUARE above, which is a quarter of
# the 183 mm Nature *page* rather than of stylia's canvas. The waffle and the two stacked subtask
# panels are sized on this one.
QUARTER_WIDTH = 1.5

# Height of the technical box row AND of the Biomedical Area strip, so the two sit in one band.
# 30 mm: for the box row that is 3 task rows at a ~5.9 mm pitch; for the strip, 5 groups at 3.64 mm
# (it was 4 groups at 4.55 mm before Antifungal was split out on 2026-08-07).
_BOX_ROW_HEIGHT_CELLS = 1.0

# Biomedical Area: a 25 mm strip, the narrowest panel in the repo. There is no room for a tick-label
# gutter (a fixed ~14 mm whatever the width), so the category names go INSIDE the bars and the y axis
# goes away — see hbar(inside_labels=True).
NARROW_STRIP_MM = 25.0
_STRIP_CROP_PAD_MM = 1.10      # measured; smaller than the box row's because there is no y label

# Bar thickness as a fraction of the row pitch. Lower than matplotlib's 0.8 on purpose: the strip's
# height is fixed by the box row (below), so with only a handful of groups the pitch is wide, and 0.8
# would give a 3.6 mm bar with a 0.9 mm gap — a near-solid block. At 0.70 and the four groups this was
# calibrated on, the bar was 3.2 mm (close to the 2.9 mm mark weight used across this figure, and
# enough to seat a 5 pt label) with a 1.4 mm gap. At five groups it is 2.55 mm with a 1.09 mm gap —
# still seats the label, but this fraction is the knob to revisit if a sixth group ever appears.
_STRIP_BAR_FRACTION = 0.70

# The three donut panels (licence, architecture, biomedical area) are one family: same width, same
# ring, same legend-beneath treatment, so they read as a set. 25 mm matches the Biomedical Area strip
# and puts three of them in 75 mm.
#
# The ring is the widest thing in each panel — the legends run 21-23 mm — so the crop width is the
# ring's own diameter plus a measured tight-bbox pad. NOTE calibrate that pad against a FULL script
# run, never an isolated one: stylia wipes matplotlib's font cache on import, so the first figure of a
# process draws with fallback text metrics and measures ~0.6 mm small.
DONUT_WIDTH_MM = 25.0
_DONUT_CROP_PAD_MM = 1.23
_DONUT_LEGEND_ROW_MM = 2.5     # per legend entry, measured


#: Rendered outer diameter of every donut's ring, in mm. Pinned rather than inherited, because
#: otherwise the ring is a side effect of how long the legend text happens to be: ``tight_layout``
#: shrinks the axes to fit a legend wider than it, so the two-row architecture key (19.94 mm) squeezed
#: its ring to 18.18 mm against the other two panels' 19.81 mm — an 8 % difference driven by nothing
#: but the word "AMD64+ARM64". Must stay at or below the narrowest natural axes width, or the ring
#: overflows its axes and is clipped. See :meth:`DonutPlot.pin_ring`.
DONUT_RING_MM = 19.6


def _donut_cells(n_legend_rows):
    """``(rows, cols)`` for a donut panel with ``n_legend_rows`` entries beneath the ring.

    Height tracks the legend rather than being fixed: the ring is square and width-bound, so the panel
    is ring + legend band. **A two-entry donut therefore comes out shorter than a four-entry one** —
    the architecture panel is ~5 mm shorter than the other two, so align the set on the rings, not on
    the panel boxes. Declaring roughly the right height matters because the legend is anchored below
    the axes: over-declare and ``tight_layout`` opens a gap between ring and key.
    """
    cols = (DONUT_WIDTH_MM - _DONUT_CROP_PAD_MM) / CELL_MM
    rows = (DONUT_WIDTH_MM + n_legend_rows * _DONUT_LEGEND_ROW_MM) / CELL_MM
    return (rows, cols)


def _narrow_strip_cells(n_bars, width_mm=NARROW_STRIP_MM):
    """``(rows, cols)`` for the Biomedical Area strip: ``width_mm`` wide once cropped, and **the same
    height as the technical box row**, so the two can sit side by side in one band.

    Height is therefore NOT derived from ``n_bars`` — the bar count only decides how the fixed axes
    height is divided. At five groups that is an 18.19 mm axes over 5 rows, a 3.64 mm pitch.
    """
    cols = (width_mm - _STRIP_CROP_PAD_MM) / CELL_MM
    return (_BOX_ROW_HEIGHT_CELLS, cols)


class FieldBarPlot(BasePlot):
    """Horizontal bar chart of a metadata field's value counts.

    Parameters
    ----------
    counts : DataFrame with columns ``value`` and ``count`` (already sorted).
    title  : panel title.
    colors : optional list of per-bar colours (same order as ``counts``); defaults to
             a single ``BAR_DEFAULT`` colour for every bar.
    color_fn : optional ``values -> list of colours`` callable, applied AFTER the top-N cap. Use
             this rather than ``colors`` when the colour depends on the value (e.g.
             :func:`plotting_colors.catchall_colors`), so the list can never fall out of step with
             the rows actually drawn.
    n      : optional top-N cap; when set, only the first ``n`` rows are shown and the
             title gets a "(top n)" suffix.
    cells  : footprint on the reference grid as ``(rows, cols)`` — taller for panels with
             more bars (see ``save_metadata_figures``).
    inside_labels : draw the category names inside the bars instead of as y tick labels, and drop the
             y axis. For panels too narrow to afford a tick-label gutter — see :func:`hbar`, which
             requires a pale bar colour for it to read. Sets ``self.label_texts`` for measurement.
    xlabel : count-axis label. Shortened on narrow panels: at 25 mm "Number of models" is 17.6 mm
             against an 18.7 mm axes, close enough to overhanging that it would set the crop width.
    """

    def __init__(self, ax=None, counts=None, title="", colors=None, color_fn=None,
                 n=None, cells=(3, 3), inside_labels=False,
                 xlabel="Number of models"):
        BasePlot.__init__(self, ax=ax, cells=cells)
        self.name = title.lower().replace(" ", "_")
        if n:
            counts = counts.head(n)
            title = f"{title} (top {n})"
        if color_fn is not None:
            colors = color_fn(counts["value"])
        self.label_values = counts["count"].tolist()
        self.label_texts = hbar(
            self.ax, counts["value"].tolist(), self.label_values,
            colors=colors if colors is not None else BAR_DEFAULT,
            inside_labels=inside_labels,
            bar_fraction=_STRIP_BAR_FRACTION if inside_labels else 0.8)
        self.label(xlabel=xlabel, ylabel="" if inside_labels else " ", title=title)

    def measure_labels(self):
        """Place the inside-bar labels, then measure them. Call BEFORE ``save`` closes the figure.

        Placement runs here rather than at draw time because it weighs a label width fixed in points
        against an axes width ``tight_layout`` is still free to change, so layout must go first. An
        inside label is also not clipped, so one hanging off a spine would spill outside the plot area
        and *set* the crop width — the failure mode a fixed-width panel cannot absorb. Sets
        ``self.n_inside`` / ``self.n_after`` and returns
        ``(axes_mm, widest_label_mm, widest_label_text)``.
        """
        import matplotlib.pyplot as plt

        plt.figure(self.fig.number)
        plt.tight_layout()
        self.n_inside, self.n_after = place_inside_labels(
            self.ax, self.label_texts, self.label_values)
        self.fig.canvas.draw()
        r = self.fig.canvas.get_renderer()
        self.axes_mm = self.ax.get_window_extent().width / self.fig.dpi * 25.4
        widest = max(self.label_texts,
                     key=lambda t: t.get_window_extent(r).width, default=None)
        self.label_mm = (widest.get_window_extent(r).width / self.fig.dpi * 25.4
                         if widest is not None else 0.0)
        self.widest_label = widest.get_text() if widest is not None else ""
        return self.axes_mm, self.label_mm, self.widest_label


#: Layout dpi for a standalone StackedFieldBarPlot, raised from matplotlib's 100.
#:
#: Nothing to do with output resolution — the PDF is vector and stylia sets the PNG's dpi on save.
#: This is about the LAYOUT GRID: matplotlib's canvas is a whole number of pixels, so every size the
#: layout can express is a multiple of one pixel — **0.254 mm at dpi 100**. Two panels sized to draw
#: bars of the same thickness (see :func:`_stack_axes_heights`) need a specific, non-round split of
#: the height between them, and on a 0.254 mm grid the nearest achievable split leaves the two bars
#: **1.8% apart**, which no choice of constants can improve. At 600 the grid is 0.042 mm and the same
#: solve lands inside 0.3%.
#:
#: Set on this class only, so the rest of the figure's carefully calibrated crops are untouched — a
#: global dpi change would shift every measured panel size in the repo by up to a quarter of a mm.
_LAYOUT_DPI = 600


class StackedFieldBarPlot(BasePlot):
    """A metadata field's value counts as horizontal bars, SEGMENTED by Subtask.

    A bar's total length is exactly the field's own count, while the segments say which subtasks make
    up that total — so Source Type and Output stop being independent of the task breakdown and the
    reader gets the joint distribution from one panel instead of cross-referencing two.

    Segment colour is the subtask palette, NOT a palette of its own: the field itself is encoded by
    the bar's position/label, so spending a second colour dimension on it would be redundant.

    Parameters
    ----------
    table  : wide DataFrame of model counts. Index = field values, ordered top-to-bottom as
             drawn; columns = segment labels in stacking order (left to right).
    colors : dict segment label -> colour; must cover every column of ``table``.
    name   : output file stem (e.g. ``"source_type_by_subtask"``).
    legend_kw : kwargs forwarded to the swatch legend, or ``None`` for **no legend** — which is what
             both panels here pass, since ``task_subtask`` and the waffle are the shared subtask key
             and a 6-entry legend does not fit on a ~52 x 25 mm panel. ``{}`` gives the primitive's
             defaults.
    xlim   : count-axis limits, or ``None`` to autoscale (the normal case here — the two panels are
             NOT on a shared scale, so each one's axis runs to its own longest bar).
    show_xlabel : ``False`` drops the ``"Number of models"`` axis title while keeping the tick labels,
             so a stacked pair states the quantity once instead of twice. The panel keeps a fully
             readable count axis of its own, unlike hiding the ticks as well: it just borrows the
             wording from its partner. Worth ~3.3 mm of the fixed tick+label band.
    cells  : footprint on the reference grid as ``(rows, cols)``.
    """

    def __init__(self, ax=None, table=None, colors=None, name="stacked_field",
                 legend_kw=None, xlim=None, cells=(2, 2), show_xlabel=True):
        # The layout dpi has to be in force when the figure is CREATED — it is what sets the pixel
        # grid every later size is quantised onto. Setting it afterwards does not re-quantise.
        with matplotlib.rc_context({"figure.dpi": _LAYOUT_DPI}):
            BasePlot.__init__(self, ax=ax, cells=cells)
        self.name = name
        segments = list(table.columns)
        # stacked_hbar wants one {segment: value} dict per row. Its xlim default of (0, 1) is for
        # fraction stacks, so pass the count axis through explicitly (None = autoscale).
        rows = [row.to_dict() for _, row in table.iterrows()]
        stacked_hbar(self.ax, list(table.index), rows, segments, colors, xlim=xlim)
        if legend_kw is not None:
            self.legend({s: colors[s] for s in segments}, **legend_kw)
        # Tick labels always stay — each panel autoscales to its own longest bar, so it needs its own
        # readable numbers. Only the axis TITLE is optional, so a stacked pair says what is being
        # counted once rather than twice.
        self.label(xlabel="Number of models" if show_xlabel else "", ylabel=" ", title=name)
        self.n_bars = len(table)
        self.show_xlabel = show_xlabel
        self.bar_mm = self.pitch_mm = self.axes_h_mm = self.fig_h_mm = None

    def measure_geometry(self):
        """Drawn bar thickness, row pitch and axes height in mm. Call BEFORE ``save`` closes the figure.

        Measured off the *rendered* axes rather than computed from the sizing constants, because that
        is the only check that catches those constants drifting: a stacked pair sized for equal bar
        thickness (see :func:`_stack_cells`) is only equal if the fixed furniture bands the solver
        subtracts are still what the renderer actually spends. Pitch is one data unit through
        ``transData`` and thickness is a segment patch's own height. ``fig_h_mm - axes_h_mm`` is what
        re-calibrates ``_AXES_BAND_MM``, and it is measured against the FIGURE, not the footprint: the
        figure is the footprint quantised onto the canvas pixel grid, so the footprint is the one size
        in the chain that no band is constant against.

        Sets ``self.pitch_mm`` / ``self.bar_mm`` / ``self.axes_h_mm`` / ``self.fig_h_mm``. The panel's
        page size is NOT measured here — that is ``plotting_base.pdf_page_mm`` on the saved file, since
        the pre-save estimate is ~1.2 mm short on width.
        """
        import matplotlib.pyplot as plt

        plt.figure(self.fig.number)
        plt.tight_layout()
        self.fig.canvas.draw()
        px_mm = 25.4 / self.fig.dpi
        t = self.ax.transData
        self.pitch_mm = abs(t.transform((0, 1))[1] - t.transform((0, 0))[1]) * px_mm
        self.bar_mm = (self.ax.patches[0].get_window_extent().height * px_mm
                       if self.ax.patches else 0.0)
        self.axes_h_mm = self.ax.get_window_extent().height * px_mm
        self.fig_h_mm = self.fig.get_size_inches()[1] * 25.4
        return self.pitch_mm, self.bar_mm, self.axes_h_mm


class TaskSubtaskBarPlot(BasePlot):
    """Combined Task + Subtask panel: one horizontal bar per subtask, in that subtask's own
    shade of its parent task's hue — the SAME ``SUBTASK_COLORS`` the stacked ``*_by_subtask``
    panels use. Bars are grouped by task and ordered by count within it, so the task grouping
    reads off the hue runs while each bar is named on the axis.

    Deliberately carries **no legend**: every bar is already labelled with its subtask name next
    to its own colour, so this panel *is* the subtask colour key, and placing it beside the
    ``*_by_subtask`` panels lets those drop their legends too (they have no room for one at
    ~52 x 25 mm). Per-task totals live in ``task_counts.csv``. ``task_subtask_waffle`` is the other
    panel that can serve as the key — it carries its own legend, on abbreviated labels, so the full
    subtask names reach the reader through THIS panel's tick labels.

    Parameters
    ----------
    sub   : subtask counts DataFrame with columns ``value`` (short display label), ``count`` and
            ``parent`` (already ordered by parent task, then count within task).
    cells : footprint on the reference grid as ``(rows, cols)``.
    """

    def __init__(self, ax=None, sub=None, cells=(2, 2)):
        BasePlot.__init__(self, ax=ax, cells=cells)
        self.name = "task_subtask"
        colors = [SUBTASK_COLORS[v] for v in sub["value"]]
        hbar(self.ax, sub["value"].tolist(), sub["count"].tolist(), colors=colors)
        self.label(xlabel="Number of models", ylabel=" ", title="Tasks & subtasks")


class _HorizontalTaskPanel(BasePlot):
    """Shared scaffolding for a panel in the horizontal task row (see ``_box_row_cells``).

    Owns the two things every panel in that row must do identically, so a panel cannot drift out of
    register with its neighbours: the **task axis** (tasks down the y axis, first on top, tick labels
    only on the leftmost panel) and the **label-fit measurement** the entry point reports.

    Subclasses draw their marks against integer y positions ``0..len(tasks)-1`` and then call
    :meth:`_task_axis`.
    """

    def _task_axis(self, tasks, show_y, xlabel, log=True):
        ax = self.ax
        if log:
            ax.set_xscale("log")
        ax.set_ylim(len(tasks) - 0.5, -0.5)      # inverted: first task on top
        ax.set_yticks(range(len(tasks)))
        ax.set_yticklabels(tasks if show_y else [""] * len(tasks))
        self.label(xlabel=xlabel, ylabel="", title=self.name)

    def measure_fit(self):
        """Measure the drawn axes and xlabel widths in mm; call BEFORE ``save`` closes the figure.

        The metric label is centred under the axes, so one wider than the axes overhangs it and
        *sets* the tight-crop width — silently blowing the row's width budget, which is split on the
        assumption that the axes is the widest thing in each panel. Sets ``self.axes_mm`` and
        ``self.xlabel_mm`` so the entry point can report the margin and the check stays visible on
        every run instead of being rediscovered from a broken layout.
        """
        import matplotlib.pyplot as plt

        plt.figure(self.fig.number)
        plt.tight_layout()
        self.fig.canvas.draw()
        r = self.fig.canvas.get_renderer()
        self.axes_mm = self.ax.get_window_extent().width / self.fig.dpi * 25.4
        self.xlabel_mm = self.ax.xaxis.get_label().get_window_extent(r).width / self.fig.dpi * 25.4
        return self.axes_mm, self.xlabel_mm


class TaskMetricBoxPlot(_HorizontalTaskPanel):
    """One HORIZONTAL box-with-dots per Task for a per-model technical metric.

    Tasks run down the y axis (first on top, as everywhere else in the repo) and the metric runs
    along a log x axis. Three of these sit side by side as one row — runtime, image size and output
    dimension — sharing the task axis, so only the leftmost panel draws its tick labels
    (``show_y``). See ``_box_row_cells`` for how the row's widths are split.

    Values ``<= 0`` are treated as NOT MEASURED, because Airtable stores a ``-1`` sentinel for a
    benchmark that was never run. They are skipped rather than imputed. ``self.coverage`` records
    ``{task: (n_measured, n_total)}``, which the caller prints and writes to
    ``technical_metrics_summary.csv``.

    **The per-task n is NOT in the tick labels**, unlike the earlier vertical version of this panel.
    A shared axis can carry only one label set, and coverage differs per metric (runtime is
    131/59/13, image size and output dimension are both 133/60/25), so a single labelled ``n`` would
    be wrong for two panels out of three. It lives in the summary CSV and the run log instead — **a
    caption must state that the runtime box for Sampling rests on 13 of 25 models.**

    A Task with **no** measurements still gets its slot, marked "not measured" in the neutral colour
    instead of being dropped: an absent box is information, an absent category is a misreading.

    The metric axis is logarithmic by default — all three metrics span more than a decade (runtime
    16-1626 s, image size 290-10242 MB, output dimension 1-5000) and a linear axis flattens the bulk
    of each distribution against the floor.

    Parameters
    ----------
    df       : the (Status == "Ready") metadata DataFrame.
    column   : the metric column to draw.
    name     : output file stem.
    xlabel   : metric-axis label, including units.
    log      : logarithmic metric axis (default True).
    show_y   : ``False`` hides the task tick labels (keeping the tick marks) for the second and third
               panels of the row, which read off the first. Such a panel **cannot be placed on its
               own** — it has no category axis.
    cells    : footprint on the reference grid as ``(rows, cols)``.
    """

    #: Box body and swarm geometry, in category-axis data units (1.0 = one task row). The swarm band
    #: is as wide as the box on purpose: these panels are ~6 mm per row, and at the old 0.10 jitter
    #: 133 Annotation points piled into a ~1 mm line. The box is unfilled (house style), so the
    #: swarm is what carries the distribution and the box only marks the quartiles.
    _BOX_WIDTH = 0.5
    _JITTER = 0.20

    def __init__(self, ax=None, df=None, column=None, name="task_metric", xlabel="",
                 log=True, show_y=True, cells=SMALL_SQUARE):
        BasePlot.__init__(self, ax=ax, cells=cells)
        self.name = name
        ax = self.ax
        rng = np.random.default_rng(RANDOM_SEED)

        tasks = [t for t in TASK_COLORS if (df["Task"] == t).any()]
        self.coverage = {}
        for i, task in enumerate(tasks):
            sub = df.loc[df["Task"] == task, column]
            vals = sub[sub > 0].to_numpy(dtype=float)
            self.coverage[task] = (len(vals), len(sub))
            if len(vals):
                box_with_jitter(ax, vals, i, TASK_COLORS[task], vert=False,
                                width=self._BOX_WIDTH, jitter_width=self._JITTER,
                                rng=rng, point_size=3, point_alpha=0.55)
            else:
                # Keep the slot. Placed mid-width in axes coordinates (get_yaxis_transform: axes x,
                # data y) so it sits where the box would have been, needing no x value on a log scale.
                ax.text(0.5, i, "not measured", ha="center", va="center",
                        transform=ax.get_yaxis_transform(), fontsize=stylia.FONTSIZE_SMALL,
                        color=REFERENCE_LINE)

        self._task_axis(tasks, show_y, xlabel, log=log)


class TaskOutputDimensionCirclesPlot(_HorizontalTaskPanel):
    """Output Dimension per Task as **decade-binned, area-proportional circles**.

    A drop-in replacement for the box-and-swarm version of this panel: same task axis, same
    horizontal shape, same footprint. The swarm was the wrong mark for this column — Output Dimension
    is heavily **tied** (68 of 133 Annotation models output a single value and 102 of them fall in the
    1-9 bin; 100 and 1000 recur across Representation and Sampling), so jittered points piled onto a
    handful of x positions and the visual density said more about the jitter than about the data.

    Binning by decade replaces overplotting with an explicit count: one circle per (task, decade),
    **area proportional to the number of models**, drawn at the bin's *geometric centre* so it sits
    between the two decade ticks that bound it rather than on top of one of them — a circle on the
    10² tick would read as "exactly 100", which is a real value in this column.

    *Area, not diameter, carries the count* — the perceptually correct choice for a count, and the
    same convention as every other sized mark in the repo. The 10 non-empty bins span 2 to 102
    models, a 51x range, which is a 7.1x range in diameter; ``_MAX_DIAMETER_MM`` pins the largest
    circle just inside the row pitch, leaving the smallest at ~0.7 mm. There is **no in-panel size
    key** — at 34 x 31 mm there is nowhere to put one — so exact counts come from
    ``output_dimension_bins.csv``; read the panel for the pattern (Annotation outputs one value,
    Representation spreads across three decades, Sampling sits at 10²-10³).

    ``self.bins`` holds the ``{task: {decade: count}}`` cross-tab for the caller to write out.

    Parameters
    ----------
    df     : the (Status == "Ready") metadata DataFrame.
    column : the metric column to bin (``"Output Dimension"``).
    name   : output file stem.
    xlabel : metric-axis label.
    show_y : ``False`` hides the task tick labels; see :class:`_HorizontalTaskPanel`.
    cells  : footprint on the reference grid as ``(rows, cols)``.
    """

    #: Diameter of the largest circle. The row pitch at this footprint is ~5.9 mm, so 4.6 mm keeps
    #: the biggest bin inside its own row with a visible gap to its neighbours.
    _MAX_DIAMETER_MM = 4.6

    def __init__(self, ax=None, df=None, column="Output Dimension", name="output_dimension",
                 xlabel="Output dimension", show_y=True, cells=SMALL_SQUARE):
        BasePlot.__init__(self, ax=ax, cells=cells)
        self.name = name
        ax = self.ax

        tasks = [t for t in TASK_COLORS if (df["Task"] == t).any()]
        vals = df[column]
        # Values <= 0 would have no decade; this column has none (all 218 models are >= 1), but the
        # guard keeps the panel honest if a -1 sentinel ever appears, as in the runtime column.
        ok = df[vals > 0]
        self.n_skipped = len(df) - len(ok)
        decades = np.floor(np.log10(ok[column].to_numpy(dtype=float))).astype(int)
        lo, hi = int(decades.min()), int(decades.max())

        self.bins = {}
        counts, xs, ys, colors = [], [], [], []
        for i, task in enumerate(tasks):
            m = (ok["Task"] == task).to_numpy()
            per = {k: int((decades[m] == k).sum()) for k in range(lo, hi + 1)}
            self.bins[task] = per
            for k, c in per.items():
                if not c:
                    continue
                counts.append(c)
                xs.append(10.0 ** (k + 0.5))     # geometric centre of [10^k, 10^(k+1))
                ys.append(i)
                colors.append(TASK_COLORS[task])

        # matplotlib's scatter `s` is an area in points squared, so s proportional to count IS area
        # proportional to count; the largest bin is pinned to _MAX_DIAMETER_MM and the rest follow.
        d_max_pt = self._MAX_DIAMETER_MM / 25.4 * 72.0
        s_max = d_max_pt ** 2
        cmax = max(counts)
        ax.scatter(xs, ys, s=[s_max * c / cmax for c in counts], color=colors,
                   edgecolors="none", zorder=3)

        # Task axis FIRST: it applies the log scale, which installs a LogLocator and would discard
        # any ticks set before it. Then ticks on the decade EDGES, so the circles visibly sit inside
        # the intervals they summarise rather than on the boundaries.
        self._task_axis(tasks, show_y, xlabel, log=True)
        ax.set_xlim(10.0 ** (lo - 0.15), 10.0 ** (hi + 1.15))
        ax.set_xticks([10.0 ** k for k in range(lo, hi + 2)])


def _donut(ax, labels, values, colors, *, ring=0.42, center_value=None, hatches=None):
    """House-style donut: ring clockwise from 12 o'clock, **no labels on the wedges**.

    Labels never go on the wedges because there is nowhere to put them. At the width this is used
    for (30 mm, one cell) labelling around the ring is not merely tight but geometrically impossible:
    ``"Non-commercial 4 (2%)"`` is 20.07 mm of text and even the name alone is 13.46 mm, so labels on
    two sides would consume ~27 mm of the 30 mm before the circle got any. The caller therefore supplies
    a legend beneath, which is also what rescues the 4-of-218 wedge — at 6.6 degrees it can carry no
    label of its own however much room the panel has.

    The hole is not decoration: it holds the **total**, which a pie has nowhere to put, so the panel
    states its own n instead of leaving it to the caption. Returns the wedge list.
    """
    values = np.asarray(values, dtype=float)
    wedges, _ = ax.pie(
        values, radius=1.0, startangle=90, counterclock=False, colors=colors,
        wedgeprops=dict(width=ring, edgecolor="white", linewidth=0.6))
    if hatches:
        # matplotlib draws a hatch in the patch's EDGE colour, which is already white here, so a
        # patterned wedge reads as white-on-hue and still registers as its colour at a glance.
        for w, h in zip(wedges, hatches):
            if h:
                w.set_hatch(h)
    ax.set_aspect("equal")
    # ``ax.pie`` always sets the limits to +/-1.25 to leave room for the outside labels it normally
    # draws. With no such labels that is 20 % of dead margin on every side, which both shrinks the ring
    # and opens a ~3.5 mm gap above a legend anchored below the axes. Pull them in to hug the circle.
    ax.set_xlim(-1.02, 1.02)
    ax.set_ylim(-1.02, 1.02)
    if center_value is not None:
        # The bare total, centred. No unit word beneath it: at this ring size a second line crowds the
        # hole, and every legend row underneath already says what is being counted.
        ax.text(0, 0, f"{int(center_value):,}", ha="center", va="center",
                fontsize=stylia.FONTSIZE_SMALL + 2)
    return wedges


class DonutPlot(BasePlot):
    """A composition as a ring with its key beneath — the shared base for this figure's three donuts.

    One implementation for `license_class_donut`, `docker_architecture` and `biomedical_area_donut`,
    so the three read as a set: same 25 mm width, same ring thickness, same legend treatment, same
    total in the hole. A subclass only prepares ``(labels, values, colors)``.

    **Why a legend rather than labels around the ring.** At 25 mm, labelling around the circle is
    geometrically impossible rather than merely tight: ``"Non-commercial 1 (<1%)"`` is 20.07 mm of text
    and the name alone is 13.46 mm, so labels on two sides would take more than the whole panel before
    the ring got any. The legend sits beneath at one entry per row, carrying name, count and share, so
    every number is in type rather than estimated off a wedge — which is also the only way a 1.7 degree
    wedge gets named at all.

    **The hole carries the total**, the one thing a pie has nowhere to put, so each panel states its
    own n instead of delegating it to the caption.

    Parameters
    ----------
    labels  : category names, in draw order (clockwise from 12 o'clock).
    values  : counts, same order. The legend shows count and share; the hole shows the sum.
    colors  : one colour per wedge.
    name    : output file stem.
    hatches : optional list of matplotlib hatch strings, one per wedge, for a panel that distinguishes
              categories by **pattern** instead of by hue (see :class:`BiomedicalAreaDonutPlot`). The
              legend swatches carry the same patterns.
    """

    def __init__(self, labels, values, colors, *, name, hatches=None, total=None,
                 ax=None, cells=None):
        values = np.asarray(values, dtype=float)
        BasePlot.__init__(self, ax=ax, cells=cells or _donut_cells(len(labels)))
        self.name = name
        self.counts = {l: int(v) for l, v in zip(labels, values)}
        # ``total`` overrides the wedge sum, and it matters wherever a model can belong to more than
        # one wedge: the Biomedical Area groups sum to 98 across 93 models, and a hole reading 98 would
        # state a model count that does not exist. Pass the distinct count; the discrepancy against the
        # legend rows is real and belongs in the caption, not hidden by quietly showing the sum.
        self.total = int(values.sum() if total is None else total)
        _donut(self.ax, labels, values, colors, center_value=self.total, hatches=hatches)
        # Legend rows are in wedge order, so scanning the key top-to-bottom walks the ring clockwise
        # from 12 o'clock, and they carry the count — the number the ring cannot show.
        #
        # NO percentage here, and that is what forces the panel size: the legend must
        # stay narrower than the ring or tight_layout shrinks the axes to fit it and the ring collapses,
        # and at 25 mm "AMD64 + ARM64 129 (62%)" is 25.94 mm against a ~20 mm budget. Dropping it costs
        # little — the share is exactly what the ring already encodes, and the hole gives the total to
        # divide by, so the panel would have been saying the same thing three ways.
        rows = [f"{l} {int(v)}" for l, v in zip(labels, values)]
        # Handle spacing is tighter than the repo default: the legend has to stay narrower than the
        # ring (see pin_ring), and at 25 mm the architecture key had no room to spare.
        self.legend({r: c for r, c in zip(rows, colors)},
                    hatches=({r: h for r, h in zip(rows, hatches)} if hatches else None),
                    ncol=1, loc="upper center", bbox_to_anchor=(0.5, -0.02),
                    handlelength=0.9, handletextpad=0.4, labelspacing=0.35,
                    borderpad=0.1)
        # A pie axis has no x or y, but stylia's preset leaves "X-axis / Units" placeholders on any
        # axis whose labels are never set — so this call is what clears them, not decoration.
        self.label(xlabel="", ylabel="")

    def pin_ring(self):
        """Force the ring to render at exactly ``DONUT_RING_MM`` across. Call BEFORE ``save``.

        Without this the ring's size is decided by ``tight_layout``, which shrinks the axes to fit a
        legend wider than it — so a panel's ring shrank or grew according to the length of its longest
        label, and the three donuts came out visibly different sizes. Pinning makes the family strictly
        uniform, and because a pie axis draws no frame it also makes the panels' crop widths uniform:
        the crop is then ring + savefig pad in every case.

        The ring is drawn at radius 1 and the *limits* are scaled instead of the radius, so the wedge
        geometry, the hatching and the hole all scale together. Runs after ``tight_layout`` for the
        same reason the bar-label placement does: it needs the axes' final size in millimetres.
        Sets ``self.axes_mm`` (the natural axes width, i.e. the headroom this panel had).
        """
        import matplotlib.pyplot as plt

        plt.figure(self.fig.number)
        plt.tight_layout()
        self.fig.canvas.draw()
        self.axes_mm = self.ax.get_window_extent().width / self.fig.dpi * 25.4
        # A radius-1 circle spans 2 data units; over an axes `axes_mm` wide with limits +/-lim it
        # renders at axes_mm / lim across, so lim = axes_mm / target.
        lim = self.axes_mm / DONUT_RING_MM
        self.ax.set_xlim(-lim, lim)
        self.ax.set_ylim(-lim, lim)
        return self.axes_mm


class LicenseClassDonutPlot(DonutPlot):
    """Licence composition as a **donut** — the compositional counterpart to the ``license`` bar panel.

    Slices are the four **reuse classes**, not the twelve individual licences. That is a deliberate
    limit, not a shortcut: five of the twelve licences cover exactly one model each, which on 218
    models is a **1.7 degree wedge** — invisible, unlabellable, and impossible to tell apart from the
    other four. The per-licence detail is what the bar panel is for; a circle can only carry the
    four-way split. Even here Non-commercial is 4 models (6.6 degrees), a sliver that exists on the
    ring only to be accounted for — it is the legend that actually tells a reader it is there.

    Ordered by count, so the smallest wedge finishes at the 12 o'clock start line rather than being
    buried between two large ones, where it would read as a rendering artefact.

    Parameters
    ----------
    counts : the grouped licence DataFrame (``value``, ``count``, ``class``) from the script.
    """

    def __init__(self, ax=None, counts=None, cells=None):
        by_class = counts.groupby("class")["count"].sum()
        # The reindex below silently drops any class with no LICENSE_CLASS_COLORS entry, which would
        # take its models out of the donut without touching the total in the hole. Counterpart to the
        # unmapped-licence guard in 01_ersilia_metadata.py: a new reuse class needs a colour here.
        _uncoloured = sorted(set(by_class.index) - set(LICENSE_CLASS_COLORS))
        if _uncoloured:
            raise KeyError(f"Licence reuse classes with no LICENSE_CLASS_COLORS entry: {_uncoloured}")
        by_class = (by_class.reindex(list(LICENSE_CLASS_COLORS))
                    .dropna()
                    .sort_values(ascending=False))
        DonutPlot.__init__(self, list(by_class.index), by_class.to_numpy(dtype=float),
                           [LICENSE_CLASS_COLORS[k] for k in by_class.index],
                           name="license_class_donut", ax=ax, cells=cells)


class ArchitectureDonutPlot(DonutPlot):
    """Share of models whose container is also built for ARM64, as a donut.

    The underlying field only ever holds ``AMD64`` or ``AMD64,ARM64`` — there is **no ARM-only build**
    — so this is "also built for ARM" versus "x86 only", not a three-way split.

    This is a **snapshot, not a trend**: dual-arch is 45 % among models incorporated in 2021 and 77 %
    among 2026 ones, so the headline share is accumulated stock rather than current build practice, and
    the incorporation-date range of the source metadata belongs in any caption.

    Colours are cobalt / tangerine rather than the turquoise + periwinkle this panel used as a pie:
    those two now belong to the licence donut beside it, and no hue may mean two different things
    across a set read together. See ``ARCH_COLORS`` — tangerine is a plain categorical hue here, not a
    warning about x86-only builds.

    Parameters
    ----------
    df : the (Status == "Ready") metadata DataFrame.
    """

    def __init__(self, ax=None, df=None, cells=None):
        counts = (df["Docker Architecture"].map(ARCH_DISPLAY)
                  .value_counts()
                  .reindex(list(ARCH_COLORS))
                  .dropna())
        DonutPlot.__init__(self, list(counts.index), counts.to_numpy(dtype=float),
                           [ARCH_COLORS[k] for k in counts.index],
                           name="docker_architecture", ax=ax, cells=cells)


class BiomedicalAreaDonutPlot(DonutPlot):
    """Biomedical Area groups as a donut — **one hue, five fill patterns**.

    The alternative to the `biomedical_area` bar strip on the same five groups and the same 95 Activity
    prediction models; both are rendered and one is picked at layout time.

    Every wedge is the Annotation crimson, because every model in this panel *is* an Annotation model —
    colour would be encoding nothing. The categories are separated by **pattern** instead
    (``BIOAREA_GROUP_HATCH``): solid for the largest group, then progressively lighter-inked patterns,
    so the ink ordering matches the size ordering, with the cross-hatch on the catch-all where it reads
    as "mixed". Patterns are white over the crimson, so a wedge still registers as red at a glance.
    ``Antifungal`` is the one exception to the ink ordering — it mirrors ADMET's stroke to read as a
    pair, since it was split out of ``Antimicrobial`` on 2026-08-07.

    That makes this the one panel in the set whose key is **not** redundant with a colour a reader
    could name: the legend swatches carry the patterns, so it cannot describe a mark the ring does not
    draw.

    Parameters
    ----------
    counts : the grouped Biomedical Area DataFrame (``value``, ``count``) from the script.
    """

    def __init__(self, ax=None, counts=None, cells=None):
        labels = counts["value"].tolist()
        DonutPlot.__init__(self, labels, counts["count"].to_numpy(dtype=float),
                           [TASK_COLORS["Annotation"]] * len(labels),
                           hatches=[BIOAREA_GROUP_HATCH.get(l, "") for l in labels],
                           # Distinct models, set by the script — NOT the wedge sum, which is 98
                           # because three models carry areas in two different groups.
                           total=counts.attrs.get("n_models"),
                           name="biomedical_area_donut", ax=ax, cells=cells)


class TaskSubtaskWafflePlot(BasePlot):
    """One square per model, coloured by subtask — the *unit* alternative to
    :class:`TaskSubtaskBarPlot`, kept alongside it so the two can be compared.

    Squares are laid out in reading order through the subtasks (which are already grouped by parent
    task), so each subtask is a contiguous run and each task a contiguous band of one hue. What this
    buys over a bar chart is that the hub total is *shown* rather than stated: n is countable, every
    model is one mark, and the six subtasks are compared as areas of a whole. What it costs is
    precision — reading 52 vs 39 off a waffle means counting.

    Unlike the other subtask panels this one carries its **own legend**, so it is self-contained and
    does not have to travel next to ``task_subtask``. Legend labels carry each subtask's count, which
    is exactly what the waffle itself conveys badly.

    The legend is also what makes the panel square. The 16-wide grid alone crops landscape (1.23 at
    60 mm, measured when 208 models filled exactly 13 rows); the legend band adds height below it and
    squares the panel up. At 218 models the grid runs to a 14th row, so it is marginally taller than
    those measurements and the legend has correspondingly less work to do. Going the other way is
    counter-productive: a portrait grid is narrower than the legend, which then sets the width while
    the extra rows add height (measured at 60 mm: 13 cols 0.78, 14 cols 0.83, 15 cols 0.94, **16 cols
    1.05**, 18 cols 1.22).

    **At quarter width the legend is what governs, so its labels are abbreviated**
    (``_LEGEND_ABBREV``). With the full subtask names the two-column key measures 48.3 mm against a
    45.0 mm panel, so the *legend* sets the crop: the grid gets squeezed to 35.9 mm inside a 47.7 mm
    page with dead space either side, and the squares fall to 1.99 mm. Abbreviated, the key comes to
    38.0 mm — narrower than the grid — so the grid sets the width instead and expands to 41.2 mm
    (**2.57 mm squares**), and the panel crops to 46.2 x 46.1 mm, square to within 0.2 %. The counts
    are kept: they are the whole point of this legend, and they are not what overflowed.

    ``self.blank`` records how many trailing cells of the grid went unfilled, so a caller can flag a
    ragged last row. **A ragged row is now unavoidable**: 16 columns divided 208 exactly (13 full
    rows), but 218 = 2 x 109 has no divisor in the usable 13-18 range, so every column count leaves a
    remainder. 16 is kept for the aspect ratio measured above and leaves 6 blank cells in row 14.

    Parameters
    ----------
    sub   : subtask counts DataFrame with columns ``value`` (short display label) and ``count``,
            already ordered by parent task then count within task.
    cols  : squares per row. The default 16 leaves ``self.blank`` = 6 cells empty in the last row at
            the current 218 models; it divided the earlier 208 into exactly 13 full rows.
    cells : footprint on the reference grid as ``(rows, cols)``.
    """

    #: Legend-only abbreviations of the subtask display labels, applied so the two-column key fits
    #: the quarter-width panel (see the class docstring for the measurements). Deliberately NOT
    #: folded into ``SUBTASK_DISPLAY``: the full names still reach the ``task_subtask`` tick labels
    #: and the stacked panels, and this key sits next to them.
    _LEGEND_ABBREV = {
        "Activity prediction": "Activity pred.",
        "Property prediction": "Property pred.",
        "Similarity search": "Similarity",
    }

    _COLS = 16
    _PAD = 0.16       # white gap between squares, in square-widths

    def __init__(self, ax=None, sub=None, cols=None, cells=(QUARTER_WIDTH, QUARTER_WIDTH)):
        BasePlot.__init__(self, ax=ax, cells=cells)
        self.name = "task_subtask_waffle"
        ax = self.ax
        cols = cols or self._COLS

        units = []
        for v, c in zip(sub["value"], sub["count"]):
            units.extend([v] * int(c))
        rows = math.ceil(len(units) / cols)
        self.blank = rows * cols - len(units)

        for i, v in enumerate(units):
            r, c = divmod(i, cols)
            ax.add_patch(mpatches.Rectangle(
                (c + self._PAD / 2, -r + self._PAD / 2), 1 - self._PAD, 1 - self._PAD,
                facecolor=SUBTASK_COLORS[v], edgecolor="none"))

        ax.set_xlim(-0.15, cols + 0.15)
        ax.set_ylim(-(rows - 1) - 0.15, 1.15)
        ax.set_aspect("equal")
        ax.set_axis_off()

        # Keyed in fill order (grouped by parent task), not palette order, so scanning the legend
        # top-to-bottom walks the waffle top-to-bottom. Two columns, on abbreviated labels: three
        # columns would overrun even the 60 mm panel, and at quarter width the unabbreviated two
        # columns already do (48.3 mm of key against a 45.0 mm panel). Spacing is tightened from the
        # matplotlib defaults on purpose: the legend band is the only thing adding height below a
        # landscape grid, so its height is what tunes the panel's aspect.
        # Left-aligned (``loc="upper left"`` anchored at x=0), not centred: the two columns hold labels
        # of different lengths, so centring the block leaves its left edge floating away from the
        # grid's own left edge above it. Aligning both to x=0 lets the key read as part of the panel.
        self.legend({f"{self._LEGEND_ABBREV.get(v, v)} ({int(c)})": SUBTASK_COLORS[v]
                     for v, c in zip(sub["value"], sub["count"])},
                    ncol=2, loc="upper left", bbox_to_anchor=(0.0, 0.0),
                    handlelength=1.1, handletextpad=0.5, labelspacing=0.35,
                    columnspacing=1.0, borderpad=0.1)
        self.label(title="Tasks & subtasks")


# --------------------------------------------------------------------------------------
# Pathogen circle-treemap internals
# --------------------------------------------------------------------------------------
# Area encodes a model's training-set size as n ** AREA_EXPONENT. Raw sizes span 96 to 5,000,000
# compounds, so a linear encoding buries the small datasets, but the earlier 1 + ln(n) compression
# had a worse flaw: every model contributed at least ~5.5 to its group regardless of size, so a
# group's total tracked how MANY models it held rather than how much data — E. coli (12 models,
# 2.4M) outranked P. falciparum (7 models, 6.0M). A power transform makes a tiny model contribute
# tiny area, so model count stops inflating groups, and the exponent dials between data volume
# (1.0) and pure count (towards 0). At 0.5 the pathogen ordering follows data volume and the
# smallest Voronoi cell stays ~9x above the readable floor.
AREA_EXPONENT = 0.5
_PARENT_PAD = 0.5      # grey ring left around a pathogen's packed dots, in dot-radius units
_LABEL_PAD = 0.4       # gap between a circle rim and its label, in font heights
_LABEL_STEP = 0.7      # outward nudge once no angle around the circle is free, in font heights
_INSIDE_FIT = 1.45     # label goes inside the circle when its width <= this x the radius
# Angles tried around a circle, as offsets from the outward radial direction. Sliding the
# label around its own circle keeps it touching the circle it names; pushing it further out
# (the _LABEL_STEP fallback) is what breaks the association, so it is the last resort.
_LABEL_ANGLES = [0, 20, -20, 40, -40, 60, -60, 80, -80, 100, -100, 120, -120, 150, -150, 180]
# Round decade references for the dot-size legend, spanning the observed 96 - 5,000,000 range.
_SIZE_LEGEND_KEYS = [100, 10_000, 1_000_000]
# Axes-fraction footprint of the Source Type legend (lower right), reserved up front so a genus
# label is never placed underneath it. Measured from the rendered panel with slack, since a
# legend's true extent is only knowable after a draw while label positions must be fixed before.
# The size legend needs no entry: it is anchored *below* the axes, clear of the packing.
_LEGEND_BLOCKS = [(0.76, 0.00, 1.00, 0.34)]


def _size_datum(n):
    """Area-proportional datum for a training-set size (a mark's AREA ~ this value)."""
    return float(n) ** AREA_EXPONENT


def _nested_pack(groups):
    """Two-level circle packing with a single global dot scale.

    ``groups`` is ``[(key, [(item_id, datum), ...]), ...]``. Returns
    ``{key: (cx, cy, r_parent, [(x, y, r, item_id), ...])}`` where **every** dot radius is
    exactly ``sqrt(datum)`` regardless of which group it sits in, so dot areas are
    comparable across pathogens.

    Nesting via ``circlify``'s own hierarchy support would not give that: it rescales each
    parent's children to fill the parent, and the packing efficiency differs per group, so
    absolute child sizes drift between parents. Instead each group is packed on its own
    (``circlify`` returns radii ``s * sqrt(datum)`` inside a unit enclosure, ``s`` constant
    per call), the group is scaled by ``1 / s`` to normalise the dots, and the parents are
    then packed by the radius each group needs.
    """
    inner = {}
    for key, items in groups:
        cs = circlify.circlify([{"id": i, "datum": d} for i, d in items],
                               show_enclosure=False)
        s = cs[0].r / math.sqrt(cs[0].ex["datum"])
        unit = 1.0 / s                 # radius the packed dots occupy once normalised
        inner[key] = (unit, unit + _PARENT_PAD, cs)

    # Parents: circlify radius ~ sqrt(datum), so datum = r**2 yields radius = C * r for one
    # global C; dividing through by C lands each parent on the radius its dots require.
    tops = circlify.circlify([{"id": k, "datum": v[1] ** 2} for k, v in inner.items()],
                             show_enclosure=False)
    scale = tops[0].r / inner[tops[0].ex["id"]][1]

    packed = {}
    for t in tops:
        key = t.ex["id"]
        unit, rho, cs = inner[key]
        cx, cy = t.x / scale, t.y / scale
        dots = [(cx + c.x * unit, cy + c.y * unit, c.r * unit, c.ex["id"]) for c in cs]
        packed[key] = (cx, cy, rho, dots)
    return packed


def _pts_per_data(ax, fig):
    """Typographic points per data unit, read off the live ``transData``.

    Requires the figure to have been drawn at the current limits. Taken from the real
    transform rather than computed from ``get_position()`` and the limits: under
    ``aspect="equal"`` the drawn axes is not the nominal rect (here it comes out ~1.2x
    larger), so the analytic estimate silently mis-sizes both the label boxes and the
    size-legend markers.
    """
    (x0, _), (x1, _) = ax.transData.transform([(0, 0), (1, 0)])
    return (x1 - x0) * 72.0 / fig.dpi


def _text_extents(ax, fig, texts, fontsize):
    """Measure each string's rendered ``(width, height)`` in DATA units.

    Real extents from the renderer, not ``len(text) * some_aspect``: a character-count estimate
    is wrong per-string (glyph widths vary, and "S. pneumoniae" is far wider than its count
    suggests next to "E. coli"), and underestimating it is what let neighbouring genus labels
    touch even though the collision test reported them clear. Requires a prior draw.
    """
    renderer = fig.canvas.get_renderer()
    inv = ax.transData.inverted()
    out = {}
    for key, text in texts.items():
        artist = ax.text(0, 0, text, fontsize=fontsize)
        bb = artist.get_window_extent(renderer=renderer)
        artist.remove()
        (x0, y0), (x1, y1) = inv.transform([(bb.x0, bb.y0), (bb.x1, bb.y1)])
        out[key] = (abs(x1 - x0), abs(y1 - y0))
    return out


def _text_box(x, y, w, h, ha, va):
    """Bounding box of a text of size ``(w, h)`` anchored at ``(x, y)`` with ``ha``/``va``."""
    x0 = x if ha == "left" else x - w if ha == "right" else x - w / 2
    y0 = y if va == "bottom" else y - h if va == "top" else y - h / 2
    return (x0, y0, x0 + w, y0 + h)


def _boxes_overlap(a, b):
    return not (a[2] <= b[0] or b[2] <= a[0] or a[3] <= b[1] or b[3] <= a[1])


def _box_hits_circle(box, cx, cy, r):
    nx = min(max(cx, box[0]), box[2])       # box point closest to the circle centre
    ny = min(max(cy, box[1]), box[3])
    return (nx - cx) ** 2 + (ny - cy) ** 2 < r * r


def _nearest_is_own(box, own_key, circles):
    """Whether the circle nearest this label box is the one it names.

    Distance is measured to the rim (``|centre - c| - r``), which is what a reader's eye
    follows. This is the constraint that stops a label being parked next to a neighbour: an
    interior circle simply has no outward direction that satisfies it, and falls back to
    being labelled inside.
    """
    px, py = (box[0] + box[2]) / 2, (box[1] + box[3]) / 2
    nearest = min(circles, key=lambda c: math.hypot(px - c[0], py - c[1]) - c[2])
    return nearest[3] == own_key


def _align_for_angle(deg):
    """Alignment for a label anchored in direction ``deg``, so text grows away from the circle."""
    deg = ((deg + 180) % 360) - 180
    if -45 <= deg < 45:
        return "left", "center"
    if 45 <= deg < 135:
        return "center", "bottom"
    if deg >= 135 or deg < -135:
        return "right", "center"
    return "center", "top"


def _legend_boxes(xlim, ylim):
    """``_LEGEND_BLOCKS`` mapped from axes fractions into data coordinates."""
    (xa, xb), (ya, yb) = xlim, ylim
    return [(xa + f0 * (xb - xa), ya + g0 * (yb - ya),
             xa + f1 * (xb - xa), ya + g1 * (yb - ya))
            for f0, g0, f1, g1 in _LEGEND_BLOCKS]


def _place_labels(circles, labels, extents, line_height, xlim, ylim):
    """Collision-free label placement for ``circles`` = ``[(cx, cy, r, key)]``.

    A circle wide enough to carry its own name (``_INSIDE_FIT``) is labelled *inside*, just
    below the top rim. Everything else is anchored outside: starting from the outward radial
    direction, the label slides around its own rim through ``_LABEL_ANGLES``, stepping
    further out only once a whole lap has failed. Candidates are tried in two tiers —

    1. clear of every placed label AND of every circle, and nearest to its own circle;
    2. clear of every placed label and still nearest to its own circle, but allowed to sit
       over a neighbouring circle.

    Tier 2 exists because in a tight packing some circles have no fully clear slot, and a
    label resting on a neighbour's grey disc (while unambiguously nearest its own) reads far
    better than one crammed inside a circle too small for it. Labels are drawn with a white
    stroke, so overlap stays legible. Only if both tiers fail does the label go inside
    regardless of fit.

    Largest circles go first, since they are what a reader looks for. Returns
    ``([(x, y, ha, va, text, inside)], [box])``.
    """
    lh = line_height

    def inside(cx, cy, r, w, h):
        x, y = cx, cy + r - 1.15 * h
        return x, y, "center", "center", _text_box(x, y, w, h, "center", "center")

    # Seeding with the legend footprints makes them behave exactly like already-placed labels.
    blocked = _legend_boxes(xlim, ylim)
    placed, out = list(blocked), []
    for cx, cy, r, key in sorted(circles, key=lambda c: -c[2]):
        text = labels[key]
        w, h = extents[key]

        chosen = inside(cx, cy, r, w, h) if w <= _INSIDE_FIT * r else None
        is_inside = chosen is not None

        base = math.degrees(math.atan2(cy, cx)) if math.hypot(cx, cy) > 1e-9 else -90.0
        if chosen is None:
            candidates = []
            for step in range(3):
                off = (_LABEL_PAD + step * _LABEL_STEP) * lh
                for delta in _LABEL_ANGLES:
                    deg = base + delta
                    a = math.radians(deg)
                    x, y = cx + math.cos(a) * (r + off), cy + math.sin(a) * (r + off)
                    ha, va = _align_for_angle(deg)
                    box = _text_box(x, y, w, h, ha, va)
                    if any(_boxes_overlap(box, b) for b in placed):
                        continue
                    if not _nearest_is_own(box, key, circles):
                        continue
                    clear = not any(_box_hits_circle(box, *c[:3]) for c in circles)
                    candidates.append((0 if clear else 1, step, abs(delta), x, y, ha, va, box))
            if candidates:
                _, _, _, x, y, ha, va, box = min(candidates, key=lambda c: c[:3])
                chosen = (x, y, ha, va, box)

        mode = "inside" if is_inside else "outside"

        if chosen is None:
            # Nothing unambiguous adjacent to the circle. Forcing the label inside is wrong here
            # — these are the *small* circles, six times too narrow for their own name, so the
            # text overhangs and collides with whatever else was forced inside nearby (which is
            # exactly how Enterobacter and Campylobacter ended up on top of each other). Instead
            # reach further out to genuinely free space and connect it with a leader line, which
            # makes ownership explicit rather than relying on proximity.
            # Only other labels are treated as obstacles here. The circles are not: these are
            # interior circles with no clear radial escape, and since the leader line already
            # establishes which circle the label belongs to, resting on a grey disc costs
            # nothing (the white stroke keeps it readable). Requiring clear space as well is
            # what left Enterobacter and Campylobacter with nowhere to go.
            for step in range(3, 14):
                off = (_LABEL_PAD + step * _LABEL_STEP) * lh
                for delta in _LABEL_ANGLES:
                    deg = base + delta
                    a = math.radians(deg)
                    x, y = cx + math.cos(a) * (r + off), cy + math.sin(a) * (r + off)
                    ha, va = _align_for_angle(deg)
                    box = _text_box(x, y, w, h, ha, va)
                    if not any(_boxes_overlap(box, b) for b in placed):
                        chosen, mode = (x, y, ha, va, box), "leader"
                        break
                if chosen:
                    break

        if chosen is None:              # give up: overhang the circle as a last resort
            chosen, mode = inside(cx, cy, r, w, h), "inside"

        x, y, ha, va, box = chosen
        placed.append(box)
        out.append((x, y, ha, va, text, mode, (cx, cy, r)))
    # Only the label boxes go back: the legend footprints are already inside the limits, so
    # returning them would stop the axis ever fitting tighter than its provisional frame.
    return out, placed[len(blocked):]


class PathogenTreemapPlot(BasePlot):
    """Circle-treemap of models per priority pathogen.

    One circle per pathogen, packed with ``circlify``, holding one dot per model. Dot **area**
    encodes that model's training-set size as ``n ** AREA_EXPONENT`` (see ``_size_datum``) on a single
    global scale, so a dot means the same thing in every circle; a pathogen's circle is only an
    ENCLOSURE, sized to whatever radius the packing needs plus a constant ring — its area encodes
    nothing, and because a packed enclosure grows with the number of dots inside it, it tracks
    model count no matter what the dots encode. Dot colour encodes Source Type. Genus labels sit outside their circle, pushed radially outward
    until they clear every other label and circle.

    Which models appear is driven by ``training_sizes_path``, not by the Airtable Target
    Organism field: that config file is the curated model-to-pathogen mapping (rows whose
    organism assignment was too generic to be meaningful were removed by hand). It therefore
    shows fewer models than the Target Organism bar panel, which still counts every Airtable
    annotation.

    Parameters
    ----------
    df                  : the (Status == "Ready") metadata DataFrame — supplies Source Type.
    pathogens_path      : CSV with ``pathogen`` / ``code`` columns (priority pathogens); sets
                          the display label and which codes are eligible.
    training_sizes_path : CSV with ``eosid`` / ``pathogen`` (code) / ``training_size`` — one
                          row per dot.
    cells               : footprint on the reference grid as ``(rows, cols)``.
    """

    def __init__(self, ax=None, df=None, pathogens_path=None, training_sizes_path=None,
                 cells=(3, 3)):
        BasePlot.__init__(self, ax=ax, cells=cells)
        self.name = "pathogen_circles"
        import pandas as pd

        ax = self.ax
        pathogens = pd.read_csv(pathogens_path)
        sizes = pd.read_csv(training_sizes_path)

        source_type = dict(zip(df["Identifier"], df["Source Type"].fillna("External")))
        display = {r["code"]: abbrev(r["pathogen"]) for _, r in pathogens.iterrows()}

        # One group per pathogen, one item per model row. Ordered by model count so the
        # circlify layout is stable across runs.
        groups = []
        for code in pathogens["code"]:
            g = sizes[sizes["pathogen"] == code]
            if g.empty:
                continue
            groups.append((code, [(r["eosid"], _size_datum(r["training_size"]))
                                  for _, r in g.iterrows()]))
        groups.sort(key=lambda kv: -len(kv[1]))

        packed = _nested_pack(groups)
        circles = [(cx, cy, rho, code) for code, (cx, cy, rho, _) in packed.items()]

        for code, (cx, cy, rho, dots) in packed.items():
            ax.add_patch(mpatches.Circle((cx, cy), rho, facecolor=hue("silver", lighten=0.15),
                                         edgecolor="white", zorder=1))
            for dx, dy, dr, eosid in dots:
                color = SOURCE_TYPE_COLORS.get(source_type.get(eosid), BAR_DEFAULT)
                ax.add_patch(mpatches.Circle((dx, dy), dr, facecolor=color,
                                             edgecolor="white", linewidth=0.3, zorder=3))

        # Labels sit outside the circles, so the axis extent depends on where they land while
        # their size in data units depends on the extent. Iterate: apply the limits, draw to
        # read the true data-to-points scale, place, refit. Three passes is ample — the extent
        # stops growing once the labels stop getting wider.
        x0 = min(c[0] - c[2] for c in circles)
        x1 = max(c[0] + c[2] for c in circles)
        y0 = min(c[1] - c[2] for c in circles)
        y1 = max(c[1] + c[2] for c in circles)
        slack = 0.10 * max(x1 - x0, y1 - y0)
        xlim, ylim = (x0 - slack, x1 + slack), (y0 - slack, y1 + slack)
        # Centred anchor on purpose: the packing is square, so anchor="N" would leave no
        # room under it and the size legend would land on the bottom circles. Centring
        # splits the slack, putting half of it just below the packing.
        ax.set_aspect("equal")
        ax.set_axis_off()

        # The size-legend band is reserved inside ylim from the start. Adding it afterwards
        # would rescale the axis after the labels were placed, so every measured text extent
        # would understate the rendered width and neighbouring genus labels would touch.
        rmax = max(math.sqrt(_size_datum(k)) for k in _SIZE_LEGEND_KEYS)

        for _ in range(4):
            ax.set_xlim(*xlim)
            ax.set_ylim(*ylim)
            self.fig.canvas.draw()
            per_data = _pts_per_data(ax, self.fig)
            extents = _text_extents(ax, self.fig, display, stylia.FONTSIZE_SMALL)
            _texts, boxes = _place_labels(circles, display, extents,
                                          stylia.FONTSIZE_SMALL / per_data, xlim, ylim)
            bx0 = min([x0] + [b[0] for b in boxes])
            bx1 = max([x1] + [b[2] for b in boxes])
            by0 = min([y0] + [b[1] for b in boxes])
            by1 = max([y1] + [b[3] for b in boxes])
            margin = 0.02 * max(bx1 - bx0, by1 - by0)
            band = 2 * rmax + 3.6 * (stylia.FONTSIZE_SMALL / per_data)
            xlim = (bx0 - margin, bx1 + margin)
            ylim = (by0 - margin - band, by1 + margin)

        # Final placement AT the converged limits, keeping them. The loop above refits the
        # limits after each placement, so its last result was decided against stale limits —
        # and the reserved legend footprint, being an axes fraction, would have moved with them.
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        self.fig.canvas.draw()
        per_data = _pts_per_data(ax, self.fig)
        extents = _text_extents(ax, self.fig, display, stylia.FONTSIZE_SMALL)
        texts, _boxes = _place_labels(circles, display, extents,
                                      stylia.FONTSIZE_SMALL / per_data, xlim, ylim)

        # Kept for inspection: overlapping genus labels are the recurring failure mode here.
        self.labels = texts
        # A white stroke keeps every label legible where it crosses a dot or a neighbouring
        # disc, without punching an opaque hole in the shapes underneath.
        for x, y, ha, va, text, mode, (ccx, ccy, cr) in texts:
            if mode == "leader":
                # Thin line from the label back to its circle's rim, so a label parked in
                # distant free space still reads as belonging to that circle.
                dx, dy = ccx - x, ccy - y
                d = math.hypot(dx, dy)
                if d > cr:
                    ax.plot([x, ccx - dx / d * cr], [y, ccy - dy / d * cr],
                            color=hue("silver"), linewidth=0.5, zorder=3)
            ax.text(x, y, text, ha=ha, va=va, fontsize=stylia.FONTSIZE_SMALL, zorder=4,
                    path_effects=[patheffects.withStroke(linewidth=1.2, foreground="white")])

        # Two legends. Colour = Source Type, inside at lower right over empty space (its
        # footprint is reserved in _LEGEND_BLOCKS so no genus label lands under it).
        leg = self.legend(SOURCE_TYPE_COLORS)
        ax.add_artist(leg)

        # Size = training-set size, drawn as nested circles in a band below the packing. Drawing
        # them in DATA coordinates makes them exactly to scale with the dots by construction — no
        # points conversion to get wrong. Nested (shared bottom tangent) rather than in a row
        # because under a power transform the keys span a 10x radius range, and a row of
        # to-scale markers that wide takes over the panel.
        h = stylia.FONTSIZE_SMALL / per_data
        radii = [math.sqrt(_size_datum(k)) for k in _SIZE_LEGEND_KEYS]
        floor = ylim[0] + 2.2 * h                # shared bottom tangent, inside the band
        cx = xlim[0] + rmax + 0.5 * h
        for key, r in sorted(zip(_SIZE_LEGEND_KEYS, radii), key=lambda t: -t[1]):
            ax.add_patch(mpatches.Circle((cx, floor + r), r, facecolor="none",
                                         edgecolor=hue("silver"), linewidth=0.7, zorder=3))
            top = floor + 2 * r
            ax.plot([cx, cx + rmax + 0.6 * h], [top, top], color=hue("silver"),
                    linewidth=0.5, zorder=3)
            ax.text(cx + rmax + 0.9 * h, top, f"{key:,}", ha="left", va="center",
                    fontsize=stylia.FONTSIZE_SMALL, zorder=4)
        ax.text(cx - rmax, ylim[0] + 0.4 * h, "training compounds", ha="left", va="bottom",
                fontsize=stylia.FONTSIZE_SMALL, zorder=4)

        self.label(title="Models targeting priority pathogens")


# --------------------------------------------------------------------------------------
# Sizing for the two stacked subtask panels (see _stack_cells and save_metadata_figures)
# --------------------------------------------------------------------------------------
#: The two cross-tab fields drawn as one vertical unit, TOP FIRST. Output (4 bars) goes above
#: Source Type (3 bars), and the order is what decides which panel drops its axis title: only the
#: upper one does, so the unit reads "Number of models" once, under the whole block, where a reader
#: expects the shared unit of a stacked pair to be. Any cross-tab field NOT listed here gets a
#: default footprint and its own axis instead of being sized into the unit.
_SUBTASK_STACK_ORDER = ("Output", "Source Type")

#: HEIGHT band of a StackedFieldBarPlot, per ``show_xlabel``: figure height minus AXES height, in mm.
#: Everything ``tight_layout`` spends on x tick labels and the optional axis title.
#:
#: MEASURED with ``tools/probe_stack_geometry.py`` and **constant to four decimals** over 19-23 mm and
#: over both fields — which is the whole reason the solver works in axes height. Re-run the probe if
#: the font sizes or the axis labels change. Dropping the axis title is worth 3.40 mm of it.
_AXES_BAND_MM = {True: 11.6634, False: 8.2659}

#: Saved page size minus FIGURE size, in mm — savefig's tight-bbox pad, measured with the same probe.
#: Height is the same either way. Width differs by 0.03 mm between the two panels because
#: ``tight_layout`` sizes each one's y tick-label column to its own content; the value here is the
#: WIDER panel's, so sizing on it keeps NEITHER over the width budget.
_PAGE_PAD_MM = 1.271
_PAGE_PAD_W_MM = 1.247

#: One canvas pixel, in mm — the finest size difference the layout can express (see ``_LAYOUT_DPI``).
#: The figure size is the requested footprint rounded to the NEAREST pixel, not floored, so a panel can
#: come out up to HALF a pixel bigger than asked for. Since the budget is a ceiling, every solved size
#: is backed off by that half pixel: without it the width lands 0.02 mm over, which is nothing on the
#: page but does mean the panel no longer honours a number someone else's layout is built on.
_PIXEL_MM = 25.4 / _LAYOUT_DPI

#: TARGET CROP size of the stacked subtask unit, in mm — the space the pair occupies on the page,
#: which is what a layout budget is actually about. Width is per panel; height is the two panels
#: TOGETHER. Requested 2026-08-05.
#:
#: NOTE these are saved-page sizes, not ``cells`` footprints. The tick labels and the axis title are
#: drawn outside the figure canvas, so ``bbox_inches="tight"`` writes a page LARGER than the
#: footprint — by ~1.27 mm in each direction here. ``_stack_cells`` works backwards from these targets
#: to the footprints, and the script re-reads the saved PDFs' own ``/MediaBox`` afterwards
#: (``plotting_base.pdf_page_mm``) and flags any overrun.
_SUBTASK_STACK_WIDTH_MM = 52.0    # saved page width of EACH panel
_SUBTASK_STACK_HEIGHT_MM = 50.0   # saved page heights of the two panels SUMMED

#: Bar geometry inside a StackedFieldBarPlot, needed to solve for equal bar thickness across two
#: panels with different bar counts. ``stacked_hbar`` leaves the y axis autoscaled, so the view
#: spans the bars' own extent — ``n - 1`` between the first and last centre plus half a bar
#: overhanging at each end — inflated by matplotlib's 5% margin top and bottom. Thickness is
#: therefore ``0.8 x axes_height / _y_span(n)``, NOT ``0.8 x axes_height / n``: at 3 bars the overhang
#: and the margins are worth 10% of the axes, and using the bar count instead leaves the two panels'
#: bars 3% apart — a small error, but in the one quantity this pair is being sized to equalise.
_BAR_HEIGHT_DATA = 0.8        # matplotlib's barh default, which stacked_hbar does not override
_Y_MARGIN = 0.05              # rcParams["axes.ymargin"]; stylia leaves it at the matplotlib default


def _y_span(n_bars):
    """A StackedFieldBarPlot's y view height in DATA units, for ``n_bars`` bars."""
    return (n_bars - 1 + _BAR_HEIGHT_DATA) * (1 + 2 * _Y_MARGIN)


def _stack_axes_heights(n_top, n_bot):
    """Axes heights in mm for the two stacked subtask panels, as ``(top, bottom)``.

    The whole sizing problem, solved in the one quantity that behaves. Page height is figure height
    plus a constant pad, and figure height is axes height plus a constant band, so the height budget
    buys ``_SUBTASK_STACK_HEIGHT_MM`` minus two pads and two bands of axes — and that is what gets
    divided between the panels.

    It is divided in the ratio of their **y spans**, because drawn bar thickness is
    ``_BAR_HEIGHT_DATA x axes_height / _y_span(n)`` — so equal thickness means equal
    ``axes_height / _y_span(n)``, nothing more. Note the split is by span, NOT by bar count: the
    half-bar overhang and the 5% margins make a 3-bar panel's span 2.8 units rather than 3, and
    ignoring that is a 3% error in the very thing being equalised.

    Only the LOWER panel carries the axis title, so it spends 3.40 mm more of its figure on the band;
    the two footprints therefore come out nearly equal even though the *axes* split is 4:3.
    """
    spans = (_y_span(n_top), _y_span(n_bot))
    axes_total = (_SUBTASK_STACK_HEIGHT_MM - 2 * _PAGE_PAD_MM - _PIXEL_MM
                  - _AXES_BAND_MM[False] - _AXES_BAND_MM[True])
    return tuple(axes_total * s / sum(spans) for s in spans)


def _stack_cells(n_top, n_bot):
    """Footprints for the two stacked subtask panels, as ``((rows, cols), (rows, cols))``.

    The solved axes heights (:func:`_stack_axes_heights`) plus each panel's own band give the figure
    height, which is what ``cells`` expresses — so the footprint is not a fudge factor here, it is the
    figure. Width is the page budget less the wider panel's measured pad, the same value for both so
    their axes frames stay registered.

    Returns fractional cells — the documented off-grid exception, unavoidable here since the sizes
    are physical millimetres rather than grid multiples.
    """
    cols = (_SUBTASK_STACK_WIDTH_MM - _PAGE_PAD_W_MM - _PIXEL_MM / 2) / CELL_MM
    return tuple(((axes_mm + _AXES_BAND_MM[show_xlabel]) / CELL_MM, cols)
                 for axes_mm, show_xlabel in zip(_stack_axes_heights(n_top, n_bot), (False, True)))


def _stack_thickness(n_top, n_bot):
    """Predicted drawn bar thickness in mm for the two panels, as ``(top, bottom)``.

    Equal by construction — :func:`_stack_axes_heights` solves for exactly that — so this is the
    *predicted* pair the script checks the *measured* pair
    (``StackedFieldBarPlot.measure_geometry``) against. A divergence means ``_AXES_BAND_MM`` no longer
    describes what the renderer does, and the probe needs re-running.
    """
    return tuple(_BAR_HEIGHT_DATA * axes_mm / _y_span(n)
                 for axes_mm, n in zip(_stack_axes_heights(n_top, n_bot), (n_top, n_bot)))


#: The horizontal box row occupies 4/6 of the page width. Its height is declared at the top of this
#: module (``_BOX_ROW_HEIGHT_CELLS``) because the Biomedical Area strip shares it.
_BOX_ROW_WIDTH_CELLS = 4          # 4/6 of the page width = 120 mm for the whole row

# Measured WIDTH furniture for a horizontal TaskMetricBoxPlot (mm), the counterpart of the height
# constants above: the task tick-label column costs 14.2 mm, and it is fixed in points, so it does not
# scale with the footprint. Verified linear to within 0.5 mm over declared widths of 33-54 mm.
_BOX_FURN_Y_MM = 19.43            # with task tick labels (show_y=True)
_BOX_FURN_NO_Y_MM = 5.24          # without (show_y=False)
# Crop pad differs between the two: the labelled panel's tick-label column sits slightly outside what
# tight_layout reserves for it, so its tight bbox runs ~1.4 mm wider than a bare panel's.
_BOX_CROP_PAD_Y_MM = 2.57
_BOX_CROP_PAD_MM = 1.15


def _box_row_cells(total_mm, n_panels, height_cells=_BOX_ROW_HEIGHT_CELLS):
    """Footprints for a row of ``n_panels`` horizontal box panels sharing one category axis.

    Sized so the panels' tight crops **sum to** ``total_mm`` and every panel gets the **same drawn
    axes width**, which is the point of the exercise: the leftmost panel has to be wider than the
    others by exactly its tick-label column, or the three metric axes come out different lengths and
    the row reads as though the metrics had been scaled differently.

    Returns ``(footprints, axes_width_mm)``. The second value is what the caller checks its axis
    labels against — an xlabel wider than the axes overhangs it and *sets* the crop width, silently
    blowing the row's budget (the runtime label is 26.2 mm against a 28.9 mm axes, the tightest of
    the three).
    """
    furniture = (_BOX_FURN_Y_MM + _BOX_CROP_PAD_Y_MM
                 + (n_panels - 1) * (_BOX_FURN_NO_Y_MM + _BOX_CROP_PAD_MM))
    axes_w = (total_mm - furniture) / n_panels
    first = (axes_w + _BOX_FURN_Y_MM) / CELL_MM
    rest = (axes_w + _BOX_FURN_NO_Y_MM) / CELL_MM
    return ([(height_cells, first)]
            + [(height_cells, rest)] * (n_panels - 1)), axes_w


def save_metadata_figures(counts, df, pathogens_path, training_sizes_path, output_dir,
                          top_n=None, cross_tabs=None):
    """Render every metadata panel as its own Nature-sized figure and save each one.

    Each panel is built standalone (``ax=None``), so it sizes itself from its ``cells``
    footprint (3 cm square grid; see ``plotting_base``) and is written as both a raster PNG
    (``output_dir/png/<name>.png``) and a vector PDF (``output_dir/pdf/<name>.pdf``) ready for
    Illustrator. Panel footprints are also recorded in ``output_dir/figure_cells.json`` so the
    intended grid layout survives the tight-crop applied on save.

    No composite figure and no A/B/C panel letters — final placement/ordering happens in
    Illustrator.

    Parameters
    ----------
    counts              : dict field -> counts DataFrame. ``counts["Subtask"]`` must also carry
                          a ``parent`` column (see the script's ordering step).
    df                  : the (Status == "Ready") metadata DataFrame (for the treemap).
    pathogens_path      : path to the priority-pathogens CSV.
    training_sizes_path : path to the curated per-model training-size CSV (drives the treemap).
    output_dir          : directory to write ``png/``, ``pdf/`` and ``figure_cells.json`` into.
    top_n               : optional dict field -> cap (e.g. {"Target Organism": 10}).
    cross_tabs          : optional dict ``field -> wide count DataFrame`` of that field crossed with
                          Subtask (index = field values in draw order, columns = subtasks in
                          stacking order). Each entry adds one :class:`StackedFieldBarPlot` named
                          ``<field>_by_subtask``. The two fields in ``_SUBTASK_STACK_ORDER`` are
                          sized as one stacked unit by ``_stack_cells``; any other field gets a
                          default footprint and its own count axis.
    """
    top_n = top_n or {}
    cross_tabs = cross_tabs or {}

    # Biomedical Area as FOUR GROUPS over Activity prediction models only — the classification and
    # the reasoning behind every membership live in ``BIOAREA_GROUP`` (default.py), and the script
    # does the counting (distinct models per group, not area assignments). Collapsing to four bars
    # also removes the top-N question entirely: there is no tail left to cut, so no tie to break.
    _bioarea = counts["Biomedical Area grouped"]

    # One entry per panel. Footprints (rows, cols) in 3 cm cells, sized for a 183 x 170 mm Nature
    # page: most panels are 2x2 (60 x 60 mm), the two ten-category fields are quarter-width squares,
    # the technical boxes are 45 mm squares, and the pathogen packing is 2.5x2.5 — wide enough for its
    # two legends and 15 genus labels.
    plots = [
        TaskSubtaskBarPlot(sub=counts["Subtask"], cells=(2, 2)),
        # Unit alternative to the bars, on the same SUBTASK_COLORS. Reports any unfilled trailing
        # cells, which the script prints. Quarter width (create_figure width=0.25 = 45.0 mm), with
        # the two stacked subtask panels sized to sit beside it in the other half of the same
        # 90 mm block — see _stack_cells.
        TaskSubtaskWafflePlot(sub=counts["Subtask"],
                              cells=(QUARTER_WIDTH, QUARTER_WIDTH)),
        # Biomedical Area: a 25 mm strip with the group names centred inside the bars and no y axis,
        # the same height as the technical box row so the two form one band. One flat FULL-STRENGTH
        # crimson rather than `catchall_colors`: every bar here is Annotation / Activity prediction,
        # so colour distinguishes nothing WITHIN the panel and is only a cross-panel cue to the
        # Annotation task. Full strength rather than a pale tint because the labels carry a white
        # outline (see hbar) — legibility is the halo's job, not the bar's, so the bar can match the
        # figure's other crimson marks instead of being washed out to accommodate text.
        FieldBarPlot(counts=_bioarea, title="Biomedical Area",
                     cells=_narrow_strip_cells(len(_bioarea)),
                     colors=TASK_COLORS["Annotation"],
                     inside_labels=True, xlabel="Models"),
        FieldBarPlot(counts=counts["Target Organism"], title="Target Organism",
                     n=top_n.get("Target Organism"), cells=QUARTER_SQUARE,
                     color_fn=catchall_colors),
        # The three donuts. One family: 25 mm wide each (75 mm for the set), ring + legend beneath,
        # total in the hole. Each carries its own numbers, so they need no shared key and can be placed
        # apart — but they are designed to be read together, which is why no hue means two different
        # things across them and why the biomedical one separates its groups by PATTERN rather than
        # spending a fourth palette on categories that are all the same task.
        LicenseClassDonutPlot(counts=counts["License grouped"]),
        ArchitectureDonutPlot(df=df),
        BiomedicalAreaDonutPlot(counts=_bioarea),
        # Footprint: the crop width is set by where the iterative label placement lands, not by the
        # packing, so it is NOT monotonic in `cells` — measured crops are 2.3 -> 65.6 mm, 2.4 -> 67.7,
        # 2.5 -> 63.1, 2.6 -> 67.4. 2.5 is the value that comes in under the 65 mm ceiling; do not
        # assume a smaller footprint gives a smaller panel here, re-measure.
        PathogenTreemapPlot(df=df, pathogens_path=pathogens_path,
                            training_sizes_path=training_sizes_path, cells=(2.5, 2.5)),
    ]

    # Technical box row: three HORIZONTAL box-with-swarm panels sharing one task axis, occupying
    # 4/6 of the page width between them. Only the first draws the task tick labels, so it is wider
    # than the other two by exactly its label column and all three metric axes come out the same
    # length. Names and labels derive from RUNTIME_BATCH so they cannot drift from the column
    # plotted. The -1 sentinels are skipped, never imputed, and any task left with nothing keeps a
    # "not measured" slot (the class reports coverage, which the script prints).
    # The third panel is decade-binned circles rather than a box: Output Dimension is heavily tied,
    # so a swarm overplots onto a few x positions and shows jitter instead of data. Same task axis
    # and footprint, so it drops into the row unchanged.
    box_row = [
        (TaskMetricBoxPlot, RUNTIME_COLUMN, f"runtime_{RUNTIME_BATCH}",
         f"Runtime for {RUNTIME_BATCH:,} molecules (s)"),
        (TaskMetricBoxPlot, "Image Size", "image_size", "Docker image size (MB)"),
        (TaskOutputDimensionCirclesPlot, "Output Dimension", "output_dimension",
         "Output dimension"),
    ]
    box_cells, box_axes_mm = _box_row_cells(
        _BOX_ROW_WIDTH_CELLS / CELLS_PER_WIDTH * stylia.SIZE * 25.4, len(box_row))
    for (cls, column, name, xlabel), cells in zip(box_row, box_cells):
        plots.append(cls(df=df, column=column, name=name, xlabel=xlabel,
                         show_y=name == box_row[0][2], cells=cells))

    # Stacked variants: one bar per field value, segmented by subtask, so a single panel carries the
    # joint distribution. No legend — `task_subtask` and the waffle are the shared subtask key.
    #
    # The two named in _SUBTASK_STACK_ORDER are drawn as one vertical unit, SOLVED by _stack_cells
    # from the crop budget (_SUBTASK_STACK_WIDTH_MM x _SUBTASK_STACK_HEIGHT_MM) so the pair fills the
    # space it is allotted and both panels draw bars of the same thickness. Each keeps its OWN
    # autoscaled count axis and its own tick labels; only the UPPER panel's axis TITLE is dropped, so
    # the unit says "Number of models" once, beneath the block.
    #
    # Footprints are what we ask for; crops are what gets written. measure_geometry() reports the
    # real crop and the real bar thickness for each panel below, so neither the budget nor the
    # equal-thickness claim rests on the constants alone.
    unit = [f for f in _SUBTASK_STACK_ORDER if f in cross_tabs]
    stack_cells = {}
    if len(unit) == 2:
        n_bars = [len(cross_tabs[f]) for f in unit]
        top, bot = _stack_cells(*n_bars)
        stack_cells = {unit[0]: top, unit[1]: bot}
        want = dict(zip(unit, _stack_thickness(*n_bars)))
        print(f"\n[subtask stack] target page {_SUBTASK_STACK_WIDTH_MM:g} mm wide x "
              f"{_SUBTASK_STACK_HEIGHT_MM:g} mm for the pair; declared "
              + ", ".join(f"{f} {stack_cells[f][1] * CELL_MM:.2f} x "
                          f"{stack_cells[f][0] * CELL_MM:.2f} mm" for f in unit)
              + f"; axes autoscaled per panel, xlabel only on {unit[1]}")
        print(f"{'':16s}predicted bar thickness "
              + ", ".join(f"{f} {want[f]:.2f} mm ({n} bars)" for f, n in zip(unit, n_bars)))

    for field, table in cross_tabs.items():
        plots.append(StackedFieldBarPlot(
            table=table, colors=SUBTASK_COLORS, legend_kw=None,
            cells=stack_cells.get(field, (1, 2)),
            show_xlabel=field != unit[0] if len(unit) == 2 else True,
            name=f"{field.lower().replace(' ', '_')}_by_subtask"))

    footprints = {}
    measured_stacks = {}
    for p in plots:
        if isinstance(p, _HorizontalTaskPanel):
            p.measure_fit()          # must happen before save() closes the figure
        if isinstance(p, FieldBarPlot) and p.label_texts:
            p.measure_labels()       # likewise
        if isinstance(p, DonutPlot):
            p.pin_ring()             # likewise: needs the laid-out axes, before save closes it
        if isinstance(p, StackedFieldBarPlot):
            p.measure_geometry()     # likewise
            measured_stacks[p.name] = p
        p.save(output_dir)
        footprints[p.name] = list(p.cells)
        if isinstance(p, FieldBarPlot) and p.label_texts:
            print(f"\n[{p.name}] {len(p.label_texts)} bars, names on the bars (no y axis): "
                  f"{p.n_inside} centred inside, {p.n_after} set after a too-short bar; "
                  f"axes {p.axes_mm:.2f} mm vs widest label {p.label_mm:.2f} mm "
                  f"('{p.widest_label}')"
                  + ("" if p.label_mm <= p.axes_mm else
                     "  <-- label OVERHANGS: it sets the crop width, so the strip is no longer "
                     f"{NARROW_STRIP_MM:g} mm. Abbreviate it in BIOAREA_DISPLAY."))
        if isinstance(p, TaskSubtaskWafflePlot):
            print(f"\n[{p.name}] unfilled trailing cells in the last row: {p.blank}")
        if isinstance(p, TaskMetricBoxPlot):
            gaps = {t: c for t, c in p.coverage.items() if c[0] < c[1]}
            print(f"\n[{p.name}] measured per task: "
                  + ", ".join(f"{t} {n}/{tot}" for t, (n, tot) in p.coverage.items()))
            if gaps:
                print(f"{'':{len(p.name) + 3}}unmeasured models skipped (Airtable -1 sentinel), "
                      f"NOT imputed: " + ", ".join(f"{t} {tot - n}" for t, (n, tot) in gaps.items()))
        if isinstance(p, TaskOutputDimensionCirclesPlot):
            print(f"\n[{p.name}] decade bins per task (area-proportional circles): "
                  + "; ".join(f"{t} " + " ".join(f"1e{k}:{c}" for k, c in b.items() if c)
                              for t, b in p.bins.items()))
            if p.n_skipped:
                print(f"{'':{len(p.name) + 3}}{p.n_skipped} non-positive values skipped "
                      f"(no decade), NOT imputed")
        if isinstance(p, _HorizontalTaskPanel):
            fits = p.xlabel_mm <= p.axes_mm
            print(f"{'':{len(p.name) + 3}}axes {p.axes_mm:.2f} mm wide vs xlabel {p.xlabel_mm:.2f} mm"
                  + ("" if fits else "  <-- xlabel OVERHANGS: it now sets the crop width, so the "
                                     "box row no longer sums to its budget. Shorten it."))
        if isinstance(p, DonutPlot):
            fits = p.axes_mm >= DONUT_RING_MM
            print(f"\n[{p.name}] ring pinned to {DONUT_RING_MM:g} mm; this panel's axes had "
                  f"{p.axes_mm:.2f} mm"
                  + ("" if fits else f"  <-- AXES NARROWER THAN THE RING: its legend is squeezing the "
                                     f"axes, so the ring is being clipped. Shorten the longest label."))
        if isinstance(p, ArchitectureDonutPlot):
            print(f"{'':{len(p.name) + 3}}{p.counts} — snapshot, not a trend "
                  f"(45% dual-arch among 2021 models vs 77% among 2026)")

    # What the pair ACTUALLY came out as. Both checks matter and neither is visible from the
    # footprints: the crop is the space the panels take on the page, and the bar thickness is whether
    # the two still read as one unit. A drift in either means the measured furniture constants above
    # have gone stale.
    if len(unit) == 2:
        panels = [measured_stacks.get(f"{f.lower().replace(' ', '_')}_by_subtask") for f in unit]
        if all(panels):
            pages = [pdf_page_mm(os.path.join(output_dir, "pdf", p.name + ".pdf")) for p in panels]
            page_h = sum(h for _, h in pages)
            page_w = max(w for w, _ in pages)
            print(f"\n[subtask stack] saved page "
                  + ", ".join(f"{p.name} {w:.2f} x {h:.2f} mm"
                              for p, (w, h) in zip(panels, pages))
                  + f"; pair {page_w:.2f} x {page_h:.2f} mm vs budget "
                  f"{_SUBTASK_STACK_WIDTH_MM:g} x {_SUBTASK_STACK_HEIGHT_MM:g}")
            for what, got, want, const in (("width", page_w, _SUBTASK_STACK_WIDTH_MM, "_PAGE_PAD_W_MM"),
                                           ("height", page_h, _SUBTASK_STACK_HEIGHT_MM, "_PAGE_PAD_MM")):
                if got > want + 0.1:
                    print(f"{'':16s}<-- OVER BUDGET on {what} by {got - want:.2f} mm: raise "
                          f"{const} and re-run")
            # The band is the figure height the axes did NOT get. The solver divides the budget by
            # these two bands, so a drift here is exactly what makes the bars unequal.
            print(f"{'':16s}axes band measured "
                  + ", ".join(f"{p.name.split('_by_')[0]} {p.fig_h_mm - p.axes_h_mm:.4f} mm"
                              for p in panels)
                  + f" vs constants "
                  + " / ".join(f"{_AXES_BAND_MM[p.show_xlabel]:.4f}" for p in panels))
            bars = [p.bar_mm for p in panels]
            ratio = max(bars) / min(bars)
            print(f"{'':16s}measured bar thickness "
                  + ", ".join(f"{p.name.split('_by_')[0]} {p.bar_mm:.3f} mm "
                              f"(pitch {p.pitch_mm:.2f}, {p.n_bars} bars)" for p in panels)
                  + (f" — EQUAL to {(ratio - 1) * 100:.1f}%" if ratio < 1.01 else
                     f"  <-- UNEQUAL by {(ratio - 1) * 100:.1f}%: re-calibrate _AXES_BAND_MM from "
                     f"the line above, which is what the solver splits the budget by"))

    with open(os.path.join(output_dir, "figure_cells.json"), "w") as f:
        json.dump(footprints, f, indent=2)
