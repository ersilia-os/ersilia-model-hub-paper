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
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.patheffects as patheffects
import stylia

from plotting_base import BasePlot, CELLS_PER_WIDTH
from plotting_colors import (TASK_COLORS, SUBTASK_COLORS, SOURCE_TYPE_COLORS, BAR_DEFAULT,
                             ARCH_COLORS, ARCH_DISPLAY, LICENSE_CLASS_COLORS, REFERENCE_LINE,
                             catchall_colors, license_colors, ordinal_shades, hue)
from plotting_utils import abbrev, box_with_jitter, hbar, stacked_hbar
from voronoi_treemap import polygon_area, polygon_centroid, voronoi_treemap
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
    legend_map : optional ``{label: colour}`` key, for panels whose bar colour encodes something
             the axis labels do not say (e.g. the licence reuse classes). Omit where colour is
             decorative or already implied by the tick labels.
    n      : optional top-N cap; when set, only the first ``n`` rows are shown and the
             title gets a "(top n)" suffix.
    cells  : footprint on the reference grid as ``(rows, cols)`` — taller for panels with
             more bars (see ``save_metadata_figures``).
    """

    def __init__(self, ax=None, counts=None, title="", colors=None, color_fn=None,
                 legend_map=None, n=None, cells=(3, 3)):
        BasePlot.__init__(self, ax=ax, cells=cells)
        self.name = title.lower().replace(" ", "_")
        if n:
            counts = counts.head(n)
            title = f"{title} (top {n})"
        if color_fn is not None:
            colors = color_fn(counts["value"])
        hbar(self.ax, counts["value"].tolist(), counts["count"].tolist(),
             colors=colors if colors is not None else BAR_DEFAULT)
        if legend_map:
            self.legend(legend_map)
        self.label(xlabel="Number of models", ylabel=" ", title=title)


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
             and a 6-entry legend does not fit on a 60 x 30 mm panel. ``{}`` gives the primitive's
             defaults.
    xlim   : count-axis limits, or ``None`` to autoscale.
    cells  : footprint on the reference grid as ``(rows, cols)``.
    """

    def __init__(self, ax=None, table=None, colors=None, name="stacked_field",
                 legend_kw=None, xlim=None, cells=(2, 2)):
        BasePlot.__init__(self, ax=ax, cells=cells)
        self.name = name
        segments = list(table.columns)
        # stacked_hbar wants one {segment: value} dict per row. Its xlim default of (0, 1) is for
        # fraction stacks, so pass the count axis through explicitly (None = autoscale).
        rows = [row.to_dict() for _, row in table.iterrows()]
        stacked_hbar(self.ax, list(table.index), rows, segments, colors, xlim=xlim)
        if legend_kw is not None:
            self.legend({s: colors[s] for s in segments}, **legend_kw)
        self.label(xlabel="Number of models", ylabel=" ", title=name)


class TaskSubtaskBarPlot(BasePlot):
    """Combined Task + Subtask panel: one horizontal bar per subtask, in that subtask's own
    shade of its parent task's hue — the SAME ``SUBTASK_COLORS`` the stacked ``*_by_subtask``
    panels use. Bars are grouped by task and ordered by count within it, so the task grouping
    reads off the hue runs while each bar is named on the axis.

    Deliberately carries **no legend**: every bar is already labelled with its subtask name next
    to its own colour, so this panel *is* the subtask colour key, and placing it beside the
    ``*_by_subtask`` panels lets those drop their legends too (they have no room for one at
    60 x 30 mm). Per-task totals live in ``task_counts.csv``. ``task_subtask_waffle`` is the other
    panel that can serve as the key — it carries its own legend.

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


class TaskMetricBoxPlot(BasePlot):
    """One box-with-dots per Task for a per-model technical metric (runtime, image size, ...).

    Values ``<= 0`` are treated as NOT MEASURED, because Airtable stores a ``-1`` sentinel for a
    benchmark that was never run. They are skipped rather than imputed, and every box's tick label
    carries its own ``n`` so a reader can never mistake a thin sample for a complete one.

    A Task with **no** measurements still gets its slot, marked "not measured" in the neutral colour
    instead of being dropped. That matters for the runtime panel: none of the 19 Sampling models have
    a 1000-molecule benchmark, and silently showing a two-box panel would read as though the hub had
    only two kinds of model. ``self.coverage`` records ``{task: (n_measured, n_total)}`` for the
    caller to report.

    The y axis is logarithmic by default: both metrics used here span more than a decade (runtime
    17.5-1491 s, image size 291-10242 MB) and a linear axis flattens the bulk of the distribution
    against the floor.

    Tick labels are single-line and **rotated 30 degrees**. Stacked horizontally they do not fit: at
    6 pt "Representation" alone is ~20 mm wide, so three of them need the full 60 mm of panel and
    already touched at that size. Rotating is what lets the panel shrink to 45 mm square.

    Parameters
    ----------
    df       : the (Status == "Ready") metadata DataFrame.
    column   : the metric column to draw.
    name     : output file stem.
    ylabel   : y-axis label, including units.
    log      : logarithmic y axis (default True).
    rotation : x tick label rotation in degrees; 0 draws them horizontally (needs a wider panel).
    cells    : footprint on the reference grid as ``(rows, cols)``.
    """

    def __init__(self, ax=None, df=None, column=None, name="task_metric", ylabel="",
                 log=True, rotation=30, cells=SMALL_SQUARE):
        BasePlot.__init__(self, ax=ax, cells=cells)
        self.name = name
        ax = self.ax
        rng = np.random.default_rng(RANDOM_SEED)

        tasks = [t for t in TASK_COLORS if (df["Task"] == t).any()]
        self.coverage = {}
        labels = []
        for i, task in enumerate(tasks):
            sub = df.loc[df["Task"] == task, column]
            vals = sub[sub > 0].to_numpy(dtype=float)
            self.coverage[task] = (len(vals), len(sub))
            if len(vals):
                box_with_jitter(ax, vals, i, TASK_COLORS[task], rng=rng, point_size=3,
                                point_alpha=0.55, jitter_width=0.10)
                labels.append(f"{task} ({len(vals)})")
            else:
                # Keep the slot. An absent box is information; an absent category is a misreading.
                # Placed mid-height in axes coordinates (get_xaxis_transform: data x, axes y) so it
                # sits where the box would have been, and needs no y value on a log scale.
                ax.text(i, 0.5, "not\nmeasured", ha="center", va="center", linespacing=1.15,
                        transform=ax.get_xaxis_transform(), fontsize=stylia.FONTSIZE_SMALL,
                        color=REFERENCE_LINE)
                labels.append(f"{task} (0/{len(sub)})")

        if log:
            ax.set_yscale("log")
        ax.set_xlim(-0.6, len(tasks) - 0.4)
        ax.set_xticks(range(len(tasks)))
        ax.set_xticklabels(labels, rotation=rotation or 0,
                           ha="right" if rotation else "center")
        self.label(ylabel=ylabel, title=name)


def _spread(ys, gap):
    """Push a column of label centres apart so consecutive ones are at least ``gap`` apart.

    ``ys`` are the natural positions for ONE side of a pie, sorted top-down. Returns the adjusted
    positions, keeping the column centred where it started so it does not drift off the panel.
    """
    out = list(ys)
    for i in range(1, len(out)):                       # cascade downwards
        out[i] = min(out[i], out[i - 1] - gap)
    for i in range(len(out) - 2, -1, -1):              # and back up, in case the cascade ran long
        out[i] = max(out[i], out[i + 1] + gap)
    shift = (ys[0] + ys[-1]) / 2 - (out[0] + out[-1]) / 2
    return [y + shift for y in out]


def _decollide(ax, texts):
    """Separate outside pie labels vertically, per side, once they have been drawn.

    A wedge's natural label sits at its mid-angle, which collides as soon as two slices are narrow
    and adjacent: the single non-commercial licence is a 1.7 degree wedge butted against a 13 % one,
    so their labels land ~0.1 data units apart against a ~0.15 unit line height. Sides are handled
    independently, since a left and a right label may share a height harmlessly.
    """
    fig = ax.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    inv = ax.transData.inverted()
    for sign in (1, -1):
        side = [t for t in texts if (1 if t.get_position()[0] >= 0 else -1) == sign]
        if len(side) < 2:
            continue
        side.sort(key=lambda t: -t.get_position()[1])
        heights = []
        for t in side:
            bb = t.get_window_extent(renderer)
            (_, y0), (_, y1) = inv.transform([(bb.x0, bb.y0), (bb.x1, bb.y1)])
            heights.append(abs(y1 - y0))
        ys = _spread([t.get_position()[1] for t in side], max(heights) * 1.05)
        for t, y in zip(side, ys):
            t.set_position((t.get_position()[0], y))


def _pie(ax, labels, values, colors):
    """House-style pie: slices clockwise from 12 o'clock, labelled OUTSIDE with count and share.

    Labelling in place rather than with a legend keeps the panel self-contained, but it only works
    while every slice is wide enough to own a label — so a pie built through this helper should have
    a handful of comparable slices, not a long tail (a 1-of-208 slice is a 1.7 degree wedge that can
    be neither seen nor labelled). Push detail like that to a bar panel instead.
    """
    values = np.asarray(values, dtype=float)
    total = values.sum()
    # A share under 0.5 % rounds to "0 %", which reads as none at all — the single non-commercial
    # licence is 1/208 = 0.5 %. Show "<1 %" instead so a real slice is never labelled zero.
    shares = ["<1%" if 0 < v / total < 0.005 else f"{v / total:.0%}" for v in values]
    _wedges, texts = ax.pie(
        values, radius=1.0, startangle=90, counterclock=False, colors=colors,
        labels=[f"{l}\n{int(v)} ({s})" for l, v, s in zip(labels, values, shares)],
        labeldistance=1.12,
        textprops=dict(fontsize=stylia.FONTSIZE_SMALL, linespacing=1.2, ha="center"),
        wedgeprops=dict(edgecolor="white", linewidth=0.6))
    ax.set_aspect("equal")
    _decollide(ax, texts)


class ArchitecturePiePlot(BasePlot):
    """Share of models whose container is built for ARM64 as well as x86-64.

    A two-slice pie, labelled outside with the count and percentage, so it needs no legend. Note the
    underlying field only ever holds ``AMD64`` or ``AMD64,ARM64`` — there is no ARM-only build, so
    this is "also built for ARM" versus "x86 only", not a three-way split.

    This is a **snapshot**, not a trend: it reflects the metadata as staged, and the share has moved
    a long way (45 % dual-arch among 2021 models, 77 % among 2026 ones), so the incorporation-date
    range of the source data belongs in any caption.

    Parameters
    ----------
    df    : the (Status == "Ready") metadata DataFrame.
    cells : footprint on the reference grid as ``(rows, cols)``.
    """

    def __init__(self, ax=None, df=None, cells=SMALL_SQUARE):
        BasePlot.__init__(self, ax=ax, cells=cells)
        self.name = "docker_architecture"

        counts = (df["Docker Architecture"].map(ARCH_DISPLAY)
                  .value_counts()
                  .reindex(list(ARCH_COLORS))
                  .dropna())
        self.counts = counts.astype(int).to_dict()
        _pie(self.ax, list(counts.index), counts.to_numpy(dtype=float),
             [ARCH_COLORS[k] for k in counts.index])
        self.label(title="Docker build architecture")


class LicenseClassPiePlot(BasePlot):
    """Licence composition as a pie — the compositional counterpart to the ``license`` bar panel.

    Slices are the four **reuse classes**, not the ten individual licences. That is a deliberate
    limit, not a shortcut: four of the ten licences cover exactly one model each, which on 208 models
    is a **1.7 degree wedge** — invisible, unlabellable, and impossible to tell apart from the other
    three. The per-licence detail is what the bar panel is for; a pie can only carry the four-way
    split. Even here Non-commercial is a single model (0.5 %), so its slice is a hairline whose label
    is the only thing a reader will actually see.

    Ordered by count so the hairline slice ends up adjacent to the 12 o'clock start rather than
    wedged between two large ones.

    Parameters
    ----------
    counts : the grouped licence DataFrame (``value``, ``count``, ``class``) from the script.
    cells  : footprint on the reference grid as ``(rows, cols)``.
    """

    def __init__(self, ax=None, counts=None, cells=SMALL_SQUARE):
        BasePlot.__init__(self, ax=ax, cells=cells)
        self.name = "license_class_pie"

        by_class = (counts.groupby("class")["count"].sum()
                    .reindex(list(LICENSE_CLASS_COLORS))
                    .dropna()
                    .sort_values(ascending=False))
        self.counts = by_class.astype(int).to_dict()
        _pie(self.ax, list(by_class.index), by_class.to_numpy(dtype=float),
             [LICENSE_CLASS_COLORS[k] for k in by_class.index])
        self.label(title="License reuse class")


class TagCloudPlot(BasePlot):
    """Word cloud of the free-text ``Tag`` field — what the hub actually covers, at finer grain
    than the curated Biomedical Area / Target Organism fields.

    Font size encodes how many models carry the tag and colour is an ordinal periwinkle shade of the
    same quantity (:func:`plotting_colors.ordinal_shades`), so colour is redundant with size by
    design — it makes the ranking scannable without adding a second claim.

    **This is the one raster panel in the repo.** ``wordcloud`` renders to a bitmap, which goes into
    the axes through ``imshow``, so this panel's PDF embeds an image instead of editable vector text —
    the one documented exception to the vector-PDF rule. It is generated at ``_PX`` pixels square so
    it stays crisp at the 600 dpi the PNG is saved at.

    **A word cloud is not a quantitative encoding.** A word's inked area depends on how many
    characters it has and whether it was rotated, so "Antimicrobial activity" looks far more than 49x
    "GPCR" (1). Read ranks from it, never ratios; the exact counts live in ``tag_counts.csv``. Layout
    is seeded with ``RANDOM_SEED`` so successive runs are identical.

    Two settings exist purely to keep the long tail legible, and both cost encoding fidelity —
    accepted, because the panel is ordinal anyway. Frequencies are fed in as **sqrt(count)**, which
    compresses the 49:1 count range to 7:1, and ``relative_scaling`` stays at wordcloud's 0.5, which
    blends frequency-proportional sizing with rank-only sizing and lifts the smallest words further.
    Together they make the size scale **ordinal, not proportional** — do not read a ratio off it.

    The numbers behind those choices, all with 59 tags placed: raw counts at 60 mm put **30 tags
    below 5 pt**, with the ten single-model tags at **2.9 pt**; sqrt at 60 mm still leaves 11 below;
    sqrt at 90 mm puts the smallest at **5.5 pt** and none below. (``relative_scaling=1.0`` would
    restore font ~ sqrt(count) and hence area ~ count, but its 7:1 font range drops the tail back to
    3.6 pt even at 90 mm.) ``self.min_font_pt`` / ``self.below_floor`` re-measure this every run
    rather than trusting it, and ``self.n_placed`` catches wordcloud silently dropping words it
    cannot fit.

    Parameters
    ----------
    counts : DataFrame with columns ``value`` (tag) and ``count``, already sorted descending.
    cells  : footprint on the reference grid as ``(rows, cols)``. Do not shrink below (3, 3) without
             re-checking ``below_floor`` — legibility here is set by the panel's physical size.
    """

    _PX = 1400              # raster side in pixels (600 dpi over a 60 mm panel, more over a larger one)
    #: wordcloud's own default: a compromise between font size tracking frequency (1.0) and tracking
    #: rank only (0.0). Kept at 0.5 for the tail legibility documented above.
    _RELATIVE_SCALING = 0.5
    #: Smallest type this repo will call legible in print.
    _PRINT_FLOOR_PT = 5.0

    def __init__(self, ax=None, counts=None, cells=(3, 3)):
        BasePlot.__init__(self, ax=ax, cells=cells)
        self.name = "tag_cloud"
        from matplotlib import font_manager
        from wordcloud import WordCloud

        freq = {v: float(np.sqrt(n)) for v, n in zip(counts["value"], counts["count"])}
        # Colour by rank, matching the bar-panel gradient: darkest = most models.
        shades = dict(zip(counts["value"], ordinal_shades(len(counts))))
        self.n_tags = len(freq)

        # Same typeface as every other panel. wordcloud needs a file path, not a family name, so the
        # family stylia selected is resolved through matplotlib's font manager rather than hardcoded.
        font_path = font_manager.findfont(
            font_manager.FontProperties(family=plt_rc_font_family()))

        wc = WordCloud(
            width=self._PX, height=self._PX, background_color="white", mode="RGB",
            font_path=font_path, prefer_horizontal=0.9, max_words=len(freq),
            relative_scaling=self._RELATIVE_SCALING, random_state=RANDOM_SEED,
            color_func=lambda word, **kw: mcolors.to_hex(shades[word]),
        ).generate_from_frequencies(freq)

        # Font sizes come back in canvas pixels; convert to points at this panel's physical width.
        panel_mm = cells[1] / CELLS_PER_WIDTH * 180.0
        pt_per_px = panel_mm / self._PX / 25.4 * 72.0
        sizes = [fs * pt_per_px for _entry, fs, _pos, _orient, _c in wc.layout_]
        self.n_placed = len(sizes)
        self.min_font_pt = min(sizes) if sizes else 0.0
        self.below_floor = sum(1 for s in sizes if s < self._PRINT_FLOOR_PT)

        self.ax.imshow(wc.to_array(), interpolation="bilinear")
        self.ax.set_axis_off()
        self.label(title="Model tags")


def plt_rc_font_family():
    """The concrete font family stylia's style selected (e.g. ``Arial``), for ``wordcloud``."""
    import matplotlib
    family = matplotlib.rcParams["font.family"][0]
    if family in ("sans-serif", "serif", "monospace", "cursive", "fantasy"):
        return matplotlib.rcParams[f"font.{family}"][0]
    return family


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

    The legend is also what makes the panel square. The 16 x 13 grid alone crops to 1.23 (wide and
    short); the legend band adds height below it and lands the whole panel at **1.05** — so the grid
    stays landscape. Going the other way is counter-productive here: 13 x 16 (also exact) drops to
    0.78, because a portrait grid is narrower than the legend, which then sets the width while the
    extra rows add height. Measured crop aspects: 13 cols 0.78, 14 cols 0.83, 15 cols 0.94, **16
    cols 1.05**, 18 cols 1.22.

    ``self.blank`` records how many trailing cells of the grid went unfilled, so a caller can flag a
    ragged last row. 16 columns is also one of the few counts that divides 208 exactly (13 full
    rows); 14, 15, 17 and 18 all leave a ragged row.

    Parameters
    ----------
    sub   : subtask counts DataFrame with columns ``value`` (short display label) and ``count``,
            already ordered by parent task then count within task.
    cols  : squares per row. The default 16 divides the current 208 models into exactly 13 full
            rows; a different total leaves ``self.blank`` cells empty in the last row.
    cells : footprint on the reference grid as ``(rows, cols)``.
    """

    _COLS = 16
    _PAD = 0.16       # white gap between squares, in square-widths

    def __init__(self, ax=None, sub=None, cols=None, cells=(2, 2)):
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
        # top-to-bottom walks the waffle top-to-bottom. Two columns: at this font "Property
        # prediction (39)" is ~25 mm wide, so three columns would overrun the 60 mm panel.
        # Spacing is tightened from the matplotlib defaults on purpose: the legend band is the only
        # thing adding height below a landscape grid, so its height is what tunes the panel's aspect.
        self.legend({f"{v} ({int(c)})": SUBTASK_COLORS[v]
                     for v, c in zip(sub["value"], sub["count"])},
                    ncol=2, loc="upper center", bbox_to_anchor=(0.5, 0.0),
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


class PathogenVoronoiPlot(BasePlot):
    """Voronoi treemap of models per priority pathogen — the space-filling alternative to
    :class:`PathogenTreemapPlot`, kept alongside it so the two can be compared.

    Two levels, both area-accurate: the panel rectangle is tessellated into one region per
    pathogen with area proportional to that pathogen's total training compounds, and each region
    is then tessellated into one cell per model with area proportional to that model's training
    set. Unlike circle packing a treemap wastes no space, so every cell area is literally its
    share of the panel and the pathogen regions carry a real quantity rather than being mere
    enclosures. Cell colour encodes Source Type.

    ``self.area_error`` records the largest relative area error over every tessellation, and
    ``self.tiny`` lists cells whose final area is too small to read; callers should surface both
    rather than imply the areas are exact.

    Parameters
    ----------
    df                  : the (Status == "Ready") metadata DataFrame — supplies Source Type.
    pathogens_path      : CSV with ``pathogen`` / ``code`` columns.
    training_sizes_path : CSV with ``eosid`` / ``pathogen`` (code) / ``training_size``.
    cells               : footprint on the reference grid as ``(rows, cols)``.
    """

    #: Cells below this fraction of the panel are reported in ``self.tiny`` as unreadable.
    MIN_READABLE_FRACTION = 1e-4

    def __init__(self, ax=None, df=None, pathogens_path=None, training_sizes_path=None,
                 cells=(3, 3), area_metric="power"):
        BasePlot.__init__(self, ax=ax, cells=cells)
        self.name = "pathogen_voronoi" if area_metric == "power" else "pathogen_voronoi_linear"
        import pandas as pd

        if area_metric not in ("linear", "power"):
            raise ValueError(f"area_metric must be 'linear' or 'power', got {area_metric!r}")
        # "linear" is faithful to compound count but illegible here (two models hold ~70% of all
        # training data); "power" applies the same n ** AREA_EXPONENT the circle panel uses.
        weigh = (lambda v: float(v)) if area_metric == "linear" else _size_datum

        ax = self.ax
        pathogens = pd.read_csv(pathogens_path)
        sizes = pd.read_csv(training_sizes_path)

        source_type = dict(zip(df["Identifier"], df["Source Type"].fillna("External")))
        display = {r["code"]: abbrev(r["pathogen"]) for _, r in pathogens.iterrows()}

        groups = []
        for code in pathogens["code"]:
            g = sizes[sizes["pathogen"] == code]
            if g.empty:
                continue
            groups.append((code, [(r["eosid"], weigh(r["training_size"]))
                                  for _, r in g.iterrows()]))
        groups.sort(key=lambda kv: -sum(v for _, v in kv[1]))

        boundary = np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]])
        errors = []

        # Level 1 — one region per pathogen, area proportional to its total training compounds.
        regions, err = voronoi_treemap([sum(v for _, v in items) for _, items in groups],
                                       boundary, seed=RANDOM_SEED, iterations=600)
        errors.append(err)

        self.tiny = []
        for (code, items), region in zip(groups, regions):
            if len(region) < 3:
                continue
            # Level 2 — one cell per model inside that pathogen's region.
            sub, sub_err = voronoi_treemap([v for _, v in items], region,
                                           seed=RANDOM_SEED, iterations=600)
            errors.append(sub_err)
            for (eosid, value), cell in zip(items, sub):
                if len(cell) < 3:
                    self.tiny.append((code, eosid, 0.0))
                    continue
                frac = polygon_area(cell)
                if frac < self.MIN_READABLE_FRACTION:
                    self.tiny.append((code, eosid, frac))
                ax.add_patch(mpatches.Polygon(
                    cell, closed=True, linewidth=0.3, edgecolor="white",
                    facecolor=SOURCE_TYPE_COLORS.get(source_type.get(eosid), BAR_DEFAULT),
                    zorder=2))
            # Region outline on top, so pathogen grouping reads over the per-model cells.
            ax.add_patch(mpatches.Polygon(region, closed=True, fill=False, linewidth=1.1,
                                          edgecolor="white", zorder=3))

        self.area_error = max(errors)

        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_aspect("equal")
        ax.set_axis_off()
        # Draw once so _text_extents below can measure against the live transform.
        self.fig.canvas.draw()

        # Region labels sit at the centroid, but only where the region is actually wide enough
        # to hold the text — a treemap fills the panel, so there is no outside to escape to.
        # Anything skipped is recorded in self.unlabelled for the caller to report.
        self.unlabelled = []
        extents = _text_extents(ax, self.fig, display, stylia.FONTSIZE_SMALL)
        for (code, _items), region in zip(groups, regions):
            if len(region) < 3:
                self.unlabelled.append(code)
                continue
            text = display[code]
            w, h = extents[code]
            span = region[:, 0].max() - region[:, 0].min()
            tall = region[:, 1].max() - region[:, 1].min()
            if w > span or h > tall:
                self.unlabelled.append(code)
                continue
            cx, cy = polygon_centroid(region)
            ax.text(cx, cy, text, ha="center", va="center", fontsize=stylia.FONTSIZE_SMALL,
                    zorder=4,
                    path_effects=[patheffects.withStroke(linewidth=1.2, foreground="white")])

        # The rectangle fills the panel, so the colour key goes below it (the tight crop on
        # save keeps it). No size key: in a treemap the area *is* the value.
        self.legend(SOURCE_TYPE_COLORS, loc="upper left", bbox_to_anchor=(0.0, -0.01),
                    ncol=len(SOURCE_TYPE_COLORS))
        self.label(title="Models targeting priority pathogens")


#: Footprint of each stacked subtask panel, keyed by field. Both are 60 mm wide and short. Source
#: Type's 3 bars at (1, 2) draw 4.9 mm bars; Output holds 4 bars, so at the same footprint its bars
#: would come out a quarter thinner. Note bar thickness is 0.8 x (axes height) / n_bars and the axes
#: is the figure MINUS a FIXED ~11.7 mm tick+xlabel band, so thickness does not scale with the
#: footprint: exactly matching 4.9 mm needs (1 - 11.66/30) x 4/3 + 11.66/30 = 1.204 cells (36.1 mm),
#: not the 4/3 x 30 = 40 mm a naive ratio suggests. 1.15 is a deliberate compromise chosen for layout
#: — taller than Source Type, shorter than a full thickness match — leaving Output's bars ~7% thinner
#: (4.53 mm vs 4.88 mm). Raise to 1.204 to match exactly. Fractional cells are the documented
#: off-grid exception.
_SUBTASK_STACK_CELLS = {
    "Source Type": (1, 2),
    "Output": (1.15, 2),
}


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
                          ``<field>_by_subtask``, sized from ``_SUBTASK_STACK_CELLS``.
    """
    top_n = top_n or {}
    cross_tabs = cross_tabs or {}

    # One entry per panel. Footprints (rows, cols) in 3 cm cells, sized for a 183 x 170 mm Nature
    # page: most panels are 2x2 (60 x 60 mm), the two ten-category fields are quarter-width squares,
    # the technical boxes are 45 mm squares, and the tag cloud and pathogen panels are 3x3 (90 mm) —
    # the cloud for legibility, the pathogen panels for their two legends and 15 genus labels.
    plots = [
        TaskSubtaskBarPlot(sub=counts["Subtask"], cells=(2, 2)),
        # Unit alternative to the bars, on the same SUBTASK_COLORS. Reports any unfilled trailing
        # cells, which the script prints.
        TaskSubtaskWafflePlot(sub=counts["Subtask"], cells=(2, 2)),
        # Quarter-width squares (45.75 mm = 183/4), so the two sit side by side in half a page row.
        # Every named area/organism is Annotation-only, so the bars carry the Annotation hue and the
        # "Any" catch-all the neutral one — see `catchall_colors` for the one exception.
        FieldBarPlot(counts=counts["Biomedical Area"], title="Biomedical Area",
                     n=top_n.get("Biomedical Area"), cells=QUARTER_SQUARE,
                     color_fn=catchall_colors),
        FieldBarPlot(counts=counts["Target Organism"], title="Target Organism",
                     n=top_n.get("Target Organism"), cells=QUARTER_SQUARE,
                     color_fn=catchall_colors),
        # Licence: 60 mm rather than quarter-width because "CC-BY-NC-ND-4.0" is a long tick label
        # and the reuse-class key needs somewhere to sit. The bottom five bars are 1-4 models, so the
        # legend goes in the empty lower right of the axes rather than below it.
        FieldBarPlot(counts=counts["License grouped"], title="License",
                     cells=(2, 2), color_fn=license_colors,
                     legend_map=LICENSE_CLASS_COLORS),
        # Pie counterpart to the bar above, at reuse-class granularity — see the class docstring for
        # why the ten individual licences cannot be a pie. Both pies sit in the 45 mm tier: with only
        # two to four slices they need no more room, and their labels sit outside the circle anyway.
        LicenseClassPiePlot(counts=counts["License grouped"], cells=SMALL_SQUARE),
        # 90 mm square, not 60: at 60 mm the ten single-model tags render at 2.9 pt. See the class
        # docstring — the panel prints its own measured minimum type size so this stays checkable.
        TagCloudPlot(counts=counts["Tag"], cells=(3, 3)),
        # Technical trio. Name and label both derive from RUNTIME_BATCH so they cannot drift from the
        # column actually plotted. The -1 sentinels are skipped, never imputed, and any task left with
        # nothing keeps a "not measured" slot (the class reports coverage, which the script prints).
        TaskMetricBoxPlot(df=df, column=RUNTIME_COLUMN, name=f"runtime_{RUNTIME_BATCH}",
                          ylabel=f"Runtime for {RUNTIME_BATCH:,} molecules (s)",
                          cells=SMALL_SQUARE),
        TaskMetricBoxPlot(df=df, column="Image Size", name="image_size",
                          ylabel="Docker image size (MB)", cells=SMALL_SQUARE),
        ArchitecturePiePlot(df=df, cells=SMALL_SQUARE),
        PathogenTreemapPlot(df=df, pathogens_path=pathogens_path,
                            training_sizes_path=training_sizes_path, cells=(3, 3)),
        # Space-filling alternative to the circle version, on the same n ** AREA_EXPONENT areas
        # so the two panels are directly comparable. Pass area_metric="linear" for areas strictly
        # proportional to compound count — faithful, but illegible with this data (see README).
        # Reports its own area-fit quality and legibility, which the script prints.
        PathogenVoronoiPlot(df=df, pathogens_path=pathogens_path,
                            training_sizes_path=training_sizes_path, cells=(3, 3),
                            area_metric="power"),
    ]

    # Stacked variants: one bar per field value, segmented by subtask, so a single panel carries the
    # joint distribution. No legend — `task_subtask` and the waffle are the shared subtask key.
    for field, table in cross_tabs.items():
        plots.append(StackedFieldBarPlot(
            table=table, colors=SUBTASK_COLORS, legend_kw=None,
            cells=_SUBTASK_STACK_CELLS[field],
            name=f"{field.lower().replace(' ', '_')}_by_subtask"))

    footprints = {}
    for p in plots:
        p.save(output_dir)
        footprints[p.name] = list(p.cells)
        if isinstance(p, TaskSubtaskWafflePlot):
            print(f"\n[{p.name}] unfilled trailing cells in the last row: {p.blank}")
        if isinstance(p, TaskMetricBoxPlot):
            gaps = {t: c for t, c in p.coverage.items() if c[0] < c[1]}
            print(f"\n[{p.name}] measured per task: "
                  + ", ".join(f"{t} {n}/{tot}" for t, (n, tot) in p.coverage.items()))
            if gaps:
                print(f"{'':{len(p.name) + 3}}unmeasured models skipped (Airtable -1 sentinel), "
                      f"NOT imputed: " + ", ".join(f"{t} {tot - n}" for t, (n, tot) in gaps.items()))
        if isinstance(p, ArchitecturePiePlot):
            print(f"\n[{p.name}] {p.counts} — snapshot, not a trend "
                  f"(45% dual-arch among 2021 models vs 77% among 2026)")
        if isinstance(p, TagCloudPlot):
            tag = f"[{p.name}]"
            print(f"\n{tag} {p.n_placed}/{p.n_tags} tags placed"
                  + ("" if p.n_placed == p.n_tags
                     else f"  <-- wordcloud DROPPED {p.n_tags - p.n_placed}"))
            print(f"{tag} smallest type {p.min_font_pt:.2f} pt at {p.cells[1] / CELLS_PER_WIDTH * 180:.0f}"
                  f" mm wide; below the {TagCloudPlot._PRINT_FLOOR_PT:g} pt print floor: {p.below_floor}")
            print(f"{tag} RASTER panel — the one exception to the vector-PDF rule. Read ranks, "
                  f"not ratios; exact counts are in tag_counts.csv")
        if isinstance(p, PathogenVoronoiPlot):
            tag = f"[{p.name}]"
            print(f"\n{tag} max relative area error: {p.area_error:.4f}")
            print(f"{tag} cells below {PathogenVoronoiPlot.MIN_READABLE_FRACTION:g} of the "
                  f"panel (unreadable at print size): {len(p.tiny)}")
            for code, eosid, frac in sorted(p.tiny, key=lambda t: t[2]):
                print(f"           {code:14s} {eosid}  {frac:.2e}")
            print(f"{tag} regions too small to label: {len(p.unlabelled)}"
                  + (f" ({', '.join(p.unlabelled)})" if p.unlabelled else ""))

    with open(os.path.join(output_dir, "figure_cells.json"), "w") as f:
        json.dump(footprints, f, indent=2)
