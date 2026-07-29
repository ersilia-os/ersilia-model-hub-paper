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
import matplotlib.patches as mpatches
import matplotlib.patheffects as patheffects
import stylia

from plotting_base import BasePlot
from plotting_colors import TASK_COLORS, SOURCE_TYPE_COLORS, BAR_DEFAULT, output_colors, hue
from plotting_utils import abbrev, hbar
from voronoi_treemap import polygon_area, polygon_centroid, voronoi_treemap
from default import RANDOM_SEED


class FieldBarPlot(BasePlot):
    """Horizontal bar chart of a metadata field's value counts.

    Parameters
    ----------
    counts : DataFrame with columns ``value`` and ``count`` (already sorted).
    title  : panel title.
    colors : optional list of per-bar colours (same order as ``counts``); defaults to
             a single ``BAR_DEFAULT`` colour for every bar.
    n      : optional top-N cap; when set, only the first ``n`` rows are shown and the
             title gets a "(top n)" suffix.
    cells  : footprint on the reference grid as ``(rows, cols)`` — taller for panels with
             more bars (see ``save_metadata_figures``).
    """

    def __init__(self, ax=None, counts=None, title="", colors=None, n=None, cells=(3, 3)):
        BasePlot.__init__(self, ax=ax, cells=cells)
        self.name = title.lower().replace(" ", "_")
        if n:
            counts = counts.head(n)
            title = f"{title} (top {n})"
        hbar(self.ax, counts["value"].tolist(), counts["count"].tolist(),
             colors=colors if colors is not None else BAR_DEFAULT)
        self.label(xlabel="Number of models", ylabel=" ", title=title)


class TaskSubtaskBarPlot(BasePlot):
    """Combined Task + Subtask panel: one horizontal bar per subtask, coloured by its
    parent task, with a legend giving the per-task totals. This merges the two former
    panels into one — subtask breakdown from the bars, task totals from the legend.

    Parameters
    ----------
    sub         : subtask counts DataFrame with columns ``value``, ``count``, ``parent``
                  (already ordered by parent task, then count within task).
    task_counts : Task counts DataFrame (``value``, ``count``) for the legend totals.
    cells       : footprint on the reference grid as ``(rows, cols)``.
    """

    def __init__(self, ax=None, sub=None, task_counts=None, cells=(3, 3)):
        BasePlot.__init__(self, ax=ax, cells=cells)
        self.name = "task_subtask"
        colors = [TASK_COLORS[p] for p in sub["parent"]]
        hbar(self.ax, sub["value"].tolist(), sub["count"].tolist(), colors=colors)
        self.legend({f"{t} (n={n})": TASK_COLORS[t]
                     for t, n in zip(task_counts["value"], task_counts["count"])})
        self.label(xlabel="Number of models", ylabel=" ", title="Tasks & subtasks")


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
        self.fig.canvas.draw()
        per_data = _pts_per_data(ax, self.fig)

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


def save_metadata_figures(counts, df, pathogens_path, training_sizes_path, output_dir,
                          top_n=None):
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
    """
    top_n = top_n or {}

    # One entry per panel. Footprints (rows, cols) in 3 cm cells: taller for panels with more
    # bars, square for the treemap. Tune here to change a panel's physical size.
    plots = [
        TaskSubtaskBarPlot(sub=counts["Subtask"], task_counts=counts["Task"], cells=(3, 3)),
        FieldBarPlot(counts=counts["Source Type"], title="Source Type", cells=(2, 3),
                     colors=[SOURCE_TYPE_COLORS[v] for v in counts["Source Type"]["value"]]),
        FieldBarPlot(counts=counts["Output"], title="Output", cells=(3, 3),
                     colors=output_colors(len(counts["Output"]))),
        FieldBarPlot(counts=counts["Biomedical Area"], title="Biomedical Area",
                     n=top_n.get("Biomedical Area"), cells=(4, 3)),
        FieldBarPlot(counts=counts["Target Organism"], title="Target Organism",
                     n=top_n.get("Target Organism"), cells=(4, 3)),
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

    footprints = {}
    for p in plots:
        p.save(output_dir)
        footprints[p.name] = list(p.cells)
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
