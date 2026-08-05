"""Shared, ax-based plotting primitives — the reusable toolbox every panel builds on.

The concrete panels in the ``plots_*`` modules are thin: they prepare data and then call
these primitives, so the *house style* (colours, reference lines, box/jitter, ROC panels,
legends, tick handling) is defined ONCE here and every figure inherits it. Nothing in this
module instantiates its own palette — all colours flow through :mod:`plotting_colors`
(``hue`` / ``REFERENCE_LINE`` / ``INK``), never ``stylia.NamedColors`` or literal hex/grey.

Each primitive takes a matplotlib ``ax`` so it works both for a standalone :class:`BasePlot`
panel and for one axis of a multi-panel figure.
"""

import matplotlib.colors as mcolors
import matplotlib.patheffects as patheffects
import numpy as np
import stylia as st
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
from matplotlib.path import Path

from plotting_colors import INK, REFERENCE_LINE, hue

# Max individual points drawn in a jittered box overlay (dense distributions are subsampled).
CAP_POINTS = 1000

# House-style legend frame: a semi-transparent white background so legends stay readable when
# they sit over the plotted data. Every legend (primitive or hand-built) uses these kwargs.
LEGEND_KW = dict(frameon=True, facecolor="white", framealpha=0.7, edgecolor="none")


# --------------------------------------------------------------------------- #
# Text                                                                         #
# --------------------------------------------------------------------------- #
def abbrev(name, name_map=None):
    """Abbreviate a genus, optionally after mapping a code to a name first.

    ``abbrev("Mycobacterium tuberculosis")`` → ``"M. tuberculosis"``.
    ``abbrev("mtuberculosis", {"mtuberculosis": "Mycobacterium tuberculosis"})`` → same.
    A single-word or unmapped value is returned unchanged. This is the ONE abbreviation
    helper (it absorbed the former ``plots_chembl_performance._pathogen_label``).
    """
    if name_map:
        name = name_map.get(name, name)
    parts = str(name).split()
    return f"{parts[0][0]}. {' '.join(parts[1:])}" if len(parts) > 1 else name


# --------------------------------------------------------------------------- #
# Reference lines & tick labels                                                #
# --------------------------------------------------------------------------- #
def sentence_case(text):
    """Capitalise an axis label's first letter, leaving the rest of the string untouched.

    House style: an axis label reads as a sentence — "Own-assay AUROC", not "own-assay AUROC".
    Applied centrally (``plotting_base`` wraps ``stylia.label`` with it), so a panel may write its
    label in whatever case reads best in code and the figures stay consistent regardless.

    Two guards keep it from damaging scientific text, which is why this is a function and not a
    ``.capitalize()`` call:

    - Only the FIRST character changes, and only if it is a lowercase letter. ``str.capitalize``
      would lowercase everything after it and turn "own-assay AUROC" into "Own-assay auroc".
    - The whole string is left alone when its first WORD already carries an internal capital. That
      protects lowercase-leading tokens whose case is meaningful: "pIC50 (nM)", "cLogP",
      "mRNA count" must never become "PIC50", "CLogP", "MRNA count".

    Strings starting with punctuation or a digit are unchanged, so "-log10(FDR-adjusted q-value)",
    "|Spearman rho|" and "12 output columns" pass through. ``None`` and "" pass through too, since
    ``stylia.label`` treats ``None`` as "leave this axis alone".
    """
    if not text:
        return text
    words = text.split()
    head = words[0] if words else ""
    if any(ch.isupper() for ch in head[1:]):
        return text
    return text[0].upper() + text[1:] if text[0].islower() else text


def ref_line(ax, value, axis="x", *, linestyle="--", linewidth=0.8, color=None, **kw):
    """One dashed reference line (chance level, threshold, baseline) in the neutral colour."""
    color = REFERENCE_LINE if color is None else color
    draw = ax.axvline if axis == "x" else ax.axhline
    return draw(value, color=color, linestyle=linestyle, linewidth=linewidth, **kw)


def abbrev_ticks(ax, labels, axis="y", *, rotation=None, name_map=None, fontsize=None):
    """Set categorical tick labels, abbreviating genus names through :func:`abbrev`."""
    labs = [abbrev(l, name_map) for l in labels]
    ticks = range(len(labels))
    if axis == "y":
        ax.set_yticks(ticks)
        ax.set_yticklabels(labs, fontsize=fontsize)
    else:
        ax.set_xticks(ticks)
        ha = "right" if rotation not in (None, 0, 90) else "center"
        ax.set_xticklabels(labs, rotation=rotation or 0, ha=ha, fontsize=fontsize)


# --------------------------------------------------------------------------- #
# Bars                                                                         #
# --------------------------------------------------------------------------- #
#: Colour of a bar label according to what it sits on: **white** where it is over its bar, the default
#: text colour where it is over the page. No backing box and no stroked outline — both were tried and
#: both are worse at this type size: a chip reads as a sticker pasted over the mark, and a stroke thick
#: enough to separate dark text from saturated crimson floods the counters of O, e and A so the glyphs
#: stop reading as type. Reversing the text out of the bar needs no extra ink at all.
LABEL_ON_BAR_COLOR = "white"


def place_inside_labels(ax, texts, values, *, pad_frac=0.025, gap_frac=0.03):
    """Position and colour bar labels adaptively. Call AFTER ``tight_layout``.

    Two placements, chosen per bar by whether the label actually fits, and the colour follows the
    placement so no backing box is ever needed:

    * **fits inside** → centred on its bar, reversed out in ``LABEL_ON_BAR_COLOR`` (white);
    * **does not fit** → immediately after the bar's end, in the default text colour (black), which is
      what it is sitting on: the page.

    Placement is adaptive because centring unconditionally destroys the mark it labels — at these group
    sizes *Antiviral* is 6.3 mm of text on a 2.4 mm bar, so a centred label swamps the bar entirely. A
    bar is never obscured by its own label here, and no label ever straddles the bar's edge, which is
    what makes a single colour per label sufficient.

    Must run after ``tight_layout``: it weighs a label width fixed in *points* against an axes width
    layout is still free to change. A label with room neither inside its bar nor after it is tucked
    inside the right spine and stays **black**, since in that case the bar is short and the space it is
    tucked into is page, not bar.

    Returns ``(n_inside, n_after)``.
    """
    fig = ax.figure
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    inv = ax.transAxes.inverted()
    lo, hi = ax.get_xlim()
    span = (hi - lo) or 1.0
    n_inside = n_after = 0
    for t, v in zip(texts, values):
        bb = inv.transform(t.get_window_extent(r))
        w = bb[1][0] - bb[0][0]
        end = (float(v) - lo) / span          # bar's end, in axes fraction
        if w + 2 * pad_frac <= end:
            t.set_ha("center")
            t.set_x(end / 2.0)
            t.set_color(LABEL_ON_BAR_COLOR)   # reversed out of the bar it sits on
            n_inside += 1
        elif end + gap_frac + w <= 1.0 - pad_frac:
            t.set_ha("left")
            t.set_x(end + gap_frac)
            n_after += 1
        else:
            t.set_ha("right")                 # no room after the bar; tuck it inside the right spine
            t.set_x(1.0 - pad_frac)           # short bar, so this lands on the page: stays black
            n_after += 1
    return n_inside, n_after


def hbar(ax, labels, values, *, colors=None, abbreviate=False, name_map=None,
         ref=None, xlim=None, inside_labels=False, label_color=None, label_size=None,
         bar_fraction=0.8):
    """Plain horizontal bars, one per (label, value), drawn top-to-bottom in the given order.

    The FIRST element sits at the top (natural reading order), so callers pass data already
    ordered the way they want it displayed. ``colors`` is a single colour or a per-bar list.

    ``inside_labels=True`` moves the category names **off the y axis and onto the bars**: the y ticks
    go away entirely and each name is drawn in front of its own bar. For a very narrow panel this is
    the difference between a chart and a column of text — a tick-label gutter is a fixed ~14 mm
    whatever the panel's width, which at 25 mm would leave almost nothing for the bars.

    A label's colour follows what it sits on: **white** where it is centred on its bar, the default
    black where it is set after a bar too short to hold it (see :func:`place_inside_labels`). No
    backing box and no outline — reversing the text out of the bar costs no extra ink, and it lets the
    bar stay full strength rather than being washed out to keep dark text legible over it.

    The caller **must** run :func:`place_inside_labels` after ``tight_layout`` — that is what decides
    per bar whether the label is centred on it or set just after it, and it cannot be decided at draw
    time because the label width is fixed in points while the axes width is not. What is drawn here is
    only a provisional centred position. Returns the list of label ``Text`` artists.
    """
    y = np.arange(len(labels))
    ax.barh(y, values, color=colors, height=bar_fraction)
    texts = []
    if inside_labels:
        ax.set_yticks([])
        # Axes-fraction x, data y — so a label's horizontal position can be clamped against the axes
        # box later without having to know the count axis' limits.
        lo, hi = ax.get_xlim()
        span = (hi - lo) or 1.0
        for yi, lab, v in zip(y, labels, values):
            texts.append(ax.text((float(v) / 2.0 - lo) / span, yi, str(lab),
                                 transform=ax.get_yaxis_transform(),
                                 ha="center", va="center", zorder=5,
                                 color=label_color,      # None -> rcParams text.color (black)
                                 fontsize=label_size or st.FONTSIZE_SMALL))
    elif abbreviate:
        abbrev_ticks(ax, labels, axis="y", name_map=name_map)
    else:
        ax.set_yticks(y)
        ax.set_yticklabels(list(labels))
    ax.invert_yaxis()  # first element on top
    if ref is not None:
        ref_line(ax, ref, axis="x")
    if xlim is not None:
        ax.set_xlim(*xlim)
    return texts


def _wedge_path(frac, *, n=64):
    """Unit wedge from 12 o'clock clockwise over ``frac`` of a turn, as a marker Path."""
    a = np.pi / 2 - np.linspace(0.0, 2 * np.pi * float(frac), n)
    verts = np.concatenate([[[0.0, 0.0]], np.column_stack([np.cos(a), np.sin(a)]), [[0.0, 0.0]]])
    return Path(verts, closed=True)


def pie_scatter(ax, x, y, frac, colors, *, s=26, edgecolor="white", linewidth=0.3, zorder=3):
    """Scatter where every point is a two-slice pie: ``frac`` of it in the first colour.

    Drawn as matplotlib **marker paths**, not patches, which is what keeps every pie the same
    physical size and perfectly round on an axis with a log scale, a categorical x and no equal
    aspect — a ``Wedge`` patch would have to be sized in data units and would come out elliptical and
    varying. Marker paths are sized in points via ``s`` and are immune to the axis transform.
    Matplotlib scales a custom marker path by a flat 0.5 rather than normalising it to its bounding
    box, so a 10 % wedge and a 90 % wedge share one radius (checked: ``max|coord| == 0.5`` for every
    fraction) — without that, slice size would leak into apparent circle size.

    The "rest" slice is one full-circle scatter call for all points, with the ``frac`` wedges drawn
    over it, so N pies cost N+1 draw calls rather than 2N.

    **Both slices are stroked** in ``edgecolor``. The wedge sits on top of the circle, so leaving it
    unstroked punches a gap in the circle's own outline wherever the wedge meets the rim — which is
    exactly where two overlapping pies need an edge to stay separable. Stroking the wedge also draws
    its two radii, giving the usual pie slice separators.

    ``colors`` is a ``(part, rest)`` pair. ``s`` is the marker area in points squared — a scalar
    for uniform pies, or an array to size each pie individually (the wedge always takes ``s`` from
    its own point, so slice and radius stay independent encodings).
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    sizes = np.broadcast_to(np.asarray(s, dtype=float), x.shape)
    part, rest = colors
    ax.scatter(x, y, s=sizes, marker="o", facecolor=rest, edgecolors=edgecolor,
               linewidths=linewidth, zorder=zorder)
    for xi, yi, fi, si in zip(x, y, np.asarray(frac, dtype=float), sizes):
        if not np.isfinite(fi) or fi <= 0:
            continue
        # A full pie is drawn as a plain circle, not a 360-degree wedge: the wedge path closes
        # through the centre, so its two radii coincide at 12 o'clock and stroke a seam across an
        # otherwise solid disc.
        marker = "o" if fi >= 1.0 else _wedge_path(fi)
        ax.scatter([xi], [yi], s=si, marker=marker, facecolor=part,
                   edgecolors=edgecolor, linewidths=linewidth, zorder=zorder + 1)


def stacked_hbar(ax, labels, fracs, seg_order, seg_colors, *, xlim=(0, 1),
                 always_show=(), default_color=None):
    """Horizontal stacked bars over ``seg_order`` — one stack per label, first label on top.

    ``fracs`` is a list (one dict per label) of ``{segment: value}``. Segments whose values
    are all zero are skipped unless named in ``always_show``. ``seg_colors`` is a dict (or any
    object with ``.get``); missing segments fall back to ``default_color`` (neutral).
    """
    default_color = REFERENCE_LINE if default_color is None else default_color
    y = np.arange(len(labels))
    lefts = np.zeros(len(labels))
    for seg in seg_order:
        vals = np.array([f.get(seg, 0.0) for f in fracs], dtype=float)
        if seg not in always_show and vals.sum() == 0:
            continue
        ax.barh(y, vals, left=lefts, color=seg_colors.get(seg, default_color), label=seg)
        lefts += vals
    ax.set_yticks(y)
    ax.set_yticklabels(list(labels))
    ax.invert_yaxis()
    if xlim is not None:
        ax.set_xlim(*xlim)


def grouped_hbar(ax, y_labels, series, colors, order, *, xlabel="", title="", ref=0.5,
                 xlim=(0, 1), counts=None, label_fn=None, abbreviate=True, name_map=None,
                 legend=True):
    """Grouped horizontal bars: one row per label, ``len(order)`` bars per row.

    ``series`` maps ``(label, group) -> value``; ``colors`` maps ``group -> colour``.
    ``counts`` (optional) maps ``(label, group) -> n`` and annotates each bar with its sample
    size at the bar's left end, so a high value on very few points is not read as solid.
    ``label_fn`` prettifies the legend group labels. Pass ``legend=False`` on a panel too small to
    hold one — the key then has to be supplied separately, e.g. as a standalone key panel.
    """
    label_fn = label_fn or (lambda g: g)
    n = len(y_labels)
    ng = max(len(order), 1)
    h = 0.8 / ng
    ys = np.arange(n)
    for gi, g in enumerate(order):
        offset = (gi - (ng - 1) / 2) * h
        vals = [series.get((lbl, g), np.nan) for lbl in y_labels]
        ax.barh(ys + offset, vals, height=h, color=colors[g], label=label_fn(g))
        if counts is not None:
            for yi, (lbl, v) in enumerate(zip(y_labels, vals)):
                c = counts.get((lbl, g))
                if c is not None and np.isfinite(v):
                    ax.text(0.015, ys[yi] + offset, f"n={int(c)}", ha="left", va="center",
                            fontsize=5, color="white")
    if ref is not None:
        ref_line(ax, ref, axis="x")
    if abbreviate:
        abbrev_ticks(ax, y_labels, axis="y", name_map=name_map)
    else:
        ax.set_yticks(ys)
        ax.set_yticklabels(list(y_labels))
    if xlim is not None:
        ax.set_xlim(*xlim)
    if legend:
        swatch_legend(ax, {label_fn(g): colors[g] for g in order})
    st.label(ax, xlabel=xlabel, ylabel="", title=title)


# --------------------------------------------------------------------------- #
# Boxes with jittered points                                                   #
# --------------------------------------------------------------------------- #
# House line weight for distribution boxes. matplotlib's boxplot defaults are 1.0 pt — TWICE stylia's
# ``lines.linewidth`` of 0.5 and twice the axes spines — which made every box the heaviest mark in its
# panel and buried the swarm it is drawn over. Boxes now sit at the house weight, with the median one
# step up so it still reads as the summary statistic rather than as another whisker.
BOX_LINEWIDTH = 0.5
MEDIAN_LINEWIDTH = 1.4


def box_with_jitter(ax, values, position, color, *, face=None, vert=True, width=0.34,
                    filled=True, median_color=None, line_color=None, jitter=True,
                    jitter_width=0.12, cap=CAP_POINTS, point_size=6, point_alpha=0.5,
                    rng=None, showfliers=False, label=None):
    """One box (with optional jittered points) — the single house style for distribution boxes.

    House style: an **unfilled** box outlined in ``color`` at the house line weight, so the jittered
    swarm underneath stays the thing you read; median in the same colour but drawn heavier; points in
    ``color``, subsampled to ``cap``. Pass ``face`` explicitly for a tinted body — panels that draw
    marks *inside* the box need one for contrast — and the median
    then reverts to :data:`plotting_colors.INK` so it survives on an opaque body.
    Set ``filled=False`` for an outline-only box in ``color``, or ``jitter=False`` to omit the
    point overlay. Returns the matplotlib boxplot dict.
    """
    vals = np.asarray([v for v in np.asarray(values, dtype=float) if np.isfinite(v)],
                      dtype=float)
    if not len(vals):
        return None
    bp = ax.boxplot([vals], positions=[position], widths=width, vert=vert,
                    showfliers=showfliers, patch_artist=filled, manage_ticks=False)
    _style_box(bp, color, face=face, filled=filled, line_color=line_color,
               median_color=median_color, label=label)
    if jitter:
        _jitter_points(ax, vals, position, color, vert=vert, jitter_width=jitter_width,
                       cap=cap, point_size=point_size, point_alpha=point_alpha, rng=rng)
    return bp


def box_from_stats(ax, stats, position, color, *, face=None, vert=True, width=0.34,
                   filled=True, median_color=None, line_color=None, points=None,
                   jitter_width=0.12, cap=CAP_POINTS, point_size=6, point_alpha=0.5,
                   rng=None, label=None):
    """A house-style box drawn from PRECOMPUTED statistics instead of raw values.

    Same look as :func:`box_with_jitter` (that function's twin), for panels whose distribution
    is too large to ship per-molecule: ``stats`` is a mapping with ``median``, ``q1``, ``q3``,
    ``whisker_lo``, ``whisker_hi`` (the summary CSV's columns). Pass ``points`` to overlay
    jittered individual values where they *are* available and few enough to be readable.
    Returns the matplotlib boxplot dict, or None if the stats are incomplete.
    """
    keys = ("median", "q1", "q3", "whisker_lo", "whisker_hi")
    if any(stats.get(k) is None or not np.isfinite(float(stats.get(k, np.nan))) for k in keys):
        return None
    bp = ax.bxp([{
        "med": float(stats["median"]), "q1": float(stats["q1"]), "q3": float(stats["q3"]),
        "whislo": float(stats["whisker_lo"]), "whishi": float(stats["whisker_hi"]),
        "fliers": [],
    }], positions=[position], widths=width, vert=vert, showfliers=False,
        patch_artist=filled, manage_ticks=False)
    _style_box(bp, color, face=face, filled=filled, line_color=line_color,
               median_color=median_color, label=label)
    if points is not None and len(points):
        _jitter_points(ax, np.asarray(points, dtype=float), position, color, vert=vert,
                       jitter_width=jitter_width, cap=cap, point_size=point_size,
                       point_alpha=point_alpha, rng=rng)
    return bp


def _style_box(bp, color, *, face, filled, line_color, median_color, label=None):
    """Apply the house box style to a boxplot dict: transparent body, ``color`` outline and median.

    The outline (box edge, whiskers, caps) takes the **category colour**, not INK, so a box belongs
    visibly to its series. The median follows suit **on an unfilled box** and is simply drawn heavier
    (``MEDIAN_LINEWIDTH``) so it still reads as the summary statistic — on a 3 mm box a dark median
    was the heaviest mark in the panel and looked like a defect.

    On a **filled** box it falls back to INK, which is the only colour guaranteed to read on an
    opaque body: several panels pass ``face=color``, where a same-hue median would vanish outright.
    Override either with ``line_color`` / ``median_color``.

    ``face=None`` leaves the body transparent (the default). Line weights come from
    ``BOX_LINEWIDTH`` / ``MEDIAN_LINEWIDTH`` rather than matplotlib's heavier boxplot defaults.
    """
    line_color = color if line_color is None else line_color
    if median_color is None:
        median_color = INK if (filled and face is not None) else line_color
    if label is not None:
        bp["boxes"][0].set_label(label)
    if filled:
        # "none" (not None) is matplotlib's transparent facecolor; None would mean "the default".
        bp["boxes"][0].set_facecolor(face if face is not None else "none")
        bp["boxes"][0].set_edgecolor(line_color)
    else:
        bp["boxes"][0].set_color(line_color)
    bp["boxes"][0].set_linewidth(BOX_LINEWIDTH)
    for element in ("whiskers", "caps"):
        for line in bp[element]:
            line.set_color(line_color)
            line.set_linewidth(BOX_LINEWIDTH)
    for line in bp["medians"]:
        line.set_color(median_color)
        line.set_linewidth(MEDIAN_LINEWIDTH)


def _jitter_points(ax, vals, position, color, *, vert, jitter_width, cap, point_size,
                   point_alpha, rng):
    """Scatter ``vals`` around ``position``, subsampled to ``cap`` when an rng is given."""
    pts = vals if len(vals) <= cap or rng is None else rng.choice(vals, cap, replace=False)
    jit = rng.uniform(-jitter_width, jitter_width, len(pts)) if rng is not None \
        else np.zeros(len(pts))
    along = position + jit
    xs, ys = (along, pts) if vert else (pts, along)
    ax.scatter(xs, ys, s=point_size, alpha=point_alpha, color=color,
               edgecolors="none", zorder=3)


# --------------------------------------------------------------------------- #
# ROC panels                                                                   #
# --------------------------------------------------------------------------- #
def roc_panel(ax, fpr, tpr, auroc, n_pos, n_neg, color, *, xlabel="", ylabel="", title=""):
    """One square ROC panel with the chance diagonal, shaded curve and an AUC/count label.

    Pass ``fpr=None`` for the degenerate single-class case (draws a "no positives" note).
    Shared by every ROC grid so all ROC panels read identically.
    """
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xticks([0, 0.5, 1])
    ax.set_yticks([0, 0.5, 1])
    ax.set_aspect("equal", adjustable="box")
    if fpr is None:
        ax.text(0.5, 0.5, "no positives", ha="center", va="center",
                fontsize=st.FONTSIZE_SMALL, color=REFERENCE_LINE)
    else:
        fpr = np.asarray(fpr, dtype=float)
        tpr = np.asarray(tpr, dtype=float)
        ax.plot([0, 1], [0, 1], "--", color=REFERENCE_LINE, linewidth=0.5)
        ax.fill_between(fpr, tpr, alpha=0.15, color=color, linewidth=0)
        ax.plot(fpr, tpr, color=color)
        ax.text(0.96, 0.06, f"AUC {auroc:.2f}\n{n_pos}+ / {n_neg}−",
                ha="right", va="bottom", fontsize=st.FONTSIZE_SMALL, linespacing=1.1)
    st.label(ax, xlabel=xlabel, ylabel=ylabel, title=title)
    ax.title.set_fontsize(st.FONTSIZE_SMALL)


# --------------------------------------------------------------------------- #
# Heatmap                                                                      #
# --------------------------------------------------------------------------- #
def diverging_cmap(low="crimson", mid="white", high="cobalt"):
    """A crimson→white→cobalt diverging colormap anchored to ArticleColors hues."""
    lo = "white" if low == "white" else hue(low)
    hi = "white" if high == "white" else hue(high)
    md = "white" if mid == "white" else hue(mid)
    return mcolors.LinearSegmentedColormap.from_list("diverging", [lo, md, hi])


def sequential_cmap(name="cobalt"):
    """A white→hue single-hue sequential colormap (for 0→1 quantities with no chance-centre)."""
    return mcolors.LinearSegmentedColormap.from_list("sequential", ["white", hue(name)])


def heatmap(ax, matrix, *, cmap, norm, annotate=True, value_fmt="{:.2f}",
            text_light_when=None, nan_color=None, highlight=None, highlight_color=None,
            x_rotation=45, row_labels=None, col_labels=None, colorbar=False,
            annot_fontsize=None):
    """Annotated heatmap of a DataFrame: shared by the AUROC and correlation matrices.

    ``cmap``/``norm`` set the colour mapping; NaN cells render in ``nan_color`` (neutral).
    ``text_light_when(v) -> bool`` picks white vs INK text per cell; ``highlight`` outlines the
    listed ``(row, col)`` cells (e.g. a matrix diagonal).
    """
    nan_color = REFERENCE_LINE if nan_color is None else nan_color
    highlight_color = INK if highlight_color is None else highlight_color
    data = matrix.values.astype(float)
    nrows, ncols = data.shape
    cmap = cmap.copy()
    cmap.set_bad(nan_color)
    im = ax.imshow(np.ma.masked_invalid(data), cmap=cmap, norm=norm, aspect="auto")
    if colorbar:
        ax.figure.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    if highlight:
        for (ri, ci) in highlight:
            ax.add_patch(Rectangle((ci - 0.5, ri - 0.5), 1, 1, fill=False,
                                   edgecolor=highlight_color, linewidth=2))
    if annotate:
        for i in range(nrows):
            for j in range(ncols):
                v = data[i, j]
                if np.isfinite(v):
                    light = bool(text_light_when(v)) if text_light_when else False
                    ax.text(j, i, value_fmt.format(v), ha="center", va="center",
                            fontsize=annot_fontsize or st.FONTSIZE_SMALL,
                            color="white" if light else INK)
    rl = list(matrix.index) if row_labels is None else row_labels
    cl = list(matrix.columns) if col_labels is None else col_labels
    ax.set_xticks(range(ncols))
    ha = "right" if x_rotation not in (None, 0, 90) else "center"
    ax.set_xticklabels(cl, rotation=x_rotation or 0, ha=ha)
    ax.set_yticks(range(nrows))
    ax.set_yticklabels(rl)
    return im


# --------------------------------------------------------------------------- #
# Legends                                                                      #
# --------------------------------------------------------------------------- #
def swatch_legend(ax, mapping, *, loc="lower right", hatches=None, **kw):
    """A legend of colour swatches from a ``{label: colour}`` mapping (semi-transparent bg).

    ``hatches`` optionally maps the same labels to matplotlib hatch strings, for a panel that
    distinguishes categories by fill pattern rather than by hue. The swatch then carries the pattern
    as well as the colour, so the key cannot describe a mark the panel does not draw.
    """
    hatches = hatches or {}
    handles = [Patch(facecolor=c, edgecolor="white", hatch=hatches.get(l) or None, label=l)
               for l, c in mapping.items()]
    return ax.legend(handles=handles, loc=loc, fontsize=st.FONTSIZE_SMALL,
                     **LEGEND_KW, **kw)


def nested_size_legend(ax, keys, areas, *, x, y_base, color=None, label_fmt="{:,}", title=None,
                       fontsize=None, linewidth=0.7, zorder=6):
    """Size key drawn as **nested circles sharing a bottom tangent**, one per key.

    Nested rather than a row or column of separate markers: when the keys span decades the largest is
    many times the smallest, so a row of to-scale markers takes over the panel while a stacked legend
    wastes the space between them. Sharing a tangent also lets a reader compare each ring directly
    against the one inside it.

    ``areas`` are the matching scatter ``s`` values (points squared) — pass the *same* function the
    data uses, so the key cannot drift from the marks. A scatter marker's path radius is
    ``sqrt(s) / 2`` points, which is what is converted here.

    Pass ``title`` rather than adding your own text: the label rows are spread at draw time, so only
    this function knows where the block actually ends, and a caller guessing an offset will collide
    with the lowest label.

    Circles are drawn as ``Ellipse`` patches with the point radius converted separately per axis, so
    they render round whatever the axes aspect or scale. Requires the axes limits to be final: the
    conversion reads the live transform, so call this last.
    """
    from matplotlib.patches import Ellipse

    color = REFERENCE_LINE if color is None else color
    fig = ax.figure
    fig.canvas.draw()
    bb = ax.get_window_extent()
    w_pts = bb.width * 72.0 / fig.dpi
    h_pts = bb.height * 72.0 / fig.dpi
    (x0, x1), (y0, y1) = ax.get_xlim(), ax.get_ylim()
    per_pt_x = (x1 - x0) / w_pts
    per_pt_y = (y1 - y0) / h_pts

    fontsize = fontsize or st.FONTSIZE_SMALL
    radii = [float(np.sqrt(a)) / 2.0 for a in areas]
    ordered = sorted(zip(keys, radii), key=lambda t: -t[1])   # largest ring first (drawn behind)
    rmax = max(radii)

    # Label rows must be de-collided. Nesting packs the ring tops within 2*rmax of each other, which
    # for a decade key is a few millimetres, while each label needs its own line height — so anchoring
    # a label at its own ring's top guarantees overlap. Spread them evenly instead and reach back to
    # each ring with a leader.
    line_h = 1.25 * fontsize * per_pt_y
    tops = [y_base + 2 * r * per_pt_y for _k, r in ordered]
    span = max(tops[0] - tops[-1], line_h * (len(ordered) - 1))
    centre = (tops[0] + tops[-1]) / 2
    label_ys = [centre + span / 2 - i * span / max(len(ordered) - 1, 1)
                for i in range(len(ordered))]

    elbow = x + rmax * per_pt_x * 1.4
    for (key, r_pts), top, ly in zip(ordered, tops, label_ys):
        rx, ry = r_pts * per_pt_x, r_pts * per_pt_y
        ax.add_patch(Ellipse((x, y_base + ry), width=2 * rx, height=2 * ry, facecolor="none",
                             edgecolor=color, linewidth=linewidth, zorder=zorder))
        ax.plot([x, elbow, elbow + rmax * per_pt_x * 0.4], [top, ly, ly],
                color=color, linewidth=0.5, solid_joinstyle="round", zorder=zorder)
        ax.text(elbow + rmax * per_pt_x * 0.7, ly, label_fmt.format(key), ha="left", va="center",
                fontsize=fontsize, zorder=zorder)

    if title:
        ax.text(x, min(label_ys) - line_h * 1.15, title, ha="center", va="top",
                fontsize=fontsize, zorder=zorder)


def marker_legend(ax, entries, *, loc="lower right", **kw):
    """A legend of point markers (semi-transparent white background).

    ``entries`` is a list of dicts: ``{"label", "color", "marker"='o', "markersize"=None,
    "linestyle"='none'}``.
    """
    handles = [
        Line2D([], [], linestyle=e.get("linestyle", "none"), marker=e.get("marker", "o"),
               markerfacecolor=e["color"], markeredgecolor="none",
               markersize=e["markersize"] if e.get("markersize") else 6,
               color=e["color"], label=e["label"])
        for e in entries
    ]
    return ax.legend(handles=handles, loc=loc, fontsize=st.FONTSIZE_SMALL,
                     **LEGEND_KW, **kw)


# --------------------------------------------------------------------------- #
# Specialised bar helper still called by a single panel                        #
# --------------------------------------------------------------------------- #
def specificity_bars(ax, df, *, title="", order=None):
    """Diverging horizontal bars of a per-pathogen specificity index (positive vs negative).

    ``order`` is an explicit pathogen sequence, FIRST ELEMENT AT THE TOP (this draws via :func:`hbar`,
    which inverts the y axis). Pass it to tie the rows to another panel's order; omit it to sort by
    the index itself, strongest first. Pathogens absent from ``order`` are dropped, so build it from
    the frame's own values.
    """
    df = df.dropna(subset=["specificity_index"])
    if order:
        df = df.set_index("pathogen").reindex([p for p in order if p in set(df["pathogen"])]) \
               .reset_index()
    else:
        df = df.sort_values("specificity_index", ascending=False)
    values = df["specificity_index"].values
    colors = [hue("turquoise") if v >= 0 else hue("crimson") for v in values]
    hbar(ax, df["pathogen"].tolist(), values, colors=colors, abbreviate=True, ref=0.0)
    st.label(ax, xlabel="Specificity index (same − mean cross-pathogen AUROC)",
             ylabel="", title=title)
