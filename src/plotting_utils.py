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
def hbar(ax, labels, values, *, colors=None, abbreviate=False, name_map=None,
         ref=None, xlim=None):
    """Plain horizontal bars, one per (label, value), drawn top-to-bottom in the given order.

    The FIRST element sits at the top (natural reading order), so callers pass data already
    ordered the way they want it displayed. ``colors`` is a single colour or a per-bar list.
    """
    y = np.arange(len(labels))
    ax.barh(y, values, color=colors)
    if abbreviate:
        abbrev_ticks(ax, labels, axis="y", name_map=name_map)
    else:
        ax.set_yticks(y)
        ax.set_yticklabels(list(labels))
    ax.invert_yaxis()  # first element on top
    if ref is not None:
        ref_line(ax, ref, axis="x")
    if xlim is not None:
        ax.set_xlim(*xlim)


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

    ``colors`` is a ``(part, rest)`` pair. ``s`` is the marker area in points squared.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    part, rest = colors
    ax.scatter(x, y, s=s, marker="o", facecolor=rest, edgecolors=edgecolor,
               linewidths=linewidth, zorder=zorder)
    for xi, yi, fi in zip(x, y, np.asarray(frac, dtype=float)):
        if not np.isfinite(fi) or fi <= 0:
            continue
        ax.scatter([xi], [yi], s=s, marker=_wedge_path(min(fi, 1.0)), facecolor=part,
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
                 xlim=(0, 1), counts=None, label_fn=None, abbreviate=True, name_map=None):
    """Grouped horizontal bars: one row per label, ``len(order)`` bars per row.

    ``series`` maps ``(label, group) -> value``; ``colors`` maps ``group -> colour``.
    ``counts`` (optional) maps ``(label, group) -> n`` and annotates each bar with its sample
    size at the bar's left end, so a high value on very few points is not read as solid.
    ``label_fn`` prettifies the legend group labels.
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
    swatch_legend(ax, {label_fn(g): colors[g] for g in order})
    st.label(ax, xlabel=xlabel, ylabel="", title=title)


# --------------------------------------------------------------------------- #
# Boxes with jittered points                                                   #
# --------------------------------------------------------------------------- #
def box_with_jitter(ax, values, position, color, *, face=None, vert=True, width=0.34,
                    filled=True, median_color=None, line_color=None, jitter=True,
                    jitter_width=0.12, cap=CAP_POINTS, point_size=6, point_alpha=0.5,
                    rng=None, showfliers=False, label=None):
    """One box (with optional jittered points) — the single house style for distribution boxes.

    House style: filled box tinted by ``face`` (defaults to ``color``), whiskers/caps/median
    in :data:`plotting_colors.INK`, jittered points in ``color`` (subsampled to ``cap``).
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
    """Apply the house box style (filled body, INK outlines/median) to a boxplot dict."""
    median_color = INK if median_color is None else median_color
    line_color = INK if line_color is None else line_color
    if label is not None:
        bp["boxes"][0].set_label(label)
    if filled:
        bp["boxes"][0].set_facecolor(face if face is not None else color)
        bp["boxes"][0].set_edgecolor(line_color)
    else:
        bp["boxes"][0].set_color(color)
    for element in ("whiskers", "caps"):
        for line in bp[element]:
            line.set_color(color if not filled else line_color)
    for line in bp["medians"]:
        line.set_color(median_color)


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
def swatch_legend(ax, mapping, *, loc="lower right", **kw):
    """A legend of colour swatches from a ``{label: colour}`` mapping (semi-transparent bg)."""
    handles = [Patch(color=c, label=l) for l, c in mapping.items()]
    return ax.legend(handles=handles, loc=loc, fontsize=st.FONTSIZE_SMALL,
                     **LEGEND_KW, **kw)


def nested_size_legend(ax, keys, areas, *, x, y_base, color=None, label_fmt="{:,}",
                       fontsize=None, linewidth=0.7, zorder=6):
    """Size key drawn as **nested circles sharing a bottom tangent**, one per key.

    Nested rather than a row or column of separate markers: when the keys span decades the largest is
    many times the smallest, so a row of to-scale markers takes over the panel while a stacked legend
    wastes the space between them. Sharing a tangent also lets a reader compare each ring directly
    against the one inside it.

    ``areas`` are the matching scatter ``s`` values (points squared) — pass the *same* function the
    data uses, so the key cannot drift from the marks. A scatter marker's path radius is
    ``sqrt(s) / 2`` points, which is what is converted here.

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
def specificity_bars(ax, df, *, title=""):
    """Diverging horizontal bars of a per-pathogen specificity index (positive vs negative)."""
    df = df.dropna(subset=["specificity_index"]).sort_values(
        "specificity_index", ascending=False)
    values = df["specificity_index"].values
    colors = [hue("turquoise") if v >= 0 else hue("crimson") for v in values]
    hbar(ax, df["pathogen"].tolist(), values, colors=colors, abbreviate=True, ref=0.0)
    st.label(ax, xlabel="Specificity index (same − mean cross-pathogen AUROC)",
             ylabel="", title=title)
