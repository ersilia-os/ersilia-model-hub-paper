"""Concrete plot classes and the figure entry point for script 01 (model metadata).

Mirrors ``zairachem/report/plots.py`` (one class per plot type, all subclassing
:class:`plotting_base.BasePlot`). ``save_metadata_figures`` renders each panel as its own
Nature-sized figure and saves it as both PNG and vector PDF (see the function docstring).
"""

import json
import os

import circlify
import numpy as np
import matplotlib.patches as mpatches
import stylia

from plotting_base import BasePlot
from plotting_colors import TASK_COLORS, SOURCE_TYPE_COLORS, BAR_DEFAULT, output_colors
from plotting_utils import abbrev


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
        c = colors if colors is not None else [BAR_DEFAULT] * len(counts)
        self.ax.barh(counts["value"][::-1], counts["count"][::-1], color=c[::-1])
        stylia.label(self.ax, xlabel="Number of models", ylabel=" ", title=title)


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
        a = self.ax
        colors = [TASK_COLORS[p] for p in sub["parent"]]
        a.barh(sub["value"][::-1], sub["count"][::-1], color=colors[::-1])
        handles = [mpatches.Patch(color=TASK_COLORS[t], label=f"{t} (n={n})")
                   for t, n in zip(task_counts["value"], task_counts["count"])]
        a.legend(handles=handles, loc="lower right")
        stylia.label(a, xlabel="Number of models", ylabel=" ", title="Tasks & subtasks")


def _packed_dots(cx, cy, r, n, gap=0.94):
    """Lay out ``n`` non-overlapping dots that fill a circle of radius ``r``.

    Points follow a golden-angle (sunflower) spiral; the dot radius is set to half the
    minimum pairwise spacing of that layout (times ``gap`` for a hairline separation), so
    dots pack the available surface without colliding. The layout is scaled so the outermost
    dot just touches the circle rim. Returns ``(xs, ys, dot_r)``.
    """
    if n == 1:
        return [cx], [cy], 0.9 * r
    golden = np.pi * (3 - np.sqrt(5))
    idx = np.arange(n)
    rad = np.sqrt((idx + 0.5) / n)          # unit-disk radii (<= 1)
    ang = idx * golden
    ux, uy = rad * np.cos(ang), rad * np.sin(ang)
    pts = np.column_stack([ux, uy])
    d = np.sqrt(((pts[:, None, :] - pts[None, :, :]) ** 2).sum(-1))
    np.fill_diagonal(d, np.inf)
    rho = 0.5 * d.min()                      # touching dot radius in unit space
    r_out = rad.max()
    scale = r / (r_out + rho)                # so the outermost dot just reaches the rim
    dot_r = gap * rho * scale
    return list(cx + ux * scale), list(cy + uy * scale), dot_r


class PathogenTreemapPlot(BasePlot):
    """Circle-treemap of models per priority pathogen.

    One circle per pathogen (area ~ model count), packed with ``circlify``. Inside each
    circle, one dot per model laid out on a golden-angle spiral. Every model dot is the same
    absolute size across pathogens (because circle area ~ n, so r/sqrt(n) is constant), and
    dot colour encodes Source Type.

    Parameters
    ----------
    df             : the (Status == "Ready") metadata DataFrame.
    pathogens_path : CSV with a ``pathogen`` column (priority pathogens).
    cells          : footprint on the reference grid as ``(rows, cols)`` — square by default.
    """

    def __init__(self, ax=None, df=None, pathogens_path=None, cells=(3, 3)):
        BasePlot.__init__(self, ax=ax, cells=cells)
        self.name = "pathogen_circles"
        import pandas as pd

        ax = self.ax
        pathogens = pd.read_csv(pathogens_path)

        # Match each pathogen to models via the multi-value Target Organism field.
        df2 = df.copy()
        df2["orgs"] = df2["Target Organism"].fillna("").str.split(",").apply(
            lambda xs: [x.strip() for x in xs])

        rows = []
        for _, p in pathogens.iterrows():
            mask = df2["orgs"].apply(
                lambda orgs: any(p["pathogen"].lower() in o.lower() for o in orgs))
            rows.append((p["pathogen"], df2[mask].copy()))

        # Sort descending by count; minimum area 0.3 keeps empty circles visible.
        rows.sort(key=lambda x: len(x[1]), reverse=True)
        values = [max(len(models), 0.3) for _, models in rows]

        circles = circlify.circlify(values, show_enclosure=False)
        circles_sorted = sorted(circles, key=lambda c: c.r, reverse=True)

        # One global dot radius so EVERY model dot is the same absolute size. Take the smallest
        # packed radius across multi-model circles: drawing every dot at this size guarantees no
        # overlap anywhere, and single-model pathogens get a small circle sized to one such dot.
        ref = [_packed_dots(c.x, c.y, c.r, len(m))[2]
               for c, (_, m) in zip(circles_sorted, rows) if len(m) >= 2]
        d0 = min(ref) if ref else 0.12

        for circ, (name, models) in zip(circles_sorted, rows):
            cx, cy, r = circ.x, circ.y, circ.r
            n = len(models)
            colors = [SOURCE_TYPE_COLORS.get(v, BAR_DEFAULT)
                      for v in models["Source Type"].fillna("External")]

            if n == 0:
                r_draw = r
            elif n == 1:
                r_draw = d0 / 0.82  # shrink the circle to just hold its single equal-size dot
            else:
                xs_d, ys_d, _ = _packed_dots(cx, cy, r, n)
                r_draw = r

            facecolor = "#e8e8e8" if n == 0 else "#f0f0f0"
            ax.add_patch(mpatches.Circle((cx, cy), r_draw,
                                         facecolor=facecolor, edgecolor="white", zorder=1))
            if n == 1:
                ax.add_patch(mpatches.Circle((cx, cy), d0, facecolor=colors[0],
                                             edgecolor="white", linewidth=0.3, zorder=3))
            elif n > 1:
                for xd, yd, col in zip(xs_d, ys_d, colors):
                    ax.add_patch(mpatches.Circle((xd, yd), d0, facecolor=col,
                                                 edgecolor="white", linewidth=0.3, zorder=3))

            ax.text(cx, cy - r_draw - 0.03, abbrev(name),
                    ha="center", va="top", fontsize=stylia.FONTSIZE_SMALL, zorder=4)

        # Fit the axis tightly to the packed circles (small slack below for the genus labels)
        # so the treemap fills its panel.
        xmin = min(c.x - c.r for c in circles_sorted)
        xmax = max(c.x + c.r for c in circles_sorted)
        ymin = min(c.y - c.r for c in circles_sorted)
        ymax = max(c.y + c.r for c in circles_sorted)
        ax.set_xlim(xmin - 0.06, xmax + 0.06)
        ax.set_ylim(ymin - 0.13, ymax + 0.03)
        ax.set_aspect("equal", anchor="N")
        ax.set_axis_off()

        legend_handles = [mpatches.Patch(color=c, label=lbl)
                          for lbl, c in SOURCE_TYPE_COLORS.items()]
        ax.legend(handles=legend_handles, loc="lower right")
        stylia.label(ax, xlabel="", ylabel="", title="Models targeting priority pathogens")


def save_metadata_figures(counts, df, pathogens_path, output_dir, top_n=None):
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
    counts         : dict field -> counts DataFrame. ``counts["Subtask"]`` must also carry a
                     ``parent`` column (see the script's ordering step).
    df             : the (Status == "Ready") metadata DataFrame (for the treemap).
    pathogens_path : path to the priority-pathogens CSV.
    output_dir     : directory to write ``png/``, ``pdf/`` and ``figure_cells.json`` into.
    top_n          : optional dict field -> cap (e.g. {"Target Organism": 10}).
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
        PathogenTreemapPlot(df=df, pathogens_path=pathogens_path, cells=(3, 3)),
    ]

    footprints = {}
    for p in plots:
        p.save(output_dir)
        footprints[p.name] = list(p.cells)

    with open(os.path.join(output_dir, "figure_cells.json"), "w") as f:
        json.dump(footprints, f, indent=2)
