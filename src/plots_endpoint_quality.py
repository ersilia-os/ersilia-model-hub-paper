"""Figures for the per-endpoint agreement audit (step 09).

Five panels off the tables written by :mod:`eval_endpoint_quality` and
:func:`eval_predictor_performance.run_activity_self_performance`, all bioactivity-only (no property
data enters any of them):

  - :class:`EndpointUprankingPlot` — one box per pathogen over its endpoints' median AUROC against
    their own peers, every endpoint drawn as a point. The headline panel: it shows both which
    pathogens have a weak block and which individual endpoints sit below chance inside it.
  - :class:`EndpointSpecificityPlot` — same-pathogen against different-pathogen AUROC, one point per
    endpoint. Separates the two ways an endpoint can be uninformative: low on both axes is
    idiosyncratic, high on both is non-specific.
  - :func:`endpoint_ranked_tail_figure` — the weakest endpoints, named. Following the precedent of
    :func:`plots_matrix_analyses.pathogen_jaccard_figure`, this one is a **plain-matplotlib
    diagnostic rather than a paper panel**: it carries a label per endpoint, which goes illegible at
    the 180 mm print width. It still draws only through :mod:`plotting_colors` hues and
    :func:`plotting_utils.sentence_case`, under the stylia print/article style, and writes PNG
    **and** PDF.
  - :class:`ActivitySelfPerformancePlot` / :class:`PathogenSubsetPerformancePlot` — moved here from
    the property-predictor figures (former step 13) since their data, activity-endpoints-predicting-
    each-other, is bioactivity-only and now comes from step 09 itself.

Colour encodes ``organism_class`` in the first three, and is a **secondary** encoding everywhere: the
first two panels name every pathogen on the axis and the third names every endpoint. That choice is
forced — :func:`plotting_colors.distinct_colors` carries only 9 substantive hues and explicitly
warns against using colour as the sole encoding beyond that, and there are 12 pathogens in scope.
The 6 organism classes fit comfortably, and using one shared key across the three panels means a
reader learns it once.

Consensus endpoints are drawn as a distinct marker rather than a distinct colour. They are not fair
comparators — a ``consensus_score`` column is its model's aggregate over the very peers it is scored
against, so its agreement is inflated by construction — and the panels have to let a reader discount
them without hiding them.
"""

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import stylia as st
from matplotlib.lines import Line2D

import plotting_base  # noqa: F401  (applies the stylia print/article style on import)
from default import ORGANISM_CLASS_ORDER, PREDICTOR_CHANCE_LEVEL, PREDICTOR_METRICS, RANDOM_SEED
from plotting_base import BasePlot
from plotting_colors import INK, REFERENCE_LINE, distinct_colors, hue
from plotting_utils import (_jitter_points, abbrev, box_with_jitter, marker_legend,
                            merge_figure_cells, ref_line, save_diagnostic_figure, sentence_case)

#: Organism class -> colour, in the canonical class order so the key reads the same in every panel.
#: Six substantive categories, comfortably inside the 9-hue palette.
CLASS_COLORS = dict(zip(ORGANISM_CLASS_ORDER, distinct_colors(len(ORGANISM_CLASS_ORDER))))

#: Marker per endpoint kind. A consensus column aggregates the peers it is scored against, so it is
#: shown but must be visually separable from the assay-derived endpoints around it.
CONSENSUS_MARKER = "D"
ENDPOINT_MARKER = "o"

POINT_SIZE = 9
CONSENSUS_POINT_SIZE = 14
POINT_ALPHA = 0.8
BOX_WIDTH = 0.6
JITTER_WIDTH = 0.18


def _class_color(organism_class):
    """Colour for an organism class, falling back to the neutral hue for an unlisted one."""
    return CLASS_COLORS.get(organism_class, REFERENCE_LINE)


def _class_legend(ax, table, **kw):
    """Key for the organism classes present, plus the consensus marker when one is drawn.

    The consensus entry is conditional rather than always present: a key must not describe a mark
    the panel does not contain, and the ranked-tail diagnostic legitimately shows no consensus
    column at all (they rank near the top of their pathogens, not in the weak tail).
    """
    entries = [{"label": c, "color": CLASS_COLORS[c], "marker": ENDPOINT_MARKER, "markersize": 3}
               for c in ORGANISM_CLASS_ORDER if (table["organism_class"] == c).any()]
    if table["is_consensus"].any():
        entries.append({"label": "consensus column", "color": INK,
                        "marker": CONSENSUS_MARKER, "markersize": 3})
    return marker_legend(ax, entries, **kw)


class EndpointUprankingPlot(BasePlot):
    """One box per pathogen over its endpoints' median same-pathogen AUROC; one point per endpoint.

    The y value of a point is that endpoint's MEDIAN AUROC against the other endpoints of its own
    pathogen — "how well does this endpoint rank the compounds its siblings call active". The dashed
    line is chance. A point below it is an endpoint that ranks its peers' actives below its own
    inactives, which is a statement about direction, not a threshold.

    Boxes are ordered by median so the panel reads as a ranking of the pathogens' blocks, not as
    config order. Each pathogen's endpoint count is on the axis, because a 6-endpoint box and a
    64-endpoint box are not equally informative and the box alone does not say which is which.
    """

    #: 12 pathogens over 4/6 of the print width leaves ~10 mm per box, which holds a rotated genus
    #: label and a 64-point swarm without either crowding the other.
    def __init__(self, table, ax=None, cells=(2, 4), seed=RANDOM_SEED):
        super().__init__(ax=ax, cells=cells)
        self.name = "09_endpoint_upranking"

        block = table.dropna(subset=["auroc_out_same_median"])
        if not len(block):
            self._unavailable()
            return

        rng = np.random.default_rng(seed)  # jitter is stochastic; seeded so the figure is stable
        order = (block.groupby("pathogen")["auroc_out_same_median"].median()
                 .sort_values(ascending=False).index.tolist())

        labels = []
        for i, pathogen in enumerate(order):
            g = block[block["pathogen"] == pathogen]
            color = _class_color(g["organism_class"].iloc[0])
            # Box over ALL the pathogen's endpoints; the points carry the individual endpoints, and
            # the consensus ones are re-drawn on top so they can be discounted by eye.
            box_with_jitter(self.ax, g["auroc_out_same_median"], i, color, width=BOX_WIDTH,
                            filled=False, jitter_width=JITTER_WIDTH, point_size=POINT_SIZE,
                            point_alpha=POINT_ALPHA, rng=rng)
            cons = g[g["is_consensus"]]
            if len(cons):
                self.ax.scatter(np.full(len(cons), i), cons["auroc_out_same_median"],
                                s=CONSENSUS_POINT_SIZE, marker=CONSENSUS_MARKER, facecolor="none",
                                edgecolors=INK, linewidths=0.6, zorder=4)
            labels.append(f"{abbrev(g['organism'].iloc[0])}\nn={len(g)}")

        self.ax.set_xticks(range(len(order)))
        self.ax.set_xticklabels(labels, rotation=90, ha="center", fontsize=st.FONTSIZE_SMALL)
        self.ax.set_xlim(-0.8, len(order) - 0.2)
        self.ref_line(PREDICTOR_CHANCE_LEVEL, axis="y")
        self.label(ylabel="median AUROC against own-pathogen endpoints")
        _class_legend(self.ax, block, loc="lower left")


class EndpointSpecificityPlot(BasePlot):
    """Same-pathogen against different-pathogen median AUROC, one point per endpoint.

    The identity diagonal is the reference that matters, not either axis alone. An endpoint far to
    the RIGHT of it ranks its own pathogen's actives better than other pathogens' — the behaviour
    the endpoint exists to have. An endpoint ON the diagonal ranks every pathogen's actives equally
    well, so whatever it is measuring is not specific to its own organism. An endpoint low on BOTH
    axes ranks nothing well.

    The two chance lines are drawn as well, since the bottom-left quadrant they cut off — below
    chance in both directions — is where an inverted or degenerate endpoint lands.
    """

    #: Square footprint: both axes are the same quantity on the same scale, so an unequal aspect
    #: would tilt the identity diagonal and make the comparison the panel exists for harder to read.
    def __init__(self, table, ax=None, cells=(3, 3)):
        super().__init__(ax=ax, cells=cells)
        self.name = "09_endpoint_specificity"

        block = table.dropna(subset=["auroc_out_same_median", "auroc_out_diff_median"])
        if not len(block):
            self._unavailable()
            return

        lo = float(min(block["auroc_out_same_median"].min(),
                       block["auroc_out_diff_median"].min())) - 0.03
        hi = float(max(block["auroc_out_same_median"].max(),
                       block["auroc_out_diff_median"].max())) + 0.03
        self.ax.plot([lo, hi], [lo, hi], linestyle="--", linewidth=0.8, color=REFERENCE_LINE,
                     zorder=1)
        self.ref_line(PREDICTOR_CHANCE_LEVEL, axis="x")
        self.ref_line(PREDICTOR_CHANCE_LEVEL, axis="y")

        for consensus, marker, size in ((False, ENDPOINT_MARKER, POINT_SIZE),
                                        (True, CONSENSUS_MARKER, CONSENSUS_POINT_SIZE)):
            g = block[block["is_consensus"] == consensus]
            if not len(g):
                continue
            self.ax.scatter(g["auroc_out_same_median"], g["auroc_out_diff_median"],
                            s=size, marker=marker, alpha=POINT_ALPHA, zorder=3,
                            color=[_class_color(c) for c in g["organism_class"]],
                            edgecolors=INK if consensus else "none",
                            linewidths=0.5 if consensus else 0)

        self.ax.set_xlim(lo, hi)
        self.ax.set_ylim(lo, hi)
        self.ax.set_aspect("equal", adjustable="box")
        self.label(xlabel="median AUROC, own pathogen",
                   ylabel="median AUROC, other pathogens")
        _class_legend(self.ax, block, loc="upper left")


#: Endpoints named in the ranked-tail diagnostic. A DISPLAY limit only — every endpoint is in
#: ``09_endpoint_quality.csv``, nothing is filtered out of any statistic, and the figure's own axis
#: label states how many of how many are shown. The full 249 would need a ~100-inch page at a
#: legible row height.
TAIL_N = 40


def endpoint_ranked_tail_figure(table, *, name, output_dir, n=TAIL_N, chance=PREDICTOR_CHANCE_LEVEL):
    """The weakest endpoints by median same-pathogen AUROC, one labelled row each.

    Drawn as a deficit from chance — a stem running from the chance line to the endpoint's value —
    rather than as a bar from zero, so the length of the mark is the quantity a reader cares about
    (how far short of chance it falls) instead of being dominated by the 0-to-0.5 span every
    endpoint shares.

    A diagnostic, not a paper panel: it carries one label per endpoint, which cannot stay legible at
    page width. Returns ``(png_path, pdf_path)``.
    """
    block = table.dropna(subset=["auroc_out_same_median"]).head(int(n))
    rows = len(block)

    fig, ax = plt.subplots(figsize=(7.2, max(3.5, rows * 0.16)))
    y = np.arange(rows)
    values = block["auroc_out_same_median"].to_numpy(dtype=float)
    colors = [_class_color(c) for c in block["organism_class"]]

    ax.hlines(y, chance, values, colors=colors, linewidth=1.0, zorder=2)
    for consensus, marker, size in ((False, ENDPOINT_MARKER, 14), (True, CONSENSUS_MARKER, 20)):
        m = (block["is_consensus"] == consensus).to_numpy()
        if not m.any():
            continue
        ax.scatter(values[m], y[m], s=size, marker=marker, zorder=3,
                   color=[c for c, keep in zip(colors, m) if keep],
                   edgecolors=INK if consensus else "none", linewidths=0.5 if consensus else 0)

    ref_line(ax, chance, axis="x")
    labels = [f"{abbrev(o)} · {c}  ({m}, {int(p)} peers)"
              for o, c, m, p in zip(block["organism"], block["column_name"],
                                    block["model_id"], block["n_peers"])]
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=5)
    ax.set_ylim(-0.8, rows - 0.2)
    ax.invert_yaxis()
    ax.set_xlabel(sentence_case(
        f"median AUROC against own-pathogen endpoints  —  weakest {rows} of "
        f"{len(table)} endpoints (full ranking in the CSV)"))
    # Upper right: the rows are sorted ascending, so the top ones are the shortest stems and leave
    # that corner empty, while the bottom rows run out to the right edge.
    handles = [Line2D([], [], linestyle="none", marker=ENDPOINT_MARKER, color=CLASS_COLORS[c],
                      markersize=4, label=c)
               for c in ORGANISM_CLASS_ORDER if (block["organism_class"] == c).any()]
    if block["is_consensus"].any():
        handles.append(Line2D([], [], linestyle="none", marker=CONSENSUS_MARKER,
                              markerfacecolor="none", markeredgecolor=INK, markersize=4,
                              label="consensus column"))
    ax.legend(handles=handles, loc="upper right", fontsize=5, frameon=True, facecolor="white",
              framealpha=0.85)
    fig.tight_layout()
    return save_diagnostic_figure(fig, name, output_dir)


#: Pair kind -> colour for the activity-vs-activity figure. Two endpoints of the SAME pathogen
#: should agree; separating them is what keeps "this model is self-consistent" from being read as
#: "this model predicts everything".
PAIR_COLORS = {True: hue("crimson"), False: hue("cobalt")}
PAIR_LABELS = {True: "same pathogen", False: "cross pathogen"}

#: Points drawn per pathogen box before subsampling, before the 2026-09-02 restriction to the 15
#: pathogens of interest this reached ~16,500 for P. falciparum's box alone (64 endpoints x 306
#: full-library targets) — far past what a 3 mm box can show. The consensus-collapsed, 15-pathogen
#: population both panels now share is far smaller, so this cap is a safety margin rather than a
#: binding constraint in practice, but it costs nothing to keep.
SELF_POINT_CAP = 400

#: Scoped to the two performance panels below, distinct from EndpointUprankingPlot/
#: EndpointSpecificityPlot's own POINT_SIZE/POINT_ALPHA/JITTER_WIDTH (tuned for a dozen boxes of a
#: few endpoints each): these two panels pool hundreds to thousands of pairs per box, so their
#: points must be smaller and more transparent to stay legible. Unchanged from the former step-13
#: figures these were moved from.
PERF_POINT_SIZE = 1.5
PERF_POINT_ALPHA = 0.25
PERF_JITTER_WIDTH = 0.22


class ActivitySelfPerformancePlot(BasePlot):
    """Activity endpoints predicting each other, grouped by the predictor endpoint's pathogen,
    restricted to the 15 pathogens of interest on BOTH sides (2026-09-02).

    One box per pathogen over its (consensus-collapsed) endpoints' AUROCs against every OTHER
    endpoint's binarized version, with points coloured by whether the target endpoint belongs to
    the same pathogen. Reads the SAME ``09_pathogen_subset_self_performance.csv`` as
    :class:`PathogenSubsetPerformancePlot` below — both restrict predictor AND target to the 15
    pathogens, with each ChEMBL model collapsed to its consensus column
    (:func:`eval_predictor_performance.pathogen_subset_endpoints`), so "cross-pathogen" here means
    "one of the other 14 priority pathogens", not any of the other 41 organisms in the full
    selection. That match matters: without it, one pathogen's box could be dominated by a single
    model's dozens of correlated sub-assays, exactly the count-domination problem the consensus
    collapse exists to avoid — see :class:`PathogenSubsetPerformancePlot`'s docstring.

    Self-pairs are excluded: an endpoint against its own binarization is 1.0 by construction and
    would add one guaranteed-perfect point to every box.
    """

    def __init__(self, subset, ax=None, cells=(2, 6), seed=RANDOM_SEED):
        super().__init__(ax=ax, cells=cells)
        self.name = "09_performance_activity_by_organism"

        block = subset[~subset["self_pair"]]
        if not len(block):
            self._unavailable()
            return

        rng = np.random.default_rng(seed)  # jitter + point subsampling are stochastic
        order = (block.groupby("predictor_pathogen")["value"].median()
                 .sort_values(ascending=False).index.tolist())

        for i, pathogen in enumerate(order):
            g = block[block["predictor_pathogen"] == pathogen]
            # Neutral box over ALL of this pathogen's pairs; the two colours ride on the points, so
            # the box must not claim either of them.
            box_with_jitter(self.ax, g["value"], i, INK, width=BOX_WIDTH, filled=False,
                            jitter=False)
            for same in (False, True):
                vals = g.loc[g["same_organism"] == same, "value"].to_numpy(dtype=float)
                vals = vals[np.isfinite(vals)]
                if not len(vals):
                    continue
                _jitter_points(self.ax, vals, i, PAIR_COLORS[same], vert=True,
                               jitter_width=PERF_JITTER_WIDTH, cap=SELF_POINT_CAP,
                               point_size=PERF_POINT_SIZE, point_alpha=PERF_POINT_ALPHA, rng=rng)

        self.ax.set_xticks(range(len(order)))
        self.ax.set_xticklabels([abbrev(o) for o in order], rotation=90, ha="center")
        self.ax.set_xlim(-0.8, len(order) - 0.2)
        self.ref_line(PREDICTOR_CHANCE_LEVEL, axis="y")
        self.label(ylabel=f"{PREDICTOR_METRICS['continuous'].upper()}")
        self.legend({PAIR_LABELS[k]: PAIR_COLORS[k] for k in (True, False)})


class PathogenSubsetPerformancePlot(BasePlot):
    """One box per activity endpoint of the 15 pathogens of interest, grouped by pathogen.

    The consensus collapse (see ``eval_predictor_performance.pathogen_subset_endpoints``) cuts the
    endpoint count far enough that each one gets its own box, which the pooled-by-pathogen
    :class:`ActivitySelfPerformancePlot` above cannot show: a pathogen's consensus score and its
    individual assay endpoints sit side by side.

    Boxes are ordered by pathogen (pathogens by median, endpoints by median within), and coloured by
    pathogen as a SECONDARY cue only — every box is identified on the axis as
    ``{pathogen} - {endpoint}``, so nothing depends on telling 15 hues apart, which
    ``plotting_colors.distinct_colors`` explicitly warns against as a sole encoding.

    Self-pairs are excluded for the same reason as the full figure: 1.0 by construction.
    """

    #: Taller than the other step-09 panels (4 rows, not 2): the tick labels carry both the pathogen
    #: and the endpoint, and at 59 categories they need roughly half the panel height on their own.
    #: At (2, 6) the boxes were squeezed into a thin strip along the top.
    def __init__(self, subset, ax=None, cells=(4, 6), seed=RANDOM_SEED):
        super().__init__(ax=ax, cells=cells)
        self.name = "09_performance_pathogen_subset"

        block = subset[~subset["self_pair"]]
        if not len(block):
            self._unavailable()
            return

        rng = np.random.default_rng(seed)
        by_pathogen = (block.groupby("predictor_pathogen")["value"].median()
                       .sort_values(ascending=False).index.tolist())
        by_endpoint = block.groupby("predictor_endpoint")["value"].median()

        order, colors = [], {}
        palette = distinct_colors(len(by_pathogen))
        for pathogen, color in zip(by_pathogen, palette):
            eps = block.loc[block["predictor_pathogen"] == pathogen,
                            "predictor_endpoint"].unique().tolist()
            for e in sorted(eps, key=lambda x: -by_endpoint[x]):
                order.append((pathogen, e))
                colors[e] = color

        for i, (_, endpoint) in enumerate(order):
            g = block[block["predictor_endpoint"] == endpoint]
            box_with_jitter(self.ax, g["value"], i, colors[endpoint], width=BOX_WIDTH,
                            filled=False, jitter=False)
            for same in (False, True):
                vals = g.loc[g["same_organism"] == same, "value"].to_numpy(dtype=float)
                vals = vals[np.isfinite(vals)]
                if not len(vals):
                    continue
                # Same/cross organism stays a shape distinction here, since colour is already
                # carrying pathogen identity: cross-organism points are the pale majority, and
                # same-organism ones are drawn darker and larger on top.
                _jitter_points(self.ax, vals, i, colors[endpoint], vert=True,
                               jitter_width=PERF_JITTER_WIDTH, cap=SELF_POINT_CAP,
                               point_size=PERF_POINT_SIZE * (3 if same else 1),
                               point_alpha=0.75 if same else PERF_POINT_ALPHA, rng=rng)

        labels = [f"{abbrev(p)} - {e.split(':')[-1]}" for p, e in order]
        self.ax.set_xticks(range(len(order)))
        self.ax.set_xticklabels(labels, rotation=90, ha="center")
        self.ax.set_xlim(-0.8, len(order) - 0.2)
        self.ref_line(PREDICTOR_CHANCE_LEVEL, axis="y")
        self.label(ylabel=f"{PREDICTOR_METRICS['continuous'].upper()}")


# --------------------------------------------------------------------------- #
# Entry point                                                                  #
# --------------------------------------------------------------------------- #
def save_activity_self_figure(output_dir, subset=None):
    """The activity-vs-activity figure, restricted to the 15 pathogens of interest on both sides
    (2026-09-02), merged into ``figure_cells.json``. Reads the SAME
    ``09_pathogen_subset_self_performance.csv`` as :func:`save_pathogen_subset_figure`.
    """
    if subset is None:
        subset = pd.read_csv(os.path.join(output_dir, "09_pathogen_subset_self_performance.csv"))
    plot = ActivitySelfPerformancePlot(subset)
    new = {}
    if plot.is_available:
        plot.save(output_dir)
        new[plot.name] = list(plot.cells)
        print(f"  figure: {plot.name}")
    else:
        print("  [skip figure] 09_performance_activity_by_organism: no rows")
    return merge_figure_cells(output_dir, new)


def save_pathogen_subset_figure(output_dir, subset=None):
    """The 15-pathogen, consensus-collapsed figure, merged into ``figure_cells.json``."""
    if subset is None:
        subset = pd.read_csv(os.path.join(output_dir, "09_pathogen_subset_self_performance.csv"))
    plot = PathogenSubsetPerformancePlot(subset)
    new = {}
    if plot.is_available:
        plot.save(output_dir)
        new[plot.name] = list(plot.cells)
        print(f"  figure: {plot.name}")
    else:
        print("  [skip figure] 09_performance_pathogen_subset: no rows")
    return merge_figure_cells(output_dir, new)


def save_endpoint_quality_figures(output_dir, table=None):
    """Build all three per-endpoint-quality panels and record the grid footprints in
    ``figure_cells.json``.

    Only the two :class:`plotting_base.BasePlot` panels have a cell footprint; the ranked-tail
    diagnostic is sized in inches by its row count, exactly as step 09's per-pathogen figure is, so
    it has no entry in the manifest. The two performance panels (:func:`save_activity_self_figure`,
    :func:`save_pathogen_subset_figure`) are built separately, since they read a different table.
    """
    if table is None:
        table = pd.read_csv(os.path.join(output_dir, "09_endpoint_quality.csv"))

    footprints = {}
    for plot in (EndpointUprankingPlot(table), EndpointSpecificityPlot(table)):
        if plot.is_available:
            plot.save(output_dir)
            footprints[plot.name] = list(plot.cells)
            print(f"  figure: {plot.name}")
        else:
            print(f"  [skip figure] {plot.name}: no rows")

    png_path, _ = endpoint_ranked_tail_figure(table, name="09_endpoint_ranked_tail",
                                              output_dir=output_dir)
    print(f"  figure: {os.path.basename(png_path)[:-4]} (diagnostic, no cell footprint)")

    merge_figure_cells(output_dir, footprints)
    return footprints
