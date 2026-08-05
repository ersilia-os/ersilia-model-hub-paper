"""Step 05 figures — EU OpenScreen validation of the ChEMBL pathogen models.

Each panel is a standalone ``BasePlot`` (one chart per file, no A/B/C letters), sized on the
3 cm cell grid and saved as PNG + PDF with a ``figure_cells.json`` footprint manifest. Panels
read ONLY the small summary CSVs written by ``eval_euopenscreen.run_all`` — never per-molecule
data.

Every panel is expressed through the shared primitives in :mod:`plotting_utils`
(``grouped_hbar``, ``roc_panel`` via :class:`GridPlot`, ``heatmap``, ``specificity_bars``,
``marker_legend``) and pulls colours only from :mod:`plotting_colors`, so the look is identical
to the other steps' figures. The reusable :class:`MetricByOrganismPlot` base is also used by the
CoAdd figures in :mod:`plots_coadd`.
"""

import json
import math
import os

import numpy as np
import pandas as pd
from matplotlib.colors import Normalize, TwoSlopeNorm
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator

from plotting_base import BasePlot, GridPlot
from plotting_colors import (
    INK,
    REFERENCE_LINE,
    SHARED_ORGANISM_COLORS,
    count_shades,
    hue,
    pathogen_activity_colors,
)
import stylia

from plotting_utils import (
    LEGEND_KW,
    abbrev,
    abbrev_ticks,
    box_from_stats,
    diverging_cmap,
    grouped_hbar,
    heatmap,
    marker_legend,
    pie_scatter,
    roc_panel,
    sentence_case,
    sequential_cmap,
    specificity_bars,
)
from default import (
    ACTIVITY_CLASSES,
    DEDUP_SUBDIR,
    FULL_SUBDIR,
    HIT_CLASSES,
    RANDOM_SEED,
    SHARED_ORGANISMS,
)

# Footprint for the panels that go into the paper's two 45 mm rows. Three fit a 183 mm row with
# gutter; the sixth panel of a 6-panel set would need a second row (6 x 45 = 270 mm > 183 mm).
# At this size no in-axes legend fits, so these panels are built with `legend=False` and their keys
# are emitted as standalone `*_key` panels to place once in Illustrator.
SMALL_SQUARE = (1.5, 1.5)

# 40 mm square, for hit_exclusivity_events. Cells are 3 cm, so mm / 30 — 40 / 30 = 4/3. Smaller
# than SMALL_SQUARE because this panel spends ~13% of its width on the AUROC column to its right,
# so at 45 mm the plot plus its numbers overran the slot the AUROC panel occupies.
EVENT_SQUARE = (4 / 3, 4 / 3)

# The activity box panel: 40 mm tall, deliberately narrow. Two categories need no width, and a
# square panel spent most of it on white space between two boxes. 1.05 cols is the FLOOR while the
# x tick labels carry their class sizes — "(n = 97,162)" is 10.4 mm wide on its own, so the two
# labels collide below this (measured +1.5 mm clearance here, -0.3 mm at 0.9 cols). Dropping the
# counts into the caption, as step 01's TaskMetricBoxPlot does, would allow ~24 mm; they are kept
# on the panel instead because the 428-vs-97,162 imbalance is the context for everything in it.
BOX_NARROW = (4 / 3, 1.05)          # -> crop 31.3 x 41.1 mm

# ---------------------------------------------------------------------------------------------- #
# THE SIX MATCHED LANDSCAPE PANELS                                                                 #
# ---------------------------------------------------------------------------------------------- #
# One rule governs all six constants below: **every panel's AXES HEIGHT equals `euos_overlap`'s
# 25.81 mm.** These panels are laid out together in Illustrator, so what has to line up is the
# PLOTTING RECTANGLE, not the file. The crops are deliberately unequal (38.0 to 44.4 mm wide, 38.9 to
# 44.3 mm tall) because each panel spends a different amount of itself on chrome — a second x axis, an
# AUROC header, 45-degree class labels — and equalising the crops would push the data rectangles out
# of register, which is the opposite of what the matching is for. Do not "fix" the crop spread.
#
# `euos_overlap` is the reference because it is the constrained panel: it spends its height on chrome
# twice over (two x axes, each with ticks AND a label), so 30 mm of height left a 9.5 mm axes for 7
# genus rows that need 17.8 mm at the house 6 pt and the labels collided outright. It cannot give any
# height back, so everything else comes to it. Its width carries a deliberate +2 mm (`+ 2 / 30`, 1
# cell = 30 mm) over the tuned 45.0 mm crop, asked for so the panel breathes.
#
# WIDTHS are matched in pairs, not globally — panels sit side by side, so only their heights have to
# agree:
#     euos_overlap  <->  submodel_auroc_summary                     23.01 / 22.93 mm
#     hit_promiscuity  <->  exclusive_hit_model_rank_..._dedup      26.74 / 26.74 mm
#     hit_exclusivity_events  <->  consensus_sum_by_hit_class_dedup 24.05 / 24.00 mm
#
# ACHIEVED heights: 25.81 (overlap), 25.88 (promiscuity), 25.85 (rank), 25.71 (events), 25.87
# (submodel), 25.73 (consensus sum) — a 0.17 mm spread, i.e. +-0.09 mm about the reference. That is
# the floor, not a lack of effort: the tight bbox quantises to text extents so axes height moves in
# ~0.25 mm steps, and `hit_exclusivity_events` and `consensus_sum_..._dedup` were both searched at a
# 0.001-cell (0.03 mm) grid with no closer value on offer.
#
# HOW TO RE-TUNE, if any panel's chrome changes (a tick label, an axis label, a rotation, a legend):
# the relation is `axes = cells/6 * 180 mm - chrome` with chrome near-constant, so correct
# `cols += dw/30`, `rows += dh/30`, iterate 4-6 times, then local-search the last ~0.1 mm to beat the
# quantisation. Measure the AXES via `ax.get_position()` scaled by the figure size, NOT the crop —
# `crop_size_mm` predicts height exactly but runs ~1.2 mm under on width, and `pdf_page_mm` measures
# the page, not the plotting area. Do NOT shrink much further either way: below ~(0.6, 0.9)
# tight_layout gives up and the crop starts GROWING as the footprint falls (see the gotchas section of
# docs/figure_conventions.md).
WIDE_OVERLAP = (1.294, 1.395 + 2 / 30)   # reference -> crop 45.0 x 39.9, axes 23.01 x 25.81 mm
WIDE_PROMISCUITY = (1.2723, 1.3298)      # -> crop 41.1 x 39.4, axes 26.74 x 25.88 mm
WIDE_RANK = (1.2961, 1.2869)             # -> crop 39.9 x 40.3, axes 26.74 x 25.85 mm
WIDE_SUBMODEL = (1.2547, 1.4413)         # -> crop 44.4 x 38.9, axes 22.93 x 25.87 mm

# `hit_exclusivity_events` drops its own genus tick labels and is placed against `euos_overlap`'s, so
# its rows must land on that panel's rows: the axes height AND the y limits (-0.6, n - 0.4) both have
# to match, not just the footprint. Changing one without the other silently relabels every row.
WIDE_EVENTS = (1.3626, 1.3293)           # -> crop 41.1 x 42.2, axes 24.05 x 25.71 mm
BOX_EVENTS = (1.3205, 1.2116)            # -> crop 38.0 x 44.3, axes 24.00 x 25.73 mm

# SMALL_SQUARE plus 4 mm of width (1 cell = 30 mm), for `hit_exclusivity_auroc`. Its own constant
# rather than a wider SMALL_SQUARE: that footprint is shared by five other panels across steps 05 and
# 06 that are sized for the paper's 45 mm rows and must not move. Widening this one takes the paper's
# row B (`consensus_max_percentile_by_activity_dedup` + this + `submodel_auroc_summary`) from 137.3 to
# ~141.3 mm, still inside a 183 mm row.
EXCLUSIVITY_BARS = (1.5, 1.5 + 4 / 30)   # -> ~49 mm wide

# exclusive vs non-exclusive(shared) hits. Periwinkle for shared rather than the turquoise this pair
# used originally: turquoise is the repo's default/positive hue and was carrying FOUR unrelated
# meanings inside this one step — "shared hits" here, "novel to model" in the overlap panels (see
# LEAKAGE_COLORS), "EF @ 5%" in the enrichment panel and "narrow-spectrum" in the hit-class boxes. The
# first two moved to periwinkle; the last two still hold it, so it is down to one meaning per panel.
EXCLUSIVITY_COLORS = {"exclusive": hue("amber"), "nonexclusive": hue("periwinkle")}

# Training-set leakage, for the twin-axis overlap panels and their shared key. One definition because
# `EuosOverlapTwinPlot` and `euos_overlap_handles` each used to spell the pair out, and a key that
# disagrees with its bars is worse than no key.
LEAKAGE_COLORS = {"novel": hue("periwinkle"), "in_training": hue("crimson")}
# hit classes for the score-distribution box: the inactive background, organism-specific hits, and
# the shared hits split into narrow- (2-3 pathogens) and broad-spectrum (>3).
HIT_CLASS_COLORS = {
    # Silver, and unfilled: the inactive class is the BACKGROUND against which the hit classes
    # are read, not a finding of its own, and silver is the palette's reserved neutral. In cobalt
    # it carried the visual weight of a substantive category — on the activity panel it is 97,162
    # of 97,590 compounds, so a saturated box for it dominates a panel that exists to show the 428.
    "inactive": hue("silver"),
    "exclusive": hue("amber"),
    "narrow": hue("turquoise"),
    # Periwinkle, not crimson. Broad-spectrum IS the pan-active class, and periwinkle already carries
    # "shared / pan-active" in `hit_exclusivity_events` (see EXCLUSIVITY_COLORS) — the two panels are
    # read together, so the same idea must not change hue between them. Applied to BOTH the full/ and
    # deduplicated/ twins: a colour that means one thing in one variant and another in the other is a
    # caption hazard, since the panels differ only in leakage filtering.
    "broad": hue("periwinkle"),
    "active": hue("crimson"),   # the four active classes collapsed to one (activity view)
}
# enrichment factor top-fractions.
EF_COLORS = {"ef_1pct": hue("cobalt"), "ef_5pct": hue("turquoise")}

_PRETTY = {
    "raw": "raw", "dedup": "dedup (no leakage)",
    "exclusive": "exclusive", "nonexclusive": "shared",
    "inactive": "inactive", "narrow": "narrow\n(2-3)", "broad": "broad\n(>3)",
    "active": "active",
    "inhib_50": "single-point (inhib_50)", "mic_10": "MIC (mic_10)",
    "ef_1pct": "EF @ 1%", "ef_5pct": "EF @ 5%",
}


def _pretty(key):
    return _PRETTY.get(key, key)


def overlap_row_order(leak_df):
    """Organisms in ``euos_overlap`` row order, bottom row first: fewest EU OpenScreen actives up.

    THE single definition of that order, because two panels depend on it. ``euos_overlap`` carries
    the genus tick labels and ``hit_exclusivity_events`` has its own labels suppressed so the two
    can be placed side by side against one set of names — which is only honest if the rows are the
    same organisms in the same sequence. Sorting independently in each panel would work until the
    active counts shifted, and then silently relabel every row of one panel.

    Returns None when the leakage report is unusable, which the caller must treat as "cannot drop
    the labels".
    """
    if leak_df is None or leak_df.empty or "n_active" not in leak_df.columns:
        return None
    d = leak_df[leak_df["n_active"] > 0].sort_values("n_active", ascending=True)
    return d["pathogen"].tolist() or None


def _apply_row_order(present, order):
    """Reorder ``present`` pathogens into ``euos_overlap`` row order, bottom row first.

    Every panel with a pathogen y axis takes its row sequence from here so that the whole step reads
    as one figure: a reader who has learned the row order on one panel keeps it on all of them, and a
    row can be traced across panels without re-reading seven genus names each time. The cost is
    deliberate — these panels used to sort by their own metric (best on top), which is easier to rank
    within a single panel but makes any cross-panel comparison a lookup. The metric is still on the
    axis, so nothing is lost that the reader cannot see.

    Organisms unknown to ``order`` are APPENDED (they become the top rows) rather than dropped, so a
    panel can never silently lose a row to an incomplete leakage report. Returns None when there is no
    usable order, which every caller reads as "keep your own sort".

    NOTE the convention: bottom row FIRST, matching ``grouped_hbar`` and bare ``ax.barh``/``scatter``
    on ``arange(n)``. :func:`plotting_utils.hbar` is the exception — it calls ``invert_yaxis`` so its
    first element lands on TOP, and callers using it must reverse this list (see
    :class:`SpecificityIndexPlot`).
    """
    if not order:
        return None
    present = list(present)
    seen, known = set(present), set(order)
    return [p for p in order if p in seen] + [p for p in present if p not in known]


def _prefer(df, group_col):
    """Keep the dedup row per (identity, group) when present, else the raw row."""
    if df.empty or "set" not in df.columns:
        return df
    key_cols = [c for c in ("pathogen", "code", group_col) if c in df.columns]
    df = df.copy()
    df["_rank"] = (df["set"] == "dedup").astype(int)  # dedup preferred
    df = df.sort_values("_rank").drop_duplicates(key_cols, keep="last")
    return df.drop(columns="_rank")


# --------------------------------------------------------------------------- #
# Panels                                                                       #
# --------------------------------------------------------------------------- #
class MetricByOrganismPlot(BasePlot):
    """Grouped horizontal bars of one metric per organism, split by a grouping column.

    A reusable panel type: the EU OpenScreen hit-exclusivity figure and the CoAdd own-strain
    figure (:mod:`plots_coadd`) both build on it.
    """

    def __init__(self, df, group_col, group_colors, group_order, metric, xlabel, title,
                 name, ref=0.5, xlim=(0, 1), prefer_best=True, count_col=None,
                 cells=(3, 3), legend=True, ax=None, row_order=None):
        BasePlot.__init__(self, ax=ax, cells=cells)
        self.name = name
        if df is None or df.empty:
            self._unavailable()
            return
        d = _prefer(df, group_col) if prefer_best else df
        d = d[d[group_col].isin(group_order)]
        if d.empty:
            self._unavailable()
            return
        # Rows follow the shared `euos_overlap` order when the caller supplies one, so this panel can
        # be read against the others; otherwise fall back to this panel's own metric sort (ascending →
        # best on top). The fallback is what the CoAdd subclasses in `plots_coadd` use — they have a
        # different organism set and no EU OpenScreen leakage report to order by.
        y_labels = _apply_row_order(d["pathogen"].unique(), row_order)
        if y_labels is None:
            order_key = d.groupby("pathogen")[metric].max().sort_values(ascending=True)
            y_labels = order_key.index.tolist()
        series = {(r["pathogen"], r[group_col]): r[metric] for _, r in d.iterrows()}
        counts = None
        if count_col is not None and count_col in d.columns:
            counts = {(r["pathogen"], r[group_col]): r[count_col] for _, r in d.iterrows()}
        grouped_hbar(self.ax, y_labels, series, group_colors, group_order,
                     xlabel=xlabel, title=title, ref=ref, xlim=xlim, counts=counts,
                     label_fn=_pretty, legend=legend)


class EuosSharedEnrichmentPlot(BasePlot):
    """Analysis 1 (companion) — EU OpenScreen enrichment factor (top 1% / 5%) per organism,
    deduplicated. Enrichment = hit rate in the top k% over the base rate (1 = no enrichment)."""

    def __init__(self, own_df, ax=None, cells=(3, 3), row_order=None):
        BasePlot.__init__(self, ax=ax, cells=cells)
        self.name = "euos_shared_enrichment"
        if own_df is None or own_df.empty or "set" not in own_df.columns:
            self._unavailable()
            return
        # one row per organism, preferring the deduplicated set
        d = own_df.copy()
        d["_r"] = (d["set"] == "dedup").astype(int)
        d = d.sort_values("_r").drop_duplicates(["pathogen", "code"], keep="last")
        if d.empty:
            self._unavailable()
            return
        # Shared `euos_overlap` row order when supplied, else this panel's own EF@5% sort.
        y_labels = _apply_row_order(d["pathogen"], row_order)
        if y_labels is None:
            d = d.sort_values("ef_5pct", ascending=True)
            y_labels = d["pathogen"].tolist()
        series = {}
        for _, r in d.iterrows():
            for ef in ("ef_1pct", "ef_5pct"):
                series[(r["pathogen"], ef)] = r[ef]
        grouped_hbar(self.ax, y_labels, series, EF_COLORS, ["ef_1pct", "ef_5pct"],
                     xlabel="enrichment factor", title="", ref=1.0, xlim=None,
                     label_fn=_pretty)


class EuosRocGridPlot(GridPlot):
    """Analysis 1 (ROC view) — a grid of small ROC curves, one per organism (dedup), each shaded
    by its AUROC on a cobalt fading colormap (pale = near chance, saturated = strong ranking).
    Curves come from the precomputed ROC summary; panels are sorted best-AUROC first."""

    def __init__(self, roc_df, set_name="dedup", cols=3, name="euos_roc_grid"):
        items = self._items(roc_df, set_name)
        self.build_grid(items, cols=cols, name=name, panel_fn=self._panel,
                        color_fn=lambda it: it["color"],
                        edge_xlabel="FPR", edge_ylabel="TPR")

    @staticmethod
    def _items(roc_df, set_name):
        from plotting_colors import auroc_shades
        if roc_df is None or roc_df.empty or "set" not in roc_df.columns:
            return []
        d = roc_df[roc_df["set"] == set_name]
        if d.empty:
            return []
        aurocs = d.groupby("pathogen")["auroc"].first().sort_values(ascending=False)
        colors = auroc_shades(aurocs.values)
        items = []
        for org, color in zip(aurocs.index, colors):
            sub = d[d["pathogen"] == org].sort_values(["fpr", "tpr"])
            items.append(dict(
                fpr=sub["fpr"].values, tpr=sub["tpr"].values,
                auroc=float(sub["auroc"].iloc[0]),
                n_pos=int(sub["n_pos"].iloc[0]), n_neg=int(sub["n_neg"].iloc[0]),
                title=abbrev(org), color=color))
        return items

    @staticmethod
    def _panel(ax, item, color, xlabel, ylabel):
        roc_panel(ax, item["fpr"], item["tpr"], item["auroc"], item["n_pos"], item["n_neg"],
                  color, xlabel=xlabel, ylabel=ylabel, title=item["title"])


class EuosOverlapTwinPlot(BasePlot):
    """Training-set overlap for the whole EU OpenScreen library and for its actives in one panel,
    via twin x-axes (the two quantities differ ~300x). Per organism: an upper, hatched bar for
    the library (top axis, compounds) and a lower, solid bar for the actives (bottom axis), each
    stacked novel (periwinkle) vs in-training (crimson) — see :data:`LEAKAGE_COLORS`.

    **No gridlines on either axis.** With two x scales on one panel a vertical guide is ambiguous by
    construction — the reader cannot tell whether it belongs to the library scale or the actives one,
    and a grey vertical crossing a stacked bar additionally reads as a division in the data. Both
    grids are switched off explicitly rather than left to the rcParam default (``axes.grid = True``),
    so re-enabling one is a deliberate act.
    """

    #: Hatch for the library bars. Many short dashes rather than `//`'s few long ones: at this
    #: footprint a bar is ~1.5 mm tall, so a sparse hatch puts one or two strokes on a bar and reads
    #: as a defect rather than a texture. The stroke is thinned to match — at the default 1.0 pt a
    #: hatch this dense closes up into a solid block and the fill colour underneath is lost.
    HATCH = "//////"
    HATCH_WIDTH = 0.4
    #: Bar thickness and the pair's offset from the row centre, in y-data units (row pitch = 1). The
    #: pair therefore spans ``OFFSET +/- HEIGHT/2``, leaving ``2 * OFFSET - HEIGHT`` between the
    #: library bar and its actives bar and ``1 - 2 * OFFSET - HEIGHT`` between neighbouring organisms.
    #: At the original 0.38 the within-pair gap was 0.04 units = **0.15 mm** at this footprint, which
    #: prints as one merged block — and the library bars carry a white edge for the hatch, so the two
    #: read as a single hatched-then-solid shape. 0.34 gives 0.08 units (~0.30 mm) within the pair
    #: against 0.28 (~1.03 mm) between organisms, a ~3.5x ratio, so the pairing is unambiguous.
    BAR_HEIGHT = 0.34
    BAR_OFFSET = 0.21
    #: Tick-count target for the actives axis. The default locator offered only 0 and 200 across a
    #: 0-397 span, too coarse to read a bar against; 4 bins crowd once the genus labels have taken
    #: their ~14 mm out of a 40 mm panel, so the numbers run together.
    ACTIVES_TICK_BINS = 3
    # No tick-label size override here. The genus labels ran at 5 pt for a while, to buy height back
    # when the panel was fighting for it, which left them a size apart from the numbers on the other
    # axis — stylia puts every label (tick, axis, title, legend) at one size, and a panel with two
    # label sizes reads as a mistake. Moving the axis names inline freed enough room that the override
    # is not needed: the axes gives each of the 7 rows 3.7 mm against a 2.5 mm line height at 6 pt.
    #: Axis names, set inline on their own tick row (see :meth:`_inline_axis_labels`).
    BOTTOM_LABEL = "Actives"
    TOP_LABEL = "Full library"
    #: x of those names in axes coordinates — left of the axes, right-aligned, so they run back into
    #: the margin the genus tick labels already occupy rather than widening the crop. Far enough out
    #: to clear the "0" tick, which is CENTRED on the axes' left edge and so reaches into negative
    #: axes coordinates itself.
    LABEL_X = -0.10

    def __init__(self, leak_df, ax=None, cells=(3, 4), legend=True):
        BasePlot.__init__(self, ax=ax, cells=cells)
        self.name = "euos_overlap"
        need = {"n_active", "n_eval_conclusive", "n_overlap", "n_overlap_active"}
        if leak_df is None or leak_df.empty or not need.issubset(leak_df.columns):
            self._unavailable()
            return
        order = overlap_row_order(leak_df)
        d = (leak_df.set_index("pathogen").loc[order].reset_index()
             if order else leak_df.iloc[0:0])
        if d.empty:
            self._unavailable()
            return
        novel, in_training = LEAKAGE_COLORS["novel"], LEAKAGE_COLORS["in_training"]
        y = np.arange(len(d))
        h = self.BAR_HEIGHT
        off = self.BAR_OFFSET
        lib_novel = (d["n_eval_conclusive"] - d["n_overlap"]).values
        lib_ovl = d["n_overlap"].values
        act_novel = (d["n_active"] - d["n_overlap_active"]).values
        act_ovl = d["n_overlap_active"].values

        # top axis = library compounds; bottom axis (self.ax) = actives (twiny shares the y-axis)
        axt = self.ax.twiny()
        lib_bars = [
            axt.barh(y + off, lib_novel, height=h, color=novel, edgecolor="white",
                     hatch=self.HATCH),
            axt.barh(y + off, lib_ovl, left=lib_novel, height=h, color=in_training,
                     edgecolor="white", hatch=self.HATCH),
        ]
        for container in lib_bars:
            for patch in container:
                patch.set_hatch_linewidth(self.HATCH_WIDTH)
        self.ax.barh(y - off, act_novel, height=h, color=novel)
        self.ax.barh(y - off, act_ovl, left=act_novel, height=h, color=in_training)

        self.ax.set_ylim(-0.6, len(d) - 0.4)
        self.ax.set_yticks(y)
        self.ax.set_yticklabels([abbrev(p) for p in d["pathogen"]])
        self.ax.set_xlim(0, float((act_novel + act_ovl).max()) * 1.05)
        axt.set_xlim(0, float(d["n_eval_conclusive"].max()) * 1.05)
        self.ax.xaxis.set_major_locator(MaxNLocator(nbins=self.ACTIVES_TICK_BINS, integer=True))
        self.ax.grid(False)
        axt.grid(False)

        # Colour = novel/in-training, texture = which axis. The legend hangs outside the right
        # edge, which at 45 mm would take more width than the plot — hence the standalone
        # EuosOverlapKeyPlot.
        if legend:
            self.ax.legend(handles=euos_overlap_handles(), fontsize=5, loc="center left",
                           bbox_to_anchor=(1.01, 0.5), **LEGEND_KW)

        # Clears stylia's placeholder axis text. Without it the panel ships with a literal
        # "Y-axis / Units" down its left edge — the x names below are set on the axis objects
        # directly and never go through `label`, so nothing else overrides the default.
        self.label(xlabel="", ylabel="", title="")
        self._inline_axis_labels(axt)

    def _inline_axis_labels(self, axt):
        """Put each x axis's name at the LEFT END of its own tick-number row, not centred beyond it.

        A centred label costs a whole text row per axis, and this panel has two x axes — ~8 mm of the
        height went to naming them. Set on the tick row instead, the names occupy the empty corner to
        the left of the numbers, under the genus labels' column, which is space the panel was already
        spending. That is what buys the plotting area its extra height at an unchanged 40 mm crop.

        The row's y position is *measured*, not assumed: the label is pinned in axes coordinates to
        the mid-height of the drawn tick labels, so it stays on their line whatever the tick padding
        or font size. Requires layout to be final, hence the ``tight_layout`` + ``draw`` first.
        """
        import matplotlib.pyplot as plt

        plt.figure(self.fig.number)
        plt.tight_layout()
        self.fig.canvas.draw()
        renderer = self.fig.canvas.get_renderer()
        for ax, text in ((self.ax, self.BOTTOM_LABEL), (axt, self.TOP_LABEL)):
            boxes = [t.get_window_extent(renderer) for t in ax.get_xticklabels() if t.get_text()]
            ax.set_xlabel(sentence_case(text))
            if boxes:
                inv = ax.transAxes.inverted()
                y = float(np.mean([inv.transform((0, (b.y0 + b.y1) / 2))[1] for b in boxes]))
                ax.xaxis.set_label_coords(self.LABEL_X, y)
                ax.xaxis.label.set_va("center")
            ax.xaxis.label.set_ha("right")


class HitExclusivityEventPlot(BasePlot):
    """Where every hit lands in its own model's ranking — one tick per compound.

    The companion to :class:`HitExclusivityPlot`, showing the distribution the AUROC summarises
    instead of the summary. Two lanes per organism: shared hits (periwinkle, upper) and exclusive
    hits (amber, lower), each drawn as an event plot — one vertical line per hit, positioned at
    that compound's percentile in the model's ranking of the whole evaluated set. Hits massed at
    the right are upranked; a bar at 0.5 that looks respectable can turn out to be a few hits at
    the top over a flat spread, and only this panel shows the difference.

    **Lines, never a density.** Several subsets are tiny — *P. aeruginosa* exclusive is a single
    compound, *E. faecium* four — and a smoothed curve over one point draws a shape that is not in
    the data. One tick per compound cannot overstate the evidence: a one-hit lane is visibly one
    line. It also removes any bandwidth or bin-width choice.

    **The right-hand axis is the same quantity as the ticks.** AUROC is the mean rank percentile
    of the positives (exactly, ``(auroc * n_neg + (n_pos + 1) / 2) / n``; at this prevalence, the
    mean of the plotted ticks to four decimals), so the printed value is where each lane's centre
    of mass sits. It is a second reading of the lane, not a second variable —
    ``eval_euopenscreen._check_percentiles_match_auroc`` asserts the identity when the data is
    written.
    """

    #: Half-height of one lane, as a fraction of the row pitch, and the offset of each lane from
    #: the row centre. Two lanes per row with a hairline between them.
    _LANE_HALF = 0.19
    _LANE_OFFSET = 0.21

    #: Lanes are drawn shared-first so the row reads periwinkle-over-amber, matching the bar order
    #: of the AUROC panel it accompanies — and matching the ``shared/exclusive`` reading order of
    #: the paired AUROC values on the right.
    _LANE_ORDER = ("nonexclusive", "exclusive")

    #: Clearance in POINTS between the axes' right spine and the first digit of the AUROC column,
    #: and the gap either side of the slash. The column's distance from the axes is DERIVED from
    #: the rendered width of the widest left-hand value, not fixed: every point spent here is a
    #: point the tick strip does not get, and on a 40 mm panel whose axes is only ~8 mm wide that
    #: is the dominant layout cost. A hand-set offset left 11.4 pt (~4 mm) of dead space, because
    #: it had to be generous enough to clear a text width nobody had measured.
    #:
    #: Points, not axes fractions, so the clearance survives a change of footprint.
    #:
    #: The pair is drawn as three separate texts anchored right / centre / left rather than one
    #: string, which is what lets each value keep its lane's colour — a single ``"0.70/0.98"``
    #: string can only be one colour, and leaving the reader to infer which value is which from
    #: position is exactly the error the colours prevent.
    #: Line width of one hit tick. Thin: in the dense lanes (S. aureus exclusive, 176 hits in a
    #: ~12 mm strip) neighbouring ticks touch and the lane reads as a solid block, which hides the
    #: very clustering the panel exists to show.
    _TICK_WIDTH = 0.4

    #: An AUROC over this many hits or fewer is marked. Four of the fourteen qualify, and two of
    #: those are n = 1 and n = 4 — at that size a single compound moves the value by tenths, so an
    #: unmarked 0.98 beside a 0.86 computed over 176 compounds invites a comparison that the data
    #: does not support. The mark is a flag, not a filter: nothing is hidden or dropped.
    _SMALL_N = 10
    _SMALL_N_MARK = "*"

    _AUROC_CLEAR_PT = 3.0
    _AUROC_GAP_PT = 1.5
    _AUROC_HEADER_PT = 6.0

    @classmethod
    def _lane_y(cls, row, lane):
        """y of one lane within an organism row. y increases UPWARD, so lane 0 sits above.

        Deliberately the single source of this arithmetic: the ticks and the right-hand AUROC
        labels are positioned independently, and when this was written out twice, flipping the
        sign in one place silently swapped every printed AUROC onto the wrong lane. Nothing
        errors when that happens — the panel just lies.
        """
        return row + (1 if lane == 0 else -1) * cls._LANE_OFFSET

    def __init__(self, pct_df, excl_df, leak_df=None, ax=None, cells=EVENT_SQUARE,
                 auroc_axis=True):
        BasePlot.__init__(self, ax=ax, cells=cells)
        self.name = "hit_exclusivity_events"
        if pct_df is None or pct_df.empty or excl_df is None or excl_df.empty:
            self._unavailable()
            return

        excl = _prefer(excl_df, "subset")
        excl = excl[excl["subset"].isin(self._LANE_ORDER)]

        # Row order and whether the genus labels are drawn are ONE decision. Given the leakage
        # report, the rows take euos_overlap's order and the labels are dropped, so this panel can
        # sit beside it against that panel's single set of names. Without it there is nothing to
        # match, so the panel falls back to its own AUROC order and keeps its own labels — never
        # unlabelled rows in an unrelated order, which is the one combination that misreads.
        shared_order = overlap_row_order(leak_df)
        available = set(excl["pathogen"])
        if shared_order and available.issubset(shared_order):
            organisms = [o for o in shared_order if o in available]
            self.borrows_labels = True
        else:
            organisms = excl.groupby("pathogen")["auroc"].max() \
                .sort_values(ascending=True).index.tolist()
            self.borrows_labels = False

        auroc = {(r["pathogen"], r["subset"]): r["auroc"] for _, r in excl.iterrows()}
        positions, offsets, colors, self.counts = [], [], [], {}
        for i, org in enumerate(organisms):
            for lane, subset in enumerate(self._LANE_ORDER):
                vals = pct_df[(pct_df["pathogen"] == org)
                              & (pct_df["subset"] == subset)]["percentile"].values
                positions.append(np.asarray(vals, dtype=float))
                offsets.append(self._lane_y(i, lane))
                colors.append(EXCLUSIVITY_COLORS[subset])
                self.counts[(org, subset)] = len(vals)

        self.ax.eventplot(positions, lineoffsets=offsets,
                          linelengths=2 * self._LANE_HALF, colors=colors,
                          linewidths=self._TICK_WIDTH)
        self.ax.set_xlim(0, 1)
        # Index 0 at the BOTTOM — y increasing upward, not the inverted axis ``hbar`` uses.
        # Getting this backwards silently mirrors the panel against the one it pairs with.
        #
        # The limits are euos_overlap's (-0.6, n - 0.4), not a symmetric (-0.5, n - 0.5). The span
        # differs (7.2 vs 7.0), so with equal axes heights the row CENTRES would drift ~0.4 mm
        # apart by the top of the panel — enough to read as misaligned against shared labels.
        self.ax.set_ylim(-0.6, len(organisms) - 0.4)
        self.ax.set_yticks(range(len(organisms)))
        if self.borrows_labels:
            self.ax.set_yticklabels([])
        else:
            abbrev_ticks(self.ax, organisms, axis="y")
        # No chance line at 0.5. On the AUROC panel a 0.5 rule marks a meaningful threshold for a
        # single bar; here every lane is a spread of ticks that straddles it anyway, so the rule
        # adds a vertical stripe through all fourteen lanes without telling the reader anything the
        # ticks do not. The value is still on the x axis if it is wanted.
        self.label(xlabel="Rank percentile", ylabel="")

        self.auroc_labels = None
        if auroc_axis:
            self._auroc_axis(organisms, auroc)

    def _text_width_pt(self, text):
        """Rendered width of ``text`` in points, at the panel's small font size."""
        self.fig.canvas.draw()
        probe = self.ax.text(0, 0, text, fontsize=stylia.FONTSIZE_SMALL)
        width = probe.get_window_extent(
            renderer=self.fig.canvas.get_renderer()).width / self.fig.dpi * 72
        probe.remove()
        return width

    def _auroc_axis(self, organisms, auroc):
        """Print each organism's two AUROCs to the right of the plot as ``shared/exclusive``.

        Values rather than a scale: a second 0-1 axis would need its own ticks and reference line
        for fourteen numbers, and at 40 mm there is no room.

        Drawn as free text on a blended transform (x in axes fractions, y in data) rather than as
        twin-axis tick labels, because the pair has to share ONE line while each value keeps its
        lane's colour. A tick label is a single text object and so a single colour; three texts —
        value, slash, value — anchored right / centre / left give the same line with two colours
        and stay aligned whatever the glyph widths. ``clip_on=False`` because all of it sits
        outside the axes.
        """
        t = self.ax.get_yaxis_transform()          # x: axes fraction, y: data coordinates
        shared_c = EXCLUSIVITY_COLORS[self._LANE_ORDER[0]]
        excl_c = EXCLUSIVITY_COLORS[self._LANE_ORDER[1]]

        # Format every value first, THEN measure — the small-n asterisk widens the string, so a
        # pad measured on the bare number would let a marked value creep back over the plot.
        def fmt(org, subset):
            v = auroc.get((org, subset))
            if v is None:
                return None
            n = self.counts.get((org, subset), 0)
            return f"{v:.2f}" + (self._SMALL_N_MARK if n <= self._SMALL_N else "")

        pairs = [(o, fmt(o, self._LANE_ORDER[0]), fmt(o, self._LANE_ORDER[1])) for o in organisms]
        self.n_marked = sum(t is not None and t.endswith(self._SMALL_N_MARK)
                            for _, a, b in pairs for t in (a, b))

        # Where the slash goes, measured rather than assumed: the widest left-hand string has to
        # fit between the spine and the slash, so the offset is clearance + that width + gap.
        widest = max((self._text_width_pt(a) for _, a, _ in pairs if a), default=0.0)
        pad = self._AUROC_CLEAR_PT + widest + self._AUROC_GAP_PT

        def at(dx_pt, y, text, ha, color, dy_pt=0.0, va="center"):
            """One text pinned to the axes' right edge, offset by ``dx_pt`` points."""
            return self.ax.annotate(
                text, xy=(1.0, y), xycoords=t, textcoords="offset points",
                xytext=(dx_pt, dy_pt), ha=ha, va=va, color=color,
                fontsize=stylia.FONTSIZE_SMALL, annotation_clip=False)

        self.auroc_labels = []
        for i, (org, a, b) in enumerate(pairs):
            if a is None and b is None:
                continue
            at(pad, i, "/", "center", INK)
            if a is not None:
                at(pad - self._AUROC_GAP_PT, i, a, "right", shared_c)
            if b is not None:
                at(pad + self._AUROC_GAP_PT, i, b, "left", excl_c)
            self.auroc_labels.append((org, a, b))

        # Header over the column instead of a rotated axis label: rotated, "AUROC" is nearly as
        # tall as the 40 mm panel, and a vertical label beside a column of numbers reads as an axis
        # that is not there.
        at(pad, len(organisms) - 0.5, "AUROC", "center", INK,
           dy_pt=self._AUROC_HEADER_PT, va="bottom")


class HitExclusivityPlot(MetricByOrganismPlot):
    """Analysis 3 — exclusive vs shared (non-exclusive) hit AUROC per organism (dedup)."""

    def __init__(self, excl_df, ax=None, cells=(3, 3), legend=True, row_order=None):
        super().__init__(
            excl_df, group_col="subset", group_colors=EXCLUSIVITY_COLORS,
            group_order=["exclusive", "nonexclusive"], metric="auroc", xlabel="AUROC",
            title="", name="hit_exclusivity_auroc", count_col="n_active",
            ref=0.5, xlim=(0, 1), prefer_best=True, cells=cells, legend=legend, ax=ax,
            row_order=row_order)


class ActiveOverlapHeatmapPlot(BasePlot):
    """Pairwise Jaccard overlap between the 7 EU OpenScreen primary-assay active sets (label-only).
    Shows how much the organisms' hit-sets coincide — the backdrop for reading the cross-organism
    AUROCs. The self-overlap diagonal (=1) is blanked so the colour scale spans the informative
    off-diagonal overlaps."""

    def __init__(self, overlap_df, ax=None, cells=(4, 4)):
        BasePlot.__init__(self, ax=ax, cells=cells)
        self.name = "active_overlap_jaccard"
        if overlap_df is None or overlap_df.empty:
            self._unavailable()
            return
        codes = [c for c in SHARED_ORGANISMS if c in set(overlap_df["code_a"])]
        name_by_code = dict(zip(overlap_df["code_a"], overlap_df["pathogen_a"]))
        mat = overlap_df.pivot_table(index="code_a", columns="code_b", values="jaccard",
                                     aggfunc="first").reindex(index=codes, columns=codes)
        off = mat.where(~np.eye(len(codes), dtype=bool))  # blank the self-overlap diagonal
        vmax = float(np.nanmax(off.values)) if np.isfinite(off.values).any() else 1.0
        labels = [abbrev(name_by_code.get(c, c)) for c in codes]
        heatmap(self.ax, off, cmap=sequential_cmap("cobalt"),
                norm=Normalize(0.0, vmax), value_fmt="{:.2f}",
                text_light_when=lambda v: v > 0.6 * vmax,
                x_rotation=45, row_labels=labels, col_labels=labels)
        self.label(title="")


class ActiveOverlapPiePlot(BasePlot):
    """Directional overlap between the 7 EU OpenScreen primary-assay active sets, as a matrix of
    pies (label-only). The circle view of :class:`ActiveOverlapHeatmapPlot`.

    The **coloured wedge** is ``containment`` = |A ∩ B| / |A|: the share of the ROW organism's actives
    that are also active against the COLUMN organism.

    The matrix is therefore **not symmetric** and both triangles must be read: *P. aeruginosa* shares
    13 of its 14 actives with *A. baumannii* (93% of its row) while the reverse cell is 13/57 = 23%.
    That asymmetry is exactly what a Jaccard heatmap hides, since a union denominator lets the larger
    set dominate — the pair scores 0.22 either way. The diagonal needs no special case: containment
    with itself is 1, so it renders as a full circle.

    **Every pie is the same size, and each organism's active count is in its y tick label instead.**
    An earlier version scaled circle area by the row's active count. That works at 120 mm, but this
    panel is a 45 mm paper panel where the tick labels take ~60 % of the width and the cell pitch is
    2.4 mm: a linear area scale then put *S. aureus* (378 actives) at 2.1 mm and *P. aeruginosa* (14)
    at **0.40 mm**, which is a dot with no readable wedge. Since the wedge is the measurement and the
    size was context, the size encoding was dropped and the counts became text — exact, rather than
    estimated off an area a reader cannot resolve. ``self.set_sizes`` still publishes the counts.

    Each ROW carries its pathogen's hue, matching ``pathogen_activity_ratios`` in step 03 via
    :func:`plotting_colors.pathogen_activity_colors` — colour is redundant with the row label, exactly
    as it is there, and exists so a reader can follow one pathogen across the two steps' panels. Note
    those hues do NOT agree with ``SHARED_ORGANISM_COLORS`` in general (only the
    ``PATHOGEN_HUE_SWAPS`` pair is reconciled), so the step-05 panels that break down BY pathogen
    (``ExclusiveHitModelRankPlot``) take the same step-03 hues rather than letting two palettes
    disagree inside one figure.

    Pie radii are sized in POINTS from the measured axes box, not in data units, so the marks stay
    round whatever ``tight_layout`` does to the axes. This requires layout to be final before the
    marks are drawn — hence the ``tight_layout`` + ``draw`` before ``pie_scatter``.
    """

    #: Pie diameter as a fraction of the matrix cell pitch. Not larger: the centre gridlines run
    #: UNDER the pies, so they are only visible in the gap between neighbours, and at 0.86 that gap
    #: closes to ~0.3 mm and the grid disappears.
    MAX_DIAM_FRAC = 0.76
    #: Stroke around each slice. The "rest" slice is white, so a white edge (``pie_scatter``'s
    #: default) would erase the circle's own outline against the page.
    EDGE_WIDTH = 0.4

    def __init__(self, overlap_df, row_colors=None, ax=None, cells=SMALL_SQUARE):
        BasePlot.__init__(self, ax=ax, cells=cells)
        self.name = "active_overlap_containment"
        self.pie_area = None
        self.set_sizes = {}
        need = {"code_a", "code_b", "n_a", "n_intersect"}
        if overlap_df is None or overlap_df.empty or not need.issubset(overlap_df.columns):
            self._unavailable()
            return
        codes = [c for c in SHARED_ORGANISMS if c in set(overlap_df["code_a"])]
        if not codes:
            self._unavailable()
            return
        name_by_code = dict(zip(overlap_df["code_a"], overlap_df["pathogen_a"]))
        d = overlap_df.copy()
        # Derived rather than required: the column was added to the summary CSV alongside `jaccard`,
        # but it is a pure function of columns that were always shipped, so older CSVs still plot.
        if "containment" not in d.columns:
            d["containment"] = d["n_intersect"] / d["n_a"]
        frac = d.pivot_table(index="code_a", columns="code_b", values="containment",
                             aggfunc="first").reindex(index=codes, columns=codes)
        n_by_code = d.groupby("code_a")["n_a"].first()
        n_row = np.array([float(n_by_code[c]) for c in codes])
        self.set_sizes = {c: int(v) for c, v in zip(codes, n_row)}

        k = len(codes)
        labels = [abbrev(name_by_code.get(c, c)) for c in codes]
        self.ax.set_xticks(range(k))
        # 90 deg, not the usual 45: rotated labels on a shared baseline clear each other only when
        # `pitch * sin(angle)` beats their line height, and at this footprint the pitch is ~2.4 mm —
        # 45 deg yields 1.7 mm against a ~1.9 mm line, so the genus names overlapped.
        self.ax.set_xticklabels(labels, rotation=90, ha="center")
        self.ax.set_yticks(range(k))
        # The row's active count rides on the y label — it is the denominator of every wedge in that
        # row, and it no longer has a size encoding to live in.
        self.ax.set_yticklabels([f"{lab} ({int(n)})" for lab, n in zip(labels, n_row)])
        self.ax.set_xlim(-0.5, k - 0.5)
        self.ax.set_ylim(k - 0.5, -0.5)   # first organism on top, as in the heatmap panels
        self.ax.set_aspect("equal")        # square cells, so one pitch governs both axes
        # Faint guides: with a sparse circle matrix there is otherwise nothing to carry the eye from a
        # pie back to its row and column labels. Drawn on the MAJOR ticks, so a line runs through the
        # circle centres and points straight at the label — not on cell boundaries, which would box
        # every mark in and add a second grid of edges competing with the circles' own outlines.
        # `set_axisbelow` keeps them under the pies; a line over a white slice would read as a radius.
        self.ax.grid(which="major", color=REFERENCE_LINE, linewidth=0.3)
        self.ax.set_axisbelow(True)
        self.label(xlabel="Also active against", ylabel="Actives of")

        # Layout must be final before radii are fixed in points (see the class docstring).
        import matplotlib.pyplot as plt
        plt.figure(self.fig.number)
        plt.tight_layout()
        self.fig.canvas.draw()
        bb = self.ax.get_window_extent()
        pitch = min(bb.width, bb.height) * 72.0 / self.fig.dpi / k
        self.pie_area = (self.MAX_DIAM_FRAC * pitch) ** 2   # one marker area (pt^2) for every pie

        # One pie_scatter call per row, because the wedge colour is per row: `pie_scatter` takes one
        # (part, rest) pair per call.
        for i, code in enumerate(codes):
            xs, fs = [], []
            for j in range(k):
                v = frac.iat[i, j]
                if np.isfinite(v):
                    xs.append(j)
                    fs.append(float(v))
            if not xs:
                continue
            part = (row_colors or {}).get(code) or hue("cobalt")
            pie_scatter(self.ax, xs, [i] * len(xs), fs, (part, "white"),
                        s=self.pie_area, edgecolor=INK, linewidth=self.EDGE_WIDTH)


class ActiveOverlapKeyPlot(BasePlot):
    """Standalone key for :class:`ActiveOverlapPiePlot` — the wedge scale, drawn once as a strip.

    Only one scale to explain, since the pies are all one size: what a given wedge angle means. The
    swatches go through :func:`pie_scatter`, the same code path as the marks, so the key cannot drift
    from what it explains. (It carried a ``nested_size_legend`` block too while circle area encoded
    the row's active count; that encoding is gone — see :class:`ActiveOverlapPiePlot`.)

    Swatches are drawn NEUTRAL, not in any pathogen's hue: in the panel colour identifies the row's
    organism, so a key swatch in one organism's colour would read as a statement about that organism
    rather than about the angle scale.

    Swatches are deliberately LARGER than the panel's own 2.1 mm pies. A key is read once, up close,
    to calibrate the eye; reproducing the mark at its true size would make the 25 %/50 % distinction
    as hard to see here as it is in the matrix, which defeats the point of having a key.
    """

    #: Containment fractions shown in the strip.
    WEDGE_FRACTIONS = (0.25, 0.5, 0.75, 1.0)
    #: Marker area (pt^2) of the swatches — fixed, and independent of the panel's pie size.
    WEDGE_AREA = 90.0

    def __init__(self, ax=None, cells=(0.5, 1.5)):
        BasePlot.__init__(self, ax=ax, cells=cells)
        self.name = "active_overlap_containment_key"
        self.ax.set_xlim(0, 1)
        self.ax.set_ylim(0, 1)
        self.ax.set_axis_off()

        fracs = list(self.WEDGE_FRACTIONS)
        xs = np.linspace(0.10, 0.90, len(fracs))
        pie_scatter(self.ax, xs, [0.62] * len(fracs), fracs,
                    (REFERENCE_LINE, "white"), s=self.WEDGE_AREA, edgecolor=INK,
                    linewidth=ActiveOverlapPiePlot.EDGE_WIDTH)
        for x, f in zip(xs, fracs):
            self.ax.text(x, 0.26, f"{int(round(f * 100))}%", ha="center", va="top", fontsize=5)
        self.ax.text(0.5, 1.0, "share of the row's actives", ha="center", va="top", fontsize=5)


class HitPromiscuityPlot(BasePlot):
    """Hit promiscuity (label-only) — how many EU OpenScreen actives are hits in 1, 2, ... 7 of
    the primary assays. The per-compound counterpart of the Jaccard overlap heatmap: most actives
    are organism-specific singletons, while a small tail is pan-active across the panel.

    The y axis is log-scaled because the distribution spans two orders of magnitude (hundreds of
    singletons vs a handful of 7-pathogen hits); every bar is annotated with its exact count.

    Bars fade with **promiscuity itself** (:func:`plotting_colors.count_shades` over ``n_pathogens``)
    — palest for the organism-specific singletons, darkest for the pan-active 7-pathogen bin, so the
    gradient runs the same direction as the panel's subject. Driving it from the bar's own *count*
    was tried and reads backwards: it makes the singleton bin the darkest thing on the panel, and it
    doubles the y axis rather than adding anything. Colour is never the sole encoding here — the exact
    count is annotated on every bar, so the fade is a redundant, scannable cue.

    The fit is linear, not log: 1-7 is a short ordinal range, where a log fit would bunch the high
    (darkest, and most interesting) end together.

    Bars are drawn from a **finite positive baseline** (``Y_FLOOR``), never matplotlib's default
    ``bottom=0``: on a log axis 0 maps to negative infinity, so the default emitted paths reaching
    y = -28,828 pt on a 130 pt page. The PNG clips them away, which is what makes this bug invisible
    in the raster preview and visible only in the vector PDF. See the log-axis gotcha in
    ``docs/figure_conventions.md`` and ``PipelineFunnelPlot``.
    """

    #: Finite positive baseline for the bars, and the bottom of the log y axis.
    Y_FLOOR = 0.7

    def __init__(self, prom_df, ax=None, cells=(3, 3)):
        BasePlot.__init__(self, ax=ax, cells=cells)
        self.name = "hit_promiscuity"
        if prom_df is None or prom_df.empty or "n_molecules" not in prom_df.columns:
            self._unavailable()
            return
        d = prom_df.sort_values("n_pathogens")
        x = d["n_pathogens"].to_numpy()
        y = d["n_molecules"].to_numpy(dtype=float)
        self.ax.bar(x, np.maximum(y - self.Y_FLOOR, 0.0), bottom=self.Y_FLOOR,
                    color=count_shades(x))
        self.ax.set_yscale("log")
        self.ax.set_ylim(self.Y_FLOOR, float(y.max()) * 3.0 if y.max() > 0 else 10.0)
        self.ax.set_xticks(x)
        for xi, yi in zip(x, y):
            if yi <= 0:
                continue
            self.ax.annotate(f"{int(yi)}", (xi, yi), ha="center", va="bottom",
                             fontsize=5, xytext=(0, 1.5), textcoords="offset points")
        total = int(d["n_molecules"].sum())
        # Terse axis labels: at 45 mm a sentence-length label doubles the panel's crop. The full
        # phrasing belongs in the caption.
        self.label(xlabel="Pathogens hit", ylabel=f"Actives (n = {total:,})", title="")


class ExclusiveHitModelRankPlot(BasePlot):
    """For each exclusive hit (active in exactly 1 of the 7 primary assays), where its own
    pathogen's model ranks it among all 7 models: rank 1 = its own pathogen scores it highest,
    rank 2 = one other pathogen's model ranks it higher, and so on to rank 7.

    One panel per ranking mode — ``raw`` consensus scores as-is, or ``percentile`` (each score first
    converted to its percentile within that model's own library distribution, which puts the models
    on a common scale; the raw scores are not calibrated across models).

    Bar height is the molecule count, split into segments coloured by the pathogen each hit belongs
    to (from the ``n_<code>`` columns of the summary CSV). Segment sizes therefore track how many
    exclusive hits each pathogen has in the first place — the legend gives those totals so a large
    segment isn't misread as a per-pathogen effect.

    **No chance line.** A dashed neutral rule used to mark ``n_chance`` (the count expected if the
    own-model rank were random, n_total / 7). It was removed: nothing on the panel said what it was,
    a bare grey rule invites being read as a threshold or an axis break, and it ran straight through
    the bar-top numbers. The value is still carried in the summary CSV's ``n_chance`` column, so it
    belongs in the caption, where it can be named.

    ``pathogen_colors`` overrides the default ``SHARED_ORGANISM_COLORS``. It is passed the step-03
    ``pathogen_activity_ratios`` hues so that this panel and
    :class:`ActiveOverlapPiePlot` — the two step-05 panels that break a total down BY pathogen, and
    which share a page — use ONE colour per pathogen. The two palettes disagree outright (step 03
    gives *S. aureus* the lime the shared palette gives *E. coli*, and vice versa), so letting each
    panel keep its own would put the same pathogen in two colours and two pathogens in each other's."""

    #: Gap between a bar's top and its count, in points.
    LABEL_PAD = 3.0

    def __init__(self, rank_df, ranking="raw", name_by_code=None, dedup=False,
                 pathogen_colors=None, ax=None, cells=(3, 3), legend=True):
        BasePlot.__init__(self, ax=ax, cells=cells)
        self.name = f"exclusive_hit_model_rank_{ranking}" + ("_dedup" if dedup else "")
        self.colors = dict(SHARED_ORGANISM_COLORS)
        self.colors.update(pathogen_colors or {})
        if rank_df is None or rank_df.empty or "rank" not in rank_df.columns:
            self._unavailable()
            return
        d = rank_df[rank_df["ranking"] == ranking].sort_values("rank") \
            if "ranking" in rank_df.columns else rank_df.sort_values("rank")
        if d.empty:
            self._unavailable()
            return
        x = d["rank"].to_numpy()
        y = d["n_molecules"].to_numpy(dtype=float)
        # stack by the hit's own pathogen when the breakdown columns are present
        codes = [c for c in SHARED_ORGANISMS if f"n_{c}" in d.columns and d[f"n_{c}"].sum() > 0]
        if codes:
            bottom = np.zeros(len(d))
            for code in codes:
                seg = d[f"n_{code}"].to_numpy(dtype=float)
                self.ax.bar(x, seg, bottom=bottom, color=self.colors[code])
                bottom += seg
            self.segment_totals = {c: int(d[f"n_{c}"].sum()) for c in codes}
            if legend:
                self.legend({f"{abbrev(c, name_by_code)} ({self.segment_totals[c]})":
                             self.colors[c] for c in codes},
                            loc="upper right", ncol=2)
        else:
            self.ax.bar(x, y, color=hue("cobalt", lighten=0.55))
        # No grid. Every bar carries its exact count, so the gridlines were already redundant — and
        # keeping them means a bar whose top lands near a gridline puts its number straight onto it
        # (37 sat on the 40 line). Padding cannot fix that in general: the collision depends on where
        # the counts happen to fall relative to the tick locator, so any fixed gap is one data refresh
        # away from breaking. Removing the grid makes the labels legible by construction.
        self.ax.grid(False)
        self.ax.set_xticks(x)
        for xi, yi in zip(x, y):
            if yi > 0:
                self.ax.annotate(f"{int(yi)}", (xi, yi), ha="center", va="bottom",
                                 fontsize=5, xytext=(0, self.LABEL_PAD),
                                 textcoords="offset points")
        total = int(d["n_total"].iloc[0]) if "n_total" in d.columns else int(y.sum())
        suffix = ""   # the dedup status is in the folder name and the caption, not the axis
        n_ranked = int(d["n_models_ranked"].iloc[0]) if "n_models_ranked" in d.columns else len(x)
        self.label(xlabel=f"Own-model rank (of {n_ranked})",
                   ylabel=f"Exclusive hits (n = {total}{suffix})", title="")


class ScoreByHitClassPlot(BasePlot):
    """Boxes of a per-compound aggregated consensus score, one box per EU OpenScreen hit class.

    A reusable panel type (both score-distribution figures below build on it). Boxes are drawn from
    the PRECOMPUTED statistics — the inactive class has ~10^5 compounds and is never shipped
    per-molecule — while the active classes, whose individual values the summary CSV does carry,
    also get a jittered point overlay. Raw scores: training-set compounds are deliberately
    included, so these panels describe the score distribution, not out-of-sample performance.
    """

    #: Box body and swarm geometry, matching the step-01 technical boxes (``TaskMetricBoxPlot``):
    #: a swarm band as wide as the box, so a few hundred points spread instead of piling into a
    #: line. The defaults in ``box_from_stats`` (0.34 / 0.12) are sized for panels with far fewer
    #: points than the 428 actives here.
    _BOX_WIDTH = 0.5
    _JITTER = 0.20
    #: Swarm point area (pt^2), as ``box_from_stats``' ``point_size``. The base keeps the helper's
    #: default; subclasses that need a lighter swarm override it (see
    #: :class:`ConsensusSumByHitClassPlot`).
    _POINT_SIZE = 6

    def __init__(self, stats_df, actives_df, classes, score_col, ylabel, name,
                 ax=None, cells=(3, 3), show_counts=True, label_rotation=0,
                 ylabel_right=False):
        BasePlot.__init__(self, ax=ax, cells=cells)
        self.name = name
        if stats_df is None or stats_df.empty or "hit_class" not in stats_df.columns:
            self._unavailable()
            return
        order = [c for c in classes if c in set(stats_df["hit_class"])]
        if not order:
            self._unavailable()
            return
        by_class = stats_df.set_index("hit_class")
        rng = np.random.default_rng(RANDOM_SEED)
        for i, hit_class in enumerate(order):
            points = None
            if (actives_df is not None and not actives_df.empty
                    and hit_class != "inactive" and score_col in actives_df.columns):
                points = actives_df.loc[actives_df["hit_class"] == hit_class,
                                        score_col].to_numpy()
            # No fills at all — the house default. This panel used to fill every class, which hid
            # all 428 active points behind their own box (only those spilling past the quartiles
            # were visible) and forced the median to INK, since a same-hue median cannot be seen
            # on an opaque body. Unfilled, the swarm carries the distribution and the median takes
            # the class colour. The inactive box has no points behind it and is simply an outline,
            # which is the right weight for a background class.
            box_from_stats(self.ax, by_class.loc[hit_class], i, HIT_CLASS_COLORS[hit_class],
                           face=None, points=points, rng=rng,
                           width=self._BOX_WIDTH, jitter_width=self._JITTER,
                           point_size=self._POINT_SIZE)
        n_models = int(stats_df["n_models_aggregated"].iloc[0]) \
            if "n_models_aggregated" in stats_df.columns else len(SHARED_ORGANISMS)
        # ``show_counts=False`` drops the "(n = ...)" line from the tick labels. It is a SIZE
        # decision, not a preference: the four class labels carrying their counts are 32.3 mm of
        # type against a ~33 mm axes at 40 mm, so they overlap (-1.7 mm between the first two).
        # Without the counts the measured gap is +0.2 mm, but that measurement is a LIE: tick
        # labels' window extents do not reflect their final placement, and the panel rendered with
        # "inactive" and "exclusive" run together. Verified by looking at the output, not by
        # measuring it. ``label_rotation=45`` is what actually clears them at this width.
        #
        # The counts do not disappear: they stay in the box-stats CSV, are exposed here as
        # ``self.class_counts`` for the entry point to print, and **a caption must carry them** —
        # the broad class of the dedup twin rests on 14 compounds, which no reader should have to
        # go looking for. Same trade-off, and the same resolution, as step 01's TaskMetricBoxPlot.
        self.class_counts = {c: int(by_class.loc[c, "n"]) for c in order}
        self.ax.set_xticks(range(len(order)))
        self.ax.set_xticklabels(
            [f"{_pretty(c)}\n(n = {self.class_counts[c]:,})" if show_counts else _pretty(c)
             for c in order],
            rotation=label_rotation, ha="right" if label_rotation else "center",
            rotation_mode="anchor" if label_rotation else None)
        self.label(xlabel="", ylabel=ylabel.format(n_models=n_models), title="")
        # ``ylabel_right`` moves the scale (ticks, tick labels AND axis label) to the right spine.
        # A placement decision for the Illustrator layout, not a style one: this panel sits to the
        # RIGHT of `hit_exclusivity_events`, whose own AUROC column runs down its right edge, so a
        # left-hand scale here would put two label stacks back to back in the gutter between them.
        # Set after ``label`` because the label position has to be moved after stylia writes it.
        if ylabel_right:
            self.ax.yaxis.tick_right()
            self.ax.yaxis.set_label_position("right")


class ConsensusMaxRocPlot(BasePlot):
    """ROC of active vs inactive under the best-of-7 score — the curve behind the box panel.

    Exactly the same numbers as :class:`ConsensusMaxByActivityPlot`, read the other way: the box
    shows where each class sits, this shows the retrieval trade-off that separation buys. A pair
    of boxes cannot answer "how much of the library must I screen to recover half the hits"; a ROC
    cannot show that a quarter of the actives score below the inactive median. Neither replaces
    the other.

    Drawn with the shared :func:`plotting_utils.roc_panel`, so it reads identically to the
    per-organism ROC grid — chance diagonal, shaded curve, AUC and class counts in the corner —
    and in the active class's colour, tying it to the box panel it accompanies.

    The curve comes from the precomputed summary: the ~10^5 inactive scores it is built from are
    never shipped per-molecule, so it cannot be recomputed at figure time (see
    ``eval_euopenscreen.run_consensus_max_by_activity``).
    """

    def __init__(self, roc_df, normalized=False, dedup=False, ax=None, cells=(3, 3)):
        BasePlot.__init__(self, ax=ax, cells=cells)
        name = "consensus_max_percentile_roc" if normalized else "consensus_max_roc"
        self.name = name + ("_dedup" if dedup else "")
        if roc_df is None or roc_df.empty or "fpr" not in roc_df.columns:
            self._unavailable()
            return
        d = roc_df.sort_values("fpr")
        self.auroc = float(d["auroc"].iloc[0])
        roc_panel(self.ax, d["fpr"].to_numpy(), d["tpr"].to_numpy(), self.auroc,
                  int(d["n_pos"].iloc[0]), int(d["n_neg"].iloc[0]),
                  HIT_CLASS_COLORS["active"], xlabel="FPR", ylabel="TPR")


class ConsensusSumByHitClassPlot(ScoreByHitClassPlot):
    """Summed consensus score (over the 7 shared-organism models) by hit class: inactive everywhere
    / hit in exactly 1 pathogen (exclusive) / narrow-spectrum (2-3) / broad-spectrum (>3).

    ``dedup=True`` reads the leakage-filtered twin, where every compound in ANY of the 7 models'
    ChEMBL training sets is dropped from all four classes. **The two are not the same picture at a
    different scale.** Leakage rises steeply with promiscuity — 17% of exclusive hits are in a
    training set (68 of 390), against 40% of narrow (61 of 153) and 69% of broad ones (31 of 45) —
    so the dedup panel loses most of its broad class and comparatively little of its exclusive one.
    Any read of the class trend has to say which variant it is describing. Counts are the ones the
    box-stats CSVs carry; re-derive them from there rather than trusting this comment.
    """

    #: Lighter swarm than the base's 6 pt^2. Four classes share the width here and the exclusive one
    #: carries ~400 points, so at 6 pt^2 (a 0.97 mm dot) its swarm read as a solid mass and the box
    #: outline it sits behind was lost in it. At 3 pt^2 (0.69 mm) the individual compounds separate
    #: and the quartiles stay legible. Applied to both twins so the two panels stay comparable.
    _POINT_SIZE = 3

    def __init__(self, stats_df, actives_df=None, dedup=False, ax=None, cells=(3, 3),
                 show_counts=True, label_rotation=0, ylabel_right=False):
        super().__init__(stats_df, actives_df, classes=HIT_CLASSES,
                         score_col="consensus_sum",
                         ylabel="summed consensus score ({n_models} models)",
                         name="consensus_sum_by_hit_class" + ("_dedup" if dedup else ""),
                         ax=ax, cells=cells, show_counts=show_counts,
                         label_rotation=label_rotation, ylabel_right=ylabel_right)


class ConsensusMaxByActivityPlot(ScoreByHitClassPlot):
    """Maximum score across the 7 shared-organism models — how confident the single most confident
    model is — for compounds inactive in every assay they were tested in vs compounds that are a hit
    in one or more pathogens (regardless of how many).

    ``normalized=True`` reads the variant whose per-model scores were converted to within-model
    library percentiles before taking the max (so the max is not biased towards whichever model
    outputs the highest values), with that maximum then re-ranked over the library so the axis is a
    plain library percentile rather than a best-of-7 value crowded against 1.0. Re-ranking is
    monotone, so this panel and the raw one differ only in axis, not in ordering."""

    def __init__(self, stats_df, actives_df=None, normalized=False, dedup=False,
                 ax=None, cells=(3, 3)):
        # Short forms: spelled out these ran to ~70 characters, and a rotated y label that long makes
        # the panel's crop taller than its footprint. Model count and dedup status go in the caption.
        ylabel = "best-model percentile" if normalized else "max consensus score"
        name = "consensus_max_percentile_by_activity" if normalized \
            else "consensus_max_by_activity"
        if dedup:
            name += "_dedup"
        super().__init__(stats_df, actives_df, classes=ACTIVITY_CLASSES,
                         score_col="consensus_max", ylabel=ylabel, name=name,
                         ax=ax, cells=cells)
        # The separation as one number. It cannot be recomputed here — the inactive scores are
        # never shipped — so it rides along in the stats CSV; see run_consensus_max_by_activity.
        self.auroc = None
        if self.is_available and stats_df is not None and "auroc" in stats_df.columns:
            value = stats_df["auroc"].dropna()
            if not value.empty:
                self.auroc = float(value.iloc[0])
                self.ax.text(0.03, 0.03, f"AUROC {self.auroc:.2f}", transform=self.ax.transAxes,
                             ha="left", va="bottom", fontsize=stylia.FONTSIZE_SMALL,
                             # Same intent as LEGEND_KW ("legends over data stay readable":
                             # white, 70% opaque, no edge) but spelled for a patch — a bbox
                             # rejects the legend-only frameon/framealpha keys.
                             bbox=dict(boxstyle="square,pad=0.25", facecolor="white",
                                       alpha=0.7, edgecolor="none"))


class CrossOrganismHeatmapPlot(BasePlot):
    """Analysis 4 — model x EU OpenScreen assay AUROC matrix (off-diagonal = cross-organism)."""

    def __init__(self, cross_df, ax=None, cells=(4, 4)):
        BasePlot.__init__(self, ax=ax, cells=cells)
        self.name = "cross_organism_heatmap"
        if cross_df is None or cross_df.empty:
            self._unavailable()
            return
        use_set = "dedup" if "dedup" in set(cross_df["set"]) else "raw"
        d = cross_df[cross_df["set"] == use_set]
        if d.empty:
            self._unavailable()
            return
        # Rows: shared organisms first (diagonal block), then non-shared, in a stable order.
        name_by_code = dict(zip(d["model_code"], d["model_pathogen"]))
        assay_name = dict(zip(d["assay_code"], d["assay_pathogen"]))
        col_codes = [c for c in SHARED_ORGANISMS if c in set(d["assay_code"])]
        present = [c for c in d["model_code"].unique()]
        row_codes = ([c for c in SHARED_ORGANISMS if c in present]
                     + [c for c in present if c not in SHARED_ORGANISMS])
        mat = d.pivot_table(index="model_code", columns="assay_code",
                            values="auroc", aggfunc="first")
        mat = mat.reindex(index=row_codes, columns=col_codes)
        mat.index = [name_by_code.get(c, c) for c in mat.index]
        mat.columns = [abbrev(assay_name.get(c, c)) for c in col_codes]
        # highlight the diagonal (model organism == assay organism)
        highlight = [(ri, ci) for ri, rc in enumerate(row_codes)
                     for ci, cc in enumerate(col_codes) if rc == cc]
        heatmap(self.ax, mat, cmap=diverging_cmap(),
                norm=TwoSlopeNorm(0.5, vmin=0.0, vmax=1.0),
                text_light_when=lambda v: v > 0.75 or v < 0.25,
                highlight=highlight, x_rotation=45,
                row_labels=[abbrev(p) for p in mat.index])
        self.label(title="")


class SpecificityIndexPlot(BasePlot):
    """Analysis 4 (companion) — per-model specificity index (own − mean cross AUROC)."""

    def __init__(self, spec_df, ax=None, cells=(3, 3), row_order=None):
        BasePlot.__init__(self, ax=ax, cells=cells)
        self.name = "specificity_index"
        if spec_df is None or spec_df.empty or spec_df["specificity_index"].dropna().empty:
            self._unavailable()
            return
        # REVERSED: `specificity_bars` draws through `plotting_utils.hbar`, which inverts the y axis so
        # its first element lands on TOP. Every other panel here is bottom-row-first. Passing the
        # shared order unreversed would put this panel's rows in the exact opposite sequence to the
        # rest — the one failure mode this whole change exists to prevent.
        order = _apply_row_order(spec_df["pathogen"].dropna().unique(), row_order)
        specificity_bars(self.ax, spec_df, title="",
                         order=list(reversed(order)) if order else None)


class SubmodelAurocPlot(BasePlot):
    """Per-pathogen: AUROC of every sub-model output on the pathogen's own EU OpenScreen assay
    (dedup), one dot per sub-model, ``consensus_score`` marked with a star. A wide horizontal
    spread means the ensemble members disagree in quality."""

    def __init__(self, grp, ax=None):
        d = grp.copy()
        d["_r"] = (d["set"] == "dedup").astype(int)  # prefer dedup per feature
        d = d.sort_values("_r").drop_duplicates(["feature"], keep="last")
        d = d.sort_values("auroc", ascending=True).reset_index(drop=True)
        n = len(d)
        BasePlot.__init__(self, ax=ax, cells=(max(2, math.ceil(n / 3)), 3))
        self.name = f"{grp['code'].iloc[0]}_submodel_auroc"
        if n == 0:
            self._unavailable()
            return
        y = np.arange(n)
        is_cons = (d["feature"] == "consensus_score").values
        self.ax.scatter(d.loc[~is_cons, "auroc"], y[~is_cons], color=hue("cobalt"),
                        zorder=3, label="sub-model")
        if is_cons.any():
            self.ax.scatter(d.loc[is_cons, "auroc"], y[is_cons], color=hue("crimson"),
                            marker="*", s=140, zorder=4, label="consensus")
        self.ref_line(0.5, axis="x")
        self.ax.set_yticks(y)
        self.ax.set_yticklabels(d["feature"].tolist(), fontsize=5)
        self.ax.set_xlim(0, 1)
        self.ax.legend(fontsize=5, loc="lower right", **LEGEND_KW)
        self.label(xlabel="AUROC", title=abbrev(grp["pathogen"].iloc[0]))


class SubmodelCorrPlot(BasePlot):
    """Per-pathogen: pairwise Spearman correlation between sub-model scores over the library.
    Cobalt = strongly correlated (rank alike), crimson = anti-correlated, white ≈ 0. Low
    off-diagonal blocks flag sub-models that rank compounds differently."""

    def __init__(self, corr, code, pathogen, ax=None):
        n = 0 if corr is None else corr.shape[0]
        side = min(6, max(3, round(n * 0.35))) if n else 3
        BasePlot.__init__(self, ax=ax, cells=(side, side))
        self.name = f"{code}_submodel_corr"
        if corr is None or corr.empty:
            self._unavailable()
            return
        labels = list(corr.columns)
        heatmap(self.ax, corr, cmap=diverging_cmap(), norm=Normalize(-1, 1),
                annotate=(n <= 12), text_light_when=lambda v: abs(v) > 0.6,
                value_fmt="{:.2f}", annot_fontsize=3.5, x_rotation=90, colorbar=True,
                row_labels=labels, col_labels=labels)
        self.ax.tick_params(labelsize=4)
        self.label(title=abbrev(pathogen))


class SubmodelAurocSummaryPlot(BasePlot):
    """Cross-pathogen summary of sub-models vs the consensus. One row per shared organism: a small
    dot per sub-model output's own-assay AUROC (dedup), a thin guide line across that organism's
    min–max, and the ``consensus_score`` drawn as a larger bubble — shows at a glance where the
    consensus sits within the spread of its ensemble members across all organisms.

    **Deliberately styled as the twin of step 03's**
    :class:`plots_chembl_performance.PathogenConsensusAurocPlot`: same marks (periwinkle member
    dots, min–max guide line in the neutral, crimson summary dot), same in-axes key. The two are
    the paper's two AUROC dot plots — one on held-out ChEMBL CV, one on EU OpenScreen — and a
    reader moving between them should not have to relearn the encoding. What the summary dot means
    is the one difference: there it is the *mean* of the CV AUROCs, here it is the actual
    ``consensus_score`` output, which is a prediction in its own right and need not sit at the
    members' centre. Sizes stay at this panel's own values (step 03 scales its mean dot by dataset
    size; there is no such quantity here), and the x axis keeps its 0.2 floor because sub-models do
    reach below chance — step 03's 0.68 floor would clip them.

    Kept at ``WIDE_SUBMODEL`` rather than step 03's 62 mm square: this panel is one of the paper's
    45 mm row-B panels and growing it would break that row. That footprint is also why the key is
    NOT in-axes as step 03's is, the one place the twinning stops: with 7 genus tick labels taking
    half of the 45 mm crop the axes is 21.7 mm wide, and the two-entry key measures 14.4 mm against
    a 7.9 mm clear corner (the bottom two rows floor at AUROC 0.490). Measured, not estimated — no
    anchor placement recovers an 84 % overrun, so the key stays a standalone ``KeyPanel`` built from
    the same :data:`SUBMODEL_KEY_ENTRIES`. ``legend=True`` still draws it in-axes for anyone
    rendering this panel larger; the placement constants below are kept for that path.
    """

    #: Weight of the min-max guide line per row, matching step 03's ``_RANGE_LINEWIDTH``. At the
    #: house 0.5 pt it stays a guide — the eye reads the dots and the consensus, with the line only
    #: tying a row together.
    _RANGE_LINEWIDTH = 0.5
    #: Key row spacing, in font-size units, as step 03's. Well above matplotlib's 0.5 default
    #: because the rows have to clear the consensus swatch, not just the text. Unlike step 03 the
    #: swatches keep their marks' TRUE sizes (see :data:`SUBMODEL_KEY_ENTRIES`) — step 03 equalises
    #: them only because its summary dot is size-encoded, so no one size would be honest there.
    _LEGEND_LABELSPACING = 0.9
    #: Key corner, and its anchor in AXES fractions. **LOWER left, where step 03 uses upper left** —
    #: the one place the two panels deliberately differ, because the clear corner is a property of
    #: the data, not of the style. Step 03's rows are dense and its members never reach its 0.68
    #: floor, so its top-left is empty. Here the sub-models spread far below chance and the top rows
    #: are the ones that reach leftmost (the best organism's members go down to 0.441), while the
    #: bottom two rows floor at 0.490 — so the empty corner is the bottom one, with ~36 % of the
    #: width against the top corner's ~30 %. Anchored on the spine, not inset, because 36 % of a
    #: ~30 mm axes is ~11 mm and "consensus" needs nearly all of it.
    _LEGEND_LOC = "lower left"
    _LEGEND_ANCHOR = (0.0, 0.0)

    def __init__(self, auroc_df, ax=None, cells=(3, 3), legend=True, row_order=None):
        BasePlot.__init__(self, ax=ax, cells=cells)
        self.name = "submodel_auroc_summary"
        if auroc_df is None or auroc_df.empty:
            self._unavailable()
            return
        d = auroc_df.copy()
        d["_r"] = (d["set"] == "dedup").astype(int)  # prefer dedup per (pathogen, feature)
        d = d.sort_values("_r").drop_duplicates(["pathogen", "feature"], keep="last")
        # Rows follow the shared `euos_overlap` order when supplied. Otherwise order by consensus
        # AUROC (best on top), falling back to the mean when there is no consensus row.
        cons = d[d["feature"] == "consensus_score"].set_index("pathogen")["auroc"]
        order = _apply_row_order(d["pathogen"].unique(), row_order)
        if order is None:
            if cons.empty:
                order = d.groupby("pathogen")["auroc"].mean() \
                    .sort_values(ascending=True).index.tolist()
            else:
                order = cons.sort_values(ascending=True).index.tolist()
        idx = {p: i for i, p in enumerate(order)}
        sub = d[d["feature"] != "consensus_score"]
        con = d[d["feature"] == "consensus_score"]
        # Min-max guide line per row, spanning the SUB-MODELS only — the spread the panel exists to
        # show. Step 03's line spans the values its small dots mark and its mean dot cannot fall
        # outside them; here the consensus is a separate model output that CAN sit beyond its
        # members' range, and when it does that is a finding, not an error to be hidden by
        # stretching the line to cover it.
        for p, i in idx.items():
            vals = sub.loc[sub["pathogen"] == p, "auroc"].to_numpy(dtype=float)
            if len(vals):
                self.ax.plot([vals.min(), vals.max()], [i, i], color=REFERENCE_LINE,
                             linewidth=self._RANGE_LINEWIDTH, zorder=1)
        # Marker areas are absolute, so they must come down with the footprint. At WIDE_SUBMODEL the
        # 7 rows share ~29 mm of axes, a 4.2 mm pitch, against which the s=40 consensus dot is 2.2 mm
        # — the old s=90 was 3.3 mm and crowded its neighbours' rows, burying the sub-model dots.
        # Re-check these if the footprint shrinks again.
        self.ax.scatter(sub["auroc"], [idx[p] for p in sub["pathogen"]], color=hue("periwinkle"),
                        s=8, alpha=0.6, edgecolors="none", zorder=2)
        self.ax.scatter(con["auroc"], [idx[p] for p in con["pathogen"]], color=hue("crimson"),
                        s=40, alpha=0.85, edgecolor="white", linewidth=0.5, zorder=3)
        self.ref_line(0.5, axis="x")
        self.ax.set_yticks(range(len(order)))
        self.ax.set_yticklabels([abbrev(p) for p in order])
        self.ax.set_xlim(0.2, 1)
        if legend:
            marker_legend(self.ax, SUBMODEL_KEY_ENTRIES, loc=self._LEGEND_LOC,
                          bbox_to_anchor=self._LEGEND_ANCHOR, borderpad=0.3,
                          handletextpad=0.4, labelspacing=self._LEGEND_LABELSPACING)
        self.label(xlabel="own-assay AUROC", title="")


# --------------------------------------------------------------------------- #
# Standalone key panels                                                        #
# --------------------------------------------------------------------------- #
# At SMALL_SQUARE none of these panels can host their own legend, so each palette is emitted once as
# a key panel to place by hand. Keys and panels share these definitions, so a key can never drift
# from the marks it explains.

#: Marker entries for the sub-model vs consensus key (SubmodelAurocSummaryPlot). Markersizes are the
#: sqrt of the panel's scatter areas, so the swatches are true to the marks they explain.
#: Periwinkle for the sub-models, matching step 03's member dots — see the class docstring.
SUBMODEL_KEY_ENTRIES = [
    {"label": "sub-model", "color": hue("periwinkle"), "markersize": 8 ** 0.5},
    {"label": "consensus", "color": hue("crimson"), "markersize": 40 ** 0.5},
]


def euos_overlap_handles():
    """Patch handles for the EU OpenScreen overlap key: colour = leakage, hatch = which axis.

    The hatch is taken from :class:`EuosOverlapTwinPlot`, stroke width included, so the key cannot
    end up showing a coarser texture than the bars it explains.
    """
    novel, in_training = LEAKAGE_COLORS["novel"], LEAKAGE_COLORS["in_training"]
    library = Patch(facecolor=novel, edgecolor="white",
                    hatch=EuosOverlapTwinPlot.HATCH, label="library")
    library.set_hatch_linewidth(EuosOverlapTwinPlot.HATCH_WIDTH)
    return [
        Patch(facecolor=novel, label="novel to model"),
        Patch(facecolor=in_training, label="in training set"),
        library,
        Patch(facecolor=novel, label="actives"),
    ]


class KeyPanel(BasePlot):
    """A legend with no chart — one palette, drawn centred on a blank axis.

    Follows the ``curation_outcome_legend`` pattern from script 02: where several panels share a
    palette, or a panel is too small to host its own legend, the key is rendered once at a known size
    and placed by hand rather than repeated per panel.

    Parameters
    ----------
    name    : output file stem (always ends ``_key`` so keys sort together).
    mapping : ``{label: colour}`` for a swatch key, or ``None`` when ``handles`` / ``entries`` is used.
    handles : explicit matplotlib handles (for hatches and other non-swatch marks).
    entries : :func:`plotting_utils.marker_legend` entries, for point-marker keys.
    ncol    : columns; a wide strip reads better than a tall column above a row of panels.
    cells   : footprint, defaulting to a 45 x 15 mm strip.
    """

    def __init__(self, name, mapping=None, handles=None, entries=None, ncol=1,
                 ax=None, cells=(0.5, 1.5)):
        BasePlot.__init__(self, ax=ax, cells=cells)
        self.name = name
        if handles is not None:
            self.ax.legend(handles=handles, loc="center", ncol=ncol, fontsize=5, **LEGEND_KW)
        elif entries is not None:
            marker_legend(self.ax, entries, loc="center", ncol=ncol)
        else:
            self.legend(mapping, loc="center", ncol=ncol)
        self.ax.set_axis_off()


# --------------------------------------------------------------------------- #
# Entry point                                                                  #
# --------------------------------------------------------------------------- #
def _read(output_dir, fname):
    path = os.path.join(output_dir, fname)
    return pd.read_csv(path) if os.path.exists(path) else pd.DataFrame()


def _step03_pathogen_colors(output_dir):
    """The step-03 ``pathogen_activity_ratios`` hue per pathogen, or ``None`` if step 03 has not run.

    Reads the *dataset counts* out of step 03's small ``dataset_sizes.csv`` summary and hands them to
    :func:`plotting_colors.pathogen_activity_colors`, which reproduces that panel's positional
    ranking. A sibling-output dependency, and the only one in this module — it is what makes the two
    steps' panels share a colour code. Absent or unreadable, the caller falls back to a single hue.

    Takes that function's default ACCENT level, not step 03's pale ``(0.62, 0.38)`` dot fills. The
    pale level was tried on the rank panel's bar segments and reverted: same hue either way, but the
    saturated level is what every step-05 panel now uses, so one pathogen reads as one colour across
    the step rather than one colour at two weights.
    """
    path = os.path.join(output_dir, "..", "03_chembl_models_performance", "dataset_sizes.csv")
    if not os.path.exists(path):
        return None
    d = pd.read_csv(path)
    if d.empty or "pathogen" not in d.columns:
        return None
    return pathogen_activity_colors(d.groupby("pathogen").size().to_dict())


def save_individual_performance_figures(indiv_dir):
    """Build the per-pathogen sub-model panels (AUROC spread + score-correlation heatmap) from
    the CSVs in ``indiv_dir`` and record their footprints in ``figure_cells.json``."""
    auroc = _read(indiv_dir, "05_submodel_auroc.csv")
    footprints = {}
    if not auroc.empty:
        for code, grp in auroc.groupby("code"):
            p = SubmodelAurocPlot(grp)
            if p.is_available:
                p.save(indiv_dir)
                footprints[p.name] = list(p.cells)
                print(f"  figure: individual_performance/{p.name}")
        for code in auroc["code"].unique():
            path = os.path.join(indiv_dir, f"{code}_submodel_corr.csv")
            if not os.path.exists(path):
                continue
            corr = pd.read_csv(path, index_col=0)
            pathogen = auroc[auroc["code"] == code]["pathogen"].iloc[0]
            p = SubmodelCorrPlot(corr, code, pathogen)
            if p.is_available:
                p.save(indiv_dir)
                footprints[p.name] = list(p.cells)
                print(f"  figure: individual_performance/{p.name}")
    with open(os.path.join(indiv_dir, "figure_cells.json"), "w") as f:
        json.dump(footprints, f, indent=2)


def save_euopenscreen_figures(output_dir):
    """Build every EU OpenScreen panel from the step-05 summary CSVs and save the available ones.

    Panels are filed the same three ways as the CSVs that feed them (see
    :func:`eval_euopenscreen.run_all`): the ``full/`` analysis keeps the compounds the models were
    trained on, ``deduplicated/`` removes them, and the top level holds the panels with no leakage
    dimension (label-only). Each of the three dirs gets its own ``png/``, ``pdf/`` and
    ``figure_cells.json``. Panels with no data are skipped (logged).

    Note the AUROC-family panels live under ``deduplicated/`` because they plot the leakage-filtered
    values (``EuosRocGridPlot(set_name="dedup")``, ``_prefer`` picking dedup, and the dedup-derived
    specificity index) — their ``full/`` counterparts exist as data, not as figures.
    """
    full_dir = os.path.join(output_dir, FULL_SUBDIR)
    dedup_dir = os.path.join(output_dir, DEDUP_SUBDIR)

    # no leakage dimension → top level
    leak = _read(output_dir, "05_leakage_report.csv")
    overlap = _read(output_dir, "05_active_overlap.csv")
    promiscuity = _read(output_dir, "05_hit_promiscuity.csv")
    # full: training-set compounds kept
    sum_stats = _read(full_dir, "05_consensus_sum_boxstats.csv")
    sum_actives = _read(full_dir, "05_consensus_sum_actives.csv")
    max_stats = _read(full_dir, "05_consensus_max_boxstats.csv")
    max_actives = _read(full_dir, "05_consensus_max_actives.csv")
    maxp_stats = _read(full_dir, "05_consensus_max_percentile_boxstats.csv")
    maxp_actives = _read(full_dir, "05_consensus_max_percentile_actives.csv")
    model_rank = _read(full_dir, "05_exclusive_hit_model_rank.csv")
    rank_compounds = _read(full_dir, "05_exclusive_hit_model_rank_compounds.csv")
    # deduplicated: training-set compounds removed
    own = _read(dedup_dir, "05_euopenscreen_auroc.csv")
    roc = _read(dedup_dir, "05_euopenscreen_roc.csv")
    excl = _read(dedup_dir, "05_hit_exclusivity.csv")
    excl_pct = _read(dedup_dir, "05_hit_exclusivity_percentiles.csv")
    maxd_roc = _read(dedup_dir, "05_consensus_max_percentile_dedup_roc.csv")
    sumd_stats = _read(dedup_dir, "05_consensus_sum_dedup_boxstats.csv")
    sumd_actives = _read(dedup_dir, "05_consensus_sum_dedup_actives.csv")
    max_roc = _read(full_dir, "05_consensus_max_roc.csv")
    maxp_roc = _read(full_dir, "05_consensus_max_percentile_roc.csv")
    cross = _read(dedup_dir, "05_cross_organism_euos.csv")
    spec = _read(dedup_dir, "05_specificity_index.csv")
    maxd_stats = _read(dedup_dir, "05_consensus_max_percentile_dedup_boxstats.csv")
    maxd_actives = _read(dedup_dir, "05_consensus_max_percentile_dedup_actives.csv")
    model_rank_dedup = _read(dedup_dir, "05_exclusive_hit_model_rank_dedup.csv")
    # code -> binomial name, for the per-pathogen legend of the rank panels
    rank_names = dict(zip(rank_compounds["code"], rank_compounds["pathogen"])) \
        if not rank_compounds.empty else None
    # its own analysis family, untouched by the full/dedup split
    submodel = _read(os.path.join(output_dir, "individual_performance"),
                     "05_submodel_auroc.csv")

    # The six panels destined for the paper's two 45 mm rows are built at SMALL_SQUARE with
    # `legend=False`; their keys follow as standalone `*_key` panels. Everything else keeps the
    # roomier footprint it had, since those are inspection figures rather than paper panels.
    # Built ahead of the group list because its key needs the per-organism totals the plot computes,
    # so the key cannot be written out by hand without risking a drift from the segments.
    # One pathogen -> one colour across BOTH steps: the step-03 `pathogen_activity_ratios` hues drive
    # every step-05 panel that breaks a total down by pathogen. Resolved once here so the panels and
    # their keys cannot end up on different palettes.
    pathogen_colors = _step03_pathogen_colors(output_dir)
    # ONE row order for every panel with a pathogen y axis, taken from `euos_overlap` (fewest EU
    # OpenScreen actives at the bottom). Resolved once here so no two panels can end up on different
    # sequences — see `_apply_row_order` for why the per-panel metric sorts were given up.
    row_order = overlap_row_order(leak)
    rank_dedup = ExclusiveHitModelRankPlot(model_rank_dedup, ranking="percentile",
                                           name_by_code=rank_names, dedup=True,
                                           pathogen_colors=pathogen_colors,
                                           cells=WIDE_RANK, legend=False)
    rank_key = [
        KeyPanel("exclusive_hit_model_rank_key",
                 {f"{abbrev(c, rank_names)} ({n})": rank_dedup.colors[c]
                  for c, n in rank_dedup.segment_totals.items()}, ncol=2, cells=(0.9, 1.5)),
    ] if getattr(rank_dedup, "segment_totals", None) else []

    overlap_pies = ActiveOverlapPiePlot(overlap, row_colors=pathogen_colors, cells=SMALL_SQUARE)
    overlap_pie_key = [ActiveOverlapKeyPlot()] if overlap_pies.is_available else []

    groups = [
        (output_dir, "", [
            EuosOverlapTwinPlot(leak, cells=WIDE_OVERLAP, legend=False),
            KeyPanel("euos_overlap_key", handles=euos_overlap_handles(), ncol=2,
                     cells=(0.55, 1.5)),
            # Both views of the same active sets, kept side by side: the symmetric Jaccard heatmap
            # and the directional containment pie matrix.
            ActiveOverlapHeatmapPlot(overlap, cells=(4, 4)),
            overlap_pies, *overlap_pie_key,
            HitPromiscuityPlot(promiscuity, cells=WIDE_PROMISCUITY),
        ]),
        (full_dir, FULL_SUBDIR, [
            ConsensusSumByHitClassPlot(sum_stats, sum_actives, cells=(3, 3)),
            ConsensusMaxByActivityPlot(max_stats, max_actives, cells=(3, 3)),
            ConsensusMaxByActivityPlot(maxp_stats, maxp_actives, normalized=True, cells=(3, 3)),
            # The same separation as a retrieval curve; see ConsensusMaxRocPlot.
            ConsensusMaxRocPlot(max_roc, cells=(3, 3)),
            ConsensusMaxRocPlot(maxp_roc, normalized=True, cells=(3, 3)),
            ExclusiveHitModelRankPlot(model_rank, ranking="raw", name_by_code=rank_names,
                                      pathogen_colors=pathogen_colors, cells=(3, 3)),
            ExclusiveHitModelRankPlot(model_rank, ranking="percentile", name_by_code=rank_names,
                                      pathogen_colors=pathogen_colors, cells=(3, 3)),
        ]),
        (dedup_dir, DEDUP_SUBDIR, [
            EuosRocGridPlot(roc, set_name="dedup", cols=3),
            EuosSharedEnrichmentPlot(own, cells=(3, 3), row_order=row_order),
            HitExclusivityPlot(excl, cells=EXCLUSIVITY_BARS, legend=False,
                               row_order=row_order),
            # The distribution behind the bars above: one tick per hit at its rank percentile,
            # with each lane's AUROC printed on the right. Kept ALONGSIDE the AUROC panel rather
            # than replacing it — the two share the exclusive/shared key, so choose one in
            # Illustrator and drop the other, as with pathogen_circles vs pathogen_voronoi.
            HitExclusivityEventPlot(excl_pct, excl, leak_df=leak, cells=WIDE_EVENTS),
            KeyPanel("hit_exclusivity_auroc_key",
                     {_pretty(g): EXCLUSIVITY_COLORS[g]
                      for g in ("exclusive", "nonexclusive")}, ncol=2, cells=(0.4, 1.5)),
            CrossOrganismHeatmapPlot(cross, cells=(4, 4)),
            SpecificityIndexPlot(spec, cells=(3, 3), row_order=row_order),
            # Styled as step 03's twin (see the class docstring) but WITHOUT its in-axes key: at
            # WIDE_SUBMODEL the axes is 21.7 mm wide and the key measures 14.4 mm against 7.9 mm of
            # clear corner, so it stays a standalone panel. Both read SUBMODEL_KEY_ENTRIES, so the
            # key cannot drift from the marks.
            SubmodelAurocSummaryPlot(submodel, cells=WIDE_SUBMODEL, legend=False,
                                     row_order=row_order),
            KeyPanel("submodel_auroc_summary_key", entries=SUBMODEL_KEY_ENTRIES, ncol=2,
                     cells=(0.4, 1.5)),
            # Matched to hit_exclusivity_events on AXES HEIGHT (see BOX_EVENTS), which forces the
            # class counts out of the tick labels — see ScoreByHitClassPlot. Its scale goes on the
            # RIGHT spine, away from the event panel's AUROC column. The full/ twin keeps the
            # counts, the left-hand scale and its 90 mm footprint.
            ConsensusSumByHitClassPlot(sumd_stats, sumd_actives, dedup=True,
                                       cells=BOX_EVENTS, show_counts=False,
                                       label_rotation=45, ylabel_right=True),
            ConsensusMaxRocPlot(maxd_roc, normalized=True, dedup=True, cells=SMALL_SQUARE),
            ConsensusMaxByActivityPlot(maxd_stats, maxd_actives, normalized=True, dedup=True,
                                       cells=BOX_NARROW),
            rank_dedup, *rank_key,
        ]),
    ]
    for target_dir, label, plots in groups:
        os.makedirs(target_dir, exist_ok=True)
        _save_group(plots, target_dir, label)


def _save_group(plots, target_dir, label):
    """Save a group of panels into ``target_dir`` and write that dir's ``figure_cells.json``."""
    prefix = f"{label}/" if label else ""
    footprints = {}
    for p in plots:
        if p.is_available:
            p.save(target_dir)
            footprints[p.name] = list(p.cells)
            print(f"  figure: {prefix}{p.name}")
        else:
            print(f"  [skip figure] {prefix}{p.name}: no data")

    with open(os.path.join(target_dir, "figure_cells.json"), "w") as f:
        json.dump(footprints, f, indent=2)
