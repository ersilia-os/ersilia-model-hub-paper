"""Publication panels for the ChEMBL antimicrobial model performance (script 03).

Two panel types, one pair per pathogen:

- :class:`RocGridPlot` — small multiples, one ROC curve per model, laid out on the
  3 cm cell grid at 6 columns (full 180 mm page width) so every ROC panel is a true
  square. Rows grow with the model count; trailing empty cells are hidden.
- :class:`RankBoxplotPlot` — paired horizontal boxplots of the out-of-fold predicted
  rank of actives vs inactives, one pair per model.

Both follow ``docs/figure_conventions.md``: ``BasePlot`` subclasses, sized in cells,
saved as PNG + vector PDF, colours only from :mod:`plotting_colors`. The ROC grid is built
on the shared :class:`GridPlot` scaffold and the shared ``roc_panel`` / ``box_with_jitter``
primitives, so it reads identically to the step-05 ROC grid.
"""

import json
import math
import os

import numpy as np
import stylia
from sklearn.metrics import auc, roc_curve

from plotting_base import BasePlot, GridPlot
from plotting_colors import (ACTIVE_INACTIVE_COLORS, BAR_DEFAULT, REFERENCE_LINE, auroc_shades,
                             distinct_colors, hue)
from plotting_utils import (abbrev, box_with_jitter, marker_legend, nested_size_legend, pie_scatter,
                            ref_line, roc_panel, swatch_legend)
from default import RANDOM_SEED

# Full page width: 6 columns of 3 cm. Chosen over 4 columns so the model-rich
# pathogens (P. falciparum 52 models, M. tuberculosis 34) stay within a
# supplementary-page height while each ROC panel keeps its 3 cm square cell.
ROC_GRID_COLS = 6

# Vertical density of the rank boxplots: how many model rows (= one actives +
# one inactives box) share a single 3 cm cell.
MODELS_PER_CELL = 4


class PathogenDatasetSizesPlot(BasePlot):
    """Size of every modelled dataset, per pathogen, with its active/inactive split.

    One box per pathogen over the sizes of its datasets, and one small **pie** per dataset drawn on
    top: the pie's filled share is that dataset's active fraction. So the box gives the pathogen's
    size distribution and each circle gives one dataset's balance, in the same panel.

    **Added negatives are excluded**, from both the size and the ratio. Airtable's ``n_compounds``
    counts them, so size is ``n_compounds - n_added_negatives - n_added_decoys`` and the active
    fraction is ``n_positives / size``; that reproduces the curation pipeline's own ``n_mol_after``
    and ``ar_after`` exactly (checked for all 151 models that join to ``25_pool_summary``). It matters
    for the 54 of 193 datasets that got negatives added: ``mtuberculosis/DR_0012`` goes from 2450
    compounds at a 0.50 active ratio to 1411 at 0.87. Decoys were never used — the column is zero for
    every model — but it is subtracted anyway so the definition does not silently depend on that.

    Pie diameter is FIXED. The y position already encodes size on a log axis, so scaling the circles
    would double-encode it, and with a 96-to-334,766 range the largest would swamp the panel while
    the smallest vanished.

    Full page width (180 x 90 mm) is a legibility requirement, not a layout preference: 15 pathogens
    share the x axis and *P. falciparum* alone holds 51 datasets, so the columns need ~12 mm each for
    the pies to separate under horizontal jitter. At half width they overlap into a smear.

    A two-slice pie at this size reads as a gestalt, not a measurement — "mostly active" versus
    "mostly inactive" is clear, 45 % versus 50 % is not. Exact ratios are in the summary CSV.

    Parameters
    ----------
    datasets : DataFrame with columns ``pathogen``, ``size`` and ``active_fraction``, one row per
               modelled dataset (see the script's derivation).
    pathogen_names : optional dict code -> full species name, for genus-abbreviated tick labels.
    cells    : footprint on the reference grid as ``(rows, cols)``.
    """

    #: Marker area in points squared. ~2.5 mm across at 180 mm width — big enough to read a slice,
    #: small enough that 51 of them in a 12 mm column still separate under jitter.
    PIE_AREA = 22
    #: Horizontal jitter as a fraction of the column pitch, so dense columns spread out.
    JITTER = 0.30

    def __init__(self, ax=None, datasets=None, pathogen_names=None, cells=(3, 6)):
        BasePlot.__init__(self, ax=ax, cells=cells)
        self.name = "pathogen_dataset_sizes"
        self.is_available = datasets is not None and len(datasets) > 0
        if not self.is_available:
            self._unavailable()
            return
        ax = self.ax
        rng = np.random.default_rng(RANDOM_SEED)

        # Most datasets first, so the eye starts where the evidence is.
        order = (datasets.groupby("pathogen")["size"].size()
                 .sort_values(ascending=False).index.tolist())
        self.counts = {p: int((datasets["pathogen"] == p).sum()) for p in order}

        c_act = ACTIVE_INACTIVE_COLORS["active"]
        c_ina = ACTIVE_INACTIVE_COLORS["inactive"]
        for i, pathogen in enumerate(order):
            g = datasets[datasets["pathogen"] == pathogen]
            sizes = g["size"].to_numpy(dtype=float)
            # Box only — the pies ARE the point overlay, so box_with_jitter's own scatter is off.
            # Pale fill on purpose: the pies sit inside the box, and at full-strength cobalt the
            # crimson/silver slices lose contrast against it.
            box_with_jitter(ax, sizes, i, BAR_DEFAULT, face=hue("cobalt", lighten=0.22),
                            jitter=False, width=0.55)
            jitter = rng.uniform(-self.JITTER, self.JITTER, len(g))
            pie_scatter(ax, i + jitter, sizes, g["active_fraction"].to_numpy(dtype=float),
                        (c_act, c_ina), s=self.PIE_AREA)

        ax.set_yscale("log")
        ax.set_xlim(-0.7, len(order) - 0.3)
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels([f"{abbrev(p, pathogen_names)}\n({self.counts[p]})" for p in order],
                           rotation=30, ha="right", style="italic")
        swatch_legend(ax, {"Active": c_act, "Inactive": c_ina}, loc="upper right")
        self.label(ylabel="Compounds per dataset", title="Modelled dataset sizes")


class PathogenActivityRatiosPlot(BasePlot):
    """Active fraction of every modelled dataset, per pathogen, with dot area encoding its size.

    The complement of :class:`PathogenDatasetSizesPlot`: that panel puts size on the y axis and the
    balance inside each mark, this one swaps them. Balance is the quantity a reader can actually read
    off a position scale here, which is the point — a pie at 2.5 mm cannot distinguish 45 % from 50 %,
    a y position can. Size drops to the secondary encoding.

    A dashed reference line marks the **0.5 balance point**: above it a dataset has more actives than
    inactives, which is true of 54 of the 193 once added negatives are excluded. Sizes and ratios use
    the same added-negative-free definition as the sibling panel.

    **Dot area is affine in sqrt(size), not proportional to size.** Sizes span 96 to 334,766 — a
    3,500x range — so area-proportional dots would put a 59x radius ratio on the panel: the biggest
    would blot out its column and the smallest would be invisible. The sqrt transform is the same
    compromise ``AREA_EXPONENT`` makes in the script-01 pathogen treemap. Read sizes off the legend,
    which is generated by the same function, rather than comparing areas by eye.

    A short horizontal bar per column marks that pathogen's **mean** active fraction, in the
    full-strength hue against the pale dots. Note it is the unweighted mean over datasets, so every
    dataset counts once regardless of size — and because the very large datasets sit near 0 % active,
    a size-weighted mean would fall well below this bar. It answers "what does a typical dataset for
    this pathogen look like", not "what fraction of its compounds are active".

    **Colour is redundant with the x axis** — one hue per pathogen, but the pathogen is already named
    on the tick. That is deliberate: it makes a column scannable and lets a reader track a pathogen
    across panels, and it means the palette does not have to carry 15 unambiguous hues (only 9 exist,
    so 6 are tints — see :func:`plotting_colors.distinct_colors`). There is no colour legend; 15
    swatches would not be readable at panel size and the axis is the real key. Dots take the pale
    ``FILL_LEVELS`` palette so heavy overlap stays readable, outlined and topped by the darker
    ``ACCENT_LEVELS`` of the same hue; both calls index hues identically, so a dot and its mean bar
    read as one category.

    Parameters
    ----------
    datasets : DataFrame with ``pathogen``, ``size`` and ``active_fraction``, one row per dataset.
    pathogen_names : optional dict code -> full species name, for genus-abbreviated tick labels.
    cells    : footprint on the reference grid as ``(rows, cols)``.
    """

    #: Marker area in points squared at ``SIZE_REF`` compounds; every other dot scales as
    #: ``sqrt(size / SIZE_REF)`` from here. 10,000 compounds -> 30 pt^2 puts the 96-compound floor at
    #: ~0.7 mm across and the 334,766-compound ceiling at ~4.9 mm.
    SIZE_REF = 10_000
    SIZE_REF_AREA = 30.0
    #: Decades shown in the size key.
    SIZE_LEGEND_KEYS = (100, 1_000, 10_000, 100_000)
    #: Horizontal jitter as a fraction of the column pitch.
    JITTER = 0.30
    #: Pale palette for the dot fills, and the darker one for their outlines and mean bars. Two
    #: lighten levels each because ``distinct_colors`` uses the second pass to get past the 9-hue
    #: limit — a single level would make pathogen 10 identical to pathogen 1.
    FILL_LEVELS = (0.62, 0.38)
    ACCENT_LEVELS = (None, 0.55)
    #: Headroom above 1 and below 0, as a fraction of the axis. The active fraction is bounded, so the
    #: ticks stop at 0 and 1; this only stops the largest dots (~2.5 mm radius) sitting on the spines.
    Y_PAD = 0.07
    #: Where the nested size key sits, in data coordinates.
    LEGEND_X = 13.6
    LEGEND_Y = 0.80

    def _area(self, size):
        return self.SIZE_REF_AREA * np.sqrt(np.asarray(size, dtype=float) / self.SIZE_REF)

    def __init__(self, ax=None, datasets=None, pathogen_names=None, cells=(3, 6)):
        BasePlot.__init__(self, ax=ax, cells=cells)
        self.name = "pathogen_activity_ratios"
        self.is_available = datasets is not None and len(datasets) > 0
        if not self.is_available:
            self._unavailable()
            return
        ax = self.ax
        rng = np.random.default_rng(RANDOM_SEED)

        order = (datasets.groupby("pathogen")["size"].size()
                 .sort_values(ascending=False).index.tolist())
        self.counts = {p: int((datasets["pathogen"] == p).sum()) for p in order}
        fills = distinct_colors(len(order), levels=self.FILL_LEVELS)
        accents = distinct_colors(len(order), levels=self.ACCENT_LEVELS)
        self.means = {}

        # Balance line first, so it sits under the data.
        ref_line(ax, 0.5, axis="y")
        for i, (pathogen, fill, accent) in enumerate(zip(order, fills, accents)):
            g = datasets[datasets["pathogen"] == pathogen]
            ratios = g["active_fraction"].to_numpy(dtype=float)
            jitter = rng.uniform(-self.JITTER, self.JITTER, len(g))
            ax.scatter(i + jitter, ratios, s=self._area(g["size"]), facecolor=fill,
                       edgecolors=accent, linewidths=0.4, zorder=3)
            self.means[pathogen] = float(np.mean(ratios))
            ax.hlines(self.means[pathogen], i - 0.42, i + 0.42, color=accent, linewidth=1.4,
                      zorder=4)

        ax.set_ylim(-self.Y_PAD, 1 + self.Y_PAD)
        ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
        ax.set_xlim(-0.8, len(order) - 0.2)
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels([f"{abbrev(p, pathogen_names)}\n({self.counts[p]})" for p in order],
                           rotation=30, ha="right", style="italic")
        # Size key as nested circles, generated from the same _area function as the dots. Drawn last:
        # it reads the live transform, so the limits above must already be final.
        nested_size_legend(ax, self.SIZE_LEGEND_KEYS,
                           [self._area(k) for k in self.SIZE_LEGEND_KEYS],
                           x=self.LEGEND_X, y_base=self.LEGEND_Y, color=hue("silver"))
        ax.text(self.LEGEND_X, self.LEGEND_Y - 0.05, "Compounds", ha="center", va="top",
                fontsize=stylia.FONTSIZE_SMALL)
        self.label(ylabel="Active fraction of dataset", title="Dataset activity ratios")


class RocGridPlot(GridPlot):
    """Grid of per-model ROC curves for one pathogen, sorted by mean AUROC.

    Parameters
    ----------
    pathogen : str
        Pathogen code, used for the output file stem.
    models : list of dict
        Each with ``name``, ``mean_auroc``, ``y_true``, ``y_pred``.
    cols : int
        Grid columns; the footprint is ``(ceil(n / cols), cols)`` cells.
    """

    def __init__(self, pathogen, models, cols=ROC_GRID_COLS, pathogen_names=None):
        items = self._items(models)
        self.build_grid(items, cols=cols, name=f"{pathogen}_roc_curves",
                        panel_fn=self._panel, color_fn=lambda it: it["color"],
                        edge_xlabel="FPR", edge_ylabel="TPR")

    @staticmethod
    def _items(models):
        ordered = sorted(models, key=lambda m: m["mean_auroc"], reverse=True)
        colors = auroc_shades([m["mean_auroc"] for m in ordered]) if ordered else []
        items = []
        for m, color in zip(ordered, colors):
            y_true = np.asarray(m["y_true"])
            y_pred = np.asarray(m["y_pred"], dtype=float)
            if len(np.unique(y_true)) < 2:
                items.append(dict(fpr=None, tpr=None, auroc=None, n_pos=0, n_neg=0,
                                  title=m["name"], color=color))
                continue
            fpr, tpr, _ = roc_curve(y_true, y_pred)
            n_pos = int(y_true.sum())
            items.append(dict(fpr=fpr, tpr=tpr, auroc=auc(fpr, tpr), n_pos=n_pos,
                              n_neg=int(len(y_true) - n_pos), title=m["name"], color=color))
        return items

    @staticmethod
    def _panel(ax, item, color, xlabel, ylabel):
        roc_panel(ax, item["fpr"], item["tpr"], item["auroc"], item["n_pos"], item["n_neg"],
                  color, xlabel=xlabel, ylabel=ylabel, title=item["title"])


class RankBoxplotPlot(BasePlot):
    """Paired rank boxplots (actives vs inactives) for every model of one pathogen.

    Ranks are out-of-fold and pooled across the 5 folds, so each compound is scored
    exactly once. Models are ordered by mean AUROC, best at the top.
    """

    def __init__(self, pathogen, models, cols=3, pathogen_names=None):
        n = len(models)
        rows = max(2, math.ceil(n / MODELS_PER_CELL))
        BasePlot.__init__(self, ax=None, cells=(rows, cols))
        self.name = f"{pathogen}_rank_boxplots"
        if n == 0:
            self.cells = (2, cols)
            self._unavailable()
            return

        ax = self.ax
        ordered = sorted(models, key=lambda m: m["mean_auroc"])  # best last = top
        c_act = ACTIVE_INACTIVE_COLORS["active"]
        c_ina = ACTIVE_INACTIVE_COLORS["inactive"]

        for i, m in enumerate(ordered):
            box_with_jitter(ax, m["rank_actives"], i * 2 + 0.3, c_act, vert=False,
                            width=0.5, jitter=False)
            box_with_jitter(ax, m["rank_inactives"], i * 2 - 0.3, c_ina, vert=False,
                            width=0.5, jitter=False)

        ax.set_xlim(0, 1)
        ax.set_ylim(-1.2, (n - 1) * 2 + 1.2)
        ax.set_yticks([i * 2 for i in range(n)])
        ax.set_yticklabels([m["name"] for m in ordered])

        # Above the axes: every row of the plot area is occupied by boxes, so an
        # in-axes legend would sit on top of the worst-performing model.
        self.legend({"actives": c_act, "inactives": c_ina}, loc="lower right",
                    bbox_to_anchor=(1.0, 1.0), ncol=2, handlelength=1.2, borderpad=0,
                    columnspacing=1.0)
        self.label(xlabel="Out-of-fold predicted rank",
                   title=abbrev(pathogen, pathogen_names))


class PathogenConsensusAurocPlot(BasePlot):
    """Condensed cross-pathogen summary of ChEMBL model CV performance.

    One row per pathogen: a small dot per retained model's mean 5-fold CV AUROC, a thin guide
    line across that pathogen's min–max, and a large crimson dot at the **mean** — the
    representative ('consensus') per-pathogen value. The mean dot's **area encodes the number of
    unique molecules** in that pathogen's ChEMBL dataset (area interpolated between a min and max
    marker size so small datasets stay visible; see the size key). Rows sorted by mean (best on
    top); the dashed line marks the step-10 retention floor (mean AUROC 0.7).

    NOTE: this summarises per-model CV AUROCs — it is NOT a per-compound consensus prediction
    (the CV pools share no held-out set, so a true ensemble AUROC is not computable here).
    """

    _SIZE_MIN = 15.0    # marker area (pt^2) for the smallest dataset
    _SIZE_MAX = 300.0   # ... and the largest

    def _area(self, c):
        """Marker area for a molecule count, sqrt-scaled between the size range (perceptual)."""
        if not np.isfinite(c) or self._cmax <= self._cmin:
            return (self._SIZE_MIN + self._SIZE_MAX) / 2
        t = (np.sqrt(c) - np.sqrt(self._cmin)) / (np.sqrt(self._cmax) - np.sqrt(self._cmin))
        t = min(1.0, max(0.0, t))
        return self._SIZE_MIN + t * (self._SIZE_MAX - self._SIZE_MIN)

    def __init__(self, entries, ax=None, cells=(4, 3)):
        # entries: {"code", "pathogen", "aurocs": [...], "n_molecules": int}
        BasePlot.__init__(self, ax=ax, cells=cells)
        self.name = "pathogen_consensus_auroc"
        entries = [e for e in entries if e.get("aurocs")]
        if not entries:
            self._unavailable()
            return
        entries = sorted(entries, key=lambda e: float(np.mean(e["aurocs"])))
        counts = [e.get("n_molecules") for e in entries if e.get("n_molecules")]
        self._cmin = float(min(counts)) if counts else 1.0
        self._cmax = float(max(counts)) if counts else 1.0
        y = np.arange(len(entries))
        crimson = ACTIVE_INACTIVE_COLORS["active"]
        for i, e in enumerate(entries):
            vals = np.asarray(e["aurocs"], dtype=float)
            self.ax.plot([vals.min(), vals.max()], [i, i], color=REFERENCE_LINE,
                         linewidth=0.8, zorder=1)
            self.ax.scatter(vals, np.full(len(vals), i), color=BAR_DEFAULT,
                            s=10, alpha=0.5, zorder=2)
            self.ax.scatter([vals.mean()], [i], color=crimson, alpha=0.75,
                            s=self._area(e.get("n_molecules", np.nan)), zorder=3)
        self.ref_line(0.7, axis="x")
        self.ax.set_yticks(y)
        self.ax.set_yticklabels([abbrev(e["code"], {e["code"]: e["pathogen"]}) for e in entries])
        self.ax.set_xlim(0.6, 1.0)
        self.ax.set_ylim(-0.8, len(entries) - 0.2)

        # Two legends (both over the empty <0.7 AUROC region, so never read as data):
        #   upper-left = colour key (model vs mean); lower-left = molecule-count size key.
        leg1 = marker_legend(
            self.ax,
            [{"label": "model", "color": BAR_DEFAULT}, {"label": "mean", "color": crimson}],
            loc="upper left")
        self.ax.add_artist(leg1)
        if counts:
            # Round reference values spanning the observed range (~150 to ~500k molecules).
            keys = [150, 10_000, 100_000]
            marker_legend(
                self.ax,
                [{"label": f"{kc:,}", "color": crimson,
                  "markersize": float(np.sqrt(self._area(kc)))} for kc in keys],
                loc="lower left", title="unique molecules",
                labelspacing=1.6, borderpad=1.0, handletextpad=1.2)

        self.label(xlabel="cross-validation AUROC", title="")


def save_consensus_figure(entries, output_dir):
    """Build and save the cross-pathogen consensus-AUROC summary; return its footprint dict."""
    p = PathogenConsensusAurocPlot(entries)
    footprints = {}
    if p.is_available:
        p.save(output_dir)
        footprints[p.name] = list(p.cells)
    return footprints


def save_dataset_sizes_figure(datasets, output_dir, pathogen_names=None):
    """Save both per-pathogen dataset panels (size on y, and ratio on y); return their footprints.

    The two are alternatives on the same data, kept side by side so the better one can be chosen at
    layout time: ``pathogen_dataset_sizes`` puts size on the position scale and balance in a pie,
    ``pathogen_activity_ratios`` swaps them.
    """
    footprints = {}
    for p in (PathogenDatasetSizesPlot(datasets=datasets, pathogen_names=pathogen_names),
              PathogenActivityRatiosPlot(datasets=datasets, pathogen_names=pathogen_names)):
        if not p.is_available:
            continue
        p.save(output_dir)
        footprints[p.name] = list(p.cells)
        print(f"[{p.name}] {sum(p.counts.values())} datasets over {len(p.counts)} pathogens; "
              f"densest column {max(p.counts.values())}")
    return footprints


def save_performance_figures(pathogen, models, output_dir, pathogen_names=None):
    """Build and save both panels for one pathogen; return their cell footprints."""
    plots = [
        RocGridPlot(pathogen, models, pathogen_names=pathogen_names),
        RankBoxplotPlot(pathogen, models, pathogen_names=pathogen_names),
    ]
    footprints = {}
    for p in plots:
        if p.is_available:
            p.save(output_dir)
            footprints[p.name] = list(p.cells)
    return footprints


def write_figure_cells(footprints, output_dir):
    """Write the accumulated ``figure_cells.json`` manifest."""
    path = os.path.join(output_dir, "figure_cells.json")
    with open(path, "w") as fh:
        json.dump(footprints, fh, indent=2)
    return path
