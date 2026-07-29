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
from sklearn.metrics import auc, roc_curve

from plotting_base import BasePlot, GridPlot
from plotting_colors import ACTIVE_INACTIVE_COLORS, BAR_DEFAULT, REFERENCE_LINE, auroc_shades
from plotting_utils import abbrev, box_with_jitter, marker_legend, roc_panel

# Full page width: 6 columns of 3 cm. Chosen over 4 columns so the model-rich
# pathogens (P. falciparum 52 models, M. tuberculosis 34) stay within a
# supplementary-page height while each ROC panel keeps its 3 cm square cell.
ROC_GRID_COLS = 6

# Vertical density of the rank boxplots: how many model rows (= one actives +
# one inactives box) share a single 3 cm cell.
MODELS_PER_CELL = 4


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
