"""Step 14 figures — how well each property column predicts pathogen activity.

Reads ONLY the long per-pair CSV written by :mod:`eval_predictor_performance`. One figure per
predictor family (:data:`default.PREDICTOR_FAMILIES`): a box per predictor over its distribution of
performance across all activity targets, jittered points behind it, predictors sorted by median.

101 predictors on one shared x-axis would leave ~1.8 mm per category at the 180 mm print width,
hence the split by family.

Box colour encodes the predictor's VALUE TYPE, which is also what selected its metric — continuous
predictors are scored by AUROC and binary ones by balanced accuracy. Both share a 0.5 chance
baseline (drawn), which is the only reason the two can share a y-axis; the colour is what stops the
`abx` figure, where both types occur, from reading as one homogeneous measure.
"""

import os

import numpy as np
import pandas as pd

from default import (CURATED_FAMILY_HUES, PREDICTOR_CHANCE_LEVEL, PREDICTOR_FAMILIES,
                     PREDICTOR_METRICS, RANDOM_SEED)
from plotting_base import BasePlot
from plotting_colors import hue
from plotting_utils import box_with_jitter, merge_figure_cells

#: Predictor value type -> colour. Two clearly separated hues, since this doubles as the metric key.
TYPE_COLORS = {"continuous": hue("cobalt"), "binary": hue("amber")}

#: Type -> the legend wording, naming the metric rather than the type: the metric is what a reader
#: needs in order to know what the y-axis means for that box.
TYPE_LABELS = {
    "continuous": f"continuous predictor ({PREDICTOR_METRICS['continuous'].upper()})",
    "binary": "binary predictor (balanced accuracy)",
}

#: Cells (rows, cols) per family, sized by category count at ~3 mm minimum per box.
FAMILY_CELLS = {"physchem": (2, 4), "cytotox": (2, 4), "abx": (2, 6)}
DEFAULT_CELLS = (2, 6)

POINT_SIZE = 1.5
POINT_ALPHA = 0.25
BOX_WIDTH = 0.6
#: Half-width of the point scatter around each box centre, so the swarm spans 2x this. 0.22 against
#: a 0.6 box keeps every point inside its own box and well clear of the neighbouring position (1.0
#: apart). Shared by all four figure families so the panels look alike.
#: NOTE: `box_with_jitter` silently draws NO jitter unless an `rng` is also passed —
#: `_jitter_points` falls back to `np.zeros` when `rng is None`, so `jitter_width` alone does
#: nothing. Every call below therefore passes both.
JITTER_WIDTH = 0.22


class PredictorPerformancePlot(BasePlot):
    """One family's predictors, each a box over its performance across all activity targets."""

    def __init__(self, family, perf, ax=None, cells=None, seed=RANDOM_SEED):
        super().__init__(ax=ax, cells=cells or FAMILY_CELLS.get(family, DEFAULT_CELLS))
        self.name = f"14_performance_{family}"

        block = perf[perf["family"] == family]
        if not len(block):
            self._unavailable()
            return

        # Sort by median performance so the figure reads as a ranking, not as config file order.
        order = (block.groupby("predictor")["value"].median()
                 .sort_values(ascending=False).index.tolist())

        rng = np.random.default_rng(seed)  # jitter is stochastic; seeded so the figure is stable
        types = block.drop_duplicates("predictor").set_index("predictor")["predictor_type"]
        for i, predictor in enumerate(order):
            vals = block.loc[block["predictor"] == predictor, "value"]
            box_with_jitter(
                self.ax, vals, i, TYPE_COLORS[types[predictor]],
                width=BOX_WIDTH, filled=False, jitter_width=JITTER_WIDTH,
                point_size=POINT_SIZE, point_alpha=POINT_ALPHA, rng=rng)

        self.ax.set_xticks(range(len(order)))
        # Only the column name is shown: the family is the figure, and the model ID is in the CSV.
        self.ax.set_xticklabels([p.split("__")[-1] for p in order],
                                rotation=90, ha="center")
        self.ax.set_xlim(-0.8, len(order) - 0.2)
        self.ref_line(PREDICTOR_CHANCE_LEVEL, axis="y")
        self.label(ylabel="Performance")

        present = [t for t in ("continuous", "binary") if (types == t).any()]
        if len(present) > 1:
            self.legend({TYPE_LABELS[t]: TYPE_COLORS[t] for t in present})


class CuratedPredictorPlot(BasePlot):
    """The curated 12 property predictors against the pathogen-subset activity targets, one panel.

    All three families share one axis because every curated predictor is continuous, so the whole
    panel is AUROC on a single scale — unlike the per-family figures, where the abx block mixes
    AUROC and balanced accuracy and the two could not be pooled.

    With only three groups, colour is a genuine primary encoding (a 3-entry legend reads at panel
    size), and families are kept in contiguous blocks so the comparison a reader makes first is
    within-family. Within each family, predictors are ordered by median.
    """

    def __init__(self, curated, ax=None, cells=(3, 4), family_hues=CURATED_FAMILY_HUES,
                 seed=RANDOM_SEED):
        super().__init__(ax=ax, cells=cells)
        self.name = "14_performance_curated_predictors"
        rng = np.random.default_rng(seed)  # jitter is stochastic; seeded so the figure is stable

        if not len(curated):
            self._unavailable()
            return

        medians = curated.groupby("predictor")["value"].median()
        order, colors = [], {}
        # Families in the order declared in CURATED_PREDICTORS, each block sorted by median.
        for family in family_hues:
            preds = curated.loc[curated["family"] == family, "predictor"].unique().tolist()
            for p in sorted(preds, key=lambda x: -medians[x]):
                order.append(p)
                colors[p] = hue(family_hues[family])

        for i, predictor in enumerate(order):
            vals = curated.loc[curated["predictor"] == predictor, "value"]
            box_with_jitter(self.ax, vals, i, colors[predictor], width=BOX_WIDTH, filled=False,
                            jitter_width=JITTER_WIDTH, point_size=POINT_SIZE * 2,
                            point_alpha=0.35, rng=rng)

        self.ax.set_xticks(range(len(order)))
        self.ax.set_xticklabels([p.split("__")[-1] for p in order], rotation=90, ha="center")
        self.ax.set_xlim(-0.8, len(order) - 0.2)
        self.ref_line(PREDICTOR_CHANCE_LEVEL, axis="y")
        self.label(ylabel=f"{PREDICTOR_METRICS['continuous'].upper()}")
        self.legend({f: hue(h) for f, h in family_hues.items()})


# --------------------------------------------------------------------------- #
# Entry point                                                                  #
# --------------------------------------------------------------------------- #
def save_curated_predictor_figure(output_dir, curated=None):
    """The curated 12-predictor figure, merged into ``figure_cells.json``."""
    if curated is None:
        curated = pd.read_csv(os.path.join(output_dir, "14_curated_predictor_performance.csv"))
    plot = CuratedPredictorPlot(curated)
    new = {}
    if plot.is_available:
        plot.save(output_dir)
        new[plot.name] = list(plot.cells)
        print(f"  figure: {plot.name}")
    else:
        print("  [skip figure] 14_performance_curated_predictors: no rows")
    return merge_figure_cells(output_dir, new)


def save_predictor_performance_figures(output_dir, perf=None, families=PREDICTOR_FAMILIES):
    """Build one figure per predictor family and merge their footprints into ``figure_cells.json``."""
    if perf is None:
        perf = pd.read_csv(os.path.join(output_dir, "14_predictor_performance.csv"))

    new = {}
    for family in families:
        plot = PredictorPerformancePlot(family, perf)
        if plot.is_available:
            plot.save(output_dir)
            new[plot.name] = list(plot.cells)
            print(f"  figure: {plot.name}")
        else:
            print(f"  [skip figure] 14_performance_{family}: no rows")
    return merge_figure_cells(output_dir, new)
