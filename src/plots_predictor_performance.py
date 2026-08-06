"""Step 15 figures — how well each property column predicts pathogen activity.

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

import json
import os

import numpy as np
import pandas as pd

from default import (CURATED_FAMILY_HUES, PREDICTOR_CHANCE_LEVEL, PREDICTOR_FAMILIES,
                     PREDICTOR_METRICS, RANDOM_SEED)
from plotting_base import BasePlot
from plotting_colors import INK, distinct_colors, hue
from plotting_utils import _jitter_points, abbrev, box_with_jitter

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


class PredictorPerformancePlot(BasePlot):
    """One family's predictors, each a box over its performance across all activity targets."""

    def __init__(self, family, perf, ax=None, cells=None):
        super().__init__(ax=ax, cells=cells or FAMILY_CELLS.get(family, DEFAULT_CELLS))
        self.name = f"13_performance_{family}"

        block = perf[perf["family"] == family]
        if not len(block):
            self._unavailable()
            return

        # Sort by median performance so the figure reads as a ranking, not as config file order.
        order = (block.groupby("predictor")["value"].median()
                 .sort_values(ascending=False).index.tolist())

        types = block.drop_duplicates("predictor").set_index("predictor")["predictor_type"]
        for i, predictor in enumerate(order):
            vals = block.loc[block["predictor"] == predictor, "value"]
            box_with_jitter(
                self.ax, vals, i, TYPE_COLORS[types[predictor]],
                width=BOX_WIDTH, filled=False,
                point_size=POINT_SIZE, point_alpha=POINT_ALPHA)

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


#: Pair kind -> colour for the activity-vs-activity figure. Two endpoints of the SAME pathogen
#: should agree; separating them is what keeps "this model is self-consistent" from being read as
#: "this model predicts everything".
PAIR_COLORS = {True: hue("crimson"), False: hue("cobalt")}
PAIR_LABELS = {True: "same organism", False: "cross organism"}

#: Points drawn per organism box before subsampling. Each box pools (its endpoints x 259 targets),
#: which reaches ~16,500 points for P. falciparum — far past what a 3 mm box can show.
SELF_POINT_CAP = 400
SELF_JITTER_WIDTH = 0.22


class ActivitySelfPerformancePlot(BasePlot):
    """Activity endpoints predicting each other, grouped by the predictor endpoint's organism.

    One box per organism over all its endpoints' AUROCs against every OTHER endpoint's binarized
    version, with points coloured by whether the target endpoint belongs to the same organism.

    Self-pairs are excluded: an endpoint against its own binarization is 1.0 by construction and
    would add one guaranteed-perfect point to every box.
    """

    def __init__(self, self_perf, ax=None, cells=(2, 6), seed=RANDOM_SEED):
        super().__init__(ax=ax, cells=cells)
        self.name = "13_performance_activity_by_organism"

        block = self_perf[~self_perf["self_pair"]]
        if not len(block):
            self._unavailable()
            return

        rng = np.random.default_rng(seed)  # jitter + point subsampling are stochastic
        order = (block.groupby("predictor_organism")["value"].median()
                 .sort_values(ascending=False).index.tolist())

        for i, organism in enumerate(order):
            g = block[block["predictor_organism"] == organism]
            # Neutral box over ALL of this organism's pairs; the two colours ride on the points, so
            # the box must not claim either of them.
            box_with_jitter(self.ax, g["value"], i, INK, width=BOX_WIDTH, filled=False,
                            jitter=False)
            for same in (False, True):
                vals = g.loc[g["same_organism"] == same, "value"].to_numpy(dtype=float)
                vals = vals[np.isfinite(vals)]
                if not len(vals):
                    continue
                _jitter_points(self.ax, vals, i, PAIR_COLORS[same], vert=True,
                               jitter_width=SELF_JITTER_WIDTH, cap=SELF_POINT_CAP,
                               point_size=POINT_SIZE, point_alpha=POINT_ALPHA, rng=rng)

        self.ax.set_xticks(range(len(order)))
        self.ax.set_xticklabels([abbrev(o) for o in order], rotation=90, ha="center")
        self.ax.set_xlim(-0.8, len(order) - 0.2)
        self.ref_line(PREDICTOR_CHANCE_LEVEL, axis="y")
        self.label(ylabel=f"{PREDICTOR_METRICS['continuous'].upper()}")
        self.legend({PAIR_LABELS[k]: PAIR_COLORS[k] for k in (True, False)})


class CuratedPredictorPlot(BasePlot):
    """The curated 12 property predictors against the pathogen-subset activity targets, one panel.

    All three families share one axis because every curated predictor is continuous, so the whole
    panel is AUROC on a single scale — unlike the per-family figures, where the abx block mixes
    AUROC and balanced accuracy and the two could not be pooled.

    With only three groups, colour is a genuine primary encoding (a 3-entry legend reads at panel
    size), and families are kept in contiguous blocks so the comparison a reader makes first is
    within-family. Within each family, predictors are ordered by median.
    """

    def __init__(self, curated, ax=None, cells=(3, 4), family_hues=CURATED_FAMILY_HUES):
        super().__init__(ax=ax, cells=cells)
        self.name = "13_performance_curated_predictors"

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
                            point_size=POINT_SIZE * 2, point_alpha=0.35)

        self.ax.set_xticks(range(len(order)))
        self.ax.set_xticklabels([p.split("__")[-1] for p in order], rotation=90, ha="center")
        self.ax.set_xlim(-0.8, len(order) - 0.2)
        self.ref_line(PREDICTOR_CHANCE_LEVEL, axis="y")
        self.label(ylabel=f"{PREDICTOR_METRICS['continuous'].upper()}")
        self.legend({f: hue(h) for f, h in family_hues.items()})


class PathogenSubsetPerformancePlot(BasePlot):
    """One box per activity endpoint of the 15 pathogens of interest, grouped by pathogen.

    The consensus collapse (see ``eval_predictor_performance.pathogen_subset_endpoints``) cuts the
    endpoint count far enough that each one gets its own box, which the pooled per-organism figure
    cannot show: a pathogen's consensus score and its individual assay endpoints sit side by side.

    Boxes are ordered by pathogen (pathogens by median, endpoints by median within), and coloured by
    pathogen as a SECONDARY cue only — every box is identified on the axis as
    ``{pathogen} - {endpoint}``, so nothing depends on telling 15 hues apart, which
    ``plotting_colors.distinct_colors`` explicitly warns against as a sole encoding.

    Self-pairs are excluded for the same reason as the full figure: 1.0 by construction.
    """

    #: Taller than the other step-13 panels (4 rows, not 2): the tick labels carry both the pathogen
    #: and the endpoint, and at 59 categories they need roughly half the panel height on their own.
    #: At (2, 6) the boxes were squeezed into a thin strip along the top.
    def __init__(self, subset, ax=None, cells=(4, 6), seed=RANDOM_SEED):
        super().__init__(ax=ax, cells=cells)
        self.name = "13_performance_pathogen_subset"

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
                               jitter_width=SELF_JITTER_WIDTH, cap=SELF_POINT_CAP,
                               point_size=POINT_SIZE * (3 if same else 1),
                               point_alpha=0.75 if same else POINT_ALPHA, rng=rng)

        labels = [f"{abbrev(p)} - {e.split(':')[-1]}" for p, e in order]
        self.ax.set_xticks(range(len(order)))
        self.ax.set_xticklabels(labels, rotation=90, ha="center")
        self.ax.set_xlim(-0.8, len(order) - 0.2)
        self.ref_line(PREDICTOR_CHANCE_LEVEL, axis="y")
        self.label(ylabel=f"{PREDICTOR_METRICS['continuous'].upper()}")


# --------------------------------------------------------------------------- #
# Entry point                                                                  #
# --------------------------------------------------------------------------- #
def save_curated_predictor_figure(output_dir, curated=None, footprints=None):
    """The curated 12-predictor figure, appended to ``figure_cells.json``."""
    if curated is None:
        curated = pd.read_csv(os.path.join(output_dir, "13_curated_predictor_performance.csv"))
    plot = CuratedPredictorPlot(curated)
    footprints = {} if footprints is None else footprints
    if plot.is_available:
        plot.save(output_dir)
        footprints[plot.name] = list(plot.cells)
        print(f"  figure: {plot.name}")
    else:
        print("  [skip figure] 13_performance_curated_predictors: no rows")
    with open(os.path.join(output_dir, "figure_cells.json"), "w") as f:
        json.dump(footprints, f, indent=2)
    return footprints


def save_pathogen_subset_figure(output_dir, subset=None, footprints=None):
    """The 15-pathogen, consensus-collapsed figure, appended to ``figure_cells.json``."""
    if subset is None:
        subset = pd.read_csv(os.path.join(output_dir, "13_pathogen_subset_self_performance.csv"))
    plot = PathogenSubsetPerformancePlot(subset)
    footprints = {} if footprints is None else footprints
    if plot.is_available:
        plot.save(output_dir)
        footprints[plot.name] = list(plot.cells)
        print(f"  figure: {plot.name}")
    else:
        print("  [skip figure] 13_performance_pathogen_subset: no rows")
    with open(os.path.join(output_dir, "figure_cells.json"), "w") as f:
        json.dump(footprints, f, indent=2)
    return footprints


def save_activity_self_figure(output_dir, self_perf=None, footprints=None):
    """The activity-vs-activity figure, appended to ``figure_cells.json``."""
    if self_perf is None:
        self_perf = pd.read_csv(os.path.join(output_dir, "13_activity_self_performance.csv"))
    plot = ActivitySelfPerformancePlot(self_perf)
    footprints = {} if footprints is None else footprints
    if plot.is_available:
        plot.save(output_dir)
        footprints[plot.name] = list(plot.cells)
        print(f"  figure: {plot.name}")
    else:
        print("  [skip figure] 13_performance_activity_by_organism: no rows")
    with open(os.path.join(output_dir, "figure_cells.json"), "w") as f:
        json.dump(footprints, f, indent=2)
    return footprints


def save_predictor_performance_figures(output_dir, perf=None, families=PREDICTOR_FAMILIES):
    """Build one figure per predictor family and record their footprints in ``figure_cells.json``."""
    if perf is None:
        perf = pd.read_csv(os.path.join(output_dir, "13_predictor_performance.csv"))

    footprints = {}
    for family in families:
        plot = PredictorPerformancePlot(family, perf)
        if plot.is_available:
            plot.save(output_dir)
            footprints[plot.name] = list(plot.cells)
            print(f"  figure: {plot.name}")
        else:
            print(f"  [skip figure] 13_performance_{family}: no rows")
    with open(os.path.join(output_dir, "figure_cells.json"), "w") as f:
        json.dump(footprints, f, indent=2)
    # Returned so save_activity_self_figure can append to the same manifest rather than
    # overwriting it with a single-entry one.
    return footprints
