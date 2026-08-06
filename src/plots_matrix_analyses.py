"""Figures for the score-matrix analyses (steps 07–08).

Two panels, both downstream of the named full-library score matrix built by step 07:

  - :func:`pathogen_jaccard_figure` (step 08) — per-pathogen top-N Jaccard, same pathogen vs.
    different pathogen. A diagnostic rather than a paper panel, so it is plain matplotlib rather than
    a :class:`plotting_base.BasePlot` on the 3 cm cell grid: the y-axis carries per-pathogen endpoint
    and pair counts that go illegible at page width. It still draws only through
    :mod:`plotting_colors` hues, :func:`plotting_utils.box_with_jitter` and
    :func:`plotting_utils.sentence_case`, under the stylia print/article style, and writes PNG **and**
    PDF — the vector copy is the readable one.
  - :class:`MeanRankDistributionPlot` (step 07) — distribution of each compound's mean percentile rank
    across every selected endpoint. A single distribution fits the page grid comfortably, so this one
    follows the standard publication convention in ``docs/figure_conventions.md``.
"""

import os

import matplotlib.pyplot as plt
import numpy as np
import stylia as st
from matplotlib.patches import Patch

import plotting_base  # noqa: F401  (applies the stylia print/article style on import)
from default import RANDOM_SEED
from plotting_base import BasePlot
from plotting_colors import INK, hue
from plotting_utils import box_with_jitter, sentence_case

# --------------------------------------------------------------------------- #
# Step 08 — per-pathogen Jaccard                                               #
# --------------------------------------------------------------------------- #
#: Rendered points per box. Every box's TRUE pair count is annotated because the cap makes a
#: 15-pair box and a 9,000-pair box look similarly dense otherwise.
POINT_CAP = 150


def pathogen_jaccard_figure(boxes, summary, *, cutoff, matrix_label, name, output_dir):
    """Per-pathogen same/different-pathogen Jaccard figure — one pair of boxes per pathogen.

    Turquoise = every unordered pair of that pathogen's own endpoint columns, crimson = its columns
    against every other pathogen's. Rows ordered by specificity (same-median minus diff-median)
    descending, matching :func:`eval_correlations.pathogen_metric_summary`.

    Pathogens with fewer than the required number of endpoints are expected to have been removed
    upstream by :func:`eval_correlations.multi_column_pathogen_nodes` — a pathogen with one endpoint
    has no same-pathogen pair at all.

    **Linear x-axis, deliberately** (user-directed): top-N Jaccard values are small and bunch near
    zero, so a linear axis compresses them — but it renders exact-zero pairs, which a log axis
    silently drops. Nothing is filtered out of the boxes.
    """
    order = list(summary["pathogen"])
    counts = summary.set_index("pathogen")
    c_same, c_diff = hue("turquoise"), hue("crimson")
    rng = np.random.default_rng(RANDOM_SEED)

    n = len(order)
    fig, ax = plt.subplots(figsize=(7.2, max(3.5, n * 0.42)))
    labels = []
    for i, pathogen in enumerate(order):
        b = boxes[pathogen]
        if len(b["same"]):
            box_with_jitter(ax, b["same"], i + 0.2, c_same, vert=False, width=0.34,
                            jitter_width=0.1, point_size=5, point_alpha=0.45, cap=POINT_CAP,
                            rng=rng)
        if len(b["diff"]):
            box_with_jitter(ax, b["diff"], i - 0.2, c_diff, vert=False, width=0.34,
                            jitter_width=0.1, point_size=5, point_alpha=0.45, cap=POINT_CAP,
                            rng=rng)
        r = counts.loc[pathogen]
        labels.append(f"{pathogen}\n{int(r.n_columns)} cols · n={int(r.n_same_pairs)}/"
                      f"{int(r.n_diff_pairs)}")

    ax.set_ylim(-0.8, n - 0.2)
    ax.set_yticks(range(n))
    ax.set_yticklabels(labels, fontsize=5)
    ax.invert_yaxis()
    ax.set_xlim(left=0)
    ax.set_xlabel(sentence_case(f"top-{cutoff} Jaccard overlap  —  {matrix_label}"))
    ax.legend(handles=[Patch(facecolor=c_same, label="same pathogen"),
                       Patch(facecolor=c_diff, label="different pathogen")],
              loc="lower right", fontsize=6, frameon=True, facecolor="white", framealpha=0.85)
    fig.tight_layout()

    png_dir = os.path.join(output_dir, "png")
    pdf_dir = os.path.join(output_dir, "pdf")
    os.makedirs(png_dir, exist_ok=True)
    os.makedirs(pdf_dir, exist_ok=True)
    png_path = os.path.join(png_dir, name + ".png")
    pdf_path = os.path.join(pdf_dir, name + ".pdf")
    fig.savefig(png_path, dpi=300)
    fig.savefig(pdf_path)
    plt.close(fig)
    return png_path, pdf_path


# --------------------------------------------------------------------------- #
# Step 07 — mean percentile rank distribution (merged in from the former step 09, 2026-08-06)                                  #
# --------------------------------------------------------------------------- #
#: Histogram bins over the mean-rank range. 1.35M values into 120 bins leaves ~11k per bin on
#: average — fine enough to show the tails without turning the outline into noise.
N_BINS = 120


class MeanRankDistributionPlot(BasePlot):
    """Histogram of per-compound mean percentile rank across every selected endpoint.

    Each compound contributes one value: the mean of its per-column percentile ranks (a column's
    percentile rank is that compound's standing *within that endpoint* over the whole library). So
    0.5 is "average across the board", high values are compounds ranking near the top of many
    endpoints at once, low values compounds ranking near the bottom.

    A dashed reference line marks 0.5 — the value the distribution is pinned to **by construction**
    (each percentile-rank column has mean ≈ 0.5, so the grand mean of the row means is ≈ 0.5 too).
    What is informative is the SPREAD and SHAPE around it, not its centre.

    The panel carries no title, per the figure conventions; contextual counts go into the axis label.
    """

    def __init__(self, values, *, name="07_mean_rank_distribution", label_note="", ax=None,
                 cells=(3, 4)):
        super().__init__(ax=ax, cells=cells)
        self.name = name
        vals = np.asarray(values, dtype=float)
        vals = vals[np.isfinite(vals)]
        if not len(vals):
            self._unavailable()
            return
        self.ax.hist(vals, bins=N_BINS, color=hue("cobalt"), edgecolor="none")
        self.ref_line(0.5, axis="x")
        self.ax.set_xlim(0, 1)
        xlabel = "mean percentile rank across endpoints"
        if label_note:
            xlabel = f"{xlabel}\n{label_note}"
        self.label(xlabel=xlabel, ylabel="compounds")
        self.ax.text(0.02, 0.97,
                     f"n = {len(vals):,}\nmean {vals.mean():.3f}   SD {vals.std(ddof=0):.3f}\n"
                     f"min {vals.min():.3f}   max {vals.max():.3f}",
                     transform=self.ax.transAxes, va="top", ha="left",
                     fontsize=st.FONTSIZE_SMALL, color=INK)
