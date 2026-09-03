"""Step 12 figure — the eos19mt Fisher-enrichment matrix.

One cell per (eos19mt antibiotic class, pathogen) pair, from :mod:`eval_abx_enrichment`'s wide
odds-ratio/p-value matrices. Colour is a diverging scale on log2(odds ratio), centred on 0 (odds
ratio == 1, no association) — cool where a class is under-represented in the pathogen's top-N,
warm where it is over-represented, matching the AUROC matrix's cool-low/warm-high convention.
Some pairs have an odds ratio of exactly 0 or infinity (e.g. a class with a single non-zero
compound library-wide, which lands either entirely inside or entirely outside a pathogen's top-N)
— the colour scale is clipped at :data:`default.ABX_ENRICHMENT_LOG2OR_CAP` so these cells still
render at the scale's extreme rather than breaking the colormap, but the printed cell text is
always the TRUE odds ratio (``0``/``inf``), never the clipped colour value.

A trailing ``*``/``**``/``***`` marks the raw (unadjusted) per-cell p-value at
:data:`default.ABX_ENRICHMENT_SIG_THRESHOLDS` — a display convenience, not a multiple-testing
claim (the long-format CSV carries a Benjamini-Hochberg-adjusted column for that).
"""

import matplotlib.colors as mcolors
import numpy as np
import stylia as st
from mpl_toolkits.axes_grid1 import make_axes_locatable

from default import ABX_ENRICHMENT_LOG2OR_CAP, ABX_ENRICHMENT_SIG_THRESHOLDS, ORGANISM_CLASS_ORDER
from plotting_colors import INK, REFERENCE_LINE, distinct_colors
from plotting_utils import diverging_cmap, heatmap, merge_figure_cells
from plots_auroc_matrix import _MatrixPlotBase, _organism_label

#: Near-full-page footprint: 38 rows need a lot of vertical room, and the cell text (odds ratios up
#: to 5 digits plus significance asterisks) needs the full print width's worth of column room too —
#: tuned visually against the rendered figure (an initial (7.2, 3) let neighbouring cells' text run
#: into each other), not derived from the row/column counts.
CELLS = (7.2, 6)


def _or_text(v):
    """The true odds ratio as display text — ``inf``/``0`` shown as such, never clipped."""
    if np.isnan(v):
        return ""
    if np.isinf(v):
        return "inf" if v > 0 else "0"
    if v >= 100:
        return f"{v:.0f}"
    if v >= 10:
        return f"{v:.1f}"
    return f"{v:.2f}"


def _sig_stars(p, thresholds=ABX_ENRICHMENT_SIG_THRESHOLDS):
    """``*``/``**``/``***`` for the raw p-value at each successive threshold, ``""`` if not
    significant or missing (degenerate class).
    """
    if not np.isfinite(p):
        return ""
    return "*" * sum(p < t for t in thresholds)


class AbxEnrichmentMatrixPlot(_MatrixPlotBase):
    """The 38 (eos19mt class) x 15 (pathogen) Fisher-enrichment matrix, with one top track for
    each pathogen's organism class (:data:`default.ORGANISM_CLASS_ORDER`) — reusing
    :class:`plots_auroc_matrix._MatrixPlotBase`'s track/colorbar chrome, on the same page-budget
    3 cm cell grid as every other matrix figure in this repo.
    """

    def __init__(self, odds_ratio, p_value, pathogens, ax=None, cells=CELLS, name=None,
                cap=ABX_ENRICHMENT_LOG2OR_CAP):
        super().__init__(ax=ax, cells=cells)
        self.name = name or "12_abx_enrichment_matrix"
        if not odds_ratio.size:
            self._unavailable()
            return

        with np.errstate(divide="ignore", invalid="ignore"):
            log2or = np.log2(odds_ratio.astype(float))
        clipped = log2or.clip(-cap, cap)

        cmap = diverging_cmap(low="cobalt", mid="white", high="crimson").copy()
        cmap.set_bad(REFERENCE_LINE)
        norm = mcolors.Normalize(vmin=-cap, vmax=cap)

        col_labels = [_organism_label(p) for p in pathogens["pathogen"]]
        heatmap(self.ax, clipped, cmap=cmap, norm=norm, annotate=False, nan_color=REFERENCE_LINE,
               x_rotation=90, row_labels=list(odds_ratio.index), col_labels=col_labels,
               annot_fontsize=st.FONTSIZE_SMALL, aspect="auto")

        clip_arr = clipped.to_numpy()
        or_arr = odds_ratio.to_numpy()
        p_arr = p_value.to_numpy()
        for i in range(clip_arr.shape[0]):
            for j in range(clip_arr.shape[1]):
                if np.isnan(clip_arr[i, j]):
                    continue  # degenerate class: no test, nothing to print — colour marks it NaN
                light = abs(clip_arr[i, j]) >= cap * 0.5
                text = _or_text(or_arr[i, j]) + _sig_stars(p_arr[i, j])
                self.ax.text(j, i, text, ha="center", va="center",
                            fontsize=st.FONTSIZE_SMALL, color="white" if light else INK)

        self.ax.tick_params(labelsize=st.FONTSIZE_SMALL, length=0)
        self.ax.grid(False)
        self.ax.set_axisbelow(False)
        self.ax.set_xlabel("")
        self.ax.set_ylabel("")

        divider = make_axes_locatable(self.ax)
        class_colors = dict(zip(ORGANISM_CLASS_ORDER,
                                distinct_colors(len(ORGANISM_CLASS_ORDER), levels=(0.55, 0.3))))
        self._track(divider, "top", pathogens["organism_class"].tolist(), class_colors,
                   "organism_class", vertical=False)
        self._colorbar(cmap, norm, divider, label="log2(odds ratio)", continuous=True)


def save_enrichment_figure(output_dir, odds_ratio, p_value, pathogens, name=None):
    """Build the Fisher-enrichment matrix figure and merge its footprint into
    ``figure_cells.json``.
    """
    plot = AbxEnrichmentMatrixPlot(odds_ratio, p_value, pathogens, name=name)
    footprints = {}
    if plot.is_available:
        plot.save(output_dir)
        footprints[plot.name] = list(plot.cells)
        print(f"  figure: {plot.name}")
    else:
        print(f"  [skip figure] {plot.name}: empty matrix")
    return merge_figure_cells(output_dir, footprints)
