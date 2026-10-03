"""xx_eosquality figure — per-pathogen eosquality scores, top-1000 hits vs the random null control.

Reads ONLY the merged summary CSV written by ``scripts/xx_eosquality.py``
(``xx_all_pathogens_summary.csv``: one row per pathogen x query set x score, with min / q25 /
median / q75 / max), per the repo's "feed figures from summary CSVs" rule — never the
per-compound score files.

One panel, four columns (typicality, extremity, support, consistency) sharing one pathogen axis.
Each row carries two horizontal boxes: the pathogen's top-1000 (crimson) and the in-sample
random-10,000 control (silver). Boxes are drawn from the precomputed quartiles; **whiskers span
min-max** (the summary carries no per-compound values to compute 1.5 x IQR whiskers from). Every
calibrated score is a CDF lookup against the reference's own distribution, so the control centres
on 0.5 by construction — marked by the dashed reference line.
"""

import os

import stylia as st

from plotting_base import MultiPanelPlot
from plotting_colors import REFERENCE_LINE, hue
from plotting_utils import (
    abbrev_ticks,
    box_from_stats,
    merge_figure_cells,
    ref_line,
    swatch_legend,
)

#: Column order: the two output-space, per-column scores first, then the two neighbourhood scores.
SCORE_ORDER = ["typicality", "extremity", "support", "consistency"]

#: Query sets drawn per row, as ``(query_set, legend label, colour, vertical offset)``.
QUERY_SETS = [
    ("top1000", "Top-1000", hue("crimson"), -0.18),
    ("random10000", "Random control", REFERENCE_LINE, 0.18),
]

BOX_WIDTH = 0.3


class EosqualityPathogenPlot(MultiPanelPlot):
    """Four score columns x 15 pathogen rows of top-1000 vs control boxes.

    ``summary`` is the merged summary table; ``order`` is a DataFrame with ``code`` and
    ``pathogen`` (full organism name) columns giving the row order, top to bottom.
    """

    def __init__(self, summary, order, cells=(3, 6)):
        axs = self._new_figure(1, len(SCORE_ORDER), cells, "xx_eosquality_pathogen_quality")
        codes = list(order["code"])
        names = dict(zip(order["code"], order["pathogen"]))
        data = summary[summary["pathogen"].isin(codes) & summary["score"].isin(SCORE_ORDER)]
        if data.empty:
            self._unavailable()
            return
        for i, score in enumerate(SCORE_ORDER):
            ax = axs.next()
            self._panel(ax, data[data["score"] == score], codes, names, first=(i == 0))
            st.label(ax, xlabel=score, ylabel="")
        # Above the last column: every in-axes region holds whiskers, so any inside placement hides data.
        swatch_legend(ax, {label: color for _, label, color, _ in QUERY_SETS}, loc="lower right",
                      bbox_to_anchor=(1.0, 1.0), ncol=len(QUERY_SETS))

    def _panel(self, ax, sub, codes, names, first):
        for yi, code in enumerate(codes):
            for query_set, _, color, dy in QUERY_SETS:
                row = sub[(sub["pathogen"] == code) & (sub["query_set"] == query_set)]
                if row.empty:
                    continue
                r = row.iloc[0]
                stats = {"median": r["median"], "q1": r["q25"], "q3": r["q75"],
                         "whisker_lo": r["min"], "whisker_hi": r["max"]}
                box_from_stats(ax, stats, yi + dy, color, vert=False, width=BOX_WIDTH)
        ref_line(ax, 0.5, axis="x")
        ax.set_xlim(0, 1)
        ax.set_ylim(len(codes) - 0.5, -0.5)
        if first:
            abbrev_ticks(ax, codes, axis="y", name_map=names)
        else:
            ax.set_yticks(range(len(codes)))
            ax.set_yticklabels([])


def save_eosquality_figure(summary, order, output_dir):
    """Build, save (PNG + PDF) and record the footprint of the per-pathogen panel."""
    plot = EosqualityPathogenPlot(summary, order)
    footprints = {}
    if plot.is_available:
        plot.save(output_dir)
        footprints[plot.name] = list(plot.cells)
        print(f"  figure: {plot.name}  ({plot.cells[0]} x {plot.cells[1]} cells)"
              f" -> {os.path.join(output_dir, 'pdf')}")
    else:
        print("  [skip figure] no eosquality summary rows to plot")
    merge_figure_cells(output_dir, footprints)
    return plot
