"""Step 07 figures — inter-model prediction correlation panels.

Reads ONLY the small summary CSVs written by :func:`eval_correlations.run_analyze`
(``07_spearman_corr.csv``, ``07_topn_overlap_N*.csv``, ``07_column_index.csv``,
``07_group_assignments.csv``, ``07_group_correlation_summary.csv``) — never the per-compound score
matrix. Every panel is an individual :class:`BasePlot`, saved as PNG + PDF with its footprint in
``figure_cells.json``. Colours come only from :mod:`plotting_colors`.
"""

import json
import os

import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
import stylia as st

from default import TOPN_CUTOFFS
from eval_correlations import build_groups
from plotting_base import BasePlot
from plotting_colors import hue
from plotting_utils import diverging_cmap, heatmap, sequential_cmap, swatch_legend

CORR_NORM = mcolors.Normalize(vmin=-1.0, vmax=1.0)


def _clamp_cells(n, lo=3, hi=6, per=4):
    """Square cell footprint scaled to node count (kept within the page grid)."""
    return (max(lo, min(hi, round(n / per) or lo)),) * 2


def _short(node):
    """Compact node label ``model:output`` for a labelled heatmap axis."""
    model_id, _, col = node.partition(":")
    return f"{model_id}·{col}" if col else model_id


# --------------------------------------------------------------------------- #
# Global overview                                                              #
# --------------------------------------------------------------------------- #
class GlobalCorrPlot(BasePlot):
    """Full node × node Spearman heatmap, hierarchically ordered, unlabelled (block overview)."""

    def __init__(self, corr, ax=None, cells=(6, 6)):
        super().__init__(ax=ax, cells=cells)
        self.name = "07_corr_global"
        if corr.empty:
            self._unavailable()
            return
        cmap = diverging_cmap("crimson", "white", "cobalt")
        heatmap(self.ax, corr, cmap=cmap, norm=CORR_NORM, annotate=False,
                row_labels=[""] * len(corr), col_labels=[""] * len(corr))
        self.ax.set_xticks([])
        self.ax.set_yticks([])
        self.label(xlabel=f"{len(corr)} model output columns",
                   ylabel=f"{len(corr)} model output columns")


class CorrDistributionPlot(BasePlot):
    """Histogram of pairwise |Spearman|: within-focus-group vs cross-group."""

    def __init__(self, corr, col_index, group_df, ax=None, cells=(3, 4)):
        super().__init__(ax=ax, cells=cells)
        self.name = "07_corr_distribution"
        if corr.empty:
            self._unavailable()
            return
        groups = build_groups(col_index, group_df)
        member_nodes = sorted({n for nodes in groups.values() for n in nodes if n in corr.index})
        abs_corr = corr.abs()
        iu = np.triu_indices(len(corr), k=1)
        all_vals = abs_corr.to_numpy()[iu]
        within = []
        for nodes in groups.values():
            nodes = [n for n in nodes if n in abs_corr.index]
            if len(nodes) < 2:
                continue
            sub = abs_corr.loc[nodes, nodes].to_numpy()
            within.extend(sub[np.triu_indices(len(nodes), k=1)].tolist())
        bins = np.linspace(0, 1, 26)
        self.ax.hist(all_vals, bins=bins, density=True, color=hue("turquoise"),
                     alpha=0.55, label="all pairs")
        if within:
            self.ax.hist(within, bins=bins, density=True, color=hue("cobalt"),
                         alpha=0.55, label="within focus group")
        swatch_legend(self.ax, {"all pairs": hue("turquoise"),
                                "within focus group": hue("cobalt")}, loc="upper right")
        self.label(xlabel="|Spearman rho|", ylabel="density")


# --------------------------------------------------------------------------- #
# Focus-group heatmaps (labelled)                                              #
# --------------------------------------------------------------------------- #
class FocusHeatmapPlot(BasePlot):
    """Labelled square heatmap of a matrix restricted to one focus group's nodes.

    Used for both the Spearman focus panels (``kind='corr'``) and the top-N overlap panels
    (``kind='overlap'``); the two differ only in colour mapping and value scale.
    """

    def __init__(self, matrix, nodes, name, *, kind="corr", ax=None):
        nodes = [n for n in nodes if n in matrix.index and n in matrix.columns]
        cells = _clamp_cells(len(nodes))
        super().__init__(ax=ax, cells=cells)
        self.name = name
        if len(nodes) < 2:
            self._unavailable()
            return
        sub = matrix.loc[nodes, nodes]
        labels = [_short(n) for n in nodes]
        if kind == "corr":
            cmap, norm = diverging_cmap("crimson", "white", "cobalt"), CORR_NORM
            light = lambda v: abs(v) > 0.6  # noqa: E731
        else:
            vmax = float(np.nanmax(sub.to_numpy())) or 1.0
            cmap, norm = sequential_cmap("cobalt"), mcolors.Normalize(0.0, vmax)
            light = lambda v: v > 0.6 * vmax  # noqa: E731
        annotate = len(nodes) <= 20
        labelled = len(nodes) <= 30  # beyond this, ticks are unreadable — show block structure only
        heatmap(self.ax, sub, cmap=cmap, norm=norm, annotate=annotate, value_fmt="{:.2f}",
                text_light_when=light, row_labels=labels if labelled else [""] * len(nodes),
                col_labels=labels if labelled else [""] * len(nodes), colorbar=True)
        self.label(xlabel="", ylabel="")  # node names live on the ticks; clear stylia placeholders
        self.ax.set_xlabel("")
        self.ax.set_ylabel("")
        if not labelled:
            self.ax.set_xticks([])
            self.ax.set_yticks([])
            self.label(xlabel=f"{len(nodes)} output columns", ylabel=f"{len(nodes)} output columns")
        for lbl in self.ax.get_xticklabels():
            lbl.set_fontsize(st.FONTSIZE_SMALL)
        for lbl in self.ax.get_yticklabels():
            lbl.set_fontsize(st.FONTSIZE_SMALL)


# --------------------------------------------------------------------------- #
# Entry point                                                                  #
# --------------------------------------------------------------------------- #
def _read(output_dir, fname, **kw):
    path = os.path.join(output_dir, fname)
    return pd.read_csv(path, **kw) if os.path.exists(path) else pd.DataFrame()


def _org_slug(label):
    """`organism:Escherichia coli` -> `ecoli` for a filesystem-safe figure name."""
    org = label.split(":", 1)[1]
    return "".join(ch for ch in org.lower() if ch.isalnum())


def save_correlation_figures(output_dir):
    """Build every step-07 panel from the summary CSVs in ``output_dir`` and record footprints."""
    corr = _read(output_dir, "07_spearman_corr.csv", index_col=0)
    col_index = _read(output_dir, "07_column_index.csv")
    group_df = _read(output_dir, "07_group_assignments.csv")
    if group_df.empty:
        group_df = pd.DataFrame(columns=["model_id", "is_cytotox", "organisms"])
    group_df = group_df.fillna({"organisms": "", "cytotox_evidence": ""})

    plots = [
        GlobalCorrPlot(corr),
        CorrDistributionPlot(corr, col_index, group_df),
    ]

    # Focus-group Spearman panels (cytotoxicity + each qualifying organism).
    groups = build_groups(col_index, group_df) if not corr.empty else {}
    for label, nodes in groups.items():
        slug = "cytotoxicity" if label == "cytotoxicity" else _org_slug(label)
        plots.append(FocusHeatmapPlot(corr, nodes, f"07_corr_{slug}", kind="corr"))

    # Top-N overlap: global overview (unlabelled via the focus helper only for the cytotox group).
    for n in TOPN_CUTOFFS:
        jac = _read(output_dir, f"07_topn_overlap_N{n}.csv", index_col=0)
        if jac.empty:
            continue
        cyto = groups.get("cytotoxicity", [])
        plots.append(FocusHeatmapPlot(jac, cyto, f"07_overlap_N{n}_cytotoxicity", kind="overlap"))

    footprints = {}
    for p in plots:
        if p.is_available:
            p.save(output_dir)
            footprints[p.name] = list(p.cells)
            print(f"  figure: {p.name}")
        else:
            print(f"  [skip figure] {p.name}: insufficient data")
    with open(os.path.join(output_dir, "figure_cells.json"), "w") as f:
        json.dump(footprints, f, indent=2)
