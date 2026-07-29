"""Figures for the CoAdd-model (eos3dys) EU OpenScreen validation (``eos3dys_validation`` subfolder
of step 05).

Panels (all dedup-preferred — CoAdd training compounds removed):
- :class:`Eos3dysHeatmapPlot` — endpoint x EU OpenScreen assay AUROC matrix, diagonal (endpoint's
  own organism) highlighted. From ``eos3dys_euos_auroc.csv``.
- :class:`Eos3dysSpecificitySwarmPlot` — same-organism vs different-organism AUROC distributions.
  From ``eos3dys_euos_auroc.csv``.
- :class:`Eos3dysHitExclusivityPlot` — AUROC on exclusive vs shared hits, per organism x endpoint.
  From ``eos3dys_hit_exclusivity.csv``; the twin of the ChEMBL ``hit_exclusivity_auroc`` panel.
- ``eos3dys_overlap_{metric}`` — the shared :class:`plots_euopenscreen.EuosOverlapTwinPlot`, reused
  once per endpoint metric on ``eos3dys_overlap_report.csv``: how much of the EU OpenScreen library
  and of its actives eos3dys already saw in CoAdd training.
- ``eos3dys_exclusive_rank_{raw,percentile}`` — the shared
  :class:`plots_euopenscreen.ExclusiveHitModelRankPlot` on ``eos3dys_exclusive_rank.csv``: for each
  exclusive hit, where its own organism's combined (inhib_50 + mic_25) score ranks among the 6
  matched organisms, stacked by pathogen.
- ``eos3dys_roc_grid_{metric}`` and ``..._{metric}_exclusive`` — the shared
  :class:`plots_euopenscreen.EuosRocGridPlot` on ``eos3dys_roc.csv``, one grid per endpoint metric
  over all actives and over exclusive actives only.
- ``eos3dys_consensus_max_by_activity`` — the shared
  :class:`plots_euopenscreen.ConsensusMaxByActivityPlot` on ``eos3dys_consensus_max_boxstats.csv``:
  max endpoint probability across the 12 matched endpoints, active vs inactive, on the RAW scale
  (all endpoints top out at ~1.0, so no percentile normalisation is needed here).

Every panel reuses a shared class rather than reimplementing it; the eos3dys-specific names are set
on the instance in :func:`save_eos3dys_figures`.

Reuses the shared primitives (``heatmap``/``diverging_cmap``/``box_with_jitter``) so the look
matches every other step; colours only from :mod:`plotting_colors`.
"""

import json
import os

import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm

from plotting_base import BasePlot
from plotting_colors import hue
from plotting_utils import abbrev, box_with_jitter, diverging_cmap, grouped_hbar, heatmap
from plots_euopenscreen import (
    ConsensusMaxByActivityPlot,
    EuosOverlapTwinPlot,
    EuosRocGridPlot,
    ExclusiveHitModelRankPlot,
)
from default import RANDOM_SEED, SHARED_ORGANISMS


def _read(subdir, fname):
    path = os.path.join(subdir, fname)
    return pd.read_csv(path) if os.path.exists(path) else pd.DataFrame()


def _prefer_set(df):
    """Keep the dedup row per (endpoint, assay) when present, else raw.

    Preferring per pair (not globally) matters: the non-organism endpoints
    (cytotoxicity/hemolytic) have no CoAdd training file, so they are raw-only and a global
    "dedup if any exists" filter would silently drop them from the matrix.
    """
    d = df.copy()
    d["_r"] = (d["set"] == "dedup").astype(int)
    d = d.sort_values("_r").drop_duplicates(["endpoint", "assay_code"], keep="last")
    return d.drop(columns="_r")


class Eos3dysHeatmapPlot(BasePlot):
    """endpoint x EU OpenScreen assay AUROC matrix for the CoAdd model (dedup). Rows are the
    model's endpoints, columns the 7 shared organisms' primary assays; the diagonal cell of each
    organism-specific endpoint (endpoint organism == assay organism) is outlined."""

    def __init__(self, df, ax=None, cells=(6, 4)):
        BasePlot.__init__(self, ax=ax, cells=cells)
        self.name = "eos3dys_euos_heatmap"
        if df is None or df.empty:
            self._unavailable()
            return
        d = _prefer_set(df)
        if d.empty:
            self._unavailable()
            return
        assay_name = dict(zip(d["assay_code"], d["assay_pathogen"]))
        ep_org = dict(zip(d["endpoint"], d["endpoint_organism"]))
        col_codes = [c for c in SHARED_ORGANISMS if c in set(d["assay_code"])]
        endpoints = sorted(d["endpoint"].unique())  # alphabetical → groups by organism prefix
        mat = d.pivot_table(index="endpoint", columns="assay_code", values="auroc",
                            aggfunc="first").reindex(index=endpoints, columns=col_codes)
        mat.columns = [abbrev(assay_name.get(c, c)) for c in col_codes]
        highlight = [(ri, ci) for ri, ep in enumerate(endpoints)
                     for ci, cc in enumerate(col_codes) if ep_org.get(ep) == cc]
        heatmap(self.ax, mat, cmap=diverging_cmap(),
                norm=TwoSlopeNorm(0.5, vmin=0.0, vmax=1.0),
                text_light_when=lambda v: v > 0.75 or v < 0.25,
                highlight=highlight, x_rotation=45, annot_fontsize=4)
        self.ax.tick_params(labelsize=5)
        self.label(title="")


class Eos3dysSpecificitySwarmPlot(BasePlot):
    """Same-organism vs different-organism AUROC for the CoAdd model's endpoints (dedup). Only
    endpoints whose organism has an EU OpenScreen assay are included (so "same" is defined); the
    two non-organism endpoints appear in the heatmap but not here. A drop toward 0.5 for the
    different-organism group would indicate organism specificity."""

    def __init__(self, df, ax=None, cells=(3, 3)):
        BasePlot.__init__(self, ax=ax, cells=cells)
        self.name = "eos3dys_same_vs_diff"
        if df is None or df.empty:
            self._unavailable()
            return
        d = _prefer_set(df)
        d = d[d["endpoint_organism"].isin(SHARED_ORGANISMS)]
        if d.empty:
            self._unavailable()
            return
        same = d[d["same_organism"]]["auroc"].dropna().to_numpy(float)
        diff = d[~d["same_organism"]]["auroc"].dropna().to_numpy(float)
        rng = np.random.default_rng(RANDOM_SEED)
        box_with_jitter(self.ax, same, 0, hue("turquoise"), vert=True, width=0.5,
                        jitter_width=0.15, point_size=12, point_alpha=0.6, rng=rng)
        box_with_jitter(self.ax, diff, 1, hue("crimson"), vert=True, width=0.5,
                        jitter_width=0.15, point_size=12, point_alpha=0.5, rng=rng)
        self.ref_line(0.5, axis="y")
        self.ax.set_xticks([0, 1])
        self.ax.set_xticklabels(["same\norganism", "different\norganism"])
        self.ax.set_xlim(-0.5, 1.5)
        self.ax.set_ylim(0, 1)
        self.label(ylabel="AUROC", title="")


class Eos3dysHitExclusivityPlot(BasePlot):
    """eos3dys AUROC on EXCLUSIVE vs SHARED EU OpenScreen hits, per organism and endpoint (dedup).

    The eos3dys twin of the ChEMBL ``hit_exclusivity_auroc`` panel: for each of the organisms that
    eos3dys and EU OpenScreen both cover, four bars — the single-point (``inhib_50``) and MIC
    (``mic_25``) endpoints, each on the organism's exclusive actives (a hit in only 1 assay) and its
    shared actives (a hit in >= 2), all against the same primary inactives. Subset is encoded by hue
    (matching that panel's amber/turquoise), endpoint by saturation. Bars carry their active count,
    which is small for the exclusive subsets — read them with the n= labels."""

    def __init__(self, excl_df, ax=None, cells=(4, 3)):
        BasePlot.__init__(self, ax=ax, cells=cells)
        self.name = "eos3dys_hit_exclusivity"
        need = {"pathogen", "metric", "subset", "auroc"}
        if excl_df is None or excl_df.empty or not need.issubset(excl_df.columns):
            self._unavailable()
            return
        d = excl_df.copy()
        d["_r"] = (d["set"] == "dedup").astype(int)          # prefer dedup per bar
        d = d.sort_values("_r").drop_duplicates(
            ["pathogen", "metric", "subset"], keep="last").drop(columns="_r")
        d["group"] = d["metric"] + " " + d["subset"]
        order = [f"{m} {s}" for m in ("inhib_50", "mic_25")
                 for s in ("exclusive", "nonexclusive")]
        order = [g for g in order if g in set(d["group"])]
        if not order:
            self._unavailable()
            return
        colors = {}
        for metric, light in (("inhib_50", None), ("mic_25", 0.55)):
            for subset, name in (("exclusive", "amber"), ("nonexclusive", "turquoise")):
                colors[f"{metric} {subset}"] = hue(name, lighten=light)
        y_labels = d.groupby("pathogen")["auroc"].max().sort_values(ascending=True).index.tolist()
        series = {(r["pathogen"], r["group"]): r["auroc"] for _, r in d.iterrows()}
        counts = {(r["pathogen"], r["group"]): r["n_active"] for _, r in d.iterrows()} \
            if "n_active" in d.columns else None
        grouped_hbar(self.ax, y_labels, series, colors, order, xlabel="AUROC", title="",
                     ref=0.5, xlim=(0, 1), counts=counts,
                     label_fn=lambda g: g.replace("nonexclusive", "shared"))


def _dedup_first(df, col="set"):
    """The dedup rows when present, else raw — the eos3dys panels all report leakage-filtered."""
    if df.empty or col not in df.columns:
        return df
    return df[df[col] == "dedup"] if "dedup" in set(df[col]) else df[df[col] == "raw"]


def save_eos3dys_figures(subdir):
    """Build the eos3dys-on-EU-OpenScreen panels from the CSVs in ``subdir`` and record their
    footprints in ``figure_cells.json``."""
    df = _read(subdir, "eos3dys_euos_auroc.csv")
    overlap = _read(subdir, "eos3dys_overlap_report.csv")
    excl = _read(subdir, "eos3dys_hit_exclusivity.csv")
    rank = _read(subdir, "eos3dys_exclusive_rank.csv")
    rank_compounds = _read(subdir, "eos3dys_exclusive_rank_compounds.csv")
    roc = _read(subdir, "eos3dys_roc.csv")
    max_stats = _read(subdir, "eos3dys_consensus_max_boxstats.csv")
    max_actives = _read(subdir, "eos3dys_consensus_max_actives.csv")
    rank_names = dict(zip(rank_compounds["code"], rank_compounds["pathogen"])) \
        if not rank_compounds.empty else None

    plots = [
        Eos3dysHeatmapPlot(df, cells=(6, 4)),
        Eos3dysSpecificitySwarmPlot(df, cells=(3, 3)),
        Eos3dysHitExclusivityPlot(excl, cells=(4, 3)),
    ]
    # one training-overlap panel per endpoint metric (the twin-axis panel takes one bar pair per row)
    if not overlap.empty and "metric" in overlap.columns:
        for metric in sorted(overlap["metric"].unique()):
            p = EuosOverlapTwinPlot(overlap[overlap["metric"] == metric], cells=(3, 4))
            p.name = f"eos3dys_overlap_{metric}"
            plots.append(p)
    # own-organism rank for exclusive hits (dedup), both rankings, stacked by pathogen
    if not rank.empty:
        for ranking in ("raw", "percentile"):
            p = ExclusiveHitModelRankPlot(_dedup_first(rank, "leakage"), ranking=ranking,
                                          name_by_code=rank_names, cells=(3, 3))
            p.name = f"eos3dys_exclusive_rank_{ranking}"
            plots.append(p)
    # ROC grids: one per endpoint metric, for all actives and for exclusive actives only
    if not roc.empty and {"metric", "subset"}.issubset(roc.columns):
        for metric in sorted(roc["metric"].unique()):
            for subset in ("all", "exclusive"):
                sub = roc[(roc["metric"] == metric) & (roc["subset"] == subset)]
                if sub.empty:
                    continue
                name = f"eos3dys_roc_grid_{metric}" + ("_exclusive" if subset == "exclusive" else "")
                p = EuosRocGridPlot(sub, set_name="dedup", cols=3, name=name)
                plots.append(p)
    # maximum endpoint probability, active vs inactive (raw scale, no normalisation needed)
    if not max_stats.empty:
        p = ConsensusMaxByActivityPlot(_dedup_first(max_stats), _dedup_first(max_actives),
                                       cells=(3, 3))
        p.name = "eos3dys_consensus_max_by_activity"
        plots.append(p)
    footprints = {}
    for p in plots:
        if p.is_available:
            p.save(subdir)
            footprints[p.name] = list(p.cells)
            print(f"  figure: eos3dys_validation/{p.name}")
        else:
            print(f"  [skip figure] {p.name}: no data")
    with open(os.path.join(subdir, "figure_cells.json"), "w") as f:
        json.dump(footprints, f, indent=2)
