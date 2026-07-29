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

from plotting_base import BasePlot, GridPlot
from plotting_colors import SHARED_ORGANISM_COLORS, hue
from plotting_utils import (
    LEGEND_KW,
    abbrev,
    box_from_stats,
    diverging_cmap,
    grouped_hbar,
    heatmap,
    marker_legend,
    roc_panel,
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

# exclusive vs non-exclusive(shared) hits.
EXCLUSIVITY_COLORS = {"exclusive": hue("amber"), "nonexclusive": hue("turquoise")}
# hit classes for the score-distribution box: the inactive background, organism-specific hits, and
# the shared hits split into narrow- (2-3 pathogens) and broad-spectrum (>3).
HIT_CLASS_COLORS = {
    "inactive": hue("cobalt"),
    "exclusive": hue("amber"),
    "narrow": hue("turquoise"),
    "broad": hue("crimson"),
    "active": hue("crimson"),   # the four active classes collapsed to one (activity view)
}
# enrichment factor top-fractions.
EF_COLORS = {"ef_1pct": hue("cobalt"), "ef_5pct": hue("turquoise")}
# primary screen vs merged secondary (confirmatory) assays.
SCREEN_COLORS = {"primary": hue("cobalt"), "secondary": hue("tangerine")}

_PRETTY = {
    "raw": "raw", "dedup": "dedup (no leakage)",
    "exclusive": "exclusive", "nonexclusive": "shared",
    "inactive": "inactive", "narrow": "narrow\n(2-3)", "broad": "broad\n(>3)",
    "active": "active\n(>=1)",
    "inhib_50": "single-point (inhib_50)", "mic_10": "MIC (mic_10)",
    "ef_1pct": "EF @ 1%", "ef_5pct": "EF @ 5%",
}


def _pretty(key):
    return _PRETTY.get(key, key)


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
                 cells=(3, 3), ax=None):
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
        # order organisms by their best available metric value (ascending → best on top)
        order_key = d.groupby("pathogen")[metric].max().sort_values(ascending=True)
        y_labels = order_key.index.tolist()
        series = {(r["pathogen"], r[group_col]): r[metric] for _, r in d.iterrows()}
        counts = None
        if count_col is not None and count_col in d.columns:
            counts = {(r["pathogen"], r[group_col]): r[count_col] for _, r in d.iterrows()}
        grouped_hbar(self.ax, y_labels, series, group_colors, group_order,
                     xlabel=xlabel, title=title, ref=ref, xlim=xlim, counts=counts,
                     label_fn=_pretty)


class EuosSharedEnrichmentPlot(BasePlot):
    """Analysis 1 (companion) — EU OpenScreen enrichment factor (top 1% / 5%) per organism,
    deduplicated. Enrichment = hit rate in the top k% over the base rate (1 = no enrichment)."""

    def __init__(self, own_df, ax=None, cells=(3, 3)):
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
    stacked novel (turquoise) vs in-training (crimson)."""

    def __init__(self, leak_df, ax=None, cells=(3, 4)):
        BasePlot.__init__(self, ax=ax, cells=cells)
        self.name = "euos_overlap"
        need = {"n_active", "n_eval_conclusive", "n_overlap", "n_overlap_active"}
        if leak_df is None or leak_df.empty or not need.issubset(leak_df.columns):
            self._unavailable()
            return
        d = leak_df[leak_df["n_active"] > 0].sort_values(
            "n_active", ascending=True).reset_index(drop=True)
        if d.empty:
            self._unavailable()
            return
        novel, in_training = hue("turquoise"), hue("crimson")
        y = np.arange(len(d))
        h = 0.38
        lib_novel = (d["n_eval_conclusive"] - d["n_overlap"]).values
        lib_ovl = d["n_overlap"].values
        act_novel = (d["n_active"] - d["n_overlap_active"]).values
        act_ovl = d["n_overlap_active"].values

        # top axis = library compounds; bottom axis (self.ax) = actives (twiny shares the y-axis)
        axt = self.ax.twiny()
        axt.barh(y + 0.21, lib_novel, height=h, color=novel, edgecolor="white", hatch="//")
        axt.barh(y + 0.21, lib_ovl, left=lib_novel, height=h, color=in_training,
                 edgecolor="white", hatch="//")
        self.ax.barh(y - 0.21, act_novel, height=h, color=novel)
        self.ax.barh(y - 0.21, act_ovl, left=act_novel, height=h, color=in_training)

        self.ax.set_ylim(-0.6, len(d) - 0.4)
        self.ax.set_yticks(y)
        self.ax.set_yticklabels([abbrev(p) for p in d["pathogen"]])
        self.ax.set_xlim(0, float((act_novel + act_ovl).max()) * 1.05)
        axt.set_xlim(0, float(d["n_eval_conclusive"].max()) * 1.05)

        # one combined legend outside the plot: colour = novel/in-training, texture = which axis
        handles = [
            Patch(facecolor=novel, label="novel to model"),
            Patch(facecolor=in_training, label="in training set"),
            Patch(facecolor=novel, edgecolor="white", hatch="//", label="library"),
            Patch(facecolor=novel, label="actives"),
        ]
        self.ax.legend(handles=handles, fontsize=5, loc="center left",
                       bbox_to_anchor=(1.01, 0.5), **LEGEND_KW)

        self.label(xlabel="actives", title="")
        axt.set_xlabel("library compounds")


class HitExclusivityPlot(MetricByOrganismPlot):
    """Analysis 3 — exclusive vs shared (non-exclusive) hit AUROC per organism (dedup)."""

    def __init__(self, excl_df, ax=None, cells=(3, 3)):
        super().__init__(
            excl_df, group_col="subset", group_colors=EXCLUSIVITY_COLORS,
            group_order=["exclusive", "nonexclusive"], metric="auroc", xlabel="AUROC",
            title="", name="hit_exclusivity_auroc", count_col="n_active",
            ref=0.5, xlim=(0, 1), prefer_best=True, cells=cells, ax=ax)


class PrimaryVsSecondaryAurocPlot(MetricByOrganismPlot):
    """Analysis 1b — own-assay AUROC on the primary vs the merged secondary EU OpenScreen
    screen, per organism (dedup preferred). Whether the model's ranking holds on the smaller
    confirmatory / dose-response assays."""

    def __init__(self, combined_df, ax=None, cells=(3, 3)):
        super().__init__(
            combined_df, group_col="screen", group_colors=SCREEN_COLORS,
            group_order=["primary", "secondary"], metric="auroc", xlabel="AUROC",
            title="", name="primary_vs_secondary_auroc",
            ref=0.5, xlim=(0, 1), prefer_best=True, cells=cells, ax=ax)


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


class HitPromiscuityPlot(BasePlot):
    """Hit promiscuity (label-only) — how many EU OpenScreen actives are hits in 1, 2, ... 7 of
    the primary assays. The per-compound counterpart of the Jaccard overlap heatmap: most actives
    are organism-specific singletons, while a small tail is pan-active across the panel.

    The y axis is log-scaled because the distribution spans two orders of magnitude (hundreds of
    singletons vs a handful of 7-pathogen hits); every bar is annotated with its exact count.
    """

    def __init__(self, prom_df, ax=None, cells=(3, 3)):
        BasePlot.__init__(self, ax=ax, cells=cells)
        self.name = "hit_promiscuity"
        if prom_df is None or prom_df.empty or "n_molecules" not in prom_df.columns:
            self._unavailable()
            return
        d = prom_df.sort_values("n_pathogens")
        x = d["n_pathogens"].to_numpy()
        y = d["n_molecules"].to_numpy(dtype=float)
        self.ax.bar(x, y, color=hue("turquoise"))
        self.ax.set_yscale("log")
        self.ax.set_ylim(0.7, float(y.max()) * 3.0 if y.max() > 0 else 10.0)
        self.ax.set_xticks(x)
        for xi, yi in zip(x, y):
            if yi <= 0:
                continue
            self.ax.annotate(f"{int(yi)}", (xi, yi), ha="center", va="bottom",
                             fontsize=5, xytext=(0, 1.5), textcoords="offset points")
        total = int(d["n_molecules"].sum())
        self.label(xlabel="number of pathogens the compound is a hit for",
                   ylabel=f"EU OpenScreen actives (n = {total})", title="")


class ExclusiveHitModelRankPlot(BasePlot):
    """For each exclusive hit (active in exactly 1 of the 7 primary assays), where its own
    pathogen's model ranks it among all 7 models: rank 1 = its own pathogen scores it highest,
    rank 2 = one other pathogen's model ranks it higher, and so on to rank 7.

    One panel per ranking mode — ``raw`` consensus scores as-is, or ``percentile`` (each score first
    converted to its percentile within that model's own library distribution, which puts the models
    on a common scale; the raw scores are not calibrated across models). The dashed line is the
    chance level, n_hits / 7.

    Bar height is the molecule count, split into segments coloured by the pathogen each hit belongs
    to (from the ``n_<code>`` columns of the summary CSV). Segment sizes therefore track how many
    exclusive hits each pathogen has in the first place — the legend gives those totals so a large
    segment isn't misread as a per-pathogen effect."""

    def __init__(self, rank_df, ranking="raw", name_by_code=None, dedup=False,
                 ax=None, cells=(3, 3)):
        BasePlot.__init__(self, ax=ax, cells=cells)
        self.name = f"exclusive_hit_model_rank_{ranking}" + ("_dedup" if dedup else "")
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
                self.ax.bar(x, seg, bottom=bottom, color=SHARED_ORGANISM_COLORS[code])
                bottom += seg
            self.legend({f"{abbrev(c, name_by_code)} ({int(d[f'n_{c}'].sum())})":
                         SHARED_ORGANISM_COLORS[c] for c in codes},
                        loc="upper right", ncol=2)
        else:
            self.ax.bar(x, y, color=hue("cobalt", lighten=0.55))
        if "n_chance" in d.columns:
            self.ref_line(float(d["n_chance"].iloc[0]), axis="y")
        self.ax.set_xticks(x)
        for xi, yi in zip(x, y):
            if yi > 0:
                self.ax.annotate(f"{int(yi)}", (xi, yi), ha="center", va="bottom",
                                 fontsize=5, xytext=(0, 1.5), textcoords="offset points")
        total = int(d["n_total"].iloc[0]) if "n_total" in d.columns else int(y.sum())
        suffix = ", no training-set hits" if dedup else ""
        n_ranked = int(d["n_models_ranked"].iloc[0]) if "n_models_ranked" in d.columns else len(x)
        self.label(xlabel=f"rank of the hit pathogen's own model (of {n_ranked})",
                   ylabel=f"exclusive hits (n = {total}{suffix})", title="")


class ScoreByHitClassPlot(BasePlot):
    """Boxes of a per-compound aggregated consensus score, one box per EU OpenScreen hit class.

    A reusable panel type (both score-distribution figures below build on it). Boxes are drawn from
    the PRECOMPUTED statistics — the inactive class has ~10^5 compounds and is never shipped
    per-molecule — while the active classes, whose individual values the summary CSV does carry,
    also get a jittered point overlay. Raw scores: training-set compounds are deliberately
    included, so these panels describe the score distribution, not out-of-sample performance.
    """

    def __init__(self, stats_df, actives_df, classes, score_col, ylabel, name,
                 ax=None, cells=(3, 3)):
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
            box_from_stats(self.ax, by_class.loc[hit_class], i, HIT_CLASS_COLORS[hit_class],
                           face=HIT_CLASS_COLORS[hit_class], points=points, rng=rng)
        n_models = int(stats_df["n_models_aggregated"].iloc[0]) \
            if "n_models_aggregated" in stats_df.columns else len(SHARED_ORGANISMS)
        self.ax.set_xticks(range(len(order)))
        self.ax.set_xticklabels([f"{_pretty(c)}\n(n = {int(by_class.loc[c, 'n']):,})"
                                 for c in order])
        self.label(xlabel="", ylabel=ylabel.format(n_models=n_models), title="")


class ConsensusSumByHitClassPlot(ScoreByHitClassPlot):
    """Summed consensus score (over the 7 shared-organism models) by hit class: inactive everywhere
    / hit in exactly 1 pathogen (exclusive) / narrow-spectrum (2-3) / broad-spectrum (>3)."""

    def __init__(self, stats_df, actives_df=None, ax=None, cells=(3, 3)):
        super().__init__(stats_df, actives_df, classes=HIT_CLASSES,
                         score_col="consensus_sum",
                         ylabel="summed consensus score ({n_models} models)",
                         name="consensus_sum_by_hit_class", ax=ax, cells=cells)


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
        ylabel = "library percentile of best-model score ({n_models} models)" if normalized \
            else "maximum consensus score ({n_models} models)"
        name = "consensus_max_percentile_by_activity" if normalized \
            else "consensus_max_by_activity"
        if dedup:
            name += "_dedup"
            ylabel += ", no training-set compounds"
        super().__init__(stats_df, actives_df, classes=ACTIVITY_CLASSES,
                         score_col="consensus_max", ylabel=ylabel, name=name,
                         ax=ax, cells=cells)


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

    def __init__(self, spec_df, ax=None, cells=(3, 3)):
        BasePlot.__init__(self, ax=ax, cells=cells)
        self.name = "specificity_index"
        if spec_df is None or spec_df.empty or spec_df["specificity_index"].dropna().empty:
            self._unavailable()
            return
        specificity_bars(self.ax, spec_df, title="")


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
    dot per sub-model output's own-assay AUROC (dedup) and the ``consensus_score`` drawn as a
    larger bubble — shows at a glance where the consensus sits within the spread of its ensemble
    members across all organisms. Mirrors the old ``xx_euopenscreen_preds`` AUROC dot plot."""

    def __init__(self, auroc_df, ax=None, cells=(3, 3)):
        BasePlot.__init__(self, ax=ax, cells=cells)
        self.name = "submodel_auroc_summary"
        if auroc_df is None or auroc_df.empty:
            self._unavailable()
            return
        d = auroc_df.copy()
        d["_r"] = (d["set"] == "dedup").astype(int)  # prefer dedup per (pathogen, feature)
        d = d.sort_values("_r").drop_duplicates(["pathogen", "feature"], keep="last")
        # order pathogens by consensus AUROC (best on top); fall back to mean if no consensus
        cons = d[d["feature"] == "consensus_score"].set_index("pathogen")["auroc"]
        if cons.empty:
            order = d.groupby("pathogen")["auroc"].mean().sort_values(ascending=True).index.tolist()
        else:
            order = cons.sort_values(ascending=True).index.tolist()
        idx = {p: i for i, p in enumerate(order)}
        sub = d[d["feature"] != "consensus_score"]
        con = d[d["feature"] == "consensus_score"]
        self.ax.scatter(sub["auroc"], [idx[p] for p in sub["pathogen"]], color=hue("cobalt"),
                        s=16, alpha=0.6, edgecolors="none", zorder=2)
        self.ax.scatter(con["auroc"], [idx[p] for p in con["pathogen"]], color=hue("crimson"),
                        s=90, alpha=0.85, edgecolor="white", linewidth=0.5, zorder=3)
        self.ref_line(0.5, axis="x")
        self.ax.set_yticks(range(len(order)))
        self.ax.set_yticklabels([abbrev(p) for p in order])
        self.ax.set_xlim(0.2, 1)
        marker_legend(self.ax, [
            {"label": "sub-model", "color": hue("cobalt"), "markersize": 5},
            {"label": "consensus", "color": hue("crimson"), "markersize": 8},
        ], loc="lower left")
        self.label(xlabel="AUROC on own EU OpenScreen assay (dedup)", title="")


# --------------------------------------------------------------------------- #
# Entry point                                                                  #
# --------------------------------------------------------------------------- #
def _read(output_dir, fname):
    path = os.path.join(output_dir, fname)
    return pd.read_csv(path) if os.path.exists(path) else pd.DataFrame()


def _combine_screens(primary, secondary):
    """Stack the primary + secondary own-assay records into one frame with a ``screen`` column,
    for the primary-vs-secondary comparison panel."""
    frames = []
    if not primary.empty:
        p = primary.copy(); p["screen"] = "primary"; frames.append(p)
    if not secondary.empty:
        s = secondary.copy(); s["screen"] = "secondary"; frames.append(s)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


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
    sec = _read(dedup_dir, "05_euopenscreen_secondary_auroc.csv")
    roc = _read(dedup_dir, "05_euopenscreen_roc.csv")
    excl = _read(dedup_dir, "05_hit_exclusivity.csv")
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

    groups = [
        (output_dir, "", [
            EuosOverlapTwinPlot(leak, cells=(3, 4)),
            ActiveOverlapHeatmapPlot(overlap, cells=(4, 4)),
            HitPromiscuityPlot(promiscuity, cells=(3, 3)),
        ]),
        (full_dir, FULL_SUBDIR, [
            ConsensusSumByHitClassPlot(sum_stats, sum_actives, cells=(3, 3)),
            ConsensusMaxByActivityPlot(max_stats, max_actives, cells=(3, 3)),
            ConsensusMaxByActivityPlot(maxp_stats, maxp_actives, normalized=True, cells=(3, 3)),
            ExclusiveHitModelRankPlot(model_rank, ranking="raw", name_by_code=rank_names,
                                      cells=(3, 3)),
            ExclusiveHitModelRankPlot(model_rank, ranking="percentile", name_by_code=rank_names,
                                      cells=(3, 3)),
        ]),
        (dedup_dir, DEDUP_SUBDIR, [
            EuosRocGridPlot(roc, set_name="dedup", cols=3),
            EuosSharedEnrichmentPlot(own, cells=(3, 3)),
            PrimaryVsSecondaryAurocPlot(_combine_screens(own, sec), cells=(3, 3)),
            HitExclusivityPlot(excl, cells=(3, 3)),
            CrossOrganismHeatmapPlot(cross, cells=(4, 4)),
            SpecificityIndexPlot(spec, cells=(3, 3)),
            SubmodelAurocSummaryPlot(submodel, cells=(3, 3)),
            ConsensusMaxByActivityPlot(maxd_stats, maxd_actives, normalized=True, dedup=True,
                                       cells=(3, 3)),
            ExclusiveHitModelRankPlot(model_rank_dedup, ranking="percentile",
                                      name_by_code=rank_names, dedup=True, cells=(3, 3)),
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
