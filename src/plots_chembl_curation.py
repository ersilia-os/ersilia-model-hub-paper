"""Concrete plot classes and the figure entry point for the ChEMBL data-curation figures
(scripts/02_chembl_data_curation.py).

Each panel reproduces a step-27 figure from ``chembl-antimicrobial-tasks``
(``scripts/27_general_plots.py``) but is rebuilt here **purely from the per-pathogen and
aggregate SUMMARY CSVs** staged into ``data/raw/chembl_curation/`` — no full molecule-level
datasets. The aggregation/binning of every panel is ported verbatim from the named upstream
function so the numbers match; the plotting is re-expressed on the paper's stylia stack
(individual figures, 3 cm cell grid, PNG + vector PDF, no A/B/C letters).

Figures whose values are NOT in the summary CSVs (chemical-space overlap/coverage, embedding
scatters, molecule-level conflict panels) are intentionally not built here.

Follows the shared plotting layer: panels subclass :class:`plotting_base.BasePlot` /
:class:`plotting_base.MultiPanelPlot` and draw through the :mod:`plotting_utils` primitives
(``stacked_hbar``, ``box_with_jitter``, ``ref_line``, ``swatch_legend``) so their look matches
every other step's figures. Colours come only from :mod:`plotting_colors`, and pathogen ticks
are abbreviated genus names (``A. baumannii``) via a code→name map.
"""

import json
import os

import numpy as np
import pandas as pd
import stylia

from plotting_base import BasePlot, MultiPanelPlot
from plotting_colors import AUROC_PASS_COLORS, CATEGORY_COLORS, INK, REFERENCE_LINE, hue, shades
from plotting_utils import (
    LEGEND_KW,
    abbrev,
    box_with_jitter,
    ref_line,
    stacked_hbar,
    swatch_legend,
)

# --------------------------------------------------------------------------------------------
# Constants ported verbatim from the upstream curation pipeline (plot_utils.py / 27_general_plots.py).
# These encode inherited (not newly chosen) binning and ordering — see scripts/README.md.
# --------------------------------------------------------------------------------------------

# Curation-outcome taxonomy shared by the discard-reason and chemical-space attrition figures so
# the two read as a pair. Derived from the ``discard_reason`` TEXT written in 21_curation_summary
# (NOT the numeric discard_step: the codes in the staged summaries were produced by a pipeline
# version whose numbering differs from the upstream STEP_LABELS constant, so keying on codes
# mislabels segments and drops the largest bucket — ≤5 molecules. Keying on the text is robust).
CURATION_OUTCOME_ORDER = [
    "retained", "not whole-cell", "qualitative-only", "≤5 molecules",
    "uncategorized", "multi-unit", "Co-ADD", "superseded by PubChem", "manually excluded",
]
CURATION_OUTCOME_COLOR = {
    "retained": hue("cobalt"),
    "not whole-cell": hue("amber"),
    "qualitative-only": hue("tangerine"),
    "≤5 molecules": hue("lime"),
    "uncategorized": hue("turquoise"),
    "multi-unit": hue("silver"),
    "Co-ADD": hue("orchid"),
    "superseded by PubChem": hue("periwinkle"),
    "manually excluded": hue("crimson"),
}
# Verbose discard_reason (as written in 21_curation_summary.csv) -> short outcome label.
_DISCARD_REASON_SHORT = {
    "manually excluded (counter-screen / flagged)": "manually excluded",
    "not ORGANISM (SINGLE PROTEIN)": "not whole-cell",
    "not ORGANISM (DISCARDED)": "not whole-cell",
    "no quantitative data (qualitative-only)": "qualitative-only",
    "multiple units (incompatible, not pooled)": "multi-unit",
    "uncategorized ((activity_type, unit) not in DR/SP)": "uncategorized",
    "Co-ADD": "Co-ADD",
    "superseded by PubChem": "superseded by PubChem",
    "<= 5 molecules": "≤5 molecules",
}

# Dataset-size bins (step 21 whole-cell sizes).
SIZE_EDGES = [5, 20, 50, 100, 1000]
SIZE_LABELS = ["≤5", "5–20", "20–50", "50–100", "100–1000", ">1000"]

# Active-ratio funnel (activity ratio by stage 22 -> 26).
STAGE_RATIO_SPECS = [  # (stage label, filename, ratio column, optional row filter)
    ("22", "22_binarisation_summary.csv", "active_ratio", None),
    ("23", "23_pool_summary.csv", "active_ratio", None),
    ("24", "24_cv_summary.csv", "active_ratio", lambda d: d[d["reason"] == "modelled"]),
    ("25", "25_pool_summary.csv", "ar_after", None),
    ("26", "26_cv_summary.csv", "active_ratio", None),
]

# Cutoff-sensitivity unit panels (one per measurement unit).
SENS_UNIT_PANELS = [
    ("umol.L-1", "Dose response (µM)"),
    ("%", "Single point (% effect)"),
    ("mm", "Single point (zone, mm)"),
]
CUTOFF_TIERS = ("low", "middle", "high")

CATS = ["DR", "SP"]
CAT_DISPLAY = {"DR": "Dose Response", "SP": "Single Point"}
CATCHALL_MIN_AUROC = 0.70   # step-26 evaluation reference line (NOT enforced)

_DR = CATEGORY_COLORS["DR"]
_SP = CATEGORY_COLORS["SP"]
_CAT_HUE = {"DR": "cobalt", "SP": "tangerine"}  # ArticleColors hue name per category
_CAT_LEGEND = {"Dose Response": _DR, "Single Point": _SP}


def _size_bins(values):
    """Bin index 0..5 for each value ('≤edge' semantics), matching plot_utils.size_bins."""
    return np.searchsorted(SIZE_EDGES, np.asarray(values, dtype=float), side="left")


def _sequential(base_name, n):
    """``n`` shades of an ArticleColors hue from a visible light tint to the full colour.

    The lightest shade is capped well above white (``lighten=0.35``, not ``0.0``) so no ordered
    bin renders as an invisible near-white block on the white page (see ``wholecell_sizes``).
    Reverse of ``plotting_colors.shades`` (full hue first, lightening towards a floor); this
    module wants the opposite direction, light-to-dark, without duplicating the ramp math."""
    if int(n) == 1:
        return [hue(base_name, lighten=0.35)]
    return list(reversed(shades(base_name, n, floor=0.35)))


# --------------------------------------------------------------------------------------------
# Pathogen labelling and house-style DR/SP box (set once, reused by every panel).
# --------------------------------------------------------------------------------------------

_SEED = None   # set in save_curation_figures from default.RANDOM_SEED
_NAMES = None  # code -> binomial name, set in save_curation_figures (abbrev is a no-op if None)


def _plabel(code):
    return abbrev(code, _NAMES)


def _labels(codes):
    return [_plabel(c) for c in codes]


def _cat_box(ax, values, position, cat, *, vert=True, width=0.34, jitter=True,
             jitter_width=0.12, point_size=6, point_alpha=0.5, rng=None):
    """House-style DR/SP distribution box: lightened category fill, INK lines, colour jitter."""
    return box_with_jitter(ax, values, position, CATEGORY_COLORS[cat],
                           face=hue(_CAT_HUE[cat], lighten=0.55), vert=vert, width=width,
                           jitter=jitter, jitter_width=jitter_width, point_size=point_size,
                           point_alpha=point_alpha, rng=rng)


# --------------------------------------------------------------------------------------------
# Data access helpers (mirror _read1 / _collect / available_pathogens from 27_general_plots.py,
# but read the staged data/raw/chembl_curation/ tree instead of the upstream stage4 tree).
# --------------------------------------------------------------------------------------------

def _read1(data_dir, code, fname):
    p = os.path.join(data_dir, code, fname)
    return pd.read_csv(p) if os.path.exists(p) else None


def _collect(data_dir, codes, fname):
    """Concatenate a per-pathogen CSV across the codes that have it."""
    frames = []
    for c in codes:
        d = _read1(data_dir, code=c, fname=fname)
        if d is not None:
            frames.append(d)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def available_pathogens(data_dir):
    """Pathogen codes present in the staged tree (a subdir with 22_binarisation_summary.csv),
    sorted for deterministic ordering."""
    codes = []
    for name in sorted(os.listdir(data_dir)):
        if name == "general":
            continue
        if os.path.exists(os.path.join(data_dir, name, "22_binarisation_summary.csv")):
            codes.append(name)
    return codes


# --------------------------------------------------------------------------------------------
# Curation attrition
# --------------------------------------------------------------------------------------------

def _stacked_outcome_barh(ax, codes, fracs, xlabel, title, legend=True):
    """Shared renderer: horizontal stacked bars per pathogen over CURATION_OUTCOME_ORDER (plus any
    unmapped outcomes appended), with the shared palette. ``fracs`` is a list (one dict per code)
    of {outcome_label: fraction}. ``legend=False`` omits the legend (use the standalone
    CurationOutcomeLegendPlot when panels are packed into a tight row)."""
    extra = [k for f in fracs for k in f
             if k not in CURATION_OUTCOME_ORDER and k not in CURATION_OUTCOME_COLOR]
    segs = CURATION_OUTCOME_ORDER + sorted(set(extra))
    stacked_hbar(ax, _labels(codes), fracs, segs, CURATION_OUTCOME_COLOR,
                 always_show=("retained",))
    if legend:
        ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5),
                  fontsize=stylia.FONTSIZE_SMALL, **LEGEND_KW)
    stylia.label(ax, xlabel=xlabel, ylabel="", title=title)


class CurationDiscardPlot(BasePlot):
    """Horizontal stacked bars of assays retained vs discarded per pathogen, by reason (fractions
    summing to 1). Outcomes come from the ``discard_reason`` text in ``21_curation_summary.csv``
    (retained = ``discarded`` is False), mapped to the shared curation-outcome taxonomy so this
    figure pairs with :class:`ChemspaceAttritionPlot`."""

    def __init__(self, data_dir, codes, ax=None, cells=(4, 4), legend=True):
        BasePlot.__init__(self, ax=ax, cells=cells)
        self.name = "curation_discard"
        fracs = []
        for c in codes:
            d = _read1(data_dir, c, "21_curation_summary.csv")
            if d is None or d.empty:
                fracs.append({})
                continue
            discarded = d["discarded"].astype(bool)
            reason = d["discard_reason"].map(_DISCARD_REASON_SHORT).fillna(
                d["discard_reason"].astype(str))
            outcome = reason.where(discarded, "retained")
            fracs.append((outcome.value_counts() / len(d)).to_dict())
        _stacked_outcome_barh(self.ax, codes, fracs, "Fraction of assays", "",
                              legend=legend)


# Chemical-space attrition. The cumulative curation stages in 21_curation_stats.csv (funnel order);
# the loss between consecutive stages is attributed to the filter applied there, mapped to the
# shared curation-outcome label so this figure uses the same labels/colours as the discard figure.
_CHEMSPACE_STAGE_ORDER = ["All", "Whole-cell", "Quantitative", "Single-unit",
                          "DR/SP", "Co-ADD", "PubChem", "≥6 mols"]
_CHEMSPACE_STAGE_LABEL = {  # destination stage -> shared curation-outcome label (loss at that stage)
    "Whole-cell": "not whole-cell",
    "Quantitative": "qualitative-only",
    "Single-unit": "multi-unit",
    "DR/SP": "uncategorized",
    "Co-ADD": "Co-ADD",
    "PubChem": "superseded by PubChem",
    "≥6 mols": "≤5 molecules",
}


class ChemspaceAttritionPlot(BasePlot):
    """Horizontal stacked bars of unique-molecule (chemical-space) attrition per pathogen: the
    fraction of the starting chemical space retained, plus the fraction removed at each curation
    stage. The chemical-space analogue of :class:`CurationDiscardPlot` (which counts assays, not
    unique molecules), sharing its outcome labels and palette. Reads the per-stage
    ``chemical_space`` column of ``21_curation_stats.csv``."""

    def __init__(self, data_dir, codes, ax=None, cells=(4, 4), legend=True):
        BasePlot.__init__(self, ax=ax, cells=cells)
        self.name = "chemspace_attrition"
        fracs = []
        for c in codes:
            d = _read1(data_dir, c, "21_curation_stats.csv")
            f = {}
            if d is not None and not d.empty:
                cs = d.set_index("stage")["chemical_space"]
                present = [s for s in _CHEMSPACE_STAGE_ORDER if s in cs.index]
                if present and float(cs[present[0]]) > 0:
                    total = float(cs[present[0]])
                    f["retained"] = float(cs[present[-1]]) / total
                    for a, b in zip(present[:-1], present[1:]):
                        lab = _CHEMSPACE_STAGE_LABEL.get(b)
                        if lab:
                            f[lab] = f.get(lab, 0.0) + (float(cs[a]) - float(cs[b])) / total
            fracs.append(f)
        _stacked_outcome_barh(self.ax, codes, fracs, "Fraction of chemical space", "",
                              legend=legend)


class CurationOutcomeLegendPlot(BasePlot):
    """Standalone shared legend for the paired curation-outcome figures (``curation_discard`` +
    ``chemspace_attrition``): one swatch per outcome in CURATION_OUTCOME_ORDER. Rendered as its
    own small panel so it can be placed once when the packed row omits per-panel legends."""

    def __init__(self, ax=None, cells=(2, 2)):
        BasePlot.__init__(self, ax=ax, cells=cells)
        self.name = "curation_outcome_legend"
        self.legend({o: CURATION_OUTCOME_COLOR[o] for o in CURATION_OUTCOME_ORDER},
                    loc="center")
        self.ax.set_axis_off()


class WholecellSizesPlot(BasePlot):
    """Horizontal stacked bars of whole-cell dataset sizes by size bin (fraction of assays),
    per pathogen. Ported from ``plot_size_stack``; sizes = ``n_mol_dataset`` for ORGANISM assays
    with numeric measurements, from ``21_curation_summary.csv``."""

    def __init__(self, data_dir, codes, ax=None, cells=(4, 4)):
        BasePlot.__init__(self, ax=ax, cells=cells)
        self.name = "wholecell_sizes"
        nb = len(SIZE_LABELS)
        mat = np.zeros((len(codes), nb))
        for i, c in enumerate(codes):
            d = _read1(data_dir, c, "21_curation_summary.csv")
            if d is None:
                continue
            wc = d[(d["target_type_curated_extra"] == "ORGANISM") & (d["n_numeric"] > 0)]
            sizes = np.asarray([v for v in wc["n_mol_dataset"].tolist()
                                if v and not pd.isna(v)], dtype=float)
            sizes = sizes[sizes > 0]
            idx = _size_bins(sizes)
            for k in range(nb):
                mat[i, k] = int((idx == k).sum())
        tot = mat.sum(axis=1, keepdims=True)
        tot[tot == 0] = 1
        mat = mat / tot
        seg_colors = dict(zip(SIZE_LABELS, _sequential("cobalt", nb)))
        fracs = [{SIZE_LABELS[k]: mat[i, k] for k in range(nb)} for i in range(len(codes))]
        stacked_hbar(self.ax, _labels(codes), fracs, SIZE_LABELS, seg_colors,
                     always_show=tuple(SIZE_LABELS))
        swatch_legend(self.ax, seg_colors, loc="lower right", title="molecules")
        self.label(xlabel="Fraction of assays", title="Whole-cell datasets by size")


# --------------------------------------------------------------------------------------------
# Binarisation & activity ratios
# --------------------------------------------------------------------------------------------

class BinarisationActiveRatioPlot(BasePlot):
    """Active-ratio boxplots (DR top · SP bottom) per pathogen, with jittered points. Ported
    from ``plot_ratio_box_grouped``; ``active_ratio`` from ``22_binarisation_summary.csv``."""

    def __init__(self, data_dir, codes, ax=None, cells=(4, 4)):
        BasePlot.__init__(self, ax=ax, cells=cells)
        self.name = "binarisation_active_ratio"
        rng = np.random.default_rng(_SEED)
        pos = np.arange(len(codes))
        for cat, dy in (("DR", -0.2), ("SP", 0.2)):
            col = "DR" if cat == "DR" else "SP"
            for i, c in enumerate(codes):
                d = _read1(data_dir, c, "22_binarisation_summary.csv")
                if d is None:
                    continue
                vals = d.loc[d["category"] == cat, "active_ratio"].dropna().tolist()
                _cat_box(self.ax, vals, pos[i] + dy, col, vert=False, width=0.32,
                         jitter_width=0.12, point_size=3, point_alpha=0.2, rng=rng)
        self.ax.set_xlim(-0.02, 1.02)
        self.ax.set_yticks(pos)
        self.ax.set_yticklabels(_labels(codes))
        self.ax.invert_yaxis()
        self.legend(_CAT_LEGEND, loc="lower right")
        self.label(xlabel="Active ratio", title="Active ratio (DR top · SP bottom)")


class ActivityRatioFlowPlot(MultiPanelPlot):
    """Active-ratio distribution at each pipeline stage (22->26), DR and SP panels, all pathogens
    pooled. Ported from ``make_activity_ratio`` (``plot_stage_ratio``); reads the 22/23/24/25/26
    summaries."""

    def __init__(self, data_dir, codes, cells=(3, 6)):
        axs = self._new_figure(1, 2, cells, "activity_ratio_flow")
        d22 = _collect(data_dir, codes, "22_binarisation_summary.csv")
        if not d22.empty:
            d22 = d22[d22["category"].isin(CATS)]
        p23 = _collect(data_dir, codes, "23_pool_summary.csv")
        c24 = _collect(data_dir, codes, "24_cv_summary.csv")
        p25 = _collect(data_dir, codes, "25_pool_summary.csv")
        c26 = _collect(data_dir, codes, "26_cv_summary.csv")
        for cat in CATS:
            self._panel(axs.next(), cat, d22, p23, c24, p25, c26)

    def _panel(self, ax, category, d22, p23, c24, p25, c26):
        def stage_col(df, value_col, reason=None):
            if "category" not in df.columns or value_col not in df.columns:
                return pd.Series([], dtype=float)
            sel = df["category"] == category
            if reason is not None:
                if "reason" not in df.columns:
                    return pd.Series([], dtype=float)
                sel = sel & (df["reason"] == reason)
            return df.loc[sel, value_col]

        stages = [
            ("22 datasets", stage_col(d22, "active_ratio")),
            ("23 pools", stage_col(p23, "active_ratio")),
            ("24 modelled", stage_col(c24, "active_ratio", reason="modelled")),
            ("25 final pools", stage_col(p25, "ar_after")),
            ("26 catch-all", stage_col(c26, "active_ratio")),
        ]
        rng = np.random.default_rng(_SEED)
        for i, (_lab, s) in enumerate(stages):
            v = np.asarray(s.dropna(), dtype=float)
            v = v[np.isfinite(v)]
            if not len(v):
                continue
            _cat_box(ax, v, i, category, vert=True, width=0.6, jitter_width=0.18,
                     point_size=8, point_alpha=0.35, rng=rng)
        ref_line(ax, 0.5, axis="y")
        ax.set_ylim(-0.02, 1.02)
        ax.set_xticks(range(len(stages)))
        ax.set_xticklabels([lab for lab, _ in stages], rotation=30, ha="right")
        stylia.label(ax, xlabel="", ylabel="active ratio",
                     title=f"Active ratio by stage — {CAT_DISPLAY[category]}")


class ActivityRatioPerPathogenPlot(MultiPanelPlot):
    """Small multiples (one panel per pathogen) of the active-ratio flow across stages 22->26,
    DR vs SP grouped boxes. Ported from ``make_activity_ratio_per_pathogen``."""

    def __init__(self, data_dir, codes, cells=(6, 6)):
        ncols = 3
        nrows = int(np.ceil(len(codes) / ncols))
        axs = self._new_figure(nrows, ncols, cells, "activity_ratio_per_pathogen")
        long = self._stage_ratio_long(data_dir, codes)
        stages = [s for s, _, _, _ in STAGE_RATIO_SPECS]
        panels = [axs.next() for _ in range(nrows * ncols)]
        for j, code in enumerate(codes):
            ax = panels[j]
            sub = long[long["pathogen"] == code]
            for cat, dx in (("DR", -0.18), ("SP", 0.18)):
                for i, st in enumerate(stages):
                    v = sub.loc[(sub["category"] == cat) & (sub["stage"] == st), "ratio"].to_numpy()
                    if not len(v):
                        continue
                    _cat_box(ax, v, i + dx, cat, vert=True, width=0.3, jitter=False)
            ref_line(ax, 0.5, axis="y", linewidth=0.6)
            ax.set_ylim(-0.02, 1.02)
            ax.set_xlim(-0.5, len(stages) - 0.5)
            ax.set_xticks(range(len(stages)))
            ax.set_xticklabels(stages)
            stylia.label(ax, xlabel="", ylabel="active ratio" if j % ncols == 0 else "",
                         title=_plabel(code))
        for ax in panels[len(codes):]:
            ax.axis("off")

    @staticmethod
    def _stage_ratio_long(data_dir, codes):
        rows = []
        for code in codes:
            for stage, fname, col, filt in STAGE_RATIO_SPECS:
                df = _read1(data_dir, code, fname)
                if df is None:
                    continue
                if filt is not None:
                    df = filt(df)
                if "category" not in df.columns or col not in df.columns:
                    continue
                for cat in CATS:
                    for v in df.loc[df["category"] == cat, col].dropna():
                        rows.append((code, cat, stage, float(v)))
        return pd.DataFrame(rows, columns=["pathogen", "category", "stage", "ratio"])


class CutoffSensitivityPlot(MultiPanelPlot):
    """One panel per measurement unit: pathogens on y, active ratio on x, three tier markers
    (low/middle/high candidate cutoff) joined by a line; the default tier is starred. Ported
    from ``plot_cutoff_sensitivity_panel``; reads ``general/27_cutoff_sensitivity.csv``."""

    def __init__(self, data_dir, codes, cells=(3, 6)):
        axs = self._new_figure(1, len(SENS_UNIT_PANELS), cells, "cutoff_sensitivity")
        df = _read1(data_dir, "general", "27_cutoff_sensitivity.csv")
        if df is None or df.empty:
            self._unavailable()
            return
        tier_colors = dict(zip(CUTOFF_TIERS, _sequential("cobalt", len(CUTOFF_TIERS))[::-1]))
        for unit, title in SENS_UNIT_PANELS:
            self._panel(axs.next(), df[df["unit"] == unit], codes, tier_colors, title)

    def _panel(self, ax, df_unit, codes, tier_colors, title):
        pos = np.arange(len(codes))
        for i, c in enumerate(codes):
            sub = df_unit[df_unit["pathogen"] == c]
            pts = []
            for t in CUTOFF_TIERS:
                row = sub[sub["tier"] == t]
                if len(row):
                    pts.append((t, float(row["active_ratio"].iloc[0]), bool(row["is_default"].iloc[0])))
            if not pts:
                continue
            ax.plot([v for _, v, _ in pts], [pos[i]] * len(pts), color=REFERENCE_LINE, zorder=1)
            for t, v, is_def in pts:
                ax.scatter([v], [pos[i]], color=tier_colors[t], zorder=3)
                if is_def:
                    ax.scatter([v], [pos[i]], color=INK, marker="*", zorder=4)
        ax.set_xlim(-0.02, 1.02)
        ax.set_yticks(pos)
        ax.set_yticklabels(_labels(codes))
        ax.invert_yaxis()
        stylia.label(ax, xlabel="Active ratio", ylabel="", title=title)


# --------------------------------------------------------------------------------------------
# Pooling & modelling
# --------------------------------------------------------------------------------------------

class PoolPartitionPlot(MultiPanelPlot):
    """Stacked bar per pathogen: molecule fraction in multi-dataset pool / singleton / not
    considered, DR and SP panels. Ported from ``plot_chemspace``; reads
    ``23_chemspace_partition.csv``."""

    def __init__(self, data_dir, codes, cells=(3, 6)):
        axs = self._new_figure(1, 2, cells, "pool_partition")
        part = _collect(data_dir, codes, "23_chemspace_partition.csv")
        if part.empty:
            self._unavailable()
            return
        for cat in CATS:
            self._panel(axs.next(), part, codes, cat)

    def _panel(self, ax, part, codes, category):
        sub = part[part["category"] == category].set_index("pathogen")
        x = np.arange(len(codes))
        f_multi, f_single, f_not = [], [], []
        for c in codes:
            if c in sub.index:
                r = sub.loc[c]
                tot = float(r["n_total"])
                if tot > 0:
                    f_multi.append(r["n_multi"] / tot)
                    f_single.append(r["n_singleton"] / tot)
                    f_not.append(r["n_not_considered"] / tot)
                    continue
            f_multi.append(0.0); f_single.append(0.0); f_not.append(1.0)
        f_multi, f_single, f_not = np.array(f_multi), np.array(f_single), np.array(f_not)
        seg = {"multi-dataset pool": hue("cobalt"), "singleton": hue("periwinkle"),
               "not considered": REFERENCE_LINE}
        ax.bar(x, f_multi, color=seg["multi-dataset pool"], label="multi-dataset pool")
        ax.bar(x, f_single, bottom=f_multi, color=seg["singleton"], label="singleton")
        ax.bar(x, f_not, bottom=f_multi + f_single, color=seg["not considered"],
               label="not considered")
        ax.set_xticks(x)
        ax.set_xticklabels(_labels(codes), rotation=90)
        ax.set_ylim(0, 1)
        swatch_legend(ax, seg, loc="upper right")
        stylia.label(ax, xlabel="", ylabel="fraction of molecules",
                     title=f"Step-23 pool partition — {CAT_DISPLAY[category]}")


def _final_pools(data_dir, codes):
    """Final pools as a long frame [pathogen, category, size, ratio, auroc] — the union of
    step-25 grown pools (size = n_mol_after, ratio = ar_after, auroc = grown_auroc) and step-26
    catch-all pools (size = n_molecules, ratio = active_ratio, auroc = auroc)."""
    frames = []
    p25 = _collect(data_dir, codes, "25_pool_summary.csv")
    if not p25.empty and {"n_mol_after", "ar_after", "grown_auroc"}.issubset(p25.columns):
        frames.append(p25[["pathogen", "category", "n_mol_after", "ar_after", "grown_auroc"]]
                      .rename(columns={"n_mol_after": "size", "ar_after": "ratio",
                                       "grown_auroc": "auroc"}))
    c26 = _collect(data_dir, codes, "26_cv_summary.csv")
    if not c26.empty and {"n_molecules", "active_ratio", "auroc"}.issubset(c26.columns):
        frames.append(c26[["pathogen", "category", "n_molecules", "active_ratio", "auroc"]]
                      .rename(columns={"n_molecules": "size", "active_ratio": "ratio"}))
    if not frames:
        return pd.DataFrame(columns=["pathogen", "category", "size", "ratio", "auroc"])
    return pd.concat(frames, ignore_index=True)


def _final_pool_codes(final, codes):
    return [c for c in codes if (final["pathogen"] == c).any()]


class PoolActiveRatiosPlot(BasePlot):
    """One dot per FINAL pool (y = active ratio, area ~ molecules), DR and SP paired per pathogen
    in a single panel (DR = cobalt, SP = tangerine) with per-group median lines. Final pools =
    step-25 grown ∪ step-26 catch-all. Single-panel, combined-DR/SP adaptation of the step-23
    ``plot_active_ratios`` bubble plot."""

    _SIZE_MIN, _SIZE_MAX = 15.0, 600.0

    def __init__(self, data_dir, codes, ax=None, cells=(3, 3)):
        BasePlot.__init__(self, ax=ax, cells=cells)
        self.name = "pool_active_ratios"
        final = _final_pools(data_dir, codes).dropna(subset=["size", "ratio"])
        final = final[final["size"] > 0]
        if final.empty:
            self._unavailable()
            return
        nmax = float(final["size"].max())
        rng = np.random.default_rng(_SEED)
        cc = _final_pool_codes(final, codes)
        for xi, c in enumerate(cc):
            for cat, dx in (("DR", -0.2), ("SP", 0.2)):
                g = final[(final["pathogen"] == c) & (final["category"] == cat)].sort_values(
                    "size", ascending=False)
                if g.empty:
                    continue
                y = g["ratio"].to_numpy(float)
                x = xi + dx + rng.uniform(-0.08, 0.08, len(y))
                self.ax.scatter(x, y, s=self._area_for(g["size"].to_numpy(), nmax),
                                color=CATEGORY_COLORS[cat], alpha=0.55, edgecolor="white",
                                linewidth=0.5, zorder=2)
                self.ax.hlines(float(np.median(y)), xi + dx - 0.16, xi + dx + 0.16,
                               color=INK, zorder=3)
        self.ax.set_ylim(-0.08, 1.08)
        self.ax.set_yticks(np.linspace(0, 1, 6))
        self.ax.set_xticks(np.arange(len(cc)))
        self.ax.set_xticklabels(_labels(cc), rotation=90)
        self.legend(_CAT_LEGEND, loc="upper right")
        # Pool-size key: drawn manually in axes-fraction coords rather than via ax.legend() —
        # matplotlib's Legend cannot reliably vertically center a label against a marker this
        # much larger than a normal legend glyph, so dot and label are placed at the same
        # explicit y instead.
        for v, y in zip((1_000, 50_000), (0.94, 0.87)):
            self.ax.scatter([0.06], [y], s=self._area_for(v, nmax), transform=self.ax.transAxes,
                            color=REFERENCE_LINE, alpha=0.55, edgecolor="white", linewidth=0.5,
                            clip_on=False)
            self.ax.text(0.11, y, f"{v:,} molecules", transform=self.ax.transAxes,
                        va="center", fontsize=stylia.FONTSIZE_SMALL)
        self.label(xlabel="", ylabel="active ratio", title="Final pool active ratio")

    def _area_for(self, n, nmax):
        return self._SIZE_MIN + (np.asarray(n, dtype=float) / nmax) * (self._SIZE_MAX - self._SIZE_MIN)


class PoolCvAurocPlot(BasePlot):
    """CV AUROC of the FINAL pools per pathogen (grown step-25 ∪ catch-all step-26), DR and SP
    paired per pathogen in a single panel (box + jittered dots, 0.70 reference line). Final-pool
    counterpart of the step-24 ``plot_auroc_box``; pairs with :class:`PoolActiveRatiosPlot`."""

    def __init__(self, data_dir, codes, ax=None, cells=(3, 3)):
        BasePlot.__init__(self, ax=ax, cells=cells)
        self.name = "pool_cv_auroc"
        final = _final_pools(data_dir, codes)
        if final.empty or final["auroc"].dropna().empty:
            self._unavailable()
            return
        rng = np.random.default_rng(_SEED)
        cc = _final_pool_codes(final, codes)
        for xi, c in enumerate(cc):
            for cat, dx in (("DR", -0.2), ("SP", 0.2)):
                v = final[(final["pathogen"] == c) & (final["category"] == cat)][
                    "auroc"].dropna().to_numpy(float)
                _cat_box(self.ax, v, xi + dx, cat, vert=True, width=0.34, jitter_width=0.1,
                         point_size=6, point_alpha=0.6, rng=rng)
        self.ref_line(0.7, axis="y")
        self.ax.set_ylim(0.4, 1.02)
        self.ax.set_xticks(np.arange(len(cc)))
        self.ax.set_xticklabels(_labels(cc), rotation=90)
        self.legend(_CAT_LEGEND, loc="lower right")
        self.label(xlabel="", ylabel="final pool CV AUROC", title="Final pool modellability")


class MergeAurocPlot(MultiPanelPlot):
    """Grouped box per pathogen: baseline (pre-merge) vs grown (post-merge) pool CV AUROC,
    DR and SP panels (0.70 guide line). Ported from ``plot_auroc_lift``; reads
    ``25_pool_summary.csv``."""

    def __init__(self, data_dir, codes, cells=(3, 6)):
        axs = self._new_figure(1, 2, cells, "merge_auroc")
        ps = _collect(data_dir, codes, "25_pool_summary.csv")
        if ps.empty:
            self._unavailable()
            return
        for cat in CATS:
            self._panel(axs.next(), ps, codes, cat)

    def _panel(self, ax, ps, codes, category):
        base_fill = hue("turquoise", lighten=0.55)
        grown_fill = hue("turquoise")
        sub = ps[ps["category"] == category]
        for xi, c in enumerate(codes):
            g = sub[sub["pathogen"] == c]
            for dx, col, fill in [(-0.2, "baseline_auroc", base_fill),
                                  (0.2, "grown_auroc", grown_fill)]:
                v = g[col].dropna().to_numpy(float)
                box_with_jitter(ax, v, xi + dx, fill, face=fill, vert=True, width=0.34,
                                jitter=False)
        ref_line(ax, 0.7, axis="y")
        ax.set_ylim(0.4, 1.02)
        ax.set_xticks(np.arange(len(codes)))
        ax.set_xticklabels(_labels(codes), rotation=90)
        stylia.label(ax, xlabel="", ylabel="pool CV AUROC",
                     title=f"Step-25 baseline → grown AUROC — {CAT_DISPLAY[category]}")


class LowDataAurocPlot(BasePlot):
    """Catch-all CV AUROC per deferred (pathogen, category), sorted, coloured by the 0.70
    reference (teal ≥0.70 / amber <0.70, not enforced). Ported from ``make_lowdata_overview``;
    reads ``26_cv_summary.csv``."""

    def __init__(self, data_dir, codes, ax=None, cells=(3, 4)):
        BasePlot.__init__(self, ax=ax, cells=cells)
        self.name = "lowdata_auroc"
        rows = []
        for c in codes:
            d = _read1(data_dir, c, "26_cv_summary.csv")
            if d is not None:
                rows += d.to_dict("records")
        df = pd.DataFrame(rows)
        if df.empty or "auroc" not in df.columns:
            self._unavailable()
            return
        df = df[df["auroc"].notna()].copy()
        if df.empty:
            self._unavailable()
            return
        df["label"] = df["pathogen"].map(_plabel) + " " + df["category"]
        df = df.sort_values("auroc", ascending=False)
        colors = [AUROC_PASS_COLORS["pass"] if a >= CATCHALL_MIN_AUROC else AUROC_PASS_COLORS["fail"]
                  for a in df["auroc"]]
        self.ax.bar(np.arange(len(df)), df["auroc"].to_numpy(), color=colors)
        self.ref_line(CATCHALL_MIN_AUROC, axis="y")
        self.ax.set_xticks(np.arange(len(df)))
        self.ax.set_xticklabels(df["label"], rotation=90)
        self.ax.set_ylim(0.5, 1.0)
        self.label(xlabel="", ylabel="catch-all CV AUROC",
                   title="Low-data catch-all pools (teal ≥0.70 / amber <0.70, reference)")


class ChemblCoveragePlot(BasePlot):
    """Donut: fraction of the bioactive ChEMBL chemical space covered by the union of all
    pathogens' cleaned molecules. Ported from ``plot_chembl_coverage``; reads the small
    ``general/27_chembl_coverage.csv`` exported by step 27 (no full ChEMBL data needed)."""

    def __init__(self, data_dir, codes, ax=None, cells=(3, 3)):
        BasePlot.__init__(self, ax=ax, cells=cells)
        self.name = "chembl_coverage"
        df = _read1(data_dir, "general", "27_chembl_coverage.csv")
        if df is None or df.empty:
            self._unavailable()
            return
        total = int(df["total_bioactive_compounds"].iloc[0])
        related = int(df.loc[df["is_union"], "n_cleaned_inchikeys"].iloc[0])
        if not total:
            self._unavailable()
            return
        rest = max(total - related, 0)
        self.ax.pie([related, rest], colors=[hue("crimson"), hue("silver", lighten=0.6)],
                    startangle=90, counterclock=False, wedgeprops={"width": 0.42})
        self.ax.text(0, 0, f"{100 * related / total:.1f}%", ha="center", va="center")
        # Short single-line title — the panel is only ~45 mm and gets overlapped with the funnel.
        self.label(title="ChEMBL coverage")


class PipelineFunnelPlot(BasePlot):
    """Consolidation funnel: datasets → folded nodes → step-23 pools → final pools (log scale
    with value labels). Ported from ``make_pipeline_summary`` panel A; reads
    ``general/27_master_table.csv`` and ``23_first_pass.csv``."""

    def __init__(self, data_dir, codes, ax=None, cells=(3, 3)):
        BasePlot.__init__(self, ax=ax, cells=cells)
        self.name = "pipeline_funnel"
        m = _read1(data_dir, "general", "27_master_table.csv")
        if m is None or m.empty:
            self._unavailable()
            return
        n_nodes = len(_collect(data_dir, codes, "23_first_pass.csv"))
        n_pool23 = (m.dropna(subset=["pool_id_23"])
                    .drop_duplicates(["pathogen", "pool_id_23"]).shape[0])
        pools = (m[m["final_status"] != "uncovered"]
                 .dropna(subset=["final_pool", "final_auroc"])
                 .drop_duplicates(["pathogen", "category", "final_pool"]))
        stages = [("Datasets", len(m)), ("Folded\nnodes", n_nodes),
                  ("Pools\n(step 23)", n_pool23), ("Final\npools", len(pools))]
        x = np.arange(len(stages))
        vals = np.array([v for _, v in stages], dtype=float)
        self.ax.set_yscale("log")
        # Draw bars from a finite positive baseline, NOT the default 0: on a log axis 0 maps to
        # -inf, so bar rectangles extend infinitely below the axes. The raster PNG clips that, but
        # the vector PDF keeps the full path — hence the huge overflowing bars seen in Illustrator.
        # Crimson to match the ChEMBL-coverage donut's covered slice (the two are overlapped).
        floor = max(vals.min() / 3.0, 1.0)
        self.ax.set_ylim(floor, vals.max() * 1.8)
        self.ax.bar(x, vals - floor, bottom=floor, color=hue("crimson"))
        for xi, v in zip(x, vals):
            self.ax.text(xi, v, f"{int(v):,}", ha="center", va="bottom",
                         fontsize=stylia.FONTSIZE_SMALL)
        self.ax.set_xticks(x)
        self.ax.set_xticklabels([s for s, _ in stages])
        self.label(xlabel="", ylabel="count (log)", title="Consolidation")


# --------------------------------------------------------------------------------------------
# Entry point
# --------------------------------------------------------------------------------------------

def save_curation_figures(data_dir, output_dir, random_seed=42, name_map=None):
    """Build every curation panel as its own Nature-sized figure and save it.

    Each panel is standalone (self-sized from its ``cells`` footprint) and written as PNG +
    vector PDF under ``output_dir``; footprints are recorded in ``figure_cells.json``. Panels
    with no available data set ``is_available = False`` and are skipped (reported at the end).

    Parameters
    ----------
    data_dir    : staged curation tree (data/raw/chembl_curation) with <pathogen>/ + general/.
    output_dir  : directory to write png/, pdf/ and figure_cells.json into.
    random_seed : seed for jitter/sampling (default RANDOM_SEED).
    name_map    : optional pathogen code -> full binomial name, so axis ticks show abbreviated
                  genus names ("A. baumannii"). When None, the raw codes are shown.
    """
    global _SEED, _NAMES
    _SEED = random_seed
    _NAMES = name_map

    codes = available_pathogens(data_dir)
    print(f"Pathogens: {len(codes)} — {', '.join(codes)}")

    plots = [
        # Row 1 — curation_discard + chemspace_attrition (60 mm each), sharing the standalone
        # legend (per-panel legends omitted).
        CurationDiscardPlot(data_dir, codes, cells=(2, 2), legend=False),
        ChemspaceAttritionPlot(data_dir, codes, cells=(2, 2), legend=False),
        CurationOutcomeLegendPlot(cells=(2, 2)),
        # Coverage donut + consolidation funnel: both ~45 mm (off-grid) and both crimson, to be
        # overlapped by hand into a single slot in Illustrator.
        ChemblCoveragePlot(data_dir, codes, cells=(1.5, 1.5)),
        PipelineFunnelPlot(data_dir, codes, cells=(1.5, 1.5)),
        WholecellSizesPlot(data_dir, codes),
        BinarisationActiveRatioPlot(data_dir, codes),
        ActivityRatioFlowPlot(data_dir, codes),
        ActivityRatioPerPathogenPlot(data_dir, codes),
        CutoffSensitivityPlot(data_dir, codes),
        PoolPartitionPlot(data_dir, codes),
        # Row 2 of the assembled figure — two 90 mm (3-cell) panels across the 180 mm width,
        # DR/SP combined per panel; final pools = step-25 grown ∪ step-26 catch-all.
        PoolActiveRatiosPlot(data_dir, codes),
        PoolCvAurocPlot(data_dir, codes),
        MergeAurocPlot(data_dir, codes),
        LowDataAurocPlot(data_dir, codes),
    ]

    footprints = {}
    saved, skipped = [], []
    for p in plots:
        if p.is_available:
            p.save(output_dir)
            footprints[p.name] = list(p.cells)
            saved.append(p.name)
        else:
            skipped.append(p.name)

    with open(os.path.join(output_dir, "figure_cells.json"), "w") as f:
        json.dump(footprints, f, indent=2)

    print(f"Saved {len(saved)} figures: {', '.join(saved)}")
    if skipped:
        print(f"Skipped (no data in staged summaries): {', '.join(skipped)}")
