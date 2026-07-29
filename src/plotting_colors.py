"""Single source of truth for the metadata figures' semantic colours.

Anchored to stylia's ``ArticleColors`` (the non-branded NPG palette used with
``stylia.set_style("article")``), mirroring ``zairachem/report/colors.py``. Every plot
pulls its colours from here so the palette can't drift across panels.

Convention (kept from the original script 01):
- Tasks: Annotation = cobalt, Representation = turquoise, Sampling = tangerine.
- Source Type: External = crimson, Internal = lime, Replicated = amber.
- Default bar colour: cobalt.
"""

import numpy as np
from stylia.colors import ArticleColors

_AC = ArticleColors()


def hue(name, lighten=None):
    """The single accessor for a raw ArticleColors hue (optionally lightened).

    Every module pulls colours from here (semantic dicts below, or ``hue`` for a
    one-off) so nothing instantiates its own ``ArticleColors``/``NamedColors`` and the
    palette can't drift. ``lighten`` follows stylia's convention: ``None`` (or 1.0) is the
    base hue, smaller values lighten towards white.
    """
    return _AC.get(name, lighten=lighten)

# Three distinct colour sets so the groupings never share hues:
#   - Task   : cool trio (blue / teal / orange)
#   - Source : warm trio (red / green / gold), reused in the treemap dots
#   - Output : single-hue (periwinkle) gradient, shaded by rank
# Colours anchor to stylia ArticleColors (NPG). Two palette rules:
#   - silver is the NEUTRAL hue: reference marks (chance diagonals, baselines) and neutral / "other"
#     categorical buckets (e.g. TARGET_TYPE_COLORS["Other"], the multi-unit curation outcome). It is
#     never used for a substantive category that carries its own meaning.
#   - fuchsia comes LAST in the pick order — it is far stronger than the rest of the palette and
#     reads as emphasis, so reach for the other hues first. It is NOT off-limits: using it is fine
#     once the other hues are taken (a categorical set needing them all), or when asked for.

# Task -> colour. Subtasks inherit their parent task's colour (see SUBTASK_PARENT in default.py).
TASK_COLORS = {
    "Annotation": _AC.cobalt,
    "Representation": _AC.turquoise,
    "Sampling": _AC.tangerine,
}

# Source Type -> colour (used for the Source Type bars, the treemap dots, and its legend).
SOURCE_TYPE_COLORS = {
    "External": _AC.crimson,
    "Internal": _AC.lime,
    "Replicated": _AC.amber,
}

# Default colour for a plain bar chart.
BAR_DEFAULT = _AC.cobalt


def output_colors(n):
    """``n`` periwinkle shades, darkest first — a single-hue gradient for the Output bars
    (which are sorted by count, so darkest = most models). Periwinkle rather than fuchsia (too
    strong, deprioritised) and it collides with neither the Task nor the Source trio."""
    lightens = np.linspace(1.0, 0.4, n)
    return [_AC.get("periwinkle", lighten=float(l)) for l in lightens]


# ---------------------------------------------------------------------------
# ChEMBL data-curation figures (02_chembl_data_curation.py)
# ---------------------------------------------------------------------------
# Semantic keys ported from ``../chembl-antimicrobial-tasks/src/plot_colors.py`` but anchored
# to stylia ArticleColors and kept within the repo convention (turquoise default, silver for
# reference/neutral, no plum/grey). DR = dose-response, SP = single-point assays.

CATEGORY_COLORS = {
    "DR": _AC.cobalt,       # dose-response datasets
    "SP": _AC.tangerine,    # single-point datasets
}

# Whether a dataset / pool survives curation or is dropped along the way.
CURATION_STATUS_COLORS = {
    "survivor": _AC.cobalt,
    "lost": _AC.amber,
}

# Assay target type.
TARGET_TYPE_COLORS = {
    "Whole-cell": _AC.turquoise,
    "Protein": _AC.amber,
    "Other": _AC.silver,
}

# Binary activity label.
ACTIVE_INACTIVE_COLORS = {
    "active": _AC.crimson,
    "inactive": _AC.silver,
}

# Pass/fail against the (inherited) 0.70 AUROC bar used in the low-data catch-all figure.
AUROC_PASS_COLORS = {
    "pass": _AC.turquoise,
    "fail": _AC.amber,
}

# One distinct hue per "shared" organism (the 7 with an EU OpenScreen primary assay), for figures
# that break a total down BY pathogen. Ordered as SHARED_ORGANISMS in default.py. Silver is excluded
# because every organism here is a substantive category (silver is the neutral hue). Fuchsia is not
# needed either — 7 organisms fit without reaching that far down the pick order, so periwinkle takes
# the seventh slot.
SHARED_ORGANISM_COLORS = {
    "abaumannii": _AC.cobalt,
    "calbicans": _AC.turquoise,
    "ecoli": _AC.lime,
    "efaecium": _AC.amber,
    "kpneumoniae": _AC.tangerine,
    "paeruginosa": _AC.crimson,
    "saureus": _AC.periwinkle,
}

# Neutral colour for reference marks (chance diagonals, baselines, gridline emphasis).
REFERENCE_LINE = _AC.silver

# Structural ink for box/whisker outlines, median lines and marker emphasis.
INK = _AC.black


# ---------------------------------------------------------------------------
# ChEMBL model-performance figures (03_chembl_models_performance.py)
# ---------------------------------------------------------------------------


def auroc_shades(values, lo=0.35, hi=1.0):
    """Cobalt shades encoding AUROC, fitted to the ``[lo, hi]`` range.

    Used to colour the per-model ROC curves so a whole grid reads at a glance:
    pale = near chance, saturated = strong ranking. ``lo`` sits below the 0.5
    chance level on purpose — fitted at exactly 0.5 a chance-level curve comes out
    white and disappears against the panel.
    """
    from stylia import FadingColormap

    cm = FadingColormap("cobalt")
    cm.fit(np.array([lo, hi]))
    return cm.transform(np.asarray(values, dtype=float))
