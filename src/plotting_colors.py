"""Single source of truth for the metadata figures' semantic colours.

Anchored to stylia's ``ArticleColors`` (the non-branded NPG palette used with
``stylia.set_style("article")``), mirroring ``zairachem/report/colors.py``. Every plot
pulls its colours from here so the palette can't drift across panels.

Convention:
- Tasks: Annotation = crimson, Representation = amber, Sampling = lime.
- Subtasks: shades of their parent task's hue (see ``SUBTASK_COLORS``).
- Source Type: External = crimson, Internal = lime, Replicated = amber.
- Default bar colour: cobalt.

Note the task and Source Type trios draw on the SAME three hues. That is deliberate but only safe
because the two are never encoded by colour in the same panel: in the stacked Source Type panels
the source type is encoded by bar position and only the task/subtask segments carry colour.
"""

import numpy as np
from stylia.colors import ArticleColors

from default import SUBTASK_DISPLAY, SUBTASK_PARENT

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

# Task -> stylia hue NAME. Both TASK_COLORS and SUBTASK_COLORS derive from this single mapping,
# so a task and its subtasks can never drift onto different hues.
TASK_HUES = {
    "Annotation": "crimson",
    "Representation": "amber",
    "Sampling": "lime",
}

# Task -> colour. Subtasks inherit their parent task's colour (see SUBTASK_PARENT in default.py).
TASK_COLORS = {t: _AC.get(h, lighten=None) for t, h in TASK_HUES.items()}

# Within-task display order for the subtask segments of a stacked bar: largest subtask first, so
# the darkest shade is also the biggest slice. Keys are the RAW subtask names (as in default.py);
# SUBTASK_COLORS below is re-keyed to the shorter display labels that reach the axes and legends.
SUBTASK_ORDER = [
    "Activity prediction", "Property calculation or prediction",   # Annotation
    "Featurization", "Projection",                                 # Representation
    "Similarity search", "Generation",                             # Sampling
]

# Lightest tint allowed for a subtask shade. Capped well short of white so the smallest segment of
# a stacked bar stays visible on the page (the "no invisible colours" rule in the figure
# conventions) and still reads as the same hue as its parent task.
SUBTASK_LIGHTEN_FLOOR = 0.5


def _subtask_colors():
    """Subtask (display label) -> a shade of its parent task's hue.

    A bar segmented by subtask should still read as its TASK at a glance, so each task's hue is
    split into as many shades as it has subtasks (base hue for the largest, lightening towards
    ``SUBTASK_LIGHTEN_FLOOR``) instead of spending six unrelated hues on six subtasks — which
    would collide with the Source Type and Output palettes and lose the task grouping entirely.
    """
    out = {}
    for task, hue_name in TASK_HUES.items():
        subs = [s for s in SUBTASK_ORDER if SUBTASK_PARENT[s] == task]
        lightens = (np.linspace(1.0, SUBTASK_LIGHTEN_FLOOR, len(subs)) if len(subs) > 1
                    else [1.0])
        for s, l in zip(subs, lightens):
            lighten = None if float(l) == 1.0 else float(l)
            out[SUBTASK_DISPLAY.get(s, s)] = _AC.get(hue_name, lighten=lighten)
    return out


# Subtask -> colour, keyed by the SHORT display label (SUBTASK_DISPLAY), in SUBTASK_ORDER order so
# a stacked bar keeps same-task shades adjacent.
SUBTASK_COLORS = _subtask_colors()

# Source Type -> colour (used for the Source Type bars, the treemap dots, and its legend).
SOURCE_TYPE_COLORS = {
    "External": _AC.crimson,
    "Internal": _AC.lime,
    "Replicated": _AC.amber,
}

# Default colour for a plain bar chart.
BAR_DEFAULT = _AC.cobalt


# The catch-all value both Biomedical Area and Target Organism use for "not specific to one area /
# organism". It is the only value in either field that is not made up purely of Annotation models:
# Biomedical Area "Any" is 26 Annotation / 58 Representation / 19 Sampling, Target Organism "Any" is
# 40 / 57 / 19, while every *named* area and organism is 100% Annotation — with one exception,
# Homo sapiens, which carries 1 Representation model among its 36 Annotation ones (2.7%).
CATCH_ALL_LABEL = "Any"


def catchall_colors(values):
    """Annotation hue per value, with the ``Any`` catch-all row in the neutral hue.

    Because the named rows of these two fields are (bar one model) all Annotation, the bars carry no
    task information of their own — so they take the Annotation hue rather than a palette of their
    own, and a reader who has learnt the task colours from the task/subtask panels reads these two
    for free. ``Any`` is a catch-all bucket, which is exactly what silver is reserved for.
    """
    return [_AC.silver if v == CATCH_ALL_LABEL else TASK_COLORS["Annotation"] for v in values]


# Categorical hues in pick order, arranged so CONSECUTIVE entries are far apart in hue — a
# categorical axis usually orders its groups meaningfully, so neighbours are what a reader compares.
# black/silver/white are excluded: silver is the reserved neutral and black is structural ink.
_CATEGORICAL_HUES = ["cobalt", "crimson", "lime", "tangerine", "turquoise",
                     "periwinkle", "amber", "fuchsia", "orchid"]


def distinct_colors(n, *, levels=(None, 0.55)):
    """``n`` categorical colours, as visually distinct as ArticleColors allows.

    There are only **9** substantive hues in the palette, so beyond that this cycles back through them
    at the next entry in ``levels``, giving ``9 x len(levels)`` before anything repeats. A hue and its
    tint are the most confusable pair produced, and they land 9 apart to keep them off each other.

    ``levels`` is the lighten value per pass (``None`` = the base hue), and it is what lets a caller
    ask for a *lighter* palette without collapsing it: passing a single value for everything would
    make entry 10 identical to entry 1. Two distinct levels are always needed, e.g. ``(0.62, 0.38)``
    for a pale fill palette against ``(None, 0.55)`` for its matching darker accents — same hue per
    index, so a fill and its accent read as one category.

    Because of the 9-hue limit, only use this where colour is a **secondary, scannable** encoding
    backed by a real key, such as a labelled categorical axis. It is not adequate as the sole encoding
    for more than 9 categories, and a legend of 15 swatches would not be readable at panel size.
    """
    hues = _CATEGORICAL_HUES
    out = []
    for i in range(int(n)):
        level = levels[(i // len(hues)) % len(levels)]
        out.append(_AC.get(hues[i % len(hues)], lighten=level))
    return out


def ordinal_shades(n, lo=0.4):
    """``n`` periwinkle shades, darkest first — this figure's ORDINAL single-hue gradient.

    Used by the tag cloud, so a gradient always reads as "rank/count within one field". Note
    periwinkle now does **double duty**: it is also the flat categorical hue for licence Copyleft and
    for x86-only builds. The two never meet in a panel and a gradient does not read as a category, but
    if that ever needs separating, cobalt is free and is the obvious hue to move this to. ``lo`` caps
    the lightest tint short of white.
    """
    lightens = np.linspace(1.0, lo, n)
    return [_AC.get("periwinkle", lighten=float(l)) for l in lightens]


# ---------------------------------------------------------------------------
# Licence reuse classes (01_ersilia_metadata.py)
# ---------------------------------------------------------------------------
# Value used for models with no licence recorded in Airtable. Not a licence — the terms are simply
# unknown, which for a reuser is its own (worse) category, so it is never folded into a real one.
LICENSE_MISSING = "Not recorded"

# Licence -> reuse class. Keys are the SIMPLIFIED identifiers (the "-or-later" / "-only" suffixes are
# collapsed upstream in the script), so GPL-3.0-or-later and GPL-3.0-only both arrive as "GPL-3.0".
# Collapsing loses a real legal distinction and is a presentation simplification only; the raw values
# stay in data/raw and in license_counts.csv.
#
# The permissive/copyleft split is a coarse reuse classification, NOT legal advice: CC-BY-4.0 and
# CC0-1.0 are grouped as permissive because they place no share-alike duty on a reuser, though
# neither is an OSI-approved software licence and both are poor fits for code.
LICENSE_CLASS = {
    "MIT": "Permissive",
    "Apache-2.0": "Permissive",
    "BSD-3-Clause": "Permissive",
    "CC0-1.0": "Permissive",
    "CC-BY-4.0": "Permissive",
    "GPL-3.0": "Copyleft",
    "AGPL-3.0": "Copyleft",
    "LGPL-3.0": "Copyleft",
    "CC-BY-NC-ND-4.0": "Non-commercial",
    LICENSE_MISSING: LICENSE_MISSING,
}

# Turquoise = the repo's default/positive hue, so it takes the class a reuser wants; periwinkle for
# the share-alike obligation; silver for the unknown-terms bucket, per the neutral-hue convention.
#
# Non-commercial gets FUCHSIA, the one place in this figure that hue is used. The convention
# deprioritises it because it reads as emphasis — which is exactly what is needed here: the class holds
# a single model (1/208 = a 1.7 degree wedge in the pie, a hairline bar in the chart), so a
# well-behaved hue would simply disappear. Emphasis is also editorially correct, since it is the only
# licence in the hub that forbids commercial reuse.
LICENSE_CLASS_COLORS = {
    "Permissive": _AC.turquoise,
    "Copyleft": _AC.periwinkle,
    "Non-commercial": _AC.fuchsia,
    LICENSE_MISSING: _AC.silver,
}


def license_colors(values):
    """Reuse-class colour per simplified licence identifier (see :data:`LICENSE_CLASS`)."""
    return [LICENSE_CLASS_COLORS[LICENSE_CLASS[v]] for v in values]


# ---------------------------------------------------------------------------
# Docker build architecture
# ---------------------------------------------------------------------------
# Airtable stores "AMD64" or "AMD64,ARM64" — there is no ARM-only build, so the split is really
# "also built for ARM" vs "x86 only". Turquoise for the dual build (the repo's positive hue, and the
# slice the panel exists to show) and periwinkle for the x86-only one: two substantive hues, since both
# are real build targets rather than one being a residual bucket — silver would have implied the
# latter. Periwinkle rather than a warning hue so the panel does not read as flagging x86-only builds
# as a problem, and it pairs with turquoise without competing.
ARCH_DISPLAY = {"AMD64,ARM64": "AMD64 + ARM64", "AMD64": "AMD64 only"}
ARCH_COLORS = {
    "AMD64 + ARM64": _AC.turquoise,
    "AMD64 only": _AC.periwinkle,
}


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
