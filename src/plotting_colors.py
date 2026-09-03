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


# Lightest tint allowed for a generic shade ramp (see ``shades``). Sits between
# SUBTASK_LIGHTEN_FLOOR (0.5) and the 0.35 used by plots_chembl_curation._sequential, both bounded
# by the same "no invisible colours" rule in docs/figure_conventions.md.
SHADE_LIGHTEN_FLOOR = 0.4


def shades(name, n, *, floor=SHADE_LIGHTEN_FLOOR):
    """``n`` shades of ONE hue: the base hue first, lightening towards ``floor``.

    Use this where the things being coloured are **variants of one category** and should still read
    as that category at a glance — several assay readouts of one organism, nested series of one
    track, segments of one parent task. The tint carries "same thing, other variant", and the ramp
    direction carries an order, so the sequence must be meaningful rather than arbitrary.

    Do NOT use it for set membership or for genuinely unrelated categories: a tint of a hue reads as
    a weaker version of that hue, so e.g. an intersection layer drawn as pale crimson would read as
    "some of the crimson set" rather than as its own group. Spend a distinct hue there instead.

    ``floor`` is capped short of white on purpose; see :data:`SHADE_LIGHTEN_FLOOR`.
    """
    if int(n) <= 1:
        return [hue(name)]
    return [hue(name, lighten=None if float(l) == 1.0 else float(l))
            for l in np.linspace(1.0, floor, int(n))]

# Two distinct colour sets so the groupings never share hues:
#   - Task   : cool trio (blue / teal / orange)
#   - Source : warm trio (red / green / gold), reused in the treemap dots
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
        for s, colour in zip(subs, shades(hue_name, len(subs), floor=SUBTASK_LIGHTEN_FLOOR)):
            out[SUBTASK_DISPLAY.get(s, s)] = colour
    return out


# Subtask -> colour, keyed by the SHORT display label (SUBTASK_DISPLAY), in SUBTASK_ORDER order so
# a stacked bar keeps same-task shades adjacent.
SUBTASK_COLORS = _subtask_colors()

# Source Type -> colour (used for the Source Type bars, the treemap dots, and its legend).
# Used ONLY by the two pathogen panels (circles and voronoi), so the hues are chosen against what a
# reader sees there rather than against the hub-wide Source Type split. Ordered by dot count in those
# panels — Internal 33, External 16, Replicated 2 — with the smallest category taking lime, which reads
# as an accent against the two heavier hues and keeps a 2-dot category from disappearing. Crimson is
# deliberately NOT here any more: it is the Annotation task hue, and these panels sit near the
# task-coloured ones.
# NOTE the hub-wide ordering is the other way round (External 156 > Internal 45 > Replicated 7), so do
# not read the hue ranking as a statement about the hub — only about these panels' dots.
SOURCE_TYPE_COLORS = {
    "Internal": _AC.periwinkle,
    "External": _AC.amber,
    "Replicated": _AC.lime,
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
    # NCSA = the University of Illinois/NCSA Open Source License, an OSI-approved MIT/BSD-style
    # licence with no share-alike duty. Added 2026-08-14 with the manual metadata revision.
    "NCSA": "Permissive",
    "GPL-3.0": "Copyleft",
    "AGPL-3.0": "Copyleft",
    "LGPL-3.0": "Copyleft",
    "CC-BY-NC-ND-4.0": "Non-commercial",
    # Both added 2026-08-14 with the manual metadata revision. `Non-commercial` is a bare Airtable
    # value, not an SPDX identifier: the upstream terms forbid commercial reuse without naming a
    # standard licence. It maps to the class of the same name, which is the honest reading — a
    # reuser learns exactly as much from it as the metadata records.
    "Non-commercial": "Non-commercial",
    "CC-BY-NC-SA-4.0": "Non-commercial",
    LICENSE_MISSING: LICENSE_MISSING,
}

# Turquoise = the repo's default/positive hue, so it takes the class a reuser wants; periwinkle for
# the share-alike obligation; silver for the unknown-terms bucket, per the neutral-hue convention.
#
# Non-commercial gets FUCHSIA, the one place in this figure that hue is used. The convention
# deprioritises it because it reads as emphasis — which is exactly what is needed here: the class holds
# 4 of 218 models (a 6.6 degree wedge in the donut ring), so a well-behaved hue would simply
# disappear. Emphasis is also editorially correct, since these are the only licences in the hub that
# forbid commercial reuse. It was a single model (1.7 degrees) until the 2026-08-14 manual metadata
# revision added `Non-commercial` x2 and `CC-BY-NC-SA-4.0`; still small enough that the argument for
# fuchsia holds.
LICENSE_CLASS_COLORS = {
    "Permissive": _AC.turquoise,
    "Copyleft": _AC.periwinkle,
    "Non-commercial": _AC.fuchsia,
    LICENSE_MISSING: _AC.silver,
}


# ---------------------------------------------------------------------------
# Docker build architecture
# ---------------------------------------------------------------------------
# Airtable stores "AMD64" or "AMD64,ARM64" — there is no ARM-only build, so the split is really
# "also built for ARM" vs "x86 only". Turquoise for the dual build (the repo's positive hue, and the
# slice the panel exists to show) and periwinkle for the x86-only one: two substantive hues, since both
# are real build targets rather than one being a residual bucket — silver would have implied the
# latter. Periwinkle rather than a warning hue so the panel does not read as flagging x86-only builds
# as a problem, and it pairs with turquoise without competing.
# Shortened from AMD64 / ARM64. As a donut legend row with its count, "AMD64+ARM64 129" is the widest
# label of the three donut panels and squeezes its own axes: it left the ring 18.7 mm against the other
# two panels' 19.8 mm, i.e. the ring's size was being set by the length of a word. Dropping the "64"
# clears it with room to spare. The full names are the Docker platform identifiers (`linux/amd64`,
# `linux/arm64`) and belong in a caption if precision matters there.
ARCH_DISPLAY = {"AMD64,ARM64": "AMD + ARM", "AMD64": "AMD only"}
# Cobalt / tangerine rather than the turquoise + periwinkle this panel used as a pie: those two are
# now the licence donut's, and the three donuts are meant to be read side by side, so no hue may mean
# two different things across them. Both are substantive hues — silver would have cast x86-only as a
# residual bucket rather than a real build target. Tangerine is used here as a plain categorical hue
# (as it is for Single Point in script 02), NOT as a warning: a caption must not read it as flagging
# x86-only builds as a problem.
# Ordered base-capability first, so the ring reads "x86 only, then also ARM" rather than by size —
# this dict's order is what drives the wedge and legend order.
ARCH_COLORS = {
    "AMD only": _AC.tangerine,
    "AMD + ARM": _AC.cobalt,
}

# Biomedical Area groups: ONE hue (the Annotation crimson, which is what every model in that panel is)
# differentiated by fill pattern instead of by colour. Solid for the largest group and progressively
# lighter-inked patterns after it, so the ink ordering matches the size ordering; the catch-all takes
# the cross-hatch, which reads as "mixed". Patterns are drawn in white over the crimson (matplotlib
# hatches use the patch edge colour), so every wedge still reads as red at a glance.
# Repeat counts are the density knob: in matplotlib a repeated hatch character packs the motif tighter,
# so these are deliberately long. At three or four repeats the ring showed finger-thick stripes and
# dots the size of the ring's own thickness, which read as damage rather than as a fill. Paired with
# the thin ``hatch.linewidth`` set in plotting_base.
#
# The ink ordering holds EXCEPT for the mirrored pair (2026-08-07). ``Antifungal`` was split out of
# ``Antimicrobial`` and takes the backslash, mirroring ADMET's forward slash, so the two read as a
# related pair rather than as a new independent category. Its ink therefore matches ADMET's despite
# it being the smallest substantive group (3 models vs 29) — a deliberate exception, chosen because
# every pattern lighter than the dots aliased against the ring edge on a wedge that thin.
BIOAREA_GROUP_HATCH = {
    "Antimicrobial": "",
    "ADMET": "/////////",
    "Antifungal": "\\\\\\\\\\\\\\\\\\",
    "Antiviral": ".........",
    "Other": "xxxxxxxxx",
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


#: Pathogen pairs whose positional hues are exchanged after the fact. The positional palette (see
#: :func:`pathogen_activity_colors`) happened to give *S. aureus* the lime and *E. coli* the
#: periwinkle that :data:`SHARED_ORGANISM_COLORS` assigns the other way round — the swap makes those
#: two organisms read as ONE colour across step 03 and step 05 instead of trading hues between them,
#: which is the single most confusing way two palettes can disagree. Applied to BOTH the panel and
#: :func:`pathogen_activity_colors` through :func:`swap_pathogen_hues`, so the two cannot drift.
PATHOGEN_HUE_SWAPS = (("ecoli", "saureus"),)


def swap_pathogen_hues(order, colors):
    """``colors`` with each :data:`PATHOGEN_HUE_SWAPS` pair's two entries exchanged.

    ``order`` is the positional pathogen list the colours were generated for; a pair with either
    member absent is skipped, so this is a no-op on a subset that does not contain both.
    """
    out = list(colors)
    index = {p: i for i, p in enumerate(order)}
    for a, b in PATHOGEN_HUE_SWAPS:
        ia, ib = index.get(a), index.get(b)
        if ia is not None and ib is not None:
            out[ia], out[ib] = out[ib], out[ia]
    return out


def pathogen_activity_colors(dataset_counts, *, levels=(None, 0.78)):
    """The hue ``pathogen_activity_ratios`` (step 03) gives each pathogen, as a ``{pathogen: colour}``.

    That panel colours **positionally** — ``distinct_colors(n)`` over pathogens ranked by how many
    ChEMBL datasets each has, descending, then :func:`swap_pathogen_hues` — so there is no fixed hue
    per pathogen to import, and a panel in another step can only match it by reproducing the ranking.
    This does that, from the ``{pathogen: n_datasets}`` mapping the step-03 summary CSV
    (``dataset_sizes.csv``) yields.

    Derived, not frozen, on purpose: the ranking moves whenever the ChEMBL curation adds or drops a
    dataset, and a hardcoded copy would silently stop matching the panel it exists to match.

    ``levels`` defaults to the panel's ACCENT palette (the base hue, used there for dot outlines and
    mean bars) rather than its pale ``(0.62, 0.38)`` dot fills — a 0.62-lightened hue is too weak to
    carry a mark on white.

    **These hues still do not agree with** :data:`SHARED_ORGANISM_COLORS` in general, and cannot: that
    dict is a fixed 7-organism palette, this one is a position in a 15-pathogen ranking. Only the
    :data:`PATHOGEN_HUE_SWAPS` pair is reconciled (*E. coli* lime, *S. aureus* periwinkle in both);
    every other organism may still differ, so do not mix the two palettes in one figure.
    """
    order = sorted(dataset_counts, key=lambda p: (-dataset_counts[p], p))
    return dict(zip(order, swap_pathogen_hues(order, distinct_colors(len(order), levels=levels))))

# Neutral colour for reference marks (chance diagonals, baselines, gridline emphasis).
REFERENCE_LINE = _AC.silver

# Structural ink for box/whisker outlines, median lines and marker emphasis.
INK = _AC.black


# ---------------------------------------------------------------------------
# Hub timeline (01b_community_stats.py)
# ---------------------------------------------------------------------------
# One hue per track of the five-track shared-axis timeline. Here colour is DECORATIVE, not
# semantic: each track is already named by its y label and no two tracks share a scale, so the
# hue only helps the eye keep its place while reading down a 15 mm band.
#
# That is exactly why this is its own palette rather than a reuse of TASK_COLORS. Those encode
# categories; borrowing them here would imply that the "models" track and the Annotation task are
# the same thing — a relationship that does not exist across tracks.
#
# Hue NAMES rather than resolved colours, because the People track needs a second, lighter weight
# of its own hue and ``hue()`` is the only sanctioned way to ask for one — ``distinct_colors()``
# hands back RGB tuples that cannot be lightened without reimplementing the blend.
#
# Assigned explicitly rather than by zipping _CATEGORICAL_HUES: People takes turquoise, which the
# pick order would have given to a fifth track. Turquoise is the repo's default/positive hue, and
# spending it on a decorative track is only acceptable because no track here encodes a category —
# see the note above.
TIMELINE_TRACK_ORDER = ["Models", "People", "Commits", "Issues"]
TIMELINE_HUES = {
    "Models": "cobalt",
    "People": "turquoise",
    "Commits": "lime",
    "Issues": "tangerine",
}
TIMELINE_COLORS = {t: _AC.get(h) for t, h in TIMELINE_HUES.items()}

#: Lighter tint of a track's own hue, for a track carrying more than one nested series. Unused
#: while every track is a single series, but the mechanism has to exist for the moment one is not:
#: nested areas need two weights of one hue, never two hues.
TIMELINE_SECONDARY_LIGHTEN = 0.45

#: Fill opacity under each track's line. 1.0 — the fills are the FULL hue, not a wash. At four
#: tracks of ~5 mm the bands are small enough that saturated colour reads as four clean stripes
#: rather than as noise, and a tinted fill at this size is barely distinguishable from white.
TIMELINE_FILL_ALPHA = 1.0


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


def count_shades(values, name="turquoise", *, log=False, headroom=0.15):
    """Fading shades over an ordinal scale — pale = low, saturated = high.

    The ordinal counterpart of :func:`auroc_shades`. Pass whichever variable the gradient is meant to
    encode; the ramp always runs pale at the minimum to saturated at the maximum, so reversing the
    direction is a matter of what you hand it, not a flag here.

    ``log=True`` fits on ``log10`` instead of the raw values, for quantities running over decades — a
    linear fit there puts everything below the largest value into the palest step or two and the
    colour stops distinguishing anything. Values are clamped at 1 first, so a zero shades as the
    floor rather than raising. Leave it off for a short ordinal range such as 1-7, where log would
    bunch the top of the scale together instead.

    ``headroom`` extends the fit floor below the smallest value by that fraction of the observed
    span, for the same reason ``auroc_shades`` anchors below chance: fitted exactly at the minimum,
    the smallest mark comes out white and disappears against the panel. A fraction rather than an
    absolute so it means the same thing on both scales. It is a legibility margin, not a data
    threshold.

    Colour is ordinal here, so use it only where the position scale already carries the value — a
    redundant, scannable encoding, never the sole one.
    """
    from stylia import FadingColormap

    v = np.asarray(values, dtype=float)
    if log:
        v = np.log10(np.maximum(v, 1.0))
    lo, hi = float(v.min()), float(v.max())
    span = (hi - lo) or 1.0
    cm = FadingColormap(name)
    cm.fit(np.array([lo - headroom * span, hi]))
    return cm.transform(v)
