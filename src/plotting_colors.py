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

# Three distinct colour sets so the groupings never share hues:
#   - Task   : cool trio (blue / teal / orange)
#   - Source : warm trio (red / green / gold), reused in the treemap dots
#   - Output : single-hue (fuchsia) gradient, shaded by rank
# Colours anchor to stylia ArticleColors (NPG). Purple/plum and grey are avoided by convention
# (grey/silver is reserved for reference lines).

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
    """``n`` fuchsia shades, darkest first — a single-hue gradient for the Output bars
    (which are sorted by count, so darkest = most models)."""
    lightens = np.linspace(1.0, 0.4, n)
    return [_AC.get("fuchsia", lighten=float(l)) for l in lightens]


# ---------------------------------------------------------------------------
# ChEMBL data-curation figures (xx_chembl_data_curation.py)
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
