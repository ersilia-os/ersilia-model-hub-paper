"""Step 06 figures — CoAdd validation of the public ChEMBL pathogen models.

Four panels, all built on the step-05 classes in :mod:`plots_euopenscreen` so the CoAdd figures
are visually interchangeable with their EU OpenScreen counterparts — only the input frame and the
output name differ:

  - ``coadd_shared_auroc``          own-strain AUROC per organism, both endpoints
  - ``coadd_overlap``               training-set overlap, library + actives (analog of ``euos_overlap``)
  - ``coadd_hit_promiscuity``       actives by number of reference strains hit
  - ``coadd_hit_exclusivity_auroc`` exclusive vs shared hit AUROC, deduplicated

The three cross-organism panels read the ``COADD_HITSET_ENDPOINT`` (inhib_50) slice only; see that
constant for why. Unlike step 05 this step keeps a flat output directory, so the dedup-only
exclusivity panel sits alongside the rest — its CSV carries ``set=dedup`` on every row.
Colours come only from :mod:`plotting_colors`.
"""

import json
import os

import pandas as pd

from default import COADD_HITSET_ENDPOINT
from plotting_colors import CATEGORY_COLORS
from plots_euopenscreen import (
    EXCLUSIVITY_COLORS,
    SMALL_SQUARE,
    EuosOverlapTwinPlot,
    HitExclusivityPlot,
    HitPromiscuityPlot,
    KeyPanel,
    MetricByOrganismPlot,
    _pretty,
    euos_overlap_handles,
)

# CoAdd endpoints: mic_10 ~ dose-response (DR), inhib_50 ~ single-point (SP).
ENDPOINT_COLORS = {"mic_10": CATEGORY_COLORS["DR"], "inhib_50": CATEGORY_COLORS["SP"]}


def _hitset_slice(df):
    """The ``COADD_HITSET_ENDPOINT`` rows of a frame that carries an ``endpoint`` column.

    The leakage report has one row per (organism, endpoint); the cross-organism panels are
    single-endpoint, so they would otherwise draw each organism twice.
    """
    if df is None or df.empty or "endpoint" not in df.columns:
        return df
    return df[df["endpoint"] == COADD_HITSET_ENDPOINT]


class CoaddSharedAurocPlot(MetricByOrganismPlot):
    """CoAdd own-strain AUROC per shared organism, both endpoints (dedup)."""

    def __init__(self, coadd_df, ax=None, cells=(3, 3)):
        super().__init__(
            coadd_df, group_col="endpoint", group_colors=ENDPOINT_COLORS,
            group_order=["inhib_50", "mic_10"], metric="auroc", xlabel="AUROC",
            title="", name="coadd_shared_auroc",
            ref=0.5, xlim=(0, 1), prefer_best=True, cells=cells, ax=ax)


class CoaddOverlapTwinPlot(EuosOverlapTwinPlot):
    """Training-set overlap for the CoAdd reference-strain library and its actives.

    Same twin-axis construction as ``euos_overlap`` (hatched upper bar = library on the top axis,
    solid lower bar = actives on the bottom axis, each stacked novel vs in-training), restricted to
    the single hit-set endpoint so each organism appears once.
    """

    def __init__(self, leak_df, ax=None, cells=(3, 4), legend=True):
        super().__init__(_hitset_slice(leak_df), ax=ax, cells=cells, legend=legend)
        self.name = "coadd_overlap"


class CoaddHitPromiscuityPlot(HitPromiscuityPlot):
    """CoAdd actives by how many reference strains they hit (label-only)."""

    def __init__(self, prom_df, ax=None, cells=(3, 3)):
        super().__init__(prom_df, ax=ax, cells=cells)
        self.name = "coadd_hit_promiscuity"


class CoaddHitExclusivityPlot(HitExclusivityPlot):
    """Exclusive vs shared hit AUROC per organism on CoAdd (dedup)."""

    def __init__(self, excl_df, ax=None, cells=(3, 3), legend=True):
        super().__init__(excl_df, ax=ax, cells=cells, legend=legend)
        self.name = "coadd_hit_exclusivity_auroc"


def _read(output_dir, fname):
    path = os.path.join(output_dir, fname)
    return pd.read_csv(path) if os.path.exists(path) else pd.DataFrame()


def save_coadd_figures(output_dir):
    """Build the step-06 panels from the summary CSVs in ``output_dir`` and record their
    footprints in ``figure_cells.json``.

    The two paper panels that carry a legend are built at ``SMALL_SQUARE`` with ``legend=False``
    and followed by a standalone ``*_key``, mirroring how step 05 lays out its 45 mm row.
    """
    auroc = _read(output_dir, "06_coadd_auroc.csv")
    leak = _read(output_dir, "06_coadd_leakage_report.csv")
    promiscuity = _read(output_dir, "06_coadd_hit_promiscuity.csv")
    exclusivity = _read(output_dir, "06_coadd_hit_exclusivity.csv")

    plots = [
        CoaddSharedAurocPlot(auroc, cells=(3, 3)),
        CoaddOverlapTwinPlot(leak, cells=SMALL_SQUARE, legend=False),
        KeyPanel("coadd_overlap_key", handles=euos_overlap_handles(), ncol=2, cells=(0.55, 1.5)),
        CoaddHitPromiscuityPlot(promiscuity, cells=SMALL_SQUARE),
        CoaddHitExclusivityPlot(exclusivity, cells=SMALL_SQUARE, legend=False),
        KeyPanel("coadd_hit_exclusivity_auroc_key",
                 {_pretty(g): EXCLUSIVITY_COLORS[g] for g in ("exclusive", "nonexclusive")},
                 ncol=2, cells=(0.4, 1.5)),
    ]

    footprints = {}
    for p in plots:
        if p.is_available:
            p.save(output_dir)
            footprints[p.name] = list(p.cells)
            print(f"  figure: {p.name}")
        else:
            print(f"  [skip figure] {p.name}: no data")
    with open(os.path.join(output_dir, "figure_cells.json"), "w") as f:
        json.dump(footprints, f, indent=2)
