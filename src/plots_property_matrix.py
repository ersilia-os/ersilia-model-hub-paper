"""Figures for the per-family property matrices (step 10).

:class:`PropertyDistributionsPlot` — one histogram per endpoint over the full reference library, as a
small-multiples grid on the 3 cm cell grid. Built for step 10's 22 physicochemical descriptors, but
family-agnostic: it takes a matrix and the endpoint names to draw.

**Why distributions and not a projection.** The abx and cytotox blocks are drawn as UMAP highlights
because "the top N most antibiotic-like compounds" is a meaningful set to locate in chemical space.
"The top N compounds by molecular weight" is not the same kind of claim — it is an arbitrary slice of
a continuous descriptor — so this block is summarised by what each descriptor's distribution over the
library actually looks like.

**Nothing is filtered.** Every selected endpoint gets a panel, including near-constant ones such as
``n_radical_electrons`` (0 for all but ~1 compound in 5,000); a flat panel is the honest rendering of
a flat column. NaNs are excluded from the histogram and counted in the panel annotation, never
imputed.
"""

import json
import os

import numpy as np

import plotting_base  # noqa: F401  (applies the stylia print/article style on import)
from plotting_base import GridPlot
from plotting_colors import INK, hue
from plotting_utils import sentence_case

#: Bins per panel. Lower than the single-panel default (120) because each panel here is roughly a
#: third of a 3 cm cell wide, where 120 bins render as noise rather than a shape.
N_BINS = 40

#: Panels per row. 5 over the 180 mm print width gives ~36 mm panels, which is the narrowest at
#: which a descriptor name still fits on one line.
COLS = 5


class PropertyDistributionsPlot(GridPlot):
    """Small-multiples grid: one full-library histogram per property endpoint.

    Each panel is titled with its ``column_name`` (not the full three-part endpoint name, which does
    not fit) and annotated with its median and, where non-zero, its NaN count. The model ID is
    dropped from the panel title deliberately — a single-model family would repeat it 22 times — so
    a caption must name the source model, or read it from the endpoint stats CSV.
    """

    def __init__(self, matrix, endpoint_names, *, name="10_physchem_distributions", cols=COLS):
        super().__init__()

        def panel(ax, endpoint, color, xlabel, ylabel):
            v = matrix[endpoint].to_numpy(dtype=float)
            finite = v[np.isfinite(v)]
            n_nan = int(len(v) - len(finite))
            if not len(finite):
                ax.text(0.5, 0.5, "no data", transform=ax.transAxes, ha="center", va="center",
                        fontsize=6, color=INK)
                ax.set_xticks([])
                ax.set_yticks([])
                ax.set_xlabel("")
                ax.set_ylabel("")
                return
            ax.hist(finite, bins=N_BINS, color=hue("cobalt"), edgecolor="none")
            title = sentence_case(endpoint.split("__", 2)[2].replace("_", " "))
            ax.set_title(title, fontsize=6, color=INK)
            note = f"med {np.median(finite):.3g}"
            if n_nan:
                note += f"\n{n_nan} NaN"
            ax.text(0.97, 0.95, note, transform=ax.transAxes, ha="right", va="top",
                    fontsize=5, color=INK)
            ax.set_yticks([])
            # Set BOTH labels unconditionally, even when empty: stylia's style seeds every axis
            # with placeholder "X-axis / Units" / "Y-axis / Units" text, so a non-edge panel that is
            # merely skipped keeps the placeholder instead of staying blank.
            ax.set_xlabel(xlabel, fontsize=6)
            ax.set_ylabel(ylabel, fontsize=6)

        self.build_grid(endpoint_names, cols=cols, name=name, panel_fn=panel,
                        edge_xlabel="value", edge_ylabel="compounds")


def save_property_distribution_figure(output_dir, matrix, endpoint_names,
                                      name="10_physchem_distributions"):
    """Draw and save the distributions grid, updating ``figure_cells.json`` in ``output_dir``."""
    if not endpoint_names:
        print(f"  [skip figure] {name}: no endpoints to draw")
        return {}
    plot = PropertyDistributionsPlot(matrix, endpoint_names, name=name)
    plot.save(output_dir)
    cells_path = os.path.join(output_dir, "figure_cells.json")
    footprints = {}
    if os.path.exists(cells_path):
        with open(cells_path) as f:
            footprints = json.load(f)
    footprints[plot.name] = list(plot.cells)
    with open(cells_path, "w") as f:
        json.dump(footprints, f, indent=2)
    print(f"  figure: {plot.name} {tuple(plot.cells)} cells")
    return footprints
