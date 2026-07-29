# Figure conventions

How publication figures are built in this repo. The goal: every chart is an **individual,
Illustrator-ready panel** sized on a fixed grid, so panels drop into a Nature-style page without
rescaling. Reference implementation: **`src/plotting_base.py`** (the `BasePlot` base class and grid
constants). Concrete examples: `src/plots_metadata.py`, `src/plots_chembl_curation.py`.

## Core rules

1. **stylia only.** All figures use [stylia](https://github.com/ersilia-os/stylia)
   (`stylia.create_figure`, `stylia.label`, `stylia.save_figure`) — never raw matplotlib figure
   constructors or `plt.savefig`. See the `/stylia-plotting` skill for the API. Presets are applied
   once at import in `plotting_base.py`: `set_format("print")` + `set_style("article")`.
2. **Individual panels.** One chart per file. No composite "overview" dashboards, and **no A/B/C
   panel letters** — final ordering and lettering happen in Illustrator.
3. **Both PNG and vector PDF.** Every figure is saved as a raster `png/<name>.png` (preview) *and* a
   vector `pdf/<name>.pdf` (for Illustrator), into `png/` and `pdf/` subfolders of the output dir.
   *(This overrides the generic stylia skill's "PNG by default" — here PDF is always produced.)*
4. **Cell-grid sizing.** Panels are sized on a 3 cm square-cell grid so they compose predictably.
5. **`figure_cells.json` manifest.** Each figure's `(rows, cols)` footprint is written to
   `output/<script>/figure_cells.json` so the intended grid size survives the tight-crop on save.
6. **Summary CSVs, not raw data.** A figure's data comes from a small pre-aggregated CSV. If the
   summary doesn't exist upstream, create it in the source repo and stage only that through
   `scripts/00_download_data.py` — never copy full per-molecule datasets here.

## The cell grid

stylia's `print` format is `SIZE = 7.09" = 180 mm` (Nature two-column width), and
`create_figure(width=w, height=h)` treats `w`/`h` as **fractions of SIZE**. Splitting that width
into `CELLS_PER_WIDTH = 6` gives true **3 cm square cells** (`plotting_base.py`):

```
width_mm  = cols / 6 × 180 mm
height_mm = rows / 6 × 180 mm
```

So a panel declares a footprint `cells=(rows, cols)` and renders at that physical size. A full
Nature page (180 × 215 mm) is **6 cells wide × ~7.2 cells tall**.

| `cells` (rows, cols) | size (mm) | typical use |
|---|---|---|
| (2, 2) | 60 × 60 | small square panel (donut, compact bar chart); 3 across a row |
| (3, 3) | 90 × 90 | square panel; 2 across a row |
| (4, 3) | 120 × 90 | tall bar chart (many categories) |
| (3, 6) / (2, 6) | full-width | wide single figure |

Rows must sum to **≤ 6 cells wide**: e.g. three 60 mm panels (2+2+2) or two 90 mm panels (3+3).
Fractional cells are allowed for off-grid sizes (e.g. `cells=(1.5, 1.5)` → 45 mm) but are an
explicit exception — prefer integer cells.

## The BasePlot pattern

Each chart is a `BasePlot` subclass. `BasePlot.__init__(ax=None, cells=(rows, cols))`:
- with `ax=None` it builds its own correctly-sized standalone figure (the normal case);
- with an `ax` supplied it draws into that axis (for composing, `cells` ignored).

`save(output_dir)` writes `png/<name>.png` (600 dpi) + `pdf/<name>.pdf` (vector, fonts embedded) and
closes the figure, targeting the plot's own figure explicitly.

```python
class MyBarPlot(BasePlot):
    def __init__(self, data, ax=None, cells=(4, 3)):
        BasePlot.__init__(self, ax=ax, cells=cells)
        self.name = "my_bar"              # output file stem
        self.is_available = len(data) > 0 # skipped (no-op save) if False
        self.ax.barh(data["label"], data["value"], color=BAR_DEFAULT)
        stylia.label(self.ax, xlabel="count", ylabel="", title="My bar")
```

An entry point builds the panels and saves each, recording footprints:

```python
def save_my_figures(data, output_dir):
    plots = [MyBarPlot(data, cells=(4, 3)), OtherPlot(data, cells=(2, 2))]
    footprints = {}
    for p in plots:
        if p.is_available:
            p.save(output_dir)
            footprints[p.name] = list(p.cells)
    json.dump(footprints, open(os.path.join(output_dir, "figure_cells.json"), "w"), indent=2)
```

**Colours** come from `src/plotting_colors.py` (semantic dicts anchored to stylia `ArticleColors`) —
never hardcode hex and never instantiate `ArticleColors`/`NamedColors` in a plot module. For a
one-off hue use the `hue(name, lighten=...)` accessor; for reference lines use `REFERENCE_LINE`
(silver); for structural box/median/marker lines use `INK`. Convention: turquoise = default/positive,
silver = reference/neutral lines, crimson = the "selected/kept" highlight; avoid plum/purple/grey as
data colours.

**Multi-panel in one file:** subclass `MultiPanelPlot` (in `plotting_base.py`) and call
`self._new_figure(nrows, ncols, cells, name)` to build a correctly-sized multi-axis figure; draw into
each axis. For a small-multiples grid where every cell is the same chart, subclass `GridPlot` and call
`self.build_grid(items, cols=..., name=..., panel_fn=..., color_fn=...)` — it owns the columns/rows
maths, edge-only axis labels and hiding of trailing empty cells.

## The shared plotting layer

A concrete panel should be **data-prep + declarative calls**, not hand-rolled matplotlib. Three shared
layers back every panel so all figures read as one system:

- **`plotting_colors.py`** — the single colour source (semantic dicts, `hue()`, `REFERENCE_LINE`, `INK`).
- **`plotting_base.py`** — `BasePlot` (+ helpers `self.label(...)`, `self.ref_line(v, axis)`,
  `self.legend({label: colour})`, `self._unavailable()`), plus the `MultiPanelPlot` / `GridPlot` bases.
- **`plotting_utils.py`** — ax-based primitives that carry the house style:
  `hbar`, `stacked_hbar`, `grouped_hbar` (bars); `box_with_jitter` (distribution boxes);
  `roc_panel` (one ROC cell); `heatmap` + `diverging_cmap`; `ref_line`; `swatch_legend` /
  `marker_legend`; `abbrev` / `abbrev_ticks` (genus abbreviation, e.g. `M. tuberculosis`).

**House style (applied uniformly):**
- **No panel titles.** Standalone panels carry no title (`self.label` accepts a `title` but ignores
  it) — titles and lettering are added in Illustrator. The *only* text kept is the per-cell
  identifier inside small-multiples (the model/organism name over each ROC or per-pathogen cell, and
  DR/SP subplot labels), set via `stylia.label`/`roc_panel` on the sub-axis.
- **Legends over data stay readable.** Every legend uses a semi-transparent white background
  (`LEGEND_KW` = `frameon=True, facecolor="white", framealpha=0.7, edgecolor="none"`), applied by
  `swatch_legend`/`marker_legend` and by any hand-built `ax.legend`.
- **No invisible colours.** Sequential/ordered palettes must stay visible on the white page — never
  let a shade fade to near-white/transparent (`_sequential` caps the lightest tint at `lighten=0.35`).
- All colour via `plotting_colors`; dashed silver reference lines (`ref_line`); horizontal bars
  ordered first-on-top (`hbar`); distribution boxes filled + tinted with `INK` median/whiskers and
  colour-matched jittered points (`box_with_jitter`); abbreviated genus tick labels everywhere.

When adding a panel, reach for an existing primitive first; only drop to raw matplotlib for genuinely
unique geometry (treemap, funnel, donut, twin-axis), and even then route colour/labels/legend through
the shared layer.

## Assembling a Nature figure

Panels are composed **in Illustrator**, not in code. Because `save_figure` uses
`bbox_inches="tight"`, the saved PDF is content-cropped and its page size deviates slightly from the
exact footprint — the authoritative footprint is in `figure_cells.json`, so place each panel to its
recorded cell box. Group panels into rows whose widths sum to ≤ 6 cells. Where two panels share a
legend, drop the per-panel legends and render one standalone legend panel to place once.

## Per-figure layouts

Script-specific footprints and page arrangements. Footprints are also written to each script's
`figure_cells.json`.

### 01_ersilia_metadata.py (`save_metadata_figures`, `src/plots_metadata.py`)

Six panels on the 3 cm grid: Tasks & subtasks `(3,3)`, Source Type `(2,3)`, Output `(3,3)`,
Biomedical Area `(4,3)`, Target Organism `(4,3)`, pathogen treemap `(3,3)`.

### 02_chembl_data_curation.py (`save_curation_figures`, `src/plots_chembl_curation.py`)

- **Row 1** — `curation_discard` + `chemspace_attrition` at **2×2 cells (60 mm)** each, per-panel
  legends omitted; `curation_outcome_legend` (a standalone 60 mm swatch key for the shared
  curation-outcome taxonomy) is placed once for both. `chemspace_attrition` hides any outcome that
  removes zero unique molecules across all pathogens (e.g. Co-ADD); it stays in the shared legend.
- **Row 2** — `pool_active_ratios` + `pool_cv_auroc` at **3×3 cells (90 mm)** each. DR and SP are
  combined per panel (DR = cobalt, SP = tangerine); both show **final pools** (step-25 grown ∪
  step-26 catch-all).
- **Coverage + funnel overlap unit** — `chembl_coverage` and `pipeline_funnel` are each **~45 mm
  (off-grid, footprint `[1.5, 1.5]`)** and both drawn in **crimson** (donut covered slice and funnel
  bars encode the same "kept" quantity), intended to be overlapped by hand into one slot.
- Remaining panels (`wholecell_sizes`, `binarisation_active_ratio`, `activity_ratio_flow`,
  `activity_ratio_per_pathogen`, `cutoff_sensitivity`, `pool_partition`, `merge_auroc`,
  `lowdata_auroc`) keep their default footprints, unassigned to a row yet.

### 03_chembl_models_performance.py (`save_performance_figures`, `src/plots_chembl_performance.py`)

Two panels per pathogen (15 pathogens, **196 models**), written to
`output/03_chembl_models_performance/individual_plots/{png,pdf}/` with their own
`figure_cells.json`. These 30 panels are **intermediate results, not paper figures** — they exist to
inspect every model, and the condensed cross-pathogen figures live at the top level of the output
dir. Note the grids show all of step 09, which includes 3 models step 10 discarded for AUROC < 0.7 —
they are drawn unmarked, so any caption claiming "the 193 hub models" would be wrong; the `retained`
column of the summary CSVs is the source of truth. See `scripts/README.md`.

- **`{pathogen}_roc_curves`** — small multiples, one ROC per model, **6 columns** (full 180 mm
  width) × `ceil(n/6)` rows, footprint `(rows, 6)`. Each panel sets `set_aspect("equal")` on a
  0–1 box so the ROC is square; without it `create_figure(nrows, ncols)` with no `width`/`height`
  squashes every panel into one default-height figure. Trailing empty cells are `axis("off")`.
  Axis labels only on the grid edges (TPR on column 0, FPR on the last occupied panel of each
  column). Curve colour = `auroc_shades(mean_auroc)`, cobalt fading colormap fitted to
  **[0.35, 1.0]** — the low anchor is below 0.5 on purpose, since fitting at chance level renders
  a chance-level curve white.
- **`{pathogen}_rank_boxplots`** — actives vs inactives out-of-fold rank, one box pair per model,
  footprint `(max(2, ceil(n/4)), 3)` — i.e. `MODELS_PER_CELL = 4`. Colours from
  `ACTIVE_INACTIVE_COLORS` (crimson / silver); legend sits **above** the axes because every row of
  the plot area is occupied by boxes.

## Gotchas

- **Log-axis bar charts:** draw bars from a *finite positive* baseline (`bottom=floor`,
  `height=value−floor`), not the default `bottom=0` — on a log axis 0 maps to −∞, so bars render
  as huge paths overflowing the axes in the vector PDF (invisible in the clipped PNG). See
  `PipelineFunnelPlot`.
- **Run in the stylia env:** figures must run in a conda env with stylia installed (e.g. `paper`),
  not `base`.
- **Small multiples need explicit sizing:** `stylia.create_figure(nrows, ncols)` sizes the *whole
  figure* to the format default, so an N×M grid gets N×M squashed panels. Always pass
  `width=cols/6, height=rows/6` and, for aspect-critical charts (ROC, PR, calibration), also
  `ax.set_aspect("equal", adjustable="box")`.
- **Data provenance:** record the source snapshot/version (e.g. ChEMBL release) in the script
  docstring; stage summaries via `scripts/00_download_data.py`, never ad-hoc copies.
