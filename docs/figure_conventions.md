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
never hardcode hex. Convention: turquoise = default/positive, silver = reference/neutral lines,
crimson = the "selected/kept" highlight; avoid plum/purple/grey as data colours.

**Multi-panel in one file:** a figure may hold >1 axis (e.g. small multiples) — create the figure
with `stylia.create_figure(nrows, ncols, width=cols/6, height=rows/6)`, store it as `self.fig`, and
draw into each axis; `save()` only needs `self.fig`/`self.name`/`self.cells`/`self.is_available`.

## Assembling a Nature figure

Panels are composed **in Illustrator**, not in code. Because `save_figure` uses
`bbox_inches="tight"`, the saved PDF is content-cropped and its page size deviates slightly from the
exact footprint — the authoritative footprint is in `figure_cells.json`, so place each panel to its
recorded cell box. Group panels into rows whose widths sum to ≤ 6 cells. Where two panels share a
legend, drop the per-panel legends and render one standalone legend panel to place once.

## Gotchas

- **Log-axis bar charts:** draw bars from a *finite positive* baseline (`bottom=floor`,
  `height=value−floor`), not the default `bottom=0` — on a log axis 0 maps to −∞, so bars render
  as huge paths overflowing the axes in the vector PDF (invisible in the clipped PNG). See
  `PipelineFunnelPlot`.
- **Run in the stylia env:** figures must run in a conda env with stylia installed (e.g. `paper`),
  not `base`.
- **Data provenance:** record the source snapshot/version (e.g. ChEMBL release) in the script
  docstring; stage summaries via `scripts/00_download_data.py`, never ad-hoc copies.
