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
  `roc_panel` (one ROC cell); `pie_scatter` (points drawn as two-slice pies); `heatmap` +
  `diverging_cmap`; `ref_line`; `swatch_legend` /
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

Fourteen panels, sized for a **183 × 170 mm** page laid out as rows of three 60 mm panels
(3 across × ~2.8 rows), with the two pathogen panels forming a 90 mm row of two.

| panel | `cells` | mm |
|---|---|---|
| `task_subtask`, `task_subtask_waffle` | (2,2) | 60 × 60 |
| `source_type_by_subtask` | (1,2) | 60 × 30 |
| `output_by_subtask` | (1.15,2) | 60 × 34.5 |
| `license` | (2,2) | 60 × 60 |
| `runtime_100`, `image_size`, `license_class_pie`, `docker_architecture` | (1.5,1.5) | 45 × 45 |
| `biomedical_area`, `target_organism` (10 categories each) | (1.525,1.525) | 45.75 × 45.75 |
| `pathogen_circles`, `pathogen_voronoi`, `tag_cloud` | (3,3) | 90 × 90 |

The pathogen panels stay 90 mm square because each carries two legends and 15 genus labels.

**Quarter-width squares (`QUARTER_SQUARE`).** `biomedical_area` and `target_organism` are a quarter
of the page width square — 183/4 = 45.75 mm = **1.525 cells** — so the two sit side by side in half a
page row. Measured at that size: 34.06 mm of axes height over 10 bars = **3.41 mm row pitch** and
2.72 mm bars against a 2.12 mm tick label, i.e. ~1.3 mm clearance, so they compress without going
illegible. The long tick labels do eat into it: the drawn axes comes out only 21.8 / 22.8 mm wide out
of 45.7 mm, the rest being category names.

**Two-colour bars in those panels (`catchall_colors`).** Every *named* Biomedical Area and Target
Organism value is made up purely of Annotation models — Biomedical Area `Any` is 26 Annotation / 58
Representation / 19 Sampling and Target Organism `Any` is 40 / 57 / 19, while all 16 named areas and
all 65 named organisms are 100 % Annotation. **One apparent exception is a known upstream
mis-annotation, not real:** `eos93h2` (`image-mol-gpcr`) is recorded in Airtable as Task =
Representation / Subtask = Featurization with Target Organism = *Homo sapiens*, but it should be an
Annotation model; it is being corrected at source. It is drawn with the Annotation rows, which is
what the corrected data will say. Since the named bars therefore carry no task information of their
own, they take the
**Annotation hue (crimson)** rather than a palette of their own, so a reader who has learnt the task
colours from the task/subtask panels reads these two panels for free. `Any` takes **silver**, which is
the repo's reserved neutral for a catch-all bucket. Note neither panel carries a key for this, so a
caption must state what the two colours mean.

**When `eos93h2` is re-annotated upstream, change Task AND Subtask together.** The figures read task
two ways: `task_subtask`, the waffle and the `*_by_subtask` panels colour by Subtask through
`SUBTASK_PARENT`, while `runtime_100` / `image_size` and the Task counts read the `Task` column
directly. Setting `Task = Annotation` while leaving `Subtask = Featurization` would make those two
groups disagree **silently** — `task_subtask` would still draw the model as Representation-amber.
Moving it to a Subtask absent from `SUBTASK_PARENT` / `SUBTASK_ORDER` fails loudly instead
(`KeyError` in `SUBTASK_COLORS`), which is the safer failure. Expected deltas once fixed: Task
Annotation 131→132 and Representation 58→57; Subtask Featurization 52→51 with +1 on whichever
Annotation subtask it lands in.

**Bar thickness across the two subtask stacks.** Thickness is `0.8 × (axes height) / n_bars`, and
the axes is the figure minus a **fixed ~11.7 mm** tick+xlabel band that does *not* scale with figure
height. So it does not follow the footprint: for the 4-bar `output_by_subtask` to match the 3-bar
`source_type_by_subtask` (4.88 mm bars off an 18.31 mm axes) it would need
`(1 − 11.66/30) × 4/3 + 11.66/30` = **1.204 cells (36.1 mm)** — not the 40 mm a naive 4/3 ratio
suggests, which would overshoot to 5.7 mm. It is set to **1.15 (34.5 mm)** instead, a deliberate
layout compromise: taller than Source Type, shorter than a full thickness match, leaving its bars
**4.53 mm vs 4.88 mm (7% thinner)**. Raise the `_PANEL_OVERRIDES` entry to 1.204 to match exactly.
Both panels' *drawn axes* come out 41.1 vs 40.9 mm wide, so their plot areas align when stacked.
Fractional cells are the documented off-grid exception.

**Task / subtask colours.** `TASK_HUES` (in `plotting_colors.py`) is the single source: Annotation
= crimson, Representation = amber, Sampling = lime. `TASK_COLORS` takes the base hues and
`SUBTASK_COLORS` splits each into one shade per subtask (base hue for the largest, lightening to
`SUBTASK_LIGHTEN_FLOOR = 0.5`), so a task and its subtasks can never drift onto different hues.
Note these are the *same three hues* as `SOURCE_TYPE_COLORS` — safe only because the two are never
both colour-encoded in one panel: in the stacked Source Type panel, source type is encoded by bar
position and only the segments carry colour. `SOURCE_TYPE_COLORS` now survives solely for the two
pathogen panels' dot/cell colours.

**Stacked variants (`StackedFieldBarPlot`).** Source Type and Output are each drawn as one stacked
bar panel whose segments give the **subtask** breakdown of the field's own total, so a single panel
carries the joint distribution. Segments reuse the subtask palette rather than one of their own,
because the field is already encoded by the bar's position — a second colour dimension for it would
be redundant. Data comes from the `*_by_subtask_counts.csv` cross-tabs written by the script.

**Waffle (`TaskSubtaskWafflePlot`).** `task_subtask_waffle` draws one square per model in reading
order through the subtasks, so each subtask is a contiguous run and each task a contiguous band of
one hue. 16 columns divides the current 208 models into exactly 13 full rows; a different total
leaves `self.blank` trailing cells empty, which the script prints. It *shows* n instead of stating
it — every model is one mark — at the cost of precision: 52 vs 39 means counting.

This is the **one subtask panel that carries its own legend** (2 columns beneath the grid, labels
including each subtask's count — precisely what the waffle conveys badly), so it is self-contained
and does not have to travel next to `task_subtask`. The legend is also what squares the panel: the
grid alone crops to 1.23, and the legend band beneath it brings the whole panel to **1.05** on a
(2,2) footprint. Counter-intuitively, *reducing* the columns makes it worse, not better — a portrait
grid is narrower than the legend, which then fixes the content width while the extra rows add height.
Measured crop aspects: 13 cols 0.78, 14 cols 0.83, 15 cols 0.94, **16 cols 1.05**, 18 cols 1.22.
Only 13 and 16 columns divide 208 without a ragged last row. Legend spacing (`handlelength`,
`labelspacing`, `columnspacing`, `borderpad`) is tightened from the matplotlib defaults because the
band's height is the knob that tunes the panel's aspect.

**Licence (`license`).** Bars are the *simplified* licence identifiers — the script collapses
`-or-later` / `-only`, so GPL-3.0-or-later (51) and GPL-3.0-only (20) become one **GPL-3.0 bar (71),
which overtakes MIT (69)**. That drops a real legal distinction, so the ungrouped `license_counts.csv`
is kept as the record alongside `license_grouped_counts.csv`. The 27 models with no licence are
labelled **Not recorded**, never dropped — unknown terms are their own, worse, category for a reuser.
Colour is the reuse class (`LICENSE_CLASS` / `LICENSE_CLASS_COLORS`): turquoise Permissive (104),
periwinkle Copyleft (76), **fuchsia** Non-commercial (1, the sole CC-BY-NC-ND-4.0), silver Not
recorded. Fuchsia is used here and nowhere else in the repo: the convention deprioritises it because
it reads as emphasis, which is exactly what a 1-of-208 category needs — a hairline bar and a 1.7°
wedge in a well-behaved hue simply disappear — and emphasis is editorially right for the only licence
in the hub that forbids commercial reuse.
That split is a coarse classification, **not legal advice** — CC-BY-4.0 and CC0-1.0 are counted
permissive because they impose no share-alike duty, though neither is OSI-approved for software. The
4 AGPL-3.0 models are all the same upstream project (`molgrad`), whose repo does declare AGPL-3.0.
60 mm rather than quarter-width because `CC-BY-NC-ND-4.0` is a long tick label and the 4-entry key
needs the empty lower right of the axes.

**Tag cloud (`tag_cloud`) — the one raster panel in the repo.** `wordcloud` renders to a bitmap that
enters the axes via `imshow`, so this panel's PDF embeds an image rather than editable vector text:
the single documented exception to the vector-PDF rule. Generated at 1400 px square (600 dpi over
60 mm, more over 90) and seeded with `RANDOM_SEED` so layout is reproducible. Colour is an ordinal
periwinkle shade of the same count the font size encodes — redundant by design, to make the ranking
scannable without adding a claim.

*Size here is ordinal, never proportional.* Two settings buy legibility at the cost of fidelity:
frequencies go in as **sqrt(count)** (compressing the 49:1 range to 7:1) and `relative_scaling` stays
at 0.5 (blending frequency-proportional with rank-only sizing). Result: 5.47–26.79 pt across a 1–49
count range, so *Descriptor* (34) is 24.6 pt against *Antimicrobial activity* (49) at 26.8 pt — a
2 pt gap for a 15-model difference. On top of that a word's inked area depends on its character
count, so *Gram-negative bacteria* (13) occupies more of the panel than *Descriptor* (34). **Read
ranks only, and only coarsely; exact counts are in `tag_counts.csv`.**

Panel size is set by legibility, not layout. Measured with all 59 tags placed: raw counts at 60 mm put
**30 tags below 5 pt** with the ten single-model tags at **2.9 pt**; sqrt at 60 mm still leaves 11
below; sqrt at **90 mm** puts the smallest at **5.83 pt with none below**. `relative_scaling=1.0`
would restore font ~ sqrt(count) (hence area ~ count) but its 7:1 font range drops the tail to 3.6 pt
even at 90 mm. `TagCloudPlot` re-measures `min_font_pt`, `below_floor` and `n_placed` every run and the
script prints them, so shrinking the panel or adding tags cannot silently break legibility. **Do not
go below (3,3) without re-reading that output.**

**Technical panels (`TaskMetricBoxPlot`, `ArchitecturePiePlot`).** Three panels describing the
containers rather than the science: `runtime_100`, `image_size` and `docker_architecture`. The two box
panels are boxes-with-dots split by Task on a **log y axis** (both metrics span more than a decade —
runtime 16–1626 s, image size 291–10242 MB — and a linear axis crushes the bulk against the floor).
Per-task quartiles are written to `technical_metrics_summary.csv`.

*Footprint `SMALL_SQUARE = (1.5, 1.5)` = 45 × 45 mm, the smallest square this figure uses; four fit
across a page row.* Reaching it required **rotating the x tick labels 30°** and putting `n` inline
rather than on a second line: at 6 pt "Representation" alone is ~20 mm wide, so three horizontal
labels need the full 60 mm and already touched at that size. Even rotated, the axis furniture does
not shrink — the drawn axes is **31.9 × 28.3 mm inside a 45 mm panel**, so roughly 55 % of the panel
height is the chart and the rest is the y label, log ticks and rotated categories. That is the floor:
below 45 mm this stops being a small chart and becomes furniture with a chart in the corner.

*The `-1` sentinel is "never benchmarked", not zero.* Airtable stores `-1` where a runtime benchmark
was not run, so `TaskMetricBoxPlot` skips non-positive values and **never imputes them**. Any task
left with nothing keeps its slot, drawn empty and labelled "not measured" rather than omitted — a
missing box is information, a missing category is a misreading. Every tick label carries its own `n`,
and the script prints per-task coverage each run.

*Why the batch size is 100 (`RUNTIME_BATCH` in `default.py`).* Coverage collapses as the batch grows:

| batch | Annotation | Representation | Sampling |
|---|---|---|---|
| 100 | 129/131 | 57/58 | **10/19** |
| 1,000 | 123/131 | 55/58 | **0/19** |
| 10,000 | 95/131 | 47/58 | **0/19** |

100 is chosen so generative models appear at all. **The trade-off is what the number then means:** at
this batch size the median CP3/CP1 ratio is **0.86** — running 100 molecules takes no longer than
running 1 — so for annotation and representation models the value is dominated by container startup,
not per-molecule work. It is the wall-clock a user waits for, **not a throughput measure**; do not
divide it by 100. (At 1,000 molecules the ratio rises to 3.29, which *is* an amortisation result, but
costs every Sampling model.)

Medians at 100 molecules: Annotation **34.1 s** (IQR 22.9–61.8), Representation **28.6 s**
(IQR 21.5–50.7), Sampling **481.4 s** (IQR 353.7–524.1) — generative models are ~14× slower, which is
the finding this batch size buys. Note Sampling rests on **10 of 19** models, so treat that box as
indicative. Image size medians: Annotation **4.5 GB**, Representation **1.8 GB**, Sampling **613 MB**.

**Pies (`_pie`).** Both pie panels share one helper: slices clockwise from 12 o'clock, labelled
*outside* with count and share so no legend is needed. A share under 0.5 % is printed as `<1%`, never
as the `0%` that `.0%` formatting would give — the single non-commercial licence is 1/208 and must not
be labelled as none. In-place labelling only works while every slice can own a label, so a pie built
through this helper needs a handful of comparable slices, not a long tail; push detail to a bar panel.

**`license_class_pie` shows the four reuse classes, not the ten licences** — Permissive 104 (50 %),
Copyleft 76 (37 %), Not recorded 27 (13 %), Non-commercial 1 (<1 %), colours matching the `license`
bar panel. That is a hard limit, not a shortcut: four of the ten licences cover exactly one model
each, which on 208 models is a **1.7° wedge**, invisible and impossible to label or tell apart from
the other three. Per-licence detail is the bar panel's job. Slices run in descending order so the
Non-commercial hairline sits next to the 12 o'clock start rather than between two large slices, and
even so its label is the only part of it a reader will see.

`docker_architecture` is a two-slice pie labelled outside with count and percentage, so it needs no
legend: **AMD64 + ARM64 129 (62%)** in turquoise vs **AMD64 only 79 (38%)** in periwinkle. Two
substantive hues, not silver — both are real build targets, and silver would have cast x86-only as a
residual bucket; periwinkle rather than a warning hue so the panel does not read as flagging those
builds as a problem.
The field only ever holds `AMD64` or `AMD64,ARM64` — there is no ARM-only build, so this is "also built
for ARM" versus "x86 only". It is a **snapshot, not a trend**: dual-arch is 45% among models
incorporated in 2021 and 77% among 2026 ones, so the 62% is accumulated stock rather than current
build practice, and a caption must give the metadata snapshot date.

**Legends across this figure's subtask panels.** `task_subtask` uses `SUBTASK_COLORS` (not
the parent-task base hue), so every subtask appears there as a labelled bar in its own colour —
that panel *is* the subtask key, and it therefore carries no legend of its own. The two `_by_subtask`
stacks carry none either: at 60 × 30 mm a 6-entry key would be as tall as the plot, and they are meant
to sit beside `task_subtask`. The **waffle** has its own (2 columns beneath the grid), as does the
`license` bar (4 reuse classes, lower right). **Consequence for layout: if you use a `_by_subtask`
stack, `task_subtask` or `task_subtask_waffle` must travel with it, or render a standalone subtask
key.**

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

**Top level of the output dir** holds the condensed cross-pathogen figures: `pathogen_dataset_sizes`
(3,6) and `pathogen_consensus_auroc` (4,3). Per-pathogen panels are intermediate results and live in
`individual_plots/`.

Two panels per pathogen (15 pathogens, **196 models**), written to
`output/03_chembl_models_performance/individual_plots/{png,pdf}/` with their own
`figure_cells.json`. These 30 panels are **intermediate results, not paper figures** — they exist to
inspect every model, and the condensed cross-pathogen figures live at the top level of the output
dir. Note the grids show all of step 09, which includes 3 models step 10 discarded for AUROC < 0.7 —
they are drawn unmarked, so any caption claiming "the 193 hub models" would be wrong; the `retained`
column of the summary CSVs is the source of truth. See `scripts/README.md`.

**`pathogen_dataset_sizes`** — footprint `(3,6)`, **full page width (180 × 90 mm)**. One box per
pathogen over the sizes of its modelled datasets on a **log y axis** (range 96 → 334,766 compounds),
with one small **pie per dataset** drawn on top whose filled share is that dataset's active fraction.
Box = the pathogen's size distribution, circle = one dataset's balance, in one panel. Tick labels
carry the dataset count and are genus-abbreviated + italic; ordered by dataset count descending.
Data: `dataset_sizes.csv`, derived from `10_reports/10_reports.csv` (the 193 step-10 keeps).

*Added negatives are excluded from both the size and the ratio.* `n_compounds` counts them, so size
is `n_compounds − n_added_negatives − n_added_decoys` and the active fraction is `n_positives / size`.
That reproduces the curation pipeline's own `n_mol_after` / `ar_after` **exactly** — verified for all
151 models that join to `25_pool_summary`. Only 54 of 193 models were given negatives, so the medians
barely move (1115 → 1098 compounds; active ratio unchanged at 0.365), but for those 54 the change is
large: `mtuberculosis/DR_0012` goes from 2450 compounds at 0.50 active to **1411 at 0.87**. Decoys are
zero for every model, but are subtracted anyway so the definition does not silently depend on that.

*Fixed pie diameter.* The y position already encodes size, so scaling the circles would double-encode
it, and across a 3,500× range the largest would swamp the panel while the smallest vanished. The full
page width is a **legibility requirement, not a preference**: 15 pathogens share the x axis and
*P. falciparum* alone holds 51 datasets, so columns need ~12 mm each for the pies to separate under
horizontal jitter (`JITTER = 0.30` of the column pitch, seeded by `RANDOM_SEED`). The box is filled at
`lighten=0.22` because the pies sit inside it and lose contrast against full-strength cobalt.

**Read a pie as a gestalt, not a measurement.** At ~2.5 mm "mostly active" vs "mostly inactive" is
clear; 45 % vs 50 % is not. Exact ratios are in `dataset_sizes.csv`. Note **54 of 193 datasets are
majority-active** once added negatives come out, and none are 0 % active. `pathogen_activity_ratios`
below exists precisely because this panel cannot resolve that.

**`pathogen_activity_ratios`** — footprint `(3,6)`, the **alternative** to the panel above on the same
`dataset_sizes.csv`: it swaps the encodings, putting the active fraction on the y axis (0–1, linear)
where a reader can actually resolve it, and demoting size to **dot area**. A dashed silver line marks
the 0.5 balance point. Both panels are rendered; pick one at layout time.

*Dot area is affine in √size, not proportional to size.* `SIZE_REF_AREA = 30` pt² at
`SIZE_REF = 10,000` compounds, everything else scaling as `√(size / SIZE_REF)`. Area-proportional dots
over a 3,500× range would mean a 59× radius ratio — the largest would blot out its column and the
smallest vanish. This is the same compromise `AREA_EXPONENT` makes in the script-01 pathogen treemap.
The size key (100 / 1,000 / 10,000 / 100,000) is generated by the same function, so read sizes off it
rather than comparing areas by eye.

*Mean bars.* A short horizontal bar per column marks the pathogen's **unweighted mean** active
fraction — every dataset counts once regardless of size. Because the very large datasets sit near 0 %
active, a size-weighted mean would fall well below the bar; the bar answers "what does a typical
dataset for this pathogen look like", not "what fraction of its compounds are active". Drawn in the
full-strength hue against the pale dots.

*Y headroom.* `Y_PAD = 0.07` above 1 and below 0, with ticks stopping at 0 and 1 so the padding never
implies out-of-range values. Without it the largest dots (~2.5 mm radius) sit on the spines.

*Palette.* Dots take the pale `FILL_LEVELS = (0.62, 0.38)` palette so heavy overlap stays readable,
outlined and topped by the darker `ACCENT_LEVELS = (None, 0.55)` of the same hue. **Both levels are
needed in each call**: `distinct_colors` uses its second pass to get past the 9-hue limit, so a single
lighten value would make pathogen 10 identical to pathogen 1. The two calls index hues identically, so
a dot and its mean bar read as one category.

*Size key (`plotting_utils.nested_size_legend`).* Nested circles sharing a bottom tangent, drawn as
`Ellipse` patches with the point radius converted **separately per axis** so they stay round whatever
the aspect or scale. It must be called last — it reads the live transform, so the axis limits have to
be final. Label rows are **de-collided and reached with leaders**: nesting packs the ring tops within
`2 × rmax` of each other, about 3.4 mm for this decade key, while four labels need ~8 mm of line
height, so anchoring each label at its own ring's top guarantees overlap.

*Colour is redundant with the x axis* — one hue per pathogen via `plotting_colors.distinct_colors`,
while the tick already names it. Deliberate: it makes a column scannable and lets a reader follow a
pathogen across panels, and it means the palette need not supply 15 unambiguous hues. **Only 9
substantive hues exist in ArticleColors**, so `distinct_colors` cycles back through them at
`lighten=0.55`, placing a hue and its tint 9 apart. There is no colour legend — 15 swatches would not
be readable at panel size, and the labelled axis is the real key. Do not use `distinct_colors` as the
*sole* encoding for more than 9 categories.

This panel surfaces what the pie version hides: the **largest datasets sit at near-zero active
fraction** while the small ones are balanced or active-heavy.

*Implementation (`plotting_utils.pie_scatter`).* Pies are matplotlib **marker paths**, not `Wedge`
patches — that is what keeps every pie round and the same physical size on an axis with a log scale, a
categorical x and no equal aspect, where a patch sized in data units would come out elliptical and
varying. Matplotlib scales a custom marker path by a flat 0.5 instead of normalising it to its
bounding box, so a 10 % wedge and a 90 % wedge share one radius (verified: `max|coord| == 0.5` at
every fraction) — without that, slice size would leak into apparent circle size. The "inactive" slice
is one full-circle scatter call for all points with the active wedges drawn over it, so N pies cost
N+1 draw calls, not 2N.

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
