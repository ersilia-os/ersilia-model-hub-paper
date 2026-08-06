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
(silver); `INK` is for structural marks that must read on an opaque body (the median of a *filled*
box — see the box style below). Convention: turquoise = default/positive,
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
  `hbar` (+ `inside_labels=True` / `place_inside_labels`, names on the bars and no y axis, for panels
  too narrow to afford a tick-label gutter), `stacked_hbar`, `grouped_hbar` (bars);
  `box_with_jitter` (distribution boxes);
  `roc_panel` (one ROC cell); `pie_scatter` (points drawn as two-slice pies); `heatmap` +
  `diverging_cmap`; `ref_line`; `swatch_legend` /
  `marker_legend`; `abbrev` / `abbrev_ticks` (genus abbreviation, e.g. `M. tuberculosis`).

**House style (applied uniformly):**
- **No panel titles.** Standalone panels carry no title (`self.label` accepts a `title` but ignores
  it) — titles and lettering are added in Illustrator. The *only* text kept is the per-cell
  identifier inside small-multiples (the model/organism name over each ROC or per-pathogen cell, and
  DR/SP subplot labels), set via `stylia.label`/`roc_panel` on the sub-axis.
- **Axis labels are sentence case**, always — "Own-assay AUROC", not "own-assay AUROC". You do not have
  to remember this: `plotting_base` wraps `stylia.label` once at import with
  `plotting_utils.sentence_case`, so a panel may write its label in any case and the figure comes out
  right. Only the first character is touched, and only when the first *word* has no internal capital —
  that is what keeps `pIC50 (nM)`, `cLogP`, `-log10(q)` and `|Spearman rho|` intact, and what makes this
  a helper rather than a `.capitalize()` call (which would also lowercase `AUROC`). Titles are *not*
  transformed: standalone panels have none, and small-multiple cell identifiers are genus abbreviations
  already cased correctly. A panel that bypasses stylia for `ax.set_xlabel` must call `sentence_case`
  itself — three do (`euos_overlap`'s inline twin-axis labels, the two `plots_column_jaccard` charts).
- **Legends over data stay readable.** Every legend uses a semi-transparent white background
  (`LEGEND_KW` = `frameon=True, facecolor="white", framealpha=0.7, edgecolor="none"`), applied by
  `swatch_legend`/`marker_legend` and by any hand-built `ax.legend`.
- **No invisible colours.** Sequential/ordered palettes must stay visible on the white page — never
  let a shade fade to near-white/transparent (`_sequential` caps the lightest tint at `lighten=0.35`).
- **Distribution boxes defer to their swarm.** `box_with_jitter` / `box_from_stats` draw an
  **unfilled** box outlined in the category colour at `BOX_LINEWIDTH = 0.5` — the house line weight,
  matching the axes spines. matplotlib's boxplot default is `1.0`, i.e. *twice* stylia's
  `lines.linewidth`, which made every box the heaviest mark in its panel and buried the jittered
  points it sits over. The median is the same colour drawn heavier (`MEDIAN_LINEWIDTH = 1.4`) so it
  still reads as the summary statistic. **Pass `face=` only when the box has something to hide or
  nothing to reveal:** a body is right where marks are drawn *inside* it (`PathogenDatasetSizesPlot`'s
  pies) or where there is no swarm at all and the fill carries the category (`rank_boxplots`). On a
  filled box the median reverts to `INK`, the one colour that survives an opaque body — several
  panels pass `face=color`, where a same-hue median would vanish.
- All colour via `plotting_colors`; dashed silver reference lines (`ref_line`); horizontal bars
  ordered first-on-top (`hbar`); abbreviated genus tick labels everywhere.

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

Thirteen panels, sized for a **183 × 170 mm** page laid out as rows of three 60 mm panels
(3 across × ~2.8 rows), with `pathogen_circles` taking a 63 mm block of its own.

| panel | `cells` | mm |
|---|---|---|
| `task_subtask` | (2,2) | 60 × 60 |
| `task_subtask_waffle` | (1.5,1.5) | 45 × 45 |
| `output_by_subtask` | (0.803,1.691) | 50.7 × 24.1 |
| `source_type_by_subtask` | (0.778,1.691) | 50.7 × 23.3 |

*(both are figure sizes; the saved pages are ~52 × 25 mm — see the stack note below)*
| `runtime_100` | (1,1.595) | 47.9 × 30 |
| `image_size`, `output_dimension` | (1,1.122) | 33.7 × 30 |
| `license_class_donut`, `biomedical_area_donut` | (1.167,0.792) | **24.9** × 35 |
| `docker_architecture` (2 legend rows, so shorter) | (1.0,0.792) | **24.9** × 30.5 |
| `biomedical_area` (4 groups, names on the bars) | (1,0.797) | **25** × 30 |
| `target_organism` (10 categories) | (1.525,1.525) | 45.75 × 45.75 |
| `pathogen_circles` | (2.5,2.5) | crop **63 × 76** |

**`pathogen_circles` comes in under a 65 mm width ceiling** at `cells=(2.5, 2.5)`, cropping to
**63.05 × 76.20 mm** — a packed layout is mostly whitespace between circles, so it absorbs the
reduction a space-filling layout could not.

**Its crop width is not monotonic in `cells`, so never infer it — measure it.** The width is set by
where the iterative label placement lands, not by the packing, and measured crops go 2.3 → 65.6 mm,
2.4 → 67.7, **2.5 → 63.1**, 2.55 → 66.2, 2.6 → 67.4. A *smaller* footprint can give a *wider* panel.
Everything in it is sized in points and does not scale with the footprint, so re-check the labels
before going lower.

**Its genus labels already collide at *both* sizes.** `H. pylori`, `S. aureus`, `Campylobacter` and
`Enterobacter` sit on top of circles rather than clear of them at 90 mm as well as at 75 mm, so this is
a pre-existing limit of the radial push-out, not something the reduction introduced. The white stroke
behind each label keeps them readable; fixing it properly means a real label-placement pass.

**`biomedical_area` is a 25 mm strip with the labels on the bars** (`_narrow_strip_cells`,
`hbar(inside_labels=True)`) — the narrowest panel in the repo, `create_figure(width=0.1325)`, cropping
to **25.23 mm**. At that width a tick-label gutter is impossible: the gutter is a **fixed ~14 mm**
whatever the panel's width, which would leave ~4 mm of bar and turn the panel into a column of text
with a hairline chart attached. So the y axis goes away entirely and each category name is drawn in
front of its own bar, in the **default black** every other label uses.

*Height is shared with the technical box row* (`_BOX_ROW_HEIGHT_CELLS`, 30 mm declared → **31.24 mm**
cropped against the row's 31.17–31.27), so the strip and the three box panels form one 145 mm band.
Height therefore does **not** follow the bar count — the four groups only decide how the fixed
18.19 mm axes is divided, giving a 4.55 mm pitch. `_STRIP_BAR_FRACTION = 0.70` rather than
matplotlib's 0.8 because at that pitch 0.8 gives a 3.6 mm bar with a 0.9 mm gap, a near-solid block;
0.70 leaves a 3.2 mm bar (close to the figure's 2.9 mm mark weight, and enough to seat a 5 pt label)
with a 1.4 mm gap.

**Label placement is adaptive, and that is the whole trick** (`plotting_utils.place_inside_labels`).
Per bar: if the label fits inside, it is **centred on the bar** and sits on a semi-transparent white
chip (`LABEL_CHIP_KW`, the same device `LEGEND_KW` uses over data); if it does not fit, it is set
**immediately after the bar's end** with the chip removed, since over the page a white chip is just a
box. Bars can therefore stay **full-strength crimson** — legibility is the chip's job, not the bar's —
instead of being washed out to a pale tint to keep dark text readable.

Two things were tried and rejected. **Unconditional centring destroys the mark it labels:** *Antiviral*
is 6.3 mm of text on a 2.4 mm bar, and a chip centred there covers the bar completely, leaving a white
pill with a red rim where the data should be. **A stroked outline** (`patheffects.withStroke`, as the
pathogen treemaps use) fails at this type size: a stroke thick enough to separate black text from
saturated crimson floods the counters of O, e and A, and the glyphs stop reading as type. It works on
the treemaps because their labels sit well inside large filled cells.

Placement must run **after `tight_layout`** — it weighs a label width fixed in *points* against an axes
width layout is still free to change — so `FieldBarPlot.measure_labels` does it just before saving. In-bar
labels are also **not clipped**, so one overrunning a spine would spill outside the plot area and *set*
the crop width, breaking the fixed 25 mm. The script prints the placement split and the widest label
against the axes every run (currently 2 centred inside, 2 after their bars; 17.74 mm axes vs 10.92 mm
for *Antimicrobial*). The count-axis label is shortened to **"Models"** for the same reason —
"Number of models" is 17.6 mm and would overhang.

**Four groups, over Activity prediction models only** (`BIOAREA_GROUP` / `ACTIVITY_SUBTASK` in
`default.py`, signed off 2026-08-02): **Antimicrobial 50, ADMET 29, Antiviral 7, Other 8** of
**92** models. `default.py` holds the classification and the evidence for every membership; the
per-area breakdown is written to `biomedical_area_groups.csv`. Substantive groups run by size
descending with the catch-all pinned **last** however large it grows, so `Other` (8) sits below
`Antiviral` (7) — the residual is always the bottom bar, even when it is not the shortest.

*Restricting to Activity prediction is what keeps `Other` honest.* `Any` (no area declared) is **22 of
the 39** "Property calculation or prediction" models — generic property predictors like logP and
solubility have no area to declare — but only **4 of the 92** Activity prediction ones. So `Any` folds
into `Other` without turning it into a catch-all, and the four bars account for *every* model in the
subtask: nothing dropped, nothing unexplained. Across all of Annotation `Other` would have been 30 and
outranked Antiviral 4:1.
**The cost is ADMET: 44 → 29**, because ADMET splits 29 Activity / 15 Property across the two
Annotation subtasks and that split is not semantically clean. ADMET therefore falls behind
Antimicrobial here, and the bar means "ADMET models annotated as Activity prediction", not "ADMET
models". The two antimicrobial models excluded are genuinely not activity endpoints — `eos4n4d`
(Gram-negative *accumulation*, a permeability property) and `eos2xeq` (antibiotic *downselection*, a
novelty/synthesizability filter).

**Biomedical Area is multi-value, so counts are of distinct MODELS, not area assignments.** Grouping
absorbs most of the multiplicity: 10 of the 12 multi-area Annotation models have all their areas inside
one group (AMR+Pneumonia, AMR+Diarrhoea, Gonorrhea+AMR). **Two models still span two groups and are
counted in both** — `eos2zmb` (Cancer + AIDS) and `eos7kpb` (ADMET + Malaria + Tuberculosis) — so the
four bars **sum to 94 against 92 models**. That is the metadata's own claim, left unresolved rather
than forced into one bucket; the script prints it every run and a caption should note it.

**Two deliberate stretches in "Antimicrobial".** *Malaria* (8) is *Plasmodium falciparum*, a
**protozoan**, inside "antimicrobial" only on the broad clinical definition; *Schistosomiasis* (2) is
*Schistosoma mansoni*, a multicellular **helminth**, not a microorganism at all. Both are kept because
Ersilia's own naming already treats them that way — the S. mansoni model's slug is literally
`antimicrobial-activity-smansoni`. **A caption saying "antimicrobial" therefore covers antibacterial,
antifungal, antiprotozoal *and* antihelminthic activity.** Rename the group to "Anti-infective" if that
overclaims for a given venue. Memberships were checked against each model's Target Organism rather than
inferred from the area name: *Peptic ulcer disease* is `eos9eyo` / *H. pylori* (bacterial, not the NSAID
aetiology), *Diarrheal diseases* is all *Campylobacter* / *E. coli*, and *Candidiasis* / *Mycetoma* are
fungal.

**An unmapped area fails loudly.** The counting step raises `KeyError` on any Biomedical Area value
missing from `BIOAREA_GROUP`, so a new area added upstream cannot silently vanish from the figure.

**Bar colour.** Every bar here is Annotation / Activity prediction, so colour distinguishes nothing
*within* the panel; the crimson survives purely as a cross-panel cue to the Annotation task.
`catchall_colors` is therefore no longer used by this panel (`target_organism` still uses it).
**`target_organism` is untouched** — 45.75 mm, tick labels, ungrouped, `Any` included — so the two
panels no longer match and cannot be presented as a pair without saying so. It also shows exactly why
`Any` had to be handled: at 116 models it sets the count axis and leaves every named organism a
hairline. Every *named* Biomedical Area and Target Organism value is made up purely of Annotation
models. On the **Jul 30 2026 snapshot** Biomedical Area `Any` was 26 Annotation / 58 Representation /
19 Sampling and Target Organism `Any` was 40 / 57 / 19, with all 16 named areas and all 65 named
organisms 100 % Annotation.

**RESOLVED 2026-08-06 — `eos93h2` is no longer an exception.** It was the one model breaking the rule:
recorded in Airtable as Task = Representation / Subtask = Featurization while carrying Target Organism =
*Homo sapiens*, when it is really a GPCR activity predictor. The upstream correction has landed, and it
changed **both** fields as this section required:

```
Task:     Representation  ->  Annotation
Subtask:  Featurization   ->  Activity prediction
```

Nothing else in the record moved (Status still Ready, Target Organism still *Homo sapiens*, Output
Dimension still 10). Script 01 already drew it with the Annotation rows, so **no panel changes** —
which is what this section predicted. Re-verified on the fresh metadata: **90 Ready models carry a
named Target Organism and 0 are non-Annotation; 103 carry a named Biomedical Area and 0 are
non-Annotation.** The invariant now holds with no exceptions, so the named-value bars taking the
Annotation hue is unconditionally correct rather than correct-in-advance.

`Any` takes **silver** in `target_organism`, the repo's reserved neutral for a catch-all bucket; that
panel carries no key, so a caption must state what its two colours mean. On the fresh snapshot
Target Organism `Any` is 40 Annotation / 61 Representation / 23 Sampling, Biomedical Area `Any` is
27 / 61 / 23.

**The Task/Subtask coupling still matters for the next such fix.** The figures read task two ways:
`task_subtask`, the waffle and the `*_by_subtask` panels colour by Subtask through `SUBTASK_PARENT`,
while `runtime_100` / `image_size` and the Task counts read the `Task` column directly. Setting
`Task = Annotation` while leaving `Subtask = Featurization` would make those two groups disagree
**silently** — `task_subtask` would still draw the model as Representation-amber. Moving a model to a
Subtask absent from `SUBTASK_PARENT` / `SUBTASK_ORDER` fails loudly instead (`KeyError` in
`SUBTASK_COLORS`), which is the safer failure.

**Do not read the observed count deltas as this fix's effect.** This section previously predicted Task
Annotation 131→132, Representation 58→57 and Subtask Featurization 52→51 for the fix *in isolation*,
and that arithmetic was right. The measured Ready-only movement between the two snapshots is different
because six models were added and four changed Status in the same pull: Annotation **131 → 130**
(+1 from `eos93h2`, −1 each for `eos18ie` and `eos1lb5` going to *In maintenance*), Representation
**58 → 61** (−1 `eos93h2`, +3 new featurizers, +1 new projector), Sampling **19 → 23**, total Ready
**208 → 214**. By subtask: Activity prediction 92 → 91, Featurization 52 → 54, Generation 8 → 12,
Projection 6 → 7, Property calculation 39 → 39, Similarity search 11 → 11.

#### The two `*_by_subtask` stacks are sized to a page budget (2026-08-05)

**They fill 52 × 50 mm between them, and their bars are the same thickness.** Those two
requirements are the whole specification, and everything below follows from them. `Output` (4 bars)
goes on **top**, `Source Type` (3 bars) **below** it carrying the `"Number of models"` axis title, so
the unit states its quantity once, under the block, where a reader expects the shared unit of a
stacked pair. Delivered: pages **51.93 × 25.36** and **51.96 × 24.60 mm** (pair 51.96 × 49.96) with
bars of **3.028 and 3.029 mm** — equal to 0.0 %, both dimensions inside the budget. The script prints
all of it every run. Bar thickness therefore sits just above the **2.9 mm** mark weight the technical
box row was built around, where the earlier 44.5 mm budget put it just below at 2.42 mm.

**The budget is in saved-page millimetres, not `cells`.** `_SUBTASK_STACK_WIDTH_MM = 52.0` (per
panel) and `_SUBTASK_STACK_HEIGHT_MM = 50.0` (the two summed) are what the panels *measure on the
page* — resized 2026-08-05 from 47.75 × 44.5, by editing those two numbers and nothing else; `_stack_cells` works backwards from them to footprints. Check the result with
`plotting_base.pdf_page_mm`, which reads the PDF's own `/MediaBox` — **not** `crop_size_mm`, whose
width runs ~1.2 mm under what `savefig` writes and which would silently put a calibrated panel over
budget.

**The size chain has three links, and only two of them are constant.** Getting this wrong is what
made two earlier attempts miss:

```
footprint (cells)  --quantised-->  figure  --band-->  axes  --pad-->  saved page
                    canvas pixels          CONSTANT          CONSTANT
```

- `figure height − axes height` is the **band**: `_AXES_BAND_MM = {True: 11.6634, False: 8.2659}`,
  everything `tight_layout` spends on x tick labels plus the optional axis title. Constant to four
  decimals across 19–23 mm and both fields. Dropping the title is worth **3.40 mm**.
- `page − figure` is savefig's tight-bbox **pad**: `_PAGE_PAD_MM = 1.271` in height, `1.27` in width.
- **The footprint is NOT in a fixed relationship to either.** matplotlib's canvas is a whole number
  of pixels, so the figure is the footprint *floored onto a pixel grid*. A band measured against the
  declared footprint therefore drifts, and an earlier version of this note recorded a bogus
  "axes height is linear in declared height with slope 1.016" — that was **aliasing**, from sweeping
  the height in 3 mm steps, which is 11.81 quanta apart. Sample finer than the quantum or you will
  fit a slope to a staircase.

**The layout dpi is raised to 600 for these panels, and that is a sizing decision, not a resolution
one** (`_LAYOUT_DPI` in `plots_metadata.py`, applied through `rc_context` at figure creation — set
afterwards it does not re-quantise). The pixel grid *is* the set of sizes the layout can express:
**0.254 mm at matplotlib's default dpi 100**. Equal bar thickness needs a specific, non-round split of
the height between the two panels, and on a 0.254 mm grid the nearest achievable split leaves the bars
**1.8 % apart** — a floor no choice of constants can beat. At dpi 600 the grid is 0.042 mm and the same
solve lands inside **0.2 %**. It is set on `StackedFieldBarPlot` alone: a global dpi change would move
every measured panel size in the repo by up to a quarter of a millimetre, and this figure has a dozen
panels calibrated to tenths.

**Equal thickness is a ratio of y SPANS, not of bar counts.** Thickness is
`0.8 × axes_height / _y_span(n)`, so equal thickness means equal `axes_height / _y_span(n)` and the
axes height is divided between the panels in the ratio of their spans. The span is **not** `n`:
`stacked_hbar` leaves the y axis autoscaled, so the view covers `n − 1` between the first and last bar
centre, plus half a bar overhanging at each end, all inflated by matplotlib's 5 % margin top and bottom
— `(n − 1 + 0.8) × 1.1`, i.e. **3.08 units for 3 bars, not 3**. Splitting by bar count instead is a
3 % error in the very quantity being equalised. Because only the lower panel carries the axis title it
spends 3.40 mm more on its band. At 44.5 mm that happened to make the two *footprints* nearly equal;
at 50 mm the band no longer cancels the 4 : 3 axes split and the upper panel is the taller footprint
(24.09 vs 23.32 mm). Neither is a target — the footprint is an output of the solve.

**Every solved size is backed off by half a canvas pixel** (`_PIXEL_MM`). The figure size is the
requested footprint rounded to the **nearest** pixel, not floored, so a panel can come out up to half a
pixel *larger* than asked — which put the width 0.02 mm over the ceiling before the back-off. Nothing
on the page, but a budget is a ceiling and someone else's layout is built on the number.

**Each stack still keeps its own count axis.** They are **not** on a shared scale — each autoscales to
its own longest bar (Source Type reaches 156, Output 106) — so **do not compare a bar in one against a
bar in the other by length**. `show_xlabel=False` on the upper panel drops only the axis *title*; its
tick labels stay, so it remains readable on its own rather than depending on its partner for a scale.

**Known residual: the axes frames are ~0.35 mm out of register.** Measured off the saved files, the
axes run **18.73 → 49.39 mm** on `output_by_subtask` and **18.37 → 49.47** on `source_type_by_subtask`
(measured from the page's own left edge) — left edges 0.37 mm apart, right edges 0.07 mm, drawn axes
widths 30.66 vs 31.10 mm. The left-edge error is **fixed in mm, not proportional**: it is a difference
of tick-label widths, so it does not shrink when the panels grow. `tight_layout`
sizes each panel's y tick-label column to its own content and `Compound` is wider than `Replicated`.
So **align the pair on their axes frames, not their page boxes**, and treat the 0.7 mm axes-width
difference as one more reason not to read bar lengths across the pair. Making the *axes* widths equal
instead is possible — give each panel a declared width compensating for its own label column, at the
cost of slightly different page widths — but it needs two more calibration constants tied to specific
tick-label strings, which break when a category is renamed upstream. Fixing it properly means drawing
both into one figure with hand-placed axes; **`01b_community_stats.py` is the worked example**, four
tracks with left/right axis edges identical to six decimal places. Use it as the template whenever
panels must be read against a common axis — this pair need not be, since their scales differ anyway.

**The waffle no longer travels with them as a flush block.** `task_subtask_waffle` is a
**quarter-width square** (`create_figure(width=0.25)` = **1.5 cells = 45.0 mm**, a quarter of stylia's
180.09 mm canvas, so 0.75 mm narrower than the `QUARTER_SQUARE = 1.525` used by `biomedical_area` /
`target_organism`, which is a quarter of the 183 mm *page*; both constants exist in
`plots_metadata.py` as `QUARTER_WIDTH` / `QUARTER_SQUARE` — don't unify them without deciding which
reference you want). The stacks are now sized to their own 44.5 mm budget rather than to the waffle's
46.09 mm crop, so at the 50 mm budget the column runs **3.9 mm taller** than the waffle and the two blocks are
aligned by hand at layout time. The waffle is still the pair's **subtask colour key** (see the legend note
below), so one of it or `task_subtask` must travel with them.

**Task / subtask colours.** `TASK_HUES` (in `plotting_colors.py`) is the single source: Annotation
= crimson, Representation = amber, Sampling = lime. `TASK_COLORS` takes the base hues and
`SUBTASK_COLORS` splits each into one shade per subtask (base hue for the largest, lightening to
`SUBTASK_LIGHTEN_FLOOR = 0.5`), so a task and its subtasks can never drift onto different hues.
**`SOURCE_TYPE_COLORS` survives solely for the two pathogen panels' dot/cell colours**, and its hues
are chosen against what a reader sees *there*, not against the hub-wide Source Type split:
**Internal = periwinkle (33 dots), External = amber (16), Replicated = lime (2)** — ordered by dot count
in those panels, with the smallest category taking lime, which reads as an accent and keeps a two-dot
category from vanishing. Crimson is deliberately no longer in this dict: it is the Annotation task hue,
and these panels sit near the task-coloured ones. **The hub-wide ordering is the reverse**
(External 156 > Internal 45 > Replicated 7), so the hue ranking says nothing about the hub — only about
these panels' dots, and a caption should not imply otherwise.

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
grid alone crops landscape (1.23 at 60 mm), and the legend band beneath it brings the whole panel
back to square. Counter-intuitively, *reducing* the columns makes it worse, not better — a portrait
grid is narrower than the legend, which then fixes the content width while the extra rows add height.
Measured crop aspects at 60 mm: 13 cols 0.78, 14 cols 0.83, 15 cols 0.94, **16 cols 1.05**, 18 cols
1.22. Only 13 and 16 columns divide 208 without a ragged last row. Legend spacing (`handlelength`,
`labelspacing`, `columnspacing`, `borderpad`) is tightened from the matplotlib defaults because the
band's height is the knob that tunes the panel's aspect.

The key is **left-aligned** (`loc="upper left"` at `x=0`), not centred: its two columns hold labels of
different lengths, so a centred block leaves its left edge floating away from the grid's own left edge
directly above it.

**At quarter width the legend governs, so its labels are abbreviated** (`_LEGEND_ABBREV`: *Activity
prediction* → *Activity pred.*, *Property prediction* → *Property pred.*, *Similarity search* →
*Similarity*). With the full names the 2-column key measures **48.3 mm against a 45.0 mm panel**, so
the *legend* set the crop: the grid was squeezed to 35.9 mm inside a 47.7 mm page with dead space
either side, and the squares fell to **1.99 mm**. Abbreviated the key is 38.0 mm — narrower than the
grid — so the grid sets the width instead, expanding to 41.2 mm (**2.57 mm squares**), and the panel
crops to **46.23 × 46.09 mm**, square to 0.3%. The counts stay: they are the point of this legend and
not what overflowed. Measured alternatives, current spacing: full names 48.3 mm (over), counts
without parentheses 45.2 mm (over), abbreviating only `prediction` 44.2 mm (fits, but the key still
sets the width), both abbreviations **38.0 mm**. The abbreviations are legend-only — the full subtask
names still reach `task_subtask`'s tick labels, which sits beside this panel.

**The technical box row (`_HorizontalTaskPanel`, `_box_row_cells`).** Three **horizontal** panels
sharing one task axis, occupying **4/6 of the page width = 120 mm** between them and 30 mm tall:
`runtime_100`, `image_size` (both `TaskMetricBoxPlot`) and `output_dimension`
(`TaskOutputDimensionCirclesPlot`). Tasks run down the y axis, first on top; the metric runs along a
**log x axis** in all three (runtime 16–1626 s, image size 291–10242 MB, output dimension 1–5000 —
every one spans more than a decade, and a linear axis crushes the bulk against the floor). Per-task
quartiles go to `technical_metrics_summary.csv`. `docker_architecture` is not part of this
group; it is one of the three 25 mm donuts.

*Only the leftmost panel draws the task tick labels* (`show_y`), so the other two **cannot be placed
on their own** — they have no category axis. The label column is a fixed **14.2 mm**, so `runtime_100`
is wider than its neighbours by exactly that much and all three metric axes come out ~28 mm: equal
plot widths, which is the point. Sizing the row any other way (three equal panels) would leave the
first axis a third shorter than the other two and make the row read as though the metrics had been
scaled differently. Delivered crops **50.44 + 34.91 + 34.91 = 120.27 mm** against the 120.06 mm
target. `_box_row_cells` returns the axes width alongside the footprints so the caller can check its
labels against it, and the script prints `axes N mm vs xlabel N mm` for every panel each run — **a
metric label wider than its axes overhangs and *sets* the crop width, silently blowing the row
budget.** The runtime label is the tight one: 26.2 mm against 28.5 mm.

*Height 30 mm = 1 cell*, giving 17.7 mm of axes over 3 tasks — a **5.9 mm row pitch**, and at
`_BOX_WIDTH = 0.5` a **2.9 mm box**, which is deliberately the same mark weight as the 2.98 mm bars in
the subtask stacks next door. The swarm band (`_JITTER = 0.20`) is as wide as the box: the old vertical
panels used 0.10, which at this pitch would pile 131 Annotation points into a ~1 mm line.

*The `-1` sentinel is "never benchmarked", not zero.* Airtable stores `-1` where a runtime benchmark
was not run, so non-positive values are skipped and **never imputed**. Any task left with nothing
keeps its slot, labelled "not measured" rather than omitted — a missing box is information, a missing
category is a misreading.

***`n` is no longer in the tick labels.*** The old vertical panels put each task's `n` in its own tick
label; a shared axis can carry only one label set, and coverage differs per metric — runtime is
**129/57/10**, image size and output dimension are both **131/58/19**. One labelled `n` would
therefore be wrong for two panels out of three. It moved to `technical_metrics_summary.csv` and the
run log, which means **a caption must state that the runtime box for Sampling rests on 10 of 19
models** — that is a real transfer of load-bearing information out of the figure.

**`output_dimension` is circles, not a box** (`TaskOutputDimensionCirclesPlot`). Output Dimension is
heavily **tied** — 68 of 131 Annotation models output a single value (100 of them fall in the 1-9
bin, which is the circle the panel draws), and 100 and 1000 recur across Representation and
Sampling — so a swarm piled onto a handful of x positions and its visual density
described the jitter rather than the data. Instead: one circle per (task, decade), **area proportional
to the number of models**, at the bin's *geometric centre* so it sits between the two decade ticks
that bound it. A circle centred on the 10² tick would read as "exactly 100", which is a real and
common value in this column.

The 9 non-empty bins span **2 to 100 models** — a 50× range, so 7.1× in diameter.
`_MAX_DIAMETER_MM = 4.6` pins the largest circle just inside the 5.9 mm row pitch, leaving the
smallest at ~0.7 mm. **There is no size key** — at 34 × 31 mm there is nowhere to put one — so the
panel is for the pattern (Annotation outputs one value; Representation spreads across three decades;
Sampling sits at 10²–10³) and the exact counts live in `output_dimension_bins.csv`. It is a drop-in
replacement: same task axis, same footprint, same orientation, so it swaps with the box version
without touching the row.

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

**Circles: `_donut`.** The one helper. It puts **nothing on the wedges** and leaves labelling to a
legend the caller supplies, because every panel that uses it is too narrow to label in place — at
25 mm, `"Non-commercial 1 (<1%)"` is 20.07 mm of text and the name alone is 13.46 mm, so labels on two
sides would consume the whole panel before the ring got any. (An outside-labelled `_pie` variant with
a `_share_labels` helper existed for the retired `license_class_pie` / architecture-pie panels; both
went with those panels.)

`_donut` also **pulls the axis limits in to ±1.02**. `ax.pie` always sets them to ±1.25 to leave room
for the outside labels it normally draws; with no such labels that is 20 % of dead margin on every
side, which both shrinks the ring and opens a ~3.5 mm gap above a legend anchored below the axes.

**The three donuts are one family** (`DonutPlot`): `license_class_donut`, `docker_architecture` and
`biomedical_area_donut`, each **25 mm wide** — 75 mm for the set — with the ring, the total in the
hole, and the key beneath. One implementation, so they cannot drift apart; a subclass only prepares
`(labels, values, colors)`.

**The ring is pinned, not inherited** (`DONUT_RING_MM = 19.6`, `DonutPlot.pin_ring`). Left to itself the
ring's size is a side effect of label length: `tight_layout` shrinks the axes to fit a legend wider than
it, so the two-row architecture key squeezed its own ring to **18.18 mm against the other two panels'
19.81 mm** — an 8 % difference driven by nothing but the string `"AMD64+ARM64"`. `pin_ring` runs after
`tight_layout` (it needs the axes' final size in mm) and scales the *axis limits*, not the radius, so
wedges, hatching and hole all scale together. Because a pie axis draws no frame, this also makes the
**crop widths uniform** — the crop is ring + savefig pad in every case, **24.89 mm** for all three.

`DONUT_RING_MM` must stay **at or below the narrowest natural axes width**, or the ring overflows its
axes and is clipped. The script prints each panel's natural axes width against the pin every run and
says so loudly if one is short — that check is what caught the architecture panel twice, first at
18.18 mm and again at 18.71 mm after the legend spacing was tightened. It was fixed properly by
shortening the labels to **`AMD only` / `AMD + ARM`**, which restores the full 19.81 mm; the full
`AMD64` / `ARM64` are the Docker platform identifiers and belong in a caption if precision matters.

*Heights still differ and that is structural.* The panel is ring + legend band, so the two-entry
architecture donut comes out **24.89 × 30.45 mm** against the others' **24.89 × 35.00 mm**. Widths and
rings now match exactly, so **align the set on the rings or the left edges; the panel boxes differ in
height only.**

*No hue means two things across the set.* Licence holds turquoise / periwinkle / silver / fuchsia, so
architecture moved off turquoise + periwinkle to **cobalt (AMD + ARM) / tangerine (AMD only)** —
both substantive, since silver would cast x86-only as a residual rather than a real build target.
Tangerine is a plain categorical hue here, as it is for Single Point in script 02; **a caption must not
read it as flagging x86-only builds as a problem.** Architecture is ordered base-capability first
("x86 only, then also ARM"), not by size.

*The biomedical donut spends no colour at all.* Every wedge is the Annotation crimson, because every
model in that panel *is* an Annotation model — colour would encode nothing. Groups separate by
**pattern** (`BIOAREA_GROUP_HATCH`): solid for the largest, then progressively lighter-inked patterns
so ink ordering matches size ordering, with the cross-hatch on the catch-all where it reads as "mixed".
Patterns are white over the crimson (matplotlib hatches take the patch edge colour), so a wedge still
registers as red at a glance, and the legend swatches carry the same patterns so the key cannot
describe a mark the ring does not draw. **Hatch density is the repeat count** — nine repeats, because at
three the ring showed finger-thick stripes and dots as thick as the ring itself — paired with
`hatch.linewidth = 0.35` set in `plotting_base` (matplotlib's 1.0 default is the weight of a bar
outline and turns a fill into damage).

*Legend rows carry the count but **not** the share.* That is what sets the panel size: the legend must
stay narrower than the ring or `tight_layout` shrinks the axes to fit it and the ring collapses, and at
25 mm `"AMD64 + ARM64 129 (62%)"` was 25.94 mm against a ~20 mm budget. Dropping the share costs little —
it is exactly what the ring already encodes, with the hole giving the total to divide by — where
dropping the count would lose the only number nothing else states. `ARCH_DISPLAY` also drops the spaces
around the `+` for the same reason.

*The hole shows a **model** count, not always the wedge sum.* `DonutPlot(total=...)` overrides it, and
`biomedical_area_donut` uses it: its groups sum to **94 across 92 models** because two models carry
areas in two groups, so the hole reads 92 while the legend rows add to 94. The discrepancy is real and
belongs in the caption rather than being hidden by quietly showing the sum.

**`license_class_donut` shows the four reuse classes, not the ten licences** — Permissive 104 (50 %),
Copyleft 76 (37 %), Not recorded 27 (13 %), Non-commercial 1 (<1 %), colours matching the `license`
bar panel. That is a hard limit, not a shortcut: four of the ten licences cover exactly one model
each, which on 208 models is a **1.7° wedge**, invisible and impossible to label or tell apart from
the other three. Per-licence detail is the bar panel's job. Wedges run in descending order so the
Non-commercial hairline finishes at the 12 o'clock start line rather than sitting between two large
wedges, where it would read as a rendering artefact.

*Why a legend rather than labels around the ring.* At 25 mm, labelling around the circle is
geometrically impossible, not merely tight: `"Non-commercial 1 (<1%)"` is **20.07 mm** of text and the
name alone is **13.46 mm**, so labels on two sides would take more than the whole panel before the ring
got any. It is also what rescues the 1.7° wedge, which can carry no label of its own at any panel size.

*Calibrating the width needs a full script run.* stylia wipes matplotlib's font cache on import, so the
**first figure of a process draws with fallback text metrics** and measures ~0.6 mm small — an isolated
sweep of this panel gave 29.88 mm where the real run gave 30.48 mm at the same footprint.
`_DONUT_CROP_PAD_MM` is calibrated against the script. The same trap applies to any other panel
calibrated by measurement.

`docker_architecture` is the two-wedge member of the donut set: **AMD only 79 (38%)** in tangerine,
**AMD + ARM 129 (62%)** in cobalt.
The field only ever holds `AMD64` or `AMD64,ARM64` — there is no ARM-only build, so this is "also built
for ARM" versus "x86 only". It is a **snapshot, not a trend**: dual-arch is 45% among models
incorporated in 2021 and 77% among 2026 ones, so the 62% is accumulated stock rather than current
build practice, and a caption must give the metadata snapshot date.

**Legends across this figure's subtask panels.** `task_subtask` uses `SUBTASK_COLORS` (not
the parent-task base hue), so every subtask appears there as a labelled bar in its own colour —
that panel *is* the subtask key, and it therefore carries no legend of its own. The two `_by_subtask`
stacks carry none either: at ~52 × 25 mm each a 6-entry key would be taller than the plot,
and they now sit directly beside the waffle, which supplies the key. The **waffle** has its own
(2 columns beneath the grid, abbreviated labels — see above). **Consequence for layout: if you use a `_by_subtask`
stack, `task_subtask` or `task_subtask_waffle` must travel with it, or render a standalone subtask
key.**

### 01b_community_stats.py (`save_timeline_figure`, `src/plots_timeline.py`)

**One** panel, `hub_timeline`, `cells=(1.05, 6)` = **180 × 31.5 mm**, containing four stacked
tracks on a single shared year axis. The only figure in the repo that is deliberately not one
chart per file — see the shared-axis note above, which prescribes exactly this. Measured page
after tight crop: **181.6 × 32.8 mm**.

| track | y measure | hue | source |
|---|---|---|---|
| Models | cumulative | cobalt | read from `output/01_models_metadata/models_over_time_by_task.csv` |
| People | cumulative | turquoise | aggregated in-script from `commit_weeks.csv` |
| Commits | per month | lime | aggregated in-script from `commit_weeks.csv` |
| Issues | per month | tangerine | aggregated in-script from `org_participation.csv` |

Only the Models track is read; the other three are aggregated by the script itself from
`data/raw/github_stats/`, scoped to the **247** repos that are the Model Hub (`repo_set.csv`) and
with `default.GITHUB_BOT_ACCOUNTS` excluded. This step is the merge of the former 08 and 09,
reduced to what the timeline needs — the eight community panels and `src/plots_community.py` went
on 2026-08-04.

**Stocks and flows share the figure, and the caption must say so.** Models and People are
cumulative — the height *is* the size of the hub. Commits and issues are per-month counts — the
height is a rate. **Vertical comparison across tracks is meaningless**; horizontal comparison is
the entire point of the shared axis.

**"People" is distinct COMMIT authors (107).** The issue/PR-author series (335) was dropped from
this figure; it survives in `01b_timeline_series.csv` and `01b_snapshot.txt`. A caption saying
"people" here is claiming people who wrote code, which is about a third of the people who took
part — worth being precise about.

**Alignment is the reason this is one file.** All four axes come out at the same figure-fraction
left and right edge, on one `xlim`. Four separately saved panels cannot achieve that:
`bbox_inches="tight"` sizes each file's tick-label column to its own content, which is what puts
the two subtask stacks 0.25 mm out of register.

**No workflow-runs track.** GitHub deletes workflow-run records after ~13-14 months, so that
series only starts 2025-06; on a six-year shared axis it is empty for three quarters of its width.
The collector no longer fetches them, so the series exists nowhere. `HubTimelinePlot` still
supports a `no_data_before` key that shades the unavailable stretch and labels it *"no data
retained"* — currently unexercised, kept because any truncated series added later needs it. An
empty stretch of track is not the same as a stretch of no activity.

#### Everything is sized around a 4.4 mm track

At 180 × 31.5 mm the shared x band costs ~12 mm once for the whole figure and the four tracks
split what is left, giving a **measured 4.4 mm** of drawing height each — roughly a quarter of the
smallest panel anywhere else in the repo. Five consequences, all forced:

- **Axes are placed by hand with `add_axes`, not as subplots.** `stylia.save_figure` calls
  `plt.tight_layout()` unconditionally; tight_layout separates gridspec rows by `h_pad` (default
  1.08 font-size units, ~2.3 mm here) and recomputes from the gridspec every time, so
  `subplots_adjust(hspace=...)` is silently discarded and stylia's wrapper takes no argument to
  change it. Axes made with `add_axes` have no subplotspec and tight_layout leaves them alone.
  That is the only way to set the inter-track gap exactly — and it returns the padding
  tight_layout was spending, which is worth ~2.3 mm per boundary.
- **`_TRACK_GAP = 0.06`**, sized by measurement, not by eye. With the tracks touching, a track's
  `0` label and the top label of the track below are centred on the same spine and **overlap by
  2.4 pt**. Sweeping the gap and measuring the label bounding boxes: 0.03 still overlaps
  (+0.3 pt), 0.05 just clears (−1.1 pt), 0.06 clears with margin (−1.7 pt). Every mm of gap comes
  straight out of the tracks — the crop stays 32.3 mm regardless — so it is kept at the minimum
  the labels allow.
- **Two y ticks per track, and the intermediate ones are removed rather than blanked.**
  `MaxNLocator(3)` still runs, but only to choose a round top value (160 rather than 215, 1,000
  rather than 1,366); `_thin_y_labels` then keeps the first and last positions. **Only ticks
  inside the view can be kept** — `MaxNLocator` routinely places one above the data (0/80/160/240
  for a 215 maximum), and taking the raw last entry once labelled a tick that is never drawn,
  which made the top label vanish.
- **Horizontal track labels**, against the house default of rotated ones (`rotation=0`,
  `set_label_coords(-0.055, 0.5)`). Rotated, "Commits" is ~7 mm of type standing in a 4.4 mm
  track: the four labels overlapped into a single unreadable column reading
  `IssuesCommitsPeopleModels`. Laid flat each is ~2.1 mm tall, at the cost of ~8 mm of the 180 mm
  width. Because every track uses the same `_YLABEL_X` they are left-aligned by construction;
  `fig.align_ylabels` is kept as well since it costs nothing. **Without it the Commits label
  sticks out**, because each label is otherwise placed against its own tick-label column and
  `1,000` is wider than `200`.
- **Within-track breakdowns are summed away, and now discarded upstream.** The sources carry
  category splits (task; core/model repo group); a 2-3 segment stack inside 4.4 mm is unreadable.
  The repo-group split is no longer computed at all — the script reads `repo_set.csv` for
  membership only. The Models track's task split has no time-series panel of its own — step 01
  shows the task composition as a snapshot (`task_subtask`) — so its CSV is the only place it
  survives.

**Vertical gridlines only** (`ax.yaxis.grid(False)`, `ax.xaxis.grid(True)`). The year rules run
down the whole figure and are how a feature in one track is located against another — they are
what the shared axis is *for*. Horizontal rules cross them every couple of mm in a track this
short, and the resulting mesh reads as texture rather than as a scale; the two y labels carry the
range instead.

**X ticks on the bottom track only** — marks as well as labels
(`tick_params(bottom=False, labelbottom=False)`). With a ~1.8 mm gap an upper track's tick marks
would hang most of the way into the track below and read as marks on *its* data. This is stricter
than `StackedFieldBarPlot.show_x`, which keeps its marks because its panels sit further apart.

**Tune the height against the CROP, not the nominal `cells`.** The nominal figure is 29.4 mm but
crops to 32.3 mm, because the outermost tick labels sit outside the axes rectangle. Choosing
`rows` from the nominal number overshot the target by 10% on the first attempt. Use
`plotting_base.crop_size_mm`, whose height is exact.

**This height is near the tight-layout floor.** With every gridline labelled, matplotlib refused to
run tight_layout at all ("cannot make Axes height small enough to accommodate all Axes
decorations"). `save_timeline_figure` catches that message and reports it as one plain line — a
regression guard. It deliberately does **not** report "not compatible with tight_layout", which
`add_axes` always raises and which is the mechanism working as intended.

**One hue per track, at full strength, and it is decorative.** `TIMELINE_HUES` /
`TIMELINE_COLORS` in `plotting_colors.py` — cobalt, turquoise, lime, tangerine. Assigned
**explicitly**, not by zipping `_CATEGORICAL_HUES`: People takes turquoise, which the pick order
would have given to a fifth track. Spending the repo's default/positive hue on a decorative track
is only acceptable because no track here encodes a category — each is named by its label and no
two share a scale. `TIMELINE_FILL_ALPHA = 1.0`: the fills are the **full hue, not a wash**, since
at 3.4 mm a tinted band is barely distinguishable from white. Hues are stored as **names**, not
resolved colours, so a nested second series can be asked for via `hue(name, lighten=...)`;
`distinct_colors()` returns tuples that cannot be lightened.

**Year ticks every year** (`YearLocator(1)`). At 180 mm all seven labels fit with room to spare;
a 45 mm panel would need step 2.
Only the bottom track keeps its tick *labels*; the other three keep the ticks so they register
against the year grid below, and drop the text — the `StackedFieldBarPlot.show_x` idiom.

**`MultiPanelPlot` does not set `self.ax`.** `BasePlot.label` / `.legend` / `.ref_line` are
therefore unavailable in this module; labelling goes through `stylia.label(ax, ...)` directly, as
in the five `plots_chembl_curation.py` subclasses.

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

**Top level of the output dir** holds the condensed cross-pathogen figures:
`pathogen_activity_ratios` (2,4) and `pathogen_consensus_auroc` (2.04,1.94). Per-pathogen panels are
intermediate results and live in `individual_plots/`.

Two panels per pathogen (15 pathogens, **196 models**), written to
`output/03_chembl_models_performance/individual_plots/{png,pdf}/` with their own
`figure_cells.json`. These 30 panels are **intermediate results, not paper figures** — they exist to
inspect every model, and the condensed cross-pathogen figures live at the top level of the output
dir. Note the grids show all of step 09, which includes 3 models step 10 discarded for AUROC < 0.7 —
they are drawn unmarked, so any caption claiming "the 193 hub models" would be wrong; the `retained`
column of the summary CSVs is the source of truth. See `scripts/README.md`.

**`pathogen_activity_ratios`** — footprint `(2,4)` = **120 × 60 mm**. One dot per modelled dataset,
with the active fraction on the y axis (0–1, linear) where a reader can actually resolve it and size
demoted to **dot area**. A dashed silver line marks the 0.5 balance point. Tick labels are
genus-abbreviated, single-line and upright; columns are ordered by dataset count descending.
Data: `dataset_sizes.csv`, derived from `10_reports/10_reports.csv` (the 193 step-10 keeps).

A superseded sibling, `pathogen_dataset_sizes` (3,6), swapped the two encodings — size on a log y axis
with a per-dataset pie carrying the balance. It was dropped because a 2.5 mm two-slice pie cannot
resolve 45 % from 50 %, which is exactly the comparison the panel had to support, and because it cost
the full page width to keep *P. falciparum*'s 51 pies separable.

*Added negatives are excluded from both the size and the ratio.* `n_compounds` counts them, so size
is `n_compounds − n_added_negatives − n_added_decoys` and the active fraction is `n_positives / size`.
That reproduces the curation pipeline's own `n_mol_after` / `ar_after` **exactly** — verified for all
151 models that join to `25_pool_summary`. Only 54 of 193 models were given negatives, so the medians
barely move (1115 → 1098 compounds; active ratio unchanged at 0.365), but for those 54 the change is
large: `mtuberculosis/DR_0012` goes from 2450 compounds at 0.50 active to **1411 at 0.87**. Decoys are
zero for every model, but are subtracted anyway so the definition does not silently depend on that.
Note **54 of 193 datasets are majority-active** once added negatives come out, and none are 0 % active.

*Dot area is affine in √size, not proportional to size.* `SIZE_REF_AREA = 22` pt² at
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

*Variable column widths.* Widths are proportional to `n_datasets ** 0.25`, floored at `MIN_WEIGHT =
0.55` of the mean. Equal columns waste the panel — three pathogens have a single dataset and need no
spread at all, while *P. falciparum* needs every millimetre for its 51. The fourth root rather than the
square root because width only has to *separate* the dots, not encode the count: √n gave the singletons
a 2.5 mm slot, too narrow for a tick label. At 120 × 60 mm — **44 % of the original 180 × 90 mm area** —
*P. falciparum* still gets **10.9 mm** against the 11 mm it had at full page width, and the singletons
4.1 mm. Jitter and mean-bar length are fractions of each column's own half-width.

*Dot sizes are tied to the footprint.* Marker areas are absolute (points²), so condensing the panel
without shrinking them proportionally just crowds the dense clusters — `SIZE_REF_AREA` went 30 → 22
alongside the 180 → 120 mm narrowing. Rescale it whenever `cells` changes.

*Dots are drawn largest-first, in a single scatter call.* Every dot is opaque; overlap is resolved by
draw order, not by alpha (stacked pale dots under transparency read as a darker, i.e. different, hue).
One call across all pathogens rather than one per column, because matplotlib takes `zorder` per artist,
not per point — a 96-compound dot can otherwise disappear under a 300,000-compound dot from the
neighbouring column. Sorting is `np.argsort(-area, kind="stable")`, so the small dots land on top.

*Tick labels are upright, not italic.* They are single-line (no
`(n)`); the per-pathogen count lives in `dataset_sizes.csv` and the
number of dots shows it anyway. Rotated labels on a shared baseline clear each other when
`column_width × sin(angle)` exceeds the label line height, so dropping the second line halved the
clearance needed and is what let the columns narrow. At `LABEL_ROTATION = 40°` the tightest column
gives 3.4 mm of perpendicular clearance against a 2.1 mm line height.

*Mean bars sit UNDER the dots* (`zorder=2.5`). Where a pathogen has one dataset the mean lands exactly
on its own dot, and a 1.4 pt bar drawn on top swallowed it whole — *Campylobacter*'s single dataset is
96 compounds, a 0.6 mm dot.

*Palette.* Dots take the pale `FILL_LEVELS = (0.62, 0.38)` palette so heavy overlap stays readable,
outlined and topped by the darker `ACCENT_LEVELS = (None, 0.78)` of the same hue. The second accent
level sits close to the base so fill-vs-accent contrast is comparable in both passes (0.62 against
1.0, 0.38 against 0.78); that makes the two accent levels similar to each other, which is fine because
the *fills* carry the hue-repeat distinction and accents only outline. **Both levels are
needed in each call**: `distinct_colors` uses its second pass to get past the 9-hue limit, so a single
lighten value would make pathogen 10 identical to pathogen 1. The two calls index hues identically, so
a dot and its mean bar read as one category.

*Size key (`plotting_utils.nested_size_legend`).* Nested circles sharing a bottom tangent, drawn as
`Ellipse` patches with the point radius converted **separately per axis** so they stay round whatever
the aspect or scale. It must be called last — it reads the live transform, so the axis limits have to
be final. Label rows are **de-collided and reached with leaders**: nesting packs the ring tops within
`2 × rmax` of each other, about 3.4 mm for this decade key, while four labels need ~8 mm of line
height, so anchoring each label at its own ring's top guarantees overlap. Pass the block's caption via
the `title` argument rather than adding your own text: only this function knows where the spread block
ends, so a caller guessing an offset collides with the lowest label.

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

`s` takes a scalar for uniform pies **or an array to size each pie individually** — each wedge reads
`s` from its own point, so slice and radius stay independent encodings. A **full pie is drawn as a
plain `"o"` marker**, not a 360° wedge: the wedge path closes through the centre, so at `frac = 1` its
two radii coincide at 12 o'clock and stroke a visible seam across an otherwise solid disc. That case
is common wherever a matrix has a diagonal (`active_overlap_containment` in step 05).

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
  `PipelineFunnelPlot` and `HitPromiscuityPlot.Y_FLOOR`.
  **This bug is invisible in review.** `HitPromiscuityPlot` shipped with `bottom=0` for a while: its
  PNG looked perfect because the raster clips, while the PDF carried path vertices at
  **y = −28,828 pt on a 131 pt page**. Checking the PNG proves nothing.

  **Audit at the source, not the PDF.** The reliable check is a property of the code: *a `bar`/`barh`
  call on a log-scaled axis with no explicit `bottom`/`left`*. Grep for `set_yscale("log")` /
  `set_xscale("log")` and inspect the bar calls near each one. Auditing the rendered PDF instead is a
  trap — a hand-rolled scan of path coordinates gives **false positives**, because matplotlib wraps
  rotated text and grouped artists in `q … cm … Q` blocks, so a raw coordinate is in local space, not
  page space. A 90° rotation turns a raw y into a page x, which made the perfectly healthy
  `02/lowdata_auroc` look like it had geometry 121 pt below its page. Doing it properly needs a real
  parser (`pikepdf`), which is not currently a dependency.

  As of the last audit there are **four log axes** in `src/`: `PipelineFunnelPlot` and
  `HitPromiscuityPlot` (bars — both use a finite floor) and `TaskMetricBoxPlot` /
  `PathogenDatasetSizesPlot` (boxes and scatter, which have no zero baseline and cannot be
  affected). Only bar marks are exposed. (A fifth, `PublicationLagPlot`, went with the retired
  growth strip; it was verified clean by decompressing its PDF's *content* stream and checking no
  path coordinate exceeded the MediaBox. Note a naive byte scan of every stream in a PDF reports
  false hits — the embedded font subsets are binary and full of large integers — so restrict any
  such check to the ASCII content stream.)

- **Run in the stylia env:** figures must run in a conda env with stylia installed (e.g. `paper`),
  not `base`.
- **Small multiples need explicit sizing:** `stylia.create_figure(nrows, ncols)` sizes the *whole
  figure* to the format default, so an N×M grid gets N×M squashed panels. Always pass
  `width=cols/6, height=rows/6` and, for aspect-critical charts (ROC, PR, calibration), also
  `ax.set_aspect("equal", adjustable="box")`.
- **Data provenance:** record the source snapshot/version (e.g. ChEMBL release) in the script
  docstring; stage summaries via `scripts/00_download_data.py`, never ad-hoc copies.

### 05_euopenscreen_validation.py (`save_euopenscreen_figures`, `src/plots_euopenscreen.py`)

Six panels are sized for the paper figure at **`SMALL_SQUARE = (1.5, 1.5)` = 45 × 45 mm**, arranged as
**two rows of three** below the step-03 row. Everything else in this script keeps its roomier
footprint — those are inspection figures, not paper panels.

| row | panels | crop widths |
|---|---|---|
| A | `hit_promiscuity`, `euos_overlap`, `exclusive_hit_model_rank_percentile_dedup` | 41.1 + 45.0 + 39.9 = 126.0 mm |
| B | `consensus_max_percentile_by_activity_dedup`, `hit_exclusivity_auroc`, `submodel_auroc_summary` | 33.4 + 50.0 + 44.4 = 127.8 mm |

Both rows still clear a 183 mm page row. Re-measure them from the saved PDFs whenever a footprint
moves — this table has gone stale twice already.

#### The six landscape panels share one axes height (2026-08-05)

**Every one of the six is tuned so its AXES HEIGHT equals `euos_overlap`'s 25.81 mm.** They are laid
out together, so what has to line up is the **plotting rectangle**, not the file. Achieved: 25.81
(overlap), 25.88 (promiscuity), 25.85 (rank), 25.87 (submodel), 25.73 (consensus sum), 25.71 (events)
— a **0.17 mm spread**, i.e. ±0.09 mm about the reference.

| panel | constant | cells | axes w × h (mm) | crop (mm) |
|---|---|---|---|---|
| `euos_overlap` *(reference)* | `WIDE_OVERLAP` | (1.2940, 1.4617) | 23.01 × 25.81 | 45.0 × 39.9 |
| `hit_promiscuity` | `WIDE_PROMISCUITY` | (1.2723, 1.3298) | 26.74 × 25.88 | 41.1 × 39.4 |
| `exclusive_hit_model_rank_percentile_dedup` | `WIDE_RANK` | (1.2961, 1.2869) | 26.74 × 25.85 | 39.9 × 40.3 |
| `submodel_auroc_summary` | `WIDE_SUBMODEL` | (1.2547, 1.4413) | 22.93 × 25.87 | 44.4 × 38.9 |
| `hit_exclusivity_events` | `WIDE_EVENTS` | (1.3626, 1.3293) | 24.05 × 25.71 | 41.1 × 42.2 |
| `consensus_sum_by_hit_class_dedup` | `BOX_EVENTS` | (1.3205, 1.2116) | 24.00 × 25.73 | 38.0 × 44.3 |

**The crops are deliberately unequal** — 38.0 to 45.0 mm wide, 38.9 to 44.3 mm tall — because each
panel spends a different amount of itself on chrome (a second x axis, an AUROC header, 45° class
labels). Equalising the crops would push the data rectangles *out* of register, which is the opposite
of the point. Do not "fix" the crop spread.

**Widths are matched in pairs, not globally**, since panels sit side by side and only heights have to
agree: `euos_overlap` ↔ `submodel_auroc_summary` (23.01 / 22.93), `hit_promiscuity` ↔
`exclusive_hit_model_rank_..._dedup` (26.74 / 26.74), `hit_exclusivity_events` ↔
`consensus_sum_..._dedup` (24.05 / 24.00).

**`euos_overlap` is the reference because it is the constrained panel.** 40 × 30 mm was tried and only
works for `hit_promiscuity`: `euos_overlap` spends its height on chrome twice over — *two* x axes, each
with ticks and a label — so at a 30 mm crop the axes came out **9.5 mm** tall for 7 genus rows that
need 17.8 mm at the house 6 pt, and the labels collided outright. It cannot give height back, so
everything else comes to it. The general rule: a twin-axis panel pays its chrome twice, so its height
floor is set by the *rows*, not by the footprint you would like. Its width also carries a deliberate
**+2 mm** (written `1.395 + 2 / 30`, 1 cell = 30 mm) over the tuned 45.0 mm crop so the panel breathes;
at 40 mm the library axis loses its `50000` tick and degrades to endpoints only.

**How to re-tune** when a panel's chrome changes (a tick label, an axis label, a rotation, a legend):
the relation is `axes = cells/6 × 180 mm − chrome` with chrome near-constant, so correct
`cols += dw/30`, `rows += dh/30`, iterate 4–6 times, then **local-search the last ~0.1 mm** — the
analytic solve lands ~0.2 mm out because the tight bbox quantises to text extents (~0.25 mm steps).
Measure the **axes** via `ax.get_position()` scaled by the figure size, *not* the crop: `crop_size_mm`
predicts height exactly but runs ~1.2 mm under on width, and `pdf_page_mm` measures the page, not the
plotting area. `hit_exclusivity_events` and `consensus_sum_..._dedup` are both at a quantisation floor
— searched at a 0.001-cell (0.03 mm) grid with no closer value on offer — which is why they sit
0.1 mm low rather than exact.

**Two older sizing flavours are now retired** and should not be reintroduced for these six: matching
another panel's *crop* (`WIDE_RANK` used to match `hit_promiscuity`'s), and fitting under a *bound*
(`WIDE_SUBMODEL` used to be "at most 45 × 40"). The shared-axes-height rule supersedes both; the
submodel panel satisfies its old bound anyway at 44.4 × 38.9 mm.

**Three per row, not six.** Six at 45 mm is 270 mm against a 183 mm row — only four fit. Two rows of
three leaves ~42 mm of gutter to distribute and puts the page at 62 (step-03 row) + 45 + 45 = **152 mm
of the 170 mm budget**. Squeezing six into one row would mean 30.5 mm each, below the repo's
documented 45 mm floor, and these panels carry 7 genus tick labels.

**Legends are stripped and emitted as standalone `*_key` panels** (`KeyPanel`, following script 02's
`curation_outcome_legend`): `euos_overlap_key` (4 entries, colour = leakage and hatch = which axis),
`exclusive_hit_model_rank_key` (7 organisms with their hit totals), `hit_exclusivity_auroc_key` and
`submodel_auroc_summary_key` (2 entries each). Place each once in Illustrator. Panels take
`legend=False`; `grouped_hbar` gained the same switch. Keys and panels share their definitions
(`SUBMODEL_KEY_ENTRIES`, `euos_overlap_handles()`, and the rank key built from the plot's own
`segment_totals`), so a key cannot drift from the marks it explains.

**Axis labels had to become terse.** At 45 mm a sentence-length label doubles the panel's crop:
`consensus_max_percentile_by_activity_dedup` came out **77 mm tall** with its original
*"library percentile of best-model score (7 models), no training-set compounds"* y label, because a
rotated label that long sets the crop height. The short forms below carry the model count, dedup
status and full phrasing **into the caption instead** — a real transfer of information out of the
panels, so the caption is now load-bearing.

| panel | was | now |
|---|---|---|
| `consensus_max_percentile_…` | "library percentile of best-model score (7 models), no training-set compounds" | "Best-model percentile" |
| `hit_promiscuity` | "number of pathogens the compound is a hit for" | "Pathogens hit" |
| `submodel_auroc_summary` | "AUROC on own EU OpenScreen assay (dedup)" | "Own-assay AUROC" |
| `exclusive_hit_model_rank_…` | "rank of the hit pathogen's own model (of 7)" | "Own-model rank (of 7)" |

The *now* column shows the rendered label. Panels may write these in lower case in code — the
sentence-case wrap in `plotting_base` capitalises them (see the house-style list above).

**`exclusive_hit_model_rank_*` carries no grid and no chance line.** The dashed neutral rule marking
`n_chance` (the count expected if the own-model rank were random, `n_total / 7`) was removed: nothing
on the panel said what it was, a bare grey rule reads as a threshold or an axis break, and it ran
straight through the bar-top numbers. The value is still in the summary CSV's `n_chance` column, so it
belongs in the caption where it can be named. The y grid went with it, for a structural reason worth
generalising: **a bar whose top lands near a gridline puts its value label onto that line**, and no
fixed label padding fixes it, because whether it collides depends on where the counts happen to fall
against the tick locator — 3 pt of padding cleared `19` and `41` but left `37` sitting exactly on the
40 line. Since every bar is annotated with its exact count, the grid was redundant anyway; dropping it
makes the labels legible by construction rather than by luck.

**`euos_overlap` has no gridlines on either axis**, set explicitly rather than left to the rcParam
default (`axes.grid = True`), so re-enabling one is a deliberate act. With two x scales on one panel a
vertical guide is ambiguous by construction — nothing tells the reader whether it belongs to the
library scale or the actives one — and a grey vertical crossing a stacked bar additionally reads as a
division in the data.

Getting there took two passes worth recording. The grid was first kept and merely reordered, which
turned out to be harder than it looks: `twiny()` axes draw *after* the parent, so the library axis's
gridlines sat **on top of the actives bars**, and `set_axisbelow` cannot fix that — it orders a grid
against its own axes' artists, and those bars belong to the other axes. The working fix was to switch
both grids off and re-draw the library-tick guides on the *bottom* axes at `zorder=0`. That rendered
correctly and was still dropped, because thin and correctly-layered did not cure the ambiguity. Bars
this legible do not need a grid at all.

The actives axis keeps `MaxNLocator(nbins=4)` — the default locator offered just 0 and 200 over a
0–397 span, too coarse to read a bar against.

**Axis names sit inline, at the left end of their own tick row** (`_inline_axis_labels`), not centred
beyond it — `Actives` on the bottom row, `Full library` on the top. A centred label costs a whole text
row *per axis*, and this panel has two, so ~8 mm of its height went to naming them; set on the tick
row they occupy the empty corner left of the numbers, under the column the genus labels already claim.
That single change took the plotting area from ~19 mm to 25.8 mm at an unchanged 40 mm crop — by far
the cheapest height available on a twin-axis panel.

Two details make it work. The row's y is **measured**, not assumed: the label is pinned in axes
coordinates to the mid-height of the drawn tick labels, so it stays on their line whatever the tick
padding or font size (hence `tight_layout` + `draw` before pinning). And `LABEL_X = -0.10`, well left
of the -0.02 tried first: the `0` tick is *centred* on the axes' left edge, so it reaches into
negative axes coordinates itself and a name set too close lands on top of it.

Dropping the centred labels means `label()` is called with everything empty purely to clear stylia's
placeholders — omit it and the panel ships with a literal **"Y-axis / Units"** down its left edge,
since the inline names are set on the axis objects directly and never go through `label()`.

`ACTIVES_TICK_BINS` came back down to **3**: once the genus labels take their ~14 mm out of a 40 mm
panel, 4 tick numbers run together (`0 100200300`).

**Never set a per-panel label size.** stylia puts *every* label at one size — `axes.labelsize`,
`xtick.labelsize`, `ytick.labelsize`, `axes.titlesize` and `legend.fontsize` are all 6.0 — and a panel
carrying two label sizes reads as a mistake, not as emphasis. This panel's genus ticks ran at 5 pt for
a while, to buy height back when it was short of it, which left them visibly smaller than the numbers
on the other axis. Moving the axis names inline freed enough room to drop the override: the axes gives
each of the 7 rows 3.7 mm against a 2.5 mm line height at 6 pt.

If a panel seems to need smaller labels, it needs a different layout — shorter label text, a bigger
footprint, or chrome removed (as here). `st.FONTSIZE_SMALL` (5) remains correct for *annotations* —
bar-value callouts, key text, in-axes legends — which are a different class of text from labels.

*Hatch density is tied to bar height.* The library bars carry `HATCH = "//////"` at
`HATCH_WIDTH = 0.4`, not matplotlib's `"//"` at the default 1.0 pt stroke. At 45 mm a bar is ~1.5 mm
tall, and a sparse hatch lands one or two strokes on it — that reads as a printing defect, not a
texture. The stroke has to be thinned as the hatch tightens, or the slashes merge into a solid block
and the fill colour underneath is lost. `set_hatch_linewidth` is per-patch (matplotlib ≥ 3.10), so no
global rcParam is touched. `euos_overlap_handles()` pulls both constants off the class, so the key
swatch cannot show a coarser texture than the bars.

**`hit_promiscuity` bars fade with `n_pathogens` — the x axis, not the bar height**
(`plotting_colors.count_shades`, the ordinal counterpart of `auroc_shades`). Palest for the
organism-specific singletons, darkest for the pan-active 7-pathogen bin, so the gradient runs the same
direction as the thing the panel is about. Shading by the bar's own *count* was tried first and reads
backwards — it makes the singleton bin the darkest mark on the panel, and it merely restates the y
axis. **Pick the variable the gradient is supposed to mean, not the one the mark already shows.**
The exact count is annotated on every bar, so colour is never the sole encoding.

`count_shades` always ramps pale-at-minimum to saturated-at-maximum, so reversing a gradient is a
matter of which variable you hand it — there is no direction flag. `log=True` fits on log10, for
quantities spanning decades where a linear fit would drop everything below the largest value into the
palest step or two; leave it **off** for a short ordinal range like 1–7, where log instead bunches the
high (darkest, most interesting) end together. `headroom` extends the fit floor below the minimum by a
*fraction of the observed span* — a fraction so it means the same on both scales — for the reason
`auroc_shades` anchors below chance: fitted exactly at the minimum, the smallest mark renders white
and vanishes. A legibility margin, not a data threshold.

**Step-05 pathogen colour comes from step 03, not from `SHARED_ORGANISM_COLORS`.** Both panels that
break a total down BY pathogen — `active_overlap_containment` and `exclusive_hit_model_rank_*` (with
its key) — take the hues `pathogen_activity_ratios` assigns, via
`plotting_colors.pathogen_activity_colors`. One pathogen therefore has one colour across steps 03 and
05. This is not cosmetic: the two palettes still **disagree** for most organisms, so leaving each panel
on its own palette would put the same pathogen in two colours inside one printed figure.

The worst case of that disagreement has since been fixed at the source. The positional palette used to
*swap* two organisms against the shared one — step 03 gave *S. aureus* the lime the shared palette
gives *E. coli*, and *E. coli* the periwinkle it gives *S. aureus*, which is the single most confusing
way two palettes can differ. `plotting_colors.PATHOGEN_HUE_SWAPS` now exchanges that one pair after the
positional assignment (via `swap_pathogen_hues`), so **E. coli is lime and S. aureus periwinkle in both
palettes**. The swap is applied inside both `PathogenActivityRatiosPlot` and
`pathogen_activity_colors`, so the panel and its step-05 consumers cannot drift apart. Every other
organism may still differ between the two palettes — do not mix them in one figure.

Those hues are **derived, not frozen**: step 03 colours *positionally*, `distinct_colors(n)` over
pathogens ranked by dataset count descending, so there is no fixed hue per pathogen to import and a
hardcoded copy would silently stop matching as soon as the ChEMBL curation adds or drops a dataset.
`plots_euopenscreen._step03_pathogen_colors` reproduces the ranking from step 03's own
`dataset_sizes.csv` — the module's **only** sibling-output dependency, and it degrades to a single hue
if step 03 has not run.

The **accent** level `(None, 0.78)` is used throughout, not step 03's pale `(0.62, 0.38)` dot fills.
The pale level was tried on `exclusive_hit_model_rank_*`'s bar segments — matching what step 03 puts
in its dots — and reverted. Same hue either way, but splitting the two step-05 panels across two
weights of it means a pathogen no longer reads as *one* colour across the step, and the pale level is
in any case too weak for `active_overlap_containment`'s ~2 mm pies. One level, everywhere.

**Marker areas are absolute and must be rescaled with the footprint.**
`submodel_auroc_summary`'s consensus dot went 90 → 40 pt² and its sub-model dots 16 → 8: at 45 mm the
7 rows share ~33 mm of axes (a 4.7 mm pitch), and a 90 pt² dot is 3.3 mm across — it crowded its
neighbours' rows and buried the sub-model spread. Same coupling as `SIZE_REF_AREA` in step 03.

**`active_overlap_containment` — a pie matrix at `SMALL_SQUARE`, cropping to 48.3 × 42.6 mm.** The
circle view of the same active sets as `active_overlap_jaccard`, which stays alongside it at `(4, 4)`.
The **wedge is containment |A ∩ B| / |A|**: the share of the ROW organism's actives also active against
the COLUMN organism. Directional, so the matrix is **not symmetric** and both triangles carry
information; a Jaccard heatmap cannot show this, because a union denominator lets the larger set
dominate (paeruginosa/abaumannii is 0.22 read either way, while containment is 93% one way and 23% the
other). **The diagonal needs no special case** — containment with itself is 1, so it draws as a full
circle. This is the one panel in the repo whose caption *must* state the direction, or the reader will
take it for symmetric.

*The size encoding was tried and dropped.* Circle area first carried the row's active count. That
works at 120 mm, but as a 45 mm paper panel the tick labels take ~60 % of the width and the pitch
falls to 2.4 mm: a linear area scale then put *S. aureus* (378) at 2.1 mm and *P. aeruginosa* (14) at
**0.40 mm** — a dot whose wedge cannot be read, losing the measurement in the two smallest rows to an
encoding that was only context. **All pies are now one size and the counts ride on the y tick labels**
(`A. baumannii (57)`), where they are exact rather than estimated off an area. General lesson: an area
scale spanning 27× does not survive a 45 mm footprint; put the quantity in text and keep the pitch for
the thing being measured.

*Radii are sized in points from the measured axes box, not in data units,* which forces the ordering:
ticks, limits, `set_aspect("equal")`, then `tight_layout()` and `canvas.draw()`, and only then the
marks. Getting it backwards sizes the pies against an axes box layout is still free to change.
`set_aspect("equal")` is what makes one pitch govern both axes.

*Grid on the MAJOR ticks, not cell boundaries* — a line through the circle centres points straight at
the label, where a boundary grid would box every mark in and set a second grid of edges against the
circles' own outlines. It is drawn under the marks (`set_axisbelow`), so it shows only in the gaps
between neighbours: that is what caps `MAX_DIAM_FRAC` at **0.76**, since at 0.86 the gap closes to
~0.3 mm and the grid vanishes.

*x tick labels rotate 90°, not the usual 45°.* Rotated labels on a shared baseline clear each other
only when `pitch × sin(angle)` beats their line height; at a 2.4 mm pitch, 45° gives 1.7 mm against a
~1.9 mm line and the genus names overlap. Same trade-off as `LABEL_ROTATION` in step 03, resolved the
other way because this pitch is far tighter.

*Its key is `active_overlap_containment_key`* `(0.5, 1.5)` = 45 × 15 mm — a wedge-scale strip only,
since there is one scale left to explain. Swatches go through `pie_scatter`, the same code path as the
marks, and are drawn **larger** than the panel's own pies: a key is read once up close to calibrate the
eye, so reproducing the mark at its true size would make 25 % vs 50 % as hard to see here as in the
matrix. They are also drawn **neutral** — colour identifies the row's organism in the panel, so a
swatch in one organism's hue would read as a claim about that organism.

#### Score-distribution boxes (`ScoreByHitClassPlot` and its two subclasses)

**No box is filled, and `inactive` is silver.** The class is the background the hit classes are
read against, not a finding: on the activity panel it is 97,162 of 97,590 compounds, so in cobalt
it carried the weight of a substantive category and dominated a panel that exists to show the 428.
Silver is the palette's reserved neutral. Its box is an outline with nothing behind it, which is
the right weight for a background class.

**`consensus_max_percentile_by_activity_dedup` is `BOX_NARROW` = 40 mm tall x 31.5 mm wide**
(crop 31.3 x 41.1). Two categories need no width. **1.05 cols is the floor while the tick labels
carry their class sizes** — `(n = 97,162)` is 10.4 mm on its own, so the two labels collide below
this (+1.5 mm clearance here, -0.3 mm at 0.9 cols). Moving the counts to the caption, as step 01's
`TaskMetricBoxPlot` does, would reach ~24 mm; they are kept because the 428-vs-97,162 imbalance is
the context for everything else in the panel.

**The active tick label is just "active".** The `(>=1)` qualifier came off — what it meant (a hit
in one or more of the 7 primary assays) belongs in the caption, and at this width it cost a line.

**Historically: fill only a box with nothing to reveal.** These panels used to pass `face=` for every class,
which hid all 428 active points behind their own box — only the ones spilling past the quartiles
were visible — and forced the median to INK, since a same-hue median cannot be seen on an opaque
body. `box_from_stats` already implements the house style; the panel was overriding it. Now the
`inactive` class keeps its fill (~10^5 compounds, genuinely never shipped per-molecule) and every
class that ships points is unfilled, so the swarm carries the distribution and the median takes the
class colour. Box width 0.5 and jitter 0.20 match the step-01 technical boxes
(`TaskMetricBoxPlot`); the `box_from_stats` defaults (0.34 / 0.12) are sized for far fewer points.
This affects all four box panels — `consensus_sum_by_hit_class`,
`consensus_max_by_activity`, `consensus_max_percentile_by_activity` and the `_dedup` twin.

**`consensus_max_*` panels print their AUROC**, carried in a new `auroc` column on every row of the
box-stats CSV. **It cannot be derived downstream, and the obvious shortcut is wrong.** AUROC is the
mean rank percentile of the positives only when the percentile is ranked over the *evaluated* set;
here the score is ranked over the whole scored library (106,290) and the AUROC is over the smaller
evaluated subset (97,590 after dropping unlabelled, then training-set, compounds). Those sets
differ non-randomly — dedup preferentially removes high scorers, because promiscuous hits are the
ones ChEMBL trained on — so the shortcut gives **0.7511 against a true 0.7574**. The inactive
scores are never shipped, so the column is the only route to the number. (The identity *does* hold
for `hit_exclusivity_events`, whose percentiles are ranked within their own evaluated set.)

**`consensus_max_*_roc` is the same numbers read the other way.** One panel per box panel
(`consensus_max_roc`, `consensus_max_percentile_roc`, `consensus_max_percentile_roc_dedup`),
drawn with the shared `plotting_utils.roc_panel` so it reads identically to the per-organism ROC
grid, in the active class's colour. Neither view replaces the other: a pair of boxes cannot answer
"how much of the library must I screen to recover half the hits", and a ROC cannot show that a
quarter of the actives score below the inactive median. The curve is exported thinned
(`_thin_curve`, 800 vertices) for the same reason the AUROC is exported — the ~10^5 inactive scores
never leave `run_consensus_max_by_activity`, so it cannot be rebuilt at figure time.

**Within-model normalisation is not just an axis change.** Raw max scores AUROC 0.7335; converting
each model to its within-model library percentile before taking the max gives 0.7943 (both `full`).
Only the final *re-ranking* is monotone — the per-model normalisation changes which model wins, and
is worth ~0.06 AUROC.

#### `hit_exclusivity_events` — the distribution behind the AUROC bars

A companion to `hit_exclusivity_auroc`, at **`EVENT_SQUARE` = 40 × 40 mm** (cells are 3 cm, so
mm/30 → 4/3), showing what the bars summarise. Smaller than the ~49 mm `EXCLUSIVITY_BARS` of the panel
it accompanies because ~13 % of its width goes to the AUROC column outside the axes; measured crop
is 40.4 × 41.0 mm. Two lanes per organism — shared hits (periwinkle, upper) and exclusive hits (amber,
lower) — drawn as a **matplotlib event plot**: one vertical line per hit, at that compound's
percentile in its model's ranking of the whole evaluated set. Mass to the right = hits upranked.
Organisms are ordered by the same rule as the AUROC panel (best on top) so the two read together.

**The right-hand axis is the same quantity as the ticks, not a second variable.** AUROC is the
mean rank percentile of the positives — exactly `(auroc * n_neg + (n_pos + 1) / 2) / n`, which at
this prevalence is the mean of the plotted ticks to four decimals. The printed value is therefore
where its lane's centre of mass sits. `eval_euopenscreen._check_percentiles_match_auroc` asserts
the identity when the CSV is written, so a mis-joined export fails loudly instead of drawing a
plausible-looking wrong distribution.

**Lines, never a density.** *P. aeruginosa* exclusive is one compound and *E. faecium* four; a
smoothed curve over one point draws a shape that is not in the data. One tick per compound cannot
overstate the evidence, and it removes any bandwidth or bin-width choice.

**The AUROC pair shares one line per organism, `shared / exclusive`, and keeps both colours.**
A tick label is one text object and therefore one colour, so this is drawn as three free texts —
value, slash, value — anchored right / centre / left, which stays aligned whatever the glyph
widths. `AUROC` is a header above the column, not a rotated axis label: rotated it is nearly as
tall as the whole 40 mm panel, and a vertical label beside a column of numbers reads as an axis
that is not there.

**The column's offset is MEASURED, in points, not guessed.** `_text_width_pt` renders the widest
left-hand value and the slash is placed at `clearance + that width + gap`, so the only free
parameter is a 3 pt clearance. Two earlier versions got this wrong in opposite directions: an
axes-fraction offset laid the numbers on top of the plot, and the hand-set point offset that
replaced it left **11.4 pt (~4 mm) of dead space** because it had to be generous enough to clear a
width nobody had measured. On a 40 mm panel that is the dominant layout cost — removing it took the
axes from **8.3 mm to 12.0 mm wide, +45 % of plotting area** at the same footprint. Points rather
than fractions so the clearance survives a change of footprint.

**It carries no genus labels, and that is tied to its row order.** `WIDE_EVENTS = (1.3626, 1.3293)`
is tuned by the shared-axes-height method above so its axes lands at **25.71 mm against
euos_overlap's 25.81 mm**, and it adopts that panel's y limits **(-0.6, n - 0.4)** rather than a
symmetric (-0.5, n - 0.5). The limits matter as much as the height: the spans differ (7.2 vs 7.0), so
equal heights with unequal spans drift the row centres ~0.4 mm apart by the top of the panel. The two
are placed side by side against **one** set of names, which frees the ~14 mm the labels used to take.

**Active counts are compound-level and reconcile with the exclusivity split (2026-08-05).**
`load_euos_primary` collapses duplicate SMILES in the primary-assay labels, keeping the highest bin
("Active prevails", the same rule upstream uses for `02_merged`). Before that the labels were
measurement-level — one row per assay well — so structures the library registered twice were counted
twice, and the own-assay active count sat up to 1 above `exclusive + shared` for 4 of 7 organisms.
The two now agree exactly in both `full/` and `deduplicated/`, which matters here because
`hit_exclusivity_auroc` and `hit_exclusivity_events` are read against the own-assay panels. Row order
is unaffected — `n_active` moved by ≤1 and crossed no boundary — so none of the six matched footprints
needed re-tuning.

**One helper, `overlap_row_order`, defines the row order for EVERY panel with a pathogen y axis** —
`n_active` ascending from the leakage report, bottom row first. `_apply_row_order` adapts it per panel
(unknown organisms are appended, never dropped) and `save_euopenscreen_figures` resolves it once and
threads it in as `row_order=`. Panels on it: `euos_overlap`, `hit_exclusivity_events`,
`euos_shared_enrichment`, `hit_exclusivity_auroc`, `specificity_index`, `submodel_auroc_summary`.
Sorting independently in each would work until the active counts moved and then silently relabel every
row of one panel. Suppressing the labels and adopting the order are also ONE decision in the code:
without a `leak_df` the panel keeps its own metric order AND its own labels, so it can never end up
with unlabelled rows in an unrelated sequence.

**The trade-off is explicit:** these panels used to sort by their own metric, best on top, which is
easier to rank *within* one panel. A shared order makes any cross-panel comparison a glance instead of
a lookup, and the metric is still on the axis, so nothing is hidden. The matrix panels
(`active_overlap_jaccard`, `active_overlap_containment`, `cross_organism_heatmap`) are **not** on this
order — they carry pathogens on *both* axes, so a reorder has to move rows and columns together to keep
the diagonal meaningful, and they are already mutually consistent on `SHARED_ORGANISMS`.

**`*` marks an AUROC computed over `_SMALL_N` = 10 hits or fewer** (four of fourteen; two are
n = 1 and n = 4). A flag, not a filter — nothing is dropped. **The caption must say what it means.**
It also has to be part of the width measurement: the pad is measured on the FORMATTED string, so a
bare-number measurement would let a marked value creep back over the plot.

**Tick width is 0.4 pt.** At 0.6 the dense lanes (S. aureus exclusive, 176 hits in ~12 mm) closed
into a solid block and hid the clustering the panel exists to show.

**No chance line at 0.5.** On the AUROC panel a 0.5 rule marks a threshold for a single bar; here
every lane is a spread that straddles it anyway, so the rule only draws a stripe through all
fourteen lanes.

**Two traps, both of which produced a silently wrong panel during development:**
- **Know which bar helper inverts the y axis, because only one of them does.** `hbar` calls
  `invert_yaxis()`, so its **first element lands on top**; `grouped_hbar` and bare
  `ax.barh`/`scatter` on `arange(n)` leave the axis alone, so their **first element is at the
  bottom**. `MetricByOrganismPlot` uses `grouped_hbar` (bottom-first); `specificity_bars` goes
  through `hbar` (top-first) and is the one panel that must be handed the shared row order
  **reversed** — `SpecificityIndexPlot` does exactly that. Get it wrong and the panel is a mirror of
  the ones it is meant to be read against, with nothing erroring. This panel places its own lanes, so
  it needs `set_ylim(-0.5, n-0.5)` to put index 0 at the bottom.
- **Verify row order by screen position, not tick index.** `get_yticklabels()` returns labels in
  ascending *data* order, which on an inverted axis is top-to-bottom — so a naive readout reports a
  correct panel as reversed. Sort the labels by `get_window_extent().y0` instead.
- **Lane geometry is computed once, in `_lane_y`.** The ticks and the right-hand labels are
  positioned independently; when that arithmetic was written out twice, fixing the sign in one
  place swapped every printed AUROC onto the wrong lane while the figure still looked entirely
  reasonable. Never inline the offset.

**The twin axis must have its grid switched off.** It inherits the style's grid and its ticks are
at the 14 lane positions, so left on it rules a grey line through every lane and the tick marks
read as data.

**Data.** `05_hit_exclusivity_percentiles.csv` (732 rows for dedup), one row per active:
`pathogen, code, eosid, subset, set, inchikey, smiles, percentile`. Written by
`eval_euopenscreen.run_euopenscreen` alongside the AUROC record it unpacks, and split into
`full/` and `deduplicated/` by the same `set` column as the other long-form tables.

