# Scripts

Numbered scripts run in order. Figure and plotting conventions (sizing, formats, per-figure
layouts) live in [`docs/figure_conventions.md`](../docs/figure_conventions.md).

## 00_download_data.py
Stages all input data in four sections: companion repos / eosvc (EU OpenScreen tasks, CoAdd data,
ChEMBL model reports and curation summaries), public GitHub files (Ersilia reference library,
DrugBank), Airtable model metadata, and Isaura precalc predictions for Ready annotation models.
Skip-if-exists, so it is safe to re-run. `--skip-isaura` skips the slow Isaura section; `--eosvc`
pulls Section 1 from eosvc instead of the companion repos.
Section 4 also fetches `eos1klk` (2D projector, Task=Representation/Projection) explicitly — its
Task isn't "Annotation", so the automatic loop above never sees it — into
`data/processed/eos1klk_projection/eos1klk_v1.csv`, feeding step 09.
**Key decision:** ChEMBL model performance metrics are recomputed downstream from the per-fold
`09_reports`; the staged `10_reports/` CSVs are kept only for fields that cannot be reconstructed
(quality weights, decision cutoffs, discarded-model reasons).

## 01_ersilia_metadata.py
Analyses the Airtable model metadata (Ready models only, n = 208) and renders 13 panels — the
task/subtask breakdown, the source-type and output composition, biomedical area, target organism,
three composition donuts, the container-metrics row and the pathogen circle-treemap.

> **Note (2026-08-04):** this section was rewritten from the code after an accidental deletion of the
> previous text. It describes what the code does today; any decision recorded only in the old prose
> is not here. Review it before trusting it as the decision log.

**Shown vs counted:** value counts for every counted field are written to `*_counts.csv`, but not all
are plotted. **Counted only, no panel:** Publication Type, Tag and the ungrouped License field (the
licence composition reaches the figure as the reuse-class donut, not per identifier). Target Organism
is capped at the **top 10** categories in its panel; Biomedical Area is drawn as **four groups over
Activity prediction models only**, not as its 16 raw areas. Full counts remain in the CSVs either way.

**Stacked panels:** Source Type and Output are each drawn as one stacked bar panel segmented by
Subtask, so a single panel carries the joint distribution instead of two. Segments use the subtask
palette (shades of the parent task's hue), never a palette of their own — the field is already encoded
by the bar position. Cross-tabs are written as `<field>_by_subtask_counts.csv`.

**Task hues:** Annotation = crimson, Representation = amber, Sampling = lime (`TASK_HUES`, the one
place to change them). `SUBTASK_COLORS` derives from the same mapping, so a task and its subtasks can
never drift onto different hues. `SOURCE_TYPE_COLORS` now survives only for the pathogen panel's dots
— the two sets never carry colour in the same panel.

**Task/subtask alternatives** — two ways to draw the same 208 models, on the same colours:
- `task_subtask` — bars, one per subtask, grouped by task and ordered by count within it. Best for
  reading counts. No legend: each bar is named on the axis next to its own colour, which makes this
  panel the subtask colour key.
- `task_subtask_waffle` — one square per model, **16 × 13 = 208 exactly** (16 is one of the few column
  counts that divides 208 without a ragged last row; the script prints `blank`, currently 0). Shows n
  rather than stating it; reading 52 vs 39 means counting. Carries its own legend with counts, so it
  is self-contained and can serve as the key too. That legend band is also what squares the panel —
  the grid alone crops landscape (1.23 at 60 mm), and at quarter width the legend labels have to be
  abbreviated (`_LEGEND_ABBREV`) or the key sets the crop and the squares fall to 1.99 mm.

**Legends:** the two `*_by_subtask` stacks carry none — no room at ~52 × 25 mm each, where a 6-entry
key would be several times taller than the plot, and they are meant to sit beside a key panel. **So if either goes into the figure, `task_subtask` or
`task_subtask_waffle` must travel with it.** The waffle and each of the three donuts carry their own.

**Biomedical Area is four groups, not 16 areas** (`BIOAREA_GROUP` in `src/default.py`, signed off
2026-08-02), over **Activity prediction models only**: Antimicrobial 50, ADMET 29, Antiviral 7, Other
8. Two things the counting has to get right, both printed on every run:
- Biomedical Area is **multi-value**, so these are counts of *distinct models* per group, not of area
  assignments. Grouping absorbs most of the multiplicity.
- Two models (`eos2zmb`, `eos7kpb`) carry areas in two different groups and are **counted in both**, so
  the bars sum to 94 over 92 models. That is the metadata's own claim, left unresolved rather than
  silently reassigned. The donut shows the 92 in its hole, so a caption must explain the discrepancy.
An area with no `BIOAREA_GROUP` entry **raises** rather than vanishing from the figure. Per-area
membership goes to `biomedical_area_groups.csv`.

**Biomedical Area / Target Organism colours:** every *named* area and organism consists purely of
Annotation models (`Any` is where the Representation and Sampling models sit), so named bars take the
**Annotation hue** — no separate palette — and `Any` takes **silver**, the reserved catch-all neutral
(`catchall_colors`). The Biomedical Area strip is one flat full-strength crimson instead, since every
bar in it is already Annotation / Activity prediction and colour distinguishes nothing within the
panel. Neither panel carries a key, so the caption must explain the colours.

**Pathogen circles.** One circle per priority pathogen, packed with `circlify`, holding one dot per
model; dot **area** encodes that model's training-set size on a single global scale, dot colour encodes
Source Type. A pathogen's circle is only an **enclosure** — its area encodes nothing and tracks model
count whatever the dots hold. Which models appear is driven by `config/model_training_sizes.csv` (the
curated model-to-pathogen mapping, rows with too-generic organism assignments removed by hand), **not**
by the Airtable Target Organism field, so it shows fewer models than the Target Organism bar panel.

**Known upstream mis-annotation:** `eos93h2` (`image-mol-gpcr`) is recorded in Airtable as
Representation / Featurization with Target Organism = *Homo sapiens*, but it should be an Annotation
model — being fixed at source. It is the only thing that breaks the "named organisms are
Annotation-only" rule, and it is already drawn as an Annotation row, so these two panels need no
change when the fix lands. When it does, **Task and Subtask must be changed together** — see the
warning in `docs/figure_conventions.md`, since a Task-only edit would leave the subtask-coloured
panels disagreeing with the task-coloured ones without erroring.

**Technical box row.** Three **horizontal** panels per Task on a log x axis, sharing one task axis and
occupying **4/6 of the page width (120 mm) × 30 mm** between them: `runtime_100`, `image_size` and
`output_dimension`. Only the leftmost draws the task tick labels, so the other two **cannot be placed
on their own**; `runtime_100` is wider by exactly its 14.2 mm label column so all three metric axes
come out the same ~28 mm length. Per-task quartiles go to `technical_metrics_summary.csv`.

**`n` is not in the tick labels** — coverage differs per metric (runtime 129/57/10, the other two
131/58/19), so a shared label set cannot carry it. **A caption must state that the runtime box for
Sampling rests on 10 of 19 models.**

**`output_dimension` is decade-binned circles, not a box.** The column is heavily tied (68 of 131
Annotation models output a single value; 100 of them fall in the 1-9 bin, which is the circle the
panel draws), so a swarm just overplots. One circle per (task, decade),
area ∝ model count, placed at the bin's geometric centre so it reads as an interval rather than as an
exact value. No size key fits at 34 mm — exact counts are in `output_dimension_bins.csv`. It is a
drop-in replacement for the box version: same axis, footprint and orientation.

**Box style (repo-wide).** Distribution boxes are **unfilled**, outlined in the category colour at
0.5 pt — the house line weight — with the median in the same colour drawn heavier. matplotlib's
boxplot default is 1.0 pt, twice stylia's `lines.linewidth`, which made every box the heaviest mark in
its panel and buried the swarm underneath. Pass `face=` only where the box has marks drawn inside it
or no swarm at all; on a filled box the median reverts to `INK` so it survives an opaque body.
**Scripts 02, 03 and 05 need a re-run to pick this up.**

**Three donuts, one family** (`DonutPlot`), **25 mm wide each** (75 mm for the set): ring, total in the
hole, key beneath with names and counts. Labels cannot go around the ring at this width — it is
geometrically impossible, not just tight ("Non-commercial 1 (<1%)" alone is 20 mm of text) — and the
legend is also the only way a 1.7° wedge gets named at all. **Legend rows carry the count but not the
share**, because the legend has to stay narrower than the ring or it collapses it, and the share is
what the ring already encodes. **The ring is pinned to a fixed 19.6 mm in all three** (`DonutPlot.pin_ring`) — otherwise
`tight_layout` sizes it around the legend and a longer label silently shrinks the ring; crop widths are
an identical 24.89 mm. **Heights still differ** (the two-entry architecture donut is 30.45 mm tall, the
others 35.00 mm), so align the set on the rings or the left edges.

- `license_class_donut` — Permissive 104 turquoise, Copyleft 76 periwinkle, Not recorded 27 silver,
  Non-commercial 1 **fuchsia** (the one use of that hue in the repo, because a 1-of-208 wedge in any
  calmer colour is invisible). Four **reuse classes, not the ten licences**: four licences have exactly
  one model each, a 1.7° wedge that cannot be seen or labelled, so per-licence detail stays in
  `license_grouped_counts.csv` (and the ungrouped `license_counts.csv` beside it, which keeps the
  `-or-later` / `-only` distinction the grouped file collapses).
- `docker_architecture` — **AMD only 79 tangerine, AMD + ARM 129 cobalt**, ordered base-capability
  first. Moved off turquoise/periwinkle so no hue means two things across the set; tangerine is a plain
  categorical hue here, **not** a warning about x86-only builds.
- `biomedical_area_donut` — the alternative to the bar strip, on the same four groups. **One hue, four
  patterns**: every model in it is an Annotation model, so colour would encode nothing; groups separate
  by hatch instead (solid → diagonal → dots → cross-hatch for the catch-all). **The hole reads 92 while
  the legend rows add to 94** — two models carry areas in two groups. That is real; a caption should say
  so.

**Runtime batch size: 100 molecules** (`RUNTIME_BATCH` in `src/default.py`; the five
`Computational Performance` columns are 1/10/100/1,000/10,000 molecules). A `-1` in those columns means
the benchmark was **never run**, not zero — those models are skipped, never imputed. 100 is chosen
because it is the largest batch where generative models still have data: coverage is **129/131, 57/58,
10/19** at 100 molecules but **123/131, 55/58, 0/19** at 1,000. The cost of that choice: at 100
molecules the median runtime is no higher than at 1 molecule (CP3/CP1 ratio 0.86), so for annotation
and representation models the number is container startup, **not throughput — do not divide it by
100**. What it buys is the Sampling box: median **481 s** against 34 s (Annotation) and 29 s
(Representation), i.e. generative models are ~14× slower. That box rests on 10 of 19 models, so it is
indicative only; every tick label carries its `n` and the script prints coverage on each run.

The architecture donut is a **snapshot, not a trend** (45% dual-arch among 2021 models vs 77% among
2026 — the 62% is accumulated stock, not current practice), so a caption needs the metadata snapshot
date.

**Panel sizes:** `task_subtask` is 60 × 60 mm (`cells=(2,2)`); the waffle and the two box panels are
45 × 45 mm; the two ten-category fields are 45.75 × 45.75 mm (a quarter of the page width, two side by
side); the three donuts are 25 mm wide; `pathogen_circles` crops to **63 × 76 mm** — on a
183 × 170 mm page. The
`*_by_subtask` stacks are sized as a **pair to a page budget of 52 × 50 mm** (2026-08-05), with
`Output` (4 bars) above `Source Type` (3 bars). Two properties are solved for and printed on every
run: the pair's saved pages fill the budget without exceeding it (**51.96 × 49.96 mm**), and both
panels draw **bars of the same thickness** (3.03 mm, equal to 0.0%) despite holding different numbers
of them. Re-sizing the pair means editing those two budget numbers only — the footprints, the height
split and the bar thickness are all solved from them. Only `Source Type` carries the `"Number of models"` axis title, so the unit states its
quantity once, beneath the block.

**Cutoff:** equal bar thickness costs a raised *layout* dpi (600, not matplotlib's 100). The layout
grid is one canvas pixel, so at the default dpi the nearest achievable split of the height leaves the
bars 1.8% apart — a floor no constant can beat. It affects sizing precision only, not output
resolution, and is set on that one plot class so the figure's other calibrated panels are untouched.
Re-measure the three sizing constants with `tools/probe_stack_geometry.py` if the fonts or labels
change. Full derivation, plus a ~0.35 mm axis-registration residual between the two panels, in
`docs/figure_conventions.md`.

**`pathogen_activity_ratios`.** One dot per modelled dataset: active fraction on the y axis where it
can be resolved, dataset size as **dot area** (affine in √size, keyed at 100/1,000/10,000/100,000 —
not area-proportional, since the range is 3,500×), one colour per pathogen. Colour is redundant with
the x axis on purpose: only 9 substantive hues exist, so `distinct_colors` reuses them as tints beyond
that, and the labelled axis is the real key. The pattern to read: the largest datasets sit at
near-zero active fraction while small ones are balanced or active-heavy.
**Added negatives are excluded from both the size and the ratio** — size is
`n_compounds − n_added_negatives − n_added_decoys`, which reproduces the curation pipeline's own
`n_mol_after`/`ar_after` exactly. Only 54 of 193 models had negatives added, so the medians barely
move, but for those 54 it matters a lot (`mtuberculosis/DR_0012`: 2450 compounds at 0.50 active →
1411 at 0.87). Decoys are zero for every model. 54 of 193 datasets are majority-active once added
negatives come out. Derived table: `dataset_sizes.csv`.
A short bar per column marks the pathogen's **unweighted mean** active
fraction (size-weighted would be much lower, since the huge datasets are near 0% active). Size key is
nested circles with leader-lined labels. **120 × 60 mm** with column widths proportional to
`n_datasets ** 0.25`: the three single-dataset pathogens need no spread, so *P. falciparum* still gets
10.9 mm in a panel with **44% of the original area**. Tick labels are single-line; counts are in the
CSV. Dot sizes scale with the footprint (`SIZE_REF_AREA`) — rescale both together.

**Area encoding:** a dot's *area* encodes its model's training-set size
as `n ** AREA_EXPONENT` with **`AREA_EXPONENT = 0.5`**, i.e. area ∝ √n, set in
`src/plots_metadata.py`. Sizes span 96 to 5,000,000 compounds, so a linear encoding buries the
small datasets. The exponent is the one knob: 1.0 is true data volume, →0 is pure model count.

**Why not `1 + ln(n)`** (the earlier choice): every model then contributed at least ~5.5 to its
pathogen's total regardless of size, so a pathogen's total tracked how *many* models it held
rather than how much data — *E. coli* (12 models, 2.4M compounds) came out **larger** than
*P. falciparum* (7 models, 6.0M). A power transform makes a small model contribute small area, so
count stops inflating groups. At 0.5 the ordering follows data volume (*Pf* 33.6% vs *Ec* 21.5%
of the training compounds) and the smallest mark stays ~9x above the readable floor. Note the count
correlation cannot be driven to zero — pathogens with more models genuinely do hold more data
(r = 0.62 even for strictly linear areas); what the transform removes is count *mechanically*
inflating size.

**Circle panel enclosures:** the grey circle is an enclosure, sized to whatever radius circlify
needs plus a constant ring, so its dots fill only 43–75% of it. Under √n it nonetheless orders
correctly (*Pf* r=77.2 > *Ec* r=62.9; enclosure area correlates 0.945 with data volume, 0.778
with model count) — but it is still not an exact quantity, so read model count by **counting
dots**, not from circle size.

**Size legend (circle panel):** decade references (100 / 10,000 / 1,000,000 compounds) drawn as
*nested* circles sharing a bottom tangent, in a band reserved inside the axis limits below the
packing. Drawn in **data coordinates**, so they are exactly to scale with the dots by
construction. Nested rather than in a row because under √n the keys span a 10x radius range and a
row that wide takes over the panel. The band is reserved *before* labels are placed — adding it
afterwards rescales the axis and every measured text extent then understates the rendered width.

**Genus labels** use real rendered text extents (`_text_extents`), not a character-count
estimate, which is what stopped neighbouring labels touching. Placement is: inside the circle if
it fits, else adjacent to its own rim, else parked in free space with a **leader line**. Forcing
an over-wide label inside a small circle is deliberately the last resort — it is what previously
put *Enterobacter* and *Campylobacter* on top of each other.

**Pathogen treemap — which models appear:** driven by `config/model_training_sizes.csv`, *not* by
the Airtable Target Organism field. That file is the hand-curated model-to-pathogen mapping (one
row per dot); three models whose organism annotation was too generic to be meaningful (eos2gth,
eos2xeq, eos74km) were removed from it by hand and so do not appear. **The treemap therefore shows
fewer models than the Target Organism bar panel**, which still counts every Airtable annotation —
e.g. *E. coli* is 12 dots here against 13 bars there. Training sizes for the 15 shipped
ChEMBL/PubChem pathogen models are unique-compound counts over each model's retained endpoints;
eos3dys is counted per organism over its CoAdd strain files; the rest were curated by hand from
their publications.

### Cumulative models over time (`models_over_time_by_task.csv`)

No panel — this is a **summary CSV only**, kept because `scripts/01b_community_stats.py` reads it for
the Models track of the hub timeline. It is the only consumer of Airtable's `Incorporation Date` column;
every panel in this script is a snapshot. (Four time-axis panels — cumulative models by Task and by
Source Type, cumulative contributors, and the publication→incorporation gap — were dropped from the
figure on 2026-08-04 along with `src/plots_growth.py`.)

**Denominator: 215 models, not 208.** The series runs on the **unfiltered** metadata, unlike every
panel here, which uses `Status == "Ready"` (208). A model currently in maintenance was still
incorporated on its date and still counts towards how the hub grew, so **a caption using this series
must give 215, not the 208 the rest of the figure states.**

**Excluded: 4 models with no `Incorporation Date`.** The date is the x axis, so there is nowhere to
put them. They are excluded, never imputed, and the script prints the count each run.

## 01b_community_stats.py
Draws the whole history of the hub as **four thin tracks on one shared year axis**: models
(cumulative), people (cumulative), commits and issues (per month). This is the repo's only
time-series figure. Outputs to `output/01b_community_stats/`, plus `01b_timeline_series.csv`
carrying every series on one month index so the figure is reproducible from a single file, and
`01b_snapshot.txt` with the GitHub collection date.

**Numbered 01b because it reads step 01's output.** The Models track is
`output/01_models_metadata/models_over_time_by_task.csv`; the other three tracks are aggregated
here from `data/raw/github_stats/`. It fails with a message naming the script to run if 01 has not
been run. This script is the merge of the former steps 08 (community aggregation) and 09 (the
timeline figure), reduced to only what the timeline needs — the eight community panels, the
contributor region map and the workflow-runs series were dropped on 2026-08-04.

**Hub-scoped, not org-wide: 7 core + 240 model = 247 repos.** Core is `ersilia`, `ersilia-pack`,
`ersilia-pack-utils`, `ersilia-maintenance`, `ersilia-model-workflows`, `eos-template`,
`ersilia-apptainer`. Models are matched by `^eos[0-9a-z]{4}$` — **an `eos` prefix is not enough**,
since it also catches `eosvc`, `eosbench`, `eosframes`, `eosquality`, `eosdev` and `eos-demo`,
which are not models. `eos-template` counts as core. The `ersilia-os` org has ~421 repos in total
including websites, grant repos and capstones; **none of those are in these three tracks**, so a
caption must not describe them as organisation-wide.

**Only the repo *list* is used, not the core/model split.** `repo_set.csv` records a group per
repo, but every track here is a total, so the group is never read. (The former step 08 stacked its
commit and issue panels by group; those panels are gone, along with `REPO_GROUP_COLORS`.)

**Stocks and flows share one figure, and the caption must say so.** Models and People are
cumulative — the height *is* the size of the hub. Commits and Issues are per-month counts — the
height is a rate. **Comparing heights across tracks is meaningless.** What the figure is for is
comparing *positions*: the Outreachy cohorts, the model-incorporation bursts and the commit spikes
all line up vertically, which is the entire point of the shared axis.

**Headline numbers** (2026-08-02 snapshot, as plotted): models reach **215**, commit authors
**107**, **15,828 commits** and **1,469 issues** in total.

**Denominator: 215 models, not 208.** The Models series runs on the **unfiltered** metadata,
unlike every panel in step 01, which uses `Status == "Ready"` (208). A model in maintenance was
still incorporated on its date and still counts towards how the hub grew.

**"People" means distinct commit authors (107), not all participants (335).** The issue/PR-author
series was dropped from the figure at the author's request; it is still in
`01b_timeline_series.csv` and in `01b_snapshot.txt`. **Captions should not say "contributors"
without qualifying it** — three times as many people take part as write code.

**Issues only, PRs excluded.** `org_participation.csv` carries both; a PR is a different act and
the track is labelled Issues. Bots are excluded via `default.GITHUB_BOT_ACCOUNTS` — an explicit
list, because GitHub types them all as "User" and a substring match on `bot` catches real people.

**Commit months come from week starts.** `stats/contributors` reports whole weeks, so a week is
assigned to the month its Monday falls in and commits near a boundary can land one month early.
Immaterial monthly; splitting a week across two months would invent precision the endpoint does
not have. Weekly rows carry a commit *count* and are expanded to one row per commit before
binning, or the monthly figure would count weeks.

**One file, not four, because of alignment.** All four axes render at identical left/right
figure-fraction positions on one `xlim`. Separately-saved panels cannot do this — the tight crop
sizes each file's tick-label column to its own content, which is what leaves the two subtask
stacks in step 01 0.25 mm out of register (see `docs/figure_conventions.md`).

**Size: 180 × 31.5 mm nominal, 32.3 mm as cropped** — the cropped number is the one tuned, since
the outermost tick labels sit outside the axes box. Each track gets **4.4 mm** of drawing height.
That budget is what forces hand-placed axes, a measured inter-track gap, two y ticks per track,
horizontal rather than rotated track labels, x ticks on the bottom track only, and summing away
any category breakdown. The figure sits near the height at which matplotlib's `tight_layout`
refuses to run; the script warns if a future change crosses that line.

**Vertical gridlines only.** The year rules are what the shared axis is for. Horizontal rules
cross them every couple of mm at this track height and read as texture rather than as a scale.

**No workflow-runs track.** GitHub deletes run records after ~13-14 months, so that series would
start only in 2025-06 and be empty across three quarters of a six-year axis. Verified: `eos4e40`
was created 2020-11 but its oldest surviving run was 2025-07-08 — the floor tracked retention, not
repo age. The collector no longer fetches them.

**Partial month trimmed from the plot, kept in the data.** The collection date fell 2 days into
2026-08; drawn as a full month it is a drop to near-zero at the right edge of three tracks, which
reads as the project stopping. The figure stops at the last complete month. **Nothing is dropped**
— `01b_timeline_series.csv` keeps the partial month and the script prints what was held back
(0 commits, 1 issue). `01b_snapshot.txt` gives all-time *and* as-plotted totals, because they
differ (336 vs 335 issue/PR authors, 1,470 vs 1,469 issues).

**NaN handling.** The cumulative source series start in different months, so the joined frame has
leading NaNs in the people columns. Those are filled with 0 and made monotonic with `cummax()`,
never interpolated: they are integer counts of people, and a straight line between two months
would invent arrivals. A **trailing** NaN is left as NaN and means something else — that series
ends before the snapshot month. `models_cumulative` is blank for the partial month because step
01's series stops at the last month carrying an Incorporation Date, which is not the same claim as
"no models were added".

**EXPLORATORY — outside the `00_download_data.py` convention.** The GitHub input is collected by
the standalone `tools/fetch_github_stats.py` rather than by step 00, because this figure may not
survive review. If it is kept, fold the collector into a "Section 5 — GitHub" of step 00 and
delete the tool.

**Snapshot-dependent.** Unlike every other step, re-running this on a different day gives
different numbers. The collection date is in `data/raw/github_stats/snapshot.json` and copied to
`output/01b_community_stats/01b_snapshot.txt`. **Captions must carry it.**

## 02_chembl_data_curation.py
Reproduces the ChEMBL data-curation story (upstream step 27 of `chembl-antimicrobial-tasks`) as 15
individual panels, rebuilt **entirely from the staged summary CSVs** in `data/raw/chembl_curation/`
— no molecule-level data needed. Outputs to `output/02_chembl_data_curation/`.
**Snapshot:** ChEMBL `chembl_36` (read from `general/27_chembl_space.json`).
**Caveat:** the staged summaries' numeric `discard_step` codes do not match the upstream
`STEP_LABELS` constant, so the curation-outcome figures key on the `discard_reason` **text**, not the
step number (otherwise the largest bucket, `≤5 molecules`, is mislabelled and dropped).
**Not reproduced (by design):** chemical-space overlap heatmaps, binarisation agreement, embeddings,
and molecule-level panels — their values are not written to the summary CSVs.

## 03_chembl_models_performance.py
For each pathogen in `data/raw/chembl_model_reports/`, loads the 5-fold cross-validation reports,
computes mean ± std AUROC per model, and renders a per-model ROC-curve grid and paired rank
boxplots (plot classes in `src/plots_chembl_performance.py`). Outputs one summary CSV per pathogen
to `output/03_chembl_models_performance/`, plus PNG + vector PDF panels and a `figure_cells.json`
footprint manifest to its `individual_plots/` subfolder.
**`individual_plots/` are intermediate results.** 30 per-pathogen panels, far too many for the
paper; they exist to inspect every model. The condensed cross-pathogen figures that go in the paper
land at the top level of the output dir.
**Cross-pathogen consensus figure** (`pathogen_consensus_auroc`, top level): one row per pathogen —
a small **periwinkle** dot per retained model's mean CV AUROC, a grey min–max guide line at the house
0.5 pt, and a translucent red dot at the **mean** (the representative per-pathogen value), best on top.
The mean dot's
**area encodes the number of unique molecules** in the pathogen's ChEMBL dataset (`n_cleaned_inchikeys`
from the staged step-02 `chembl_curation/general/27_chembl_coverage.csv`), area-scaled between a min
and max marker so small datasets stay visible. **No size key fits at this footprint** — the range
(currently 170 → 500,199 molecules) is printed on every run and belongs in the caption, or the size
encoding is decorative. The two-entry colour key sits inside the axes, upper left: rows run
worst-to-best upward so the top rows' marks are pushed right, leaving that corner the one region the
data does not reach. Dashed line = the step-10 retention floor (mean AUROC 0.7).
**This is a summary of per-model CV AUROCs, not a per-compound consensus prediction** — the CV pools
share no held-out set, so a true ensemble AUROC is not computable here.
**Scope: step 09, not step 10.** The staged per-pathogen folders are the pipeline's `09_reports`
(**196 trained models**), a superset of the **193** retained by step 10. Three models were trained
and reported and *then* discarded for mean AUROC < 0.7 — `calbicans/588506` (0.578),
`hpylori/SP_catchall` (0.656), `pfalciparum/743093_merged2` (0.533). **All 196 are plotted, with no
visual distinction**; the verdict is instead recorded per model in the summary CSVs as `retained`
(bool) + `discard_reason`, joined from `10_reports/`. The other three step-10 discards are
`untrainable: min class 0 < 5 folds` and never produced a step-09 report, so they cannot appear.
Note *H. pylori* has only 2 step-09 models, one of which is discarded. The script prints a
196 = 193 + 3 reconciliation and warns if the staged `10_reports/` drifts out of sync with
`09_reports/`.
**Join key:** the report **file stem**, forced to `str` on both sides — 26 models are named with
digits only (e.g. `1242`), so joining on the in-file `model_name` column would type them as `int64`
and silently miss them.
**ROC grid: 6 columns.** Full 180 mm page width, so each ROC panel is a true 3 cm square cell and
the model-rich pathogens stay within a supplementary-page height (*P. falciparum* 52 models → 9
rows / 27 cm; *M. tuberculosis* 34 → 6 rows). No models are dropped from any figure.
**Rank boxplots: 4 models per 3 cm cell** (minimum 2 cells tall), which keeps ~7.5 mm of vertical
space per model — enough for the actives/inactives box pair to stay legible at every model count.
**AUROC colour scale fitted to [0.35, 1.0].** Curves are shaded by mean AUROC on a cobalt fading
colormap; the low anchor sits below the 0.5 chance level because fitting at exactly 0.5 renders a
chance-level curve white and invisible.

## 04_ersilia_predictions.sh
First live **Ersilia CLI** step (all earlier predictions came from the Isaura precalc cache). Runs
`ersilia fetch → serve → run → close` per model, writing one bare `{eosid}.csv` per model per library
to `output/04_ersilia_predictions/{euopenscreen,coadd}/` (consolidated inputs in `inputs/`).
**Library-major ordering:** all EU OpenScreen predictions (every model) run first, then all CoAdd —
EU OpenScreen has priority. Idempotent (skip-if-exists), with an append-only `_failures.log`.
**Models (16):** the 15 in `config/pathogens_of_interest.csv` plus the CoAdd model `eos3dys`; every
model predicts both libraries.
**Libraries:** EU OpenScreen reuses `data/raw/euopenscreen_data/02_merged/02_only_smiles.csv`
(106,290 compounds); CoAdd uses the `std_smiles` column of `data/raw/coadd_data/00_smiles_info.csv`,
dropping 199 rows with empty `std_smiles` (failed upstream standardization) and deduplicating →
100,005. Both are ~100k compounds. Blank SMILES are filtered from both inputs so Ersilia is never
fed an empty row.
**Environment:** requires `ersilia` in a conda env (`ERSILIA_ENV`, default `ersilia`; deliberately
not a `requirements.txt` dependency).
**Compute:** ~16 × ~206k predictions — an overnight job. `SMOKE=1` runs one model over
~1,500-compound subsets for a fast end-to-end check.

## 05_euopenscreen_validation.py
Evaluates how well the ChEMBL pathogen models predict the **EU OpenScreen** screening library from
the step-04 predictions, and renders the figures. Analyses, all joined on SMILES:
(1) own-assay AUROC for the 7 organisms with an EU OpenScreen primary assay; (3) exclusive vs
shared-hit AUROC; (4) a model × EU OpenScreen assay AUROC matrix (off-diagonal = a model predicting
a *different* organism's data) plus a per-model specificity index; (4c) the label-only hit-promiscuity
distribution (how many actives are hits in 1, 2, … 7 pathogens); (4d) the summed consensus score
across the 7 models by hit class (inactive / exclusive / narrow / broad), the maximum score by plain
activity (raw and percentile-normalised), and the rank each exclusive hit's own model gives it among
the 7, each of the last two also with training-set compounds removed. Writes the summary CSVs
(including
per-organism ROC-curve points and a leakage/overlap report) + `BasePlot` figures (png + pdf +
`figure_cells.json`) to `output/05_euopenscreen_validation/`. Analysis and figures are split, with
shared IO/metric primitives factored out: `src/eval_common.py` (load predictions, training keys,
merge + metrics), `src/eval_euopenscreen.py` (the EU OpenScreen analyses), `src/plots_euopenscreen.py`
(reads only the summaries).
**CoAdd is a separate step.** This script is EU OpenScreen only. The CoAdd validation (each model on
its CoAdd reference strain, both endpoints) lives in `src/eval_coadd.py` + `src/plots_coadd.py` and is
run by `06_coadd_validation.py` (documented below).
**Metrics:** AUROC, AUPRC, BEDROC(α=20), enrichment factor EF@1%/EF@5% — the ranking metrics matter
because external hit rates are very low (~0.01–0.4% active). Defined in `src/metrics.py` (ported from
the `new-modelling` repo).
**"Shared" organisms (7):** those with an EU OpenScreen primary assay — abaumannii, calbicans, ecoli,
efaecium, kpneumoniae, paeruginosa, saureus (`SHARED_ORGANISMS` in `src/default.py`, from
`primary_assays_manual.csv`). The other 8 models are ChEMBL-only and appear only as extra rows in the
cross-organism matrix (evaluated against the 7 assays, no own-assay diagonal).
**Labels — the primary assay only, one row per compound** (`load_euos_primary`). Read from
`02_binarised_assays/{primary_assay_id}.csv`, conclusive rows only (`bin ∈ {0,1}`), then collapsed to
one row per SMILES keeping the **highest bin ("Active prevails")**.
- **Not `02_merged/02_{code}.csv`.** That file is the union of *all* 5–6 live EU OpenScreen assays for the
  pathogen. The non-primary ones are small (~5,300 compounds vs ~101,000) but hit-**enriched 12×–60×**
  — they include academic compound batches and dose-response confirmation sets that run >50% active.
  Using it would multiply actives by 1.4×–3.5× (e.g. *A. baumannii* 57 → 199) and destroy the ~1e-4
  prevalence these ranking metrics are read against. The primary screen is the only uniformly
  screened, unbiased set. `02_merged` is used *only* for the InChIKey lookup here, and as
  `02_only_smiles.csv` for the step-04 prediction input.
- **Why the collapse.** The per-assay file is *measurement*-level (one row per assay well) and the
  library registers a few structures more than once as separate compounds, so those recur as repeated
  SMILES. Uncollapsed they were counted twice in `n_eval`/`n_active` and carried double weight in
  every metric, which put this loader ≤1 active out of step with the *compound*-level exclusivity
  subsets. Collapses are logged per assay: 2, 2, 0, 2, 7, 0, 2 rows for abaumannii, calbicans, ecoli,
  efaecium, kpneumoniae, paeruginosa, saureus.
- **The rule matches upstream deliberately** —
  `eu-openscreen-antimicrobial-tasks/scripts/02_binarise_and_merge.py::merge_pathogen_rows` uses the
  same max-bin rule when building `02_merged`; it just never applies it per assay. It decides one real
  case: *A. baumannii* has a **discordant** pair, two separately registered library compounds
  (EOS101879 @ 93.05 → active, EOS17004 @ 11.03 → inactive) sharing one structure and InChIKey. A
  replicate disagreement, now resolved to active by rule rather than by row order.
**Exclusivity (analysis 3):** uses the precomputed EU OpenScreen subsets in
`06_subset_data/exclusivity/` — "exclusive" = active in only 1 of the 7 primary assays (i.e.
organism-specific), "shared" (non-exclusive) = active in ≥2 (i.e. pan-active) — each scored against
the same primary inactives, with `n=` labels. This *is* the pan-active-vs-specific specificity test
(a drop toward 0.5 on the exclusive/specific bars = the model captures generic, not organism-specific,
activity); an on-the-fly recompute was checked to match these upstream subsets to ≤0.0002 AUROC, so
only this one is kept.
**Paper-figure panels (two rows of three).** `hit_promiscuity` (40 × 39.3 mm) and `euos_overlap`
(47 × 39.8 mm — wider so its library axis keeps a middle tick, plus a later +2 mm) have their own footprints
(`WIDE_PROMISCUITY` / `WIDE_OVERLAP`) and, more importantly, are matched on
**plotting-area height** (~25.8 mm) so their data rectangles line up when placed side by side — equal
crops alone would not do that, since `euos_overlap` spends more of itself on chrome. That is why
`hit_promiscuity`'s crop is the shorter of the two. 40 × 30 mm was tried and only works for
`hit_promiscuity`: `euos_overlap` has two x axes, whose ticks and labels eat its height twice, leaving
a 9.5 mm plotting area in which the 7 genus labels collided. Its axis names were then moved **inline**,
to the left end of each tick row (`Actives`, `Full library`), which cost nothing and took the plotting
area from ~19 to 25.8 mm. `hit_exclusivity_auroc` is at `EXCLUSIVITY_BARS` = `SMALL_SQUARE` + 4 mm of
width (~50 × 46 mm measured); the remaining three panels are at
`SMALL_SQUARE = (1.5, 1.5)` = 45 × 45 mm, below the step-03 row:
row A `hit_promiscuity`, `euos_overlap`, `exclusive_hit_model_rank_percentile_dedup`;
row B `consensus_max_percentile_by_activity_dedup`, `hit_exclusivity_auroc`, `submodel_auroc_summary`.
Three per row, not six — six at 45 mm is 270 mm against a 183 mm row. Measured crop widths: row A
127.1 mm, row B 128.3 mm, both clear of the 183 mm row. Stacked under the step-03 row (182.8 × 62.4 mm)
the three rows come to **148.4 mm of the 170 mm** page height (tallest panel per row: 62.4 + 39.8 + 46.2).
Everything else in this script keeps its larger footprint (inspection figures, not paper panels).

**Their legends are separate `*_key` panels** to place once in Illustrator: `euos_overlap_key`,
`exclusive_hit_model_rank_key`, `hit_exclusivity_auroc_key`, `submodel_auroc_summary_key`. Keys and
panels share their definitions so they cannot drift apart.

**Axis labels on those six are deliberately terse** — at 45 mm a sentence-length y label made
`consensus_max_percentile_by_activity_dedup` 77 mm tall. Model counts, dedup status and full phrasing
now live **in the caption**, so the caption is load-bearing for these panels.

**Active-set overlap (analysis 4b, `05_active_overlap.csv`):** label-only pairwise overlap between
the 7 assays' active sets. The hit-sets overlap substantially (most actives are pan-active), which is
*why* the cross-organism AUROCs (analysis 4) are high — high off-diagonal AUROC reflects shared
promiscuous hits, not necessarily cross-organism prediction. Read the two together.
**Two measures, two panels**, because the active sets are very unequal (14 for *P. aeruginosa* up to
378 for *S. aureus*) and a single number cannot serve both readings:
- `jaccard` = |A ∩ B| / |A ∪ B| → `active_overlap_jaccard`, a heatmap. Symmetric, but the union
  denominator lets the larger set dominate: *P. aeruginosa* shares 13 of its 14 actives with
  *A. baumannii* and still scores only 0.22.
- `containment` = |A ∩ B| / |A| → `active_overlap_containment`, a 45 mm matrix of pies where the
  **coloured wedge is the share of the row organism's actives that are also active against the column
  organism**. **Directional, so the matrix is not symmetric and both triangles must be read**
  (paeruginosa→abaumannii = 93%, abaumannii→paeruginosa = 23%) — the caption must say the row is the
  denominator or the panel will be misread. The diagonal is 1 by construction and renders as a full
  circle. Its key is the standalone `active_overlap_containment_key` strip.
  **Decision:** all pies are drawn at one size and each organism's active count is printed in its y
  tick label. Circle area encoded that count in an earlier version, but at 45 mm the cell pitch is
  2.4 mm, which put *P. aeruginosa* (14 actives) at 0.40 mm across — too small to read the wedge,
  which is the actual measurement. Row colour is the pathogen's step-03 `pathogen_activity_ratios`
  hue, shared with `exclusive_hit_model_rank_*` so one pathogen has one colour across both steps.
Compounds are matched on the **raw SMILES string** (consistent with analysis 4c; the leakage report
matches on InChIKey instead).
**Hit promiscuity (analysis 4c, `05_hit_promiscuity.csv` + `05_promiscuous_hits.csv` +
`hit_promiscuity`):** label-only, per-compound counterpart of the Jaccard view — for every distinct
EU OpenScreen active, in how many of the 7 primary assays it is a hit. **Decisions:** hits are
counted from the **primary** assays only (one assay per pathogen, so the axis runs 1–7); compounds
are matched on **SMILES**, consistent with the overlap analysis; and **all 588 distinct actives are
counted** (the union across assays), with no testing-coverage filter. 52 of them have a conclusive
result in fewer than 7 assays, so their count is a lower bound — nothing is dropped, but those
compounds are tallied per bin in `n_incomplete_coverage` and each compound carries its own
`n_assays_tested`. Result: 390 / 108 / 45 / 20 / 11 / 12 / 2 molecules hit 1…7 pathogens, i.e. two
thirds of hits are organism-specific singletons and only 2 compounds hit the whole panel. The
figure uses a log y axis (the distribution spans two orders of magnitude) with counts annotated.
**Summed consensus score by hit class (analysis 4d, `05_consensus_sum_boxstats.csv` +
`05_consensus_sum_actives.csv` + `consensus_sum_by_hit_class`):** for every library compound the 7
shared-organism models' `consensus_score` values are summed (range 0–7; all 7 models score the same
106,290 compounds with no missing values, so the sums are comparable), then compared across three
classes taken from the **primary** assays only: `inactive` = `bin == 0` in every assay where it has
a conclusive result and active in none, `exclusive` = hit in exactly 1, `narrow` = hit in 2–3, and
`broad` = hit in more than 3 (i.e. the shared/non-exclusive hits of analysis 3 split by breadth).
**Cutoff:** the narrow/broad boundary is `NARROW_MAX_PATHOGENS = 3` in `src/default.py` — a
requested classification, not a fitted value; change it there and both the class labels and the
figure follow. **Training-set compounds are deliberately kept** (raw only, no dedup), so this
describes the score distribution, *not* an out-of-sample performance estimate. 101,024 compounds
classify (100,436 / 390 / 153 / 45); the 5,266 compounds with no conclusive primary result cannot be
classified and are reported as `n_unlabelled_excluded`. Medians rise monotonically
3.41 → 3.76 → 4.31 → 4.68. The figure draws the boxes from the precomputed statistics (median/IQR,
1.5×IQR whiskers, fliers hidden — the inactive class is never shipped per-molecule) and overlays
the 588 individual actives as jittered points, seeded with `RANDOM_SEED`.
**Maximum consensus score, active vs inactive (`05_consensus_max_boxstats.csv` +
`05_consensus_max_actives.csv` + `consensus_max_by_activity`):** the same compounds and the same
conventions, but aggregating each compound's 7 scores with `max` (how confident the single most
confident model is, in [0, 1]) and collapsing to two classes — `inactive` as above vs `active` = a
hit in ≥1 pathogen. Medians 0.637 (n = 100,436) vs 0.712 (n = 588). Shares the analysis core
(`_consensus_scores_with_hit_counts`) and the panel base (`ScoreByHitClassPlot`) with the summed
variant, so the two figures are directly comparable. A **normalised twin** is written alongside it
(`05_consensus_max_percentile_{boxstats,actives}.csv` + `consensus_max_percentile_by_activity`,
`normalize=True`): each model's score is first replaced by its percentile within that model's own
library distribution, so the max is not biased towards whichever model happens to output the highest
values — the same transform as the `percentile` ranking below. That maximum is then **re-ranked over
the library** (`rerank=True`), because the maximum of 7 values is high by construction (for 7
independent uniforms E[max] = 7/8 = 0.875; observed library median 0.830, below 0.875 only because
the models correlate, mean pairwise Spearman ρ = 0.461). Re-ranking is monotone, so it leaves every
ordering — and therefore every ranking metric — bit-identical, and only makes the axis interpretable:
"this compound's best-model percentile beats that fraction of the library's best-model percentiles".
Medians on that axis are 0.491 (inactive, ≈0.5 by construction since inactives are 99.4% of the
library) vs 0.919 (active, Q1 0.648). The `score_scale` column in each boxstats CSV records which
scale a row came from (`raw` / `percentile_reranked`).
**Aggregate choice.** `max` was compared against mean, median, 2nd-highest, top-2/top-3 mean and
max−mean on active-vs-inactive separation (descriptive AUROC, raw scores): max-percentile 0.794 and
top-2 mean 0.795 lead, plain `max` on raw scores gets 0.734, and averaging over more models loses
separation (mean-percentile 0.741) because an organism-specific hit has one high score and six
mediocre ones. Not aggregating at all — one observation per (compound, pathogen) pair scored by that
pathogen's own model — separates best (0.839, 706,804 pairs) but changes the unit of observation, so
it is recorded here as an option rather than plotted.
**Own-model rank for exclusive hits (analysis 4e, `05_exclusive_hit_model_rank.csv` +
`05_exclusive_hit_model_rank_compounds.csv` + `exclusive_hit_model_rank_{raw,percentile}`):** for
each of the 390 exclusive hits (active in exactly 1 primary assay, so exactly one of the 7 models is
the "right" one), the 7 models' scores are ranked best-first and the position of the hit pathogen's
own model is recorded — rank 1 means its own pathogen scores it highest, rank 2 that one other
pathogen's model ranks it above its own, and so on to 7. The dashed line on the figure is chance
(390/7 = 55.7). **Two rankings, because the models' raw `consensus_score` values are not calibrated
to a common scale** (library medians range 0.353 for *E. faecium* to 0.606 for *P. aeruginosa*):
`raw` ranks the scores as-is and gives rank-1 = 80/390, but *C. albicans* then takes rank 1 for 185
of the 390 compounds and *E. faecium* for none — largely a per-model offset effect. `percentile`
converts each score to its percentile within that model's own library distribution first, which
evens the top-1 model distribution out and gives rank-1 = 114/390. Ties are resolved with
`method="min"` (best case for the true pathogen) and counted in `n_tied_scores` (2 raw, 1
percentile). Both panels are written; read the percentile one as the scale-corrected version.
**Output layout: `full/` vs `deduplicated/`.** Results are filed by leakage status so the two can
never be confused (`FULL_SUBDIR` / `DEDUP_SUBDIR` in `src/default.py`). `full/` holds the analyses
that KEEP the compounds the models were trained on (15 CSVs, 7 figures); `deduplicated/` holds those
with every training-set compound REMOVED (13 CSVs, 14 figures); the top level keeps only what has no
leakage dimension — the label-only tables (`05_active_overlap`, `05_hit_promiscuity`,
`05_promiscuous_hits`) and the leakage audit itself (`05_leakage_report`) — 4 CSVs, 6 figures. Each
dir has its own `png/`, `pdf/` and `figure_cells.json`. The five long-form metric tables
(`05_euopenscreen_auroc`, `_roc`, `05_hit_exclusivity`, `_percentiles`, `05_cross_organism_euos`)
carry both variants in a `set` column, so they are written into **both** subfolders filtered to that
folder's variant — each subfolder is self-contained. The AUROC-family figures live under
`deduplicated/` because they plot the leakage-filtered values (`EuosRocGridPlot(set_name="dedup")`,
`_prefer` picking dedup, the dedup-derived specificity index); their `full/` counterparts exist as
data, not as figures. `individual_performance/` and `eos3dys_validation/` are separate analysis
families and keep their own layout, with both variants in their `set` column.
**Leakage-filtered twins (`*_dedup*`).** The max-consensus and rank figures each get a second
version with every compound present in **any of the 7 models' ChEMBL training sets** removed
(`shared_training_inchikeys` + `euos_inchikeys`; all 101,024 classified compounds carry an InChIKey,
so coverage is complete, and a compound without one would be kept and logged rather than dropped on
missing information). **Any-model, not own-model,** filtering: these two figures compare all 7 models
against each other, so leakage in a *non-own* model inflates that model's score and pushes the own
model *down* the ranking. `05_consensus_max_percentile_dedup_*` filters **both** classes (inactives
100,436 → 97,162, actives 588 → 428; 3,434 dropped in total) so the comparison stays like-for-like —
medians move 0.491 → 0.492 (inactive) and 0.919 → 0.867 (active), i.e. most but not all of the
active-class advantage survives. `05_exclusive_hit_model_rank_dedup*` drops 68 of the 390 exclusive
hits → 322, and rank-1 goes from 114/390 (29.2%) to 88/322 (27.3%) against a chance of 46.0. A
`leakage` column (`raw` / `dedup`) tags every row of all four files.
**Bars are stacked by the hit's own pathogen** (`n_<code>` columns in the summary CSV, colours from
`SHARED_ORGANISM_COLORS` in `src/plotting_colors.py` — one hue per shared organism, the 7 usable
ArticleColors categoricals). Bar height stays the molecule count so the chance line keeps its
meaning. **Read segments with care:** they scale with how many exclusive hits each pathogen has
(*S. aureus* 198, *C. albicans* 88, *K. pneumoniae* 66, *E. coli* 23, *A. baumannii* 10,
*E. faecium* 4, *P. aeruginosa* 1 — totals repeated in the legend), so a big segment is not a
per-pathogen effect. On the percentile ranking the per-pathogen rank-1 rates are *S. aureus* 67/198,
*C. albicans* 28/88, *K. pneumoniae* 14/66, *E. coli* 3/23, *A. baumannii* 1/10,
*P. aeruginosa* 1/1, *E. faecium* 0/4 — available per compound in
`05_exclusive_hit_model_rank_compounds.csv`.
**Leakage:** every metric is reported both `raw` and InChIKey-`dedup` (removing compounds in the
model's ChEMBL training set, read from `chembl-antimicrobial-models/output/07_datasets/{code}/*.csv`).
If that sibling repo is absent, dedup silently becomes a no-op and only `raw` rows are written
(logged). Dedup AUROCs run slightly below raw, as expected.
**Headline score:** `consensus_score` where present; single-dataset pathogen models (campylobacter,
hpylori, ngonorrhoeae) have no consensus and expose one output column, which is used as the headline
(logged). Multiple outputs with no consensus would be skipped, never silently averaged.
**Partial data:** missing step-04 model files are skipped with a logged message. The EU OpenScreen
predictions are complete (16/16 models); the script re-runs cleanly if any are re-generated.
**Individual sub-model performance** (`individual_performance/` subfolder): per shared pathogen, looks
at every sub-model output column (not just `consensus_score`) to check whether the ensemble members
agree. Two panels per pathogen: (a) each sub-model's AUROC on the pathogen's own assay (dedup, dot plot,
consensus starred) — a wide spread means members differ in quality; (b) a Spearman correlation heatmap
between sub-model scores over the library — low off-diagonal values mean members rank compounds
differently. Correlation is label-free, so computed on all rows; AUROC uses the dedup set.
A **cross-pathogen summary** (`submodel_auroc_summary`, top-level figure) overlays, per organism,
every sub-model's own-assay AUROC (small dots) with the `consensus_score` as a larger bubble —
the whole-hub view of where the consensus sits within its ensemble's spread.
**CoAdd model on EU OpenScreen** (`eos3dys_validation/` subfolder): the CoAdd model `eos3dys` emits
many endpoints (`{organism}_{strain}_{inhib_50|mic_25}` + cytotoxicity/hemolytic); each is scored
against every EU OpenScreen primary assay → `eos3dys_euos_auroc.csv`, drawn as an endpoint×assay
heatmap (diagonal = matching organism) plus a same-vs-different-organism AUROC swarm. This is the
*opposite direction* from the ChEMBL-models-on-CoAdd step (`src/eval_coadd.py`): it tests the CoAdd
model's generalization to EU OpenScreen. Dedup here removes each endpoint's own CoAdd-training
compounds (`data/raw/coadd_data/{03_binarised_inhibition,05_binarised_mic}/{strain}.csv`). Endpoints
are heterogeneous (strains, inhibition vs MIC), so read the matrix as indicative. Code:
`src/eval_eos3dys.py` + `src/plots_eos3dys.py`.
**Two own-assay eos3dys analyses mirroring the ChEMBL ones**, over the **6 organisms with both an
EU OpenScreen assay and an eos3dys endpoint** — abaumannii, calbicans, ecoli, kpneumoniae,
paeruginosa, saureus (`matched_endpoints`; *E. faecium* has an EUOS assay but no eos3dys endpoint and
is logged as skipped). Labels are always EU OpenScreen — the out-of-sample direction for this model —
and CoAdd serves only as the dedup source; scoring eos3dys against CoAdd labels would be in-sample.
- **Training overlap** (`eos3dys_overlap_report.csv` + `eos3dys_overlap_{inhib_50,mic_25}`, the shared
  `EuosOverlapTwinPlot` reused once per endpoint metric): how much of the EU OpenScreen library, and
  of its actives, eos3dys already saw in CoAdd training — matched on **InChIKey**, because the two
  sources standardise SMILES differently. Only ~450 of the ~101k EUOS compounds are in the
  `inhib_50` training set (~95 for `mic_25`), but they are enriched in actives: e.g. 31 of
  *S. aureus*' 378 actives and 18 of *C. albicans*' 171, and 3 of *P. aeruginosa*' 14. So dedup
  barely changes the library but removes a real slice of the positives.
- **Exclusive vs shared hits** (`eos3dys_hit_exclusivity.csv` + `eos3dys_hit_exclusivity`): the twin
  of the ChEMBL `hit_exclusivity_auroc` panel, using the **same upstream exclusivity subsets** so the
  two compare directly (those subsets are defined over all 7 assays, including the *E. faecium* one
  eos3dys cannot score). Four bars per organism — `inhib_50`/`mic_25` × exclusive/shared, subset by
  hue and endpoint by saturation, with `n=` labels. Dedup AUROCs: shared hits score 0.80–0.95
  everywhere, exclusive hits drop sharply for the Gram-negatives (*E. coli* 0.48, *K. pneumoniae*
  0.56 on `inhib_50`) but hold up for *S. aureus* (0.85) and *C. albicans* (0.88). Several exclusive
  subsets are tiny (*P. aeruginosa* n=1, *A. baumannii* n=9–10) — read those bars off the counts.
**Strain pinning (bug fixed 2026-07-29).** Several organisms have more than one strain endpoint per
metric (`ecoli` has ATCC25922/lpxC/tolC, `paeruginosa` has ATCC27853/PAO1/PAO397). `matched_endpoints`
now pins the strain to `COADD_REF_STRAINS` — the same reference strain the CoAdd step uses — and logs
any fallback. Before the fix it silently kept whichever strain came last, so *E. coli* was scored with
the efflux-deficient `tolC` mutant (73,391 training compounds) and *P. aeruginosa* with `PAO397`
(14,373) instead of their wild-type references (~81,600 each). Mutant strains are not interchangeable
with the reference, so any earlier *E. coli* / *P. aeruginosa* eos3dys number is superseded.
**Three further eos3dys analyses, all over the same 6 matched organisms and all EU OpenScreen-labelled:**
- **Own-organism rank for exclusive hits** (`eos3dys_exclusive_rank{,_compounds}.csv` +
  `eos3dys_exclusive_rank_{raw,percentile}`): each organism gets ONE score per compound —
  `inhib_50 + mic_25` summed, since both are probabilities from one model (summed rather than
  max'd because `mic_25` runs systematically above `inhib_50`, 50–86% of compounds depending on
  organism, so a max would just return the MIC head). The 6 organisms are ranked and the position of
  the hit's own organism recorded. **Exclusivity is recounted over the 6 matched organisms**, not read
  from the upstream 7-assay subsets, so the ranking candidates and the exclusivity definition span the
  same set — a compound hit in *E. faecium* plus one matched organism counts as exclusive here but
  shared in the ChEMBL panels. Dedup drops hits in ANY matched endpoint's CoAdd training set (25 of
  438 → 413). Both rankings are written, as for the ChEMBL models: rank-1 = 99/413 raw vs **137/413
  percentile** (chance 68.8), the gap again reflecting per-organism score-scale differences.
- **Max endpoint probability, active vs inactive** (`eos3dys_consensus_max_{boxstats,actives}.csv` +
  `eos3dys_consensus_max_by_activity`): per compound the highest of the 12 matched endpoint
  probabilities. **No percentile normalisation**, unlike the ChEMBL twin: the 7 ChEMBL models' maxima
  ran 0.753–0.974 so low-scoring ones could never win a max, whereas all 12 eos3dys endpoints top out
  at ~1.0. Dedup medians 0.734 (inactive, n=100,007) vs 0.865 (active, n=535). Note the max of 12
  correlated probabilities is high by construction, so read the gap, not the absolute level.

### Score-distribution boxes

**No fills, and `inactive` is silver** — a background class of 97,162 against 428 should not
carry the visual weight of a finding. The activity panel is 40 x 31.5 mm; 1.05 cols is the floor
while the tick labels keep their `(n = ...)` counts, which are kept deliberately.

**Boxes with points are unfilled.** Only `inactive` keeps a fill — it has ~10^5 compounds and is
never shipped per-molecule, so there is nothing to reveal behind it. Every class that ships points
shows its swarm through an unfilled box, matching the step-01 technical boxes. Previously all
classes were filled, hiding the 428 active points behind their own box.

**The `consensus_max_*` box stats now carry an `auroc` column** and the panels print it: **0.7574**
(percentile max, dedup), 0.7943 (percentile max, full), 0.7335 (raw max, full).
**It cannot be derived from the shipped actives.** The "AUROC = mean percentile of the positives"
identity needs the percentile ranked over the evaluated set, but this score is ranked over the
whole library and then filtered — and the filtering is not random, since dedup drops high scorers
preferentially. The shortcut gives 0.7511 against a true 0.7574.

**ROC curves accompany the box panels** (`05_consensus_max_roc.csv`,
`05_consensus_max_percentile_roc.csv`, `05_consensus_max_percentile_dedup_roc.csv`; 800 thinned
vertices each). Same data as the boxes, read as a retrieval trade-off rather than a pair of
distributions — the dedup curve reaches TPR 0.5 at roughly FPR 0.12, i.e. half the hits are
recovered in the top ~12% of the library.

**Normalisation is worth ~0.06 AUROC** (0.7335 raw max → 0.7943 percentile max, full). Only the
final re-rank is monotone; putting the 7 models on a common scale before the max genuinely
improves separation.

### Hit-exclusivity event plot (`hit_exclusivity_events`)

The same analysis as the `hit_exclusivity_auroc` bars, shown as a distribution instead of a
summary: one vertical line per hit at its percentile in that organism's model ranking, periwinkle
lane for shared hits and amber for exclusive, with the pair of AUROCs printed to the right as
`shared / exclusive`, each value in its lane's colour. 40 × 40 mm (`EVENT_SQUARE`), smaller than
the ~49 mm AUROC panel because the value column sits outside the axes. Reading: mass pushed to the right = hits upranked. A respectable-looking bar
can be a handful of hits at the very top over an otherwise flat spread, and only this panel
distinguishes the two.

**The right axis is not a second variable.** AUROC is the mean rank percentile of the positives,
so the printed number is where the lane's centre of mass sits. The export asserts that identity
(`_check_percentiles_match_auroc`) rather than trusting it.

**Lines, not a density — deliberately.** *P. aeruginosa* exclusive has n=1 and *E. faecium* n=4.
A smoothed curve over one point invents a shape; one tick per compound cannot. It also avoids
picking a bandwidth or a bin width.

**No genus labels — it is placed against `euos_overlap`'s.** Both panels take their row order from
one `overlap_row_order` helper (`n_active` ascending), and this panel's axes height and y limits
are matched to that panel's so the rows line up. Label suppression and row order are a single
decision in the code, so unlabelled rows in the wrong order cannot happen.

**`*` = AUROC over 10 hits or fewer** (P. aeruginosa shared n=3 and exclusive n=1, E. faecium
exclusive n=4, A. baumannii exclusive n=8). Flagged, not filtered — **the caption must define the
mark.**

**Layout is measured, not hand-tuned.** The AUROC column's distance from the plot is derived from
the rendered width of the widest value, leaving a 3 pt clearance. The hand-set offset it replaced
wasted ~4 mm, which on a 40 mm panel cost 45% of the plotting area.

**New summary CSV:** `05_hit_exclusivity_percentiles.csv`, one row per active
(`pathogen, code, eosid, subset, set, inchikey, smiles, percentile`), 732 rows for dedup. Written
next to the AUROC record it unpacks.

**Tolerance note:** the AUROC-identity check runs at `tol=1e-4` and cannot be tightened —
`metrics.compute_metrics` stores `round(auroc, 4)`, so the comparison is against a value already
quantised to 1e-4 and the measured worst deviation is 4.86e-05, which is the rounding rather than
a defect.


## 06_coadd_validation.py
The mirror of step 05 on **CoAdd**: the public ChEMBL pathogen models scored against CoAdd
growth-inhibition / MIC labels, from the same step-04 predictions. Writes to
`output/06_coadd_validation/`; analysis in `src/eval_coadd.py`, figure in `src/plots_coadd.py`.
CoAdd is richer (many strains/cutoffs per organism), so this first version is deliberately narrow —
**shared organisms on their reference strain only**.
**Organisms (8, `COADD_REF_STRAINS`):** every organism with both a ChEMBL model and a CoAdd
reference strain — the 7 EU OpenScreen-shared ones **plus S. pneumoniae** (`ATCC700677`). efaecium and
spneumoniae are **MIC-only** (no single-point file, so `inhib_50` is skipped and logged).
**Endpoints:** the reference strain at two headline cutoffs — `inhib_50` (≈ ChEMBL single-point 50%)
and `mic_10` (≈ ChEMBL dose-response 10 µM) → own-strain AUROC bars per organism, grouped by endpoint
(`06_coadd_auroc.csv` + the `coadd_shared_auroc` figure).
**Leakage** (`06_coadd_leakage_report.csv`): overlap between each model's ChEMBL training InChIKeys
and its CoAdd reference-strain data; metrics reported `raw` + `dedup` as in step 05.
**Not evaluated here:** the CoAdd model `eos3dys` (it is validated the *other* way — on EU OpenScreen —
in step 05's `eos3dys_validation/`). **Partial data:** missing step-04 CoAdd predictions are skipped
with a logged message; re-runs cleanly as step 04 completes. **Deferred:** the full multi-strain /
multi-cutoff CoAdd matrix, cross-organism/specificity, ROC grid and per-submodel breakdown.
## 07_score_matrices.py
Builds the foundation of the correlation analysis: the full-library prediction matrices. All annotation
models were run by `00_download_data.py` on the same ~1.35M-compound reference library
(`data/processed/annotation_preds_ref_library/`), aligned on the `key` hash — a clean rectangular
matrix. Produces **five** matrices, each 1,355,109 x 260, plus the parquet cache behind them. Engine in
`src/eval_correlations.py`; writes to `output/07_score_matrices/`.

| matrix | transform |
|---|---|
| named | raw scores, columns renamed `{pathogen_code}__{model_id}__{column_name}` |
| z-score | `(x - mean) / std` per column |
| rank-percentile | percentile rank within each column, bounded [0, 1] |
| z-score + L2 row-norm | each compound's profile divided by its Euclidean norm |
| rank-percentile + L1 row-norm | each compound's profile divided by its absolute sum (rows sum to 1) |

**Key decisions, for review:**
- **Which endpoints:** only `selected == Yes` rows of the manually curated
  `config/08_endpoint_selection.csv`. That file, not any automatic organism/model filter, is the single
  source of truth. Template generator: `tools/build_endpoint_selection_template.py`.
- **Pathogen codes** come from `config/pathogens_of_interest.csv` where available (12 of 56 organisms);
  the remaining 44 (mostly single-model gut-microbiome species) get a mechanical fallback — first letter
  of genus + species epithet, lowercased (`Bacteroides caccae` -> `bcaccae`). Every fallback is printed
  for audit, and a code collision raises rather than silently merging two organisms.
- **Column scaling vs. row normalization are independent transforms**, composed in that order. The two
  row norms are deliberately *different*: L1 suits the non-negative percentile matrix (giving a
  compositional vector — each endpoint's relative share of that compound's activity), L2 suits the
  signed, mean-centred z-score matrix (the norm cosine similarity is built on). L1 on signed data has
  no equally clean reading.
- **Caveat, not a cutoff:** row normalization amplifies a compound whose profile is near-flat/near-zero
  everywhere (small vector / small norm). Expected for the technique, but such a profile must not be
  read as strong endpoint preference. Rows with a degenerate norm (< `ROW_NORM_EPS` = 1e-12) divide to
  inf/NaN and are counted and reported — never dropped or clipped. None occur in the current selection.
- **Pre-existing missing values are kept, not imputed.** 15 of 1,355,109 compounds lack one endpoint
  (`pfalciparum__eos4zfy__maip_score`); each is normalized over its 259 available endpoints and the
  missing cell stays NaN. Counted and printed; nothing filled or dropped.
- **CSV, not parquet**, at explicit user request: ~26 GB of text across the five versus ~7 GB as float32
  parquet, and noticeably slower to write and read. A real tradeoff, accepted deliberately. The raw
  prediction CSVs (~15 GB) are read at most **once**, into `07_score_matrix_full.parquet`.

## 08_pathogen_jaccard.py
Asks whether endpoints of the SAME pathogen pick out the same compounds more than endpoints of
DIFFERENT pathogens do. For every pair of the 260 endpoint columns, the Jaccard overlap of their
top-1000 highest-scoring compounds, aggregated per pathogen into a same-vs-different box pair.
Figure in `src/plots_matrix_analyses.py`; writes the 260x260 Jaccard matrices (reused downstream),
a per-pathogen summary CSV, and PNG + PDF to `output/08_pathogen_jaccard/`.
**Key decisions, for review:**
- **Cutoff: top 1000** of 1,355,109 compounds per column — user-directed.
- **Minimum endpoints per pathogen: 5** (default; a command-line argument, so several thresholds can be
  run and compared rather than one overwriting another — outputs are suffixed `min<K>`). A pathogen must
  have at least K endpoints **of its own**; those below are removed from the analysis *entirely*, ceasing
  to be different-pathogen partners too, not merely losing their own box. **User-directed, not fitted.**
  At `min5`: 209 of 260 columns, 11 of 56 pathogens. At `min2` (the weakest defensible bar, 2 endpoints
  being the minimum for any same-pathogen pair to exist): 219 columns, 15 pathogens — but `bfragilis` and
  `cneoformans` rest on a single pair each and `ngonorrhoeae`/`lmajor` on three, which is why 5 is the
  more defensible default.
- **Three matrix variants, not five.** Top-N Jaccard depends only on each column's own internal ranking,
  and both column scalings are strictly increasing per column — so z-scoring and rank-percentiling cannot
  change any column's top-1000 set. The unscaled and both scaled matrices give the same result, computed
  once; only the row-normalized matrices change rankings. The script **asserts** this at runtime rather
  than assuming it, and the assertion earns its keep: baseline == rank-percentiled holds exactly, but
  baseline == z-scored comes back **False** — `(x - mean) / std` in float32 reorders near-tied values in
  one column of 260 (`lmajor__eos60mw__leishmania_*`, 155 of its top-1000 shift), touching 138 of 67,600
  cells at ~0.001.
- **Same-model pairs are included** (the literal "each column against all others"). This matters at
  pathogen level: at `min5`, **3 of the 11 surviving pathogens** (`espp`, `spneumoniae`, `efaecium`) have
  `n_same_pairs_excl_same_model == 0` — every same-pathogen pair they have comes from one model's multiple
  output columns, so their box is that model agreeing with itself and says nothing about cross-model
  specificity. `espp` tops the `min5` baseline ranking on exactly that basis. The summary CSV carries
  `same_median_excl_same_model` alongside `same_median`; read them together.
- **Linear x-axis** (user-directed). Values bunch near zero so the different-pathogen boxes render as thin
  slivers, but exact-zero pairs are shown rather than silently dropped by a log axis. Nothing is filtered.
- **Plain matplotlib, not the 3 cm-cell publication grid** — a deliberate, user-approved departure from
  `docs/figure_conventions.md`, since the per-pathogen endpoint and pair counts on the y-axis go illegible
  at page width. PNG *and* PDF are still written; the vector copy is the readable one.

## 09_mean_rank_distribution.py
Collapses the rank-percentile matrix along its columns to one number per compound — that compound's
**mean percentile rank** across the 260 endpoints — and shows the distribution over the full library as a
histogram. Answers "how highly does this molecule rank, on average, across everything we predict".
Figure in `src/plots_matrix_analyses.py`; writes to `output/09_mean_rank_distribution/`.
**Key decisions, for review:**
- **Mean, not sum.** The request was phrased as a column-wise sum; the reported quantity is the mean
  (sum / 260), which is what "average rank" means and keeps the value on the interpretable [0, 1]
  percentile scale. They differ only by that constant factor — multiply by `n_endpoints` for the sum.
- **The centre of this distribution is fixed by construction, not a finding.** Each percentile-rank column
  has mean ~= 0.5, so the grand mean of the row means is ~= 0.5 necessarily (observed 0.5000). Only the
  **spread and shape** carry information: observed SD 0.133 against the ~0.018 that 260 mutually
  independent endpoints would give, and a range of 0.140-0.869 (no compound ranks in the top decile of
  everything, nor the bottom decile of everything). What that ~7x inflation means substantively is not
  something these numbers settle.
- **No rows dropped for missing values** — the 15 partially-scored compounds are averaged over 259
  endpoints, with `n_endpoints` per row in the CSV so they stay visible.
- Follows the standard publication figure convention (`BasePlot`, PNG + PDF + `figure_cells.json`) — a
  single distribution fits the page grid, unlike the step 08 diagnostic.

## 10_reference_library_projection.py
Visualises where the reference library sits in chemical space, highlighting each pathogen's
highest-predicted-activity compounds. `eos1klk` (fetched by `00_download_data.py`) computes four 2D
projections of the library — PCA, UMAP, t-SNE, TMAP; each gets its own figure: a silver background of the
full library's density, with each of the 15 pathogens in `config/pathogens_of_interest.csv` getting its own
panel showing its `PROJECTION_TOP_N` highest-scoring compounds (by `consensus_score`) in crimson. This is a
**rank cutoff (top-N), never a score threshold** — no `consensus_score` cutoff value is chosen or reviewed,
only a compound count. Analysis in `src/eval_projection.py`, figures in `src/plots_projection.py`; writes to
`output/10_reference_library_projection/`.
**Memory:** only one pathogen's `key` + `consensus_score` columns (of each ~424 MB prediction file) are read
at a time, immediately reduced to its `PROJECTION_TOP_N` highest rows and discarded — no more than one
pathogen's raw scores, and never all 15 prediction files, are held in memory together. The score ranking
does not depend on projection method, so each pathogen's top-N is computed once and reused across all four.
**Key decisions, for review:**
- **Top N per pathogen** (`PROJECTION_TOP_N`, in `src/default.py`) is **1000** — user-directed.
- **Grid resolution** (`PROJECTION_BINS`) for the silver background density is **60x60** cells per method —
  sized to stay legible at the ~30x30 mm panel size in a 15-pathogen small-multiples grid, not a fitted
  value.
