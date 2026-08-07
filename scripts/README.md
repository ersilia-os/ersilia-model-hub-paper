# Scripts

Numbered scripts run in order. Figure and plotting conventions (sizing, formats, per-figure
layouts) live in [`docs/figure_conventions.md`](../docs/figure_conventions.md).

## 00_download_data.py
Stages all input data in four sections: companion repos / eosvc (EU OpenScreen tasks, CoAdd data,
ChEMBL model reports and curation summaries), public GitHub files (Ersilia reference library,
DrugBank), Airtable model metadata, and Isaura precalc predictions for Ready annotation models.
Skip-if-exists, so it is safe to re-run. `--skip-isaura` skips the slow Isaura section; `--eosvc`
pulls Section 1 from eosvc instead of the companion repos.
Section 4 also fetches two Representation models explicitly, since the automatic loop is filtered to
`Task == "Annotation"` and never sees them: `eos1klk` (2D projector, Subtask=Projection) into
`data/processed/eos1klk_projection/eos1klk_v1.csv`, feeding step 10; and `eos4djh`
(datamol-basic-descriptors, Subtask=Featurization) into
`data/processed/annotation_preds_ref_library/eos4djh_v1.csv`, feeding step 10. Both pin the major
version only (`v1`), matching the loop. **`eos4djh` lands in the annotation folder despite not being an
annotation model** — the property steps read every property model from that one directory, and the model's real
Task is recorded on `PHYSCHEM_MODEL_ID` in `src/default.py` rather than in the path.
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

**Upstream mis-annotation — FIXED 2026-08-06.** `eos93h2` (`image-mol-gpcr`) was recorded in Airtable as
Representation / Featurization with Target Organism = *Homo sapiens*, when it is really a GPCR activity
predictor. It was the only thing breaking the "named organisms are Annotation-only" rule that the
Biomedical Area and Target Organism panel colours depend on. The correction has landed and it moved
**both** fields together, as required: `Task` Representation → Annotation and `Subtask` Featurization →
Activity prediction. Nothing else in the record changed. Script 01 already drew it as an Annotation row,
so **no panel needs changing**. Re-verified on the fresh metadata: 90 Ready models carry a named Target
Organism and 0 are non-Annotation; 103 carry a named Biomedical Area and 0 are non-Annotation — the rule
now holds with no exceptions. The Task/Subtask coupling warning in `docs/figure_conventions.md` still
applies to any future re-annotation, since a Task-only edit leaves the subtask-coloured panels
disagreeing with the task-coloured ones **without erroring**.

**`eos93h2` is now visible to `00_download_data.py`'s Section 4 loop**, which filters
`Task == "Annotation" & Status == "Ready"` — a test it previously failed. It has `Release = v1.0.0` and
no staged prediction file, so a re-run will try to fetch it from Isaura. It is a **human GPCR** model
with no antimicrobial relevance and appears in no config file, so its 10 endpoints do not enter
`config/08_endpoint_selection.csv` without a deliberate decision.

**Metadata snapshot refreshed 2026-08-06: 219 → 225 models, Ready 208 → 214.** Every "n = 208" below
predates it. The previous snapshot is preserved at `tmp/airtable_metadata_20260730_backup.csv`. Six
models were added — `eos19dk` (molcompass, Projection), `eos5g6m` (glacier-embeddings), `eos5mnx`
(sand-shape-descriptor), `eos6pj2` (nafm-embeddings), `eos8zvb` (pymolgen) and `eos84nf`
(genmol-scaffold-decoration, *In progress*) — **none of them an Annotation model**, so none adds
bioactivity endpoints. Four models changed Status: `eos8vud`, `eos69e6` and `eos4qda` became Ready,
while **`eos18ie` (antibiotics-ai-saureus) and `eos1lb5` went to *In maintenance***. Ready-only Task
counts: Annotation 131 → 130, Representation 58 → 61, Sampling 19 → 23.
**`eos18ie` is a selected endpoint whose model is now in maintenance** — it is row 1 of
`config/08_endpoint_selection.csv` (`saureus_inhibition_probability`, `Yes`) and is in the step-07
cache, so it stays in the score matrices while dropping out of every `Status == "Ready"` metadata panel.
That inconsistency is unresolved and needs a decision.

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
matrix. **Defines** five matrices, one row per compound and one column per `selected == Yes` endpoint
(**300** as the config stands), but by default **writes only the parquet cache** behind them plus the
mean-rank outputs below. Engine in `src/eval_correlations.py`; writes to `output/07_score_matrices/`.

**The five CSVs are opt-in (`--write-matrix-csvs`), and off by default.** They are ~23.5 GB of text and
~1 h 54 m of pure serialisation, and **no code in this repo reads them** — every downstream step
(08–14) re-derives what it needs from `07_score_matrix_full.parquet`, which is columnar and ~15x
smaller. `src/eval_auroc_matrix.py` says so explicitly at its top. They exist as human-readable exports
only, so a rebuild no longer pays for them unless asked:

```
python 07_score_matrices.py                      # parquet cache + mean-rank figure  (~5 min)
python 07_score_matrices.py --write-matrix-csvs  # also the five CSV exports        (+~1 h 54 m)
```

**Measured, so the tradeoff is reviewable rather than asserted** (2026-08-07, this hardware): building
the base matrix from the parquet takes **0.74 s**, deriving z-score + L2 row-norm **3.5 s** and
rank-percentile + L1 row-norm **51.1 s** — **55 s** for all three. Reading the three equivalent CSVs
would take **~84 s** at the measured 0.181 GB/s. So deriving is *modestly* faster per run, not
dramatically; what actually decides it is the **one-off ~1 h 54 m write** to save ~29 s per run, plus
the drift risk — a materialised copy can fall out of sync with the config, and did (the CSVs sat at 260
columns while the selection said 300). The parquet *is* the saved intermediate; the CSVs were a second,
redundant copy of it in a slower format.

**Only 3 of the 5 are ever consumed**, and step 08 is the only consumer — it re-derives them in memory
rather than reading any CSV. The other two are exports for a human, not pipeline inputs:

| matrix | transform | used downstream? |
|---|---|---|
| named | raw scores, columns renamed `{pathogen_code}__{model_id}__{column_name}` | **yes** — step 08 `baseline` |
| z-score | `(x - mean) / std` per column | no |
| rank-percentile | percentile rank within each column, bounded [0, 1] | no |
| z-score + L2 row-norm | each compound's profile divided by its Euclidean norm | **yes** — step 08 `zscore_l2rownorm` |
| rank-percentile + L1 row-norm | each compound's profile divided by its absolute sum (rows sum to 1) | **yes** — step 08 `rankpct_l1rownorm` |

**Why the two column-scaled-only variants are unused:** top-N Jaccard depends only on each column's
internal ranking, and both column scalings are strictly increasing *per column*, so neither can change
a column's top-1000 set — their Jaccard matrices are identical to `baseline`'s. Row normalisation is
different: it mixes values *across* a row, so it genuinely reorders each column, which is why the two
`*_rownorm` variants are real inputs. Step 08 asserts this rather than assuming it (see its section).
Keeping the two in `VARIANTS` is deliberate: a plain z-score or percentile table is the interpretable
thing to hand someone, which is why they are behind the flag rather than deleted.

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
- **CSV, not parquet**, for those exports when they are requested, at explicit user request: ~23.5 GB of
  text across the five versus ~7 GB as float32 parquet, and noticeably slower to write and read. A real
  tradeoff, accepted deliberately — and since 2026-08-06 one paid only on demand, because the format
  argument applies to a human reading the file, not to the pipeline, which uses the parquet. The raw
  prediction CSVs (~15 GB) are read at most **once**, into `07_score_matrix_full.parquet`.
- **The `260`-based numbers throughout the step 08–14 sections below were measured under the earlier
  260-endpoint selection.** Steps 07–14 were rebuilt from scratch at 300 endpoints on 2026-08-06; where
  a re-measured value is known it is given alongside, but treat any bare `260`/`67,600`/`26,260` as
  historical.

### Mean percentile rank per compound (merged in from the former step 09, 2026-08-06)

Collapses the rank-percentile matrix along its columns to one number per compound — its **mean
percentile rank** across all selected endpoints — and shows the distribution over the full library as a
histogram. Answers "how highly does this molecule rank, on average, across everything we predict".
Writes `07_mean_rank_per_compound.csv`, `07_mean_rank_quantiles.csv` and the
`07_mean_rank_distribution` panel (PNG + PDF + `figure_cells.json`) into `output/07_score_matrices/`.
Figure class in `src/plots_matrix_analyses.py`.

**Why it lives here now:** it needs the named matrix and the `rank_pct` scaling, both of which this
script already has, and it consumed nothing else. As a separate step it rebuilt the same matrix from
the same parquet a second time. `base` is now built once and serves both sections; the mean-rank block
is **skip-if-exists** on `07_mean_rank_per_compound.csv`, matching the five matrix CSVs. The former
`scripts/09_mean_rank_distribution.py` was removed from the tree (recoverable via
`git show <rev>:scripts/09_mean_rank_distribution.py`), and its stale 260-endpoint
`output/09_mean_rank_distribution/` was deleted.

**Key decisions, for review:**
- **Mean, not sum.** The request was phrased as a column-wise sum; the reported quantity is the mean,
  which is what "average rank" means and keeps the value on the interpretable [0, 1] percentile scale.
  They differ only by a constant factor — multiply by `n_endpoints` for the sum.
- **The centre of this distribution is fixed by construction, not a finding.** Each percentile-rank
  column has mean ~= 0.5, so the grand mean of the row means is ~= 0.5 necessarily (observed **0.5000**
  at 300 endpoints). Only the **spread and shape** carry information: observed **SD 0.1289** against the
  ~0.0167 that 300 mutually independent endpoints would give — a ~7.7x inflation — over a range of
  **0.1518–0.8616**, i.e. no compound ranks in the top decile of everything nor the bottom decile of
  everything. What that inflation means substantively (shared chemistry, correlated training data,
  genuine broad-spectrum compounds) is **not** something these numbers settle.
- **Quantiles** at 300 endpoints: 0.1973 / 0.2379 / 0.2887 / 0.4004 / 0.5023 / 0.5997 / 0.7059 / 0.7580
  / 0.8017 for the 0.1st, 1st, 5th, 25th, 50th, 75th, 95th, 99th and 99.9th percentiles.
- **No rows dropped for missing values** — the 15 partially-scored compounds
  (`pfalciparum__eos4zfy__maip_score`) are averaged over the 299 endpoints they have, with
  `n_endpoints` per row in the CSV so they stay visible. NaN-skipping mean, never zero-filled.
- Follows the standard publication figure convention (`BasePlot`, PNG + PDF + `figure_cells.json`).
  **This is the only figure step 07 produces**, so `figure_cells.json` in that folder holds exactly one
  entry.
- Previous 260-endpoint values, for comparison: mean 0.5000, SD 0.133, range 0.140–0.869.

## 08_pathogen_jaccard.py
Asks whether endpoints of the SAME pathogen pick out the same compounds more than endpoints of
DIFFERENT pathogens do. For every pair of the 260 endpoint columns, the Jaccard overlap of their
top-1000 highest-scoring compounds, aggregated per pathogen into a same-vs-different box pair.
Figure in `src/plots_matrix_analyses.py`; writes the 260x260 Jaccard matrices (reused downstream),
a per-pathogen summary CSV, and PNG + PDF to `output/08_pathogen_jaccard/`.
**Key decisions, for review:**
- **Cutoff: top 1000** of 1,355,109 compounds per column — user-directed.
- **Scope: the 15 curated pathogens of interest** (`config/pathogens_of_interest.csv`), replacing the
  former `min<K>`-endpoint thresholds on 2026-08-07 (user-directed). **254 of 300 columns, 15 of 57
  pathogens**; the other 42 pathogens / 46 columns are removed *entirely*, ceasing to be
  different-pathogen partners too, not merely losing their own box.
  **This narrows the comparator and a caption must say so:** a pathogen's crimson box is now against
  the other 14 priority pathogens only, **not** against all 57, so "specific to this pathogen" here
  means "relative to the other priority pathogens". The old `min5` figure compared against every
  pathogen that cleared the bar, including the gut-microbiome organisms.
  **Two of the 15 have a single endpoint** (*Campylobacter*, *H. pylori*), so no same-pathogen pair
  exists and they carry a crimson box only, with `same_median = NaN`. They are **kept, not dropped**:
  for a pathogen on the priority list, having too little in the hub to assess is itself the result,
  and hiding the row would hide it. The run prints them by name every time.
  Superseded values, for reference: `min5` was 249 of 300 columns / 12 of 57 pathogens, `min2` 259 /
  16, and at `min2` `bfragilis` topped the ranking at 0.5026 on a *single pair*.
- **Pathogen codes: the two configs were aligned at source (2026-08-07).**
  `config/08_endpoint_selection.csv` spelled two organisms `Campylobacter spp` / `Enterobacter spp`
  while `config/pathogens_of_interest.csv` said `Campylobacter` / `Enterobacter`. The name lookup
  therefore missed and `_pathogen_code` fell through to its mechanical fallback, coding them **`cspp`**
  and **`espp`** — which both mislabelled two of the 15 in this figure and would have made a filter on
  the config's `code` column silently drop them. The selection config was edited to match (7 rows), so
  the codes are now `campylobacter` and `enterobacter` everywhere and every lookup is an **exact
  match**. `default.PATHOGEN_ORGANISM_ALIASES`, which existed only to paper over this, was removed
  along with its use in `src/eval_predictor_performance.py`.
  **Do not reintroduce substring matching** if they ever diverge again: `Candida albicans` would
  capture *C. glabrata*, and `Streptococcus pneumoniae` would capture *S. parasanguinis* and
  *S. salivarius*. Fix the spelling instead. `pathogens_of_interest_nodes` **raises** if any of the 15
  fails to resolve, so a future rename fails loudly rather than quietly shrinking the figure.
- **Three matrix variants, not five — and this step is the only consumer of any of them.** The three are
  `baseline` (step 07's *named*), `zscore_l2rownorm` and `rankpct_l1rownorm`; step 07's plain *z-score*
  and *rank-percentile* variants are never used by anything. Top-N Jaccard depends only on each column's
  own internal ranking, and both column scalings are strictly increasing per column — so z-scoring and
  rank-percentiling cannot change any column's top-1000 set. The unscaled and both scaled matrices give
  the same result, computed once; only the row-normalized matrices change rankings, because row
  normalisation mixes values *across* a row. The script **asserts** this at runtime rather
  than assuming it, and the assertion earns its keep: baseline == rank-percentiled holds exactly, but
  baseline == z-scored comes back **False** — `(x - mean) / std` in float32 reorders near-tied values in
  exactly **one column of 300** (`lmajor__eos60mw__leishmania_mlp`, 155 of its top-1000 members shift),
  touching **156 of 90,000 cells**, max absolute difference 0.001032 and mean 0.000599 over those cells.
  Re-verified on the 300-endpoint rebuild: still that one column and no other, and
  baseline == rank-percentiled differs in **0 of 90,000** cells. (At 260 endpoints it was 138 of 67,600.)
  Note the difference floor is ~1/2000 = 5e-4, the granularity of a Jaccard over two ~1000-element sets,
  so any comparison reporting differences below that is reading float noise — e.g. from a CSV round-trip
  rather than the in-memory arrays the assertion uses.
- **Same-model pairs are included** (the literal "each column against all others"). This matters at
  pathogen level: **3 of the 15** (`enterobacter`, `spneumoniae`, `efaecium`) have
  `n_same_pairs_excl_same_model == 0` — every same-pathogen pair they have comes from one model's multiple
  output columns, so their box is that model agreeing with itself and says nothing about cross-model
  specificity. `enterobacter` tops the baseline ranking (0.1105) on exactly that basis, so **the top row
  of the figure is an artifact, not a finding** — the first row with genuine cross-model signal is
  `ecoli` at 0.0471 with 167 cross-model pairs. The summary CSV carries
  `same_median_excl_same_model` alongside `same_median`; read them together.
- **Linear x-axis** (user-directed). Values bunch near zero so the different-pathogen boxes render as thin
  slivers, but exact-zero pairs are shown rather than silently dropped by a log axis. Nothing is filtered.
- **Plain matplotlib, not the 3 cm-cell publication grid** — a deliberate, user-approved departure from
  `docs/figure_conventions.md`, since the per-pathogen endpoint and pair counts on the y-axis go illegible
  at page width. PNG *and* PDF are still written; the vector copy is the readable one.

## 09_reference_library_projection.py
Visualises where the reference library sits in chemical space, highlighting each pathogen's
highest-predicted-activity compounds. `eos1klk` (fetched by `00_download_data.py`) computes four 2D
projections of the library — PCA, UMAP, t-SNE, TMAP; each gets its own figure: a silver background of the
full library's density, with each of the 15 pathogens in `config/pathogens_of_interest.csv` getting its own
panel showing its `PROJECTION_TOP_N` highest-scoring compounds (by `consensus_score`) in crimson. This is a
**rank cutoff (top-N), never a score threshold** — no `consensus_score` cutoff value is chosen or reviewed,
only a compound count. Analysis in `src/eval_projection.py`, figures in `src/plots_projection.py`; writes to
`output/09_reference_library_projection/`.
**Memory:** only one pathogen's `key` + `consensus_score` columns (of each ~424 MB prediction file) are read
at a time, immediately reduced to its `PROJECTION_TOP_N` highest rows and discarded — no more than one
pathogen's raw scores, and never all 15 prediction files, are held in memory together. The score ranking
does not depend on projection method, so each pathogen's top-N is computed once and reused across all four.
**Key decisions, for review:**
- **Top N per pathogen** (`PROJECTION_TOP_N`, in `src/default.py`) is **1000** — user-directed.
- **Grid resolution** (`PROJECTION_BINS`) for the silver background density is **60x60** cells per method —
  sized to stay legible at the ~30x30 mm panel size in a 15-pathogen small-multiples grid, not a fitted
  value.

## 10_physchem_matrix.py
Builds the **physchem** block of the reference library from `config/physchem_models.csv` — all 22
descriptors of `eos4djh` (datamol-basic-descriptors) — one row per compound, plus a per-endpoint stats
CSV and a 22-panel distributions figure. Writes `output/10_physchem_matrix/`. Collection logic in
`src/eval_property_matrix.py` (shared with step 12), figure in `src/plots_property_matrix.py`.
**Split out of the former `11_additional_properties.py` (2026-08-06)**, which built the physchem and
cytotox blocks in one 46-column table. Each family now builds its own matrix, so a consumer opens only
the blocks it needs and the two kinds of number cannot be conflated. The superseded script is kept at
`tmp/superseded_scripts/11_additional_properties.py` — it was never committed, so there is no git
history to recover it from.
**Column naming:** `physchem__{model_id}__{column_name}`, the same three-part shape as the pathogen
matrix (step 07, prefix = pathogen code), the abx block (step 11) and the cytotox block (step 12), so
the families join column-wise on `key`.
**These are calculations, not predictions.** `eos4djh` wraps deterministic RDKit-family arithmetic, so
this is the one property block with no model error, no training set and no leakage dimension. A caption
must not describe it as "predicted", and steps 13/14 must not read it as interchangeable with the
predicted blocks.
**All 22 descriptors are kept even though 5 are near-redundant** (a signed-off request, not an
oversight). Verified on 200k rows: `n_rings`, `n_aliphatic_rings`, `n_aromatic_rings` and
`n_saturated_rings` are each **exactly** the sum of their two carbocycle/heterocycle components on 100%
of rows, and `n_radical_electrons` is 0 for all but 38 of 200,000 compounds (the run reports it as the
one endpoint non-zero in <1% of the library). Anything fitting a model or drawing a correlation matrix
off this block should drop those five first — the independent set is the other 17.
(`n_saturated_rings` is *not* a duplicate of `n_aliphatic_rings`; they differ on 34% of rows.)
**Observed descriptor ranges are the library's filter, not the model's.** MW 198.9–699.5, cLogP
−2.0–7.0, `n_rings` ≤ 6 — the Ersilia reference library is already drug-like-filtered, so these bounds
describe the input set. A caption must not read them as a property of `eos4djh`.
**Upstream spelling kept:** `n_aliphatic_heterocyles`, `n_aromatic_heterocyles` and
`n_saturated_heterocyles` are misspelled in Datamol itself. The config matches the real CSV header
rather than correcting it, or the lookup would miss.
**Figure: distributions, not a UMAP.** The abx and cytotox blocks are drawn as UMAP highlights because
"the top N most antibiotic-like compounds" is a set worth locating in chemical space; "the top N
compounds by molecular weight" is an arbitrary slice of a continuous descriptor. So this block is
summarised by each descriptor's full-library distribution — a 5-column small-multiples grid, **(5, 5)
cells**. Every endpoint gets a panel including the near-constant `n_radical_electrons`; a flat panel is
the honest rendering of a flat column.
**Un-normalized only**, matching the other two blocks: choosing a transform is a decision for after the
blocks are joined, and `10_physchem_endpoint_stats.csv` is what that decision should be made from.
**No rows dropped**, nothing imputed; 0 missing values across all 22 endpoints.

## 11_abx_resemblance_matrix.py
The antibiotic-resemblance score matrix and its endpoints on the library UMAP, from
`config/antibiotic_resemblance.csv` (55 selected endpoints across 4 models). **Renamed from
`12_antibiotic_resemblance_matrix.py` (2026-08-06)**; outputs moved to
`output/11_abx_resemblance_matrix/` with the `11_abx_*` prefix, and the figure key in
`figure_cells.json` became `11_umap_abx_endpoints_max1000`. Engine `src/eval_abx_matrix.py`, figure
`src/plots_abx_projection.py`.
**Column naming:** `abx__{model_id}__{column_name}` — a constant group code in the pathogen slot.
**Un-normalized only**, by request, for the same reason as the other two blocks.
**The highlight rule is not a rank cutoff.** 44 of the 55 endpoints are binary or small integer counts,
so each panel shows every compound with a value > 0, capped at `PROJECTION_TOP_N` — never padded. **19
endpoints hit that cap**, so their panels show an arbitrary key-ordered subset, not the full flagged
set; read each panel's `n_shown/n_nonzero` annotation.
**4 endpoints are constant zero library-wide** (`arsenic_cpds`, `b_lactamase_inhibitors`,
`lipopeptides`, `polypeptides`) and are omitted from the **figure only** — they remain in the matrix and
the stats CSV. 32 of 55 have fewer than 1000 non-zero compounds. 66 missing values, kept.
**Reuses step 09's `09_umap_background.csv`** rather than recomputing the density, so its panels are
directly comparable to the pathogen ones.

### Pathogen x abx overlap (merged in from the former step 14, 2026-08-06)
Outputs renamed `14_*` -> `11_*`; the 15 figure keys are now `11_umap_abx_overlap_{pathogen}`.
Merged because the abx side of the intersection is the highlights table this step writes a few lines
earlier — as a separate step it re-read it from disk. The superseded script is at
`tmp/superseded_scripts/14_pathogen_abx_overlap.py` (never committed), and
`output/14_pathogen_abx_overlap/` was deleted after its two CSVs were verified byte-identical to the
regenerated `11_*` ones.

**`figure_cells.json` had to become append-not-replace for this merge.** Step 11 now writes two figure
families into one dir, and all three `save_*_figures` entry points previously opened the manifest with
`"w"` from an empty dict — so whichever ran last truncated it to its own entries, silently losing the
abx grid's footprint. They now route through `plotting_utils.merge_figure_cells`. Verified: step 11's
manifest holds **16** entries (1 abx grid + 15 overlap).

Overlays each pathogen's predicted hits with each antibiotic-resemblance endpoint's compounds on the
`eos1klk` UMAP, on the same silver full-library background used by steps 10, 12 and 13. **One figure per
pathogen** (15), each a 3x3 grid with one panel per abx endpoint and three point groups: crimson =
pathogen only, cobalt = antibiotic-like only, lime = **both** (drawn last and larger, since the
intersection is the subject). Analysis in `src/eval_pathogen_abx_overlap.py`, figures in
`src/plots_pathogen_abx_overlap.py`; writes into `output/11_abx_resemblance_matrix/`.
**Pure set arithmetic — no new scoring and no threshold.** It reads two summary CSVs that already exist
(`10_top1000_per_pathogen.csv` and `12_abx_umap_highlights.csv`) and never touches a prediction file or a
full matrix.
**Endpoints** are `ABX_OVERLAP_ENDPOINTS` in `src/default.py` — 9 of the 55 selected in
`config/antibiotic_resemblance.csv`, user-directed, spanning one continuous learned score (`abx_score`),
two AntibioticDB similarity counts, and six substructure/class flags. Two were requested under names that
do not exist in the config and were corrected to its spelling: `betalactan_motif` -> `betalactam_motif`,
`b_lactam_all` -> `b_lactams_all`.
**The two sides are NOT selected the same way, and a panel cannot be read without knowing this.** The
pathogen side is a rank cutoff — the top `PROJECTION_TOP_N` = 1000 by `consensus_score`, so every pathogen
contributes exactly 1000. The abx side is not: only `abx_score` is continuous, so step 11's rule (every
compound with value > 0, capped at 1000) applies. Three of the nine have fewer positives library-wide than
the cap and are drawn exhaustively (`carbepenem_motif` 346, `ansamycins_rifamycins_macrolides` 577,
`b_lactams_all` 733); the other six hit the cap and are an **arbitrary** subset of their positives.
Panels mark a capped endpoint with a trailing `*` on the `n=` annotation, and `14_overlap_counts.csv`
carries `n_abx`/`abx_capped` per row — so an intersection of zero can be read as "no overlap" rather than
confused with "this endpoint barely has any hits at all".
**`14_overlap_counts.csv`** holds `n_pathogen`, `n_abx`, `n_both` and the Jaccard index for all 135
(pathogen x endpoint) pairs — the summary to read before interpreting any single panel.

## 12_cytotox_matrix.py
Builds the **cytotox** block from `config/cytotoxicity_models.csv` (24 selected endpoints across 4
models: `eos42ez`, `eos7m30`, `eos3le9`, `eos3dys`), plus a per-endpoint stats CSV. Writes
`output/12_cytotox_matrix/`. Shares `src/eval_property_matrix.py` with step 10.
**Split out of the former `11_additional_properties.py`** alongside step 10, same rationale.
**Every column here is a prediction**, unlike step 10's deterministic descriptors. The two must not be
treated interchangeably.
**The cytotox figure is the toxicity projection below**, merged in from the former step 13 — it needs
this matrix plus step 09's UMAP background. This step also supplies the matrix that steps 13 and 14
read.
**`eos7m30`'s own 8 physicochemical columns stay `No`** in the config, so `eos4djh` remains the single
physchem source and no two blocks can disagree about molecular weight or logP.
**Carrying the model ID matters here:** two models score HepG2 —
`cytotox__eos42ez__cytotoxicity_hepg2` vs `cytotox__eos3le9__ic50_hepg2_72h_5um` — and the column name
is what keeps them apart.
**Endpoint selection is not made here** — it is the manually curated `selected` column. No threshold or
cutoff is applied; raw outputs pass through unchanged. 0 missing values across all 24 endpoints.
**Alignment:** prediction files are concatenated column-wise only after their `key` order is verified
against the first file; any file whose order differs is reindexed on `key` rather than concatenated
blindly.

### Toxicity projection (merged in from the former step 13, 2026-08-06)
Outputs renamed `13_*` -> `12_*`; the figure key in `figure_cells.json` is now
`12_umap_top1000_toxicity`. Merged because it consumes exactly this step's matrix and nothing else.
The superseded script is at `tmp/superseded_scripts/13_toxicity_projection.py` (never committed, so
no git history to recover it from), and `output/13_toxicity_projection/` was deleted after its four
CSVs were verified byte-identical to the regenerated `12_*` ones.

The toxicity counterpart of step 10, on the same `eos1klk` chemical-space layout: a silver
full-library density background with, in crimson, each toxicity endpoint's most toxic compounds. One
small-multiples figure (`13_umap_top1000_toxicity`), 24 panels — one per `selected == "Yes"` endpoint in
`config/cytotoxicity_models.csv` — scored from the step-12 cytotox matrix. Analysis in
`src/eval_tox_projection.py`, figure in `src/plots_tox_projection.py`; writes to
`output/12_cytotox_matrix/`. Panels are labelled with both the endpoint and its model ID, since
several endpoints read overlapping biology from different models (two HepG2 readouts, three cytotoxicity
readouts).
**Top N per endpoint** is `PROJECTION_TOP_N` = **1000**, shared with step 09 — a **rank cutoff, never a
score threshold**, so no score value is chosen or reviewed, only a compound count.
**Projection method:** only UMAP (`TOX_PROJECTION_METHOD`) of the four `PROJECTION_METHODS` is drawn —
24 endpoints x 4 methods would be 96 panels, and UMAP is the layout step 09's pathogen figures are read
from, so the two are directly comparable. The engine takes the method as an argument if others are wanted.
**Direction — which end is "most toxic":** all 24 endpoints are ranked **highest-first**
(`TOX_RANK_DESCENDING`). For the 23 classification endpoints the score is the predicted probability of
the toxic/active class, so this is direct. `ld50_zhu` is the one regression endpoint and the one genuine
ambiguity: eos7m30 emits it in `log(1/(mol/kg))`, where the reciprocal inverts the dose scale, so a
higher value is a *lower* LD50 and therefore more acutely toxic. Note that the eos7m30 column metadata
labels it `direction: low`, which describes the underlying LD50 dose rather than the transformed value
actually emitted — **the two disagree, and this was resolved in favour of highest-first by user sign-off**,
corroborated by its positive Spearman correlation (+0.21 to +0.30) with all six independently-trained
cytotoxicity models (eos42ez x3, eos3le9 x2, eos3dys) on a 150k-compound sample.
**Score ranges of each top 1000 are printed to stdout** — worth reading before interpreting a panel, as
the endpoints are very unevenly saturated: `dili`'s top 1000 spans 0.996-0.999 and `sr_are`'s 0.964-1.000
(the cutoff separates almost nothing), whereas `nr_ar_lbd`'s spans 0.393-0.863.
**Memory:** the ~720 MB step-12 cytotox matrix is streamed in chunks reading `key` + the 24 score columns only
(never the `input` SMILES column, which is most of the file), with each endpoint's running top-N reduced
after every chunk.

## 13_curated_predictions.py
Asks whether the property/resemblance columns carry any signal about pathogen activity. Treats every
column of the step-10 (physchem), step-11 (abx) and step-12 (cytotox) blocks as a **predictor**, and every curated
activity endpoint as a binary **target**, giving one performance value per (predictor, target) pair —
**101 x 300 = 30,300** as re-run on 2026-08-06 (was 101 x 260 = 26,260 at 260 endpoints). Predictors by
family: abx 55, cytotox 24, physchem 22 — now read from three separate per-family matrices (steps 10,
11, 12) rather than two. Three figures, one per predictor family, each a box-with-jitter per predictor over its
distribution across all 260 targets, sorted by median, with a chance line at 0.5. Analysis in
`src/eval_predictor_performance.py`, figures in `src/plots_predictor_performance.py`; writes to
`output/13_curated_predictions/`.
**This is a descriptive association measure, NOT a trained model** — nothing is fitted, nothing is split,
no random seed is involved. Every value is a rank statistic over the full library.
**Target binarization:** `ACTIVITY_BINARIZE_TOP_N` = **1000** highest-scoring compounds = positive class,
the remaining 1,354,109 = negative. A user-directed **rank cutoff on a fixed count, never a score
threshold**, so prevalence is a constant 0.0738% across all targets. All 260 selected endpoints are
asserted to be `direction == higher` — the engine raises rather than silently inverting a `lower` one.
**Metric choice** is driven by the predictor's own value type, resolved on the **full column, never a
subsample**: continuous -> AUROC, binary -> balanced accuracy. Both share a 0.5 chance baseline, which is
the only reason they can share a y-axis; box colour encodes which is which. The full-column rule is not
pedantry — `n_radical_electrons` is 0 for all but ~1 compound in 5,000 and classifies as *binary* on any
sample, but as *continuous* (correctly) on the full column.
**AUROC is reported RAW and may fall below 0.5**, so anti-correlation stays visible rather than being
folded to its mirror image. This is load-bearing: `qed` has a median AUROC of **0.31**, i.e. it is a
strong *inverse* predictor of activity, which folding would have disguised as moderate signal.
**AUROC is computed by the Mann-Whitney rank-sum identity**, not `sklearn.roc_auc_score`: the predictor is
ranked once (ties averaged, so the result is exact) and each target is then a 1000-element gather. 30,300
direct sklearn calls over 1.35M rows would run for hours. Verified identical to sklearn (<= 1e-16) on ten
real (predictor, target) pairs plus synthetic no-tie/heavy-tie cases.
**Missing values — nothing imputed, nothing dropped library-wide.** Two independent cases, both signed
off and both recorded in the outputs: (1) one *target*, `eos4zfy:maip_score`, has 15 unscored compounds —
an unscored compound cannot be claimed to be in the top 1000, so it is ineligible for the positive class
and stays negative; (2) eleven *predictors* (all `eos6ojg`) share the same 6 unscored compounds and are
evaluated **pairwise-complete** over 1,355,103, with `n_compounds`/`n_unscored` in the summary CSV.
Filling 0 was rejected — these are similarity counts, where 0 asserts "no similar antibiotic found".
**Degenerate predictors give NaN, never 0.5:** four `eos19mt` flags (`arsenic_cpds`,
`b_lactamase_inhibitors`, `lipopeptides`, `polypeptides`) are all-zero library-wide, so balanced accuracy
is undefined; they account for all 1,040 NaN values and appear as empty slots on the abx figure rather
than as chance-level boxes.
**`same_model` flag:** `eos3dys` supplies 2 cytotoxicity predictors *and* 20 of the 260 targets, so those
40 pairs share training data and are not independent evidence. They are kept and flagged in
`13_predictor_performance.csv` (exactly 40 rows) for exclusion downstream.
**Fourth figure — the activity endpoints as predictors of each other**
(`13_performance_activity_by_organism`, from `13_activity_self_performance.csv`). Same machinery with
activity on both sides: each endpoint's **raw, un-binarized** score as the predictor against every
endpoint's top-1000 binarized version as the target — 260 x 260 = 67,600 AUROCs. All activity endpoints
are continuous, so AUROC applies throughout with no metric selection. The x-axis groups by the
**predictor** endpoint's organism (56 organisms), and points are coloured by whether the target endpoint
belongs to the same organism.
**Self-pairs are excluded from the figure**, kept in the CSV under a `self_pair` flag. An endpoint against
its own binarization is 1.0 by construction — the top 1000 of a score *are* its 1000 highest values — so
it measures nothing and would add one guaranteed-perfect point to every box. It is retained as a
correctness check: the 260 self-pairs come back at 0.999952-1.000000, the sub-1.0 values being score ties
at the top-1000 boundary (tied compounds straddling the cutoff), not an error.
**Same-organism vs cross-organism is the distinction the figure exists to make** (6,334 vs 61,006 pairs
after removing self-pairs). Two endpoints of one pathogen agreeing is model self-consistency; an endpoint
scoring highly against *unrelated* pathogens is a different claim entirely. Note that the **41 organisms
with a single selected endpoint have no same-organism pairs at all** once self-pairs are removed, so
their boxes are cross-organism only and are not comparable to the multi-endpoint organisms on that axis.
**Point subsampling:** each organism box pools (its endpoints x 259 targets), reaching ~16,500 points for
*P. falciparum*, so the jittered overlay is capped at 400 points per colour, seeded with `RANDOM_SEED`.
The box itself is computed from **all** values, not the subsample.
**Fifth figure — the 15 pathogens of interest, consensus models collapsed**
(`13_performance_pathogen_subset`, from `13_pathogen_subset_self_performance.csv`). The same
activity-vs-activity block restricted to `config/pathogens_of_interest.csv` and with each ChEMBL
antimicrobial model reduced to its single `consensus_score` column: **215 endpoints -> 59** (12 consensus,
47 single) across all 15 pathogens, so 59 x 59 = 3,481 pairs. One box per endpoint, ordered by pathogen
(pathogens by median, endpoints by median within), labelled `{pathogen} - {endpoint}`.
**Why collapse:** a single ChEMBL model contributes up to 52 highly correlated sub-endpoints (one per
source assay) — *P. falciparum*'s eos4an7 alone had 52 of the 260 — which dominate any pooled view by
count alone. Where a model publishes a consensus column that column *is* its headline score; where it
never had more than one sub-model there was nothing to take a consensus over (`eos7iak`, `eos9eyo`,
`eos5qya`), so its single endpoint is kept as-is.
**The rule is applied per (model_id, organism), not per model** — `eos3dys` spans six organisms and has
no consensus column, so each of its organisms keeps its own endpoints.
**Organism matching is an explicit alias map** (`PATHOGEN_ORGANISM_ALIASES`), never a genus substring:
the two configs spelled `Campylobacter`/`Campylobacter spp` and `Enterobacter`/`Enterobacter spp` (aligned 2026-08-07)
differently, while substring matching would wrongly capture *C. glabrata* for *C. albicans* and
*S. parasanguinis*/*S. salivarius* for *S. pneumoniae* — all distinct organisms in the curation.
**Colour is a secondary cue only.** Every box is identified on the axis, because
`plotting_colors.distinct_colors` supplies only 9 substantive hues and explicitly warns that a
15-swatch legend is unreadable at panel size. Same-organism points are drawn larger and more opaque than
cross-organism ones, since colour is already carrying pathogen identity here.
**Panel height is 4 cells, not 2** (the other step-13 panels): the two-part tick labels need roughly half
the panel height at 59 categories, and at 2 cells the boxes collapsed into a strip along the top.

**Sixth figure — the curated 12 predictors, three families in one panel**
(`13_performance_curated_predictors`, from `13_curated_predictor_performance.csv`). A user-directed
shortlist of the 101 property columns (`CURATED_PREDICTORS` in `src/default.py`) scored against the same
59 consensus-collapsed pathogen endpoints as the fifth figure — 12 x 59 = 708 pairs. One box per
predictor, families in contiguous blocks and sorted by median within each, coloured physchem /
cytotox / abx.
**Predictors:** physchem `mw, tpsa, clogp`; cytotox `cytotoxicity_hepg2, cytotoxicity_hskmc,
cytotoxicity_imr90, cytotoxicity_ic50, ic50_hepg2_72h_5um, ic50_hepg2_72h_10um`; abx `abx_score,
num_sim_0_5_all, num_sim_0_5_subset`. Declared as bare column names and resolved against the
`{family}__{model_id}__{column_name}` names at read time, so a model version bump needs no edit here.
**All twelve are continuous, so the whole panel is one AUROC scale** — the engine asserts this rather
than assuming it, and raises if a binary column is ever added to the shortlist. This is what the
per-family figures cannot do: their abx block mixes AUROC and balanced accuracy, which must not be
pooled on one axis. It is also why colour can be a primary encoding here (3 families, a legible
3-entry legend) rather than the secondary cue it is forced to be at 15 pathogen groups.
**`same_model` flag:** `cytotox__eos3dys__cytotoxicity_ic50` shares a model with the eos3dys activity
targets, giving 16 non-independent pairs, flagged in the CSV.

### Config changes: *M. tuberculosis* endpoint selection

*M. tuberculosis* is one of the 15 pathogens of interest but at the last commit had **all 35 of its
`eos43d6` columns marked `No`**, plus every other MTB bioactivity endpoint — so it was absent from every
activity figure. It was also the only ChEMBL pathogen model not fully selected (1/35, against 52/52 for
`eos4an7`, 19/19 for `eos8lcw`, 13/13 for `eos5eya` and so on). Fixed in two user-directed steps:

1. `eos43d6,...,consensus_score` flipped `No -> Yes`, taking the selection **260 -> 261**. A minimum fix
   to get MTB into the figures at all.
2. **All remaining MTB `bioactivity` endpoints selected, with two exceptions** (2026-08-06),
   taking the selection **261 -> 300**. MTB now carries **40 endpoints**, second only to *P. falciparum*
   (64). Per model:

| model | endpoints | selected | rationale |
|---|---|---|---|
| `eos43d6` antimicrobial-activity-mtuberculosis | 35 | **35/35** | Now consistent with every other pathogen model, which are all 100% selected. |
| `eos9ivc` anti-mtb-seattle | 3 | **3/3** | Whole-cell MIC50/MIC90/WCS endpoints. |
| `eos46ev` chemtb | 1 | **1/1** | Whole-cell activity probability. |
| `eos7kpb` h3d-virtual-screening-cascade-light | 2 | **1/2** — `mtb_norm` only | The normalised column is kept, the raw `mtb` dropped. |
| `eos24jm` qcrb-tb | 1 | **0/1** | **Target-based** (QcrB inhibition), not whole-cell activity — a different claim from every other endpoint in the block. |

**`eos7kpb` now uses opposite conventions across organisms, and this is unresolved.** For MTB the
normalised column is selected and the raw one is not; for *P. falciparum* the raw `pf_k1` and `pf_nf54`
are selected and their `_norm` twins are not. Same model, same config, two rules. Either the MTB rows or
the *P. falciparum* rows should change so the model is read one way — flagged, not silently fixed, since
which convention is right is a scientific call.

**Permeability endpoints stay `No`, and that is not an MTB-specific decision.** MTB's 9 permeability
endpoints (`eos1lb5` x6, `eos3ujl`, `eos5jv3`, `eos8d8a`) are excluded under the config-wide rule that
only `bioactivity` rows are ever selected: permeability is 0/22 across all organisms, as are ADME 0/14,
toxicity 0/6, class_prediction 0/8, structural_alert 0/7 and physicochemical 0/2.

**`eos9ivc` (anti-mtb-seattle) arrived after the step-07 cache was built.** Its predictions were
missing when the cache was written and were downloaded on 2026-08-06 as
`annotation_preds_ref_library/eos9ivc_v2.csv` (1,355,109 rows, all 3 columns continuous in [0, 1], 0
missing). It needs no registration beyond this config — *M. tuberculosis* is already in
`config/pathogens_of_interest.csv`, so its columns name themselves `mtuberculosis__eos9ivc__*`.

**Downstream effect — step 07 DOES need a rebuild, and 08/09/15/16 after it.** The cache's
"contains every referenced endpoint, `Yes` and `No`" guarantee only covers models whose prediction
files existed at build time. Verified against `07_score_matrix_full.parquet` (393 columns): of the 300
now-selected endpoints, **297 are cached and 3 are missing — exactly `eos9ivc`'s**
`mic50_10um`, `mic90_10um`, `wcs_70percent`. So the eos43d6 flip alone would have needed no rebuild, but
adding eos9ivc does. Step 07 skips on file existence, not on content, so it will **not** notice by
itself: the 1.5 GB parquet and the five CSVs (25 GB total, ~2 h to regenerate) have to be removed for it
to rebuild.

**Obsolete as of the 2026-08-06 rebuild:** steps 07-14 were all re-run at 300 endpoints, so the counts
in their sections are current unless explicitly labelled historical. Any bare `260` / `67,600` /
`26,260` still in this file is a superseded value kept for comparison, not a live one.

## 14_auroc_matrix.py
Collapses each organism's activity endpoints into **one score per organism**, then draws the AUROC
matrix: **15 organism rows x 27 columns** (the same 15 organism aggregates, then cytotoxicity 6, abx
resemblance 3, physchem 3 = 405 cells). Each cell carries its printed AUROC. Assembly and computation in
`src/eval_auroc_matrix.py`, figure in `src/plots_auroc_matrix.py`; writes to `output/14_auroc_matrix/`.
Replaces an earlier per-endpoint 61 x 73 version, which was too dense to read: a pathogen contributed up
to 13 correlated endpoints, so its band said more about how many assays it has than about the pathogen.
**The merge:** an organism's endpoint columns are scaled to **percentile rank** within the full library
(`ORGANISM_MERGE_METHOD = "rank_pct"`) and then **averaged** (`ORGANISM_MERGE_AGG = "mean"`). Scaling
*before* aggregating is the point — raw scores from different models sit on unrelated ranges, and
averaging them directly would weight whichever endpoint has the widest one. Endpoints are the
consensus-collapsed 15-pathogen set (61 endpoints; 1-13 per organism).
**Two properties of the merge, recorded rather than corrected** — both follow directly from merging the
endpoints as selected: (1) **five organisms have exactly one endpoint** (Campylobacter, Enterobacter,
*E. faecium*, *H. pylori*, *S. pneumoniae*), so nothing is merged and their score IS that endpoint's
percentile rank — their row is not the same kind of quantity as *E. coli*'s 11-endpoint mean; (2) a ChEMBL
`consensus_score` is itself an aggregate over sub-models, so averaging it with individual assay endpoints
gives it equal weight to a single assay.
**Source:** step 07's **parquet cache**, re-scaled here — not `07_score_matrix_named_rankpct.csv`, which
is stale (260 columns, predating two config changes) and lacks several of these endpoints.
Step 07's own mean-rank section rebuilds from the parquet for the same reason.
**Colour: discrete 0.1-wide bins from 0.2 to 1.0, on a DIVERGING scale pinned to chance**
(`DivergingColormap("crimson_cobalt")`, reversed to cool-low / warm-high). Cool = below 0.5, near-white
at 0.5, warm = above. A diverging map is the right family here **because 0.5 is a real neutral with data
on both sides**: 44 of the 405 cells (10.9%) fall below chance, down to 0.262, meaning those predictors
rank actives *below* inactives — signal with a direction, not absence of signal. It would be the wrong
choice for a quantity with no natural centre, where the pale middle would land on an arbitrary value.
**The arms are unequal** (0.24 below chance, 0.50 above), so bin midpoints are mapped through a
two-slope normalisation that pins 0.5 to the colormap's centre. Sampling linearly across the whole range
instead would put the white point at 0.6 and quietly assert that 0.6 is chance.
**Nothing is clipped** — the 0.2 floor sits below the matrix minimum of 0.262. This reverses an earlier
clip-at-0.5 setting, which flattened all 44 below-chance cells into one grey bin.
*Colour history, since the choice was iterated:* `SpectralColormap("npg")` was tried first and rejected —
its hues carry no rank, so magenta (1.0) vs blue (0.8) vs red (0.6) gave a reader no way to tell high
from low, and red read as "hot" while sitting near the bottom. A monotonic `ContinuousColormap("plum")`
fixed the ordering but was visually poor. The diverging scale solves both and adds the below-chance arm.
**Gridlines are off.** stylia's article style draws a grid; over a filled heatmap it lands on top of the
cells and strikes through the printed values.
**Peripheral tracks:** block category, organism class and predictor model on the top edge; organism class
on the left. No swatch legend (user-directed) — the block track labels itself in place and the organisms
are named on both axes, but note the narrow class bands (Fungi, Protozoa, Helminths, Mycobacteria) are
unkeyed colour only. Track axes **share x/y with the grid**, which is what guarantees cell alignment.
### Second view: shared actives (`14_overlap_matrix`)
The same 15 x 27 axes and annotation tracks, a different quantity: **how many of the row organism's
1000 actives fall in the column's own top 1000** (0-1000). AUROC asks "does this predictor RANK the
row's actives highly across the whole library"; this asks how many of the very same molecules it puts
at the top. A predictor can do the first well without doing the second. Values in
`14_overlap_matrix.csv`.
**A raw count, not Jaccard.** Both sets have exactly 1000 members, so `J = i / (2000 - i)` is a
monotone re-expression of the same number and orders the matrix identically — the count is the one a
reader can act on ("724 of the 1000 shared"). **The measure is symmetric**, so the bioactivity block is
a symmetric matrix, unlike AUROC's.
**Bins are NON-UNIFORM** (`OVERLAP_MATRIX_BINS` = 0, 1, 10, 25, 50, 100, 200, 400, 750), tuned to a
heavily skewed distribution: off-diagonal counts run 0-724 with a **median of 3**, and two random
1000-compound sets out of 1,355,109 would share ~0.7 by chance. The boundary at 1 is the informative
one — it separates the 33% of cells sharing **no** compound at all from those sharing some. No bin
holds more than a third. The colourbar uses **uniform** spacing so the narrow low bins stay readable.
**Colour is SEQUENTIAL cobalt**, not diverging: a count's neutral is 0, at the end of the scale rather
than in the middle, so there is nothing to diverge around. Matches the overlap heatmaps in step 08 and
the EU OpenScreen validation.

### Both figures carry a mean row, and both blank their diagonal
**Mean row.** A separate one-row band under each grid gives the **column-wise mean**, i.e. how a
predictor does averaged over all 15 organisms. It is kept OUTSIDE the main axes rather than added as a
16th row, so it can never be mistaken for another organism, and it shares the x axis so it stays
column-aligned. It also carries the column tick labels, since appending a band below the grid occupies
exactly where those labels would otherwise sit.
**The diagonal is excluded from the mean.** Those cells are an entity against itself — 1.00 AUROC and
1000/1000 shared, by construction — so including them would inflate every bioactivity column by a
guaranteed maximum.
**The diagonal is blanked in both grids**, drawn as a **dashed cell**, which distinguishes "deliberately
not shown" from "no value". For the overlap matrix this also frees the scale: the largest real overlap
is 724, and keeping the diagonal at 1000 would spend colour range on a cell that cannot be anything
else. Follows the convention in step 08 and `ActiveOverlapHeatmapPlot`. **The true values, diagonals
included, are still in both CSVs** — blanking is a display choice, and the build-time diagonal assert
still runs on the unmodified matrix.

**Three build-time checks**, each of which would otherwise fail silently and still look plausible: every
aggregate column's mean must be ~0.5 (a mean of percentile ranks must be; catches scaling along the wrong
axis), the diagonal must be 1.0 across all 15 (an organism against its own binarization; doubles as a
transpose test), and no cell may be missing. The script exits rather than drawing if any fails.
**Nothing is excluded any more.** This previously dropped 3 `eos9ivc` *M. tuberculosis* endpoints that
were absent from a step-07 cache built before their prediction file was staged; the rebuilt parquet
carries all three, and step 14 now reports 405 cells with none missing. The guard that caused it is
still live and still worth knowing about: step 14 intersects the selection with the columns actually in
the cache (`available=cached`) and **silently proceeds with fewer endpoints** rather than failing, so a
stale cache shows up only as a lower endpoint count in the log.

### Rebuild status

`config/08_endpoint_selection.csv` selects **300** endpoints (up from 260: all 35 eos43d6
*M. tuberculosis* rows plus eos46ev/eos7kpb/eos9ivc).

**Steps 07-14 were deleted and rebuilt from scratch on 2026-08-06** (27 GB -> 3.1 GB, 11 m 20 s,
all eight exit 0). Every step's own internal assertion passed and every value reproduced the
incremental runs, so the from-scratch result is consistent, not merely successful:

| step | time | verification |
|---|---|---|
| 07 | 4m32s | parquet 1,355,109 x 395, **0 models skipped**; mean rank 0.5000 / SD 0.1289 |
| 08 | 2m07s | 3 Jaccard matrices at 300x300 (re-run 2026-08-07: 15 pathogens of interest, 254/300 columns) |
| 09 | 35s | 4 background grids; 15 pathogens at exactly 1000/1000 |
| 10 | 18s | 1,355,109 x 22, 0 missing |
| 11 | 62s | 55 endpoints; **16 `figure_cells` entries** (1 grid + 15 overlap) |
| 12 | 43s | 24 endpoints + merged toxicity projection |
| 13 | 89s | key order verified across 4 sources; 300 targets x 101 predictors = 30,300 pairs |
| 14 | 34s | **diagonal exactly 1.0**, organism-score column means 0.5000, 15 x 27 |

Two things only a from-scratch run surfaced:
- **`13_excluded_targets.csv` is correctly absent** — it is written only when a selected endpoint is
  missing from the parquet, which was true when the cache lacked `eos9ivc`. The file had survived
  every rename since as a stale artifact. Step 13 writes **5** CSVs, not 6.
- **1,200 of the 30,300 pairs are undefined** — the 4 constant-zero abx endpoints x 300 targets. A
  constant column has no ranking, so no AUROC exists. Recorded, never imputed.

**Not regenerated, by decision:** the five scaled matrix CSVs. Nothing reads them (only 3 of the 5
would ever be consumed, all by step 08, which derives them in memory instead) and they cost ~1 h 54 m.
Run `python 07_score_matrices.py --write-matrix-csvs` if they are wanted as exports.

**Still stale:** step 01 only, for an unrelated reason — the metadata refresh (Ready 208 -> 214).
