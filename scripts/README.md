# Scripts

Numbered scripts run in order. Figure and plotting conventions (sizing, formats, per-figure
layouts) live in [`docs/figure_conventions.md`](../docs/figure_conventions.md).

## 00_download_data.py
Stages all input data in four sections: companion repos / eosvc (EU OpenScreen tasks, CoAdd data,
ChEMBL model reports and curation summaries), public GitHub files (Ersilia reference library,
DrugBank), Airtable model metadata, and Isaura precalc predictions for Ready annotation models.
Skip-if-exists, so it is safe to re-run. `--skip-isaura` skips the slow Isaura section; `--eosvc`
pulls Section 1 from eosvc instead of the companion repos.
**Key decision:** ChEMBL model performance metrics are recomputed downstream from the per-fold
`09_reports`; the staged `10_reports/` CSVs are kept only for fields that cannot be reconstructed
(quality weights, decision cutoffs, discarded-model reasons).

## 01_ersilia_metadata.py
Analyses the Airtable model metadata (Ready models only) and renders 14 panels — the task/subtask
breakdown, source type and output composition, biomedical area, target organism, licensing, tags,
container metrics and two pathogen treemaps.
**Shown vs counted:** value counts for every counted field are written to `*_counts.csv`, but not all
are plotted — Publication Type is counted only. Biomedical Area and Target Organism are capped at the
top 10 categories in the figure (full counts remain in the CSVs).

**Stacked panels:** Source Type and Output are each drawn as one stacked bar panel segmented by
Subtask, so a single panel carries the joint distribution instead of two. Segments use the subtask
palette (shades of the parent task's hue), never a palette of their own — the field is already encoded
by the bar position. Cross-tabs are written as `<field>_by_subtask_counts.csv`.

**Task hues:** Annotation = crimson, Representation = amber, Sampling = lime (`TASK_HUES`, the one
place to change them). These are the same three hues as `SOURCE_TYPE_COLORS`, which now survives only
for the two pathogen panels' dots/cells — the two never carry colour in the same panel.

**Task/subtask alternatives** — two ways to draw the same 208 models, on the same colours:
- `task_subtask` — bars. Best for reading counts. No legend: each bar is labelled next to its own
  colour, which makes this panel the subtask colour key.
- `task_subtask_waffle` — one square per model, 16 × 13 = 208 exactly. Shows n rather than stating it;
  reading 52 vs 39 means counting. Carries its own legend (with counts), so it is self-contained and
  can serve as the key too. That legend band is also what squares the panel — the grid alone crops to
  1.23, with the legend it lands at 1.05.

**Legends:** the two `*_by_subtask` stacks carry none (no room at 60 × 30 mm, and they are meant to
sit beside a key panel), so if one goes into the figure, `task_subtask` or `task_subtask_waffle` must
travel with it. The waffle and the `license` bar have their own.

**Biomedical Area / Target Organism colours:** every *named* area and organism consists purely of
Annotation models (`Any` is where all the Representation and Sampling models sit), so the named bars
are drawn in the **Annotation hue** — no separate palette — and `Any` in **silver**, the reserved
neutral for a catch-all. Neither panel carries a key, so the caption must explain the two colours.

**Known upstream mis-annotation:** `eos93h2` (`image-mol-gpcr`) is recorded in Airtable as
Representation / Featurization with Target Organism = *Homo sapiens*, but it should be an Annotation
model — being fixed at source. It is the only thing that breaks the "named organisms are
Annotation-only" rule, and it is already drawn as an Annotation row, so these two panels need no
change when the fix lands. When it does, **Task and Subtask must be changed together** — see the
warning in `docs/figure_conventions.md`, since a Task-only edit would leave the subtask-coloured
panels disagreeing with the task-coloured ones without erroring.

**License panel.** Bars use *simplified* identifiers: `-or-later` / `-only` are collapsed, so
GPL-3.0-or-later (51) + GPL-3.0-only (20) = one **GPL-3.0 bar (71), which overtakes MIT (69)**. This
loses a real legal distinction, so the ungrouped `license_counts.csv` is kept next to
`license_grouped_counts.csv`. The 27 models with no license are shown as **Not recorded**, not
dropped. Colour = reuse class: Permissive 104 (turquoise), Copyleft 76 (periwinkle), Non-commercial 1
(**fuchsia**, so a single-model bar is still visible), Not recorded 27 (silver).
**That split is a coarse classification, not legal advice** —
CC-BY-4.0 and CC0-1.0 are counted permissive despite not being OSI-approved for software. The 4
AGPL-3.0 models are all the same upstream project (`molgrad`), which does declare AGPL-3.0 on GitHub.

**Tag cloud.** All 59 tags from the free-text `Tag` field, sized by model count. **The only raster
panel in the repo** — `wordcloud` renders a bitmap, so this PDF embeds an image rather than vector
text. **Size is ordinal only:** frequencies are fed in as sqrt(count) and `relative_scaling=0.5`, so
the 1–49 count range compresses to 5.5–26.8 pt, and a word's area also depends on its character count
(*Gram-negative bacteria*, 13 models, covers more panel than *Descriptor*, 34). Read ranks coarsely;
exact counts are in `tag_counts.csv`. **Its 90 mm size is a legibility requirement, not a layout
choice** — at 60 mm the ten single-model tags fall to 2.9 pt. The script prints the measured smallest
type size and how many tags fall below the 5 pt print floor on every run; keep that at 0.

**Technical panels.** `runtime_100` and `image_size` are boxes-with-dots per Task on a log y axis, at
**45 × 45 mm** — the smallest square in the figure, which needs the category labels rotated 30° with
`n` inline (three horizontal labels need the full 60 mm). `docker_architecture` is a two-slice pie, also
45 × 45 mm (**AMD64 + ARM64 129 / 62%** turquoise vs **AMD64 only 79 / 38%** periwinkle — there is no
ARM-only build). Per-task quartiles go to `technical_metrics_summary.csv`.

**License pie.** `license_class_pie` is the compositional counterpart to the `license` bar: Permissive
104 (50%) turquoise, Copyleft 76 (37%) periwinkle, Not recorded 27 (13%) silver, Non-commercial 1
(<1%) **fuchsia** — the one use of that hue in the repo, because a 1-of-208 slice in any calmer colour
is invisible. It shows the four **reuse classes, not the ten licences**: four licences have exactly one
model each, a 1.7° wedge that cannot be seen or labelled, so per-licence detail stays in the bar
panel.

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

The architecture pie is a **snapshot, not a trend** (45% dual-arch among 2021 models vs 77% among
2026 — the 62% is accumulated stock, not current practice), so a caption needs the metadata snapshot
date.

**Panel sizes:** `task_subtask`, the waffle and the `license` bar are 60 × 60 mm (`cells=(2,2)`); the
two box panels and the two pies are 45 × 45 mm; the two ten-category fields are 45.75 × 45.75 mm (a
quarter of the page width, two side by side); the tag cloud and the pathogen panels are 90 × 90 mm —
on a 183 × 170 mm page. The `*_by_subtask` stacks are 60 mm wide and shorter: 30 mm for Source Type's 3
bars, 34.5 mm for Output's 4. That leaves Output's bars ~7% thinner (4.53 vs 4.88 mm) — an exact
thickness match needs 36.1 mm, a chosen trade-off documented in `docs/figure_conventions.md`.

**`pathogen_dataset_sizes` (full page width).** One box per pathogen over its modelled dataset sizes
(log axis, 96–334,766 compounds), with a small pie per dataset whose filled share is that dataset's
active fraction. **Added negatives are excluded from both the size and the ratio** — size is
`n_compounds − n_added_negatives − n_added_decoys`, which reproduces the curation pipeline's own
`n_mol_after`/`ar_after` exactly. Only 54 of 193 models had negatives added, so the medians barely
move, but for those 54 it matters a lot (`mtuberculosis/DR_0012`: 2450 compounds at 0.50 active →
1411 at 0.87). Decoys are zero for every model. Derived table: `dataset_sizes.csv`.
**Read the pies as a gestalt, not a measurement** — at ~2.5 mm, "mostly active" vs "mostly inactive"
is clear but 45% vs 50% is not; exact ratios are in the CSV. 54 of 193 datasets are majority-active
once added negatives come out. Full page width is a legibility requirement: *P. falciparum* alone has
51 datasets in one column.

**`pathogen_activity_ratios` (alternative to the above).** Same data, encodings swapped: active
fraction on the y axis where it can be resolved, dataset size as **dot area** (affine in √size, keyed
at 100/1,000/10,000/100,000 — not area-proportional, since the range is 3,500×), one colour per
pathogen. Colour is redundant with the x axis on purpose: only 9 substantive hues exist, so
`distinct_colors` reuses them as tints beyond that, and the labelled axis is the real key. Shows what
the pie version cannot: the largest datasets sit at near-zero active fraction while small ones are
balanced or active-heavy. A short bar per column marks the pathogen's **unweighted mean** active
fraction (size-weighted would be much lower, since the huge datasets are near 0% active). Size key is
nested circles with leader-lined labels. Both panels are rendered — pick one at layout time.

**Area encoding (both pathogen panels):** a mark's *area* encodes its model's training-set size
as `n ** AREA_EXPONENT` with **`AREA_EXPONENT = 0.5`**, i.e. area ∝ √n, set in
`src/plots_metadata.py`. Sizes span 96 to 5,000,000 compounds, so a linear encoding buries the
small datasets. The exponent is the one knob: 1.0 is true data volume, →0 is pure model count.

**Why not `1 + ln(n)`** (the earlier choice): every model then contributed at least ~5.5 to its
pathogen's total regardless of size, so a pathogen's total tracked how *many* models it held
rather than how much data — *E. coli* (12 models, 2.4M compounds) came out **larger** than
*P. falciparum* (7 models, 6.0M). A power transform makes a small model contribute small area, so
count stops inflating groups. At 0.5 the ordering follows data volume (*Pf* 33.6% vs *Ec* 21.5%
of the Voronoi panel) and the smallest cell stays ~9x above the readable floor. Note the count
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

**Pathogen Voronoi treemap (`pathogen_voronoi`)** — a space-filling alternative to the circle
panel, rendered alongside it so the two can be compared; neither has been chosen yet. Two
area-accurate levels: the panel rectangle is tessellated into one region per pathogen (area =
that pathogen's share of all training compounds), each region into one cell per model. Geometry
is in `src/voronoi_treemap.py` (additively weighted power diagram, weights fitted by gradient
ascent on the optimal-transport dual — no new dependencies). The script prints the achieved
max relative area error, which is **0.9%**, plus any cells too small to read.

**Area metric:** defaults to `power` (√n, shared with the circle panel). Passing
`area_metric="linear"` gives areas strictly proportional to compound count — faithful, but it
does **not** produce a readable figure with this data: 9 of 15 regions are then too small to hold
their own label and 11 of 51 cells fall below 1e-4 of the panel, because two models (eos4zfy 5M,
eos5nqn 2M) carry ~70% of all training data. Under √n, 0 cells are unreadable and 5 regions are
still too small to label. The figure prints its area error and both legibility counts each run.

**Once a panel is chosen, delete the other** rather than leaving both in `output/`.

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
a small blue dot per retained model's mean CV AUROC, a grey min–max guide line, and a translucent
red dot at the **mean** (the representative per-pathogen value), sorted best-first. The mean dot's
**area encodes the number of unique molecules** in the pathogen's ChEMBL dataset (`n_cleaned_inchikeys`
from the staged step-02 `chembl_curation/general/27_chembl_coverage.csv`), area-scaled between a min
and max marker so small datasets stay visible; a framed size legend gives round references
(150 / 10,000 / 100,000 molecules). Dashed line = the step-10 retention floor (mean AUROC 0.7).
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
(1) own-assay AUROC for the 7 organisms with an EU OpenScreen primary assay; (1b) the same models on
the merged secondary (confirmatory) assays, for a primary-vs-secondary comparison; (3) exclusive vs
shared-hit AUROC; (4) a model × EU OpenScreen assay AUROC matrix (off-diagonal = a model predicting
a *different* organism's data) plus a per-model specificity index; (4c) the label-only hit-promiscuity
distribution (how many actives are hits in 1, 2, … 7 pathogens); (4d) the summed consensus score
across the 7 models by hit class (inactive / exclusive / narrow / broad), the maximum score by plain
activity (raw and percentile-normalised), and the rank each exclusive hit's own model gives it among
the 7, each of the last two also with training-set compounds removed. Writes twenty-two summary CSVs
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
**Secondary (analysis 1b):** own-assay AUROC on `06_subset_data/secondary/{code}_secondary.csv` —
the upstream merge of every **non-primary** assay for the organism (academic sub-screens +
dose-response/IC50), deduplicated active-prevails, from the `eu-openscreen-antimicrobial-tasks`
repo. Written to `05_euopenscreen_secondary_auroc.csv` and drawn as a primary-vs-secondary grouped
bar panel. Note the secondary set pools assays of mixed type/size/cutoff, so it is not a like-for-like
single screen — treat the comparison as indicative.
**Exclusivity (analysis 3):** uses the precomputed EU OpenScreen subsets in
`06_subset_data/exclusivity/` — "exclusive" = active in only 1 of the 7 primary assays (i.e.
organism-specific), "shared" (non-exclusive) = active in ≥2 (i.e. pan-active) — each scored against
the same primary inactives, with `n=` labels. This *is* the pan-active-vs-specific specificity test
(a drop toward 0.5 on the exclusive/specific bars = the model captures generic, not organism-specific,
activity); an on-the-fly recompute was checked to match these upstream subsets to ≤0.0002 AUROC, so
only this one is kept.
**Active-set overlap (analysis 4b, `05_active_overlap.csv` + `active_overlap_jaccard`):** label-only
pairwise Jaccard between the 7 assays' active sets. The hit-sets overlap substantially (most actives
are pan-active), which is *why* the cross-organism AUROCs (analysis 4) are high — high off-diagonal
AUROC reflects shared promiscuous hits, not necessarily cross-organism prediction. Read the two together.
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
that KEEP the compounds the models were trained on (13 CSVs, 5 figures); `deduplicated/` holds those
with every training-set compound REMOVED (10 CSVs, 9 figures); the top level keeps only what has no
leakage dimension — the label-only tables (`05_active_overlap`, `05_hit_promiscuity`,
`05_promiscuous_hits`) and the leakage audit itself (`05_leakage_report`) — 4 CSVs, 3 figures. Each
dir has its own `png/`, `pdf/` and `figure_cells.json`. The five long-form metric tables
(`05_euopenscreen_auroc`, `_secondary_auroc`, `_roc`, `05_hit_exclusivity`, `05_cross_organism_euos`)
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
  *S. aureus*' 379 actives and 18 of *C. albicans*' 172, and 3 of *P. aeruginosa*' 14. So dedup
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
- **ROC grids** (`eos3dys_roc.csv` + `eos3dys_roc_grid_{inhib_50,mic_25}{,_exclusive}`): one grid per
  endpoint metric — 6 organisms, dedup — over all actives and over exclusive actives only. Dedup
  AUROC all → exclusive: *C. albicans* 0.90→0.88, *S. aureus* 0.89→0.86 and *A. baumannii*
  0.88→0.82 hold up on `inhib_50`, while *E. coli* 0.78→0.48 and *K. pneumoniae* 0.74→0.56 collapse.
  `mic_25` is weaker throughout. *P. aeruginosa*'s exclusive panel rests on a single active.
- **Max endpoint probability, active vs inactive** (`eos3dys_consensus_max_{boxstats,actives}.csv` +
  `eos3dys_consensus_max_by_activity`): per compound the highest of the 12 matched endpoint
  probabilities. **No percentile normalisation**, unlike the ChEMBL twin: the 7 ChEMBL models' maxima
  ran 0.753–0.974 so low-scoring ones could never win a max, whereas all 12 eos3dys endpoints top out
  at ~1.0. Dedup medians 0.734 (inactive, n=100,007) vs 0.865 (active, n=535). Note the max of 12
  correlated probabilities is high by construction, so read the gap, not the absolute level.

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

## 07_prediction_correlations.py
An **exploratory** assessment of how much the annotation models AGREE. All ~130 models were run by
`00_download_data.py` on the same ~1.35M-compound reference library
(`data/processed/annotation_preds_ref_library/`), aligned on the `key` hash — a clean rectangular
matrix. Works at the **output-column level** (each of the ~1500 output columns is a node), computes a
library-wide Spearman correlation matrix and top-N hit-overlap, and highlights two questions: do
different cytotoxicity models correlate, and do same-organism models correlate? Analysis in
`src/eval_correlations.py`, figures in `src/plots_correlations.py`; writes to
`output/07_prediction_correlations/`. Runs in two stages — `python 07_prediction_correlations.py`
(build the matrix + auto groups) then, **after reviewing the auto groups**,
`python 07_prediction_correlations.py --analyze` (correlations + figures); `--all` skips the gate.
**Sample:** correlation is computed on a fixed-seed (`RANDOM_SEED=42`) **200,000-compound** sample
(`CORR_SAMPLE_N`) so the ~1500-node matrix is tractable; Spearman ρ is essentially stable at this
size. Raise `CORR_SAMPLE_N` (or set it to `None` for the full library) for the final, narrowed
analysis. The sampled `key × node` matrix is cached as `07_score_matrix.parquet` and reused on re-run.
**Correlation:** **Spearman** (rank-based) so heterogeneous output scales (probabilities, regressions,
categoricals) are comparable; imputes any NaN prediction cells to the column mean rank (count logged).
**Top-N overlap:** Jaccard of the top-**100** and top-**500** highest-scoring compounds
(`TOPN_CUTOFFS`), computed **only for probability-type output columns** (values in `[0,1]`, a known
"higher = more" direction) — regression/categorical columns have no meaningful "top" so they enter the
Spearman matrix but not the overlap. **Value-type tagging:** a column is `categorical` if integer-like
with ≤ `CATEGORICAL_MAX_UNIQUE` (=10) distinct values, `probability` if all values ∈ `[0,1]`, else
`continuous` — a display heuristic, revisable. **Grouping needs sign-off:** organism groups come from
splitting Airtable `Target Organism`; cytotoxicity from a regex (`CYTOTOX_REGEX`) over
Tag/Title/Interpretation/Description. The auto assignment is written to `07_group_assignments.csv`
**for the user to review and edit before the `--analyze` stage draws focus-group figures.**
**Deferred:** narrowed follow-up analyses guided by this first pass (e.g. full-library correlation,
per-organism deep dives, redundancy pruning).
