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
`data/processed/annotation_preds_ref_library/eos4djh_v1.csv`, feeding step 08. Both pin the major
version only (`v1`), matching the loop. **`eos4djh` lands in the annotation folder despite not being an
annotation model** — the property steps read every property model from that one directory, and the model's real
Task is recorded on `PHYSCHEM_MODEL_ID` in `src/default.py` rather than in the path.
**Key decision:** ChEMBL model performance metrics are recomputed downstream from the per-fold
`09_reports`; the staged `10_reports/` CSVs are kept only for fields that cannot be reconstructed
(quality weights, decision cutoffs, discarded-model reasons).

## 01_ersilia_metadata.py
Analyses the Airtable model metadata (Ready models only, n = 218 on the hand-revised
2026-08-14 export) and renders 13 panels — the
task/subtask breakdown, the source-type and output composition, biomedical area, target organism,
three composition donuts, the container-metrics row and the pathogen circle-treemap.

> **Note (2026-08-04):** this section was rewritten from the code after an accidental deletion of the
> previous text. It describes what the code does today; any decision recorded only in the old prose
> is not here. Review it before trusting it as the decision log.

**Shown vs counted:** value counts for every counted field are written to `*_counts.csv`, but not all
are plotted. **Counted only, no panel:** Publication Type, Tag and the ungrouped License field (the
licence composition reaches the figure as the reuse-class donut, not per identifier). Target Organism
is capped at the **top 10** categories in its panel; Biomedical Area is drawn as **five groups over
Activity prediction models only**, not as its 18 raw areas. Full counts remain in the CSVs either way.

**Stacked panels:** Source Type and Output are each drawn as one stacked bar panel segmented by
Subtask, so a single panel carries the joint distribution instead of two. Segments use the subtask
palette (shades of the parent task's hue), never a palette of their own — the field is already encoded
by the bar position. Cross-tabs are written as `<field>_by_subtask_counts.csv`.

**Task hues:** Annotation = crimson, Representation = amber, Sampling = lime (`TASK_HUES`, the one
place to change them). `SUBTASK_COLORS` derives from the same mapping, so a task and its subtasks can
never drift onto different hues. `SOURCE_TYPE_COLORS` now survives only for the pathogen panel's dots
— the two sets never carry colour in the same panel.

**Task/subtask alternatives** — two ways to draw the same 218 models, on the same colours:
- `task_subtask` — bars, one per subtask, grouped by task and ordered by count within it. Best for
  reading counts. No legend: each bar is named on the axis next to its own colour, which makes this
  panel the subtask colour key.
- `task_subtask_waffle` — one square per model, **16 wide × 14 rows** (the script prints `blank`,
  currently **6**). A ragged last row is unavoidable at 218 = 2 × 109, which has no divisor in the
  usable 13–18 column range; 16 divided the earlier 208 exactly and is kept for the aspect ratio. Shows n
  rather than stating it; reading 52 vs 39 means counting. Carries its own legend with counts, so it
  is self-contained and can serve as the key too. That legend band is also what squares the panel —
  the grid alone crops landscape (1.23 at 60 mm), and at quarter width the legend labels have to be
  abbreviated (`_LEGEND_ABBREV`) or the key sets the crop and the squares fall to 1.99 mm.

**Legends:** the two `*_by_subtask` stacks carry none — no room at ~52 × 25 mm each, where a 6-entry
key would be several times taller than the plot, and they are meant to sit beside a key panel. **So if either goes into the figure, `task_subtask` or
`task_subtask_waffle` must travel with it.** The waffle and each of the three donuts carry their own.

**Biomedical Area is five groups, not 18 areas** (`BIOAREA_GROUP` in `src/default.py`, signed off
2026-08-02, **Antifungal split out of Antimicrobial 2026-08-07**), over **Activity prediction models
only**: Antimicrobial 51, ADMET 27, Other 10, Antiviral 7, Antifungal 3. Two things the counting has to get right, both printed on every run:
- Biomedical Area is **multi-value**, so these are counts of *distinct models* per group, not of area
  assignments. Grouping absorbs most of the multiplicity.
- Five models (`eos2zmb`, `eos3f8h`, `eos60mw`, `eos7kpb`, `eos8jx6`) carry areas in two different
  groups and are **counted in both**, so the bars sum to 98 over 93 models. That is the metadata's own
  claim, left unresolved rather than silently reassigned. The donut shows the 93 in its hole, so a
  caption must explain the discrepancy.
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
`Task == "Annotation" & Status == "Ready"` — a test it previously failed. It is a **human GPCR** model
with no antimicrobial relevance and appears in no config file, so its 10 endpoints do not enter
`config/08_endpoint_selection.csv` without a deliberate decision.

**Metadata snapshot refreshed 2026-08-06: 219 → 225 models, Ready 208 → 214.** Superseded by the
08-07 freeze below; kept as the record of that move. Six models were added — `eos19dk` (molcompass,
Projection), `eos5g6m` (glacier-embeddings), `eos5mnx` (sand-shape-descriptor), `eos6pj2`
(nafm-embeddings), `eos8zvb` (pymolgen) and `eos84nf` (genmol-scaffold-decoration, *In progress*) —
**none of them an Annotation model**, so none adds bioactivity endpoints. Four models changed Status:
`eos8vud`, `eos69e6` and `eos4qda` became Ready, while `eos18ie` (antibiotics-ai-saureus) and
`eos1lb5` went to *In maintenance*. Ready-only Task counts: Annotation 131 → 130, Representation
58 → 61, Sampling 19 → 23.

### Metadata HAND-REVISED 2026-08-14 — 218 models, 218 Ready

**The catalogue was revised by hand offline and the revisions are NOT yet in Airtable.** The corrected
export is `data/raw/airtable_metadata_manual.csv`, and `src/default.py` now pins
`AIRTABLE_METADATA_FILE` to that literal filename instead of deriving it from
`AIRTABLE_SNAPSHOT_DATE`. **Re-downloading would undo the revisions**, which is why the derived name
is commented out rather than the date bumped. `AIRTABLE_SNAPSHOT_DATE` stays at `2026-08-12` as the
record of the last genuine Airtable pull. **To revert once the sync lands:** bump the date, restore
the f-string, re-run `00_download_data.py`. Both `00` and `01` now print the *filename* alongside the
date, so no run can be mistaken for a dated snapshot.

Schema is identical (same 44 columns). The four `In progress` models (`eos1tt2`, `eos48ue`, `eos6ru5`,
`eos6wdw`) are gone, so the total falls 222 → 218 while **Ready is unchanged at 218** — no panel loses
a model. Every model's Description and Interpretation was rewritten; `Image Size` and the five
`Computational Performance` columns differ on most rows but only in the second decimal (rounding, not
re-measurement).

**Substantive field changes:**
- **Task/Subtask, 4 models.** `eos1vms` (chembl-multitask-descriptor) Representation/Featurization →
  **Annotation/Activity prediction**; `eos526j` (aizynthfinder) → **Annotation/Property calculation**;
  `eos4k4f` (standardization) Annotation → Representation/Featurization; `eos935d` (meta-trans)
  Annotation/Activity prediction → **Sampling/Generation**. Ready Task counts: Annotation 133 → 133,
  Representation 61 → **60**, Sampling 24 → **25**; Featurization 54 → 53, Generation 13 → 14.
- **Biomedical Area, 9 models.** Six pathogen models gained `Antimicrobial resistance` (28 → **35**,
  the largest single move in that figure); `eos1amr` gained `Alzheimer`; `eos60mw` (cidalsdb) gained
  **`Leishmaniasis`, a new vocabulary value** that had no `BIOAREA_GROUP` entry and therefore *raised*.
  **Decision (2026-08-14): `Leishmaniasis` → `Antimicrobial`**, on the same reasoning already recorded
  for Malaria — *Leishmania major* is a protozoan, inside "antimicrobial" only on the broad clinical
  definition. Grouped counts move Antimicrobial 49 → **51**, ADMET 28 → 27, Other 9 → **10**, and the
  models spanning two groups go 3 → **5** (adding `eos60mw` and `eos8jx6`).
- **License, 22 models**, including **three vocabulary values with no `LICENSE_CLASS` entry** —
  `Non-commercial` (×2), `CC-BY-NC-SA-4.0`, `NCSA` — which the donut was **silently dropping**, since
  it groups by class and drops NaN. **Decision (2026-08-14):** `Non-commercial` and `CC-BY-NC-SA-4.0`
  → **Non-commercial**; `NCSA` (the University of Illinois/NCSA licence, MIT/BSD-style with no
  share-alike duty) → **Permissive**. A **raise-on-unmapped guard** was added to `01` so this cannot
  recur silently. Classes move Permissive 111 → **116**, Copyleft 79 → **67**, Non-commercial 1 → **4**,
  Not recorded 27 → **31**. `LGPL-3.0-only` leaves the vocabulary.
- **Target Organism, 5 models.** `eos3ev6`, `eos5jz9`, `eos7nno`, `eos935d`: `Any` → `Homo sapiens`
  (`Any` 125 → 121, `Homo sapiens` 37 → **41**). `eos74km` reordered only.
- Publication Type 9 models, Publication Year 28, Output Consistency 4 (`Variable` → `Fixed`).

**⚠️ OPEN: the "named organisms are Annotation-only" colour rule now has one exception.** `eos935d`
(meta-trans) moved to **Sampling/Generation** while keeping `Target Organism = Homo sapiens` and
`Biomedical Area = ADMET`. The `target_organism` panel colours named bars with the **Annotation hue**
on the strength of that rule, so its `Homo sapiens` bar (41 models) now contains 1 Sampling model. The
Biomedical Area panels are unaffected — they are restricted to Activity prediction, which excludes it.
Nothing raises; the rule is prose, not an assertion. **Needs a decision:** accept it as a 1-in-41
approximation and say so in the caption, correct the model's Airtable record, or drop the rule and
colour those panels differently.

**Endpoint selection unaffected — verified.** No model in `config/08_endpoint_selection.csv` changed
`Task`, and the only `Target Organism` change inside it (`eos74km`) is a reordering of the same set.
All four reclassified models have **zero rows** in the config, as do all four `Any → Homo sapiens`
models. The two models newly visible to `00`'s Section 4 loop get **no endpoint rows**: `eos526j` is
retrosynthesis planning, not a bioactivity endpoint, and `eos1vms` declares `Target Organism = Any` —
so neither enters without a deliberate decision, the same rule already applied to `eos93h2` and
`eos3f8h`. **Steps 07–14 need no re-run.**

### Metadata snapshot RE-FROZEN 2026-08-12 — 222 models, 218 Ready

**Taken to repair a data-destroying bug in this script, not to refresh the catalogue.** The download
round-tripped the Airtable export through pandas without `keep_default_na=False`. `None` is a
legitimate value of the License vocabulary *and* a member of pandas' default NA set, so every model
that correctly declares "this upstream repo has no LICENSE file" arrived in the snapshot
indistinguishable from a model that declares nothing. 28 models were affected, and the loss
manufactured 28 spurious repo-vs-Airtable drift findings in `tools/audit_model_metadata.py`. Both
reads now pass `keep_default_na=False`; the 2026-08-07 and earlier files still carry the defect and
**must not be used to judge whether a model declares a licence**.

Two models left the shared view between the two snapshots (`eos2gth`, `eos7d58`), so Ready goes
220 → 218. Anything downstream that quotes a model count needs re-running.

### Metadata snapshot FROZEN 2026-08-07 — 224 models, 220 Ready, 0 in maintenance

**The Airtable metadata is no longer refreshed implicitly.** The snapshot date is a constant in
`src/default.py`:

```python
AIRTABLE_SNAPSHOT_DATE = "2026-08-07"
AIRTABLE_METADATA_FILE = f"airtable_metadata_{AIRTABLE_SNAPSHOT_DATE.replace('-', '')}.csv"
```

Every consumer — `00_download_data.py`, `01_ersilia_metadata.py` and
`tools/build_endpoint_selection_template.py` — reads the **dated** file, at that point
`data/raw/airtable_metadata_20260807.csv`. There is deliberately **no generic
`airtable_metadata.csv`** any more, so a script that misses the constant fails loudly instead of
silently reading a stale file. The previous snapshot stays on disk as
`data/raw/airtable_metadata_20260806.csv`. **To take a new snapshot, bump the date in
`src/default.py` and re-run `00_download_data.py`** — it downloads the new dated file and leaves the
old ones untouched. Both `00` and `01` print the snapshot date on every run.

> **Superseded 2026-08-14** (see the section at the top of this log): `AIRTABLE_METADATA_FILE` is now
> pinned to the literal `airtable_metadata_manual.csv` rather than derived from the date, and the two
> scripts print the filename alongside the date. The date-derived mechanism above is what to restore
> once the hand revisions reach Airtable.

**Every "n = 208" and "n = 214" above this line predates the freeze.** Ready-only Task counts:
Annotation 130 → **135**, Representation 61 → **61**, Sampling 23 → **24**. By subtask: Activity
prediction 91 → 95, Property calculation 39 → 40, Generation 12 → 13; Featurization (54), Projection
(7) and Similarity search (11) unchanged.

**Five models were repaired and returned to Ready** — this snapshot was taken specifically to capture
that: `eos18ie` (antibiotics-ai-saureus), `eos1lb5` (mycobacterium-permeability), `eos9ivc`
(anti-mtb-seattle) and `eos9tyg` (ncats-pampa74) all moved *In maintenance* → **Ready**, and `eos84nf`
(genmol-scaffold-decoration) moved *In progress* → **Ready**. Four of the five were also repackaged
(`Last Packaging Date` = 2026-08-07) and three gained an ARM64 build, so their container metrics moved
too — `eos1lb5` and `eos9ivc` roughly tripled in image size (~2.0 → ~5.8 GB) on the rebuild, and
`eos18ie` gained its first `Computational Performance 3` measurement (was the `-1` sentinel, now
1492.7). **The technical box row therefore shifts for reasons unrelated to model count.**

**`eos18ie`'s maintenance inconsistency is RESOLVED.** It was row 1 of
`config/08_endpoint_selection.csv` (`saureus_inhibition_probability`, `Yes`) and in the step-07 cache
while dropping out of every `Status == "Ready"` panel. Now that it is Ready again the score matrices
and the metadata panels agree, and no decision is outstanding.

**Two models left the dataset entirely — `eos3nn9` (mpro-covid19) and `eos7asg` (padel).** They were
*In maintenance* on 08-06 and are now **`Archived`**. This is the one snapshot-diff trap worth
knowing: **the shared view filters `Archived` out**, so an archived model does not appear with a
changed Status — it vanishes from the CSV export, and a naive row-count diff reads it as a deletion.
Verified against the base directly (both records still exist, Status `Archived`). Neither is
referenced by any config file or any script, so **they are dropped from all figures**, per decision on
2026-08-07. This is why the total falls 225 → 224 while Ready rises 214 → 220.

**One model was added: `eos3f8h` (eu-openscreen-hts)**, Ready, Annotation / Activity prediction — the
only new record. It is the first *added* Annotation model in either recent snapshot, so unlike the
08-06 refresh **this one does add candidate bioactivity endpoints**. Its endpoints are **not** in
`config/08_endpoint_selection.csv`; adding them is a scientific selection decision and has not been
made. Note its direct relevance to step 05 (EU OpenScreen validation).

**New Biomedical Area → Antifungal is now its own group (decided 2026-08-07).** `eos3f8h` carries
`Fungal infections`, an area no previous snapshot had, and `01`'s unmapped-area guard raised on it
exactly as designed. Rather than folding it into Antimicrobial, the fungal areas were **split out**:
`Candidiasis`, `Mycetoma` and `Fungal infections` now map to a fifth group, **Antifungal** (3 models —
`eos3f8h`, `eos8jx6`, `eos4f95`). Consequences to know:
- **A caption saying "antimicrobial" no longer covers antifungal activity** — it is reported
  separately. The protozoal/helminth stretches (Malaria, Schistosomiasis) are unchanged.
- The donut's hatch scheme gained a fifth pattern. `Antifungal` takes the **backslash, mirroring
  ADMET's forward slash**, so the two read as a pair; this is a **deliberate exception to the
  ink-ordering rule** (its ink matches ADMET's despite being the smallest group), documented in
  `src/plotting_colors.py`. Every pattern lighter than the dots aliased on a wedge that thin.
- **At n = 3 the Antifungal wedge is a sliver whose hatch is not legible in the ring** — it reads only
  in the legend swatch. The bar strip shows it cleanly; prefer that panel if the group matters.
- The strip's bar pitch drops 4.55 → 3.64 mm and the bar 3.2 → 2.55 mm. Still seats the 5 pt label,
  but `_STRIP_BAR_FRACTION` is the knob to revisit if a sixth group ever appears.

**Isaura re-run (Section 4) on the new snapshot: 3 fetched, 3 still missing.** Downloaded on 08-07 —
`eos9tyg` (ncats-pampa74, 116 MB) and `eos93h2` (image-mol-gpcr, 352 MB, first fetch since the 08-06
Task correction made it visible to the loop). **`eos3f8h` (eu-openscreen-hts) was not yet precalculated
on 08-07 and was fetched on 2026-08-11** once it landed in Isaura: 288 MB, 1,355,109 rows covering the
full reference library, **7 continuous per-pathogen columns** (`abaumannii`, `calbicans`, `ecoli`,
`efaecalis`, `kpneumoniae`, `paeruginosa`, `saureus`), all in (0, 1) with **no nulls**. Its endpoints
are therefore now available to the score matrices, though they are still **not** in
`config/08_endpoint_selection.csv` — that selection decision is open.

Still absent from Isaura: `eos3wzy`, `eos935d`, `eos2b6f` — **pre-existing gaps**, already reported
missing by the 07-30 run and unrelated to this snapshot. All other prediction files cover the full
1,355,109-compound reference library.

### eos3f8h endpoints added to the selection — 2026-08-11 (selection 300 → 307)

**All 7 of `eos3f8h`'s columns were added to `config/08_endpoint_selection.csv` as `selected = Yes`**
(user-directed), taking the file 395 → 402 rows and the selection **300 → 307 endpoints**. Field values
are taken from the model's own Airtable record, not inferred: `direction = higher` (Interpretation:
*"higher values indicate greater predicted probability of growth inhibition"*), `assay_type =
bioactivity` (single-point growth inhibition, 50–70 % cut-offs at 41.7–50 µM), `sensitivity =
wild-type` (no resistant or sensitized strains in this assay set).

**Steps 07, 08, 09, 13 and 14 were rebuilt at 307 endpoints on 2026-08-11** (table at the end of this
file). Steps 10, 11 and 12 were deliberately **not** re-run: none reads the endpoint selection, and
step 11's parquet is its own separate cache, so their outputs are unaffected and step 14 reads them
as-is. Note that step 07 **skips on file existence, not content**, so it does not notice a changed
selection on its own — the previous outputs had to be moved aside to force the rebuild. That
300-endpoint backup was **deleted with `tmp/` on 2026-08-12** (user-directed), so the only record of
those numbers is the tables in this file; the build itself is reproducible only by reverting the
selection and re-running.

**Six of the seven organisms already have endpoints in the selection**, so these are largely
*additional* endpoints for existing pathogens rather than new coverage: *S. aureus* 22 → 23,
*E. coli* 23 → 24, *C. albicans* 20 → 21, *P. aeruginosa* 19 → 20, *A. baumannii* 13 → 14,
*K. pneumoniae* 13 → 14. Expect them to be **correlated** with the endpoints already there, which
matters for any statistic that treats endpoints as independent.

**`Enterococcus faecalis` is the exception and needs a decision.** It had **zero** endpoints before, so
`eos3f8h` is its only one — and it is **not in `config/pathogens_of_interest.csv`** (that file lists
*E. faecium*, a different species). Consequences: it does not become a row in the pathogen matrices
without being added there, and as a single-endpoint organism its aggregate score would simply *be* that
endpoint's percentile rank, not a merge — the same caveat step 14 already records for Campylobacter,
Enterobacter, *E. faecium*, *H. pylori* and *S. pneumoniae*.

**Circularity warning for any EU OpenScreen evaluation.** `eos3f8h` was trained by Ersilia on the EU
OpenScreen HTS data (source: `ersilia-os/eu-openscreen-antimicrobial-tasks`), which is the **same
companion repo that supplies step 05's validation tasks**. Step 05 does not read this config, so
nothing is currently circular. But **evaluating `eos3f8h` against EU OpenScreen ground truth would be
scoring it on its own training data** — its reported mean AUROC of 0.94 is 5-fold cross-validation on
that set, not external validation. Keep it out of any EU-OpenScreen-based benchmark.

**Technical box row.** Three **horizontal** panels per Task on a log x axis, sharing one task axis and
occupying **4/6 of the page width (120 mm) × 30 mm** between them: `runtime_100`, `image_size` and
`output_dimension`. Only the leftmost draws the task tick labels, so the other two **cannot be placed
on their own**; `runtime_100` is wider by exactly its 14.2 mm label column so all three metric axes
come out the same ~28 mm length. Per-task quartiles go to `technical_metrics_summary.csv`.

**`n` is not in the tick labels** — coverage differs per metric (runtime 131/59/13, the other two
133/60/25), so a shared label set cannot carry it. **A caption must state that the runtime box for
Sampling rests on 13 of 25 models.**

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

- `license_class_donut` — Permissive 116 turquoise, Copyleft 67 periwinkle, Not recorded 31 silver,
  Non-commercial 4 **fuchsia** (the one use of that hue in the repo, because a 4-of-218 wedge in any
  calmer colour is invisible). Four **reuse classes, not the twelve licences**: five licences have
  exactly one model each, a 1.7° wedge that cannot be seen or labelled, so per-licence detail stays in
  `license_grouped_counts.csv` (and the ungrouped `license_counts.csv` beside it, which keeps the
  `-or-later` / `-only` distinction the grouped file collapses). A licence with no `LICENSE_CLASS`
  entry **raises** (added 2026-08-14) rather than being dropped from the donut — the same guard the
  Biomedical Area grouping has.
- `docker_architecture` — **AMD only 77 tangerine, AMD + ARM 141 cobalt**, ordered base-capability
  first. Moved off turquoise/periwinkle so no hue means two things across the set; tangerine is a plain
  categorical hue here, **not** a warning about x86-only builds.
- `biomedical_area_donut` — the alternative to the bar strip, on the same five groups. **One hue, five
  patterns**: every model in it is an Annotation model, so colour would encode nothing; groups separate
  by hatch instead (solid → diagonal → dots → cross-hatch for the catch-all, with Antifungal on the
  mirrored backslash). **The hole reads 95 while the legend rows add to 98** — three models carry areas
  in two groups. That is real; a caption should say so. **At n = 3 the Antifungal wedge is a sliver and
  its hatch is not legible in the ring itself** — only in the legend swatch. If that group matters to
  the argument, prefer the bar strip, where it reads cleanly.

**Runtime batch size: 100 molecules** (`RUNTIME_BATCH` in `src/default.py`; the five
`Computational Performance` columns are 1/10/100/1,000/10,000 molecules). A `-1` in those columns means
the benchmark was **never run**, not zero — those models are skipped, never imputed. 100 is chosen
because it is the largest batch where generative models still have data: coverage is **131/133, 59/60,
13/25** at 100 molecules but **125/133, 58/60, 1/25** at 1,000. The cost of that choice: at 100
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

**Denominator: the same 218 as every panel, since 2026-08-14.** The series runs on the **unfiltered**
metadata while every panel here uses `Status == "Ready"`, and the two used to differ — a model in
maintenance was still incorporated on its date and still counts towards how the hub grew. The
hand-revised export contains **only Ready models**, so unfiltered and Ready are both 218 and the
caveat is currently moot. **It returns the moment a non-Ready model reappears**, so keep reading the
denominator off the run log rather than assuming the two agree.

**Excluded: models with no `Incorporation Date`.** The date is the x axis, so there is nowhere to
put them. They are excluded, never imputed, and the script prints the count each run — **currently
0**, since the four undated models were the `In progress` rows the 2026-08-14 revision removed.

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

**Headline numbers** (2026-08-02 GitHub snapshot, as plotted): models reach **211**, commit authors
**107**, **15,828 commits** and **1,469 issues** in total. The Models track stops at 211 rather than
218 because the partial month 2026-08 is excluded from the plotted series (it stays in the CSV), and
the seven most recent models were incorporated in it.

**Denominator: the same 218 as step 01, since 2026-08-14.** The Models series runs on the
**unfiltered** metadata while every panel in step 01 uses `Status == "Ready"`. The hand-revised
export contains only Ready models, so the two coincide at 218; the distinction returns as soon as a
non-Ready model does.

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
Builds the foundation of the correlation analysis: the full-library prediction matrix cache. All
annotation models were run by `00_download_data.py` on the same ~1.35M-compound reference library
(`data/processed/annotation_preds_ref_library/`), aligned on the `key` hash — a clean rectangular
matrix. **Defines the named matrix**, one row per compound and one column per `selected == Yes`
endpoint (**307** as the config stands, though the built cache is still 300 — see the eos3f8h note
below), and caches it as `07_score_matrix_full.parquet`. Every downstream step re-derives whatever
scaled/row-normalized variant it needs from that cache in memory (column z-score or rank-percentile,
then row L2 or L1 normalization — `src/eval_correlations.py`'s `scale_matrix`/`row_normalize`); none
is written to disk as a CSV. Writes to `output/07_score_matrices/`.

```
python 07_score_matrices.py   # parquet cache + mean-rank figure  (~5 min)
```

**A `--write-matrix-csvs` flag used to export five scaled/normalized matrices as CSVs (~23.5 GB,
~1 h 54 m) — removed 2026-09-01** (see "Pipeline reorganization" below), after confirming it had
gone unused since being added: only 3 of the 5 variants were ever consumed by anything, even
in-memory, and none was ever read back from disk by any script. Kept for the record, since the
reasoning still applies to any future "materialise a variant" request: `named`, `zscore_l2rownorm`
and `rankpct_l1rownorm` fed step 09's Jaccard matrices; plain `zscore` and plain `rankpct` (no row
normalization) never had a consumer at all, because top-N Jaccard depends only on each column's own
internal ranking, and both column scalings are strictly increasing *per column* — neither can change
a column's top-1000 set, so their Jaccard matrices are provably identical to the unscaled baseline's
(step 09 asserts this at runtime rather than assuming it). Deriving the 3 real variants from the
parquet in memory measured at **55 s** total (2026-08-07, this hardware) against **~84 s** to read
the equivalent CSVs back — cheaper every run, with no risk of the materialised copy drifting out of
sync with the config (which happened once: the CSVs sat at 260 columns while the selection said 300).

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
- The raw prediction CSVs (~15 GB) are read at most **once**, into `07_score_matrix_full.parquet`.
- **The `260`-based numbers throughout the step 09–14 sections below were measured under the earlier
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

**Why it lives here now:** it needs the named matrix and the `rank_pct` scaling, which this script
already builds from the parquet cache, and it consumed nothing else. As a separate step it rebuilt
the same matrix from the same parquet a second time. The mean-rank block is **skip-if-exists** on
`07_mean_rank_per_compound.csv`. The former `scripts/09_mean_rank_distribution.py` was removed from
the tree (recoverable via `git show <rev>:scripts/09_mean_rank_distribution.py`), and its stale
260-endpoint `output/09_mean_rank_distribution/` was deleted.

**Key decisions, for review:**
- **Mean, not sum.** The request was phrased as a column-wise sum; the reported quantity is the mean,
  which is what "average rank" means and keeps the value on the interpretable [0, 1] percentile scale.
  They differ only by a constant factor — multiply by `n_endpoints` for the sum.
- **The centre of this distribution is fixed by construction, not a finding.** Each percentile-rank
  column has mean ~= 0.5, so the grand mean of the row means is ~= 0.5 necessarily (observed **0.5000**
  at 307 endpoints). Only the **spread and shape** carry information: observed **SD 0.1305** against the
  ~0.0165 that 307 mutually independent endpoints would give — a ~7.9x inflation — over a range of
  **0.1495–0.8645**, i.e. no compound ranks in the top decile of everything nor the bottom decile of
  everything. (At 300 endpoints: SD 0.1289 against ~0.0167, a ~7.7x inflation, range 0.1518–0.8616 —
  adding seven correlated endpoints widened the spread slightly rather than averaging it away.) What that inflation means substantively (shared chemistry, correlated training data,
  genuine broad-spectrum compounds) is **not** something these numbers settle.
- **Quantiles** at 307 endpoints: 0.1945 / 0.2351 / 0.2860 / 0.3990 / 0.5026 / 0.6012 / 0.7080 / 0.7605
  / 0.8047 for the 0.1st, 1st, 5th, 25th, 50th, 75th, 95th, 99th and 99.9th percentiles.
  (At 300: 0.1973 / 0.2379 / 0.2887 / 0.4004 / 0.5023 / 0.5997 / 0.7059 / 0.7580 / 0.8017.)
- **No rows dropped for missing values** — the 15 partially-scored compounds
  (`pfalciparum__eos4zfy__maip_score`) are averaged over the 306 endpoints they have, with
  `n_endpoints` per row in the CSV so they stay visible. NaN-skipping mean, never zero-filled.
- Follows the standard publication figure convention (`BasePlot`, PNG + PDF + `figure_cells.json`).
  **This is the only figure step 07 produces**, so `figure_cells.json` in that folder holds exactly one
  entry.
- Previous 260-endpoint values, for comparison: mean 0.5000, SD 0.133, range 0.140–0.869.

### Pipeline reorganization — 2026-09-01

Steps 07-15 were renumbered and one script split off, to fix a real circular dependency and match
the numbering to the order the scripts actually need to run in (previously: 08 pathogen Jaccard, 09
projection, 10 physchem matrix, 11 abx matrix, 12 cytotox matrix, 13 curated predictions, 14 AUROC
matrix, 15 the per-endpoint audit — now: 08 property matrices *(new)*, 09 pathogen Jaccard, 10
projection, 11 abx figures, 12 cytotox figures, 13 curated predictions (now including physchem's
stats/figures), 14 AUROC matrix, 15 the per-endpoint audit unchanged).

**The cycle, found and verified two ways** (reading `10_reference_library_projection.py`'s old guard
directly, and cross-checking output file timestamps, which showed live staleness — step 11 had been
rerun after step 14 last consumed an older version of its output): the old step 09 (projection) hard-
exited unless the old step 14's output already existed, checked *before* it computed anything — but
step 14 needed the old steps 11/12's matrices, and those steps' own scripts needed step 09's
chemical-space background to draw their figures, in the same script that built their matrix. So
`09 → 14 → {11, 12} → 9`, a genuine cycle: no valid run order existed for a from-scratch build.

**The fix has two parts:**
1. **`08_property_matrices.py` is new** — it pulls the matrix-building (and only the matrix-building)
   out of the old steps 10/11/12 into one foundational step, the same "build once, cache, everyone
   re-slices it" pattern step 07 already applies to the bioactivity matrix. Steps 11-13 and 14 now
   read a cache here instead of rebuilding a matrix themselves, so 13 and 14 no longer need any
   figure-drawing script (10/11/12) to have run at all — only 07 and 08.
2. **Step 10's (was 09's) hard exit on step 14's output moved from the top of the script to just
   before the one panel that needs it** (the matrix-14-matched UMAP figure) — every other output of
   that script has no dependency on step 14 and is produced regardless; the one panel now falls back
   to the script's own ranking with a printed note if step 14 hasn't run, rather than blocking
   everything. See that step's section below.

**True run order now: 07 and 08 first (either order) → 09, 10 → 11, 12, 13 → 14 → [rerun 10 for its
bonus panel, optional] → 15.** Nothing is circular any more.

**Folded in at the same time, per explicit decision:** the former standalone `10_physchem_matrix.py`
is gone — its matrix moved into `08_property_matrices.py` and its stats/figures moved into
`13_curated_predictions.py`, since physchem's matrix now feeds only that one downstream analysis
(physchem was dropped from the AUROC matrix on 2026-09-02). `scripts/xx_group_jaccard.py` and
`scripts/xx_non_abx_matrix.py` — confirmed genuine members of this family, reading step 08/09's
outputs — keep their `xx_` prefix for now (not renumbered in this pass) but had their internal path
references updated so they still run.

**Also fixed, while in this code:** a real bug in `figure_cells.json` handling —
`plots_predictor_performance.py` (4 entry points, used by step 13) and `plots_auroc_matrix.
save_overlap_matrix_figure` (step 14) each opened `figure_cells.json` with `"w"` rather than merging
with what was already on disk, which can silently truncate a manifest that another figure in the
same output directory already wrote (this is what step 13's new physchem figures would otherwise
have clobbered, since they now land in the same folder). Both now route through the existing
`plotting_utils.merge_figure_cells` helper, already used correctly elsewhere (e.g.
`save_auroc_matrix_figure`). Also fixed: `plots_abx_projection.py` redefined `_pivot`/`_grid_extent`
instead of importing them from `plots_projection.py` (now imports); two stale docstring headers
mislabelling `eval_predictor_performance.py` and `plots_auroc_matrix.py` as "Step 15"/"Step 16" when
they back steps 13/14; and `plots_tox_projection.py`'s header, stale since cytotox absorbed the
toxicity-projection step on 2026-08-06 but still saying "Step 13" (it backs step 12).

**Step 07's `--write-matrix-csvs` flag removed, same day, on review.** The flag (and the `VARIANTS`
dict / CSV-writing loop behind it) was confirmed unused since being added — no code in the repo ever
read one of the five exported CSVs back, and only 3 of the 5 variants had ever been consumed even
in-memory. See step 07's section above for what was removed and why the reasoning still holds for
any future request to materialise a variant.

## 08_property_matrices.py

The physchem, abx and cytotox full-library score matrices, built once and cached — the same
"read raw prediction CSVs once, cache as a matrix, everyone downstream re-slices it" pattern step 07
already applies to the bioactivity matrix, extended to the other three model families
(`config/physchem_models.csv`, `config/antibiotic_resemblance.csv`, `config/cytotoxicity_models.csv`).
No figures, no chemical-space background — pure matrix-building, engine functions reused unchanged
from `src/eval_property_matrix.py` (physchem, cytotox) and `src/eval_abx_matrix.py` (abx). Writes
`08_physchem_matrix_named.csv` + `08_physchem_endpoint_stats.csv`, `08_abx_matrix_full.parquet` +
`08_abx_matrix_named.csv` + `08_abx_endpoint_stats.csv`, and `08_cytotox_matrix_named.csv` +
`08_cytotox_endpoint_stats.csv` to `output/08_property_matrices/`. Ran in **~13 s** on a warm abx
parquet cache (2026-09-01, this hardware): abx reuses its cache in 0.2 s, physchem and cytotox each
rebuild from their (much smaller than the pathogen library's) raw prediction files in a few seconds.

**Added 2026-09-01 to fix a real circular dependency.** Before this step existed, steps 10, 11 and 12
each built their own matrix from raw predictions *and* drew a figure that needed step 9's chemical-space
background in the same script. That made steps 11/12 (and therefore step 14, which read their matrix
output) depend on step 9 having run — and step 9 itself depends on step 14's output for one panel
(the matrix-14-matched UMAP). The result was a genuine cycle (`9 → 14 → {11,12} → 9`), verified two
ways: reading `10_reference_library_projection.py`'s old guard directly, and cross-checking output
timestamps, which showed live staleness (step 11 had been rerun after step 14 last consumed an older
version of its output). Splitting matrix-building out of 10/11/12 into this one foundational step
breaks the cycle completely: steps 13 and 14 now read straight from here and no longer need any
figure-drawing script (10, 11, 12) to have run first. See the `10_reference_library_projection.py`
section below for the other half of the fix (the one remaining forward reference, now non-blocking).

**Column naming**: `{family}__{model_id}__{column_name}` — physchem, abx, cytotox — the same
three-part shape as the pathogen matrix (step 07, prefix = pathogen code), so all four families join
column-wise on `key`. **Un-normalized only** in all three: no scaling or row normalization, since
that is a decision for whoever joins the blocks downstream (step 13), not this step. Nothing is
dropped or imputed; missing values are counted in each family's `*_endpoint_stats.csv`.

**Physchem (`eos4djh`, 22 descriptors)**: these are *calculations, not predictions* — deterministic
RDKit-family arithmetic, so unlike the other two families there is no model error, no training set
and no leakage dimension; a caption must not describe them as "predicted". All 22 are kept even
though 5 are near-redundant (signed off, not an oversight): `n_rings`, `n_aliphatic_rings`,
`n_aromatic_rings` and `n_saturated_rings` are each exactly the sum of their two
carbocycle/heterocycle components on 100% of rows, and `n_radical_electrons` is 0 for all but ~1
compound in 5,000. Anything fitting a model or reading a correlation matrix off this block should
drop those five first. Upstream spelling kept: `n_aliphatic_heterocyles`/`n_aromatic_heterocyles`/
`n_saturated_heterocyles` are misspelled in Datamol itself; the config matches the real CSV header
rather than correcting it, or the lookup would miss.

**Abx (55 endpoints across 4 models)**: same layout and machinery as the pathogen matrix builder,
factored into its own module (`src/eval_abx_matrix.py`) because it has no organism dimension (a
constant `abx` group code fills the pathogen slot) and its endpoints are almost all discrete (54 of
55 are binary flags or small integer counts; only `abx_score` is continuous). Cached as its own
parquet (`08_abx_matrix_full.parquet`) via `build_full_library_matrix`, same mechanism as step 07's
cache.

**Per-endpoint stats: one shared function for all three families, since 2026-09-01.**
`eval_abx_matrix.endpoint_stats` and `eval_property_matrix.property_endpoint_stats` were near-
duplicate implementations of the same full-library summary (n_unique/min/max/mean/n_nonzero/
pct_nonzero/n_nan), differing only in how they were called (a `sel` DataFrame vs. a bare endpoint-
name list) and that the abx version was missing `median`/`std`. Unified onto
`property_endpoint_stats` — the abx callers (this step and step 11) now pass
`list(matrix.columns)` + `config_csv=<selection path>` instead of a `sel` frame, and
`eval_abx_matrix.endpoint_stats` was deleted. **`08_abx_endpoint_stats.csv` and
`11_abx_endpoint_stats.csv` gained `median` and `std` columns** as a result — additive, verified by
rerunning both steps and confirming every other value unchanged (same 19-endpoint cap count, same
Fisher enrichment numbers). Also fixed while in this code: `08_property_matrices.py` and
`12_cytotox_matrix.py` now import `default.TOX_PREFIX` (`"cytotox"`) instead of hardcoding the
literal string a constant already existed for.

**Cytotox (24 endpoints across 4 models: `eos42ez`, `eos7m30`, `eos3le9`, `eos3dys`)**: every column
here *is* a prediction, unlike physchem — the two must not be treated interchangeably. `eos7m30`'s
own 8 physicochemical columns stay `No` in the config so `eos4djh` remains the single physchem
source and no two blocks can disagree about molecular weight or logP. Carrying the model ID matters:
two models score HepG2 (`cytotox__eos42ez__cytotoxicity_hepg2` vs.
`cytotox__eos3le9__ic50_hepg2_72h_5um`) and the column name is what keeps them distinguishable.


### Pipeline reorganization — 2026-09-02

A second renumbering, on top of the 2026-09-01 split above: **09 becomes bioactivity-only, and the
AUROC matrix moves from slot 14 to slot 10** (previously: 09 pathogen Jaccard, 10 projection, 11 abx
figures, 12 cytotox figures, 13 curated predictions, 14 AUROC matrix, 15 the per-endpoint audit — now:
09 pathogen Jaccard *and* the per-endpoint audit (merged), 10 AUROC matrix, 11 projection, 12 abx
figures, 13 cytotox figures, 14 curated predictions *and* the endpoint confounder check (merged)).
**07, 08, 09, 10, 11, 12, 13, 14 — 8 scripts, down from 9** (the former standalone
`15_bioactivity_endpoints.py` no longer exists as its own file).

**The insight, found by reading the actual code rather than assumed:** the former step 15's headline
metric — each activity endpoint's AUROC against its own pathogen's peer endpoints, i.e. "can this
endpoint uprank the compounds its siblings call active" — is computed by
`eval_predictor_performance.activity_self_performance`, which reads **only the step-07 bioactivity
matrix**. It never touches step 08's physchem/abx/cytotox data at all. It only used to live inside the
(property-predictor) step 13 because it shared target-binarization code with that analysis, not
because of a real data dependency. Splitting it out means step 09 now answers "do a pathogen's
endpoints agree with each other" **twice, at two granularities** (per-pathogen Jaccard, unchanged from
before; per-endpoint Jaccard *and* AUROC, absorbed from the former step 15) using **nothing but step
07** — no property data enters step 09 at all.

**The AUROC matrix (old step 14), by contrast, is *not* bioactivity-only** — it deliberately merges in
2 confounder columns (`cytotox_rank_sum`, `abx_rank_sum`) from step 08 by design, so it doesn't belong
inside a bioactivity-only step either. But it only ever needed steps 07 and 08, so it can move much
earlier than its old slot 14 with no new dependency. Moving it to slot 10 — right after 09, right
before the projection step — means the projection step's one AUROC-matched UMAP panel (see that
step's section) now has its dependency satisfied on the very first pass through the pipeline; no more
optional rerun.

**The confounder check (also former step 15)** — which endpoints' "hits" are really tracking a
physchem/abx/cytotox property rather than genuine pathogen-specific bioactivity — genuinely does need
property-predictor data, so it moved to step 14 instead, run after that step's own predictor-
performance analysis has produced `14_predictor_performance.csv`. It is **not** re-merged into step
09's per-endpoint ranking table this pass ("for the moment", per the same review that decided this
split) — it is written as its own `14_endpoint_confounders.csv` and printed alongside step 09's
ranking for the weakest endpoints, joined at read time rather than at write time.

**True run order now: 07 and 08 first (either order) → 09, 10 → 11 → 12, 13, 14.** Every step needs at
most 07 and 08, plus, for step 11's one bonus panel, step 10 (which now always runs immediately
before it). Nothing is optional or circular.

**Verified at every step, by re-running and diffing against the pre-reorg output before deleting it:**
step 09's per-pathogen Jaccard numbers unchanged; its newly-absorbed per-endpoint Jaccard/AUROC
specificity table (`09_endpoint_quality.csv`, 255 rows) numerically identical on every shared column
to the old `15_endpoint_quality.csv`, including the MTB known-answer check; step 09's directed pairs
CSV (64,770 rows) and per-pathogen summary (12 rows) identical; step 10's `10_auroc_matrix_phylo.csv`
byte-identical to the old `14_auroc_matrix_phylo.csv`; step 11's AUROC-matched UMAP panel confirmed to
complete on the first pass (no "step 10 missing" fallback message); step 14's
`14_predictor_performance.csv` (31,007 rows) and `14_curated_predictor_performance.csv` (840 rows)
byte-identical to the old step 13's; step 14's new `14_endpoint_confounders.csv` numerically identical
on all 255 endpoints (predictor, value and deviation) to the `confounder_*` columns the old
`15_endpoint_quality.csv` carried. `xx_group_jaccard.py` and `xx_non_abx_matrix.py` both re-run clean
against the new paths (the latter's dependency on the AUROC matrix, `output/10_auroc_matrix/
10_auroc_matrix_phylo.csv`, updated accordingly).

**A real, pre-existing bug found and fixed while renumbering:** `src/plots_property_matrix.py`'s
physchem-UMAP figure name was hardcoded `f"13_{method}_top{top_n}_physchem"` — a stale prefix from
before either reorg — so it silently wrote `13_umap_top1000_physchem.png` into step 14's own output
directory instead of `14_umap_top1000_physchem.png`, undetected because nothing checked the figure
name against its own step number. Fixed to `f"14_{method}_top{top_n}_physchem"` and re-verified.

## 09_bioactivity_endpoints.py

### Part 1 — per-pathogen top-1000 Jaccard

Asks whether endpoints of the SAME pathogen pick out the same compounds more than endpoints of
DIFFERENT pathogens do. For every pair of the 260 endpoint columns, the Jaccard overlap of their
top-1000 highest-scoring compounds, aggregated per pathogen into a same-vs-different box pair.
Figure in `src/plots_matrix_analyses.py`; writes the 307x307 baseline Jaccard matrix (reused
downstream), a per-pathogen summary CSV, and PNG + PDF to `output/09_bioactivity_endpoints/`.
**Key decisions, for review:**
- **Cutoff: top 1000** of 1,355,109 compounds per column — user-directed.
- **Scope: the 15 curated pathogens of interest** (`config/pathogens_of_interest.csv`), replacing the
  former `min<K>`-endpoint thresholds on 2026-08-07 (user-directed). **260 of 307 columns, 15 of 58
  pathogens**; the other 43 pathogens / 47 columns are removed *entirely*, ceasing to be
  different-pathogen partners too, not merely losing their own box. The 58th pathogen is
  *E. faecalis*, new with `eos3f8h` on 2026-08-11 and **excluded by design** — it is absent from
  `config/pathogens_of_interest.csv`, so its one column is dropped here (user-confirmed).
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
- **Only the baseline variant is computed and published — trimmed from 3 to 1 on 2026-09-01.**
  `zscore_l2rownorm` and `rankpct_l1rownorm` (row-normalized on top of z-score/rank-percentile
  columns) used to get their own full Jaccard computation, summary CSV and figure alongside
  `baseline`, but nothing downstream ever read either one — only `baseline` (this step's own
  per-organism-class aggregation, Part 2 below). They were a genuine robustness check (does row-normalizing change which
  pathogen looks "specific"?), not dead weight, but at ~2x the compute of baseline alone
  (**~2m07s → ~74s** on a from-scratch build, this hardware) for 2 figures nobody consumed, they were
  cut on review. Recoverable from `git log` if that question needs revisiting.
- **A cheaper, separate check is still run and kept: does *column* scaling alone ever change a top-N
  Jaccard set?** Unrelated to the row-normalized variants above — it never row-normalizes anything.
  Top-N Jaccard depends only on each column's own internal ranking, and both of step 07's column
  scalings (z-score, rank-percentile) are strictly increasing per column, so in exact arithmetic
  neither could change any column's top-1000 set. The script **asserts** this at runtime rather than
  assuming it, over the same cached matrix it already holds (no extra I/O), and the assertion earns
  its keep: baseline == rank-percentiled holds exactly, but baseline == z-scored comes back
  **False** — `(x - mean) / std` in float32 reorders near-tied values in exactly **one column of 300**
  (`lmajor__eos60mw__leishmania_mlp`, 155 of its top-1000 members shift), touching **156 of 90,000
  cells**, max absolute difference 0.001032 and mean 0.000599 over those cells. Re-verified on the
  300-endpoint rebuild: still that one column and no other, and baseline == rank-percentiled differs
  in **0 of 90,000** cells. (At 260 endpoints it was 138 of 67,600.) **On the 307-endpoint rebuild the
  two check verdicts are unchanged** (z-scored `False`, rank-percentiled `True`, now over
  307x307 = 94,249 cells), but the script prints only those verdicts, so the per-cell counts above
  were **not** re-quantified at 307 and still describe the 300 build. Note the difference floor is
  ~1/2000 = 5e-4, the granularity of a Jaccard over two ~1000-element sets, so any comparison
  reporting differences below that is reading float noise — e.g. from a CSV round-trip rather than
  the in-memory arrays the assertion uses.
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

### Part 2 — per-organism-class Jaccard (absorbed from the former standalone `xx_group_jaccard.py`, 2026-09-02)

Part 1 one level up: the same top-1000 Jaccard question asked of an **organism class** rather than a
pathogen. Do all the Gram-negative endpoints agree with each other more than with the Gram-positives,
the fungi and the rest? That is the level at which cross-organism transfer would show up. A pure
re-aggregation of Part 1's own in-memory Jaccard matrix and summary — nothing is recomputed, nothing
is re-read from disk, and it runs in seconds.

**Scope: the 15 pathogens of interest** (user-directed, matching Part 1), **260 endpoints in 6
classes**.

**Four of the six classes are a single organism, and that has to reach every caption.** Mycobacteria
is *M. tuberculosis* alone, Helminths *S. mansoni* alone, Fungi *C. albicans* alone (21 of 21
endpoints), Protozoa *P. falciparum* alone (64 of 64). For those four, "same class" and "same
pathogen" are the same set of pairs, so their `same_median` is a verbatim copy of Part 1's and their
same-class-different-pathogen box is empty. **Only Gram-negative (8 organisms) and Gram-positive (3)
carry information the per-pathogen figure does not already show.** Rather than caveat this in prose,
the script *asserts* the four copies against Part 1's own in-memory summary — if the class join were
wrong those medians would stop matching, so the degeneracy doubles as the correctness check.

**Three boxes per class, not two.** `same` is dominated by within-pathogen pairs that Part 1 already
published, so the box carrying the new information is **`same_excl_same_organism`** — same class,
different pathogen. `same_excl_same_model` is carried too, for continuity with Part 1 and because it
matters *more* here: a class pools several models, so its internal agreement can be one model
agreeing with itself across two organisms.

**Result:** Gram-negative same-class median 0.0116, falling to **0.0106** once within-pathogen
pairs are removed, against a different-class median of 0.0000 — its endpoints genuinely agree across
species. Gram-positive is weaker (0.0020 -> 0.0015 vs 0.0005).

**No threshold** — sorted by specificity, never filtered, as in Part 1. The figure is plain
matplotlib, matching `pathogen_jaccard_figure`, so the two read as one diagnostic family; a class with
one organism simply draws no middle box, which is the visual form of the degeneracy. Engine
`src/eval_group_jaccard.py`, figure `group_jaccard_figure` in `src/plots_matrix_analyses.py`. Writes
`09_group_jaccard_top1000_summary.csv` (6 rows), `09_group_jaccard_pairs.csv` (directed pairs), and
`09_group_jaccard_top1000` (png/pdf, diagnostic, no `figure_cells.json` entry).

**Shared-code change (from when this was `xx_group_jaccard.py`):** `eval_correlations.
pathogen_metric_boxes` / `pathogen_metric_summary` gained an optional `exclusions` mapping (default
`{"same_excl_same_model": "same_model"}` = Part 1's behaviour exactly), which `eval_group_jaccard`'s
class-level aggregators reuse rather than duplicate. Verified as a no-op on Part 1's own output.

### Part 3 — per-endpoint agreement audit (absorbed from the former standalone `15_bioactivity_endpoints.py`, 2026-09-02)

Per-**endpoint** agreement audit: can each endpoint uprank the compounds its own pathogen's other
endpoints call active? The per-pathogen Jaccard figure above (Part 1) and this step's own AUROC
self-performance figures (below) both report **one number per pathogen**, so an individual endpoint
that behaves badly is averaged into its pathogen's box. This section regroups the same pair
statistics **by endpoint** instead.

**AUROC self-performance is computed fresh here, bioactivity-only.** Each activity endpoint's
**raw, un-binarized** score is used as a predictor against every OTHER endpoint's top-1000
binarized version as the target — 307 x 307 AUROCs, via
`eval_predictor_performance.run_activity_self_performance`, which reads **nothing but the step-07
matrix**. This is the key architectural fact behind the 2026-09-02 reorg (see the "Pipeline
reorganization" note above): the self-performance engine never touches step 08's physchem/abx/cytotox
data, so it belongs here, in a bioactivity-only step, rather than beside the property-predictor
analysis in step 14. Writes `09_activity_self_performance.csv` (all 307 targets) and
`09_pathogen_subset_self_performance.csv` (the 15 pathogens of interest, consensus models collapsed
to one column per model — 260 endpoints -> 70, 12 consensus / 58 single). Self-pairs (an endpoint
against its own binarization) are 1.0 by construction — verified min 0.999952, max 1.000000 over all
307 — kept in the CSVs as a correctness check, excluded from both figures below.

**Two more figures, moved here from the former property-predictor step (old step 13) since their
data is bioactivity-only. Both now read the SAME 15-pathogen, consensus-collapsed subset (70
endpoints, 4,900 pairs, 4,830 after removing self-pairs) — 2026-09-02, user-directed, narrowed from
the original all-56-organism version so "cross-pathogen" means one of the other 14 priority
pathogens, not any of the other 41 organisms in the full selection, and so a pathogen's box can never
be dominated by one model's dozens of correlated sub-assays (the same reason the consensus collapse
exists at all — see `pathogen_subset_endpoints`):**
- **`09_performance_activity_by_organism`** — one box per pathogen (15) over its endpoints' AUROCs
  against every OTHER endpoint's binarized version, points coloured by whether the target belongs to
  the same pathogen (460 same-pathogen vs 4,370 cross-pathogen pairs). Points are subsampled at 400
  per colour per box (seeded with `RANDOM_SEED`); the box itself uses all values.
- **`09_performance_pathogen_subset`** — the same 70-endpoint pool, one box per **endpoint** rather
  than per pathogen, ordered by pathogen then by median within it, labelled `{pathogen} - {endpoint}`
  — the finer-grained view that lets a pathogen's consensus score and its individual assay endpoints
  sit side by side rather than being pooled into one box.

**Minimum endpoints: >5** (`MIN_ENDPOINTS = 6`). A pathogen with 5 or fewer endpoints cannot support
a per-endpoint peer distribution — with 3 endpoints each one is judged on 2 peers. Keeps **12 of the
15** pathogens of interest / **255 endpoints**; drops *N. gonorrhoeae* (3), *Campylobacter* (1) and
*H. pylori* (1). Per the contract of `pathogens_of_interest_nodes` and `multi_column_pathogen_nodes`,
this removes those three **entirely**, so they stop being cross-pathogen partners too: the
`diff_median` here is against the other **11** pathogens. That is narrower than Part 1's other-14,
which is itself narrower than the all-57 comparator the pre-2026-08-07 figures used. A caption must
say so — "specific to this pathogen" here means "relative to the other 11 priority pathogens".

**Consensus columns are kept, and are not fair comparators.** A `consensus_score` column aggregates
the very sub-endpoints it is scored against, so its agreement with them is inflated by construction.
It is kept (flagged `is_consensus`) as its pathogen's within-block ceiling rather than dropped or
collapsed. This differs deliberately from `pathogen_subset_endpoints` above, which collapses each
model to its consensus column for the two performance figures: that reduction is what this analysis
exists to look inside. The behaviour is as expected — **11 of the 12 consensus columns are the
single strongest endpoint of their pathogen** — which is also the check that the ranking is oriented
correctly. The exception is reported: `eos8v1a` (*S. mansoni*) ranks 4th weakest of 11.

**Same-model peers are INCLUDED** in every headline statistic, matching Part 1's "each column
against all the others". Where one model supplies most of a pathogen's endpoints this measures a
model agreeing largely with itself, so every row carries `n_same_model_peers` and the pathogen
summary carries `n_single_model` (endpoints with no cross-model peer at all). Two blocks are
single-model outright — Enterobacter (6 endpoints, all `eos9bpi`) and *E. faecium* / *S. pneumoniae*
(9 each) — and Enterobacter topping the specificity ranking is the same artifact already flagged for
Part 1's figure, not a finding. The directed pairs CSV keeps its `same_model` flag, so the
same-model-excluded view can be recovered without recomputing either matrix.

**No threshold, no flag column, no filtering.** The outputs are rankings; deciding which endpoints
are genuinely off is a scientific judgement. The one count reported, `n_below_chance_peers`, is a
definition — AUROC < 0.5 means the peer's actives were ranked below its inactives — not a chosen
cutoff. `TAIL_N = 40` in `plots_endpoint_quality` is a **display** limit on the diagnostic figure
only: the axis label states "weakest 40 of 255" and the full ranking is in the CSV.

**Median, not mean**, throughout, following Part 1's `same_median` / `diff_median` convention. The
two differ enough to matter: MTB's `chembl_single_point_7` is 0.4605 as a mean and 0.4663 as a
median, and `eos43d6:consensus_score` is 0.7951 vs 0.8471. Both below-chance MTB values are asserted
in the script so the wrong statistic cannot be reported silently.

**Checks, all of which would otherwise fail silently rather than raise:**
- **Endpoint counts** are asserted against `config/08_endpoint_selection.csv` joined to
  `config/pathogens_of_interest.csv`, not against a hardcoded list.
- **Known answer:** every one of MTB's 40 endpoints must match a direct `groupby` on this step's own
  `09_activity_self_performance.csv` using *its own* organism columns — an independent path from the
  config-organism → code → named column → endpoint-key mapping this step builds. Currently exact
  (max deviation 0.0).
- **Cross-metric:** Jaccard and AUROC specificity must be positively rank-correlated, since they
  measure the same thing by different means. Spearman ρ = **0.482** (p = 3.0e-16). A near-zero value
  would mean the join between the two naming conventions
  (`{pathogen}__{model}__{column}` for Part 1's Jaccard matrix, `{model}:{column}` for the
  self-performance frame) is wrong and every endpoint has someone else's Jaccard.

**What it currently shows:**

**8 of 255 endpoints rank their own pathogen's actives below chance.** MTB's two known columns are
there, but so are six others, and the weakest is not MTB's: `eos21dr:chembl_dose_response_6`
(*A. baumannii*) at 0.345. The others are `eos6wb7:chembl_single_point_4` (0.443),
`eos8jx6:chembl_dose_response_4` / `_3` / `chembl_single_point_4` (*C. albicans*, 0.481 / 0.489 /
0.496) and `eos81zy:chembl_single_point_0` (0.498). Every one is a `chembl_*` sub-assay column, none
is a consensus column, and all eight sit in blocks dominated by same-model peers.

At pathogen level the ranking is *C. albicans* (0.615) < *S. aureus* (0.624) < **MTB (0.653)** <
*P. aeruginosa* < ... < *S. mansoni* (0.922). So MTB is weak but **not uniquely so** — the finding
that motivated this analysis generalises to at least two other pathogens, and MTB's position looks
worse on Jaccard (where it is last) than on AUROC (where it is third from last). The two metrics
disagree because Jaccard is a hard top-1000 set overlap and AUROC uses the whole ranking; MTB's
40 endpoints rank compounds consistently enough to beat chance while almost never agreeing on which
exact 1000 are the best.

**47 endpoints have a negative AUROC specificity** and 37 a non-positive Jaccard specificity — they
agree with other pathogens' endpoints at least as much as with their own. That is reported, not
filtered, and is the second axis of the specificity scatter panel.

**Outputs:** `09_endpoint_quality.csv` (255 rows, weakest first — the confounder columns the former
step 15 carried here are gone; step 14 now writes its own confounder table, see that step's section),
`09_endpoint_pairs.csv` (64,770 directed rows carrying both metrics and the `same_model` flag),
`09_pathogen_endpoint_summary.csv` (12 rows), and five figures in `png/` + `pdf/`
(`09_endpoint_upranking`, `09_endpoint_specificity`, `09_performance_activity_by_organism`,
`09_performance_pathogen_subset` — all four `BasePlot`s with a `figure_cells.json` footprint — plus
`09_endpoint_ranked_tail`, a plain-matplotlib **diagnostic** sized in inches by its row count, since
it carries one label per endpoint and cannot stay legible at page width).

## 10_auroc_matrix.py
Collapses each organism's activity endpoints into **one score per organism**, then draws the AUROC
matrix: **15 organism rows x 17 columns** (the same 15 organism aggregates, then the two merged
property rank-sums, cytotox and abx — see "Merged predictors, dropped physchem" below for how the
matrix got here from an earlier 27-column/405-cell version). Each cell carries its printed AUROC.
Assembly and computation in `src/eval_auroc_matrix.py`, figure in `src/plots_auroc_matrix.py`; writes
to `output/10_auroc_matrix/`. Replaces an earlier per-endpoint 61 x 73 version, which was too dense
to read: a pathogen contributed up to 13 correlated endpoints, so its band said more about how many
assays it has than about the pathogen.
**Depends on nothing but steps 07 and 08** (2026-09-01 — see the "Pipeline reorganization" note
after step 07): the abx and cytotox columns are read straight from `08_abx_matrix_named.csv` /
`08_cytotox_matrix_named.csv`, not from steps 12/13's own output, which now only draw figures. This
is what breaks a real circular dependency the pipeline used to have (step 11 needs this step's
output for one panel; before this change, this step needed steps 12/13, which needed step 11's
background — a cycle).
**The merge:** an organism's endpoint columns are scaled to **percentile rank** within the full library
(`ORGANISM_MERGE_METHOD = "rank_pct"`) and then **averaged** (`ORGANISM_MERGE_AGG = "mean"`). Scaling
*before* aggregating is the point — raw scores from different models sit on unrelated ranges, and
averaging them directly would weight whichever endpoint has the widest one. Endpoints are the
consensus-collapsed 15-pathogen set (61 endpoints; 1-13 per organism).
**Two properties of the merge, recorded rather than corrected** — both follow directly from merging the
endpoints as selected: (1) **five organisms have exactly one endpoint** (Campylobacter, Enterobacter,
*E. faecium*, *H. pylori*, *S. pneumoniae*), so nothing is merged and their score IS that endpoint's
percentile rank — their row is not the same kind of quantity as *E. coli*'s 12-endpoint mean; (2) a ChEMBL
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
### Second view: shared actives (`10_overlap_matrix_top{N}`)
The same 15 x 27 axes and annotation tracks, a different quantity: **how many of the row organism's
top-N actives fall in the column's own top N**. AUROC asks "does this predictor RANK the
row's actives highly across the whole library"; this asks how many of the very same molecules it puts
at the top. A predictor can do the first well without doing the second.

> **Superseded 2026-09-01:** the 100 and 10000 panels described below were dropped — the pipeline now
> draws only the top-1000 cutoff, and `OVERLAP_MATRIX_CUTOFFS` no longer exists in `src/default.py`.
> Kept here for the reasoning behind the non-uniform bins and the multi-cutoff comparison itself.

**Drawn at three cutoffs** (`OVERLAP_MATRIX_CUTOFFS` = 100 / 1000 / 10000), each its own CSV and
figure, because where the line is drawn changes the picture completely — this is the finding, not
noise. Over the 210 off-diagonal bioactivity cells:

| cutoff | off-diagonal range | median | cells sharing nothing | median as % of N | chance (`N²/n_total`) |
|---|---|---|---|---|---|
| top-100 | 0-59 | 0 | **164 of 210** | 0.00% | 0.007 |
| top-1000 | 0-724 | 3 | 58 of 210 | 0.30% | 0.7 |
| top-10000 | 7-4119 | 410 | **0 of 210** | 4.10% | 73.8 |

At top-100 most organism pairs share nothing at all; by top-10000 every pair shares something and the
median pair shares 410 compounds. **The three panels are not interchangeable**: chance overlap grows
as N², so at top-10000 a cell showing a handful of shared actives is *below* chance, while the same
count at top-100 is far above it. A caption must say which cutoff it is showing and what chance is
there.

Bin boundaries are rescaled per cutoff by `plots_auroc_matrix.overlap_bins(top_n)`, which scales
`OVERLAP_MATRIX_BINS` by `top_n / 1000` so each bin stays the same **fraction** of the cutoff and the
three panels stay comparable. Two boundaries are pinned rather than scaled: `0`, and `1` — the
sharing-nothing/sharing-something boundary, which is meaningful at every cutoff and would otherwise
be scaled away (at top-100 the second bin would land on 1 and collide with it). After rounding the
remaining boundaries are forced strictly increasing, since the naive scale collides at top-100
(10 -> 1, 25 -> 2.5). `overlap_bins(1000)` returns `OVERLAP_MATRIX_BINS` unchanged, and the top-1000
CSV is **byte-identical** to the single-cutoff version it replaces.

**The AUROC matrix is not triplicated.** It uses the whole ranking, so its cutoff only decides which
compounds count as positives; `ACTIVITY_BINARIZE_TOP_N` remains its single value.

**Cost:** the property CSVs are read **once** for all three cutoffs (`eval_auroc_matrix.predictor_tops`
takes a list and computes every cutoff from the column already in memory). Reading them per cutoff
would have tripled the step's I/O over several GB for no new data.

**Renamed outputs:** the single-cutoff `10_overlap_matrix.csv` / `.png` / `.pdf` were deleted on
2026-08-10 after confirmation. They are fully reproduced by the `_top1000` files — the CSV and PNG
were byte-identical, the PDF differed only in its embedded creation timestamp. Anything referring to
the old names needs updating to `10_overlap_matrix_top1000`.
**A raw count, not Jaccard.** Both sets have exactly 1000 members, so `J = i / (2000 - i)` is a
monotone re-expression of the same number and orders the matrix identically — the count is the one a
reader can act on ("724 of the 1000 shared"). **The measure is symmetric**, so the bioactivity block is
a symmetric matrix, unlike AUROC's.
**Bins are NON-UNIFORM** (`OVERLAP_MATRIX_BINS` = 0, 1, 10, 25, 50, 100, 200, 400, 750 at top-1000,
rescaled per cutoff as above), tuned to a heavily skewed distribution: at top-1000 the off-diagonal
counts run 0-724 with a **median of 3**, and two random 1000-compound sets out of 1,355,109 would
share ~0.7 by chance. Uniform bins would put most of the matrix in one class. No bin holds more than a
third. The colourbar uses **uniform** spacing so the narrow low bins stay readable.
**Colour was SEQUENTIAL cobalt** (a single hue), matching the overlap heatmaps in step 09 and the EU
OpenScreen validation, until the 2026-08-31 change below. Still sequential, not diverging, either
way: a count's neutral is 0, at the end of the scale rather than in the middle, so there is nothing
to diverge around.

**Multi-hue spectrum — 2026-08-31 (user-directed, iterated twice).** The original ramp was a plain
white → `cobalt` 2-stop gradient (`plotting_utils.sequential_cmap`), which never gets darker than the
named hue itself (`#457B9D`) — over an 8-bin scale the mid-to-high bins came out as similar
pale-to-medium blues, hard to tell apart. A single-hue white→cobalt→dark-navy 3-stop version was
tried first and rejected (**not what was wanted** — "I want you to use more than one color, it is a
spectrum"). **Now uses `plotting_utils.spectrum_cmap`** over `default.OVERLAP_MATRIX_SPECTRUM =
("amber", "turquoise", "cobalt", "periwinkle")`: white → amber → turquoise → cobalt → periwinkle, a
genuine multi-hue ramp. The four hues are ordered by **decreasing perceptual luminance** (amber 0.77
→ turquoise 0.64 → cobalt 0.45 → periwinkle 0.41), so the scale still darkens monotonically and reads
unambiguously as a magnitude gradient even though hue changes at each stop — the same principle
behind well-known multi-hue sequential palettes (e.g. `YlGnBu`), and NOT the same thing as a
qualitative/rainbow map: stylia's `SpectralColormap` cycles hues with no such ordering, which is
exactly why it was rejected for step 10's (diverging) AUROC scale — see that figure's colour history
above. **This is a deliberate departure from "every set-overlap figure in the paper reads on one
hue"** — step 10's overlap matrices no longer match step 09/05's plain cobalt sequential style, and a
caption should say so if the two are shown together.

**Sampled from 0.0, not the small legibility-margin offset other sequential scales in this repo use**
(`plotting_colors.count_shades`'s `headroom`): the 0-shared-actives bin is meant to read as true,
blank white, since 0 is a common, meaningful value here, not an edge case to keep visible — 76.2% of
bioactivity off-diagonal cells are exactly 0 at top-100, 26.7% at top-1000, 0% at top-10000 (so the
white bin is invisible in that panel by construction, not a bug). All three cutoffs share this
colormap function, so the change applies to all three, not only top-1000.

**Linear colour scale for step 10 only — 2026-09-02 (user-directed).** The discrete, non-uniform bins
above are drawn with `spacing="uniform"` — every bin gets equal visual height regardless of its
numeric width — which is exactly what a LOG-scale colourbar looks like (values 1/10/100/750 evenly
spaced) even though nothing is actually log-transformed. Flagged as confusing, so **step 10's
`10_overlap_matrix_top1000_{hclust,phylo}` now use a plain continuous LINEAR scale instead**:
`OverlapMatrixPlot(..., continuous_color=True)` swaps `discrete_overlap_cmap`'s `(ListedColormap,
BoundaryNorm)` for `(plotting_utils.spectrum_cmap(OVERLAP_MATRIX_SPECTRUM),
matplotlib.colors.Normalize(0, top_n))` — same four-hue palette, same 0-is-white anchor, but colour
now varies smoothly and proportionally with the actual value, and the colourbar ticks are matplotlib's
own evenly-value-spaced defaults rather than the skewed bin boundaries.
**The non-abx section below (see "Non-abx robustness check") now uses this same `continuous_color=True`
scale too, as of the 2026-09-02 fold that absorbed it into this step** — it used to keep its own three
cutoffs at the discrete default (`continuous_color=False`), on the reasoning that its non-uniform bins
were the deliberate fix for this same skew and not what was raised here. Folding it in to draw the
exact same plots as this section overrode that: both now share the class
(`plots_auroc_matrix.OverlapMatrixPlot`) *and* the scale.
**Consequence, expected rather than a bug:** with a median of 2-3 shared actives (of 1000) and a real
max of 724, a strictly linear 0-1000 scale renders almost every off-diagonal cell pale — only the
tightly-correlated Enterobacterales block (400-724) reads as strongly coloured. That IS what "normal"
(linear) looks like on data this skewed; the discrete bins existed specifically to counteract it.
**Cell-annotation text colour** (white vs ink) is picked from the colormap's own relative luminance
crossing 0.5 (`plots_auroc_matrix._dark_threshold`, ≈682 of 1000 for this palette+cutoff) rather than
a guessed fraction of `top_n` — `spectrum_cmap` mixes hues of different intrinsic lightness, so
"darker" and "higher value" don't track each other as simply as they would for a single-hue ramp.

**The two Δ-marked property columns — SUPERSEDED 2026-09-02.** `Δmw`, `Δclogp`, `Δtpsa` used to flag
the three physchem columns whose own "top N" was defined as furthest-from-median rather than
highest-value. Physchem was dropped from the matrix entirely (see "Merged predictors, dropped
physchem" below), so this marker no longer applies to anything and the code that drew it has been
deleted, not just left unused. The two property columns that remain (`cytotox`, `abx sim`) are both
one-directional highest-value rank-sums and need no such marker.

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
else. Follows the convention in step 09 and `ActiveOverlapHeatmapPlot`. **The true values, diagonals
included, are still in both CSVs** — blanking is a display choice, and the build-time diagonal assert
still runs on the unmodified matrix.

**Three build-time checks**, each of which would otherwise fail silently and still look plausible: every
aggregate column's mean must be ~0.5 (a mean of percentile ranks must be; catches scaling along the wrong
axis), the diagonal must be 1.0 across all 15 (an organism against its own binarization; doubles as a
transpose test), and no cell may be missing. The script exits rather than drawing if any fails.
**Nothing is excluded any more.** This previously dropped 3 `eos9ivc` *M. tuberculosis* endpoints that
were absent from a step-07 cache built before their prediction file was staged; the rebuilt parquet
carries all three, and step 10 now reports 405 cells with none missing. The guard that caused it is
still live and still worth knowing about: step 10 intersects the selection with the columns actually in
the cache (`available=cached`) and **silently proceeds with fewer endpoints** rather than failing, so a
stale cache shows up only as a lower endpoint count in the log.

### Row-order comparison — 2026-08-31, simplified 2026-09-01, simplified again 2026-09-02

**Current state (2026-09-02): `hclust` is GONE too — `phylo` is the only order the script produces.**
The `hclust`/`phylo` comparison below (kept from the 2026-09-01 simplification) did its own job in
turn: `phylo` is what was kept. What the script draws now: **2 figures total** —
`10_auroc_matrix_phylo` and `10_overlap_matrix_top1000_phylo` — plus the `10_phylo_dendrogram`
diagnostic (`10_hclust_dendrogram` was dropped alongside `hclust`). `10_row_order_comparison.csv`
still carries both a `baseline_position` and an `hclust_position` column for reference (cheap, no
heatmap drawn from either), so both retired orders remain inspectable. See "Merged predictors, dropped
physchem — 2026-09-02" below for the column-count change (27 → 17) that ran alongside this. The rest
of this subsection is the **history of how the `phylo`-only choice was reached** — it describes the
fuller comparisons that ran on 2026-08-31 and 2026-09-01 and are not what the script produces today.

**`config/organism_taxonomy.csv` verified against live NCBI Taxonomy — 2026-09-01.** The table was
originally hand-compiled from general microbiology knowledge, not queried from a database, and was
flagged at the time as needing a spot-check. It has now been checked directly against NCBI's
E-utilities (`esearch` + `efetch`, `db=taxonomy`, one query per organism, taxids: *A. baumannii* 470,
*Campylobacter* 194, *E. coli* 562, *Enterobacter* 547, *H. pylori* 210, *K. pneumoniae* 573,
*N. gonorrhoeae* 485, *P. aeruginosa* 287, *E. faecium* 1352, *S. aureus* 1280, *S. pneumoniae* 1313,
*M. tuberculosis* 1773, *C. albicans* 5476, *P. falciparum* 5833, *S. mansoni* 6183).

**7 of the 75 cells (15 organisms x 5 ranks) were wrong and have been corrected:**

| organism | rank | had | corrected to |
|---|---|---|---|
| *Campylobacter* | taxonomic_class | Campylobacteria | **Epsilonproteobacteria** |
| *H. pylori* | taxonomic_class | Campylobacteria | **Epsilonproteobacteria** |
| *S. aureus* | order | Bacillales | **Caryophanales** |
| *M. tuberculosis* | taxonomic_class | Actinomycetia | **Actinomycetes** |
| *C. albicans* | taxonomic_class | Saccharomycetes | **Pichiomycetes** |
| *C. albicans* | order | Saccharomycetales | **Serinales** |
| *S. mansoni* | order | Diplostomida | **Strigeidida** |

`Caryophanales`/`Pichiomycetes`/`Serinales` are genuinely recent NCBI reclassifications (Bacillales
and Saccharomycetes/Saccharomycetales are still the names in most textbooks and much current
literature) — this is exactly the *C. albicans* taxonomy-in-flux caveat flagged when the table was
first written, now confirmed rather than merely suspected.

**None of the 7 corrections changes any figure.** Verified directly: `phylogeny_organism_order` and
`phylogeny_class_linkages`'s linkage matrices are **byte-identical** before and after the fix, and a
full pipeline re-run produced a byte-identical `10_auroc_matrix_phylo.csv`. This holds for a
structural reason, not by luck: *M. tuberculosis*, *C. albicans* and *S. mansoni* are singletons in
their `organism_class` (nothing to order against), and the two corrections inside the 8-member
Gram-negative group and the 3-member Gram-positive group both replaced one "different from its
neighbours" label with another — *Campylobacter*/*H. pylori*'s wrong class still differed from the
Gammaproteobacteria/Betaproteobacteria of the rest of the group either way, and *S. aureus*'s wrong
order still differed from the other two organisms' `Lactobacillales` either way. The table is now
correct **and** the dendrogram/heatmap it drives were never actually wrong.

Rows kept the same `ORGANISM_CLASS_ORDER` grouping throughout; two alternative WITHIN-class orders
were added for comparison against the (then-)existing alphabetical default, applied to **both** the
AUROC matrix and all three overlap-matrix cutoffs (12 figures total: 3 orders x (1 AUROC + 3
overlap)). None of the three was meant to be a "final" order on its own — all three were written and
drawn so they could be read side by side, using the SAME `hclust`/`phylo` organism order for the
AUROC matrix and every overlap cutoff, so a given organism sat in the same row across all four matrix
families.

- **`10_auroc_matrix_hclust.*`** — hierarchical clustering (`scipy.cluster.hierarchy`, no new
  dependency) on each organism's full 27-column row profile, correlation distance (1 − Pearson r,
  `ORGANISM_HCLUST_METRIC`) with average linkage (`ORGANISM_HCLUST_LINKAGE`), leaf order tidied by
  `optimal_leaf_ordering`. Clusters by profile SHAPE — which predictors rank an organism high vs
  low — not by overall AUROC level. **Visibly recovers a block**: *E. coli*, *K. pneumoniae*,
  *A. baumannii*, *Enterobacter*, *P. aeruginosa* and *N. gonorrhoeae* land adjacent at mutual AUROC
  0.94–1.00, with *H. pylori* and *Campylobacter* (the two Campylobacterota, weakest links to the
  rest at 0.47–0.88) pushed to the group's edges.
- **`10_auroc_matrix_phylo.*`** — within each class, organisms are ordered by the tuple (phylum,
  taxonomic class, order, family, genus) from `config/organism_taxonomy.csv`, so organisms sharing a
  family or genus end up adjacent. **Originally described here as "a plain lineage sort, not a
  computed tree" — superseded 2026-08-31**, the same day, once the `10_phylo_dendrogram` diagnostic
  was added: `phylogeny_organism_order` now derives from `eval_auroc_matrix.phylogeny_class_linkages`
  (a real ultrametric tree built from the taxonomy table, see the row-order-comparison section
  below), not an independent sort — the two cannot disagree. Two divergences from a naive
  alphabetical/genus read are worth knowing: *Campylobacter* and *H. pylori* are placed together
  (both Campylobacterota, a phylum split out of Proteobacteria) ahead of the rest of the
  Gram-negatives, and *N. gonorrhoeae* (Betaproteobacteria) sits apart from the six
  Gammaproteobacteria that make up the rest of that class.
  **`config/organism_taxonomy.csv` was verified against live NCBI Taxonomy on 2026-09-01** — see the
  dedicated note in the row-order-comparison section below for the method and what it found.
- **Only Gram-negative (8 organisms) and Gram-positive (3) actually reorder.** The other four
  classes are singletons (Mycobacteria, Fungi, Protozoa, Helminths) and are — and must be — in the
  identical row position across all three CSVs; `10_row_order_comparison.csv` confirms this
  (`baseline_position == hclust_position == phylo_position` for all four).
- **`10_row_order_comparison.csv`** carries one row per organism with its class and its 0-based
  position under each of the three orders (computed once, on the AUROC matrix's row set, and reused
  for every overlap cutoff), so the three don't need opening side by side to compare.
- **Overlap matrices**: `10_overlap_matrix_top{100,1000,10000}_{hclust,phylo}.csv` reuse
  `eval_auroc_matrix.reorder_bioactivity_axes` unchanged — the overlap matrix has the identical
  `(organisms x organisms+predictors)` shape as the AUROC matrix, so no new reordering code was
  needed, only re-running the existing per-cutoff invariants (diagonal == `top_n`, bioactivity block
  symmetric) after reordering to confirm nothing was corrupted by the permutation.
- Every reordered matrix re-runs its baseline invariant checks before drawing (`diagonal_check` for
  the AUROC matrix, the diagonal/symmetry asserts for each overlap cutoff), and all of them route
  through the now-shared `plotting_utils.merge_figure_cells` so `figure_cells.json` accumulates all
  12 entries instead of the last writer truncating the manifest — the same bug class already fixed
  for step 12's abx grid, now caught here before it could bite (`save_auroc_matrix_figure` previously
  wrote the manifest with `"w"` on every call).
- **`10_hclust_dendrogram.*` and `10_phylo_dendrogram.*`** (2026-08-31) — diagnostics, not
  `BasePlot`s on the 3 cm cell grid (same departure as step 09's `pathogen_jaccard_figure`), one
  horizontal dendrogram per class with 2+ organisms (only Gram-negative and Gram-positive draw
  anything), sized in inches by leaf count, drawn by the shared `plots_auroc_matrix.
  save_dendrogram_figure`. Neither is entered into `figure_cells.json`.
  - **`10_hclust_dendrogram`**: x-axis = 1 − Pearson correlation distance. Reads the exact same `Z`
    linkage that ordered the `hclust` heatmap rows (`eval_auroc_matrix.hierarchical_class_linkages`,
    factored out of `hierarchical_organism_order`), so the tree shown IS the tree the heatmap was
    permuted by.
  - **`10_phylo_dendrogram`**: x-axis = taxonomic rank distance, an integer 1–5 = how many of the 5
    ranks (phylum, taxonomic_class, order, family, genus) two organisms do NOT share — sharing
    everything down to family merges at 1, sharing nothing (different phylum) merges at 5. This is
    an exact ultrametric by construction (a real nested classification, not a fitted distance), so
    `eval_auroc_matrix.phylogeny_class_linkages` assembles the one tree it implies directly
    (`_cartesian_linkage_from_order`) rather than approximating it with `scipy.cluster.hierarchy.
    linkage` on a pairwise matrix, which would face unnecessary tie-breaking — every pair among the
    three Enterobacteriaceae genera here (*E. coli*, *Enterobacter*, *K. pneumoniae*) is exactly
    equidistant. **`phylogeny_organism_order` was refactored to derive from this same linkage**
    (previously an independent tuple-sort with identical output), so the `phylo` heatmap and its
    dendrogram are now guaranteed to agree — verified: the row order is byte-identical to the
    already-documented one above (Gram-negative: Campylobacter, H. pylori, N. gonorrhoeae,
    Enterobacter, E. coli, K. pneumoniae, A. baumannii, P. aeruginosa; Gram-positive: S. aureus,
    E. faecium, S. pneumoniae).
  - **A real bug caught by cross-checking the two**: `scipy.cluster.hierarchy.dendrogram` with
    `orientation="left"` draws leaf 0 at the BOTTOM and leaf n−1 at the TOP — the opposite of the
    corresponding heatmap's top-to-bottom row order. Both figures now call `ax.invert_yaxis()` so a
    reader can line up the dendrogram against its heatmap by eye; without it the two read in
    opposite vertical directions with no visual cue that anything was flipped.

### Merged predictors, dropped physchem, simplified layout — 2026-09-02

**Property predictors are two MERGED rank-sum columns, not 12 raw ones (user-directed).** Cytotoxicity's
6 columns and abx-resemblance's 3 are each percentile-ranked (`rank_pct`, the same scaling
`organism_scores` already uses) then SUMMED per compound — one column per family
(`eval_auroc_matrix.merged_predictor_scores`, written to `10_merged_predictors.csv` /
`xx_merged_predictors.csv`). **Physchem (`mw`, `clogp`, `tpsa`) was dropped, not merged** — those
three measure unrelated quantities from each other and from "how cytotoxic/antibiotic-like", so
summing their ranks would not be a meaningful number the way a within-family sum is. The matrix goes
from **15 x 27 to 15 x 17** (15 organisms + cytotox + abx). This is a **shared-constant change**
(`default.AUROC_MATRIX_BLOCKS`), confirmed with the user to apply to both this analysis and the
non-abx robustness check — see the "Non-abx robustness check" section below for what changed there.

**Sum, not mean — a display-magnitude choice, not a scientific one.** Both are a monotone rescale of
each other (`mean = sum / n_columns` for a fixed `n_columns`), so the ranking, AUROC and top-N
selection built from the merged column are byte-identical either way; nothing downstream depends on
which was used. `default.PREDICTOR_MERGE_AGG = "sum"` records the choice for review.

**Bonus finding, not requested but a direct consequence:** the non-abx section's old
`NON_ABX_MATRIX_BLOCKS` workaround — cutting the abx block down to `abx_score` alone because its
filter forces `num_sim_0_5_all`/`_subset` to a constant 0, which used to make `classify_predictor`
call them binary and `aggregated_matrix` raise — is no longer needed. `DataFrame.rank(pct=True)`
gives every tied (constant-0) row the same 0.5, a pure constant folded into every compound's sum, so
the merged abx column stays well-defined and continuous. **Verified on the actual re-run**: on the
filtered subset, `abx__merged__rank_sum` ranges exactly `[1.000, 2.000]` — two constant columns at
0.5 each plus `abx_score`'s own real `[0, 1]` rank_pct, exactly as the arithmetic predicts. The two
scripts now use the identical `AUROC_MATRIX_BLOCKS`.

**Missing-value handling in the merge**: `merged_predictor_scores` sums with `skipna=False`, so a row
missing any of a family's raw columns gets NaN in the merged column rather than a silently-partial
sum over fewer inputs than the rest of the library (66 rows have a missing `eos6ojg` value in the
full abx block; only 6 of those fall inside the specific 3 columns merged here). `aggregated_matrix`
and `_top_indices` already handle a NaN-bearing property column exactly this way (pairwise-complete
AUROC, NaN excluded from ever entering a top-N set) — no special-casing was needed.

**Column labels, top tracks, and panel shape (further user-directed cosmetic pass, same day):**
- **`cytotox_rank_sum`/`abx_rank_sum` were tried and rejected** in favour of plain **`cytotox`** /
  **`abx sim`** (`plots_auroc_matrix.PREDICTOR_FAMILY_LABELS`) — both merged columns are literally
  named `"rank_sum"` (see `AUROC_MATRIX_BLOCKS`), so *some* disambiguation was necessary; the family
  name alone reads cleaner than repeating the column name. "abx sim" rather than bare "abx": the
  merged column folds `abx_score` in alongside the two similarity counts, so "sim" (matching the
  raw `num_sim_0_5_*` column names) reads more accurately than the bare family key would.
- **The top model/class/block tracks are GONE** from both `AurocMatrixPlot` and `OverlapMatrixPlot`
  (user-directed: "not so informative now"). With only 2 property columns left, "model" was always
  "merged", "class" over the property columns was always the neutral/NA colour, and "block" only ever
  separated 15 bioactivity columns from 2 property ones — a distinction the `PREDICTOR_FAMILY_LABELS`
  column labels already make directly. `BLOCK_HUES` and the associated `col_block`/`col_class`
  (top)/`col_model`/`model_ids`/`model_colors` variables were deleted, not just left unused. Only the
  row-side (left) organism-class band remains.
- **Cells are now literally square**, not stretched rectangles: `heatmap()` gained an `aspect`
  parameter (default `"auto"`, unchanged for every other caller — `plots_eos3dys.py` and
  `plots_euopenscreen.py` also share this function) and the two matrix classes pass `aspect="equal"`.
  Panel footprint narrowed from the 27-column-era `cells=(3, 6)` (180 x 90 mm landscape) to
  `cells=(3, 3.6)` (108 x 90 mm), sized for the new 17-column width; `aspect="equal"` makes the cells
  exactly square regardless of any residual mismatch between this footprint and the actual data
  aspect, rather than depending on getting the mm math exactly right.

### Physchem overlap redefined as "top extreme" — 2026-08-31 — **SUPERSEDED 2026-09-02**

**Physchem was dropped from the matrix entirely on 2026-09-02** (see the section above) — this whole
mechanism (`_top_extreme_indices`, `OVERLAP_EXTREME_FAMILIES`, the `Δ`-prefixed column labels) no
longer applies to anything and has been deleted from the code, not just left unused. Kept below as
the historical record of why it existed and what it found while physchem was still in the matrix.

**For the overlap matrix's physchem block only** (`mw`, `clogp`, `tpsa` — the 3 columns in
`AUROC_MATRIX_BLOCKS`'s `"physchem"` family), a property's own "top N" is now the N compounds
FURTHEST FROM THE COLUMN MEDIAN in either direction (`eval_auroc_matrix._top_extreme_indices`,
gated by `default.OVERLAP_EXTREME_FAMILIES = {"physchem"}`), not the N compounds with the highest raw
value (`_top_indices`, still used everywhere else, including the abx-resemblance and cytotoxicity
property blocks). User-directed: MW/logP/TPSA can associate with pathogen activity by being unusually
far from typical in EITHER direction — an unusually small OR large molecule — not only by being
heavy, lipophilic or polar, so "the top 1000 MW compounds" was a one-sided question where a two-sided
one is more informative. **This changes only the overlap matrices' physchem columns**
(`10_overlap_matrix_top{100,1000,10000}*.csv`); the AUROC matrix's physchem columns are unaffected —
they were never a top-N selection to begin with, only a continuous Mann-Whitney rank against each
organism's own top-N actives.

**A caveat worth carrying alongside this change:** the AUROC block still assumes physchem associates
*monotonically* with activity (`rankdata` over the raw column, direction-signed but one-directional).
If a property's real association is genuinely two-sided — extreme in *either* direction predicts
activity — the AUROC matrix would still read it as chance (~0.5) even though the overlap matrix now
surfaces it. The two views are not measuring the same thing for physchem any more, and that should be
stated wherever both are read together. No AUROC-side change has been made; flagged for a decision if
a folded/extreme AUROC is wanted too.

**Verified against the current library** (`10_physchem_matrix_named.csv`, medians:
mw 374.09, clogp 3.32, tpsa 74.68): the switch is a near no-op for two of the three columns and a
complete flip for the third.
- **`clogp`: 0 of 1000 compounds in common** with the old top-1000-highest set — every one of the
  1000 most-extreme compounds is now BELOW the median, because the library's low-logP tail is
  further from the median (median − min = 5.32) than its high-logP tail (max − median = 3.68). "Top
  clogp" used to mean "most lipophilic"; the overlap matrix's `clogp` column now means "most
  hydrophilic", the opposite claim.
- **`mw` and `tpsa`: 1000/1000 and 999/1000 unchanged** — for both, the upper tail already happens to
  be the more extreme one in this library, so "top extreme" and "top highest" pick almost the same
  compounds. The redefinition is not vacuous for these two, just currently non-binding; it would bite
  if the reference library's distribution shifted.
- The top-10000 overlap matrix's off-diagonal minimum moved 7 → 0 between the two runs, a direct,
  verified consequence of the `clogp` flip (no other column changed).

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
- **`14_excluded_targets.csv` is correctly absent** — it is written only when a selected endpoint is
  missing from the parquet, which was true when the cache lacked `eos9ivc`. The file had survived
  every rename since as a stale artifact. Step 14 writes **5** CSVs, not 6.
- **1,200 of the 30,300 pairs are undefined** — the 4 constant-zero abx endpoints x 300 targets. A
  constant column has no ranking, so no AUROC exists. Recorded, never imputed.

**Not regenerated, by decision:** the five scaled matrix CSVs. Nothing reads them (only 3 of the 5
would ever be consumed, all by step 09, which derives them in memory instead) and they cost ~1 h 54 m.

### Partial rebuild at 307 endpoints — 2026-08-11

Triggered by adding `eos3f8h`'s 7 endpoints to the selection (300 -> 307). **Only the steps that read
the endpoint selection or step 07's parquet were re-run** — 07, 08, 09, 10, 14, all exit 0. Steps 11,
12 and 13 read neither (step 12 has its own separate parquet cache), so they were left alone and step
10 consumed their existing outputs. The prior outputs were moved aside to force the rebuild and then
**deleted with `tmp/` on 2026-08-12** (user-directed), so the 300-endpoint values survive only as the
figures quoted in this file.

| step | verification |
|---|---|
| 07 | parquet 1,355,109 x **402** (all Yes+No columns), **0 models skipped**; named matrix 1,355,109 x **307**; mean rank 0.5000 / **SD 0.1305** |
| 08 | 3 Jaccard matrices at **307x307**; kept **260 of 307 columns, 15 of 58 pathogens**; both identity checks unchanged |
| 09 | **byte-identical to the 300-endpoint build** — it draws only `eos3dys` CoAdd endpoints (20 across 9 organisms), so `eos3f8h` never enters it |
| 13 | 101 x 307 = **31,007** pairs (**1,228 undefined**); activity-self **307x307 = 94,249**; pathogen subset **70x70 = 4,900** |
| 14 | **15 x 27**, 405 cells none missing, **diagonal exactly 1.0**, organism-score column means 0.5000 |

**Six organism aggregates each gained exactly one endpoint**, and none gained more: *A. baumannii*
4->5, *C. albicans* 3->4, *E. coli* 11->12, *K. pneumoniae* 3->4, *P. aeruginosa* 5->6, *S. aureus*
4->5. *P. falciparum* stays at 13 and the five single-endpoint organisms (Campylobacter, Enterobacter,
*E. faecium*, *H. pylori*, *S. pneumoniae*) are unchanged, so step 10 is still **15 rows**.

**The AUROC matrix moved much more than "+7 of 307 endpoints" suggests: 210 of 405 cells changed**,
max absolute delta **0.115**, mean **0.0098**. This is mechanical, not a new finding — a changed
organism aggregate re-binarizes that organism's top-1000 actives, which shifts every cell in its row
*and* every cell in its column, so six affected organisms touch ~210 cells. **Any number quoted from
the 300-endpoint AUROC matrix must be re-read from the new CSV**, not carried over.

**Weight warning for the small aggregates.** For organisms with few endpoints, one `eos3f8h` column is
a large share of the merged score — a third of *C. albicans* and a quarter of *K. pneumoniae*. Those
rows now lean substantially on a single model, and on one whose endpoints are expected to correlate
with the others present. Whether that is acceptable weighting is a judgement call, not settled here.
(The CSV-export flag this note originally pointed to was removed 2026-09-01 as unused — see the
"Pipeline reorganization" note above; the matrices themselves are unaffected, only that export path.)

**Still stale:** step 01 only, for an unrelated reason — the metadata refresh (Ready 208 -> 214).

## 11_reference_library_projection.py
Visualises where the reference library sits in chemical space, highlighting each pathogen's
highest-predicted-activity compounds. `eos1klk` (fetched by `00_download_data.py`) computes four 2D
projections of the library — PCA, UMAP, t-SNE, TMAP; each gets its own figure: a silver background of the
full library's density, with each of the 15 pathogens in `config/pathogens_of_interest.csv` getting its own
panel showing its `PROJECTION_TOP_N` highest-scoring compounds (by `consensus_score`) in crimson. This is a
**rank cutoff (top-N), never a score threshold** — no `consensus_score` cutoff value is chosen or reviewed,
only a compound count. Analysis in `src/eval_projection.py`, figures in `src/plots_projection.py`; writes to
`output/11_reference_library_projection/`.
**The UMAP panel is redrawn to match step 10's AUROC matrix exactly, so it now differs from the PCA/t-SNE/
TMAP panels described above** (`Matrix14OrderedUmapGridPlot` in `src/plots_projection.py`,
`matrix10_top_n_per_organism`/`run_matrix10_umap` in `src/eval_projection.py`):
- **Row order** matches step 10's phylogeny-within-organism-class order exactly (read from
  `output/10_auroc_matrix/10_row_order_comparison.csv`'s `phylo_position`), not the alphabetical order the
  other three methods still use.
- **Top-1000 selection** is step 10's own aggregate per-organism score (rank-pct-scaled endpoints averaged
  per organism, read straight from `output/10_auroc_matrix/10_organism_scores.parquet`), not this script's
  own per-pathogen `consensus_score`. **These are different compound sets** for any organism with 2+ merged
  endpoints — only the five single-endpoint organisms (Campylobacter, Enterobacter, E. faecium, H. pylori,
  S. pneumoniae) are guaranteed identical to the PCA/t-SNE/TMAP top-1000.
- **Optional, not fatal, if `10_auroc_matrix.py` has not been run yet** (changed 2026-09-01 — see the
  "Pipeline reorganization" note after step 07 for why). Every other output of this script — the
  background grids, the per-pathogen PCA/t-SNE/TMAP panels, the coadd panels — has no dependency on
  step 10 and is produced regardless. Only this one panel needs it; if step 10 hasn't run, the UMAP
  panel falls back to this script's own per-pathogen `consensus_score` ranking (the same one PCA/
  t-SNE/TMAP use) with a printed note, and rerunning this script after step 10 upgrades it to the
  matched panel. Before this fix, the check ran at the top of the script and blocked EVERY output on
  step 10 — which was one leg of a real circular dependency (this script's background also gates
  steps 12/13/14's figures, and those steps' matrices used to gate step 10 in turn).
- **Density-shaded scatter, one shared colour scale** — user-directed, revised 2026-09-01 (the original
  4-hue `spectrum_cmap(OVERLAP_MATRIX_SPECTRUM)` read as busy/uneven at this scale): each panel's top-1000
  points are shaded by local KDE density (`scipy.stats.gaussian_kde`) through
  `stylia.DivergingColormap("crimson_cobalt").cmap.reversed()` — the house "coolwarm" equivalent (also the
  hue pair matrix 14's own AUROC matrix uses), cool cobalt = sparse, warm crimson = dense, one shared scale
  across all 15 panels. Points are drawn sparsest-first so the densest cluster in each panel is always on
  top, never partly painted over by sparser points. Marker size is smaller than the PCA/t-SNE/TMAP panels'
  (`MATRIX10_UMAP_POINT_SIZE = 2.5` vs `POINT_SIZE = 4`), sized for the smaller panel below.
- **Footprint: `cells=(1.8, 3.0)` = 90 x 54 mm** — user-directed, "half of stylia's two-column print width"
  (180 mm / 2 = 90 mm) for a 3-row x 5-col grid of strictly square 18 mm panels (90 mm / 5 cols = 18 mm;
  3 rows x 18 mm = 54 mm). A fractional-cells figure per the sanctioned exception in
  `docs/figure_conventions.md`.
- **No frame, no ticks, no "UMAP 1"/"UMAP 2" axis labels, no marker legend, no class-colour swatch**
  (the last two dropped 2026-09-01, user-directed) — the only panel chrome is the pathogen name, drawn as a
  real `loc="left"` axis title (the same mechanism `plotting_utils.roc_panel` uses for its own per-cell
  identifier) so it sits in the space above the axes and can never overlap a scatter point, however the
  data happens to be laid out.
- **Fixed axis range `(-1, 1)` both axes** — user-directed. The real UMAP extent runs slightly wider on both
  axes (roughly `-1.17` to `1.09`, per `11_umap_background.csv`), so a thin sliver of the outer background —
  and possibly a handful of top-1000 points — falls outside this fixed frame and is not drawn.
**Memory:** only one pathogen's `key` + `consensus_score` columns (of each ~424 MB prediction file) are read
at a time, immediately reduced to its `PROJECTION_TOP_N` highest rows and discarded — no more than one
pathogen's raw scores, and never all 15 prediction files, are held in memory together. The score ranking
does not depend on projection method, so each pathogen's top-N is computed once and reused across all four.
**Second figure family — the same view from `eos3dys` (CoAdd), UMAP only.** Alongside the 15 ChEMBL-derived
pathogen models, the script draws one more grid driven by a single independent predictor: `COADD_MODEL_ID`
(`eos3dys`, `coadd-antimicrobial-activity`, an array of LazyQSAR models trained on CoADD screening data),
which covers **9** of the same organisms. Same library, same layout, same background grid (reused, not
recomputed) and the same top-1000 rank cutoff, so the two families can be read side by side: agreement is
evidence from an independently trained model, disagreement is a finding.
- **`eos3dys` has no `consensus_score`** — its 22 outputs are independent per-strain, per-assay endpoints —
  so its top 1000 is ranked **per endpoint**, on that endpoint's own score, not per organism. Still a pure
  rank cutoff.
- **Endpoints are grouped onto organisms, one panel each** (9 panels, 3x3 grid). Seven organisms carry 1–2
  endpoints; **E. coli and P. aeruginosa carry 4 each** (extra strains: `lpxC`/`tolC`, `PAO397`). Endpoints
  sharing a panel are drawn as shades of one hue, with the **darkest shade the wild-type reference strain**
  — that ordering comes from `config/08_endpoint_selection.csv` row order, which happens to list wild-type
  before sensitised/resistant. Each endpoint's name is printed in its own shade, so the annotation is the
  key (a shared legend cannot name endpoints that differ per panel).
- **Endpoint set** = the `eos3dys` rows of `config/08_endpoint_selection.csv` with
  `assay_type == "bioactivity"` → **20** endpoints. This excludes its two `Homo sapiens` columns
  (`cytotoxicity_ic50`, `hemolitic_activity`), which are selected in `config/cytotoxicity_models.csv` and
  belong to step 13 instead.
- **UMAP only** (`COADD_PROJECTION_METHOD`) — 9 organism panels at four methods each is 4x the panels for
  no added insight, and UMAP is the layout the pathogen figures are read from. Same rationale as
  `TOX_PROJECTION_METHOD`; the engine takes the method as an argument, so this is revisable.
**Key decisions, for review:**
- **Top N per pathogen** (`PROJECTION_TOP_N`, in `src/default.py`) is **1000** — user-directed. The same
  count is used per endpoint for the `eos3dys` family.
- **Grid resolution** (`PROJECTION_BINS`) for the silver background density is **60x60** cells per method —
  sized to stay legible at the ~30x30 mm panel size in a 15-pathogen small-multiples grid, not a fitted
  value.
**Known metadata gap (not fixed here):** `eos3dys`'s Airtable `Target Organism` field omits *Escherichia
coli* although the model has four E. coli endpoints. This figure reads organisms from
`config/08_endpoint_selection.csv`, so it is unaffected, but the Airtable record should be corrected.

### Former `10_physchem_matrix.py` — folded away 2026-09-01, its stats/figures removed 2026-09-02

This script no longer exists as its own file/folder, and its history has two stages. **First
(2026-09-01):** its matrix-building half (raw physchem calculation, column naming, the "all 22 kept
though 5 are near-redundant" and "upstream spelling kept" notes) moved to `08_property_matrices.py`
above, and its stats/figures half (endpoint stats CSV, the 22-panel distributions grid, the
3-descriptor UMAP panel) moved to `14_confounder_analysis.py`, next to the predictor-performance
analysis that was physchem's only remaining consumer. **Second (2026-09-02, user-directed):** that
stats/figures half was removed from step 14 entirely, not moved again — it tested nothing about
activity, and step 14's whole point is the confound check, so a purely descriptive block of physchem
figures didn't belong there. Physchem itself is unaffected: its matrix (step 08) still feeds step
14's predictor-performance analysis as one of the three property families, unchanged; only its
standalone per-descriptor stats/distributions/UMAP are gone. **All prior findings quoted under this
heading are historical** (from before the 2026-09-02 removal) and recoverable from `git log` if the
per-descriptor stats/distributions are ever needed again: 5 near-redundant descriptors (verified on
200k rows: `n_rings`/`n_aliphatic_rings`/`n_aromatic_rings`/`n_saturated_rings` each exactly the sum
of their two carbocycle/heterocycle components; `n_radical_electrons` 0 for all but ~38 of 200,000),
observed ranges MW 198.9–699.5 / cLogP −2.0–7.0 / `n_rings` ≤ 6 (the library's filter, not the
model's), and the three misspelled Datamol column names kept verbatim to match the real CSV header.

## 12_abx_projection.py
The antibiotic-resemblance score matrix's endpoints on the library UMAP, plus their pathogen
enrichment. **Matrix-building moved to `08_property_matrices.py` on 2026-09-01** (see the "Pipeline
reorganization" note after step 07) — this step now reads `08_abx_matrix_named.csv` directly rather
than rebuilding it from `config/antibiotic_resemblance.csv`'s 55 selected endpoints across 4 models,
so it no longer needs a raw prediction file open. Everything below is otherwise unchanged: same
column naming (`abx__{model_id}__{column_name}`, a constant group code in the pathogen slot), same
un-normalized convention. Engine `src/eval_abx_matrix.py`, figure `src/plots_abx_projection.py`;
writes to `output/12_abx_projection/`.
**The highlight rule is not a rank cutoff.** 44 of the 55 endpoints are binary or small integer counts,
so each panel shows every compound with a value > 0, capped at `PROJECTION_TOP_N` — never padded. **19
endpoints hit that cap**, so their panels show an arbitrary key-ordered subset, not the full flagged
set; read each panel's `n_shown/n_nonzero` annotation.
**4 endpoints are constant zero library-wide** (`arsenic_cpds`, `b_lactamase_inhibitors`,
`lipopeptides`, `polypeptides`) and are omitted from the **figure only** — they remain in the matrix and
the stats CSV. 32 of 55 have fewer than 1000 non-zero compounds. 66 missing values, kept.
**Reuses step 11's `11_umap_background.csv`** rather than recomputing the density, so its panels are
directly comparable to the pathogen ones.

### Pathogen top-1000 x eos19mt antibiotic-class Fisher enrichment (replaces the former Jaccard/UMAP
overlap analysis, 2026-09-01)
For every pathogen and every one of eos19mt's 38 antibiotic structural classes (not the mixed
9-endpoint, 4-model set the earlier overlap version used — "belonging to specific antibiotic classes
as defined by eos19mt" is the question actually asked), tests whether being in the pathogen's top
`PROJECTION_TOP_N` (by `consensus_score`, from step 11) is associated with belonging to the class
(value > 0), via a one-sided (`alternative="greater"`) Fisher exact test over the full ~1,355,109
compound reference library. Analysis in `src/eval_abx_enrichment.py`, figure in
`src/plots_abx_enrichment.py`; both replace and delete `src/eval_pathogen_abx_overlap.py` /
`src/plots_pathogen_abx_overlap.py` and their outputs (user-directed).

**Contingency table**, per (pathogen, class) pair, over the full library:

```
                      in top-N     not in top-N
    class member         a              c
    not class member     b              d
```

Compounds with a missing value in that class's column are **excluded from the pair's table**
(never imputed to 0) — the class's own missingness. `scipy.stats.fisher_exact` on `[[a, b], [c, d]]`
gives `(odds_ratio, p_value)`. A class with **no variance** over the resulting universe (e.g. the
four constant-zero eos19mt columns — `arsenic_cpds`, `b_lactamase_inhibitors`, `lipopeptides`,
`polypeptides`) has nothing to test: its `odds_ratio`/`p_value` are `NaN`, not a
computed-but-meaningless value.

**Row order** is eos19mt's 38 classes in `config/antibiotic_resemblance.csv`'s file order — verified
to match the model's own native output order (`eos19mt_v2.csv`'s header). **Column order** is the 15
pathogens in step 10's phylogenetic order, reusing `eval_auroc_matrix.bioactivity_order` /
`organism_order` / `phylogeny_organism_order` rather than re-deriving it, so this matrix's columns can
never disagree with the AUROC matrix's rows.

**`12_abx_enrichment_long.csv`** — 570 rows (15 pathogens x 38 classes): `pathogen_code`, `pathogen`,
`column_name`, contingency counts `a,b,c,d`, `n_class_nan`, `odds_ratio`, `p_value`, and a
Benjamini-Hochberg-adjusted `p_value_fdr` (for reference only — the wide matrices below and the
figure's asterisks always use the raw, unadjusted `p_value`, since multiple-testing correction was
not part of what was asked).

**`12_abx_enrichment_odds_ratio.csv` / `12_abx_enrichment_pvalue.csv`** — the two 38 x 15 wide
matrices (rows = classes, columns = pathogens, both orders above), raw values, one metric each.

**Figure** (`12_abx_enrichment_matrix`): colour = `log2(odds ratio)`, diverging around 0, clipped at
`ABX_ENRICHMENT_LOG2OR_CAP` so the 0/infinite odds ratios that occur (e.g. `glycopeptides` has a
single non-zero compound library-wide, so it lands either entirely inside or entirely outside a
pathogen's top-N) still render at the scale's extreme rather than breaking the colormap — the printed
cell text is always the TRUE odds ratio, never the clipped colour value. A trailing `*`/`**`/`***`
marks the raw p-value at `ABX_ENRICHMENT_SIG_THRESHOLDS` (0.05/0.01/0.001) — a display convenience,
not a multiple-testing claim. One top track shows each pathogen's `organism_class`.

## 13_toxicity_projection.py
Endpoint stats and the toxicity projection for the **cytotox** block (`config/cytotoxicity_models.csv`,
24 selected endpoints across 4 models: `eos42ez`, `eos7m30`, `eos3le9`, `eos3dys`).
**Matrix-building moved to `08_property_matrices.py` on 2026-09-01** (see the "Pipeline
reorganization" note after step 07) — this step reads `08_cytotox_matrix_named.csv` directly rather
than rebuilding it, so it no longer opens a raw prediction file for the matrix itself. Writes
`output/13_toxicity_projection/`. Endpoint-stats logic shared with step 08 via `src/eval_property_matrix.py`.
**Every column here is a prediction**, unlike the physchem block's deterministic descriptors. The two
must not be treated interchangeably.
**`eos7m30`'s own 8 physicochemical columns stay `No`** in the config, so `eos4djh` remains the single
physchem source and no two blocks can disagree about molecular weight or logP.
**Carrying the model ID matters here:** two models score HepG2 —
`cytotox__eos42ez__cytotoxicity_hepg2` vs `cytotox__eos3le9__ic50_hepg2_72h_5um` — and the column name
is what keeps them apart.
**Endpoint selection is not made here** — it is the manually curated `selected` column. No threshold or
cutoff is applied; raw outputs pass through unchanged. 0 missing values across all 24 endpoints.

### Toxicity projection (merged in from the former step 13, 2026-08-06)
The toxicity counterpart of step 11, on the same `eos1klk` chemical-space layout: a silver
full-library density background with, in crimson, each toxicity endpoint's most toxic compounds. One
small-multiples figure (`13_umap_top1000_toxicity`), 24 panels — one per `selected == "Yes"` endpoint in
`config/cytotoxicity_models.csv` — scored from step 08's cytotox matrix. Analysis in
`src/eval_tox_projection.py`, figure in `src/plots_tox_projection.py`; writes to
`output/13_toxicity_projection/`. Panels are labelled with both the endpoint and its model ID, since
several endpoints read overlapping biology from different models (two HepG2 readouts, three cytotoxicity
readouts).
**Top N per endpoint** is `PROJECTION_TOP_N` = **1000**, shared with step 11 — a **rank cutoff, never a
score threshold**, so no score value is chosen or reviewed, only a compound count.
**Projection method:** only UMAP (`TOX_PROJECTION_METHOD`) of the four `PROJECTION_METHODS` is drawn —
24 endpoints x 4 methods would be 96 panels, and UMAP is the layout step 11's pathogen figures are read
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
**Memory:** the ~720 MB step-13 cytotox matrix is streamed in chunks reading `key` + the 24 score columns only
(never the `input` SMILES column, which is most of the file), with each endpoint's running top-N reduced
after every chunk.

## 14_confounder_analysis.py

**Physchem's own endpoint stats/distributions/UMAP panel — gained 2026-09-01, removed 2026-09-02.**
This step briefly carried physchem's descriptive figures (moved here when the former standalone
`10_physchem_matrix.py` was folded away), then had them removed again — user-directed, on review —
since they tested nothing about activity and this step's whole point is the confound check below; a
purely descriptive block of physchem figures didn't belong in it. See the "Former
`10_physchem_matrix.py`" note (in step 11's section) for the two-stage history and the historical
findings. Physchem itself is unaffected — its raw matrix (step 08) still feeds the predictor
performance analysis below as one of the three property families — only its standalone figures are
gone, and this step no longer depends on step 11 (the projection background) at all.

The script asks whether the property/resemblance columns carry any signal about pathogen activity.
Treats every column of step 08's physchem, abx and cytotox blocks as a
**predictor**, and every curated activity endpoint as a binary **target**, giving one performance
value per (predictor, target) pair — **101 x 307 = 31,007** as re-run on 2026-08-11 (was
101 x 300 = 30,300 on 2026-08-06 and 101 x 260 = 26,260 at 260 endpoints). Predictors by family: abx
55, cytotox 24, physchem 22 — all three read from step 08's cache directly, so this step needs
nothing but steps 07, 08 and (for the confounder check below) 09 — no figure-drawing script (11, 12,
13) has to run first. Three figures, one per predictor family, each a
box-with-jitter per predictor over its distribution across all 260 targets, sorted by median, with a
chance line at 0.5. Analysis in `src/eval_predictor_performance.py`, figures in
`src/plots_predictor_performance.py`; writes to `output/14_confounder_analysis/`.
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
ranked once (ties averaged, so the result is exact) and each target is then a 1000-element gather. 31,007
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
`14_predictor_performance.csv` (exactly 40 rows) for exclusion downstream.
**Note (2026-09-02): this step used to draw two more figures here** — activity endpoints as
predictors of each other, and the same restricted to the 15 pathogens of interest — using data that
turned out to be bioactivity-only (no property data entered either computation). Both moved to step 09
in the second reorg pass; see that step's "Part 3" section for the current `09_activity_self_performance.csv`
/ `09_pathogen_subset_self_performance.csv` and their figures.

**The curated 12 predictors, three families in one panel**
(`14_performance_curated_predictors`, from `14_curated_predictor_performance.csv`). A user-directed
shortlist of the 101 property columns (`CURATED_PREDICTORS` in `src/default.py`) scored against the same
59 consensus-collapsed pathogen endpoints step 09's `09_pathogen_subset_self_performance.csv` uses
(`pathogen_subset_endpoints`) — 12 x 59 = 708 pairs. One box per
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

### Endpoint confounder check (absorbed from the former standalone `15_bioactivity_endpoints.py`, 2026-09-02)

For every endpoint in step 09's `09_endpoint_quality.csv` ranking (255 endpoints, the same >5-endpoint
scope Part 3 of that step uses), finds the strongest non-same-model physchem/abx/cytotox predictor of
its own top-1000, using `eval_endpoint_quality.confounder_stats` against this step's own
`14_predictor_performance.csv`. The question is different from step 09's peer-agreement one: an
endpoint whose top-1000 is recovered almost perfectly by a single property column is, whatever its
peer agreement, largely reproducing a property filter rather than a bioactivity ranking.

**Strength is `|value - chance|`**, so a strong INVERSE predictor counts as strongly as a direct one;
the signed value is kept in `confounder_value` so the direction stays readable. `same_model` rows are
excluded — `eos3dys` appears in both `config/cytotoxicity_models.csv` and the activity endpoint
selection, so leaving them in would let one of its columns "predict" another of its own columns and
be reported as an external confounder.

**Written as its own output for now, not merged into step 09's table** ("for the moment", per the
review that decided this split): `14_endpoint_confounders.csv` (`endpoint`, `confounder_predictor`,
`confounder_family`, `confounder_metric`, `confounder_value`, `confounder_abs_dev`), joined against
step 09's ranking only for the printed report below, at read time rather than at write time — step 09
runs before this step, so it cannot carry this step's confounder columns when its own table is built.

**All 255 in-scope endpoints have a non-same-model confounder candidate** — physchem is the
strongest predictor for 111 of them, cytotox for 89, abx-resemblance for 55.

**Among the 40 weakest endpoints by step 09's `auroc_out_same_median` ranking, 5 have a confounder at
`|AUROC - 0.5| >= 0.4`:**

| endpoint | confounder predictor | confounder value |
|---|---|---|
| `eos8jx6:chembl_dose_response_7` | `physchem__eos4djh__n_hetero_atoms` | 0.908 |
| `eos4an7:chembl_single_point_5` | `physchem__eos4djh__sas` | 0.946 |
| `eos8lcw:chembl_single_point_2` | `abx__eos11sm__abx_score` | 0.928 |
| `eos21dr:chembl_dose_response_7` | `physchem__eos4djh__tpsa` | 0.056 |
| `eos5eya:chembl_single_point_5` | `physchem__eos4djh__n_hetero_atoms` | 0.073 |

e.g. `eos8lcw:chembl_single_point_2`'s (*S. aureus*) top-1000 is recovered at AUROC 0.928 by
`abx_score` alone, so its "actives" track antibiotic-resemblance at least as much as a genuine
*S. aureus*-specific signal. The two lowest-value rows above are equally strong confounders read the
other way: a LOW value means the predictor's own top-1000 (highest TPSA / most heteroatoms) is
enriched for the endpoint's *inactives*, which is exactly as informative as a high value once the
sign is read correctly — this is why strength is measured as absolute deviation from chance rather
than the raw value.

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

### Non-abx robustness check (absorbed from the former standalone `xx_non_abx_matrix.py`, 2026-09-02)

Repeats the analysis above — the exact same two figures, same `phylo` row order, same
`continuous_color` overlap scale — on a library purged of antibiotic-like compounds. If the
organism-to-organism structure above were driven by known antibiotics converging on the same
scaffolds, the matrix would be reporting the library's composition rather than the models' behaviour.

**Two behaviour changes came with the fold, both required by "draw the exact same plots as the main
analysis":** the 100/10000 overlap cutoffs the standalone script also drew are gone (the main
analysis dropped those on 2026-09-01, this section now matches it — only top-1000 is drawn), and the
overlap figure now uses the same `continuous_color=True` linear scale the main analysis switched to
on 2026-09-02, not the discrete non-uniform-bin scale the standalone script still used. Neither is a
new decision — both simplifications already happened to the *main* analysis; this fold just brings
the non-abx section into line with it. The `phylo` row order itself is not recomputed — it depends
only on the fixed 15-organism list and `config/organism_taxonomy.csv`, neither of which changes
under the mask — and the dendrogram diagnostic is not redrawn (it would be identical).

**Filter** (user-directed; a compound is kept only if all hold), from step 08's abx parquet:
`eos11sm:abx_score < 0.5`; `eos2xeq:` `is_sim_known_ab`, `nitrofuran_motif`, `fluoroquinolone_motif`,
`carbepenem_motif`, `betalactam_motif` all `== 0`; `eos6ojg:` `num_sim_0_5_all`, `num_sim_0_5_subset`
both `== 0`. (`carbepenem_motif` is spelled with an "e" upstream — verified, not a typo. The unrelated
eos19mt column `b_lactams_carbapenems` uses an "a".) Per-filter pass counts are written to
`10_nonabx_filter_summary.csv`. Kept exactly as-is for now — not revisited in this pass.

**Cutoffs reviewed and kept on 2026-08-11.** The filter removes only 2.38%, so the alternatives were
quantified before deciding to keep it. The levers and what each would have cost:

- **`abx_score < 0.5` sits at the 99th percentile** (median 0.123, q90 0.219, q99 0.513), so it
  removes just 1.07% and is close to a no-op. Tightening to `< 0.3` / `< 0.25` / `< 0.2` would have
  removed 4.04% / 6.80% / 13.23% on that term alone.
- **Tanimoto 0.5 is a tight similarity radius** — only 1.45% of the library has any known-antibiotic
  neighbour there. `num_sim_0_3` flags **28.81%** and is by far the strongest available lever.
- **`eos19mt`'s 38 antibiotic-class structural flags are in `config/antibiotic_resemblance.csv` but
  deliberately not used** — they would flag 11.19% of the library (or 5.92% restricted to 22 core
  antibiotic classes, though 76,473 of those 80,234 hits are `sulfonamides` alone, a motif ubiquitous
  outside antibiotics).

The strictest combination considered (motifs + eos19mt + `num_sim_0_3` + `abx_score < 0.15`) would
have retained 51%. **The original cutoffs were kept**, so this analysis is a narrow test of *close
antibiotic analogues* rather than of antibiotic-like chemistry broadly — which is what its
conclusions can support, and no more.

**Retains 1,322,835 of 1,355,109 (97.62%)**; 32,274 removed. **6 compounds have a NaN in an eos6ojg
column and are excluded**, since NaN fails `== 0` — the conservative reading, reported in the log and
the summary CSV rather than silently applied.

**Ranks are recomputed within the subset** (user-directed): the mask is applied *before*
`scale_matrix`, so percentile ranks are relative to the 1.32M survivors. This is the counterfactual
"what if the library had never contained these compounds", not "which of the full library's top-N are
antibiotic-like". It is also what keeps the column-mean ≈ 0.5 assertion meaningful.

**The abx block is cut to one column, giving 15 x 25 — SUPERSEDED 2026-09-02.** Kept as history: the
filter used to force `num_sim_0_5_all`/`num_sim_0_5_subset` to a constant 0 for every retained
compound; `classify_predictor` then called them binary and `aggregated_matrix` raised, because the
matrix puts every column on one AUROC scale. The workaround was `NON_ABX_MATRIX_BLOCKS`, cutting the
abx block down to `abx_score` alone.

**Uses the same merged cytotox/abx rank-sum columns as the main analysis (15 x 17), and needs no
workaround.** `NON_ABX_MATRIX_BLOCKS`, an earlier abx-block-to-one-column workaround, has been
deleted — this section uses the identical `AUROC_MATRIX_BLOCKS` the main analysis does.
`DataFrame.rank(pct=True)` gives every tied (constant-0) row the same 0.5, a pure constant folded
into the sum, so the merged abx column stays well-defined and continuous without needing to drop
anything. **Verified on the actual re-run**: `abx__merged__rank_sum` on the filtered subset ranges
exactly `[1.000, 2.000]` — two constant columns at 0.5 each plus `abx_score`'s own real `[0, 1]`
rank_pct, exactly as the arithmetic predicts. See the "Merged predictors, dropped physchem" section
above for the full method and the sum-vs-mean rationale
(`eval_auroc_matrix.merged_predictor_scores`, `PREDICTOR_MERGE_AGG = "sum"`). Ranks are still
recomputed WITHIN this section's masked subset before merging, same rule as everything else here.

**Row alignment is the failure mode that would not raise.** The main analysis above indexes its
property CSVs *positionally* from positions derived from the step-07 parquet and never reads their
`key`; a filtered score matrix against unfiltered property CSVs would give a full, plausible, wrong
matrix, and the diagonal check would still pass because it only touches the bioactivity block. Three
guards: `_assert_key_alignment` runs before any compute here (the main analysis above does not call
it — masking is this section's own risk to guard against), the abx parquet's key order is checked
against step 07's, and one mask is threaded through every reader with
`eval_auroc_matrix._apply_row_mask` re-checking its length at each.

**Results (375-cell, 15 x 27/25 era) — SUPERSEDED 2026-09-02**, kept as history since the column set
that produced them (raw physchem `clogp`, per-model cytotoxicity columns) no longer exists:

| block | median &#124;delta&#124; vs step 10 | max |
|---|---|---|
| all 375 cells | 0.0290 | 0.5737 |
| bioactivity 15 x 15 | 0.0253 | 0.3194 |
| property 15 x 10 | 0.0507 | 0.5737 |

39 of 375 cells crossed the 0.5 chance line, the largest against `clogp` (*S. aureus* 0.297 -> 0.685)
and raw cytotoxicity columns. Cross-checked at the time by an independent key-based join on one cell
(returned 0.6848, matching the positional pipeline).

**Re-run on the merged 15 x 17 columns, 2026-09-02:** over 255 cells (17 shared columns), max
`|delta|` **0.5759**, median **0.0250** — a similar median to before, on a fifth as many cells. The
largest shifts are all in the merged **abx** column, not cytotox: *Enterobacter* 0.9345 -> 0.3587
(+0.5759), *S. pneumoniae* 0.7685 -> 0.3774 (+0.3912), *S. aureus* 0.9794 -> 0.6477 (+0.3317),
*A. baumannii* (abx) 0.9963 -> 0.6712 (+0.3251) and (cytotox) 0.5115 -> 0.8392 (+0.3277). This is
expected, not surprising: the abx block is now dominated by `abx_score`'s own truncated-to-(0.004,
0.5) range on this subset (see above), so an organism whose full-library abx AUROC leaned heavily on
`abx_score` specifically will move the most under the filter. **No cell needs re-verification by an
independent join this time** — the mechanism producing the largest shifts (a truncated single column
now carrying nearly all of the merged predictor's signal) is understood directly from the merge
arithmetic, unlike the 2026-08-11 `clogp` shift, which needed the independent check to rule out a
masking bug. The delta table is written to `10_nonabx_auroc_delta_vs_full.csv` and reported, never
asserted against a bound — **this is an observation, not a conclusion**; the interpretation is a
scientific call.
