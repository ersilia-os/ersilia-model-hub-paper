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
Analyses the Airtable model metadata (Ready models only) and renders each summary field as its own
figure, plus a pathogen circle-treemap. Value counts for all eight metadata fields are written as
`*_counts.csv` to `output/01_models_metadata/`.
**Shown vs counted:** 6 of the 8 counted fields are plotted (Task and Subtask merged into one
panel); License and Publication Type are counted but not plotted. Biomedical Area and Target
Organism are capped at the top 10 categories in the figure (full counts remain in the CSVs).

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
boxplots. Outputs one summary CSV and two figures per pathogen to
`output/03_chembl_models_performance/`.

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
