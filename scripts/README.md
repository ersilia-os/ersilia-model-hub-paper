# Scripts

## 00_download_data.py
Downloads all input data required by the pipeline: the Ersilia reference library, DrugBank SMILES, Airtable model metadata, and precomputed model outputs from Isaura for all Ready annotation models against the reference library. All 129 annotation model CSVs are saved to `data/processed/annotation_preds_ref_library/`. Downloads are skipped if the output file already exists, so the script is safe to re-run.

## 02_chembl_models_performance.py
For each pathogen in `data/raw/chembl_model_reports/`, loads all 5-fold cross-validation report CSVs, computes mean ± std AUROC per model, and generates two plots: a subplot grid with one ROC curve per model (coloured by mean AUROC using FadingColormap cobalt fitted [0.5, 1.0], best model top-left) and horizontal paired rank boxplots (fold 0; actives in turquoise, inactives in crimson). Outputs one summary CSV and two PNGs per pathogen to `output/02_chembl_models_performance/`.

## 02_euopenscreen_auroc.py
Converts EU OpenScreen prediction H5 files to flat CSVs (one per model, saved in `data/processed/02_euopenscreen_preds/`) and computes per-feature AUROCs against experimental activity data from the `eu-openscreen-antimicrobial-tasks` repository. Only the 6–7 pathogens with available ground-truth task files are evaluated. Outputs a dot plot (`auroc.png`) and a scores table (`auroc_scores.csv`) in `output/02_euopenscreen_preds/`.

## 03_coadd_benchmark.py
Evaluates the CoAdd model (eos3dys) against EU OpenScreen experimental data. For each of the 22 CoAdd endpoints (pathogen/strain/condition-specific predictions), computes AUROC against every available EU OpenScreen binary task, excluding compounds that appeared in the CoAdd training set for that endpoint. Outputs an AUROC matrix CSV and a heatmap in `output/03_coadd_benchmark/`.

**Training data:** per-strain binarised MIC files from the `coadd-binary-tasks` sibling repo (`data/processed/coadd/05_binarised_mic/`), staged by `00_download_data.py`. Without this data, leakage removal is skipped and a warning is printed.

## 04_crossactivity_analysis.py
Assesses whether ChEMBL and CoAdd models have pathogen-specific discriminatory power, or whether they learn general antimicrobial features. Three analyses: (1) Jaccard overlap between active compound sets across the 7 EU OpenScreen tasks; (2) specificity index per ChEMBL model (same-pathogen AUROC minus mean cross-pathogen AUROC); (3) pan-active vs specific-active AUROC split — actives shared across ≥2 tasks vs actives unique to a single pathogen — using the matched ChEMBL `consensus_score` and the best-matching CoAdd endpoint. Outputs to `output/04_crossactivity_analysis/`.

**Threshold:** compounds with fewer than 5 specific-actives in a given task are flagged and their specific-active AUROC is set to NaN (insufficient data for a reliable estimate).

## 01_admet_properties.py
Downloads precomputed outputs for model `eos74km` (Antimicrobial class specificity prediction) from the Isaura public bucket and plots the distribution of each predicted property across the Ersilia reference library.

## filter_compound_families.py
Annotates the Ersilia reference library and DrugBank SMILES with known antimicrobial chemical family membership using SMARTS substructure searches. Each compound receives a boolean column per family; assignment is non-exclusive (a compound may belong to multiple families).

**Families covered (12):** beta-lactams, tetracyclines, fluoroquinolones, sulfonamides, oxazolidinones, nitroimidazoles, rifamycins, phenicols, quinolines, nitrofurans, macrolides, diaminopyrimidines.

**SMARTS notes:**
- Fluoroquinolones require both the 4-oxo-3-carboxylic acid quinolone core AND a fluorine atom.
- Macrolides are detected programmatically (SMARTS cannot constrain ring size): a compound qualifies if it contains a lactone within a ring of size ≥ 12.
- Nitroimidazole and nitrofuran patterns use branch notation `[N+](=O)([O-])` (not chain) to correctly attach the nitro group to the ring carbon.
- Rifamycins are matched via the aminonaphthalenediol core (`Nc1cc(O)c2ccccc2c1O`) shared by all clinically used rifamycins; the complex ansa-bridge is not required for matching.