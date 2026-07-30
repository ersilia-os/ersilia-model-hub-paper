
# Project-wide random seed for any stochastic step (sampling, jitter, splits).
RANDOM_SEED = 42

AIRTABLE_SHARE_URL = (
    "https://airtable.com/appR6ZwgLgG8RTdoU/shr7scXQV3UYqnM6Q/tblAfOWRbA7bI1VTB"
)
AIRTABLE_BASE_ID = "appR6ZwgLgG8RTdoU"
AIRTABLE_VIEW_ID = "viwy0inGR1xv2Tpfe"

REFERENCE_LIBRARY_URL = (
    "https://raw.githubusercontent.com/ersilia-os/"
    "ersilia-model-hub-maintained-inputs/main/inputs/reference_library_smiles.csv"
)

DRUGBANK_URL = (
    "https://raw.githubusercontent.com/ersilia-os/"
    "sars-cov-2-chemspace/refs/heads/main/data/drugbank_smiles.csv"
)

COADD_MODEL_ID = "eos3dys"

# Maps each model Subtask to its parent Task, so subtasks can inherit the parent
# task's colour in the metadata figures.
SUBTASK_PARENT = {
    "Activity prediction": "Annotation",
    "Property calculation or prediction": "Annotation",
    "Featurization": "Representation",
    "Projection": "Representation",
    "Similarity search": "Sampling",
    "Generation": "Sampling",
}

# Airtable records each model's container runtime benchmark in five columns, one per input batch
# size. A value of -1 means the benchmark was never run for that model — NOT that it took no time —
# so those rows must be skipped, never treated as a measurement of zero.
RUNTIME_COLUMNS = {
    1: "Computational Performance 1",
    10: "Computational Performance 2",
    100: "Computational Performance 3",
    1_000: "Computational Performance 4",
    10_000: "Computational Performance 5",
}
RUNTIME_NOT_MEASURED = -1
#: Batch size reported in the paper figure. Changing this changes which models the runtime panel can
#: show at all, because coverage collapses as the batch grows: at 100 molecules it is 129/131
#: Annotation, 57/58 Representation and 10/19 Sampling; at 1,000 it is 123/131, 55/58 and **0/19**.
#: 100 is chosen so generative models appear, and that is the whole reason — at this batch size the
#: median CP3/CP1 ratio is 0.86, i.e. running 100 molecules takes no longer than running 1, so for
#: annotation and representation models the number is dominated by container startup rather than
#: per-molecule work. It is still the wall-clock a user waits for; it is NOT a throughput measure.
RUNTIME_BATCH = 100
RUNTIME_COLUMN = RUNTIME_COLUMNS[RUNTIME_BATCH]

# Shorter display labels for crowded metadata figure axes.
SUBTASK_DISPLAY = {
    "Property calculation or prediction": "Property prediction",
}
BIOAREA_DISPLAY = {
    "Antimicrobial resistance": "AMR",
    "Diarrheal diseases": "Diarrhoea",
}

ERSILIA_MODEL_IDS = {
    "abaumannii":"eos21dr",
    "calbicans":"eos8jx6",
    "campylobacter":"eos7iak",
    "ecoli":"eos5eya",
    "efaecium":"eos81zy",
    "enterobacter":"eos9bpi",
    "hpylori":"eos9eyo",
    "kpneumoniae":"eos6wb7",
    "mtuberculosis":"eos43d6",
    "ngonorrhoeae":"eos5qya",
    "paeruginosa":"eos2e3s",
    "pfalciparum":"eos4an7",
    "saureus":"eos8lcw",
    "smansoni":"eos8v1a",
    "spneumoniae":"eos5q52",
}

# ---------------------------------------------------------------------------
# External-validation step (05): EU OpenScreen + CoAdd
# ---------------------------------------------------------------------------
# Headline prediction column shared by every step-04 Ersilia model output CSV.
SCORE_COL = "consensus_score"

# "Shared" organisms = the intersection of our ChEMBL pathogen models with the
# external screening data (organisms that HAVE an EU OpenScreen primary assay and
# a CoAdd reference strain). The remaining ChEMBL-only pathogens are addressed by the
# cross-organism analyses, which read the full config table rather than a fixed list.
SHARED_ORGANISMS = [
    "abaumannii", "calbicans", "ecoli", "efaecium",
    "kpneumoniae", "paeruginosa", "saureus",
]

# Hit classes for an EU OpenScreen compound, from how many of the 7 primary assays it hits:
# inactive (0, in every assay where it has a conclusive result), exclusive (1),
# narrow-spectrum (2 to NARROW_MAX_PATHOGENS) and broad-spectrum (more than that).
# The 3-pathogen narrow/broad boundary is a user-set choice, not a fitted value.
NARROW_MAX_PATHOGENS = 3
HIT_CLASSES = ["inactive", "exclusive", "narrow", "broad"]
# The same compounds collapsed to plain activity: active = a hit in one or more of the 7 assays.
ACTIVITY_CLASSES = ["inactive", "active"]

# Step-05 output layout. Results are filed by leakage status so the two are never confused:
# FULL keeps the compounds the models were trained on, DEDUP removes every compound present in any
# of the 7 models' ChEMBL training sets. Outputs with no leakage dimension (label-only tables, the
# leakage audit) stay at the top level of the step's output dir.
FULL_SUBDIR = "full"
DEDUP_SUBDIR = "deduplicated"

# CoAdd reference strain per organism that has both a ChEMBL model and CoAdd data — the CoAdd
# analogue of the EU OpenScreen "shared" set, and the organisms step 06 evaluates. This is the
# strain whose binarised file is used as ground truth. Ported from new-modelling/src/default.py
# (coadd config/manual/coadd_strains.csv), plus spneumoniae.
# NOTE: efaecium and spneumoniae have CoAdd data only for MIC (no single-point inhibition file),
# so their inhib_50 evaluation is skipped and logged.
COADD_REF_STRAINS = {
    "abaumannii":  "ATCC19606",
    "calbicans":   "ATCC90028",
    "ecoli":       "ATCC25922",
    "efaecium":    "ATCC51559",   # MIC only
    "kpneumoniae": "ATCC700603",
    "paeruginosa": "ATCC27853",
    "saureus":     "ATCC43300",
    "spneumoniae": "ATCC700677",  # MIC only
}

# CoAdd ground-truth label columns (precomputed upstream; we only select which
# cutoff column to read, we do NOT re-binarise).
#   inhib_50 = active if >=50% single-point inhibition (matches ChEMBL SP=50%)
#   mic_10   = active if MIC <=10 uM  (matches ChEMBL DR=10 uM)
COADD_INHIB_COL = "inhib_50"
COADD_MIC_COL   = "mic_10"
# Endpoint label -> (binarised subdir, label column). Both are evaluated.
COADD_ENDPOINTS = {
    "inhib_50": ("03_binarised_inhibition", COADD_INHIB_COL),
    "mic_10":   ("05_binarised_mic",        COADD_MIC_COL),
}

# ---------------------------------------------------------------------------
# Inter-model correlation step (07)
# ---------------------------------------------------------------------------
# Annotation-model predictions on the shared reference library, staged by
# 00_download_data.py as data/processed/<subdir>/{model_id}_{version}.csv. Every file shares the
# join columns `key` (compound hash) and `input` (SMILES); output columns are model-specific.
ANNOTATION_PREDS_SUBDIR = "annotation_preds_ref_library"

# Correlation is computed on a fixed-seed sample of the ~1.35M-compound library (RANDOM_SEED) so
# the ~1500-node matrix is tractable. Spearman rho is essentially stable at this size; raise
# CORR_SAMPLE_N (or set it to None to use the full library) for the final, narrowed analysis.
CORR_SAMPLE_N = 200_000
# Rows read per chunk when filtering the large prediction files (one is 15.9 GB / 619 columns).
CORR_CHUNK_SIZE = 200_000

# Top-N hit-overlap depths (user-chosen). Overlap = Jaccard of the top-N highest-scoring compounds.
TOPN_CUTOFFS = (100, 500)

# Per-output-column value-type tagging (computed on the sample):
#   categorical  : integer-like with <= CATEGORICAL_MAX_UNIQUE distinct values
#   probability  : all values within [0, 1]  (a known "higher = more" direction)
#   continuous   : anything else (e.g. molecular weight, logP, free energy)
# Only probability columns receive top-N overlap; all columns enter the Spearman matrix.
# CATEGORICAL_MAX_UNIQUE is a display heuristic, not a data filter — revisable.
CATEGORICAL_MAX_UNIQUE = 10

# Case-insensitive regex that ADVISORILY flags a model as toxicity-related in the review CSV
# (07_group_assignments.csv, column `is_cytotox`), matched against the Airtable
# Tag / Title / Interpretation / Description fields. Broad on purpose — the precise cytotoxicity
# FOCUS GROUP is built by the two regexes below, not by this flag.
CYTOTOX_REGEX = r"cytotox|hepg2|cardiotox|herg|tox21|toxcast|\bdili\b|\btoxic(ity)?\b"

# Precise cytotoxicity focus-group membership (output-column level). A node joins the group when
# EITHER its OUTPUT-COLUMN name matches TOX_COLUMN_REGEX, OR its model TITLE matches the narrower
# TOX_TITLE_REGEX. The title rule rescues dedicated cell/organ-tox models whose single output is
# generically named (e.g. the hERG models eos30gr:activity_80, eos43at:pic50); the narrow title
# terms deliberately EXCLUDE broad "toxicity" panels (Tox21 nuclear-receptor assays, adverse-drug-
# reaction organ classes) and metabolic-stability / activity models that only mention toxicity.
TOX_COLUMN_REGEX = (
    r"cytotox|hepg2|hepatotox|cardiotox|herg|mitotox|toxicity|\btox\b|tox_|_tox|\bdili\b|dili_|_dili"
)
TOX_TITLE_REGEX = r"cytotox|hepg2|hepatotox|cardiotox|herg|\bdili\b"

# Organisms excluded from the per-organism focus heatmaps. 'Homo sapiens' is a catch-all of
# unrelated human endpoints (hERG, logP, DILI, permeability), not a meaningful same-organism set;
# its models still appear in the cytotoxicity group and the global heatmap.
ORGANISM_EXCLUDE = {"Homo sapiens"}

# Minimum number of models sharing an organism for that organism to get its own focus heatmap
# (a correlation needs at least two nodes). Display minimum, not a data filter.
ORGANISM_MIN_MODELS = 2