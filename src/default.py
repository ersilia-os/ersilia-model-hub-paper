
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

# Automation accounts to exclude from every community count. GitHub's ``type`` field reports all of
# these as "User", so they cannot be filtered programmatically and have to be listed. Note the
# converse trap: matching on the substring "bot" also catches real people (e.g. MuoboTone).
GITHUB_BOT_ACCOUNTS = {
    "dependabot[bot]", "github-actions[bot]", "model-request-bot[bot]",
    "vercel[bot]", "bitnami-bot",
}

# Known duplicate GitHub handles in the Airtable ``Contributor`` field, mapped to the account the
# person actually uses. Applied before any contributor is counted, so the same person is never
# counted twice. This is a data correction, not a display relabelling: without it the distinct
# contributor count is 31 rather than 30.
#   Richioo -> Richiio : one model (2023-12-03) filed under a handle that differs from the same
#                        contributor's other five (2023-08 to 2024-02) by a single character.
# The Airtable record should be corrected upstream; this mapping is the interim fix.
# NOTE (2026-08-04): no current consumer. Script 01 stopped counting Airtable contributors when the
# contributors-over-time panel was dropped, and step 08 counts people from the GitHub API instead.
# Kept as the record of the upstream duplicate, which is still uncorrected; re-apply it in any script
# that counts the Airtable ``Contributor`` field.
CONTRIBUTOR_ALIASES = {
    "Richioo": "Richiio",
}

# Shorter display labels for crowded metadata figure axes.
SUBTASK_DISPLAY = {
    "Property calculation or prediction": "Property prediction",
}
BIOAREA_DISPLAY = {
    "Antimicrobial resistance": "AMR",
    "Diarrheal diseases": "Diarrhoea",
}

# The Annotation subtask the grouped Biomedical Area panel is built from. Restricting to it is what
# keeps the "Other" bucket honest: `Any` (no disease area) is 22 of the 39 "Property calculation or
# prediction" models but only 4 of the 92 Activity prediction ones, because generic property
# predictors (logP, solubility) have no area to declare.
ACTIVITY_SUBTASK = "Activity prediction"

# Biomedical Area -> one of four groups, for the compact grouped panel. Signed off 2026-08-02.
#
# Every raw Airtable value must appear here: the counting step raises on an unmapped value rather than
# dropping it, so a new area added upstream fails loudly instead of vanishing from the figure.
#
# Membership was checked against each model's Target Organism, not inferred from the area name:
#   Peptic ulcer disease  -> eos9eyo, Helicobacter pylori (bacterial, not the NSAID aetiology)
#   Diarrheal diseases    -> all Campylobacter / E. coli, i.e. Gram-negative bacteria
#   Candidiasis, Mycetoma -> Candida albicans, Madurella mycetomatis (fungal; antifungal counts)
#
# TWO DELIBERATE STRETCHES, kept because Ersilia's own naming already treats them this way (the
# S. mansoni model's slug is literally `antimicrobial-activity-smansoni`):
#   Malaria         -> Plasmodium falciparum is a PROTOZOAN, inside "antimicrobial" only on the broad
#                      clinical definition.
#   Schistosomiasis -> Schistosoma mansoni is a multicellular HELMINTH, not a microorganism at all.
# A caption that says "antimicrobial" therefore covers antibacterial, antifungal, antiprotozoal and
# antihelminthic activity. Rename the group to "Anti-infective" if that overclaims for a given venue.
#
# `Any` (no area declared) goes to Other rather than being dropped, so the four groups account for
# every Activity prediction model. At this subtask it is only 4 models, so Other stays a genuine
# residual instead of the 26-model catch-all it would be across all of Annotation.
BIOAREA_GROUP = {
    "ADMET": "ADMET",
    "Antimicrobial resistance": "Antimicrobial",
    "Tuberculosis": "Antimicrobial",
    "Malaria": "Antimicrobial",
    "Pneumonia": "Antimicrobial",
    "Diarrheal diseases": "Antimicrobial",
    "Gonorrhea": "Antimicrobial",
    "Schistosomiasis": "Antimicrobial",
    "Candidiasis": "Antimicrobial",
    "Mycetoma": "Antimicrobial",
    "Peptic ulcer disease": "Antimicrobial",
    "COVID-19": "Antiviral",
    "AIDS": "Antiviral",
    "Hepatitis B": "Antiviral",
    "Cancer": "Other",
    "Alzheimer": "Other",
    "Any": "Other",
}

#: Catch-all group, pinned last in the panel however large it gets.
BIOAREA_GROUP_OTHER = "Other"

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

# Endpoint backing the cross-organism hit-set analyses of step 06 (overlap, promiscuity,
# exclusivity). inhib_50 is the CoAdd counterpart of the EU OpenScreen primary screen: ~81,600
# compounds tested per organism against ~4,500 for mic_10, so it is the only endpoint with enough
# actives to split by promiscuity. Only the 6 organisms with an inhibition file take part —
# efaecium and spneumoniae are MIC-only (see COADD_REF_STRAINS).
COADD_HITSET_ENDPOINT = "inhib_50"

# ---------------------------------------------------------------------------
# Score matrices and their analyses (steps 07-09)
# ---------------------------------------------------------------------------
# Annotation-model predictions on the shared reference library, staged by
# 00_download_data.py as data/processed/<subdir>/{model_id}_{version}.csv. Every file shares the
# join columns `key` (compound hash) and `input` (SMILES); output columns are model-specific.
ANNOTATION_PREDS_SUBDIR = "annotation_preds_ref_library"

# Rows read per chunk when filtering the large prediction files (one is 15.9 GB / 619 columns).
CORR_CHUNK_SIZE = 200_000

# Row-wise normalization (eval_correlations.row_normalize): a compound whose profile norm falls
# below this divides to inf/NaN. Not a data filter — such rows are counted and reported, never
# dropped or clipped; the value only sets what counts as "degenerate" for that warning.
ROW_NORM_EPS = 1e-12

# Endpoint-type classification used by tools/build_endpoint_selection_template.py when generating
# the reviewable template for config/08_endpoint_selection.csv, stratifying by whether an output
# column measures membrane PERMEABILITY/uptake/efflux/accumulation rather than direct
# antimicrobial/antiparasitic ACTIVITY (potency). The Airtable "Permeability" Tag alone under-counts
# this (only 5 of 11 permeability-type models carry that exact tag — e.g. "Efflux susceptibility in
# gram-negative bacteria" and "Entry-way rules (eNTRy) for gram-negative bacteria" are tagged only
# "Antimicrobial activity"), so this case-insensitive regex is matched over Title + Tag +
# Interpretation together: a hit means "permeability", otherwise "activity". Verified to reproduce
# the same 11 permeability / 25 activity split as a manual title-by-title read, with the matched
# evidence written out for review. Advisory only — the curated CSV is the source of truth.
ENDPOINT_TYPE_REGEX = r"permea|retain|accumulat|entry|efflux"

# ---------------------------------------------------------------------------
# Reference-library projection step (10)
# ---------------------------------------------------------------------------
# eos1klk (Task=Representation, Subtask=Projection) computes four 2D layouts of the SAME
# ~1.35M-compound reference library staged by 00_download_data.py — one output column pair
# ({method}_x, {method}_y) per method. Its Task is not "Annotation", so it is invisible to
# 00_download_data.py's automatic Section 4 loop and is fetched by an explicit call instead.
PROJECTION_MODEL_ID = "eos1klk"
PROJECTION_PREDS_SUBDIR = "eos1klk_projection"
PROJECTION_METHODS = ["pca", "umap", "tsne", "tmap"]

# Grid resolution (bins x bins) each projection is reduced to before plotting — a display
# resolution, not a data filter. 60 keeps a 15-pathogen small-multiples panel (~30x30 mm each)
# legible without over-resolving past what that panel size can show. Revisable.
PROJECTION_BINS = 60

# Number of highest-scoring compounds highlighted per pathogen, on top of the silver full-library
# density background. A user-directed value (not a fitted cutoff): the figure shows each
# pathogen's top PROJECTION_TOP_N compounds by consensus_score, ranked, rather than a score
# threshold — so it never needs a threshold sign-off, only a count.
PROJECTION_TOP_N = 1000