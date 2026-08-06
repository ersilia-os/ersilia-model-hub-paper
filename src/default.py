
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

# ---------------------------------------------------------------------------
# Basic physicochemical descriptors (step 11)
# ---------------------------------------------------------------------------
# eos4djh (datamol-basic-descriptors, Task=Representation, Subtask=Featurization) computes 22
# basic descriptors — molecular weight, cLogP, H-bond donors/acceptors, ring and atom counts —
# over the same reference library. These are deterministic RDKit-family calculations wrapped in
# the Datamol API, not predictions, so they carry no training set and no leakage dimension.
# Like PROJECTION_MODEL_ID above, its Task is not "Annotation", so 00_download_data.py's
# automatic Section 4 loop never sees it and it is fetched by an explicit call instead. Its
# predictions are nonetheless filed alongside the annotation models, since step 11 reads them
# from that one folder.
PHYSCHEM_MODEL_ID = "eos4djh"

# Grid resolution (bins x bins) each projection is reduced to before plotting — a display
# resolution, not a data filter. 60 keeps a 15-pathogen small-multiples panel (~30x30 mm each)
# legible without over-resolving past what that panel size can show. Revisable.
PROJECTION_BINS = 60

# Number of highest-scoring compounds highlighted per pathogen, on top of the silver full-library
# density background. A user-directed value (not a fitted cutoff): the figure shows each
# pathogen's top PROJECTION_TOP_N compounds by consensus_score, ranked, rather than a score
# threshold — so it never needs a threshold sign-off, only a count.
PROJECTION_TOP_N = 1000

# ---------------------------------------------------------------------------
# Step 13 — toxicity endpoints projected onto the same reference-library layout.
# Only UMAP of the four PROJECTION_METHODS is drawn: 24 endpoint panels at four methods each
# would be 96 panels for no added insight, and UMAP is the layout the pathogen figures are read
# from. Revisable — the step-13 engine takes the method as an argument.
TOX_PROJECTION_METHOD = "umap"

# Column-name prefix for the step-11 cytotoxicity table, matching the
# {prefix}__{model_id}__{column_name} convention shared with the bioactivity matrix
# (build_named_library_matrix in eval_correlations) and antibiotic resemblance.
TOX_PREFIX = "cytotox"

# Every selected endpoint in config/cytotoxicity_models.csv is ranked HIGHEST-FIRST to define
# "most toxic". For the 23 classification endpoints this is the predicted probability of the
# toxic/active class, so it is direct. ld50_zhu is the one regression endpoint: eos7m30 emits it
# in log(1/(mol/kg)) units, where the reciprocal inverts the dose scale, so a HIGHER value is a
# LOWER LD50 and therefore MORE acutely toxic. The eos7m30 column metadata labels this endpoint
# `direction: low`, which describes the underlying LD50 dose rather than the transformed value it
# actually emits; user-confirmed to rank it highest-first like the rest. Corroborated by its
# positive Spearman correlation (+0.21 to +0.30) with all six independently-trained cytotoxicity
# models (eos42ez x3, eos3le9 x2, eos3dys) on a 150k-compound sample.
TOX_RANK_DESCENDING = True

# ---------------------------------------------------------------------------
# Step 14 — pathogen hits vs antibiotic-resemblance hits on the library UMAP.
# The subset of config/antibiotic_resemblance.csv's 55 selected endpoints drawn against every
# pathogen: nine spanning the three kinds of evidence — one continuous learned score, two
# similarity counts against AntibioticDB, and six substructure/class flags. User-directed.
#
# Two names differ from how they were requested and are corrected here to the config's spelling:
#   betalactan_motif -> betalactam_motif (eos2xeq)      b_lactam_all -> b_lactams_all (eos19mt)
#
# Only abx_score is continuous. The rest are binary flags or small integer counts, so step 12's
# highlight rule (every compound > 0, capped at PROJECTION_TOP_N) is what defines the abx side —
# NOT a rank cutoff. Three of these nine have fewer than PROJECTION_TOP_N positives library-wide
# (carbepenem_motif 346, ansamycins_rifamycins_macrolides 577, b_lactams_all 733) and so are drawn
# exhaustively; the rest hit the cap and are an arbitrary subset of their positives.
ABX_OVERLAP_ENDPOINTS = [
    "abx_score",
    "num_sim_0_5_all",
    "num_sim_0_5_subset",
    "fluoroquinolone_motif",
    "nitrofuran_motif",
    "betalactam_motif",
    "b_lactams_all",
    "carbepenem_motif",
    "ansamycins_rifamycins_macrolides",
]
# ---------------------------------------------------------------------------
# Step 13 — property/resemblance columns as predictors of pathogen activity.
#
# Each activity endpoint (the selected == "Yes" rows of config/08_endpoint_selection.csv, all of
# them direction == "higher") is binarized into a two-class target: its ACTIVITY_BINARIZE_TOP_N
# highest-scoring compounds are the positive class, every remaining compound is negative. A
# user-directed RANK cutoff on a fixed count, not a score threshold — so no score value is ever
# chosen or reviewed, and prevalence is a constant across all targets.
#
# Deliberately a SEPARATE constant from PROJECTION_TOP_N despite sharing the value 1000: that one
# picks how many points a projection figure draws, this one defines the positive class of a
# classification target. They answer to different reviews and should be free to diverge.
ACTIVITY_BINARIZE_TOP_N = 1000

# The three predictor blocks, in figure order. Each is one column-name prefix of the
# {family}__{model_id}__{column_name} convention, and each gets its own step-13 figure — 101
# predictors on a single shared x-axis would leave ~1.8 mm per category at the 180 mm print width.
PREDICTOR_FAMILIES = ["physchem", "cytotox", "abx"]

# Metric per predictor, selected from the predictor's own value type (resolved on the FULL column,
# never a subsample):
#   continuous -> AUROC             ranks the predictor against the binary target
#   binary     -> balanced accuracy AUROC of a two-valued score is a degenerate single-operating-
#                 point measure, so the standard imbalanced-data pairing is used instead
# Both share a 0.5 chance baseline, which is the only reason the two can sit on one y-axis.
PREDICTOR_METRICS = {"continuous": "auroc", "binary": "balanced_accuracy"}
PREDICTOR_CHANCE_LEVEL = 0.5

# ---------------------------------------------------------------------------
# Step 15, pathogen-subset variant — the activity endpoints of the 15 pathogens of interest only,
# with each ChEMBL antimicrobial model reduced to its single consensus column.
#
# Those models contribute up to 52 highly correlated sub-endpoints each (one per source ChEMBL/
# PubChem assay), which dominate any pooled view by sheer count. Where the model publishes a
# CONSENSUS_COLUMN, that one column is its headline score and the sub-endpoints are collapsed into
# it; where it does not (a model with a single sub-model never had anything to take a consensus
# over — campylobacter/eos7iak, hpylori/eos9eyo, ngonorrhoeae/eos5qya), its endpoints are kept as
# they are. The rule is applied per (model_id, organism), not per model, because eos3dys spans six
# organisms. This cuts 214 endpoints to 59 across the 15 pathogens.
CONSENSUS_COLUMN = "consensus_score"

# config/pathogens_of_interest.csv and config/08_endpoint_selection.csv spell two organisms
# differently. An explicit alias map, never genus-substring matching: "Candida albicans" would
# otherwise capture "Candida glabrata", and "Streptococcus pneumoniae" would capture
# S. parasanguinis and S. salivarius — all three distinct organisms in the curation.
PATHOGEN_ORGANISM_ALIASES = {
    "Campylobacter": "Campylobacter spp",
    "Enterobacter": "Enterobacter spp",
}

# The curated 12-predictor subset scored against the pathogen-subset targets (step 13's sixth
# figure): a hand-picked, user-directed shortlist of the 101 property columns, three families in one
# panel. All twelve are continuous, so the whole panel is AUROC on one scale — unlike the per-family
# figures, where the abx block mixes AUROC and balanced accuracy.
#   physchem : the three descriptors a medicinal chemist reads first
#   cytotox  : the six direct cell-viability readouts (the Tox21/ADMET panels are left out)
#   abx      : the learned antibiotic-likeness score plus the two 0.5-cutoff similarity counts
# Keyed by family so the figure can colour by it; values are the bare column_name, resolved against
# the {family}__{model_id}__{column_name} names at read time so a model version bump needs no edit.
CURATED_PREDICTORS = {
    "physchem": ["mw", "tpsa", "clogp"],
    "cytotox": ["cytotoxicity_hepg2", "cytotoxicity_hskmc", "cytotoxicity_imr90",
                "cytotoxicity_ic50", "ic50_hepg2_72h_5um", "ic50_hepg2_72h_10um"],
    "abx": ["abx_score", "num_sim_0_5_all", "num_sim_0_5_subset"],
}

# Family -> hue for that panel. Three well-separated hues; the legend is only three entries, so
# colour is a genuine primary encoding here rather than the secondary cue it has to be at 15 groups.
CURATED_FAMILY_HUES = {"physchem": "turquoise", "cytotox": "crimson", "abx": "cobalt"}

# ---------------------------------------------------------------------------
# Step 14 — the AUROC matrix (predictors x pathogen activity endpoints).
#
# organism_class has no natural order, so the axis order is DECLARED here: the three bacterial
# classes together first, then the eukaryotic pathogens roughly by cell complexity. This is a
# display choice and not a claim about the organisms — change the list to reorder the figure.
# Every organism_class value present in the curation must appear, or the ordering step raises
# rather than silently dropping a block of endpoints off the end of the axis.
ORGANISM_CLASS_ORDER = [
    "Gram-negative bacteria",
    "Gram-positive bacteria",
    "Mycobacteria (acid-fast)",
    "Fungi",
    "Protozoa",
    "Helminths",
]

# Discrete AUROC bins, 0.1 wide, running 0.2 -> 1.0 and CENTRED ON 0.5.
#
# AUROC has a principled neutral: 0.5 is chance. 44 of the 405 cells (10.9%) fall below it, down to
# 0.262 — those predictors rank actives BELOW inactives, which is real signal with a direction, not
# absence of signal. An earlier version clipped everything under 0.5 into one grey bin, which erased
# that. The floor is 0.2 because the matrix minimum is 0.262, so nothing is actually clipped.
AUROC_MATRIX_BINS = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

# The value the diverging colour scale pivots on. Cells near it render near-white, so "no better than
# chance" reads as visually empty, and the two directions get opposite hues.
AUROC_MATRIX_CENTER = 0.5

# stylia's NPG DIVERGING colormap, reversed so it runs cool -> warm: cobalt (blue) below chance,
# near-white at chance, crimson (red) above. A diverging map is the right family here precisely
# BECAUSE 0.5 is a meaningful midpoint with data on both sides — it is the wrong choice for a
# quantity with no natural centre, where the pale middle would land on an arbitrary value.
#
# This replaced, in order: a SpectralColormap("npg") ramp (hues carry no rank — magenta vs blue vs
# red gives a reader no way to tell high from low) and a monotonic ContinuousColormap("plum") ramp
# (correctly ordered but visually poor). The arms are asymmetric (0.24 below chance, 0.50 above), so
# the bins are mapped through a two-slope normalisation that pins 0.5 to the colormap's centre —
# without it the white point would drift off chance and the figure would imply a neutral at 0.6.
AUROC_MATRIX_CMAP = "crimson_cobalt"

# Column blocks of the matrix, in x-axis order: (block label, family, [column_name, ...]).
# The bioactivity block is not listed — it is the 59 activity endpoints, ordered separately.
# Order within each block is user-directed: cytotoxicity by model_id, the other two as given.
AUROC_MATRIX_BLOCKS = [
    ("cytotoxicity", "cytotox", ["cytotoxicity_ic50",            # eos3dys
                                 "ic50_hepg2_72h_5um", "ic50_hepg2_72h_10um",   # eos3le9
                                 "cytotoxicity_hepg2", "cytotoxicity_hskmc",
                                 "cytotoxicity_imr90"]),         # eos42ez
    ("abx resemblance", "abx", ["abx_score", "num_sim_0_5_all", "num_sim_0_5_subset"]),
    ("physchem", "physchem", ["mw", "clogp", "tpsa"]),
]

# Step 16's per-organism merge: each organism's endpoints are collapsed into ONE score per compound
# before any AUROC is computed, so a pathogen's weight in the figure comes from the pathogen and not
# from how many assays it happens to have (P. falciparum has 13 endpoints, H. pylori has 1).
#
# rank_pct, not zscore: percentile ranks are bounded [0, 1] and robust to outliers and to endpoints
# on unrelated native scales, so no single endpoint dominates an organism's mean. Unbounded z-scores
# would let one long-tailed endpoint carry the 11-endpoint E. coli and 13-endpoint P. falciparum
# aggregates. This is also the scaling step 07's mean-rank section already averages over.
ORGANISM_MERGE_METHOD = "rank_pct"

# Mean, not max or median: a compound must rank well ACROSS an organism's assays rather than in any
# one of them. Max would let a single noisy endpoint set the organism's whole column.
#
# Two consequences that are recorded rather than corrected, both following directly from "merge the
# endpoints currently in the grid":
#   - 5 organisms (Campylobacter, Enterobacter, E. faecium, H. pylori, S. pneumoniae) have exactly
#     ONE endpoint, so nothing is merged and their score IS that endpoint's percentile rank. Their
#     row is not the same kind of quantity as E. coli's.
#   - A ChEMBL `consensus_score` is itself an aggregate over sub-models, so averaging it with
#     individual assay endpoints gives it equal weight to a single assay.
ORGANISM_MERGE_AGG = "mean"

# ---------------------------------------------------------------------------
# Step 16's second view: how many of the row's actives fall in the column's top N.
#
# Different question from the AUROC matrix. AUROC asks "does this predictor RANK the row organism's
# actives highly, across the whole library"; this asks "how many of the row's 1000 actives are among
# the column's own top 1000". A predictor can do the first well without doing the second.
#
# The raw intersection count (0-1000), not Jaccard. Both sets have exactly ACTIVITY_BINARIZE_TOP_N
# members, so Jaccard = i / (2N - i) is a monotone re-expression of the same number and orders the
# matrix identically — the count is simply the one a reader can act on ("724 of the 1000 shared").
# Note the measure is SYMMETRIC, so the bioactivity block is a symmetric matrix, unlike AUROC's.
#
# NON-UNIFORM bins tuned to the distribution: off-diagonal counts run 0-724 with a MEDIAN of 3, and
# two random 1000-compound sets out of 1,355,109 would share ~0.7 by chance. Uniform bins would put
# most of the matrix in one class. The boundary at 1 is the informative one — it separates the 33% of
# cells with NO shared compound at all from those with some. No bin holds more than a third.
OVERLAP_MATRIX_BINS = [0, 1, 10, 25, 50, 100, 200, 400, 750]

# SEQUENTIAL, not diverging: unlike AUROC, an overlap count has no meaningful midpoint — its neutral
# is 0, at the end of the scale. Cobalt matches the overlap heatmaps in step 08 and the EU OpenScreen
# validation, so every set-overlap figure in the paper reads on one hue.
OVERLAP_MATRIX_HUE = "cobalt"

# The self-overlap diagonal is 1000/1000 by construction while the largest real value is 724, so it is
# BLANKED in the figure (dashed cell) and excluded from the column means. Follows the convention in
# step 08 and ActiveOverlapHeatmapPlot. True values stay in the CSV — blanking is a display choice.
OVERLAP_BLANK_DIAGONAL = True
