"""One-off generator for config/08_endpoint_selection.csv — the manually-curated master list of
which (model_id, column_name) endpoints enter the step-08 same-pathogen vs. different-pathogen
comparison.

This REPLACED an earlier automatic filter (minimum-models-per-organism plus single-organism-only,
both since removed) with full manual control: every endpoint whose model
targets at least one real pathogen (any Target Organism other than Homo sapiens / Rattus norvegicus
/ Mus musculus) gets a row, with organism / organism_class / sensitivity / assay_type / direction
best-effort INFERRED from Airtable metadata (Title, Tag, Interpretation) and, for a handful of
multi-organism models, from the column-name convention itself. Everything here is a starting point
for manual review, not a final answer — see the column-by-column notes below.

Run once: `python tools/build_endpoint_selection_template.py`. Refuses to overwrite an existing
config/08_endpoint_selection.csv (it is meant to be hand-edited after the first generation) — delete
it first if you deliberately want to regenerate from scratch.

Columns
-------
model_id       Ersilia eos identifier.
slug           Airtable "Slug" (human-readable model name), for readability when scanning the CSV.
column_name    The model's output column.
direction      "higher" (default) or "lower" — which direction means "more" for this column. Only
               one model in the whole catalog states an explicit inversion in its Interpretation
               text (eos5jv3: "lower values indicate higher permeability"); everything else defaults
               to "higher" since practically all outputs here are probabilities of activity/effect.
organism       Target organism for this specific COLUMN (not just the model — see below).
organism_class One of: Gram-negative bacteria, Gram-positive bacteria, Mycobacteria (acid-fast),
               Protozoa, Helminths, Fungi, Viruses, or blank for non-pathogen / not-organism-specific
               columns. Fungi and Viruses aren't in the 5 classes named in chat but the catalog does
               contain fungal (Candida, Cryptococcus, Madurella) and viral (SARS-CoV-2, HIV,
               Hepatitis B) pathogens, so they get their own class rather than being forced in.
sensitivity    "wild-type" (default) or "resistant", inferred from Interpretation/Description text
               or, for a few models, from the reference-strain code named in the column itself (see
               STRAIN_SENSITIVITY below — e.g. S. aureus ATCC 43300 is a standard MRSA reference
               strain, P. falciparum K1 is a standard chloroquine-resistant reference strain). A
               THIRD value, "sensitized", is used for the couple of models built on deliberately
               permeability/efflux-compromised strains (E. coli lpxC/tolC knockouts, a "sensitized"
               GNEProp strain) — biologically these are neither wild-type nor drug-resistant, so
               forcing them into the binary would misrepresent them; re-tag by hand if you'd rather
               collapse this into one of the two.
assay_type     Best-effort guess: bioactivity, permeability, toxicity, ADME, physicochemical,
               class_prediction, structural_alert. Built from ENDPOINT_TYPE_REGEX (src/default.py)
               plus explicit overrides for named non-bioactivity columns (cytotoxicity, hemolysis,
               clearance, CYP, solubility, Caco-2, structural alerts, the eos74km/eos6m2k
               class-level outputs).
selected       "Yes" only when organism/organism_class/sensitivity/assay_type are all unambiguous
               AND assay_type == "bioactivity" AND organism isn't a non-pathogen host. Everything
               else — including permeability endpoints, which are NOT ambiguous, just a different
               question — defaults to "No" per instruction: doubt (of any kind, including "this is a
               deliberately different assay type") means unselected until reviewed by hand.
"""

import os
import re
import sys

import pandas as pd

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))

from default import AIRTABLE_METADATA_FILE, ENDPOINT_TYPE_REGEX  # noqa: E402

meta_path = os.path.join(root, "..", "data", "raw", AIRTABLE_METADATA_FILE)
col_index_path = os.path.join(root, "..", "output", "07_prediction_correlations", "07_column_index.csv")
config_dir = os.path.join(root, "..", "config")
out_path = os.path.join(config_dir, "08_endpoint_selection.csv")

NON_PATHOGEN_HOSTS = {"Homo sapiens", "Rattus norvegicus", "Mus musculus"}

# Standard microbiological classification (textbook taxonomy) for every pathogen appearing anywhere
# in the catalog's Target Organism field. Gram stain only applies to true bacteria; M. tuberculosis
# is acid-fast (atypical mycolic-acid cell wall, conventionally its own category, not gram+/-).
ORGANISM_CLASS = {
    # Gram-negative bacteria
    "Escherichia coli": "Gram-negative bacteria", "Pseudomonas aeruginosa": "Gram-negative bacteria",
    "Acinetobacter baumannii": "Gram-negative bacteria", "Neisseria gonorrhoeae": "Gram-negative bacteria",
    "Klebsiella pneumoniae": "Gram-negative bacteria", "Helicobacter pylori": "Gram-negative bacteria",
    "Campylobacter spp": "Gram-negative bacteria", "Enterobacter spp": "Gram-negative bacteria",
    "Burkholderia cenocepacia": "Gram-negative bacteria", "Fusobacterium nucleatum": "Gram-negative bacteria",
    "Bacteroides caccae": "Gram-negative bacteria", "Bacteroides fragilis": "Gram-negative bacteria",
    "Bacteroides ovatus": "Gram-negative bacteria", "Bacteroides thetaiotaomicron": "Gram-negative bacteria",
    "Bacteroides uniformis": "Gram-negative bacteria", "Bacteroides vulgatus": "Gram-negative bacteria",
    "Bacteroides xylanisolvens": "Gram-negative bacteria", "Prevotella copri": "Gram-negative bacteria",
    "Parabacteroides distasonis": "Gram-negative bacteria", "Parabacteroides merdae": "Gram-negative bacteria",
    "Bilophila wadsworthia": "Gram-negative bacteria", "Odoribacter splanchnicus": "Gram-negative bacteria",
    "Veillonella parvula": "Gram-negative bacteria", "Akkermansia muciniphila": "Gram-negative bacteria",
    # Gram-positive bacteria
    "Staphylococcus aureus": "Gram-positive bacteria", "Streptococcus pneumoniae": "Gram-positive bacteria",
    "Streptococcus parasanguinis": "Gram-positive bacteria", "Streptococcus salivarius": "Gram-positive bacteria",
    "Enterococcus faecium": "Gram-positive bacteria", "Enterococcus faecalis": "Gram-positive bacteria",
    "Bifidobacterium adolescentis": "Gram-positive bacteria", "Bifidobacterium longum": "Gram-positive bacteria",
    "Blautia obeum": "Gram-positive bacteria", "Clostridium bolteae": "Gram-positive bacteria",
    "Clostridium difficile": "Gram-positive bacteria", "Clostridium perfringens": "Gram-positive bacteria",
    "Clostridium ramosum": "Gram-positive bacteria", "Clostridium saccharolyticum": "Gram-positive bacteria",
    "Collinsella aerofaciens": "Gram-positive bacteria", "Coprococcus comes": "Gram-positive bacteria",
    "Dorea formicigenerans": "Gram-positive bacteria", "Eggerthella lenta": "Gram-positive bacteria",
    "Eubacterium eligens": "Gram-positive bacteria", "Eubacterium rectale": "Gram-positive bacteria",
    "Lactobacillus paracasei": "Gram-positive bacteria", "Roseburia hominis": "Gram-positive bacteria",
    "Roseburia intestinalis": "Gram-positive bacteria", "Ruminococcus bromii": "Gram-positive bacteria",
    "Ruminococcus gnavus": "Gram-positive bacteria", "Ruminococcus torques": "Gram-positive bacteria",
    # Other classes
    "Mycobacterium tuberculosis": "Mycobacteria (acid-fast)",
    "Plasmodium falciparum": "Protozoa", "Leishmania major": "Protozoa",
    "Schistosoma mansoni": "Helminths",
    "Candida albicans": "Fungi", "Candida glabrata": "Fungi", "Cryptococcus neoformans": "Fungi",
    "Cryptococcus deuterogattii": "Fungi", "Madurella mycetomatis": "Fungi",
    "SARS-CoV-2": "Viruses", "HIV": "Viruses", "Hepatitis B virus": "Viruses",
}

# Reference-strain codes with a well-documented resistance/sensitization phenotype, keyed by
# substring match on the column name (case-insensitive). Everything not listed defaults to
# "wild-type". Sources: ATCC 43300 (S. aureus, MRSA reference), P. falciparum K1 (chloroquine-
# resistant reference strain) are standard, widely-cited antimicrobial-screening reference strains;
# lpxC/tolC E. coli knockouts and PAO397 are efflux/outer-membrane mutants used specifically to
# remove the permeability barrier (hypersusceptible "sensitized" strains, not drug-resistant ones).
STRAIN_SENSITIVITY = {
    "atcc43300": "resistant",      # S. aureus ATCC 43300 = MRSA reference strain
    "pf_k1": "resistant",          # P. falciparum K1 = chloroquine-resistant reference strain
    "ecoli_lpxc": "sensitized",    # outer-membrane permeability mutant
    "ecoli_tolc": "sensitized",    # efflux-deficient mutant
    "pao397": "sensitized",        # P. aeruginosa efflux mutant
}
SENSITIZED_MODEL_KEYWORDS = re.compile(r"sensitiz", re.IGNORECASE)  # e.g. eos5nqn "sensitized strain"
RESISTANT_MODEL_KEYWORDS = re.compile(r"drug[- ]resistant|resistant bacteria", re.IGNORECASE)  # eos5xng

# Non-bioactivity assay-type overrides by exact column name (host-toxicity / ADME / physicochemical
# / class-level / structural-alert columns that would otherwise default to "bioactivity").
ASSAY_TYPE_COLUMN_OVERRIDE = {
    "cytotoxicity_ic50": "toxicity", "hemolitic_activity": "toxicity",
    "cho": "toxicity", "cho_norm": "toxicity", "hepg2": "toxicity", "hepg2_norm": "toxicity",
    "clint_h": "ADME", "clint_m": "ADME", "clint_r": "ADME",
    "clint_h_norm": "ADME", "clint_m_norm": "ADME", "clint_r_norm": "ADME",
    "cyp2c9": "ADME", "cyp2c19": "ADME", "cyp3a4": "ADME", "cyp2d6": "ADME",
    "cyp2c9_norm": "ADME", "cyp2c19_norm": "ADME", "cyp3a4_norm": "ADME", "cyp2d6_norm": "ADME",
    "caco_2": "permeability", "caco_2_norm": "permeability",  # Caco-2 IS a permeability assay
    "aq_sol": "physicochemical", "aq_sol_norm": "physicochemical",
    "has_pains": "structural_alert", "has_brenk": "structural_alert",
    "is_sim_known_ab": "structural_alert", "nitrofuran_motif": "structural_alert",
    "fluoroquinolone_motif": "structural_alert", "carbepenem_motif": "structural_alert",
    "betalactam_motif": "structural_alert",
    "apscore_total": "class_prediction", "apscore_gpositive": "class_prediction",
    "apscore_gnegative": "class_prediction",
    "acid_fast": "class_prediction", "fungi": "class_prediction", "gram_negative": "class_prediction",
    "gram_positive": "class_prediction", "inactive": "class_prediction",
}

# Explicit per-column organism overrides for the 7 models whose Target Organism field lists more
# than one organism, so a genuinely multi-organism MODEL still gets an unambiguous organism per
# COLUMN wherever the column name/semantics make that clear. Built by reading each model's column
# list against its title/interpretation (see the tool's docstring); anything genuinely ambiguous
# (e.g. eos2xeq's structural-alert flags, eos74km's class-probability outputs) is left with organism
# "" (blank) rather than guessed.
MULTI_ORGANISM_OVERRIDE = {
    "eos2xeq": {c: "" for c in [  # structural-alert flags, not organism-specific despite the tag
        "has_pains", "has_brenk", "is_sim_known_ab", "nitrofuran_motif",
        "fluoroquinolone_motif", "carbepenem_motif", "betalactam_motif"]},
    "eos3dys": {
        "abaumannii_ATCC19606_inhib_50": "Acinetobacter baumannii",
        "abaumannii_ATCC19606_mic_25": "Acinetobacter baumannii",
        "calbicans_ATCC90028_inhib_50": "Candida albicans", "calbicans_ATCC90028_mic_25": "Candida albicans",
        "cdeuterogattii_CBS7750_mic_25": "Cryptococcus deuterogattii",
        "cglabrata_ATCC90030_mic_25": "Candida glabrata",
        "cneoformans_H99_inhib_50": "Cryptococcus neoformans", "cneoformans_H99_mic_25": "Cryptococcus neoformans",
        "ecoli_ATCC25922_inhib_50": "Escherichia coli", "ecoli_ATCC25922_mic_25": "Escherichia coli",
        "ecoli_lpxC_inhib_50": "Escherichia coli", "ecoli_tolC_inhib_50": "Escherichia coli",
        "kpneumoniae_ATCC700603_inhib_50": "Klebsiella pneumoniae",
        "kpneumoniae_ATCC700603_mic_25": "Klebsiella pneumoniae",
        "paeruginosa_ATCC27853_inhib_50": "Pseudomonas aeruginosa",
        "paeruginosa_ATCC27853_mic_25": "Pseudomonas aeruginosa",
        "paeruginosa_PAO397_inhib_50": "Pseudomonas aeruginosa", "paeruginosa_PAO397_mic_25": "Pseudomonas aeruginosa",
        "saureus_ATCC43300_inhib_50": "Staphylococcus aureus", "saureus_ATCC43300_mic_25": "Staphylococcus aureus",
        "cytotoxicity_ic50": "Homo sapiens", "hemolitic_activity": "Homo sapiens",
    },
    "eos60mw": {
        "leishmania_rf": "Leishmania major", "leishmania_mlp": "Leishmania major",
        "leishmania_chemberta": "Leishmania major",
        "coronavirus_gcn": "SARS-CoV-2", "coronavirus_gb": "SARS-CoV-2", "coronavirus_chemberta": "SARS-CoV-2",
    },
    "eos6m2k": {
        "apscore_total": "", "apscore_gpositive": "", "apscore_gnegative": "",
        "akkermansia_muciniphila_nt5021": "Akkermansia muciniphila",
        "bacteroides_caccae_nt5050": "Bacteroides caccae",
        "bacteroides_fragilis_et_nt5033": "Bacteroides fragilis", "bacteroides_fragilis_nt_nt5003": "Bacteroides fragilis",
        "bacteroides_ovatus_nt5054": "Bacteroides ovatus",
        "bacteroides_thetaiotaomicron_nt5004": "Bacteroides thetaiotaomicron",
        "bacteroides_uniformis_nt5002": "Bacteroides uniformis", "bacteroides_vulgatus_nt5001": "Bacteroides vulgatus",
        "bacteroides_xylanisolvens_nt5064": "Bacteroides xylanisolvens",
        "bifidobacterium_adolescentis_nt5022": "Bifidobacterium adolescentis",
        "bifidobacterium_longum_nt5028": "Bifidobacterium longum",
        "bilophila_wadsworthia_nt5036": "Bilophila wadsworthia", "blautia_obeum_nt5069": "Blautia obeum",
        "clostridium_bolteae_nt5026": "Clostridium bolteae", "clostridium_difficile_nt5083": "Clostridium difficile",
        "clostridium_perfringens_nt5032": "Clostridium perfringens",
        "clostridium_ramosum_nt5006": "Clostridium ramosum",
        "clostridium_saccharolyticum_nt5037": "Clostridium saccharolyticum",
        "collinsella_aerofaciens_nt5073": "Collinsella aerofaciens", "coprococcus_comes_nt5048": "Coprococcus comes",
        "dorea_formicigenerans_nt5076": "Dorea formicigenerans", "eggerthella_lenta_nt5024": "Eggerthella lenta",
        "escherichia_coli_ed1a_nt5078": "Escherichia coli", "escherichia_coli_iai1_nt5077": "Escherichia coli",
        "eubacterium_eligens_nt5075": "Eubacterium eligens", "eubacterium_rectale_nt5009": "Eubacterium rectale",
        "fusobacterium_nucleatum_nt5025": "Fusobacterium nucleatum",
        "lactobacillus_paracasei_nt5042": "Lactobacillus paracasei",
        "odoribacter_splanchnicus_nt5081": "Odoribacter splanchnicus",
        "parabacteroides_distasonis_nt5074": "Parabacteroides distasonis",
        "parabacteroides_merdae_nt5071": "Parabacteroides merdae", "prevotella_copri_nt5019": "Prevotella copri",
        "roseburia_hominis_nt5079": "Roseburia hominis", "roseburia_intestinalis_nt5011": "Roseburia intestinalis",
        "ruminococcus_bromii_nt5045": "Ruminococcus bromii", "ruminococcus_gnavus_nt5046": "Ruminococcus gnavus",
        "ruminococcus_torques_nt5047": "Ruminococcus torques",
        "streptococcus_parasanguinis_nt5072": "Streptococcus parasanguinis",
        "streptococcus_salivarius_nt5038": "Streptococcus salivarius",
        "veillonella_parvula_nt5017": "Veillonella parvula",
    },
    "eos74km": {c: "" for c in ["acid_fast", "fungi", "gram_negative", "gram_positive", "inactive"]},
    "eos7kpb": {
        "pf_nf54": "Plasmodium falciparum", "pf_nf54_norm": "Plasmodium falciparum",
        "pf_k1": "Plasmodium falciparum", "pf_k1_norm": "Plasmodium falciparum",
        "mtb": "Mycobacterium tuberculosis", "mtb_norm": "Mycobacterium tuberculosis",
        "cho": "Homo sapiens", "cho_norm": "Homo sapiens", "hepg2": "Homo sapiens", "hepg2_norm": "Homo sapiens",
        "clint_h": "Homo sapiens", "clint_h_norm": "Homo sapiens",
        "clint_m": "Mus musculus", "clint_m_norm": "Mus musculus",
        "clint_r": "Rattus norvegicus", "clint_r_norm": "Rattus norvegicus",
        "caco_2": "Homo sapiens", "caco_2_norm": "Homo sapiens",
        "aq_sol": "", "aq_sol_norm": "",
        "cyp2c9": "Homo sapiens", "cyp2c19": "Homo sapiens", "cyp3a4": "Homo sapiens", "cyp2d6": "Homo sapiens",
        "cyp2c9_norm": "Homo sapiens", "cyp2c19_norm": "Homo sapiens", "cyp3a4_norm": "Homo sapiens",
        "cyp2d6_norm": "Homo sapiens",
    },
    "eos9x3z": {"gn_activity": "Escherichia coli"},  # "Gram-negative activity...proxy"; tagged E.coli+S.aureus
}
MULTI_ORGANISM_ORGANISM_CLASS_OVERRIDE = {
    ("eos6m2k", "apscore_gpositive"): "Gram-positive bacteria",
    ("eos6m2k", "apscore_gnegative"): "Gram-negative bacteria",
    ("eos74km", "acid_fast"): "Mycobacteria (acid-fast)",
    ("eos74km", "fungi"): "Fungi",
    ("eos74km", "gram_negative"): "Gram-negative bacteria",
    ("eos74km", "gram_positive"): "Gram-positive bacteria",
}

PERM_REGEX = re.compile(ENDPOINT_TYPE_REGEX, re.IGNORECASE)
LOWER_DIRECTION_REGEX = re.compile(r"lower value.*(higher|more|greater)", re.IGNORECASE)


def split_organisms(raw):
    if not isinstance(raw, str) or not raw.strip():
        return []
    return [o.strip() for o in raw.split(",") if o.strip() and o.strip().lower() != "any"]


def infer_sensitivity(model_id, column, blob):
    key = column.lower()
    for code, val in STRAIN_SENSITIVITY.items():
        if code in key:
            return val
    if SENSITIZED_MODEL_KEYWORDS.search(blob):
        return "sensitized"
    if RESISTANT_MODEL_KEYWORDS.search(blob):
        return "resistant"
    return "wild-type"


def infer_assay_type(column, blob):
    if column in ASSAY_TYPE_COLUMN_OVERRIDE:
        return ASSAY_TYPE_COLUMN_OVERRIDE[column]
    if PERM_REGEX.search(blob):
        return "permeability"
    return "bioactivity"


def infer_direction(interpretation):
    if isinstance(interpretation, str) and LOWER_DIRECTION_REGEX.search(interpretation):
        return "lower"
    return "higher"


def main():
    if os.path.exists(out_path):
        print(f"[skip] {out_path} already exists — delete it first to regenerate from scratch.")
        return

    meta = pd.read_csv(meta_path).fillna("")
    col_index = pd.read_csv(col_index_path)
    meta_by_id = meta.set_index("Identifier")

    rows = []
    for model_id, model_cols in col_index.groupby("model_id"):
        if model_id not in meta_by_id.index:
            continue
        m = meta_by_id.loc[model_id]
        organisms = split_organisms(m["Target Organism"])
        pathogens = [o for o in organisms if o not in NON_PATHOGEN_HOSTS]
        if not pathogens:
            continue  # no real pathogen target at all — out of scope for this template

        slug = m["Slug"]
        blob_model = " ".join([str(m["Title"]), str(m["Tag"]), str(m["Interpretation"])])
        direction = infer_direction(m["Interpretation"])
        is_multi = model_id in MULTI_ORGANISM_OVERRIDE

        for _, r in model_cols.iterrows():
            column = r["output_col"]
            if is_multi:
                organism = MULTI_ORGANISM_OVERRIDE[model_id].get(column, "")
            else:
                organism = pathogens[0] if len(pathogens) == 1 else ""

            organism_class = MULTI_ORGANISM_ORGANISM_CLASS_OVERRIDE.get(
                (model_id, column), ORGANISM_CLASS.get(organism, ""))
            sensitivity = infer_sensitivity(model_id, column, blob_model) if organism \
                and organism not in NON_PATHOGEN_HOSTS else ""
            assay_type = infer_assay_type(column, blob_model)

            is_clean_pathogen = bool(organism) and organism not in NON_PATHOGEN_HOSTS
            selected = "Yes" if (is_clean_pathogen and organism_class and assay_type == "bioactivity"
                                 and sensitivity in ("wild-type", "")) else "No"

            rows.append({
                "model_id": model_id, "slug": slug, "column_name": column, "direction": direction,
                "organism": organism, "organism_class": organism_class, "sensitivity": sensitivity,
                "assay_type": assay_type, "selected": selected,
            })

    df = pd.DataFrame(rows, columns=["model_id", "slug", "column_name", "direction", "organism",
                                     "organism_class", "sensitivity", "assay_type", "selected"])
    df = df.sort_values(["organism_class", "organism", "model_id", "column_name"],
                        na_position="last").reset_index(drop=True)
    os.makedirs(config_dir, exist_ok=True)
    df.to_csv(out_path, index=False)
    print(f"[build] wrote {len(df)} endpoint rows across {df['model_id'].nunique()} models -> {out_path}")
    print(f"[build] selected=Yes by default: {(df['selected'] == 'Yes').sum()}")
    print(df["organism_class"].value_counts(dropna=False))
    print(df["assay_type"].value_counts(dropna=False))


if __name__ == "__main__":
    main()
