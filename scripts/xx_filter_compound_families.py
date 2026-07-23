"""
Annotate the Ersilia reference library and DrugBank SMILES with known antimicrobial
chemical family membership using SMARTS substructure searches.

Each compound receives a boolean column per family (non-exclusive — a compound may
belong to multiple families).

Families (12):
  beta_lactam, tetracycline, fluoroquinolone, sulfonamide, oxazolidinone,
  nitroimidazole, rifamycin, phenicol, quinoline, nitrofuran,
  macrolide, diaminopyrimidine

Macrolides are detected programmatically (RDKit SMARTS cannot constrain ring size):
a compound is classified as a macrolide if it contains an ester (C(=O)O) where both
the carbonyl carbon and the ester oxygen belong to the same ring of size >= 12.

Requires data/compound_lists/reference_library_smiles.csv and
data/compound_lists/drugbank_smiles.csv (run 00_download_data.py first).

Outputs:
  - output/xx_compound_families/reference_library_families.csv
  - output/xx_compound_families/drugbank_families.csv
"""

import os
import sys

import pandas as pd
from rdkit import Chem

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))

data_dir = os.path.join(root, "..", "data","raw", "compound_lists")
output_dir = os.path.join(root, "..", "output", "xx_compound_families")
os.makedirs(output_dir, exist_ok=True)

# ---------------------------------------------------------------------------
# SMARTS patterns — one entry per family (None = handled programmatically)
# ---------------------------------------------------------------------------

FAMILY_SMARTS = {
    # 4-membered beta-lactam ring (azetidinone)
    "beta_lactam": "[#6]1(=O)[#7][#6][#6]1",
    # keto-enol-amide triad of the A-ring (common to all tetracyclines/glycylcyclines)
    "tetracycline": "C(O)=C(C(N)=O)C(=O)",
    # 4-oxo-3-carboxylic acid quinolone (F presence enforced separately in classify())
    "fluoroquinolone": "O=c1c(C(=O)O)cnc2ccccc21",
    # sulfonamide SO2N core
    "sulfonamide": "[#16](=O)(=O)[#7]",
    # 2-oxazolidinone 5-membered ring: N-C(=O)-O-C-C
    "oxazolidinone": "[#7]1[#6](=O)[#8][#6][#6]1",
    # nitro group directly on imidazole ring (branch syntax required)
    "nitroimidazole": "[N+](=O)([O-])c1cnc[n]1",
    # aminonaphthalenediol core shared by all rifamycins
    "rifamycin": "Nc1cc(O)c2ccccc2c1O",
    # dichloroacetamide — defines the entire phenicol class
    "phenicol": "ClC(Cl)C(=O)[NH]",
    # quinoline ring system
    "quinoline": "c1ccc2ncccc2c1",
    # nitro group on furan ring (branch syntax required)
    "nitrofuran": "[N+](=O)([O-])c1ccco1",
    # programmatic — macrolactone ring size >= 12; see is_macrolide()
    "macrolide": None,
    # 2,4-diaminopyrimidine (trimethoprim scaffold)
    "diaminopyrimidine": "Nc1nccc(N)n1",
}

# Pre-compile all SMARTS patterns once
_COMPILED = {}
for _family, _smarts in FAMILY_SMARTS.items():
    if _smarts is not None:
        _pat = Chem.MolFromSmarts(_smarts)
        if _pat is None:
            raise ValueError(f"Invalid SMARTS for {_family}: {_smarts}")
        _COMPILED[_family] = _pat

# SMARTS for ester group (used in macrolide detection)
_ESTER_SMARTS = Chem.MolFromSmarts("C(=O)O")
# SMARTS for fluorine (used in fluoroquinolone check)
_F_SMARTS = Chem.MolFromSmarts("[F]")


def is_macrolide(mol):
    """Return True if the molecule contains a macrolactone ring (size >= 12)."""
    if mol is None:
        return False
    ester_matches = mol.GetSubstructMatches(_ESTER_SMARTS)
    if not ester_matches:
        return False
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    for match in ester_matches:
        # match = (carbonyl_C, oxo_O, ester_O)
        c_idx, _, o_idx = match
        for ring in rings:
            if len(ring) >= 12 and c_idx in ring and o_idx in ring:
                return True
    return False


def classify(smiles_series):
    """
    Given a Series of SMILES strings, return a DataFrame of boolean columns,
    one per antimicrobial family.
    """
    results = {family: [] for family in FAMILY_SMARTS}

    for smi in smiles_series:
        mol = Chem.MolFromSmiles(str(smi)) if pd.notna(smi) else None

        for family, pat in _COMPILED.items():
            if mol is None:
                results[family].append(False)
            elif family == "fluoroquinolone":
                results[family].append(
                    mol.HasSubstructMatch(pat) and mol.HasSubstructMatch(_F_SMARTS)
                )
            else:
                results[family].append(mol.HasSubstructMatch(pat))

        results["macrolide"].append(is_macrolide(mol))

    return pd.DataFrame(results)


# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------

print("Loading reference library...")
ref_df = pd.read_csv(os.path.join(data_dir, "reference_library_smiles.csv"))

print("Loading DrugBank...")
db_df = pd.read_csv(os.path.join(data_dir, "drugbank_smiles.csv"))

# ---------------------------------------------------------------------------
# Classify
# ---------------------------------------------------------------------------

print(f"\nClassifying {len(ref_df)} reference library compounds...")
ref_families = classify(ref_df["input"])
ref_out = pd.concat([ref_df, ref_families], axis=1)
ref_out_path = os.path.join(output_dir, "reference_library_families.csv")
ref_out.to_csv(ref_out_path, index=False)
print(f"  Saved to {ref_out_path}")

print(f"\nClassifying {len(db_df)} DrugBank compounds...")
db_families = classify(db_df["smiles"])
db_out = pd.concat([db_df, db_families], axis=1)
db_out_path = os.path.join(output_dir, "drugbank_families.csv")
db_out.to_csv(db_out_path, index=False)
print(f"  Saved to {db_out_path}")

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

families = list(FAMILY_SMARTS.keys())
summary = pd.DataFrame(
    {
        "family": families,
        "reference_library": [int(ref_families[f].sum()) for f in families],
        "drugbank": [int(db_families[f].sum()) for f in families],
    }
)
print("\n--- Compound family counts ---")
print(summary.to_string(index=False))
