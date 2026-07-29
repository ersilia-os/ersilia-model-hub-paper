"""Reproduce the ChEMBL data-curation figures as Nature-ready panels.

These figures re-tell the curation story produced upstream by
``chembl-antimicrobial-tasks/scripts/27_general_plots.py`` (step 27), but are rebuilt here
purely from the small per-pathogen and aggregate SUMMARY CSVs staged by 00_download_data.py
into data/raw/chembl_curation/ — no full molecule-level datasets are copied. Figures whose
values are not present in those summaries (chemical-space overlap/coverage, embedding
scatters, molecule-level conflict panels) are intentionally skipped; save_curation_figures
prints which ones.

ChEMBL snapshot: chembl_36 (recorded in data/raw/chembl_curation/general/27_chembl_space.json).

Each panel is an individual figure saved as PNG + vector PDF on the 3 cm cell grid, with a
figure_cells.json footprint manifest. No A/B/C panel letters (ordering happens in Illustrator).

Output
------
    output/02_chembl_data_curation/png/<panel>.png
    output/02_chembl_data_curation/pdf/<panel>.pdf
    output/02_chembl_data_curation/figure_cells.json
"""

import json
import os
import sys

import pandas as pd

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))

from plots_chembl_curation import save_curation_figures  # noqa: E402
from default import RANDOM_SEED  # noqa: E402

config_dir = os.path.join(root, "..", "config")
data_dir = os.path.join(root, "..", "data", "raw", "chembl_curation")
output_dir = os.path.join(root, "..", "output", "02_chembl_data_curation")
os.makedirs(output_dir, exist_ok=True)

# Pathogen code -> full binomial, so the figures show abbreviated genus names ("A. baumannii")
# consistently with the other steps instead of the raw squashed codes.
pathogen_names = (
    pd.read_csv(os.path.join(config_dir, "pathogens_of_interest.csv"))
    .set_index("code")["pathogen"]
    .to_dict()
)

# Provenance: record the ChEMBL snapshot the curation ran against.
space_path = os.path.join(data_dir, "general", "27_chembl_space.json")
if os.path.exists(space_path):
    with open(space_path) as f:
        space = json.load(f)
    print(f"ChEMBL snapshot: {space.get('chembl_db')} "
          f"({space.get('bioactive_compounds')} bioactive compounds)")

save_curation_figures(data_dir=data_dir, output_dir=output_dir, random_seed=RANDOM_SEED,
                      name_map=pathogen_names)
