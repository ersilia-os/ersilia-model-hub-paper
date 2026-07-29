"""Analyse Ersilia Model Hub metadata and plot summary statistics.

Requires data/raw/airtable_metadata.csv (run 00_download_data.py first),
config/pathogens_of_interest.csv and config/model_training_sizes.csv.

Counts are written per field; each metadata panel (Tasks & subtasks, Source Type,
Output, Biomedical Area, Target Organism, pathogen circle-treemap) is saved as its own
Nature-sized figure — raster PNG + vector PDF ready for Illustrator — using the reusable
plotting stack in src/ (plotting_base, plotting_colors, plots_metadata).

Output
------
    output/01_models_metadata/*_counts.csv
    output/01_models_metadata/png/<panel>.png
    output/01_models_metadata/pdf/<panel>.pdf
    output/01_models_metadata/figure_cells.json   # {panel: [rows, cols]} grid footprints
"""

import os
import sys
import pandas as pd

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))

from plots_metadata import save_metadata_figures  # noqa: E402
from plotting_utils import abbrev  # noqa: E402
from default import SUBTASK_PARENT, SUBTASK_DISPLAY, BIOAREA_DISPLAY  # noqa: E402

outpath = os.path.join(root, "..", "output", "01_models_metadata")
os.makedirs(outpath, exist_ok=True)

pathogens_path = os.path.join(root, "..", "config", "pathogens_of_interest.csv")
training_sizes_path = os.path.join(root, "..", "config", "model_training_sizes.csv")

df = pd.read_csv(os.path.join(root, "..", "data", "raw", "airtable_metadata.csv"),
                 encoding="utf-8-sig")
print(len(df))
df = df[df["Status"] == "Ready"]
print(len(df))

# Fields with comma-separated multiple values are exploded before counting.
multi_value_fields = {"Output", "Biomedical Area", "Target Organism"}
fields = ["Output", "License", "Task", "Subtask", "Source Type",
          "Biomedical Area", "Target Organism", "Publication Type"]

# Biomedical Area and Target Organism have many categories; cap at top 10 in the figure.
top_n = {"Biomedical Area": 10, "Target Organism": 10}

print(f"Total models: {len(df)}\n")

# Compute counts for every field, keep them in memory, and persist each as a CSV.
counts = {}
for field in fields:
    if field in multi_value_fields:
        values = df[field].dropna().str.split(",").explode().str.strip()
    else:
        values = df[field].dropna()
    c = values.value_counts().reset_index()
    c.columns = ["value", "count"]
    counts[field] = c
    fname = field.lower().replace(" ", "_") + "_counts.csv"
    c.to_csv(os.path.join(outpath, fname), index=False)

# Subtask: attach parent task, order by parent (matching Task order) then count within task,
# and apply the shorter display label. The parent column feeds the combined Task/Subtask panel.
task_order = {t: i for i, t in enumerate(counts["Task"]["value"])}
sub = counts["Subtask"].copy()
sub["parent"] = sub["value"].map(SUBTASK_PARENT)
sub["task_rank"] = sub["parent"].map(task_order)
sub = (sub.sort_values(["task_rank", "count"], ascending=[True, False])
       .drop(columns=["task_rank"])
       .reset_index(drop=True))
sub["value"] = sub["value"].replace(SUBTASK_DISPLAY)
counts["Subtask"] = sub

# Shorten crowded axis labels: Biomedical Area abbreviations and genus-abbreviated organisms.
counts["Biomedical Area"]["value"] = counts["Biomedical Area"]["value"].replace(BIOAREA_DISPLAY)
counts["Target Organism"]["value"] = counts["Target Organism"]["value"].map(abbrev)

# Render each panel as its own PNG + PDF figure.
save_metadata_figures(
    counts=counts,
    df=df,
    pathogens_path=pathogens_path,
    training_sizes_path=training_sizes_path,
    output_dir=outpath,
    top_n=top_n,
)
