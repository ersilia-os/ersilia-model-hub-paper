"""Analyse Ersilia Model Hub metadata and plot summary statistics.

Requires data/raw/airtable_metadata.csv (run 00_download_data.py first),
config/pathogens_of_interest.csv and config/model_training_sizes.csv.

Counts are written per field; each metadata panel (Tasks & subtasks, Source Type,
Output, Biomedical Area, Target Organism, pathogen circle-treemap) is saved as its own
Nature-sized figure — raster PNG + vector PDF ready for Illustrator — using the reusable
plotting stack in src/ (plotting_base, plotting_colors, plots_metadata).

Source Type and Output are ALSO rendered as stacked bars segmented by Task and by Subtask
(four extra panels), so the joint distribution reads from one panel instead of two. The plain
single-colour versions are kept alongside them; which pair goes into the paper figure is a
layout choice made in Illustrator.

Output
------
    output/01_models_metadata/*_counts.csv                    # one per field
    output/01_models_metadata/*_by_{task,subtask}_counts.csv   # field x task cross-tabs
    output/01_models_metadata/license_grouped_counts.csv       # simplified licences + reuse class
    output/01_models_metadata/technical_metrics_summary.csv     # runtime + image size per task
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
from plotting_colors import LICENSE_CLASS, LICENSE_MISSING, SUBTASK_ORDER  # noqa: E402
from plotting_utils import abbrev  # noqa: E402
from default import (SUBTASK_PARENT, SUBTASK_DISPLAY, BIOAREA_DISPLAY,  # noqa: E402
                     RUNTIME_BATCH, RUNTIME_COLUMN)

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
multi_value_fields = {"Output", "Biomedical Area", "Target Organism", "Tag"}
fields = ["Output", "License", "Task", "Subtask", "Source Type",
          "Biomedical Area", "Target Organism", "Publication Type", "Tag",
          "Docker Architecture"]

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

# License, simplified for the figure. Two steps, both presentation-only:
#   1. Collapse the "-or-later" / "-only" suffixes, so GPL-3.0-or-later (51) and GPL-3.0-only (20)
#      become one GPL-3.0 bar (71). This drops a real legal distinction, which is why the ungrouped
#      license_counts.csv written above is kept as the record.
#   2. Label models with no license as LICENSE_MISSING rather than dropping them — unknown terms are
#      their own (worse) category for anyone reusing a model, not an absence of data to ignore.
# Ties are broken alphabetically so the four single-model licenses come out in a stable order.
lic = (df["License"].fillna(LICENSE_MISSING)
       .str.replace(r"-(or-later|only)$", "", regex=True))
lic_counts = lic.value_counts().reset_index()
lic_counts.columns = ["value", "count"]
lic_counts["class"] = lic_counts["value"].map(LICENSE_CLASS)
lic_counts = (lic_counts.sort_values(["count", "value"], ascending=[False, True])
              .reset_index(drop=True))
counts["License grouped"] = lic_counts
lic_counts.to_csv(os.path.join(outpath, "license_grouped_counts.csv"), index=False)

# Per-task summary of the two technical metrics behind the box panels, so the quartiles the figure
# shows are recoverable without re-reading the raw metadata. n_not_measured counts the Airtable -1
# sentinel (benchmark never run) — those rows are excluded from the statistics, never imputed, which
# is why the count is reported alongside.
tech_rows = []
for column, unit in [(RUNTIME_COLUMN, f"seconds for {RUNTIME_BATCH:,} molecules"),
                     ("Image Size", "MB")]:
    for task, g in df.groupby("Task"):
        v = g[column][g[column] > 0]
        tech_rows.append({
            "metric": column, "unit": unit, "task": task,
            "n_measured": len(v), "n_not_measured": int((g[column] <= 0).sum()),
            "median": v.median(), "q1": v.quantile(0.25), "q3": v.quantile(0.75),
            "min": v.min(), "max": v.max(),
        })
pd.DataFrame(tech_rows).to_csv(
    os.path.join(outpath, "technical_metrics_summary.csv"), index=False)

# Shorten crowded axis labels: Biomedical Area abbreviations and genus-abbreviated organisms.
counts["Biomedical Area"]["value"] = counts["Biomedical Area"]["value"].replace(BIOAREA_DISPLAY)
counts["Target Organism"]["value"] = counts["Target Organism"]["value"].map(abbrev)


# Field x Subtask model counts, feeding the stacked Source Type and Output panels. Rows keep the
# descending order the field's own counts use; columns follow SUBTASK_ORDER (grouped by parent task)
# so same-task shades stay adjacent in the stack. Each cross-tab is also written out as a CSV
# alongside the one-dimensional counts.
subtask_order = [SUBTASK_DISPLAY.get(s, s) for s in SUBTASK_ORDER]


def cross_tab(field):
    """Wide count table of ``field`` (rows) by Subtask (columns), in display order."""
    d = df[[field, "Subtask"]].copy()
    if field in multi_value_fields:
        d[field] = d[field].str.split(",")
        d = d.explode(field)
        d[field] = d[field].str.strip()
    d["Subtask"] = d["Subtask"].replace(SUBTASK_DISPLAY)
    table = pd.crosstab(d[field], d["Subtask"])
    rows = [v for v in counts[field]["value"] if v in table.index]
    cols = [c for c in subtask_order if c in table.columns]
    return table.loc[rows, cols]


cross_tabs = {}
for field in ["Source Type", "Output"]:
    table = cross_tab(field)
    cross_tabs[field] = table
    stem = f"{field.lower().replace(' ', '_')}_by_subtask_counts.csv"
    table.to_csv(os.path.join(outpath, stem))

# Render each panel as its own PNG + PDF figure.
save_metadata_figures(
    counts=counts,
    df=df,
    pathogens_path=pathogens_path,
    training_sizes_path=training_sizes_path,
    output_dir=outpath,
    top_n=top_n,
    cross_tabs=cross_tabs,
)
