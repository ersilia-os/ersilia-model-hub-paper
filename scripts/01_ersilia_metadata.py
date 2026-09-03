"""Analyse Ersilia Model Hub metadata and plot summary statistics.

Requires the frozen Airtable snapshot named by AIRTABLE_METADATA_FILE in src/default.py
(run 00_download_data.py first), config/pathogens_of_interest.csv and
config/model_training_sizes.csv.

Counts are written per field; each metadata panel (Tasks & subtasks, Biomedical Area,
Target Organism, the three donuts, the technical box row, pathogen circle-treemap) is saved as its
own Nature-sized figure — raster PNG + vector PDF ready for Illustrator — using the reusable
plotting stack in src/ (plotting_base, plotting_colors, plots_metadata).

Source Type and Output are ALSO rendered as stacked bars segmented by Subtask, so the joint
distribution reads from one panel instead of two.

Output
------
    output/01_models_metadata/*_counts.csv                       # one per field
    output/01_models_metadata/*_by_subtask_counts.csv            # field x subtask cross-tabs
    output/01_models_metadata/license_grouped_counts.csv         # simplified licences + reuse class
    output/01_models_metadata/technical_metrics_summary.csv      # runtime, image size, output dim per task
    output/01_models_metadata/output_dimension_bins.csv          # task x decade counts (circle panel)
    output/01_models_metadata/biomedical_area_groups.csv         # 4-group classification, per area
    output/01_models_metadata/models_over_time_by_task.csv       # cumulative series for step 01b
    output/01_models_metadata/png/<panel>.png
    output/01_models_metadata/pdf/<panel>.pdf
    output/01_models_metadata/figure_cells.json   # {panel: [rows, cols]} grid footprints
"""

import os
import sys
import numpy as np
import pandas as pd

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))

from plots_metadata import save_metadata_figures  # noqa: E402
from plotting_colors import (LICENSE_CLASS, LICENSE_MISSING, SUBTASK_ORDER,  # noqa: E402
                             TASK_COLORS)
from plotting_utils import abbrev  # noqa: E402
from default import (SUBTASK_PARENT, SUBTASK_DISPLAY, BIOAREA_DISPLAY,  # noqa: E402
                     BIOAREA_GROUP, BIOAREA_GROUP_OTHER, ACTIVITY_SUBTASK,
                     RUNTIME_BATCH, RUNTIME_COLUMN,
                     AIRTABLE_METADATA_FILE, AIRTABLE_SNAPSHOT_DATE)

outpath = os.path.join(root, "..", "output", "01_models_metadata")
os.makedirs(outpath, exist_ok=True)

pathogens_path = os.path.join(root, "..", "config", "pathogens_of_interest.csv")
training_sizes_path = os.path.join(root, "..", "config", "model_training_sizes.csv")

# `keep_default_na=False` matches 00_download_data.py's read of the same file: `None` is a
# legitimate License value ("repo checked, confirmed no LICENSE file"), and pandas' default NA set
# would otherwise coerce it to NaN, indistinguishable from a genuinely blank cell.
df_all = pd.read_csv(os.path.join(root, "..", "data", "raw", AIRTABLE_METADATA_FILE),
                     encoding="utf-8-sig", keep_default_na=False)
print(f"Metadata: {AIRTABLE_METADATA_FILE} "
      f"(last Airtable pull {AIRTABLE_SNAPSHOT_DATE})")
print(len(df_all))
# The cumulative series below needs the unfiltered frame (every model that carries a date, whatever
# its Status), so df_all is kept alongside the Ready-only frame the snapshot panels use.
df = df_all[df_all["Status"] == "Ready"]
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
#      their own (worse) category for anyone reusing a model, not an absence of data to ignore. The
#      literal Airtable value "None" (repo checked, confirmed no LICENSE file) is folded into the
#      same LICENSE_MISSING bucket as a genuinely blank cell — per project convention, Ersilia never
#      leaves the license check undone, so "None" and "not recorded" describe the same reuser-facing
#      outcome and are kept as one field rather than two.
# Ties are broken alphabetically so the four single-model licenses come out in a stable order.
lic = (df["License"].replace({"": LICENSE_MISSING, "None": LICENSE_MISSING})
       .str.replace(r"-(or-later|only)$", "", regex=True))
lic_counts = lic.value_counts().reset_index()
lic_counts.columns = ["value", "count"]
lic_counts["class"] = lic_counts["value"].map(LICENSE_CLASS)
# Same guard as BIOAREA_GROUP below, and for the same reason. LicenseClassDonutPlot groups by the
# class and drops NaN, so a licence with no LICENSE_CLASS entry would silently remove its models
# from the donut — a figure quietly stating a smaller n than the rest of the panel set. Added
# 2026-08-14, when the manual metadata revision introduced three unmapped values at once
# (Non-commercial, CC-BY-NC-SA-4.0, NCSA) covering 4 models.
_unmapped_lic = sorted(lic_counts.loc[lic_counts["class"].isna(), "value"])
if _unmapped_lic:
    raise KeyError(f"License values with no LICENSE_CLASS entry: {_unmapped_lic}. "
                   "Add them to src/plotting_colors.py — an unmapped licence must not vanish "
                   "from the donut.")
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
                     ("Image Size", "MB"),
                     ("Output Dimension", "values per prediction")]:
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

# Output Dimension is drawn as decade-binned circles rather than a box-and-swarm (the column is
# heavily tied — 68 of 133 Annotation models output a single value, and 102 of them fall in the 1-9
# bin), and a circle's area is the only thing carrying its count. This is the record of the exact
# numbers behind those nine circles.
od = df[df["Output Dimension"] > 0].copy()
od["decade"] = np.floor(np.log10(od["Output Dimension"])).astype(int)
od_bins = pd.crosstab(od["Task"], od["decade"])
od_bins.columns = [f"{10 ** c}-{10 ** (c + 1) - 1}" for c in od_bins.columns]
od_bins.to_csv(os.path.join(outpath, "output_dimension_bins.csv"))

# Biomedical Area collapsed into five groups (BIOAREA_GROUP, signed off 2026-08-02, Antifungal split
# out 2026-08-07, Leishmaniasis added 2026-08-14), over ACTIVITY
# PREDICTION models only. Two things this counting has to get right:
#   1. Biomedical Area is MULTI-VALUE, so counts are of DISTINCT MODELS per group, not of area
#      assignments. Grouping absorbs most of the multiplicity — 16 of the 22 multi-area Annotation
#      models have all their areas inside one group (AMR+Pneumonia, AMR+Diarrhoea, Gonorrhea+AMR).
#   2. Models that still span two groups are counted in BOTH, so the bars sum to more than n. That is
#      the metadata's own claim and is not silently resolved; the count is printed below.
act = df[df["Subtask"] == ACTIVITY_SUBTASK]
ba = act[["Identifier", "Biomedical Area"]].dropna(subset=["Biomedical Area"]).copy()
ba["Biomedical Area"] = ba["Biomedical Area"].str.split(",")
ba = ba.explode("Biomedical Area")
ba["Biomedical Area"] = ba["Biomedical Area"].str.strip()
ba["group"] = ba["Biomedical Area"].map(BIOAREA_GROUP)
_unmapped = sorted(ba.loc[ba["group"].isna(), "Biomedical Area"].unique())
if _unmapped:
    raise KeyError(f"Biomedical Area values with no BIOAREA_GROUP entry: {_unmapped}. "
                   "Add them to src/default.py — an unmapped area must not vanish from the figure.")

# Substantive groups by size descending, with the catch-all pinned last however large it grows.
_g = ba.groupby("group")["Identifier"].nunique().sort_values(ascending=False)
_order = [g for g in _g.index if g != BIOAREA_GROUP_OTHER] + \
         ([BIOAREA_GROUP_OTHER] if BIOAREA_GROUP_OTHER in _g.index else [])
counts["Biomedical Area grouped"] = pd.DataFrame(
    {"value": _order, "count": [int(_g[g]) for g in _order]})
# Distinct models behind the five groups, which is NOT their sum: five models carry areas in two
# different groups, so the counts add to 98 over 93 models. The donut shows this in its hole, and
# would otherwise state a model count that does not exist.
counts["Biomedical Area grouped"].attrs["n_models"] = int(ba["Identifier"].nunique())
(ba.groupby(["group", "Biomedical Area"])["Identifier"].nunique()
   .rename("n_models").reset_index()
   .to_csv(os.path.join(outpath, "biomedical_area_groups.csv"), index=False))

_spanning = ba.groupby("Identifier")["group"].nunique()
_spanning = sorted(_spanning[_spanning > 1].index)
print(f"\n[biomedical_area] {ACTIVITY_SUBTASK} models: {ba['Identifier'].nunique()} of {len(act)}; "
      + ", ".join(f"{g} {int(_g[g])}" for g in _order))
print(f"{'':18s}bars sum to {int(_g.sum())} because {len(_spanning)} model(s) span two groups "
      f"and are counted in both: {', '.join(_spanning)}")

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

# ---------------------------------------------------------------------------
# Cumulative models over time (summary CSV only — no panel in this script)
# ---------------------------------------------------------------------------
# Every panel above is a snapshot. This block adds no figure of its own: it is the only consumer of
# Airtable's date columns, and it exists to write the one month-indexed summary CSV that
# scripts/01b_community_stats.py reads for its Models track. Keep it in step 01 rather than moving
# it into 01b — that script reads this one series as a pre-aggregated summary, never raw Airtable.
#
# It runs on df_all, NOT the Ready-only frame: a model that is currently in maintenance was still
# incorporated on its date and still counts towards how the hub grew.
#
# Models with no Incorporation Date are excluded — the date is the axis, so there is nowhere to
# put them. They are never imputed, and the count is printed below for the caption.
dated = df_all[df_all["Incorporation Date"].notna()].copy()
dated["month"] = pd.to_datetime(dated["Incorporation Date"]).dt.to_period("M")
months = pd.period_range(dated["month"].min(), dated["month"].max(), freq="M")
# Month ends, so the last point sits at the real end of the series rather than on its first day.
month_ends = months.to_timestamp(how="end")

n_undated = int(df_all["Incorporation Date"].isna().sum())
print(f"\n[growth] {len(dated)} models carry an Incorporation Date "
      f"({dated['Incorporation Date'].min()} to {dated['Incorporation Date'].max()}); "
      f"{n_undated} have none and are excluded, not imputed")

# Reindexed onto the full month range before the cumulative sum, so months with no new model still
# carry a point and the series is flat there rather than interpolated across the gap. Columns follow
# TASK_COLORS so the task order matches every other panel in the figure.
by_task = pd.crosstab(dated["month"], dated["Task"]).reindex(months, fill_value=0)
by_task = by_task[[c for c in TASK_COLORS if c in by_task.columns]].cumsum()
by_task.index = month_ends
by_task.to_csv(os.path.join(outpath, "models_over_time_by_task.csv"))

print(f"[growth] cumulative total at {month_ends[-1].date()}: "
      f"{int(by_task.iloc[-1].sum())} models (all dated) "
      f"-> models_over_time_by_task.csv, read by scripts/01b_community_stats.py")
