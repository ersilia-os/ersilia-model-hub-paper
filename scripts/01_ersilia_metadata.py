"""Analyse Ersilia Model Hub metadata and plot summary statistics.

Requires data/raw/airtable_metadata.csv (run 00_download_data.py first).

Output
------
    output/01_models_metadata/*_counts.csv
    output/01_models_metadata/*.png
"""

import os
import sys
import numpy as np
import pandas as pd
import circlify
import stylia
import matplotlib.patches as mpatches

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))

outpath = os.path.join(root, "..", "output", "01_models_metadata")
os.makedirs(outpath, exist_ok=True)

df = pd.read_csv(os.path.join(root, "..", "data", "raw", "airtable_metadata.csv"), encoding="utf-8-sig")
print(len(df))
df = df[df["Status"] == "Ready"]
print(len(df))

# Fields with comma-separated multiple values are exploded before counting
multi_value_fields = {"Output", "Biomedical Area", "Target Organism"}
fields = ["Output", "License", "Task", "Subtask", "Source Type", "Biomedical Area", "Target Organism", "Publication Type"]

print(f"Total models: {len(df)}\n")
for field in fields:
    if field in multi_value_fields:
        values = df[field].dropna().str.split(",").explode().str.strip()
    else:
        values = df[field].dropna()
    counts = values.value_counts().reset_index()
    counts.columns = ["value", "count"]
    fname = field.lower().replace(" ", "_") + "_counts.csv"
    counts.to_csv(os.path.join(outpath, fname), index=False)

# --- Plots ---
stylia.set_format("print")
stylia.set_style("article")

nc = stylia.NamedColors()

# Target Organism has 59 categories; cap at top 15 for readability
top_n = {"Target Organism": 15}

# Task→colour mapping; subtask inherits from its parent task
task_color = {
    "Annotation":     nc.cobalt,
    "Representation": nc.turquoise,
    "Sampling":       nc.tangerine,
}
subtask_parent = {
    "Activity prediction":                "Annotation",
    "Property calculation or prediction": "Annotation",
    "Featurization":                      "Representation",
    "Projection":                         "Representation",
    "Similarity search":                  "Sampling",
    "Generation":                         "Sampling",
}


def plot_bar(ax, counts, title, colors=None, n=None):
    if n:
        counts = counts.head(n)
        title = f"{title} (top {n})"
    c = colors if colors is not None else [nc.cobalt] * len(counts)
    ax.barh(counts["value"][::-1], counts["count"][::-1], color=c[::-1])
    stylia.label(ax, xlabel="Number of models", ylabel=" ", title=title)


# Task + Subtask: single figure with linked colours
task_counts = pd.read_csv(os.path.join(outpath, "task_counts.csv"))
subtask_counts = pd.read_csv(os.path.join(outpath, "subtask_counts.csv"))

# Reorder subtasks grouped by parent task (matching task bar order), then by count within each group
subtask_counts["parent"] = subtask_counts["value"].map(subtask_parent)
task_order = {t: i for i, t in enumerate(task_counts["value"])}
subtask_counts = (subtask_counts
                  .assign(task_rank=subtask_counts["parent"].map(task_order))
                  .sort_values(["task_rank", "count"], ascending=[True, False])
                  .drop(columns=["parent", "task_rank"])
                  .reset_index(drop=True))

fig, axs = stylia.create_figure(1, 2, width=1)
plot_bar(axs.next(), task_counts, "Task",
         colors=[task_color[v] for v in task_counts["value"]])
plot_bar(axs.next(), subtask_counts, "Subtask",
         colors=[task_color[subtask_parent[v]] for v in subtask_counts["value"]])
stylia.save_figure(os.path.join(outpath, "task_subtask.png"))

# All other fields — individual figures
other_fields = [f for f in fields if f not in {"Task", "Subtask"}]
for field in other_fields:
    fname = field.lower().replace(" ", "_") + "_counts.csv"
    counts = pd.read_csv(os.path.join(outpath, fname))
    fig, axs = stylia.create_figure(1, 1)
    plot_bar(axs.next(), counts, field, n=top_n.get(field))
    stylia.save_figure(os.path.join(outpath, field.lower().replace(" ", "_") + ".png"))


# --- Circle treemap: models per priority pathogen ---
pathogens_path = os.path.join(root, "..", "data", "raw", "pathogens_of_interest.csv")

source_color = {
    "External":   nc.cobalt,
    "Internal":   nc.turquoise,
    "Replicated": nc.tangerine,
}


def abbrev(name):
    """Abbreviate genus: 'Mycobacterium tuberculosis' → 'M. tuberculosis'."""
    parts = name.split()
    return f"{parts[0][0]}. {' '.join(parts[1:])}" if len(parts) > 1 else name


def sunflower(n, cx, cy, r, pad=0.18):
    """Place n dots inside a circle using the golden-angle sunflower spiral."""
    golden = np.pi * (3 - np.sqrt(5))
    rmax = r * (1 - pad)
    xs = [cx + rmax * np.sqrt((i + 0.5) / n) * np.cos(i * golden) for i in range(n)]
    ys = [cy + rmax * np.sqrt((i + 0.5) / n) * np.sin(i * golden) for i in range(n)]
    return xs, ys


def plot_circle_treemap(ax, df, pathogens_path):
    pathogens = pd.read_csv(pathogens_path)

    # Match each pathogen to models via multi-value Target Organism field
    df2 = df.copy()
    df2["orgs"] = df2["Target Organism"].fillna("").str.split(",").apply(
        lambda xs: [x.strip() for x in xs])

    rows = []
    for _, p in pathogens.iterrows():
        mask = df2["orgs"].apply(
            lambda orgs: any(p["pathogen"].lower() in o.lower() for o in orgs))
        rows.append((p["pathogen"], df2[mask].copy()))

    # Sort descending by count; use minimum area 0.3 so empty circles stay visible
    rows.sort(key=lambda x: len(x[1]), reverse=True)
    values = [max(len(models), 0.3) for _, models in rows]

    circles = circlify.circlify(values, show_enclosure=False)
    # Match largest circle to largest count
    circles_sorted = sorted(circles, key=lambda c: c.r, reverse=True)

    # Global dot-size scaling across all models (sqrt of output dimension)
    sqrt_dim = np.sqrt(df["Output Dimension"].fillna(1).clip(lower=1))
    lo, hi = sqrt_dim.min(), sqrt_dim.max()

    for circ, (name, models) in zip(circles_sorted, rows):
        cx, cy, r = circ.x, circ.y, circ.r
        n = len(models)

        facecolor = "#e8e8e8" if n == 0 else "#f0f0f0"
        ax.add_patch(mpatches.Circle((cx, cy), r,
                                      facecolor=facecolor, edgecolor="white", zorder=1))
        if n > 0:
            xs_d, ys_d = sunflower(n, cx, cy, r)
            sdim = np.sqrt(models["Output Dimension"].fillna(1).clip(lower=1))
            s = 15 + 150 * (sdim - lo) / (hi - lo + 1e-9)
            c = [source_color.get(v, nc.cobalt) for v in models["Source Type"].fillna("External")]
            ax.scatter(xs_d, ys_d, s=s.values, c=c, zorder=3)

        ax.text(cx, cy - r - 0.04, abbrev(name),
                ha="center", va="top", fontsize=stylia.FONTSIZE_SMALL, zorder=4)

    ax.set_xlim(-1.25, 1.25)
    ax.set_ylim(-1.25, 1.25)
    ax.set_aspect("equal")
    ax.set_axis_off()

    legend_handles = [mpatches.Patch(color=c, label=lbl) for lbl, c in source_color.items()]
    ax.legend(handles=legend_handles, loc="lower right")
    stylia.label(ax, title="Models targeting priority pathogens")


fig, axs = stylia.create_figure(1, 1, width=0.5, height=0.5)
plot_circle_treemap(axs.next(), df, pathogens_path)
stylia.save_figure(os.path.join(outpath, "pathogen_circles.png"))
