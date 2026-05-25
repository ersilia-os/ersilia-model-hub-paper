"""Plot output distributions for eos7d58 (ADMET properties) across the reference library.

Requires data/processed/annotation_preds_ref_library/eos7d58_v1.0.0.csv
(run 00_download_data.py first).

Output
------
    output/xx_admet_properties/eos7d58_distributions.png
"""

import math
import os
import sys

import pandas as pd
import stylia

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))

annotation_dir = os.path.join(root, "..", "data", "processed", "annotation_preds_ref_library")
output_dir = os.path.join(root, "..", "output", "xx_admet_properties")
os.makedirs(output_dir, exist_ok=True)

MODEL_ID = "eos7d58"
MODEL_VERSION = "v1.0.0"
input_csv = os.path.join(annotation_dir, f"{MODEL_ID}_{MODEL_VERSION}.csv")

if not os.path.exists(input_csv):
    raise FileNotFoundError(
        f"{input_csv} not found — run 00_download_data.py first."
    )

# --- Load outputs ---
df = pd.read_csv(input_csv)
print(f"Loaded {len(df)} rows, columns: {df.columns.tolist()}")

output_cols = [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c])]
print(f"Plotting distributions for: {output_cols}")

# --- Plot ---
stylia.set_format("print")
stylia.set_style("article")

nc = stylia.NamedColors()

ncols = min(len(output_cols), 3)
nrows = math.ceil(len(output_cols) / ncols)


def plot_histogram(ax, values, col_name):
    ax.hist(values.dropna(), bins=40, color=nc.cobalt)
    stylia.label(ax, xlabel=col_name, ylabel=" ")


fig, axs = stylia.create_figure(nrows, ncols, width=0.5)
for col in output_cols:
    plot_histogram(axs.next(), df[col], col)

fig_path = os.path.join(output_dir, f"{MODEL_ID}_distributions.png")
stylia.save_figure(fig_path)
print(f"  -> {fig_path}")
