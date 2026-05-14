"""Download eos74km (Antimicrobial class specificity prediction) outputs for the
Ersilia reference library and plot the distribution of each output.

Data sources
------------
- Reference library: downloaded from ersilia-model-hub-maintained-inputs via src/utils.py
- Model outputs: fetched from Isaura public bucket via src/isaura_utils.py
  Model: eos74km  Version: v1.0.1  Bucket: isaura-public

Output
------
    data/raw/compounds/reference_library_smiles.csv
    output/01_admet_properties/eos74km_v1.0.1.csv
    output/01_admet_properties/eos74km_distributions.png
"""

import math
import os
import sys

import pandas as pd
import stylia

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))

from isaura_utils import download_from_isaura
from utils import download_reference_library

# --- Paths ---
data_dir = os.path.join(root, "..", "data", "raw", "compounds")
output_dir = os.path.join(root, "..", "output", "01_admet_properties")
os.makedirs(data_dir, exist_ok=True)
os.makedirs(output_dir, exist_ok=True)

MODEL_ID = "eos74km"
MODEL_VERSION = "v1"
reference_csv = os.path.join(data_dir, "reference_library_smiles.csv")
output_csv = os.path.join(output_dir, f"{MODEL_ID}_{MODEL_VERSION}.csv")

# --- Download reference library if needed ---
if not os.path.exists(reference_csv):
    print("Downloading reference library...")
    df_ref = download_reference_library()
    df_ref.to_csv(reference_csv, index=False)
    print(f"  -> {reference_csv} ({len(df_ref)} compounds)")

# --- Download model outputs from Isaura if needed ---
if not os.path.exists(output_csv):
    print(f"Downloading {MODEL_ID} outputs from Isaura...")
    download_from_isaura(
        model_id=MODEL_ID,
        model_version=MODEL_VERSION,
        input_csv=reference_csv,
        output_path=output_dir,
    )
    print(f"  -> {output_csv}")

# --- Load outputs ---
df = pd.read_csv(output_csv)
print(f"Loaded {len(df)} rows, columns: {df.columns.tolist()}")

# Output columns: all numeric columns except the input key
output_cols = [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c])]
print(f"Plotting distributions for: {output_cols}")

# --- Plot ---
# Format: print | Style: article — change with stylia.set_format() / stylia.set_style()
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

stylia.save_figure(os.path.join(output_dir, "eos74km_distributions.png"))
print(f"  -> {os.path.join(output_dir, 'eos74km_distributions.png')}")
