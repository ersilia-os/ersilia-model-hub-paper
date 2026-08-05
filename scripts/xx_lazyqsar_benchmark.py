"""Plot LazyQSAR v3.4.2 benchmark results (AUROC / AUPR vs Polaris LeaderBoard Best).

Reproduces, unchanged, the two-panel plot built in merge_plot_results.py of the
ersilia-ml-benchmark companion repo (cloned at the same level as this repo).

Requires data/raw/lazyqsar_benchmark/all_results.csv (run 00_download_data.py first)
and the ersilia-ml-benchmark companion repo checked out alongside this one, for its
src/common/defaults.py (benchmarks reference scores + official metric per benchmark).

Output
------
    output/xx_lazyqsar_benchmark/Figure_SX.png
    output/xx_lazyqsar_benchmark/Figure_SX.pdf
"""

import os
import sys

import numpy as np
import pandas as pd
import stylia
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

root = os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.abspath(os.path.join(root, ".."))
data_csv = os.path.join(repo_root, "data", "raw", "lazyqsar_benchmark", "all_results.csv")
output_dir = os.path.join(repo_root, "output", "xx_lazyqsar_benchmark")
os.makedirs(output_dir, exist_ok=True)

sys.path.append(os.path.join(repo_root, "..", "ersilia-ml-benchmark", "src", "common"))
from defaults import benchmarks, metrics

if not os.path.exists(data_csv):
    raise FileNotFoundError(f"{data_csv} not found — run 00_download_data.py first.")

df = pd.read_csv(data_csv)

MODE_LABELS = {"none": "RF (ECFP)", "fast": "LazyQSAR (fast)", "slow": "LazyQSAR (slow)"}

# --- Plot ---------------------------------------------------------------
# Two stacked panels: AUROC on top, AUPR at the bottom. A benchmark only
# gets a Polaris LeaderBoard Best reference line on the panel matching its
# official metric (roc_auc XOR pr_auc, per defaults.metrics), drawn as a
# horizontal black line spanning that benchmark's bar group. The AUPR panel
# additionally shows, per benchmark, the "expected" AUPR of a no-skill
# classifier (the positive-class prevalence, i.e. ratio_test) as a dashed
# black line.
benchmark_names = sorted(benchmarks.keys())
short_names = [b.split("/")[1].replace("-", "\n") for b in benchmark_names]
modes = ["none", "fast", "slow"]
PANELS = [
    ("auroc", "roc_auc", "AUROC", (0.5, 1.0), 0.1),
    ("aupr", "pr_auc", "AUPR", (0.1, 1.0), 0.2),
]
expected_aupr = {b: df.loc[df["benchmark"] == b, "ratio_test"].iloc[0] for b in benchmark_names}

stylia.set_format("print")
stylia.set_style("article")
colors = stylia.ArticleColors()
mode_colors = {"none": colors.silver, "fast": colors.cobalt, "slow": colors.crimson}

fig, axs = stylia.create_figure(nrows=2, ncols=1, width=1.0, height=0.6)

n_modes = len(modes)
bar_width = 0.8 / n_modes
group_half_width = 0.8 / 2
x = np.arange(len(benchmark_names))

# Capture the real Axes objects once; re-indexing through stylia's
# AxisManager (e.g. axs[0] again later) re-applies default axis styling
# and would silently wipe out labels/ticks set below.
panel_axes = [axs[panel_i] for panel_i in range(len(PANELS))]

for panel_i, (score_col, metric_name, ylabel, ylim, ytick_step) in enumerate(PANELS):
    ax = panel_axes[panel_i]
    for i, mode in enumerate(modes):
        scores = [
            df.loc[(df["benchmark"] == b) & (df["lq_mode"] == mode), score_col].iloc[0]
            for b in benchmark_names
        ]
        offset = (i - (n_modes - 1) / 2) * bar_width
        ax.bar(x + offset, scores, width=bar_width, color=mode_colors[mode], zorder=3)

    for xi, b in zip(x, benchmark_names):
        if metrics[b] == metric_name:
            ax.hlines(
                benchmarks[b],
                xi - group_half_width,
                xi + group_half_width,
                color=colors.black,
                linewidth=stylia.LINEWIDTH_THICK,
                zorder=5,
            )

    if score_col == "aupr":
        for xi, b in zip(x, benchmark_names):
            ax.hlines(
                expected_aupr[b],
                xi - group_half_width,
                xi + group_half_width,
                color=colors.black,
                linestyle="--",
                linewidth=stylia.LINEWIDTH_THICK,
                zorder=5,
            )

    ax.set_xticks(x)
    if panel_i == len(PANELS) - 1:
        ax.set_xticklabels(short_names, rotation=0, ha="center")
    else:
        ax.set_xticklabels([])
    ax.set_ylim(*ylim)
    ytick_start = np.ceil(ylim[0] / ytick_step) * ytick_step
    ax.set_yticks(np.round(np.arange(ytick_start, ylim[1] + 1e-9, ytick_step), 2))
    stylia.label(ax, xlabel="", ylabel=ylabel)

legend_handles = [Patch(color=mode_colors[mode], label=MODE_LABELS[mode]) for mode in modes]
legend_handles.append(Line2D([0], [0], color=colors.black, linewidth=stylia.LINEWIDTH_THICK, label="Polaris LeaderBoard Best"))
legend_handles.append(Line2D([0], [0], color=colors.black, linewidth=stylia.LINEWIDTH_THICK, linestyle="--", label="Expected AUPR (prevalence)"))
fig.suptitle("LazyQSAR v3.4.2 benchmark results", y=1.06, fontsize=stylia.FONTSIZE_BIG)
fig.legend(
    handles=legend_handles,
    loc="upper center",
    bbox_to_anchor=(0.5, 1.0),
    ncol=len(legend_handles),
    frameon=False,
)

stylia.save_figure(os.path.join(output_dir, "Figure_SX.png"))
stylia.save_figure(os.path.join(output_dir, "Figure_SX.pdf"))
