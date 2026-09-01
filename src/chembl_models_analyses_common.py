"""Shared, library-independent helpers for the xx_chembl_models_*.py
family of scripts (currently: _drugbank.py, _reflib.py) — the applicability-
domain check for the ChEMBL/PubChem pathogen models, run against different
compound libraries. Anything here is a property of the trained models
themselves, not of whichever library is being scored against them: the
consensus-threshold reconstruction, SMILES validity filtering, the
cache-merge pattern, and the 3-panel comparison figure. Per-library data
loading and the nearest-neighbour similarity search are NOT here — they
differ enough (data source shape, similarity-search architecture) between
libraries that sharing them would cost more in indirection than it saves.
"""

import os

import numpy as np
import pandas as pd
import stylia
from rdkit import Chem

W_COLS = ["w1", "w2", "w3", "w4", "w5", "w6"]


def require(path, hint):
    if not os.path.exists(path):
        raise FileNotFoundError(f"{path} not found.\n  {hint}")
    return path


def valid_smiles_idx(smiles_list):
    """Indices of smiles_list RDKit can parse."""
    return [i for i, s in enumerate(smiles_list) if Chem.MolFromSmiles(s) is not None]


def consensus_threshold(endpoints, k_star):
    """The published consensus_score threshold, reconstructed rather than fetched:
    apply consensus.py's W1-W6 (quality) + W7 (per-compound cutoff ramp) weighting
    formula to every sub-model's own decision_cutoff_rank, at the boundary where
    every sub-model's prob_rank sits exactly at its own cutoff (so W7=0 for all of
    them by construction), then the pathogen's tanh transform with its fixed k_star
    from 12b_k_star.json. Verified against the live eos21dr model: reproduces its
    published abaumannii consensus_score threshold (0.846) to within rounding
    (0.8459 computed). Library-independent: a property of the trained models.
    """
    cutoffs = endpoints["decision_cutoff_rank"].to_numpy(dtype=float)
    w_quality = endpoints[W_COLS].to_numpy(dtype=float)
    w7 = np.zeros(len(endpoints))
    w_eff = np.column_stack([w_quality, w7]).mean(axis=1)
    raw = (cutoffs * w_eff).sum() / w_eff.sum()
    return float(0.5 + 0.5 * np.tanh(k_star * (raw - 0.5)) / np.tanh(k_star / 2))


def merge_into_combined_summary(new_results, combined_path):
    """Merge new_results into whatever's already on disk at combined_path rather
    than overwriting it: rows for pathogens in new_results replace their own
    prior rows (in case of a re-run); every other pathogen's previously saved
    rows are kept untouched — so pathogens can be run one at a time,
    incrementally, without earlier runs' results disappearing."""
    if os.path.exists(combined_path):
        existing = pd.read_csv(combined_path)
        existing = existing[~existing["pathogen"].isin(new_results["pathogen"].unique())]
        combined = pd.concat([existing, new_results], ignore_index=True)
    else:
        combined = new_results
    combined.to_csv(combined_path, index=False)
    return combined


# ---------------------------------------------------------------------------
# Figure: 1 row x 3 columns (panel A sized to match B+C combined).
#   A: jittered strip by pathogen, one point per endpoint (sub-models + each
#      pathogen's consensus, where available), y = % of the library above
#      threshold.
#   B: % of the library in applicability domain (y, shared with C) vs
#      corr(rank, similarity to training actives) (x).
#   C: same shared y vs % of the library above threshold (x).
# Colors are pathogen-dependent and consistent across all three panels.
# ---------------------------------------------------------------------------

def plot_strip(ax, panel1_df, pathogen_color, all_pathogens, rng, ylabel):
    present = [p for p in all_pathogens if p in set(panel1_df["pathogen"])]
    for i, pathogen in enumerate(present):
        sub = panel1_df[panel1_df["pathogen"] == pathogen]
        jitter = rng.uniform(-0.15, 0.15, size=len(sub))
        is_consensus = (sub["endpoint"] == "consensus_score").to_numpy()
        x = np.full(len(sub), i) + jitter
        color = pathogen_color[pathogen]
        ax.scatter(x[~is_consensus], sub["pct_above"][~is_consensus], color=color, s=20, marker="o")
        ax.scatter(
            x[is_consensus], sub["pct_above"][is_consensus],
            color=color, s=35, marker="*", linewidth=0.5, edgecolor="black", zorder=3,
        )
    ax.set_xticks(range(len(present)))
    ax.set_xticklabels(present, rotation=90, ha="center")
    stylia.label(ax, xlabel="", ylabel=ylabel)

    dot_handle = ax.scatter([], [], color="black", s=20, marker="o", label="Model score")
    star_handle = ax.scatter(
        [], [], color="black", s=35, marker="*", linewidth=0.5, edgecolor="black", label="Consensus score"
    )
    ax.legend(handles=[dot_handle, star_handle], handletextpad=0.1)


def plot_indomain_vs_corr(ax, df, pathogen_color, ylim, in_domain_col, ylabel):
    ax.scatter(df["corr_rank_vs_sim_to_actives"], df[in_domain_col],
               color=df["pathogen"].map(pathogen_color).tolist())
    ax.axvline(0, color=stylia.NamedColors().silver, linewidth=0.5)
    ax.set_ylim(ylim)
    stylia.label(ax, xlabel="Correlation between rank\nand sim. to actives", ylabel=ylabel)


def plot_indomain_vs_hitrate(ax, df, pathogen_color, ylim, above_col, in_domain_col, xlabel, ylabel):
    ax.scatter(df[above_col], df[in_domain_col],
               color=df["pathogen"].map(pathogen_color).tolist())
    ax.set_ylim(ylim)
    stylia.label(ax, xlabel=xlabel, ylabel=ylabel)
