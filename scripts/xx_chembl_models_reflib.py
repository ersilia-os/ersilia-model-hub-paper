"""Applicability-domain check for the ChEMBL/PubChem pathogen models — reference
library variant. Same question as xx_chembl_models_drugbank.py (is a
sub-model's hit rate driven by genuine chemical resemblance to its own training
actives, or extrapolation noise from a small/narrow training pool?), against the
Ersilia reference library (1,355,109 compounds,
data/raw/compound_lists/reference_library_smiles.csv) instead of DrugBank
(11,347) — in practice a random LIBRARY_SAMPLE_SIZE=250,000-compound sample of
it (see the constant's comment below for why: a full-library pilot run was
timed and found intractable, and library-side subsampling turns out to be
statistically unbiased here in a way training-pool subsampling would not be).
Two things are genuinely different here, not just "a bigger file":

1. DATA SOURCE. chembl-antimicrobial-models' internal pipeline has no
   reference-library predictions (only DrugBank/CoAdd), so this script reads
   Isaura's precalculated PACKAGED-model output instead — public column names
   (consensus_score, chembl_dose_response_0, ...), already rank-transformed
   (confirmed this session by reading lazyqsar's predict_rank() source), no
   pre-rank "raw" score. Bridging Isaura's public names to
   data/raw/chembl_model_reports/10_reports/10_reports.csv's internal names
   (DR_0001, SP_catchall, ...) needs each eos repo's own
   model/checkpoints/reports.csv (model_name=public, original_name=internal) —
   fetched once per pathogen, cached alongside the Isaura files.

2. SIMILARITY-SEARCH ARCHITECTURE. At 11,347 DrugBank queries, indexing each
   endpoint's (small) training pool and looping over the query library was
   fine. At 1.36M queries, doing that per endpoint would mean paying FPSim2's
   per-query fixed overhead well over a million times per endpoint. So this
   script INVERTS the loop: the reference library is indexed ONCE, globally
   (build_or_load_reference_index(), persisted to
   output/xx_chembl_models_reflib/.cache/reference_library_sample250000.h5), and
   queried with each endpoint's training compounds instead — Sigma(n_train)
   is ~3.05M across all 193 endpoints TOTAL (a fixed number, independent of
   library size), vs. up to 1.36M x 193 = 262M calls under the drugbank
   script's design. Each training-compound query returns similarity to every
   library compound; results are folded into a running elementwise max
   (max_sim_train, and max_sim_active from active training compounds only)
   rather than stored, to avoid materializing millions of 1.36M-length
   vectors. Parallelized across N_WORKERS processes, chunking the (per
   endpoint, usually far smaller than 1.36M) training compound list.
   Consequently there's no RDKit/FPSim2 size-based dispatch here (unlike the
   drugbank script) — the indexed side is always the full library, so FPSim2
   wins essentially everywhere.

Library-independent pieces (consensus-threshold reconstruction, SMILES
validity filtering, cache-merge pattern, the 3-panel figure) are shared with
the drugbank script via src/chembl_models_analyses_common.py — see that
module's docstring for why the rest isn't shared.

Requires:
  data/raw/compound_lists/reference_library_smiles.csv        (staged by 00_download_data.py)
  data/raw/chembl_model_reports/10_reports/10_reports.csv     (staged by 00_download_data.py)
  data/processed/annotation_preds_ref_library/{eosid}_{version}.csv
      — reference-library Isaura predictions for the 15 pathogen models.
        Run 00_download_data.py to fetch (routed through it per this repo's
        data-download convention — Isaura is category 4 there).
  ../chembl-antimicrobial-models/output/07_datasets/{pathogen}/{name}.csv
  ../chembl-antimicrobial-models/output/09_reports/{pathogen}/{name}_folds.json
      (also staged locally at data/raw/chembl_model_reports/{pathogen}/{name}_folds.json)
  ../chembl-antimicrobial-models/output/12_drugbank/12b_k_star.json  (for consensus thresholds)

Usage:
    python xx_chembl_models_reflib.py                              # all 15 pathogens (default)
    python xx_chembl_models_reflib.py --pathogen mtuberculosis
    python xx_chembl_models_reflib.py --pathogen mtuberculosis --endpoint 449762_merged3
        # single-endpoint pilot run — see the plan's mandatory-pilot step before a full run
    python xx_chembl_models_reflib.py --pathogen mtuberculosis --force  # redo despite cache
    python xx_chembl_models_reflib.py --compare  # also save the DrugBank-vs-reflib comparison figure
        # (needs xx_chembl_models_drugbank.py's applicability_domain_summary.csv to already exist)

Outputs:
    output/xx_chembl_models_reflib/.cache/reference_library_sample250000.h5   (persistent FPSim2 index, built once)
    output/xx_chembl_models_reflib/{pathogen}_applicability_domain.csv
    output/xx_chembl_models_reflib/applicability_domain_summary.csv   (all pathogens combined)
    output/xx_chembl_models_reflib/pathogen_endpoint_hit_rate.png / .pdf
    output/xx_chembl_models_reflib/drugbank_vs_reflib_comparison.png / .pdf   (only with --compare)
"""

import argparse
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
import requests
import stylia
from FPSim2 import FPSim2Engine
from FPSim2.io import create_db_file
from rdkit import RDLogger
from scipy import stats

RDLogger.DisableLog("rdApp.*")

root = os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.abspath(os.path.join(root, ".."))
sys.path.append(os.path.join(root, "..", "src"))
from default import RANDOM_SEED
from chembl_models_analyses_common import (
    require, valid_smiles_idx, consensus_threshold, merge_into_combined_summary,
    plot_strip, plot_indomain_vs_corr, plot_indomain_vs_hitrate,
)

CAM_ROOT = os.path.join(repo_root, "..", "chembl-antimicrobial-models")
LOCAL_FOLDS_DIR = os.path.join(repo_root, "data", "raw", "chembl_model_reports")
REPORTS_10_PATH = os.path.join(LOCAL_FOLDS_DIR, "10_reports", "10_reports.csv")
K_STAR_PATH = os.path.join(CAM_ROOT, "output", "12_drugbank", "12b_k_star.json")
REFLIB_DIR = os.path.join(repo_root, "data", "processed", "annotation_preds_ref_library")
REF_LIBRARY_SMILES_PATH = os.path.join(repo_root, "data", "raw", "compound_lists", "reference_library_smiles.csv")
PATHOGENS_CONFIG_PATH = os.path.join(repo_root, "config", "pathogens_of_interest.csv")

output_dir = os.path.join(repo_root, "output", "xx_chembl_models_reflib")
INDEX_DIR = os.path.join(output_dir, ".cache")
# LIBRARY_SAMPLE_SIZE: a full 1.36M-compound run was piloted on mtuberculosis's
# single largest endpoint (329,692 training compounds) and killed after 85+ min
# still running — extrapolated to the full Sigma(n_train)~=3.05M across all 193
# endpoints, that's many hours to a day+, not tractable. Fix: subsample the
# LIBRARY (the FPSim2-indexed side under this script's inverted architecture),
# not the training pool. This is statistically unbiased in a way that subsampling
# the training pool would not be: pct_library_above_threshold/pct_library_in_domain
# are proportions over the library, and random sampling gives an unbiased, tight
# estimate of them (~0.1-0.3% standard error at n=100k for the proportions in
# play here); each sampled library compound's max_sim_train/max_sim_active is
# still computed EXACTLY against the FULL, un-subsampled training set — no
# approximation on the nearest-neighbour side, so none of the extreme-value bias
# that subsampling the training pool itself would introduce. One fixed, seeded
# sample is drawn once and reused for every pathogen/endpoint, so results stay
# comparable across the whole analysis.
LIBRARY_SAMPLE_SIZE = 250_000
INDEX_PATH = os.path.join(INDEX_DIR, f"reference_library_sample{LIBRARY_SAMPLE_SIZE}.h5")
LIBRARY_SMILES_CACHE = os.path.join(INDEX_DIR, f"library_smiles_sample{LIBRARY_SAMPLE_SIZE}.csv")
os.makedirs(INDEX_DIR, exist_ok=True)

SIM_CUTOFF = 0.4
N_WORKERS = 8  # matches the drugbank script — this machine has 16 cores, 8 leaves headroom.

pathogens_df = pd.read_csv(PATHOGENS_CONFIG_PATH)
PATHOGEN_TO_EOSID = dict(zip(pathogens_df["code"], pathogens_df["eosid"]))
ALL_PATHOGENS = list(PATHOGEN_TO_EOSID.keys())

parser = argparse.ArgumentParser(
    description="Applicability-domain check for the ChEMBL/PubChem pathogen models' reference-library "
                 "predictions. See the top of this file for the full methodology."
)
parser.add_argument("--pathogen", default="all", help="Pathogen code, or 'all' for every pathogen (default: all).")
parser.add_argument("--endpoint", default=None,
                     help="Restrict to a single internal endpoint name (e.g. 449762_merged3) — for pilot runs.")
parser.add_argument("--force", action="store_true", help="Recompute even if a cached *_applicability_domain.csv exists.")
parser.add_argument("--compare", action="store_true",
                     help="Also save a DrugBank-vs-reference-library comparison figure "
                          "(needs xx_chembl_models_drugbank.py's summary to already exist).")
args = parser.parse_args()
pathogens = ALL_PATHOGENS if args.pathogen == "all" else [args.pathogen]


# ---------------------------------------------------------------------------
# One-time, persistent reference-library FPSim2 index
# ---------------------------------------------------------------------------

def build_or_load_reference_index():
    if os.path.exists(INDEX_PATH) and os.path.exists(LIBRARY_SMILES_CACHE):
        library_smiles = pd.read_csv(LIBRARY_SMILES_CACHE)["smiles"].tolist()
        print(f"Reference library index already built: {INDEX_PATH} ({len(library_smiles)} compounds)")
        return library_smiles

    print(f"Building reference library FPSim2 index (one-time; {LIBRARY_SAMPLE_SIZE}-compound "
          f"random sample, seed={RANDOM_SEED})...")
    t0 = time.time()
    ref_df = pd.read_csv(REF_LIBRARY_SMILES_PATH)
    all_smiles = ref_df["input"].tolist()
    valid_idx = valid_smiles_idx(all_smiles)
    print(f"  {len(valid_idx)}/{len(all_smiles)} valid ({time.time()-t0:.1f}s)")

    rng = np.random.default_rng(RANDOM_SEED)
    sample_idx = rng.choice(valid_idx, size=min(LIBRARY_SAMPLE_SIZE, len(valid_idx)), replace=False)
    library_smiles = [all_smiles[i] for i in sample_idx]
    print(f"  Sampled {len(library_smiles)}/{len(valid_idx)} valid library compounds "
          f"(random, seed={RANDOM_SEED})")

    create_db_file(
        mols_source=[(s, i) for i, s in enumerate(library_smiles)],
        filename=INDEX_PATH,
        mol_format="smiles",
        fp_type="Morgan",
        fp_params={"radius": 2, "fpSize": 2048},
    )
    pd.DataFrame({"smiles": library_smiles}).to_csv(LIBRARY_SMILES_CACHE, index=False)
    print(f"  Index built: {INDEX_PATH} ({time.time()-t0:.1f}s total)")
    return library_smiles


# ---------------------------------------------------------------------------
# Inverted nearest-neighbour search: index the library once (above), query
# with each endpoint's (much smaller) training compounds, chunked across
# N_WORKERS processes. Each worker accumulates its own running elementwise
# max over its chunk of training-compound queries; the main process combines
# the (<=N_WORKERS) partial arrays with a final elementwise max.
# ---------------------------------------------------------------------------

_worker_engine = None


def _init_reflib_worker(index_path):
    global _worker_engine
    _worker_engine = FPSim2Engine(index_path)


def _query_train_chunk(chunk):
    train_chunk_smiles, active_flags, n_library = chunk
    max_train = np.zeros(n_library)
    max_active = np.zeros(n_library)
    for smi, is_active in zip(train_chunk_smiles, active_flags):
        res = _worker_engine.similarity(smi, threshold=0.0, n_workers=1)
        if len(res) == 0:
            continue
        mol_id, coeff = res["mol_id"], res["coeff"]
        max_train[mol_id] = np.maximum(max_train[mol_id], coeff)
        if is_active:
            max_active[mol_id] = np.maximum(max_active[mol_id], coeff)
    return max_train, max_active


def nearest_neighbour_similarity_reflib(train_smiles, active_mask, n_library):
    has_actives = active_mask is not None and active_mask.any()
    active_flags = active_mask if active_mask is not None else np.zeros(len(train_smiles), dtype=bool)

    n_workers = min(N_WORKERS, len(train_smiles))
    idx_splits = [s for s in np.array_split(np.arange(len(train_smiles)), n_workers) if len(s)]
    chunks = [
        ([train_smiles[i] for i in idxs], active_flags[idxs], n_library)
        for idxs in idx_splits
    ]
    with ProcessPoolExecutor(
        max_workers=len(chunks), initializer=_init_reflib_worker, initargs=(INDEX_PATH,)
    ) as ex:
        results = list(ex.map(_query_train_chunk, chunks))

    max_train = np.maximum.reduce([r[0] for r in results])
    max_active = np.maximum.reduce([r[1] for r in results]) if has_actives else np.full(n_library, np.nan)
    return max_train, max_active


# ---------------------------------------------------------------------------
# Isaura data + internal<->public name mapping
# ---------------------------------------------------------------------------

def load_isaura_predictions(eosid):
    matches = [f for f in os.listdir(REFLIB_DIR) if f.startswith(f"{eosid}_v") and f.endswith(".csv")]
    if not matches:
        raise FileNotFoundError(
            f"No {eosid}_v*.csv in {REFLIB_DIR}.\n"
            f"  Run: python 00_download_data.py (fetches reference-library Isaura predictions)"
        )
    return pd.read_csv(os.path.join(REFLIB_DIR, sorted(matches)[-1]))


def fetch_internal_to_public_mapping(eosid):
    """internal name (DR_0001, ...) -> public name (chembl_dose_response_0, ...),
    from the live eos repo's own model/checkpoints/reports.csv. Cached alongside
    the Isaura predictions (same directory, skip-if-exists)."""
    cache_path = os.path.join(REFLIB_DIR, f"{eosid}_reports.csv")
    if not os.path.exists(cache_path):
        url = f"https://raw.githubusercontent.com/ersilia-os/{eosid}/main/model/checkpoints/reports.csv"
        resp = requests.get(url, timeout=60)
        resp.raise_for_status()
        with open(cache_path, "wb") as f:
            f.write(resp.content)
    reports = pd.read_csv(cache_path)
    return dict(zip(reports["original_name"], reports["model_name"]))


def compute_consensus_row_reflib(pathogen, endpoints, k_star_all, isaura_df, library_smiles):
    """Isaura already has consensus_score precomputed for multi-submodel pathogens
    (absent entirely for the 3 single-submodel ones) — no Tanimoto/FPSim2 work
    needed here, only the threshold (library-independent, reconstructed the same
    way as the drugbank script)."""
    if "consensus_score" not in isaura_df.columns:
        print(f"  No consensus for {pathogen} (single sub-model model — Isaura has no consensus_score column).")
        return None
    k_star_entry = k_star_all.get(pathogen)
    if not k_star_entry:
        print(f"  No consensus for {pathogen}: missing 12b_k_star.json entry.")
        return None

    thr = round(consensus_threshold(endpoints, k_star_entry["k_star"]), 4)
    cons = (
        pd.DataFrame({"smiles": library_smiles})
        .merge(isaura_df[["input", "consensus_score"]].rename(columns={"input": "smiles"}), on="smiles", how="left")
        ["consensus_score"].to_numpy()
    )
    n_missing = np.isnan(cons).sum()
    if n_missing:
        print(f"  NOTE: {n_missing}/{len(cons)} library compounds missing from Isaura consensus_score "
              f"— excluded from the hit-rate.")
    pct_above = 100 * np.nanmean(cons > thr)
    print(f"  consensus_score (M={k_star_entry['M']}, k_star={k_star_entry['k_star']}): "
          f"threshold={thr} above={pct_above:.2f}%")
    return {
        "pathogen": pathogen, "endpoint": "consensus_score", "threshold": thr,
        "n_train": np.nan, "n_train_active": np.nan, "pct_train_active": np.nan,
        "oof_active_mean": np.nan, "oof_inactive_mean": np.nan,
        "library_score_mean": round(float(np.nanmean(cons)), 4),
        "pct_library_above_threshold": round(pct_above, 4),
        "pct_library_in_domain": np.nan,
        "corr_rank_vs_sim_to_actives": np.nan,
        "similarity_trend_top_minus_bottom_decile": np.nan,
    }


# ---------------------------------------------------------------------------
# Per pathogen
# ---------------------------------------------------------------------------

def run_pathogen(pathogen, reports_10, k_star_all, library_smiles, force=False):
    eosid = PATHOGEN_TO_EOSID[pathogen]
    cache_path = os.path.join(output_dir, f"{pathogen}_applicability_domain.csv")
    endpoints = reports_10[reports_10["pathogen"] == pathogen]
    if args.endpoint:
        endpoints = endpoints[endpoints["model_name"] == args.endpoint]

    if os.path.exists(cache_path) and not force and not args.endpoint:
        print(f"\n{pathogen}: cached at {cache_path} — loading without recomputing (--force to redo).")
        return pd.read_csv(cache_path)

    print(f"\n{'='*70}\n{pathogen} ({eosid})\n{'='*70}")
    if endpoints.empty:
        print(f"  No retained endpoints for {pathogen} in 10_reports.csv — skipping.")
        return None

    isaura_df = load_isaura_predictions(eosid)
    internal_to_public = fetch_internal_to_public_mapping(eosid)
    n_library = len(library_smiles)

    rows = []
    for i, row in enumerate(endpoints.itertuples()):
        internal_name = row.model_name
        t0 = time.time()

        public_name = internal_to_public.get(internal_name)
        if public_name is None or public_name not in isaura_df.columns:
            print(f"  [{i+1}/{len(endpoints)}] SKIP {internal_name}: no public-name mapping or not in Isaura file")
            continue

        dataset_path = os.path.join(CAM_ROOT, "output", "07_datasets", pathogen, f"{internal_name}.csv")
        local_folds = os.path.join(LOCAL_FOLDS_DIR, pathogen, f"{internal_name}_folds.json")
        cam_folds = os.path.join(CAM_ROOT, "output", "09_reports", pathogen, f"{internal_name}_folds.json")
        folds_path = local_folds if os.path.exists(local_folds) else cam_folds
        if not (os.path.exists(dataset_path) and os.path.exists(folds_path)):
            print(f"  [{i+1}/{len(endpoints)}] SKIP {internal_name}: missing training set or OOF folds")
            continue

        train = pd.read_csv(dataset_path)
        with open(folds_path) as f:
            folds = json.load(f)
        oof_y, oof_yhat = [], []
        for fold in folds.values():
            oof_y.extend(fold["y_true"])
            oof_yhat.extend(fold["y_hat"])
        oof_y, oof_yhat = np.array(oof_y), np.array(oof_yhat)
        oof_active, oof_inactive = oof_yhat[oof_y == 1], oof_yhat[oof_y == 0]

        train_smiles_list = train["smiles"].tolist()
        train_valid = valid_smiles_idx(train_smiles_list)
        train_smiles = [train_smiles_list[i] for i in train_valid]
        train_labels = train["bin"].to_numpy()[train_valid]
        active_mask = train_labels == 1
        has_actives = active_mask.any()

        max_sim_train, max_sim_active = nearest_neighbour_similarity_reflib(train_smiles, active_mask, n_library)

        lib_rank = (
            pd.DataFrame({"smiles": library_smiles})
            .merge(isaura_df[["input", public_name]].rename(columns={"input": "smiles"}), on="smiles", how="left")
            [public_name].to_numpy()
        )

        threshold = round(float(row.decision_cutoff_rank), 4)
        pct_above = 100 * np.nanmean(lib_rank > threshold)
        pct_in_domain = 100 * (max_sim_train >= SIM_CUTOFF).mean()

        valid_pair = ~np.isnan(lib_rank)
        corr = (
            np.corrcoef(lib_rank[valid_pair], max_sim_active[valid_pair])[0, 1]
            if has_actives and valid_pair.sum() > 1 else np.nan
        )

        decile_trend = np.nan
        if has_actives:
            d = pd.DataFrame({"rank": lib_rank[valid_pair], "max_sim_active": max_sim_active[valid_pair]})
            d["decile"] = pd.qcut(d["rank"], 10, labels=False, duplicates="drop")
            decile_sim = d.groupby("decile")["max_sim_active"].mean()
            if len(decile_sim) > 1:
                decile_trend = decile_sim.iloc[-1] - decile_sim.iloc[0]

        rows.append({
            "pathogen": pathogen,
            "endpoint": internal_name,
            "threshold": threshold,
            "n_train": len(train),
            "n_train_active": int(train_labels.sum()),
            "pct_train_active": round(100 * train_labels.mean(), 2),
            "oof_active_mean": round(float(oof_active.mean()), 4) if len(oof_active) else np.nan,
            "oof_inactive_mean": round(float(oof_inactive.mean()), 4) if len(oof_inactive) else np.nan,
            "library_score_mean": round(float(np.nanmean(lib_rank)), 4),
            "pct_library_above_threshold": round(pct_above, 4),
            "pct_library_in_domain": round(pct_in_domain, 4),
            "corr_rank_vs_sim_to_actives": round(float(corr), 4) if not np.isnan(corr) else np.nan,
            "similarity_trend_top_minus_bottom_decile": (
                round(float(decile_trend), 4) if not np.isnan(decile_trend) else np.nan
            ),
        })
        print(f"  [{i+1}/{len(endpoints)}] {internal_name} (n_train={len(train)}): "
              f"above={pct_above:.2f}% in_domain={pct_in_domain:.2f}% corr={corr:.3f} "
              f"({time.time()-t0:.1f}s)")
        # Written after every endpoint, not just at the end: at this runtime scale a
        # multi-endpoint pathogen can run for hours, and losing all of it to an
        # interruption partway through would be expensive to redo.
        pd.DataFrame(rows).to_csv(cache_path, index=False)

    if not args.endpoint:
        cons_row = compute_consensus_row_reflib(pathogen, endpoints, k_star_all, isaura_df, library_smiles)
        if cons_row is not None:
            rows.append(cons_row)
            pd.DataFrame(rows).to_csv(cache_path, index=False)

    pathogen_df = pd.DataFrame(rows)
    print(f"  Saved: {cache_path}  ({len(pathogen_df)} endpoints)")
    return pathogen_df


# ---------------------------------------------------------------------------
# --compare: DrugBank-vs-reference-library comparison figure
# ---------------------------------------------------------------------------

def plot_drugbank_vs_reflib_comparison(reflib_df, drugbank_summary_path, pathogen_color, output_dir):
    drugbank_summary_path = require(
        drugbank_summary_path,
        "Run xx_chembl_models_drugbank.py first to generate the DrugBank-side summary.",
    )
    db = pd.read_csv(drugbank_summary_path)
    merged = db.merge(reflib_df, on=["pathogen", "endpoint"], suffixes=("_db", "_rl"))
    merged["avg_in_domain"] = merged[["pct_drugbank_in_domain", "pct_library_in_domain"]].mean(axis=1)
    nc = stylia.NamedColors()

    def plot_comparison(ax, xcol, ycol, xlabel, ylabel, size_col=None):
        cols = [xcol, ycol] + ([size_col] if size_col else [])
        df = merged.dropna(subset=cols)
        ax.plot([0, 100], [0, 100], linestyle="--", color=nc.silver, zorder=1, linewidth=stylia.LINEWIDTH)
        for pathogen in pathogen_color:
            sub = df[df["pathogen"] == pathogen]
            kwargs = {"s": 2 + sub[size_col] * 2} if size_col else {}
            ax.scatter(sub[xcol], sub[ycol], color=pathogen_color[pathogen], zorder=2, **kwargs)
        pear = stats.pearsonr(df[xcol], df[ycol])
        title = f"Pearson r = {pear.statistic:.3f}"
        if size_col:
            title += " (size = avg % in domain)"
            for v in [1, 10, 50]:
                ax.scatter([], [], s=2 + v * 2, color=nc.silver, label=f"{v}%")
            ax.legend(title="avg % in domain", loc="upper left")
        stylia.label(ax, xlabel=xlabel, ylabel=ylabel, title=title)

    fig, axs = stylia.create_figure(1, 3)
    plot_comparison(axs.next(), "pct_drugbank_above_threshold", "pct_library_above_threshold",
                     "% DrugBank above threshold", "% Reference library above threshold")
    plot_comparison(axs.next(), "pct_drugbank_in_domain", "pct_library_in_domain",
                     "% DrugBank in applicability domain", "% Reference library in applicability domain")
    plot_comparison(axs.next(), "pct_drugbank_above_threshold", "pct_library_above_threshold",
                     "% DrugBank above threshold", "% Reference library above threshold",
                     size_col="avg_in_domain")

    fig_path_png = os.path.join(output_dir, "drugbank_vs_reflib_comparison.png")
    fig_path_pdf = os.path.join(output_dir, "drugbank_vs_reflib_comparison.pdf")
    stylia.save_figure(fig_path_png)
    stylia.save_figure(fig_path_pdf)
    print(f"Saved: {fig_path_png} / .pdf")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

library_smiles = build_or_load_reference_index()

reports_10 = pd.read_csv(REPORTS_10_PATH)
if os.path.exists(K_STAR_PATH):
    with open(K_STAR_PATH) as f:
        k_star_all = json.load(f)
else:
    print(f"NOTE: {K_STAR_PATH} not found — no consensus rows will be computed for any pathogen.")
    k_star_all = {}

t_all = time.time()
results = [run_pathogen(p, reports_10, k_star_all, library_smiles, force=args.force) for p in pathogens]
results = [r for r in results if r is not None and len(r)]

if not results:
    print("\nNo results produced.")
    sys.exit(0)

new_results = pd.concat(results, ignore_index=True)
combined_path = os.path.join(output_dir, "applicability_domain_summary.csv")
combined = merge_into_combined_summary(new_results, combined_path)
print(f"\nTotal elapsed: {time.time()-t_all:.1f}s. Combined summary: {combined_path} "
      f"({len(combined)} endpoints across {combined['pathogen'].nunique()} pathogen(s) total; "
      f"{len(new_results)} from this run)")

if args.endpoint:
    print("Single-endpoint pilot run — skipping the combined figure.")
    sys.exit(0)

if combined["pathogen"].nunique() < 2:
    print("\nFewer than 2 pathogens in the combined summary — skipping the 3-panel figure.")
    sys.exit(0)

panel1_df = combined.rename(columns={"pct_library_above_threshold": "pct_above"})
bc_df = combined[combined["endpoint"] != "consensus_score"]

# Format: print | Style: article — change with stylia.set_format() / stylia.set_style()
stylia.set_format("print")
stylia.set_style("article")

rng = np.random.default_rng(RANDOM_SEED)
pal = stylia.CategoricalPalette("npg")
pathogen_color = dict(zip(ALL_PATHOGENS, pal.get(len(ALL_PATHOGENS))))

y_margin = 0.05 * (bc_df["pct_library_in_domain"].max() - bc_df["pct_library_in_domain"].min())
shared_ylim = (
    bc_df["pct_library_in_domain"].min() - y_margin,
    bc_df["pct_library_in_domain"].max() + y_margin,
)

fig, axs = stylia.create_figure(1, 3, width_ratios=[2, 1, 1], height=0.35)
plot_strip(axs.next(), panel1_df, pathogen_color, ALL_PATHOGENS, rng, ylabel="% reference library above threshold")
ax_b, ax_c = axs.next(), axs.next()
plot_indomain_vs_corr(ax_b, bc_df, pathogen_color, shared_ylim, "pct_library_in_domain",
                       ylabel="% reference library in applicability domain")
plot_indomain_vs_hitrate(ax_c, bc_df, pathogen_color, shared_ylim, "pct_library_above_threshold",
                          "pct_library_in_domain", xlabel="% reference library above threshold", ylabel="")
ax_c.sharey(ax_b)

fig_path_png = os.path.join(output_dir, "pathogen_endpoint_hit_rate.png")
fig_path_pdf = os.path.join(output_dir, "pathogen_endpoint_hit_rate.pdf")
stylia.save_figure(fig_path_png)
stylia.save_figure(fig_path_pdf)
print(f"Saved: {fig_path_png} / .pdf")

if args.compare:
    drugbank_summary_path = os.path.join(repo_root, "output", "xx_chembl_models_drugbank", "applicability_domain_summary.csv")
    plot_drugbank_vs_reflib_comparison(combined, drugbank_summary_path, pathogen_color, output_dir)
