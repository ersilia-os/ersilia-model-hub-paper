"""Applicability-domain check for the ChEMBL/PubChem pathogen models' DrugBank
predictions: is a sub-model's DrugBank "hit rate" (compounds scoring above its
recommended threshold) driven by genuine chemical resemblance to its own
training actives, or is it extrapolation noise from a small/narrow training pool?

For every retained sub-model of a pathogen (one row per endpoint in
data/raw/chembl_model_reports/10_reports/10_reports.csv), this script:
  1. Loads the endpoint's training compounds + labels and its per-fold
     out-of-fold (OOF) scores.
  2. Loads DrugBank's raw score and rank-transformed score on that endpoint.
  3. Computes, for every DrugBank compound, its nearest-neighbour Tanimoto
     similarity (Morgan r=2, 2048 bit) to the endpoint's training actives and
     to the full training pool (the applicability-domain check), via FPSim2
     (bit-packed/SIMD Tanimoto search — verified to reproduce RDKit's
     BulkTanimotoSimilarity bit-for-bit) with the query list split across
     N_WORKERS processes for endpoints large enough to benefit.
  4. Reports: % of DrugBank above the endpoint's recommended threshold, % of
     DrugBank inside a conventional applicability domain (max similarity to
     any training compound >= 0.4), and the correlation between predicted
     rank and similarity to training actives — a near-zero or negative
     correlation on a low-domain-coverage endpoint is the signature found for
     several pfalciparum/mtuberculosis endpoints in this repo's exploratory
     analysis (see conversation history for the original exploratory work).

A pathogen already computed (output/xx_chembl_models_drugbank/{pathogen}_
applicability_domain.csv exists) is loaded from that cache instead of
recomputed — pass --force to redo it anyway. This is what makes `--pathogen
all` cheap to re-run after the first full pass: only new/forced pathogens
pay the Tanimoto cost.

Per-endpoint per-molecule data (training sets, OOF folds, DrugBank score/rank)
is read directly from the companion chembl-antimicrobial-models repo, cloned
at the same level as this one — per this repo's convention, full per-molecule
datasets are never copied into ersilia-model-hub-paper's data/ tree, only the
small pre-aggregated 10_reports.csv summary (already staged by
00_download_data.py) and this script's own endpoint-level output.

Each pathogen's consensus_score (12 of 15 pathogens have one — the other 3
ship a single sub-model, so there is no consensus to compute) is not fetched
from anywhere: its recommended threshold is reconstructed locally from
10_reports.csv's per-sub-model W1-W6/decision_cutoff_rank plus the
companion repo's 12b_k_star.json, applying the same weighting+tanh formula
consensus.py itself uses (see consensus_threshold() below) — verified to
reproduce the published eos21dr/abaumannii threshold (0.846) to within
rounding. This is what lets panel A mark each pathogen's consensus with a
star without depending on anything outside this script + the companion repo.

Requires, in ../chembl-antimicrobial-models (pull with `eosvc download --path
<subpath>` there if missing):
    output/07_datasets/{pathogen}/{name}.csv
    output/09_reports/{pathogen}/{name}_folds.json   (also staged locally at
        data/raw/chembl_model_reports/{pathogen}/{name}_folds.json — used
        preferentially since it is already in this repo)
    output/12_drugbank/score/{pathogen}.csv
    output/12_drugbank/rank/{pathogen}.csv
    output/12_drugbank/12b_k_star.json                  (for consensus rows)
    output/14_consensus/{pathogen}_transformed.csv      (for consensus rows)

Nothing under this repo's tmp/ is required — this script is fully
self-contained given the companion repo + its own cached output/.

Usage:
    python xx_chembl_models_drugbank.py                            # all 15 pathogens (default); cached ones are free
    python xx_chembl_models_drugbank.py --pathogen mtuberculosis
    python xx_chembl_models_drugbank.py --pathogen mtuberculosis --force   # redo despite cache

Outputs:
    output/xx_chembl_models_drugbank/{pathogen}_applicability_domain.csv
    output/xx_chembl_models_drugbank/applicability_domain_summary.csv   (all pathogens combined)
    output/xx_chembl_models_drugbank/pathogen_endpoint_hit_rate.png / .pdf   (panels A/B/C, see below)
"""

import argparse
import json
import os
import sys
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
import stylia
from FPSim2 import FPSim2Engine
from FPSim2.io import create_db_file
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem

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
# Library-specific paths, centralized here (not spread across the script) so that
# pointing this whole script at a different compound library later is a matter of
# editing these few lines. (The reference-library variant, xx_chembl_models_
# reflib.py, needed enough more than a path swap — different data
# source shape, inverted similarity-search architecture — that it's a separate
# script rather than a parameterization of this one; see its own docstring.)
K_STAR_PATH = os.path.join(CAM_ROOT, "output", "12_drugbank", "12b_k_star.json")
SCORE_DIR = os.path.join(CAM_ROOT, "output", "12_drugbank", "score")
RANK_DIR = os.path.join(CAM_ROOT, "output", "12_drugbank", "rank")
CONSENSUS_DIR = os.path.join(CAM_ROOT, "output", "14_consensus")
output_dir = os.path.join(repo_root, "output", "xx_chembl_models_drugbank")
os.makedirs(output_dir, exist_ok=True)

SIM_CUTOFF = 0.4

# Nearest-neighbour Tanimoto search tuning (see nearest_neighbour_similarity()).
# FPSim2's bit-packed/SIMD search has more fixed per-query overhead than RDKit's
# BulkTanimotoSimilarity (HDF5 index open, per-query sanitization, structured-array
# construction), so it only pays for itself once the training pool is large enough
# that its per-comparison speed advantage outweighs that overhead. Benchmarked
# directly against RDKit on real endpoints (mtuberculosis, DrugBank-scale queries):
#   n_train=1,937   -> FPSim2 roughly on par with / slightly worse than RDKit
#   n_train=23,462  -> FPSim2 serial already 1.5x faster (34.0s -> 22.1s),
#                      ~5x faster with 8-way process parallelism (-> ~6.7s)
#   n_train=329,692 -> FPSim2 serial 2.2x faster (520s -> 276s incl. index build)
# LARGE_TRAIN_THRESHOLD is set well inside the confirmed-faster region so RDKit
# (already optimal for small pools) is never made slower by this change.
LARGE_TRAIN_THRESHOLD = 10_000
N_WORKERS = 8  # user-specified: this machine has 16 cores, 8 leaves headroom.
# Below this many (query x train) comparisons, skip the process pool even for a
# large-enough training pool — fork overhead (tens of ms x 8) would exceed the
# time saved on a small job.
PARALLEL_MIN_COMPARISONS = 20_000_000
ALL_PATHOGENS = [
    "abaumannii", "calbicans", "campylobacter", "ecoli", "efaecium", "enterobacter",
    "hpylori", "kpneumoniae", "mtuberculosis", "ngonorrhoeae", "paeruginosa",
    "pfalciparum", "saureus", "smansoni", "spneumoniae",
]

parser = argparse.ArgumentParser(
    description="Applicability-domain check for the ChEMBL/PubChem pathogen models' DrugBank "
                 "predictions. See the top of this file for the full methodology."
)
parser.add_argument("--pathogen", default="all", help="Pathogen code, or 'all' for every pathogen (default: all).")
parser.add_argument("--force", action="store_true", help="Recompute even if a cached *_applicability_domain.csv exists.")
args = parser.parse_args()
pathogens = ALL_PATHOGENS if args.pathogen == "all" else [args.pathogen]


# Populated once per worker process by _init_worker, not once per query.
_worker_engine = None
_worker_active_ids = None


def _init_worker(fp_path, active_ids):
    global _worker_engine, _worker_active_ids
    _worker_engine = FPSim2Engine(fp_path)
    _worker_active_ids = active_ids


def _query_chunk(smiles_chunk):
    n = len(smiles_chunk)
    max_train = np.full(n, np.nan)
    max_active = np.full(n, np.nan)
    has_actives = _worker_active_ids is not None and len(_worker_active_ids) > 0
    for i, smi in enumerate(smiles_chunk):
        res = _worker_engine.similarity(smi, threshold=0.0, n_workers=1)
        if len(res) == 0:
            continue
        max_train[i] = res["coeff"].max()
        if has_actives:
            mask = np.isin(res["mol_id"], _worker_active_ids)
            if mask.any():
                max_active[i] = res["coeff"][mask].max()
    return max_train, max_active


def _rdkit_nearest_neighbour(query_fps, train_smiles, active_mask):
    """Original single-process RDKit BulkTanimotoSimilarity approach — still the
    faster option below LARGE_TRAIN_THRESHOLD (see the benchmark note above).
    Takes pre-computed query fingerprints (query_fps) rather than SMILES: the
    query set (DrugBank) is the same across every endpoint of a pathogen, so
    run_pathogen() fingerprints it once and reuses it — re-fingerprinting per
    endpoint here would redo that work once per endpoint for no reason."""
    train_fps = [AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(s), 2, nBits=2048) for s in train_smiles]
    has_actives = active_mask is not None and active_mask.any()
    max_train = np.empty(len(query_fps))
    max_active = np.full(len(query_fps), np.nan)
    for j, fp in enumerate(query_fps):
        sims = np.array(DataStructs.BulkTanimotoSimilarity(fp, train_fps))
        max_train[j] = sims.max()
        if has_actives:
            max_active[j] = sims[active_mask].max()
    return max_train, max_active


def _fpsim2_nearest_neighbour(query_smiles, train_smiles, active_mask):
    """FPSim2-backed search, optionally parallelized across N_WORKERS processes —
    the faster option at/above LARGE_TRAIN_THRESHOLD (see the benchmark note above).
    One FPSim2 index is built from train_smiles, then query_smiles is searched
    against it — in parallel once the workload is large enough for that to pay for
    its own fork overhead (see PARALLEL_MIN_COMPARISONS), single-process otherwise.
    """
    active_ids = (
        np.where(active_mask)[0] if active_mask is not None and active_mask.any() else np.array([], dtype=int)
    )
    n_query, n_train = len(query_smiles), len(train_smiles)

    with tempfile.TemporaryDirectory() as d:
        fp_path = os.path.join(d, "index.h5")
        create_db_file(
            mols_source=[(s, i) for i, s in enumerate(train_smiles)],
            filename=fp_path,
            mol_format="smiles",
            fp_type="Morgan",
            fp_params={"radius": 2, "fpSize": 2048},
        )

        if n_query * n_train < PARALLEL_MIN_COMPARISONS or N_WORKERS <= 1:
            _init_worker(fp_path, active_ids)
            max_train, max_active = _query_chunk(query_smiles)
        else:
            chunks = [c.tolist() for c in np.array_split(np.array(query_smiles, dtype=object), N_WORKERS)]
            with ProcessPoolExecutor(
                max_workers=N_WORKERS, initializer=_init_worker, initargs=(fp_path, active_ids)
            ) as ex:
                chunk_results = list(ex.map(_query_chunk, chunks))
            max_train = np.concatenate([r[0] for r in chunk_results])
            max_active = np.concatenate([r[1] for r in chunk_results])

    return max_train, max_active


def nearest_neighbour_similarity(query_smiles, query_fps, train_smiles, active_mask):
    """For every compound in query_smiles: its max Tanimoto (Morgan r=2, 2048-bit)
    similarity to train_smiles overall (max_sim_train), and to just the
    active_mask subset of it (max_sim_active). Dispatches to whichever backend
    benchmarked faster for this training-pool size (see LARGE_TRAIN_THRESHOLD).
    query_fps are query_smiles pre-fingerprinted once per pathogen (RDKit path
    only — FPSim2 fingerprints queries internally, so query_smiles is enough there).
    """
    if len(train_smiles) >= LARGE_TRAIN_THRESHOLD:
        return _fpsim2_nearest_neighbour(query_smiles, train_smiles, active_mask)
    return _rdkit_nearest_neighbour(query_fps, train_smiles, active_mask)


def compute_consensus_row(pathogen, endpoints, k_star_all, db_smiles):
    """The consensus_score row for panel A — needs only db_smiles (no fingerprints,
    no Tanimoto), so it's cheap enough to backfill into an existing cache without
    redoing that pathogen's expensive applicability-domain pass. Returns None for
    the 3 pathogens with a single sub-model (no consensus exists) or if the
    companion repo's k_star/consensus files aren't available for this pathogen.
    """
    cons_path = os.path.join(CONSENSUS_DIR, f"{pathogen}_transformed.csv")
    k_star_entry = k_star_all.get(pathogen)
    if not (k_star_entry and k_star_entry["M"] > 1 and os.path.exists(cons_path)):
        print(f"  No consensus for {pathogen} (single sub-model, or missing "
              f"12b_k_star.json / 14_consensus/{pathogen}_transformed.csv).")
        return None

    thr = round(consensus_threshold(endpoints, k_star_entry["k_star"]), 4)
    cons_scores = (
        pd.DataFrame({"smiles": db_smiles})
        .merge(pd.read_csv(cons_path)[["smiles", "consensus_score"]], on="smiles", how="left")
        ["consensus_score"].to_numpy()
    )
    n_missing = np.isnan(cons_scores).sum()
    if n_missing:
        print(f"  NOTE: {n_missing}/{len(cons_scores)} DrugBank compounds missing from "
              f"{cons_path} — excluded from the consensus hit-rate.")
    pct_above_cons = 100 * np.nanmean(cons_scores > thr)
    print(f"  consensus_score (M={k_star_entry['M']}, k_star={k_star_entry['k_star']}): "
          f"threshold={thr} above={pct_above_cons:.1f}%")
    return {
        "pathogen": pathogen, "endpoint": "consensus_score", "threshold": thr,
        "n_train": np.nan, "n_train_active": np.nan, "pct_train_active": np.nan,
        "oof_active_mean": np.nan, "oof_inactive_mean": np.nan,
        "drugbank_raw_mean": round(float(np.nanmean(cons_scores)), 4),
        "pct_drugbank_above_threshold": round(pct_above_cons, 2),
        "pct_drugbank_in_domain": np.nan,
        "corr_rank_vs_sim_to_actives": np.nan,
        "similarity_trend_top_minus_bottom_decile": np.nan,
    }


def run_pathogen(pathogen, reports_10, k_star_all, force=False):
    cache_path = os.path.join(output_dir, f"{pathogen}_applicability_domain.csv")
    endpoints = reports_10[reports_10["pathogen"] == pathogen]

    if os.path.exists(cache_path) and not force:
        cached = pd.read_csv(cache_path)
        needs_consensus = (
            not endpoints.empty
            and not (cached["endpoint"] == "consensus_score").any()
            and k_star_all.get(pathogen, {}).get("M", 1) > 1
        )
        if not needs_consensus:
            print(f"\n{pathogen}: cached at {cache_path} — loading without recomputing (--force to redo).")
            return cached
        print(f"\n{pathogen}: cached sub-models found but consensus row missing — "
              f"backfilling just that (cheap, no Tanimoto recompute).")
        score_path = require(
            os.path.join(SCORE_DIR, f"{pathogen}.csv"),
            f"Run: cd ../chembl-antimicrobial-models && eosvc download --path output/12_drugbank/score/{pathogen}.csv",
        )
        db_smiles = pd.read_csv(score_path)["smiles"].to_numpy()
        cons_row = compute_consensus_row(pathogen, endpoints, k_star_all, db_smiles)
        if cons_row is None:
            return cached
        backfilled = pd.concat([cached, pd.DataFrame([cons_row])], ignore_index=True)
        backfilled.to_csv(cache_path, index=False)
        print(f"  Saved: {cache_path}  ({len(backfilled)} endpoints, consensus row backfilled)")
        return backfilled

    print(f"\n{'='*70}\n{pathogen}\n{'='*70}")
    if endpoints.empty:
        print(f"  No retained endpoints for {pathogen} in 10_reports.csv — skipping.")
        return None

    score_path = require(
        os.path.join(SCORE_DIR, f"{pathogen}.csv"),
        f"Run: cd ../chembl-antimicrobial-models && eosvc download --path output/12_drugbank/score/{pathogen}.csv",
    )
    rank_path = require(
        os.path.join(RANK_DIR, f"{pathogen}.csv"),
        f"Run: cd ../chembl-antimicrobial-models && eosvc download --path output/12_drugbank/rank/{pathogen}.csv",
    )
    db_score_all = pd.read_csv(score_path)
    db_rank_all = pd.read_csv(rank_path)

    t_start = time.time()
    db_smiles_list = db_score_all["smiles"].tolist()
    db_valid_master = valid_smiles_idx(db_smiles_list)
    db_smiles_master = [db_smiles_list[i] for i in db_valid_master]
    # Fingerprinted once per pathogen and reused across every endpoint's RDKit-path
    # call below — the query set (DrugBank) doesn't change per endpoint, only the
    # (usually much smaller) training set does.
    db_fps_master = [AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(s), 2, nBits=2048) for s in db_smiles_master]
    print(f"  DrugBank compounds: {len(db_smiles_master)}/{len(db_score_all)} valid "
          f"({time.time()-t_start:.1f}s)")

    rows = []
    for i, row in enumerate(endpoints.itertuples()):
        internal_name = row.model_name
        t0 = time.time()

        dataset_path = os.path.join(CAM_ROOT, "output", "07_datasets", pathogen, f"{internal_name}.csv")
        local_folds = os.path.join(LOCAL_FOLDS_DIR, pathogen, f"{internal_name}_folds.json")
        cam_folds = os.path.join(CAM_ROOT, "output", "09_reports", pathogen, f"{internal_name}_folds.json")
        folds_path = local_folds if os.path.exists(local_folds) else cam_folds

        if not (os.path.exists(dataset_path) and os.path.exists(folds_path)):
            print(f"  [{i+1}/{len(endpoints)}] SKIP {internal_name}: missing training set or OOF folds")
            continue
        if internal_name not in db_score_all.columns:
            print(f"  [{i+1}/{len(endpoints)}] SKIP {internal_name}: not a column in DrugBank score file")
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

        db_raw = db_score_all[internal_name].to_numpy()[db_valid_master]
        db_rank = db_rank_all[internal_name].to_numpy()[db_valid_master]

        train_smiles_list = train["smiles"].tolist()
        train_valid = valid_smiles_idx(train_smiles_list)
        train_smiles = [train_smiles_list[i] for i in train_valid]
        train_labels = train["bin"].to_numpy()[train_valid]
        active_mask = train_labels == 1
        has_actives = active_mask.any()

        max_sim_train, max_sim_active = nearest_neighbour_similarity(
            db_smiles_master, db_fps_master, train_smiles, active_mask
        )

        threshold = round(float(row.decision_cutoff_rank), 4)
        pct_above = 100 * (db_rank > threshold).mean()
        pct_in_domain = 100 * (max_sim_train >= SIM_CUTOFF).mean()
        corr = np.corrcoef(db_rank, max_sim_active)[0, 1] if has_actives else np.nan

        d = pd.DataFrame({"rank": db_rank, "max_sim_active": max_sim_active})
        decile_trend = np.nan
        if has_actives:
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
            "drugbank_raw_mean": round(float(np.nanmean(db_raw)), 4),
            "pct_drugbank_above_threshold": round(pct_above, 2),
            "pct_drugbank_in_domain": round(pct_in_domain, 2),
            "corr_rank_vs_sim_to_actives": round(float(corr), 4) if not np.isnan(corr) else np.nan,
            "similarity_trend_top_minus_bottom_decile": (
                round(float(decile_trend), 4) if not np.isnan(decile_trend) else np.nan
            ),
        })
        print(f"  [{i+1}/{len(endpoints)}] {internal_name} (n_train={len(train)}): "
              f"above={pct_above:.1f}% in_domain={pct_in_domain:.1f}% corr={corr:.3f} "
              f"({time.time()-t0:.1f}s)")

    # No single training pool backs a consensus, so this row gets no
    # applicability-domain columns (n_train, corr, etc. stay NaN) — it's here
    # purely so panel A can mark it with a star.
    db_smiles_valid = db_score_all["smiles"].to_numpy()[db_valid_master]
    cons_row = compute_consensus_row(pathogen, endpoints, k_star_all, db_smiles_valid)
    if cons_row is not None:
        rows.append(cons_row)

    pathogen_df = pd.DataFrame(rows)
    pathogen_df.to_csv(cache_path, index=False)
    print(f"  Saved: {cache_path}  ({len(pathogen_df)} endpoints, {time.time()-t_start:.1f}s total)")
    return pathogen_df


reports_10 = pd.read_csv(REPORTS_10_PATH)
if os.path.exists(K_STAR_PATH):
    with open(K_STAR_PATH) as f:
        k_star_all = json.load(f)
else:
    print(f"NOTE: {K_STAR_PATH} not found — no consensus rows will be computed for any pathogen.\n"
          f"  Run: cd ../chembl-antimicrobial-models && eosvc download --path output/12_drugbank/12b_k_star.json")
    k_star_all = {}

t_all = time.time()
results = [run_pathogen(p, reports_10, k_star_all, force=args.force) for p in pathogens]
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

# ---------------------------------------------------------------------------
# Figure: 1 row x 3 columns (panel A sized to match B+C combined).
#   A: jittered strip by pathogen, one point per endpoint (sub-models + each
#      pathogen's consensus, where available), y = % DrugBank above threshold.
#   B: % DrugBank in applicability domain (y, shared with C) vs
#      corr(rank, similarity to training actives) (x).
#   C: same shared y vs % DrugBank above threshold (x).
# Colors are pathogen-dependent and consistent across all three panels.
# ---------------------------------------------------------------------------
if combined["pathogen"].nunique() < 2:
    print("\nFewer than 2 pathogens in the combined summary — skipping the 3-panel figure "
          "(needs at least 2 to be a meaningful comparison). Re-run with more pathogens or --pathogen all.")
    sys.exit(0)

panel1_df = combined.rename(columns={"pct_drugbank_above_threshold": "pct_above"})
bc_df = combined[combined["endpoint"] != "consensus_score"]  # B/C are applicability-domain-only; NaN for consensus

# Format: print | Style: article — change with stylia.set_format() / stylia.set_style()
stylia.set_format("print")
stylia.set_style("article")

rng = np.random.default_rng(RANDOM_SEED)
pal = stylia.CategoricalPalette("npg")
palette_colors = pal.get(len(ALL_PATHOGENS))
pathogen_color = dict(zip(ALL_PATHOGENS, palette_colors))


y_margin = 0.05 * (bc_df["pct_drugbank_in_domain"].max() - bc_df["pct_drugbank_in_domain"].min())
shared_ylim = (
    bc_df["pct_drugbank_in_domain"].min() - y_margin,
    bc_df["pct_drugbank_in_domain"].max() + y_margin,
)

fig, axs = stylia.create_figure(1, 3, width_ratios=[2, 1, 1], height=0.35)
plot_strip(axs.next(), panel1_df, pathogen_color, ALL_PATHOGENS, rng, ylabel="% DrugBank above threshold")
ax_b, ax_c = axs.next(), axs.next()
plot_indomain_vs_corr(ax_b, bc_df, pathogen_color, shared_ylim, "pct_drugbank_in_domain",
                       ylabel="% DrugBank in applicability domain")
plot_indomain_vs_hitrate(ax_c, bc_df, pathogen_color, shared_ylim, "pct_drugbank_above_threshold",
                          "pct_drugbank_in_domain", xlabel="% DrugBank above threshold",
                          ylabel="% DrugBank in applicability domain")
ax_c.sharey(ax_b)

fig_path_png = os.path.join(output_dir, "pathogen_endpoint_hit_rate.png")
fig_path_pdf = os.path.join(output_dir, "pathogen_endpoint_hit_rate.pdf")
stylia.save_figure(fig_path_png)
stylia.save_figure(fig_path_pdf)
print(f"Saved: {fig_path_png} / .pdf")
