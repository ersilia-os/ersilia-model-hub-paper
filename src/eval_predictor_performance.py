"""Steps 09/14 analysis engine — activity endpoints as predictors of each other, and property/
resemblance columns as predictors of pathogen activity.

Two independent orchestrators share the machinery below, split along a real data boundary: step 09's
:func:`run_activity_self_performance` reads ONLY the step-07 bioactivity matrix (no property data
enters it at all), while step 14's :func:`run_predictor_performance` additionally reads the step-08
physchem/abx/cytotox property blocks and treats every one of their columns as a PREDICTOR, with every
curated activity endpoint as a binary TARGET, giving one performance value per (predictor, target)
pair.

This is a descriptive association measure, not a trained model: nothing is fitted, nothing is split,
and no seed is involved. Each number is a rank statistic over the full 1.35M-compound library.

Targets
-------
The ``selected == "Yes"`` rows of ``config/08_endpoint_selection.csv``, read from step 07's cached
matrix. All are ``direction == "higher"`` (asserted), so binarization is unambiguous: the
:data:`default.ACTIVITY_BINARIZE_TOP_N` highest-scoring compounds are the positive class and every
remaining compound is negative — a rank cutoff on a fixed count, never a score threshold.

Metrics
-------
Chosen per predictor from its own value type, resolved on the FULL column (see
:func:`classify_predictor` for why a subsample is not good enough):

*   **continuous -> AUROC**, via the Mann-Whitney rank-sum identity rather than
    ``sklearn.roc_auc_score``. The predictor is ranked ONCE (ties averaged, so the result is exact);
    each target is then a 1000-element gather and sum. The direct route would be ~26k AUROC calls
    over 1.35M rows each — hours of work for numerically identical answers. Verified against sklearn
    in the step-14 script's spot-check.
*   **binary -> balanced accuracy.** AUROC of a two-valued score collapses to a single operating
    point, so the standard imbalanced-data pairing is used instead.

Both share a 0.5 chance baseline, which is the only reason the two can share a y-axis.

AUROC is reported RAW and is free to fall below 0.5 — a predictor that anti-predicts activity stays
visible as such rather than being folded to its mirror image (user-directed).

Missing values
--------------
Nothing is imputed and nothing is dropped library-wide. Two independent cases, both user-directed
and both recorded in the outputs rather than absorbed silently:

*   A TARGET with unscored compounds (currently ``eos4zfy:maip_score``, 15 of 1,355,109): an
    unscored compound cannot be claimed to be in the top N, so it is ineligible for the positive
    class and stays negative.
*   A PREDICTOR with unscored compounds (currently the 11 ``eos6ojg`` columns, the same 6 compounds
    in each): evaluated pairwise-complete over its scored subset, with ``n_compounds`` and
    ``n_unscored`` carried in the summary CSV. Filling 0 was rejected — these are similarity counts,
    where 0 means "no similar antibiotic found", which missing data does not support.
"""

import os

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

from default import (ACTIVITY_BINARIZE_TOP_N, COADD_MODEL_ID, CONSENSUS_COLUMN,
                     CURATED_PREDICTORS, PREDICTOR_CHANCE_LEVEL,
                     PREDICTOR_FAMILIES, PREDICTOR_METRICS)

#: Parquet columns read at once when reducing targets to their top-N indices. 20 columns of 1.35M
#: float64 is ~215 MB in flight; each is reduced to 1000 int32 indices and dropped immediately.
TARGET_BATCH = 20

KEY_COL = "key"


# --------------------------------------------------------------------------- #
# Targets                                                                      #
# --------------------------------------------------------------------------- #
def selected_targets(selection_csv):
    """The ``selected == "Yes"`` activity endpoints, with their step-07 matrix column names.

    Raises if any selected endpoint is not ``direction == "higher"``: the whole binarization assumes
    "top N = most active", and a ``lower`` endpoint would silently invert. Better to fail here than
    to publish an inverted target.
    """
    sel = pd.read_csv(selection_csv)
    sel = sel[sel["selected"] == "Yes"].copy()

    wrong = sel[sel["direction"] != "higher"]
    if len(wrong):
        raise ValueError(
            f"{len(wrong)} selected endpoint(s) in {os.path.basename(selection_csv)} are not "
            f"direction=higher: {wrong['column_name'].tolist()}. Binarization assumes top-N is "
            "most active — resolve the direction before running step 09 or 14.")

    sel["target"] = sel["model_id"] + ":" + sel["column_name"]
    return sel[["model_id", "column_name", "organism", "target"]].reset_index(drop=True)


def available_targets(parquet_path, targets):
    """Drop selected endpoints whose column is absent from the step-07 cache, reporting them.

    The cache is built once from the raw prediction files and holds every endpoint the selection CSV
    referenced AT THAT TIME. A model staged later (or an endpoint re-enabled after the cache was
    written) will be selected in the config but missing from the cache.

    Excluding rather than raising is user-directed, and the excluded list is both printed and written
    to ``{step}_excluded_targets.csv`` by the caller (step 09 or step 14) — a narrowed analysis must
    never be silent. The fix is to remove ``07_score_matrix_full.parquet`` and re-run
    ``07_score_matrices.py``, which re-extracts from the raw files.
    """
    available = set(pq.ParquetFile(parquet_path).schema.names)
    keep = targets["target"].isin(available)
    missing = targets[~keep]
    if len(missing):
        print(f"  [!] {len(missing)} of {len(targets)} selected endpoint(s) are ABSENT from "
              f"{os.path.basename(parquet_path)} and are EXCLUDED from this step:")
        for r in missing.itertuples():
            print(f"        {r.target}  ({r.organism})")
        print("      The cache predates their prediction file being staged. To include them, "
              "remove the parquet and re-run 07_score_matrices.py.")
    return targets[keep].reset_index(drop=True), missing.reset_index(drop=True)


def target_top_indices(parquet_path, targets, top_n=ACTIVITY_BINARIZE_TOP_N, batch=TARGET_BATCH):
    """``{target: int64[top_n]}`` of ROW POSITIONS of each target's ``top_n`` highest scores.

    Reads the step-07 matrix in column batches and reduces each column immediately, so only ~1 MB of
    indices survives rather than a 1.35M x 260 float matrix (~2.8 GB). ``np.argpartition`` is used
    rather than a full sort: only membership of the top-N matters, not the order within it.

    Row POSITIONS (not keys) are returned because every source is verified to share one key order,
    which makes the whole step positional and index-free.
    """
    available = set(pq.ParquetFile(parquet_path).schema.names)
    missing = [t for t in targets["target"] if t not in available]
    if missing:
        raise ValueError(f"{len(missing)} target column(s) absent from "
                         f"{os.path.basename(parquet_path)}: {missing[:5]}... "
                         "Re-run 07_score_matrices.py.")

    names = targets["target"].tolist()
    out, n_nan = {}, {}
    for i in range(0, len(names), batch):
        chunk = names[i:i + batch]
        df = pd.read_parquet(parquet_path, columns=chunk)
        for col in chunk:
            v = df[col].to_numpy(dtype=float)
            nan = int(np.isnan(v).sum())
            if nan:
                n_nan[col] = nan
            # An unscored compound cannot be claimed to be in the top N, so -inf makes it
            # ineligible for the positive class; it stays in the negative class rather than being
            # removed from the library. User-directed (see the note below). NaN would otherwise
            # sort HIGH under argpartition and silently become a positive.
            v = np.where(np.isnan(v), -np.inf, v)
            idx = np.argpartition(-v, top_n - 1)[:top_n]
            out[col] = np.sort(idx)
        del df
        print(f"  [targets] {min(i + batch, len(names))}/{len(names)} reduced to top {top_n}")

    # Reported, never silently absorbed: currently only eos4zfy:maip_score, 15 of 1,355,109
    # compounds — the same partially-scored compounds step 09 keeps and tracks. Signed off to stay
    # in the negative class. A column missing enough scores to matter would show up here loudly.
    if n_nan:
        total = sum(n_nan.values())
        print(f"  [targets] {len(n_nan)} column(s) have unscored compounds ({total} values total): "
              f"{n_nan} — kept in the negative class, none dropped")
    return out


# --------------------------------------------------------------------------- #
# Predictors                                                                   #
# --------------------------------------------------------------------------- #
def predictor_index(property_csvs, families=PREDICTOR_FAMILIES):
    """Every predictor column across the property blocks, as a metadata DataFrame.

    ``property_csvs`` is an iterable of per-family matrix CSVs (all step 08). The family is
    read from each column's own ``{family}__`` prefix, not from which file it came from, so the
    split into one CSV per family cannot silently mislabel a predictor.

    Columns: ``predictor`` (the full ``{family}__{model_id}__{column}`` name), ``family``,
    ``model_id``, ``column_name``, ``source`` (which CSV to read it from).
    """
    rows = []
    for path in property_csvs:
        header = pd.read_csv(path, nrows=0).columns
        for c in header:
            parts = c.split("__")
            if len(parts) != 3 or parts[0] not in families:
                continue  # key / input / anything not following the convention
            rows.append({"predictor": c, "family": parts[0], "model_id": parts[1],
                         "column_name": parts[2], "source": path})
    index = pd.DataFrame(rows)
    order = {f: i for i, f in enumerate(families)}
    index["_o"] = index["family"].map(order)
    return index.sort_values(["_o", "predictor"]).drop(columns="_o").reset_index(drop=True)


def classify_predictor(values):
    """``"binary"`` if the column's distinct values are a subset of {0, 1}, else ``"continuous"``.

    Deliberately evaluated on the FULL column. ``physchem__eos4djh__n_radical_electrons`` is 0 for
    all but roughly 1 compound in 5,000, so it presents as binary in any subsample while being a
    count; only the full column settles it. The resolved type is written to the output CSV so the
    call is reviewable rather than implicit.
    """
    u = pd.unique(values[~np.isnan(values)])
    return "binary" if set(np.asarray(u, dtype=float)) <= {0.0, 1.0} else "continuous"


# --------------------------------------------------------------------------- #
# Metric kernels                                                               #
# --------------------------------------------------------------------------- #
def auroc_from_ranks(ranks, top_idx, n_total):
    """AUROC of a continuous predictor against a top-N binary target, from precomputed ranks.

    The Mann-Whitney identity: with tie-averaged ranks over all ``n_total`` compounds,
    ``AUROC = (R_pos - n_pos(n_pos+1)/2) / (n_pos * n_neg)``, where ``R_pos`` is the rank sum of the
    positive class. Exact under ties, and identical to ``sklearn.roc_auc_score``.
    """
    n_pos = len(top_idx)
    n_neg = n_total - n_pos
    if n_pos == 0 or n_neg == 0:
        return np.nan
    r_pos = float(ranks[top_idx].sum())
    return (r_pos - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)


def balanced_accuracy_from_mask(mask, n_mask_true, top_idx, n_total):
    """Balanced accuracy of a binary predictor against a top-N binary target.

    ``mask`` is the predictor's boolean "predicted positive" vector and ``n_mask_true`` its total,
    both computed once per predictor; only the ``top_idx`` gather is per-target.
    Returns NaN for a constant predictor (all-0 or all-1), where the measure is undefined rather
    than 0.5 — a degenerate column must not masquerade as a chance-level one.
    """
    n_pos = len(top_idx)
    n_neg = n_total - n_pos
    if n_pos == 0 or n_neg == 0 or n_mask_true == 0 or n_mask_true == n_total:
        return np.nan
    tp = int(mask[top_idx].sum())
    fp = n_mask_true - tp
    tn = n_neg - fp
    return 0.5 * (tp / n_pos + tn / n_neg)


# --------------------------------------------------------------------------- #
# Activity endpoints as predictors of each other                               #
# --------------------------------------------------------------------------- #
def pathogen_subset_endpoints(targets, pathogens_csv, consensus_col=CONSENSUS_COLUMN):
    """The activity endpoints of the 15 pathogens of interest, with consensus models collapsed.

    Two reductions, in order:

    1.  Keep only endpoints whose organism is one of ``config/pathogens_of_interest.csv``, matched
        EXACTLY. The two configs were aligned on 2026-08-07 (they previously spelled Campylobacter
        and Enterobacter differently and needed an alias map). Never match on a genus substring if
        they diverge again: that would wrongly capture *C. glabrata* for *C. albicans* and two
        non-pneumoniae streptococci for *S. pneumoniae* — fix the spelling instead.
    2.  Per (model_id, organism): if that model publishes ``consensus_col`` for the organism, keep
        ONLY that column; otherwise keep all of its endpoints. Grouping by (model_id, organism)
        rather than model_id matters — eos3dys spans six organisms and has no consensus column, so
        each of its organisms keeps its own endpoints.

    Returns the same shape as :func:`selected_targets` plus a ``pathogen`` column (the canonical
    name from the pathogen config) and ``is_consensus``.
    """
    pathogens = pd.read_csv(pathogens_csv)
    wanted = {p: p for p in pathogens["pathogen"]}

    sub = targets[targets["organism"].isin(wanted)].copy()
    sub["pathogen"] = sub["organism"].map(wanted)

    absent = sorted(set(pathogens["pathogen"]) - set(sub["pathogen"]))
    if absent:
        print(f"  [pathogen-subset] {len(absent)} pathogen(s) have no selected endpoint and are "
              f"absent from this figure: {absent}")

    keep = []
    for (_, _), g in sub.groupby(["model_id", "organism"], sort=False):
        has = (g["column_name"] == consensus_col).any()
        keep.append(g[g["column_name"] == consensus_col] if has else g)
    out = pd.concat(keep).copy()
    out["is_consensus"] = out["column_name"] == consensus_col
    out = out.sort_values(["pathogen", "model_id", "column_name"]).reset_index(drop=True)

    print(f"  [pathogen-subset] {len(sub)} endpoints -> {len(out)} after collapsing consensus "
          f"models, across {out['pathogen'].nunique()} pathogens "
          f"({int(out['is_consensus'].sum())} consensus, {int((~out['is_consensus']).sum())} single)")
    return out


def activity_self_performance(parquet_path, targets, tops, n_total, batch=TARGET_BATCH):
    """The 260 x 260 block: each activity endpoint's RAW score against every endpoint's top-N
    binarization, as AUROC.

    Same machinery as the property predictors, with the activity endpoints on both sides — the
    x-axis entity is the un-binarized score, the y-axis target is the binarized one. Every activity
    endpoint is continuous, so AUROC applies throughout and no metric selection is needed.

    Three pair kinds are labelled rather than filtered, since each answers a different question:

    *   ``self_pair`` — an endpoint against its OWN binarization. AUROC is 1.0 by construction (the
        top N of a score are exactly its N highest values), so it measures nothing; it is retained
        only as a correctness check that the pipeline returns 1.0 where it must, and is excluded
        from the figure.
    *   ``same_organism`` — two endpoints of the same pathogen. Expected to agree; this is model
        self-consistency, not cross-pathogen predictive power. Note that the 41 organisms with a
        single selected endpoint have NO same-organism pairs once self-pairs are removed.
    *   ``same_model`` — two columns of one model, which share training data.
    """
    from scipy.stats import rankdata

    meta = targets.set_index("target")
    names = targets["target"].tolist()
    records = []

    for i in range(0, len(names), batch):
        chunk = names[i:i + batch]
        df = pd.read_parquet(parquet_path, columns=chunk)
        for pcol in chunk:
            v = df[pcol].to_numpy(dtype=float)
            n_nan = int(np.isnan(v).sum())
            # Same pairwise-complete rule as the property predictors (see run_predictor_performance).
            n_used, remap = n_total, None
            if n_nan:
                valid = ~np.isnan(v)
                n_used = int(valid.sum())
                remap = np.full(n_total, -1, dtype=np.int64)
                remap[valid] = np.arange(n_used)
                v = v[valid]
            ranks = rankdata(v, method="average")
            pm = meta.loc[pcol]

            for tcol, idx in tops.items():
                if remap is not None:
                    idx = remap[idx]
                    idx = idx[idx >= 0]
                tm = meta.loc[tcol]
                rec = {
                    "predictor_endpoint": pcol, "predictor_model_id": pm["model_id"],
                    "predictor_column": pm["column_name"], "predictor_organism": pm["organism"],
                    "target_endpoint": tcol, "target_model_id": tm["model_id"],
                    "target_column": tm["column_name"], "target_organism": tm["organism"],
                    "metric": PREDICTOR_METRICS["continuous"],
                    "value": auroc_from_ranks(ranks, idx, n_used),
                    "self_pair": pcol == tcol,
                    "same_organism": pm["organism"] == tm["organism"],
                    "same_model": pm["model_id"] == tm["model_id"],
                }
                # Carried through when the caller passed the pathogen-subset table, which adds the
                # canonical pathogen name and the consensus flag on top of selected_targets' columns.
                for extra in ("pathogen", "is_consensus"):
                    if extra in meta.columns:
                        rec[f"predictor_{extra}"] = pm[extra]
                        rec[f"target_{extra}"] = tm[extra]
                records.append(rec)
        del df
        print(f"  [activity-self] {min(i + batch, len(names))}/{len(names)} endpoints ranked")

    return pd.DataFrame(records)


def curated_predictor_performance(property_csvs, targets, tops, n_total,
                                  curated=CURATED_PREDICTORS):
    """The hand-picked :data:`default.CURATED_PREDICTORS` scored against a given target set.

    Same kernels as :func:`run_predictor_performance`, but over a shortlist rather than all 101
    property columns, and
    against the consensus-collapsed pathogen targets rather than every selected endpoint. Every
    curated predictor is continuous, so the whole block is AUROC — verified here rather than
    assumed, since a binary column slipping in would silently be scored with the wrong metric.
    """
    from scipy.stats import rankdata

    index = predictor_index(property_csvs)
    wanted = {(f, c) for f, cols in curated.items() for c in cols}
    sel = index[[(r.family, r.column_name) in wanted for r in index.itertuples()]].copy()

    missing = wanted - {(r.family, r.column_name) for r in sel.itertuples()}
    if missing:
        raise ValueError(f"Curated predictor(s) not found in the property blocks: {sorted(missing)}")

    meta = targets.set_index("target")
    records = []
    for source, group in sel.groupby("source", sort=False):
        block = pd.read_csv(source, usecols=group["predictor"].tolist())
        for r in group.itertuples():
            v = block[r.predictor].to_numpy(dtype=float)
            n_nan = int(np.isnan(v).sum())
            n_used, remap = n_total, None
            if n_nan:  # pairwise-complete, as in run_predictor_performance
                valid = ~np.isnan(v)
                n_used = int(valid.sum())
                remap = np.full(n_total, -1, dtype=np.int64)
                remap[valid] = np.arange(n_used)
                v = v[valid]

            ptype = classify_predictor(v)
            if ptype != "continuous":
                raise ValueError(
                    f"Curated predictor {r.predictor} is {ptype}, not continuous — this figure "
                    "puts every predictor on one AUROC axis and cannot mix in a binary column.")

            ranks = rankdata(v, method="average")
            for target, idx in tops.items():
                if remap is not None:
                    idx = remap[idx]
                    idx = idx[idx >= 0]
                tm = meta.loc[target]
                records.append({
                    "predictor": r.predictor, "family": r.family,
                    "predictor_model_id": r.model_id, "predictor_column": r.column_name,
                    "metric": PREDICTOR_METRICS["continuous"],
                    "target_endpoint": target, "target_model_id": tm["model_id"],
                    "target_column": tm["column_name"], "organism": tm["organism"],
                    "pathogen": tm.get("pathogen"), "is_consensus": tm.get("is_consensus"),
                    "value": auroc_from_ranks(ranks, idx, n_used),
                    # cytotox__eos3dys__cytotoxicity_ic50 shares a model with the eos3dys targets.
                    "same_model": r.model_id == tm["model_id"],
                })
            print(f"  {r.predictor}: median AUROC "
                  f"{np.nanmedian([x['value'] for x in records if x['predictor'] == r.predictor]):.3f}")
        del block
    return pd.DataFrame(records)


# --------------------------------------------------------------------------- #
# Orchestrator                                                                 #
# --------------------------------------------------------------------------- #
def _bioactivity_total(parquet_path):
    """Row count of the step-07 matrix. No property-CSV alignment check is needed here: step 09's
    orchestrator never reads property data, so there is nothing to align it against.
    """
    n = pq.ParquetFile(parquet_path).metadata.num_rows
    print(f"[activity-self] {n} compounds in {os.path.basename(parquet_path)}")
    return n


def _assert_key_alignment(property_csvs, parquet_path):
    """Every source must share one key ORDER, since the whole step is positional."""
    # Step 07's matrix carries `key` as the parquet INDEX, not as a column.
    k07 = pd.read_parquet(parquet_path, columns=[]).index.to_series().reset_index(drop=True)
    for path in property_csvs:
        k = pd.read_csv(path, usecols=[KEY_COL])[KEY_COL]
        if not k.equals(k07):
            raise ValueError(
                f"Key order in {os.path.basename(path)} differs from the step-07 matrix; step 14 "
                "compares them positionally. Rebuild them from the same reference library before "
                "continuing.")
    print(f"[predictor-perf] key order verified across {len(property_csvs) + 1} sources "
          f"({len(k07)} compounds)")
    return len(k07)


def run_activity_self_performance(parquet_path, selection_csv, pathogens_csv, output_dir,
                                  top_n=ACTIVITY_BINARIZE_TOP_N):
    """Step 09 orchestrator — activity endpoints as predictors of each other.

    Bioactivity-only: reads nothing but the step-07 matrix, no property/resemblance data enters this
    analysis at all, which is exactly what lets it live in step 09 rather than alongside the
    property-predictor analysis in step 14. Writes ``09_activity_self_performance.csv`` (all
    selected targets) and ``09_pathogen_subset_self_performance.csv`` (restricted to the 15 pathogens
    of interest, consensus models collapsed).
    """
    n_total = _bioactivity_total(parquet_path)

    targets, excluded = available_targets(parquet_path, selected_targets(selection_csv))
    if len(excluded):
        excluded.to_csv(os.path.join(output_dir, "09_excluded_targets.csv"), index=False)
        print(f"      -> 09_excluded_targets.csv ({len(excluded)} rows)")
    print(f"[activity-self] {len(targets)} activity targets, binarized at top {top_n} "
          f"(prevalence {top_n / n_total:.5%})")
    tops = target_top_indices(parquet_path, targets, top_n=top_n)

    print(f"\n[activity-self] {len(targets)} x {len(targets)} endpoint block")
    self_perf = activity_self_performance(parquet_path, targets, tops, n_total)
    self_path = os.path.join(output_dir, "09_activity_self_performance.csv")
    self_perf.to_csv(self_path, index=False)

    # The diagonal is 1.0 by construction; anything else means the ranking and the top-N selection
    # have drifted apart, so it is checked rather than assumed.
    diag = self_perf.loc[self_perf["self_pair"], "value"]
    print(f"  self-pair AUROC: min {diag.min():.6f}, max {diag.max():.6f} "
          f"({len(diag)} pairs — must all be 1.0)")
    off = self_perf[~self_perf["self_pair"]]
    print(f"  same-organism pairs: {int(off['same_organism'].sum())}, "
          f"cross-organism: {int((~off['same_organism']).sum())}")
    print(f"  -> {os.path.basename(self_path)} ({len(self_perf)} rows)")

    print(f"\n[pathogen-subset] activity endpoints of the pathogens of interest")
    sub_targets = pathogen_subset_endpoints(targets, pathogens_csv)
    sub_tops = {t: tops[t] for t in sub_targets["target"]}
    subset = activity_self_performance(parquet_path, sub_targets, sub_tops, n_total)
    sub_path = os.path.join(output_dir, "09_pathogen_subset_self_performance.csv")
    subset.to_csv(sub_path, index=False)
    sub_diag = subset.loc[subset["self_pair"], "value"]
    print(f"  self-pair AUROC: min {sub_diag.min():.6f}, max {sub_diag.max():.6f} "
          f"({len(sub_diag)} pairs — must all be 1.0)")
    print(f"  -> {os.path.basename(sub_path)} ({len(subset)} rows)")

    return targets, tops, self_perf, subset


def run_predictor_performance(property_csvs, parquet_path, selection_csv, output_dir,
                              pathogens_csv=None, top_n=ACTIVITY_BINARIZE_TOP_N,
                              coadd_model=COADD_MODEL_ID):
    """Step 14 orchestrator. Writes the long per-pair CSV and the per-predictor summary."""
    from scipy.stats import rankdata

    n_total = _assert_key_alignment(property_csvs, parquet_path)

    targets, excluded = available_targets(parquet_path, selected_targets(selection_csv))
    if len(excluded):
        excluded.to_csv(os.path.join(output_dir, "14_excluded_targets.csv"), index=False)
        print(f"      -> 14_excluded_targets.csv ({len(excluded)} rows)")
    print(f"[predictor-perf] {len(targets)} activity targets, binarized at top {top_n} "
          f"(prevalence {top_n / n_total:.5%})")
    tops = target_top_indices(parquet_path, targets, top_n=top_n)

    predictors = predictor_index(property_csvs)
    print(f"[predictor-perf] {len(predictors)} predictors: "
          f"{predictors['family'].value_counts().to_dict()}")

    target_meta = targets.set_index("target")
    records, summaries = [], []

    for source, group in predictors.groupby("source", sort=False):
        cols = group["predictor"].tolist()
        block = pd.read_csv(source, usecols=cols)
        for r in group.itertuples():
            v = block[r.predictor].to_numpy(dtype=float)
            n_nan = int(np.isnan(v).sum())

            # PAIRWISE-COMPLETE (user-directed): a predictor with unscored compounds is evaluated
            # over its scored subset only. No value is imputed — these are similarity counts, where
            # filling 0 would assert "no similar antibiotic found", a claim missing data does not
            # support — and no compound is removed from any other predictor's evaluation.
            # `remap` renumbers full-library row positions onto the compacted subset so each
            # target's positive indices still line up; positives that fall on an unscored compound
            # simply drop out of that pair.
            n_used, remap = n_total, None
            if n_nan:
                valid = ~np.isnan(v)
                n_used = int(valid.sum())
                remap = np.full(n_total, -1, dtype=np.int64)
                remap[valid] = np.arange(n_used)
                v = v[valid]
                print(f"  [{r.predictor}] {n_nan} unscored compound(s) excluded "
                      f"— evaluated over {n_used}")

            ptype = classify_predictor(v)
            metric = PREDICTOR_METRICS[ptype]
            if ptype == "continuous":
                ranks = rankdata(v, method="average")
                score_fn = lambda idx, _r=ranks, _n=n_used: auroc_from_ranks(_r, idx, _n)
            else:
                mask = v == 1
                n_mask = int(mask.sum())
                score_fn = lambda idx, _m=mask, _k=n_mask, _n=n_used: \
                    balanced_accuracy_from_mask(_m, _k, idx, _n)

            vals = []
            for target, idx in tops.items():
                if remap is not None:
                    idx = remap[idx]
                    idx = idx[idx >= 0]
                value = score_fn(idx)
                vals.append(value)
                meta = target_meta.loc[target]
                records.append({
                    "predictor": r.predictor, "family": r.family,
                    "predictor_model_id": r.model_id, "predictor_column": r.column_name,
                    "predictor_type": ptype, "metric": metric,
                    "target_endpoint": target, "target_model_id": meta["model_id"],
                    "target_column": meta["column_name"], "organism": meta["organism"],
                    "value": value,
                    # Predictor and target from the SAME model share training data, so these pairs
                    # are not independent evidence. Flagged, never silently dropped.
                    "same_model": r.model_id == meta["model_id"],
                })
            a = np.asarray(vals, dtype=float)
            finite = a[np.isfinite(a)]
            summaries.append({
                "predictor": r.predictor, "family": r.family, "model_id": r.model_id,
                "column_name": r.column_name, "predictor_type": ptype, "metric": metric,
                "n_compounds": n_used, "n_unscored": n_nan,
                "n_targets": len(a), "n_defined": len(finite),
                "median": np.median(finite) if len(finite) else np.nan,
                "q1": np.percentile(finite, 25) if len(finite) else np.nan,
                "q3": np.percentile(finite, 75) if len(finite) else np.nan,
                "min": finite.min() if len(finite) else np.nan,
                "max": finite.max() if len(finite) else np.nan,
            })
            flag = "" if len(finite) == len(a) else f"  [{len(a) - len(finite)} undefined]"
            print(f"  {r.predictor} ({ptype}, {metric}): "
                  f"median {summaries[-1]['median']:.3f}{flag}")
        del block

    perf = pd.DataFrame(records)
    summary = pd.DataFrame(summaries)

    perf_path = os.path.join(output_dir, "14_predictor_performance.csv")
    summary_path = os.path.join(output_dir, "14_predictor_summary.csv")
    perf.to_csv(perf_path, index=False)
    summary.to_csv(summary_path, index=False)

    n_same = int(perf["same_model"].sum())
    n_undef = int(perf["value"].isna().sum())
    print(f"\n  -> {os.path.basename(perf_path)} ({len(perf)} rows; "
          f"{n_same} same-model pairs flagged, {n_undef} undefined)")
    print(f"  -> {os.path.basename(summary_path)} ({len(summary)} rows)")
    print(f"  chance level for both metrics: {PREDICTOR_CHANCE_LEVEL}")

    # --- The curated 12-predictor shortlist, restricted to the 15 pathogens of interest with
    # consensus models collapsed. (Their bioactivity-only self-performance is step 09's job; this
    # is the property-predictor shortlist scored against the same target set.) ---
    curated = None
    if pathogens_csv:
        sub_targets = pathogen_subset_endpoints(targets, pathogens_csv)
        sub_tops = {t: tops[t] for t in sub_targets["target"]}
        print(f"\n[curated] {sum(len(v) for v in CURATED_PREDICTORS.values())} predictors "
              f"x {len(sub_targets)} pathogen endpoints")
        curated = curated_predictor_performance(property_csvs, sub_targets, sub_tops,
                                                n_total)
        cur_path = os.path.join(output_dir, "14_curated_predictor_performance.csv")
        curated.to_csv(cur_path, index=False)
        print(f"  -> {os.path.basename(cur_path)} ({len(curated)} rows; "
              f"{int(curated['same_model'].sum())} same-model pairs flagged)")

    return perf, summary, curated
