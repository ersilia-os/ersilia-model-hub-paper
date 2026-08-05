"""Step 05 analysis engine — evaluate the ChEMBL pathogen models on EU OpenScreen.

Consumes the step-04 Ersilia prediction CSVs
(``output/04_ersilia_predictions/euopenscreen/{eosid}.csv``, headline column ``consensus_score``,
SMILES column ``input``) plus the EU OpenScreen ground-truth label files, and writes small
pre-aggregated summary CSVs into ``output/05_euopenscreen_validation/`` that the figures
(``plots_euopenscreen``) consume. The figures never touch per-molecule data.

Analyses:
  1. Own primary assay — each of the 7 organisms with an EU OpenScreen primary assay, scored by
     its own model (AUROC/AUPRC/BEDROC/EF + ROC curves).
  3. Shared vs exclusive hits — precomputed EU OpenScreen exclusivity subsets.
  4. Cross-organism — every model x every EU OpenScreen assay (off-diagonal = a model predicting
     a DIFFERENT organism's data) + per-model specificity index.
  4b. Active-set overlap + hit promiscuity — label-only views of the 7 assays' hit sets: pairwise
     Jaccard, and how many actives are hits in 1, 2, ... 7 pathogens (which compounds are the
     promiscuous, pan-active ones).
  + Per-submodel breakdown (all output columns, not just consensus) and a training-set leakage
    report.

Every metric is reported both ``raw`` and ``dedup`` (InChIKey-deduplicated against the model's
ChEMBL training set); dedup degrades gracefully to raw-only if the training repo is absent.
Missing step-04 model files are skipped with a logged message (no silent caps). The shared
IO/merge/metric primitives live in :mod:`eval_common`; CoAdd validation lives in
:mod:`eval_coadd`.
"""

import os

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, roc_curve

from default import (
    ACTIVITY_CLASSES,
    DEDUP_SUBDIR,
    FULL_SUBDIR,
    HIT_CLASSES,
    NARROW_MAX_PATHOGENS,
    SHARED_ORGANISMS,
)
from metrics import compute_metrics
from eval_common import (
    METRIC_KEYS,
    evaluate,
    load_predictions,
    merged_variants,
    training_inchikeys,
)


# --------------------------------------------------------------------------- #
# Loaders                                                                      #
# --------------------------------------------------------------------------- #
def load_euos_primary(euos_root, code, assay_id):
    """EU OpenScreen primary-assay labels as DataFrame[smiles, bin, inchikey].

    ONE ROW PER COMPOUND. Reads ``02_binarised_assays/{assay_id}.csv`` (smiles, bin), keeps only
    conclusive rows (bin in {0, 1}), collapses duplicate SMILES keeping the HIGHEST bin, and enriches
    inchikey from ``02_merged/02_{code}.csv``.

    **Why the primary assay and not ``02_merged``.** ``02_merged/02_{code}.csv`` is the union of ALL
    5-6 live EU OpenScreen assays for the pathogen, not its primary screen. The non-primary ones are small
    (~5,300 compounds against ~101,000) but hit-ENRICHED 12x to 60x, because they include academic
    compound batches and dose-response confirmation sets that run >50% active. Pooling them multiplies
    the active count by 1.4x-3.5x and destroys the ~1e-4 prevalence that AUROC, BEDROC and EF@k are
    read against. The primary assay is the only uniformly screened, unbiased set, so it is the source
    of truth for "did this compound inhibit this organism". ``02_merged`` is used here ONLY for the
    InChIKey lookup (it is a superset of every primary assay's SMILES, so the join is complete), and
    elsewhere only as ``02_only_smiles.csv``, the step-04 prediction input.

    **Why the de-duplication.** The per-assay file is MEASUREMENT-level: one row per assay well. The
    EU OpenScreen library registers some structures more than once as separate compounds, so those
    appear as repeated SMILES. Left alone they are counted twice in ``n_eval``/``n_active`` and carry
    double weight in every metric, and they made this loader disagree with the compound-level
    exclusivity subsets in ``06_subset_data/exclusivity/`` by up to one active per organism.

    The rule is **max bin, "Active prevails"** — deliberately the same rule upstream applies in
    ``eu-openscreen-antimicrobial-tasks/scripts/02_binarise_and_merge.py::merge_pathogen_rows`` when
    it builds ``02_merged``. Upstream simply does not apply it to the per-assay files. Keeping the
    rules identical is the point: the two repos must not resolve the same conflict differently.

    It matters for real data, not just for counting. *A. baumannii* has a DISCORDANT pair — two
    separately registered library compounds, EOS101879 (93.05, active) and EOS17004 (11.03,
    inactive), with the same structure and InChIKey. That is a replicate disagreement, and this rule
    resolves it to active. Before this, it was resolved by whichever row happened to come first.
    """
    path = os.path.join(euos_root, "02_binarised_assays", f"{assay_id}.csv")
    if not os.path.exists(path):
        print(f"  [skip] binarised assay {assay_id}.csv not found")
        return None
    labels = pd.read_csv(path)
    labels = labels[labels["bin"].isin([0, 1])].copy()
    # Collapse measurement rows to compounds. Logged rather than silent: if upstream ever starts
    # shipping compound-level per-assay files this count goes to zero, and if the library grows new
    # duplicate registrations it goes up — either way the change should be visible in the run log.
    n_before = len(labels)
    # Counted BEFORE collapsing, while both rows of a conflicting pair are still present.
    n_disc = int(labels.groupby("smiles")["bin"].nunique().gt(1).sum())
    labels = (labels.sort_values("bin", ascending=False)
                    .drop_duplicates("smiles", keep="first"))
    n_collapsed = n_before - len(labels)
    if n_collapsed:
        disc_note = f", {n_disc} discordant (resolved active)" if n_disc else ""
        print(f"  [dedup] {assay_id}: collapsed {n_collapsed} duplicate label row(s) to "
              f"{len(labels)} compounds{disc_note}")
    merged_path = os.path.join(euos_root, "02_merged", f"02_{code}.csv")
    if os.path.exists(merged_path):
        inchi = pd.read_csv(merged_path, usecols=["smiles", "inchikey"]).drop_duplicates("smiles")
        labels = labels.merge(inchi, on="smiles", how="left")
    else:
        labels["inchikey"] = np.nan
    return labels[["smiles", "bin", "inchikey"]]


# --------------------------------------------------------------------------- #
# ROC-curve helpers (EU OpenScreen own-assay only)                             #
# --------------------------------------------------------------------------- #
def _thin_curve(fpr, tpr, cap=800):
    """Down-sample a ROC curve to <= ``cap`` points, keeping shape-defining vertices.

    Always keeps the endpoints and every TPR-jump vertex (where a positive is recovered) so
    the few-positives staircase is preserved exactly; the remaining budget is filled with a
    uniform subsample of the (dense, negative-driven) FPR steps."""
    n = len(fpr)
    if n <= cap:
        return fpr, tpr
    keep = np.zeros(n, dtype=bool)
    keep[0] = keep[-1] = True
    keep[1:] |= (tpr[1:] != tpr[:-1])          # TPR jumps = positives recovered
    budget = cap - int(keep.sum())
    if budget > 0:
        rest = np.where(~keep)[0]
        keep[rest[np.linspace(0, len(rest) - 1, budget).round().astype(int)]] = True
    return fpr[keep], tpr[keep]


def _roc_records(pred, labels, train_keys, pathogen, code):
    """ROC-curve points per set for one own-assay evaluation.

    Uses sklearn ``roc_curve`` (drop_intermediate=True → compact staircase), so the
    few-positives shape is preserved. Returns a list of {pathogen, code, set, fpr, tpr,
    n_pos, n_neg, auroc} rows (one per curve vertex)."""
    out = []
    for set_name, sub in merged_variants(pred, labels, train_keys):
        if sub["bin"].nunique() < 2:
            continue
        y, s = sub["bin"].values, sub["score"].values
        fpr, tpr, _ = roc_curve(y, s)
        n_pos, n_neg = int(y.sum()), int((y == 0).sum())
        auroc = round(roc_auc_score(y, s), 4)
        fpr, tpr = _thin_curve(fpr, tpr)
        for x, t in zip(fpr, tpr):
            out.append({"pathogen": pathogen, "code": code, "set": set_name,
                        "fpr": round(float(x), 5), "tpr": round(float(t), 5),
                        "n_pos": n_pos, "n_neg": n_neg, "auroc": auroc})
    return out


# --------------------------------------------------------------------------- #
# Analyses                                                                     #
# --------------------------------------------------------------------------- #
def run_euopenscreen(pred_dir, euos_root, config, train_cache):
    """Analyses 1 (own-assay), 3 (exclusivity) and 4 (cross-organism) for EU OpenScreen.

    Returns (own_records, exclusivity_records, cross_records, roc_records,
    exclusivity_pct_records). ``config`` is the
    pathogen table (pathogen, code, eosid); ``train_cache`` is a {code: set|None} of keys.
    """
    primary = pd.read_csv(os.path.join(euos_root, "primary_assays_manual.csv"))
    code_to_assay = dict(zip(primary["pathogen_code"], primary["assay_eos_id"]))
    name_by_code = dict(zip(config["code"], config["pathogen"]))

    # Pre-load the 7 shared primary assays once (reused by own + cross analyses).
    assays = {}
    for code in SHARED_ORGANISMS:
        assay_id = code_to_assay.get(code)
        if assay_id is None:
            continue
        lab = load_euos_primary(euos_root, code, assay_id)
        if lab is not None:
            assays[code] = lab

    own_records, cross_records, exclusivity_records, roc_records = [], [], [], []
    exclusivity_pct_records = []

    for _, row in config.iterrows():
        code, eosid, pathogen = row["code"], row["eosid"], row["pathogen"]
        pred = load_predictions(pred_dir, "euopenscreen", eosid)
        if pred is None:
            continue
        train_keys = train_cache.get(code)

        # Analysis 4: this model vs every shared assay (diagonal = own organism).
        for assay_code, lab in assays.items():
            base = {
                "model_pathogen": pathogen, "model_code": code, "model_eosid": eosid,
                "assay_pathogen": name_by_code.get(assay_code, assay_code),
                "assay_code": assay_code,
            }
            cross_records.extend(evaluate(pred, lab, train_keys, base))

        # Analysis 1: own primary assay (only for shared organisms).
        if code in assays:
            base = {"pathogen": pathogen, "code": code, "eosid": eosid}
            own_records.extend(evaluate(pred, assays[code], train_keys, base))
            roc_records.extend(_roc_records(pred, assays[code], train_keys, pathogen, code))

            # Analysis 3: exclusive vs non-exclusive (shared) hits vs primary inactives.
            inactives = assays[code][assays[code]["bin"] == 0].copy()
            for mode in ("exclusive", "nonexclusive"):
                sub = _load_exclusivity_task(euos_root, code, mode, inactives)
                if sub is None:
                    continue
                base = {"pathogen": pathogen, "code": code, "eosid": eosid, "subset": mode}
                exclusivity_records.extend(evaluate(pred, sub, train_keys, base))
                exclusivity_pct_records.extend(
                    _exclusivity_percentiles(pred, sub, train_keys, base))

    return (own_records, exclusivity_records, cross_records, roc_records,
            exclusivity_pct_records)


def _exclusivity_percentiles(pred, labels, train_keys, base):
    """Where each active lands in the model's own ranking, as a percentile in [0, 1].

    The exclusivity AUROC says how well a model separates one subset of hits from the primary
    inactives, but a single number cannot show *how*: whether the hits pile up at the top, or
    whether a good score is a handful of hits at the very top dragging a long tail. This exports
    the distribution behind the number, one row per active, for the event-plot panel.

    Percentile is the active's rank among ALL evaluated compounds (its subset's actives plus the
    same primary inactives the AUROC uses), so 1.0 = top-scoring compound in the library and 0.5 =
    the middle. Ties take the average rank, matching how AUROC treats them.

    The two are the same quantity: with 1-based ranks,
    ``AUROC = (mean_rank_pos - (n_pos + 1) / 2) / n_neg``, so the mean of these percentiles is
    ``(auroc * n_neg + (n_pos + 1) / 2) / n``. At this prevalence (~1e-4) that is the AUROC to
    four decimal places. :func:`_check_percentiles_match_auroc` asserts the exact identity — the
    ridge's centre of mass IS the bar it replaces, and the export is wrong if it is not.
    """
    rows = []
    for set_name, sub in merged_variants(pred, labels, train_keys):
        if sub["bin"].nunique() < 2:
            continue
        pct = sub["score"].rank(pct=True, method="average")
        active = sub["bin"].values == 1
        for key, smi, p in zip(sub.loc[active, "inchikey"], sub.loc[active, "smiles"],
                               pct.values[active]):
            rows.append({**base, "set": set_name, "inchikey": key, "smiles": smi,
                         "percentile": float(p)})
    return rows


def _check_percentiles_match_auroc(pct_df, excl_df, tol=1e-4):
    """Assert the exported percentiles reproduce the AUROC they are meant to unpack.

    Guards the export against a silent mis-join: a percentile taken over the wrong compound set,
    or actives matched to the wrong subset, would still look like a plausible distribution but
    would no longer integrate to the published AUROC. Reports the worst deviation.

    ``tol`` cannot go below 1e-4. ``metrics.compute_metrics`` stores ``round(auroc, 4)``, so the
    comparison is against a value already quantised to 1e-4 and a deviation of up to 5e-5 is the
    rounding, not a defect — measured worst case is 4.86e-05. Tightening this further would fail
    on correct data; catching a real mis-join does not need more resolution than this, since a
    wrong compound set moves the mean percentile by orders of magnitude more.
    """
    if pct_df.empty or excl_df.empty:
        return None
    keys = ["pathogen", "code", "subset", "set"]
    obs = pct_df.groupby(keys)["percentile"].agg(["mean", "size"]).reset_index()
    ref = excl_df[keys + ["auroc", "n_eval", "n_active"]]
    both = obs.merge(ref, on=keys, how="inner")
    n_neg = both["n_eval"] - both["n_active"]
    expected = (both["auroc"] * n_neg + (both["n_active"] + 1) / 2) / both["n_eval"]
    worst = float((both["mean"] - expected).abs().max())
    if worst > tol:
        raise AssertionError(
            f"exclusivity percentiles do not reproduce the AUROC (worst deviation {worst:.2e} "
            f"> {tol:.0e}); the export is joined to the wrong compound set")
    return worst


def _load_exclusivity_task(euos_root, code, mode, inactives):
    """Subset actives (bin=1) + primary inactives (bin=0) for one exclusivity mode.

    ``mode`` in {"exclusive", "nonexclusive"}; reads
    ``06_subset_data/exclusivity/{code}_{mode}.csv`` (smiles, inchikey — actives only).
    Returns DataFrame[smiles, bin, inchikey] or None if the file is missing/empty.
    """
    path = os.path.join(euos_root, "06_subset_data", "exclusivity", f"{code}_{mode}.csv")
    if not os.path.exists(path):
        return None
    actives = pd.read_csv(path)
    if actives.empty:
        return None
    actives = actives.assign(bin=1)
    cols = ["smiles", "bin", "inchikey"]
    inactives = inactives[[c for c in cols if c in inactives.columns]]
    actives = actives[[c for c in cols if c in actives.columns]]
    return pd.concat([actives, inactives], ignore_index=True)


def run_individual_performance(pred_dir, euos_root, config, train_cache):
    """Per shared pathogen, look at ALL sub-model outputs (not just consensus_score).

    Returns (auroc_records, corr_mats):
      - auroc_records: per (pathogen, feature, set) AUROC + metrics on the pathogen's own
        EU OpenScreen assay — to see whether the ensemble members differ in quality.
      - corr_mats: {code: DataFrame} pairwise Spearman correlation between the sub-model
        score columns over the full prediction library — to see whether they rank compounds
        differently (label-free, so computed on all rows).
    """
    primary = pd.read_csv(os.path.join(euos_root, "primary_assays_manual.csv"))
    code_to_assay = dict(zip(primary["pathogen_code"], primary["assay_eos_id"]))
    name_by_code = dict(zip(config["code"], config["pathogen"]))
    eosid_by_code = dict(zip(config["code"], config["eosid"]))

    auroc_records, corr_mats = [], {}
    for code in SHARED_ORGANISMS:
        assay_id, eosid = code_to_assay.get(code), eosid_by_code.get(code)
        if assay_id is None or eosid is None:
            continue
        path = os.path.join(pred_dir, "euopenscreen", f"{eosid}.csv")
        if not os.path.exists(path):
            print(f"  [skip indiv] euopenscreen/{eosid}.csv not present yet")
            continue
        pred = pd.read_csv(path)
        if "input" not in pred.columns:
            continue
        pred = pred.rename(columns={"input": "smiles"})
        feature_cols = [c for c in pred.columns if c not in ("key", "smiles")]
        if len(feature_cols) < 2:  # nothing to compare
            continue

        # score-ranking agreement: Spearman corr between sub-model scores over the library
        corr_mats[code] = pred[feature_cols].corr(method="spearman").round(3)

        # per-feature AUROC on the pathogen's own assay (raw + dedup)
        lab = load_euos_primary(euos_root, code, assay_id)
        if lab is None:
            continue
        train_keys = train_cache.get(code)
        pathogen = name_by_code.get(code, code)
        for set_name, sub in merged_variants(pred, lab, train_keys):
            for feat in feature_cols:
                s = sub[["bin", feat]].dropna()
                if s["bin"].nunique() < 2:
                    continue
                m = compute_metrics(s["bin"].values, s[feat].values)
                rec = {"pathogen": pathogen, "code": code, "eosid": eosid,
                       "feature": feat, "set": set_name}
                rec.update({k: m[k] for k in METRIC_KEYS})
                auroc_records.append(rec)
    return auroc_records, corr_mats


def build_specificity_index(cross_df):
    """Per-model specificity index = own-organism AUROC - mean cross-organism AUROC.

    Uses the deduplicated cross matrix when available, else raw. Only defined for models
    whose organism has an own (diagonal) assay, i.e. the shared organisms.
    """
    if cross_df.empty:
        return pd.DataFrame(columns=["pathogen", "code", "same_pathogen_auroc",
                                     "mean_cross_auroc", "specificity_index"])
    use_set = "dedup" if "dedup" in set(cross_df["set"]) else "raw"
    df = cross_df[cross_df["set"] == use_set]
    rows = []
    for code, grp in df.groupby("model_code"):
        same = grp[grp["assay_code"] == code]["auroc"]
        cross = grp[grp["assay_code"] != code]["auroc"].dropna()
        if same.empty or same.isna().all():
            continue
        same_auroc = float(same.dropna().iloc[0]) if not same.dropna().empty else np.nan
        mean_cross = float(cross.mean()) if not cross.empty else np.nan
        if np.isnan(same_auroc) or np.isnan(mean_cross):
            continue
        rows.append({
            "pathogen": grp["model_pathogen"].iloc[0],
            "code": code,
            "same_pathogen_auroc": round(same_auroc, 4),
            "mean_cross_auroc": round(mean_cross, 4),
            "specificity_index": round(same_auroc - mean_cross, 4),
        })
    return pd.DataFrame(rows)


def build_leakage_report(config, euos_root, train_cache):
    """Overlap between each shared model's training InChIKeys and its EU OpenScreen assay."""
    primary = pd.read_csv(os.path.join(euos_root, "primary_assays_manual.csv"))
    code_to_assay = dict(zip(primary["pathogen_code"], primary["assay_eos_id"]))
    rows = []
    for _, row in config.iterrows():
        code, pathogen = row["code"], row["pathogen"]
        if code not in SHARED_ORGANISMS:
            continue
        train_keys = train_cache.get(code)
        assay_id = code_to_assay.get(code)
        lab = load_euos_primary(euos_root, code, assay_id) if assay_id else None
        n_train = len(train_keys) if train_keys else 0
        if lab is None:
            rows.append({"pathogen": pathogen, "code": code, "n_train": n_train,
                         "n_eval_conclusive": 0, "n_active": 0, "n_inactive": 0,
                         "n_overlap": 0, "n_overlap_active": 0, "n_overlap_inactive": 0})
            continue
        eval_keys = lab["inchikey"].dropna().astype(str)
        overlap = set(eval_keys) & (train_keys or set())
        ov = lab[lab["inchikey"].astype(str).isin(overlap)]
        rows.append({
            "pathogen": pathogen, "code": code, "n_train": n_train,
            "n_eval_conclusive": int(len(lab)),
            "n_active": int((lab["bin"] == 1).sum()),
            "n_inactive": int((lab["bin"] == 0).sum()),
            "n_overlap": int(len(overlap)),
            "n_overlap_active": int((ov["bin"] == 1).sum()),
            "n_overlap_inactive": int((ov["bin"] == 0).sum()),
        })
    return pd.DataFrame(rows)


def run_active_overlap(euos_root, config):
    """Pairwise active-compound overlap between the 7 EU OpenScreen primary assays.

    Label-only (no model): the interpretive backdrop for the cross-organism AUROCs — a compound
    active against one organism is often active against others, so high off-diagonal AUROC need
    not mean cross-organism prediction. Returns long-form rows
    (code_a, pathogen_a, code_b, pathogen_b, n_a, n_b, n_intersect, jaccard, containment).

    Two measures, because the active sets are very unequal in size (14 for P. aeruginosa vs 378
    for S. aureus) and they answer different questions:

    ``jaccard``     |A ∩ B| / |A ∪ B| — symmetric, but dominated by the larger set: P. aeruginosa
                    shares 13 of its 14 actives with A. baumannii and still scores only 0.22.
    ``containment`` |A ∩ B| / |A| — DIRECTIONAL: the share of row organism A's actives that are
                    also active against column organism B. The matrix is therefore NOT symmetric
                    (paeruginosa→abaumannii = 0.93, abaumannii→paeruginosa = 0.23), and the
                    diagonal is 1 by construction.

    Compounds are matched on the raw SMILES string (consistent with :func:`run_hit_promiscuity`;
    the leakage report matches on InChIKey instead).
    """
    primary = pd.read_csv(os.path.join(euos_root, "primary_assays_manual.csv"))
    code_to_assay = dict(zip(primary["pathogen_code"], primary["assay_eos_id"]))
    name_by_code = dict(zip(config["code"], config["pathogen"]))
    actives = {}
    for code in SHARED_ORGANISMS:
        assay_id = code_to_assay.get(code)
        if assay_id is None:
            continue
        lab = load_euos_primary(euos_root, code, assay_id)
        if lab is not None:
            actives[code] = set(lab.loc[lab["bin"] == 1, "smiles"])
    codes = [c for c in SHARED_ORGANISMS if c in actives]
    rows = []
    for a in codes:
        for b in codes:
            inter = len(actives[a] & actives[b])
            union = len(actives[a] | actives[b])
            rows.append({
                "code_a": a, "pathogen_a": name_by_code.get(a, a),
                "code_b": b, "pathogen_b": name_by_code.get(b, b),
                "n_a": len(actives[a]), "n_b": len(actives[b]),
                "n_intersect": inter,
                "jaccard": round(inter / union, 4) if union else 0.0,
                "containment": round(inter / len(actives[a]), 4) if actives[a] else 0.0,
            })
    return rows


def run_hit_promiscuity(euos_root, config):
    """Hit promiscuity — how many EU OpenScreen actives are hits in 1, 2, ... 7 pathogens.

    Label-only (no model), the per-compound counterpart of :func:`run_active_overlap`: instead of
    pairwise overlap it counts, for every distinct active, in how many of the 7 primary assays it
    is a hit (``bin == 1``). Compounds are matched on SMILES, consistent with the overlap
    analysis. ALL distinct actives are counted (the union across assays), including the few that
    were not tested in every assay — for those the count is a lower bound, so the number of such
    compounds is reported per bin as ``n_incomplete_coverage`` rather than dropped.

    Returns (dist_rows, hit_rows):
      - dist_rows: the aggregated distribution (n_pathogens, n_molecules, frac_molecules,
        n_molecules_ge, n_incomplete_coverage) — this is what the figure reads.
      - hit_rows: one row per distinct active (smiles, inchikey, n_pathogens, pathogens,
        n_assays_tested), most promiscuous first — the compound-level table.
    """
    primary = pd.read_csv(os.path.join(euos_root, "primary_assays_manual.csv"))
    code_to_assay = dict(zip(primary["pathogen_code"], primary["assay_eos_id"]))
    name_by_code = dict(zip(config["code"], config["pathogen"]))

    actives, tested, keys = {}, {}, {}
    for code in SHARED_ORGANISMS:
        assay_id = code_to_assay.get(code)
        if assay_id is None:
            continue
        lab = load_euos_primary(euos_root, code, assay_id)
        if lab is None:
            continue
        actives[code] = set(lab.loc[lab["bin"] == 1, "smiles"])
        tested[code] = set(lab["smiles"])                      # conclusive rows only
        keys.update(lab.dropna(subset=["inchikey"])
                       .set_index("smiles")["inchikey"].astype(str).to_dict())
    codes = [c for c in SHARED_ORGANISMS if c in actives]
    if not codes:
        print("  [skip] no primary assays available for hit promiscuity")
        return [], []

    n_assays = len(codes)
    hit_rows = []
    for smiles in sorted(set.union(*[actives[c] for c in codes])):
        hit_codes = [c for c in codes if smiles in actives[c]]
        hit_rows.append({
            "smiles": smiles,
            "inchikey": keys.get(smiles, ""),
            "n_pathogens": len(hit_codes),
            "pathogens": ";".join(name_by_code.get(c, c) for c in hit_codes),
            "n_assays_tested": sum(smiles in tested[c] for c in codes),
        })
    hit_rows.sort(key=lambda r: (-r["n_pathogens"], r["smiles"]))

    total = len(hit_rows)
    counts = {k: 0 for k in range(1, n_assays + 1)}
    incomplete = {k: 0 for k in range(1, n_assays + 1)}
    for r in hit_rows:
        counts[r["n_pathogens"]] += 1
        if r["n_assays_tested"] < n_assays:
            incomplete[r["n_pathogens"]] += 1
    dist_rows = []
    for k in range(1, n_assays + 1):
        dist_rows.append({
            "n_pathogens": k,
            "n_molecules": counts[k],
            "frac_molecules": round(counts[k] / total, 5) if total else 0.0,
            "n_molecules_ge": sum(counts[j] for j in range(k, n_assays + 1)),
            "n_incomplete_coverage": incomplete[k],
        })
    n_incomplete = sum(incomplete.values())
    print(f"  {total} distinct actives across {n_assays} primary assays; "
          f"{n_incomplete} not tested in all {n_assays} (their n_pathogens is a lower bound)")
    return dist_rows, hit_rows


def shared_training_inchikeys(models_root):
    """Union of the 7 shared-organism models' ChEMBL training InChIKeys, or None if unavailable.

    The leakage filter for the multi-model analyses (summed/max consensus, own-model rank). Those
    compare all 7 models against each other, so a compound is only leakage-free if NO model was
    trained on it — dropping the union is the analogue of :func:`eval_common.merged_variants`'
    per-model ``dedup`` for a single-model evaluation. Returns None (logged) when no training set is
    readable, in which case the filtered variants are skipped rather than silently equal to raw.
    """
    keys = set()
    for code in SHARED_ORGANISMS:
        k = training_inchikeys(models_root, code)
        if k:
            keys.update(k)
    if not keys:
        print("  [note] no ChEMBL training sets found — skipping the leakage-filtered variants")
        return None
    return keys


def euos_inchikeys(euos_root, config):
    """smiles -> InChIKey over the 7 primary assays, for leakage filtering by InChIKey."""
    primary = pd.read_csv(os.path.join(euos_root, "primary_assays_manual.csv"))
    code_to_assay = dict(zip(primary["pathogen_code"], primary["assay_eos_id"]))
    keys = {}
    for code in SHARED_ORGANISMS:
        assay_id = code_to_assay.get(code)
        lab = load_euos_primary(euos_root, code, assay_id) if assay_id else None
        if lab is None:
            continue
        keys.update(lab.dropna(subset=["inchikey"])
                       .drop_duplicates("smiles")
                       .set_index("smiles")["inchikey"].astype(str).to_dict())
    return keys


def _drop_leaked(df, smiles_to_key, train_keys, label):
    """Drop rows whose InChIKey is in ``train_keys`` (no-op when either input is missing).

    A compound with no InChIKey cannot be checked, so it is KEPT and counted in the log — never
    dropped on the basis of missing information.
    """
    if not train_keys or not smiles_to_key:
        return df, 0
    keyed = pd.Series(df.index.map(lambda s: smiles_to_key.get(s)), index=df.index)
    leaked = keyed.isin(train_keys)
    n_unknown = int(keyed.isna().sum())
    if n_unknown:
        print(f"  [note] {n_unknown} {label} compounds have no InChIKey — kept (cannot be checked)")
    print(f"  [dedup] dropped {int(leaked.sum())} of {len(df)} {label} compounds present in a "
          "ChEMBL training set")
    return df[~leaked], int(leaked.sum())


def _consensus_scores_with_hit_counts(pred_dir, euos_root, config, agg, normalize=False,
                                      rerank=False, train_keys=None, smiles_to_key=None):
    """Per-compound aggregated consensus score + primary-assay hit counts.

    The shared core of the score-distribution analyses. Aggregates the 7 shared-organism models'
    ``consensus_score`` per compound with ``agg`` ("sum" → total in [0, 7]; "max" → the single
    most confident model's score in [0, 1]), and counts, from the PRIMARY assays only, how many
    gave a conclusive result (``n_assays_tested``) and how many called the compound a hit
    (``n_pathogens``). No leakage filtering — training-set compounds are deliberately kept.

    With ``normalize=True`` each model's score is first replaced by its percentile within that
    model's OWN library distribution (over all scored compounds, before any label filtering), which
    puts the 7 models on a common scale before aggregating — the models' raw scores are not
    calibrated to each other (library medians range from ~0.35 to ~0.61), so a raw ``max`` partly
    reflects which model happens to output the highest values. Same transform as the ``percentile``
    ranking in :func:`run_exclusive_hit_model_rank`.

    With ``rerank=True`` the aggregated value is finally re-expressed as its own percentile within
    the library. This matters for ``agg="max"``: the maximum of 7 values is high by construction
    (for 7 independent uniforms the expected maximum is 7/8), so an un-re-ranked max-percentile axis
    is compressed into its top decile and its baseline is not interpretable. Re-ranking is a
    monotone transform — it leaves every ranking metric (AUROC, EF, ...) bit-identical and only
    makes the axis readable: "this compound's best-model percentile beats that fraction of the
    library's best-model percentiles".

    Only compounds scored by all 7 models are returned, so the aggregate is comparable across
    compounds; any incomplete ones are counted and logged. Returns (DataFrame indexed by smiles
    with [score, n_pathogens, n_assays_tested], n_unlabelled) where ``n_unlabelled`` counts the
    compounds with no conclusive primary result at all (already removed — they cannot be
    classified), or (None, 0) if any model's predictions are missing.
    """
    primary = pd.read_csv(os.path.join(euos_root, "primary_assays_manual.csv"))
    code_to_assay = dict(zip(primary["pathogen_code"], primary["assay_eos_id"]))
    eosid_by_code = dict(zip(config["code"], config["eosid"]))

    scores = {}
    for code in SHARED_ORGANISMS:
        eosid = eosid_by_code.get(code)
        pred = load_predictions(pred_dir, "euopenscreen", eosid) if eosid else None
        if pred is None:
            print(f"  [skip] no predictions for {code} — cannot aggregate over all 7 models")
            return None, 0
        scores[code] = pred.drop_duplicates("smiles").set_index("smiles")["score"]
    mat = pd.DataFrame(scores)
    complete = mat.notna().all(axis=1)
    if not complete.all():
        print(f"  [note] {int((~complete).sum())} compounds lack a score from all 7 models "
              f"— excluded ({agg} would not be comparable across compounds)")
    mat = mat[complete]
    if normalize:
        mat = mat.rank(pct=True)          # within-model percentile, over the whole library
    score = mat.sum(axis=1) if agg == "sum" else mat.max(axis=1)
    if rerank:
        score = score.rank(pct=True)      # library percentile OF the aggregate

    n_active = pd.Series(0, index=score.index, dtype=int)
    n_tested = pd.Series(0, index=score.index, dtype=int)
    for code in SHARED_ORGANISMS:
        assay_id = code_to_assay.get(code)
        lab = load_euos_primary(euos_root, code, assay_id) if assay_id else None
        if lab is None:
            continue
        # `load_euos_primary` already returns one row per compound, resolved by max bin. This used to
        # carry its own `drop_duplicates(keep="first")` — a SECOND, row-order-dependent rule that
        # agreed with the first only by luck (the active row precedes the inactive one in the
        # A. baumannii file). No guard is needed in its place: the `reindex` below raises outright on
        # a duplicate index, so a regression in the loader fails loudly instead of double-counting.
        lab = lab.set_index("smiles")["bin"]
        aligned = lab.reindex(score.index)
        n_tested += aligned.notna().astype(int)
        n_active += (aligned == 1).astype(int)

    df = pd.DataFrame({"score": score, "n_pathogens": n_active, "n_assays_tested": n_tested})
    unlabelled = int((df["n_assays_tested"] == 0).sum())
    df = df[df["n_assays_tested"] > 0]  # no conclusive primary result → unclassifiable
    leak_note = "no leakage filtering"
    if train_keys:
        df, n_leaked = _drop_leaked(df, smiles_to_key, train_keys, "classified")
        leak_note = f"{n_leaked} training-set compounds removed"
    print(f"  {len(df)} classified compounds ({unlabelled} with no conclusive primary result "
          f"are unlabelled); {agg} over 7 models, {leak_note}")
    return df, unlabelled


def _class_box_stats(df, classes, score_col, unlabelled, n_models, score_scale="raw",
                     leakage="raw"):
    """Box statistics (n, median, q1, q3, 1.5xIQR whiskers, mean, min, max) per hit class.

    Only summary numbers leave this function — the inactive class has ~10^5 compounds and is never
    shipped per-molecule to the figures.
    """
    rows = []
    for hit_class in classes:
        v = df.loc[df["hit_class"] == hit_class, score_col]
        if v.empty:
            continue
        q1, med, q3 = (float(v.quantile(q)) for q in (0.25, 0.5, 0.75))
        iqr = q3 - q1
        inside = v[(v >= q1 - 1.5 * iqr) & (v <= q3 + 1.5 * iqr)]
        rows.append({
            "hit_class": hit_class, "n": int(len(v)),
            "median": round(med, 4), "q1": round(q1, 4), "q3": round(q3, 4),
            "whisker_lo": round(float(inside.min()), 4),
            "whisker_hi": round(float(inside.max()), 4),
            "mean": round(float(v.mean()), 4),
            "min": round(float(v.min()), 4), "max": round(float(v.max()), 4),
            "n_models_aggregated": n_models, "n_unlabelled_excluded": unlabelled,
            "score_scale": score_scale, "leakage": leakage,
        })
    return rows


def _active_rows(df, score_col):
    """Per-compound rows for the active classes only (a few hundred), for the jitter overlay."""
    actives = df[df["hit_class"] != "inactive"].copy()
    actives = actives.reset_index().rename(columns={"index": "smiles"})
    return actives.sort_values(["n_pathogens", score_col], ascending=[False, False])[
        ["smiles", "hit_class", "n_pathogens", "n_assays_tested", score_col]
    ].round({score_col: 4}).to_dict("records")


def run_consensus_sum_by_hit_class(pred_dir, euos_root, config,
                                   train_keys=None, smiles_to_key=None):
    """Summed consensus score across the 7 pathogen models, split by EU OpenScreen hit class.

    Each compound's 7 model scores are summed (one score in [0, 7]) and the compound is assigned
    to one of four classes from the PRIMARY assay labels only:

      - ``inactive``  — ``bin == 0`` in every primary assay where it has a conclusive result,
                        and active in none (compounds tested in only some assays still count,
                        per the definition "inactive across all tested datasets").
      - ``exclusive`` — a hit in exactly 1 of the 7 primary assays.
      - ``narrow``    — a hit in 2-3 primary assays.
      - ``broad``     — a hit in more than 3 primary assays (4-7), i.e. broad-spectrum.

    The last two are the split of the "shared" (non-exclusive) hits used by analysis 3. Leakage
    filtering is opt-in, exactly as for :func:`run_consensus_max_by_activity`: pass ``train_keys``
    + ``smiles_to_key`` to drop every compound in ANY of the 7 models' ChEMBL training sets, from
    all four classes alike. Omit them for the raw variant, which deliberately keeps training-set
    compounds and so describes the score distribution rather than out-of-sample performance.

    Note the dedup twin's class balance is NOT the raw one rescaled: leakage rises steeply with
    promiscuity (17% of exclusive hits are in a training set — 68 of 390 — against 40% of narrow,
    61 of 153, and 69% of broad, 31 of 45), so the broad class thins out far more than the
    exclusive one.

    Returns (stat_rows, active_rows): the per-class box statistics the figure reads, and the
    per-compound summed score for the three active classes.
    """
    df, unlabelled = _consensus_scores_with_hit_counts(
        pred_dir, euos_root, config, agg="sum",
        train_keys=train_keys, smiles_to_key=smiles_to_key)
    if df is None:
        return [], []
    df = df.rename(columns={"score": "consensus_sum"})
    df["hit_class"] = np.select(
        [df["n_pathogens"] > NARROW_MAX_PATHOGENS,      # 4-7 pathogens
         df["n_pathogens"] >= 2,                        # 2-3 pathogens
         df["n_pathogens"] == 1],
        ["broad", "narrow", "exclusive"], default="inactive")
    return (_class_box_stats(df, HIT_CLASSES, "consensus_sum", unlabelled,
                             len(SHARED_ORGANISMS), score_scale="raw",
                             leakage="dedup" if train_keys else "raw"),
            _active_rows(df, "consensus_sum"))


def run_consensus_max_by_activity(pred_dir, euos_root, config, normalize=False,
                                  train_keys=None, smiles_to_key=None):
    """MAXIMUM consensus score across the 7 pathogen models, active vs inactive.

    The simplest read of the same data: per compound, the single highest of the 7 models' scores —
    i.e. how confident the most confident model is — split into just two classes from the PRIMARY
    assays:

      - ``inactive`` — ``bin == 0`` in every primary assay where it has a conclusive result.
      - ``active``   — a hit in one or more of the 7 primary assays, regardless of how many.

    ``normalize=False`` takes the max of the raw ``consensus_score`` values (in [0, 1]);
    ``normalize=True`` first converts each model's score to its within-model library percentile (so
    the max is not biased towards whichever model outputs the highest values) AND re-ranks the
    resulting maximum over the library, which removes the best-of-7 baseline from the axis without
    changing any ordering — see :func:`_consensus_scores_with_hit_counts`. Both variants are
    written, as two panels; the normalised one is the readable version of the same content.

    Same conventions as the summed variant: primary assays only, compounds with no conclusive primary
    result reported as unlabelled. Leakage filtering is opt-in: pass ``train_keys`` (from
    :func:`shared_training_inchikeys`) + ``smiles_to_key`` (from :func:`euos_inchikeys`) to drop
    every compound present in ANY of the 7 models' ChEMBL training sets. BOTH classes are filtered
    the same way — filtering only the actives would strip their trained-on (highest-scoring) members
    while leaving the inactive class's leaked compounds in, biasing the gap downward. Omit them for
    the raw variant, which deliberately keeps training-set compounds. The ``leakage`` column of the
    returned stats records which was used.

    Returns (stat_rows, active_rows, roc_rows): the per-class box statistics the figure reads, the
    per-compound maximum score for the active class (which keeps its ``n_pathogens`` count), and
    the thinned ROC curve of active-vs-inactive under this score.
    """
    df, unlabelled = _consensus_scores_with_hit_counts(
        pred_dir, euos_root, config, agg="max", normalize=normalize, rerank=normalize,
        train_keys=train_keys, smiles_to_key=smiles_to_key)
    if df is None:
        return [], []
    df = df.rename(columns={"score": "consensus_max"})
    df["hit_class"] = np.where(df["n_pathogens"] >= 1, "active", "inactive")
    scale = "percentile_reranked" if normalize else "raw"
    stats = _class_box_stats(df, ACTIVITY_CLASSES, "consensus_max", unlabelled,
                             len(SHARED_ORGANISMS), score_scale=scale,
                             leakage="dedup" if train_keys else "raw")

    # Separation as ONE number, on every stats row so the figure can annotate it without a join.
    #
    # It has to be computed here and cannot be recovered downstream. The obvious shortcut — AUROC
    # is the mean rank percentile of the positives — does NOT hold for this score, because the
    # percentile is ranked over the whole scored library while the AUROC is over the smaller
    # evaluated subset (unlabelled compounds removed, then training-set compounds removed). The two
    # sets differ non-randomly: dedup preferentially drops high scorers, since promiscuous hits are
    # the ones ChEMBL trained on. Deriving it from the shipped actives that way gives 0.7511 against
    # a true 0.7574. The inactive scores are never shipped per-molecule (~10^5 rows), so this column
    # is the only route to the number.
    y = (df["hit_class"] == "active").to_numpy(dtype=int)
    score = df["consensus_max"].to_numpy()
    two_class = 0 < y.sum() < len(y)
    auroc = round(float(roc_auc_score(y, score)), 4) if two_class else np.nan
    for row in stats:
        row["auroc"] = auroc

    # The curve behind that number. Same construction as the per-organism ROCs
    # (:func:`_roc_records`): sklearn's compact staircase, then :func:`_thin_curve` to cap the
    # vertex count while keeping every TPR jump, so the shape survives at ~10^5 negatives against
    # a few hundred positives. Computed here for the same reason the AUROC is — the inactive
    # scores never leave this function.
    roc_rows = []
    if two_class:
        fpr, tpr, _ = roc_curve(y, score)
        fpr, tpr = _thin_curve(fpr, tpr)
        n_pos, n_neg = int(y.sum()), int((y == 0).sum())
        roc_rows = [{"score_scale": scale, "leakage": "dedup" if train_keys else "raw",
                     "fpr": round(float(x), 5), "tpr": round(float(t), 5),
                     "n_pos": n_pos, "n_neg": n_neg, "auroc": auroc}
                    for x, t in zip(fpr, tpr)]
    return stats, _active_rows(df, "consensus_max"), roc_rows


def run_exclusive_hit_model_rank(pred_dir, euos_root, config, train_keys=None,
                                 smiles_to_key=None):
    """For each EXCLUSIVE hit, where does its own pathogen's model rank it among the 7 models?

    An exclusive hit is active in exactly 1 of the 7 primary assays, so exactly one of the 7 models
    is the "right" one. The 7 models' scores for that compound are ranked best-first and we record
    the position of the right model: rank 1 = the compound's own pathogen scores it highest, rank 2
    = one other pathogen's model ranks it above its own, and so on to rank 7. The figure histograms
    those ranks (chance level = n_hits / 7).

    Two rankings are computed, because the models' raw ``consensus_score`` values are NOT calibrated
    to a common scale (per-library medians range from ~0.35 to ~0.61), so a raw ranking partly
    reflects each model's output offset rather than the compound:

      - ``raw``        — rank the raw consensus scores, as-is.
      - ``percentile`` — rank each model's score converted to its percentile within that model's own
                         library distribution first, which puts the 7 models on a common footing.

    Ties are resolved with ``method="min"`` (the best-case rank for the true pathogen) and counted.
    Primary assays only. Leakage filtering is opt-in: pass ``train_keys`` (from
    :func:`shared_training_inchikeys`) + ``smiles_to_key`` (from :func:`euos_inchikeys`) to drop the
    exclusive hits present in ANY of the 7 models' ChEMBL training sets. Any-model rather than
    own-model filtering, because leakage in a NON-own model inflates that model's score and pushes
    the own model down the ranking — a spurious effect in the opposite direction; dropping
    any-leaked compounds leaves all 7 scores leakage-free for every compound that remains.

    Returns (rank_rows, compound_rows):
      - rank_rows: per (ranking, rank) the number of compounds, plus fraction, cumulative count and
        the chance baseline — what the figure reads. Also carries one ``n_<code>`` column per
        organism with that rank bin's breakdown by the hit's own pathogen, so the bars can be
        stacked by pathogen without touching per-compound data.
      - compound_rows: per compound and ranking, its pathogen, the rank of its own model and the
        pathogen whose model ranked it top (a few hundred rows per ranking).
    """
    primary = pd.read_csv(os.path.join(euos_root, "primary_assays_manual.csv"))
    code_to_assay = dict(zip(primary["pathogen_code"], primary["assay_eos_id"]))
    name_by_code = dict(zip(config["code"], config["pathogen"]))
    eosid_by_code = dict(zip(config["code"], config["eosid"]))

    scores = {}
    for code in SHARED_ORGANISMS:
        eosid = eosid_by_code.get(code)
        pred = load_predictions(pred_dir, "euopenscreen", eosid) if eosid else None
        if pred is None:
            print(f"  [skip] no predictions for {code} — cannot rank the 7 models")
            return [], []
        scores[code] = pred.drop_duplicates("smiles").set_index("smiles")["score"]
    mat = pd.DataFrame(scores)
    mat = mat[mat.notna().all(axis=1)]

    # Exclusive hits: active in exactly one primary assay → that assay's code is the true pathogen.
    hit_codes = {}
    for code in SHARED_ORGANISMS:
        assay_id = code_to_assay.get(code)
        lab = load_euos_primary(euos_root, code, assay_id) if assay_id else None
        if lab is None:
            continue
        for smiles in lab.loc[lab["bin"] == 1, "smiles"]:
            hit_codes.setdefault(smiles, []).append(code)
    truth = {s: c[0] for s, c in hit_codes.items() if len(c) == 1 and s in mat.index}
    if train_keys and smiles_to_key:
        leaked = {s for s in truth if smiles_to_key.get(s) in train_keys}
        unknown = sum(1 for s in truth if s not in smiles_to_key)
        if unknown:
            print(f"  [note] {unknown} exclusive hits have no InChIKey — kept (cannot be checked)")
        print(f"  [dedup] dropped {len(leaked)} of {len(truth)} exclusive hits present in a "
              "ChEMBL training set (any of the 7 models)")
        truth = {s: c for s, c in truth.items() if s not in leaked}
        if not truth:
            print("  [skip] no exclusive hits left after leakage filtering")
            return [], []
    n_missing = sum(1 for s, c in hit_codes.items() if len(c) == 1 and s not in mat.index)
    if n_missing:
        print(f"  [note] {n_missing} exclusive hits have no complete set of 7 model scores "
              "— excluded from the ranking")
    if not truth:
        print("  [skip] no exclusive hits with predictions from all 7 models")
        return [], []
    order = sorted(truth)
    n_models = mat.shape[1]

    rank_rows, compound_rows = [], []
    for ranking in ("raw", "percentile"):
        # percentile: each model's score → its rank within that model's own library distribution
        m = mat if ranking == "raw" else mat.rank(pct=True)
        sub = m.loc[order]
        ranks = sub.rank(axis=1, ascending=False, method="min")
        n_ties = int((sub.apply(lambda r: r.duplicated().any(), axis=1)).sum())
        true_rank = pd.Series([int(ranks.at[s, truth[s]]) for s in order], index=order)
        top_code = sub.idxmax(axis=1)
        true_code = pd.Series([truth[s] for s in order], index=order)
        counts = true_rank.value_counts()
        # rank x own-pathogen breakdown, so the bars can be stacked by pathogen
        by_pathogen = pd.crosstab(true_rank, true_code)
        total = len(true_rank)
        for k in range(1, n_models + 1):
            n = int(counts.get(k, 0))
            row = {
                "ranking": ranking, "rank": k, "n_molecules": n,
                "frac_molecules": round(n / total, 5) if total else 0.0,
                "n_molecules_le": int(sum(int(counts.get(j, 0)) for j in range(1, k + 1))),
                "n_total": total, "n_chance": round(total / n_models, 2),
                "n_models_ranked": n_models, "n_tied_scores": n_ties,
                "leakage": "dedup" if train_keys else "raw",
            }
            for code in SHARED_ORGANISMS:
                row[f"n_{code}"] = int(by_pathogen.at[k, code]) \
                    if k in by_pathogen.index and code in by_pathogen.columns else 0
            rank_rows.append(row)
        for s in order:
            compound_rows.append({
                "smiles": s, "ranking": ranking,
                "pathogen": name_by_code.get(truth[s], truth[s]), "code": truth[s],
                "true_pathogen_rank": int(true_rank[s]),
                "top_ranked_pathogen": name_by_code.get(top_code[s], top_code[s]),
                "leakage": "dedup" if train_keys else "raw",
            })
        print(f"  {ranking}: rank-1 = {int(counts.get(1, 0))}/{total} "
              f"(chance {total / n_models:.1f}); {n_ties} compounds had tied scores")
    return rank_rows, compound_rows


# --------------------------------------------------------------------------- #
# Orchestration                                                               #
# --------------------------------------------------------------------------- #
def _write_csvs(frames, target_dir, label):
    """Write ``{filename: DataFrame}`` into ``target_dir``, returning them keyed by relative path."""
    os.makedirs(target_dir, exist_ok=True)
    out = {}
    for fname, df in frames.items():
        df.to_csv(os.path.join(target_dir, fname), index=False)
        rel = os.path.join(label, fname) if label else fname
        out[rel] = df
        print(f"  wrote {rel}: {len(df)} rows")
    return out


def _filter_set(df, set_name, fname):
    """One variant of a long-form metric table (``set`` column of ``raw`` / ``dedup``).

    A frame with no ``set`` column is returned whole (it has no leakage dimension to filter on).
    An empty result is written anyway and logged — an empty ``deduplicated/`` table is the honest
    signal that no training sets were available, not something to hide by falling back to raw.
    """
    if df.empty or "set" not in df.columns:
        return df
    sub = df[df["set"] == set_name]
    if sub.empty:
        print(f"  [note] no '{set_name}' rows for {fname} — writing an empty table")
    return sub


def run_all(pred_dir, euos_root, models_root, config_path, output_dir):
    """Run the EU OpenScreen analyses and write the summary CSVs into ``output_dir``.

    Outputs are split three ways so the leakage status of every number is unambiguous: the FULL
    analysis (training-set compounds kept) under ``full/``, the DEDUPLICATED one (leakage-filtered)
    under ``deduplicated/``, and anything with no leakage dimension — label-only tables, and the
    leakage audit itself — at the top level. ``individual_performance/`` and (in the caller)
    ``eos3dys_validation/`` are separate analysis families and keep their own layout.

    Returns a dict of the written DataFrames keyed by path relative to ``output_dir`` (also useful
    for testing). Prints a short coverage summary; nothing here mutates git-tracked or step-04 data.
    """
    config = pd.read_csv(config_path)  # columns: pathogen, code, eosid

    # Cache training InChIKeys per pathogen once (dedup source; may be None if absent).
    train_cache = {code: training_inchikeys(models_root, code) for code in config["code"]}
    if all(v is None for v in train_cache.values()):
        print("[leakage] no training sets found under "
              f"{os.path.join(models_root, 'output', '07_datasets')} — reporting RAW only")

    print("[EU OpenScreen] evaluating own-assay, exclusivity and cross-organism ...")
    own, excl, cross, roc, excl_pct = run_euopenscreen(pred_dir, euos_root, config, train_cache)
    print("[EU OpenScreen] active-set overlap ...")
    overlap = run_active_overlap(euos_root, config)
    print("[EU OpenScreen] hit promiscuity (actives shared across pathogens) ...")
    prom_dist, prom_hits = run_hit_promiscuity(euos_root, config)
    print("[EU OpenScreen] summed consensus score by hit class ...")
    sum_stats, sum_actives = run_consensus_sum_by_hit_class(pred_dir, euos_root, config)
    print("[EU OpenScreen] maximum consensus score, active vs inactive ...")
    max_stats, max_actives, max_roc = run_consensus_max_by_activity(
        pred_dir, euos_root, config)
    print("[EU OpenScreen] maximum within-model percentile, active vs inactive ...")
    maxp_stats, maxp_actives, maxp_roc = run_consensus_max_by_activity(
        pred_dir, euos_root, config, normalize=True)
    print("[EU OpenScreen] own-model rank for exclusive hits ...")
    rank_dist, rank_compounds = run_exclusive_hit_model_rank(pred_dir, euos_root, config)

    # Leakage-filtered twins of the two multi-model figures: every compound present in ANY of the 7
    # models' ChEMBL training sets is dropped (both classes, so the comparison stays like-for-like).
    print("[EU OpenScreen] leakage-filtered (dedup) twins ...")
    shared_train = shared_training_inchikeys(models_root)
    key_map = euos_inchikeys(euos_root, config) if shared_train else None
    maxd_stats, maxd_actives, maxd_roc = [], [], []
    sumd_stats, sumd_actives = [], []
    rankd_dist, rankd_compounds = [], []
    if shared_train:
        sumd_stats, sumd_actives = run_consensus_sum_by_hit_class(
            pred_dir, euos_root, config, train_keys=shared_train, smiles_to_key=key_map)
        maxd_stats, maxd_actives, maxd_roc = run_consensus_max_by_activity(
            pred_dir, euos_root, config, normalize=True,
            train_keys=shared_train, smiles_to_key=key_map)
        rankd_dist, rankd_compounds = run_exclusive_hit_model_rank(
            pred_dir, euos_root, config, train_keys=shared_train, smiles_to_key=key_map)

    own_df = pd.DataFrame(own)
    excl_df = pd.DataFrame(excl)
    excl_pct_df = pd.DataFrame(excl_pct)
    # The exported distribution must integrate back to the published AUROC — see the function.
    worst = _check_percentiles_match_auroc(excl_pct_df, excl_df)
    if worst is not None:
        print(f"[EU OpenScreen] exclusivity percentiles reproduce AUROC to {worst:.2e} "
              f"({len(excl_pct_df)} actives)")
    cross_df = pd.DataFrame(cross)
    roc_df = pd.DataFrame(roc)
    spec_df = build_specificity_index(cross_df)
    leak_df = build_leakage_report(config, euos_root, train_cache)
    overlap_df = pd.DataFrame(overlap)

    # Outputs are grouped by whether they describe the FULL (training-set compounds kept) or the
    # DEDUPLICATED (leakage-filtered) analysis, so the two can never be confused:
    #   <output_dir>/                — no leakage dimension at all (label-only, or about leakage)
    #   <output_dir>/full/           — training-set compounds INCLUDED
    #   <output_dir>/deduplicated/   — training-set compounds REMOVED
    # The long-form metric tables carry both variants in a ``set`` column, so they are written into
    # BOTH subfolders, each filtered to that folder's variant (see SPLIT_BY_SET).
    main_outputs = {
        "05_leakage_report.csv": leak_df,          # the leakage audit itself
        "05_active_overlap.csv": overlap_df,       # label-only, no model involved
        "05_hit_promiscuity.csv": pd.DataFrame(prom_dist),
        "05_promiscuous_hits.csv": pd.DataFrame(prom_hits),
    }
    split_outputs = {                              # both variants, filtered per folder
        "05_euopenscreen_auroc.csv": own_df,
        "05_euopenscreen_roc.csv": roc_df,
        "05_hit_exclusivity.csv": excl_df,
        "05_hit_exclusivity_percentiles.csv": excl_pct_df,
        "05_cross_organism_euos.csv": cross_df,
    }
    full_outputs = {
        "05_consensus_sum_boxstats.csv": pd.DataFrame(sum_stats),
        "05_consensus_sum_actives.csv": pd.DataFrame(sum_actives),
        "05_consensus_max_boxstats.csv": pd.DataFrame(max_stats),
        "05_consensus_max_actives.csv": pd.DataFrame(max_actives),
        "05_consensus_max_percentile_boxstats.csv": pd.DataFrame(maxp_stats),
        "05_consensus_max_percentile_actives.csv": pd.DataFrame(maxp_actives),
        "05_consensus_max_roc.csv": pd.DataFrame(max_roc),
        "05_consensus_max_percentile_roc.csv": pd.DataFrame(maxp_roc),
        "05_exclusive_hit_model_rank.csv": pd.DataFrame(rank_dist),
        "05_exclusive_hit_model_rank_compounds.csv": pd.DataFrame(rank_compounds),
    }
    dedup_outputs = {
        # built from the dedup cross matrix when available (see build_specificity_index)
        "05_specificity_index.csv": spec_df,
        "05_consensus_sum_dedup_boxstats.csv": pd.DataFrame(sumd_stats),
        "05_consensus_sum_dedup_actives.csv": pd.DataFrame(sumd_actives),
        "05_consensus_max_percentile_dedup_boxstats.csv": pd.DataFrame(maxd_stats),
        "05_consensus_max_percentile_dedup_roc.csv": pd.DataFrame(maxd_roc),
        "05_consensus_max_percentile_dedup_actives.csv": pd.DataFrame(maxd_actives),
        "05_exclusive_hit_model_rank_dedup.csv": pd.DataFrame(rankd_dist),
        "05_exclusive_hit_model_rank_dedup_compounds.csv": pd.DataFrame(rankd_compounds),
    }

    outputs = {}
    outputs.update(_write_csvs(main_outputs, output_dir, ""))
    outputs.update(_write_csvs(full_outputs, os.path.join(output_dir, FULL_SUBDIR), FULL_SUBDIR))
    outputs.update(_write_csvs(dedup_outputs, os.path.join(output_dir, DEDUP_SUBDIR), DEDUP_SUBDIR))
    for subdir, set_name in ((FULL_SUBDIR, "raw"), (DEDUP_SUBDIR, "dedup")):
        filtered = {f: _filter_set(df, set_name, f) for f, df in split_outputs.items()}
        outputs.update(_write_csvs(filtered, os.path.join(output_dir, subdir), subdir))

    # Individual sub-model performance, per shared pathogen → subfolder.
    print("[individual] per-pathogen sub-model AUROC + score correlations ...")
    indiv_dir = os.path.join(output_dir, "individual_performance")
    os.makedirs(indiv_dir, exist_ok=True)
    ip_auroc, ip_corr = run_individual_performance(pred_dir, euos_root, config, train_cache)
    pd.DataFrame(ip_auroc).to_csv(
        os.path.join(indiv_dir, "05_submodel_auroc.csv"), index=False)
    for code, cm in ip_corr.items():
        cm.to_csv(os.path.join(indiv_dir, f"{code}_submodel_corr.csv"))
    print(f"  wrote individual_performance/ ({len(ip_auroc)} auroc rows, "
          f"{len(ip_corr)} corr matrices)")

    return outputs
