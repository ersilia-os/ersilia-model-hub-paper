"""CoAdd-model (eos3dys) validation on EU OpenScreen — the ``eos3dys_validation`` sub-analysis
of step 05.

The CoAdd model ``eos3dys`` emits many endpoint columns (``{organism}_{strain}_{inhib_50|mic_25}``
for its reference strains, plus ``cytotoxicity_ic50`` and ``hemolitic_activity``). This module
scores **each endpoint** against **each EU OpenScreen primary assay** (the 7 shared organisms),
so the diagonal is "this endpoint's organism vs its own EU OpenScreen assay" and the off-diagonal
is cross-organism — a generalization + organism-specificity test for the CoAdd model.

Alongside that endpoint x assay matrix, two own-assay analyses mirror what the ChEMBL models get on
EU OpenScreen, for the organisms eos3dys and EU OpenScreen BOTH cover (see :func:`matched_endpoints`
— *E. faecium* has an EU OpenScreen assay but no eos3dys endpoint, so 6 of the 7 qualify):

  - training-set overlap (:func:`build_eos3dys_overlap_report`) — how much of the evaluation set,
    and of its actives, eos3dys already saw in CoAdd training; the twin of the ChEMBL leakage report.
  - exclusive vs shared hit AUROC (:func:`run_eos3dys_exclusivity`) — the twin of analysis 3, using
    the same upstream exclusivity subsets so the two figures compare directly.

This is the *opposite direction* from :mod:`eval_coadd` (which scores our ChEMBL models against
CoAdd labels). Metrics are reported ``raw`` + ``dedup``; here dedup removes compounds in the
endpoint's own CoAdd training set (its binarised strain file), i.e. compounds the model has seen.
CoAdd is used ONLY as that dedup source — the labels are always EU OpenScreen, which is the
out-of-sample direction for this model. Reuses the shared primitives in :mod:`eval_common` and the
EU OpenScreen loaders in :mod:`eval_euopenscreen`.
"""

import os

import numpy as np
import pandas as pd

from default import ACTIVITY_CLASSES, COADD_MODEL_ID, COADD_REF_STRAINS, SHARED_ORGANISMS
from eval_common import evaluate
from eval_euopenscreen import _load_exclusivity_task, load_euos_primary

# eos3dys endpoint metric suffix -> CoAdd binarised subdir holding that endpoint's labels/training set.
_METRIC_SUBDIR = {
    "inhib_50": "03_binarised_inhibition",
    "mic_25": "05_binarised_mic",
}


def _endpoint_to_strain(endpoint):
    """Parse an eos3dys endpoint into ``(strain, subdir)`` for its CoAdd binarised file.

    ``abaumannii_ATCC19606_inhib_50`` -> ``("abaumannii_ATCC19606", "03_binarised_inhibition")``.
    Returns None for endpoints with no organism strain (``cytotoxicity_ic50``,
    ``hemolitic_activity``) or an unrecognised metric suffix.
    """
    for metric, subdir in _METRIC_SUBDIR.items():
        suffix = "_" + metric
        if endpoint.endswith(suffix):
            return endpoint[: -len(suffix)], subdir
    return None


def _endpoint_organism(endpoint):
    """Organism code an endpoint belongs to (the first ``_``-delimited token)."""
    return endpoint.split("_")[0]


def coadd_training_inchikeys(coadd_root, endpoint):
    """InChIKeys the CoAdd model saw for this endpoint (its binarised strain file).

    Used to remove already-seen compounds from the EU OpenScreen evaluation. Returns a set, or
    None (dedup becomes a no-op for this endpoint) when the endpoint has no strain file.
    """
    parsed = _endpoint_to_strain(endpoint)
    if parsed is None:
        return None
    strain, subdir = parsed
    path = os.path.join(coadd_root, subdir, f"{strain}.csv")
    if not os.path.exists(path):
        return None
    try:
        keys = pd.read_csv(path, usecols=["inchikey"])["inchikey"]
    except (ValueError, KeyError):
        return None
    keys = set(keys.dropna().astype(str))
    return keys or None


def run_eos3dys_euopenscreen(pred_dir, euos_root, coadd_root, config, eos_id=COADD_MODEL_ID):
    """Every eos3dys endpoint scored against every EU OpenScreen primary assay (raw + dedup).

    Returns a list of metric records with fields: endpoint, endpoint_organism, assay_code,
    assay_pathogen, same_organism, set, + metrics. ``config`` is the pathogen table (for the
    assay pathogen names).
    """
    primary = pd.read_csv(os.path.join(euos_root, "primary_assays_manual.csv"))
    code_to_assay = dict(zip(primary["pathogen_code"], primary["assay_eos_id"]))
    name_by_code = dict(zip(config["code"], config["pathogen"]))

    # The 7 shared EU OpenScreen primary assays = the matrix columns.
    assays = {}
    for code in SHARED_ORGANISMS:
        assay_id = code_to_assay.get(code)
        if assay_id is None:
            continue
        lab = load_euos_primary(euos_root, code, assay_id)
        if lab is not None:
            assays[code] = lab

    path = os.path.join(pred_dir, "euopenscreen", f"{eos_id}.csv")
    if not os.path.exists(path):
        print(f"  [skip] euopenscreen/{eos_id}.csv not present")
        return []
    pred = pd.read_csv(path)
    if "input" not in pred.columns:
        print(f"  [skip] euopenscreen/{eos_id}.csv missing 'input' column")
        return []
    pred = pred.rename(columns={"input": "smiles"})
    endpoints = [c for c in pred.columns if c not in ("key", "smiles")]

    records = []
    for ep in endpoints:
        ep_pred = pred[["smiles", ep]].rename(columns={ep: "score"})
        train_keys = coadd_training_inchikeys(coadd_root, ep)
        ep_org = _endpoint_organism(ep)
        for assay_code, lab in assays.items():
            base = {
                "endpoint": ep, "endpoint_organism": ep_org,
                "assay_code": assay_code,
                "assay_pathogen": name_by_code.get(assay_code, assay_code),
                "same_organism": ep_org == assay_code,
            }
            records.extend(evaluate(ep_pred, lab, train_keys, base))
    return records


def eos3dys_endpoints(pred_dir, eos_id=COADD_MODEL_ID):
    """The eos3dys output columns present in the step-04 EU OpenScreen prediction file."""
    path = os.path.join(pred_dir, "euopenscreen", f"{eos_id}.csv")
    if not os.path.exists(path):
        return []
    header = pd.read_csv(path, nrows=0).columns.tolist()
    return [c for c in header if c not in ("key", "smiles", "input")]


def matched_endpoints(pred_dir, eos_id=COADD_MODEL_ID):
    """{organism code: {metric: endpoint}} for organisms eos3dys AND EU OpenScreen both cover.

    The organisms this sub-analysis can evaluate own-assay: they need an EU OpenScreen primary
    assay (``SHARED_ORGANISMS``) and an eos3dys endpoint of their own. *E. faecium* has an
    EU OpenScreen assay but no eos3dys endpoint, so it drops out (logged by the callers).

    Several organisms have MORE than one strain endpoint per metric (``ecoli`` has ATCC25922, lpxC
    and tolC; ``paeruginosa`` has ATCC27853, PAO1 and PAO397). The strain is therefore pinned to
    :data:`default.COADD_REF_STRAINS` — the same reference strain the CoAdd evaluation step uses —
    so the two steps describe the same organism. A non-reference strain is used only if the
    reference one is missing, and that substitution is logged; mutant strains (tolC is
    efflux-deficient, PAO397 is a small subset) are not interchangeable with the wild-type
    reference, so picking one silently would quietly change what the figure means.
    """
    matched = {}
    for ep in eos3dys_endpoints(pred_dir, eos_id):
        parsed = _endpoint_to_strain(ep)
        if parsed is None:
            continue                                   # cytotoxicity / hemolytic: no organism
        strain_key, _ = parsed
        code = _endpoint_organism(ep)
        if code not in SHARED_ORGANISMS:
            continue                                   # no EU OpenScreen assay for this organism
        metric = ep[len(strain_key) + 1:]
        ref = COADD_REF_STRAINS.get(code)
        is_ref = strain_key == f"{code}_{ref}" if ref else False
        current = matched.setdefault(code, {}).get(metric)
        if current is None:
            matched[code][metric] = ep
        elif is_ref:
            print(f"  [note] {code} {metric}: using reference strain {ep} (not {current})")
            matched[code][metric] = ep
    # Report any organism/metric left on a non-reference strain.
    for code, by_metric in matched.items():
        ref = COADD_REF_STRAINS.get(code)
        for metric, ep in by_metric.items():
            if ref and not ep.startswith(f"{code}_{ref}_"):
                print(f"  [note] {code} {metric}: reference strain {ref} has no endpoint — "
                      f"falling back to {ep}")
    return matched


def build_eos3dys_overlap_report(pred_dir, euos_root, coadd_root, config, eos_id=COADD_MODEL_ID):
    """Overlap between each endpoint's CoAdd training set and its own EU OpenScreen assay.

    The eos3dys twin of :func:`eval_euopenscreen.build_leakage_report`: per (organism, endpoint),
    how many of the EU OpenScreen compounds — and specifically how many of its ACTIVES — are
    compounds eos3dys already saw in CoAdd training, i.e. how much of the evaluation is genuinely
    novel. This is what the ``dedup`` metrics then remove. Matched on InChIKey (the two sources
    standardise SMILES differently, so a SMILES join would undercount).

    Returns rows with the same schema as the ChEMBL leakage report (pathogen, code, n_train,
    n_eval_conclusive, n_active, n_inactive, n_overlap, n_overlap_active, n_overlap_inactive) plus
    ``endpoint`` and ``metric``, so the existing overlap panel can be reused per endpoint.
    """
    primary = pd.read_csv(os.path.join(euos_root, "primary_assays_manual.csv"))
    code_to_assay = dict(zip(primary["pathogen_code"], primary["assay_eos_id"]))
    name_by_code = dict(zip(config["code"], config["pathogen"]))

    rows = []
    for code, by_metric in sorted(matched_endpoints(pred_dir, eos_id).items()):
        assay_id = code_to_assay.get(code)
        lab = load_euos_primary(euos_root, code, assay_id) if assay_id else None
        if lab is None:
            continue
        eval_keys = set(lab["inchikey"].dropna().astype(str))
        for metric, ep in sorted(by_metric.items()):
            train_keys = coadd_training_inchikeys(coadd_root, ep)
            overlap = eval_keys & (train_keys or set())
            ov = lab[lab["inchikey"].astype(str).isin(overlap)]
            rows.append({
                "pathogen": name_by_code.get(code, code), "code": code,
                "endpoint": ep, "metric": metric,
                "n_train": len(train_keys) if train_keys else 0,
                "n_eval_conclusive": int(len(lab)),
                "n_active": int((lab["bin"] == 1).sum()),
                "n_inactive": int((lab["bin"] == 0).sum()),
                "n_overlap": int(len(overlap)),
                "n_overlap_active": int((ov["bin"] == 1).sum()),
                "n_overlap_inactive": int((ov["bin"] == 0).sum()),
            })
    return rows


def run_eos3dys_exclusivity(pred_dir, euos_root, coadd_root, config, eos_id=COADD_MODEL_ID):
    """eos3dys AUROC on EXCLUSIVE vs SHARED EU OpenScreen hits, per organism and endpoint.

    The eos3dys twin of analysis 3 (``hit_exclusivity``): each organism's exclusive actives (a hit
    in only 1 of the EU OpenScreen primary assays, i.e. organism-specific) and its shared actives
    (a hit in >= 2, i.e. pan-active) are each scored against the SAME primary inactives, by that
    organism's own eos3dys endpoint. A drop towards 0.5 on the exclusive bars means the endpoint is
    capturing generic rather than organism-specific activity.

    Uses the upstream precomputed subsets in ``06_subset_data/exclusivity/`` — the same files the
    ChEMBL analysis uses, so the two figures are directly comparable. Those subsets are defined
    over all 7 EU OpenScreen assays (including *E. faecium*, which eos3dys cannot score); only the
    6 organisms with an eos3dys endpoint are evaluated. Metrics are ``raw`` + ``dedup`` (dedup
    removing that endpoint's CoAdd training compounds).
    """
    primary = pd.read_csv(os.path.join(euos_root, "primary_assays_manual.csv"))
    code_to_assay = dict(zip(primary["pathogen_code"], primary["assay_eos_id"]))
    name_by_code = dict(zip(config["code"], config["pathogen"]))

    path = os.path.join(pred_dir, "euopenscreen", f"{eos_id}.csv")
    if not os.path.exists(path):
        print(f"  [skip] euopenscreen/{eos_id}.csv not present")
        return []
    pred = pd.read_csv(path).rename(columns={"input": "smiles"})

    records = []
    for code, by_metric in sorted(matched_endpoints(pred_dir, eos_id).items()):
        assay_id = code_to_assay.get(code)
        lab = load_euos_primary(euos_root, code, assay_id) if assay_id else None
        if lab is None:
            continue
        inactives = lab[lab["bin"] == 0].copy()
        for metric, ep in sorted(by_metric.items()):
            ep_pred = pred[["smiles", ep]].rename(columns={ep: "score"})
            train_keys = coadd_training_inchikeys(coadd_root, ep)
            for mode in ("exclusive", "nonexclusive"):
                sub = _load_exclusivity_task(euos_root, code, mode, inactives)
                if sub is None:
                    print(f"  [skip] no {mode} subset for {code}")
                    continue
                base = {"pathogen": name_by_code.get(code, code), "code": code,
                        "endpoint": ep, "metric": metric, "subset": mode}
                records.extend(evaluate(ep_pred, sub, train_keys, base))
    return records


def _matched_labels(euos_root, config, codes):
    """{code: DataFrame[smiles, bin, inchikey]} — the own primary assay of each matched organism."""
    primary = pd.read_csv(os.path.join(euos_root, "primary_assays_manual.csv"))
    code_to_assay = dict(zip(primary["pathogen_code"], primary["assay_eos_id"]))
    labels = {}
    for code in codes:
        assay_id = code_to_assay.get(code)
        lab = load_euos_primary(euos_root, code, assay_id) if assay_id else None
        if lab is not None:
            labels[code] = lab
    return labels


def _hit_counts(labels, index):
    """(n_pathogens, n_assays_tested) over the MATCHED organisms only, aligned to ``index``.

    Exclusivity for the eos3dys analyses is counted over just the organisms eos3dys can score, so
    the ranking candidates and the exclusivity definition span the same set. This differs on purpose
    from the upstream 7-assay subsets used by :func:`run_eos3dys_exclusivity` and the ChEMBL panels:
    a compound hit in *E. faecium* plus one matched organism is "shared" there but "exclusive" here.
    """
    n_active = pd.Series(0, index=index, dtype=int)
    n_tested = pd.Series(0, index=index, dtype=int)
    for lab in labels.values():
        aligned = lab.drop_duplicates("smiles").set_index("smiles")["bin"].reindex(index)
        n_tested += aligned.notna().astype(int)
        n_active += (aligned == 1).astype(int)
    return n_active, n_tested


def _combined_scores(pred, matched):
    """DataFrame[smiles x organism] of each organism's inhib_50 + mic_25 summed.

    One score per organism per compound, so the rank analysis ranks the 6 organisms rather than 12
    endpoints. Summing the two heads (rather than taking the better one) keeps both pieces of
    evidence: ``mic_25`` runs systematically above ``inhib_50`` for the same organism, so a max
    would return the MIC head for most compounds. Both are probabilities from one model, so they are
    summed directly with no rescaling. Organisms missing a head contribute the head they have.
    """
    out = {}
    for code, by_metric in matched.items():
        cols = [ep for ep in by_metric.values() if ep in pred.columns]
        if cols:
            out[code] = pred[cols].sum(axis=1)
    return pd.DataFrame(out, index=pred.index)


def endpoint_training_union(coadd_root, matched):
    """Union of the matched endpoints' CoAdd training InChIKeys (the multi-endpoint dedup source).

    The consensus-max panel mixes all 12 endpoints, so a compound is only leakage-free if NO
    endpoint was trained on it — the eos3dys analogue of
    :func:`eval_euopenscreen.shared_training_inchikeys`. Returns None (logged) if nothing is
    readable, so the caller reports raw only rather than a silently-equal "dedup".
    """
    keys = set()
    for by_metric in matched.values():
        for ep in by_metric.values():
            k = coadd_training_inchikeys(coadd_root, ep)
            if k:
                keys.update(k)
    if not keys:
        print("  [note] no CoAdd training sets readable — reporting raw only")
        return None
    return keys


def run_eos3dys_exclusive_rank(pred_dir, euos_root, coadd_root, config, eos_id=COADD_MODEL_ID):
    """For each EXCLUSIVE hit, where does its own organism rank among the matched organisms?

    The eos3dys twin of :func:`eval_euopenscreen.run_exclusive_hit_model_rank`. Exclusive = a hit in
    exactly 1 of the matched organisms' EU OpenScreen assays (see :func:`_hit_counts`); each
    organism is scored by its ``inhib_50 + mic_25`` combined score (see :func:`_combined_scores`),
    the organisms are ranked best-first, and the position of the hit's own organism is recorded —
    rank 1 means its own organism scores it highest.

    Two rankings, as for the ChEMBL models: ``raw`` combined scores, and ``percentile`` (each
    organism's combined score converted to its within-library percentile first). The organisms' raw
    score distributions differ — library medians span ~0.93 to ~1.14 — so a raw ranking partly
    reflects which organism's endpoints output higher values.

    Returns (rank_rows, compound_rows); rank_rows carry per-organism ``n_<code>`` breakdown columns
    so the bars can be stacked by pathogen, matching the ChEMBL panel.
    """
    matched = matched_endpoints(pred_dir, eos_id)
    labels = _matched_labels(euos_root, config, matched)
    name_by_code = dict(zip(config["code"], config["pathogen"]))
    path = os.path.join(pred_dir, "euopenscreen", f"{eos_id}.csv")
    if not os.path.exists(path) or not labels:
        print(f"  [skip] {eos_id} predictions or matched assays unavailable")
        return [], []
    pred = pd.read_csv(path).rename(columns={"input": "smiles"}).drop_duplicates("smiles")
    pred = pred.set_index("smiles")
    mat = _combined_scores(pred, matched)
    mat = mat[mat.notna().all(axis=1)]

    n_active, _ = _hit_counts(labels, mat.index)
    exclusive = n_active[n_active == 1].index
    truth = {}
    for code, lab in labels.items():
        act = set(lab.loc[lab["bin"] == 1, "smiles"])
        for s in exclusive:
            if s in act:
                truth[s] = code
    key_map = {}
    for lab in labels.values():
        key_map.update(lab.dropna(subset=["inchikey"]).drop_duplicates("smiles")
                          .set_index("smiles")["inchikey"].astype(str).to_dict())
    train_union = endpoint_training_union(coadd_root, matched)

    codes = sorted(mat.columns)
    n_models = len(codes)
    rank_rows, compound_rows = [], []
    for leakage in ("raw", "dedup"):
        keep = sorted(truth)
        if leakage == "dedup":
            if not train_union:
                continue
            keep = [s for s in keep if key_map.get(s) not in train_union]
            print(f"  [dedup] {len(sorted(truth)) - len(keep)} of {len(truth)} exclusive hits are "
                  "CoAdd training compounds — dropped")
        if not keep:
            continue
        for ranking in ("raw", "percentile"):
            m = mat if ranking == "raw" else mat.rank(pct=True)
            sub = m.loc[keep, codes]
            ranks = sub.rank(axis=1, ascending=False, method="min")
            n_ties = int(sub.apply(lambda r: r.duplicated().any(), axis=1).sum())
            true_rank = pd.Series([int(ranks.at[s, truth[s]]) for s in keep], index=keep)
            top_code = sub.idxmax(axis=1)
            counts = true_rank.value_counts()
            by_pathogen = pd.crosstab(true_rank, pd.Series([truth[s] for s in keep], index=keep))
            total = len(true_rank)
            for k in range(1, n_models + 1):
                n = int(counts.get(k, 0))
                row = {
                    "ranking": ranking, "leakage": leakage, "rank": k, "n_molecules": n,
                    "frac_molecules": round(n / total, 5) if total else 0.0,
                    "n_molecules_le": int(sum(int(counts.get(j, 0)) for j in range(1, k + 1))),
                    "n_total": total, "n_chance": round(total / n_models, 2),
                    "n_models_ranked": n_models, "n_tied_scores": n_ties,
                }
                for code in codes:
                    row[f"n_{code}"] = int(by_pathogen.at[k, code]) \
                        if k in by_pathogen.index and code in by_pathogen.columns else 0
                rank_rows.append(row)
            for s in keep:
                compound_rows.append({
                    "smiles": s, "ranking": ranking, "leakage": leakage,
                    "pathogen": name_by_code.get(truth[s], truth[s]), "code": truth[s],
                    "true_pathogen_rank": int(true_rank[s]),
                    "top_ranked_pathogen": name_by_code.get(top_code[s], top_code[s]),
                })
            print(f"  {leakage}/{ranking}: rank-1 = {int(counts.get(1, 0))}/{total} "
                  f"(chance {total / n_models:.1f}); {n_ties} tied")
    return rank_rows, compound_rows


def run_eos3dys_consensus_max(pred_dir, euos_root, coadd_root, config, eos_id=COADD_MODEL_ID):
    """MAXIMUM endpoint probability across the matched endpoints, active vs inactive.

    The eos3dys twin of the ChEMBL ``consensus_max_by_activity`` panel: per compound, the single
    highest of the matched organisms' endpoint probabilities — how confident the most confident
    endpoint is. NO percentile normalisation, and none needed: unlike the 7 ChEMBL models (whose
    maxima ran 0.753-0.974, so the low-scoring ones could never win a max), all 12 eos3dys endpoints
    top out at ~1.0, so the scales are aligned where a max reads them.

    Classes, over the matched organisms only: ``active`` = a hit in >= 1 of their EU OpenScreen
    assays, ``inactive`` = ``bin == 0`` in every one it was tested in. Reported ``raw`` and
    ``dedup`` (dedup dropping compounds in ANY matched endpoint's CoAdd training set, both classes,
    so the comparison stays like-for-like).

    Returns (stat_rows, active_rows).
    """
    matched = matched_endpoints(pred_dir, eos_id)
    labels = _matched_labels(euos_root, config, matched)
    path = os.path.join(pred_dir, "euopenscreen", f"{eos_id}.csv")
    if not os.path.exists(path) or not labels:
        print(f"  [skip] {eos_id} predictions or matched assays unavailable")
        return [], []
    pred = pd.read_csv(path).rename(columns={"input": "smiles"})
    pred = pred.drop_duplicates("smiles").set_index("smiles")
    eps = [ep for by_metric in matched.values() for ep in by_metric.values()
           if ep in pred.columns]
    scores = pred[eps]
    scores = scores[scores.notna().all(axis=1)]
    n_active, n_tested = _hit_counts(labels, scores.index)

    df = pd.DataFrame({"consensus_max": scores.max(axis=1), "n_pathogens": n_active,
                       "n_assays_tested": n_tested})
    unlabelled = int((df["n_assays_tested"] == 0).sum())
    df = df[df["n_assays_tested"] > 0]
    df["hit_class"] = np.where(df["n_pathogens"] >= 1, "active", "inactive")

    key_map = {}
    for lab in labels.values():
        key_map.update(lab.dropna(subset=["inchikey"]).drop_duplicates("smiles")
                          .set_index("smiles")["inchikey"].astype(str).to_dict())
    train_union = endpoint_training_union(coadd_root, matched)

    stat_rows, active_rows = [], []
    for leakage in ("raw", "dedup"):
        d = df
        if leakage == "dedup":
            if not train_union:
                continue
            leaked = pd.Series(d.index.map(lambda s: key_map.get(s)), index=d.index) \
                .isin(train_union)
            print(f"  [dedup] dropped {int(leaked.sum())} of {len(d)} compounds present in a "
                  "CoAdd training set")
            d = d[~leaked]
        for hit_class in ACTIVITY_CLASSES:
            v = d.loc[d["hit_class"] == hit_class, "consensus_max"]
            if v.empty:
                continue
            q1, med, q3 = (float(v.quantile(q)) for q in (0.25, 0.5, 0.75))
            iqr = q3 - q1
            inside = v[(v >= q1 - 1.5 * iqr) & (v <= q3 + 1.5 * iqr)]
            stat_rows.append({
                "hit_class": hit_class, "set": leakage, "n": int(len(v)),
                "median": round(med, 4), "q1": round(q1, 4), "q3": round(q3, 4),
                "whisker_lo": round(float(inside.min()), 4),
                "whisker_hi": round(float(inside.max()), 4),
                "mean": round(float(v.mean()), 4),
                "min": round(float(v.min()), 4), "max": round(float(v.max()), 4),
                "n_models_aggregated": len(eps), "n_unlabelled_excluded": unlabelled,
                "score_scale": "raw",
            })
        act = d[d["hit_class"] == "active"].reset_index().rename(columns={"index": "smiles"})
        active_rows.extend(act.assign(set=leakage)[
            ["smiles", "set", "hit_class", "n_pathogens", "n_assays_tested", "consensus_max"]
        ].round({"consensus_max": 4}).to_dict("records"))
    return stat_rows, active_rows


def run_all(pred_dir, euos_root, coadd_root, config_path, output_dir, eos_id=COADD_MODEL_ID):
    """Run the eos3dys-on-EU-OpenScreen analyses and write their summary CSVs."""
    config = pd.read_csv(config_path)  # columns: pathogen, code, eosid
    print(f"[eos3dys] evaluating {eos_id} endpoints on the EU OpenScreen assays ...")
    df = pd.DataFrame(run_eos3dys_euopenscreen(pred_dir, euos_root, coadd_root, config, eos_id))

    matched = matched_endpoints(pred_dir, eos_id)
    skipped = [c for c in SHARED_ORGANISMS if c not in matched]
    if skipped:
        print(f"  [note] no eos3dys endpoint for {', '.join(skipped)} — own-assay analyses cover "
              f"the {len(matched)} organisms that have both")
    print("[eos3dys] training-set overlap with the EU OpenScreen assays ...")
    overlap_df = pd.DataFrame(
        build_eos3dys_overlap_report(pred_dir, euos_root, coadd_root, config, eos_id))
    print("[eos3dys] exclusive vs shared hit AUROC ...")
    excl_df = pd.DataFrame(
        run_eos3dys_exclusivity(pred_dir, euos_root, coadd_root, config, eos_id))
    print("[eos3dys] own-organism rank for exclusive hits ...")
    rank_rows, rank_compounds = run_eos3dys_exclusive_rank(
        pred_dir, euos_root, coadd_root, config, eos_id)
    print("[eos3dys] maximum endpoint probability, active vs inactive ...")
    max_stats, max_actives = run_eos3dys_consensus_max(
        pred_dir, euos_root, coadd_root, config, eos_id)

    os.makedirs(output_dir, exist_ok=True)
    outputs = {
        "eos3dys_euos_auroc.csv": df,
        "eos3dys_overlap_report.csv": overlap_df,
        "eos3dys_hit_exclusivity.csv": excl_df,
        "eos3dys_exclusive_rank.csv": pd.DataFrame(rank_rows),
        "eos3dys_exclusive_rank_compounds.csv": pd.DataFrame(rank_compounds),
        "eos3dys_consensus_max_boxstats.csv": pd.DataFrame(max_stats),
        "eos3dys_consensus_max_actives.csv": pd.DataFrame(max_actives),
    }
    for fname, out in outputs.items():
        out.to_csv(os.path.join(output_dir, fname), index=False)
        print(f"  wrote {fname}: {len(out)} rows")
    return outputs
