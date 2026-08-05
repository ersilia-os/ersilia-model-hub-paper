"""Step 06 analysis engine — evaluate the public ChEMBL pathogen models on CoAdd.

Each organism that has both a ChEMBL model and a CoAdd reference strain (``COADD_REF_STRAINS``,
8 organisms) is scored on that reference strain, for both endpoints (single-point ``inhib_50``
and MIC ``mic_10``); metrics are reported ``raw`` + ``dedup`` (dedup removes the model's ChEMBL
training compounds). Writes small summary CSVs into ``output/06_coadd_validation/`` that
``plots_coadd`` consumes.

This is the mirror of step 05 (EU OpenScreen) on the CoAdd dataset. The CoAdd model ``eos3dys``
is NOT evaluated here — it is validated the other way round (on EU OpenScreen) in step 05's
``eos3dys_validation/``. CoAdd is richer (many strains/cutoffs per organism); this first version
stays on the reference strain + two headline cutoffs, deferring the full multi-strain matrix.

Reuses the shared IO/merge/metric primitives in :mod:`eval_common`.
"""

import os

import pandas as pd

from default import COADD_ENDPOINTS, COADD_HITSET_ENDPOINT, COADD_REF_STRAINS
from eval_common import evaluate, load_predictions, training_inchikeys


def load_coadd_labels(coadd_root, code, subdir, label_col):
    """CoAdd labels for one organism's reference strain as DataFrame[smiles, bin, inchikey].

    File is ``{subdir}/{code}_{strain}.csv`` with SMILES in ``std_smiles`` (renamed to
    ``smiles`` for the join). Returns None (logged) if the strain file is absent — e.g. efaecium
    and spneumoniae are MIC-only, so their single-point inhibition file does not exist.
    """
    strain = COADD_REF_STRAINS[code]
    path = os.path.join(coadd_root, subdir, f"{code}_{strain}.csv")
    if not os.path.exists(path):
        print(f"  [skip] CoAdd {subdir}/{code}_{strain}.csv not found")
        return None
    df = pd.read_csv(path)
    if label_col not in df.columns:
        print(f"  [skip] CoAdd {code}_{strain}.csv missing '{label_col}'")
        return None
    df = df.rename(columns={"std_smiles": "smiles", label_col: "bin"})
    df = df[df["bin"].isin([0, 1])]
    return df[["smiles", "bin", "inchikey"]]


def run_coadd(pred_dir, coadd_root, config, train_cache):
    """Each ChEMBL model on its organism's CoAdd reference strain, both endpoints.

    Iterates ``COADD_REF_STRAINS`` (the organisms with a designated CoAdd reference strain — the
    CoAdd analogue of the EU OpenScreen shared set). Returns a list of metric records with an
    ``endpoint`` field (inhib_50 / mic_10), raw + dedup. Models with no step-04 CoAdd prediction
    yet are skipped (logged).
    """
    name_by_code = dict(zip(config["code"], config["pathogen"]))
    eosid_by_code = dict(zip(config["code"], config["eosid"]))
    records = []
    for code in COADD_REF_STRAINS:
        eosid = eosid_by_code.get(code)
        if eosid is None:
            continue
        pred = load_predictions(pred_dir, "coadd", eosid)
        if pred is None:
            continue
        train_keys = train_cache.get(code)
        for endpoint, (subdir, label_col) in COADD_ENDPOINTS.items():
            lab = load_coadd_labels(coadd_root, code, subdir, label_col)
            if lab is None:
                continue
            base = {
                "pathogen": name_by_code.get(code, code), "code": code, "eosid": eosid,
                "strain": COADD_REF_STRAINS[code], "endpoint": endpoint,
            }
            records.extend(evaluate(pred, lab, train_keys, base))
    return records


def build_coadd_leakage_report(config, coadd_root, train_cache):
    """Overlap between each model's ChEMBL training InChIKeys and its CoAdd reference-strain data.

    One row per (organism, endpoint) with a strain file — independent of whether the step-04
    prediction exists yet. Mirrors ``eval_euopenscreen.build_leakage_report``.
    """
    name_by_code = dict(zip(config["code"], config["pathogen"]))
    modelled = set(config["code"])
    rows = []
    for code in COADD_REF_STRAINS:
        if code not in modelled:
            continue
        train_keys = train_cache.get(code)
        n_train = len(train_keys) if train_keys else 0
        strain = COADD_REF_STRAINS[code]
        for endpoint, (subdir, label_col) in COADD_ENDPOINTS.items():
            lab = load_coadd_labels(coadd_root, code, subdir, label_col)
            if lab is None:
                continue
            eval_keys = set(lab["inchikey"].dropna().astype(str))
            overlap = eval_keys & (train_keys or set())
            ov = lab[lab["inchikey"].astype(str).isin(overlap)]
            rows.append({
                "pathogen": name_by_code.get(code, code), "code": code,
                "strain": strain, "endpoint": endpoint, "n_train": n_train,
                "n_eval_conclusive": int(len(lab)),
                "n_active": int((lab["bin"] == 1).sum()),
                "n_inactive": int((lab["bin"] == 0).sum()),
                "n_overlap": int(len(overlap)),
                "n_overlap_active": int((ov["bin"] == 1).sum()),
                "n_overlap_inactive": int((ov["bin"] == 0).sum()),
            })
    return rows


def _load_hitset_labels(coadd_root, config):
    """Reference-strain label frames for the hit-set endpoint, one per organism that has a file.

    Shared front-end for the three cross-organism analyses (promiscuity, exclusivity, and the
    overlap panel's organism set). Returns ``(labels, codes, name_by_code)`` where ``labels`` maps
    code -> DataFrame[smiles, bin, inchikey] and ``codes`` preserves ``COADD_REF_STRAINS`` order.
    Organisms with no file for this endpoint (efaecium, spneumoniae) are absent, logged by
    :func:`load_coadd_labels`.
    """
    subdir, label_col = COADD_ENDPOINTS[COADD_HITSET_ENDPOINT]
    name_by_code = dict(zip(config["code"], config["pathogen"]))
    modelled = set(config["code"])
    labels = {}
    for code in COADD_REF_STRAINS:
        if code not in modelled:
            continue
        lab = load_coadd_labels(coadd_root, code, subdir, label_col)
        if lab is not None:
            labels[code] = lab
    codes = [c for c in COADD_REF_STRAINS if c in labels]
    return labels, codes, name_by_code


def _hit_counts(labels, codes):
    """Map each distinct active SMILES to the list of organism codes it is a hit in.

    Matching is on the standardised SMILES string, consistent with the rest of step 06 (the
    predictions are joined to labels on ``smiles`` too). Also returns, per SMILES, how many
    organisms actually tested it, so partial coverage can be reported rather than assumed.
    """
    actives = {c: set(labels[c].loc[labels[c]["bin"] == 1, "smiles"]) for c in codes}
    tested = {c: set(labels[c]["smiles"]) for c in codes}
    hit_codes, n_tested = {}, {}
    for smiles in set().union(*[actives[c] for c in codes]) if codes else set():
        hit_codes[smiles] = [c for c in codes if smiles in actives[c]]
        n_tested[smiles] = sum(smiles in tested[c] for c in codes)
    return hit_codes, n_tested


def run_coadd_hit_promiscuity(coadd_root, config):
    """Hit promiscuity — how many CoAdd actives are hits in 1, 2, ... N reference strains.

    The CoAdd mirror of :func:`eval_euopenscreen.run_hit_promiscuity`: label-only (no model), on
    the ``COADD_HITSET_ENDPOINT`` reference strains. ALL distinct actives are counted (the union
    across organisms), including those not tested everywhere — for those the count is a lower
    bound, so they are reported per bin as ``n_incomplete_coverage`` rather than dropped.

    Returns ``(dist_rows, hit_rows)``: the aggregated distribution the figure reads, and the
    per-compound table.
    """
    labels, codes, name_by_code = _load_hitset_labels(coadd_root, config)
    if not codes:
        print("  [skip] no CoAdd reference strains available for hit promiscuity")
        return [], []
    hit_codes, n_tested = _hit_counts(labels, codes)
    keys = {}
    for code in codes:
        lab = labels[code].dropna(subset=["inchikey"])
        keys.update(lab.set_index("smiles")["inchikey"].astype(str).to_dict())

    n_org = len(codes)
    hit_rows = [{
        "smiles": smiles,
        "inchikey": keys.get(smiles, ""),
        "n_pathogens": len(hits),
        "pathogens": ";".join(name_by_code.get(c, c) for c in hits),
        "n_assays_tested": n_tested[smiles],
    } for smiles, hits in hit_codes.items()]
    hit_rows.sort(key=lambda r: (-r["n_pathogens"], r["smiles"]))

    total = len(hit_rows)
    counts = {k: 0 for k in range(1, n_org + 1)}
    incomplete = {k: 0 for k in range(1, n_org + 1)}
    for r in hit_rows:
        counts[r["n_pathogens"]] += 1
        if r["n_assays_tested"] < n_org:
            incomplete[r["n_pathogens"]] += 1
    dist_rows = [{
        "n_pathogens": k,
        "n_molecules": counts[k],
        "frac_molecules": round(counts[k] / total, 5) if total else 0.0,
        "n_molecules_ge": sum(counts[j] for j in range(k, n_org + 1)),
        "n_incomplete_coverage": incomplete[k],
    } for k in range(1, n_org + 1)]
    n_incomplete = sum(incomplete.values())
    print(f"  {total} distinct actives across {n_org} reference strains ({COADD_HITSET_ENDPOINT}); "
          f"{n_incomplete} not tested in all {n_org} (their n_pathogens is a lower bound)")
    return dist_rows, hit_rows


def run_coadd_exclusivity(pred_dir, coadd_root, config, train_cache):
    """Exclusive vs shared hit AUROC per organism, on the CoAdd reference strains.

    The CoAdd mirror of step 05's analysis 3. EU OpenScreen reads precomputed exclusivity subsets
    from ``06_subset_data/exclusivity/``; CoAdd has no such files, so the split is derived here
    from the labels themselves: an active is ``exclusive`` when it is a hit in exactly one of the
    reference strains and ``nonexclusive`` (shared) when it hits two or more. Negatives are the
    organism's own reference-strain inactives, matching how ``_load_exclusivity_task`` builds the
    step-05 task so the two figures stay comparable.

    Only ``dedup`` records are returned — this analysis is reported leakage-filtered only. An
    organism whose subset collapses to a single class (too few actives after filtering) is skipped
    by :func:`evaluate` and logged here, rather than silently omitted.
    """
    labels, codes, name_by_code = _load_hitset_labels(coadd_root, config)
    if not codes:
        print("  [skip] no CoAdd reference strains available for hit exclusivity")
        return []
    hit_codes, _ = _hit_counts(labels, codes)
    eosid_by_code = dict(zip(config["code"], config["eosid"]))

    records = []
    for code in codes:
        eosid = eosid_by_code.get(code)
        if eosid is None:
            continue
        pred = load_predictions(pred_dir, "coadd", eosid)
        if pred is None:
            continue
        lab = labels[code]
        inactives = lab[lab["bin"] == 0]
        actives = lab[lab["bin"] == 1]
        for mode, keep in (("exclusive", lambda n: n == 1), ("nonexclusive", lambda n: n >= 2)):
            sub_act = actives[actives["smiles"].map(lambda s: keep(len(hit_codes[s])))]
            if sub_act.empty:
                print(f"  [skip] {code} {mode}: no actives in this subset")
                continue
            task = pd.concat([sub_act, inactives], ignore_index=True)
            base = {
                "pathogen": name_by_code.get(code, code), "code": code, "eosid": eosid,
                "strain": COADD_REF_STRAINS[code], "endpoint": COADD_HITSET_ENDPOINT,
                "subset": mode,
            }
            got = [r for r in evaluate(pred, task, train_cache.get(code), base)
                   if r["set"] == "dedup"]
            if not got:
                print(f"  [skip] {code} {mode}: no dedup record "
                      f"(no training set, or <2 classes after leakage filtering)")
            records.extend(got)
    return records


def run_all(pred_dir, coadd_root, models_root, config_path, output_dir):
    """Run the CoAdd analyses and write the summary CSVs into ``output_dir``.

    A thin orchestrator, mirroring ``eval_euopenscreen.run_all``. Returns a dict of the written
    DataFrames.
    """
    config = pd.read_csv(config_path)  # columns: pathogen, code, eosid
    train_cache = {code: training_inchikeys(models_root, code) for code in config["code"]}
    if all(v is None for v in train_cache.values()):
        print("[leakage] no training sets found under "
              f"{os.path.join(models_root, 'output', '07_datasets')} — reporting RAW only")

    print("[CoAdd] evaluating own-strain AUROC (inhib_50 + mic_10) ...")
    coadd_df = pd.DataFrame(run_coadd(pred_dir, coadd_root, config, train_cache))
    leak_df = pd.DataFrame(build_coadd_leakage_report(config, coadd_root, train_cache))

    print(f"[CoAdd] hit promiscuity across reference strains ({COADD_HITSET_ENDPOINT}) ...")
    prom_rows, hit_rows = run_coadd_hit_promiscuity(coadd_root, config)

    print(f"[CoAdd] exclusive vs shared hit AUROC, dedup only ({COADD_HITSET_ENDPOINT}) ...")
    excl_df = pd.DataFrame(run_coadd_exclusivity(pred_dir, coadd_root, config, train_cache))

    os.makedirs(output_dir, exist_ok=True)
    outputs = {
        "06_coadd_auroc.csv": coadd_df,
        "06_coadd_leakage_report.csv": leak_df,
        "06_coadd_hit_promiscuity.csv": pd.DataFrame(prom_rows),
        "06_coadd_promiscuous_hits.csv": pd.DataFrame(hit_rows),
        "06_coadd_hit_exclusivity.csv": excl_df,
    }
    for fname, df in outputs.items():
        df.to_csv(os.path.join(output_dir, fname), index=False)
        print(f"  wrote {fname}: {len(df)} rows")
    return outputs
