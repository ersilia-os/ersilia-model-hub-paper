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

from default import COADD_ENDPOINTS, COADD_REF_STRAINS
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

    os.makedirs(output_dir, exist_ok=True)
    outputs = {
        "06_coadd_auroc.csv": coadd_df,
        "06_coadd_leakage_report.csv": leak_df,
    }
    for fname, df in outputs.items():
        df.to_csv(os.path.join(output_dir, fname), index=False)
        print(f"  wrote {fname}: {len(df)} rows")
    return outputs
