"""Shared external-validation primitives, reused by both the EU OpenScreen and CoAdd steps.

These are the dataset-agnostic building blocks — loading a step-04 prediction CSV, loading a
model's ChEMBL training InChIKeys for leakage removal, and the raw/dedup merge + metric
computation — that every external-validation analysis reuses. Dataset-specific loaders and
analyses live in :mod:`eval_euopenscreen` and :mod:`eval_coadd`.
"""

import glob
import os

import pandas as pd

from default import SCORE_COL
from metrics import compute_metrics

METRIC_KEYS = ["n_eval", "n_active", "prevalence",
               "auroc", "auprc", "bedroc", "ef_1pct", "ef_5pct"]


def load_predictions(pred_dir, library, eosid, score_col=SCORE_COL):
    """Load one model's predictions for a library as DataFrame[smiles, score].

    ``library`` in {"euopenscreen", "coadd"}. The step-04 SMILES column is ``input``;
    it is renamed to ``smiles`` for a uniform join key. The headline score is
    ``consensus_score`` when present; single-dataset pathogen models (e.g. campylobacter,
    hpylori, ngonorrhoeae) have no consensus to aggregate and expose their sole output
    column instead, which is then used as the headline (logged). Returns None (logged) if
    the file is missing, has no ``input`` column, or has multiple outputs but no
    ``consensus_score`` (ambiguous — never silently aggregated).
    """
    path = os.path.join(pred_dir, library, f"{eosid}.csv")
    if not os.path.exists(path):
        print(f"  [skip] {library}/{eosid}.csv not present yet")
        return None
    header = pd.read_csv(path, nrows=0).columns.tolist()
    if "input" not in header:
        print(f"  [skip] {library}/{eosid}.csv missing 'input' column")
        return None
    feature_cols = [c for c in header if c not in ("key", "input")]
    if score_col in header:
        use = score_col
    elif len(feature_cols) == 1:
        use = feature_cols[0]
        print(f"  [note] {library}/{eosid}.csv has no '{score_col}'; using sole output '{use}'")
    else:
        print(f"  [skip] {library}/{eosid}.csv has no '{score_col}' and "
              f"{len(feature_cols)} outputs — ambiguous, not aggregating")
        return None
    df = pd.read_csv(path, usecols=["input", use])
    return df.rename(columns={"input": "smiles", use: "score"})[["smiles", "score"]]


def training_inchikeys(models_root, code):
    """Union of training InChIKeys for a pathogen from the ChEMBL models repo.

    Reads ``output/07_datasets/{code}/*.csv`` (columns inchikey, smiles, bin). Returns a
    set, or None (logged once) when the training repo/path is absent — in which case
    de-duplication becomes a no-op and only ``raw`` metrics are reported.
    """
    dpath = os.path.join(models_root, "output", "07_datasets", code)
    if not os.path.isdir(dpath):
        return None
    keys = set()
    for f in sorted(glob.glob(os.path.join(dpath, "*.csv"))):
        try:
            col = pd.read_csv(f, usecols=["inchikey"])["inchikey"]
        except (ValueError, KeyError):
            continue
        keys.update(col.dropna().astype(str).tolist())
    return keys or None


def merged_variants(pred, labels, train_keys):
    """Return [(set_name, merged_df)]: always ("raw", ...), plus ("dedup", ...) when
    training keys are available. The join is on smiles; dedup drops rows whose InChIKey
    is in the model's ChEMBL training set. Single source of the merge/dedup logic."""
    merged = labels.merge(pred, on="smiles", how="inner")
    variants = [("raw", merged)]
    if train_keys and "inchikey" in merged.columns:
        variants.append(("dedup", merged[~merged["inchikey"].astype(str).isin(train_keys)]))
    return variants


def evaluate(pred, labels, train_keys, base_record):
    """Join predictions to labels on smiles and return raw (+ dedup) metric records.

    ``pred`` has [smiles, score]; ``labels`` has [smiles, bin, inchikey]. Emits one record
    per available ``set`` ("raw", and "dedup" when train_keys is not None). Skips (logged)
    when a merged set has fewer than two classes.
    """
    records = []
    for set_name, sub in merged_variants(pred, labels, train_keys):
        if sub["bin"].nunique() < 2:
            continue
        m = compute_metrics(sub["bin"].values, sub["score"].values)
        rec = dict(base_record)
        rec["set"] = set_name
        rec.update({k: m[k] for k in METRIC_KEYS})
        records.append(rec)
    return records
