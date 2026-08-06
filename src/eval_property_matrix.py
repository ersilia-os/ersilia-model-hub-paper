"""Shared builder for the per-family property matrices (steps 10 and 12).

One function assembles a ``key x endpoint`` block over the full reference library from a curated
selection config, so each property family gets its own step without duplicating the collection
logic:

  - step 10 ``config/physchem_models.csv``      -> ``physchem__{model_id}__{column}``  (22 columns)
  - step 12 ``config/cytotoxicity_models.csv``  -> ``cytotox__{model_id}__{column}``   (24 columns)

The ``{prefix}__{model_id}__{column_name}`` naming is the same three-part shape used by the pathogen
matrix (:func:`eval_correlations.build_named_score_matrix`, where the prefix is the pathogen code)
and by the abx block (:mod:`eval_abx_matrix`, where it is a constant group code). Carrying the model
ID keeps provenance readable where two models cover overlapping biology.

**Un-normalized only**, matching the abx block: no scaling and no row normalization. Choosing a
transform is a decision for after the family blocks are joined, and the stats these functions write
are what that decision should be made from.

**Nothing is dropped.** Missing values are counted and reported, never imputed or filtered, and a
constant column appears in the stats with ``n_unique == 1`` rather than being removed.
"""

import glob
import os

import numpy as np
import pandas as pd

KEY_COL = "key"
SMILES_COL = "input"


def resolve_prediction_file(pred_dir, model_id):
    """Highest-version Isaura prediction file for a model, e.g. ``eos4djh_v1.csv``."""
    matches = sorted(glob.glob(os.path.join(pred_dir, f"{model_id}_v*.csv")))
    if not matches:
        raise FileNotFoundError(
            f"No prediction file for {model_id} in {pred_dir} — run 00_download_data.py first."
        )
    return matches[-1]


def selected_endpoints(config_csv):
    """``model_id -> [column_name, ...]`` for the ``selected == "Yes"`` rows, order preserved."""
    config = pd.read_csv(config_csv)
    sel = config[config["selected"] == "Yes"]
    endpoints = {}
    for model_id, group in sel.groupby("model_id", sort=False):
        endpoints[model_id] = group["column_name"].tolist()
    return endpoints, len(config), len(sel)


def build_property_matrix(pred_dir, config_csv, prefix, verbose=True):
    """The ``selected == "Yes"`` block for one property family, as a wide DataFrame.

    Returns ``(matrix, endpoint_names)``. ``matrix`` carries ``key``, ``input`` and one column per
    endpoint named ``{prefix}__{model_id}__{column_name}``; ``endpoint_names`` is those column
    names in config order.

    Prediction files are concatenated column-wise on the assumption that they share a row order,
    which is **verified** per model against the first file's keys; any file whose order differs is
    reindexed on ``key`` rather than concatenated blindly.
    """
    endpoints, n_config, n_selected = selected_endpoints(config_csv)
    if verbose:
        print(f"{os.path.basename(config_csv)}: {n_selected} of {n_config} endpoints selected "
              f"-> '{prefix}' columns")

    rename_map = {(model_id, col): f"{prefix}__{model_id}__{col}"
                  for model_id, cols in endpoints.items() for col in cols}
    endpoint_names = list(rename_map.values())
    duplicates = sorted({c for c in endpoint_names if endpoint_names.count(c) > 1})
    if duplicates:
        raise ValueError(f"Endpoint names collide: {duplicates}")

    frames = []
    reference_keys = None
    reference_model = None
    for model_id, cols in endpoints.items():
        path = resolve_prediction_file(pred_dir, model_id)
        if verbose:
            print(f"  reading {os.path.basename(path)} ({len(cols)} columns)...")
        df = pd.read_csv(path, usecols=[KEY_COL, SMILES_COL] + cols)
        if reference_keys is None:
            reference_keys = df[KEY_COL]
            reference_model = model_id
            frames.append(df[[KEY_COL, SMILES_COL]])
        elif not df[KEY_COL].equals(reference_keys):
            if verbose:
                print(f"    key order differs from {reference_model}; aligning on '{KEY_COL}'.")
            df = df.set_index(KEY_COL).reindex(reference_keys).reset_index()
        frames.append(df[cols].rename(columns={c: rename_map[(model_id, c)] for c in cols}))

    matrix = pd.concat(frames, axis=1)
    if verbose:
        print(f"  matrix: {len(matrix):,} rows x {len(endpoint_names)} endpoints")
    return matrix, endpoint_names


def property_endpoint_stats(matrix, endpoint_names, config_csv=None):
    """Per-endpoint summary over the FULL library — one row per endpoint, nothing excluded.

    Mirrors :func:`eval_abx_matrix.endpoint_stats` but does not require a ``direction`` column,
    which the physchem config deliberately omits (no descriptor has a "better" end). If the config
    does carry one it is joined on, so the cytotox stats keep it.

    A constant column appears with ``n_unique == 1`` rather than being filtered out — that is what
    makes it auditable, since a constant column cannot be z-scored (std = 0) and has no top-N.
    """
    rows = []
    for name in endpoint_names:
        _, model_id, column_name = name.split("__", 2)
        v = matrix[name]
        finite = v[np.isfinite(v)]
        rows.append({
            "endpoint": name, "model_id": model_id, "column_name": column_name,
            "n_unique": int(v.nunique()),
            "min": float(finite.min()) if len(finite) else np.nan,
            "max": float(finite.max()) if len(finite) else np.nan,
            "mean": float(finite.mean()) if len(finite) else np.nan,
            "median": float(finite.median()) if len(finite) else np.nan,
            "std": float(finite.std()) if len(finite) else np.nan,
            "n_nonzero": int((v > 0).sum()),
            "n_nan": int(v.isna().sum()),
        })
    out = pd.DataFrame(rows)
    out["pct_nonzero"] = 100.0 * out["n_nonzero"] / len(matrix)
    if config_csv is not None:
        config = pd.read_csv(config_csv)
        if "direction" in config.columns:
            out = out.merge(config[["model_id", "column_name", "direction"]],
                            on=["model_id", "column_name"], how="left")
    return out


def report_missing(matrix, endpoint_names):
    """Print missing-value coverage per endpoint. Reporting only — nothing is dropped."""
    print("\nMissing values per endpoint:")
    n = len(matrix)
    for name in endpoint_names:
        n_missing = int(matrix[name].isna().sum())
        print(f"  {name}: {n_missing} missing ({100 * n_missing / n:.2f}%)")
