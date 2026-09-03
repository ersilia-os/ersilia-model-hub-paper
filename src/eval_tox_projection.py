"""Step 13 analysis engine — reference-library projection coloured by predicted toxicity.

Same layout and same background as :mod:`eval_projection` (step 11), but the highlighted overlay
is per *toxicity endpoint* rather than per pathogen: for each endpoint selected in
``config/cytotoxicity_models.csv`` this finds the ``PROJECTION_TOP_N`` most toxic compounds in the
step-08 cytotox matrix — a rank cutoff, never a score threshold — and attaches their
:data:`default.TOX_PROJECTION_METHOD` coordinates.

"Most toxic" is the HIGHEST value for every endpoint (:data:`default.TOX_RANK_DESCENDING`); see
that constant for why this also holds for ``ld50_zhu``, the one regression endpoint, whose
reciprocal-log units invert the dose scale.

Memory: the step-08 cytotox matrix is ~690 MB, most of it the ``input`` SMILES column, which is
never read here. The 24 score columns are streamed in chunks and each endpoint's running top-N is
reduced after every chunk, so no more than one chunk plus 24 x N rows is ever held.
"""

import os

import pandas as pd

from default import (CORR_CHUNK_SIZE, PROJECTION_TOP_N, TOX_PREFIX,
                     TOX_PROJECTION_METHOD, TOX_RANK_DESCENDING)
from eval_projection import load_projection


def selected_endpoints(config_csv, prefix=TOX_PREFIX):
    """The ``selected == "Yes"`` rows of the endpoint config, with their step-08 column names.

    Returns a DataFrame with ``model_id``, ``column_name`` and ``endpoint`` (the prefixed
    ``{prefix}__{model_id}__{column_name}`` name as written by step 08), in config file order so
    the figure's panels stay grouped by model.
    """
    config = pd.read_csv(config_csv)
    sel = config[config["selected"] == "Yes"].copy()
    sel["endpoint"] = [f"{prefix}__{r.model_id}__{r.column_name}" for r in sel.itertuples()]
    return sel[["model_id", "slug", "column_name", "endpoint"]].reset_index(drop=True)


def endpoint_top_n(properties_csv, endpoints, n=PROJECTION_TOP_N, chunk_size=CORR_CHUNK_SIZE,
                   descending=TOX_RANK_DESCENDING):
    """Each endpoint's ``n`` most toxic compounds, as ``{endpoint: DataFrame[key, score]}``.

    Streams ``key`` + the endpoint columns only (never ``input``) and re-reduces every endpoint's
    running top-N after each chunk, so memory stays flat regardless of library size. ``nlargest``
    is used because every endpoint ranks most-toxic-highest; ``descending=False`` flips to
    ``nsmallest`` for a hypothetical future endpoint that does not.
    """
    take = (lambda df, col: df.nlargest(n, col)) if descending else \
           (lambda df, col: df.nsmallest(n, col))
    tops = {e: None for e in endpoints}
    n_rows = 0
    for chunk in pd.read_csv(properties_csv, usecols=["key"] + list(endpoints),
                             chunksize=chunk_size):
        n_rows += len(chunk)
        for e in endpoints:
            part = chunk[["key", e]]
            tops[e] = part if tops[e] is None else pd.concat([tops[e], part], ignore_index=True)
            tops[e] = take(tops[e], e)
    print(f"[tox-projection] ranked {n_rows} compounds across {len(endpoints)} endpoints")
    return {e: df.rename(columns={e: "score"}).reset_index(drop=True) for e, df in tops.items()}


def attach_coordinates(tops, proj, method=TOX_PROJECTION_METHOD):
    """Join every endpoint's top-N keys to their ``{method}_x/_y`` coordinates, one long table."""
    coord_cols = [f"{method}_{ax}" for ax in ("x", "y")]
    frames = []
    for endpoint, df in tops.items():
        merged = df.merge(proj[["key"] + coord_cols], on="key", how="left")
        merged.insert(0, "endpoint", endpoint)
        frames.append(merged)
    return pd.concat(frames, ignore_index=True)


def run_all(projection_file, properties_csv, config_csv, output_dir, background_path,
            method=TOX_PROJECTION_METHOD, top_n=PROJECTION_TOP_N,
            chunk_size=CORR_CHUNK_SIZE):
    """Step 13's toxicity projection: each endpoint's top-``top_n`` on the shared background grid.

    Writes ``13_top{top_n}_per_endpoint.csv`` (one row per (endpoint, compound)).

    **The background grid is REUSED, not recomputed.** ``background_path`` points at step 11's
    ``11_{method}_background.csv`` — the same full-library density over the same ``eos1klk`` layout,
    so the toxicity panels sit on a background identical to the pathogen ones (step 11) and the abx
    ones (step 12, which already reused it). Recomputing it here produced a byte-identical file, so
    the only thing the duplicate bought was a second chance to drift.
    """
    endpoints_df = selected_endpoints(config_csv)
    endpoints = endpoints_df["endpoint"].tolist()
    print(f"[tox-projection] {len(endpoints)} selected endpoints from "
          f"{os.path.basename(config_csv)}")

    proj = load_projection(projection_file)
    print(f"[tox-projection] loaded {os.path.basename(projection_file)} for {len(proj)} compounds")

    if not os.path.exists(background_path):
        raise FileNotFoundError(
            f"Missing {background_path}. Run `python 11_reference_library_projection.py` first — "
            "step 13 reuses its background grid rather than recomputing it.")
    print(f"  [{method}] reusing background grid {os.path.basename(background_path)}")

    tops = endpoint_top_n(properties_csv, endpoints, n=top_n, chunk_size=chunk_size)
    table = attach_coordinates(tops, proj, method=method)

    # Report where each endpoint's rank cutoff landed — no threshold is applied, but a top-N whose
    # scores are all near-saturated (e.g. dili, whose library median is already ~0.78) means the
    # cutoff is separating far less than it does for a sharply-peaked endpoint.
    for r in endpoints_df.itertuples():
        g = table[table["endpoint"] == r.endpoint]
        missing = int(g[f"{method}_x"].isna().sum())
        note = f", {missing} without coordinates" if missing else ""
        print(f"  {r.model_id}/{r.column_name}: top {len(g)}/{top_n} "
              f"(score {g['score'].min():.3f}-{g['score'].max():.3f}){note}")

    out_path = os.path.join(output_dir, f"13_top{top_n}_per_endpoint.csv")
    table.to_csv(out_path, index=False)
    print(f"  -> {os.path.basename(out_path)} ({len(table)} rows)")
    return endpoints_df
