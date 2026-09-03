"""Steps 08/12 analysis engine — the antibiotic-resemblance score matrix and its UMAP highlights.

``build_abx_named_matrix`` (matrix-building) is called from step 08; ``endpoint_highlights``/
``load_umap`` (UMAP highlight selection) are called from step 12, which reads the matrix step 08
already cached rather than calling ``build_abx_named_matrix`` itself. Per-endpoint summary stats are
computed by :func:`eval_property_matrix.property_endpoint_stats`, shared with the physchem/cytotox
blocks rather than duplicated here.

Same layout and machinery as the pathogen matrices in :mod:`eval_correlations`: every model was run
on the SAME ~1.35M-compound reference library, staged by ``00_download_data.py`` as
``data/processed/annotation_preds_ref_library/{model_id}_{version}.csv`` and aligned on ``key``, so
the predictions form a clean rectangular matrix with no compound alignment needed.

Which endpoints enter is decided by ONE file: ``config/antibiotic_resemblance.csv``,
``selected == Yes`` rows only.

Two things differ from :func:`eval_correlations.build_named_score_matrix`, and both are why this
lives in its own module rather than extending that one:

  1. **No organism dimension.** Antibiotic resemblance is not pathogen-specific, so the config has
     no ``organism`` column and there is no pathogen code to look up. Columns are named
     ``{group}__{model_id}__{column_name}`` with a constant group code (``abx``), keeping the
     three-part shape :func:`eval_correlations.parse_named_column` expects so the blocks can be
     concatenated later.
  2. **The endpoints are almost all discrete.** 54 of the 55 selected columns are binary flags or
     small integer counts; only ``abx_score`` is continuous. That makes a plain top-N rank cutoff
     — which is what step 11 uses on continuous ``consensus_score`` — meaningless here, and drives
     the highlight rule in :func:`endpoint_highlights`.

Nothing is dropped, filtered or thresholded: every selected column reaches the matrix and the stats
table, including the four that are constant zero over the whole library.
"""

import os

import pandas as pd
import pyarrow.parquet as pq

from default import PROJECTION_TOP_N, RANDOM_SEED
from eval_correlations import build_full_library_matrix

#: Group code occupying the first slot of every column name, where the pathogen matrices carry a
#: pathogen code. Makes each column self-identifying once the blocks are joined.
GROUP_CODE = "abx"


def named_column(model_id, column_name, group=GROUP_CODE):
    """``{group}__{model_id}__{column_name}`` — the three-part name used across the joined blocks."""
    return f"{group}__{model_id}__{column_name}"


def build_abx_named_matrix(pred_dir, selection_path, full_matrix_cache_path=None,
                           group=GROUP_CODE, seed=RANDOM_SEED):
    """The ``selected == Yes`` full-library matrix, columns renamed via :func:`named_column`.

    If ``full_matrix_cache_path`` already holds every needed column, only those columns are read
    from it — checked against the parquet's own schema, so no full-file load is needed to decide.
    Otherwise a full :func:`eval_correlations.build_full_library_matrix` build runs and covers every
    endpoint the selection CSV references (``Yes`` AND ``No``), so re-enabling a ``No`` row later
    needs no re-read of the raw prediction files.
    """
    sel_all = pd.read_csv(selection_path)
    sel = sel_all[sel_all["selected"] == "Yes"].copy()
    old_cols = [f"{r.model_id}:{r.column_name}" for r in sel.itertuples()]
    rename_map = {f"{r.model_id}:{r.column_name}": named_column(r.model_id, r.column_name, group)
                  for r in sel.itertuples()}

    if full_matrix_cache_path and os.path.exists(full_matrix_cache_path):
        schema_cols = set(pq.ParquetFile(full_matrix_cache_path).schema.names)
        if set(old_cols) <= schema_cols:
            print(f"[abx-matrix] reusing cached {os.path.basename(full_matrix_cache_path)} "
                  f"({len(old_cols)} of its columns)")
            matrix = pd.read_parquet(full_matrix_cache_path, columns=old_cols)
            return matrix.rename(columns=rename_map), sel
        print(f"[abx-matrix] cache is missing {len(set(old_cols) - schema_cols)} needed column(s) "
              "— building from raw files")

    full = build_full_library_matrix(pred_dir, selection_path, seed=seed)
    if full_matrix_cache_path:
        full.to_parquet(full_matrix_cache_path)
        print(f"[abx-matrix] wrote cache {os.path.basename(full_matrix_cache_path)} {full.shape}")
    return full[old_cols].rename(columns=rename_map), sel


def endpoint_highlights(matrix, stats, proj, top_n=PROJECTION_TOP_N, group=GROUP_CODE):
    """Compounds to highlight per endpoint, with UMAP coordinates attached.

    **The rule (user-directed): every compound with a value > 0, capped at ``top_n``, highest value
    first.** It is deliberately NOT step 11's "top ``top_n`` by score":

      - Only ``abx_score`` is continuous. For the other 54 columns a plain top-N would pad the set
        with arbitrarily chosen ZERO-valued compounds — e.g. ``glycopeptides`` has 1 non-zero
        compound in 1,355,109, so a top-1000 would be 1 real hit and 999 tie-broken zeros drawn as
        a cloud that looks like a finding.
      - Capping never pads: an endpoint with 3 hits contributes 3 rows, not ``top_n``.
      - Where the cap binds on a binary column (e.g. ``sulfonamides``, 76,473 flagged) the 1000
        kept are an ARBITRARY subset — ties are broken by ``key`` order, which is deterministic and
        reproducible but carries no meaning. The ``n_shown``/``n_nonzero`` pair in ``stats`` is what
        makes that visible; do not read a capped panel as exhaustive.

    Endpoints with ``n_nonzero == 0`` yield no rows at all (nothing to draw). They stay in ``stats``.

    Returns ``(highlights, stats)`` — ``stats`` gains an ``n_shown`` column.
    """
    coords = proj.set_index("key")[["umap_x", "umap_y"]]
    parts, shown = [], {}
    for r in stats.itertuples():
        v = matrix[r.endpoint]
        nz = v[v > 0]
        if not len(nz):
            shown[r.endpoint] = 0
            continue
        top = nz.nlargest(top_n)
        block = pd.DataFrame({"endpoint": r.endpoint, "model_id": r.model_id,
                              "column_name": r.column_name, "key": top.index,
                              "value": top.to_numpy()})
        block = block.join(coords, on="key")
        parts.append(block)
        shown[r.endpoint] = len(block)
    stats = stats.copy()
    stats["n_shown"] = stats["endpoint"].map(shown).astype(int)
    highlights = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(
        columns=["endpoint", "model_id", "column_name", "key", "value", "umap_x", "umap_y"])
    return highlights, stats


def load_umap(projection_file):
    """eos1klk's ``key`` + UMAP coordinate pair only (2 of its 8 coordinate columns)."""
    df = pd.read_csv(projection_file, usecols=["key", "umap_x", "umap_y"])
    for c in ("umap_x", "umap_y"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df
