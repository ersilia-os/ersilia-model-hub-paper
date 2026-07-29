"""Step 07 analysis engine — inter-model prediction correlations on the reference library.

The annotation models were all run on the SAME ~1.35M-compound reference library and staged by
``00_download_data.py`` as ``data/processed/annotation_preds_ref_library/{model_id}_{version}.csv``,
each with join columns ``key`` (compound hash) + ``input`` (SMILES) and model-specific output
columns. Because every file is aligned on ``key`` this is a clean rectangular matrix — no compound
alignment needed.

This engine, at the **output-column** level (each output column is a node):
  1. samples ``CORR_SAMPLE_N`` compounds (``RANDOM_SEED``) and builds a ``key × node`` score matrix,
     caching it as ``07_score_matrix.parquet``;
  2. tags every node's value type (probability / continuous / categorical);
  3. auto-assigns grouping metadata (Target Organism, cytotoxicity) from Airtable — written for
     USER REVIEW before figures;
  4. computes a library-wide **Spearman** correlation matrix across all nodes;
  5. computes **top-N Jaccard overlap** (N in ``TOPN_CUTOFFS``) for probability nodes only;
  6. summarises within-group vs cross-group correlation for the review groups.

Steps 1–3 are the ``build`` stage (stop here to review ``07_group_assignments.csv``); steps 4–6 are
the ``analyze`` stage. Depends only on pandas / numpy / scipy — no project eval primitives, because
the data layout (``key``/``input``) differs from the ``smiles``/``bin`` contract in ``eval_common``.
"""

import glob
import os
import re

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import squareform

from default import (
    CATEGORICAL_MAX_UNIQUE,
    CORR_CHUNK_SIZE,
    CORR_SAMPLE_N,
    CYTOTOX_REGEX,
    ORGANISM_EXCLUDE,
    ORGANISM_MIN_MODELS,
    RANDOM_SEED,
    TOPN_CUTOFFS,
    TOX_COLUMN_REGEX,
    TOX_TITLE_REGEX,
)

JOIN_COLS = ("key", "input")
MATRIX_FILE = "07_score_matrix.parquet"


# --------------------------------------------------------------------------- #
# File selection                                                               #
# --------------------------------------------------------------------------- #
def latest_version_files(pred_dir):
    """Map each model_id to its highest-version prediction file.

    Files are ``{model_id}_v{n}.csv``; some models ship several versions. Only one file per model
    becomes a set of nodes (two versions of one model would correlate ~1.0 and inflate the matrix).
    Returns ``{model_id: (version, path)}`` keeping the numerically largest version.
    """
    chosen = {}
    for path in sorted(glob.glob(os.path.join(pred_dir, "*.csv"))):
        stem = os.path.basename(path)[:-4]
        model_id, _, ver = stem.rpartition("_")
        if not model_id:
            continue
        vnum = int(re.sub(r"\D", "", ver) or 0)
        if model_id not in chosen or vnum > chosen[model_id][0]:
            chosen[model_id] = (vnum, path)
    return {m: (f"v{v}", p) for m, (v, p) in chosen.items()}


# --------------------------------------------------------------------------- #
# Score matrix                                                                 #
# --------------------------------------------------------------------------- #
def sample_keys(reference_file, n, seed):
    """Draw ``n`` compound keys (in a fixed order) from one prediction file's ``key`` column.

    Every file is aligned on ``key``, so a sample drawn from any one file indexes all of them.
    ``n=None`` returns the full key set. Uses ``numpy.random.RandomState(seed)`` for reproducibility.
    """
    keys = pd.read_csv(reference_file, usecols=["key"])["key"]
    if n is None or n >= len(keys):
        return keys.to_numpy()
    rng = np.random.RandomState(seed)
    idx = np.sort(rng.choice(len(keys), size=n, replace=False))
    return keys.to_numpy()[idx]


def _read_model_columns(path, key_order, chunk_size):
    """Read one prediction file's output columns, filtered to ``key_order``, in chunks.

    Returns a DataFrame indexed by ``key`` (reindexed to ``key_order``, so all models align) with
    one float column per model output. Chunked so the 15.9 GB / 619-column file never loads whole.
    """
    key_set = set(key_order)
    parts = []
    for chunk in pd.read_csv(path, chunksize=chunk_size):
        keep = chunk[chunk["key"].isin(key_set)]
        if len(keep):
            parts.append(keep)
    if not parts:
        return None
    df = pd.concat(parts, ignore_index=True).drop_duplicates(subset=["key"])
    out_cols = [c for c in df.columns if c not in JOIN_COLS]
    df = df.set_index("key")[out_cols].apply(pd.to_numeric, errors="coerce")
    return df.reindex(key_order).astype(np.float32)


def build_score_matrix(files, key_order, chunk_size=CORR_CHUNK_SIZE, verbose=True):
    """Assemble the ``key × node`` score matrix from every model file.

    ``files`` is ``{model_id: (version, path)}``. Node columns are named ``"{model_id}:{output}"``.
    Non-numeric / failed outputs become NaN. Returns a float32 DataFrame indexed by ``key``.
    """
    blocks = []
    for i, (model_id, (ver, path)) in enumerate(sorted(files.items()), 1):
        if verbose:
            size_gb = os.path.getsize(path) / 1e9
            print(f"  [{i}/{len(files)}] {model_id} {ver} ({size_gb:.1f} GB) ...", flush=True)
        df = _read_model_columns(path, key_order, chunk_size)
        if df is None or df.shape[1] == 0:
            print(f"    [skip] {model_id}: no usable output columns")
            continue
        df.columns = [f"{model_id}:{c}" for c in df.columns]
        blocks.append(df)
    matrix = pd.concat(blocks, axis=1)
    matrix.index.name = "key"
    return matrix


# --------------------------------------------------------------------------- #
# Column classification                                                        #
# --------------------------------------------------------------------------- #
def classify_columns(matrix, model_of):
    """Tag each node's value type from its sampled values.

    Returns a DataFrame ``[node, model_id, output_col, n_nonnull, n_unique, vmin, vmax,
    value_type]`` with ``value_type`` in {probability, continuous, categorical}. Only probability
    columns are eligible for top-N overlap.
    """
    rows = []
    for node in matrix.columns:
        s = matrix[node].dropna()
        n_unique = int(s.nunique())
        vmin = float(s.min()) if len(s) else np.nan
        vmax = float(s.max()) if len(s) else np.nan
        integer_like = len(s) and np.allclose(s.values, np.round(s.values), equal_nan=False)
        if n_unique <= CATEGORICAL_MAX_UNIQUE and integer_like:
            vtype = "categorical"
        elif len(s) and vmin >= 0.0 and vmax <= 1.0:
            vtype = "probability"
        else:
            vtype = "continuous"
        model_id = node.split(":", 1)[0]
        rows.append({
            "node": node, "model_id": model_id, "output_col": node.split(":", 1)[1],
            "n_nonnull": int(len(s)), "n_unique": n_unique, "vmin": vmin, "vmax": vmax,
            "value_type": vtype,
        })
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# Grouping metadata (auto — reviewed by the user before figures)              #
# --------------------------------------------------------------------------- #
def _split_multi(value):
    """Split a comma-separated Airtable multi-value cell into a clean list (drops 'Any'/blanks)."""
    if pd.isna(value):
        return []
    return [p.strip() for p in str(value).split(",") if p.strip() and p.strip().lower() != "any"]


def assign_groups(meta, model_ids):
    """Auto group assignments per model for review: organism membership + cytotoxicity flag.

    Organism = split ``Target Organism`` (multi-value, 'Any' dropped). Cytotoxicity = case-insensitive
    ``CYTOTOX_REGEX`` over Tag / Title / Interpretation / Description, with the matched text recorded
    as evidence so the user can audit each flag. Returns one row per model.
    """
    pat = re.compile(CYTOTOX_REGEX, re.IGNORECASE)
    by_id = meta.set_index("Identifier")
    rows = []
    for model_id in sorted(model_ids):
        m = by_id.loc[model_id] if model_id in by_id.index else None
        title = "" if m is None else str(m.get("Title", ""))
        organisms = [] if m is None else _split_multi(m.get("Target Organism"))
        blob_fields = [] if m is None else [
            str(m.get(c, "")) for c in ("Tag", "Title", "Interpretation", "Description")
        ]
        hits = sorted({h.group(0) for f in blob_fields for h in pat.finditer(f)})
        rows.append({
            "model_id": model_id,
            "title": title,
            "target_organism_raw": "" if m is None else str(m.get("Target Organism", "")),
            "organisms": ";".join(organisms),
            "biomedical_area": "" if m is None else str(m.get("Biomedical Area", "")),
            "subtask": "" if m is None else str(m.get("Subtask", "")),
            "is_cytotox": bool(hits),
            "cytotox_evidence": ";".join(hits),
        })
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# Correlation & overlap                                                        #
# --------------------------------------------------------------------------- #
def spearman_matrix(matrix):
    """Library-wide Spearman correlation across all nodes (rank then Pearson, BLAS fast path).

    Ranks every column then takes ``numpy.corrcoef`` on the ranks (Spearman = Pearson-of-ranks),
    which is a single BLAS call — orders of magnitude faster than pandas pairwise for ~1500 nodes.
    Any NaN prediction cells are imputed to the column's mean rank first (count logged); a
    constant / all-NaN column yields NaN correlations (rendered as blank cells downstream).
    Returns a square DataFrame indexed by node.
    """
    ranked = matrix.rank()
    n_nan = int(ranked.isna().to_numpy().sum())
    if n_nan:
        print(f"  [spearman] imputed {n_nan} NaN rank cells to column mean before correlation")
        ranked = ranked.fillna(ranked.mean())
    corr = np.corrcoef(ranked.to_numpy(dtype=np.float64), rowvar=False)
    return pd.DataFrame(corr, index=matrix.columns, columns=matrix.columns)


def topn_overlap_matrix(matrix, nodes, n):
    """Pairwise Jaccard overlap of the top-``n`` highest-scoring compounds among ``nodes``.

    Builds a boolean top-``n`` membership matrix (compounds × nodes); intersection counts come from
    ``Mᵀ·M`` and union from ``|A|+|B|−|A∩B|``. Returns a square DataFrame of Jaccard values.
    """
    sub = matrix[nodes].to_numpy(dtype=np.float64)
    n = min(n, sub.shape[0])
    # argsort descending; NaN sorts to the end so it never enters a top-n set.
    order = np.argsort(-sub, axis=0)[:n, :]
    mask = np.zeros(sub.shape, dtype=bool)
    for j in range(sub.shape[1]):
        col = sub[:, j]
        valid = order[:, j][np.isfinite(col[order[:, j]])]
        mask[valid, j] = True
    inter = mask.T.astype(np.int64) @ mask.astype(np.int64)
    sizes = mask.sum(axis=0)
    union = sizes[:, None] + sizes[None, :] - inter
    with np.errstate(divide="ignore", invalid="ignore"):
        jac = np.where(union > 0, inter / union, np.nan)
    return pd.DataFrame(jac, index=nodes, columns=nodes)


def cluster_order(corr):
    """Hierarchical-clustering leaf order for a correlation matrix (average linkage on 1−corr).

    Returns the node labels in clustered order so the global heatmap shows block structure. Falls
    back to the original order if the matrix is too small or degenerate.
    """
    if corr.shape[0] < 3:
        return list(corr.columns)
    d = 1.0 - corr.to_numpy(dtype=np.float64)
    d = np.nan_to_num((d + d.T) / 2.0, nan=1.0)
    np.fill_diagonal(d, 0.0)
    try:
        order = leaves_list(linkage(squareform(d, checks=False), method="average"))
        return [corr.columns[i] for i in order]
    except Exception as exc:  # pragma: no cover - defensive
        print(f"  [cluster] falling back to input order: {exc}")
        return list(corr.columns)


def group_correlation_summary(corr, col_index, groups):
    """Within-group vs cross-group mean/median |Spearman| for each review group.

    ``groups`` maps a group label to the set of member nodes. For each group, compares the mean and
    median absolute correlation among its own nodes against the mean absolute correlation of those
    nodes with everything outside the group. Returns a tidy DataFrame.
    """
    abs_corr = corr.abs()
    all_nodes = list(corr.columns)
    rows = []
    for label, members in groups.items():
        members = [m for m in members if m in abs_corr.index]
        if len(members) < 2:
            continue
        within = abs_corr.loc[members, members].to_numpy()
        iu = np.triu_indices(len(members), k=1)
        within_vals = within[iu]
        others = [n for n in all_nodes if n not in set(members)]
        cross_vals = abs_corr.loc[members, others].to_numpy().ravel() if others else np.array([])
        rows.append({
            "group": label, "n_nodes": len(members),
            "within_mean_abs_rho": float(np.nanmean(within_vals)),
            "within_median_abs_rho": float(np.nanmedian(within_vals)),
            "cross_mean_abs_rho": float(np.nanmean(cross_vals)) if len(cross_vals) else np.nan,
        })
    return pd.DataFrame(rows)


def cytotox_nodes(col_index, group_df):
    """Nodes in the precise cytotoxicity focus group (see TOX_COLUMN_REGEX / TOX_TITLE_REGEX).

    A node joins if its output-column name matches ``TOX_COLUMN_REGEX`` (cell/organ-tox endpoints
    across any model) or its model title matches the narrower ``TOX_TITLE_REGEX`` (rescues dedicated
    hERG/cardiotox/DILI models with generically-named outputs). Returns the node list.
    """
    col_pat = re.compile(TOX_COLUMN_REGEX, re.IGNORECASE)
    title_pat = re.compile(TOX_TITLE_REGEX, re.IGNORECASE)
    titles = group_df["title"].fillna("") if "title" in group_df.columns \
        else pd.Series("", index=group_df.index)
    title_by_model = dict(zip(group_df["model_id"], titles))
    nodes = []
    for r in col_index.itertuples():
        col_hit = bool(col_pat.search(str(r.output_col)))
        title_hit = bool(title_pat.search(str(title_by_model.get(r.model_id, ""))))
        if col_hit or title_hit:
            nodes.append(r.node)
    return nodes


def build_groups(col_index, group_df):
    """Assemble ``{group_label: [nodes]}`` for the cytotoxicity and per-organism focus groups.

    Cytotoxicity membership is the refined column/title rule (:func:`cytotox_nodes`). Per-organism
    groups come from each model's ``organisms`` list; ``ORGANISM_EXCLUDE`` organisms are dropped, as
    are organisms with fewer than ``ORGANISM_MIN_MODELS`` distinct models (a correlation needs two).
    """
    org_map = {r.model_id: r.organisms.split(";") if r.organisms else []
               for r in group_df.itertuples()}
    groups = {}
    cyto = cytotox_nodes(col_index, group_df)
    if cyto:
        groups["cytotoxicity"] = cyto
    node_model = dict(zip(col_index["node"], col_index["model_id"]))
    org_to_nodes = {}
    for node, model_id in node_model.items():
        for org in org_map.get(model_id, []):
            if org in ORGANISM_EXCLUDE:
                continue
            org_to_nodes.setdefault(org, []).append(node)
    for org, nodes in org_to_nodes.items():
        n_models = len({node_model[n] for n in nodes})
        if n_models >= ORGANISM_MIN_MODELS:
            groups[f"organism:{org}"] = nodes
    return groups


# --------------------------------------------------------------------------- #
# Stage orchestrators                                                          #
# --------------------------------------------------------------------------- #
def run_build(pred_dir, meta_path, output_dir, sample_n=CORR_SAMPLE_N, seed=RANDOM_SEED):
    """Build stage: score matrix + column index + auto group assignments (then STOP for review).

    Idempotent: if ``07_score_matrix.parquet`` already exists it is reused. Writes the matrix,
    ``07_column_index.csv`` and ``07_group_assignments.csv`` into ``output_dir``.
    """
    os.makedirs(output_dir, exist_ok=True)
    files = latest_version_files(pred_dir)
    print(f"[build] {len(files)} models (latest version each)")

    matrix_path = os.path.join(output_dir, MATRIX_FILE)
    if os.path.exists(matrix_path):
        print(f"[build] reusing cached matrix {matrix_path}")
        matrix = pd.read_parquet(matrix_path)
    else:
        ref_file = next(iter(files.values()))[1]
        key_order = sample_keys(ref_file, sample_n, seed)
        print(f"[build] sampled {len(key_order)} compounds (seed={seed})")
        matrix = build_score_matrix(files, key_order)
        matrix.to_parquet(matrix_path)
        print(f"[build] wrote {matrix_path}  shape={matrix.shape}")

    meta = pd.read_csv(meta_path)
    col_index = classify_columns(matrix, model_of=None)
    group_df = assign_groups(meta, {n.split(':', 1)[0] for n in matrix.columns})
    # Fold value-type + group metadata into the column index for downstream convenience.
    cyto_models = set(group_df.loc[group_df["is_cytotox"], "model_id"])
    org_map = dict(zip(group_df["model_id"], group_df["organisms"]))
    col_index["is_cytotox"] = col_index["model_id"].isin(cyto_models)
    col_index["organisms"] = col_index["model_id"].map(org_map).fillna("")

    col_index.to_csv(os.path.join(output_dir, "07_column_index.csv"), index=False)
    group_df.to_csv(os.path.join(output_dir, "07_group_assignments.csv"), index=False)
    print(f"[build] {len(col_index)} nodes | "
          f"{(col_index['value_type'] == 'probability').sum()} probability, "
          f"{(col_index['value_type'] == 'continuous').sum()} continuous, "
          f"{(col_index['value_type'] == 'categorical').sum()} categorical")
    print(f"[build] {group_df['is_cytotox'].sum()} models auto-flagged cytotoxicity")
    print("\n[build] REVIEW 07_group_assignments.csv, then run with --analyze")
    return {"matrix": matrix, "col_index": col_index, "group_assignments": group_df}


def run_analyze(output_dir):
    """Analyze stage: Spearman matrix, top-N overlap, group summary (reads reviewed groups).

    Requires the build stage's ``07_score_matrix.parquet``, ``07_column_index.csv`` and the
    (reviewed) ``07_group_assignments.csv`` in ``output_dir``. Writes all correlation/overlap CSVs.
    """
    matrix = pd.read_parquet(os.path.join(output_dir, MATRIX_FILE))
    col_index = pd.read_csv(os.path.join(output_dir, "07_column_index.csv"))
    group_df = pd.read_csv(os.path.join(output_dir, "07_group_assignments.csv")).fillna(
        {"organisms": "", "cytotox_evidence": ""})

    print("[analyze] Spearman correlation ...")
    corr = spearman_matrix(matrix)
    order = cluster_order(corr)
    corr = corr.loc[order, order]
    corr.to_csv(os.path.join(output_dir, "07_spearman_corr.csv"))

    groups = build_groups(col_index, group_df)
    # Transparent membership: exactly which nodes landed in each focus group, for audit.
    membership = pd.DataFrame(
        [{"group": g, "node": n} for g, nodes in groups.items() for n in nodes])
    membership.to_csv(os.path.join(output_dir, "07_group_membership.csv"), index=False)
    summary = group_correlation_summary(corr, col_index, groups)
    summary.to_csv(os.path.join(output_dir, "07_group_correlation_summary.csv"), index=False)

    prob_nodes = col_index.loc[col_index["value_type"] == "probability", "node"].tolist()
    prob_nodes = [n for n in prob_nodes if n in matrix.columns]
    print(f"[analyze] top-N overlap over {len(prob_nodes)} probability nodes ...")
    for n in TOPN_CUTOFFS:
        jac = topn_overlap_matrix(matrix, prob_nodes, n)
        jac = jac.loc[[c for c in order if c in prob_nodes], [c for c in order if c in prob_nodes]]
        jac.to_csv(os.path.join(output_dir, f"07_topn_overlap_N{n}.csv"))
        print(f"  wrote 07_topn_overlap_N{n}.csv  ({len(prob_nodes)} nodes)")

    print(f"[analyze] wrote correlation + overlap summaries to {output_dir}")
    return {"corr": corr, "group_summary": summary}
