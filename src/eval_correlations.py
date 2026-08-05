"""Score-matrix engine — the full-library prediction matrices and the metrics computed on them.

Every annotation model was run on the SAME ~1.35M-compound reference library and staged by
``00_download_data.py`` as ``data/processed/annotation_preds_ref_library/{model_id}_{version}.csv``,
each with join columns ``key`` (compound hash) + ``input`` (SMILES) and model-specific output
columns. Because every file is aligned on ``key`` this is a clean rectangular matrix — no compound
alignment needed.

Which endpoints enter the analysis is decided by ONE file: the manually curated
``config/08_endpoint_selection.csv``, ``selected == Yes`` rows only. That is the single source of
truth — there is no automatic organism/model filter.

What this module provides, at the **output-column** level (each output column is a node):

  1. **matrix construction** — :func:`build_full_library_matrix` reads the raw prediction CSVs once
     (~15 GB) into a cacheable parquet; :func:`build_named_score_matrix` re-slices that cache and
     renames columns ``{pathogen_code}__{model_id}__{column_name}`` (step 07).
  2. **normalization** — :func:`scale_matrix` per column (z-score / rank-percentile) and
     :func:`row_normalize` per compound profile (L1 / L2). Independent and composable (step 07).
  3. **top-N Jaccard** — :func:`topn_jaccard_matrix`, full-library scale (step 08).
  4. **pathogen aggregation** — :func:`column_metric_pairs`, :func:`multi_column_pathogen_nodes`,
     :func:`pathogen_metric_boxes` and :func:`pathogen_metric_summary` turn a node x node metric
     matrix into same-pathogen vs. different-pathogen distributions (step 08).

Depends only on pandas / numpy — no project eval primitives, because the data layout
(``key``/``input``) differs from the ``smiles``/``bin`` contract in ``eval_common``.
"""

import glob
import os
import re

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

from default import CORR_CHUNK_SIZE, RANDOM_SEED, ROW_NORM_EPS

JOIN_COLS = ("key", "input")

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


def _read_model_columns(path, key_order, chunk_size, usecols=None):
    """Read one prediction file's output columns, filtered to ``key_order``, in chunks.

    ``usecols`` (optional) restricts which OUTPUT columns are parsed at all — worth passing when
    only a handful of a wide model's columns are needed (e.g. one file has 619), since the C parser
    then skips tokenizing/converting the rest. ``key`` is always read regardless. Returns a
    DataFrame indexed by ``key`` (reindexed to ``key_order``, so all models align) with one float
    column per requested output. Chunked so the 15.9 GB / 619-column file never loads whole.
    """
    key_set = set(key_order)
    read_cols = None if usecols is None else ["key", *usecols]
    parts = []
    for chunk in pd.read_csv(path, chunksize=chunk_size, usecols=read_cols):
        keep = chunk[chunk["key"].isin(key_set)]
        if len(keep):
            parts.append(keep)
    if not parts:
        return None
    df = pd.concat(parts, ignore_index=True).drop_duplicates(subset=["key"])
    out_cols = [c for c in df.columns if c not in JOIN_COLS]
    df = df.set_index("key")[out_cols].apply(pd.to_numeric, errors="coerce")
    return df.reindex(key_order).astype(np.float32)


def build_score_matrix(files, key_order, chunk_size=CORR_CHUNK_SIZE, verbose=True,
                       usecols_by_model=None):
    """Assemble the ``key × node`` score matrix from every model file.

    ``files`` is ``{model_id: (version, path)}``. Node columns are named ``"{model_id}:{output}"``.
    ``usecols_by_model`` (optional) is ``{model_id: [output_col, ...]}`` — when a model_id has an
    entry, only those output columns are read (see :func:`_read_model_columns`); models absent from
    the dict are read in full, so passing ``None`` keeps the original all-columns behaviour. Worth
    using when only a curated subset of a wide model's columns is needed. Non-numeric / failed
    outputs become NaN. Returns a float32 DataFrame indexed by ``key``.
    """
    usecols_by_model = usecols_by_model or {}
    blocks = []
    for i, (model_id, (ver, path)) in enumerate(sorted(files.items()), 1):
        if verbose:
            size_gb = os.path.getsize(path) / 1e9
            print(f"  [{i}/{len(files)}] {model_id} {ver} ({size_gb:.1f} GB) ...", flush=True)
        df = _read_model_columns(path, key_order, chunk_size,
                                 usecols=usecols_by_model.get(model_id))
        if df is None or df.shape[1] == 0:
            print(f"    [skip] {model_id}: no usable output columns")
            continue
        df.columns = [f"{model_id}:{c}" for c in df.columns]
        blocks.append(df)
    matrix = pd.concat(blocks, axis=1)
    matrix.index.name = "key"
    return matrix


def build_full_library_matrix(pred_dir, endpoint_selection_path, seed=RANDOM_SEED):
    """Build the FULL ~1.35M-compound ``key x node`` matrix, restricted to every (model_id, column)
    endpoint referenced anywhere in ``endpoint_selection_path`` — both ``selected == Yes`` AND
    ``No`` rows, so the cache stays valid if the curation changes later (e.g. Mtb re-enabled).

    This reads the FULL key set (:func:`sample_keys` with ``n=None``) from whichever needed
    model's file is available first,
    and only the specific columns each model needs (:func:`build_score_matrix`'s
    ``usecols_by_model``) rather than every column of every file — required at this scale, since
    some of these models have dozens of columns and the raw files run into the GB range each.

    A model referenced in the CSV with no raw prediction file in ``pred_dir`` (e.g. Mtb's
    ``eos9ivc``, not yet downloaded) is skipped with a logged message, exactly like
    :func:`build_score_matrix` already does for unusable columns — never an error, since the
    curated CSV may reference models whose data isn't staged yet.
    """
    sel = pd.read_csv(endpoint_selection_path)
    needed = sel.groupby("model_id")["column_name"].apply(list).to_dict()
    files = latest_version_files(pred_dir)

    missing_models = [m for m in needed if m not in files]
    for m in missing_models:
        print(f"  [full-library] skip {m}: no raw prediction file in {pred_dir}")
    usable = {m: files[m] for m in needed if m in files}
    if not usable:
        raise FileNotFoundError(f"None of the {len(needed)} models in {endpoint_selection_path} "
                                f"have a raw prediction file in {pred_dir}")

    ref_path = next(iter(usable.values()))[1]
    print(f"  [full-library] drawing full key set from {ref_path} ...")
    key_order = sample_keys(ref_path, n=None, seed=seed)
    print(f"  [full-library] {len(key_order)} compounds x {len(usable)} models "
          f"({sum(len(c) for c in needed.values())} requested columns, "
          f"{len(missing_models)} model(s) skipped)")

    matrix = build_score_matrix(usable, key_order, usecols_by_model=needed)
    return matrix


def _pathogen_code(organism, org_to_code):
    """Short pathogen code for a column name: the curated code if ``organism`` is in
    ``org_to_code`` (from ``config/pathogens_of_interest.csv``), otherwise a mechanical fallback —
    first letter of genus + species epithet, lowercased (``"Bacteroides caccae"`` -> ``"bcaccae"``)
    — matching how the curated codes already look. Single-word names (e.g. ``"Campylobacter spp"``
    is two words, but a hypothetical one-word organism) are lowercased as-is.
    """
    if organism in org_to_code:
        return org_to_code[organism]
    parts = organism.replace(".", "").split()
    return parts[0].lower() if len(parts) == 1 else (parts[0][0] + parts[-1]).lower()


def build_named_score_matrix(pred_dir, endpoint_selection_path, pathogens_of_interest_path,
                             full_matrix_cache_path=None, seed=RANDOM_SEED):
    """The ``selected == Yes`` full-library matrix, columns renamed
    ``"{pathogen_code}__{model_id}__{column_name}"`` so the pathogen and source model are readable
    straight off the column name, not just recoverable by joining back to
    ``endpoint_selection_path``.

    Pathogen codes come from ``pathogens_of_interest_path`` (``config/pathogens_of_interest.csv``,
    the existing curated organism -> code mapping) where available, and a mechanical fallback
    (:func:`_pathogen_code`) otherwise — most of the ``selected == Yes`` set is single-model
    gut-microbiome organisms outside the 15-pathogen priority list, so most codes ARE
    fallback-generated; every fallback is printed for audit, and a code collision raises rather
    than silently overwriting a column.

    If ``full_matrix_cache_path`` (typically ``07_score_matrix_full.parquet`` from
    :func:`build_full_library_matrix`, which already covers every ``selected`` AND ``No`` row) has
    every needed column, this reads just those columns from the cache instead of re-reading the
    raw prediction files — checked via the parquet's own schema (no full-file load needed to check).
    Falls back to a full :func:`build_full_library_matrix` build otherwise.
    """
    sel_all = pd.read_csv(endpoint_selection_path)
    sel = sel_all[sel_all["selected"] == "Yes"].copy()

    pathogens = pd.read_csv(pathogens_of_interest_path)
    org_to_code = dict(zip(pathogens["pathogen"], pathogens["code"]))

    unique_orgs = sorted(sel["organism"].unique())
    codes = {org: _pathogen_code(org, org_to_code) for org in unique_orgs}
    fallback_used = [(org, c) for org, c in codes.items() if org not in org_to_code]
    code_counts = pd.Series(list(codes.values())).value_counts()
    collisions = code_counts[code_counts > 1]
    if len(collisions):
        colliding_orgs = [o for o, c in codes.items() if c in collisions.index]
        raise ValueError(f"Pathogen code collision(s) {list(collisions.index)} across organisms "
                         f"{colliding_orgs} — resolve in {pathogens_of_interest_path} or "
                         "_pathogen_code's fallback before proceeding.")
    print(f"[named-matrix] {len(unique_orgs)} organisms: {len(unique_orgs) - len(fallback_used)} "
          f"from {pathogens_of_interest_path}, {len(fallback_used)} fallback-generated")
    for org, c in fallback_used:
        print(f"    fallback: {org!r} -> {c!r}")

    old_cols = [f"{r.model_id}:{r.column_name}" for r in sel.itertuples()]
    rename_map = {f"{r.model_id}:{r.column_name}": f"{codes[r.organism]}__{r.model_id}__{r.column_name}"
                 for r in sel.itertuples()}

    if full_matrix_cache_path and os.path.exists(full_matrix_cache_path):
        schema_cols = set(pq.ParquetFile(full_matrix_cache_path).schema.names)
        if set(old_cols) <= schema_cols:
            print(f"[named-matrix] reusing cached {full_matrix_cache_path} "
                 f"({len(old_cols)} of its columns)")
            matrix = pd.read_parquet(full_matrix_cache_path, columns=old_cols)
            return matrix.rename(columns=rename_map)
        print(f"[named-matrix] cache at {full_matrix_cache_path} is missing "
             f"{len(set(old_cols) - schema_cols)} needed column(s) — building from raw files")

    full = build_full_library_matrix(pred_dir, endpoint_selection_path, seed=seed)
    return full[old_cols].rename(columns=rename_map)


#: Column-wise scaling methods for :func:`scale_matrix`.
SCALE_METHODS = ("zscore", "rank_pct")


def scale_matrix(matrix, method):
    """Column-wise scaling of a ``key x node`` score matrix, computed over every row given (the
    full library, if that's what ``matrix`` is) — same shape/index/columns as the input.

    ``method="zscore"``: ``(x - mean) / std`` per column (``ddof=0``, population std). A
    zero-variance column (constant value) would divide by zero; these are detected up front and
    printed as a warning (their output is left as NaN, not silently zeroed or dropped) rather than
    failing partway through.

    ``method="rank_pct"``: each value replaced by its percentile rank within its own column
    (``pandas.DataFrame.rank(pct=True)``, average rank for ties) — bounded ``[0, 1]``, robust to
    outliers and to columns on very different native scales, the same rank-based idea behind
    every other metric in this analysis (Spearman, AUROC, Jaccard top/bottom-X).
    """
    if method not in SCALE_METHODS:
        raise ValueError(f"method must be one of {SCALE_METHODS}, got {method!r}")
    if method == "zscore":
        std = matrix.std(axis=0, ddof=0)
        zero_var = std[std == 0].index.tolist()
        if zero_var:
            print(f"[scale] WARNING: {len(zero_var)} zero-variance column(s) -> NaN after "
                 f"z-scoring: {zero_var}")
        return (matrix - matrix.mean(axis=0)) / std
    return matrix.rank(pct=True)


ROW_NORMS = ("l1", "l2")


def row_normalize(matrix, norm):
    """Row-wise normalization of a ``key x node`` score matrix — same shape/index/columns.

    Complementary to (and composable with) `scale_matrix`: that puts every *column* on a common
    footing, this rescales each *compound's* profile across nodes so profiles are compared by
    shape rather than by overall magnitude.

    ``norm="l1"``: divide each row by the sum of its absolute values, so every row sums to 1.
    Best suited to non-negative input (e.g. a ``rank_pct``-scaled matrix), where the result reads
    as a compositional vector — each node's relative share of that compound's total activity.

    ``norm="l2"``: divide each row by its Euclidean norm, so every row lands on the unit sphere.
    The conventional choice for signed, mean-centred input (e.g. a ``zscore``-scaled matrix), and
    the normalization cosine similarity between two profiles is built on.

    Rows whose norm is zero or near-zero (``< ROW_NORM_EPS``) are counted and printed as a
    warning: they divide to inf/NaN rather than being silently dropped or clipped. Note also that
    a compound with a near-flat profile is *amplified* by row normalization — expected behaviour,
    but worth remembering when reading such rows.

    Pre-existing NaNs in the input (a compound a model failed to score) are likewise counted and
    reported, not imputed: the norm is taken over that row's non-null values, so the row is unit
    norm across the endpoints it *does* have and the missing cell stays NaN.
    """
    if norm not in ROW_NORMS:
        raise ValueError(f"norm must be one of {ROW_NORMS}, got {norm!r}")
    n_missing_by_col = matrix.isna().sum()
    n_missing = int(n_missing_by_col.sum())
    if n_missing:
        cols = n_missing_by_col[n_missing_by_col > 0]
        print(f"[rownorm] {n_missing} pre-existing NaN cell(s) in "
              f"{len(cols)} column(s) — norms taken over each row's non-null values, NaNs kept: "
              f"{cols.to_dict()}")
    if norm == "l1":
        norms = matrix.abs().sum(axis=1)
    else:
        norms = np.sqrt((matrix ** 2).sum(axis=1))
    n_degenerate = int((norms < ROW_NORM_EPS).sum())
    if n_degenerate:
        print(f"[rownorm] WARNING: {n_degenerate} row(s) with {norm} norm < {ROW_NORM_EPS:g} "
              "-> inf/NaN after normalization (left as-is, not dropped)")
    return matrix.div(norms, axis=0)


# --------------------------------------------------------------------------- #
# Top-N overlap                                                                #
# --------------------------------------------------------------------------- #
def topn_jaccard_matrix(matrix, n, verbose=False):
    """Pairwise Jaccard overlap of each column's top-``n`` compounds — full-library scale.

    The membership mask is filled one column at a time with ``np.argpartition`` (an O(rows) partial
    selection) on a float32 view, rather than materialising a float64 copy plus a complete
    ``np.argsort`` of every column. At 1.35M x 260 that is ~1.8 GB peak versus ~9 GB.

    This is the plain top-vs-top Jaccard the name says: no maximisation over bottom-``n`` or
    cross (top-vs-bottom) set combinations, and no direction flipping.

    NaN values never enter a top-``n`` set (a compound a model failed to score is not a hit). A
    column with fewer than ``n`` finite values therefore contributes a smaller set, and the union
    term ``|A| + |B| - |A n B|`` accounts for that; pairs with an empty union give NaN.
    """
    nodes = list(matrix.columns)
    values = matrix.to_numpy(dtype=np.float32, copy=False)
    n_rows = values.shape[0]
    n = min(n, n_rows)
    mask = np.zeros(values.shape, dtype=bool)
    for j in range(values.shape[1]):
        col = values[:, j]
        finite = np.isfinite(col)
        if finite.all():
            idx = np.argpartition(-col, n - 1)[:n]
        else:
            # Restrict to finite entries first, then map back to original row positions.
            where = np.flatnonzero(finite)
            sub = col[where]
            k = min(n, sub.shape[0])
            idx = where[np.argpartition(-sub, k - 1)[:k]] if k else where[:0]
        mask[idx, j] = True
        if verbose and (j + 1) % 50 == 0:
            print(f"  [topn-jaccard] top-{n} mask: {j + 1}/{values.shape[1]} columns")

    # One float32 view of the mask, reused for both matmul operands (two separate `astype` calls
    # would hold two full-size copies at once).
    fmask = mask.astype(np.float32)
    inter = fmask.T @ fmask
    sizes = mask.sum(axis=0).astype(np.float64)
    del fmask
    union = sizes[:, None] + sizes[None, :] - inter
    with np.errstate(divide="ignore", invalid="ignore"):
        jac = np.where(union > 0, inter / union, np.nan)
    return pd.DataFrame(jac, index=nodes, columns=nodes)

def parse_named_column(name):
    """Split a ``"{pathogen}__{model_id}__{column_name}"`` column name into its three parts.

    The naming convention written by :func:`build_named_score_matrix` — the pathogen and source
    model are readable straight off the column name, so no metadata join is needed here.
    """
    parts = name.split("__", 2)
    if len(parts) != 3:
        raise ValueError(f"column {name!r} is not '{{pathogen}}__{{model_id}}__{{column_name}}'")
    return parts[0], parts[1], parts[2]


def column_metric_pairs(metric_matrix):
    """Long-format, DIRECTED view of a square node x node metric matrix.

    One row per (node, partner) ordered pair, so every node carries its complete set of partners
    and per-node statistics need no symmetry bookkeeping. Self-pairs are dropped. Columns:
    ``node, partner, pathogen, partner_pathogen, model_id, partner_model_id, category, same_model,
    value`` — ``category`` is ``same_pathogen`` / ``different_pathogen`` read off the column-name
    prefixes.

    ``same_model`` is carried but NOT used to filter: two columns of one model agreeing says little
    about pathogen specificity, so the flag lets a caller report the same-pathogen statistic both
    with and without those pairs instead of silently picking one.
    """
    nodes = list(metric_matrix.columns)
    parsed = {n: parse_named_column(n) for n in nodes}
    values = metric_matrix.to_numpy(dtype=np.float64)
    idx = np.arange(len(nodes))
    a, b = np.meshgrid(idx, idx, indexing="ij")
    off = a != b
    a, b = a[off], b[off]
    pathogen = np.array([parsed[n][0] for n in nodes])
    model_id = np.array([parsed[n][1] for n in nodes])
    same_pathogen = pathogen[a] == pathogen[b]
    return pd.DataFrame({
        "node": np.array(nodes)[a], "partner": np.array(nodes)[b],
        "pathogen": pathogen[a], "partner_pathogen": pathogen[b],
        "model_id": model_id[a], "partner_model_id": model_id[b],
        "category": np.where(same_pathogen, "same_pathogen", "different_pathogen"),
        "same_model": model_id[a] == model_id[b],
        "value": values[a, b],
    })

def multi_column_pathogen_nodes(metric_matrix, min_columns=2):
    """Nodes of the pathogens contributing at least ``min_columns`` output columns.

    A pathogen represented by a single column has NO same-pathogen partner, so it can never carry
    the same-vs-different comparison. Dropping such pathogens removes them from the analysis
    ENTIRELY — they stop being different-pathogen partners for the pathogens that remain, not just
    stop having a box of their own.
    """
    pathogen = pd.Series({n: parse_named_column(n)[0] for n in metric_matrix.columns})
    keep = pathogen.value_counts()
    keep = set(keep[keep >= min_columns].index)
    return [n for n in metric_matrix.columns if pathogen[n] in keep]


def pathogen_metric_boxes(pairs):
    """Per-PATHOGEN same/different value arrays, aggregated from a :func:`column_metric_pairs` frame.

    For pathogen P: ``same`` is every unordered pair of P's own columns (deduplicated — the pairs
    frame is directed, so each within-pathogen pair appears twice and would otherwise be
    double-counted), ``diff`` is every pair from one of P's columns to a column of some OTHER
    pathogen (already counted once each, since only P-side rows are taken).

    Returns ``{pathogen: {"same": array, "diff": array, "same_excl_same_model": array}}``.
    ``same_excl_same_model`` repeats ``same`` without pairs whose two columns come from the same
    model — the check on whether a pathogen's internal agreement is really one model agreeing with
    itself.
    """
    out = {}
    for pathogen, d in pairs.groupby("pathogen", sort=False):
        same = d[d["category"] == "same_pathogen"]
        same = same[same["node"] < same["partner"]]  # one row per unordered pair
        diff = d[d["category"] == "different_pathogen"]
        out[pathogen] = {
            "same": same["value"].to_numpy(float),
            "same_excl_same_model": same.loc[~same["same_model"], "value"].to_numpy(float),
            "diff": diff["value"].to_numpy(float),
        }
    return out


def pathogen_metric_summary(boxes, n_columns):
    """Tidy per-pathogen summary of :func:`pathogen_metric_boxes`, ordered by specificity.

    ``specificity`` is ``same_median - diff_median`` — the same definition step 08's breakdown uses —
    and the frame is sorted by it descending, which is also the row order used by the figure. A
    pathogen whose specificity is negative agrees MORE with other pathogens' columns than with its
    own; that is reported, not filtered.
    """
    rows = []
    for pathogen, b in boxes.items():
        rows.append({
            "pathogen": pathogen, "n_columns": int(n_columns[pathogen]),
            "n_same_pairs": len(b["same"]),
            "n_same_pairs_excl_same_model": len(b["same_excl_same_model"]),
            "n_diff_pairs": len(b["diff"]),
            "same_median": np.median(b["same"]) if len(b["same"]) else np.nan,
            "same_median_excl_same_model": (np.median(b["same_excl_same_model"])
                                            if len(b["same_excl_same_model"]) else np.nan),
            "diff_median": np.median(b["diff"]) if len(b["diff"]) else np.nan,
        })
    out = pd.DataFrame(rows)
    out["specificity"] = out["same_median"] - out["diff_median"]
    return out.sort_values("specificity", ascending=False).reset_index(drop=True)
