"""Step 14 — one score per organism, then the aggregated AUROC matrix.

Each organism's activity endpoints are collapsed into a SINGLE per-compound score before any AUROC is
computed, so a pathogen's weight in the figure comes from the pathogen rather than from how many
assays it happens to have (P. falciparum has 13 endpoints in this set, H. pylori has 1).

    rows (y)  15 organisms, each binarized at ACTIVITY_BINARIZE_TOP_N on its aggregate score
    cols (x)  the same 15 aggregates, then cytotoxicity (6), abx resemblance (3), physchem (3) = 27

The merge (see :data:`default.ORGANISM_MERGE_METHOD` / :data:`default.ORGANISM_MERGE_AGG`): the
organism's endpoint columns are scaled to percentile ranks within the full library, then averaged.
Five organisms have exactly one endpoint, so nothing is merged and their score IS that endpoint's
percentile rank — their row is not the same kind of quantity as E. coli's 11-endpoint mean.

Source scores come from step 07's parquet CACHE, re-scaled here, not from
``07_score_matrix_named_rankpct.csv``: that 6.8 GB CSV is stale (260 columns, predating two config
changes) and lacks several of these endpoints. Step 09 rebuilds from the parquet for the same reason.
"""

import os

import numpy as np
import pandas as pd

from default import (ACTIVITY_BINARIZE_TOP_N, AUROC_MATRIX_BLOCKS, ORGANISM_CLASS_ORDER,
                     ORGANISM_MERGE_AGG, ORGANISM_MERGE_METHOD)
from eval_correlations import scale_matrix
from eval_predictor_performance import (auroc_from_ranks, classify_predictor,
                                        pathogen_subset_endpoints, predictor_index,
                                        selected_targets)

KEY_COL = "key"


def bioactivity_order(selection_csv, pathogens_csv, available=None,
                      class_order=ORGANISM_CLASS_ORDER):
    """The activity endpoints that feed the merge, ordered, with their organism metadata.

    Same consensus-collapsed 15-pathogen set the per-endpoint version used; the ordering matters less
    now (endpoints are merged away) but the ``organism`` / ``organism_class`` columns drive both the
    grouping and the axis order.
    """
    targets = selected_targets(selection_csv)
    if available is not None:
        dropped = targets[~targets["target"].isin(available)]
        if len(dropped):
            print(f"  [auroc-matrix] {len(dropped)} selected endpoint(s) are absent from the "
                  f"step-07 cache and excluded: {dropped['target'].tolist()}")
        targets = targets[targets["target"].isin(available)].reset_index(drop=True)
    subset = pathogen_subset_endpoints(targets, pathogens_csv)

    sel = pd.read_csv(selection_csv)
    sel["target"] = sel["model_id"] + ":" + sel["column_name"]
    classes = sel.drop_duplicates("target").set_index("target")["organism_class"]
    subset = subset.copy()
    subset["organism_class"] = subset["target"].map(classes)

    unknown = sorted(set(subset["organism_class"].dropna()) - set(class_order))
    if unknown:
        raise ValueError(
            f"organism_class value(s) {unknown} are not in ORGANISM_CLASS_ORDER — add them there "
            "so their endpoints keep a defined position on the axis.")

    rank = {c: i for i, c in enumerate(class_order)}
    subset["_class"] = subset["organism_class"].map(rank)
    subset["_consensus"] = ~subset["is_consensus"]
    return (subset.sort_values(["_class", "organism", "_consensus", "model_id", "column_name"])
            .drop(columns=["_class", "_consensus"]).reset_index(drop=True))


def organism_order(endpoints, class_order=ORGANISM_CLASS_ORDER):
    """One row per organism: ``organism``, ``organism_class``, ``n_endpoints``, in axis order."""
    rank = {c: i for i, c in enumerate(class_order)}
    out = (endpoints.groupby(["organism_class", "organism"], sort=False)
           .size().reset_index(name="n_endpoints"))
    out["_class"] = out["organism_class"].map(rank)
    return out.sort_values(["_class", "organism"]).drop(columns="_class").reset_index(drop=True)


def organism_scores(parquet_path, endpoints, method=ORGANISM_MERGE_METHOD,
                    agg=ORGANISM_MERGE_AGG):
    """``(n_compounds x n_organisms)`` aggregate score matrix.

    Reads only this set's endpoint columns, scales them column-wise (percentile rank within the full
    library), then aggregates across each organism's columns. Scaling BEFORE aggregating is the whole
    point: raw scores from different models are on unrelated scales and averaging them directly would
    weight whichever endpoint happens to have the widest range.
    """
    cols = endpoints["target"].tolist()
    print(f"[organism-scores] reading {len(cols)} endpoint columns from "
          f"{os.path.basename(parquet_path)}")
    raw = pd.read_parquet(parquet_path, columns=cols)

    print(f"[organism-scores] scaling column-wise ({method})")
    scaled = scale_matrix(raw, method)
    del raw

    by_organism = endpoints.groupby("organism")["target"].apply(list)
    out = {}
    for organism, targets in by_organism.items():
        block = scaled[targets]
        out[organism] = getattr(block, agg)(axis=1)
        note = " (single endpoint — nothing merged)" if len(targets) == 1 else ""
        print(f"  {organism}: {agg} of {len(targets)} endpoint(s){note}")
    scores = pd.DataFrame(out, index=scaled.index)
    del scaled
    return scores


def _top_indices(values, top_n=ACTIVITY_BINARIZE_TOP_N):
    """Row positions of the ``top_n`` highest values; NaN can never enter the positive class."""
    v = np.where(np.isnan(values), -np.inf, values)
    return np.sort(np.argpartition(-v, top_n - 1)[:top_n])


def predictor_order(blocks=AUROC_MATRIX_BLOCKS):
    """The property columns as an ordered DataFrame: ``block``, ``family``, ``column_name``."""
    return pd.DataFrame([{"block": block, "family": family, "column_name": col}
                         for block, family, cols in blocks for col in cols])


def aggregated_matrix(scores, organisms, cols, property_csvs,
                      top_n=ACTIVITY_BINARIZE_TOP_N):
    """The ``n_organisms x (n_organisms + n_property)`` AUROC matrix, plus the resolved column axis.

    Rows are organisms binarized at ``top_n`` on their aggregate score; columns are the aggregate
    scores themselves followed by the property predictors. Uses the same Mann-Whitney rank-sum kernel
    as step 13 — rank each predictor once, then a ``top_n``-element gather per target.
    """
    from scipy.stats import rankdata

    order = organisms["organism"].tolist()
    tops = {o: _top_indices(scores[o].to_numpy(dtype=float), top_n) for o in order}
    n_total = len(scores)

    # --- Bioactivity block: each organism's aggregate as a predictor ---
    matrix = pd.DataFrame(index=order, dtype=float)
    for predictor in order:
        ranks = rankdata(scores[predictor].to_numpy(dtype=float), method="average")
        matrix[predictor] = [auroc_from_ranks(ranks, tops[t], n_total) for t in order]

    # --- Property blocks ---
    index = predictor_index(property_csvs)
    lookup = {(r.family, r.column_name): r for r in index.itertuples()}
    resolved, model_ids, sources = [], [], []
    for r in cols.itertuples():
        key = (r.family, r.column_name)
        if key not in lookup:
            raise ValueError(f"Predictor {r.family}/{r.column_name} is not in the step-11/12 "
                             "property blocks.")
        resolved.append(lookup[key].predictor)
        model_ids.append(lookup[key].model_id)
        sources.append(lookup[key].source)
    cols = cols.copy()
    cols["predictor"] = resolved
    cols["model_id"] = model_ids
    cols["source"] = sources

    for source, group in cols.groupby("source", sort=False):
        block = pd.read_csv(source, usecols=group["predictor"].tolist())
        for predictor in group["predictor"]:
            v = block[predictor].to_numpy(dtype=float)
            n_nan = int(np.isnan(v).sum())
            n_used, remap = n_total, None
            if n_nan:  # pairwise-complete, as in step 13
                valid = ~np.isnan(v)
                n_used = int(valid.sum())
                remap = np.full(n_total, -1, dtype=np.int64)
                remap[valid] = np.arange(n_used)
                v = v[valid]
            ptype = classify_predictor(v)
            if ptype != "continuous":
                raise ValueError(f"Property predictor {predictor} is {ptype}, not continuous — "
                                 "this matrix puts every column on one AUROC scale.")
            ranks = rankdata(v, method="average")
            vals = []
            for t in order:
                idx = tops[t]
                if remap is not None:
                    idx = remap[idx]
                    idx = idx[idx >= 0]
                vals.append(auroc_from_ranks(ranks, idx, n_used))
            matrix[predictor] = vals
        del block

    n_missing = int(matrix.isna().sum().sum())
    if n_missing:
        raise ValueError(f"{n_missing} matrix cell(s) have no value.")
    print(f"[auroc-matrix] {matrix.shape[0]} rows x {matrix.shape[1]} columns "
          f"({matrix.size} cells, none missing)")
    return matrix, cols


def _predictor_tops(cols, property_csvs, n_total, top_n=ACTIVITY_BINARIZE_TOP_N):
    """``{predictor: int64[top_n]}`` row positions of each property predictor's highest values."""
    out = {}
    for source, group in cols.groupby("source", sort=False):
        block = pd.read_csv(source, usecols=group["predictor"].tolist())
        for predictor in group["predictor"]:
            out[predictor] = _top_indices(block[predictor].to_numpy(dtype=float), top_n)
        del block
    return out


def overlap_matrix(scores, organisms, cols, property_csvs,
                   top_n=ACTIVITY_BINARIZE_TOP_N):
    """How many of the ROW organism's ``top_n`` actives fall in the COLUMN's own ``top_n``.

    A different question from the AUROC matrix — that asks whether a predictor RANKS an organism's
    actives highly across the whole library, this asks how many of the very same molecules it puts at
    the top. A predictor can do the first well without doing the second.

    The raw intersection count, not Jaccard: both sets have exactly ``top_n`` members, so
    ``J = i / (2 * top_n - i)`` is a monotone re-expression of the same number and orders the matrix
    identically. The count is the one a reader can act on.

    The measure is SYMMETRIC, so the bioactivity block is a symmetric matrix (the AUROC block is not),
    and the diagonal is ``top_n`` by construction. The figure blanks it; these values do not.
    """
    order = organisms["organism"].tolist()
    n_total = len(scores)
    tops = {o: set(_top_indices(scores[o].to_numpy(dtype=float), top_n)) for o in order}
    tops.update({p: set(v) for p, v in
                 _predictor_tops(cols, property_csvs, n_total, top_n).items()})

    column_keys = order + cols["predictor"].tolist()
    matrix = pd.DataFrame(index=order, columns=column_keys, dtype=float)
    for r in order:
        for c in column_keys:
            matrix.loc[r, c] = len(tops[r] & tops[c])

    off = matrix.values[~np.eye(len(order), M=len(column_keys), dtype=bool)]
    print(f"[overlap] {matrix.shape[0]} x {matrix.shape[1]}; off-diagonal {int(off.min())}-"
          f"{int(off.max())} of {top_n} (median {int(np.median(off))}; "
          f"chance for two random top-{top_n} sets ~ {top_n ** 2 / n_total:.1f})")
    return matrix.astype(int)


def axes_table(organisms, cols):
    """Row and column axis metadata, written alongside the matrix so the ordering and the annotation
    tracks are reproducible without re-deriving them."""
    r = organisms.copy()
    r.insert(0, "axis", "row")
    r.insert(1, "position", range(len(r)))
    r = r.rename(columns={"organism": "key"})
    r["block"] = "bioactivity"

    ca = r.copy()
    ca["axis"] = "column"
    cb = cols[["predictor", "block", "family", "column_name", "model_id"]].copy()
    cb.insert(0, "axis", "column")
    cb = cb.rename(columns={"predictor": "key"})
    cb.insert(1, "position", range(len(ca), len(ca) + len(cb)))
    return pd.concat([r, ca, cb], ignore_index=True)


def diagonal_check(matrix, organisms):
    """Each organism's aggregate against its OWN binarization — 1.0 only if rows and columns are
    aligned the same way, so this doubles as a transpose test."""
    return np.array([matrix.loc[o, o] for o in organisms], dtype=float)
