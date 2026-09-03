"""Step 10 — one score per organism, then the aggregated AUROC matrix.

Each organism's activity endpoints are collapsed into a SINGLE per-compound score before any AUROC is
computed, so a pathogen's weight in the figure comes from the pathogen rather than from how many
assays it happens to have (P. falciparum has 13 endpoints in this set, H. pylori has 1).

    rows (y)  15 organisms, each binarized at ACTIVITY_BINARIZE_TOP_N on its aggregate score
    cols (x)  the same 15 aggregates, then cytotoxicity (1) and abx resemblance (1) = 17

Cytotoxicity and abx resemblance are each a single MERGED rank-sum predictor (2026-09-02,
user-directed — see :func:`merged_predictor_scores`): every raw column in the family is percentile-
ranked, then summed per compound, collapsing 6 cytotoxicity columns and 3 abx-resemblance columns to
one each. Physchem (mw/clogp/tpsa) was dropped rather than merged — those three measure unrelated
quantities from each other and from "how cytotoxic/antibiotic-like", so summing their ranks together
would not be a meaningful number the way a within-family sum is.

The merge (see :data:`default.ORGANISM_MERGE_METHOD` / :data:`default.ORGANISM_MERGE_AGG`): the
organism's endpoint columns are scaled to percentile ranks within the full library, then averaged.
Five organisms have exactly one endpoint, so nothing is merged and their score IS that endpoint's
percentile rank — their row is not the same kind of quantity as E. coli's 12-endpoint mean.

Source scores come from step 07's parquet CACHE, re-scaled here, not from
``07_score_matrix_named_rankpct.csv``: that 6.8 GB CSV is stale (260 columns, predating two config
changes) and lacks several of these endpoints. Step 09 rebuilds from the parquet for the same reason.
"""

import os

import numpy as np
import pandas as pd

from default import (ACTIVITY_BINARIZE_TOP_N, AUROC_MATRIX_BLOCKS, MERGED_PREDICTOR_GROUPS,
                     ORGANISM_CLASS_ORDER, ORGANISM_MERGE_AGG, ORGANISM_MERGE_METHOD,
                     PREDICTOR_MERGE_AGG)
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


#: Coarsest-to-finest rank columns read from ``config/organism_taxonomy.csv``. The number of these
#: (5) is also the maximum possible :func:`phylogeny_class_linkages` merge height — two organisms
#: sharing none of them (different phylum) merge at height 5, the tree's root.
TAXONOMY_RANK_COLUMNS = ["phylum", "taxonomic_class", "order", "family", "genus"]


def _cartesian_linkage_from_order(gaps):
    """A valid scipy linkage matrix for ``len(gaps) + 1`` leaves, given ONLY in ``ORDER``.

    ``gaps[i]`` is the merge height between leaves ``i`` and ``i + 1`` — implicitly, the height at
    which the whole run between two leaves would be joined, valid only when ``gaps`` comes from a
    genuine ultrametric (see :func:`phylogeny_class_linkages`: two organisms at the same rank-depth
    are always exactly equidistant, so the "gap to the right neighbour" fully determines every
    cross-cluster distance — no pairwise matrix or clustering ambiguity to resolve).

    Repeatedly merges the two ADJACENT clusters with the smallest remaining gap — never a distant
    pair — so the leaf order given in is preserved exactly, with no tie-breaking freedom of the kind
    ``scipy.cluster.hierarchy.linkage`` has on a matrix with many exactly-equal distances (e.g. every
    pair among 3 organisms that share a family but not a genus). Always picking the global minimum
    next also means merges are emitted in non-decreasing height order, which is what scipy expects.
    """
    n = len(gaps) + 1
    nodes = [(i, 1) for i in range(n)]     # (cluster id, size), in current left-to-right order
    gaps = list(gaps)
    next_id = n
    z = []
    while len(nodes) > 1:
        i = min(range(len(gaps)), key=lambda k: gaps[k])
        (id1, s1), (id2, s2) = nodes[i], nodes[i + 1]
        z.append([id1, id2, gaps[i], s1 + s2])
        nodes[i:i + 2] = [(next_id, s1 + s2)]
        next_id += 1
        del gaps[i]
    return np.array(z, dtype=float)


#: Tied-sibling display order overrides, as ``(first, second)`` — ``first`` should appear immediately
#: before ``second``. A pair belongs here only when the taxonomy itself gives NO reason to prefer
#: either order: the two merge with each other before either merges with anything else, so swapping
#: them changes no tree height or topology, only which one is drawn "first". Same convention as
#: ``plotting_colors.PATHOGEN_HUE_SWAPS`` (documented, minimal, applied after the general algorithm).
#:
#: (2026-09-02, user-directed) *S. pneumoniae* / *E. faecium*: both Bacillota/Bacilli/Lactobacillales,
#: different family, so they merge with each other at height 2 — and *S. aureus* (Bacillales) is
#: EXACTLY as distant from each of them (height 3 to both), confirmed by direct computation, so its
#: position is unaffected by which of the two sits "on top". The plain alphabetical tie-break (by
#: family: Enterococcaceae < Streptococcaceae) put *E. faecium* first; swapped here for readability.
PHYLO_SIBLING_SWAPS = (("Streptococcus pneumoniae", "Enterococcus faecium"),)


def _apply_sibling_swaps(members, swaps=PHYLO_SIBLING_SWAPS):
    """``members`` with each :data:`PHYLO_SIBLING_SWAPS` pair placed ``(first, second)`` if both are
    present and currently adjacent in the other order — a no-op otherwise, so a pair missing from
    this class or already in the requested order is left alone.
    """
    members = list(members)
    for first, second in swaps:
        if first in members and second in members:
            i, j = members.index(first), members.index(second)
            if j == i - 1:
                members[j], members[i] = members[i], members[j]
    return members


def phylogeny_class_linkages(organisms, taxonomy_csv, class_order=ORGANISM_CLASS_ORDER,
                             ranks=TAXONOMY_RANK_COLUMNS):
    """``{organism_class: (Z, members)}`` for every class with 2+ organisms, built directly from the
    real NCBI-taxonomy lineage in ``taxonomy_csv`` — not fitted to any AUROC value (the alternative
    "cluster on the AUROC row profile" ordering was tried and dropped 2026-09-02 in favour of this
    one).

    ``members`` is sorted by the lineage tuple ``(phylum, taxonomic_class, order, family, genus,
    organism)`` — the SAME sort :func:`phylogeny_organism_order` uses — so the dendrogram
    (:func:`plots_auroc_matrix.save_dendrogram_figure`) and the ``phylo`` heatmap's row order can
    never disagree on which organism sits where.

    Merge height = ``len(ranks)`` minus how many of the 5 ranks the two sides share, i.e. how many
    levels up you must climb to find a common ancestor: sharing everything down to family merges at
    height 1, sharing only phylum merges at height 4, sharing nothing merges at height 5. This is an
    exact ultrametric by construction — a real nested classification, not a fitted distance — so
    :func:`_cartesian_linkage_from_order` assembles the one tree it implies directly, rather than
    approximating it with ``scipy.cluster.hierarchy.linkage`` on a full pairwise matrix, which would
    face unnecessary tie-breaking on the many exactly-equal distances a coarse, 5-rank taxonomy
    produces (e.g. every pair among the three Enterobacteriaceae genera here is equidistant).
    """
    tax = pd.read_csv(taxonomy_csv).set_index("organism")
    missing = sorted(set(organisms["organism"]) - set(tax.index))
    if missing:
        raise ValueError(f"organism(s) {missing} are not in {taxonomy_csv} — add their lineage so "
                          "the phylogeny-based order can place them.")

    by_class = organisms.groupby("organism_class", sort=False)["organism"].apply(list)
    out = {}
    for cls in class_order:
        members = sorted(by_class.get(cls, []),
                         key=lambda o: tuple(tax.loc[o, ranks]) + (o,))
        members = _apply_sibling_swaps(members)
        if len(members) < 2:
            continue
        gaps = []
        for a, b in zip(members[:-1], members[1:]):
            shared = 0
            for col in ranks:
                if tax.loc[a, col] != tax.loc[b, col]:
                    break
                shared += 1
            gaps.append(len(ranks) - shared)
        out[cls] = (_cartesian_linkage_from_order(gaps), members)
    return out


def phylogeny_organism_order(organisms, taxonomy_csv, class_order=ORGANISM_CLASS_ORDER):
    """Row order: within each class, sort by real NCBI-taxonomy lineage.

    Class blocks are concatenated in ``class_order`` sequence — only the sequence WITHIN each class
    changes. Derived from :func:`phylogeny_class_linkages` so this order and the dendrogram
    diagnostic read the same tree; a class with fewer than 2 members has no tree and keeps its
    (single-member, trivial) order.
    """
    linkages = phylogeny_class_linkages(organisms, taxonomy_csv, class_order)
    by_class = organisms.groupby("organism_class", sort=False)["organism"].apply(list)
    out = []
    for cls in class_order:
        if cls not in linkages:
            out.extend(by_class.get(cls, []))
            continue
        _, members = linkages[cls]
        out.extend(members)
    return out


def reorder_bioactivity_axes(matrix, organisms, new_order):
    """``matrix``/``organisms`` re-indexed to ``new_order`` — a permutation of the same organisms.

    Splits the columns into the bioactivity block (the organism names themselves, reindexed to
    ``new_order``) and the property blocks (everything else, left in their original order, since
    they are not organism-linked). The diagonal stays each organism against its own binarization, so
    :func:`diagonal_check` still passes unchanged.
    """
    bio_cols = organisms["organism"].tolist()
    prop_cols = [c for c in matrix.columns if c not in set(bio_cols)]
    new_matrix = matrix.loc[new_order, new_order + prop_cols]
    new_organisms = organisms.set_index("organism").loc[new_order].reset_index()
    return new_matrix, new_organisms


def organism_scores(parquet_path, endpoints, method=ORGANISM_MERGE_METHOD,
                    agg=ORGANISM_MERGE_AGG, row_mask=None):
    """``(n_compounds x n_organisms)`` aggregate score matrix.

    Reads only this set's endpoint columns, scales them column-wise (percentile rank within the
    library), then aggregates across each organism's columns. Scaling BEFORE aggregating is the whole
    point: raw scores from different models are on unrelated scales and averaging them directly would
    weight whichever endpoint happens to have the widest range.

    ``row_mask`` is an optional boolean array over the FULL library, for running the matrix on a
    subset (see step 10's non-abx section). It is applied **before** the scaling, so the
    percentile ranks are recomputed within the subset rather than inherited from the full library —
    the two are different analyses, and this is the one that asks "what would the matrix look like if
    the library had never contained these compounds". ``scale_matrix`` ranks over every row it is
    given, so the choice is made here and nowhere else.

    **Every other reader of a positional source must be given the same mask**, or the property block
    will be indexed with positions that no longer mean what they did — see :func:`aggregated_matrix`.
    """
    cols = endpoints["target"].tolist()
    print(f"[organism-scores] reading {len(cols)} endpoint columns from "
          f"{os.path.basename(parquet_path)}")
    raw = pd.read_parquet(parquet_path, columns=cols)
    if row_mask is not None:
        raw = raw[row_mask]
        print(f"[organism-scores] row mask applied BEFORE scaling: {len(raw):,} of "
              f"{len(row_mask):,} compounds ({100 * len(raw) / len(row_mask):.2f}%)")

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


def _apply_row_mask(block, row_mask, source):
    """``block`` restricted to ``row_mask``, with the length checked first.

    The length check is the whole value of this helper. Every consumer of a property CSV in this
    module indexes it by row POSITION, so a mask built against a different row count would not raise
    — it would quietly select the wrong compounds and produce a full, plausible, wrong matrix.
    """
    if row_mask is None:
        return block
    if len(block) != len(row_mask):
        raise ValueError(
            f"{os.path.basename(source)} has {len(block)} rows but the row mask covers "
            f"{len(row_mask)} — they are compared positionally, so a mismatch would silently "
            "select the wrong compounds. Rebuild both from the same reference library.")
    return block[row_mask].reset_index(drop=True)


def _top_indices(values, top_n=ACTIVITY_BINARIZE_TOP_N):
    """Row positions of the ``top_n`` highest values; NaN can never enter the positive class."""
    v = np.where(np.isnan(values), -np.inf, values)
    if top_n > len(v):
        raise ValueError(f"top_n={top_n} exceeds the {len(v)} rows available.")
    return np.sort(np.argpartition(-v, top_n - 1)[:top_n])


def predictor_order(blocks=AUROC_MATRIX_BLOCKS):
    """The property columns as an ordered DataFrame: ``block``, ``family``, ``column_name``."""
    return pd.DataFrame([{"block": block, "family": family, "column_name": col}
                         for block, family, cols in blocks for col in cols])


def merged_predictor_scores(property_csvs, groups=MERGED_PREDICTOR_GROUPS, agg=PREDICTOR_MERGE_AGG,
                            row_mask=None):
    """One merged predictor column per family in ``groups`` (2026-09-02, user-directed): each raw
    column is percentile-ranked, then the ranks are aggregated (SUM by default) across that family's
    raw columns — 6 cytotoxicity columns and 3 abx-resemblance columns collapse to 1 each.

    ``row_mask`` is applied BEFORE ranking, same convention as :func:`organism_scores` and
    :func:`aggregated_matrix`: for a masked run (step 10's non-abx section) each merged
    predictor's ranks are relative to the surviving subset, not inherited from the full library.

    ``property_csvs`` must be the RAW per-family matrices (step 08) — this is what the merge is
    computed FROM, not the merged output itself.

    Returns a DataFrame with one column per family, named ``"{family}__merged__rank_sum"`` — fits
    :func:`eval_predictor_performance.predictor_index`'s ``{family}__{model_id}__{column}`` parsing
    unchanged, with ``"merged"`` standing in for a model_id since the column is no longer tied to one
    model. **Summed with ``skipna=False``**: a row missing any of a family's raw columns (66 of
    1,355,109 do, in the abx block — see step 08) gets NaN in the merged column rather than a
    silently-partial sum over fewer inputs than the rest of the library — the same "never impute,
    never silently drop" rule as everywhere else in this pipeline. Downstream, :func:`aggregated_matrix`
    and :func:`_top_indices` already handle a NaN-bearing property column exactly this way (pairwise-
    complete AUROC, NaN excluded from ever entering a top-N set), so no special-casing is needed here.
    """
    index = predictor_index(property_csvs)
    lookup = {(r.family, r.column_name): r for r in index.itertuples()}
    out = {}
    for family, raw_columns in groups.items():
        by_source = {}
        for col in raw_columns:
            key = (family, col)
            if key not in lookup:
                raise ValueError(f"Merge input {family}/{col} not found in {property_csvs}.")
            r = lookup[key]
            by_source.setdefault(r.source, []).append(r.predictor)
        ranked = []
        for source, predictors in by_source.items():
            block = pd.read_csv(source, usecols=predictors)
            block = _apply_row_mask(block, row_mask, source)
            ranked.append(scale_matrix(block, "rank_pct"))
        combined = pd.concat(ranked, axis=1)
        merged = getattr(combined, agg)(axis=1, skipna=False)
        n_nan = int(merged.isna().sum())
        out[f"{family}__merged__rank_sum"] = merged
        print(f"[merged-predictors] {family}: {agg} of {len(raw_columns)} rank_pct-scaled column(s) "
              f"-> {family}__merged__rank_sum, range [{merged.min():.3f}, {merged.max():.3f}]"
              + (f", {n_nan} row(s) NaN (missing input)" if n_nan else ""))
    return pd.DataFrame(out)


def aggregated_matrix(scores, organisms, cols, property_csvs,
                      top_n=ACTIVITY_BINARIZE_TOP_N, row_mask=None, return_top_indices=False):
    """The ``n_organisms x (n_organisms + n_property)`` AUROC matrix, plus the resolved column axis.

    Rows are organisms binarized at ``top_n`` on their aggregate score; columns are the aggregate
    scores themselves followed by the property predictors. Uses the same Mann-Whitney rank-sum kernel
    as step 14 — rank each predictor once, then a ``top_n``-element gather per target.

    ``row_mask`` MUST be the same boolean array that was passed to :func:`organism_scores`. The
    property CSVs are indexed **positionally** by row positions derived from ``scores``; if ``scores``
    is a masked subset and the CSVs are not, every property-block AUROC is silently wrong — and the
    diagonal check downstream would still pass, because it only touches the bioactivity block. That
    is the one corruption none of this module's guards catches, which is why the mask is a parameter
    here rather than something the caller applies beforehand.

    ``return_top_indices=True`` additionally returns ``{predictor: positions}`` — the same
    ``top_n`` top-index computation :func:`predictor_tops` would do, piggybacked onto the property
    block this function already reads, so a caller that also needs :func:`overlap_matrix` for the
    SAME ``top_n`` and property CSVs does not have to read them a second time.
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
            raise ValueError(f"Predictor {r.family}/{r.column_name} is not in the step-08 "
                             "property blocks.")
        resolved.append(lookup[key].predictor)
        model_ids.append(lookup[key].model_id)
        sources.append(lookup[key].source)
    cols = cols.copy()
    cols["predictor"] = resolved
    cols["model_id"] = model_ids
    cols["source"] = sources

    predictor_top_indices = {}
    for source, group in cols.groupby("source", sort=False):
        block = pd.read_csv(source, usecols=group["predictor"].tolist())
        block = _apply_row_mask(block, row_mask, source)
        for predictor in group["predictor"]:
            v = block[predictor].to_numpy(dtype=float)
            if return_top_indices:
                # Same computation predictor_tops would do on this column — taken here, while the
                # block is already in memory, rather than re-reading property_csvs for it.
                predictor_top_indices[predictor] = _top_indices(v, top_n)
            n_nan = int(np.isnan(v).sum())
            n_used, remap = n_total, None
            if n_nan:  # pairwise-complete, as in step 14
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
    if return_top_indices:
        return matrix, cols, predictor_top_indices
    return matrix, cols


def predictor_tops(cols, property_csvs, cutoffs=(ACTIVITY_BINARIZE_TOP_N,), row_mask=None):
    """``{top_n: {predictor: int64[top_n]}}`` row positions of each property predictor's top values.

    Takes a LIST of cutoffs and reads each property CSV exactly once, computing every cutoff from
    the column already in memory. The overlap view is drawn at three cutoffs, and the property CSVs
    total several GB — calling this once per cutoff would triple the step's I/O for no new data,
    since ``_top_indices`` is a cheap ``argpartition`` over an array that has already been read.

    Every predictor uses the plain highest-value :func:`_top_indices`: since 2026-09-02 the only
    property predictors are the merged cytotox/abx rank-sums (physchem, the one family that needed a
    two-sided "furthest from median" top-N, was dropped rather than merged), and higher rank-sum is
    unambiguously "more cytotoxic-like" / "more antibiotic-like" — a one-directional quantity like
    every other predictor here.

    ``row_mask`` must match the one given to :func:`organism_scores` — same reason as
    :func:`aggregated_matrix`.
    """
    out = {n: {} for n in cutoffs}
    for source, group in cols.groupby("source", sort=False):
        block = pd.read_csv(source, usecols=group["predictor"].tolist())
        block = _apply_row_mask(block, row_mask, source)
        for predictor in group["predictor"]:
            v = block[predictor].to_numpy(dtype=float)
            for n in cutoffs:
                out[n][predictor] = _top_indices(v, n)
        del block
    return out


def overlap_matrix(scores, organisms, cols, property_csvs,
                   top_n=ACTIVITY_BINARIZE_TOP_N, predictor_top_indices=None, row_mask=None):
    """How many of the ROW organism's ``top_n`` actives fall in the COLUMN's own ``top_n``.

    A different question from the AUROC matrix — that asks whether a predictor RANKS an organism's
    actives highly across the whole library, this asks how many of the very same molecules it puts at
    the top. A predictor can do the first well without doing the second.

    The raw intersection count, not Jaccard: both sets have exactly ``top_n`` members, so
    ``J = i / (2 * top_n - i)`` is a monotone re-expression of the same number and orders the matrix
    identically. The count is the one a reader can act on.

    The measure is SYMMETRIC, so the bioactivity block is a symmetric matrix (the AUROC block is not),
    and the diagonal is ``top_n`` by construction. The figure blanks it; these values do not.

    ``predictor_top_indices`` accepts the ``{predictor: positions}`` mapping for THIS cutoff, already
    computed by :func:`predictor_tops`. Pass it when drawing several cutoffs so the property CSVs are
    read once for all of them; omit it and this reads them itself.
    """
    order = organisms["organism"].tolist()
    n_total = len(scores)
    if predictor_top_indices is None:
        predictor_top_indices = predictor_tops(cols, property_csvs, (top_n,), row_mask)[top_n]
    tops = {o: set(_top_indices(scores[o].to_numpy(dtype=float), top_n)) for o in order}
    tops.update({p: set(v) for p, v in predictor_top_indices.items()})

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
