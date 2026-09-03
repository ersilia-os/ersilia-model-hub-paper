"""Per-ENDPOINT agreement statistics — the endpoint-level view of step 09's two pathogen-level
agreement summaries, plus a confounder check used by step 14.

Step 09 answers "do a pathogen's endpoints agree with each other more than with other pathogens'
endpoints?" TWICE — once via top-N Jaccard overlap, once via AUROC self-performance — and both of
its own figures roll the answer up to **one number per pathogen**. A single endpoint that behaves
badly is therefore invisible in those panels: it is averaged into its pathogen's box. This module
regroups exactly the same pair statistics **by endpoint** instead of by pathogen, so the question
becomes "can THIS endpoint uprank the compounds its own pathogen's other endpoints call active?".

Two statistics per endpoint, reported side by side, each against the same two comparators:

- **AUROC (upranking)** — the endpoint's raw score used to rank the library, scored against a peer
  endpoint's top-``ACTIVITY_BINARIZE_TOP_N`` compounds as the positive class. Read straight from
  step 09's ``09_activity_self_performance.csv``, which already holds the full directed 300x300
  grid. This is the "can it uprank" number: 0.5 is chance, and **below 0.5 means the endpoint ranks
  its peers' actives BELOW its own inactives**.
- **Top-N Jaccard (overlap)** — the overlap of the two endpoints' top-1000 compound sets, read from
  step 09's ``09_jaccard_top1000_baseline_matrix.csv``. Symmetric, so it has no direction.

A third, optional statistic — the **confounder check** — is added by step 14, after its own
property-predictor analysis has run: the strongest physchem/abx/cytotox predictor of each endpoint's
own top-1000, read from step 14's ``14_predictor_performance.csv``. It is not available at step-09
time (step 14 has not run yet), so :func:`endpoint_quality_table` treats it as optional.

Each is summarised as a ``same_median`` (over the endpoint's own pathogen's other endpoints), a
``diff_median`` (over every endpoint of the OTHER pathogens in scope), and their difference, which
follows step 09's naming and definition exactly::

    specificity = same_median - diff_median

Three properties of the scope are load-bearing and must reach any caption:

1.  **Restricting the node set removes a pathogen entirely**, per the contract of
    :func:`eval_correlations.pathogens_of_interest_nodes` and
    :func:`eval_correlations.multi_column_pathogen_nodes`. Applying the ">5 endpoints" rule on top
    of the 15 pathogens of interest therefore drops 3 pathogens as *partners* as well as as
    subjects, so ``diff_median`` here is against the other **11** pathogens — a narrower comparator
    than step 09's other-14, which is itself narrower than the all-57 comparator the pre-2026-08-07
    figures used.
2.  **``consensus_score`` columns are kept and analysed, but they are not fair comparators.** A
    consensus column is its model's aggregate over the very sub-endpoints it is being scored
    against, so its agreement with them is inflated by construction. It is flagged ``is_consensus``
    and is best read as that pathogen's within-block ceiling.
3.  **Same-model peers are INCLUDED** in every headline statistic — the literal "this endpoint
    against all the others of its pathogen". Where one model supplies most of a pathogen's
    endpoints this measures a model agreeing largely with itself, so ``n_same_model_peers`` is
    carried on every row as the caveat, and the directed pairs frame keeps its ``same_model`` flag
    so the excluded view can be recovered without recomputing anything.

Nothing here applies a threshold or flags an endpoint as bad. Every function returns a **ranking**;
which endpoints are genuinely off is a scientific judgement left to the reader. The one count that
is reported, ``n_below_chance_peers``, is a definition (AUROC < 0.5 is inverse ranking) rather than
a chosen cutoff.

Everything is a re-slice of CSVs steps 09 and 14 already wrote — no pass over the 1.5 GB score
matrix, so the whole module runs in seconds.
"""

import numpy as np
import pandas as pd

from default import CONSENSUS_COLUMN, PREDICTOR_CHANCE_LEVEL
from eval_correlations import (column_metric_pairs, multi_column_pathogen_nodes,
                               parse_named_column, pathogens_of_interest_nodes)


def _endpoint_key(model_id, column_name):
    """The ``"{model_id}:{column_name}"`` endpoint identifier used by every self-performance /
    predictor-performance CSV (steps 09 and 14).

    The bridge between the repo's two naming conventions: step 09's matrices are keyed
    ``"{pathogen}__{model_id}__{column_name}"`` (built by
    :func:`eval_correlations.build_named_score_matrix`), while the self-performance (step 09) and
    predictor-performance (step 14) long frames are keyed this way. Joining the two metrics is the
    only place both appear, so the conversion lives here rather than being spelled out at each call
    site.
    """
    return f"{model_id}:{column_name}"


def endpoint_nodes(metric_matrix, pathogens_of_interest_path, endpoint_selection_path,
                   min_endpoints, consensus_col=CONSENSUS_COLUMN):
    """The in-scope endpoint columns and their metadata.

    Two reductions, composed from the existing node selectors rather than reimplemented:

    1.  :func:`eval_correlations.pathogens_of_interest_nodes` — the 15 curated pathogens, matched on
        the organism name EXACTLY (it raises if one cannot be resolved, so an upstream rename fails
        loudly instead of quietly shrinking the analysis).
    2.  :func:`eval_correlations.multi_column_pathogen_nodes` — of those, the pathogens carrying at
        least ``min_endpoints`` columns.

    Both remove an excluded pathogen ENTIRELY, so the pathogens dropped here also stop being
    different-pathogen partners for the ones that remain. The dropped counts are returned so the
    caller can report them rather than have them vanish silently.

    Returns ``(nodes, meta, dropped)`` where ``nodes`` is the list of named columns, ``meta`` is a
    frame with one row per node (``named_column``, ``endpoint``, ``pathogen``, ``organism``,
    ``organism_class``, ``model_id``, ``column_name``, ``is_consensus``), and ``dropped`` is a
    Series of ``pathogen -> n_columns`` for the pathogens of interest that failed the
    ``min_endpoints`` rule.
    """
    poi_nodes = pathogens_of_interest_nodes(metric_matrix, pathogens_of_interest_path,
                                            endpoint_selection_path)
    nodes = multi_column_pathogen_nodes(metric_matrix.loc[poi_nodes, poi_nodes],
                                        min_columns=min_endpoints)

    poi_counts = pd.Series([parse_named_column(n)[0] for n in poi_nodes]).value_counts()
    kept = {parse_named_column(n)[0] for n in nodes}
    dropped = poi_counts[~poi_counts.index.isin(kept)]

    sel = pd.read_csv(endpoint_selection_path)
    sel = sel[sel["selected"] == "Yes"]
    org = sel.set_index(["model_id", "column_name"])[["organism", "organism_class"]]

    rows = []
    for name in nodes:
        pathogen, model_id, column_name = parse_named_column(name)
        organism, organism_class = org.loc[(model_id, column_name)]
        rows.append({
            "endpoint": _endpoint_key(model_id, column_name),
            "named_column": name,
            "pathogen": pathogen,
            "organism": organism,
            "organism_class": organism_class,
            "model_id": model_id,
            "column_name": column_name,
            "is_consensus": column_name == consensus_col,
        })
    return nodes, pd.DataFrame(rows), dropped


def jaccard_endpoint_pairs(jaccard_matrix, nodes):
    """Directed long-format Jaccard pairs over ``nodes``, keyed by endpoint identifier.

    Thin wrapper over :func:`eval_correlations.column_metric_pairs` — the same directed frame step
    09 aggregates by pathogen — with the ``"{model_id}:{column_name}"`` keys added so it can be
    joined to step 09's AUROC frame. Self-pairs are already dropped by ``column_metric_pairs``.
    """
    pairs = column_metric_pairs(jaccard_matrix.loc[nodes, nodes])
    pairs = pairs.rename(columns={"value": "jaccard"})
    for side, col in (("node", "endpoint"), ("partner", "peer")):
        parsed = pairs[side].map(parse_named_column)
        pairs[col] = [_endpoint_key(p[1], p[2]) for p in parsed]
    return pairs


def jaccard_endpoint_stats(pairs):
    """Per-endpoint same/different-pathogen Jaccard medians and their difference.

    ``specificity`` is ``same_median - diff_median``, the same definition and sign convention as
    :func:`eval_correlations.pathogen_metric_summary`. A negative value means the endpoint's top-1000
    overlaps other pathogens' endpoints MORE than its own pathogen's — reported, never filtered.

    ``n_same_model_peers`` counts how many of the endpoint's same-pathogen peers come from its own
    model. It is not used to filter anything; it is the caveat a reader needs in order to know
    whether a high ``same_median`` is cross-model agreement or one model agreeing with itself.
    """
    same = pairs[pairs["category"] == "same_pathogen"]
    diff = pairs[pairs["category"] == "different_pathogen"]
    out = pd.DataFrame({
        "jac_same_median": same.groupby("endpoint")["jaccard"].median(),
        "jac_diff_median": diff.groupby("endpoint")["jaccard"].median(),
        "n_peers": same.groupby("endpoint")["jaccard"].size(),
        "n_same_model_peers": same.groupby("endpoint")["same_model"].sum(),
    })
    out["jac_specificity"] = out["jac_same_median"] - out["jac_diff_median"]
    return out


def auroc_endpoint_pairs(self_perf, meta):
    """Directed AUROC pairs from step 09's ``09_activity_self_performance.csv``, restricted to scope.

    Filters both sides to the in-scope endpoints and drops the self-pairs, which are 1.0 by
    construction (an endpoint against its own binarization) and would otherwise put one guaranteed
    perfect value into every endpoint's same-pathogen distribution.

    ``same_pathogen`` is recomputed from ``meta`` rather than taken from the file's ``same_organism``
    column: that column was computed against all 57 organisms, and the comparator here is the
    in-scope pathogens only.
    """
    pathogen = meta.set_index("endpoint")["pathogen"]
    keep = set(meta["endpoint"])
    block = self_perf[self_perf["predictor_endpoint"].isin(keep)
                      & self_perf["target_endpoint"].isin(keep)
                      & ~self_perf["self_pair"]].copy()
    block = block.rename(columns={"predictor_endpoint": "endpoint",
                                  "target_endpoint": "peer", "value": "auroc"})
    block["pathogen"] = block["endpoint"].map(pathogen)
    block["peer_pathogen"] = block["peer"].map(pathogen)
    block["category"] = np.where(block["pathogen"] == block["peer_pathogen"],
                                 "same_pathogen", "different_pathogen")
    return block[["endpoint", "peer", "pathogen", "peer_pathogen", "category",
                  "same_model", "auroc"]]


def auroc_endpoint_stats(pairs, chance=PREDICTOR_CHANCE_LEVEL):
    """Per-endpoint AUROC agreement, in both directions.

    AUROC is **not** symmetric — an endpoint can be well predicted by its peers without being able
    to predict them — so both directions are summarised and they answer different questions:

    - **outgoing** (``auroc_out_*``): this endpoint's raw score ranking its peers' actives. The
      primary "can it uprank" statistic.
    - **incoming** (``auroc_in_*``): its peers' scores ranking THIS endpoint's actives, i.e. how
      predictable its own top-1000 is. Secondary and diagnostic: an endpoint that is neither
      predictive nor predictable is idiosyncratic, while one that is predictable but not predictive
      is more consistent with a noisy or degenerate score.

    ``n_below_chance_peers`` counts the same-pathogen peers this endpoint ranks below chance. It is
    a definition, not a threshold: AUROC < 0.5 means the peer's actives were ranked below its
    inactives, which is a statement about direction rather than a chosen cutoff.
    """
    same = pairs[pairs["category"] == "same_pathogen"]
    diff = pairs[pairs["category"] == "different_pathogen"]
    out = pd.DataFrame({
        "auroc_out_same_median": same.groupby("endpoint")["auroc"].median(),
        "auroc_out_diff_median": diff.groupby("endpoint")["auroc"].median(),
        "auroc_out_min": same.groupby("endpoint")["auroc"].min(),
        "n_below_chance_peers": same.assign(low=same["auroc"] < chance)
                                    .groupby("endpoint")["low"].sum(),
        "auroc_in_same_median": same.groupby("peer")["auroc"].median(),
        "auroc_in_diff_median": diff.groupby("peer")["auroc"].median(),
    })
    out["auroc_out_specificity"] = out["auroc_out_same_median"] - out["auroc_out_diff_median"]
    return out


def confounder_stats(predictor_perf, meta, chance=PREDICTOR_CHANCE_LEVEL):
    """The strongest physchem / cytotox / abx predictor of each endpoint's own top-1000.

    Read from step 14's ``14_predictor_performance.csv`` (101 property predictors x 300 activity
    targets). The question it answers is different from the peer-agreement one: an endpoint whose
    top-1000 is recovered almost perfectly by a single physicochemical descriptor is, whatever its
    peer agreement, largely reproducing a property filter rather than a bioactivity ranking.

    Strength is ``|value - chance|``, so a strong INVERSE predictor counts as strongly as a direct
    one — ``qed`` and ``clogp`` predict several endpoints well below 0.5, and calling those weak
    would invert the finding. The signed ``confounder_value`` is kept alongside so the direction
    stays readable.

    ``same_model`` rows are excluded: ``eos3dys`` appears in both ``config/cytotoxicity_models.csv``
    and the endpoint selection, so leaving them in would let one of its columns "predict" another of
    its own columns and be reported as an external confounder.
    """
    block = predictor_perf[predictor_perf["target_endpoint"].isin(set(meta["endpoint"]))
                           & ~predictor_perf["same_model"]].copy()
    block["abs_dev"] = (block["value"] - chance).abs()
    block = block.dropna(subset=["abs_dev"])
    if not len(block):
        return pd.DataFrame(columns=["confounder_predictor", "confounder_family",
                                     "confounder_metric", "confounder_value",
                                     "confounder_abs_dev"])
    best = block.loc[block.groupby("target_endpoint")["abs_dev"].idxmax()]
    return (best.set_index("target_endpoint")[["predictor", "family", "metric", "value", "abs_dev"]]
                .rename(columns={"predictor": "confounder_predictor",
                                 "family": "confounder_family",
                                 "metric": "confounder_metric",
                                 "value": "confounder_value",
                                 "abs_dev": "confounder_abs_dev"})
                .rename_axis("endpoint"))


def endpoint_quality_table(meta, jac_stats, auroc_stats, conf_stats=None):
    """The master per-endpoint frame: metadata + both agreement metrics + the confounder check.

    ``conf_stats`` is optional (default ``None``): step 09 builds this table before step 14's
    property-predictor analysis has run, so it calls this with the Jaccard/AUROC statistics only.
    Step 14 re-derives the same table with ``conf_stats`` supplied, adding the confounder columns.

    Sorted by ``auroc_out_same_median`` ASCENDING, so the endpoints least able to uprank their own
    pathogen's peers are at the top. ``rank_overall`` and ``rank_within_pathogen`` follow the same
    order (rank 1 = weakest), the second one because a pathogen whose whole block is weak — MTB is
    the clearest case — would otherwise bury its own relative outliers at the bottom of a global
    ranking.

    No row is ever dropped and no flag column is added: the table is a ranking, and deciding which
    endpoints are genuinely off is left to the reader.
    """
    parts = [jac_stats, auroc_stats] + ([conf_stats] if conf_stats is not None else [])
    out = meta.set_index("endpoint").join(parts).reset_index()
    out = out.sort_values("auroc_out_same_median", ascending=True).reset_index(drop=True)
    out["rank_overall"] = out["auroc_out_same_median"].rank(method="min").astype("Int64")
    out["rank_within_pathogen"] = (out.groupby("pathogen")["auroc_out_same_median"]
                                      .rank(method="min").astype("Int64"))
    return out


def pathogen_endpoint_summary(table, chance=PREDICTOR_CHANCE_LEVEL):
    """Per-pathogen roll-up of the endpoint table, weakest pathogen first.

    ``n_below_chance`` counts endpoints whose median AUROC against their own peers is below chance —
    the pathogen-level tally of the same definition used per endpoint. ``n_single_model`` reports
    how many of the pathogen's endpoints have no cross-model peer at all, which is the condition
    under which its within-pathogen agreement is one model agreeing with itself.
    """
    rows = []
    for pathogen, g in table.groupby("pathogen", sort=False):
        rows.append({
            "pathogen": pathogen,
            "organism": g["organism"].iloc[0],
            "organism_class": g["organism_class"].iloc[0],
            "n_endpoints": len(g),
            "n_models": g["model_id"].nunique(),
            "n_consensus": int(g["is_consensus"].sum()),
            "auroc_out_same_median": g["auroc_out_same_median"].median(),
            "auroc_out_specificity_median": g["auroc_out_specificity"].median(),
            "jac_same_median": g["jac_same_median"].median(),
            "jac_specificity_median": g["jac_specificity"].median(),
            "n_below_chance": int((g["auroc_out_same_median"] < chance).sum()),
            "n_single_model": int((g["n_same_model_peers"] == g["n_peers"]).sum()),
        })
    return (pd.DataFrame(rows).sort_values("auroc_out_same_median", ascending=True)
                              .reset_index(drop=True))
