"""Top-N Jaccard aggregated by ORGANISM CLASS rather than by pathogen — step 09 one level up.

Step 08 asks whether a *pathogen's* endpoints pick out the same compounds as each other more than as
other pathogens' endpoints. This module asks the same question of an organism CLASS: do all the
Gram-negative endpoints agree with each other more than with the Gram-positives, the fungi and the
rest? That is the level at which cross-organism transfer would show up — a compound series active
across the Gram-negatives is a different claim from one active across everything.

It is a pure re-aggregation of `09_jaccard_top1000_baseline_matrix.csv`, which step 09 already wrote
over all 57 organisms. Nothing is recomputed.

**Three boxes per class, not two.** The class-level ``same`` box is dominated by pairs that are
already same-*pathogen*, which step 09 published — so the box that carries the new information is
``same_excl_same_organism``: same class, different pathogen. It is computed through the
``exclusions`` hook of :func:`eval_correlations.pathogen_metric_boxes`, alongside the
``same_excl_same_model`` box that step 09 has always carried.

**The scope makes four of the six classes degenerate, and that is reported rather than hidden.**
Restricted to the 15 pathogens of interest (the chosen scoping, matching step 09), Mycobacteria is
*M. tuberculosis* alone, Helminths is *S. mansoni* alone, Fungi is *C. albicans* alone and Protozoa
is *P. falciparum* alone. For those four, ``same`` is a verbatim copy of that pathogen's step-09
same-pathogen box and ``same_excl_same_organism`` is EMPTY. Only Gram-negative (8 organisms) and
Gram-positive (3) carry information the per-pathogen figure does not already show. Every summary row
therefore carries ``n_organisms``, the figure draws no third box where there is nothing to draw, and
the driver asserts the four copies against step 09's own summary — turning the degeneracy into a
correctness check instead of a silent caveat.
"""

import numpy as np
import pandas as pd

from eval_correlations import (column_metric_pairs, pathogen_metric_boxes,
                               pathogen_metric_summary, parse_named_column)

#: Both exclusion boxes, in the order they should read on the figure and in the summary.
#: ``same_excl_same_organism`` is the one this analysis exists for; ``same_excl_same_model`` is
#: carried for continuity with step 09, where it guards against a class's agreement being one model
#: agreeing with itself — which matters more here, not less, because a class pools several models.
CLASS_EXCLUSIONS = {
    "same_excl_same_organism": "same_organism",
    "same_excl_same_model": "same_model",
}


def node_organism_classes(nodes, endpoint_selection_path):
    """``{named column: organism_class}`` for every node, joined from the endpoint selection.

    The organism class is **not** encoded in the column name — ``build_named_score_matrix`` writes
    ``{pathogen_code}__{model_id}__{column_name}`` and nothing more — so a join is unavoidable. It
    goes through ``(model_id, column_name)`` rather than through the pathogen code, because that
    pair is the selection file's own key and cannot be affected by the mechanical code fallback in
    :func:`eval_correlations._pathogen_code`.

    Raises if any node has no row in the selection, rather than silently producing a NaN class that
    would then form its own group.
    """
    sel = pd.read_csv(endpoint_selection_path)
    sel = sel[sel["selected"] == "Yes"]
    classes = sel.set_index(["model_id", "column_name"])["organism_class"]

    out, missing = {}, []
    for name in nodes:
        _, model_id, column_name = parse_named_column(name)
        try:
            out[name] = classes.loc[(model_id, column_name)]
        except KeyError:
            missing.append(name)
    if missing:
        raise ValueError(
            f"{len(missing)} node(s) have no selected row in {endpoint_selection_path}: "
            f"{missing[:5]}{' ...' if len(missing) > 5 else ''}")
    if pd.isna(pd.Series(out)).any():
        blank = sorted(n for n, c in out.items() if pd.isna(c))
        raise ValueError(f"{len(blank)} node(s) have a blank organism_class: {blank[:5]}")
    return out


def class_metric_pairs(metric_matrix, nodes, endpoint_selection_path):
    """Directed pair frame regrouped by organism class.

    Reuses :func:`eval_correlations.column_metric_pairs` verbatim, then re-labels it: the original
    per-pathogen codes are preserved as ``organism`` / ``partner_organism`` with a ``same_organism``
    flag, and ``pathogen`` / ``category`` are overwritten with the class-level values. The
    downstream aggregators key off ``pathogen`` and ``category`` and never read a column name, so
    they group by class unchanged.

    Keeping the pathogen codes is what makes ``same_excl_same_organism`` possible, and it is also
    what lets the pairs CSV be re-sliced later without recomputing the Jaccard matrix.
    """
    node_class = node_organism_classes(nodes, endpoint_selection_path)

    pairs = column_metric_pairs(metric_matrix.loc[nodes, nodes])
    pairs = pairs.rename(columns={"pathogen": "organism", "partner_pathogen": "partner_organism"})
    pairs["same_organism"] = pairs["organism"] == pairs["partner_organism"]
    pairs["pathogen"] = pairs["node"].map(node_class)
    pairs["partner_pathogen"] = pairs["partner"].map(node_class)
    pairs["category"] = np.where(pairs["pathogen"] == pairs["partner_pathogen"],
                                 "same_pathogen", "different_pathogen")
    return pairs, pd.Series(node_class, name="organism_class")


def class_metric_boxes(pairs):
    """Per-class ``same`` / ``same_excl_same_organism`` / ``same_excl_same_model`` / ``diff`` boxes."""
    return pathogen_metric_boxes(pairs, exclusions=CLASS_EXCLUSIONS)


def class_metric_summary(boxes, n_columns, n_organisms):
    """Tidy per-class summary, ordered by specificity, with the organism count carried on every row.

    ``n_organisms`` is not decoration: a class of one organism has an empty
    ``same_excl_same_organism`` box and a ``same_median`` that merely restates step 09, so the count
    is what tells a reader which rows are informative. It sits immediately after ``n_columns`` so the
    two are read together.
    """
    out = pathogen_metric_summary(boxes, n_columns, exclusions=CLASS_EXCLUSIONS)
    out = out.rename(columns={"pathogen": "organism_class"})
    out.insert(2, "n_organisms", out["organism_class"].map(n_organisms).astype(int))
    return out
