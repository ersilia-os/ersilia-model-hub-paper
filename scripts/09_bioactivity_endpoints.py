"""Step 09 — bioactivity-only endpoint agreement: per-pathogen top-1000 Jaccard, the same one level
up by organism class (absorbed from the former standalone ``xx_group_jaccard.py``, 2026-09-02), then
the per-endpoint audit (Jaccard specificity, AUROC self-performance), absorbed from the former
standalone step 15 (2026-09-02). Depends on nothing but step 07's cached bioactivity matrix — no
physchem/abx/cytotox property data enters this step at all, which is exactly what lets the AUROC
self-performance analysis live here rather than in step 14 alongside the property-predictor analysis.

Part 1 asks whether endpoints of the SAME pathogen pick out the same compounds more than endpoints of
DIFFERENT pathogens do. For every pair of the selected endpoint columns, the Jaccard overlap of their
top-1000 highest-scoring compounds (out of 1,355,109), then aggregated per pathogen: turquoise =
every unordered pair of that pathogen's own columns, crimson = its columns against every other
pathogen's.

**Scope: the 15 curated pathogens of interest** (``config/pathogens_of_interest.csv``), replacing the
former ``min<K>``-endpoint thresholds (2026-08-07). Restricting the node set removes every other
pathogen ENTIRELY — the ~40 gut-microbiome and other organisms stop being different-pathogen partners
too — so each pathogen's crimson box is against the other 14 of interest only, **not** against all 57.
That is a deliberately narrower comparator than the old ``min5`` figure used, and a caption must say
so: "specific to this pathogen" here means "relative to the other priority pathogens".

**Two of the 15 have a single endpoint** (*Campylobacter*, *H. pylori*), so they have no
same-pathogen pair and no turquoise box — only a crimson one, and a NaN ``same_median``. They are
kept rather than dropped: for a pathogen on the priority list, having too little in the hub to assess
is itself the result, and hiding the row would hide it.

**Only the unscaled (baseline) matrix is computed and published — trimmed from 3 variants to 1 on
2026-09-01.** The other two (z-score + L2 row-norm, rank-percentile + L1 row-norm) used to get their
own full Jaccard computation, summary CSV and diagnostic figure, but nothing downstream ever read
either one — only ``baseline`` (this step's own per-organism-class aggregation, part 2 below). They
were a genuine robustness check (does row-normalizing change which pathogen looks "specific"?) rather
than dead weight, but at ~2x the compute of baseline alone for 2 figures nobody consumed, they were
cut on review; re-add them from ``git log`` if that robustness question needs revisiting.

**A cheaper, different check is still run and kept: does *column* scaling alone ever change a top-N
Jaccard set?** Top-N Jaccard depends only on each column's own internal ranking, and both of step 07's
column scalings (z-score, rank-percentile) are strictly increasing per column — so in exact arithmetic
neither could change any column's top-1000 set, and the baseline matrix should equal what either
scaling would give. The identity is **asserted at runtime**, not assumed, over the SAME cached
``base`` matrix this step already holds (no extra I/O) — and it earns its keep: baseline ==
rank-percentiled holds exactly, but baseline == z-scored comes back FALSE, because ``(x - mean) / std``
in float32 reorders near-tied values. That per-cell effect was quantified at 300 endpoints (one column,
~0.001 on 156 of 90,000 Jaccard cells); the script prints only the True/False check, so it was NOT
re-quantified on the 307-endpoint rebuild, where the check still returns the same verdicts. This check
is unrelated to the row-normalized variants dropped above — it never row-normalizes anything.

The full 307x307 Jaccard matrices are written out for reuse by later analysis — they are computed over
ALL pathogens, and only the per-pathogen aggregation is restricted to the 15.

Part 2 (absorbed from the former standalone ``xx_group_jaccard.py``, 2026-09-02) is part 1 one level
up: the SAME baseline Jaccard matrix, aggregated by ORGANISM CLASS instead of by pathogen. Do all the
Gram-negative endpoints agree with each other more than with the Gram-positives, the fungi and the
rest? That is the level at which cross-organism transfer would show up. Nothing is recomputed — it
re-slices part 1's own in-memory matrix and summary — and it runs in seconds. **Four of the six
classes are a single organism** under the 15-pathogen scoping (Mycobacteria = *M. tuberculosis*,
Helminths = *S. mansoni*, Fungi = *C. albicans*, Protozoa = *P. falciparum*), so for those four the
same-class box is a verbatim copy of that pathogen's part-1 box — asserted against part 1's own
summary as a correctness check, not just a caveat. Only Gram-negative (8 organisms) and Gram-positive
(3) carry information part 1 does not already show. See ``src/eval_group_jaccard.py`` for the full
three-box (``same`` / ``same_excl_same_organism`` / ``diff``) contract.

Part 3 (absorbed from the former standalone step 15) regroups the SAME per-pathogen statistics **by
endpoint** instead: can THIS endpoint uprank the compounds its own pathogen's other endpoints call
active? A single endpoint that behaves badly is invisible in part 1's pathogen-level boxes, averaged
into its pathogen's box. Two statistics per endpoint: AUROC upranking (computed fresh here, over the
step-07 matrix only, via :func:`eval_predictor_performance.run_activity_self_performance`) and
top-1000 Jaccard overlap (reused from part 1's own baseline matrix). Scope: pathogens of interest with
**>5 endpoints** (``MIN_ENDPOINTS``) — narrower than part 1's 15, since a pathogen with 5 or fewer
endpoints cannot support a per-endpoint peer distribution. See ``src/eval_endpoint_quality.py`` for
the full statistics contract (consensus columns, same-model peers, no threshold applied). The
confounder check that used to run alongside this in step 15 (is an endpoint's top-1000 really
tracking a physchem/abx/cytotox property?) now runs in step 14, once its property-predictor data
exists.

Reporting choices worth knowing:

  - **Same-model pairs are INCLUDED** in the boxes — the literal "each column against all others". Two
    output columns of one model agreeing says little about pathogen specificity, so the summary CSV
    carries ``same_median_excl_same_model`` alongside ``same_median``. Read them together: several
    pathogens have NO cross-model same-pathogen pair at all, so their turquoise box is one model
    agreeing with itself.
  - **Linear x-axis** — values bunch near zero, but exact-zero pairs render instead of being silently
    dropped by a log axis. Nothing is filtered.

  This is a DIAGNOSTIC, not a paper panel: plain matplotlib rather than the 3 cm cell grid, because
  the y axis carries per-pathogen endpoint and pair counts that go illegible at page width.

    python 09_bioactivity_endpoints.py

Outputs
-------
    output/09_bioactivity_endpoints/09_jaccard_top1000_baseline_matrix.csv          (307x307, reused)
    output/09_bioactivity_endpoints/09_pathogen_jaccard_top1000_baseline_summary.csv
    output/09_bioactivity_endpoints/png/09_pathogen_jaccard_top1000_baseline.png
    output/09_bioactivity_endpoints/pdf/09_pathogen_jaccard_top1000_baseline.pdf
    output/09_bioactivity_endpoints/09_group_jaccard_top1000_summary.csv        (6 organism classes)
    output/09_bioactivity_endpoints/09_group_jaccard_pairs.csv                  (directed pairs)
    output/09_bioactivity_endpoints/png/09_group_jaccard_top1000.png
    output/09_bioactivity_endpoints/pdf/09_group_jaccard_top1000.pdf
    output/09_bioactivity_endpoints/09_activity_self_performance.csv            (all selected targets)
    output/09_bioactivity_endpoints/09_pathogen_subset_self_performance.csv     (15 pathogens, consensus collapsed)
    output/09_bioactivity_endpoints/09_endpoint_quality.csv                     (one row per endpoint)
    output/09_bioactivity_endpoints/09_endpoint_pairs.csv                       (directed, Jaccard + AUROC)
    output/09_bioactivity_endpoints/09_pathogen_endpoint_summary.csv
    output/09_bioactivity_endpoints/png|pdf/09_endpoint_*.{png,pdf}
    output/09_bioactivity_endpoints/png|pdf/09_performance_activity_by_organism.{png,pdf}
    output/09_bioactivity_endpoints/png|pdf/09_performance_pathogen_subset.{png,pdf}
    output/09_bioactivity_endpoints/figure_cells.json
"""

import os
import sys
import time

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))

from default import (  # noqa: E402
    ACTIVITY_BINARIZE_TOP_N, ANNOTATION_PREDS_SUBDIR, ORGANISM_CLASS_ORDER, PREDICTOR_CHANCE_LEVEL,
)
from eval_correlations import (  # noqa: E402
    build_named_score_matrix,
    column_metric_pairs,
    pathogens_of_interest_nodes,
    parse_named_column,
    pathogen_metric_boxes,
    pathogen_metric_summary,
    scale_matrix,
    topn_jaccard_matrix,
)
from eval_endpoint_quality import (  # noqa: E402
    auroc_endpoint_pairs,
    auroc_endpoint_stats,
    endpoint_nodes,
    endpoint_quality_table,
    jaccard_endpoint_pairs,
    jaccard_endpoint_stats,
    pathogen_endpoint_summary,
)
from eval_group_jaccard import class_metric_boxes, class_metric_pairs, class_metric_summary  # noqa: E402
from eval_predictor_performance import run_activity_self_performance  # noqa: E402
from plots_endpoint_quality import (  # noqa: E402
    save_activity_self_figure, save_endpoint_quality_figures, save_pathogen_subset_figure,
)
from plots_matrix_analyses import group_jaccard_figure, pathogen_jaccard_figure  # noqa: E402

#: ">5 endpoints" — a pathogen needs at least this many for each of its endpoints to be judged
#: against a peer distribution rather than against one or two other columns.
MIN_ENDPOINTS = 6

CUTOFF = 1000

pred_dir = os.path.join(root, "..", "data", "processed", ANNOTATION_PREDS_SUBDIR)
matrix_dir = os.path.join(root, "..", "output", "07_score_matrices")
output_dir = os.path.join(root, "..", "output", "09_bioactivity_endpoints")
config_dir = os.path.join(root, "..", "config")
os.makedirs(output_dir, exist_ok=True)

endpoint_selection_path = os.path.join(config_dir, "08_endpoint_selection.csv")
pathogens_of_interest_path = os.path.join(config_dir, "pathogens_of_interest.csv")
full_matrix_cache_path = os.path.join(matrix_dir, "07_score_matrix_full.parquet")

# (slug, human label, how to derive the matrix from the base named matrix). A single-entry list
# rather than unrolled code: the z-score+L2rownorm / rank-percentile+L1rownorm variants that used to
# sit alongside "baseline" here were dropped 2026-09-01 (see the docstring), and this shape is kept
# so a future robustness variant can be added back without restructuring the two loops below.
VARIANTS = [
    ("baseline", "unscaled scores (= rank-percentiled; = z-scored bar one column)", lambda m: m),
]

# ----------------------------------------------------------------------------- #
# 1. The 307x307 baseline Jaccard matrix (cached)                                #
# ----------------------------------------------------------------------------- #
jaccards = {}
missing = [v for v in VARIANTS
           if not os.path.exists(os.path.join(output_dir, f"09_jaccard_top{CUTOFF}_{v[0]}_matrix.csv"))]
base = None
if missing:
    t0 = time.time()
    base = build_named_score_matrix(
        pred_dir=pred_dir, endpoint_selection_path=endpoint_selection_path,
        pathogens_of_interest_path=pathogens_of_interest_path,
        full_matrix_cache_path=full_matrix_cache_path)
    print(f"[pathogen-jaccard] base named matrix {base.shape} ready in {time.time() - t0:.1f}s")

for slug, label, derive in VARIANTS:
    jac_path = os.path.join(output_dir, f"09_jaccard_top{CUTOFF}_{slug}_matrix.csv")
    if os.path.exists(jac_path):
        jaccards[slug] = pd.read_csv(jac_path, index_col=0)
        print(f"[pathogen-jaccard] {slug}: reusing cached {os.path.basename(jac_path)}")
        continue
    t0 = time.time()
    variant = derive(base)
    jac = topn_jaccard_matrix(variant, CUTOFF)
    del variant
    jac.to_csv(jac_path)
    jaccards[slug] = jac
    print(f"[pathogen-jaccard] {slug}: top-{CUTOFF} Jaccard {jac.shape} in "
          f"{time.time() - t0:.1f}s -> {os.path.basename(jac_path)}")

# The scaling-invariance claim is asserted, not assumed: if a future change to scale_matrix broke
# monotonicity, the "baseline covers all three" label would silently become a lie.
if base is not None:
    for name, method in (("z-scored", "zscore"), ("rank-percentiled", "rank_pct")):
        other = topn_jaccard_matrix(scale_matrix(base, method), CUTOFF)
        same = np.allclose(jaccards["baseline"].to_numpy(), other.to_numpy(), equal_nan=True)
        print(f"[pathogen-jaccard] check: baseline == {name} top-{CUTOFF} Jaccard -> {same}")
    del base

# ----------------------------------------------------------------------------- #
# 2. Per-pathogen aggregation and figure, per threshold
# ----------------------------------------------------------------------------- #
for slug, label, _ in VARIANTS:
    jac = jaccards[slug]
    nodes = pathogens_of_interest_nodes(jac, pathogens_of_interest_path, endpoint_selection_path)
    all_pathogens = pd.Series([parse_named_column(c)[0] for c in jac.columns]).value_counts()
    n_columns = pd.Series([parse_named_column(n)[0] for n in nodes]).value_counts()
    dropped = all_pathogens.drop(index=n_columns.index)

    pairs = column_metric_pairs(jac.loc[nodes, nodes])
    boxes = pathogen_metric_boxes(pairs)
    summary = pathogen_metric_summary(boxes, n_columns)
    stem = f"09_pathogen_jaccard_top{CUTOFF}_{slug}"
    summary.to_csv(os.path.join(output_dir, f"{stem}_summary.csv"), index=False)

    png_path, pdf_path = pathogen_jaccard_figure(
        boxes, summary, cutoff=CUTOFF,
        matrix_label=f"{label}  (15 pathogens of interest)",
        name=stem, output_dir=output_dir)

    print(f"[pathogen-jaccard] {slug}")
    print(f"    kept {len(nodes)} of {len(jac.columns)} columns, "
          f"{len(summary)} of {len(all_pathogens)} pathogens (the curated 15)")
    print(f"    excluded {len(dropped)} other pathogens / {int(dropped.sum())} columns — they are "
          "no longer different-pathogen partners either")
    # A priority pathogen with one endpoint has no same-pathogen pair. It is KEPT, with a diff box
    # and a NaN same_median, because the gap is the finding for a pathogen on the priority list.
    no_same = summary[summary["n_same_pairs"] == 0]
    if len(no_same):
        print(f"    {len(no_same)} pathogen(s) with a single endpoint, so no same-pathogen pair "
              f"and no same box: {', '.join(no_same['pathogen'])}")
    print(f"    wrote {os.path.basename(png_path)} + {os.path.basename(pdf_path)}")
    print(summary.to_string(index=False, float_format=lambda v: f"{v:.4f}"))
    print()

# ----------------------------------------------------------------------------- #
# 3. Per-organism-class aggregation and figure (absorbed from the former         #
#    standalone xx_group_jaccard.py, 2026-09-02)                                 #
# ----------------------------------------------------------------------------- #
# Re-slices `jac`/`nodes`/`summary`, still in memory from section 2's single-entry VARIANTS loop —
# no new I/O, no new Jaccard computation, just a different grouping of the same 307x307 matrix.
class_pairs, node_class = class_metric_pairs(jac, nodes, endpoint_selection_path)

n_columns_by_class = node_class.value_counts()
organisms = pd.Series([parse_named_column(n)[0] for n in nodes], index=nodes)
n_organisms = pd.DataFrame({"organism_class": node_class, "organism": organisms}) \
                .groupby("organism_class")["organism"].nunique()

print(f"[group-jaccard] {len(nodes)} endpoints of the 15 pathogens of interest -> "
      f"{n_columns_by_class.size} organism classes")
for cls in ORGANISM_CLASS_ORDER:
    if cls in n_columns_by_class.index:
        note = "  <- single organism: same-class == same-pathogen" if n_organisms[cls] == 1 else ""
        print(f"    {cls:<26} {int(n_columns_by_class[cls]):>3} endpoints, "
              f"{int(n_organisms[cls]):>2} organism(s){note}")

# The class assignment is the one thing that could be wrong without raising, so it is checked
# against the config rather than against a hardcoded list.
expected_class = (pd.read_csv(endpoint_selection_path).query("selected == 'Yes'")
              .merge(pd.read_csv(pathogens_of_interest_path), left_on="organism",
                     right_on="pathogen")["organism_class"].value_counts())
assert n_columns_by_class.sort_index().equals(expected_class.sort_index()), (
    f"class endpoint counts disagree with {endpoint_selection_path}:\n"
    f"{pd.DataFrame({'nodes': n_columns_by_class, 'config': expected_class}).to_string()}")
unknown = set(n_columns_by_class.index) - set(ORGANISM_CLASS_ORDER)
assert not unknown, f"organism class(es) absent from ORGANISM_CLASS_ORDER: {sorted(unknown)}"

class_boxes = class_metric_boxes(class_pairs)
class_summary = class_metric_summary(class_boxes, n_columns_by_class, n_organisms)
class_summary.to_csv(os.path.join(output_dir, f"09_group_jaccard_top{CUTOFF}_summary.csv"),
                     index=False)
class_pairs.to_csv(os.path.join(output_dir, "09_group_jaccard_pairs.csv"), index=False)

# The degeneracy, as a correctness check: for a class of one organism, "same class" and "same
# pathogen" are the same set of pairs, so the class median MUST equal that pathogen's part-1 median
# exactly and the different-pathogen box must be empty. Checked against `summary` (part 1's own,
# still in memory above), not re-read from disk. If the class join were wrong, this is where it
# would show.
step09_by_pathogen = summary.set_index("pathogen")
single = [c for c in class_summary["organism_class"] if n_organisms[c] == 1]
code_of = {node_class[n]: parse_named_column(n)[0] for n in nodes}
for cls in single:
    row = class_summary.set_index("organism_class").loc[cls]
    code = code_of[cls]
    ours, theirs = float(row["same_median"]), float(step09_by_pathogen.loc[code, "same_median"])
    assert abs(ours - theirs) < 1e-12, (
        f"{cls} is a single organism ({code}) so its same-class median must equal part 1's "
        f"same-pathogen median: {ours} vs {theirs}")
    assert not len(class_boxes[cls]["same_excl_same_organism"]), (
        f"{cls} has one organism but a non-empty same-class different-pathogen box")
print(f"[group-jaccard] {len(single)} single-organism class(es) reproduce part 1 exactly "
      f"({', '.join(single)}) and have an empty different-pathogen box")

informative = [c for c in class_summary["organism_class"] if n_organisms[c] > 1]
print(f"[group-jaccard] {len(informative)} class(es) carry new information: "
      f"{', '.join(informative)}")

print("\n[group-jaccard] per-class summary (by specificity, descending):")
print(class_summary.to_string(index=False, float_format=lambda v: f"{v:.4f}"))

for cls in informative:
    b = class_boxes[cls]
    same, excl = np.median(b["same"]), np.median(b["same_excl_same_organism"])
    print(f"\n[group-jaccard] {cls}: same-class median {same:.4f} -> {excl:.4f} once "
          f"within-pathogen pairs are removed ({len(b['same'])} -> "
          f"{len(b['same_excl_same_organism'])} pairs), against a different-class median of "
          f"{np.median(b['diff']):.4f}")

png_path, pdf_path = group_jaccard_figure(
    class_boxes, class_summary, cutoff=CUTOFF,
    matrix_label="unscaled scores  (15 pathogens of interest, by organism class)",
    name=f"09_group_jaccard_top{CUTOFF}", output_dir=output_dir)
print(f"\n  figure: {os.path.basename(png_path)[:-4]} (diagnostic, no cell footprint)")
print()

assert CUTOFF == ACTIVITY_BINARIZE_TOP_N, (
    f"Jaccard cutoff {CUTOFF} != ACTIVITY_BINARIZE_TOP_N {ACTIVITY_BINARIZE_TOP_N}: the two "
    "metrics would be computed over different compound sets and could not be compared.")

# ----------------------------------------------------------------------------- #
# 4. Per-endpoint agreement audit (former step 15)                               #
# ----------------------------------------------------------------------------- #
# 4a. AUROC self-performance — bioactivity-only, computed fresh from the step-07 matrix.
targets, tops, self_perf, subset = run_activity_self_performance(
    parquet_path=full_matrix_cache_path, selection_csv=endpoint_selection_path,
    pathogens_csv=pathogens_of_interest_path, output_dir=output_dir)
save_activity_self_figure(output_dir, subset)
save_pathogen_subset_figure(output_dir, subset)

# 4b. Scope: pathogens of interest with > MIN_ENDPOINTS - 1 endpoints, narrower than part 1's 15.
jaccard = jaccards["baseline"]
nodes, meta, dropped_endpoints = endpoint_nodes(jaccard, pathogens_of_interest_path,
                                                endpoint_selection_path, MIN_ENDPOINTS)
counts = meta["pathogen"].value_counts()
print(f"\n[endpoint-quality] {len(nodes)} endpoints across {counts.size} pathogens "
      f"(>{MIN_ENDPOINTS - 1} endpoints each), {int(meta['is_consensus'].sum())} of them consensus "
      "columns")
if len(dropped_endpoints):
    print(f"    dropped {len(dropped_endpoints)} pathogen(s) of interest below the endpoint "
          "minimum — they are no longer cross-pathogen partners either: "
          + ", ".join(f"{p} ({n})" for p, n in dropped_endpoints.items()))

# The whole analysis is a per-pathogen peer comparison, so a wrong node set would silently change
# every statistic rather than raise. Checked against the config, not against a hardcoded list.
expected = (pd.read_csv(endpoint_selection_path)
              .query("selected == 'Yes'")
              .merge(pd.read_csv(pathogens_of_interest_path), left_on="organism",
                     right_on="pathogen")["code"].value_counts())
expected = expected[expected >= MIN_ENDPOINTS]
assert counts.sort_index().equals(expected.sort_index()), (
    f"endpoint counts disagree with {endpoint_selection_path}:\n"
    f"{pd.DataFrame({'nodes': counts, 'config': expected}).to_string()}")

# 4c. The two metrics, per endpoint. No confounder check here — that needs step 14's property-
# predictor data, which does not exist yet; it is added by step 14 once it does.
jac_pairs = jaccard_endpoint_pairs(jaccard, nodes)
auroc_pairs = auroc_endpoint_pairs(self_perf, meta)
jac_stats = jaccard_endpoint_stats(jac_pairs)
auroc_stats = auroc_endpoint_stats(auroc_pairs)
print(f"[endpoint-quality] {len(jac_pairs):,} Jaccard pairs, {len(auroc_pairs):,} AUROC pairs")

table = endpoint_quality_table(meta, jac_stats, auroc_stats)
endpoint_summary = pathogen_endpoint_summary(table)

# One directed row per (endpoint, peer) carrying both metrics, so the same-model-excluded view and
# any per-pair follow-up can be recovered without recomputing either matrix.
endpoint_pairs = auroc_pairs.merge(jac_pairs[["endpoint", "peer", "jaccard"]],
                                   on=["endpoint", "peer"], how="outer")
assert len(endpoint_pairs) == len(auroc_pairs) == len(jac_pairs), (
    f"pair frames do not align: {len(jac_pairs)} Jaccard, {len(auroc_pairs)} AUROC, "
    f"{len(endpoint_pairs)} merged — the two naming conventions disagree on some endpoint.")

table.to_csv(os.path.join(output_dir, "09_endpoint_quality.csv"), index=False)
endpoint_pairs.to_csv(os.path.join(output_dir, "09_endpoint_pairs.csv"), index=False)
endpoint_summary.to_csv(os.path.join(output_dir, "09_pathogen_endpoint_summary.csv"), index=False)

# 4d. Checks that would otherwise fail silently.
# The join between the two naming conventions is the one place this could be wrong without raising:
# a mismatched key would give every endpoint someone else's Jaccard. The two metrics measure the
# same thing by different means, so their specificities must be positively correlated.
rho, p = spearmanr(table["jac_specificity"], table["auroc_out_specificity"], nan_policy="omit")
print(f"[endpoint-quality] Jaccard vs AUROC specificity: Spearman rho = {rho:.3f} (p = {p:.2e})")
assert rho > 0, ("Jaccard and AUROC specificity are not positively correlated — the endpoint key "
                 "join is likely wrong.")

# Known answer, taken over the whole MTB block rather than a few hand-copied values. The reference
# is computed straight off the fresh self_perf frame using ITS OWN organism columns, which is an
# independent path from the pathogen mapping this script builds (config organism -> code -> named
# column -> endpoint key): if either the node selection or the same/different split were wrong, the
# two would disagree. MTB is the reference block because it is the largest same-organism set (40
# endpoints, 39 peers each) and the one the analysis was motivated by.
observed = table.set_index("endpoint")["auroc_out_same_median"]
mtb = self_perf[(self_perf["predictor_organism"] == "Mycobacterium tuberculosis")
                & (self_perf["target_organism"] == "Mycobacterium tuberculosis")
                & ~self_perf["self_pair"]]
reference = mtb.groupby("predictor_endpoint")["value"].median()
assert set(reference.index) == set(table.loc[table["pathogen"] == "mtuberculosis", "endpoint"]), (
    "the MTB endpoints in scope differ from those self_perf labelled M. tuberculosis")
delta = (observed.reindex(reference.index) - reference).abs().max()
assert delta < 1e-9, f"MTB medians differ from a direct groupby by up to {delta:.2e}"
print(f"[endpoint-quality] known-answer check passed on all {len(reference)} MTB endpoints "
      f"(max deviation {delta:.1e} from a direct groupby)")

# 4e. Report.
below = table[table["auroc_out_same_median"] < PREDICTOR_CHANCE_LEVEL]
print(f"\n[endpoint-quality] {len(below)} endpoint(s) rank their own pathogen's actives below "
      f"chance (median AUROC < {PREDICTOR_CHANCE_LEVEL}):")
for _, r in below.iterrows():
    print(f"    {r['organism']:<28} {r['model_id']}:{r['column_name']:<26} "
          f"AUROC {r['auroc_out_same_median']:.4f}  over {int(r['n_peers'])} peers "
          f"({int(r['n_same_model_peers'])} same-model)")

cols = ["organism", "model_id", "column_name", "is_consensus", "n_peers", "n_same_model_peers",
        "auroc_out_same_median", "auroc_out_specificity", "jac_specificity"]
print("\n[endpoint-quality] weakest 15 endpoints by median AUROC against their own peers:")
print(table[cols].head(15).to_string(index=False, float_format=lambda v: f"{v:.4f}"))

print("\n[endpoint-quality] per-pathogen summary (weakest first):")
print(endpoint_summary.to_string(index=False, float_format=lambda v: f"{v:.4f}"))

# A consensus column should be at or near the top of its own pathogen's ranking; where it is not,
# that is a finding about the model's aggregation, not a bug in this script.
cons = table[table["is_consensus"]]
# rank 1 is the WEAKEST, and a pathogen has n_peers + 1 endpoints, so this is "in the bottom half
# of its own pathogen".
odd = cons[cons["rank_within_pathogen"] <= (cons["n_peers"] + 1) / 2]
if len(odd):
    print(f"\n[endpoint-quality] {len(odd)} consensus column(s) in the WEAKER half of their own "
          "pathogen's ranking — worth a look, since a consensus aggregates the very peers it is "
          "scored against:")
    print(odd[["organism", "model_id", "rank_within_pathogen", "n_peers",
               "auroc_out_same_median"]].to_string(index=False,
                                                   float_format=lambda v: f"{v:.4f}"))

print()
save_endpoint_quality_figures(output_dir, table)

print(f"\nDone → {output_dir}")
