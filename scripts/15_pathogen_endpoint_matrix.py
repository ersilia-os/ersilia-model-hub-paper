"""Step 15 — one pathogen's own endpoint x endpoint top-1000 shared-actives matrix.

Step 09 computes and caches the full 307x307 top-1000 Jaccard matrix across every selected
bioactivity endpoint, but only ever consumes it aggregated into same-pathogen / different-pathogen
box plots. This step draws the raw endpoint x endpoint matrix for a SINGLE pathogen at a time — the
detailed, per-endpoint view that a pathogen-level box plot averages away.

**Scope: `RUNS` below, not all 15 pathogens of interest.** This is the start of a "detailed analysis
of enrichment for the endpoints of one single organism" — E. coli was chosen as the first worked
example (24 selected bioactivity endpoints across 7 models, the widest spread among the pathogens of
interest with more than a handful of endpoints), then P. falciparum (malaria) was added the same
way. Looping this over the remaining pathogens (`config/pathogens_of_interest.csv`) is a natural
follow-up once the shape of the analysis is agreed, not done here.

**Cell value: the RAW top-1000 shared-actives COUNT, not Jaccard (user-directed, 2026-09-03).** How
many of the two endpoints' top-1000 highest-scoring compounds are literally the same compounds — the
same quantity step 10's own overlap matrix draws (`eval_auroc_matrix.overlap_matrix`), here computed
generically per column pair via :func:`eval_correlations.topn_overlap_matrix`. Most of the
bioactivity endpoints need no new computation for the UNDERLYING top-1000 sets: this reads step 09's
already-cached, already-validated Jaccard matrix and slices it down to one pathogen's own columns —
but since that cache stores the Jaccard RATIO, not the count, getting the count means one fresh
top-1000 computation over the pathogen's raw columns (read from step 07's parquet). The result is
cross-checked against the step-09 cache: raw counts are converted back to Jaccard
(`jaccard = inter / (2 * CUTOFF - inter)`, exact here since every endpoint has far more than 1000
scored compounds) and compared to the cached ratio before anything is trusted.

**An extra predictor node, one per pathogen, per `RUNS[*]["predictor_family"]` (user-directed):
`abx` (resemblance-to-known-antibiotics) for E. coli, `cytotox` (cytotoxicity) for P. falciparum.**
Neither is in step 09's cache — step 09 is bioactivity-only by design — so this reuses step 10's own
merge (`eval_auroc_matrix.merged_predictor_scores`, the same rank-summed
`{family}__merged__rank_sum` column step 10's AUROC/overlap matrices carry) against the same raw
endpoint columns, folded into the one fresh top-1000 computation above.

**P. falciparum's `consensus_only` collapse (user-directed):** P. falciparum has 64 selected
bioactivity endpoints — `eos4an7` alone contributes 24 raw ChEMBL/PubChem sub-assay columns plus its
own `consensus_score`. At the fixed 180x180 mm page-width footprint this step draws at, a 65x65
matrix would put each cell at ~2.8 mm, too small to keep the printed values legible. Rather than
resize the figure, any model with a `consensus_score` among its selected columns is collapsed to
just that one column (`_select_endpoints`) — `eos4an7` goes from 24 columns to 1, bringing
P. falciparum to 13 bioactivity nodes + 1 predictor = 14, comparable in scale to E. coli's 25. This
rule is scoped to whichever runs opt into it via `consensus_only=True` in `RUNS`; E. coli's `False`
entry draws every one of its 24 bioactivity endpoints.

    python 07_score_matrices.py           # writes the parquet cache (if not already run)
    python 08_property_matrices.py        # abx + cytotox matrices (if not already run)
    python 09_bioactivity_endpoints.py    # writes the cached 307x307 Jaccard matrix (if not already run)
    python 15_pathogen_endpoint_matrix.py

**Chord (circos-style) diagrams, one per `RUNS[*]["chord"]` entry (currently 2, both E. coli only),
are ALTERNATIVE views of the exact same matrix — no new data.** The 24 bioactivity endpoints sit on
a ring; a ribbon between two endpoints is coloured/widened by their shared-actives count on the
identical `Normalize(0, CUTOFF)` scale the matrix's own colorbar uses. The predictor node is NOT a
25th ribbon node — it is a radial bar track around the ring instead (one bar per endpoint), so its
relationship to every endpoint reads as a single glance rather than one more crowded ribbon. Nothing
is filtered: a zero-count pair is simply not drawn (invisible either way), so an endpoint with
nothing to show shows exactly that.

**Two ring-ordering modes (user-directed), since the matrix's own model-grouped order made the first
version look like an arbitrary mandala rather than showing structure:** `"cluster"` (hierarchical
clustering / optimal leaf ordering on the shared-actives counts, so strongly-linked endpoints sit
next to each other) and `"abx"` (a direct descending rank by overlap with the predictor node,
clockwise from 12 o'clock). See `src/plots_pathogen_endpoint_matrix.py`'s `PathogenChordPlot` and
`_ordered_labels` for the full layout and ordering logic.

**The chord figures ALSO collapse each model to its consensus score (user-directed) —
`_collapse_to_consensus`, applied ONLY to `chord_sub`/`chord_model_id`, never to `sub` itself, so the
matrix keeps every endpoint.** E. coli's `eos5eya` (11 raw ChEMBL/PubChem sub-columns + its own
`consensus_score`) collapses to 1 node instead of 13, since a reader comparing 24 crowded ribbons
gets more signal from one ribbon per model than from the same model's sub-assays talking mostly to
each other. A model with no consensus column (and the predictor node, which never has one) is
unaffected. For P. falciparum this would be a no-op if its chord list were ever populated — its
`consensus_only` matrix-level collapse already did the same thing to `eos4an7`.

Outputs
-------
    output/15_pathogen_endpoint_matrix/15_{code}_overlap_top1000_matrix.csv           (per RUNS entry)
    output/15_pathogen_endpoint_matrix/png|pdf/15_{code}_endpoint_overlap_top1000_matrix.*
    output/15_pathogen_endpoint_matrix/png|pdf/15_{code}_endpoint_chord_top1000_{order}.*  (per chord entry)
    output/15_pathogen_endpoint_matrix/figure_cells.json
"""

import os
import sys

import numpy as np
import pandas as pd

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))

from default import ACTIVITY_BINARIZE_TOP_N, MERGED_PREDICTOR_GROUPS  # noqa: E402
from eval_auroc_matrix import merged_predictor_scores  # noqa: E402
from eval_correlations import parse_named_column, topn_overlap_matrix  # noqa: E402
from eval_predictor_performance import _assert_key_alignment  # noqa: E402
from plots_pathogen_endpoint_matrix import (  # noqa: E402
    save_pathogen_chord_figure, save_pathogen_endpoint_matrix_figure,
)

#: One entry per pathogen to draw. See the module docstring for why each is in scope, what
#: `consensus_only` does, why the two pathogens carry different `predictor_family` values, and why
#: only E. coli also gets chord diagrams. `chord` is a list of ring-ordering modes (empty = no chord
#: figure) — see `plots_pathogen_endpoint_matrix._ordered_labels` for what each mode means. Order is
#: the order figures are drawn/printed in.
RUNS = [
    {"code": "ecoli", "consensus_only": False, "predictor_family": "abx",
     "chord": ["cluster", "abx"]},
    {"code": "pfalciparum", "consensus_only": True, "predictor_family": "cytotox", "chord": []},
]

CUTOFF = 1000

config_dir = os.path.join(root, "..", "config")
step07_dir = os.path.join(root, "..", "output", "07_score_matrices")
step08_dir = os.path.join(root, "..", "output", "08_property_matrices")
step09_dir = os.path.join(root, "..", "output", "09_bioactivity_endpoints")
output_dir = os.path.join(root, "..", "output", "15_pathogen_endpoint_matrix")
os.makedirs(output_dir, exist_ok=True)

#: Where each predictor family's raw per-compound matrix lives (step 08). Only the families
#: actually referenced in `RUNS` are checked for existence / key-alignment below — a family not in
#: `MERGED_PREDICTOR_GROUPS` (`default.py`) would raise a plain `KeyError` from that dict, same as
#: any other config typo in this script.
PREDICTOR_CSV = {
    "abx": os.path.join(step08_dir, "08_abx_matrix_named.csv"),
    "cytotox": os.path.join(step08_dir, "08_cytotox_matrix_named.csv"),
}

endpoint_selection_path = os.path.join(config_dir, "08_endpoint_selection.csv")
pathogens_of_interest_path = os.path.join(config_dir, "pathogens_of_interest.csv")
jaccard_cache_path = os.path.join(step09_dir, f"09_jaccard_top{CUTOFF}_baseline_matrix.csv")
parquet_path = os.path.join(step07_dir, "07_score_matrix_full.parquet")

assert CUTOFF == ACTIVITY_BINARIZE_TOP_N, (
    f"Cutoff {CUTOFF} != ACTIVITY_BINARIZE_TOP_N {ACTIVITY_BINARIZE_TOP_N}: this step must read "
    "the same cutoff step 09 cached, or the two would not be comparable.")

#: Only the predictor families actually used by `RUNS`, deduplicated — so adding a pathogen that
#: reuses a family already in scope doesn't re-check or re-verify that family's CSV twice.
needed_families = list(dict.fromkeys(r["predictor_family"] for r in RUNS))
needed_csvs = [PREDICTOR_CSV[f] for f in needed_families]

for path, step in [(jaccard_cache_path, "09_bioactivity_endpoints.py"),
                   (parquet_path, "07_score_matrices.py"),
                   *[(p, "08_property_matrices.py") for p in needed_csvs]]:
    if not os.path.exists(path):
        sys.exit(f"Missing {path}. Run `python {step}` first.")

jaccard = pd.read_csv(jaccard_cache_path, index_col=0)
print(f"[pathogen-endpoint-matrix] loaded cached {jaccard.shape} Jaccard matrix from "
      f"{os.path.basename(jaccard_cache_path)} (used only to cross-check the fresh overlap counts)")

pathogens = pd.read_csv(pathogens_of_interest_path)
selection = pd.read_csv(endpoint_selection_path)

for family, path in zip(needed_families, needed_csvs):
    _assert_key_alignment([path], parquet_path)
    print(f"[pathogen-endpoint-matrix] {family}: key order verified against {parquet_path}")


def _select_endpoints(expected, consensus_only):
    """``expected`` unchanged, unless ``consensus_only``: any model with a ``consensus_score``
    among its selected columns is collapsed to just that one row; a model with no consensus column
    keeps every row it had. See the module docstring for why (P. falciparum's node-count blowup).
    """
    if not consensus_only:
        return expected
    rows = []
    for _, group in expected.groupby("model_id", sort=False):
        if "consensus_score" in set(group["column_name"]):
            group = group[group["column_name"] == "consensus_score"]
        rows.append(group)
    return pd.concat(rows, ignore_index=True)


def _collapse_to_consensus(sub, model_id):
    """Row/column-filtered copy of ``sub`` (plus the matching ``model_id``): any model with a
    ``"{model}:consensus_score"`` label among its columns is collapsed to just that one label; a
    model with no consensus column (including the predictor node, which never has one) keeps every
    label it had. Same rule as :func:`_select_endpoints`'s ``consensus_only``, applied post-hoc here
    to the CHORD figures only (user-directed) — the matrix keeps every endpoint regardless.
    """
    by_model = {}
    for lab, m in zip(sub.columns, model_id):
        by_model.setdefault(m, []).append(lab)
    keep = set()
    for labs in by_model.values():
        consensus = [lab for lab in labs if lab.endswith(":consensus_score")]
        keep.update(consensus if consensus else labs)
    kept_labels = [lab for lab in sub.columns if lab in keep]
    label_model = dict(zip(sub.columns, model_id))
    kept_model_id = [label_model[lab] for lab in kept_labels]
    return sub.loc[kept_labels, kept_labels], kept_model_id


def run_pathogen(code, consensus_only, predictor_family, chord=()):
    pathogen_name = pathogens.set_index("code").loc[code, "pathogen"]
    predictor_csv = PREDICTOR_CSV[predictor_family]
    predictor_label = f"{predictor_family}:rank_sum"

    # Cross-checked against config, not a hardcoded count — the same "assert against config"
    # pattern steps 09/10 use throughout, so an upstream selection change fails loudly here rather
    # than silently drawing a matrix of the wrong size. Checked on the FULL selected set, before
    # any consensus_only display-subsetting, since this validates the cache itself.
    nodes_all = [c for c in jaccard.columns if parse_named_column(c)[0] == code]
    expected_all = selection[(selection["selected"] == "Yes") & (selection["organism"] == pathogen_name)]
    assert len(nodes_all) == len(expected_all), (
        f"{code}: {len(nodes_all)} node(s) in the cached matrix vs {len(expected_all)} selected "
        f"endpoint(s) in {endpoint_selection_path} for {pathogen_name!r} — selection and cache "
        "disagree.")

    expected = _select_endpoints(expected_all, consensus_only)
    if consensus_only:
        print(f"[pathogen-endpoint-matrix] {code}: consensus_only collapsed "
              f"{len(expected_all)} -> {len(expected)} bioactivity endpoints")

    nodes = [f"{code}__{r.model_id}__{r.column_name}" for r in expected.itertuples()]
    model_id = expected["model_id"].tolist()
    column_name = expected["column_name"].tolist()
    labels = [f"{m}:{c}" for m, c in zip(model_id, column_name)]
    jaccard_ref = jaccard.loc[nodes, nodes]
    jaccard_ref.index = labels
    jaccard_ref.columns = labels

    meta = expected.set_index(["model_id", "column_name"])[["sensitivity"]]
    sensitivity = [meta.loc[(m, c), "sensitivity"] for m, c in zip(model_id, column_name)]

    # --- Add the pathogen's own aggregate predictor as an extra node (user-directed), and compute
    # RAW top-1000 shared-actives counts (user-directed) instead of reusing the step-09 Jaccard
    # cache. Key order must match across the parquet and the predictor CSV, or the two would be
    # combined positionally-wrong — the same guard step 10 runs before doing exactly this same
    # merge (checked once per family, above, since it does not depend on the pathogen).
    old_cols = [f"{m}:{c}" for m, c in zip(model_id, column_name)]
    raw = pd.read_parquet(parquet_path, columns=old_cols).reset_index(drop=True)
    assert list(raw.columns) == labels, "raw parquet column order does not match the cached labels"
    raw.columns = labels

    # Only `predictor_family` is merged — this pathogen's `RUNS` entry picks exactly one family.
    merged = merged_predictor_scores([predictor_csv],
                                     groups={predictor_family: MERGED_PREDICTOR_GROUPS[predictor_family]})
    merged = merged.rename(columns={f"{predictor_family}__merged__rank_sum": predictor_label})

    combined = pd.concat([raw, merged], axis=1)
    del raw, merged
    full_overlap = topn_overlap_matrix(combined, CUTOFF)
    del combined

    # The bioactivity block of this fresh computation must reproduce step 09's cache exactly, once
    # converted back to a Jaccard ratio — checked, not assumed, before the extra row/column is
    # trusted alongside it. Exact here (not merely close) because every endpoint in this pipeline
    # has far more than CUTOFF scored compounds, so both top-CUTOFF sets always have exactly CUTOFF
    # members and `jaccard = inter / (2 * CUTOFF - inter)` holds without approximation.
    inter = full_overlap.loc[labels, labels].to_numpy()
    with np.errstate(divide="ignore", invalid="ignore"):
        computed_jaccard = inter / (2 * CUTOFF - inter)
    delta = computed_jaccard - jaccard_ref.to_numpy()
    assert np.nanmax(np.abs(delta)) < 1e-6, (
        f"Freshly computed {code} block disagrees with the step-09 cache by up to "
        f"{np.nanmax(np.abs(delta)):.2e} (Jaccard-equivalent) — the parquet read is not "
        "reproducing the cached matrix.")
    print(f"[pathogen-endpoint-matrix] {code} block cross-checked against the step-09 cache "
          f"(max |delta| {np.nanmax(np.abs(delta)):.2e}, Jaccard-equivalent)")

    labels = labels + [predictor_label]
    model_id = model_id + [predictor_family]
    sensitivity = sensitivity + ["n/a (predictor, not an assay)"]
    sub = full_overlap.loc[labels, labels]

    csv_path = os.path.join(output_dir, f"15_{code}_overlap_top{CUTOFF}_matrix.csv")
    sub.to_csv(csv_path)
    n_bio = sub.shape[0] - 1
    n_models = len(set(model_id) - {predictor_family})
    print(f"[pathogen-endpoint-matrix] {pathogen_name} ({code}): {sub.shape[0]} nodes "
          f"({n_bio} bioactivity endpoints across {n_models} models + 1 aggregate "
          f"{predictor_family} predictor) -> {os.path.basename(csv_path)}")

    save_pathogen_endpoint_matrix_figure(
        output_dir, sub, model_id, name=f"15_{code}_endpoint_overlap_top{CUTOFF}_matrix",
        top_n=CUTOFF)

    # --- Alternative view(s) of the SAME matrix (user-directed) — no new data, see the module
    # docstring for why the predictor node becomes a bar track instead of a 25th ribbon node, why
    # the chord figures collapse each model to its consensus score, and for what each ring-ordering
    # mode in `chord` means.
    if chord:
        chord_sub, chord_model_id = _collapse_to_consensus(sub, model_id)
        if len(chord_sub) < len(sub):
            print(f"[pathogen-endpoint-matrix] {code}: chord figures collapsed "
                  f"{len(sub)} -> {len(chord_sub)} nodes (consensus-only)")
        for order in chord:
            save_pathogen_chord_figure(
                output_dir, chord_sub, chord_model_id, predictor_family,
                name=f"15_{code}_endpoint_chord_top{CUTOFF}_{order}", top_n=CUTOFF, order=order)

    # --- Console summary: reported, not filtered or interpreted (per CLAUDE.md's sign-off rules). ---
    diag_mask = pd.DataFrame(np.eye(len(sub), dtype=bool), index=sub.index, columns=sub.columns)
    off_diag = sub.mask(diag_mask)
    stacked = off_diag.stack()
    stacked.index.names = ["endpoint", "peer"]
    stacked = stacked[stacked.index.get_level_values(0) < stacked.index.get_level_values(1)]

    print(f"\n[pathogen-endpoint-matrix] {len(stacked)} unordered pairs, "
          f"median {stacked.median():.0f} of {CUTOFF} shared, "
          f"range {stacked.min():.0f}-{stacked.max():.0f}")

    print("\n  strongest 5 pairs:")
    for (a, b), v in stacked.sort_values(ascending=False).head(5).items():
        print(f"    {a:<45} {b:<45} {v:.0f}")

    print("\n  weakest 5 pairs:")
    for (a, b), v in stacked.sort_values().head(5).items():
        print(f"    {a:<45} {b:<45} {v:.0f}")

    per_model = pd.Series(model_id).value_counts()
    print(f"\n  {len(per_model)} group(s) ({n_models} bioactivity model(s) + {predictor_family}): "
          + ", ".join(f"{m} ({n})" for m, n in per_model.items()))

    per_sensitivity = pd.Series(sensitivity).value_counts()
    print("  sensitivity: " + ", ".join(f"{s} ({n})" for s, n in per_sensitivity.items()))

    # A node (endpoint or the predictor_family node) whose median shared-actives count against
    # every OTHER node sits near zero is flagged here for visibility, not excluded or explained —
    # interpreting why is a call for the user. The matrix is symmetric, so a row median of the
    # diagonal-blanked matrix already IS the median against every other node, in both directions.
    # NEAR_ZERO_COUNT (5% of CUTOFF) is an informal report threshold, not a scientific cutoff —
    # nothing is dropped or filtered on it, only named in this printout — chosen to read at roughly
    # the same selectivity as the 0.05 Jaccard threshold this step used before switching to counts.
    NEAR_ZERO_COUNT = round(0.05 * CUTOFF)
    own_median = off_diag.median(axis=1, skipna=True)
    near_zero = own_median[own_median < NEAR_ZERO_COUNT].sort_values()
    if len(near_zero):
        print(f"\n  {len(near_zero)} node(s) with median shared-actives count < {NEAR_ZERO_COUNT} "
              f"(of {CUTOFF}) against every other {pathogen_name} node in this set (near-isolated, "
              "not filtered):")
        for ep, v in near_zero.items():
            print(f"    {ep:<45} median {v:.0f}")

    print()


for run in RUNS:
    run_pathogen(run["code"], run["consensus_only"], run["predictor_family"], run["chord"])

print(f"Done -> {output_dir}")
