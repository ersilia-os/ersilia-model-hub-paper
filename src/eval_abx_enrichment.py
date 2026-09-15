"""Step 12 analysis engine — Fisher enrichment of pathogen top-1000 sets against eos19mt's
antibiotic structural classes.

Pure statistics over two things script 12 already has: the named abx matrix built in its section 1
(:func:`eval_abx_matrix.build_abx_named_matrix`) and each pathogen's top-N compounds from step 11
(``11_top{PROJECTION_TOP_N}_per_pathogen.csv``). For every (pathogen, eos19mt class) pair, tests
whether being in the pathogen's top-N is associated with belonging to the class, via a one-sided
Fisher exact test over the full reference library.

Rows are eos19mt's 38 classes in the model's own (== the config file's) order; columns are the 15
pathogens in step 10's phylogenetic order — reusing :func:`eval_auroc_matrix.bioactivity_order` /
:func:`eval_auroc_matrix.organism_order` / :func:`eval_auroc_matrix.phylogeny_organism_order`
rather than re-deriving it, so this matrix's column order can never drift from the AUROC matrix's
row order.

Contingency table, per pair, over the full library:

                      in top-N     not in top-N
    class member         a              c
    not class member     b              d

Compounds with a missing value in that class's column are excluded from the pair's table entirely
(never imputed to 0 — the class's own missingness). A class with zero variance over the resulting
universe (e.g. one of the four constant-zero eos19mt columns) has nothing to test — its
``odds_ratio``/``p_value`` are ``NaN`` rather than a computed-but-meaningless value.
"""

import os

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

from default import ABX_ENRICHMENT_ALTERNATIVE, ABX_ENRICHMENT_MODEL_ID, ABX_ENRICHMENT_SIG_THRESHOLDS
from eval_abx_matrix import GROUP_CODE, named_column
from eval_auroc_matrix import bioactivity_order, organism_order, phylogeny_organism_order


def eos19mt_classes(sel, model_id=ABX_ENRICHMENT_MODEL_ID, group=GROUP_CODE):
    """eos19mt's selected rows of ``sel`` (as returned by
    :func:`eval_abx_matrix.build_abx_named_matrix`), in file/native order.

    File order matches the model's own output order — confirmed against
    ``data/processed/annotation_preds_ref_library/eos19mt_v2.csv``'s header — so this IS "the same
    order as in the eos19mt model", not an editorial re-sort.
    """
    classes = sel[sel["model_id"] == model_id].reset_index(drop=True).copy()
    if not len(classes):
        raise ValueError(f"No selected rows for model_id={model_id!r} in the abx selection config.")
    classes["endpoint"] = [named_column(r.model_id, r.column_name, group)
                           for r in classes.itertuples()]
    return classes


def pathogen_phylo_order(selection_csv, pathogens_csv, taxonomy_csv):
    """The 15 pathogens, with ``code``/``pathogen``/``organism_class``, in step 10's phylogenetic
    order — the same three ``eval_auroc_matrix`` calls ``scripts/10_auroc_matrix.py`` makes, so this
    matrix's column order can never disagree with the AUROC matrix's row order. ``available=None``
    (no step-07 parquet cache filter): every pathogen is one of the fixed 15 regardless of which
    endpoints happen to be cached, so no dependency on that cache is needed here.
    """
    endpoints = bioactivity_order(selection_csv, pathogens_csv, available=None)
    organisms = organism_order(endpoints)
    order = phylogeny_organism_order(organisms, taxonomy_csv)

    codes = pd.read_csv(pathogens_csv).set_index("pathogen")["code"]
    class_of = organisms.set_index("organism")["organism_class"]
    rows = [{"pathogen": o, "code": codes[o], "organism_class": class_of[o]} for o in order]
    return pd.DataFrame(rows)


def _fisher_cell(col, top_keys, alternative=ABX_ENRICHMENT_ALTERNATIVE):
    """One (pathogen, class) contingency table + Fisher test, from a class's full-library Series
    (indexed by ``key``) and a pathogen's top-N ``key`` set.

    NaN compounds are dropped from this pair's table before anything is counted. Returns the four
    contingency counts plus ``odds_ratio``/``p_value`` — both ``NaN`` when the class has no
    variance over the resulting universe (constant zero, or constant nonzero).
    """
    present = col.dropna()
    n_nan = int(col.isna().sum())
    in_top = present.index.isin(top_keys)
    is_member = present.to_numpy() > 0

    a = int(np.sum(in_top & is_member))
    b = int(np.sum(in_top & ~is_member))
    c = int(np.sum(~in_top & is_member))
    d = int(np.sum(~in_top & ~is_member))

    n_member = a + c
    if n_member == 0 or n_member == len(present):
        odds_ratio, p_value = np.nan, np.nan
    else:
        odds_ratio, p_value = fisher_exact([[a, b], [c, d]], alternative=alternative)
    return {"a": a, "b": b, "c": c, "d": d, "n_class_nan": n_nan,
            "odds_ratio": odds_ratio, "p_value": p_value}


def fisher_table(matrix, classes, pathogen_top, pathogens, alternative=ABX_ENRICHMENT_ALTERNATIVE):
    """Long-format result: one row per (pathogen, class) pair.

    ``matrix`` is the full-library named abx matrix (script 12 section 1, indexed by ``key``).
    ``pathogen_top`` is step 11's ``11_top{N}_per_pathogen.csv``. ``pathogens``/``classes`` fix the
    row/column order this table is later pivoted into.
    """
    rows = []
    for p in pathogens.itertuples():
        top_keys = set(pathogen_top.loc[pathogen_top["pathogen_code"] == p.code, "key"])
        for cls in classes.itertuples():
            cell = _fisher_cell(matrix[cls.endpoint], top_keys, alternative=alternative)
            rows.append({"pathogen_code": p.code, "pathogen": p.pathogen,
                        "model_id": cls.model_id, "column_name": cls.column_name,
                        "endpoint": cls.endpoint, "n_pathogen_top": len(top_keys), **cell})
    return pd.DataFrame(rows)


def benjamini_hochberg(p_values):
    """FDR-adjusted p-values (Benjamini-Hochberg, 1995). NaN entries stay NaN and do not count
    toward ``m``. No new dependency: ``scipy.stats`` has no BH implementation and ``statsmodels``
    is not a project dependency.
    """
    p = np.asarray(p_values, dtype=float)
    out = np.full_like(p, np.nan)
    finite = np.isfinite(p)
    m = int(finite.sum())
    if m == 0:
        return out
    idx = np.flatnonzero(finite)
    order = idx[np.argsort(p[idx])]
    ranked = p[order] * m / (np.arange(m) + 1)
    adjusted = np.minimum.accumulate(ranked[::-1])[::-1]  # enforce monotonicity, largest p first
    out[order] = np.minimum(adjusted, 1.0)
    return out


def enrichment_matrices(long_df, classes, pathogens):
    """``(odds_ratio_df, p_value_df)`` — wide pivots, rows in ``classes`` order (``column_name``),
    columns in ``pathogens`` order (``code``).
    """
    def _pivot(value_col):
        wide = long_df.pivot(index="column_name", columns="pathogen_code", values=value_col)
        return wide.reindex(index=classes["column_name"], columns=pathogens["code"])
    return _pivot("odds_ratio"), _pivot("p_value")


def run_all(matrix, sel, pathogen_top_csv, selection_csv, pathogens_csv, taxonomy_csv, output_dir,
           model_id=ABX_ENRICHMENT_MODEL_ID, alternative=ABX_ENRICHMENT_ALTERNATIVE):
    """Orchestrator (step 12). Writes ``12_abx_enrichment_long.csv``,
    ``12_abx_enrichment_odds_ratio.csv`` and ``12_abx_enrichment_pvalue.csv``.
    """
    classes = eos19mt_classes(sel, model_id=model_id)
    pathogens = pathogen_phylo_order(selection_csv, pathogens_csv, taxonomy_csv)
    pathogen_top = pd.read_csv(pathogen_top_csv)

    print(f"[abx-enrichment] {len(pathogens)} pathogens x {len(classes)} {model_id} classes "
          f"(Fisher exact, alternative={alternative!r})")
    long_df = fisher_table(matrix, classes, pathogen_top, pathogens, alternative=alternative)
    long_df["p_value_fdr"] = benjamini_hochberg(long_df["p_value"])

    odds_ratio, p_value = enrichment_matrices(long_df, classes, pathogens)

    long_path = os.path.join(output_dir, "12_abx_enrichment_long.csv")
    or_path = os.path.join(output_dir, "12_abx_enrichment_odds_ratio.csv")
    p_path = os.path.join(output_dir, "12_abx_enrichment_pvalue.csv")
    long_df.to_csv(long_path, index=False)
    odds_ratio.to_csv(or_path)
    p_value.to_csv(p_path)
    print(f"  -> {os.path.basename(long_path)} ({len(long_df)} rows)")
    print(f"  -> {os.path.basename(or_path)} / {os.path.basename(p_path)} {odds_ratio.shape}")

    degenerate = sorted(long_df.loc[long_df['odds_ratio'].isna(), 'column_name'].unique())
    if degenerate:
        print(f"  {len(degenerate)} class(es) have no variance library-wide (NaN throughout): "
              f"{degenerate}")

    sig_at = ABX_ENRICHMENT_SIG_THRESHOLDS[0]
    sig = long_df.dropna(subset=["p_value"])
    sig = sig[sig["p_value"] < sig_at]
    print(f"  {len(sig)}/{len(long_df)} pairs significant at p<{sig_at} (raw, uncorrected)")

    top = long_df.replace([np.inf], np.nan).nlargest(
        10, "odds_ratio")[["pathogen_code", "column_name", "odds_ratio", "p_value", "a", "c"]]
    print("\n  Largest finite odds ratios:")
    print("   " + top.to_string(index=False).replace("\n", "\n   "))
    return classes, pathogens, long_df, odds_ratio, p_value
