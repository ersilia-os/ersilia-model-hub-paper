"""Step 14 analysis engine — pathogen hits vs antibiotic-resemblance hits on the library UMAP.

Pure set arithmetic over two summary CSVs that already exist; nothing here reads a prediction file
or the 1.35M-row matrices:

    step 10  ``10_top1000_per_pathogen.csv``      each pathogen's top-N by ``consensus_score``
    step 12  ``12_abx_umap_highlights.csv``       each abx endpoint's highlighted compounds

For every (pathogen, abx endpoint) pair the two ``key`` sets are split three ways — pathogen only,
abx only, and the intersection — and each compound keeps the UMAP coordinates already carried by
the source tables. The intersection is the point of the figure: compounds a pathogen model ranks
highly AND that resemble a known antibiotic class.

The two sides are NOT selected the same way, and the difference matters when reading a panel:

*   The pathogen side IS a rank cutoff — the top ``PROJECTION_TOP_N`` by ``consensus_score``, so
    every pathogen contributes exactly ``PROJECTION_TOP_N`` compounds.
*   The abx side is NOT. Only ``abx_score`` is continuous; the rest are binary flags or small
    integer counts, so step 12 takes every compound with a value > 0, capped at
    ``PROJECTION_TOP_N``. An endpoint with 346 positives library-wide contributes 346 compounds,
    and one whose cap binds contributes an arbitrary 1000 of its positives.

So an intersection of zero can mean "no overlap" or simply "this endpoint only has 3 hits at all".
:func:`overlap_counts` therefore carries ``n_abx`` and ``abx_capped`` on every row.
"""

import os

import pandas as pd

from default import ABX_OVERLAP_ENDPOINTS, PROJECTION_TOP_N
from eval_abx_matrix import GROUP_CODE as ABX_PREFIX

#: Column pair carrying the shared UMAP coordinates in both source tables.
COORD_COLS = ("umap_x", "umap_y")


def selected_abx_endpoints(config_csv, columns=ABX_OVERLAP_ENDPOINTS, prefix=ABX_PREFIX):
    """Resolve the requested abx column names against ``config/antibiotic_resemblance.csv``.

    Returns a DataFrame with ``model_id``, ``column_name`` and ``endpoint`` (the prefixed
    ``{prefix}__{model_id}__{column_name}`` name step 12 writes), in the order given by
    ``columns``. Raises on a name that is not a selected endpoint in the config, rather than
    silently dropping a requested panel — a typo'd endpoint must fail loudly.
    """
    config = pd.read_csv(config_csv)
    sel = config[config["selected"] == "Yes"]
    by_column = {r.column_name: r for r in sel.itertuples()}

    unknown = [c for c in columns if c not in by_column]
    if unknown:
        raise ValueError(
            f"Requested abx endpoint(s) {unknown} are not selected rows of "
            f"{os.path.basename(config_csv)}. Available: {sorted(by_column)}")

    rows = [{"model_id": by_column[c].model_id, "slug": by_column[c].slug, "column_name": c,
             "endpoint": f"{prefix}__{by_column[c].model_id}__{c}"} for c in columns]
    return pd.DataFrame(rows)


def _keys(df):
    """The ``key`` column as a set, tolerating an empty frame."""
    return set(df["key"]) if len(df) else set()


def overlap_table(pathogen_top, abx_highlights, pathogens, endpoints, coord_cols=COORD_COLS):
    """One row per (pathogen, endpoint, compound) with a ``group`` of pathogen/abx/both.

    Coordinates are taken from whichever source table the compound came from; for a ``both``
    compound the pathogen table's are used (they are the same eos1klk layout, so this is a choice
    of provenance, not of position).
    """
    xcol, ycol = coord_cols
    frames = []
    for p in pathogens.itertuples():
        pat = pathogen_top[pathogen_top["pathogen_code"] == p.code]
        pat_keys = _keys(pat)
        for e in endpoints.itertuples():
            abx = abx_highlights[abx_highlights["endpoint"] == e.endpoint]
            abx_keys = _keys(abx)
            both = pat_keys & abx_keys

            parts = [
                pat[pat["key"].isin(pat_keys - abx_keys)].assign(group="pathogen"),
                abx[abx["key"].isin(abx_keys - pat_keys)].assign(group="abx"),
                pat[pat["key"].isin(both)].assign(group="both"),
            ]
            block = pd.concat([q[["key", xcol, ycol, "group"]] for q in parts], ignore_index=True)
            block.insert(0, "pathogen_code", p.code)
            block.insert(1, "endpoint", e.endpoint)
            frames.append(block)
    return pd.concat(frames, ignore_index=True)


def overlap_counts(table, pathogen_top, abx_highlights, pathogens, endpoints,
                   top_n=PROJECTION_TOP_N):
    """Per-(pathogen, endpoint) set sizes, plus the Jaccard index of the two sets.

    ``abx_capped`` flags an endpoint whose highlight list hit step 12's cap, i.e. one whose drawn
    compounds are an arbitrary subset of its positives — the panels where a small intersection is
    least interpretable.
    """
    rows = []
    for p in pathogens.itertuples():
        n_pat = int((pathogen_top["pathogen_code"] == p.code).sum())
        for e in endpoints.itertuples():
            n_abx = int((abx_highlights["endpoint"] == e.endpoint).sum())
            g = table[(table["pathogen_code"] == p.code) & (table["endpoint"] == e.endpoint)]
            n_both = int((g["group"] == "both").sum())
            union = n_pat + n_abx - n_both
            rows.append({
                "pathogen_code": p.code, "pathogen": p.pathogen,
                "endpoint": e.endpoint, "model_id": e.model_id, "column_name": e.column_name,
                "n_pathogen": n_pat, "n_abx": n_abx, "n_both": n_both,
                "jaccard": n_both / union if union else float("nan"),
                "abx_capped": n_abx >= top_n,
            })
    return pd.DataFrame(rows)


def run_all(pathogen_top_csv, abx_highlights_csv, pathogens_csv, abx_config_csv, output_dir,
            top_n=PROJECTION_TOP_N):
    """Overlap orchestrator (step 11). Writes ``11_overlap_points.csv`` and ``11_overlap_counts.csv``."""
    pathogens = pd.read_csv(pathogens_csv)
    endpoints = selected_abx_endpoints(abx_config_csv)
    pathogen_top = pd.read_csv(pathogen_top_csv)
    abx_highlights = pd.read_csv(abx_highlights_csv)

    print(f"[overlap] {len(pathogens)} pathogens x {len(endpoints)} abx endpoints")
    for e in endpoints.itertuples():
        n = int((abx_highlights["endpoint"] == e.endpoint).sum())
        cap = " (cap binds — arbitrary subset of its positives)" if n >= top_n else ""
        print(f"  {e.model_id}/{e.column_name}: {n} compounds{cap}")

    missing = [e.endpoint for e in endpoints.itertuples()
               if not (abx_highlights["endpoint"] == e.endpoint).any()]
    if missing:
        raise ValueError(f"No rows in {os.path.basename(abx_highlights_csv)} for {missing} — "
                         "re-run 11_abx_resemblance_matrix.py.")

    table = overlap_table(pathogen_top, abx_highlights, pathogens, endpoints)
    counts = overlap_counts(table, pathogen_top, abx_highlights, pathogens, endpoints, top_n=top_n)

    points_path = os.path.join(output_dir, "11_overlap_points.csv")
    counts_path = os.path.join(output_dir, "11_overlap_counts.csv")
    table.to_csv(points_path, index=False)
    counts.to_csv(counts_path, index=False)
    print(f"  -> {os.path.basename(points_path)} ({len(table)} rows)")
    print(f"  -> {os.path.basename(counts_path)} ({len(counts)} rows)")

    top = counts.nlargest(10, "n_both")[["pathogen_code", "column_name", "n_both", "n_abx",
                                         "jaccard", "abx_capped"]]
    print("\n  Largest intersections (n_both):")
    print("   " + top.to_string(index=False).replace("\n", "\n   "))
    return pathogens, endpoints, counts
