"""Step 11 analysis engine — reference-library chemical-space projection, coloured by pathogen
activity.

``eos1klk`` computes four 2D projections (PCA, UMAP, t-SNE, TMAP) of the SAME ~1.35M-compound
reference library staged by ``00_download_data.py``, one output column pair (``{method}_x``,
``{method}_y``) per method. For each pathogen in ``config/pathogens_of_interest.csv`` this engine
finds its ``PROJECTION_TOP_N`` highest-scoring compounds by ``consensus_score`` — a rank cutoff,
never a score threshold — and, separately, a full-library density grid per method for the silver
background layer every pathogen's highlighted points sit on.

Memory: only one pathogen's ``key`` + ``consensus_score`` columns (of a ~424 MB / ~14-column
prediction file) are ever read at a time, via ``usecols``, immediately reduced to its top-N rows
and discarded. All 15 pathogens' raw scores, or all 15 prediction files, are never held in memory
together — only the small aggregated background grids and the (15 x PROJECTION_TOP_N)-row top-N
table are written out, per the repo's "feed figures from summary CSVs" rule.

This module also drives step 11's SECOND figure family (see :func:`run_coadd`): the same layout and
the same background grids, but highlighted by ``COADD_MODEL_ID`` (eos3dys, CoAdd-trained) rather
than the per-pathogen ChEMBL models. Because that model publishes no ``consensus_score``, its top-N
is ranked per ENDPOINT rather than per pathogen; the endpoints are then grouped onto organisms for
plotting.
"""

import os

import numpy as np
import pandas as pd

from default import (COADD_MODEL_ID, COADD_PROJECTION_METHOD, CORR_CHUNK_SIZE, PROJECTION_BINS,
                     PROJECTION_METHODS, PROJECTION_TOP_N)
from eval_correlations import latest_version_files

JOIN_COLS = ("key", "input")


# --------------------------------------------------------------------------- #
# Projection table                                                             #
# --------------------------------------------------------------------------- #
def load_projection(projection_file):
    """Load eos1klk's ``key`` + the 8 ``{method}_x/_y`` columns once.

    Small (~1.35M rows x 9 numeric columns), so this is read whole — unlike the pathogen
    prediction files, which are read column-restricted and in chunks (see
    :func:`_read_pathogen_scores`).
    """
    coord_cols = [f"{m}_{ax}" for m in PROJECTION_METHODS for ax in ("x", "y")]
    df = pd.read_csv(projection_file, usecols=["key"] + coord_cols)
    for c in coord_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def method_extent(proj, method):
    """``(xmin, xmax, ymin, ymax)`` for one projection method, ignoring NaN coordinates.

    Computed once per method and shared by the background grid and every pathogen's grid, so all
    16 grids for a method (1 background + 15 pathogens) align on the same cells and can be
    overlaid directly.
    """
    x, y = proj[f"{method}_x"], proj[f"{method}_y"]
    return float(x.min()), float(x.max()), float(y.min()), float(y.max())


# --------------------------------------------------------------------------- #
# Gridding                                                                     #
# --------------------------------------------------------------------------- #
def _tidy_grid(x_edges, y_edges, **value_grids):
    """Flatten one or more ``(bins, bins)`` histogram grids into a tidy DataFrame.

    Every cell is included (even empty ones), so every pathogen's grid for a method has the same
    ``bin_i``/``bin_j`` index as that method's background grid.
    """
    x_centers = (x_edges[:-1] + x_edges[1:]) / 2
    y_centers = (y_edges[:-1] + y_edges[1:]) / 2
    bin_i, bin_j = np.meshgrid(np.arange(len(x_centers)), np.arange(len(y_centers)), indexing="ij")
    xc, yc = np.meshgrid(x_centers, y_centers, indexing="ij")
    out = {"bin_i": bin_i.ravel(), "bin_j": bin_j.ravel(),
           "x_center": xc.ravel(), "y_center": yc.ravel()}
    for name, grid in value_grids.items():
        out[name] = grid.ravel()
    return pd.DataFrame(out)


def background_density(proj, method, bins, extent):
    """Full-library compound count per grid cell — the shared grey base layer for one method."""
    x, y = proj[f"{method}_x"].to_numpy(), proj[f"{method}_y"].to_numpy()
    valid = np.isfinite(x) & np.isfinite(y)
    xmin, xmax, ymin, ymax = extent
    count, x_edges, y_edges = np.histogram2d(
        x[valid], y[valid], bins=bins, range=[[xmin, xmax], [ymin, ymax]])
    return _tidy_grid(x_edges, y_edges, n_compounds=count)


def _score_column(pred_path):
    """The column holding this pathogen file's headline score.

    Almost every prediction file has a synthesized ``consensus_score`` (Isaura's aggregate over
    several ChEMBL sub-models). Three pathogens — campylobacter (eos7iak), hpylori (eos9eyo),
    ngonorrhoeae (eos5qya) — have only ONE ChEMBL sub-model, so there was never more than one
    column to take a consensus OVER and it is missing; for those, the file's single score column
    (e.g. ``chembl_dose_response_0``) IS the headline score and is used directly.
    """
    header = pd.read_csv(pred_path, nrows=0).columns
    if "consensus_score" in header:
        return "consensus_score"
    candidates = [c for c in header if c not in JOIN_COLS]
    if len(candidates) != 1:
        raise ValueError(
            f"{pred_path}: no consensus_score and {len(candidates)} candidate score columns "
            f"{candidates} — ambiguous, needs review.")
    print(f"    [note] {os.path.basename(pred_path)} has no consensus_score — "
          f"using its only score column '{candidates[0]}' instead")
    return candidates[0]


def _read_pathogen_scores(pred_path, chunk_size):
    """Chunked read of one pathogen prediction file's ``key`` + headline score column ONLY.

    ``usecols`` already keeps this to 2 of the file's ~14 columns; ``chunksize`` adds the same
    row-batching safety margin used for the much wider files in ``eval_correlations``. The score
    column (see :func:`_score_column`) is renamed to ``consensus_score`` so every pathogen's grid
    is built the same way regardless of which column it came from.
    """
    score_col = _score_column(pred_path)
    parts = [chunk.rename(columns={score_col: "consensus_score"}) for chunk in
              pd.read_csv(pred_path, usecols=["key", score_col], chunksize=chunk_size)]
    return pd.concat(parts, ignore_index=True)


def pathogen_top_n(proj, pred_path, n, chunk_size=CORR_CHUNK_SIZE):
    """This pathogen's top-``n`` highest ``consensus_score`` compounds, with EVERY projection
    method's coordinates attached.

    The score ranking does not depend on projection method, so this is computed ONCE per pathogen
    (not once per method): reads only ``key`` + the headline score column from ``pred_path``,
    reduces immediately to its ``n`` highest rows, then looks those ``n`` keys up against every
    ``{method}_x/_y`` column pair in one merge — never a method x pathogen size table.
    """
    scores = _read_pathogen_scores(pred_path, chunk_size).nlargest(n, "consensus_score")
    coord_cols = [f"{m}_{ax}" for m in PROJECTION_METHODS for ax in ("x", "y")]
    return scores.merge(proj[["key"] + coord_cols], on="key", how="left").reset_index(drop=True)


# --------------------------------------------------------------------------- #
# AUROC-matrix-aligned UMAP variant                                           #
# --------------------------------------------------------------------------- #
def auroc_matched_top_n_per_organism(scores_parquet, row_order_csv, proj, top_n=PROJECTION_TOP_N,
                                method="umap"):
    """Each organism's top-``n`` compounds by step 10's OWN aggregate score, in step 10's OWN row
    order, with ``method``'s projection coordinates attached.

    Reads ``10_organism_scores.parquet`` (the 1.35M x 15 rank-pct-averaged, cross-endpoint-merged
    aggregate score matrix step 10 writes, plus its ``key`` column) and
    ``10_row_order_comparison.csv`` (step 10's phylogeny-within-class row order and organism_class)
    directly, rather than re-deriving the merge — this guarantees the SAME top-``n`` compounds and
    the SAME row order the AUROC matrix draws, at the cost of requiring ``10_auroc_matrix.py`` to
    have already been run. This is the one panel of step 11 that depends on another step; the caller
    (``11_reference_library_projection.py``) treats it as optional and skips it, not the whole
    script, when step 10 hasn't run yet — see that script for the guard (in the current numbering,
    step 10 always runs first, so this is a fallback rather than the expected path).

    Unlike :func:`pathogen_top_n` (per-pathogen raw ``consensus_score``, used by the PCA/t-SNE/TMAP
    panels), this ranks each organism's ``organism_scores`` aggregate — a genuinely different set of
    compounds for any organism with 2+ merged endpoints (see ``eval_auroc_matrix.organism_scores``).
    """
    order = pd.read_csv(row_order_csv).sort_values("phylo_position")
    # Step 10 saves this parquet with its compound "key" as the DataFrame INDEX (not a column) —
    # reset it so every organism column can be selected alongside "key" by name below.
    scores = pd.read_parquet(scores_parquet).reset_index()
    coord_cols = [f"{method}_{ax}" for ax in ("x", "y")]

    frames = []
    for row in order.itertuples():
        top = scores.nlargest(top_n, row.organism)[[row.organism, "key"]].rename(
            columns={row.organism: "score"})
        top.insert(0, "phylo_position", row.phylo_position)
        top.insert(0, "organism_class", row.organism_class)
        top.insert(0, "organism", row.organism)
        frames.append(top)
    table = pd.concat(frames, ignore_index=True)
    return table.merge(proj[["key"] + coord_cols], on="key", how="left").reset_index(drop=True)


def run_auroc_matched_umap(proj, scores_parquet, row_order_csv, top_n=PROJECTION_TOP_N, method="umap"):
    """The AUROC-matrix-aligned UMAP top-N table, for :class:`plots_projection.AurocMatchedUmapGridPlot`.

    Returned in memory only — no CSV is written. Nothing reads this table back from disk (the figure
    is built from this same in-memory return value), and it used to be written as a CSV purely as an
    export; removed 2026-09-02 as unused. Contrast with ``11_top{top_n}_per_pathogen.csv``
    (:func:`run_all`), which the PCA/t-SNE/TMAP panels and step 12's Fisher enrichment both read back
    and so must stay a real file.
    """
    table = auroc_matched_top_n_per_organism(scores_parquet, row_order_csv, proj, top_n=top_n,
                                        method=method)
    missing = int(table[f"{method}_x"].isna().sum())
    note = f", {missing} without coordinates" if missing else ""
    print(f"  [auroc-matched-umap] {len(table)} rows, {table['organism'].nunique()} organisms{note}")
    return table


# --------------------------------------------------------------------------- #
# Orchestrator                                                                #
# --------------------------------------------------------------------------- #
def run_all(projection_file, pred_dir, pathogens_csv, output_dir, bins=PROJECTION_BINS,
            top_n=PROJECTION_TOP_N, chunk_size=CORR_CHUNK_SIZE):
    """Step 11 orchestrator: one background density grid per method, plus each pathogen's
    top-``top_n`` compounds (with every method's coordinates) computed once.

    Writes ``11_{method}_background.csv`` (one row per grid cell: ``bin_i, bin_j, x_center,
    y_center, n_compounds``) for each of :data:`default.PROJECTION_METHODS`, and one
    ``11_top{top_n}_per_pathogen.csv`` (one row per (pathogen, compound): ``pathogen_code, key,
    consensus_score`` + every method's ``{method}_x/_y``).

    Returns the loaded projection table so :func:`run_coadd` can reuse it: it is ~1.35M rows x 9
    columns and both families need it, so the caller threads one copy through rather than reading
    the same file twice per run.
    """
    pathogens = pd.read_csv(pathogens_csv)
    proj = load_projection(projection_file)
    print(f"[projection] loaded {os.path.basename(projection_file)} for {len(proj)} compounds")

    for method in PROJECTION_METHODS:
        extent = method_extent(proj, method)
        background = background_density(proj, method, bins, extent)
        background.to_csv(os.path.join(output_dir, f"11_{method}_background.csv"), index=False)
        print(f"  [{method}] background grid ({bins}x{bins}) -> 11_{method}_background.csv")

    files = latest_version_files(pred_dir)
    tops = []
    for _, row in pathogens.iterrows():
        code, eosid = row["code"], row["eosid"]
        entry = files.get(eosid)
        if entry is None:
            print(f"  [skip] {code} ({eosid}): no prediction file in {pred_dir}")
            continue
        _, pred_path = entry
        top = pathogen_top_n(proj, pred_path, top_n, chunk_size)
        top.insert(0, "pathogen_code", code)
        tops.append(top)
        print(f"  {code}: top {len(top)}/{top_n} "
              f"(score {top['consensus_score'].min():.3f}-{top['consensus_score'].max():.3f})")
    if tops:
        out_path = os.path.join(output_dir, f"11_top{top_n}_per_pathogen.csv")
        pd.concat(tops, ignore_index=True).to_csv(out_path, index=False)
        print(f"  -> {os.path.basename(out_path)} ({len(tops)} pathogens)")
    return proj


# --------------------------------------------------------------------------- #
# eos3dys (CoAdd) variant — same layout, an independent predictor                #
# --------------------------------------------------------------------------- #
def coadd_endpoints(selection_csv, model_id=COADD_MODEL_ID):
    """The model's bioactivity endpoints from ``config/08_endpoint_selection.csv``, in config order.

    Filters on ``assay_type == "bioactivity"``, which drops eos3dys's two ``Homo sapiens`` columns
    (``cytotoxicity_ic50``, ``hemolitic_activity``) and keeps its 20 pathogen endpoints. For eos3dys
    this is currently the same set as ``selected == "Yes"``, but ``assay_type`` is the filter that
    states the actual decision — pathogen endpoint, not human toxicity — and those two columns ARE
    selected in ``config/cytotoxicity_models.csv``, where step 13 picks them up instead.

    Adds two derived columns used by the figure: ``organism_token`` (the ``column_name`` prefix
    before the first underscore, e.g. ``paeruginosa``) and ``label`` (the remainder, e.g.
    ``PAO397_mic_25``), which is what a ~30 mm panel has room to print. Config order is preserved
    because it determines the shade ramp: it happens to place wild-type strains before sensitised
    ones, so the darkest shade in a panel is the wild-type reference strain.

    Raises if ``organism_token`` -> ``organism`` is not 1:1 — that would be a config error which
    would otherwise silently mislabel or split a panel.
    """
    config = pd.read_csv(selection_csv)
    sel = config[(config["model_id"] == model_id)
                 & (config["assay_type"] == "bioactivity")].copy()
    if not len(sel):
        raise ValueError(
            f"{selection_csv}: no bioactivity endpoints for {model_id} — nothing to plot.")
    sel["organism_token"] = [c.split("_", 1)[0] for c in sel["column_name"]]
    sel["label"] = [c.split("_", 1)[1] if "_" in c else c for c in sel["column_name"]]
    clashes = {t: sorted(g["organism"].unique())
               for t, g in sel.groupby("organism_token") if g["organism"].nunique() > 1}
    if clashes:
        raise ValueError(
            f"{selection_csv}: column-name prefix maps to more than one organism for {model_id}: "
            f"{clashes} — fix the spelling in the config, do not guess here.")
    return sel[["model_id", "column_name", "organism", "organism_token",
                "label"]].reset_index(drop=True)


def coadd_endpoint_top_n(pred_path, endpoint_cols, n=PROJECTION_TOP_N, chunk_size=CORR_CHUNK_SIZE):
    """Each endpoint's ``n`` highest-scoring compounds, as ``{column_name: DataFrame[key, score]}``.

    eos3dys has no ``consensus_score`` to rank an organism by (see :func:`_score_column` for the
    same situation among the pathogen models), so every endpoint is ranked on its own column. All 20
    are ranked in ONE pass: streams ``key`` + the 20 endpoint columns and re-reduces each endpoint's
    running top-N after every chunk, so memory stays flat at one chunk plus 20 x n rows. ``input``
    is never read — that SMILES column is most of the file's 671 MB.

    This mirrors :func:`eval_tox_projection.endpoint_top_n` and is deliberately NOT imported from
    it: that module already imports :func:`load_projection` from here, so importing back would be a
    circular import.
    """
    tops = {c: None for c in endpoint_cols}
    n_rows = 0
    for chunk in pd.read_csv(pred_path, usecols=["key"] + list(endpoint_cols),
                             chunksize=chunk_size):
        n_rows += len(chunk)
        for c in endpoint_cols:
            part = chunk[["key", c]]
            tops[c] = part if tops[c] is None else pd.concat([tops[c], part], ignore_index=True)
            tops[c] = tops[c].nlargest(n, c)
    print(f"[coadd-projection] ranked {n_rows} compounds across {len(endpoint_cols)} endpoints")
    return {c: df.rename(columns={c: "score"}).reset_index(drop=True) for c, df in tops.items()}


def run_coadd(proj, pred_path, selection_csv, output_dir, method=COADD_PROJECTION_METHOD,
              top_n=PROJECTION_TOP_N, chunk_size=CORR_CHUNK_SIZE):
    """Step 11's second family: :data:`default.COADD_MODEL_ID`'s top-``top_n`` per endpoint.

    Writes ``11_coadd_top{top_n}_per_endpoint.csv`` (one row per (endpoint, compound):
    ``organism, column_name, key, score, {method}_x, {method}_y``) and returns the endpoint frame
    from :func:`coadd_endpoints` for the figure to group into panels.

    ``proj`` is :func:`run_all`'s already-loaded projection table, and the background grids it wrote
    are reused unchanged — the point of this family is that it sits on a background identical to the
    pathogen panels', so the two are directly comparable.
    """
    endpoints = coadd_endpoints(selection_csv)
    cols = endpoints["column_name"].tolist()
    print(f"[coadd-projection] {len(cols)} bioactivity endpoints across "
          f"{endpoints['organism'].nunique()} organisms from {os.path.basename(selection_csv)}")

    tops = coadd_endpoint_top_n(pred_path, cols, n=top_n, chunk_size=chunk_size)

    coord_cols = [f"{method}_{ax}" for ax in ("x", "y")]
    frames = []
    for r in endpoints.itertuples():
        merged = tops[r.column_name].merge(proj[["key"] + coord_cols], on="key", how="left")
        merged.insert(0, "column_name", r.column_name)
        merged.insert(0, "organism", r.organism)
        frames.append(merged)
        # Report where each endpoint's rank cutoff landed. No threshold is applied, but a top-N
        # whose scores are all near-saturated means the cutoff is separating far less than it does
        # for a sharply-peaked endpoint — and the 20 panels are only comparable if this is checked.
        missing = int(merged[f"{method}_x"].isna().sum())
        note = f", {missing} without coordinates" if missing else ""
        print(f"  {r.column_name}: top {len(merged)}/{top_n} "
              f"(score {merged['score'].min():.3f}-{merged['score'].max():.3f}){note}")

    table = pd.concat(frames, ignore_index=True)
    out_path = os.path.join(output_dir, f"11_coadd_top{top_n}_per_endpoint.csv")
    table.to_csv(out_path, index=False)
    print(f"  -> {os.path.basename(out_path)} ({len(table)} rows)")
    return endpoints
