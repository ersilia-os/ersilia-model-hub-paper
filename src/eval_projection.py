"""Step 09 analysis engine — reference-library chemical-space projection, coloured by pathogen
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
"""

import os

import numpy as np
import pandas as pd

from default import CORR_CHUNK_SIZE, PROJECTION_BINS, PROJECTION_METHODS, PROJECTION_TOP_N
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
# Orchestrator                                                                #
# --------------------------------------------------------------------------- #
def run_all(projection_file, pred_dir, pathogens_csv, output_dir, bins=PROJECTION_BINS,
            top_n=PROJECTION_TOP_N, chunk_size=CORR_CHUNK_SIZE):
    """Step 09 orchestrator: one background density grid per method, plus each pathogen's
    top-``top_n`` compounds (with every method's coordinates) computed once.

    Writes ``09_{method}_background.csv`` (one row per grid cell: ``bin_i, bin_j, x_center,
    y_center, n_compounds``) for each of :data:`default.PROJECTION_METHODS`, and one
    ``09_top{top_n}_per_pathogen.csv`` (one row per (pathogen, compound): ``pathogen_code, key,
    consensus_score`` + every method's ``{method}_x/_y``).
    """
    pathogens = pd.read_csv(pathogens_csv)
    proj = load_projection(projection_file)
    print(f"[projection] loaded {os.path.basename(projection_file)} for {len(proj)} compounds")

    for method in PROJECTION_METHODS:
        extent = method_extent(proj, method)
        background = background_density(proj, method, bins, extent)
        background.to_csv(os.path.join(output_dir, f"09_{method}_background.csv"), index=False)
        print(f"  [{method}] background grid ({bins}x{bins}) -> 09_{method}_background.csv")

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
        out_path = os.path.join(output_dir, f"09_top{top_n}_per_pathogen.csv")
        pd.concat(tops, ignore_index=True).to_csv(out_path, index=False)
        print(f"  -> {os.path.basename(out_path)} ({len(tops)} pathogens)")
