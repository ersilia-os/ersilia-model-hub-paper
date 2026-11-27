"""Load the frozen Airtable metadata snapshot used throughout the analysis."""

import os

import pandas as pd

from default import AIRTABLE_METADATA_FILE, AIRTABLE_STATUS_REMAP

root = os.path.dirname(os.path.abspath(__file__))
AIRTABLE_METADATA_PATH = os.path.join(root, "..", "data", "raw", AIRTABLE_METADATA_FILE)


def load_airtable_metadata(path: str = AIRTABLE_METADATA_PATH, verbose: bool = True) -> pd.DataFrame:
    """Read the Airtable snapshot and apply the project's Status remapping.

    The raw file in ``data/raw/`` is never modified; the remap (``AIRTABLE_STATUS_REMAP`` in
    ``src/default.py``) is applied only to the returned DataFrame, so every script sees the same
    statuses.

    Parameters
    ----------
    path : str
        CSV to read. Defaults to the snapshot named by ``AIRTABLE_METADATA_FILE``.
    verbose : bool
        If True, print how many models had their Status remapped.

    Returns
    -------
    pandas.DataFrame
        One row per model, with the remapped ``Status`` column.

    Raises
    ------
    FileNotFoundError
        If the snapshot has not been downloaded yet (run ``scripts/00_download_data.py``).
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"{path} not found — run scripts/00_download_data.py first.")
    # `keep_default_na=False` is load-bearing: `None` is a legitimate License value ("repo has no
    # LICENSE file"), and pandas' default NA set would turn it into NaN, indistinguishable from blank.
    df = pd.read_csv(path, encoding="utf-8-sig", keep_default_na=False)
    # The side effect of keep_default_na=False: a numeric column with any blank cell (e.g. runtime
    # or image size of a model never benchmarked) is read as text. Restore numbers wherever every
    # non-blank value parses, so blanks become NaN as they would with a default read.
    for col in df.columns[df.dtypes == object]:
        filled = df[col].str.strip() != ""
        parsed = pd.to_numeric(df.loc[filled, col], errors="coerce")
        if filled.any() and parsed.notna().all():
            df[col] = pd.to_numeric(df[col].where(filled), errors="coerce")
    for old, new in AIRTABLE_STATUS_REMAP.items():
        mask = df["Status"] == old
        if verbose and mask.any():
            print(f"  Status remap: {mask.sum()} '{old}' -> '{new}' "
                  f"({', '.join(sorted(df.loc[mask, 'Identifier']))})")
        df.loc[mask, "Status"] = new
    return df


def prediction_version(model_id: str, meta: pd.DataFrame | None = None) -> str:
    """Major version of a model's Isaura prediction file, from the snapshot's ``Release``.

    ``v2.0.0`` -> ``"v2"``, matching the ``{model_id}_{version}.csv`` names written by
    ``00_download_data.py``.

    Raises
    ------
    ValueError
        If the model is not in the snapshot or has no release.
    """
    meta = load_airtable_metadata(verbose=False) if meta is None else meta
    release = meta.loc[meta["Identifier"] == model_id, "Release"]
    if release.empty or not str(release.iloc[0]).strip():
        raise ValueError(f"{model_id} has no Release in {AIRTABLE_METADATA_FILE}.")
    return str(release.iloc[0]).strip().split(".")[0]
