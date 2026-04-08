import gzip
import os
import shutil
import urllib.request

import pandas as pd

from isaura.const import (
    MINIO_CLOUD_AK,
    MINIO_CLOUD_SK,
    MINIO_ENDPOINT,
    MINIO_ENDPOINT_CLOUD,
    MINIO_LOCAL_AK,
    MINIO_LOCAL_SK,
)
from isaura.manage import IsauraReader

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_ISAURA_DIR = os.path.join(ROOT, "data", "raw", "isaura")
DEFAULT_COMPOUNDS_DIR = os.path.join(ROOT, "data", "raw", "compounds")
DEFAULT_METADATA_DIR = os.path.join(ROOT, "data", "raw")

_EMH_PREFIX = "emh_paper"


def download_from_isaura(
    model_id: str,
    model_version: str,
    input_csv: str,
    bucket: str = "isaura-public",
    output_csv: str = None,
    approximate: bool = False,
    compress: bool = False,
    cloud: bool = True,
) -> str:
    """Download precalculated model outputs from Isaura into data/raw/isaura/.

    Output files follow the eosframes naming convention with a repo-specific
    prefix: emh_paper_<model_id>_<version>.csv (e.g. emh_paper_eos8a4x_v1.csv).

    Parameters
    ----------
    model_id : str
        Ersilia model identifier (e.g. 'eos8a4x').
    model_version : str
        Model version string (e.g. 'v1').
    input_csv : str
        Path to CSV file containing the input compounds (SMILES or InChIKey).
    bucket : str
        Isaura bucket name. Defaults to 'isaura-public'.
    output_csv : str, optional
        Destination path for the output CSV. Defaults to
        data/raw/isaura/emh_paper_<model_id>_<version>.csv(.gz if compress=True).
    approximate : bool
        Whether to use approximate nearest-neighbour search for unseen inputs.
    compress : bool
        If True, gzip-compress the output file. The path returned will end in .gz.
    cloud : bool
        If True (default), connect to the cloud MinIO endpoint. If False, use
        the local MinIO instance.

    Returns
    -------
    str
        Absolute path to the written output file.
    """
    if output_csv is None:
        os.makedirs(DEFAULT_ISAURA_DIR, exist_ok=True)
        filename = f"{_EMH_PREFIX}_{model_id}_{model_version}.csv"
        output_csv = os.path.join(DEFAULT_ISAURA_DIR, filename)

    if cloud:
        endpoint, access_key, secrete = MINIO_ENDPOINT_CLOUD, MINIO_CLOUD_AK, MINIO_CLOUD_SK
    else:
        endpoint, access_key, secrete = MINIO_ENDPOINT, MINIO_LOCAL_AK, MINIO_LOCAL_SK

    reader = IsauraReader(
        model_id=model_id,
        model_version=model_version,
        bucket=bucket,
        input_csv=input_csv,
        approximate=approximate,
        endpoint=endpoint,
        access_key=access_key,
        secrete=secrete,
    )
    reader.read(output_csv=output_csv)

    if not os.path.exists(output_csv):
        raise FileNotFoundError(
            f"Isaura returned no results for {model_id} {model_version} — output file was not written."
        )

    if compress:
        compressed_path = output_csv + ".gz"
        with open(output_csv, "rb") as f_in, gzip.open(compressed_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
        os.remove(output_csv)
        return compressed_path

    return output_csv


_REFERENCE_LIBRARY_URL = (
    "https://raw.githubusercontent.com/ersilia-os/"
    "ersilia-model-hub-maintained-inputs/main/inputs/reference_library_smiles.csv"
)


def download_reference_library(
    output_csv: str = None,
) -> str:
    """Download the Ersilia reference library of compounds.

    Downloads the canonical reference set of SMILES from the
    ersilia-model-hub-maintained-inputs repository.

    Parameters
    ----------
    output_csv : str, optional
        Destination path for the CSV file. Defaults to
        data/raw/compounds/reference_library_smiles.csv.

    Returns
    -------
    str
        Absolute path to the written output CSV.
    """
    if output_csv is None:
        os.makedirs(DEFAULT_COMPOUNDS_DIR, exist_ok=True)
        output_csv = os.path.join(DEFAULT_COMPOUNDS_DIR, "reference_library_smiles.csv")

    tmp = output_csv + ".tmp"
    urllib.request.urlretrieve(_REFERENCE_LIBRARY_URL, tmp)
    df = pd.read_csv(tmp)
    df = df.rename(columns={"standardized_smiles": "input"})
    df.to_csv(output_csv, index=False)
    os.remove(tmp)
    return output_csv


_AIRTABLE_SHARE_URL = (
    "https://airtable.com/appR6ZwgLgG8RTdoU/shr7scXQV3UYqnM6Q/tblAfOWRbA7bI1VTB"
)
_AIRTABLE_VIEW_ID = "viwy0inGR1xv2Tpfe"


def download_ersilia_metadata(output_csv: str = None) -> str:
    """Download Ersilia model metadata from the public Airtable share.

    Fetches all records from the Ersilia model hub Airtable base at
    https://airtable.com/appR6ZwgLgG8RTdoU/shr7scXQV3UYqnM6Q/tblAfOWRbA7bI1VTB
    and saves them as a CSV file.

    Parameters
    ----------
    output_csv : str, optional
        Destination path for the CSV file. Defaults to
        data/raw/ersilia_metadata.csv.

    Returns
    -------
    str
        Absolute path to the written CSV file.
    """
    import io
    import json
    import re

    import requests

    if output_csv is None:
        os.makedirs(DEFAULT_METADATA_DIR, exist_ok=True)
        output_csv = os.path.join(DEFAULT_METADATA_DIR, "ersilia_metadata.csv")

    session = requests.Session()
    session.headers.update({
        "User-Agent": (
            "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 "
            "Chrome/120.0.0.0 Safari/537.36"
        )
    })

    # Step 1: load the share page to obtain the signed accessPolicy and CSRF token
    resp = session.get(_AIRTABLE_SHARE_URL, timeout=30)
    resp.raise_for_status()

    m = re.search(r"initData\s*=\s*(\{.+?\})\s*;\s*\n", resp.text, re.DOTALL)
    if not m:
        raise RuntimeError("Could not parse Airtable share page — page structure may have changed.")
    init = json.loads(m.group(1))
    csrf = init["csrfToken"]
    access_policy = init["accessPolicy"]

    # Step 2: download the CSV using the signed accessPolicy
    csv_url = f"https://airtable.com/v0.3/view/{_AIRTABLE_VIEW_ID}/downloadCsv"
    csv_resp = session.get(
        csv_url,
        params={"accessPolicy": access_policy},
        headers={
            "Accept": "*/*",
            "x-csrf-token": csrf,
            "x-airtable-application-id": "appR6ZwgLgG8RTdoU",
            "x-requested-with": "XMLHttpRequest",
            "x-time-zone": "UTC",
            "x-user-locale": "en-US",
            "Referer": _AIRTABLE_SHARE_URL,
        },
        timeout=60,
    )
    csv_resp.raise_for_status()

    df = pd.read_csv(io.StringIO(csv_resp.text))
    df.to_csv(output_csv, index=False)
    return output_csv
