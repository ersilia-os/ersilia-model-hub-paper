import gzip
import os
import shutil
import urllib.request

from isaura.manage import IsauraReader

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_ISAURA_DIR = os.path.join(ROOT, "data", "raw", "isaura")
DEFAULT_COMPOUNDS_DIR = os.path.join(ROOT, "data", "raw", "compounds")

_EMH_PREFIX = "emh_paper"


def download_from_isaura(
    model_id: str,
    model_version: str,
    input_csv: str,
    bucket: str = "isaura-public",
    output_csv: str = None,
    approximate: bool = False,
    compress: bool = False,
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

    Returns
    -------
    str
        Absolute path to the written output file.
    """
    if output_csv is None:
        os.makedirs(DEFAULT_ISAURA_DIR, exist_ok=True)
        filename = f"{_EMH_PREFIX}_{model_id}_{model_version}.csv"
        output_csv = os.path.join(DEFAULT_ISAURA_DIR, filename)

    reader = IsauraReader(
        model_id=model_id,
        model_version=model_version,
        bucket=bucket,
        input_csv=input_csv,
        approximate=approximate,
    )
    reader.read(output_csv=output_csv)

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

    urllib.request.urlretrieve(_REFERENCE_LIBRARY_URL, output_csv)
    return output_csv
