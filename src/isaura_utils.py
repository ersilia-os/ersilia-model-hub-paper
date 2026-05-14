import gzip
import shutil
import os
import sys

from isaura.const import (
    MINIO_CLOUD_AK,
    MINIO_CLOUD_SK,
    MINIO_ENDPOINT,
    MINIO_ENDPOINT_CLOUD,
    MINIO_LOCAL_AK,
    MINIO_LOCAL_SK,
)
from isaura.manage import IsauraReader

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root))

def download_from_isaura(
    model_id: str,
    model_version: str,
    input_csv: str,
    output_path: str,
    bucket: str = "isaura-public",
    approximate: bool = False,
    compress: bool = False,
    cloud: bool = True,
) -> str:
    """Download precalculated model outputs from Isaura into data/raw/isaura/.

    Output files follow the eosframes naming convention
    prefix: <model_id>_<version>.csv (e.g. emh_paper_eos8a4x_v1.csv).

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
        Destination path for the output CSV.(.gz if compress=True).
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
    filename = f"{model_id}_{model_version}.csv"
    output_csv = os.path.join(output_path, filename)

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
