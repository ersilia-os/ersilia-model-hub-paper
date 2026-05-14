import io
import json
import os
import re
import sys

import pandas as pd
import requests

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root))

from default import AIRTABLE_BASE_ID, AIRTABLE_SHARE_URL, AIRTABLE_VIEW_ID, REFERENCE_LIBRARY_URL

def download_reference_library() -> pd.DataFrame:
    """Download the canonical reference set of SMILES from the ersilia-model-hub-maintained-inputs repository."""
    resp = requests.get(REFERENCE_LIBRARY_URL, timeout=60)
    resp.raise_for_status()
    df = pd.read_csv(io.StringIO(resp.text))
    df = df.rename(columns={"standardized_smiles": "input"})
    n_before = len(df)
    df = df.drop_duplicates(subset=["input"])
    n_dropped = n_before - len(df)
    if n_dropped:
        print(f"  Removed {n_dropped} duplicate SMILES from reference library.")
    return df

def download_airtable_metadata() -> pd.DataFrame:
    """Fetch all model records from the Ersilia Model Hub public Airtable share and output a dataframe"""
    session = requests.Session()
    session.headers.update({
        "User-Agent": (
            "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 "
            "Chrome/120.0.0.0 Safari/537.36"
        )
    })
    resp = session.get(AIRTABLE_SHARE_URL, timeout=30)
    resp.raise_for_status()
    m = re.search(r"initData\s*=\s*(\{)", resp.text)
    if not m:
        raise RuntimeError(
            "Could not parse Airtable share page — page structure may have changed."
        )
    init, _ = json.JSONDecoder().raw_decode(resp.text, m.start(1))
    csrf_token = init["csrfToken"]
    access_policy = init["accessPolicy"]

    # Download the full table as CSV using the signed credentials
    csv_url = f"https://airtable.com/v0.3/view/{AIRTABLE_VIEW_ID}/downloadCsv"
    csv_resp = session.get(
        csv_url,
        params={"accessPolicy": access_policy},
        headers={
            "Accept": "*/*",
            "x-csrf-token": csrf_token,
            "x-airtable-application-id": AIRTABLE_BASE_ID,
            "x-requested-with": "XMLHttpRequest",
            "x-time-zone": "UTC",
            "x-user-locale": "en-US",
            "Referer": AIRTABLE_SHARE_URL,
        },
        timeout=60,
    )
    csv_resp.raise_for_status()

    df = pd.read_csv(io.StringIO(csv_resp.text))
    return df