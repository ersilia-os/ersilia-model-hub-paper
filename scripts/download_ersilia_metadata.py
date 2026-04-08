"""Download Ersilia model metadata from the public Airtable share.

Usage
-----
    python scripts/download_ersilia_metadata.py
    python scripts/download_ersilia_metadata.py --output path/to/file.csv
    python scripts/download_ersilia_metadata.py --force
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from src.utils import DEFAULT_METADATA_DIR, download_ersilia_metadata

DEFAULT_OUTPUT = os.path.join(DEFAULT_METADATA_DIR, "ersilia_metadata.csv")


def main():
    parser = argparse.ArgumentParser(description="Download Ersilia model metadata from Airtable.")
    parser.add_argument("--output", default=DEFAULT_OUTPUT, help="Output CSV path")
    parser.add_argument("--force", action="store_true", help="Re-download even if file already exists")
    args = parser.parse_args()

    if os.path.exists(args.output) and not args.force:
        print(f"Already cached: {args.output}")
        return

    print("Downloading Ersilia metadata from Airtable...")
    path = download_ersilia_metadata(output_csv=args.output)
    print(f"  -> {path}")


if __name__ == "__main__":
    main()
