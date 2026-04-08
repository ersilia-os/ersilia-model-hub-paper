"""Download precalculated Isaura outputs for the Ersilia reference library.

Usage
-----
    python scripts/download_reference.py --model eos42ez
    python scripts/download_reference.py --model eos42ez --version v1
    python scripts/download_reference.py --model eos42ez --force
"""

import argparse
import os
import sys

# Allow imports from src/ regardless of working directory
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from src.utils import DEFAULT_COMPOUNDS_DIR, download_from_isaura, download_reference_library

REFERENCE_LIBRARY_CSV = os.path.join(DEFAULT_COMPOUNDS_DIR, "reference_library_smiles.csv")


def main():
    parser = argparse.ArgumentParser(
        description="Download Isaura precalculations for the Ersilia reference library."
    )
    parser.add_argument("--model", required=True, help="Ersilia model ID (e.g. eos42ez)")
    parser.add_argument("--version", default="v1", help="Model version (default: v1)")
    parser.add_argument("--bucket", default="isaura-public", help="Isaura bucket (default: isaura-public)")
    parser.add_argument("--force", action="store_true", help="Re-download reference library even if cached")
    parser.add_argument("--local", action="store_true", help="Use local MinIO instead of cloud (default: cloud)")
    args = parser.parse_args()

    # Step 1: reference library
    if os.path.exists(REFERENCE_LIBRARY_CSV) and not args.force:
        print(f"Reference library already cached: {REFERENCE_LIBRARY_CSV}")
    else:
        print("Downloading reference library...")
        download_reference_library()
        print(f"  -> {REFERENCE_LIBRARY_CSV}")

    # Step 2: Isaura precalculations
    print(f"Downloading Isaura outputs for {args.model} {args.version}...")
    output_csv = download_from_isaura(
        model_id=args.model,
        model_version=args.version,
        input_csv=REFERENCE_LIBRARY_CSV,
        bucket=args.bucket,
        cloud=not args.local,
    )
    print(f"  -> {output_csv}")


if __name__ == "__main__":
    main()
