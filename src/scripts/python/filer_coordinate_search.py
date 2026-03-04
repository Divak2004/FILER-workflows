#!/usr/bin/env python3
"""
Recipe 3 — Search for overlapping tracks by genomic coordinates.

Queries the FILER discovery endpoint to find all tracks in a genome build
that overlap a specific region, with optional metadata filtering.

Usage:
  python src/scripts/python/filer_coordinate_search.py \
    --region "chr1:100000-200000" \
    --genome-build hg38 \
    --filter-string '.metadata."Data source" == "ENCODE"' \
    --out search_results.tsv
"""

import argparse
import os
import sys
import requests
import pandas as pd

# Update this to your actual endpoint filename
ENDPOINT = "https://tf.lisanwanglab.org/FILER2/get_overlapping_tracks_by_coord.php"

def coordinate_search(region: str, build: str, filters: dict) -> pd.DataFrame:
    payload = {
        "region": region,
        "genomeBuild": build,
        "outputFormat": "json",
        **filters
    }

    try:
        # Using GET to handle potentially complex filter strings
        r = requests.get(ENDPOINT, params=payload, timeout=300)
        r.raise_for_status()
        
        # The script may return an error string instead of JSON
        if r.text.startswith("ERROR"):
            print(f"API {r.text.strip()}", file=sys.stderr)
            sys.exit(1)

        data = r.json()
        return pd.DataFrame(data) if data else pd.DataFrame()

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

def main():
    p = argparse.ArgumentParser(description="Global search for FILER tracks by coordinates.")
    p.add_argument("--region", required=True, help='e.g., "chr1:100-200"')
    p.add_argument("--genome-build", required=True, choices=["hg19", "hg38"])
    p.add_argument("--filter-string", default=".", help='jq-style filter string')
    p.add_argument("--full-metadata", action="store_true", help="Return all metadata fields")
    p.add_argument("--count-only", type=int, choices=[0, 1], default=1, 
                   help="1: Return IDs and hit counts only. 0: Return detailed info.")
    p.add_argument("--out", default="search_results.tsv", help="Output TSV path")
    
    args = p.parse_args()

    filters = {
        "filterString": args.filter_string,
        "fullMetadata": 1 if (args.full_metadata or args.filter_string != ".") else 0,
        "countOnly": args.count_only
    }

    print(f"[recipe03] Searching {args.genome_build} for overlaps in {args.region}...", file=sys.stderr)

    df = coordinate_search(args.region, args.genome_build, filters)

    if df.empty:
        print("[recipe03] No overlapping tracks found.", file=sys.stderr)
        sys.exit(0)

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    df.to_csv(args.out, sep="\t", index=False)
    print(f"[recipe03] Found {len(df)} tracks → {args.out}", file=sys.stderr)

if __name__ == "__main__":
    main()