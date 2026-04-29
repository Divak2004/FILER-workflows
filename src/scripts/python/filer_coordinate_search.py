#!/usr/bin/env python3
"""
Recipe 3 — Search for overlapping tracks by genomic coordinates.

Queries the FILER discovery endpoint to find all tracks in a genome build
that overlap a specific region, with optional metadata filtering.

Usage:
  # Convenience flags (same style as recipe 1) — combined with `and` into a jq filter
  python src/scripts/python/filer_coordinate_search.py \
    --region "chr1:100000-200000" \
    --genome-build hg38 \
    --assay "ATAC-seq" \
    --data-source "ENCODE" \
    --out search_results.tsv

  # Raw jq filter string (overrides convenience flags)
  python src/scripts/python/filer_coordinate_search.py \
    --region "chr1:100000-200000" \
    --genome-build hg38 \
    --filter-string '.data_source == "ENCODE" and .assay == "ATAC-seq"' \
    --out search_results.tsv
"""

import argparse
import os
import sys
import requests
import pandas as pd

# Update this to your actual endpoint filename
ENDPOINT = "https://tf.lisanwanglab.org/FILER2/get_overlapping_tracks_by_coord.php"

# Map user-friendly CLI names → metadata field names used in jq filter expressions
JQ_FIELD_MAP = {
    "assay":           ".assay",
    "cell_type":       ".cell_type",
    "tissue_category": ".tissue_category",
    "data_source":     ".data_source",
    "track_id":        ".identifier",
}


def build_filter_string(named: dict) -> str:
    """Combine named filters into a jq boolean expression joined with `and`."""
    clauses = [
        f'{JQ_FIELD_MAP[k]} == "{v}"'
        for k, v in named.items()
        if v is not None
    ]
    return " and ".join(clauses) if clauses else "."


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
        if not data:
            return pd.DataFrame()

        df = pd.DataFrame(data)
        df.insert(df.columns.get_loc("identifier") + 1, "queryRegion", region)
        return df

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

def main():
    p = argparse.ArgumentParser(description="Global search for FILER tracks by coordinates.")
    p.add_argument("--region", required=True, help='e.g., "chr1:100-200"')
    p.add_argument("--genome-build", required=True, choices=["hg19", "hg38"])
    p.add_argument("--assay",           help='e.g. "ATAC-seq"')
    p.add_argument("--cell-type",       help='e.g. "CD14+ monocyte"')
    p.add_argument("--tissue-category", help='e.g. "Blood"')
    p.add_argument("--data-source",     help='e.g. "ENCODE"')
    p.add_argument("--track-id",        help="Specific track identifier")
    p.add_argument("--filter-string", default=None,
                   help='Raw jq-style filter string. Overrides the convenience flags above.')
    p.add_argument("--full-metadata", action="store_true", help="Return all metadata fields")
    p.add_argument("--count-only", type=int, choices=[0, 1], default=1,
                   help="1: Return IDs and hit counts only. 0: Return detailed info.")
    p.add_argument("--out", default="search_results.tsv", help="Output TSV path")

    args = p.parse_args()

    if args.filter_string is not None:
        filter_string = args.filter_string
    else:
        filter_string = build_filter_string({
            "assay":           args.assay,
            "cell_type":       args.cell_type,
            "tissue_category": args.tissue_category,
            "data_source":     args.data_source,
            "track_id":        args.track_id,
        })

    has_filter = filter_string != "."
    filters = {
        "filterString": filter_string,
        "fullMetadata": 1 if (args.full_metadata or has_filter) else 0,
        "countOnly": args.count_only
    }

    print(f"[recipe03] Searching {args.genome_build} for overlaps in {args.region}...", file=sys.stderr)
    if has_filter:
        print(f"[recipe03] filterString: {filter_string}", file=sys.stderr)

    df = coordinate_search(args.region, args.genome_build, filters)

    if df.empty:
        print("[recipe03] No overlapping tracks found.", file=sys.stderr)
        sys.exit(0)

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    df.to_csv(args.out, sep="\t", index=False)
    print(f"[recipe03] Found {len(df)} tracks → {args.out}", file=sys.stderr)

if __name__ == "__main__":
    main()
