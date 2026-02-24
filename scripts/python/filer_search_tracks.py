#!/usr/bin/env python3
"""
Recipe 1 — Find FILER tracks by metadata filters.

Usage:
  python filer_search_tracks.py \
    --genome-build hg38 \
    --assay "ATAC-seq" \
    --cell-type "CD14+ monocyte" \
    --data-source ENCODE \
    --out tracks.tsv
"""
import argparse
import sys
import requests
import pandas as pd

BASE_URL = "https://tf.lisanwanglab.org/FILER2"

REQUIRED_COLS = [
    "track_id", "genome_build", "assay",
    "cell_type", "data_source", "bed_gz_url",
]


def search_tracks(genome_build: str, filters: dict) -> pd.DataFrame:
    params = {"genome_build": genome_build, **filters}
    r = requests.get(f"{BASE_URL}/api/tracks", params=params, timeout=60)
    r.raise_for_status()
    data = r.json()

    # Adjust key to match actual API response structure
    records = data.get("results") or data.get("tracks") or data
    df = pd.DataFrame(records)

    # Reorder: required columns first, then any extras
    available_required = [c for c in REQUIRED_COLS if c in df.columns]
    extras = [c for c in df.columns if c not in REQUIRED_COLS]
    return df[available_required + extras]


def main():
    p = argparse.ArgumentParser(description="Search FILER tracks by metadata filters.")
    p.add_argument("--genome-build", default="hg38", choices=["hg19", "hg38"])
    p.add_argument("--assay", help='e.g. "ATAC-seq"')
    p.add_argument("--cell-type", help='e.g. "CD14+ monocyte"')
    p.add_argument("--biosample-type", help='e.g. "primary cell"')
    p.add_argument("--life-stage", help='e.g. "adult"')
    p.add_argument("--data-source", help='e.g. "ENCODE"')
    p.add_argument("--out", default="tracks.tsv", help="Output TSV path")
    p.add_argument("--json", action="store_true", help="Also write tracks.json")
    args = p.parse_args()

    filters = {k: v for k, v in {
        "assay":          args.assay,
        "cell_type":      args.cell_type,
        "biosample_type": args.biosample_type,
        "life_stage":     args.life_stage,
        "data_source":    args.data_source,
    }.items() if v is not None}

    df = search_tracks(args.genome_build, filters)
    df.to_csv(args.out, sep="\t", index=False)
    print(f"[recipe01] Wrote {len(df)} tracks → {args.out}", file=sys.stderr)

    if args.json:
        json_out = args.out.replace(".tsv", ".json")
        df.to_json(json_out, orient="records", indent=2)
        print(f"[recipe01] Wrote {json_out}", file=sys.stderr)


if __name__ == "__main__":
    main()