#!/usr/bin/env python3
"""
Recipe 1 — Find FILER tracks by metadata filters.

Queries the FILER metadata endpoint and returns a table of matching tracks.

Usage:
  python src/scripts/python/filer_search_tracks.py \\
    --genome-build hg38 \\
    --assay "ATAC-seq" \\
    --cell-type "CD14+ monocyte" \\
    --data-source ENCODE \\
    --out tracks.tsv \\
    --json

  # Advanced: raw jq filter string
  python src/scripts/python/filer_search_tracks.py \\
    --genome-build hg38 \\
    --filter-string '.data_source == "ENCODE" and .assay == "ATAC-seq"' \\
    --out tracks.tsv
"""
import argparse
import sys
import requests
import pandas as pd

ENDPOINT = "https://tf.lisanwanglab.org/FILER2/get_metadata.php"

# Columns to include in output (in order); others are appended after
STANDARD_COLS = [
    "identifier",
    "genome_build",
    "assay",
    "cell_type",
    "biosample_type",
    "tissue_category",
    "system_category",
    "life_stage",
    "data_source",
    "data_category",
    "classification",
    "output_type",
    "track_name",
    "processed_file_download_url",
    "tabix_file_url",
]

# Map user-friendly CLI names → PHP parameter names
PARAM_MAP = {
    "assay":           "assayType",
    "cell_type":       "cellType",
    "tissue_category": "tissueCategory",
    "data_source":     "dataSource",
    "track_id":        "trackID",
}


def search_tracks(genome_build: str, filters: dict) -> pd.DataFrame:
    """
    Query FILER metadata endpoint and return results as a DataFrame.

    Parameters
    ----------
    genome_build : str
        Genome build, e.g. "hg38" or "hg19".
    filters : dict
        Keys are PHP parameter names (e.g. "assayType", "cellType") or
        the special key "filterString" for a raw jq expression.

    Returns
    -------
    pd.DataFrame
        One row per matching track. Key columns: identifier, genome_build,
        assay, cell_type, tissue_category, data_source, tabix_file_url.
    """
    params = {
        "genomeBuild":  genome_build,
        "outputFormat": "json",
        **filters,
    }

    r = requests.get(ENDPOINT, params=params, timeout=60)
    r.raise_for_status()

    data = r.json()
    if not data:
        return pd.DataFrame(columns=STANDARD_COLS)

    df = pd.DataFrame(data)

    # Reorder: standard columns first (where present), then any extras
    present_std = [c for c in STANDARD_COLS if c in df.columns]
    extras = [c for c in df.columns if c not in STANDARD_COLS]
    return df[present_std + extras]


def main():
    p = argparse.ArgumentParser(
        description="Search FILER tracks by metadata filters.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument(
        "--genome-build", default="hg38", choices=["hg19", "hg38"],
        help="Genome build (default: hg38)",
    )
    p.add_argument("--assay",           help='e.g. "ATAC-seq"')
    p.add_argument("--cell-type",       help='e.g. "CD14+ monocyte"')
    p.add_argument("--tissue-category", help='e.g. "blood"')
    p.add_argument("--data-source",     help='e.g. "ENCODE"')
    p.add_argument("--track-id",        help="Specific track identifier")
    p.add_argument(
        "--filter-string",
        help='Raw jq filter, e.g. \'.data_source == "ENCODE"\'. '
             "Overrides other filter flags.",
    )
    p.add_argument("--out", default="tracks.tsv", help="Output TSV path (default: tracks.tsv)")
    p.add_argument("--json", action="store_true", help="Also write <out>.json")
    args = p.parse_args()

    # Build filters dict using PHP param names
    filters = {}

    if args.filter_string:
        # Raw jq filter takes precedence over named params
        filters["filterString"] = args.filter_string
    else:
        named = {
            "assay":           args.assay,
            "cell_type":       args.cell_type,
            "tissue_category": args.tissue_category,
            "data_source":     args.data_source,
            "track_id":        args.track_id,
        }
        for friendly, val in named.items():
            if val is not None:
                filters[PARAM_MAP[friendly]] = val

    df = search_tracks(args.genome_build, filters)

    if df.empty:
        print("[recipe01] No tracks matched the given filters.", file=sys.stderr)
        sys.exit(0)

    df.to_csv(args.out, sep="\t", index=False)
    print(f"[recipe01] {len(df)} tracks → {args.out}", file=sys.stderr)

    if args.json:
        json_out = args.out.replace(".tsv", "") + ".json"
        df.to_json(json_out, orient="records", indent=2)
        print(f"[recipe01] Wrote {json_out}", file=sys.stderr)


if __name__ == "__main__":
    main()