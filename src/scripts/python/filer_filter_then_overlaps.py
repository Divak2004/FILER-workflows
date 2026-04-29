#!/usr/bin/env python3
"""
Workflow — Coordinate Search → Get Overlaps → Joined Table

Runs Recipe 3 (coordinate search with optional metadata filter) followed by
Recipe 2 (overlap extraction) and joins their outputs into a single
row-per-interval table.

Steps:
  1. Recipe 3 — query coordinate endpoint with the metadata filter to find
     tracks both matching biological criteria AND overlapping the region.
  2. Rank by num_overlaps and keep the top N tracks.
  3. Recipe 2 — fetch actual overlapping intervals from those top N tracks.
  4. Join interval table with track metadata into final output.

Usage:
  python filer_filter_then_overlaps.py \\
    --genome-build hg38 \\
    --region "chr1:100000-200000" \\
    --assay "ATAC-seq" \\
    --tissue-category "Blood" \\
    --data-source "ENCODE" \\
    --top 100 \\
    --out output/10-filter-then-overlaps/results.tsv
"""

import argparse
import os
import sys
import requests
import pandas as pd

# ── Endpoints ────────────────────────────────────────────────────────────────
COORD_ENDPOINT    = "https://tf.lisanwanglab.org/FILER2/get_overlapping_tracks_by_coord.php"
OVERLAPS_ENDPOINT = "https://tf.lisanwanglab.org/FILER2/get_overlaps.php"

HIT_SEP = "@@@"

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


# ── Recipe 3 — coordinate search ─────────────────────────────────────────────

def recipe3_coordinate_search(region: str, genome_build: str,
                               filter_string: str = ".") -> pd.DataFrame:
    """
    Query FILER coordinate endpoint. Always uses fullMetadata=1 so that
    num_overlaps and all metadata fields (including wget_command, file_size,
    processed_file_md5, etc.) are present for ranking and downstream use.
    """
    print(f"[workflow] Step 1 — querying coordinate endpoint for {region}...", file=sys.stderr)
    if filter_string != ".":
        print(f"[workflow]   filterString: {filter_string}", file=sys.stderr)
    params = {
        "region":       region,
        "genomeBuild":  genome_build,
        "outputFormat": "json",
        "filterString": filter_string,
        "fullMetadata": 1,
        "countOnly":    0,
    }

    r = requests.get(COORD_ENDPOINT, params=params, timeout=300)
    r.raise_for_status()

    if r.text.startswith("ERROR"):
        print(f"[workflow] Coordinate endpoint error: {r.text.strip()}", file=sys.stderr)
        sys.exit(1)

    data = r.json()
    return pd.DataFrame(data) if data else pd.DataFrame()


# ── Rank and select top N ─────────────────────────────────────────────────────

def rank_and_select(df: pd.DataFrame, top_n: int) -> pd.DataFrame:
    """Sort by num_overlaps descending and return the top N rows."""
    if "num_overlaps" not in df.columns:
        print("[workflow]   Warning: num_overlaps not found; ordering may be arbitrary.",
              file=sys.stderr)
        return df.head(top_n).reset_index(drop=True)

    ranked = df.sort_values("num_overlaps", ascending=False)
    print(f"[workflow] Step 2 — ranking {len(ranked)} tracks; selecting top {top_n}.", file=sys.stderr)
    return ranked.head(top_n).reset_index(drop=True)


# ── Recipe 2 — fetch overlapping intervals ────────────────────────────────────

def _fetch_chunk(region: str, track_ids_str: str) -> pd.DataFrame:
    """GET one batch of track IDs from the overlaps endpoint."""
    r = requests.get(
        OVERLAPS_ENDPOINT,
        params={"region": region, "trackIDs": track_ids_str},
        timeout=300,
    )
    r.raise_for_status()

    text = r.text.strip()
    if text.startswith("ERROR:"):
        print(f"[workflow] Overlaps endpoint error: {text}", file=sys.stderr)
        sys.exit(1)

    data = r.json()
    return pd.DataFrame(data) if isinstance(data, list) else pd.DataFrame()


def recipe2_get_overlaps(region: str, track_ids: list,
                          chunk_size: int = 250) -> pd.DataFrame:
    """
    Fetch overlapping intervals for a list of track IDs in batches.
    Returns one row per interval with identifier, queryRegion, and hitString.
    """
    batches = [track_ids[i:i+chunk_size] for i in range(0, len(track_ids), chunk_size)]
    print(
        f"[workflow] Step 3 — fetching intervals: {len(track_ids)} tracks "
        f"in {len(batches)} batch(es)...",
        file=sys.stderr,
    )

    rows: list[dict] = []
    for i, batch in enumerate(batches, 1):
        print(f"[workflow]   batch {i}/{len(batches)} ({len(batch)} IDs)...", file=sys.stderr)
        try:
            df_batch = _fetch_chunk(region, ",".join(batch))
        except Exception as e:
            print(f"[workflow]   Warning: batch {i} failed — {e}", file=sys.stderr)
            continue

        for _, row in df_batch.iterrows():
            identifier   = row.get("Identifier")
            query_region = row.get("queryRegion")
            features     = row.get("features", [])
            if not isinstance(features, list):
                continue
            for feature in features:
                if not isinstance(feature, dict):
                    continue
                hit_string = HIT_SEP.join("" if v is None else str(v) for v in feature.values())
                rows.append({
                    "identifier":  identifier,
                    "queryRegion": query_region,
                    "hitString":   hit_string,
                })

    return pd.DataFrame(rows)


# ── Final join ────────────────────────────────────────────────────────────────

def join_results(r2: pd.DataFrame, top_tracks: pd.DataFrame) -> pd.DataFrame:
    """
    Left-join Recipe 2 interval rows with the top-N track metadata.
    Result is one row per interval with full track metadata attached.
    Column order: identifier → queryRegion → hitString → metadata columns.
    """
    final = r2.merge(top_tracks, on="identifier", how="left")

    front_cols = ["identifier", "queryRegion", "hitString"]
    meta_cols  = [c for c in top_tracks.columns if c != "identifier"]
    ordered    = front_cols + [c for c in meta_cols if c in final.columns]
    extra      = [c for c in final.columns if c not in ordered]
    return final[ordered + extra]


# ── CLI ───────────────────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser(
        description="FILER workflow: coordinate search → overlaps → joined table.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    p.add_argument("--genome-build", required=True, choices=["hg19", "hg38"])
    p.add_argument("--region",       required=True, help='e.g. "chr1:100000-200000"')
    p.add_argument("--top",          type=int, default=100,
                   help="Number of top tracks by num_overlaps to query in Recipe 2 (default: 100)")
    p.add_argument("--chunk-size",   type=int, default=250,
                   help="Batch size for Recipe 2 requests (default: 250)")
    p.add_argument("--out",          default="output/10-filter-then-overlaps/results.tsv",
                   help="Final output TSV (default: output/10-filter-then-overlaps/results.tsv)")

    f = p.add_argument_group("Metadata filters (applied server-side via jq filterString)")
    f.add_argument("--assay",           help='e.g. "ATAC-seq"')
    f.add_argument("--cell-type",       help='e.g. "CD14+ monocyte"')
    f.add_argument("--tissue-category", help='e.g. "Blood"')
    f.add_argument("--data-source",     help='e.g. "ENCODE"')
    f.add_argument("--track-id",        help="Specific track identifier")
    f.add_argument("--filter-string",   default=None,
                   help="Raw jq filter string. Overrides the convenience flags above.")

    args = p.parse_args()

    # ── Build filter string ───────────────────────────────────────────────────
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

    # ── Step 1: Recipe 3 ──────────────────────────────────────────────────────
    df_r3 = recipe3_coordinate_search(args.region, args.genome_build, filter_string)
    if df_r3.empty:
        print("[workflow] Recipe 3 returned no overlapping tracks. Exiting.", file=sys.stderr)
        sys.exit(0)
    print(f"[workflow]   {len(df_r3)} tracks from Recipe 3.", file=sys.stderr)

    # ── Step 2: Rank and select top N ─────────────────────────────────────────
    top_tracks = rank_and_select(df_r3, args.top)

    # ── Step 3: Recipe 2 ──────────────────────────────────────────────────────
    track_ids = top_tracks["identifier"].tolist()
    df_r2 = recipe2_get_overlaps(args.region, track_ids, chunk_size=args.chunk_size)
    if df_r2.empty:
        print("[workflow] Recipe 2 returned no intervals. Exiting.", file=sys.stderr)
        sys.exit(0)
    print(f"[workflow]   {len(df_r2)} intervals from Recipe 2.", file=sys.stderr)

    # ── Step 4: Join ──────────────────────────────────────────────────────────
    print(f"[workflow] Step 4 — joining interval and metadata tables...", file=sys.stderr)
    final = join_results(df_r2, top_tracks)

    out_dir = os.path.dirname(args.out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    final.to_csv(args.out, sep="\t", index=False)
    print(f"[workflow] Done — {len(final)} rows written to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
