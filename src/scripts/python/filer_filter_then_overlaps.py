#!/usr/bin/env python3
"""
Workflow — Metadata Filter → Coordinate Search → Get Overlaps → Joined Table

Runs all three FILER recipes in sequence and joins their outputs into a single
row-per-interval table.

Steps:
  1. Recipe 1 — query metadata endpoint to find tracks matching filters
  2. Recipe 3 — query coordinate endpoint to find tracks overlapping a region
  3. Intersect results from steps 1 and 2 on identifier, rank by num_overlaps
     (union of all columns from both recipes is preserved)
  4. Recipe 2 — fetch actual overlapping intervals from the top N tracks
  5. Join interval table with track metadata into final output

Usage:
  python filer_filter_then_overlaps.py \\
    --genome-build hg38 \\
    --region "chr1:100000-200000" \\
    --assay "ATAC-seq" \\
    --tissue-category "Blood" \\
    --data-source "ENCODE" \\
    --top 100 \\
    --out output/04-workflow/results.tsv
"""

import argparse
import os
import sys
import requests
import pandas as pd

# ── Endpoints ────────────────────────────────────────────────────────────────
METADATA_ENDPOINT = "https://tf.lisanwanglab.org/FILER2/get_metadata.php"
COORD_ENDPOINT    = "https://tf.lisanwanglab.org/FILER2/get_overlapping_tracks_by_coord.php"
OVERLAPS_ENDPOINT = "https://tf.lisanwanglab.org/FILER2/get_overlaps.php"

HIT_SEP = "@@@"

# ── Recipe 1 — metadata search ───────────────────────────────────────────────
PARAM_MAP = {
    "assay":           "assayType",
    "cell_type":       "cellType",
    "tissue_category": "tissueCategory",
    "data_source":     "dataSource",
    "track_id":        "trackID",
}

STANDARD_COLS = [
    "identifier", "genome_build", "assay", "cell_type", "biosample_type",
    "tissue_category", "system_category", "life_stage", "data_source",
    "data_category", "classification", "output_type", "track_name",
    "processed_file_download_url", "tabix_file_url",
]


def recipe1_search_tracks(genome_build: str, filters: dict) -> pd.DataFrame:
    """Query FILER metadata endpoint. Returns one row per matching track."""
    params = {"genomeBuild": genome_build, "outputFormat": "json", **filters}
    print(f"[workflow] Step 1 — querying metadata endpoint...", file=sys.stderr)

    r = requests.get(METADATA_ENDPOINT, params=params, timeout=60)
    r.raise_for_status()
    data = r.json()

    if not data:
        return pd.DataFrame(columns=STANDARD_COLS)

    df = pd.DataFrame(data)
    present_std = [c for c in STANDARD_COLS if c in df.columns]
    extras = [c for c in df.columns if c not in STANDARD_COLS]
    return df[present_std + extras]


# ── Recipe 3 — coordinate search ─────────────────────────────────────────────

def recipe3_coordinate_search(region: str, genome_build: str,
                               filter_string: str = ".") -> pd.DataFrame:
    """
    Query FILER coordinate endpoint. Always uses fullMetadata=1 so that
    num_overlaps and all metadata fields (including wget_command, file_size,
    processed_file_md5, etc.) are present for the intersection step.
    """
    print(f"[workflow] Step 2 — querying coordinate endpoint for {region}...", file=sys.stderr)
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
        f"[workflow] Step 4 — fetching intervals: {len(track_ids)} tracks "
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


# ── Intersection + ranking ────────────────────────────────────────────────────

def intersect_and_rank(r1: pd.DataFrame, r3: pd.DataFrame, top_n: int) -> pd.DataFrame:
    """
    Inner-join Recipe 1 and Recipe 3 on identifier, taking the UNION of all
    columns from both. Where a column exists in both, Recipe 1's value is
    preferred and Recipe 3's value fills in any nulls. Rows are ranked by
    num_overlaps (from Recipe 3) and the top N are returned.
    """
    print(
        f"[workflow] Step 3 — intersecting {len(r1)} (R1) × {len(r3)} (R3) tracks...",
        file=sys.stderr,
    )

    merged = r1.merge(r3, on="identifier", how="inner", suffixes=("", "_r3"))

    # Coalesce duplicate columns: prefer R1 value, fall back to R3
    for col in r3.columns:
        r3_col = f"{col}_r3"
        if r3_col in merged.columns:
            merged[col] = merged[col].combine_first(merged[r3_col])
            merged.drop(columns=[r3_col], inplace=True)

    # Ensure num_overlaps is present (comes from R3)
    if "num_overlaps" not in merged.columns:
        print("[workflow]   Warning: num_overlaps not found after merge; ordering may be arbitrary.",
              file=sys.stderr)
    else:
        merged = merged.sort_values("num_overlaps", ascending=False)

    print(f"[workflow]   {len(merged)} tracks in common; selecting top {top_n}.", file=sys.stderr)
    return merged.head(top_n).reset_index(drop=True)


# ── Final join ────────────────────────────────────────────────────────────────

def join_results(r2: pd.DataFrame, top_tracks: pd.DataFrame) -> pd.DataFrame:
    """
    Left-join Recipe 2 interval rows with ranked track metadata.
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
        description="FILER end-to-end workflow: metadata → coordinate search → overlaps → joined table.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    # Shared
    p.add_argument("--genome-build", required=True, choices=["hg19", "hg38"])
    p.add_argument("--region",       required=True, help='e.g. "chr1:100000-200000"')
    p.add_argument("--top",          type=int, default=100,
                   help="Number of top tracks by num_overlaps to query in Recipe 2 (default: 100)")
    p.add_argument("--chunk-size",   type=int, default=250,
                   help="Batch size for Recipe 2 requests (default: 250)")
    p.add_argument("--out",          default="output/04-workflow/results.tsv",
                   help="Final output TSV (default: output/04-workflow/results.tsv)")

    # Recipe 1 filters
    r1 = p.add_argument_group("Recipe 1 filters (metadata search)")
    r1.add_argument("--assay")
    r1.add_argument("--cell-type")
    r1.add_argument("--tissue-category")
    r1.add_argument("--data-source")
    r1.add_argument("--track-id")
    r1.add_argument("--filter-string-r1", dest="filter_string_r1",
                    help="Raw jq filter for Recipe 1 (overrides named filters)")

    # Recipe 3 filters
    r3 = p.add_argument_group("Recipe 3 filters (coordinate search)")
    r3.add_argument("--filter-string-r3", dest="filter_string_r3", default=".",
                    help="Raw jq filter for Recipe 3 coordinate search (default: none)")

    args = p.parse_args()

    # ── Step 1: Recipe 1 ──────────────────────────────────────────────────────
    if args.filter_string_r1:
        r1_filters = {"filterString": args.filter_string_r1}
    else:
        r1_filters = {}
        for friendly, php_param in PARAM_MAP.items():
            val = getattr(args, friendly.replace("-", "_"), None)
            if val:
                r1_filters[php_param] = val

    df_r1 = recipe1_search_tracks(args.genome_build, r1_filters)
    if df_r1.empty:
        print("[workflow] Recipe 1 returned no tracks. Exiting.", file=sys.stderr)
        sys.exit(0)
    print(f"[workflow]   {len(df_r1)} tracks from Recipe 1.", file=sys.stderr)

    # ── Step 2: Recipe 3 ──────────────────────────────────────────────────────
    df_r3 = recipe3_coordinate_search(args.region, args.genome_build, args.filter_string_r3)
    if df_r3.empty:
        print("[workflow] Recipe 3 returned no overlapping tracks. Exiting.", file=sys.stderr)
        sys.exit(0)
    print(f"[workflow]   {len(df_r3)} tracks from Recipe 3.", file=sys.stderr)

    # ── Step 3: Intersect and rank (union of columns) ─────────────────────────
    top_tracks = intersect_and_rank(df_r1, df_r3, args.top)
    if top_tracks.empty:
        print("[workflow] No tracks in common between Recipe 1 and Recipe 3. Exiting.", file=sys.stderr)
        sys.exit(0)

    # ── Step 4: Recipe 2 ──────────────────────────────────────────────────────
    track_ids = top_tracks["identifier"].tolist()
    df_r2 = recipe2_get_overlaps(args.region, track_ids, chunk_size=args.chunk_size)
    if df_r2.empty:
        print("[workflow] Recipe 2 returned no intervals. Exiting.", file=sys.stderr)
        sys.exit(0)
    print(f"[workflow]   {len(df_r2)} intervals from Recipe 2.", file=sys.stderr)

    # ── Step 5: Join ──────────────────────────────────────────────────────────
    print(f"[workflow] Step 5 — joining interval and metadata tables...", file=sys.stderr)
    final = join_results(df_r2, top_tracks)

    out_dir = os.path.dirname(args.out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    final.to_csv(args.out, sep="\t", index=False)
    print(f"[workflow] Done — {len(final)} rows written to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()