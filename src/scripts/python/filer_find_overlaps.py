#!/usr/bin/env python3
"""
Recipe 2 — Extract genomic overlaps from specific tracks.

Queries the FILER overlaps endpoint to retrieve data intervals from 
specific tracks overlapping a given genomic region.

Usage:
  # Manual IDs
  python filer_find_overlaps.py --region "chr1:100000-200000" --track-ids "NG123,NG129"
  
  # From a TSV file (output of Recipe 1)
  python filer_find_overlaps.py --region "chr1:100000-200000" --file tracks.tsv --id-col identifier
"""

import argparse
import sys
import os
import requests
import pandas as pd
import json

ENDPOINT = "https://tf.lisanwanglab.org/FILER2/get_overlaps.php"

# The server only reads $_GET, so both params go in the query string.
# Keep chunks small enough to stay under the ~8KB URL limit.
# ~250 IDs × ~6 chars each + region + overhead ≈ safe margin.
CHUNK_SIZE = 250


def _get_chunk(region: str, track_ids_str: str) -> pd.DataFrame:
    """GET one batch — both params in the query string, which is all PHP reads."""
    try:
        r = requests.get(
            ENDPOINT,
            params={"region": region, "trackIDs": track_ids_str},
            timeout=300,
        )
        r.raise_for_status()
    except requests.exceptions.RequestException as e:
        body = r.text.strip() if "r" in locals() else str(e)
        print(f"[recipe02] Request failed: {body}", file=sys.stderr)
        sys.exit(1)

    text = r.text.strip()
    if text.startswith("ERROR:"):
        print(f"[recipe02] Server error: {text}", file=sys.stderr)
        sys.exit(1)

    try:
        data = r.json()
    except json.JSONDecodeError:
        print(f"[recipe02] Could not parse JSON response:\n{text[:500]}", file=sys.stderr)
        sys.exit(1)

    if isinstance(data, list):
        return pd.DataFrame(data)
    elif isinstance(data, dict):
        for key in ("data", "results", "records", "overlaps"):
            if key in data and isinstance(data[key], list):
                return pd.DataFrame(data[key])
        return pd.DataFrame([data])
    else:
        print(f"[recipe02] Unexpected response type: {type(data)}", file=sys.stderr)
        sys.exit(1)


def fetch_overlaps(region: str, ids: list, chunk_size: int = CHUNK_SIZE) -> pd.DataFrame:
    chunks = [ids[i : i + chunk_size] for i in range(0, len(ids), chunk_size)]
    n_chunks = len(chunks)

    print(
        f"[recipe02] Splitting {len(ids)} IDs into {n_chunks} batches of ≤{chunk_size}...",
        file=sys.stderr,
    )

    frames = []
    for i, chunk in enumerate(chunks, 1):
        print(f"[recipe02]   batch {i}/{n_chunks} ({len(chunk)} IDs)...", file=sys.stderr)
        df = _get_chunk(region, ",".join(map(str, chunk)))
        if not df.empty:
            frames.append(df)

    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def load_ids_from_file(path: str, id_col: str) -> list:
    if not os.path.exists(path):
        print(f"[recipe02] Error: File '{path}' not found.", file=sys.stderr)
        sys.exit(1)

    if path.endswith(".json"):
        with open(path, "r") as f:
            data = json.load(f)
        df = pd.DataFrame(data)
    else:
        df = pd.read_csv(path, sep="\t")

    if id_col not in df.columns:
        available = ", ".join(df.columns.tolist())
        print(
            f"[recipe02] Error: Column '{id_col}' not found in '{path}'.\n"
            f"           Available columns: {available}",
            file=sys.stderr,
        )
        sys.exit(1)

    return df[id_col].dropna().unique().tolist()


def main():
    p = argparse.ArgumentParser(description="Extract genomic overlaps from FILER tracks.")
    p.add_argument("--region", required=True, help='Genomic region, e.g. "chr1:100000-200000"')
    p.add_argument("--track-ids", help="Comma-separated track IDs, e.g. NG123,NG129")
    p.add_argument("--file", help="Path to TSV or JSON file containing track IDs")
    p.add_argument(
        "--id-col",
        default="identifier",
        help="Column name for track IDs in --file (default: identifier)",
    )
    p.add_argument("--out", default="overlaps.tsv", help="Output TSV path (default: overlaps.tsv)")
    p.add_argument("--json", action="store_true", help="Also write JSON alongside the TSV output")
    p.add_argument(
        "--chunk-size",
        type=int,
        default=CHUNK_SIZE,
        help=f"Max track IDs per GET request (default: {CHUNK_SIZE})",
    )

    args = p.parse_args()

    if args.track_ids and args.file:
        print("[recipe02] Warning: both --track-ids and --file given; using --track-ids.", file=sys.stderr)

    if args.track_ids:
        ids = [t.strip() for t in args.track_ids.split(",") if t.strip()]
    elif args.file:
        ids = load_ids_from_file(args.file, args.id_col)
    else:
        print("[recipe02] Error: must provide either --track-ids or --file.", file=sys.stderr)
        p.print_help(sys.stderr)
        sys.exit(1)

    if not ids:
        print("[recipe02] Error: no track IDs resolved.", file=sys.stderr)
        sys.exit(1)

    print(f"[recipe02] Querying region '{args.region}' against {len(ids)} track(s)...", file=sys.stderr)

    df = fetch_overlaps(args.region, ids, chunk_size=args.chunk_size)

    if df.empty:
        print("[recipe02] No overlaps found.", file=sys.stderr)
        sys.exit(0)

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    df.to_csv(args.out, sep="\t", index=False)
    print(f"[recipe02] {len(df)} interval(s) written to '{args.out}'", file=sys.stderr)

    if args.json:
        json_out = (args.out[: -len(".tsv")] if args.out.endswith(".tsv") else args.out) + ".json"
        df.to_json(json_out, orient="records", indent=2)
        print(f"[recipe02] JSON written to '{json_out}'", file=sys.stderr)


if __name__ == "__main__":
    main()