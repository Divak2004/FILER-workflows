#!/usr/bin/env python3
"""
Recipe 04 — Map rsIDs to Genomic Regions

Converts one or more dbSNP rsIDs to genomic coordinates using the Ensembl
REST API and writes results as a TSV with one row per rsID–mapping.

The `region` column is formatted as chrN:start-end (0-based start) and is
ready to pass directly to filer_coordinate_search.py (Recipe 3) or
filer_filter_then_overlaps.py (Recipe 10).

Usage:
  # Single rsID
  python src/scripts/python/rsid_to_positions.py \\
    --rsid rs699 \\
    --genome-build GRCh38 \\
    --out output/11-rsid-to-region/results.tsv

  # Batch from file (one rsID per line, # lines ignored)
  python src/scripts/python/rsid_to_positions.py \\
    --file input/rsids.txt \\
    --genome-build GRCh38 \\
    --out output/11-rsid-to-region/results.tsv
"""

import argparse
import os
import sys
import time
from typing import Union

import pandas as pd
import requests

SERVER  = "https://rest.ensembl.org"
HEADERS = {"Content-Type": "application/json", "Accept": "application/json"}

OUTPUT_COLS = ["rsid", "chromosome", "start", "end", "region",
               "allele_string", "strand", "assembly"]


# ── Core lookup helpers ───────────────────────────────────────────────────────

def _parse_mappings(rsid: str, mappings: list, genome_build: str) -> list[dict]:
    """Convert a list of Ensembl mapping dicts into output rows."""
    rows = []
    for m in mappings:
        assembly = m.get("assembly_name", "")
        if genome_build and genome_build not in assembly:
            continue

        chrom  = m.get("seq_region_name", "?")
        start  = m.get("start")
        end    = m.get("end")
        # 0-based start to match UCSC / FILER chrN:start-end convention
        region = f"chr{chrom}:{start - 1}-{end}"

        rows.append({
            "rsid":          rsid,
            "chromosome":    f"chr{chrom}",
            "start":         start,
            "end":           end,
            "region":        region,
            "allele_string": m.get("allele_string", ""),
            "strand":        m.get("strand", ""),
            "assembly":      assembly,
        })
    return rows


def map_rsid_single(rsid: str, genome_build: str = "GRCh38") -> Union[list[dict], str]:
    """
    Map one rsID via GET /variation/human/{id}.

    Returns a list of mapping dicts on success, or an error string on failure.
    """
    url = f"{SERVER}/variation/human/{rsid}"
    try:
        r = requests.get(url, headers=HEADERS, timeout=30)
        r.raise_for_status()
    except requests.exceptions.HTTPError:
        if r.status_code == 400:
            return f"Error: '{rsid}' not recognised as a valid rsID."
        if r.status_code == 404:
            return f"Error: rsID '{rsid}' not found in Ensembl."
        return f"HTTP {r.status_code} for '{rsid}'."
    except requests.exceptions.RequestException as e:
        return f"Connection error for '{rsid}': {e}"

    mappings = r.json().get("mappings", [])
    if not mappings:
        return f"Error: no genomic mappings found for '{rsid}'."

    rows = _parse_mappings(rsid, mappings, genome_build)
    if not rows:
        return f"Error: no mappings for '{rsid}' matching assembly '{genome_build}'."
    return rows


def map_rsids_batch(rsids: list[str], genome_build: str = "GRCh38",
                    chunk_size: int = 200) -> dict:
    """
    Map a list of rsIDs via POST /variation/homo_sapiens.

    Automatically splits into chunks and pauses between them to respect
    the Ensembl 15 req/s rate limit.

    Returns a dict keyed by rsID; values are lists of mapping dicts or
    error strings.
    """
    url     = f"{SERVER}/variation/homo_sapiens"
    results = {}
    chunks  = [rsids[i:i + chunk_size] for i in range(0, len(rsids), chunk_size)]

    for i, chunk in enumerate(chunks):
        if i > 0:
            time.sleep(0.5)

        print(f"  batch {i + 1}/{len(chunks)} ({len(chunk)} IDs)…", file=sys.stderr)

        try:
            r = requests.post(url, headers=HEADERS, json={"ids": chunk}, timeout=60)
            r.raise_for_status()
        except requests.exceptions.RequestException as e:
            for rsid in chunk:
                results[rsid] = f"Request error: {e}"
            continue

        batch_data = r.json()

        for rsid in chunk:
            entry    = batch_data.get(rsid)
            if entry is None:
                results[rsid] = f"Error: '{rsid}' not found in Ensembl."
                continue

            mappings = entry.get("mappings", [])
            if not mappings:
                results[rsid] = f"Error: no genomic mappings found for '{rsid}'."
                continue

            rows = _parse_mappings(rsid, mappings, genome_build)
            results[rsid] = rows if rows else \
                f"Error: no mappings for '{rsid}' matching assembly '{genome_build}'."

    return results


# ── Input loading ─────────────────────────────────────────────────────────────

def load_rsids_from_file(path: str) -> list[str]:
    """Read rsIDs from a plain-text file, one per line. # lines are ignored."""
    if not os.path.exists(path):
        raise FileNotFoundError(f"Input file not found: '{path}'")

    rsids = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line and not line.startswith("#"):
                rsids.append(line)

    if not rsids:
        raise ValueError(f"No rsIDs found in '{path}'.")

    return rsids


# ── CLI ───────────────────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser(
        description="Map rsIDs to genomic regions using the Ensembl REST API.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--rsid",  help="Single rsID, e.g. rs699")
    src.add_argument("--file",  help="Text file with one rsID per line (# lines ignored)")

    p.add_argument("--genome-build", default="GRCh38",
                   help="Assembly to filter on: GRCh38 (default) or GRCh37. "
                        "Pass '' to return all assemblies.")
    p.add_argument("--chunk-size",   type=int, default=200,
                   help="rsIDs per POST request for batch mode (default: 200, max: 200)")
    p.add_argument("--out", default="output/11-rsid-to-region/results.tsv",
                   help="Output TSV path (default: output/11-rsid-to-region/results.tsv)")

    args = p.parse_args()

    # ── Collect rsIDs ─────────────────────────────────────────────────────────
    if args.rsid:
        rsids = [args.rsid]
    else:
        try:
            rsids = load_rsids_from_file(args.file)
        except (FileNotFoundError, ValueError) as e:
            print(f"[rsid_to_positions] {e}", file=sys.stderr)
            sys.exit(1)

    print(f"[rsid_to_positions] Mapping {len(rsids)} rsID(s) to {args.genome_build}…",
          file=sys.stderr)

    # ── Resolve ───────────────────────────────────────────────────────────────
    if len(rsids) == 1:
        result = map_rsid_single(rsids[0], genome_build=args.genome_build)
        if isinstance(result, str):
            print(f"[rsid_to_positions] {result}", file=sys.stderr)
            sys.exit(1)
        all_rows = result
        errors   = []
    else:
        batch    = map_rsids_batch(rsids, genome_build=args.genome_build,
                                   chunk_size=args.chunk_size)
        all_rows = []
        errors   = []
        for rsid, val in batch.items():
            if isinstance(val, str):
                errors.append(f"  {rsid}: {val}")
            else:
                all_rows.extend(val)

    # ── Report errors ─────────────────────────────────────────────────────────
    if errors:
        print(f"[rsid_to_positions] {len(errors)} rsID(s) could not be resolved:",
              file=sys.stderr)
        for e in errors:
            print(e, file=sys.stderr)

    if not all_rows:
        print("[rsid_to_positions] No mappings resolved. Exiting.", file=sys.stderr)
        sys.exit(1)

    # ── Write output ──────────────────────────────────────────────────────────
    df = pd.DataFrame(all_rows, columns=OUTPUT_COLS)

    out_dir = os.path.dirname(args.out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    df.to_csv(args.out, sep="\t", index=False)
    print(f"[rsid_to_positions] {len(df)} mapping(s) written to {args.out}",
          file=sys.stderr)

    # Print region strings to stdout for easy shell piping
    for region in df["region"].unique():
        print(region)


if __name__ == "__main__":
    main()