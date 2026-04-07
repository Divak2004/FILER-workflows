#!/usr/bin/env python3
"""
Recipe 05 — Map gene symbols to Genomic Regions

Converts one or more gene symbols to genomic coordinates using the Ensembl
REST API and writes results as a TSV with one row per gene.

The `region` column is formatted as chrN:start-end (0-based start) and is
ready to pass directly to filer_coordinate_search.py (Recipe 3).

Usage:
  # Single gene
  python src/scripts/python/gene_to_positions.py \
    --gene AGT \
    --genome-build GRCh38 \
    --out output/5-gene-to-positions/results.tsv

  # Batch from file (one gene symbol per line, # lines ignored)
  python src/scripts/python/gene_to_positions.py \
    --file input/genes.txt \
    --genome-build GRCh38 \
    --out output/5-gene-to-positions/results.tsv
"""

import argparse
import os
import sys
import time

import pandas as pd
import requests

SERVER  = "https://rest.ensembl.org"
HEADERS = {"Content-Type": "application/json", "Accept": "application/json"}

OUTPUT_COLS = ["gene", "ensembl_id", "chromosome", "start", "end", "region",
               "strand", "assembly"]


# ── Core lookup helpers ───────────────────────────────────────────────────────

def _parse_gene(gene: str, data: dict) -> dict:
    """Convert an Ensembl lookup response into an output row."""
    chrom  = data.get("seq_region_name", "?")
    start  = data.get("start")
    end    = data.get("end")
    region = f"chr{chrom}:{start - 1}-{end}"

    return {
        "gene":        gene,
        "ensembl_id":  data.get("id", ""),
        "chromosome":  f"chr{chrom}",
        "start":       start,
        "end":         end,
        "region":      region,
        "strand":      data.get("strand", ""),
        "assembly":    data.get("assembly_name", ""),
    }


def map_gene_single(gene: str) -> dict | str:
    """
    Map one gene symbol via GET /lookup/symbol/homo_sapiens/{gene}.

    Returns a row dict on success, or an error string on failure.
    """
    url = f"{SERVER}/lookup/symbol/homo_sapiens/{gene}"
    try:
        r = requests.get(url, headers=HEADERS, params={"expand": 0}, timeout=30)
        r.raise_for_status()
    except requests.exceptions.HTTPError:
        if r.status_code == 400:
            return f"Error: '{gene}' not recognised as a valid gene symbol."
        if r.status_code == 404:
            return f"Error: gene '{gene}' not found in Ensembl."
        return f"HTTP {r.status_code} for '{gene}'."
    except requests.exceptions.RequestException as e:
        return f"Connection error for '{gene}': {e}"

    return _parse_gene(gene, r.json())


def map_genes_batch(genes: list[str], chunk_size: int = 1000) -> dict:
    """
    Map a list of gene symbols via POST /lookup/symbol/homo_sapiens.

    Returns a dict keyed by gene symbol; values are row dicts or error strings.
    """
    url     = f"{SERVER}/lookup/symbol/homo_sapiens"
    results = {}
    chunks  = [genes[i:i + chunk_size] for i in range(0, len(genes), chunk_size)]

    for i, chunk in enumerate(chunks):
        if i > 0:
            time.sleep(0.5)

        print(f"  batch {i + 1}/{len(chunks)} ({len(chunk)} genes)…", file=sys.stderr)

        try:
            r = requests.post(url, headers=HEADERS, json={"symbols": chunk}, timeout=60)
            r.raise_for_status()
        except requests.exceptions.RequestException as e:
            for gene in chunk:
                results[gene] = f"Request error: {e}"
            continue

        batch_data = r.json()

        for gene in chunk:
            entry = batch_data.get(gene)
            if entry is None:
                results[gene] = f"Error: '{gene}' not found in Ensembl."
            else:
                results[gene] = _parse_gene(gene, entry)

    return results


# ── Input loading ─────────────────────────────────────────────────────────────

def load_genes_from_file(path: str) -> list[str]:
    """Read gene symbols from a plain-text file, one per line. # lines ignored."""
    if not os.path.exists(path):
        raise FileNotFoundError(f"Input file not found: '{path}'")

    genes = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line and not line.startswith("#"):
                genes.append(line)

    if not genes:
        raise ValueError(f"No gene symbols found in '{path}'.")

    return genes


# ── CLI ───────────────────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser(
        description="Map gene symbols to genomic regions using the Ensembl REST API.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--gene", help="Single gene symbol, e.g. AGT")
    src.add_argument("--file", help="Text file with one gene symbol per line (# lines ignored)")

    p.add_argument("--genome-build", default="GRCh38",
                   help="Assembly to filter on: GRCh38 (default) or GRCh37.")
    p.add_argument("--chunk-size", type=int, default=1000,
                   help="Genes per POST request for batch mode (default: 1000)")
    p.add_argument("--out", default="output/5-gene-to-positions/results.tsv",
                   help="Output TSV path")

    args = p.parse_args()

    if args.gene:
        genes = [args.gene]
    else:
        try:
            genes = load_genes_from_file(args.file)
        except (FileNotFoundError, ValueError) as e:
            print(f"[gene_to_positions] {e}", file=sys.stderr)
            sys.exit(1)

    print(f"[gene_to_positions] Mapping {len(genes)} gene(s)…", file=sys.stderr)

    if len(genes) == 1:
        result = map_gene_single(genes[0])
        if isinstance(result, str):
            print(f"[gene_to_positions] {result}", file=sys.stderr)
            sys.exit(1)
        all_rows = [result]
        errors   = []
    else:
        batch    = map_genes_batch(genes, chunk_size=args.chunk_size)
        all_rows = []
        errors   = []
        for gene, val in batch.items():
            if isinstance(val, str):
                errors.append(f"  {gene}: {val}")
            else:
                all_rows.append(val)

    if errors:
        print(f"[gene_to_positions] {len(errors)} gene(s) could not be resolved:",
              file=sys.stderr)
        for e in errors:
            print(e, file=sys.stderr)

    if not all_rows:
        print("[gene_to_positions] No mappings resolved. Exiting.", file=sys.stderr)
        sys.exit(1)

    df = pd.DataFrame(all_rows, columns=OUTPUT_COLS)

    out_dir = os.path.dirname(args.out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    df.to_csv(args.out, sep="\t", index=False)
    print(f"[gene_to_positions] {len(df)} gene(s) written to {args.out}", file=sys.stderr)

    # Print region strings to stdout for easy shell piping
    for region in df["region"].unique():
        print(region)


if __name__ == "__main__":
    main()