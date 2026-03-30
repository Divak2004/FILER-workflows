import requests
import sys
import time
from typing import Union

SERVER = "https://rest.ensembl.org"
HEADERS = {"Content-Type": "application/json", "Accept": "application/json"}


def map_rsid_to_region(rsid: str, genome_build: str = "GRCh38") -> Union[list[dict], str]:
    """
    Map a single rsID to its genomic position(s) using the Ensembl REST API.

    Uses the /variation/human/{id} endpoint, which returns structured mappings
    with chromosome, start, and end directly — no HGVS string parsing needed.

    A single rsID can map to multiple locations (e.g. positions in assembly
    exceptions), so a list of dicts is returned, one per mapping.

    Parameters
    ----------
    rsid : str
        dbSNP rsID, e.g. "rs699".
    genome_build : str
        "GRCh38" (default) or "GRCh37". Used to filter mappings by assembly.

    Returns
    -------
    list[dict] | str
        On success: list of dicts with keys:
            rsid, chromosome, start, end, region, allele_string, strand
        On failure: error string.
    """
    url = f"{SERVER}/variation/human/{rsid}"

    try:
        r = requests.get(url, headers=HEADERS, timeout=30)
        r.raise_for_status()
    except requests.exceptions.HTTPError as e:
        if r.status_code == 400:
            return f"Error: '{rsid}' not recognised as a valid rsID."
        if r.status_code == 404:
            return f"Error: rsID '{rsid}' not found in Ensembl."
        return f"HTTP error for '{rsid}': {e}"
    except requests.exceptions.RequestException as e:
        return f"Connection error for '{rsid}': {e}"

    data = r.json()
    mappings = data.get("mappings", [])

    if not mappings:
        return f"Error: no genomic mappings found for '{rsid}'."

    results = []
    for m in mappings:
        # Filter by assembly if specified
        assembly = m.get("assembly_name", "")
        if genome_build and genome_build not in assembly:
            continue

        chrom = m.get("seq_region_name", "?")
        start = m.get("start")
        end   = m.get("end")

        # Format as UCSC/FILER-style region string (chr-prefixed, 0-based start)
        region = f"chr{chrom}:{start - 1}-{end}"

        results.append({
            "rsid":          rsid,
            "chromosome":    f"chr{chrom}",
            "start":         start,          # 1-based (Ensembl convention)
            "end":           end,
            "region":        region,         # 0-based start, matches FILER format
            "allele_string": m.get("allele_string", ""),
            "strand":        m.get("strand", ""),
            "assembly":      assembly,
        })

    if not results:
        return f"Error: no mappings for '{rsid}' matching assembly '{genome_build}'."

    return results


def map_rsids_batch(rsids: list[str], genome_build: str = "GRCh38",
                    chunk_size: int = 200) -> dict:
    """
    Map a list of rsIDs to genomic positions using the POST batch endpoint.

    Splits into chunks of chunk_size (Ensembl limit is ~200 per request) and
    pauses briefly between chunks to respect the 15 req/s rate limit.

    Parameters
    ----------
    rsids : list[str]
        List of rsIDs, e.g. ["rs699", "rs1800562"].
    genome_build : str
        "GRCh38" (default) or "GRCh37".
    chunk_size : int
        IDs per POST request. Max 200.

    Returns
    -------
    dict
        Keys are rsIDs. Values are either a list of mapping dicts (see
        map_rsid_to_region) or an error string.
    """
    url = f"{SERVER}/variation/homo_sapiens"
    results = {}

    chunks = [rsids[i:i + chunk_size] for i in range(0, len(rsids), chunk_size)]

    for i, chunk in enumerate(chunks):
        if i > 0:
            time.sleep(0.5)  # stay well under the 15 req/s rate limit

        try:
            r = requests.post(
                url,
                headers=HEADERS,
                json={"ids": chunk},
                timeout=60,
            )
            r.raise_for_status()
        except requests.exceptions.RequestException as e:
            for rsid in chunk:
                results[rsid] = f"Request error: {e}"
            continue

        batch_data = r.json()  # dict keyed by rsID

        for rsid in chunk:
            entry = batch_data.get(rsid)
            if entry is None:
                results[rsid] = f"Error: '{rsid}' not found in Ensembl."
                continue

            mappings = entry.get("mappings", [])
            if not mappings:
                results[rsid] = f"Error: no genomic mappings found for '{rsid}'."
                continue

            hits = []
            for m in mappings:
                assembly = m.get("assembly_name", "")
                if genome_build and genome_build not in assembly:
                    continue

                chrom  = m.get("seq_region_name", "?")
                start  = m.get("start")
                end    = m.get("end")
                region = f"chr{chrom}:{start - 1}-{end}"

                hits.append({
                    "rsid":          rsid,
                    "chromosome":    f"chr{chrom}",
                    "start":         start,
                    "end":           end,
                    "region":        region,
                    "allele_string": m.get("allele_string", ""),
                    "strand":        m.get("strand", ""),
                    "assembly":      assembly,
                })

            results[rsid] = hits if hits else f"Error: no mappings for '{rsid}' matching assembly '{genome_build}'."

    return results


# ── Example usage ─────────────────────────────────────────────────────────────

if __name__ == "__main__":
    # Single rsID
    print("=== Single lookup ===")
    result = map_rsid_to_region("rs699")
    if isinstance(result, str):
        print(result)  # error
    else:
        for hit in result:
            print(f"{hit['rsid']}  →  {hit['region']}  ({hit['allele_string']})")

    print()

    # Batch lookup
    print("=== Batch lookup ===")
    batch = map_rsids_batch(["rs699", "rs1800562", "rs7903146"])
    for rsid, hits in batch.items():
        if isinstance(hits, str):
            print(f"{rsid}: {hits}")
        else:
            for hit in hits:
                print(f"{hit['rsid']}  →  {hit['region']}  ({hit['allele_string']})")