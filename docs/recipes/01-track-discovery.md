# Recipe 1 — List Tracks Matching Metadata Filters

## What this does

Given a set of standardized metadata filters (cell type, assay, biosample type, life stage, data source), query the FILER API and return a table of matching tracks saved as `tracks.tsv` (and optionally `tracks.json`).

**Use when:** you want to assemble a relevant set of tracks before querying overlaps. Almost every downstream workflow starts here.

---

## Inputs

| Parameter | Example values | Required |
|---|---|---|
| `genome_build` | `hg19`, `hg38` | Yes |
| `assay` | `ATAC-seq`, `ChIP-seq`, `DNase-seq` | No |
| `cell_type` | `CD14+ monocyte`, `T cell` | No |
| `biosample_type` | `primary cell`, `tissue` | No |
| `life_stage` | `adult`, `embryonic` | No |
| `data_source` | `ENCODE`, `Roadmap` | No |

At least one filter beyond `genome_build` is recommended to keep results manageable.

## Outputs

| File | Description |
|---|---|
| `tracks.tsv` | Tab-separated table of matching tracks |
| `tracks.json` | Same content as JSON array (optional) |

**Required columns in `tracks.tsv`:**
`track_id`, `genome_build`, `assay`, `cell_type`, `data_source`, `bed_gz_url`

**Optional columns:** `biosample_type`, `life_stage`, `tissue`, `biosample`, `experiment_id`

**Time:** Typically seconds for a filtered query; up to ~30s for broad/unfiltered queries returning thousands of tracks.

---

## Prerequisites

```bash
# Python
pip install requests pandas

# R
install.packages(c("httr", "jsonlite", "dplyr"))
```

Access to the FILER API endpoint (no authentication required for public data).

---

## Python (recommended)

### Quick one-liner (script)

```bash
python scripts/python/filer_search_tracks.py \
  --genome-build hg38 \
  --assay "ATAC-seq" \
  --cell-type "CD14+ monocyte" \
  --data-source ENCODE \
  --out tracks.tsv \
  --json   # also write tracks.json
```

### Full script: `scripts/python/filer_search_tracks.py`

```python
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
```

### Notebook / inline snippet

```python
import requests
import pandas as pd

BASE_URL = "https://tf.lisanwanglab.org/FILER2"

params = {
    "genome_build": "hg38",
    "assay":        "ATAC-seq",
    "cell_type":    "CD14+ monocyte",
    "data_source":  "ENCODE",
}

r = requests.get(f"{BASE_URL}/api/tracks", params=params, timeout=60)
r.raise_for_status()

df = pd.DataFrame(r.json()["results"])  # adjust key if needed
print(df.shape)
df.head()
```

---

## Bash (minimal)

```bash
BASE="https://tf.lisanwanglab.org/FILER2"

curl -sG "${BASE}/api/tracks" \
  --data-urlencode "genome_build=hg38" \
  --data-urlencode "assay=ATAC-seq" \
  --data-urlencode "cell_type=CD14+ monocyte" \
  --data-urlencode "data_source=ENCODE" \
  > tracks.json

# Convert to TSV (requires jq)
jq -r '
  ["track_id","genome_build","assay","cell_type","data_source","bed_gz_url"],
  (.results[] | [.track_id,.genome_build,.assay,.cell_type,.data_source,.bed_gz_url])
  | @tsv
' tracks.json > tracks.tsv
```

---

## R (secondary)

```r
library(httr)
library(jsonlite)
library(dplyr)

BASE_URL <- "https://tf.lisanwanglab.org/FILER2"

search_tracks <- function(genome_build = "hg38", ...) {
  filters <- Filter(Negate(is.null), list(genome_build = genome_build, ...))
  resp    <- GET(paste0(BASE_URL, "/api/tracks"), query = filters)
  stop_for_status(resp)
  as_tibble(fromJSON(content(resp, as = "text", encoding = "UTF-8"))$results)
}

tracks <- search_tracks(
  genome_build = "hg38",
  assay        = "ATAC-seq",
  cell_type    = "CD14+ monocyte",
  data_source  = "ENCODE"
)

write.table(tracks, "tracks.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
print(tracks)
```

---

## Expected output

`tracks.tsv` (first 3 rows, columns abbreviated):

```
track_id          genome_build  assay     cell_type         data_source  bed_gz_url
ENCFF001ABC       hg38          ATAC-seq  CD14+ monocyte    ENCODE       https://...
ENCFF002DEF       hg38          ATAC-seq  CD14+ monocyte    ENCODE       https://...
ENCFF003GHI       hg38          ATAC-seq  CD14+ monocyte    ENCODE       https://...
```

The fixture at `examples/expected/recipe01_tracks.tsv` contains a known-good 20-row result you can diff against.

---

## Discover available filter values

Not sure which values are valid for a given field? Use the metadata values endpoint:

```bash
# List all available assay types
curl -sG "${BASE}/api/metadata/values" \
  --data-urlencode "field=assay" \
  --data-urlencode "include_counts=true" \
  | jq '.values[:10]'

# Filter cell types relevant to a specific assay context
curl -sG "${BASE}/api/metadata/values" \
  --data-urlencode "field=cell_type" \
  --data-urlencode "q=mono" \
  --data-urlencode "filters=assay:ATAC-seq" \
  --data-urlencode "sort=count" \
  | jq '.values'
```

```python
# Python equivalent
r = requests.get(f"{BASE_URL}/api/metadata/values", params={
    "field": "cell_type",
    "q": "mono",
    "filters": "assay:ATAC-seq",
    "sort": "count",
    "include_counts": "true",
}, timeout=30)
print(r.json()["values"])
```

---

## Scale-up notes

**Large result sets:** If a broad query returns thousands of tracks, check whether the endpoint supports `limit` / `offset` or `page` / `page_size` pagination and loop accordingly.

**Caching:** Save `tracks.tsv` with a datestamp (e.g., `tracks_20250224.tsv`) so runs are reproducible. The Track Set manifest in Recipe 1.3 formalizes this pattern.

**Batching multiple filter combinations:** Loop over a list of filter dicts and concatenate results into a single DataFrame, then deduplicate on `track_id`.

```python
filter_sets = [
    {"assay": "ATAC-seq",  "cell_type": "CD14+ monocyte"},
    {"assay": "DNase-seq", "cell_type": "CD14+ monocyte"},
]
df = pd.concat([search_tracks("hg38", f) for f in filter_sets]).drop_duplicates("track_id")
```

---

## Next steps

- **Recipe 1.3** — Turn this track list into a reusable Track Set manifest (with provenance, date, filters)
- **Recipe 3.1** — Query overlaps for a genomic region against this track set
- **Recipe 4.1** — Summarize and plot overlap counts by assay / cell type