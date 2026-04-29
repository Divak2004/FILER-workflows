# Recipe 10 — Filter Tracks by Metadata, Find Regional Overlaps, Extract Intervals

## What this does

Runs two FILER endpoints in sequence and joins their outputs into a single
row-per-interval table:

1. **Recipe 3** — query the coordinate endpoint with a metadata filter to find
   tracks that both match biological criteria (assay, tissue, data source, etc.)
   and overlap a genomic region of interest.
2. **Rank** — sort by `num_overlaps`, keep the top N.
3. **Recipe 2** — fetch the actual overlapping intervals from those top N tracks.
4. **Join** — attach full track metadata to every interval row.

**Use when:** you have both a biological question (e.g. "ENCODE ATAC-seq blood tracks")
and a genomic locus of interest (e.g. a GWAS hit), and you want the actual peak intervals
from the most relevant tracks at that locus.

> **Note:** Earlier versions of this recipe ran the metadata endpoint (Recipe 1) and the
> coordinate endpoint (Recipe 3) separately and intersected the results. With
> `fullMetadata=1`, Recipe 3 already returns the full metadata for each overlapping
> track, so the separate Recipe 1 step is redundant — applying the same filter to
> Recipe 3 produces an identical track set. Recipe 1 is still useful on its own when
> you want the full track universe regardless of region (e.g. for [Recipe 11](./11-filer-selective-install.md)).

---

## Inputs

### Shared parameters

| Parameter | Example values | Required |
|---|---|---|
| `--genome-build` | `hg19`, `hg38` | Yes |
| `--region` | `chr1:100000-200000` | Yes |
| `--top` | `100` (default) | No |
| `--chunk-size` | `250` (default) | No |
| `--out` | `output/10-filter-then-overlaps/results.tsv` | No |

### Metadata filters

These flags mirror Recipe 3 and are joined with `and` into a single jq `filterString`
applied server-side by the coordinate endpoint.

| Parameter | Example values | Notes |
|---|---|---|
| `--assay` | `ATAC-seq`, `ChIP-seq` | |
| `--cell-type` | `CD14+ monocyte` | |
| `--tissue-category` | `Blood`, `Brain` | |
| `--data-source` | `ENCODE`, `Blueprint` | |
| `--track-id` | `NGBLPL2W2SM2WC` | |
| `--filter-string` | `.data_source == "ENCODE" and .assay == "ATAC-seq"` | Raw jq, overrides named filters |

---

## Outputs

| File | Description |
|---|---|
| `output/10-filter-then-overlaps/results.tsv` | One row per overlapping interval, with full track metadata attached |

**Key columns in the final output:**

| Column | Source | Example value |
|---|---|---|
| `identifier` | Recipe 3 | `NGBLPL2W2SM2WC` |
| `queryRegion` | Recipe 2 | `chr1:100000-200000` |
| `hitString` | Recipe 2 | `chr19@@@44905261@@@44906731@@@Peak_166852...` |
| `num_overlaps` | Recipe 3 | `5` |
| `assay` | Recipe 3 | `ATAC-seq` |
| `tissue_category` | Recipe 3 | `Blood` |
| `data_source` | Recipe 3 | `ENCODE` |
| `track_name` | Recipe 3 | `ENCODE K562 ATAC-seq peaks` |
| `processed_file_download_url` | Recipe 3 | `https://tf.lisanwanglab.org/GADB/…bed.gz` |
| `processed_file_md5` | Recipe 3 | `3b5013b3bed9eb54...` |
| `wget_command` | Recipe 3 | `wget https://… -P FILER2/Annotationtracks/…` |

**Time:** Dominated by Recipe 3 (5–60 s) and Recipe 2 (~1 min per 250 tracks).
Expect 1–3 minutes end-to-end for `--top 100`.

---

## Prerequisites
```bash
# 1. Ensure you are in the FILER-Workflows folder

# 2. Ensure your venv is active
source venv/bin/activate

# 3. Install the project and its core dependencies
pip install -e .
```

No authentication required for public FILER data.

---

## Python

### Example scripts

Basic usage — ENCODE ATAC-seq blood tracks overlapping a region:
```bash
python src/scripts/python/filer_filter_then_overlaps.py \
  --genome-build hg38 \
  --region "chr1:100000-200000" \
  --assay "ATAC-seq" \
  --tissue-category "Blood" \
  --data-source "ENCODE" \
  --top 100 \
  --out output/10-filter-then-overlaps/results.tsv
```

Using a raw jq filter for finer control (e.g. OR logic across data sources):
```bash
python src/scripts/python/filer_filter_then_overlaps.py \
  --genome-build hg38 \
  --region "chr1:100000-200000" \
  --filter-string '(.data_source == "ENCODE" or .data_source == "Blueprint") and .assay == "ATAC-seq"' \
  --top 50 \
  --out output/10-filter-then-overlaps/results.tsv
```

Querying a smaller region with a tighter filter:
```bash
python src/scripts/python/filer_filter_then_overlaps.py \
  --genome-build hg38 \
  --region "chr19:44905791-44909393" \
  --filter-string '.data_source == "ENCODE" and .assay == "ATAC-seq" and .life_stage == "Adult"' \
  --top 100 \
  --out output/10-filter-then-overlaps/results.tsv
```

### Full script: [`src/scripts/python/filer_filter_then_overlaps.py`](../../src/scripts/python/filer_filter_then_overlaps.py)

---

## How the steps connect
```
Recipe 3 (coordinate search with metadata filter, fullMetadata=1)
  └─ N tracks matching biological criteria AND overlapping the region
          │
          ▼
    Sort by num_overlaps, keep top N
          │
          ▼
Recipe 2 (get intervals)
  └─ one row per overlapping interval
          │
          ▼
    Left-join with track metadata
          │
          ▼
    Final table: one row per interval + full metadata
```

---

## Notes

**`--filter-string` requires `fullMetadata=1`:** the script always sets this on the
coordinate endpoint so that filters can be evaluated server-side and so that the
returned metadata is rich enough for downstream use (e.g. [Recipe 11](./11-filer-selective-install.md)
needs `processed_file_md5`, `file_size`, etc.).

**Hit string:** the final table includes `hitString`, a FILER-formatted string built
from the underlying feature dict values joined with `@@@`. The set/order of fields
(and therefore the number of `@@@`-separated tokens) varies by track type.

**Tracks with biology match but zero overlap:** because we use Recipe 3 with the
filter applied server-side, tracks that match the biological filter but have no
overlap in the region simply don't appear in the output. If you want to know the
size of the full biological universe regardless of region, run [Recipe 1](./01-track-discovery.md)
separately with the same filter.

**Upstream recipes:** the individual steps can also be run on their own. See
[Recipe 2](./02-track-overlaps.md) and [Recipe 3](./03-coordinate-search.md) for details.
