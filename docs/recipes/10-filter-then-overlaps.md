# Recipe 10 — Filter Tracks by Metadata, Find Regional Overlaps, Extract Intervals

## What this does

Runs three FILER endpoints in sequence and joins their outputs into a single
row-per-interval table:

1. **Recipe 1** — query the metadata endpoint to find tracks matching biological
   filters (assay, tissue, data source, etc.)
2. **Recipe 3** — query the coordinate endpoint to find all tracks overlapping a
   genomic region of interest
3. **Intersect** — take the tracks common to both, rank by `num_overlaps`, keep top N
4. **Recipe 2** — fetch the actual overlapping intervals from those top N tracks
5. **Join** — attach full track metadata to every interval row

**Use when:** you have both a biological question (e.g. "ENCODE ATAC-seq blood tracks")
and a genomic locus of interest (e.g. a GWAS hit), and you want the actual peak intervals
from the most relevant tracks at that locus.

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

### Recipe 1 filters (metadata search)

| Parameter | Example values | Notes |
|---|---|---|
| `--assay` | `ATAC-seq`, `ChIP-seq` | |
| `--cell-type` | `CD14+ monocyte` | |
| `--tissue-category` | `Blood`, `Brain` | |
| `--data-source` | `ENCODE`, `Blueprint` | |
| `--track-id` | `NGBLPL2W2SM2WC` | |
| `--filter-string-r1` | `.data_source == "ENCODE" and .assay == "ATAC-seq"` | Overrides named filters |

### Recipe 3 filters (coordinate search)

| Parameter | Example values | Notes |
|---|---|---|
| `--filter-string-r3` | `.tissue_category == "Blood"` | jq-style, default: none |

---

## Outputs

| File | Description |
|---|---|
| `output/10-filter-then-overlaps/results.tsv` | One row per overlapping interval, with full track metadata attached |

**Key columns in the final output:**

| Column | Source | Example value |
|---|---|---|
| `identifier` | All recipes | `NGBLPL2W2SM2WC` |
| `queryRegion` | Recipe 2 | `chr1:100000-200000` |
| `hitString` | Recipe 2 | `chr19@@@44905261@@@44906731@@@Peak_166852...` |
| `num_overlaps` | Recipe 3 | `5` |
| `assay` | Recipe 1 | `ATAC-seq` |
| `tissue_category` | Recipe 1 | `Blood` |
| `data_source` | Recipe 1 | `ENCODE` |
| `track_name` | Recipe 1 | `ENCODE K562 ATAC-seq peaks` |

**Time:** Dominated by Recipe 3 (5–60 s) and Recipe 2 (~1 min per 250 tracks).
Expect 2–5 minutes end-to-end for `--top 100`.

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

Using raw jq filters for finer control over each step:
```bash
python src/scripts/python/filer_filter_then_overlaps.py \
  --genome-build hg38 \
  --region "chr1:100000-200000" \
  --filter-string-r1 '.data_source == "ENCODE" and .assay == "ATAC-seq"' \
  --filter-string-r3 '.tissue_category == "Blood" and .life_stage == "Adult"' \
  --top 50 \
  --out output/10-filter-then-overlaps/results.tsv
```

Querying a smaller region with multiple data sources:
```bash
python src/scripts/python/filer_filter_then_overlaps.py \
  --genome-build hg38 \
  --region "chr19:44905791-44909393" \
  --filter-string-r1 '(.data_source == "ENCODE" or .data_source == "Blueprint") and .assay == "ATAC-seq"' \
  --top 100 \
  --out output/10-filter-then-overlaps/results.tsv
```

### Full script: [`src/scripts/python/filer_filter_then_overlaps.py`](../../src/scripts/python/filer_filter_then_overlaps.py)

---

## How the steps connect
```
Recipe 1 (metadata filter)
  └─ N tracks matching biological criteria
                    ╲
                     ╠═ intersect on identifier, rank by num_overlaps, top N
                    ╱
Recipe 3 (coordinate search)
  └─ M tracks overlapping region
          │
          ▼
    Top N track IDs
          │
          ▼
Recipe 2 (get intervals)
  └─ one row per overlapping interval
          │
          ▼
    Left-join with Recipe 1 metadata
          │
          ▼
    Final table: one row per interval + full metadata
```

---

## Notes

**Intersection size:** if Recipe 1 and Recipe 3 share very few tracks in common
(e.g. a highly specific tissue filter against a region with mostly different data sources),
`--top N` may be smaller than N. The script reports the actual intersection size before
querying Recipe 2.

**`--filter-string-r1` vs `--filter-string-r3`:** these are independent filters applied
at different stages. R1 filters the FILER metadata table directly. R3 filters are applied
server-side after the giggle overlap search, and require `fullMetadata=1` to work — the
script handles this automatically.

:**Hit string:** the final table includes `hitString`, a FILER-formatted string built
from the underlying feature dict values joined with `@@@`. The set/order of fields
(and therefore the number of `@@@`-separated tokens) varies by track type.

**Upstream recipes:** this workflow can also be run step-by-step using the individual
recipe scripts. See [Recipe 1](../01-track-discovery/README.md),
[Recipe 2](../02-track-overlaps/README.md), and
[Recipe 3](../03-coordinate-search/README.md) for details.