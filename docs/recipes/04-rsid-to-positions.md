# Recipe 4 — Map rsIDs to Genomic Regions

## What this does

Converts one or more dbSNP rsIDs to genomic coordinates using the Ensembl REST API,
returning results in the `chrN:start-end` format accepted by all other FILER recipes.

Two modes are supported:

- **Single lookup** — one rsID via `GET /variation/human/{id}`
- **Batch lookup** — a list of rsIDs via `POST /variation/homo_sapiens` (up to 200 per request, automatically chunked)

The output is a TSV with one row per rsID–mapping pair, ready to be piped directly into
Recipe 3 (coordinate search) or Recipe 10 (full end-to-end workflow).

---

## Inputs

### Parameters

| Parameter | Example values | Required |
|---|---|---|
| `--rsid` | `rs699` | One of `--rsid` or `--file` |
| `--file` | `input/rsids.txt` | One of `--rsid` or `--file` |
| `--genome-build` | `GRCh38` (default), `GRCh37` | No |
| `--chunk-size` | `200` (default) | No |
| `--out` | `output/11-rsid-to-positions/results.tsv` | No |

### Input file format (for `--file`)

Plain text, one rsID per line, `#` lines ignored:
```
# GWAS hits for my locus
rs699
rs1800562
rs7903146
```

---

## Outputs

| File | Description |
|---|---|
| `output/11-rsid-to-positions/results.tsv` | One row per rsID–mapping, with region string and allele info |

**Key columns:**

| Column | Example value | Notes |
|---|---|---|
| `rsid` | `rs699` | Input rsID |
| `chromosome` | `chr1` | chr-prefixed |
| `start` | `230710048` | 1-based (Ensembl convention) |
| `end` | `230710048` | 1-based |
| `region` | `chr1:230710047-230710048` | 0-based start; ready for FILER |
| `allele_string` | `A/G` | ref/alt |
| `strand` | `1` | |
| `assembly` | `GRCh38` | |

> ℹ️ **Coordinate conventions:** Ensembl returns 1-based positions. The `region` column
> uses a 0-based start (`start - 1`) to match the UCSC/FILER `chrN:start-end` format
> expected by Recipes 2, 3, and 10.

> ℹ️ **Multiple mappings:** some variants map to more than one location (e.g. PAR regions,
> assembly exceptions). Each mapping becomes its own row; most rsIDs will have exactly one.

---

## Prerequisites
```bash
# 1. Ensure you are in the FILER-Workflows folder

# 2. Ensure your venv is active
source venv/bin/activate

# 3. Install the project and its core dependencies
pip install -e .
```

No authentication required. The Ensembl public REST API is free and requires no API key.
The public server is rate-limited to **15 requests/second**; the script handles this
automatically with a brief sleep between batch chunks.

---

## Python

### Example scripts

Single rsID lookup:
```bash
python src/scripts/python/rsid_to_positions.py \
  --rsid rs699 \
  --genome-build GRCh38 \
  --out output/11-rsid-to-positions/results.tsv
```

Batch lookup from a file:
```bash
python src/scripts/python/rsid_to_positions.py \
  --file input/rsids.txt \
  --genome-build GRCh38 \
  --out output/11-rsid-to-positions/results.tsv
```

Use GRCh37 (hg19) coordinates:
```bash
python src/scripts/python/rsid_to_positions.py \
  --file input/rsids.txt \
  --genome-build GRCh37 \
  --out output/11-rsid-to-positions/results.tsv
```

### Full script: [`src/scripts/python/rsid_to_positions.py`](../../src/scripts/python/rsid_to_positions.py)

---

## Chaining into downstream recipes

The `region` column in the output is formatted as `chrN:start-end` and can be passed
directly to Recipe 3 or Recipe 10.

### Pipe a single rsID into Recipe 3

```bash
# Step 1: resolve rsID to region
python src/scripts/python/rsid_to_positions.py \
  --rsid rs699 \
  --out output/11-rsid-to-positions/results.tsv

# Step 2: extract the region string and pass to Recipe 3
REGION=$(awk 'NR==2 {print $5}' output/11-rsid-to-positions/results.tsv)

python src/scripts/python/filer_coordinate_search.py \
  --region "$REGION" \
  --genome-build hg38 \
  --count-only 0 \
  --full-metadata \
  --out output/03-coordinate-search/search_results.tsv
```

### Pipe a single rsID into Recipe 2

```bash
  # Pipeline — rsID → Recipe 3 → Recipe 2 (coordinate search → overlapping intervals):
  # Step 1: resolve rsID to region
  python src/scripts/python/rsid_to_positions.py \
    --rsid rs699 \
    --out output/11-rsid-to-positions/results.tsv

  # Step 2: extract the region string and pass to Recipe 3
  REGION=$(awk 'NR==2 {print $5}' output/11-rsid-to-positions/results.tsv)

  python src/scripts/python/filer_coordinate_search.py \
    --region "$REGION" \
    --genome-build hg38 \
    --count-only 0 \
    --out output/03-coordinate-search/search_results.tsv

  # Step 3: extract overlapping intervals from discovered tracks
  python src/scripts/python/filer_find_overlaps.py \
    --region "$REGION" \
    --file output/03-coordinate-search/search_results.tsv \
    --id-col identifier \
    --out output/02-track-overlaps/overlaps.tsv
```

### Pipe a batch into Recipe 10

```bash
# Step 1: resolve all rsIDs
python src/scripts/python/rsid_to_positions.py \
  --file input/rsids.txt \
  --out output/11-rsid-to-positions/results.tsv

# Step 2: loop over each resolved region
tail -n +2 output/11-rsid-to-positions/results.tsv | cut -f5 | while read REGION; do
  python src/scripts/python/filer_filter_then_overlaps.py \
    --genome-build hg38 \
    --region "$REGION" \
    --assay "ATAC-seq" \
    --tissue-category "Blood" \
    --top 100 \
    --out "output/10-filter-then-overlaps/${REGION//[^a-zA-Z0-9]/_}.tsv"
done
```

---

## Notes

**GRCh37 vs GRCh37.p13:** the `--genome-build` filter matches on substring, so `GRCh37`
matches both `GRCh37` and `GRCh37.p13`. Pass `""` to return all assemblies unfiltered.

**Upstream / downstream recipes:**

| Recipe | Description |
|---|---|
| [Recipe 3](../03-coordinate-search/README.md) | Find FILER tracks overlapping a region |
| [Recipe 10](../10-filter-then-overlaps/README.md) | Full end-to-end: metadata filter → region search → intervals |