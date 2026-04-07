# Recipe 5 — Map Gene Symbols to Genomic Regions

## What this does

Converts one or more gene symbols to genomic coordinates using the Ensembl REST API,
returning results in the `chrN:start-end` format accepted by all other FILER recipes.

Two modes are supported:

- **Single lookup** — one gene symbol via `GET /lookup/symbol/homo_sapiens/{gene}`
- **Batch lookup** — a list of gene symbols via `POST /lookup/symbol/homo_sapiens` (up to 1000 per request, automatically chunked)

The output is a TSV with one row per gene, ready to be piped directly into
Recipe 3 (coordinate search) or Recipe 10 (full end-to-end workflow).

---

## Inputs

### Parameters

| Parameter | Example values | Required |
|---|---|---|
| `--gene` | `AGT` | One of `--gene` or `--file` |
| `--file` | `input/genes.txt` | One of `--gene` or `--file` |
| `--genome-build` | `GRCh38` (default), `GRCh37` | No |
| `--chunk-size` | `1000` (default) | No |
| `--out` | `output/05-gene-to-positions/results.tsv` | No |

### Input file format (for `--file`)

Plain text, one gene symbol per line, `#` lines ignored:
```
# Candidate genes from GWAS
AGT
BRCA1
TP53
```

---

## Outputs

| File | Description |
|---|---|
| `output/05-gene-to-positions/results.tsv` | One row per gene, with region string and Ensembl metadata |

**Key columns:**

| Column | Example value | Notes |
|---|---|---|
| `gene` | `AGT` | Input gene symbol |
| `ensembl_id` | `ENSG00000135744` | Ensembl stable gene ID |
| `chromosome` | `chr1` | chr-prefixed |
| `start` | `230710048` | 1-based (Ensembl convention) |
| `end` | `230723266` | 1-based |
| `region` | `chr1:230710047-230723266` | 0-based start; ready for FILER |
| `strand` | `1` | |
| `assembly` | `GRCh38` | |

> ℹ️ **Coordinate conventions:** Ensembl returns 1-based positions. The `region` column
> uses a 0-based start (`start - 1`) to match the UCSC/FILER `chrN:start-end` format
> expected by Recipes 2, 3, and 10.

> ℹ️ **Gene body size:** gene regions can span tens to hundreds of kilobases (e.g. BRCA1
> spans ~80kb), so Recipe 3 may return significantly more overlapping tracks than a
> point variant query would. Consider using `--limit` in Recipe 2 to cap the number of
> tracks queried downstream.

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

Single gene lookup:
```bash
python src/scripts/python/gene_to_positions.py \
  --gene AGT \
  --genome-build GRCh38 \
  --out output/05-gene-to-positions/results.tsv
```

Batch lookup from a file:
```bash
python src/scripts/python/gene_to_positions.py \
  --file input/genes.txt \
  --genome-build GRCh38 \
  --out output/05-gene-to-positions/results.tsv
```

Use GRCh37 (hg19) coordinates:
```bash
python src/scripts/python/gene_to_positions.py \
  --file input/genes.txt \
  --genome-build GRCh37 \
  --out output/05-gene-to-positions/results.tsv
```

### Full script: [`src/scripts/python/gene_to_positions.py`](../../src/scripts/python/gene_to_positions.py)

---

## Chaining into downstream recipes

The `region` column in the output is formatted as `chrN:start-end` and can be passed
directly to Recipe 3 or Recipe 10.

### Pipe a single gene into Recipe 3

```bash
# Step 1: resolve gene symbol to region
python src/scripts/python/gene_to_positions.py \
  --gene AGT \
  --out output/05-gene-to-positions/results.tsv

# Step 2: extract the region string and pass to Recipe 3
REGION=$(awk 'NR==2 {print $6}' output/05-gene-to-positions/results.tsv)

python src/scripts/python/filer_coordinate_search.py \
  --region "$REGION" \
  --genome-build hg38 \
  --count-only 0 \
  --full-metadata \
  --out output/03-coordinate-search/search_results.tsv
```

### Pipe a single gene into Recipe 2

```bash
# Pipeline — gene → Recipe 3 → Recipe 2 (coordinate search → overlapping intervals):
# Step 1: resolve gene symbol to region
python src/scripts/python/gene_to_positions.py \
  --gene AGT \
  --out output/05-gene-to-positions/results.tsv

# Step 2: extract the region string and pass to Recipe 3
REGION=$(awk 'NR==2 {print $6}' output/05-gene-to-positions/results.tsv)

python src/scripts/python/filer_coordinate_search.py \
  --region "$REGION" \
  --genome-build hg38 \
  --count-only 0 \
  --out output/03-coordinate-search/search_results.tsv

# Step 3: extract overlapping intervals from discovered tracks. Be warned that this may not finish running since gene regions can produce a large number of overlapping tracks.
python src/scripts/python/filer_find_overlaps.py \
  --region "$REGION" \
  --file output/03-coordinate-search/search_results.tsv \
  --id-col identifier \
  --limit 100 \
  --out output/02-track-overlaps/overlaps.tsv
```

### Pipe a batch into Recipe 10

```bash
# Step 1: resolve all gene symbols
python src/scripts/python/gene_to_positions.py \
  --file input/genes.txt \
  --out output/05-gene-to-positions/results.tsv

# Step 2: loop over each resolved region
tail -n +2 output/05-gene-to-positions/results.tsv | cut -f6 | while read REGION; do
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

**Gene symbols vs Ensembl IDs:** this script accepts HGNC gene symbols (e.g. `AGT`,
`BRCA1`). If you have Ensembl IDs instead (e.g. `ENSG00000135744`), use the
`/lookup/id` endpoint directly or convert symbols first via a tool like MyGene.info.

**Upstream / downstream recipes:**

| Recipe | Description |
|---|---|
| [Recipe 3](../03-coordinate-search/README.md) | Find FILER tracks overlapping a region |
| [Recipe 10](../10-filter-then-overlaps/README.md) | Full end-to-end: metadata filter → region search → intervals |