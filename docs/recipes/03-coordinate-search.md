# Recipe 3 — Find Tracks Overlapping a Genomic Region

## What this does

Query the FILER coordinate-overlap endpoint with a genomic region and optional metadata
filters to return all tracks whose intervals intersect that region, saved as
`search_results.tsv`.

**Use when:** you have a specific locus of interest (e.g., a GWAS hit, a candidate enhancer)
and want to discover which FILER tracks have signal there. Pair with Recipe 1 to pre-filter
the track universe before querying.

---

## Inputs

| Parameter | PHP param | Example values | Required |
|---|---|---|---|
| `region` | `region` | `chr1:100000-200000` | Yes |
| `genome_build` | `genomeBuild` | `hg19`, `hg38` | Yes |
| `filter_string` | `filterString` | `.data_source == "ENCODE"` | No |
| `full_metadata` | `fullMetadata` | `0` (default), `1` | No |
| `count_only` | `countOnly` | `1` (default), `0` | No |
| `output_format` | `outputFormat` | `json` (default), `html` | No |

**`count_only` and `full_metadata` interaction:**

| `count_only` | `full_metadata` | What you get |
|---|---|---|
| `1` | `0` | `identifier` + `num_overlaps` only (default, fastest) |
| `0` | `0` | `identifier` + `processed_file_download_url` + `num_overlaps` |
| `0` | `1` | All 36 metadata columns + `num_overlaps` |

> ⚠️ `filterString` requires `full_metadata=1` to work. The filter is applied server-side
> after metadata is fetched, so fields like `.data_source` must be present in the response
> to match against. Without `full_metadata=1`, every filter evaluates to `null == "value"`
> and returns empty results. The Python script handles this automatically.

---

## Outputs

| File | Description |
|---|---|
| `search_results.tsv` | Tab-separated table of overlapping tracks |

**Key columns returned by the API (when `count_only=0, full_metadata=1`):**

| Column | Example value |
|---|---|
| `identifier` | `NGBLPL2W2SM2WC` |
| `genome_build` | `hg38` |
| `assay` | `ATAC-seq` |
| `cell_type` | `K562` |
| `tissue_category` | `Blood` |
| `system_category` | `Immune` |
| `life_stage` | `Adult` |
| `data_source` | `ENCODE` |
| `classification` | `ATAC-seq peaks` |
| `track_name` | `ENCODE K562 ATAC-seq peaks` |
| `processed_file_download_url` | `https://tf.lisanwanglab.org/GADB/…bed.gz` |
| `num_overlaps` | `3` |

When `count_only=1`, only `identifier` and `num_overlaps` are returned.

**Time:** Typically 5–60 seconds depending on region size and filter breadth.

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

## Python (recommended)

### Example scripts

Search all hg38 tracks overlapping a region, returning counts only (default):
```bash
python src/scripts/python/filer_coordinate_search.py \
  --region "chr1:100000-200000" \
  --genome-build hg38 \
  --out output/03-coordinate-search/search_results.tsv
```

Return all metadata fields for overlapping tracks:
```bash
python src/scripts/python/filer_coordinate_search.py \
  --region "chr1:100000-200000" \
  --genome-build hg38 \
  --count-only 0 \
  --full-metadata \
  --out output/03-coordinate-search/search_results.tsv
```

Filter to ENCODE ATAC-seq tracks from adult blood (single data source):
```bash
python src/scripts/python/filer_coordinate_search.py \
  --region "chr1:100000-200000" \
  --genome-build hg38 \
  --filter-string '.data_source == "ENCODE" and .assay == "ATAC-seq" and .tissue_category == "Blood" and .life_stage == "Adult"' \
  --count-only 0 \
  --out output/03-coordinate-search/search_results.tsv
```

Filter across multiple data sources using OR:
```bash
python src/scripts/python/filer_coordinate_search.py \
  --region "chr1:100000-200000" \
  --genome-build hg38 \
  --filter-string '(.data_source == "ENCODE" or .data_source == "Blueprint") and .output_type == "peaks"' \
  --count-only 0 \
  --out output/03-coordinate-search/search_results.tsv
```

### Full script: [`src/scripts/python/filer_coordinate_search.py`](../../src/scripts/python/filer_coordinate_search.py)

---

## Bash (minimal)

Make sure that `jq` is installed on your computer.
```bash
BASE="https://tf.lisanwanglab.org/FILER2"

curl -sG "${BASE}/get_overlapping_tracks_by_coord.php" \
  --data-urlencode "region=chr1:100000-200000" \
  --data-urlencode "genomeBuild=hg38" \
  --data-urlencode "outputFormat=json" \
  --data-urlencode "countOnly=1" \
  > output/03-coordinate-search/search_results.json

jq -r '
  ["identifier","num_overlaps"],
  (.[] | [.identifier, .num_overlaps])
  | @tsv
' output/03-coordinate-search/search_results.json \
  > output/03-coordinate-search/search_results.tsv
```

For full metadata, add `fullMetadata=1` and `countOnly=0` and extend the `jq` field list:
```bash
curl -sG "${BASE}/get_overlapping_tracks_by_coord.php" \
  --data-urlencode "region=chr1:100000-200000" \
  --data-urlencode "genomeBuild=hg38" \
  --data-urlencode "outputFormat=json" \
  --data-urlencode "countOnly=0" \
  --data-urlencode "fullMetadata=1" \
  > output/03-coordinate-search/search_results.json

jq -r '
  ["identifier","data_source","assay","cell_type","tissue_category","life_stage","track_name","num_overlaps"],
  (.[] | [.identifier,.data_source,.assay,.cell_type,.tissue_category,.life_stage,.track_name,.num_overlaps])
  | @tsv
' output/03-coordinate-search/search_results.json \
  > output/03-coordinate-search/search_results.tsv
```

---

## Advanced: metadata filtering with `filterString`

The `filterString` parameter accepts a jq-style boolean expression evaluated server-side
against each track's metadata. **Requires `fullMetadata=1`** — without it, metadata fields
are not present in the response and every filter returns empty results.

Single condition:
```bash
curl -sG "${BASE}/get_overlapping_tracks_by_coord.php" \
  --data-urlencode "region=chr1:100000-200000" \
  --data-urlencode "genomeBuild=hg38" \
  --data-urlencode 'filterString=.data_source == "ENCODE"' \
  --data-urlencode "countOnly=0" \
  --data-urlencode "fullMetadata=1" \
  --data-urlencode "outputFormat=json"
```

Multiple conditions:
```bash
curl -sG "${BASE}/get_overlapping_tracks_by_coord.php" \
  --data-urlencode "region=chr1:100000-200000" \
  --data-urlencode "genomeBuild=hg38" \
  --data-urlencode 'filterString=.data_source == "ENCODE" and .assay == "ATAC-seq" and .tissue_category == "Blood"' \
  --data-urlencode "countOnly=0" \
  --data-urlencode "fullMetadata=1" \
  --data-urlencode "outputFormat=json"
```
```python
payload = {
    "region":       "chr1:100000-200000",
    "genomeBuild":  "hg38",
    "filterString": '.data_source == "ENCODE" and .tissue_category == "Blood"',
    "outputFormat": "json",
    "fullMetadata": 1,
    "countOnly":    0,
}
```

Pass `filterString=.` (the default) to apply no filter and return all overlapping tracks.