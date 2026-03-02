# Recipe 1 — List Tracks Matching Metadata Filters

## What this does

Query the FILER metadata endpoint with standardized filters (cell type, assay, tissue category,
data source) and return a table of matching tracks saved as `tracks.tsv` (and optionally
`tracks.json`).

**Use when:** you want to assemble a relevant set of tracks before querying overlaps. Almost
every downstream workflow starts here.

---

## Inputs

| Parameter | PHP param | Example values | Required |
|---|---|---|---|
| `genome_build` | `genomeBuild` | `hg19`, `hg38` | Yes |
| `assay` | `assayType` | `ATAC-seq`, `ChIP-seq`, `DNase-seq` | No |
| `cell_type` | `cellType` | `CD14+ monocyte`, `T cell` | No |
| `tissue_category` | `tissueCategory` | `blood`, `brain` | No |
| `data_source` | `dataSource` | `ENCODE`, `Roadmap`, `DASHR2` | No |
| `track_id` | `trackID` | `ENCFF001ABC` | No |
| `output_format` | `outputFormat` | `json` (default), `tsv` | No |

At least one filter beyond `genome_build` is recommended to keep results manageable.

## Outputs

| File | Description |
|---|---|
| `tracks.tsv` | Tab-separated table of matching tracks |
| `tracks.json` | Same content as a JSON array (optional) |

**Key columns returned by the API:**

| Column | Example value |
|---|---|
| `identifier` | `NGBLPL2W2SM2WC` |
| `genome_build` | `hg38` |
| `assay` | `WGB-Seq`, `ATAC-seq` |
| `cell_type` | `B cell`, `CD14+ monocyte` |
| `biosample_type` | `Primary cell` |
| `tissue_category` | `Blood` |
| `system_category` | `Cardiovascular`, `Immune` |
| `life_stage` | `Child`, `Adult` |
| `data_source` | `Blueprint`, `ENCODE` |
| `data_category` | `Methylation`, `Chromatin accessibility` |
| `classification` | `WGB-Seq peaks` |
| `output_type` | `peaks` |
| `track_name` | `Blueprint B cell WGB-Seq peaks (bed4) [Life stage: Child]` |
| `processed_file_download_url` | `https://tf.lisanwanglab.org/GADB/…bed.gz` |
| `tabix_file_url` | `https://tf.lisanwanglab.org/GADB/…bed.gz.tbi` |

Additional columns also returned: `file_name`, `number_of_intervals`, `bp_covered`,
`file_size`, `file_format`, `encode_experiment_id`, `biological_replicate`,
`technical_replicate`, `antibody`, `downloaded_date`, `release_date`,
`date_added_to_filer`, `processed_file_md5`, `link_out_url`, `track_description`,
`biosamples_term_id`.

**Time:** Typically 2–10 seconds depending on filter breadth.

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

This script demonstrates all parameters that you can use, with the exception of track-id since it would return at most 1 track back. Feel free to add, remove, or change any lines to suit your needs.

```bash
python src/scripts/python/filer_search_tracks.py \
  --genome-build hg38 \
  --assay "ATAC-seq" \
  --data-source "ENCODE" \
  --cell-type "K562" \
  --tissue-category "Blood" \
  --out output/01-track-discovery/tracks.tsv \
  --json
```

This searches for tracks that pertain to brain tissues and use Hi-C assays.

```bash
python src/scripts/python/filer_search_tracks.py \
  --genome-build hg38 \
  --assay "Hi-C" \
  --tissue-category "Brain" \
  --out output/01-track-discovery/tracks.tsv \
  --json
```

This searches up a track by its ID.
```bash
python src/scripts/python/filer_search_tracks.py \
  --genome-build hg38 \
  --track-id "NGBLPL2W2SM2WC" \
  --out output/01-track-discovery/tracks.tsv \
  --json
```

### Full script: [`src/scripts/python/filer_search_tracks.py`](../../src/scripts/python/filer_search_tracks.py).

---

## Bash (minimal) 
Make sure that jq is installed on your computer.

```bash
BASE="https://tf.lisanwanglab.org/FILER2"

curl -sG "${BASE}/get_metadata.php" \
  --data-urlencode "genomeBuild=hg38" \
  --data-urlencode "assayType=ATAC-seq" \
  --data-urlencode "dataSource=ENCODE" \
  --data-urlencode "outputFormat=json" \
  > output/01-track-discovery/tracks.json

jq -r '
  ["identifier","genome_build","assay","cell_type","biosample_type",
   "tissue_category","life_stage","data_source","track_name","processed_file_download_url","tabix_file_url"],
  (.[] | [
    .identifier,.genome_build,.assay,.cell_type,.biosample_type,
    .tissue_category,.life_stage,.data_source,.track_name,.processed_file_download_url,.tabix_file_url
  ])
  | @tsv
' output/01-track-discovery/tracks.json > output/01-track-discovery/tracks.tsv
```

---

## Advanced: raw `filterString` (jq syntax)

For power users, you can pass a raw jq filter string directly. This bypasses the named
parameters and gives you full flexibility:

```bash
curl -sG "${BASE}/get_metadata.php" \
  --data-urlencode 'filterString=.data_source == "ENCODE" and .assay == "ATAC-seq"' \
  --data-urlencode "genomeBuild=hg38" \
  --data-urlencode "outputFormat=json" \
  > output/01-track-discovery/tracks.json
```

```python
params = {
    "genomeBuild":    "hg38",
    "filterString":   '.data_source == "ENCODE" and .assay == "ATAC-seq"',
    "outputFormat":   "json",
}
```