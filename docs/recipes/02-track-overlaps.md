# Recipe 2 — Extract Overlaps from Specific Tracks

## What this does

Query the FILER overlaps endpoint with a genomic region and a list of track IDs, and return all data intervals from those tracks that overlap the region, saved as `overlaps.tsv` (and optionally `overlaps.json`).

**Use when:** you have a **small** set of tracks from Recipe 1 and want to pull the actual genomic features (peaks, intervals, etc.) that fall within a region of interest. For context, it takes approximately a minute to process 250 tracks.

---

## Inputs

| Parameter | PHP param | Example values | Required |
|---|---|---|---|
| `region` | `region` | `chr19:44905791-44909393` | Yes |
| `track_ids` | `trackIDs` | `NGBLPL2W2SM2WC,NG123` | Yes (or `--file`) |
| `--file` | _(reads `trackIDs` from a TSV/JSON file)_ | `output/01-track-discovery/tracks.tsv` | Yes (or `--track-ids`) |
| `--id-col` | _(column name in `--file`)_ | `identifier` (default) | No |

The region must be in `chrN:start-end` format (e.g. `chr1:100000-200000`).

Track IDs can be supplied directly via `--track-ids` or loaded from a file produced by
Recipe 1 via `--file`. If both are given, `--track-ids` takes precedence.

## Outputs

| File | Description |
|---|---|
| `overlaps.tsv` | Tab-separated table of overlapping intervals, one row per feature |
| `overlaps.json` | Same content as a JSON array (optional) |

**Output columns from this recipe:**

| Column | Example value |
|---|---|
| `Identifier` | `NGBLPL2W2SM2WC` |
| `queryRegion` | `chr19:44905791-44909393` |
| `hitString` | Feature fields joined by `@@@`, using the API-provided field order for that feature |

The `hitString` field set/order (and thus the number of tokens) varies by track type.
Each token corresponds to one value from the underlying feature record returned by FILER.
See the schema file under the res folder.

**Time:** Varies with the number of tracks and interval density. Expect roughly 1 minute per batch of 250 tracks

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

Query a region against a list of track IDs supplied directly.

```bash
python src/scripts/python/filer_find_overlaps.py \
  --region "chr19:44905791-44909393" \
  --track-ids "NGBLPL2W2SM2WC,NGENCC4CIT5FBQ" \
  --out output/02-track-overlaps/overlaps.tsv \
  --json
```

Query a region against all tracks returned by Recipe 1 (recommended for large searches). Please be warned that larger numbers of tracks take a longer time. 

```bash
python src/scripts/python/filer_find_overlaps.py \
  --region "chr19:44905791-44909393" \
  --file output/01-track-discovery/tracks.tsv \
  --id-col identifier \
  --out output/02-track-overlaps/overlaps.tsv \
  --json
```

Use a smaller chunk size if the server returns errors on batches (default is 250).

```bash
python src/scripts/python/filer_find_overlaps.py \
  --region "chr19:44905791-44909393" \
  --file output/01-track-discovery/tracks.tsv \
  --id-col identifier \
  --out output/02-track-overlaps/overlaps.tsv \
  --chunk-size 100 \
  --json
```

### Full script: [`src/scripts/python/filer_find_overlaps.py`](../../src/scripts/python/filer_find_overlaps.py).

---

## Bash (minimal)
Make sure that `jq` is installed on your computer.

```bash
BASE="https://tf.lisanwanglab.org/FILER2"

curl -sG "${BASE}/get_overlaps.php" \
  --data-urlencode "region=chr19:44905791-44909393" \
  --data-urlencode "trackIDs=NGBLPL2W2SM2WC,NG123" \
  > output/02-track-overlaps/overlaps.json

jq -r '
  ["Identifier","queryRegion","hitString"],
  (.[] | .Identifier as $id | .queryRegion as $qr |
    .features[] | [$id, $qr, (to_entries | map(.value|tostring) | join("@@@"))]
  )
  | @tsv
' output/02-track-overlaps/overlaps.json > output/02-track-overlaps/overlaps.tsv
```

> **Note:** The `curl` approach works for small numbers of track IDs. For large lists
> (hundreds or thousands of IDs) the URL will exceed the server's length limit — use the
> Python script instead, which handles chunking automatically.

---

## Notes

**Chunking:** The server only reads query string parameters (`$_GET`), so both `region` and
`trackIDs` must be passed as GET params. To stay within the server's URL length limit, the
Python script automatically splits large ID lists into batches of `--chunk-size` (default 250) and merges the results.
Use `--chunk-size` to tune this if needed.

**Track IDs not found:** If a track ID does not exist in FILER, the API returns an error
object inside `features` for that track rather than failing the whole request. These rows will
appear in the output with an `ERROR` token inside `hitString` and can be filtered out downstream.

**Upstream dependency:** This recipe is designed to run after Recipe 1. The `--file` +
`--id-col identifier` pattern directly accepts the `tracks.tsv` output from Recipe 1 without
any transformation.