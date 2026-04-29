# Recipe 11 — Filter Tracks by Metadata and Region, Then Install Locally

## What this does

Run **Recipe 3** (coordinate search with a metadata filter) to find tracks overlapping
your locus, then download and index those tracks directly to a local directory using
`filer_install.py install --from-tracks`.

This is the recommended path when you want a **biologically scoped local FILER
installation** — one that contains only the tracks relevant to a specific locus and
biological context — rather than installing the full FILER data archive.

**Use when:** you have a genomic region of interest (e.g. a GWAS hit, a candidate
enhancer) and want a lightweight local copy of the tracks that overlap it, ready
for offline querying with Giggle and tabix.

> **Why Recipe 3 and not Recipe 10?** Recipe 3 with `--full-metadata` already returns
> every column the installer needs (`identifier`, `file_name`, `file_size`,
> `processed_file_md5`, `processed_file_download_url`/`wget_command`) plus
> `num_overlaps`. The installer accepts a `--top N` flag that ranks by `num_overlaps`
> and caps the install size, so you do not need Recipe 10 just to limit the number of
> tracks. Recipe 10 wraps Recipe 3 with an additional Recipe 2 call to fetch
> *interval-level* data (`hitString`), which the installer ignores. **Use Recipe 10 →
> install only when you also want the per-interval `hitString` data preserved alongside
> the install** (e.g. for downstream analysis of which peaks fell where). For pure
> "find and download," Recipe 3 → install is faster and simpler.

---

## How the steps connect

```
filer_coordinate_search.py --count-only 0 --full-metadata
  └─ search_results.tsv   (one row per overlapping track, full metadata)
          │
          ▼
filer_install.py install --from-tracks search_results.tsv
  ├─ deduplicate on identifier
  ├─ check disk space
  ├─ download each .bed.gz / .vcf.gz with size + MD5 verification
  ├─ giggle index  (per directory)
  └─ tabix index   (per file)
          │
          ▼
    Local FILER installation — mirroring the server's directory structure
```

---

## Inputs

### `filer_coordinate_search.py` parameters

| Parameter | Example values | Required |
|---|---|---|
| `--genome-build` | `hg19`, `hg38` | Yes |
| `--region` | `chr1:100000-200000` | Yes |
| `--count-only` | `0` (must be `0` to get download URLs) | Yes |
| `--full-metadata` | _(flag — required for installer)_ | Yes |
| `--out` | `output/03-coordinate-search/search_results.tsv` | No |

#### Metadata filters (applied server-side via jq filterString)

| Parameter | Example values | Notes |
|---|---|---|
| `--assay` | `ATAC-seq`, `ChIP-seq` | |
| `--cell-type` | `CD14+ monocyte` | |
| `--tissue-category` | `Blood`, `Brain` | |
| `--data-source` | `ENCODE`, `Blueprint` | |
| `--track-id` | `NGBLPL2W2SM2WC` | |
| `--filter-string` | `.data_source == "ENCODE" and .assay == "ATAC-seq"` | Raw jq, overrides named filters |

### `filer_install.py install --from-tracks` parameters

| Parameter | Example values | Required |
|---|---|---|
| `--from-tracks` | `output/03-coordinate-search/search_results.tsv` | Yes |
| `--target-dir` | `FILER_data` | Yes |
| `--giggle` | `giggle`, `/usr/local/bin/giggle` | Yes |
| `--tabix` | `tabix`, `/usr/local/bin/tabix` | Yes |
| `--top` | `100` (rank by `num_overlaps` desc, install only top N) | No |
| `--skip-index` | _(flag)_ | No |
| `--keep-going` | _(flag)_ | No |
| `--verbose` | _(flag)_ | No |

---

## Outputs

### Step 1 — `filer_coordinate_search.py`

| File | Description |
|---|---|
| `search_results.tsv` | One row per overlapping track with full metadata (all columns required by the installer plus `num_overlaps`) |

**Key columns in `search_results.tsv`:**

| Column | Example value |
|---|---|
| `identifier` | `NGBLPL2W2SM2WC` |
| `num_overlaps` | `5` |
| `assay` | `ATAC-seq` |
| `tissue_category` | `Blood` |
| `data_source` | `ENCODE` |
| `track_name` | `ENCODE K562 ATAC-seq peaks` |
| `file_name` | `ENCFF006OFA.bed.gz` |
| `file_size` | `5064124` |
| `processed_file_md5` | `3b5013b3bed9eb54...` |
| `wget_command` | `wget https://… -P FILER2/Annotationtracks/…` |
| `processed_file_download_url` | `https://tf.lisanwanglab.org/GADB/…bed.gz` |
| `tabix_index_download` | `wget https://…bed.gz.tbi -P …` |

### Step 2 — `filer_install.py install --from-tracks`

| Output | Description |
|---|---|
| `<target-dir>/…/*.bed.gz` | Downloaded annotation track files, mirroring the server's relative directory structure |
| `<target-dir>/…/giggle_index/` | Giggle index for each directory of `.bed.gz` files |
| `<target-dir>/…/*.bed.gz.tbi` | Tabix index for each individual track file |
| `<target-dir>/track_metadata.tsv` | One row per downloaded track with all metadata columns plus `local_path` |

**Time:**
- `filer_coordinate_search.py`: 5–60 seconds depending on region size and filter breadth.
- `filer_install.py install --from-tracks`: depends on total file size and network speed;
  MD5 verification is performed on every file after download.

---

## Prerequisites

### Python

Python 3.8 or later is required. `filer_install.py` uses only the Python standard
library (`argparse`, `csv`, `hashlib`, `json`, `re`, `shutil`, `subprocess`, `sys`,
`time`, `urllib`), so no additional pip packages are needed beyond what the project
itself declares.

```bash
# Confirm your Python version
python --version   # must be 3.8+

# Create and activate a virtual environment (recommended)
python -m venv venv
source venv/bin/activate   # Windows: venv\Scripts\activate

# Install the project and its declared dependencies
pip install -e .
```

### System utilities

The following standard Unix tools must be present on your `PATH`. They are available
by default on macOS and most Linux distributions; on Windows use WSL.

| Tool | Used by | Notes |
|---|---|---|
| `wget` | `install_filer.sh` and `filer_install.py` (template mode) | Install via `apt install wget` or `brew install wget` |
| `awk` | `install_filer.sh` | Part of GNU coreutils / macOS base system |
| `envsubst` | `install_filer.sh` | Part of `gettext`; install via `apt install gettext` or `brew install gettext` |
| `md5sum` | `install_filer.sh` | Linux: coreutils. macOS: available as `md5`; install GNU coreutils via `brew install coreutils` |
| `bash` 4+ | `install_filer.sh` | macOS ships with bash 3; upgrade via `brew install bash` |

```bash
# Quick check — all should print a path or version
which wget awk envsubst md5sum bash
bash --version   # should be 4.x or 5.x
```

### Bioinformatics tools

#### tabix (htslib)

tabix is used to create per-file indexes (`.tbi`) for `.bed.gz` and `.vcf.gz` tracks.

```bash
# Option A — conda (recommended)
conda install -c bioconda htslib

# Option B — apt (Debian / Ubuntu)
sudo apt install tabix

# Option C — build from source
git clone --recurse-submodules https://github.com/samtools/htslib.git
cd htslib && autoreconf -i && ./configure && make && sudo make install

# Option D - homebrew
brew install htslib

# Option E - install conda first
brew install miniforge
conda init zsh
# restart your terminal, then:
conda install -c bioconda htslib

# Verify
tabix --version
```

#### giggle (FILER fork)

The FILER workflow requires the FILER-specific fork of giggle, **not** the upstream
version. Build it from source:

```bash
git clone https://github.com/pkuksa/FILER_giggle
cd FILER_giggle
make

# On mac, you may need to run
make HTSLIB=/opt/homebrew/opt/htslib

### WARNING: You may need to modify some of the files (in particular the MakeFiles) so that make runs properly

# Add to PATH for the current session
export PATH="$PWD:$PATH"

# To persist across sessions, add the export line to ~/.bashrc or ~/.zshrc

# Verify
giggle --version
```

> **Important:** the standard `giggle` package available through conda/apt is a
> different build and is not compatible with FILER indexes. Always use the FILER fork
> above.

### Full prerequisite check

Run this block to verify everything is in place before starting:

```bash
echo "--- Python ---"
python --version

echo "--- System utilities ---"
for cmd in wget awk envsubst md5sum bash; do
  command -v "$cmd" && $cmd --version 2>&1 | head -1 || echo "MISSING: $cmd"
done

echo "--- Bioinformatics tools ---"
tabix  --version 2>&1 | head -1 || echo "MISSING: tabix"
giggle | head -1 || echo "MISSING: giggle"
```

No authentication is required for public FILER data.

---

## Python

### Example scripts

Basic usage — ENCODE ATAC-seq blood tracks overlapping a locus, then install locally:

```bash
# Step 1: find tracks at the locus matching the biological filter
python src/scripts/python/filer_coordinate_search.py \
  --genome-build hg38 \
  --region "chr1:100000-200000" \
  --assay "ATAC-seq" \
  --tissue-category "Blood" \
  --data-source "ENCODE" \
  --count-only 0 \
  --full-metadata \
  --out output/03-coordinate-search/search_results.tsv

# Step 2: install only those tracks
python src/scripts/python/filer_install.py install \
  --from-tracks output/03-coordinate-search/search_results.tsv \
  --target-dir FILER_data \
  --giggle giggle \
  --tabix tabix \
  --top 100
```

Using a raw jq filter for finer control (e.g. OR logic across data sources):

```bash
python src/scripts/python/filer_coordinate_search.py \
  --genome-build hg38 \
  --region "chr19:44905791-44909393" \
  --filter-string '.data_source == "ENCODE" and .assay == "ATAC-seq" and .life_stage == "Adult"' \
  --count-only 0 \
  --full-metadata \
  --out output/03-coordinate-search/search_results.tsv

python src/scripts/python/filer_install.py install \
  --from-tracks output/03-coordinate-search/search_results.tsv \
  --target-dir FILER_data \
  --giggle giggle \
  --tabix tabix
```

Download only, skip indexing (useful for staging files before indexing separately):

```bash
python src/scripts/python/filer_install.py install \
  --from-tracks output/03-coordinate-search/search_results.tsv \
  --target-dir FILER_data \
  --giggle giggle \
  --tabix tabix \
  --skip-index
```

Continue past individual download failures rather than stopping on the first error:

```bash
python src/scripts/python/filer_install.py install \
  --from-tracks output/03-coordinate-search/search_results.tsv \
  --target-dir FILER_data \
  --giggle giggle \
  --tabix tabix \
  --keep-going \
  --verbose
```

Cap the install at the top 100 most-overlapping tracks (useful when a broad filter
returns many tracks):

```bash
python src/scripts/python/filer_install.py install \
  --from-tracks output/03-coordinate-search/search_results.tsv \
  --target-dir FILER_data \
  --giggle giggle \
  --tabix tabix \
  --top 100
```

The installer ranks rows by `num_overlaps` descending after deduplication, then keeps
the top N. If the input TSV has no `num_overlaps` column (e.g. Recipe 1 output) the
installer takes the first N rows and prints a warning.

#### Alternative: Recipe 10 → install (when you need per-interval `hitString` data)

Use Recipe 10 instead of Recipe 3 when you want the per-interval `hitString` data
preserved alongside the install for downstream analysis (e.g. of which peaks fell
where). For ranking/top-N alone, prefer the installer's `--top` flag above — it
avoids the extra Recipe 2 API call.

```bash
# Step 1: Recipe 10 fetches both ranking and per-interval data
python src/scripts/python/filer_filter_then_overlaps.py \
  --genome-build hg38 \
  --region "chr1:100000-200000" \
  --assay "ATAC-seq" \
  --tissue-category "Blood" \
  --data-source "ENCODE" \
  --top 100 \
  --out output/10-filter-then-overlaps/results.tsv

# Step 2: install — dedupes per-interval rows on identifier automatically
python src/scripts/python/filer_install.py install \
  --from-tracks output/10-filter-then-overlaps/results.tsv \
  --target-dir FILER_data \
  --giggle giggle \
  --tabix tabix
```

### Integration examples

These examples show how to chain other recipes into the Recipe 11 install + query workflow.

#### rsID → install → query (Recipe 4 → 3 → 11)

Start from a GWAS variant, resolve it to a genomic region, find and install overlapping
ATAC-seq tracks, then query the local index:

```bash
# Step 1: resolve rsID to a genomic region (Recipe 4)
python src/scripts/python/rsid_to_positions.py \
  --rsid rs699 \
  --genome-build GRCh38 \
  --out output/4-rsid-to-positions/results.tsv

# Step 2: extract the region string
REGION=$(awk 'NR==2 {print $5}' output/4-rsid-to-positions/results.tsv)
echo "Resolved region: $REGION"

# Step 3: find tracks overlapping the region matching the biological filter (Recipe 3)
python src/scripts/python/filer_coordinate_search.py \
  --genome-build hg38 \
  --region "$REGION" \
  --assay "ATAC-seq" \
  --tissue-category "Blood" \
  --data-source "ENCODE" \
  --count-only 0 \
  --full-metadata \
  --out output/03-coordinate-search/search_results.tsv

# Step 4: install the tracks locally
python src/scripts/python/filer_install.py install \
  --from-tracks output/03-coordinate-search/search_results.tsv \
  --target-dir FILER_data \
  --giggle giggle \
  --tabix tabix

# Step 5: find the giggle index(es) created during install and query them
#         (the exact path depends on which tracks the API returned)
find FILER_data -type d -name giggle_index | while read IDX; do
  echo "=== $IDX ==="
  giggle search -i "$IDX" -r "$REGION" -c
done
```

#### Gene symbol → install → query (Recipe 5 → 3 → 11)

Start from a gene of interest, install overlapping tracks for the gene body, then
drill into a specific track:

```bash
# Step 1: resolve gene symbol to a genomic region (Recipe 5)
python src/scripts/python/gene_to_positions.py \
  --gene BRCA1 \
  --genome-build GRCh38 \
  --out output/05-gene-to-positions/results.tsv

# Step 2: extract the region string
REGION=$(awk 'NR==2 {print $6}' output/05-gene-to-positions/results.tsv)
echo "Resolved region: $REGION"

# Step 3: find tracks overlapping the gene body (Recipe 3)
python src/scripts/python/filer_coordinate_search.py \
  --genome-build hg38 \
  --region "$REGION" \
  --assay "ChIP-seq" \
  --tissue-category "Breast" \
  --count-only 0 \
  --full-metadata \
  --out output/03-coordinate-search/search_results.tsv

# Step 4: install
python src/scripts/python/filer_install.py install \
  --from-tracks output/03-coordinate-search/search_results.tsv \
  --target-dir FILER_data \
  --giggle giggle \
  --tabix tabix

# Step 5: find the giggle index(es) and query them
find FILER_data -type d -name giggle_index | while read IDX; do
  echo "=== $IDX ==="
  giggle search -i "$IDX" -r "$REGION" -c
done

# Step 6: pick a track from the metadata file, then inspect its intervals with tabix
#         (local_path column gives you the exact file location)
TRACK_PATH=$(awk -F'\t' 'NR==2 {print $2}' FILER_data/track_metadata.tsv)
echo "Inspecting: $TRACK_PATH"
tabix "$TRACK_PATH" "$REGION"
```

#### Batch rsIDs → install all into a shared index (Recipe 4 → 3 → 11)

Process a list of GWAS variants, collect results into a shared install directory,
then query across all of them:

```bash
# Step 1: resolve all rsIDs to regions (Recipe 4)
python src/scripts/python/rsid_to_positions.py \
  --file input/rsids.txt \
  --genome-build GRCh38 \
  --out output/4-rsid-to-positions/results.tsv

# Step 2: run Recipe 3 for each region, collecting results
mkdir -p output/11-selective-install/batch
tail -n +2 output/4-rsid-to-positions/results.tsv | while IFS=$'\t' read -r RSID CHROM START END REGION ALLELE STRAND ASM; do
  SAFE_NAME="${REGION//[^a-zA-Z0-9]/_}"
  echo "Processing $RSID → $REGION"
  python src/scripts/python/filer_coordinate_search.py \
    --genome-build hg38 \
    --region "$REGION" \
    --assay "ATAC-seq" \
    --data-source "ENCODE" \
    --count-only 0 \
    --full-metadata \
    --out "output/11-selective-install/batch/${SAFE_NAME}.tsv"
done

# Step 3: merge all batch results into a single TSV
head -1 output/11-selective-install/batch/*.tsv | head -1 > output/11-selective-install/merged.tsv
tail -q -n +2 output/11-selective-install/batch/*.tsv >> output/11-selective-install/merged.tsv

# Step 4: install everything into one shared directory (deduplicates automatically)
python src/scripts/python/filer_install.py install \
  --from-tracks output/11-selective-install/merged.tsv \
  --target-dir FILER_data \
  --giggle giggle \
  --tabix tabix \
  --keep-going

# Step 5: query the shared index(es) with all variant positions at once
tail -n +2 output/4-rsid-to-positions/results.tsv \
  | awk -F'\t' '{printf "%s\t%s\t%s\n", $2, $3-1, $4}' \
  | sort -k1,1 -k2,2n \
  | bgzip > query_variants.bed.gz

find FILER_data -type d -name giggle_index | while read IDX; do
  echo "=== $IDX ==="
  giggle search -i "$IDX" -q query_variants.bed.gz -o -b
done
```

#### Batch genes → install and query (Recipe 5 → 3 → 11)

Same pattern as above but starting from gene symbols:

```bash
# Step 1: resolve all gene symbols to regions (Recipe 5)
python src/scripts/python/gene_to_positions.py \
  --file input/genes.txt \
  --genome-build GRCh38 \
  --out output/05-gene-to-positions/results.tsv

# Step 2: run Recipe 3 for each gene region
mkdir -p output/11-selective-install/batch
tail -n +2 output/05-gene-to-positions/results.tsv | while IFS=$'\t' read -r GENE ENSID CHROM START END REGION STRAND ASM; do
  SAFE_NAME="${GENE}"
  echo "Processing $GENE → $REGION"
  python src/scripts/python/filer_coordinate_search.py \
    --genome-build hg38 \
    --region "$REGION" \
    --assay "DNase-seq" \
    --tissue-category "Brain" \
    --count-only 0 \
    --full-metadata \
    --out "output/11-selective-install/batch/${SAFE_NAME}.tsv"
done

# Step 3: merge and install
head -1 output/11-selective-install/batch/*.tsv | head -1 > output/11-selective-install/merged.tsv
tail -q -n +2 output/11-selective-install/batch/*.tsv >> output/11-selective-install/merged.tsv

python src/scripts/python/filer_install.py install \
  --from-tracks output/11-selective-install/merged.tsv \
  --target-dir FILER_data \
  --giggle giggle \
  --tabix tabix \
  --keep-going

# Step 4: build a query file from the gene regions and search
tail -n +2 output/05-gene-to-positions/results.tsv \
  | awk -F'\t' '{split($6, a, /[:-]/); printf "%s\t%s\t%s\n", a[1], a[2], a[3]}' \
  | sort -k1,1 -k2,2n \
  | bgzip > query_genes.bed.gz

find FILER_data -type d -name giggle_index | while read IDX; do
  echo "=== $IDX ==="
  giggle search -i "$IDX" -q query_genes.bed.gz -o -b
done
```

#### Track discovery → install → query (Recipe 1 → 11)
Start from a metadata search, install all matching tracks, then query the local index:

```bash
# Step 1: search for tracks by metadata (Recipe 1)
python src/scripts/python/filer_install.py search \
  --genome-build hg38 \
  --assay "ATAC-seq" \
  --tissue-category "Blood" \
  --data-source "ENCODE" \
  --out output/01-track-discovery/tracks.tsv

# Step 2: install those tracks locally
python src/scripts/python/filer_install.py install \
  --from-tracks output/01-track-discovery/tracks.tsv \
  --target-dir FILER_data \
  --giggle giggle \
  --tabix tabix

# Step 3: find the giggle index(es) and query them
find FILER_data -type d -name giggle_index | while read IDX; do
  echo "=== $IDX ==="
  giggle search -i "$IDX" -r "chr1:100000-200000" -c
done

# Step 4: inspect a specific track with tabix
TRACK_PATH=$(awk -F'\t' 'NR==2 {print $2}' FILER_data/track_metadata.tsv)
echo "Inspecting: $TRACK_PATH"
tabix "$TRACK_PATH" "chr1:100000-200000"
```

### Full scripts

- [`filer_coordinate_search.py`](../../src/scripts/python/filer_coordinate_search.py) — Recipe 3 (default upstream for this recipe)
- [`filer_filter_then_overlaps.py`](../../src/scripts/python/filer_filter_then_overlaps.py) — Recipe 10 (alternative upstream when ranking or interval data is needed)
- [`filer_install.py`](../../src/scripts/python/filer_install.py) — Unified FILER CLI (`search` and `install` subcommands)

---

## Notes

**Supported TSV sources:** `--from-tracks` accepts Recipe 3 (with `--count-only 0
--full-metadata`), Recipe 10, and Recipe 1 outputs. The required columns are
`identifier`, `file_name`, `file_size`, `processed_file_md5`, and either `wget_command`
or `processed_file_download_url` — all present in Recipe 3 output by default and in
Recipe 10 output (which embeds the same Recipe 3 columns).

**Directory structure:** downloaded files mirror the server's relative path (e.g.
`FILER2/Annotationtracks/ENCODE/data/ATAC-seq/narrowpeak/hg38/`), so a selective
installation is structurally compatible with a full FILER installation in the same
target directory.

**Deduplication:** Recipe 3 output is one row per track, so no dedup is needed.
Recipe 10 output is one row per overlapping *interval*, so a single track file appears
in multiple rows; the install step automatically deduplicates on `identifier` before
downloading, so each file is downloaded exactly once either way.

**MD5 verification:** every file is verified against `processed_file_md5` after
download. If verification fails, the file is deleted and reported as an error. Use
`--keep-going` to continue with remaining tracks rather than aborting.

**Re-running:** if a file already exists at the expected path with the correct size and
MD5, it is skipped without re-downloading. This makes it safe to re-run the install
step after a partial failure.

**Upstream recipes:** see [Recipe 3](./03-coordinate-search.md) for the default upstream
step, [Recipe 10](./10-filter-then-overlaps.md) for ranking/interval-data variants, and
[Recipe 1](./01-track-discovery.md) for the metadata-only path.

---

## Querying Downloaded FILER Tracks with Giggle and Tabix

After installing tracks with `filer_install.py`, your `FILER_data/` directory
will look something like this:

```
FILER_data/
├── GADB/Annotationtracks/ENCODE/data/ATAC-seq/narrowpeak/hg38/
│   ├── giggle_index/          ← Giggle index built over all .bed.gz in this directory
│   ├── ENCFF008ECD.bed.gz     ← Downloaded track (bgzipped BED)
│   ├── ENCFF008ECD.bed.gz.tbi ← Tabix index for that track
│   ├── ENCFF124SPE.bed.gz
│   ├── ENCFF124SPE.bed.gz.tbi
│   └── ...
└── track_metadata.tsv         ← One row per track with identifier, local_path, and all metadata
```

The `giggle_index/` directory and the `.tbi` files are created automatically
during installation.  The `track_metadata.tsv` file maps every track identifier
to its local path and full metadata (assay, cell type, tissue, data source,
etc.), which is useful for interpreting results from either tool.

---

### Giggle — search across many tracks at once

Giggle is designed for asking "which of my downloaded tracks overlap a region
(or set of regions)?" It operates on the `giggle_index/` directory.

> **Note:** For the following queries, please replace the path of the giggle_index folder if it differs

#### List indexed files

To see which tracks are in the index:

```bash
giggle search \
  -i FILER_data/GADB/Annotationtracks/ENCODE/data/ATAC-seq/narrowpeak/hg38/giggle_index \
  -l
```

This prints one filename per line — each corresponds to a downloaded `.bed.gz`
track.

#### Search a single region

To find which tracks overlap a specific genomic region:

```bash
giggle search \
  -i FILER_data/GADB/Annotationtracks/ENCODE/data/ATAC-seq/narrowpeak/hg38/giggle_index \
  -r "chr1:100000-200000"
```

This returns overlap counts per indexed file.  To also see the actual
overlapping intervals in BED format, add `-o -b`:

```bash
giggle search \
  -i FILER_data/GADB/Annotationtracks/ENCODE/data/ATAC-seq/narrowpeak/hg38/giggle_index \
  -r "chr1:100000-200000" \
  -o -b
```

#### Search multiple regions with a query file

If you have a BED file of regions you want to test (e.g. a set of GWAS loci or
peaks from your own experiment), pass it with `-q`:

```bash
giggle search \
  -i FILER_data/GADB/Annotationtracks/ENCODE/data/ATAC-seq/narrowpeak/hg38/giggle_index \
  -q my_regions.bed.gz \
  -c
```

> **Note:** The query file must be bgzipped and sorted (`sort -k1,1 -k2,2n
> my_regions.bed | bgzip > my_regions.bed.gz`).

The `-c` flag gives overlap counts per indexed file across all query regions.
Replace `-c` with `-o -b` to get the per-record overlapping intervals in BED
format instead.

#### Filter results by filename pattern

If you only care about results from specific tracks, use `-f` with a regex:

```bash
giggle search \
  -i FILER_data/GADB/Annotationtracks/ENCODE/data/ATAC-seq/narrowpeak/hg38/giggle_index \
  -r "chr1:100000-200000" \
  -o -b \
  -f "ENCFF008ECD"
```

#### Searching multiple regions inline (CSV) — caveat

Giggle accepts comma-separated regions via `-r`:

```bash
giggle search \
  -i FILER_data/GADB/Annotationtracks/ENCODE/data/ATAC-seq/narrowpeak/hg38/giggle_index \
  -r "chr1:100000-200000,chr2:500000-600000,chr3:1000000-1100000" \
  -o -b
```

> **Caveat:** When using `-r` with multiple regions, the output's first column
> (which is supposed to identify which query region produced each hit) will
> show `(null)` instead of the actual region.  This makes it impossible to
> tell which hits came from which query region.
>
> **For multi-region searches, prefer a query file (`-q`) instead** — it
> correctly labels each hit with the originating query region.  See
> [Search multiple regions with a query file](#search-multiple-regions-with-a-query-file)
> above.

If you have a few regions and want a quick one-liner without creating a file:

```bash
printf "chr1\t100000\t200000\nchr2\t500000\t600000\nchr3\t1000000\t1100000\n" \
  | sort -k1,1 -k2,2n \
  | bgzip > query_regions.bed.gz

giggle search \
  -i FILER_data/GADB/Annotationtracks/ENCODE/data/ATAC-seq/narrowpeak/hg38/giggle_index \
  -q query_regions.bed.gz \
  -o -b
```

---

### Tabix — query a single track by region

Tabix is designed for fast random access into an individual `.bed.gz` file.
Use it when you already know which track you're interested in and want to
extract intervals from a specific region.

> **Note:** For the following queries, please replace the path of the file if it differs

#### Query a region

```bash
tabix FILER_data/GADB/Annotationtracks/ENCODE/data/ATAC-seq/narrowpeak/hg38/ENCFF008ECD.bed.gz \
  chr1:100000-200000
```

This prints all intervals in that file overlapping `chr1:100000-200000`, one
per line in BED format.

#### Query multiple regions

You can pass several regions on the command line:

```bash
tabix FILER_data/GADB/Annotationtracks/ENCODE/data/ATAC-seq/narrowpeak/hg38/ENCFF008ECD.bed.gz \
  chr1:100000-200000 chr2:500000-600000
```

Or pipe from a BED file of regions:

```bash
tabix -R my_regions.bed \
  FILER_data/GADB/Annotationtracks/ENCODE/data/ATAC-seq/narrowpeak/hg38/ENCFF008ECD.bed.gz
```

#### View the full file header

Some tracks include a header line.  To see it:

```bash
tabix -H FILER_data/GADB/Annotationtracks/ENCODE/data/ATAC-seq/narrowpeak/hg38/ENCFF008ECD.bed.gz
```

---

### Typical workflow: Giggle to discover, Tabix to inspect

A common pattern is to use Giggle first to find which tracks are interesting,
then drill into specific ones with Tabix:

```bash
# 1. Which tracks overlap my region?
giggle search \
  -i FILER_data/GADB/Annotationtracks/ENCODE/data/ATAC-seq/narrowpeak/hg38/giggle_index \
  -r "chr1:100000-200000" \
  -c

# 2. Look up the top hit in track_metadata.tsv to understand what it is
#    (e.g. grep for the filename to find the assay, cell type, etc.)
grep "ENCFF008ECD" FILER_data/track_metadata.tsv

# 3. Pull the actual intervals from that track
tabix FILER_data/GADB/Annotationtracks/ENCODE/data/ATAC-seq/narrowpeak/hg38/ENCFF008ECD.bed.gz \
  chr1:100000-200000
```

---

### Quick reference

| Task | Tool | Key flags |
|---|---|---|
| Which tracks overlap a region? | `giggle search -i <index> -r <region>` | `-c` counts, `-o -b` BED output |
| Which tracks overlap regions in a file? | `giggle search -i <index> -q <file.bed.gz>` | `-c` counts, `-o -b` BED output |
| List files in an index | `giggle search -i <index> -l` | |
| Filter by track name | `giggle search ... -f <regex>` | |
| Extract intervals from one track | `tabix <file.bed.gz> <region>` | |
| Extract using a BED of regions | `tabix -R <regions.bed> <file.bed.gz>` | |
| View file header | `tabix -H <file.bed.gz>` | |

---

## To-Do
Provide instructions for modifying the GitHub for FILER_giggle