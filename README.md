# FILER-workflows

Python scripts, notebooks, and a thin client library for working with [FILER](https://tf.lisanwanglab.org/FILER2) — a functional genomics database of regulatory annotation tracks from ENCODE, Blueprint, Roadmap, and other consortia.

---

## What's in this repo

| Path | Description |
|---|---|
| `src/filerpy/` | `filerpy` — thin Python client for the FILER HTTP API |
| `src/scripts/python/` | Standalone CLI scripts (one per recipe) |
| `notebooks/python/` | Jupyter notebooks mirroring the CLI recipes |
| `docs/recipes/` | Detailed documentation for each recipe |
| `input/` | Sample input files (`rsids.txt`, `genes.txt`) |
| `res/` | FILER schema reference (`filer-schema.tsv`) |

---

## Installation

Python 3.9+ is required.

```bash
# Create and activate a virtual environment (recommended)
python3 -m venv venv
source venv/bin/activate        # Windows: venv\Scripts\activate

# Install the project and all dependencies
pip install -e ".[dev]"
```

No API key or authentication is required for public FILER data.

---

## Recipes

Each recipe is a self-contained workflow. They are designed to be run independently or chained together.

| Recipe | Script | Description |
|---|---|---|
| [01 — Track discovery](docs/recipes/01-track-discovery.md) | `filer_search_tracks.py` | Search tracks by assay, cell type, tissue, data source |
| [02 — Track overlaps](docs/recipes/02-track-overlaps.md) | `filer_find_overlaps.py` | Fetch overlapping intervals from a set of track IDs |
| [03 — Coordinate search](docs/recipes/03-coordinate-search.md) | `filer_coordinate_search.py` | Find all tracks with signal in a genomic region |
| [04 — rsID → positions](docs/recipes/04-rsid-to-positions.md) | `rsid_to_positions.py` | Convert dbSNP rsIDs to `chrN:start-end` via Ensembl |
| [05 — Gene → positions](docs/recipes/05-gene-to-positions.md) | `gene_to_positions.py` | Convert gene symbols to `chrN:start-end` via Ensembl |
| [10 — Filter then overlaps](docs/recipes/10-filter-then-overlaps.md) | `filer_filter_then_overlaps.py` | End-to-end: metadata filter → coordinate search → interval extraction |
| [11 — Selective install](docs/recipes/11-filer-selective-install.md) | `filer_install.py` | Run various recipes then download and index matching tracks locally |

### Common starting points

**Find ATAC-seq blood tracks in hg38:**
```bash
python src/scripts/python/filer_search_tracks.py \
  --genome-build hg38 \
  --assay "ATAC-seq" \
  --tissue-category "Blood" \
  --out output/01-track-discovery/tracks.tsv
```

**Find all tracks with signal at a locus:**
```bash
python src/scripts/python/filer_coordinate_search.py \
  --region "chr1:100000-200000" \
  --genome-build hg38 \
  --count-only 0 \
  --full-metadata \
  --out output/03-coordinate-search/results.tsv
```

**End-to-end from a GWAS variant:**
```bash
# Resolve rsID to coordinates
python src/scripts/python/rsid_to_positions.py \
  --rsid rs699 \
  --out output/04-rsid-to-positions/results.tsv

# Extract the region and run the full workflow
REGION=$(awk 'NR==2 {print $5}' output/04-rsid-to-positions/results.tsv)
python src/scripts/python/filer_filter_then_overlaps.py \
  --genome-build hg38 \
  --region "$REGION" \
  --assay "ATAC-seq" \
  --tissue-category "Blood" \
  --top 100 \
  --out output/10-filter-then-overlaps/results.tsv
```

---

## filerpy library

`filerpy` provides a Python API that wraps the same FILER endpoints used by the CLI scripts.

```python
from filerpy.client import search_tracks, get_overlapping_tracks, fetch_overlaps
from filerpy.trackset import make_trackset, save_trackset, load_trackset

# Search for tracks
df = search_tracks("hg38", assayType="ATAC-seq", cellType="CD14+ monocyte")

# Find tracks overlapping a region
hits = get_overlapping_tracks("chr1:100000-200000", "hg38", fullMetadata=1)

# Fetch interval data from specific tracks
overlaps = fetch_overlaps("chr19:44905791-44909393", df["identifier"].tolist())

# Save a pinned, versioned track set for reproducibility
ts = make_trackset(df, "hg38", {"assayType": "ATAC-seq"}, name="atac_blood")
save_trackset(ts, out_dir="output/tracksets", df=df)
```

### Key modules

| Module | Purpose |
|---|---|
| `filerpy.client` | HTTP requests to FILER endpoints (`search_tracks`, `get_overlapping_tracks`, `fetch_overlaps`, `get_metadata`) |
| `filerpy.trackset` | Read/write Track Set manifests (`.trackset.json` / `.trackset.tsv`) for reproducible analyses |

---

## Notebooks

Interactive Jupyter notebooks are available in `notebooks/python/` for Recipes 01, 02, 03, and 10. They mirror the CLI scripts and can be used as a starting point for custom analyses.

---

## Recipe chaining

The recipes are designed to compose. The typical patterns are:

```
rsID / gene symbol
       │
  Recipe 4/5 (Ensembl lookup → region)
       │
  Recipe 3 (coordinate search → overlapping track IDs)
       │
  Recipe 1 (metadata filter → biological context)
       │
  Recipe 2 (fetch actual intervals)
       │
  Recipe 10 (all of the above in one script)
       │
  Recipe 11 (+ download and index locally)
```

See the individual recipe docs in `docs/recipes/` for detailed parameter tables, output column descriptions, and copy-paste examples.

---

## Prerequisites for Recipe 11 (local install)

Recipe 11 downloads tracks and builds Giggle and tabix indexes locally. It requires:

- **tabix** (htslib): `conda install -c bioconda htslib` or `brew install htslib`
- **giggle** (FILER fork — the standard conda/apt version is not compatible):

```bash
git clone https://github.com/pkuksa/FILER_giggle
cd FILER_giggle
make   # on macOS: make HTSLIB=/opt/homebrew/opt/htslib
export PATH="$PWD:$PATH"
```

### macOS Makefile modifications

If building on macOS with Homebrew-installed htslib fails due to missing dependencies or linker errors, replace the two Makefiles below with the full contents shown here. They add an `HTSLIB` switch so you can point at the Homebrew install (`make HTSLIB=/opt/homebrew/opt/htslib`) and link the extra dependencies Homebrew's htslib requires.

**Replace `FILER_giggle/Makefile` with:**

```makefile
BIN=bin
OBJ=obj

ifdef HTSLIB
all:
else
all: htslib
endif
	@mkdir -p $(OBJ)
	@mkdir -p $(BIN)
	cd src; $(MAKE)

htslib:
	$(shell cd lib/htslib && autoreconf --install && autoreconf) # for Mac OS may need to use specific autoconf version, e.g., 2.69
	cd lib/htslib; ./configure --disable-bz2 --disable-lzma --disable-libcurl --host=x86_64 # may need to explicitly specify host if configure fails, e.g., --host=x86_64-apple-darwin20.3.0
	$(MAKE) -C lib/htslib

server:
	@mkdir -p $(OBJ)
	@mkdir -p $(BIN)
	cd src; $(MAKE) server

clean:
	rm -rf $(BIN)/*
	rm -rf $(OBJ)/*
	cd lib/htslib && $(MAKE) clean
```

**Replace `FILER_giggle/src/Makefile` with:**

```makefile
BIN=../bin
OBJ=../obj
LIBD=../lib

ifdef HTSLIB
HTS_ROOT = $(HTSLIB)
HTS_LIB = $(HTS_ROOT)/lib/libhts.a
HTS_INCLUDE = -I$(HTS_ROOT)/include
else
HTS_ROOT = $(LIBD)/htslib
HTS_LIB = $(HTS_ROOT)/libhts.a
HTS_INCLUDE = -I$(HTS_ROOT)
endif

LIBMHD_INCLUDES=$(HOME)/usr/local/include
LIBMHD_LIBS=$(HOME)/usr/local/lib
CFLAGS=-O2 -D_FILE_OFFSET_BITS=64 -Werror -Wuninitialized
#CFLAGS=-g -D_FILE_OFFSET_BITS=64 -Werror -Wuninitialized
CFLAGS+=-DBLOCK_STORE
CFLAGS+=-fcommon
INCLUDES=$(HTS_INCLUDE) -L/opt/homebrew/opt/openssl@3/lib -L/opt/homebrew/opt/libdeflate/lib -L/opt/homebrew/opt/xz/lib
LIBS=-ldl -lz -lm -pthread -lcurl
LIBS+=-lcrypto
LIBS+=-lbz2 -llzma -ldeflate
CC=gcc


LIB=bpt.o \
    disk_store.o \
    cache.o \
    data_reg.o \
    file_read.o \
    giggle_index.o \
    lists.o \
    ll.o \
    util.o \
    timer.o \
    wah.o \
    kfunc.o \
    index.o \
    offset_index.o \
    search.o \
    leaf.o \
    pq.o \
    jsw_avltree.o \
    fastlz.o


.SUFFIXES:

OBJS=$(addprefix $(OBJ)/, $(LIB))

.SECONDARY: $(OBJS)


PROG=giggle \
     api_test \
     offset_idx_lookup

LIST=$(addprefix $(BIN)/, $(PROG))


all: check-env $(LIST) library

server: INCLUDES+=-I$(LIBMHD_INCLUDES) -L$(LIBMHD_LIBS)
server: LIBS+=-lmicrohttpd
server: $(BIN)/server_enrichment

library: $(OBJS)
	ar -cvq $(LIBD)/libgiggle.a $(OBJS)

$(OBJ)/%.o: %.c
	$(CC) -c $(CFLAGS) -o $@ $< \
		$(HTS_INCLUDE)


$(BIN)/%: %.c $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ \
		$(INCLUDES) \
		-DSAMTOOLS=1 \
		$(HTS_LIB) \
		$(LIBS)


check-env:
ifdef TRAVIS_BUILD_DIR
HTS_ROOT=$(TRAVIS_BUILD_DIR)/htslib
endif
```

> **Tab characters required.** Makefile recipe lines must start with a real tab, not spaces. If you copy these blocks from a rendered view that converts tabs, re-indent the recipe lines (the ones under each target) with tabs before running `make`.

See [Recipe 11](docs/recipes/11-filer-selective-install.md) for the full prerequisites checklist and usage examples.
