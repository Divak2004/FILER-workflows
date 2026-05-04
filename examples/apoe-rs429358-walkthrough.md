# Does the APOE-ε4 variant rs429358 sit in active brain chromatin?

## The question

> **The APOE-ε4 variant rs429358 is the strongest common genetic risk factor for
> late-onset Alzheimer's disease. Does it sit inside active regulatory chromatin in
> brain tissue, and if so, which peaks call it?**

rs429358 (chr19, APOE locus) is famous as a missense SNP — Cys→Arg in APOE — but a
secondary question is whether it *also* disrupts a regulatory element. The FILER
pipeline answers this end-to-end without leaving the command line.

## Recipe chain

```
Recipe 4 (rsid_to_positions.py)
  rs429358 → chr19:44908684-44908685
        │
        ▼
Recipe 10 (filer_filter_then_overlaps.py)
  Find brain ATAC-seq / ChIP-seq tracks at that position AND pull
  their peak intervals — one row per overlap, with full track metadata
        │
        ▼
Recipe 11 (filer_install.py install --from-tracks)
  (optional) Install those tracks locally for offline tabix / giggle queries
```

## The commands

### Step 1 — Resolve the variant to a region (Recipe 4)

```bash
python src/scripts/python/rsid_to_positions.py \
  --rsid rs429358 \
  --genome-build GRCh38 \
  --out output/4-rsid-to-positions/results.tsv

REGION=$(awk 'NR==2 {print $5}' output/4-rsid-to-positions/results.tsv)
echo "$REGION"   # → chr19:44908684-44908685
```

### Step 2 — Find brain regulatory tracks AND pull their peak intervals (Recipe 10)

Recipe 10 chains Recipe 3 (coordinate search with metadata filter) and Recipe 2
(interval extraction) into one call. Output is one row per overlapping interval
with full track metadata attached.

```bash
python src/scripts/python/filer_filter_then_overlaps.py \
  --genome-build hg38 \
  --region "$REGION" \
  --filter-string '.tissue_category == "Brain" and (.assay == "ATAC-seq" or .assay == "ChIP-seq")' \
  --top 200 \
  --out output/10-filter-then-overlaps/brain_apoe.tsv
```

### Step 3 (optional) — Install the relevant tracks locally (Recipe 11)

Recipe 11 dedupes the input on `identifier`, so it consumes Recipe 10's
interval-level output directly without any reshape.

```bash
python src/scripts/python/filer_install.py install \
  --from-tracks output/10-filter-then-overlaps/brain_apoe.tsv \
  --target-dir FILER_data \
  --giggle giggle \
  --tabix tabix \
  --top 50
```

After this step, every downloaded `.bed.gz` has a sibling `.tbi`, and each
leaf directory containing `.bed.gz` files has its own `giggle_index/` (so a
50-track install produces several indices, one per data source / peak format —
e.g. RoadMap narrowpeak, RoadMap broadpeak, IHEC bed5, ATACdb bed13, …).
The resolved local path of each track is stored in the `local_path` column
of `FILER_data/track_metadata.tsv`. You can now query the local copy without
touching the network.

### Step 4 (optional) — Query the local indices

#### Tabix — pull peaks overlapping rs429358 from a single track

```bash
# local_path is column 2 of track_metadata.tsv
TRACK=$(awk -F'\t' 'NR==2 {print $2}' FILER_data/track_metadata.tsv)

tabix "$TRACK" chr19:44908684-44908685
```

#### Tabix — sweep every installed track for the variant position

```bash
find FILER_data -name '*.bed.gz' | while read -r f; do
  hits=$(tabix "$f" chr19:44908684-44908685)
  if [ -n "$hits" ]; then
    echo "=== $(basename "$f") ==="
    echo "$hits"
  fi
done
```

#### Giggle — score one region against every installed index

`filer_install.py` creates one `giggle_index/` per leaf folder, so you query
each in turn. `-s` prints per-track overlap counts and Fisher enrichment
(combo score in column 8).

```bash
# Build a one-line query BED, bgzip it (giggle requires bgzipped input)
printf 'chr19\t44908684\t44908685\trs429358\n' \
  | bgzip -c > query.bed.gz

find FILER_data -type d -name giggle_index | while read -r idx; do
  echo "=== $idx ==="
  giggle search -i "$idx" -q query.bed.gz -s
done
```

#### Giggle — return the actual overlapping intervals (`-v`)

```bash
find FILER_data -type d -name giggle_index | while read -r idx; do
  giggle search -i "$idx" -q query.bed.gz -v
done
```

The last column of each row is the source track path, so you can pipe the
combined output into `awk -F'\t' '{print $NF}' | sort -u` to get the list of
tracks whose peaks call rs429358.