# How active is the CACNA1C gene body across brain ATAC-seq?

## The question

> **CACNA1C is the strongest common GWAS gene shared across schizophrenia and
> bipolar disorder. The risk signals are noncoding and scattered across the
> ~600 kb gene body — so any effect must go through regulatory chromatin
> somewhere inside the locus. Of all brain ATAC-seq tracks in FILER,
> what fraction call at least one peak inside CACNA1C, and which tracks
> have the densest signal?**

For a missense variant ("does it disrupt the protein?") a single overlap is
enough to tell a story. For a noncoding GWAS gene, a single overlap means
nothing — you need a **rate**. A gene that is open in 47/300 brain
ATAC-seq tracks is meaningfully active; 5/300 is background.

This workflow pairs **Recipe 1** (universe → denominator) with **Recipe 3 with
`--full-metadata`** (filtered, region-restricted track list with `num_overlaps`
per track → numerator + ranking).

## Recipe chain

```
Recipe 5 (gene_to_positions.py)
  CACNA1C → chr12:1970771-2697950
        │
        ├──────────────────────────────────────┐
        ▼                                      ▼
Recipe 1 (filer_search_tracks.py)        Recipe 3 (filer_coordinate_search.py)
  Full universe: every brain               Same biology + region constraint,
  ATAC-seq track in FILER                  one row per overlapping track with
  → denominator                            num_overlaps + full metadata
                                           → numerator + ranking
        │                                      │
        └──────────────────┬───────────────────┘
                           ▼
                rate = numerator / denominator
                top-N tracks in the gene, ranked by num_overlaps
                           │
                           ▼ (optional)
              Recipe 2 (filer_find_overlaps.py)
              Pull the actual peak intervals from the
              top-ranked tracks for closer inspection
```

## The commands

### Step 1 — Resolve the gene to a region (Recipe 5)

```bash
python src/scripts/python/gene_to_positions.py \
  --gene CACNA1C \
  --genome-build GRCh38 \
  --out output/05-gene-to-positions/cacna1c.tsv

# Recipe 5 emits the FILER-ready region string in column 6.
REGION=$(awk 'NR==2 {print $6}' output/05-gene-to-positions/cacna1c.tsv)
echo "$REGION"   # e.g. chr12:1970771-2697950
```

The gene body is the region of interest — no expansion window needed. CACNA1C
is unusually large (~660 kb), so expect Recipe 3 to return more tracks here
than a point-variant query would.

### Step 2 — Universe of brain ATAC-seq tracks (Recipe 1, denominator)

```bash
python src/scripts/python/filer_search_tracks.py \
  --genome-build hg38 \
  --assay "ATAC-seq" \
  --tissue-category "Brain" \
  --out output/01-track-discovery/brain_atac.tsv

# subtract the header to count tracks
DENOM=$(($(wc -l < output/01-track-discovery/brain_atac.tsv) - 1))
echo "Brain ATAC-seq tracks in FILER: $DENOM"
```

### Step 3 — Same biology, restricted to the gene, ranked by num_overlaps (Recipe 3, numerator)

```bash
python src/scripts/python/filer_coordinate_search.py \
  --genome-build hg38 \
  --region "$REGION" \
  --filter-string '.assay == "ATAC-seq" and .tissue_category == "Brain"' \
  --full-metadata \
  --count-only 0 \
  --out output/03-coordinate-search/cacna1c_brain_atac.tsv
```

The same `--filter-string` is reused verbatim — that's what guarantees the
numerator is a strict subset of the denominator. If the two filters drift,
the rate is meaningless.

### Step 4 — Compute the rate and the top tracks

Recipe 3 returns one row per overlapping track (no dedup needed). The relevant
columns are `identifier` (col 1), `num_overlaps` (col 38, the last column),
`tissue_category` (col 12), `cell_type` (col 9), and `track_name` (col 37).

```bash
NUMER=$(($(wc -l < output/03-coordinate-search/cacna1c_brain_atac.tsv) - 1))

echo "Tracks with a peak inside CACNA1C : $NUMER / $DENOM"
echo "Rate                              : $(awk -v n="$NUMER" -v d="$DENOM" 'BEGIN{printf "%.1f%%\n", 100*n/d}')"
```

```bash
# Top 20 tracks at the locus, ranked by num_overlaps (descending)
{ head -1 output/03-coordinate-search/cacna1c_brain_atac.tsv ;
  tail -n +2 output/03-coordinate-search/cacna1c_brain_atac.tsv \
    | sort -t$'\t' -k38,38nr ; } \
  | head -21 \
  | cut -f1,38,12,9,37   # identifier, num_overlaps, tissue_category, cell_type, track_name
```

### Step 5 (optional) — Browse the actual peaks in the top tracks (Recipe 2)

The rate and ranking from Step 4 tell you *which* tracks are most active in
CACNA1C, but not *where* inside the ~660 kb gene body the peaks fall. Recipe 2
pulls the underlying interval rows so you can inspect peak coordinates,
signal values, and summit offsets.

Order is preserved from the input file, so reusing the sorted Recipe 3 output
keeps the highest-`num_overlaps` tracks first. `--top 50` caps the run at the
50 most relevant tracks (~12 s at the 1 min / 250 tracks rule of thumb).

```bash
# Save the sorted Recipe 3 output so Recipe 2 sees tracks in rank order.
{ head -1 output/03-coordinate-search/cacna1c_brain_atac.tsv ;
  tail -n +2 output/03-coordinate-search/cacna1c_brain_atac.tsv \
    | sort -t$'\t' -k38,38nr ; } \
  > output/03-coordinate-search/cacna1c_brain_atac_ranked.tsv

python src/scripts/python/filer_find_overlaps.py \
  --region "$REGION" \
  --file output/03-coordinate-search/cacna1c_brain_atac_ranked.tsv \
  --id-col identifier \
  --top 50 \
  --out output/02-track-overlaps/cacna1c_brain_atac_peaks.tsv \
  --json
```

The output is one row per peak per track, with `hitString` carrying the
per-feature fields (peak coordinates, `signalValue`, `pValue`, summit offset,
etc.). Skip this step if you only need the rate and the track ranking — it's
purely for users who want to drill into the peaks themselves.

## How to read the result

A high rate (e.g. 50%+) at the gene level says CACNA1C is **constitutively
accessible in brain** — regulatory elements inside it are real and
broadly used. A moderate rate (10–30%) concentrated in specific cell types
(look at `cell_type` in the top-N output) suggests **cell-type-restricted**
regulation, which is often what GWAS fine-mapping cares about most. A
near-zero rate is a negative result: the locus probably acts through a
different tissue or mechanism. If you ran the optional Step 5, the peak
coordinates show *where inside the gene* the activity concentrates — useful
for cross-referencing GWAS fine-mapped credible sets.