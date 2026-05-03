# Worked Example — How active is the schizophrenia GWAS variant rs1006737 across adult-brain ATAC-seq?

## The question

> **rs1006737 (CACNA1C intron 3) is the strongest common GWAS hit shared across
> schizophrenia and bipolar disorder. It is noncoding — so any effect must go
> through regulatory chromatin. Of all adult-brain ATAC-seq tracks in FILER,
> what fraction actually call a peak there, and which tracks have the strongest
> signal?**

For a missense variant ("does it disrupt the protein?") a single overlap is
enough to tell a story. For a noncoding GWAS hit, a single overlap means
nothing — you need a **rate**. A locus that is open in 47/300 adult-brain
ATAC-seq tracks is meaningfully active; 5/300 is background.

This recipe pairs **Recipe 1** (universe → denominator) with **Recipe 10**
(filtered, locus-restricted, with peak intervals → numerator + evidence).

## Recipe chain

```
Recipe 4 (rsid_to_positions.py)
  rs1006737 → chr12:2236129-2236130
        │
        ├──────────────────────────────────────┐
        ▼                                      ▼
Recipe 1 (filer_search_tracks.py)        Recipe 10 (filer_filter_then_overlaps.py)
  Full universe: every adult-brain         Same biology + region constraint,
  ATAC-seq track in FILER                  with peak intervals attached
  → denominator                            → numerator + ranked tracks
        │                                      │
        └──────────────────┬───────────────────┘
                           ▼
                rate = numerator / denominator
                top-N ranked tracks at the locus
```

## The commands

### Step 1 — Resolve the variant (Recipe 4)

```bash
python src/scripts/python/rsid_to_positions.py \
  --rsid rs1006737 \
  --genome-build GRCh38 \
  --out output/4-rsid-to-positions/cacna1c.tsv

POINT=$(awk 'NR==2 {print $5}' output/4-rsid-to-positions/cacna1c.tsv)
# Expand to a small window so we can rank tracks by num_overlaps in the locus,
# not just at the single base. 5 kb on each side is a reasonable peak-scale window.
CHR=$(echo "$POINT" | cut -d: -f1)
POS=$(echo "$POINT" | cut -d: -f2 | cut -d- -f1)
REGION="${CHR}:$((POS-5000))-$((POS+5000))"
echo "$REGION"   # e.g. chr12:2231129-2241129
```

### Step 2 — Universe of adult-brain ATAC-seq tracks (Recipe 1, denominator)

```bash
python src/scripts/python/filer_search_tracks.py \
  --genome-build hg38 \
  --assay "ATAC-seq" \
  --tissue-category "Brain" \
  --out output/01-track-discovery/brain_atac.tsv

# subtract the header to count tracks
DENOM=$(($(wc -l < output/01-track-discovery/brain_atac.tsv) - 1))
echo "Adult-brain ATAC-seq tracks in FILER: $DENOM"
```

If you also want to restrict to adult life stage (Recipe 1 doesn't expose
`life_stage` as a named flag, so use `--filter-string` for parity with Step 3):

```bash
python src/scripts/python/filer_search_tracks.py \
  --genome-build hg38 \
  --filter-string '.assay == "ATAC-seq" and .tissue_category == "Brain" and .life_stage == "Adult"' \
  --out output/01-track-discovery/brain_atac_adult.tsv

DENOM=$(($(wc -l < output/01-track-discovery/brain_atac_adult.tsv) - 1))
```

### Step 3 — Same biology, restricted to the locus, with peaks attached (Recipe 10, numerator)

```bash
python src/scripts/python/filer_filter_then_overlaps.py \
  --genome-build hg38 \
  --region "$REGION" \
  --filter-string '.assay == "ATAC-seq" and .tissue_category == "Brain" and .life_stage == "Adult"' \
  --top 200 \
  --out output/10-filter-then-overlaps/cacna1c_brain_atac.tsv
```

The same `--filter-string` is reused verbatim — that's what guarantees the
numerator is a strict subset of the denominator. If the two filters drift,
the rate is meaningless.

### Step 4 — Compute the rate and the top tracks

```bash
# unique tracks calling a peak in the locus (one row per interval, so dedupe on identifier)
NUMER=$(awk -F'\t' 'NR>1 {print $1}' output/10-filter-then-overlaps/cacna1c_brain_atac.tsv \
  | sort -u | wc -l | tr -d ' ')

echo "Tracks with a peak at rs1006737 ± 5 kb : $NUMER / $DENOM"
echo "Rate                                  : $(awk -v n="$NUMER" -v d="$DENOM" 'BEGIN{printf "%.1f%%\n", 100*n/d}')"
```

```bash
# Top 20 tracks at the locus, by num_overlaps (one row per track)
awk -F'\t' 'NR==1 || !seen[$1]++' \
    output/10-filter-then-overlaps/cacna1c_brain_atac.tsv \
  | sort -t$'\t' -k4,4nr \
  | head -21 \
  | cut -f1,4,7,8,11   # identifier, num_overlaps, tissue_category, life_stage, track_name
```

## How to read the result

A high rate (e.g. 50%+) at an adult-brain locus says the variant sits in
chromatin that is **constitutively open in adult brain** — the regulatory
element is real and broad. A moderate rate (10–30%) concentrated in specific
cell types (look at `cell_type` in the top-N output) suggests a **cell-type-
restricted** element, which is often what GWAS fine-mapping cares about most.
A near-zero rate is a negative result: the variant probably acts through a
different tissue, a different mechanism, or a different SNP in LD with it.