# Calculate Motif Enrichment via Fisher's Exact Test

Compares motif hit frequencies between foreground (anchor regions) and
background (random loops) sequences.

## Usage

``` r
.calc_motif_enrichment(fg_gr, bg_gr, genome_obj, pval_thresh, species_id)
```

## Arguments

- fg_gr:

  GRanges object. Foreground genomic regions.

- bg_gr:

  GRanges object. Background genomic regions.

- genome_obj:

  BSgenome object. Reference genome sequence.

- pval_thresh:

  Numeric. Cutoff for calling a motif match.

- species_id:

  Numeric. Taxonomy ID for JASPAR filtering (e.g., 9606).

## Value

Data frame of motif enrichment statistics (Pvalue, FDR, OddsRatio).
