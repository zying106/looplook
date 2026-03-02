# Run Dual Motif Analysis for Loop Anchors

Extracts genomic coordinates for proximal (gene-centric) and distal
(enhancer-centric) anchors of loops connecting to target genes. Scans
sequences against JASPAR core motifs to identify enriched TFBS.

## Usage

``` r
run_distal_motif_analysis(
  target_genes,
  loop_df,
  genome_id,
  pval_thresh,
  current_proj_name,
  out_dir,
  top_n = 5
)
```

## Arguments

- target_genes:

  Character vector of target gene symbols.

- loop_df:

  Data frame. Loop annotation table containing coordinates
  (chr1/start1/end1) and gene assignments.

- genome_id:

  Character. Genome assembly (e.g., "hg19", "mm10") for sequence
  extraction.

- pval_thresh:

  Numeric. P-value cutoff for motifmatchr scanning.

- current_proj_name:

  Character. Project prefix.

- out_dir:

  Character. Output directory.

- top_n:

  Integer. Number of top enriched motifs to output as SeqLogos.

## Value

Invisible `NULL`. Triggers internal plotting functions.
