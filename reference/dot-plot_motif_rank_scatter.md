# Plot Motif Rank Scatter

Creates a scatter plot ranking motifs by FDR, where point size
represents Odds Ratio and color represents the TF Family.

## Usage

``` r
.plot_motif_rank_scatter(res_df, out_dir, prefix, fdr_thresh = 0.05)
```

## Arguments

- res_df:

  Data frame of motif enrichment results.

- out_dir:

  Character. Output directory.

- prefix:

  Character. File prefix.

- fdr_thresh:

  Numeric. FDR threshold for significance coloration.

## Value

Invisible `NULL`.
