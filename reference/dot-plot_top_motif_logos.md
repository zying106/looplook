# Plot Top Motif Sequence Logos

Retrieves Position Frequency Matrices (PFM) from JASPAR and plots
sequence logos for the top N enriched motifs.

## Usage

``` r
.plot_top_motif_logos(res_df, out_dir, prefix, top_n)
```

## Arguments

- res_df:

  Data frame of motif enrichment results.

- out_dir:

  Character. Output directory.

- prefix:

  Character. File prefix.

- top_n:

  Integer. Number of top motifs to plot.

## Value

Invisible `NULL`.
