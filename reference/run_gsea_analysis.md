# Run Custom Gene Set Enrichment Analysis (GSEA)

Evaluates whether target genes are significantly enriched at the
extremes (up/down-regulated) of the global ranked expression list.

## Usage

``` r
run_gsea_analysis(
  target_genes,
  global_glist,
  gsea_ntop,
  current_proj_name,
  out_dir
)
```

## Arguments

- target_genes:

  Character vector. Target gene symbols forming the custom gene set.

- global_glist:

  Named numeric vector. Globally ranked LFC values.

- gsea_ntop:

  Integer. Maximum number of target genes to use (randomly downsampled
  if exceeded).

- current_proj_name:

  Character. Prefix for plot and table.

- out_dir:

  Character. Output directory.

## Value

GSEA result object. Saves a composite GSEA plot (ES curve, barcode, LFC
heatmap) to PDF.
