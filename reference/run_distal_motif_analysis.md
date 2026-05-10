# Run Dual Motif Analysis for Loop Anchors

Run Dual Motif Analysis for Loop Anchors

## Usage

``` r
run_distal_motif_analysis(
  target_genes,
  loop_df,
  genome_id,
  pval_thresh,
  current_proj_name,
  top_n = 5,
  jaspar_db = JASPAR2020::JASPAR2020,
  jaspar_collection = "CORE"
)
```

## Arguments

- jaspar_db:

  A JASPAR database object (e.g., `JASPAR2020::JASPAR2020` or
  `JASPAR2024::JASPAR2024`). Default: `JASPAR2020::JASPAR2020`.

- jaspar_collection:

  Character. JASPAR collection to query (e.g., `"CORE"`, `"CNE"`).
  Default: `"CORE"`.

## Value

A named list containing motif enrichment results and plot objects.
