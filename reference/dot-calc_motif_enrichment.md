# Calculate Motif Enrichment via Fisher's Exact Test

Calculate Motif Enrichment via Fisher's Exact Test

## Usage

``` r
.calc_motif_enrichment(
  fg_gr,
  bg_gr,
  genome_obj,
  pval_thresh,
  species_id,
  jaspar_db = JASPAR2020::JASPAR2020,
  jaspar_collection = "CORE"
)
```

## Value

A data frame of enrichment results, or `NULL` if input is empty.
