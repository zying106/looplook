# Perform GO Enrichment and Generate Network Plot

Perform GO Enrichment and Generate Network Plot

## Usage

``` r
run_go_enrichment(
  genes,
  org_db,
  universe_genes,
  cnet_nSample = 50,
  project_name = "Analysis"
)
```

## Value

A list with `result` (data frame) and `plot` (ggplot) elements.
