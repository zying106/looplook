# Construct and Visualize STRING PPI Network

Construct and Visualize STRING PPI Network

## Usage

``` r
run_ppi_analysis(
  target_genes,
  global_glist,
  org_db,
  ppi_score,
  ppi_ntop,
  current_proj_name
)
```

## Value

A `ggplot` object representing the PPI network, or `NULL` if no
interactions found.
