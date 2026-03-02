# Construct and Visualize STRING PPI Network

Maps target genes to STRING database, extracts the high-confidence
protein-protein interaction subnetwork, removes isolated nodes, and
visualizes the network colored by LFC and sized by degree.

## Usage

``` r
run_ppi_analysis(
  target_genes,
  global_glist,
  org_db,
  ppi_score,
  ppi_ntop,
  current_proj_name,
  out_dir
)
```

## Arguments

- target_genes:

  Character vector of gene symbols.

- global_glist:

  Named numeric vector of LFC values.

- org_db:

  Character. Organism database to automatically infer STRING species ID.

- ppi_score:

  Numeric. Minimum combined score threshold for STRING edges (e.g.,
  400).

- ppi_ntop:

  Integer. Maximum number of genes to include (prioritized by highest
  absolute LFC).

- current_proj_name:

  Character. Project prefix.

- out_dir:

  Character. Output directory.

## Value

Invisible `NULL`. Saves a `ggraph` network PDF.
