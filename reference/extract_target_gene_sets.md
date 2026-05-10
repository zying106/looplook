# Extract Target Gene Sets from Annotation Results

Parses loop and target annotations to extract valid gene lists.

## Usage

``` r
extract_target_gene_sets(
  annotation_res,
  src,
  active_loop_types = NULL,
  include_Filled = TRUE,
  use_nearest_gene = FALSE,
  target_mapping_mode = "all"
)
```

## Value

A named list of character vectors, each containing target gene symbols.
