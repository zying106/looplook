# Extract Target Gene Sets from Annotation Results

Parses loop and target annotations using strict column matching to
extract valid gene lists for downstream profiling.

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

## Arguments

- annotation_res:

  List output from `annotate_peaks_and_loops`.

- src:

  Character. Source to extract from ("loops" or "targets").

- active_loop_types:

  Character vector. Which loop types to extract (if src is "loops").

- include_Filled:

  Logical. Whether to use `_Filled` columns containing merged
  host/target genes.

- use_nearest_gene:

  Logical. Whether to strictly use the 1D nearest gene (SYMBOL) instead
  of 3D loops.

- target_mapping_mode:

  Character. Mapping mode for targets ("all" or "promoter").

## Value

Named list of gene character vectors.
