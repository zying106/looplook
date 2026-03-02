# Generate Expression Heatmap and Connectivity Raincloud Plots

Plots Z-score normalized expression heatmaps for target genes. Also
generates scatter and Raincloud plots correlating 3D connectivity
(Degree) with expression levels and LFC, inheriting Hub classifications
from upstream.

## Usage

``` r
run_heatmap_and_connectivity(
  target_genes,
  tpm_mat_raw,
  meta_raw,
  loop_stats_df,
  global_glist,
  heatmap_ntop,
  cor_method,
  current_proj_name,
  out_dir,
  source_type,
  target_col = NULL,
  skip_heatmap = FALSE
)
```

## Arguments

- target_genes:

  Character vector of gene symbols.

- tpm_mat_raw:

  Data frame. Raw TPM/FPKM expression matrix.

- meta_raw:

  Data frame. Sample metadata.

- loop_stats_df:

  Data frame. Promoter or Distal statistics containing node degree and
  Hub labels.

- global_glist:

  Named numeric vector. Global LFC values.

- heatmap_ntop:

  Integer. Max highly variable genes for the heatmap.

- cor_method:

  Character. Correlation method.

- current_proj_name:

  Character. Project prefix.

- out_dir:

  Character. Output directory.

- source_type:

  Character. Source type ("loops" or "targets") for plot subtitling.

- target_col:

  Character or NULL. Specific column in `loop_stats_df` to use as
  connectivity degree (e.g., "n_Linked_Distal").

- skip_heatmap:

  Logical. If `TRUE`, skips drawing the heatmap and only draws
  connectivity plots.

## Value

Invisible `NULL`. Saves multiple PDF plots.
