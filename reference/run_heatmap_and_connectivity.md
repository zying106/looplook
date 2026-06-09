# Generate Expression Heatmap and Connectivity Plots

Generate Expression Heatmap and Connectivity Plots

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
  source_type,
  target_col = NULL,
  skip_heatmap = FALSE
)
```

## Value

A named list of plot objects (Heatmap, Scatter, Raincloud_LFC,
Raincloud_Expr).
