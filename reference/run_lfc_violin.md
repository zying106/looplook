# Generate LFC Violin and Boxplot

Generate LFC Violin and Boxplot

## Usage

``` r
run_lfc_violin(
  target_genes,
  global_glist,
  stat_test = "wilcox.test",
  project_name
)
```

## Value

A `ggplot` object, or `NULL` if fewer than 3 valid targets.
