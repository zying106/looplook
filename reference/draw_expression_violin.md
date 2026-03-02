# Draw Expression Violin Plot

Creates a violin + boxplot showing log2-transformed gene expression
grouped by loop type (e.g., promoter-enhancer, enhancer-enhancer).

## Usage

``` r
draw_expression_violin(
  plot_data,
  project_name,
  filename,
  unit_type,
  group_colors
)
```

## Arguments

- plot_data:

  (data.frame) Must contain columns: `loop_type`, `expression_value`.

- project_name:

  (character) Project or sample name for plot title.

- filename:

  (character) Output file path (e.g., "expr_violin.pdf").

- unit_type:

  (character) Expression unit (e.g., "TPM", "FPKM").

- group_colors:

  (character vector) Named or ordered colors for each `loop_type`.

## Value

A `ggplot` object of the violin + boxplot.
