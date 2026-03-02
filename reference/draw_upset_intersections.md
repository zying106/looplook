# Generate UpSet Plot for Gene Set Intersections

Visualizes intersections among multiple gene sets using the classic
UpSetR package. Uses grid graphics capture to ensure plot generation in
all environments.

## Usage

``` r
draw_upset_intersections(
  gene_lists,
  project_name,
  filename,
  group_colors = NULL
)
```

## Arguments

- gene_lists:

  A named list of character vectors of gene identifiers.

- project_name:

  Character. Used for the file title.

- filename:

  Character. Output file path (must end in .pdf or .png).

- group_colors:

  Optional named character vector. Not used in this version.

## Value

Invisibly returns `NULL`. Saves the plot to `filename`.

## Examples

``` r
gene_sets <- list(
  Upregulated = c("A", "B", "C", "D"),
  Downregulated = c("C", "D", "E", "F"),
  Bound_by_TF = c("B", "D", "F", "G")
)

tf <- tempfile(fileext = ".pdf")

if (requireNamespace("UpSetR", quietly = TRUE)) {
  draw_upset_intersections(
    gene_lists = gene_sets,
    project_name = "Transcriptional Regulation",
    filename = tf
  )
}
#> Warning: `aes_string()` was deprecated in ggplot2 3.0.0.
#> ℹ Please use tidy evaluation idioms with `aes()`.
#> ℹ See also `vignette("ggplot2-in-packages")` for more information.
#> ℹ The deprecated feature was likely used in the UpSetR package.
#>   Please report the issue to the authors.
#> Warning: The `size` argument of `element_line()` is deprecated as of ggplot2 3.4.0.
#> ℹ Please use the `linewidth` argument instead.
#> ℹ The deprecated feature was likely used in the UpSetR package.
#>   Please report the issue to the authors.
```
