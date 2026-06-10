# Draw Simplified Flower Plot for Core vs. Unique Genes

Creates a circular "flower" diagram where each petal represents a gene
set, showing the number of genes unique to that set. The center displays
the size of the core intersection (genes shared by all sets). Designed
for intuitive comparison of shared vs. condition-specific genes across
2-6 groups.

## Usage

``` r
draw_flower_simplified(gene_lists, project_name, group_colors)
```

## Arguments

- gene_lists:

  A named list of character vectors containing gene identifiers.

- project_name:

  Character. Prefix for the plot title.

- group_colors:

  Named character vector for specific group mappings.

## Value

Invisibly returns the `ggplot` object.

## Examples

``` r
gene_sets <- list(
    Control = c("TP53", "BRCA1", "MYC", "EGFR"),
    Treated = c("BRCA1", "MYC", "EGFR", "KRAS"),
    Resistant = c("MYC", "EGFR", "KRAS", "BRAF")
)
draw_flower_simplified(
    gene_lists = gene_sets,
    project_name = "Drug Response",
    group_colors = c(Control = "#E41A1C", Treated = "#377EB8", Resistant = "#4DAF4A")
)
```
