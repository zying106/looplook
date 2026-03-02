# Draw Simplified Flower Plot for Core vs. Unique Genes

Creates a circular "flower" diagram where each petal represents a gene
set, showing the number of genes unique to that set. The center displays
the size of the core intersection (genes shared by all sets). Designed
for intuitive comparison of shared vs. condition-specific genes across
2–6 groups.

## Usage

``` r
draw_flower_simplified(gene_lists, project_name, filename, group_colors)
```

## Arguments

- gene_lists:

  A named list of character vectors containing gene identifiers (e.g.,
  symbols or Entrez IDs).

- project_name:

  Character. Prefix for the plot title (e.g., "Differential Loops").

- filename:

  Character. Output file path (e.g., "flower.png" or "flower.pdf").

- group_colors:

  Character vector. Colors for each group (named or in same order as
  `gene_lists`).

## Value

Invisibly returns `NULL`. Saves the plot to `filename`.

## Examples

``` r
# 1. Create dummy gene sets with some overlap
gene_sets <- list(
  Control = c("TP53", "BRCA1", "MYC"),
  Treated = c("BRCA1", "MYC", "EGFR"),
  Resistant = c("MYC", "EGFR", "KRAS")
)
# 2. Draw the flower plot
draw_flower_simplified(
  gene_lists = gene_sets,
  project_name = "Drug Response Study",
  filename = tempfile(fileext = ".png"),
  group_colors = c(Control = "#E41A1C", Treated = "#377EB8", Resistant = "#4DAF4A")
)
#>     Saved (Simplified Flower Plot with inner counts): /tmp/RtmpkfVz6j/file2312463bce2b.png
```
