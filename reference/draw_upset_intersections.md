# Generate UpSet Plot for Gene Set Intersections

Visualizes intersections among multiple gene sets using the classic
UpSetR package. Uses grid graphics capture to ensure plot generation in
all environments.

## Usage

``` r
draw_upset_intersections(
  gene_lists,
  project_name = "UpSet Plot",
  group_colors = NULL
)
```

## Arguments

- gene_lists:

  A named list of character vectors of gene identifiers.

- project_name:

  Character. Used for the file title.

- group_colors:

  Reserved for future use. Currently ignored.

## Value

Invisibly returns the `grob` object.

## Examples

``` r
gene_sets <- list(
  Upregulated = c("TP53", "BRCA1", "MYC", "EGFR"),
  Downregulated = c("BRCA1", "MYC", "CDKN1A", "BAX"),
  Bound_by_TF = c("MYC", "EGFR", "CDKN1A", "KRAS")
)
draw_upset_intersections(
  gene_lists = gene_sets,
  project_name = "Transcriptional Regulation"
)
```
