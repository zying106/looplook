# Generate LFC Violin and Boxplot

Visualizes the Log2 Fold Change distribution of target genes compared to
the genomic background.

## Usage

``` r
run_lfc_violin(
  target_genes,
  global_glist,
  stat_test = "wilcox.test",
  project_name,
  out_dir
)
```

## Arguments

- target_genes:

  Character vector of target gene symbols.

- global_glist:

  Named numeric vector of global log fold changes (names = gene
  symbols).

- stat_test:

  Character. Statistical test ("t.test" or "wilcox.test").

- project_name:

  Character. Project-specific prefix for output.

- out_dir:

  Character. Output directory.

## Value

Invisible `NULL`. Saves a narrow-format PDF plot.
