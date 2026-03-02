# Generate Summary GO Lollipop Facet Plot

Combines GO enrichment results from multiple loop types (e.g., E-P, P-P)
and plots them side-by-side in a faceted lollipop chart, categorized by
GO Ontology (BP, CC, MF).

## Usage

``` r
plot_summary_go_lollipop(all_go_results, base_project_name, out_dir)
```

## Arguments

- all_go_results:

  List of data frames containing GO results.

- base_project_name:

  Character. Prefix for output.

- out_dir:

  Character. Output directory.

## Value

Invisible `NULL`. Saves a combined PDF plot.
