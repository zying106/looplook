# Draw Enhancer Source Distribution

Bar plot showing the origin of functional enhancer anchors (e.g., native
vs. promoter-derived).

## Usage

``` r
draw_enhancer_source_distribution(loop_data, project_name, filename)
```

## Arguments

- loop_data:

  (data.frame) Must contain:

  - `functional_anchor1_type`, `anchor1_source`

  - `functional_anchor2_type`, `anchor2_source`

- project_name:

  (character) Project name for title. Only anchors where type == "E" are
  considered.

- filename:

  (character) Output file path (e.g., "enh_sources.pdf").

## Value

A `ggplot` object of the bar plot showing enhancer source distribution.
