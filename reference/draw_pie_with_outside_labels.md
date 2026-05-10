# Internal: Draw Pie Chart with Outside Labels

Simplified pie chart with labels placed outside the slices, using
RColorBrewer palettes for genomic annotation categories.

## Usage

``` r
draw_pie_with_outside_labels(data_df, group_col, title, palette)
```

## Arguments

- data_df:

  Data frame with an annotation column.

- group_col:

  Character. Column name for grouping.

- title:

  Character. Plot title.

- palette:

  Character. RColorBrewer palette name.

## Value

A ggplot object, or NULL if data is empty.
