# Internal: Draw Circular Bar Plot (Gene Counts)

Optimized for strictly vertical label alignment using y=0 anchor.

## Usage

``` r
draw_circular_bar_plot(data_df, project_name, filename, color_vec)
```

## Arguments

- data_df:

  Data frame containing loop and gene information.

- project_name:

  Character string for the project title.

- filename:

  Character string for the output filename.

- color_vec:

  Named character vector for colors.

## Value

A ggplot object representing the circular bar plot.
