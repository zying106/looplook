# Internal: Draw Comparison Bar Chart

Visualizes the comparison of loop counts between original and filtered
datasets.

## Usage

``` r
draw_comparison_bar(original_df, filtered_df, filename, color_vec)
```

## Arguments

- original_df:

  Data frame containing the original loop data.

- filtered_df:

  Data frame containing the filtered loop data.

- filename:

  Character string for the output filename.

- color_vec:

  Named character vector for colors.

## Value

A ggplot object visualizing the comparison between groups.
