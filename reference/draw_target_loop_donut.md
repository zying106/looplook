# Internal: Draw Target-Associated Loop Distribution (Donut Chart)

Visualizes the distribution of loop types that are connected to target
peaks.

## Usage

``` r
draw_target_loop_donut(loop_data, project_name, filename, color_vec)
```

## Arguments

- loop_data:

  Data frame containing loop information.

- project_name:

  Character string for the project title.

- filename:

  Character string for the output filename.

- color_vec:

  Named character vector for colors.

## Value

A ggplot object representing the donut chart of loop type distribution.
