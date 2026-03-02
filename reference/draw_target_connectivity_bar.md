# Internal: Draw Target-Loop Connectivity

Visualizes how many target peaks overlap with loops, and breakdown by
Loop Type.

## Usage

``` r
draw_target_connectivity_bar(
  bed_info,
  cluster_info,
  project_name,
  filename,
  color_palette = "Set2"
)
```

## Arguments

- bed_info:

  Data frame containing peak information.

- cluster_info:

  Data frame containing cluster information (currently unused in plot
  but kept for consistency).

- project_name:

  Character string for the project title.

- filename:

  Character string for the output filename.

- color_palette:

  Character string for the color palette (default: "Set2").

## Value

A ggplot object representing the connectivity bar chart.
