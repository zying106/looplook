# Internal: Draw Target Genomic Annotation Distribution (Pie Chart)

Optimized: Labels only show Count (%), and hides labels for small slices
(\<2%) to avoid overlap.

## Usage

``` r
draw_target_annotation_pie(
  bed_info,
  project_name,
  filename,
  color_palette = "Set3"
)
```

## Arguments

- bed_info:

  Data frame containing annotation information.

- project_name:

  Character string for the project title.

- filename:

  Character string for the output filename.

- color_palette:

  Character string for the color palette (default: "Set3").

## Value

A ggplot object representing the annotation pie chart.
