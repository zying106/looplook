# Internal: Build Refinement Visualization Suite

Generates diagnostic plots for expression-aware refinement: dumbbell
comparison, donut, Sankey, karyotype heatmaps, and rose plot.

## Usage

``` r
build_refinement_plots(
  original_loop_df,
  loop_df,
  bed_info,
  whitelist,
  project_name,
  karyo_bin_size,
  species
)
```

## Arguments

- original_loop_df:

  Loop annotation before refinement.

- loop_df:

  Loop annotation after refinement.

- bed_info:

  Target annotation data frame (optional).

- whitelist:

  Character vector of active gene symbols.

- project_name:

  Character. Project prefix for plot titles.

- karyo_bin_size:

  Integer. Bin size for karyotype heatmaps.

- species:

  Character. Genome assembly.

## Value

A named list of ggplot / htmlwidget / grob objects.
