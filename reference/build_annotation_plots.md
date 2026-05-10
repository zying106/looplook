# Internal: Build Annotation Visualization Suite

Generates the complete diagnostic plot collection (donut, rose,
karyotype heatmaps, flower plot, pie charts) for the annotation results.

## Usage

``` r
build_annotation_plots(
  plot_df,
  bed_info,
  cluster_info,
  target_connected_loops,
  txdb_obj,
  org_db_pkg,
  species,
  project_name,
  color_palette,
  karyo_bin_size
)
```

## Arguments

- plot_df:

  Loop annotation data frame with loop_type and All_Anchor_Genes.

- bed_info:

  Target annotation data frame (optional).

- cluster_info:

  Anchor annotation data frame.

- target_connected_loops:

  Data frame of target-connected loops (optional).

- txdb_obj:

  TxDb object for gene coordinate lookup.

- org_db_pkg:

  Character. Organism database package name.

- species:

  Character. Genome assembly.

- project_name:

  Character. Project prefix for plot titles.

- color_palette:

  Character. RColorBrewer palette name.

- karyo_bin_size:

  Integer. Bin size for karyotype heatmaps.

## Value

A named list of ggplot / grob objects.
