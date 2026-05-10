# Internal: Draw Karyotype Heatmap

Creates a genome-wide heatmap of genomic feature density (e.g., loops)
across chromosomes, binned by a fixed window size, and rendered as a
karyotype plot.

## Usage

``` r
draw_karyo_heatmap_internal(
  gr_data,
  title_prefix,
  bin_size,
  sat_level,
  ref_txdb,
  plot_species,
  unit_label,
  custom_colors = NULL
)
```

## Arguments

- gr_data:

  (GRanges) Genomic ranges to visualize (e.g., loop anchors).

- title_prefix:

  (character) Subtitle descriptor (e.g., sample name).

- bin_size:

  (integer) Bin width in base pairs (e.g., 1e6 for 1 Mb).

- sat_level:

  (numeric) Quantile (0–1) for color saturation (e.g., 0.95).

- ref_txdb:

  (TxDb or similar) Reference genome annotation for chromosome lengths.

- plot_species:

  (character) Genome build/species code (e.g., "hg38", "mm10").

- unit_label:

  (character) Unit for load annotation (e.g., "loops").

- custom_colors:

  (character vector) Optional custom color palette.

## Value

A `looplook_karyo` object wrapping a rendered PNG payload. Use
[`print()`](https://rdrr.io/r/base/print.html) to display.
