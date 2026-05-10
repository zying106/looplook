# Integrative visualization of 3D chromatin loops and genomic features

Generates an integrative genomic track plot displaying chromatin loops
as arcs, loop anchors as rectangles, optional overlapping features
(e.g., ChIP-seq peaks), and annotated genes. Loop arcs can be colored or
sized by interaction score (7th column in BEDPE).

## Usage

``` r
plot_peaks_interactions(
  bedpe_file,
  target_bed = NULL,
  chr = NULL,
  from = NULL,
  to = NULL,
  species = "hg38",
  max_levels = 10,
  base_anchor_height = 0.05,
  loop_color = "#5D6D7E",
  anchor_color = "#3498DB",
  overlap_color = "#02ABB4",
  exon_color = "#2C3E50",
  intron_color = "black",
  score_to_alpha = TRUE,
  min_score = NULL,
  save_file = NULL
)
```

## Arguments

- bedpe_file:

  Character. Path to a BEDPE file (at least 6 columns; 7th column used
  as score if present).

- target_bed:

  Optional character. Path to a BED file (e.g., peaks) to overlay below
  the loop track.

- chr:

  Character. Chromosome name (e.g., "chr8"). If NULL, inferred from the
  most frequent chromosome in the BEDPE.

- from:

  Numeric. Start coordinate of the region to plot.

- to:

  Numeric. End coordinate of the region to plot.

- species:

  Character. Genome assembly: "hg38", "hg19", "mm10", or "mm9".

- max_levels:

  Integer. Maximum number of vertical levels for arc stacking (default:
  10).

- base_anchor_height:

  Numeric. Height of anchor rectangles (default: 0.05).

- loop_color:

  Character. Default color for arcs when no score is provided (default:
  "#5D6D7E").

- anchor_color:

  Character. Color for loop anchor rectangles (default: "#3498DB").

- overlap_color:

  Character. Color for overlap track (default: "#02ABB4").

- exon_color:

  Character. Gene exon fill color (default: "#2C3E50").

- intron_color:

  Character. Gene intron line color (default: "black").

- score_to_alpha:

  Logical. Whether to map interaction scores to arc transparency.

- min_score:

  Optional numeric. Floor value for score mapping.

- save_file:

  Character. Optional path to save the plot via
  [`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html).
  When set, the plot is written to this file and returned invisibly.

## Value

A `ggplot` object (invisibly when `save_file` is non-NULL).

## Examples

``` r
if (requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE) &&
  requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  bedpe_path <- system.file("extdata", "example_loops_1.bedpe", package = "looplook")
  bed_path <- system.file("extdata", "example_peaks.bed", package = "looplook")
  p <- plot_peaks_interactions(
    bedpe_file = bedpe_path,
    target_bed = bed_path,
    chr = "chr1",
    from = 11884299,
    to   = 12106581,
    species = "hg38"
  )
  print(p)
}
```
