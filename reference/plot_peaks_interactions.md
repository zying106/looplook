# Publication-ready visualization toolkit for rendering genomic track data and statistical summaries related to 3D chromatin interactions

Generates an integrative genomic track plot resembling IGV, displaying
chromatin loops as arcs, loop anchors as rectangles, optional
overlapping features (e.g., ChIP-seq peaks), and annotated genes. Loop
arcs can be colored or sized by interaction score (7th column in BEDPE).

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

  Character. Chromosome name (e.g., "chr8"). If NULL, inferred from most
  frequent chromosome in BEDPE.

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

  Numeric. Height of anchor rectangles (default: 0.1).

- loop_color:

  Character. Default color for arcs when no score is provided (e.g.,
  "#5D6D7E").

- anchor_color:

  Character. Color for loop anchor rectangles (default: "#3498DB").

- overlap_color:

  Character. Color for overlap track (default: "#E74C3C").

- exon_color:

  Character. Gene exon fill color (default: "#2C3E50").

- intron_color:

  Character. Gene intron line color (default: "black").

- score_to_alpha:

  Logical. Whether to map interaction scores to arc transparency.

- min_score:

  Logical. If TRUE, use score to control arc line width instead of color
  (not yet implemented in current version; future extension).

- save_file:

  Optional character. File path to save the plot (e.g.,
  "region_plot.pdf").

## Value

A `ggplot` object.

## Examples

``` r
# 1. Get paths to example files included in the package
bedpe_path <- system.file("extdata", "example_loops_1.bedpe", package = "looplook")
bed_path <- system.file("extdata", "example_peaks.bed", package = "looplook")

# 2. Run plotting (requires TxDb package for gene annotation)
if (bedpe_path != "" &&
  requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE) &&
  requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  # Example A: Basic plot with loops only
  # Plotting the region around the example loop (chr1:10,000-50,000)
  p1 <- plot_peaks_interactions(
    bedpe_file = bedpe_path,
    target_bed = bed_path,
    chr = "chr1",
    from = 11884299,
    to = 12106581,
    species = "hg38"
  )
  print(p1)

  # Example B: Integrative plot with overlapping peaks and output to file
  p2 <- plot_peaks_interactions(
    bedpe_file = bedpe_path,
    target_bed = bed_path,
    chr = "chr1",
    from = 11884299,
    to = 12106581,
    species = "hg38",
    save_file = tempfile(fileext = ".pdf")
  )
}
#>   2169 genes were dropped because they have exons located on both strands of
#>   the same reference sequence or on more than one reference sequence, so cannot
#>   be represented by a single genomic range.
#>   Use 'single.strand.genes.only=FALSE' to get all the genes in a GRangesList
#>   object, or use suppressMessages() to suppress this message.
#> 'select()' returned 1:1 mapping between keys and columns
#> 'select()' returned 1:1 mapping between keys and columns
#> Warning: The `size` argument of `element_rect()` is deprecated as of ggplot2 3.4.0.
#> ℹ Please use the `linewidth` argument instead.
#> ℹ The deprecated feature was likely used in the looplook package.
#>   Please report the issue at <https://github.com/zying106/looplook/issues>.

#>   2169 genes were dropped because they have exons located on both strands of
#>   the same reference sequence or on more than one reference sequence, so cannot
#>   be represented by a single genomic range.
#>   Use 'single.strand.genes.only=FALSE' to get all the genes in a GRangesList
#>   object, or use suppressMessages() to suppress this message.
#> 'select()' returned 1:1 mapping between keys and columns
#> 'select()' returned 1:1 mapping between keys and columns
```
