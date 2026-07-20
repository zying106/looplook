# Filter Chromatin Loops by Blacklist and/or Region of Interest

Applies optional blacklist and region-of-interest (ROI) filtering to a
single
[`GInteractions`](https://rdrr.io/pkg/InteractionSet/man/GInteractions-class.html)
object. This is the single-sample counterpart of
[`consolidate_chromatin_loops`](https://zying106.github.io/looplook/reference/consolidate_chromatin_loops.md)
without the merge/consensus step.

## Usage

``` r
filter_chromatin_loops(
  gi,
  blacklist_species = NULL,
  region_of_interest = NULL,
  roi_mode = c("any", "both"),
  min_score = NULL,
  quiet = FALSE
)
```

## Arguments

- gi:

  A
  [`GInteractions`](https://rdrr.io/pkg/InteractionSet/man/GInteractions-class.html)
  object.

- blacklist_species:

  Character. Species/build for built-in ENCODE blacklist (`"hg38"`,
  `"hg19"`, `"mm10"`, `"mm9"`), a path to a custom BED file, or `NULL`
  to skip blacklist filtering.

- region_of_interest:

  Character. Path to a BED file defining regions of interest, or `NULL`
  to skip ROI filtering.

- roi_mode:

  Character. `"any"` (default): keep loops where *either* anchor
  overlaps the ROI. `"both"`: keep loops where *both* anchors overlap
  the ROI.

- min_score:

  Numeric or `NULL`. Minimum interaction score to retain a loop. Loops
  with `score < min_score` are discarded. Assumes higher scores = better
  interactions. Default: `NULL` (no score filtering).

- quiet:

  Logical. If `TRUE`, suppress progress messages. Default: `FALSE`.

## Value

A filtered
[`GInteractions`](https://rdrr.io/pkg/InteractionSet/man/GInteractions-class.html)
object.

## Examples

``` r
gi <- bedpe_to_gi(system.file("extdata", "example_loops_1.bedpe", package = "looplook"))

# Blacklist only
gi_clean <- filter_chromatin_loops(gi, blacklist_species = "hg38", quiet = TRUE)

# ROI (any anchor) with custom BED
roi_path <- system.file("extdata", "example_k27ac_peaks.bed", package = "looplook")
gi_roi <- filter_chromatin_loops(gi, region_of_interest = roi_path, roi_mode = "any", quiet = TRUE)

# Score + blacklist + ROI (both)
gi_strict <- filter_chromatin_loops(gi,
  min_score = 5,
  blacklist_species = "hg38", region_of_interest = roi_path,
  roi_mode = "both", quiet = TRUE
)
```
