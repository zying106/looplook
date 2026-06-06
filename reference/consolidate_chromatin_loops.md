# Consolidate and Integrate Chromatin Loops from Replicates or Multiple Sources

This function consolidates chromatin loops from multiple BEDPE files. It
is designed for two main purposes:

1.  **Replicate Consolidation**: Merging biological or technical
    replicates to identify high-confidence, reproducible loops (e.g., 3
    replicates of H3K27ac HiChIP).

2.  **Multi-Omics Integration**: The framework can be used to identify
    multi-source consensus by integrating datasets from various
    experimental designs, such as HiChIP assays targeting different
    factors (e.g., integrating **H3K27ac** and **H3K4me3**, or
    overlapping Hi-C with ChIA-PET).

The function supports three modes:

- `"consensus"`: Implements graph-based connected component analysis to
  cluster spatially proximal anchors across samples. Only retains
  clusters detected in \>= min_consensus biological replicates.

- `"intersect"`: Enforces strict reference-based filtering, retaining
  loops that show full genomic overlap with the reference file (File 1).

- `"union"`: Retains all chromatin interactions across the entire
  cohort, ideal for exploratory pan-tissue analyses.

**Connected-component chaining**: Graph-based clustering may
transitively chain loci (A–B and B–C merge, pulling A and C into the
same cluster even if they are far apart). Inspect `n_members` in the
output and reduce `gap` if clusters appear inflated.

It also supports a **two-stage filtering strategy** to maximize
signal-to-noise ratio:

- **Pre-filtering** (`min_raw_score`): Removes low-confidence noise
  (e.g., singleton reads) from raw files *before* merging.

- **Post-filtering** (`min_score`): Filters the final consensus loops
  based on their aggregated score.

- **Replicate-Balanced Aggregation**: In `"consensus"` and `"union"`
  modes, each cluster score is computed as the mean of per-source mean
  scores, so one replicate with many fragmented calls cannot dominate
  the final representative score.

## Usage

``` r
consolidate_chromatin_loops(
  files = NULL,
  gap = 1000,
  mode = c("consensus", "intersect", "union"),
  min_consensus = NULL,
  score_col = NULL,
  min_raw_score = NULL,
  min_score = NULL,
  blacklist_species = NULL,
  region_of_interest = NULL,
  roi_mode = c("any", "both"),
  out_file = NULL,
  write_output = TRUE,
  quiet = FALSE
)
```

## Arguments

- files:

  Character vector. Paths to BEDPE files (at least two).

- gap:

  Numeric. Distance (bp) to consider loops as overlapping. Default 1000.

- mode:

  Character. Choose one of the following: "consensus", "intersect",
  "union". Merge strategy:

  - `"intersect"`: Strict reference-based filtering (keeps loops in File
    1 supported by ALL other files).

  - `"union"`: Merges all detected loops into a comprehensive map.

  - `"consensus"`: Graph-based clustering to find a consensus set
    supported by a majority of samples. (Formerly "reproducible").

- min_consensus:

  Integer. Minimum number of replicates a loop must appear in (only
  effective when `mode = "consensus"`). If `NULL` (default), the
  threshold is automatically calculated:

  - For 2–3 replicates: Requires all (100\\

  - For \\N \ge 4\\: Requires \\\lfloor 0.75N \rfloor + 1\\ (e.g., 3 for
    N=4, 4 for N=5, 7 for N=8).

- score_col:

  Integer or `NULL`. Column index to use as interaction score when
  reading BEDPE files. Passed to
  [`bedpe_to_gi`](https://zying106.github.io/looplook/reference/bedpe_to_gi.md).
  If `NULL` (default), auto-detection is used (see
  [`bedpe_to_gi`](https://zying106.github.io/looplook/reference/bedpe_to_gi.md)).

- min_raw_score:

  Numeric. **Pre-filtering threshold**. Loops with a raw score (e.g.,
  read count) below this value in individual files will be discarded
  **before** any merging or intersection.

  - Recommended value: `2` (to remove singleton noise loops with
    count=1).

  - Default: `NULL` (no pre-filtering).

- min_score:

  Numeric. **Post-filtering threshold**. Minimum score to keep a
  consolidated loop **after** merging.

  - For `"consensus"` and `"union"` modes, this filters loops based on a
    replicate-balanced representative score (mean of per-source cluster
    means).

  - For `"intersect"` mode, this filters the retained File 1 loops by
    their original score.

  - Default: `NULL` (no post-filtering).

- blacklist_species:

  Character. Species/build for built-in blacklist (e.g., "hg38", "hg19",
  "mm10", "mm9"), or a path to a custom BED file.

- region_of_interest:

  Character. Path to BED file defining regions of interest (ROI). Only
  loops overlapping these regions will be kept.

- roi_mode:

  Character. How loops must overlap `region_of_interest`. `"any"`
  (default): keep loops where *either* anchor overlaps the ROI (suitable
  for promoter-centric or enhancer-gene queries). `"both"`: keep loops
  where *both* anchors overlap the ROI (suitable for TAD confinement or
  domain-internal interaction queries).

- out_file:

  Character. The file name (including the file path) for saving results
  in the extended BEDPE format.

- write_output:

  Logical. If `TRUE` (default), write the consolidated BEDPE file when
  `out_file` is provided. If `FALSE`, return the `GInteractions` object
  without creating directories or files.

- quiet:

  Logical. If `TRUE`, suppress progress messages while preserving
  warnings. Default: `FALSE`.

## Value

A filtered
[`GInteractions`](https://rdrr.io/pkg/InteractionSet/man/GInteractions-class.html)
object with metadata columns:

- `score`:

  Replicate-balanced consensus score.

- `n_members`:

  Number of raw loops merged into this entry (1 for intersect mode where
  no coordinate merging occurs).

- `n_reps`:

  Number of input files that support this entry.

When `write_output = TRUE` and `out_file` is provided, an extended BEDPE
file is written with the additional columns `n_members` and `n_reps`
appended after the standard BEDPE fields.

## Examples

``` r
# 1. Get paths to example BEDPE files included in the package
f1 <- system.file("extdata", "example_loops_1.bedpe", package = "looplook")
f2 <- system.file("extdata", "example_loops_2.bedpe", package = "looplook")

# 2. Run consolidation (ensure files exist)
# Example A: Intersect Mode
# Only keeps loops present in f1 that are also supported by f2
res_intersect <- consolidate_chromatin_loops(
  files = c(f1, f2),
  mode = "intersect",
  gap = 1000,
  out_file = tempfile(fileext = ".bedpe")
)
#> >>> Reading BEDPE files
#>     File 1: 300 loops
#>     File 2: 300 loops
#> >>> Intersect mode: Reference-based filtering (No Coordinate Merging)
#>     Base: File 1. Criterion: Must overlap with ALL other files.
#>     Intersecting with File 2...
#> Finished! Saved to /tmp/RtmpjNEYlw/file9d2f383429c2.bedpe
#> Finished! Final loops: 12

# Example B: Consensus Mode (formerly Reproducible)
# Finds consensus loops supported by both replicates (default for N=2)
res_consensus <- consolidate_chromatin_loops(
  files = c(f1, f2),
  mode = "consensus",
  gap = 1000,
  out_file = tempfile(fileext = ".bedpe")
)
#> >>> Reading BEDPE files
#>     File 1: 300 loops
#>     File 2: 300 loops
#> >>> Clustering mode (Union/Consensus): Merging coordinates via Graph
#> >>> Consensus mode: Keeping clusters in >= 2 replicates
#> Finished! Saved to /tmp/RtmpjNEYlw/file9d2f46041b05.bedpe
#> Finished! Final loops: 12

# Example C: Union Mode
# Merges all loops into a single map
res_union <- consolidate_chromatin_loops(
  files = c(f1, f2),
  mode = "union",
  gap = 1000,
  out_file = tempfile(fileext = ".bedpe")
)
#> >>> Reading BEDPE files
#>     File 1: 300 loops
#>     File 2: 300 loops
#> >>> Clustering mode (Union/Consensus): Merging coordinates via Graph
#> >>> Union mode: Keeping all clusters
#> Finished! Saved to /tmp/RtmpjNEYlw/file9d2f2ba638d8.bedpe
#> Finished! Final loops: 586

# Example D: Dual Filtering Strategy (Recommended for HiChIP)
# 1. Pre-filter: Discard singletons (score < 2) to remove noise.
# 2. Merge: Find loops present in both replicates.
# 3. Post-filter: Keep only strong consensus loops (score > 5).
res_clean <- consolidate_chromatin_loops(
  files = c(f1, f2),
  mode = "consensus",
  min_raw_score = 2, # Pre-filter (remove noise)
  min_score = 5, # Post-filter (keep strong loops)
  gap = 1000,
  out_file = tempfile(fileext = ".bedpe")
)
#> >>> Reading BEDPE files
#>     File 1: 115 loops
#>     File 2: 100 loops
#> >>> Clustering mode (Union/Consensus): Merging coordinates via Graph
#> >>> Consensus mode: Keeping clusters in >= 2 replicates
#> Finished! Saved to /tmp/RtmpjNEYlw/file9d2f2e4fcaed.bedpe
#> Finished! Final loops: 4

# Inspect results
length(res_intersect)
#> [1] 12
length(res_clean)
#> [1] 4
```
