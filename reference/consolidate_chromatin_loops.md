# Consolidate and Integrate Chromatin Loops from Multiple Sources

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

It also supports a **two-stage filtering strategy** to maximize
signal-to-noise ratio:

- **Pre-filtering** (`min_raw_score`): Removes low-confidence noise
  (e.g., singleton reads) from raw files *before* merging.

- **Post-filtering** (`min_score`): Filters the final consensus loops
  based on their aggregated score (e.g., average intensity).

## Usage

``` r
consolidate_chromatin_loops(
  files = NULL,
  gap = 1000,
  mode = c("consensus", "intersect", "union"),
  min_consensus = NULL,
  min_raw_score = NULL,
  min_score = NULL,
  blacklist_species = NULL,
  region_of_interest = NULL,
  out_file = NULL
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

  - For 2 replicates: Requires both (2/2).

  - For \>2 replicates: Requires strict majority (\>75\\

  - (e.g., 3 for N=3, 4 for N=4, 4 for N=5).

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

  - For `"consensus"` mode, this filters the consensus loops based on
    their representative score.

  - Default: `NULL` (no post-filtering).

- blacklist_species:

  Character. Species/build for built-in blacklist (e.g., "hg38", "hg19",
  "mm10", "mm9"), or a path to a custom BED file.

- region_of_interest:

  Character. Path to BED file. Only loops overlapping these regions will
  be kept.

- out_file:

  Character. The file name (including the file path) for saving results
  in the extended BEDPE format.

## Value

A filtered
[`GInteractions`](https://rdrr.io/pkg/InteractionSet/man/GInteractions-class.html)
object.

## Examples

``` r
# 1. Get paths to example BEDPE files included in the package
f1 <- system.file("extdata", "example_loops_1.bedpe", package = "looplook")
f2 <- system.file("extdata", "example_loops_2.bedpe", package = "looplook")

# 2. Run consolidation (ensure files exist)
if (f1 != "" && f2 != "") {
  # Example A: Intersect Mode
  # Only keeps loops present in f1 that are also supported by f2
  res_intersect <- consolidate_chromatin_loops(
    files = c(f1, f2),
    mode = "intersect",
    gap = 1000,
    out_file = tempfile(fileext = ".bedpe")
  )

  # Example B: Consensus Mode (formerly Reproducible)
  # Finds consensus loops supported by both replicates (default for N=2)
  res_consensus <- consolidate_chromatin_loops(
    files = c(f1, f2),
    mode = "consensus",
    gap = 1000,
    out_file = tempfile(fileext = ".bedpe")
  )

  # Example C: Union Mode
  # Merges all loops into a single map
  res_union <- consolidate_chromatin_loops(
    files = c(f1, f2),
    mode = "union",
    gap = 1000,
    out_file = tempfile(fileext = ".bedpe")
  )

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

  # Inspect results
  length(res_intersect)
  length(res_clean)
}
#> >>> Reading BEDPE files
#>     File 1: 300 loops
#>     File 2: 300 loops
#> >>> Intersect mode: Reference-based filtering (No Coordinate Merging)
#>     Base: File 1. Criterion: Must overlap with ALL other files.
#>     Intersecting with File 2...
#> Finished! Saved to /tmp/Rtmpp3oEeP/file41ba54cefc97.bedpe
#> Finished! Final loops: 12
#> >>> Reading BEDPE files
#>     File 1: 300 loops
#>     File 2: 300 loops
#> >>> Clustering mode (Union/Consensus): Merging coordinates via Graph
#> >>> Consensus mode: Keeping clusters in >= 2 replicates
#> Finished! Saved to /tmp/Rtmpp3oEeP/file41ba6e8d0d42.bedpe
#> Finished! Final loops: 11
#> >>> Reading BEDPE files
#>     File 1: 300 loops
#>     File 2: 300 loops
#> >>> Clustering mode (Union/Consensus): Merging coordinates via Graph
#> >>> Union mode: Keeping all clusters
#> Finished! Saved to /tmp/Rtmpp3oEeP/file41ba56cb154.bedpe
#> Finished! Final loops: 589
#> >>> Reading BEDPE files
#>     File 1: 115 loops
#>     File 2: 100 loops
#> >>> Clustering mode (Union/Consensus): Merging coordinates via Graph
#> >>> Consensus mode: Keeping clusters in >= 2 replicates
#> Finished! Saved to /tmp/Rtmpp3oEeP/file41ba6c6f4ef5.bedpe
#> Finished! Final loops: 4
#> [1] 4
```
