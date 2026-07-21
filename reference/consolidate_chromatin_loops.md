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

- `"intersect"`: Reference-based filtering. Retains loops in File 1
  whose anchors overlap with loops in every other file within the
  specified `gap` tolerance. Coordinates and scores are inherited from
  File 1 without merging. **Important: Intersect mode is NOT
  symmetric.** The output depends on which file is listed first – loops
  are retained from File 1 only. Changing file order may produce
  different results.

- `"union"`: Retains all chromatin interactions across the entire
  cohort, ideal for exploratory pan-tissue analyses.

**Connected-component chaining**: Graph-based clustering may
transitively chain loci (A-B and B-C merge, pulling A and C into the
same cluster even if they are far apart). By default, a warning is
emitted when any cluster span exceeds the chaining threshold
(`max(3*gap, 5*med_width*gap/1000)`). Use `chaining_policy` to control
this behavior (`"warn"`, `"none"`, `"drop"`, or `"error"`). **Note:**
the default `"warn"` policy does *not* remove affected clusters – their
inflated coordinates (from `min(start)` to `max(end)` across all chained
members) flow into downstream analyses unchanged. When chaining is a
concern for downstream overlapping operations, use
`chaining_policy = "drop"` to exclude wide clusters, or reduce `gap` to
prevent chaining from forming.

It also supports a **two-stage filtering strategy** to improve
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
  chaining_policy = c("warn", "none", "drop", "error"),
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

  Character vector. Paths to BEDPE files (at least two). **Empty
  files:** A zero-byte BEDPE file is treated as an empty replicate and
  triggers a warning. Empty replicates remain part of the replicate set:
  they count toward the default consensus denominator (`min_consensus`)
  and are required by intersect mode (which demands support from all
  inputs). They contribute no interactions in union mode. A file must be
  truly zero-byte to qualify as an empty replicate; any non-empty file
  must be a valid BEDPE with at least six columns. Verify that an empty
  file represents a legitimate biological result rather than a failed
  loop-calling run.

- gap:

  Numeric. Maximum distance (in base pairs) between loop anchors to
  consider them as overlapping. Default: `1000` (1 kb). Typical range:
  500-5000 bp depending on data resolution. Use `gap = 0` for strict
  overlap only.

- mode:

  Character. Choose one of the following: "consensus", "intersect",
  "union". Merge strategy:

  - `"consensus"`: Graph-based clustering to find a consensus set
    supported by a majority of samples (default). Formerly
    "reproducible".

  - `"intersect"`: Strict reference-based filtering (keeps loops in File
    1 supported by ALL other files). The result is **not** symmetric:
    the first input file serves as the reference, and changing file
    order changes results.

  - `"union"`: Merges all detected loops into a comprehensive map.

- min_consensus:

  Integer. Minimum number of replicates a loop must appear in (only
  effective when `mode = "consensus"`). If `NULL` (default), the
  threshold is automatically calculated:

  - For 2-3 replicates: Requires all (100\\

  - For \\N \ge 4\\: Requires \\\ge 75\\\\ of replicates
    (`ceiling(0.75 * n_reps)`; e.g., 3 for N=4, 4 for N=5, 6 for N=8).

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

- chaining_policy:

  Character. Controls behaviour when connected-component chaining
  produces clusters with span exceeding the chaining threshold
  (`max(3*gap, 5*med_width*gap/1000)`). `"warn"` (default): emit a
  warning and retain all clusters *unchanged*, including their
  potentially inflated coordinates. `"none"`: silently accept all
  clusters. `"drop"`: remove clusters exceeding the threshold (safe for
  downstream overlapping analyses that may be affected by inflated
  spans). `"error"`: stop with an error.

- blacklist_species:

  Character. Species/build for built-in ENCODE blacklist (`"hg38"`,
  `"hg19"`, `"mm10"`, `"mm9"`), or a path to a custom BED file. When a
  species name is recognised, the bundled blacklist is used; otherwise
  the value is treated as a file path.

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
  no coordinate merging occurs). For `"intersect"` mode, results are
  **not symmetric** – changing file order changes output.

- `n_reps`:

  Number of input files that support this entry.

- `cluster_id`:

  Connected-component cluster identifier.

The returned object carries a `looplook_metadata` attribute (access via
`attr(x, "looplook_metadata")`) with package version, call parameters,
diagnostics, and database versions. When `write_output = TRUE` and
`out_file` is provided, an extended BEDPE file is written with the
additional columns `n_members` and `n_reps` appended after the standard
BEDPE fields.

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
#>     File 1: 300 loops (raw: 300)
#>     File 2: 300 loops (raw: 300)
#> >>> Intersect mode: Reference-based filtering (No Coordinate Merging)
#>     Base: File 1 (first input). Output coordinates and scores come exclusively from File 1.
#>     Intersecting with File 2...
#> Finished! Final loops: 12
#> Finished! Saved to /tmp/RtmpwD8XlJ/filea06359ef710d.bedpe

# Example B: Consensus Mode (formerly Reproducible)
# Finds consensus loops supported by both replicates (default for N=2)
res_consensus <- consolidate_chromatin_loops(
  files = c(f1, f2),
  mode = "consensus",
  gap = 1000,
  out_file = tempfile(fileext = ".bedpe")
)
#> >>> Reading BEDPE files
#>     File 1: 300 loops (raw: 300)
#>     File 2: 300 loops (raw: 300)
#> >>> Clustering mode (Union/Consensus): Merging coordinates via Graph
#> --- Gap Diagnosis ---
#>   Anchor width: median = 4,859 bp  (IQR: 3,433.25 - 7,146.25 bp)
#>   Adjacent-anchor gap: median = 0 bp  (57% within current gap = 1,000 bp)
#>   Gap / median_anchor_width ratio: 0.2x  (effective width expansion: 1.4x)
#>   Gap appears appropriate for the anchor width distribution.
#>   [i]  Super-enhancer-like anchors -- very broad regulatory domains.
#> --- End Gap Diagnosis ---
#> >>> Consensus mode: Keeping clusters in >= 2 replicates
#> --- Post-Clustering Diagnosis ---
#>   Clusters formed: 12  (from 24 loops surviving consensus)
#>   Members per cluster: median = 2, IQR = 2-2, max = 2
#>   Consensus retention: 24 / 600 input loops (4.0%)
#>   [i]  Low retention -- many loops failed consensus. Gap may be too small for reproducible calls across replicates.
#>   Cluster span: median = 6,260  |  max = 26,481  |  threshold = 5xmed_width(4,859) x gap/1000 = 24,295 bp
#>   Largest cluster spans:
#>     #1: max_span = 26,481 bp, n_members = 2, n_reps = 2
#>     #3: max_span = 22,248 bp, n_members = 2, n_reps = 2
#>     #4: max_span = 11,840 bp, n_members = 2, n_reps = 2
#>   Chaining: 0/12 above threshold -- PASS.
#> --- End Post-Clustering Diagnosis ---
#> Finished! Final loops: 12
#> Finished! Saved to /tmp/RtmpwD8XlJ/filea063758c0301.bedpe

# Example C: Union Mode
# Merges all loops into a single map
res_union <- consolidate_chromatin_loops(
  files = c(f1, f2),
  mode = "union",
  gap = 1000,
  out_file = tempfile(fileext = ".bedpe")
)
#> >>> Reading BEDPE files
#>     File 1: 300 loops (raw: 300)
#>     File 2: 300 loops (raw: 300)
#> >>> Clustering mode (Union/Consensus): Merging coordinates via Graph
#> --- Gap Diagnosis ---
#>   Anchor width: median = 4,859 bp  (IQR: 3,433.25 - 7,146.25 bp)
#>   Adjacent-anchor gap: median = 0 bp  (57% within current gap = 1,000 bp)
#>   Gap / median_anchor_width ratio: 0.2x  (effective width expansion: 1.4x)
#>   Gap appears appropriate for the anchor width distribution.
#>   [i]  Super-enhancer-like anchors -- very broad regulatory domains.
#> --- End Gap Diagnosis ---
#> >>> Union mode: Keeping all clusters
#> --- Post-Clustering Diagnosis ---
#>   Clusters formed: 586  (from 600 loops surviving consensus)
#>   Members per cluster: median = 1, IQR = 1-1, max = 2
#>   Consensus retention: 600 / 600 input loops (100.0%)
#>   Cluster span: median = 6,739  |  max = 27,921  |  threshold = 5xmed_width(4,859) x gap/1000 = 24,295 bp
#>   Largest cluster spans:
#>     #397: max_span = 27,921 bp, n_members = 1, n_reps = 1
#>     #463: max_span = 27,921 bp, n_members = 1, n_reps = 1
#>     #427: max_span = 26,620 bp, n_members = 1, n_reps = 1
#>   Chaining: 0/586 above threshold -- PASS.
#> --- End Post-Clustering Diagnosis ---
#> Finished! Final loops: 586
#> Finished! Saved to /tmp/RtmpwD8XlJ/filea06374d2327b.bedpe

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
#>     File 1: 115 loops (raw: 300)
#>     File 2: 100 loops (raw: 300)
#> >>> Clustering mode (Union/Consensus): Merging coordinates via Graph
#> --- Gap Diagnosis ---
#>   Anchor width: median = 5,403 bp  (IQR: 3,921.5 - 8,977 bp)
#>   Adjacent-anchor gap: median = 7,161.5 bp  (43% within current gap = 1,000 bp)
#>   Gap / median_anchor_width ratio: 0.2x  (effective width expansion: 1.4x)
#>   Gap appears appropriate for the anchor width distribution.
#>   [i]  Super-enhancer-like anchors -- very broad regulatory domains.
#> --- End Gap Diagnosis ---
#> >>> Consensus mode: Keeping clusters in >= 2 replicates
#> --- Post-Clustering Diagnosis ---
#>   Clusters formed: 7  (from 14 loops surviving consensus)
#>   Members per cluster: median = 2, IQR = 2-2, max = 2
#>   Consensus retention: 14 / 215 input loops (6.5%)
#>   [i]  Low retention -- many loops failed consensus. Gap may be too small for reproducible calls across replicates.
#>   Cluster span: median = 5,672  |  max = 22,248  |  threshold = 5xmed_width(5,403) x gap/1000 = 27,015 bp
#>   Largest cluster spans:
#>     #2: max_span = 22,248 bp, n_members = 2, n_reps = 2
#>     #3: max_span = 11,840 bp, n_members = 2, n_reps = 2
#>     #5: max_span = 6,849 bp, n_members = 2, n_reps = 2
#>   Chaining: 0/7 above threshold -- PASS.
#> --- End Post-Clustering Diagnosis ---
#> Finished! Final loops: 4
#> Finished! Saved to /tmp/RtmpwD8XlJ/filea0632b1d75bf.bedpe

# Inspect results
length(res_intersect)
#> [1] 12
length(res_clean)
#> [1] 4
```
