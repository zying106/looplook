# Chromatin-Aware refinement of loop anchor classification

Refines loop anchor types (P, E, G, eP, eG) using chromatin mark
evidence (H3K4me1, H3K4me3, and optionally H3K27ac, ATAC, H3K27me3).
Anchors with dual promoter-enhancer chromatin signatures are flagged as
`"dual"`. This function can be called before or after
[`refine_loop_anchors_by_expression`](https://zying106.github.io/looplook/reference/refine_loop_anchors_by_expression.md)
– the two refinements address orthogonal biological questions (chromatin
identity vs. transcriptional activity) and do not depend on ordering.

## Usage

``` r
refine_loop_anchors_by_chromatin(
  annotation_res,
  chromatin_beds = list(),
  anchor_gap = 200,
  anchor_min_overlap = 100,
  species = "hg38",
  out_dir = "./results/chromatin",
  project_name = "HiChIP",
  color_palette = "Paired",
  candidate_types = NULL,
  recompute_targets = TRUE,
  write_output = TRUE,
  quiet = FALSE,
  sankey_colors = NULL,
  chromatin_bw = NULL,
  bw_ratio_threshold = 3,
  enhancer_bed = NULL
)
```

## Arguments

- annotation_res:

  List. Output from
  [`annotate_peaks_and_loops`](https://zying106.github.io/looplook/reference/annotate_peaks_and_loops.md)
  or
  [`refine_loop_anchors_by_expression`](https://zying106.github.io/looplook/reference/refine_loop_anchors_by_expression.md).

- chromatin_beds:

  Named list of BED file paths. At minimum `H3K4me1` and `H3K4me3` are
  required for meaningful reclassification; `H3K27ac`, `ATAC`, and
  `H3K27me3` improve confidence but are optional.

- anchor_gap:

  Integer. Gap tolerance for mark overlap. Default `200`.

- anchor_min_overlap:

  Integer. Minimum overlap for mark overlap. Default `100`.

- species:

  Character. Genome assembly. Default `"hg38"`.

- out_dir:

  Character. Output directory for Excel export. Default
  `"./results/chromatin"`.

- project_name:

  Character. Project prefix. Default `"HiChIP"`.

- color_palette:

  Character. RColorBrewer palette name for the rose chart. Dumbbell and
  mark-enrichment heatmap use fixed academic palettes. Default:
  `"Paired"`.

- candidate_types:

  Character vector or `NULL`. Anchor types to validate and reclassify.
  `NULL` (default): auto-selects `c("eP","eG")` for refined input,
  `c("P","G","E")` for raw. Set explicitly to `c("P","G","E","eP","eG")`
  for full-range analysis.

- recompute_targets:

  Logical. If `TRUE` (default), re-run target gene assignment using
  updated anchor types. Requires the input `annotation_res` to contain
  the `looplook_anchor_state` attribute (present when using
  [`annotate_peaks_and_loops`](https://zying106.github.io/looplook/reference/annotate_peaks_and_loops.md)
  output). Default `FALSE`.

- write_output:

  Logical. Write Excel workbook. Default `TRUE`.

- quiet:

  Logical. Suppress messages. Default `FALSE`.

- sankey_colors:

  Named character vector or `NULL`. Override the default type-to-color
  mapping in the Sankey diagram. Names must be anchor types (`"P"`,
  `"E"`, `"G"`, `"eP"`, `"eG"`, `"dual"`), values are hex colours. When
  `NULL` (default), the colourblind-safe Wong palette is used:
  `E = "#E69F00"` (orange), `dual = "#CC0000"` (red), `P = "#0072B2"`
  (blue), `eP = "#009E73"` (bluish-green), `G = "#CC79A7"`
  (reddish-purple), `eG = "#56B4E9"` (sky blue).

- chromatin_bw:

  Named list of bigWig file paths, or `NULL`. When provided,
  H3K4me1/H3K4me3 signal ratios are computed for dual-positive anchors
  to distinguish true dual-function elements (ratio \>= 3) from
  promoter-proximal H3K4me1 shoulders (ratio \< 3, reclassified to P).
  Requires the rtracklayer package. The list must include named elements
  `"H3K4me1"` and `"H3K4me3"`. Default: `NULL` (BED-only mode — all
  H3K4me3-positive anchors are classified as P). **Strongly
  recommended** when ChIP-seq data are available: bigWig signals provide
  the quantitative H3K4me1/H3K4me3 ratio needed to distinguish true
  dual-function elements from promoter-proximal H3K4me1 shoulders. Use
  the companion script `inst/scripts/diagnose_h3k4me_ratio.R` to explore
  the ratio distribution in your data and choose an appropriate
  threshold.

- bw_ratio_threshold:

  Numeric. Minimum H3K4me1/H3K4me3 ratio to classify a dual-positive
  anchor as true `"dual"`. Anchors below this threshold are reclassified
  as `"P"` (promoter shoulder). Default: `3` (H3K4me1 signal must be at
  least 3× H3K4me3 to override promoter identity). Only used when
  `chromatin_bw` is provided.

- enhancer_bed:

  Character or `NULL`. Path to a BED file of known enhancer regions
  (e.g., FANTOM5, ENCODE cCREs). Anchors overlapping these regions
  receive high-confidence enhancer classification: `"E"` when H3K4me3 is
  absent, `"dual"` when H3K4me3 is also present. This curated evidence
  takes priority over chromatin-mark-derived confidence levels. Default:
  `NULL`.

## Value

An invisible named list with updated `loop_annotation`,
`anchor_loci_annotation`, `promoter_centric_stats`,
`distal_element_stats`, `chromatin_validation`, `plots`
(`Chromatin_Dumbbell`: anchor-type before/after comparison;
`Chromatin_UpSet`: loop-type UpSet plot (dot matrix + log10 bar chart)),
`plot_list` (alias of `plots`), `qc_summary`, and `metadata`. When
`recompute_targets = TRUE` (default), target links are rebuilt from
chromatin-updated anchor states, producing chromatin-aware
`target_annotation` and `target_gene_links`. Set to `FALSE` to preserve
pre-chromatin target assignments and use `target_source = "loops"` for
downstream profiling. producing chromatin-aware `target_annotation` and
`target_gene_links`. The input must carry the `looplook_anchor_state`
attribute (present when using
[`annotate_peaks_and_loops`](https://zying106.github.io/looplook/reference/annotate_peaks_and_loops.md)
output).

## Details

**Reclassification rules (minimum input: H3K4me1 + H3K4me3):**

- Enhancer BED overlap (highest priority): overlapped anchors become
  `"E"` (H3K4me3 absent) or `"dual"` (H3K4me3 present).

- P + H3K4me1(+) and H3K4me3(+) -\> `"dual"` (dual-function) or `"P"`
  (when H3K4me1 is likely a promoter shoulder; this is resolved by
  bigWig ratio \>= 3 when `chromatin_bw` is provided).

- eP/eG + gold_standard or high_confidence + active/primed enhancer
  chromatin -\> `"E"` (enhancer identity confirmed).

- P + H3K4me1(+) and H3K4me3(-) *and* (H3K27ac(+) or ATAC(+)) -\> `"E"`
  (conservative: requires active-mark confirmation beyond H3K4me1
  alone).

- eP/eG + gold_standard or high_confidence + active/primed enhancer
  chromatin -\> `"E"` (enhancer identity confirmed; matches the P-\>E
  rule to guarantee the same outcome whether expression or chromatin
  refinement runs first).

- eP + promoter_like (H3K4me3+, irrespective of H3K27me3 status) -\>
  `"P"` (revert).

- eG + promoter_like -\> `"G"` (revert).

- E + H3K4me1(+) and H3K4me3(+) -\> `"dual"` (dual-function).

- E + H3K4me3(+) without H3K4me1 -\> `"P"` (unannotated promoter).

- G + H3K4me1(+) and H3K4me3(+) -\> `"dual"` (gene-body dual).

- G + H3K4me3(+) without H3K4me1 -\> `"P"` (internal promoter).

- G + H3K4me1(+) and H3K4me3(-) *and* (H3K27ac(+) or ATAC(+)) -\> `"E"`
  (conservative intronic enhancer).

- All other anchors: unchanged.

**Important – H3K4me3 dependency:** Rules that test for H3K4me3
*absence* (e.g., P -\> E, G -\> E) require H3K4me3 to be provided and
explicitly called as absent at the anchor. When `H3K4me3` is not
included in `chromatin_beds`, H3K4me3 is `NA` for all anchors, and these
rules are skipped entirely (conservative: no reclassification without
explicit data). Include `H3K4me3` to enable the full set of
reclassification rules.

**Chromatin state inference:** Each anchor is also assigned a
`chromatin_state` based on its mark combination (highest-priority match
first):

- `"conflicting_marks"`: H3K27me3+ coexisting with any active mark
  (H3K4me1+, H3K27ac+, or H3K4me3+) – bivalent/poised/ambiguous
  chromatin; takes priority over enhancer-like and dual-like
  classifications.

- `"dual_like"`: H3K4me1+ H3K4me3+

- `"active_enhancer_like"`: H3K4me1+ H3K27ac+ ATAC+

- `"primed_enhancer_like"`: H3K4me1+ with H3K27ac+ or ATAC+

- `"weak_enhancer_like"`: H3K4me1+ only, or H3K27ac+ only

- `"promoter_like"`: H3K4me3+

- `"repressed"`: H3K27me3+ (no active marks)

- `"uncertain"`: marks tested but none present

- `"no_data"`: no marks tested

**Candidate anchor selection:** When `candidate_types = NULL` (default),
raw annotation from
[`annotate_peaks_and_loops`](https://zying106.github.io/looplook/reference/annotate_peaks_and_loops.md)
tests `P`, `G`, and `E` anchors. Expression-refined input tests `eP` and
`eG` anchors.

After reclassification, `loop_type` is recomputed from the updated
anchor types.

**Pipeline guidance:** The two refinement modules serve distinct roles
in a complete analysis:

- **Expression-aware refinement** evaluates transcriptional activity (is
  the gene expressed?) and produces activity-aware metrics
  (`Active_Target_Genes`, `Retained_In_Functional_Network`,
  `Refinement_Action`) required for downstream GSEA, GO enrichment, and
  differential expression profiling.

- **Chromatin-aware refinement** evaluates chromatin identity (what *is*
  this anchor?) using direct histone mark evidence, correcting anchor
  types and reclassifying eP/eG into definitive categories (E, P, G,
  dual).

When you have **both RNA-seq and ChIP-seq**, use
`refine_loop_anchors_by_expression` first, then
`refine_loop_anchors_by_chromatin` – expression marks silent anchors as
hypotheses (eP/eG), chromatin confirms or corrects them into definitive
types, and the expression-derived activity columns
(`Active_Target_Genes` etc.) pass through into downstream profiling.
When you have **only ChIP-seq**, run chromatin refinement directly on
the output of
[`annotate_peaks_and_loops`](https://zying106.github.io/looplook/reference/annotate_peaks_and_loops.md).

## See also

[`annotate_peaks_and_loops`](https://zying106.github.io/looplook/reference/annotate_peaks_and_loops.md)
for initial 3D annotation,
[`refine_loop_anchors_by_expression`](https://zying106.github.io/looplook/reference/refine_loop_anchors_by_expression.md)
for expression-aware refinement,
[`validate_epeG_by_chromatin`](https://zying106.github.io/looplook/reference/validate_epeG_by_chromatin.md)
for standalone chromatin validation.

## Examples

``` r
# 1. Get paths to the required example files in the package
rdata_path <- system.file("extdata", "analysis_results.RData", package = "looplook")
k4me1_path <- system.file("extdata", "example_h3k4me1_peaks.bed", package = "looplook")
k4me3_path <- system.file("extdata", "example_h3k4me3_peaks.bed", package = "looplook")

# 2. Load pre-computed annotation result
temp_env <- new.env()
load(rdata_path, envir = temp_env)
raw_annotation <- temp_env[[ls(temp_env)[1]]]
raw_annotation$loop_annotation <- head(raw_annotation$loop_annotation, 10)
raw_annotation$target_annotation <- head(raw_annotation$target_annotation, 3)

# 3. Run chromatin-aware refinement
res_chromatin <- refine_loop_anchors_by_chromatin(
    annotation_res = raw_annotation,
    chromatin_beds = list(
        H3K4me1 = k4me1_path,
        H3K4me3 = k4me3_path
    ),
    anchor_gap = 500,
    anchor_min_overlap = 100,
    species = "hg38",
    out_dir = tempdir(),
    project_name = "Example_Chromatin",
    write_output = FALSE,
    quiet = TRUE
)
table(res_chromatin$loop_annotation$loop_type)
#> 
#> E-E E-P G-G G-P P-P 
#>   1   2   1   1   5 
```
