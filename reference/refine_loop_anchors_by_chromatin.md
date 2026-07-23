# Chromatin-Aware refinement of loop anchor classification

Refines loop anchor types (P, E, G, eP, eG) using chromatin mark
evidence (H3K4me1, H3K4me3, and optionally H3K27ac, ATAC, H3K27me3).
Anchors with dual promoter-enhancer chromatin signatures are flagged as
`"dual"`. **Supported multi-omic workflow:** When both expression and
chromatin data are used, the required order is:


      annotate_peaks_and_loops()
      ->refine_loop_anchors_by_expression()
      ->refine_loop_anchors_by_chromatin()

Chromatin-only refinement (without prior expression refinement) is also
supported.

**Unsupported:** Running
[`refine_loop_anchors_by_expression()`](https://zying106.github.io/looplook/reference/refine_loop_anchors_by_expression.md)
on a chromatin-refined object will raise an error. The reverse order
(chromatin -\>expression) is not supported because `clean_anchor()` only
handles P and G types and cannot correctly process `dual` or other types
set by chromatin refinement.

## Usage

``` r
refine_loop_anchors_by_chromatin(
  annotation_res,
  chromatin_beds = list(),
  anchor_gap = 200,
  anchor_min_overlap = 100,
  species = "hg38",
  tss_region = NULL,
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
  enhancer_bed = NULL,
  hub_metric = c("unique_contacts", "total_loops"),
  allow_rerefine = FALSE
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

  Integer. Candidate search radius in bp between an anchor and a
  chromatin mark peak. Default `200`. Does *not* bypass the
  physical-overlap requirement enforced by `anchor_min_overlap`.

- anchor_min_overlap:

  Integer. Minimum required physical overlap (bp) between an anchor and
  a chromatin mark or enhancer BED region. Default `100`. Proximity-only
  matches (within `anchor_gap` but without actual base-pair
  intersection) are excluded. Reduce to `1--10` for narrow features
  (TFBS, DHS, SNPs).

- species:

  Character. Genome assembly. Default `"hg38"`.

- tss_region:

  Numeric vector of length 2, or `NULL`. Promoter window around the
  anchor centre in bp, used for TSS re-annotation of E-\>P / G-\>P
  anchors. Default: `NULL` (auto-inherits from upstream annotation
  metadata; falls back to `c(-2000, 2000)` if neither is available). A
  warning is emitted if an explicit value differs from the upstream
  annotation.

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
  `NULL` (default): uses all five types `c("P","G","E","eP","eG")`
  regardless of input source. Chromatin evidence can reclassify any
  anchor type (e.g. P-\>dual, E-\>P), so restricting to
  expression-silenced types alone would miss biologically important
  transitions. Set to a subset to limit reclassification scope.

- recompute_targets:

  Logical. If `TRUE` (default), re-run target gene assignment using
  updated anchor types. Requires the input `annotation_res` to contain
  the `looplook_anchor_state` attribute (present when using
  [`annotate_peaks_and_loops`](https://zying106.github.io/looplook/reference/annotate_peaks_and_loops.md)
  output). Default `TRUE`.

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
  to distinguish dual-signature elements (ratio \>= 3) from
  promoter-proximal H3K4me1 shoulders (ratio \< 3, reclassified to P).
  Requires the rtracklayer package. The list must include named elements
  `"H3K4me1"` and `"H3K4me3"`. Default: `NULL` (BED-only mode –
  H3K4me3-positive anchors lacking H3K4me1 overlap are reclassified as
  P; dual-positive (H3K4me1+/H3K4me3+) anchors retain their original
  type pending bigWig resolution). **Strongly recommended** when
  ChIP-seq data are available: bigWig signals provide the quantitative
  H3K4me1/H3K4me3 ratio needed to distinguish dual-signature elements
  from promoter-proximal H3K4me1 shoulders. Use the companion script
  `inst/scripts/diagnose_h3k4me_ratio.R` to explore the ratio
  distribution in your data and choose an appropriate threshold.
  **Important:** The H3K4me1/H3K4me3 ratio is only interpretable when
  both bigWig tracks are generated with comparable normalization (e.g.,
  both CPM, both RPGC), identical scaling, binning, and smoothing, from
  matched biological material. Do not mix tracks with different
  normalization schemes (e.g., one fold-enrichment and one raw coverage)
  or from different cell types or batches. A separate low-signal
  threshold is not applied because input peaks already pass a peak
  caller (e.g. MACS2) and represent statistically enriched regions;
  signal at called peaks is expected to be well above noise.

- bw_ratio_threshold:

  Numeric. Minimum H3K4me1/H3K4me3 ratio to classify a dual-positive
  anchor as true `"dual"`. Anchors below this threshold are reclassified
  as `"P"` (promoter shoulder). Default: `3` (H3K4me1 signal must be at
  least 3x H3K4me3 to override promoter identity). Only used when
  `chromatin_bw` is provided.

- enhancer_bed:

  Character or `NULL`. Path to a BED file of known enhancer regions
  (e.g., FANTOM5, ENCODE cCREs). Anchors overlapping these regions
  receive high-confidence enhancer classification: `"E"` when H3K4me3 is
  absent, `"dual"` when H3K4me3 is also present. This curated evidence
  takes priority over chromatin-mark-derived enhancer evidence levels.
  Default: `NULL`.

- hub_metric:

  Character. Which connectivity count to use for hub detection.
  `"unique_contacts"` (default) or `"total_loops"`.

- allow_rerefine:

  Logical. If `FALSE` (default), stops when the input has already
  undergone chromatin refinement. Repeated chromatin refinement is not
  guaranteed to be idempotent (anchor types and TSS assignments have
  already been modified). Set `TRUE` only for exploratory use and
  restart from a common upstream object for formal parameter
  comparisons.

## Value

An invisible named list with updated `loop_annotation`,
`anchor_loci_annotation`, `anchor_annotation` (alias of
`anchor_loci_annotation`), `promoter_centric_stats`,
`distal_element_stats`, `chromatin_validation` (includes
`chromatin_state`, `positional_type` and `final_type` columns per
anchor), `plots` (`Chromatin_Dumbbell`: anchor-type before/after
comparison; `Chromatin_Sankey`: anchor-type reclassification flow;
`Chromatin_MarkHeatmap`: percentage of anchors positive for each mark
per reclassification group; `Chromatin_UpSet`: loop-type UpSet plot (dot
matrix + log10 bar chart)), `plot_list` (alias of `plots`),
`qc_summary`, `metadata`, and when `recompute_targets = TRUE`,
chromatin-aware `target_annotation` and `target_gene_links`. Set
`recompute_targets = FALSE` to preserve pre-chromatin target
assignments.

## Details

**Reclassification rules (minimum input: H3K4me1 + H3K4me3):**

- Enhancer BED overlap (highest priority): overlapped anchors become
  `"E"` (H3K4me3 absent) or `"dual"` (H3K4me3 present).

- P + H3K4me1(+) and H3K4me3(+) -\> `"dual"` (dual-signature) or `"P"`
  (when H3K4me1 is likely a promoter shoulder; this is resolved by
  bigWig ratio \>= 3 when `chromatin_bw` is provided).

- eP/eG + canonical or strong + active/intermediate enhancer chromatin
  -\> `"E"` (enhancer identity confirmed).

- P + H3K4me1(+) and H3K4me3(-) *and* (H3K27ac(+) or ATAC(+)) -\> `"E"`
  (conservative: requires active-mark confirmation beyond H3K4me1
  alone).

- eP/eG + canonical or strong + active/intermediate enhancer chromatin
  -\> `"E"` (enhancer identity confirmed; matches the P-\>E rule to
  guarantee the same outcome whether expression or chromatin refinement
  runs first).

- eP/eG with H3K4me3(+) but without H3K4me1(+) -\> `"P"` (revert to
  promoter). This includes H3K4me3(+)/H3K27me3(+) bivalent anchors.

- Triple-positive H3K4me1(+)/H3K4me3(+)/H3K27me3(+) eP/eG anchors are
  classified as `conflicting_marks` and retain their previous eP/eG
  type. Peak co-occurrence alone is not sufficient evidence for
  restoration to a functional promoter state.

- E + H3K4me1(+) and H3K4me3(+) -\> resolved by bigWig ratio:
  me1-dominant -\> `"dual"`, unresolved (no bigWig) -\> keep `"E"`,
  not-me1-dominant -\> `"P"`.

- E + H3K4me3(+) without H3K4me1 -\> `"P"` (unannotated promoter).

- G + H3K4me1(+) and H3K4me3(+) -\> same three-way bigWig resolution as
  E above.

- G + H3K4me3(+) without H3K4me1 -\> `"P"` (internal promoter).

- G + H3K4me1(+) and H3K4me3(-) *and* (H3K27ac(+) or ATAC(+)) -\> `"E"`
  (conservative intronic enhancer).

- All other anchors: unchanged.

**Required inputs:** `H3K4me1` and `H3K4me3` are required. Their BED
files must exist and contain at least one peak (an empty required-mark
file cannot be interpreted as genome-wide negative evidence). Optional
marks (`H3K27ac`, `ATAC`, `H3K27me3`) that are not provided are recorded
as `NA` (not_tested), and rules that depend on their absence are skipped
conservatively.

**Chromatin state inference:** Each anchor is also assigned a
`chromatin_state` based on its mark combination (highest-priority match
first):

- `"conflicting_marks"`: H3K27me3+ coexisting with any active mark
  (H3K4me1+, H3K27ac+, H3K4me3+, or ATAC+) – bivalent/poised/ambiguous
  chromatin; takes priority over enhancer-like and dual-like
  classifications.

- `"dual_like"`: H3K4me1+ H3K4me3+

- `"active_enhancer_like"`: H3K4me1+ H3K27ac+ ATAC+

- `"intermediate_enhancer_like"`: H3K4me1+ with H3K27ac+ or ATAC+

- `"other_enhancer_like"`: H3K4me1+ only, or H3K27ac+ only

- `"promoter_like"`: H3K4me3+

- `"repressed"`: H3K27me3+ (no active marks)

- `"uncertain"`: marks tested but none present

- `"no_data"`: no marks tested

**Candidate anchor selection:** When `candidate_types = NULL` (default),
all five anchor types `c("P","G","E","eP","eG")` are tested regardless
of input source. Set `candidate_types` to a subset (e.g.
`c("eP", "eG")`) to restrict reclassification scope.

After reclassification, `loop_type` is recomputed from the updated
anchor types.

**Bivalent / conflicting-mark domains:** Anchors with H3K27me3+
coexisting alongside active marks (H3K4me3+, H3K4me1+, or H3K27ac+)
receive `chromatin_state = "conflicting_marks"`. This is common in
stem-cell and developmental contexts, where H3K4me3+/H3K27me3+ bivalent
domains poise genes for activation while maintaining transcriptional
silence. The reclassification rules distinguish between: bivalent
double-positive anchors (H3K4me3+ / H3K27me3+, H3K4me1-), which are
restored to `P` while retaining the `conflicting_marks` state; and
triple-positive anchors (H3K4me1+ / H3K4me3+ / H3K27me3+), which retain
their previous eP/eG type because the co-presence of three marks on bulk
ChIP does not resolve whether H3K4me1 and H3K4me3 come from the same
nucleosome, allele, or cell. To distinguish active from poised
promoters, run expression refinement before chromatin refinement and use
the retained expression state together with
`chromatin_state = "conflicting_marks"`. The `chromatin_state` column in
`chromatin_validation` remains available for independent filtering of
conflicting-mark anchors.

**Pipeline guidance:** The two refinement modules serve distinct roles
in a complete analysis:

- **Expression-aware refinement** evaluates transcriptional activity (is
  the gene expressed?) and produces expression-filtered
  `Putative_Target_Genes` and `Promoter_Target_Genes` plus
  `Retained_In_Functional_Network` and `Refinement_Action` for
  downstream profiling.

- **Chromatin-aware refinement** evaluates chromatin identity (what *is*
  this anchor?) using direct histone mark evidence, correcting anchor
  types and reclassifying eP/eG into definitive categories (E, P, G,
  dual).

When you have **both RNA-seq and ChIP-seq**, use
`refine_loop_anchors_by_expression` first, then
`refine_loop_anchors_by_chromatin` – expression marks silent anchors as
hypotheses (eP/eG), chromatin confirms or corrects them into definitive
types, and the expression-filtered target gene columns
(`Putative_Target_Genes`, `Promoter_Target_Genes` etc.) carry the
current target set into downstream profiling. When you have **only
ChIP-seq**, run chromatin refinement directly on the output of
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
#> Chromatin Sankey skipped: install 'networkD3' and 'htmlwidgets' packages.
#> Chromatin MarkHeatmap skipped: install 'ComplexHeatmap' and 'circlize' packages.
#>   chromatin role: 1 unchanged E anchor(s) blocked from strict target (nearest-gene not loop-supported)
#>   chromatin role: 1 anchor(s) assigned enhancer_candidate role
#>   chromatin role: 1 anchor(s) excluded from strict assignment (positional_candidate)
table(res_chromatin$loop_annotation$loop_type)
#> 
#> E-E E-P G-G G-P P-P 
#>   1   2   1   3   3 
```
