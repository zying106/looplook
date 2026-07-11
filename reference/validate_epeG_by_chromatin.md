# Orthogonal validation of eP/eG reclassification by chromatin marks

Validates the expression-aware eP/eG reclassification produced by
[`refine_loop_anchors_by_expression`](https://zying106.github.io/looplook/reference/refine_loop_anchors_by_expression.md),
or the raw P/G/E classification from
[`annotate_peaks_and_loops`](https://zying106.github.io/looplook/reference/annotate_peaks_and_loops.md),
against user-supplied chromatin mark BED files.

## Usage

``` r
validate_epeG_by_chromatin(
  annotation_res,
  chromatin_beds = list(),
  anchor_gap = 200,
  anchor_min_overlap = 100,
  candidate_types = NULL,
  species = "hg38",
  quiet = FALSE
)
```

## Arguments

- annotation_res:

  List. The raw foundational output object returned by
  [`annotate_peaks_and_loops`](https://zying106.github.io/looplook/reference/annotate_peaks_and_loops.md)
  or the refined output from
  [`refine_loop_anchors_by_expression`](https://zying106.github.io/looplook/reference/refine_loop_anchors_by_expression.md).

- chromatin_beds:

  Named list of BED file paths. Accepts up to five canonical histone
  marks: `H3K4me1`, `H3K27ac`, `ATAC`, `H3K27me3`, and `H3K4me3`.
  Unknown names are dropped with a warning; unmatched case is resolved
  to the canonical names. An empty list classifies every anchor as
  `"uncertain"`.

- anchor_gap:

  Integer. Gap tolerance for mark overlap. Default `200`.

- anchor_min_overlap:

  Integer. Minimum overlap for mark overlap. Default `100`.

- candidate_types:

  Character vector or `NULL`. Anchor types to validate. `NULL`
  (default): uses all five types `c("P","G","E","eP","eG")` regardless
  of input source.

- species:

  Character. Genome assembly string. Default: `"hg38"`. Accepts any
  assembly string; used for seqlevel harmonization.

- quiet:

  Logical. Suppress messages. Default: `FALSE`.

## Value

A data frame with columns:

- anchor_id, chr, start, end, anchor_type, anchor_gene, cluster_id:

  Anchor identifiers.

- H3K4me1, H3K27ac, ATAC, H3K27me3, H3K4me3:

  Logical or `NA`. `TRUE` = overlap, `FALSE` = tested but absent, `NA` =
  not tested.

- enhancer_evidence:

  Factor: `gold_standard` \> `high_confidence` \> `supported` \> `weak`
  \> `uncertain`.

- evidence:

  Human-readable string of the constituting marks (e.g.
  `"H3K4me1+; H3K27ac+; H3K4me3-; H3K27me3-"`). Anchors with
  tested-positive H3K4me3 at `supported` confidence are annotated with a
  `"promoter_like"` tag.

## Details

Confidence levels follow ENCODE active-enhancer criteria. Each anchor is
classified into the highest applicable level:

- gold_standard:

  All 5 marks tested: H3K4me1+, H3K27ac+, ATAC+, H3K27me3-, H3K4me3-.

- high_confidence:

  H3K4me1+ and (H3K27ac+ or ATAC+); H3K4me3 must be absent when tested.

- supported:

  Any of H3K4me1, H3K27ac, or ATAC positive. Anchors with H3K4me3+
  (regardless of H3K27me3 status) receive a `"promoter_like"` evidence
  tag.

- weak:

  Tested negative markers (H3K27me3-, H3K4me3-) present but no active
  marks.

- uncertain:

  Marks tested but none identified; or no marks tested at all.

Marks not provided are recorded as `NA` and never treated as negative
evidence. Case-insensitive name matching allows flexible input (e.g.,
`h3k4me1` normalised to `H3K4me1`).

## Examples

``` r
rdata_path <- system.file("extdata", "analysis_results.RData", package = "looplook")
tmp <- new.env()
load(rdata_path, envir = tmp)
raw_annotation <- tmp[[ls(tmp)[1]]]

# Create dummy chromatin BED files for demonstration
bed_dir <- tempdir()
writeLines("chr6\t10410000\t10413000", file.path(bed_dir, "H3K4me1.bed"))
writeLines("chr6\t10411000\t10414000", file.path(bed_dir, "H3K27ac.bed"))

# Run validation (using raw annotation; pass refined for eP/eG only)
result <- validate_epeG_by_chromatin(
    annotation_res = raw_annotation,
    chromatin_beds = list(
        H3K4me1 = file.path(bed_dir, "H3K4me1.bed"),
        H3K27ac = file.path(bed_dir, "H3K27ac.bed")
    ),
    quiet = TRUE
)
#> Warning: The 2 combined objects have no sequence levels in common. (Use
#>   suppressWarnings() to suppress this warning.)
#> Warning: The 2 combined objects have no sequence levels in common. (Use
#>   suppressWarnings() to suppress this warning.)
table(result$enhancer_evidence)
#> 
#>   gold_standard high_confidence       supported            weak       uncertain 
#>               0               0               0               0             385 
```
