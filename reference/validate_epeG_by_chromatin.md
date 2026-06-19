# Orthogonal validation of eP/eG reclassification by chromatin marks

Validates the expression-aware eP/eG reclassification produced by
[`refine_loop_anchors_by_expression`](https://zying106.github.io/looplook/reference/refine_loop_anchors_by_expression.md)
using orthogonal chromatin data (ATAC-seq, ChIP-seq). Each eP/eG anchor
is tested for overlap with user-supplied mark BED files and assigned a
confidence level based on the ENCODE active-enhancer signature.

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

  List. Output from
  [`annotate_peaks_and_loops`](https://zying106.github.io/looplook/reference/annotate_peaks_and_loops.md)
  or
  [`refine_loop_anchors_by_expression`](https://zying106.github.io/looplook/reference/refine_loop_anchors_by_expression.md).
  When the refined output is provided, all anchors classified as `eP` or
  `eG` are validated. When the raw annotation is provided, all `P`, `G`,
  and `E` anchors are evaluated (chromatin evidence may reclassify E to
  P or dual). Anchors are tested for overlap with chromatin mark BED
  files and assigned confidence levels based on ENCODE active-enhancer
  criteria.

- chromatin_beds:

  Named list of BED file paths. Names must be mark identifiers: any of
  `"H3K4me1"`, `"H3K27ac"`, `"ATAC"`, `"H3K27me3"`, `"H3K4me3"`.
  Additional names are ignored with a warning.

- anchor_gap:

  Integer. Maximum gap (bp) between an anchor and a chromatin peak for
  them to be considered overlapping. Default: `200`.

- anchor_min_overlap:

  Integer. Minimum overlap (bp) required. Default: `100`.

- candidate_types:

  Character vector or `NULL`. Candidate anchor types to validate. `NULL`
  (default): `c("eP","eG")` for refined input, `c("P","G","E")` for raw.
  Set to `c("P","G","E","eP","eG")` to cover all positional categories.

- species:

  Character. Genome assembly for seqlevel harmonization. Default:
  `"hg38"`.

- quiet:

  Logical. Suppress progress messages. Default: `FALSE`.

## Value

A data frame with one row per candidate anchor (eP/eG when the input is
refined, P/G when the input is raw from `annotate_peaks_and_loops`):

- `anchor_id`:

  Anchor identifier.

- `chr`, `start`, `end`:

  Anchor coordinates.

- `anchor_type`:

  Original type (P or G) before reclassification.

- `anchor_gene`:

  Gene symbol(s) at this anchor.

- `cluster_id`:

  Loop cluster identifier.

- `H3K4me1`, `H3K27ac`, `ATAC`:

  Logical. TRUE if the anchor overlaps the corresponding positive mark.

- `H3K27me3`, `H3K4me3`:

  Logical. TRUE if the anchor overlaps the corresponding negative mark.

- `confidence`:

  Factor with levels: `"gold_standard"`, `"high_confidence"`,
  `"supported"`, `"weak"`, `"uncertain"`.

- `evidence`:

  Human-readable summary of which marks supported the classification.

## Details

**Confidence levels (ENCODE active-enhancer criteria):**

- `"gold_standard"`: All five marks align with the canonical
  active-enhancer signature: H3K4me1(+), H3K27ac(+), ATAC(+),
  H3K4me3(-), H3K27me3(-). Requires all five marks to be provided.

- `"high_confidence"`: H3K4me1(+) and H3K27ac(+), or H3K4me1(+) and
  ATAC(+). At least two supporting marks present.

- `"supported"`: At least one enhancer-associated mark (H3K4me1,
  H3K27ac, or ATAC) is present.

- `"weak"`: Only exclusion marks are informative (H3K27me3(-) or
  H3K4me3(-)) without positive enhancer evidence.

- `"uncertain"`: No chromatin data provided or no overlaps detected.

**Mark semantics:**

- *Positive marks* (presence = supporting evidence): `H3K4me1`,
  `H3K27ac`, `ATAC`.

- *Negative marks* (absence = supporting evidence): `H3K27me3`,
  `H3K4me3`.

## Examples

``` r
# Load pre-computed annotation result
rdata_path <- system.file("extdata", "analysis_results.RData",
                          package = "looplook")
temp_env <- new.env()
load(rdata_path, envir = temp_env)
raw_annotation <- temp_env[[ls(temp_env)[1]]]

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
table(result$confidence)
#> 
#>   gold_standard high_confidence       supported            weak       uncertain 
#>               0               0               0               0             385 
```
