# Expression-Aware refinement of loop anchors and target linkages

Integrates quantitative RNA-seq data (e.g., TPM/FPKM) with 3D structural
data to reclassify regulatory elements and annotate functional status,
deriving a functionally interpretable regulatory network from physical
chromatin contacts. All structural loops are preserved; refinement
status columns indicate which loops belong to the high-confidence active
subset.

## Usage

``` r
refine_loop_anchors_by_expression(
  annotation_res,
  expr_matrix_file = NULL,
  sample_columns = NULL,
  threshold = 1,
  unit_type = "TPM",
  species = "hg38",
  out_dir = "./results/filtered",
  project_name = "HiChIP",
  color_palette = "Paired",
  karyo_bin_size = 1e+05,
  reclassify_by_expression = TRUE,
  hub_percentile = 0.95,
  chromatin_beds = list(),
  write_output = TRUE,
  quiet = FALSE
)
```

## Arguments

- annotation_res:

  List. The raw foundational output object returned by
  [`annotate_peaks_and_loops`](https://zying106.github.io/looplook/reference/annotate_peaks_and_loops.md).

- expr_matrix_file:

  Path to a normalised expression matrix (TPM/FPKM, genes x samples).
  Required for refinement. Default: `NULL`.

- sample_columns:

  Character vector or integer indices. Columns in `expr_matrix_file` to
  average. Default: `NULL`.

- threshold:

  Numeric. Minimum expression (e.g. TPM \>= 1) for a gene to be
  considered active. Default: `1`.

- unit_type:

  Character. Expression unit for plot labels (e.g., `"TPM"`). Default:
  `"TPM"`.

- species:

  Character. Genome assembly. One of `"hg38"`, `"hg19"`, `"mm10"`,
  `"mm9"`. Default: `"hg38"`. Extensible to other species via
  `species_txdb_pkg()` and related helpers.

- out_dir:

  Character. Output directory for the Excel results file. Default:
  `"./results/filtered"`.

- project_name:

  Character. Prefix for output files (automatically appends
  `"_Filtered"`). Default: `"HiChIP"`.

- color_palette:

  Character. RColorBrewer palette name for loop-type colour assignments.
  Default: `"Paired"`.

- karyo_bin_size:

  Integer. Bin width in bp for karyotype heatmaps. Default: `1e5`.

- reclassify_by_expression:

  Logical. If `TRUE` (default), silent promoters (P) and gene bodies (G)
  are reclassified as eP/eG.

- hub_percentile:

  Numeric (0-1). Node-degree quantile for hub detection. Default:
  `0.95`.

- chromatin_beds:

  Named list of BED file paths for orthogonal chromatin mark validation
  of eP/eG anchors. When non-empty, a *Chromatin Validation* sheet is
  added to the Excel workbook (see
  [`validate_epeG_by_chromatin`](https://zying106.github.io/looplook/reference/validate_epeG_by_chromatin.md)).
  Default: [`list()`](https://rdrr.io/r/base/list.html) (skip).

- write_output:

  Logical. If `TRUE` (default), write the refined Excel workbook to
  `out_dir`. If `FALSE`, return results without creating directories or
  files.

- quiet:

  Logical. If `TRUE`, suppress progress messages while preserving
  warnings. Default: `FALSE`.

## Value

An invisible named list:

- `loop_annotation` – Full refined 3D network with updated `loop_type`
  (e.g., eP-P) and two target gene columns:

  - `Active_Target_Genes`: Expression-filtered active-only targets (no
    fallback).

  - `Putative_Target_Genes`: Display column; may include linear
    nearest-gene fallback when `Active_Target_Genes` is empty.

  Refinement status columns: `Has_Active_Target`,
  `Retained_In_Functional_Network`, and `Refinement_Action`
  (`"retained_active_target"`, `"reclassified_silent_anchor"`,
  `"expression_filtered_no_active_target"`, or
  `"structural_only_no_active_target"`). All structural loops are
  preserved; filter on `Retained_In_Functional_Network` for the
  high-confidence active subset.

- `anchor_loci_annotation` – Filtered non-redundant anchor-locus
  annotations with expressed targets.

- `anchor_annotation` – Backward-compatible alias of
  `anchor_loci_annotation`.

- `promoter_centric_stats` – Gene-level connectivity statistics.

- `distal_element_stats` – Distal-element connectivity statistics.

- `target_annotation` – Target features (peaks) with gene assignments.
  Key columns include:

  - `All_Loop_Connected_Genes`: All genes from loop-connected anchors
    (P/G types).

  - `Regulated_promoter_genes`: Promoter genes supported by loop-anchor
    context.

  - `Assigned_Target_Genes`: Promoter-first 3D assignment (prioritises P
    \> G \> E).

  - `*_Filled` variants: Linear nearest-gene fallback when strict 3D
    assignments are empty.

  - `Regulated_promoter_Evidence`: Provenance of
    `Regulated_promoter_genes` (e.g., `local_promoter_overlap`,
    `direct_opposite_promoter`). **Read with**
    `Regulated_promoter_genes`; do not cross-reference with
    `Assigned_Target_Genes` or other columns.

  - `Regulated_promoter_Fallback_Evidence`: Provenance of
    `Regulated_promoter_genes_Filled`. **Read with**
    `Regulated_promoter_genes_Filled`; indicates which `*_Filled` column
    supplied the fallback gene.

- `target_gene_links` – Long-format peak-gene provenance table. Each row
  records one peak-gene linkage with full provenance. **Read**
  `evidence`, `anchor_role`, and `gene_role` **together as a group** –
  they jointly describe how each gene was assigned to each peak; do not
  interpret any one column in isolation.

  - `input_id`, `loop_ID`, `anchor_id`: Identifiers.

  - `gene`: Linked gene symbol.

  - `gene_role`: `"promoter"`, `"gene_body"`, or `"linear_annotation"`.

  - `source`: `"loop_anchor"` (3D-derived) or `"linear_annotation"`
    (nearest gene).

  - `evidence`: Provenance label – `"local_promoter_overlap"` (peak
    overlaps anchor promoter), `"direct_opposite_promoter"` (opposite
    anchor is promoter), `"gene_body_context"` (gene body linkage),
    `"expanded_promoter_loop"` (via ego-network expansion),
    `"linear_annotation"` (direct nearest gene), or `"linear_fallback"`
    (filled when 3D assignment was empty).

  - `anchor_role`: `"local_anchor"`, `"opposite_anchor"`,
    `"expanded_anchor"`, or `"linear_annotation"`.

  - `used_as_fallback`: Logical. `TRUE` when this link was added via the
    `*_Filled` linear nearest-gene fallback mechanism.

  - `in_regulated_promoter` through `in_assigned_target_filled`: Logical
    membership flags indicating which target annotation column(s) this
    gene appears in. A gene may appear in multiple columns
    simultaneously.

  - (Refine only) `Mean_Expression`: Per-gene mean expression value.

  - (Refine only) `Passes_Expression_Filter`: Logical. `TRUE` if
    `Mean_Expression >= threshold`.

- `chromatin_validation` – Data frame from
  [`validate_epeG_by_chromatin`](https://zying106.github.io/looplook/reference/validate_epeG_by_chromatin.md)
  with confidence levels and evidence strings for each eP/eG anchor.
  `NULL` when `chromatin_beds` is empty.

- `plots` – Named list of ggplot objects (dumbbell, rose, karyotype).

- `plot_list` – Backward-compatible alias of `plots`.

- `metadata` – Internal metadata list (parameters, versions, database
  versions). Not intended for direct use.

If `write_output = TRUE`, also writes `_Refined_Results.xlsx` to
`out_dir`. The workbook contains a *Chromatin Validation* sheet (when
`chromatin_beds` is provided) and a *Functional Loop Annotation* sheet
with only loops where `Retained_In_Functional_Network == TRUE`.

## Details

**Algorithmic Framework:**

- **Target Filtration:** Parses merged gene assignments (e.g.,
  `"GeneA;GeneB"`), evaluates individual genes against a defined
  expression threshold, and retains only transcriptionally active
  targets.

- **Biological Reclassification:** Reclassifies physically annotated
  promoters (`P`) and gene bodies (`G`) lacking active transcription as
  *expression-filtered silent* regulatory elements (`eP`, `eG`).
  **Important:** `eP`/`eG` labels indicate transcriptional silence at
  the reference gene – they do **not** constitute evidence of enhancer
  activity. Orthogonal chromatin data (ATAC-seq, H3K27ac, H3K4me1) are
  required for functional enhancer interpretation. The labels are
  retained for backward compatibility; interpret them as "inactive-P" /
  "inactive-G" rather than "enhancer-P" / "enhancer-G".

- **Expression-Aware Connectivity Statistics:** Recomputes
  promoter-centric and distal-element connectivity after
  expression-aware anchor refinement, while preserving all structural
  loops in the refined loop annotation. This separates the complete
  physical contact map from the high-confidence active subset.

- **External Target Refinement:** Filters auxiliary target mapping
  columns (e.g., `Assigned_Target_Genes_Filled`) based on expression
  criteria, ensuring that mapped 1D genomic features are exclusively
  linked to active genes.

- **Target Provenance Preservation:** Recomputes `*_Filled` membership
  flags in `target_gene_links` after expression filtering, retains only
  links still used by the refined target columns, and appends
  `Mean_Expression` plus `Passes_Expression_Filter`. Evidence labels
  such as `local_promoter_overlap`, `direct_opposite_promoter`, and
  `linear_fallback` are preserved.

**Design Philosophy:** This function does not discard structural loops
based on expression state. Hi-C, HiChIP, and PLAC-seq capture 3D
chromatin contacts; RNA-seq captures current transcriptional state. A
silent promoter may reflect cell-state, time-point, or technical factors
rather than absence of physical contact. All structural loops are
retained with refinement status columns, and a high-confidence
functional subset is provided via `Retained_In_Functional_Network` and
the *Functional Loop Annotation* Excel sheet.

**Interpretation of eP/eG labels:** `eP` and `eG` are
**expression-filtered silent states**, not functional enhancer
classifications. Bulk RNA-seq silence can arise from cell-type
proportions, time-point effects, sequencing depth, or promoter pausing –
none of which imply the locus has gained enhancer activity. These labels
should be read as "transcriptionally inactive P/G" and treated as
hypotheses requiring orthogonal validation (ATAC-seq, H3K27ac, H3K4me1,
or H3K27me3 depletion). Users with matched chromatin data should overlay
eP/eG loci against these tracks before interpreting them as putative
regulatory elements.

## Examples

``` r
# 1. Get paths to the required example files in the package
rdata_path <- system.file("extdata", "analysis_results.RData", package = "looplook")
expr_path <- system.file("extdata", "example_tpm.txt", package = "looplook")

# Safely load the pre-computed annotation result from RData
temp_env <- new.env()
load(rdata_path, envir = temp_env)
raw_annotation <- temp_env[[ls(temp_env)[1]]]
raw_annotation$loop_annotation <- head(raw_annotation$loop_annotation, 6)
raw_annotation$target_annotation <- head(raw_annotation$target_annotation, 3)
raw_annotation$promoter_centric_stats <- head(raw_annotation$promoter_centric_stats, 6)
raw_annotation$distal_element_stats <- head(raw_annotation$distal_element_stats, 6)

res_reclassified <- refine_loop_anchors_by_expression(
    annotation_res = raw_annotation,
    expr_matrix_file = expr_path,
    sample_columns = "con1",
    threshold = 1.0,
    unit_type = "TPM",
    species = "hg38",
    out_dir = tempdir(),
    project_name = "Example_Reclassified",
    reclassify_by_expression = TRUE,
    write_output = FALSE,
    quiet = TRUE
)
#> Warning: 100% of P/G anchors were reclassified to eP/eG. eP/eG labels indicate expression-aware enhancer-like states. Validate with orthogonal chromatin data (ATAC-seq, H3K27ac) before interpreting eP/eG anchors as functional enhancers.
print(table(res_reclassified$loop_annotation$loop_type))
#> 
#>   E-E  E-eP  G-eP  P-eG eP-eP 
#>     1     2     1     1     1 
```
