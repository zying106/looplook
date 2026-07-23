# Expression-Aware refinement of loop anchors and target linkages

Integrates quantitative transcription data with 3D structural data.
Preferred input order: CAGE-seq / PRO-seq \> GRO-seq / NET-seq \> TT-seq
\> RNA-seq (methods that directly measure promoter activity or Pol II
elongation provide the most direct evidence for P/eP classification). to
reclassify regulatory elements and annotate functional status, deriving
a functionally interpretable regulatory network from physical chromatin
contacts. All structural loop rows are preserved; anchor type and gene
assignments are reclassified in place (silent promoters become eP/eG
with positional gene preserved; expression filtering is applied to
`Putative_Target_Genes` and `Promoter_Target_Genes`, not to anchor gene
assignments). Refinement status columns indicate which loops belong to
the high-confidence active subset.

## Usage

``` r
refine_loop_anchors_by_expression(
  annotation_res,
  expr_matrix_file = NULL,
  sample_columns = NULL,
  threshold = 1,
  threshold_mode = c("absolute", "quantile"),
  unit_type = "TPM",
  species = "hg38",
  out_dir = "./results/filtered",
  project_name = "HiChIP",
  color_palette = "Paired",
  karyo_bin_size = 1e+05,
  reclassify_by_expression = TRUE,
  hub_percentile = 0.95,
  hub_metric = c("unique_contacts", "total_loops"),
  chromatin_beds = list(),
  write_output = TRUE,
  quiet = FALSE,
  allow_rerefine = FALSE
)
```

## Arguments

- annotation_res:

  List. The raw foundational output object returned by
  [`annotate_peaks_and_loops`](https://zying106.github.io/looplook/reference/annotate_peaks_and_loops.md).

- expr_matrix_file:

  Path to a normalised expression matrix (genes x samples). Accepts
  steady-state RNA-seq (TPM/FPKM/RPKM), nascent transcription data
  (NET-seq, PRO-seq, GRO-seq, TT-seq), or CAGE-seq. Required for
  refinement. Default: `NULL`.

- sample_columns:

  Character vector or integer indices. Columns in `expr_matrix_file` to
  average. Default: `NULL`.

- threshold:

  Numeric. Minimum expression value (e.g., TPM \>= 1) for a gene to be
  considered active. Default: `1`. Set to `0` to retain any detected
  expression. Use `threshold_mode = "quantile"` to specify a quantile
  instead of an absolute value. For nascent transcription data (NET-seq,
  PRO-seq), use gene-body aggregated signal (e.g., RPKM) and adjust the
  threshold accordingly – Pol II elongation signal is typically sparser
  than steady-state mRNA, so a lower absolute threshold or quantile mode
  is recommended. Gene name matching is **case-insensitive**. Gene
  identifiers are canonicalised to uppercase before matching, so
  `"TP53"`, `"Tp53"`, and `"tp53"` are treated as the same gene. This
  accommodates the common human (all-uppercase) vs mouse (Title Case)
  convention mismatch.

- threshold_mode:

  Character. How to interpret the `threshold` value. `"absolute"`
  (default): `threshold` is a direct expression cutoff (e.g., TPM \>=
  1). `"quantile"`: `threshold` is a quantile of the expression
  distribution (e.g., `0.75` means the top 25\\ highly expressed genes
  are considered active). Quantile mode is dataset-adaptive and may be
  more robust across experiments with different sequencing depths. The
  effective expression cutoff is reported in the log.

- unit_type:

  Character. Expression unit for plot labels (e.g., `"TPM"`, `"FPKM"`,
  `"RPKM"`, `"NETseq_RPKM"`, `"raw_count"`). Default: `"TPM"`.

- species:

  Character. Genome assembly string. Default: `"hg38"`. For
  non-human/non-mouse species, pass `annotation_res` produced by
  [`annotate_peaks_and_loops`](https://zying106.github.io/looplook/reference/annotate_peaks_and_loops.md)
  with custom `txdb`/`org_db`.

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

  Integer. Bin width in base pairs (bp) for karyotype heatmaps. Default:
  `1e5` (100 kb). Typical range: 5e4-5e5 depending on genome size and
  resolution.

- reclassify_by_expression:

  Logical. If `TRUE` (default), silent promoters (P) and gene bodies (G)
  are reclassified as eP/eG while positional gene attribution is
  retained. Expression filtering is applied to current target-gene
  columns, not to anchor_gene. If `FALSE`, both anchor types and gene
  symbols are preserved unchanged – genes absent from the expression
  matrix are kept (their expression state is unmeasured, not silent).
  Measured-silent genes remain in the gene list but are tracked
  separately via the `Target_Expression_State` and
  `Passes_Expression_Filter` columns in the output.

- hub_percentile:

  Numeric (0-1). Node-degree quantile for hub detection. Default:
  `0.95`.

- hub_metric:

  Character. Which connectivity count to use for hub detection.
  `"unique_contacts"` (default): counts distinct neighbour anchor IDs.
  `"total_loops"`: counts all loop rows (backward compatible).

- chromatin_beds:

  Named list of BED file paths for orthogonal chromatin mark validation
  of eP/eG anchors. When non-empty, a *Chromatin Validation* sheet is
  added to the Excel workbook (see
  [`validate_epeG_by_chromatin`](https://zying106.github.io/looplook/reference/validate_epeG_by_chromatin.md)).
  Default: [`list()`](https://rdrr.io/r/base/list.html) (skip). Note: if
  you plan to use
  [`refine_loop_anchors_by_chromatin`](https://zying106.github.io/looplook/reference/refine_loop_anchors_by_chromatin.md),
  the chromatin validation and reclassification are handled there; you
  can skip this parameter to avoid redundant BED file reads.

- write_output:

  Logical. If `TRUE` (default), write the refined Excel workbook to
  `out_dir`. If `FALSE`, return results without creating directories or
  files.

- quiet:

  Logical. If `TRUE`, suppress progress messages while preserving
  warnings. Default: `FALSE`.

- allow_rerefine:

  Logical. If `FALSE` (default), stops when the input has already been
  expression-refined. Set to `TRUE` to override (for interactive
  exploration only).

## Value

A named list:

- `loop_annotation` – Full refined 3D network with updated `loop_type`
  (e.g., eP-P) and two target gene columns:

  - `Putative_Target_Genes`: **Primary target-assignment result.**
    Candidate target genes based on final anchor type, TSS assignment,
    loop topology, and the current refinement stage. In
    expression-refined objects, all listed genes pass the expression
    threshold. **Use for:** enhancer–gene linking benchmarks, CRISPRi
    validation, nearest-gene comparisons, loop connectivity analysis,
    full network export.

  - `Promoter_Target_Genes`: **Promoter-only subset of
    `Putative_Target_Genes`.** Contains only target genes supported by
    promoter-role (P/eP/dual) anchors. Always satisfies
    `Promoter_Target_Genes subset of Putative_Target_Genes`. **Use
    for:** promoter-centric GO, KEGG, GSEA, active regulatory network
    from promoter-anchored loops.

  - `Target_Expression_State`: Per-loop summary of target gene
    expression status (`"active"`, `"active_partial"`,
    `"measured_silent"`, `"unmeasured"`, `"mixed"`, `"no_target"`,
    `"not_assessed"`). **`"unmeasured"` != `"measured_silent"`:**
    unmeasured means the gene was not in the expression input and its
    activity is unknown; measured-silent means it was measured but fell
    below the activity threshold.

  Refinement status columns: `Has_Active_Target`,
  `Retained_In_Functional_Network`, and `Refinement_Action`
  (`"retained_active_target"`, `"reclassified_silent_anchor"`,
  `"expression_filtered_no_active_target"`, or
  `"structural_only_no_active_target"`). All structural loop rows are
  preserved (anchor types/genes are reclassified in place for silent
  elements); filter on `Retained_In_Functional_Network` for the
  high-confidence active subset.

- `anchor_loci_annotation` – Filtered non-redundant anchor-locus
  annotations with expressed targets.

- `anchor_annotation` – Backward-compatible alias of
  `anchor_loci_annotation`.

- `promoter_centric_stats` – Gene-level connectivity statistics.

- `distal_element_stats` – Distal anchor connectivity (E, dual, G, eG).

- `target_annotation` – Target features (peaks) with gene assignments.
  Key columns include:

  - `All_Loop_Connected_Genes`: Inclusive provenance union of all
    loop-anchor gene links. May include strict assignment-eligible
    targets and non-strict positional/enhancer candidates. Not a
    confirmed target-gene set.

  - `Regulated_promoter_genes`: Promoter genes supported by loop-anchor
    context.

  - `Assigned_Target_Genes`: Policy-based 3D assignment (default:
    promoter-first, then shorter path wins; see `target_priority`).

  - `*_Filled` variants: Linear nearest-gene fallback when strict 3D
    assignments are empty.

  - `Regulated_promoter_Evidence`: Provenance of
    `Regulated_promoter_genes` (e.g., `local_promoter_overlap`,
    `distal_promoter`). **Read with** `Regulated_promoter_genes`; do not
    cross-reference with `Assigned_Target_Genes` or other columns.

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

  - `gene_role`: `"promoter"`, `"gene_body"`, `"enhancer_candidate"`,
    `"positional_candidate"`, or `"linear_annotation"`.

  - `source`: `"loop_anchor"` (3D-derived) or `"linear_annotation"`
    (nearest gene).

  - `evidence`: Provenance label – `"local_promoter_overlap"` (peak
    overlaps anchor promoter), `"distal_promoter"` (promoter on the
    distal loop anchor), `"gene_body_context"` /
    `"distal_gene_body_context"` (gene body linkage),
    `"local_enhancer_candidate"` / `"distal_enhancer_candidate"` /
    `"expanded_enhancer_candidate"` (enhancer-associated linkage),
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

  - (Refine only) `Mean_Expression`: Per-gene mean expression value
    (`NA` when the gene is not in the expression matrix).

  - (Refine only) `Expression_State`: Character. Per-gene expression
    status: `"active"` (measured and above threshold),
    `"measured_silent"` (measured but below threshold), `"unmeasured"`
    (gene not present in expression input), `"measured_not_assessed"`
    (measured but no threshold available). Loop-level aggregation fields
    also include `"active_partial"`, `"mixed"`, `"no_target"`, and
    `"not_assessed"`.

  - (Refine only) `Passes_Expression_Filter`: Logical with three states.
    `TRUE` = measured and active; `FALSE` = measured but below
    threshold; `NA` = unmeasured or not assessed. **`NA`
    (unmeasured/not-assessed) \\\neq\\ `FALSE` (measured-silent).**

- `chromatin_validation` – Data frame from
  [`validate_epeG_by_chromatin`](https://zying106.github.io/looplook/reference/validate_epeG_by_chromatin.md)
  with enhancer evidence levels and evidence strings for each eP/eG
  anchor. `NULL` when `chromatin_beds` is not provided or is empty
  (default).

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
  *element-Promoter like* (`eP`) and *element-Genebody like* (`eG`)
  elements. **Important:** `eP`/`eG` labels denote structural genomic
  position (near a TSS or within a gene body) – they do **not**
  constitute evidence of enhancer activity, even though the associated
  gene is transcriptionally silent. Orthogonal chromatin data (ATAC-seq,
  H3K27ac, H3K4me1) are required for functional regulatory
  interpretation. The `e` prefix stands for **element** (a structurally
  defined genomic feature akin to a cis-regulatory element), not
  enhancer.

- **Expression-Aware Connectivity Statistics:** Recomputes
  promoter-centric and distal-element connectivity after
  expression-aware anchor refinement, while preserving all structural
  loop rows in the refined loop annotation (anchor types and genes are
  reclassified in place for silent elements). This separates the
  complete physical contact map from the high-confidence active subset.

- **External Target Refinement:** Filters auxiliary target mapping
  columns (e.g., `Assigned_Target_Genes_Filled`) based on expression
  criteria, ensuring that mapped 1D genomic features are exclusively
  linked to active genes.

- **Target Provenance Preservation:** Recomputes `*_Filled` membership
  flags in `target_gene_links` after expression filtering, retains all
  structural links with expression provenance, and appends
  `Mean_Expression`, `Expression_State`, and `Passes_Expression_Filter`.
  Evidence labels such as `local_promoter_overlap`, `distal_promoter`,
  and `linear_fallback` are preserved.

**Design Philosophy:** This function does not discard structural loop
rows based on expression state. Hi-C, HiChIP, and PLAC-seq capture 3D
chromatin contacts; RNA-seq captures current transcriptional state. A
silent promoter may reflect cell-state, time-point, or technical factors
rather than absence of physical contact. However, anchor gene
assignments and types are reclassified in place:
`anchor1_gene`/`anchor2_gene` may be set to `NA` and anchor types
downgraded (e.g. `P -> eP`, `G -> eG`) for silent regulatory elements.
All rows are retained; a high-confidence functional subset is provided
via `Retained_In_Functional_Network` and the *Functional Loop
Annotation* Excel sheet.

**Interpretation of eP/eG labels:** `eP` (**e**lement-**P**romoter like)
and `eG` (**e**lement-**G**enebody like) are **structural positional
classifications**, not functional enhancer classifications. An eP anchor
resides near a TSS and is structurally promoter-like; an eG anchor lies
within a gene body and is structurally genebody-like. The associated
gene may be transcriptionally silent in the current sample – this can
arise from cell-type proportions, time-point effects, sequencing depth,
or promoter pausing – but the anchor retains its positional identity as
a putative regulatory element. These labels should be treated as
hypotheses requiring orthogonal validation (ATAC-seq, H3K27ac, H3K4me1,
or H3K27me3 depletion). Users with matched chromatin data should overlay
eP/eG loci against these tracks before interpreting them as functionally
active regulatory elements.

**Pipeline guidance:** The recommended order is **expression refinement
first, then chromatin refinement.** When matched chromatin data are
available, follow expression refinement with
[`refine_loop_anchors_by_chromatin`](https://zying106.github.io/looplook/reference/refine_loop_anchors_by_chromatin.md)
– chromatin evidence resolves eP/eG into definitive types (E, P, G,
dual), while the expression-derived columns (`Putative_Target_Genes`,
`Promoter_Target_Genes`, `Retained_In_Functional_Network`) carry forward
the filtered target set for downstream profiling. The reverse order
(chromatin -\>expression) is unsupported and will stop. and will produce
incorrect results. See
[`refine_loop_anchors_by_chromatin`](https://zying106.github.io/looplook/reference/refine_loop_anchors_by_chromatin.md)
for the full pipeline recommendation.

## See also

[`annotate_peaks_and_loops`](https://zying106.github.io/looplook/reference/annotate_peaks_and_loops.md)
for initial 3D annotation,
[`refine_loop_anchors_by_chromatin`](https://zying106.github.io/looplook/reference/refine_loop_anchors_by_chromatin.md)
for chromatin-aware reclassification.

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
#> Warning: Only 27.3% (3 / 11) of annotation gene symbols match the expression matrix row names. Low ID mapping suggests different gene identifier conventions. Check that expression matrix row names and OrgDb annotations use the same convention (e.g., both SYMBOL, or both ENSEMBL). If the expression matrix uses Ensembl IDs or Entrez IDs while annotations use SYMBOL, convert identifiers first (e.g., via AnnotationDbi::mapIds or org.*.eg.db). 
print(table(res_reclassified$loop_annotation$loop_type))
#> 
#> E-E E-P G-P P-P 
#>   1   2   2   1 
```
