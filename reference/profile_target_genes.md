# End-to-end functional annotation pipeline integrating JASPAR transcription factor motif analysis, gene ontology enrichment, and protein-protein interaction network analysis

This function serves as a fully automated, end-to-end multi-omics
analysis pipeline. It elegantly bridges 3D genomic interactions (e.g.,
Hi-C, HiChIP) with transcriptomic data (RNA-seq) to systematically
characterize the functional landscape and regulatory mechanisms
underlying the identified target genes.

The suite integrates multiple analytical modules:

- **Differential Expression (LFC) Profiling:** Compares the Log2 Fold
  Change of target genes against the genomic background using tailored
  violin/boxplots.

- **GSEA Analysis:** Evaluates the rank-distribution of custom target
  gene sets within the global differential expression landscape.

- **Expression & 3D Connectivity Dynamics:** Generates Z-score
  expression heatmaps and correlates 3D structural connectivity (e.g.,
  Hub degree) with gene expression/LFC using scatter and sophisticated
  Raincloud plots.

- **TF Motif Analysis:** Scans JASPAR core motifs across proximal and
  distal loop anchors, generating motif enrichment rank scatter plots
  and sequence logos (SeqLogos).

- **Functional Enrichment (GO):** Performs Gene Ontology (BP)
  enrichment, visualizing results through divergent concept networks
  (cnetplot) and aggregated faceted lollipop charts.

- **Interaction Networks (PPI):** Constructs and visualizes
  high-confidence Protein-Protein Interaction networks via the STRING
  database.

## Usage

``` r
profile_target_genes(
  annotation_res,
  diff_file,
  lfc_col,
  expr_matrix_file,
  metadata_file,
  target_source = c("loops", "targets"),
  target_mapping_mode = c("all", "promoter"),
  loop_types = c("E-P", "P-P"),
  include_Filled = TRUE,
  use_nearest_gene = FALSE,
  group_order = NULL,
  project_name = "Analysis",
  out_dir = "./",
  org_db = "org.Hs.eg.db",
  run_motif = FALSE,
  genome_id = "hg19",
  motif_p_thresh = 1e-04,
  motif_ntop = 5,
  run_go = TRUE,
  run_ppi = FALSE,
  ppi_score = 400,
  ppi_nSample = 400,
  heatmap_nSample = 1000,
  gsea_nSample = 99999,
  cnet_nSample = 50,
  stat_test = "wilcox.test",
  cor_method = "pearson"
)
```

## Arguments

- annotation_res:

  List. The result object returned by
  [`annotate_peaks_and_loops`](https://zying106.github.io/looplook/reference/annotate_peaks_and_loops.md)
  or
  [`refine_loop_anchors_by_expression`](https://zying106.github.io/looplook/reference/refine_loop_anchors_by_expression.md).

- diff_file:

  Character. Path to the differential expression file (CSV/TSV). Must
  contain gene IDs (rownames) and a Log Fold Change column.

- lfc_col:

  Character. The column name in `diff_file` representing Log2 Fold
  Change (e.g., "log2FoldChange").

- expr_matrix_file:

  Character. Path to the normalized expression matrix (Genes x Samples,
  e.g., TPM or FPKM).

- metadata_file:

  Character. Path to the sample metadata file. Must contain "SampleID"
  and "Group" columns.

- target_source:

  Character vector. Source of target genes to analyze. Options include
  `"loops"` (genes connected by specific loops) and `"targets"` (genes
  overlapping with external target BED).

- target_mapping_mode:

  Character. Specifies the mapping strategy for target genes. Must be
  one of:

  - `"all"` (Default): Accepts all physically connected 3D targets
    (including Gene Bodies).

  - `"promoter"`: Strictly enforces that the 3D loop must explicitly
    anchor at a Promoter region.

- loop_types:

  Character vector. The specific loop types to analyze (e.g.,
  `c("E-P", "P-P")`). Active only when `"loops"` is included in
  `target_source`.

- include_Filled:

  Logical. If `TRUE`, utilizes the comprehensively merged gene
  assignment (retaining both loop-assigned targets and local host
  genes).

- use_nearest_gene:

  Logical. If `TRUE`, bypasses 3D loop-based gene assignment and
  strictly uses the nearest 1D gene (`SYMBOL`), serving as a spatial
  baseline reference.

- group_order:

  Character vector. Optional factor levels to sort sample groups in
  visualizations.

- project_name:

  Character. Prefix for all output files and plot titles.

- out_dir:

  Character. Directory path for saving output PDF plots and CSV tables.

- org_db:

  Character. Organism annotation database (e.g., "org.Hs.eg.db" for
  human, "org.Mm.eg.db" for mouse).

- run_motif:

  Logical. Whether to perform Transcription Factor Binding Site (TFBS)
  motif analysis on proximal and distal anchors.

- genome_id:

  Character. Reference genome assembly for motif sequence extraction
  (e.g., "hg38", "hg19", "mm10").

- motif_p_thresh:

  Numeric. P-value threshold for `motifmatchr` scanning.

- motif_ntop:

  Numeric. Number of top enriched motifs to output as SeqLogos.

- run_go:

  Logical. Whether to perform Gene Ontology (GO) ORA enrichment.

- run_ppi:

  Logical. Whether to construct Protein-Protein Interaction networks via
  the STRING database (requires internet access).

- ppi_score:

  Numeric. Minimum combined confidence score for STRING interaction
  edges (0-1000).

- ppi_nSample:

  Numeric. Maximum number of top variable/LFC genes to include in the
  PPI network to prevent visual overcrowding.

- heatmap_nSample:

  Numeric. Maximum number of highly variable genes to plot in the
  expression heatmap.

- gsea_nSample:

  Numeric. Maximum number of target genes to sample for GSEA to maintain
  computational efficiency.

- cnet_nSample:

  Numeric. Number of top GO terms to display in the divergent concept
  network plot.

- stat_test:

  Character. Statistical test for LFC comparisons (`"wilcox.test"` or
  `"t.test"`).

- cor_method:

  Character. Method for sample correlation matrices (`"pearson"` or
  `"spearman"`).

## Value

An invisible nested list containing:

- `go_results`: A list of top GO enrichment tables for each analyzed
  data source.

- `target_gene_sets`: The specific lists of gene IDs utilized for each
  analysis task.

## Examples

``` r
# =========================================================================
# 1. Get paths to real example files included in the package
# =========================================================================
rdata_path <- system.file("extdata", "analysis_results.RData", package = "looplook")
expr_path <- system.file("extdata", "example_tpm.txt", package = "looplook")
diff_path <- system.file("extdata", "example_deg.txt", package = "looplook")
meta_path <- system.file("extdata", "example_coldata.txt", package = "looplook")

# Check if all necessary files exist before running
if (rdata_path != "" && expr_path != "" && diff_path != "" && meta_path != "") {
  # Safely load the pre-computed annotation result from RData
  temp_env <- new.env()
  load(rdata_path, envir = temp_env)
  raw_annotation <- temp_env[[ls(temp_env)[1]]]

  out_base <- tempdir()

  # =========================================================================
  # Scenario A: Comprehensive Integrative Profiling (Loops + Targets)
  #
  # Application:
  # Best for complete HiChIP + ATAC-seq (or ChIP-seq) joint analysis.
  # It maps genes regulated by BOTH 3D distal loops and 1D direct target overlaps.
  # =========================================================================
  res_A <- profile_target_genes(
    annotation_res = raw_annotation,
    diff_file = diff_path,
    lfc_col = "log2FoldChange",
    expr_matrix_file = expr_path,
    metadata_file = meta_path,
    target_source = c("loops", "targets"),
    target_mapping_mode = "all",
    include_Filled = TRUE,
    use_nearest_gene = FALSE,
    project_name = "Scenario_A_Integrative",
    out_dir = out_base,
    run_go = FALSE, # Set to TRUE in real analysis
    run_ppi = FALSE
  )

  # =========================================================================
  # Scenario B: Peak-Centric & Promoter-Strict Profiling
  #
  # Application:
  # Focuses ONLY on 1D peaks ("targets") and STRICTLY restricts the mapping
  # to genes where peaks fall directly on their "promoter".
  # Ideal for baseline epigenetic profiling, ignoring distal 3D regulation.
  # =========================================================================
  res_B <- profile_target_genes(
    annotation_res = raw_annotation,
    diff_file = diff_path,
    lfc_col = "log2FoldChange",
    expr_matrix_file = expr_path,
    metadata_file = meta_path,
    target_source = "targets", # Only 1D targets
    target_mapping_mode = "promoter", # Strict promoter mapping
    include_Filled = TRUE,
    use_nearest_gene = FALSE,
    project_name = "Scenario_B_PromoterOnly",
    out_dir = out_base,
    run_go = FALSE,
    run_ppi = FALSE
  )

  # =========================================================================
  # Scenario C: Highly Conservative Distal Enhancer-Targeting
  #
  # Application:
  # Focuses on 3D loops. By setting `use_nearest_gene = TRUE` and
  # `include_Filled = FALSE`, it forces the selection of strictly the nearest
  # genes to the anchors, preventing the expansion of gene lists from broader
  # annotated regions. Ensures a highly conservative gene target list.
  # =========================================================================
  res_C <- profile_target_genes(
    annotation_res = raw_annotation,
    diff_file = diff_path,
    lfc_col = "log2FoldChange",
    expr_matrix_file = expr_path,
    metadata_file = meta_path,
    target_source = "loops", # Only 3D loops
    target_mapping_mode = "all",
    include_Filled = FALSE,
    use_nearest_gene = TRUE, # Strictly nearest gene
    project_name = "Scenario_C_StrictLoops",
    out_dir = out_base,
    run_go = FALSE,
    run_ppi = FALSE
  )
}
#> >>> Analysis Init | Root Project: Scenario_A_Integrative
#> --- Reading files...
#> Loaded Promoter Centric Stats.
#> 
#> ================================================================
#> >>> Processing Source: [loops]
#>      Loops: Successfully located 'Putative_Target_Genes'
#> 
#> --- Task: EP_Genes (Valid Genes: 36) ---
#>     Saved LFC Violin Plot: /tmp/Rtmp2HPPhq/Scenario_A_Integrative_loops_EP_Genes_LFC_Violin.pdf
#>     Saved GSEA Plot (Custom): /tmp/Rtmp2HPPhq/Scenario_A_Integrative_loops_EP_Genes_GSEA.pdf
#>     [Extra Analysis] Found 'n_Linked_Distal'. Running distal connectivity analysis...
#> 
#> --- Task: PP_Genes (Valid Genes: 124) ---
#>     Saved LFC Violin Plot: /tmp/Rtmp2HPPhq/Scenario_A_Integrative_loops_PP_Genes_LFC_Violin.pdf
#>     Saved GSEA Plot (Custom): /tmp/Rtmp2HPPhq/Scenario_A_Integrative_loops_PP_Genes_GSEA.pdf
#>     [Extra Analysis] Found 'n_Linked_Distal'. Running distal connectivity analysis...
#> 
#> ================================================================
#> >>> Processing Source: [targets]
#>       Targets: Successfully located 'Assigned_Target_Genes_Filled'
#> 
#> --- Task: Target_Genes (Valid Genes: 937) ---
#>     Saved LFC Violin Plot: /tmp/Rtmp2HPPhq/Scenario_A_Integrative_targets_Target_Genes_LFC_Violin.pdf
#>     Saved GSEA Plot (Custom): /tmp/Rtmp2HPPhq/Scenario_A_Integrative_targets_Target_Genes_GSEA.pdf
#>     [Extra Analysis] Found 'n_Linked_Distal'. Running distal connectivity analysis...
#> 
#>  All analysis complete.
#> >>> Analysis Init | Root Project: Scenario_B_PromoterOnly_Promoter
#> --- Reading files...
#> Loaded Promoter Centric Stats.
#> 
#> ================================================================
#> >>> Processing Source: [targets]
#>       Targets: Successfully located 'Regulated_promoter_genes_Filled'
#> 
#> --- Task: Target_Genes (Valid Genes: 931) ---
#>     Saved LFC Violin Plot: /tmp/Rtmp2HPPhq/Scenario_B_PromoterOnly_Promoter_targets_Target_Genes_LFC_Violin.pdf
#>     Saved GSEA Plot (Custom): /tmp/Rtmp2HPPhq/Scenario_B_PromoterOnly_Promoter_targets_Target_Genes_GSEA.pdf
#>     [Extra Analysis] Found 'n_Linked_Distal'. Running distal connectivity analysis...
#> 
#>  All analysis complete.
#> >>> Analysis Init | Root Project: Scenario_C_StrictLoops
#> --- Reading files...
#> Loaded Promoter Centric Stats.
#> 
#> ================================================================
#> >>> Processing Source: [loops]
#>      Loops: Successfully located 'Putative_Target_Genes'
#> 
#> --- Task: EP_Genes (Valid Genes: 36) ---
#>     Saved LFC Violin Plot: /tmp/Rtmp2HPPhq/Scenario_C_StrictLoops_loops_EP_Genes_LFC_Violin.pdf
#>     Saved GSEA Plot (Custom): /tmp/Rtmp2HPPhq/Scenario_C_StrictLoops_loops_EP_Genes_GSEA.pdf
#>     [Extra Analysis] Found 'n_Linked_Distal'. Running distal connectivity analysis...
#> 
#> --- Task: PP_Genes (Valid Genes: 124) ---
#>     Saved LFC Violin Plot: /tmp/Rtmp2HPPhq/Scenario_C_StrictLoops_loops_PP_Genes_LFC_Violin.pdf
#>     Saved GSEA Plot (Custom): /tmp/Rtmp2HPPhq/Scenario_C_StrictLoops_loops_PP_Genes_GSEA.pdf
#>     [Extra Analysis] Found 'n_Linked_Distal'. Running distal connectivity analysis...
#> 
#>  All analysis complete.
```
