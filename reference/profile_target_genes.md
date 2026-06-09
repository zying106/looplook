# Integrative functional annotation and profiling of target genes

Integrates 3D genomic interaction data (e.g., Hi-C, HiChIP) with
transcriptomic profiles (RNA-seq) to evaluate the functional and
regulatory landscape of target genes. The workflow sequentially performs
differential expression profiling, gene set enrichment, transcription
factor motif scanning, Gene Ontology (GO) enrichment, and
protein-protein interaction (PPI) network construction.

## Usage

``` r
profile_target_genes(
  annotation_res,
  diff_file,
  lfc_col = "log2FoldChange",
  expr_matrix_file,
  metadata_file,
  target_source = c("loops", "targets"),
  target_mapping_mode = c("all", "promoter"),
  loop_types = c("E-P", "P-P"),
  include_Filled = TRUE,
  use_nearest_gene = FALSE,
  group_order = NULL,
  project_name = "Analysis",
  org_db = "org.Hs.eg.db",
  run_motif = FALSE,
  genome_id = "hg38",
  motif_p_thresh = 1e-04,
  motif_ntop = 5,
  run_go = FALSE,
  run_ppi = FALSE,
  ppi_score = 400,
  ppi_nSample = 400,
  heatmap_nSample = 99999,
  gsea_nSample = 99999,
  cnet_nSample = 50,
  stat_test = "wilcox.test",
  cor_method = "pearson",
  seed = NULL
)
```

## Arguments

- annotation_res:

  List. The result object returned by
  [`annotate_peaks_and_loops`](https://zying106.github.io/looplook/reference/annotate_peaks_and_loops.md).

- diff_file:

  Character. Path to the differential expression file (CSV/TSV).

- lfc_col:

  Character. The column name in `diff_file` representing Log2 Fold
  Change.

- expr_matrix_file:

  Character. Path to the normalized expression matrix.

- metadata_file:

  Character. Path to the sample metadata file.

- target_source:

  Character vector. Source of target genes to analyze.

- target_mapping_mode:

  Character. Specifies the mapping strategy for target genes.

- loop_types:

  Character vector. The specific loop types to analyze.

- include_Filled:

  Logical. If `TRUE`, utilizes the comprehensively merged gene
  assignment.

- use_nearest_gene:

  Logical. If `TRUE`, bypasses 3D loop-based gene assignment.

- group_order:

  Character vector. Optional factor levels to sort sample groups.

- project_name:

  Character. Prefix for all output files and plot titles.

- org_db:

  Character. Organism annotation database (e.g., "org.Hs.eg.db").

- run_motif:

  Logical. Whether to perform Transcription Factor Binding Site motif
  analysis.

- genome_id:

  Character. Reference genome assembly for motif sequence extraction.

- motif_p_thresh:

  Numeric. P-value threshold for scanning.

- motif_ntop:

  Numeric. Number of top enriched motifs to output.

- run_go:

  Logical. Whether to perform Gene Ontology (GO) enrichment.

- run_ppi:

  Logical. Whether to construct Protein-Protein Interaction networks.

- ppi_score:

  Numeric. Minimum combined confidence score for STRING edges.

- ppi_nSample:

  Numeric. Maximum number of genes to include in PPI.

- heatmap_nSample:

  Numeric. Maximum number of genes to plot in heatmap.

- gsea_nSample:

  Numeric. Maximum number of target genes to sample for GSEA.

- cnet_nSample:

  Numeric. Number of top GO terms to display in cnetplot.

- stat_test:

  Character. Statistical test for LFC comparisons.

- cor_method:

  Character. Method for sample correlation matrices.

- seed:

  Integer or NULL. Random seed for reproducible GSEA down-sampling and
  motif GC-matched background sampling. When `NULL` (default), the
  global RNG state is used; set to a positive integer for fully
  reproducible results. The seed is recorded in the result object as
  `attr(result, "seed")`.

## Value

An invisible nested list indexed by `target_source` (e.g., `"targets"`,
`"loops"`). Each element contains:

- `go_results`:

  Named list of data frames (one per gene set) containing GO enrichment
  results (if `run_go = TRUE`).

- `target_gene_sets`:

  Named list of character vectors containing target gene symbols.

- `plots`:

  Named list of ggplot objects (LFC_Violin, GSEA, Heatmap, Scatter,
  GO_Network, PPI_Network, etc.).

## Details

Two analysis steps use random sampling: GSEA target-gene down-sampling
(controlled by `gsea_nSample`, via unweighted sampling without
replacement to reduce enrichment bias) and motif background anchor
sampling (GC-matched, limited to 2 000 background regions per contrast).
GSEA tie-breaking for duplicate ranked values is deterministic
(position-based offset). For fully reproducible results, set the `seed`
parameter.

**Exploratory modules:** The GO enrichment (`run_go`), motif scanning
(`run_motif`), and PPI network (`run_ppi`) modules are *research-grade*
analyses that depend on external databases and algorithms. Results
should be treated as hypothesis-generating and validated with
independent experimental approaches. All three modules are disabled by
default.

## Note

If any downstream analysis module fails (e.g. due to missing optional
packages or network timeouts), the error propagates and stops the entire
function. To obtain partial results when some modules are unavailable,
disable problematic steps via `run_go = FALSE`, `run_ppi = FALSE`, or
`run_motif = FALSE`, or call the corresponding internal functions
individually.

## Examples

``` r
rdata_path <- system.file("extdata", "analysis_results.RData", package = "looplook")
diff_path <- system.file("extdata", "example_deg.txt", package = "looplook")
expr_path <- system.file("extdata", "example_tpm.txt", package = "looplook")
meta_path <- system.file("extdata", "example_coldata.txt", package = "looplook")
tmp <- new.env()
load(rdata_path, envir = tmp)
res <- tmp[[ls(tmp)[1]]]
profile_res <- profile_target_genes(
    annotation_res = res,
    diff_file = diff_path,
    expr_matrix_file = expr_path,
    metadata_file = meta_path,
    run_go = FALSE,
    run_ppi = FALSE,
    run_motif = FALSE,
    heatmap_nSample = 20,
    gsea_nSample = 20,
    cnet_nSample = 5
)
#> >>> Analysis Init | Root Project: Analysis
#> --- Reading files...
#> 
#> ================================================================
#> >>> Processing Source: [loops]
#> 
#> --- Task: EP_Genes (Valid Genes: 13) ---
#> 
#> Warning: qvalue::qvalue() failed, returning NA for qvalue. Error: missing values and NaN's not allowed if 'na.rm' is FALSE
#> Warning: ComplexHeatmap/circlize not installed; skipping heatmap.
#> 
#> --- Task: PP_Genes (Valid Genes: 67) ---
#> Warning: qvalue::qvalue() failed, returning NA for qvalue. Error: missing values and NaN's not allowed if 'na.rm' is FALSE
#> Warning: ComplexHeatmap/circlize not installed; skipping heatmap.
#> 
#> ================================================================
#> >>> Processing Source: [targets]
#> 
#> --- Task: Target_Genes (Valid Genes: 217) ---
#> Warning: qvalue::qvalue() failed, returning NA for qvalue. Error: missing values and NaN's not allowed if 'na.rm' is FALSE
#> Warning: ComplexHeatmap/circlize not installed; skipping heatmap.
#> 
#>  All analysis complete.
names(profile_res)
#> [1] "loops"   "targets"
```
