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
  genome_id = c("hg38", "hg19", "mm10", "mm9"),
  motif_p_thresh = 1e-04,
  motif_ntop = 5,
  motif_n_perm = 0L,
  run_go = FALSE,
  run_ppi = FALSE,
  ppi_score = 400,
  ppi_nSample = 400,
  ppi_species_id = NULL,
  heatmap_nSample = 99999,
  gsea_nSample = NULL,
  cnet_nSample = 50,
  universe_genes = NULL,
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
  Change. Default `"log2FoldChange"`.

- expr_matrix_file:

  Character. Path to the normalized expression matrix.

- metadata_file:

  Character. Path to the sample metadata file.

- target_source:

  Character vector. Source of target genes to analyze. Default
  `c("loops", "targets")`.

- target_mapping_mode:

  Character. Mapping strategy: `"all"` (all anchor-gene connections,
  `Putative_Target_Genes`) or `"promoter"` (promoter-only subset,
  `Promoter_Target_Genes`). Expression filtering is determined by the
  refinement stage of the supplied annotation object, not by this
  parameter. Default `"all"`.

- loop_types:

  Character vector. The specific loop types to analyze. Default
  `c("E-P", "P-P")`.

- include_Filled:

  Logical. If `TRUE`, utilizes the comprehensively merged gene
  assignment. Default `TRUE`.

- use_nearest_gene:

  Logical. If `TRUE`, bypasses 3D loop-based gene assignment. Default
  `FALSE`.

- group_order:

  Character vector. Optional factor levels to sort sample groups.
  Default `NULL`.

- project_name:

  Character. Prefix for all output files and plot titles. Default
  `"Analysis"`.

- org_db:

  Character. Organism annotation database (e.g., "org.Hs.eg.db").
  Default `"org.Hs.eg.db"`.

- run_motif:

  Logical. Whether to perform Transcription Factor Binding Site motif
  analysis. Default `FALSE`.

- genome_id:

  Character. Reference genome assembly for motif sequence extraction.
  Default `"hg38"`. One of `c("hg38", "hg19", "mm10", "mm9")`.

- motif_p_thresh:

  Numeric. P-value threshold for scanning. Default `1e-4`.

- motif_ntop:

  Numeric. Number of top enriched motifs to output. Default `5`.

- motif_n_perm:

  Integer. Number of within-component foreground/background label
  permutations for empirical motif P-value estimation. When positive,
  labels are shuffled within each loop component to partially account
  for anchor non-independence. `0` (default) retains Fisher exact test.
  Use `10-100` for testing. For inference, choose according to the
  number of motifs and desired P-value resolution; `1000` may be
  insufficient after multiple-testing correction across a full motif
  library. This is an exploratory calibration, not a fully
  component-level model.

- run_go:

  Logical. Whether to perform Gene Ontology (GO) enrichment. Default
  `FALSE`.

- run_ppi:

  Logical. Whether to construct Protein-Protein Interaction networks.
  Default `FALSE`.

- ppi_score:

  Numeric. Minimum combined confidence score for STRING edges. Default
  `400`.

- ppi_nSample:

  Numeric. Maximum number of genes to include in PPI. Default `400`.

- ppi_species_id:

  Integer or `NULL`. NCBI taxonomy ID for STRING PPI analysis (e.g.\\
  9606 for human, 10090 for mouse). When `NULL` (default), the ID is
  resolved from the `org_db` package name (only standard `org.*.eg.db`
  packages are recognised: human, mouse, rat, fruit fly, worm, yeast,
  zebrafish). An error is raised when the species cannot be resolved –
  no implicit guessing. Set explicitly for any other species. See
  <https://string-db.org> for a complete taxonomy list.

- heatmap_nSample:

  Numeric. Maximum number of genes to plot in heatmap. Default `99999`
  (effectively no limit). Reduce to `20-50` for readable heatmaps.

- gsea_nSample:

  Numeric or `NULL`. Maximum number of target genes to sample for GSEA.
  Default `NULL` (no down-sampling; all target genes are used).
  Down-sampling introduces Monte Carlo variance and should only be used
  when the target set is extremely large (\>1000 genes) for
  computational efficiency. For reproducible results, set the `seed`
  parameter when using down-sampling.

- cnet_nSample:

  Numeric. Number of top GO terms to display in cnetplot. Default `50`.

- universe_genes:

  Named numeric vector or `NULL`. GO background universe. Names are gene
  symbols, values are ranking metrics (e.g. LFC). When `NULL`, the full
  differential table is used as background. Default `NULL`.

- stat_test:

  Character. Statistical test for LFC comparisons. Default
  `"wilcox.test"`.

- cor_method:

  Character. Method for sample correlation matrices. Default
  `"pearson"`.

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

- `motif_results`:

  Named list of motif enrichment result tables, indexed by analysis
  task. Each task contains `proximal` and `distal` data frames when
  available (if `run_motif = TRUE`).

- `target_gene_sets`:

  Named list of character vectors containing target gene symbols.

- `plots`:

  Named list of ggplot objects (LFC_Violin, GSEA, Heatmap, Scatter,
  GO_Network, PPI_Network, etc.).

- `warnings`:

  Character vector of module-level warnings (e.g., "AnnotationDbi::GO
  failed: ..."). Empty if all modules succeeded.

## Details

GSEA operates on the full target gene set by default
(`gsea_nSample = NULL`). Down-sampling is available for very large sets
but introduces Monte Carlo variance. Motif background anchor sampling is
GC-matched, limited to 2 000 background regions per contrast. GSEA
tie-breaking for duplicate ranked values is deterministic
(position-based offset). For fully reproducible results, set the `seed`
parameter.

**Exploratory modules:** The GO enrichment (`run_go`), motif scanning
(`run_motif`), and PPI network (`run_ppi`) modules are *research-grade*
analyses that depend on external databases and algorithms. Results
should be treated as hypothesis-generating and validated with
independent experimental approaches. All three modules are disabled by
default. **Motif enrichment note:** By default, Fisher's exact test
treats each anchor as an independent observation. When
`motif_n_perm > 0`, empirical P-values are estimated by shuffling
foreground/background labels within loop components, which partially
accounts for component dependence. Pure-label components are
non-exchangeable; the procedure remains exploratory rather than a fully
component-level inferential model.

## Note

Downstream modules are run in a fail-soft mode. Module-level failures
(e.g. due to missing optional packages or network timeouts) are captured
as warnings and stored in the returned `warnings` element, allowing
other modules to complete. To proactively disable specific modules, set
`run_go = FALSE`, `run_ppi = FALSE`, or `run_motif = FALSE`.

## See also

[`annotate_peaks_and_loops`](https://zying106.github.io/looplook/reference/annotate_peaks_and_loops.md)
for initial 3D annotation,
[`refine_loop_anchors_by_expression`](https://zying106.github.io/looplook/reference/refine_loop_anchors_by_expression.md)
for expression-aware refinement.

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
#> Warning: GSEA: down-sampling 20 of 67 target genes. GSEA results represent a random subset, not the full gene set. Set gsea_nSample = NULL for full analysis. For fully reproducible results, set the `seed` parameter in profile_target_genes() (e.g., seed = 42).
#> Warning: qvalue::qvalue() failed, returning NA for qvalue. Error: missing values and NaN's not allowed if 'na.rm' is FALSE
#> Warning: ComplexHeatmap/circlize not installed; skipping heatmap.
#> 
#> ================================================================
#> >>> Processing Source: [targets]
#> 
#> --- Task: Target_Genes (Valid Genes: 217) ---
#> Warning: GSEA: down-sampling 20 of 217 target genes. GSEA results represent a random subset, not the full gene set. Set gsea_nSample = NULL for full analysis. For fully reproducible results, set the `seed` parameter in profile_target_genes() (e.g., seed = 42).
#> Warning: qvalue::qvalue() failed, returning NA for qvalue. Error: missing values and NaN's not allowed if 'na.rm' is FALSE
#> Warning: ComplexHeatmap/circlize not installed; skipping heatmap.
#> 
#>  All analysis complete.
names(profile_res)
#> [1] "loops"   "targets"
```
