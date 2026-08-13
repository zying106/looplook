# Render a Publication-Ready 3D Annotation Report

One-click parameterised R Markdown report that executes the full
looplook pipeline (annotation -\> refinement -\> profiling) and renders
an interpretation-ready HTML document suitable for sharing with
collaborators.

## Usage

``` r
looplook_report(
  bedpe_file = NULL,
  target_bed = NULL,
  expr_matrix_file = NULL,
  sample_columns = NULL,
  species = "hg38",
  project_name = "looplook Analysis",
  out_dir = "looplook_results",
  threshold = 1,
  unit_type = "TPM",
  reclassify_by_expression = TRUE,
  tss_region = c(-2000, 2000),
  neighbor_hop = 0,
  hub_percentile = 0.95,
  color_palette = "Set2",
  target_mapping_mode = "all",
  include_Filled = TRUE,
  use_nearest_gene = FALSE,
  target_source = "targets",
  loop_types = c("E-P", "P-P"),
  stat_test = "wilcox.test",
  run_go = FALSE,
  run_ppi = FALSE,
  run_motif = FALSE,
  genome_id = species,
  motif_p_thresh = 1e-04,
  motif_ntop = 5,
  motif_n_perm = 0L,
  ppi_score = 400,
  ppi_nSample = 400,
  heatmap_nSample = 99999,
  gsea_nSample = 99999,
  cnet_nSample = 50,
  karyo_bin_size = 1e+05,
  diff_file = NULL,
  lfc_col = "log2FoldChange",
  metadata_file = NULL,
  precomputed_res = NULL,
  chromatin_beds = list(),
  chromatin_bw = NULL,
  bw_ratio_threshold = 3,
  enhancer_bed = NULL,
  output_file = NULL,
  quiet = FALSE,
  seed = NULL,
  universe_genes = NULL,
  ...
)
```

## Arguments

- bedpe_file:

  Path to a BEDPE file of chromatin loops. Default `NULL` (requires
  `precomputed_res` instead).

- target_bed:

  Optional path to a BED file of genomic features.

- expr_matrix_file:

  Optional path to a normalised expression matrix.

- sample_columns:

  Sample columns in the expression matrix to average. Default `NULL`
  (all columns).

- species:

  Genome assembly (`"hg38"`, `"hg19"`, `"mm10"`, `"mm9"`).

- project_name:

  Character. Project prefix for the report title.

- out_dir:

  Output directory. Created if missing.

- threshold:

  Numeric. Expression threshold for active gene classification. Default
  `1.0`.

- unit_type:

  Character. Expression unit label for plot annotations. Default
  `"TPM"`.

- reclassify_by_expression:

  Logical. Reclassify silent promoters as eP/eG. Default `TRUE`.

- tss_region:

  Numeric vector of length 2. TSS flanking region in bp. Default
  `c(-2000, 2000)`.

- neighbor_hop:

  Integer. k-hop ego-network expansion order for loop connectivity
  analysis. `0` (default) restricts to direct loop contacts. `1`
  additionally includes 2-hop expanded targets. Values greater than `1`
  are not supported.

- hub_percentile:

  Numeric. Quantile threshold for hub classification. Default `0.95`.

- color_palette:

  Character. RColorBrewer qualitative palette name. Default `"Set2"`.

- target_mapping_mode:

  Character. Mapping strategy for target genes. Default `"all"`.

- include_Filled:

  Logical. Include comprehensively merged gene assignments. Default
  `TRUE`.

- use_nearest_gene:

  Logical. Bypass 3D loop-based assignment and use linear proximity.
  Default `FALSE`.

- target_source:

  Character vector. Source of target genes to profile. Default
  `"targets"`.

- loop_types:

  Character vector. Loop types to include in profiling. Default
  `c("E-P", "P-P")`.

- stat_test:

  Character. Statistical test for violin comparisons. Default
  `"wilcox.test"`.

- run_go:

  Logical. Run GO enrichment (requires clusterProfiler). Default
  `FALSE`.

- run_ppi:

  Logical. Run protein-protein interaction network analysis. Default
  `FALSE`.

- run_motif:

  Logical. Run transcription factor motif analysis. Default `FALSE`.

- genome_id:

  Character. Reference genome for motif scanning. Defaults to `species`.

- motif_p_thresh:

  Numeric. P-value threshold for motif enrichment. Default `1e-4`.

- motif_ntop:

  Numeric. Number of top motifs to display. Default `5`.

- motif_n_perm:

  Integer. Number of within-component label permutations for empirical
  motif P-value estimation. Passed to
  [`profile_target_genes`](https://zying106.github.io/looplook/reference/profile_target_genes.md).
  Default `0` (Fisher exact test).

- ppi_score:

  Numeric. Minimum STRING combined score. Default `400`.

- ppi_nSample:

  Numeric. Maximum genes to include in PPI. Default `400`.

- heatmap_nSample:

  Numeric. Maximum genes in expression heatmap. Default `99999`.

- gsea_nSample:

  Numeric. Maximum genes sampled for GSEA. Default `99999`.

- cnet_nSample:

  Numeric. Number of GO terms in cnetplot. Default `50`.

- karyo_bin_size:

  Numeric. Bin size for karyotype heatmaps. Default `100000`.

- diff_file:

  Optional differential expression result file.

- lfc_col:

  Column name for log2 fold change in `diff_file`.

- metadata_file:

  Optional sample metadata file.

- precomputed_res:

  Optional. Either a `.RData` file path or an in-memory list object
  returned by `annotate_peaks_and_loops`. When provided, annotation is
  skipped and refinement starts from this object.

- chromatin_beds:

  Named list of BED file paths for orthogonal chromatin mark validation
  (passed to
  [`refine_loop_anchors_by_expression`](https://zying106.github.io/looplook/reference/refine_loop_anchors_by_expression.md)).
  When non-empty, a *Chromatin Validation* section appears in the report
  with confidence-level distribution for eP/eG anchors. Default:
  [`list()`](https://rdrr.io/r/base/list.html) (skip).

- chromatin_bw:

  Named list of bigWig file paths, or `NULL`. Passed to
  [`refine_loop_anchors_by_chromatin`](https://zying106.github.io/looplook/reference/refine_loop_anchors_by_chromatin.md).
  Requires `"H3K4me1"` and `"H3K4me3"`. **Strongly recommended** for
  resolving true dual-signature elements. Default: `NULL`.

- bw_ratio_threshold:

  Numeric. Minimum H3K4me1/H3K4me3 ratio for dual classification.
  Default: `3`.

- enhancer_bed:

  Character or `NULL`. Path to a curated enhancer BED file (e.g.
  FANTOM5, ENCODE cCREs). Default: `NULL`.

- output_file:

  Character. Output HTML file name. `NULL` derives from `project_name`.

- quiet:

  Logical. Suppress rendering output. Default `FALSE`.

- seed:

  Integer or NULL. Passed to
  [`profile_target_genes`](https://zying106.github.io/looplook/reference/profile_target_genes.md)
  for reproducible GSEA and motif sampling. Default `NULL`.

- universe_genes:

  Named numeric vector or `NULL`. GO background universe. Passed to
  [`profile_target_genes`](https://zying106.github.io/looplook/reference/profile_target_genes.md).
  Default `NULL`.

- ...:

  Additional arguments passed to
  [`rmarkdown::render`](https://pkgs.rstudio.com/rmarkdown/reference/render.html).

## Value

The path to the generated HTML report (invisibly).

## Details

The profiling stage uses R's random number generator for GSEA gene-set
down-sampling (`gsea_nSample`) and motif background anchor sampling.
Call [`set.seed()`](https://rdrr.io/r/base/Random.html) before
`looplook_report()` for fully reproducible results.

## Examples

``` r
if (requireNamespace("rmarkdown", quietly = TRUE) &&
  requireNamespace("knitr", quietly = TRUE)) {
  temp_env <- new.env()
  load(system.file("extdata", "analysis_results.RData", package = "looplook"), envir = temp_env)
  precomputed_res <- temp_env[[ls(temp_env)[1]]]
  precomputed_res$loop_annotation <- head(precomputed_res$loop_annotation, 6)
  precomputed_res$target_annotation <- head(precomputed_res$target_annotation, 3)
  precomputed_res$promoter_centric_stats <- head(precomputed_res$promoter_centric_stats, 6)
  precomputed_res$distal_element_stats <- head(precomputed_res$distal_element_stats, 6)

  report_path <- looplook_report(
    precomputed_res = precomputed_res,
    project_name = "Example",
    out_dir = tempdir(),
    output_file = "looplook-example-report.html",
    quiet = TRUE,
    run_go = FALSE,
    run_ppi = FALSE,
    run_motif = FALSE
  )
  file.exists(report_path)
}
#> [1] TRUE
```
