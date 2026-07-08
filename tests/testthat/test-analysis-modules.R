# tests/testthat/test-analysis-modules.R — unit tests for analysis pipeline modules
# Profiling integration tests (GO/GSEA/PPI/motif) are computationally heavy;
# skip on Bioconductor to stay within the build time budget.
skip_on_bioc()

# Shared setup: load pre-computed annotation and expression data
rdata_path <- system.file("extdata", "analysis_results.RData", package = "looplook")
has_data <- file.exists(rdata_path) && rdata_path != ""
if (has_data) {
  temp_env <- new.env()
  load(rdata_path, envir = temp_env)
  raw_annotation <- temp_env[[ls(temp_env)[1]]]

  gene_sets <- looplook:::extract_target_gene_sets(raw_annotation,
    src = "targets",
    include_Filled = TRUE
  )
  target_genes <- head(intersect(
    gene_sets[["Target_Genes"]],
    rownames(looplook:::read_robust_general(
      system.file("extdata", "example_tpm.txt", package = "looplook"),
      header = TRUE, row_name = 1, min_cols = 1
    ))
  ), 12)

  diff_df <- looplook:::read_robust_general(
    system.file("extdata", "example_deg.txt", package = "looplook"),
    header = TRUE, row_name = 1, min_cols = 1
  )
  global_glist <- sort(setNames(
    diff_df[["log2FoldChange"]],
    rownames(diff_df)
  ), decreasing = TRUE)
}

test_that("run_lfc_violin returns ggplot with valid input", {
  skip_if_not(has_data, "Pre-computed RData not available")
  skip_if_not_installed("ggplot2")

  p <- looplook:::run_lfc_violin(target_genes, global_glist, "wilcox.test", "Test_Violin")
  expect_s3_class(p, "ggplot")
  # fewer than 3 genes returns NULL
  p_null <- looplook:::run_lfc_violin(
    head(target_genes, 2), global_glist,
    "wilcox.test", "TooFew"
  )
  expect_null(p_null)
})

test_that("run_gsea_analysis returns list with result and plot", {
  skip_if_not(has_data, "Pre-computed RData not available")
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("enrichplot")

  # Ensure .Random.seed is initialized for clusterProfiler::GSEA
  set.seed(1)
  # Use small gene set for speed
  small_targets <- head(target_genes, 10)
  out <- looplook:::run_gsea_analysis(
    small_targets, global_glist, 20,
    "Test_GSEA"
  )
  expect_type(out, "list")
  expect_true("result" %in% names(out))
  expect_true("plot" %in% names(out))
  # fewer than 2 overlap should return NULL
  out_null <- looplook:::run_gsea_analysis(
    c("FAKE1", "FAKE2"), global_glist, 50,
    "Empty"
  )
  expect_null(out_null$result)
})

test_that("run_go_enrichment returns list with result and plot", {
  skip_if_not(has_data, "Pre-computed RData not available")
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("org.Hs.eg.db")

  # Use small gene set for speed
  small_targets <- head(target_genes, 8)
  out <- looplook:::run_go_enrichment(small_targets, "org.Hs.eg.db", global_glist,
    cnet_nSample = 5, project_name = "Test_GO"
  )
  expect_type(out, "list")
  expect_true("result" %in% names(out))
  # result should be a data.frame with expected columns
  if (!is.null(out$result) && nrow(out$result) > 0) {
    expect_s3_class(out$result, "data.frame")
  }
})

test_that("run_go_enrichment with NULL universe_genes", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("org.Hs.eg.db")

  out <- looplook:::run_go_enrichment(
    c("TP53", "BRCA1", "MYC"), "org.Hs.eg.db", NULL,
    cnet_nSample = 3, project_name = "Test_null_univ"
  )
  expect_type(out, "list")
})

test_that("run_go_enrichment handles unmapped genes gracefully", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("org.Hs.eg.db")

  # Mix of real and fake gene names
  out <- looplook:::run_go_enrichment(
    c("TP53", "FAKEGENE123", "BRCA1"), "org.Hs.eg.db",
    setNames(c(1, 2, 3), c("TP53", "FAKEGENE123", "BRCA1")),
    cnet_nSample = 3, project_name = "Test_unmapped"
  )
  expect_type(out, "list")
})

test_that("run_ppi_analysis returns ggplot for valid input", {
  skip_if_not(has_data, "Pre-computed RData not available")
  skip_if_not_installed("STRINGdb")
  skip_if_not_installed("ggraph")
  run_network_tests <- identical(
    tolower(Sys.getenv("LOOPLOOK_RUN_NETWORK_TESTS", unset = "false")),
    "true"
  )
  skip_if_not(
    run_network_tests,
    "External STRINGdb integration test disabled; set LOOPLOOK_RUN_NETWORK_TESTS=true to run."
  )
  # STRINGdb requires external network access even when explicitly enabled.
  has_net <- tryCatch(
    {
      con <- url("https://string-db.org", open = "r")
      close(con)
      TRUE
    },
    error = function(e) FALSE
  )
  skip_if_not(has_net, "STRINGdb not reachable")

  p <- looplook:::run_ppi_analysis(target_genes, global_glist, "org.Hs.eg.db",
    ppi_score = 700, ppi_ntop = 20, "Test_PPI"
  )
  # may return NULL if no interactions found — not a failure
  if (!is.null(p)) {
    expect_s3_class(p, "gg")
  }
})

test_that("run_heatmap_and_connectivity returns plot list", {
  skip_if_not(has_data, "Pre-computed RData not available")
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")

  tpm_mat <- looplook:::read_robust_general(
    system.file("extdata", "example_tpm.txt", package = "looplook"),
    header = TRUE, row_name = 1, min_cols = 1
  )
  meta_raw <- data.frame(
    SampleID = colnames(tpm_mat),
    Group = rep(c("A", "B"), length.out = ncol(tpm_mat)),
    stringsAsFactors = FALSE
  )

  plots <- looplook:::run_heatmap_and_connectivity(target_genes, tpm_mat, meta_raw,
    loop_stats_df = NULL, global_glist,
    heatmap_ntop = 50, cor_method = "pearson",
    current_proj_name = "Test_Heat", source_type = "targets"
  )
  expect_type(plots, "list")
})

test_that("run_heatmap_and_connectivity with loop_stats_df and scatter", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggpointdensity")

  tpm_mat <- data.frame(
    S1 = c(10, 20, 30, 40, 50, 60, 70, 80),
    S2 = c(15, 25, 35, 45, 55, 65, 75, 85),
    row.names = paste0("Gene", 1:8)
  )
  meta_raw <- data.frame(
    SampleID = c("S1", "S2"),
    Group = c("A", "B"),
    stringsAsFactors = FALSE
  )
  loop_stats_df <- data.frame(
    Gene = paste0("Gene", 1:8),
    Total_Loops = c(10, 8, 6, 4, 3, 2, 1, 1),
    Is_High_Connectivity_Gene = c("Yes", "Yes", "No", "No", "No", "No", "No", "No"),
    Is_High_Distal_Connectivity_Gene = c("Yes", "No", "No", "No", "No", "No", "No", "No"),
    stringsAsFactors = FALSE
  )
  global_glist <- setNames(c(2, 1, 0.5, -0.5, -1, -2, 0, 0.3), paste0("Gene", 1:8))

  plots <- looplook:::run_heatmap_and_connectivity(
    paste0("Gene", 1:8), tpm_mat, meta_raw,
    loop_stats_df, global_glist,
    heatmap_ntop = 50, cor_method = "pearson",
    current_proj_name = "Test", source_type = "targets"
  )
  expect_type(plots, "list")
  expect_true("Scatter" %in% names(plots))
})

test_that("run_heatmap_and_connectivity with n_Linked_Distal target_col", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggpointdensity")

  tpm_mat <- data.frame(
    S1 = c(10, 20, 30, 40, 50, 60),
    S2 = c(15, 25, 35, 45, 55, 65),
    row.names = paste0("Gene", 1:6)
  )
  meta_raw <- data.frame(
    SampleID = c("S1", "S2"),
    Group = c("A", "B"),
    stringsAsFactors = FALSE
  )
  loop_stats_df <- data.frame(
    Gene = paste0("Gene", 1:6),
    Total_Loops = c(10, 8, 6, 4, 3, 2),
    n_Linked_Distal = c(5, 4, 3, 2, 1, 0),
    stringsAsFactors = FALSE
  )
  global_glist <- setNames(c(2, 1, 0.5, -0.5, -1, -2), paste0("Gene", 1:6))

  plots <- looplook:::run_heatmap_and_connectivity(
    paste0("Gene", 1:6), tpm_mat, meta_raw,
    loop_stats_df, global_glist,
    heatmap_ntop = 50, cor_method = "pearson",
    current_proj_name = "Test", source_type = "targets",
    target_col = "n_Linked_Distal", skip_heatmap = TRUE
  )
  expect_type(plots, "list")
  expect_true("Scatter" %in% names(plots))
})

test_that("run_heatmap_and_connectivity returns empty for no valid samples", {
  tpm_mat <- data.frame(S1 = c(10, 20), row.names = c("A", "B"))
  meta_raw <- data.frame(SampleID = c("X", "Y"), Group = c("A", "B"), stringsAsFactors = FALSE)
  plots <- looplook:::run_heatmap_and_connectivity(
    c("A", "B"), tpm_mat, meta_raw, NULL, setNames(c(1, 2), c("A", "B")),
    50, "pearson", "Test", "targets"
  )
  expect_equal(length(plots), 0)
})

test_that("run_heatmap_and_connectivity returns empty for invalid target_col", {
  tpm_mat <- data.frame(S1 = c(10, 20), row.names = c("A", "B"))
  meta_raw <- data.frame(SampleID = c("S1"), Group = c("A"), stringsAsFactors = FALSE)
  loop_stats_df <- data.frame(Gene = c("A", "B"), Total_Loops = c(1, 2), stringsAsFactors = FALSE)
  plots <- looplook:::run_heatmap_and_connectivity(
    c("A", "B"), tpm_mat, meta_raw, loop_stats_df, setNames(c(1, 2), c("A", "B")),
    50, "pearson", "Test", "targets",
    target_col = "NonExistent"
  )
  expect_equal(length(plots), 0)
})

# --- run_gsea_analysis: deterministic tie-breaking ---
test_that("run_gsea_analysis uses deterministic tie-breaking", {
  gl <- setNames(c(2.5, 2.5, 2.0, 1.0, 0.5), c("A", "B", "C", "D", "E"))
  tg <- c("A", "B")
  # Run twice with same seed; results must be deterministic (no RNG for tie-breaking).
  # GSEA internal permutation may accumulate ~1e-7 floating-point differences.
  old_seed <- tryCatch(.Random.seed, error = function(e) NULL)
  on.exit(
    {
      if (!is.null(old_seed)) assign(".Random.seed", old_seed, envir = globalenv())
    },
    add = TRUE
  )
  set.seed(42)
  r1 <- looplook:::run_gsea_analysis(tg, gl, NULL, "Test")
  set.seed(42)
  r2 <- looplook:::run_gsea_analysis(tg, gl, NULL, "Test")
  expect_equal(r1$result$NES, r2$result$NES, tolerance = 1e-6)
  expect_equal(r1$result$pvalue, r2$result$pvalue, tolerance = 1e-6)
})

# --- plot_summary_go_lollipop: handles valid and edge inputs ---
test_that("plot_summary_go_lollipop returns ggplot objects", {
  # Mock GO results mimicking clusterProfiler output
  go_df <- data.frame(
    ONTOLOGY = rep(c("BP", "MF"), each = 3),
    Description = c(
      "cell cycle", "DNA repair", "apoptosis",
      "kinase activity", "DNA binding", "ATP binding"
    ),
    pvalue = c(0.001, 0.01, 0.02, 0.005, 0.03, 0.04),
    Count = c(30, 20, 15, 25, 18, 12),
    geneID = c("A/B", "C/D", "E/F", "G/H", "I/J", "K/L"),
    CleanLoopType = "test",
    LoopType = "test",
    Source = "targets",
    stringsAsFactors = FALSE
  )
  plots <- looplook:::plot_summary_go_lollipop(list(go_df), "Test")
  expect_type(plots, "list")
  expect_gt(length(plots), 0)
  # Each element should be a ggplot
  for (p in plots) expect_s3_class(p, "ggplot")
})

test_that("plot_summary_go_lollipop handles empty input", {
  expect_equal(length(looplook:::plot_summary_go_lollipop(list(), "Empty")), 0L)
  null_df <- data.frame()
  expect_equal(length(looplook:::plot_summary_go_lollipop(list(null_df), "Null")), 0L)
})

test_that("plot_summary_go_lollipop with ggtext styling", {
  skip_if_not_installed("ggtext")
  go_df <- data.frame(
    ONTOLOGY = c("BP", "MF", "CC"),
    Description = c("cell cycle", "kinase activity", "nucleus"),
    pvalue = c(0.001, 0.01, 0.02),
    Count = c(20, 15, 10),
    geneID = c("A/B", "C/D", "E/F"),
    CleanLoopType = "test",
    LoopType = "test",
    Source = "targets",
    stringsAsFactors = FALSE
  )
  plots <- looplook:::plot_summary_go_lollipop(list(go_df), "Test_ggtext")
  expect_gt(length(plots), 0)
})

# --- .sample_gc_matched_background: edge cases ---
test_that(".sample_gc_matched_background guards: empty foreground", {
  empty_gr <- GenomicRanges::GRanges()
  bg_gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1:10, 11:20))
  # No foreground to match → returns all background as-is
  result <- looplook:::.sample_gc_matched_background(empty_gr, bg_gr, NULL,
    max_bg = 5L, gc_bins = 5L
  )
  expect_equal(length(result), 10L)
})

test_that(".sample_gc_matched_background guards: empty background", {
  fg_gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1:10, 11:20))
  empty_bg <- GenomicRanges::GRanges()
  result <- looplook:::.sample_gc_matched_background(fg_gr, empty_bg, NULL,
    max_bg = 5L, gc_bins = 5L
  )
  expect_equal(length(result), 0L)
})

# ════════════════════════════════════════════════════════════════════════════
# Additional coverage tests
# ════════════════════════════════════════════════════════════════════════════

test_that("run_gsea_analysis with gsea_nSample limiting", {
  # Test the gsea_nSample branch (lines 406-414)
  gl <- setNames(
    c(5, 4, 3, 2, 1, 0, -1, -2, -3, -4),
    c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J")
  )
  tg <- c("A", "B", "C", "D", "E", "F", "G", "H")

  set.seed(42)
  out <- looplook:::run_gsea_analysis(tg, gl, gsea_nSample = 4, "Test_ntop")
  expect_type(out, "list")
  expect_true("result" %in% names(out))
})

test_that("run_gsea_analysis with valid <= gsea_nSample", {
  # When valid genes <= gsea_nSample, should use all valid genes (line 413)
  gl <- setNames(c(5, 4, 3), c("A", "B", "C"))
  tg <- c("A", "B", "C")

  set.seed(42)
  out <- looplook:::run_gsea_analysis(tg, gl, gsea_nSample = 10, "Test_all")
  expect_type(out, "list")
})

test_that("run_gsea_analysis with all-zero weights", {
  # When all weights are 0, should use uniform sampling (line 410)
  gl <- setNames(rep(0, 10), LETTERS[1:10])
  tg <- LETTERS[1:8]

  set.seed(42)
  out <- looplook:::run_gsea_analysis(tg, gl, gsea_nSample = 4, "Test_zero_wt")
  expect_type(out, "list")
})

test_that("run_gsea_analysis with < 2 valid targets returns NULL", {
  gl <- setNames(c(5, 4, 3), c("A", "B", "C"))
  # Only 1 valid target
  out <- looplook:::run_gsea_analysis(c("A", "FAKE"), gl, NULL, "Test_few")
  expect_null(out$result)
  expect_null(out$plot)
})

test_that("run_lfc_violin with t.test", {
  gl <- setNames(rnorm(100, 0, 1), paste0("Gene", 1:100))
  tg <- names(gl)[1:10]
  p <- looplook:::run_lfc_violin(tg, gl, "t.test", "Test_ttest")
  expect_s3_class(p, "ggplot")
})

test_that("run_lfc_violin with very few genes returns NULL", {
  gl <- setNames(c(1, 2), c("A", "B"))
  p <- looplook:::run_lfc_violin("A", gl, "wilcox.test", "Test_single")
  expect_null(p)
})

test_that("plot_summary_go_lollipop handles empty ONTOLOGY", {
  go_df <- data.frame(
    pvalue = c(0.01, 0.02),
    Count = c(10, 5),
    Description = c("term1", "term2"),
    stringsAsFactors = FALSE
  )
  result <- looplook:::plot_summary_go_lollipop(list(go_df), "Test")
  # Should handle missing ONTOLOGY gracefully
  expect_type(result, "list")
})

test_that("plot_summary_go_lollipop with all-zero Count", {
  go_df <- data.frame(
    ONTOLOGY = c("BP", "BP"),
    pvalue = c(0.01, 0.02),
    Count = c(0, 0),
    Description = c("term1", "term2"),
    stringsAsFactors = FALSE
  )
  # Should not error on zero Count (scale_f protection)
  expect_no_error(
    looplook:::plot_summary_go_lollipop(list(go_df), "Test_zero")
  )
})

test_that("extract_target_gene_sets with NULL target_annotation", {
  anno <- list(
    target_annotation = NULL,
    loop_annotation = data.frame(
      loop_type = c("E-P"),
      Putative_Target_Genes = c("GENE_X"),
      stringsAsFactors = FALSE
    )
  )
  res <- looplook:::extract_target_gene_sets(anno, c("targets", "loops"))
  expect_true("EP_Genes" %in% names(res))
  expect_false("Target_Genes" %in% names(res))
})

test_that("extract_target_gene_sets error on missing column", {
  anno <- list(
    target_annotation = data.frame(wrong_col = "A", stringsAsFactors = FALSE)
  )
  expect_error(
    looplook:::extract_target_gene_sets(anno, "targets"),
    "Required column"
  )
})

test_that(".subset_motif_loop_df returns full df for non-loops source", {
  df <- data.frame(loop_type = c("E-P", "P-P"), stringsAsFactors = FALSE)
  r <- looplook:::.subset_motif_loop_df(df, "targets", "EP_Genes")
  expect_equal(nrow(r), 2)
})

test_that(".subset_motif_loop_df returns full df for no match", {
  df <- data.frame(loop_type = c("E-P", "P-P"), stringsAsFactors = FALSE)
  r <- looplook:::.subset_motif_loop_df(df, "loops", "ZZ_Genes")
  expect_equal(nrow(r), 2)
})
