# tests/testthat/test-mock-modules.R — mock data tests for uncovered code paths

test_that("plot_summary_go_lollipop returns plot list with mock GO results", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggtext")

  mock_go <- list(
    data.frame(
      Description = c(
        "immune response", "cell cycle", "DNA repair", "apoptosis",
        "signal transduction", "metabolism", "transport", "development"
      ),
      pvalue = c(1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 0.01, 0.05),
      p.adjust = c(1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 0.01, 0.05, 0.1),
      Count = c(50, 40, 30, 25, 20, 15, 10, 5),
      ONTOLOGY = rep(c("BP", "MF"), 4),
      CleanLoopType = c(rep("EP", 4), rep("PP", 4)),
      LoopType = c(rep("E-P", 4), rep("P-P", 4)),
      stringsAsFactors = FALSE
    )
  )

  plots <- looplook:::plot_summary_go_lollipop(mock_go, "TestMock")
  expect_type(plots, "list")
  expect_true(length(plots) >= 1)
  expect_s3_class(plots[[1]], "ggplot")
})

test_that("plot_summary_go_lollipop handles empty and NULL inputs", {
  expect_type(looplook:::plot_summary_go_lollipop(list(), "Empty"), "list")
  expect_type(looplook:::plot_summary_go_lollipop(
    list(data.frame(
      Description = character(0), pvalue = numeric(0),
      Count = integer(0), ONTOLOGY = character(0), CleanLoopType = character(0),
      stringsAsFactors = FALSE
    )), "NoRows"
  ), "list")
})

test_that("compute_refined_stats handles NULL upstream stats", {
  mock_loop <- data.frame(
    chr1 = "chr1", start1 = 1:6 * 1000, end1 = 1:6 * 1000 + 500,
    chr2 = "chr1", start2 = 1:6 * 2000, end2 = 1:6 * 2000 + 500,
    a1_id = paste0("A", 1:6), a2_id = paste0("B", 1:6),
    anchor1_type = c("P", "P", "E", "E", "eP", "eG"),
    anchor2_type = c("P", "E", "P", "G", "P", "P"),
    anchor1_gene = c("G1", "G2", "G3", "G4", "G5", "G6"),
    anchor2_gene = c("G7", "G8", "G1", "G9", "G10", "G11"),
    loop_type = c("P-P", "E-P", "E-P", "E-G", "eP-P", "eG-P"),
    cluster_id = paste0("C", c(1, 1, 2, 2, 3, 3)),
    Putative_Target_Genes = c("G1", "G8", "G1", "G9", "G5", "G11"),
    stringsAsFactors = FALSE
  )

  vals <- setNames(runif(11, 0, 10), paste0("G", 1:11))
  names(vals)[1] <- "G1"

  res <- looplook:::compute_refined_stats(mock_loop,
    upstream_promoter_stats = NULL,
    vals = vals, threshold = 1, hub_percentile = 0.95
  )

  expect_type(res, "list")
  expect_true("promoter_centric" %in% names(res))
  expect_true("distal_element" %in% names(res))
  if (!is.null(res$promoter_centric)) {
    expect_s3_class(res$promoter_centric, "data.frame")
  }
})

test_that("compute_refined_stats handles upstream stats merge", {
  mock_loop <- data.frame(
    chr1 = "chr1", start1 = c(1000, 2000), end1 = c(1500, 2500),
    chr2 = "chr2", start2 = c(3000, 4000), end2 = c(3500, 4500),
    a1_id = c("A1", "A2"), a2_id = c("B1", "B2"),
    anchor1_type = c("P", "P"), anchor2_type = c("E", "E"),
    anchor1_gene = c("TP53", "BRCA1"), anchor2_gene = c("MYC", "EGFR"),
    loop_type = c("E-P", "E-P"), cluster_id = c("C1", "C1"),
    Putative_Target_Genes = c("TP53", "BRCA1"),
    stringsAsFactors = FALSE
  )

  upstream_prom <- data.frame(
    Gene = c("TP53", "BRCA1"),
    Total_Loops = c(5, 3),
    n_Linked_Promoters = c(2, 1),
    n_Linked_Distal = c(3, 2),
    stringsAsFactors = FALSE
  )

  vals <- c(TP53 = 8, BRCA1 = 3, MYC = 5, EGFR = 2)
  res <- looplook:::compute_refined_stats(mock_loop,
    upstream_promoter_stats = upstream_prom,
    vals = vals, threshold = 1, hub_percentile = 0.95
  )

  expect_false(is.null(res$promoter_centric))
  expect_true("Is_Active_Gene" %in% colnames(res$promoter_centric))
})

test_that("run_lfc_violin handles t.test and edge cases", {
  global_glist <- setNames(rnorm(100), paste0("Gene", 1:100))
  targets <- names(global_glist)[1:20]

  p <- looplook:::run_lfc_violin(targets, global_glist, "t.test", "Test_T")
  expect_s3_class(p, "ggplot")

  # few targets but >=3
  p3 <- looplook:::run_lfc_violin(targets[1:3], global_glist, "wilcox.test", "Test_3")
  expect_s3_class(p3, "ggplot")

  # <3 returns NULL
  expect_null(looplook:::run_lfc_violin(targets[1:2], global_glist, "t.test", "TooFew"))

  # invalid stat test
  expect_error(looplook:::run_lfc_violin(targets, global_glist, "invalid", "Bad"))
})

test_that("build_refinement_plots returns expected plot names", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggtext")
  skip_if_not_installed("ggrepel")

  mock_orig <- data.frame(
    loop_type = c("E-P", "P-P", "E-P"),
    stringsAsFactors = FALSE
  )
  mock_loop <- data.frame(
    loop_type = c("E-P", "P-P", "E-P"),
    chr1 = "chr1", start1 = c(1, 2, 3) * 1000, end1 = c(1, 2, 3) * 1000 + 500,
    chr2 = "chr1", start2 = c(4, 5, 6) * 1000, end2 = c(4, 5, 6) * 1000 + 500,
    a1_id = c("A1", "A2", "A3"), a2_id = c("B1", "B2", "B3"),
    cluster_id = c("C1", "C2", "C3"),
    Putative_Target_Genes = c("TP53", "BRCA1", "MYC"),
    stringsAsFactors = FALSE
  )

  plots <- looplook:::build_refinement_plots(mock_orig, mock_loop,
    bed_info = NULL, whitelist = c("TP53", "BRCA1"),
    project_name = "Test", karyo_bin_size = 1e6, species = "hg38"
  )

  expect_type(plots, "list")
  expect_true("Comparison_Dumbbell" %in% names(plots))
  expect_true("Rose" %in% names(plots))
  expect_s3_class(plots$Comparison_Dumbbell, "ggplot")
})

test_that("build_annotation_plots returns core plot names", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("TxDb.Hsapiens.UCSC.hg38.knownGene")
  skip_if_not_installed("org.Hs.eg.db")

  mock_plot_df <- data.frame(
    loop_type = c("E-P", "P-P"),
    loop_genes = c("TP53;BRCA1", "MYC"),
    All_Anchor_Genes = c("TP53;BRCA1", "MYC"),
    chr1 = "chr1", start1 = c(1, 2) * 1000, end1 = c(1, 2) * 1000 + 500,
    chr2 = "chr1", start2 = c(3, 4) * 1000, end2 = c(3, 4) * 1000 + 500,
    stringsAsFactors = FALSE
  )
  mock_cluster <- data.frame(
    annotation = c("Promoter", "Distal Intergenic"),
    cluster_id = c("C1", "C2"),
    stringsAsFactors = FALSE
  )
  mock_bed <- data.frame(
    annotation = c("Promoter", "Intron"),
    Linked_Loop_IDs = c("L1", NA),
    stringsAsFactors = FALSE
  )
  mock_tgt <- data.frame(
    loop_type = c("E-P"),
    stringsAsFactors = FALSE
  )

  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene

  plots <- looplook:::build_annotation_plots(
    plot_df = mock_plot_df, bed_info = mock_bed,
    cluster_info = mock_cluster, target_connected_loops = mock_tgt,
    txdb_obj = txdb, org_db_pkg = "org.Hs.eg.db",
    species = "hg38", project_name = "TestAnno",
    color_palette = "Set2", karyo_bin_size = 1e7
  )

  expect_type(plots, "list")
  expect_true("Basic_Donut" %in% names(plots))
  expect_true("Basic_Circular" %in% names(plots))
  expect_s3_class(plots$Basic_Donut, "ggplot")
})
