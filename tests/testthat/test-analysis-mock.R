# tests/testthat/test-analysis-mock.R — mock-only tests for analysis.R internal functions

test_that(".plot_motif_rank_scatter generates scatter with mock data", {
  skip_if_not_installed("ggplot2")

  mock_res <- data.frame(
    MotifID = paste0("MA", sprintf("%04d", 1:20), ".1"),
    MotifName = paste("TF", LETTERS[1:20]),
    Pvalue = 10^-(1:20),
    FDR = c(10^-(1:10), rep(0.5, 10)),
    OddsRatio = runif(20, 0.8, 1.6),
    Family = c(
      rep("bZIP", 5), rep("Homeobox", 5), rep("Unknown", 5),
      rep("bHLH", 5)
    ),
    stringsAsFactors = FALSE
  )

  p <- looplook:::.plot_motif_rank_scatter(mock_res, "Test_Motif", fdr_thresh = 0.05)
  expect_s3_class(p, "ggplot")

  # empty data
  p_null <- looplook:::.plot_motif_rank_scatter(
    data.frame(
      MotifID = character(0), MotifName = character(0),
      Pvalue = numeric(0), FDR = numeric(0), OddsRatio = numeric(0),
      Family = character(0), stringsAsFactors = FALSE
    ),
    "Empty"
  )
  expect_null(p_null)


  # no Family column
  res_no_fam <- mock_res[, setdiff(colnames(mock_res), "Family")]
  p2 <- looplook:::.plot_motif_rank_scatter(res_no_fam, "NoFam")
  expect_s3_class(p2, "ggplot")
})

test_that(".plot_save_motif generates barplot with mock data", {
  skip_if_not_installed("ggplot2")

  mock_res <- data.frame(
    MotifID = paste0("MA", sprintf("%04d", 1:20), ".1"),
    MotifName = paste("TF", LETTERS[1:20]),
    Pvalue = 10^-(1:20),
    FDR = 10^-(1:20),
    OddsRatio = rep(1.2, 20),
    stringsAsFactors = FALSE
  )

  p <- looplook:::.plot_save_motif(mock_res, "Test_Bar")
  expect_s3_class(p, "ggplot")

  # NULL returns NULL
  expect_null(looplook:::.plot_save_motif(NULL, "Null"))

  # empty data
  expect_null(looplook:::.plot_save_motif(
    mock_res[integer(0), ], "Empty"
  ))
})

test_that("run_heatmap_and_connectivity handles distal-only mode with mock data", {
  skip_if_not_installed("ggpointdensity")
  skip_if_not_installed("viridis")

  genes <- paste0("Gene", 1:30)
  tpm_mat <- as.data.frame(matrix(rexp(30 * 4, rate = 1), 30, 4,
    dimnames = list(genes, c("s1", "s2", "s3", "s4"))
  ))
  meta_raw <- data.frame(
    SampleID = c("s1", "s2", "s3", "s4"),
    Group = c("A", "A", "B", "B"), stringsAsFactors = FALSE
  )

  loop_stats <- data.frame(
    Gene = genes[1:15],
    Total_Loops = sample(1:10, 15, replace = TRUE),
    n_Linked_Distal = sample(0:5, 15, replace = TRUE),
    stringsAsFactors = FALSE
  )
  global_glist <- setNames(rnorm(30), genes)
  targets <- genes[1:10]

  plots <- looplook:::run_heatmap_and_connectivity(
    targets, tpm_mat, meta_raw, loop_stats, global_glist,
    heatmap_ntop = 50, cor_method = "pearson",
    current_proj_name = "Test_Distal", source_type = "targets",
    target_col = "n_Linked_Distal", skip_heatmap = TRUE
  )
  expect_type(plots, "list")
})

test_that(".annotate_motif_families handles empty input", {
  res_empty <- looplook:::.annotate_motif_families(NULL)
  expect_null(res_empty)

  res_no_rows <- looplook:::.annotate_motif_families(
    data.frame(
      MotifID = character(0), MotifName = character(0),
      Pvalue = numeric(0), stringsAsFactors = FALSE
    )
  )
  expect_equal(nrow(res_no_rows), 0)
})

test_that(".subset_motif_loop_df restricts motif universe for loop-derived tasks", {
  loop_df <- data.frame(
    loop_type = c("E-P", "P-P", "E-P", "G-P"),
    stringsAsFactors = FALSE
  )

  ep_only <- looplook:::.subset_motif_loop_df(loop_df, src = "loops", task_name = "EP_Genes")
  expect_equal(nrow(ep_only), 2)
  expect_true(all(ep_only$loop_type == "E-P"))

  unchanged_targets <- looplook:::.subset_motif_loop_df(loop_df, src = "targets", task_name = "Target_Genes")
  expect_equal(nrow(unchanged_targets), nrow(loop_df))
})

test_that(".prepare_motif_anchor_sets deduplicates anchors and matches classes", {
  loop_df <- data.frame(
    chr1 = rep("chr1", 6),
    start1 = c(100, 100, 100, 700, 900, 1100),
    end1 = c(150, 150, 150, 750, 950, 1150),
    chr2 = rep("chr1", 6),
    start2 = c(300, 500, 650, 1300, 1500, 1700),
    end2 = c(350, 550, 700, 1350, 1550, 1750),
    a1_id = c("A", "A", "A", "E", "G", "I"),
    a2_id = c("B", "C", "D", "F", "H", "J"),
    anchor1_gene = c("TP53", "TP53", "TP53", "CTRL1", "CTRL2", "CTRL3"),
    anchor1_type = c("P", "P", "P", "P", "P", "P"),
    anchor2_gene = c("EnhA", "EnhB", "PartnerProm", "BgEnh", "GeneBody", "SilentProm"),
    anchor2_type = c("E", "eP", "P", "E", "G", "eG"),
    stringsAsFactors = FALSE
  )

  motif_sets <- looplook:::.prepare_motif_anchor_sets(loop_df, target_genes = "TP53")

  expect_equal(motif_sets$target_loop_n, 3)
  expect_equal(sort(names(motif_sets$proximal_fg)), "A")
  # eP/eG no longer distal-like; only E anchors are distal for motif sets
  expect_equal(sort(names(motif_sets$distal_fg)), "B")
  expect_equal(sort(names(motif_sets$proximal_bg)), c("E", "G", "I"))
  expect_equal(sort(names(motif_sets$distal_bg)), "F")
})

# --- GSEA seed reproducibility ---
test_that("profile_target_genes with fixed seed produces reproducible output", {
  rng_capture <- new.env(parent = emptyenv())
  rng_capture$draws <- list()
  testthat::local_mocked_bindings(
    run_gsea_analysis = function(...) list(result = NULL, plot = NULL),
    run_heatmap_and_connectivity = function(...) list(),
    # Mock the heavy violin resampling with a lightweight RNG draw so the
    # test verifies that `seed` actually governs downstream RNG (not just
    # the deterministic gene-set extraction).
    run_lfc_violin = function(...) {
      rng_capture$draws[[length(rng_capture$draws) + 1L]] <-
        sample.int(1e6, 1L)
      NULL
    },
    .package = "looplook"
  )
  rdata_path <- system.file("extdata", "analysis_results.RData", package = "looplook")
  diff_path <- system.file("extdata", "example_deg.txt", package = "looplook")
  expr_path <- system.file("extdata", "example_tpm.txt", package = "looplook")
  meta_path <- system.file("extdata", "example_coldata.txt", package = "looplook")
  skip_if(
    rdata_path == "" || diff_path == "" || expr_path == "" || meta_path == "",
    "Example data not available"
  )

  tmp <- new.env()
  load(rdata_path, envir = tmp)
  res <- tmp[[ls(tmp)[1]]]

  run_with_seed <- function(s) {
    looplook::profile_target_genes(
      annotation_res = res,
      diff_file = diff_path,
      expr_matrix_file = expr_path,
      metadata_file = meta_path,
      run_go = FALSE, run_ppi = FALSE, run_motif = FALSE,
      gsea_nSample = 5, seed = s
    )
  }

  set.seed(42)
  draws1 <- c()
  rng_capture$draws <- list()
  run1 <- run_with_seed(123)
  draws1 <- unlist(rng_capture$draws)

  set.seed(99)
  rng_capture$draws <- list()
  run2 <- run_with_seed(123)
  draws2 <- unlist(rng_capture$draws)

  # Same seed must govern the RNG draws made inside the pipeline, regardless
  # of the RNG state before the call.
  expect_equal(draws1, draws2,
    info = "Same seed must produce identical downstream RNG draws"
  )
  expect_equal(
    attr(run1, "seed"),
    attr(run2, "seed"),
    info = "Run metadata must record the same seed"
  )
})

test_that("profile_target_genes with NULL seed has NULL in metadata", {
  testthat::local_mocked_bindings(
    run_gsea_analysis = function(...) list(result = NULL, plot = NULL),
    run_heatmap_and_connectivity = function(...) list(),
    .package = "looplook"
  )
  rdata_path <- system.file("extdata", "analysis_results.RData", package = "looplook")
  diff_path <- system.file("extdata", "example_deg.txt", package = "looplook")
  expr_path <- system.file("extdata", "example_tpm.txt", package = "looplook")
  meta_path <- system.file("extdata", "example_coldata.txt", package = "looplook")
  skip_if(
    rdata_path == "" || diff_path == "" || expr_path == "" || meta_path == "",
    "Example data not available"
  )

  tmp <- new.env()
  load(rdata_path, envir = tmp)
  res <- tmp[[ls(tmp)[1]]]

  out <- profile_target_genes(
    annotation_res = res,
    diff_file = diff_path,
    expr_matrix_file = expr_path,
    metadata_file = meta_path,
    run_go = FALSE, run_ppi = FALSE, run_motif = FALSE,
    gsea_nSample = 5
  )

  expect_null(attr(out, "seed"),
    info = "No seed argument should record NULL in result metadata"
  )
})
