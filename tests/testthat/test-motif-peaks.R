# tests/testthat/test-motif-peaks.R — peak-overlap stratification helpers

test_that(".load_motif_peak_gr reads a narrowPeak file into GRanges", {
  pk <- tempfile(fileext = ".narrowPeak")
  writeLines(c(
    "chr1\t100\t200\tpeak1\t10\t.\t1\t2\t3\t50",
    "chr2\t300\t400\tpeak2\t20\t.\t1\t2\t3\t80"
  ), pk)

  gr <- looplook:::.load_motif_peak_gr(pk)
  expect_s4_class(gr, "GRanges")
  expect_equal(length(gr), 2)
  expect_equal(as.character(GenomicRanges::seqnames(gr)), c("chr1", "chr2"))
  expect_equal(GenomicRanges::start(gr), c(101L, 301L)) # 0-based -> 1-based
  expect_equal(GenomicRanges::end(gr), c(200L, 400L))

  # GRanges passthrough
  gr2 <- looplook:::.load_motif_peak_gr(gr)
  expect_identical(length(gr2), 2L)

  # duplicate rows collapse to unique intervals
  pk2 <- tempfile(fileext = ".bed")
  writeLines(c("chr1\t100\t200", "chr1\t100\t200"), pk2)
  gr3 <- looplook:::.load_motif_peak_gr(pk2)
  expect_equal(length(gr3), 1)

  # invalid input errors
  expect_error(looplook:::.load_motif_peak_gr(tempfile(fileext = ".txt")),
    "not found")
  bad <- tempfile(fileext = ".bed")
  writeLines("chr1\t100", bad)
  expect_error(looplook:::.load_motif_peak_gr(bad), "at least three")
  expect_error(looplook:::.load_motif_peak_gr(42), "GRanges object or a path")
})

test_that(".make_peak_overlap_fg windows anchors at overlapping peak midpoints", {
  fg <- GenomicRanges::GRanges(
    c("chr1", "chr1", "chr2"),
    IRanges::IRanges(c(1, 1000, 1), c(500, 1500, 500))
  )
  # chr1:500-600 overlaps anchor 1; chr1:1100-1200 overlaps anchor 2;
  # chr2:900-1000 overlaps nothing.
  pk <- GenomicRanges::GRanges(
    c("chr1", "chr1", "chr2"),
    IRanges::IRanges(c(500, 1100, 900), c(600, 1200, 1000))
  )

  win <- looplook:::.make_peak_overlap_fg(fg, pk)
  expect_s4_class(win, "GRanges")
  expect_equal(length(win), 2)
  # windows are centred on peak midpoints (550 and 1150), 500 bp wide
  # (even-width resize centres at midpoint - 250)
  expect_equal(GenomicRanges::start(win), c(300L, 900L))
  expect_equal(GenomicRanges::end(win), c(799L, 1399L))

  # identical windows deduplicated
  fg2 <- GenomicRanges::GRanges(c("chr1", "chr1"), IRanges::IRanges(c(1, 1), c(700, 700)))
  pk2 <- GenomicRanges::GRanges("chr1", IRanges::IRanges(500, 600))
  win2 <- looplook:::.make_peak_overlap_fg(fg2, pk2)
  expect_equal(length(win2), 1)

  # no overlap -> empty GRanges
  pk3 <- GenomicRanges::GRanges("chr9", IRanges::IRanges(1, 100))
  win3 <- suppressWarnings(looplook:::.make_peak_overlap_fg(fg, pk3))
  expect_s4_class(win3, "GRanges")
  expect_equal(length(win3), 0)

  # NULL-safe
  expect_equal(length(looplook:::.make_peak_overlap_fg(NULL, pk)), 0)
})

test_that(".make_peak_overlap_fg picks the peak with the largest overlap", {
  fg <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1000, 2000))
  # small overlap (1800-2000) and large overlap (1000-2500): the large one wins
  pk <- GenomicRanges::GRanges(
    c("chr1", "chr1"),
    IRanges::IRanges(c(1800, 1000), c(2000, 2500))
  )
  win <- looplook:::.make_peak_overlap_fg(fg, pk)
  expect_equal(length(win), 1)
  # centred on the midpoint of the 1000-2500 peak (1750)
  expect_equal(GenomicRanges::start(win), 1500L)
  expect_equal(GenomicRanges::end(win), 1999L)
})

test_that(".make_peak_overlap_fg inherits cluster_id from source anchors", {
  fg <- GenomicRanges::GRanges(
    c("chr1", "chr1"),
    IRanges::IRanges(c(1, 1000), c(500, 1500))
  )
  S4Vectors::mcols(fg)$cluster_id <- c("c1", "c2")
  pk <- GenomicRanges::GRanges(
    c("chr1", "chr1"),
    IRanges::IRanges(c(200, 1100), c(300, 1200))
  )
  win <- looplook:::.make_peak_overlap_fg(fg, pk)
  expect_equal(length(win), 2)
  expect_equal(S4Vectors::mcols(win)$cluster_id, c("c1", "c2"))

  # no cluster_id on anchors -> no cluster_id on windows
  win2 <- looplook:::.make_peak_overlap_fg(
    GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 500)), pk
  )
  expect_null(S4Vectors::mcols(win2)$cluster_id)
})

test_that("profile_target_genes accepts the motif_n_perm auto sentinel", {
  dummy <- list(
    annotation_res = list(),
    diff_file = "missing_diff.csv",
    expr_matrix_file = "missing_expr.txt",
    metadata_file = "missing_meta.txt"
  )
  # -2 (or any integer < -1) is rejected by the public guard
  expect_error(
    do.call(looplook::profile_target_genes, c(dummy, list(motif_n_perm = -2L))),
    "finite integer >= -1"
  )
  # -1 passes the guard and fails later for an unrelated reason
  err <- tryCatch(
    do.call(looplook::profile_target_genes, c(dummy, list(motif_n_perm = -1L))),
    error = function(e) conditionMessage(e)
  )
  expect_false(grepl("finite integer", err))
})

test_that(".motif_peak_qc reports overlap statistics", {
  fg <- GenomicRanges::GRanges(
    c("chr1", "chr1", "chr2"),
    IRanges::IRanges(c(1, 1000, 1), c(500, 1500, 500))
  )
  pk <- GenomicRanges::GRanges(
    c("chr1", "chr2"),
    IRanges::IRanges(c(500, 900), c(600, 1000))
  )

  qc <- looplook:::.motif_peak_qc(fg, pk, "distal")
  expect_s3_class(qc, "data.frame")
  expect_equal(qc$Set, "distal")
  expect_equal(qc$N_FG_Anchors, 3)
  expect_equal(qc$N_FG_Overlap_Peaks, 1)
  expect_equal(qc$FG_Overlap_Fraction, 1 / 3)
  expect_equal(qc$N_Peaks_Total, 2)
  expect_equal(qc$N_Peaks_Overlap_Anchors, 1)

  expect_null(looplook:::.motif_peak_qc(NULL, pk, "x"))
  expect_null(looplook:::.motif_peak_qc(fg, NULL, "x"))
})

test_that(".resolve_motif_n_perm handles the auto sentinel", {
  gr_with_ids <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 100))
  S4Vectors::mcols(gr_with_ids)$cluster_id <- "c1"
  gr_without_ids <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 100))

  sets_with <- list(proximal_fg = gr_with_ids, distal_fg = gr_without_ids)
  sets_without <- list(proximal_fg = gr_without_ids, distal_fg = gr_without_ids)

  expect_equal(looplook:::.resolve_motif_n_perm(0L, sets_with), 0L)
  expect_equal(looplook:::.resolve_motif_n_perm(50L, sets_with), 50L)
  expect_equal(looplook:::.resolve_motif_n_perm(-1L, sets_with), 100L)
  expect_equal(looplook:::.resolve_motif_n_perm(-1L, sets_without), 0L)
  expect_equal(looplook:::.resolve_motif_n_perm(NA_integer_, sets_with), 0L)
})

test_that("run_distal_motif_analysis keeps the stable schema without peaks", {
  empty_loop_df <- data.frame(
    chr1 = character(0), start1 = integer(0), end1 = integer(0),
    chr2 = character(0), start2 = integer(0), end2 = integer(0),
    stringsAsFactors = FALSE
  )
  out <- suppressWarnings(
    looplook:::run_distal_motif_analysis(
      target_genes = "TP53", loop_df = empty_loop_df,
      genome_id = "hg38", pval_thresh = 1e-4,
      current_proj_name = "Test", peak_gr = NULL
    )
  )
  expect_named(out, c("results", "plots"))
  expect_named(out$results, c("proximal", "distal"))
  expect_false("qc" %in% names(out))
})

test_that(".plot_save_motif tolerates missing Log2OddsRatio column", {
  skip_if_not_installed("ggplot2")
  mock_res <- data.frame(
    MotifID = paste0("MA", sprintf("%04d", 1:5), ".1"),
    MotifName = paste("TF", LETTERS[1:5]),
    Pvalue = 10^-(1:5),
    FDR = 10^-(1:5),
    OddsRatio = c(0.5, 1, 1.2, 2, 4),
    stringsAsFactors = FALSE
  )
  p <- looplook:::.plot_save_motif(mock_res, "Test_Bar")
  expect_s3_class(p, "ggplot")

  p2 <- looplook:::.plot_motif_rank_scatter(mock_res, "Test_Rank")
  expect_s3_class(p2, "ggplot")
})

test_that(".plot_save_motif guards FDR = 0 and NA", {
  skip_if_not_installed("ggplot2")
  mock_res <- data.frame(
    MotifID = paste0("MA", sprintf("%04d", 1:5), ".1"),
    MotifName = paste("TF", LETTERS[1:5]),
    Pvalue = c(0, 1e-10, 1e-5, 0.5, 1),
    FDR = c(0, 1e-10, NA, 0.5, 1),
    OddsRatio = c(30, 2, 1.5, 1, 0.8),
    stringsAsFactors = FALSE
  )
  p <- looplook:::.plot_save_motif(mock_res, "Test_Bar")
  expect_s3_class(p, "ggplot")
  pd <- p$data
  expect_false(any(is.na(pd$FDR_safe) | !is.finite(pd$FDR_safe)))
  expect_equal(pd$FDR_safe[!is.na(pd$FDR) & pd$FDR == 0], 1e-300)
  expect_equal(pd$FDR_safe[is.na(pd$FDR)], 1)
})

test_that(".plot_save_motif keeps the x-axis non-negative", {
  skip_if_not_installed("ggplot2")
  # degenerate case: all top bars have FDR = 1 (x = 0)
  mock_res <- data.frame(
    MotifID = paste0("MA", sprintf("%04d", 1:5), ".1"),
    MotifName = paste("TF", LETTERS[1:5]),
    Pvalue = rep(1, 5),
    FDR = rep(1, 5),
    OddsRatio = rep(1.2, 5),
    stringsAsFactors = FALSE
  )
  p <- looplook:::.plot_save_motif(mock_res, "Test_Bar")
  expect_s3_class(p, "ggplot")
  x_scale <- p$scales$get_scales("x")
  # expansion vector: c(mult_lower, add_lower, mult_upper, add_upper)
  expect_equal(x_scale$expand[1], 0)
  expect_equal(x_scale$expand[3], 0.05)
})

test_that(".plot_motif_rank_scatter keeps non-significant dots small", {
  skip_if_not_installed("ggplot2")
  mock_res <- data.frame(
    MotifID = paste0("MA", sprintf("%04d", 1:8), ".1"),
    MotifName = paste("TF", LETTERS[1:8]),
    Pvalue = 10^-(1:8),
    FDR = c(1e-10, 1e-10, 1e-4, 1e-12, 0.8, 0.9, 0.95, 1),
    OddsRatio = c(3, 1.5, 0.25, 30, 8, 2, 1, 0.5),
    Family = "Unknown",
    stringsAsFactors = FALSE
  )
  p <- looplook:::.plot_motif_rank_scatter(mock_res, "Test_Rank", fdr_thresh = 0.05)
  expect_s3_class(p, "ggplot")
  pd <- p$data
  sig <- pd$FDR < 0.05
  expect_true(all(pd$SizeVal[sig] > 0))
  expect_true(all(pd$SizeVal[!sig] == 0))
  # FDR ties (both 1e-10) are broken by |log2(OR)| descending
  expect_lt(which(pd$OddsRatio == 3), which(pd$OddsRatio == 1.5))
})

test_that(".plot_save_motif guards FDR = 0 and NA", {
  skip_if_not_installed("ggplot2")
  mock_res <- data.frame(
    MotifID = paste0("MA", sprintf("%04d", 1:5), ".1"),
    MotifName = paste("TF", LETTERS[1:5]),
    Pvalue = c(0, 1e-10, 1e-5, 0.5, 1),
    FDR = c(0, 1e-10, NA, 0.5, 1),
    OddsRatio = c(30, 2, 1.5, 1, 0.8),
    stringsAsFactors = FALSE
  )
  p <- looplook:::.plot_save_motif(mock_res, "Test_Bar")
  expect_s3_class(p, "ggplot")
  pd <- p$data
  expect_false(any(is.na(pd$FDR_safe) | !is.finite(pd$FDR_safe)))
  expect_equal(pd$FDR_safe[!is.na(pd$FDR) & pd$FDR == 0], 1e-300)
  expect_equal(pd$FDR_safe[is.na(pd$FDR)], 1)
})
