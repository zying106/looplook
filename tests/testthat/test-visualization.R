# tests/testthat/test-visualization.R

test_that("Module 5: IGV-Style Track Visualization generates plot successfully", {
  skip_if_not_installed("ggplot2")
  f1 <- tempfile(fileext = ".bedpe")
  writeLines(
    "chr1\t100\t140\tchr1\t400\t440",
    con = f1
  )

  track_plot <- looplook::plot_peaks_interactions(
    bedpe_file = f1,
    chr = "chr1",
    from = 90,
    to = 450,
    score_to_alpha = TRUE,
    show_gene_track = FALSE
  )

  expect_s3_class(track_plot, "ggplot")
})

test_that("prepare_track_data stacks overlapping arcs within max_levels", {
  bedpe_path <- tempfile(fileext = ".bedpe")
  on.exit(unlink(bedpe_path), add = TRUE)
  writeLines(
    c(
      "chr1\t100\t140\tchr1\t400\t440",
      "chr1\t120\t160\tchr1\t320\t360",
      "chr1\t150\t190\tchr1\t260\t300"
    ),
    con = bedpe_path
  )

  d <- looplook:::prepare_track_data(
    bedpe_file = bedpe_path,
    target_bed = NULL,
    chr = "chr1",
    from = 90,
    to = 450,
    species = "hg38",
    max_levels = 2,
    base_anchor_height = 0.05,
    loop_color = "#5D6D7E",
    anchor_color = "#3498DB",
    score_to_alpha = FALSE,
    min_score = NULL,
    show_gene_track = FALSE
  )

  expect_gt(nrow(d$bez_df), 0)

  arc_levels <- sort(unique(d$bez_df$arc_level))
  expect_equal(arc_levels, c(1, 2))
  expect_lte(max(arc_levels), 2)

  peak_by_loop <- stats::aggregate(y ~ loop_i, data = d$bez_df, FUN = max)
  expect_equal(length(unique(peak_by_loop$y)), nrow(peak_by_loop))
})

test_that("draw_flower_simplified returns ggplot and handles edge cases", {
  skip_if_not_installed("ggplot2")

  gene_sets <- list(
    A = c("TP53", "BRCA1", "MYC"),
    B = c("BRCA1", "MYC", "EGFR"),
    C = c("MYC", "EGFR", "KRAS")
  )
  p <- looplook:::draw_flower_simplified(
    gene_sets, "Test_Flower",
    c(A = "red", B = "blue", C = "green")
  )
  expect_s3_class(p, "ggplot")

  p_null <- looplook:::draw_flower_simplified(list(X = c("A", "B")), "Solo", NULL)
  expect_null(p_null)
})

test_that("plot_peaks_interactions with score column", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggforce")
  skip_if_not_installed("ggrepel")

  f1 <- tempfile(fileext = ".bedpe")
  writeLines(
    c(
      "chr1\t100\t140\tchr1\t400\t440\t10",
      "chr1\t150\t190\tchr1\t350\t390\t20"
    ),
    con = f1
  )

  p <- looplook::plot_peaks_interactions(
    bedpe_file = f1,
    chr = "chr1",
    from = 90,
    to = 450,
    score_to_alpha = TRUE,
    show_gene_track = FALSE
  )
  expect_s3_class(p, "ggplot")
  unlink(f1)
})

test_that("plot_peaks_interactions with target_bed overlay", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggforce")
  skip_if_not_installed("ggrepel")

  bedpe <- tempfile(fileext = ".bedpe")
  bed <- tempfile(fileext = ".bed")
  writeLines("chr1\t100\t140\tchr1\t400\t440", bedpe)
  writeLines("chr1\t110\t130", bed)

  p <- looplook::plot_peaks_interactions(
    bedpe_file = bedpe,
    target_bed = bed,
    chr = "chr1",
    from = 90,
    to = 450,
    show_gene_track = FALSE
  )
  expect_s3_class(p, "ggplot")
  unlink(c(bedpe, bed))
})

test_that("prepare_track_data with min_score filter", {
  bedpe_path <- tempfile(fileext = ".bedpe")
  writeLines(
    c(
      "chr1\t100\t140\tchr1\t400\t440\t5",
      "chr1\t150\t190\tchr1\t350\t390\t50"
    ),
    bedpe_path
  )

  d <- looplook:::prepare_track_data(
    bedpe_file = bedpe_path,
    target_bed = NULL,
    chr = "chr1",
    from = 90,
    to = 450,
    species = "hg38",
    max_levels = 10,
    base_anchor_height = 0.05,
    loop_color = "#5D6D7E",
    anchor_color = "#3498DB",
    score_to_alpha = FALSE,
    min_score = 10,
    show_gene_track = FALSE
  )
  # Only the loop with score >= 10 should remain
  expect_equal(length(unique(d$bez_df$loop_i)), 1)
  unlink(bedpe_path)
})

test_that("prepare_track_data infers chr when NULL", {
  bedpe_path <- tempfile(fileext = ".bedpe")
  writeLines(
    c(
      "chr5\t100\t140\tchr5\t400\t440",
      "chr5\t150\t190\tchr5\t350\t390"
    ),
    bedpe_path
  )

  d <- looplook:::prepare_track_data(
    bedpe_file = bedpe_path,
    target_bed = NULL,
    chr = NULL, # Should infer chr5
    from = NULL,
    to = NULL,
    species = "hg38",
    max_levels = 10,
    base_anchor_height = 0.05,
    loop_color = "#5D6D7E",
    anchor_color = "#3498DB",
    score_to_alpha = FALSE,
    min_score = NULL,
    show_gene_track = FALSE
  )
  expect_equal(d$chr, "chr5")
  unlink(bedpe_path)
})

test_that("plot_peaks_interactions: interactive mode returns girafe", {
  skip_if_not_installed("ggiraph")
  bedpe_path <- tempfile(fileext = ".bedpe")
  writeLines("chr1\t11890000\t11890500\tchr1\t11905000\t11905500", bedpe_path)
  p <- plot_peaks_interactions(
    bedpe_file = bedpe_path, chr = "chr1", from = 11884299, to = 12106581,
    show_gene_track = FALSE, interactive = TRUE
  )
  expect_s3_class(p, "girafe")
  unlink(bedpe_path)
})

test_that("plot_peaks_interactions: interactive=FALSE unchanged", {
  bedpe_path <- tempfile(fileext = ".bedpe")
  writeLines("chr1\t11890000\t11890500\tchr1\t11905000\t11905500", bedpe_path)
  p <- plot_peaks_interactions(
    bedpe_file = bedpe_path, chr = "chr1", from = 11884299, to = 12106581,
    show_gene_track = FALSE, interactive = FALSE
  )
  expect_s3_class(p, "ggplot")
  unlink(bedpe_path)
})
