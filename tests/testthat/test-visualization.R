# tests/testthat/test-visualization.R

test_that("Module 5: IGV-Style Track Visualization generates plot successfully", {
  f1 <- system.file("extdata", "example_loops_1.bedpe", package = "looplook")
  atac_path <- system.file("extdata", "example_peaks.bed", package = "looplook")

  skip_if(f1 == "" || atac_path == "")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("TxDb.Hsapiens.UCSC.hg38.knownGene")
  skip_if_not_installed("org.Hs.eg.db")

  out_base <- tempdir()
  save_path <- file.path(out_base, "Test_Locus_Track.pdf")

  track_plot <- suppressWarnings(suppressMessages(
    plot_peaks_interactions(
      bedpe_file = f1,
      target_bed = atac_path,
      species = "hg38",
      score_to_alpha = TRUE,
      save_file = save_path
    )
  ))

  expect_s3_class(track_plot, "ggplot")

  expect_true(file.exists(save_path))
})
