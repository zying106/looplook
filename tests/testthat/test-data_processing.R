# tests/testthat/test-data_processing.R

test_that("data_processing modules run successfully on example bedpe files", {
  loop1 <- system.file("extdata", "example_loops_1.bedpe", package = "looplook")
  loop2 <- system.file("extdata", "example_loops_2.bedpe", package = "looplook")
  h3k27ac_peaks <- system.file("extdata", "example_k27ac_peaks.bed", package = "looplook")
  skip_if(loop1 == "" || loop2 == "")

  res_clean <- suppressWarnings(suppressMessages(
    consolidate_chromatin_loops(
      files = c(loop1, loop2),
      mode = "consensus",
      min_raw_score = 2,
      min_score = 5,
      gap = 1000,
      blacklist_species = "hg38",
      region_of_interest = h3k27ac_peaks,
      out_file = tempfile(fileext = ".bedpe")
    )
  ))

  expect_true(!is.null(res_clean))
})
