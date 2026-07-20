# tests/testthat/test-utils_plots.R

test_that("draw_circular_bar_plot returns ggplot with valid input", {
  skip_if_not_installed("ggplot2")
  mock <- data.frame(
    loop_type = c("E-P", "P-P", "E-E"),
    loop_genes = c("GeneA;GeneB", "GeneC", "GeneD;GeneE"),
    stringsAsFactors = FALSE
  )
  p <- looplook:::draw_circular_bar_plot(mock, "Test", NULL, c("E-P" = "red", "P-P" = "blue", "E-E" = "green"))
  expect_s3_class(p, "ggplot")
})

test_that("draw_pie_with_outside_labels handles empty and valid data", {
  skip_if_not_installed("ggplot2")
  expect_null(looplook:::draw_pie_with_outside_labels(data.frame(), "anno", "Empty", "Set2"))

  df <- data.frame(
    annotation = c("Promoter (<=1kb)", "Distal Intergenic", "Intron (1-2kb)"),
    stringsAsFactors = FALSE
  )
  p <- looplook:::draw_pie_with_outside_labels(df, "annotation", "Test", "Set2")
  expect_s3_class(p, "ggplot")
})
