# tests/testthat/test-annotation.R

test_that("Module 2: annotate_peaks_and_loops runs comprehensive Example A", {
  global_out <- system.file("extdata", "example_loops_1.bedpe", package = "looplook")
  expr_path <- system.file("extdata", "example_tpm.txt", package = "looplook")
  atac_path <- system.file("extdata", "example_peaks.bed", package = "looplook")


  skip_if(global_out == "" || expr_path == "" || atac_path == "")
  skip_if_not_installed("TxDb.Hsapiens.UCSC.hg38.knownGene")
  skip_if_not_installed("org.Hs.eg.db")

  out_base <- tempdir()


  res_integrated <- suppressWarnings(suppressMessages(
    annotate_peaks_and_loops(
      bedpe_file = global_out,
      target_bed = atac_path,
      expr_matrix_file = expr_path,
      sample_columns = c("con1", "con2"),
      species = "hg38",
      tss_region = c(-2000, 2000),
      neighbor_hop = 0,
      hub_percentile = 0.95,
      out_dir = out_base,
      project_name = "Example_HiChIP_Integrative_Test"
    )
  ))


  expect_type(res_integrated, "list")
  expect_true(all(c("target_annotation", "loop_annotation", "anchor_annotation") %in% names(res_integrated)))


  expect_true("Assigned_Target_Genes_Filled" %in% colnames(res_integrated$target_annotation))
})
