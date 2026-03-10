# tests/testthat/test-refine.R

test_that("Module 3: refine_loop_anchors_by_expression runs successfully", {
  rdata_path <- system.file("extdata", "analysis_results.RData", package = "looplook")
  expr_path <- system.file("extdata", "example_tpm.txt", package = "looplook")
  skip_if(rdata_path == "" || expr_path == "")

  temp_env <- new.env()
  load(rdata_path, envir = temp_env)
  res_integrated <- temp_env[[ls(temp_env)[1]]]
  out_base <- tempdir()


  res_basic <- suppressWarnings(suppressMessages(
    refine_loop_anchors_by_expression(
      annotation_res = res_integrated,
      expr_matrix_file = expr_path,
      sample_columns = c("con1", "con2"),
      threshold = 1.0,
      unit_type = "TPM",
      reclassify_by_expression = FALSE,
      out_dir = out_base,
      project_name = "Test_Basic_Filter"
    )
  ))
  expect_type(res_basic, "list")


  res_reclass <- suppressWarnings(suppressMessages(
    refine_loop_anchors_by_expression(
      annotation_res = res_integrated,
      expr_matrix_file = expr_path,
      sample_columns = c("con1", "con2"),
      threshold = 1.0,
      unit_type = "TPM",
      reclassify_by_expression = TRUE,
      out_dir = out_base,
      project_name = "Test_Reclass_Filter"
    )
  ))
  expect_type(res_reclass, "list")
})
