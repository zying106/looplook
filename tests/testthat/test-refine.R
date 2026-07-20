# tests/testthat/test-refine.R
# Tests load pre-computed annotation from inst/extdata/analysis_results.RData.
# If tests fail with "missing required columns", the RData fixture was generated
# by an older pipeline version and needs to be regenerated.

test_that("Module 3: refine_loop_anchors_by_expression runs successfully", {
  rdata_path <- system.file("extdata", "analysis_results.RData", package = "looplook")
  expr_path <- system.file("extdata", "example_tpm.txt", package = "looplook")
  skip_if(rdata_path == "" || expr_path == "")

  temp_env <- new.env()
  load(rdata_path, envir = temp_env)
  res_integrated <- temp_env[[ls(temp_env)[1]]]
  res_integrated$loop_annotation$Promoter_Target_Genes <- NA_character_
  res_integrated$loop_annotation$Putative_Target_Genes <- NA_character_
  out_base <- tempdir()


  res_basic <- looplook:::.with_known_upstream_noise_suppressed(
    refine_loop_anchors_by_expression(
      annotation_res = res_integrated,
      expr_matrix_file = expr_path,
      sample_columns = c("con1", "con2"),
      threshold = 1.0,
      unit_type = "TPM",
      reclassify_by_expression = FALSE,
      out_dir = out_base,
      project_name = "Test_Basic_Filter",
      write_output = FALSE,
      quiet = TRUE
    )
  )
  expect_type(res_basic, "list")
  expect_true(all(c("plots", "plot_list") %in% names(res_basic)))
  expect_identical(res_basic$plots, res_basic$plot_list)


  res_reclass <- looplook:::.with_known_upstream_noise_suppressed(
    refine_loop_anchors_by_expression(
      annotation_res = res_integrated,
      expr_matrix_file = expr_path,
      sample_columns = c("con1", "con2"),
      threshold = 1.0,
      unit_type = "TPM",
      reclassify_by_expression = TRUE,
      out_dir = out_base,
      project_name = "Test_Reclass_Filter",
      write_output = FALSE,
      quiet = TRUE
    )
  )
  expect_type(res_reclass, "list")
  expect_true(all(c("plots", "plot_list") %in% names(res_reclass)))
  expect_identical(res_reclass$plots, res_reclass$plot_list)
})

test_that("refine_loop_anchors_by_expression respects quiet and write_output flags", {
  rdata_path <- system.file("extdata", "analysis_results.RData", package = "looplook")
  expr_path <- system.file("extdata", "example_tpm.txt", package = "looplook")
  skip_if(rdata_path == "" || expr_path == "")

  temp_env <- new.env()
  load(rdata_path, envir = temp_env)
  res_integrated <- temp_env[[ls(temp_env)[1]]]
  res_integrated$loop_annotation$Promoter_Target_Genes <- NA_character_
  res_integrated$loop_annotation$Putative_Target_Genes <- NA_character_
  out_base <- tempfile(pattern = "looplook_refine_nowrite_")
  unlink(out_base, recursive = TRUE, force = TRUE)
  expect_false(dir.exists(out_base))

  expect_no_message({
    res_basic <- looplook:::.with_known_upstream_noise_suppressed(
      refine_loop_anchors_by_expression(
        annotation_res = res_integrated,
        expr_matrix_file = expr_path,
        sample_columns = c("con1", "con2"),
        threshold = 1.0,
        unit_type = "TPM",
        reclassify_by_expression = FALSE,
        out_dir = out_base,
        project_name = "Test_Basic_NoWrite",
        write_output = FALSE,
        quiet = TRUE
      )
    )
  })

  expect_type(res_basic, "list")
  expect_false(dir.exists(out_base))
})

test_that("refined workbook hides internal loop ids and uses consistent distal sheet name", {
  rdata_path <- system.file("extdata", "analysis_results.RData", package = "looplook")
  expr_path <- system.file("extdata", "example_tpm.txt", package = "looplook")
  skip_if(rdata_path == "" || expr_path == "")

  temp_env <- new.env()
  load(rdata_path, envir = temp_env)
  res_integrated <- temp_env[[ls(temp_env)[1]]]
  res_integrated$loop_annotation$Promoter_Target_Genes <- NA_character_
  res_integrated$loop_annotation$Putative_Target_Genes <- NA_character_
  out_base <- tempfile(pattern = "looplook_refine_export_")
  unlink(out_base, recursive = TRUE, force = TRUE)

  looplook:::.with_known_upstream_noise_suppressed(
    refine_loop_anchors_by_expression(
      annotation_res = res_integrated,
      expr_matrix_file = expr_path,
      sample_columns = c("con1", "con2"),
      threshold = 1.0,
      unit_type = "TPM",
      reclassify_by_expression = TRUE,
      out_dir = out_base,
      project_name = "Test_Export_Filtered",
      write_output = TRUE,
      quiet = TRUE
    )
  )

  workbook_file <- file.path(out_base, "Test_Export_Filtered_Refined_Results.xlsx")
  expect_true(file.exists(workbook_file))

  sheet_names <- openxlsx::getSheetNames(workbook_file)
  expect_true("Filtered Distal Element Stats" %in% sheet_names)
  expect_false("Filtered Distal Stats" %in% sheet_names)

  loop_export <- openxlsx::read.xlsx(workbook_file, sheet = "Filtered Loop Annotation")
  expect_false(any(c("a1_id", "a2_id") %in% colnames(loop_export)))

  unlink(out_base, recursive = TRUE, force = TRUE)
})

test_that("refinement status columns are consistent and mutually exclusive", {
  rdata_path <- system.file("extdata", "analysis_results.RData", package = "looplook")
  expr_path <- system.file("extdata", "example_tpm.txt", package = "looplook")
  skip_if(rdata_path == "" || expr_path == "")

  temp_env <- new.env()
  load(rdata_path, envir = temp_env)
  res_integrated <- temp_env[[ls(temp_env)[1]]]
  res_integrated$loop_annotation$Promoter_Target_Genes <- NA_character_
  res_integrated$loop_annotation$Putative_Target_Genes <- NA_character_

  refined <- looplook:::.with_known_upstream_noise_suppressed(
    refine_loop_anchors_by_expression(
      annotation_res = res_integrated,
      expr_matrix_file = expr_path,
      sample_columns = c("con1", "con2"),
      threshold = 1.0,
      reclassify_by_expression = TRUE,
      out_dir = tempdir(),
      project_name = "Test_Status",
      write_output = FALSE,
      quiet = TRUE
    )
  )
  la <- refined$loop_annotation

  # Status columns exist (Active_Target_Genes removed in pipeline upgrade)
  expect_true("Has_Active_Target" %in% colnames(la))
  expect_true("Retained_In_Functional_Network" %in% colnames(la))
  expect_true("Refinement_Action" %in% colnames(la))

  # Has_Active_Target == Retained_In_Functional_Network (always equal)
  expect_identical(la$Has_Active_Target, la$Retained_In_Functional_Network)

  # Refinement_Action categories
  valid_actions <- c(
    "retained_active_target", "reclassified_silent_anchor",
    "expression_filtered_no_active_target", "structural_only_no_active_target",
    "expression_state_unknown"
  )
  expect_true(all(la$Refinement_Action %in% valid_actions))
  expect_false(any(is.na(la$Refinement_Action)))
})

test_that("Functional Loop Annotation sheet always exists and matches retained count", {
  rdata_path <- system.file("extdata", "analysis_results.RData", package = "looplook")
  expr_path <- system.file("extdata", "example_tpm.txt", package = "looplook")
  skip_if(rdata_path == "" || expr_path == "")

  temp_env <- new.env()
  load(rdata_path, envir = temp_env)
  res_integrated <- temp_env[[ls(temp_env)[1]]]
  res_integrated$loop_annotation$Promoter_Target_Genes <- NA_character_
  res_integrated$loop_annotation$Putative_Target_Genes <- NA_character_
  out_base <- tempfile(pattern = "looplook_func_sheet_")
  unlink(out_base, recursive = TRUE, force = TRUE)

  looplook:::.with_known_upstream_noise_suppressed(
    refine_loop_anchors_by_expression(
      annotation_res = res_integrated,
      expr_matrix_file = expr_path,
      sample_columns = c("con1", "con2"),
      threshold = 1.0,
      reclassify_by_expression = TRUE,
      out_dir = out_base,
      project_name = "Test_Func_Sheet",
      write_output = TRUE,
      quiet = TRUE
    )
  )

  wb_file <- file.path(out_base, "Test_Func_Sheet_Filtered_Refined_Results.xlsx")
  expect_true(file.exists(wb_file))

  sheet_names <- openxlsx::getSheetNames(wb_file)
  expect_true("Functional Loop Annotation" %in% sheet_names)

  all_loops <- openxlsx::read.xlsx(wb_file, sheet = "Filtered Loop Annotation")
  func_loops <- openxlsx::read.xlsx(wb_file, sheet = "Functional Loop Annotation")

  # Functional sheet row count == retained count
  retained_n <- sum(all_loops$Retained_In_Functional_Network == TRUE, na.rm = TRUE)
  expect_equal(nrow(func_loops), retained_n)

  # All rows in functional sheet have Retained == TRUE
  if (nrow(func_loops) > 0) {
    expect_true(all(func_loops$Retained_In_Functional_Network == TRUE))
  }

  # Inactive loops are NOT in functional sheet
  inactive_in_func <- func_loops$Retained_In_Functional_Network == FALSE
  expect_false(any(inactive_in_func))

  unlink(out_base, recursive = TRUE, force = TRUE)
})

test_that("refine survives missing cluster_id and gracefully handles cluster-optional plots", {
  rdata_path <- system.file("extdata", "analysis_results.RData", package = "looplook")
  expr_path <- system.file("extdata", "example_tpm.txt", package = "looplook")
  tmp <- new.env()
  load(rdata_path, envir = tmp)
  res <- tmp[[ls(tmp)[1]]]

  # Remove cluster_id to simulate a stripped annotation object
  res$loop_annotation$cluster_id <- NULL
  expect_false("cluster_id" %in% colnames(res$loop_annotation))

  refined <- looplook:::refine_loop_anchors_by_expression(
    annotation_res = res,
    expr_matrix_file = expr_path,
    sample_columns = "con1",
    threshold = 1.0,
    reclassify_by_expression = TRUE,
    out_dir = tempdir(),
    project_name = "TestNoCluster",
    write_output = FALSE,
    quiet = TRUE
  )

  # Main flow must succeed and re-add cluster_id as NA
  expect_true("cluster_id" %in% colnames(refined$loop_annotation))

  # Cluster-dependent donut: may be silently skipped or built with fallback;
  # the current implementation gracefully handles all-NA cluster_id.
  expect_true(is.null(refined$plots$Target_Loop_Donut) ||
    inherits(refined$plots$Target_Loop_Donut, "ggplot"))

  # Cluster-independent plots must still be present
  expect_true(!is.null(refined$plots$Comparison_Dumbbell))
  expect_true(!is.null(refined$plots$Rose))
})
