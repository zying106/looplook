# tests/testthat/test-pipeline.R

test_that("Module 4: profile_target_genes exhaustive branching", {
  rdata_path <- system.file("extdata", "analysis_results.RData", package = "looplook")
  expr_path <- system.file("extdata", "example_tpm.txt", package = "looplook")
  diff_path <- system.file("extdata", "example_deg.txt", package = "looplook")
  meta_path <- system.file("extdata", "example_coldata.txt", package = "looplook")
  skip_if(rdata_path == "" || expr_path == "")

  temp_env <- new.env()
  load(rdata_path, envir = temp_env)
  raw_annotation <- temp_env[[ls(temp_env)[1]]]
  out_base <- tempdir()


  res_A <- suppressWarnings(suppressMessages(profile_target_genes(
    annotation_res = raw_annotation, diff_file = diff_path, lfc_col = "log2FoldChange",
    expr_matrix_file = expr_path, metadata_file = meta_path, target_source = "targets",
    target_mapping_mode = "all", include_Filled = TRUE, project_name = "Test_A",
    out_dir = out_base, run_go = TRUE, run_ppi = FALSE, run_motif = FALSE, gsea_nSample = 50
  )))
  expect_type(res_A, "list")

  res_B <- suppressWarnings(suppressMessages(profile_target_genes(
    annotation_res = raw_annotation, diff_file = diff_path, lfc_col = "log2FoldChange",
    expr_matrix_file = expr_path, metadata_file = meta_path, target_source = "targets",
    target_mapping_mode = "promoter", include_Filled = FALSE, project_name = "Test_B",
    out_dir = out_base, run_go = FALSE, run_ppi = FALSE, run_motif = FALSE, gsea_nSample = 50
  )))
  expect_type(res_B, "list")

  res_C <- suppressWarnings(suppressMessages(profile_target_genes(
    annotation_res = raw_annotation, diff_file = diff_path, lfc_col = "log2FoldChange",
    expr_matrix_file = expr_path, metadata_file = meta_path, target_source = "targets",
    target_mapping_mode = "all", include_Filled = FALSE, use_nearest_gene = TRUE, project_name = "Test_C",
    out_dir = out_base, run_go = FALSE, run_ppi = TRUE, run_motif = TRUE, gsea_nSample = 50
  )))
  expect_type(res_C, "list")

  res_D <- suppressWarnings(suppressMessages(profile_target_genes(
    annotation_res = raw_annotation, diff_file = diff_path, lfc_col = "log2FoldChange",
    expr_matrix_file = expr_path, metadata_file = meta_path, target_source = "loops",
    target_mapping_mode = "all", include_Filled = TRUE, project_name = "Test_D",
    out_dir = out_base, run_go = FALSE, run_ppi = FALSE, run_motif = FALSE, gsea_nSample = 50
  )))
  expect_type(res_D, "list")
})
