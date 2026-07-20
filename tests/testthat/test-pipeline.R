# tests/testthat/test-pipeline.R

test_that("Module 4: profile_target_genes exhaustive branching", {
  rdata_path <- system.file("extdata", "analysis_results.RData", package = "looplook")
  expr_path <- system.file("extdata", "example_tpm.txt", package = "looplook")
  diff_path <- system.file("extdata", "example_deg.txt", package = "looplook")
  meta_path <- system.file("extdata", "example_coldata.txt", package = "looplook")
  skip_if(
    rdata_path == "" || expr_path == "" || diff_path == "" || meta_path == "",
    "Test data not available"
  )
  skip_if_not_installed("org.Hs.eg.db")
  skip_if_not_installed("clusterProfiler")

  temp_env <- new.env()
  load(rdata_path, envir = temp_env)
  raw_annotation <- temp_env[[ls(temp_env)[1]]]
  raw_annotation$loop_annotation <- head(raw_annotation$loop_annotation, 20)
  raw_annotation$target_annotation <- head(raw_annotation$target_annotation, 10)
  raw_annotation$promoter_centric_stats <- head(raw_annotation$promoter_centric_stats, 20)
  raw_annotation$distal_element_stats <- head(raw_annotation$distal_element_stats, 20)
  res_A <- looplook:::.with_messages_silenced(
    looplook::profile_target_genes(
      annotation_res = raw_annotation, diff_file = diff_path, lfc_col = "log2FoldChange",
      expr_matrix_file = expr_path, metadata_file = meta_path, target_source = "targets",
      target_mapping_mode = "all", include_Filled = TRUE, project_name = "Test_A",
      run_go = FALSE, run_ppi = FALSE, run_motif = FALSE, stat_test = "t.test", gsea_nSample = 20
    )
  )
  expect_type(res_A, "list")

  res_B <- tryCatch(
    looplook:::.with_messages_silenced(
      looplook::profile_target_genes(
        annotation_res = raw_annotation, diff_file = diff_path, lfc_col = "log2FoldChange",
        expr_matrix_file = expr_path, metadata_file = meta_path, target_source = "targets",
        target_mapping_mode = "promoter", include_Filled = FALSE, project_name = "Test_B",
        run_go = FALSE, run_ppi = FALSE, run_motif = FALSE, stat_test = "t.test", gsea_nSample = 20
      )
    ),
    error = function(e) NULL
  )
  expect_true(!is.null(res_B))

  res_C <- tryCatch(
    looplook:::.with_messages_silenced(
      looplook::profile_target_genes(
        annotation_res = raw_annotation, diff_file = diff_path, lfc_col = "log2FoldChange",
        expr_matrix_file = expr_path, metadata_file = meta_path, target_source = "targets",
        target_mapping_mode = "all", include_Filled = FALSE, use_nearest_gene = TRUE, project_name = "Test_C",
        run_go = FALSE, run_ppi = FALSE, run_motif = FALSE, stat_test = "t.test", gsea_nSample = 20
      )
    ),
    error = function(e) NULL
  )
  expect_true(!is.null(res_C))

  res_D <- looplook:::.with_messages_silenced(
    looplook::profile_target_genes(
      annotation_res = raw_annotation, diff_file = diff_path, lfc_col = "log2FoldChange",
      expr_matrix_file = expr_path, metadata_file = meta_path, target_source = "loops",
      target_mapping_mode = "all", include_Filled = TRUE, project_name = "Test_D",
      run_go = FALSE, run_ppi = FALSE, run_motif = FALSE, stat_test = "t.test", gsea_nSample = 20
    )
  )
  expect_type(res_D, "list")
})
