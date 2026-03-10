# tests/testthat/test-extract.R

test_that("extract_target_gene_sets correctly parses annotation lists", {
  mock_anno <- list(
    target_annotation = data.frame(
      SYMBOL = c("GENE_A", "GENE_B", NA),
      Assigned_Target_Genes_Filled = c("GENE_A;GENE_C", "GENE_B", "GENE_D")
    ),
    loop_annotation = data.frame(
      loop_type = c("E-P", "P-P"),
      Putative_Target_Genes = c("GENE_X;GENE_Y", "GENE_Z")
    )
  )

  res_targets <- extract_target_gene_sets(mock_anno, src = "targets", include_Filled = TRUE)
  expect_type(res_targets, "list")
  expect_true("Target_Genes" %in% names(res_targets))
  expect_setequal(res_targets$Target_Genes, c("GENE_A", "GENE_B", "GENE_C", "GENE_D"))

  res_nearest <- extract_target_gene_sets(mock_anno, src = "targets", use_nearest_gene = TRUE)
  expect_setequal(res_nearest$Target_Genes, c("GENE_A", "GENE_B"))

  res_loops <- extract_target_gene_sets(mock_anno, src = "loops", active_loop_types = c("E-P", "P-P"))
  expect_true("EP_Genes" %in% names(res_loops))
  expect_true("PP_Genes" %in% names(res_loops))
  expect_setequal(res_loops$EP_Genes, c("GENE_X", "GENE_Y"))

  bad_anno <- list(target_annotation = data.frame(wrong_col = 1:3))
  expect_error(extract_target_gene_sets(bad_anno, src = "targets"), "Required column")
})
