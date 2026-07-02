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

  res_targets <- looplook:::extract_target_gene_sets(mock_anno, src = "targets", include_Filled = TRUE)
  expect_type(res_targets, "list")
  expect_true("Target_Genes" %in% names(res_targets))
  expect_setequal(res_targets$Target_Genes, c("GENE_A", "GENE_B", "GENE_C", "GENE_D"))

  res_nearest <- looplook:::extract_target_gene_sets(mock_anno, src = "targets", use_nearest_gene = TRUE)
  expect_setequal(res_nearest$Target_Genes, c("GENE_A", "GENE_B"))

  res_loops <- looplook:::extract_target_gene_sets(mock_anno, src = "loops", active_loop_types = c("E-P", "P-P"))
  expect_true("EP_Genes" %in% names(res_loops))
  expect_true("PP_Genes" %in% names(res_loops))
  expect_setequal(res_loops$EP_Genes, c("GENE_X", "GENE_Y"))

  bad_anno <- list(target_annotation = data.frame(wrong_col = 1:3))
  expect_error(looplook:::extract_target_gene_sets(bad_anno, src = "targets"), "Required column")
})

test_that("extract_target_gene_sets 'active' mode reads Active_Target_Genes on loops branch", {
  # Active_Target_Genes is the strict whitelist-filtered column produced by
  # refine_loop_anchors_by_expression (no loop-anchor fallback).
  mock_refined <- list(
    loop_annotation = data.frame(
      loop_type = c("E-P", "P-P"),
      Putative_Target_Genes = c("GENE_X;GENE_Y;SILENT_GENE", "GENE_Z;SILENT_GENE2"),
      Active_Target_Genes    = c("GENE_X;GENE_Y",                       "GENE_Z"),
      stringsAsFactors = FALSE
    )
  )

  # active mode: only Active_Target_Genes (no SILENT_* fallback)
  res_active <- looplook:::extract_target_gene_sets(
    mock_refined, src = "loops",
    active_loop_types = c("E-P", "P-P"),
    target_mapping_mode = "active"
  )
  expect_setequal(res_active$EP_Genes, c("GENE_X", "GENE_Y"))
  expect_setequal(res_active$PP_Genes, c("GENE_Z"))

  # default mode: Putative_Target_Genes includes the silent fallback genes
  res_putative <- looplook:::extract_target_gene_sets(
    mock_refined, src = "loops",
    active_loop_types = c("E-P", "P-P"),
    target_mapping_mode = "all"
  )
  expect_true("SILENT_GENE" %in% res_putative$EP_Genes)
  expect_true("SILENT_GENE2" %in% res_putative$PP_Genes)

  # 'active' on the targets branch must error and point to the loops branch,
  # because Active_Target_Genes lives on loop_annotation, not bed_info.
  expect_error(
    looplook:::extract_target_gene_sets(
      list(target_annotation = data.frame(Assigned_Target_Genes_Filled = "G1")),
      src = "targets", target_mapping_mode = "active"
    ),
    "target_source = 'loops'"
  )

  # 'active' without the Active_Target_Genes column must error with a helpful
  # hint to run refine_loop_anchors_by_expression first.
  mock_unrefined <- list(
    loop_annotation = data.frame(
      loop_type = c("E-P"),
      Putative_Target_Genes = "G1",
      stringsAsFactors = FALSE
    )
  )
  expect_error(
    looplook:::extract_target_gene_sets(
      mock_unrefined, src = "loops", target_mapping_mode = "active"
    ),
    "refine_loop_anchors_by_expression"
  )
})
