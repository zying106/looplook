# tests/testthat/test-synthetic-invariants.R
# Minimal synthetic tests for critical invariants (2026-07-15)
# All data in-memory; no file I/O; targets <5s total runtime.
# Note: uses internal helper clean_gene_names — run via devtools::test() for namespace access

# Helper: bed_info where each input_id's gene is in Assigned_Target_Genes,
# so that .mark_target_gene_link_membership() flags them as retained.
.bed_for_links <- function(links) {
  data.frame(
    input_id = links$input_id,
    Assigned_Target_Genes = links$gene,
    stringsAsFactors = FALSE
  )
}


# ══════════════════════════════════════════════════════════════════════════
# 1. Expression state → Passes_Expression_Filter mapping
# ══════════════════════════════════════════════════════════════════════════

test_that("Expression_State: active→TRUE, measured_silent→FALSE, unmeasured→NA", {
  links <- data.frame(
    input_id = c("p1", "p2", "p3"),
    loop_ID = c("L1", "L2", "L3"),
    anchor_id = c("A1", "A2", "A3"),
    gene = c("GENE1", "GENE2", "GENE3"),
    gene_role = "promoter", source = "loop_anchor",
    evidence = "local_promoter_overlap", anchor_role = "proximal",
    path_length = 0L, Mean_Expression = c(50, 0.2, NA_real_),
    stringsAsFactors = FALSE
  )
  vals <- c(GENE1 = 50, GENE2 = 0.2)
  threshold <- 1

  res <- looplook:::.filter_refined_target_gene_links(
    links, .bed_for_links(links), vals, threshold
  )

  expect_equal(res$Expression_State[res$gene == "GENE1"], "active")
  expect_equal(res$Expression_State[res$gene == "GENE2"], "measured_silent")
  expect_equal(res$Expression_State[res$gene == "GENE3"], "unmeasured")

  expect_true(res$Passes_Expression_Filter[res$gene == "GENE1"])
  expect_false(res$Passes_Expression_Filter[res$gene == "GENE2"])
  expect_true(is.na(res$Passes_Expression_Filter[res$gene == "GENE3"]))
})


# ══════════════════════════════════════════════════════════════════════════
# 2. Has_Active_Target: unmeasured/mixed → NA, active → TRUE, silent → FALSE
# ══════════════════════════════════════════════════════════════════════════

test_that("Has_Active_Target: active→TRUE, measured_silent→FALSE, unmeasured→NA", {
  loop_df <- data.frame(
    anchor1_type = c("P", "P", "P"),
    anchor2_type = c("E", "E", "E"),
    anchor1_gene = c("G1", "G2", "G3"),
    anchor2_gene = c("", "", ""),
    Promoter_Target_Genes = c("G1", "", ""),
    stringsAsFactors = FALSE
  )
  whitelist <- c("G1")
  measured_set <- c("G1", "G2")
  orig_ptg <- c("G1", "G2", "G3")

  res <- looplook:::.refine_compute_targets(
    loop_df, orig_ptg, whitelist,
    orig_anchor1_type = loop_df$anchor1_type,
    orig_anchor2_type = loop_df$anchor2_type,
    measured_set = measured_set
  )

  expect_true(res$Has_Active_Target[1])
  expect_false(res$Has_Active_Target[2])
  expect_true(is.na(res$Has_Active_Target[3]))
})


# ══════════════════════════════════════════════════════════════════════════
# 3. Refinement_Action: unmeasured → expression_state_unknown
# ══════════════════════════════════════════════════════════════════════════

test_that("Refinement_Action: unmeasured→expression_state_unknown, active→retained", {
  loop_df <- data.frame(
    anchor1_type = c("P", "P"),
    anchor2_type = c("E", "E"),
    anchor1_gene = c("G1", "G3"),
    anchor2_gene = c("", ""),
    Promoter_Target_Genes = c("G1", ""),
    stringsAsFactors = FALSE
  )
  whitelist <- c("G1")
  measured_set <- c("G1", "G2")
  orig_ptg <- c("G1", "G3")

  res <- looplook:::.refine_compute_targets(
    loop_df, orig_ptg, whitelist,
    orig_anchor1_type = loop_df$anchor1_type,
    orig_anchor2_type = loop_df$anchor2_type,
    measured_set = measured_set
  )

  expect_equal(res$Refinement_Action[1], "retained_active_target")
  expect_equal(res$Refinement_Action[2], "expression_state_unknown")
})


# ══════════════════════════════════════════════════════════════════════════
# 4. Promoter_Target_Genes ⊆ Putative_Target_Genes invariant
# ══════════════════════════════════════════════════════════════════════════

test_that("Promoter ⊆ Putative: promoter targets subset of all targets", {
  loop_df <- data.frame(
    anchor1_type = c("P", "P", "P"),
    anchor2_type = c("E", "E", "E"),
    anchor1_gene = c("A;B", "C", "D;E"),
    anchor2_gene = c("", "", ""),
    Promoter_Target_Genes = c("A", "", "D"),
    stringsAsFactors = FALSE
  )
  wl <- toupper(c("A", "D"))
  ms <- toupper(c("A", "B", "C", "D", "E"))
  orig_ptg <- c("A;B", "C", "D;E")

  res <- looplook:::.refine_compute_targets(
    loop_df, orig_ptg, wl,
    orig_anchor1_type = loop_df$anchor1_type,
    orig_anchor2_type = loop_df$anchor2_type,
    measured_set = ms
  )

  for (i in seq_len(nrow(res))) {
    promoter_genes <- clean_gene_names(res$Promoter_Target_Genes[i])
    putative_genes <- clean_gene_names(res$Putative_Target_Genes[i])
    if (length(promoter_genes) > 0) {
      expect_true(all(promoter_genes %in% putative_genes),
        info = sprintf("Row %d: Promoter targets not in Putative", i)
      )
    }
  }
})


# ══════════════════════════════════════════════════════════════════════════
# 5. Retained_In_Functional_Network == Has_Active_Target
# ══════════════════════════════════════════════════════════════════════════

test_that("Retained_In_Functional_Network equals Has_Active_Target", {
  loop_df <- data.frame(
    anchor1_type = c("P", "P", "P"),
    anchor2_type = c("E", "E", "E"),
    anchor1_gene = c("G1", "G2", "G3"),
    anchor2_gene = c("", "", ""),
    Promoter_Target_Genes = c("G1", "", ""),
    stringsAsFactors = FALSE
  )
  wl <- c("G1")
  ms <- c("G1", "G2")
  orig_ptg <- c("G1", "G2", "G3")

  res <- looplook:::.refine_compute_targets(
    loop_df, orig_ptg, wl,
    orig_anchor1_type = loop_df$anchor1_type,
    orig_anchor2_type = loop_df$anchor2_type,
    measured_set = ms
  )

  expect_equal(res$Retained_In_Functional_Network, res$Has_Active_Target)
})


# ══════════════════════════════════════════════════════════════════════════
# 6. Expression_State in target links never empty
# ══════════════════════════════════════════════════════════════════════════

test_that("Expression_State in target links is never empty or NA", {
  links <- data.frame(
    input_id = c("p1", "p2"),
    loop_ID = c("L1", "L2"), anchor_id = c("A1", "A2"),
    gene = c("G1", "G2"), gene_role = "promoter",
    source = "loop_anchor", evidence = "local_promoter_overlap",
    anchor_role = "proximal", path_length = 0L,
    Mean_Expression = c(50, NA_real_), stringsAsFactors = FALSE
  )
  vals <- c(G1 = 50)
  threshold <- 1

  res <- looplook:::.filter_refined_target_gene_links(
    links, .bed_for_links(links), vals, threshold
  )

  expect_false(any(res$Expression_State == "" | is.na(res$Expression_State)))
})


# ══════════════════════════════════════════════════════════════════════════
# 7. Passes_Expression_Filter: unmeasured = NA, never FALSE
# ══════════════════════════════════════════════════════════════════════════

test_that("Passes_Expression_Filter: unmeasured is NA, not FALSE", {
  links <- data.frame(
    input_id = "p1", loop_ID = "L1", anchor_id = "A1",
    gene = "UNMEASURED", gene_role = "promoter",
    source = "loop_anchor", evidence = "local_promoter_overlap",
    anchor_role = "proximal", path_length = 0L,
    Mean_Expression = NA_real_, stringsAsFactors = FALSE
  )
  vals <- c(G1 = 50)
  threshold <- 1

  res <- looplook:::.filter_refined_target_gene_links(
    links, .bed_for_links(links), vals, threshold
  )

  expect_equal(res$Expression_State, "unmeasured")
  expect_true(is.na(res$Passes_Expression_Filter))
  expect_false(identical(res$Passes_Expression_Filter, FALSE))
})


# ══════════════════════════════════════════════════════════════════════════
# 8. Target_Expression_State: six loop-level classifications
# ══════════════════════════════════════════════════════════════════════════

test_that("Target_Expression_State: computed on pre-filter PTG", {
  ms <- c("A", "B", "C")
  wl <- c("A")

  state_for <- function(gene_string) {
    loop_df <- data.frame(
      anchor1_type = "P", anchor2_type = "E",
      anchor1_gene = "", anchor2_gene = "",
      Promoter_Target_Genes = "",
      stringsAsFactors = FALSE
    )
    # Target_Expression_State is now computed on the PRE-filter original_ptg
    res <- looplook:::.refine_compute_targets(
      loop_df, gene_string, wl,
      orig_anchor1_type = loop_df$anchor1_type,
      orig_anchor2_type = loop_df$anchor2_type,
      measured_set = ms
    )
    unname(res$Target_Expression_State)
  }

  expect_equal(state_for("A"), "active")
  expect_equal(state_for("A;D"), "active_partial")
  expect_equal(state_for("B"), "measured_silent")
  expect_equal(state_for("D"), "unmeasured")
  expect_equal(state_for("B;D"), "mixed")
  expect_equal(state_for(""), "no_target")
})


# ══════════════════════════════════════════════════════════════════════════
# 9. Refinement_Action: never empty or NA
# ══════════════════════════════════════════════════════════════════════════

test_that("Refinement_Action never empty or NA", {
  loop_df <- data.frame(
    anchor1_type = c("P", "E", "E", "P"),
    anchor2_type = c("E", "E", "E", "E"),
    anchor1_gene = c("G1", "", "", "G2"),
    anchor2_gene = c("", "", "", ""),
    Promoter_Target_Genes = c("G1", "", "", ""),
    stringsAsFactors = FALSE
  )
  wl <- c("G1")
  ms <- c("G1", "G2")
  orig_ptg <- c("G1", "", "G3", "G2")

  res <- looplook:::.refine_compute_targets(
    loop_df, orig_ptg, wl,
    orig_anchor1_type = loop_df$anchor1_type,
    orig_anchor2_type = loop_df$anchor2_type,
    measured_set = ms
  )

  actions <- res$Refinement_Action
  expect_false(any(is.na(actions)))
  expect_false(any(actions == ""))

  valid <- c(
    "retained_active_target", "expression_state_unknown",
    "expression_filtered_no_active_target", "structural_only_no_active_target",
    "reclassified_silent_anchor"
  )
  expect_true(all(actions %in% valid))
})


# ══════════════════════════════════════════════════════════════════════════
# 10. Has_Active_Target consistent with Target_Expression_State
# ══════════════════════════════════════════════════════════════════════════

test_that("Has_Active_Target consistent with Target_Expression_State", {
  loop_df <- data.frame(
    anchor1_type = c("P", "P", "P", "P"),
    anchor2_type = c("E", "E", "E", "E"),
    anchor1_gene = c("A", "B", "D", ""),
    anchor2_gene = c("", "", "", ""),
    Promoter_Target_Genes = c("A", "", "", ""),
    stringsAsFactors = FALSE
  )
  wl <- c("A")
  ms <- c("A", "B")
  orig_ptg <- c("A", "B", "D", "")

  res <- looplook:::.refine_compute_targets(
    loop_df, orig_ptg, wl,
    orig_anchor1_type = loop_df$anchor1_type,
    orig_anchor2_type = loop_df$anchor2_type,
    measured_set = ms
  )

  expect_equal(unname(res$Target_Expression_State[1]), "active")
  expect_true(res$Has_Active_Target[1])

  expect_equal(unname(res$Target_Expression_State[2]), "measured_silent")
  expect_false(res$Has_Active_Target[2])

  expect_equal(unname(res$Target_Expression_State[3]), "unmeasured")
  expect_true(is.na(res$Has_Active_Target[3]))

  expect_equal(unname(res$Target_Expression_State[4]), "no_target")
  expect_false(res$Has_Active_Target[4])
})


# ══════════════════════════════════════════════════════════════════════════
# 11. clean_gene_names: edge cases
# ══════════════════════════════════════════════════════════════════════════

test_that("clean_gene_names: empty, NA, whitespace, duplicates", {
  expect_equal(clean_gene_names("", ";"), character(0))
  expect_equal(clean_gene_names(NA_character_, ";"), character(0))
  expect_equal(clean_gene_names("  ", ";"), character(0))
  expect_equal(clean_gene_names("A;B;A", ";"), c("A", "B"))
  expect_equal(clean_gene_names("A; ;B", ";"), c("A", "B"))
})


# ══════════════════════════════════════════════════════════════════════════
# 12. target_gene_links schema includes Expression_State
# ══════════════════════════════════════════════════════════════════════════

test_that("target_gene_links schema includes Expression_State column", {
  empty <- looplook:::.empty_target_gene_links()
  expect_true("Expression_State" %in% colnames(empty))
  expect_true("Passes_Expression_Filter" %in% colnames(empty))
  expect_true("Mean_Expression" %in% colnames(empty))
  expect_type(empty$Expression_State, "character")
  expect_type(empty$Passes_Expression_Filter, "logical")
})


# ══════════════════════════════════════════════════════════════════════════
# 13. target_priority: promoter_then_distance vs distance_then_role
# ══════════════════════════════════════════════════════════════════════════

.make_mock_links <- function() {
  data.frame(
    input_id = c("p1", "p1", "p2", "p2"),
    gene = c("GENEA", "GENEB", "GENEC", "GENED"),
    gene_role = c("gene_body", "promoter", "promoter", "promoter"),
    evidence = c(
      "gene_body_context", "distal_promoter",
      "distal_promoter", "distal_promoter"
    ),
    source = "loop_anchor",
    path_length = c(0L, 1L, 1L, 1L),
    stringsAsFactors = FALSE
  )
}
.make_mock_bed <- function() {
  data.frame(input_id = c("p1", "p2"), stringsAsFactors = FALSE)
}

test_that("target_priority: promoter_then_distance prefers promoter over gene_body", {
  links <- .make_mock_links()
  bed <- .make_mock_bed()
  # p1 has gene_body path0 and promoter path1 -> promoter_then_distance picks promoter
  res <- looplook:::.aggregate_strict_targets(
    links, bed,
    target_priority = "promoter_then_distance"
  )
  res_p1 <- res$Assigned_Target_Genes[res$input_id == "p1"]
  # Should pick GENEB (promoter), not GENEA (gene_body)
  expect_true(grepl("GENEB", res_p1), info = "promoter_then_distance should prefer promoter over gene_body")
  expect_false(grepl("GENEA", res_p1), info = "promoter_then_distance should NOT pick gene_body over promoter")
})

test_that("target_priority: distance_then_role prefers gene_body path0 over promoter path1", {
  links <- .make_mock_links()
  bed <- .make_mock_bed()
  res <- looplook:::.aggregate_strict_targets(
    links, bed,
    target_priority = "distance_then_role"
  )
  res_p1 <- res$Assigned_Target_Genes[res$input_id == "p1"]
  expect_true(grepl("GENEA", res_p1),
    info = "distance_then_role should prefer gene_body path0 over promoter path1"
  )
})

test_that("target_priority: tied promoter genes are all retained", {
  links <- .make_mock_links()
  bed <- .make_mock_bed()
  # p2 has two promoter path1 genes -> both should be retained
  res <- looplook:::.aggregate_strict_targets(
    links, bed,
    target_priority = "promoter_then_distance"
  )
  res_p2 <- res$Assigned_Target_Genes[res$input_id == "p2"]
  expect_true(grepl("GENEC", res_p2) && grepl("GENED", res_p2),
    info = "tied promoter genes should both be retained (semicolon-delimited)"
  )
})


# ══════════════════════════════════════════════════════════════════════════
# 14. strict_assignment_eligible filtering in target aggregation
# ══════════════════════════════════════════════════════════════════════════

test_that("strict gene_body link remains assignment-eligible after chromatin refinement", {
  links <- data.frame(
    input_id = "p1", gene = "MYC",
    gene_role = "gene_body", source = "loop_anchor",
    evidence = "basic_gene_body_retained_after_chromatin",
    path_length = 1L, strict_assignment_eligible = TRUE,
    stringsAsFactors = FALSE
  )
  bed <- data.frame(input_id = "p1", stringsAsFactors = FALSE)
  res <- looplook:::.aggregate_strict_targets(links, bed, max_primary_hop = 1L)
  expect_equal(res$Assigned_Target_Genes, "MYC")
  expect_true(is.na(res$Expanded_Target_Genes) | res$Expanded_Target_Genes == "")
})

test_that("resolve_chromatin_gene_role maps final types correctly", {
  # Unified resolver: final_type + old_type + TSS determines role + strict.
  # Test final-type-only role mapping (all old_type="P", no TSS):
  r <- looplook:::.resolve_chromatin_gene_role(
    old_type   = c("P", "P", "P", "P", "P", "P", NA),
    final_type = c("P", "eP", "dual", "G", "eG", "E", NA),
    tss_supported = FALSE,
    has_gene   = TRUE
  )
  # NA old_type + NA final_type → role=NA, strict=NA (unresolved anchor)
  expect_equal(r$role,
    c("promoter", "promoter", "promoter", "gene_body", "gene_body",
      "enhancer_candidate", NA_character_))
  expect_equal(r$strict,
    c(TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, NA))
})

test_that("P0-2: E->P/dual without TSS is positional_candidate, not strict", {
  # E->P/dual anchors have promoter-like chromatin but the nearest gene
  # is not reliable.  Until TSS reannotation confirms gene identity,
  # the link stays positional_candidate and does NOT enter strict targets.
  links <- data.frame(
    input_id = "p1", gene = "SOX2",
    gene_role = "positional_candidate", source = "loop_anchor",
    evidence = "positional_candidate_after_chromatin",
    path_length = 1L, strict_assignment_eligible = FALSE,
    stringsAsFactors = FALSE
  )
  bed <- data.frame(input_id = "p1", stringsAsFactors = FALSE)
  res <- looplook:::.aggregate_strict_targets(links, bed, max_primary_hop = 1L)
  expect_true(is.na(res$Assigned_Target_Genes) | res$Assigned_Target_Genes == "")
  expect_true(is.na(res$Regulated_promoter_genes) |
    res$Regulated_promoter_genes == "")
})

test_that("P0-1: E->E unchanged nearest gene is NOT strict-eligible", {
  # Unchanged E anchors carry only ChIPseeker nearest genes.
  # These must not enter strict target columns even when the E anchor
  # participates in a loop (source == "loop_anchor").
  links <- data.frame(
    input_id = "p1", gene = "SOX2",
    gene_role = "enhancer_candidate", source = "loop_anchor",
    evidence = "local_enhancer_candidate",
    path_length = 0L, strict_assignment_eligible = FALSE,
    stringsAsFactors = FALSE
  )
  bed <- data.frame(input_id = "p1", stringsAsFactors = FALSE)
  res <- looplook:::.aggregate_strict_targets(links, bed, max_primary_hop = 1L)
  # All_Loop_Connected_Genes includes it (structural connectivity)
  expect_equal(res$All_Loop_Connected_Genes, "SOX2")
  # But strict columns exclude it (no structural gene provenance)
  expect_true(is.na(res$Assigned_Target_Genes) | res$Assigned_Target_Genes == "")
})

test_that("G->P/dual host gene enters strict target as promoter (rank 1)", {
  # G/eG anchors have structural gene-body overlap.  When chromatin promotes
  # them to P/dual, the host gene is a valid promoter-associated target.
  links <- data.frame(
    input_id = "p1", gene = "MYC",
    gene_role = "promoter", source = "loop_anchor",
    evidence = "distal_promoter",
    path_length = 1L, strict_assignment_eligible = TRUE,
    stringsAsFactors = FALSE
  )
  bed <- data.frame(input_id = "p1", stringsAsFactors = FALSE)
  res <- looplook:::.aggregate_strict_targets(links, bed, max_primary_hop = 1L)
  expect_equal(res$Assigned_Target_Genes, "MYC")
  expect_equal(res$Regulated_promoter_genes, "MYC")
})

test_that("P/eP->E structurally supported gene enters as enhancer_candidate (rank 3)", {
  # P/eP anchors have annotated promoter genes.  When chromatin downgrades
  # them to E, the gene is still structurally valid but at rank 3.
  links <- data.frame(
    input_id = "p1", gene = "TP53",
    gene_role = "enhancer_candidate", source = "loop_anchor",
    evidence = "distal_enhancer_candidate",
    path_length = 1L, strict_assignment_eligible = TRUE,
    stringsAsFactors = FALSE
  )
  bed <- data.frame(input_id = "p1", stringsAsFactors = FALSE)
  res <- looplook:::.aggregate_strict_targets(links, bed, max_primary_hop = 1L)
  expect_equal(res$Assigned_Target_Genes, "TP53")
  # Must NOT be in Regulated_promoter_genes (no longer promoter role)
  expect_true(is.na(res$Regulated_promoter_genes) |
    res$Regulated_promoter_genes == "")
})

test_that("enhancer_candidate: rank 3, loses to gene_body at rank 2", {
  links <- data.frame(
    input_id = c("p1", "p1"),
    gene = c("SOX2", "MYC"),
    gene_role = c("enhancer_candidate", "gene_body"),
    source = "loop_anchor",
    evidence = c("distal_enhancer_candidate", "gene_body_context"),
    path_length = c(1L, 1L),
    strict_assignment_eligible = TRUE,
    stringsAsFactors = FALSE
  )
  bed <- data.frame(input_id = "p1", stringsAsFactors = FALSE)
  res <- looplook:::.aggregate_strict_targets(links, bed, max_primary_hop = 1L)
  # gene_body (rank 2) beats enhancer_candidate (rank 3)
  expect_equal(res$Assigned_Target_Genes, "MYC")
  expect_true(grepl("SOX2", res$All_Loop_Connected_Genes))
})

test_that("positional_candidate excluded from strict target", {
  # Positional candidates are never strict-eligible; always excluded.
  links <- data.frame(
    input_id = "p1", gene = "SOX2",
    gene_role = "positional_candidate", source = "loop_anchor",
    evidence = "positional_candidate_after_chromatin",
    path_length = 1L, strict_assignment_eligible = FALSE,
    stringsAsFactors = FALSE
  )
  bed <- data.frame(input_id = "p1", stringsAsFactors = FALSE)
  res <- looplook:::.aggregate_strict_targets(links, bed, max_primary_hop = 1L)
  expect_true(is.na(res$Assigned_Target_Genes) | res$Assigned_Target_Genes == "")
})

test_that("TSS-supported: promoter role enters strict target", {
  links <- data.frame(
    input_id = "p1", gene = "EGFR",
    gene_role = "promoter", source = "loop_anchor",
    evidence = "txdb_promoter_overlap",
    path_length = 1L, strict_assignment_eligible = TRUE,
    stringsAsFactors = FALSE
  )
  bed <- data.frame(input_id = "p1", stringsAsFactors = FALSE)
  res <- looplook:::.aggregate_strict_targets(links, bed, max_primary_hop = 1L)
  expect_equal(res$Assigned_Target_Genes, "EGFR")
})

test_that("P/eP anchor: promoter role preserved", {
  links <- data.frame(
    input_id = "p1", gene = "TP53",
    gene_role = "promoter", source = "loop_anchor",
    evidence = "local_promoter_overlap",
    path_length = 0L, strict_assignment_eligible = TRUE,
    stringsAsFactors = FALSE
  )
  bed <- data.frame(input_id = "p1", stringsAsFactors = FALSE)
  res <- looplook:::.aggregate_strict_targets(links, bed, max_primary_hop = 1L)
  expect_equal(res$Assigned_Target_Genes, "TP53")
})

test_that("positional_candidate is in all-connected but excluded from strict and expanded", {
  links <- data.frame(
    input_id = "p1", gene = "SOX2",
    gene_role = "positional_candidate", source = "loop_anchor",
    evidence = "positional_candidate_after_chromatin",
    path_length = 2L, strict_assignment_eligible = FALSE,
    stringsAsFactors = FALSE
  )
  bed <- data.frame(input_id = "p1", stringsAsFactors = FALSE)
  res <- looplook:::.aggregate_strict_targets(links, bed, max_primary_hop = 1L)
  expect_equal(res$All_Loop_Connected_Genes, "SOX2")
  expect_true(is.na(res$Assigned_Target_Genes) | res$Assigned_Target_Genes == "")
  expect_true(is.na(res$Expanded_Target_Genes) | res$Expanded_Target_Genes == "")
})

test_that("gene_body role enters strict target at rank 2", {
  links <- data.frame(
    input_id = "p1", gene = "MYC",
    gene_role = "gene_body", source = "loop_anchor",
    evidence = "basic_annotation_gene_body",
    path_length = 1L, strict_assignment_eligible = TRUE,
    stringsAsFactors = FALSE
  )
  bed <- data.frame(input_id = "p1", stringsAsFactors = FALSE)
  res <- looplook:::.aggregate_strict_targets(links, bed, max_primary_hop = 1L)
  expect_equal(res$Assigned_Target_Genes, "MYC")
})

test_that("promoter role preserved as rank 1 in aggregation", {
  links <- data.frame(
    input_id = "p1", gene = "TP53",
    gene_role = "promoter", source = "loop_anchor",
    evidence = "basic_annotation_promoter",
    path_length = 1L, strict_assignment_eligible = TRUE,
    stringsAsFactors = FALSE
  )
  bed <- data.frame(input_id = "p1", stringsAsFactors = FALSE)
  res <- looplook:::.aggregate_strict_targets(links, bed, max_primary_hop = 1L)
  expect_equal(res$Assigned_Target_Genes, "TP53")
})

test_that("positional_candidate always excluded regardless of provenance flags", {
  links <- data.frame(
    input_id = "p1", gene = "SOX2",
    gene_role = "positional_candidate", source = "loop_anchor",
    evidence = "basic_annotation_positional_only",
    path_length = 1L, strict_assignment_eligible = FALSE,
    stringsAsFactors = FALSE
  )
  bed <- data.frame(input_id = "p1", stringsAsFactors = FALSE)
  res <- looplook:::.aggregate_strict_targets(links, bed, max_primary_hop = 1L)
  expect_true(is.na(res$Assigned_Target_Genes) | res$Assigned_Target_Genes == "")
})

test_that("per-loop gene summary not globally merged", {
  links <- data.frame(
    input_id = c("p1", "p2"),
    gene = c("TP53", "EGFR"),
    gene_role = c("promoter", "promoter"),
    source = "loop_anchor",
    evidence = "distal_promoter",
    path_length = c(1L, 1L),
    strict_assignment_eligible = TRUE,
    stringsAsFactors = FALSE
  )
  bed <- data.frame(input_id = c("p1", "p2"), stringsAsFactors = FALSE)
  res <- looplook:::.aggregate_strict_targets(links, bed, max_primary_hop = 1L)
  expect_equal(res$Assigned_Target_Genes[res$input_id == "p1"], "TP53")
  expect_equal(res$Assigned_Target_Genes[res$input_id == "p2"], "EGFR")
  expect_false(grepl("TP53;EGFR", res$Assigned_Target_Genes[1]))
})

test_that("same gene multi-evidence: one excluded link does not remove other", {
  links <- data.frame(
    input_id = c("p1", "p1"),
    gene = c("MYC", "MYC"),
    gene_role = c("positional_candidate", "promoter"),
    source = "loop_anchor",
    evidence = c("positional_candidate_after_chromatin", "txdb_promoter_overlap"),
    path_length = c(1L, 1L),
    strict_assignment_eligible = c(FALSE, TRUE),
    stringsAsFactors = FALSE
  )
  bed <- data.frame(input_id = "p1", stringsAsFactors = FALSE)
  res <- looplook:::.aggregate_strict_targets(links, bed, max_primary_hop = 1L)
  expect_equal(res$Assigned_Target_Genes, "MYC")
})


# ══════════════════════════════════════════════════════════════════════════
# 15. Old schema positional candidate exclusion (looplook28)
# ══════════════════════════════════════════════════════════════════════════

test_that("old schema positional candidate excluded from strict (no eligibility col)", {
  links <- data.frame(
    input_id = "p1",
    gene = "SOX2",
    gene_role = "positional_candidate",
    source = "loop_anchor",
    evidence = "positional_candidate_after_chromatin",
    path_length = 1L,
    stringsAsFactors = FALSE
  )
  bed <- data.frame(input_id = "p1", stringsAsFactors = FALSE)
  res <- looplook:::.aggregate_strict_targets(links, bed, max_primary_hop = 1L)
  expect_equal(res$All_Loop_Connected_Genes, "SOX2")
  expect_equal(res$Candidate_Positional_Genes, "SOX2")
  expect_true(is.na(res$Assigned_Target_Genes) || res$Assigned_Target_Genes == "")
  expect_true(is.na(res$Expanded_Target_Genes) || res$Expanded_Target_Genes == "")
})

test_that("old schema positional candidate excluded from expanded (path > 1, no eligibility col)", {
  links <- data.frame(
    input_id = "p1",
    gene = "SOX2",
    gene_role = "positional_candidate",
    source = "loop_anchor",
    evidence = "positional_candidate_after_chromatin",
    path_length = 2L,
    stringsAsFactors = FALSE
  )
  bed <- data.frame(input_id = "p1", stringsAsFactors = FALSE)
  res <- looplook:::.aggregate_strict_targets(links, bed, max_primary_hop = 1L)
  expect_equal(res$Candidate_Positional_Genes, "SOX2")
  expect_true(is.na(res$Assigned_Target_Genes) || res$Assigned_Target_Genes == "")
  expect_true(is.na(res$Expanded_Target_Genes) || res$Expanded_Target_Genes == "")
})

test_that("positional candidate is complement of strict: anti-join prevents overlap", {
  links <- data.frame(
    input_id = c("p1", "p1"),
    gene = c("SOX2", "SOX2"),
    gene_role = c("positional_candidate", "promoter"),
    source = "loop_anchor",
    evidence = c("positional_candidate_after_chromatin", "local_promoter_overlap"),
    path_length = c(1L, 1L),
    strict_assignment_eligible = c(FALSE, TRUE),
    stringsAsFactors = FALSE
  )
  bed <- data.frame(input_id = "p1", stringsAsFactors = FALSE)
  res <- looplook:::.aggregate_strict_targets(links, bed, max_primary_hop = 1L)
  expect_equal(res$Assigned_Target_Genes, "SOX2",
    info = "Gene with both positional and strict evidence enters Assigned"
  )
  expect_true(is.na(res$Candidate_Positional_Genes) || res$Candidate_Positional_Genes == "",
    info = "Candidate must not include genes already in strict set"
  )
})

test_that("NA gene_role does not crash positional eligibility enforcement", {
  links <- data.frame(
    input_id = "p1",
    gene = "GENE1",
    gene_role = NA_character_,
    source = "loop_anchor",
    evidence = "unknown",
    path_length = 1L,
    strict_assignment_eligible = TRUE,
    stringsAsFactors = FALSE
  )
  bed <- data.frame(input_id = "p1", stringsAsFactors = FALSE)
  expect_no_error(
    looplook:::.aggregate_strict_targets(links, bed, max_primary_hop = 1L)
  )
})

# ══════════════════════════════════════════════════════════════════════════
# 16. target_gene_links invariant: positional_candidate never strict-eligible
# ══════════════════════════════════════════════════════════════════════════

test_that("invariant: no positional_candidate row has strict_assignment_eligible = TRUE", {
  links <- data.frame(
    input_id = "p1", gene = "SOX2",
    gene_role = "positional_candidate",
    source = "loop_anchor",
    evidence = "positional_candidate_after_chromatin",
    path_length = 1L,
    strict_assignment_eligible = FALSE,
    stringsAsFactors = FALSE
  )
  bed <- data.frame(input_id = "p1", stringsAsFactors = FALSE)
  res <- looplook:::.aggregate_strict_targets(links, bed, max_primary_hop = 1L)
  expect_true(is.na(res$Assigned_Target_Genes) | res$Assigned_Target_Genes == "")
  expect_equal(res$Candidate_Positional_Genes, "SOX2")
})

test_that("invariant: old schema (no eligibility col) also prevents positional from strict", {
  links <- data.frame(
    input_id = "p1", gene = "SOX2",
    gene_role = "positional_candidate",
    source = "loop_anchor",
    evidence = "positional_candidate_after_chromatin",
    path_length = 1L,
    stringsAsFactors = FALSE
  )
  bed <- data.frame(input_id = "p1", stringsAsFactors = FALSE)
  res <- looplook:::.aggregate_strict_targets(links, bed, max_primary_hop = 1L)
  expect_true(is.na(res$Assigned_Target_Genes) | res$Assigned_Target_Genes == "")
  expect_equal(res$Candidate_Positional_Genes, "SOX2")
})

test_that("invariant: chromatin_annotate_links enforces positional_candidate ineligibility", {
  new_links <- data.frame(
    input_id = "p1", loop_ID = "L1", anchor_id = "A1",
    gene = "SOX2",
    gene_role = "positional_candidate",
    source = "loop_anchor",
    evidence = "positional_candidate_after_chromatin",
    anchor_role = "distal", path_length = 1L,
    anchor_type_before_chromatin = "E",
    anchor_type_after_chromatin = "P",
    TSS_supported = FALSE,
    stringsAsFactors = FALSE
  )
  reclass_map <- data.frame(
    anchor_id = "A1",
    old_type = "E", new_type = "P",
    changed = TRUE,
    stringsAsFactors = FALSE
  )
  annotation_res <- list(chromatin_validation = NULL)
  res <- looplook:::.chromatin_annotate_links(new_links, reclass_map, annotation_res)
  expect_false(any(res$gene_role == "positional_candidate" & res$strict_assignment_eligible),
    info = "positional_candidate rows must never have strict_assignment_eligible = TRUE"
  )
})

test_that("invariant: old schema input to chromatin_annotate_links (no eligibility col)", {
  new_links <- data.frame(
    input_id = "p1", loop_ID = "L1", anchor_id = "A1",
    gene = "SOX2",
    gene_role = "positional_candidate",
    source = "loop_anchor",
    evidence = "positional_candidate_after_chromatin",
    anchor_role = "distal", path_length = 1L,
    anchor_type_before_chromatin = "E",
    anchor_type_after_chromatin = "P",
    TSS_supported = FALSE,
    stringsAsFactors = FALSE
  )
  reclass_map <- data.frame(
    anchor_id = "A1",
    old_type = "E", new_type = "P",
    changed = TRUE,
    stringsAsFactors = FALSE
  )
  annotation_res <- list(chromatin_validation = NULL)
  res <- looplook:::.chromatin_annotate_links(new_links, reclass_map, annotation_res)
  expect_false(any(res$gene_role == "positional_candidate" & res$strict_assignment_eligible),
    info = "old schema: positional_candidate must still get strict_assignment_eligible = FALSE"
  )
})

test_that("invariant: positional_candidate with NA eligibility also forced to FALSE", {
  links <- data.frame(
    input_id = "p1", gene = "SOX2",
    gene_role = "positional_candidate",
    source = "loop_anchor",
    evidence = "positional_candidate_after_chromatin",
    path_length = 1L,
    strict_assignment_eligible = NA,
    stringsAsFactors = FALSE
  )
  bed <- data.frame(input_id = "p1", stringsAsFactors = FALSE)
  res <- looplook:::.aggregate_strict_targets(links, bed, max_primary_hop = 1L)
  expect_true(is.na(res$Assigned_Target_Genes) | res$Assigned_Target_Genes == "")
  expect_equal(res$Candidate_Positional_Genes, "SOX2")
})
