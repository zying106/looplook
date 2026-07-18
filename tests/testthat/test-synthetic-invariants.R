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
        gene     = c("GENEA", "GENEB", "GENEC", "GENED"),
        gene_role = c("gene_body", "promoter", "promoter", "promoter"),
        evidence = c("gene_body_context", "distal_promoter",
                     "distal_promoter", "distal_promoter"),
        source   = "loop_anchor",
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
        links, bed, target_priority = "promoter_then_distance")
    res_p1 <- res$Assigned_Target_Genes[res$input_id == "p1"]
    # Should pick GENEB (promoter), not GENEA (gene_body)
    expect_true(grepl("GENEB", res_p1), info = "promoter_then_distance should prefer promoter over gene_body")
    expect_false(grepl("GENEA", res_p1), info = "promoter_then_distance should NOT pick gene_body over promoter")
})

test_that("target_priority: distance_then_role prefers gene_body path0 over promoter path1", {
    links <- .make_mock_links()
    bed <- .make_mock_bed()
    res <- looplook:::.aggregate_strict_targets(
        links, bed, target_priority = "distance_then_role")
    res_p1 <- res$Assigned_Target_Genes[res$input_id == "p1"]
    expect_true(grepl("GENEA", res_p1),
        info = "distance_then_role should prefer gene_body path0 over promoter path1")
})

test_that("target_priority: tied promoter genes are all retained", {
    links <- .make_mock_links()
    bed <- .make_mock_bed()
    # p2 has two promoter path1 genes -> both should be retained
    res <- looplook:::.aggregate_strict_targets(
        links, bed, target_priority = "promoter_then_distance")
    res_p2 <- res$Assigned_Target_Genes[res$input_id == "p2"]
    expect_true(grepl("GENEC", res_p2) && grepl("GENED", res_p2),
        info = "tied promoter genes should both be retained (semicolon-delimited)")
})
