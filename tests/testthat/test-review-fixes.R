# tests/testthat/test-review-fixes.R
# Synthetic tests for P0/P1 fixes (2026-07-14)
# Dependencies: TxDb/OrgDb installed packages (skip_if_not_installed per-test)

# ── 1. Expression: case-insensitive matching (P1-3) ──────────────────────────

test_that("load_expression_matrix normalises gene IDs to uppercase", {
    tmp <- tempfile(fileext = ".csv")
    writeLines("Gene\tS1\tS2\ntp53\t10\t20\nBRCA1\t30\t40", tmp)
    vals <- looplook:::load_expression_matrix(tmp, c("S1", "S2"))
    expect_true("TP53" %in% names(vals))
    expect_true("BRCA1" %in% names(vals))
    expect_false("tp53" %in% names(vals))
})

test_that(".get_expr uses toupper lookup", {
    vals <- c(TP53 = 10, BRCA1 = 30)
    expect_equal(unname(looplook:::.get_expr("tp53", vals)), 10)
    expect_true(is.na(looplook:::.get_expr("NONEXIST", vals)))
})


# ── 2. Expression: duplicate gene IDs deduplicated (P1-4) ────────────────────

test_that("load_expression_matrix deduplicates gene IDs", {
    tmp <- tempfile(fileext = ".csv")
    writeLines("Gene\tS1\nGENE1\t10\nGENE1\t20\nGENE2\t5", tmp)
    expect_message(
        vals <- looplook:::load_expression_matrix(tmp, "S1"), "Removed"
    )
    expect_equal(sum(duplicated(names(vals))), 0)
    expect_true("GENE1" %in% names(vals))
})


# ── 3. Expression: missing gene returns NA, not 0 (P1-2) ─────────────────────

test_that(".get_expr returns NA for unmatched genes, keeps true zero", {
    vals <- c(TP53 = 10, BRCA1 = 0)
    expect_true(is.na(looplook:::.get_expr("NONEXIST", vals)))
    expect_equal(unname(looplook:::.get_expr("BRCA1", vals)), 0)
})


# ── 4. Chromatin: recompute_targets = FALSE preserves objects (P0-2) ─────────

test_that("recompute_targets = FALSE preserves target_annotation", {
    skip("Requires full Bioconductor environment with compatible loop/BED data")
})


# ── 5. Chromatin: Is_Active_Gene inherited, not overwritten (P0-4) ───────────

test_that("chromatin refinement inherits Is_Active_Gene as Not_assessed", {
    skip("Requires full Bioconductor environment with compatible loop/BED data")
})


# ── 6. Chromatin: promoter stats tracking columns via .build_promoter_centric_df (P0-3) ──

test_that(".build_promoter_centric_df adds was/is tracking columns", {
    raw_stats <- data.frame(
        Gene = c("GENE1"),
        Total_Loops_Filtered = c(10),
        n_Linked_Promoters_Filtered = c(3),
        n_Linked_Distal_Filtered = c(7),
        Dominant_Interaction_Filtered = c("E-P"),
        stringsAsFactors = FALSE
    )
    upstream <- data.frame(
        Gene = "GENE1", Total_Loops = 10, n_Linked_Promoters = 3,
        n_Linked_Distal = 7, Dominant_Interaction = "E-P",
        Is_High_Connectivity_Gene = "Yes", Is_High_Distal_Connectivity_Gene = "Yes",
        Is_Active_Gene = "Yes", stringsAsFactors = FALSE
    )
    vals <- c(GENE1 = 100)
    res <- looplook:::.build_promoter_centric_df(
        raw_stats, upstream, vals,
        threshold = 1, hub_percentile = 0.95
    )
    expect_true("was_promoter_before" %in% colnames(res))
    expect_true("is_promoter_after" %in% colnames(res))
})


# ── 7. Seqlevel harmonization: GInteractions works (P0-6) ────────────────────

test_that("seqlevelsStyle harmonization on GInteractions", {
    skip_if_not_installed("InteractionSet")
    tmp <- tempfile(fileext = ".bedpe")
    writeLines("1\t100\t200\t1\t300\t400", tmp)
    gi <- suppressMessages(looplook:::bedpe_to_gi(tmp))
    GenomeInfoDb::seqlevelsStyle(gi) <- "UCSC"
    ref_gr <- GenomicRanges::GRanges("chr1:1-100")
    GenomeInfoDb::seqlevelsStyle(ref_gr) <- "NCBI"
    result <- looplook:::.harmonize_seqlevels(gi, ref_gr, "Test")
    expect_equal(GenomeInfoDb::seqlevelsStyle(result)[1], "NCBI")
})


# ── 8. Evidence levels use new labels (label rename) ─────────────────────────

test_that("enhancer_evidence uses canonical/strong/limited factor levels", {
    empty_anchors <- data.frame(
        anchor_id = character(), chr = character(), start = integer(),
        end = integer(), anchor_type = character(), anchor_gene = character(),
        cluster_id = character(), H3K4me1 = logical(), H3K27ac = logical(),
        ATAC = logical(), H3K27me3 = logical(), H3K4me3 = logical(),
        stringsAsFactors = FALSE
    )
    empty_matrix <- as.data.frame(matrix(
        nrow = 0, ncol = 5,
        dimnames = list(NULL, c("H3K4me1", "H3K27ac", "ATAC", "H3K27me3", "H3K4me3"))
    ))
    res <- looplook:::.assign_chromatin_confidence(
        empty_anchors, empty_matrix, character(0),
        c("H3K4me1", "H3K27ac", "ATAC", "H3K27me3", "H3K4me3")
    )
    expect_equal(
        levels(res$enhancer_evidence),
        c("canonical", "strong", "supported", "limited", "uncertain")
    )
})


# ── 9. .is_distal_like excludes eP/eG (P1-1) ─────────────────────────────────

test_that(".is_distal_like excludes eP and eG", {
    expect_true(looplook:::.is_distal_like("E"))
    expect_true(looplook:::.is_distal_like("dual"))
    expect_false(looplook:::.is_distal_like("eP"))
    expect_false(looplook:::.is_distal_like("eG"))
    expect_false(looplook:::.is_distal_like("P"))
    expect_false(looplook:::.is_distal_like("G"))
})


# ── 10. GO: enrichGO returns p.adjust column ─────────────────────────────────

test_that("enrichGO output has p.adjust column", {
    skip_if_not_installed("org.Hs.eg.db")
    skip_if_not_installed("clusterProfiler")
    suppressMessages({
        go <- clusterProfiler::enrichGO(
            gene = c("1", "2"), OrgDb = org.Hs.eg.db::org.Hs.eg.db,
            keyType = "ENTREZID", ont = "BP",
            pvalueCutoff = 1, qvalueCutoff = 1,
            minGSSize = 100, maxGSSize = 200
        )
    })
    # Even if empty, the column should exist
    expect_true("p.adjust" %in% colnames(as.data.frame(go)))
})
