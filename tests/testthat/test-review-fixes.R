# tests/testthat/test-review-fixes.R
# Synthetic tests for P0/P1 fixes (2026-07-14)
# Dependencies: TxDb/OrgDb installed packages (skip_if_not_installed per-test)

# Shared annotation bases are provided by helper-fixtures.R (built once per
# test session). `.get_shared_*` aliases are kept for readability and map to
# the session-wide fixtures.
.get_shared_target_base <- function() get_base_annotation_target()
.get_shared_hop1_base <- function() get_base_annotation_hop1()

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


# ── 2. Expression: duplicate gene IDs rejected (P1-4) ───────────────────────

test_that("load_expression_matrix rejects duplicate gene IDs", {
  tmp <- tempfile(fileext = ".csv")
  writeLines("Gene\tS1\nGENE1\t10\nGENE1\t20\nGENE2\t5", tmp)
  expect_error(
    looplook:::load_expression_matrix(tmp, "S1"),
    "duplicated gene identifier"
  )
  unlink(tmp)
})


# ── 3. Expression: missing gene returns NA, not 0 (P1-2) ─────────────────────

test_that(".get_expr returns NA for unmatched genes, keeps true zero", {
  vals <- c(TP53 = 10, BRCA1 = 0)
  expect_true(is.na(looplook:::.get_expr("NONEXIST", vals)))
  expect_equal(unname(looplook:::.get_expr("BRCA1", vals)), 0)
})


# ── 4. Chromatin: recompute_targets = FALSE preserves objects (P0-2) ─────────

test_that("recompute_targets = FALSE preserves target_annotation", {
  skip_if_not_installed("org.Hs.eg.db")
  skip_if(is.null(get_test_txdb()), "Sample TxDb not available")
  h3k4me1 <- tempfile(fileext = ".bed")
  h3k4me3 <- tempfile(fileext = ".bed")
  writeLines("chr6\t10411900\t10412700", h3k4me1)
  writeLines("chr6\t10411900\t10412700", h3k4me3)
  base_res <- .get_shared_target_base()
  before_ta <- base_res$target_annotation
  before_tgl <- base_res$target_gene_links
  cr <- refine_loop_anchors_by_chromatin(base_res,
    chromatin_beds = list(H3K4me1 = h3k4me1, H3K4me3 = h3k4me3),
    recompute_targets = FALSE, write_output = FALSE, quiet = TRUE
  )
  expect_false(is.null(cr$target_gene_links))
  expect_equal(cr$target_annotation$Assigned_Target_Genes, before_ta$Assigned_Target_Genes)
  unlink(c(h3k4me1, h3k4me3))
})

test_that("chromatin refinement inherits Is_Active_Gene as Not_assessed", {
  skip_if_not_installed("org.Hs.eg.db")
  rdata_path <- system.file("extdata", "analysis_results.RData", package = "looplook")
  skip_if(rdata_path == "", "Pre-computed test data not available")
  tmp <- new.env()
  load(rdata_path, envir = tmp)
  base_res <- tmp[[ls(tmp)[1]]]
  # NOTE: shrink_annotation_res() was never implemented; the full
  # analysis_results.RData object is small enough to use directly.
  h3k4me1 <- tempfile(fileext = ".bed")
  h3k4me3 <- tempfile(fileext = ".bed")
  # Cover the chr1 anchor region used in analysis_results.RData
  writeLines("chr1\t109492600\t151588200", h3k4me1)
  writeLines("chr1\t109492600\t151588200", h3k4me3)
  cr <- refine_loop_anchors_by_chromatin(base_res,
    chromatin_beds = list(H3K4me1 = h3k4me1, H3K4me3 = h3k4me3),
    recompute_targets = TRUE, write_output = FALSE, quiet = TRUE
  )
  expect_true("Is_Active_Gene" %in% colnames(cr$promoter_centric_stats))
  iag <- cr$promoter_centric_stats$Is_Active_Gene
  expect_true(any(!is.na(iag)), info = "At least one gene must have an activity state")
  expect_setequal(unique(iag[!is.na(iag)]), "Not_assessed")
  unlink(c(h3k4me1, h3k4me3))
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

# ── 10. Annotation metadata: neighbor_hop recording ───────────────────────────

test_that("basic annotation records neighbor_hop=0 in metadata and state", {
  skip_if_not_installed("org.Hs.eg.db")
  skip_if(is.null(get_test_txdb()), "Sample TxDb not available")
  txdb_obj <- get_test_txdb()
  tmp_bedpe <- tempfile(fileext = ".bedpe")
  writeLines("chr6\t10412000\t10412600\tchr6\t10415000\t10415600", tmp_bedpe)
  res <- annotate_peaks_and_loops(bedpe_file = tmp_bedpe,
    txdb = txdb_obj, org_db = "org.Hs.eg.db", species = "hg19",
    neighbor_hop = 0, out_dir = tempdir(), write_output = FALSE, quiet = TRUE)
  expect_equal(res$metadata$parameters$neighbor_hop, 0L)
  expect_equal(attr(res, "looplook_anchor_state")$neighbor_hop, 0L)
  unlink(tmp_bedpe)
})


# ══════════════════════════════════════════════════════════════════════════
# 11. allow_rerefine: chromatin re-refinement stop/allow (P0-2)
# ══════════════════════════════════════════════════════════════════════════

test_that("refine_loop_anchors_by_chromatin stops on already chromatin-refined object", {
  skip_if_not_installed("org.Hs.eg.db")
  skip_if(is.null(get_test_txdb()), "Sample TxDb not available")
  h3k4me1 <- tempfile(fileext = ".bed")
  h3k4me3 <- tempfile(fileext = ".bed")
  writeLines("chr6\t10411900\t10412700", h3k4me1)
  writeLines("chr6\t10411900\t10412700", h3k4me3)
  base_res <- .get_shared_target_base()
  cr <- refine_loop_anchors_by_chromatin(base_res,
    chromatin_beds = list(H3K4me1 = h3k4me1, H3K4me3 = h3k4me3),
    recompute_targets = FALSE, write_output = FALSE, quiet = TRUE
  )
  expect_error(
    refine_loop_anchors_by_chromatin(cr,
      chromatin_beds = list(H3K4me1 = h3k4me1, H3K4me3 = h3k4me3),
      recompute_targets = FALSE, write_output = FALSE, quiet = TRUE
    ),
    "already chromatin-refined"
  )
  unlink(c(h3k4me1, h3k4me3))
})

test_that("refine_loop_anchors_by_chromatin allow_rerefine=TRUE proceeds", {
  skip_if_not_installed("org.Hs.eg.db")
  skip_if(is.null(get_test_txdb()), "Sample TxDb not available")
  h3k4me1 <- tempfile(fileext = ".bed")
  h3k4me3 <- tempfile(fileext = ".bed")
  writeLines("chr6\t10411900\t10412700", h3k4me1)
  writeLines("chr6\t10411900\t10412700", h3k4me3)
  base_res <- .get_shared_target_base()
  cr <- refine_loop_anchors_by_chromatin(base_res,
    chromatin_beds = list(H3K4me1 = h3k4me1, H3K4me3 = h3k4me3),
    recompute_targets = FALSE, write_output = FALSE, quiet = TRUE
  )
  expect_no_error(
    refine_loop_anchors_by_chromatin(cr,
      chromatin_beds = list(H3K4me1 = h3k4me1, H3K4me3 = h3k4me3),
      recompute_targets = FALSE, write_output = FALSE, quiet = TRUE,
      allow_rerefine = TRUE
    )
  )
  unlink(c(h3k4me1, h3k4me3))
})


# ══════════════════════════════════════════════════════════════════════════
# 12. gene_assignment_policy: stored in anchor_state and preserved through
#    expression refinement (conflict_min_expr vs effective_threshold)
# ══════════════════════════════════════════════════════════════════════════

test_that("gene_assignment_policy stored in anchor_state after annotation", {
  skip_if_not_installed("org.Hs.eg.db")
  skip_if(is.null(get_test_txdb()), "Sample TxDb not available")
  txdb_obj <- get_test_txdb()
  tmp_bedpe <- tempfile(fileext = ".bedpe")
  writeLines("chr6\t10412000\t10412600\tchr6\t10415000\t10415600", tmp_bedpe)
  base_res <- annotate_peaks_and_loops(
    bedpe_file = tmp_bedpe,
    txdb = txdb_obj, org_db = "org.Hs.eg.db", species = "hg19",
    min_expr = 0.5,
    conflict_strategy = "expression_first",
    co_dominance_ratio = 0.2,
    out_dir = tempdir(), write_output = FALSE, quiet = TRUE
  )
  state <- attr(base_res, "looplook_anchor_state")
  policy <- state$gene_assignment_policy
  expect_equal(policy$conflict_min_expr, 0.5)
  expect_equal(policy$conflict_strategy, "expression_first")
  expect_equal(policy$co_dominance_ratio, 0.2)
  expect_true(!is.null(policy$biotype_order))
  unlink(tmp_bedpe)
})

test_that("conflict_min_expr preserved after expression refinement (not overwritten)", {
  skip_if_not_installed("org.Hs.eg.db")
  skip_if(is.null(get_test_txdb()), "Sample TxDb not available")
  expr_path <- system.file("extdata", "example_tpm.txt", package = "looplook")
  skip_if(expr_path == "", "Example TPM file not available")
  txdb_obj <- get_test_txdb()
  tmp_bedpe <- tempfile(fileext = ".bedpe")
  writeLines("chr6\t10412000\t10412600\tchr6\t10415000\t10415600", tmp_bedpe)
  base_res <- annotate_peaks_and_loops(
    bedpe_file = tmp_bedpe,
    txdb = txdb_obj, org_db = "org.Hs.eg.db", species = "hg19",
    min_expr = 0.1, conflict_strategy = "biotype_first",
    out_dir = tempdir(), write_output = FALSE, quiet = TRUE
  )
  refined <- refine_loop_anchors_by_expression(base_res,
    expr_matrix_file = expr_path, sample_columns = c("con1", "con2"),
    threshold = 1, reclassify_by_expression = FALSE,
    out_dir = tempdir(), write_output = FALSE, quiet = TRUE
  )
  state <- attr(refined, "looplook_anchor_state")
  # conflict_min_expr must survive unchanged
  expect_equal(state$gene_assignment_policy$conflict_min_expr, 0.1,
    info = "conflict_min_expr must not be overwritten by expression threshold"
  )
  # expression_policy must record the refinement threshold
  expect_equal(state$expression_policy$effective_threshold, 1,
    info = "expression_policy$effective_threshold must match the refinement threshold"
  )
  expect_true(
    state$gene_assignment_policy$conflict_min_expr !=
      state$expression_policy$effective_threshold,
    info = "conflict_min_expr and effective_threshold must be independent"
  )
  unlink(tmp_bedpe)
})


# ══════════════════════════════════════════════════════════════════════════
# 13. .reannotate_tss_genes: anchor ID alignment via metadata column
# ══════════════════════════════════════════════════════════════════════════

test_that(".reannotate_tss_genes returns valid structure for single anchor", {
  skip_if_not_installed("org.Hs.eg.db")
  skip_if(is.null(get_test_txdb()), "Sample TxDb not available")
  txdb_obj <- get_test_txdb()
  validation <- data.frame(
    anchor_id = "A1", chr = "chr6", start = 10412000L, end = 10413000L,
    anchor_type = "E", stringsAsFactors = FALSE
  )
  result <- looplook:::.reannotate_tss_genes(
    "A1", validation, txdb_obj,
    org_db_pkg = "org.Hs.eg.db",
    tss_region = c(-2000, 2000), quiet = TRUE
  )
  # Core contract: output must have these columns with correct types
  expect_equal(result$anchor_id, "A1")
  expect_true("alignment_method" %in% colnames(result))
  expect_true("gene_after" %in% colnames(result))
  expect_true("TSS_supported" %in% colnames(result))
  expect_type(result$TSS_supported, "logical")
})

test_that(".reannotate_tss_genes handles duplicate coordinates without crash", {
  skip_if_not_installed("org.Hs.eg.db")
  skip_if(is.null(get_test_txdb()), "Sample TxDb not available")
  txdb_obj <- get_test_txdb()
  validation <- data.frame(
    anchor_id = c("A1", "A2"),
    chr = c("chr6", "chr6"),
    start = c(10412000L, 10412000L),
    end = c(10413000L, 10413000L),
    anchor_type = c("E", "E"),
    stringsAsFactors = FALSE
  )
  result <- looplook:::.reannotate_tss_genes(
    c("A1", "A2"), validation, txdb_obj,
    org_db_pkg = "org.Hs.eg.db",
    tss_region = c(-2000, 2000), quiet = TRUE
  )
  # Both anchors must return rows in order
  expect_setequal(result$anchor_id, c("A1", "A2"))
  expect_equal(nrow(result), 2L)
  # alignment_method must NOT be NA (some matching strategy was attempted)
  expect_false(any(is.na(result$alignment_method)),
    info = "alignment_method must be populated for every anchor"
  )
})

test_that(".reannotate_tss_genes returns gene_after=NA for unannotated region", {
  skip_if_not_installed("org.Hs.eg.db")
  skip_if(is.null(get_test_txdb()), "Sample TxDb not available")
  txdb_obj <- get_test_txdb()
  # Use chr6 at a coordinate unlikely to overlap any gene in the sample TxDb
  validation <- data.frame(
    anchor_id = "A_Nogene",
    chr = "chr6", start = 1L, end = 100L,
    anchor_type = "E",
    stringsAsFactors = FALSE
  )
  result <- looplook:::.reannotate_tss_genes(
    "A_Nogene", validation, txdb_obj,
    org_db_pkg = "org.Hs.eg.db",
    tss_region = c(-2000, 2000), quiet = TRUE
  )
  expect_equal(result$anchor_id, "A_Nogene")
  expect_true(is.na(result$gene_after) || result$gene_after == "")
  expect_false(result$TSS_supported)
})


# ══════════════════════════════════════════════════════════════════════════
# 14. BEDPE/BED coordinate validation (P2-2)
# ══════════════════════════════════════════════════════════════════════════

test_that("BEDPE validator rejects negative coordinates", {
  tmp <- tempfile(fileext = ".bedpe")
  writeLines("chr1\t-100\t200\tchr1\t300\t400", tmp)
  expect_error(
    looplook:::.validate_bedpe_df(data.table::fread(tmp), quiet = TRUE),
    "negative"
  )
  unlink(tmp)
})

test_that("BEDPE validator rejects non-finite coordinates", {
  tmp <- tempfile(fileext = ".bedpe")
  writeLines("chr1\tInf\t200\tchr1\t300\t400", tmp)
  expect_error(
    looplook:::.validate_bedpe_df(data.table::fread(tmp), quiet = TRUE),
    "non-finite"
  )
  unlink(tmp)
})

test_that("BEDPE validator rejects empty chromosome names", {
  tmp <- tempfile(fileext = ".bedpe")
  writeLines("\t100\t200\tchr1\t300\t400", tmp)
  expect_error(
    looplook:::.validate_bedpe_df(data.table::fread(tmp, header = FALSE), quiet = TRUE),
    "chromosome"
  )
  unlink(tmp)
})


# ══════════════════════════════════════════════════════════════════════════
# 15. empty_annotation_result schema (P2-10)
# ══════════════════════════════════════════════════════════════════════════

test_that("empty_annotation_result returns complete schema with qc_status", {
  res <- looplook:::.empty_annotation_result(
    "No valid anchors", "hg38", NULL, NULL
  )
  expect_true(is.list(res))
  expect_true("metadata" %in% names(res))
  expect_equal(res$metadata$parameters$qc_status, "empty_input")
  expect_true(nzchar(res$metadata$parameters$qc_reason))
  expect_equal(nrow(res$loop_annotation), 0L)
})


# ══════════════════════════════════════════════════════════════════════════
# 16. distal-element stats: role-based promoter counting (P1)
# ══════════════════════════════════════════════════════════════════════════

test_that("distal stats use Neighbor_Role when role columns are present", {
  loop_df <- data.frame(
    a1_id = c("E1", "E1"), a2_id = c("P1", "G1"),
    chr1 = "chr1", start1 = 100L, end1 = 200L,
    chr2 = "chr1", start2 = 300L, end2 = 400L,
    anchor1_type = c("E", "E"), anchor2_type = c("P", "P"),
    anchor1_gene = c("", ""), anchor2_gene = c("GENE_P", "GENE_G"),
    anchor1_gene_role = c(NA_character_, NA_character_),
    anchor2_gene_role = c("promoter", "gene_body"),
    loop_type = c("E-P", "E-P"),
    cluster_id = c("c1", "c1"),
    stringsAsFactors = FALSE
  )
  res <- looplook:::.build_distal_element_df(loop_df, 0.95)
  expect_equal(res$n_Linked_Promoters, 1L,
    info = "Only the promoter-role neighbor should count as promoter-linked"
  )
  expect_equal(res$Target_Genes, "GENE_P",
    info = "gene_body neighbor should not appear in Target_Genes"
  )
})

test_that("distal stats fall back to type-based when role columns are absent", {
  loop_df <- data.frame(
    a1_id = "E1", a2_id = "P1",
    chr1 = "chr1", start1 = 100L, end1 = 200L,
    chr2 = "chr1", start2 = 300L, end2 = 400L,
    anchor1_type = "E", anchor2_type = "P",
    anchor1_gene = "", anchor2_gene = "GENE1",
    loop_type = "E-P", cluster_id = "c1",
    stringsAsFactors = FALSE
  )
  res <- looplook:::.build_distal_element_df(loop_df, 0.95)
  expect_equal(res$n_Linked_Promoters, 1L,
    info = "Without role columns, type-based promoter detection should work"
  )
})


# ══════════════════════════════════════════════════════════════════════════
# 17. promoter-centric stats: enhancer-like count (P1)
# ══════════════════════════════════════════════════════════════════════════

test_that("promoter stats count pure-E neighbor as enhancer-like", {
  loop_df <- data.frame(
    a1_id = "P1", a2_id = "E1",
    anchor1_type = "P", anchor2_type = "E",
    anchor1_gene = "GENE1", anchor2_gene = "",
    anchor1_gene_role = "promoter", anchor2_gene_role = NA_character_,
    loop_type = "E-P",
    stringsAsFactors = FALSE
  )
  raw <- looplook:::.compute_raw_promoter_stats(loop_df)
  expect_equal(raw$n_Linked_EnhancerLike_Filtered, 1L,
    info = "Pure E neighbor must appear in enhancer-like count"
  )
  expect_equal(raw$n_Linked_Promoters_Filtered, 0L,
    info = "E neighbor is not a promoter connection"
  )
})


# ══════════════════════════════════════════════════════════════════════════
# 18. hub_percentile propagation (P1-3)
# ══════════════════════════════════════════════════════════════════════════

test_that("promoter stats respect non-default hub_percentile", {
  # 5 hub genes (4 loops each), 5 non-hub (1 loop each).
  # 50th percentile = 1, so hub genes should be > floor(3) and > p50.
  genes <- c(
    rep(c("HUB1", "HUB2", "HUB3", "HUB4", "HUB5"), each = 4),
    "LOW1", "LOW2", "LOW3", "LOW4", "LOW5"
  )
  loop_df <- data.frame(
    a1_id = paste0("A", seq_along(genes)),
    a2_id = paste0("B", seq_along(genes)),
    anchor1_type = rep("P", length(genes)),
    anchor2_type = rep("E", length(genes)),
    anchor1_gene = genes, anchor2_gene = "",
    anchor1_gene_role = rep("promoter", length(genes)),
    anchor2_gene_role = rep(NA_character_, length(genes)),
    loop_type = rep("E-P", length(genes)),
    stringsAsFactors = FALSE
  )
  raw <- looplook:::.compute_raw_promoter_stats(loop_df)
  vals <- setNames(rep(10, length(unique(genes))), unique(genes))
  res_50 <- looplook:::.build_promoter_centric_df(raw, NULL,
    vals = vals, threshold = 0, hub_percentile = 0.50
  )
  res_95 <- looplook:::.build_promoter_centric_df(raw, NULL,
    vals = vals, threshold = 0, hub_percentile = 0.95
  )
  n_hub_50 <- sum(res_50$Is_Regulatory_Hub == "Yes")
  n_hub_95 <- sum(res_95$Is_Regulatory_Hub == "Yes")
  expect_true(n_hub_50 >= n_hub_95,
    info = "Lower percentile should produce at least as many hubs"
  )
})


# ══════════════════════════════════════════════════════════════════════════
# 19. DAG max-support shortest path (P2-1 fix)
# ══════════════════════════════════════════════════════════════════════════

test_that("DAG path: vertex names are matched correctly (not as.integer)", {
  skip_if_not_installed("igraph")
  # 6-node ring: A1→A5 has shortest path A1-A6-A5 (2 hops) via the back
  g <- igraph::make_ring(6, directed = FALSE)
  igraph::V(g)$name <- c("A1", "A2", "A3", "A4", "A5", "A6")
  igraph::E(g)$n_support <- rep(1, 6)
  p <- looplook:::.dag_max_support_path(g, "A1", "A4", 3)
  expect_equal(igraph::as_ids(p), c("A1", "A2", "A3", "A4"))
})

test_that("DAG path: edge input order does not affect result", {
  skip_if_not_installed("igraph")
  # Two identical graphs with reversed edge-input order.
  edges_fwd <- matrix(c("A1", "A2", "A2", "A3", "A3", "A4"), ncol = 2, byrow = TRUE)
  edges_rev <- matrix(c("A4", "A3", "A3", "A2", "A2", "A1"), ncol = 2, byrow = TRUE)
  g_fwd <- igraph::graph_from_edgelist(edges_fwd, directed = FALSE)
  igraph::E(g_fwd)$n_support <- c(2, 1, 3)
  g_rev <- igraph::graph_from_edgelist(edges_rev, directed = FALSE)
  igraph::E(g_rev)$n_support <- c(3, 1, 2)
  p_fwd <- igraph::as_ids(looplook:::.dag_max_support_path(g_fwd, "A1", "A4", 3))
  p_rev <- igraph::as_ids(looplook:::.dag_max_support_path(g_rev, "A1", "A4", 3))
  expect_equal(p_fwd, c("A1", "A2", "A3", "A4"))
  expect_equal(p_rev, c("A1", "A2", "A3", "A4"))
})

test_that("DAG path: prefers higher support among equal-hop paths", {
  skip_if_not_installed("igraph")
  # 6-node ring, 3-hop: A1-A2-A3-A4 (support sum 1+5+3=9) vs A1-A6-A5-A4 (1+1+1=3)
  g <- igraph::make_ring(6, directed = FALSE)
  igraph::V(g)$name <- c("A1", "A2", "A3", "A4", "A5", "A6")
  igraph::E(g)$n_support <- c(1, 5, 3, 1, 1, 1)
  p <- looplook:::.dag_max_support_path(g, "A1", "A4", 3)
  # Should pick A1-A2-A3-A4 (support 9) over A1-A6-A5-A4 (support 3)
  expect_equal(igraph::as_ids(p), c("A1", "A2", "A3", "A4"))
})

test_that("DAG path: minimum hop is first objective (3-hop low > 4-hop high)", {
  skip_if_not_installed("igraph")
  # 3-hop path: A1-B1-C1-A4, support = 1,1,1
  # 4-hop path: A1-D1-E1-F1-A4, support = 100,100,100,100
  edges <- matrix(c(
    "A1", "B1", "B1", "C1", "C1", "A4",
    "A1", "D1", "D1", "E1", "E1", "F1", "F1", "A4"
  ), ncol = 2, byrow = TRUE)
  g <- igraph::graph_from_edgelist(edges, directed = FALSE)
  igraph::E(g)$n_support <- c(1, 1, 1, 100, 100, 100, 100)
  p <- looplook:::.dag_max_support_path(g, "A1", "A4", 3)
  expect_equal(igraph::as_ids(p), c("A1", "B1", "C1", "A4"),
    info = "Must choose 3-hop path even though 4-hop has higher support"
  )
})

test_that("DAG path: returns NULL when no path exists (disconnected)", {
  skip_if_not_installed("igraph")
  g <- igraph::make_ring(3, directed = FALSE)
  igraph::V(g)$name <- c("A1", "A2", "A3")
  igraph::E(g)$n_support <- c(1, 1, 1)
  g <- igraph::add_vertices(g, 2, name = c("B1", "B2"))
  g <- igraph::add_edges(g, c("B1", "B2"))
  igraph::E(g)$n_support[is.na(igraph::E(g)$n_support)] <- 1
  expect_null(looplook:::.dag_max_support_path(g, "A1", "B2", 1))
})


# ══════════════════════════════════════════════════════════════════════════
# 21. d≥3 expanded path → strict target propagation invariant
# ══════════════════════════════════════════════════════════════════════════

test_that("3-hop promoter gene enters Expanded but NOT Assigned (primary guard)", {
  links <- data.frame(
    input_id = "peak1",
    gene = "GENE3",
    gene_role = "promoter",
    source = "loop_anchor",
    evidence = "expanded_promoter_loop",
    anchor_role = "expanded_anchor",
    path_length = 3L,
    strict_assignment_eligible = TRUE,
    stringsAsFactors = FALSE
  )
  bed <- data.frame(input_id = "peak1", stringsAsFactors = FALSE)
  res <- looplook:::.aggregate_strict_targets(links, bed, max_primary_hop = 1L)
  expect_equal(res$Expanded_Target_Genes, "GENE3",
    info = "3-hop gene must enter Expanded"
  )
  expect_true(
    is.na(res$Assigned_Target_Genes) ||
      res$Assigned_Target_Genes == "",
    info = "3-hop gene must NOT enter Assigned (primary)"
  )
})

test_that("3-hop gene-body enters Expanded not Assigned", {
  links <- data.frame(
    input_id = "peak1",
    gene = "GENE3",
    gene_role = "gene_body",
    source = "loop_anchor",
    evidence = "expanded_gene_body_context",
    anchor_role = "expanded_anchor",
    path_length = 3L,
    strict_assignment_eligible = TRUE,
    stringsAsFactors = FALSE
  )
  bed <- data.frame(input_id = "peak1", stringsAsFactors = FALSE)
  res <- looplook:::.aggregate_strict_targets(links, bed, max_primary_hop = 1L)
  expect_equal(res$Expanded_Target_Genes, "GENE3")
  expect_true(is.na(res$Assigned_Target_Genes) || res$Assigned_Target_Genes == "")
})

test_that("3-hop positional_candidate enters neither Expanded nor Assigned", {
  links <- data.frame(
    input_id = "peak1",
    gene = "GENE3",
    gene_role = "positional_candidate",
    source = "loop_anchor",
    evidence = "positional_candidate_after_chromatin",
    anchor_role = "expanded_anchor",
    path_length = 3L,
    strict_assignment_eligible = FALSE,
    stringsAsFactors = FALSE
  )
  bed <- data.frame(input_id = "peak1", stringsAsFactors = FALSE)
  res <- looplook:::.aggregate_strict_targets(links, bed, max_primary_hop = 1L)
  expect_equal(res$All_Loop_Connected_Genes, "GENE3",
    info = "3-hop candidate still in all-connected"
  )
  expect_true(is.na(res$Expanded_Target_Genes) || res$Expanded_Target_Genes == "",
    info = "positional_candidate must not enter Expanded"
  )
  expect_true(is.na(res$Assigned_Target_Genes) || res$Assigned_Target_Genes == "",
    info = "positional_candidate must not enter Assigned"
  )
})

test_that("1-hop promoter with path_length=0 correctly enters Assigned", {
  links <- data.frame(
    input_id = "peak1",
    gene = "GENE1",
    gene_role = "promoter",
    source = "loop_anchor",
    evidence = "local_promoter_overlap",
    anchor_role = "local_anchor",
    path_length = 0L,
    strict_assignment_eligible = TRUE,
    stringsAsFactors = FALSE
  )
  bed <- data.frame(input_id = "peak1", stringsAsFactors = FALSE)
  res <- looplook:::.aggregate_strict_targets(links, bed, max_primary_hop = 1L)
  expect_equal(res$Assigned_Target_Genes, "GENE1",
    info = "Direct promoter gene must enter Assigned"
  )
  expect_true(is.na(res$Expanded_Target_Genes) || res$Expanded_Target_Genes == "",
    info = "Direct promoter must not appear in Expanded"
  )
})

test_that("DAG path: returns NULL for non-existent vertex name", {
  skip_if_not_installed("igraph")
  g <- igraph::make_ring(3, directed = FALSE)
  igraph::V(g)$name <- c("A1", "A2", "A3")
  igraph::E(g)$n_support <- c(1, 1, 1)
  expect_null(looplook:::.dag_max_support_path(g, "A1", "NOPE", 1))
})


# ══════════════════════════════════════════════════════════════════════════
# 20. canonical anchor registry + distal stats (P1-3 fix)
# ══════════════════════════════════════════════════════════════════════════

test_that("distal stats with anchor_registry uses canonical coordinates (not raw)", {
  # Deliberately make raw loop coords different from canonical merged coords.
  gr <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1"),
    ranges = IRanges::IRanges(start = c(100, 280), end = c(200, 420)),
    anchor_id = c("A1", "A2"),
    cluster_id = c("canonical_cluster", "canonical_cluster")
  )
  names(gr) <- c("A1", "A2")
  loop_df <- data.frame(
    a1_id = "A1", a2_id = "A2",
    chr1 = "chr1", start1 = 100L, end1 = 200L,
    chr2 = "chr1", start2 = 300L, end2 = 350L,
    anchor1_type = "P", anchor2_type = "E",
    anchor1_gene = "GENE1", anchor2_gene = "",
    anchor1_gene_role = "promoter",
    anchor2_gene_role = NA_character_,
    loop_type = "E-P", cluster_id = "raw_cluster",
    stringsAsFactors = FALSE
  )
  res <- looplook:::.build_distal_element_df(loop_df, 0.95,
    anchor_registry = gr
  )
  # Must use canonical registry, not loop_df raw coords
  expect_equal(res$start, 280L) # canonical: 280, raw: 300
  expect_equal(res$end, 420L) # canonical: 420, raw: 350
  expect_equal(res$cluster_id, "canonical_cluster")
  expect_true("cluster_id" %in% colnames(res))
})

test_that("distal stats with registry: row-order invariance (2 rows)", {
  gr <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1", "chr1"),
    ranges = IRanges::IRanges(start = c(100, 200, 280), end = c(150, 250, 420)),
    anchor_id = c("P1", "P2", "E1"),
    cluster_id = c("c1", "c1", "c1")
  )
  names(gr) <- c("P1", "P2", "E1")
  # Two loop rows both reference distal anchor E1, different raw coords
  base <- data.frame(
    a1_id = c("P1", "P2"),
    a2_id = c("E1", "E1"),
    chr1 = c("chr1", "chr1"),
    start1 = c(100L, 200L),
    end1 = c(150L, 250L),
    chr2 = c("chr1", "chr1"),
    start2 = c(300L, 340L),
    end2 = c(350L, 390L),
    anchor1_type = c("P", "P"),
    anchor2_type = c("E", "E"),
    anchor1_gene = c("G1", "G2"),
    anchor2_gene = c("", ""),
    anchor1_gene_role = c("promoter", "promoter"),
    anchor2_gene_role = c(NA_character_, NA_character_),
    loop_type = c("E-P", "E-P"),
    cluster_id = c("raw", "raw"),
    stringsAsFactors = FALSE
  )
  res1 <- looplook:::.build_distal_element_df(base, 0.95,
    anchor_registry = gr
  )
  res2 <- looplook:::.build_distal_element_df(base[2:1, ], 0.95,
    anchor_registry = gr
  )
  expect_equal(res1$start, 280L)
  expect_equal(res1$end, 420L)
  expect_equal(res1$cluster_id, "c1")
  expect_equal(res1$start, res2$start)
  expect_equal(res1$end, res2$end)
  expect_equal(res1$cluster_id, res2$cluster_id)
  expect_equal(res1$Total_Loops, res2$Total_Loops)
})

test_that("distal stats without anchor_registry still works (backward compat)", {
  loop_df <- data.frame(
    a1_id = "A1", a2_id = "A2",
    chr1 = "chr1", start1 = 100L, end1 = 200L,
    chr2 = "chr1", start2 = 300L, end2 = 400L,
    anchor1_type = "P", anchor2_type = "E",
    anchor1_gene = "GENE1", anchor2_gene = "",
    loop_type = "E-P", cluster_id = "c1",
    stringsAsFactors = FALSE
  )
  res <- looplook:::.build_distal_element_df(loop_df, 0.95)
  expect_s3_class(res, "data.frame")
  expect_true(nrow(res) > 0)
})


# ══════════════════════════════════════════════════════════════════════════
# 22. DAG multi-component: no coercion warning (P1)
# ══════════════════════════════════════════════════════════════════════════

test_that("DAG path: ignores unrelated disconnected component without warning", {
  skip_if_not_installed("igraph")
  # A1-A2-A3-A4 is the target path; B1-B2 is an unrelated component.
  edges <- matrix(c(
    "A1", "A2",
    "A2", "A3",
    "A3", "A4",
    "B1", "B2"
  ), ncol = 2, byrow = TRUE)
  g <- igraph::graph_from_edgelist(edges, directed = FALSE)
  igraph::E(g)$n_support <- c(1, 1, 1, 1)
  expect_no_warning({
    p <- looplook:::.dag_max_support_path(g, "A1", "A4", 3)
  })
  expect_equal(igraph::as_ids(p), c("A1", "A2", "A3", "A4"))
})


# ══════════════════════════════════════════════════════════════════════════
# 23. End-to-end: graph path → target gene links propagation (P1)
# ══════════════════════════════════════════════════════════════════════════

test_that("end-to-end: 3-hop graph path sets path_length=3 and anchor_role=expanded_anchor", {
  skip_if_not_installed("igraph")
  # 3-hop chain A1-A2-A3-A4; A4 is a promoter with gene GENE4
  edges <- matrix(c(
    "A1", "A2",
    "A2", "A3",
    "A3", "A4"
  ), ncol = 2, byrow = TRUE)
  g <- igraph::graph_from_edgelist(edges, directed = FALSE)
  igraph::E(g)$n_support <- c(1, 1, 1)

  hit_df <- data.frame(
    qid = 1L, sid = 1L, anchor_id = "A1",
    stringsAsFactors = FALSE
  )
  loop_annotation_final <- data.frame(
    loop_ID = c("L1", "L2", "L3"),
    a1_id = c("A1", "A2", "A3"),
    a2_id = c("A2", "A3", "A4"),
    stringsAsFactors = FALSE
  )
  map_info <- data.frame(
    anchor_id = c("A1", "A2", "A3", "A4"),
    type_code = c("E", "E", "E", "P"),
    SYMBOL = c("", "", "", "GENE4"),
    stringsAsFactors = FALSE
  )
  ego_list_target <- list(A1 = igraph::V(g)[c("A1", "A2", "A3", "A4")])
  bed_info <- data.frame(input_id = "Peak_1", stringsAsFactors = FALSE)

  result <- looplook:::.build_target_gene_links(
    hit_df = hit_df,
    bed_info = bed_info,
    loop_annotation_final = loop_annotation_final,
    map_info = map_info,
    ego_list_target = ego_list_target,
    g = g
  )

  gene4_rows <- result[result$gene == "GENE4" &
    result$anchor_role == "expanded_anchor", ]
  expect_equal(nrow(gene4_rows), 1L,
    info = "GENE4 should appear once as expanded_anchor"
  )
  expect_equal(gene4_rows$path_length, 3L,
    info = "3-hop graph path must produce path_length = 3, not 0 or NA"
  )
  expect_equal(gene4_rows$evidence, "expanded_promoter_loop",
    info = "promoter at expanded distance must get expanded_promoter_loop evidence"
  )
})


# ══════════════════════════════════════════════════════════════════════════
# 24. Chromatin shared input validator (P2)
# ══════════════════════════════════════════════════════════════════════════

test_that("chromatin validator rejects unnamed chromatin_beds list", {
  expect_error(
    looplook:::.validate_chromatin_overlap_inputs(
      chromatin_beds = list("a.bed", "b.bed"),
      candidate_types = NULL,
      anchor_gap = 200,
      anchor_min_overlap = 100
    ),
    "named list"
  )
})

# ══════════════════════════════════════════════════════════════════════════
# Metadata inheritance: neighbor_hop across refinement stages
# ══════════════════════════════════════════════════════════════════════════

test_that("basic annotation records neighbor_hop in metadata and state", {
  skip_if_not_installed("org.Hs.eg.db")
  skip_if(is.null(get_test_txdb()), "Sample TxDb not available")

  res <- .get_shared_hop1_base()

  expect_equal(res$metadata$parameters$neighbor_hop, 1L)
  expect_equal(attr(res, "looplook_anchor_state")$neighbor_hop, 1L)
})

test_that("metadata: expression refinement inherits neighbor_hop (state >= metadata)", {
  skip_if_not_installed("org.Hs.eg.db")
  skip_if(is.null(get_test_txdb()), "Sample TxDb not available")

  res <- .get_shared_hop1_base()

  expr_path <- system.file("extdata", "example_tpm.txt", package = "looplook")
  skip_if(expr_path == "", "Example TPM file not available")
  refined <- refine_loop_anchors_by_expression(res,
    expr_matrix_file = expr_path, sample_columns = c("con1", "con2"),
    threshold = 1, out_dir = tempdir(), write_output = FALSE, quiet = TRUE)

  expect_equal(refined$metadata$parameters$neighbor_hop, 1L)
  expect_equal(attr(refined, "looplook_anchor_state")$neighbor_hop, 1L)
})

# ══════════════════════════════════════════════════════════════════════════
# Regulated promoter evidence strict filtering tests
# ══════════════════════════════════════════════════════════════════════════

test_that("regulated evidence: strict promoter retained", {
  links <- data.frame(
    input_id = "P1", source = "loop_anchor",
    gene = "G1", gene_role = "promoter",
    strict_assignment_eligible = TRUE,
    path_length = 0L, evidence = "local_promoter_overlap",
    stringsAsFactors = FALSE
  )
  res <- looplook:::.summarise_regulated_promoter_evidence(links)
  expect_false(res$Regulated_promoter_Evidence[1] == "none")
})

test_that("regulated evidence: non-strict promoter excluded", {
  links <- data.frame(
    input_id = "P1", source = "loop_anchor",
    gene = "G1", gene_role = "promoter",
    strict_assignment_eligible = FALSE,
    path_length = 0L, evidence = "local_promoter_overlap",
    stringsAsFactors = FALSE
  )
  res <- looplook:::.summarise_regulated_promoter_evidence(links)
  expect_equal(nrow(res), 0L)
})

test_that("regulated evidence: old schema promoter (no strict col) retained", {
  links <- data.frame(
    input_id = "P1", source = "loop_anchor",
    gene = "G1", gene_role = "promoter",
    path_length = 0L, evidence = "local_promoter_overlap",
    stringsAsFactors = FALSE
  )
  res <- looplook:::.summarise_regulated_promoter_evidence(links)
  expect_false(res$Regulated_promoter_Evidence[1] == "none")
})

test_that("chromatin validator rejects NA or empty mark names", {
  expect_error(
    looplook:::.validate_chromatin_overlap_inputs(
      chromatin_beds = setNames(list("a.bed"), NA_character_),
      candidate_types = NULL,
      anchor_gap = 200,
      anchor_min_overlap = 100
    ),
    "names must not"
  )
  expect_error(
    looplook:::.validate_chromatin_overlap_inputs(
      chromatin_beds = setNames(list("a.bed"), ""),
      candidate_types = NULL,
      anchor_gap = 200,
      anchor_min_overlap = 100
    ),
    "names must not"
  )
})

test_that("chromatin validator rejects unknown candidate_type", {
  expect_error(
    looplook:::.validate_chromatin_overlap_inputs(
      chromatin_beds = list(H3K4me1 = "x.bed"),
      candidate_types = "EP",
      anchor_gap = 200,
      anchor_min_overlap = 100
    ),
    "Unknown candidate_type"
  )
})

test_that("chromatin validator rejects duplicate candidate_types", {
  expect_error(
    looplook:::.validate_chromatin_overlap_inputs(
      chromatin_beds = list(H3K4me1 = "x.bed"),
      candidate_types = c("E", "E"),
      anchor_gap = 200,
      anchor_min_overlap = 100
    ),
    "unique"
  )
})

test_that("chromatin validator rejects non-integer or non-finite anchor_gap", {
  expect_error(
    looplook:::.validate_chromatin_overlap_inputs(
      chromatin_beds = list(H3K4me1 = "x.bed"),
      candidate_types = NULL,
      anchor_gap = 1.5,
      anchor_min_overlap = 100
    ),
    "finite integer"
  )
  expect_error(
    looplook:::.validate_chromatin_overlap_inputs(
      chromatin_beds = list(H3K4me1 = "x.bed"),
      candidate_types = NULL,
      anchor_gap = Inf,
      anchor_min_overlap = 100
    ),
    "finite integer"
  )
})

test_that("chromatin validator rejects non-integer anchor_min_overlap", {
  expect_error(
    looplook:::.validate_chromatin_overlap_inputs(
      chromatin_beds = list(H3K4me1 = "x.bed"),
      candidate_types = NULL,
      anchor_gap = 200,
      anchor_min_overlap = 1.5
    ),
    "finite positive integer"
  )
  expect_error(
    looplook:::.validate_chromatin_overlap_inputs(
      chromatin_beds = list(H3K4me1 = "x.bed"),
      candidate_types = NULL,
      anchor_gap = 200,
      anchor_min_overlap = 0
    ),
    "finite positive integer"
  )
})


# ══════════════════════════════════════════════════════════════════════════
# 25. Consolidation: reject duplicate replicate file paths (P1)
# ══════════════════════════════════════════════════════════════════════════

test_that("consolidation rejects duplicate replicate file paths", {
  f1 <- tempfile(fileext = ".bedpe")
  writeLines("chr1\t100\t200\tchr1\t500\t600", f1)
  expect_error(
    consolidate_chromatin_loops(
      files = c(f1, f1),
      write_output = FALSE
    ),
    "Duplicate BEDPE replicate"
  )
  unlink(f1)
})

test_that("consolidation rejects replicate paths that resolve to same file", {
  f1 <- tempfile(fileext = ".bedpe")
  writeLines("chr1\t100\t200\tchr1\t500\t600", f1)
  f1_rel <- basename(f1)
  old <- setwd(dirname(f1))
  on.exit(setwd(old))
  expect_error(
    consolidate_chromatin_loops(
      files = c(f1_rel, f1),
      write_output = FALSE
    ),
    "Duplicate BEDPE replicate"
  )
  unlink(f1)
})


# ══════════════════════════════════════════════════════════════════════════
# 26. Consolidation QC status and metadata (P2)
# ══════════════════════════════════════════════════════════════════════════

test_that("consolidation empty_after_raw_score_filter has correct metadata", {
  f1 <- tempfile(fileext = ".bedpe")
  f2 <- tempfile(fileext = ".bedpe")
  writeLines("chr1\t100\t200\tchr1\t300\t400\t1", f1)
  writeLines("chr1\t150\t250\tchr1\t350\t450\t2", f2)
  res <- consolidate_chromatin_loops(
    files = c(f1, f2), mode = "consensus",
    min_raw_score = 10, quiet = TRUE
  )
  meta <- attr(res, "looplook_metadata")
  expect_equal(meta$parameters$qc_status, "empty_after_raw_score_filter")
  expect_match(meta$parameters$qc_reason, "min_raw_score")
  diag <- meta$diagnostics
  expect_true("input_counts" %in% names(diag))
  expect_equal(diag$input_counts$raw_count, c(1L, 1L))
  expect_equal(diag$input_counts$retained_count, c(0L, 0L))
  unlink(c(f1, f2))
})

test_that("consolidation records effective min_consensus for N=4", {
  fs <- replicate(4, {
    f <- tempfile(fileext = ".bedpe")
    writeLines("chr1\t100\t200\tchr1\t300\t400", f)
    f
  })
  on.exit(unlink(fs))
  # N=4, min_consensus=NULL → ceiling(0.75*4) = 3
  res <- consolidate_chromatin_loops(
    files = fs, mode = "consensus",
    gap = 0, quiet = TRUE
  )
  meta <- attr(res, "looplook_metadata")
  expect_equal(meta$parameters$min_consensus_effective, 3L)
  expect_null(meta$parameters$min_consensus_requested)
})

test_that("consolidation input_counts in metadata for non-empty result", {
  f1 <- tempfile(fileext = ".bedpe")
  f2 <- tempfile(fileext = ".bedpe")
  writeLines("chr1\t100\t200\tchr1\t300\t400", f1)
  writeLines("chr1\t100\t200\tchr1\t300\t400", f2)
  res <- consolidate_chromatin_loops(
    files = c(f1, f2), mode = "consensus",
    gap = 0, quiet = TRUE
  )
  meta <- attr(res, "looplook_metadata")
  expect_equal(meta$parameters$qc_status, "ok")
  diag <- meta$diagnostics
  expect_true("input_counts" %in% names(diag))
  # Both files have 1 loop each, both survive (no min_raw_score filtering)
  expect_equal(diag$input_counts$raw_count, c(1L, 1L))
  expect_equal(diag$input_counts$retained_count, c(1L, 1L))
  unlink(c(f1, f2))
})

test_that("consolidation empty_after_merge when no cluster meets consensus threshold", {
  f1 <- tempfile(fileext = ".bedpe")
  f2 <- tempfile(fileext = ".bedpe")
  # Two loops on different chromosomes – cannot merge, each cluster n_reps=1
  writeLines("chr1\t100\t200\tchr1\t300\t400", f1)
  writeLines("chr2\t100\t200\tchr2\t300\t400", f2)
  res <- consolidate_chromatin_loops(
    files = c(f1, f2), mode = "consensus",
    gap = 0, quiet = TRUE
  )
  meta <- attr(res, "looplook_metadata")
  expect_equal(meta$parameters$qc_status, "empty_after_merge")
  expect_match(meta$parameters$qc_reason, "consensus")
  unlink(c(f1, f2))
})

test_that("consolidation empty_after_postfilter when min_score removes everything", {
  f1 <- tempfile(fileext = ".bedpe")
  f2 <- tempfile(fileext = ".bedpe")
  # Identical loops on both files → they merge with sufficient consensus.
  # post-filter min_score=99 removes them all.
  writeLines("chr1\t100\t200\tchr1\t300\t400", f1)
  writeLines("chr1\t100\t200\tchr1\t300\t400", f2)
  res <- consolidate_chromatin_loops(
    files = c(f1, f2), mode = "consensus",
    gap = 0, min_score = 99, quiet = TRUE
  )
  meta <- attr(res, "looplook_metadata")
  expect_equal(meta$parameters$qc_status, "empty_after_postfilter")
  expect_match(meta$parameters$qc_reason, "min_score")
  unlink(c(f1, f2))
})


# ══════════════════════════════════════════════════════════════════════════
# 27. Empty replicate file policy (P1)
# ══════════════════════════════════════════════════════════════════════════

test_that("consolidation: all files empty → empty_input with warning", {
  f1 <- tempfile(fileext = ".bedpe")
  f2 <- tempfile(fileext = ".bedpe")
  file.create(f1)
  file.create(f2)
  expect_warning(
    res <- consolidate_chromatin_loops(
      files = c(f1, f2), mode = "consensus", quiet = TRUE
    ),
    "empty.*0 bytes"
  )
  meta <- attr(res, "looplook_metadata")
  expect_equal(meta$parameters$qc_status, "empty_input")
  expect_equal(meta$diagnostics$input_counts$raw_count, c(0L, 0L))
  expect_equal(meta$diagnostics$input_counts$retained_count, c(0L, 0L))
  unlink(c(f1, f2))
})

test_that("consolidation: one empty + one normal → empty_after_merge with warning", {
  f1 <- tempfile(fileext = ".bedpe")
  f2 <- tempfile(fileext = ".bedpe")
  file.create(f1)
  writeLines("chr1\t100\t200\tchr1\t300\t400", f2)
  expect_warning(
    res <- consolidate_chromatin_loops(
      files = c(f1, f2), mode = "consensus",
      gap = 0, quiet = TRUE
    ),
    "empty.*0 bytes"
  )
  meta <- attr(res, "looplook_metadata")
  expect_equal(meta$parameters$qc_status, "empty_after_merge")
  expect_equal(meta$parameters$min_consensus_effective, 2L)
  unlink(c(f1, f2))
})

test_that("consolidation: N=4 with 1 empty + 3 overlapping files retains the loop", {
  f_empty <- tempfile(fileext = ".bedpe")
  f1 <- tempfile(fileext = ".bedpe")
  f2 <- tempfile(fileext = ".bedpe")
  f3 <- tempfile(fileext = ".bedpe")
  file.create(f_empty)
  writeLines("chr1\t100\t200\tchr1\t300\t400", f1)
  writeLines("chr1\t100\t200\tchr1\t300\t400", f2)
  writeLines("chr1\t100\t200\tchr1\t300\t400", f3)
  # N=4, min_consensus=NULL → ceiling(0.75*4) = 3
  # 3 normal files support the loop → n_reps=3 >= 3 → retained
  expect_warning(
    res <- consolidate_chromatin_loops(
      files = c(f_empty, f1, f2, f3),
      mode = "consensus", gap = 0, quiet = TRUE
    ),
    "empty.*0 bytes"
  )
  meta <- attr(res, "looplook_metadata")
  expect_equal(meta$parameters$min_consensus_effective, 3L)
  expect_equal(meta$parameters$qc_status, "ok")
  expect_equal(length(res), 1L)
  unlink(c(f_empty, f1, f2, f3))
})


# ══════════════════════════════════════════════════════════════════════════
# 28. Effective consensus threshold: metadata matches actual retention (P1)
# ══════════════════════════════════════════════════════════════════════════

test_that("N=4: 3-rep cluster retained, 2-rep cluster dropped at effective min_consensus=3", {
  # Files 1,2,3 have chr1 loop → n_reps=3 for chr1 cluster
  # Files 3,4 have chr2 loop     → n_reps=2 for chr2 cluster
  f1 <- tempfile(fileext = ".bedpe")
  f2 <- tempfile(fileext = ".bedpe")
  f3 <- tempfile(fileext = ".bedpe")
  f4 <- tempfile(fileext = ".bedpe")
  writeLines("chr1\t100\t200\tchr1\t300\t400", f1)
  writeLines("chr1\t100\t200\tchr1\t300\t400", f2)
  writeLines("chr1\t100\t200\tchr1\t300\t400\nchr2\t100\t200\tchr2\t300\t400", f3)
  writeLines("chr2\t100\t200\tchr2\t300\t400", f4)
  res <- consolidate_chromatin_loops(
    files = c(f1, f2, f3, f4),
    mode = "consensus", gap = 0, quiet = TRUE
  )
  meta <- attr(res, "looplook_metadata")
  expect_equal(meta$parameters$min_consensus_effective, 3L)
  expect_equal(length(res), 1L,
    info = "chr1 cluster (3/4 reps) should be retained"
  )
  # Verify the retained loop is the chr1 cluster
  a1 <- InteractionSet::anchors(res, "first")
  expect_equal(as.character(GenomicRanges::seqnames(a1)), "chr1")
  unlink(c(f1, f2, f3, f4))
})


# ══════════════════════════════════════════════════════════════════════════
# 29. Union / intersect empty replicate cross-mode tests (P1)
# ══════════════════════════════════════════════════════════════════════════

test_that("union mode: one empty + one normal retains the normal loops", {
  f_empty <- tempfile(fileext = ".bedpe")
  f_norm <- tempfile(fileext = ".bedpe")
  file.create(f_empty)
  writeLines("chr1\t100\t200\tchr1\t300\t400", f_norm)
  expect_warning(
    res <- consolidate_chromatin_loops(
      files = c(f_empty, f_norm),
      mode = "union", gap = 0, quiet = TRUE
    ),
    "empty.*0 bytes"
  )
  meta <- attr(res, "looplook_metadata")
  expect_equal(meta$parameters$qc_status, "ok")
  expect_equal(length(res), 1L)
  expect_equal(meta$diagnostics$input_counts$raw_count, c(0L, 1L))
  expect_equal(S4Vectors::mcols(res)$n_reps, 1L)
  unlink(c(f_empty, f_norm))
})

test_that("union mode: two normal + one empty retains all normal loops", {
  f_empty <- tempfile(fileext = ".bedpe")
  f1 <- tempfile(fileext = ".bedpe")
  f2 <- tempfile(fileext = ".bedpe")
  file.create(f_empty)
  writeLines("chr1\t100\t200\tchr1\t300\t400", f1)
  writeLines("chr1\t100\t200\tchr1\t300\t500", f2)
  expect_warning(
    res <- consolidate_chromatin_loops(
      files = c(f_empty, f1, f2),
      mode = "union", gap = 0, quiet = TRUE
    ),
    "empty.*0 bytes"
  )
  # The two loops overlap at both anchors (chr1:100-200 matches,
  # chr1:300-400 overlaps chr1:300-500) → merged into one union cluster.
  meta <- attr(res, "looplook_metadata")
  expect_equal(meta$parameters$qc_status, "ok")
  expect_equal(meta$diagnostics$input_counts$raw_count, c(0L, 1L, 1L))
  expect_equal(length(res), 1L)
  expect_equal(S4Vectors::mcols(res)$n_reps, 2L)
  unlink(c(f_empty, f1, f2))
})

test_that("intersect mode: first file empty → empty_after_merge with warning", {
  f_empty <- tempfile(fileext = ".bedpe")
  f_norm <- tempfile(fileext = ".bedpe")
  file.create(f_empty)
  writeLines("chr1\t100\t200\tchr1\t300\t400", f_norm)
  expect_warning(
    res <- consolidate_chromatin_loops(
      files = c(f_empty, f_norm),
      mode = "intersect", gap = 0, quiet = TRUE
    ),
    "empty.*0 bytes"
  )
  meta <- attr(res, "looplook_metadata")
  expect_equal(meta$parameters$qc_status, "empty_after_merge")
  expect_equal(length(res), 0L)
  expect_equal(meta$diagnostics$input_counts$raw_count, c(0L, 1L))
  unlink(c(f_empty, f_norm))
})

test_that("intersect mode: non-reference file empty → empty_after_merge with warning", {
  f_norm <- tempfile(fileext = ".bedpe")
  f_empty <- tempfile(fileext = ".bedpe")
  writeLines("chr1\t100\t200\tchr1\t300\t400", f_norm)
  file.create(f_empty)
  expect_warning(
    res <- consolidate_chromatin_loops(
      files = c(f_norm, f_empty),
      mode = "intersect", gap = 0, quiet = TRUE
    ),
    "empty.*0 bytes"
  )
  meta <- attr(res, "looplook_metadata")
  expect_equal(meta$parameters$qc_status, "empty_after_merge")
  expect_equal(length(res), 0L)
  unlink(c(f_norm, f_empty))
})

test_that("intersect mode: zero-length result has correct mcols schema", {
  f1 <- tempfile(fileext = ".bedpe")
  f2 <- tempfile(fileext = ".bedpe")
  # Two non-overlapping loops → intersect produces 0 result
  writeLines("chr1\t100\t200\tchr1\t300\t400", f1)
  writeLines("chr2\t100\t200\tchr2\t300\t400", f2)
  res <- consolidate_chromatin_loops(
    files = c(f1, f2),
    mode = "intersect", gap = 0, quiet = TRUE
  )
  expect_equal(length(res), 0L)
  expect_true("n_reps" %in% colnames(S4Vectors::mcols(res)))
  expect_true("n_members" %in% colnames(S4Vectors::mcols(res)))
  expect_true("cluster_id" %in% colnames(S4Vectors::mcols(res)))
  unlink(c(f1, f2))
})


# ══════════════════════════════════════════════════════════════════════════
# 30. target_connected_loops: direct endpoints vs whole component (P1)
# ══════════════════════════════════════════════════════════════════════════

test_that(".build_target_hit_linkage returns only direct endpoint loops, not whole component", {
  skip_if_not_installed("GenomicRanges")
  # Component: A1-A2, A2-A3 (target peak hits A1 only, A3 is a branch)
  gr_anchors <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1", "chr1"),
    ranges = IRanges::IRanges(start = c(100, 200, 300), end = c(150, 250, 350)),
    anchor_id = c("A1", "A2", "A3"),
    cluster_id = c("c1", "c1", "c1")
  )
  names(gr_anchors) <- c("A1", "A2", "A3")
  hits <- GenomicRanges::findOverlaps(
    GenomicRanges::GRanges("chr1", IRanges::IRanges(100, 101)),
    gr_anchors
  )
  loop_df <- data.frame(
    loop_ID = c("L1", "L2"),
    a1_id = c("A1", "A2"), a2_id = c("A2", "A3"),
    cluster_id = c("c1", "c1"),
    stringsAsFactors = FALSE
  )
  result <- looplook:::.build_target_hit_linkage(
    hits = hits, gr_anchors = gr_anchors, anchor_topo_map = data.frame(anchor_id = character(), stringsAsFactors = FALSE),
    loop_annotation_final = loop_df, bed_info = data.frame(input_id = character(), stringsAsFactors = FALSE),
    map_info = data.frame(
      anchor_id = c("A1", "A2", "A3"),
      type_code = c("E", "E", "E"), SYMBOL = c("", "", ""),
      stringsAsFactors = FALSE
    ),
    ego_list_target = list()
  )
  tcl <- result$target_connected_loops
  expect_false(is.null(tcl), info = "Should have target-connected loops")
  # Only L1 (A1-A2) is direct; L2 (A2-A3) is a branch
  expect_equal(sort(tcl$loop_ID), "L1",
    info = "Only direct endpoint loops, not whole component branches"
  )
})

test_that(".build_target_hit_linkage same component: both endpoints hit → both loops included", {
  skip_if_not_installed("GenomicRanges")
  gr_anchors <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1", "chr1"),
    ranges = IRanges::IRanges(start = c(100, 200, 300), end = c(150, 250, 350)),
    anchor_id = c("A1", "A2", "A3"),
    cluster_id = c("c1", "c1", "c1")
  )
  names(gr_anchors) <- c("A1", "A2", "A3")
  # Target peaks hit BOTH A1 and A2
  hits <- GenomicRanges::findOverlaps(
    GenomicRanges::GRanges("chr1", IRanges::IRanges(c(100, 200), c(101, 201))),
    gr_anchors
  )
  loop_df <- data.frame(
    loop_ID = c("L1", "L2"),
    a1_id = c("A1", "A2"), a2_id = c("A2", "A3"),
    cluster_id = c("c1", "c1"),
    stringsAsFactors = FALSE
  )
  result <- looplook:::.build_target_hit_linkage(
    hits = hits, gr_anchors = gr_anchors, anchor_topo_map = data.frame(anchor_id = character(), stringsAsFactors = FALSE),
    loop_annotation_final = loop_df, bed_info = data.frame(input_id = character(), stringsAsFactors = FALSE),
    map_info = data.frame(
      anchor_id = c("A1", "A2", "A3"),
      type_code = c("E", "E", "E"), SYMBOL = c("", "", ""),
      stringsAsFactors = FALSE
    ),
    ego_list_target = list()
  )
  tcl <- result$target_connected_loops
  # Both A1 and A2 are hit, so L1 (A1-A2) is direct; L2 is only attached to A2
  expect_equal(sort(tcl$loop_ID), c("L1", "L2"),
    info = "Both incident loops should be included when both endpoints are hit"
  )
})


# ══════════════════════════════════════════════════════════════════════════
# 31. Expression case collision: core loader and profile both reject (P1)
# ══════════════════════════════════════════════════════════════════════════

test_that("load_expression_matrix rejects case-colliding gene IDs", {
  tmp <- tempfile(fileext = ".csv")
  writeLines("Gene\tS1\nGENE1\t10\ngene1\t20", tmp)
  expect_error(
    looplook:::load_expression_matrix(tmp, "S1"),
    "collide after case-insensitive matching"
  )
  unlink(tmp)
})


# ══════════════════════════════════════════════════════════════════════════
# 32. Refinement donut uses Linked_Loop_IDs from main annotation (P2)
# ══════════════════════════════════════════════════════════════════════════

test_that(".build_refinement_donut uses Linked_Loop_IDs, not raw re-overlap", {
  skip_if_not_installed("ggplot2")
  loop_df <- data.frame(
    loop_ID = c("L1", "L2"), loop_type = c("E-P", "P-P"),
    a1_id = c("A1", "A2"), a2_id = c("A2", "A3"),
    chr1 = c("chr1", "chr1"), start1 = c(100L, 200L), end1 = c(150L, 250L),
    chr2 = c("chr1", "chr1"), start2 = c(200L, 300L), end2 = c(250L, 350L),
    cluster_id = c("c1", "c1"),
    stringsAsFactors = FALSE
  )
  # bed_info with Linked_Loop_IDs from main annotation pipeline:
  # only L1 is linked to the target peak, L2 is a component branch
  bed_info <- data.frame(
    input_id = "Peak_1",
    Linked_Loop_IDs = "L1",
    stringsAsFactors = FALSE
  )
  p <- looplook:::.build_refinement_donut(
    bed_info = bed_info, loop_df = loop_df,
    custom_colors = c("E-P" = "orange", "P-P" = "blue"),
    project_name = "Test"
  )
  expect_false(is.null(p),
    info = "Should produce a donut for accepted Linked_Loop_IDs"
  )
})

test_that(".build_refinement_donut returns NULL when Linked_Loop_IDs is empty", {
  skip_if_not_installed("ggplot2")
  loop_df <- data.frame(
    loop_ID = "L1", loop_type = "E-P",
    a1_id = "A1", a2_id = "A2",
    cluster_id = "c1",
    stringsAsFactors = FALSE
  )
  bed_info <- data.frame(
    input_id = "Peak_1",
    Linked_Loop_IDs = NA_character_,
    stringsAsFactors = FALSE
  )
  p <- looplook:::.build_refinement_donut(
    bed_info = bed_info, loop_df = loop_df,
    custom_colors = c("E-P" = "orange"),
    project_name = "Test"
  )
  expect_null(p,
    info = "Should return NULL when no Linked_Loop_IDs"
  )
})


# ══════════════════════════════════════════════════════════════════════════
# 22. Promoter activity: unmeasured / silent / active (looplook68-69 fix)
# ══════════════════════════════════════════════════════════════════════════

test_that(".build_promoter_centric_df distinguishes unmeasured, silent, active", {
  raw <- data.frame(
    Gene = c("G_UNMEASURED", "G_SILENT", "G_ACTIVE"),
    Total_Loops_Filtered = c(1L, 1L, 1L),
    n_Unique_Contacts_Filtered = c(1L, 1L, 1L),
    n_Linked_Promoters_Filtered = c(0L, 0L, 0L),
    n_Linked_Distal_Filtered = c(1L, 1L, 1L),
    n_Linked_EnhancerLike_Filtered = c(1L, 1L, 1L),
    Dominant_Interaction_Filtered = c("E-P", "E-P", "E-P"),
    stringsAsFactors = FALSE
  )
  vals <- c(G_SILENT = 0, G_ACTIVE = 10)
  out <- looplook:::.build_promoter_centric_df(
    raw_stats_df = raw, upstream_promoter_stats = NULL,
    vals = vals, threshold = 1, hub_percentile = 0.95
  )
  state <- setNames(out$Is_Active_Gene, out$Gene)
  expect_equal(state[["G_UNMEASURED"]], "Not_assessed",
    info = "Gene not in expression matrix must be Not_assessed"
  )
  expect_equal(state[["G_SILENT"]], "No",
    info = "Gene with low expression must be No"
  )
  expect_equal(state[["G_ACTIVE"]], "Yes",
    info = "Gene above threshold must be Yes"
  )
})

test_that(".build_promoter_centric_df upstream branch also distinguishes three states", {
  raw <- data.frame(
    Gene = c("G_UNMEASURED", "G_SILENT", "G_ACTIVE"),
    Total_Loops_Filtered = c(1L, 1L, 1L),
    n_Unique_Contacts_Filtered = c(1L, 1L, 1L),
    n_Linked_Promoters_Filtered = c(0L, 0L, 0L),
    n_Linked_Distal_Filtered = c(1L, 1L, 1L),
    n_Linked_EnhancerLike_Filtered = c(1L, 1L, 1L),
    Dominant_Interaction_Filtered = c("E-P", "E-P", "E-P"),
    stringsAsFactors = FALSE
  )
  upstream <- data.frame(
    Gene = c("G_UNMEASURED", "G_SILENT", "G_ACTIVE"),
    Total_Loops = c(2L, 2L, 2L), n_Linked_Promoters = c(0L, 0L, 0L),
    n_Linked_Distal = c(2L, 2L, 2L),
    Dominant_Interaction = c("E-P", "E-P", "E-P"),
    Is_High_Connectivity_Gene = "No", Is_High_Distal_Connectivity_Gene = "No",
    Is_Active_Gene = "Not_assessed", stringsAsFactors = FALSE
  )
  vals <- c(G_SILENT = 0, G_ACTIVE = 10)
  out <- looplook:::.build_promoter_centric_df(
    raw_stats_df = raw, upstream_promoter_stats = upstream,
    vals = vals, threshold = 1, hub_percentile = 0.95
  )
  state <- setNames(out$Is_Active_Gene, out$Gene)
  expect_equal(state[["G_UNMEASURED"]], "Not_assessed")
  expect_equal(state[["G_SILENT"]], "No")
  expect_equal(state[["G_ACTIVE"]], "Yes")
})

# ══════════════════════════════════════════════════════════════════════════
# Chromatin BED boundary tests (required/optional empty, shared seqlevel,
# zero overlap, invalid path)
# ══════════════════════════════════════════════════════════════════════════

test_that("optional empty chromatin BED is non-fatal (FALSE/not-called)", {
  anchors <- data.frame(
    anchor_id = "A1", chr = "chr1", start = 1L, end = 100L,
    anchor_type = "E", anchor_gene = "", cluster_id = "C1",
    stringsAsFactors = FALSE
  )
  empty_bed <- tempfile(fileext = ".bed")
  file.create(empty_bed)
  known_marks <- c("H3K4me1", "H3K4me3", "H3K27ac")
  res <- looplook:::.overlap_chromatin_marks(
    anchors, chromatin_beds = list(H3K27ac = empty_bed),
    provided_marks = "H3K27ac", known_marks = known_marks,
    anchor_gap = 200L, anchor_min_overlap = 1L,
    log_message = function(...) {}
  )
  expect_true(all(res$mark_matrix$H3K27ac == FALSE))
  expect_equal(res$valid_provided_marks, "H3K27ac")
  unlink(empty_bed)
})

test_that("required empty chromatin BED is fatal", {
  anchors <- data.frame(
    anchor_id = "A1", chr = "chr1", start = 1L, end = 100L,
    anchor_type = "E", anchor_gene = "", cluster_id = "C1",
    stringsAsFactors = FALSE
  )
  empty_bed <- tempfile(fileext = ".bed")
  file.create(empty_bed)
  known_marks <- c("H3K4me1", "H3K4me3")
  expect_error(
    looplook:::.overlap_chromatin_marks(
      anchors, chromatin_beds = list(H3K4me1 = empty_bed),
      provided_marks = "H3K4me1", known_marks = known_marks,
      anchor_gap = 200L, anchor_min_overlap = 1L,
      log_message = function(...) {}
    ),
    "zero peaks"
  )
  unlink(empty_bed)
})

test_that("non-empty required BED with shared seqlevels but zero overlap warns", {
  anchors <- data.frame(
    anchor_id = "A1", chr = "chr1", start = 1L, end = 100L,
    anchor_type = "E", anchor_gene = "", cluster_id = "C1",
    stringsAsFactors = FALSE
  )
  far_bed <- tempfile(fileext = ".bed")
  writeLines("chr1\t1000000\t1000100", far_bed)
  known_marks <- c("H3K4me1", "H3K4me3")
  expect_warning(
    res <- looplook:::.overlap_chromatin_marks(
      anchors, chromatin_beds = list(H3K4me1 = far_bed),
      provided_marks = "H3K4me1", known_marks = known_marks,
      anchor_gap = 200L, anchor_min_overlap = 1L,
      log_message = function(...) {}
    ),
    "tested-but-not-called"
  )
  expect_true(all(res$mark_matrix$H3K4me1 == FALSE))
  unlink(far_bed)
})

test_that("non-empty BED with no shared seqlevels is fatal", {
  anchors <- data.frame(
    anchor_id = "A1", chr = "chr1", start = 1L, end = 100L,
    anchor_type = "E", anchor_gene = "", cluster_id = "C1",
    stringsAsFactors = FALSE
  )
  mismatch_bed <- tempfile(fileext = ".bed")
  writeLines("chrUn_gl000220\t1000\t2000", mismatch_bed)
  known_marks <- c("H3K4me1", "H3K4me3")
  expect_error(
    looplook:::.overlap_chromatin_marks(
      anchors, chromatin_beds = list(H3K4me1 = mismatch_bed),
      provided_marks = "H3K4me1", known_marks = known_marks,
      anchor_gap = 200L, anchor_min_overlap = 1L,
      log_message = function(...) {}
    ),
    "no seqlevels in common"
  )
  unlink(mismatch_bed)
})

test_that("invalid optional BED path warns and skips", {
  anchors <- data.frame(
    anchor_id = "A1", chr = "chr1", start = 1L, end = 100L,
    anchor_type = "E", anchor_gene = "", cluster_id = "C1",
    stringsAsFactors = FALSE
  )
  known_marks <- c("H3K4me1", "H3K4me3")
  expect_warning(
    res <- looplook:::.overlap_chromatin_marks(
      anchors, chromatin_beds = list(H3K27ac = "/no/such/file.bed"),
      provided_marks = "H3K27ac", known_marks = known_marks,
      anchor_gap = 200L, anchor_min_overlap = 1L,
      log_message = function(...) {}
    ),
    "not found"
  )
})

# ══════════════════════════════════════════════════════════════════════════
# End-to-end neighbor_hop = 1 test (2-hop target expansion)
# ══════════════════════════════════════════════════════════════════════════

test_that("annotate_peaks_and_loops with neighbor_hop=1 validates 2-hop results", {
  skip_if_not_installed("org.Hs.eg.db")
  skip_if(is.null(get_test_txdb()), "Sample TxDb not available")
  txdb_obj <- get_test_txdb()

  # A1(100-600) —L1— A2(1100-1600) —L2— A3(10415000-10415600,gene)
  #                                 \_L3— A4(200000-200600,no gene)
  # neighbor_hop=1 → reaches A3 (2-hop) but not A4 (3-hop)
  tmp_bedpe <- tempfile(fileext = ".bedpe")
  writeLines(c(
    "chr6\t100\t600\tchr6\t1100\t1600",               # L1
    "chr6\t1100\t1600\tchr6\t10415000\t10415600",     # L2: A3 has TFAP2A-AS1
    "chr6\t10415000\t10415600\tchr6\t200000\t200600"  # L3: A4 no gene
  ), tmp_bedpe)
  tmp_target <- tempfile(fileext = ".bed")
  writeLines("chr6\t100\t600", tmp_target)

  res <- annotate_peaks_and_loops(
    bedpe_file = tmp_bedpe, target_bed = tmp_target,
    txdb = txdb_obj, org_db = "org.Hs.eg.db", species = "hg19",
    neighbor_hop = 1,
    out_dir = tempdir(), write_output = FALSE, quiet = TRUE
  )

  links <- res$target_gene_links
  expect_true(!is.null(links) && nrow(links) > 0)

  # Primary targets must have path_length <= 1
  primary <- links[links$anchor_role %in% c("local_anchor", "opposite_anchor"), ]
  expect_true(all(primary$path_length <= 1L))

  # 2-hop expanded targets should exist (A3 is 2 hops from A1)
  expanded <- links[links$anchor_role == "expanded_anchor", ]
  expect_gt(nrow(expanded), 0L)
  expect_equal(unique(expanded$path_length), 2L,
    info = "expanded targets must have path_length 2")

  # Primary targets invariant between hop=0 and hop=1
  res0 <- annotate_peaks_and_loops(
    bedpe_file = tmp_bedpe, target_bed = tmp_target,
    txdb = txdb_obj, org_db = "org.Hs.eg.db", species = "hg19",
    neighbor_hop = 0,
    out_dir = tempdir(), write_output = FALSE, quiet = TRUE
  )
  expect_equal(res0$target_annotation$Assigned_Target_Genes,
    res$target_annotation$Assigned_Target_Genes)

  unlink(c(tmp_bedpe, tmp_target))
})

# ══════════════════════════════════════════════════════════════════════════
# Lazy ego vs full ego equivalence test
# ══════════════════════════════════════════════════════════════════════════

test_that("lazy target ego produces same links as full ego", {
  skip_if_not_installed("igraph")
  # Chain: A1–A2–A3 (2 edges, no ring).  A3 is promoter with GENE3.
  # distance(A1, A3) = 2 → GENE3 is expanded_anchor.
  g <- igraph::graph_from_edgelist(
    matrix(c("A1","A2", "A2","A3"), ncol = 2, byrow = TRUE),
    directed = FALSE
  )
  igraph::E(g)$n_support <- c(1, 1)

  hit_df <- data.frame(qid = 1L, sid = 1L, anchor_id = "A1",
    stringsAsFactors = FALSE)
  bed_info <- data.frame(input_id = "Peak_1", stringsAsFactors = FALSE)
  la <- data.frame(
    loop_ID = c("L1", "L2"),
    a1_id = c("A1", "A2"),
    a2_id = c("A2", "A3"),
    stringsAsFactors = FALSE
  )
  map_info <- data.frame(
    anchor_id = c("A1", "A2", "A3"),
    type_code = c("E", "E", "P"),
    SYMBOL = c("", "", "GENE3"),
    stringsAsFactors = FALSE
  )

  # Full ego
  nodes_in_graph <- igraph::V(g)$name
  full_ego <- igraph::ego(g, order = 2L, nodes = nodes_in_graph, mode = "all")
  names(full_ego) <- nodes_in_graph

  links_full <- looplook:::.build_target_gene_links(
    hit_df = hit_df, bed_info = bed_info,
    loop_annotation_final = la, map_info = map_info,
    ego_list_target = full_ego, g = g, neighbor_hop = 1L
  )

  # Lazy ego (empty cache, computed on demand)
  links_lazy <- looplook:::.build_target_gene_links(
    hit_df = hit_df, bed_info = bed_info,
    loop_annotation_final = la, map_info = map_info,
    ego_list_target = list(), g = g, neighbor_hop = 1L
  )

  # GENE3 must be 2-hop expanded in both
  g3_full <- links_full[links_full$gene == "GENE3" &
    links_full$anchor_role == "expanded_anchor", ]
  g3_lazy <- links_lazy[links_lazy$gene == "GENE3" &
    links_lazy$anchor_role == "expanded_anchor", ]
  expect_equal(nrow(g3_full), 1L)
  expect_equal(nrow(g3_lazy), 1L)
  expect_equal(g3_full$path_length, 2L)
  expect_equal(g3_lazy$path_length, 2L)

  cols <- c("input_id", "anchor_id", "gene", "path_length",
    "anchor_role", "evidence")
  expect_equal(
    links_full[order(links_full$gene), cols],
    links_lazy[order(links_lazy$gene), cols],
    ignore_attr = TRUE
  )
})

# ══════════════════════════════════════════════════════════════════════════
# Pair broadcast test: multiple peaks hitting same anchor
# ══════════════════════════════════════════════════════════════════════════

test_that("pair dedup + join preserves per-peak expanded rows", {
  skip_if_not_installed("igraph")
  # Chain: A1–A2–A3.  Two peaks both hit A1 → both see GENE3.
  g <- igraph::graph_from_edgelist(
    matrix(c("A1","A2", "A2","A3"), ncol = 2, byrow = TRUE),
    directed = FALSE
  )
  igraph::E(g)$n_support <- c(1, 1)

  hit_df <- data.frame(
    qid = c(1L, 2L), sid = c(1L, 1L),
    anchor_id = c("A1", "A1"),
    stringsAsFactors = FALSE
  )
  bed_info <- data.frame(
    input_id = c("Peak_1", "Peak_2"),
    stringsAsFactors = FALSE
  )
  la <- data.frame(
    loop_ID = c("L1", "L2"),
    a1_id = c("A1", "A2"),
    a2_id = c("A2", "A3"),
    stringsAsFactors = FALSE
  )
  map_info <- data.frame(
    anchor_id = c("A1", "A2", "A3"),
    type_code = c("E", "E", "P"),
    SYMBOL = c("", "", "GENE3"),
    stringsAsFactors = FALSE
  )

  links <- looplook:::.build_target_gene_links(
    hit_df = hit_df, bed_info = bed_info,
    loop_annotation_final = la, map_info = map_info,
    ego_list_target = list(), g = g, neighbor_hop = 1L
  )

  g3 <- links[links$gene == "GENE3" &
    links$anchor_role == "expanded_anchor", ]
  expect_equal(sort(g3$input_id), c("Peak_1", "Peak_2"),
    info = "both peaks must see the 2-hop gene GENE3 as expanded")
  expect_true(all(g3$path_length == 2L))
})

test_that("neighbor_hop > 1 is rejected", {
  expect_error(
    annotate_peaks_and_loops(bedpe_file = "no", neighbor_hop = 2),
    "not supported"
  )
})

# ══════════════════════════════════════════════════════════════════════════
# .resolve_chromatin_gene_role() matrix tests
# ══════════════════════════════════════════════════════════════════════════

test_that("resolver: P/P retains promoter strict", {
  r <- looplook:::.resolve_chromatin_gene_role("P", "P")
  expect_equal(r$role, "promoter")
  expect_equal(r$strict, TRUE)
})

test_that("resolver: P/E becomes enhancer_candidate strict", {
  r <- looplook:::.resolve_chromatin_gene_role("P", "E")
  expect_equal(r$role, "enhancer_candidate")
  expect_equal(r$strict, TRUE)
})

test_that("resolver: G/P host gene becomes promoter strict", {
  r <- looplook:::.resolve_chromatin_gene_role("G", "P")
  expect_equal(r$role, "promoter")
  expect_equal(r$strict, TRUE)
})

test_that("resolver: G/E becomes enhancer_candidate strict", {
  r <- looplook:::.resolve_chromatin_gene_role("G", "E")
  expect_equal(r$role, "enhancer_candidate")
  expect_equal(r$strict, TRUE)
})

test_that("resolver: E/E becomes enhancer_candidate NOT strict", {
  r <- looplook:::.resolve_chromatin_gene_role("E", "E")
  expect_equal(r$role, "enhancer_candidate")
  expect_equal(r$strict, FALSE)
})

test_that("resolver: E/P without TSS becomes positional_candidate NOT strict", {
  r <- looplook:::.resolve_chromatin_gene_role("E", "P", tss_supported = FALSE)
  expect_equal(r$role, "positional_candidate")
  expect_equal(r$strict, FALSE)
})

test_that("resolver: E/P with TSS becomes promoter strict", {
  r <- looplook:::.resolve_chromatin_gene_role("E", "P", tss_supported = TRUE)
  expect_equal(r$role, "promoter")
  expect_equal(r$strict, TRUE)
})

test_that("resolver: E/dual with TSS becomes promoter strict", {
  r <- looplook:::.resolve_chromatin_gene_role("E", "dual", tss_supported = TRUE)
  expect_equal(r$role, "promoter")
  expect_equal(r$strict, TRUE)
})

test_that("resolver: has_gene=FALSE forces strict=FALSE regardless of type", {
  r <- looplook:::.resolve_chromatin_gene_role("P", "P", has_gene = FALSE)
  expect_equal(r$role, "promoter")
  expect_equal(r$strict, FALSE)
})

test_that("resolver: eP/eG treated same as P/G for role", {
  for (tt in c("eP", "eG")) {
    r <- looplook:::.resolve_chromatin_gene_role(tt, tt)
    expect_true(r$role %in% c("promoter", "gene_body"))
    expect_true(r$strict)
  }
})

test_that("resolver: vectorised form works correctly", {
  old <- c("P", "G", "E", "E", "P")
  final <- c("E", "P", "E", "P", "P")
  tss <- c(FALSE, FALSE, FALSE, TRUE, TRUE)
  r <- looplook:::.resolve_chromatin_gene_role(old, final, tss_supported = tss)
  expect_equal(r$role, c("enhancer_candidate", "promoter",
    "enhancer_candidate", "promoter", "promoter"))
  expect_equal(r$strict, c(TRUE, TRUE, FALSE, TRUE, TRUE))
})

# ══════════════════════════════════════════════════════════════════════════
# Layer consistency: anchor → loop → target must share role/strict
# ══════════════════════════════════════════════════════════════════════════

test_that("chromatin refinement propagates resolver result across layers", {
  skip_if_not_installed("org.Hs.eg.db")
  skip_if(is.null(get_test_txdb()), "Sample TxDb not available")
  txdb_obj <- get_test_txdb()

  tmp_bedpe <- tempfile(fileext = ".bedpe")
  writeLines("chr6\t10412000\t10412600\tchr6\t10415000\t10415600", tmp_bedpe)
  tmp_target <- tempfile(fileext = ".bed")
  writeLines("chr6\t10412000\t10412600", tmp_target)

  res <- annotate_peaks_and_loops(
    bedpe_file = tmp_bedpe, target_bed = tmp_target,
    txdb = txdb_obj,
    org_db = "org.Hs.eg.db", species = "hg19",
    out_dir = tempdir(), write_output = FALSE, quiet = TRUE
  )

  h3k4me1 <- tempfile(fileext = ".bed")
  h3k4me3 <- tempfile(fileext = ".bed")
  writeLines("chr6\t10411900\t10412700", h3k4me1)
  writeLines("chr6\t10411900\t10412700", h3k4me3)

  cr <- refine_loop_anchors_by_chromatin(res,
    chromatin_beds = list(H3K4me1 = h3k4me1, H3K4me3 = h3k4me3),
    recompute_targets = TRUE, write_output = FALSE, quiet = TRUE
  )

  st <- attr(cr, "looplook_anchor_state")
  mi <- st$map_info
  la <- cr$loop_annotation

  # Structural assertions: required columns must exist
  expect_true(all(c("anchor_id", "effective_gene_role",
    "strict_assignment_eligible") %in% colnames(mi)))
  expect_true(all(c("a1_id", "a2_id", "anchor1_gene_role",
    "anchor2_gene_role", "anchor1_strict_eligible",
    "anchor2_strict_eligible") %in% colnames(la)))

  role_map <- setNames(mi$effective_gene_role, mi$anchor_id)
  strict_map <- setNames(mi$strict_assignment_eligible, mi$anchor_id)

  # a1/a2 role/strict match between map_info and loop_annotation
  expect_equal(unname(role_map[la$a1_id]), la$anchor1_gene_role)
  expect_equal(unname(strict_map[la$a1_id]), la$anchor1_strict_eligible)
  expect_equal(unname(role_map[la$a2_id]), la$anchor2_gene_role)
  expect_equal(unname(strict_map[la$a2_id]), la$anchor2_strict_eligible)

  # Check target_gene_links also shares the same resolver output
  links <- cr$target_gene_links
  expect_true(!is.null(links) && nrow(links) > 0,
    info = "target_gene_links must be non-empty when target_bed is provided")
  # Only compare loop-anchor links (exclude linear_annotation rows with NA anchor_id)
  loop_links <- links[links$source == "loop_anchor" &
    !is.na(links$anchor_id) & links$anchor_id != "", ]
  expect_gt(nrow(loop_links), 0)
  expect_equal(unname(role_map[loop_links$anchor_id]), loop_links$gene_role)
  expect_equal(unname(strict_map[loop_links$anchor_id]),
    loop_links$strict_assignment_eligible)

  unlink(c(tmp_bedpe, tmp_target, h3k4me1, h3k4me3))
})

# ══════════════════════════════════════════════════════════════════════════
# Old-schema fallback tests: enhancer_candidate without strict column
# ══════════════════════════════════════════════════════════════════════════

test_that("old-schema enhancer_candidate without strict column defaults to FALSE", {
  mi <- data.frame(
    anchor_id = c("A1", "A2"),
    SYMBOL = c("SOX2", "MYC"),
    type_code = c("E", "E"),
    effective_gene_role = c("enhancer_candidate", "enhancer_candidate"),
    stringsAsFactors = FALSE
  )
  # No strict_assignment_eligible column — simulates pre-87 RData
  gm <- looplook:::.target_anchor_gene_map(mi)
  expect_true(all(gm$strict_assignment_eligible == FALSE))
  expect_true(all(gm$gene_role == "enhancer_candidate"))
})

test_that("old-schema promoter/gene_body without strict column defaults to TRUE", {
  mi <- data.frame(
    anchor_id = c("A1", "A2"),
    SYMBOL = c("SOX2", "MYC"),
    type_code = c("P", "G"),
    effective_gene_role = c("promoter", "gene_body"),
    stringsAsFactors = FALSE
  )
  gm <- looplook:::.target_anchor_gene_map(mi)
  expect_true(all(gm$strict_assignment_eligible == TRUE))
})

test_that("old-schema enhancer_candidate with strict=NA resolved to FALSE in aggregator", {
  links <- data.frame(
    input_id = "Peak_1",
    loop_ID = "L1",
    anchor_id = "A1",
    gene = "GENE1",
    gene_role = "enhancer_candidate",
    source = "loop_anchor",
    evidence = "test",
    anchor_role = "local_anchor",
    path_length = 0L,
    strict_assignment_eligible = NA,
    stringsAsFactors = FALSE
  )
  bed_info <- data.frame(input_id = "Peak_1", stringsAsFactors = FALSE)
  result <- looplook:::.aggregate_strict_targets(links, bed_info)
  # Enhancer_candidate with NA strict should NOT enter strict target summary
  expect_true(is.na(result$Assigned_Target_Genes[1]) ||
    result$Assigned_Target_Genes[1] == "")
})

# ══════════════════════════════════════════════════════════════════════════
# Direct P–P co-assignment tests
# ══════════════════════════════════════════════════════════════════════════

.make_pp_links <- function(genes, paths, anchors, loop = "L1", input = "P1",
                            roles = rep("promoter", length(genes)),
                            strict = rep(TRUE, length(genes))) {
  data.frame(
    input_id = input, loop_ID = loop, anchor_id = anchors,
    gene = genes, gene_role = roles,
    strict_assignment_eligible = strict,
    path_length = paths,
    source = "loop_anchor", evidence = "test",
    anchor_role = ifelse(paths == 0, "local_anchor", "opposite_anchor"),
    stringsAsFactors = FALSE
  )
}

test_that("P-P: direct promoter pair co-assigned (promoter_then_distance)", {
  links <- .make_pp_links(c("A","B"), c(0L,1L), c("X","Y"))
  bed <- data.frame(input_id = "P1", stringsAsFactors = FALSE)
  r <- looplook:::.aggregate_strict_targets(links, bed,
    target_priority = "promoter_then_distance")
  expect_equal(r$Assigned_Target_Genes[1], "A;B")
  expect_equal(r$Regulated_promoter_genes[1], "A;B")
})

test_that("P-P: distance_then_role keeps local-only", {
  links <- .make_pp_links(c("A","B"), c(0L,1L), c("X","Y"))
  bed <- data.frame(input_id = "P1", stringsAsFactors = FALSE)
  r <- looplook:::.aggregate_strict_targets(links, bed,
    target_priority = "distance_then_role")
  expect_equal(r$Assigned_Target_Genes[1], "A")
  expect_equal(r$Regulated_promoter_genes[1], "A;B")
})

test_that("P-P: different loops are not paired", {
  links <- rbind(
    .make_pp_links("A", 0L, "X", loop = "L1"),
    .make_pp_links("B", 1L, "Y", loop = "L2")
  )
  bed <- data.frame(input_id = "P1", stringsAsFactors = FALSE)
  r <- looplook:::.aggregate_strict_targets(links, bed,
    target_priority = "promoter_then_distance")
  expect_equal(r$Assigned_Target_Genes[1], "A")
})

test_that("P-P: same anchor not counted as pair", {
  links <- .make_pp_links(c("A","A"), c(0L,1L), c("X","X"))
  bed <- data.frame(input_id = "P1", stringsAsFactors = FALSE)
  r <- looplook:::.aggregate_strict_targets(links, bed,
    target_priority = "promoter_then_distance")
  expect_equal(r$Assigned_Target_Genes[1], "A")
})

test_that("P-P: one end not strict → not co-assigned", {
  links <- .make_pp_links(c("A","B"), c(0L,1L), c("X","Y"),
    strict = c(TRUE, FALSE))
  bed <- data.frame(input_id = "P1", stringsAsFactors = FALSE)
  r <- looplook:::.aggregate_strict_targets(links, bed,
    target_priority = "promoter_then_distance")
  expect_equal(r$Assigned_Target_Genes[1], "A")
})

test_that("P-P: positional_candidate cannot co-assign", {
  links <- .make_pp_links(c("A","B"), c(0L,1L), c("X","Y"),
    roles = c("promoter", "positional_candidate"),
    strict = c(TRUE, FALSE))
  bed <- data.frame(input_id = "P1", stringsAsFactors = FALSE)
  r <- looplook:::.aggregate_strict_targets(links, bed,
    target_priority = "promoter_then_distance")
  expect_equal(r$Assigned_Target_Genes[1], "A")
})

test_that("P-P: multi-loop union preserves all genes", {
  links <- rbind(
    .make_pp_links(c("A","B"), c(0L,1L), c("X","Y"), loop = "L1"),
    .make_pp_links(c("A","C"), c(0L,1L), c("X","Z"), loop = "L2")
  )
  bed <- data.frame(input_id = "P1", stringsAsFactors = FALSE)
  r <- looplook:::.aggregate_strict_targets(links, bed,
    target_priority = "promoter_then_distance")
  expect_equal(r$Assigned_Target_Genes[1], "A;B;C")
})

test_that("P-P: same gene on both anchors deduplicated", {
  links <- .make_pp_links(c("A","A"), c(0L,1L), c("X","Y"))
  bed <- data.frame(input_id = "P1", stringsAsFactors = FALSE)
  r <- looplook:::.aggregate_strict_targets(links, bed,
    target_priority = "promoter_then_distance")
  expect_equal(r$Assigned_Target_Genes[1], "A")
})

test_that("P-P: path2 promoter excluded", {
  links <- rbind(
    .make_pp_links("A", 0L, "X", loop = "L1"),
    .make_pp_links("B", 2L, "Y", loop = "L1")
  )
  bed <- data.frame(input_id = "P1", stringsAsFactors = FALSE)
  r <- looplook:::.aggregate_strict_targets(links, bed,
    target_priority = "promoter_then_distance")
  expect_equal(r$Assigned_Target_Genes[1], "A")
})

test_that("P-P: union preserves already-assigned genes from other loops", {
  # A,C are each path0 on different loops → normal priority assigns A;C
  # B is path1 on L1 with A → P-P co-assigns B
  # Result: A;B;C (union)
  links <- rbind(
    .make_pp_links(c("A","B"), c(0L,1L), c("X","Y"), loop = "L1"),
    .make_pp_links("C", 0L, "Z", loop = "L2")
  )
  bed <- data.frame(input_id = "P1", stringsAsFactors = FALSE)
  r <- looplook:::.aggregate_strict_targets(links, bed,
    target_priority = "promoter_then_distance")
  expect_equal(r$Assigned_Target_Genes[1], "A;B;C")
})

test_that("P-P: empty anchor_id does not trigger co-assignment", {
  links <- .make_pp_links(c("A","B"), c(0L,1L), c("X", NA_character_),
    loop = "L1", input = "P1", strict = c(TRUE, TRUE))
  links$anchor_id[2] <- NA
  bed <- data.frame(input_id = "P1", stringsAsFactors = FALSE)
  r <- looplook:::.aggregate_strict_targets(links, bed,
    target_priority = "promoter_then_distance")
  expect_equal(r$Assigned_Target_Genes[1], "A")
})

test_that("P-P: different inputs not merged", {
  links <- rbind(
    .make_pp_links(c("A","B"), c(0L,1L), c("X","Y"), loop = "L1", input = "P1"),
    .make_pp_links("C", 0L, "Z", loop = "L1", input = "P2")
  )
  bed <- data.frame(input_id = c("P1","P2"), stringsAsFactors = FALSE)
  r <- looplook:::.aggregate_strict_targets(links, bed,
    target_priority = "promoter_then_distance")
  expect_equal(r$Assigned_Target_Genes[1], "A;B")
  expect_equal(r$Assigned_Target_Genes[2], "C")
})

test_that("P-P: empty string anchor_id does not trigger co-assignment", {
  links <- .make_pp_links(c("A","B"), c(0L,1L), c("X","Y"))
  links$anchor_id[2] <- ""
  bed <- data.frame(input_id = "P1", stringsAsFactors = FALSE)
  r <- looplook:::.aggregate_strict_targets(links, bed,
    target_priority = "promoter_then_distance")
  expect_equal(r$Assigned_Target_Genes[1], "A")
})

test_that("P-P: max_primary_hop=0 excludes path1 from primary", {
  links <- .make_pp_links(c("A","B"), c(0L,1L), c("X","Y"))
  bed <- data.frame(input_id = "P1", stringsAsFactors = FALSE)
  r <- looplook:::.aggregate_strict_targets(links, bed,
    target_priority = "promoter_then_distance", max_primary_hop = 0L)
  expect_equal(r$Assigned_Target_Genes[1], "A")
})

test_that("dual is rejected as explicit chromatin candidate_type", {
  expect_error(
    looplook:::.validate_chromatin_overlap_inputs(
      chromatin_beds = list(H3K4me1 = "x.bed"),
      candidate_types = "dual",
      anchor_gap = 200,
      anchor_min_overlap = 100
    ),
    "Unknown candidate_type"
  )
})
