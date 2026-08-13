# tests/testthat/test-annotation.R
# Shared annotation bases (TxDb + tiny example) are provided by
# helper-fixtures.R so they are built once per session across test files.

test_that("packaged annotation example keeps the expected output contract", {
  rdata_path <- system.file("extdata", "analysis_results.RData", package = "looplook")
  skip_if(rdata_path == "", "Pre-computed test data not available")

  temp_env <- new.env()
  load(rdata_path, envir = temp_env)
  res_integrated <- temp_env[[ls(temp_env)[1]]]

  expect_type(res_integrated, "list")
  expect_true(all(c("target_annotation", "loop_annotation", "anchor_loci_annotation", "anchor_annotation") %in% names(res_integrated)))
  # The packaged fixture intentionally excludes the large pre-computed
  # `plots`/`plot_list` objects (they are re-generated on demand by the
  # pipeline functions); only the annotation tables are retained.
  expect_false("plots" %in% names(res_integrated))
  expect_gt(nrow(res_integrated$loop_annotation), 0)
  expect_gt(nrow(res_integrated$target_annotation), 0)
  expect_true("Assigned_Target_Genes_Filled" %in% colnames(res_integrated$target_annotation))
})

test_that("annotate_peaks_and_loops shared setup: flags + anchor_gap proximity", {
  skip_if_not_installed("org.Hs.eg.db")
  skip_if(is.null(get_test_txdb()), "Sample TxDb not available")
  txdb_obj <- get_test_txdb()

  out_base <- tempfile(pattern = "looplook_anno_nowrite_")
  unlink(out_base, recursive = TRUE, force = TRUE)
  expect_false(dir.exists(out_base))

  # Test 1: quiet + write_output = FALSE
  res <- suppressPackageStartupMessages(
    looplook:::.with_known_upstream_noise_suppressed(
      looplook::annotate_peaks_and_loops(
        bedpe_file = system.file("extdata", "example_loops_1.bedpe", package = "looplook"),
        txdb = txdb_obj, org_db = "org.Hs.eg.db", species = "hg19",
        out_dir = out_base, write_output = FALSE, quiet = TRUE
      )
    )
  )
  expect_type(res, "list")
  expect_false(dir.exists(out_base), info = "write_output=FALSE must not create directory")

  # Test 2: anchor_gap proximity linking
  target_bed <- tempfile(fileext = ".bed")
  loop_bedpe <- tempfile(fileext = ".bedpe")
  # Peak 200bp away from anchor, 0bp actual overlap
  writeLines("chr6\t10412800\t10413000", target_bed)
  writeLines("chr6\t10412000\t10412600\tchr6\t10415000\t10415600", loop_bedpe)

  # anchor_gap=200 + anchor_min_overlap=1: 0bp actual overlap → should NOT link
  # anchor_min_overlap >= 1L filter is always applied: proximity-only without
  # physical overlap is excluded.
  res_gap <- annotate_peaks_and_loops(
    bedpe_file = loop_bedpe, target_bed = target_bed,
    txdb = txdb_obj, org_db = "org.Hs.eg.db", species = "hg19",
    anchor_gap = 200L, anchor_min_overlap = 1L,
    out_dir = tempdir(), write_output = FALSE, quiet = TRUE
  )
  has_loop <- !is.na(res_gap$target_annotation$Linked_Loop_IDs) &
    res_gap$target_annotation$Linked_Loop_IDs != ""
  expect_false(any(has_loop),
    info = "Peak within gap but without physical overlap should NOT link"
  )

  # Default (anchor_gap=-1L): strict overlap only, should NOT link
  res_strict <- annotate_peaks_and_loops(
    bedpe_file = loop_bedpe, target_bed = target_bed,
    txdb = txdb_obj, org_db = "org.Hs.eg.db", species = "hg19",
    out_dir = tempdir(), write_output = FALSE, quiet = TRUE
  )
  has_loop2 <- !is.na(res_strict$target_annotation$Linked_Loop_IDs) &
    res_strict$target_annotation$Linked_Loop_IDs != ""
  expect_false(any(has_loop2),
    info = "Peak without physical overlap should NOT be linked with default strict mode"
  )
  unlink(c(target_bed, loop_bedpe))
})

# --- validate_epeG_by_chromatin: evidence/missing-data logic ---
# Uses pre-computed annotation results with guaranteed P/G anchors.
test_that("validate_epeG_by_chromatin: no marks → all uncertain", {
  rdata_path <- system.file("extdata", "analysis_results.RData", package = "looplook")
  skip_if(rdata_path == "", "Pre-computed test data not available")
  tmp <- new.env()
  load(rdata_path, envir = tmp)
  res <- tmp[[ls(tmp)[1]]]

  val <- validate_epeG_by_chromatin(res, chromatin_beds = list(), quiet = TRUE)
  skip_if(nrow(val) == 0, "No P/G anchors in pre-computed data")
  expect_true(all(val$enhancer_evidence == "uncertain"))
  expect_true(all(is.na(val$H3K4me1) & is.na(val$H3K27ac) &
    is.na(val$ATAC) & is.na(val$H3K27me3) & is.na(val$H3K4me3)))
})

test_that("validate_epeG_by_chromatin: missing negative marks → not weak", {
  rdata_path <- system.file("extdata", "analysis_results.RData", package = "looplook")
  skip_if(rdata_path == "", "Pre-computed test data not available")
  tmp <- new.env()
  load(rdata_path, envir = tmp)
  res <- tmp[[ls(tmp)[1]]]

  first_anchor <- res$loop_annotation[1, ]
  h3k4me1 <- tempfile(fileext = ".bed")
  writeLines(sprintf(
    "%s\t%d\t%d", first_anchor$chr1, first_anchor$start1,
    first_anchor$end1
  ), h3k4me1)

  val <- validate_epeG_by_chromatin(res, chromatin_beds = list(
    H3K4me1 = h3k4me1
  ), quiet = TRUE)
  skip_if(nrow(val) == 0, "No P/G anchors in pre-computed data")
  expect_true(all(is.na(val$H3K27me3)))
  expect_true(all(is.na(val$H3K4me3)))
  expect_false(any(val$enhancer_evidence == "limited"),
    info = "Missing negative marks should not produce 'weak' classification"
  )
  unlink(h3k4me1)
})

# --- .assign_chromatin_confidence: all enhancer evidence levels ---
test_that(".assign_chromatin_confidence: canonical (all 5 marks aligned)", {
  anchors <- data.frame(
    anchor_id = "A1", chr = "chr1", start = 1L, end = 100L,
    anchor_type = "eP", anchor_gene = "GENE1", cluster_id = "C1",
    stringsAsFactors = FALSE
  )
  known_marks <- c("H3K4me1", "H3K27ac", "ATAC", "H3K27me3", "H3K4me3")
  mm <- as.data.frame(matrix(NA, nrow = 1, ncol = 5, dimnames = list(NULL, known_marks)))
  mm$H3K4me1 <- TRUE
  mm$H3K27ac <- TRUE
  mm$ATAC <- TRUE
  mm$H3K27me3 <- FALSE
  mm$H3K4me3 <- FALSE

  res <- looplook:::.assign_chromatin_confidence(anchors, mm, known_marks, known_marks)
  expect_equal(as.character(res$enhancer_evidence), "canonical")
  expect_false(any(is.na(res[, known_marks])))
})

test_that(".assign_chromatin_confidence: strong (H3K4me1 + H3K27ac)", {
  anchors <- data.frame(
    anchor_id = "A1", chr = "chr1", start = 1L, end = 100L,
    anchor_type = "eP", anchor_gene = "GENE1", cluster_id = "C1",
    stringsAsFactors = FALSE
  )
  known_marks <- c("H3K4me1", "H3K27ac", "ATAC", "H3K27me3", "H3K4me3")
  mm <- as.data.frame(matrix(NA, nrow = 1, ncol = 5, dimnames = list(NULL, known_marks)))
  mm$H3K4me1 <- TRUE
  mm$H3K27ac <- TRUE
  # Only 2 marks provided → all_five = FALSE → canonical impossible
  # H3K4me1+ and H3K27ac+ → strong

  res <- looplook:::.assign_chromatin_confidence(
    anchors, mm,
    c("H3K4me1", "H3K27ac"), known_marks
  )
  expect_equal(as.character(res$enhancer_evidence), "strong")
})

test_that(".assign_chromatin_confidence: supported (one positive mark)", {
  anchors <- data.frame(
    anchor_id = "A1", chr = "chr1", start = 1L, end = 100L,
    anchor_type = "eP", anchor_gene = "GENE1", cluster_id = "C1",
    stringsAsFactors = FALSE
  )
  known_marks <- c("H3K4me1", "H3K27ac", "ATAC", "H3K27me3", "H3K4me3")
  mm <- as.data.frame(matrix(NA, nrow = 1, ncol = 5, dimnames = list(NULL, known_marks)))
  mm$ATAC <- TRUE # only ATAC positive, no H3K4me1

  res <- looplook:::.assign_chromatin_confidence(anchors, mm, c("ATAC"), known_marks)
  expect_equal(as.character(res$enhancer_evidence), "supported")
})

test_that(".assign_chromatin_confidence: weak (negative marks tested, no positives)", {
  anchors <- data.frame(
    anchor_id = "A1", chr = "chr1", start = 1L, end = 100L,
    anchor_type = "eP", anchor_gene = "GENE1", cluster_id = "C1",
    stringsAsFactors = FALSE
  )
  known_marks <- c("H3K4me1", "H3K27ac", "ATAC", "H3K27me3", "H3K4me3")
  mm <- as.data.frame(matrix(NA, nrow = 1, ncol = 5, dimnames = list(NULL, known_marks)))
  mm$H3K27me3 <- FALSE
  mm$H3K4me3 <- FALSE # negative marks tested, absent

  res <- looplook:::.assign_chromatin_confidence(
    anchors, mm,
    c("H3K27me3", "H3K4me3"), known_marks
  )
  expect_equal(as.character(res$enhancer_evidence), "limited")
})

test_that(".assign_chromatin_confidence: uncertain (all marks absent, NA cols OK)", {
  anchors <- data.frame(
    anchor_id = "A1", chr = "chr1", start = 1L, end = 100L,
    anchor_type = "eP", anchor_gene = "GENE1", cluster_id = "C1",
    stringsAsFactors = FALSE
  )
  known_marks <- c("H3K4me1", "H3K27ac", "ATAC", "H3K27me3", "H3K4me3")
  mm <- as.data.frame(matrix(NA, nrow = 1, ncol = 5, dimnames = list(NULL, known_marks)))

  res <- looplook:::.assign_chromatin_confidence(anchors, mm, character(0), known_marks)
  expect_equal(as.character(res$enhancer_evidence), "uncertain")
  # All mark columns should be NA
  expect_true(all(is.na(res[, known_marks])))
})

test_that(".assign_chromatin_confidence: NA negative marks do not produce weak", {
  anchors <- data.frame(
    anchor_id = "A1", chr = "chr1", start = 1L, end = 100L,
    anchor_type = "eP", anchor_gene = "GENE1", cluster_id = "C1",
    stringsAsFactors = FALSE
  )
  known_marks <- c("H3K4me1", "H3K27ac", "ATAC", "H3K27me3", "H3K4me3")
  mm <- as.data.frame(matrix(NA, nrow = 1, ncol = 5, dimnames = list(NULL, known_marks)))
  mm$H3K4me1 <- FALSE
  mm$H3K27ac <- FALSE
  mm$ATAC <- FALSE # all active marks tested & absent
  # negative marks are NA → should NOT trigger "limited"

  res <- looplook:::.assign_chromatin_confidence(
    anchors, mm,
    c("H3K4me1", "H3K27ac", "ATAC"), known_marks
  )
  expect_equal(as.character(res$enhancer_evidence), "uncertain")
  expect_true(all(is.na(res$H3K27me3)))
  expect_true(all(is.na(res$H3K4me3)))
})

test_that(".assign_chromatin_confidence: canonical fails with missing mark", {
  anchors <- data.frame(
    anchor_id = "A1", chr = "chr1", start = 1L, end = 100L,
    anchor_type = "eP", anchor_gene = "GENE1", cluster_id = "C1",
    stringsAsFactors = FALSE
  )
  known_marks <- c("H3K4me1", "H3K27ac", "ATAC", "H3K27me3", "H3K4me3")
  mm <- as.data.frame(matrix(NA, nrow = 1, ncol = 5, dimnames = list(NULL, known_marks)))
  mm$H3K4me1 <- TRUE
  mm$H3K27ac <- TRUE
  mm$ATAC <- TRUE
  mm$H3K27me3 <- FALSE
  mm$H3K4me3 <- NA # H3K4me3 NOT tested → all_five fails

  res <- looplook:::.assign_chromatin_confidence(
    anchors, mm,
    c("H3K4me1", "H3K27ac", "ATAC", "H3K27me3"), known_marks
  )
  # Not canonical (only 4 marks provided)
  expect_equal(as.character(res$enhancer_evidence), "strong")
})

# --- .record_database_versions ---
test_that(".record_database_versions returns expected structure", {
  dbv <- looplook:::.record_database_versions("hg38")
  expect_type(dbv, "list")
  expect_true(all(c(
    "TxDb", "OrgDb", "BSgenome", "JASPAR", "clusterProfiler",
    "txdb_pkg", "orgdb_pkg"
  ) %in% names(dbv)))
})

test_that(".record_database_versions handles NULL species", {
  dbv <- looplook:::.record_database_versions(NULL)
  expect_type(dbv, "list")
  expect_true(all(is.na(dbv$TxDb) & is.na(dbv$OrgDb) & is.na(dbv$BSgenome)))
})

# --- refine_loop_anchors_by_chromatin: E anchors included + target NULL ---
test_that("refine_loop_anchors_by_chromatin: recompute_targets=FALSE preserves upstream targets", {
  rdata_path <- system.file("extdata", "analysis_results.RData", package = "looplook")
  skip_if(rdata_path == "", "Pre-computed test data not available")
  tmp <- new.env()
  load(rdata_path, envir = tmp)
  res <- tmp[[ls(tmp)[1]]]
  skip_if(is.null(res$loop_annotation), "No loop_annotation in test data")

  h3k4me1 <- tempfile(fileext = ".bed")
  h3k4me3 <- tempfile(fileext = ".bed")
  writeLines(sprintf(
    "%s\t%d\t%d", res$loop_annotation$chr1[1],
    res$loop_annotation$start1[1], res$loop_annotation$end1[1]
  ), h3k4me1)
  writeLines(sprintf(
    "%s\t%d\t%d", res$loop_annotation$chr1[1],
    res$loop_annotation$start1[1], res$loop_annotation$end1[1]
  ), h3k4me3)

  out <- refine_loop_anchors_by_chromatin(res,
    chromatin_beds = list(H3K4me1 = h3k4me1, H3K4me3 = h3k4me3),
    write_output = FALSE, quiet = TRUE, recompute_targets = FALSE
  )
  # recompute_targets = FALSE preserves upstream target_annotation/target_gene_links
  # (not NULL — they are inherited from the input annotation_res)
  expect_false(
    is.null(out$target_annotation),
    "recompute_targets=FALSE should preserve upstream target_annotation"
  )
  expect_false(
    is.null(out$target_gene_links),
    "recompute_targets=FALSE should preserve upstream target_gene_links"
  )
  expect_s3_class(out$loop_annotation, "data.frame")
  expect_true("chromatin_validation" %in% names(out))
  unlink(c(h3k4me1, h3k4me3))
})

# --- recompute_targets: chromatin-aware target gene links ---
test_that("recompute_targets=TRUE returns non-NULL target tables", {
  rdata_path <- system.file("extdata", "analysis_results.RData", package = "looplook")
  skip_if(rdata_path == "", "Pre-computed test data not available")
  tmp <- new.env()
  load(rdata_path, envir = tmp)
  res <- tmp[[ls(tmp)[1]]]
  skip_if(
    !is.list(attr(res, "looplook_anchor_state")),
    "No anchor_state attribute in test data"
  )

  h3k4me1 <- tempfile(fileext = ".bed")
  h3k4me3 <- tempfile(fileext = ".bed")
  h3k27ac <- tempfile(fileext = ".bed")
  first_anchor <- res$loop_annotation[1, ]
  writeLines(sprintf(
    "%s\t%d\t%d", first_anchor$chr1,
    first_anchor$start1, first_anchor$end1
  ), h3k4me1)
  writeLines(sprintf(
    "%s\t%d\t%d", first_anchor$chr1,
    first_anchor$start1, first_anchor$end1
  ), h3k4me3)
  writeLines(sprintf(
    "%s\t%d\t%d", first_anchor$chr1,
    first_anchor$start1, first_anchor$end1
  ), h3k27ac)

  cr <- refine_loop_anchors_by_chromatin(res,
    chromatin_beds = list(
      H3K4me1 = h3k4me1, H3K4me3 = h3k4me3,
      H3K27ac = h3k27ac
    ),
    recompute_targets = TRUE, write_output = FALSE, quiet = TRUE
  )
  expect_true(!is.null(cr$target_annotation) && !is.null(cr$target_gene_links))
  unlink(c(h3k4me1, h3k4me3, h3k27ac))
})

test_that("recompute_targets=FALSE preserves upstream target tables", {
  rdata_path <- system.file("extdata", "analysis_results.RData", package = "looplook")
  skip_if(rdata_path == "", "Pre-computed test data not available")
  tmp <- new.env()
  load(rdata_path, envir = tmp)
  res <- tmp[[ls(tmp)[1]]]

  h3k4me1 <- tempfile(fileext = ".bed")
  h3k4me3 <- tempfile(fileext = ".bed")
  first_anchor <- res$loop_annotation[1, ]
  writeLines(sprintf(
    "%s\t%d\t%d", first_anchor$chr1,
    first_anchor$start1, first_anchor$end1
  ), h3k4me1)
  writeLines(sprintf(
    "%s\t%d\t%d", first_anchor$chr1,
    first_anchor$start1, first_anchor$end1
  ), h3k4me3)

  cr <- refine_loop_anchors_by_chromatin(res,
    chromatin_beds = list(H3K4me1 = h3k4me1, H3K4me3 = h3k4me3),
    recompute_targets = FALSE, write_output = FALSE, quiet = TRUE
  )
  expect_false(
    is.null(cr$target_annotation),
    "recompute_targets=FALSE should preserve upstream target_annotation"
  )
  expect_false(
    is.null(cr$target_gene_links),
    "recompute_targets=FALSE should preserve upstream target_gene_links"
  )
  unlink(c(h3k4me1, h3k4me3))
})

# --- chromatin-aware target remapping synthetic tests ---
# Use .assign_chromatin_confidence + .chromatin_reclassify directly
# to test reclassification rules without needing real TxDb/OrgDb.

.make_toy_validation <- function(old_type, h3k4me1, h3k4me3, h3k27ac = FALSE,
                                 atac = FALSE, h3k27me3 = FALSE) {
  known_marks <- c("H3K4me1", "H3K27ac", "ATAC", "H3K27me3", "H3K4me3")
  anchors <- data.frame(
    anchor_id = "A1", chr = "chr1", start = 1L, end = 100L,
    anchor_type = old_type, anchor_gene = "GENE1", cluster_id = "C1",
    stringsAsFactors = FALSE
  )
  mm <- as.data.frame(matrix(NA, 1, 5, dimnames = list(NULL, known_marks)))
  mm$H3K4me1 <- h3k4me1
  mm$H3K4me3 <- h3k4me3
  mm$H3K27ac <- h3k27ac
  mm$ATAC <- atac
  mm$H3K27me3 <- h3k27me3
  provided <- known_marks[!is.na(as.vector(mm[1, ]))]
  list(
    validation = looplook:::.assign_chromatin_confidence(anchors, mm, provided, known_marks),
    anchors = anchors,
    provided_marks = provided
  )
}

test_that("P + H3K4me1+ H3K4me3+ H3K27ac+ H3K27me3- (no bw) → P", {
  tv <- .make_toy_validation("P", TRUE, TRUE, TRUE, TRUE, FALSE)
  rm <- looplook:::.chromatin_reclassify(tv$validation, provided_marks = tv$provided_marks)
  expect_equal(rm$new_type[1], "P")
})

test_that("P + H3K4me1+ H3K4me3- H3K27ac+ → E (conservative)", {
  tv <- .make_toy_validation("P", TRUE, FALSE, TRUE, FALSE, FALSE)
  rm <- looplook:::.chromatin_reclassify(tv$validation, provided_marks = tv$provided_marks)
  expect_equal(rm$new_type[1], "E")
  expect_true(rm$changed[1])
})

test_that("P + H3K4me1+ H3K4me3- without H3K27ac/ATAC → stays P", {
  tv <- .make_toy_validation("P", TRUE, FALSE, FALSE, FALSE, FALSE)
  rm <- looplook:::.chromatin_reclassify(tv$validation, provided_marks = tv$provided_marks)
  expect_equal(rm$new_type[1], "P")
  expect_false(rm$changed[1])
})

test_that("E + H3K4me3+ H3K4me1- → P (unannotated promoter)", {
  tv <- .make_toy_validation("E", FALSE, TRUE, FALSE, FALSE, FALSE)
  rm <- looplook:::.chromatin_reclassify(tv$validation, provided_marks = tv$provided_marks)
  expect_equal(rm$new_type[1], "P")
  expect_true(rm$changed[1])
})

test_that("E + H3K4me1+ H3K4me3+ (no bw) → keeps E (unresolved dual)", {
  tv <- .make_toy_validation("E", TRUE, TRUE, TRUE, TRUE, FALSE)
  rm <- looplook:::.chromatin_reclassify(tv$validation, provided_marks = tv$provided_marks)
  expect_equal(rm$new_type[1], "E") # no bigWig → unresolved → keep old type
})

test_that("G + H3K4me1+ H3K4me3- H3K27ac+ → E (conservative intronic enhancer)", {
  tv <- .make_toy_validation("G", TRUE, FALSE, TRUE, FALSE, FALSE)
  rm <- looplook:::.chromatin_reclassify(tv$validation, provided_marks = tv$provided_marks)
  expect_equal(rm$new_type[1], "E")
  expect_true(rm$changed[1])
})

test_that("G + H3K4me1+ only → stays G (not enough evidence)", {
  tv <- .make_toy_validation("G", TRUE, FALSE, FALSE, FALSE, FALSE)
  rm <- looplook:::.chromatin_reclassify(tv$validation, provided_marks = tv$provided_marks)
  expect_equal(rm$new_type[1], "G")
  expect_false(rm$changed[1])
})

test_that("eP + dual_like (no bw) → eP (unresolved keeps old type)", {
  tv <- .make_toy_validation("eP", TRUE, TRUE, TRUE, TRUE, FALSE)
  rm <- looplook:::.chromatin_reclassify(tv$validation, provided_marks = tv$provided_marks)
  expect_equal(rm$new_type[1], "eP")
})

test_that("eP + H3K4me3+ only (promoter_like) → P", {
  tv <- .make_toy_validation("eP", FALSE, TRUE, FALSE, FALSE, FALSE)
  rm <- looplook:::.chromatin_reclassify(tv$validation, provided_marks = tv$provided_marks)
  # H3K4me3+ without H3K4me1 (promoter_like, no dual_like) → P via rule 5
  expect_equal(rm$new_type[1], "P")
})

test_that("eP + dual_like + me1_dominant → dual (bigWig resolved)", {
  tv <- .make_toy_validation("eP", TRUE, TRUE, FALSE, FALSE, FALSE)
  bw <- stats::setNames(2, tv$validation$anchor_id[1]) # log2(4) > log2(3)
  rm <- looplook:::.chromatin_reclassify(tv$validation,
    bw_ratio = bw,
    provided_marks = tv$provided_marks
  )
  # Rule 2: eP + dual_like + me1_dominant → dual
  expect_equal(rm$new_type[1], "dual")
})

test_that("eP + dual_like + not_me1_dominant → P (K4me3 dominates)", {
  tv <- .make_toy_validation("eP", TRUE, TRUE, FALSE, FALSE, FALSE)
  bw <- stats::setNames(1, tv$validation$anchor_id[1]) # log2(2) < log2(3)
  rm <- looplook:::.chromatin_reclassify(tv$validation,
    bw_ratio = bw,
    provided_marks = tv$provided_marks
  )
  # Rule 4: eP + dual_like + not_me1_dominant → P
  expect_equal(rm$new_type[1], "P")
})

test_that("eP + H3K4me3+ H3K27me3+ (bivalent) → P, not eP", {
  tv <- .make_toy_validation("eP", FALSE, TRUE, FALSE, FALSE, TRUE)
  rm <- looplook:::.chromatin_reclassify(tv$validation, provided_marks = tv$provided_marks)
  # Rule 5: eP + H3K4me3+ without H3K4me1 → P, even when H3K27me3 present
  expect_equal(rm$new_type[1], "P")
})

test_that("enhancer BED + eP + H3K4me3 → dual (highest priority)", {
  tv <- .make_toy_validation("eP", FALSE, TRUE, FALSE, FALSE, FALSE)
  rm <- looplook:::.chromatin_reclassify(tv$validation,
    enhancer_anchors = tv$validation$anchor_id[1],
    provided_marks = tv$provided_marks
  )
  # Rule 1: enhancer BED + H3K4me3 → dual (overrides promoter_like)
  expect_equal(rm$new_type[1], "dual")
})

# ══════════════════════════════════════════════════════════════════════════
# Triple-positive eP/eG regression tests (H3K4me1+/H3K4me3+/H3K27me3+)
# ══════════════════════════════════════════════════════════════════════════

test_that("triple-positive eP retains eP (conflicting_marks)", {
  tv <- .make_toy_validation("eP", h3k4me1 = TRUE, h3k4me3 = TRUE, h3k27me3 = TRUE)
  rm <- looplook:::.chromatin_reclassify(tv$validation, provided_marks = tv$provided_marks)
  expect_equal(rm$chromatin_state[1], "conflicting_marks")
  expect_equal(rm$new_type[1], "eP")
  expect_false(rm$changed[1])
})

test_that("triple-positive eG retains eG (conflicting_marks)", {
  tv <- .make_toy_validation("eG", h3k4me1 = TRUE, h3k4me3 = TRUE, h3k27me3 = TRUE)
  rm <- looplook:::.chromatin_reclassify(tv$validation, provided_marks = tv$provided_marks)
  expect_equal(rm$chromatin_state[1], "conflicting_marks")
  expect_equal(rm$new_type[1], "eG")
  expect_false(rm$changed[1])
})

test_that("bivalent eP without H3K4me1 reverts to P", {
  tv <- .make_toy_validation("eP", h3k4me1 = FALSE, h3k4me3 = TRUE, h3k27me3 = TRUE)
  rm <- looplook:::.chromatin_reclassify(tv$validation, provided_marks = tv$provided_marks)
  expect_equal(rm$chromatin_state[1], "conflicting_marks")
  expect_equal(rm$new_type[1], "P")
  expect_true(rm$changed[1])
})

# --- annotate_peaks_and_loops output contract: a1_id / a2_id + anchor_state ---
test_that("annotate_peaks_and_loops returns loop_annotation with a1_id / a2_id", {
  skip_if_not_installed("org.Hs.eg.db")

  res <- get_base_annotation_target()

  expect_true("a1_id" %in% colnames(res$loop_annotation),
    info = "loop_annotation must carry a1_id for downstream chromatin refinement"
  )
  expect_true("a2_id" %in% colnames(res$loop_annotation),
    info = "loop_annotation must carry a2_id for downstream chromatin refinement"
  )

  as <- attr(res, "looplook_anchor_state")
  expect_true(!is.null(as), info = "anchor_state attribute must be present")
  expect_true(!is.null(as$target_hit_df))
  expect_true(nrow(as$target_hit_df) > 0)
  expect_true("anchor_id" %in% colnames(as$target_hit_df))
})

test_that("annotate_peaks_and_loops without target_bed has anchor_state but NULL hit_df", {
  skip_if_not_installed("org.Hs.eg.db")
  txdb_obj <- get_test_txdb()

  tiny_bedpe <- tempfile(fileext = ".bedpe")
  writeLines("chr6\t10412000\t10412600\tchr6\t10415000\t10415600", tiny_bedpe)

  res <- annotate_peaks_and_loops(
    bedpe_file = tiny_bedpe,
    txdb = txdb_obj, org_db = "org.Hs.eg.db", species = "hg19",
    out_dir = tempdir(), write_output = FALSE, quiet = TRUE
  )

  as <- attr(res, "looplook_anchor_state")
  expect_true(!is.null(as))
  expect_true(is.null(as$target_hit_df) || nrow(as$target_hit_df) == 0,
    info = "target_hit_df should be empty/NULL when no target_bed provided"
  )

  unlink(tiny_bedpe)
})

# --- chromatin_target_action in recomputed target_gene_links ---
test_that("recompute_targets=TRUE populates chromatin_target_action", {
  skip_if_not_installed("org.Hs.eg.db")

  h3k4me1 <- tempfile(fileext = ".bed")
  h3k4me3 <- tempfile(fileext = ".bed")
  writeLines("chr6\t10411900\t10412700", h3k4me1)
  writeLines("chr6\t10411900\t10412700", h3k4me3)

  base_res <- get_base_annotation_target()

  cr <- refine_loop_anchors_by_chromatin(base_res,
    chromatin_beds = list(H3K4me1 = h3k4me1, H3K4me3 = h3k4me3),
    recompute_targets = TRUE, write_output = FALSE, quiet = TRUE
  )

  expect_true(!is.null(cr$target_gene_links),
    info = "recompute_targets=TRUE must return target_gene_links"
  )
  if (nrow(cr$target_gene_links) > 0) {
    expect_true("chromatin_target_action" %in% colnames(cr$target_gene_links),
      info = "target_gene_links must carry chromatin_target_action after recompute"
    )
  }

  unlink(c(h3k4me1, h3k4me3))
})

# --- recompute_targets=FALSE leaves target tables untouched ---
test_that("recompute_targets=FALSE preserves upstream target_annotation", {
  skip_if_not_installed("org.Hs.eg.db")

  h3k4me1 <- tempfile(fileext = ".bed")
  h3k4me3 <- tempfile(fileext = ".bed")
  writeLines("chr6\t10411900\t10412700", h3k4me1)
  writeLines("chr6\t10411900\t10412700", h3k4me3)

  base_res <- get_base_annotation_target()

  cr <- refine_loop_anchors_by_chromatin(base_res,
    chromatin_beds = list(H3K4me1 = h3k4me1, H3K4me3 = h3k4me3),
    recompute_targets = FALSE, write_output = FALSE, quiet = TRUE
  )

  expect_false(is.null(cr$target_annotation),
    info = "recompute_targets=FALSE: target_annotation should be preserved, not NULL"
  )
  expect_false(is.null(cr$target_gene_links),
    info = "recompute_targets=FALSE: target_gene_links should be preserved, not NULL"
  )

  unlink(c(h3k4me1, h3k4me3))
})

# --- recompute_targets: biological consequence checks ---
test_that("chromatin reclassification propagates anchor type changes into loop_type", {
  skip_if_not_installed("org.Hs.eg.db")

  # Anchor A1 (chr6:10412000-10412600) overlaps both a promoter and the
  # H3K4me1/H3K27ac marks but NOT H3K4me3 -> P + H3K4me1+ H3K4me3- H3K27ac+
  # must be reclassified to E. A2 receives no marks and stays P.
  h3k4me1 <- tempfile(fileext = ".bed")
  h3k4me3 <- tempfile(fileext = ".bed")
  h3k27ac <- tempfile(fileext = ".bed")
  writeLines("chr6\t10411900\t10412700", h3k4me1)
  writeLines("chr6\t20000000\t20001000", h3k4me3)
  writeLines("chr6\t10411900\t10412700", h3k27ac)

  base_res <- get_base_annotation_target()

  cr <- refine_loop_anchors_by_chromatin(base_res,
    chromatin_beds = list(H3K4me1 = h3k4me1, H3K4me3 = h3k4me3, H3K27ac = h3k27ac),
    recompute_targets = TRUE, write_output = FALSE, quiet = TRUE
  )

  cv <- cr$chromatin_validation
  pos_type <- as.character(cv$positional_type)
  fin_type <- as.character(cv$final_type)
  changed <- pos_type != fin_type
  expect_true(any(changed)) # the crafted marks must actually reclassify an anchor
  changed_ids <- cv$anchor_id[changed]
  expect_true(length(changed_ids) > 0)
  expect_true(all(pos_type[changed] != fin_type[changed]))
  # The reclassified anchor type propagates into the loop topology.
  expect_true("E" %in% c(cr$loop_annotation$anchor1_type, cr$loop_annotation$anchor2_type))
  expect_true("E-P" %in% cr$loop_annotation$loop_type)

  unlink(c(h3k4me1, h3k4me3, h3k27ac))
})

test_that("recomputed target_gene_links carry chromatin provenance columns", {
  skip_if_not_installed("org.Hs.eg.db")

  h3k4me1 <- tempfile(fileext = ".bed")
  h3k4me3 <- tempfile(fileext = ".bed")
  writeLines("chr6\t10411900\t10412700", h3k4me1)
  writeLines("chr6\t10411900\t10412700", h3k4me3)

  base_res <- get_base_annotation_target()

  cr <- refine_loop_anchors_by_chromatin(base_res,
    chromatin_beds = list(H3K4me1 = h3k4me1, H3K4me3 = h3k4me3),
    recompute_targets = TRUE, write_output = FALSE, quiet = TRUE
  )

  links <- cr$target_gene_links
  expect_true(!is.null(links) && nrow(links) > 0)
  expect_true("anchor_type_before_chromatin" %in% colnames(links))
  expect_true("anchor_type_after_chromatin" %in% colnames(links))
  expect_true("chromatin_confidence" %in% colnames(links))
  expect_true("chromatin_evidence" %in% colnames(links))

  unlink(c(h3k4me1, h3k4me3))
})

test_that("recomputed *_Filled columns consistent with strict 3D columns", {
  skip_if_not_installed("org.Hs.eg.db")

  h3k4me1 <- tempfile(fileext = ".bed")
  h3k4me3 <- tempfile(fileext = ".bed")
  writeLines("chr6\t10411900\t10412700", h3k4me1)
  writeLines("chr6\t10411900\t10412700", h3k4me3)

  base_res <- get_base_annotation_target()

  cr <- refine_loop_anchors_by_chromatin(base_res,
    chromatin_beds = list(H3K4me1 = h3k4me1, H3K4me3 = h3k4me3),
    recompute_targets = TRUE, write_output = FALSE, quiet = TRUE
  )

  ta <- cr$target_annotation
  expect_true("Regulated_promoter_genes_Filled" %in% colnames(ta))
  expect_true("Assigned_Target_Genes_Filled" %in% colnames(ta))
  expect_true("Regulated_promoter_Evidence" %in% colnames(ta))
  expect_true("Regulated_promoter_Fallback_Evidence" %in% colnames(ta))

  reg_filled <- ta$Regulated_promoter_genes_Filled
  reg_strict <- ta$Regulated_promoter_genes
  reg_evidence <- ta$Regulated_promoter_Fallback_Evidence

  has_strict <- !is.na(reg_strict) & nzchar(reg_strict)
  has_filled <- !is.na(reg_filled) & nzchar(reg_filled)

  if (any(has_strict & has_filled)) {
    idx <- which(has_strict & has_filled)
    expect_equal(reg_filled[idx], reg_strict[idx])
  }
  if (any(has_strict)) {
    expect_true(all(reg_evidence[has_strict] == "none"))
  }

  unlink(c(h3k4me1, h3k4me3))
})

# --- P0-1: expression refinement preserves and updates anchor_state ---
test_that("expression refinement preserves and updates looplook_anchor_state", {
  skip_if_not_installed("org.Hs.eg.db")
  txdb_obj <- get_test_txdb()

  tiny_bedpe <- tempfile(fileext = ".bedpe")
  target_bed <- tempfile(fileext = ".bed")
  expr_file <- tempfile(fileext = ".tsv")
  writeLines("chr6\t10412000\t10412600\tchr6\t10415000\t10415600", tiny_bedpe)
  writeLines("chr6\t10412000\t10413000", target_bed)
  writeLines("Gene\tS1\nGENE6\t0", expr_file)

  raw_res <- annotate_peaks_and_loops(
    bedpe_file = tiny_bedpe, target_bed = target_bed,
    txdb = txdb_obj, org_db = "org.Hs.eg.db", species = "hg19",
    expr_matrix_file = expr_file, sample_columns = "S1",
    out_dir = tempdir(), write_output = FALSE, quiet = TRUE
  )

  expr_res <- refine_loop_anchors_by_expression(raw_res,
    expr_matrix_file = expr_file, sample_columns = "S1",
    threshold = 0.5, species = "hg19",
    out_dir = tempdir(), write_output = FALSE, quiet = TRUE
  )

  st <- attr(expr_res, "looplook_anchor_state", exact = TRUE)
  expect_true(!is.null(st))
  expect_true("target_hit_df" %in% names(st))
  expect_true("map_info" %in% names(st))

  unlink(c(tiny_bedpe, target_bed, expr_file))
})
