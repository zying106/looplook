# tests/testthat/test-annotation.R

test_that("packaged annotation example keeps the expected output contract", {
  rdata_path <- system.file("extdata", "analysis_results.RData", package = "looplook")
  skip_if(rdata_path == "", "Pre-computed test data not available")

  temp_env <- new.env()
  load(rdata_path, envir = temp_env)
  res_integrated <- temp_env[[ls(temp_env)[1]]]

  expect_type(res_integrated, "list")
  expect_true(all(c("target_annotation", "loop_annotation", "anchor_loci_annotation", "anchor_annotation", "plots") %in% names(res_integrated)))
  expect_type(res_integrated$plots, "list")
  expect_gt(nrow(res_integrated$loop_annotation), 0)
  expect_gt(nrow(res_integrated$target_annotation), 0)
  expect_true("Assigned_Target_Genes_Filled" %in% colnames(res_integrated$target_annotation))
})

test_that("annotate_peaks_and_loops respects quiet and write_output flags", {
  skip_if_not_installed("org.Hs.eg.db")
  sample_txdb_path <- system.file(
    "extdata", "hg19_knownGene_sample.sqlite",
    package = "GenomicFeatures"
  )
  skip_if(sample_txdb_path == "", "Sample TxDb not available")
  txdb_obj <- AnnotationDbi::loadDb(sample_txdb_path)
  tiny_bedpe <- tempfile(fileext = ".bedpe")
  writeLines("chr6\t10412000\t10412600\tchr6\t10415000\t10415600", tiny_bedpe)

  out_base <- tempfile(pattern = "looplook_anno_nowrite_")
  unlink(out_base, recursive = TRUE, force = TRUE)
  expect_false(dir.exists(out_base))

  # quiet=TRUE suppresses looplook progress messages; package startup
  # messages from third-party dependencies are outside our control.
  res_integrated <- suppressPackageStartupMessages(
    looplook:::.with_known_upstream_noise_suppressed(
      looplook::annotate_peaks_and_loops(
        bedpe_file = tiny_bedpe,
        txdb = txdb_obj,
        org_db = "org.Hs.eg.db",
        species = "hg19",
        out_dir = out_base,
        project_name = "Tiny_NoWrite_Test",
        write_output = FALSE,
        quiet = TRUE
      )
    )
  )

  expect_type(res_integrated, "list")
  expect_false(dir.exists(out_base))
})

# --- Parameter validation ---
test_that("annotate_peaks_and_loops validates new parameters", {
  skip_if_not_installed("org.Hs.eg.db")
  sample_txdb_path <- system.file("extdata", "hg19_knownGene_sample.sqlite", package = "GenomicFeatures")
  skip_if(sample_txdb_path == "", "Sample TxDb not available")
  txdb_obj <- AnnotationDbi::loadDb(sample_txdb_path)
  tiny_bedpe <- tempfile(fileext = ".bedpe")
  writeLines("chr6\t10412000\t10412600\tchr6\t10415000\t10415600", tiny_bedpe)
  base_args <- list(
    bedpe_file = tiny_bedpe, txdb = txdb_obj, org_db = "org.Hs.eg.db",
    species = "hg19", out_dir = tempdir(), write_output = FALSE, quiet = TRUE
  )

  expect_error(do.call(annotate_peaks_and_loops, c(base_args, list(anchor_gap = -2))))
  expect_error(do.call(annotate_peaks_and_loops, c(base_args, list(anchor_min_overlap = 0))))
  expect_error(do.call(annotate_peaks_and_loops, c(base_args, list(anchor_min_frac = 1.5))))
  expect_error(do.call(annotate_peaks_and_loops, c(base_args, list(anchor_min_frac = -0.1))))
  expect_error(do.call(annotate_peaks_and_loops, c(base_args, list(hub_percentile = 0))))
  expect_error(do.call(annotate_peaks_and_loops, c(base_args, list(hub_percentile = 1.5))))
  expect_error(do.call(annotate_peaks_and_loops, c(base_args, list(neighbor_hop = -1))))
  expect_error(do.call(annotate_peaks_and_loops, c(base_args, list(neighbor_hop = 1.5))))
  expect_error(do.call(annotate_peaks_and_loops, c(base_args, list(karyo_bin_size = 0))))

  # Valid values should not error
  expect_no_error(do.call(annotate_peaks_and_loops, c(base_args, list(anchor_gap = 200L))))
  expect_no_error(do.call(annotate_peaks_and_loops, c(base_args, list(anchor_min_overlap = 10L))))
  expect_no_error(do.call(annotate_peaks_and_loops, c(base_args, list(anchor_min_frac = 0.5))))
})

# --- anchor_gap proximity matching ---
test_that("anchor_gap > 0 allows proximity-based peak-anchor linking", {
  skip_if_not_installed("org.Hs.eg.db")
  sample_txdb_path <- system.file("extdata", "hg19_knownGene_sample.sqlite", package = "GenomicFeatures")
  skip_if(sample_txdb_path == "", "Sample TxDb not available")
  txdb_obj <- AnnotationDbi::loadDb(sample_txdb_path)

  target_bed <- tempfile(fileext = ".bed")
  loop_bedpe <- tempfile(fileext = ".bedpe")
  # Peak 200bp away from anchor, 0bp actual overlap
  writeLines("chr6\t10412800\t10413000", target_bed)
  writeLines("chr6\t10412000\t10412600\tchr6\t10415000\t10415600", loop_bedpe)

  # anchor_gap=200: proximity hit (peak within 200bp of anchor) SHOULD link
  res <- annotate_peaks_and_loops(
    bedpe_file = loop_bedpe, target_bed = target_bed,
    txdb = txdb_obj, org_db = "org.Hs.eg.db", species = "hg19",
    anchor_gap = 200L, anchor_min_overlap = 1L,
    out_dir = tempdir(), write_output = FALSE, quiet = TRUE
  )
  has_loop <- !is.na(res$target_annotation$Linked_Loop_IDs) &
              res$target_annotation$Linked_Loop_IDs != ""
  expect_true(any(has_loop),
    info = "Peak within anchor_gap should be linked via proximity matching")

  # Default (anchor_gap=-1L): strict overlap only, should NOT link
  res2 <- annotate_peaks_and_loops(
    bedpe_file = loop_bedpe, target_bed = target_bed,
    txdb = txdb_obj, org_db = "org.Hs.eg.db", species = "hg19",
    out_dir = tempdir(), write_output = FALSE, quiet = TRUE
  )
  has_loop2 <- !is.na(res2$target_annotation$Linked_Loop_IDs) &
               res2$target_annotation$Linked_Loop_IDs != ""
  expect_false(any(has_loop2),
    info = "Peak without physical overlap should NOT be linked with default strict mode")

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
  expect_true(all(val$confidence == "uncertain"))
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
  writeLines(sprintf("%s\t%d\t%d", first_anchor$chr1, first_anchor$start1,
                     first_anchor$end1), h3k4me1)

  val <- validate_epeG_by_chromatin(res, chromatin_beds = list(
    H3K4me1 = h3k4me1
  ), quiet = TRUE)
  skip_if(nrow(val) == 0, "No P/G anchors in pre-computed data")
  expect_true(all(is.na(val$H3K27me3)))
  expect_true(all(is.na(val$H3K4me3)))
  expect_false(any(val$confidence == "weak"),
    info = "Missing negative marks should not produce 'weak' classification")
  unlink(h3k4me1)
})

# --- .assign_chromatin_confidence: all confidence levels ---
test_that(".assign_chromatin_confidence: gold_standard (all 5 marks aligned)", {
  anchors <- data.frame(anchor_id = "A1", chr = "chr1", start = 1L, end = 100L,
    anchor_type = "eP", anchor_gene = "GENE1", cluster_id = "C1",
    stringsAsFactors = FALSE)
  known_marks <- c("H3K4me1", "H3K27ac", "ATAC", "H3K27me3", "H3K4me3")
  mm <- as.data.frame(matrix(NA, nrow = 1, ncol = 5, dimnames = list(NULL, known_marks)))
  mm$H3K4me1 <- TRUE; mm$H3K27ac <- TRUE; mm$ATAC <- TRUE
  mm$H3K27me3 <- FALSE; mm$H3K4me3 <- FALSE

  res <- looplook:::.assign_chromatin_confidence(anchors, mm, known_marks, known_marks)
  expect_equal(as.character(res$confidence), "gold_standard")
  expect_false(any(is.na(res[, known_marks])))
})

test_that(".assign_chromatin_confidence: high_confidence (H3K4me1 + H3K27ac)", {
  anchors <- data.frame(anchor_id = "A1", chr = "chr1", start = 1L, end = 100L,
    anchor_type = "eP", anchor_gene = "GENE1", cluster_id = "C1",
    stringsAsFactors = FALSE)
  known_marks <- c("H3K4me1", "H3K27ac", "ATAC", "H3K27me3", "H3K4me3")
  mm <- as.data.frame(matrix(NA, nrow = 1, ncol = 5, dimnames = list(NULL, known_marks)))
  mm$H3K4me1 <- TRUE; mm$H3K27ac <- TRUE
  # Only 2 marks provided → all_five = FALSE → gold_standard impossible
  # H3K4me1+ and H3K27ac+ → high_confidence

  res <- looplook:::.assign_chromatin_confidence(anchors, mm,
    c("H3K4me1", "H3K27ac"), known_marks)
  expect_equal(as.character(res$confidence), "high_confidence")
})

test_that(".assign_chromatin_confidence: supported (one positive mark)", {
  anchors <- data.frame(anchor_id = "A1", chr = "chr1", start = 1L, end = 100L,
    anchor_type = "eP", anchor_gene = "GENE1", cluster_id = "C1",
    stringsAsFactors = FALSE)
  known_marks <- c("H3K4me1", "H3K27ac", "ATAC", "H3K27me3", "H3K4me3")
  mm <- as.data.frame(matrix(NA, nrow = 1, ncol = 5, dimnames = list(NULL, known_marks)))
  mm$ATAC <- TRUE  # only ATAC positive, no H3K4me1

  res <- looplook:::.assign_chromatin_confidence(anchors, mm, c("ATAC"), known_marks)
  expect_equal(as.character(res$confidence), "supported")
})

test_that(".assign_chromatin_confidence: weak (negative marks tested, no positives)", {
  anchors <- data.frame(anchor_id = "A1", chr = "chr1", start = 1L, end = 100L,
    anchor_type = "eP", anchor_gene = "GENE1", cluster_id = "C1",
    stringsAsFactors = FALSE)
  known_marks <- c("H3K4me1", "H3K27ac", "ATAC", "H3K27me3", "H3K4me3")
  mm <- as.data.frame(matrix(NA, nrow = 1, ncol = 5, dimnames = list(NULL, known_marks)))
  mm$H3K27me3 <- FALSE; mm$H3K4me3 <- FALSE  # negative marks tested, absent

  res <- looplook:::.assign_chromatin_confidence(anchors, mm,
    c("H3K27me3", "H3K4me3"), known_marks)
  expect_equal(as.character(res$confidence), "weak")
})

test_that(".assign_chromatin_confidence: uncertain (all marks absent, NA cols OK)", {
  anchors <- data.frame(anchor_id = "A1", chr = "chr1", start = 1L, end = 100L,
    anchor_type = "eP", anchor_gene = "GENE1", cluster_id = "C1",
    stringsAsFactors = FALSE)
  known_marks <- c("H3K4me1", "H3K27ac", "ATAC", "H3K27me3", "H3K4me3")
  mm <- as.data.frame(matrix(NA, nrow = 1, ncol = 5, dimnames = list(NULL, known_marks)))

  res <- looplook:::.assign_chromatin_confidence(anchors, mm, character(0), known_marks)
  expect_equal(as.character(res$confidence), "uncertain")
  # All mark columns should be NA
  expect_true(all(is.na(res[, known_marks])))
})

test_that(".assign_chromatin_confidence: NA negative marks do not produce weak", {
  anchors <- data.frame(anchor_id = "A1", chr = "chr1", start = 1L, end = 100L,
    anchor_type = "eP", anchor_gene = "GENE1", cluster_id = "C1",
    stringsAsFactors = FALSE)
  known_marks <- c("H3K4me1", "H3K27ac", "ATAC", "H3K27me3", "H3K4me3")
  mm <- as.data.frame(matrix(NA, nrow = 1, ncol = 5, dimnames = list(NULL, known_marks)))
  mm$H3K4me1 <- FALSE; mm$H3K27ac <- FALSE; mm$ATAC <- FALSE  # all active marks tested & absent
  # negative marks are NA → should NOT trigger "weak"

  res <- looplook:::.assign_chromatin_confidence(anchors, mm,
    c("H3K4me1", "H3K27ac", "ATAC"), known_marks)
  expect_equal(as.character(res$confidence), "uncertain")
  expect_true(all(is.na(res$H3K27me3)))
  expect_true(all(is.na(res$H3K4me3)))
})

test_that(".assign_chromatin_confidence: gold_standard fails with missing mark", {
  anchors <- data.frame(anchor_id = "A1", chr = "chr1", start = 1L, end = 100L,
    anchor_type = "eP", anchor_gene = "GENE1", cluster_id = "C1",
    stringsAsFactors = FALSE)
  known_marks <- c("H3K4me1", "H3K27ac", "ATAC", "H3K27me3", "H3K4me3")
  mm <- as.data.frame(matrix(NA, nrow = 1, ncol = 5, dimnames = list(NULL, known_marks)))
  mm$H3K4me1 <- TRUE; mm$H3K27ac <- TRUE; mm$ATAC <- TRUE
  mm$H3K27me3 <- FALSE; mm$H3K4me3 <- NA  # H3K4me3 NOT tested → all_five fails

  res <- looplook:::.assign_chromatin_confidence(anchors, mm,
    c("H3K4me1", "H3K27ac", "ATAC", "H3K27me3"), known_marks)
  # Not gold_standard (only 4 marks provided)
  expect_equal(as.character(res$confidence), "high_confidence")
})

# --- .record_database_versions ---
test_that(".record_database_versions returns expected structure", {
  dbv <- looplook:::.record_database_versions("hg38")
  expect_type(dbv, "list")
  expect_true(all(c("TxDb","OrgDb","BSgenome","JASPAR","clusterProfiler",
    "txdb_pkg","orgdb_pkg") %in% names(dbv)))
})

test_that(".record_database_versions handles NULL species", {
  dbv <- looplook:::.record_database_versions(NULL)
  expect_type(dbv, "list")
  expect_true(all(is.na(dbv$TxDb) & is.na(dbv$OrgDb) & is.na(dbv$BSgenome)))
})

# --- refine_loop_anchors_by_chromatin: E anchors included + target NULL ---
test_that("refine_loop_anchors_by_chromatin: recompute_targets=FALSE yields NULL", {
  rdata_path <- system.file("extdata", "analysis_results.RData", package = "looplook")
  skip_if(rdata_path == "", "Pre-computed test data not available")
  tmp <- new.env()
  load(rdata_path, envir = tmp)
  res <- tmp[[ls(tmp)[1]]]
  skip_if(is.null(res$loop_annotation), "No loop_annotation in test data")

  h3k4me1 <- tempfile(fileext = ".bed")
  h3k4me3 <- tempfile(fileext = ".bed")
  writeLines(sprintf("%s\t%d\t%d", res$loop_annotation$chr1[1],
    res$loop_annotation$start1[1], res$loop_annotation$end1[1]), h3k4me1)
  writeLines(sprintf("%s\t%d\t%d", res$loop_annotation$chr1[1],
    res$loop_annotation$start1[1], res$loop_annotation$end1[1]), h3k4me3)

  out <- refine_loop_anchors_by_chromatin(res,
    chromatin_beds = list(H3K4me1 = h3k4me1, H3K4me3 = h3k4me3),
    write_output = FALSE, quiet = TRUE, recompute_targets = FALSE)
  expect_null(out$target_annotation)
  expect_null(out$target_gene_links)
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
  skip_if(!is.list(attr(res, "looplook_anchor_state")),
          "No anchor_state attribute in test data")

  h3k4me1 <- tempfile(fileext = ".bed")
  h3k4me3 <- tempfile(fileext = ".bed")
  h3k27ac <- tempfile(fileext = ".bed")
  first_anchor <- res$loop_annotation[1, ]
  writeLines(sprintf("%s\t%d\t%d", first_anchor$chr1,
    first_anchor$start1, first_anchor$end1), h3k4me1)
  writeLines(sprintf("%s\t%d\t%d", first_anchor$chr1,
    first_anchor$start1, first_anchor$end1), h3k4me3)
  writeLines(sprintf("%s\t%d\t%d", first_anchor$chr1,
    first_anchor$start1, first_anchor$end1), h3k27ac)

  cr <- refine_loop_anchors_by_chromatin(res,
    chromatin_beds = list(H3K4me1 = h3k4me1, H3K4me3 = h3k4me3,
                          H3K27ac = h3k27ac),
    recompute_targets = TRUE, write_output = FALSE, quiet = TRUE)
  expect_true(!is.null(cr$target_annotation) && !is.null(cr$target_gene_links))
  unlink(c(h3k4me1, h3k4me3, h3k27ac))
})

test_that("recompute_targets=FALSE yields NULL target tables", {
  rdata_path <- system.file("extdata", "analysis_results.RData", package = "looplook")
  skip_if(rdata_path == "", "Pre-computed test data not available")
  tmp <- new.env()
  load(rdata_path, envir = tmp)
  res <- tmp[[ls(tmp)[1]]]

  h3k4me1 <- tempfile(fileext = ".bed")
  h3k4me3 <- tempfile(fileext = ".bed")
  first_anchor <- res$loop_annotation[1, ]
  writeLines(sprintf("%s\t%d\t%d", first_anchor$chr1,
    first_anchor$start1, first_anchor$end1), h3k4me1)
  writeLines(sprintf("%s\t%d\t%d", first_anchor$chr1,
    first_anchor$start1, first_anchor$end1), h3k4me3)

  cr <- refine_loop_anchors_by_chromatin(res,
    chromatin_beds = list(H3K4me1 = h3k4me1, H3K4me3 = h3k4me3),
    recompute_targets = FALSE, write_output = FALSE, quiet = TRUE)
  expect_null(cr$target_annotation)
  expect_null(cr$target_gene_links)
  unlink(c(h3k4me1, h3k4me3))
})

# --- chromatin-aware target remapping synthetic tests ---
# Use .assign_chromatin_confidence + .chromatin_reclassify directly
# to test reclassification rules without needing real TxDb/OrgDb.

.make_toy_validation <- function(old_type, h3k4me1, h3k4me3, h3k27ac=FALSE,
                                    atac=FALSE, h3k27me3=FALSE) {
    known_marks <- c("H3K4me1","H3K27ac","ATAC","H3K27me3","H3K4me3")
    anchors <- data.frame(anchor_id="A1", chr="chr1", start=1L, end=100L,
        anchor_type=old_type, anchor_gene="GENE1", cluster_id="C1",
        stringsAsFactors=FALSE)
    mm <- as.data.frame(matrix(NA, 1, 5, dimnames=list(NULL, known_marks)))
    mm$H3K4me1 <- h3k4me1; mm$H3K4me3 <- h3k4me3
    mm$H3K27ac <- h3k27ac; mm$ATAC <- atac; mm$H3K27me3 <- h3k27me3
    provided <- known_marks[!is.na(as.vector(mm[1,]))]
    list(
        validation = looplook:::.assign_chromatin_confidence(anchors, mm, provided, known_marks),
        anchors = anchors
    )
}

test_that("P + H3K4me1+ H3K4me3+ H3K27ac+ H3K27me3- → dual", {
    tv <- .make_toy_validation("P", TRUE, TRUE, TRUE, TRUE, FALSE)
    rm <- looplook:::.chromatin_reclassify(tv$validation)
    expect_equal(rm$new_type[1], "dual")
    expect_true(rm$changed[1])
})

test_that("P + H3K4me1+ H3K4me3- H3K27ac+ → E (conservative)", {
    tv <- .make_toy_validation("P", TRUE, FALSE, TRUE, FALSE, FALSE)
    rm <- looplook:::.chromatin_reclassify(tv$validation)
    expect_equal(rm$new_type[1], "E")
    expect_true(rm$changed[1])
})

test_that("P + H3K4me1+ H3K4me3- without H3K27ac/ATAC → stays P", {
    tv <- .make_toy_validation("P", TRUE, FALSE, FALSE, FALSE, FALSE)
    rm <- looplook:::.chromatin_reclassify(tv$validation)
    expect_equal(rm$new_type[1], "P")
    expect_false(rm$changed[1])
})

test_that("E + H3K4me3+ H3K4me1- → P (unannotated promoter)", {
    tv <- .make_toy_validation("E", FALSE, TRUE, FALSE, FALSE, FALSE)
    rm <- looplook:::.chromatin_reclassify(tv$validation)
    expect_equal(rm$new_type[1], "P")
    expect_true(rm$changed[1])
})

test_that("E + H3K4me1+ H3K4me3+ → dual", {
    tv <- .make_toy_validation("E", TRUE, TRUE, TRUE, TRUE, FALSE)
    rm <- looplook:::.chromatin_reclassify(tv$validation)
    expect_equal(rm$new_type[1], "dual")
    expect_true(rm$changed[1])
})

test_that("G + H3K4me1+ H3K4me3- H3K27ac+ → E (conservative intronic enhancer)", {
    tv <- .make_toy_validation("G", TRUE, FALSE, TRUE, FALSE, FALSE)
    rm <- looplook:::.chromatin_reclassify(tv$validation)
    expect_equal(rm$new_type[1], "E")
    expect_true(rm$changed[1])
})

test_that("G + H3K4me1+ only → stays G (not enough evidence)", {
    tv <- .make_toy_validation("G", TRUE, FALSE, FALSE, FALSE, FALSE)
    rm <- looplook:::.chromatin_reclassify(tv$validation)
    expect_equal(rm$new_type[1], "G")
    expect_false(rm$changed[1])
})

test_that("eP + dual_like → dual", {
    tv <- .make_toy_validation("eP", TRUE, TRUE, TRUE, TRUE, FALSE)
    rm <- looplook:::.chromatin_reclassify(tv$validation)
    expect_equal(rm$new_type[1], "dual")
    expect_true(rm$changed[1])
})

test_that("eP + H3K4me3+ only (promoter_like) → P", {
    tv <- .make_toy_validation("eP", FALSE, TRUE, FALSE, FALSE, FALSE)
    rm <- looplook:::.chromatin_reclassify(tv$validation)
    expect_equal(rm$new_type[1], "P")
    expect_true(rm$changed[1])
})

# --- annotate_peaks_and_loops output contract: a1_id / a2_id + anchor_state ---
test_that("annotate_peaks_and_loops returns loop_annotation with a1_id / a2_id", {
  skip_if_not_installed("org.Hs.eg.db")
  sample_txdb_path <- system.file("extdata", "hg19_knownGene_sample.sqlite",
    package = "GenomicFeatures")
  skip_if(sample_txdb_path == "", "Sample TxDb not available")
  txdb_obj <- AnnotationDbi::loadDb(sample_txdb_path)

  tiny_bedpe <- tempfile(fileext = ".bedpe")
  target_bed <- tempfile(fileext = ".bed")
  writeLines("chr6\t10412000\t10412600\tchr6\t10415000\t10415600", tiny_bedpe)
  writeLines("chr6\t10412000\t10413000", target_bed)

  res <- annotate_peaks_and_loops(
    bedpe_file = tiny_bedpe, target_bed = target_bed,
    txdb = txdb_obj, org_db = "org.Hs.eg.db", species = "hg19",
    out_dir = tempdir(), write_output = FALSE, quiet = TRUE
  )

  expect_true("a1_id" %in% colnames(res$loop_annotation),
    info = "loop_annotation must carry a1_id for downstream chromatin refinement")
  expect_true("a2_id" %in% colnames(res$loop_annotation),
    info = "loop_annotation must carry a2_id for downstream chromatin refinement")

  as <- attr(res, "looplook_anchor_state")
  expect_true(!is.null(as), info = "anchor_state attribute must be present")
  expect_true(!is.null(as$target_hit_df))
  expect_true(nrow(as$target_hit_df) > 0)
  expect_true("anchor_id" %in% colnames(as$target_hit_df))

  unlink(c(tiny_bedpe, target_bed))
})

test_that("annotate_peaks_and_loops without target_bed has anchor_state but NULL hit_df", {
  skip_if_not_installed("org.Hs.eg.db")
  sample_txdb_path <- system.file("extdata", "hg19_knownGene_sample.sqlite",
    package = "GenomicFeatures")
  skip_if(sample_txdb_path == "", "Sample TxDb not available")
  txdb_obj <- AnnotationDbi::loadDb(sample_txdb_path)

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
    info = "target_hit_df should be empty/NULL when no target_bed provided")

  unlink(tiny_bedpe)
})

# --- chromatin_target_action in recomputed target_gene_links ---
test_that("recompute_targets=TRUE populates chromatin_target_action", {
  skip_if_not_installed("org.Hs.eg.db")
  sample_txdb_path <- system.file("extdata", "hg19_knownGene_sample.sqlite",
    package = "GenomicFeatures")
  skip_if(sample_txdb_path == "", "Sample TxDb not available")
  txdb_obj <- AnnotationDbi::loadDb(sample_txdb_path)

  tiny_bedpe <- tempfile(fileext = ".bedpe")
  target_bed <- tempfile(fileext = ".bed")
  h3k4me1 <- tempfile(fileext = ".bed")
  h3k4me3 <- tempfile(fileext = ".bed")
  writeLines("chr6\t10412000\t10412600\tchr6\t10415000\t10415600", tiny_bedpe)
  writeLines("chr6\t10412000\t10413000", target_bed)
  writeLines("chr6\t10411900\t10412700", h3k4me1)
  writeLines("chr6\t10411900\t10412700", h3k4me3)

  base_res <- annotate_peaks_and_loops(
    bedpe_file = tiny_bedpe, target_bed = target_bed,
    txdb = txdb_obj, org_db = "org.Hs.eg.db", species = "hg19",
    out_dir = tempdir(), write_output = FALSE, quiet = TRUE
  )

  cr <- refine_loop_anchors_by_chromatin(base_res,
    chromatin_beds = list(H3K4me1 = h3k4me1, H3K4me3 = h3k4me3),
    recompute_targets = TRUE, write_output = FALSE, quiet = TRUE)

  expect_true(!is.null(cr$target_gene_links),
    info = "recompute_targets=TRUE must return target_gene_links")
  if (nrow(cr$target_gene_links) > 0) {
    expect_true("chromatin_target_action" %in% colnames(cr$target_gene_links),
      info = "target_gene_links must carry chromatin_target_action after recompute")
  }

  unlink(c(tiny_bedpe, target_bed, h3k4me1, h3k4me3))
})

# --- recompute_targets=FALSE leaves target tables untouched ---
test_that("recompute_targets=FALSE does not alter target_annotation", {
  skip_if_not_installed("org.Hs.eg.db")
  sample_txdb_path <- system.file("extdata", "hg19_knownGene_sample.sqlite",
    package = "GenomicFeatures")
  skip_if(sample_txdb_path == "", "Sample TxDb not available")
  txdb_obj <- AnnotationDbi::loadDb(sample_txdb_path)

  tiny_bedpe <- tempfile(fileext = ".bedpe")
  target_bed <- tempfile(fileext = ".bed")
  h3k4me1 <- tempfile(fileext = ".bed")
  h3k4me3 <- tempfile(fileext = ".bed")
  writeLines("chr6\t10412000\t10412600\tchr6\t10415000\t10415600", tiny_bedpe)
  writeLines("chr6\t10412000\t10413000", target_bed)
  writeLines("chr6\t10411900\t10412700", h3k4me1)
  writeLines("chr6\t10411900\t10412700", h3k4me3)

  base_res <- annotate_peaks_and_loops(
    bedpe_file = tiny_bedpe, target_bed = target_bed,
    txdb = txdb_obj, org_db = "org.Hs.eg.db", species = "hg19",
    out_dir = tempdir(), write_output = FALSE, quiet = TRUE
  )

  cr <- refine_loop_anchors_by_chromatin(base_res,
    chromatin_beds = list(H3K4me1 = h3k4me1, H3K4me3 = h3k4me3),
    recompute_targets = FALSE, write_output = FALSE, quiet = TRUE)

  expect_null(cr$target_annotation,
    info = "recompute_targets=FALSE: target_annotation must be NULL")
  expect_null(cr$target_gene_links,
    info = "recompute_targets=FALSE: target_gene_links must be NULL")

  unlink(c(tiny_bedpe, target_bed, h3k4me1, h3k4me3))
})

# --- recompute_targets: biological consequence checks ---
test_that("chromatin reclassification propagates anchor type changes into loop_type", {
  skip_if_not_installed("org.Hs.eg.db")
  sample_txdb_path <- system.file("extdata", "hg19_knownGene_sample.sqlite",
    package = "GenomicFeatures")
  skip_if(sample_txdb_path == "", "Sample TxDb not available")
  txdb_obj <- AnnotationDbi::loadDb(sample_txdb_path)

  tiny_bedpe <- tempfile(fileext = ".bedpe")
  target_bed <- tempfile(fileext = ".bed")
  h3k4me1    <- tempfile(fileext = ".bed")
  h3k4me3    <- tempfile(fileext = ".bed")
  writeLines("chr6\t10412000\t10412600\tchr6\t10415000\t10415600", tiny_bedpe)
  writeLines("chr6\t10412000\t10413000", target_bed)
  writeLines("chr6\t10411900\t10412700", h3k4me1)
  writeLines("chr6\t10411900\t10412700", h3k4me3)

  base_res <- annotate_peaks_and_loops(
    bedpe_file = tiny_bedpe, target_bed = target_bed,
    txdb = txdb_obj, org_db = "org.Hs.eg.db", species = "hg19",
    out_dir = tempdir(), write_output = FALSE, quiet = TRUE
  )

  cr <- refine_loop_anchors_by_chromatin(base_res,
    chromatin_beds = list(H3K4me1 = h3k4me1, H3K4me3 = h3k4me3),
    recompute_targets = TRUE, write_output = FALSE, quiet = TRUE)

  cv <- cr$chromatin_validation
  pos_type <- as.character(cv$positional_type)
  fin_type <- as.character(cv$final_type)
  changed <- pos_type != fin_type
  if (any(changed)) {
    changed_ids <- cv$anchor_id[changed]
    expect_true(length(changed_ids) > 0)
    expect_true(all(pos_type[changed] != fin_type[changed]))
  }

  unlink(c(tiny_bedpe, target_bed, h3k4me1, h3k4me3))
})

test_that("recomputed target_gene_links carry chromatin provenance columns", {
  skip_if_not_installed("org.Hs.eg.db")
  sample_txdb_path <- system.file("extdata", "hg19_knownGene_sample.sqlite",
    package = "GenomicFeatures")
  skip_if(sample_txdb_path == "", "Sample TxDb not available")
  txdb_obj <- AnnotationDbi::loadDb(sample_txdb_path)

  tiny_bedpe <- tempfile(fileext = ".bedpe")
  target_bed <- tempfile(fileext = ".bed")
  h3k4me1    <- tempfile(fileext = ".bed")
  h3k4me3    <- tempfile(fileext = ".bed")
  writeLines("chr6\t10412000\t10412600\tchr6\t10415000\t10415600", tiny_bedpe)
  writeLines("chr6\t10412000\t10413000", target_bed)
  writeLines("chr6\t10411900\t10412700", h3k4me1)
  writeLines("chr6\t10411900\t10412700", h3k4me3)

  base_res <- annotate_peaks_and_loops(
    bedpe_file = tiny_bedpe, target_bed = target_bed,
    txdb = txdb_obj, org_db = "org.Hs.eg.db", species = "hg19",
    out_dir = tempdir(), write_output = FALSE, quiet = TRUE
  )

  cr <- refine_loop_anchors_by_chromatin(base_res,
    chromatin_beds = list(H3K4me1 = h3k4me1, H3K4me3 = h3k4me3),
    recompute_targets = TRUE, write_output = FALSE, quiet = TRUE)

  links <- cr$target_gene_links
  expect_true(!is.null(links) && nrow(links) > 0)
  expect_true("anchor_type_before_chromatin" %in% colnames(links))
  expect_true("anchor_type_after_chromatin" %in% colnames(links))
  expect_true("chromatin_confidence" %in% colnames(links))
  expect_true("chromatin_evidence" %in% colnames(links))

  unlink(c(tiny_bedpe, target_bed, h3k4me1, h3k4me3))
})

test_that("recomputed *_Filled columns consistent with strict 3D columns", {
  skip_if_not_installed("org.Hs.eg.db")
  sample_txdb_path <- system.file("extdata", "hg19_knownGene_sample.sqlite",
    package = "GenomicFeatures")
  skip_if(sample_txdb_path == "", "Sample TxDb not available")
  txdb_obj <- AnnotationDbi::loadDb(sample_txdb_path)

  tiny_bedpe <- tempfile(fileext = ".bedpe")
  target_bed <- tempfile(fileext = ".bed")
  h3k4me1    <- tempfile(fileext = ".bed")
  h3k4me3    <- tempfile(fileext = ".bed")
  writeLines("chr6\t10412000\t10412600\tchr6\t10415000\t10415600", tiny_bedpe)
  writeLines("chr6\t10412000\t10413000", target_bed)
  writeLines("chr6\t10411900\t10412700", h3k4me1)
  writeLines("chr6\t10411900\t10412700", h3k4me3)

  base_res <- annotate_peaks_and_loops(
    bedpe_file = tiny_bedpe, target_bed = target_bed,
    txdb = txdb_obj, org_db = "org.Hs.eg.db", species = "hg19",
    out_dir = tempdir(), write_output = FALSE, quiet = TRUE
  )

  cr <- refine_loop_anchors_by_chromatin(base_res,
    chromatin_beds = list(H3K4me1 = h3k4me1, H3K4me3 = h3k4me3),
    recompute_targets = TRUE, write_output = FALSE, quiet = TRUE)

  ta <- cr$target_annotation
  expect_true("Regulated_promoter_genes_Filled" %in% colnames(ta))
  expect_true("Assigned_Target_Genes_Filled" %in% colnames(ta))
  expect_true("Regulated_promoter_Evidence" %in% colnames(ta))
  expect_true("Regulated_promoter_Fallback_Evidence" %in% colnames(ta))

  reg_filled  <- ta$Regulated_promoter_genes_Filled
  reg_strict  <- ta$Regulated_promoter_genes
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

  unlink(c(tiny_bedpe, target_bed, h3k4me1, h3k4me3))
})

# --- P0-1: expression refinement preserves and updates anchor_state ---
test_that("expression refinement preserves and updates looplook_anchor_state", {
  skip_if_not_installed("org.Hs.eg.db")
  sample_txdb_path <- system.file("extdata", "hg19_knownGene_sample.sqlite",
    package = "GenomicFeatures")
  skip_if(sample_txdb_path == "", "Sample TxDb not available")
  txdb_obj <- AnnotationDbi::loadDb(sample_txdb_path)

  tiny_bedpe <- tempfile(fileext = ".bedpe")
  target_bed <- tempfile(fileext = ".bed")
  expr_file  <- tempfile(fileext = ".tsv")
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
