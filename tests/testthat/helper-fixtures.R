# tests/testthat/helper-fixtures.R
# Session-wide test fixtures shared across test files. The hg19 sample TxDb
# from GenomicFeatures and the tiny example annotation bases are expensive to
# load/compute, so they are built once per test session and reused by every
# test file that needs them (annotate, refine, chromatin tests).

.fixture_cache <- new.env(parent = emptyenv())

.fixture_txdb_path <- function() {
  system.file("extdata", "hg19_knownGene_sample.sqlite", package = "GenomicFeatures")
}

#' Get the shared hg19 sample TxDb (loaded once per session)
get_test_txdb <- function() {
  if (is.null(.fixture_cache$txdb)) {
    p <- .fixture_txdb_path()
    if (p == "") {
      return(NULL)
    }
    .fixture_cache$txdb <- AnnotationDbi::loadDb(p)
  }
  .fixture_cache$txdb
}

#' Shared annotate base: default parameters with a target BED
#' (chr6:10412000-10412600 / 10415000-10415600 loops, target chr6:10412000-10413000).
get_base_annotation_target <- function() {
  if (is.null(.fixture_cache$target)) {
    txdb_obj <- get_test_txdb()
    tiny_bedpe <- tempfile(fileext = ".bedpe")
    target_bed <- tempfile(fileext = ".bed")
    writeLines("chr6\t10412000\t10412600\tchr6\t10415000\t10415600", tiny_bedpe)
    writeLines("chr6\t10412000\t10413000", target_bed)
    .fixture_cache$target <- annotate_peaks_and_loops(
      bedpe_file = tiny_bedpe, target_bed = target_bed,
      txdb = txdb_obj, org_db = "org.Hs.eg.db", species = "hg19",
      out_dir = tempdir(), write_output = FALSE, quiet = TRUE
    )
  }
  .fixture_cache$target
}

#' Shared annotate base: neighbor_hop = 1 without a target BED.
get_base_annotation_hop1 <- function() {
  if (is.null(.fixture_cache$hop1)) {
    txdb_obj <- get_test_txdb()
    tiny_bedpe <- tempfile(fileext = ".bedpe")
    writeLines("chr6\t10412000\t10412600\tchr6\t10415000\t10415600", tiny_bedpe)
    .fixture_cache$hop1 <- annotate_peaks_and_loops(
      bedpe_file = tiny_bedpe,
      txdb = txdb_obj, org_db = "org.Hs.eg.db", species = "hg19",
      neighbor_hop = 1, out_dir = tempdir(), write_output = FALSE, quiet = TRUE
    )
  }
  .fixture_cache$hop1
}
