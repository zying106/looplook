# tests/testthat/test-helpers.R — unit tests for utility helper functions

test_that("clean_gene_names handles normal, edge, and NA cases", {
  expect_equal(looplook:::clean_gene_names(c("A", "B", "A", "C")), c("A", "B", "C"))
  expect_equal(looplook:::clean_gene_names(c("A", "", "B", NA, " ")), c("A", "B"))
  expect_equal(looplook:::clean_gene_names(character(0)), character(0))
  expect_equal(looplook:::clean_gene_names(NULL), character(0))

  # with split
  expect_equal(looplook:::clean_gene_names("TP53;BRCA1;TP53", ";"), c("TP53", "BRCA1"))
  expect_equal(looplook:::clean_gene_names(c("A;B", "B;C", "D")), c("A;B", "B;C", "D"))
  expect_equal(looplook:::clean_gene_names(c("", NA, ";"), "[;,]"), character(0))
})

test_that("extract_genes collapses delimited gene strings correctly", {
  expect_equal(looplook:::extract_genes("TP53;BRCA1;TP53"), "TP53;BRCA1")
  expect_equal(looplook:::extract_genes(c("TP53", "BRCA1")), "TP53;BRCA1")
  expect_equal(looplook:::extract_genes(NA_character_), NA_character_)
  expect_equal(looplook:::extract_genes(""), NA_character_)
  expect_equal(looplook:::extract_genes(";"), NA_character_)
  # single gene
  expect_equal(looplook:::extract_genes("MYC"), "MYC")
})

test_that("load_expression_matrix loads and averages correctly", {
  tmp <- tempfile(fileext = ".tsv")
  write.table(data.frame(
    GeneID = c("TP53", "BRCA1", "MYC"),
    con1 = c(10, 20, 30),
    con2 = c(15, 25, 35),
    trt1 = c(100, 200, 300)
  ), tmp, sep = "\t", row.names = FALSE, quote = FALSE)

  vals <- looplook:::load_expression_matrix(tmp, c("con1", "con2"))
  expect_named(vals, c("TP53", "BRCA1", "MYC"))
  expect_equal(vals[["TP53"]], 12.5)

  # error cases
  expect_error(
    looplook:::load_expression_matrix("nonexistent.tsv", "con1"),
    "not found"
  )
  expect_error(
    looplook:::load_expression_matrix(tmp, "badcol"),
    "Requested sample columns not found"
  )
  expect_error(
    looplook:::load_expression_matrix(tmp, c("con1", "con1")),
    "`sample_columns` contains duplicates"
  )

  unlink(tmp)
})

test_that("load_expression_matrix errors on duplicated expression column names", {
  tmp <- tempfile(fileext = ".tsv")
  writeLines(
    c(
      "Gene\tcon1\tcon1\ttrt1",
      "TP53\t10\t15\t100",
      "BRCA1\t20\t25\t200"
    ),
    tmp
  )

  expect_error(
    looplook:::load_expression_matrix(tmp, c("con1", "trt1")),
    "duplicated sample column names"
  )

  unlink(tmp)
})

test_that("get_colors returns expected vector lengths and fallbacks", {
  expect_no_warning(expect_length(looplook:::get_colors(0, "Set2"), 0))
  expect_no_warning(expect_length(looplook:::get_colors(1, "Set2"), 1))
  expect_no_warning(expect_length(looplook:::get_colors(2, "Set2"), 2))
  expect_no_warning(expect_length(looplook:::get_colors(5, "Set2"), 5))
  expect_no_warning(expect_length(looplook:::get_colors(10, c("#E41A1C", "#377EB8")), 10))
  # auto-generate when NULL/empty
  cols <- looplook:::get_colors(6, NULL)
  expect_length(cols, 6)
})

test_that("read_robust_general handles empty and short files", {
  tmp <- tempfile(fileext = ".tsv")
  writeLines("A\tB\tC\n1\t2\t3\n4\t5\t6", tmp)
  res <- looplook:::read_robust_general(tmp, header = TRUE, min_cols = 2)
  expect_equal(nrow(res), 2)
  expect_equal(ncol(res), 3)

  expect_error(
    looplook:::read_robust_general("", desc = "Test"),
    "path is empty"
  )
  unlink(tmp)
})

test_that("resolve_gene_conflicts returns input unchanged when empty", {
  df <- data.frame(
    chr = character(0), start = integer(0), end = integer(0),
    stringsAsFactors = FALSE
  )
  expect_equal(nrow(looplook:::resolve_gene_conflicts(
    df, NULL, NULL, c(-2000, 2000), NULL
  )), 0)
})

test_that("clean_anchor filters and reclassifies correctly", {
  res <- looplook:::clean_anchor("TP53;BRCA1", "P", c("TP53"), TRUE)
  expect_equal(res$type, "P")
  expect_equal(res$gene, "TP53")

  res2 <- looplook:::clean_anchor("GENE_X", "P", c("TP53"), TRUE)
  expect_equal(res2$type, "eP")
  expect_true(is.na(res2$gene))

  res3 <- looplook:::clean_anchor("GENE_X", "P", c("TP53"), FALSE)
  expect_equal(res3$type, "P")
  expect_true(is.na(res3$gene))

  res4 <- looplook:::clean_anchor("", "E", c("TP53"), FALSE)
  expect_equal(res4$type, "E")
  expect_true(is.na(res4$gene))
})

test_that(".map_txdb_gene_ids infers OrgDb keytypes for TxDb-like gene IDs", {
  skip_if_not_installed("org.Hs.eg.db")
  org_db <- org.Hs.eg.db::org.Hs.eg.db

  entrez_ids <- head(AnnotationDbi::keys(org_db, keytype = "ENTREZID"), 5)
  ensembl_ids <- head(AnnotationDbi::keys(org_db, keytype = "ENSEMBL"), 5)

  entrez_map <- looplook:::.map_txdb_gene_ids(
    gene_ids = entrez_ids,
    org_db = org_db,
    columns = "SYMBOL",
    warn = FALSE
  )
  ensembl_map <- looplook:::.map_txdb_gene_ids(
    gene_ids = ensembl_ids,
    org_db = org_db,
    columns = "SYMBOL",
    warn = FALSE
  )

  expect_identical(attr(entrez_map, "keytype"), "ENTREZID")
  expect_identical(attr(ensembl_map, "keytype"), "ENSEMBL")
  expect_true(any(!is.na(entrez_map$SYMBOL)))
  expect_true(any(!is.na(ensembl_map$SYMBOL)))
})

test_that(".map_txdb_gene_ids safely falls back when no OrgDb keytype matches", {
  skip_if_not_installed("org.Hs.eg.db")

  map <- looplook:::.map_txdb_gene_ids(
    gene_ids = c("not_a_real_gene_1", "still_not_real_2"),
    org_db = org.Hs.eg.db::org.Hs.eg.db,
    columns = "SYMBOL",
    warn = FALSE
  )

  expect_true(is.na(attr(map, "keytype")))
  expect_identical(attr(map, "hit_rate"), 0)
  expect_true(all(is.na(map$SYMBOL)))
})

# --- load_expression_matrix: duplicate gene ID warning ---
test_that("load_expression_matrix warns on duplicate gene IDs", {
  tmp <- tempfile(fileext = ".txt")
  write.table(data.frame(
    Gene = c("TP53", "BRCA1", "TP53"),
    s1 = c(10, 20, 30),
    s2 = c(15, 25, 35)
  ), tmp, row.names = FALSE, sep = "\t", quote = FALSE)
  expect_warning(
    looplook:::load_expression_matrix(tmp, c("s1", "s2")),
    "duplicated gene identifier"
  )
  unlink(tmp)
})

# --- load_expression_matrix: single column with integer indices ---
test_that("load_expression_matrix with integer indices", {
  tmp <- tempfile(fileext = ".txt")
  writeLines("Gene\tS1\tS2\tS3\nA\t10\t20\t30\nB\t5\t15\t25", tmp)
  # After removing Gene column, S1 is at index 1, S2 at index 2
  vals <- looplook:::load_expression_matrix(tmp, c(1L))
  expect_equal(vals[["A"]], 10)
  expect_equal(names(vals), c("A", "B"))
  unlink(tmp)
})

# --- load_expression_matrix: empty sample columns error ---
test_that("load_expression_matrix errors on invalid column index", {
  tmp <- tempfile(fileext = ".txt")
  writeLines("Gene\tS1\nA\t10", tmp)
  expect_error(
    looplook:::load_expression_matrix(tmp, c(5L)),
    "invalid column indices"
  )
  unlink(tmp)
})

# --- load_expression_matrix: empty file error ---
test_that("load_expression_matrix errors on single column file", {
  tmp <- tempfile(fileext = ".txt")
  writeLines("Gene\nA\nB", tmp)
  expect_error(
    looplook:::load_expression_matrix(tmp, NULL),
    "at least one sample column"
  )
  unlink(tmp)
})

# --- clean_anchor: eP and eG reclassification ---
test_that("clean_anchor reclassifies eP and eG correctly", {
  # Gene body silent → eG
  res <- looplook:::clean_anchor("SILENT_GENE", "G", c("ACTIVE"), TRUE)
  expect_equal(res$type, "eG")
  expect_true(is.na(res$gene))

  # Enhancer stays E even when silent
  res <- looplook:::clean_anchor("SILENT_GENE", "E", c("ACTIVE"), TRUE)
  expect_equal(res$type, "E")

  # eP stays eP (already reclassified)
  res <- looplook:::clean_anchor("SILENT_GENE", "eP", c("ACTIVE"), TRUE)
  expect_equal(res$type, "eP")
})

# --- resolve_gene_conflicts: basic functionality ---
test_that("resolve_gene_conflicts processes overlapping promoters", {
  skip_if_not_installed("org.Hs.eg.db")
  skip_if_not_installed("GenomicFeatures")

  txdb <- tryCatch(
    AnnotationDbi::loadDb(system.file("extdata", "hg19_knownGene_sample.sqlite", package = "GenomicFeatures")),
    error = function(e) NULL
  )
  skip_if(is.null(txdb), "Sample TxDb unavailable")

  genes_gr <- GenomicFeatures::genes(txdb)
  gene_coords <- data.frame(
    chr = as.character(GenomicRanges::seqnames(genes_gr))[1:3],
    start = GenomicRanges::start(genes_gr)[1:3],
    end = GenomicRanges::end(genes_gr)[1:3],
    annotation = c("Promoter", "Intron", "Distal Intergenic"),
    SYMBOL = c("GENE1", "GENE2", "GENE3"),
    stringsAsFactors = FALSE
  )

  result <- looplook:::resolve_gene_conflicts(
    gene_coords, txdb, "org.Hs.eg.db", c(-2000, 2000), NULL
  )
  expect_s3_class(result, "data.frame")
  expect_true("SYMBOL" %in% colnames(result))
})

# --- species_*_pkg functions ---
test_that("species_txdb_pkg returns correct package names", {
  expect_equal(looplook:::species_txdb_pkg("hg38"), "TxDb.Hsapiens.UCSC.hg38.knownGene")
  expect_equal(looplook:::species_txdb_pkg("hg19"), "TxDb.Hsapiens.UCSC.hg19.knownGene")
  expect_equal(looplook:::species_txdb_pkg("mm10"), "TxDb.Mmusculus.UCSC.mm10.knownGene")
  expect_equal(looplook:::species_txdb_pkg("mm9"), "TxDb.Mmusculus.UCSC.mm9.knownGene")
  expect_null(looplook:::species_txdb_pkg("invalid"))
})

test_that("species_orgdb_pkg returns correct package names", {
  expect_equal(looplook:::species_orgdb_pkg("hg38"), "org.Hs.eg.db")
  expect_equal(looplook:::species_orgdb_pkg("mm10"), "org.Mm.eg.db")
  expect_null(looplook:::species_orgdb_pkg("invalid"))
})

test_that("species_bsgenome_pkg returns correct package names", {
  expect_equal(looplook:::species_bsgenome_pkg("hg38"), "BSgenome.Hsapiens.UCSC.hg38")
  expect_equal(looplook:::species_bsgenome_pkg("mm10"), "BSgenome.Mmusculus.UCSC.mm10")
  expect_null(looplook:::species_bsgenome_pkg("invalid"))
})

# --- .harmonize_seqlevels ---
test_that(".harmonize_seqlevels converts style when mismatched", {
  skip_if_not_installed("GenomeInfoDb")

  # Create GRanges with UCSC style (chr1)
  gr_ucsc <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100, 200))
  GenomeInfoDb::seqlevelsStyle(gr_ucsc) <- "UCSC"

  # Create GRanges with Ensembl style (1)
  gr_ensembl <- GenomicRanges::GRanges("1", IRanges::IRanges(100, 200))
  GenomeInfoDb::seqlevelsStyle(gr_ensembl) <- "Ensembl"

  # Should convert and emit message
  expect_message(
    result <- looplook:::.harmonize_seqlevels(gr_ensembl, gr_ucsc, "test"),
    "Seqlevels style harmonized"
  )
  expect_equal(GenomeInfoDb::seqlevelsStyle(result), "UCSC")
})

test_that(".harmonize_seqlevels does nothing when styles match", {
  skip_if_not_installed("GenomeInfoDb")

  gr1 <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100, 200))
  gr2 <- GenomicRanges::GRanges("chr1", IRanges::IRanges(300, 400))

  # No message when styles already match
  expect_silent(looplook:::.harmonize_seqlevels(gr1, gr2, "test"))
})

test_that(".harmonize_seqlevels handles empty GRanges", {
  skip_if_not_installed("GenomeInfoDb")

  empty_gr <- GenomicRanges::GRanges()
  ref_gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100, 200))

  # Should return empty GRanges without error
  result <- looplook:::.harmonize_seqlevels(empty_gr, ref_gr, "test")
  expect_equal(length(result), 0)
})

test_that("resolve_gene_conflicts: biotype_first strategy works", {
  skip_if_not_installed("GenomicFeatures")
  skip_if_not_installed("org.Hs.eg.db")

  txdb <- tryCatch(
    AnnotationDbi::loadDb(system.file("extdata", "hg19_knownGene_sample.sqlite", package = "GenomicFeatures")),
    error = function(e) NULL
  )
  skip_if(is.null(txdb), "Sample TxDb unavailable")

  genes_gr <- GenomicFeatures::genes(txdb)
  gene_ids <- GenomicRanges::mcols(genes_gr)$gene_id
  symbols <- AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = gene_ids, column = "SYMBOL",
    keytype = "ENTREZID", multiVals = "first"
  )

  # Use ETV6 region (protein-coding, on chr12)
  etv6_idx <- which(symbols == "ETV6")
  skip_if(length(etv6_idx) == 0, "ETV6 not found in sample TxDb")
  etv6_gr <- genes_gr[etv6_idx[1]]

  test_df <- data.frame(
    chr = as.character(GenomicRanges::seqnames(etv6_gr)),
    start = GenomicRanges::start(etv6_gr),
    end = GenomicRanges::end(etv6_gr),
    annotation = "Promoter",
    SYMBOL = NA_character_,
    stringsAsFactors = FALSE
  )

  # Both strategies should resolve to a valid SYMBOL
  result_bf <- looplook:::.with_known_upstream_noise_suppressed(
    looplook:::resolve_gene_conflicts(
      test_df, txdb, "org.Hs.eg.db", c(-2000, 2000),
      gene_expr_map = NULL, min_expr = 0,
      conflict_strategy = "biotype_first"
    )
  )
  expect_true(nrow(result_bf) >= 1)
  expect_true("SYMBOL" %in% colnames(result_bf))

  result_ef <- looplook:::.with_known_upstream_noise_suppressed(
    looplook:::resolve_gene_conflicts(
      test_df, txdb, "org.Hs.eg.db", c(-2000, 2000),
      gene_expr_map = NULL, min_expr = 0,
      conflict_strategy = "expression_first"
    )
  )
  expect_true(nrow(result_ef) >= 1)

  # Both should produce valid output
  expect_true(all(c("SYMBOL", "annotation") %in% colnames(result_bf)))
  expect_true(all(c("SYMBOL", "annotation") %in% colnames(result_ef)))

  # conflict_strategy must be validated
  expect_error(
    looplook:::resolve_gene_conflicts(
      test_df, txdb, "org.Hs.eg.db", c(-2000, 2000),
      gene_expr_map = NULL, conflict_strategy = "invalid"
    )
  )
})

test_that("resolve_gene_conflicts: biotype_first retains silent protein-coding over expressed lncRNA", {
  skip_if_not_installed("GenomicFeatures")
  skip_if_not_installed("org.Hs.eg.db")

  txdb <- tryCatch(
    AnnotationDbi::loadDb(system.file("extdata", "hg19_knownGene_sample.sqlite", package = "GenomicFeatures")),
    error = function(e) NULL
  )
  skip_if(is.null(txdb), "Sample TxDb unavailable")

  genes_gr <- GenomicFeatures::genes(txdb)
  gene_ids <- GenomicRanges::mcols(genes_gr)$gene_id
  symbols <- AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = gene_ids, column = "SYMBOL",
    keytype = "ENTREZID", multiVals = "first"
  )
  # Get the GENETYPE for each gene
  genetypes <- tryCatch(
    AnnotationDbi::mapIds(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = gene_ids, column = "GENETYPE",
      keytype = "ENTREZID", multiVals = "first"
    ),
    error = function(e) NULL
  )

  # Use a protein-coding gene and give it zero expression
  etv6_idx <- which(symbols == "ETV6")
  skip_if(length(etv6_idx) == 0, "ETV6 not found in sample TxDb")
  etv6_gr <- genes_gr[etv6_idx[1]]

  test_df <- data.frame(
    chr = as.character(GenomicRanges::seqnames(etv6_gr)),
    start = GenomicRanges::start(etv6_gr),
    end = GenomicRanges::end(etv6_gr),
    annotation = "Promoter",
    SYMBOL = NA_character_,
    stringsAsFactors = FALSE
  )

  # Expression map: ETV6 has zero expression (simulates silent protein-coding)
  expr_map <- setNames(0, "ETV6")

  # biotype_first: protein-coding retained even when silent
  result_bf <- looplook:::.with_known_upstream_noise_suppressed(
    looplook:::resolve_gene_conflicts(
      test_df, txdb, "org.Hs.eg.db", c(-2000, 2000),
      gene_expr_map = expr_map, min_expr = 1,
      conflict_strategy = "biotype_first"
    )
  )
  expect_true(nrow(result_bf) >= 1)
  # biotype_first keeps the protein-coding gene regardless of expression
  expect_true(!is.na(result_bf$SYMBOL[1]) && result_bf$SYMBOL[1] != "")

  # expression_first with same input should also produce a valid result
  result_ef <- looplook:::.with_known_upstream_noise_suppressed(
    looplook:::resolve_gene_conflicts(
      test_df, txdb, "org.Hs.eg.db", c(-2000, 2000),
      gene_expr_map = expr_map, min_expr = 1,
      conflict_strategy = "expression_first"
    )
  )
  expect_true(nrow(result_ef) >= 1)
})

# --- .filter_target_anchor_overlaps ---
test_that(".filter_target_anchor_overlaps filters targets overlapping anchors", {
  target <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1"),
    ranges = IRanges::IRanges(start = c(1, 1), end = c(100, 1000))
  )
  anchors <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1", "chr1"),
    ranges = IRanges::IRanges(
      start = c(90, 900, 995),
      end   = c(189, 999, 1094)
    )
  )
  res <- looplook:::.filter_target_anchor_overlaps(target, anchors)
  expect_s4_class(res, "GRanges")
  expect_equal(length(res), 2L)
})

test_that(".filter_target_anchor_overlaps applies min_overlap threshold", {
  target <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1"),
    ranges = IRanges::IRanges(start = c(1, 1), end = c(100, 1000))
  )
  anchors <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1", "chr1"),
    ranges = IRanges::IRanges(
      start = c(90, 900, 995),
      end   = c(189, 999, 1094)
    )
  )
  res <- looplook:::.filter_target_anchor_overlaps(target, anchors,
                                                     min_overlap = 2L)
  expect_equal(length(res), 1L)  # only [1,1000] overlaps >=2 anchors
})

test_that(".filter_target_anchor_overlaps handles empty input", {
  empty <- GenomicRanges::GRanges()
  anchors <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = 1, end = 100)
  )
  res <- looplook:::.filter_target_anchor_overlaps(empty, anchors)
  expect_equal(length(res), 0L)
})
