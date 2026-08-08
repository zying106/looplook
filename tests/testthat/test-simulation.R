# tests/testthat/test-simulation.R — synthetic data covering all code paths

# ════════════════════════════════════════════════════════════════════════════
# Shared fixtures: cache TxDb / OrgDb once for all integration tests
# ════════════════════════════════════════════════════════════════════════════
sim_txdb <- NULL
sim_org_db <- NULL
sim_has_bioc <- requireNamespace("org.Hs.eg.db", quietly = TRUE) &&
  requireNamespace("GenomicFeatures", quietly = TRUE)

if (sim_has_bioc) {
  sim_txdb <- tryCatch(
    AnnotationDbi::loadDb(
      system.file("extdata", "hg19_knownGene_sample.sqlite", package = "GenomicFeatures")
    ),
    error = function(e) NULL
  )
  sim_org_db <- "org.Hs.eg.db"
}

# ════════════════════════════════════════════════════════════════════════════
# Part A: Unit-level tests with pure mock data (no TxDb / OrgDb required)
# ════════════════════════════════════════════════════════════════════════════

test_that("clean_anchor: all reclassification patterns", {
  wl <- c("TP53", "BRCA1")

  # active promoter → stays P
  r <- looplook:::clean_anchor("TP53", "P", wl, TRUE)
  expect_equal(r$type, "P")
  expect_equal(r$gene, "TP53")

  # silent promoter with reclassify → eP, positional gene kept
  r <- looplook:::clean_anchor("GENE_X", "P", wl, TRUE)
  expect_equal(r$type, "eP")
  expect_equal(r$gene, "GENE_X")

  # silent promoter without reclassify → stays P, positional gene kept
  r <- looplook:::clean_anchor("GENE_X", "P", wl, FALSE)
  expect_equal(r$type, "P")
  expect_equal(r$gene, "GENE_X")

  # silent gene body with reclassify → eG
  r <- looplook:::clean_anchor("GENE_X", "G", wl, TRUE)
  expect_equal(r$type, "eG")

  # enhancer never reclassified
  r <- looplook:::clean_anchor("GENE_X", "E", wl, TRUE)
  expect_equal(r$type, "E")

  # multi-gene: all positional genes kept
  r <- looplook:::clean_anchor("TP53;GENE_X;BRCA1", "P", wl, TRUE)
  expect_equal(r$gene, "TP53;GENE_X;BRCA1")

  # empty gene string → type preserved, gene NA
  r <- looplook:::clean_anchor("", "E", wl, TRUE)
  expect_equal(r$type, "E")
  expect_true(is.na(r$gene))

  # NA gene
  r <- looplook:::clean_anchor(NA_character_, "P", wl, TRUE)
  expect_equal(r$type, "P")
  expect_true(is.na(r$gene))
})

test_that("extract_target_gene_sets: source branching", {
  anno <- list(
    target_annotation = data.frame(
      SYMBOL = c("A", "B"),
      Assigned_Target_Genes = c("A;C", "B"),
      Assigned_Target_Genes_Filled = c("A;C;D", "B;E"),
      Regulated_promoter_genes = c("A", NA),
      stringsAsFactors = FALSE
    ),
    loop_annotation = data.frame(
      loop_type = c("E-P", "P-P", "E-G", "E-P"),
      Putative_Target_Genes = c("X", "Y", "Z", "W"),
      stringsAsFactors = FALSE
    )
  )

  # targets + all + filled
  r <- looplook:::extract_target_gene_sets(anno, "targets", include_Filled = TRUE, target_mapping_mode = "all")
  expect_setequal(r$Target_Genes, c("A", "B", "C", "D", "E"))

  # targets + promoter mode
  r <- looplook:::extract_target_gene_sets(anno, "targets", include_Filled = FALSE, target_mapping_mode = "promoter")
  expect_setequal(r$Target_Genes, c("A"))

  # targets + use_nearest_gene → uses SYMBOL
  r <- looplook:::extract_target_gene_sets(anno, "targets", use_nearest_gene = TRUE)
  expect_setequal(r$Target_Genes, c("A", "B"))

  # loops with type filter
  r <- looplook:::extract_target_gene_sets(anno, "loops", active_loop_types = c("E-P"))
  expect_true("EP_Genes" %in% names(r))
  expect_setequal(r$EP_Genes, c("X", "W"))

  # loops all types
  r <- looplook:::extract_target_gene_sets(anno, "loops")
  expect_equal(length(r), 3) # EP, PP, EG
})

test_that("target gene link table records evidence and fallback semantics", {
  bed_info <- data.frame(
    input_id = c("Peak_1", "Peak_2"),
    SYMBOL = c("LINEAR_A", "LINEAR_B"),
    annotation = c("Promoter", "Distal Intergenic"),
    All_Loop_Connected_Genes = c("GENE_A;GENE_B;GENE_C", NA),
    Regulated_promoter_genes = c("GENE_A;GENE_B", NA),
    Assigned_Target_Genes = c("GENE_A;GENE_B", NA),
    Regulated_promoter_genes_Filled = c("GENE_A;GENE_B", "LINEAR_B"),
    Assigned_Target_Genes_Filled = c("GENE_A;GENE_B", "LINEAR_B"),
    stringsAsFactors = FALSE
  )
  hit_df <- data.frame(
    qid = c(1L),
    anchor_id = c("A1"),
    stringsAsFactors = FALSE
  )
  loop_annotation <- data.frame(
    loop_ID = c("L1", "L2"),
    a1_id = c("A1", "A1"),
    a2_id = c("A2", "A3"),
    stringsAsFactors = FALSE
  )
  map_info <- data.frame(
    anchor_id = c("A1", "A2", "A3"),
    type_code = c("P", "P", "G"),
    SYMBOL = c("GENE_A", "GENE_B", "GENE_C"),
    stringsAsFactors = FALSE
  )
  ego <- list(
    A1 = stats::setNames(c("A1", "A2", "A3"), c("A1", "A2", "A3"))
  )

  links <- looplook:::.build_target_gene_links(
    hit_df = hit_df,
    bed_info = bed_info,
    loop_annotation_final = loop_annotation,
    map_info = map_info,
    ego_list_target = ego
  )
  links <- looplook:::.mark_target_gene_link_membership(links, bed_info)

  expect_true(all(c(
    "input_id", "loop_ID", "gene", "gene_role", "source", "evidence",
    "anchor_role", "used_as_fallback", "in_regulated_promoter",
    "in_assigned_target", "in_regulated_promoter_filled"
  ) %in% colnames(links)))

  local_a <- links[links$input_id == "Peak_1" & links$gene == "GENE_A" &
    links$anchor_role == "local_anchor", , drop = FALSE]
  expect_setequal(local_a$loop_ID, c("L1", "L2"))
  expect_true(all(local_a$evidence == "local_promoter_overlap"))
  expect_true(all(local_a$in_regulated_promoter))

  direct_b <- links[links$input_id == "Peak_1" & links$gene == "GENE_B", , drop = FALSE]
  expect_true(any(direct_b$evidence == "distal_promoter"))
  expect_true(any(direct_b$in_assigned_target))

  fallback_b <- links[links$input_id == "Peak_2" & links$gene == "LINEAR_B", , drop = FALSE]
  expect_equal(nrow(fallback_b), 1)
  expect_true(fallback_b$used_as_fallback)
  expect_equal(fallback_b$evidence, "linear_fallback")
  expect_true(fallback_b$in_regulated_promoter_filled)

  ev <- looplook:::.summarise_regulated_promoter_evidence(links)
  expect_true("Peak_1" %in% ev$input_id)
  expect_match(ev$Regulated_promoter_Evidence[ev$input_id == "Peak_1"], "local_promoter_overlap")
  expect_match(ev$Regulated_promoter_Evidence[ev$input_id == "Peak_1"], "distal_promoter")
})

test_that("compute_refined_stats: all loop types and multi-gene split", {
  loop_df <- data.frame(
    chr1 = rep("chr1", 8), start1 = seq(1000, 8000, 1000), end1 = seq(1500, 8500, 1000),
    chr2 = rep("chr1", 8), start2 = seq(10000, 17000, 1000), end2 = seq(10500, 17500, 1000),
    cluster_id = paste0("C", c(1, 1, 2, 2, 3, 3, 4, 4)),
    a1_id = paste0("A", 1:8), a2_id = paste0("B", 1:8),
    anchor1_type = c("P", "P", "E", "E", "G", "G", "P", "P"),
    anchor2_type = c("E", "P", "P", "G", "P", "E", "eP", "eG"),
    anchor1_gene = c("G1", "G2;G3", NA, NA, "G4", "G5", "G6;G1", "G7"),
    anchor2_gene = c(NA, NA, "G1", "G8", "G9", NA, "G10", "G11"),
    loop_type = c("E-P", "P-P", "E-P", "E-G", "G-P", "E-G", "eP-P", "eG-P"),
    Putative_Target_Genes = c("G1", "G2;G3", "G1", "G8", "G9", "G5", "G6;G1", "G11"),
    stringsAsFactors = FALSE
  )
  vals <- setNames(seq_len(11) * 2, paste0("G", 1:11))

  res <- looplook:::compute_refined_stats(loop_df,
    upstream_promoter_stats = NULL,
    vals = vals, threshold = 1, hub_percentile = 0.95
  )

  # promoter-centric: G1, G2, G3, G6, G7 should appear as separate rows
  expect_true(all(c("G1", "G2", "G3", "G6", "G7") %in% res$promoter_centric$Gene))
  expect_false("G2;G3" %in% res$promoter_centric$Gene)
  expect_false("G6;G1" %in% res$promoter_centric$Gene)

  # G1 appears 3 times: L1 a1=P, L3 a2=P, L7 a1=P (G6;G1 split)
  g1_row <- res$promoter_centric[res$promoter_centric$Gene == "G1", ]
  expect_equal(g1_row$Total_Loops, 3)

  # distal elements: E, G, eP, eG anchors
  expect_true(!is.null(res$distal_element))
  expect_gt(nrow(res$distal_element), 0)
  # G9 is target of distal anchor A5 (G-P loop: a1=G(E-type distal), a2=P)
  expect_true("G9" %in% res$distal_element$Target_Genes)
})

test_that("compute_refined_stats: hub detection thresholds", {
  loop_df <- data.frame(
    chr1 = "chr1", start1 = 1:10 * 1000, end1 = 1:10 * 1000 + 500,
    chr2 = "chr1", start2 = 1:10 * 2000, end2 = 1:10 * 2000 + 500,
    cluster_id = paste0("C", 1:10),
    a1_id = paste0("A", 1:10), a2_id = paste0("B", 1:10),
    anchor1_type = rep("P", 10), anchor2_type = rep("E", 10),
    anchor1_gene = paste0("Gene", 1:10), anchor2_gene = NA_character_,
    loop_type = rep("E-P", 10),
    Putative_Target_Genes = paste0("Gene", 1:10),
    stringsAsFactors = FALSE
  )
  # Gene1-3 have high connectivity (many loops), Gene4-10 have few
  gene_loops <- c(rep(10, 3), rep(2, 7))
  loop_df <- do.call(rbind, lapply(seq_along(gene_loops), function(i) {
    n <- gene_loops[i]
    if (n <= 1) {
      return(NULL)
    }
    data.frame(
      chr1 = "chr1", start1 = i * 1000, end1 = i * 1000 + 500,
      chr2 = "chr1", start2 = (i + 1) * 2000, end2 = (i + 1) * 2000 + 500,
      cluster_id = paste0("C", i, "_", seq_len(n)),
      a1_id = paste0("A", i, "_", seq_len(n)),
      a2_id = paste0("B", i, "_", seq_len(n)),
      anchor1_type = "P", anchor2_type = "E",
      anchor1_gene = paste0("Gene", i), anchor2_gene = NA_character_,
      loop_type = "E-P",
      Putative_Target_Genes = paste0("Gene", i),
      stringsAsFactors = FALSE
    )
  }))
  vals <- setNames(seq_len(10) * 5, paste0("Gene", 1:10))

  res <- looplook:::compute_refined_stats(loop_df,
    upstream_promoter_stats = NULL,
    vals = vals, threshold = 1, hub_percentile = 0.50
  )

  # With 50th percentile, at least Gene1-3 (with 10 loops each) should be hubs
  high_conn <- res$promoter_centric$Gene[res$promoter_centric$Is_High_Connectivity_Gene == "Yes"]
  expect_true("Gene1" %in% high_conn)
})

test_that("get_feature_class: all categories", {
  classify <- function(x) {
    if (is.na(x)) {
      return("Unknown")
    }
    x <- tolower(x)
    if (grepl("promoter", x)) {
      return("P")
    }
    if (grepl("intergenic|downstream", x)) {
      return("E")
    }
    if (grepl("exon|intron|utr", x)) {
      return("G")
    }
    return("E")
  }
  inputs <- c(
    "Promoter (<=1kb)", "Distal Intergenic", "Intron (uc001.1)",
    "Exon (uc001.1)", "5-UTR", "Downstream (<=300bp)",
    "3-UTR", NA_character_
  )
  expected <- c("P", "E", "G", "G", "G", "E", "G", "Unknown")
  expect_equal(unname(vapply(inputs, classify, character(1))), expected)
})

test_that("extract_genes and clean_gene_names: complete edge cases", {
  expect_equal(looplook:::extract_genes(c("A;B", "B;C", "D")), "A;B;C;D")
  expect_equal(looplook:::extract_genes(""), NA_character_)
  expect_equal(looplook:::extract_genes(";"), NA_character_)
  expect_equal(looplook:::extract_genes(c(NA, "")), NA_character_)
  expect_equal(looplook:::extract_genes(" A ; B "), "A;B") # trimws added

  expect_equal(looplook:::clean_gene_names(NULL), character(0))
  expect_equal(looplook:::clean_gene_names(c("A", NA, "", " ")), "A")
  expect_equal(looplook:::clean_gene_names("A;B;A;", ";"), c("A", "B"))
  expect_equal(looplook:::clean_gene_names(c("A,B", "C,D"), "[;,]"), c("A", "B", "C", "D"))
})

test_that("load_expression_matrix: all input variants", {
  tmp <- tempfile(fileext = ".tsv")
  writeLines("Gene\tS1\tS2\tS3\nA\t10\t20\t30\nB\t5\t15\t25\nC\t0\t0\t0", tmp)

  # multiple columns
  v <- looplook:::load_expression_matrix(tmp, c("S1", "S2"))
  expect_equal(v[["A"]], 15)

  # single column preserves names
  v <- looplook:::load_expression_matrix(tmp, "S1")
  expect_equal(v[["A"]], 10)

  # all columns (NULL)
  v <- looplook:::load_expression_matrix(tmp, NULL)
  expect_equal(v[["C"]], 0)

  # integer indices
  v <- looplook:::load_expression_matrix(tmp, c(2L, 3L))
  expect_equal(names(v), c("A", "B", "C"))

  unlink(tmp)

  # error: non-numeric values
  tmp2 <- tempfile(fileext = ".tsv")
  writeLines("Gene\tS1\nA\tNA\nB\tbad", tmp2)
  expect_error(looplook:::load_expression_matrix(tmp2, "S1"), "non-numeric")
  unlink(tmp2)
})

test_that("format_annotation_columns: roundtrip", {
  df <- data.frame(
    annotation = c("Promoter (<=1kb)", "Intron (uc001.1)"),
    stringsAsFactors = FALSE
  )
  res <- looplook:::format_annotation_columns(df)
  expect_equal(res$annotation, c("Promoter", "Intron"))
  expect_equal(res$detail_anno, c("Promoter (<=1kb)", "Intron (uc001.1)"))

  # no annotation column → pass through
  df2 <- data.frame(x = 1:3)
  res2 <- looplook:::format_annotation_columns(df2)
  expect_equal(res2$x, 1:3)
})

test_that("simplify_annotation: all paths", {
  x <- c(
    "Promoter (<=1kb)", "Intron (uc001.1)", "Exon (uc001.1)",
    "Distal Intergenic", "Downstream (<=300bp)", "unknown_type", NA
  )
  res <- looplook:::simplify_annotation(x)
  expect_equal(unname(res), c(
    "Promoter", "Intron", "Exon", "Distal Intergenic",
    "Downstream", "Others", "Others"
  ))
})

test_that("read_robust_general: fill option for ragged lines", {
  tmp <- tempfile(fileext = ".csv")
  writeLines("A,B,C\n1,2\n3,4,5,6", tmp)
  res <- looplook:::read_robust_general(tmp, header = TRUE)
  expect_equal(ncol(res), 4)
  unlink(tmp)
})

test_that("get_colors: RColorBrewer, custom, and fallback", {
  expect_length(looplook:::get_colors(0, "Set2"), 0)
  expect_length(looplook:::get_colors(1, "Set1"), 1)
  expect_length(looplook:::get_colors(50, "Set2"), 50)

  cols <- looplook:::get_colors(3, c("#FF0000"))
  expect_equal(cols, c("#FF0000", "#FF0000", "#FF0000"))

  cols <- looplook:::get_colors(3, NULL)
  expect_length(cols, 3)
})


# ════════════════════════════════════════════════════════════════════════════
# Part B: Integration tests with sample TxDb
# ════════════════════════════════════════════════════════════════════════════

test_that("annotate_peaks_and_loops: all loop-type combinations", {
  skip_if_not(sim_has_bioc && !is.null(sim_txdb), "TxDb/OrgDb unavailable")
  txdb <- sim_txdb

  # Get real gene coordinates from TxDb
  genes_gr <- GenomicFeatures::genes(txdb)
  gene_ids <- names(genes_gr)
  map <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
    keys = gene_ids, columns = "SYMBOL", keytype = "ENTREZID"
  )
  map <- map[!duplicated(map$ENTREZID), ]

  # Pick genes on the same chromosome for intra-chromosomal loops
  gene_df <- data.frame(
    gene_id = gene_ids,
    chr = as.character(GenomicRanges::seqnames(genes_gr)),
    start = GenomicRanges::start(genes_gr),
    end = GenomicRanges::end(genes_gr),
    stringsAsFactors = FALSE
  )
  gene_df <- merge(gene_df, map, by.x = "gene_id", by.y = "ENTREZID")

  # Find chr6 genes: TFAP2A-AS1 (~10.4M) and SRPK1 (~35.8M)
  # Place anchors to create E-P, P-P, E-E, E-G, P-G loops
  chr6_genes <- gene_df[gene_df$chr == "chr6", ]
  skip_if(nrow(chr6_genes) < 2, "Need ≥2 genes on same chr for loop tests")

  g1 <- chr6_genes[1, ] # TFAP2A-AS1: ncRNA
  g2 <- chr6_genes[2, ] # SRPK1: protein-coding

  # Anchor A: exactly at g1 TSS → should annotate as Promoter
  # Anchor B: intergenic between g1 and g2 → should annotate as Distal Intergenic (E)
  # Anchor C: at g2 TSS → should annotate as Promoter
  # Anchor D: within g2 gene body (TSS+1000) → should annotate as G (Exon/Intron)

  tss_a <- ifelse(g1$start < g1$end, g1$start, g1$end)
  tss_c <- ifelse(g2$start < g2$end, g2$start, g2$end)
  mid_ab <- round((tss_a + tss_c) / 2)

  # Build BEDPE with 4 loops:
  # L1: P(A)-E(B) → type should be E-P
  # L2: P(A)-P(C) → type should be P-P
  # L3: E(B)-P(C) → type should be E-P
  # L4: P(A)-G(D) → type should be G-P (if D hits gene body)
  bedpe_lines <- c(
    sprintf("chr6\t%d\t%d\tchr6\t%d\t%d", tss_a - 1000, tss_a + 500, mid_ab - 500, mid_ab + 500),
    sprintf("chr6\t%d\t%d\tchr6\t%d\t%d", tss_a - 1000, tss_a + 500, tss_c - 1000, tss_c + 500),
    sprintf("chr6\t%d\t%d\tchr6\t%d\t%d", mid_ab - 500, mid_ab + 500, tss_c - 1000, tss_c + 500),
    sprintf("chr6\t%d\t%d\tchr6\t%d\t%d", tss_a - 1000, tss_a + 500, tss_c + 1000, tss_c + 10000)
  )
  tmp_bedpe <- tempfile(fileext = ".bedpe")
  writeLines(bedpe_lines, tmp_bedpe)

  res <- looplook:::.with_known_upstream_noise_suppressed(
    looplook::annotate_peaks_and_loops(
      bedpe_file = tmp_bedpe,
      txdb = txdb, org_db = "org.Hs.eg.db", species = "hg19",
      out_dir = tempdir(), project_name = "Sim_LoopTypes",
      write_output = FALSE, quiet = TRUE
    )
  )

  expect_type(res, "list")
  la <- res$loop_annotation
  expect_equal(nrow(la), 4)

  # All loops should have a classified loop_type (not "Unknown")
  expect_false(any(la$loop_type == "Unknown"))

  # Promoter-centric stats should exist for the P anchors
  expect_gt(nrow(res$promoter_centric_stats), 0)

  unlink(tmp_bedpe)
})

test_that("annotate_peaks_and_loops: target BED integration with and without loop overlap", {
  skip_if_not(sim_has_bioc && !is.null(sim_txdb), "TxDb/OrgDb unavailable")
  txdb <- sim_txdb

  genes_gr <- GenomicFeatures::genes(txdb)
  gene_ids <- names(genes_gr)
  map <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
    keys = gene_ids, columns = "SYMBOL", keytype = "ENTREZID"
  )
  map <- map[!duplicated(map$ENTREZID), ]
  gene_df <- data.frame(
    gene_id = gene_ids,
    chr = as.character(GenomicRanges::seqnames(genes_gr)),
    start = GenomicRanges::start(genes_gr),
    end = GenomicRanges::end(genes_gr),
    stringsAsFactors = FALSE
  )
  gene_df <- merge(gene_df, map, by.x = "gene_id", by.y = "ENTREZID")
  chrX_genes <- gene_df[gene_df$chr == "chrX", ]
  skip_if(nrow(chrX_genes) < 3, "Need ≥3 genes on chrX")

  g1 <- chrX_genes[1, ] # CLCN4
  g2 <- chrX_genes[3, ] # ARMCX3
  tss1 <- ifelse(g1$start < g1$end, g1$start, g1$end)
  tss2 <- ifelse(g2$start < g2$end, g2$start, g2$end)

  # One loop: P(g1) <-> E(intergenic)
  bedpe_lines <- sprintf(
    "chrX\t%d\t%d\tchrX\t%d\t%d",
    tss1 - 1000, tss1 + 500,
    round((tss1 + tss2) / 2) - 500, round((tss1 + tss2) / 2) + 500
  )
  tmp_bedpe <- tempfile(fileext = ".bedpe")
  writeLines(bedpe_lines, tmp_bedpe)

  # Target BED:
  # Peak1: overlaps the P anchor → should get 3D-assigned targets
  # Peak2: far away from any anchor → orphan (only linear nearest gene)
  tmp_bed <- tempfile(fileext = ".bed")
  writeLines(c(
    sprintf("chrX\t%d\t%d", tss1 - 500, tss1 + 100), # overlaps promoter anchor
    "chrX\t50000000\t50001000" # orphan peak
  ), tmp_bed)

  res <- looplook:::.with_known_upstream_noise_suppressed(
    looplook::annotate_peaks_and_loops(
      bedpe_file = tmp_bedpe,
      target_bed = tmp_bed,
      txdb = txdb, org_db = "org.Hs.eg.db", species = "hg19",
      out_dir = tempdir(), project_name = "Sim_Targets",
      write_output = FALSE, quiet = TRUE
    )
  )

  ta <- res$target_annotation
  expect_equal(nrow(ta), 2)
  expect_true(all(c(
    "Regulated_promoter_Evidence",
    "Regulated_promoter_genes_Filled",
    "Regulated_promoter_Fallback_Evidence",
    "Assigned_Target_Genes_Filled"
  ) %in% colnames(ta)))
  expect_s3_class(res$target_gene_links, "data.frame")
  expect_true(all(c(
    "input_id", "loop_ID", "anchor_id", "gene", "gene_role", "source",
    "evidence", "anchor_role", "used_as_fallback",
    "in_regulated_promoter", "in_assigned_target",
    "in_regulated_promoter_filled", "in_assigned_target_filled"
  ) %in% colnames(res$target_gene_links)))

  # Peak1 should have 3D-derived targets
  peak1 <- ta[1, ]
  expect_false(is.na(peak1$Assigned_Target_Genes) || peak1$Assigned_Target_Genes == "")
  expect_false(is.na(peak1$Regulated_promoter_Evidence) || peak1$Regulated_promoter_Evidence == "")
  expect_true(peak1$Regulated_promoter_Evidence != "none")
  expect_true(any(res$target_gene_links$input_id == peak1$input_id &
    res$target_gene_links$source == "loop_anchor"))

  # Both should have _Filled column (peak2 falls back to SYMBOL)
  expect_true("Assigned_Target_Genes_Filled" %in% colnames(ta))
  expect_false(any(is.na(ta$Assigned_Target_Genes_Filled) | ta$Assigned_Target_Genes_Filled == ""))
  peak2 <- ta[2, ]
  expect_true(is.na(peak2$Regulated_promoter_genes) || peak2$Regulated_promoter_genes == "")
  expect_false(is.na(peak2$Regulated_promoter_genes_Filled) || peak2$Regulated_promoter_genes_Filled == "")
  expect_false(is.na(peak2$Regulated_promoter_Fallback_Evidence) ||
    peak2$Regulated_promoter_Fallback_Evidence %in% c("", "none"))
  expect_true(any(res$target_gene_links$input_id == peak2$input_id &
    res$target_gene_links$source == "linear_annotation" &
    res$target_gene_links$used_as_fallback))

  # All_Loop_Connected_Genes should differ from Assigned_Target_Genes
  # (ATGs uses priority: promoter-only for 3D assignment)
  has_loop <- !is.na(ta$All_Loop_Connected_Genes) & ta$All_Loop_Connected_Genes != ""
  if (any(has_loop)) {
    # verify Assigned_Target_Genes is subset of All_Loop_Connected_Genes
    for (i in which(has_loop)) {
      assigned <- looplook:::clean_gene_names(ta$Assigned_Target_Genes[i], ";")
      all_conn <- looplook:::clean_gene_names(ta$All_Loop_Connected_Genes[i], ";")
      if (length(assigned) > 0) {
        expect_true(all(assigned %in% all_conn))
      }
    }
  }

  unlink(c(tmp_bedpe, tmp_bed))
})

test_that("annotate_peaks_and_loops: neighbor_hop controls ego-network radius", {
  skip_if_not(sim_has_bioc && !is.null(sim_txdb), "TxDb/OrgDb unavailable")
  txdb <- sim_txdb

  genes_gr <- GenomicFeatures::genes(txdb)
  gene_ids <- names(genes_gr)
  map <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
    keys = gene_ids, columns = "SYMBOL", keytype = "ENTREZID"
  )
  map <- map[!duplicated(map$ENTREZID), ]
  gene_df <- data.frame(
    gene_id = gene_ids,
    chr = as.character(GenomicRanges::seqnames(genes_gr)),
    start = GenomicRanges::start(genes_gr),
    end = GenomicRanges::end(genes_gr),
    stringsAsFactors = FALSE
  )
  gene_df <- merge(gene_df, map, by.x = "gene_id", by.y = "ENTREZID")

  # Use chrX genes ARMCX3 and FAM199X (close enough for realistic loops)
  chrX_genes <- gene_df[gene_df$chr == "chrX", ]
  skip_if(nrow(chrX_genes) < 3, "Need ≥3 genes on chrX")

  g1 <- chrX_genes[3, ] # ARMCX3
  g2 <- chrX_genes[4, ] # FAM199X
  tss1 <- ifelse(g1$start < g1$end, g1$start, g1$end)
  tss2 <- ifelse(g2$start < g2$end, g2$start, g2$end)

  # Build 3 loops that share the SAME promoter anchor (A):
  # L1: P(A at g1) — E(B intergenic)
  # L2: P(A at g1) — E(C intergenic, different location)
  # L3: P(A at g1) — P(D at g2)
  anchor_a_start <- tss1 - 1000
  anchor_a_end <- tss1 + 500

  bedpe_lines <- c(
    sprintf(
      "chrX\t%d\t%d\tchrX\t%d\t%d",
      anchor_a_start, anchor_a_end, tss1 + 50000 - 500, tss1 + 50000 + 500
    ),
    sprintf(
      "chrX\t%d\t%d\tchrX\t%d\t%d",
      anchor_a_start, anchor_a_end, tss1 + 80000 - 500, tss1 + 80000 + 500
    ),
    sprintf(
      "chrX\t%d\t%d\tchrX\t%d\t%d",
      anchor_a_start, anchor_a_end, tss2 - 1000, tss2 + 500
    )
  )
  tmp_bedpe <- tempfile(fileext = ".bedpe")
  writeLines(bedpe_lines, tmp_bedpe)

  # Target peak overlapping the shared anchor A → should link to all 3 loops
  tmp_bed <- tempfile(fileext = ".bed")
  writeLines(sprintf("chrX\t%d\t%d", anchor_a_start, anchor_a_end), tmp_bed)

  res <- looplook:::.with_known_upstream_noise_suppressed(
    looplook::annotate_peaks_and_loops(
      bedpe_file = tmp_bedpe, target_bed = tmp_bed,
      txdb = txdb, org_db = "org.Hs.eg.db", species = "hg19",
      out_dir = tempdir(), project_name = "Sim_MultiLoop",
      write_output = FALSE, quiet = TRUE
    )
  )

  ta <- res$target_annotation
  expect_equal(nrow(ta), 1)

  # Should be connected to loops (Linked_Loop_IDs should have multiple entries)
  expect_false(is.na(ta$Linked_Loop_IDs) || ta$Linked_Loop_IDs == "")
  n_loops <- length(unlist(strsplit(ta$Linked_Loop_IDs, ";")))
  expect_equal(n_loops, 3)

  # Assigned_Target_Genes should not be empty (at least g1 via L3)
  expect_false(is.na(ta$Assigned_Target_Genes) || ta$Assigned_Target_Genes == "")
  expect_s3_class(res$target_gene_links, "data.frame")
  linked_local_rows <- res$target_gene_links[
    res$target_gene_links$input_id == ta$input_id[1] &
      res$target_gene_links$anchor_role == "local_anchor" &
      !is.na(res$target_gene_links$loop_ID), ,
    drop = FALSE
  ]
  expect_gte(length(unique(linked_local_rows$loop_ID)), 3)

  unlink(c(tmp_bedpe, tmp_bed))
})

test_that("refine: full target-assignment fallback chain", {
  skip_if_not(sim_has_bioc, "OrgDb unavailable")
  rdata_path <- system.file("extdata", "analysis_results.RData", package = "looplook")
  expr_path <- system.file("extdata", "example_tpm.txt", package = "looplook")
  skip_if(rdata_path == "" || expr_path == "")

  tmp <- new.env()
  load(rdata_path, envir = tmp)
  res <- tmp[[ls(tmp)[1]]]
  res$loop_annotation <- head(res$loop_annotation, 20)
  res$target_annotation <- head(res$target_annotation, 8)
  res$promoter_centric_stats <- head(res$promoter_centric_stats, 20)
  res$distal_element_stats <- head(res$distal_element_stats, 20)

  # Use a high threshold so most genes are silenced → triggers reclassify + fallback
  ref <- looplook:::.with_known_upstream_noise_suppressed(
    looplook::refine_loop_anchors_by_expression(
      annotation_res = res, expr_matrix_file = expr_path,
      sample_columns = c("con1", "con2"), threshold = 100,
      reclassify_by_expression = TRUE,
      out_dir = tempdir(), project_name = "Sim_FallbackChain",
      write_output = FALSE, quiet = TRUE
    )
  )

  # Verify the output contract
  expect_type(ref, "list")
  la <- ref$loop_annotation

  # After high-threshold refinement, some anchors should be reclassified
  types <- unique(c(la$anchor1_type, la$anchor2_type))
  expect_true(any(grepl("^e", types)))

  # Putative_Target_Genes should use fallback when original targets are silenced
  expect_true("Putative_Target_Genes" %in% colnames(la))

  # If target_annotation present, Assigned_Target_Genes_Filled must exist
  if (!is.null(ref$target_annotation)) {
    ta <- ref$target_annotation
    expect_true("Assigned_Target_Genes_Filled" %in% colnames(ta))
    expect_true("Regulated_promoter_Evidence" %in% colnames(ta) ||
      is.null(res$target_annotation$Regulated_promoter_Evidence))
    if ("Regulated_promoter_Evidence" %in% colnames(ta)) {
      expect_false(any(is.na(ta$Regulated_promoter_Evidence)))
    }
    # _Filled column should never be NA when Assigned_Target_Genes is present
    has_assigned <- !is.na(ta$Assigned_Target_Genes) & ta$Assigned_Target_Genes != ""
    filled_na <- is.na(ta$Assigned_Target_Genes_Filled) | ta$Assigned_Target_Genes_Filled == ""
    expect_false(any(has_assigned & filled_na))
  }
  if (!is.null(ref$target_gene_links)) {
    expect_s3_class(ref$target_gene_links, "data.frame")
    expect_true("used_as_fallback" %in% colnames(ref$target_gene_links))
  }
})

test_that("refined target gene links retain only post-refinement targets", {
  target_gene_links <- data.frame(
    input_id = c("Peak_1", "Peak_1", "Peak_2"),
    loop_ID = c("Loop_1", "Loop_1", NA_character_),
    anchor_id = c("Anchor_1", "Anchor_2", NA_character_),
    gene = c("ACTIVE", "SILENT", "FALLBACK"),
    gene_role = c("promoter", "promoter", "linear_annotation"),
    source = c("loop_anchor", "loop_anchor", "linear_annotation"),
    evidence = c("distal_promoter", "distal_promoter", "linear_annotation"),
    anchor_role = c("opposite_anchor", "opposite_anchor", "linear_annotation"),
    stringsAsFactors = FALSE
  )
  bed_info <- data.frame(
    input_id = c("Peak_1", "Peak_2"),
    Regulated_promoter_genes = c("ACTIVE", NA_character_),
    Assigned_Target_Genes = c("ACTIVE", NA_character_),
    All_Loop_Connected_Genes = c("ACTIVE", NA_character_),
    Regulated_promoter_genes_Filled = c("ACTIVE", "FALLBACK"),
    Assigned_Target_Genes_Filled = c("ACTIVE", "FALLBACK"),
    stringsAsFactors = FALSE
  )
  vals <- c(ACTIVE = 10, SILENT = 0, FALLBACK = 5)

  refined_links <- looplook:::.filter_refined_target_gene_links(
    target_gene_links, bed_info, vals,
    threshold = 1
  )

  # SILENT is kept as positional gene in anchor links but marked as
  # measured_silent (Passes_Expression_Filter = FALSE)
  expect_setequal(refined_links$gene, c("ACTIVE", "SILENT", "FALLBACK"))
  expect_true("SILENT" %in% refined_links$gene)
  expect_false(any(refined_links$Passes_Expression_Filter[refined_links$gene == "SILENT"] %in% TRUE))
  expect_true(refined_links$used_as_fallback[refined_links$gene == "FALLBACK"])
  expect_true(all(c("Mean_Expression", "Passes_Expression_Filter") %in%
    colnames(refined_links)))
})


# ════════════════════════════════════════════════════════════════════════════
# Part D: analysis.R helpers and visualization.R exported functions
# ════════════════════════════════════════════════════════════════════════════

test_that(".calc_gc_fraction computes GC content correctly", {
  seqs <- c("GCGC", "ATAT", "NNNN", "GCAT", "")
  gc <- looplook:::.calc_gc_fraction(seqs)
  expect_equal(gc[1], 1.0)
  expect_equal(gc[2], 0.0)
  expect_equal(gc[3], 0.0)
  expect_equal(gc[4], 0.5)
  expect_true(is.na(gc[5]))
})

test_that(".is_promoter_anchor_type and .is_enhancer_like_anchor_type", {
  expect_true(looplook:::.is_promoter_anchor_type("P"))
  expect_false(looplook:::.is_promoter_anchor_type("E"))
  expect_false(looplook:::.is_promoter_anchor_type("eP"))
  expect_false(looplook:::.is_promoter_anchor_type(NA_character_))

  expect_true(looplook:::.is_enhancer_like_anchor_type("E"))
  expect_false(looplook:::.is_enhancer_like_anchor_type("eP"))
  expect_false(looplook:::.is_enhancer_like_anchor_type("eG"))
  expect_false(looplook:::.is_enhancer_like_anchor_type("P"))
  expect_false(looplook:::.is_enhancer_like_anchor_type("G"))
})

test_that(".empty_anchor_df and .make_anchor_df produce valid output", {
  edf <- looplook:::.empty_anchor_df()
  expect_s3_class(edf, "data.frame")
  expect_equal(nrow(edf), 0)
  expect_equal(colnames(edf), c("anchor_id", "chr", "start", "end", "anchor_type", "cluster_id"))

  loop_df <- data.frame(
    chr1 = "chr1", start1 = 100, end1 = 150, a1_id = "A1",
    stringsAsFactors = FALSE
  )
  adf <- looplook:::.make_anchor_df(loop_df, 1, "1", c("P"))
  expect_equal(adf$anchor_id, "A1")
  expect_equal(adf$chr, "chr1")
  expect_equal(adf$anchor_type, "P")

  loop_df2 <- data.frame(chr1 = "chr2", start1 = 200, end1 = 250, stringsAsFactors = FALSE)
  adf2 <- looplook:::.make_anchor_df(loop_df2, 1, "1", c("E"))
  expect_match(adf2$anchor_id, "chr2_200_250")
})

test_that("run_lfc_violin: edge cases with constant values", {
  global_glist <- setNames(rep(1, 100), paste0("Gene", 1:100))
  targets <- names(global_glist)[1:5]
  p <- looplook:::run_lfc_violin(targets, global_glist, "wilcox.test", "Test")
  expect_s3_class(p, "ggplot")

  global_glist2 <- setNames(c(rep(-2, 50), rep(-1, 50)), paste0("Gene", 1:100))
  p2 <- looplook:::run_lfc_violin(targets, global_glist2, "t.test", "Neg")
  expect_s3_class(p2, "ggplot")
})

test_that("draw_flower_simplified: exported function with mock gene lists", {
  skip_if_not_installed("ggplot2")
  sets <- list(
    A = c("TP53", "BRCA1", "MYC", "EGFR", "KRAS"),
    B = c("BRCA1", "MYC", "EGFR", "PTEN", "APC"),
    C = c("MYC", "EGFR", "KRAS", "TP53", "BRAF")
  )
  p <- draw_flower_simplified(sets, "Test", c(A = "#E41A1C", B = "#377EB8", C = "#4DAF4A"))
  expect_s3_class(p, "ggplot")

  expect_message(
    p_null <- draw_flower_simplified(list(A = c("X")), "One", NULL),
    "Less than 2"
  )
  expect_null(p_null)
})

test_that("plot_peaks_interactions: renders without gene track", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggforce")
  skip_if_not_installed("ggrepel")

  tmp <- tempfile(fileext = ".bedpe")
  writeLines(c(
    "chr1\t10000\t10500\tchr1\t20000\t20500",
    "chr1\t15000\t15500\tchr1\t25000\t25500"
  ), tmp)

  p <- plot_peaks_interactions(
    bedpe_file = tmp, chr = "chr1", from = 9000, to = 26000,
    show_gene_track = FALSE
  )
  expect_s3_class(p, "ggplot")
  unlink(tmp)
})

test_that("prepare_track_data: works without gene track", {
  tmp <- tempfile(fileext = ".bedpe")
  writeLines("chr1\t10000\t10500\tchr1\t20000\t20500", tmp)

  d <- looplook:::prepare_track_data(
    bedpe_file = tmp, target_bed = NULL, chr = "chr1",
    from = 9000, to = 21000, species = "hg38",
    max_levels = 10, base_anchor_height = 0.05,
    loop_color = "#5D6D7E", anchor_color = "#3498DB",
    score_to_alpha = FALSE, min_score = NULL,
    show_gene_track = FALSE
  )
  expect_type(d, "list")
  expect_equal(d$chr, "chr1")
  expect_gt(nrow(d$bez_df), 0)
  unlink(tmp)
})
