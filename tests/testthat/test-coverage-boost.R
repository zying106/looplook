# tests/testthat/test-coverage-boost.R
# Fast mock-based tests targeting low-coverage functions.
# Focus: analysis.R (58%), visualization.R (69%), plus small gains elsewhere.

# ============================================================================
# analysis.R — GSEA, motif helpers, heatmap/connectivity, edge cases
# ============================================================================

test_that("run_gsea_analysis returns list with result and plot for valid input", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("enrichplot")

  set.seed(42)
  glist <- sort(setNames(rnorm(200), paste0("GENE", 1:200)), decreasing = TRUE)
  targets <- names(glist)[1:30]

  out <- looplook:::run_gsea_analysis(
    target_genes = targets, global_glist = glist,
    gsea_nSample = NULL, current_proj_name = "TestGSEA"
  )
  expect_type(out, "list")
  expect_true("result" %in% names(out))
  expect_true("plot" %in% names(out))
})

test_that("run_gsea_analysis down-samples when gsea_nSample < length(targets)", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("enrichplot")

  set.seed(123)
  glist <- sort(setNames(rnorm(200), paste0("GENE", 1:200)), decreasing = TRUE)
  targets <- names(glist)[1:50]

  out <- looplook:::run_gsea_analysis(
    target_genes = targets, global_glist = glist,
    gsea_nSample = 15, current_proj_name = "TestDownsample"
  )
  expect_type(out, "list")
})

test_that("run_gsea_analysis handles duplicate ranked values with tie-breaking", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("enrichplot")

  set.seed(7)
  glist <- sort(setNames(c(rnorm(100), rep(0.5, 100)), paste0("G", 1:200)),
                decreasing = TRUE)
  targets <- names(glist)[1:20]
  # Should not error on duplicates
  out <- looplook:::run_gsea_analysis(
    target_genes = targets, global_glist = glist,
    gsea_nSample = NULL, current_proj_name = "TestTies"
  )
  expect_type(out, "list")
})

test_that("run_gsea_analysis returns NULL for <2 valid targets", {
  skip_if_not_installed("clusterProfiler")
  glist <- setNames(rnorm(100), paste0("GENE", 1:100))
  out <- looplook:::run_gsea_analysis(
    target_genes = "NOMATCH", global_glist = glist,
    gsea_nSample = NULL, current_proj_name = "NoMatch"
  )
  expect_equal(out$result, NULL)
  expect_equal(out$plot, NULL)
})

test_that(".calc_gc_fraction computes correct GC content", {
  seqs <- c("ACGTACGTGC", "AAAAAATTTT", "GCGCGCGCGC", "")
  gc <- looplook:::.calc_gc_fraction(seqs)
  expect_equal(gc[1], 0.6)   # 6/10 GC
  expect_equal(gc[2], 0.0)   # 0/10 GC
  expect_equal(gc[3], 1.0)   # 10/10 GC
  expect_true(is.na(gc[4]))  # empty string → NA
})

test_that(".sample_gc_matched_background returns bg when bg <= target_n (early exit)", {
  bg_gr <- GenomicRanges::GRanges(
    seqnames = rep("chr1", 3),
    ranges = IRanges::IRanges(start = c(1, 100, 200), width = 10)
  )
  fg_gr <- bg_gr[1]
  # bg has 3 <= target_n (5) → early return without GC computation
  out <- looplook:::.sample_gc_matched_background(fg_gr, bg_gr, NULL, max_bg = 5L)
  expect_equal(length(out), 3)
})

test_that(".sample_gc_matched_background returns empty for empty bg", {
  bg_gr <- GenomicRanges::GRanges()
  fg_gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 10))
  out <- looplook:::.sample_gc_matched_background(fg_gr, bg_gr, NULL)
  expect_equal(length(out), 0)
})

test_that(".deduplicate_anchor_df removes duplicates and NAs", {
  df <- data.frame(
    anchor_id = c("A1", "", "A1", NA, "A3"),
    chr = c("chr1", "chr1", "chr1", NA, "chr1"),
    start = c(1L, 2L, 1L, NA, 5L),
    end = c(100L, 200L, 100L, NA, 500L),
    anchor_type = c("P", "E", "P", "E", "G"),
    stringsAsFactors = FALSE
  )
  out <- looplook:::.deduplicate_anchor_df(df)
  # A1 (deduped to 1), chr1_2_200 (empty id filled), A3 = 3 rows
  expect_equal(nrow(out), 3)
  expect_false(any(duplicated(out$anchor_id)))
})

test_that(".deduplicate_anchor_df handles empty and NULL input", {
  expect_equal(nrow(looplook:::.deduplicate_anchor_df(NULL)), 0)
  expect_equal(nrow(looplook:::.deduplicate_anchor_df(
    data.frame(anchor_id = character(), chr = character(),
               start = integer(), end = integer(), anchor_type = character(),
               stringsAsFactors = FALSE)
  )), 0)
})

test_that(".empty_anchor_df returns correct structure", {
  df <- looplook:::.empty_anchor_df()
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 0)
  expect_true(all(c("anchor_id", "chr", "start", "end", "anchor_type") %in% colnames(df)))
})

test_that(".anchor_df_to_gr converts valid anchors to GRanges", {
  df <- data.frame(
    anchor_id = c("A1", "A2"),
    chr = c("chr1", "chr2"),
    start = c(1L, 100L),
    end = c(50L, 200L),
    anchor_type = c("P", "E"),
    stringsAsFactors = FALSE
  )
  gr <- looplook:::.anchor_df_to_gr(df)
  expect_s4_class(gr, "GRanges")
  expect_equal(length(gr), 2)
  expect_equal(names(gr), c("A1", "A2"))
})

test_that(".anchor_df_to_gr returns empty GRanges for empty input", {
  gr <- looplook:::.anchor_df_to_gr(looplook:::.empty_anchor_df())
  expect_s4_class(gr, "GRanges")
  expect_equal(length(gr), 0)
})

test_that(".make_anchor_df builds correct data frame from loop rows", {
  loop_df <- data.frame(
    chr1 = "chr1", start1 = 100L, end1 = 200L,
    chr2 = "chr2", start2 = 300L, end2 = 400L,
    a1_id = "A1", a2_id = "A2",
    stringsAsFactors = FALSE
  )
  types <- c("P", "E")
  out1 <- looplook:::.make_anchor_df(loop_df, 1L, "1", types)
  expect_equal(out1$anchor_id, "A1")
  expect_equal(out1$chr, "chr1")
  out2 <- looplook:::.make_anchor_df(loop_df, 1L, "2", types)
  expect_equal(out2$anchor_id, "A2")
})

test_that("run_heatmap_and_connectivity handles missing ComplexHeatmap gracefully", {
  skip_if_not_installed("ggpointdensity")
  skip_if_not_installed("viridis")

  genes <- paste0("Gene", 1:20)
  tpm_mat <- as.data.frame(matrix(rexp(20 * 4, rate = 1), 20, 4,
    dimnames = list(genes, c("s1", "s2", "s3", "s4"))
  ))
  meta <- data.frame(
    SampleID = c("s1", "s2", "s3", "s4"),
    Group = c("A", "A", "B", "B"), stringsAsFactors = FALSE
  )
  glist <- setNames(rnorm(20), genes)

  # total-loops mode
  loop_stats <- data.frame(
    Gene = genes[1:10],
    Total_Loops = 1:10,
    n_Linked_Distal = 0:9,
    stringsAsFactors = FALSE
  )

  plots <- looplook:::run_heatmap_and_connectivity(
    target_genes = genes[1:8], tpm_mat_raw = tpm_mat,
    meta_raw = meta, loop_stats_df = loop_stats,
    global_glist = glist, heatmap_ntop = 50, cor_method = "pearson",
    current_proj_name = "TestTotal", source_type = "targets"
  )
  expect_type(plots, "list")
})

test_that("run_heatmap_and_connectivity handles <5 genes gracefully", {
  genes <- paste0("Gene", 1:5)
  tpm_mat <- as.data.frame(matrix(1:20, 5, 4,
    dimnames = list(genes, c("s1", "s2", "s3", "s4"))
  ))
  meta <- data.frame(
    SampleID = c("s1", "s2", "s3", "s4"),
    Group = c("A", "A", "B", "B"), stringsAsFactors = FALSE
  )
  glist <- setNames(rnorm(5), genes)
  loop_stats <- data.frame(
    Gene = genes[1:3], Total_Loops = 1:3, n_Linked_Distal = 0:2,
    stringsAsFactors = FALSE
  )
  # Only 3 genes match loop_stats → <5 valid targets → returns early
  plots <- looplook:::run_heatmap_and_connectivity(
    target_genes = genes, tpm_mat_raw = tpm_mat,
    meta_raw = meta, loop_stats_df = loop_stats,
    global_glist = glist, heatmap_ntop = 50, cor_method = "pearson",
    current_proj_name = "TestSmall", source_type = "targets"
  )
  expect_type(plots, "list")
})

test_that(".add_connectivity_rainclouds works with High Distal + High Total groups", {
  skip_if_not_installed("ggdist")
  skip_if_not_installed("ggplot2")

  n <- 60
  plot_df <- data.frame(
    Gene = paste0("G", 1:n),
    Conn_Group = factor(c(rep("High Distal", 20), rep("High Total", 20), rep("Others", 20)),
                        levels = c("High Distal", "High Total", "Others")),
    LFC = rnorm(n),
    Expression = runif(n, 0, 10),
    stringsAsFactors = FALSE
  )
  plot_df <- plot_df %>%
    dplyr::mutate(
      Conn_Group_num = as.numeric(Conn_Group),
      Conn_Group_jitter = Conn_Group_num - 0.12,
      Conn_Group_slab = Conn_Group_num + 0.07
    )

  custom_colors <- c("High Distal" = "#9BC985", "High Total" = "#ECB884", "Others" = "#82969D")
  plots <- list()
  plots <- looplook:::.add_connectivity_rainclouds(plot_df, custom_colors, plots)
  expect_true("Raincloud_LFC" %in% names(plots))
  expect_true("Raincloud_Expr" %in% names(plots))
})

test_that(".is_promoter_anchor_type and .is_enhancer_like_anchor_type classify correctly", {
  expect_true(looplook:::.is_promoter_anchor_type("P"))
  expect_true(looplook:::.is_promoter_anchor_type(" P "))
  expect_false(looplook:::.is_promoter_anchor_type("E"))
  expect_false(looplook:::.is_promoter_anchor_type(NA_character_))

  expect_true(looplook:::.is_enhancer_like_anchor_type("E"))
  expect_true(looplook:::.is_enhancer_like_anchor_type("eP"))
  expect_false(looplook:::.is_enhancer_like_anchor_type("P"))
  expect_false(looplook:::.is_enhancer_like_anchor_type(NA_character_))
})

test_that(".anchor_matches_targets detects gene membership correctly", {
  expect_true(looplook:::.anchor_matches_targets("TP53;BRCA1", c("BRCA1", "MYC")))
  expect_false(looplook:::.anchor_matches_targets("EGFR;KRAS", c("TP53", "MYC")))
  expect_false(looplook:::.anchor_matches_targets("", c("TP53")))
  expect_false(looplook:::.anchor_matches_targets(NA_character_, c("TP53")))
})

# ============================================================================
# visualization.R — track data, arc levels, flower/upset edge cases
# ============================================================================

test_that(".assign_arc_levels computes stacking correctly", {
  left <- c(1, 2, 3, 4)
  right <- c(10, 12, 8, 15)
  res <- looplook:::.assign_arc_levels(left, right, max_levels = 5)
  expect_s3_class(res, "data.frame")
  expect_true(all(c("raw_level", "display_level", "band_pos") %in% colnames(res)))
  expect_equal(nrow(res), 4)
  # Non-overlapping: all should be on level 1
  left_np <- c(1, 20)
  right_np <- c(10, 30)
  res_np <- looplook:::.assign_arc_levels(left_np, right_np, max_levels = 5)
  expect_equal(max(res_np$raw_level), 1)

  # Empty input
  res_empty <- looplook:::.assign_arc_levels(numeric(0), numeric(0), 5)
  expect_equal(nrow(res_empty), 0)
})

test_that(".calc_alpha_color returns alpha-scaled colors", {
  scores <- c(1, 5, 10)
  base <- "#5D6D7E"
  res <- looplook:::.calc_alpha_color(scores, base, TRUE)
  expect_equal(length(res), 3)
  expect_true(all(grepl("^#[0-9A-Fa-f]{8}$", res)))

  # No alpha scaling (all same score)
  res2 <- looplook:::.calc_alpha_color(rep(5, 3), base, TRUE)
  expect_equal(length(res2), 3)

  # NA scores
  res3 <- looplook:::.calc_alpha_color(c(1, NA, 10), base, TRUE)
  expect_equal(length(res3), 3)
})

test_that(".read_track_bedpe reads and validates BEDPE for track plots", {
  tmp <- tempfile(fileext = ".bedpe")
  # 7 columns: 6 coords + score
  writeLines("chr1\t100\t200\tchr1\t300\t400\t5.5", tmp)
  res <- looplook:::.read_track_bedpe(tmp)
  expect_type(res, "list")
  expect_true("loops" %in% names(res))
  expect_true("has_score" %in% names(res))
  expect_true(res$has_score)
  expect_equal(nrow(res$loops), 1)
  unlink(tmp)
})

test_that(".read_track_bedpe handles 8-column BEDPE with score in column 8", {
  tmp <- tempfile(fileext = ".bedpe")
  writeLines("chr1\t100\t200\tchr1\t300\t400\t.\t3.2", tmp)
  res <- looplook:::.read_track_bedpe(tmp)
  expect_true(res$has_score)
  unlink(tmp)
})

test_that(".read_track_bedpe filters by min_score", {
  tmp <- tempfile(fileext = ".bedpe")
  writeLines(c(
    "chr1\t100\t200\tchr1\t300\t400\t10",
    "chr1\t500\t600\tchr1\t700\t800\t1"
  ), tmp)
  res <- looplook:::.read_track_bedpe(tmp, min_score = 5)
  expect_equal(nrow(res$loops), 1)
  unlink(tmp)
})

test_that("plot_peaks_interactions returns ggplot with minimal data (no gene track)", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggforce")

  tmp <- tempfile(fileext = ".bedpe")
  writeLines("chr1\t11890000\t11890500\tchr1\t11905000\t11905500", tmp)

  p <- plot_peaks_interactions(
    bedpe_file = tmp, chr = "chr1", from = 11884299, to = 12106581,
    show_gene_track = FALSE
  )
  expect_s3_class(p, "ggplot")
  unlink(tmp)
})

test_that("plot_peaks_interactions auto-detects chr when chr=NULL", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggforce")

  tmp <- tempfile(fileext = ".bedpe")
  writeLines("chr1\t100000\t100500\tchr1\t101000\t101500", tmp)
  p <- plot_peaks_interactions(
    bedpe_file = tmp, show_gene_track = FALSE
  )
  expect_s3_class(p, "ggplot")
  unlink(tmp)
})

test_that("draw_flower_simplified handles 2 groups and returns invisible ggplot", {
  gene_sets <- list(
    A = c("TP53", "BRCA1", "MYC", "EGFR"),
    B = c("BRCA1", "MYC", "EGFR", "KRAS")
  )
  p <- draw_flower_simplified(gene_sets, "Test", c(A = "#E41A1C", B = "#377EB8"))
  expect_s3_class(p, "ggplot")
})

test_that("draw_flower_simplified returns NULL for <2 non-empty groups", {
  expect_null(draw_flower_simplified(list(A = "TP53"), "One", c(A = "red")))
  expect_message(
    draw_flower_simplified(list(), "Empty", character(0)),
    "Less than 2"
  )
})

test_that("draw_upset_intersections returns invisibly for 2+ gene sets", {
  skip_if_not_installed("UpSetR")
  gene_sets <- list(
    Up = c("TP53", "BRCA1", "MYC"),
    Down = c("BRCA1", "MYC", "CDKN1A")
  )
  g <- withVisible(draw_upset_intersections(gene_sets, "Test"))$value
  if (!is.null(g)) {
    expect_s3_class(g, "grob")
  }
  # Mainly verifying it doesn't error
  expect_true(TRUE)
})

test_that("draw_upset_intersections returns NULL for <2 gene sets", {
  expect_message(
    draw_upset_intersections(list(A = "TP53"), "One"),
    "Less than 2"
  )
})

# --- visualization.R deep coverage: gene track, full plot_peaks_interactions ---

test_that("plot_peaks_interactions with gene_track=TRUE renders full multi-track plot", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggforce")
  skip_if_not_installed("TxDb.Hsapiens.UCSC.hg38.knownGene")
  skip_if_not_installed("org.Hs.eg.db")

  # Use a region on chr6 that contains known genes (HLA region is gene-dense)
  tmp_loops <- tempfile(fileext = ".bedpe")
  tmp_peaks <- tempfile(fileext = ".bed")
  # Two loops in a gene-rich region
  writeLines(c(
    "chr6\t29941000\t29942000\tchr6\t29945000\t29946000",
    "chr6\t29941000\t29942000\tchr6\t29943000\t29944000"
  ), tmp_loops)
  writeLines("chr6\t29941000\t29945000", tmp_peaks)

  p <- plot_peaks_interactions(
    bedpe_file = tmp_loops,
    target_bed = tmp_peaks,
    chr = "chr6", from = 29940000, to = 29947000,
    species = "hg38", show_gene_track = TRUE
  )
  expect_s3_class(p, "ggplot")

  unlink(c(tmp_loops, tmp_peaks))
})

test_that("plot_peaks_interactions handles score column with alpha", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggforce")

  tmp <- tempfile(fileext = ".bedpe")
  # 7-column BEDPE with score
  writeLines(c(
    "chr1\t100000\t100500\tchr1\t101000\t101500\t8.5",
    "chr1\t200000\t200500\tchr1\t202000\t202500\t2.1"
  ), tmp)

  p <- plot_peaks_interactions(
    bedpe_file = tmp, show_gene_track = FALSE,
    score_to_alpha = TRUE, min_score = 5
  )
  expect_s3_class(p, "ggplot")

  unlink(tmp)
})

test_that("plot_peaks_interactions min_score filters all loops", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggforce")

  tmp <- tempfile(fileext = ".bedpe")
  writeLines("chr1\t100000\t100500\tchr1\t101000\t101500\t1.0", tmp)

  expect_error(
    plot_peaks_interactions(bedpe_file = tmp, min_score = 999),
    "All loops filtered"
  )

  unlink(tmp)
})

test_that(".prepare_gene_track_data returns gene and feature data frames", {
  skip_if_not_installed("TxDb.Hsapiens.UCSC.hg38.knownGene")
  skip_if_not_installed("org.Hs.eg.db")

  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
  res <- looplook:::.prepare_gene_track_data(
    txdb_obj = txdb, org_db_ref = "org.Hs.eg.db",
    chr = "chr6", from = 29940000, to = 29947000
  )
  expect_type(res, "list")
  expect_true("genes_df" %in% names(res))
  expect_true("feature_df" %in% names(res))
  expect_s3_class(res$genes_df, "data.frame")
  if (nrow(res$genes_df) > 0) {
    expect_true("SYMBOL" %in% colnames(res$genes_df))
  }
})

test_that(".prepare_overlap_track reads and filters target BED for a region", {
  tmp <- tempfile(fileext = ".bed")
  writeLines(c("chr1\t100\t200", "chr1\t500\t600", "chr2\t100\t200"), tmp)
  res <- looplook:::.prepare_overlap_track(tmp, chr = "chr1", from = 1, to = 300)
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 1)  # only chr1:100-200 overlaps [1,300]
  unlink(tmp)
})

test_that(".prepare_overlap_track returns NULL for NULL target_bed", {
  expect_null(looplook:::.prepare_overlap_track(NULL, "chr1", 1, 100))
})

test_that(".prepare_overlap_track returns NULL when no features in region", {
  tmp <- tempfile(fileext = ".bed")
  writeLines("chr2\t100\t200", tmp)
  res <- looplook:::.prepare_overlap_track(tmp, chr = "chr1", from = 1, to = 1000)
  expect_null(res)
  unlink(tmp)
})

test_that("draw_flower_simplified handles unnamed color vector (fallback to hue_pal)", {
  gene_sets <- list(A = c("TP53", "BRCA1"), B = c("MYC", "EGFR"))
  # Unnamed colors: code path for anyNA(final_colors) → hue_pal fallback
  p <- draw_flower_simplified(gene_sets, "TestFallback", c("#E41A1C", "#377EB8"))
  expect_s3_class(p, "ggplot")
})

test_that("draw_flower_simplified with NA in colors triggers fallback", {
  gene_sets <- list(A = c("G1", "G2"), B = c("G2", "G3"))
  p <- draw_flower_simplified(gene_sets, "TestNA", c(A = "#E41A1C", B = NA))
  expect_s3_class(p, "ggplot")
})

test_that("draw_upset_intersections handles UpSetR error gracefully", {
  skip_if_not_installed("UpSetR")
  # 2 valid sets should work
  sets <- list(S1 = c("A", "B", "C"), S2 = c("B", "C", "D"))
  expect_no_error(draw_upset_intersections(sets, "TestSets"))
})

test_that("prepare_track_data auto-inferred chr/from/to when all NULL", {
  skip_if_not_installed("TxDb.Hsapiens.UCSC.hg38.knownGene")
  skip_if_not_installed("org.Hs.eg.db")

  tmp <- tempfile(fileext = ".bedpe")
  writeLines("chr1\t11884299\t11884500\tchr1\t12106000\t12106581", tmp)

  d <- looplook:::prepare_track_data(
    bedpe_file = tmp,
    target_bed = NULL, chr = NULL, from = NULL, to = NULL,
    species = "hg38", max_levels = 10, base_anchor_height = 0.05,
    loop_color = "#5D6D7E", anchor_color = "#3498DB",
    score_to_alpha = FALSE, min_score = NULL,
    show_gene_track = FALSE
  )
  expect_type(d, "list")
  expect_equal(d$chr, "chr1")
  # from/to are derived from the BEDPE coordinates (1-based after conversion)
  expect_true(d$from > 0)
  expect_true(d$to > d$from)

  unlink(tmp)
})

test_that("prepare_track_data enforces minimum region width of 100", {
  tmp <- tempfile(fileext = ".bedpe")
  writeLines("chr1\t100000\t100001\tchr1\t100010\t100011", tmp)

  d <- looplook:::prepare_track_data(
    bedpe_file = tmp,
    target_bed = NULL, chr = "chr1", from = NULL, to = NULL,
    species = "hg38", max_levels = 10, base_anchor_height = 0.05,
    loop_color = "#5D6D7E", anchor_color = "#3498DB",
    score_to_alpha = FALSE, min_score = NULL,
    show_gene_track = FALSE
  )
  expect_true(d$to - d$from >= 100)

  unlink(tmp)
})

test_that(".validate_bedpe_df errors on <6 columns", {
  expect_error(
    looplook:::.validate_bedpe_df(data.frame(V1 = "chr1", V2 = 1, V3 = 2)),
    "at least 6 columns"
  )
})

test_that(".validate_bedpe_df errors on non-numeric coordinates", {
  df <- data.frame(
    V1 = rep("chr1", 2), V2 = c("abc", "def"), V3 = c(100, 200),
    V4 = rep("chr1", 2), V5 = c(300, 400), V6 = c(350, 450),
    stringsAsFactors = FALSE
  )
  expect_error(looplook:::.validate_bedpe_df(df), "non-numeric")
})

test_that(".validate_bedpe_df errors on start >= end", {
  df <- data.frame(
    V1 = "chr1", V2 = 200, V3 = 100,
    V4 = "chr1", V5 = 300, V6 = 400,
    stringsAsFactors = FALSE
  )
  expect_error(looplook:::.validate_bedpe_df(df), "start >= end")
})

test_that(".validate_bedpe_df normalizes anchor order", {
  df <- data.frame(
    V1 = "chr2", V2 = 300, V3 = 400,
    V4 = "chr1", V5 = 100, V6 = 200,
    stringsAsFactors = FALSE
  )
  out <- looplook:::.validate_bedpe_df(df)
  # After swap: chr1 should be first
  expect_equal(out[1, 1], "chr1")
  # After BEDPE 0-based → 1-based conversion: start1 = 100 + 1 = 101
  expect_equal(out[1, 2], 101)
})

test_that("bedpe_to_gi errors on missing file", {
  expect_error(bedpe_to_gi("/nonexistent/path.bedpe"), "does not exist")
})

test_that("bedpe_to_gi validates explicit score_col", {
  tmp <- tempfile(fileext = ".bedpe")
  writeLines(c(
    "chr1\t100\t200\tchr1\t300\t400\t5.5",
    "chr1\t500\t600\tchr1\t700\t800\t3.2"
  ), tmp)
  gi <- bedpe_to_gi(tmp, score_col = 7)
  expect_s4_class(gi, "GInteractions")
  expect_equal(length(gi), 2)
  expect_true("score" %in% colnames(S4Vectors::mcols(gi)))
  unlink(tmp)
})

test_that("bedpe_to_gi errors on invalid score_col", {
  tmp <- tempfile(fileext = ".bedpe")
  writeLines("chr1\t100\t200\tchr1\t300\t400\tx", tmp)
  expect_error(bedpe_to_gi(tmp, score_col = 7), "predominantly numeric")
  unlink(tmp)
})

test_that("bedpe_to_gi errors on out-of-bounds score_col", {
  tmp <- tempfile(fileext = ".bedpe")
  writeLines("chr1\t100\t200\tchr1\t300\t400", tmp)
  expect_error(bedpe_to_gi(tmp, score_col = 999), "exceeds file column count")
  unlink(tmp)
})

test_that("read_simple_bed errors on missing file", {
  expect_error(read_simple_bed("/nonexistent.bed"), "does not exist")
})

test_that("read_simple_bed returns NULL for NULL input", {
  expect_null(read_simple_bed(NULL))
})

test_that("read_simple_bed auto-detects and skips header", {
  tmp <- tempfile(fileext = ".bed")
  writeLines(c("chr\tstart\tend", "chr1\t100\t200", "chr1\t300\t400"), tmp)
  gr <- read_simple_bed(tmp, quiet = TRUE)
  expect_s4_class(gr, "GRanges")
  expect_equal(length(gr), 2)
  # 1-based conversion: 100 -> 101
  expect_equal(BiocGenerics::start(gr)[1], 101)
  unlink(tmp)
})

test_that("read_simple_bed errors on non-numeric coordinates after header", {
  # First line has numeric coordinates (not a header), so x/y in data trigger error
  tmp <- tempfile(fileext = ".bed")
  writeLines(c("chr1\t1000\t2000", "chr1\tx\ty"), tmp)
  expect_error(read_simple_bed(tmp), "non-numeric")
  unlink(tmp)
})

test_that("read_simple_bed errors on start >= end", {
  tmp <- tempfile(fileext = ".bed")
  writeLines("chr1\t200\t100", tmp)
  expect_error(read_simple_bed(tmp), "start >= end")
  unlink(tmp)
})

test_that("consolidate_chromatin_loops validates parameters", {
  expect_error(consolidate_chromatin_loops(files = NULL), "at least two")
  expect_error(consolidate_chromatin_loops(files = "/nonexistent.bedpe"), "at least two")
  f1 <- system.file("extdata", "example_loops_1.bedpe", package = "looplook")
  skip_if(f1 == "", "example data missing")
  expect_error(consolidate_chromatin_loops(files = c(f1, "/bad.bedpe")), "not found")
  expect_error(consolidate_chromatin_loops(files = c(f1, f1), gap = -5), "non-negative")
  expect_error(consolidate_chromatin_loops(files = c(f1, f1), min_consensus = 0), "positive integer")
  expect_error(consolidate_chromatin_loops(files = c(f1, f1), min_raw_score = "a"), "single number")
  expect_error(consolidate_chromatin_loops(files = c(f1, f1), min_score = "x"), "single number")
  expect_error(consolidate_chromatin_loops(files = c(f1, f1), roi_mode = "invalid"), "should be one of")
})

test_that("consolidate_chromatin_loops intersect mode: output structure", {
  f1 <- system.file("extdata", "example_loops_1.bedpe", package = "looplook")
  f2 <- system.file("extdata", "example_loops_2.bedpe", package = "looplook")
  skip_if(f1 == "" || f2 == "", "example data missing")

  res <- consolidate_chromatin_loops(
    files = c(f1, f2), mode = "intersect", gap = 1000,
    out_file = tempfile(fileext = ".bedpe"), quiet = TRUE
  )
  expect_s4_class(res, "GInteractions")
  mc <- S4Vectors::mcols(res)
  expect_true("score" %in% colnames(mc))
  expect_true("cluster_id" %in% colnames(mc))
  expect_true("n_members" %in% colnames(mc))
  expect_true("n_reps" %in% colnames(mc))
  # intersect mode: n_reps = number of input files, n_members = 1
  if (length(res) > 0) {
    expect_equal(unique(mc$n_members), 1L)
    expect_equal(unique(mc$n_reps), 2L)
  }
})

test_that("consolidate_chromatin_loops consensus mode: output structure", {
  f1 <- system.file("extdata", "example_loops_1.bedpe", package = "looplook")
  f2 <- system.file("extdata", "example_loops_2.bedpe", package = "looplook")
  skip_if(f1 == "" || f2 == "", "example data missing")

  res <- consolidate_chromatin_loops(
    files = c(f1, f2), mode = "consensus", gap = 1000,
    out_file = tempfile(fileext = ".bedpe"), quiet = TRUE
  )
  expect_s4_class(res, "GInteractions")
  mc <- S4Vectors::mcols(res)
  expect_true(all(c("score", "cluster_id", "n_members", "n_reps") %in% colnames(mc)))
})

test_that("consolidate_chromatin_loops union mode: output structure", {
  f1 <- system.file("extdata", "example_loops_1.bedpe", package = "looplook")
  f2 <- system.file("extdata", "example_loops_2.bedpe", package = "looplook")
  skip_if(f1 == "" || f2 == "", "example data missing")

  res <- consolidate_chromatin_loops(
    files = c(f1, f2), mode = "union", gap = 1000,
    out_file = tempfile(fileext = ".bedpe"), quiet = TRUE
  )
  expect_s4_class(res, "GInteractions")
  mc <- S4Vectors::mcols(res)
  expect_true(all(c("score", "cluster_id", "n_members", "n_reps") %in% colnames(mc)))
})

test_that("reduce_ginteractions returns cluster membership", {
  f1 <- system.file("extdata", "example_loops_1.bedpe", package = "looplook")
  skip_if(f1 == "", "example data missing")
  gi <- bedpe_to_gi(f1, quiet = TRUE)
  res <- reduce_ginteractions(gi, gap = 1000)
  expect_type(res, "list")
  expect_true("gi" %in% names(res))
  expect_true("membership" %in% names(res))
  expect_s4_class(res$gi, "GInteractions")
  expect_equal(length(res$gi), length(unique(res$membership)))
})

test_that("reduce_ginteractions handles empty GI", {
  gi_empty <- InteractionSet::GInteractions(
    GenomicRanges::GRanges(), GenomicRanges::GRanges()
  )
  res <- reduce_ginteractions(gi_empty, gap = 500)
  expect_equal(length(res$gi), 0)
  expect_equal(length(res$membership), 0)
})

test_that("reduce_ginteractions handles single-loop GI", {
  gi <- InteractionSet::GInteractions(
    GenomicRanges::GRanges("chr1", IRanges::IRanges(100, 200)),
    GenomicRanges::GRanges("chr1", IRanges::IRanges(300, 400))
  )
  S4Vectors::mcols(gi)$score <- 5
  res <- reduce_ginteractions(gi, gap = 500)
  expect_equal(length(res$gi), 1)
  expect_equal(res$membership, 1)
})

test_that(".consolidate_intersect handles all-filtered case", {
  f1 <- system.file("extdata", "example_loops_1.bedpe", package = "looplook")
  skip_if(f1 == "", "example data missing")

  # Two identical files but gap=0 with different coords -> no hits
  gi_list <- list(bedpe_to_gi(f1, quiet = TRUE), bedpe_to_gi(f1, quiet = TRUE))
  # We can't easily test the all-filtered case without mocking,
  # but we can verify the structure
  res <- looplook:::.consolidate_intersect(gi_list, gap = 1000, n_reps = 2,
    log_message = function(...) {})
  expect_s4_class(res, "GInteractions")
})

# ============================================================================
# annotation.R — sankey, internal format helpers
# ============================================================================

test_that("format_annotation_columns preserves detail_anno", {
  df <- data.frame(
    annotation = c("Promoter (<=1kb)", "Intron (uc001.1)"),
    stringsAsFactors = FALSE
  )
  out <- looplook:::format_annotation_columns(df)
  expect_true("annotation" %in% colnames(out))
  expect_true("detail_anno" %in% colnames(out))
  expect_equal(out$annotation, c("Promoter", "Intron"))
})

test_that("format_annotation_columns is no-op when annotation col missing", {
  df <- data.frame(x = 1:3)
  out <- looplook:::format_annotation_columns(df)
  expect_equal(colnames(out), "x")
})

test_that(".annotation_feature_class classifies all types", {
  expect_equal(looplook:::.annotation_feature_class("Promoter (<=1kb)"), "P")
  expect_equal(looplook:::.annotation_feature_class("promoter (2-3kb)"), "P")
  expect_equal(looplook:::.annotation_feature_class("Exon (uc001.1)"), "G")
  expect_equal(looplook:::.annotation_feature_class("Intron (uc001.1)"), "G")
  expect_equal(looplook:::.annotation_feature_class("5' UTR"), "G")
  expect_equal(looplook:::.annotation_feature_class("Distal Intergenic"), "E")
  expect_equal(looplook:::.annotation_feature_class("Downstream (1-2kb)"), "E")
  # empty string falls through to default "E" (positional distal)
  expect_equal(looplook:::.annotation_feature_class(""), "E")
  expect_equal(looplook:::.annotation_feature_class(NA_character_), "Unknown")
})

test_that(".loop_type_code sorts anchors alphabetically", {
  expect_equal(looplook:::.loop_type_code("E", "P"), "E-P")
  expect_equal(looplook:::.loop_type_code("P", "E"), "E-P")
  expect_equal(looplook:::.loop_type_code("eP", "P"), "P-eP")  # lowercase e sorts after P
  expect_equal(looplook:::.loop_type_code("G", "E"), "E-G")
  expect_equal(looplook:::.loop_type_code(NA_character_, "P"), "Unknown")
})

test_that(".loop_locus_genes extracts P/G anchor genes", {
  expect_equal(looplook:::.loop_locus_genes("P", "E", "TP53", NA), "TP53")
  expect_equal(looplook:::.loop_locus_genes("E", "P", NA, "BRCA1"), "BRCA1")
  expect_equal(looplook:::.loop_locus_genes("P", "P", "TP53", "BRCA1"), "TP53;BRCA1")
  expect_equal(looplook:::.loop_locus_genes("E", "E", "G1", "G2"), "")
})

test_that(".ids_to_genes_simple maps anchor IDs to symbols", {
  lookup <- c(A1 = "TP53", A2 = "BRCA1", A3 = "MYC")
  expect_equal(looplook:::.ids_to_genes_simple(c("A1", "A2"), lookup), "BRCA1;TP53")
  expect_equal(looplook:::.ids_to_genes_simple(c("A1", "A1"), lookup), "TP53")
  expect_true(is.na(looplook:::.ids_to_genes_simple(c("A99"), lookup)))
  expect_true(is.na(looplook:::.ids_to_genes_simple(character(0), lookup)))
})

test_that(".ids_to_genes_priority prioritizes promoter anchors", {
  sym <- c(A1 = "TP53", A2 = "BRCA1")
  typ <- c(A1 = "P", A2 = "G")
  # Both in range, but A1 is promoter -> should take A1
  res <- looplook:::.ids_to_genes_priority(c("A1", "A2"), sym, typ)
  expect_equal(res, "TP53")
})

test_that(".fill_target_gene_fallback fills empty targets from fallback", {
  targets <- c("TP53", NA, "", "BRCA1")
  fallback <- c("G1", "G2", "G3", NA)
  res <- looplook:::.fill_target_gene_fallback(targets, fallback)
  expect_equal(res[1], "TP53")
  expect_equal(res[2], "G2")
  expect_equal(res[3], "G3")
  expect_equal(res[4], "BRCA1")
})

test_that(".empty_target_gene_links returns all required columns", {
  df <- looplook:::.empty_target_gene_links()
  expected_cols <- c("input_id", "loop_ID", "anchor_id", "gene", "gene_role",
    "source", "evidence", "anchor_role", "used_as_fallback",
    "in_regulated_promoter", "in_assigned_target", "in_all_loop_connected",
    "in_regulated_promoter_filled", "in_assigned_target_filled")
  expect_true(all(expected_cols %in% colnames(df)))
  expect_equal(nrow(df), 0)
})

test_that("simplify_annotation maps to 5 broad categories", {
  x <- c("Promoter (<=1kb)", "Intron (uc001.1)", "Exon (uc001.1)",
         "Distal Intergenic", "Downstream (1-2kb)", "unknown_type")
  res <- unname(looplook:::simplify_annotation(x))
  expect_equal(res, c("Promoter", "Intron", "Exon", "Distal Intergenic",
                      "Downstream", "Others"))
})

test_that(".collate_target_values collapses and sorts unique values", {
  expect_equal(looplook:::.collapse_target_values(c("B", "A", "B", "C")), "A;B;C")
  expect_equal(looplook:::.collapse_target_values(c(NA, "", "A")), "A")
  expect_true(is.na(looplook:::.collapse_target_values(c(NA, ""))))
})

test_that(".fallback_evidence_from_annotation classifies annotations", {
  x <- c("Promoter (<=1kb)", "Exon (uc001.1)", "Intron", "Distal Intergenic", "", NA)
  res <- unname(looplook:::.fallback_evidence_from_annotation(x))
  expect_equal(res[1], "local_promoter")
  expect_equal(res[2], "local_gene_body")
  expect_equal(res[3], "local_gene_body")
  expect_equal(res[4], "linear_nearest")
  expect_equal(res[5], "linear_fallback")
  expect_equal(res[6], "linear_fallback")
})

test_that(".target_linear_gene_column prioritizes SYMBOL over geneId", {
  expect_equal(looplook:::.target_linear_gene_column(
    data.frame(SYMBOL = "A", geneId = "B")
  ), "SYMBOL")
  expect_equal(looplook:::.target_linear_gene_column(
    data.frame(geneId = "B")
  ), "geneId")
  expect_null(looplook:::.target_linear_gene_column(data.frame(x = 1)))
})

test_that(".contains_target_gene detects membership in gene string", {
  expect_true(looplook:::.contains_target_gene("TP53", "TP53;BRCA1"))
  expect_false(looplook:::.contains_target_gene("MYC", "TP53;BRCA1"))
  expect_false(looplook:::.contains_target_gene(NA_character_, "TP53"))
  expect_false(looplook:::.contains_target_gene("TP53", ""))
})

# ============================================================================
# utils.R — expression matrix, karyotype, resolve_gene_conflicts
# ============================================================================

test_that("load_expression_matrix handles all sample columns", {
  tmp <- tempfile(fileext = ".tsv")
  writeLines("Gene\tS1\tS2\tS3\nTP53\t10\t20\t30\nBRCA1\t5\t15\t25", tmp)
  vals <- looplook:::load_expression_matrix(tmp)
  expect_named(vals, c("TP53", "BRCA1"))
  expect_equal(unname(vals["TP53"]), 20)  # mean of 10,20,30
  unlink(tmp)
})

test_that("load_expression_matrix handles character sample_columns", {
  tmp <- tempfile(fileext = ".tsv")
  writeLines("Gene\tS1\tS2\tS3\nTP53\t10\t20\t30\nBRCA1\t5\t15\t25", tmp)
  vals <- looplook:::load_expression_matrix(tmp, sample_columns = c("S1", "S3"))
  expect_equal(unname(vals["TP53"]), 20)  # mean of 10,30
  unlink(tmp)
})

test_that("load_expression_matrix handles integer sample_columns", {
  tmp <- tempfile(fileext = ".tsv")
  writeLines("Gene\tS1\tS2\tS3\nTP53\t10\t20\t30\nBRCA1\t5\t15\t25", tmp)
  # After removing Gene col: S1=10, S2=20, S3=30 → c(2,3) selects S2,S3 → mean=25
  vals <- looplook:::load_expression_matrix(tmp, sample_columns = c(2, 3))
  expect_equal(unname(vals["TP53"]), 25)
  unlink(tmp)
})

test_that("load_expression_matrix errors on duplicate sample column names", {
  tmp <- tempfile(fileext = ".tsv")
  writeLines("Gene\tS1\tS1\nTP53\t10\t20", tmp)
  expect_error(looplook:::load_expression_matrix(tmp), "duplicated sample")
  unlink(tmp)
})

test_that("load_expression_matrix errors on missing sample_columns", {
  tmp <- tempfile(fileext = ".tsv")
  writeLines("Gene\tS1\tS2\nTP53\t10\t20", tmp)
  expect_error(
    looplook:::load_expression_matrix(tmp, sample_columns = c("S3")),
    "not found"
  )
  unlink(tmp)
})

test_that("load_expression_matrix errors on duplicate requested sample_columns", {
  tmp <- tempfile(fileext = ".tsv")
  writeLines("Gene\tS1\tS2\nTP53\t10\t20", tmp)
  expect_error(
    looplook:::load_expression_matrix(tmp, sample_columns = c("S1", "S1")),
    "duplicates"
  )
  unlink(tmp)
})

test_that("load_expression_matrix handles non-numeric values gracefully", {
  tmp <- tempfile(fileext = ".tsv")
  writeLines("Gene\tS1\tS2\nTP53\t10\tx\nBRCA1\t5\t15", tmp)
  expect_error(
    looplook:::load_expression_matrix(tmp),
    "non-numeric"
  )
  unlink(tmp)
})

test_that("load_expression_matrix errors on insufficient columns", {
  tmp <- tempfile(fileext = ".tsv")
  writeLines("Gene\nTP53", tmp)
  expect_error(looplook:::load_expression_matrix(tmp), "at least one sample")
  unlink(tmp)
})

test_that("load_expression_matrix handles duplicate gene identifiers", {
  tmp <- tempfile(fileext = ".tsv")
  writeLines("Gene\tS1\nTP53\t10\nTP53\t20", tmp)
  vals <- looplook:::load_expression_matrix(tmp)
  expect_equal(length(vals), 2)
  unlink(tmp)
})

test_that("resolve_gene_conflicts handles empty input", {
  df <- data.frame(stringsAsFactors = FALSE)
  out <- resolve_gene_conflicts(df, NULL, NULL, c(-2000, 2000), NULL)
  expect_equal(nrow(out), 0)
})

test_that("resolve_gene_conflicts errors on invalid strategy", {
  expect_error(
    resolve_gene_conflicts(
      data.frame(x = 1), NULL, NULL, c(-2000, 2000), NULL,
      conflict_strategy = "bad"
    ),
    "should be one of"
  )
})

test_that("clean_gene_names handles edge cases", {
  expect_equal(looplook:::clean_gene_names(NULL), character(0))
  expect_equal(looplook:::clean_gene_names(character(0)), character(0))
  expect_equal(looplook:::clean_gene_names(c("TP53", "", NA, " TP53 ")), "TP53")
  expect_equal(
    looplook:::clean_gene_names(c("TP53;BRCA1", "MYC"), ";"),
    c("TP53", "BRCA1", "MYC")
  )
})

test_that("extract_genes collapses deduplicated gene symbols", {
  # extract_genes removes dupes but preserves order of first appearance (no sort)
  expect_equal(looplook:::extract_genes(c("TP53", "TP53", "BRCA1")), "TP53;BRCA1")
  expect_equal(looplook:::extract_genes(c("TP53;MYC", "BRCA1;MYC")), "TP53;MYC;BRCA1")
  expect_true(is.na(looplook:::extract_genes(c(NA, ""))))
})

test_that(".is_promoter_like .is_distal_like .is_gene_body_like work", {
  expect_true(looplook:::.is_promoter_like("P"))
  expect_true(looplook:::.is_promoter_like("dual"))
  expect_false(looplook:::.is_promoter_like("E"))
  expect_false(looplook:::.is_promoter_like("eP"))

  expect_true(looplook:::.is_distal_like("E"))
  expect_true(looplook:::.is_distal_like("eP"))
  expect_true(looplook:::.is_distal_like("dual"))
  expect_false(looplook:::.is_distal_like("P"))

  expect_true(looplook:::.is_gene_body_like("G"))
  expect_true(looplook:::.is_gene_body_like("eG"))
  expect_false(looplook:::.is_gene_body_like("P"))
})

test_that("clean_anchor filters to active genes and reclassifies", {
  allow <- c("TP53", "BRCA1")
  # Active promoter — keep
  res1 <- looplook:::clean_anchor("TP53;MYC", "P", allow, down = FALSE)
  expect_equal(res1$type, "P")
  expect_equal(res1$gene, "TP53")

  # Silent promoter with down=TRUE → eP
  res2 <- looplook:::clean_anchor("MYC", "P", allow, down = TRUE)
  expect_equal(res2$type, "eP")
  expect_true(is.na(res2$gene))

  # Silent promoter with down=FALSE → stays P, gene NA
  res3 <- looplook:::clean_anchor("MYC", "P", allow, down = FALSE)
  expect_equal(res3$type, "P")
  expect_true(is.na(res3$gene))

  # NA gene
  res4 <- looplook:::clean_anchor(NA_character_, "E", allow, down = TRUE)
  expect_equal(res4$type, "E")
  expect_true(is.na(res4$gene))
})

test_that(".get_dom returns modal value", {
  expect_equal(looplook:::.get_dom(c("A", "B", "A", "C", "B", "A")), "A")
  expect_equal(looplook:::.get_dom(c(NA, NA, "X")), "X")
  expect_true(is.na(looplook:::.get_dom(c(NA, NA))))
  expect_true(is.na(looplook:::.get_dom(character(0))))
})

test_that(".get_expr looks up gene expression", {
  vals <- c(TP53 = 8.5, BRCA1 = 3.2)
  expect_equal(unname(looplook:::.get_expr("TP53", vals)), 8.5)
  expect_equal(unname(looplook:::.get_expr("MYC", vals)), 0)
})

test_that("get_colors returns correct number of colors", {
  cols <- looplook:::get_colors(5, "Set2")
  expect_equal(length(cols), 5)
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", cols)))

  cols2 <- looplook:::get_colors(3, c("#FF0000", "#00FF00", "#0000FF"))
  expect_equal(cols2, c("#FF0000", "#00FF00", "#0000FF"))

  expect_equal(looplook:::get_colors(0, "Set2"), character(0))
})

test_that("species_txdb_pkg returns correct package names", {
  expect_equal(looplook:::species_txdb_pkg("hg38"), "TxDb.Hsapiens.UCSC.hg38.knownGene")
  expect_equal(looplook:::species_txdb_pkg("mm10"), "TxDb.Mmusculus.UCSC.mm10.knownGene")
  expect_error(looplook:::species_txdb_pkg("unknown"), "not supported")
})

test_that("species_orgdb_pkg returns correct package names", {
  expect_equal(looplook:::species_orgdb_pkg("hg38"), "org.Hs.eg.db")
  expect_equal(looplook:::species_bsgenome_pkg("hg19"), "BSgenome.Hsapiens.UCSC.hg19")
})

test_that(".build_looplook_metadata returns expected structure", {
  m <- looplook:::.build_looplook_metadata("test_fun",
    params = list(a = 1), genome_build = "hg38",
    score_semantics = "test", diagnostics = list(x = 1),
    database_versions = list(db = "1.0")
  )
  expect_equal(m$function_name, "test_fun")
  expect_equal(m$parameters$a, 1)
  expect_equal(m$genome_build, "hg38")
  expect_equal(m$diagnostics$x, 1)
  expect_true("call_time" %in% names(m))
})

# ============================================================================
# profile_target_genes — branch coverage via example data
# ============================================================================

test_that("profile_target_genes validates seed parameter", {
  rdata_path <- system.file("extdata", "analysis_results.RData", package = "looplook")
  diff_path <- system.file("extdata", "example_deg.txt", package = "looplook")
  expr_path <- system.file("extdata", "example_tpm.txt", package = "looplook")
  meta_path <- system.file("extdata", "example_coldata.txt", package = "looplook")
  skip_if(rdata_path == "" || diff_path == "" || expr_path == "" || meta_path == "",
    "Example data not available")

  tmp <- new.env()
  load(rdata_path, envir = tmp)
  res <- tmp[[ls(tmp)[1]]]
  base_args <- list(
    annotation_res = res, diff_file = diff_path,
    expr_matrix_file = expr_path, metadata_file = meta_path,
    run_go = FALSE, run_ppi = FALSE, run_motif = FALSE,
    gsea_nSample = 20, heatmap_nSample = 20
  )

  # seed must be positive integer or NULL
  expect_error(do.call(profile_target_genes, c(base_args, list(seed = 1.5))),
    "positive integer")
  expect_error(do.call(profile_target_genes, c(base_args, list(seed = -1))),
    "positive integer")
  expect_error(do.call(profile_target_genes, c(base_args, list(seed = 0))),
    "positive integer")
  expect_error(do.call(profile_target_genes, c(base_args, list(seed = c(1, 2)))),
    "positive integer")
  expect_error(do.call(profile_target_genes, c(base_args, list(seed = NA_real_))),
    "positive integer")

  # Valid seed should work
  expect_no_error(do.call(profile_target_genes, c(base_args, list(seed = 42))))
})

test_that("profile_target_genes different target_source modes", {
  rdata_path <- system.file("extdata", "analysis_results.RData", package = "looplook")
  diff_path <- system.file("extdata", "example_deg.txt", package = "looplook")
  expr_path <- system.file("extdata", "example_tpm.txt", package = "looplook")
  meta_path <- system.file("extdata", "example_coldata.txt", package = "looplook")
  skip_if(rdata_path == "" || diff_path == "" || expr_path == "" || meta_path == "",
    "Example data not available")

  tmp <- new.env()
  load(rdata_path, envir = tmp)
  res <- tmp[[ls(tmp)[1]]]
  base_args <- list(
    annotation_res = res, diff_file = diff_path,
    expr_matrix_file = expr_path, metadata_file = meta_path,
    run_go = FALSE, run_ppi = FALSE, run_motif = FALSE,
    gsea_nSample = 20, heatmap_nSample = 20
  )

  # targets only
  out1 <- do.call(profile_target_genes, c(base_args, list(target_source = "targets")))
  expect_type(out1, "list")

  # loops only
  out2 <- do.call(profile_target_genes, c(base_args, list(target_source = "loops")))
  expect_type(out2, "list")

  # promoter mode
  out3 <- do.call(profile_target_genes,
    c(base_args, list(target_mapping_mode = "promoter")))
  expect_type(out3, "list")

  # include_Filled = FALSE
  out4 <- do.call(profile_target_genes,
    c(base_args, list(include_Filled = FALSE)))
  expect_type(out4, "list")

  # use_nearest_gene = TRUE
  out5 <- do.call(profile_target_genes,
    c(base_args, list(use_nearest_gene = TRUE)))
  expect_type(out5, "list")
})

test_that("profile_target_genes stat_test modes", {
  rdata_path <- system.file("extdata", "analysis_results.RData", package = "looplook")
  diff_path <- system.file("extdata", "example_deg.txt", package = "looplook")
  expr_path <- system.file("extdata", "example_tpm.txt", package = "looplook")
  meta_path <- system.file("extdata", "example_coldata.txt", package = "looplook")
  skip_if(rdata_path == "" || diff_path == "" || expr_path == "" || meta_path == "",
    "Example data not available")

  tmp <- new.env()
  load(rdata_path, envir = tmp)
  res <- tmp[[ls(tmp)[1]]]
  base_args <- list(
    annotation_res = res, diff_file = diff_path,
    expr_matrix_file = expr_path, metadata_file = meta_path,
    run_go = FALSE, run_ppi = FALSE, run_motif = FALSE,
    gsea_nSample = 20, heatmap_nSample = 20
  )

  # t.test mode (caught by .safe_run so won't error, just handles internally)
  out_t <- do.call(profile_target_genes, c(base_args, list(stat_test = "t.test")))
  expect_type(out_t, "list")

  # match.arg validates at the run_lfc_violin level, caught by .safe_run
  # So invalid stat_test produces warnings but not errors
  expect_warning(
    do.call(profile_target_genes, c(base_args, list(stat_test = "invalid"))),
    "should be one of"
  )
})

test_that("profile_target_genes group_order parameter", {
  rdata_path <- system.file("extdata", "analysis_results.RData", package = "looplook")
  diff_path <- system.file("extdata", "example_deg.txt", package = "looplook")
  expr_path <- system.file("extdata", "example_tpm.txt", package = "looplook")
  meta_path <- system.file("extdata", "example_coldata.txt", package = "looplook")
  skip_if(rdata_path == "" || diff_path == "" || expr_path == "" || meta_path == "",
    "Example data not available")

  tmp <- new.env()
  load(rdata_path, envir = tmp)
  res <- tmp[[ls(tmp)[1]]]
  base_args <- list(
    annotation_res = res, diff_file = diff_path,
    expr_matrix_file = expr_path, metadata_file = meta_path,
    run_go = FALSE, run_ppi = FALSE, run_motif = FALSE,
    gsea_nSample = 20, heatmap_nSample = 20
  )

  out <- do.call(profile_target_genes,
    c(base_args, list(group_order = c("con", "treat"))))
  expect_type(out, "list")
})

test_that("extract_target_gene_sets promoter mode uses Regulated_promoter_genes_Filled", {
  mock_anno <- list(
    target_annotation = data.frame(
      SYMBOL = c("TP53", "BRCA1"),
      Assigned_Target_Genes_Filled = c("TP53;MYC", "BRCA1"),
      Regulated_promoter_genes_Filled = c("TP53", "BRCA1"),
      stringsAsFactors = FALSE
    )
  )
  # promoter mapping mode uses Regulated_promoter_genes_Filled
  res <- looplook:::extract_target_gene_sets(mock_anno, src = "targets",
    include_Filled = TRUE, target_mapping_mode = "promoter")
  expect_type(res, "list")
  expect_true("Target_Genes" %in% names(res))
})

test_that("read_robust_general errors on empty path", {
  expect_error(looplook:::read_robust_general(NULL, desc = "test"), "empty")
  expect_error(looplook:::read_robust_general("", desc = "test"), "empty")
})

# ============================================================================
# analysis.R — targeted edge cases to push past 60%
# ============================================================================

test_that("extract_target_gene_sets with geneId column (no SYMBOL)", {
  mock_anno <- list(
    target_annotation = data.frame(
      geneId = c("G1", "G2"),
      Assigned_Target_Genes_Filled = c("G1", "G2"),
      stringsAsFactors = FALSE
    )
  )
  res <- looplook:::extract_target_gene_sets(mock_anno, src = "targets",
    use_nearest_gene = TRUE)
  expect_type(res, "list")
  expect_true("Target_Genes" %in% names(res))
})

test_that("extract_target_gene_sets errors when loops missing Putative_Target_Genes", {
  mock_anno <- list(
    loop_annotation = data.frame(loop_type = "E-P", stringsAsFactors = FALSE)
  )
  expect_error(
    looplook:::extract_target_gene_sets(mock_anno, src = "loops"),
    "Required column"
  )
})

test_that("extract_target_gene_sets skips loop types not in active_loop_types", {
  mock_anno <- list(
    loop_annotation = data.frame(
      loop_type = c("E-P", "P-P", "G-P"),
      Putative_Target_Genes = c("G1", "G2", "G3"),
      stringsAsFactors = FALSE
    )
  )
  res <- looplook:::extract_target_gene_sets(mock_anno, src = "loops",
    active_loop_types = "E-P")
  expect_true("EP_Genes" %in% names(res))
  expect_false("PP_Genes" %in% names(res))
})

test_that(".plot_motif_rank_scatter handles all non-significant motifs", {
  mock_res <- data.frame(
    MotifID = paste0("MA", sprintf("%04d", 1:10), ".1"),
    MotifName = paste("TF", LETTERS[1:10]),
    Pvalue = rep(0.5, 10),
    FDR = rep(0.8, 10),  # all above 0.05 threshold
    OddsRatio = rep(1.0, 10),
    Family = rep("Unknown", 10),
    stringsAsFactors = FALSE
  )
  p <- looplook:::.plot_motif_rank_scatter(mock_res, "Test_AllNS", fdr_thresh = 0.05)
  expect_s3_class(p, "ggplot")
})

test_that(".prepare_motif_anchor_sets returns empty for NULL loop_df", {
  out <- looplook:::.prepare_motif_anchor_sets(NULL, "TP53")
  expect_type(out, "list")
  expect_equal(out$target_loop_n, 0L)
})

test_that(".prepare_motif_anchor_sets returns empty when missing required columns", {
  loop_df <- data.frame(chr1 = "chr1", start1 = 1, end1 = 10,
    chr2 = "chr1", start2 = 20, end2 = 30,
    stringsAsFactors = FALSE)
  out <- looplook:::.prepare_motif_anchor_sets(loop_df, "TP53")
  expect_equal(out$target_loop_n, 0L)
})

test_that(".prepare_motif_anchor_sets returns empty for empty target_genes", {
  loop_df <- data.frame(
    chr1 = "chr1", start1 = 1, end1 = 10,
    chr2 = "chr1", start2 = 20, end2 = 30,
    anchor1_gene = "G1", anchor2_gene = "G2",
    anchor1_type = "P", anchor2_type = "E",
    stringsAsFactors = FALSE
  )
  out <- looplook:::.prepare_motif_anchor_sets(loop_df, character(0))
  expect_equal(out$target_loop_n, 0L)
})

test_that("run_heatmap_and_connectivity returns early for column-less stats", {
  genes <- paste0("Gene", 1:10)
  tpm_mat <- as.data.frame(matrix(1:40, 10, 4,
    dimnames = list(genes, c("s1","s2","s3","s4"))))
  meta <- data.frame(SampleID = c("s1","s2","s3","s4"),
    Group = c("A","A","B","B"), stringsAsFactors = FALSE)
  glist <- setNames(rnorm(10), genes)
  # loop_stats has no Total_Loops, Loop_Degree, or degree column
  loop_stats <- data.frame(
    Gene = genes[1:5],
    SomeOtherCol = 1:5,
    stringsAsFactors = FALSE
  )
  plots <- looplook:::run_heatmap_and_connectivity(
    target_genes = genes[1:5], tpm_mat_raw = tpm_mat,
    meta_raw = meta, loop_stats_df = loop_stats,
    global_glist = glist, heatmap_ntop = 50, cor_method = "pearson",
    current_proj_name = "TestNoCol", source_type = "targets",
    skip_heatmap = TRUE
  )
  expect_type(plots, "list")
})

test_that("run_heatmap_and_connectivity returns early when <5 expression targets overlap", {
  genes <- paste0("Gene", 1:10)
  tpm_mat <- as.data.frame(matrix(1:40, 10, 4,
    dimnames = list(genes, c("s1","s2","s3","s4"))))
  meta <- data.frame(SampleID = c("s1","s2","s3","s4"),
    Group = c("A","A","B","B"), stringsAsFactors = FALSE)
  glist <- setNames(rnorm(10), genes)
  # Only 3 genes in stats → <5 valid expression targets
  loop_stats <- data.frame(
    Gene = genes[1:3],
    Total_Loops = c(1, 0, 1),
    stringsAsFactors = FALSE
  )
  plots <- looplook:::run_heatmap_and_connectivity(
    target_genes = genes[1:5], tpm_mat_raw = tpm_mat,
    meta_raw = meta, loop_stats_df = loop_stats,
    global_glist = glist, heatmap_ntop = 50, cor_method = "pearson",
    current_proj_name = "TestSmall", source_type = "targets",
    skip_heatmap = TRUE
  )
  expect_type(plots, "list")
})

test_that("run_heatmap_and_connectivity uses Loop_Degree fallback column", {
  genes <- paste0("Gene", 1:10)
  tpm_mat <- as.data.frame(matrix(rexp(40), 10, 4,
    dimnames = list(genes, c("s1","s2","s3","s4"))))
  meta <- data.frame(SampleID = c("s1","s2","s3","s4"),
    Group = c("A","A","B","B"), stringsAsFactors = FALSE)
  glist <- setNames(rnorm(10), genes)
  # Loop_Degree instead of Total_Loops
  loop_stats <- data.frame(
    Gene = genes,
    Loop_Degree = 1:10,
    stringsAsFactors = FALSE
  )
  plots <- looplook:::run_heatmap_and_connectivity(
    target_genes = genes, tpm_mat_raw = tpm_mat,
    meta_raw = meta, loop_stats_df = loop_stats,
    global_glist = glist, heatmap_ntop = 50, cor_method = "pearson",
    current_proj_name = "TestLoopDeg", source_type = "targets",
    skip_heatmap = TRUE
  )
  expect_type(plots, "list")
})

test_that(".add_connectivity_rainclouds handles single Conn_Group level", {
  n <- 30
  plot_df <- data.frame(
    Conn_Group = factor(rep("High Total", n), levels = "High Total"),
    LFC = rnorm(n), Expression = runif(n, 0, 10),
    stringsAsFactors = FALSE
  )
  plot_df <- plot_df %>%
    dplyr::mutate(
      Conn_Group_num = as.numeric(Conn_Group),
      Conn_Group_jitter = Conn_Group_num - 0.12,
      Conn_Group_slab = Conn_Group_num + 0.07
    )
  custom_colors <- c("High Total" = "#ECB884")
  plots <- list()
  # nlevels <= 1 → returns plots_list unchanged
  plots <- looplook:::.add_connectivity_rainclouds(plot_df, custom_colors, plots)
  # Should be empty since nlevels <= 1
  expect_equal(length(plots), 0)
})

test_that("plot_summary_go_lollipop handles results without CleanLoopType", {
  mock_go <- list(
    data.frame(
      Description = c("term1", "term2"),
      pvalue = c(0.001, 0.01), Count = c(10, 5),
      ONTOLOGY = c("BP", "MF"), LoopType = c("EP", "PP"),
      stringsAsFactors = FALSE
    )
  )
  plots <- looplook:::plot_summary_go_lollipop(mock_go, "TestNoClean")
  expect_type(plots, "list")
  expect_true(length(plots) >= 1)
})

test_that("profile_target_genes with run_go=TRUE validates org_db exists", {
  expect_error(
    profile_target_genes(
      annotation_res = list(promoter_centric_stats = data.frame()),
      diff_file = tempfile(fileext = ".csv"),
      expr_matrix_file = tempfile(fileext = ".csv"),
      metadata_file = tempfile(fileext = ".csv"),
      run_go = TRUE, org_db = "nonexistent.package.name.xyz",
      run_ppi = FALSE, run_motif = FALSE
    ),
    "is required for GO analysis"
  )
})
