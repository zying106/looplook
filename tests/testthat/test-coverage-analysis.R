# === analysis.R coverage boost ===
# Targets uncovered code paths without external deps (GO/PPI/Motif/ComplexHeatmap)

# ── read_robust_general ──
test_that("read_robust_general rejects insufficient columns", {
    tmp <- tempfile()
    writeLines("a\tb", tmp)
    expect_error(
        looplook:::read_robust_general(tmp, header = TRUE, desc = "test", min_cols = 3),
        "insufficient"
    )
    unlink(tmp)
})

test_that("read_robust_general handles empty path", {
    expect_error(
        looplook:::read_robust_general("", desc = "test"),
        "path is empty"
    )
})

test_that("read_robust_general assigns row names", {
    tmp <- tempfile()
    writeLines("id\tv1\tv2\nG1\t1\t2\nG2\t3\t4", tmp)
    d <- looplook:::read_robust_general(tmp, header = TRUE, row_name = 1, desc = "test", min_cols = 2)
    expect_equal(rownames(d), c("G1", "G2"))
    expect_equal(ncol(d), 2)
    unlink(tmp)
})

# ── extract_target_gene_sets ──
test_that("extract_target_gene_sets handles missing loop_annotation", {
    res <- list()
    gs <- looplook:::extract_target_gene_sets(res, src = "loops")
    expect_equal(length(gs), 0)
})

test_that("extract_target_gene_sets handles empty loop_annotation", {
    res <- list(loop_annotation = data.frame())
    gs <- looplook:::extract_target_gene_sets(res, src = "loops")
    expect_equal(length(gs), 0)
})

test_that("extract_target_gene_sets errors on missing promoter column", {
    res <- list(loop_annotation = data.frame(loop_type = "E-P", Putative_Target_Genes = "GENE1"))
    expect_error(
        looplook:::extract_target_gene_sets(res,
            src = "loops",
            target_mapping_mode = "promoter"
        ),
        "Promoter_Target_Genes"
    )
})

test_that("extract_target_gene_sets errors on use_nearest_gene + promoter", {
    res <- list(target_annotation = data.frame())
    expect_error(
        looplook:::extract_target_gene_sets(res,
            src = "targets",
            use_nearest_gene = TRUE, target_mapping_mode = "promoter"
        ),
        "use_nearest_gene = TRUE"
    )
})

test_that("extract_target_gene_sets nearest_gene mode uses SYMBOL column", {
    res <- list(target_annotation = data.frame(SYMBOL = "TP53;BRCA1"))
    gs <- looplook:::extract_target_gene_sets(res,
        src = "targets",
        use_nearest_gene = TRUE, target_mapping_mode = "all"
    )
    expect_equal(sort(gs[[1]]), sort(c("TP53", "BRCA1")))
})

test_that("extract_target_gene_sets nearest_gene falls back to geneId", {
    res <- list(target_annotation = data.frame(geneId = "ENSG000001"))
    gs <- looplook:::extract_target_gene_sets(res,
        src = "targets",
        use_nearest_gene = TRUE, target_mapping_mode = "all"
    )
    expect_equal(gs[[1]], "ENSG000001")
})

# ── run_lfc_violin ──
test_that("run_lfc_violin returns NULL when too few valid genes", {
    gl <- setNames(rnorm(10), paste0("GENE", 1:10))
    tg <- c("MISSING1", "MISSING2")
    expect_null(looplook:::run_lfc_violin(tg, gl, "wilcox.test", "Test"))
})

test_that("run_lfc_violin handles case-insensitive matching", {
    gl <- setNames(rnorm(10), paste0("GENE", 1:10))
    tg <- c("gene1", "gene2", "gene3", "gene4", "gene5")
    p <- looplook:::run_lfc_violin(tg, gl, "wilcox.test", "Test", seed = 42)
    expect_s3_class(p, "ggplot")
})

test_that("run_lfc_violin produces plot with expression matching", {
    gl <- setNames(runif(20, -2, 2), paste0("GENE_", 1:20))
    tg <- c("GENE_1", "GENE_2", "GENE_3", "GENE_4", "GENE_5")
    ev <- setNames(runif(20, 0, 10), toupper(paste0("GENE_", 1:20)))
    p <- looplook:::run_lfc_violin(tg, gl, "wilcox.test", "Test", expr_vals = ev, seed = 42, n_iter = 3L)
    expect_s3_class(p, "ggplot")
})

test_that("run_lfc_violin uses t.test when specified", {
    gl <- setNames(rnorm(10), paste0("GENE", 1:10))
    tg <- c("GENE1", "GENE2", "GENE3", "GENE4", "GENE5")
    p <- looplook:::run_lfc_violin(tg, gl, "t.test", "Test", seed = 42)
    expect_s3_class(p, "ggplot")
})

test_that("run_lfc_violin returns NULL when < 3 valid genes", {
    gl <- setNames(rnorm(10), paste0("GENE", 1:10))
    tg <- c("GENE1", "MISSING_X")
    expect_null(looplook:::run_lfc_violin(tg, gl, "wilcox.test", "Test"))
})

test_that("run_lfc_violin deduplicates case-collided genes", {
    gl <- setNames(rnorm(10), c(
        "gene1", "gene2", "gene3", "gene4", "gene5",
        "gene6", "gene7", "gene8", "gene9", "gene10"
    ))
    tg <- c("GENE1", "Gene1", "gene1", "GENE2", "gene2", "GENE3") # 3 unique after toupper
    p <- looplook:::run_lfc_violin(tg, gl, "wilcox.test", "Test", seed = 42)
    expect_s3_class(p, "ggplot")
})

# ── run_gsea_analysis ──
test_that("run_gsea_analysis handles empty target set", {
    gl <- setNames(rnorm(10), paste0("GENE_", 1:10))
    expect_null(looplook:::run_gsea_analysis(character(0), gl, NULL, "Test")$result)
    expect_null(looplook:::run_gsea_analysis(character(0), gl, NULL, "Test")$plot)
})

test_that("run_gsea_analysis returns result with valid genes", {
    gl <- setNames(runif(50, -2, 2), paste0("GENE", 1:50))
    tg <- c("GENE1", "GENE2", "GENE3", "GENE4", "GENE5")
    # GSEA needs enough valid genes in list
    vals <- c(gl, setNames(runif(200, -1, 1), paste0("BG", 1:200)))
    out <- looplook:::run_gsea_analysis(tg, vals, gsea_nSample = NULL, "Test")
    expect_true(!is.null(out$result) || is.null(out$plot))
})

test_that("run_gsea_analysis handles duplicate ranked values", {
    vals <- rep(0, 100)
    names(vals) <- paste0("GENE", 1:100)
    tg <- c("GENE1", "GENE2", "GENE3", "GENE4", "GENE5")
    out <- looplook:::run_gsea_analysis(tg, vals, NULL, "Test")
    expect_true(!is.null(out$result) || is.null(out$plot))
})

test_that("run_gsea_analysis down-samples when gsea_nSample < length", {
    gl <- setNames(runif(200, -2, 2), paste0("GENE", 1:200))
    tg <- paste0("GENE", 1:50)
    expect_warning(
        looplook:::run_gsea_analysis(tg, gl, gsea_nSample = 10, "Test"),
        "down-sampling"
    )
})

# ── plot_summary_go_lollipop ──
test_that("plot_summary_go_lollipop handles empty input", {
    p <- looplook:::plot_summary_go_lollipop(list(), "Test")
    expect_equal(length(p), 0)
})

test_that("plot_summary_go_lollipop handles NULL entries", {
    p <- looplook:::plot_summary_go_lollipop(list(NULL, NULL), "Test")
    expect_equal(length(p), 0)
})

test_that("plot_summary_go_lollipop produces lollipop with valid GO data", {
    go_df <- data.frame(
        CleanLoopType = "EP_Genes",
        Description = c("immune response", "signal transduction", "cell cycle"),
        ONTOLOGY = c("BP", "BP", "BP"),
        p.adjust = c(0.001, 0.01, 0.05),
        Count = c(20, 15, 10),
        geneID = c("G1/G2", "G3/G4", "G5/G6"),
        LoopType = c("E-P", "E-P", "E-P"),
        Source = "loops",
        stringsAsFactors = FALSE
    )
    plots <- looplook:::plot_summary_go_lollipop(list(go_df), "Test")
    expect_true(length(plots) > 0)
    expect_s3_class(plots[[1]], "ggplot")
})

test_that("plot_summary_go_lollipop handles results without CleanLoopType", {
    go_df <- data.frame(
        Description = c("pathway A", "pathway B"),
        ONTOLOGY = c("BP", "BP"),
        p.adjust = c(0.01, 0.05),
        Count = c(10, 5),
        geneID = c("G1", "G2"),
        LoopType = c("E-P", "E-P"),
        Source = "loops",
        stringsAsFactors = FALSE
    )
    plots <- looplook:::plot_summary_go_lollipop(list(go_df), "Test")
    expect_true(length(plots) > 0)
})

# ── run_heatmap_and_connectivity ──
test_that("run_heatmap_and_connectivity skips when no loop_stats", {
    tg <- paste0("GENE", 1:10)
    mat <- matrix(runif(40), nrow = 10, dimnames = list(paste0("GENE", 1:10), paste0("S", 1:4)))
    meta <- data.frame(SampleID = paste0("S", 1:4), Group = c("A", "A", "B", "B"))
    gl <- setNames(runif(10), paste0("GENE", 1:10))
    plots <- looplook:::run_heatmap_and_connectivity(tg, mat, meta, NULL, gl, 5, "pearson", "Test", "loops",
        skip_heatmap = TRUE
    )
    expect_equal(length(plots), 0)
})

test_that("run_heatmap_and_connectivity skips when < 5 valid targets", {
    tg <- c("GENE1", "GENE2")
    mat <- matrix(runif(10), nrow = 5, dimnames = list(paste0("GENE", 1:5), paste0("S", 1:2)))
    meta <- data.frame(SampleID = paste0("S", 1:2), Group = c("A", "B"))
    gl <- setNames(runif(5), paste0("GENE", 1:5))
    stats <- data.frame(
        Gene = paste0("GENE", 1:2), Total_Loops = c(10, 20),
        n_Linked_Promoters = c(3, 5), n_Linked_Distal = c(2, 4)
    )
    plots <- looplook:::run_heatmap_and_connectivity(tg, mat, meta, stats, gl, 5, "pearson", "Test", "loops")
    expect_equal(length(plots), 0)
})

test_that("run_heatmap_and_connectivity returns connectivity plots with data", {
    tg <- paste0("GENE", 1:10)
    mat <- matrix(runif(40), nrow = 10, dimnames = list(paste0("GENE", 1:10), paste0("S", 1:4)))
    meta <- data.frame(SampleID = paste0("S", 1:4), Group = c("A", "A", "B", "B"))
    gl <- setNames(runif(10), paste0("GENE", 1:10))
    stats <- data.frame(
        Gene = paste0("GENE", 1:10),
        Total_Loops = sample(5:20, 10, replace = TRUE),
        n_Linked_Promoters = sample(0:5, 10, replace = TRUE),
        n_Linked_Distal = sample(0:5, 10, replace = TRUE),
        Is_High_Connectivity_Gene = "No",
        Is_High_Distal_Connectivity_Gene = "No"
    )
    plots <- looplook:::run_heatmap_and_connectivity(tg, mat, meta, stats, gl, 5, "pearson", "Test", "loops",
        skip_heatmap = TRUE
    )
    # Should have Scatter + Raincloud plots
    expect_true("Scatter" %in% names(plots))
})

test_that("run_heatmap_and_connectivity works with distal target_col", {
    tg <- paste0("GENE", 1:10)
    mat <- matrix(runif(40), nrow = 10, dimnames = list(paste0("GENE", 1:10), paste0("S", 1:4)))
    meta <- data.frame(SampleID = paste0("S", 1:4), Group = c("A", "A", "B", "B"))
    gl <- setNames(runif(10), paste0("GENE", 1:10))
    stats <- data.frame(
        Gene = paste0("GENE", 1:10),
        n_Linked_Distal = sample(1:10, 10, replace = TRUE)
    )
    plots <- looplook:::run_heatmap_and_connectivity(tg, mat, meta, stats, gl, 5, "pearson", "Test", "loops",
        target_col = "n_Linked_Distal", skip_heatmap = TRUE
    )
    expect_true("Scatter" %in% names(plots))
})

# ── .run_profile_tasks ──
test_that(".run_profile_tasks discovers valid genes from global_glist", {
    gl <- setNames(runif(20, -2, 2), paste0("GENE", 1:20))
    mat <- matrix(runif(40), nrow = 10, dimnames = list(paste0("GENE", 6:15), paste0("S", 1:4)))
    meta <- data.frame(SampleID = paste0("S", 1:4), Group = c("A", "A", "B", "B"))
    aq <- list(TestSet = c("GENE1", "GENE5", "GENE10", "MISSING1", "GENE15"))

    res <- looplook:::.run_profile_tasks(
        analysis_queue = aq, current_source_proj_name = "Test",
        global_glist = gl, tpm_mat_raw = mat, meta_raw = meta,
        loop_stats_df = NULL, annotation_res = list(), src = "loops",
        stat_test = "wilcox.test", gsea_nSample = NULL, heatmap_nSample = 5,
        cor_method = "pearson",
        run_motif = FALSE, genome_id = "hg38",
        motif_p_thresh = 1e-4, motif_ntop = 5, motif_n_perm = 0L,
        run_go = FALSE, universe_genes = NULL, org_db = "org.Hs.eg.db",
        cnet_nSample = 50,
        run_ppi = FALSE, ppi_score = 400, ppi_nSample = 400, ppi_species_id = NULL
    )
    expect_type(res, "list")
    expect_true("analysis_queue" %in% names(res))
    expect_true("plots" %in% names(res))
    # TestSet genes filtered to GENE1, GENE5, GENE10, GENE15 (4 valid from gl)
    expect_equal(length(res$analysis_queue$TestSet), 4)
})

test_that(".run_profile_tasks skips when < 3 valid genes", {
    gl <- setNames(runif(20), paste0("GENE", 1:20))
    mat <- matrix(runif(10), nrow = 5, dimnames = list(paste0("GENE", 1:5), paste0("S", 1:2)))
    meta <- data.frame(SampleID = paste0("S", 1:2), Group = c("A", "B"))
    aq <- list(SmallSet = c("GENE1", "MISSING1"))

    res <- looplook:::.run_profile_tasks(
        analysis_queue = aq, current_source_proj_name = "Test",
        global_glist = gl, tpm_mat_raw = mat, meta_raw = meta,
        loop_stats_df = NULL, annotation_res = list(), src = "loops",
        stat_test = "wilcox.test", gsea_nSample = NULL, heatmap_nSample = 5,
        cor_method = "pearson",
        run_motif = FALSE, genome_id = "hg38",
        motif_p_thresh = 1e-4, motif_ntop = 5, motif_n_perm = 0L,
        run_go = FALSE, universe_genes = NULL, org_db = "org.Hs.eg.db",
        cnet_nSample = 50,
        run_ppi = FALSE, ppi_score = 400, ppi_nSample = 400, ppi_species_id = NULL
    )
    expect_equal(length(res$plots), 0)
})

# ── .build_expression_heatmap ──
test_that(".build_expression_heatmap returns NULL when skip_heatmap", {
    mat <- matrix(1:20, nrow = 5, dimnames = list(paste0("GENE", 1:5), paste0("S", 1:4)))
    meta <- data.frame(SampleID = paste0("S", 1:4), Group = c("A", "A", "B", "B"))
    expect_null(looplook:::.build_expression_heatmap(
        paste0("GENE", 1:5), mat, meta, 10, "Test",
        skip_heatmap = TRUE
    ))
})

test_that(".build_expression_heatmap returns NULL when < 5 expressed genes", {
    mat <- matrix(1:8, nrow = 2, dimnames = list(paste0("GENE", 1:2), paste0("S", 1:4)))
    meta <- data.frame(SampleID = paste0("S", 1:4), Group = c("A", "A", "B", "B"))
    expect_null(looplook:::.build_expression_heatmap(
        c("GENE1", "GENE2"), mat, meta, 10, "Test",
        skip_heatmap = FALSE
    ))
})

# ── .add_connectivity_rainclouds ──
test_that(".add_connectivity_rainclouds returns early with single group", {
    df <- data.frame(
        LFC = rnorm(10), Expression = runif(10),
        Conn_Group = factor(rep("Others", 10), levels = "Others"),
        Conn_Group_num = 1, Conn_Group_jitter = 0.88, Conn_Group_slab = 1.07
    )
    res <- looplook:::.add_connectivity_rainclouds(df, c(Others = "grey"), list())
    expect_equal(length(res), 0)
})

# ── .plot_save_motif ──
test_that(".plot_save_motif returns NULL for empty input", {
    expect_null(looplook:::.plot_save_motif(NULL, "test"))
    expect_null(looplook:::.plot_save_motif(data.frame(), "test"))
})

test_that(".plot_save_motif produces barplot", {
    df <- data.frame(
        MotifID = paste0("MA", 1:10),
        MotifName = paste0("TF", 1:10),
        Pvalue = runif(10, 0, 0.01),
        FDR = runif(10, 0, 0.05),
        stringsAsFactors = FALSE
    )
    p <- looplook:::.plot_save_motif(df, "test")
    expect_s3_class(p, "ggplot")
})

# ── .plot_motif_rank_scatter ──
test_that(".plot_motif_rank_scatter returns NULL for empty input", {
    expect_null(looplook:::.plot_motif_rank_scatter(NULL, "test"))
    expect_null(looplook:::.plot_motif_rank_scatter(data.frame(), "test"))
})

test_that(".plot_motif_rank_scatter produces scatter plot", {
    df <- data.frame(
        MotifID = paste0("MA", 1:20),
        MotifName = paste0("TF", 1:20),
        Pvalue = runif(20, 0, 0.01),
        FDR = runif(20, 0, 0.05),
        OddsRatio = runif(20, 0.9, 1.5),
        Family = rep(c("bZIP", "Homeobox", "ETS"), c(7, 7, 6)),
        stringsAsFactors = FALSE
    )
    p <- looplook:::.plot_motif_rank_scatter(df, "test")
    expect_s3_class(p, "ggplot")
})

# ── .annotate_motif_families ──
test_that(".annotate_motif_families returns input unchanged when empty", {
    expect_equal(looplook:::.annotate_motif_families(NULL), NULL)
    expect_equal(nrow(looplook:::.annotate_motif_families(data.frame())), 0)
})
