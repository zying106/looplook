#' @title Integrative functional annotation and profiling of target genes
#'
#' @description
#' Integrates 3D genomic interaction data (e.g., Hi-C, HiChIP) with transcriptomic
#' profiles (RNA-seq) to evaluate the functional and regulatory landscape of target genes.
#' The workflow sequentially performs differential expression profiling, gene set enrichment,
#' transcription factor motif scanning, Gene Ontology (GO) enrichment, and
#' protein-protein interaction (PPI) network construction.
#'
#' @details
#' GSEA operates on the full target gene set by default (\code{gsea_nSample = NULL}).
#' Down-sampling is available for very large sets but introduces Monte Carlo
#' variance. Motif background anchor sampling is GC-matched, limited to 2 000
#' background regions per contrast. GSEA tie-breaking for duplicate ranked
#' values is deterministic (position-based offset). For fully reproducible
#' results, set the \code{seed}
#' parameter.
#'
#' \strong{Exploratory modules:} The GO enrichment (\code{run_go}), motif scanning
#' (\code{run_motif}), and PPI network (\code{run_ppi}) modules are
#' \emph{research-grade} analyses that depend on external databases and algorithms.
#' Results should be treated as hypothesis-generating and validated with
#' independent experimental approaches. All three modules are disabled by default.
#' \strong{Motif enrichment note:} By default, Fisher's exact test treats each
#' anchor as an independent observation. When \code{motif_n_perm > 0},
#' empirical P-values are estimated by shuffling foreground/background labels
#' within loop components, which partially accounts for component dependence.
#' Pure-label components are non-exchangeable; the procedure remains exploratory
#' rather than a fully component-level inferential model.
#'
#' @param annotation_res List. The result object returned by \code{\link{annotate_peaks_and_loops}}.
#' @param diff_file Character. Path to the differential expression file (CSV/TSV).
#' @param lfc_col Character. The column name in \code{diff_file} representing Log2 Fold Change.
#'   Default \code{"log2FoldChange"}.
#' @param expr_matrix_file Character. Path to the normalized expression matrix.
#' @param metadata_file Character. Path to the sample metadata file.
#' @param target_source Character vector. Source of target genes to analyze.
#'   Default \code{c("loops", "targets")}.
#' @param target_mapping_mode Character. Mapping strategy: \code{"all"}
#'   (all anchor-gene connections, \code{Putative_Target_Genes}) or
#'   \code{"promoter"} (promoter-only subset, \code{Promoter_Target_Genes}).
#'   Expression filtering is determined by the refinement stage of the
#'   supplied annotation object, not by this parameter.
#'   Default \code{"all"}.
#' @param loop_types Character vector. The specific loop types to analyze.
#'   Default \code{c("E-P", "P-P")}.
#' @param include_Filled Logical. If \code{TRUE}, utilizes the comprehensively merged gene assignment.
#'   Default \code{TRUE}.
#' @param use_nearest_gene Logical. If \code{TRUE}, bypasses 3D loop-based gene assignment.
#'   Default \code{FALSE}.
#' @param group_order Character vector. Optional factor levels to sort sample groups.
#'   Default \code{NULL}.
#' @param project_name Character. Prefix for all output files and plot titles.
#'   Default \code{"Analysis"}.
#' @param org_db Character. Organism annotation database (e.g., "org.Hs.eg.db").
#'   Default \code{"org.Hs.eg.db"}.
#' @param run_motif Logical. Whether to perform Transcription Factor Binding Site motif analysis.
#'   Default \code{FALSE}.
#' @param genome_id Character. Reference genome assembly for motif sequence extraction.
#'   Default \code{"hg38"}.
#' @param motif_p_thresh Numeric. P-value threshold for scanning.
#'   Default \code{1e-4}.
#' @param motif_ntop Numeric. Number of top enriched motifs to output.
#'   Default \code{5}.
#' @param motif_n_perm Integer. Number of within-component foreground/background
#'   label permutations for empirical motif P-value estimation. When positive,
#'   labels are shuffled within each loop component to partially account for
#'   anchor non-independence. \code{0} (default) retains Fisher exact test.
#'   Use \code{10-100} for testing. For inference, choose according to the
#'   number of motifs and desired P-value resolution; \code{1000} may be
#'   insufficient after multiple-testing correction across a full motif library.
#'   This is an exploratory calibration, not a fully component-level model.
#' @param run_go Logical. Whether to perform Gene Ontology (GO) enrichment.
#'   Default \code{FALSE}.
#' @param universe_genes Named numeric vector or \code{NULL}. GO background
#'   universe.  Names are gene symbols, values are ranking metrics (e.g. LFC).
#'   When \code{NULL}, the full differential table is used as background.
#'   Default \code{NULL}.
#' @param run_ppi Logical. Whether to construct Protein-Protein Interaction networks.
#'   Default \code{FALSE}.
#' @param ppi_score Numeric. Minimum combined confidence score for STRING edges.
#'   Default \code{400}.
#' @param ppi_nSample Numeric. Maximum number of genes to include in PPI.
#'   Default \code{400}.
#' @param ppi_species_id Integer or \code{NULL}. NCBI taxonomy ID for STRING
#'   PPI analysis (e.g.\ 9606 for human, 10090 for mouse). \code{NULL} (default)
#'   attempts automatic inference from the OrgDb package name. Set explicitly
#'   for species whose OrgDb naming is not recognised. See
#'   \url{https://string-db.org} for a complete taxonomy list.
#' @param heatmap_nSample Numeric. Maximum number of genes to plot in heatmap.
#'   Default \code{99999} (effectively no limit). Reduce to \code{20-50} for
#'   readable heatmaps.
#' @param gsea_nSample Numeric or \code{NULL}. Maximum number of target genes
#'   to sample for GSEA. Default \code{NULL} (no down-sampling; all target
#'   genes are used). Down-sampling introduces Monte Carlo variance and
#'   should only be used when the target set is extremely large (>1000 genes)
#'   for computational efficiency. For reproducible results, set the
#'   \code{seed} parameter when using down-sampling.
#' @param cnet_nSample Numeric. Number of top GO terms to display in cnetplot.
#'   Default \code{50}.
#' @param stat_test Character. Statistical test for LFC comparisons.
#'   Default \code{"wilcox.test"}.
#' @param cor_method Character. Method for sample correlation matrices.
#'   Default \code{"pearson"}.
#' @param seed Integer or NULL. Random seed for reproducible GSEA down-sampling
#'   and motif GC-matched background sampling. When \code{NULL} (default),
#'   the global RNG state is used; set to a positive integer for fully
#'   reproducible results. The seed is recorded in the result object as
#'   \code{attr(result, "seed")}.
#'
#' @return An invisible nested list indexed by \code{target_source} (e.g., \code{"targets"}, \code{"loops"}).
#'   Each element contains:
#'   \describe{
#'     \item{\code{go_results}}{Named list of data frames (one per gene set) containing GO enrichment results (if \code{run_go = TRUE}).}
#'     \item{\code{motif_results}}{Named list of motif enrichment result tables, indexed by analysis task. Each task contains \code{proximal} and \code{distal} data frames when available (if \code{run_motif = TRUE}).}
#'     \item{\code{target_gene_sets}}{Named list of character vectors containing target gene symbols.}
#'     \item{\code{plots}}{Named list of ggplot objects (LFC_Violin, GSEA, Heatmap, Scatter, GO_Network, PPI_Network, etc.).}
#'     \item{\code{warnings}}{Character vector of module-level warnings (e.g., "[GO] failed: ..."). Empty if all modules succeeded.}
#'   }
#'
#' @note Downstream modules are run in a fail-soft mode. Module-level failures
#'   (e.g. due to missing optional packages or network timeouts) are captured as
#'   warnings and stored in the returned \code{warnings} element, allowing other
#'   modules to complete. To proactively disable specific modules, set
#'   \code{run_go = FALSE}, \code{run_ppi = FALSE}, or \code{run_motif = FALSE}.
#'
#' @examples
#' rdata_path <- system.file("extdata", "analysis_results.RData", package = "looplook")
#' diff_path <- system.file("extdata", "example_deg.txt", package = "looplook")
#' expr_path <- system.file("extdata", "example_tpm.txt", package = "looplook")
#' meta_path <- system.file("extdata", "example_coldata.txt", package = "looplook")
#' tmp <- new.env()
#' load(rdata_path, envir = tmp)
#' res <- tmp[[ls(tmp)[1]]]
#' profile_res <- profile_target_genes(
#'     annotation_res = res,
#'     diff_file = diff_path,
#'     expr_matrix_file = expr_path,
#'     metadata_file = meta_path,
#'     run_go = FALSE,
#'     run_ppi = FALSE,
#'     run_motif = FALSE,
#'     heatmap_nSample = 20,
#'     gsea_nSample = 20,
#'     cnet_nSample = 5
#' )
#' names(profile_res)
#'
#' @seealso \code{\link{annotate_peaks_and_loops}} for initial 3D annotation,
#'   \code{\link{refine_loop_anchors_by_expression}} for expression-aware refinement.
#'
#' @export
profile_target_genes <- function(
  annotation_res,
  diff_file,
  lfc_col = "log2FoldChange",
  expr_matrix_file,
  metadata_file,
  target_source = c("loops", "targets"),
  target_mapping_mode = c("all", "promoter"),
  loop_types = c("E-P", "P-P"),
  include_Filled = TRUE,
  use_nearest_gene = FALSE,
  group_order = NULL,
  project_name = "Analysis",
  org_db = "org.Hs.eg.db",
  run_motif = FALSE,
  genome_id = "hg38",
  motif_p_thresh = 1e-4,
  motif_ntop = 5,
  motif_n_perm = 0L,
  run_go = FALSE,
  run_ppi = FALSE,
  ppi_score = 400,
  ppi_nSample = 400,
  ppi_species_id = NULL,
  heatmap_nSample = 99999,
  gsea_nSample = NULL,
  cnet_nSample = 50,
  universe_genes = NULL,
  stat_test = "wilcox.test",
  cor_method = "pearson",
  seed = NULL
) {
    target_source <- match.arg(target_source, several.ok = TRUE)
    target_mapping_mode <- match.arg(target_mapping_mode)
    genome_id <- match.arg(genome_id, c("hg38", "hg19", "mm10", "mm9"))

    if (!is.numeric(motif_n_perm) || length(motif_n_perm) != 1L ||
        is.na(motif_n_perm) || !is.finite(motif_n_perm) ||
        motif_n_perm < 0 || motif_n_perm != floor(motif_n_perm)) {
        stop("`motif_n_perm` must be a single finite non-negative integer.",
            call. = FALSE
        )
    }

    # Seed management: withr::local_seed provides a local RNG context
    # without leaking .Random.seed into the global environment.
    if (!is.null(seed)) {
        if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) ||
            seed != as.integer(seed) || seed < 1L) {
            stop("`seed` must be a single positive integer or NULL.", call. = FALSE)
        }
        withr::local_seed(seed)
    }
    used_seed <- if (!is.null(seed)) seed else NULL

    root_project_name <- project_name
    if (target_mapping_mode == "promoter") root_project_name <- paste0(root_project_name, "_Promoter")
    if (use_nearest_gene && "targets" %in% target_source) {
        root_project_name <- paste0(root_project_name, "_RefNearest")
    } else if ("targets" %in% target_source) {
        if (!include_Filled) root_project_name <- paste0(root_project_name, "_LoopOnly")
    }

    message(">>> Analysis Init | Root Project: ", root_project_name)

    if (run_go) {
        if (!requireNamespace(org_db, quietly = TRUE)) {
            stop("Package '", org_db, "' is required for GO analysis. Please install it.")
        }
    }

    if (run_motif) {
        bs_pkg <- species_bsgenome_pkg(genome_id)
        if (is.null(bs_pkg) || !requireNamespace(bs_pkg, quietly = TRUE)) {
            warning("Unsupported genome or missing BSgenome package. Disabling Motif Analysis.")
            run_motif <- FALSE
        }
    }

    if (run_ppi) {
        if (!requireNamespace("STRINGdb", quietly = TRUE)) {
            stop("Package 'STRINGdb' is required for PPI analysis. Please install it.")
        }
        if (!requireNamespace("ggraph", quietly = TRUE)) {
            stop("Package 'ggraph' is required for PPI analysis. Please install it.")
        }
    }

    message("--- Reading files...")
    diff_df_raw <- read_robust_general(diff_file, header = TRUE, row_name = 1, desc = "Diff", min_cols = 1)
    tpm_mat_raw <- read_robust_general(expr_matrix_file, header = TRUE, row_name = 1, desc = "Expr", min_cols = 1)

    # Validate expression matrix -- same checks as load_expression_matrix()
    gene_ids <- rownames(tpm_mat_raw)
    if (is.null(gene_ids) || length(gene_ids) == 0) {
        stop("Expression matrix has no row names (gene identifiers).")
    }
    dup_ids <- unique(gene_ids[duplicated(gene_ids)])
    if (length(dup_ids) > 0) {
        warning(sprintf(
            "Expression matrix contains %d duplicated gene identifier(s). The first occurrence is retained; subsequent duplicates are removed. Example(s): %s",
            length(dup_ids), paste(head(dup_ids, 3), collapse = ", ")
        ))
        tpm_mat_raw <- tpm_mat_raw[!duplicated(gene_ids), , drop = FALSE]
        gene_ids <- rownames(tpm_mat_raw)
    }
    # Check for non-numeric content
    numeric_cols <- vapply(tpm_mat_raw, is.numeric, logical(1))
    if (!all(numeric_cols)) {
        stop(
            "Expression matrix contains non-numeric column(s): ",
            paste(colnames(tpm_mat_raw)[!numeric_cols], collapse = ", ")
        )
    }
    # Warn about case collisions
    id_upper <- toupper(gene_ids)
    case_dup <- unique(id_upper[duplicated(id_upper)])
    if (length(case_dup) > 0) {
        warning(sprintf(
            "Expression matrix contains %d gene identifier(s) that collide after case-insensitive matching. Example(s): %s",
            length(case_dup), paste(head(case_dup, 3), collapse = ", ")
        ))
    }

    # Validate universe_genes if provided
    if (!is.null(universe_genes)) {
        if (!is.numeric(universe_genes) || is.null(names(universe_genes))) {
            stop("`universe_genes` must be a named numeric vector: ",
                "names = gene IDs, values = ranking metric (e.g. LFC).",
                call. = FALSE
            )
        }
        keep <- !is.na(names(universe_genes)) &
            nzchar(trimws(names(universe_genes))) &
            is.finite(universe_genes)
        universe_genes <- universe_genes[keep]
        names(universe_genes) <- trimws(names(universe_genes))
        if (anyDuplicated(toupper(names(universe_genes)))) {
            stop("`universe_genes` contains duplicate or case-colliding gene IDs.",
                call. = FALSE
            )
        }
        if (length(universe_genes) == 0L) {
            stop("`universe_genes` contains no valid entries after validation ",
                "(all were non-finite, NA-named, or blank).",
                call. = FALSE
            )
        }
    }

    expr_vals <- rowMeans(tpm_mat_raw, na.rm = TRUE)
    names(expr_vals) <- toupper(rownames(tpm_mat_raw))
    meta_raw <- read_robust_general(metadata_file, header = TRUE, row_name = NULL, desc = "Meta", min_cols = 2)
    colnames(meta_raw)[c(1, 2)] <- c("SampleID", "Group")
    meta_raw$SampleID <- trimws(as.character(meta_raw$SampleID))
    if (!is.null(group_order)) meta_raw$Group <- factor(meta_raw$Group, levels = group_order)
    meta_raw <- meta_raw %>% dplyr::arrange(Group)

    if (!lfc_col %in% colnames(diff_df_raw)) stop("LFC column ", lfc_col, " not found")
    clean_diff <- diff_df_raw[!is.na(diff_df_raw[[lfc_col]]) & is.finite(diff_df_raw[[lfc_col]]), , drop = FALSE]
    global_glist <- sort(setNames(clean_diff[[lfc_col]], rownames(clean_diff)), decreasing = TRUE)

    # Warn if the LFC-ranked list (used as GO universe / GSEA background)
    # appears to be a DEG-only subset rather than a full tested-gene table.
    # A small universe relative to the expression matrix biases enrichment
    # statistics toward the pre-selected genes.
    n_expr_genes <- length(expr_vals)
    n_universe <- length(global_glist)
    if (n_universe > 0 && n_expr_genes > 0 &&
        n_universe < 500 && n_universe / n_expr_genes < 0.10) {
        warning(
            "The LFC file contains only ", n_universe, " genes with finite LFC, ",
            "compared to ", n_expr_genes, " genes in the expression matrix (",
            round(100 * n_universe / n_expr_genes, 1), "%). ",
            "If this is a DEG-only file (not a complete DESeq2/edgeR tested-gene ",
            "table), GO enrichment and GSEA backgrounds will be biased toward the ",
            "pre-selected gene set.  For proper enrichment statistics, provide the ",
            "full tested-gene table with all genes (including non-significant ones).",
            call. = FALSE
        )
    }

    loop_stats_df <- annotation_res$promoter_centric_stats
    final_master_list <- list()

    for (src in target_source) {
        current_source_proj_name <- paste0(root_project_name, "_", src)
        message("\n================================================================")
        message(">>> Processing Source: [", src, "]")

        active_loop_types <- if (src == "loops") loop_types else NULL
        raw_gene_sets <- extract_target_gene_sets(annotation_res, src, active_loop_types, include_Filled, use_nearest_gene, target_mapping_mode)

        if (length(raw_gene_sets) == 0) {
            warning("No gene sets found. Skipping.")
            next
        }

        analysis_queue <- raw_gene_sets
        task_results <- .run_profile_tasks(
            analysis_queue = analysis_queue,
            current_source_proj_name = current_source_proj_name,
            global_glist = global_glist,
            tpm_mat_raw = tpm_mat_raw,
            meta_raw = meta_raw,
            loop_stats_df = loop_stats_df,
            annotation_res = annotation_res,
            src = src,
            expr_vals = expr_vals,
            stat_test = stat_test,
            gsea_nSample = gsea_nSample,
            heatmap_nSample = heatmap_nSample,
            cor_method = cor_method,
            run_motif = run_motif,
            genome_id = genome_id,
            motif_p_thresh = motif_p_thresh,
            motif_ntop = motif_ntop,
            motif_n_perm = motif_n_perm,
            run_go = run_go,
            universe_genes = universe_genes,
            org_db = org_db,
            cnet_nSample = cnet_nSample,
            run_ppi = run_ppi,
            ppi_score = ppi_score,
            ppi_nSample = ppi_nSample,
            ppi_species_id = ppi_species_id
        )
        analysis_queue <- task_results$analysis_queue
        source_go_results <- task_results$go_results
        source_motif_results <- task_results$motif_results
        source_plots <- task_results$plots

        if (run_go && length(source_go_results) > 0) {
            p_go_sum <- plot_summary_go_lollipop(source_go_results, current_source_proj_name)
            if (length(p_go_sum) > 0) source_plots$Summary_GO <- p_go_sum
        }
        final_master_list[[src]] <- list(
            go_results = source_go_results,
            motif_results = source_motif_results,
            target_gene_sets = analysis_queue, plots = source_plots,
            warnings = task_results$warnings
        )
    }
    message("\n All analysis complete.")
    attr(final_master_list, "seed") <- used_seed
    return(invisible(final_master_list))
}

#' @title Robust Data Reader
#' @description Safely reads standard genomic formats utilizing `data.table::fread` for intelligent format inference.
#' @param f Character. Path to the input file.
#' @param header Logical. Whether the file contains a header row.
#' @param row_name Integer or NULL. Column index to be used as row names.
#' @param desc Character. Short description for error logging.
#' @param min_cols Integer. Minimum number of columns required.
#' @return A data frame.
#' @keywords internal
#' @noRd
read_robust_general <- function(f, header = FALSE, row_name = NULL, desc = "file", min_cols = 3) {
    if (is.null(f) || length(f) == 0 || f == "") stop(desc, " path is empty.")
    if (!file.exists(f)) stop(desc, " not found: ", f)

    d_dt <- data.table::fread(f, header = header, data.table = FALSE, showProgress = FALSE, fill = Inf)

    if (!is.null(row_name) && ncol(d_dt) > 1) {
        rownames(d_dt) <- d_dt[, row_name]
        d_dt <- d_dt[, -row_name, drop = FALSE]
    }

    if (ncol(d_dt) < min_cols) {
        stop(desc, " has insufficient columns (found ", ncol(d_dt), ", required ", min_cols, ").")
    }
    return(d_dt)
}

#' @title Extract Target Gene Sets from Annotation Results
#' @description Parses loop and target annotations to extract valid gene lists.
#' @param annotation_res List. Annotation result from
#'   \code{annotate_peaks_and_loops} or
#'   \code{refine_loop_anchors_by_expression}.
#' @param src Character vector. One or both of \code{"targets"} and
#'   \code{"loops"} to select the annotation source.
#' @param active_loop_types Character vector or \code{NULL}. Loop types to
#'   include when \code{src} includes \code{"loops"}. \code{NULL} uses all.
#' @param include_Filled Logical. Whether to use the \code{*_Filled} (linear
#'   nearest-gene fallback-augmented) columns for the \code{"targets"} branch.
#'   Default \code{TRUE}.
#' @param use_nearest_gene Logical. If \code{TRUE}, bypasses 3D-loop
#'   assignment and uses the linear nearest-gene column directly.
#'   Default \code{FALSE}.
#' @param target_mapping_mode Character. \code{"all"} or
#'   \code{"promoter"}. Default \code{"all"}.
#' @return A named list of character vectors, each containing target gene symbols.
#' @keywords internal
#' @noRd
extract_target_gene_sets <- function(annotation_res, src, active_loop_types = NULL, include_Filled = TRUE, use_nearest_gene = FALSE, target_mapping_mode = "all") {
    raw_gene_sets <- list()

    if ("targets" %in% src && !is.null(annotation_res$target_annotation)) {
        bed_info <- annotation_res$target_annotation
        target_col <- NULL
        if (use_nearest_gene && target_mapping_mode == "promoter") {
            stop(
                "`use_nearest_gene = TRUE` bypasses promoter-specific ",
                "loop mapping. Use `target_mapping_mode = 'all'` or ",
                "disable nearest-gene mode.",
                call. = FALSE
            )
        }
        if (use_nearest_gene) {
            if ("SYMBOL" %in% colnames(bed_info)) {
                target_col <- "SYMBOL"
            } else if ("geneId" %in% colnames(bed_info)) target_col <- "geneId"
            if (is.null(target_col)) stop("Targets: 'SYMBOL' or 'geneId' required when use_nearest_gene is TRUE.")
        } else {
            base_col <- if (target_mapping_mode == "promoter") "Regulated_promoter_genes" else "Assigned_Target_Genes"
            desired_col <- if (include_Filled) paste0(base_col, "_Filled") else base_col
            if (desired_col %in% colnames(bed_info)) {
                target_col <- desired_col
            } else {
                stop("Targets: Required column '", desired_col, "' not found.")
            }
        }
        if (!is.null(target_col)) {
            gs <- clean_gene_names(bed_info[[target_col]], "[;,]")
            if (length(gs) > 0) raw_gene_sets[["Target_Genes"]] <- gs
        }
    }

    if ("loops" %in% src && !is.null(annotation_res$loop_annotation) &&
        nrow(annotation_res$loop_annotation) > 0) {
        loop_df <- annotation_res$loop_annotation
        # Map: all ->Putative_Target_Genes, promoter ->Promoter_Target_Genes.
        # Expression filtering is encoded by the refinement stage of the
        # supplied annotation object, not by a separate mode.
        gene_col <- switch(target_mapping_mode,
            all = "Putative_Target_Genes",
            promoter = "Promoter_Target_Genes"
        )
        if (!gene_col %in% colnames(loop_df)) {
            stop("Loops: Required column '", gene_col, "' not found. ",
                "Ensure the annotation object was produced by the current ",
                "pipeline version.",
                call. = FALSE
            )
        }
        use_types <- if (is.null(active_loop_types)) unique(loop_df$loop_type) else intersect(active_loop_types, unique(loop_df$loop_type))
        if (length(use_types) > 0) {
            for (lt in use_types) {
                sub_df <- loop_df[loop_df$loop_type == lt, ]
                if (nrow(sub_df) > 0) {
                    gs <- clean_gene_names(sub_df[[gene_col]], "[;,]")
                    if (length(gs) > 0) {
                        safe_name <- paste0(gsub("-", "", lt, fixed = TRUE), "_Genes")
                        raw_gene_sets[[safe_name]] <- gs
                    }
                }
            }
        }
    }
    return(raw_gene_sets)
}

#' @title Run per-task profiling pipeline (violin, GSEA, heatmap, motif, GO, PPI)
#' @return A list with \code{analysis_queue}, \code{go_results},
#'   \code{motif_results}, \code{plots}, and \code{warnings}.
#' @keywords internal
#' @noRd
.run_profile_tasks <- function(
  analysis_queue, current_source_proj_name, global_glist,
  tpm_mat_raw, meta_raw, loop_stats_df, annotation_res, src,
  expr_vals = NULL,
  stat_test, gsea_nSample, heatmap_nSample, cor_method,
  run_motif, genome_id, motif_p_thresh, motif_ntop, motif_n_perm,
  run_go, universe_genes, org_db, cnet_nSample,
  run_ppi, ppi_score, ppi_nSample, ppi_species_id
) {
    go_results <- list()
    motif_results <- list()
    plots <- list()
    warn_env <- new.env(parent = emptyenv())
    warn_env$warnings <- character()

    .safe_run <- function(module_name, expr) {
        out <- tryCatch(
            list(result = expr, warning = NULL),
            error = function(e) {
                w_msg <- paste0("[", module_name, "] failed: ", conditionMessage(e))
                warning(w_msg, call. = FALSE)
                list(result = NULL, warning = w_msg)
            }
        )
        if (!is.null(out$warning)) {
            warn_env$warnings <- c(warn_env$warnings, out$warning)
        }
        out$result
    }

    for (task_name in names(analysis_queue)) {
        target_genes <- analysis_queue[[task_name]]
        current_proj_name <- paste0(current_source_proj_name, "_", task_name)
        idx <- match(toupper(target_genes), toupper(names(global_glist)))
        target_genes <- unique(names(global_glist)[idx[!is.na(idx)]])
        analysis_queue[[task_name]] <- target_genes

        message("\n--- Task: ", task_name, " (Valid Genes: ", length(target_genes), ") ---")
        if (length(target_genes) < 3) {
            message("  Too few genes (<3), skipping.")
            next
        }

        mapped_upper <- toupper(target_genes)
        if (any(duplicated(mapped_upper))) {
            dup_upper <- unique(mapped_upper[duplicated(mapped_upper)])
            collided <- unique(target_genes[toupper(target_genes) %in% dup_upper])
            if (length(collided) > 0) {
                warning(
                    "Task '", task_name, "': ", length(dup_upper),
                    " gene symbol(s) collide after case-insensitive matching (",
                    "e.g., ", paste(head(collided, 6), collapse = ", "), "). ",
                    "Only the first matching symbol from the ranked list is used. ",
                    "Consider normalising gene identifier case before analysis.",
                    call. = FALSE
                )
            }
        }
        task_plots <- list()

        # Violin
        p_vio <- .safe_run(
            "Violin",
            run_lfc_violin(target_genes, global_glist, stat_test, current_proj_name,
                expr_vals = expr_vals
            )
        )
        if (!is.null(p_vio)) task_plots$LFC_Violin <- p_vio

        # GSEA
        gsea_out <- .safe_run(
            "GSEA",
            run_gsea_analysis(target_genes, global_glist, gsea_nSample, current_proj_name)
        )
        if (!is.null(gsea_out) && !is.null(gsea_out$plot)) task_plots$GSEA <- gsea_out$plot

        # Connectivity heatmap (total loops)
        heat_plots <- .safe_run(
            "Heatmap",
            run_heatmap_and_connectivity(
                target_genes, tpm_mat_raw, meta_raw, loop_stats_df,
                global_glist, heatmap_nSample, cor_method, current_proj_name,
                source_type = src, target_col = NULL
            )
        )
        if (!is.null(heat_plots) && length(heat_plots) > 0) task_plots <- c(task_plots, heat_plots)

        # Connectivity (distal loops)
        if (!is.null(loop_stats_df) && "n_Linked_Distal" %in% colnames(loop_stats_df)) {
            dist_plots <- .safe_run(
                "DistalHeatmap",
                run_heatmap_and_connectivity(
                    target_genes, tpm_mat_raw, meta_raw, loop_stats_df,
                    global_glist, heatmap_nSample, cor_method, current_proj_name,
                    source_type = src, target_col = "n_Linked_Distal", skip_heatmap = TRUE
                )
            )
            if (!is.null(dist_plots) && length(dist_plots) > 0) task_plots <- c(task_plots, dist_plots)
        }

        # Motif
        if (run_motif) {
            motif_loop_df <- .subset_motif_loop_df(annotation_res$loop_annotation, src, task_name)
            motif_out <- .safe_run(
                "Motif",
                run_distal_motif_analysis(
                    target_genes, motif_loop_df,
                    genome_id, motif_p_thresh, current_proj_name, motif_ntop,
                    motif_n_perm = motif_n_perm,
                    anchor_registry = .get_anchor_registry(annotation_res)
                )
            )
            if (!is.null(motif_out)) {
                if (!is.null(motif_out$plots) && length(motif_out$plots) > 0) {
                    task_plots <- c(task_plots, motif_out$plots)
                }
                if (!is.null(motif_out$results)) {
                    motif_results[[task_name]] <- motif_out$results
                }
            }
        }

        # GO
        if (run_go) {
            go_out <- .safe_run(
                "GO",
                run_go_enrichment(
                    target_genes, org_db,
                    universe_genes = if (!is.null(universe_genes)) universe_genes else global_glist,
                    cnet_nSample = cnet_nSample,
                    project_name = current_proj_name
                )
            )
            if (!is.null(go_out) && !is.null(go_out$result) && nrow(go_out$result) > 0) {
                top_go <- if ("ONTOLOGY" %in% colnames(go_out$result)) {
                    go_out$result %>%
                        dplyr::group_by(ONTOLOGY) %>%
                        dplyr::arrange(p.adjust) %>%
                        dplyr::slice_head(n = 5) %>%
                        dplyr::ungroup()
                } else {
                    head(go_out$result[order(go_out$result$p.adjust), ], 15)
                }
                top_go$CleanLoopType <- task_name
                top_go$LoopType <- if ("ONTOLOGY" %in% colnames(go_out$result)) {
                    paste0(task_name, "\n(", top_go$ONTOLOGY, ")")
                } else {
                    task_name
                }
                top_go$Source <- src
                go_results[[length(go_results) + 1]] <- top_go
            }
            if (!is.null(go_out) && !is.null(go_out$plot)) task_plots$GO_Network <- go_out$plot
        }

        # PPI
        if (run_ppi) {
            p_ppi <- .safe_run(
                "PPI",
                run_ppi_analysis(
                    target_genes, global_glist, org_db, ppi_score,
                    ppi_nSample, current_proj_name, ppi_species_id
                )
            )
            if (!is.null(p_ppi)) task_plots$PPI_Network <- p_ppi
        }

        plots[[task_name]] <- task_plots
    }

    list(
        analysis_queue = analysis_queue, go_results = go_results,
        motif_results = motif_results, plots = plots,
        warnings = warn_env$warnings
    )
}

#' @title Generate LFC Violin and Boxplot
#' @return A \code{ggplot} object, or \code{NULL} if fewer than 3 valid targets.
#' @keywords internal
#' @noRd
run_lfc_violin <- function(target_genes, global_glist, stat_test = c("wilcox.test", "t.test"), project_name, expr_vals = NULL, seed = NULL, n_iter = 100L) {
    stat_test <- match.arg(stat_test)
    # Build case-insensitive key for matching: human SYMBOL is typically
    # uppercase, but mouse Title Case and mixed conventions also appear.
    # Canonicalise both sides to uppercase for intersection, then map back.
    lfc_names_upper <- toupper(names(global_glist))
    target_upper <- toupper(target_genes)
    valid_idx <- which(target_upper %in% lfc_names_upper & !duplicated(target_upper))
    valid_targets <- target_genes[valid_idx]
    if (length(valid_targets) < 3) {
        return(NULL)
    }

    # Map uppercase back to original LFC names for value lookup
    lfc_lookup <- setNames(names(global_glist), lfc_names_upper)
    target_lfc <- global_glist[lfc_lookup[toupper(valid_targets)]]

    # Cliff's delta: non-parametric effect size, [-1, 1]
    # |d| >= 0.474 = large, >= 0.33 = medium, >= 0.147 = small
    cliff_d <- function(x, y) {
        r <- rank(c(x, y))
        R1 <- sum(r[seq_along(x)])
        U1 <- R1 - length(x) * (length(x) + 1) / 2
        2 * U1 / (length(x) * length(y)) - 1
    }

    # --- Helper: draw one matched background ---
    # RNG state is managed externally via withr::local_seed()
    .draw_matched_background <- function() {
        if (!is.null(expr_vals) && length(expr_vals) >= 10) {
            expr_names_upper <- toupper(names(expr_vals))
            expr_lookup <- setNames(names(expr_vals), expr_names_upper)
            common_upper <- base::intersect(lfc_names_upper, expr_names_upper)
            common <- lfc_lookup[common_upper]
            if (length(common) > length(valid_targets) * 2) {
                expr_common <- expr_vals[expr_lookup[common_upper]]
                names(expr_common) <- common
                decile_breaks <- stats::quantile(expr_common,
                    probs = seq(0, 1, 0.1), na.rm = TRUE
                )
                decile_breaks <- unique(decile_breaks)
                if (length(decile_breaks) > 1) {
                    target_decile <- as.integer(cut(
                        expr_vals[expr_lookup[toupper(valid_targets)]],
                        breaks = decile_breaks, include.lowest = TRUE
                    ))
                    all_decile <- as.integer(cut(
                        expr_common,
                        breaks = decile_breaks, include.lowest = TRUE
                    ))
                    k_per_target <- 10L
                    valid_target_keys <- toupper(valid_targets)
                    matched_bg <- unlist(lapply(
                        sort(unique(stats::na.omit(target_decile))),
                        function(d) {
                            n_target_d <- sum(target_decile == d, na.rm = TRUE)
                            pool_upper <- setdiff(common_upper[all_decile == d], valid_target_keys)
                            if (length(pool_upper) == 0) {
                                return(character(0))
                            }
                            pool <- lfc_lookup[pool_upper]
                            n_sample <- min(k_per_target * n_target_d, length(pool))
                            sample(pool, size = n_sample, replace = FALSE)
                        }
                    ))
                    if (length(matched_bg) >= length(valid_targets)) {
                        return(global_glist[matched_bg])
                    }
                }
            }
        }
        # Fallback: full background (no expression matching possible)
        global_glist[setdiff(names(global_glist), lfc_lookup[toupper(valid_targets)])]
    }

    # --- Iterate matched sampling to account for matching uncertainty ---
    # Pre-generate deterministic seeds so that the representative background
    # draw for the violin plot exactly reproduces the best iteration.
    base_seed <- if (!is.null(seed)) {
        as.integer(seed)
    } else {
        sample.int(.Machine$integer.max, 1L)
    }
    draw_seeds <- as.integer(
        (as.double(base_seed) + seq_len(n_iter) - 1L) %%
            (.Machine$integer.max - 1) + 1
    )

    if (n_iter > 1L && !is.null(expr_vals) && length(expr_vals) >= 10) {
        iter_d <- numeric(n_iter)
        iter_md <- numeric(n_iter)
        for (i in seq_len(n_iter)) {
            withr::local_seed(draw_seeds[i])
            other_lfc <- .draw_matched_background()
            iter_d[i] <- cliff_d(target_lfc, other_lfc)
            iter_md[i] <- median(target_lfc) - median(other_lfc)
        }
        d <- stats::median(iter_d, na.rm = TRUE)
        ci_d <- stats::quantile(iter_d, c(0.025, 0.975), na.rm = TRUE)
        med_diff <- stats::median(iter_md, na.rm = TRUE)
        ci_m <- stats::quantile(iter_md, c(0.025, 0.975), na.rm = TRUE)
        n_iter_msg <- sprintf(" (n_iter=%d)", n_iter)
    } else {
        # Single draw (backward compatible or no expression data)
        other_lfc <- .draw_matched_background()
        d <- cliff_d(target_lfc, other_lfc)
        ci_d <- c(NA_real_, NA_real_)
        med_diff <- median(target_lfc) - median(other_lfc)
        ci_m <- c(NA_real_, NA_real_)
        n_iter_msg <- ""
    }

    # --- Representative background for the violin plot ---
    # Use the iteration whose Cliff's delta is closest to the median, with
    # the exact same seed so the visual distribution matches the effect size.
    if (n_iter > 1L && !is.null(expr_vals) && length(expr_vals) >= 10) {
        best_i <- which.min(abs(iter_d - d))[1]
        withr::local_seed(draw_seeds[best_i])
        plot_bg <- .draw_matched_background()
    } else {
        plot_bg <- .draw_matched_background()
    }
    plot_data <- data.frame(
        LFC = c(target_lfc, plot_bg),
        Group = factor(c(rep("Target", length(target_lfc)), rep("Background", length(plot_bg))), levels = c("Target", "Background"))
    )

    n_target <- length(target_lfc)
    n_back <- length(plot_bg)

    p_val <- if (stat_test == "wilcox.test") stats::wilcox.test(target_lfc, plot_bg)$p.value else stats::t.test(target_lfc, plot_bg)$p.value

    # Bootstrap CI for single-draw mode; iteration-based CI for multi-draw mode
    if (n_iter <= 1L) {
        if (!is.null(seed)) withr::local_seed(seed + 1L)
        n_boot <- 1000L
        boot_d <- replicate(n_boot, cliff_d(
            sample(target_lfc, replace = TRUE),
            sample(plot_bg, replace = TRUE)
        ))
        boot_med <- replicate(
            n_boot,
            median(sample(target_lfc, replace = TRUE)) -
                median(sample(plot_bg, replace = TRUE))
        )
        ci_d <- stats::quantile(boot_d, c(0.025, 0.975), na.rm = TRUE)
        ci_m <- stats::quantile(boot_med, c(0.025, 0.975), na.rm = TRUE)
    }

    p_label <- sprintf(
        "Cliff's Delta = %.3f [%.3f, %.3f]  DeltaMedian = %.3f [%.3f, %.3f]  P = %s%s",
        d, ci_d[1], ci_d[2],
        med_diff, ci_m[1], ci_m[2],
        formatC(p_val, format = "e", digits = 2),
        n_iter_msg
    )
    x_labels <- c("Target" = paste0("Target\n(n=", n_target, ")"), "Background" = paste0("Background\n(n=", n_back, ")"))

    y_min <- quantile(plot_data$LFC, 0.001, na.rm = TRUE)
    y_max <- quantile(plot_data$LFC, 0.999, na.rm = TRUE)
    y_pad <- (y_max - y_min) * 0.03
    cols <- c("Target" = "#E41A1C", "Background" = "#999999")

    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Group, y = LFC, fill = Group)) +
        ggplot2::geom_violin(trim = TRUE, alpha = 0.5, color = NA, scale = "width", width = 0.55, adjust = 2.5) +
        ggplot2::geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.9, color = "black") +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey30", linewidth = 0.5) +
        ggplot2::scale_x_discrete(labels = x_labels) +
        ggplot2::scale_y_continuous(expand = c(0, 0)) +
        ggplot2::coord_cartesian(ylim = c(y_min - y_pad, y_max + y_pad)) +
        ggplot2::scale_fill_manual(values = cols) +
        ggplot2::labs(title = project_name, subtitle = paste0("Stat: ", stat_test, ", P: ", p_label), y = "Log2 Fold Change", x = NULL) +
        ggplot2::theme_classic() +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 12), plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10), legend.position = "none")

    return(p)
}

#' @title Run Custom Gene Set Enrichment Analysis (GSEA)
#' @param target_genes Character vector of target gene symbols.
#' @param global_glist Named numeric vector (gene-level ranking metric,
#'   e.g. log2 fold change). Names are gene symbols; values are the
#'   ranking statistic sorted descending.
#' @param gsea_nSample Integer or \code{NULL}. Maximum number of target
#'   genes to sample for enrichment. \code{NULL} uses all targets.
#' @param current_proj_name Character. Prefix for plot titles and
#'   file names.
#' @return A list with \code{result} (data frame) and \code{plot} (ggplot) elements.
#' @keywords internal
#' @noRd
run_gsea_analysis <- function(target_genes, global_glist, gsea_nSample, current_proj_name) {
    curr_glist <- global_glist
    if (any(duplicated(curr_glist))) {
        curr_glist <- curr_glist + seq_along(curr_glist) * 1e-12
        curr_glist <- sort(curr_glist, decreasing = TRUE)
    }
    # Case-insensitive matching: detect collisions before normalising
    upper_gl <- toupper(names(curr_glist))
    upper_tg <- toupper(target_genes)
    # Deduplicate gene names in ranked list: collisions after case
    # normalisation keep the entry with the larger absolute value.
    names(curr_glist) <- upper_gl
    if (any(duplicated(names(curr_glist)))) {
        dup_n <- sum(duplicated(names(curr_glist)))
        warning(
            "GSEA ranked list: ", dup_n,
            " gene symbol(s) collide after toupper(). Keeping the entry with ",
            "larger absolute value per gene. Consider normalising gene ",
            "identifier case before analysis.",
            call. = FALSE
        )
        curr_glist <- curr_glist[order(abs(curr_glist), decreasing = TRUE)]
        curr_glist <- curr_glist[!duplicated(names(curr_glist))]
        curr_glist <- sort(curr_glist, decreasing = TRUE)
    }
    curr_targets <- unique(upper_tg)
    curr_targets <- intersect(curr_targets, names(curr_glist))

    if (!is.null(gsea_nSample) && length(curr_targets) > gsea_nSample) {
        warning("GSEA: down-sampling ", gsea_nSample, " of ", length(curr_targets),
            " target genes. GSEA results represent a random subset, ",
            "not the full gene set. Set gsea_nSample = NULL for full analysis. ",
            "For fully reproducible results, set the `seed` parameter in ",
            "profile_target_genes() (e.g., seed = 42).",
            call. = FALSE
        )
        curr_targets <- sample(curr_targets, size = gsea_nSample, replace = FALSE)
    }
    if (length(curr_targets) < 2) {
        return(list(result = NULL, plot = NULL))
    }

    term_df <- data.frame(
        term = current_proj_name,
        gene = curr_targets,
        stringsAsFactors = FALSE
    )
    gsea_res <- tryCatch(
        # GSEA in profile_target_genes is exploratory: pvalueCutoff=1.1 and
        # minGSSize=10 ensure permissive but stable enrichment.
        clusterProfiler::GSEA(curr_glist,
            TERM2GENE = term_df,
            pvalueCutoff = 1.1, minGSSize = 10, maxGSSize = 50000,
            verbose = FALSE, seed = TRUE
        ),
        error = function(e) {
            warning("GSEA failed for ", current_proj_name, ": ",
                conditionMessage(e),
                call. = FALSE
            )
            return(NULL)
        }
    )

    if (is.null(gsea_res) || nrow(as.data.frame(gsea_res)) == 0) {
        return(list(result = NULL, plot = NULL))
    }

    p_out <- NULL
    p_temp <- tryCatch(
        .with_known_upstream_noise_suppressed(
            enrichplot::gseaplot2(gsea_res, geneSetID = 1, subplots = 1)
        ),
        error = function(e) NULL
    )
    d <- NULL
    if (inherits(p_temp, "ggplot")) {
        d <- p_temp$data
    } else if (inherits(p_temp, "aplot") || inherits(p_temp, "gglist") || is.list(p_temp)) {
        for (sub_p in p_temp) {
            if (inherits(sub_p, "ggplot") && !is.null(sub_p$data) && "runningScore" %in% colnames(sub_p$data)) {
                d <- sub_p$data
                break
            }
        }
    } else if (!is.null(p_temp$data)) d <- p_temp$data

    if (!is.null(d) && is.data.frame(d) && "runningScore" %in% colnames(d)) {
        if (!"geneList" %in% colnames(d)) d$geneList <- curr_glist[d$x]
        if (!"position" %in% colnames(d)) d$position <- as.numeric(names(curr_glist)[d$x] %in% curr_targets)
        max_rank <- max(d$x)
        res_df <- as.data.frame(gsea_res)
        nes_val <- res_df$NES[1]
        pval_val <- res_df$pvalue[1]
        main_col <- if (!is.na(nes_val) && nes_val >= 0) "#E41A1C" else "#377EB8"

        p1 <- ggplot2::ggplot(d, ggplot2::aes(x = x, y = runningScore)) +
            ggplot2::geom_line(color = main_col, linewidth = 1) +
            ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
            ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(0, max_rank)) +
            ggplot2::theme_bw() +
            ggplot2::labs(x = NULL, y = "ES", title = paste0(current_proj_name, "\nNES: ", round(nes_val, 3), "  nominal P: ", formatC(pval_val, format = "e", digits = 2))) +
            ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(), panel.grid = ggplot2::element_blank(), plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 10))

        hit_data <- d[d$position == 1, ]
        p2 <- ggplot2::ggplot(hit_data, ggplot2::aes(x = x, y = 1)) +
            ggplot2::geom_segment(ggplot2::aes(xend = x, yend = 0), color = "black", alpha = 0.6) +
            ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(0, max_rank)) +
            ggplot2::scale_y_continuous(expand = c(0, 0)) +
            ggplot2::theme_void() +
            ggplot2::theme(panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5))

        p3 <- ggplot2::ggplot(d, ggplot2::aes(x = x, y = geneList)) +
            ggplot2::geom_segment(ggplot2::aes(xend = x, yend = 0, color = geneList)) +
            ggplot2::scale_color_gradient2(low = "#1B7837", mid = "white", high = "#762A83", midpoint = 0, limits = c(-3, 3), oob = scales::squish) +
            ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(0, max_rank)) +
            ggplot2::coord_cartesian(ylim = c(quantile(d$geneList, 0.005, na.rm = TRUE), quantile(d$geneList, 0.995, na.rm = TRUE))) +
            ggplot2::theme_classic() +
            ggplot2::labs(x = "Rank", y = "LFC") +
            ggplot2::theme(legend.position = "none", axis.text.y = ggplot2::element_text(size = 8))

        if (requireNamespace("aplot", quietly = TRUE)) {
            p_out <- aplot::plot_list(p1, p2, p3, ncol = 1, heights = c(2, 0.5, 1.5))
        } else {
            p_out <- .with_known_upstream_noise_suppressed(
                enrichplot::gseaplot2(gsea_res, geneSetID = 1, title = as.data.frame(gsea_res)$Description[1])
            )
        }
    }
    return(list(result = as.data.frame(gsea_res), plot = p_out))
}

#' @title Perform GO Enrichment and Generate Network Plot
#' @param genes Character vector of gene symbols (or ENTREZ IDs on
#'   second attempt after mapping failure) to test for GO over-representation.
#' @param org_db Character. Organism annotation database package name
#'   (e.g. \code{"org.Hs.eg.db"}).
#' @param universe_genes Named numeric vector (gene-level ranking metric)
#'   used as the background universe. Names must match the \code{genes}
#'   identifier convention.
#' @param cnet_nSample Integer. Number of top GO terms to display in the
#'   cnetplot network. Default \code{50}.
#' @param project_name Character. Prefix for plot titles. Default
#'   \code{"Analysis"}.
#' @importFrom methods slot<-
#' @return A list with \code{result} (data frame) and \code{plot} (ggplot) elements.
#' @keywords internal
#' @noRd
run_go_enrichment <- function(genes, org_db, universe_genes, cnet_nSample = 50, project_name = "Analysis") {
    clean_genes <- clean_gene_names(genes)
    org_db_obj <- .get_org_db_obj(org_db)
    valid_keys <- AnnotationDbi::keytypes(org_db_obj)
    primary_key <- if ("ENTREZID" %in% valid_keys) "ENTREZID" else valid_keys[1]
    symbol_key <- if ("SYMBOL" %in% valid_keys) "SYMBOL" else valid_keys[1]

    gene_entrez <- .with_known_upstream_noise_suppressed(AnnotationDbi::mapIds(
        org_db_obj,
        keys = clean_genes,
        column = primary_key,
        keytype = symbol_key,
        multiVals = "first"
    ))
    # Count genes dropped due to multi-mapping collisions
    gene_entrez_list <- .with_known_upstream_noise_suppressed(AnnotationDbi::mapIds(
        org_db_obj,
        keys = clean_genes, column = primary_key,
        keytype = symbol_key, multiVals = "list"
    ))
    n_multi <- sum(vapply(gene_entrez_list, length, integer(1)) > 1, na.rm = TRUE)
    if (n_multi > 0) {
        message(sprintf("GO: %d gene(s) have multiple ENTREZID mappings (first kept).", n_multi))
    }
    valid_entrez <- na.omit(gene_entrez)
    n_dropped <- length(clean_genes) - length(valid_entrez)
    if (n_dropped > 0) {
        message(sprintf(
            "GO: %d of %d genes (%.1f%%) could not be mapped to ENTREZID.",
            n_dropped, length(clean_genes),
            100 * n_dropped / length(clean_genes)
        ))
    }

    use_symbol_mode <- length(valid_entrez) < 5 || (length(valid_entrez) / length(clean_genes) < 0.1)

    if (use_symbol_mode) {
        final_genes <- clean_genes
        final_keytype <- symbol_key
    } else {
        final_genes <- valid_entrez
        final_keytype <- primary_key
    }

    final_universe <- NULL
    if (!is.null(universe_genes)) {
        if (use_symbol_mode) {
            final_universe <- names(universe_genes)
        } else {
            univ_entrez <- .with_known_upstream_noise_suppressed(AnnotationDbi::mapIds(
                org_db_obj,
                keys = names(universe_genes),
                column = primary_key,
                keytype = symbol_key,
                multiVals = "first"
            ))
            final_universe <- unique(na.omit(univ_entrez))
            if (length(universe_genes) > 0 &&
                length(final_universe) / length(universe_genes) < 0.5) {
                warning(
                    "Only ", round(length(final_universe) / length(universe_genes) * 100, 1),
                    "% of background genes mapped to ENTREZID. ",
                    "GO enrichment background may be incomplete.",
                    call. = FALSE
                )
            }
        }
    }

    ego <- tryCatch(
        clusterProfiler::enrichGO(gene = final_genes, universe = final_universe, OrgDb = org_db_obj, keyType = final_keytype, ont = "ALL", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.1, minGSSize = 10, maxGSSize = 800, readable = (final_keytype == primary_key)),
        error = function(e) {
            warning("clusterProfiler::enrichGO failed: ", conditionMessage(e), call. = FALSE)
            return(NULL)
        }
    )
    if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
        return(list(result = NULL, plot = NULL))
    }

    p_cnet <- NULL
    top_n <- if (!is.null(cnet_nSample)) {
        max(1L, min(as.integer(cnet_nSample), nrow(as.data.frame(ego))))
    } else {
        min(5L, nrow(as.data.frame(ego)))
    }
    fc_vec <- universe_genes
    # Convert fc_vec names to ENTREZ only when the GO result is in ENTREZ
    # (not readable mode).  When readable=TRUE, GO geneID is already SYMBOL
    # and fc_vec names should remain SYMBOL for intersect() to work.
    readable_mode <- final_keytype == primary_key # enrichGO readable=TRUE
    if (!use_symbol_mode && !readable_mode &&
        exists("univ_entrez", inherits = FALSE)) {
        valid_map <- univ_entrez[!is.na(univ_entrez)]
        name_idx <- match(names(fc_vec), names(valid_map))
        has_map <- !is.na(name_idx)
        names(fc_vec)[has_map] <- as.character(valid_map[name_idx[has_map]])
    }

    genes_to_label <- c()
    top_df <- head(as.data.frame(ego), top_n)

    gene_to_pathways <- list()
    for (i in seq_len(nrow(top_df))) {
        gs <- unlist(strsplit(top_df$geneID[i], "/", fixed = TRUE))
        for (g in gs) gene_to_pathways[[g]] <- c(gene_to_pathways[[g]], top_df$ID[i])
    }

    for (i in seq_len(nrow(top_df))) {
        gs <- unlist(strsplit(top_df$geneID[i], "/", fixed = TRUE))
        valid_g <- intersect(gs, names(fc_vec))
        if (length(valid_g) > 0) genes_to_label <- c(genes_to_label, head(valid_g[order(abs(fc_vec[valid_g]), decreasing = TRUE)], 3))
    }

    hub_genes <- names(gene_to_pathways)[lengths(gene_to_pathways) >= 2]
    valid_hub <- intersect(hub_genes, names(fc_vec))
    if (length(valid_hub) > 0) genes_to_label <- c(genes_to_label, head(valid_hub[order(abs(fc_vec[valid_hub]), decreasing = TRUE)], 5))

    genes_to_label <- unique(genes_to_label)

    ego_df <- as.data.frame(ego)
    ego_df$Description <- vapply(ego_df$Description, function(x) paste(strwrap(x, width = 35), collapse = "\n"), FUN.VALUE = character(1))
    slot(ego, "result") <- ego_df

    old_ggrepel <- getOption("ggrepel.max.overlaps", 10)
    options(ggrepel.max.overlaps = 100)
    on.exit(options(ggrepel.max.overlaps = old_ggrepel), add = TRUE)
    p_cnet <- .with_known_upstream_noise_suppressed(
        enrichplot::cnetplot(ego, foldChange = fc_vec, showCategory = top_n, node_label = "category")
    )

    if (length(genes_to_label) > 0 && requireNamespace("ggraph", quietly = TRUE)) {
        p_cnet <- p_cnet + ggraph::geom_node_text(ggplot2::aes(filter = name %in% genes_to_label, label = name), repel = TRUE, size = 3.5, fontface = "bold.italic", bg.color = "white", bg.r = 0.15, max.overlaps = Inf)
    }

    p_cnet <- .with_known_upstream_noise_suppressed(
        p_cnet + ggplot2::scale_color_distiller(palette = "PuOr", name = "Log2FC") + ggplot2::labs(title = paste0("GO Network: ", project_name)) + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))
    )

    return(list(result = as.data.frame(ego), plot = p_cnet))
}

#' @title Construct and Visualize STRING PPI Network
#' @param target_genes Character vector of target gene symbols for PPI
#'   network construction.
#' @param global_glist Named numeric vector (gene-level ranking metric).
#'   Used to colour nodes by LFC and rank by connectivity.
#' @param org_db Character. Organism annotation database package name
#'   (e.g. \code{"org.Hs.eg.db"}). Used to resolve STRING species ID.
#' @param ppi_score Numeric. Minimum STRING combined confidence score
#'   for edge inclusion. Default \code{400}.
#' @param ppi_ntop Integer. Maximum number of top-ranked target genes
#'   to include in the PPI construction. Default \code{400}.
#' @param current_proj_name Character. Prefix for plot titles.
#' @details
#' Nodes are the top \code{ppi_ntop} target genes by \code{abs(LFC)}; edges
#' are derived from the STRING knowledge base (combined score >=
#' \code{ppi_score}). \strong{Interpretation note:} node degree in the PPI
#' network reflects known protein-protein interactions from STRING, not
#' direct evidence of 3D chromatin connectivity. Hub genes in this network
#' should not be conflated with regulatory hubs in the loop-anchor graph.
#' @importFrom utils capture.output
#' @return A \code{ggplot} object representing the PPI network, or \code{NULL} if no interactions found.
#' @keywords internal
#' @noRd
run_ppi_analysis <- function(target_genes, global_glist, org_db, ppi_score, ppi_ntop, current_proj_name, ppi_species_id = NULL) {
    # Resolve STRING species ID from OrgDb package name.
    # If ppi_species_id is explicitly provided (NCBI taxonomy ID), use it
    # directly. Otherwise, try inference from the OrgDb package name.
    if (!is.null(ppi_species_id)) {
        if (!is.numeric(ppi_species_id) || length(ppi_species_id) != 1L || is.na(ppi_species_id)) {
            stop("`ppi_species_id` must be a single NCBI taxonomy ID integer", call. = FALSE)
        }
        species_id <- as.integer(ppi_species_id)
    } else {
        species_map <- c(
            "org.Hs.eg.db" = 9606L, # human
            "org.Mm.eg.db" = 10090L, # mouse
            "org.Rn.eg.db" = 10116L, # rat
            "org.Dm.eg.db" = 7227L, # fruit fly
            "org.Ce.eg.db" = 6239L, # worm
            "org.Sc.eg.db" = 4932L, # yeast
            "org.Dr.eg.db" = 7955L # zebrafish
        )
        org_pkg_name <- if (is.character(org_db) && length(org_db) == 1L && nzchar(org_db)) {
            org_db
        } else {
            attr(org_db, "package")
        }
        species_id <- species_map[[org_pkg_name]]
        if (is.null(species_id)) {
            # Pattern-based detection for known-but-unlisted OrgDbs
            species_id <- if (grepl("\\.[Mm]m\\.", org_pkg_name)) {
                10090L
            } else if (grepl("\\.[Rr]n\\.", org_pkg_name)) {
                10116L
            } else if (grepl("\\.[Dd]m\\.", org_pkg_name)) {
                7227L
            } else {
                NULL
            }
            if (is.null(species_id)) {
                warning("Could not resolve STRING species ID for OrgDb '", org_pkg_name,
                    "'. PPI analysis skipped. Supported OrgDbs include: ",
                    "org.Hs.eg.db, org.Mm.eg.db, org.Rn.eg.db, org.Dm.eg.db.",
                    call. = FALSE
                )
                return(NULL)
            }
            warning("Inferred STRING species ID ", species_id, " for '", org_pkg_name,
                "'. If this is incorrect, verify the OrgDb package.",
                call. = FALSE
            )
        }
    } # end else (ppi_species_id not provided)

    string_db_obj <- tryCatch(
        suppressMessages(
            STRINGdb::STRINGdb$new(
                species = species_id,
                score_threshold = ppi_score,
                version = "12.0",
                input_directory = tempdir()
            )
        ),
        error = function(e) {
            warning("STRINGdb initialisation failed: ", conditionMessage(e),
                ". Skipping PPI analysis.",
                call. = FALSE
            )
            return(NULL)
        }
    )
    if (is.null(string_db_obj)) {
        return(NULL)
    }

    ppi_genes <- target_genes
    if (!is.null(ppi_ntop) && length(ppi_genes) > ppi_ntop) {
        valid_in_lfc <- intersect(ppi_genes, names(global_glist))
        if (length(valid_in_lfc) > 0) {
            ppi_genes <- head(valid_in_lfc[order(abs(global_glist[valid_in_lfc]), decreasing = TRUE)], ppi_ntop)
        } else {
            ppi_genes <- head(ppi_genes, ppi_ntop)
        }
    }

    # Capture STRINGdb mapping output to prevent raw text in the report
    targets_mapped <- tryCatch(
        suppressMessages({
            capture.output(
                out <- string_db_obj$map(
                    data.frame(gene = ppi_genes), "gene",
                    removeUnmappedRows = TRUE
                )
            )
            out
        }),
        error = function(e) {
            warning("STRINGdb gene mapping failed: ", conditionMessage(e),
                ". Skipping PPI analysis.",
                call. = FALSE
            )
            return(data.frame())
        }
    )
    if (nrow(targets_mapped) == 0) {
        return(NULL)
    }

    hits <- targets_mapped$STRING_id
    if (length(hits) <= 1) {
        return(NULL)
    }

    g_string <- tryCatch(
        string_db_obj$get_subnetwork(hits),
        error = function(e) {
            warning("STRINGdb network construction failed: ",
                conditionMessage(e),
                ". Skipping PPI analysis.",
                call. = FALSE
            )
            return(NULL)
        }
    )
    if (is.null(g_string) || igraph::vcount(g_string) == 0) {
        return(NULL)
    }
    g_string <- igraph::delete_vertices(g_string, igraph::V(g_string)[igraph::degree(g_string) == 0])
    if (igraph::vcount(g_string) == 0) {
        return(NULL)
    }

    map_df <- targets_mapped[targets_mapped$STRING_id %in% igraph::V(g_string)$name, ]
    map_df <- map_df[!duplicated(map_df$STRING_id), ]
    symbol_map <- setNames(map_df$gene, map_df$STRING_id)
    igraph::V(g_string)$symbol <- symbol_map[igraph::V(g_string)$name]

    lfc_vals <- setNames(as.numeric(global_glist), names(global_glist))[igraph::V(g_string)$symbol]
    lfc_vals[is.na(lfc_vals)] <- 0
    igraph::V(g_string)$lfc <- as.numeric(lfc_vals)
    igraph::V(g_string)$deg <- as.numeric(igraph::degree(g_string))

    if (is.null(igraph::E(g_string)$combined_score)) igraph::E(g_string)$combined_score <- ppi_score
    igraph::E(g_string)$combined_score <- as.numeric(igraph::E(g_string)$combined_score)

    top_n_labels <- 25
    num_nodes <- length(igraph::V(g_string)$deg)
    threshold_deg <- if (num_nodes > top_n_labels) sort(igraph::V(g_string)$deg, decreasing = TRUE)[top_n_labels] else 0
    igraph::V(g_string)$label_text <- ifelse(igraph::V(g_string)$deg >= threshold_deg, igraph::V(g_string)$symbol, NA)

    p_ppi <- ggraph::ggraph(g_string, layout = "fr") +
        ggraph::geom_edge_link(ggplot2::aes(alpha = combined_score), color = "grey60", edge_width = 0.5, show.legend = FALSE) +
        ggraph::geom_node_point(ggplot2::aes(color = lfc, size = deg), stroke = 0.5) +
        ggraph::geom_node_text(ggplot2::aes(label = label_text), repel = TRUE, size = 3.5, max.overlaps = Inf, fontface = "bold", bg.color = "white", bg.r = 0.1) +
        ggplot2::scale_color_distiller(palette = "PuOr", name = "LFC") +
        ggplot2::scale_size_continuous(range = c(2, 8), guide = "none") +
        ggraph::scale_edge_alpha_continuous(range = c(0.4, 0.9)) +
        ggraph::theme_graph(base_family = "sans", background = "white") +
        ggplot2::labs(title = paste0("PPI Network: ", current_proj_name), subtitle = paste0("Interacting Nodes: ", num_nodes, " | Score: ", min(igraph::E(g_string)$combined_score, na.rm = TRUE), "-", max(igraph::E(g_string)$combined_score, na.rm = TRUE))) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))

    return(p_ppi)
}

#' @title Generate Summary GO Lollipop Facet Plot
#' @return A list of \code{ggplot} objects, one per ontology facet.
#' @keywords internal
#' @noRd
plot_summary_go_lollipop <- function(all_go_results, base_project_name) {
    plot_list <- list()
    valid_results <- Filter(function(x) !is.null(x) && is.data.frame(x) && nrow(x) > 0, all_go_results)
    if (length(valid_results) == 0) {
        return(plot_list)
    }

    final_go_df <- do.call(rbind, valid_results)
    if (is.null(final_go_df)) {
        return(plot_list)
    }

    if (!"CleanLoopType" %in% colnames(final_go_df)) final_go_df$CleanLoopType <- final_go_df$LoopType
    use_ggtext <- requireNamespace("ggtext", quietly = TRUE)

    for (ltype in unique(final_go_df$CleanLoopType)) {
        sub_df <- final_go_df %>% dplyr::filter(CleanLoopType == ltype)
        if (nrow(sub_df) == 0) next

        sub_df$logP <- -log10(sub_df$p.adjust)
        max_count <- max(sub_df$Count, na.rm = TRUE)
        if (!is.finite(max_count) || max_count == 0) max_count <- 1
        scale_f <- max(sub_df$logP, na.rm = TRUE) / max_count

        sub_df <- sub_df %>%
            dplyr::group_by(ONTOLOGY) %>%
            dplyr::arrange(logP) %>%
            dplyr::ungroup() %>%
            dplyr::mutate(Description_unique = factor(Description, levels = unique(Description)))

        onto_levels <- sort(unique(sub_df$ONTOLOGY))
        sub_df$ONTOLOGY <- factor(sub_df$ONTOLOGY, levels = onto_levels)
        onto_colors <- RColorBrewer::brewer.pal(max(3, length(onto_levels)), "Dark2")[seq_len(length(onto_levels))]

        if (use_ggtext) {
            sub_df$ONTOLOGY_Plot <- factor(sub_df$ONTOLOGY, levels = onto_levels, labels = paste0("<span style='color:", onto_colors, "'>", onto_levels, "</span>"))
        } else {
            sub_df$ONTOLOGY_Plot <- sub_df$ONTOLOGY
        }

        p_go <- ggplot2::ggplot(sub_df, ggplot2::aes(y = Description_unique)) +
            ggplot2::geom_segment(ggplot2::aes(x = 0, xend = logP, yend = Description_unique, color = ONTOLOGY), linewidth = 3) +
            ggplot2::geom_point(ggplot2::aes(x = logP, color = ONTOLOGY), size = 5) +
            ggplot2::geom_path(ggplot2::aes(x = Count * scale_f, group = 1), color = "grey60", linewidth = 1.5, linetype = "11") +
            ggplot2::geom_point(ggplot2::aes(x = Count * scale_f), color = "grey60", size = 4, shape = 17) +
            ggplot2::scale_x_continuous(name = expression(-log[10](FDR)), expand = ggplot2::expansion(mult = c(0, 0.6)), sec.axis = ggplot2::sec_axis(~ . / scale_f, name = "Gene Counts")) +
            ggplot2::facet_grid(ONTOLOGY_Plot ~ ., scales = "free_y", space = "free_y", switch = "y") +
            ggplot2::labs(y = NULL, title = paste0("GO Enrichment: ", ltype), subtitle = "Colored Dot: Significance | Grey Triangle: Gene Count") +
            ggplot2::theme_bw() +
            ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10, color = "black"), strip.placement = "outside", strip.background = ggplot2::element_blank(), panel.grid.major.y = ggplot2::element_line(color = "grey95", linetype = "dashed"), legend.position = "none") +
            ggplot2::scale_color_brewer(palette = "Dark2")

        if (use_ggtext) {
            p_go <- p_go + ggplot2::theme(strip.text.y.left = ggtext::element_markdown(angle = 0, face = "bold", size = 12))
        } else {
            p_go <- p_go + ggplot2::theme(strip.text.y.left = ggplot2::element_text(angle = 0, face = "bold", size = 12, color = onto_colors))
        }

        plot_list[[ltype]] <- p_go
    }
    return(plot_list)
}

#' Internal: Build expression heatmap sub-plot
#' @keywords internal
#' @noRd
.build_expression_heatmap <- function(target_genes, curr_mat, curr_meta,
                                      heatmap_ntop, current_proj_name,
                                      skip_heatmap) {
    if (skip_heatmap) {
        return(NULL)
    }
    # Case-insensitive matching using canonical uppercase keys.
    # After match, use uppercase keys throughout so that downstream
    # indexing (e.g. mouse Trp53 vs human TP53) is consistent.
    target_key <- toupper(target_genes)
    expr_key <- toupper(rownames(curr_mat))
    expr_idx <- match(target_key, expr_key)
    matched <- !is.na(expr_idx)
    idx <- expr_idx[matched]
    matched_keys <- target_key[matched]
    # Deduplicate: multiple target genes may map to the same uppercase key
    keep <- !duplicated(matched_keys)
    idx <- idx[keep]
    matched_keys <- matched_keys[keep]
    curr_mat <- curr_mat[idx, , drop = FALSE]
    rownames(curr_mat) <- matched_keys
    expr_genes <- matched_keys
    if (length(expr_genes) < 5) {
        return(NULL)
    }
    mat_plot <- log2(curr_mat[expr_genes, , drop = FALSE] + 1)
    if (!is.null(heatmap_ntop) && nrow(mat_plot) > heatmap_ntop) {
        row_vars <- apply(mat_plot, 1, var, na.rm = TRUE)
        mat_plot <- mat_plot[head(names(sort(row_vars, decreasing = TRUE)), heatmap_ntop), , drop = FALSE]
    }
    mat_scaled <- t(scale(t(mat_plot)))
    mat_scaled[mat_scaled > 2] <- 2
    mat_scaled[mat_scaled < -2] <- -2
    mat_scaled[is.na(mat_scaled)] <- 0
    col_fun <- circlize::colorRamp2(c(-2, 0, 2), c("#2BB2D1", "white", "#FF8181"))
    groups <- unique(curr_meta$Group)
    n_groups <- length(groups)
    cols <- if (n_groups <= 8) {
        RColorBrewer::brewer.pal(max(3, n_groups), "Set2")[seq_len(n_groups)]
    } else {
        grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(n_groups)
    }
    names(cols) <- as.character(groups)
    ha <- ComplexHeatmap::HeatmapAnnotation(
        Group = curr_meta$Group,
        col = list(Group = cols),
        simple_anno_size = unit(0.3, "cm")
    )
    ComplexHeatmap::Heatmap(mat_scaled,
        name = "Z-score", col = col_fun,
        cluster_columns = FALSE, show_row_names = (nrow(mat_scaled) <= 80),
        top_annotation = ha, border = TRUE,
        column_title = paste0("Expression Heatmap\n", current_proj_name),
        use_raster = FALSE
    )
}

#' Internal: Add raincloud plots to connectivity plot list
#' @keywords internal
#' @noRd
.add_connectivity_rainclouds <- function(plot_df_rc, custom_colors, plots_list) {
    if (nlevels(plot_df_rc$Conn_Group) <= 1) {
        return(plots_list)
    }
    clean_theme <- ggplot2::theme_classic() + ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 11, color = "black"),
        legend.position = "none",
        axis.text.x = ggplot2::element_text(angle = 20, hjust = 1, size = 10, color = "black"),
        axis.text.y = ggplot2::element_text(size = 10, color = "black"),
        axis.title.y = ggplot2::element_text(size = 12, face = "bold"),
        axis.line = ggplot2::element_line(color = "black", linewidth = 0.6),
        axis.ticks = ggplot2::element_line(color = "black"),
        panel.grid = ggplot2::element_blank()
    )
    get_pval_str <- function(val_col) {
        lvls <- levels(plot_df_rc$Conn_Group)
        if (!"Others" %in% lvls) {
            return("Wilcox P: NA (No 'Others' group)")
        }
        res <- character()
        if ("High Distal" %in% lvls) {
            p <- tryCatch(
                stats::wilcox.test(
                    plot_df_rc[[val_col]][plot_df_rc$Conn_Group == "High Distal"],
                    plot_df_rc[[val_col]][plot_df_rc$Conn_Group == "Others"]
                )$p.value,
                error = function(e) NA_real_
            )
            if (!is.na(p)) {
                res <- c(res, paste0(
                    "Distal=", signif(p, 3), " (",
                    dplyr::case_when(
                        p < 0.001 ~ "***", p < 0.01 ~ "**",
                        p < 0.05 ~ "*", TRUE ~ "ns"
                    ), ")"
                ))
            }
        }
        if ("High Total" %in% lvls) {
            p <- tryCatch(
                stats::wilcox.test(
                    plot_df_rc[[val_col]][plot_df_rc$Conn_Group == "High Total"],
                    plot_df_rc[[val_col]][plot_df_rc$Conn_Group == "Others"]
                )$p.value,
                error = function(e) NA_real_
            )
            if (!is.na(p)) {
                res <- c(res, paste0(
                    "Total=", signif(p, 3), " (",
                    dplyr::case_when(
                        p < 0.001 ~ "***", p < 0.01 ~ "**",
                        p < 0.05 ~ "*", TRUE ~ "ns"
                    ), ")"
                ))
            }
        }
        if (length(res) == 0) {
            return("Wilcox P: NA")
        }
        paste0("Wilcox P (vs Others):\n", paste(res, collapse = " | "))
    }
    base_box <- function(y_var, y_lab, density_adjust = 0.5) {
        ggplot2::ggplot(plot_df_rc, ggplot2::aes(fill = Conn_Group)) +
            ggplot2::geom_jitter(
                ggplot2::aes(
                    x = .data$Conn_Group_jitter,
                    y = .data[[y_var]], color = .data$Conn_Group
                ),
                shape = 16, width = 0.03, height = 0, alpha = 0.6, size = 0.8, stroke = 0
            ) +
            ggplot2::stat_boxplot(
                ggplot2::aes(
                    x = .data$Conn_Group_num,
                    y = .data[[y_var]], color = .data$Conn_Group
                ),
                geom = "errorbar", width = 0.05, linewidth = 0.5
            ) +
            ggplot2::geom_boxplot(
                ggplot2::aes(
                    x = .data$Conn_Group_num,
                    y = .data[[y_var]], color = .data$Conn_Group
                ),
                width = 0.12, notch = TRUE, outlier.shape = NA, alpha = 1, linewidth = 0.5
            ) +
            ggplot2::stat_summary(
                ggplot2::aes(
                    x = .data$Conn_Group_num,
                    y = .data[[y_var]]
                ),
                fun = median, fun.min = median, fun.max = median,
                geom = "crossbar", width = 0.1, color = "black", linewidth = 0.4
            ) +
            {
                if (requireNamespace("ggdist", quietly = TRUE)) {
                    ggdist::stat_slab(
                        ggplot2::aes(
                            x = .data$Conn_Group_slab,
                            y = .data[[y_var]], fill = .data$Conn_Group
                        ),
                        adjust = density_adjust, width = 0.35, justification = 0, normalize = "groups", alpha = 0.3, color = NA
                    )
                }
            } +
            {
                if (requireNamespace("ggdist", quietly = TRUE)) {
                    ggdist::stat_slab(
                        ggplot2::aes(
                            x = .data$Conn_Group_slab,
                            y = .data[[y_var]], color = .data$Conn_Group
                        ),
                        adjust = density_adjust, width = 0.35, justification = 0, normalize = "groups",
                        fill = NA, alpha = 0.5, linewidth = 0.4
                    )
                }
            } +
            ggplot2::scale_x_continuous(
                breaks = seq_along(levels(plot_df_rc$Conn_Group)),
                labels = levels(plot_df_rc$Conn_Group)
            ) +
            ggplot2::scale_fill_manual(values = custom_colors) +
            ggplot2::scale_color_manual(values = custom_colors) +
            ggplot2::labs(
                title = "Regulation: High_connectivity vs Others",
                subtitle = get_pval_str(y_var), x = NULL, y = y_lab
            ) +
            clean_theme
    }
    plots_list$Raincloud_LFC <- base_box("LFC", "Log2 Fold Change (LFC)", density_adjust = 1.5) +
        ggplot2::coord_cartesian(
            xlim = c(0.75, length(levels(plot_df_rc$Conn_Group)) + 0.6),
            ylim = c(-4, 4)
        ) +
        ggplot2::geom_hline(
            yintercept = 0, linetype = "dashed",
            color = "grey45", linewidth = 0.6
        )
    plots_list$Raincloud_Expr <- base_box("Expression", "Log2(Mean Expression + 1)") +
        ggplot2::coord_cartesian(
            xlim = c(0.75, length(levels(plot_df_rc$Conn_Group)) + 0.6)
        )
    plots_list
}

#' @title Generate Expression Heatmap and Connectivity Plots
#' @param target_genes Character vector of target gene symbols to
#'   include in the heatmap and connectivity analysis.
#' @param tpm_mat_raw Numeric matrix of normalised expression values
#'   (genes x samples) with gene symbols as row names.
#' @param meta_raw Data frame with columns \code{SampleID} and
#'   \code{Group} defining sample groups.
#' @param loop_stats_df Data frame. Promoter-centric statistics from
#'   \code{refine_loop_anchors_by_expression} or
#'   \code{annotate_peaks_and_loops}.
#' @param global_glist Named numeric vector (gene-level ranking metric)
#'   used to sort genes for heatmap display.
#' @param heatmap_ntop Integer. Maximum number of top genes to include
#'   in the heatmap.
#' @param cor_method Character. Correlation method for the heatmap
#'   sample annotation (passed to \code{cor}).
#' @param current_proj_name Character. Prefix for plot titles.
#' @param source_type Character. Source label (e.g. \code{"loops"} or
#'   \code{"targets"}) for legend titles.
#' @param target_col Character or \code{NULL}. Column name in
#'   \code{loop_stats_df} to use as connectivity degree. \code{NULL}
#'   auto-detects \code{Total_Loops}, \code{Loop_Degree}, or
#'   \code{degree}.
#' @param skip_heatmap Logical. If \code{TRUE}, skip the ComplexHeatmap
#'   rendering and only produce connectivity plots. Default \code{FALSE}.
#' @return A named list of plot objects (Heatmap, Scatter, Raincloud_LFC, Raincloud_Expr).
#' @keywords internal
#' @noRd
run_heatmap_and_connectivity <- function(target_genes, tpm_mat_raw, meta_raw, loop_stats_df, global_glist, heatmap_ntop, cor_method, current_proj_name, source_type, target_col = NULL, skip_heatmap = FALSE) {
    plots_list <- list()
    if (!skip_heatmap && (!requireNamespace("ComplexHeatmap", quietly = TRUE) ||
        !requireNamespace("circlize", quietly = TRUE))) {
        warning("ComplexHeatmap/circlize not installed; skipping heatmap.", call. = FALSE)
        skip_heatmap <- TRUE
    }
    colnames(tpm_mat_raw) <- trimws(colnames(tpm_mat_raw))
    valid_s <- intersect(meta_raw$SampleID, colnames(tpm_mat_raw))
    if (length(valid_s) == 0) {
        return(plots_list)
    }

    curr_mat <- tpm_mat_raw[, valid_s, drop = FALSE]
    curr_meta <- meta_raw %>% dplyr::filter(SampleID %in% valid_s)

    plots_list$Heatmap <- .build_expression_heatmap(
        target_genes, curr_mat, curr_meta, heatmap_ntop, current_proj_name, skip_heatmap
    )

    if (is.null(loop_stats_df)) {
        return(plots_list)
    }
    use_col <- NULL
    display_tag <- "Total Loops"

    if (!is.null(target_col)) {
        if (target_col %in% colnames(loop_stats_df)) {
            use_col <- target_col
            display_tag <- target_col
        } else {
            return(plots_list)
        }
    } else {
        if ("Total_Loops" %in% colnames(loop_stats_df)) {
            use_col <- "Total_Loops"
        } else if ("Loop_Degree" %in% colnames(loop_stats_df)) {
            use_col <- "Loop_Degree"
        } else if ("degree" %in% colnames(loop_stats_df)) use_col <- "degree"
    }

    if (is.null(use_col)) {
        return(plots_list)
    }

    gene_col_name <- colnames(loop_stats_df)[1]
    valid_targets <- intersect(target_genes, loop_stats_df[[gene_col_name]])
    if (length(valid_targets) < 5) {
        return(plots_list)
    }

    cols_to_extract <- unique(c(gene_col_name, use_col, intersect(colnames(loop_stats_df), c("Is_High_Connectivity_Gene", "Is_High_Distal_Connectivity_Gene", "High_Connectivity_Gene"))))
    stats_subset <- loop_stats_df[loop_stats_df[[gene_col_name]] %in% valid_targets, cols_to_extract]
    colnames(stats_subset)[which(colnames(stats_subset) == use_col)] <- "Degree"
    colnames(stats_subset)[1] <- "Gene"

    valid_expr_targets <- intersect(stats_subset$Gene, rownames(curr_mat))
    if (length(valid_expr_targets) < 5) {
        return(plots_list)
    }

    stats_subset <- stats_subset[stats_subset$Gene %in% valid_expr_targets, ]
    plot_df <- stats_subset %>%
        dplyr::mutate(Expression = as.numeric(log2(rowMeans(curr_mat[stats_subset$Gene, , drop = FALSE], na.rm = TRUE) + 1)), LFC = as.numeric(global_glist[stats_subset$Gene]), Log10Degree = log10(Degree)) %>%
        dplyr::filter(!is.na(Expression), !is.na(LFC), Degree >= 1)

    if (nrow(plot_df) < 5) {
        return(plots_list)
    }

    full_title_suffix <- paste0(if (source_type == "loops") "Looped (Specified Types)" else "Looped Targets", " | ", display_tag)

    plots_list$Scatter <- if (requireNamespace("ggpointdensity", quietly = TRUE) &&
        requireNamespace("viridis", quietly = TRUE) &&
        requireNamespace("ggpubr", quietly = TRUE)) {
        ggplot2::ggplot(plot_df, ggplot2::aes(x = Log10Degree, y = Expression)) +
            ggpointdensity::geom_pointdensity(alpha = 0.6, size = 1.5) +
            viridis::scale_color_viridis(option = "D", name = "Density") +
            ggplot2::geom_smooth(method = "lm", formula = y ~ x, color = "black", se = TRUE, linewidth = 0.8) +
            ggpubr::stat_cor(method = cor_method, label.x.npc = "left", label.y.npc = "top", size = 4) +
            ggplot2::theme_bw() +
            ggplot2::theme(panel.grid = ggplot2::element_blank(), panel.background = ggplot2::element_rect(color = "black", fill = "transparent"), legend.key = ggplot2::element_rect(fill = "transparent"), plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")) +
            ggplot2::labs(title = "Connectivity vs Expression (Scatter)", subtitle = paste0(full_title_suffix, "\nGenes: ", nrow(plot_df)), x = paste0("Log10 (", display_tag, ")"), y = "Log2(Mean Expression + 1)")
    } else {
        NULL
    }

    if ("Is_High_Distal_Connectivity_Gene" %in% colnames(plot_df) && "Is_High_Connectivity_Gene" %in% colnames(plot_df)) {
        plot_df_rc <- plot_df %>%
            dplyr::mutate(
                high_distal = Is_High_Distal_Connectivity_Gene %in% c("Yes", "TRUE", TRUE, 1),
                high_total = Is_High_Connectivity_Gene %in% c("Yes", "TRUE", TRUE, 1)
            )
        plot_df_rc$Conn_Group <- ifelse(plot_df_rc$high_distal, "High Distal", ifelse(plot_df_rc$high_total, "High Total", "Others"))
        plot_df_rc$high_distal <- NULL
        plot_df_rc$high_total <- NULL
        plot_df_rc$Conn_Group <- factor(plot_df_rc$Conn_Group, levels = c("High Distal", "High Total", "Others"))
        custom_colors <- c("High Distal" = "#9BC985", "High Total" = "#ECB884", "Others" = "#82969D")
    } else {
        deg_thresh <- max(quantile(plot_df$Degree, 0.75, na.rm = TRUE), 2)
        plot_df_rc <- plot_df %>% dplyr::mutate(Conn_Group = factor(ifelse(Degree >= deg_thresh, "High Total", "Others"), levels = c("High Total", "Others")))
        custom_colors <- c("High Total" = "#ECB884", "Others" = "#82969D")
    }

    plot_df_rc <- plot_df_rc %>%
        dplyr::filter(!is.na(Conn_Group)) %>%
        droplevels() %>%
        dplyr::mutate(
            Conn_Group_num = as.numeric(Conn_Group)
        ) %>%
        dplyr::mutate(
            Conn_Group_jitter = .data$Conn_Group_num - 0.12,
            Conn_Group_slab = .data$Conn_Group_num + 0.07
        )

    plots_list <- .add_connectivity_rainclouds(
        plot_df_rc, custom_colors, plots_list
    )
    return(plots_list)
}

.anchor_matches_targets <- function(gene_string, target_genes) {
    genes <- clean_gene_names(gene_string, ";")
    length(genes) > 0 && any(genes %in% target_genes)
}

.subset_motif_loop_df <- function(loop_df, src, task_name) {
    if (!is.data.frame(loop_df) || !identical(src, "loops") || !"loop_type" %in% colnames(loop_df)) {
        return(loop_df)
    }

    loop_types <- unique(as.character(loop_df$loop_type))
    task_map <- paste0(gsub("-", "", loop_types, fixed = TRUE), "_Genes")
    matched_types <- loop_types[task_map == task_name]
    if (length(matched_types) == 0) {
        return(loop_df)
    }
    loop_df[loop_df$loop_type %in% matched_types, , drop = FALSE]
}

.is_promoter_anchor_type <- function(anchor_type) {
    anchor_type <- trimws(as.character(anchor_type))
    !is.na(anchor_type) & .is_promoter_like(anchor_type)
}

.is_enhancer_like_anchor_type <- function(anchor_type) {
    anchor_type <- trimws(as.character(anchor_type))
    !is.na(anchor_type) & .is_distal_like(anchor_type)
}

.empty_anchor_df <- function() {
    data.frame(
        anchor_id = character(),
        chr = character(),
        start = integer(),
        end = integer(),
        anchor_type = character(),
        cluster_id = character(),
        stringsAsFactors = FALSE
    )
}

.deduplicate_anchor_df <- function(anchor_df) {
    if (is.null(anchor_df) || nrow(anchor_df) == 0) {
        return(.empty_anchor_df())
    }

    anchor_df <- anchor_df[!is.na(anchor_df$chr) & nzchar(anchor_df$chr), , drop = FALSE]
    anchor_df <- anchor_df[!is.na(anchor_df$start) & !is.na(anchor_df$end), , drop = FALSE]
    if (nrow(anchor_df) == 0) {
        return(.empty_anchor_df())
    }

    anchor_df$anchor_id <- ifelse(
        is.na(anchor_df$anchor_id) | !nzchar(anchor_df$anchor_id),
        paste(anchor_df$chr, anchor_df$start, anchor_df$end, sep = "_"),
        as.character(anchor_df$anchor_id)
    )
    anchor_df <- anchor_df[!duplicated(anchor_df$anchor_id), , drop = FALSE]
    rownames(anchor_df) <- NULL
    anchor_df
}

.anchor_df_to_gr <- function(anchor_df) {
    anchor_df <- .deduplicate_anchor_df(anchor_df)
    if (nrow(anchor_df) == 0) {
        return(.with_known_upstream_noise_suppressed(GenomicRanges::GRanges()))
    }

    gr <- .with_known_upstream_noise_suppressed(GenomicRanges::GRanges(
        seqnames = anchor_df$chr,
        ranges = IRanges::IRanges(start = anchor_df$start, end = anchor_df$end)
    ))
    names(gr) <- anchor_df$anchor_id
    S4Vectors::mcols(gr)$anchor_id <- anchor_df$anchor_id
    S4Vectors::mcols(gr)$anchor_type <- anchor_df$anchor_type
    if ("cluster_id" %in% colnames(anchor_df)) {
        S4Vectors::mcols(gr)$cluster_id <- anchor_df$cluster_id
    }
    gr
}

.make_anchor_df <- function(loop_df, idx, side, anchor_types) {
    if (length(idx) == 0) {
        return(.empty_anchor_df())
    }

    chr_col <- paste0("chr", side)
    start_col <- paste0("start", side)
    end_col <- paste0("end", side)
    id_col <- paste0("a", side, "_id")
    anchor_id <- if (id_col %in% colnames(loop_df)) {
        loop_df[[id_col]][idx]
    } else {
        paste(loop_df[[chr_col]][idx], loop_df[[start_col]][idx], loop_df[[end_col]][idx], sep = "_")
    }
    cluster_id <- if ("cluster_id" %in% colnames(loop_df)) {
        as.character(loop_df$cluster_id[idx])
    } else {
        rep(NA_character_, length(idx)) # no real component -- do not fabricate
    }

    data.frame(
        anchor_id = as.character(anchor_id),
        chr = as.character(loop_df[[chr_col]][idx]),
        start = as.integer(loop_df[[start_col]][idx]),
        end = as.integer(loop_df[[end_col]][idx]),
        anchor_type = as.character(anchor_types[idx]),
        cluster_id = cluster_id,
        stringsAsFactors = FALSE
    )
}

.prepare_motif_anchor_sets <- function(
  loop_df, target_genes,
  anchor_registry = NULL
) {
    empty_gr <- .with_known_upstream_noise_suppressed(GenomicRanges::GRanges())
    empty_sets <- list(
        target_loop_n = 0L,
        proximal_fg = empty_gr,
        distal_fg = empty_gr,
        proximal_bg = empty_gr,
        distal_bg = empty_gr
    )
    if (!is.data.frame(loop_df) || nrow(loop_df) == 0) {
        return(empty_sets)
    }

    col_g1 <- intersect(c("anchor1_gene", "Anchor1_Gene", "gene_name_1", "Gene_Name_1", "Symbol_1", "nearest_gene_1", "gene1"), colnames(loop_df))[1]
    col_g2 <- intersect(c("anchor2_gene", "Anchor2_Gene", "gene_name_2", "Gene_Name_2", "Symbol_2", "nearest_gene_2", "gene2"), colnames(loop_df))[1]
    col_t1 <- intersect(c("anchor1_type", "Anchor1_Type", "type1"), colnames(loop_df))[1]
    col_t2 <- intersect(c("anchor2_type", "Anchor2_Type", "type2"), colnames(loop_df))[1]
    required_cols <- c("chr1", "start1", "end1", "chr2", "start2", "end2")
    if (any(is.na(c(col_g1, col_g2, col_t1, col_t2))) || !all(required_cols %in% colnames(loop_df))) {
        return(empty_sets)
    }

    target_genes <- clean_gene_names(target_genes)
    if (length(target_genes) == 0) {
        return(empty_sets)
    }

    a1_hits <- vapply(loop_df[[col_g1]], .anchor_matches_targets,
        target_genes = target_genes, FUN.VALUE = logical(1)
    )
    a2_hits <- vapply(loop_df[[col_g2]], .anchor_matches_targets,
        target_genes = target_genes, FUN.VALUE = logical(1)
    )
    t1 <- trimws(as.character(loop_df[[col_t1]]))
    t2 <- trimws(as.character(loop_df[[col_t2]]))

    a1_target_promoter <- a1_hits & .is_promoter_anchor_type(t1)
    a2_target_promoter <- a2_hits & .is_promoter_anchor_type(t2)
    target_loop_idx <- which(a1_target_promoter | a2_target_promoter)
    is_bg_loop <- !(seq_len(nrow(loop_df)) %in% target_loop_idx)

    proximal_fg_df <- rbind(
        .make_anchor_df(loop_df, which(a1_target_promoter), "1", t1),
        .make_anchor_df(loop_df, which(a2_target_promoter), "2", t2)
    )
    distal_fg_df <- rbind(
        .make_anchor_df(loop_df, which(a1_target_promoter & .is_enhancer_like_anchor_type(t2)), "2", t2),
        .make_anchor_df(loop_df, which(a2_target_promoter & .is_enhancer_like_anchor_type(t1)), "1", t1)
    )
    proximal_bg_df <- rbind(
        .make_anchor_df(loop_df, which(is_bg_loop & .is_promoter_anchor_type(t1)), "1", t1),
        .make_anchor_df(loop_df, which(is_bg_loop & .is_promoter_anchor_type(t2)), "2", t2)
    )
    distal_bg_df <- rbind(
        .make_anchor_df(loop_df, which(is_bg_loop & .is_enhancer_like_anchor_type(t1)), "1", t1),
        .make_anchor_df(loop_df, which(is_bg_loop & .is_enhancer_like_anchor_type(t2)), "2", t2)
    )

    # Exclude anchors that appear in both foreground and background.
    # A single anchor can participate in multiple loops -- a target loop
    # and a non-target loop -- and would otherwise violate the independence
    # assumption of Fisher's test and permutation.
    fg_ids <- unique(c(proximal_fg_df$anchor_id, distal_fg_df$anchor_id))
    n_bg_before <- nrow(proximal_bg_df) + nrow(distal_bg_df)
    proximal_bg_df <- proximal_bg_df[!proximal_bg_df$anchor_id %in% fg_ids, , drop = FALSE]
    distal_bg_df <- distal_bg_df[!distal_bg_df$anchor_id %in% fg_ids, , drop = FALSE]
    n_bg_after <- nrow(proximal_bg_df) + nrow(distal_bg_df)
    n_overlap <- n_bg_before - n_bg_after
    if (n_overlap > 0) {
        message(sprintf(
            "Motif: %d anchor(s) present in both FG and BG excluded from background (%d remaining).",
            n_overlap, n_bg_after
        ))
    }

    # When a canonical anchor registry is available, remap all anchor
    # coordinates to the merged canonical positions.  This ensures that
    # motif results at anchor_merge_gap > 0 are consistent with chromatin
    # refinement and downstream stats.
    sets <- list(
        proximal_fg = .anchor_df_to_gr(proximal_fg_df),
        distal_fg = .anchor_df_to_gr(distal_fg_df),
        proximal_bg = .anchor_df_to_gr(proximal_bg_df),
        distal_bg = .anchor_df_to_gr(distal_bg_df)
    )
    if (!is.null(anchor_registry) && length(anchor_registry) > 0) {
        reg_lookup <- GenomicRanges::GRanges(
            seqnames = GenomicRanges::seqnames(anchor_registry),
            ranges   = GenomicRanges::ranges(anchor_registry)
        )
        names(reg_lookup) <- names(anchor_registry)
        for (nm in names(sets)) {
            gr <- sets[[nm]]
            if (length(gr) == 0) next
            ids <- names(gr)
            idx <- match(ids, names(reg_lookup))
            n_missing <- sum(is.na(idx))
            if (n_missing > 0) {
                stop(
                    "Motif anchor registry mismatch: ", n_missing, " of ",
                    length(ids), " ", nm, " anchor(s) were not found in the ",
                    "canonical anchor registry. This may occur when the ",
                    "annotation object was not produced by the current ",
                    "pipeline version. Please re-run annotate_peaks_and_loops().",
                    call. = FALSE
                )
            }
            canonical_gr <- reg_lookup[idx]
            names(canonical_gr) <- ids
            # Preserve metadata from original GRanges
            for (mc in names(S4Vectors::mcols(gr))) {
                S4Vectors::mcols(canonical_gr)[[mc]] <-
                    S4Vectors::mcols(gr)[[mc]]
            }
            sets[[nm]] <- canonical_gr
        }
    }

    c(list(target_loop_n = length(target_loop_idx)), sets)
}

.calc_gc_fraction <- function(seq_set) {
    seq_chr <- as.character(seq_set)
    seq_len <- nchar(seq_chr)
    gc_len <- nchar(gsub("[^GCgc]", "", seq_chr))
    out <- rep(NA_real_, length(seq_chr))
    keep <- seq_len > 0
    out[keep] <- gc_len[keep] / seq_len[keep]
    out
}

.sample_gc_matched_background <- function(fg_gr, bg_gr, genome_obj, max_bg = 2000L, gc_bins = 5L) {
    if (length(bg_gr) == 0) {
        return(bg_gr)
    }

    target_n <- min(length(bg_gr), as.integer(max_bg))
    if (target_n <= 0L) {
        return(bg_gr[0])
    }
    if (length(bg_gr) <= target_n || length(fg_gr) == 0) {
        return(bg_gr)
    }

    fg_gc <- .calc_gc_fraction(BSgenome::getSeq(genome_obj, fg_gr))
    bg_gc <- .calc_gc_fraction(BSgenome::getSeq(genome_obj, bg_gr))

    # Detect foreground anchors with non-finite GC (zero-length sequences
    # from chromosome-boundary trimming or gap regions). These cannot be
    # GC-matched; they are included in the final set but excluded from the
    # GC bin proportion calculation.  Warn if they are a material fraction.
    n_na_fg <- sum(!is.finite(fg_gc))
    if (n_na_fg > 0) {
        na_frac <- n_na_fg / length(fg_gc)
        if (na_frac >= 0.01) {
            warning(
                n_na_fg, " of ", length(fg_gc), " foreground anchors (",
                round(na_frac * 100, 1), "%) have non-finite GC content ",
                "and cannot be GC-matched. GC bin proportions are derived ",
                "from the remaining ", length(fg_gc) - n_na_fg,
                " finite-GC anchors. If this fraction is high, inspect your ",
                "anchor regions for zero-width sequences (e.g. at chromosome ",
                "boundaries or assembly gaps).",
                call. = FALSE
            )
        }
    }

    finite_gc <- c(fg_gc[is.finite(fg_gc)], bg_gc[is.finite(bg_gc)])
    if (length(finite_gc) < 2 || length(unique(finite_gc)) < 2) {
        return(bg_gr[sample(seq_along(bg_gr), target_n)])
    }

    gc_breaks <- unique(stats::quantile(
        finite_gc,
        probs = seq(0, 1, length.out = gc_bins + 1),
        na.rm = TRUE,
        names = FALSE,
        type = 8
    ))
    if (length(gc_breaks) < 2) {
        return(bg_gr[sample(seq_along(bg_gr), target_n)])
    }

    fg_bin <- cut(fg_gc, breaks = gc_breaks, include.lowest = TRUE, labels = FALSE)
    bg_bin <- cut(bg_gc, breaks = gc_breaks, include.lowest = TRUE, labels = FALSE)
    fg_tab <- table(fg_bin)
    if (length(fg_tab) == 0 || sum(fg_tab) == 0) {
        return(bg_gr[sample(seq_along(bg_gr), target_n)])
    }

    fg_prop <- as.numeric(fg_tab) / sum(fg_tab)
    desired <- floor(target_n * fg_prop)
    remainder <- target_n - sum(desired)
    if (remainder > 0) {
        bump_idx <- rep(seq_along(desired), length.out = remainder)
        desired[bump_idx] <- desired[bump_idx] + 1L
    }

    bin_ids <- as.integer(names(fg_tab))
    selected <- integer()
    n_matched <- 0L # track GC-matched count
    for (i in seq_along(bin_ids)) {
        candidates <- which(bg_bin == bin_ids[i])
        take_n <- min(length(candidates), desired[i])
        if (take_n > 0) {
            selected <- c(selected, sample(candidates, take_n))
            n_matched <- n_matched + take_n
        }
    }
    selected <- unique(selected)
    n_fallback <- 0L
    if (length(selected) < target_n) {
        remaining <- setdiff(seq_along(bg_gr), selected)
        if (length(remaining) > 0) {
            n_fallback <- min(length(remaining), target_n - length(selected))
            selected <- c(selected, sample(remaining, n_fallback))
        }
    }

    # Warn if a substantial fraction of background couldn't be GC-matched
    matched_frac <- n_matched / max(length(selected), 1L)
    if (matched_frac < 0.9) {
        warning(
            "Only ", round(matched_frac * 100), "% of background sequences ",
            "were GC-matched to foreground (", n_matched, "/", length(selected),
            "). The remaining ", n_fallback,
            " were randomly sampled, which may bias motif enrichment. ",
            "Consider increasing the background pool or reducing gc_bins.",
            call. = FALSE
        )
    }

    bg_gr[sort(unique(selected))]
}

.empty_motif_output <- function() {
    list(results = list(proximal = NULL, distal = NULL), plots = list())
}
#' @title Run Dual Motif Analysis for Loop Anchors
#' @param target_genes Character vector of target gene symbols to anchor
#'   the foreground set.
#' @param loop_df Data frame. Loop annotation from
#'   \code{refine_loop_anchors_by_expression} or
#'   \code{annotate_peaks_and_loops}.
#' @param genome_id Character. Genome assembly ID (e.g. \code{"hg38"},
#'   \code{"mm10"}).
#' @param pval_thresh Numeric. P-value cutoff for
#'   \code{motifmatchr::matchMotifs}. Default \code{1e-4}.
#' @param current_proj_name Character. Prefix for plot titles and file names.
#' @param top_n Integer. Number of top enriched motifs to display as
#'   sequence logos. Default \code{5}.
#' @param jaspar_db A JASPAR database object (e.g., \code{JASPAR2020::JASPAR2020} or \code{JASPAR2024::JASPAR2024}). Default: \code{NULL} (auto-resolves to \code{JASPAR2020::JASPAR2020} if installed).
#' @param jaspar_collection Character. JASPAR collection to query (e.g., \code{"CORE"}, \code{"CNE"}). Default: \code{"CORE"}.
#' @param motif_max_bg Integer passed to \code{\link{.sample_gc_matched_background}}. Default \code{2000L}.
#' @param motif_gc_bins Integer. Number of GC-content bins for background matching. Default \code{5L}. Increase for regions with highly skewed GC content (e.g., CpG islands).
#' @return A named list containing motif enrichment results and plot objects.
#' @keywords internal
#' @noRd

run_distal_motif_analysis <- function(
  target_genes, loop_df, genome_id, pval_thresh,
  current_proj_name, top_n = 5, jaspar_db = NULL,
  jaspar_collection = "CORE", motif_max_bg = 2000L, motif_gc_bins = 5L,
  motif_n_perm = 0L,
  anchor_registry = NULL
) {
    if (!requireNamespace("motifmatchr", quietly = TRUE) ||
        !requireNamespace("TFBSTools", quietly = TRUE)) {
        warning("Packages 'motifmatchr' and 'TFBSTools' are required for motif analysis. Skipping.", call. = FALSE)
        return(.empty_motif_output())
    }
    if (is.null(jaspar_db)) {
        if (!requireNamespace("JASPAR2020", quietly = TRUE)) {
            warning("Package 'JASPAR2020' is required for motif analysis. Skipping.", call. = FALSE)
            return(.empty_motif_output())
        }
        jaspar_db <- JASPAR2020::JASPAR2020
    }
    bs_pkg <- species_bsgenome_pkg(genome_id)
    if (is.null(bs_pkg)) stop("Unsupported genome: ", genome_id)
    species_id <- if (grepl("mm", genome_id)) 10090 else if (grepl("hg", genome_id)) 9606 else NULL
    genome_obj <- get0(bs_pkg, envir = asNamespace(bs_pkg))
    if (!is.data.frame(loop_df)) loop_df <- as.data.frame(loop_df)

    motif_sets <- .prepare_motif_anchor_sets(loop_df, target_genes,
        anchor_registry = anchor_registry
    )
    has_proximal <- length(motif_sets$proximal_fg) >= 5 && length(motif_sets$proximal_bg) >= 5
    has_distal <- length(motif_sets$distal_fg) >= 5 && length(motif_sets$distal_bg) >= 5
    if (!has_proximal && !has_distal) {
        return(.empty_motif_output())
    }

    plots_list <- list()
    result_tables <- list(proximal = NULL, distal = NULL)
    if (has_proximal) {
        enrich_prox <- .calc_motif_enrichment(
            motif_sets$proximal_fg, motif_sets$proximal_bg,
            genome_obj, pval_thresh, species_id, jaspar_db, jaspar_collection,
            max_bg = motif_max_bg, gc_bins = motif_gc_bins, n_perm = motif_n_perm
        )
        res_prox <- .annotate_motif_families(enrich_prox, jaspar_db, jaspar_collection)
        result_tables$proximal <- res_prox
        plots_list$Proximal_Motif_Bar <- .plot_save_motif(res_prox, paste0(current_proj_name, "_Motif_Proximal"))
        plots_list$Proximal_Motif_Logos <- .plot_top_motif_logos(res_prox, top_n, jaspar_db)
        plots_list$Proximal_Motif_Rank <- .plot_motif_rank_scatter(res_prox, paste0(current_proj_name, "_Motif_Proximal"))
    }

    if (has_distal) {
        enrich_dist <- .calc_motif_enrichment(
            motif_sets$distal_fg, motif_sets$distal_bg,
            genome_obj, pval_thresh, species_id, jaspar_db, jaspar_collection,
            max_bg = motif_max_bg, gc_bins = motif_gc_bins, n_perm = motif_n_perm
        )
        res_dist <- .annotate_motif_families(enrich_dist, jaspar_db, jaspar_collection)
        result_tables$distal <- res_dist
        plots_list$Distal_Motif_Bar <- .plot_save_motif(res_dist, paste0(current_proj_name, "_Motif_Distal"))
        plots_list$Distal_Motif_Logos <- .plot_top_motif_logos(res_dist, top_n, jaspar_db)
        plots_list$Distal_Motif_Rank <- .plot_motif_rank_scatter(res_dist, paste0(current_proj_name, "_Motif_Distal"))
    }

    return(list(results = result_tables, plots = plots_list))
}

#' @title Plot Motif Rank Scatter
#' @return A \code{ggplot} object, or \code{NULL} if input is empty.
#' @keywords internal
#' @noRd
.plot_motif_rank_scatter <- function(res_df, prefix, fdr_thresh = 0.05) {
    if (is.null(res_df) || nrow(res_df) == 0) {
        return(NULL)
    }
    if (!"Family" %in% colnames(res_df)) res_df$Family <- "Unknown"

    plot_df <- res_df %>%
        dplyr::mutate(FDR = ifelse(is.na(FDR), 1, FDR), LogFDR = -log10(FDR + 1e-300), Is_Sig = FDR < fdr_thresh, OddsRatio = pmax(0.8, pmin(1.6, as.numeric(OddsRatio)))) %>%
        dplyr::arrange(dplyr::desc(LogFDR)) %>%
        dplyr::mutate(Rank = dplyr::row_number())

    if (sum(plot_df$Is_Sig, na.rm = TRUE) == 0) {
        plot_df$PlotFamily <- "Not Significant"
    } else {
        fam_counts <- table(plot_df$Family[plot_df$Is_Sig & !plot_df$Family %in% c("Unknown", "", NA)])
        plot_df$PlotFamily <- ifelse(plot_df$Is_Sig, ifelse(plot_df$Family %in% names(sort(fam_counts, decreasing = TRUE))[seq_len(min(10, length(fam_counts)))], plot_df$Family, "Others"), "Not Significant")
    }

    plot_df$PlotFamily <- factor(plot_df$PlotFamily, levels = unique(c(setdiff(unique(plot_df$PlotFamily), c("Others", "Not Significant")), "Others", "Not Significant")))
    color_map <- setNames(c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A")[seq_len(length(setdiff(levels(plot_df$PlotFamily), c("Others", "Not Significant"))))], setdiff(levels(plot_df$PlotFamily), c("Others", "Not Significant")))
    color_map["Others"] <- "black"
    color_map["Not Significant"] <- "grey85"

    return(ggplot2::ggplot(plot_df, ggplot2::aes(x = Rank, y = LogFDR, size = OddsRatio, color = PlotFamily)) +
        ggplot2::geom_point(alpha = 0.8) +
        ggplot2::scale_color_manual(values = color_map, name = "TF Family (Top 10)") +
        ggplot2::scale_radius(name = "Odds Ratio", range = c(1, 7), breaks = c(0.8, 1.0, 1.2, 1.4, 1.6)) +
        ggplot2::geom_hline(yintercept = -log10(fdr_thresh), linetype = "dashed", color = "black", alpha = 0.5) +
        ggplot2::labs(title = paste0("Motif Enrichment Rank: ", basename(prefix)), x = "Rank", y = "-log10(FDR)") +
        ggplot2::theme_classic() +
        ggplot2::theme(legend.position = "right", legend.text = ggplot2::element_text(size = 8), legend.title = ggplot2::element_text(size = 9, face = "bold")))
}

#' @title Calculate Motif Enrichment via Fisher's Exact Test
#' @return A data frame of enrichment results, or \code{NULL} if input is empty.
#' @keywords internal
#' @noRd
.calc_motif_enrichment <- function(
  fg_gr, bg_gr, genome_obj, pval_thresh, species_id,
  jaspar_db = NULL, jaspar_collection = "CORE",
  max_bg = 2000L, gc_bins = 5L, n_perm = 0L
) {
    if (is.null(jaspar_db)) {
        if (!requireNamespace("JASPAR2020", quietly = TRUE)) {
            warning("Package 'JASPAR2020' is required for motif enrichment. Skipping.", call. = FALSE)
            return(NULL)
        }
        jaspar_db <- JASPAR2020::JASPAR2020
    }
    fg_gr <- GenomicRanges::resize(fg_gr[GenomicRanges::start(fg_gr) > 0], width = 500, fix = "center")
    bg_gr <- GenomicRanges::resize(bg_gr[GenomicRanges::start(bg_gr) > 0], width = 500, fix = "center")
    sl <- GenomeInfoDb::seqinfo(genome_obj)
    GenomeInfoDb::seqinfo(fg_gr) <- sl[GenomeInfoDb::seqlevels(fg_gr)]
    fg_gr <- GenomicRanges::trim(fg_gr)
    GenomeInfoDb::seqinfo(bg_gr) <- sl[GenomeInfoDb::seqlevels(bg_gr)]
    bg_gr <- GenomicRanges::trim(bg_gr)
    bg_gr <- .sample_gc_matched_background(fg_gr, bg_gr, genome_obj, max_bg = max_bg, gc_bins = gc_bins)
    if (length(fg_gr) == 0 || length(bg_gr) == 0) {
        message("Motif enrichment skipped: no valid foreground or background regions after trimming.")
        return(NULL)
    }
    n_fg <- length(fg_gr)
    n_bg <- length(bg_gr)

    fg_seq <- BSgenome::getSeq(genome_obj, fg_gr)
    bg_seq <- BSgenome::getSeq(genome_obj, bg_gr)
    pfm_list <- if (!is.null(species_id)) {
        TFBSTools::getMatrixSet(jaspar_db, list(species = species_id, collection = jaspar_collection))
    } else {
        TFBSTools::getMatrixSet(jaspar_db, list(collection = jaspar_collection))
    }
    if (length(pfm_list) == 0 && !is.null(species_id)) {
        pfm_list <- TFBSTools::getMatrixSet(jaspar_db, list(collection = jaspar_collection))
    }
    if (length(pfm_list) == 0) {
        message("Motif enrichment skipped: no motifs found in JASPAR database.")
        return(NULL)
    }

    # --- Motif scanning (run once, shared by observed + all permutations) ---
    fg_counts <- colSums(as.matrix(motifmatchr::motifMatches(motifmatchr::matchMotifs(pfm_list, fg_seq, out = "matches", p.cutoff = pval_thresh))))
    bg_counts <- colSums(as.matrix(motifmatchr::motifMatches(motifmatchr::matchMotifs(pfm_list, bg_seq, out = "matches", p.cutoff = pval_thresh))))
    all_motif_ids <- union(names(fg_counts), names(bg_counts))

    # --- Observed enrichment (Fisher test) ---
    .motif_log_or <- function(a, b, n_fg, n_bg) {
        log(((a + 0.5) * (n_bg - b + 0.5)) /
            ((n_fg - a + 0.5) * (b + 0.5)))
    }
    .build_result_row <- function(m, fg_hits, bg_hits, pval, log_or) {
        if (is.na(pval)) pval <- 1
        if (is.na(log_or)) log_or <- 0
        data.frame(
            MotifID = m,
            MotifName = if (m %in% names(pfm_list)) TFBSTools::name(pfm_list[[m]]) else m,
            Pvalue = pval,
            OddsRatio = exp(log_or),
            LogOddsRatio = log_or,
            FG_Hits = fg_hits, FG_Total = n_fg,
            BG_Hits = bg_hits, BG_Total = n_bg,
            stringsAsFactors = FALSE
        )
    }
    observed <- lapply(all_motif_ids, function(m) {
        a <- if (m %in% names(fg_counts)) fg_counts[[m]] else 0
        b <- if (m %in% names(bg_counts)) bg_counts[[m]] else 0
        ft <- if (a > 0) fisher.test(matrix(c(a, b, n_fg - a, n_bg - b), nrow = 2), alternative = "greater") else NULL
        pval <- if (!is.null(ft)) ft$p.value else 1
        .build_result_row(m, a, b, pval, .motif_log_or(a, b, n_fg, n_bg))
    })
    res_df <- do.call(rbind, observed)
    # Fixed schema: initialize all provenance columns
    res_df$Fisher_Pvalue <- res_df$Pvalue
    res_df$Empirical_Pvalue <- NA_real_
    res_df$Pvalue_Used <- res_df$Fisher_Pvalue
    res_df$Pvalue_Method <- "fisher_exact"
    res_df$Permutation_N <- 0L
    res_df$N_Components <- NA_integer_
    res_df$N_Mixed_Components <- NA_integer_

    # --- Component-block permutation ---
    # Anchors within the same loop component share connectivity structure;
    # treating them as independent inflates significance.  We shuffle fg/bg
    # labels within each component, keeping the component structure intact.
    has_clusters <- FALSE
    permutation_ran <- FALSE
    if (n_perm > 0L) {
        fg_cluster <- S4Vectors::mcols(fg_gr)$cluster_id
        bg_cluster <- S4Vectors::mcols(bg_gr)$cluster_id
        has_clusters <- !is.null(fg_cluster) && !is.null(bg_cluster) &&
            length(c(fg_cluster, bg_cluster)) == n_fg + n_bg &&
            all(!is.na(fg_cluster) & nzchar(as.character(fg_cluster))) &&
            all(!is.na(bg_cluster) & nzchar(as.character(bg_cluster)))
        if (has_clusters) {
            fg_labels <- rep(TRUE, n_fg)
            bg_labels <- rep(FALSE, n_bg)
            all_labels <- c(fg_labels, bg_labels)
            all_clusters <- c(fg_cluster, bg_cluster)
            fg_mat <- as.matrix(motifmatchr::motifMatches(
                motifmatchr::matchMotifs(pfm_list, fg_seq, out = "matches", p.cutoff = pval_thresh)
            ))
            bg_mat <- as.matrix(motifmatchr::motifMatches(
                motifmatchr::matchMotifs(pfm_list, bg_seq, out = "matches", p.cutoff = pval_thresh)
            ))
            motif_ids <- all_motif_ids
            all_mat <- rbind(
                fg_mat[, motif_ids, drop = FALSE],
                bg_mat[, motif_ids, drop = FALSE]
            )
            clust_list <- split(seq_along(all_clusters), all_clusters)

            # Detect pure-label components: no FG/BG mixing -> permutation non-informative
            n_mixed <- sum(vapply(clust_list, function(cl) {
                length(unique(all_labels[cl])) > 1L
            }, integer(1)))
            if (n_mixed == 0L) {
                message("    Motif permutation skipped: no component contains both FG and BG anchors; within-component label permutation is non-informative. Retaining Fisher exact test.")
            } else {
                permutation_ran <- TRUE
                null_log_or <- matrix(NA_real_,
                    nrow = n_perm, ncol = length(motif_ids),
                    dimnames = list(NULL, motif_ids)
                )
                for (p in seq_len(n_perm)) {
                    perm_labels <- all_labels
                    for (cl in clust_list) {
                        perm_labels[cl] <- sample(perm_labels[cl])
                    }
                    p_fg <- which(perm_labels)
                    p_bg <- which(!perm_labels)
                    p_fg_hits <- colSums(all_mat[p_fg, , drop = FALSE])
                    p_bg_hits <- colSums(all_mat[p_bg, , drop = FALSE])
                    for (j in seq_along(motif_ids)) {
                        a <- p_fg_hits[j]
                        b <- p_bg_hits[j]
                        null_log_or[p, j] <- .motif_log_or(a, b, n_fg, n_bg)
                    }
                }
                # Empirical P-value: fraction of null log-OR >= observed log-OR
                res_df$Empirical_Pvalue <- res_df$Fisher_Pvalue
                for (j in seq_along(motif_ids)) {
                    obs_log_or <- res_df$LogOddsRatio[res_df$MotifID == motif_ids[j]]
                    if (length(obs_log_or) == 0L || !is.finite(obs_log_or)) next
                    n_ge <- sum(null_log_or[, j] >= obs_log_or, na.rm = TRUE)
                    res_df$Empirical_Pvalue[res_df$MotifID == motif_ids[j]] <-
                        (n_ge + 1) / (n_perm + 1)
                }
                res_df$Pvalue <- res_df$Empirical_Pvalue
                message(sprintf(
                    "    Motif component-block permutation: %d iterations, %d components (%d mixed), %d motifs",
                    n_perm, length(clust_list), n_mixed, length(motif_ids)
                ))
            }
        } else {
            message("    Motif permutation skipped: no cluster_id in anchor metadata.")
        }
    }

    if (!is.null(res_df) && nrow(res_df) > 0) {
        if (permutation_ran) {
            res_df$Pvalue_Used <- res_df$Empirical_Pvalue
            res_df$Pvalue <- res_df$Pvalue_Used
            res_df$Pvalue_Method <- "within_component_permutation"
            res_df$Permutation_N <- n_perm
            res_df$N_Components <- length(clust_list)
            res_df$N_Mixed_Components <- n_mixed
        }
        res_df$FDR <- p.adjust(res_df$Pvalue, method = "BH")
        res_df <- res_df[order(res_df$FDR), ]
    }
    return(res_df)
}

#' @title Plot and Save Motif Results (Barplot)
#' @return A \code{ggplot} barplot, or \code{NULL} if input is empty.
#' @keywords internal
#' @noRd
.plot_save_motif <- function(res_df, prefix) {
    if (is.null(res_df) || nrow(res_df) == 0) {
        return(NULL)
    }
    top_df <- head(res_df, 15)
    top_df$MotifLabel <- factor(paste0(top_df$MotifName, " (", top_df$MotifID, ")"), levels = rev(paste0(top_df$MotifName, " (", top_df$MotifID, ")")))
    return(ggplot2::ggplot(top_df, ggplot2::aes(x = -log10(.data$FDR), y = MotifLabel)) +
        ggplot2::geom_col(fill = "#E7298A", width = 0.7) +
        ggplot2::labs(title = paste0("Motif Enrichment: ", basename(prefix)), x = "-log10(FDR)", y = NULL) +
        ggplot2::theme_classic())
}

#' @title Plot Top Motif Sequence Logos
#' @return A list of sequence logo plots, or \code{NULL} if input is empty.
#' @keywords internal
#' @noRd
.plot_top_motif_logos <- function(
  res_df, top_n,
  jaspar_db = NULL
) {
    if (is.null(res_df) || nrow(res_df) == 0) {
        return(NULL)
    }
    if (!requireNamespace("TFBSTools", quietly = TRUE)) {
        message("Motif family annotation skipped: install 'TFBSTools' package.")
        return(NULL)
    }
    if (is.null(jaspar_db)) {
        if (!requireNamespace("JASPAR2020", quietly = TRUE)) {
            message("Motif family annotation skipped: install 'JASPAR2020' package.")
            return(NULL)
        }
        jaspar_db <- JASPAR2020::JASPAR2020
    }
    top_df <- head(res_df[order(res_df$FDR), ], top_n)
    pfm_list <- TFBSTools::getMatrixSet(jaspar_db, opts = list(ID = top_df$MotifID))
    if (length(pfm_list) == 0) {
        message("Motif logo plot skipped: no PFM matrices found for top motifs.")
        return(NULL)
    }

    plot_list <- list()
    for (i in seq_along(top_df$MotifID)) if (top_df$MotifID[i] %in% names(pfm_list)) plot_list[[paste0(top_df$MotifName[i], " (", top_df$MotifID[i], ")")]] <- TFBSTools::Matrix(pfm_list[[top_df$MotifID[i]]])
    if (length(plot_list) == 0) {
        message("Motif logo plot skipped: no matching PFM matrices for top motifs.")
        return(NULL)
    }

    if (!requireNamespace("ggseqlogo", quietly = TRUE)) {
        message("Motif logo plot skipped: install 'ggseqlogo' package.")
        return(NULL)
    }
    return(ggseqlogo::ggseqlogo(plot_list, ncol = 1) + ggplot2::theme_classic() + ggplot2::theme(axis.text.x = ggplot2::element_blank(), strip.text = ggplot2::element_text(size = 10, face = "bold", hjust = 0), strip.background = ggplot2::element_rect(fill = "grey95", color = NA)) + ggplot2::labs(y = "Bits", title = paste0("Top ", top_n, " Enriched Motifs (SeqLogo)")))
}

#' @title Annotate Motif Families
#' @return A data frame with added \code{Family} column, sorted by P-value.
#' @keywords internal
#' @noRd
.annotate_motif_families <- function(
  res_df,
  jaspar_db = NULL, jaspar_collection = "CORE"
) {
    if (is.null(res_df) || nrow(res_df) == 0) {
        return(res_df)
    }
    if (!requireNamespace("TFBSTools", quietly = TRUE)) {
        message("Motif family annotation skipped: install 'TFBSTools' package.")
        return(res_df)
    }
    if (is.null(jaspar_db)) {
        if (!requireNamespace("JASPAR2020", quietly = TRUE)) {
            message("Motif family annotation skipped: install 'JASPAR2020' package.")
            return(res_df)
        }
        jaspar_db <- JASPAR2020::JASPAR2020
    }
    meta_df <- do.call(rbind, lapply(TFBSTools::getMatrixSet(jaspar_db, list(collection = jaspar_collection)), function(x) data.frame(MotifID = TFBSTools::ID(x), Family = paste(if (is.null(TFBSTools::tags(x)$family)) "Unknown" else TFBSTools::tags(x)$family, collapse = "; "), stringsAsFactors = FALSE)))
    res_df <- merge(res_df[, !colnames(res_df) %in% c("Family", "Class")], meta_df, by = "MotifID", all.x = TRUE)
    return(res_df[order(res_df$FDR), c(c("MotifID", "MotifName", "Family", "Pvalue", "FDR", "OddsRatio", "FG_Hits", "FG_Total", "BG_Hits", "BG_Total"), setdiff(colnames(res_df), c("MotifID", "MotifName", "Family", "Pvalue", "FDR", "OddsRatio", "FG_Hits", "FG_Total", "BG_Hits", "BG_Total")))])
}


#' Render a Publication-Ready 3D Annotation Report
#'
#' One-click parameterised R Markdown report that executes the full looplook
#' pipeline (annotation -> refinement -> profiling) and renders an
#' interpretation-ready HTML document suitable for sharing with collaborators.
#'
#' @details
#' The profiling stage uses R's random number generator for GSEA gene-set
#' down-sampling (\code{gsea_nSample}) and motif background anchor sampling. Call
#' \code{set.seed()} before \code{looplook_report()} for fully reproducible
#' results.
#'
#' @param bedpe_file Path to a BEDPE file of chromatin loops. Default \code{NULL}
#'   (requires \code{precomputed_res} instead).
#' @param target_bed Optional path to a BED file of genomic features.
#' @param expr_matrix_file Optional path to a normalised expression matrix.
#' @param sample_columns Sample columns in the expression matrix to average.
#'   Default \code{NULL} (all columns).
#' @param species Genome assembly (\code{"hg38"}, \code{"hg19"}, \code{"mm10"}, \code{"mm9"}).
#' @param project_name Character. Project prefix for the report title.
#' @param out_dir Output directory. Created if missing.
#' @param threshold Numeric. Expression threshold for active gene classification.
#'   Default \code{1.0}.
#' @param reclassify_by_expression Logical. Reclassify silent promoters as eP/eG.
#'   Default \code{TRUE}.
#' @param run_go Logical. Run GO enrichment (requires clusterProfiler).
#'   Default \code{FALSE}.
#' @param universe_genes Named numeric vector or \code{NULL}. GO background
#'   universe. Passed to \code{\link{profile_target_genes}}. Default \code{NULL}.
#' @param diff_file Optional differential expression result file.
#' @param lfc_col Column name for log2 fold change in \code{diff_file}.
#' @param metadata_file Optional sample metadata file.
#' @param precomputed_res Optional. Either a \code{.RData} file path or an
#'   in-memory list object returned by \code{annotate_peaks_and_loops}.
#'   When provided, annotation is skipped and refinement starts from this object.
#' @param chromatin_beds Named list of BED file paths for orthogonal chromatin
#'   mark validation (passed to \code{\link{refine_loop_anchors_by_expression}}).
#'   When non-empty, a \emph{Chromatin Validation} section appears in the report
#'   with confidence-level distribution for eP/eG anchors. Default: \code{list()} (skip).
#' @param chromatin_bw Named list of bigWig file paths, or \code{NULL}. Passed
#'   to \code{\link{refine_loop_anchors_by_chromatin}}. Requires \code{"H3K4me1"}
#'   and \code{"H3K4me3"}. \strong{Strongly recommended} for resolving true
#'   dual-signature elements. Default: \code{NULL}.
#' @param bw_ratio_threshold Numeric. Minimum H3K4me1/H3K4me3 ratio for dual
#'   classification. Default: \code{3}.
#' @param enhancer_bed Character or \code{NULL}. Path to a curated enhancer
#'   BED file (e.g. FANTOM5, ENCODE cCREs). Default: \code{NULL}.
#' @param unit_type Character. Expression unit label for plot annotations. Default \code{"TPM"}.
#' @param tss_region Numeric vector of length 2. TSS flanking region in bp. Default \code{c(-2000, 2000)}.
#' @param neighbor_hop Integer. k-hop ego-network expansion order for loop connectivity analysis. Default \code{0}.
#' @param hub_percentile Numeric. Quantile threshold for hub classification. Default \code{0.95}.
#' @param color_palette Character. RColorBrewer qualitative palette name. Default \code{"Set2"}.
#' @param target_mapping_mode Character. Mapping strategy for target genes. Default \code{"all"}.
#' @param include_Filled Logical. Include comprehensively merged gene assignments. Default \code{TRUE}.
#' @param use_nearest_gene Logical. Bypass 3D loop-based assignment and use linear proximity. Default \code{FALSE}.
#' @param target_source Character vector. Source of target genes to profile. Default \code{"targets"}.
#' @param loop_types Character vector. Loop types to include in profiling. Default \code{c("E-P", "P-P")}.
#' @param stat_test Character. Statistical test for violin comparisons. Default \code{"wilcox.test"}.
#' @param run_ppi Logical. Run protein-protein interaction network analysis. Default \code{FALSE}.
#' @param run_motif Logical. Run transcription factor motif analysis. Default \code{FALSE}.
#' @param genome_id Character. Reference genome for motif scanning. Defaults to \code{species}.
#' @param motif_p_thresh Numeric. P-value threshold for motif enrichment. Default \code{1e-4}.
#' @param motif_ntop Numeric. Number of top motifs to display. Default \code{5}.
#' @param motif_n_perm Integer. Number of within-component label permutations
#'   for empirical motif P-value estimation. Passed to
#'   \code{\link{profile_target_genes}}. Default \code{0} (Fisher exact test).
#' @param ppi_score Numeric. Minimum STRING combined score. Default \code{400}.
#' @param ppi_nSample Numeric. Maximum genes to include in PPI. Default \code{400}.
#' @param heatmap_nSample Numeric. Maximum genes in expression heatmap. Default \code{99999}.
#' @param gsea_nSample Numeric. Maximum genes sampled for GSEA. Default \code{99999}.
#' @param cnet_nSample Numeric. Number of GO terms in cnetplot. Default \code{50}.
#' @param karyo_bin_size Numeric. Bin size for karyotype heatmaps. Default \code{100000}.
#' @param output_file Character. Output HTML file name. \code{NULL} derives from \code{project_name}.
#' @param quiet Logical. Suppress rendering output. Default \code{FALSE}.
#' @param seed Integer or NULL. Passed to \code{\link{profile_target_genes}} for
#'   reproducible GSEA and motif sampling. Default \code{NULL}.
#' @param ... Additional arguments passed to \code{rmarkdown::render}.
#'
#' @return The path to the generated HTML report (invisibly).
#'
#' @export
#'
#' @examples
#' if (requireNamespace("rmarkdown", quietly = TRUE) &&
#'     requireNamespace("knitr", quietly = TRUE)) {
#'     temp_env <- new.env()
#'     load(system.file("extdata", "analysis_results.RData", package = "looplook"), envir = temp_env)
#'     precomputed_res <- temp_env[[ls(temp_env)[1]]]
#'     precomputed_res$loop_annotation <- head(precomputed_res$loop_annotation, 6)
#'     precomputed_res$target_annotation <- head(precomputed_res$target_annotation, 3)
#'     precomputed_res$promoter_centric_stats <- head(precomputed_res$promoter_centric_stats, 6)
#'     precomputed_res$distal_element_stats <- head(precomputed_res$distal_element_stats, 6)
#'
#'     report_path <- looplook_report(
#'         precomputed_res = precomputed_res,
#'         project_name = "Example",
#'         out_dir = tempdir(),
#'         output_file = "looplook-example-report.html",
#'         quiet = TRUE,
#'         run_go = FALSE,
#'         run_ppi = FALSE,
#'         run_motif = FALSE
#'     )
#'     file.exists(report_path)
#' }
looplook_report <- function(
  bedpe_file = NULL,
  target_bed = NULL,
  expr_matrix_file = NULL,
  sample_columns = NULL,
  species = "hg38",
  project_name = "looplook Analysis",
  out_dir = "looplook_results",
  threshold = 1.0,
  unit_type = "TPM",
  reclassify_by_expression = TRUE,
  tss_region = c(-2000, 2000),
  neighbor_hop = 0,
  hub_percentile = 0.95,
  color_palette = "Set2",
  target_mapping_mode = "all",
  include_Filled = TRUE,
  use_nearest_gene = FALSE,
  target_source = "targets",
  loop_types = c("E-P", "P-P"),
  stat_test = "wilcox.test",
  run_go = FALSE,
  run_ppi = FALSE,
  run_motif = FALSE,
  genome_id = species,
  motif_p_thresh = 1e-4,
  motif_ntop = 5,
  motif_n_perm = 0L,
  ppi_score = 400,
  ppi_nSample = 400,
  heatmap_nSample = 99999,
  gsea_nSample = 99999,
  cnet_nSample = 50,
  karyo_bin_size = 1e5,
  diff_file = NULL,
  lfc_col = "log2FoldChange",
  metadata_file = NULL,
  precomputed_res = NULL,
  chromatin_beds = list(),
  chromatin_bw = NULL,
  bw_ratio_threshold = 3,
  enhancer_bed = NULL,
  output_file = NULL,
  quiet = FALSE,
  seed = NULL,
  universe_genes = NULL,
  ...
) {
    if (!is.character(species) || length(species) != 1L || is.na(species) || !nzchar(species)) {
        stop("`species` must be a single non-empty string", call. = FALSE)
    }
    normalize_report_path <- function(path) {
        if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
            return(path)
        }
        if (!file.exists(path)) {
            return(path)
        }
        normalizePath(path, mustWork = TRUE)
    }

    # Locate template
    template <- system.file("rmarkdown", "templates", "looplook-report",
        "skeleton", "skeleton.Rmd",
        package = "looplook"
    )
    if (!nzchar(template)) stop("Report template not found. Reinstall looplook.")
    if (!requireNamespace("rmarkdown", quietly = TRUE)) {
        stop("Package 'rmarkdown' is required for looplook_report(). Install it with: install.packages('rmarkdown')", call. = FALSE)
    }

    # Create output directory
    if (!dir.exists(out_dir)) {
        dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    }
    out_dir <- normalizePath(out_dir, mustWork = TRUE)

    bedpe_file <- normalize_report_path(bedpe_file)
    target_bed <- normalize_report_path(target_bed)
    expr_matrix_file <- normalize_report_path(expr_matrix_file)
    diff_file <- normalize_report_path(diff_file)
    metadata_file <- normalize_report_path(metadata_file)
    if (is.character(precomputed_res)) {
        precomputed_res <- normalize_report_path(precomputed_res)
    }

    # Prepare output filename
    if (is.null(output_file)) {
        output_file <- paste0(gsub("[[:space:]]+", "_", project_name), "_Report.html")
    }

    if (!quiet) message(">> Generating looplook report: ", output_file)

    # Render report
    out_path <- rmarkdown::render(
        input = template,
        params = list(
            bedpe_file = bedpe_file,
            target_bed = target_bed,
            expr_matrix_file = expr_matrix_file,
            sample_columns = sample_columns,
            tss_region = tss_region,
            species = species,
            project_name = project_name,
            out_dir = out_dir,
            threshold = threshold,
            unit_type = unit_type,
            reclassify_by_expression = reclassify_by_expression,
            neighbor_hop = neighbor_hop,
            hub_percentile = hub_percentile,
            color_palette = color_palette,
            target_mapping_mode = target_mapping_mode,
            include_Filled = include_Filled,
            use_nearest_gene = use_nearest_gene,
            target_source = target_source,
            loop_types = loop_types,
            stat_test = stat_test,
            run_go = run_go,
            universe_genes = universe_genes,
            run_ppi = run_ppi,
            run_motif = run_motif,
            genome_id = genome_id,
            motif_p_thresh = motif_p_thresh,
            motif_ntop = motif_ntop,
            motif_n_perm = motif_n_perm,
            ppi_score = ppi_score,
            ppi_nSample = ppi_nSample,
            heatmap_nSample = heatmap_nSample,
            gsea_nSample = gsea_nSample,
            cnet_nSample = cnet_nSample,
            karyo_bin_size = karyo_bin_size,
            diff_file = diff_file,
            lfc_col = lfc_col,
            metadata_file = metadata_file,
            precomputed_res = precomputed_res,
            chromatin_beds = chromatin_beds,
            chromatin_bw = chromatin_bw,
            bw_ratio_threshold = bw_ratio_threshold,
            enhancer_bed = enhancer_bed, seed = seed
        ),
        output_dir = out_dir,
        output_file = output_file,
        quiet = quiet,
        ...
    )

    if (!quiet) message(">> Report saved to: ", out_path)
    return(invisible(out_path))
}
