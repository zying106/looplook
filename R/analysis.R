#' @title End-to-end functional annotation pipeline integrating JASPAR transcription factor motif analysis, gene ontology enrichment, and protein-protein interaction network analysis
#'
#' @description
#' This function serves as a fully automated, end-to-end multi-omics analysis pipeline.
#' It elegantly bridges 3D genomic interactions (e.g., Hi-C, HiChIP) with transcriptomic data (RNA-seq)
#' to systematically characterize the functional landscape and regulatory mechanisms underlying the identified target genes.
#'
#' The suite integrates multiple analytical modules:
#' \itemize{
#'   \item \strong{Differential Expression (LFC) Profiling:} Compares the Log2 Fold Change of target genes against the genomic background using tailored violin/boxplots.
#'   \item \strong{GSEA Analysis:} Evaluates the rank-distribution of custom target gene sets within the global differential expression landscape.
#'   \item \strong{Expression & 3D Connectivity Dynamics:} Generates Z-score expression heatmaps and correlates 3D structural connectivity (e.g., Hub degree) with gene expression/LFC using scatter and sophisticated Raincloud plots.
#'   \item \strong{TF Motif Analysis:} Scans JASPAR core motifs across proximal and distal loop anchors, generating motif enrichment rank scatter plots and sequence logos (SeqLogos).
#'   \item \strong{Functional Enrichment (GO):} Performs Gene Ontology (BP) enrichment, visualizing results through divergent concept networks (cnetplot) and aggregated faceted lollipop charts.
#'   \item \strong{Interaction Networks (PPI):} Constructs and visualizes high-confidence Protein-Protein Interaction networks via the STRING database.
#' }
#'
#'
#' @param annotation_res List. The result object returned by \code{\link{annotate_peaks_and_loops}} or \code{\link{refine_loop_anchors_by_expression}}.
#' @param diff_file Character. Path to the differential expression file (CSV/TSV). Must contain gene IDs (rownames) and a Log Fold Change column.
#' @param lfc_col Character. The column name in \code{diff_file} representing Log2 Fold Change (e.g., "log2FoldChange").
#' @param expr_matrix_file Character. Path to the normalized expression matrix (Genes x Samples, e.g., TPM or FPKM).
#' @param metadata_file Character. Path to the sample metadata file. Must contain "SampleID" and "Group" columns.
#' @param target_source Character vector. Source of target genes to analyze. Options include \code{"loops"} (genes connected by specific loops) and \code{"targets"} (genes overlapping with external target BED).
#' @param target_mapping_mode Character. Specifies the mapping strategy for target genes. Must be one of:
#'   \itemize{
#'     \item \code{"all"} (Default): Accepts all physically connected 3D targets (including Gene Bodies).
#'     \item \code{"promoter"}: Strictly enforces that the 3D loop must explicitly anchor at a Promoter region.
#'   }
#' @param loop_types Character vector. The specific loop types to analyze (e.g., \code{c("E-P", "P-P")}). Active only when \code{"loops"} is included in \code{target_source}.
#' @param include_Filled Logical. If \code{TRUE}, utilizes the comprehensively merged gene assignment (retaining both loop-assigned targets and local host genes).
#' @param use_nearest_gene Logical. If \code{TRUE}, bypasses 3D loop-based gene assignment and strictly uses the nearest 1D gene (\code{SYMBOL}), serving as a spatial baseline reference.
#' @param group_order Character vector. Optional factor levels to sort sample groups in visualizations.
#' @param project_name Character. Prefix for all output files and plot titles.
#' @param out_dir Character. Directory path for saving output PDF plots and CSV tables.
#' @param org_db Character. Organism annotation database (e.g., "org.Hs.eg.db" for human, "org.Mm.eg.db" for mouse).
#' @param run_motif Logical. Whether to perform Transcription Factor Binding Site (TFBS) motif analysis on proximal and distal anchors.
#' @param genome_id Character. Reference genome assembly for motif sequence extraction (e.g., "hg38", "hg19", "mm10").
#' @param motif_p_thresh Numeric. P-value threshold for \code{motifmatchr} scanning.
#' @param motif_ntop Numeric. Number of top enriched motifs to output as SeqLogos.
#' @param run_go Logical. Whether to perform Gene Ontology (GO) ORA enrichment.
#' @param run_ppi Logical. Whether to construct Protein-Protein Interaction networks via the STRING database (requires internet access).
#' @param ppi_score Numeric. Minimum combined confidence score for STRING interaction edges (0-1000).
#' @param ppi_nSample Numeric. Maximum number of top variable/LFC genes to include in the PPI network to prevent visual overcrowding.
#' @param heatmap_nSample Numeric. Maximum number of highly variable genes to plot in the expression heatmap.
#' @param gsea_nSample Numeric. Maximum number of target genes to sample for GSEA to maintain computational efficiency.
#' @param cnet_nSample Numeric. Number of top GO terms to display in the divergent concept network plot.
#' @param stat_test Character. Statistical test for LFC comparisons (\code{"wilcox.test"} or \code{"t.test"}).
#' @param cor_method Character. Method for sample correlation matrices (\code{"pearson"} or \code{"spearman"}).
#'
#' @return An invisible nested list containing:
#' \itemize{
#'   \item \code{go_results}: A list of top GO enrichment tables for each analyzed data source.
#'   \item \code{target_gene_sets}: The specific lists of gene IDs utilized for each analysis task.
#' }
#' @export
#'
#' @examples
#' # =========================================================================
#' # 1. Get paths to real example files included in the package
#' # =========================================================================
#' rdata_path <- system.file("extdata", "analysis_results.RData", package = "looplook")
#' expr_path <- system.file("extdata", "example_tpm.txt", package = "looplook")
#' diff_path <- system.file("extdata", "example_deg.txt", package = "looplook")
#' meta_path <- system.file("extdata", "example_coldata.txt", package = "looplook")
#'
#' # Check if all necessary files exist before running
#' if (rdata_path != "" && expr_path != "" && diff_path != "" && meta_path != "") {
#'   # Safely load the pre-computed annotation result from RData
#'   temp_env <- new.env()
#'   load(rdata_path, envir = temp_env)
#'   raw_annotation <- temp_env[[ls(temp_env)[1]]]
#'
#'   out_base <- tempdir()
#'
#'   # =========================================================================
#'   # Scenario A: Comprehensive Integrative Profiling (Loops + Targets)
#'   #
#'   # Application:
#'   # Best for complete HiChIP + ATAC-seq (or ChIP-seq) joint analysis.
#'   # It maps genes regulated by BOTH 3D distal loops and 1D direct target overlaps.
#'   # =========================================================================
#'   res_A <- profile_target_genes(
#'     annotation_res = raw_annotation,
#'     diff_file = diff_path,
#'     lfc_col = "log2FoldChange",
#'     expr_matrix_file = expr_path,
#'     metadata_file = meta_path,
#'     target_source = c("loops", "targets"),
#'     target_mapping_mode = "all",
#'     include_Filled = TRUE,
#'     use_nearest_gene = FALSE,
#'     project_name = "Scenario_A_Integrative",
#'     out_dir = out_base,
#'     run_go = FALSE, # Set to TRUE in real analysis
#'     run_ppi = FALSE
#'   )
#'
#'   # =========================================================================
#'   # Scenario B: Peak-Centric & Promoter-Strict Profiling
#'   #
#'   # Application:
#'   # Focuses ONLY on 1D peaks ("targets") and STRICTLY restricts the mapping
#'   # to genes where peaks fall directly on their "promoter".
#'   # Ideal for baseline epigenetic profiling, ignoring distal 3D regulation.
#'   # =========================================================================
#'   res_B <- profile_target_genes(
#'     annotation_res = raw_annotation,
#'     diff_file = diff_path,
#'     lfc_col = "log2FoldChange",
#'     expr_matrix_file = expr_path,
#'     metadata_file = meta_path,
#'     target_source = "targets", # Only 1D targets
#'     target_mapping_mode = "promoter", # Strict promoter mapping
#'     include_Filled = TRUE,
#'     use_nearest_gene = FALSE,
#'     project_name = "Scenario_B_PromoterOnly",
#'     out_dir = out_base,
#'     run_go = FALSE,
#'     run_ppi = FALSE
#'   )
#'
#'   # =========================================================================
#'   # Scenario C: Highly Conservative Distal Enhancer-Targeting
#'   #
#'   # Application:
#'   # Focuses on 3D loops. By setting `use_nearest_gene = TRUE` and
#'   # `include_Filled = FALSE`, it forces the selection of strictly the nearest
#'   # genes to the anchors, preventing the expansion of gene lists from broader
#'   # annotated regions. Ensures a highly conservative gene target list.
#'   # =========================================================================
#'   res_C <- profile_target_genes(
#'     annotation_res = raw_annotation,
#'     diff_file = diff_path,
#'     lfc_col = "log2FoldChange",
#'     expr_matrix_file = expr_path,
#'     metadata_file = meta_path,
#'     target_source = "loops", # Only 3D loops
#'     target_mapping_mode = "all",
#'     include_Filled = FALSE,
#'     use_nearest_gene = TRUE, # Strictly nearest gene
#'     project_name = "Scenario_C_StrictLoops",
#'     out_dir = out_base,
#'     run_go = FALSE,
#'     run_ppi = FALSE
#'   )
#' }
profile_target_genes <- function(
  annotation_res,
  diff_file,
  lfc_col,
  expr_matrix_file,
  metadata_file,
  target_source = c("loops", "targets"),
  target_mapping_mode = c("all", "promoter"),
  loop_types = c("E-P", "P-P"),
  include_Filled = TRUE,
  use_nearest_gene = FALSE,
  group_order = NULL,
  project_name = "Analysis",
  out_dir = "./",
  org_db = "org.Hs.eg.db",
  run_motif = FALSE,
  genome_id = "hg19",
  motif_p_thresh = 1e-4,
  motif_ntop = 5,
  run_go = TRUE,
  run_ppi = FALSE,
  ppi_score = 400,
  ppi_nSample = 400,
  heatmap_nSample = 1000,
  gsea_nSample = 99999,
  cnet_nSample = 50,
  stat_test = "wilcox.test",
  cor_method = "pearson"
) {
  # --- Environment Setup ---
  target_source <- match.arg(target_source, choices = c("loops", "targets"), several.ok = TRUE)
  target_mapping_mode <- match.arg(target_mapping_mode)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  root_project_name <- project_name


  if (use_nearest_gene && "targets" %in% target_source) {
    message(">>> [Mode] Reference Mode Active: Using Nearest Genes (SYMBOL).")
    root_project_name <- paste0(root_project_name, "_RefNearest")
  } else if ("targets" %in% target_source) {
    if (target_mapping_mode == "promoter") {
      root_project_name <- paste0(root_project_name, "_Promoter")
    }
    if (!include_Filled) {
      root_project_name <- paste0(root_project_name, "_LoopOnly")
    }
  }

  message(">>> Analysis Init | Root Project: ", root_project_name)

  # Check Packages
  req_pkgs <- c("RColorBrewer", "ggplot2", "clusterProfiler", "AnnotationDbi", "dplyr", "ComplexHeatmap", "circlize")
  if (!requireNamespace(org_db, quietly = TRUE)) stop("Package ", org_db, " missing. Please install it.")
  for (pkg in req_pkgs) if (!requireNamespace(pkg, quietly = TRUE)) stop("Package ", pkg, " missing")

  # Check Motif Packages
  if (run_motif) {
    motif_pkgs <- c("motifmatchr", "JASPAR2020", "TFBSTools", "BSgenome", "GenomicRanges", "ggseqlogo")
    bs_pkg <- switch(genome_id,
      "hg19" = "BSgenome.Hsapiens.UCSC.hg19",
      "hg38" = "BSgenome.Hsapiens.UCSC.hg38",
      "mm9"  = "BSgenome.Mmusculus.UCSC.mm9",
      "mm10" = "BSgenome.Mmusculus.UCSC.mm10",
      NULL
    )

    if (is.null(bs_pkg)) {
      warning(" Unsupported genome_id: '", genome_id, "'. Supported: hg19, hg38, mm9, mm10. Disabling Motif Analysis.")
      run_motif <- FALSE
    } else {
      motif_pkgs <- c(motif_pkgs, bs_pkg)
      missing_motif <- motif_pkgs[!vapply(motif_pkgs, requireNamespace, logical(1), quietly = TRUE)]
      if (length(missing_motif) > 0) {
        warning(" Missing packages for Motif analysis: ", toString(missing_motif), ". Skipping Motif step.")
        run_motif <- FALSE
      } else {
        message("Motif Analysis Enabled using ", bs_pkg)
      }
    }
  }

  if ("loops" %in% target_source && is.null(annotation_res$loop_annotation)) stop("Missing loop_annotation in input object")

  # --- Data Loading ---
  message("--- Reading files...")
  diff_df_raw <- read_robust_general(diff_file, header = TRUE, row_name = 1, desc = "Diff", min_cols = 1)
  tpm_mat_raw <- read_robust_general(expr_matrix_file, header = TRUE, row_name = 1, desc = "Expr", min_cols = 1)
  meta_raw <- read_robust_general(metadata_file, header = TRUE, row_name = NULL, desc = "Meta", min_cols = 1)

  if (ncol(meta_raw) >= 2) colnames(meta_raw)[c(1, 2)] <- c("SampleID", "Group")
  meta_raw$SampleID <- trimws(as.character(meta_raw$SampleID))
  if (!is.null(group_order)) meta_raw$Group <- factor(meta_raw$Group, levels = group_order)
  meta_raw <- meta_raw %>% dplyr::arrange(Group)

  if (!lfc_col %in% colnames(diff_df_raw)) stop("LFC column", lfc_col, "not found")
  clean_diff <- diff_df_raw[!is.na(diff_df_raw[[lfc_col]]) & is.finite(diff_df_raw[[lfc_col]]), ]
  global_glist <- sort(setNames(clean_diff[[lfc_col]], rownames(clean_diff)), decreasing = TRUE)

  loop_stats_df <- annotation_res$promoter_centric_stats
  if (!is.null(loop_stats_df)) message("Loaded Promoter Centric Stats.")

  final_master_list <- list()

  # --- Main Loop ---
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
    source_go_results <- list()

    for (task_name in names(analysis_queue)) {
      target_genes <- analysis_queue[[task_name]]
      current_proj_name <- paste0(current_source_proj_name, "_", task_name)
      target_genes <- intersect(target_genes, names(global_glist))

      message("\n--- Task: ", task_name, " (Valid Genes: ", length(target_genes), ") ---")
      if (length(target_genes) < 3) {
        message("  Too few genes (<3), skipping.")
        next
      }

      # 3.1 Violin
      tryCatch(
        {
          run_lfc_violin(target_genes, global_glist, stat_test, current_proj_name, out_dir)
        },
        error = function(e) warning("Violin failed: ", e$message)
      )

      # 3.2 GSEA
      tryCatch(
        {
          gsea_res <- run_gsea_analysis(target_genes, global_glist, gsea_nSample, current_proj_name, out_dir)
          if (!is.null(gsea_res)) write.csv(as.data.frame(gsea_res), file.path(out_dir, paste0(current_proj_name, "_GSEA_Results.csv")), row.names = FALSE)
        },
        error = function(e) warning("GSEA failed: ", e$message)
      )

      # 3.3 Visualizations (Total Loops)
      tryCatch(
        {
          run_heatmap_and_connectivity(
            target_genes, tpm_mat_raw, meta_raw, loop_stats_df, global_glist,
            heatmap_nSample, cor_method,
            current_proj_name, out_dir,
            source_type = src,
            target_col = NULL
          )
        },
        error = function(e) warning("Vis (Total): ", e$message)
      )

      # 3.3b Visualizations (Distal Loops)
      if (!is.null(loop_stats_df) && "n_Linked_Distal" %in% colnames(loop_stats_df)) {
        message("    [Extra Analysis] Found 'n_Linked_Distal'. Running distal connectivity analysis...")
        tryCatch(
          {
            run_heatmap_and_connectivity(
              target_genes, tpm_mat_raw, meta_raw, loop_stats_df, global_glist,
              heatmap_nSample, cor_method,
              current_proj_name, out_dir,
              source_type = src,
              target_col = "n_Linked_Distal",
              skip_heatmap = TRUE
            )
          },
          error = function(e) warning("Vis (Distal): ", e$message)
        )
      }

      # 3.4 Motif Analysis (Proximal & Distal)
      if (run_motif) {
        message("    [Motif] Running Motif Analysis (Proximal vs Distal)...")
        tryCatch(
          {
            run_distal_motif_analysis(
              target_genes = target_genes,
              loop_df = annotation_res$loop_annotation,
              genome_id = genome_id,
              pval_thresh = motif_p_thresh,
              current_proj_name = current_proj_name,
              out_dir = out_dir,
              top_n = motif_ntop
            )
          },
          error = function(e) warning("    [Motif] Failed: ", e$message)
        )
      }

      # 3.5 GO
      tryCatch(
        {
          if (run_go) {
            go_res <- run_go_enrichment(target_genes, org_db, global_glist, cnet_nSample, current_proj_name, out_dir)
            if (!is.null(go_res) && nrow(go_res) > 0) {
              write.csv(go_res, file.path(out_dir, paste0(current_proj_name, "_GO_Enrichment_All.csv")), row.names = FALSE)
              if ("ONTOLOGY" %in% colnames(go_res)) {
                top_go <- go_res %>%
                  dplyr::group_by(ONTOLOGY) %>%
                  dplyr::arrange(pvalue) %>%
                  dplyr::slice_head(n = 5) %>%
                  dplyr::ungroup()
                top_go$CleanLoopType <- task_name
                top_go$LoopType <- paste0(task_name, "\n(", top_go$ONTOLOGY, ")")
              } else {
                top_go <- head(go_res[order(go_res$pvalue), ], 15)
                top_go$CleanLoopType <- task_name
                top_go$LoopType <- task_name
              }
              if (nrow(top_go) > 0) {
                top_go$Source <- src
                source_go_results[[length(source_go_results) + 1]] <- top_go
              }
            }
          }
        },
        error = function(e) warning("GO : ", e$message)
      )

      # 3.6 PPI
      tryCatch(
        {
          if (run_ppi) run_ppi_analysis(target_genes, global_glist, org_db, ppi_score, ppi_nSample, current_proj_name, out_dir)
        },
        error = function(e) warning("PPI failed: ", e$message)
      )
    }

    if (run_go && length(source_go_results) > 0) {
      tryCatch(
        {
          plot_summary_go_lollipop(source_go_results, current_source_proj_name, out_dir)
        },
        error = function(e) warning("    Plot Failed: ", e$message)
      )
    }
    final_master_list[[src]] <- list(go_results = source_go_results, target_gene_sets = analysis_queue)
  }
  message("\n All analysis complete.")
  invisible(final_master_list)
}


#' @title Robust Data Reader
#' @description Safely reads standard genomic formats (BED, BEDPE, CSV, TSV) using `data.table::fread` (if available) with progressive fallbacks to base R. Automatically handles delimiters and validates minimum column counts.
#' @param f Character. Path to the input file.
#' @param header Logical. Whether the file contains a header row.
#' @param row_name Integer or NULL. Column index to be used as row names.
#' @param desc Character. Short description for error logging (e.g., "BEDPE").
#' @param min_cols Integer. Minimum number of columns required to pass validation.
#' @return A data frame.
#' @keywords internal
read_robust_general <- function(f, header = FALSE, row_name = NULL, desc = "file", min_cols = 3) {
  if (is.null(f) || length(f) == 0 || f == "") stop(desc, " path is empty.")
  if (!file.exists(f)) stop(desc, " not found: ", f)

  d <- NULL
  rn <- if (is.null(row_name)) NULL else 1

  if (requireNamespace("data.table", quietly = TRUE)) {
    try(
      {
        d_dt <- data.table::fread(f, header = header, data.table = FALSE, showProgress = FALSE)


        if (!is.null(row_name) && ncol(d_dt) > 1) {
          rownames(d_dt) <- d_dt[, 1]
          d_dt <- d_dt[, -1, drop = FALSE]
        }


        if (ncol(d_dt) >= min_cols) d <- d_dt
      },
      silent = TRUE
    )
  }


  if (is.null(d)) {
    d <- tryCatch(
      utils::read.table(f,
        header = header, sep = "\t", row.names = rn,
        check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = ""
      ),
      error = function(e) NULL
    )
  }


  if (is.null(d) || ncol(d) < min_cols) {
    d <- tryCatch(
      utils::read.table(f,
        header = header, sep = ",", row.names = rn,
        check.names = FALSE, stringsAsFactors = FALSE, quote = "\"", comment.char = ""
      ),
      error = function(e) NULL
    )
  }


  if (is.null(d) || ncol(d) < min_cols) {
    d <- tryCatch(
      utils::read.table(f,
        header = header, sep = "", row.names = rn,
        check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = ""
      ),
      error = function(e) NULL
    )
  }


  if (is.null(d)) {
    stop("failed to read ", desc, ". Format not recognized.")
  }

  if (ncol(d) < min_cols) {
    stop(desc, " has insufficient columns (found ", ncol(d), ", required ", min_cols, ").")
  }

  return(d)
}

#' @title Extract Target Gene Sets from Annotation Results
#' @description Parses loop and target annotations using strict column matching to extract valid gene lists for downstream profiling.
#' @param annotation_res List output from `annotate_peaks_and_loops`.
#' @param src Character. Source to extract from ("loops" or "targets").
#' @param active_loop_types Character vector. Which loop types to extract (if src is "loops").
#' @param include_Filled Logical. Whether to use `_Filled` columns containing merged host/target genes.
#' @param use_nearest_gene Logical. Whether to strictly use the 1D nearest gene (SYMBOL) instead of 3D loops.
#' @param target_mapping_mode Character. Mapping mode for targets ("all" or "promoter").
#' @return Named list of gene character vectors.
#' @keywords internal
extract_target_gene_sets <- function(annotation_res, src, active_loop_types = NULL, include_Filled = TRUE, use_nearest_gene = FALSE, target_mapping_mode = "all") {
  raw_gene_sets <- list()

  .clean_split_genes <- function(gene_vec) {
    if (is.null(gene_vec)) {
      return(character(0))
    }
    genes <- unlist(strsplit(as.character(gene_vec), "[;,]"))
    genes <- unique(trimws(genes))
    genes <- genes[genes != "" & genes != "NA" & !is.na(genes)]
    return(genes)
  }

  # --- Source: Targets ---
  if ("targets" %in% src && !is.null(annotation_res$target_annotation)) {
    bed_info <- annotation_res$target_annotation
    target_col <- NULL

    if (use_nearest_gene) {
      if ("SYMBOL" %in% colnames(bed_info)) {
        target_col <- "SYMBOL"
      } else if ("geneId" %in% colnames(bed_info)) target_col <- "geneId"

      if (!is.null(target_col)) {
        message("      Targets: Using Nearest Gene '", target_col, "'")
      } else {
        stop("       Targets: 'SYMBOL' or 'geneId' column required when use_nearest_gene is TRUE.")
      }
    } else {
      base_col <- if (target_mapping_mode == "promoter") "Regulated_promoter_genes" else "Assigned_Target_Genes"
      desired_col <- if (include_Filled) paste0(base_col, "_Filled") else base_col

      if (desired_col %in% colnames(bed_info)) {
        target_col <- desired_col
        message("      Targets: Successfully located '", target_col, "'")
      } else {
        stop("     Targets: Required column '", desired_col, "' not found in target_annotation. Please ensure the annotation step was run correctly.")
      }
    }

    if (!is.null(target_col)) {
      gs <- .clean_split_genes(bed_info[[target_col]])
      if (length(gs) > 0) raw_gene_sets[["Target_Genes"]] <- gs
    }
  }

  # --- Source: Loops ---
  if ("loops" %in% src && !is.null(annotation_res$loop_annotation)) {
    loop_df <- annotation_res$loop_annotation
    gene_col <- "Putative_Target_Genes"

    if (gene_col %in% colnames(loop_df)) {
      message("     Loops: Successfully located '", gene_col, "'")
    } else {
      stop("     Loops: Required column '", gene_col, "' not found in loop_annotation. Please ensure the loop processing step was run correctly.")
    }

    use_types <- if (is.null(active_loop_types)) unique(loop_df$loop_type) else intersect(active_loop_types, unique(loop_df$loop_type))

    if (length(use_types) > 0) {
      for (lt in use_types) {
        sub_df <- loop_df[loop_df$loop_type == lt, ]
        if (nrow(sub_df) > 0) {
          gs <- .clean_split_genes(sub_df[[gene_col]])
          if (length(gs) > 0) {
            safe_name <- paste0(gsub("-", "", lt), "_Genes")
            raw_gene_sets[[safe_name]] <- gs
          }
        }
      }
    } else {
      warning("   No loop types matched. Requested: ", toString(active_loop_types))
    }
  }

  return(raw_gene_sets)
}

#' @title Generate LFC Violin and Boxplot
#' @description Visualizes the Log2 Fold Change distribution of target genes compared to the genomic background.
#' @param target_genes Character vector of target gene symbols.
#' @param global_glist Named numeric vector of global log fold changes (names = gene symbols).
#' @param stat_test Character. Statistical test ("t.test" or "wilcox.test").
#' @param project_name Character. Project-specific prefix for output.
#' @param out_dir Character. Output directory.
#' @return Invisible `NULL`. Saves a narrow-format PDF plot.
#' @keywords internal
run_lfc_violin <- function(target_genes, global_glist, stat_test = "wilcox.test", project_name, out_dir) {
  valid_targets <- intersect(target_genes, names(global_glist))

  if (length(valid_targets) < 3) {
    warning("[LFC Violin] Too few valid targets (<3). Skipping.")
    return(NULL)
  }

  target_lfc <- global_glist[valid_targets]
  other_genes <- setdiff(names(global_glist), valid_targets)
  other_lfc <- global_glist[other_genes]

  plot_data <- data.frame(
    LFC = c(target_lfc, other_lfc),
    Group = factor(
      c(
        rep("Target", length(target_lfc)),
        rep("Background", length(other_lfc))
      ),
      levels = c("Target", "Background")
    )
  )


  n_target <- length(target_lfc)
  n_back <- length(other_lfc)

  p_val <- tryCatch(
    {
      if (stat_test == "wilcox.test") {
        wilcox.test(target_lfc, other_lfc)$p.value
      } else {
        t.test(target_lfc, other_lfc)$p.value
      }
    },
    error = function(e) NA
  )

  p_label <- if (is.na(p_val)) "NA" else formatC(p_val, format = "e", digits = 2)

  x_labels <- c(
    "Target" = paste0("Target\n(n=", n_target, ")"),
    "Background" = paste0("Background\n(n=", n_back, ")")
  )


  y_min <- quantile(plot_data$LFC, 0.001, na.rm = TRUE)
  y_max <- quantile(plot_data$LFC, 0.999, na.rm = TRUE)
  y_pad <- (y_max - y_min) * 0.03


  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    return(NULL)
  }
  cols <- c("Target" = "#E41A1C", "Background" = "#999999")

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Group, y = LFC, fill = Group)) +
    ggplot2::geom_violin(trim = TRUE, alpha = 0.5, color = NA) +
    ggplot2::geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.9, color = "black") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
    ggplot2::scale_x_discrete(labels = x_labels) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::coord_cartesian(ylim = c(y_min - y_pad, y_max + y_pad)) +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::labs(
      title = project_name,
      subtitle = paste0("Stat: ", stat_test, ", P: ", p_label),
      y = "Log2 Fold Change",
      x = NULL
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 12),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10),
      legend.position = "none"
    )


  out_file <- file.path(out_dir, paste0(project_name, "_LFC_Violin.pdf"))
  tryCatch(
    {
      ggplot2::ggsave(out_file, p, width = 3.0, height = 6.0)
      message("    Saved LFC Violin Plot: ", out_file)
    },
    error = function(e) message("    Save Failed: ", e$message)
  )
}

#' @title Run Custom Gene Set Enrichment Analysis (GSEA)
#' @description Evaluates whether target genes are significantly enriched at the extremes (up/down-regulated) of the global ranked expression list.
#' @param target_genes Character vector. Target gene symbols forming the custom gene set.
#' @param global_glist Named numeric vector. Globally ranked LFC values.
#' @param gsea_ntop Integer. Maximum number of target genes to use (randomly downsampled if exceeded).
#' @param current_proj_name Character. Prefix for plot and table.
#' @param out_dir Character. Output directory.
#' @return GSEA result object. Saves a composite GSEA plot (ES curve, barcode, LFC heatmap) to PDF.
#' @keywords internal
run_gsea_analysis <- function(target_genes, global_glist, gsea_ntop, current_proj_name, out_dir) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    return(NULL)
  }
  if (!requireNamespace("enrichplot", quietly = TRUE)) {
    return(NULL)
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    return(NULL)
  }


  curr_glist <- global_glist

  if (any(duplicated(curr_glist))) {
    curr_glist <- curr_glist + runif(length(curr_glist), 0, 1e-6)
    curr_glist <- sort(curr_glist, decreasing = TRUE)
  }

  names(curr_glist) <- toupper(names(curr_glist))


  curr_targets <- toupper(target_genes)


  if (!is.null(gsea_ntop) && length(curr_targets) > gsea_ntop) {
    message(sprintf("    [GSEA Info] Downsampling targets: %d -> %d", length(curr_targets), gsea_ntop))

    curr_targets <- sample(curr_targets, gsea_ntop)
  }


  overlap_count <- length(intersect(curr_targets, names(curr_glist)))
  if (overlap_count < 2) {
    message("    [GSEA Info] Skipped: Overlap too small (<2).")
    return(NULL)
  }


  term_df <- data.frame(term = current_proj_name, gene = curr_targets, stringsAsFactors = FALSE)

  gsea_res <- tryCatch(
    {
      clusterProfiler::GSEA(
        curr_glist,
        TERM2GENE = term_df,
        pvalueCutoff = 1.1,
        minGSSize = 2,
        maxGSSize = 50000,
        verbose = FALSE,
        seed = 123
      )
    },
    error = function(e) {
      warning("    Calculation failed: ", e$message)
      return(NULL)
    }
  )


  if (!is.null(gsea_res) && nrow(gsea_res@result) > 0) {
    out_file <- file.path(out_dir, paste0(current_proj_name, "_GSEA.pdf"))

    tryCatch(
      {
        p_temp <- enrichplot::gseaplot2(gsea_res, geneSetID = 1, subplots = 1)


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
        } else if (!is.null(p_temp$data)) {
          d <- p_temp$data
        }


        if (is.null(d) || !is.data.frame(d) || !("runningScore" %in% colnames(d))) {
          stop("Data extraction from enrichplot::gseaplot2 failed due to package version mismatch.")
        }

        if (!"geneList" %in% colnames(d)) d$geneList <- curr_glist[d$x]
        if (!"position" %in% colnames(d)) {
          gene_at_rank <- names(curr_glist)[d$x]
          d$position <- as.numeric(gene_at_rank %in% curr_targets)
        }


        max_rank <- max(d$x)
        nes_val <- gsea_res@result$NES[1]
        pval_val <- gsea_res@result$pvalue[1]
        main_col <- if (!is.na(nes_val) && nes_val >= 0) "#E41A1C" else "#377EB8"

        # p1: ES Curve
        p1 <- ggplot2::ggplot(d, ggplot2::aes(x = x, y = runningScore)) +
          ggplot2::geom_line(color = main_col, linewidth = 1) +
          ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
          ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(0, max_rank)) +
          ggplot2::theme_bw() +
          ggplot2::labs(
            x = NULL, y = "ES",
            title = paste0(current_proj_name, "\nNES: ", round(nes_val, 3), "  P: ", formatC(pval_val, format = "e", digits = 2))
          ) +
          ggplot2::theme(
            axis.text.x = ggplot2::element_blank(),
            axis.ticks.x = ggplot2::element_blank(),
            panel.grid = ggplot2::element_blank(),
            plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 10)
          )

        # p2: Barcode
        hit_data <- d[d$position == 1, ]
        p2 <- ggplot2::ggplot(hit_data, ggplot2::aes(x = x, y = 1)) +
          ggplot2::geom_segment(ggplot2::aes(xend = x, yend = 0), color = "black", alpha = 0.6) +
          ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(0, max_rank)) +
          ggplot2::scale_y_continuous(expand = c(0, 0)) +
          ggplot2::theme_void() +
          ggplot2::theme(panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5))

        # p3: LFC Colorbar
        lfc_min_bound <- quantile(d$geneList, 0.005, na.rm = TRUE)
        lfc_max_bound <- quantile(d$geneList, 0.995, na.rm = TRUE)

        p3 <- ggplot2::ggplot(d, ggplot2::aes(x = x, y = geneList)) +
          ggplot2::geom_segment(ggplot2::aes(xend = x, yend = 0, color = geneList)) +
          ggplot2::scale_color_gradient2(
            low = "#1B7837",
            mid = "white",
            high = "#762A83",
            midpoint = 0
          ) +
          ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(0, max_rank)) +
          ggplot2::coord_cartesian(ylim = c(lfc_min_bound, lfc_max_bound)) +
          ggplot2::theme_classic() +
          ggplot2::labs(x = "Rank", y = "LFC") +
          ggplot2::theme(legend.position = "none", axis.text.y = ggplot2::element_text(size = 8))

        if (requireNamespace("aplot", quietly = TRUE)) {
          final_p <- aplot::plot_list(p1, p2, p3, ncol = 1, heights = c(2, 0.5, 1.5))
          ggplot2::ggsave(out_file, final_p, width = 6, height = 6)
          message("    Saved GSEA Plot (Custom): ", out_file)
        } else {
          p_fallback <- enrichplot::gseaplot2(gsea_res, geneSetID = 1, title = gsea_res@result$Description[1])
          ggplot2::ggsave(out_file, p_fallback, width = 7, height = 6)
          message("    Saved GSEA Plot (Fallback): ", out_file)
        }
      },
      error = function(e) {
        warning("    [GSEA Plot ] ", e$message)
      }
    )
  }

  return(gsea_res)
}

#' @title Perform GO Enrichment and Generate Network Plot
#' @description Runs Biological Process GO enrichment via `clusterProfiler` and visualizes the top pathways and their core Hub genes in a highly readable, non-overlapping divergent concept network (Cnetplot).
#' @param genes Character vector of target gene symbols.
#' @param org_db Character. Organism database (e.g., "org.Hs.eg.db").
#' @param universe_genes Named numeric vector (LFC) used as background universe and color scale.
#' @param cnet_nSample Integer. Number of top GO pathways to display in the network.
#' @param project_name Character. Prefix for outputs.
#' @param out_dir Character. Output directory.
#' @return Data frame of raw GO enrichment results. Saves a Cnetplot PDF.
#' @keywords internal
run_go_enrichment <- function(genes, org_db, universe_genes, cnet_nSample = 50, project_name = "Analysis", out_dir = "./") {
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    return(NULL)
  }


  clean_genes <- unique(trimws(as.character(genes)))
  clean_genes <- clean_genes[clean_genes != "" & clean_genes != "NA"]

  gene_entrez <- suppressMessages(tryCatch(
    {
      AnnotationDbi::mapIds(get(org_db), keys = clean_genes, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
    },
    error = function(e) NULL
  ))

  valid_entrez <- na.omit(gene_entrez)

  use_symbol_mode <- FALSE
  final_genes <- valid_entrez
  final_keytype <- "ENTREZID"

  if (length(valid_entrez) < 5 || (length(valid_entrez) / length(clean_genes) < 0.1)) {
    message("    [GO] Low mapping. Switching to SYMBOL mode.")
    use_symbol_mode <- TRUE
    final_genes <- clean_genes
    final_keytype <- "SYMBOL"
  }

  final_universe <- NULL
  if (!is.null(universe_genes)) {
    if (use_symbol_mode) {
      final_universe <- names(universe_genes)
    } else {
      univ_entrez <- suppressMessages(tryCatch(
        {
          AnnotationDbi::mapIds(get(org_db), keys = names(universe_genes), column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
        },
        error = function(e) NULL
      ))
      final_universe <- na.omit(univ_entrez)
    }
  }

  run_enrich <- function(input_genes, univ, key_type, p_val) {
    clusterProfiler::enrichGO(
      gene          = input_genes,
      universe      = univ,
      OrgDb         = get(org_db),
      keyType       = key_type,
      ont           = "ALL",
      pAdjustMethod = "BH",
      pvalueCutoff  = p_val,
      qvalueCutoff  = 1.0,
      minGSSize     = 5,
      maxGSSize     = 800,
      readable      = (key_type == "ENTREZID")
    )
  }

  ego <- tryCatch(
    {
      run_enrich(final_genes, final_universe, final_keytype, 0.05)
    },
    error = function(e) NULL
  )
  if (is.null(ego) || nrow(ego) == 0) {
    ego <- tryCatch(
      {
        run_enrich(final_genes, NULL, final_keytype, 0.2)
      },
      error = function(e) NULL
    )
  }
  if (is.null(ego) || nrow(ego) == 0) {
    ego <- tryCatch(
      {
        run_enrich(final_genes, NULL, final_keytype, 0.8)
      },
      error = function(e) NULL
    )
  }

  if (is.null(ego) || nrow(ego) == 0) {
    return(NULL)
  }


  tryCatch(
    {
      if (requireNamespace("enrichplot", quietly = TRUE) && requireNamespace("ggplot2", quietly = TRUE)) {
        top_n <- if (!is.null(cnet_nSample)) min(cnet_nSample, 5) else 5
        fc_vec <- universe_genes


        genes_to_label <- c()
        if (!is.null(ego) && nrow(ego@result) > 0) {
          top_df <- head(ego@result, top_n)

          gene_to_pathways <- list()
          for (i in seq_len(nrow(top_df))) {
            genes <- unlist(strsplit(top_df$geneID[i], "/"))
            for (g in genes) {
              gene_to_pathways[[g]] <- c(gene_to_pathways[[g]], top_df$ID[i])
            }
          }


          for (i in seq_len(nrow(top_df))) {
            genes <- unlist(strsplit(top_df$geneID[i], "/"))
            valid_g <- intersect(genes, names(fc_vec))
            if (length(valid_g) > 0) {
              g_sorted <- valid_g[order(abs(fc_vec[valid_g]), decreasing = TRUE)]
              genes_to_label <- c(genes_to_label, head(g_sorted, 3))
            }
          }

          hub_genes <- names(gene_to_pathways)[vapply(gene_to_pathways, length, integer(1)) >= 2]
          valid_hub <- intersect(hub_genes, names(fc_vec))
          if (length(valid_hub) > 0) {
            hub_sorted <- valid_hub[order(abs(fc_vec[valid_hub]), decreasing = TRUE)]
            genes_to_label <- c(genes_to_label, head(hub_sorted, 5))
          }


          genes_to_label <- unique(genes_to_label)
          message(sprintf("    [GO Plot] Intelligently selected %d key genes to highlight.", length(genes_to_label)))

          ego@result$Description <- vapply(ego@result$Description, function(x) {
            paste(strwrap(x, width = 35), collapse = "\n")
          }, FUN.VALUE = character(1))
        }


        options(ggrepel.max.overlaps = 100)

        p_cnet <- enrichplot::cnetplot(
          ego,
          foldChange = fc_vec,
          circular = FALSE,
          colorEdge = TRUE,
          showCategory = top_n,
          node_label = "category",
          repel = TRUE,
          cex_label_category = 1.3
        )


        if (length(genes_to_label) > 0 && requireNamespace("ggraph", quietly = TRUE)) {
          p_cnet <- p_cnet +
            ggraph::geom_node_text(
              ggplot2::aes(filter = name %in% genes_to_label, label = name),
              repel = TRUE,
              size = 3.5,
              fontface = "bold.italic",
              bg.color = "white",
              bg.r = 0.15,
              max.overlaps = Inf
            )
        }

        p_cnet <- p_cnet +
          ggplot2::scale_color_gradient2(name = "Log2FC", low = "#57992B", mid = "white", high = "#CD2E91", midpoint = 0) +
          ggplot2::labs(title = paste0("GO Network: ", project_name)) +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))

        ggplot2::ggsave(file.path(out_dir, paste0(project_name, "_GO_Network.pdf")), p_cnet, width = 12, height = 9)
      }
    },
    error = function(e) message("    [GO Plot] Cnet plot skipped/failed: ", e$message)
  )

  return(as.data.frame(ego))
}

#' @title Construct and Visualize STRING PPI Network
#' @description Maps target genes to STRING database, extracts the high-confidence protein-protein interaction subnetwork, removes isolated nodes, and visualizes the network colored by LFC and sized by degree.
#' @param target_genes Character vector of gene symbols.
#' @param global_glist Named numeric vector of LFC values.
#' @param org_db Character. Organism database to automatically infer STRING species ID.
#' @param ppi_score Numeric. Minimum combined score threshold for STRING edges (e.g., 400).
#' @param ppi_ntop Integer. Maximum number of genes to include (prioritized by highest absolute LFC).
#' @param current_proj_name Character. Project prefix.
#' @param out_dir Character. Output directory.
#' @return Invisible `NULL`. Saves a `ggraph` network PDF.
#' @keywords internal
run_ppi_analysis <- function(target_genes, global_glist, org_db, ppi_score, ppi_ntop, current_proj_name, out_dir) {
  message("    >>> [PPI Init] Threshold = ", ppi_score)


  species_id <- if (grepl("Mm", org_db, ignore.case = TRUE)) 10090 else 9606
  string_db_obj <- tryCatch(
    {
      STRINGdb::STRINGdb$new(version = "11.5", species = species_id, score_threshold = ppi_score, input_directory = "")
    },
    error = function(e) {
      message("    [PPI] Init failed: ", e$message)
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
    message(sprintf("    [PPI] Downsampled to top %d genes by LFC.", length(ppi_genes)))
  }

  # 3. Map & Get Network
  targets_mapped <- string_db_obj$map(data.frame(gene = ppi_genes), "gene", removeUnmappedRows = TRUE)
  if (nrow(targets_mapped) == 0) {
    return()
  }

  hits <- targets_mapped$STRING_id
  if (length(hits) <= 1) {
    return()
  }

  g_string <- string_db_obj$get_subnetwork(hits)

  original_node_count <- igraph::vcount(g_string)
  g_string <- igraph::delete_vertices(g_string, igraph::V(g_string)[igraph::degree(g_string) == 0])
  filtered_node_count <- igraph::vcount(g_string)

  if (filtered_node_count == 0) {
    message("    [PPI] No interacting nodes left after removing isolated points.")
    return()
  }

  if (original_node_count > filtered_node_count) {
    message(sprintf(
      "    [PPI] Removed %d isolated nodes. Remaining nodes: %d",
      original_node_count - filtered_node_count, filtered_node_count
    ))
  }

  map_df <- targets_mapped[targets_mapped$STRING_id %in% igraph::V(g_string)$name, ]
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
    ggraph::geom_edge_link(ggplot2::aes(alpha = combined_score), color = "grey60", width = 0.5, show.legend = FALSE) +
    ggraph::geom_node_point(ggplot2::aes(color = lfc, size = deg), stroke = 0.5) +
    ggraph::geom_node_text(ggplot2::aes(label = label_text), repel = TRUE, size = 3.5, max.overlaps = Inf, fontface = "bold", bg.color = "white", bg.r = 0.1) +
    ggplot2::scale_color_gradient2(low = "#57992B", mid = "white", high = "#CD2E91", midpoint = 0, name = "LFC") +
    ggplot2::scale_size_continuous(range = c(2, 8), guide = "none") +
    ggraph::scale_edge_alpha_continuous(range = c(0.4, 0.9)) +
    ggraph::theme_graph(base_family = "sans", background = "white") +
    ggplot2::labs(
      title = paste0("PPI Network: ", current_proj_name),
      subtitle = paste0("Interacting Nodes: ", num_nodes, " | Score: ", min(igraph::E(g_string)$combined_score, na.rm = TRUE), "-", max(igraph::E(g_string)$combined_score, na.rm = TRUE))
    ) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))

  ggplot2::ggsave(file.path(out_dir, paste0(current_proj_name, "_PPI_Network_ggraph.pdf")), p_ppi, width = 9, height = 9)
}


#' @title Generate Summary GO Lollipop Facet Plot
#' @description Combines GO enrichment results from multiple loop types (e.g., E-P, P-P) and plots them side-by-side in a faceted lollipop chart, categorized by GO Ontology (BP, CC, MF).
#' @param all_go_results List of data frames containing GO results.
#' @param base_project_name Character. Prefix for output.
#' @param out_dir Character. Output directory.
#' @return Invisible `NULL`. Saves a combined PDF plot.
#' @keywords internal
plot_summary_go_lollipop <- function(all_go_results, base_project_name, out_dir) {
  message("\n>>> Plotting Summary GO Lollipop...")
  valid_results <- Filter(function(x) !is.null(x) && is.data.frame(x) && nrow(x) > 0, all_go_results)

  if (length(valid_results) == 0) {
    message("    [Info] No valid GO results to plot.")
    return()
  }

  final_go_df <- do.call(rbind, valid_results)
  if (is.null(final_go_df)) {
    return()
  }

  if (!"CleanLoopType" %in% colnames(final_go_df)) {
    final_go_df$CleanLoopType <- final_go_df$LoopType
  }

  use_ggtext <- requireNamespace("ggtext", quietly = TRUE)
  if (!use_ggtext) {
    message("      [Info] The 'ggtext' package is not installed. Colored facet labels will be disabled.")
    message("             To enable this feature, please run: install.packages('ggtext')")
  }

  # --- Loop through each unique LoopType ---
  unique_types <- unique(final_go_df$CleanLoopType)

  for (ltype in unique_types) {
    # 1. Filter data for current task
    sub_df <- final_go_df %>% dplyr::filter(CleanLoopType == ltype)
    if (nrow(sub_df) == 0) next

    # 2. Calc Stats
    sub_df$logP <- -log10(sub_df$pvalue)
    max_logp <- max(sub_df$logP, na.rm = TRUE)
    max_count <- max(sub_df$Count, na.rm = TRUE)
    scale_f <- max_logp / max_count * 1.0

    # 3. Factor Reordering
    sub_df <- sub_df %>%
      dplyr::group_by(ONTOLOGY) %>%
      dplyr::arrange(logP) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(Description_unique = factor(Description, levels = unique(Description)))

    onto_levels <- sort(unique(sub_df$ONTOLOGY))
    sub_df$ONTOLOGY <- factor(sub_df$ONTOLOGY, levels = onto_levels)
    n_onto <- length(onto_levels)
    onto_colors <- RColorBrewer::brewer.pal(max(3, n_onto), "Dark2")[seq_len(n_onto)]

    if (use_ggtext) {
      colored_labels <- paste0("<span style='color:", onto_colors, "'>", onto_levels, "</span>")
      sub_df$ONTOLOGY_Plot <- factor(sub_df$ONTOLOGY, levels = onto_levels, labels = colored_labels)
    } else {
      sub_df$ONTOLOGY_Plot <- sub_df$ONTOLOGY
    }

    # 4. Base Plot
    p_go <- ggplot2::ggplot(sub_df, ggplot2::aes(y = Description_unique)) +
      # Lollipop Stick & Head
      ggplot2::geom_segment(ggplot2::aes(x = 0, xend = logP, yend = Description_unique, color = ONTOLOGY), size = 3) +
      ggplot2::geom_point(ggplot2::aes(x = logP, color = ONTOLOGY), size = 5) +

      # Gene Count Line & Triangle
      ggplot2::geom_path(ggplot2::aes(x = Count * scale_f, group = 1), color = "grey60", size = 1.5, linetype = "11") +
      ggplot2::geom_point(ggplot2::aes(x = Count * scale_f), color = "grey60", size = 4, shape = 17) +

      # Axis
      ggplot2::scale_x_continuous(
        name = expression(-log[10](p - value)),
        expand = ggplot2::expansion(mult = c(0, 0.6)),
        sec.axis = ggplot2::sec_axis(~ . / scale_f, name = "Gene Counts")
      ) +
      ggplot2::facet_grid(ONTOLOGY_Plot ~ ., scales = "free_y", space = "free_y", switch = "y") +
      ggplot2::labs(
        y = NULL,
        title = paste0("GO Enrichment: ", ltype),
        subtitle = "Colored Dot: Significance | Grey Triangle: Gene Count"
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.text.y = ggplot2::element_text(size = 10, color = "black"),
        strip.placement = "outside",
        strip.background = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_line(color = "grey95", linetype = "dashed"),
        legend.position = "none"
      ) +
      ggplot2::scale_color_brewer(palette = "Dark2")

    if (use_ggtext) {
      p_go <- p_go + ggplot2::theme(strip.text.y.left = ggtext::element_markdown(angle = 0, face = "bold", size = 12))
    } else {
      p_go <- p_go + ggplot2::theme(strip.text.y.left = ggplot2::element_text(angle = 0, face = "bold", size = 12, color = onto_colors))
    }

    # 5. Save
    clean_ltype_name <- gsub("[^A-Za-z0-9]", "", ltype)
    out_name <- file.path(out_dir, paste0(base_project_name, "_GO_Lollipop_", clean_ltype_name, ".pdf"))

    n_terms <- nrow(sub_df)
    n_groups <- length(unique(sub_df$ONTOLOGY))
    plot_height <- max(5, (n_terms * 0.4) + (n_groups * 0.5))

    tryCatch(
      {
        ggplot2::ggsave(out_name, p_go, width = 9, height = plot_height)
        message("    Saved Separate GO Plot (Densest Dashed): ", out_name)
      },
      error = function(e) message("    Save Failed: ", e$message)
    )
  }
}


#' @title Generate Expression Heatmap and Connectivity Raincloud Plots
#' @description Plots Z-score normalized expression heatmaps for target genes. Also generates scatter and Raincloud plots correlating 3D connectivity (Degree) with expression levels and LFC, inheriting Hub classifications from upstream.
#' @param target_genes Character vector of gene symbols.
#' @param tpm_mat_raw Data frame. Raw TPM/FPKM expression matrix.
#' @param meta_raw Data frame. Sample metadata.
#' @param loop_stats_df Data frame. Promoter or Distal statistics containing node degree and Hub labels.
#' @param global_glist Named numeric vector. Global LFC values.
#' @param heatmap_ntop Integer. Max highly variable genes for the heatmap.
#' @param cor_method Character. Correlation method.
#' @param current_proj_name Character. Project prefix.
#' @param out_dir Character. Output directory.
#' @param source_type Character. Source type ("loops" or "targets") for plot subtitling.
#' @param target_col Character or NULL. Specific column in `loop_stats_df` to use as connectivity degree (e.g., "n_Linked_Distal").
#' @param skip_heatmap Logical. If `TRUE`, skips drawing the heatmap and only draws connectivity plots.
#' @return Invisible `NULL`. Saves multiple PDF plots.
#' @keywords internal
run_heatmap_and_connectivity <- function(target_genes, tpm_mat_raw, meta_raw, loop_stats_df, global_glist, heatmap_ntop, cor_method, current_proj_name, out_dir, source_type, target_col = NULL, skip_heatmap = FALSE) {
  req_pkgs <- c("ggpointdensity", "viridis", "ggpubr", "ggdist")
  for (pkg in req_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      warning("    Missing package '", pkg, "'. Please install it via install.packages('", pkg, "')")
      return()
    }
  }

  colnames(tpm_mat_raw) <- trimws(colnames(tpm_mat_raw))
  valid_s <- intersect(meta_raw$SampleID, colnames(tpm_mat_raw))
  if (length(valid_s) == 0) {
    return()
  }
  curr_mat <- tpm_mat_raw[, valid_s, drop = FALSE]
  curr_meta <- meta_raw %>% dplyr::filter(SampleID %in% valid_s)


  if (!skip_heatmap) {
    expr_genes <- intersect(target_genes, rownames(curr_mat))
    if (length(expr_genes) >= 5) {
      tryCatch(
        {
          mat_plot <- log2(curr_mat[expr_genes, , drop = FALSE] + 1)
          if (!is.null(heatmap_ntop) && nrow(mat_plot) > heatmap_ntop) {
            row_vars <- apply(mat_plot, 1, var, na.rm = TRUE)
            top_genes <- head(names(sort(row_vars, decreasing = TRUE)), heatmap_ntop)
            mat_plot <- mat_plot[top_genes, , drop = FALSE]
          }
          mat_scaled <- t(scale(t(mat_plot)))
          mat_scaled[mat_scaled > 2] <- 2
          mat_scaled[mat_scaled < -2] <- -2
          mat_scaled[is.na(mat_scaled)] <- 0

          if (requireNamespace("ComplexHeatmap", quietly = TRUE)) {
            col_fun <- circlize::colorRamp2(c(-2, 0, 2), c("#2BB2D1", "white", "#FF8181"))
            groups <- unique(curr_meta$Group)
            n_groups <- length(groups)
            if (n_groups <= 9) {
              cols <- RColorBrewer::brewer.pal(max(3, n_groups), "Set1")[seq_len(n_groups)]
            } else {
              cols <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(n_groups)
            }
            group_cols <- setNames(cols, groups)
            ha <- ComplexHeatmap::HeatmapAnnotation(Group = curr_meta$Group, col = list(Group = group_cols), simple_anno_size = unit(0.3, "cm"))
            hm <- ComplexHeatmap::Heatmap(
              mat_scaled,
              name = "Z-score", col = col_fun, cluster_columns = FALSE,
              show_row_names = (nrow(mat_scaled) <= 80), top_annotation = ha, border = TRUE,
              column_title = paste0("Expression Heatmap\n", current_proj_name), use_raster = FALSE
            )
            out_file_hm <- file.path(out_dir, paste0(current_proj_name, "_Expression_Heatmap.pdf"))
            pdf(out_file_hm, width = 7, height = 8)
            tryCatch({
              ComplexHeatmap::draw(hm)
            }, finally = {
              dev.off()
            })
            message("    Saved Heatmap: ", out_file_hm)
          }
        },
        error = function(e) warning("    [Heatmap] ", e$message)
      )
    }
  }


  if (is.null(loop_stats_df)) {
    return()
  }
  use_col <- NULL
  file_tag <- "Connectivity"
  display_tag <- "Total Loops"

  if (!is.null(target_col)) {
    if (target_col %in% colnames(loop_stats_df)) {
      use_col <- target_col
      file_tag <- paste0("Connectivity_", target_col)
      display_tag <- target_col
    } else {
      warning("    Requested column '", target_col, "' not found.")
      return()
    }
  } else {
    if ("Total_Loops" %in% colnames(loop_stats_df)) {
      use_col <- "Total_Loops"
    } else if ("Loop_Degree" %in% colnames(loop_stats_df)) {
      use_col <- "Loop_Degree"
    } else if ("degree" %in% colnames(loop_stats_df)) use_col <- "degree"
  }
  if (is.null(use_col)) {
    return()
  }

  gene_col_name <- colnames(loop_stats_df)[1]
  valid_targets <- intersect(target_genes, loop_stats_df[[gene_col_name]])
  if (length(valid_targets) < 5) {
    return()
  }

  cols_to_extract <- c(gene_col_name, use_col)
  if ("Is_High_Connectivity_Gene" %in% colnames(loop_stats_df)) cols_to_extract <- c(cols_to_extract, "Is_High_Connectivity_Gene")
  if ("Is_High_Distal_Connectivity_Gene" %in% colnames(loop_stats_df)) cols_to_extract <- c(cols_to_extract, "Is_High_Distal_Connectivity_Gene")
  if ("High_Connectivity_Gene" %in% colnames(loop_stats_df)) cols_to_extract <- c(cols_to_extract, "High_Connectivity_Gene")
  cols_to_extract <- unique(cols_to_extract)

  stats_subset <- loop_stats_df[loop_stats_df[[gene_col_name]] %in% valid_targets, cols_to_extract]

  colnames(stats_subset)[which(colnames(stats_subset) == use_col)] <- "Degree"
  colnames(stats_subset)[1] <- "Gene"

  valid_expr_targets <- intersect(stats_subset$Gene, rownames(curr_mat))
  if (length(valid_expr_targets) < 5) {
    return()
  }
  stats_subset <- stats_subset[stats_subset$Gene %in% valid_expr_targets, ]

  mean_expr <- rowMeans(curr_mat[stats_subset$Gene, , drop = FALSE], na.rm = TRUE)
  expr_vals <- log2(mean_expr + 1)
  lfc_vals <- global_glist[stats_subset$Gene]

  plot_df <- stats_subset %>%
    dplyr::mutate(
      Expression = as.numeric(expr_vals),
      LFC = as.numeric(lfc_vals),
      Log10Degree = log10(Degree)
    ) %>%
    dplyr::filter(!is.na(Expression), !is.na(LFC)) %>%
    dplyr::filter(Degree >= 1)

  if (nrow(plot_df) < 5) {
    return()
  }

  title_suffix <- if (source_type == "loops") "Looped (Specified Types)" else "Looped Targets"
  full_title_suffix <- paste0(title_suffix, " | ", display_tag)


  tryCatch(
    {
      p_scatter <- ggplot2::ggplot(plot_df, ggplot2::aes(x = Log10Degree, y = Expression))
      if (requireNamespace("ggpointdensity", quietly = TRUE) && requireNamespace("viridis", quietly = TRUE)) {
        p_scatter <- p_scatter + ggpointdensity::geom_pointdensity(alpha = 0.6, size = 1.5) + viridis::scale_color_viridis(option = "D", name = "Density")
      } else {
        p_scatter <- p_scatter + ggplot2::geom_point(alpha = 0.5, color = "#2c7bb6")
      }
      p_scatter <- p_scatter + ggplot2::geom_smooth(method = "lm", formula = y ~ x, color = "black", se = TRUE, size = 0.8)
      if (requireNamespace("ggpubr", quietly = TRUE)) {
        p_scatter <- p_scatter + ggpubr::stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top", size = 4)
      }
      p_scatter <- p_scatter + ggplot2::theme_bw() + ggplot2::theme(panel.grid = ggplot2::element_blank(), panel.background = ggplot2::element_rect(color = "black", fill = "transparent"), legend.key = ggplot2::element_rect(fill = "transparent"), plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")) +
        ggplot2::labs(title = paste0("Connectivity vs Expression (Scatter)"), subtitle = paste0(full_title_suffix, "\nGenes: ", nrow(plot_df)), x = paste0("Log10 (", display_tag, ")"), y = "Log2(Mean Expression + 1)")
      out_scatter <- file.path(out_dir, paste0(current_proj_name, "_", file_tag, "_Scatter.pdf"))
      ggplot2::ggsave(out_scatter, p_scatter, width = 6, height = 5)
      message("    Saved Scatter Plot (", display_tag, "): ", out_scatter)
    },
    error = function(e) message("    Scatter Plot Failed: ", e$message)
  )


  tryCatch(
    {
      has_distal <- "Is_High_Distal_Connectivity_Gene" %in% colnames(plot_df)
      has_total <- "Is_High_Connectivity_Gene" %in% colnames(plot_df)
      has_old_hc <- "High_Connectivity_Gene" %in% colnames(plot_df)

      if (has_distal && has_total) {
        df_distal <- plot_df %>%
          dplyr::filter(Is_High_Distal_Connectivity_Gene %in% c("Yes", "TRUE", TRUE, 1)) %>%
          dplyr::mutate(Conn_Group = "High Distal")
        df_total <- plot_df %>%
          dplyr::filter(Is_High_Connectivity_Gene %in% c("Yes", "TRUE", TRUE, 1)) %>%
          dplyr::mutate(Conn_Group = "High Total")
        df_others <- plot_df %>%
          dplyr::filter(!(Is_High_Distal_Connectivity_Gene %in% c("Yes", "TRUE", TRUE, 1)) & !(Is_High_Connectivity_Gene %in% c("Yes", "TRUE", TRUE, 1))) %>%
          dplyr::mutate(Conn_Group = "Others")
        plot_df_rc <- dplyr::bind_rows(df_distal, df_total, df_others)
        plot_df_rc$Conn_Group <- factor(plot_df_rc$Conn_Group, levels = c("High Distal", "High Total", "Others"))

        file_tag_rc <- "Combined_Connectivity"
        custom_colors <- c("High Distal" = "#9BC985", "High Total" = "#ECB884", "Others" = "#82969D")
      } else if (has_total) {
        plot_df_rc <- plot_df %>% dplyr::mutate(Conn_Group = ifelse(Is_High_Connectivity_Gene %in% c("Yes", "TRUE", TRUE, 1), "High Total", "Others"))
        plot_df_rc$Conn_Group <- factor(plot_df_rc$Conn_Group, levels = c("High Total", "Others"))
        file_tag_rc <- file_tag
        custom_colors <- c("High Total" = "#ECB884", "Others" = "#82969D")
      } else if (has_old_hc) {
        plot_df_rc <- plot_df %>% dplyr::mutate(Conn_Group = ifelse(High_Connectivity_Gene %in% c("Yes", "TRUE", TRUE, 1), "High Total", "Others"))
        plot_df_rc$Conn_Group <- factor(plot_df_rc$Conn_Group, levels = c("High Total", "Others"))
        file_tag_rc <- file_tag
        custom_colors <- c("High Total" = "#ECB884", "Others" = "#82969D")
      } else {
        message("    [Vis Info] Using dynamic threshold for grouping as upstream labels are missing.")
        deg_thresh <- quantile(plot_df$Degree, 0.75, na.rm = TRUE)
        if (deg_thresh < 2) deg_thresh <- 2
        plot_df_rc <- plot_df %>% dplyr::mutate(Conn_Group = ifelse(Degree >= deg_thresh, "High Total", "Others"))
        plot_df_rc$Conn_Group <- factor(plot_df_rc$Conn_Group, levels = c("High Total", "Others"))
        file_tag_rc <- file_tag
        custom_colors <- c("High Total" = "#ECB884", "Others" = "#82969D")
      }

      plot_df_rc <- plot_df_rc %>% dplyr::filter(!is.na(Conn_Group))
      plot_df_rc$Conn_Group <- droplevels(plot_df_rc$Conn_Group)

      if (nlevels(plot_df_rc$Conn_Group) > 1) {
        clean_theme <- ggplot2::theme_classic() +
          ggplot2::theme(
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
          if ("Others" %in% lvls) {
            res <- c()
            if ("High Distal" %in% lvls) {
              p <- tryCatch(wilcox.test(plot_df_rc[[val_col]][plot_df_rc$Conn_Group == "High Distal"], plot_df_rc[[val_col]][plot_df_rc$Conn_Group == "Others"])$p.value, error = function(e) NA)
              if (!is.na(p)) {
                s <- dplyr::case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "ns")
                res <- c(res, paste0("Distal=", signif(p, 3), " (", s, ")"))
              }
            }
            if ("High Total" %in% lvls) {
              p <- tryCatch(wilcox.test(plot_df_rc[[val_col]][plot_df_rc$Conn_Group == "High Total"], plot_df_rc[[val_col]][plot_df_rc$Conn_Group == "Others"])$p.value, error = function(e) NA)
              if (!is.na(p)) {
                s <- dplyr::case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "ns")
                res <- c(res, paste0("Total=", signif(p, 3), " (", s, ")"))
              }
            }
            if (length(res) > 0) {
              return(paste0("Wilcox P (vs Others):\n", paste(res, collapse = " | ")))
            }
            return("Wilcox P: NA")
          } else {
            return("Wilcox P: NA (No 'Others' group)")
          }
        }


        plot_width <- if (nlevels(plot_df_rc$Conn_Group) == 3) 4.2 else 3.2
        plot_height <- 5.5


        subtitle_text_lfc <- get_pval_str("LFC")

        p_box_lfc <- ggplot2::ggplot(plot_df_rc, ggplot2::aes(fill = Conn_Group)) +
          ggplot2::geom_jitter(
            ggplot2::aes(x = as.numeric(Conn_Group) - 0.12, y = LFC, color = Conn_Group),
            shape = 16, width = 0.03, height = 0, alpha = 0.6, size = 0.8
          ) +
          ggplot2::stat_boxplot(
            ggplot2::aes(x = as.numeric(Conn_Group), y = LFC, color = Conn_Group),
            geom = "errorbar", width = 0.05, linewidth = 0.5
          ) +
          ggplot2::geom_boxplot(
            ggplot2::aes(x = as.numeric(Conn_Group), y = LFC, color = Conn_Group),
            width = 0.12, notch = TRUE, outlier.shape = NA, alpha = 1, linewidth = 0.5
          ) +
          ggplot2::stat_summary(
            ggplot2::aes(x = as.numeric(Conn_Group), y = LFC),
            fun = median, fun.min = median, fun.max = median,
            geom = "crossbar", width = 0.1, color = "black", linewidth = 0.4
          ) +
          ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.6) +
          ggdist::stat_slab(
            ggplot2::aes(x = as.numeric(Conn_Group) + 0.07, y = LFC, fill = Conn_Group),
            adjust = 0.5, width = 0.35, justification = 0,
            alpha = 0.3, color = NA
          ) +
          ggdist::stat_slab(
            ggplot2::aes(x = as.numeric(Conn_Group) + 0.07, y = LFC, color = Conn_Group),
            adjust = 0.5, width = 0.35, justification = 0,
            fill = NA, alpha = 0.5, linewidth = 0.4
          ) +
          ggplot2::scale_x_continuous(
            breaks = seq_along(levels(plot_df_rc$Conn_Group)),
            labels = levels(plot_df_rc$Conn_Group)
          ) +
          ggplot2::coord_cartesian(xlim = c(0.75, length(levels(plot_df_rc$Conn_Group)) + 0.6)) +
          ggplot2::scale_fill_manual(values = custom_colors) +
          ggplot2::scale_color_manual(values = custom_colors) +
          ggplot2::labs(
            title = "Regulation: High_connectivity vs Others",
            subtitle = subtitle_text_lfc,
            x = NULL, y = "Log2 Fold Change (LFC)"
          ) +
          clean_theme

        out_box_lfc <- file.path(out_dir, paste0(current_proj_name, "_", file_tag_rc, "_Raincloud_LFC.pdf"))
        suppressWarnings(ggplot2::ggsave(out_box_lfc, p_box_lfc, width = plot_width, height = plot_height))
        message("    Saved Raincloud Plot (LFC): ", out_box_lfc)

        subtitle_text_expr <- get_pval_str("Expression")

        p_box_expr <- ggplot2::ggplot(plot_df_rc, ggplot2::aes(fill = Conn_Group)) +
          ggplot2::geom_jitter(
            ggplot2::aes(x = as.numeric(Conn_Group) - 0.12, y = Expression, color = Conn_Group),
            shape = 16, width = 0.03, height = 0, alpha = 0.6, size = 0.8
          ) +
          ggplot2::stat_boxplot(
            ggplot2::aes(x = as.numeric(Conn_Group), y = Expression, color = Conn_Group),
            geom = "errorbar", width = 0.05, linewidth = 0.5
          ) +
          ggplot2::geom_boxplot(
            ggplot2::aes(x = as.numeric(Conn_Group), y = Expression, color = Conn_Group),
            width = 0.12, notch = TRUE, outlier.shape = NA, alpha = 1, linewidth = 0.5
          ) +
          ggplot2::stat_summary(
            ggplot2::aes(x = as.numeric(Conn_Group), y = Expression),
            fun = median, fun.min = median, fun.max = median,
            geom = "crossbar", width = 0.1, color = "black", linewidth = 0.4
          ) +
          ggdist::stat_slab(
            ggplot2::aes(x = as.numeric(Conn_Group) + 0.07, y = Expression, fill = Conn_Group),
            adjust = 0.5, width = 0.35, justification = 0,
            alpha = 0.3, color = NA
          ) +
          ggdist::stat_slab(
            ggplot2::aes(x = as.numeric(Conn_Group) + 0.07, y = Expression, color = Conn_Group),
            adjust = 0.5, width = 0.35, justification = 0,
            fill = NA, alpha = 0.5, linewidth = 0.4
          ) +
          ggplot2::scale_x_continuous(
            breaks = seq_along(levels(plot_df_rc$Conn_Group)),
            labels = levels(plot_df_rc$Conn_Group)
          ) +
          ggplot2::coord_cartesian(xlim = c(0.75, length(levels(plot_df_rc$Conn_Group)) + 0.6)) +
          ggplot2::scale_fill_manual(values = custom_colors) +
          ggplot2::scale_color_manual(values = custom_colors) +
          ggplot2::labs(
            title = "Regulation: High_connectivity vs Others",
            subtitle = subtitle_text_expr,
            x = NULL, y = "Log2(Mean Expression + 1)"
          ) +
          clean_theme

        out_box_expr <- file.path(out_dir, paste0(current_proj_name, "_", file_tag_rc, "_Raincloud_Expr.pdf"))
        suppressWarnings(ggplot2::ggsave(out_box_expr, p_box_expr, width = plot_width, height = plot_height))
        message("    Saved Raincloud Plot (Expr): ", out_box_expr)
      } else {
        warning("   Not enough variance to split groups.")
      }
    },
    error = function(e) warning("    Raincloud Plot Failed: ", e$message)
  )
}


#' @title Run Dual Motif Analysis for Loop Anchors
#' @description Extracts genomic coordinates for proximal (gene-centric) and distal (enhancer-centric) anchors of loops connecting to target genes. Scans sequences against JASPAR core motifs to identify enriched TFBS.
#' @param target_genes Character vector of target gene symbols.
#' @param loop_df Data frame. Loop annotation table containing coordinates (chr1/start1/end1) and gene assignments.
#' @param genome_id Character. Genome assembly (e.g., "hg19", "mm10") for sequence extraction.
#' @param pval_thresh Numeric. P-value cutoff for motifmatchr scanning.
#' @param current_proj_name Character. Project prefix.
#' @param out_dir Character. Output directory.
#' @param top_n Integer. Number of top enriched motifs to output as SeqLogos.
#' @return Invisible `NULL`. Triggers internal plotting functions.
#' @keywords internal
run_distal_motif_analysis <- function(target_genes, loop_df, genome_id, pval_thresh, current_proj_name, out_dir, top_n = 5) {
  requireNamespace("motifmatchr")
  requireNamespace("JASPAR2020")
  requireNamespace("TFBSTools")
  requireNamespace("GenomicRanges")
  requireNamespace("ggseqlogo")
  bs_pkg <- switch(genome_id,
    "hg19" = "BSgenome.Hsapiens.UCSC.hg19",
    "hg38" = "BSgenome.Hsapiens.UCSC.hg38",
    "mm9"  = "BSgenome.Mmusculus.UCSC.mm9",
    "mm10" = "BSgenome.Mmusculus.UCSC.mm10"
  )
  species_id <- if (grepl("mm", genome_id)) 10090 else 9606
  if (!requireNamespace(bs_pkg, quietly = TRUE)) {
    stop("Please install ", bs_pkg, " via BiocManager::install('", bs_pkg, "')")
  }
  genome_obj <- get0(bs_pkg, envir = asNamespace(bs_pkg))
  if (!is.data.frame(loop_df)) loop_df <- as.data.frame(loop_df)

  possibilities_1 <- c("anchor1_gene", "Anchor1_Gene", "gene_name_1", "Gene_Name_1", "Symbol_1", "nearest_gene_1", "gene1")
  possibilities_2 <- c("anchor2_gene", "Anchor2_Gene", "gene_name_2", "Gene_Name_2", "Symbol_2", "nearest_gene_2", "gene2")
  col_g1 <- intersect(possibilities_1, colnames(loop_df))[1]
  col_g2 <- intersect(possibilities_2, colnames(loop_df))[1]
  if (is.na(col_g1) || is.na(col_g2)) {
    warning("    Could not find gene name columns.")
    return()
  }
  loop_df$gene_name_1 <- loop_df[[col_g1]]
  loop_df$gene_name_2 <- loop_df[[col_g2]]
  warning("     Using loop columns: ", col_g1, " & ", col_g2)

  idx_1 <- which(loop_df$gene_name_1 %in% target_genes)
  idx_2 <- which(loop_df$gene_name_2 %in% target_genes)
  target_indices <- unique(c(idx_1, idx_2))
  warning("    [Motif Debug] Found ", length(target_indices), " loops overlapping with target genes.")
  if (length(target_indices) < 5) {
    warning("     Too few loops associated with target genes (<5). Skipping.")
    return()
  }

  target_loops <- loop_df[target_indices, ]
  proximal_gr_list <- list()
  distal_gr_list <- list()

  for (i in seq_len(nrow(target_loops))) {
    row <- target_loops[i, ]
    g1 <- row$gene_name_1
    g2 <- row$gene_name_2
    is_a1_target <- g1 %in% target_genes
    is_a2_target <- g2 %in% target_genes
    if (is_a1_target) {
      proximal_gr_list[[length(proximal_gr_list) + 1]] <- GenomicRanges::GRanges(seqnames = row$chr1, ranges = IRanges::IRanges(start = row$start1, end = row$end1))
      distal_gr_list[[length(distal_gr_list) + 1]] <- GenomicRanges::GRanges(seqnames = row$chr2, ranges = IRanges::IRanges(start = row$start2, end = row$end2))
    }
    if (is_a2_target) {
      proximal_gr_list[[length(proximal_gr_list) + 1]] <- GenomicRanges::GRanges(seqnames = row$chr2, ranges = IRanges::IRanges(start = row$start2, end = row$end2))
      distal_gr_list[[length(distal_gr_list) + 1]] <- GenomicRanges::GRanges(seqnames = row$chr1, ranges = IRanges::IRanges(start = row$start1, end = row$end1))
    }
  }
  proximal_gr <- do.call(c, proximal_gr_list)
  distal_gr <- do.call(c, distal_gr_list)

  bg_indices <- setdiff(seq_len(nrow(loop_df)), target_indices)
  if (length(bg_indices) > 2000) bg_indices <- sample(bg_indices, 2000)
  bg_loops <- loop_df[bg_indices, ]
  bg_gr <- c(GenomicRanges::GRanges(seqnames = bg_loops$chr1, ranges = IRanges::IRanges(start = bg_loops$start1, end = bg_loops$end1)), GenomicRanges::GRanges(seqnames = bg_loops$chr2, ranges = IRanges::IRanges(start = bg_loops$start2, end = bg_loops$end2)))

  # --- Proximal Analysis ---
  message("    [Motif] analyzing Proximal (Target) Anchors (Species: ", species_id, ")...")
  res_prox <- .calc_motif_enrichment(proximal_gr, bg_gr, genome_obj, pval_thresh, species_id)
  res_prox <- .annotate_motif_families(res_prox)
  .plot_save_motif(res_prox, out_dir, paste0(current_proj_name, "_Motif_Proximal"))
  .plot_top_motif_logos(res_prox, out_dir, paste0(current_proj_name, "_Motif_Proximal"), top_n)
  .plot_motif_rank_scatter(res_prox, out_dir, paste0(current_proj_name, "_Motif_Proximal"))

  # --- Distal Analysis ---
  message("    [Motif] analyzing Distal (Opposite) Anchors...")
  res_dist <- .calc_motif_enrichment(distal_gr, bg_gr, genome_obj, pval_thresh, species_id)
  res_dist <- .annotate_motif_families(res_dist)
  .plot_save_motif(res_dist, out_dir, paste0(current_proj_name, "_Motif_Distal"))
  .plot_top_motif_logos(res_dist, out_dir, paste0(current_proj_name, "_Motif_Distal"), top_n)
  .plot_motif_rank_scatter(res_dist, out_dir, paste0(current_proj_name, "_Motif_Distal"))
}


#' @title Plot Motif Rank Scatter
#' @description Creates a scatter plot ranking motifs by FDR, where point size represents Odds Ratio and color represents the TF Family.
#' @param res_df Data frame of motif enrichment results.
#' @param out_dir Character. Output directory.
#' @param prefix Character. File prefix.
#' @param fdr_thresh Numeric. FDR threshold for significance coloration.
#' @return Invisible `NULL`.
#' @keywords internal
.plot_motif_rank_scatter <- function(res_df, out_dir, prefix, fdr_thresh = 0.05) {
  if (is.null(res_df) || nrow(res_df) == 0) {
    return()
  }

  req_pkgs <- c("ggplot2", "dplyr")
  for (pkg in req_pkgs) if (!requireNamespace(pkg, quietly = TRUE)) {
    return()
  }

  # Check if Family column exists
  if (!"Family" %in% colnames(res_df)) res_df$Family <- "Unknown"

  # 1. Prepare Data & Calculate Ranks
  plot_df <- res_df %>%
    dplyr::mutate(
      FDR = ifelse(is.na(FDR), 1, FDR),
      LogFDR = -log10(FDR + 1e-300),
      Is_Sig = FDR < fdr_thresh,
      OddsRatio = as.numeric(OddsRatio),
      OddsRatio = ifelse(is.infinite(OddsRatio) | OddsRatio > 1.6, 1.6, OddsRatio),
      OddsRatio = ifelse(is.infinite(OddsRatio) | OddsRatio < 0.8, 0.8, OddsRatio)
    ) %>%
    dplyr::arrange(desc(LogFDR)) %>%
    dplyr::mutate(Rank = dplyr::row_number())

  # 2. Process Families (Keep Top 10 frequent among significant ones)
  if (sum(plot_df$Is_Sig, na.rm = TRUE) == 0) {
    plot_df$PlotFamily <- "Not Significant"
  } else {
    fam_counts <- table(plot_df$Family[plot_df$Is_Sig & !plot_df$Family %in% c("Unknown", "", NA)])
    if (length(fam_counts) > 0) {
      top_fams <- names(sort(fam_counts, decreasing = TRUE))[seq_len(min(10, length(fam_counts)))]
      plot_df$PlotFamily <- ifelse(plot_df$Family %in% top_fams, plot_df$Family, "Others")
    } else {
      plot_df$PlotFamily <- "Others"
    }
    plot_df$PlotFamily <- ifelse(plot_df$Is_Sig, plot_df$PlotFamily, "Not Significant")
  }

  # 3. Order factors for plotting
  all_plot_fams <- unique(plot_df$PlotFamily)
  sig_fams <- setdiff(all_plot_fams, c("Others", "Not Significant"))

  factor_levels <- sig_fams
  if ("Others" %in% all_plot_fams) factor_levels <- c(factor_levels, "Others")
  if ("Not Significant" %in% all_plot_fams) factor_levels <- c(factor_levels, "Not Significant")

  plot_df$PlotFamily <- factor(plot_df$PlotFamily, levels = factor_levels)

  # 4. Color Palette (Custom distinct colors excluding grey)
  distinct_cols <- c(
    "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
    "#FFFF33", "#A65628", "#F781BF", "#1B9E77", "#D95F02"
  )

  color_map <- c()
  if (length(sig_fams) > 0) {
    color_map <- setNames(distinct_cols[seq_len(length(sig_fams))], sig_fams)
  }

  if ("Others" %in% all_plot_fams) {
    color_map["Others"] <- "black"
  }
  if ("Not Significant" %in% all_plot_fams) {
    color_map["Not Significant"] <- "grey85"
  }

  # 5. Plotting
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = Rank, y = LogFDR, size = OddsRatio, color = PlotFamily)) +
    ggplot2::geom_point(alpha = 0.8) +
    ggplot2::scale_color_manual(values = color_map, name = "TF Family (Top 10)") +
    ggplot2::scale_radius(name = "Odds Ratio", range = c(1, 7), breaks = c(0.8, 1.0, 1.2, 1.4, 1.6)) +
    ggplot2::geom_hline(yintercept = -log10(fdr_thresh), linetype = "dashed", color = "black", alpha = 0.5) +
    ggplot2::labs(
      title = paste0("Motif Enrichment Rank: ", basename(prefix)),
      x = "Rank",
      y = "-log10(FDR)"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position = "right",
      legend.text = ggplot2::element_text(size = 8),
      legend.title = ggplot2::element_text(size = 9, face = "bold")
    )

  out_pdf <- file.path(out_dir, paste0(prefix, "_RankScatter.pdf"))
  ggplot2::ggsave(out_pdf, p, width = 8, height = 6)
  message("      Saved Rank Scatter Plot: ", out_pdf)
}


#' @title Calculate Motif Enrichment via Fisher's Exact Test
#' @description Compares motif hit frequencies between foreground (anchor regions) and background (random loops) sequences.
#' @param fg_gr GRanges object. Foreground genomic regions.
#' @param bg_gr GRanges object. Background genomic regions.
#' @param genome_obj BSgenome object. Reference genome sequence.
#' @param pval_thresh Numeric. Cutoff for calling a motif match.
#' @param species_id Numeric. Taxonomy ID for JASPAR filtering (e.g., 9606).
#' @return Data frame of motif enrichment statistics (Pvalue, FDR, OddsRatio).
#' @keywords internal
.calc_motif_enrichment <- function(fg_gr, bg_gr, genome_obj, pval_thresh, species_id) {
  fg_gr <- GenomicRanges::resize(fg_gr, width = 500, fix = "center")
  bg_gr <- GenomicRanges::resize(bg_gr, width = 500, fix = "center")
  fg_gr <- fg_gr[GenomicRanges::start(fg_gr) > 0]
  bg_gr <- bg_gr[GenomicRanges::start(bg_gr) > 0]
  if (length(fg_gr) == 0 || length(bg_gr) == 0) {
    return(NULL)
  }
  fg_seq <- BSgenome::getSeq(genome_obj, fg_gr)
  bg_seq <- BSgenome::getSeq(genome_obj, bg_gr)

  opts <- list(species = species_id, collection = "CORE")
  pfm_list <- tryCatch(
    {
      TFBSTools::getMatrixSet(JASPAR2020::JASPAR2020, opts)
    },
    error = function(e) {
      TFBSTools::getMatrixSet(JASPAR2020::JASPAR2020, list(collection = "CORE"))
    }
  )
  if (length(pfm_list) == 0) {
    return(NULL)
  }

  fg_hits <- motifmatchr::matchMotifs(pfm_list, fg_seq, out = "matches", p.cutoff = pval_thresh)
  bg_hits <- motifmatchr::matchMotifs(pfm_list, bg_seq, out = "matches", p.cutoff = pval_thresh)

  fg_mat <- as.matrix(motifmatchr::motifMatches(fg_hits))
  bg_mat <- as.matrix(motifmatchr::motifMatches(bg_hits))
  fg_counts <- colSums(fg_mat)
  bg_counts <- colSums(bg_mat)
  n_fg <- length(fg_seq)
  n_bg <- length(bg_seq)

  results_list <- list()
  all_motifs <- union(names(fg_counts), names(bg_counts))
  for (motif_name in all_motifs) {
    a <- if (motif_name %in% names(fg_counts)) fg_counts[[motif_name]] else 0
    b <- if (motif_name %in% names(bg_counts)) bg_counts[[motif_name]] else 0
    c <- n_fg - a
    d <- n_bg - b
    if (a > 0) {
      ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")
      results_list[[motif_name]] <- data.frame(MotifID = motif_name, MotifName = TFBSTools::name(pfm_list[[motif_name]]), Pvalue = ft$p.value, OddsRatio = ft$estimate, FG_Hits = a, FG_Total = n_fg, BG_Hits = b, BG_Total = n_bg)
    }
  }
  if (length(results_list) == 0) {
    return(NULL)
  }
  res_df <- do.call(rbind, results_list)
  if (!is.null(res_df) && nrow(res_df) > 0) {
    res_df$FDR <- p.adjust(res_df$Pvalue, method = "BH")
    res_df <- res_df[order(res_df$Pvalue), ]
  }
  return(res_df)
}


#' @title Plot and Save Motif Results (Barplot)
#' @description Plots the top 15 enriched motifs ranked by -log10(P-value).
#' @param res_df Data frame of motif enrichment results.
#' @param out_dir Character. Output directory.
#' @param prefix Character. File prefix.
#' @return Invisible `NULL`.
#' @keywords internal
.plot_save_motif <- function(res_df, out_dir, prefix) {
  if (is.null(res_df) || nrow(res_df) == 0) {
    return()
  }
  csv_file <- file.path(out_dir, paste0(prefix, ".csv"))
  write.csv(res_df, csv_file, row.names = FALSE)
  top_df <- head(res_df, 15)
  top_df$MotifLabel <- paste0(top_df$MotifName, " (", top_df$MotifID, ")")
  top_df$LogP <- -log10(top_df$Pvalue)
  top_df$MotifLabel <- factor(top_df$MotifLabel, levels = rev(top_df$MotifLabel))
  p <- ggplot2::ggplot(top_df, ggplot2::aes(x = LogP, y = MotifLabel)) +
    ggplot2::geom_col(fill = "#E7298A", width = 0.7) +
    ggplot2::labs(title = paste0("Motif Enrichment: ", basename(prefix)), x = "-log10(P-value)", y = NULL) +
    ggplot2::theme_classic()
  pdf_file <- file.path(out_dir, paste0(prefix, ".pdf"))
  ggplot2::ggsave(pdf_file, p, width = 6, height = 6)
  message("      Saved Motif Plot: ", pdf_file)
}


#' @title Plot Top Motif Sequence Logos
#' @description Retrieves Position Frequency Matrices (PFM) from JASPAR and plots sequence logos for the top N enriched motifs.
#' @param res_df Data frame of motif enrichment results.
#' @param out_dir Character. Output directory.
#' @param prefix Character. File prefix.
#' @param top_n Integer. Number of top motifs to plot.
#' @return Invisible `NULL`.
#' @keywords internal
.plot_top_motif_logos <- function(res_df, out_dir, prefix, top_n) {
  if (is.null(res_df) || nrow(res_df) == 0) {
    return()
  }

  # Filter top N by P-value
  top_df <- head(res_df[order(res_df$Pvalue), ], top_n)
  motif_ids <- top_df$MotifID
  motif_names <- top_df$MotifName

  # Fetch matrices from JASPAR
  pfm_list <- tryCatch(
    {
      TFBSTools::getMatrixSet(JASPAR2020::JASPAR2020, opts = list(ID = motif_ids))
    },
    error = function(e) {
      return(NULL)
    }
  )

  if (is.null(pfm_list) || length(pfm_list) == 0) {
    return()
  }

  # Convert to PCM/ICM list for plotting
  plot_list <- list()
  for (i in seq_along(motif_ids)) {
    mid <- motif_ids[i]
    if (mid %in% names(pfm_list)) {
      mat <- TFBSTools::Matrix(pfm_list[[mid]])
      label <- paste0(motif_names[i], " (", mid, ")")
      plot_list[[label]] <- mat
    }
  }

  if (length(plot_list) == 0) {
    return()
  }

  # Plot with ggseqlogo
  p <- ggseqlogo::ggseqlogo(plot_list, ncol = 1) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = 10, face = "bold", hjust = 0),
      strip.background = ggplot2::element_rect(fill = "grey95", color = NA)
    ) +
    ggplot2::labs(y = "Bits", title = paste0("Top ", top_n, " Enriched Motifs (SeqLogo)"))

  out_pdf <- file.path(out_dir, paste0(prefix, "_Logos.pdf"))
  ggplot2::ggsave(out_pdf, p, width = 6, height = 2 * length(plot_list), limitsize = FALSE)
  message("      Saved Logo Plot: ", out_pdf)
}


#' @title Annotate Motif Families
#' @description Appends JASPAR TF family annotations to the motif enrichment results.
#' @param res_df Data frame. Motif enrichment results from `.calc_motif_enrichment`.
#' @return Data frame with an appended `Family` column.
#' @keywords internal
.annotate_motif_families <- function(res_df) {
  if (is.null(res_df) || nrow(res_df) == 0) {
    return(res_df)
  }

  requireNamespace("JASPAR2020")
  requireNamespace("TFBSTools")

  message("      [Motif Info] Annotating motif families (Family Only)...")

  tryCatch(
    {
      # 1. Fetch all matrix objects
      db_motifs <- TFBSTools::getMatrixSet(JASPAR2020::JASPAR2020, list(collection = "CORE"))

      # 2. Extract metadata (Using tags$family which worked previously)
      meta_list <- lapply(db_motifs, function(x) {
        tags <- TFBSTools::tags(x)

        # Extract Family only
        fam <- tags$family
        if (is.null(fam)) fam <- "Unknown"

        data.frame(
          MotifID = TFBSTools::ID(x),
          Family = paste(fam, collapse = "; "),
          stringsAsFactors = FALSE
        )
      })
      meta_df <- do.call(rbind, meta_list)

      # 3. Merge with result
      # Remove old Family/Class columns if they exist
      res_df <- res_df[, !colnames(res_df) %in% c("Family", "Class")]

      # Left join
      res_df <- merge(res_df, meta_df, by = "MotifID", all.x = TRUE)

      # 4. Reorder Columns (No Class)
      cols <- c("MotifID", "MotifName", "Family", "Pvalue", "FDR", "OddsRatio", "FG_Hits", "FG_Total", "BG_Hits", "BG_Total")
      extra_cols <- setdiff(colnames(res_df), cols)
      res_df <- res_df[, c(cols, extra_cols)]

      res_df <- res_df[order(res_df$Pvalue), ]
    },
    error = function(e) {
      warning("      [Motif] Family annotation failed: ", e$message)
      return(res_df)
    }
  )

  return(res_df)
}
