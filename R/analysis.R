#' Comprehensive Functional Profiling of Loop-Associated Target Genes
#'
#' Performs an integrated downstream analysis to characterize the functional and expression
#' landscapes of genes associated with specific chromatin loops.
#'
#' @description
#' This function acts as an automated pipeline that bridges 3D genomic interactions (loops)
#' with transcriptomic data (RNA-seq). It extracts target genes linked to loops (e.g.,
#' Enhancer-Promoter or Promoter-Promoter) and performs a comprehensive multi-omics
#' profiling suite, including:
#' \itemize{
#'   \item \strong{Differential Expression Analysis:} Compares the Log2 Fold Change (LFC) of loop-target genes against the genomic background to assess regulatory impact.
#'   \item \strong{Functional Enrichment:} Performs Gene Ontology (GO) ORA and Gene Set Enrichment Analysis (GSEA) to identify biological pathways.
#'   \item \strong{Expression Clustering:} Generates TPM heatmaps and sample correlation matrices to visualize expression patterns across groups.
#'   \item \strong{Interaction Networks:} (Optional) Constructs Protein-Protein Interaction (PPI) networks using the STRING database.
#' }
#'
#' @details
#' The analysis workflow proceeds in the following steps:
#' \enumerate{
#'   \item \strong{Target Extraction:} Identifies target genes based on the provided \code{annotation_res}, filtering by loop type (e.g., "EP") and method ("overlap" or "cluster").
#'   \item \strong{LFC Profiling:} Checks if target genes show significant global shifts in expression (e.g., are EP-loop targets generally upregulated?).
#'   \item \strong{Enrichment Analysis:} Uses \code{clusterProfiler} to run GO and GSEA. GSEA requires a ranked gene list derived from the provided \code{diff_file}.
#'   \item \strong{Visualization:} Produces PDF plots for all analyses in the specified \code{out_dir}.
#' }
#' Robust error handling (`tryCatch`) is implemented to ensure that a failure in one module (e.g., PPI) does not stop the entire pipeline.
#'
#' @param annotation_res List. The result object returned by \code{\link{annotate_peaks_and_loops}}.
#' @param use_filtered Logical. If \code{TRUE}, uses the refined/filtered loops; otherwise uses raw loops.
#' @param diff_file Character. Path to the differential expression file (CSV/TSV). Must contain gene IDs (rownames) and a Log Fold Change column.
#' @param lfc_col Character. The column name in \code{diff_file} representing Log2 Fold Change (e.g., "log2FoldChange").
#' @param tpm_file Character. Path to the TPM expression matrix (Genes x Samples).
#' @param metadata_file Character. Path to the sample metadata file. Must contain "SampleID" and "Group" columns.
#' @param target_source Character vector. Source of target genes:
#'   \itemize{
#'     \item \code{"overlap"}: Genes physically overlapping loop anchors.
#'     \item \code{"cluster"}: Genes associated with loop clusters (network awareness).
#'   }
#' @param loop_types Character vector. The types of loops to analyze (e.g., \code{c("EP", "PP")}).
#' @param include_direct Logical. Whether to include genes directly overlapping peaks (even without loops) in the analysis.
#' @param group_order Character vector. Optional factor levels to order the groups in plots (e.g., \code{c("Ctrl", "Treat")}).
#' @param project_name Character. Prefix for all output files.
#' @param out_dir Character. Directory to save output plots and tables.
#' @param org_db Character. Annotation database for enrichment (e.g., "org.Hs.eg.db" for human, "org.Mm.eg.db" for mouse).
#' @param run_go Logical. Whether to run Gene Ontology enrichment.
#' @param run_ppi Logical. Whether to run PPI network analysis (requires internet).
#' @param ppi_score Numeric. Minimum confidence score for STRING PPI interactions (default: 400).
#' @param ppi_ntop Numeric. Max number of top LFC genes to include in the PPI network.
#' @param heatmap_ntop Numeric. Max number of variable genes to plot in the heatmap.
#' @param gsea_ntop Numeric. Max number of genes to use for GSEA visualization.
#' @param cnet_ntop Numeric. Number of top enriched terms to show in the GO concept network plot.
#' @param stat_test Character. Statistical test for LFC comparison ("t.test" or "wilcox.test").
#' @param cor_method Character. Correlation method for sample clustering ("pearson", "spearman").
#'
#' @return An invisible list containing:
#' \itemize{
#'   \item \code{go_results}: A list of top GO enrichment tables for each analyzed loop type.
#'   \item \code{target_gene_sets}: The specific lists of gene IDs used for each analysis task.
#' }
#' @export
#' @examples
#' # =========================================================================
#' # 1. Prepare Mock Input Files (Simulating user data)
#' # =========================================================================
#' tmp_dir <- tempdir()
#'
#' # A. Mock TPM Matrix (Genes x Samples)
#' tpm_df <- data.frame(
#'   GeneID = c("TP53", "EGFR", "MYC", "GAPDH", "ACTB", "CD4"),
#'   S1 = c(10, 20, 5, 100, 50, 0),
#'   S2 = c(12, 18, 6, 105, 52, 0),
#'   S3 = c(50, 5, 50, 98, 51, 20),
#'   S4 = c(55, 4, 52, 100, 53, 22)
#' )
#' f_tpm <- file.path(tmp_dir, "mock_tpm.csv")
#' write.csv(tpm_df, f_tpm, row.names = FALSE)
#'
#' # B. Mock Differential Expression Table
#' diff_df <- data.frame(
#'   GeneID = c("TP53", "EGFR", "MYC", "GAPDH", "ACTB", "CD4"),
#'   log2FoldChange = c(2.5, -2.0, 3.1, 0.1, 0.05, 5.0),
#'   pvalue = c(0.001, 0.001, 0.001, 0.9, 0.8, 0.0001)
#' )
#' f_diff <- file.path(tmp_dir, "mock_diff.csv")
#' write.csv(diff_df, f_diff, row.names = FALSE)
#'
#' # C. Mock Metadata
#' meta_df <- data.frame(
#'   SampleID = c("S1", "S2", "S3", "S4"),
#'   Group = c("Ctrl", "Ctrl", "Treat", "Treat")
#' )
#' f_meta <- file.path(tmp_dir, "mock_meta.csv")
#' write.csv(meta_df, f_meta, row.names = FALSE)
#'
#' # =========================================================================
#' # 2. Prepare Mock Annotation Result (KITCHEN SINK VERSION)
#' # =========================================================================
#' # Create a robust dataframe with ALL column variations
#' loop_data <- data.frame(
#'   loop_id = c("L1", "L2", "L3"),
#'   # Support both "loop_type" and "type"
#'   loop_type = c("EP", "EP", "PP"),
#'   type = c("EP", "EP", "PP"),
#'
#'   # Support all naming conventions for gene 1
#'   gene_name_1 = c("TP53", "MYC", "GAPDH"),
#'   gene_id_1 = c("TP53", "MYC", "GAPDH"),
#'   anchor1_gene = c("TP53", "MYC", "GAPDH"),
#'
#'   # Support all naming conventions for gene 2
#'   gene_name_2 = c("Enhancer1", "Enhancer2", "Promoter2"),
#'   gene_id_2 = c("Enhancer1", "Enhancer2", "Promoter2"),
#'   anchor2_gene = c("Enhancer1", "Enhancer2", "Promoter2"),
#'   cluster_id = c("C1", "C2", "C3"),
#'   stringsAsFactors = FALSE
#' )
#'
#' # [CRITICAL FIX] Populate BOTH 'loop_annotation' (Raw) and 'filtered_loops' (Filtered)
#' # This ensures extract_target_gene_sets() finds data regardless of its default usage.
#' mock_res <- list(
#'   loop_annotation = loop_data, # For use_filtered = FALSE
#'   filtered_loops = loop_data, # For use_filtered = TRUE (or default)
#'   cluster_annotation = data.frame(
#'     cluster_id = c("C1", "C2", "C3"),
#'     expressed_genes_Total = c("TP53", "MYC", "GAPDH"),
#'     stringsAsFactors = FALSE
#'   )
#' )
#'
#' # =========================================================================
#' # 3. Run the Profiling Function
#' # =========================================================================
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   res <- profile_target_genes(
#'     annotation_res = mock_res,
#'     use_filtered = FALSE,
#'     diff_file = f_diff,
#'     lfc_col = "log2FoldChange",
#'     tpm_file = f_tpm,
#'     metadata_file = f_meta,
#'     target_source = "overlap",
#'     loop_types = c("EP"),
#'     project_name = "Example_Run",
#'     out_dir = file.path(tmp_dir, "output"),
#'     run_go = FALSE,
#'     run_ppi = FALSE,
#'     stat_test = "t.test"
#'   )
#'
#'   if (!is.null(res)) {
#'     message("✅ Results generated successfully!")
#'     print(names(res))
#'     # print(list.files(file.path(tmp_dir, "output"))) # Optional: check files
#'   }
#' }
profile_target_genes <- function(
  annotation_res,
  use_filtered = TRUE,
  diff_file,
  lfc_col,
  tpm_file,
  metadata_file,
  target_source = c("overlap", "cluster"),
  loop_types = c("EP", "PP"),
  include_direct = FALSE,
  group_order = NULL,
  project_name = "Analysis",
  out_dir = "./",
  # --- Params ---
  org_db = "org.Hs.eg.db",
  run_go = TRUE,
  run_ppi = FALSE,
  ppi_score = 400,
  ppi_ntop = 400,
  heatmap_ntop = 1000,
  gsea_ntop = 99999,
  cnet_ntop = 50,
  stat_test = "wilcox.test",
  cor_method = "pearson"
) {
  target_source <- match.arg(target_source, several.ok = TRUE)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  suffix <- if (use_filtered) "_Filtered" else "_Raw"
  base_project_name <- paste0(project_name, suffix)
  message(">>> Analysis Init | Base Project: ", base_project_name)

  # === 0. Dependency Check ===
  if (!requireNamespace("RColorBrewer", quietly = TRUE)) stop("Please install 'RColorBrewer'")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Please install 'ggplot2'")

  # Load data
  message("--- Reading files...")
  # [CRITICAL FIX] Ensure min_cols=1 is used for diff and tpm to avoid 'insufficient columns' error
  diff_df_raw <- read_robust_general(diff_file, header = TRUE, row_name = 1, desc = "Diff Table", min_cols = 1)
  tpm_mat_raw <- read_robust_general(tpm_file, header = TRUE, row_name = 1, desc = "TPM Matrix", min_cols = 1)
  meta_raw <- read_robust_general(metadata_file, header = TRUE, row_name = NULL, desc = "Metadata", min_cols = 1)

  # Check duplicates
  .check_duplicate_rownames <- function(mat, name) {
    rn <- rownames(mat)
    if (is.null(rn)) {
      return()
    }
    if (anyDuplicated(rn)) {
      dup_genes <- unique(rn[duplicated(rn) | duplicated(rn, fromLast = TRUE)])
      stop(sprintf("❌ %s contains duplicate rownames (e.g., %s).", name, paste(head(dup_genes, 3), collapse = ", ")))
    }
  }
  .check_duplicate_rownames(diff_df_raw, "Differential expression file")
  .check_duplicate_rownames(tpm_mat_raw, "TPM matrix")

  # Process metadata
  if (ncol(meta_raw) >= 2) colnames(meta_raw)[c(1, 2)] <- c("SampleID", "Group")
  meta_raw$SampleID <- trimws(as.character(meta_raw$SampleID))
  if (!is.null(group_order)) meta_raw$Group <- factor(meta_raw$Group, levels = group_order)
  meta_raw <- meta_raw %>% dplyr::arrange(Group)

  # Global LFC list
  if (!lfc_col %in% colnames(diff_df_raw)) stop("❌ LFC column not found")
  clean_diff <- diff_df_raw[!is.na(diff_df_raw[[lfc_col]]) & is.finite(diff_df_raw[[lfc_col]]), ]
  global_glist <- sort(setNames(clean_diff[[lfc_col]], rownames(clean_diff)), decreasing = TRUE)

  # Step 1: Extract target gene sets
  message("---Extracting target gene sets...")
  raw_gene_sets <- extract_target_gene_sets(annotation_res, target_source, loop_types, include_direct)

  if (length(raw_gene_sets) == 0) {
    warning("raw_gene_sets is empty! Skipping analysis.")
    return(NULL)
  }

  # Step 2: Build analysis queue
  analysis_queue <- build_analysis_queue(raw_gene_sets, loop_types)
  if (length(analysis_queue) == 0) {
    warning("analysis_queue is empty!")
    return(NULL)
  }

  all_go_results <- list()

  # Step 3: Run per-task analysis
  for (task_name in names(analysis_queue)) {
    target_genes <- analysis_queue[[task_name]]
    current_proj_name <- paste0(base_project_name, "_", task_name)

    message("\n============== Analysis: ", task_name, " (Genes: ", length(target_genes), ") ==============")

    if (length(target_genes) < 3) {
      message(" Too few genes (<3), skipping.")
      next
    }

    # === 3.1 LFC Boxplot (TryCatch) ===
    tryCatch(
      {
        run_lfc_boxplot(target_genes, global_glist, stat_test, current_proj_name, out_dir)
      },
      error = function(e) {
        warning("  LFC Boxplot failed for ", task_name, ": ", e$message)
      }
    )

    # === 3.2 GSEA (TryCatch) ===
    tryCatch(
      {
        # [REMOVED] seed argument
        run_gsea_analysis(target_genes, global_glist, gsea_ntop, current_proj_name, out_dir)
      },
      error = function(e) {
        warning(" GSEA failed for ", task_name, " (likely no significant pathways): ", e$message)
      }
    )

    # === 3.3 Heatmap (TryCatch) ===
    tryCatch(
      {
        run_heatmap_and_correlation(target_genes, tpm_mat_raw, meta_raw, heatmap_ntop, cor_method, current_proj_name, out_dir)
      },
      error = function(e) {
        warning(" Heatmap failed for ", task_name, ": ", e$message)
      }
    )

    # === 3.4 GO Enrichment (TryCatch) ===
    tryCatch(
      {
        if (run_go) {
          go_result <- run_go_enrichment(target_genes, org_db, global_glist, cnet_ntop, current_proj_name, out_dir)
          if (!is.null(go_result)) {
            write.csv(go_result, file.path(out_dir, paste0(current_proj_name, "_GO_Enrichment_All.csv")), row.names = FALSE)
            top_go <- head(go_result[order(go_result$pvalue), ], 10)
            if (nrow(top_go) > 0) {
              top_go$LoopType <- task_name
              all_go_results[[task_name]] <- top_go
            }
          }
        }
      },
      error = function(e) {
        warning(" GO Analysis failed for ", task_name, ": ", e$message)
      }
    )

    # === 3.5 PPI (TryCatch) ===
    tryCatch(
      {
        if (run_ppi) {
          run_ppi_analysis(target_genes, global_glist, org_db, ppi_score, ppi_ntop, current_proj_name, out_dir)
        }
      },
      error = function(e) {
        warning("  PPI Analysis failed for ", task_name, ": ", e$message)
      }
    )
  }

  # Step 4: Summary GO Plot
  tryCatch(
    {
      if (run_go && length(all_go_results) > 0) {
        plot_summary_go_lollipop(all_go_results, base_project_name, out_dir)
      }
    },
    error = function(e) {
      warning(" Summary Plot failed: ", e$message)
    }
  )

  message("\n All analysis complete.")
  invisible(list(go_results = all_go_results, target_gene_sets = analysis_queue))
}


#' Internal: Robust Data Reader
#'
#' Safely reads standard genomic formats (BED, BEDPE, CSV, TSV) using data.table (if avail)
#' or base R. Automatically handles delimiters and checks column counts.
#'
#' @param f Character. Path to input file.
#' @param header Logical. Whether the file contains a header. Default: FALSE (standard for BED/BEDPE).
#' @param row_name Either \code{NULL} or integer column index for row names.
#' @param desc Character. Description for error messages.
#' @param min_cols Integer. Minimum required columns.
#'
#' @return A data frame.
#'
#' @importFrom utils read.table
#' @keywords internal
read_robust_general <- function(f, header = FALSE, row_name = NULL, desc = "file", min_cols = 3) {
  # 1. 检查文件是否存在
  if (is.null(f) || length(f) == 0 || f == "") stop(esc, " path is empty.")
  if (!file.exists(f)) stop(desc, " not found: ", f)

  d <- NULL
  rn <- if (is.null(row_name)) NULL else 1

  # 2. 优先尝试 data.table::fread (速度快)
  if (requireNamespace("data.table", quietly = TRUE)) {
    try(
      {
        # fread 自动检测分隔符
        d_dt <- data.table::fread(f, header = header, data.table = FALSE, showProgress = FALSE)

        # 手动处理行名 (fread 不直接支持 row.names 参数)
        if (!is.null(row_name) && ncol(d_dt) > 1) {
          rownames(d_dt) <- d_dt[, 1]
          d_dt <- d_dt[, -1, drop = FALSE]
        }

        # 检查列数是否达标
        if (ncol(d_dt) >= min_cols) d <- d_dt
      },
      silent = TRUE
    )
  }

  # 3. 降级方案 A: 制表符 (TSV/BED)
  if (is.null(d)) {
    d <- tryCatch(
      utils::read.table(f,
        header = header, sep = "\t", row.names = rn,
        check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = ""
      ),
      error = function(e) NULL
    )
  }

  # 4. 降级方案 B: 逗号 (CSV)
  if (is.null(d) || ncol(d) < min_cols) {
    d <- tryCatch(
      utils::read.table(f,
        header = header, sep = ",", row.names = rn,
        check.names = FALSE, stringsAsFactors = FALSE, quote = "\"", comment.char = ""
      ),
      error = function(e) NULL
    )
  }

  # 5. 降级方案 C: 任意空白符
  if (is.null(d) || ncol(d) < min_cols) {
    d <- tryCatch(
      utils::read.table(f,
        header = header, sep = "", row.names = rn,
        check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = ""
      ),
      error = function(e) NULL
    )
  }

  # 6. 最终检查
  if (is.null(d)) {
    stop("failed to read ", desc, ". Format not recognized.")
  }

  if (ncol(d) < min_cols) {
    stop("❌ ", desc, " has insufficient columns (found ", ncol(d), ", required ", min_cols, ").")
  }

  return(d)
}

#' Internal: Extract Target Gene Sets from Annotation Results
#'
#' Parses loop and cluster annotations to extract target genes based on specified sources
#' ("overlap" from loop-gene overlap, "cluster" from peak-associated expressed genes).
#'
#' @param annotation_res List output from upstream annotation functions, expected to contain
#'   \code{loop_annotation} and/or \code{cluster_annotation}.
#' @param target_source Character vector: "overlap", "cluster", or both.
#' @param loop_types Character vector of loop types to consider (e.g., "EP", "PP").
#' @param include_direct Logical. If \code{TRUE} and "cluster" is in \code{target_source},
#'   include direct peak-linked genes under key "Direct".
#'
#' @return Named list of gene vectors, keyed by loop type or "Direct".
#'
#' @importFrom dplyr filter pull
#' @importFrom tidyr separate_rows
#' @keywords internal
extract_target_gene_sets <- function(annotation_res, target_source, loop_types, include_direct) {
  raw_gene_sets <- list()
  if ("overlap" %in% target_source && !is.null(annotation_res$loop_annotation)) {
    loop_df <- annotation_res$loop_annotation
    if (all(c("loop_type", "loop_genes") %in% colnames(loop_df))) {
      for (lt in loop_types) {
        sub_df <- loop_df %>%
          dplyr::filter(loop_type == lt) %>%
          dplyr::filter(!is.na(loop_genes) & loop_genes != "")
        if (nrow(sub_df) > 0) {
          gs <- unique(trimws(unlist(strsplit(as.character(sub_df$loop_genes), ";"))))
          gs <- gs[gs != ""]
          if (length(gs) > 0) raw_gene_sets[[lt]] <- gs
        }
      }
    }
  }
  if ("cluster" %in% target_source && include_direct && !is.null(annotation_res$cluster_annotation)) {
    clust_df <- annotation_res$cluster_annotation
    sym_col <- if ("expressed_symbol" %in% colnames(clust_df)) "expressed_symbol" else "SYMBOL"
    if (sym_col %in% colnames(clust_df)) {
      direct_genes <- unique(as.character(na.omit(clust_df[[sym_col]])))
      if (length(direct_genes) > 0) raw_gene_sets[["Direct"]] <- direct_genes
    }
  }
  return(raw_gene_sets)
}

#' Internal: Build Analysis Queue from Raw Gene Sets
#'
#' Constructs a named list of gene sets to analyze, including:
#' - Individual loop types (e.g., "EP", "PP")
#' - Combined set (e.g., "EP+PP")
#' - "Direct" if present
#' - "All_Combined" if multiple sets exist
#'
#' @param raw_gene_sets Named list. Output from \code{extract_target_gene_sets()}.
#' @param loop_types Character vector. Loop types originally requested (e.g., c("EP", "PP")).
#'
#' @return Named list of gene vectors for downstream analysis.
#'
#' @keywords internal
build_analysis_queue <- function(raw_gene_sets, loop_types) {
  analysis_queue <- list()
  for (lt in loop_types) {
    if (!is.null(raw_gene_sets[[lt]])) analysis_queue[[lt]] <- raw_gene_sets[[lt]]
  }
  valid_types <- intersect(loop_types, names(raw_gene_sets))
  if (length(valid_types) > 1) {
    combined_name <- paste(valid_types, collapse = "+")
    analysis_queue[[combined_name]] <- unique(unlist(raw_gene_sets[valid_types]))
  }
  if ("Direct" %in% names(raw_gene_sets)) {
    analysis_queue[["Direct"]] <- raw_gene_sets[["Direct"]]
    if (length(analysis_queue) > 1) {
      analysis_queue[["All_Combined"]] <- unique(unlist(raw_gene_sets))
    }
  }
  return(analysis_queue)
}

#' Internal: Generate LFC Boxplot Comparing Target vs Background Genes
#'
#' Creates a violin + boxplot showing log fold change distribution for target genes
#' versus all other genes, with statistical test annotation.
#'
#' @param target_genes Character vector of target gene symbols.
#' @param global_glist Named numeric vector of log fold changes (names = gene symbols).
#' @param stat_test Statistical test: "t.test" or "wilcox.test".
#' @param current_proj_name Project-specific prefix for output filename.
#' @param out_dir Output directory.
#'
#' @return Invisible \code{NULL}; saves PDF plot.
#'
#' @keywords internal
run_lfc_boxplot <- function(target_genes, global_glist, stat_test, current_proj_name, out_dir) {
  pdata <- data.frame(
    Gene = names(global_glist), LFC = as.numeric(global_glist),
    Type = factor(ifelse(names(global_glist) %in% target_genes, "Target", "Background"), levels = c("Target", "Background"))
  )
  if (sum(pdata$Type == "Target") == 0) {
    return()
  }

  test_res <- NULL
  test_name <- ""
  if (grepl("t\\.test", stat_test, ignore.case = TRUE) || grepl("^t$", stat_test, ignore.case = TRUE)) {
    test_res <- t.test(LFC ~ Type, data = pdata)
    test_name <- "T-test"
  } else {
    test_res <- wilcox.test(LFC ~ Type, data = pdata)
    test_name <- "Wilcoxon"
  }
  p_val_num <- test_res$p.value
  p_txt <- if (p_val_num < 2.2e-16) "p < 2.2e-16" else paste0("p = ", format(p_val_num, digits = 3))
  label_text <- paste0(test_name, "\n", p_txt)
  ylim_range <- quantile(pdata$LFC, c(0.01, 0.99))
  p_box <- ggplot2::ggplot(pdata, aes(x = Type, y = LFC, fill = Type)) +
    geom_violin(alpha = 0.4, color = NA, trim = TRUE) +
    geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
    coord_cartesian(ylim = ylim_range) +
    scale_fill_manual(values = c("Background" = "grey70", "Target" = "#E41A1C")) +
    theme_classic() +
    labs(title = paste0("LFC: ", current_proj_name)) +
    annotate("text", x = Inf, y = Inf, label = label_text, hjust = 1.1, vjust = 1.5, size = 3.5, fontface = "italic")
  ggsave(file.path(out_dir, paste0(current_proj_name, "_LFC_Boxplot.pdf")), p_box, width = 4, height = 5)
}

#' Run GSEA Analysis (Internal)
#'
#' @param target_genes Character vector.
#' @param global_glist Named numeric vector (Ranked LFC).
#' @param gsea_ntop Integer.
#' @param current_proj_name Character.
#' @param out_dir Character.
#' @return A GSEA result object or NULL if analysis fails.
#' @importFrom clusterProfiler GSEA
#' @importFrom enrichplot gseaplot2
#' @importFrom ggplot2 ggplot geom_line geom_hline geom_segment scale_x_continuous scale_y_continuous scale_color_gradientn theme_bw theme_classic theme_void labs theme element_blank element_text element_rect ggsave
#' @importFrom aplot plot_list
#' @keywords internal
run_gsea_analysis <- function(target_genes, global_glist, gsea_ntop, current_proj_name, out_dir) {
  # 1. 准备排序列表 (Ranked List)
  curr_glist <- global_glist

  # 处理重复值：添加微小扰动以打破 Ties
  # 注意：由于移除了 seed，这里的 runif 每次运行结果会微小不同，属于预期行为
  if (any(duplicated(curr_glist))) {
    curr_glist <- curr_glist + runif(length(curr_glist), 0, 1e-6)
    curr_glist <- sort(curr_glist, decreasing = TRUE)
  }

  # 2. 准备目标基因集 (Target Gene Set)
  curr_targets <- target_genes
  if (length(curr_targets) > gsea_ntop) {
    message("[GSEA] Downsampling targets > ", gsea_ntop)
    curr_targets <- sample(curr_targets, gsea_ntop)
  }

  names(curr_glist) <- toupper(names(curr_glist))
  curr_targets <- toupper(curr_targets)

  # 检查重叠基因数量
  if (length(intersect(curr_targets, names(curr_glist))) <= 5) {
    message("[GSEA] Skipped: Too few targets found in global gene list.")
    return(NULL)
  }

  # 3. 运行 GSEA
  # 构建自定义基因集 Data Frame
  term_df <- data.frame(term = current_proj_name, gene = curr_targets, stringsAsFactors = FALSE)

  gsea_res <- tryCatch(
    {
      # [CHANGE] Removed 'seed' argument
      clusterProfiler::GSEA(
        curr_glist,
        TERM2GENE = term_df,
        pvalueCutoff = 1.0, # 返回所有结果以便绘图
        minGSSize = 3,
        maxGSSize = 50000,
        verbose = FALSE
      )
    },
    error = function(e) {
      message("[GSEA] Calculation failed: ", e$message)
      return(NULL)
    }
  )

  # 4. 绘图
  if (!is.null(gsea_res) && nrow(gsea_res@result) > 0) {
    # === [核心修复] 获取绘图数据 ===
    # 使用 gseaplot2 生成对象并提取数据，规避使用 ::: 调用内部函数
    gsdata <- tryCatch(
      {
        # 生成 ES 曲线图 (subplots=1)
        p_temp <- enrichplot::gseaplot2(gsea_res, geneSetID = 1, subplots = 1)
        d <- p_temp$data

        # [Robustness] 确保 rank (x) 对应的数值存在
        # d$x 通常是排序后的索引 (1, 2, 3...)

        # 1. 补全 LFC 值 (用于热图颜色)
        if (!"geneList" %in% colnames(d)) {
          d$geneList <- curr_glist[d$x]
        }

        # 2. 补全 Position 值 (用于 Barcode)
        # 手动计算：当前位置的基因名是否在目标集合中
        if (!"position" %in% colnames(d)) {
          gene_at_rank <- names(curr_glist)[d$x]
          d$position <- as.numeric(gene_at_rank %in% curr_targets)
        }

        d # 返回处理后的数据框
      },
      error = function(e) {
        message("[GSEA] Plot data extraction failed: ", e$message)
        return(NULL)
      }
    )

    if (is.null(gsdata)) {
      return(NULL)
    }

    # 提取统计值用于标题
    max_rank <- max(gsdata$x)
    nes_val <- gsea_res@result$NES[1]
    p_val <- gsea_res@result$pvalue[1]

    # 设定主色调 (上调红，下调蓝)
    main_color <- if (!is.na(nes_val) && nes_val >= 0) "#E41A1C" else "#377EB8"

    # --- 子图 1: Enrichment Score (ES) ---
    p1 <- ggplot(gsdata, aes(x = x, y = runningScore)) +
      geom_line(color = main_color, linewidth = 1) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      scale_x_continuous(expand = c(0, 0), limits = c(0, max_rank)) +
      theme_bw() +
      labs(
        x = NULL, y = "ES",
        title = paste0(current_proj_name, "\nNES: ", round(nes_val, 3), "  P-val: ", format(p_val, digits = 3, scientific = TRUE))
      ) +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold")
      )

    # --- 子图 2: Hit Position (Barcode) ---
    # 筛选出 Hit 的位置
    hit_data <- gsdata[gsdata$position == 1, ]

    if (nrow(hit_data) > 0) {
      p2 <- ggplot(hit_data, aes(x = x, y = 1)) +
        geom_segment(aes(xend = x, yend = 0), color = "black", linewidth = 0.2, alpha = 0.8) +
        scale_y_continuous(expand = c(0, 0)) +
        scale_x_continuous(expand = c(0, 0), limits = c(0, max_rank)) +
        theme_void() +
        theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))
    } else {
      p2 <- ggplot() +
        theme_void() # 空图兜底
    }

    # --- 子图 3: Ranked List Metric (LFC) ---
    # 动态加载配色
    grad_colors <- if (requireNamespace("RColorBrewer", quietly = TRUE)) rev(RColorBrewer::brewer.pal(11, "PRGn")) else c("green", "white", "purple")

    p3 <- ggplot(gsdata, aes(x = x, y = geneList)) +
      geom_segment(aes(xend = x, yend = 0, color = geneList), linewidth = 0.5) +
      scale_color_gradientn(colors = grad_colors) +
      scale_x_continuous(expand = c(0, 0), limits = c(0, max_rank)) +
      theme_classic() +
      labs(x = "Rank", y = "LFC") +
      theme(legend.position = "none")

    # --- 拼图并保存 ---
    if (requireNamespace("aplot", quietly = TRUE)) {
      # 使用 aplot 对齐坐标轴
      final_plot <- aplot::plot_list(p1, p2, p3, ncol = 1, heights = c(2, 0.5, 1.5))

      out_file <- file.path(out_dir, paste0(current_proj_name, "_GSEA.pdf"))
      tryCatch(
        {
          ggsave(out_file, final_plot, width = 6, height = 6)
          message("    Saved (GSEA): ", out_file)
        },
        error = function(e) {
          message("[GSEA] Save failed: ", e$message)
        }
      )
    }
  }
}

#' Internal: Generate Expression Heatmap and Sample Correlation Plot
#'
#' Produces a z-score scaled heatmap of target genes across samples, annotated by group,
#' plus a sample-sample correlation plot.
#'
#' @param target_genes Character vector of gene symbols.
#' @param tpm_mat_raw TPM expression matrix (genes x samples).
#' @param meta_raw Metadata data frame with "SampleID" and "Group".
#' @param heatmap_ntop Max number of genes to plot (selected by highest variance).
#' @param cor_method Correlation method ("pearson", "spearman", etc.).
#' @param current_proj_name Project-specific prefix.
#' @param out_dir Output directory.
#'
#' @return Invisible \code{NULL}; saves two PDF files.
#'
#' @keywords internal
run_heatmap_and_correlation <- function(target_genes, tpm_mat_raw, meta_raw, heatmap_ntop, cor_method, current_proj_name, out_dir) {
  colnames(tpm_mat_raw) <- trimws(colnames(tpm_mat_raw))
  valid_s <- intersect(meta_raw$SampleID, colnames(tpm_mat_raw))
  if (length(valid_s) == 0) {
    return()
  }
  curr_meta <- meta_raw %>% dplyr::filter(SampleID %in% valid_s)
  valid_t <- intersect(target_genes, rownames(tpm_mat_raw))
  if (length(valid_t) <= 5) {
    return()
  }

  mat_plot <- log2(tpm_mat_raw[valid_t, as.character(curr_meta$SampleID), drop = FALSE] + 1)
  if (nrow(mat_plot) > heatmap_ntop) {
    row_vars <- apply(mat_plot, 1, var, na.rm = TRUE)
    top_genes <- head(names(sort(row_vars, decreasing = TRUE)), heatmap_ntop)
    mat_plot <- mat_plot[top_genes, , drop = FALSE]
  }
  mat_scaled <- t(scale(t(mat_plot)))
  mat_scaled[mat_scaled > 2] <- 2
  mat_scaled[mat_scaled < -2] <- -2
  mat_scaled[is.na(mat_scaled)] <- 0

  col_fun <- circlize::colorRamp2(c(-2, 0, 2), c("#313695", "white", "#A50026"))
  groups <- unique(curr_meta$Group)
  group_cols <- setNames(
    brewer.pal(max(3, length(groups)), "Set1")[seq_along(groups)],
    groups
  )
  ha <- ComplexHeatmap::HeatmapAnnotation(Group = curr_meta$Group, col = list(Group = group_cols), simple_anno_size = unit(0.3, "cm"))
  hm <- ComplexHeatmap::Heatmap(mat_scaled,
    name = "Z-score", col = col_fun, cluster_columns = FALSE,
    show_row_names = (nrow(mat_scaled) <= 80), top_annotation = ha, border = TRUE,
    column_title = paste0("Heatmap: ", current_proj_name)
  )
  pdf(file.path(out_dir, paste0(current_proj_name, "_Heatmap.pdf")), width = 7, height = 8)
  ComplexHeatmap::draw(hm)
  dev.off()

  tryCatch(
    {
      cor_mat <- cor(t(mat_plot), method = cor_method)
      method_title <- paste0(toupper(substr(cor_method, 1, 1)), substr(cor_method, 2, nchar(cor_method)))
      pdf(file.path(out_dir, paste0(current_proj_name, "_Correlation.pdf")), width = 8, height = 8)
      n_dims <- ncol(cor_mat)
      corrplot::corrplot(cor_mat,
        method = "color", order = "hclust", hclust.method = "ward.D2", addrect = 3, rect.col = "black", rect.lwd = 2,
        type = "full", tl.pos = "n", col = corrplot::COL2("PRGn", 200), cl.pos = "r", cl.ratio = 0.15, cl.align.text = "l", cl.offset = 0.5,
        outline = FALSE, title = paste0("Correlation Analysis: ", current_proj_name), mar = c(1, 1, 3, 4)
      )
      mtext(paste0(method_title, "\nIndex"), side = 4, line = 0.2, at = n_dims + 3, las = 1, adj = 0.5, cex = 0.9, font = 2)
      dev.off()
    },
    error = function(e) message("   ⚠️ Correlation plot failed: ", e$message)
  )
}

#' Internal: Perform GO Enrichment and Generate Circular Network Plot
#'
#' Runs GO enrichment (Biological Process) and creates a circular gene-pathway network
#' highlighting top pathways and most differentially expressed genes.
#'
#' @param target_genes Character vector of gene symbols.
#' @param org_db Organism database (e.g., "org.Hs.eg.db").
#' @param global_glist Named log fold change vector for coloring genes (names = symbols).
#' @param cnet_ntop Number of top pathways to display in the plot.
#' @param current_proj_name Project-specific prefix.
#' @param out_dir Output directory.
#'
#' @return Data frame of GO results (or \code{NULL} if failed).
#'
#' @importFrom clusterProfiler bitr enrichGO
#' @importFrom enrichplot cnetplot
#' @importFrom ggplot2 ggsave theme element_text
#' @keywords internal
run_go_enrichment <- function(target_genes, org_db, global_glist, cnet_ntop, current_proj_name, out_dir) {
  eg <- tryCatch(
    {
      clusterProfiler::bitr(target_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org_db)
    },
    error = function(e) NULL
  )
  if (is.null(eg) || nrow(eg) == 0) {
    return(NULL)
  }

  go_res <- tryCatch(
    {
      clusterProfiler::enrichGO(
        gene = eg$ENTREZID, OrgDb = org_db, ont = "BP", pAdjustMethod = "BH",
        pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE
      )
    },
    error = function(e) NULL
  )
  if (is.null(go_res) || nrow(go_res) == 0) {
    return(NULL)
  }

  # Circular network plot
  tryCatch(
    {
      top_n_cats <- 5
      go_sub <- go_res
      go_sub@result <- head(go_res@result, top_n_cats)
      all_genes_in_top <- unique(unlist(strsplit(go_sub@result$geneID, "/")))
      valid_genes <- intersect(all_genes_in_top, names(global_glist))
      if (length(valid_genes) > 0) {
        lfc_in_top <- global_glist[valid_genes]
        lfc_in_top <- lfc_in_top[order(abs(lfc_in_top), decreasing = TRUE)]
        keep_genes <- names(head(lfc_in_top, cnet_ntop))
        for (i in seq_len(nrow(go_sub@result))) {
          g_in_path <- unlist(strsplit(go_sub@result$geneID[i], "/"))
          g_keep <- intersect(g_in_path, keep_genes)
          go_sub@result$geneID[i] <- paste(g_keep, collapse = "/")
        }
        p_cnet <- enrichplot::cnetplot(go_sub,
          foldChange = global_glist, circular = TRUE, colorEdge = TRUE,
          showCategory = top_n_cats, cex_label_category = 1.5, cex_label_gene = 0.8
        ) +
          scale_color_gradient2(name = "Log2FC", low = "#313695", mid = "white", high = "#A50026", midpoint = 0) +
          labs(
            title = paste0("GO Circular Network: ", current_proj_name),
            subtitle = paste0("Top ", top_n_cats, " Pathways | Top ", length(keep_genes), " Genes (by |LFC|)")
          ) +
          theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "right")
        ggsave(file.path(out_dir, paste0(current_proj_name, "_GO_Circular_Network_TopGenes.pdf")), p_cnet, width = 10, height = 9)
      }
    },
    error = function(e) message("GO Circular Plot failed: ", e$message)
  )

  return(as.data.frame(go_res))
}

#' Internal: Construct and Visualize PPI Network
#'
#' Maps target genes to STRING, builds a subnetwork, and visualizes it with node size
#' reflecting degree and color reflecting log fold change.
#'
#' @param target_genes Character vector of gene symbols.
#' @param global_glist Named log fold change vector.
#' @param org_db Organism database to infer species (e.g., "org.Hs.eg.db" -> 9606).
#' @param ppi_score Minimum STRING interaction score (0–1000).
#' @param ppi_ntop Max number of genes to include (prioritized by |LFC|).
#' @param current_proj_name Project-specific prefix.
#' @param out_dir Output directory.
#'
#' @return Invisible \code{NULL}; saves PDF plot if successful.
#'
#' @keywords internal
run_ppi_analysis <- function(target_genes, global_glist, org_db, ppi_score, ppi_ntop, current_proj_name, out_dir) {
  message("   >>> [PPI Init] Threshold = ", ppi_score)
  species_id <- if (grepl("Mm", org_db, ignore.case = TRUE)) 10090 else 9606
  string_db_obj <- STRINGdb::STRINGdb$new(version = "11.5", species = species_id, score_threshold = ppi_score, input_directory = "")
  ppi_genes <- target_genes
  if (length(ppi_genes) > ppi_ntop) {
    valid_in_lfc <- intersect(ppi_genes, names(global_glist))
    if (length(valid_in_lfc) > 0) {
      ppi_genes <- head(valid_in_lfc[order(abs(global_glist[valid_in_lfc]), decreasing = TRUE)], ppi_ntop)
    } else {
      ppi_genes <- head(ppi_genes, ppi_ntop)
    }
  }
  targets_mapped <- string_db_obj$map(data.frame(gene = ppi_genes), "gene", removeUnmappedRows = TRUE)
  if (nrow(targets_mapped) == 0) {
    return()
  }

  hits <- targets_mapped$STRING_id
  if (length(hits) <= 1) {
    return()
  }

  g_string <- string_db_obj$get_subnetwork(hits)
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
    ggraph::geom_edge_link(aes(alpha = combined_score), color = "grey60", width = 0.5, show.legend = FALSE) +
    ggraph::geom_node_point(aes(color = lfc, size = deg), stroke = 0.5) +
    ggraph::geom_node_text(aes(label = label_text),
      repel = TRUE, size = 3.5, max.overlaps = Inf,
      fontface = "bold", bg.color = "white", bg.r = 0.1
    ) +
    scale_color_gradient2(low = "#313695", mid = "white", high = "#A50026", midpoint = 0, name = "LFC") +
    scale_size_continuous(range = c(2, 8), guide = "none") +
    ggraph::scale_edge_alpha_continuous(range = c(0.4, 0.9)) +
    ggraph::theme_graph(base_family = "sans", background = "white") +
    labs(
      title = paste0("PPI Network: ", current_proj_name),
      subtitle = paste0("Nodes: ", num_nodes, " | Score: ", min(igraph::E(g_string)$combined_score, na.rm = TRUE), "-", max(igraph::E(g_string)$combined_score, na.rm = TRUE))
    ) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  ggsave(file.path(out_dir, paste0(current_proj_name, "_PPI_Network_ggraph.pdf")), p_ppi, width = 9, height = 9)
}

#' Internal: Generate Summary Lollipop Plot of GO Results Across Loop Types
#'
#' Combines GO enrichment results from multiple loop types into a summary lollipop plot.
#' Visualizes the top enriched terms for each loop type side-by-side.
#'
#' @param all_go_results Named list of GO result data frames (from \code{perform_go_enrichment}).
#' @param base_project_name Base project name for output file.
#' @param out_dir Output directory.
#'
#' @return Invisible \code{NULL}; saves PDF plot.
#'
#' @importFrom dplyr bind_rows slice_min mutate arrange group_by
#' @importFrom ggplot2 ggplot aes geom_segment geom_point facet_grid theme_bw labs scale_color_brewer ggsave
#' @keywords internal
plot_summary_go_lollipop <- function(all_go_results, base_project_name, out_dir) {
  message("\n>>> Plotting Summary GO Lollipop...")
  final_go_df <- do.call(rbind, all_go_results)
  final_go_df$logP <- -log10(final_go_df$pvalue)

  max_logp <- max(final_go_df$logP, na.rm = TRUE)
  max_count <- max(final_go_df$Count, na.rm = TRUE)
  scale_f <- max_logp / max_count * 1.0

  final_go_df <- final_go_df %>%
    dplyr::group_by(LoopType) %>%
    dplyr::arrange(logP) %>%
    dplyr::mutate(Description_unique = paste(Description, LoopType, sep = "___")) %>%
    dplyr::ungroup()
  final_go_df$Description_unique <- factor(final_go_df$Description_unique, levels = final_go_df$Description_unique)

  p_go_dual <- ggplot(final_go_df, aes(y = Description_unique)) +
    geom_segment(aes(x = 0, xend = logP, yend = Description_unique, color = LoopType), size = 5) +
    geom_point(aes(x = logP, color = LoopType), size = 6) +
    geom_path(aes(x = Count * scale_f, group = LoopType), color = "grey40", size = 1) +
    geom_point(aes(x = Count * scale_f), color = "grey40", size = 3) +
    scale_x_continuous(
      name = expression(-log[10](p - value)), expand = expansion(mult = c(0, 0.1)),
      sec.axis = sec_axis(~ . / scale_f, name = "Gene Counts")
    ) +
    scale_y_discrete(labels = function(x) sub("___.*", "", x)) +
    scale_color_brewer(palette = "Set2") +
    labs(y = NULL, title = "GO Enrichment Analysis (BP)", subtitle = "Lollipop: Significance (-logP) | Line: Gene Counts") +
    theme_classic() +
    theme(
      axis.text.y = element_text(size = 10, color = "black"),
      axis.text.x.bottom = element_text(color = "black"),
      axis.text.x.top = element_text(color = "grey40", face = "bold"),
      axis.title.x.top = element_text(color = "grey40", face = "bold"),
      axis.line = element_line(size = 0.8),
      panel.grid.major.y = element_line(color = "grey95", linetype = "dashed"),
      legend.position = "none"
    ) +
    facet_wrap(~LoopType, scales = "free_y", ncol = 1)

  plot_height <- max(6, length(unique(final_go_df$Description_unique)) * 0.35)
  ggsave(file.path(out_dir, paste0(base_project_name, "_GO_DualAxis_Lollipop_Thick.pdf")), p_go_dual, width = 9, height = plot_height)
}
