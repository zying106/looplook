# Classify genomic feature based on ChIPseeker annotation string.
# Returns "P" (Promoter), "E" (Distal Intergenic/Enhancer), or "G" (Gene body).
get_feature_class <- function(anno_str) {
  if (grepl("Promoter", anno_str, ignore.case = TRUE)) {
    return("P")
  }
  if (grepl("Distal Intergenic", anno_str, ignore.case = TRUE)) {
    return("E")
  }
  return("G")
}

# Standardize ChIPseeker output: extract broad annotation and keep detailed version.
format_annotation_columns <- function(df) {
  if ("annotation" %in% colnames(df)) {
    df <- df %>%
      dplyr::rename(detail_anno = annotation) %>%
      dplyr::mutate(annotation = gsub(" \\(.*", "", detail_anno)) %>%
      dplyr::relocate(annotation, .before = detail_anno)
  }
  return(df)
}


#' @title Annotate Chromatin Loops, Identify Clusters, and Integrate Auxiliary BED Features
#'
#' @description
#' Performs comprehensive genomic annotation of chromatin loop anchors (e.g., from HiChIP, ChIA-PET, or Hi-C)
#' and classifies loops into biologically meaningful types (e.g., Enhancer-Promoter (EP), Promoter-Promoter (PP)).
#'
#' \strong{Key Feature - Auxiliary Annotation:}
#' This function explicitly supports **integrative analysis** by accepting an auxiliary BED file
#' (via \code{target_bed}). This allows users to annotate external genomic features (e.g., ATAC-seq peaks,
#' ChIP-seq binding sites, or GWAS variants) and map them directly to chromatin loop clusters to identify
#' feature-specific connectivity and target genes.
#'
#' @details
#' This function implements a streamlined "Fast & Clean" workflow:
#' \enumerate{
#'   \item \strong{Loop Annotation}: Annotates loop anchors using `ChIPseeker` to determine genomic context.
#'   \item \strong{Classification}: Classifies loops (EP, PP, etc.) based on the genomic features of both anchors.
#'   \item \strong{Clustering}: Aggregates spatially connected anchors into clusters using graph theory to identify regulatory hubs.
#'   \item \strong{Target Integration}: If \code{target_bed} is provided, the function annotates these external peaks
#'   and calculates their overlap with loop clusters. This links the auxiliary features (e.g., Open Chromatin)
#'   to the 3D chromatin architecture, outputting a specific \code{target_annotation} dataset.
#' }
#'
#'
#' @param bedpe_file Character. Path to the BEDPE file containing loop coordinates.
#'   Must contain at least 6 columns: chr1, start1, end1, chr2, start2, end2.
#' @param target_bed Optional character. Path to an **auxiliary BED file** (e.g., ATAC-seq peaks, enhancer lists)
#'   to be annotated and integrated with the loop data. If provided, the output will include an
#'   \code{overlap_annotation} component linking these features to loop clusters.
#' @param species Character. Genome build for annotation. Supported: "hg38", "hg19", "mm10", "mm9".
#' @param tss_region Numeric vector of length 2. Region around TSS to be considered as promoter.
#'   Default: \code{c(-2000, 2000)}.
#' @param out_dir Character. Output directory for intermediate files (default: current directory).
#' @param project_name Character. Prefix for internal labeling of results (default: "HiChIP").
#' @param color_palette Character. Name of color palette (reserved for future plotting consistency).
#' @param karyo_bin_size Numeric. Bin size for karyotype visualization (reserved for downstream pipelines).
#'
#' @return A list with three components:
#' \describe{
#'   \item{\code{loop_annotation}}{Data frame with per-loop details (type, assigned genes, anchor info).}
#'   \item{\code{cluster_annotation}}{Aggregated cluster-level annotations and gene summaries.}
#'   \item{\code{target_annotation}}{Annotation of the \code{target_bed} features (if provided), including their connectivity to loop clusters and target genes.}
#' }
#'
#' @importFrom ChIPseeker annotatePeak
#' @importFrom GenomicRanges makeGRangesFromDataFrame findOverlaps reduce
#' @importFrom igraph graph_from_data_frame components
#' @importFrom dplyr %>%
#'
#' @export
#'
#' @examples
#' # 1. Get paths to example files included in the package
#' bedpe_path <- system.file("extdata", "example_loops_1.bedpe", package = "looplook")
#' bed_path <- system.file("extdata", "example_peaks.bed", package = "looplook")
#'
#' # 2. Check if files and required annotation database exist
#' # (ChIPseeker requires the TxDb package for the specific genome)
#' if (bedpe_path != "" &&
#'   requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE)) {
#'   # 3. Run the function with an auxiliary BED file (e.g., ATAC-seq)
#'   res <- annotate_peaks_and_loops(
#'     bedpe_file = bedpe_path,
#'     target_bed = bed_path, # <--- Integrating auxiliary data here
#'     species = "hg38",
#'     tss_region = c(-2000, 2000),
#'     out_dir = tempdir(),
#'     project_name = "Example_HiChIP"
#'   )
#'
#'   # 4. Access results
#'   head(res$loop_annotation)
#'
#'   # Access the annotation of the auxiliary BED file
#'   if (!is.null(res$target_annotation)) {
#'     head(res$target_annotation)
#'   }
#' }
annotate_peaks_and_loops <- function(
  bedpe_file,
  target_bed = NULL,
  species = "hg38",
  tss_region = c(-2000, 2000),
  out_dir = "./results",
  project_name = "HiChIP",
  color_palette = "Set2",
  karyo_bin_size = 1e5
) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  # Load DBs
  if (species == "hg38") {
    txdb_pkg <- "TxDb.Hsapiens.UCSC.hg38.knownGene"
    org_db_pkg <- "org.Hs.eg.db"
  } else if (species == "hg19") {
    txdb_pkg <- "TxDb.Hsapiens.UCSC.hg19.knownGene"
    org_db_pkg <- "org.Hs.eg.db"
  } else if (species == "mm10") {
    txdb_pkg <- "TxDb.Mmusculus.UCSC.mm10.knownGene"
    org_db_pkg <- "org.Mm.eg.db"
  } else if (species == "mm9") {
    txdb_pkg <- "TxDb.Mmusculus.UCSC.mm9.knownGene"
    org_db_pkg <- "org.Mm.eg.db"
  } else {
    stop("Unsupported species.")
  }

  txdb_obj <- get(txdb_pkg)

  # Step 1: Read loops
  message("Step 1: Reading BEDPE file...")
  loops <- read_robust_general(bedpe_file, min_cols = 6, desc = "BEDPE")
  colnames(loops)[seq_len(6)] <- c("chr1", "start1", "end1", "chr2", "start2", "end2")
  n_input_loops <- nrow(loops)
  message("    >>> Raw Input Loops: ", prettyNum(n_input_loops, big.mark = ","))

  # Construct Anchors & Clustering
  anchors <- dplyr::bind_rows(loops %>% dplyr::select(chr = chr1, start = start1, end = end1), loops %>% dplyr::select(chr = chr2, start = start2, end = end2)) %>%
    dplyr::distinct() %>%
    dplyr::mutate(anchor_id = paste0("A", seq_len(dplyr::n())))
  loops <- loops %>%
    dplyr::left_join(anchors %>% dplyr::select(chr, start, end, a1_id = anchor_id), by = c("chr1" = "chr", "start1" = "start", "end1" = "end")) %>%
    dplyr::left_join(anchors %>% dplyr::select(chr, start, end, a2_id = anchor_id), by = c("chr2" = "chr", "start2" = "start", "end2" = "end"))

  message("Step 2: Clustering loops...")
  valid_loops <- loops %>% dplyr::filter(!is.na(a1_id) & !is.na(a2_id))
  g <- igraph::graph_from_data_frame(valid_loops[, c("a1_id", "a2_id")], directed = FALSE)
  comp <- igraph::components(g)
  anchors$cluster_id <- NA
  comm <- intersect(anchors$anchor_id, names(comp$membership))
  if (length(comm) > 0) anchors$cluster_id[match(comm, anchors$anchor_id)] <- as.character(comp$membership[comm])
  anchors <- anchors %>% dplyr::filter(!is.na(cluster_id))
  loops <- loops %>% dplyr::left_join(anchors %>% dplyr::select(anchor_id, cluster_id), by = c("a1_id" = "anchor_id"))

  # Force unique names for Cluster Regions
  gr_anchors <- GenomicRanges::makeGRangesFromDataFrame(anchors, keep.extra.columns = TRUE)
  gr_list <- GenomicRanges::GRangesList(split(gr_anchors, gr_anchors$cluster_id))
  cluster_regions <- unlist(GenomicRanges::reduce(gr_list))
  cluster_regions$cluster_id <- names(cluster_regions)
  names(cluster_regions) <- paste0("peak_", seq_along(cluster_regions))

  # Step 3: Classification
  message("Step 3: Biological Classification...")
  anchor_anno <- ChIPseeker::annotatePeak(gr_anchors, TxDb = txdb_obj, tssRegion = tss_region, annoDb = org_db_pkg, verbose = FALSE)
  anchor_anno_df <- format_annotation_columns(as.data.frame(anchor_anno))
  anchor_anno_df$type_code <- vapply(
    anchor_anno_df$annotation,
    get_feature_class,
    FUN.VALUE = character(1)
  )
  map_info <- anchor_anno_df %>% dplyr::select(anchor_id, type_code, SYMBOL)
  loops_annotated <- loops %>%
    dplyr::left_join(map_info %>% dplyr::rename(t1 = type_code, s1 = SYMBOL), by = c("a1_id" = "anchor_id")) %>%
    dplyr::left_join(map_info %>% dplyr::rename(t2 = type_code, s2 = SYMBOL), by = c("a2_id" = "anchor_id"))

  classify_and_target <- function(t1, t2, s1, s2) {
    if (is.na(t1) || is.na(t2)) {
      return(list(type = "Unknown", gene = NA))
    }
    types <- sort(c(t1, t2))
    code <- paste(types, collapse = "-")
    genes <- c()
    if (code == "P-P") genes <- c(s1, s2) else if (code == "E-P" || code == "G-P") genes <- if (t1 == "P") s1 else s2 else if (code == "E-G") genes <- if (t1 == "G") s1 else s2 else genes <- c(s1, s2)
    return(list(type = code, gene = paste(unique(na.omit(genes)), collapse = ";")))
  }
  res_list <- apply(loops_annotated, 1, function(r) classify_and_target(r["t1"], r["t2"], r["s1"], r["s2"]))
  loops_annotated$loop_type <- vapply(
    res_list,
    function(x) x$type,
    FUN.VALUE = character(1)
  )
  loops_annotated$loop_gene <- vapply(
    res_list,
    function(x) x$gene,
    FUN.VALUE = character(1)
  )

  loop_annotation_final <- loops_annotated %>% dplyr::select(chr1, start1, end1, chr2, start2, end2, cluster_id, loop_type, loop_genes = loop_gene, anchor1_gene = s1, anchor1_type = t1, anchor2_gene = s2, anchor2_type = t2)

  # Step 4: Aggregation
  message("Step 4: Aggregating Cluster Annotations...")
  extract_genes <- function(genes_vec) {
    paste(unique(na.omit(unlist(strsplit(genes_vec, ";")))), collapse = ";")
  }
  agg_cluster <- loops_annotated %>%
    dplyr::filter(!is.na(cluster_id)) %>%
    dplyr::group_by(cluster_id) %>%
    dplyr::summarise(loop_count = dplyr::n(), all_loop_genes = extract_genes(loop_gene)) %>%
    dplyr::ungroup()

  gene_annot <- ChIPseeker::annotatePeak(cluster_regions, TxDb = txdb_obj, tssRegion = tss_region, annoDb = org_db_pkg, verbose = FALSE)
  cluster_info <- format_annotation_columns(as.data.frame(gene_annot))

  # [Renaming 1] GENENAME -> Gene_description
  if ("GENENAME" %in% colnames(cluster_info)) cluster_info <- cluster_info %>% dplyr::rename(Gene_description = GENENAME)

  cluster_info$cluster_id <- as.character(cluster_info$cluster_id)
  cluster_info <- dplyr::left_join(cluster_info, agg_cluster, by = "cluster_id")

  # =========================================================
  # Step 5: Target Analysis (Lightweight - NO LINKS TABLE)
  # =========================================================
  bed_info <- NULL
  target_connected_loops <- NULL

  if (!is.null(target_bed)) {
    message("Step 5: Integrating Target Annotations...")

    bed_target <- read_robust_general(target_bed, min_cols = 3, desc = "Target BED")
    colnames(bed_target)[c(1, 2, 3)] <- c("chr", "start", "end")
    gr_bed <- GenomicRanges::makeGRangesFromDataFrame(bed_target)
    gr_bed$input_id <- paste0("Peak_", seq_along(gr_bed))
    names(gr_bed) <- gr_bed$input_id

    # Annotate Peaks
    bed_annot <- ChIPseeker::annotatePeak(gr_bed, TxDb = txdb_obj, tssRegion = tss_region, annoDb = org_db_pkg, verbose = FALSE)
    bed_info <- format_annotation_columns(as.data.frame(bed_annot))

    # [Renaming 2] GENENAME -> Gene_description
    if ("GENENAME" %in% colnames(bed_info)) bed_info <- bed_info %>% dplyr::rename(Gene_description = GENENAME)

    # Find Overlaps
    hits <- GenomicRanges::findOverlaps(gr_bed, cluster_regions)

    if (length(hits) > 0) {
      hit_cluster_ids <- unique(cluster_regions$cluster_id[S4Vectors::subjectHits(hits)])
      message("    >>> Found overlaps with ", length(hit_cluster_ids), " unique loop clusters.")

      # 仅提取 Loop 用于画图，不生成详细连接表，不写 CSV
      target_connected_loops <- loop_annotation_final %>% dplyr::filter(cluster_id %in% hit_cluster_ids)
    }

    # Update Bed Info (Summary only)
    hit_df <- data.frame(qid = S4Vectors::queryHits(hits), sid = S4Vectors::subjectHits(hits))
    cluster_map <- data.frame(sid = seq_along(cluster_regions), cluster_id = as.character(cluster_regions$cluster_id)) %>% dplyr::left_join(agg_cluster, by = "cluster_id")
    hit_df <- hit_df %>% dplyr::left_join(cluster_map, by = "sid")

    # [Renaming 3] loop_genes_Total -> all_annoted_gene
    summary_df <- hit_df %>%
      dplyr::group_by(qid) %>%
      dplyr::summarise(all_annoted_gene = extract_genes(all_loop_genes), .groups = "drop") %>%
      dplyr::mutate(join_id = paste0("Peak_", qid))

    bed_info <- dplyr::left_join(bed_info, summary_df, by = c("input_id" = "join_id")) %>%
      dplyr::select(-any_of(c("join_id", "qid")))
  }

  # Step 6: Visualization
  message("Step 6: Generating Visualizations...")
  plot_df <- loop_annotation_final
  loop_types_sorted <- sort(unique(plot_df$loop_type))
  if (!exists("get_colors")) stop("Helper 'get_colors' missing.")
  custom_colors <- get_colors(length(loop_types_sorted), color_palette)
  names(custom_colors) <- loop_types_sorted

  if (exists("draw_rose_plot")) try(draw_rose_plot(plot_df, project_name, file.path(out_dir, paste0(project_name, "_Basic_Rose.pdf")), custom_colors), silent = TRUE)
  if (exists("draw_circular_bar_plot")) try(draw_circular_bar_plot(plot_df, project_name, file.path(out_dir, paste0(project_name, "_Basic_Circular.pdf")), custom_colors), silent = TRUE)

  # Karyotype
  unique_loop_genes <- plot_df %>%
    dplyr::filter(!is.na(loop_genes) & loop_genes != "") %>%
    tidyr::separate_rows(loop_genes, sep = ";") %>%
    dplyr::pull(loop_genes) %>%
    unique() %>%
    trimws()
  unique_loop_genes <- unique_loop_genes[unique_loop_genes != ""]
  if (length(unique_loop_genes) > 0 && exists("draw_karyo_heatmap_internal")) {
    all_genes_gr <- GenomicFeatures::genes(txdb_obj)
    map <- AnnotationDbi::select(if (grepl("hg", species)) get("org.Hs.eg.db") else get("org.Mm.eg.db"), keys = as.character(S4Vectors::mcols(all_genes_gr)$gene_id), columns = "SYMBOL", keytype = "ENTREZID")
    S4Vectors::mcols(all_genes_gr)$SYMBOL <- map$SYMBOL[match(S4Vectors::mcols(all_genes_gr)$gene_id, map$ENTREZID)]
    target_genes_gr <- all_genes_gr[S4Vectors::mcols(all_genes_gr)$SYMBOL %in% unique_loop_genes]
    try(draw_karyo_heatmap_internal(target_genes_gr, "Gene Load (Basic)", file.path(out_dir, paste0(project_name, "_Basic_Karyo_Genes.pdf")), karyo_bin_size, 0.99, txdb_obj, species, "Genes"), silent = TRUE)
  }

  # Anchor Load
  a1_coords <- plot_df %>% dplyr::select(chr = chr1, start = start1, end = end1)
  a2_coords <- plot_df %>% dplyr::select(chr = chr2, start = start2, end = end2)
  all_anchors <- dplyr::bind_rows(a1_coords, a2_coords) %>% dplyr::distinct()
  if (nrow(all_anchors) > 0 && exists("draw_karyo_heatmap_internal")) {
    gr_anchors_plot <- GenomicRanges::makeGRangesFromDataFrame(all_anchors)
    try(draw_karyo_heatmap_internal(gr_anchors_plot, "Loop Anchor Load (Basic)", file.path(out_dir, paste0(project_name, "_Basic_Karyo_Anchors.pdf")), karyo_bin_size, 0.99, txdb_obj, species, "Anchors"), silent = TRUE)
  }

  # Flower
  temp_df_flower <- plot_df %>%
    dplyr::filter(!is.na(loop_genes) & loop_genes != "") %>%
    tidyr::separate_rows(loop_genes, sep = ";") %>%
    dplyr::mutate(loop_genes = trimws(loop_genes)) %>%
    dplyr::filter(loop_genes != "")
  gene_sets <- split(temp_df_flower$loop_genes, temp_df_flower$loop_type)
  gene_sets <- lapply(gene_sets, unique)
  if (length(gene_sets) > 1 && exists("draw_flower_simplified")) try(draw_flower_simplified(gene_sets, project_name, file.path(out_dir, paste0(project_name, "_Basic_Flower.pdf")), custom_colors), silent = TRUE)

  # Target Visualization
  if (!is.null(bed_info)) {
    message("    Plotting Target Visualizations...")
    # 1. Pie Chart
    if (exists("draw_target_annotation_pie")) try(draw_target_annotation_pie(bed_info, project_name, file.path(out_dir, paste0(project_name, "_Target_Genomic_Distribution.pdf")), color_palette), silent = TRUE)

    # 2. Connectivity Bar
    if (exists("draw_target_connectivity_bar")) try(draw_target_connectivity_bar(bed_info, cluster_info, project_name, file.path(out_dir, paste0(project_name, "_Target_Connectivity.pdf")), color_palette), silent = TRUE)

    # 3. Donut Chart (Uses the lightweight subset)
    if (!is.null(target_connected_loops) && nrow(target_connected_loops) > 0 && exists("draw_target_loop_donut")) {
      try(draw_target_loop_donut(target_connected_loops, project_name, file.path(out_dir, paste0(project_name, "_Target_Loop_Donut.pdf")), custom_colors), silent = TRUE)
    }
  }

  # Step 7: Export
  message("Step 7: Exporting to Excel...")
  tryCatch(
    {
      wb <- openxlsx::createWorkbook()
      openxlsx::addWorksheet(wb, "Loop Annotation")
      openxlsx::writeData(wb, "Loop Annotation", loop_annotation_final)
      openxlsx::addWorksheet(wb, "Cluster Annotation")
      openxlsx::writeData(wb, "Cluster Annotation", cluster_info)
      if (!is.null(bed_info)) {
        openxlsx::addWorksheet(wb, "Target Annotation")
        openxlsx::writeData(wb, "Target Annotation", bed_info)
      }
      # No Link Sheet Generated
      openxlsx::saveWorkbook(wb, file.path(out_dir, paste0(project_name, "_Basic_Results.xlsx")), overwrite = TRUE)
      message("    Excel file saved.")
    },
    error = function(e) message("Export failed: ", e$message)
  )

  message("✅ Analysis Complete. Results saved to: ", out_dir)
  return(list(loop_annotation = loop_annotation_final, cluster_annotation = cluster_info, target_annotation = bed_info))
}

#' Refine Loop Classification Using Gene Expression Data
#'
#' Reclassifies loop anchors based on RNA-seq expression (e.g., TPM) and integrates auxiliary target features.
#'
#' @description
#' This function filters false-positive regulatory interactions by leveraging gene expression data.
#' It downgrades "inactive" Promoters to Enhancers, re-evaluates loop types (e.g., `P-P` might become `E-P`),
#' and optionally maps these refined interactions to an auxiliary BED file (e.g., ATAC-seq).
#'
#' @details
#' \strong{Refinement Logic:}
#' \itemize{
#'   \item **Active Gene Detection:** Genes with expression > \code{threshold} are considered active.
#'   \item **Anchor Downgrading:** If an anchor is annotated as a Promoter ("P") but the gene is inactive, it is reclassified as an Enhancer ("E") if \code{downgrade_nonexpressed = TRUE}.
#'   \item **Loop Reclassification:** Loop types are updated based on the new anchor types (e.g., an inactive `P-P` loop becomes `E-E`).
#' }
#'
#' \strong{Target Integration:}
#' If \code{target_bed} is provided, the function maps the \strong{refined, expressed genes} to these external peaks
#' via physical overlap with loop clusters. This allows you to see which peaks are connected to \emph{active} genes.
#'
#' @param annotation_res List. Output from \code{\link{annotate_peaks_and_loops}} (The "Before Refinement" data).
#' @param expr_matrix_file Character. Path to expression matrix (rows = genes, columns = samples). Supports CSV or TSV.
#' @param sample_columns Vector (Integer or Character). Columns in expression matrix to average (e.g., c("Ctrl1", "Ctrl2")).
#' @param target_bed Optional character. Path to an auxiliary BED file (e.g., ATAC-seq). If provided, connectivity analysis is updated based on expression.
#' @param threshold Numeric. Minimum expression (e.g., TPM > 1) to consider a gene "active" (default: 1).
#' @param unit_type Character. Expression unit for labeling (e.g., "TPM"; default: "TPM").
#' @param draw_original_violin Logical. If TRUE, draws QC violin plots.
#' @param species Character. Genome build ("hg38", "hg19", "mm10", "mm9").
#' @param out_dir Character. Directory to save all PDF plots and Excel file.
#' @param project_name Character. Prefix for output files (default: "HiChIP_Filtered").
#' @param color_palette Character. Color palette for plots (default: "Set2").
#' @param karyo_bin_size Numeric. Bin size for karyotype heatmap (default: 1e5).
#' @param color_saturation Numeric. Saturation level for karyotype colors (default: 0.99).
#' @param downgrade_nonexpressed Logical. If TRUE, non-expressed promoters become enhancers (default: TRUE).
#'
#' @return A list with refined annotations:
#'   \itemize{
#'     \item \code{loop_annotation}: Updated loop types and functional genes.
#'     \item \code{cluster_annotation}: Cluster info with updated expressed genes.
#'     \item \code{target_annotation}: (If \code{target_bed} provided) Annotation of external peaks linked to active genes.
#'   }
#'
#' @importFrom dplyr %>% filter group_by summarise ungroup mutate select rename left_join inner_join arrange desc
#' @importFrom ggplot2 ggplot aes geom_col geom_text scale_fill_brewer theme_classic labs scale_y_continuous expansion ggsave geom_rect coord_polar xlim scale_fill_manual theme_void annotate
#' @importFrom scales percent
#' @importFrom utils read.table write.table head
#' @importFrom GenomicRanges makeGRangesFromDataFrame findOverlaps
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom openxlsx createWorkbook addWorksheet writeData saveWorkbook
#' @export
#' @examples
#' # 1. Prepare prerequisites
#' bedpe_path <- system.file("extdata", "example_loops_1.bedpe", package = "looplook")
#' bed_path <- system.file("extdata", "example_peaks.bed", package = "looplook")
#'
#' # Check environment requirements
#' if (bedpe_path != "" &&
#'   requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE)) {
#'   # Step A: Run basic annotation first
#'   # Note: We suppress messages to keep the example clean
#'   res_raw <- suppressMessages(annotate_peaks_and_loops(
#'     bedpe_file = bedpe_path,
#'     species = "hg38",
#'     out_dir = tempdir()
#'   ))
#'
#'   # 2. Prepare a dummy Expression Matrix
#'   expr_df <- data.frame(
#'     GeneID = c("TP53", "MYC", "GAPDH", "ACTB"),
#'     Sample1 = c(10.5, 0.1, 50.2, 0.0),
#'     Sample2 = c(11.0, 0.2, 51.0, 0.1)
#'   )
#'   expr_file <- tempfile(fileext = ".txt")
#'   write.table(expr_df, expr_file, sep = "\t", quote = FALSE, row.names = FALSE)
#'
#'   # 3. Run the Refinement Function with Target Integration
#'   res_refined <- refine_loop_anchors_by_expression(
#'     annotation_res = res_raw,
#'     expr_matrix_file = expr_file,
#'     sample_columns = c("Sample1", "Sample2"),
#'     target_bed = bed_path, # <--- Integrate ATAC-seq here
#'     threshold = 1.0,
#'     species = "hg38",
#'     out_dir = tempdir(),
#'     project_name = "Example_Refined"
#'   )
#'
#'   # 4. Check results
#'   head(res_refined$loop_annotation)
#'   if (!is.null(res_refined$target_annotation)) {
#'     head(res_refined$target_annotation)
#'   }
#' }
refine_loop_anchors_by_expression <- function(
    annotation_res,
    expr_matrix_file,
    sample_columns,
    target_bed = NULL,
    threshold = 1,
    unit_type = "TPM",
    draw_original_violin = TRUE,
    species = "hg38",
    out_dir = "./results/filtered",
    project_name = "HiChIP_Filtered",
    color_palette = "Set2",
    karyo_bin_size = 1e5,
    color_saturation = 0.99,
    downgrade_nonexpressed = TRUE
) {
  # --- 1. 加载数据 ---
  message(">>> [Step 1] Loading Data...")

  # 读取 Loop 数据
  if (is.null(annotation_res$loop_annotation)) stop("'loop_annotation' missing.")
  loop_df <- annotation_res$loop_annotation

  # 读取表达矩阵
  if (!file.exists(expr_matrix_file)) stop("Expression file not found.")
  d <- tryCatch(read.table(expr_matrix_file, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE), error = function(e) NULL)
  if (is.null(d)) d <- tryCatch(read.table(expr_matrix_file, header = TRUE, sep = ",", row.names = 1, check.names = FALSE), error = function(e) NULL)
  expr_mat <- d

  # 处理样本列
  if (is.character(sample_columns)) {
    sub_mat <- expr_mat[, intersect(sample_columns, colnames(expr_mat)), drop = FALSE]
  } else {
    sub_mat <- expr_mat[, sample_columns, drop = FALSE]
  }

  # 计算活跃基因白名单
  vals <- if (ncol(sub_mat) > 1) rowMeans(sub_mat, na.rm = TRUE) else sub_mat[, 1]
  whitelist <- names(vals)[vals > threshold & !is.na(vals) & names(vals) != ""]
  message("    >>> Active Genes: ", length(whitelist))

  # --- 2. 核心重分类 (Reclassification) ---
  message(">>> [Step 2] Reclassifying Loops...")

  check_anc <- function(g, t, allow = whitelist, down = downgrade_nonexpressed) {
    g_char <- as.character(g)
    t_char <- as.character(t)
    gs <- unlist(strsplit(g_char, ";"))

    if (any(gs %in% allow)) {
      return(list(type = t_char, gene = g_char))
    }
    # [Fix] Use NA_character_ for strict type safety
    return(list(type = if (down) "E" else t_char, gene = NA_character_))
  }

  a1 <- mapply(check_anc, loop_df$anchor1_gene, loop_df$anchor1_type, SIMPLIFY = FALSE)
  a2 <- mapply(check_anc, loop_df$anchor2_gene, loop_df$anchor2_type, SIMPLIFY = FALSE)

  loop_df$functional_anchor1_type <- vapply(a1, function(x) x$type, FUN.VALUE = character(1))
  loop_df$functional_anchor2_type <- vapply(a2, function(x) x$type, FUN.VALUE = character(1))

  reclass <- function(t1, t2, s1, s2) {
    code <- paste(sort(c(t1, t2)), collapse = "-")
    genes <- c()
    if (code == "P-P") genes <- c(s1, s2) else if (grepl("P", code)) genes <- if (t1 == "P") s1 else s2 else genes <- c(s1, s2)

    # Ensure character return for vapply safety
    gene_res <- paste(unique(na.omit(genes)), collapse = ";")
    if (gene_res == "") gene_res <- NA_character_
    list(type = code, gene = gene_res)
  }

  res_list <- mapply(
    reclass,
    loop_df$functional_anchor1_type,
    loop_df$functional_anchor2_type,
    vapply(a1, function(x) x$gene, FUN.VALUE = character(1)),
    vapply(a2, function(x) x$gene, FUN.VALUE = character(1)),
    SIMPLIFY = FALSE
  )
  loop_df$functional_type <- vapply(res_list, function(x) x$type, FUN.VALUE = character(1))
  loop_df$functional_genes <- vapply(res_list, function(x) x$gene, FUN.VALUE = character(1))

  # --- 3. 聚合 Cluster 信息 ---
  message(">>> [Step 3] Updating Clusters...")

  agg_cluster <- loop_df %>%
    filter(!is.na(cluster_id)) %>%
    group_by(cluster_id) %>%
    summarise(
      expressed_genes_Total = paste(unique(na.omit(unlist(strsplit(as.character(functional_genes), ";")))), collapse = ";")
    ) %>%
    ungroup()

  clust_info <- annotation_res$cluster_annotation
  if ("GENENAME" %in% colnames(clust_info)) clust_info <- rename(clust_info, Gene_description = GENENAME)
  clust_info <- left_join(clust_info, agg_cluster, by = "cluster_id")

  # --- 3.5 Target Integration ---
  message(">>> [Step 3.5] Integrating Targets (ID-based)...")

  bed_info <- NULL
  target_connected_loops <- NULL

  # Define DB names for later use
  txdb_pkg_name <- if (grepl("mm", species)) "TxDb.Mmusculus.UCSC.mm10.knownGene" else "TxDb.Hsapiens.UCSC.hg38.knownGene"
  org_db_name <- if (grepl("mm", species)) "org.Mm.eg.db" else "org.Hs.eg.db"

  if (!is.null(target_bed) && file.exists(target_bed)) {
    bed_raw <- read.table(target_bed, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    if(ncol(bed_raw) >= 3) {
      colnames(bed_raw)[c(1,2,3)] <- c("chr", "start", "end")
    } else {
      stop("Target BED file must have at least 3 columns.")
    }

    bed_raw$input_id <- paste0("Peak_", seq_len(nrow(bed_raw)))
    gr_target <- makeGRangesFromDataFrame(bed_raw, keep.extra.columns = TRUE)

    # Load Annotation DB
    if(!requireNamespace(txdb_pkg_name, quietly=TRUE)) {
      message("TxDb package ", txdb_pkg_name, " not installed. Skipping detailed annotation.")
      bed_annot <- as.data.frame(gr_target)
    } else {
      bed_annot <- ChIPseeker::annotatePeak(gr_target, TxDb = get(txdb_pkg_name), verbose = FALSE)
      bed_annot <- as.data.frame(bed_annot)
    }

    bed_info <- format_annotation_columns(bed_annot)
    if ("GENENAME" %in% colnames(bed_info)) bed_info <- rename(bed_info, Gene_description = GENENAME)

    if (!is.null(annotation_res$cluster_annotation)) {
      gr_clusters <- makeGRangesFromDataFrame(annotation_res$cluster_annotation, keep.extra.columns = TRUE)
      hits <- findOverlaps(gr_target, gr_clusters)

      if (length(hits) > 0) {
        hit_cids <- unique(gr_clusters$cluster_id[subjectHits(hits)])
        target_connected_loops <- loop_df %>%
          filter(cluster_id %in% hit_cids) %>%
          mutate(loop_type = functional_type)

        hit_df <- data.frame(qid = queryHits(hits), cluster_id = gr_clusters$cluster_id[subjectHits(hits)])
        hit_df$input_id <- gr_target$input_id[hit_df$qid]

        peak_genes <- hit_df %>%
          inner_join(agg_cluster %>% select(cluster_id, expressed_genes_Total), by = "cluster_id") %>%
          filter(expressed_genes_Total != "") %>%
          group_by(input_id) %>%
          summarise(
            all_annotated_gene = paste(unique(na.omit(unlist(strsplit(as.character(expressed_genes_Total), ";")))), collapse = ";")
          )

        message("   [DEBUG] Peaks connected to genes: ", nrow(peak_genes))
        if ("all_annotated_gene" %in% colnames(bed_info)) bed_info$all_annotated_gene <- NULL
        bed_info <- bed_info %>% left_join(peak_genes, by = "input_id")
      } else {
        message("    [INFO] No overlaps found.")
        bed_info$all_annotated_gene <- NA
      }
    }
  } else {
    if(!is.null(target_bed)) warning("    Target BED file path invalid.")
    bed_info <- annotation_res$target_annotation
    if (!is.null(bed_info)) bed_info$all_annotated_gene <- NA
  }

  # --- 4. 可视化 ---
  message(">>> [Step 4] Visualization...")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE) # Ensure directory exists

  plot_df <- loop_df %>% mutate(loop_type = functional_type, loop_genes = functional_genes)

  if (!exists("get_colors")) stop("Helper 'get_colors' missing.")
  colors <- get_colors(length(unique(plot_df$loop_type)), color_palette)
  names(colors) <- sort(unique(plot_df$loop_type))

  if (exists("draw_comparison_bar")) try(draw_comparison_bar(annotation_res$loop_annotation, plot_df, file.path(out_dir, paste0(project_name, "_Comparison_Bar.pdf")), colors), silent = TRUE)
  if (exists("draw_rose_plot")) try(draw_rose_plot(plot_df, project_name, file.path(out_dir, paste0(project_name, "_Filtered_Rose.pdf")), colors), silent = TRUE)
  if (exists("draw_circular_bar_plot")) try(draw_circular_bar_plot(plot_df, project_name, file.path(out_dir, paste0(project_name, "_Filtered_Circular.pdf")), colors), silent = TRUE)

  # === [NEW] KaryoploteR Visualizations ===
  message("    Plotting KaryoploteR Maps...")

  # Ensure TxDb object is loaded for Karyo
  if(requireNamespace(txdb_pkg_name, quietly=TRUE)) {
    txdb_obj <- get(txdb_pkg_name)

    # 1. Basic Gene Load (Karyo)
    unique_loop_genes <- plot_df %>%
      dplyr::filter(!is.na(loop_genes) & loop_genes != "") %>%
      tidyr::separate_rows(loop_genes, sep = ";") %>%
      dplyr::pull(loop_genes) %>%
      unique() %>%
      trimws()
    unique_loop_genes <- unique_loop_genes[unique_loop_genes != ""]

    if (length(unique_loop_genes) > 0 && exists("draw_karyo_heatmap_internal")) {
      all_genes_gr <- GenomicFeatures::genes(txdb_obj)
      # Map Symbols
      if(requireNamespace(org_db_name, quietly=TRUE)) {
        map <- AnnotationDbi::select(get(org_db_name), keys = as.character(S4Vectors::mcols(all_genes_gr)$gene_id), columns = "SYMBOL", keytype = "ENTREZID")
        S4Vectors::mcols(all_genes_gr)$SYMBOL <- map$SYMBOL[match(S4Vectors::mcols(all_genes_gr)$gene_id, map$ENTREZID)]
        target_genes_gr <- all_genes_gr[S4Vectors::mcols(all_genes_gr)$SYMBOL %in% unique_loop_genes]

        try(draw_karyo_heatmap_internal(target_genes_gr, "Gene Load (Basic)", file.path(out_dir, paste0(project_name, "_Basic_Karyo_Genes.pdf")), karyo_bin_size, 0.99, txdb_obj, species, "Genes"), silent = TRUE)
      }
    }

    # 2. Basic Anchor Load (Karyo)
    a1_coords <- plot_df %>% dplyr::select(chr = chr1, start = start1, end = end1)
    a2_coords <- plot_df %>% dplyr::select(chr = chr2, start = start2, end = end2)
    all_anchors <- dplyr::bind_rows(a1_coords, a2_coords) %>% dplyr::distinct()

    if (nrow(all_anchors) > 0 && exists("draw_karyo_heatmap_internal")) {
      gr_anchors_plot <- GenomicRanges::makeGRangesFromDataFrame(all_anchors)
      try(draw_karyo_heatmap_internal(gr_anchors_plot, "Loop Anchor Load (Basic)", file.path(out_dir, paste0(project_name, "_Basic_Karyo_Anchors.pdf")), karyo_bin_size, 0.99, txdb_obj, species, "Anchors"), silent = TRUE)
    }
  } else {
    message("    [INFO] TxDb package missing, skipping KaryoploteR plots.")
  }

  # === Target Plots (Inline Logic) ===
  if (!is.null(bed_info)) {
    message("    Plotting Target Connectivity...")

    try({
      t_col <- if ("all_annotated_gene" %in% colnames(bed_info)) "all_annotated_gene" else "loop_genes_Total"
      if (t_col %in% colnames(bed_info)) {
        p_data <- bed_info
        p_data[[t_col]] <- as.character(p_data[[t_col]])
        p_data[[t_col]][is.na(p_data[[t_col]])] <- ""

        summ <- p_data %>%
          mutate(is_con = nchar(.data[[t_col]]) > 0) %>%
          mutate(Status = ifelse(is_con, "Connected", "Orphan")) %>%
          group_by(Status) %>%
          summarise(Count = n()) %>%
          mutate(Pct = Count / sum(Count))

        p <- ggplot(summ, aes(x = Status, y = Count, fill = Status)) +
          geom_col(width = 0.6, color = "black") +
          geom_text(aes(label = paste0(Count, "\n(", scales::percent(Pct, 0.1), ")")), vjust = -0.5, fontface = "bold") +
          scale_fill_brewer(palette = "Pastel1") +
          theme_classic() +
          labs(title = paste0(project_name, ": Connectivity")) +
          scale_y_continuous(expand = expansion(mult = c(0, 0.2)))
        ggsave(file.path(out_dir, paste0(project_name, "_Refined_Target_Connectivity.pdf")), p, width = 5, height = 6)
      }
    })

    if (!is.null(target_connected_loops) && nrow(target_connected_loops) > 0) {
      message("    Plotting Target Donut...")
      try({
        df <- target_connected_loops %>%
          group_by(loop_type) %>%
          summarise(Count = n()) %>%
          mutate(Pct = Count / sum(Count)) %>%
          arrange(desc(loop_type))
        df$ymax <- cumsum(df$Pct)
        df$ymin <- c(0, head(df$ymax, n = -1))
        df$labelPos <- (df$ymax + df$ymin) / 2
        df$labelTxt <- paste0(df$loop_type, "\n", df$Count, "\n(", scales::percent(df$Pct, 0.1), ")")

        p <- ggplot(df, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3, fill = loop_type)) +
          geom_rect(color = "white") +
          geom_text(x = 3.5, aes(y = labelPos, label = labelTxt), size = 3, fontface = "bold") +
          coord_polar(theta = "y") +
          xlim(c(2, 4.5)) +
          scale_fill_manual(values = colors) +
          theme_void() +
          annotate("text", x = 0, y = 0, label = paste0("Total\n", sum(df$Count)), size = 5, fontface = "bold")
        ggsave(file.path(out_dir, paste0(project_name, "_Refined_Target_Loop_Donut.pdf")), p, width = 6, height = 6)
      })
    }
  }


  message(">>> [Step 5] Exporting...")
  tryCatch({
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, "Filtered Loop Annotation")
    openxlsx::writeData(wb, "Filtered Loop Annotation", plot_df)
    openxlsx::addWorksheet(wb, "Cluster Annotation")
    openxlsx::writeData(wb, "Cluster Annotation", clust_info)
    if (!is.null(bed_info)) {
      openxlsx::addWorksheet(wb, "Target Annotation")
      openxlsx::writeData(wb, "Target Annotation", bed_info)
    }
    openxlsx::saveWorkbook(wb, file.path(out_dir, paste0(project_name, "_Filtered_Results.xlsx")), overwrite = TRUE)
    message("    Excel saved.")
  },
  error = function(e) message("    Export Failed: ", e$message)
  )

  message("Done. Output: ", out_dir)
  return(list(loop_annotation = plot_df, cluster_annotation = clust_info, target_annotation = bed_info))
}
