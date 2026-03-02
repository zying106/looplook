#' Internal Package Imports
#'
#' @name looplook_imports
#' @noRd
#' @importFrom grDevices col2rgb rgb
#' @importFrom stats fisher.test median na.omit p.adjust quantile reorder runif setNames t.test var wilcox.test
#' @importFrom utils head install.packages read.table write.csv
NULL

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    ".", "All_Anchor_Genes", "All_Loop_Connected_Genes", "Assigned_Target_Genes",
    "CleanLoopType", "Conn_Group", "Count", "Degree", "Description",
    "Description_unique", "Distal_Anchor_ID", "Dominant_Interaction",
    "Dominant_Interaction_Filtered", "Expression", "Expression_Status", "FDR",
    "Feature", "Filtered", "Final_Label", "Fraction", "GENENAME", "Gene",
    "Genomic_Distribution", "Group", "High_Connectivity_Gene",
    "Is_Active_Gene", "Is_High_Connectivity_Distal_Element",
    "Is_High_Connectivity_Gene", "Is_High_Distal_Connectivity_Gene", "L1_Raw",
    "L2_Raw", "L3_Raw", "LFC", "Label", "LabelText", "Label_Text",
    "Linked_Loop_IDs", "Log10Degree", "LogFDR", "LogP", "Loop_Type",
    "Mean_Expression_Temp", "MotifLabel", "ONTOLOGY", "OddsRatio", "Original",
    "Percentage", "PlotFamily", "Putative_Target_Genes", "Rank",
    "Regulated_promoter_genes", "SANKEY_RAW_GENES", "SYMBOL", "SampleID",
    "Simplified", "Source", "Stage", "Status", "Target_Genes",
    "Target_Genes_Filtered", "Total_Loops", "Total_Loops_Filtered",
    "Unique_Gene_Count", "a1_id", "a2_id", "all_cluster_loop_genes", "all_of",
    "anchor1_gene", "anchor1_source", "anchor1_type", "anchor2_gene",
    "anchor2_source", "anchor2_type", "anchor_id", "annotation", "chr",
    "cluster_id", "col2rgb", "combined_score", "count", "deg", "detail_anno",
    "elementNROWS", "everything", "expansion", "final_color", "final_fill",
    "final_label", "final_symbols", "fisher.test", "fraction",
    "functional_anchor1_type", "functional_anchor2_type", "geneList",
    "gene_id", "gene_level", "geom_hline", "group", "has_active", "head",
    "hjust", "install.packages", "is_e_type", "is_lower_e", "label",
    "labelPosition", "label_text", "label_x", "len", "lfc", "linked_loops",
    "log_expr", "logP", "loop_ID", "loop_genes", "loop_genes_Total", "loop_i",
    "loop_type", "median", "mid1", "mid2", "n", "n_Linked_Distal",
    "n_Linked_Distal_Filtered", "n_Linked_Promoters",
    "n_Linked_Promoters_Filtered", "na.omit", "name", "p.adjust", "plot_label",
    "prop", "proximate_loop_gene", "pvalue", "qid", "quantile", "query_idx",
    "read.table", "reg_loop_genes", "reorder", "rgb", "runif", "runningScore",
    "scale_color_identity", "scale_fill_identity", "setNames",
    "single_loop_genes", "strand", "t.test", "t1", "t2", "tgt_genes_p",
    "tgt_genes_pg", "tgt_genes_prio", "topo_genes_p", "topo_genes_pg", "tpm",
    "tx_id", "type", "type_code", "type_rank", "valid_genes", "valid_tpms",
    "var", "width", "wilcox.test", "write.csv", "y_mid", "ymax", "ymin", "ypos", ":=",
    "Loop_Connection", "Neighbor_Gene", "Neighbor_Type", "s1", "s2", "x", "y"
  ))
}

#' Internal: Generate Colors
#'
#' Helper to generate a vector of n colors.
#'
#' @param n Integer. Number of colors.
#' @param palette_input Character. Palette name or custom colors.
#' @return Hex color codes.
#'
#' @importFrom RColorBrewer brewer.pal.info brewer.pal
#' @importFrom grDevices colorRampPalette
#' @importFrom scales hue_pal
#' @keywords internal
get_colors <- function(n, palette_input) {
  if (length(palette_input) == 1 && palette_input %in% row.names(RColorBrewer::brewer.pal.info)) {
    max_avail <- RColorBrewer::brewer.pal.info[palette_input, "maxcolors"]
    pal <- RColorBrewer::brewer.pal(min(n, max_avail), palette_input)
    return(grDevices::colorRampPalette(pal)(n))
  } else if (length(palette_input) >= 1) {
    if (length(palette_input) < n) {
      return(rep(palette_input, length.out = n))
    }
    return(palette_input[seq_len(n)])
  } else {
    return(scales::hue_pal()(n))
  }
}


#' Internal: Draw Karyotype Heatmap
#'
#' Creates a genome-wide heatmap of genomic feature density (e.g., loops) across chromosomes,
#' binned by a fixed window size, and rendered as a karyotype plot.
#'
#' @param gr_data (GRanges) Genomic ranges to visualize (e.g., loop anchors).
#' @param title_prefix (character) Subtitle descriptor (e.g., sample name).
#' @param filename (character) Output PDF path.
#' @param bin_size (integer) Bin width in base pairs (e.g., 1e6 for 1 Mb).
#' @param sat_level (numeric) Quantile (0–1) for color saturation (e.g., 0.95).
#' @param ref_txdb (TxDb or similar) Reference genome annotation for chromosome lengths.
#' @param plot_species (character) Genome build/species code (e.g., "hg38", "mm10").
#' @param unit_label (character) Unit for load annotation (e.g., "loops").
#' @keywords internal
#' @importFrom GenomeInfoDb seqinfo seqlevelsStyle keepSeqlevels seqlengths seqlevels
#' @importFrom GenomicRanges GRanges tileGenome countOverlaps
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors mcols
#' @importFrom grDevices pdf dev.off colorRampPalette
#' @importFrom karyoploteR plotKaryotype kpRect getDefaultPlotParams
#' @importFrom fields image.plot
#' @return Invisibly returns `NULL`. Side effect: saves a PDF karyotype heatmap to `filename`.
draw_karyo_heatmap_internal <- function(gr_data, title_prefix, filename, bin_size, sat_level, ref_txdb, plot_species, unit_label, custom_colors = NULL) {
  standard_chroms <- paste0("chr", c(seq_len(22), "X", "Y"))
  if (grepl("mm", plot_species)) standard_chroms <- paste0("chr", c(seq_len(19), "X", "Y"))

  std_seqinfo <- GenomeInfoDb::seqinfo(ref_txdb)
  try(
    {
      GenomeInfoDb::seqlevelsStyle(gr_data) <- "UCSC"
    },
    silent = TRUE
  )

  existing <- intersect(GenomeInfoDb::seqlevels(gr_data), standard_chroms)
  if (length(existing) == 0) {
    return(invisible(NULL))
  }

  gr_data <- GenomeInfoDb::keepSeqlevels(gr_data, existing, pruning.mode = "coarse")
  GenomeInfoDb::seqlevels(gr_data) <- standard_chroms
  common <- intersect(GenomeInfoDb::seqlevels(gr_data), GenomeInfoDb::seqlevels(std_seqinfo))
  GenomeInfoDb::seqlengths(gr_data)[common] <- GenomeInfoDb::seqlengths(std_seqinfo)[common]
  valid_chroms <- intersect(standard_chroms, GenomeInfoDb::seqlevels(std_seqinfo))

  if (length(valid_chroms) > 0) {
    full_genome_gr <- GenomicRanges::GRanges(seqnames = valid_chroms, ranges = IRanges::IRanges(start = 1, end = GenomeInfoDb::seqlengths(std_seqinfo)[valid_chroms]))
    GenomeInfoDb::seqinfo(full_genome_gr) <- std_seqinfo[valid_chroms]
    tiles <- GenomicRanges::tileGenome(GenomeInfoDb::seqinfo(full_genome_gr), tilewidth = bin_size, cut.last.tile.in.chrom = TRUE)
    hits <- GenomicRanges::countOverlaps(tiles, gr_data)

    bin_size_mb <- bin_size / 1e6
    median_val <- median(hits[hits > 0], na.rm = TRUE)
    if (is.na(median_val)) median_val <- 0


    if (is.null(custom_colors)) {
      heatmap_colors <- c("#FFFFFF", "#FFFFCC", "#FFEDA0", "#FED976", "#FEB24C", "#FD8D3C", "#FC4E2A", "#E31A1C", "#BD0026", "#800026", "#000000")
    } else {
      heatmap_colors <- custom_colors
    }

    if (max(hits) == 0) {
      max_load <- 0
      S4Vectors::mcols(tiles)$color <- "white"
      cols <- c("white")
    } else {
      cutoff <- as.numeric(quantile(hits[hits > 0], probs = sat_level, names = FALSE))
      if (is.na(cutoff) || cutoff < 1) cutoff <- max(hits)
      max_load <- round(cutoff / bin_size_mb, 1)
      capped <- ifelse(hits > cutoff, cutoff, hits)

      col_func <- grDevices::colorRampPalette(heatmap_colors)
      cols <- col_func(100)

      idx <- ceiling((capped / cutoff) * 99) + 1
      idx[hits == 0] <- 1
      S4Vectors::mcols(tiles)$color <- cols[idx]
    }

    grDevices::pdf(filename, width = 10, height = 12)
    graphics::par(oma = c(2, 2, 6, 2))
    pp <- karyoploteR::getDefaultPlotParams(plot.type = 1)
    pp$leftmargin <- 0.08
    pp$rightmargin <- 0.08
    pp$data1height <- 100

    kp <- karyoploteR::plotKaryotype(genome = plot_species, plot.type = 1, chromosomes = valid_chroms, plot.params = pp, main = NULL)
    karyoploteR::kpRect(kp, data = tiles, y0 = 0, y1 = 1, col = S4Vectors::mcols(tiles)$color, border = NA)

    main_title <- paste0("Loop Analysis: ", title_prefix, "\n(Genomic Load: Median ~", round(median_val / bin_size_mb, 1), " ", unit_label, "/MB)")
    graphics::mtext(main_title, side = 3, line = 1, outer = TRUE, cex = 1.2, font = 2)

    if (requireNamespace("fields", quietly = TRUE)) {
      fields::image.plot(
        legend.only = TRUE, zlim = c(0, max_load), col = cols,
        legend.lab = paste0("Load (", unit_label, "/MB)"), legend.mar = 4.5, smallplot = c(0.88, 0.91, 0.3, 0.7)
      )
    }
    grDevices::dev.off()
    message("    Saved Heatmap: ", filename)
  }
}

#' Draw Expression Violin Plot
#'
#' Creates a violin + boxplot showing log2-transformed gene expression
#' grouped by loop type (e.g., promoter-enhancer, enhancer-enhancer).
#'
#' @param plot_data (data.frame) Must contain columns: `loop_type`, `expression_value`.
#' @param project_name (character) Project or sample name for plot title.
#' @param filename (character) Output file path (e.g., "expr_violin.pdf").
#' @param unit_type (character) Expression unit (e.g., "TPM", "FPKM").
#' @param group_colors (character vector) Named or ordered colors for each `loop_type`.
#' @keywords internal
#' @importFrom ggplot2 ggplot aes geom_violin geom_boxplot scale_fill_manual theme_minimal theme element_text element_blank labs ggsave
#' @return A `ggplot` object of the violin + boxplot.
draw_expression_violin <- function(plot_data, project_name, filename, unit_type, group_colors) {
  # 1. Log transform
  plot_data$log_expr <- log2(plot_data$expression_value + 0.01)


  plot_data_unique <- plot_data %>%
    dplyr::select(loop_type, loop_genes, log_expr) %>%
    dplyr::distinct(loop_type, loop_genes, .keep_all = TRUE)


  count_df <- plot_data_unique %>%
    dplyr::group_by(loop_type) %>%
    dplyr::summarise(n = dplyr::n()) %>%
    dplyr::ungroup()


  new_labels <- stats::setNames(
    paste0(count_df$loop_type, "\n(n=", count_df$n, ")"),
    count_df$loop_type
  )


  p <- ggplot2::ggplot(plot_data_unique, ggplot2::aes(x = loop_type, y = log_expr, fill = loop_type)) +
    ggplot2::geom_violin(scale = "width", trim = FALSE, alpha = 0.7, color = "grey30", size = 0.3) +
    ggplot2::geom_boxplot(width = 0.15, fill = "white", color = "black", outlier.shape = NA, alpha = 0.9) +
    ggplot2::scale_fill_manual(values = group_colors) +
    ggplot2::scale_x_discrete(labels = new_labels) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5, face = "bold", color = "black", size = 11),
      panel.grid.major.x = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = "none"
    ) +
    ggplot2::labs(
      title = paste0(project_name, ": Gene Expression by Loop Type"),
      x = NULL,
      y = paste0("Expression (Log2 ", unit_type, ")")
    )

  ggplot2::ggsave(filename, p, width = 8, height = 6)
  message("    Saved (Violin): ", filename)
  return(p)
}

#' Draw Enhancer Source Distribution
#'
#' Bar plot showing the origin of functional enhancer anchors
#' (e.g., native vs. promoter-derived).
#'
#' @param loop_data (data.frame) Must contain:
#'   - `functional_anchor1_type`, `anchor1_source`
#'   - `functional_anchor2_type`, `anchor2_source`
#' @param project_name (character) Project name for title.
#'   Only anchors where type == "E" are considered.
#' @param filename (character) Output file path (e.g., "enh_sources.pdf").
#' @keywords internal
#' @importFrom dplyr select filter group_by summarise ungroup mutate arrange bind_rows
#' @importFrom ggplot2 ggplot aes geom_bar geom_text scale_fill_manual scale_y_continuous theme_classic theme element_text element_blank labs ggsave
#' @importFrom scales percent
#' @return A `ggplot` object of the bar plot showing enhancer source distribution.
draw_enhancer_source_distribution <- function(loop_data, project_name, filename) {
  a1 <- loop_data %>% dplyr::select(type = functional_anchor1_type, source = anchor1_source)
  a2 <- loop_data %>% dplyr::select(type = functional_anchor2_type, source = anchor2_source)
  plot_data <- dplyr::bind_rows(a1, a2) %>%
    dplyr::filter(type == "E") %>%
    dplyr::group_by(source) %>%
    dplyr::summarise(Count = dplyr::n(), .groups = "drop") %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Percentage = Count / sum(Count)) %>%
    dplyr::arrange(dplyr::desc(Count))

  source_colors <- c(
    "Native" = "#E0E0E0",
    "Promoter-derived enhancer" = "#D95F02",
    "Gene_body-derived enhancer" = "#66A61E"
  )
  avail_sources <- unique(plot_data$source)
  final_colors <- source_colors[avail_sources]
  if (any(is.na(final_colors))) {
    final_colors <- get_colors(length(avail_sources), "Set2")
    names(final_colors) <- avail_sources
  }

  plot_data$Label <- paste0(plot_data$Count, "\n(", scales::percent(plot_data$Percentage, accuracy = 0.1), ")")

  p <- ggplot(plot_data, aes(x = reorder(source, -Count), y = Count, fill = source)) +
    geom_bar(stat = "identity", width = 0.7, color = "black", alpha = 0.8) +
    geom_text(aes(label = Label), vjust = -0.5, fontface = "bold", size = 4) +
    scale_fill_manual(values = final_colors) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    theme_classic(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1, face = "bold", color = "black"),
      axis.title.x = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "none"
    ) +
    labs(
      title = paste0(project_name, ": Origin of Functional Enhancers"),
      y = "Number of Anchors",
      subtitle = "Breakdown of anchors classified as Enhancers after filtering"
    )
  ggsave(filename, p, width = 7, height = 6)
  message("    Saved (Enhancer Sources): ", filename)
  return(p)
}


#' Internal: Draw Rose Chart (Loop Counts)
#'
#' Visualizes the proportion of loop types based on LOOP COUNTS.
#'
#' @param data_df Data frame containing loop information.
#' @param project_name Character string for the project title.
#' @param filename Character string for the output filename.
#' @param color_vec Named character vector for colors.
#'
#' @return A ggplot object representing the rose plot.
#'
#' @keywords internal
draw_rose_plot <- function(data_df, project_name, filename, color_vec) {
  rose_data <- data_df %>%
    dplyr::group_by(loop_type) %>%
    dplyr::summarise(Count = dplyr::n()) %>%
    dplyr::mutate(Percentage = Count / sum(Count)) %>%
    dplyr::arrange(dplyr::desc(Count))
  rose_data$Label <- paste0(rose_data$Count, "\n(", scales::percent(rose_data$Percentage, 0.1), ")")

  p <- ggplot2::ggplot(rose_data, ggplot2::aes(x = reorder(loop_type, -Count), y = Count, fill = loop_type)) +
    ggplot2::geom_bar(stat = "identity", width = 1, color = "white") +
    ggplot2::coord_polar(theta = "x", start = 0) +
    ggplot2::scale_fill_manual(values = color_vec) +
    ggplot2::geom_text(ggplot2::aes(y = Count, label = Label), size = 3.5, color = "black", vjust = -0.5) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(axis.title = ggplot2::element_blank(), axis.text = ggplot2::element_blank(), panel.grid = ggplot2::element_blank(), legend.position = "right") +
    ggplot2::labs(title = paste0(project_name, ": Loop Proportion (By Count)"))
  ggplot2::ggsave(filename, p, width = 8, height = 7)
  message("    Saved: ", filename)
}

#' Internal: Draw Circular Bar Plot (Gene Counts)
#'
#' Optimized for strictly vertical label alignment using y=0 anchor.
#'
#' @param data_df Data frame containing loop and gene information.
#' @param project_name Character string for the project title.
#' @param filename Character string for the output filename.
#' @param color_vec Named character vector for colors.
#'
#' @return A ggplot object representing the circular bar plot.
#'
#' @keywords internal
draw_circular_bar_plot <- function(data_df, project_name, filename, color_vec) {
  circ_data <- data_df %>%
    dplyr::filter(!is.na(loop_genes) & loop_genes != "") %>%
    tidyr::separate_rows(loop_genes, sep = ";") %>%
    dplyr::group_by(loop_type) %>%
    dplyr::summarise(Unique_Gene_Count = dplyr::n_distinct(trimws(loop_genes))) %>%
    dplyr::arrange(Unique_Gene_Count) %>%
    dplyr::mutate(Label_Text = paste0(loop_type, " : ", Unique_Gene_Count))

  circ_data$loop_type <- factor(circ_data$loop_type, levels = circ_data$loop_type)
  if (nrow(circ_data) == 0) {
    return(NULL)
  }
  max_gene_count <- max(circ_data$Unique_Gene_Count, na.rm = TRUE)

  p <- ggplot2::ggplot(circ_data, ggplot2::aes(x = loop_type, fill = loop_type)) +
    ggplot2::geom_col(ggplot2::aes(y = max_gene_count), width = 0.05, fill = "grey92", color = NA) +
    ggplot2::geom_col(ggplot2::aes(y = Unique_Gene_Count), width = 0.8, color = "white", size = 0.2) +
    ggplot2::geom_text(ggplot2::aes(y = Unique_Gene_Count + max_gene_count * 0.02, label = Label_Text), hjust = 0, size = 3.5, fontface = "bold") +
    ggplot2::coord_polar(theta = "y", start = 0, clip = "off") +
    ggplot2::scale_fill_manual(values = color_vec, name = "Loop Type") +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(axis.title = ggplot2::element_blank(), axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), panel.grid = ggplot2::element_blank(), legend.position = "right", plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"), plot.margin = ggplot2::margin(t = 20, r = 100, b = 20, l = 20, unit = "pt")) +
    ggplot2::scale_y_continuous(limits = c(-max_gene_count * 0.4, max_gene_count * 1.3)) +
    ggplot2::labs(title = paste0(project_name, ": Unique Target Genes (Ascending)"))

  ggplot2::ggsave(filename, p, width = 10, height = 7)
  message(" Saved: ", filename)
}

#' Internal: Draw Comparison Bar Chart
#'
#' Visualizes the comparison of loop counts between original and filtered datasets.
#'
#' @param original_df Data frame containing the original loop data.
#' @param filtered_df Data frame containing the filtered loop data.
#' @param filename Character string for the output filename.
#' @param color_vec Named character vector for colors.
#'
#' @return A ggplot object visualizing the comparison between groups.
#'
#' @keywords internal
draw_comparison_bar <- function(original_df, filtered_df, filename, color_vec) {
  fmt_type <- function(x) {
    return(x)
  }
  df_compare <- dplyr::bind_rows(
    original_df %>% dplyr::mutate(loop_type = fmt_type(loop_type)) %>% dplyr::group_by(loop_type) %>% dplyr::summarise(Count = dplyr::n()) %>% dplyr::mutate(Stage = "Original"),
    filtered_df %>% dplyr::group_by(loop_type) %>% dplyr::summarise(Count = dplyr::n()) %>% dplyr::mutate(Stage = "Filtered")
  ) %>% dplyr::mutate(Stage = factor(Stage, levels = c("Original", "Filtered")))

  p <- ggplot2::ggplot(df_compare, ggplot2::aes(x = loop_type, y = Count, fill = loop_type, alpha = Stage)) +
    ggplot2::geom_bar(stat = "identity", position = "dodge", color = "black") +
    ggplot2::scale_fill_manual(values = color_vec) +
    ggplot2::scale_alpha_manual(values = c(0.4, 1.0)) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::labs(title = "Loop Classification Changes", y = "Number of Loops", x = "Loop Type")
  ggplot2::ggsave(filename, p, width = 8, height = 6)
  message("    Saved: ", filename)
}

#' Internal: Draw Target Genomic Annotation Distribution (Pie Chart)
#'
#' Optimized: Labels only show Count (%), and hides labels for small slices (<2%) to avoid overlap.
#'
#' @param bed_info Data frame containing annotation information.
#' @param project_name Character string for the project title.
#' @param filename Character string for the output filename.
#' @param color_palette Character string for the color palette (default: "Set3").
#'
#' @return A ggplot object representing the annotation pie chart.
#'
#' @keywords internal
draw_target_annotation_pie <- function(bed_info, project_name, filename, color_palette = "Set3") {
  plot_data <- bed_info %>%
    dplyr::mutate(Feature = gsub(" \\(.*", "", annotation)) %>%
    dplyr::mutate(Feature = ifelse(grepl("Promoter", Feature), "Promoter", Feature)) %>%
    dplyr::group_by(Feature) %>%
    dplyr::summarise(Count = dplyr::n()) %>%
    dplyr::mutate(Percentage = Count / sum(Count)) %>%
    dplyr::arrange(dplyr::desc(Count))

  plot_data <- plot_data %>%
    dplyr::mutate(
      Label_Text = paste0(Count, "\n(", scales::percent(Percentage, 0.1), ")"),
      Final_Label = ifelse(Percentage >= 0.02, Label_Text, "")
    )


  n_groups <- nrow(plot_data)
  if (!exists("get_colors")) custom_colors <- scales::hue_pal()(n_groups) else custom_colors <- get_colors(n_groups, color_palette)


  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = "", y = Count, fill = Feature)) +
    ggplot2::geom_bar(stat = "identity", width = 1, color = "white", size = 0.5) +
    ggplot2::coord_polar("y", start = 0) +
    ggplot2::geom_text(
      ggplot2::aes(label = Final_Label),
      position = ggplot2::position_stack(vjust = 0.5),
      size = 3.5,
      fontface = "bold",
      color = "black"
    ) +
    ggplot2::scale_fill_manual(values = custom_colors, name = "Genomic Feature") +
    ggplot2::theme_void(base_size = 14) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", margin = ggplot2::margin(b = 10)),
      legend.position = "right",
      legend.title = ggplot2::element_text(face = "bold")
    ) +
    ggplot2::labs(title = paste0(project_name, ": Target Peak Genomic Distribution"))

  ggplot2::ggsave(filename, p, width = 9, height = 6)
  message("    Saved (Target Pie): ", filename)
}

#' Internal: Draw Target-Loop Connectivity
#'
#' Visualizes how many target peaks overlap with loops, and breakdown by Loop Type.
#'
#' @param bed_info Data frame containing peak information.
#' @param cluster_info Data frame containing cluster information (currently unused in plot but kept for consistency).
#' @param project_name Character string for the project title.
#' @param filename Character string for the output filename.
#' @param color_palette Character string for the color palette (default: "Set2").
#'
#' @return A ggplot object representing the connectivity bar chart.
#'
#' @keywords internal
draw_target_connectivity_bar <- function(bed_info, cluster_info, project_name, filename, color_palette = "Set2") {
  summary_df <- bed_info %>%
    dplyr::mutate(Status = ifelse(!is.na(loop_genes_Total) & loop_genes_Total != "", "Connected (Has Loop)", "Orphan (No Loop)")) %>%
    dplyr::group_by(Status) %>%
    dplyr::summarise(Count = dplyr::n()) %>%
    dplyr::mutate(Percentage = Count / sum(Count))

  p1 <- ggplot2::ggplot(summary_df, ggplot2::aes(x = Status, y = Count, fill = Status)) +
    ggplot2::geom_col(width = 0.6, color = "black") +
    ggplot2::geom_text(ggplot2::aes(label = Count), vjust = -0.5, fontface = "bold") +
    ggplot2::scale_fill_brewer(palette = "Pastel1") +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(legend.position = "none", axis.title.x = ggplot2::element_blank()) +
    ggplot2::labs(title = "Target Connectivity", y = "Number of Peaks")

  ggplot2::ggsave(filename, p1, width = 6, height = 6)
  message("    Saved (Target Connectivity): ", filename)
}


#' Internal: Draw Target-Associated Loop Distribution (Donut Chart)
#'
#' Visualizes the distribution of loop types that are connected to target peaks.
#'
#' @param loop_data Data frame containing loop information.
#' @param project_name Character string for the project title.
#' @param filename Character string for the output filename.
#' @param color_vec Named character vector for colors.
#'
#' @return A ggplot object representing the donut chart of loop type distribution.
#'
#' @keywords internal
draw_target_loop_donut <- function(loop_data, project_name, filename, color_vec) {
  plot_data <- loop_data %>%
    dplyr::group_by(loop_type) %>%
    dplyr::summarise(Count = dplyr::n()) %>%
    dplyr::mutate(Percentage = Count / sum(Count)) %>%
    dplyr::arrange(dplyr::desc(Count))

  plot_data$Label <- paste0(
    plot_data$loop_type, "\n",
    plot_data$Count, " (", scales::percent(plot_data$Percentage, 0.1), ")"
  )


  plot_data$ymax <- cumsum(plot_data$Percentage)
  plot_data$ymin <- c(0, head(plot_data$ymax, n = -1))
  plot_data$labelPosition <- (plot_data$ymax + plot_data$ymin) / 2


  p <- ggplot2::ggplot(plot_data, ggplot2::aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3, fill = loop_type)) +
    ggplot2::geom_rect(color = "white") +
    ggplot2::geom_text(x = 4.2, ggplot2::aes(y = labelPosition, label = Label), size = 3.5, color = "black") +
    ggplot2::coord_polar(theta = "y") +
    ggplot2::scale_fill_manual(values = color_vec) +
    ggplot2::xlim(c(2, 4.5)) +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "none", plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")) +
    ggplot2::labs(title = paste0(project_name, ": Loops Connected to Targets"))

  ggplot2::ggsave(filename, p, width = 7, height = 7)
  message("    Saved (Target Donut): ", filename)
}
