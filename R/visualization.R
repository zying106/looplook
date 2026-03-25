#' Publication-ready visualization toolkit for rendering genomic track data and statistical summaries related to 3D chromatin interactions
#'
#' Generates an integrative genomic track plot resembling IGV, displaying chromatin loops as arcs,
#' loop anchors as rectangles, optional overlapping features (e.g., ChIP-seq peaks), and annotated genes.
#' Loop arcs can be colored or sized by interaction score (7th column in BEDPE).
#'
#' @param bedpe_file Character. Path to a BEDPE file (at least 6 columns; 7th column used as score if present).
#' @param target_bed Optional character. Path to a BED file (e.g., peaks) to overlay below the loop track.
#' @param chr Character. Chromosome name (e.g., "chr8"). If NULL, inferred from most frequent chromosome in BEDPE.
#' @param from Numeric. Start coordinate of the region to plot.
#' @param to Numeric. End coordinate of the region to plot.
#' @param species Character. Genome assembly: "hg38", "hg19", "mm10", or "mm9".
#' @param max_levels Integer. Maximum number of vertical levels for arc stacking (default: 10).
#' @param base_anchor_height Numeric. Height of anchor rectangles (default: 0.1).
#' @param loop_color Character. Default color for arcs when no score is provided (e.g., "#5D6D7E").
#' @param anchor_color Character. Color for loop anchor rectangles (default: "#3498DB").
#' @param overlap_color Character. Color for overlap track (default: "#E74C3C").
#' @param exon_color Character. Gene exon fill color (default: "#2C3E50").
#' @param intron_color Character. Gene intron line color (default: "black").
#' @param score_to_alpha Logical. Whether to map interaction scores to arc transparency.
#' @param min_score Logical. If TRUE, use score to control arc line width instead of color (not yet implemented in current version; future extension).
#' @param save_file Optional character. File path to save the plot (e.g., "region_plot.pdf").
#' @return A \code{ggplot} object.
#' @importFrom dplyr %>%
#' @importFrom ggplot2 ggplot geom_rect geom_segment annotate coord_cartesian scale_x_continuous theme_classic labs ggsave arrow unit theme element_blank element_rect element_text margin
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggforce geom_bezier
#' @importFrom scales rescale comma
#' @importFrom GenomicFeatures genes exonsBy
#' @importFrom AnnotationDbi keys keytypes
#' @importFrom IRanges disjointBins IRanges
#' @export
#' @examples
#' # 1. Get paths to example files included in the package
#' bedpe_path <- system.file("extdata", "example_loops_1.bedpe", package = "looplook")
#' bed_path <- system.file("extdata", "example_peaks.bed", package = "looplook")
#'
#' # 2. Run plotting (requires TxDb package for gene annotation)
#' if (bedpe_path != "" &&
#'   requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE) &&
#'   requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
#'   # Example : Integrative plot with overlapping peaks and output to file
#'   p2 <- plot_peaks_interactions(
#'     bedpe_file = bedpe_path,
#'     target_bed = bed_path,
#'     chr = "chr1",
#'     from = 11884299,
#'     to = 12106581,
#'     species = "hg38",
#'     save_file = tempfile(fileext = ".pdf")
#'   )
#' }
plot_peaks_interactions <- function(
  bedpe_file,
  target_bed = NULL,
  chr = NULL,
  from = NULL,
  to = NULL,
  species = "hg38",
  max_levels = 10,
  base_anchor_height = 0.05,
  loop_color = "#5D6D7E",
  anchor_color = "#3498DB",
  overlap_color = "#02ABB4",
  exon_color = "#2C3E50",
  intron_color = "black",
  score_to_alpha = TRUE,
  min_score = NULL,
  save_file = NULL
) {
  # === 0. Configuration & Checks ===
  if (!requireNamespace("ggforce", quietly = TRUE)) stop("Please install 'ggforce'")
  if (!requireNamespace("GenomicFeatures", quietly = TRUE)) stop("Please install 'GenomicFeatures'")
  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) stop("Please install 'AnnotationDbi'")
  if (!requireNamespace("scales", quietly = TRUE)) stop("Please install 'scales'")

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
    stop("Species not supported.")
  }

  if (!requireNamespace(txdb_pkg, quietly = TRUE)) stop("Please install ", txdb_pkg)
  if (!requireNamespace(org_db_pkg, quietly = TRUE)) stop("Please install ", org_db_pkg)

  if (!requireNamespace(txdb_pkg, quietly = TRUE)) {
    stop("Package ", txdb_pkg, " is not installed.")
  }
  txdb <- utils::getFromNamespace(txdb_pkg, txdb_pkg)

  # === 1. Data Reading & Preprocessing ===
  loops_raw <- read.table(bedpe_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  colnames(loops_raw)[seq_len(6)] <- c("chr1", "start1", "end1", "chr2", "start2", "end2")

  # Score Auto-detect
  has_score <- FALSE
  if (ncol(loops_raw) >= 7) {
    is_numeric_col <- function(x) {
      nums <- as.numeric(x)
      return(sum(!is.na(nums)) / length(nums) > 0.9 && !all(is.na(nums)))
    }
    if (ncol(loops_raw) >= 8 && is_numeric_col(loops_raw[[8]])) {
      colnames(loops_raw)[8] <- "score"
      has_score <- TRUE
    } else if (is_numeric_col(loops_raw[[7]])) {
      colnames(loops_raw)[7] <- "score"
      has_score <- TRUE
    }
  }

  loops <- loops_raw %>%
    dplyr::mutate(dplyr::across(c(start1, end1, start2, end2), as.numeric),
      chr1 = as.character(chr1), chr2 = as.character(chr2)
    )

  # Score Filtering
  if (has_score) {
    loops$score <- as.numeric(loops$score)
    if (!is.null(min_score)) {
      loops <- loops %>% dplyr::filter(score >= min_score)
      if (nrow(loops) == 0) stop("All loops were filtered out by min_score.")
    }
  }

  if (is.null(chr)) {
    all_chr <- c(loops$chr1, loops$chr2)
    chr <- names(sort(table(all_chr), decreasing = TRUE))[1]
  }
  loops_chr <- loops %>% dplyr::filter(chr1 == chr2 & chr1 == chr)

  if (is.null(from)) from <- min(c(loops_chr$start1, loops_chr$start2))
  if (is.null(to)) to <- max(c(loops_chr$end1, loops_chr$end2))

  # Viewport Filter: Keep if AT LEAST ONE anchor is in window
  loops_view <- loops_chr %>%
    dplyr::filter(
      (end1 >= from & start1 <= to) | (end2 >= from & start2 <= to)
    )

  text_indent <- from + ((to - from) * 0.005)

  # === 2. Loop Data ===
  anchor_ymin <- 0
  anchor_ymax <- base_anchor_height

  bez_df <- data.frame()
  anchors <- data.frame()
  plot_ymax <- anchor_ymax + 0.5

  if (nrow(loops_view) > 0) {
    anchors <- dplyr::bind_rows(
      loops_view %>% dplyr::transmute(chr = chr1, start = start1, end = end1, loop_i = dplyr::row_number(), score = if (has_score) score else 1),
      loops_view %>% dplyr::transmute(chr = chr2, start = start2, end = end2, loop_i = dplyr::row_number(), score = if (has_score) score else 1)
    ) %>% dplyr::mutate(ymin = anchor_ymin, ymax = anchor_ymax)

    loops_calc <- loops_view %>%
      dplyr::mutate(mid1 = (start1 + end1) / 2, mid2 = (start2 + end2) / 2, loop_i = dplyr::row_number(), center = (mid1 + mid2) / 2, peak = anchor_ymax + (max_levels * 0.1) + 0.1)

    bez_list <- lapply(seq_len(nrow(loops_calc)), function(i) {
      d <- loops_calc[i, ]
      data.frame(loop_i = d$loop_i, x = c(d$mid1, d$center, d$mid2), y = c(anchor_ymax, d$peak, anchor_ymax), score = if (has_score) d$score else 1, stringsAsFactors = FALSE)
    })
    bez_df <- do.call(rbind, bez_list)

    max_loop_y <- max(bez_df$y)
    plot_ymax <- max_loop_y * 1.05

    calc_alpha_color <- function(scores, base_col, use_alpha) {
      if (any(is.na(scores))) scores[is.na(scores)] <- min(scores, na.rm = TRUE)
      if (use_alpha) {
        if (max(scores) == min(scores)) {
          alphas <- rep(0.8, length(scores))
        } else {
          alphas <- scales::rescale(scores, to = c(0.1, 1.0))
        }
      } else {
        alphas <- rep(0.8, length(scores))
      }
      rgb_val <- col2rgb(base_col)
      rgb(rgb_val[1], rgb_val[2], rgb_val[3], alpha = alphas * 255, maxColorValue = 255)
    }

    do_map <- has_score && score_to_alpha
    bez_df$final_color <- calc_alpha_color(bez_df$score, loop_color, do_map)
    anchors$final_fill <- calc_alpha_color(anchors$score, anchor_color, do_map)
  }

  # === 3. Overlap BED Data ===
  overlap_df_plot <- NULL
  if (!is.null(target_bed)) {
    ob <- read.table(target_bed, header = FALSE, sep = "\t")[c(1, 2, 3)]
    colnames(ob) <- c("chr", "start", "end")
    overlap_df_plot <- ob %>% dplyr::filter(chr == chr, !(end < from | start > to))
    if (nrow(overlap_df_plot) > 0) {
      overlap_df_plot <- overlap_df_plot %>% dplyr::mutate(ymin = -0.15, ymax = -0.10)
    }
  }

  # === 4. Gene Data ===
  genes_gr <- GenomicFeatures::genes(txdb)
  genes_df <- as.data.frame(genes_gr) %>% dplyr::filter(as.character(seqnames) == chr, end > from, start < to)

  feature_df <- data.frame()
  if (nrow(genes_df) > 0) {
    if (!requireNamespace(org_db_pkg, quietly = TRUE)) {
      stop("Package ", org_db_pkg, " is not installed.")
    }
    org_db <- utils::getFromNamespace(org_db_pkg, org_db_pkg)
    try(
      {
        symbol_map <- AnnotationDbi::select(org_db, keys = unique(as.character(genes_df$gene_id)), columns = "SYMBOL", keytype = "ENTREZID")
        symbol_map <- symbol_map[!duplicated(symbol_map$ENTREZID), ]
        genes_df <- dplyr::left_join(genes_df, symbol_map, by = c("gene_id" = "ENTREZID"))
      },
      silent = TRUE
    )

    if (!"SYMBOL" %in% colnames(genes_df)) genes_df$SYMBOL <- genes_df$gene_id
    genes_df <- genes_df %>% dplyr::mutate(final_label = ifelse(is.na(SYMBOL), gene_id, SYMBOL), label_x = pmax(from, pmin(to, ifelse(strand == "+", start, end))))

    bins <- IRanges::disjointBins(IRanges::IRanges(genes_df$start, genes_df$end))
    genes_df$gene_level <- bins

    tx_keytype <- if ("TXNAME" %in% AnnotationDbi::keytypes(txdb)) "TXNAME" else "TXID"
    tx2gene <- AnnotationDbi::select(txdb, keys = AnnotationDbi::keys(txdb, tx_keytype), columns = c("GENEID"), keytype = tx_keytype)
    colnames(tx2gene) <- c("tx_id", "gene_id")

    exons_list <- GenomicFeatures::exonsBy(txdb, "tx", use.names = TRUE)
    exons_gr <- unlist(exons_list)
    names(exons_gr) <- NULL
    exons_flat <- as.data.frame(exons_gr)
    exons_flat$tx_id <- rep(names(exons_list), times = S4Vectors::elementNROWS(exons_list))
    exons_joined <- exons_flat %>%
      dplyr::left_join(tx2gene, by = "tx_id") %>%
      dplyr::filter(gene_id %in% genes_df$gene_id) %>%
      dplyr::filter(start < to & end > from)

    if (nrow(exons_joined) > 0) {
      longest_tx <- exons_joined %>%
        dplyr::group_by(gene_id, tx_id) %>%
        dplyr::summarise(len = sum(width), .groups = "drop") %>%
        dplyr::arrange(desc(len)) %>%
        dplyr::group_by(gene_id) %>%
        dplyr::slice(1)
      feature_df <- exons_joined %>% dplyr::filter(tx_id %in% longest_tx$tx_id)
      feature_df <- dplyr::left_join(feature_df, genes_df[, c("gene_id", "gene_level")], by = "gene_id")
    }
  }

  # === 5. Gene Y Coordinates ===
  gene_start_y <- -0.25
  row_height <- 0.12
  plot_ymin <- gene_start_y

  if (nrow(genes_df) > 0) {
    genes_df <- genes_df %>% dplyr::mutate(y_mid = gene_start_y - (gene_level - 1) * row_height)
    plot_ymin <- min(genes_df$y_mid) - 0.2
    if (nrow(feature_df) > 0) {
      feature_df <- feature_df %>% dplyr::mutate(y_mid = gene_start_y - (gene_level - 1) * row_height, ymin = y_mid - 0.025, ymax = y_mid + 0.025)
    }
  }

  # === 6. Plotting ===
  p <- ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = -0.04, linetype = "dashed", color = "grey85", size = 0.5) +
    ggplot2::geom_hline(yintercept = -0.18, linetype = "dashed", color = "grey85", size = 0.5) +
    ggplot2::annotate("text", x = text_indent, y = plot_ymax - 0.01, label = "loop track", hjust = 0, vjust = 1, size = 4, fontface = "bold", color = "black")

  if (!is.null(overlap_df_plot)) {
    p <- p + ggplot2::annotate("text", x = text_indent, y = -0.08, label = "target track", hjust = 0, vjust = 0, size = 4, fontface = "bold", color = "black")
  }

  p <- p + ggplot2::annotate("text", x = text_indent, y = -0.21, label = "gene track", hjust = 0, vjust = 0, size = 4, fontface = "bold", color = "black")


  if (!is.null(overlap_df_plot) && nrow(overlap_df_plot) > 0) {
    p <- p + ggplot2::geom_rect(data = overlap_df_plot, ggplot2::aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax), fill = overlap_color, alpha = 1)
  }

  if (nrow(genes_df) > 0) {
    p <- p + ggplot2::geom_segment(data = genes_df, ggplot2::aes(x = pmax(start, from), xend = pmin(end, to), y = y_mid, yend = y_mid), color = intron_color, size = 0.5) +
      ggplot2::geom_segment(data = genes_df, ggplot2::aes(x = ifelse(strand == "+", pmin(end, to), pmax(start, from)), xend = ifelse(strand == "+", pmin(end, to), pmax(start, from)), y = y_mid, yend = y_mid), arrow = ggplot2::arrow(length = ggplot2::unit(0.15, "cm"), type = "open"), color = intron_color, size = 0.5)

    if (nrow(feature_df) > 0) {
      p <- p + ggplot2::geom_rect(data = feature_df, ggplot2::aes(xmin = pmax(start, from), xmax = pmin(end, to), ymin = ymin, ymax = ymax), fill = exon_color, color = NA)
    } else {
      p <- p + ggplot2::geom_rect(data = genes_df, ggplot2::aes(xmin = pmax(start, from), xmax = pmin(end, to), ymin = y_mid - 0.025, ymax = y_mid + 0.025), fill = exon_color)
    }

    p <- p + ggrepel::geom_text_repel(data = genes_df, ggplot2::aes(x = label_x, y = y_mid, label = final_label), nudge_y = -0.05, direction = "x", force = 1, size = 3, segment.size = 0.3, segment.color = "grey60", segment.linetype = "dashed", min.segment.length = 0)
  }

  if (nrow(bez_df) > 0) {
    p <- p +
      ggplot2::geom_rect(data = anchors, ggplot2::aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax, fill = final_fill), color = NA) +
      ggplot2::scale_fill_identity() +
      ggforce::geom_bezier(data = bez_df, ggplot2::aes(x = x, y = y, group = loop_i, color = final_color), size = 0.6) +
      ggplot2::scale_color_identity()
  }

  p <- p +
    ggplot2::coord_cartesian(xlim = c(from, to), ylim = c(plot_ymin, plot_ymax), expand = FALSE) +
    ggplot2::scale_x_continuous(labels = scales::comma) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.line.y = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(), axis.title.y = ggplot2::element_blank(), panel.border = ggplot2::element_rect(color = "black", fill = NA, size = 1), plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::labs(x = paste0(chr, ":", from, "-", to), title = "Loops Integrative View")

  if (!is.null(save_file)) ggplot2::ggsave(save_file, p, width = 12, height = 6)
  return(p)
}

#' Draw Simplified Flower Plot for Core vs. Unique Genes
#'
#' Creates a circular "flower" diagram where each petal represents a gene set,
#' showing the number of genes unique to that set. The center displays the size
#' of the core intersection (genes shared by all sets). Designed for intuitive
#' comparison of shared vs. condition-specific genes across 2–6 groups.
#'
#' @param gene_lists A named list of character vectors containing gene identifiers (e.g., symbols or Entrez IDs).
#' @param project_name Character. Prefix for the plot title (e.g., "Differential Loops").
#' @param filename Character. Output file path (e.g., "flower.png" or "flower.pdf").
#' @param group_colors Character vector. Colors for each group (named or in same order as \code{gene_lists}).
#' @return Invisibly returns \code{NULL}. Saves the plot to \code{filename}.
#' @importFrom ggplot2 ggplot geom_polygon annotate coord_fixed theme_void labs ggsave theme element_text margin
#' @importFrom ggforce geom_circle
#' @importFrom scales hue_pal
#' @export
#' @examples
#' # 1. Create dummy gene sets with some overlap
#' gene_sets <- list(
#'   Control = c("TP53", "BRCA1", "MYC"),
#'   Treated = c("BRCA1", "MYC", "EGFR"),
#'   Resistant = c("MYC", "EGFR", "KRAS")
#' )
#' # 2. Draw the flower plot
#' draw_flower_simplified(
#'   gene_lists = gene_sets,
#'   project_name = "Drug Response Study",
#'   filename = tempfile(fileext = ".png"),
#'   group_colors = c(Control = "#E41A1C", Treated = "#377EB8", Resistant = "#4DAF4A")
#' )
draw_flower_simplified <- function(gene_lists, project_name, filename, group_colors) {
  gene_lists <- gene_lists[vapply(gene_lists, length, FUN.VALUE = integer(1)) > 0]
  n_groups <- length(gene_lists)
  if (n_groups < 2) {
    message("Less than 2 non-empty gene lists; skipping flower plot.")
    return(invisible(NULL))
  }
  group_names <- names(gene_lists)

  final_colors <- if (!is.null(names(group_colors))) {
    group_colors[names(gene_lists)]
  } else {
    group_colors[seq_along(gene_lists)]
  }
  if (any(is.na(final_colors))) final_colors <- scales::hue_pal()(n_groups)

  core_genes <- Reduce(intersect, gene_lists)
  core_count <- length(core_genes)

  petal_counts <- vapply(group_names, function(g) {
    others_union <- unique(unlist(gene_lists[setdiff(group_names, g)]))
    unique_to_g <- setdiff(gene_lists[[g]], others_union)
    length(unique_to_g)
  }, FUN.VALUE = integer(1))

  center_x <- 0
  center_y <- 0
  ellipse_a <- 3.8
  ellipse_b <- 1.6
  r_offset <- 1.1
  core_radius <- 1.4

  get_ellipse <- function(angle_deg, cx, cy, a, b, offset, group_lbl) {
    t <- seq(0, 2 * pi, length.out = 100)
    rad <- angle_deg * pi / 180
    x <- a * cos(t)
    y <- b * sin(t)
    x_rot <- x * cos(rad) - y * sin(rad)
    y_rot <- x * sin(rad) + y * cos(rad)
    data.frame(x = x_rot + offset * cos(rad) + cx, y = y_rot + offset * sin(rad) + cy, group = group_lbl)
  }

  plot_df <- data.frame()
  group_label_df <- data.frame()
  count_label_df <- data.frame()
  angles <- seq(90, 90 + 360 * (n_groups - 1) / n_groups, length.out = n_groups)

  for (i in seq_len(n_groups)) {
    nm <- group_names[i]
    ang <- angles[i]
    coords <- get_ellipse(ang, center_x, center_y, ellipse_a, ellipse_b, r_offset, nm)
    plot_df <- rbind(plot_df, coords)

    group_lab_r <- r_offset + ellipse_a + 0.6
    group_label_df <- rbind(group_label_df, data.frame(x = center_x + group_lab_r * cos(ang * pi / 180), y = center_y + group_lab_r * sin(ang * pi / 180), label = nm, group = nm))

    count_lab_r <- r_offset + ellipse_a * 0.65
    count_label_df <- rbind(count_label_df, data.frame(x = center_x + count_lab_r * cos(ang * pi / 180), y = center_y + count_lab_r * sin(ang * pi / 180), label = as.character(petal_counts[i]), group = nm))
  }

  p <- ggplot() +
    geom_polygon(data = plot_df, aes(x = x, y = y, fill = group), alpha = 0.6, color = "white", size = 0.7) +
    ggforce::geom_circle(aes(x0 = 0, y0 = 0, r = core_radius), fill = "white", color = "grey70", size = 1) +
    annotate("text", x = 0, y = 0.35, label = "Core", color = "grey50", size = 5.5, fontface = "bold") +
    annotate("text", x = 0, y = -0.3, label = core_count, color = "black", size = 9, fontface = "bold") +
    geom_text(data = count_label_df, aes(x = x, y = y, label = label), color = "black", size = 7, fontface = "bold") +
    geom_text(data = group_label_df, aes(x = x, y = y, label = label, color = group), size = 6, fontface = "bold") +
    scale_fill_manual(values = final_colors) +
    scale_color_manual(values = final_colors) +
    coord_fixed() +
    theme_void() +
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 16, face = "bold", margin = margin(b = 20)), plot.margin = margin(30, 30, 30, 30)) +
    labs(title = paste0(project_name, ": Simplified Flower Plot\n(Core vs. Unique)"))

  ggsave(filename, p, width = 9, height = 9, bg = "white")
  message("    Saved (Simplified Flower Plot with inner counts): ", filename)
}


#' Generate UpSet Plot for Gene Set Intersections
#'
#' Visualizes intersections among multiple gene sets using the classic UpSetR package.
#' Uses grid graphics capture to ensure plot generation in all environments.
#'
#' @param gene_lists A named list of character vectors of gene identifiers.
#' @param project_name Character. Used for the file title.
#' @param filename Character. Output file path (must end in .pdf or .png).
#' @param group_colors Optional named character vector. Not used in this version.
#' @return Invisibly returns \code{NULL}. Saves the plot to \code{filename}.
#'
#' @importFrom UpSetR fromList upset
#' @importFrom grDevices pdf png dev.off
#' @importFrom grid grid.grabExpr grid.draw
#' @export
#' @examples
#' gene_sets <- list(
#'   Upregulated = c("A", "B", "C", "D"),
#'   Downregulated = c("C", "D", "E", "F"),
#'   Bound_by_TF = c("B", "D", "F", "G")
#' )
#'
#' tf <- tempfile(fileext = ".pdf")
#'
#' if (requireNamespace("UpSetR", quietly = TRUE)) {
#'   draw_upset_intersections(
#'     gene_lists = gene_sets,
#'     project_name = "Transcriptional Regulation",
#'     filename = tf
#'   )
#' }
draw_upset_intersections <- function(gene_lists, project_name, filename, group_colors = NULL) {
  # 1. Check Dependency
  if (!requireNamespace("UpSetR", quietly = TRUE)) {
    stop("Package 'UpSetR' required for this function. Please install it.")
  }

  # 2. Data Preparation
  gene_lists <- gene_lists[vapply(gene_lists, length, FUN.VALUE = integer(1)) > 0]
  if (length(gene_lists) < 2) {
    message("Less than 2 non-empty gene lists; skipping UpSet plot.")
    return(invisible(NULL))
  }

  # Convert list to UpSetR binary matrix format
  input_data <- UpSetR::fromList(gene_lists)

  # Check if data is valid
  if (nrow(input_data) == 0 || ncol(input_data) == 0) {
    warning("Input data for UpSet plot is empty. Skipping.")
    return(invisible(NULL))
  }

  # 3. Plotting (Capture Phase)
  # UpSetR draws as a side-effect. We use grid.grabExpr to capture
  # everything it draws into a single "Grob" (Graphical Object).
  # This prevents the "Blank PDF" issue caused by device desync.

  upset_grob <- tryCatch(
    {
      grid::grid.grabExpr({
        UpSetR::upset(
          input_data,
          nsets = length(gene_lists),
          nintersects = 40,
          mb.ratio = c(0.55, 0.45),
          order.by = "freq",
          mainbar.y.label = "Gene Intersection Size",
          sets.x.label = "Set Size",
          text.scale = c(1.3, 1.3, 1, 1, 1.3, 1)
        )
      })
    },
    error = function(e) {
      warning("UpSetR plotting failed: ", e$message)
      return(NULL)
    }
  )

  if (is.null(upset_grob)) {
    return(invisible(NULL))
  }

  # 4. Saving Phase (Draw Phase)
  # Now we simply draw the captured object to the file device.
  tryCatch(
    {
      is_pdf <- grepl("\\.pdf$", filename, ignore.case = TRUE)

      if (is_pdf) {
        grDevices::pdf(file = filename, width = 10, height = 6)
      } else {
        grDevices::png(file = filename, width = 10, height = 6, units = "in", res = 300)
      }

      # Explicitly draw the captured plot
      grid::grid.draw(upset_grob)

      grDevices::dev.off()
      message("    Saved (UpSet Plot): ", filename)
    },
    error = function(e) {
      try(grDevices::dev.off(), silent = TRUE)
      message("     Save failed: ", e$message)
    }
  )

  return(invisible(NULL))
}
