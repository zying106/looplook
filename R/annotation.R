#' @title Advanced 3D Chromatin Loop Annotation & Multi-Omics Feature Mapping Engine
#'
#' @description
#' A highly sophisticated, dual-purpose 3D genomic engine designed to:
#' \enumerate{
#'   \item \strong{Annotate Chromatin Loops:} Accurately classify 3D spatial interactions (e.g., Enhancer-Promoter, Promoter-Promoter) using a strict structural hierarchy.
#'   \item \strong{Map External Features (genomic features to 3D):} Act as a \strong{"3D Spatial Bridge"} to link auxiliary 1D genomic features (e.g., GWAS risk SNPs, ATAC-seq peaks, ChIP-seq binding sites) to their genuine regulated target genes, fundamentally outperforming simple "nearest gene" heuristics.
#' }
#'
#' \strong{Core Philosophy: Bridging the genomic features-3D Gap}
#' Traditional annotations typically assign non-coding variants to the nearest linear gene. This engine prioritizes \strong{physical 3D chromatin contacts}.
#' If a GWAS SNP or enhancer lands in an anchor looping to a distal gene, the distal gene is accurately assigned. If no spatial loop exists, the engine intelligently falls back to the nearest active local gene (Smart Fallback), ensuring comprehensive and gapless annotation coverage.
#'
#' \strong{Key Algorithmic Innovations:}
#' \itemize{
#'   \item \strong{Biotype & Expression-Aware Conflict Resolution:}
#'     When an anchor overlaps multiple promoters (e.g., dense gene loci or bidirectional promoters), it executes a rigorous 3-step resolution:
#'     \enumerate{
#'       \item \emph{Expression Pre-filter:} Eliminates transcriptionally silent genes based on the user-provided expression matrix.
#'       \item \emph{Functional Biotype Prioritization::} Prioritizes the remaining candidates by functional class in the following order: \code{Protein Coding > Antisense > lncRNA > Pseudogene}.
#'       \item \emph{Dominant Expression Tiebreaker:} Designates the gene with the highest transcriptional abundance as the target gene to further resolve any remaining mapping ambiguities.
#'     }
#'   \item \strong{Dynamic Topology Control (\code{neighbor_hop}):}
#'     Controls network diffusion depth. Allows signal propagation from a SNP -> Enhancer -> Hub -> Target Gene, which is ideal for uncovering multi-way super-enhancer cliques.
#'   \item \strong{Hub Detection:}
#'     Calculates precise node degrees to identify high-connectivity Promoter and Distal hubs based on the \code{hub_percentile}.
#' }
#'
#' @param bedpe_file Character. Path to the BEDPE file containing 3D loop coordinates (The "3D Bridge").
#' @param target_bed Optional character. Path to an auxiliary BED file (e.g., GWAS summary stats, eQTLs, ChIP-seq/ATAC-seq peaks). If provided, maps these 1D regulatory regions to 3D target genes.
#' @param species Character. Reference genome build. Supported: \code{"hg38"}, \code{"hg19"}, \code{"mm10"}, \code{"mm9"}.
#' @param tss_region Numeric vector. Defines the promoter window around the TSS. Default: \code{c(-2000, 2000)}.
#' @param out_dir Character. Directory for saving generated PDF plots and Excel results.
#' @param expr_matrix_file Optional character. Path to the normalized RNA-seq matrix (e.g., TPM/FPKM). \strong{Highly recommended} as it fuels the intelligent conflict resolution and filters out silent elements.
#' @param sample_columns Character vector or Integer indices. Specifies which columns in \code{expr_matrix_file} to average for the baseline expression reference.
#' @param project_name Character. Prefix for all output files and plot titles.
#' @param color_palette Character. Color brewer palette name for visualizations (e.g., \code{"Set2"}).
#' @param karyo_bin_size Numeric. Bin size for Karyotype genomic density heatmaps. Default: \code{1e5} (100kb).
#' @param neighbor_hop Integer. \strong{Network Diffusion Depth}. \code{0} (Default) captures direct physical connections only; \code{1} captures indirect 1-hop hub-mediated connections.
#' @param hub_percentile Numeric (0-1). Top percentile threshold for defining high-connectivity "Regulatory Hubs" (default: \code{0.95}).
#'
#' @return An invisible list of comprehensive data frames (also auto-saved as a multi-sheet \code{.xlsx}):
#' \itemize{
#'   \item \code{target_annotation}:
#'     Detailed annotation of \code{target_bed}. Features the important column: \code{Assigned_Target_Genes_Filled} (Loop-prioritized target, gracefully falling back to the local nearest gene if unlooped).
#'   \item \code{loop_annotation}:
#'     The fully annotated 3D network. Features \code{Putative_Target_Genes} capturing precisely resolved regulatory targets.
#'   \item \code{promoter_centric_stats}:
#'     Gene-level topological statistics. Contains structural node degrees and identifies \code{Is_High_Connectivity_Gene}.
#'   \item \code{distal_element_stats}:
#'     Enhancer-level topological statistics, highlighting critical regulatory anchors orchestrating multiple promoters.
#'   \item \code{anchor_annotation}:
#'     Locus-level summarization of all 1D anchor footprints.
#' }
#'
#' @export
#'
#' @examples
#' # 1. Get paths to example files included in the package
#' bedpe_path <- system.file("extdata", "example_loops_1.bedpe", package = "looplook")
#' bed_path <- system.file("extdata", "example_peaks.bed", package = "looplook")
#' expr_path <- system.file("extdata", "example_tpm.txt", package = "looplook")
#'
#' # 2. Check if files and required annotation databases exist
#' if (bedpe_path != "" && bed_path != "" && expr_path != "" &&
#'   requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE) &&
#'   requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
#'   # =========================================================================
#'   # Example A: Integrative Analysis (Loops + GWAS/ChIP-seq + Expression)
#'   # =========================================================================
#'   # Here, `target_bed` can represent GWAS risk SNPs, ATAC-seq peaks, or
#'   # ChIP-seq peaks. The engine bridges these 1D regions to their 3D genes.
#'   res_integrated <- annotate_peaks_and_loops(
#'     bedpe_file = bedpe_path,
#'     target_bed = bed_path,
#'     expr_matrix_file = expr_path,
#'     sample_columns = c("con1", "con2"),
#'     species = "hg38",
#'     tss_region = c(-2000, 2000),
#'     out_dir = tempdir(),
#'     color_palette = "Set2",
#'     karyo_bin_size = 1e5,
#'     neighbor_hop = 0,
#'     hub_percentile = 0.95,
#'     project_name = "Example_HiChIP_Integrative"
#'   )
#'
#'   # View integrated annotations linking external peaks/SNPs to target genes
#'   head(res_integrated$target_annotation)
#'
#'   # =========================================================================
#'   # Example B: Deep Analysis of Loops ONLY (No auxiliary peaks)
#'   # =========================================================================
#'   res_loops_only <- annotate_peaks_and_loops(
#'     bedpe_file = bedpe_path,
#'     target_bed = NULL,
#'     species = "hg38",
#'     expr_matrix_file = expr_path,
#'     sample_columns = c("con1", "con2"),
#'     tss_region = c(-2000, 2000),
#'     out_dir = tempdir(),
#'     color_palette = "Set1",
#'     karyo_bin_size = 1e5,
#'     neighbor_hop = 0,
#'     hub_percentile = 0.95,
#'     project_name = "Example_Loops_Only"
#'   )
#'
#'   # View standalone loop annotations
#'   head(res_loops_only$loop_annotation)
#' }
annotate_peaks_and_loops <- function(
  bedpe_file,
  target_bed = NULL,
  species = "hg38",
  tss_region = c(-2000, 2000),
  out_dir = "./results",
  expr_matrix_file = NULL,
  sample_columns = NULL,
  project_name = "HiChIP",
  color_palette = "Set2",
  karyo_bin_size = 1e5,
  neighbor_hop = 0,
  hub_percentile = 0.95
) {
  # Ensure output directory exists
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  # =========================================================
  # 0. Load database & define helper functions
  extract_genes <- function(genes_vec) {
    paste(unique(na.omit(unlist(strsplit(as.character(genes_vec), ";")))), collapse = ";")
  }
  extract_ids <- function(id_vec) {
    paste(unique(na.omit(as.character(id_vec))), collapse = ";")
  }

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

  if (!requireNamespace(txdb_pkg, quietly = TRUE)) stop("Package", txdb_pkg, "missing")
  if (!requireNamespace(org_db_pkg, quietly = TRUE)) stop("Package", org_db_pkg, "missing")
  txdb_obj <- utils::getFromNamespace(txdb_pkg, txdb_pkg)

  gene_expr_map <- NULL
  if (!is.null(expr_matrix_file) && !is.null(sample_columns)) {
    message("Step 0: Loading expression data...")
    d <- tryCatch(read.table(expr_matrix_file, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE), error = function(e) NULL)
    if (is.null(d)) d <- tryCatch(read.table(expr_matrix_file, header = TRUE, sep = ",", row.names = 1, check.names = FALSE), error = function(e) NULL)
    if (!is.null(d)) {
      if (is.character(sample_columns)) sub_mat <- d[, intersect(sample_columns, colnames(d)), drop = FALSE] else sub_mat <- d[, sample_columns, drop = FALSE]
      if (ncol(sub_mat) > 0) {
        vals <- if (ncol(sub_mat) > 1) rowMeans(sub_mat, na.rm = TRUE) else sub_mat[, 1]
        gene_expr_map <- vals
        message("    >>> Expression loaded for ", length(vals), " genes.")
      }
    }
  }


  resolve_gene_conflicts <- function(gr_input, current_anno_df) {
    if (length(gr_input) == 0) {
      return(current_anno_df)
    }
    tryCatch(
      {
        all_genes <- GenomicFeatures::genes(txdb_obj)
        hits <- GenomicRanges::findOverlaps(gr_input, GenomicFeatures::promoters(all_genes, upstream = abs(tss_region[1]), downstream = abs(tss_region[2])))

        if (length(hits) > 0) {
          candidates <- data.frame(query_idx = S4Vectors::queryHits(hits), gene_id = names(all_genes)[S4Vectors::subjectHits(hits)], stringsAsFactors = FALSE)


          gene_map <- suppressMessages(AnnotationDbi::select(utils::getFromNamespace(org_db_pkg, org_db_pkg), keys = unique(candidates$gene_id), columns = c("SYMBOL", "GENETYPE"), keytype = "ENTREZID"))


          gene_map$tpm <- if (!is.null(gene_expr_map)) ifelse(is.na(gene_expr_map[gene_map$SYMBOL]), 0, gene_expr_map[gene_map$SYMBOL]) else 0


          gene_map <- gene_map %>%
            dplyr::mutate(
              type_rank = dplyr::case_when(
                grepl("protein", GENETYPE, ignore.case = TRUE) ~ 1,
                grepl("antisense", GENETYPE, ignore.case = TRUE) ~ 2,
                grepl("lncRNA|ncrna", GENETYPE, ignore.case = TRUE) ~ 3,
                grepl("pseudo", GENETYPE, ignore.case = TRUE) ~ 4,
                TRUE ~ 5
              )
            )

          resolved_candidates <- candidates %>%
            dplyr::left_join(gene_map, by = c("gene_id" = "ENTREZID")) %>%
            dplyr::group_by(query_idx) %>%
            dplyr::mutate(has_active = any(tpm > 0)) %>%
            dplyr::filter(!has_active | tpm > 0) %>%
            dplyr::filter(type_rank == min(type_rank, na.rm = TRUE)) %>%
            dplyr::summarise(
              valid_genes = list(SYMBOL[!is.na(SYMBOL) & SYMBOL != ""]),
              valid_tpms = list(tpm[!is.na(SYMBOL) & SYMBOL != ""]),
              .groups = "drop"
            ) %>%
            dplyr::rowwise() %>%
            dplyr::mutate(
              final_symbols = {
                genes <- unlist(valid_genes)
                tpms <- unlist(valid_tpms)

                if (length(genes) == 0) {
                  NA_character_
                } else if (length(genes) == 1) {
                  genes[1]
                } else {
                  max_tpm <- max(tpms, na.rm = TRUE)
                  if (max_tpm <= 0) {
                    paste(sort(unique(genes)), collapse = ";")
                  } else {
                    threshold <- max_tpm * 0.1
                    active_genes <- genes[tpms >= threshold]
                    paste(sort(unique(active_genes)), collapse = ";")
                  }
                }
              }
            ) %>%
            dplyr::ungroup() %>%
            dplyr::filter(!is.na(final_symbols))


          match_idx <- match(resolved_candidates$query_idx, seq_len(nrow(current_anno_df)))
          current_anno_df$SYMBOL[match_idx] <- resolved_candidates$final_symbols
          current_anno_df$annotation[match_idx] <- "Promoter"
        }
        return(current_anno_df)
      },
      error = function(e) {
        message("    [Conflict Resolution] ", e$message)
        return(current_anno_df)
      }
    )
  }


  get_feature_class <- function(anno_str) {
    if (is.na(anno_str)) {
      return("Unknown")
    }
    anno_str <- tolower(anno_str)

    if (grepl("promoter", anno_str)) {
      return("P")
    }

    if (grepl("intergenic|downstream", anno_str)) {
      return("E")
    }
    if (grepl("exon|intron|utr", anno_str)) {
      return("G")
    }

    return("E")
  }

  message("Step 1: Reading BEDPE file...")
  loops <- read_robust_general(bedpe_file, min_cols = 6, desc = "BEDPE")
  colnames(loops)[seq_len(6)] <- c("chr1", "start1", "end1", "chr2", "start2", "end2")
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
  gr_anchors <- GenomicRanges::makeGRangesFromDataFrame(anchors, keep.extra.columns = TRUE)
  gr_anchors$anchor_id <- anchors$anchor_id
  gr_list <- GenomicRanges::GRangesList(split(gr_anchors, gr_anchors$cluster_id))
  cluster_regions <- unlist(GenomicRanges::reduce(gr_list))
  cluster_regions$cluster_id <- names(cluster_regions)
  names(cluster_regions) <- paste0("peak_", seq_along(cluster_regions))

  message("Step 3: Biological Classification & Topology...")
  anchor_anno <- ChIPseeker::annotatePeak(gr_anchors, TxDb = txdb_obj, tssRegion = tss_region, annoDb = org_db_pkg, verbose = FALSE)
  anchor_anno_df <- format_annotation_columns(as.data.frame(anchor_anno))
  anchor_anno_df <- resolve_gene_conflicts(gr_anchors, anchor_anno_df)
  anchor_anno_df$type_code <- vapply(anchor_anno_df$annotation, get_feature_class, FUN.VALUE = character(1))
  map_info <- anchor_anno_df %>% dplyr::select(anchor_id, type_code, SYMBOL)
  loops_annotated <- loops %>%
    dplyr::left_join(map_info %>% dplyr::rename(t1 = type_code, s1 = SYMBOL), by = c("a1_id" = "anchor_id")) %>%
    dplyr::left_join(map_info %>% dplyr::rename(t2 = type_code, s2 = SYMBOL), by = c("a2_id" = "anchor_id"))

  get_locus_genes <- function(t1, t2, s1, s2) {
    genes <- c()
    if (!is.na(t1) && t1 %in% c("P", "G")) genes <- c(genes, s1)
    if (!is.na(t2) && t2 %in% c("P", "G")) genes <- c(genes, s2)
    paste(unique(na.omit(genes)), collapse = ";")
  }
  get_type_code <- function(t1, t2) {
    if (is.na(t1) || is.na(t2)) {
      return("Unknown")
    }
    paste(sort(c(t1, t2)), collapse = "-")
  }
  loops_annotated$loop_type <- apply(loops_annotated, 1, function(r) get_type_code(r["t1"], r["t2"]))
  loops_annotated$single_loop_genes <- apply(loops_annotated, 1, function(r) get_locus_genes(r["t1"], r["t2"], r["s1"], r["s2"]))
  loops_annotated$reg_loop_genes <- apply(loops_annotated, 1, function(r) get_locus_genes(r["t1"], r["t2"], r["s1"], r["s2"]))

  message("    Calculating Topology (Hops)...")
  map_info$SYMBOL <- trimws(map_info$SYMBOL)
  valid_pg_nodes <- map_info %>% dplyr::filter(type_code %in% c("P", "G") & !is.na(SYMBOL) & SYMBOL != "")
  lookup_pg_symbol <- valid_pg_nodes$SYMBOL
  names(lookup_pg_symbol) <- valid_pg_nodes$anchor_id
  lookup_pg_type <- valid_pg_nodes$type_code
  names(lookup_pg_type) <- valid_pg_nodes$anchor_id
  lookup_p_symbol <- map_info %>%
    dplyr::filter(type_code == "P" & !is.na(SYMBOL) & SYMBOL != "") %>%
    dplyr::pull(SYMBOL)
  names(lookup_p_symbol) <- map_info %>%
    dplyr::filter(type_code == "P" & !is.na(SYMBOL) & SYMBOL != "") %>%
    dplyr::pull(anchor_id)
  nodes_in_graph <- igraph::V(g)$name
  ids_to_genes_simple <- function(ids, lookup) {
    valid <- intersect(ids, names(lookup))
    if (length(valid) == 0) {
      return(NA_character_)
    }
    paste(sort(unique(lookup[valid])), collapse = ";")
  }

  ids_to_genes_priority <- function(ids, lookup_sym, lookup_typ) {
    valid <- intersect(ids, names(lookup_sym))
    if (length(valid) == 0) {
      return(NA_character_)
    }
    genes_present <- lookup_sym[valid]
    paste(sort(unique(genes_present)), collapse = ";")
  }

  input_hop <- if (is.null(neighbor_hop)) 0 else neighbor_hop
  ego_list_loop <- igraph::ego(g, order = input_hop, nodes = nodes_in_graph, mode = "all")
  names(ego_list_loop) <- nodes_in_graph
  ego_list_target <- igraph::ego(g, order = input_hop + 1, nodes = nodes_in_graph, mode = "all")
  names(ego_list_target) <- nodes_in_graph
  anchor_topo_map <- data.frame(anchor_id = nodes_in_graph, topo_genes_p = vapply(ego_list_loop, function(x) ids_to_genes_simple(names(x), lookup_p_symbol), character(1)), topo_genes_pg = vapply(ego_list_loop, function(x) ids_to_genes_simple(names(x), lookup_pg_symbol), character(1)), tgt_genes_pg = vapply(ego_list_target, function(x) ids_to_genes_simple(names(x), lookup_pg_symbol), character(1)), tgt_genes_p = vapply(ego_list_target, function(x) ids_to_genes_simple(names(x), lookup_p_symbol), character(1)), tgt_genes_prio = vapply(ego_list_target, function(x) ids_to_genes_priority(names(x), lookup_pg_symbol, lookup_pg_type), character(1)), stringsAsFactors = FALSE)
  anchor_topo_map[is.na(anchor_topo_map)] <- NA_character_

  message("Step 4: Constructing Loop Tables...")
  loops_annotated <- loops_annotated %>%
    dplyr::left_join(anchor_topo_map %>% dplyr::select(anchor_id, pg1 = topo_genes_pg, p1 = topo_genes_p), by = c("a1_id" = "anchor_id")) %>%
    dplyr::left_join(anchor_topo_map %>% dplyr::select(anchor_id, pg2 = topo_genes_pg, p2 = topo_genes_p), by = c("a2_id" = "anchor_id")) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(proximate_loop_gene = dplyr::case_when((!is.na(t1) & t1 == "G" & !is.na(t2) & t2 == "P") ~ extract_genes(pg2), (!is.na(t1) & t1 == "P" & !is.na(t2) & t2 == "G") ~ extract_genes(pg1), TRUE ~ extract_genes(c(pg1, pg2)))) %>%
    dplyr::ungroup()
  clust_vec <- setNames(anchors$cluster_id, anchors$anchor_id)
  loops_annotated$cluster_id <- clust_vec[loops_annotated$a1_id]
  agg_cluster_reg <- loops_annotated %>%
    dplyr::filter(!is.na(cluster_id)) %>%
    dplyr::group_by(cluster_id) %>%
    dplyr::summarise(all_cluster_loop_genes = extract_genes(reg_loop_genes), .groups = "drop")
  loop_annotation_final <- loops_annotated %>%
    dplyr::left_join(agg_cluster_reg, by = "cluster_id") %>%
    dplyr::mutate(loop_ID = paste0("L", seq_len(dplyr::n()))) %>%
    dplyr::select(loop_ID, chr1, start1, end1, chr2, start2, end2, cluster_id, loop_type, anchor1_gene = s1, anchor1_type = t1, anchor2_gene = s2, anchor2_type = t2, all_cluster_loop_genes, single_loop_genes, proximate_loop_gene, a1_id, a2_id) %>%
    dplyr::rename(All_Anchor_Genes = single_loop_genes, Putative_Target_Genes = proximate_loop_gene, Cluster_All_Genes = all_cluster_loop_genes)
  agg_cluster_locus <- loop_annotation_final %>%
    dplyr::filter(!is.na(cluster_id)) %>%
    dplyr::group_by(cluster_id) %>%
    dplyr::summarise(Cluster_Locus_Genes = extract_genes(All_Anchor_Genes), .groups = "drop")
  gene_annot <- ChIPseeker::annotatePeak(cluster_regions, TxDb = txdb_obj, tssRegion = tss_region, annoDb = org_db_pkg, verbose = FALSE)
  cluster_info <- format_annotation_columns(as.data.frame(gene_annot))
  if ("GENENAME" %in% colnames(cluster_info)) cluster_info <- cluster_info %>% dplyr::rename(Gene_description = GENENAME)
  cluster_info$cluster_id <- as.character(cluster_info$cluster_id)
  cluster_info <- cluster_info %>% dplyr::left_join(agg_cluster_locus, by = "cluster_id")

  message("    Generating Promoter Centric Stats...")
  raw_stats_df <- dplyr::bind_rows(
    loop_annotation_final %>% dplyr::filter(anchor1_type == "P" & !is.na(anchor1_gene)) %>% dplyr::select(Gene = anchor1_gene, Neighbor_Type = anchor2_type, Loop_Type = loop_type),
    loop_annotation_final %>% dplyr::filter(anchor2_type == "P" & !is.na(anchor2_gene)) %>% dplyr::select(Gene = anchor2_gene, Neighbor_Type = anchor1_type, Loop_Type = loop_type)
  ) %>%
    tidyr::separate_rows(Gene, sep = ";") %>%
    dplyr::mutate(Gene = trimws(Gene)) %>%
    dplyr::filter(Gene != "" & !is.na(Gene)) %>%
    dplyr::group_by(Gene) %>%
    dplyr::summarise(
      Total_Loops = dplyr::n(),
      n_Linked_Promoters = sum(Neighbor_Type == "P", na.rm = TRUE),
      n_Linked_Distal = sum(Neighbor_Type %in% c("E", "G"), na.rm = TRUE),
      Dominant_Interaction = names(which.max(table(Loop_Type))),
      .groups = "drop"
    )

  final_cutoff <- max(quantile(raw_stats_df$Total_Loops, hub_percentile, na.rm = TRUE), 3)
  distal_cutoff <- max(quantile(raw_stats_df$n_Linked_Distal, hub_percentile, na.rm = TRUE), 2)

  promoter_centric_df <- raw_stats_df %>%
    dplyr::mutate(
      Is_High_Connectivity_Gene = dplyr::if_else(Total_Loops >= final_cutoff, "Yes", "No"),
      Is_High_Distal_Connectivity_Gene = dplyr::if_else(n_Linked_Distal >= distal_cutoff, "Yes", "No")
    ) %>%
    dplyr::arrange(dplyr::desc(n_Linked_Distal))

  message("    Generating Distal Element Stats...")
  distal_raw_df <- dplyr::bind_rows(
    loop_annotation_final %>% dplyr::filter(anchor1_type %in% c("E", "G")) %>% dplyr::select(Distal_Anchor_ID = a1_id, Distal_Type = anchor1_type, Neighbor_Gene = anchor2_gene, Neighbor_Type = anchor2_type, Loop_Type = loop_type),
    loop_annotation_final %>% dplyr::filter(anchor2_type %in% c("E", "G")) %>% dplyr::select(Distal_Anchor_ID = a2_id, Distal_Type = anchor2_type, Neighbor_Gene = anchor1_gene, Neighbor_Type = anchor1_type, Loop_Type = loop_type)
  ) %>%
    dplyr::group_by(Distal_Anchor_ID) %>%
    dplyr::summarise(
      Total_Loops = dplyr::n(),
      n_Linked_Promoters = sum(Neighbor_Type == "P", na.rm = TRUE),
      n_Linked_Distal = sum(Neighbor_Type %in% c("E", "G"), na.rm = TRUE),
      Dominant_Interaction = names(which.max(table(Loop_Type))),
      Target_Genes = extract_genes(Neighbor_Gene[Neighbor_Type == "P"]),
      .groups = "drop"
    )

  anchor_coords_map <- anchors %>%
    dplyr::select(anchor_id, chr, start, end, cluster_id) %>%
    dplyr::distinct()

  if (nrow(distal_raw_df) > 0) {
    final_cutoff_dist <- max(quantile(distal_raw_df$Total_Loops, hub_percentile, na.rm = TRUE), 3)
    distal_element_df <- distal_raw_df %>%
      dplyr::left_join(anchor_coords_map, by = c("Distal_Anchor_ID" = "anchor_id")) %>%
      dplyr::mutate(Is_High_Connectivity_Distal_Element = dplyr::if_else(Total_Loops >= final_cutoff_dist, "Yes", "No")) %>%
      dplyr::select(chr, start, end, cluster_id, Total_Loops, n_Linked_Promoters, n_Linked_Distal, Dominant_Interaction, Is_High_Connectivity_Distal_Element, Target_Genes) %>%
      dplyr::arrange(dplyr::desc(n_Linked_Promoters))
  } else {
    distal_element_df <- NULL
  }

  bed_info <- NULL
  target_connected_loops <- NULL
  if (!is.null(target_bed)) {
    message("Step 5: Integrating Target Annotations...")
    bed_target <- read_robust_general(target_bed, min_cols = 3, desc = "Target BED")
    colnames(bed_target)[c(1, 2, 3)] <- c("chr", "start", "end")
    gr_bed <- GenomicRanges::makeGRangesFromDataFrame(bed_target)
    gr_bed$input_id <- paste0("Peak_", seq_len(nrow(bed_target)))
    names(gr_bed) <- gr_bed$input_id
    bed_annot <- ChIPseeker::annotatePeak(gr_bed, TxDb = txdb_obj, tssRegion = tss_region, annoDb = org_db_pkg, verbose = FALSE)
    bed_info <- format_annotation_columns(as.data.frame(bed_annot))
    if ("GENENAME" %in% colnames(bed_info)) bed_info <- bed_info %>% dplyr::rename(Gene_description = GENENAME)
    message("    Refining Target annotation...")
    bed_info <- resolve_gene_conflicts(gr_bed, bed_info)
    hits <- GenomicRanges::findOverlaps(gr_bed, gr_anchors)
    if (length(hits) > 0) {
      target_connected_loops <- loop_annotation_final %>% dplyr::filter(cluster_id %in% unique(gr_anchors$cluster_id[S4Vectors::subjectHits(hits)]))
      hit_df <- data.frame(qid = S4Vectors::queryHits(hits), sid = S4Vectors::subjectHits(hits))
      hit_df$anchor_id <- gr_anchors$anchor_id[hit_df$sid]
      hit_df <- hit_df %>% dplyr::left_join(anchor_topo_map, by = "anchor_id")
      anchor_loop_agg <- dplyr::bind_rows(loop_annotation_final %>% dplyr::select(anchor_id = a1_id, loop_ID), loop_annotation_final %>% dplyr::select(anchor_id = a2_id, loop_ID)) %>%
        dplyr::distinct() %>%
        dplyr::group_by(anchor_id) %>%
        dplyr::summarise(linked_loops = extract_ids(loop_ID), .groups = "drop")
      hit_df <- hit_df %>% dplyr::left_join(anchor_loop_agg, by = "anchor_id")
      summary_df <- hit_df %>%
        dplyr::group_by(qid) %>%
        dplyr::summarise(All_Loop_Connected_Genes = extract_genes(tgt_genes_pg), Regulated_promoter_genes = extract_genes(tgt_genes_p), Assigned_Target_Genes = extract_genes(tgt_genes_prio), Linked_Loop_IDs = extract_ids(linked_loops), .groups = "drop") %>%
        dplyr::mutate(join_id = paste0("Peak_", qid))
      bed_info <- dplyr::left_join(bed_info, summary_df, by = c("input_id" = "join_id")) %>% dplyr::select(-any_of(c("join_id", "qid")))
    } else {
      bed_info$All_Loop_Connected_Genes <- NA
      bed_info$Regulated_promoter_genes <- NA
      bed_info$Assigned_Target_Genes <- NA
      bed_info$Linked_Loop_IDs <- NA
    }

    fill_logic <- function(target, fallback) dplyr::case_when(!is.na(target) & target != "" ~ target, !is.na(fallback) & fallback != "" ~ fallback, TRUE ~ NA_character_)

    bed_info <- bed_info %>% dplyr::mutate(
      All_Loop_Connected_Genes_Filled = fill_logic(All_Loop_Connected_Genes, SYMBOL),
      Regulated_promoter_genes_Filled = fill_logic(Regulated_promoter_genes, SYMBOL),
      Assigned_Target_Genes_Filled = fill_logic(Assigned_Target_Genes, SYMBOL)
    )

    if ("Linked_Loop_IDs" %in% colnames(bed_info)) {
      target_col <- if ("Gene_description" %in% colnames(bed_info)) "Gene_description" else "SYMBOL"
      if (target_col %in% colnames(bed_info)) bed_info <- bed_info %>% dplyr::relocate(Linked_Loop_IDs, .after = all_of(target_col))
    }
  }


  # Step 6: Visualization
  # =========================================================
  message("Step 6: Generating Visualizations...")
  plot_df <- loop_annotation_final
  plot_df$loop_genes <- plot_df$All_Anchor_Genes

  loop_types_sorted <- sort(unique(plot_df$loop_type))
  if (!exists("get_colors")) stop("Helper 'get_colors' missing.")
  custom_colors <- get_colors(length(loop_types_sorted), color_palette)
  names(custom_colors) <- loop_types_sorted

  red_palette <- c("#FFFFFF", "#FFFFCC", "#FFEDA0", "#FED976", "#FEB24C", "#FD8D3C", "#FC4E2A", "#E31A1C", "#BD0026", "#800026", "#000000")
  blue_palette <- c("#FFFFFF", "#E1F5FE", "#B3E5FC", "#4FC3F7", "#039BE5", "#0277BD", "#01579B", "#000000")
  purple_palette <- c("#FFFFFF", "#F3E5F5", "#E1BEE7", "#BA68C8", "#9C27B0", "#7B1FA2", "#4A148C", "#000000")

  # 1. Basic Charts
  if (requireNamespace("ggplot2", quietly = TRUE) && requireNamespace("dplyr", quietly = TRUE)) {
    tryCatch(
      {
        donut_data <- plot_df %>%
          dplyr::group_by(loop_type) %>%
          dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
          dplyr::mutate(
            prop = n / sum(n),
            legend_label = paste0(loop_type, " (n=", n, ", ", round(prop * 100, 1), "%)")
          ) %>%
          dplyr::arrange(dplyr::desc(n))

        donut_data$loop_type <- factor(donut_data$loop_type, levels = donut_data$loop_type)

        p_donut <- ggplot2::ggplot(donut_data, ggplot2::aes(x = 2, y = n, fill = loop_type)) +
          ggplot2::geom_bar(stat = "identity", color = "white") +
          ggplot2::coord_polar(theta = "y") +
          ggplot2::xlim(0.5, 2.9) +
          ggplot2::geom_text(
            ggplot2::aes(x = 2.65, label = loop_type),
            position = ggplot2::position_stack(vjust = 0.5),
            size = 2.5
          ) +
          ggplot2::scale_fill_manual(values = custom_colors, labels = setNames(donut_data$legend_label, donut_data$loop_type)) +
          ggplot2::theme_void() +
          ggplot2::labs(title = paste0(project_name, ": Loop Type Distribution")) +
          ggplot2::theme(
            legend.position = "right",
            plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
          )

        ggplot2::ggsave(file.path(out_dir, paste0(project_name, "_Basic_Donut.pdf")), p_donut, width = 8, height = 6)
      },
      error = function(e) warning("Basic Donut failed: ", e$message)
    )
  }

  if (exists("draw_circular_bar_plot")) {
    try(draw_circular_bar_plot(plot_df, project_name, file.path(out_dir, paste0(project_name, "_Basic_Circular.pdf")), custom_colors), silent = TRUE)
  }

  if ("All_Anchor_Genes" %in% colnames(plot_df) && exists("draw_karyo_heatmap_internal")) {
    genes_loop <- unique(trimws(unlist(strsplit(as.character(plot_df$All_Anchor_Genes), ";"))))
    genes_loop <- genes_loop[genes_loop != "" & !is.na(genes_loop)]
    if (length(genes_loop) > 0) {
      all_genes_gr <- GenomicFeatures::genes(txdb_obj)
      map <- AnnotationDbi::select(utils::getFromNamespace(org_db_pkg, org_db_pkg), keys = as.character(S4Vectors::mcols(all_genes_gr)$gene_id), columns = "SYMBOL", keytype = "ENTREZID")
      S4Vectors::mcols(all_genes_gr)$SYMBOL <- map$SYMBOL[match(S4Vectors::mcols(all_genes_gr)$gene_id, map$ENTREZID)]
      target_genes_gr <- all_genes_gr[S4Vectors::mcols(all_genes_gr)$SYMBOL %in% genes_loop]
      try(draw_karyo_heatmap_internal(target_genes_gr, "Loop Genes Distribution", file.path(out_dir, paste0(project_name, "_Basic_Karyo_LoopGenes.pdf")), karyo_bin_size, 0.99, txdb_obj, species, "Genes", custom_colors = red_palette), silent = TRUE)
    }
  }

  all_anchors <- dplyr::bind_rows(plot_df %>% dplyr::select(chr = chr1, start = start1, end = end1), plot_df %>% dplyr::select(chr = chr2, start = start2, end = end2)) %>% dplyr::distinct()
  if (nrow(all_anchors) > 0 && exists("draw_karyo_heatmap_internal")) {
    try(draw_karyo_heatmap_internal(GenomicRanges::makeGRangesFromDataFrame(all_anchors), "Loop Anchor Load", file.path(out_dir, paste0(project_name, "_Basic_Karyo_Anchors.pdf")), karyo_bin_size, 0.99, txdb_obj, species, "Anchors", custom_colors = blue_palette), silent = TRUE)
  }

  if (exists("draw_flower_simplified")) {
    temp_df_flower <- plot_df %>%
      dplyr::filter(!is.na(All_Anchor_Genes) & All_Anchor_Genes != "") %>%
      tidyr::separate_rows(All_Anchor_Genes, sep = ";") %>%
      dplyr::mutate(All_Anchor_Genes = trimws(All_Anchor_Genes)) %>%
      dplyr::filter(All_Anchor_Genes != "")
    gene_sets <- split(temp_df_flower$All_Anchor_Genes, temp_df_flower$loop_type)
    gene_sets <- lapply(gene_sets, unique)
    if (length(gene_sets) > 1) try(draw_flower_simplified(gene_sets, project_name, file.path(out_dir, paste0(project_name, "_Basic_Flower.pdf")), custom_colors), silent = TRUE)
  }

  if (!is.null(bed_info)) {
    message("    Plotting Target Visualizations...")
    bed_info_for_plot <- bed_info

    if ("Assigned_Target_Genes_Filled" %in% colnames(bed_info) && exists("draw_karyo_heatmap_internal")) {
      genes_target <- unique(trimws(unlist(strsplit(as.character(bed_info$Assigned_Target_Genes_Filled), ";"))))
      genes_target <- genes_target[genes_target != "" & !is.na(genes_target)]
      if (length(genes_target) > 0) {
        all_genes_gr <- GenomicFeatures::genes(txdb_obj)
        map <- AnnotationDbi::select(utils::getFromNamespace(org_db_pkg, org_db_pkg), keys = as.character(S4Vectors::mcols(all_genes_gr)$gene_id), columns = "SYMBOL", keytype = "ENTREZID")
        S4Vectors::mcols(all_genes_gr)$SYMBOL <- map$SYMBOL[match(S4Vectors::mcols(all_genes_gr)$gene_id, map$ENTREZID)]
        target_genes_gr <- all_genes_gr[S4Vectors::mcols(all_genes_gr)$SYMBOL %in% genes_target]
        try(draw_karyo_heatmap_internal(target_genes_gr, "Target Genes (Assigned+Local)", file.path(out_dir, paste0(project_name, "_Basic_Karyo_TargetGenes.pdf")), karyo_bin_size, 0.99, txdb_obj, species, "Genes", custom_colors = purple_palette), silent = TRUE)
      }
    }

    if (!is.null(target_connected_loops) && nrow(target_connected_loops) > 0 && requireNamespace("ggplot2", quietly = TRUE) && requireNamespace("dplyr", quietly = TRUE)) {
      tryCatch(
        {
          target_rose_data <- target_connected_loops %>%
            dplyr::group_by(loop_type) %>%
            dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
            dplyr::mutate(
              prop = n / sum(n),
              legend_label = paste0(loop_type, " (n=", n, ", ", round(prop * 100, 1), "%)")
            ) %>%
            dplyr::arrange(dplyr::desc(n))

          target_rose_data$loop_type <- factor(target_rose_data$loop_type, levels = target_rose_data$loop_type)

          p_target_rose <- ggplot2::ggplot(target_rose_data, ggplot2::aes(x = loop_type, y = n, fill = loop_type)) +
            ggplot2::geom_bar(stat = "identity", width = 1, color = "white") +
            ggplot2::coord_polar(theta = "x") +
            ggplot2::scale_fill_manual(values = custom_colors, labels = setNames(target_rose_data$legend_label, target_rose_data$loop_type)) +
            ggplot2::theme_void() +
            ggplot2::theme(
              legend.position = "right",
              plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
            ) +
            ggplot2::labs(title = paste0(project_name, ": Target Connected Loops (Rose)"))

          ggplot2::ggsave(file.path(out_dir, paste0(project_name, "_Target_Rose.pdf")), p_target_rose, width = 8, height = 6)
        },
        error = function(e) warning("Target Rose failed: ", e$message)
      )
    }


    draw_pie_with_outside_labels <- function(data_df, group_col, title, output_path, palette) {
      if (!requireNamespace("ggplot2", quietly = TRUE)) {
        return()
      }
      simplify_anno <- function(x) {
        if (grepl("Promoter", x, ignore.case = TRUE)) {
          return("Promoter")
        }
        if (grepl("Intron", x, ignore.case = TRUE)) {
          return("Intron")
        }
        if (grepl("Exon", x, ignore.case = TRUE)) {
          return("Exon")
        }
        if (grepl("Intergenic", x, ignore.case = TRUE)) {
          return("Distal Intergenic")
        }
        if (grepl("Downstream", x, ignore.case = TRUE)) {
          return("Downstream")
        }
        return("Others")
      }
      plot_data <- data_df
      plot_data$Simplified <- vapply(plot_data[[group_col]], simplify_anno, FUN.VALUE = character(1))

      stats <- plot_data %>%
        dplyr::group_by(Simplified) %>%
        dplyr::summarise(Count = dplyr::n(), .groups = "drop") %>%
        dplyr::mutate(
          Fraction = Count / sum(Count),
          LabelText = ifelse(Fraction >= 0.01, paste0(Count, " (", round(Fraction * 100, 1), "%)"), "")
        ) %>%
        dplyr::arrange(dplyr::desc(Simplified)) %>%
        dplyr::mutate(
          ymax = cumsum(Fraction),
          ymin = c(0, head(ymax, n = -1)),
          ypos = (ymax + ymin) / 2,
          hjust = ifelse(ypos < 0.5, 0, 1)
        )

      p <- ggplot2::ggplot(stats, ggplot2::aes(y = Fraction, fill = Simplified)) +
        ggplot2::geom_bar(ggplot2::aes(x = 1), width = 1, stat = "identity", color = "white") +
        ggplot2::geom_segment(
          data = subset(stats, LabelText != ""),
          ggplot2::aes(x = 1.51, xend = 1.62, y = ypos, yend = ypos), color = "grey50", size = 0.5
        ) +
        ggplot2::geom_text(
          data = subset(stats, LabelText != ""),
          ggplot2::aes(x = 1.65, y = ypos, label = LabelText, hjust = hjust), size = 3.5, fontface = "bold", check_overlap = FALSE
        ) +
        ggplot2::coord_polar("y", start = 0) +
        ggplot2::xlim(0.5, 2.5) +
        ggplot2::scale_fill_brewer(palette = palette, name = "Genomic Feature") +
        ggplot2::theme_void() +
        ggplot2::labs(title = title) +
        ggplot2::theme(
          legend.position = "bottom",
          legend.box.spacing = ggplot2::unit(-2, "pt"),
          legend.margin = ggplot2::margin(t = -2, b = 0),
          plot.margin = ggplot2::margin(t = 5, r = 5, b = 2, l = 5),
          plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 12)
        ) +
        ggplot2::guides(fill = ggplot2::guide_legend(nrow = 2, byrow = TRUE))

      ggplot2::ggsave(output_path, p, width = 7, height = 5)
    }

    message("    Generating Pie Chart 0: All Anchors Genomic Distribution...")
    if (nrow(cluster_info) > 0) {
      try(draw_pie_with_outside_labels(cluster_info, "annotation", paste0(project_name, ": All Anchors Genomic Distribution"), file.path(out_dir, paste0(project_name, "_Basic_Anchor_Genomic_Distribution.pdf")), color_palette), silent = TRUE)
    }

    message("    Generating Pie Chart 1: All Targets Distribution...")
    try(draw_pie_with_outside_labels(bed_info, "annotation", paste0(project_name, ": All Targets Genomic Distribution"), file.path(out_dir, paste0(project_name, "_Target_Genomic_Distribution.pdf")), color_palette), silent = TRUE)

    message("    Generating Pie Chart 2: Loop-Connected Targets Distribution...")
    linked_bed <- bed_info %>% dplyr::filter(!is.na(Linked_Loop_IDs) & Linked_Loop_IDs != "" & Linked_Loop_IDs != "NA")
    if (nrow(linked_bed) > 0) {
      try(draw_pie_with_outside_labels(linked_bed, "annotation", paste0(project_name, ": Loop-Connected Targets Distribution"), file.path(out_dir, paste0(project_name, "_Target_Loop_Genomic_Distribution.pdf")), color_palette), silent = TRUE)
    } else {
      message("    No linked targets found, skipping pie chart 2.")
    }
  }

  message("Step 7: Exporting to Excel...")
  tryCatch(
    {
      loop_annotation_clean <- loop_annotation_final %>% dplyr::select(-any_of(c("a1_id", "a2_id")))
      wb <- openxlsx::createWorkbook()
      openxlsx::addWorksheet(wb, "Loop Annotation")
      openxlsx::writeData(wb, "Loop Annotation", loop_annotation_clean)
      openxlsx::addWorksheet(wb, "Anchor Annotation")
      openxlsx::writeData(wb, "Anchor Annotation", cluster_info)
      openxlsx::addWorksheet(wb, "Promoter Stats")
      openxlsx::writeData(wb, "Promoter Stats", promoter_centric_df)
      if (!is.null(distal_element_df)) {
        openxlsx::addWorksheet(wb, "Distal Element Stats")
        openxlsx::writeData(wb, "Distal Element Stats", distal_element_df)
      }
      if (!is.null(bed_info)) {
        openxlsx::addWorksheet(wb, "Target Annotation")
        openxlsx::writeData(wb, "Target Annotation", bed_info)
      }
      openxlsx::saveWorkbook(wb, file.path(out_dir, paste0(project_name, "_Basic_Results.xlsx")), overwrite = TRUE)
      message("    Excel file saved.")
    },
    error = function(e) warning("Export failed: ", e$message)
  )

  message("Analysis Complete.")
  return(list(anchor_annotation = cluster_info, loop_annotation = if (exists("loop_annotation_clean")) loop_annotation_clean else loop_annotation_final, promoter_centric_stats = promoter_centric_df, distal_element_stats = distal_element_df, target_annotation = bed_info))
}

#' @title Transcriptome-Guided 3D Network Refinement & Anchor Reclassification
#'
#' @description
#' A critical filtration engine that integrates quantitative RNA-seq data (e.g., TPM/FPKM)
#' with 3D structural data to transition from a "Physical Contact Network" to a
#' "Functionally Active Regulatory Network".
#'
#' It systematically evaluates every loop anchor, purging transcriptional noise,
#' splitting conjoined multi-gene loci, and biologically reclassifying silent promoters into enhancer-like regulatory elements.
#'
#' @details
#' \strong{Core Algorithmic Innovations:}
#' \itemize{
#'   \item \strong{Multi-Gene Parsing & Precision Filtration:}
#'     Flawlessly handles merged gene assignments (e.g., \code{"GeneA;GeneB"}) generated by upstream biotype-aware annotations. It splits, individually evaluates against the expression threshold, and securely recombines the surviving active targets.
#'   \item \strong{Biologically-Driven Downgrading (eP / eG):}
#'     A purely physical Promoter (\code{P}) loop anchor that exhibits no active transcription is biologically reclassified as an enhancer-like element (\code{eP}). Gene bodies without transcription become \code{eG}. This corrects the regulatory syntax (e.g., turning a silent \code{P-P} loop into a functionally accurate \code{eP-P} enhancer-promoter loop).
#'   \item \strong{Topological Hub Inheritance (Crucial):}
#'     Unlike naive filters that recalculate node degrees from scratch, this engine \strong{inherits} the foundational 3D Hub classifications (e.g., \code{Is_High_Connectivity_Gene}) from the raw physical data. A structural hub remains structurally important; this function simply annotates its functional activation state, preventing the network from collapsing due to tissue-specific silencing.
#'   \item \strong{Automated Target Purification:}
#'     Intelligently detects external mapping columns (e.g., \code{Assigned_Target_Genes_Filled}) and filters them based on expression, ensuring that auxiliary BED targets are strictly linked only to functionally verified active genes.
#' }
#'
#' @param annotation_res List. The raw foundational output object returned by \code{\link{annotate_peaks_and_loops}}.
#' @param expr_matrix_file Character. Path to the normalized expression matrix (rows = genes, columns = samples). Supports CSV or TSV.
#' @param sample_columns Character vector or Integer indices. The specific sample columns in the matrix to average for calculating the baseline expression.
#' @param threshold Numeric. The minimum expression value (e.g., TPM > 1) required to consider a gene "biologically active" (default: \code{1}).
#' @param unit_type Character. Expression unit used for plot labeling (e.g., \code{"TPM"}; default: \code{"TPM"}).
#' @param species Character. Reference genome build. Supported: \code{"hg38"}, \code{"hg19"}, \code{"mm10"}, \code{"mm9"}.
#' @param out_dir Character. Directory path to save the generated PDF Karyotypes, Rose plots, and Excel results.
#' @param project_name Character. Prefix for all output files (automatically appends \code{"_Filtered"} if not present).
#' @param color_palette Character. Color brewer palette name for global visualization (default: \code{"Set2"}).
#' @param karyo_bin_size Numeric. Bin size for Karyotype genomic density heatmaps (default: \code{1e5}).
#' @param reclassify_by_expression Logical. If \code{TRUE} (default), mathematically silent Promoters (\code{P}) and Gene Bodies (\code{G}) are downgraded to \code{eP} and \code{eG}, preserving the 3D loop but correcting its topological class.
#' @param hub_percentile Numeric (0-1). Top percentile threshold used as a fallback for calculating Hub metrics if upstream inheritance is missing (default: \code{0.95}).
#'
#' @return An invisible list containing the refined multi-omics datasets (also exported as \code{_Refined_Results.xlsx}):
#' \itemize{
#'   \item \code{loop_annotation}: The fully filtered 3D network with updated \code{loop_type} (e.g., eP-P) and active \code{Putative_Target_Genes}.
#'   \item \code{anchor_annotation}: Cluster annotations updated with expressed localized targets.
#'   \item \code{promoter_centric_stats}: Inherited and updated topological statistics for regulatory targets.
#'   \item \code{distal_element_stats}: Enhancer-level statistics preserving high-connectivity enhancer status.
#'   \item \code{target_annotation}: External features strictly linked to verified active loop components.
#' }
#'
#' @importFrom dplyr %>% filter group_by summarise ungroup mutate select rename left_join full_join arrange desc case_when rowwise coalesce any_of distinct pull
#' @importFrom ggplot2 ggplot aes geom_bar geom_segment geom_point geom_text scale_color_manual scale_fill_manual theme_minimal theme_void labs coord_polar xlim ggsave
#' @importFrom tidyr pivot_longer separate_rows
#' @importFrom GenomicRanges makeGRangesFromDataFrame findOverlaps
#' @importFrom openxlsx createWorkbook addWorksheet writeData saveWorkbook
#' @export
#'
#' @examples
#' # 1. Get paths to the required example files in the package
#' rdata_path <- system.file("extdata", "analysis_results.RData", package = "looplook")
#' expr_path <- system.file("extdata", "example_tpm.txt", package = "looplook")
#'
#' # 2. Check if files exist before running
#' if (rdata_path != "" && expr_path != "") {
#'   # Safely load the pre-computed annotation result from RData
#'   temp_env <- new.env()
#'   load(rdata_path, envir = temp_env)
#'   # Extract the first object found in the RData file
#'   raw_annotation <- temp_env[[ls(temp_env)[1]]]
#'
#'   # =========================================================================
#'   # Example A: Standard filtering WITHOUT loop reclassification
#'   # =========================================================================
#'   res_basic <- refine_loop_anchors_by_expression(
#'     annotation_res = raw_annotation,
#'     expr_matrix_file = expr_path,
#'     sample_columns = c("con1", "con2"),
#'     threshold = 1.0,
#'     unit_type = "TPM",
#'     species = "hg38",
#'     out_dir = tempdir(),
#'     project_name = "Example_Basic",
#'     reclassify_by_expression = FALSE
#'   )
#'
#'   print(table(res_basic$loop_annotation$loop_type))
#'
#'   # =========================================================================
#'   # Example B: Advanced filtering WITH Transcriptome-Guided Reclassification
#'   # =========================================================================
#'   res_reclassified <- refine_loop_anchors_by_expression(
#'     annotation_res = raw_annotation,
#'     expr_matrix_file = expr_path,
#'     sample_columns = c("con1", "con2"),
#'     threshold = 1.0,
#'     unit_type = "TPM",
#'     species = "hg38",
#'     out_dir = tempdir(),
#'     project_name = "Example_Reclassified",
#'     reclassify_by_expression = TRUE
#'   )
#'
#'   # View the biologically corrected loop types (e.g., transition from P-P to eP-P)
#'   print(table(res_reclassified$loop_annotation$loop_type))
#' }
refine_loop_anchors_by_expression <- function(
  annotation_res,
  expr_matrix_file,
  sample_columns,
  threshold = 1,
  unit_type = "TPM",
  species = "hg38",
  out_dir = "./results/filtered",
  project_name = "HiChIP",
  color_palette = "Set2",
  karyo_bin_size = 1e5,
  reclassify_by_expression = TRUE,
  hub_percentile = 0.95
) {
  # --- 0. Setup ---
  if (!grepl("_Filtered$", project_name)) project_name <- paste0(project_name, "_Filtered")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  message(">>> [Refinement] Project Name: ", project_name)

  generate_colors <- function(n, pal_input) {
    if (length(pal_input) == 1 && pal_input %in% rownames(RColorBrewer::brewer.pal.info)) {
      max_pal <- RColorBrewer::brewer.pal.info[pal_input, "maxcolors"]
      base_cols <- RColorBrewer::brewer.pal(max(3, min(max_pal, n)), pal_input)
      if (n > length(base_cols)) {
        return(colorRampPalette(base_cols)(n))
      } else {
        return(base_cols[seq_len(n)])
      }
    } else {
      if (length(pal_input) < n) {
        return(colorRampPalette(pal_input)(n))
      } else {
        return(pal_input[seq_len(n)])
      }
    }
  }

  # --- 1. Load Data ---
  message(">>> [Step 1] Loading Data & Expression Matrix...")
  if (is.null(annotation_res$loop_annotation)) stop("'loop_annotation' missing.")

  original_loop_df <- annotation_res$loop_annotation
  loop_df <- annotation_res$loop_annotation
  clust_info <- annotation_res$anchor_annotation
  bed_info <- annotation_res$target_annotation

  # Check and reconstruct IDs if missing
  if (!"a1_id" %in% colnames(loop_df) || !"a2_id" %in% colnames(loop_df)) {
    message("    [Info] 'a1_id'/'a2_id' columns missing. Reconstructing from coordinates...")
    loop_df <- loop_df %>%
      dplyr::mutate(
        a1_id = paste(chr1, start1, end1, sep = "_"),
        a2_id = paste(chr2, start2, end2, sep = "_")
      )
  }

  # Load upstream stats (This holds the ORIGINAL structural connectivity)
  upstream_promoter_stats <- annotation_res$promoter_centric_stats
  upstream_distal_stats <- annotation_res$distal_element_stats

  if (!file.exists(expr_matrix_file)) stop("Expression file not found.")
  d <- tryCatch(read.table(expr_matrix_file, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE), error = function(e) NULL)
  if (is.null(d)) d <- tryCatch(read.table(expr_matrix_file, header = TRUE, sep = ",", row.names = 1, check.names = FALSE), error = function(e) NULL)
  if (is.null(d)) stop("Could not read expression matrix.")

  if (is.character(sample_columns)) sub_mat <- d[, intersect(sample_columns, colnames(d)), drop = FALSE] else sub_mat <- d[, sample_columns, drop = FALSE]
  if (ncol(sub_mat) == 0) stop("No valid sample columns found.")

  vals <- if (ncol(sub_mat) > 1) rowMeans(sub_mat, na.rm = TRUE) else sub_mat[, 1]
  whitelist <- names(vals)[vals > threshold & !is.na(vals) & names(vals) != ""]
  message(sprintf("    >>> Active Genes (> %s %s): %d", threshold, unit_type, length(whitelist)))

  # --- 2. Update Anchors & Loops ---
  message(">>> [Step 2] Updating Anchors & Loops...")

  clean_anchor <- function(g, t, allow = whitelist, down = reclassify_by_expression) {
    g_char <- as.character(g)
    t_char <- as.character(t)
    if (is.na(g_char) || g_char == "") {
      return(list(type = t_char, gene = NA_character_))
    }
    gs <- unlist(strsplit(g_char, ";"))
    active_gs <- gs[trimws(gs) %in% allow]
    if (length(active_gs) > 0) {
      return(list(type = t_char, gene = paste(unique(active_gs), collapse = ";")))
    } else {
      if (down) {
        new_type <- dplyr::case_when(t_char == "P" ~ "eP", t_char == "G" ~ "eG", TRUE ~ t_char)
        return(list(type = new_type, gene = NA_character_))
      } else {
        return(list(type = t_char, gene = NA_character_))
      }
    }
  }

  a1_res <- mapply(clean_anchor, loop_df$anchor1_gene, loop_df$anchor1_type, SIMPLIFY = FALSE)
  a2_res <- mapply(clean_anchor, loop_df$anchor2_gene, loop_df$anchor2_type, SIMPLIFY = FALSE)
  loop_df$anchor1_type <- vapply(a1_res, function(x) x$type, character(1))
  loop_df$anchor1_gene <- vapply(a1_res, function(x) x$gene, character(1))
  loop_df$anchor2_type <- vapply(a2_res, function(x) x$type, character(1))
  loop_df$anchor2_gene <- vapply(a2_res, function(x) x$gene, character(1))
  loop_df$loop_type <- mapply(function(t1, t2) paste(sort(c(t1, t2)), collapse = "-"), loop_df$anchor1_type, loop_df$anchor2_type)

  extract_genes <- function(genes_vec) {
    res <- unique(na.omit(unlist(strsplit(as.character(genes_vec), ";"))))
    if (length(res) == 0) {
      return(NA_character_)
    }
    paste(res, collapse = ";")
  }
  is_enh_like <- function(t) {
    t %in% c("E", "eP", "eG")
  }
  is_promoter <- function(t) {
    t == "P"
  }

  loop_df <- loop_df %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      Putative_Target_Genes = dplyr::case_when(
        (is_promoter(anchor1_type) & is_enh_like(anchor2_type)) ~ extract_genes(anchor1_gene),
        (is_enh_like(anchor1_type) & is_promoter(anchor2_type)) ~ extract_genes(anchor2_gene),
        (is_promoter(anchor1_type) & is_promoter(anchor2_type)) ~ extract_genes(c(anchor1_gene, anchor2_gene)),
        TRUE ~ extract_genes(c(anchor1_gene, anchor2_gene))
      )
    ) %>%
    dplyr::ungroup()

  # --- 3. Stats Update ---
  message(">>> [Step 3] Updating Stats...")

  # 3.1 Update Clusters
  agg_cluster <- loop_df %>%
    dplyr::filter(!is.na(cluster_id)) %>%
    dplyr::group_by(cluster_id) %>%
    dplyr::summarise(Cluster_All_Genes = extract_genes(Putative_Target_Genes), .groups = "drop")
  loop_df <- loop_df %>%
    dplyr::select(-any_of("Cluster_All_Genes")) %>%
    dplyr::left_join(agg_cluster, by = "cluster_id")
  if (!is.null(clust_info)) {
    clust_info <- clust_info %>%
      dplyr::select(-any_of("Cluster_All_Genes")) %>%
      dplyr::left_join(agg_cluster, by = "cluster_id")
  }

  get_dom <- function(x) {
    if (length(x) == 0) {
      return(NA_character_)
    }
    names(which.max(table(x)))
  }
  get_expr <- function(g) {
    e <- vals[g]
    e[is.na(e)] <- 0
    return(e)
  }

  # 3.2 Update Promoter Stats (MERGED Enhancers/GeneBodies -> n_Linked_Distal)
  raw_stats_df <- dplyr::bind_rows(
    loop_df %>% dplyr::filter(anchor1_type == "P" & !is.na(anchor1_gene)) %>% dplyr::select(Gene = anchor1_gene, Neighbor_Type = anchor2_type, Loop_Type = loop_type),
    loop_df %>% dplyr::filter(anchor2_type == "P" & !is.na(anchor2_gene)) %>% dplyr::select(Gene = anchor2_gene, Neighbor_Type = anchor1_type, Loop_Type = loop_type)
  ) %>%
    dplyr::group_by(Gene) %>%
    dplyr::summarise(
      Total_Loops_Filtered = dplyr::n(),
      n_Linked_Promoters_Filtered = sum(Neighbor_Type == "P", na.rm = TRUE),
      n_Linked_Distal_Filtered = sum(Neighbor_Type %in% c("E", "eP", "eG", "G"), na.rm = TRUE),
      Dominant_Interaction_Filtered = get_dom(Loop_Type),
      .groups = "drop"
    )

  promoter_centric_df <- NULL
  if (nrow(raw_stats_df) > 0) {
    if (!is.null(upstream_promoter_stats)) {
      promoter_centric_df <- upstream_promoter_stats %>%
        dplyr::left_join(raw_stats_df, by = "Gene") %>%
        dplyr::mutate(
          Total_Loops = coalesce(Total_Loops_Filtered, 0),
          n_Linked_Promoters = coalesce(n_Linked_Promoters_Filtered, 0),
          n_Linked_Distal = coalesce(n_Linked_Distal_Filtered, 0),
          Dominant_Interaction = coalesce(Dominant_Interaction_Filtered, "None"),
          Mean_Expression_Temp = get_expr(Gene),
          Is_Active_Gene = dplyr::if_else(Mean_Expression_Temp > threshold, "Yes", "No")
        ) %>%
        dplyr::select(Gene, Total_Loops, n_Linked_Promoters, n_Linked_Distal, Dominant_Interaction, any_of(c("Is_High_Connectivity_Gene", "Is_High_Distal_Connectivity_Gene")), Is_Active_Gene, everything()) %>%
        dplyr::select(-any_of(c(
          "Total_Loops_Filtered", "n_Linked_Promoters_Filtered",
          "n_Linked_Distal_Filtered", "Dominant_Interaction_Filtered",
          "Is_Regulatory_Hub", "Mean_Expression_Temp",
          "n_Linked_Enhancers", "n_Linked_GeneBodies"
        )))
    } else {
      final_cutoff <- max(quantile(raw_stats_df$Total_Loops_Filtered, hub_percentile, na.rm = TRUE), 3)
      distal_cutoff <- max(quantile(raw_stats_df$n_Linked_Distal_Filtered, hub_percentile, na.rm = TRUE), 2)

      promoter_centric_df <- raw_stats_df %>%
        dplyr::rename(
          Total_Loops = Total_Loops_Filtered,
          n_Linked_Promoters = n_Linked_Promoters_Filtered,
          n_Linked_Distal = n_Linked_Distal_Filtered,
          Dominant_Interaction = Dominant_Interaction_Filtered
        ) %>%
        dplyr::mutate(
          Mean_Expression_Temp = get_expr(Gene),
          Is_Active_Gene = dplyr::if_else(Mean_Expression_Temp > threshold, "Yes", "No"),
          Is_High_Connectivity_Gene = dplyr::if_else(Total_Loops >= final_cutoff, "Yes", "No"),
          Is_High_Distal_Connectivity_Gene = dplyr::if_else(n_Linked_Distal >= distal_cutoff, "Yes", "No")
        ) %>%
        dplyr::select(Gene, Total_Loops, n_Linked_Promoters, n_Linked_Distal, Dominant_Interaction, Is_High_Connectivity_Gene, Is_High_Distal_Connectivity_Gene, Is_Active_Gene, everything()) %>%
        dplyr::select(-any_of("Mean_Expression_Temp"))
    }
    promoter_centric_df <- promoter_centric_df %>% dplyr::arrange(dplyr::desc(n_Linked_Distal))
  }


  distal_element_df <- NULL
  if ("a1_id" %in% colnames(loop_df)) {
    distal_raw_df <- dplyr::bind_rows(
      loop_df %>% dplyr::filter(anchor1_type %in% c("E", "eP", "eG")) %>% dplyr::select(Distal_Anchor_ID = a1_id, Neighbor_Type = anchor2_type, Loop_Type = loop_type, Neighbor_Gene = anchor2_gene),
      loop_df %>% dplyr::filter(anchor2_type %in% c("E", "eP", "eG")) %>% dplyr::select(Distal_Anchor_ID = a2_id, Neighbor_Type = anchor1_type, Loop_Type = loop_type, Neighbor_Gene = anchor1_gene)
    ) %>%
      dplyr::group_by(Distal_Anchor_ID) %>%
      dplyr::summarise(
        Total_Loops_Filtered = dplyr::n(),
        n_Linked_Distal_Filtered = sum(Neighbor_Type %in% c("E", "eP", "eG", "G"), na.rm = TRUE),
        n_Linked_Promoters_Filtered = sum(Neighbor_Type == "P", na.rm = TRUE),
        Dominant_Interaction_Filtered = get_dom(Loop_Type),
        Target_Genes_Filtered = extract_genes(Neighbor_Gene[Neighbor_Type == "P"]),
        .groups = "drop"
      )

    anchor_map <- dplyr::bind_rows(
      loop_df %>% dplyr::select(anchor_id = a1_id, chr = chr1, start = start1, end = end1, cluster_id),
      loop_df %>% dplyr::select(anchor_id = a2_id, chr = chr2, start = start2, end = end2, cluster_id)
    ) %>% dplyr::distinct()

    if (nrow(distal_raw_df) > 0) {
      if (!is.null(upstream_distal_stats) && "Distal_Anchor_ID" %in% colnames(upstream_distal_stats)) {
        common_ids <- intersect(upstream_distal_stats$Distal_Anchor_ID, distal_raw_df$Distal_Anchor_ID)
        if (length(common_ids) > 0) {
          temp_df <- upstream_distal_stats %>%
            dplyr::left_join(distal_raw_df, by = "Distal_Anchor_ID") %>%
            dplyr::mutate(
              Total_Loops = coalesce(Total_Loops_Filtered, 0),
              n_Linked_Distal = coalesce(n_Linked_Distal_Filtered, 0),
              n_Linked_Promoters = coalesce(n_Linked_Promoters_Filtered, 0),
              Dominant_Interaction = coalesce(Dominant_Interaction_Filtered, "None"),
              Target_Genes = coalesce(Target_Genes_Filtered, "")
            )
        } else {
          final_cutoff_dist <- max(quantile(distal_raw_df$Total_Loops_Filtered, hub_percentile, na.rm = TRUE), 3)
          temp_df <- distal_raw_df %>%
            dplyr::rename(
              Total_Loops = Total_Loops_Filtered,
              n_Linked_Distal = n_Linked_Distal_Filtered,
              n_Linked_Promoters = n_Linked_Promoters_Filtered,
              Dominant_Interaction = Dominant_Interaction_Filtered,
              Target_Genes = Target_Genes_Filtered
            ) %>%
            dplyr::mutate(Is_High_Connectivity_Distal_Element = dplyr::if_else(Total_Loops >= final_cutoff_dist, "Yes", "No"))
        }
      } else {
        final_cutoff_dist <- max(quantile(distal_raw_df$Total_Loops_Filtered, hub_percentile, na.rm = TRUE), 3)
        temp_df <- distal_raw_df %>%
          dplyr::rename(
            Total_Loops = Total_Loops_Filtered,
            n_Linked_Distal = n_Linked_Distal_Filtered,
            n_Linked_Promoters = n_Linked_Promoters_Filtered,
            Dominant_Interaction = Dominant_Interaction_Filtered,
            Target_Genes = Target_Genes_Filtered
          ) %>%
          dplyr::mutate(Is_High_Connectivity_Distal_Element = dplyr::if_else(Total_Loops >= final_cutoff_dist, "Yes", "No"))
      }

      temp_df <- temp_df %>% dplyr::select(-any_of(c("chr", "start", "end", "cluster_id", "Distal_Type", "Distal_Type_Filtered", "Total_Loops_Filtered", "Target_Genes_Filtered", "n_Linked_Distal_Filtered", "n_Linked_Promoters_Filtered", "Dominant_Interaction_Filtered")))

      distal_element_df <- temp_df %>%
        dplyr::left_join(anchor_map, by = c("Distal_Anchor_ID" = "anchor_id")) %>%
        dplyr::select(chr, start, end, cluster_id, Total_Loops, n_Linked_Promoters, n_Linked_Distal, Dominant_Interaction, any_of("Is_High_Connectivity_Distal_Element"), Target_Genes) %>%
        dplyr::filter(Total_Loops > 0) %>%
        dplyr::arrange(dplyr::desc(Total_Loops))
    }
  }


  message(">>> [Step 4] Refining Target Annotations...")
  if (!is.null(bed_info)) {
    cols_to_clean <- grep("Priority|Strict|Physical|Loop_Genes|Promoters|Filled|Target_Genes|Assigned", colnames(bed_info), value = TRUE)
    raw_tgt_col <- "Target_Priority_Genes_Filled"
    if (!raw_tgt_col %in% colnames(bed_info)) {
      raw_tgt_col <- "Assigned_Target_Genes_Filled"
      if (!raw_tgt_col %in% colnames(bed_info)) {
        raw_tgt_col <- grep("Filled", cols_to_clean, value = TRUE)[1]
      }
    }
    if (!is.na(raw_tgt_col) && raw_tgt_col %in% colnames(bed_info)) {
      bed_info$SANKEY_RAW_GENES <- bed_info[[raw_tgt_col]]
    }
    for (col in cols_to_clean) {
      if (col %in% colnames(bed_info)) {
        bed_info[[col]] <- vapply(as.character(bed_info[[col]]), function(x) {
          if (is.na(x) || x == "" || x == "NA") {
            return(NA_character_)
          }
          gs <- unlist(strsplit(x, ";"))
          gs_active <- gs[trimws(gs) %in% whitelist]
          if (length(gs_active) == 0) {
            return(NA_character_)
          }
          return(paste(unique(sort(trimws(gs_active))), collapse = ";"))
        }, FUN.VALUE = character(1))
      }
    }
  }

  # --- 5. Visualization ---
  message(">>> [Step 5] Generating Visualizations...")
  red_palette <- c("#FFFFFF", "#FFFFCC", "#FFEDA0", "#FED976", "#FEB24C", "#FD8D3C", "#FC4E2A", "#E31A1C", "#BD0026", "#800026", "#000000")
  purple_palette <- c("#FFFFFF", "#F3E5F5", "#E1BEE7", "#BA68C8", "#9C27B0", "#7B1FA2", "#4A148C", "#000000")

  all_types <- sort(unique(c(unique(original_loop_df$loop_type), unique(loop_df$loop_type))))
  custom_colors <- generate_colors(length(all_types), color_palette)
  names(custom_colors) <- all_types

  # 5.1 Comparison Dumbbell
  if (requireNamespace("ggplot2", quietly = TRUE) && requireNamespace("dplyr", quietly = TRUE) && requireNamespace("tidyr", quietly = TRUE)) {
    tryCatch(
      {
        df_orig <- original_loop_df %>%
          dplyr::group_by(loop_type) %>%
          dplyr::summarise(Original = dplyr::n(), .groups = "drop")
        df_filt <- loop_df %>%
          dplyr::group_by(loop_type) %>%
          dplyr::summarise(Filtered = dplyr::n(), .groups = "drop")
        df_dumbbell <- dplyr::full_join(df_orig, df_filt, by = "loop_type") %>%
          dplyr::mutate(Original = ifelse(is.na(Original), 0, Original), Filtered = ifelse(is.na(Filtered), 0, Filtered), is_e_type = grepl("e", loop_type)) %>%
          dplyr::arrange(is_e_type, dplyr::desc(Original))
        df_dumbbell$loop_type <- factor(df_dumbbell$loop_type, levels = rev(df_dumbbell$loop_type))
        df_long <- df_dumbbell %>% tidyr::pivot_longer(cols = c("Original", "Filtered"), names_to = "Source", values_to = "Count")

        p_dumb <- ggplot2::ggplot() +
          ggplot2::geom_segment(data = df_dumbbell, ggplot2::aes(y = loop_type, yend = loop_type, x = Original, xend = Filtered), color = "#b2b2b2", size = 0.8) +
          ggplot2::geom_point(data = df_long, ggplot2::aes(x = Count, y = loop_type, color = Source), size = 3) +
          ggplot2::scale_color_manual(values = c("Original" = "#999999", "Filtered" = "#E69F00")) +
          ggplot2::theme_minimal() +
          ggplot2::labs(title = paste0(project_name, ": Filtration Effect (Dumbbell)"), x = "Number of Loops", y = "Loop Type") +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"), legend.position = "top")
        ggplot2::ggsave(file.path(out_dir, paste0(project_name, "_Comparison_Dumbbell.pdf")), p_dumb, width = 8, height = 8)
      },
      error = function(e) warning("     Comparison Dumbbell failed: ", e$message)
    )
  }

  # 5.2 Target Loop Donut
  if (!is.null(bed_info) && requireNamespace("ggplot2", quietly = TRUE) && requireNamespace("dplyr", quietly = TRUE)) {
    tryCatch(
      {
        gr_bed <- GenomicRanges::makeGRangesFromDataFrame(bed_info, keep.extra.columns = TRUE)
        active_anc <- dplyr::bind_rows(loop_df %>% dplyr::select(chr = chr1, start = start1, end = end1, cluster_id), loop_df %>% dplyr::select(chr = chr2, start = start2, end = end2, cluster_id)) %>% dplyr::distinct()
        if (nrow(active_anc) > 0) {
          gr_anc <- GenomicRanges::makeGRangesFromDataFrame(active_anc, keep.extra.columns = TRUE)
          hits <- GenomicRanges::findOverlaps(gr_bed, gr_anc)
          if (length(hits) > 0) {
            hit_ids <- unique(gr_anc$cluster_id[S4Vectors::subjectHits(hits)])
            tgt_loops <- loop_df %>% dplyr::filter(cluster_id %in% hit_ids)
            if (nrow(tgt_loops) > 0) {
              donut_data <- tgt_loops %>%
                dplyr::group_by(loop_type) %>%
                dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
                dplyr::mutate(fraction = count / sum(count), legend_label = paste0(loop_type, " (n=", count, ", ", round(fraction * 100, 1), "%)"), plot_label = loop_type, is_lower_e = grepl("^e", loop_type)) %>%
                dplyr::arrange(is_lower_e, dplyr::desc(count)) %>%
                dplyr::mutate(loop_type = factor(loop_type, levels = loop_type))
              p_donut <- ggplot2::ggplot(donut_data, ggplot2::aes(x = 2, y = count, fill = loop_type)) +
                ggplot2::geom_bar(stat = "identity", width = 1, color = "white") +
                ggplot2::coord_polar(theta = "y") +
                ggplot2::xlim(0.5, 2.9) +
                ggplot2::geom_text(ggplot2::aes(x = 2.8, label = plot_label), position = ggplot2::position_stack(vjust = 0.5), size = 3) +
                ggplot2::scale_fill_manual(values = custom_colors, labels = setNames(donut_data$legend_label, as.character(donut_data$loop_type)), name = "Loop Type") +
                ggplot2::theme_void() +
                ggplot2::labs(title = paste0(project_name, ": Loops Connected to Targets")) +
                ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14), legend.position = "right", legend.text = ggplot2::element_text(size = 10))
              ggplot2::ggsave(file.path(out_dir, paste0(project_name, "_Target_Loop_Donut.pdf")), p_donut, width = 8, height = 6)
            }
          }
        }
      },
      error = function(e) warning("    Donut Plot failed: ", e$message)
    )
  }

  # 5.3 Sankey Diagram
  if (!is.null(bed_info) && "SANKEY_RAW_GENES" %in% colnames(bed_info)) {
    message("    Drawing Target Sankey Diagram...")
    tryCatch(
      {
        if (requireNamespace("networkD3", quietly = TRUE)) {
          raw_bed <- bed_info
          total_targets <- nrow(raw_bed)
          get_label_mapping <- function(vec) {
            tbl <- table(vec)
            tbl <- tbl[tbl > 0]
            if (length(tbl) == 0) {
              return(character(0))
            }
            labels <- paste0(names(tbl), " (n=", tbl, ", ", round(as.numeric(tbl) / total_targets * 100, 1), "%)")
            names(labels) <- names(tbl)
            return(labels)
          }
          check_status_strict <- function(g_str) {
            if (is.na(g_str) || g_str == "" || g_str == "NA") {
              return("No Gene Assigned")
            }
            gs <- unlist(strsplit(as.character(g_str), ";"))
            gs <- trimws(gs)
            if (any(gs %in% whitelist)) {
              return("Active")
            }
            return("Inactive")
          }
          sankey_data_raw <- raw_bed %>%
            dplyr::mutate(L1_Raw = dplyr::case_when(grepl("Distal", annotation, ignore.case = TRUE) ~ "Distal Intergenic", grepl("Exon", annotation, ignore.case = TRUE) ~ "Exon", grepl("Intron", annotation, ignore.case = TRUE) ~ "Intron", grepl("Promoter", annotation, ignore.case = TRUE) ~ "Promoter", TRUE ~ "Others"), L2_Raw = ifelse(!is.na(Linked_Loop_IDs) & Linked_Loop_IDs != "" & Linked_Loop_IDs != "NA", "Connected", "Unconnected"), L3_Raw = vapply(SANKEY_RAW_GENES, check_status_strict, FUN.VALUE = character(1))) %>%
            dplyr::filter(L3_Raw != "No Gene Assigned")
          if (nrow(sankey_data_raw) > 0) {
            l1_map <- get_label_mapping(sankey_data_raw$L1_Raw)
            l2_map <- get_label_mapping(sankey_data_raw$L2_Raw)
            l3_map <- get_label_mapping(sankey_data_raw$L3_Raw)
            sankey_data_ready <- sankey_data_raw %>%
              dplyr::mutate(Genomic_Distribution = l1_map[L1_Raw], Loop_Connection = l2_map[L2_Raw], Expression_Status = l3_map[L3_Raw]) %>%
              dplyr::filter(!is.na(Genomic_Distribution) & !is.na(Loop_Connection) & !is.na(Expression_Status))
            if (nrow(sankey_data_ready) > 0) {
              links <- dplyr::bind_rows(sankey_data_ready %>% dplyr::group_by(source = Genomic_Distribution, target = Loop_Connection) %>% dplyr::summarise(value = dplyr::n(), .groups = "drop"), sankey_data_ready %>% dplyr::group_by(source = Loop_Connection, target = Expression_Status) %>% dplyr::summarise(value = dplyr::n(), .groups = "drop"))
              nodes_vec <- unique(c(links$source, links$target))
              nodes <- data.frame(name = nodes_vec, stringsAsFactors = FALSE)
              links$IDsource <- match(links$source, nodes$name) - 1
              links$IDtarget <- match(links$target, nodes$name) - 1
              sankey_colors <- generate_colors(length(nodes_vec), color_palette)
              colourScale <- sprintf('d3.scaleOrdinal().range(["%s"])', paste(sankey_colors, collapse = '","'))
              sn <- networkD3::sankeyNetwork(Links = links, Nodes = nodes, Source = "IDsource", Target = "IDtarget", Value = "value", NodeID = "name", units = "Targets", fontSize = 14, fontFamily = "Arial", nodeWidth = 15, nodePadding = 15, iterations = 0, height = 600, width = 900, colourScale = colourScale, sinksRight = FALSE)
              js_gradient <- 'function(el, x) { var svg = d3.select(el).select("svg"); function createValidID(name) { if (!name) return "unknown"; return name.replace(/[^a-zA-Z0-9-]/g, "_"); } svg.selectAll(".link").each(function(d) { var gradientID = "gradient-" + createValidID(d.source.name) + "-" + createValidID(d.target.name); if (svg.select("#" + gradientID).empty()) { var gradient = svg.append("defs").append("linearGradient").attr("id", gradientID).attr("gradientUnits", "userSpaceOnUse").attr("x1", d.source.x + d.source.dx / 2).attr("y1", d.source.y + d.source.dy / 2).attr("x2", d.target.x + d.target.dx / 2).attr("y2", d.target.y + d.target.dy / 2); var sourceColor = d3.select(el).selectAll(".node").filter(function(node) { return node.name === d.source.name; }).select("rect").style("fill"); var targetColor = d3.select(el).selectAll(".node").filter(function(node) { return node.name === d.target.name; }).select("rect").style("fill"); gradient.append("stop").attr("offset", "0%").attr("stop-color", sourceColor); gradient.append("stop").attr("offset", "100%").attr("stop-color", targetColor); } d3.select(this).style("stroke", "url(#" + gradientID + ")"); }); svg.selectAll(".node rect").style("stroke", "black").style("stroke-width", "1px"); }'
              sn_gradient <- htmlwidgets::onRender(sn, js_gradient)
              html_file <- file.path(out_dir, paste0(project_name, "_Target_Sankey.html"))
              htmlwidgets::saveWidget(sn_gradient, file = html_file, selfcontained = TRUE)
              if (requireNamespace("webshot2", quietly = TRUE) && Sys.getenv("_R_CHECK_PACKAGE_NAME_") == "") {
                webshot2::webshot(html_file, file = file.path(out_dir, paste0(project_name, "_Target_Sankey.pdf")), vwidth = 1100, vheight = 800, delay = 1)
              }
            }
          }
        }
      },
      error = function(e) warning("    Sankey Plot Failed: ", e$message)
    )
  }

  # 5.4 Karyotype
  if (exists("draw_karyo_heatmap_internal")) {
    tryCatch(
      {
        txdb_pkg <- if (species == "hg38") "TxDb.Hsapiens.UCSC.hg38.knownGene" else if (species == "mm10") "TxDb.Mmusculus.UCSC.mm10.knownGene" else NULL
        org_db <- if (species == "hg38") "org.Hs.eg.db" else if (grepl("mm", species)) "org.Mm.eg.db" else NULL
        if (!is.null(txdb_pkg) && !is.null(org_db)) {
          txdb_obj <- utils::getFromNamespace(txdb_pkg, txdb_pkg)
          all_genes <- GenomicFeatures::genes(txdb_obj)

          # === 修复核心：安全提取 org.db 对象 ===
          org_db_obj <- utils::getFromNamespace(org_db, org_db)
          map <- AnnotationDbi::select(org_db_obj, keys = as.character(all_genes$gene_id), columns = "SYMBOL", keytype = "ENTREZID")
          # ======================================

          all_genes$SYMBOL <- map$SYMBOL[match(all_genes$gene_id, map$ENTREZID)]
          g_active <- unique(trimws(unlist(strsplit(as.character(loop_df$Putative_Target_Genes), ";"))))
          g_active <- g_active[!is.na(g_active) & g_active != ""]
          if (length(g_active) > 0) draw_karyo_heatmap_internal(all_genes[all_genes$SYMBOL %in% g_active], "Refined Active Genes", file.path(out_dir, paste0(project_name, "_Refined_Karyo_Active.pdf")), karyo_bin_size, 0.99, txdb_obj, species, "Genes", custom_colors = red_palette)
          clean_tgt_col <- "Target_Priority_Genes_Filled"
          if (!clean_tgt_col %in% colnames(bed_info)) clean_tgt_col <- "Assigned_Target_Genes_Filled"
          if (!clean_tgt_col %in% colnames(bed_info)) clean_tgt_col <- grep("Filled", colnames(bed_info), value = TRUE)[1]
          if (!is.null(bed_info) && !is.na(clean_tgt_col)) {
            g_tgt <- unique(trimws(unlist(strsplit(as.character(bed_info[[clean_tgt_col]]), ";"))))
            g_tgt <- g_tgt[!is.na(g_tgt) & g_tgt != ""]
            if (length(g_tgt) > 0) draw_karyo_heatmap_internal(all_genes[all_genes$SYMBOL %in% g_tgt], "Refined Target Genes", file.path(out_dir, paste0(project_name, "_Refined_Karyo_TargetGenes.pdf")), karyo_bin_size, 0.99, txdb_obj, species, "Genes", custom_colors = purple_palette)
          }
        }
      },
      error = function(e) warning("    Karyo Plot failed: ", e$message)
    )
  }

  # 5.5 Rose Plot
  if (requireNamespace("ggplot2", quietly = TRUE) && requireNamespace("dplyr", quietly = TRUE)) {
    tryCatch(
      {
        rose_data <- loop_df %>%
          dplyr::group_by(loop_type) %>%
          dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
          dplyr::mutate(fraction = count / sum(count), legend_label = paste0(loop_type, " (n=", count, ", ", round(fraction * 100, 1), "%)"), is_lower_e = grepl("^e", loop_type)) %>%
          dplyr::arrange(dplyr::desc(count))
        plot_order <- rose_data$loop_type
        rose_data$loop_type <- factor(rose_data$loop_type, levels = plot_order)
        legend_order <- rose_data %>%
          dplyr::arrange(is_lower_e, dplyr::desc(count)) %>%
          dplyr::pull(loop_type)
        p_rose <- ggplot2::ggplot(rose_data, ggplot2::aes(x = loop_type, y = count, fill = loop_type)) +
          ggplot2::geom_bar(stat = "identity", width = 1, color = "white") +
          ggplot2::coord_polar(theta = "x") +
          ggplot2::scale_fill_manual(values = custom_colors, labels = setNames(rose_data$legend_label, as.character(rose_data$loop_type)), breaks = legend_order, name = "Loop Type") +
          ggplot2::theme_void() +
          ggplot2::labs(title = paste0(project_name, ": Loop Proportion (By Count)")) +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14), legend.position = "right", legend.text = ggplot2::element_text(size = 10))
        ggplot2::ggsave(file.path(out_dir, paste0(project_name, "_Rose.pdf")), p_rose, width = 8, height = 6)
      },
      error = function(e) warning("     Rose Plot failed: ", e$message)
    )
  }

  # --- 6. Export ---
  message(">>> [Step 6] Exporting Refined Results...")
  tryCatch(
    {
      wb <- openxlsx::createWorkbook()
      loop_export <- loop_df %>% dplyr::select(-any_of(c("loop_genes", "single_loop_genes", "proximate_loop_gene")))
      openxlsx::addWorksheet(wb, "Filtered Loop Annotation")
      openxlsx::writeData(wb, "Filtered Loop Annotation", loop_export)
      openxlsx::addWorksheet(wb, "Filtered Anchor Annotation")
      openxlsx::writeData(wb, "Filtered Anchor Annotation", clust_info)

      if (!is.null(promoter_centric_df)) {
        openxlsx::addWorksheet(wb, "Filtered Promoter Stats")
        openxlsx::writeData(wb, "Filtered Promoter Stats", promoter_centric_df)
      }
      if (!is.null(distal_element_df)) {
        openxlsx::addWorksheet(wb, "Filtered Distal Stats")
        openxlsx::writeData(wb, "Filtered Distal Stats", distal_element_df)
      }

      bed_export <- bed_info %>% dplyr::select(-any_of("SANKEY_RAW_GENES"))
      if (!is.null(bed_info)) {
        openxlsx::addWorksheet(wb, "Filtered Target Annotation")
        openxlsx::writeData(wb, "Filtered Target Annotation", bed_export)
      }

      openxlsx::saveWorkbook(wb, file.path(out_dir, paste0(project_name, "_Refined_Results.xlsx")), overwrite = TRUE)
      message("    Excel saved.")
    },
    error = function(e) warning("    Export Failed: ", e$message)
  )

  message("Refinement Complete.")
  return(list(
    loop_annotation = loop_df,
    anchor_annotation = clust_info,
    promoter_centric_stats = promoter_centric_df,
    distal_element_stats = distal_element_df,
    target_annotation = bed_info
  ))
}


#' @title Standardize and Clean ChIPseeker Annotations
#'
#' @description
#' An internal helper function that parses the verbose annotation strings generated by
#' \code{ChIPseeker::annotatePeak}. It extracts the broad genomic feature category
#' while preserving the exact spatial details.
#'
#' @details
#' \code{ChIPseeker} often outputs annotations with highly specific distance or transcript
#' information, such as \code{"Promoter (<=1kb)"} or \code{"Intron (uc001.1/exon 1)"}.
#' This function creates a clean, categorical \code{annotation} column (e.g., \code{"Promoter"}, \code{"Intron"})
#' which is strictly required for robust downstream regular expression matching and Pie/Donut chart visualizations,
#' while safely moving the original verbose string to a new \code{detail_anno} column.
#'
#' @param df A data frame representation of a \code{csAnno} object (generated by \code{as.data.frame(annotatePeak(...))}).
#'
#' @return A modified data frame where:
#' \itemize{
#'   \item \code{annotation} contains the broad feature class.
#'   \item \code{detail_anno} contains the original verbose string.
#' }
#'
#' @keywords internal
format_annotation_columns <- function(df) {
  if ("annotation" %in% colnames(df)) {
    df <- df %>%
      dplyr::rename(detail_anno = annotation) %>%
      dplyr::mutate(annotation = gsub(" \\(.*", "", detail_anno)) %>%
      dplyr::relocate(annotation, .before = detail_anno)
  }
  return(df)
}
