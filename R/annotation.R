#' Internal: Collapse IDs
#' @return A semicolon-delimited character string of unique, non-NA IDs.
#' @keywords internal
#' @noRd
.annotation_extract_ids <- function(id_vec) {
    paste(unique(na.omit(as.character(id_vec))), collapse = ";")
}

#' Internal: Resolve AnnotationDb Package Name
#' @keywords internal
#' @noRd
.pkg_from_annotation_db <- function(x) {
    db_path <- tryCatch(AnnotationDbi::dbfile(x), error = function(e) NULL)
    if (!is.null(db_path) && length(db_path) == 1L && nzchar(db_path)) {
        return(basename(dirname(dirname(db_path))))
    }
    pkg_attr <- attr(x, "package")
    if (!is.null(pkg_attr) && length(pkg_attr) == 1L && nzchar(pkg_attr) &&
        pkg_attr != "AnnotationDbi") {
        return(pkg_attr)
    }
    NULL
}

#' Internal: Resolve Annotation Resource
#' @keywords internal
#' @noRd
.resolve_annotation_resource <- function(arg, type, desc, species) {
    if (any(inherits(arg, c("TxDb", "OrgDb", "AnnotationDb")))) {
        return(list(obj = arg, pkg = .pkg_from_annotation_db(arg)))
    }
    if (is.character(arg) && nzchar(arg)) {
        if (!requireNamespace(arg, quietly = TRUE)) stop(desc, " '", arg, "' not installed")
        return(list(obj = utils::getFromNamespace(arg, arg), pkg = arg))
    }
    pkg <- if (type == "txdb") species_txdb_pkg(species) else species_orgdb_pkg(species)
    if (!requireNamespace(pkg, quietly = TRUE)) stop(desc, " '", pkg, "' not installed")
    list(obj = utils::getFromNamespace(pkg, pkg), pkg = pkg)
}

#' Internal: Convert ChIPseeker Annotation to Anchor Class
#' @details
#' \strong{Important:} The classification below is purely \emph{positional}
#' (where the anchor sits relative to known gene annotations), not
#' \emph{functional}.  The label \code{"E"} means "distal / intergenic" and
#' does \strong{not} imply enhancer activity. Orthogonal evidence (ATAC-seq,
#' H3K27ac, EP300) is required for functional enhancer classification.
#' @return A single character: \code{"P"} (promoter), \code{"G"} (gene body),
#'   \code{"E"} (distal/intergenic -- positional, not functional enhancer),
#'   or \code{"Unknown"}.
#' @keywords internal
#' @noRd
.annotation_feature_class <- function(anno_str) {
    # Positional classification only: P = promoter, G = gene body,
    # E = distal/intergenic (NOT a functional enhancer).
    # Fallback "E" means the region does not overlap annotated genes.
    if (length(anno_str) == 0 || is.na(anno_str)) {
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

    "E"  # positional default for unrecognised annotation types;
         # ChIPseeker changes to annotation strings could silently
         # route novel types here.  Review after ChIPseeker upgrades.
}

#' Internal: Extract Loop Locus Genes
#' @return A semicolon-delimited character string of gene symbols from P/G anchors, or an empty string.
#' @keywords internal
#' @noRd
.loop_locus_genes <- function(t1, t2, s1, s2) {
    genes <- c()
    if (!is.na(t1) && (.is_promoter_like(t1) || .is_gene_body_like(t1))) genes <- c(genes, s1)
    if (!is.na(t2) && (.is_promoter_like(t2) || .is_gene_body_like(t2))) genes <- c(genes, s2)
    paste(unique(na.omit(genes)), collapse = ";")
}

#' Internal: Build Loop Type Code
#' @details Anchors are sorted alphabetically so that lowercase prefixes
#'   (eP, eG) sort before uppercase (E, G, P).  This gives a canonical
#'   ordering where expression-filtered anchors appear first (e.g. eP-P
#'   rather than P-eP).
#' @return A two-letter code string (e.g. \code{"E-P"}, \code{"P-P"}) with anchors sorted alphabetically, or \code{"Unknown"}.
#' @keywords internal
#' @noRd
.loop_type_code <- function(t1, t2) {
    if (length(t1) == 0 || length(t2) == 0 || is.na(t1) || is.na(t2)) {
        return("Unknown")
    }
    paste(sort(c(t1, t2)), collapse = "-")
}

#' Internal: Convert Anchor IDs to Genes
#' @return A semicolon-delimited gene symbol string, or \code{NA_character_} if no valid mapping exists.
#' @keywords internal
#' @noRd
.ids_to_genes_simple <- function(ids, lookup) {
    valid <- intersect(ids, names(lookup))
    if (length(valid) == 0) {
        return(NA_character_)
    }
    genes <- lookup[valid]
    genes <- genes[!is.na(genes)]
    if (length(genes) == 0) {
        return(NA_character_)
    }
    paste(sort(unique(genes)), collapse = ";")
}

#' Internal: Convert Anchor IDs to Promoter-Priority Genes
#' @keywords internal
#' @noRd
.ids_to_genes_priority <- function(ids, lookup_sym, lookup_typ) {
    valid <- intersect(ids, names(lookup_sym))
    if (length(valid) == 0) {
        return(NA_character_)
    }
    promoter_ids <- valid[!is.na(lookup_typ[valid]) & .is_promoter_like(lookup_typ[valid])]
    use_ids <- if (length(promoter_ids) > 0) promoter_ids else valid
    genes_present <- lookup_sym[use_ids]
    genes_present <- genes_present[!is.na(genes_present)]
    if (length(genes_present) == 0) {
        return(NA_character_)
    }
    paste(sort(unique(genes_present)), collapse = ";")
}

#' Internal: Fill Empty Target Gene Assignments
#' @keywords internal
#' @noRd
.fill_target_gene_fallback <- function(target, fallback) {
    dplyr::case_when(
        !is.na(target) & target != "" ~ target,
        !is.na(fallback) & fallback != "" ~ fallback,
        TRUE ~ NA_character_
    )
}

#' Internal: Integrate Optional Target BED Annotation
#' @keywords internal
#' @noRd
.annotate_target_bed <- function(
  target_bed, txdb_obj, org_db_pkg, tss_region, gene_expr_map, min_expr,
  conflict_strategy,
  gr_anchors, anchor_topo_map, loop_annotation_final, map_info, ego_list_target,
  log_message,
  anchor_gap = -1L, anchor_min_overlap = 1L, anchor_min_frac = 0,
  co_dominance_ratio = 0.1
) {
    bed_target <- read_robust_general(target_bed, min_cols = 3, desc = "Target BED")
    colnames(bed_target)[c(1, 2, 3)] <- c("chr", "start", "end")
    if (nrow(bed_target) == 0) {
        warning("Target BED contains no features; skipping target annotation.")
        return(list(
            bed_info = NULL,
            target_connected_loops = NULL,
            target_gene_links = NULL
        ))
    }

    bed_target$start <- bed_target$start + 1 # BED is 0-based; GRanges is 1-based
    gr_bed <- .with_known_upstream_noise_suppressed({
        gr_bed <- GenomicRanges::makeGRangesFromDataFrame(bed_target)
        gr_bed$input_id <- paste0("Peak_", seq_len(nrow(bed_target)))
        names(gr_bed) <- gr_bed$input_id
        gr_bed
    })
    bed_annot <- .with_known_upstream_noise_suppressed(
        ChIPseeker::annotatePeak(gr_bed, TxDb = txdb_obj, tssRegion = tss_region, annoDb = org_db_pkg, verbose = FALSE)
    )
    bed_info <- format_annotation_columns(as.data.frame(bed_annot))
    if ("GENENAME" %in% colnames(bed_info)) bed_info <- bed_info %>% dplyr::rename(Gene_description = GENENAME)
    log_message("    Refining Target annotation...")
    bed_info <- resolve_gene_conflicts(bed_info, txdb_obj, org_db_pkg, tss_region, gene_expr_map, min_expr = min_expr, conflict_strategy = conflict_strategy, co_dominance_ratio = co_dominance_ratio)
    gr_bed <- .harmonize_seqlevels(gr_bed, gr_anchors, "target BED")
    n_peaks <- length(gr_bed)
    n_anchors <- length(gr_anchors)

    # anchor_gap = -1L: use findOverlaps default (strict physical overlap).
    # anchor_gap >= 0:  proximity-based matching -- intervals within `anchor_gap`
    #   bp of each other are considered hits even without physical overlap.
    #   This is the intended behavior for cross-experiment integration.
    #   The pintersect post-filter only applies when anchor_min_overlap > 1L
    #   (user explicitly requests a stricter bp-level filter).
    if (anchor_gap >= 0L) {
        hits <- GenomicRanges::findOverlaps(
            gr_bed, gr_anchors,
            maxgap = anchor_gap
        )
    } else {
        hits <- GenomicRanges::findOverlaps(gr_bed, gr_anchors)
    }

    # Post-filter by minimum physical overlap -- only when user explicitly
    # requests stricter bp-level filtering. Proximity-only hits (no actual
    # overlap but within anchor_gap) are NOT filtered here.
    n_hits_raw <- length(hits)
    if (anchor_min_overlap > 1L && n_hits_raw > 0) {
        q_gr <- gr_bed[S4Vectors::queryHits(hits)]
        s_gr <- gr_anchors[S4Vectors::subjectHits(hits)]
        overlap_w <- GenomicRanges::width(
            GenomicRanges::pintersect(q_gr, s_gr)
        )
        hits <- hits[overlap_w >= anchor_min_overlap]
    }

    # Fraction-based filtering: overlap must cover at least anchor_min_frac of the peak
    if (anchor_min_frac > 0 && length(hits) > 0) {
        q_gr <- gr_bed[S4Vectors::queryHits(hits)]
        s_gr <- gr_anchors[S4Vectors::subjectHits(hits)]
        overlap_w <- GenomicRanges::width(
            GenomicRanges::pintersect(q_gr, s_gr)
        )
        peak_w <- GenomicRanges::width(q_gr)
        frac <- overlap_w / peak_w
        hits <- hits[frac >= anchor_min_frac]
    }

    # Diagnostic: peak-to-anchor overlap summary
    n_peaks_hit <- length(unique(S4Vectors::queryHits(hits)))
    frac_info <- if (anchor_min_frac > 0)
        sprintf(", min_frac=%.2f", anchor_min_frac) else ""
    bp_info <- if (anchor_min_overlap > 1L)
        sprintf(", min_overlap=%s bp", format(anchor_min_overlap, big.mark = ",")) else ""
    log_message(sprintf(
        "    Peak-anchor overlap: %s / %s peaks (%.1f%%) hit anchors (gap=%s bp%s%s%s)",
        format(n_peaks_hit, big.mark = ","),
        format(n_peaks, big.mark = ","),
        n_peaks_hit / max(n_peaks, 1) * 100,
        format(anchor_gap, big.mark = ","),
        bp_info,
        frac_info,
        if (n_hits_raw != length(hits))
            sprintf(" [%s raw hits before filtering]", format(n_hits_raw, big.mark = ","))
        else ""
    ))
    if (n_peaks_hit == 0) {
        log_message("    No peaks overlapped loop anchors. Check that target BED and loop BEDPE use the same genome build and seqlevel style.")
    }

    target_connected_loops <- NULL
    target_gene_links <- NULL
    if (length(hits) > 0) {
        target_connected_loops <- loop_annotation_final %>%
            dplyr::filter(cluster_id %in% unique(gr_anchors$cluster_id[S4Vectors::subjectHits(hits)]))
        hit_df <- data.frame(qid = S4Vectors::queryHits(hits), sid = S4Vectors::subjectHits(hits))
        hit_df$anchor_id <- gr_anchors$anchor_id[hit_df$sid]
        hit_df <- hit_df %>% dplyr::left_join(anchor_topo_map, by = "anchor_id")
        anchor_loop_agg <- dplyr::bind_rows(
            loop_annotation_final %>% dplyr::select(anchor_id = a1_id, loop_ID),
            loop_annotation_final %>% dplyr::select(anchor_id = a2_id, loop_ID)
        ) %>%
            dplyr::distinct() %>%
            dplyr::group_by(anchor_id) %>%
            dplyr::summarise(linked_loops = .annotation_extract_ids(loop_ID), .groups = "drop")
        hit_df <- hit_df %>% dplyr::left_join(anchor_loop_agg, by = "anchor_id")
        summary_df <- hit_df %>%
            dplyr::group_by(qid) %>%
            dplyr::summarise(
                All_Loop_Connected_Genes = extract_genes(tgt_genes_pg),
                Regulated_promoter_genes = extract_genes(tgt_genes_p),
                Assigned_Target_Genes = extract_genes(tgt_genes_prio),
                Linked_Loop_IDs = .annotation_extract_ids(linked_loops),
                .groups = "drop"
            ) %>%
            dplyr::mutate(join_id = paste0("Peak_", qid))
        bed_info <- dplyr::left_join(bed_info, summary_df, by = c("input_id" = "join_id")) %>%
            dplyr::select(-any_of(c("join_id", "qid")))
        target_gene_links <- .build_target_gene_links(
            hit_df = hit_df,
            bed_info = bed_info,
            loop_annotation_final = loop_annotation_final,
            map_info = map_info,
            ego_list_target = ego_list_target
        )
    } else {
        bed_info$All_Loop_Connected_Genes <- NA
        bed_info$Regulated_promoter_genes <- NA
        bed_info$Assigned_Target_Genes <- NA
        bed_info$Linked_Loop_IDs <- NA
    }

    if (is.null(target_gene_links)) {
        target_gene_links <- .build_target_gene_links(
            hit_df = data.frame(qid = integer(0), anchor_id = character(0)),
            bed_info = bed_info,
            loop_annotation_final = loop_annotation_final,
            map_info = map_info,
            ego_list_target = ego_list_target
        )
    }

    evidence_df <- .summarise_regulated_promoter_evidence(target_gene_links)
    bed_info <- dplyr::left_join(bed_info, evidence_df, by = "input_id")
    bed_info$Regulated_promoter_Evidence <- ifelse(
        is.na(bed_info$Regulated_promoter_Evidence) |
            bed_info$Regulated_promoter_Evidence == "",
        "none",
        bed_info$Regulated_promoter_Evidence
    )

    fallback_col <- .target_linear_gene_column(bed_info)
    fallback_vec <- if (!is.null(fallback_col)) bed_info[[fallback_col]] else rep(NA_character_, nrow(bed_info))
    ann_vec <- if ("annotation" %in% colnames(bed_info)) {
        bed_info$annotation
    } else {
        rep(NA_character_, nrow(bed_info))
    }
    fallback_evidence <- .fallback_evidence_from_annotation(ann_vec)

    bed_info <- bed_info %>% dplyr::mutate(
        All_Loop_Connected_Genes_Filled = .fill_target_gene_fallback(All_Loop_Connected_Genes, fallback_vec),
        Regulated_promoter_genes_Filled = .fill_target_gene_fallback(Regulated_promoter_genes, fallback_vec),
        Assigned_Target_Genes_Filled = .fill_target_gene_fallback(Assigned_Target_Genes, fallback_vec),
        Regulated_promoter_Fallback_Evidence = dplyr::case_when(
            !is.na(Regulated_promoter_genes) & Regulated_promoter_genes != "" ~ "none",
            !is.na(Regulated_promoter_genes_Filled) & Regulated_promoter_genes_Filled != "" ~
                fallback_evidence,
            TRUE ~ "none"
        )
    )
    target_gene_links <- .mark_target_gene_link_membership(target_gene_links, bed_info)

    if ("Linked_Loop_IDs" %in% colnames(bed_info)) {
        target_col <- if ("Gene_description" %in% colnames(bed_info)) "Gene_description" else "SYMBOL"
        if (target_col %in% colnames(bed_info)) {
            bed_info <- bed_info %>% dplyr::relocate(Linked_Loop_IDs, .after = dplyr::all_of(target_col))
        }
    }

    # Return the real peak-to-anchor hit_df (with populated anchor_ids) for
    # downstream chromatin-aware target remapping
    stored_hit_df <- if (exists("hit_df", inherits = FALSE) && !is.null(hit_df) && nrow(hit_df) > 0) {
        hit_df %>%
            dplyr::select(dplyr::any_of(c("qid", "sid", "anchor_id"))) %>%
            dplyr::distinct()
    } else {
        data.frame(qid = integer(0), sid = integer(0),
                   anchor_id = character(0),
                   stringsAsFactors = FALSE)
    }

    list(
        bed_info = bed_info,
        target_connected_loops = target_connected_loops,
        target_gene_links = target_gene_links,
        hit_df = stored_hit_df
    )
}

#' @title Spatial mapping of 1D genomic features to 3D chromatin interaction targets
#'
#' @description
#' A function for loop annotation and target mapping, designed to integrate 1D genomic features with 3D chromatin architecture.
#' \enumerate{
#'   \item \strong{Loop Annotation:} Classifies 3D spatial interactions (e.g., distal-to-promoter, promoter-to-promoter) using positional anchor labels relative to gene annotations. \strong{Important:} anchor type \code{"E"} denotes a \emph{positional} distal/intergenic classification -- it does \strong{not} imply functional enhancer activity. Orthogonal chromatin data are required for functional interpretation.
#'   \item \strong{Feature-to-Target Mapping:} Links 1D genomic features (e.g., GWAS risk SNPs, ATAC-seq peaks, ChIP-seq binding sites) to putative target genes via 3D chromatin contacts, providing a loop-based alternative to linear proximity-based assignments.
#' }
#'
#' @details
#' \strong{Mapping Strategy and Fallback Mechanism}
#' The method prioritizes physical 3D chromatin contacts while keeping strict
#' and fallback semantics separate. \code{Regulated_promoter_genes} reports
#' promoter genes supported by loop-anchor context, \code{Assigned_Target_Genes}
#' preserves the historical promoter-first 3D assignment, and \code{*_Filled}
#' columns add a linear nearest-gene fallback only when strict 3D assignments are
#' empty. \code{Regulated_promoter_Evidence},
#' \code{Regulated_promoter_Fallback_Evidence}, and \code{target_gene_links}
#' provide row-level provenance for these decisions.
#'
#' \strong{Hierarchical Conflict Resolution}
#' To address complex loci where a single anchor overlaps multiple promoters (e.g., dense gene clusters or bidirectional promoters), the function executes a 3-step resolution:
#' \enumerate{
#'   \item \emph{Biotype Prioritization:} Selects the highest-priority candidates by functional class: \code{Protein Coding > small-ncRNA (miRNA, snoRNA, snRNA, rRNA, scaRNA) > Antisense > lncRNA/ncRNA > Pseudogene}.
#'   \item \emph{Expression Filter:} Within the selected biotype tier, excludes transcriptionally silent genes using a user-provided expression matrix. If no gene in the tier is expressed, all candidates in that tier are retained.
#'   \item \emph{Expression Tiebreaker:} Among remaining candidates of equal biotype priority, retains all genes whose expression >= \code{co_dominance_ratio} x the group maximum. Default \code{0.1} (one order of magnitude; ~3.3-fold in log2 space). This co-dominant rule preserves functionally redundant candidates such as bidirectional promoter pairs, where co-expressed partners typically fall within 2-3 fold of each other. Edge case: when all candidates in the best biotype tier have TPM = 0, all are retained (the tiebreaker cannot distinguish them), which may produce multi-gene assignments at transcriptionally silent loci.
#' }
#'
#' \strong{Network Topology Analysis}
#' \itemize{
#'   \item \strong{Ego-Network Expansion (\code{neighbor_hop}):} Implements k-hop neighborhood expansion via \code{igraph::ego()}. A value of \code{0} restricts to direct contacts, while \code{1} includes secondary contacts to capture broader regulatory cliques.
#'   \item \strong{Hub Detection:} Utilizes a node-degree quantile threshold (\code{hub_percentile}) to identify highly connected regulatory elements.
#' }
#'
#' \strong{Peak-Anchor Overlap Control}
#' When \code{target_bed} is provided, three parameters control how target peaks
#' are matched to loop anchors. They act as a cascade of increasingly stringent
#' filters:
#' \enumerate{
#'   \item \code{anchor_gap} -- expands the search radius around each anchor.
#'   \item \code{anchor_min_overlap} -- requires a minimum physical overlap in bp.
#'   \item \code{anchor_min_frac} -- requires the overlap to cover a minimum fraction of the peak.
#' }
#' The table below lists suggested starting points for common experimental designs.
#'
#' \tabular{llll}{
#'   \strong{Experimental design} \tab \code{anchor_gap} \tab \code{anchor_min_overlap} \tab \code{anchor_min_frac} \cr
#'   Same-experiment HiChIP / ChIA-PET (default) \tab \code{0} \tab \code{1} \tab \code{0} \cr
#'   Cross-experiment (e.g. ATAC-seq peaks x HiChIP loops) \tab \code{200-500} \tab \code{10} \tab \code{0} \cr
#'   Broad histone-mark peaks (H3K27ac, 2-5 kb) \tab \code{0} \tab \code{100} \tab \code{0.1} \cr
#'   Super-enhancers / wide domains (20-80 kb) \tab \code{0} \tab \code{500} \tab \code{0.05} \cr
#'   Point features (GWAS SNPs, eQTLs, 1 bp) \tab \code{0} \tab \code{1} \tab \code{0} \cr
#'   Stringent (high-confidence only) \tab \code{0} \tab \code{100} \tab \code{0.5}
#' }
#'
#' When \code{quiet = FALSE} (default), the function prints a diagnostic line
#' reporting how many peaks overlapped loop anchors after filtering, helping you
#' tune these thresholds for your data.
#'
#' @param bedpe_file Character. Path to a BEDPE file (at least 6 columns: chr1, start1, end1, chr2, start2, end2).
#'   Additional columns beyond 6 are retained in the output; anchor swapping only affects columns 1-6.
#' @param target_bed Optional path to a BED file of genomic features (e.g., ChIP-seq peaks, GWAS SNPs). When provided, these 1D regions are mapped to 3D target genes. Default: \code{NULL}.
#' @param txdb A \code{\link[GenomicFeatures]{TxDb}} object, a package name string, or \code{NULL} to auto-resolve from \code{species}. Default: \code{NULL}.
#' @param org_db An \code{OrgDb} object, a package name string, or \code{NULL} to auto-resolve from \code{species}. Default: \code{NULL}.
#' @param species Character. Genome assembly used when \code{txdb} and \code{org_db} are \code{NULL}. One of \code{"hg38"}, \code{"hg19"}, \code{"mm10"}, \code{"mm9"}. Default: \code{"hg38"}. The package is currently tested on human and mouse; the architecture supports extension to other species by adding entries to \code{species_txdb_pkg()}, \code{species_orgdb_pkg()}, and \code{species_bsgenome_pkg()}.
#' @param tss_region Numeric vector of length 2. Promoter window around the TSS in bp. Default: \code{c(-2000, 2000)} (typical for mammalian protein-coding genes; may need widening for broad domains like HOX clusters, or narrowing for compact genomes).
#' @param out_dir Character. Output directory for the Excel results file. Default: \code{"./results"}.
#' @param expr_matrix_file Optional path to a normalised expression matrix (TPM/FPKM, genes x samples). Enables expression-aware conflict resolution. Default: \code{NULL}.
#' @param sample_columns Character vector or integer indices. Columns in \code{expr_matrix_file} to average for baseline expression. Default: \code{NULL}.
#' @param min_expr Numeric. Minimum expression value for a gene to be considered
#'   active during anchor-level conflict resolution. Used only when
#'   \code{expr_matrix_file} is provided. Default: \code{0} (any non-zero
#'   expression, i.e. TPM > 0). Increase to \code{1} or higher to require
#'   stronger evidence. Note: when \code{min_expr = 0}, the code uses a strict
#'   greater-than comparison (\code{> 0}) to exclude truly undetected genes;
#'   when \code{min_expr >= 1}, it uses \code{>= min_expr}.
#'   See \code{\link{refine_loop_anchors_by_expression}} for the
#'   separate \code{threshold} parameter that controls promoter reclassification.
#' @param conflict_strategy Character. Conflict resolution order for
#'   overlapping gene assignments. \code{"biotype_first"} (default): select
#'   the best biotype tier first, then apply expression filtering within that
#'   tier. \code{"expression_first"}: apply expression filtering across all
#'   biotypes first, then pick the best biotype among expressed candidates.
#' @param co_dominance_ratio Numeric (0-1). In the expression tiebreaker step,
#'   genes with expression >= \code{co_dominance_ratio x max(expression)} in the
#'   group are retained together. Default: \code{0.1} (i.e. within one order of
#'   magnitude). Lower values (e.g. \code{0.01}) retain more co-dominant
#'   candidates; higher values (e.g. \code{0.5}) are more stringent.
#' @param project_name Character. Prefix for output files and plot titles. Default: \code{"HiChIP"}.
#' @param color_palette Character. RColorBrewer palette name. Default: \code{"Set2"}.
#' @param karyo_bin_size Integer. Bin width in bp for karyotype heatmaps. Default: \code{1e5}.
#' @param neighbor_hop Integer. k-hop ego-network expansion order via \code{igraph::ego()}. \code{0} restricts to direct contacts. Default: \code{0}. Target gene assignment searches one additional hop (\code{neighbor_hop + 1}) to capture genes at the opposite anchor of directly connected loops.
#' @param anchor_gap Integer. Search radius: how far apart (bp) can a peak and loop anchor be for the peak to be considered "near" the anchor? \code{-1L} (default): strict physical overlap required (GenomicRanges default -- peak and anchor must share at least 1 bp). \code{0L}: adjacent intervals (peak end == anchor start) also count. \code{>0}: explicit gap tolerance (e.g. \code{200} for cross-experiment integration). When \code{>= 0}, the result includes both physically overlapping pairs AND proximity-only pairs (within gap but no actual overlap). Use \code{anchor_min_overlap > 1} to require actual physical overlap among these candidates.
#' @param anchor_min_overlap Integer. After candidate pairs are found (via \code{anchor_gap}), how many base pairs of actual physical overlap are required? Default \code{1L}: any touch counts (including proximity-only hits when \code{anchor_gap >= 0}). Increase to \code{10-100} to filter out spurious boundary overlaps. Setting this \code{> 1} with \code{anchor_gap >= 0} ensures that even with gap tolerance, only pairs with genuine physical overlap are retained.
#' @param anchor_min_frac Numeric (0-1). After the first two filters, what fraction of the \emph{peak} width must physically overlap the anchor? Default \code{0}: any fraction accepted. Set to \code{0.1-0.5} when peaks are broad (e.g. H3K27ac domains, 2-5 kb) so a 1 bp overlap does not link the entire broad peak. Ignored for point features (SNPs, eQTLs). Applied last, only to pairs that passed \code{anchor_gap} and \code{anchor_min_overlap}.
#' @param hub_percentile Numeric (0-1). Loop-count quantile for hub detection. Default: \code{0.95}. Genes or distal elements with connectivity at or above this quantile are flagged as hubs. A minimum floor of 3 (promoter-centric) or 2 (distal) is applied to avoid false hubs in small datasets.
#' @param write_output Logical. If \code{TRUE} (default), write the Excel workbook to \code{out_dir}. If \code{FALSE}, return results without creating directories or files.
#' @param quiet Logical. If \code{TRUE}, suppress progress messages while preserving warnings. Default: \code{FALSE}.
#'
#' @return An invisible named list:
#' \itemize{
#'   \item \code{target_annotation} -- Target features (peaks) with gene assignments.
#'     Key columns include:
#'     \itemize{
#'       \item \code{All_Loop_Connected_Genes}: All genes from loop-connected anchors (P/G types).
#'       \item \code{Regulated_promoter_genes}: Promoter genes supported by loop-anchor context.
#'       \item \code{Assigned_Target_Genes}: Promoter-first 3D assignment (prioritises P > G > E).
#'       \item \code{*_Filled} variants: Linear nearest-gene fallback when strict 3D assignments are empty.
#'       \item \code{Regulated_promoter_Evidence}: Provenance of \code{Regulated_promoter_genes}
#'         (e.g., \code{local_promoter_overlap}, \code{direct_opposite_promoter}).
#'         \strong{Read with} \code{Regulated_promoter_genes}; do not cross-reference
#'         with \code{Assigned_Target_Genes} or other columns.
#'       \item \code{Regulated_promoter_Fallback_Evidence}: Provenance of
#'         \code{Regulated_promoter_genes_Filled}.
#'         \strong{Read with} \code{Regulated_promoter_genes_Filled}; indicates
#'         which \code{*_Filled} column supplied the fallback gene.
#'     }
#'   \item \code{target_gene_links} -- Long-format peak-gene provenance table.
#'     Each row records one peak-gene linkage with full provenance.
#'     \strong{Read} \code{evidence}, \code{anchor_role}, and \code{gene_role}
#'     \strong{together as a group} -- they jointly describe how each gene was
#'     assigned to each peak; do not interpret any one column in isolation.
#'     \itemize{
#'       \item \code{input_id}, \code{loop_ID}, \code{anchor_id}: Identifiers.
#'       \item \code{gene}: Linked gene symbol.
#'       \item \code{gene_role}: \code{"promoter"}, \code{"gene_body"}, or \code{"linear_annotation"}.
#'       \item \code{source}: \code{"loop_anchor"} (3D-derived) or \code{"linear_annotation"} (nearest gene).
#'       \item \code{evidence}: Provenance label --
#'         \code{"local_promoter_overlap"} (peak overlaps anchor promoter),
#'         \code{"direct_opposite_promoter"} (opposite anchor is promoter),
#'         \code{"gene_body_context"} (gene body linkage),
#'         \code{"expanded_promoter_loop"} (via ego-network expansion),
#'         \code{"linear_annotation"} (direct nearest gene),
#'         or \code{"linear_fallback"} (filled when 3D assignment was empty).
#'       \item \code{anchor_role}: \code{"local_anchor"}, \code{"opposite_anchor"},
#'         \code{"expanded_anchor"}, or \code{"linear_annotation"}.
#'       \item \code{used_as_fallback}: Logical. \code{TRUE} when this link was added
#'         via the \code{*_Filled} linear nearest-gene fallback mechanism.
#'       \item \code{in_regulated_promoter} through \code{in_assigned_target_filled}:
#'         Logical membership flags indicating which target annotation column(s)
#'         this gene appears in.
#'     }
#'   \item \code{loop_annotation} -- Annotated 3D interactome with \code{Putative_Target_Genes}.
#'   \item \code{anchor_loci_annotation} -- Non-redundant anchor-locus genomic classifications after within-cluster interval reduction.
#'   \item \code{anchor_annotation} -- Backward-compatible alias of \code{anchor_loci_annotation}.
#'   \item \code{promoter_centric_stats} -- Gene-level connectivity statistics.
#'   \item \code{distal_element_stats} -- Distal-element connectivity statistics.
#'   \item \code{plots} -- Named list of ggplot/grob objects:
#'     \code{Basic_Donut}, \code{Basic_Circular}, \code{Basic_Flower},
#'     \code{Karyo_LoopGenes}, \code{Karyo_Anchors},
#'     \code{Anchor_Genomic_Distribution}, and (when \code{target_bed} is
#'     provided) \code{Karyo_TargetGenes}, \code{Target_Rose},
#'     \code{Target_Genomic_Distribution}, \code{Target_Loop_Genomic_Distribution}.
#'   \item \code{plot_list} -- Backward-compatible alias of \code{plots}.
#'   \item \code{metadata} -- Internal metadata list (parameters, versions). Not intended for direct use.
#' }
#' If \code{write_output = TRUE}, also writes a multi-sheet Excel workbook to \code{out_dir}.
#'
#' @export
#'
#' @examples
#' # Minimal runnable example for package checks
#' if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
#'     txdb_example <- AnnotationDbi::loadDb(
#'         system.file("extdata", "hg19_knownGene_sample.sqlite", package = "GenomicFeatures")
#'     )
#'     bedpe_path <- tempfile(fileext = ".bedpe")
#'     writeLines(
#'         "chr6\t10412000\t10412600\tchr6\t10415000\t10415600",
#'         bedpe_path
#'     )
#'
#'     res <- annotate_peaks_and_loops(
#'         bedpe_file = bedpe_path,
#'         txdb = txdb_example,
#'         org_db = "org.Hs.eg.db",
#'         species = "hg19",
#'         out_dir = tempdir(),
#'         project_name = "Quick_Example",
#'         write_output = FALSE,
#'         quiet = TRUE
#'     )
#'     head(res$loop_annotation, 1)
#' }

annotate_peaks_and_loops <- function(
  bedpe_file,
  target_bed = NULL,
  txdb = NULL,
  org_db = NULL,
  species = "hg38",
  tss_region = c(-2000, 2000),
  out_dir = "./results",
  expr_matrix_file = NULL,
  sample_columns = NULL,
  project_name = "HiChIP",
  color_palette = "Set2",
  karyo_bin_size = 1e5,
  neighbor_hop = 0,
  hub_percentile = 0.95,
  min_expr = 0,
  conflict_strategy = c("biotype_first", "expression_first"),
  co_dominance_ratio = 0.1,
  anchor_gap = -1L,
  anchor_min_overlap = 1L,
  anchor_min_frac = 0,
  write_output = TRUE,
  quiet = FALSE
) {
    species <- match.arg(species, c("hg38", "hg19", "mm10", "mm9"))
    stopifnot(length(tss_region) == 2L)
    conflict_strategy <- match.arg(conflict_strategy)
    .validate_annotation_params(anchor_gap, anchor_min_overlap, anchor_min_frac,
                                hub_percentile, neighbor_hop, karyo_bin_size)
    log_message <- function(...) {
        if (!quiet) message(...)
    }

    if (write_output && !dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

    tx_res <- .resolve_annotation_resource(txdb, "txdb", "TxDb", species)
    org_res <- .resolve_annotation_resource(org_db, "orgdb", "OrgDb", species)
    txdb_obj <- tx_res$obj
    org_db_pkg <- org_res$pkg
    if (is.null(org_db_pkg) || !nzchar(org_db_pkg)) {
        stop(
            "Unable to resolve OrgDb package name from supplied object. ",
            "Pass a package name string such as 'org.Hs.eg.db', or an installed OrgDb package object."
        )
    }

    gene_expr_map <- NULL
    if (!is.null(expr_matrix_file) && !is.null(sample_columns)) {
        log_message("Step 0: Loading expression data...")
        gene_expr_map <- load_expression_matrix(expr_matrix_file, sample_columns)
        log_message("    >>> Expression loaded for ", length(gene_expr_map), " genes.")
    } else if (!is.null(expr_matrix_file) && is.null(sample_columns)) {
        warning("expr_matrix_file provided but sample_columns is NULL; expression data will not be used for conflict resolution.", call. = FALSE)
    }

    # --- Step 1 & 2: Read BEDPE and cluster anchors ---
    bedpe_data <- .read_and_cluster_bedpe(bedpe_file, log_message, quiet)

    # --- Step 3: Biological Classification & Topology ---
    if (length(bedpe_data$gr_anchors) == 0) {
        warning("No valid loop anchors found; returning empty annotation.")
        return(list(
            anchor_loci_annotation = data.frame(), anchor_annotation = data.frame(),
            loop_annotation = data.frame(),
            promoter_centric_stats = data.frame(), distal_element_stats = data.frame(),
            target_annotation = NULL, target_gene_links = NULL,
            plots = list(), plot_list = list()
        ))
    }
    class_data <- .classify_anchors_and_topology(
        gr_anchors = bedpe_data$gr_anchors,
        loops = bedpe_data$loops,
        g = bedpe_data$g,
        txdb_obj = txdb_obj,
        org_db_pkg = org_db_pkg,
        tss_region = tss_region,
        gene_expr_map = gene_expr_map,
        min_expr = min_expr,
        conflict_strategy = conflict_strategy,
        neighbor_hop = neighbor_hop,
        log_message = log_message,
        co_dominance_ratio = co_dominance_ratio
    )

    # --- Step 4: Build annotation tables & stats ---
    tables <- .build_annotation_stats_tables(
        loops_annotated = class_data$loops_annotated,
        anchors = bedpe_data$anchors,
        map_info = class_data$map_info,
        anchor_topo_map = class_data$anchor_topo_map,
        cluster_regions = bedpe_data$cluster_regions,
        txdb_obj = txdb_obj,
        org_db_pkg = org_db_pkg,
        tss_region = tss_region,
        hub_percentile = hub_percentile,
        log_message = log_message
    )
    loop_annotation_final <- tables$loop_annotation_final
    cluster_info <- tables$cluster_info
    promoter_centric_df <- tables$promoter_centric_df
    distal_element_df <- tables$distal_element_df

    # --- Step 5: Integrate target annotations ---
    bed_info <- NULL
    target_connected_loops <- NULL
    target_gene_links <- NULL
    if (!is.null(target_bed)) {
        log_message("Step 5: Integrating Target Annotations...")
        target_res <- .annotate_target_bed(
            target_bed = target_bed,
            txdb_obj = txdb_obj,
            org_db_pkg = org_db_pkg,
            tss_region = tss_region,
            gene_expr_map = gene_expr_map,
            min_expr = min_expr,
            conflict_strategy = conflict_strategy,
            gr_anchors = bedpe_data$gr_anchors,
            anchor_topo_map = class_data$anchor_topo_map,
            loop_annotation_final = loop_annotation_final,
            map_info = class_data$map_info,
            ego_list_target = class_data$ego_list_target,
            log_message = log_message,
            anchor_gap = anchor_gap,
            anchor_min_overlap = anchor_min_overlap,
            anchor_min_frac = anchor_min_frac,
            co_dominance_ratio = co_dominance_ratio
        )
        bed_info <- target_res$bed_info
        target_connected_loops <- target_res$target_connected_loops
        target_gene_links <- target_res$target_gene_links
    }

    # --- Step 6: Visualization ---
    log_message("Step 6: Generating Visualizations (Returning plot objects)...")
    plot_list <- .generate_annotation_plots(
        loop_annotation_final, bed_info, cluster_info, target_connected_loops,
        txdb_obj, org_db_pkg, species, project_name, color_palette,
        karyo_bin_size, quiet
    )

    # --- Step 7: Export ---
    # Keep a1_id/a2_id in the R object for downstream functions
    # (chromatin refinement needs them).  Only strip for Excel export.
    loop_annotation_export <- loop_annotation_final %>%
        dplyr::select(-any_of(c("a1_id", "a2_id")))
    if (write_output) {
        log_message("Step 7: Exporting to Excel...")
        .export_annotation_excel(
            loop_annotation_clean = loop_annotation_export,
            cluster_info = cluster_info,
            promoter_centric_df = promoter_centric_df,
            distal_element_df = distal_element_df,
            bed_info = bed_info,
            target_gene_links = target_gene_links,
            out_dir = out_dir,
            project_name = project_name
        )
        log_message("    Excel file saved.")
    }

    log_message("Analysis Complete.")

    # Save intermediate anchor state for downstream chromatin-aware recomputation
    anchor_state <- if (!is.null(class_data$map_info)) {
        list(
            map_info         = class_data$map_info,
            anchor_topo_map  = class_data$anchor_topo_map,
            gr_anchors       = bedpe_data$gr_anchors,
            ego_list_target  = class_data$ego_list_target,
            target_hit_df    = if (exists("target_res", inherits = FALSE) && !is.null(target_res)) target_res$hit_df else NULL
        )
    } else NULL

    result <- list(
        anchor_loci_annotation = cluster_info,
        anchor_annotation = cluster_info,
        loop_annotation = loop_annotation_final,
        promoter_centric_stats = promoter_centric_df,
        distal_element_stats = distal_element_df,
        target_annotation = bed_info,
        target_gene_links = target_gene_links,
        plots = plot_list,
        plot_list = plot_list,
        metadata = .build_looplook_metadata(
            fun = "annotate_peaks_and_loops",
            params = list(
                species = species,
                tss_region = tss_region,
                neighbor_hop = neighbor_hop,
                hub_percentile = hub_percentile,
                min_expr = min_expr,
                conflict_strategy = conflict_strategy,
                co_dominance_ratio = co_dominance_ratio,
                has_target_bed = !is.null(target_bed),
                has_expression = !is.null(expr_matrix_file)
            ),
            genome_build = species,
            score_semantics = "loop type classification; higher n_members = more evidence",
            database_versions = .record_database_versions(species)
        )
    )
    attr(result, "looplook_anchor_state") <- anchor_state
    return(result)
}

# --- Internal helpers extracted from annotate_peaks_and_loops ---

#' Internal: Validate annotation parameters
#' @keywords internal
#' @noRd
.validate_annotation_params <- function(anchor_gap, anchor_min_overlap,
                                         anchor_min_frac, hub_percentile,
                                         neighbor_hop, karyo_bin_size) {
    if (!is.numeric(anchor_gap) || length(anchor_gap) != 1L ||
        is.na(anchor_gap) || anchor_gap < -1L)
        stop("anchor_gap must be >= -1 (-1 = strict overlap, 0 = adjacent, >0 = proximity)", call. = FALSE)
    if (!is.numeric(anchor_min_overlap) || length(anchor_min_overlap) != 1L ||
        is.na(anchor_min_overlap) || anchor_min_overlap < 1L)
        stop("anchor_min_overlap must be a positive integer (>= 1)", call. = FALSE)
    if (!is.numeric(anchor_min_frac) || length(anchor_min_frac) != 1L ||
        is.na(anchor_min_frac) || anchor_min_frac < 0 || anchor_min_frac > 1)
        stop("anchor_min_frac must be in [0, 1]", call. = FALSE)
    if (!is.numeric(hub_percentile) || length(hub_percentile) != 1L ||
        is.na(hub_percentile) || hub_percentile <= 0 || hub_percentile > 1)
        stop("hub_percentile must be in (0, 1]", call. = FALSE)
    if (!is.numeric(neighbor_hop) || length(neighbor_hop) != 1L ||
        is.na(neighbor_hop) || neighbor_hop < 0 || neighbor_hop != floor(neighbor_hop))
        stop("neighbor_hop must be a non-negative integer", call. = FALSE)
    if (!is.numeric(karyo_bin_size) || length(karyo_bin_size) != 1L ||
        is.na(karyo_bin_size) || karyo_bin_size < 1)
        stop("karyo_bin_size must be a positive number", call. = FALSE)
}

#' Internal: Generate annotation plots with optional message silencing
#' @keywords internal
#' @noRd
.generate_annotation_plots <- function(loop_annotation_final, bed_info,
    cluster_info, target_connected_loops, txdb_obj, org_db_pkg, species,
    project_name, color_palette, karyo_bin_size, quiet) {
    plot_df <- loop_annotation_final
    plot_df$loop_genes <- plot_df$All_Anchor_Genes
    if (quiet) {
        .with_messages_silenced(
            .with_known_upstream_noise_suppressed(
                build_annotation_plots(
                    plot_df = plot_df, bed_info = bed_info,
                    cluster_info = cluster_info,
                    target_connected_loops = target_connected_loops,
                    txdb_obj = txdb_obj, org_db_pkg = org_db_pkg,
                    species = species, project_name = project_name,
                    color_palette = color_palette, karyo_bin_size = karyo_bin_size
                )
            )
        )
    } else {
        .with_known_upstream_noise_suppressed(
            build_annotation_plots(
                plot_df = plot_df, bed_info = bed_info,
                cluster_info = cluster_info,
                target_connected_loops = target_connected_loops,
                txdb_obj = txdb_obj, org_db_pkg = org_db_pkg,
                species = species, project_name = project_name,
                color_palette = color_palette, karyo_bin_size = karyo_bin_size
            )
        )
    }
}

#' Internal: Read BEDPE and cluster anchors
#' @return A list with \code{loops}, \code{anchors}, \code{gr_anchors},
#'   \code{cluster_regions}, and \code{g} (igraph object).
#' @keywords internal
#' @noRd
.read_and_cluster_bedpe <- function(bedpe_file, log_message, quiet) {
    log_message("Step 1: Reading BEDPE file...")
    loops <- data.table::fread(bedpe_file, header = FALSE)
    ncol_loops <- ncol(loops)
    colnames(loops)[seq_len(6)] <- c("chr1", "start1", "end1", "chr2", "start2", "end2")
    if (ncol_loops > 6 && !quiet) {
        message(
            "    BEDPE contains ", ncol_loops, " columns; columns beyond 6 ",
            "are retained as-is in the output. ",
            "If column 7+ represent strand or per-anchor metadata, ",
            "note that anchor swapping only affects columns 1-6."
        )
    }
    loops <- .validate_bedpe_df(loops, quiet = quiet)
    anchors <- data.table::rbindlist(list(
        loops %>% dplyr::select(chr = chr1, start = start1, end = end1),
        loops %>% dplyr::select(chr = chr2, start = start2, end = end2)
    )) %>%
        unique() %>%
        dplyr::mutate(anchor_id = paste0("A", seq_len(dplyr::n())))
    loops <- loops %>%
        dplyr::left_join(
            anchors %>% dplyr::select(chr, start, end, a1_id = anchor_id),
            by = c("chr1" = "chr", "start1" = "start", "end1" = "end")
        ) %>%
        dplyr::left_join(
            anchors %>% dplyr::select(chr, start, end, a2_id = anchor_id),
            by = c("chr2" = "chr", "start2" = "start", "end2" = "end")
        )

    log_message("Step 2: Clustering loops...")
    valid_loops <- loops %>% dplyr::filter(!is.na(a1_id) & !is.na(a2_id))
    g <- igraph::graph_from_data_frame(
        valid_loops[, c("a1_id", "a2_id")], directed = FALSE
    )
    comp <- igraph::components(g)
    anchors$cluster_id <- NA
    comm <- intersect(anchors$anchor_id, names(comp$membership))
    if (length(comm) > 0)
        anchors$cluster_id[match(comm, anchors$anchor_id)] <-
            as.character(comp$membership[comm])
    anchors <- anchors %>% dplyr::filter(!is.na(cluster_id))
    loops <- loops %>%
        dplyr::left_join(
            anchors %>% dplyr::select(anchor_id, cluster_id),
            by = c("a1_id" = "anchor_id")
        )
    gr_anchors <- .with_known_upstream_noise_suppressed({
        gr_anchors <- GenomicRanges::makeGRangesFromDataFrame(
            anchors, keep.extra.columns = TRUE
        )
        gr_anchors$anchor_id <- anchors$anchor_id
        gr_anchors
    })
    cluster_regions <- .with_known_upstream_noise_suppressed({
        gr_list <- GenomicRanges::GRangesList(
            split(gr_anchors, gr_anchors$cluster_id)
        )
        cluster_regions <- unlist(GenomicRanges::reduce(gr_list))
        cluster_regions$cluster_id <- names(cluster_regions)
        names(cluster_regions) <- paste0("peak_", seq_along(cluster_regions))
        cluster_regions
    })

    list(
        loops = loops, anchors = anchors,
        gr_anchors = gr_anchors, cluster_regions = cluster_regions, g = g
    )
}

#' Internal: Classify anchors and compute network topology
#' @return A list with \code{map_info}, \code{loops_annotated},
#'   \code{anchor_topo_map}, and \code{ego_list_target}.
#' @keywords internal
#' @noRd
.classify_anchors_and_topology <- function(
    gr_anchors, loops, g, txdb_obj, org_db_pkg, tss_region,
    gene_expr_map, min_expr, conflict_strategy, neighbor_hop, log_message,
    co_dominance_ratio = 0.1
) {
    log_message("Step 3: Biological Classification & Topology...")
    anchor_anno <- .with_known_upstream_noise_suppressed(
        ChIPseeker::annotatePeak(
            gr_anchors, TxDb = txdb_obj, tssRegion = tss_region,
            annoDb = org_db_pkg, verbose = FALSE
        )
    )
    anchor_anno_df <- format_annotation_columns(as.data.frame(anchor_anno))
    anchor_anno_df <- resolve_gene_conflicts(
        anchor_anno_df, txdb_obj, org_db_pkg, tss_region,
        gene_expr_map, min_expr = min_expr,
        conflict_strategy = conflict_strategy,
        co_dominance_ratio = co_dominance_ratio
    )
    anchor_anno_df$type_code <- vapply(
        anchor_anno_df$annotation, .annotation_feature_class,
        FUN.VALUE = character(1)
    )
    map_info <- anchor_anno_df %>%
        dplyr::select(anchor_id, type_code, SYMBOL)
    loops_annotated <- loops %>%
        dplyr::left_join(
            map_info %>% dplyr::rename(t1 = type_code, s1 = SYMBOL),
            by = c("a1_id" = "anchor_id")
        ) %>%
        dplyr::left_join(
            map_info %>% dplyr::rename(t2 = type_code, s2 = SYMBOL),
            by = c("a2_id" = "anchor_id")
        )

    loops_annotated$loop_type <- unlist(
        Map(.loop_type_code, loops_annotated$t1, loops_annotated$t2),
        use.names = FALSE
    )
    locus_genes <- unlist(
        Map(.loop_locus_genes, loops_annotated$t1, loops_annotated$t2,
            loops_annotated$s1, loops_annotated$s2),
        use.names = FALSE
    )
    loops_annotated$single_loop_genes <- locus_genes
    loops_annotated$reg_loop_genes <- locus_genes

    log_message("    Calculating Topology (Hops, neighbor_hop = ", neighbor_hop, ")...")
    map_info$SYMBOL <- trimws(map_info$SYMBOL)
    valid_pg_nodes <- map_info %>%
        dplyr::filter((.is_promoter_like(type_code) | .is_gene_body_like(type_code)) & !is.na(SYMBOL) & SYMBOL != "")
    lookup_pg_symbol <- valid_pg_nodes$SYMBOL
    names(lookup_pg_symbol) <- valid_pg_nodes$anchor_id
    lookup_pg_type <- valid_pg_nodes$type_code
    names(lookup_pg_type) <- valid_pg_nodes$anchor_id
    lookup_p_symbol <- map_info %>%
        dplyr::filter(.is_promoter_like(type_code) & !is.na(SYMBOL) & SYMBOL != "") %>%
        dplyr::pull(SYMBOL)
    names(lookup_p_symbol) <- map_info %>%
        dplyr::filter(.is_promoter_like(type_code) & !is.na(SYMBOL) & SYMBOL != "") %>%
        dplyr::pull(anchor_id)
    nodes_in_graph <- igraph::V(g)$name
    input_hop <- if (is.null(neighbor_hop)) 0 else neighbor_hop
    ego_list_loop <- igraph::ego(
        g, order = input_hop, nodes = nodes_in_graph, mode = "all"
    )
    names(ego_list_loop) <- nodes_in_graph
    ego_list_target <- igraph::ego(
        g, order = input_hop + 1, nodes = nodes_in_graph, mode = "all"
    )
    names(ego_list_target) <- nodes_in_graph
    anchor_topo_map <- data.frame(
        anchor_id = nodes_in_graph,
        topo_genes_p = vapply(ego_list_loop, function(x)
            .ids_to_genes_simple(names(x), lookup_p_symbol), character(1)),
        topo_genes_pg = vapply(ego_list_loop, function(x)
            .ids_to_genes_simple(names(x), lookup_pg_symbol), character(1)),
        tgt_genes_pg = vapply(ego_list_target, function(x)
            .ids_to_genes_simple(names(x), lookup_pg_symbol), character(1)),
        tgt_genes_p = vapply(ego_list_target, function(x)
            .ids_to_genes_simple(names(x), lookup_p_symbol), character(1)),
        tgt_genes_prio = vapply(ego_list_target, function(x)
            .ids_to_genes_priority(names(x), lookup_pg_symbol, lookup_pg_type),
            character(1)),
        stringsAsFactors = FALSE
    )
    anchor_topo_map[is.na(anchor_topo_map)] <- NA_character_

    list(
        map_info = map_info,
        loops_annotated = loops_annotated,
        anchor_topo_map = anchor_topo_map,
        ego_list_target = ego_list_target
    )
}

#' Internal: Build loop annotation tables and compute connectivity stats
#' @return A list with \code{loop_annotation_final}, \code{cluster_info},
#'   \code{promoter_centric_df}, and \code{distal_element_df}.
#' @keywords internal
#' @noRd
.build_annotation_stats_tables <- function(
    loops_annotated, anchors, map_info, anchor_topo_map,
    cluster_regions, txdb_obj, org_db_pkg, tss_region,
    hub_percentile, log_message
) {
    log_message("Step 4: Constructing Loop Tables...")
    loops_annotated <- loops_annotated %>%
        dplyr::left_join(
            anchor_topo_map %>% dplyr::select(anchor_id, pg1 = topo_genes_pg),
            by = c("a1_id" = "anchor_id")
        ) %>%
        dplyr::left_join(
            anchor_topo_map %>% dplyr::select(anchor_id, pg2 = topo_genes_pg),
            by = c("a2_id" = "anchor_id")
        ) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(proximate_loop_gene = dplyr::case_when(
            (!is.na(t1) & t1 == "G" & !is.na(t2) & t2 == "P") ~
                extract_genes(pg2),
            (!is.na(t1) & t1 == "P" & !is.na(t2) & t2 == "G") ~
                extract_genes(pg1),
            TRUE ~ extract_genes(c(pg1, pg2))
        )) %>%
        dplyr::ungroup()
    clust_vec <- setNames(anchors$cluster_id, anchors$anchor_id)
    loops_annotated$cluster_id <- clust_vec[loops_annotated$a1_id]
    agg_cluster_reg <- loops_annotated %>%
        dplyr::filter(!is.na(cluster_id)) %>%
        dplyr::group_by(cluster_id) %>%
        dplyr::summarise(
            all_cluster_loop_genes = extract_genes(reg_loop_genes),
            .groups = "drop"
        )
    loop_annotation_final <- loops_annotated %>%
        dplyr::left_join(agg_cluster_reg, by = "cluster_id") %>%
        dplyr::mutate(loop_ID = paste0("L", seq_len(dplyr::n()))) %>%
        dplyr::select(
            loop_ID, chr1, start1, end1, chr2, start2, end2,
            cluster_id, loop_type,
            anchor1_gene = s1, anchor1_type = t1,
            anchor2_gene = s2, anchor2_type = t2,
            all_cluster_loop_genes, single_loop_genes, proximate_loop_gene,
            a1_id, a2_id
        ) %>%
        dplyr::rename(
            All_Anchor_Genes = single_loop_genes,
            Putative_Target_Genes = proximate_loop_gene,
            Cluster_All_Genes = all_cluster_loop_genes
        )
    agg_cluster_locus <- loop_annotation_final %>%
        dplyr::filter(!is.na(cluster_id)) %>%
        dplyr::group_by(cluster_id) %>%
        dplyr::summarise(
            Cluster_Locus_Genes = extract_genes(All_Anchor_Genes),
            .groups = "drop"
        )
    if (length(cluster_regions) == 0) {
        warning("No cluster regions formed from loop anchors.")
        cluster_info <- data.frame()
    } else {
        gene_annot <- .with_known_upstream_noise_suppressed(
            ChIPseeker::annotatePeak(
                cluster_regions, TxDb = txdb_obj,
                tssRegion = tss_region, annoDb = org_db_pkg, verbose = FALSE
            )
        )
        cluster_info <- format_annotation_columns(as.data.frame(gene_annot))
    }
    if ("GENENAME" %in% colnames(cluster_info))
        cluster_info <- cluster_info %>% dplyr::rename(Gene_description = GENENAME)
    cluster_info$cluster_id <- as.character(cluster_info$cluster_id)
    cluster_info <- cluster_info %>%
        dplyr::left_join(agg_cluster_locus, by = "cluster_id")

    # --- Promoter-centric stats ---
    log_message("    Generating Promoter Centric Stats...")
    raw_stats_df <- dplyr::bind_rows(
        loop_annotation_final %>%
            dplyr::filter(.is_promoter_like(anchor1_type) & !is.na(anchor1_gene)) %>%
            dplyr::select(
                Gene = anchor1_gene, Neighbor_Type = anchor2_type,
                Loop_Type = loop_type
            ) %>%
            dplyr::mutate(Gene = as.character(Gene)),
        loop_annotation_final %>%
            dplyr::filter(.is_promoter_like(anchor2_type) & !is.na(anchor2_gene)) %>%
            dplyr::select(
                Gene = anchor2_gene, Neighbor_Type = anchor1_type,
                Loop_Type = loop_type
            ) %>%
            dplyr::mutate(Gene = as.character(Gene))
    ) %>%
        tidyr::separate_rows(Gene, sep = ";") %>%
        dplyr::mutate(Gene = trimws(Gene)) %>%
        dplyr::filter(Gene != "" & !is.na(Gene)) %>%
        dplyr::group_by(Gene) %>%
        dplyr::summarise(
            Total_Loops = dplyr::n(),
            n_Linked_Promoters = sum(.is_promoter_like(Neighbor_Type), na.rm = TRUE),
            n_Linked_Distal = sum((.is_distal_like(Neighbor_Type) | .is_gene_body_like(Neighbor_Type)), na.rm = TRUE),
            Dominant_Interaction = names(which.max(table(Loop_Type))),
            .groups = "drop"
        )

    if (nrow(raw_stats_df) == 0) {
        promoter_centric_df <- data.frame(
            Gene = character(), Total_Loops = integer(),
            n_Linked_Promoters = integer(), n_Linked_Distal = integer(),
            Dominant_Interaction = character(),
            Is_High_Connectivity_Gene = character(),
            Is_High_Distal_Connectivity_Gene = character(),
            stringsAsFactors = FALSE
        )
    } else {
        final_cutoff <- max(
            quantile(raw_stats_df$Total_Loops, hub_percentile, na.rm = TRUE), 3
        )
        distal_cutoff <- max(
            quantile(raw_stats_df$n_Linked_Distal, hub_percentile, na.rm = TRUE), 2
        )
        promoter_centric_df <- raw_stats_df %>%
            dplyr::mutate(
                Is_High_Connectivity_Gene = dplyr::if_else(
                    Total_Loops >= final_cutoff, "Yes", "No"
                ),
                Is_High_Distal_Connectivity_Gene = dplyr::if_else(
                    n_Linked_Distal >= distal_cutoff, "Yes", "No"
                )
            ) %>%
            dplyr::arrange(dplyr::desc(n_Linked_Distal))
    }

    # --- Distal element stats ---
    log_message("    Generating Distal Element Stats...")
    distal_element_df <- .build_distal_element_stats(
        loop_annotation_final, anchors, hub_percentile
    )

    list(
        loop_annotation_final = loop_annotation_final,
        cluster_info = cluster_info,
        promoter_centric_df = promoter_centric_df,
        distal_element_df = distal_element_df
    )
}

#' Internal: Export annotation results to multi-sheet Excel workbook
#' @keywords internal
#' @noRd
.export_annotation_excel <- function(
    loop_annotation_clean, cluster_info, promoter_centric_df,
    distal_element_df, bed_info, target_gene_links, out_dir, project_name
) {
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, "Loop Annotation")
    openxlsx::writeData(wb, "Loop Annotation", loop_annotation_clean)
    openxlsx::addWorksheet(wb, "Anchor Loci Annotation")
    openxlsx::writeData(wb, "Anchor Loci Annotation", cluster_info)
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
    if (!is.null(target_gene_links) && nrow(target_gene_links) > 0) {
        openxlsx::addWorksheet(wb, "Target Gene Links")
        openxlsx::writeData(wb, "Target Gene Links", target_gene_links)
    }
    tryCatch(
        openxlsx::saveWorkbook(
            wb,
            file.path(out_dir, paste0(project_name, "_Basic_Results.xlsx")),
            overwrite = TRUE
        ),
        error = function(e)
            warning("Failed to save Excel workbook: ", conditionMessage(e),
                    call. = FALSE)
    )
}

#' Internal: Build distal element connectivity data frame
#' @keywords internal
#' @noRd
.build_distal_element_stats <- function(loop_annotation_final, anchors,
                                          hub_percentile) {
    distal_raw_df <- dplyr::bind_rows(
        loop_annotation_final %>%
            dplyr::filter(anchor1_type %in% c("E", "G")) %>%
            dplyr::select(Distal_Anchor_ID = a1_id, Distal_Type = anchor1_type,
                Neighbor_Gene = anchor2_gene, Neighbor_Type = anchor2_type,
                Loop_Type = loop_type),
        loop_annotation_final %>%
            dplyr::filter(anchor2_type %in% c("E", "G")) %>%
            dplyr::select(Distal_Anchor_ID = a2_id, Distal_Type = anchor2_type,
                Neighbor_Gene = anchor1_gene, Neighbor_Type = anchor1_type,
                Loop_Type = loop_type)
    ) %>%
        dplyr::group_by(Distal_Anchor_ID) %>%
        dplyr::summarise(
            Total_Loops = dplyr::n(),
            n_Linked_Promoters = sum(.is_promoter_like(Neighbor_Type), na.rm = TRUE),
            n_Linked_Distal = sum((.is_distal_like(Neighbor_Type) | .is_gene_body_like(Neighbor_Type)), na.rm = TRUE),
            Dominant_Interaction = names(which.max(table(Loop_Type))),
            Target_Genes = extract_genes(Neighbor_Gene[.is_promoter_like(Neighbor_Type)]),
            .groups = "drop"
        )
    anchor_coords_map <- anchors %>%
        dplyr::select(anchor_id, chr, start, end, cluster_id) %>%
        dplyr::distinct()
    if (nrow(distal_raw_df) == 0) {
        return(data.frame(
            chr = character(), start = integer(), end = integer(),
            cluster_id = character(), Total_Loops = integer(),
            n_Linked_Promoters = integer(), n_Linked_Distal = integer(),
            Dominant_Interaction = character(),
            Is_High_Connectivity_Distal_Element = character(),
            Target_Genes = character(), stringsAsFactors = FALSE
        ))
    }
    final_cutoff <- max(quantile(distal_raw_df$Total_Loops, hub_percentile,
        na.rm = TRUE), 3)
    distal_raw_df %>%
        dplyr::left_join(anchor_coords_map,
            by = c("Distal_Anchor_ID" = "anchor_id")) %>%
        dplyr::mutate(Is_High_Connectivity_Distal_Element =
            dplyr::if_else(Total_Loops >= final_cutoff, "Yes", "No")) %>%
        dplyr::select(chr, start, end, cluster_id, Total_Loops,
            n_Linked_Promoters, n_Linked_Distal, Dominant_Interaction,
            Is_High_Connectivity_Distal_Element, Target_Genes) %>%
        dplyr::arrange(dplyr::desc(Total_Loops))
}

.empty_target_gene_links <- function() {
    data.frame(
        input_id = character(),
        loop_ID = character(),
        anchor_id = character(),
        gene = character(),
        gene_role = character(),
        source = character(),
        evidence = character(),
        anchor_role = character(),
        used_as_fallback = logical(),
        in_regulated_promoter = logical(),
        in_assigned_target = logical(),
        in_all_loop_connected = logical(),
        in_regulated_promoter_filled = logical(),
        in_assigned_target_filled = logical(),
        stringsAsFactors = FALSE
    )
}

.target_gene_link_flags <- function(n) {
    data.frame(
        used_as_fallback = rep(FALSE, n),
        in_regulated_promoter = rep(FALSE, n),
        in_assigned_target = rep(FALSE, n),
        in_all_loop_connected = rep(FALSE, n),
        in_regulated_promoter_filled = rep(FALSE, n),
        in_assigned_target_filled = rep(FALSE, n),
        stringsAsFactors = FALSE
    )
}

.target_linear_gene_column <- function(bed_info) {
    if ("SYMBOL" %in% colnames(bed_info)) {
        return("SYMBOL")
    }
    if ("geneId" %in% colnames(bed_info)) {
        return("geneId")
    }
    NULL
}

.collapse_target_values <- function(x, default = NA_character_) {
    x <- unique(trimws(as.character(na.omit(x))))
    x <- x[nzchar(x)]
    if (length(x) == 0) {
        return(default)
    }
    paste(sort(x), collapse = ";")
}

.target_anchor_gene_map <- function(map_info) {
    if (is.null(map_info) || nrow(map_info) == 0 ||
        !all(c("anchor_id", "type_code", "SYMBOL") %in% colnames(map_info))) {
        return(data.frame(
            anchor_id = character(), gene = character(),
            gene_role = character(), stringsAsFactors = FALSE
        ))
    }

    map_info %>%
        dplyr::filter((.is_promoter_like(type_code) | .is_gene_body_like(type_code)), !is.na(SYMBOL), SYMBOL != "") %>%
        dplyr::select(anchor_id, type_code, SYMBOL) %>%
        dplyr::mutate(SYMBOL = as.character(SYMBOL)) %>%
        tidyr::separate_rows(SYMBOL, sep = ";") %>%
        dplyr::mutate(gene = trimws(SYMBOL)) %>%
        dplyr::filter(!is.na(gene), gene != "") %>%
        dplyr::transmute(
            anchor_id = as.character(anchor_id),
            gene = gene,
            gene_role = dplyr::if_else(.is_promoter_like(type_code), "promoter", "gene_body")
        ) %>%
        dplyr::distinct()
}

.linear_target_gene_links <- function(bed_info) {
    empty <- .empty_target_gene_links()[, c(
        "input_id", "loop_ID", "anchor_id", "gene", "gene_role",
        "source", "evidence", "anchor_role"
    )]
    linear_col <- .target_linear_gene_column(bed_info)
    if (is.null(linear_col) || is.null(bed_info) || nrow(bed_info) == 0) {
        return(empty)
    }

    input_ids <- if ("input_id" %in% colnames(bed_info)) {
        as.character(bed_info$input_id)
    } else {
        paste0("Peak_", seq_len(nrow(bed_info)))
    }

    rows <- lapply(seq_len(nrow(bed_info)), function(i) {
        genes <- clean_gene_names(bed_info[[linear_col]][i], ";")
        if (length(genes) == 0) {
            return(NULL)
        }
        data.frame(
            input_id = input_ids[[i]],
            loop_ID = NA_character_,
            anchor_id = NA_character_,
            gene = genes,
            gene_role = "linear_annotation",
            source = "linear_annotation",
            evidence = "linear_annotation",
            anchor_role = "linear_annotation",
            stringsAsFactors = FALSE
        )
    })

    out <- do.call(rbind, Filter(Negate(is.null), rows))
    if (is.null(out)) {
        return(empty)
    }
    out
}

.build_target_gene_links <- function(
  hit_df, bed_info, loop_annotation_final, map_info, ego_list_target
) {
    base_cols <- c(
        "input_id", "loop_ID", "anchor_id", "gene", "gene_role",
        "source", "evidence", "anchor_role"
    )
    empty <- .empty_target_gene_links()[, base_cols]
    gene_map <- .target_anchor_gene_map(map_info)
    linear_rows <- .linear_target_gene_links(bed_info)

    if (is.null(hit_df) || nrow(hit_df) == 0 || nrow(gene_map) == 0) {
        rows <- dplyr::bind_rows(linear_rows)
        if (nrow(rows) == 0) {
            return(.empty_target_gene_links())
        }
        rows <- dplyr::distinct(rows)
        return(cbind(rows, .target_gene_link_flags(nrow(rows))))
    }

    hit_base <- hit_df %>%
        dplyr::filter(!is.na(anchor_id), anchor_id != "") %>%
        dplyr::transmute(
            input_id = paste0("Peak_", qid),
            local_anchor_id = as.character(anchor_id)
        ) %>%
        dplyr::distinct()

    edge_long <- dplyr::bind_rows(
        loop_annotation_final %>% dplyr::select(loop_ID, local_anchor_id = a1_id, opposite_anchor_id = a2_id),
        loop_annotation_final %>% dplyr::select(loop_ID, local_anchor_id = a2_id, opposite_anchor_id = a1_id)
    ) %>%
        dplyr::filter(!is.na(local_anchor_id), !is.na(opposite_anchor_id)) %>%
        dplyr::mutate(
            local_anchor_id = as.character(local_anchor_id),
            opposite_anchor_id = as.character(opposite_anchor_id)
        ) %>%
        dplyr::distinct()

    local_seed <- hit_base %>%
        dplyr::inner_join(edge_long %>% dplyr::select(local_anchor_id, loop_ID) %>% dplyr::distinct(),
            by = "local_anchor_id"
        )
    if (nrow(local_seed) == 0) {
        local_seed <- hit_base %>% dplyr::mutate(loop_ID = NA_character_)
    }

    local_rows <- local_seed %>%
        dplyr::left_join(gene_map, by = c("local_anchor_id" = "anchor_id")) %>%
        dplyr::filter(!is.na(gene), gene != "") %>%
        dplyr::transmute(
            input_id,
            loop_ID,
            anchor_id = local_anchor_id,
            gene,
            gene_role,
            source = "loop_anchor",
            evidence = dplyr::if_else(gene_role == "promoter", "local_promoter_overlap", "gene_body_context"),
            anchor_role = "local_anchor"
        )

    direct_rows <- hit_base %>%
        dplyr::inner_join(edge_long, by = "local_anchor_id") %>%
        dplyr::left_join(gene_map, by = c("opposite_anchor_id" = "anchor_id")) %>%
        dplyr::filter(!is.na(gene), gene != "") %>%
        dplyr::transmute(
            input_id,
            loop_ID,
            anchor_id = opposite_anchor_id,
            gene,
            gene_role,
            source = "loop_anchor",
            evidence = dplyr::if_else(gene_role == "promoter", "direct_opposite_promoter", "gene_body_context"),
            anchor_role = "opposite_anchor"
        )

    expanded_seed <- do.call(rbind, lapply(seq_len(nrow(hit_base)), function(i) {
        local_id <- hit_base$local_anchor_id[[i]]
        ego_ids <- if (local_id %in% names(ego_list_target)) {
            names(ego_list_target[[local_id]])
        } else {
            character(0)
        }
        direct_ids <- edge_long$opposite_anchor_id[edge_long$local_anchor_id == local_id]
        expanded_ids <- setdiff(ego_ids, c(local_id, direct_ids))
        if (length(expanded_ids) == 0) {
            return(NULL)
        }
        data.frame(
            input_id = hit_base$input_id[[i]],
            anchor_id = expanded_ids,
            stringsAsFactors = FALSE
        )
    }))
    expanded_rows <- empty
    if (!is.null(expanded_seed) && nrow(expanded_seed) > 0) {
        expanded_rows <- expanded_seed %>%
            dplyr::left_join(gene_map, by = "anchor_id") %>%
            dplyr::filter(!is.na(gene), gene != "") %>%
            dplyr::transmute(
                input_id,
                loop_ID = NA_character_,
                anchor_id,
                gene,
                gene_role,
                source = "loop_anchor",
                evidence = dplyr::if_else(gene_role == "promoter", "expanded_promoter_loop", "gene_body_context"),
                anchor_role = "expanded_anchor"
            )
    }

    rows <- dplyr::bind_rows(local_rows, direct_rows, expanded_rows, linear_rows) %>%
        dplyr::distinct()
    if (nrow(rows) == 0) {
        return(.empty_target_gene_links())
    }

    cbind(rows, .target_gene_link_flags(nrow(rows)))
}

.summarise_regulated_promoter_evidence <- function(target_gene_links) {
    if (is.null(target_gene_links) || nrow(target_gene_links) == 0) {
        return(data.frame(input_id = character(), Regulated_promoter_Evidence = character()))
    }
    target_gene_links %>%
        dplyr::filter(source == "loop_anchor", gene_role == "promoter") %>%
        dplyr::group_by(input_id) %>%
        dplyr::summarise(
            Regulated_promoter_Evidence = .collapse_target_values(evidence, default = "none"),
            .groups = "drop"
        )
}

.fallback_evidence_from_annotation <- function(annotation) {
    vapply(annotation, function(x) {
        if (is.na(x) || x == "") {
            return("linear_fallback")
        }
        x <- tolower(as.character(x))
        if (grepl("promoter", x)) {
            return("local_promoter")
        }
        if (grepl("exon|intron|utr", x)) {
            return("local_gene_body")
        }
        "linear_nearest"
    }, FUN.VALUE = character(1))
}

.contains_target_gene <- function(gene, gene_string) {
    if (is.na(gene) || gene == "" || is.na(gene_string) || gene_string == "") {
        return(FALSE)
    }
    gene %in% clean_gene_names(gene_string, ";")
}

.mark_target_gene_link_membership <- function(target_gene_links, bed_info) {
    if (is.null(target_gene_links) || nrow(target_gene_links) == 0) {
        return(.empty_target_gene_links())
    }

    target_cols <- c(
        "Regulated_promoter_genes", "Assigned_Target_Genes",
        "All_Loop_Connected_Genes", "Regulated_promoter_genes_Filled",
        "Assigned_Target_Genes_Filled"
    )
    lookup_cols <- intersect(c("input_id", target_cols), colnames(bed_info))
    link_df <- target_gene_links %>%
        dplyr::select(-dplyr::any_of(c(
            "used_as_fallback", "in_regulated_promoter", "in_assigned_target",
            "in_all_loop_connected", "in_regulated_promoter_filled",
            "in_assigned_target_filled"
        ))) %>%
        dplyr::left_join(bed_info[, lookup_cols, drop = FALSE], by = "input_id")

    for (col in target_cols) {
        if (!col %in% colnames(link_df)) {
            link_df[[col]] <- NA_character_
        }
    }

    link_df$in_regulated_promoter <- unlist(Map(.contains_target_gene, link_df$gene, link_df$Regulated_promoter_genes), use.names = FALSE)
    link_df$in_assigned_target <- unlist(Map(.contains_target_gene, link_df$gene, link_df$Assigned_Target_Genes), use.names = FALSE)
    link_df$in_all_loop_connected <- unlist(Map(.contains_target_gene, link_df$gene, link_df$All_Loop_Connected_Genes), use.names = FALSE)
    link_df$in_regulated_promoter_filled <- unlist(Map(.contains_target_gene, link_df$gene, link_df$Regulated_promoter_genes_Filled), use.names = FALSE)
    link_df$in_assigned_target_filled <- unlist(Map(.contains_target_gene, link_df$gene, link_df$Assigned_Target_Genes_Filled), use.names = FALSE)

    link_df$evidence[link_df$source == "linear_annotation"] <- "linear_annotation"
    reg_empty <- is.na(link_df$Regulated_promoter_genes) | link_df$Regulated_promoter_genes == ""
    assigned_empty <- is.na(link_df$Assigned_Target_Genes) | link_df$Assigned_Target_Genes == ""
    link_df$used_as_fallback <- link_df$source == "linear_annotation" &
        ((reg_empty & link_df$in_regulated_promoter_filled) |
            (assigned_empty & link_df$in_assigned_target_filled))
    link_df$evidence[link_df$used_as_fallback] <- "linear_fallback"

    link_df %>%
        dplyr::select(
            dplyr::any_of(c(
                "input_id", "loop_ID", "anchor_id", "gene", "gene_role",
                "source", "evidence", "anchor_role",
                "anchor_type_before_chromatin", "anchor_type_after_chromatin",
                "chromatin_confidence", "chromatin_evidence",
                "chromatin_target_action",
                "used_as_fallback", "in_regulated_promoter",
                "in_assigned_target", "in_all_loop_connected",
                "in_regulated_promoter_filled", "in_assigned_target_filled"
            ))
        ) %>%
        dplyr::distinct()
}

#' Internal: Filter Target Gene Links After Expression-Aware Refinement
#'
#' Re-marks membership flags against expression-filtered target columns,
#' appends \code{Mean_Expression} and \code{Passes_Expression_Filter},
#' and retains only rows where at least one membership flag is still
#' \code{TRUE}. Evidence labels such as \code{local_promoter_overlap}
#' and \code{linear_fallback} are preserved.
#'
#' @param target_gene_links Data frame from
#'   \code{\link{.build_target_gene_links}}.
#' @param bed_info Target annotation data frame after expression filtering.
#' @param vals Named numeric vector of per-gene mean expression.
#' @param threshold Numeric. Expression threshold.
#' @return A data frame with columns: \code{input_id}, \code{loop_ID},
#'   \code{anchor_id}, \code{gene}, \code{Mean_Expression},
#'   \code{Passes_Expression_Filter}, \code{gene_role}, \code{source},
#'   \code{evidence}, \code{anchor_role},
#'   \code{anchor_type_before_chromatin},
#'   \code{anchor_type_after_chromatin},
#'   \code{chromatin_confidence},
#'   \code{chromatin_evidence},
#'   \code{chromatin_target_action},
#'   \code{used_as_fallback},
#'   \code{in_regulated_promoter}, \code{in_assigned_target},
#'   \code{in_all_loop_connected}, \code{in_regulated_promoter_filled},
#'   \code{in_assigned_target_filled}.
#' @keywords internal
#' @noRd
.filter_refined_target_gene_links <- function(target_gene_links, bed_info, vals, threshold) {
    refined_cols <- c(
        "input_id", "loop_ID", "anchor_id", "gene", "Mean_Expression",
        "Passes_Expression_Filter", "gene_role", "source", "evidence",
        "anchor_role",
        # Preserve chromatin provenance when chromatin refinement ran
        # before expression refinement (both pipeline orders supported).
        "anchor_type_before_chromatin",
        "anchor_type_after_chromatin",
        "chromatin_confidence",
        "chromatin_evidence",
        "chromatin_target_action",
        "used_as_fallback",
        "in_regulated_promoter",
        "in_assigned_target",
        "in_all_loop_connected",
        "in_regulated_promoter_filled",
        "in_assigned_target_filled"
    )
    empty_refined <- function() {
        out <- .empty_target_gene_links()
        out$Mean_Expression <- numeric()
        out$Passes_Expression_Filter <- logical()
        out[, base::intersect(refined_cols, colnames(out)), drop = FALSE]
    }

    link_df <- .mark_target_gene_link_membership(target_gene_links, bed_info)
    if (nrow(link_df) == 0) {
        return(empty_refined())
    }

    expr <- unname(vals[link_df$gene])
    link_df$Mean_Expression <- as.numeric(expr)
    link_df$Passes_Expression_Filter <- !is.na(link_df$Mean_Expression) &
        link_df$Mean_Expression >= threshold
    link_df$retained_after_refinement <- link_df$in_regulated_promoter |
        link_df$in_assigned_target |
        link_df$in_all_loop_connected |
        link_df$in_regulated_promoter_filled |
        link_df$in_assigned_target_filled

    link_df %>%
        dplyr::filter(retained_after_refinement) %>%
        dplyr::select(dplyr::any_of(refined_cols)) %>%
        dplyr::distinct()
}

#' Internal: Build Karyotype Gene GRanges
#'
#' Shared helper that loads the full gene catalog, maps gene IDs to SYMBOLs,
#' and subsets to a given gene list. Used by karyotype heatmap sections in
#' both \code{\link{build_annotation_plots}} and related refinement plots.
#'
#' @param gene_symbols Character vector of gene symbols to retain.
#' @param txdb_obj A \code{TxDb} object for gene coordinates.
#' @param org_db_pkg Character. OrgDb package name for symbol mapping.
#' @param context Character. Diagnostic label for mapping messages.
#' @return A \code{GRanges} object subset to matching genes.
#' @keywords internal
#' @noRd
.build_karyo_gene_gr <- function(gene_symbols, txdb_obj, org_db_pkg, context,
                                 annotated_genes_gr = NULL) {
    if (!is.null(annotated_genes_gr)) {
        return(annotated_genes_gr[
            S4Vectors::mcols(annotated_genes_gr)$SYMBOL %in% gene_symbols
        ])
    }
    all_genes_gr <- .with_known_upstream_noise_suppressed(
        GenomicFeatures::genes(txdb_obj)
    )
    map <- .map_txdb_gene_ids(
        gene_ids = .extract_txdb_gene_ids(all_genes_gr),
        org_db = org_db_pkg,
        columns = "SYMBOL",
        context = context,
        warn = FALSE
    )
    S4Vectors::mcols(all_genes_gr)$SYMBOL <- map$SYMBOL[match(
        .extract_txdb_gene_ids(all_genes_gr), map$gene_id
    )]
    all_genes_gr[S4Vectors::mcols(all_genes_gr)$SYMBOL %in% gene_symbols]
}

#' Internal: Build Loop Type Donut Plot
#'
#' @param plot_df Loop annotation data frame.
#' @param custom_colors Named color vector keyed by loop_type.
#' @param project_name Character. Project prefix for the plot title.
#' @return A \code{ggplot} object.
#' @keywords internal
#' @noRd
.build_donut_plot <- function(plot_df, custom_colors, project_name) {
    donut_data <- plot_df %>%
        dplyr::group_by(loop_type) %>%
        dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
        dplyr::mutate(
            prop = n / sum(n),
            legend_label = paste0(loop_type, " (n=", n, ", ", round(prop * 100, 1), "%)")
        ) %>%
        dplyr::arrange(dplyr::desc(n))

    donut_data$loop_type <- factor(donut_data$loop_type, levels = donut_data$loop_type)

    ggplot2::ggplot(donut_data, ggplot2::aes(x = 2, y = n, fill = loop_type)) +
        ggplot2::geom_bar(stat = "identity", color = "white") +
        ggplot2::coord_polar(theta = "y") +
        ggplot2::xlim(0.5, 2.9) +
        ggplot2::geom_text(ggplot2::aes(x = 2.65, label = loop_type),
            position = ggplot2::position_stack(vjust = 0.5), size = 2.5
        ) +
        ggplot2::scale_fill_manual(
            values = rev(custom_colors),
            labels = setNames(donut_data$legend_label, donut_data$loop_type)
        ) +
        ggplot2::theme_void() +
        ggplot2::labs(title = paste0(project_name, ": Loop Type Distribution")) +
        ggplot2::theme(
            legend.position = "right",
            plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
        )
}

#' Internal: Build Karyotype Loop-Genes Heatmap
#'
#' @param plot_df Loop annotation data frame with \code{All_Anchor_Genes} column.
#' @param txdb_obj A \code{TxDb} object.
#' @param org_db_pkg Character. OrgDb package name.
#' @param species Character. Genome assembly.
#' @param karyo_bin_size Integer. Bin size for karyotype heatmaps.
#' @param red_palette Character vector. Red color palette.
#' @return A karyotype grob, or \code{NULL}.
#' @keywords internal
#' @noRd
.build_karyo_loop_genes_plot <- function(
  plot_df, txdb_obj, org_db_pkg, species, karyo_bin_size, red_palette,
  annotated_genes_gr = NULL
) {
    if (!"All_Anchor_Genes" %in% colnames(plot_df)) {
        return(NULL)
    }
    genes_loop <- clean_gene_names(plot_df$All_Anchor_Genes, ";")
    if (length(genes_loop) == 0) {
        return(NULL)
    }
    target_genes_gr <- .build_karyo_gene_gr(
        genes_loop, txdb_obj, org_db_pkg,
        "build_annotation_plots loop-gene karyotype",
        annotated_genes_gr = annotated_genes_gr
    )
    draw_karyo_heatmap_internal(
        target_genes_gr, "Loop Genes Distribution", karyo_bin_size, 0.99,
        txdb_obj, species, "Genes",
        custom_colors = red_palette
    )
}

#' Internal: Build Karyotype Anchor-Load Heatmap
#'
#' @param plot_df Loop annotation data frame with coordinate columns.
#' @param txdb_obj A \code{TxDb} object.
#' @param species Character. Genome assembly.
#' @param karyo_bin_size Integer. Bin size for karyotype heatmaps.
#' @param blue_palette Character vector. Blue color palette.
#' @return A karyotype grob, or \code{NULL}.
#' @keywords internal
#' @noRd
.build_karyo_anchors_plot <- function(
  plot_df, txdb_obj, species, karyo_bin_size, blue_palette
) {
    all_anchors <- dplyr::bind_rows(
        plot_df %>% dplyr::select(chr = chr1, start = start1, end = end1),
        plot_df %>% dplyr::select(chr = chr2, start = start2, end = end2)
    ) %>% dplyr::distinct()
    if (nrow(all_anchors) == 0) {
        return(NULL)
    }
    draw_karyo_heatmap_internal(
        .with_known_upstream_noise_suppressed(
            GenomicRanges::makeGRangesFromDataFrame(all_anchors)
        ),
        "Loop Anchor Load", karyo_bin_size, 0.99, txdb_obj, species, "Anchors",
        custom_colors = blue_palette
    )
}

#' Internal: Build Simplified Flower Plot for Annotation
#'
#' @param plot_df Loop annotation data frame.
#' @param project_name Character. Project prefix for the plot title.
#' @param custom_colors Named color vector keyed by loop_type.
#' @return A \code{ggplot} object, or \code{NULL} if fewer than 2 gene sets.
#' @keywords internal
#' @noRd
.build_flower_plot <- function(plot_df, project_name, custom_colors) {
    temp_df_flower <- plot_df %>%
        dplyr::filter(!is.na(All_Anchor_Genes) & All_Anchor_Genes != "") %>%
        tidyr::separate_rows(All_Anchor_Genes, sep = ";") %>%
        dplyr::mutate(All_Anchor_Genes = trimws(All_Anchor_Genes)) %>%
        dplyr::filter(All_Anchor_Genes != "")
    gene_sets <- split(temp_df_flower$All_Anchor_Genes, temp_df_flower$loop_type)
    gene_sets <- lapply(gene_sets, unique)
    if (length(gene_sets) <= 1) {
        return(NULL)
    }
    draw_flower_simplified(gene_sets, project_name, custom_colors)
}

#' Internal: Build Target-Connected Rose Plot
#'
#' @param target_connected_loops Data frame of target-connected loops.
#' @param custom_colors Named color vector keyed by loop_type.
#' @param project_name Character. Project prefix for the plot title.
#' @return A \code{ggplot} object, or \code{NULL}.
#' @keywords internal
#' @noRd
.build_target_rose_plot <- function(
  target_connected_loops, custom_colors, project_name
) {
    if (is.null(target_connected_loops) || nrow(target_connected_loops) == 0) {
        return(NULL)
    }
    rose_data <- target_connected_loops %>%
        dplyr::group_by(loop_type) %>%
        dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
        dplyr::mutate(
            prop = n / sum(n),
            legend_label = paste0(loop_type, " (n=", n, ", ", round(prop * 100, 1), "%)")
        ) %>%
        dplyr::arrange(dplyr::desc(n))

    rose_data$loop_type <- factor(rose_data$loop_type, levels = rose_data$loop_type)

    ggplot2::ggplot(rose_data, ggplot2::aes(x = loop_type, y = n, fill = loop_type)) +
        ggplot2::geom_bar(stat = "identity", width = 1, color = "white") +
        ggplot2::coord_polar(theta = "x") +
        ggplot2::scale_fill_manual(
            values = custom_colors,
            labels = setNames(rose_data$legend_label, rose_data$loop_type)
        ) +
        ggplot2::theme_void() +
        ggplot2::labs(title = paste0(project_name, ": Target-Linked Loop Types")) +
        ggplot2::theme(
            legend.position = "right",
            plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
        )
}

#' Internal: Build Target Genomic Distribution Pie Charts
#'
#' @param bed_info Target annotation data frame.
#' @param color_palette Character. RColorBrewer palette name.
#' @param project_name Character. Project prefix for plot titles.
#' @return A named list of \code{ggplot} objects (may be empty).
#' @keywords internal
#' @noRd
.build_target_genomic_pies <- function(bed_info, color_palette, project_name) {
    plot_list <- list()
    if (is.null(bed_info)) {
        return(plot_list)
    }
    plot_list$Target_Genomic_Distribution <- draw_pie_with_outside_labels(
        bed_info, "annotation",
        paste0(project_name, ": All Targets Genomic Distribution"), color_palette
    )
    if ("Linked_Loop_IDs" %in% colnames(bed_info)) {
        linked_rows <- bed_info[!is.na(bed_info$Linked_Loop_IDs) &
            bed_info$Linked_Loop_IDs != "", ]
        if (nrow(linked_rows) > 0) {
            plot_list$Target_Loop_Genomic_Distribution <- draw_pie_with_outside_labels(
                linked_rows, "annotation",
                paste0(project_name, ": Loop-Connected Targets Distribution"),
                color_palette
            )
        }
    }
    plot_list
}

#' Internal: Build Karyotype Target-Genes Heatmap
#'
#' @param bed_info Target annotation data frame (or NULL).
#' @param txdb_obj A \code{TxDb} object.
#' @param org_db_pkg Character. OrgDb package name.
#' @param species Character. Genome assembly.
#' @param karyo_bin_size Integer. Bin size for karyotype heatmaps.
#' @param purple_palette Character vector. Purple color palette.
#' @return A karyotype grob, or \code{NULL}.
#' @keywords internal
#' @noRd
.build_karyo_target_genes_plot <- function(
  bed_info, txdb_obj, org_db_pkg, species,
  karyo_bin_size, purple_palette, annotated_genes_gr = NULL
) {
    if (is.null(bed_info) ||
        !"Assigned_Target_Genes_Filled" %in% colnames(bed_info)) {
        return(NULL)
    }
    genes_target <- clean_gene_names(
        bed_info$Assigned_Target_Genes_Filled, ";"
    )
    if (length(genes_target) == 0) {
        return(NULL)
    }
    target_genes_gr <- .build_karyo_gene_gr(
        genes_target, txdb_obj, org_db_pkg,
        "build_annotation_plots target-gene karyotype",
        annotated_genes_gr = annotated_genes_gr
    )
    draw_karyo_heatmap_internal(
        target_genes_gr, "Target Genes (Assigned+Local)", karyo_bin_size,
        0.99, txdb_obj, species, "Genes",
        custom_colors = purple_palette
    )
}

#' Internal: Build Annotation Visualization Suite
#'
#' Generates the complete diagnostic plot collection (donut, rose, karyotype
#' heatmaps, flower plot, pie charts) for the annotation results.
#'
#' @param plot_df Loop annotation data frame with loop_type and All_Anchor_Genes.
#' @param bed_info Target annotation data frame (optional).
#' @param cluster_info Anchor annotation data frame.
#' @param target_connected_loops Data frame of target-connected loops (optional).
#' @param txdb_obj TxDb object for gene coordinate lookup.
#' @param org_db_pkg Character. Organism database package name.
#' @param species Character. Genome assembly.
#' @param project_name Character. Project prefix for plot titles.
#' @param color_palette Character. RColorBrewer palette name.
#' @param karyo_bin_size Integer. Bin size for karyotype heatmaps.
#' @return A named list of ggplot / grob objects.
#' @keywords internal
#' @noRd
build_annotation_plots <- function(
  plot_df, bed_info, cluster_info,
  target_connected_loops, txdb_obj, org_db_pkg, species, project_name,
  color_palette, karyo_bin_size
) {
    loop_types_sorted <- sort(unique(plot_df$loop_type))
    custom_colors <- get_colors(length(loop_types_sorted), color_palette)
    names(custom_colors) <- loop_types_sorted

    red_palette <- c("#FFFFFF", "#FFFFCC", "#FFEDA0", "#FED976", "#FEB24C", "#FD8D3C", "#FC4E2A", "#E31A1C", "#BD0026", "#800026", "#000000")
    blue_palette <- c("#FFFFFF", "#E1F5FE", "#B3E5FC", "#4FC3F7", "#039BE5", "#0277BD", "#01579B", "#000000")
    purple_palette <- c("#FFFFFF", "#F3E5F5", "#E1BEE7", "#BA68C8", "#9C27B0", "#7B1FA2", "#4A148C", "#000000")

    # Pre-load and annotate gene catalogue once for all karyotype plots
    need_loop_genes <- "All_Anchor_Genes" %in% colnames(plot_df) &&
        length(clean_gene_names(plot_df$All_Anchor_Genes, ";")) > 0
    need_target_genes <- !is.null(bed_info) &&
        "Assigned_Target_Genes_Filled" %in% colnames(bed_info) &&
        length(clean_gene_names(bed_info$Assigned_Target_Genes_Filled, ";")) > 0
    annotated_genes_gr <- NULL
    if (need_loop_genes || need_target_genes) {
        all_genes_gr <- .with_known_upstream_noise_suppressed(
            GenomicFeatures::genes(txdb_obj)
        )
        map <- .map_txdb_gene_ids(
            gene_ids = .extract_txdb_gene_ids(all_genes_gr),
            org_db = org_db_pkg,
            columns = "SYMBOL",
            context = "build_annotation_plots karyotype",
            warn = FALSE
        )
        S4Vectors::mcols(all_genes_gr)$SYMBOL <- map$SYMBOL[match(
            .extract_txdb_gene_ids(all_genes_gr), map$gene_id
        )]
        annotated_genes_gr <- all_genes_gr
    }

    plot_list <- list()

    plot_list$Basic_Donut <- .build_donut_plot(plot_df, custom_colors, project_name)

    plot_list$Basic_Circular <- draw_circular_bar_plot(plot_df, project_name,
        color_vec = custom_colors
    )

    plot_list$Karyo_LoopGenes <- .build_karyo_loop_genes_plot(
        plot_df, txdb_obj, org_db_pkg, species, karyo_bin_size, red_palette,
        annotated_genes_gr = annotated_genes_gr
    )

    plot_list$Karyo_Anchors <- .build_karyo_anchors_plot(
        plot_df, txdb_obj, species, karyo_bin_size, blue_palette
    )

    plot_list$Basic_Flower <- .build_flower_plot(
        plot_df, project_name, custom_colors
    )

    if (nrow(cluster_info) > 0) {
        plot_list$Anchor_Genomic_Distribution <- draw_pie_with_outside_labels(
            cluster_info, "annotation",
            paste0(project_name, ": Anchor Loci Genomic Distribution"), color_palette
        )
    }

    if (!is.null(bed_info)) {
        plot_list$Karyo_TargetGenes <- .build_karyo_target_genes_plot(
            bed_info, txdb_obj, org_db_pkg, species, karyo_bin_size, purple_palette,
            annotated_genes_gr = annotated_genes_gr
        )

        plot_list$Target_Rose <- .build_target_rose_plot(
            target_connected_loops, custom_colors, project_name
        )

        target_pies <- .build_target_genomic_pies(bed_info, color_palette, project_name)
        for (nm in names(target_pies)) plot_list[[nm]] <- target_pies[[nm]]
    }

    return(plot_list)
}

# --- Internal helpers extracted from refine_loop_anchors_by_expression ---

#' Internal: Load and validate data for expression-aware refinement
#' @return A list with loop_df, original_loop_df, clust_info, bed_info,
#'   target_gene_links, whitelist, vals, has_cluster_id, upstream_promoter_stats.
#' @keywords internal
#' @noRd
.refine_load_validate_data <- function(
    annotation_res, expr_matrix_file, sample_columns,
    threshold, unit_type, log_message
) {
    log_message(">>> [Step 1] Loading Data & Expression Matrix...")
    if (is.null(annotation_res$loop_annotation)) stop("'loop_annotation' missing.")

    original_loop_df <- annotation_res$loop_annotation
    loop_df <- annotation_res$loop_annotation
    clust_info <- annotation_res$anchor_loci_annotation
    if (is.null(clust_info)) {
        clust_info <- annotation_res$anchor_annotation
    }
    bed_info <- annotation_res$target_annotation
    target_gene_links <- annotation_res$target_gene_links

    required_cols <- c(
        "chr1", "start1", "end1", "chr2", "start2", "end2",
        "anchor1_gene", "anchor1_type", "anchor2_gene", "anchor2_type",
        "Putative_Target_Genes"
    )
    missing_cols <- setdiff(required_cols, colnames(loop_df))
    if (length(missing_cols) > 0) {
        stop(
            "annotation_res$loop_annotation is missing required columns: ",
            paste(missing_cols, collapse = ", "), ". ",
            "Ensure the input was generated by annotate_peaks_and_loops()."
        )
    }
    if (!"cluster_id" %in% colnames(loop_df)) {
        warning("'cluster_id' not found in loop_annotation; cluster-level stats and target-donut will be skipped.", call. = FALSE)
        loop_df$cluster_id <- NA_character_
        has_cluster_id <- FALSE
    } else {
        has_cluster_id <- TRUE
    }

    if (!"a1_id" %in% colnames(loop_df) || !"a2_id" %in% colnames(loop_df)) {
        log_message("    Reconstructing anchor IDs from coordinates (omitted from upstream loop_annotation output).")
        loop_df <- loop_df %>%
            dplyr::mutate(
                a1_id = paste(chr1, start1, end1, sep = "_"),
                a2_id = paste(chr2, start2, end2, sep = "_")
            )
    }

    upstream_promoter_stats <- annotation_res$promoter_centric_stats

    if (is.null(expr_matrix_file) || is.null(sample_columns)) {
        stop("Expression matrix file and sample columns are required for refinement.")
    }
    vals <- load_expression_matrix(expr_matrix_file, sample_columns)
    whitelist <- names(vals)[vals >= threshold & !is.na(vals) & names(vals) != ""]
    log_message(sprintf("    >>> Active Genes (>= %s %s): %d", threshold, unit_type, length(whitelist)))

    anno_genes <- unique(c(
        trimws(unlist(strsplit(na.omit(loop_df$anchor1_gene), ";"))),
        trimws(unlist(strsplit(na.omit(loop_df$anchor2_gene), ";")))
    ))
    anno_genes <- anno_genes[nzchar(anno_genes)]
    if (length(anno_genes) > 0) {
        overlap_rate <- length(intersect(whitelist, anno_genes)) / length(anno_genes)
        if (overlap_rate < 0.1) {
            warning(
                sprintf(
                    "Only %.1f%% of annotation gene symbols match the expression matrix row names. ",
                    overlap_rate * 100
                ),
                "Check that expression matrix row names use the same gene identifier convention (e.g., SYMBOL)."
            )
        }
    }

    list(
        loop_df = loop_df,
        original_loop_df = original_loop_df,
        clust_info = clust_info,
        bed_info = bed_info,
        target_gene_links = target_gene_links,
        whitelist = whitelist,
        vals = vals,
        has_cluster_id = has_cluster_id,
        upstream_promoter_stats = upstream_promoter_stats
    )
}

#' Internal: Reclassify anchors and compute expression-filtered targets
#' @return The updated loop_df.
#' @keywords internal
#' @noRd
.refine_apply_anchor_updates <- function(
    loop_df, whitelist, reclassify_by_expression, log_message
) {
    log_message(">>> [Step 2] Updating Anchors & Loops...")
    orig_anchor1_type <- loop_df$anchor1_type
    orig_anchor2_type <- loop_df$anchor2_type
    loop_df <- .refine_reclassify_anchors(loop_df, whitelist, reclassify_by_expression)
    original_ptg <- loop_df$Putative_Target_Genes

    loop_df <- .refine_compute_targets(
        loop_df, original_ptg, whitelist,
        orig_anchor1_type, orig_anchor2_type
    )

    log_message(sprintf(
        "    Retained: %d / %d loops",
        sum(loop_df$Has_Active_Target), nrow(loop_df)
    ))

    if (isTRUE(reclassify_by_expression)) {
        n_reclassified <- sum(
            orig_anchor1_type == "P" & loop_df$anchor1_type == "eP" |
                orig_anchor2_type == "P" & loop_df$anchor2_type == "eP" |
                orig_anchor1_type == "G" & loop_df$anchor1_type == "eG" |
                orig_anchor2_type == "G" & loop_df$anchor2_type == "eG",
            na.rm = TRUE
        )
        n_pg <- sum(
            orig_anchor1_type %in% c("P", "G") |
                orig_anchor2_type %in% c("P", "G"),
            na.rm = TRUE
        )
        if (n_pg > 0 && n_reclassified / n_pg > 0.5) {
            warning(
                round(n_reclassified / n_pg * 100),
                "% of P/G anchors were reclassified to eP/eG. ",
                "eP/eG labels indicate transcriptionally inactive promoter/gene-body states; enhancer activity requires orthogonal chromatin evidence. ",
                "Validate with orthogonal chromatin data (ATAC-seq, H3K27ac) ",
                "before interpreting eP/eG anchors as functional enhancers.",
                call. = FALSE
            )
        }
    }
    loop_df
}

# --- Internal helpers for refine_loop_anchors_by_expression ---

#' Reclassify anchor types based on expression whitelist.
#' @keywords internal
#' @noRd
.refine_reclassify_anchors <- function(loop_df, whitelist, reclassify_by_expression) {
    a1_res <- Map(
        function(g, t) clean_anchor(g, t, allow = whitelist, down = reclassify_by_expression),
        loop_df$anchor1_gene, loop_df$anchor1_type
    )
    a2_res <- Map(
        function(g, t) clean_anchor(g, t, allow = whitelist, down = reclassify_by_expression),
        loop_df$anchor2_gene, loop_df$anchor2_type
    )
    loop_df$anchor1_type <- vapply(a1_res, function(x) x$type, character(1))
    loop_df$anchor1_gene <- vapply(a1_res, function(x) x$gene, character(1))
    loop_df$anchor2_type <- vapply(a2_res, function(x) x$type, character(1))
    loop_df$anchor2_gene <- vapply(a2_res, function(x) x$gene, character(1))
    loop_df$loop_type <- unlist(Map(
        function(t1, t2) paste(sort(c(t1, t2)), collapse = "-"),
        loop_df$anchor1_type, loop_df$anchor2_type
    ), use.names = FALSE)
    loop_df
}

#' Compute active/fallback target genes and refinement status flags.
#' @keywords internal
#' @noRd
.refine_compute_targets <- function(loop_df, original_ptg, whitelist, orig_anchor1_type = NULL, orig_anchor2_type = NULL) {
    filter_genes_wl <- function(x) {
        if (is.na(x) || x == "") {
            return(NA_character_)
        }
        gs <- trimws(unlist(strsplit(as.character(x), ";")))
        gs <- unique(gs[gs %in% whitelist])
        if (length(gs) == 0) {
            return(NA_character_)
        }
        paste(sort(gs), collapse = ";")
    }
    filtered_ptg <- vapply(original_ptg, filter_genes_wl, character(1))

    is_enh_like <- .is_distal_like
    is_promoter <- .is_promoter_like
    is_gene_body <- .is_gene_body_like

    loop_df <- loop_df %>%
        dplyr::rowwise() %>%
        dplyr::mutate(
            .fallback_ptg = dplyr::case_when(
                (is_promoter(anchor1_type) & is_enh_like(anchor2_type)) ~ extract_genes(anchor1_gene),
                (is_enh_like(anchor1_type) & is_promoter(anchor2_type)) ~ extract_genes(anchor2_gene),
                (is_promoter(anchor1_type) & is_promoter(anchor2_type)) ~ extract_genes(c(anchor1_gene, anchor2_gene)),
                (is_promoter(anchor1_type) & is_gene_body(anchor2_type)) ~ extract_genes(anchor1_gene),
                (is_gene_body(anchor1_type) & is_promoter(anchor2_type)) ~ extract_genes(anchor2_gene),
                (is_gene_body(anchor1_type) & is_enh_like(anchor2_type)) ~ extract_genes(anchor1_gene),
                (is_enh_like(anchor1_type) & is_gene_body(anchor2_type)) ~ extract_genes(anchor2_gene),
                TRUE ~ extract_genes(c(anchor1_gene, anchor2_gene))
            )
        ) %>%
        dplyr::ungroup()

    loop_df$Active_Target_Genes <- filtered_ptg
    loop_df$Putative_Target_Genes <- filtered_ptg
    empty_idx <- is.na(loop_df$Putative_Target_Genes) | loop_df$Putative_Target_Genes == ""
    loop_df$Putative_Target_Genes[empty_idx] <- loop_df$.fallback_ptg[empty_idx]
    loop_df <- loop_df %>% dplyr::select(-.fallback_ptg)

    has_active_target <- !is.na(loop_df$Active_Target_Genes) & loop_df$Active_Target_Genes != ""
    if (!is.null(orig_anchor1_type) && !is.null(orig_anchor2_type)) {
        reclassified_a1 <- orig_anchor1_type != loop_df$anchor1_type & loop_df$anchor1_type %in% c("eP", "eG")
        reclassified_a2 <- orig_anchor2_type != loop_df$anchor2_type & loop_df$anchor2_type %in% c("eP", "eG")
        has_reclassified <- reclassified_a1 | reclassified_a2
    } else {
        has_reclassified <- loop_df$anchor1_type %in% c("eP", "eG") | loop_df$anchor2_type %in% c("eP", "eG")
    }
    had_original <- !is.na(original_ptg) & original_ptg != ""

    loop_df$Has_Active_Target <- has_active_target
    loop_df$Refinement_Action <- dplyr::case_when(
        has_active_target ~ "retained_active_target",
        !has_active_target & has_reclassified & had_original ~ "reclassified_silent_anchor",
        !has_active_target & had_original ~ "expression_filtered_no_active_target",
        TRUE ~ "structural_only_no_active_target"
    )
    loop_df$Retained_In_Functional_Network <- has_active_target
    loop_df
}

#' Refine target BED annotations by expression filtering.
#' @keywords internal
#' @noRd
.refine_target_annotations <- function(bed_info, loop_df, whitelist, target_gene_links, vals, threshold) {
    cols_to_clean <- grep("Strict|Physical|Loop_Genes|promoter|Filled|Target_Genes|Assigned", colnames(bed_info), value = TRUE)
    cols_to_clean <- cols_to_clean[!grepl("Evidence|Linked_Loop_IDs|SANKEY_RAW_GENES", cols_to_clean)]
    raw_tgt_col <- "Assigned_Target_Genes_Filled"
    if (!raw_tgt_col %in% colnames(bed_info)) {
        raw_tgt_col <- grep("Filled", cols_to_clean, value = TRUE)[1]
    }
    if (!is.na(raw_tgt_col) && raw_tgt_col %in% colnames(bed_info)) bed_info$SANKEY_RAW_GENES <- bed_info[[raw_tgt_col]]
    for (col in cols_to_clean) {
        if (col %in% colnames(bed_info)) {
            bed_info[[col]] <- vapply(as.character(bed_info[[col]]), function(x) {
                if (is.na(x) || x == "") {
                    return(NA_character_)
                }
                gs <- unlist(strsplit(x, ";"))
                gs_active <- gs[trimws(gs) %in% whitelist]
                if (length(gs_active) == 0) {
                    return(NA_character_)
                }
                paste(unique(sort(trimws(gs_active))), collapse = ";")
            }, FUN.VALUE = character(1))
        }
    }

    any_sym_in_wl <- function(s) {
        if (is.na(s) || s == "") {
            return(FALSE)
        }
        any(trimws(unlist(strsplit(as.character(s), ";"))) %in% whitelist)
    }
    filter_sym_expressed <- function(s) {
        if (is.na(s) || s == "") {
            return(NA_character_)
        }
        gs <- trimws(unlist(strsplit(as.character(s), ";")))
        gs <- gs[gs %in% whitelist]
        if (length(gs) == 0) {
            return(NA_character_)
        }
        paste(sort(gs), collapse = ";")
    }
    linear_col <- .target_linear_gene_column(bed_info)
    if (!is.null(linear_col)) {
        has_sym <- vapply(bed_info[[linear_col]], any_sym_in_wl, logical(1))
        sym_fill <- vapply(bed_info[[linear_col]], filter_sym_expressed, character(1))
    } else {
        has_sym <- rep(FALSE, nrow(bed_info))
        sym_fill <- rep(NA_character_, nrow(bed_info))
    }

    if ("Regulated_promoter_genes" %in% colnames(bed_info) &&
        "Regulated_promoter_genes_Filled" %in% colnames(bed_info)) {
        has_reg <- !is.na(bed_info$Regulated_promoter_genes) & bed_info$Regulated_promoter_genes != ""
        bed_info$Regulated_promoter_genes_Filled <- dplyr::case_when(
            has_reg ~ bed_info$Regulated_promoter_genes,
            has_sym ~ sym_fill,
            TRUE ~ NA_character_
        )
        if ("Regulated_promoter_Fallback_Evidence" %in% colnames(bed_info)) {
            bed_info$Regulated_promoter_Fallback_Evidence <- dplyr::case_when(
                has_reg ~ "none",
                has_sym ~ bed_info$Regulated_promoter_Fallback_Evidence,
                TRUE ~ "none"
            )
        }
    }

    if ("Linked_Loop_IDs" %in% colnames(bed_info) &&
        "loop_ID" %in% colnames(loop_df) &&
        "Active_Target_Genes" %in% colnames(loop_df)) {
        loop_tgt <- loop_df %>%
            dplyr::filter(!is.na(Active_Target_Genes) & Active_Target_Genes != "") %>%
            dplyr::select(loop_ID, Active_Target_Genes) %>%
            dplyr::distinct()
        get_loop_tgt <- function(linked) {
            if (is.na(linked) || linked == "") {
                return(NA_character_)
            }
            ids <- trimws(unlist(strsplit(as.character(linked), ";")))
            tgt <- loop_tgt$Active_Target_Genes[match(ids, loop_tgt$loop_ID)]
            genes <- unique(trimws(unlist(strsplit(na.omit(tgt), ";"))))
            genes <- genes[genes != ""]
            if (length(genes) == 0) {
                return(NA_character_)
            }
            paste(sort(genes), collapse = ";")
        }
        loop_tgt_vec <- vapply(bed_info$Linked_Loop_IDs, get_loop_tgt, FUN.VALUE = character(1))
        if ("Assigned_Target_Genes" %in% colnames(bed_info)) {
            empty <- is.na(bed_info$Assigned_Target_Genes) | bed_info$Assigned_Target_Genes == ""
            fill_ok <- !is.na(loop_tgt_vec) & loop_tgt_vec != ""
            bed_info$Assigned_Target_Genes[empty & fill_ok] <- loop_tgt_vec[empty & fill_ok]
            if ("Assigned_Target_Genes_Filled" %in% colnames(bed_info)) {
                has_tgt <- !is.na(bed_info$Assigned_Target_Genes) & bed_info$Assigned_Target_Genes != ""
                bed_info$Assigned_Target_Genes_Filled <- dplyr::case_when(
                    has_tgt ~ bed_info$Assigned_Target_Genes,
                    has_sym ~ sym_fill,
                    TRUE ~ NA_character_
                )
            }
        }
    }

    if (!is.null(target_gene_links)) {
        target_gene_links <- .filter_refined_target_gene_links(
            target_gene_links, bed_info, vals, threshold
        )
    }

    list(bed_info = bed_info, target_gene_links = target_gene_links)
}

#' Export refined results to Excel workbook.
#' @keywords internal
#' @noRd
.refine_export_workbook <- function(loop_df, clust_info, promoter_centric_df,
    distal_element_df, bed_info, target_gene_links, out_dir, project_name,
    chromatin_validation = NULL) {
    wb <- openxlsx::createWorkbook()
    loop_export <- loop_df %>%
        dplyr::select(-any_of(c("a1_id", "a2_id", "loop_genes", "single_loop_genes", "proximate_loop_gene")))

    openxlsx::addWorksheet(wb, "Filtered Loop Annotation")
    openxlsx::writeData(wb, "Filtered Loop Annotation", loop_export)

    functional_loops <- loop_export %>% dplyr::filter(Retained_In_Functional_Network == TRUE)
    openxlsx::addWorksheet(wb, "Functional Loop Annotation")
    openxlsx::writeData(wb, "Functional Loop Annotation", functional_loops)

    openxlsx::addWorksheet(wb, "Filtered Anchor Loci Annotation")
    openxlsx::writeData(wb, "Filtered Anchor Loci Annotation", clust_info)

    if (!is.null(promoter_centric_df)) {
        openxlsx::addWorksheet(wb, "Filtered Promoter Stats")
        openxlsx::writeData(wb, "Filtered Promoter Stats", promoter_centric_df)
    }
    if (!is.null(distal_element_df)) {
        openxlsx::addWorksheet(wb, "Filtered Distal Element Stats")
        openxlsx::writeData(wb, "Filtered Distal Element Stats", distal_element_df)
    }
    if (!is.null(bed_info)) {
        bed_export <- bed_info %>% dplyr::select(-any_of("SANKEY_RAW_GENES"))
        openxlsx::addWorksheet(wb, "Filtered Target Annotation")
        openxlsx::writeData(wb, "Filtered Target Annotation", bed_export)
    }
    if (!is.null(target_gene_links) && nrow(target_gene_links) > 0) {
        openxlsx::addWorksheet(wb, "Filtered Target Gene Links")
        openxlsx::writeData(wb, "Filtered Target Gene Links", target_gene_links)
    }
    if (!is.null(chromatin_validation) && nrow(chromatin_validation) > 0) {
        openxlsx::addWorksheet(wb, "Chromatin Validation")
        openxlsx::writeData(wb, "Chromatin Validation", chromatin_validation)
    }

    tryCatch(
        openxlsx::saveWorkbook(wb, file.path(out_dir, paste0(project_name, "_Refined_Results.xlsx")), overwrite = TRUE),
        error = function(e) warning("Failed to save refined Excel workbook: ", conditionMessage(e), call. = FALSE)
    )
}

#' @title Expression-Aware refinement of loop anchors and target linkages
#'
#' @description
#' Integrates quantitative RNA-seq data (e.g., TPM/FPKM) with 3D structural data
#' to reclassify regulatory elements and annotate functional status, deriving a
#' functionally interpretable regulatory network from physical chromatin contacts.
#' All structural loops are preserved; refinement status columns indicate which
#' loops belong to the high-confidence active subset.
#'
#' @details
#' \strong{Algorithmic Framework:}
#' \itemize{
#'   \item \strong{Target Filtration:} Parses merged gene assignments (e.g., \code{"GeneA;GeneB"}), evaluates individual genes against a defined expression threshold, and retains only transcriptionally active targets.
#'   \item \strong{Biological Reclassification:} Reclassifies physically annotated promoters (\code{P}) and gene bodies (\code{G}) lacking active transcription as \emph{expression-filtered silent} regulatory elements (\code{eP}, \code{eG}). \strong{Important:} \code{eP}/\code{eG} labels indicate transcriptional silence at the reference gene -- they do \strong{not} constitute evidence of enhancer activity. Orthogonal chromatin data (ATAC-seq, H3K27ac, H3K4me1) are required for functional enhancer interpretation. The labels are retained for backward compatibility; interpret them as "inactive-P" / "inactive-G" rather than "enhancer-P" / "enhancer-G".
#'   \item \strong{Expression-Aware Connectivity Statistics:} Recomputes promoter-centric and distal-element connectivity after expression-aware anchor refinement, while preserving all structural loops in the refined loop annotation. This separates the complete physical contact map from the high-confidence active subset.
#'   \item \strong{External Target Refinement:} Filters auxiliary target mapping columns (e.g., \code{Assigned_Target_Genes_Filled}) based on expression criteria, ensuring that mapped 1D genomic features are exclusively linked to active genes.
#'   \item \strong{Target Provenance Preservation:} Recomputes \code{*_Filled}
#'   membership flags in \code{target_gene_links} after expression filtering,
#'   retains only links still used by the refined target columns, and appends
#'   \code{Mean_Expression} plus \code{Passes_Expression_Filter}. Evidence labels
#'   such as \code{local_promoter_overlap}, \code{direct_opposite_promoter}, and
#'   \code{linear_fallback} are preserved.
#' }
#'
#' \strong{Design Philosophy:}
#' This function does not discard structural loops based on expression state.
#' Hi-C, HiChIP, and PLAC-seq capture 3D chromatin contacts; RNA-seq captures
#' current transcriptional state. A silent promoter may reflect cell-state,
#' time-point, or technical factors rather than absence of physical contact.
#' All structural loops are retained with refinement status columns, and a
#' high-confidence functional subset is provided via
#' \code{Retained_In_Functional_Network} and the \emph{Functional Loop Annotation}
#' Excel sheet.
#'
#' \strong{Interpretation of eP/eG labels:}
#' \code{eP} and \code{eG} are \strong{expression-filtered silent states}, not
#' functional enhancer classifications. Bulk RNA-seq silence can arise from
#' cell-type proportions, time-point effects, sequencing depth, or promoter
#' pausing -- none of which imply the locus has gained enhancer activity.
#' These labels should be read as "transcriptionally inactive P/G" and
#' treated as hypotheses requiring orthogonal validation (ATAC-seq, H3K27ac,
#' H3K4me1, or H3K27me3 depletion). Users with matched chromatin data should
#' overlay eP/eG loci against these tracks before interpreting them as
#' putative regulatory elements.
#'
#' @param annotation_res List. The raw foundational output object returned by \code{\link{annotate_peaks_and_loops}}.
#' @param expr_matrix_file Path to a normalised expression matrix (TPM/FPKM, genes x samples). Required for refinement. Default: \code{NULL}.
#' @param sample_columns Character vector or integer indices. Columns in \code{expr_matrix_file} to average. Default: \code{NULL}.
#' @param threshold Numeric. Minimum expression (e.g. TPM >= 1) for a gene to be considered active. Default: \code{1}.
#'   \strong{Gene name matching is case-sensitive.} Ensure gene identifiers in
#'   \code{expr_matrix_file} use the same case as the symbols returned by your
#'   OrgDb (typically all-uppercase for human and mouse, e.g. \code{"TP53"},
#'   not \code{"Tp53"}). Mismatched case will cause expressed genes to be
#'   misclassified as silent (eP/eG).
#' @param unit_type Character. Expression unit for plot labels (e.g., \code{"TPM"}). Default: \code{"TPM"}.
#' @param species Character. Genome assembly. One of \code{"hg38"}, \code{"hg19"}, \code{"mm10"}, \code{"mm9"}. Default: \code{"hg38"}. Extensible to other species via \code{species_txdb_pkg()} and related helpers.
#' @param out_dir Character. Output directory for the Excel results file. Default: \code{"./results/filtered"}.
#' @param project_name Character. Prefix for output files (automatically appends \code{"_Filtered"}). Default: \code{"HiChIP"}.
#' @param color_palette Character. RColorBrewer palette name for loop-type colour assignments. Default: \code{"Paired"}.
#' @param karyo_bin_size Integer. Bin width in bp for karyotype heatmaps. Default: \code{1e5}.
#' @param reclassify_by_expression Logical. If \code{TRUE} (default), silent promoters (P) and gene bodies (G) are reclassified as eP/eG.
#' @param hub_percentile Numeric (0-1). Node-degree quantile for hub detection. Default: \code{0.95}.
#' @param chromatin_beds Named list of BED file paths for orthogonal chromatin
#'   mark validation of eP/eG anchors. When non-empty, a \emph{Chromatin
#'   Validation} sheet is added to the Excel workbook (see
#'   \code{\link{validate_epeG_by_chromatin}}). Default: \code{list()} (skip).
#'   Note: if you plan to use \code{\link{refine_loop_anchors_by_chromatin}},
#'   the chromatin validation and reclassification are handled there; you can
#'   skip this parameter to avoid redundant BED file reads.
#' @param write_output Logical. If \code{TRUE} (default), write the refined Excel workbook to \code{out_dir}. If \code{FALSE}, return results without creating directories or files.
#' @param quiet Logical. If \code{TRUE}, suppress progress messages while preserving warnings. Default: \code{FALSE}.
#'
#' @return An invisible named list:
#' \itemize{
#'   \item \code{loop_annotation} -- Full refined 3D network with updated \code{loop_type}
#'     (e.g., eP-P) and two target gene columns:
#'     \itemize{
#'       \item \code{Active_Target_Genes}: Expression-filtered active-only targets (no fallback).
#'       \item \code{Putative_Target_Genes}: Display column; may include linear nearest-gene
#'         fallback when \code{Active_Target_Genes} is empty.
#'     }
#'     Refinement status columns:
#'     \code{Has_Active_Target}, \code{Retained_In_Functional_Network}, and
#'     \code{Refinement_Action} (\code{"retained_active_target"},
#'     \code{"reclassified_silent_anchor"}, \code{"expression_filtered_no_active_target"},
#'     or \code{"structural_only_no_active_target"}).
#'     All structural loops are preserved; filter on \code{Retained_In_Functional_Network}
#'     for the high-confidence active subset.
#'   \item \code{anchor_loci_annotation} -- Filtered non-redundant anchor-locus annotations with expressed targets.
#'   \item \code{anchor_annotation} -- Backward-compatible alias of \code{anchor_loci_annotation}.
#'   \item \code{promoter_centric_stats} -- Gene-level connectivity statistics.
#'   \item \code{distal_element_stats} -- Distal-element connectivity statistics.
#'   \item \code{target_annotation} -- Target features (peaks) with gene assignments.
#'     Key columns include:
#'     \itemize{
#'       \item \code{All_Loop_Connected_Genes}: All genes from loop-connected anchors (P/G types).
#'       \item \code{Regulated_promoter_genes}: Promoter genes supported by loop-anchor context.
#'       \item \code{Assigned_Target_Genes}: Promoter-first 3D assignment (prioritises P > G > E).
#'       \item \code{*_Filled} variants: Linear nearest-gene fallback when strict 3D assignments are empty.
#'       \item \code{Regulated_promoter_Evidence}: Provenance of \code{Regulated_promoter_genes}
#'         (e.g., \code{local_promoter_overlap}, \code{direct_opposite_promoter}).
#'         \strong{Read with} \code{Regulated_promoter_genes}; do not cross-reference
#'         with \code{Assigned_Target_Genes} or other columns.
#'       \item \code{Regulated_promoter_Fallback_Evidence}: Provenance of
#'         \code{Regulated_promoter_genes_Filled}.
#'         \strong{Read with} \code{Regulated_promoter_genes_Filled}; indicates
#'         which \code{*_Filled} column supplied the fallback gene.
#'     }
#'   \item \code{target_gene_links} -- Long-format peak-gene provenance table.
#'     Each row records one peak-gene linkage with full provenance.
#'     \strong{Read} \code{evidence}, \code{anchor_role}, and \code{gene_role}
#'     \strong{together as a group} -- they jointly describe how each gene was
#'     assigned to each peak; do not interpret any one column in isolation.
#'     \itemize{
#'       \item \code{input_id}, \code{loop_ID}, \code{anchor_id}: Identifiers.
#'       \item \code{gene}: Linked gene symbol.
#'       \item \code{gene_role}: \code{"promoter"}, \code{"gene_body"}, or \code{"linear_annotation"}.
#'       \item \code{source}: \code{"loop_anchor"} (3D-derived) or \code{"linear_annotation"} (nearest gene).
#'       \item \code{evidence}: Provenance label --
#'         \code{"local_promoter_overlap"} (peak overlaps anchor promoter),
#'         \code{"direct_opposite_promoter"} (opposite anchor is promoter),
#'         \code{"gene_body_context"} (gene body linkage),
#'         \code{"expanded_promoter_loop"} (via ego-network expansion),
#'         \code{"linear_annotation"} (direct nearest gene),
#'         or \code{"linear_fallback"} (filled when 3D assignment was empty).
#'       \item \code{anchor_role}: \code{"local_anchor"}, \code{"opposite_anchor"},
#'         \code{"expanded_anchor"}, or \code{"linear_annotation"}.
#'       \item \code{used_as_fallback}: Logical. \code{TRUE} when this link was added
#'         via the \code{*_Filled} linear nearest-gene fallback mechanism.
#'       \item \code{in_regulated_promoter} through \code{in_assigned_target_filled}:
#'         Logical membership flags indicating which target annotation column(s)
#'         this gene appears in. A gene may appear in multiple columns simultaneously.
#'       \item (Refine only) \code{Mean_Expression}: Per-gene mean expression value.
#'       \item (Refine only) \code{Passes_Expression_Filter}: Logical. \code{TRUE} if
#'         \code{Mean_Expression >= threshold}.
#'     }
#'   \item \code{chromatin_validation} -- Data frame from
#'     \code{\link{validate_epeG_by_chromatin}} with confidence levels and
#'     evidence strings for each eP/eG anchor. \code{NULL} when
#'     \code{chromatin_beds} is not provided or is empty (default).
#'   \item \code{plots} -- Named list of ggplot objects (dumbbell, rose, karyotype).
#'   \item \code{plot_list} -- Backward-compatible alias of \code{plots}.
#'   \item \code{metadata} -- Internal metadata list (parameters, versions, database versions). Not intended for direct use.
#' }
#' If \code{write_output = TRUE}, also writes \code{_Refined_Results.xlsx} to \code{out_dir}.
#' The workbook contains a \emph{Chromatin Validation} sheet (when
#' \code{chromatin_beds} is provided) and a \emph{Functional Loop Annotation}
#' sheet with only loops where \code{Retained_In_Functional_Network == TRUE}.
#'
#' @importFrom dplyr %>% filter group_by summarise ungroup mutate select rename left_join full_join arrange desc case_when rowwise coalesce any_of distinct pull
#' @importFrom ggplot2 ggplot aes geom_bar geom_segment geom_point geom_text scale_color_manual scale_fill_manual theme_minimal theme_void labs coord_polar xlim
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
#' # Safely load the pre-computed annotation result from RData
#' temp_env <- new.env()
#' load(rdata_path, envir = temp_env)
#' raw_annotation <- temp_env[[ls(temp_env)[1]]]
#' raw_annotation$loop_annotation <- head(raw_annotation$loop_annotation, 6)
#' raw_annotation$target_annotation <- head(raw_annotation$target_annotation, 3)
#' raw_annotation$promoter_centric_stats <- head(raw_annotation$promoter_centric_stats, 6)
#' raw_annotation$distal_element_stats <- head(raw_annotation$distal_element_stats, 6)
#'
#' res_reclassified <- refine_loop_anchors_by_expression(
#'     annotation_res = raw_annotation,
#'     expr_matrix_file = expr_path,
#'     sample_columns = "con1",
#'     threshold = 1.0,
#'     unit_type = "TPM",
#'     species = "hg38",
#'     out_dir = tempdir(),
#'     project_name = "Example_Reclassified",
#'     reclassify_by_expression = TRUE,
#'     write_output = FALSE,
#'     quiet = TRUE
#' )
#' print(table(res_reclassified$loop_annotation$loop_type))
refine_loop_anchors_by_expression <- function(
  annotation_res,
  expr_matrix_file = NULL,
  sample_columns = NULL,
  threshold = 1,
  unit_type = "TPM",
  species = "hg38",
  out_dir = "./results/filtered",
  project_name = "HiChIP",
  color_palette = "Paired",
  karyo_bin_size = 1e5,
  reclassify_by_expression = TRUE,
  hub_percentile = 0.95,
  chromatin_beds = list(),
  write_output = TRUE,
  quiet = FALSE
) {
    species <- match.arg(species, c("hg38", "hg19", "mm10", "mm9"))
    log_message <- function(...) {
        if (!quiet) message(...)
    }

    # --- 0. Setup ---
    if (!grepl("_Filtered$", project_name)) project_name <- paste0(project_name, "_Filtered")
    if (write_output && !dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    log_message(">>> [Refinement] Project Name: ", project_name)
    if (isTRUE(reclassify_by_expression)) {
        log_message(">>> NOTE: eP/eG labels indicate expression-filtered silent P/G anchors.")
        log_message(">>> These are NOT functional enhancer calls. Validate with ATAC-seq,")
        log_message(">>> H3K27ac, or other chromatin data before biological interpretation.")
    }

    # --- 1. Load & Validate ---
    refine_data <- .refine_load_validate_data(
        annotation_res, expr_matrix_file, sample_columns,
        threshold, unit_type, log_message
    )
    loop_df      <- refine_data$loop_df
    original_loop_df <- refine_data$original_loop_df
    clust_info   <- refine_data$clust_info
    bed_info     <- refine_data$bed_info
    target_gene_links <- refine_data$target_gene_links
    whitelist    <- refine_data$whitelist
    vals         <- refine_data$vals
    has_cluster_id <- refine_data$has_cluster_id
    upstream_promoter_stats <- refine_data$upstream_promoter_stats

    # --- 2. Update Anchors & Loops ---
    loop_df <- .refine_apply_anchor_updates(
        loop_df, whitelist, reclassify_by_expression, log_message
    )

    # --- 3. Stats Update ---
    log_message(">>> [Step 3] Updating Stats...")
    if (has_cluster_id) {
        agg_cluster <- loop_df %>%
            dplyr::filter(!is.na(cluster_id)) %>%
            dplyr::group_by(cluster_id) %>%
            dplyr::summarise(
                Cluster_All_Genes = extract_genes(Putative_Target_Genes),
                Cluster_Active_Target_Genes = extract_genes(Active_Target_Genes),
                .groups = "drop"
            )
        loop_df <- loop_df %>%
            dplyr::select(-any_of(c("Cluster_All_Genes", "Cluster_Active_Target_Genes"))) %>%
            dplyr::left_join(agg_cluster, by = "cluster_id")
        if (!is.null(clust_info) && "cluster_id" %in% colnames(clust_info)) {
            clust_info <- clust_info %>%
                dplyr::select(-any_of(c("Cluster_All_Genes", "Cluster_Active_Target_Genes"))) %>%
                dplyr::left_join(agg_cluster, by = "cluster_id")
        }
    }

    stats_res <- compute_refined_stats(
        loop_df = loop_df,
        upstream_promoter_stats = upstream_promoter_stats,
        vals = vals, threshold = threshold,
        hub_percentile = hub_percentile
    )
    promoter_centric_df <- stats_res$promoter_centric
    distal_element_df <- stats_res$distal_element

    # --- 4. Refine Target Annotations ---
    log_message(">>> [Step 4] Refining Target Annotations...")
    if (!is.null(bed_info)) {
        tgt_refined <- .refine_target_annotations(
            bed_info, loop_df, whitelist, target_gene_links, vals, threshold
        )
        bed_info <- tgt_refined$bed_info
        target_gene_links <- tgt_refined$target_gene_links
    }

    # --- 5. Visualization ---
    log_message(">>> [Step 5] Generating Visualizations (Returning plot objects)...")
    plot_list <- if (quiet) {
        .with_messages_silenced(
            .with_known_upstream_noise_suppressed(
                build_refinement_plots(
                    original_loop_df = original_loop_df,
                    loop_df = loop_df,
                    bed_info = bed_info,
                    whitelist = whitelist,
                    project_name = project_name,
                    karyo_bin_size = karyo_bin_size,
                    species = species,
                    color_palette = color_palette
                )
            )
        )
    } else {
        .with_known_upstream_noise_suppressed(
            build_refinement_plots(
                original_loop_df = original_loop_df,
                loop_df = loop_df,
                bed_info = bed_info,
                whitelist = whitelist,
                project_name = project_name,
                karyo_bin_size = karyo_bin_size,
                species = species,
                color_palette = color_palette
            )
        )
    }

    # --- 5.5 Orthogonal chromatin validation (optional) ---
    chromatin_validation <- NULL
    if (length(chromatin_beds) > 0) {
        log_message(">>> [Step 5.5] Orthogonal Chromatin Validation...")
        chromatin_validation <- validate_epeG_by_chromatin(
            annotation_res = list(loop_annotation = loop_df),
            chromatin_beds = chromatin_beds,
            anchor_gap = 200,
            anchor_min_overlap = 100,
            species = species,
            quiet = quiet
        )
    }

    # --- 6. Export ---
    if (write_output) {
        log_message(">>> [Step 6] Exporting Refined Results...")
        .refine_export_workbook(
            loop_df, clust_info, promoter_centric_df, distal_element_df,
            bed_info, target_gene_links, out_dir, project_name,
            chromatin_validation = chromatin_validation
        )
        log_message("    Excel saved.")
    }

    log_message("Refinement Complete.")

    out <- list(
        loop_annotation = loop_df,
        anchor_loci_annotation = clust_info,
        anchor_annotation = clust_info,
        promoter_centric_stats = promoter_centric_df,
        distal_element_stats = distal_element_df,
        target_annotation = bed_info,
        target_gene_links = target_gene_links,
        chromatin_validation = chromatin_validation,
        plots = plot_list,
        plot_list = plot_list,
        metadata = .build_looplook_metadata(
            fun = "refine_loop_anchors_by_expression",
            params = list(
                threshold = threshold,
                unit_type = unit_type,
                reclassify_by_expression = reclassify_by_expression,
                hub_percentile = hub_percentile,
                species = species
            ),
            genome_build = species,
            score_semantics = "expression-refined; eP/eG = expression-filtered silent P/G (NOT functional enhancers; validate with ATAC/H3K27ac)",
            database_versions = .record_database_versions(species)
        )
    )

    # Inherit and update anchor_state so expression -> chromatin ->
    # recompute_targets pipeline stays closed.
    anchor_state <- attr(annotation_res, "looplook_anchor_state", exact = TRUE)
    if (!is.null(anchor_state)) {
        attr(out, "looplook_anchor_state") <- .update_anchor_state_from_loop_df(
            anchor_state, loop_df
        )
    }

    return(out)
}

#' @title Chromatin-guided refinement of loop anchor classification
#'
#' @description
#' Refines loop anchor types (P, E, G, eP, eG) using chromatin mark evidence
#' (H3K4me1, H3K4me3, and optionally H3K27ac, ATAC, H3K27me3).  Anchors with
#' dual promoter-enhancer chromatin signatures are flagged as \code{"dual"}.
#' This function can be called before or after
#' \code{\link{refine_loop_anchors_by_expression}} -- the two refinements
#' address orthogonal biological questions (chromatin identity vs.
#' transcriptional activity) and do not depend on ordering.
#'
#' @details
#' \strong{Reclassification rules (minimum input: H3K4me1 + H3K4me3):}
#' \itemize{
#'   \item P + H3K4me1(+) and H3K4me3(+) -> \code{"dual"} (dual-function).
#'   \item P + H3K4me1(+) and H3K4me3(-) \emph{and} (H3K27ac(+) or ATAC(+))
#'         -> \code{"E"} (conservative: requires active-mark confirmation beyond
#'         H3K4me1 alone).
#'   \item eP/eG + gold_standard or high_confidence + active/primed enhancer
#'         chromatin -> kept (confirmed distal).
#'   \item eP + promoter_like (H3K4me3+, H3K27me3-) -> \code{"P"} (revert).
#'   \item eG + promoter_like -> \code{"G"} (revert).
#'   \item E + H3K4me1(+) and H3K4me3(+) -> \code{"dual"} (dual-function).
#'   \item E + H3K4me3(+) without H3K4me1 -> \code{"P"} (unannotated promoter).
#'   \item G + H3K4me1(+) and H3K4me3(+) -> \code{"dual"} (gene-body dual).
#'   \item G + H3K4me3(+) without H3K4me1 -> \code{"P"} (internal promoter).
#'   \item G + H3K4me1(+) and H3K4me3(-) \emph{and} (H3K27ac(+) or ATAC(+))
#'         -> \code{"E"} (conservative intronic enhancer).
#'   \item All other anchors: unchanged.
#' }
#'
#' \strong{Chromatin state inference:} Each anchor is also assigned a
#' \code{chromatin_state} based on its mark combination (highest-priority
#' match first):
#' \itemize{
#'   \item \code{"conflicting_marks"}: H3K27me3+ coexisting with any
#'         active mark (H3K4me1+, H3K27ac+, or H3K4me3+) —
#'         bivalent/poised/ambiguous chromatin; takes priority over
#'         enhancer-like and dual-like classifications.
#'   \item \code{"dual_like"}: H3K4me1+ H3K4me3+
#'   \item \code{"active_enhancer_like"}: H3K4me1+ H3K27ac+ ATAC+
#'   \item \code{"primed_enhancer_like"}: H3K4me1+ with H3K27ac+ or ATAC+
#'   \item \code{"weak_enhancer_like"}: H3K4me1+ only, or H3K27ac+ only
#'   \item \code{"promoter_like"}: H3K4me3+
#'   \item \code{"repressed"}: H3K27me3+ (no active marks)
#'   \item \code{"uncertain"}: marks tested but none present
#'   \item \code{"no_data"}: no marks tested
#' }
#'
#' \strong{Candidate anchor selection:} When \code{candidate_types = NULL}
#' (default), raw annotation from \code{\link{annotate_peaks_and_loops}} tests
#' \code{P}, \code{G}, and \code{E} anchors.  Expression-refined input tests
#' \code{eP} and \code{eG} anchors.
#'
#' After reclassification, \code{loop_type} is recomputed from the updated
#' anchor types.
#'
#' @param annotation_res List. Output from
#'   \code{\link{annotate_peaks_and_loops}} or
#'   \code{\link{refine_loop_anchors_by_expression}}.
#' @param chromatin_beds Named list of BED file paths. At minimum
#'   \code{H3K4me1} and \code{H3K4me3} are required for meaningful
#'   reclassification; \code{H3K27ac}, \code{ATAC}, and \code{H3K27me3}
#'   improve confidence but are optional.
#' @param anchor_gap Integer. Gap tolerance for mark overlap. Default \code{200}.
#' @param anchor_min_overlap Integer. Minimum overlap for mark overlap. Default \code{100}.
#' @param species Character. Genome assembly. Default \code{"hg38"}.
#' @param out_dir Character. Output directory for Excel export. Default \code{"./results/chromatin"}.
#' @param project_name Character. Project prefix. Default \code{"HiChIP"}.
#' @param candidate_types Character vector or \code{NULL}. Anchor types to
#'   validate and reclassify. \code{NULL} (default): auto-selects
#'   \code{c("eP","eG")} for refined input, \code{c("P","G","E")} for raw.
#'   Set explicitly to \code{c("P","G","E","eP","eG")} for full-range analysis.
#' @param recompute_targets Logical. If \code{TRUE}, re-run target gene
#'   assignment using updated anchor types. Requires the input
#'   \code{annotation_res} to contain the \code{looplook_anchor_state}
#'   attribute (present when using \code{\link{annotate_peaks_and_loops}}
#'   output). Default \code{FALSE}.
#' @param write_output Logical. Write Excel workbook. Default \code{TRUE}.
#' @param quiet Logical. Suppress messages. Default \code{FALSE}.
#'
#' @return An invisible named list with updated \code{loop_annotation},
#'   \code{anchor_loci_annotation}, \code{promoter_centric_stats},
#'   \code{distal_element_stats}, \code{chromatin_validation}, and
#'   \code{metadata}. When \code{recompute_targets = FALSE} (default),
#'   \code{target_annotation} and \code{target_gene_links} are \code{NULL}
#'   (pre-chromatin anchor types); downstream profiling should use
#'   \code{target_source = "loops"}. When \code{recompute_targets = TRUE},
#'   target links are rebuilt from chromatin-updated anchor states,
#'   producing chromatin-aware \code{target_annotation} and
#'   \code{target_gene_links}.  The input must carry the
#'   \code{looplook_anchor_state} attribute (present when using
#'   \code{\link{annotate_peaks_and_loops}} output).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' refined_chromatin <- refine_loop_anchors_by_chromatin(
#'   annotation_res = refined_expr,
#'   chromatin_beds = list(
#'     H3K4me1 = "H3K4me1_peaks.bed",
#'     H3K4me3 = "H3K4me3_peaks.bed",
#'     H3K27ac = "H3K27ac_peaks.bed"
#'   )
#' )
#' table(refined_chromatin$loop_annotation$loop_type)
#' }
refine_loop_anchors_by_chromatin <- function(
    annotation_res,
    chromatin_beds = list(),
    anchor_gap = 200,
    anchor_min_overlap = 100,
    species = "hg38",
    out_dir = "./results/chromatin",
    project_name = "HiChIP",
    candidate_types = NULL,
    recompute_targets = FALSE,
    write_output = TRUE,
    quiet = FALSE
) {
    species <- match.arg(species, c("hg38", "hg19", "mm10", "mm9"))
    if (!grepl("_Chromatin$", project_name)) project_name <- paste0(project_name, "_Chromatin")
    log_message <- function(...) { if (!quiet) message(...) }

    log_message(">>> [Chromatin Refinement] Project: ", project_name)
    if (length(chromatin_beds) == 0) {
        stop("chromatin_beds must contain at least H3K4me1 and H3K4me3 BED files.", call. = FALSE)
    }

    # --- 1. Load & validate ---
    loop_df <- annotation_res$loop_annotation
    if (is.null(loop_df)) stop("'annotation_res$loop_annotation' is missing.")
    known_marks <- c("H3K4me1", "H3K27ac", "ATAC", "H3K27me3", "H3K4me3")
    provided_marks <- intersect(names(chromatin_beds), known_marks)
    if (!all(c("H3K4me1", "H3K4me3") %in% provided_marks)) {
        stop("chromatin_beds must include at least 'H3K4me1' and 'H3K4me3'.", call. = FALSE)
    }

    # --- 2. Extract candidate anchors and overlap marks ---
    epeG_anchors <- .extract_epeG_anchors(annotation_res, log_message,
                                           candidate_types)
    if (inherits(epeG_anchors, "data.frame") && nrow(epeG_anchors) == 0) {
        warning("No P/G/eP/eG anchors found; returning input unchanged.", call. = FALSE)
        return(annotation_res)
    }

    mark_matrix <- .overlap_chromatin_marks(
        epeG_anchors, chromatin_beds, provided_marks, known_marks,
        anchor_gap, anchor_min_overlap, log_message
    )

    validation <- .assign_chromatin_confidence(
        epeG_anchors, mark_matrix, provided_marks, known_marks
    )

    # --- 3. Reclassify anchors based on chromatin evidence ---
    log_message(">>> Reclassifying anchors by chromatin evidence...")
    reclass_map <- .chromatin_reclassify(validation)
    log_message(sprintf("    Reclassified %d / %d anchors", sum(reclass_map$changed), nrow(reclass_map)))

    # --- 4. Apply reclassification to loop_annotation ---
    loop_df <- .chromatin_update_loops(loop_df, reclass_map, validation)

    # --- 5. Recompute loop_type and stats ---
    log_message("    Recomputing loop types and stats...")
    loop_df <- .chromatin_recompute_loop_types(loop_df)

    # --- 6. Summary ---
    log_message("--- Chromatin Refinement Summary ---")
    tab <- table(validation$confidence, useNA = "ifany")
    for (lvl in names(tab)) {
        log_message(sprintf("  %-16s: %d anchors", lvl, tab[lvl]))
    }
    n_promoter_like <- sum(grepl("promoter_like", validation$evidence, fixed = TRUE))
    if (n_promoter_like > 0) {
        log_message(sprintf("  %-16s: %d anchors", "promoter-like", n_promoter_like))
    }
    # Add positional/final type layers to validation output
    validation$positional_type <- validation$anchor_type  # pre-reclass
    final_types <- setNames(reclass_map$new_type, reclass_map$anchor_id)
    validation$final_type <- final_types[validation$anchor_id]
    validation$final_type[is.na(validation$final_type)] <- validation$positional_type[is.na(validation$final_type)]

    log_message(sprintf("  Reclassified      : %d anchors", sum(reclass_map$changed)))
    log_message("--- End Chromatin Refinement ---")

    # --- 7. Recompute stats after reclassification ---
    log_message("    Recomputing connectivity stats...")
    new_promoter_stats <- .compute_raw_promoter_stats(loop_df)
    # Build a dummy expression vector: all genes considered "active" since
    # chromatin refinement does not assess expression status.
    all_genes <- unique(new_promoter_stats$Gene)
    dummy_vals <- setNames(rep(1, length(all_genes)), all_genes)
    promoter_centric_df <- .build_promoter_centric_df(
        new_promoter_stats, annotation_res$promoter_centric_stats,
        vals = dummy_vals, threshold = 0,
        hub_percentile = 0.95
    )
    distal_element_df <- tryCatch(
        .build_distal_element_df(loop_df, 0.95),
        error = function(e) NULL
    )

    # --- 8. Recompute target links if requested ---
    ta  <- NULL; tgl <- NULL
    if (isTRUE(recompute_targets)) {
        log_message(">>> Recomputing target gene links with chromatin-aware anchors...")
        recomputed <- .recompute_targets_from_anchor_state(
            list(loop_annotation = loop_df,
                 chromatin_validation = validation), annotation_res, reclass_map
        )
        ta  <- recomputed$target_annotation
        tgl <- recomputed$target_gene_links
        log_message(sprintf("    Updated %d target gene links.", nrow(tgl)))
    }

    # --- 9. Build output ---
    # When recompute_targets = FALSE, target_annotation / target_gene_links
    # are NULL (pre-chromatin types).  Use profile_target_genes(target_source = "loops")
    # for chromatin-aware profiling.
    qc_summary <- data.frame(
        n_candidate_anchors = nrow(validation),
        n_reclassified      = sum(reclass_map$changed),
        n_gold_standard     = sum(validation$confidence == "gold_standard"),
        n_high_confidence   = sum(validation$confidence == "high_confidence"),
        n_supported         = sum(validation$confidence == "supported"),
        n_weak              = sum(validation$confidence == "weak"),
        n_dual              = sum(reclass_map$new_type == "dual"),
        n_promoter_like     = n_promoter_like,
        stringsAsFactors = FALSE
    )

    out <- list(
        loop_annotation = loop_df,
        chromatin_validation = validation,
        qc_summary = qc_summary,
        anchor_loci_annotation = annotation_res$anchor_loci_annotation,
        anchor_annotation = annotation_res$anchor_annotation,
        promoter_centric_stats = promoter_centric_df,
        distal_element_stats = distal_element_df,
        target_annotation = ta,
        target_gene_links = tgl,
        metadata = .build_looplook_metadata(
            fun = "refine_loop_anchors_by_chromatin",
            params = list(
                marks_provided = provided_marks,
                anchor_gap = anchor_gap,
                anchor_min_overlap = anchor_min_overlap,
                species = species
            ),
            genome_build = species,
            score_semantics = "chromatin-guided reclassification; dual = promoter-enhancer dual-function",
            database_versions = .record_database_versions(species)
        )
    )

    # Inherit and update anchor_state so expression -> chromatin ->
    # recompute_targets pipeline stays closed end-to-end.
    anchor_state <- attr(annotation_res, "looplook_anchor_state", exact = TRUE)
    if (!is.null(anchor_state)) {
        attr(out, "looplook_anchor_state") <- .update_anchor_state_from_loop_df(
            anchor_state, loop_df
        )
    }

    if (write_output) {
        log_message(">>> Exporting to Excel...")
        wb <- openxlsx::createWorkbook()
        openxlsx::addWorksheet(wb, "Loop Annotation")
        openxlsx::writeData(wb, "Loop Annotation",
            loop_df %>% dplyr::select(-any_of(c("a1_id", "a2_id"))))
        openxlsx::addWorksheet(wb, "Chromatin Validation")
        openxlsx::writeData(wb, "Chromatin Validation", validation)
        if (!is.null(ta) && nrow(ta) > 0) {
            openxlsx::addWorksheet(wb, "Target Annotation")
            openxlsx::writeData(wb, "Target Annotation", ta)
        }
        if (!is.null(tgl) && nrow(tgl) > 0) {
            openxlsx::addWorksheet(wb, "Target Gene Links")
            openxlsx::writeData(wb, "Target Gene Links", tgl)
        }
        tryCatch(
            openxlsx::saveWorkbook(wb,
                file.path(out_dir, paste0(project_name, "_Chromatin_Results.xlsx")),
                overwrite = TRUE),
            error = function(e)
                warning("Failed to save Excel: ", conditionMessage(e), call. = FALSE)
        )
        log_message("    Excel saved.")
    }

    log_message("Chromatin Refinement Complete.")
    return(invisible(out))
}

# --- Internal helpers for refine_loop_anchors_by_chromatin ---

#' Internal: Build reclassification map from chromatin validation
#' @keywords internal
#' @noRd
.chromatin_reclassify <- function(validation) {
    n <- nrow(validation)
    out <- data.frame(
        anchor_id        = validation$anchor_id,
        positional_type  = validation$anchor_type,  # original ChIPseeker type
        old_type         = validation$anchor_type,
        new_type         = validation$anchor_type,
        # chromatin_state is inferred per-row in the loop below
        chromatin_state  = NA_character_,
        changed          = FALSE,
        stringsAsFactors = FALSE
    )

    for (i in seq_len(n)) {
        conf <- as.character(validation$confidence[i])
        old  <- out$old_type[i]
        is_promoter_like <- grepl("promoter_like", validation$evidence[i], fixed = TRUE)
        h3k4me1_pos <- isTRUE(validation$H3K4me1[i])
        h3k4me3_pos <- isTRUE(validation$H3K4me3[i])
        h3k27ac_pos  <- isTRUE(validation$H3K27ac[i])
        h3k27me3_pos <- isTRUE(validation$H3K27me3[i])

        # Infer chromatin state from mark combination.
        # H3K27me3 (repressive) coexisting with any active mark takes
        # highest priority: bivalent/poised/conflicting chromatin should
        # not be called enhancer-like or dual-like.
        if (h3k27me3_pos && (h3k4me1_pos || h3k27ac_pos || h3k4me3_pos)) {
            out$chromatin_state[i] <- "conflicting_marks"
        } else if (h3k4me1_pos && h3k4me3_pos) {
            out$chromatin_state[i] <- "dual_like"
        } else if (h3k4me1_pos && h3k27ac_pos &&
                   isTRUE(validation$ATAC[i])) {
            out$chromatin_state[i] <- "active_enhancer_like"
        } else if (h3k4me1_pos && (h3k27ac_pos ||
                   isTRUE(validation$ATAC[i]))) {
            out$chromatin_state[i] <- "primed_enhancer_like"
        } else if (h3k4me1_pos || h3k27ac_pos) {
            out$chromatin_state[i] <- "weak_enhancer_like"
        } else if (h3k4me3_pos) {
            out$chromatin_state[i] <- "promoter_like"
        } else if (h3k27me3_pos) {
            out$chromatin_state[i] <- "repressed"
        } else if (any(!is.na(c(validation$H3K4me1[i], validation$H3K4me3[i],
                                validation$H3K27ac[i])))) {
            out$chromatin_state[i] <- "uncertain"
        } else {
            out$chromatin_state[i] <- "no_data"
        }

        # P + H3K4me1+ H3K4me3+ -> dual
        if (old == "P" && h3k4me1_pos && h3k4me3_pos) {
            out$new_type[i] <- "dual"
            out$changed[i]  <- TRUE
        }
        # P + H3K4me1+ H3K4me3- + (H3K27ac+ or ATAC+) -> E (conservative)
        # Requires active mark confirmation beyond H3K4me1 alone.
        else if (old == "P" && h3k4me1_pos &&
                 isTRUE(!is.na(validation$H3K4me3[i]) && !validation$H3K4me3[i]) &&
                 (h3k27ac_pos || isTRUE(validation$ATAC[i]))) {
            out$new_type[i] <- "E"
            out$changed[i]  <- TRUE
        }
        # eP/eG with dual_like chromatin -> reclassify as dual
        else if (old %in% c("eP", "eG") &&
                 out$chromatin_state[i] == "dual_like") {
            out$new_type[i] <- "dual"
            out$changed[i]  <- TRUE
        }
        # eP with promoter_like -> revert to P
        else if (old == "eP" && out$chromatin_state[i] == "promoter_like") {
            out$new_type[i] <- "P"
            out$changed[i]  <- TRUE
        }
        # eG with promoter_like -> revert to G
        else if (old == "eG" && out$chromatin_state[i] == "promoter_like") {
            out$new_type[i] <- "G"
            out$changed[i]  <- TRUE
        }
        # eP/eG confirmed distal (gold/high, active/primed enhancer) -> keep
        else if (old %in% c("eP", "eG") &&
                 conf %in% c("gold_standard", "high_confidence") &&
                 out$chromatin_state[i] %in% c("active_enhancer_like",
                                               "primed_enhancer_like")) {
            # keep -- confirmed distal regulatory element, no conflicting promoter marks
        }
        # eP + promoter_like -> revert to P
        else if (old == "eP" && is_promoter_like) {
            out$new_type[i] <- "P"
            out$changed[i]  <- TRUE
        }
        # eG + promoter_like -> revert to G
        else if (old == "eG" && is_promoter_like) {
            out$new_type[i] <- "G"
            out$changed[i]  <- TRUE
        }
        # E + H3K4me1+ H3K4me3+ -> dual (unannotated dual-function locus)
        else if (old == "E" && h3k4me1_pos && h3k4me3_pos) {
            out$new_type[i] <- "dual"
            out$changed[i]  <- TRUE
        }
        # E + H3K4me3+ (no H3K4me1) -> P (unannotated promoter / ncRNA gene)
        else if (old == "E" && h3k4me3_pos &&
                 isTRUE(!is.na(validation$H3K4me1[i]) && !validation$H3K4me1[i])) {
            out$new_type[i] <- "P"
            out$changed[i]  <- TRUE
        }
        # G + H3K4me1+ H3K4me3+ -> dual (gene body dual-function element)
        else if (old == "G" && h3k4me1_pos && h3k4me3_pos) {
            out$new_type[i] <- "dual"
            out$changed[i]  <- TRUE
        }
        # G + H3K4me3+ (no H3K4me1) -> P (internal promoter)
        else if (old == "G" && h3k4me3_pos &&
                 isTRUE(!is.na(validation$H3K4me1[i]) && !validation$H3K4me1[i])) {
            out$new_type[i] <- "P"
            out$changed[i]  <- TRUE
        }
        # G + H3K4me1+ H3K4me3- + (H3K27ac+ or ATAC+) -> E (conservative intronic enhancer)
        else if (old == "G" && h3k4me1_pos &&
                 isTRUE(!is.na(validation$H3K4me3[i]) && !validation$H3K4me3[i]) &&
                 (h3k27ac_pos || isTRUE(validation$ATAC[i]))) {
            out$new_type[i] <- "E"
            out$changed[i]  <- TRUE
        }
        # default: unchanged
    }
    out
}

#' Internal: Update loop_annotation anchor types from reclassification map
#' @keywords internal
#' @noRd
.chromatin_update_loops <- function(loop_df, reclass_map, validation) {
    type_map <- setNames(reclass_map$new_type, reclass_map$anchor_id)

    # Update anchor1
    if ("a1_id" %in% colnames(loop_df)) {
        idx1 <- match(loop_df$a1_id, reclass_map$anchor_id)
        hits1 <- which(!is.na(idx1) & reclass_map$changed[idx1])
        if (length(hits1) > 0) {
            loop_df$anchor1_type[hits1] <- type_map[loop_df$a1_id[hits1]]
        }
    }
    # Update anchor2
    if ("a2_id" %in% colnames(loop_df)) {
        idx2 <- match(loop_df$a2_id, reclass_map$anchor_id)
        hits2 <- which(!is.na(idx2) & reclass_map$changed[idx2])
        if (length(hits2) > 0) {
            loop_df$anchor2_type[hits2] <- type_map[loop_df$a2_id[hits2]]
        }
    }
    loop_df
}

#' Internal: Recompute loop_type from updated anchor types
#' @keywords internal
#' @noRd
.chromatin_recompute_loop_types <- function(loop_df) {
    if (!all(c("anchor1_type", "anchor2_type") %in% colnames(loop_df))) {
        return(loop_df)
    }
    loop_df$loop_type <- vapply(seq_len(nrow(loop_df)), function(i) {
        .loop_type_code(loop_df$anchor1_type[i], loop_df$anchor2_type[i])
    }, character(1))
    loop_df
}

#' Internal: Recompute target gene links from updated anchor types
#'
#' Uses the saved anchor state from annotate_peaks_and_loops to rerun
#' peak-to-gene mapping after chromatin reclassification, without
#' re-running ChIPseeker annotation of the target peaks.
#'
#' @param annotation_res The chromatin-refined annotation result.
#' @param original_res The original annotation (with anchor state).
#' @param reclass_map Reclassification data frame from .chromatin_reclassify.
#' @return Updated target_annotation and target_gene_links.
#' @keywords internal
#' @noRd
.recompute_targets_from_anchor_state <- function(annotation_res, original_res,
                                                   reclass_map) {
    anchor_state <- attr(original_res, "looplook_anchor_state")
    if (is.null(anchor_state)) {
        warning("No anchor state found in original annotation; cannot recompute target links.", call. = FALSE)
        return(list(target_annotation = NULL, target_gene_links = NULL))
    }

    map_info         <- anchor_state$map_info
    anchor_topo_map  <- anchor_state$anchor_topo_map
    gr_anchors       <- anchor_state$gr_anchors
    ego_list_target  <- anchor_state$ego_list_target

    if (is.null(map_info) || nrow(map_info) == 0) {
        return(list(target_annotation = NULL, target_gene_links = NULL))
    }

    # Update map_info type_code for reclassified anchors
    type_updates <- setNames(reclass_map$new_type, reclass_map$anchor_id)
    matched <- match(map_info$anchor_id, names(type_updates))
    update_idx <- which(!is.na(matched) & reclass_map$changed[matched])
    if (length(update_idx) > 0) {
        map_info$type_code[update_idx] <- type_updates[map_info$anchor_id[update_idx]]
    }

    # Rebuild topology maps with updated types
    map_info$SYMBOL <- trimws(map_info$SYMBOL)
    valid_pg <- map_info %>%
        dplyr::filter(
            (.is_promoter_like(type_code) | .is_gene_body_like(type_code)) &
            !is.na(SYMBOL) & SYMBOL != ""
        )
    lookup_pg_symbol <- setNames(valid_pg$SYMBOL, valid_pg$anchor_id)
    lookup_pg_type   <- setNames(valid_pg$type_code, valid_pg$anchor_id)

    lookup_p_symbol <- map_info %>%
        dplyr::filter(type_code %in% c("P", "dual") & !is.na(SYMBOL) & SYMBOL != "") %>%
        { setNames(.$SYMBOL, .$anchor_id) }

    nodes_in_graph <- names(ego_list_target)
    anchor_topo_map <- data.frame(
        anchor_id = nodes_in_graph,
        tgt_genes_pg = vapply(ego_list_target, function(x)
            .ids_to_genes_simple(names(x), lookup_pg_symbol), character(1)),
        tgt_genes_p  = vapply(ego_list_target, function(x)
            .ids_to_genes_simple(names(x), lookup_p_symbol), character(1)),
        tgt_genes_prio = vapply(ego_list_target, function(x)
            .ids_to_genes_priority(names(x), lookup_pg_symbol, lookup_pg_type),
            character(1)),
        stringsAsFactors = FALSE
    )
    anchor_topo_map[is.na(anchor_topo_map)] <- NA_character_

    # Rebuild target_gene_links using stored peak-to-anchor hits and
    # chromatin-updated topology.  If no hit_df was saved (e.g. no target_bed
    # was provided), fall back to linear-only links with a warning.
    hit_df <- anchor_state$target_hit_df
    if (is.null(hit_df) || nrow(hit_df) == 0) {
        warning("No stored peak-to-anchor hit_df found; target_gene_links will be linear-only.",
                call. = FALSE)
        hit_df <- data.frame(qid = integer(0), anchor_id = character(0),
                              stringsAsFactors = FALSE)
    }
    new_links <- .build_target_gene_links(
        hit_df    = hit_df,
        bed_info  = original_res$target_annotation,
        loop_annotation_final = annotation_res$loop_annotation,
        map_info  = map_info,
        ego_list_target = ego_list_target
    )

    # Annotate with chromatin-aware fields
    if (nrow(new_links) > 0 && "anchor_id" %in% colnames(new_links)) {
        type_before <- setNames(reclass_map$old_type, reclass_map$anchor_id)
        type_after  <- setNames(reclass_map$new_type, reclass_map$anchor_id)
        changed_map <- setNames(reclass_map$changed, reclass_map$anchor_id)
        new_links$anchor_type_before_chromatin <- type_before[new_links$anchor_id]
        new_links$anchor_type_after_chromatin  <- type_after[new_links$anchor_id]
        new_links$chromatin_target_action <- ifelse(
            !is.na(new_links$anchor_type_before_chromatin) &
            new_links$anchor_type_before_chromatin != new_links$anchor_type_after_chromatin,
            paste0(new_links$anchor_type_before_chromatin, "->",
                   new_links$anchor_type_after_chromatin),
            "unchanged"
        )
        # Attach chromatin validation per anchor (one-to-many from validation)
        if (!is.null(annotation_res$chromatin_validation)) {
            cv <- annotation_res$chromatin_validation
            conf_map <- setNames(as.character(cv$confidence), cv$anchor_id)
            evid_map <- setNames(cv$evidence, cv$anchor_id)
            new_links$chromatin_confidence <- conf_map[new_links$anchor_id]
            new_links$chromatin_evidence    <- evid_map[new_links$anchor_id]
        }
    }

    # Rebuild target_annotation columns from the chromatin-aware links,
    # aggregating from evidence/gene_role rather than from pre-existing
    # membership flags (which are all FALSE at this point).
    new_bed_info <- original_res$target_annotation
    if (!is.null(new_bed_info)) {
        # Clear old target assignment columns unconditionally so that
        # an empty new_links result does not silently retain stale values.
        strict_cols <- c(
            "All_Loop_Connected_Genes",
            "Regulated_promoter_genes",
            "Assigned_Target_Genes",
            "All_Loop_Connected_Genes_Filled",
            "Regulated_promoter_genes_Filled",
            "Assigned_Target_Genes_Filled",
            "Regulated_promoter_Evidence",
            "Regulated_promoter_Fallback_Evidence"
        )
        for (col in base::intersect(strict_cols, colnames(new_bed_info))) {
            new_bed_info[[col]] <- NA_character_
        }
    }

    if (!is.null(new_bed_info) && !is.null(new_links) && nrow(new_links) > 0) {
        loop_links <- new_links %>%
            dplyr::filter(source == "loop_anchor", !is.na(gene), gene != "")

        # All_Loop_Connected_Genes: all genes from loop anchors
        all_agg <- loop_links %>%
            dplyr::group_by(input_id) %>%
            dplyr::summarise(
                All_Loop_Connected_Genes = paste(sort(unique(gene)), collapse = ";"),
                .groups = "drop"
            )

        # Regulated_promoter_genes: promoter genes supported by loop context
        promoter_agg <- loop_links %>%
            dplyr::filter(gene_role == "promoter") %>%
            dplyr::group_by(input_id) %>%
            dplyr::summarise(
                Regulated_promoter_genes = paste(sort(unique(gene)), collapse = ";"),
                .groups = "drop"
            )

        # Assigned_Target_Genes: promoter-first priority
        assigned_agg <- loop_links %>%
            dplyr::mutate(priority = dplyr::case_when(
                gene_role == "promoter" ~ 1L,
                evidence %in% c("direct_opposite_promoter",
                                "local_promoter_overlap",
                                "expanded_promoter_loop") ~ 1L,
                gene_role == "gene_body" ~ 2L,
                TRUE ~ 3L
            )) %>%
            dplyr::group_by(input_id) %>%
            dplyr::filter(priority == min(priority, na.rm = TRUE)) %>%
            dplyr::summarise(
                Assigned_Target_Genes = paste(sort(unique(gene)), collapse = ";"),
                .groups = "drop"
            )

        # Join back into bed_info
        new_bed_info <- new_bed_info %>%
            dplyr::select(-dplyr::any_of(c(
                "All_Loop_Connected_Genes",
                "Regulated_promoter_genes",
                "Assigned_Target_Genes"
            ))) %>%
            dplyr::left_join(all_agg, by = "input_id") %>%
            dplyr::left_join(promoter_agg, by = "input_id") %>%
            dplyr::left_join(assigned_agg, by = "input_id")

        # Rebuild Filled columns and Evidence from updated target_gene_links
        evidence_df <- .summarise_regulated_promoter_evidence(new_links)
        new_bed_info <- new_bed_info %>%
            dplyr::select(-dplyr::any_of(c("Regulated_promoter_Evidence"))) %>%
            dplyr::left_join(evidence_df, by = "input_id")
        new_bed_info$Regulated_promoter_Evidence <- ifelse(
            is.na(new_bed_info$Regulated_promoter_Evidence) |
                new_bed_info$Regulated_promoter_Evidence == "",
            "none",
            new_bed_info$Regulated_promoter_Evidence
        )

        # Rebuild Filled fallback columns
        fallback_col <- .target_linear_gene_column(new_bed_info)
        fallback_vec <- if (!is.null(fallback_col)) new_bed_info[[fallback_col]] else
            rep(NA_character_, nrow(new_bed_info))
        ann_vec <- if ("annotation" %in% colnames(new_bed_info)) {
            new_bed_info$annotation
        } else {
            rep(NA_character_, nrow(new_bed_info))
        }
        fallback_evidence <- .fallback_evidence_from_annotation(ann_vec)

        new_bed_info <- new_bed_info %>% dplyr::mutate(
            All_Loop_Connected_Genes_Filled = .fill_target_gene_fallback(
                All_Loop_Connected_Genes, fallback_vec),
            Regulated_promoter_genes_Filled = .fill_target_gene_fallback(
                Regulated_promoter_genes, fallback_vec),
            Assigned_Target_Genes_Filled = .fill_target_gene_fallback(
                Assigned_Target_Genes, fallback_vec),
            Regulated_promoter_Fallback_Evidence =
                dplyr::case_when(
                    !is.na(Regulated_promoter_genes) & Regulated_promoter_genes != "" ~ "none",
                    !is.na(Regulated_promoter_genes_Filled) & Regulated_promoter_genes_Filled != "" ~
                        fallback_evidence,
                    TRUE ~ "none"
                )
        )
    }

    # Mark membership flags against the NEW target_annotation, not the old one.
    # This ensures chromatin-updated anchor types are reflected.
    if (!is.null(new_bed_info) && nrow(new_links) > 0) {
        new_links <- .mark_target_gene_link_membership(new_links, new_bed_info)
    }

    list(target_annotation = new_bed_info, target_gene_links = new_links)
}

#' Internal: Build Refinement Dumbbell Comparison Plot
#'
#' Compares loop-type counts before vs. after expression-aware filtering.
#'
#' @param original_loop_df Loop annotation before refinement.
#' @param loop_df Loop annotation after refinement.
#' @param project_name Character. Project prefix for the plot title.
#' @return A \code{ggplot} object.
#' @keywords internal
#' @noRd
.build_dumbbell_plot <- function(original_loop_df, loop_df, project_name) {
    df_orig <- original_loop_df %>%
        dplyr::group_by(loop_type) %>%
        dplyr::summarise(Original = dplyr::n(), .groups = "drop")
    df_filt <- loop_df %>%
        dplyr::group_by(loop_type) %>%
        dplyr::summarise(Refined = dplyr::n(), .groups = "drop")
    df_dumbbell <- dplyr::full_join(df_orig, df_filt, by = "loop_type") %>%
        dplyr::mutate(
            Original = ifelse(is.na(Original), 0, Original),
            Refined = ifelse(is.na(Refined), 0, Refined),
            is_e_type = grepl("e", loop_type)
        ) %>%
        dplyr::arrange(is_e_type, dplyr::desc(Original))
    df_dumbbell$loop_type <- factor(df_dumbbell$loop_type,
        levels = rev(df_dumbbell$loop_type)
    )
    df_long <- df_dumbbell %>%
        tidyr::pivot_longer(
            cols = c("Original", "Refined"),
            names_to = "Source", values_to = "Count"
        )

    ggplot2::ggplot() +
        ggplot2::geom_segment(
            data = df_dumbbell,
            ggplot2::aes(y = loop_type, yend = loop_type, x = Original, xend = Refined),
            color = "#b2b2b2", linewidth = 0.8
        ) +
        ggplot2::geom_point(
            data = df_long,
            ggplot2::aes(x = Count, y = loop_type, color = Source), size = 3
        ) +
        ggplot2::scale_color_manual(values = c("Original" = "#999999", "Refined" = "#E69F00")) +
        ggplot2::theme_minimal() +
        ggplot2::labs(
            title = paste0(project_name, ": Refinement Effect (Dumbbell)"),
            x = "Number of Loops", y = "Loop Type"
        ) +
        ggplot2::theme(
            plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
            legend.position = "top"
        )
}

#' Internal: Build Refinement Target-Loop Donut Plot
#'
#' Donut chart of loop-type distribution among target-connected refined loops.
#'
#' @param bed_info Target annotation data frame (or NULL).
#' @param loop_df Refined loop annotation.
#' @param custom_colors Named color vector keyed by loop_type.
#' @param project_name Character. Project prefix for the plot title.
#' @return A \code{ggplot} object, or \code{NULL} if no target-connected loops.
#' @keywords internal
#' @noRd
.build_refinement_donut <- function(bed_info, loop_df, custom_colors, project_name) {
    if (is.null(bed_info)) {
        return(NULL)
    }
    gr_bed <- GenomicRanges::makeGRangesFromDataFrame(bed_info,
        keep.extra.columns = TRUE
    )
    if (!"cluster_id" %in% colnames(loop_df) || all(is.na(loop_df$cluster_id))) {
        return(NULL)
    }
    active_anc <- dplyr::bind_rows(
        loop_df %>% dplyr::select(chr = chr1, start = start1, end = end1, cluster_id),
        loop_df %>% dplyr::select(chr = chr2, start = start2, end = end2, cluster_id)
    ) %>%
        dplyr::filter(!is.na(cluster_id)) %>%
        dplyr::distinct()
    if (nrow(active_anc) == 0) {
        return(NULL)
    }
    gr_anc <- GenomicRanges::makeGRangesFromDataFrame(active_anc,
        keep.extra.columns = TRUE
    )
    hits <- GenomicRanges::findOverlaps(gr_bed, gr_anc)
    if (length(hits) == 0) {
        return(NULL)
    }
    hit_ids <- unique(gr_anc$cluster_id[S4Vectors::subjectHits(hits)])
    tgt_loops <- loop_df %>% dplyr::filter(cluster_id %in% hit_ids)
    if (nrow(tgt_loops) == 0) {
        return(NULL)
    }
    donut_data <- tgt_loops %>%
        dplyr::group_by(loop_type) %>%
        dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
        dplyr::mutate(
            fraction = count / sum(count),
            legend_label = paste0(
                loop_type, " (n=", count, ", ",
                round(fraction * 100, 1), "%)"
            ),
            plot_label = loop_type, is_lower_e = grepl("^e", loop_type)
        ) %>%
        dplyr::arrange(is_lower_e, dplyr::desc(count)) %>%
        dplyr::mutate(loop_type = factor(loop_type, levels = loop_type))
    ggplot2::ggplot(
        donut_data,
        ggplot2::aes(x = 2, y = count, fill = loop_type)
    ) +
        ggplot2::geom_bar(stat = "identity", width = 1, color = "white") +
        ggplot2::coord_polar(theta = "y") +
        ggplot2::xlim(0.5, 2.9) +
        ggplot2::geom_text(ggplot2::aes(x = 2.8, label = plot_label),
            position = ggplot2::position_stack(vjust = 0.5), size = 3
        ) +
        ggplot2::scale_fill_manual(
            values = custom_colors,
            labels = setNames(donut_data$legend_label, donut_data$loop_type)
        ) +
        ggplot2::theme_void() +
        ggplot2::labs(title = paste0(project_name, ": Refined Structural Loops at Target Regions")) +
        ggplot2::theme(
            legend.position = "right",
            plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
        )
}

#' Internal: Build Multi-Omics Sankey Tracking Plot
#'
#' Three-tier Sankey diagram tracing genomic features through loop
#' connection status to expression activity.
#'
#' @param bed_info Target annotation data frame.
#' @param whitelist Character vector of active gene symbols.
#' @param project_name Character. Project prefix for the plot title.
#' @return A \code{networkD3} sankey widget, or \code{NULL} if requirements
#'   are not met.
#' @keywords internal
#' @noRd
.build_sankey_plot <- function(bed_info, whitelist, project_name) {
    if (is.null(bed_info)) {
        return(NULL)
    }
    if (!"SANKEY_RAW_GENES" %in% colnames(bed_info)) {
        return(NULL)
    }
    if (!requireNamespace("networkD3", quietly = TRUE) || !requireNamespace("htmlwidgets", quietly = TRUE)) {
        return(NULL)
    }

    get_label_mapping <- function(vec, total_targets) {
        tab <- table(vec)
        label_map <- list()
        for (i in seq_along(tab)) {
            name <- names(tab)[i]
            ct <- as.integer(tab[i])
            pct <- round(ct / total_targets * 100, 1)
            label_map[[name]] <- sprintf("%s (n=%d, %.1f%%)", name, ct, pct)
        }
        label_map
    }
    check_status_strict <- function(s, wl) {
        if (is.na(s) || s == "") {
            return("Inactive")
        }
        gs <- trimws(unlist(strsplit(as.character(s), ";")))
        gs <- gs[gs != ""]
        if (length(gs) == 0) {
            return("Inactive")
        }
        if (any(gs %in% wl)) {
            return("Active")
        }
        return("Inactive")
    }

    bed_info$L1_Raw <- vapply(bed_info$annotation, function(a) {
        if (grepl("Intergenic", a, ignore.case = TRUE)) {
            return("Distal Intergenic")
        }
        if (grepl("Promoter", a, ignore.case = TRUE)) {
            return("Promoter")
        }
        if (grepl("Exon", a, ignore.case = TRUE)) {
            return("Exon")
        }
        if (grepl("Intron", a, ignore.case = TRUE)) {
            return("Intron")
        }
        return("Others")
    }, FUN.VALUE = character(1))

    bed_info$L2_Raw <- ifelse(
        !is.na(bed_info$Linked_Loop_IDs) & bed_info$Linked_Loop_IDs != "",
        "Connected", "Unconnected"
    )

    bed_info$L3_Raw <- vapply(bed_info$SANKEY_RAW_GENES, check_status_strict,
        wl = whitelist, FUN.VALUE = character(1)
    )

    sankey_df <- bed_info %>%
        dplyr::filter(
            !is.na(L3_Raw) & L3_Raw != "" &
                !is.na(L1_Raw) & L1_Raw != "" &
                !is.na(L2_Raw) & L2_Raw != ""
        )
    if (nrow(sankey_df) == 0) {
        return(NULL)
    }

    total_targets <- nrow(sankey_df)
    label_map_l1 <- get_label_mapping(sankey_df$L1_Raw, total_targets)
    label_map_l2 <- get_label_mapping(sankey_df$L2_Raw, total_targets)
    label_map_l3 <- get_label_mapping(sankey_df$L3_Raw, total_targets)

    sankey_df$L1_Label <- label_map_l1[sankey_df$L1_Raw]
    sankey_df$L2_Label <- label_map_l2[sankey_df$L2_Raw]
    sankey_df$L3_Label <- label_map_l3[sankey_df$L3_Raw]

    sankey_df <- sankey_df %>%
        dplyr::filter(
            !is.na(.data$L1_Label) & .data$L1_Label != "" &
                !is.na(.data$L2_Label) & .data$L2_Label != "" &
                !is.na(.data$L3_Label) & .data$L3_Label != ""
        )
    if (nrow(sankey_df) == 0) {
        return(NULL)
    }

    nodes <- unique(c(
        unlist(sankey_df$L1_Label, use.names = FALSE),
        unlist(sankey_df$L2_Label, use.names = FALSE),
        unlist(sankey_df$L3_Label, use.names = FALSE)
    ))
    if (length(nodes) < 2) {
        return(NULL)
    }
    nodes <- data.frame(name = nodes, stringsAsFactors = FALSE)

    get_idx <- function(label) match(label, nodes$name) - 1
    links <- data.frame(
        source = get_idx(sankey_df$L1_Label),
        target = get_idx(sankey_df$L2_Label),
        value = 1,
        stringsAsFactors = FALSE
    )
    links2 <- data.frame(
        source = get_idx(sankey_df$L2_Label),
        target = get_idx(sankey_df$L3_Label),
        value = 1,
        stringsAsFactors = FALSE
    )
    links <- rbind(links, links2)
    links <- links %>%
        dplyr::group_by(.data$source, .data$target) %>%
        dplyr::summarise(value = dplyr::n(), .groups = "drop")

    sankey_colors <- get_colors(nrow(nodes), "Paired")
    color_scale <- paste0('d3.scaleOrdinal().range(["',
        paste(sankey_colors, collapse = '","'), '"])')
    sn <- networkD3::sankeyNetwork(
        Links = links, Nodes = nodes,
        Source = "source", Target = "target",
        Value = "value", NodeID = "name",
        units = "TWh", fontSize = 12, nodeWidth = 30,
        colourScale = networkD3::JS(color_scale),
        iterations = 0
    )

    sn$sizingPolicy$defaultWidth <- "100%"
    sn$sizingPolicy$defaultHeight <- "450px"
    sn <- htmlwidgets::onRender(sn, sprintf('
	function(el, x) {
	  var svg = d3.select(el).select("svg");
	  function createValidID(name) {
	    if (!name) return "unknown";
	    return name.replace(/[^a-zA-Z0-9-]/g, "_");
	  }
	  svg.selectAll(".link").each(function(d) {
	    var gradientID = "gradient-" + createValidID(d.source.name) +
	      "-" + createValidID(d.target.name);
	    if (svg.select("#" + gradientID).empty()) {
	      var gradient = svg.append("defs")
	        .append("linearGradient")
	        .attr("id", gradientID)
	        .attr("gradientUnits", "userSpaceOnUse")
	        .attr("x1", d.source.x + d.source.dx / 2)
	        .attr("y1", d.source.y + d.source.dy / 2)
	        .attr("x2", d.target.x + d.target.dx / 2)
	        .attr("y2", d.target.y + d.target.dy / 2);
	      var sourceColor = d3.select(el).selectAll(".node")
	        .filter(function(node) { return node.name === d.source.name; })
	        .select("rect").style("fill");
	      var targetColor = d3.select(el).selectAll(".node")
	        .filter(function(node) { return node.name === d.target.name; })
	        .select("rect").style("fill");
	      gradient.append("stop").attr("offset", "0%%")
	        .attr("stop-color", sourceColor);
	      gradient.append("stop").attr("offset", "100%%")
	        .attr("stop-color", targetColor);
	    }
	    d3.select(this).style("stroke", "url(#" + gradientID + ")")
	      .style("stroke-opacity", 0.6)
	      .style("stroke-width", function(d) { return Math.max(2, d.width); });
	  });
	  svg.selectAll(".node rect")
	    .style("stroke", "#333333")
	    .style("stroke-width", "1px");
	  svg.selectAll("text")
	    .style("font-size", "12px")
	    .style("font-weight", "bold");
	}
	'))
    sn
}

#' Internal: Build Refinement Karyotype Heatmaps
#'
#' Generates \code{Refined_Karyo_Active} and \code{Refined_Karyo_TargetGenes}
#' karyotype heatmaps from refined loop and target annotation data.
#'
#' @param loop_df Refined loop annotation.
#' @param bed_info Target annotation data frame (or NULL).
#' @param species Character. Genome assembly.
#' @param karyo_bin_size Integer. Bin size for karyotype heatmaps.
#' @param red_palette Character vector. Red color palette.
#' @param purple_palette Character vector. Purple color palette.
#' @return A named list of karyotype grob objects (may be empty).
#' @keywords internal
#' @noRd
.build_refinement_karyotypes <- function(
  loop_df, bed_info, species,
  karyo_bin_size, red_palette, purple_palette
) {
    plot_list <- list()
    txdb_pkg <- tryCatch(species_txdb_pkg(species), error = function(e) NULL)
    org_db <- tryCatch(species_orgdb_pkg(species), error = function(e) NULL)
    if (is.null(txdb_pkg) || is.null(org_db) ||
        !requireNamespace(txdb_pkg, quietly = TRUE) ||
        !requireNamespace(org_db, quietly = TRUE)) {
        return(plot_list)
    }
    txdb_obj <- utils::getFromNamespace(txdb_pkg, txdb_pkg)

    # Determine which gene sets need karyotypes
    need_active <- "Active_Target_Genes" %in% colnames(loop_df)
    need_target <- !is.null(bed_info) &&
        "Assigned_Target_Genes_Filled" %in% colnames(bed_info)
    if (!need_active && !need_target) {
        return(plot_list)
    }

    g_active <- if (need_active) {
        clean_gene_names(loop_df$Active_Target_Genes, ";")
    } else {
        character(0)
    }
    g_target <- if (need_target) {
        clean_gene_names(bed_info$Assigned_Target_Genes_Filled, ";")
    } else {
        character(0)
    }

    g_all <- unique(c(g_active, g_target))
    if (length(g_all) == 0) {
        return(plot_list)
    }

    # Load gene catalogue once and annotate with SYMBOL
    all_genes_gr <- .with_known_upstream_noise_suppressed(
        GenomicFeatures::genes(txdb_obj)
    )
    map <- .map_txdb_gene_ids(
        gene_ids = .extract_txdb_gene_ids(all_genes_gr),
        org_db = org_db,
        columns = "SYMBOL",
        context = "build_refinement_plots karyotype",
        warn = FALSE
    )
    S4Vectors::mcols(all_genes_gr)$SYMBOL <- map$SYMBOL[match(
        .extract_txdb_gene_ids(all_genes_gr), map$gene_id
    )]

    if (length(g_active) > 0) {
        target_genes_gr <- all_genes_gr[
            S4Vectors::mcols(all_genes_gr)$SYMBOL %in% g_active
        ]
        plot_list$Refined_Karyo_Active <- draw_karyo_heatmap_internal(
            target_genes_gr,
            "Refined Active Genes", karyo_bin_size, 0.99, txdb_obj, species,
            "Genes",
            custom_colors = red_palette
        )
    }

    if (length(g_target) > 0) {
        target_genes_gr <- all_genes_gr[
            S4Vectors::mcols(all_genes_gr)$SYMBOL %in% g_target
        ]
        plot_list$Refined_Karyo_TargetGenes <- draw_karyo_heatmap_internal(
            target_genes_gr,
            "Refined Target Genes", karyo_bin_size, 0.99, txdb_obj, species,
            "Genes",
            custom_colors = purple_palette
        )
    }
    plot_list
}

#' Internal: Build Refinement Rose Plot
#'
#' Polar bar chart of refined loop-type distribution.
#'
#' @param loop_df Refined loop annotation.
#' @param custom_colors Named color vector keyed by loop_type.
#' @param project_name Character. Project prefix for the plot title.
#' @return A \code{ggplot} object.
#' @keywords internal
#' @noRd
.build_rose_plot <- function(loop_df, custom_colors, project_name) {
    rose_data <- loop_df %>%
        dplyr::group_by(loop_type) %>%
        dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
        dplyr::mutate(
            fraction = count / sum(count),
            legend_label = paste0(loop_type, " (n=", count, ", ",
                round(fraction * 100, 1), "%)"),
            is_lower_e = grepl("^e", loop_type)
        ) %>%
        dplyr::arrange(dplyr::desc(count))

    plot_order <- rose_data$loop_type
    rose_data$loop_type <- factor(rose_data$loop_type, levels = plot_order)
    legend_order <- rose_data %>%
        dplyr::arrange(is_lower_e, dplyr::desc(count)) %>%
        dplyr::pull(loop_type)

    ggplot2::ggplot(rose_data, ggplot2::aes(x = loop_type, y = count, fill = loop_type)) +
        ggplot2::geom_bar(stat = "identity", width = 1, color = "white") +
        ggplot2::coord_polar(theta = "x") +
        ggplot2::scale_fill_manual(
            values = custom_colors,
            labels = setNames(rose_data$legend_label, as.character(rose_data$loop_type)),
            breaks = legend_order,
            name = "Loop Type"
        ) +
        ggplot2::theme_void() +
        ggplot2::labs(title = paste0(project_name, ": Structural Loop Types After Reclassification"), subtitle = "Full refined network (all loops). For active subset, filter on Retained_In_Functional_Network.") +
        ggplot2::theme(
            plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14),
            legend.position = "right",
            legend.text = ggplot2::element_text(size = 10)
        )
}
#' Internal: Build Refinement Visualization Suite
#'
#' Generates diagnostic plots for expression-aware refinement: dumbbell
#' comparison, donut, Sankey, karyotype heatmaps, and rose plot.
#'
#' @param original_loop_df Loop annotation before refinement.
#' @param loop_df Loop annotation after refinement.
#' @param bed_info Target annotation data frame (optional).
#' @param whitelist Character vector of active gene symbols.
#' @param project_name Character. Project prefix for plot titles.
#' @param karyo_bin_size Integer. Bin size for karyotype heatmaps.
#' @param species Character. Genome assembly.
#' @return A named list of ggplot / htmlwidget / grob objects.
#' @keywords internal
#' @noRd
build_refinement_plots <- function(
  original_loop_df, loop_df, bed_info,
  whitelist, project_name, karyo_bin_size, species,
  color_palette = "Paired"
) {
    red_palette <- c("#FFFFFF", "#FFFFCC", "#FFEDA0", "#FED976", "#FEB24C", "#FD8D3C", "#FC4E2A", "#E31A1C", "#BD0026", "#800026", "#000000")
    purple_palette <- c("#FFFFFF", "#F3E5F5", "#E1BEE7", "#BA68C8", "#9C27B0", "#7B1FA2", "#4A148C", "#000000")

    # Assign colours by descending loop-type frequency (not alphabetical)
    type_counts <- loop_df %>%
        dplyr::group_by(loop_type) %>%
        dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
        dplyr::arrange(dplyr::desc(n))
    all_types <- type_counts$loop_type
    custom_colors <- get_colors(length(all_types), color_palette)
    names(custom_colors) <- all_types

    plot_list <- list()

    plot_list$Comparison_Dumbbell <- .build_dumbbell_plot(
        original_loop_df, loop_df, project_name
    )

    plot_list$Target_Loop_Donut <- .build_refinement_donut(
        bed_info, loop_df, custom_colors, project_name
    )

    plot_list$Target_Sankey <- .build_sankey_plot(
        bed_info, whitelist, project_name
    )

    karyo_plots <- .build_refinement_karyotypes(
        loop_df, bed_info, species, karyo_bin_size,
        red_palette, purple_palette
    )
    for (nm in names(karyo_plots)) plot_list[[nm]] <- karyo_plots[[nm]]

    plot_list$Rose <- .build_rose_plot(loop_df, custom_colors, project_name)

    return(plot_list)
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
#' @noRd
format_annotation_columns <- function(df) {
    if ("annotation" %in% colnames(df)) {
        df <- df %>%
            dplyr::rename(detail_anno = annotation) %>%
            dplyr::mutate(annotation = gsub(" \\(.*", "", detail_anno)) %>%
            dplyr::relocate(annotation, .before = detail_anno)
    }
    return(df)
}

#' @title Orthogonal validation of eP/eG reclassification by chromatin marks
#'
#' @description
#' Validates the expression-aware eP/eG reclassification produced by
#' \code{\link{refine_loop_anchors_by_expression}} using orthogonal chromatin
#' data (ATAC-seq, ChIP-seq). Each eP/eG anchor is tested for overlap with
#' user-supplied mark BED files and assigned a confidence level based on the
#' ENCODE active-enhancer signature.
#'
#' @details
#' \strong{Confidence levels (ENCODE active-enhancer criteria):}
#' \itemize{
#'   \item \code{"gold_standard"}: All five marks align with the canonical
#'     active-enhancer signature: H3K4me1(+), H3K27ac(+), ATAC(+),
#'     H3K4me3(-), H3K27me3(-). Requires all five marks to be provided.
#'   \item \code{"high_confidence"}: H3K4me1(+) and H3K27ac(+), or
#'     H3K4me1(+) and ATAC(+). At least two supporting marks present.
#'   \item \code{"supported"}: At least one enhancer-associated mark
#'     (H3K4me1, H3K27ac, or ATAC) is present.
#'   \item \code{"weak"}: Only exclusion marks are informative
#'     (H3K27me3(-) or H3K4me3(-)) without positive enhancer evidence.
#'   \item \code{"uncertain"}: No chromatin data provided or no overlaps
#'     detected.
#' }
#'
#' \strong{Mark semantics:}
#' \itemize{
#'   \item \emph{Positive marks} (presence = supporting evidence):
#'     \code{H3K4me1}, \code{H3K27ac}, \code{ATAC}.
#'   \item \emph{Negative marks} (absence = supporting evidence):
#'     \code{H3K27me3}, \code{H3K4me3}.
#' }
#'
#' @param annotation_res List. Output from
#'   \code{\link{annotate_peaks_and_loops}} or
#'   \code{\link{refine_loop_anchors_by_expression}}.
#'   When the refined output is provided, all anchors classified as
#'   \code{eP} or \code{eG} are validated. When the raw annotation
#'   is provided, all \code{P}, \code{G}, and \code{E} anchors are
#'   evaluated (chromatin evidence may reclassify E to P or dual).
#'   Anchors are tested for overlap with chromatin mark BED files and
#'   assigned confidence levels based on ENCODE active-enhancer criteria.
#'   annotated as P or G are tested.
#' @param chromatin_beds Named list of BED file paths. Names must be
#'   mark identifiers: any of \code{"H3K4me1"}, \code{"H3K27ac"},
#'   \code{"ATAC"}, \code{"H3K27me3"}, \code{"H3K4me3"}.
#'   Additional names are ignored with a warning.
#' @param anchor_gap Integer. Maximum gap (bp) between an anchor and a
#'   chromatin peak for them to be considered overlapping. Default: \code{200}.
#' @param anchor_min_overlap Integer. Minimum overlap (bp) required.
#'   Default: \code{100}.
#' @param candidate_types Character vector or \code{NULL}. Candidate anchor
#'   types to validate. \code{NULL} (default): \code{c("eP","eG")} for refined
#'   input, \code{c("P","G","E")} for raw.  Set to \code{c("P","G","E","eP","eG")}
#'   to cover all positional categories.
#' @param species Character. Genome assembly for seqlevel harmonization.
#'   Default: \code{"hg38"}.
#' @param quiet Logical. Suppress progress messages. Default: \code{FALSE}.
#'
#' @return A data frame with one row per candidate anchor (eP/eG when the input
#'   is refined, P/G when the input is raw from \code{annotate_peaks_and_loops}):
#' \describe{
#'   \item{\code{anchor_id}}{Anchor identifier.}
#'   \item{\code{chr}, \code{start}, \code{end}}{Anchor coordinates.}
#'   \item{\code{anchor_type}}{Original type (P or G) before reclassification.}
#'   \item{\code{anchor_gene}}{Gene symbol(s) at this anchor.}
#'   \item{\code{cluster_id}}{Loop cluster identifier.}
#'   \item{\code{H3K4me1}, \code{H3K27ac}, \code{ATAC}}{Logical. TRUE if
#'     the anchor overlaps the corresponding positive mark.}
#'   \item{\code{H3K27me3}, \code{H3K4me3}}{Logical. TRUE if the anchor
#'     overlaps the corresponding negative mark.}
#'   \item{\code{confidence}}{Factor with levels:
#'     \code{"gold_standard"}, \code{"high_confidence"}, \code{"supported"},
#'     \code{"weak"}, \code{"uncertain"}.}
#'   \item{\code{evidence}}{Human-readable summary of which marks supported
#'     the classification.}
#' }
#'
#' @importFrom GenomicRanges GRanges findOverlaps makeGRangesFromDataFrame
#' @importFrom S4Vectors queryHits
#' @export
#'
#' @examples
#' # Load pre-computed annotation result
#' rdata_path <- system.file("extdata", "analysis_results.RData",
#'                           package = "looplook")
#' temp_env <- new.env()
#' load(rdata_path, envir = temp_env)
#' raw_annotation <- temp_env[[ls(temp_env)[1]]]
#'
#' # Create dummy chromatin BED files for demonstration
#' bed_dir <- tempdir()
#' writeLines("chr6\t10410000\t10413000", file.path(bed_dir, "H3K4me1.bed"))
#' writeLines("chr6\t10411000\t10414000", file.path(bed_dir, "H3K27ac.bed"))
#'
#' # Run validation (using raw annotation; pass refined for eP/eG only)
#' result <- validate_epeG_by_chromatin(
#'     annotation_res = raw_annotation,
#'     chromatin_beds = list(
#'         H3K4me1 = file.path(bed_dir, "H3K4me1.bed"),
#'         H3K27ac = file.path(bed_dir, "H3K27ac.bed")
#'     ),
#'     quiet = TRUE
#' )
#' table(result$confidence)
#'
validate_epeG_by_chromatin <- function(
    annotation_res,
    chromatin_beds = list(),
    anchor_gap = 200,
    anchor_min_overlap = 100,
    candidate_types = NULL,
    species = "hg38",
    quiet = FALSE
) {
    log_message <- function(...) {
        if (!quiet) message(...)
    }

    # ---- 1. Identify candidate anchors ----
    epeG_anchors <- .extract_epeG_anchors(annotation_res, log_message,
                                           candidate_types)

    # ---- 2. Validate chromatin_beds input ----
    known_marks <- c("H3K4me1", "H3K27ac", "ATAC", "H3K27me3", "H3K4me3")
    if (length(chromatin_beds) == 0) {
        warning("No chromatin_beds provided; all anchors classified as 'uncertain'.",
                call. = FALSE)
    }
    unknown <- setdiff(names(chromatin_beds), known_marks)
    if (length(unknown) > 0) {
        warning("Unknown mark names ignored: ",
                paste(unknown, collapse = ", "),
                ". Expected: ", paste(known_marks, collapse = ", "),
                call. = FALSE)
    }
    provided_marks <- intersect(names(chromatin_beds), known_marks)

    # ---- 3+4. Overlap marks with anchors ----
    mark_matrix <- .overlap_chromatin_marks(
        epeG_anchors, chromatin_beds, provided_marks, known_marks,
        anchor_gap, anchor_min_overlap, log_message
    )

    # ---- 5. Assign confidence levels and evidence ----
    result <- .assign_chromatin_confidence(
        epeG_anchors, mark_matrix, provided_marks, known_marks
    )

    # ---- 6. Summary ----
    log_message("--- Validation Summary ---")
    tab <- table(result$confidence, useNA = "ifany")
    for (lvl in names(tab)) {
        log_message(sprintf("  %-16s: %d anchors", lvl, tab[lvl]))
    }
    n_promoter_like <- sum(grepl("promoter_like", result$evidence, fixed = TRUE))
    if (n_promoter_like > 0) {
        log_message(sprintf("  %-16s: %d anchors (H3K4me3+ H3K27me3- active marks; may reflect promoter signal)", "promoter-like", n_promoter_like))
    }
    log_message("--- End Validation ---")

    out <- result %>%
        dplyr::arrange(confidence, anchor_id) %>%
        dplyr::select(
            anchor_id, chr, start, end,
            anchor_type, anchor_gene, cluster_id,
            H3K4me1, H3K27ac, ATAC, H3K27me3, H3K4me3,
            confidence, evidence
        )
    attr(out, "looplook_metadata") <- list(
        package = "looplook",
        version = as.character(utils::packageVersion("looplook")),
        function_name = "validate_epeG_by_chromatin",
        call_time = Sys.time(),
        parameters = list(
            anchor_gap = anchor_gap,
            anchor_min_overlap = anchor_min_overlap,
            species = species,
            marks_provided = names(chromatin_beds),
            mark_files = lapply(chromatin_beds, basename)
        ),
        r_version = R.version.string,
        platform = R.version$platform
    )
    out
}

#' Internal: Empty validation result data frame
#' @keywords internal
#' @noRd
.empty_validation_result <- function() {
    data.frame(
        anchor_id = character(), chr = character(),
        start = integer(), end = integer(),
        anchor_type = character(), anchor_gene = character(),
        cluster_id = character(),
        H3K4me1 = logical(), H3K27ac = logical(), ATAC = logical(),
        H3K27me3 = logical(), H3K4me3 = logical(),
        confidence = factor(levels = c("gold_standard", "high_confidence",
                                       "supported", "weak", "uncertain")),
        evidence = character(),
        stringsAsFactors = FALSE
    )
}

#' Internal: Extract eP/eG (or P/G) anchors from annotation results
#' @keywords internal
#' @noRd
.extract_epeG_anchors <- function(annotation_res, log_message, candidate_types = NULL) {
    loop_df <- annotation_res$loop_annotation
    if (is.null(loop_df)) stop("'annotation_res$loop_annotation' is missing.")
    has_refined <- "Retained_In_Functional_Network" %in% colnames(loop_df)
    has_anchor_ids <- all(c("a1_id", "a2_id") %in% colnames(loop_df))

    if (has_anchor_ids) {
        anchor_map <- dplyr::bind_rows(
            loop_df %>% dplyr::select(anchor_id = a1_id, chr = chr1,
                start = start1, end = end1, anchor_type = anchor1_type,
                anchor_gene = anchor1_gene, cluster_id),
            loop_df %>% dplyr::select(anchor_id = a2_id, chr = chr2,
                start = start2, end = end2, anchor_type = anchor2_type,
                anchor_gene = anchor2_gene, cluster_id)
        ) %>% dplyr::distinct()
    } else {
        anchor_map <- dplyr::bind_rows(
            loop_df %>% dplyr::mutate(anchor_id = paste(chr1, start1, end1,
                sep = "_")) %>% dplyr::select(anchor_id, chr = chr1,
                start = start1, end = end1, anchor_type = anchor1_type,
                anchor_gene = anchor1_gene, cluster_id),
            loop_df %>% dplyr::mutate(anchor_id = paste(chr2, start2, end2,
                sep = "_")) %>% dplyr::select(anchor_id, chr = chr2,
                start = start2, end = end2, anchor_type = anchor2_type,
                anchor_gene = anchor2_gene, cluster_id)
        ) %>% dplyr::distinct()
    }

    type_label <- if (has_refined) "eP/eG" else "P/G"
    # Use all positional categories (P/G/E) for raw input -- chromatin
    # evidence may reclassify E anchors to P or dual as well.
    type_filter <- if (!is.null(candidate_types)) {
        candidate_types
    } else if (has_refined) {
        c("eP", "eG")
    } else {
        c("P", "G", "E")
    }
    epeG_anchors <- anchor_map %>% dplyr::filter(anchor_type %in% type_filter)
    if (nrow(epeG_anchors) == 0) {
        message("No ", type_label, " anchors found. Returning empty result.")
        return(data.frame(anchor_id = character(), chr = character(),
            start = integer(), end = integer(), anchor_type = character(),
            anchor_gene = character(), cluster_id = character(),
            stringsAsFactors = FALSE))
    }
    log_message(sprintf("Validating %d %s anchors...", nrow(epeG_anchors), type_label))
    epeG_anchors
}

#' Internal: Overlap chromatin mark BEDs with anchor GRanges
#' @keywords internal
#' @noRd
.overlap_chromatin_marks <- function(epeG_anchors, chromatin_beds, provided_marks,
                                      known_marks, anchor_gap, anchor_min_overlap,
                                      log_message) {
    if (inherits(epeG_anchors, "data.frame") && nrow(epeG_anchors) == 0)
        return(data.frame())
    gr_anchors <- .with_known_upstream_noise_suppressed(
        GenomicRanges::makeGRangesFromDataFrame(epeG_anchors, keep.extra.columns = TRUE)
    )
    mark_matrix <- as.data.frame(
        matrix(NA, nrow = nrow(epeG_anchors), ncol = length(known_marks),
               dimnames = list(NULL, known_marks))
    )
    for (mark in provided_marks) {
        bed_path <- chromatin_beds[[mark]]
        if (is.null(bed_path) || !file.exists(bed_path)) {
            warning("BED file not found for ", mark, ": ", bed_path, call. = FALSE)
            next
        }
        log_message(sprintf("  Overlapping with %s ...", mark))
        mark_gr <- read_simple_bed(bed_path, quiet = TRUE)
        mark_gr <- .harmonize_seqlevels(mark_gr, gr_anchors, mark)
        hits <- GenomicRanges::findOverlaps(gr_anchors, mark_gr, maxgap = anchor_gap)
        if (anchor_min_overlap > 1L && length(hits) > 0) {
            q_gr <- gr_anchors[S4Vectors::queryHits(hits)]
            s_gr <- mark_gr[S4Vectors::subjectHits(hits)]
            overlap_w <- GenomicRanges::width(GenomicRanges::pintersect(q_gr, s_gr))
            hits <- hits[overlap_w >= anchor_min_overlap]
        }
        hit_anchors <- unique(S4Vectors::queryHits(hits))
        mark_matrix[[mark]][hit_anchors] <- TRUE
        log_message(sprintf("    %d / %d anchors overlap %s peaks",
            length(hit_anchors), nrow(epeG_anchors), mark))
    }
    # For marks that were actually provided, non-overlapping anchors are
    # tested and absent (FALSE), not untested (NA).  Marks that were never
    # provided stay NA.
    for (mark in provided_marks) {
        na_idx <- is.na(mark_matrix[[mark]])
        mark_matrix[[mark]][na_idx] <- FALSE
    }
    mark_matrix
}

#' Internal: Assign chromatin confidence levels and evidence strings
#' @keywords internal
#' @noRd
.assign_chromatin_confidence <- function(epeG_anchors, mark_matrix,
                                          provided_marks, known_marks) {
    if (nrow(mark_matrix) == 0) return(.empty_validation_result())
    result <- epeG_anchors
    result$H3K4me1  <- mark_matrix$H3K4me1
    result$H3K27ac  <- mark_matrix$H3K27ac
    result$ATAC     <- mark_matrix$ATAC
    result$H3K27me3 <- mark_matrix$H3K27me3
    result$H3K4me3  <- mark_matrix$H3K4me3

    has_positive <- vapply(seq_len(nrow(result)), function(i) {
        isTRUE(result$H3K4me1[i]) || isTRUE(result$H3K27ac[i]) || isTRUE(result$ATAC[i])
    }, logical(1))
    has_negative_excl <- vapply(seq_len(nrow(result)), function(i) {
        .is_absent <- function(x) !is.na(x) && !x
        .is_absent(result$H3K27me3[i]) && .is_absent(result$H3K4me3[i])
    }, logical(1))
    all_five <- length(provided_marks) == 5
    .is_absent <- function(x) !is.na(x) && !x

    result$confidence <- vapply(seq_len(nrow(result)), function(i) {
        h3k4me1_t <- !is.na(result$H3K4me1[i]); h3k27ac_t <- !is.na(result$H3K27ac[i])
        atac_t    <- !is.na(result$ATAC[i])
        h3k27me3_t <- !is.na(result$H3K27me3[i]); h3k4me3_t <- !is.na(result$H3K4me3[i])
        neg <- c(
            if (.is_absent(result$H3K27me3[i])) "H3K27me3-",
            if (.is_absent(result$H3K4me3[i])) "H3K4me3-"
        )
        if (all_five && h3k4me1_t && h3k27ac_t && atac_t && h3k27me3_t && h3k4me3_t &&
            result$H3K4me1[i] && result$H3K27ac[i] && result$ATAC[i] &&
            .is_absent(result$H3K27me3[i]) && .is_absent(result$H3K4me3[i]))
            return("gold_standard")
        if (h3k4me1_t && result$H3K4me1[i] &&
            ((h3k27ac_t && result$H3K27ac[i]) || (atac_t && result$ATAC[i])))
            return("high_confidence")
        if (has_positive[i]) return("supported")
        if (length(neg) > 0) return("weak")
        return("uncertain")
    }, character(1))

    result$confidence <- factor(result$confidence,
        levels = c("gold_standard", "high_confidence", "supported", "weak", "uncertain"))

    result$evidence <- vapply(seq_len(nrow(result)), function(i) {
        parts <- c()
        if (isTRUE(result$H3K4me1[i]))  parts <- c(parts, "H3K4me1+")
        if (isTRUE(result$H3K27ac[i]))  parts <- c(parts, "H3K27ac+")
        if (isTRUE(result$ATAC[i]))     parts <- c(parts, "ATAC+")
        if (isTRUE(result$H3K27me3[i])) parts <- c(parts, "H3K27me3+")
        if (isTRUE(result$H3K4me3[i]))  parts <- c(parts, "H3K4me3+")
        if (isTRUE(!is.na(result$H3K27me3[i]) && !result$H3K27me3[i]))
            parts <- c(parts, "H3K27me3-")
        if (isTRUE(!is.na(result$H3K4me3[i]) && !result$H3K4me3[i]))
            parts <- c(parts, "H3K4me3-")
        # Flag: active marks present + H3K4me3 tested & positive +
        # H3K27me3 tested & absent -> chromatin looks more like active promoter
        # than distal enhancer.  Kept as 'supported' but tagged for filtering.
        if (result$confidence[i] == "supported" &&
            isTRUE(result$H3K4me3[i]) &&
            isTRUE(!is.na(result$H3K27me3[i]) && !result$H3K27me3[i])) {
            parts <- c(parts, "promoter_like")
        }
        if (length(parts) == 0) return("no_data")
        paste(parts, collapse = "; ")
    }, character(1))
    result
}
