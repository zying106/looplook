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
  # Try the "package" attribute first -- most reliable for installed pkgs.
  pkg_attr <- attr(x, "package")
  if (!is.null(pkg_attr) && length(pkg_attr) == 1L && nzchar(pkg_attr) &&
    pkg_attr != "AnnotationDbi") {
    if (requireNamespace(pkg_attr, quietly = TRUE)) {
      return(pkg_attr)
    }
  }
  # Fall back to db path inference, but only if the inferred name
  # corresponds to a real installed package.
  db_path <- tryCatch(AnnotationDbi::dbfile(x), error = function(e) NULL)
  if (!is.null(db_path) && length(db_path) == 1L && nzchar(db_path)) {
    candidate <- basename(dirname(dirname(db_path)))
    if (requireNamespace(candidate, quietly = TRUE)) {
      return(candidate)
    }
  }
  # Cannot reliably infer a package name -- this is likely a custom
  # SQLite file or a TxDb built from an uninstalled source.
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
    .require_pkg(arg, desc, "stop")
    return(list(obj = utils::getFromNamespace(arg, arg), pkg = arg))
  }
  pkg <- if (type == "txdb") species_txdb_pkg(species) else species_orgdb_pkg(species)
  if (is.null(pkg)) {
    stop(
      desc, " package for species '", species, "' is not built-in. ",
      "Pass your own ", if (type == "txdb") "TxDb" else "OrgDb",
      " object or package name via the ",
      if (type == "txdb") "txdb" else "org_db", " parameter. ",
      "See ?annotate_peaks_and_loops for details."
    )
  }
  .require_pkg(pkg, desc, "stop")
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

  # Fallback: unrecognised ChIPseeker annotation -- do NOT assume distal.
  # Return "Unknown" so these anchors are excluded from P/G/E target
  # assignment rather than silently treated as enhancers.
  if (!is.na(anno_str) && nzchar(anno_str)) {
    warning("Unrecognised ChIPseeker annotation type '", anno_str,
      "', returning 'Unknown'. This anchor will be excluded ",
      "from P/G/E-based target assignment. ",
      "This may indicate a ChIPseeker version change; ",
      "review after package upgrades.",
      call. = FALSE
    )
  }
  "Unknown"
}

#' Internal: Extract Loop Locus Genes
#' @return A semicolon-delimited character string of gene symbols from P/G anchors, or an empty string.
#' @keywords internal
#' @noRd
.loop_locus_genes <- function(t1, t2, s1, s2) {
  genes <- c()
  if (!is.na(t1) && (.is_target_promoter_like(t1) || .is_target_gene_body_like(t1))) genes <- c(genes, s1)
  if (!is.na(t2) && (.is_target_promoter_like(t2) || .is_target_gene_body_like(t2))) genes <- c(genes, s2)
  extract_genes(genes)
}

#' Internal: Build Loop Type Code
#' @details Anchors are sorted with case-insensitive alphabetical ordering
#'   so that lowercase prefixes (eP, eG) sort before uppercase (E, G, P)
#'   when the base letter is the same. This gives a canonical ordering
#'   where expression-filtered anchors appear first (e.g. \code{"eP-P"}
#'   rather than \code{"P-eP"}).
#' @return A two-letter code string (e.g. \code{"E-P"}, \code{"P-P"}) with anchors sorted alphabetically, or \code{"Unknown"}.
#' @keywords internal
#' @noRd
.loop_type_code <- function(t1, t2) {
  if (length(t1) == 0 || length(t2) == 0 || is.na(t1) || is.na(t2)) {
    return("Unknown")
  }
  paste(c(t1, t2)[order(tolower(c(t1, t2)), c(t1, t2), method = "radix")], collapse = "-")
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
  extract_genes(lookup[valid])
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
# --- .annotate_target_bed helpers ---
#' @keywords internal
#' Internal: Read, annotate, and resolve gene conflicts for target BED

#'
#' Reads a BED file, converts to GRanges, runs ChIPseeker annotation,
#' resolves gene conflicts, and harmonises seqlevels with loop anchors.
#'
#' @return A list(gr_bed, bed_info, n_peaks), or NULL if the BED is empty.
#' @keywords internal
#' @noRd
#' @noRd
.prepare_target_bed_granges <- function(
  target_bed, txdb_obj, org_db_pkg, tss_region,
  gene_expr_map, min_expr, conflict_strategy,
  gr_anchors, co_dominance_ratio, biotype_order, log_message
) {
  bed_target <- read_robust_general(target_bed, min_cols = 3, desc = "Target BED")
  colnames(bed_target)[c(1, 2, 3)] <- c("chr", "start", "end")
  if (nrow(bed_target) == 0) {
    warning("Target BED contains no features; skipping target annotation.")
    return(NULL)
  }
  # Validate coordinates before conversion to GRanges.
  bed_target$start <- suppressWarnings(as.numeric(bed_target$start))
  bed_target$end <- suppressWarnings(as.numeric(bed_target$end))
  if (any(is.na(bed_target$start) | is.na(bed_target$end))) {
    stop("Target BED contains non-numeric start/end coordinates.", call. = FALSE)
  }
  if (any(!is.finite(bed_target$start) | !is.finite(bed_target$end))) {
    stop("Target BED contains non-finite coordinates (Inf, -Inf, NaN).", call. = FALSE)
  }
  if (any(bed_target$start < 0 | bed_target$end < 0)) {
    stop("Target BED contains negative coordinates.", call. = FALSE)
  }
  if (any(bed_target$start >= bed_target$end)) {
    stop("Target BED contains rows with start >= end (zero-width or invalid).", call. = FALSE)
  }
  if (any(bed_target$start != floor(bed_target$start) |
    bed_target$end != floor(bed_target$end))) {
    stop("Target BED contains non-integer coordinates. Genomic coordinates must be integers.",
      call. = FALSE
    )
  }
  chr_col <- trimws(as.character(bed_target[[1]]))
  if (anyNA(chr_col) || any(!nzchar(chr_col))) {
    stop("Target BED contains empty or missing chromosome names.", call. = FALSE)
  }
  bed_target[[1]] <- chr_col

  bed_target$start <- bed_target$start + 1 # BED is 0-based; GRanges is 1-based
  gr_bed <- .with_known_upstream_noise_suppressed({
    gr_bed <- GenomicRanges::makeGRangesFromDataFrame(bed_target)
    gr_bed$input_id <- paste0("Peak_", seq_len(nrow(bed_target)))
    names(gr_bed) <- gr_bed$input_id
    gr_bed
  })
  bed_annot <- .with_known_upstream_noise_suppressed(
    ChIPseeker::annotatePeak(gr_bed,
      TxDb = txdb_obj, tssRegion = tss_region,
      annoDb = org_db_pkg, verbose = FALSE
    )
  )
  bed_info <- format_annotation_columns(as.data.frame(bed_annot))
  if ("GENENAME" %in% colnames(bed_info)) {
    bed_info <- bed_info %>% dplyr::rename(Gene_description = GENENAME)
  }

  log_message("    Refining Target annotation...")
  bed_info <- resolve_gene_conflicts(
    bed_info, txdb_obj, org_db_pkg, tss_region,
    gene_expr_map,
    min_expr = min_expr,
    conflict_strategy = conflict_strategy,
    co_dominance_ratio = co_dominance_ratio,
    biotype_order = biotype_order,
    unmeasured_policy = "keep"
  )
  gr_bed <- .harmonize_seqlevels(gr_bed, gr_anchors, "target BED")

  list(gr_bed = gr_bed, bed_info = bed_info, n_peaks = length(gr_bed))
}

#' Internal: Find peak-to-anchor overlaps with cascade filtering
#'
#' Applies three-stage filtering: gap tolerance, minimum physical overlap,
#' and minimum overlap fraction. Emits diagnostic summary.
#'
#' @return A \code{\link[S4Vectors]{Hits}} object (may be empty).
#' @keywords internal
#' @noRd
.find_peak_anchor_overlaps <- function(
  gr_bed, gr_anchors,
  anchor_gap, anchor_min_overlap, anchor_min_frac,
  n_peaks, log_message
) {
  if (anchor_gap >= 0L) {
    hits <- GenomicRanges::findOverlaps(gr_bed, gr_anchors, maxgap = anchor_gap)
  } else {
    hits <- GenomicRanges::findOverlaps(gr_bed, gr_anchors)
  }

  n_hits_raw <- length(hits)
  if (anchor_min_overlap >= 1L && n_hits_raw > 0) {
    q_gr <- gr_bed[S4Vectors::queryHits(hits)]
    s_gr <- gr_anchors[S4Vectors::subjectHits(hits)]
    overlap_w <- GenomicRanges::width(GenomicRanges::pintersect(q_gr, s_gr))
    hits <- hits[overlap_w >= anchor_min_overlap]
  }

  if (anchor_min_frac > 0 && length(hits) > 0) {
    q_gr <- gr_bed[S4Vectors::queryHits(hits)]
    s_gr <- gr_anchors[S4Vectors::subjectHits(hits)]
    overlap_w <- GenomicRanges::width(GenomicRanges::pintersect(q_gr, s_gr))
    frac <- overlap_w / GenomicRanges::width(q_gr)
    hits <- hits[frac >= anchor_min_frac]
  }

  # Diagnostic summary
  n_peaks_hit <- length(unique(S4Vectors::queryHits(hits)))
  frac_info <- if (anchor_min_frac > 0) {
    sprintf(", min_frac=%.2f", anchor_min_frac)
  } else {
    ""
  }
  bp_info <- if (anchor_min_overlap >= 1L) {
    sprintf(", min_overlap=%s bp", format(anchor_min_overlap, big.mark = ","))
  } else {
    ""
  }
  log_message(sprintf(
    "    Peak-anchor overlap: %s / %s peaks (%.1f%%) hit anchors (gap=%s bp%s%s%s)",
    format(n_peaks_hit, big.mark = ","),
    format(n_peaks, big.mark = ","),
    n_peaks_hit / max(n_peaks, 1) * 100,
    format(anchor_gap, big.mark = ","),
    bp_info,
    frac_info,
    if (n_hits_raw != length(hits)) {
      sprintf(" [%s raw hits before filtering]", format(n_hits_raw, big.mark = ","))
    } else {
      ""
    }
  ))
  if (n_peaks_hit == 0) {
    log_message("    No peaks overlapped loop anchors. Check that target BED and loop BEDPE use the same genome build and seqlevel style.")
  }

  hits
}

#' Internal: Build target-gene linkage from peak-anchor hits
#'
#' Constructs hit_df, joins topology data, computes per-peak gene summaries,
#' and delegates to \code{.build_target_gene_links}.
#'
#' @return A list(bed_info, target_connected_loops, target_gene_links, hit_df).
#' @keywords internal
#' @noRd
.build_target_hit_linkage <- function(
  hits, gr_anchors, anchor_topo_map,
  loop_annotation_final, bed_info, map_info, ego_list_target,
  g = NULL, neighbor_hop = 0L
) {
  target_connected_loops <- NULL
  target_gene_links <- NULL
  hit_df <- NULL

  if (length(hits) > 0) {
    hit_anchor_ids <- gr_anchors$anchor_id[S4Vectors::subjectHits(hits)]
    target_connected_loops <- loop_annotation_final %>%
      dplyr::filter(a1_id %in% hit_anchor_ids |
        a2_id %in% hit_anchor_ids)
    hit_df <- data.frame(
      qid = S4Vectors::queryHits(hits),
      sid = S4Vectors::subjectHits(hits)
    )
    hit_df$anchor_id <- gr_anchors$anchor_id[hit_df$sid]
    hit_df <- hit_df %>% dplyr::left_join(anchor_topo_map, by = "anchor_id")

    anchor_loop_agg <- dplyr::bind_rows(
      loop_annotation_final %>% dplyr::select(anchor_id = a1_id, loop_ID),
      loop_annotation_final %>% dplyr::select(anchor_id = a2_id, loop_ID)
    ) %>%
      dplyr::distinct() %>%
      dplyr::group_by(anchor_id) %>%
      dplyr::summarise(
        linked_loops = .annotation_extract_ids(loop_ID),
        .groups = "drop"
      )
    hit_df <- hit_df %>% dplyr::left_join(anchor_loop_agg, by = "anchor_id")

    # Only aggregate Linked_Loop_IDs here; gene summary columns
    # (All_Loop_Connected_Genes, Regulated_promoter_genes,
    # Assigned_Target_Genes) are produced by the unified finalizer
    # .finalize_current_target_annotation() below.
    summary_df <- hit_df %>%
      dplyr::group_by(qid) %>%
      dplyr::summarise(
        Linked_Loop_IDs = .annotation_extract_ids(linked_loops),
        .groups = "drop"
      ) %>%
      dplyr::mutate(join_id = paste0("Peak_", qid))
    bed_info <- dplyr::left_join(bed_info, summary_df,
      by = c("input_id" = "join_id")
    ) %>%
      dplyr::select(-any_of(c("join_id", "qid")))

    target_gene_links <- .build_target_gene_links(
      hit_df = hit_df, bed_info = bed_info,
      loop_annotation_final = loop_annotation_final,
      map_info = map_info, ego_list_target = ego_list_target,
      g = g, neighbor_hop = neighbor_hop
    )
  } else {
    bed_info$All_Loop_Connected_Genes <- NA
    bed_info$Regulated_promoter_genes <- NA
    bed_info$Assigned_Target_Genes <- NA
    bed_info$Linked_Loop_IDs <- NA
  }

  # Fallback: build empty links if nothing was produced
  if (is.null(target_gene_links)) {
    target_gene_links <- .build_target_gene_links(
      hit_df = data.frame(qid = integer(0), anchor_id = character(0)),
      bed_info = bed_info,
      loop_annotation_final = loop_annotation_final,
      map_info = map_info, ego_list_target = ego_list_target,
      g = g, neighbor_hop = neighbor_hop
    )
  }

  list(
    bed_info = bed_info,
    target_connected_loops = target_connected_loops,
    target_gene_links = target_gene_links,
    hit_df = hit_df
  )
}

#' Internal: Extract vertex IDs from igraph objects
#'
#' Handles \code{igraph.vs}, numeric vertex IDs, and character vectors,
#' regardless of the \code{return.vs.es} option.
#' @param x An \code{igraph.vs}, numeric, or character vector.
#' @param graph An igraph object used to resolve numeric IDs to names.
#' @return Character vector of vertex names.
#' @keywords internal
#' @noRd
.igraph_vertex_ids <- function(x, graph) {
  if (is.null(x) || length(x) == 0L) {
    return(character())
  }
  if (inherits(x, "igraph.vs")) {
    return(as.character(igraph::as_ids(x)))
  }
  if (is.numeric(x)) {
    if (missing(graph) || is.null(graph)) {
      stop("A graph object is required to resolve numeric igraph vertex IDs.",
        call. = FALSE
      )
    }
    vn <- igraph::vertex_attr(graph, "name")
    if (is.null(vn)) {
      stop("The graph has no vertex 'name' attribute; ",
        "numeric vertex IDs cannot be resolved to anchor IDs.",
        call. = FALSE
      )
    }
    idx <- as.integer(x)
    valid <- !is.na(idx) & idx >= 1L & idx <= length(vn)
    if (any(!valid)) {
      stop("Invalid numeric igraph vertex IDs were detected; ",
        "the stored topology may be corrupted or incompatible.",
        call. = FALSE
      )
    }
    out <- rep(NA_character_, length(idx))
    out[valid] <- as.character(vn[idx[valid]])
    return(out[!is.na(out)])
  }
  as.character(x)
}

#' Internal: Empty TSS assignment provenance schema
#' @keywords internal
#' @noRd
.empty_tss_assignment_provenance <- function() {
  data.frame(
    anchor_id = character(),
    positional_type = character(),
    final_type = character(),
    gene_before = character(),
    TxDb_gene = character(),
    final_gene = character(),
    TSS_supported = logical(),
    TSS_support_status = character(),
    Gene_Assignment_Confidence = character(),
    Gene_Assignment_Evidence = character(),
    stringsAsFactors = FALSE
  )
}

#' Internal: Test whether an object is a terminal empty annotation result
#'
#' Returns TRUE if the object carries \code{qc_status = "empty_input"}.
#' Such objects cannot be refined and must be guarded by
#' \code{\link{refine_loop_anchors_by_expression}} and
#' \code{\link{refine_loop_anchors_by_chromatin}}.
#'
#' @param x Any object (NULL / non-list inputs are safe -- returns FALSE).
#' @return Logical.
#' @keywords internal
#' @noRd
.is_empty_annotation_result <- function(x) {
  is.list(x) &&
    is.list(x$metadata) &&
    is.list(x$metadata$parameters) &&
    identical(x$metadata$parameters$qc_status, "empty_input")
}

#' Internal: Terminal empty annotation result with standard slots
#'
#' Returns a minimal terminal empty result for diagnostics and reporting.
#' Refinement is not supported on this object -- callers must guard with
#' \code{\link{.is_empty_annotation_result}}.
#'
#' @param reason Character string explaining why the result is empty.
#' @param species Character. Genome assembly identifier.
#' @param expr_matrix_file Character or NULL.
#' @param target_bed Character or NULL.
#' @return A named list with all standard annotation output slots.
#' @keywords internal
#' @noRd
.empty_annotation_result <- function(reason, species, expr_matrix_file, target_bed) {
  list(
    anchor_loci_annotation = data.frame(),
    anchor_annotation = data.frame(),
    loop_annotation = data.frame(),
    promoter_centric_stats = data.frame(),
    distal_element_stats = data.frame(),
    target_annotation = NULL,
    target_gene_links = NULL,
    plots = list(),
    plot_list = list(),
    metadata = .build_looplook_metadata(
      fun = "annotate_peaks_and_loops",
      params = list(
        species = species,
        has_expression = !is.null(expr_matrix_file),
        has_target_bed = !is.null(target_bed),
        qc_status = "empty_input",
        qc_reason = reason
      ),
      genome_build = species
    )
  )
}

#' Internal: Relocate target annotation columns after finalization
#' @keywords internal
#' @noRd
.relocate_target_annotation_columns <- function(bed_info) {
  if (is.null(bed_info)) return(bed_info)

  # Core target columns in logical reading order:
  # ChIPseeker annotation -> loop connectivity -> primary targets ->
  # expanded/candidate -> linear-fallback fill -> evidence
  ordered_target_cols <- c(
    "Linked_Loop_IDs",
    "All_Loop_Connected_Genes",
    "Regulated_promoter_genes",
    "Assigned_Target_Genes",
    "Expanded_Target_Genes",
    "Candidate_Positional_Genes",
    "All_Loop_Connected_Genes_Filled",
    "Regulated_promoter_genes_Filled",
    "Assigned_Target_Genes_Filled",
    "Regulated_promoter_Evidence",
    "Regulated_promoter_Fallback_Evidence"
  )

  present_target <- intersect(ordered_target_cols, colnames(bed_info))
  other_cols <- setdiff(colnames(bed_info), present_target)

  # Explicit select() with ordered columns ensures nothing is dropped
  # and target columns appear in the declared order
  bed_info <- bed_info %>%
    dplyr::select(dplyr::all_of(other_cols), dplyr::all_of(present_target))
  bed_info
}

#' Internal: Integrate Optional Target BED Annotation
#'
#' Orchestrates the full target-BED annotation pipeline: read & annotate,
#' compute peak-anchor overlaps, build gene linkage, and finalize evidence
#' columns.  Delegates to four internal helpers.
#'
#' @keywords internal
#' @noRd
.annotate_target_bed <- function(
  target_bed, txdb_obj, org_db_pkg, tss_region, gene_expr_map, min_expr,
  conflict_strategy,
  gr_anchors, anchor_topo_map, loop_annotation_final, map_info, ego_list_target,
  log_message,
  neighbor_hop = 0L,
  anchor_gap = -1L, anchor_min_overlap = 1L, anchor_min_frac = 0,
  co_dominance_ratio = 0.1,
  biotype_order = c("protein", "small_ncRNA", "antisense", "lncRNA", "pseudogene"),
  g = NULL,
  target_priority = c("promoter_then_distance", "distance_then_role"),
  max_primary_hop = 1L
) {
  # Step 1: Read, annotate, resolve conflicts
  prep <- .prepare_target_bed_granges(
    target_bed = target_bed, txdb_obj = txdb_obj,
    org_db_pkg = org_db_pkg, tss_region = tss_region,
    gene_expr_map = gene_expr_map, min_expr = min_expr,
    conflict_strategy = conflict_strategy,
    gr_anchors = gr_anchors,
    co_dominance_ratio = co_dominance_ratio,
    biotype_order = biotype_order,
    log_message = log_message
  )
  if (is.null(prep)) {
    return(list(
      bed_info = NULL,
      target_connected_loops = NULL,
      target_gene_links = NULL,
      hit_df = data.frame(
        qid = integer(0), sid = integer(0),
        anchor_id = character(0),
        stringsAsFactors = FALSE
      )
    ))
  }
  gr_bed <- prep$gr_bed
  bed_info <- prep$bed_info
  n_peaks <- prep$n_peaks

  # Step 2: Cascade-filtered peak-to-anchor overlaps
  hits <- .find_peak_anchor_overlaps(
    gr_bed = gr_bed, gr_anchors = gr_anchors,
    anchor_gap = anchor_gap,
    anchor_min_overlap = anchor_min_overlap,
    anchor_min_frac = anchor_min_frac,
    n_peaks = n_peaks, log_message = log_message
  )

  # Step 3: Build target-gene linkage from hits
  linkage <- .build_target_hit_linkage(
    hits = hits, gr_anchors = gr_anchors,
    anchor_topo_map = anchor_topo_map,
    loop_annotation_final = loop_annotation_final,
    bed_info = bed_info, map_info = map_info,
    ego_list_target = ego_list_target,
    g = g, neighbor_hop = neighbor_hop
  )
  bed_info <- linkage$bed_info
  target_connected_loops <- linkage$target_connected_loops
  target_gene_links <- linkage$target_gene_links
  hit_df <- linkage$hit_df

  # Unified finalizer: always call so that strict/filled/evidence columns
  # are created even when target_gene_links is empty.  Schema stability
  # is required for downstream profile_target_genes() with include_Filled=TRUE.
  target_gene_links <- .ensure_target_link_schema(target_gene_links)
  final <- .finalize_current_target_annotation(
    bed_info = bed_info,
    all_target_gene_links = target_gene_links,
    has_expression = FALSE,
    target_priority = target_priority,
    max_primary_hop = max_primary_hop
  )
  bed_info <- final$target_annotation
  target_gene_links <- final$target_gene_links

  # Step 4: Column ordering (evidence/fallback/membership already done by
  # unified finalizer at Step 3b)
  bed_info <- .relocate_target_annotation_columns(bed_info)

  # Step 5: Store peak-to-anchor mapping for downstream chromatin-aware
  # target recomputation (refine_loop_anchors_by_chromatin with
  # recompute_targets = TRUE).
  stored_hit_df <- if (!is.null(hit_df) && nrow(hit_df) > 0) {
    hit_df %>%
      dplyr::select(dplyr::any_of(c("qid", "sid", "anchor_id"))) %>%
      dplyr::distinct()
  } else {
    data.frame(
      qid = integer(0), sid = integer(0),
      anchor_id = character(0),
      stringsAsFactors = FALSE
    )
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
#'   Same-experiment HiChIP / ChIA-PET (default) \tab \code{-1} \tab \code{1} \tab \code{0} \cr
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
#' @param species Character. Genome assembly string (e.g. \code{"hg38"},
#'   \code{"mm10"}, \code{"danRer11"}, \code{"dm6"}). Default: \code{"hg38"}.
#'   When \code{txdb} and \code{org_db} are \code{NULL}, auto-resolved from
#'   built-in species (hg38/hg19/mm10/mm9); for any other species you must
#'   pass \code{txdb} and \code{org_db} as objects or package name strings
#'   directly (e.g. \code{txdb = TxDb.Dmelanogaster.UCSC.dm6.ensGene},
#'   \code{org_db = "org.Dm.eg.db"}). The \code{species} string is also used
#'   for karyotype ideograms and JASPAR motif species filtering.
#' @param tss_region Numeric vector of length 2. Promoter window around the TSS in bp. Default: \code{c(-2000, 2000)} (+/-2 kb; typical for mammalian protein-coding genes; may need widening for broad domains like HOX clusters, or narrowing for compact genomes).
#' @param anchor_merge_gap Integer. Merge loop anchors within this many bp on
#'   the same chromosome before building the connectivity graph. Default
#'   \code{0} (no merging; anchors must match exactly on chr, start, end).
#'   Set to \code{50--200} when input BEDPE comes directly from a loop caller
#'   without prior replicate consolidation -- the same biological anchor may
#'   appear with slightly different boundaries in different loop rows, which
#'   would otherwise fragment the graph. Not needed when input is
#'   pre-consolidated via \code{\link{consolidate_chromatin_loops}}.
#' @param out_dir Character. Output directory for the Excel results file. Default: \code{"./results"}.
#' @param expr_matrix_file Optional path to a normalised expression matrix (genes x samples). Accepts steady-state RNA-seq (TPM/FPKM), nascent transcription data (NET-seq, PRO-seq, GRO-seq, TT-seq), or CAGE-seq. Enables expression-aware conflict resolution. Default: \code{NULL}.
#' @param sample_columns Character vector or integer indices. Columns in \code{expr_matrix_file} to average for baseline expression. Default: \code{NULL}.
#' @param min_expr Numeric. Minimum expression value for a gene to be considered
#'   active during anchor-level conflict resolution. Used only when
#'   \code{expr_matrix_file} is provided. Default: \code{0} (any non-zero
#'   expression). Increase to \code{1} or higher to require
#'   stronger evidence. Note: when \code{min_expr = 0}, the code uses a strict
#'   greater-than comparison (\code{> 0}) to exclude truly undetected genes;
#'   when \code{min_expr >= 1}, it uses \code{>= min_expr}. For nascent
#'   transcription data (NET-seq, PRO-seq), gene-body aggregated signal
#'   should be used as the quantitative input.
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
#' @param biotype_order Character vector. Custom ordering of biotype categories
#'   for gene conflict resolution (passed to \code{\link{resolve_gene_conflicts}}).
#'   Five keywords: \code{"protein"}, \code{"small_ncRNA"}, \code{"antisense"},
#'   \code{"lncRNA"}, \code{"pseudogene"}. Listed categories get top priority;
#'   unlisted keep their default relative order appended at the bottom.
#'   Default: \code{c("protein", "small_ncRNA", "antisense", "lncRNA", "pseudogene")}.
#' @param project_name Character. Prefix for output files and plot titles. Default: \code{"HiChIP"}.
#' @param color_palette Character. RColorBrewer palette name. Default: \code{"Set2"}.
#' @param karyo_bin_size Integer. Bin width in base pairs (bp) for karyotype heatmaps. Default: \code{1e5} (100 kb). Typical range: 5e4-5e5 depending on genome size and resolution.
#' @param neighbor_hop Integer. k-hop ego-network expansion order via \code{igraph::ego()}. \code{0} (default) restricts to direct loop contacts. \code{1} additionally includes 2-hop expanded targets for exploratory network analysis. Values greater than \code{1} are not supported. Target gene assignment searches one additional hop (\code{neighbor_hop + 1}) to capture genes at the opposite anchor of directly connected loops.
#' @param anchor_gap Integer. Candidate search radius in bp between a
#'   peak and a loop anchor. \code{-1L} (default): strict physical overlap
#'   required. When \code{>= 0}, expands the candidate search window
#'   (e.g. \code{200} for cross-experiment coordinate shifts).
#'   \strong{Note:} final retention always requires at least
#'   \code{anchor_min_overlap} bp of physical overlap; this parameter
#'   only controls which candidates are evaluated, not whether
#'   proximity-only pairs are retained.
#' @param anchor_min_overlap Integer. Minimum required physical overlap
#'   (bp) between a peak and an anchor. Default \code{1L}: at least 1 bp
#'   of actual sequence overlap required -- proximity-only pairs within
#'   the \code{anchor_gap} window but without physical overlap are
#'   excluded. Increase to \code{10-100} for broad peaks to avoid
#'   spurious boundary overlaps.
#' @param anchor_min_frac Numeric (0-1). After the first two filters, what fraction of the \emph{peak} width must physically overlap the anchor? Default \code{0}: any fraction accepted. Set to \code{0.1-0.5} when peaks are broad (e.g. H3K27ac domains, 2-5 kb) so a 1 bp overlap does not link the entire broad peak. Ignored for point features (SNPs, eQTLs). Applied last, only to pairs that passed \code{anchor_gap} and \code{anchor_min_overlap}.
#' @param hub_percentile Numeric (0-1). Loop-count quantile for hub detection. Default: \code{0.95}. Genes or distal elements with connectivity at or above this quantile are flagged as hubs. A minimum floor of 3 (promoter-centric) or 2 (distal) is applied to avoid false hubs in small datasets.
#' @param hub_metric Character. Which connectivity count to use for hub detection. \code{"unique_contacts"} (default): counts distinct neighbour anchor IDs, robust to duplicate/replicate loop rows. \code{"total_loops"}: counts all loop rows (backward compatible; may inflate hub calls when input contains unconsolidated replicates).
#' @param write_output Logical. If \code{TRUE}, write the Excel workbook to \code{out_dir} (default: \code{FALSE}). If \code{FALSE}, return results without creating directories or files.
#' @param quiet Logical. If \code{TRUE}, suppress progress messages while preserving warnings. Default: \code{FALSE}.
#' @param target_priority Character. How to prioritise among multiple candidate
#'   target genes per input feature. The policy applies only within primary
#'   target links (default \code{path_length <= 1}); longer paths are reported
#'   separately as \code{Expanded_Target_Genes} and do not participate in
#'   \code{Assigned_Target_Genes} selection.
#'   \code{"promoter_then_distance"} (default): within primary links,
#'   strict promoter genes outrank strict gene-body genes, which outrank
#'   structurally supported strict enhancer candidates; within each tier
#'   shorter paths win.  Exception: direct strict promoter--promoter contacts
#'   (same loop, path0 + path1) co-assign both endpoints via union,
#'   because the technical endpoint orientation does not reflect
#'   biological regulatory direction.
#'   \code{"distance_then_role"}: within primary links, path-length dominates --
#'   the closest linked gene wins; at equal distance promoter beats gene-body
#'   (legacy behaviour).
#'   The policy affects \code{Assigned_Target_Genes} only.
#'   \code{Regulated_promoter_genes} always reports all promoter-linked
#'   genes regardless of the chosen policy.
#'
#' @return A named list:
#' \itemize{
#'   \item \code{target_annotation} -- Target features (peaks) with gene assignments.
#'     Key columns include:
#'     \itemize{
#'       \item \code{All_Loop_Connected_Genes}: Inclusive provenance union of all loop-anchor gene links. May include strict assignment-eligible targets and non-strict positional/enhancer candidates. Not a confirmed target-gene set.
#'       \item \code{Regulated_promoter_genes}: Promoter genes supported by loop-anchor context.
#'       \item \code{Assigned_Target_Genes}: Policy-based 3D assignment (default: promoter-first, then shorter path wins; see \code{target_priority}).
#'       \item \code{*_Filled} variants: Linear nearest-gene fallback when strict 3D assignments are empty.
#'       \item \code{Regulated_promoter_Evidence}: Provenance of \code{Regulated_promoter_genes}
#'         (e.g., \code{local_promoter_overlap}, \code{distal_promoter}).
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
#'       \item \code{gene_role}: \code{"promoter"}, \code{"gene_body"},
#'       \code{"enhancer_candidate"}, \code{"positional_candidate"},
#'       or \code{"linear_annotation"}.
#'       \item \code{source}: \code{"loop_anchor"} (3D-derived) or \code{"linear_annotation"} (nearest gene).
#'       \item \code{evidence}: Provenance label --
#'         \code{"local_promoter_overlap"} (peak overlaps anchor promoter),
#'         \code{"distal_promoter"} (promoter on the distal loop anchor),
#'         \code{"gene_body_context"} / \code{"distal_gene_body_context"} (gene body linkage),
#'         \code{"local_enhancer_candidate"} / \code{"distal_enhancer_candidate"} /
#'         \code{"expanded_enhancer_candidate"} (enhancer-associated linkage),
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
#'   \item \code{loop_annotation} -- Annotated 3D interactome with
#'     \code{Putative_Target_Genes} (all P/G anchor genes connected
#'     through the loop network) and \code{Promoter_Target_Genes}
#'     (promoter-only subset, P-side genes).  G--P/P--G asymmetry is
#'     respected.  Does not include linear nearest-gene fallback.
#'   \item \code{anchor_loci_annotation} -- Non-redundant anchor-locus genomic classifications after within-cluster interval reduction.
#'   \item \code{anchor_annotation} -- Backward-compatible alias of \code{anchor_loci_annotation}.
#'   \item \code{promoter_centric_stats} -- Gene-level connectivity statistics.
#'   \item \code{distal_element_stats} -- Distal anchor connectivity (E, dual, G, eG).
#'   \item \code{plots} -- Named list of ggplot/grob objects:
#'     \code{Basic_Donut}, \code{Basic_Circular}, \code{Basic_Flower},
#'     \code{Karyo_LoopGenes}, \code{Karyo_Anchors},
#'     \code{Anchor_Genomic_Distribution}, and (when \code{target_bed} is
#'     provided) \code{Karyo_TargetGenes}, \code{Target_Rose},
#'     \code{Target_Genomic_Distribution}, \code{Target_Loop_Genomic_Distribution}.
#'   \item \code{plot_list} -- Backward-compatible alias of \code{plots}.
#'   \item \code{metadata} -- Internal metadata list (parameters, versions). Not intended for direct use.
#' }
#' The returned object carries a \code{looplook_anchor_state} attribute (access via
#' \code{attr(result, "looplook_anchor_state")}) containing internal anchor
#' topology data required by \code{\link{refine_loop_anchors_by_chromatin}}
#' with \code{recompute_targets = TRUE}.
#' If \code{write_output = TRUE}, also writes a multi-sheet Excel workbook to \code{out_dir}.
#'
#' @seealso \code{\link{refine_loop_anchors_by_expression}} for expression-aware refinement,
#'   \code{\link{refine_loop_anchors_by_chromatin}} for chromatin-aware reclassification,
#'   \code{\link{profile_target_genes}} for automated functional profiling.
#'
#' @export
#'
#' @examples
#' # Minimal runnable example for package checks
#' if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
#'   txdb_example <- AnnotationDbi::loadDb(
#'     system.file("extdata", "hg19_knownGene_sample.sqlite", package = "GenomicFeatures")
#'   )
#'   bedpe_path <- tempfile(fileext = ".bedpe")
#'   writeLines(
#'     "chr6\t10412000\t10412600\tchr6\t10415000\t10415600",
#'     bedpe_path
#'   )
#'
#'   res <- annotate_peaks_and_loops(
#'     bedpe_file = bedpe_path,
#'     txdb = txdb_example,
#'     org_db = "org.Hs.eg.db",
#'     species = "hg19",
#'     out_dir = tempdir(),
#'     project_name = "Quick_Example",
#'     write_output = FALSE,
#'     quiet = TRUE
#'   )
#'   head(res$loop_annotation, 1)
#' }
annotate_peaks_and_loops <- function(
  bedpe_file,
  target_bed = NULL,
  txdb = NULL,
  org_db = NULL,
  species = "hg38",
  tss_region = c(-2000, 2000),
  anchor_merge_gap = 0,
  out_dir = "./results",
  expr_matrix_file = NULL,
  sample_columns = NULL,
  project_name = "HiChIP",
  color_palette = "Set2",
  karyo_bin_size = 1e5,
  neighbor_hop = 0,
  hub_percentile = 0.95,
  hub_metric = c("unique_contacts", "total_loops"),
  min_expr = 0,
  conflict_strategy = c("biotype_first", "expression_first"),
  co_dominance_ratio = 0.1,
  biotype_order = c("protein", "small_ncRNA", "antisense", "lncRNA", "pseudogene"),
  anchor_gap = -1L,
  anchor_min_overlap = 1L,
  anchor_min_frac = 0,
  write_output = FALSE,
  quiet = FALSE,
  target_priority = c("promoter_then_distance", "distance_then_role")
) {
  .assert_nonempty_string(species, "species")
  tss_region <- .validate_tss_region(tss_region)
  conflict_strategy <- match.arg(conflict_strategy)
  target_priority <- match.arg(target_priority)
  max_primary_hop <- 1L
  .validate_annotation_params(
    anchor_gap, anchor_min_overlap, anchor_min_frac,
    hub_percentile, neighbor_hop, karyo_bin_size
  )
  .assert_scalar_count(anchor_merge_gap, "anchor_merge_gap")
  hub_metric <- match.arg(hub_metric)
  log_message <- .make_log_message(quiet)
  # ChIPseeker::annotatePeak and GenomeInfoDb::seqlevelsStyle<- emit
  # genome-info lines to stdout via cat(). Suppress when quiet=TRUE.
  if (quiet) {
    sink(file = nullfile(), type = "output")
    on.exit(sink(type = "output"), add = TRUE)
  }

  .ensure_out_dir(write_output, out_dir)

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
  if (!is.null(expr_matrix_file)) {
    log_message("Step 0: Loading expression data...")
    gene_expr_map <- load_expression_matrix(expr_matrix_file, sample_columns)
    log_message("    >>> Expression loaded for ", length(gene_expr_map), " genes.")
  }

  # --- Step 1 & 2: Read BEDPE and cluster anchors ---
  bedpe_data <- .read_and_cluster_bedpe(
    bedpe_file, log_message, quiet,
    anchor_merge_gap
  )

  # --- Step 3: Biological Classification & Topology ---
  if (length(bedpe_data$gr_anchors) == 0) {
    reason <- if (!is.null(bedpe_data$qc_reason)) {
      bedpe_data$qc_reason
    } else {
      paste(
        "No valid anchors after parsing. Input file may be empty,",
        "malformed, or all coordinates failed validation (non-numeric,",
        "start >= end, negative, non-finite)."
      )
    }
    if (is.null(bedpe_data$qc_status) ||
      bedpe_data$qc_status != "empty_input") {
      warning("No valid loop anchors found; returning empty annotation.")
    }
    return(.empty_annotation_result(
      reason,
      species = species,
      expr_matrix_file = expr_matrix_file,
      target_bed = target_bed
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
    co_dominance_ratio = co_dominance_ratio,
    biotype_order = biotype_order
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
    hub_metric = hub_metric,
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
      neighbor_hop = neighbor_hop,
      anchor_gap = anchor_gap,
      anchor_min_overlap = anchor_min_overlap,
      anchor_min_frac = anchor_min_frac,
      co_dominance_ratio = co_dominance_ratio,
      biotype_order = biotype_order,
      g = bedpe_data$g,
      target_priority = target_priority,
      max_primary_hop = max_primary_hop
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

  # Save intermediate anchor state for downstream chromatin-aware recomputation.
  # Freeze the original type_code and gene so that repeated refine runs do
  # not degrade Refinement_Action labels, and positional gene can be recovered
  # after expression/chromatin refinement clears or overwrites current gene.
  if (!is.null(class_data$map_info)) {
    class_data$map_info$original_type_code <- class_data$map_info$type_code
    class_data$map_info$positional_gene <- class_data$map_info$SYMBOL
  }
  anchor_state <- if (!is.null(class_data$map_info)) {
    list(
      map_info = class_data$map_info,
      anchor_topo_map = class_data$anchor_topo_map,
      gr_anchors = bedpe_data$gr_anchors,
      ego_list_loop = class_data$ego_list_loop,
      ego_list_target = class_data$ego_list_target,
      neighbor_hop = neighbor_hop,
      tss_region = as.numeric(tss_region),
      g = bedpe_data$g,
      txdb_obj = txdb_obj,
      org_db_pkg = org_db_pkg,
      target_hit_df = if (exists("target_res", inherits = FALSE) && !is.null(target_res)) target_res$hit_df else NULL,
      gene_assignment_policy = list(
        gene_expr_map         = gene_expr_map,
        conflict_min_expr     = min_expr,
        conflict_strategy     = conflict_strategy,
        co_dominance_ratio    = co_dominance_ratio,
        biotype_order         = biotype_order,
        unmeasured_policy     = "keep"
      )
    )
  } else {
    NULL
  }

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
        hub_metric = hub_metric,
        unmeasured_policy = "keep",
        target_priority = match.arg(target_priority),
        max_primary_hop = 1L,
        min_expr = min_expr,
        conflict_strategy = conflict_strategy,
        co_dominance_ratio = co_dominance_ratio,
        has_target_bed = !is.null(target_bed),
        expression_requested = !is.null(expr_matrix_file),
        expression_used = !is.null(gene_expr_map),
        anchor_merge_gap = bedpe_data$merge_stats$anchor_merge_gap,
        n_raw_anchors = bedpe_data$merge_stats$n_raw_anchors,
        n_merged_anchors = bedpe_data$merge_stats$n_merged_anchors,
        n_self_loops_removed = bedpe_data$merge_stats$n_self_loops,
        n_parallel_edges_collapsed = bedpe_data$merge_stats$n_parallel_collapsed,
        n_unique_edges = bedpe_data$merge_stats$n_unique_edges,
        anchor_gap = anchor_gap,
        anchor_min_overlap = anchor_min_overlap,
        anchor_min_frac = anchor_min_frac,
        resource_inherited = FALSE,
        resource_txdb_pkg = tx_res$pkg,
        resource_orgdb_pkg = org_db_pkg
      ),
      genome_build = species,
      score_semantics = "loop type classification; higher n_members = more evidence",
      database_versions = .record_database_versions(species,
        txdb_obj = txdb_obj, org_db_pkg = org_db_pkg,
        resource_reused = FALSE
      )
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
  .assert_scalar_count(anchor_gap, "anchor_gap", min = -1L)
  .assert_scalar_count(anchor_min_overlap, "anchor_min_overlap", min = 1L)
  .assert_scalar_number(anchor_min_frac, "anchor_min_frac", min = 0, max = 1)
  if (!is.numeric(hub_percentile) || length(hub_percentile) != 1L ||
    is.na(hub_percentile) || !is.finite(hub_percentile) ||
    hub_percentile <= 0 || hub_percentile > 1) {
    stop("hub_percentile must be a finite number in (0, 1]", call. = FALSE)
  }
  .assert_scalar_count(neighbor_hop, "neighbor_hop")
  if (neighbor_hop > 1L) {
    stop("neighbor_hop > 1 is not supported. Direct 1-hop loop targets are always evaluated; larger values trigger computationally intractable network expansion.", call. = FALSE)
  }
  .assert_scalar_count(karyo_bin_size, "karyo_bin_size", min = 1L)
}

#' Internal: Generate annotation plots with optional message silencing
#' @keywords internal
#' @noRd
.generate_annotation_plots <- function(
  loop_annotation_final, bed_info,
  cluster_info, target_connected_loops, txdb_obj, org_db_pkg, species,
  project_name, color_palette, karyo_bin_size, quiet
) {
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
.read_and_cluster_bedpe <- function(bedpe_file, log_message, quiet,
                                    anchor_merge_gap = 0) {
  .assert_scalar_count(anchor_merge_gap, "anchor_merge_gap")
  anchor_merge_gap <- as.integer(anchor_merge_gap)
  log_message("Step 1: Reading BEDPE file...")
  loops <- data.table::fread(bedpe_file, header = FALSE)
  ncol_loops <- ncol(loops)
  if (ncol_loops < 6L) {
    if (nrow(loops) == 0L && ncol_loops == 0L) {
      warning("BEDPE file is empty (0 rows, 0 columns).")
      return(list(
        loops = data.frame(),
        anchors = data.frame(),
        gr_anchors = GenomicRanges::GRanges(),
        cluster_regions = GenomicRanges::GRanges(),
        g = igraph::make_empty_graph(n = 0L, directed = FALSE),
        merge_stats = list(
          n_raw_anchors = 0L, n_merged_anchors = 0L,
          n_self_loops = 0L, n_parallel_collapsed = 0L,
          n_unique_edges = 0L,
          anchor_merge_gap = as.integer(anchor_merge_gap),
          merge_diagnostics = NULL,
          n_suspicious_merges = 0L,
          max_span_ratio = NA_real_
        ),
        qc_status = "empty_input",
        qc_reason = "Empty BEDPE file (0 rows, 0 columns)"
      ))
    }
    stop("BEDPE file must have at least 6 columns; found ", ncol_loops, ".",
      call. = FALSE
    )
  }
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
    unique()
  # Merge nearby anchors within anchor_merge_gap to avoid fragmenting the
  # loop graph when the same biological anchor appears with slightly different
  # coordinates across loop rows.  Default 0 = no merging (backward compatible).
  # Build a raw-coords ->merged-anchor-ID mapping via reduce + findOverlaps,
  # then use the mapping to assign anchor IDs to the original loop coordinates.
  if (anchor_merge_gap > 0) {
    gr_raw <- GenomicRanges::makeGRangesFromDataFrame(anchors)
    raw_w <- GenomicRanges::width(gr_raw)
    # reduce(min.gapwidth = k) merges ranges with gap < k.
    # We want "within X bp" ->gap <= X ->use min.gapwidth = X + 1.
    gr_merged <- GenomicRanges::reduce(gr_raw,
      min.gapwidth = anchor_merge_gap + 1L
    )
    ov <- GenomicRanges::findOverlaps(gr_raw, gr_merged)
    # each raw anchor maps to exactly one merged region (reduce is a partition)
    merged_hit <- S4Vectors::subjectHits(ov)
    merged_w <- GenomicRanges::width(gr_merged)

    # --- Merge diagnostics: detect transitive chaining ---
    # reduce() is transitive: if A is near B and B near C, A+B+C merge
    # even when A and C are far apart.  Compute per-merged-region
    # diagnostics to flag suspiciously wide merged anchors.
    n_raw_per_merged <- tabulate(merged_hit, nbins = length(gr_merged))
    raw_median_w <- tapply(raw_w[S4Vectors::queryHits(ov)], merged_hit, stats::median)
    span_ratio <- merged_w / raw_median_w
    # Max pairwise gap: true interval gap between consecutive raw
    # anchors within each merged region (next_start - prev_end - 1).
    # This is the most direct measure of chaining.
    max_gap <- vapply(seq_along(gr_merged), function(i) {
      idx <- which(merged_hit == i)
      if (length(idx) < 2L) {
        return(0L)
      }
      ord <- order(
        GenomicRanges::start(gr_raw)[idx],
        GenomicRanges::end(gr_raw)[idx]
      )
      s <- GenomicRanges::start(gr_raw)[idx[ord]]
      e <- GenomicRanges::end(gr_raw)[idx[ord]]
      max(pmax(0L, s[-1L] - e[-length(e)] - 1L))
    }, integer(1L))

    merge_diag <- data.frame(
      merged_anchor_id = paste0("A", seq_along(gr_merged)),
      merged_width = merged_w,
      n_raw_anchors = n_raw_per_merged,
      raw_median_width = raw_median_w,
      span_ratio = round(span_ratio, 2),
      max_pairwise_gap = max_gap,
      stringsAsFactors = FALSE
    )
    # Flag suspicious merges: large span ratio with many raw anchors
    # is the hallmark of transitive chaining.
    suspicious <- merge_diag$span_ratio > 10 & merge_diag$n_raw_anchors > 3
    n_suspicious <- sum(suspicious)

    # Build lookup: raw coords ->canonical merged anchor coords
    raw_to_merged <- data.frame(
      chr = anchors$chr,
      start = anchors$start,
      end = anchors$end,
      merged_chr = as.character(GenomicRanges::seqnames(gr_merged))[merged_hit],
      merged_start = GenomicRanges::start(gr_merged)[merged_hit],
      merged_end = GenomicRanges::end(gr_merged)[merged_hit],
      stringsAsFactors = FALSE
    )
    # Build merged anchor table with new IDs (one row per merged region)
    merged_anchors <- unique(raw_to_merged[, c("merged_chr", "merged_start", "merged_end")])
    colnames(merged_anchors) <- c("chr", "start", "end")
    merged_anchors$anchor_id <- paste0("A", seq_len(nrow(merged_anchors)))
    # Create two lookup tables keyed by raw coords, carrying anchor ID
    a1_lookup <- merge(raw_to_merged, merged_anchors,
      by.x = c("merged_chr", "merged_start", "merged_end"),
      by.y = c("chr", "start", "end")
    )[, c("chr", "start", "end", "anchor_id")]
    colnames(a1_lookup) <- c("chr1", "start1", "end1", "a1_id")
    a2_lookup <- a1_lookup
    colnames(a2_lookup) <- c("chr2", "start2", "end2", "a2_id")
    loops <- dplyr::left_join(loops, a1_lookup,
      by = c("chr1", "start1", "end1")
    )
    loops <- dplyr::left_join(loops, a2_lookup,
      by = c("chr2", "start2", "end2")
    )
    # Detect and remove self-loops: when both loop anchors merge into
    # the same merged anchor (a1_id == a2_id), the loop is degenerate
    # and would inflate degree / connectivity statistics.
    self_loop_mask <- !is.na(loops$a1_id) & !is.na(loops$a2_id) &
      as.character(loops$a1_id) == as.character(loops$a2_id)
    n_self <- sum(self_loop_mask, na.rm = TRUE)
    if (n_self > 0) {
      loops <- loops[!self_loop_mask, ]
      if (!quiet) {
        message(sprintf(
          "    Removed %d self-loop(s) caused by anchor merging", n_self
        ))
      }
    }
    # Set anchors to merged_anchors for downstream graph/computation
    anchors <- merged_anchors
    if (!quiet) {
      message(sprintf(
        "    Anchor merging: %d raw -> %d merged (gap = %d bp)",
        nrow(raw_to_merged), nrow(merged_anchors), anchor_merge_gap
      ))
      if (n_suspicious > 0) {
        warning(sprintf(
          "%d merged anchor(s) have span_ratio > 10 with > 3 raw anchors -- possible transitive chaining. Consider reducing `anchor_merge_gap`.",
          n_suspicious
        ), call. = FALSE)
        worst_rows <- merge_diag[suspicious, , drop = FALSE]
        worst <- worst_rows[order(
          -worst_rows$span_ratio,
          -worst_rows$n_raw_anchors
        )[seq_len(min(3, n_suspicious))], ]
        message(sprintf(
          "      Worst: %s  (width=%d bp, n_raw=%d, span_ratio=%.1f, max_gap=%d bp)",
          worst$merged_anchor_id, worst$merged_width,
          worst$n_raw_anchors, worst$span_ratio, worst$max_pairwise_gap
        ))
      }
    }
  } else {
    anchors <- anchors %>%
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
    # Self-loops can also exist in raw BEDPE without anchor merging
    # (both ends map to identical coordinates).  Remove them consistently
    # so that degree / connectivity statistics are not inflated.
    self_loop_mask <- !is.na(loops$a1_id) & !is.na(loops$a2_id) &
      as.character(loops$a1_id) == as.character(loops$a2_id)
    n_self_raw <- sum(self_loop_mask, na.rm = TRUE)
    if (n_self_raw > 0) {
      loops <- loops[!self_loop_mask, ]
      if (!quiet) {
        message(sprintf(
          "    Removed %d raw self-loop(s) from input", n_self_raw
        ))
      }
    }
  }

  log_message("Step 2: Clustering loops...")
  valid_loops <- loops %>% dplyr::filter(!is.na(a1_id) & !is.na(a2_id))
  # Aggregate parallel edges after anchor merge: when multiple raw loops
  # map to the same (a1_id, a2_id) pair, they are treated as independent
  # support for a single biological contact.  The graph topology (degree,
  # clustering, shortest-path hops) is based on unique edges; n_support
  # is stored as an edge attribute for downstream confidence weighting.
  #
  # Semantics of n_support: n_support counts raw BEDPE row multiplicity
  # per canonical edge -- i.e. how many loop rows from the input map to
  # the same (a1_id, a2_id) pair.  It does NOT automatically incorporate
  # replicate-level metadata (n_reps, n_members) from consolidated input.
  # When input has been pre-consolidated via consolidate_chromatin_loops(),
  # every consensus loop contributes one row regardless of its n_reps,
  # so all edges carry n_support = 1.  Downstream equal-hop tie-breaking
  # then defaults to a deterministic vertex-name sort order.
  n_raw_edges <- nrow(valid_loops)
  edge_df <- valid_loops %>%
    dplyr::group_by(a1_id, a2_id) %>%
    dplyr::summarise(n_support = dplyr::n(), .groups = "drop")
  n_unique_edges <- nrow(edge_df)
  g <- igraph::graph_from_data_frame(
    edge_df[, c("a1_id", "a2_id")],
    directed = FALSE
  )
  igraph::E(g)$n_support <- edge_df$n_support
  if (n_unique_edges < n_raw_edges) {
    log_message(sprintf(
      "    Parallel edges: %d raw loops collapsed to %d unique edges (%d redundant)",
      n_raw_edges, n_unique_edges, n_raw_edges - n_unique_edges
    ))
  }
  comp <- igraph::components(g)
  anchors$cluster_id <- NA
  comm <- intersect(anchors$anchor_id, names(comp$membership))
  if (length(comm) > 0) {
    anchors$cluster_id[match(comm, anchors$anchor_id)] <-
      as.character(comp$membership[comm])
  }
  anchors <- anchors %>% dplyr::filter(!is.na(cluster_id))
  loops <- loops %>%
    dplyr::left_join(
      anchors %>% dplyr::select(anchor_id, cluster_id),
      by = c("a1_id" = "anchor_id")
    )
  gr_anchors <- .with_known_upstream_noise_suppressed({
    gr_anchors <- GenomicRanges::makeGRangesFromDataFrame(
      anchors,
      keep.extra.columns = TRUE
    )
    gr_anchors$anchor_id <- anchors$anchor_id
    names(gr_anchors) <- anchors$anchor_id
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
    gr_anchors = gr_anchors, cluster_regions = cluster_regions, g = g,
    merge_stats = list(
      n_raw_anchors = if (anchor_merge_gap > 0) nrow(raw_to_merged) else nrow(anchors),
      n_merged_anchors = nrow(anchors),
      n_self_loops = if (anchor_merge_gap > 0) n_self else n_self_raw,
      n_parallel_collapsed = n_raw_edges - n_unique_edges,
      n_unique_edges = n_unique_edges,
      anchor_merge_gap = anchor_merge_gap,
      merge_diagnostics = if (anchor_merge_gap > 0) merge_diag else NULL,
      n_suspicious_merges = if (anchor_merge_gap > 0) n_suspicious else 0L,
      max_span_ratio = if (anchor_merge_gap > 0 && nrow(merge_diag) > 0) {
        max(merge_diag$span_ratio)
      } else {
        1
      }
    )
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
  co_dominance_ratio = 0.1,
  biotype_order = c("protein", "small_ncRNA", "antisense", "lncRNA", "pseudogene")
) {
  log_message("Step 3: Biological Classification & Topology...")
  anchor_anno <- .with_known_upstream_noise_suppressed(
    ChIPseeker::annotatePeak(
      gr_anchors,
      TxDb = txdb_obj, tssRegion = tss_region,
      annoDb = org_db_pkg, verbose = FALSE
    )
  )
  anchor_anno_df <- format_annotation_columns(as.data.frame(anchor_anno))
  anchor_anno_df <- resolve_gene_conflicts(
    anchor_anno_df, txdb_obj, org_db_pkg, tss_region,
    gene_expr_map,
    min_expr = min_expr,
    conflict_strategy = conflict_strategy,
    co_dominance_ratio = co_dominance_ratio,
    biotype_order = biotype_order,
    unmeasured_policy = "keep"
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
    Map(
      .loop_locus_genes, loops_annotated$t1, loops_annotated$t2,
      loops_annotated$s1, loops_annotated$s2
    ),
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
    g,
    order = input_hop, nodes = nodes_in_graph, mode = "all"
  )
  names(ego_list_loop) <- nodes_in_graph
  ego_list_target <- list()
  # Target ego networks are generated lazily in .build_target_gene_links()
  # after target-anchor overlaps are known, avoiding precomputation of
  # 2-hop neighbourhoods for every graph node.

  anchor_topo_map <- data.frame(
    anchor_id = nodes_in_graph,
    topo_genes_p = vapply(ego_list_loop, function(x) {
      .ids_to_genes_simple(.igraph_vertex_ids(x, g), lookup_p_symbol)
    }, character(1)),
    topo_genes_pg = vapply(ego_list_loop, function(x) {
      .ids_to_genes_simple(.igraph_vertex_ids(x, g), lookup_pg_symbol)
    }, character(1)),
    stringsAsFactors = FALSE
  )
  anchor_topo_map[is.na(anchor_topo_map)] <- NA_character_

  list(
    map_info = map_info,
    loops_annotated = loops_annotated,
    anchor_topo_map = anchor_topo_map,
    ego_list_loop = ego_list_loop,
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
  hub_percentile, hub_metric, log_message
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
    dplyr::left_join(
      anchor_topo_map %>% dplyr::select(anchor_id, p1 = topo_genes_p),
      by = c("a1_id" = "anchor_id")
    ) %>%
    dplyr::left_join(
      anchor_topo_map %>% dplyr::select(anchor_id, p2 = topo_genes_p),
      by = c("a2_id" = "anchor_id")
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      proximate_loop_gene = dplyr::case_when(
        (!is.na(t1) & t1 == "G" & !is.na(t2) & t2 == "P") ~
          extract_genes(pg2),
        (!is.na(t1) & t1 == "P" & !is.na(t2) & t2 == "G") ~
          extract_genes(pg1),
        TRUE ~ extract_genes(c(pg1, pg2))
      ),
      proximate_promoter_gene = dplyr::case_when(
        (!is.na(t1) & t1 == "G" & !is.na(t2) & t2 == "P") ~
          extract_genes(p2),
        (!is.na(t1) & t1 == "P" & !is.na(t2) & t2 == "G") ~
          extract_genes(p1),
        TRUE ~ extract_genes(c(p1, p2))
      )
    ) %>%
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
      all_cluster_loop_genes, single_loop_genes,
      proximate_loop_gene, proximate_promoter_gene,
      a1_id, a2_id
    ) %>%
    dplyr::rename(
      All_Anchor_Genes = single_loop_genes,
      Putative_Target_Genes = proximate_loop_gene,
      Promoter_Target_Genes = proximate_promoter_gene,
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
        cluster_regions,
        TxDb = txdb_obj,
        tssRegion = tss_region, annoDb = org_db_pkg, verbose = FALSE
      )
    )
    cluster_info <- format_annotation_columns(as.data.frame(gene_annot))
  }
  if ("GENENAME" %in% colnames(cluster_info)) {
    cluster_info <- cluster_info %>% dplyr::rename(Gene_description = GENENAME)
  }
  cluster_info$cluster_id <- as.character(cluster_info$cluster_id)
  cluster_info <- cluster_info %>%
    dplyr::left_join(agg_cluster_locus, by = "cluster_id")

  # --- Promoter-centric stats ---
  # Total_Loops counts every raw loop row (including parallel edges after
  # anchor merge).  n_Unique_Contacts counts distinct neighbour anchor IDs,
  # matching the graph degree semantics where parallel edges are collapsed.
  log_message("    Generating Promoter Centric Stats...")
  raw_stats_df <- dplyr::bind_rows(
    loop_annotation_final %>%
      dplyr::filter(.is_promoter_like(anchor1_type) & !is.na(anchor1_gene)) %>%
      dplyr::select(
        Gene = anchor1_gene, Neighbor_Type = anchor2_type,
        Loop_Type = loop_type, contact_id = a2_id
      ) %>%
      dplyr::mutate(Gene = as.character(Gene)),
    loop_annotation_final %>%
      dplyr::filter(.is_promoter_like(anchor2_type) & !is.na(anchor2_gene)) %>%
      dplyr::select(
        Gene = anchor2_gene, Neighbor_Type = anchor1_type,
        Loop_Type = loop_type, contact_id = a1_id
      ) %>%
      dplyr::mutate(Gene = as.character(Gene))
  ) %>%
    tidyr::separate_rows(Gene, sep = ";") %>%
    dplyr::mutate(Gene = trimws(Gene)) %>%
    dplyr::filter(Gene != "" & !is.na(Gene)) %>%
    dplyr::group_by(Gene) %>%
    dplyr::summarise(
      Total_Loops = dplyr::n(),
      n_Unique_Contacts = dplyr::n_distinct(contact_id),
      n_Linked_Promoters = sum(.is_promoter_like(Neighbor_Type), na.rm = TRUE),
      n_Linked_Distal = sum((.is_distal_like(Neighbor_Type) | .is_gene_body_like(Neighbor_Type)), na.rm = TRUE),
      n_Linked_EnhancerLike = sum(.is_distal_like(Neighbor_Type), na.rm = TRUE),
      Dominant_Interaction = names(which.max(table(Loop_Type))),
      .groups = "drop"
    )

  if (nrow(raw_stats_df) == 0) {
    promoter_centric_df <- data.frame(
      Gene = character(), Total_Loops = integer(),
      n_Unique_Contacts = integer(),
      n_Linked_Promoters = integer(), n_Linked_Distal = integer(),
      n_Linked_EnhancerLike = integer(),
      Dominant_Interaction = character(),
      Is_High_Connectivity_Gene = character(),
      Is_High_Distal_Connectivity_Gene = character(),
      stringsAsFactors = FALSE
    )
  } else {
    hub_col <- if (hub_metric == "unique_contacts") "n_Unique_Contacts" else "Total_Loops"
    final_cutoff <- max(
      quantile(raw_stats_df[[hub_col]], hub_percentile, na.rm = TRUE), 3
    )
    distal_cutoff <- max(
      quantile(raw_stats_df$n_Linked_Distal, hub_percentile, na.rm = TRUE), 2
    )
    promoter_centric_df <- raw_stats_df %>%
      dplyr::mutate(
        Is_High_Connectivity_Gene = dplyr::if_else(
          .data[[hub_col]] >= final_cutoff, "Yes", "No"
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
    loop_annotation_final, anchors, hub_percentile, hub_metric
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
  .add_sheet(wb, "Loop Annotation", loop_annotation_clean)
  .add_sheet(wb, "Anchor Loci Annotation", cluster_info)
  .add_sheet(wb, "Promoter Stats", promoter_centric_df)
  .add_sheet(wb, "Distal Element Stats", distal_element_df)
  .add_sheet(wb, "Target Annotation", bed_info)
  .add_sheet(wb, "Target Gene Links", target_gene_links, require_rows = TRUE)
  .save_workbook(wb, out_dir, project_name, "_Basic_Results.xlsx",
    "Failed to save Excel workbook: "
  )
}

#' Internal: Build distal element connectivity data frame
#' @keywords internal
#' @noRd
.build_distal_element_stats <- function(loop_annotation_final, anchors,
                                        hub_percentile, hub_metric) {
  # Distal anchor types: enhancer-like (E), dual-signature (dual),
  # and gene-body anchors (G, eG) after expression/chromatin refinement.
  distal_types <- c("E", "dual", "G", "eG")
  distal_raw_df <- dplyr::bind_rows(
    loop_annotation_final %>%
      dplyr::filter(anchor1_type %in% distal_types) %>%
      dplyr::select(
        Distal_Anchor_ID = a1_id, Distal_Type = anchor1_type,
        Neighbor_Gene = anchor2_gene, Neighbor_Type = anchor2_type,
        Loop_Type = loop_type, contact_id = a2_id
      ),
    loop_annotation_final %>%
      dplyr::filter(anchor2_type %in% distal_types) %>%
      dplyr::select(
        Distal_Anchor_ID = a2_id, Distal_Type = anchor2_type,
        Neighbor_Gene = anchor1_gene, Neighbor_Type = anchor1_type,
        Loop_Type = loop_type, contact_id = a1_id
      )
  ) %>%
    dplyr::group_by(Distal_Anchor_ID, Distal_Type) %>%
    dplyr::summarise(
      Total_Loops = dplyr::n(),
      n_Unique_Contacts = dplyr::n_distinct(contact_id),
      n_Linked_Promoters = sum(.is_promoter_like(Neighbor_Type), na.rm = TRUE),
      n_Linked_Distal = sum((.is_distal_like(Neighbor_Type) | .is_gene_body_like(Neighbor_Type)), na.rm = TRUE),
      n_Linked_EnhancerLike = sum(.is_distal_like(Neighbor_Type), na.rm = TRUE),
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
      cluster_id = character(), Distal_Type = character(),
      Total_Loops = integer(), n_Unique_Contacts = integer(),
      n_Linked_Promoters = integer(), n_Linked_Distal = integer(),
      n_Linked_EnhancerLike = integer(),
      Dominant_Interaction = character(),
      Is_High_Connectivity_Distal_Element = character(),
      Target_Genes = character(), stringsAsFactors = FALSE
    ))
  }
  hub_col <- if (hub_metric == "unique_contacts") "n_Unique_Contacts" else "Total_Loops"
  final_cutoff <- max(quantile(distal_raw_df[[hub_col]], hub_percentile,
    na.rm = TRUE
  ), 3)
  distal_raw_df %>%
    dplyr::left_join(anchor_coords_map,
      by = c("Distal_Anchor_ID" = "anchor_id")
    ) %>%
    dplyr::mutate(
      Is_High_Connectivity_Distal_Element =
        dplyr::if_else(.data[[hub_col]] >= final_cutoff, "Yes", "No")
    ) %>%
    dplyr::select(
      chr, start, end, cluster_id, Distal_Type,
      Total_Loops, n_Unique_Contacts,
      n_Linked_Promoters, n_Linked_Distal, n_Linked_EnhancerLike,
      Dominant_Interaction,
      Is_High_Connectivity_Distal_Element, Target_Genes
    ) %>%
    dplyr::arrange(dplyr::desc(Total_Loops))
}

# Unified target-gene-link schema is defined by .empty_target_gene_links().
# All functions use .ensure_target_link_schema() which references the
# empty prototype as the canonical column set.

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
    path_length = integer(),
    Mean_Expression = numeric(),
    Passes_Expression_Filter = logical(),
    Expression_State = character(),
    anchor_type_before_chromatin = character(),
    anchor_type_after_chromatin = character(),
    chromatin_confidence = character(),
    chromatin_evidence = character(),
    chromatin_target_action = character(),
    used_as_fallback = logical(),
    in_regulated_promoter = logical(),
    in_assigned_target = logical(),
    in_all_loop_connected = logical(),
    in_regulated_promoter_filled = logical(),
    in_assigned_target_filled = logical(),
    retained_after_refinement = logical(),
    strict_assignment_eligible = logical(),
    in_expanded_target = logical(),
    is_expanded_path = logical(),
    TSS_supported = logical(),
    TSS_support_status = character(),
    Gene_Assignment_Confidence = character(),
    Gene_Assignment_Evidence = character(),
    stringsAsFactors = FALSE
  )
}

#' Ensure a target_gene_links data frame has the canonical schema.
#' Missing columns are filled with NA of the correct type so that
#' results from different pipeline stages (raw, expression-refined,
#' chromatin-refined) have identical column structure.
#' @keywords internal
#' @noRd
.ensure_target_link_schema <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  proto <- .empty_target_gene_links()
  if (nrow(x) == 0L) {
    return(proto)
  }
  missing <- setdiff(names(proto), names(x))
  if (length(missing) == 0L) {
    return(x[, names(proto), drop = FALSE])
  }
  n <- nrow(x)
  for (m in missing) {
    col <- proto[[m]]
    if (is.character(col)) {
      x[[m]] <- rep(NA_character_, n)
    } else if (is.logical(col)) {
      x[[m]] <- rep(NA, n)
    } else if (is.integer(col)) {
      x[[m]] <- rep(NA_integer_, n)
    } else if (is.numeric(col)) {
      x[[m]] <- rep(NA_real_, n)
    } else {
      stop("Unsupported schema column type for: ", m)
    }
  }
  x[, names(proto), drop = FALSE]
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

#' Resolve chromatin-aware gene role and strict eligibility from origin + final type.
#'
#' Single source of truth for both anchor-level (map_info, loop_df) and
#' link-level (target_gene_links) role/eligibility decisions.  Prevents the
#' two-layer inconsistency where anchor-level used final-type-only logic and
#' link-level used origin-aware logic, producing contradictory roles for the
#' same anchor.
#'
#' Rules:
#'   - Chromatin state (final_type) determines the regulatory role category.
#'   - Gene-identity confidence depends on origin (old_type) + TSS validation:
#'     * P/eP/G/eG anchors carry structurally validated genes (TSS/gene-body overlap).
#'     * E anchors carry only ChIPseeker nearest genes -- NOT reliable for strict
#'       targets unless confirmed by independent TSS reannotation.
#'
#' @param old_type Character vector of pre-chromatin anchor types (P, eP, G, eG, E).
#' @param final_type Character vector of chromatin-refined anchor types.
#' @param tss_supported Logical vector: TRUE when TSS reannotation found a
#'   TxDb-supported promoter gene for this anchor.
#' @param has_gene Logical vector: TRUE when the anchor carries a non-empty gene symbol.
#' @return A list with two parallel vectors: \code{role} (character) and
#'   \code{strict} (logical).
#' @keywords internal
#' @noRd
.resolve_chromatin_gene_role <- function(old_type, final_type,
                                          tss_supported = FALSE,
                                          has_gene = TRUE) {
  # Coerce scalar or empty defaults to match input length (dplyr::case_when
  # does not auto-recycle scalars or zero-length vectors in compound & conditions).
  n <- length(old_type)
  if (length(tss_supported) == 0L) {
    tss_supported <- rep(FALSE, n)
  } else if (length(tss_supported) == 1L && n > 1L) {
    tss_supported <- rep(tss_supported, n)
  }
  if (length(has_gene) == 0L) {
    has_gene <- rep(TRUE, n)
  } else if (length(has_gene) == 1L && n > 1L) {
    has_gene <- rep(has_gene, n)
  }
  tss_pos <- !is.na(tss_supported) & tss_supported
  gene_ok <- !is.na(has_gene) & has_gene
  has_info <- !is.na(old_type) & !is.na(final_type)

  # Empty-input guard: dplyr::case_when with zero-length vectors can
  # produce recycling errors when conditions differ in length.
  if (n == 0L) {
    return(list(role = character(0), strict = logical(0)))
  }

  # --- gene_role ---
  role <- dplyr::case_when(
    !has_info ~ NA_character_,
    # E -> P/dual without TSS: chromatin state is promoter-like, but
    # gene identity is unresolved (nearest gene only).
    old_type == "E" &
      final_type %in% c("P", "dual") &
      !tss_pos ~ "positional_candidate",
    # All other P/dual: promoter
    final_type %in% c("P", "eP", "dual") ~ "promoter",
    # G/eG: gene_body
    final_type %in% c("G", "eG") ~ "gene_body",
    # E (all origins): enhancer_candidate
    final_type == "E" ~ "enhancer_candidate",
    TRUE ~ "positional_candidate"
  )

  # --- strict_assignment_eligible ---
  strict <- dplyr::case_when(
    !has_info ~ NA,
    # P/eP -> P/dual: annotated promoter gene, always strict
    old_type %in% c("P", "eP") &
      final_type %in% c("P", "eP", "dual") &
      gene_ok ~ TRUE,
    # G/eG -> P/dual: gene-body host gene, strict (moderate confidence)
    old_type %in% c("G", "eG") &
      final_type %in% c("P", "dual") &
      gene_ok ~ TRUE,
    # E -> P/dual WITH TSS: TSS-validated gene, strict
    old_type == "E" &
      final_type %in% c("P", "dual") &
      tss_pos &
      gene_ok ~ TRUE,
    # G/eG -> G/eG (unchanged): gene_body, always strict
    old_type %in% c("G", "eG") &
      final_type %in% c("G", "eG") &
      gene_ok ~ TRUE,
    # P/eP/G/eG -> E: structurally supported gene, strict at rank 3
    final_type == "E" &
      old_type %in% c("P", "eP", "G", "eG") &
      gene_ok ~ TRUE,
    # Everything else (E->E, E->P/dual without TSS, positional, no gene):
    # NOT strict
    TRUE ~ FALSE
  )

  list(role = role, strict = strict)
}

.target_anchor_gene_map <- function(map_info) {
  if (.is_null_or_empty(map_info) ||
    !all(c("anchor_id", "SYMBOL") %in% colnames(map_info))) {
    return(data.frame(
      anchor_id = character(), gene = character(),
      gene_role = character(),
      strict_assignment_eligible = logical(),
      stringsAsFactors = FALSE
    ))
  }

  if ("effective_gene_role" %in% colnames(map_info)) {
    # Use effective_gene_role when available (post-chromatin).
    strict_col <- if ("strict_assignment_eligible" %in% colnames(map_info)) {
      as.logical(map_info$strict_assignment_eligible)
    } else {
      map_info$effective_gene_role %in% c("promoter", "gene_body")
    }
    strict_col[is.na(strict_col)] <- map_info$effective_gene_role[is.na(strict_col)] %in%
      c("promoter", "gene_body")
    strict_col[map_info$effective_gene_role == "positional_candidate"] <- FALSE
    map_info %>%
      dplyr::mutate(
        SYMBOL = as.character(SYMBOL),
        gene_role = as.character(effective_gene_role),
        strict_assignment_eligible = strict_col
      ) %>%
      tidyr::separate_rows(SYMBOL, sep = ";") %>%
      dplyr::mutate(gene = trimws(SYMBOL)) %>%
      dplyr::filter(
        gene_role %in% c("promoter", "gene_body", "enhancer_candidate", "positional_candidate"),
        !is.na(gene), gene != ""
      ) %>%
      dplyr::select(anchor_id, gene, gene_role, strict_assignment_eligible) %>%
      dplyr::distinct()
  } else {
    map_info %>%
      dplyr::mutate(SYMBOL = as.character(SYMBOL)) %>%
      tidyr::separate_rows(SYMBOL, sep = ";") %>%
      dplyr::mutate(gene = trimws(SYMBOL)) %>%
      dplyr::filter(
        !is.na(gene), gene != "",
        .is_target_promoter_like(type_code) | .is_target_gene_body_like(type_code)
      ) %>%
      dplyr::transmute(
        anchor_id = as.character(anchor_id),
        gene = gene,
        gene_role = dplyr::if_else(.is_target_promoter_like(type_code), "promoter", "gene_body"),
        strict_assignment_eligible = TRUE
      ) %>%
      dplyr::distinct()
  }
}

.linear_target_gene_links <- function(bed_info) {
  empty <- .empty_target_gene_links()[, c(
    "input_id", "loop_ID", "anchor_id", "gene", "gene_role",
    "source", "evidence", "anchor_role"
  )]
  linear_col <- .target_linear_gene_column(bed_info)
  if (is.null(linear_col) || .is_null_or_empty(bed_info)) {
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
      path_length = NA_integer_,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, Filter(Negate(is.null), rows))
  if (is.null(out)) {
    return(empty)
  }
  out
}

#' Internal: Maximum-support shortest path via DAG dynamic programming
#'
#' Builds a shortest-path DAG (only edges on minimum-hop routes from source
#' to target) and uses layer-by-layer DP to select the path with maximum
#' total n_support.  O(V + E) -- avoids the combinatorial enumeration of
#' all shortest paths that \code{all_shortest_paths()} performs.
#'
#' @param g An igraph object with edge attribute \code{n_support}
#'   (raw BEDPE row multiplicity per canonical edge; see
#'   \code{.read_and_cluster_bedpe} for semantics).
#' @param from Source vertex.
#' @param to Target vertex.
#' @param d Unweighted shortest distance (pre-computed).
#' @return A vertex sequence of the selected path, or NULL.
#' @keywords internal
#' @noRd
.dag_max_support_path <- function(g, from, to, d) {
  # Defensive: n_support must be finite non-negative (normal graph has
  # positive integer row counts).  A corrupt graph could produce NA/NaN DP.
  support <- igraph::E(g)$n_support
  if (is.null(support) || anyNA(support) ||
    any(!is.finite(support)) || any(support < 0)) {
    stop("Graph edge attribute `n_support` must contain finite non-negative values.",
      call. = FALSE
    )
  }
  vertex_names <- igraph::V(g)$name
  from_idx <- match(as.character(from), vertex_names)
  to_idx <- match(as.character(to), vertex_names)
  if (is.na(from_idx) || is.na(to_idx)) {
    return(NULL)
  }

  # Keep distances as numeric (Inf for disconnected components).
  dist_from <- drop(igraph::distances(g,
    v = from, to = igraph::V(g),
    weights = NA, mode = "all"
  ))
  dist_to <- drop(igraph::distances(g,
    v = to, to = igraph::V(g),
    weights = NA, mode = "all"
  ))

  if (!is.finite(dist_from[to_idx]) || dist_from[to_idx] != d) {
    return(NULL)
  }

  # Build DAG: only edges on minimum-hop paths, oriented from lower
  # dist_from layer to higher dist_from layer.
  el <- igraph::as_edgelist(g, names = FALSE)
  u0 <- el[, 1L]
  v0 <- el[, 2L]
  uv_ok <- is.finite(dist_from[u0]) & is.finite(dist_to[v0]) &
    (dist_from[u0] + 1 + dist_to[v0] == d)
  vu_ok <- is.finite(dist_from[v0]) & is.finite(dist_to[u0]) &
    (dist_from[v0] + 1 + dist_to[u0] == d)
  keep <- uv_ok | vu_ok
  if (!any(keep)) {
    return(NULL)
  }

  dag_from <- integer(sum(keep))
  dag_to <- integer(sum(keep))
  dag_sup <- numeric(sum(keep))
  idx <- 1L
  for (i in which(keep)) {
    if (uv_ok[i]) {
      dag_from[idx] <- u0[i]
      dag_to[idx] <- v0[i]
    } else {
      dag_from[idx] <- v0[i]
      dag_to[idx] <- u0[i]
    }
    dag_sup[idx] <- igraph::E(g)$n_support[i]
    idx <- idx + 1L
  }

  # Layer-by-layer DP: best_support[node] = max total support from source
  n_v <- length(vertex_names)
  best_support <- rep(-Inf, n_v)
  best_support[from_idx] <- 0
  best_pred <- rep(NA_integer_, n_v)

  # Deterministic tie-breaking: when support is equal, prefer
  # lexicographically smaller predecessor vertex name.
  layer_order <- order(dist_from)
  for (vi in layer_order) {
    if (!is.finite(best_support[vi]) || best_support[vi] < 0) next
    if (dist_from[vi] >= d) next
    out_mask <- dag_from == vi
    if (!any(out_mask)) next
    nbrs <- dag_to[out_mask]
    sups <- dag_sup[out_mask]
    for (j in seq_along(nbrs)) {
      nbr <- nbrs[j]
      cand <- best_support[vi] + sups[j]
      better <- cand > best_support[nbr] ||
        (cand == best_support[nbr] &&
          vertex_names[vi] < vertex_names[best_pred[nbr]])
      if (better) {
        best_support[nbr] <- cand
        best_pred[nbr] <- vi
      }
    }
  }

  if (!is.finite(best_support[to_idx]) || best_support[to_idx] < 0) {
    return(NULL)
  }

  # Backtrack with cycle / broken-chain defense
  cur <- to_idx
  path_rev <- cur
  while (cur != from_idx) {
    cur <- best_pred[cur]
    if (is.na(cur)) {
      return(NULL)
    }
    path_rev <- c(path_rev, cur)
    if (length(path_rev) > d + 1L) {
      stop("DAG predecessor cycle detected.", call. = FALSE)
    }
  }
  path_nodes <- rev(path_rev)
  if (length(path_nodes) - 1L != as.integer(d)) {
    stop("Selected path does not match minimum-hop distance.", call. = FALSE)
  }
  igraph::V(g)[path_nodes]
}

.build_target_gene_links <- function(
  hit_df, bed_info, loop_annotation_final, map_info, ego_list_target,
  g = NULL, neighbor_hop = 0L
) {
  base_cols <- c(
    "input_id", "loop_ID", "anchor_id", "gene", "gene_role",
    "strict_assignment_eligible",
    "source", "evidence", "anchor_role", "path_length"
  )
  empty <- .empty_target_gene_links()[, base_cols]
  gene_map <- .target_anchor_gene_map(map_info)
  linear_rows <- .linear_target_gene_links(bed_info)

  if (.is_null_or_empty(hit_df) || nrow(gene_map) == 0) {
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

  # Lazy target ego: compute 2-hop neighbourhoods only for anchors
  # that are actually hit by the target BED, avoiding precomputation
  # for every graph node.
  if (neighbor_hop == 1L && !is.null(g) && nrow(hit_base) > 0L) {
    if (is.null(ego_list_target)) ego_list_target <- list()
    cached_ids <- names(ego_list_target)
    if (is.null(cached_ids)) cached_ids <- character(0)
    required_ids <- unique(hit_base$local_anchor_id)
    missing_ids <- setdiff(required_ids, cached_ids)
    valid_ids <- intersect(missing_ids, igraph::V(g)$name)
    if (length(valid_ids) > 0L) {
      lazy_ego <- igraph::ego(g, order = 2L, nodes = valid_ids, mode = "all")
      names(lazy_ego) <- valid_ids
      ego_list_target <- c(ego_list_target, lazy_ego)
    }
  }

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
      strict_assignment_eligible,
      source = "loop_anchor",
      evidence = dplyr::case_when(
        gene_role == "promoter"           ~ "local_promoter_overlap",
        gene_role == "gene_body"          ~ "gene_body_context",
        gene_role == "enhancer_candidate"  ~ "local_enhancer_candidate",
        TRUE                              ~ "positional_candidate"
      ),
      anchor_role = "local_anchor",
      path_length = 0L
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
      strict_assignment_eligible,
      source = "loop_anchor",
      evidence = dplyr::case_when(
        gene_role == "promoter"           ~ "distal_promoter",
        gene_role == "gene_body"          ~ "distal_gene_body_context",
        gene_role == "enhancer_candidate"  ~ "distal_enhancer_candidate",
        TRUE                              ~ "positional_candidate"
      ),
      anchor_role = "opposite_anchor",
      path_length = 1L
    )

  # Pre-build direct neighbour index: avoids scanning edge_long for each hit
  direct_index <- split(edge_long$opposite_anchor_id, edge_long$local_anchor_id)
  expanded_seed <- do.call(rbind, lapply(seq_len(nrow(hit_base)), function(i) {
    local_id <- hit_base$local_anchor_id[[i]]
    ego_ids <- if (local_id %in% names(ego_list_target)) {
      .igraph_vertex_ids(ego_list_target[[local_id]], g)
    } else {
      character(0)
    }
    direct_ids <- direct_index[[local_id]]
    if (is.null(direct_ids)) direct_ids <- character(0)
    expanded_ids <- setdiff(ego_ids, c(local_id, direct_ids))
    if (length(expanded_ids) == 0) {
      return(NULL)
    }
    data.frame(
      input_id = hit_base$input_id[[i]],
      local_anchor_id = hit_base$local_anchor_id[[i]],
      anchor_id = expanded_ids,
      stringsAsFactors = FALSE
    )
  }))
  expanded_rows <- empty
  if (!is.null(expanded_seed) && nrow(expanded_seed) > 0) {
    # Pre-filter: drop expanded anchors without gene assignments.
    # This avoids path computation for anchors that cannot produce
    # target genes, reducing the path loop size.
    gene_anchor_ids <- gene_map %>%
      dplyr::filter(!is.na(.data$gene), .data$gene != "") %>%
      dplyr::distinct(.data$anchor_id)
    expanded_seed <- expanded_seed %>%
      dplyr::semi_join(gene_anchor_ids, by = "anchor_id")
    if (nrow(expanded_seed) == 0) {
      expanded_seed <- NULL
    }
  }
  if (!is.null(expanded_seed) && nrow(expanded_seed) > 0) {
    # Compute shortest-path info BEFORE joining gene_map so that the
    # one-row-per-anchor alignment is unambiguous.  Joining gene_map
    # afterwards lets each gene row inherit its anchor's path fields.
    if (!is.null(g)) {
      # Build undirected pair ->loop_ID lookup from the edge table
      pair_lookup <- dplyr::bind_rows(
        edge_long %>%
          dplyr::transmute(
            pair = paste(local_anchor_id, opposite_anchor_id, sep = "|||"),
            loop_ID = loop_ID
          ),
        edge_long %>%
          dplyr::transmute(
            pair = paste(opposite_anchor_id, local_anchor_id, sep = "|||"),
            loop_ID = loop_ID
          )
      ) %>%
        dplyr::group_by(pair) %>%
        dplyr::summarise(
          loop_IDs = .collapse_genes(loop_ID),
          .groups = "drop"
        )
      pair_map <- setNames(pair_lookup$loop_IDs, pair_lookup$pair)

      # Lexicographic path optimisation:
      #   1. Minimise hop count (unweighted shortest distance)
      #   2. Among minimum-hop paths, prefer higher edge support.
      # For d <= 2, enumerate all shortest paths (safe in sparse
      # Hi-C graphs).  For d >= 3, use a single weighted path to
      # avoid combinatorial explosion.
      #
      # We do NOT use all_shortest_paths(weights=...) because
      # weighted shortest paths may have different hop counts;
      # hop count must be the primary objective.
      max_sp <- 200L
      h1 <- neighbor_hop == 1L
      # Deduplicate (local_anchor_id, anchor_id) pairs: multiple target
      # peaks may hit the same local anchor, producing identical graph
      # node pairs.  Compute paths once, then join back.
      unique_pairs <- expanded_seed %>%
        dplyr::distinct(.data$local_anchor_id, .data$anchor_id)
      # Cache first-order neighbours for all nodes referenced by
      # unique_pairs, avoiding repeated igraph::neighbors() calls
      # when a single node appears in many pairs.
      needed_nodes <- unique(c(unique_pairs$local_anchor_id, unique_pairs$anchor_id))
      neighbor_cache <- stats::setNames(
        lapply(needed_nodes, function(node_id) {
          .igraph_vertex_ids(igraph::neighbors(g, node_id, mode = "all"), g)
        }),
        needed_nodes
      )
      path_info <- lapply(seq_len(nrow(unique_pairs)), function(j) {
        from <- unique_pairs$local_anchor_id[j]
        to <- unique_pairs$anchor_id[j]
        # Unweighted distance: determines minimum hop count.
        # When neighbor_hop = 1, all expanded candidates are
        # exactly 2 hops away (2-hop ego minus direct neighbours).
        d <- if (h1) 2L else tryCatch(
          igraph::distances(g,
            v = from, to = to, mode = "all",
            weights = NA
          )[1, 1],
          error = function(e) Inf
        )
        if (is.infinite(d) || is.na(d) || d < 1) {
          return(c(loop_ID = NA_character_, path_length = NA_integer_))
        }
        if (d == 2L) {
          # 2-hop: each common neighbour forms one path.
          # No need to call all_shortest_paths(); O(degree)
          # is safe and avoids the enumeration risk entirely.
          nb_from <- neighbor_cache[[from]]
          nb_to <- neighbor_cache[[to]]
          common <- intersect(nb_from, nb_to)
          if (length(common) == 0L) {
            return(c(loop_ID = NA_character_, path_length = NA_integer_))
          }
          # Score ALL common neighbours by edge support.
          # O(degree) -- no need to truncate before scoring.
          .mid_support <- function(mid) {
            vp <- as.character(c(from, mid, mid, to))
            names(vp) <- NULL
            eids <- igraph::get_edge_ids(g,
              vp = vp,
              directed = FALSE, error = FALSE
            )
            if (any(eids == 0L)) {
              return(NA_real_)
            }
            sum(igraph::E(g)$n_support[eids], na.rm = TRUE)
          }
          supports <- vapply(common, .mid_support, numeric(1))
          best_sup <- max(supports, na.rm = TRUE)
          best_mids <- common[!is.na(supports) & supports == best_sup]
          # Stable deterministic order: sort by anchor ID.
          best_mids <- sort(best_mids)
          if (length(best_mids) > max_sp) {
            warning("d=2 path enumeration truncated: ",
              length(best_mids), " optimal midpoints found, ",
              "keeping first ", max_sp, ". ",
              "Results are deterministic (sorted by anchor ID).",
              call. = FALSE
            )
            best_mids <- best_mids[seq_len(max_sp)]
          }
          if (length(best_mids) == 0L) {
            best <- list(igraph::shortest_paths(
              g,
              from = from, to = to, mode = "all"
            )$vpath[[1]])
          } else {
            # Keep only the lexicographically first best mid
            # (sorted above).  Consistent with d>=3 which also
            # returns a single deterministic path.  loop_ID
            # now always represents one path's edge set.
            ids <- as.character(c(from, best_mids[1L], to))
            names(ids) <- NULL
            best <- list(igraph::V(g)[ids])
          }
        } else {
          # d >= 3: use shortest-path DAG + DP to find the
          # minimum-hop maximum-support path.  O(V+E) --
          # avoids the combinatorial enumeration of
          # all_shortest_paths().
          best_path <- .dag_max_support_path(g, from, to, d)
          if (is.null(best_path) || length(best_path) < 2L) {
            return(c(loop_ID = NA_character_, path_length = NA_integer_))
          }
          best <- list(best_path)
        }
        n_hops <- length(best[[1]]) - 1L
        if (!identical(as.integer(n_hops), as.integer(d))) {
          stop("Internal path invariant failed: selected path length (",
            n_hops, ") differs from unweighted shortest distance (",
            d, ").",
            call. = FALSE
          )
        }
        all_lids <- character(0)
        for (sp in best) {
          nodes <- .igraph_vertex_ids(sp, g)
          pair_keys <- paste(nodes[seq_len(n_hops)], nodes[-1], sep = "|||")
          all_lids <- c(all_lids, unname(pair_map[pair_keys]))
        }
        lid_str <- .collapse_genes(all_lids)
        if (is.na(lid_str) || lid_str == "") lid_str <- NA_character_
        c(loop_ID = lid_str, path_length = as.character(n_hops))
      })
      unique_pairs$loop_ID_path <- vapply(path_info, `[`, character(1), "loop_ID")
      unique_pairs$path_length <- as.integer(vapply(path_info, `[`, character(1), "path_length"))
      expanded_seed <- expanded_seed %>%
        dplyr::left_join(unique_pairs, by = c("local_anchor_id", "anchor_id"))
    } else {
      expanded_seed$loop_ID_path <- NA_character_
      expanded_seed$path_length <- NA_integer_
    }
    # Now join gene_map -- each gene row inherits its anchor's path columns
    expanded_rows <- expanded_seed %>%
      dplyr::left_join(gene_map, by = "anchor_id") %>%
      dplyr::filter(!is.na(gene), gene != "") %>%
      dplyr::transmute(
        input_id,
        loop_ID = loop_ID_path,
        anchor_id,
        gene,
        gene_role,
        strict_assignment_eligible,
        source = "loop_anchor",
        evidence = dplyr::case_when(
          gene_role == "promoter"            ~ "expanded_promoter_loop",
          gene_role == "gene_body"           ~ "expanded_gene_body_context",
          gene_role == "enhancer_candidate"   ~ "expanded_enhancer_candidate",
          gene_role == "positional_candidate" ~ "expanded_positional_candidate",
          TRUE                                ~ "expanded_anchor"
        ),
        anchor_role = "expanded_anchor",
        path_length = path_length
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
  if (.is_null_or_empty(target_gene_links)) {
    return(data.frame(input_id = character(), Regulated_promoter_Evidence = character()))
  }
  target_gene_links %>%
    dplyr::mutate(
      strict_assignment_eligible = if ("strict_assignment_eligible" %in% colnames(target_gene_links)) {
        target_gene_links$strict_assignment_eligible
      } else {
        gene_role %in% c("promoter", "gene_body")
      }
    ) %>%
    dplyr::filter(
      source == "loop_anchor",
      gene_role == "promoter",
      strict_assignment_eligible %in% TRUE,
      is.finite(path_length),
      path_length <= 1
    ) %>%
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
  if (.is_null_or_empty(target_gene_links)) {
    return(.empty_target_gene_links())
  }

  target_cols <- c(
    "Regulated_promoter_genes", "Assigned_Target_Genes",
    "All_Loop_Connected_Genes", "Regulated_promoter_genes_Filled",
    "Assigned_Target_Genes_Filled", "Expanded_Target_Genes"
  )
  lookup_cols <- intersect(c("input_id", target_cols), colnames(bed_info))
  link_df <- target_gene_links %>%
    dplyr::select(-dplyr::any_of(c(
      "used_as_fallback", "in_regulated_promoter", "in_assigned_target",
      "in_all_loop_connected", "in_regulated_promoter_filled",
      "in_assigned_target_filled", "in_expanded_target"
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
  link_df$in_expanded_target <- unlist(Map(.contains_target_gene, link_df$gene, link_df$Expanded_Target_Genes), use.names = FALSE)

  # Compute is_expanded_path from path_length for ALL pipeline stages
  if ("path_length" %in% colnames(link_df)) {
    link_df$is_expanded_path <- dplyr::case_when(
      is.na(link_df$path_length) ~ NA,
      is.finite(link_df$path_length) & link_df$path_length > 1 ~ TRUE,
      is.finite(link_df$path_length) ~ FALSE,
      TRUE ~ NA
    )
  }

  link_df$evidence[link_df$source == "linear_annotation"] <- "linear_annotation"
  reg_empty <- is.na(link_df$Regulated_promoter_genes) | link_df$Regulated_promoter_genes == ""
  assigned_empty <- is.na(link_df$Assigned_Target_Genes) | link_df$Assigned_Target_Genes == ""
  link_df$used_as_fallback <- link_df$source == "linear_annotation" &
    ((reg_empty & link_df$in_regulated_promoter_filled) |
      (assigned_empty & link_df$in_assigned_target_filled))
  link_df$evidence[link_df$used_as_fallback] <- "linear_fallback"

  .ensure_target_link_schema(link_df) %>% dplyr::distinct()
}

#' Internal: Filter Target Gene Links After Expression-Aware Refinement
#'
#' Re-marks membership flags against expression-filtered target columns,
#' appends \code{Mean_Expression} and \code{Passes_Expression_Filter},
#' and retains all structural links while marking which rows contribute
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
    "Passes_Expression_Filter", "Expression_State",
    "gene_role", "source", "evidence",
    "anchor_role", "path_length",
    # Preserve chromatin provenance when chromatin refinement ran
    # before expression refinement (preserved for forward-compatibility
    # even though the reverse order is now unsupported).
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
    "in_assigned_target_filled",
    "retained_after_refinement",
    "strict_assignment_eligible",
    "in_expanded_target",
    "is_expanded_path",
    "TSS_supported",
    "TSS_support_status",
    "Gene_Assignment_Confidence",
    "Gene_Assignment_Evidence"
  )
  empty_refined <- function() {
    out <- .empty_target_gene_links()
    out$Mean_Expression <- numeric()
    out$Passes_Expression_Filter <- logical()
    out$Expression_State <- character()
    out[, base::intersect(refined_cols, colnames(out)), drop = FALSE]
  }

  link_df <- .mark_target_gene_link_membership(target_gene_links, bed_info)
  if (nrow(link_df) == 0) {
    return(empty_refined())
  }

  # Resolve expression state per gene: active / measured_silent / unmeasured
  gene_keys <- toupper(link_df$gene)
  expr <- unname(vals[gene_keys])
  link_df$Mean_Expression <- as.numeric(expr)

  measured <- gene_keys %in% toupper(names(vals))
  passes <- .passes_expression_threshold(link_df$Mean_Expression, threshold)

  link_df$Expression_State <- dplyr::case_when(
    !measured ~ "unmeasured",
    measured & passes ~ "active",
    measured & !passes ~ "measured_silent",
    TRUE ~ "unmeasured"
  )
  # Passes_Expression_Filter: TRUE=active, FALSE=measured_silent, NA=unmeasured
  link_df$Passes_Expression_Filter <- dplyr::case_when(
    link_df$Expression_State == "active" ~ TRUE,
    link_df$Expression_State == "measured_silent" ~ FALSE,
    TRUE ~ NA
  )

  link_df$retained_after_refinement <- link_df$in_regulated_promoter |
    link_df$in_assigned_target |
    link_df$in_all_loop_connected |
    link_df$in_regulated_promoter_filled |
    link_df$in_assigned_target_filled

  # Keep ALL structural links -- silent and unmeasured links retain their
  # expression provenance.  Summary columns are filtered elsewhere.
  # The `retained_after_refinement` flag marks links that enter at least
  # one current summary target column.
  link_df %>%
    .ensure_target_link_schema() %>%
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
  if (.is_null_or_empty(target_connected_loops)) {
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
  threshold, unit_type, log_message,
  threshold_mode = c("absolute", "quantile")
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
    "Putative_Target_Genes", "Promoter_Target_Genes"
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

  if (is.null(expr_matrix_file)) {
    stop("`expr_matrix_file` is required for expression refinement.", call. = FALSE)
  }
  vals <- load_expression_matrix(expr_matrix_file, sample_columns)
  threshold_mode <- match.arg(threshold_mode)
  if (threshold_mode == "quantile") {
    if (!is.numeric(threshold) || length(threshold) != 1L ||
      !is.finite(threshold) || threshold < 0 || threshold > 1) {
      stop("For quantile mode, `threshold` must be a single finite number in [0, 1].",
        call. = FALSE
      )
    }
    finite_vals <- vals[is.finite(vals)]
    if (length(finite_vals) == 0L) {
      stop("No finite expression values available for quantile threshold calculation.",
        call. = FALSE
      )
    }
    effective_threshold <- stats::quantile(finite_vals, threshold, na.rm = TRUE, names = FALSE)
    log_message(sprintf(
      "    >>> Quantile mode: threshold=%.2f -> effective expression cutoff = %.3f %s",
      threshold, effective_threshold, unit_type
    ))
  } else {
    if (!is.numeric(threshold) || length(threshold) != 1L ||
      !is.finite(threshold) || threshold < 0) {
      stop("`threshold` must be a single finite non-negative number.",
        call. = FALSE
      )
    }
    effective_threshold <- threshold
  }
  passes <- .passes_expression_threshold(vals, effective_threshold)
  whitelist <- names(vals)[passes & names(vals) != ""]
  log_message(sprintf(
    "    >>> Active Genes (>= %s %s%s): %d",
    if (threshold_mode == "quantile") sprintf("%.3f (%.0f%%ile)", effective_threshold, threshold * 100) else as.character(effective_threshold),
    unit_type,
    if (threshold_mode == "quantile") "" else "",
    length(whitelist)
  ))

  anno_genes <- unique(c(
    trimws(unlist(strsplit(na.omit(loop_df$anchor1_gene), ";"))),
    trimws(unlist(strsplit(na.omit(loop_df$anchor2_gene), ";")))
  ))
  anno_genes <- anno_genes[nzchar(anno_genes)]
  if (length(anno_genes) > 0) {
    # ID mapping diagnostic: use ALL measured genes, not just the
    # expression-filtered whitelist.  Low mapping rate indicates a
    # genuine identifier mismatch (e.g. SYMBOL vs Ensembl), whereas
    # low active rate may simply reflect biological silencing.
    measured_genes <- names(vals)
    anno_genes_upper <- toupper(anno_genes)
    n_measured_mapped <- length(intersect(toupper(measured_genes), anno_genes_upper))
    id_mapping_rate <- n_measured_mapped / length(anno_genes)
    n_active_mapped <- length(intersect(toupper(whitelist), anno_genes_upper))
    active_fraction <- n_active_mapped / length(anno_genes)
    unmeasured_fraction <- 1 - id_mapping_rate
    log_message(sprintf(
      "    >>> ID mapping: %.1f%% (%d / %d) annotation genes found in expression matrix",
      id_mapping_rate * 100, n_measured_mapped, length(anno_genes)
    ))
    log_message(sprintf(
      "    >>> Active fraction: %.1f%% (%d / %d) annotation genes pass expression threshold",
      active_fraction * 100, n_active_mapped, length(anno_genes)
    ))
    if (unmeasured_fraction > 0) {
      log_message(sprintf(
        "    >>> Unmeasured: %.1f%% (%d / %d) annotation genes absent from expression matrix",
        unmeasured_fraction * 100, length(anno_genes) - n_measured_mapped, length(anno_genes)
      ))
    }
    if (id_mapping_rate < 0.5) {
      warning(
        sprintf(
          "Only %.1f%% (%d / %d) of annotation gene symbols match the expression matrix row names. ",
          id_mapping_rate * 100, n_measured_mapped, length(anno_genes)
        ),
        "Low ID mapping suggests different gene identifier conventions. ",
        "Check that expression matrix row names and OrgDb annotations use the same convention ",
        "(e.g., both SYMBOL, or both ENSEMBL). ",
        "If the expression matrix uses Ensembl IDs or Entrez IDs while annotations use SYMBOL, ",
        "convert identifiers first (e.g., via AnnotationDbi::mapIds or org.*.eg.db). ",
        call. = FALSE
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
    measured_set = names(vals),
    effective_threshold = effective_threshold,
    has_cluster_id = has_cluster_id,
    upstream_promoter_stats = upstream_promoter_stats
  )
}

#' Internal: Reclassify anchors and compute expression-filtered targets
#' @return The updated loop_df.
#' @keywords internal
#' @noRd
.refine_apply_anchor_updates <- function(
  loop_df, whitelist, reclassify_by_expression, log_message,
  annotation_res = NULL, measured_set = NULL
) {
  log_message(">>> [Step 2] Updating Anchors & Loops...")

  # Preserve the true original anchor types before this (or any prior)
  # refine run. Using loop_df$anchor_type directly would degrade
  # Refinement_Action labels on re-run (already-eP types would not be
  # detected as reclassified). The anchor_state stores the annotation-time
  # type permanently in original_type_code.
  anchor_state <- attr(annotation_res, "looplook_anchor_state", exact = TRUE)
  if (!is.null(anchor_state) && "map_info" %in% names(anchor_state) &&
    "original_type_code" %in% colnames(anchor_state$map_info)) {
    type_lookup <- setNames(
      anchor_state$map_info$original_type_code,
      anchor_state$map_info$anchor_id
    )
    orig_anchor1_type <- if ("a1_id" %in% colnames(loop_df)) {
      type_lookup[as.character(loop_df$a1_id)]
    } else {
      loop_df$anchor1_type
    }
    orig_anchor2_type <- if ("a2_id" %in% colnames(loop_df)) {
      type_lookup[as.character(loop_df$a2_id)]
    } else {
      loop_df$anchor2_type
    }
    # Fall back to current type for anchors not in lookup
    orig_anchor1_type[is.na(orig_anchor1_type)] <- loop_df$anchor1_type[is.na(orig_anchor1_type)]
    orig_anchor2_type[is.na(orig_anchor2_type)] <- loop_df$anchor2_type[is.na(orig_anchor2_type)]
  } else {
    orig_anchor1_type <- loop_df$anchor1_type
    orig_anchor2_type <- loop_df$anchor2_type
  }
  loop_df <- .refine_reclassify_anchors(
    loop_df, whitelist, reclassify_by_expression,
    measured_set
  )
  original_ptg <- loop_df$Putative_Target_Genes
  original_promoter_targets <- if ("Promoter_Target_Genes" %in% colnames(loop_df)) {
    loop_df$Promoter_Target_Genes
  } else {
    NULL
  }

  loop_df <- .refine_compute_targets(
    loop_df, original_ptg, whitelist,
    orig_anchor1_type, orig_anchor2_type,
    measured_set = measured_set,
    original_promoter_targets = original_promoter_targets
  )

  n_active <- sum(loop_df$Has_Active_Target %in% TRUE)
  n_unknown <- sum(is.na(loop_df$Has_Active_Target))
  log_message(sprintf(
    "    Active: %d; unknown: %d; total: %d loops",
    n_active, n_unknown, nrow(loop_df)
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
        "eP/eG = element-Promoter like / element-Genebody like: structurally defined anchors near TSS or within gene bodies. ",
        "Transcriptional silence does not imply enhancer activity. ",
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
.refine_reclassify_anchors <- function(loop_df, whitelist, reclassify_by_expression,
                                       measured_set = NULL) {
  a1_res <- Map(
    function(g, t) {
      clean_anchor(g, t,
        allow = whitelist, down = reclassify_by_expression,
        measured = measured_set
      )
    },
    loop_df$anchor1_gene, loop_df$anchor1_type
  )
  a2_res <- Map(
    function(g, t) {
      clean_anchor(g, t,
        allow = whitelist, down = reclassify_by_expression,
        measured = measured_set
      )
    },
    loop_df$anchor2_gene, loop_df$anchor2_type
  )
  loop_df$anchor1_type <- vapply(a1_res, function(x) x$type, character(1))
  loop_df$anchor1_gene <- vapply(a1_res, function(x) x$gene, character(1))
  loop_df$anchor2_type <- vapply(a2_res, function(x) x$type, character(1))
  loop_df$anchor2_gene <- vapply(a2_res, function(x) x$gene, character(1))
  loop_df$loop_type <- unlist(Map(
    function(t1, t2) paste(c(t1, t2)[order(tolower(c(t1, t2)), c(t1, t2), method = "radix")], collapse = "-"),
    loop_df$anchor1_type, loop_df$anchor2_type
  ), use.names = FALSE)
  loop_df
}

#' Compute active/fallback target genes and refinement status flags.
#' @keywords internal
#' @noRd
.refine_compute_targets <- function(loop_df, original_ptg, whitelist,
                                    orig_anchor1_type = NULL,
                                    orig_anchor2_type = NULL,
                                    measured_set = NULL,
                                    original_promoter_targets = NULL) {
  whitelist_upper <- toupper(whitelist)
  filter_genes_wl <- function(x) {
    if (is.na(x) || x == "") {
      return(NA_character_)
    }
    gs <- trimws(unlist(strsplit(as.character(x), ";")))
    gs <- unique(gs[toupper(gs) %in% whitelist_upper])
    if (length(gs) == 0) {
      return(NA_character_)
    }
    paste(sort(gs), collapse = ";")
  }
  .target_expr_state <- function(gene_string) {
    gs <- clean_gene_names(gene_string, ";")
    if (length(gs) == 0L) {
      return("no_target")
    }
    keys <- toupper(gs)
    measured <- keys %in% toupper(measured_set)
    active <- keys %in% whitelist_upper
    n_measured <- sum(measured)
    n_active <- sum(active)
    n_total <- length(gs)
    if (n_active > 0L) {
      return(if (n_measured == n_total) "active" else "active_partial")
    }
    if (n_measured == 0L) {
      return("unmeasured")
    }
    if (n_measured == n_total) {
      return("measured_silent")
    }
    return("mixed")
  }

  # Compute Target_Expression_State on PRE-filter original PTG
  if (is.null(measured_set)) {
    loop_df$Target_Expression_State <- "not_assessed"
  } else {
    # measured_set may be character(0) when no gene has any finite
    # measurement; every target is then unmeasured, not not_assessed.
    loop_df$Target_Expression_State <- vapply(
      original_ptg, .target_expr_state, character(1)
    )
  }

  # Filter current target sets by expression -- no fallback re-introduction.
  # All remaining genes must pass the expression threshold.
  loop_df$Putative_Target_Genes <- vapply(
    original_ptg, filter_genes_wl, character(1)
  )

  if (!is.null(original_promoter_targets)) {
    loop_df$Promoter_Target_Genes <- vapply(
      original_promoter_targets, filter_genes_wl, character(1)
    )
  } else if ("Promoter_Target_Genes" %in% colnames(loop_df)) {
    loop_df$Promoter_Target_Genes <- vapply(
      loop_df$Promoter_Target_Genes, filter_genes_wl, character(1)
    )
  } else {
    stop(
      "Input annotation lacks 'Promoter_Target_Genes' column. ",
      "Re-run annotate_peaks_and_loops() with the current ",
      "pipeline version to generate this column.",
      call. = FALSE
    )
  }

  has_active <- !is.na(loop_df$Putative_Target_Genes) & loop_df$Putative_Target_Genes != ""
  unknown_idx <- loop_df$Target_Expression_State %in% c("unmeasured", "mixed")
  loop_df$Has_Active_Target <- has_active
  loop_df$Has_Active_Target[unknown_idx] <- NA

  if (!is.null(orig_anchor1_type) && !is.null(orig_anchor2_type)) {
    reclassified_a1 <- orig_anchor1_type != loop_df$anchor1_type & loop_df$anchor1_type %in% c("eP", "eG")
    reclassified_a2 <- orig_anchor2_type != loop_df$anchor2_type & loop_df$anchor2_type %in% c("eP", "eG")
    has_reclassified <- reclassified_a1 | reclassified_a2
  } else {
    has_reclassified <- loop_df$anchor1_type %in% c("eP", "eG") | loop_df$anchor2_type %in% c("eP", "eG")
  }
  had_original <- !is.na(original_ptg) & original_ptg != ""

  loop_df$Refinement_Action <- dplyr::case_when(
    loop_df$Has_Active_Target %in% TRUE ~ "retained_active_target",
    is.na(loop_df$Has_Active_Target) ~ "expression_state_unknown",
    loop_df$Has_Active_Target %in% FALSE & has_reclassified & had_original ~ "reclassified_silent_anchor",
    loop_df$Has_Active_Target %in% FALSE & had_original ~ "expression_filtered_no_active_target",
    TRUE ~ "structural_only_no_active_target"
  )
  loop_df$Retained_In_Functional_Network <- loop_df$Has_Active_Target
  loop_df
}

#' Refine target BED annotations by expression filtering.
#' @keywords internal
#' @noRd
.refine_target_annotations <- function(bed_info, loop_df, whitelist, target_gene_links, vals, threshold) {
  whitelist_upper <- toupper(whitelist)
  cols_to_clean <- base::intersect(c(
    "All_Loop_Connected_Genes",
    "Regulated_promoter_genes",
    "Assigned_Target_Genes",
    "All_Loop_Connected_Genes_Filled",
    "Regulated_promoter_genes_Filled",
    "Assigned_Target_Genes_Filled"
  ), colnames(bed_info))
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
        gs_active <- gs[toupper(trimws(gs)) %in% whitelist_upper]
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
    any(toupper(trimws(unlist(strsplit(as.character(s), ";")))) %in% whitelist_upper)
  }
  filter_sym_expressed <- function(s) {
    if (is.na(s) || s == "") {
      return(NA_character_)
    }
    gs <- trimws(unlist(strsplit(as.character(s), ";")))
    gs <- gs[toupper(gs) %in% whitelist_upper]
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
      # Recomputed evidence when using linear fallback after
      # expression filtering.  The pre-expression evidence may
      # have been "none" (strict was non-empty); now that strict
      # is empty and we fall back to the linear annotation gene,
      # the evidence must reflect the annotation type.
      fallback_evidence <- .fallback_evidence_from_annotation(
        bed_info$annotation
      )
      bed_info$Regulated_promoter_Fallback_Evidence <- dplyr::case_when(
        has_reg ~ "none",
        has_sym ~ fallback_evidence,
        TRUE ~ "none"
      )
    }
  }

  # NOTE (kept for reference, currently INERT): this loop-connected-gene fill
  # block is overwritten by .finalize_current_target_annotation(), which runs
  # after .refine_target_annotations() in the expression-refinement pipeline
  # and unconditionally clears + re-aggregates all strict/Filled columns via
  # .fill_fallback_targets().  The final *Filled columns are therefore always
  # one-level fallbacks (strict -> expression-filtered linear gene), and the
  # loop-filled genes written below never survive to the output.  It is
  # retained only as documentation of the intended two-level design; remove
  # or relocate it (into .fill_fallback_targets) if two-level filling is
  # ever required.
  if ("Linked_Loop_IDs" %in% colnames(bed_info) &&
    "loop_ID" %in% colnames(loop_df) &&
    "Putative_Target_Genes" %in% colnames(loop_df)) {
    loop_tgt <- loop_df %>%
      dplyr::filter(!is.na(Putative_Target_Genes) & Putative_Target_Genes != "") %>%
      dplyr::select(loop_ID, Putative_Target_Genes) %>%
      dplyr::distinct()
    get_loop_tgt <- function(linked) {
      if (is.na(linked) || linked == "") {
        return(NA_character_)
      }
      ids <- trimws(unlist(strsplit(as.character(linked), ";")))
      tgt <- loop_tgt$Putative_Target_Genes[match(ids, loop_tgt$loop_ID)]
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

  # Rebuild All_Loop_Connected_Genes_Filled alongside the other two Filled
  # so that all three columns use the same expression-filtered fallback.
  if ("All_Loop_Connected_Genes_Filled" %in% colnames(bed_info)) {
    has_all <- !is.na(bed_info$All_Loop_Connected_Genes) & bed_info$All_Loop_Connected_Genes != ""
    bed_info$All_Loop_Connected_Genes_Filled <- dplyr::case_when(
      has_all ~ bed_info$All_Loop_Connected_Genes,
      has_sym ~ sym_fill,
      TRUE ~ NA_character_
    )
  }

  # Rebuild Regulated_promoter_Evidence from current gene set.
  # After expression filtering, evidence must only explain the genes
  # that remain.  We defer the full rebuild until after link filtering
  # so that .summarise_regulated_promoter_evidence() can use
  # Passes_Expression_Filter to select only current promoter links.
  if ("Regulated_promoter_Evidence" %in% colnames(bed_info) &&
    "Regulated_promoter_genes" %in% colnames(bed_info)) {
    has_promoter <- !is.na(bed_info$Regulated_promoter_genes) &
      bed_info$Regulated_promoter_genes != ""
    bed_info$Regulated_promoter_Evidence <- ifelse(
      has_promoter,
      bed_info$Regulated_promoter_Evidence,
      "none"
    )
  }

  if (!is.null(target_gene_links)) {
    target_gene_links <- .filter_refined_target_gene_links(
      target_gene_links, bed_info, vals, threshold
    )

    # Rebuild evidence from current (Passes==TRUE) promoter links.
    # This correctly handles partial filtering: when some promoter
    # genes are filtered out, their evidence types are also removed.
    if ("Regulated_promoter_Evidence" %in% colnames(bed_info)) {
      current_promoter_links <- target_gene_links %>%
        dplyr::filter(
          source == "loop_anchor",
          gene_role == "promoter",
          Passes_Expression_Filter %in% TRUE
        )
      if (nrow(current_promoter_links) > 0) {
        evidence_df <- .summarise_regulated_promoter_evidence(
          current_promoter_links
        )
        bed_info <- bed_info %>%
          dplyr::select(-dplyr::any_of("Regulated_promoter_Evidence")) %>%
          dplyr::left_join(evidence_df, by = "input_id")
      }
      bed_info$Regulated_promoter_Evidence <- ifelse(
        is.na(bed_info$Regulated_promoter_Evidence) |
          bed_info$Regulated_promoter_Evidence == "",
        "none",
        bed_info$Regulated_promoter_Evidence
      )
    }
  }

  list(bed_info = bed_info, target_gene_links = target_gene_links)
}

#' Export refined results to Excel workbook.
#' @keywords internal
#' @noRd
.refine_export_workbook <- function(
  loop_df, clust_info, promoter_centric_df,
  distal_element_df, bed_info, target_gene_links, out_dir, project_name,
  chromatin_validation = NULL
) {
  wb <- openxlsx::createWorkbook()
  loop_export <- loop_df %>%
    dplyr::select(-any_of(c("a1_id", "a2_id", "loop_genes", "single_loop_genes", "proximate_loop_gene")))

  .add_sheet(wb, "Filtered Loop Annotation", loop_export)

  functional_loops <- loop_export %>% dplyr::filter(Retained_In_Functional_Network)
  .add_sheet(wb, "Functional Loop Annotation", functional_loops)


  .add_sheet(wb, "Filtered Anchor Loci Annotation", clust_info)
  .add_sheet(wb, "Filtered Promoter Stats", promoter_centric_df)
  .add_sheet(wb, "Filtered Distal Element Stats", distal_element_df)
  .add_sheet(wb, "Filtered Target Annotation", bed_info, drop_cols = "SANKEY_RAW_GENES")
  .add_sheet(wb, "Filtered Target Gene Links", target_gene_links, require_rows = TRUE)
  .add_sheet(wb, "Chromatin Validation", chromatin_validation, require_rows = TRUE)

  .save_workbook(wb, out_dir, project_name, "_Refined_Results.xlsx",
    "Failed to save refined Excel workbook: "
  )
}
#' Internal: Validate expression refinement inputs and load expression data
#'
#' Handles workflow-order guard, allow_rerefine guard, project setup,
#' and delegates to .refine_load_validate_data() for expression matrix
#' loading and threshold resolution.
#'
#' @return Named list of validated inputs ready for downstream refinement steps.
#' @keywords internal
#' @noRd
.expression_validate_and_load <- function(
  annotation_res, expr_matrix_file,
  sample_columns, threshold, threshold_mode, unit_type,
  reclassify_by_expression, allow_rerefine, project_name,
  out_dir, write_output, quiet
) {
  log_message <- .make_log_message(quiet)

  if (!grepl("_Filtered$", project_name)) {
    project_name <- paste0(project_name, "_Filtered")
  }
  .ensure_out_dir(write_output, out_dir)
  log_message(">>> [Refinement] Project Name: ", project_name)
  if (isTRUE(reclassify_by_expression)) {
    log_message(">>> NOTE: eP/eG = element-Promoter like / element-Genebody like: structurally")
    log_message(">>> defined anchors near TSS or within gene bodies. The associated gene")
    log_message(">>> may be silent; validate with ATAC-seq, H3K27ac, or chromatin data")
    log_message(">>> for functional regulatory interpretation.")
  }

  up_meta <- annotation_res$metadata

  # Guard: chromatin ->expression is unsupported workflow order.
  if (!is.null(up_meta) && isTRUE(up_meta$function_name == "refine_loop_anchors_by_chromatin")) {
    stop(
      "Unsupported workflow order detected: ",
      "refine_loop_anchors_by_chromatin() was run before ",
      "refine_loop_anchors_by_expression().\n",
      "The only supported multi-omic workflow is:\n",
      "  annotate_peaks_and_loops() ->refine_loop_anchors_by_expression() ",
      "-> refine_loop_anchors_by_chromatin()\n",
      "Please restart from the original annotation object returned by ",
      "annotate_peaks_and_loops().",
      call. = FALSE
    )
  }

  # Guard against re-refinement: repeated expression refinement is not idempotent.
  if (!is.null(up_meta) &&
    isTRUE(up_meta$function_name == "refine_loop_anchors_by_expression")) {
    if (!isTRUE(allow_rerefine)) {
      stop(
        "Expression refinement cannot be applied to an already ",
        "expression-refined object. Anchor types and gene assignments ",
        "have already been modified by the first pass, and a second ",
        "pass would produce results that depend on the calling history. ",
        "For threshold sensitivity analysis, restart from the original ",
        "annotation result returned by annotate_peaks_and_loops(). ",
        "To override this guard (for interactive exploration only), ",
        "set `allow_rerefine = TRUE`.",
        call. = FALSE
      )
    } else {
      warning(
        "Re-refining an already expression-refined object. ",
        "Results depend on the calling history. ",
        "For reproducible sensitivity analysis, restart from the ",
        "original annotation object.",
        call. = FALSE
      )
    }
  }

  refine_data <- .refine_load_validate_data(
    annotation_res, expr_matrix_file, sample_columns,
    threshold, unit_type, log_message,
    threshold_mode = threshold_mode
  )

  list(
    up_meta = up_meta,
    loop_df = refine_data$loop_df,
    original_loop_df = refine_data$original_loop_df,
    clust_info = refine_data$clust_info,
    bed_info = refine_data$bed_info,
    target_gene_links = refine_data$target_gene_links,
    whitelist = refine_data$whitelist,
    vals = refine_data$vals,
    measured_set = refine_data$measured_set,
    effective_threshold = refine_data$effective_threshold,
    has_cluster_id = refine_data$has_cluster_id,
    upstream_promoter_stats = refine_data$upstream_promoter_stats,
    project_name = project_name
  )
}

#' @title Expression-Aware refinement of loop anchors and target linkages
#'
#' @description
#' Integrates quantitative transcription data with 3D structural data.
#' Preferred input order: CAGE-seq / PRO-seq > GRO-seq / NET-seq >
#' TT-seq > RNA-seq (methods that directly measure promoter activity or
#' Pol II elongation provide the most direct evidence for P/eP classification).
#' to reclassify regulatory elements and annotate functional status, deriving a
#' functionally interpretable regulatory network from physical chromatin contacts.
#' All structural loop rows are preserved; anchor type and gene assignments
#' are reclassified in place (silent promoters become eP/eG with
#' positional gene preserved; expression filtering is applied to
#' \code{Putative_Target_Genes} and \code{Promoter_Target_Genes}, not
#' to anchor gene assignments). Refinement status columns indicate which
#' loops belong to the high-confidence active subset.
#'
#' @details
#' \strong{Algorithmic Framework:}
#' \itemize{
#'   \item \strong{Target Filtration:} Parses merged gene assignments (e.g., \code{"GeneA;GeneB"}), evaluates individual genes against a defined expression threshold, and retains only transcriptionally active targets.
#'   \item \strong{Biological Reclassification:} Reclassifies physically annotated promoters (\code{P}) and gene bodies (\code{G}) lacking active transcription as \emph{element-Promoter like} (\code{eP}) and \emph{element-Genebody like} (\code{eG}) elements. \strong{Important:} \code{eP}/\code{eG} labels denote structural genomic position (near a TSS or within a gene body) -- they do \strong{not} constitute evidence of enhancer activity, even though the associated gene is transcriptionally silent. Orthogonal chromatin data (ATAC-seq, H3K27ac, H3K4me1) are required for functional regulatory interpretation. The \code{e} prefix stands for \strong{element} (a structurally defined genomic feature akin to a cis-regulatory element), not enhancer.
#'   \item \strong{Expression-Aware Connectivity Statistics:} Recomputes promoter-centric and distal-element connectivity after expression-aware anchor refinement, while preserving all structural loop rows in the refined loop annotation (anchor types and genes are reclassified in place for silent elements). This separates the complete physical contact map from the high-confidence active subset.
#'   \item \strong{External Target Refinement:} Filters auxiliary target mapping columns (e.g., \code{Assigned_Target_Genes_Filled}) based on expression criteria, ensuring that mapped 1D genomic features are exclusively linked to active genes.
#'   \item \strong{Target Provenance Preservation:} Recomputes \code{*_Filled}
#'   membership flags in \code{target_gene_links} after expression filtering,
#'   retains all structural links with expression provenance, and appends
#'   \code{Mean_Expression}, \code{Expression_State}, and
#'   \code{Passes_Expression_Filter}. Evidence labels
#'   such as \code{local_promoter_overlap}, \code{distal_promoter}, and
#'   \code{linear_fallback} are preserved.
#' }
#'
#' \strong{Design Philosophy:}
#' This function does not discard structural loop rows based on expression state.
#' Hi-C, HiChIP, and PLAC-seq capture 3D chromatin contacts; RNA-seq captures
#' current transcriptional state. A silent promoter may reflect cell-state,
#' time-point, or technical factors rather than absence of physical contact.
#' However, anchor gene assignments and types are reclassified in place:
#' \code{anchor1_gene}/\code{anchor2_gene} may be set to \code{NA} and anchor
#' types downgraded (e.g. \code{P -> eP}, \code{G -> eG}) for silent regulatory
#' elements. All rows are retained; a high-confidence functional subset is
#' provided via \code{Retained_In_Functional_Network} and the
#' \emph{Functional Loop Annotation} Excel sheet.
#'
#' \strong{Interpretation of eP/eG labels:}
#' \code{eP} (\strong{e}lement-\strong{P}romoter like) and \code{eG}
#' (\strong{e}lement-\strong{G}enebody like) are \strong{structural positional
#' classifications}, not functional enhancer classifications. An eP anchor
#' resides near a TSS and is structurally promoter-like; an eG anchor lies
#' within a gene body and is structurally genebody-like. The associated
#' gene may be transcriptionally silent in the current sample -- this can
#' arise from cell-type proportions, time-point effects, sequencing depth,
#' or promoter pausing -- but the anchor retains its positional identity
#' as a putative regulatory element. These labels should be treated as
#' hypotheses requiring orthogonal validation (ATAC-seq, H3K27ac, H3K4me1,
#' or H3K27me3 depletion). Users with matched chromatin data should
#' overlay eP/eG loci against these tracks before interpreting them as
#' functionally active regulatory elements.
#'
#' \strong{Pipeline guidance:} The recommended order is \strong{expression
#' refinement first, then chromatin refinement.} When matched chromatin data
#' are available, follow expression refinement with
#' \code{\link{refine_loop_anchors_by_chromatin}} -- chromatin evidence
#' resolves eP/eG into definitive types (E, P, G, dual), while the
#' expression-derived columns (\code{Putative_Target_Genes},
#' \code{Promoter_Target_Genes}, \code{Retained_In_Functional_Network})
#' carry forward the filtered target set for downstream profiling.
#' The reverse order (chromatin ->expression) is unsupported and will stop.
#' and will produce incorrect results. See
#' \code{\link{refine_loop_anchors_by_chromatin}} for the full
#' pipeline recommendation.
#'
#' @param annotation_res List. The raw foundational output object returned by \code{\link{annotate_peaks_and_loops}}.
#' @param expr_matrix_file Path to a normalised expression matrix (genes x samples). Accepts steady-state RNA-seq (TPM/FPKM/RPKM), nascent transcription data (NET-seq, PRO-seq, GRO-seq, TT-seq), or CAGE-seq. Required for refinement. Default: \code{NULL}.
#' @param sample_columns Character vector or integer indices. Columns in \code{expr_matrix_file} to average. Default: \code{NULL}.
#' @param threshold Numeric. Minimum expression value (e.g., TPM >= 1) for a gene to be considered active. Default: \code{1}. Set to \code{0} to retain any detected expression. Use \code{threshold_mode = "quantile"} to specify a quantile instead of an absolute value. For nascent transcription data (NET-seq, PRO-seq), use gene-body aggregated signal (e.g., RPKM) and adjust the threshold accordingly -- Pol II elongation signal is typically sparser than steady-state mRNA, so a lower absolute threshold or quantile mode is recommended.
#'   Gene name matching is \strong{case-insensitive}. Gene identifiers are
#'   canonicalised to uppercase before matching, so \code{"TP53"}, \code{"Tp53"},
#'   and \code{"tp53"} are treated as the same gene. This accommodates the
#'   common human (all-uppercase) vs mouse (Title Case) convention mismatch.
#' @param threshold_mode Character. How to interpret the \code{threshold}
#'   value. \code{"absolute"} (default): \code{threshold} is a direct expression
#'   cutoff (e.g., TPM >= 1). \code{"quantile"}: \code{threshold} is a quantile
#'   of the expression distribution (e.g., \code{0.75} means the top 25\% most
#'   highly expressed genes are considered active). Quantile mode is
#'   dataset-adaptive and may be more robust across experiments with different
#'   sequencing depths. The effective expression cutoff is reported in the log.
#' @param unit_type Character. Expression unit for plot labels (e.g., \code{"TPM"}, \code{"FPKM"}, \code{"RPKM"}, \code{"NETseq_RPKM"}, \code{"raw_count"}). Default: \code{"TPM"}.
#' @param species Character. Genome assembly string. Default: \code{"hg38"}.
#'   For non-human/non-mouse species, pass \code{annotation_res} produced by
#'   \code{\link{annotate_peaks_and_loops}} with custom \code{txdb}/\code{org_db}.
#' @param out_dir Character. Output directory for the Excel results file. Default: \code{"./results/filtered"}.
#' @param project_name Character. Prefix for output files (automatically appends \code{"_Filtered"}). Default: \code{"HiChIP"}.
#' @param color_palette Character. RColorBrewer palette name for loop-type colour assignments. Default: \code{"Paired"}.
#' @param karyo_bin_size Integer. Bin width in base pairs (bp) for karyotype heatmaps. Default: \code{1e5} (100 kb). Typical range: 5e4-5e5 depending on genome size and resolution.
#' @param allow_rerefine Logical. If \code{FALSE} (default), stops when
#'   the input has already been expression-refined. Set to \code{TRUE}
#'   to override (for interactive exploration only).
#' @param reclassify_by_expression Logical. If \code{TRUE} (default), silent
#'   promoters (P) and gene bodies (G) are reclassified as eP/eG while
#'   positional gene attribution is retained. Expression filtering is applied
#'   to current target-gene columns, not to anchor_gene. If \code{FALSE}, both anchor types and
#'   gene symbols are preserved unchanged -- genes absent from the expression
#'   matrix are kept (their expression state is unmeasured, not silent).
#'   Measured-silent genes remain in the gene list but are tracked separately
#'   via the \code{Target_Expression_State} and \code{Passes_Expression_Filter}
#'   columns in the output.
#' @param hub_percentile Numeric (0-1). Node-degree quantile for hub detection. Default: \code{0.95}.
#' @param hub_metric Character. Which connectivity count to use for hub detection. \code{"unique_contacts"} (default): counts distinct neighbour anchor IDs. \code{"total_loops"}: counts all loop rows (backward compatible).
#' @param chromatin_beds Named list of BED file paths for orthogonal chromatin
#'   mark validation of eP/eG anchors. When non-empty, a \emph{Chromatin
#'   Validation} sheet is added to the Excel workbook (see
#'   \code{\link{validate_epeG_by_chromatin}}). Default: \code{list()} (skip).
#'   Note: if you plan to use \code{\link{refine_loop_anchors_by_chromatin}},
#'   the chromatin validation and reclassification are handled there; you can
#'   skip this parameter to avoid redundant BED file reads.
#' @param write_output Logical. If \code{TRUE}, write the refined Excel workbook to \code{out_dir} (default: \code{FALSE}). If \code{FALSE}, return results without creating directories or files.
#' @param quiet Logical. If \code{TRUE}, suppress progress messages while preserving warnings. Default: \code{FALSE}.
#'
#' @return A named list:
#' \itemize{
#'   \item \code{loop_annotation} -- Full refined 3D network with updated \code{loop_type}
#'     (e.g., eP-P) and two target gene columns:
#'     \itemize{
#'       \item \code{Putative_Target_Genes}: \strong{Primary target-assignment
#'         result.} Candidate target genes based on final anchor type, TSS
#'         assignment, loop topology, and the current refinement stage.
#'         In expression-refined objects, all listed genes pass the
#'         expression threshold. \strong{Use for:} enhancer--gene linking
#'         benchmarks, CRISPRi validation, nearest-gene comparisons, loop
#'         connectivity analysis, full network export.
#'       \item \code{Promoter_Target_Genes}: \strong{Promoter-only subset
#'         of \code{Putative_Target_Genes}.} Contains only target genes
#'         supported by promoter-role (P/eP/dual) anchors. Always satisfies
#'         \code{Promoter_Target_Genes subset of Putative_Target_Genes}.
#'         \strong{Use for:} promoter-centric GO, KEGG, GSEA, active
#'         regulatory network from promoter-anchored loops.
#'       \item \code{Target_Expression_State}: Per-loop summary of target
#'         gene expression status (\code{"active"}, \code{"active_partial"},
#'         \code{"measured_silent"}, \code{"unmeasured"}, \code{"mixed"},
#'         \code{"no_target"}, \code{"not_assessed"}).
#'         \strong{\code{"unmeasured"} != \code{"measured_silent"}:}
#'         unmeasured means the gene was not in the expression input and
#'         its activity is unknown; measured-silent means it was measured
#'         but fell below the activity threshold.
#'     }
#'     Refinement status columns:
#'     \code{Has_Active_Target}, \code{Retained_In_Functional_Network}, and
#'     \code{Refinement_Action} (\code{"retained_active_target"},
#'     \code{"reclassified_silent_anchor"}, \code{"expression_filtered_no_active_target"},
#'     or \code{"structural_only_no_active_target"}).
#'     All structural loop rows are preserved (anchor types/genes are
#'     reclassified in place for silent elements); filter on
#'     \code{Retained_In_Functional_Network} for the high-confidence active subset.
#'   \item \code{anchor_loci_annotation} -- Filtered non-redundant anchor-locus annotations with expressed targets.
#'   \item \code{anchor_annotation} -- Backward-compatible alias of \code{anchor_loci_annotation}.
#'   \item \code{promoter_centric_stats} -- Gene-level connectivity statistics.
#'   \item \code{distal_element_stats} -- Distal anchor connectivity (E, dual, G, eG).
#'   \item \code{target_annotation} -- Target features (peaks) with gene assignments.
#'     Key columns include:
#'     \itemize{
#'       \item \code{All_Loop_Connected_Genes}: Inclusive provenance union of all loop-anchor gene links. May include strict assignment-eligible targets and non-strict positional/enhancer candidates. Not a confirmed target-gene set.
#'       \item \code{Regulated_promoter_genes}: Promoter genes supported by loop-anchor context.
#'       \item \code{Assigned_Target_Genes}: Policy-based 3D assignment (default: promoter-first, then shorter path wins; see \code{target_priority}).
#'       \item \code{*_Filled} variants: Linear nearest-gene fallback when strict 3D assignments are empty.
#'       \item \code{Regulated_promoter_Evidence}: Provenance of \code{Regulated_promoter_genes}
#'         (e.g., \code{local_promoter_overlap}, \code{distal_promoter}).
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
#'       \item \code{gene_role}: \code{"promoter"}, \code{"gene_body"},
#'       \code{"enhancer_candidate"}, \code{"positional_candidate"},
#'       or \code{"linear_annotation"}.
#'       \item \code{source}: \code{"loop_anchor"} (3D-derived) or \code{"linear_annotation"} (nearest gene).
#'       \item \code{evidence}: Provenance label --
#'         \code{"local_promoter_overlap"} (peak overlaps anchor promoter),
#'         \code{"distal_promoter"} (promoter on the distal loop anchor),
#'         \code{"gene_body_context"} / \code{"distal_gene_body_context"} (gene body linkage),
#'         \code{"local_enhancer_candidate"} / \code{"distal_enhancer_candidate"} /
#'         \code{"expanded_enhancer_candidate"} (enhancer-associated linkage),
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
#'       \item (Refine only) \code{Mean_Expression}: Per-gene mean expression value
#'         (\code{NA} when the gene is not in the expression matrix).
#'       \item (Refine only) \code{Expression_State}: Character. Per-gene expression
#'         status: \code{"active"} (measured and above threshold),
#'         \code{"measured_silent"} (measured but below threshold),
#'         \code{"unmeasured"} (gene not present in expression input),
#'         \code{"measured_not_assessed"} (measured but no threshold available).
#'         Loop-level aggregation fields also include \code{"active_partial"},
#'         \code{"mixed"}, \code{"no_target"}, and \code{"not_assessed"}.
#'       \item (Refine only) \code{Passes_Expression_Filter}: Logical with three
#'         states. \code{TRUE} = measured and active; \code{FALSE} = measured
#'         but below threshold; \code{NA} = unmeasured or not assessed.
#'         \strong{\code{NA} (unmeasured/not-assessed) \eqn{\neq}
#'         \code{FALSE} (measured-silent).}
#'     }
#'   \item \code{chromatin_validation} -- Data frame from
#'     \code{\link{validate_epeG_by_chromatin}} with enhancer evidence levels and
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
#'
#' @seealso \code{\link{annotate_peaks_and_loops}} for initial 3D annotation,
#'   \code{\link{refine_loop_anchors_by_chromatin}} for chromatin-aware reclassification.
#'
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
#'   annotation_res = raw_annotation,
#'   expr_matrix_file = expr_path,
#'   sample_columns = "con1",
#'   threshold = 1.0,
#'   unit_type = "TPM",
#'   species = "hg38",
#'   out_dir = tempdir(),
#'   project_name = "Example_Reclassified",
#'   reclassify_by_expression = TRUE,
#'   write_output = FALSE,
#'   quiet = TRUE
#' )
#' print(table(res_reclassified$loop_annotation$loop_type))
refine_loop_anchors_by_expression <- function(
  annotation_res,
  expr_matrix_file = NULL,
  sample_columns = NULL,
  threshold = 1,
  threshold_mode = c("absolute", "quantile"),
  unit_type = "TPM",
  species = "hg38",
  out_dir = "./results/filtered",
  project_name = "HiChIP",
  color_palette = "Paired",
  karyo_bin_size = 1e5,
  reclassify_by_expression = TRUE,
  hub_percentile = 0.95,
  hub_metric = c("unique_contacts", "total_loops"),
  chromatin_beds = list(),
  write_output = FALSE,
  quiet = FALSE,
  allow_rerefine = FALSE
) {
  threshold_mode <- match.arg(threshold_mode)
  hub_metric <- match.arg(hub_metric)
  .assert_nonempty_string(species, "species")
  log_message <- .make_log_message(quiet)

  # --- 0.  Empty annotation guard ---
  if (.is_empty_annotation_result(annotation_res)) {
    stop("Cannot refine an empty annotation result: no valid loop anchors are available.",
      call. = FALSE
    )
  }

  # --- 0-1. Setup & validate ---
  expr_validated <- .expression_validate_and_load(
    annotation_res, expr_matrix_file, sample_columns,
    threshold, threshold_mode, unit_type, reclassify_by_expression,
    allow_rerefine, project_name, out_dir, write_output, quiet
  )
  up_meta <- expr_validated$up_meta
  loop_df <- expr_validated$loop_df
  original_loop_df <- expr_validated$original_loop_df
  clust_info <- expr_validated$clust_info
  bed_info <- expr_validated$bed_info
  target_gene_links <- expr_validated$target_gene_links
  whitelist <- expr_validated$whitelist
  vals <- expr_validated$vals
  measured_set <- expr_validated$measured_set
  effective_threshold <- expr_validated$effective_threshold
  has_cluster_id <- expr_validated$has_cluster_id
  upstream_promoter_stats <- expr_validated$upstream_promoter_stats
  project_name <- expr_validated$project_name

  # --- 2. Update Anchors & Loops ---
  loop_df <- .refine_apply_anchor_updates(
    loop_df, whitelist, reclassify_by_expression, log_message,
    annotation_res = annotation_res, measured_set = measured_set
  )

  # --- 3. Stats Update ---
  log_message(">>> [Step 3] Updating Stats...")
  if (has_cluster_id) {
    agg_cluster <- loop_df %>%
      dplyr::filter(!is.na(cluster_id)) %>%
      dplyr::group_by(cluster_id) %>%
      dplyr::summarise(
        Cluster_All_Genes = extract_genes(Putative_Target_Genes),
        .groups = "drop"
      )
    loop_df <- loop_df %>%
      dplyr::select(-any_of(c("Cluster_All_Genes"))) %>%
      dplyr::left_join(agg_cluster, by = "cluster_id")
    if (!is.null(clust_info) && "cluster_id" %in% colnames(clust_info)) {
      clust_info <- clust_info %>%
        dplyr::select(-any_of(c("Cluster_All_Genes"))) %>%
        dplyr::left_join(agg_cluster, by = "cluster_id")
    }
  }

  stats_res <- compute_refined_stats(
    loop_df = loop_df,
    upstream_promoter_stats = upstream_promoter_stats,
    vals = vals, threshold = effective_threshold,
    hub_percentile = hub_percentile,
    hub_metric = hub_metric,
    anchor_registry = {
      st <- attr(annotation_res, "looplook_anchor_state", exact = TRUE)
      if (!is.null(st)) st$gr_anchors else NULL
    }
  )
  promoter_centric_df <- stats_res$promoter_centric
  distal_element_df <- stats_res$distal_element

  # --- 4. Refine Target Annotations ---
  log_message(">>> [Step 4] Refining Target Annotations...")
  if (!is.null(bed_info)) {
    tgt_refined <- .refine_target_annotations(
      bed_info, loop_df, whitelist, target_gene_links, vals, effective_threshold
    )
    bed_info <- tgt_refined$bed_info
    target_gene_links <- tgt_refined$target_gene_links

    # Re-aggregate strict target columns from current links using the
    # unified policy-controlled finalizer.  This ensures a silent promoter
    # link does not block an active gene-body link from becoming the
    # Assigned_Target_Genes entry.
    up_meta_tp <- annotation_res$metadata$parameters$target_priority
    up_meta_mph <- annotation_res$metadata$parameters$max_primary_hop
    if (!is.null(target_gene_links)) {
      final <- .finalize_current_target_annotation(
        bed_info = bed_info,
        all_target_gene_links = target_gene_links,
        has_expression = TRUE,
        target_priority = if (!is.null(up_meta_tp)) up_meta_tp else "promoter_then_distance",
        max_primary_hop = if (!is.null(up_meta_mph)) up_meta_mph else 1L
      )
      bed_info <- final$target_annotation
      target_gene_links <- final$target_gene_links
      bed_info <- .relocate_target_annotation_columns(bed_info)
    }
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

  # Capture upstream resource provenance for chain-of-custody tracking
  up_meta_params <- annotation_res$metadata$parameters

  resolved_neighbor_hop <- {
    st <- attr(annotation_res, "looplook_anchor_state", exact = TRUE)
    if (!is.null(st$neighbor_hop)) st$neighbor_hop
    else if (!is.null(up_meta_params$neighbor_hop)) up_meta_params$neighbor_hop
    else 0L
  }

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
        expression_refined = TRUE,
        threshold = threshold,
        effective_threshold = effective_threshold,
        threshold_mode = threshold_mode,
        unit_type = unit_type,
        reclassify_by_expression = reclassify_by_expression,
        hub_percentile = hub_percentile,
        hub_metric = hub_metric,
        allow_rerefine = allow_rerefine,
        # Forward upstream conflict-resolution policy so that
        # downstream chromatin TSS re-annotation inherits the
        # same biotype-priority and expression-tiebreaker rules.
        conflict_strategy = if (!is.null(up_meta_params$conflict_strategy)) {
          up_meta_params$conflict_strategy
        } else {
          "biotype_first"
        },
        co_dominance_ratio = if (!is.null(up_meta_params$co_dominance_ratio)) {
          up_meta_params$co_dominance_ratio
        } else {
          0.1
        },
        min_expr = if (!is.null(up_meta_params$min_expr)) {
          up_meta_params$min_expr
        } else {
          0
        },
        tss_region = if (!is.null(up_meta_params$tss_region)) {
          as.numeric(up_meta_params$tss_region)
        } else {
          NULL
        },
        target_priority = if (!is.null(up_meta_params$target_priority)) {
          up_meta_params$target_priority
        } else {
          "promoter_then_distance"
        },
        max_primary_hop = if (!is.null(up_meta_params$max_primary_hop)) {
          up_meta_params$max_primary_hop
        } else {
          1L
        },
        neighbor_hop = resolved_neighbor_hop,
        species = species,
        measured_gene_count = length(names(vals)),
        active_gene_count = length(whitelist),
        resource_inherited = !is.null(
          attr(annotation_res, "looplook_anchor_state", exact = TRUE)
        ),
        resource_changed = FALSE,
        upstream_resource_txdb = if (!is.null(up_meta_params)) {
          up_meta_params$resource_txdb_pkg
        } else {
          NA_character_
        },
        upstream_resource_orgdb = if (!is.null(up_meta_params)) {
          up_meta_params$resource_orgdb_pkg
        } else {
          NA_character_
        }
      ),
      genome_build = species,
      score_semantics = "expression-refined; eP/eG = element-Promoter like / element-Genebody like (structurally defined regulatory elements; validate with ATAC/H3K27ac)",
      database_versions = {
        as <- attr(annotation_res, "looplook_anchor_state", exact = TRUE)
        .record_database_versions(species,
          txdb_obj = if (!is.null(as)) as$txdb_obj else NULL,
          org_db_pkg = if (!is.null(as)) as$org_db_pkg else NULL,
          resource_reused = !is.null(as)
        )
      }
    )
  )

  # Inherit and update anchor_state so expression -> chromatin ->
  # recompute_targets pipeline stays closed.
  anchor_state <- attr(annotation_res, "looplook_anchor_state", exact = TRUE)
  if (!is.null(anchor_state)) {
    attr(out, "looplook_anchor_state") <- .update_anchor_state_from_loop_df(
      anchor_state, loop_df
    )
    attr(out, "looplook_anchor_state")$neighbor_hop <- resolved_neighbor_hop
    # Store expression data so downstream chromatin -> recompute_targets
    # can recover expression provenance for genes newly assigned to
    # E->P anchors (genes that never appeared in pre-chromatin target links).
    attr(out, "looplook_anchor_state")$expr_vals <- vals
    attr(out, "looplook_anchor_state")$expr_threshold <- effective_threshold
    attr(out, "looplook_anchor_state")$measured_set <- measured_set
    # Update gene_assignment_policy with current expression map so that
    # downstream chromatin TSS re-annotation queries the same expression
    # data.  conflict_min_expr is intentionally NOT updated -- it controls
    # the multi-candidate gene conflict resolution policy and is distinct
    # from the expression refinement threshold (see vignette for details).
    if (!is.null(attr(out, "looplook_anchor_state")$gene_assignment_policy)) {
      attr(out, "looplook_anchor_state")$gene_assignment_policy$gene_expr_map <- vals
    }
    # Store expression_policy separately so that the active/silent
    # threshold and the gene-conflict threshold are never conflated.
    attr(out, "looplook_anchor_state")$expression_policy <- list(
      threshold           = threshold,
      effective_threshold = effective_threshold,
      threshold_mode      = threshold_mode,
      unit_type           = unit_type
    )
  }

  return(out)
}

#' Internal: Validate chromatin refinement inputs and resolve tss_region
#'
#' Handles parameter validation, allow_rerefine guard, tss_region inheritance
#' from upstream annotation, and chromatin_bw/bigWig path checks.
#' Returns resolved tss_region, annotation metadata, and anchor_state.

#' Internal: Shared chromatin overlap parameter validator
#'
#' Validates chromatin_beds structure, candidate_types, anchor_gap,
#' and anchor_min_overlap.  Used by both
#' \code{\link{refine_loop_anchors_by_chromatin}} and
#' \code{\link{validate_epeG_by_chromatin}}.
#'
#' @keywords internal
#' @noRd
.validate_chromatin_overlap_inputs <- function(chromatin_beds, candidate_types,
                                               anchor_gap, anchor_min_overlap) {
  if (!is.list(chromatin_beds) ||
    (length(chromatin_beds) > 0 && is.null(names(chromatin_beds)))) {
    stop("`chromatin_beds` must be a named list of BED file paths.",
      call. = FALSE
    )
  }
  if (anyNA(names(chromatin_beds)) || any(!nzchar(names(chromatin_beds)))) {
    stop("`chromatin_beds` names must not be NA or empty strings.",
      call. = FALSE
    )
  }
  valid_types <- c("P", "G", "E", "eP", "eG")
  if (!is.null(candidate_types)) {
    if (!is.character(candidate_types) || length(candidate_types) == 0L ||
      anyNA(candidate_types) || anyDuplicated(candidate_types)) {
      stop("`candidate_types` must be a character vector of unique, non-NA anchor types.",
        call. = FALSE
      )
    }
    unknown <- setdiff(candidate_types, valid_types)
    if (length(unknown) > 0) {
      stop("Unknown candidate_type(s): ", paste(sQuote(unknown), collapse = ", "),
        ". Valid types: ", paste(valid_types, collapse = ", "), ".",
        call. = FALSE
      )
    }
  }
  if (!is.numeric(anchor_gap) || length(anchor_gap) != 1L ||
    !is.finite(anchor_gap) || anchor_gap != floor(anchor_gap) ||
    anchor_gap < -1L) {
    stop("`anchor_gap` must be a finite integer >= -1.", call. = FALSE)
  }
  if (!is.numeric(anchor_min_overlap) || length(anchor_min_overlap) != 1L ||
    !is.finite(anchor_min_overlap) ||
    anchor_min_overlap != floor(anchor_min_overlap) ||
    anchor_min_overlap < 1L) {
    stop("`anchor_min_overlap` must be a finite positive integer.", call. = FALSE)
  }
}

#' @keywords internal
#' @noRd
.chromatin_validate_inputs <- function(
  chromatin_bw, bw_ratio_threshold,
  annotation_res, allow_rerefine, tss_region, chromatin_beds, project_name,
  quiet
) {
  log_message <- .make_log_message(quiet)

  if (!is.numeric(bw_ratio_threshold) || length(bw_ratio_threshold) != 1L ||
    is.na(bw_ratio_threshold) || !is.finite(bw_ratio_threshold) ||
    bw_ratio_threshold <= 0) {
    stop("`bw_ratio_threshold` must be a single finite positive number (e.g. 3)", call. = FALSE)
  }

  if (!is.null(chromatin_bw)) {
    if (!is.list(chromatin_bw) || is.null(names(chromatin_bw))) {
      stop("`chromatin_bw` must be a named list of bigWig file paths (e.g. list(H3K4me1='path.bw', H3K4me3='path.bw'))", call. = FALSE)
    }
    if (!all(c("H3K4me1", "H3K4me3") %in% names(chromatin_bw))) {
      stop("`chromatin_bw` must include 'H3K4me1' and 'H3K4me3'", call. = FALSE)
    }
    for (nm in names(chromatin_bw)) {
      if (!file.exists(chromatin_bw[[nm]])) {
        stop("bigWig file not found for '", nm, "': ", chromatin_bw[[nm]], call. = FALSE)
      }
    }
    # Deliberately not routed through .require_pkg(): the semantics here are
    # "warn and CONTINUE in BED-only mode", which matches no on_missing mode.
    if (!requireNamespace("rtracklayer", quietly = TRUE)) {
      warning("Package 'rtracklayer' is required for bigWig processing; falling back to BED-only mode.", call. = FALSE)
    }
  }

  up_meta <- annotation_res$metadata

  # Guard: repeated chromatin refinement is not idempotent.
  if (!is.null(up_meta) && isTRUE(up_meta$function_name == "refine_loop_anchors_by_chromatin")) {
    if (!isTRUE(allow_rerefine)) {
      stop("Chromatin refinement cannot be applied to an already chromatin-refined object. ",
        "Anchor types and TSS assignments have already been modified. ",
        "Restart from the annotation- or expression-refined result. ",
        "To override (exploration only), set allow_rerefine = TRUE.",
        call. = FALSE
      )
    }
  }

  anchor_state <- attr(annotation_res, "looplook_anchor_state", exact = TRUE)
  metadata_tss <- if (!is.null(up_meta) && !is.null(up_meta$parameters)) {
    up_meta$parameters$tss_region
  } else {
    NULL
  }
  state_tss <- if (!is.null(anchor_state)) anchor_state$tss_region else NULL
  upstream_tss <- if (!is.null(metadata_tss)) metadata_tss else state_tss
  if (!is.null(upstream_tss)) {
    if (is.null(tss_region)) {
      tss_region <- upstream_tss
      log_message(
        ">>> Using tss_region from upstream annotation: ",
        paste(tss_region, collapse = ", ")
      )
    } else if (!identical(as.numeric(tss_region), as.numeric(upstream_tss))) {
      warning("Chromatin tss_region (", paste(tss_region, collapse = ","),
        ") differs from upstream annotation (",
        paste(upstream_tss, collapse = ","),
        "). This may produce inconsistent TSS re-annotation results. ",
        "Consider using the same tss_region across all refinement steps.",
        call. = FALSE
      )
    }
  }
  if (is.null(tss_region)) tss_region <- c(-2000, 2000)
  tss_region <- .validate_tss_region(tss_region)

  log_message(">>> [Chromatin Refinement] Project: ", project_name)

  # Chromatin-only workflow: inform user expression fields will be NA.
  if (!is.null(up_meta) &&
    !isTRUE(up_meta$function_name == "refine_loop_anchors_by_expression")) {
    log_message(
      "Running chromatin-only refinement. ",
      "Expression-dependent target fields (Mean_Expression, ",
      "Passes_Expression_Filter, Expression_State) will be ",
      "reported as not assessed."
    )
  }

  if (!is.null(up_meta) &&
    isTRUE(up_meta$function_name == "refine_loop_anchors_by_chromatin")) {
    log_message(
      "Re-refining an already chromatin-refined object ",
      "(allow_rerefine = TRUE). ",
      "Results depend on calling history and are not comparable ",
      "to single-pass refinement. For reproducible comparisons, ",
      "restart from a common upstream object."
    )
  }
  if (length(chromatin_beds) == 0) {
    stop("chromatin_beds must contain at least H3K4me1 and H3K4me3 BED files.", call. = FALSE)
  }

  list(tss_region = tss_region, up_meta = up_meta, anchor_state = anchor_state)
}

#' @title Chromatin-Aware refinement of loop anchor classification
#'
#' @description
#' Refines loop anchor types (P, E, G, eP, eG) using chromatin mark evidence
#' (H3K4me1, H3K4me3, and optionally H3K27ac, ATAC, H3K27me3).  Anchors with
#' dual promoter-enhancer chromatin signatures are flagged as \code{"dual"}.
#' \strong{Supported multi-omic workflow:} When both expression and
#' chromatin data are used, the required order is:\preformatted{
#'   annotate_peaks_and_loops()
#'   ->refine_loop_anchors_by_expression()
#'   ->refine_loop_anchors_by_chromatin()
#' }
#' Chromatin-only refinement (without prior expression refinement) is
#' also supported.
#'
#' \strong{Unsupported:} Running \code{refine_loop_anchors_by_expression()}
#' on a chromatin-refined object will raise an error. The reverse order
#' (chromatin ->expression) is not supported because
#' \code{clean_anchor()} only handles P and G types and cannot correctly
#' process \code{dual} or other types set by chromatin refinement.
#'
#' @details
#' \strong{Reclassification rules (minimum input: H3K4me1 + H3K4me3):}
#' \itemize{
#'   \item Enhancer BED overlap (highest priority): overlapped anchors become
#'         \code{"E"} (H3K4me3 absent) or \code{"dual"} (H3K4me3 present).
#'   \item P + H3K4me1(+) and H3K4me3(+) -> \code{"dual"} (dual-signature) or
#'         \code{"P"} (when H3K4me1 is likely a promoter shoulder; this is
#'         resolved by bigWig ratio >= 3 when \code{chromatin_bw} is provided).
#'   \item eP/eG + canonical or strong + active/intermediate enhancer
#'         chromatin -> \code{"E"} (enhancer identity confirmed).
#'   \item P + H3K4me1(+) and H3K4me3(-) \emph{and} (H3K27ac(+) or ATAC(+))
#'         -> \code{"E"} (conservative: requires active-mark confirmation beyond
#'         H3K4me1 alone).
#'   \item eP/eG + canonical or strong + active/intermediate enhancer
#'         chromatin -> \code{"E"} (enhancer identity confirmed; matches the
#'         P->E rule to guarantee the same outcome whether expression or
#'         chromatin refinement runs first).
#'   \item eP/eG with H3K4me3(+) but without H3K4me1(+) -> \code{"P"} (revert
#'         to promoter). This includes H3K4me3(+)/H3K27me3(+) bivalent anchors.
#'   \item Triple-positive H3K4me1(+)/H3K4me3(+)/H3K27me3(+) eP/eG anchors
#'         are classified as \code{conflicting_marks} and retain their
#'         previous eP/eG type. Peak co-occurrence alone is not sufficient
#'         evidence for restoration to a functional promoter state.
#'   \item E + H3K4me1(+) and H3K4me3(+) -> resolved by bigWig ratio:
#'         me1-dominant -> \code{"dual"}, unresolved (no bigWig) -> keep \code{"E"},
#'         not-me1-dominant -> \code{"P"}.
#'   \item E + H3K4me3(+) without H3K4me1 -> \code{"P"} (unannotated promoter).
#'   \item G + H3K4me1(+) and H3K4me3(+) -> same three-way bigWig resolution
#'         as E above.
#'   \item G + H3K4me3(+) without H3K4me1 -> \code{"P"} (internal promoter).
#'   \item G + H3K4me1(+) and H3K4me3(-) \emph{and} (H3K27ac(+) or ATAC(+))
#'         -> \code{"E"} (conservative intronic enhancer).
#'   \item All other anchors: unchanged.
#' }
#'
#' \strong{Required inputs:} \code{H3K4me1} and \code{H3K4me3} are required.
#' Their BED files must exist and contain at least one peak (an empty
#' required-mark file cannot be interpreted as genome-wide negative evidence).
#' Optional marks (\code{H3K27ac}, \code{ATAC}, \code{H3K27me3}) that are not
#' provided are recorded as \code{NA} (not_tested), and rules that depend on
#' their absence are skipped conservatively.
#'
#' \strong{Chromatin state inference:} Each anchor is also assigned a
#' \code{chromatin_state} based on its mark combination (highest-priority
#' match first):
#' \itemize{
#'   \item \code{"conflicting_marks"}: H3K27me3+ coexisting with any
#'         active mark (H3K4me1+, H3K27ac+, H3K4me3+, or ATAC+)  --
#'         bivalent/poised/ambiguous chromatin; takes priority over
#'         enhancer-like and dual-like classifications.
#'   \item \code{"dual_like"}: H3K4me1+ H3K4me3+
#'   \item \code{"active_enhancer_like"}: H3K4me1+ H3K27ac+ ATAC+
#'   \item \code{"intermediate_enhancer_like"}: H3K4me1+ with H3K27ac+ or ATAC+
#'   \item \code{"other_enhancer_like"}: H3K4me1+ only, or H3K27ac+ only
#'   \item \code{"promoter_like"}: H3K4me3+
#'   \item \code{"repressed"}: H3K27me3+ (no active marks)
#'   \item \code{"uncertain"}: marks tested but none present
#'   \item \code{"no_data"}: no marks tested
#' }
#'
#' \strong{Candidate anchor selection:} When \code{candidate_types = NULL}
#' (default), all five anchor types \code{c("P","G","E","eP","eG")} are tested
#' regardless of input source.  Set \code{candidate_types} to a subset
#' (e.g. \code{c("eP", "eG")}) to restrict reclassification scope.
#'
#' After reclassification, \code{loop_type} is recomputed from the updated
#' anchor types.
#'
#' \strong{Bivalent / conflicting-mark domains:} Anchors with H3K27me3+
#' coexisting alongside active marks (H3K4me3+, H3K4me1+, or H3K27ac+)
#' receive \code{chromatin_state = "conflicting_marks"}. This is common in
#' stem-cell and developmental contexts, where H3K4me3+/H3K27me3+ bivalent
#' domains poise genes for activation while maintaining transcriptional
#' silence. The reclassification rules distinguish between:
#' bivalent double-positive anchors (H3K4me3+ / H3K27me3+, H3K4me1-),
#' which are restored to \code{P} while retaining the
#' \code{conflicting_marks} state; and triple-positive anchors
#' (H3K4me1+ / H3K4me3+ / H3K27me3+), which retain their previous
#' eP/eG type because the co-presence of three marks on bulk ChIP
#' does not resolve whether H3K4me1 and H3K4me3 come from the
#' same nucleosome, allele, or cell. To distinguish active from
#' poised promoters, run expression refinement before chromatin
#' refinement and use the retained expression state together with
#' \code{chromatin_state = "conflicting_marks"}. The \code{chromatin_state} column in
#' \code{chromatin_validation} remains available for independent filtering
#' of conflicting-mark anchors.
#'
#' \strong{Pipeline guidance:}
#' The two refinement modules serve distinct roles in a complete analysis:
#' \itemize{
#'   \item \strong{Expression-aware refinement} evaluates transcriptional activity
#'         (is the gene expressed?) and produces expression-filtered
#'         \code{Putative_Target_Genes} and \code{Promoter_Target_Genes}
#'         plus \code{Retained_In_Functional_Network} and
#'         \code{Refinement_Action} for downstream profiling.
#'   \item \strong{Chromatin-aware refinement} evaluates chromatin identity
#'         (what \emph{is} this anchor?) using direct histone mark evidence,
#'         correcting anchor types and reclassifying eP/eG into definitive
#'         categories (E, P, G, dual).
#' }
#' When you have \strong{both RNA-seq and ChIP-seq}, use
#' \code{refine_loop_anchors_by_expression} first, then
#' \code{refine_loop_anchors_by_chromatin} -- expression marks silent anchors
#' as hypotheses (eP/eG), chromatin confirms or corrects them into definitive
#' types, and the expression-filtered target gene columns
#' (\code{Putative_Target_Genes}, \code{Promoter_Target_Genes} etc.)
#' carry the current target set into downstream profiling.
#' When you have \strong{only ChIP-seq}, run chromatin refinement directly on
#' the output of \code{\link{annotate_peaks_and_loops}}.
#'
#' @param annotation_res List. Output from
#'   \code{\link{annotate_peaks_and_loops}} or
#'   \code{\link{refine_loop_anchors_by_expression}}.
#' @param chromatin_beds Named list of BED file paths. At minimum
#'   \code{H3K4me1} and \code{H3K4me3} are required for meaningful
#'   reclassification; \code{H3K27ac}, \code{ATAC}, and \code{H3K27me3}
#'   improve confidence but are optional.
#' @param anchor_gap Integer. Candidate search radius in bp between an anchor
#'   and a chromatin mark peak. Default \code{200}. Does \emph{not} bypass the
#'   physical-overlap requirement enforced by \code{anchor_min_overlap}.
#' @param anchor_min_overlap Integer. Minimum required physical overlap (bp)
#'   between an anchor and a chromatin mark or enhancer BED region. Default
#'   \code{100}. Proximity-only matches (within \code{anchor_gap} but without
#'   actual base-pair intersection) are excluded. Reduce to \code{1--10}
#'   for narrow features (TFBS, DHS, SNPs).
#' @param species Character. Genome assembly. Default \code{"hg38"}.
#' @param tss_region Numeric vector of length 2, or \code{NULL}. Promoter
#'   window around the anchor centre in bp, used for TSS re-annotation of
#'   E->P / G->P anchors.  Default: \code{NULL} (auto-inherits from upstream
#'   annotation metadata; falls back to \code{c(-2000, 2000)} if neither
#'   is available).  A warning is emitted if an explicit value differs from
#'   the upstream annotation.
#' @param out_dir Character. Output directory for Excel export. Default \code{"./results/chromatin"}.
#' @param project_name Character. Project prefix. Default \code{"HiChIP"}.
#' @param color_palette Character. RColorBrewer palette name for the
#'   rose chart. Dumbbell and mark-enrichment heatmap use fixed academic
#'   palettes. Default: \code{"Paired"}.
#' @param sankey_colors Named character vector or \code{NULL}. Override the
#'   default type-to-color mapping in the Sankey diagram. Names must be
#'   anchor types (\code{"P"}, \code{"E"}, \code{"G"}, \code{"eP"},
#'   \code{"eG"}, \code{"dual"}), values are hex colours. When \code{NULL}
#'   (default), the colourblind-safe Wong palette is used:
#'   \code{E = "#E69F00"} (orange), \code{dual = "#CC0000"} (red),
#'   \code{P = "#0072B2"} (blue), \code{eP = "#009E73"} (bluish-green),
#'   \code{G = "#CC79A7"} (reddish-purple), \code{eG = "#56B4E9"} (sky blue).
#' @param chromatin_bw Named list of bigWig file paths, or \code{NULL}.
#'   When provided, H3K4me1/H3K4me3 signal ratios are computed for
#'   dual-positive anchors to distinguish dual-signature elements
#'   (ratio >= 3) from promoter-proximal H3K4me1 shoulders (ratio < 3,
#'   reclassified to P).  Requires the \pkg{rtracklayer} package.
#'   The list must include named elements \code{"H3K4me1"} and
#'   \code{"H3K4me3"}.  Default: \code{NULL} (BED-only mode --
#'   H3K4me3-positive anchors lacking H3K4me1 overlap are reclassified as P;
#'   dual-positive (H3K4me1+/H3K4me3+) anchors retain their original type
#'   pending bigWig resolution).
#'   \strong{Strongly recommended} when ChIP-seq data are available:
#'   bigWig signals provide the quantitative H3K4me1/H3K4me3 ratio needed
#'   to distinguish dual-signature elements from promoter-proximal
#'   H3K4me1 shoulders.  Use the companion script
#'   \code{inst/scripts/diagnose_h3k4me_ratio.R} to explore the ratio
#'   distribution in your data and choose an appropriate threshold.
#'   \strong{Important:} The H3K4me1/H3K4me3 ratio is only interpretable
#'   when both bigWig tracks are generated with comparable normalization
#'   (e.g., both CPM, both RPGC), identical scaling, binning, and
#'   smoothing, from matched biological material. Do not mix tracks with
#'   different normalization schemes (e.g., one fold-enrichment and one
#'   raw coverage) or from different cell types or batches. A separate
#'   low-signal threshold is not applied because input peaks already
#'   pass a peak caller (e.g. MACS2) and represent statistically
#'   enriched regions; signal at called peaks is expected to be well
#'   above noise.
#' @param bw_ratio_threshold Numeric. Minimum H3K4me1/H3K4me3 ratio to
#'   classify a dual-positive anchor as true \code{"dual"}.  Anchors
#'   below this threshold are reclassified as \code{"P"} (promoter
#'   shoulder).  Default: \code{3} (H3K4me1 signal must be at least
#'   3x H3K4me3 to override promoter identity).  Only used when
#'   \code{chromatin_bw} is provided.
#' @param enhancer_bed Character or \code{NULL}. Path to a BED file of
#'   known enhancer regions (e.g., FANTOM5, ENCODE cCREs).  Anchors
#'   overlapping these regions receive high-confidence enhancer
#'   classification: \code{"E"} when H3K4me3 is absent, \code{"dual"}
#'   when H3K4me3 is also present.  This curated evidence takes priority
#'   over chromatin-mark-derived enhancer evidence levels.  Default: \code{NULL}.
#' @param candidate_types Character vector or \code{NULL}. Anchor types to
#'   validate and reclassify. \code{NULL} (default): uses all five types
#'   \code{c("P","G","E","eP","eG")} regardless of input source. Chromatin
#'   evidence can reclassify any anchor type (e.g. P->dual, E->P), so
#'   restricting to expression-silenced types alone would miss biologically
#'   important transitions.
#'   Set to a subset to limit reclassification scope.
#' @param hub_metric Character. Which connectivity count to use for hub
#'   detection. \code{"unique_contacts"} (default) or \code{"total_loops"}.
#' @param recompute_targets Logical. If \code{TRUE} (default), re-run target gene
#'   assignment using updated anchor types. Requires the input
#'   \code{annotation_res} to contain the \code{looplook_anchor_state}
#'   attribute (present when using \code{\link{annotate_peaks_and_loops}}
#'   output). Default \code{TRUE}.
#' @param write_output Logical. Write Excel workbook. Default \code{FALSE}.
#' @param quiet Logical. Suppress messages. Default \code{FALSE}.
#' @param allow_rerefine Logical. If \code{FALSE} (default), stops when
#'   the input has already undergone chromatin refinement. Repeated chromatin
#'   refinement is not guaranteed to be idempotent (anchor types and TSS
#'   assignments have already been modified). Set \code{TRUE} only for
#'   exploratory use and restart from a common upstream object for formal
#'   parameter comparisons.
#'
#' @return An invisible named list with updated \code{loop_annotation},
#'   \code{anchor_loci_annotation}, \code{anchor_annotation} (alias of
#'   \code{anchor_loci_annotation}), \code{promoter_centric_stats},
#'   \code{distal_element_stats}, \code{chromatin_validation}
#'   (includes \code{chromatin_state}, \code{positional_type} and
#'   \code{final_type} columns per anchor),
#'   \code{plots} (\code{Chromatin_Dumbbell}: anchor-type before/after
#'   comparison; \code{Chromatin_Sankey}: anchor-type reclassification flow;
#'   \code{Chromatin_MarkHeatmap}: percentage of anchors positive for each
#'   mark per reclassification group; \code{Chromatin_UpSet}: loop-type UpSet
#'   plot (dot matrix + log10 bar chart)),
#'   \code{plot_list} (alias of \code{plots}), \code{qc_summary},
#'   \code{metadata}, and when \code{recompute_targets = TRUE},
#'   chromatin-aware \code{target_annotation} and
#'   \code{target_gene_links}. Set \code{recompute_targets = FALSE} to
#'   preserve pre-chromatin target assignments.
#'
#' @seealso \code{\link{annotate_peaks_and_loops}} for initial 3D annotation,
#'   \code{\link{refine_loop_anchors_by_expression}} for expression-aware refinement,
#'   \code{\link{validate_epeG_by_chromatin}} for standalone chromatin validation.
#'
#' @export
#'
#' @examples
#' # 1. Get paths to the required example files in the package
#' rdata_path <- system.file("extdata", "analysis_results.RData", package = "looplook")
#' k4me1_path <- system.file("extdata", "example_h3k4me1_peaks.bed", package = "looplook")
#' k4me3_path <- system.file("extdata", "example_h3k4me3_peaks.bed", package = "looplook")
#'
#' # 2. Load pre-computed annotation result
#' temp_env <- new.env()
#' load(rdata_path, envir = temp_env)
#' raw_annotation <- temp_env[[ls(temp_env)[1]]]
#' raw_annotation$loop_annotation <- head(raw_annotation$loop_annotation, 10)
#' raw_annotation$target_annotation <- head(raw_annotation$target_annotation, 3)
#'
#' # 3. Run chromatin-aware refinement
#' res_chromatin <- refine_loop_anchors_by_chromatin(
#'   annotation_res = raw_annotation,
#'   chromatin_beds = list(
#'     H3K4me1 = k4me1_path,
#'     H3K4me3 = k4me3_path
#'   ),
#'   anchor_gap = 500,
#'   anchor_min_overlap = 100,
#'   species = "hg38",
#'   out_dir = tempdir(),
#'   project_name = "Example_Chromatin",
#'   write_output = FALSE,
#'   quiet = TRUE
#' )
#' table(res_chromatin$loop_annotation$loop_type)
refine_loop_anchors_by_chromatin <- function(
  annotation_res,
  chromatin_beds = list(),
  anchor_gap = 200,
  anchor_min_overlap = 100,
  species = "hg38",
  tss_region = NULL,
  out_dir = "./results/chromatin",
  project_name = "HiChIP",
  color_palette = "Paired",
  candidate_types = NULL,
  recompute_targets = TRUE,
  write_output = FALSE,
  quiet = FALSE,
  sankey_colors = NULL,
  chromatin_bw = NULL,
  bw_ratio_threshold = 3,
  enhancer_bed = NULL,
  hub_metric = c("unique_contacts", "total_loops"),
  allow_rerefine = FALSE
) {
  hub_metric <- match.arg(hub_metric)
  .assert_nonempty_string(species, "species")
  if (!grepl("_Chromatin$", project_name)) project_name <- paste0(project_name, "_Chromatin")
  log_message <- .make_log_message(quiet)
  .ensure_out_dir(write_output, out_dir)

  # --- 0. Empty annotation guard ---
  if (.is_empty_annotation_result(annotation_res)) {
    stop("Cannot refine an empty annotation result: no valid loop anchors are available.",
      call. = FALSE
    )
  }

  # --- 0. Validate inputs & inherit tss_region ---
  validated <- .chromatin_validate_inputs(
    chromatin_bw, bw_ratio_threshold, annotation_res,
    allow_rerefine, tss_region, chromatin_beds, project_name, quiet
  )
  tss_region <- validated$tss_region
  up_meta <- validated$up_meta
  anchor_state <- validated$anchor_state
  # Resolve neighbor_hop early: old objects may carry metadata but
  # lack the anchor_state field.  Restore it before target recomputation
  # so that expanded targets use the correct hop setting.
  if (is.null(anchor_state$neighbor_hop)) {
    anchor_state$neighbor_hop <- if (!is.null(up_meta$parameters$neighbor_hop)) {
      up_meta$parameters$neighbor_hop
    } else {
      0L
    }
    attr(annotation_res, "looplook_anchor_state") <- anchor_state
  }
  .validate_chromatin_overlap_inputs(
    chromatin_beds, candidate_types,
    anchor_gap, anchor_min_overlap
  )

  # --- 1. Load & validate ---
  loop_df <- annotation_res$loop_annotation
  if (is.null(loop_df)) stop("'annotation_res$loop_annotation' is missing.")
  known_marks <- c("H3K4me1", "H3K27ac", "ATAC", "H3K27me3", "H3K4me3")
  # Case-insensitive matching of user-provided mark names against canonical
  # names (e.g., "h3k4me1"  ->  "H3K4me1"). Unmatched names are silently dropped.
  mark_lookup <- setNames(known_marks, toupper(known_marks))
  bed_names <- names(chromatin_beds)
  matched_idx <- toupper(bed_names) %in% names(mark_lookup)
  if (any(!matched_idx)) {
    warning("Unrecognised chromatin_beds name(s) ignored: ",
      paste(bed_names[!matched_idx], collapse = ", "),
      ". Expected: ", paste(known_marks, collapse = ", "),
      call. = FALSE
    )
  }
  # Normalise list element names to canonical case so downstream
  # chromatin_beds[[mark]] lookups work with canonical mark names.
  names(chromatin_beds)[matched_idx] <- unname(mark_lookup[toupper(bed_names[matched_idx])])
  if (anyDuplicated(names(chromatin_beds))) {
    stop("chromatin_beds contains duplicate mark names after case-insensitive normalisation; ",
      "ensure each mark appears only once.",
      call. = FALSE
    )
  }
  provided_marks <- intersect(names(chromatin_beds), known_marks)
  if (!all(c("H3K4me1", "H3K4me3") %in% provided_marks)) {
    stop("chromatin_beds must include at least 'H3K4me1' and 'H3K4me3'.", call. = FALSE)
  }
  valid_types <- c("P", "G", "E", "eP", "eG")
  if (!is.null(candidate_types)) {
    if (!is.character(candidate_types) || length(candidate_types) == 0L ||
      any(is.na(candidate_types)) || anyDuplicated(candidate_types)) {
      stop("`candidate_types` must be a character vector of unique, non-NA anchor types.",
        call. = FALSE
      )
    }
    unknown <- setdiff(candidate_types, valid_types)
    if (length(unknown) > 0) {
      stop("Unknown candidate_type(s): ", paste(sQuote(unknown), collapse = ", "),
        ". Valid types: ", paste(valid_types, collapse = ", "), ".",
        call. = FALSE
      )
    }
  }

  # --- 2. Extract candidate anchors and overlap marks ---
  epeG_anchors <- .extract_epeG_anchors(
    annotation_res, log_message,
    candidate_types
  )
  if (inherits(epeG_anchors, "data.frame") && nrow(epeG_anchors) == 0) {
    warning("No P/G/eP/eG anchors found; returning input unchanged.", call. = FALSE)
    return(annotation_res)
  }

  overlap_res <- .overlap_chromatin_marks(
    epeG_anchors, chromatin_beds, provided_marks, known_marks,
    anchor_gap, anchor_min_overlap, log_message
  )
  mark_matrix <- overlap_res$mark_matrix
  valid_provided_marks <- overlap_res$valid_provided_marks

  validation <- .assign_chromatin_confidence(
    epeG_anchors, mark_matrix, valid_provided_marks, known_marks
  )

  # --- 2b. Compute bigWig ratio for dual-positive anchors (optional) ---
  bw_metrics <- NULL
  bw_ratio <- NULL
  if (!is.null(chromatin_bw) && requireNamespace("rtracklayer", quietly = TRUE)) {
    log_message(">>> Computing H3K4me1 / H3K4me3 ratio from bigWig ...")
    bw_metrics <- .compute_bw_ratio(
      .with_known_upstream_noise_suppressed(
        GenomicRanges::makeGRangesFromDataFrame(epeG_anchors, keep.extra.columns = TRUE)
      ),
      mark_matrix,
      chromatin_bw[["H3K4me1"]],
      chromatin_bw[["H3K4me3"]]
    )
    bw_ratio <- setNames(bw_metrics$log2_ratio, bw_metrics$anchor_id)
    log_message(sprintf("    Ratio available for %d dual-positive anchors", nrow(bw_metrics)))
  }

  # --- 2c. External enhancer BED (optional) ---
  enhancer_anchors <- character(0)
  if (!is.null(enhancer_bed)) {
    if (!file.exists(enhancer_bed)) {
      stop("enhancer_bed file not found: ", enhancer_bed, call. = FALSE)
    }
    log_message(">>> Overlapping with external enhancer BED: ", basename(enhancer_bed))
    enh_gr <- read_simple_bed(enhancer_bed, quiet = quiet)
    enh_gr <- .harmonize_seqlevels(
      enh_gr,
      .with_known_upstream_noise_suppressed(
        GenomicRanges::makeGRangesFromDataFrame(epeG_anchors, keep.extra.columns = TRUE)
      ), "enhancer BED"
    )
    gh <- GenomicRanges::findOverlaps(
      .with_known_upstream_noise_suppressed(
        GenomicRanges::makeGRangesFromDataFrame(epeG_anchors, keep.extra.columns = TRUE)
      ), enh_gr,
      maxgap = -1L, minoverlap = anchor_min_overlap
    )
    enhancer_anchors <- epeG_anchors$anchor_id[unique(S4Vectors::queryHits(gh))]
    log_message(sprintf(
      "    %d / %d anchors overlap known enhancers",
      length(enhancer_anchors), nrow(epeG_anchors)
    ))
  }

  # --- 3. Reclassify anchors based on chromatin evidence ---
  log_message(">>> Reclassifying anchors by chromatin evidence...")
  reclass_map <- .chromatin_reclassify(validation,
    bw_ratio = bw_ratio,
    enhancer_anchors = enhancer_anchors,
    bw_ratio_threshold = bw_ratio_threshold,
    provided_marks = valid_provided_marks
  )
  log_message(sprintf("    Reclassified %d / %d anchors", sum(reclass_map$changed), nrow(reclass_map)))

  # --- 4. Apply reclassification to loop_annotation ---
  loop_df <- .chromatin_update_loops(loop_df, reclass_map, validation)

  # --- 4b. Restore gene symbols for eP/eG/E/G ->P/G/dual ---
  # Reuse the TxDb/OrgDb saved from the initial annotation so that
  # custom transcript catalogues and non-model species are respected.
  # Fall back to species-based resolution only when no saved resource
  # exists (e.g., annotation was loaded from file without anchor_state).
  anchor_state <- attr(annotation_res, "looplook_anchor_state", exact = TRUE)
  # Warn if stored topology lacks graph and ego-list elements are numeric
  # (requires graph$name to resolve).  Character anchor IDs do not need it.
  if (!is.null(anchor_state$ego_list_loop) && is.null(anchor_state$g)) {
    needs_graph <- any(vapply(
      anchor_state$ego_list_loop,
      function(x) is.numeric(x) && !inherits(x, "igraph.vs"),
      logical(1)
    ))
    if (needs_graph) {
      warning(
        "Stored topology is present but the graph object is missing. ",
        "Neighbor-hop gene reconstruction may be unavailable for ",
        "legacy or imported objects.",
        call. = FALSE
      )
    }
  }
  if (!is.null(anchor_state$txdb_obj)) {
    txdb_obj <- anchor_state$txdb_obj
    org_db_pkg <- anchor_state$org_db_pkg # may be NULL -- geneId fallback handles it
    txdb_reused <- TRUE
    # Validate TxDb connection (may be stale after saveRDS/readRDS).
    # Attempt auto-recovery from the recorded dbfile path when possible.
    up_meta <- annotation_res$metadata
    txdb_dbfile <- if (!is.null(up_meta) && !is.null(up_meta$database_versions)) {
      up_meta$database_versions$txdb_dbfile
    } else {
      NULL
    }
    txdb_obj <- .validate_or_recover_txdb(txdb_obj, txdb_dbfile)
    if (!is.null(org_db_pkg)) {
      log_message("    Using TxDb/OrgDb from initial annotation.")
    } else {
      log_message("    Using TxDb from initial annotation (no OrgDb; geneId fallback).")
    }
  } else {
    txdb_reused <- FALSE
    tx_res <- tryCatch(
      .resolve_annotation_resource(NULL, "txdb", "TxDb", species),
      error = function(e) NULL
    )
    org_res <- tryCatch(
      .resolve_annotation_resource(NULL, "orgdb", "OrgDb", species),
      error = function(e) NULL
    )
    txdb_obj <- if (!is.null(tx_res)) tx_res$obj else NULL
    org_db_pkg <- if (!is.null(org_res)) org_res$pkg else NULL
  }
  # Resolve conflict-resolution parameters for TSS re-annotation.
  # Primary source: gene_assignment_policy stored in anchor_state during
  # basic annotation (and updated by expression refinement).
  # Fallback: metadata$parameters for objects saved without the policy.
  policy <- anchor_state$gene_assignment_policy
  up_params <- if (!is.null(up_meta) && !is.null(up_meta$parameters)) {
    up_meta$parameters
  } else {
    list()
  }
  # gene_expr_map: policy first, then anchor_state$expr_vals, then NULL
  expr_map <- if (!is.null(policy$gene_expr_map)) {
    policy$gene_expr_map
  } else {
    anchor_state$expr_vals
  }
  # conflict_min_expr is the gene-conflict resolution threshold (distinct
  # from expression refinement's effective_threshold).  Fall back to the
  # deprecated 'min_expr' field for objects saved before the split.
  chr_min_expr <- if (!is.null(policy$conflict_min_expr)) {
    policy$conflict_min_expr
  } else if (!is.null(policy$min_expr)) {
    policy$min_expr
  } else if (!is.null(up_params$min_expr)) {
    up_params$min_expr
  } else {
    0
  }
  chr_conflict_strategy <- if (!is.null(policy$conflict_strategy)) {
    policy$conflict_strategy
  } else if (!is.null(up_params$conflict_strategy)) {
    up_params$conflict_strategy
  } else {
    "biotype_first"
  }
  chr_co_dominance <- if (!is.null(policy$co_dominance_ratio)) {
    policy$co_dominance_ratio
  } else if (!is.null(up_params$co_dominance_ratio)) {
    up_params$co_dominance_ratio
  } else {
    0.1
  }
  chr_biotype_order <- if (!is.null(policy$biotype_order)) {
    policy$biotype_order
  } else {
    c("protein", "small_ncRNA", "antisense", "lncRNA", "pseudogene")
  }
  restore_res <- .chromatin_restore_genes(
    loop_df, reclass_map, annotation_res, validation,
    txdb_obj = txdb_obj,
    org_db_pkg = org_db_pkg,
    tss_region = tss_region,
    gene_expr_map = expr_map,
    min_expr = chr_min_expr,
    conflict_strategy = chr_conflict_strategy,
    co_dominance_ratio = chr_co_dominance,
    biotype_order = chr_biotype_order,
    quiet = exists("quiet", inherits = TRUE) && quiet
  )
  loop_df <- restore_res$loop_df
  tss_reannotation <- restore_res$tss_reannotation
  tss_provenance_out <- restore_res$tss_assignment_provenance

  # --- 5-8. Post-reclassification pipeline (recompute types, summaries, stats, targets, plots) ---
  # Resolve hub_percentile from upstream metadata; fall back to 0.95.
  resolved_hub_percentile <- if (!is.null(up_meta) &&
    !is.null(up_meta$parameters$hub_percentile)) {
    up_meta$parameters$hub_percentile
  } else {
    0.95
  }
  pipeline <- .chromatin_post_reclassify(
    loop_df = loop_df, validation = validation, reclass_map = reclass_map,
    annotation_res = annotation_res, anchor_state = anchor_state,
    tss_provenance_out = tss_provenance_out,
    color_palette = color_palette, sankey_colors = sankey_colors,
    recompute_targets = recompute_targets, hub_metric = hub_metric,
    resolved_hub_percentile = resolved_hub_percentile,
    project_name = project_name, log_message = log_message
  )
  loop_df <- pipeline$loop_df
  validation <- pipeline$validation
  promoter_centric_df <- pipeline$promoter_centric_df
  distal_element_df <- pipeline$distal_element_df
  n_promoter_like <- pipeline$n_promoter_like
  chromatin_plots <- pipeline$chromatin_plots
  ta <- pipeline$ta
  tgl <- pipeline$tgl
  anchor_state <- pipeline$anchor_state

  # --- 9. Build output ---
  return(.chromatin_assemble_output(
    loop_df = loop_df, validation = validation,
    bw_metrics = bw_metrics, tss_reannotation = tss_reannotation,
    tss_provenance_out = tss_provenance_out,
    clust_info_in = annotation_res$anchor_loci_annotation,
    promoter_centric_df = promoter_centric_df,
    distal_element_df = distal_element_df,
    ta = ta, tgl = tgl, chromatin_plots = chromatin_plots,
    final_anchor_state = anchor_state,
    annotation_res = annotation_res, reclass_map = reclass_map,
    n_promoter_like = n_promoter_like,
    up_meta = up_meta, valid_provided_marks = valid_provided_marks,
    anchor_gap = anchor_gap, anchor_min_overlap = anchor_min_overlap,
    species = species, candidate_types = candidate_types,
    tss_region = tss_region, bw_ratio_threshold = bw_ratio_threshold,
    chromatin_bw = chromatin_bw, enhancer_bed = enhancer_bed,
    recompute_targets = recompute_targets, hub_metric = hub_metric,
    resolved_hub_percentile = resolved_hub_percentile,
    txdb_obj = txdb_obj, org_db_pkg = org_db_pkg, txdb_reused = txdb_reused,
    out_dir = out_dir, project_name = project_name,
    write_output = write_output, log_message = log_message
  ))
}

# --- Internal helpers for refine_loop_anchors_by_chromatin ---

#' Internal: Export chromatin refinement results to Excel
#' @keywords internal
#' @noRd
.chromatin_export_excel <- function(
  loop_df, validation, ta, tgl,
  tss_reannotation, tss_provenance_out, out_dir, project_name
) {
  wb <- openxlsx::createWorkbook()
  .add_sheet(wb, "Loop Annotation", loop_df, drop_cols = c("a1_id", "a2_id"))
  .add_sheet(wb, "Chromatin Validation", validation)
  .add_sheet(wb, "Target Annotation", ta, require_rows = TRUE)
  .add_sheet(wb, "Target Gene Links", tgl, require_rows = TRUE)
  .add_sheet(wb, "TSS Reannotation", tss_reannotation, require_rows = TRUE)
  .add_sheet(wb, "TSS Assignment Provenance", tss_provenance_out, require_rows = TRUE)
  .save_workbook(wb, out_dir, project_name, "_Chromatin_Results.xlsx",
    "Failed to save Excel: "
  )
}

#' Internal: Run post-reclassification chromatin pipeline (sections 5-8)
#'
#' Recomputes loop types, synchronizes anchor state, rebuilds connectivity
#' statistics, generates plots, and optionally reconstructs target gene links.
#' Returns a named list of all computed pipeline outputs.
#'
#' @keywords internal
#' @noRd

.chromatin_post_reclassify <- function(
  loop_df, validation, reclass_map,
  annotation_res, anchor_state,
  tss_provenance_out,
  color_palette, sankey_colors,
  recompute_targets, hub_metric,
  resolved_hub_percentile,
  project_name, log_message
) {
  # --- 5. Recompute loop_type, gene summaries, and stats ---
  log_message("    Recomputing loop types and gene summaries...")
  loop_df <- .chromatin_recompute_loop_types(loop_df)
  # Sync anchor_state map_info from updated loop_df so that neighbor-hop
  # topology expansion uses post-chromatin types and genes, not pre-
  # chromatin stale state. Use passed-in state unless unavailable.
  if (is.null(anchor_state)) {
    anchor_state <- attr(annotation_res, "looplook_anchor_state", exact = TRUE)
  }
  if (!is.null(anchor_state)) {
    anchor_state <- .update_anchor_state_from_loop_df(anchor_state, loop_df)
  }
  # Update chromatin validation with post-restoration gene state.
  # MUST run AFTER anchor_state sync so gene_after_chromatin reflects
  # the TSS-restored final gene, not the pre-restoration stale state.
  if (!is.null(validation) && nrow(validation) > 0) {
    if (!is.null(anchor_state) && !is.null(anchor_state$map_info)) {
      final_lookup <- setNames(
        anchor_state$map_info$SYMBOL,
        anchor_state$map_info$anchor_id
      )
      validation$gene_before_chromatin <- validation$anchor_gene
      validation$gene_after_chromatin <- final_lookup[validation$anchor_id]
    }
    # Join TSS provenance from restoration result
    if (!is.null(tss_provenance_out) && nrow(tss_provenance_out) > 0L) {
      validation <- dplyr::left_join(
        validation,
        tss_provenance_out,
        by = "anchor_id"
      )
    }
    # Join dual_ratio_state separately from reclass_map
    if ("dual_ratio_state" %in% colnames(reclass_map)) {
      ratio_prov <- reclass_map %>%
        dplyr::select(anchor_id, dual_ratio_state) %>%
        dplyr::filter(anchor_id %in% validation$anchor_id)
      if (nrow(ratio_prov) > 0) {
        validation <- dplyr::left_join(validation, ratio_prov, by = "anchor_id")
      }
    }
  }
  # --- Build complete anchor-level effective_gene_role from all anchors ---
  # Uses map_info (all anchors) as base, then overlays reclass_map + TSS
  # provenance for anchors that participated in chromatin refinement.
  # This ensures anchors outside candidate_types still have valid roles.
  if (!is.null(anchor_state) && !is.null(anchor_state$map_info) &&
    "a1_id" %in% colnames(loop_df)) {
    mi <- anchor_state$map_info
    # Preserve any existing effective_gene_role from prior refinement runs
    existing_role <- if ("effective_gene_role" %in% colnames(mi)) {
      as.character(mi$effective_gene_role)
    } else {
      rep(NA_character_, nrow(mi))
    }
    role_map <- data.frame(
      anchor_id = mi$anchor_id,
      type_code = mi$type_code,
      effective_gene_role = dplyr::case_when(
        !is.na(existing_role) & existing_role != "" ~ existing_role,
        mi$type_code %in% c("P", "eP") ~ "promoter",
        mi$type_code %in% c("G", "eG") ~ "gene_body",
        TRUE ~ "other"
      ),
      strict_assignment_eligible = rep(NA, nrow(mi)),
      stringsAsFactors = FALSE
    )
    # Overlay chromatin refinement results when available
    if (!is.null(reclass_map) && nrow(reclass_map) > 0) {
      has_tss <- !is.null(tss_provenance_out) && nrow(tss_provenance_out) > 0
      over_df <- reclass_map[, c("anchor_id", "old_type", "new_type"), drop = FALSE]
      if (has_tss) {
        over_df <- merge(over_df,
          tss_provenance_out[, c("anchor_id", "TSS_supported"), drop = FALSE],
          by = "anchor_id", all.x = TRUE, sort = FALSE
        )
      } else {
        over_df$TSS_supported <- NA
      }
      # Unified origin-aware resolver: same rules for anchor-level and
      # link-level.  TSS validates E-origin gene identity but does not
      # gate chromatin-state classification.  G/eG host genes are
      # structurally supported and do not require TSS.
      # gene_idx maps over_df (reclassified) -> mi (all anchors)
      gene_idx <- match(over_df$anchor_id, mi$anchor_id)
      gene_symbol <- as.character(mi$SYMBOL[gene_idx])
      has_gene <- !is.na(gene_idx) & !is.na(gene_symbol) &
        nzchar(trimws(gene_symbol))
      resolved <- .resolve_chromatin_gene_role(
        old_type     = over_df$old_type,
        final_type   = over_df$new_type,
        tss_supported = over_df$TSS_supported,
        has_gene     = has_gene
      )
      over_df$over_role <- resolved$role
      over_df$over_strict <- resolved$strict
      # overlay_idx maps role_map (all anchors) -> over_df (reclassified)
      overlay_idx <- match(role_map$anchor_id, over_df$anchor_id)
      ovr <- over_df$over_role[overlay_idx]
      ovs <- over_df$over_strict[overlay_idx]
      role_map$effective_gene_role[!is.na(ovr)] <- ovr[!is.na(ovr)]
      role_map$strict_assignment_eligible[!is.na(ovs)] <- ovs[!is.na(ovs)]
    }
    # Fill strict_assignment_eligible for anchors NOT in reclass_map
    # (non-chromatin-evaluated anchors): P/eP/G/eG are strict, E is not.
    na_strict <- is.na(role_map$strict_assignment_eligible)
    if (any(na_strict)) {
      role_map$strict_assignment_eligible[na_strict] <-
        role_map$effective_gene_role[na_strict] %in%
        c("promoter", "gene_body")
    }
    # Write to anchor_state$map_info
    map_idx <- match(anchor_state$map_info$anchor_id, role_map$anchor_id)
    anchor_state$map_info$effective_gene_role <- role_map$effective_gene_role[map_idx]
    anchor_state$map_info$strict_assignment_eligible <- role_map$strict_assignment_eligible[map_idx]
    # Write to loop_df
    role_lookup <- setNames(role_map$effective_gene_role, role_map$anchor_id)
    eligible_lookup <- setNames(role_map$strict_assignment_eligible, role_map$anchor_id)
    loop_df$anchor1_gene_role <- unname(role_lookup[as.character(loop_df$a1_id)])
    loop_df$anchor2_gene_role <- unname(role_lookup[as.character(loop_df$a2_id)])
    loop_df$anchor1_strict_eligible <- unname(eligible_lookup[as.character(loop_df$a1_id)])
    loop_df$anchor2_strict_eligible <- unname(eligible_lookup[as.character(loop_df$a2_id)])
    # Candidate positional genes: anchors whose own-gene link is NOT
    # strict-eligible (E->E nearest gene, E->P/dual unresolved, etc.).
    # Uses strict_assignment_eligible instead of checking gene_role alone,
    # so that enhancer_candidate with strict=FALSE (E->E) is also captured.
    loop_df <- loop_df %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        Candidate_Positional_Genes = extract_genes(c(
          if (isTRUE(!is.na(anchor1_strict_eligible) && !anchor1_strict_eligible &&
            !is.na(anchor1_gene) && nzchar(anchor1_gene))) anchor1_gene else NA_character_,
          if (isTRUE(!is.na(anchor2_strict_eligible) && !anchor2_strict_eligible &&
            !is.na(anchor2_gene) && nzchar(anchor2_gene))) anchor2_gene else NA_character_
        ))
      ) %>%
      dplyr::ungroup()
  }

  loop_df <- .rebuild_loop_gene_summaries(
    loop_df,
    expr_vals      = anchor_state$expr_vals,
    expr_threshold = anchor_state$expr_threshold,
    ego_list_loop  = anchor_state$ego_list_loop,
    map_info       = anchor_state$map_info,
    graph          = anchor_state$g
  )

  # --- 6. Summary ---
  log_message("--- Chromatin Refinement Summary ---")
  tab <- table(validation$enhancer_evidence, useNA = "ifany")
  for (lvl in names(tab)) {
    log_message(sprintf("  %-16s: %d anchors", lvl, tab[lvl]))
  }
  n_promoter_like <- sum(validation$is_promoter_like, na.rm = TRUE)
  if (n_promoter_like > 0) {
    log_message(sprintf("  %-16s: %d anchors", "promoter-like", n_promoter_like))
  }
  # Add positional/final type and chromatin_state layers to validation output
  validation$positional_type <- validation$anchor_type # pre-reclass
  final_types <- setNames(reclass_map$new_type, reclass_map$anchor_id)
  validation$final_type <- final_types[validation$anchor_id]
  validation$final_type[is.na(validation$final_type)] <- validation$positional_type[is.na(validation$final_type)]
  chromatin_states <- setNames(reclass_map$chromatin_state, reclass_map$anchor_id)
  validation$chromatin_state <- chromatin_states[validation$anchor_id]
  validation$chromatin_state[is.na(validation$chromatin_state)] <- "no_data"

  log_message(sprintf("  Reclassified      : %d anchors", sum(reclass_map$changed)))
  if (!is.null(tss_provenance_out) && nrow(tss_provenance_out) > 0) {
    n_tss_supported <- sum(tss_provenance_out$TSS_supported, na.rm = TRUE)
    n_tss_provisional <- sum(tss_provenance_out$Gene_Assignment_Confidence == "provisional", na.rm = TRUE)
    n_tss_unresolved <- sum(tss_provenance_out$Gene_Assignment_Confidence == "unresolved", na.rm = TRUE)
    n_tss_unavailable <- sum(tss_provenance_out$TSS_support_status == "TSS_validation_unavailable", na.rm = TRUE)
    log_message(sprintf("  TSS supported     : %d anchors", n_tss_supported))
    log_message(sprintf("  TSS provisional   : %d anchors", n_tss_provisional))
    log_message(sprintf("  TSS unresolved    : %d anchors", n_tss_unresolved))
    if (n_tss_unavailable > 0) {
      log_message(sprintf("  TSS not available : %d anchors", n_tss_unavailable))
    }
  }
  # --- Role transition QC: verify that chromatin reclassification is
  #     actually propagating to target gene roles.  If links are reclassified
  #     but no target summary rows change, warn loudly. ---
  n_target_reclassified <- sum(
    reclass_map$changed &
      reclass_map$old_type %in% c("P", "eP", "G", "eG", "E") &
      reclass_map$new_type %in% c("P", "eP", "dual", "G", "eG", "E"),
    na.rm = TRUE
  )
  # Count role transitions on reclassified anchors
  n_g_to_p <- sum(reclass_map$old_type %in% c("G", "eG") &
    reclass_map$new_type %in% c("P", "dual"), na.rm = TRUE)
  n_e_to_p <- sum(reclass_map$old_type == "E" &
    reclass_map$new_type %in% c("P", "dual"), na.rm = TRUE)
  n_p_to_e <- sum(reclass_map$old_type %in% c("P", "eP") &
    reclass_map$new_type == "E", na.rm = TRUE)
  n_g_to_e <- sum(reclass_map$old_type %in% c("G", "eG") &
    reclass_map$new_type == "E", na.rm = TRUE)

  log_message("--- Role Transition QC ---")
  log_message(sprintf("  G/eG -> P/dual   : %d anchors", n_g_to_p))
  log_message(sprintf("  E -> P/dual      : %d anchors", n_e_to_p))
  log_message(sprintf("  P/eP -> E        : %d anchors", n_p_to_e))
  log_message(sprintf("  G/eG -> E        : %d anchors", n_g_to_e))
  log_message(sprintf("  reclassified total : %d anchors", n_target_reclassified))

  if (n_target_reclassified > 0 && n_g_to_p + n_e_to_p == 0 && n_p_to_e == 0) {
    log_message("  (all role transitions are within same category -- target impact expected to be minimal)")
  }

  log_message("--- End Chromatin Refinement ---")

  # --- 6b. Visualization ---
  log_message("    Generating plots...")
  chromatin_plots <- list(
    Chromatin_Dumbbell = .build_chromatin_dumbbell_plot(reclass_map, project_name),
    Chromatin_Sankey = .build_chromatin_sankey_plot(reclass_map, project_name, sankey_colors = sankey_colors),
    Chromatin_MarkHeatmap = .build_chromatin_mark_heatmap(validation, reclass_map, project_name),
    Chromatin_UpSet = .build_loop_type_upset(loop_df,
      paste0(project_name, ": Loop Types After Chromatin Refinement"),
      subtitle = "Loop type distribution after chromatin-aware anchor refinement."
    )
  )

  # --- 7. Recompute stats after reclassification ---
  log_message("    Recomputing connectivity stats...")
  new_promoter_stats <- .compute_raw_promoter_stats(loop_df)
  has_expr <- !is.null(anchor_state$expr_vals) &&
    !is.null(anchor_state$expr_threshold)
  if (has_expr) {
    use_vals <- anchor_state$expr_vals
    use_thr <- anchor_state$expr_threshold
  } else {
    # No assessable expression data: all genes get Not_assessed.
    # Downstream override replaces Is_Active_Gene from upstream stats
    # where available, or sets Not_assessed for new genes.
    use_vals <- setNames(
      rep(NA_real_, length(unique(new_promoter_stats$Gene))),
      unique(new_promoter_stats$Gene)
    )
    use_thr <- 0
  }
  promoter_centric_df <- .build_promoter_centric_df(
    new_promoter_stats, annotation_res$promoter_centric_stats,
    vals = use_vals, threshold = use_thr,
    hub_percentile = resolved_hub_percentile,
    hub_metric = hub_metric
  )
  # When no expression data is available, all genes are Not_assessed.
  # If upstream stats already have Is_Active_Gene, inherit those values.
  if (!has_expr) {
    promoter_centric_df$Is_Active_Gene <- "Not_assessed"
    upstream_stats <- annotation_res$promoter_centric_stats
    if (!is.null(upstream_stats) && all(c("Gene", "Is_Active_Gene") %in% colnames(upstream_stats))) {
      upstream_active <- setNames(
        as.character(upstream_stats$Is_Active_Gene),
        as.character(upstream_stats$Gene)
      )
      matched <- promoter_centric_df$Gene %in% names(upstream_active)
      promoter_centric_df$Is_Active_Gene[matched] <-
        unname(upstream_active[promoter_centric_df$Gene[matched]])
    }
  }
  distal_element_df <- tryCatch(
    .build_distal_element_df(loop_df, resolved_hub_percentile,
      hub_metric = hub_metric,
      anchor_registry = anchor_state$gr_anchors
    ),
    error = function(e) {
      warning(
        "Failed to rebuild distal-element statistics: ",
        conditionMessage(e),
        call. = FALSE
      )
      NULL
    }
  )

  # --- 8. Recompute target links if requested ---
  # Preserve original assignments when recompute_targets = FALSE
  ta <- annotation_res$target_annotation
  tgl <- annotation_res$target_gene_links
  if (isTRUE(recompute_targets)) {
    log_message(">>> Recomputing target gene links with chromatin-aware anchors...")
    # Use the already-synchronized anchor_state from earlier in this function
    if (!is.null(anchor_state)) {
      attr(annotation_res, "looplook_anchor_state") <- .update_anchor_state_from_loop_df(
        anchor_state, loop_df
      )
    }
    recomputed <- .recompute_targets_from_anchor_state(
      list(
        loop_annotation = loop_df,
        chromatin_validation = validation
      ), annotation_res, reclass_map,
      tss_assignment_provenance = tss_provenance_out
    )
    ta <- recomputed$target_annotation
    tgl <- .ensure_target_link_schema(recomputed$target_gene_links)
    log_message(sprintf("    Updated %d target gene links.", nrow(tgl)))
  }
  list(
    loop_df = loop_df, validation = validation,
    promoter_centric_df = promoter_centric_df,
    distal_element_df = distal_element_df,
    n_promoter_like = n_promoter_like,
    chromatin_plots = chromatin_plots,
    ta = ta, tgl = tgl,
    anchor_state = anchor_state
  )
}
.chromatin_assemble_output <- function(
  loop_df, validation,
  bw_metrics, tss_reannotation, tss_provenance_out,
  clust_info_in, promoter_centric_df, distal_element_df,
  ta, tgl, chromatin_plots, final_anchor_state,
  annotation_res, reclass_map, n_promoter_like,
  up_meta, valid_provided_marks,
  anchor_gap, anchor_min_overlap, species, candidate_types, tss_region,
  bw_ratio_threshold, chromatin_bw, enhancer_bed, recompute_targets, hub_metric,
  resolved_hub_percentile,
  txdb_obj, org_db_pkg, txdb_reused,
  out_dir, project_name, write_output, log_message
) {
  n_tss_supported <- if (!is.null(tss_provenance_out)) {
    sum(tss_provenance_out$TSS_supported, na.rm = TRUE)
  } else {
    0L
  }
  n_tss_provisional <- if (!is.null(tss_provenance_out)) {
    sum(tss_provenance_out$Gene_Assignment_Confidence == "provisional", na.rm = TRUE)
  } else {
    0L
  }
  n_tss_unresolved <- if (!is.null(tss_provenance_out)) {
    sum(tss_provenance_out$Gene_Assignment_Confidence == "unresolved", na.rm = TRUE)
  } else {
    0L
  }
  n_tss_unavailable <- if (!is.null(tss_provenance_out)) {
    sum(tss_provenance_out$TSS_support_status == "TSS_validation_unavailable", na.rm = TRUE)
  } else {
    0L
  }

  qc_summary <- data.frame(
    n_candidate_anchors = nrow(validation),
    n_reclassified = sum(reclass_map$changed, na.rm = TRUE),
    n_canonical = sum(validation$enhancer_evidence == "canonical", na.rm = TRUE),
    n_strong = sum(validation$enhancer_evidence == "strong", na.rm = TRUE),
    n_supported = sum(validation$enhancer_evidence == "supported", na.rm = TRUE),
    n_limited = sum(validation$enhancer_evidence == "limited", na.rm = TRUE),
    n_dual = sum(reclass_map$new_type == "dual", na.rm = TRUE),
    n_promoter_like = n_promoter_like,
    n_TSS_supported = n_tss_supported,
    n_TSS_provisional = n_tss_provisional,
    n_TSS_unresolved = n_tss_unresolved,
    n_TSS_validation_unavailable = n_tss_unavailable,
    stringsAsFactors = FALSE
  )

  clust_info <- clust_info_in
  if (!is.null(clust_info) && nrow(clust_info) > 0 &&
    "cluster_id" %in% colnames(clust_info) &&
    "cluster_id" %in% colnames(loop_df)) {
    loop_cluster_genes <- loop_df %>%
      dplyr::filter(!is.na(cluster_id)) %>%
      dplyr::group_by(cluster_id) %>%
      dplyr::summarise(
        Cluster_Locus_Genes = extract_genes(All_Anchor_Genes),
        Cluster_All_Genes = extract_genes(Putative_Target_Genes),
        .groups = "drop"
      )
    clust_info <- clust_info %>%
      dplyr::select(-dplyr::any_of(c(
        "Cluster_Locus_Genes",
        "Cluster_All_Genes"
      ))) %>%
      dplyr::left_join(loop_cluster_genes, by = "cluster_id")
  }

  resolved_neighbor_hop <- {
    if (!is.null(final_anchor_state$neighbor_hop)) final_anchor_state$neighbor_hop
    else if (!is.null(up_meta$parameters$neighbor_hop)) up_meta$parameters$neighbor_hop
    else 0L
  }
  final_anchor_state$neighbor_hop <- resolved_neighbor_hop

  out <- list(
    loop_annotation = loop_df,
    chromatin_validation = validation,
    qc_summary = qc_summary,
    bigWig_metrics = bw_metrics,
    tss_reannotation = tss_reannotation,
    tss_assignment_provenance = tss_provenance_out,
    anchor_loci_annotation = clust_info,
    anchor_annotation = clust_info,
    promoter_centric_stats = promoter_centric_df,
    distal_element_stats = distal_element_df,
    target_annotation = ta,
    target_gene_links = tgl,
    plots = chromatin_plots,
    plot_list = chromatin_plots,
    metadata = .build_looplook_metadata(
      fun = "refine_loop_anchors_by_chromatin",
      params = list(
        chromatin_refined = TRUE,
        refinement_order = if (!is.null(up_meta) &&
          identical(up_meta$function_name, "refine_loop_anchors_by_expression")) {
          "expression_then_chromatin"
        } else {
          "chromatin_only"
        },
        expression_refined = !is.null(up_meta) &&
          identical(up_meta$function_name, "refine_loop_anchors_by_expression"),
        anchor_gap = anchor_gap,
        anchor_min_overlap = anchor_min_overlap,
        species = species,
        candidate_types = if (is.null(candidate_types)) "derived" else candidate_types,
        tss_region = tss_region,
        bw_ratio_threshold = bw_ratio_threshold,
        bw_pseudocount = if (!is.null(bw_metrics) && nrow(bw_metrics) > 0L) bw_metrics$pseudocount[[1]] else NA_real_,
        chromatin_bw_requested = !is.null(chromatin_bw),
        chromatin_bw_used = !is.null(bw_metrics) && nrow(bw_metrics) > 0L,
        has_enhancer_bed = !is.null(enhancer_bed),
        recompute_targets = recompute_targets,
        hub_percentile = resolved_hub_percentile,
        hub_metric = hub_metric,
        target_priority = if (!is.null(up_meta$parameters$target_priority)) {
          up_meta$parameters$target_priority
        } else {
          "promoter_then_distance"
        },
        max_primary_hop = if (!is.null(up_meta$parameters$max_primary_hop)) {
          up_meta$parameters$max_primary_hop
        } else {
          1L
        },
        neighbor_hop = resolved_neighbor_hop,
        dual_ratio_definition = "anchor_wide_mean_log2_ratio",
        provided_chromatin_marks = sort(valid_provided_marks),
        tss_support_policy = "provisional_on_no_annotated_TSS",
        upstream_expression = if (!is.null(up_meta) &&
          identical(up_meta$function_name, "refine_loop_anchors_by_expression")) {
          up_meta$parameters
        } else {
          list(expression_refined = FALSE)
        },
        resource_inherited = txdb_reused,
        resource_changed = FALSE,
        upstream_resource_txdb = if (!is.null(up_meta) &&
          !is.null(up_meta$parameters)) {
          pkg <- up_meta$parameters$upstream_resource_txdb
          if (is.null(pkg)) pkg <- up_meta$parameters$resource_txdb_pkg
          if (is.null(pkg) && !is.null(up_meta$database_versions)) {
            pkg <- up_meta$database_versions$txdb_pkg
          }
          pkg
        } else {
          NA_character_
        },
        upstream_resource_orgdb = if (!is.null(up_meta) &&
          !is.null(up_meta$parameters)) {
          pkg <- up_meta$parameters$upstream_resource_orgdb
          if (is.null(pkg)) pkg <- up_meta$parameters$resource_orgdb_pkg
          if (is.null(pkg) && !is.null(up_meta$database_versions)) {
            pkg <- up_meta$database_versions$orgdb_pkg
          }
          pkg
        } else {
          NA_character_
        }
      ),
      genome_build = species,
      score_semantics = "chromatin-aware reclassification; dual = promoter-enhancer dual-signature chromatin state (not functional proof)",
      database_versions = .record_database_versions(species,
        txdb_obj = txdb_obj, org_db_pkg = org_db_pkg,
        resource_reused = txdb_reused
      )
    )
  )
  attr(out, "looplook_anchor_state") <- final_anchor_state
  if (write_output) {
    log_message(">>> Exporting to Excel...")
    .chromatin_export_excel(
      loop_df, validation, ta, tgl,
      tss_reannotation, tss_provenance_out, out_dir, project_name
    )
    log_message("    Excel saved.")
  }
  log_message("Chromatin Refinement Complete.")
  return(invisible(out))
}

#' Internal: Build reclassification map from chromatin validation
#' @keywords internal
#' @noRd
.chromatin_reclassify <- function(validation, bw_ratio = NULL,
                                  enhancer_anchors = character(0),
                                  bw_ratio_threshold = 3,
                                  provided_marks = character(0)) {
  # Three-state mark status:
  #   detected   = BED provided AND peak overlaps anchor
  #   not_called = BED provided AND no peak at anchor
  #   not_tested = BED not provided
  # Only not_called provides negative evidence for exclusionary rules.
  .is_detected <- function(x) !is.na(x) & x
  .is_not_called <- function(x, mark) mark %in% provided_marks & !is.na(x) & !x
  out <- data.frame(
    anchor_id        = validation$anchor_id,
    positional_type  = validation$anchor_type,
    old_type         = validation$anchor_type,
    new_type         = validation$anchor_type,
    chromatin_state  = NA_character_,
    changed          = FALSE,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(
      h3k4me1_p = .is_detected(validation$H3K4me1),
      h3k4me3_p = .is_detected(validation$H3K4me3),
      h3k27ac_p = .is_detected(validation$H3K27ac),
      h3k27me3_p = .is_detected(validation$H3K27me3),
      atac_p = .is_detected(validation$ATAC),
      # Explicitly tested and not called (requires BED to be provided)
      h3k4me3_not_called = .is_not_called(validation$H3K4me3, "H3K4me3"),
      h3k4me1_not_called = .is_not_called(validation$H3K4me1, "H3K4me1"),
      is_promoter_like = validation$is_promoter_like,
      conf_chr = as.character(validation$enhancer_evidence),

      # --- bigWig dual resolution: three-state classification ---
      # dual_ratio_state classifies each dual-positive anchor:
      #   "me1_dominant"     = ratio data available, H3K4me1 dominates -> dual
      #   "not_me1_dominant" = ratio data available, H3K4me3 dominates -> P/G
      #   "unresolved"       = ratio data unavailable -> keep old type
      # "unresolved" != "not_me1_dominant": missing data does not
      # imply promoter dominance.
      dual_ratio_state = if (!is.null(bw_ratio)) {
        log2r <- bw_ratio[as.character(anchor_id)]
        dplyr::case_when(
          is.na(log2r) ~ "unresolved",
          log2r >= log2(bw_ratio_threshold) ~ "me1_dominant",
          TRUE ~ "not_me1_dominant"
        )
      } else {
        "unresolved"
      },

      # External enhancer BED overlap -- curated knowledge, high confidence.
      # Use vectorised `&` (not scalar `&&`) so each anchor is evaluated
      # independently against the enhancer_anchors set.
      is_enhancer_bed = length(enhancer_anchors) > 0L &
        as.character(anchor_id) %in% enhancer_anchors,

      # --- chromatin_state inference (first match wins) ---
      chromatin_state = dplyr::case_when(
        h3k27me3_p & (h3k4me1_p | h3k27ac_p | h3k4me3_p | atac_p) ~ "conflicting_marks",
        h3k4me1_p & h3k4me3_p ~ "dual_like",
        h3k4me1_p & h3k27ac_p & atac_p ~ "active_enhancer_like",
        h3k4me1_p & (h3k27ac_p | atac_p) ~ "intermediate_enhancer_like",
        h3k4me1_p | h3k27ac_p ~ "other_enhancer_like",
        h3k4me3_p ~ "promoter_like",
        h3k27me3_p ~ "repressed",
        !is.na(validation$H3K4me1) | !is.na(validation$H3K4me3) |
          !is.na(validation$H3K27ac) ~ "uncertain",
        TRUE ~ "no_data"
      ),

      # --- anchor-type reclassification (first match wins) ---
      new_type = dplyr::case_when(
        # 1. Enhancer BED -- highest priority curated knowledge
        #    (FANTOM5, ENCODE cCREs).  Overrides all signal-level
        #    rules including conflicting_marks and dual_like.
        is_enhancer_bed & h3k4me3_p ~ "dual",
        is_enhancer_bed & h3k4me3_not_called ~ "E",
        # 2. eP/eG + dual_like: resolved by bigWig ratio when available.
        #    me1-dominant -> dual-signature element.
        #    unresolved (no bigWig) -> keep old type (conservative).
        #    not_me1_dominant (K4me3 dominates) -> promoter identity.
        old_type %in% c("eP", "eG") & chromatin_state == "dual_like" &
          dual_ratio_state == "me1_dominant" ~ "dual",
        old_type %in% c("eP", "eG") & chromatin_state == "dual_like" &
          dual_ratio_state == "unresolved" ~ old_type,
        old_type %in% c("eP", "eG") & chromatin_state == "dual_like" ~ "P",
        # 3. eP/eG with H3K4me3+ but NOT dual_like (i.e. H3K4me3
        #    without H3K4me1) -> promoter identity.  This catches
        #    promoter_like anchors and bivalent (H3K4me3+H3K27me3)
        #    anchors that lack H3K4me1.  H3K4me3 defines promoter
        #    identity; conflicting H3K27me3 indicates poised status
        #    but does not negate the structural promoter identity.
        #    Note: use the original validation columns directly to
        #    avoid dplyr NSE edge cases with ! on computed columns.
        old_type == "eP" & !is.na(validation$H3K4me3) & validation$H3K4me3 &
          (is.na(validation$H3K4me1) | !validation$H3K4me1) ~ "P",
        old_type == "eG" & !is.na(validation$H3K4me3) & validation$H3K4me3 &
          (is.na(validation$H3K4me1) | !validation$H3K4me1) ~ "P",
        # 4. eP/eG + canonical/strong active-enhancer chromatin -> E
        old_type == "eP" &
          conf_chr %in% c("canonical", "strong") &
          chromatin_state %in% c("active_enhancer_like", "intermediate_enhancer_like") ~ "E",
        old_type == "eG" &
          conf_chr %in% c("canonical", "strong") &
          chromatin_state %in% c("active_enhancer_like", "intermediate_enhancer_like") ~ "E",
        # 5. eP/eG + promoter_like (H3K4me3+ only, no dual_like or
        #    active-enhancer caught above) -> P.  This is the general
        #    rule: any eP/eG with K4me3 should be a promoter.
        #    (Note: rules 2-3 already handle the specific subcases.)
        old_type == "eP" & chromatin_state == "promoter_like" ~ "P",
        old_type == "eG" & chromatin_state == "promoter_like" ~ "P",
        # 6. Conflicting active + repressive marks (H3K27me3+ with
        #    H3K4me1/me3/H3K27ac/ATAC) -> bivalent/poised chromatin.
        #    Keep original type for non-eP/eG anchors; do not output
        #    a high-confidence P/E/dual classification.
        chromatin_state == "conflicting_marks" ~ old_type,
        # 7. Dual-positive P/G/E with bigWig resolution
        old_type == "P" & h3k4me1_p & h3k4me3_p &
          dual_ratio_state == "me1_dominant" ~ "dual",
        old_type == "P" & h3k4me1_p & h3k4me3_p &
          dual_ratio_state == "unresolved" ~ "P",
        old_type == "P" & h3k4me1_p & h3k4me3_p ~ "P",
        old_type == "P" & h3k4me1_p & h3k4me3_not_called &
          (h3k27ac_p | atac_p) ~ "E",
        # 8. E/G dual-positive with bigWig resolution (non-silenced anchors)
        old_type == "E" & h3k4me1_p & h3k4me3_p &
          dual_ratio_state == "me1_dominant" ~ "dual",
        old_type == "E" & h3k4me1_p & h3k4me3_p &
          dual_ratio_state == "unresolved" ~ old_type,
        old_type == "E" & h3k4me1_p & h3k4me3_p ~ "P",
        old_type == "E" & h3k4me3_p & h3k4me1_not_called ~ "P",
        old_type == "G" & h3k4me1_p & h3k4me3_p &
          dual_ratio_state == "me1_dominant" ~ "dual",
        old_type == "G" & h3k4me1_p & h3k4me3_p &
          dual_ratio_state == "unresolved" ~ old_type,
        old_type == "G" & h3k4me1_p & h3k4me3_p ~ "P",
        old_type == "G" & h3k4me3_p & h3k4me1_not_called ~ "P",
        # 9. G -> E: intronic enhancer (H3K4me1+ H3K4me3- + H3K27ac/ATAC)
        old_type == "G" & h3k4me1_p & h3k4me3_not_called &
          (h3k27ac_p | atac_p) ~ "E",
        TRUE ~ old_type
      ),
      changed = new_type != old_type
    ) %>%
    dplyr::select(
      -h3k4me1_p, -h3k4me3_p, -h3k4me3_not_called, -h3k4me1_not_called,
      -h3k27ac_p, -h3k27me3_p, -atac_p,
      -is_promoter_like, -conf_chr
    )
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

#' Internal: Restore Gene Symbols for Anchors Reclassified to P/G/dual
#'
#' When upstream processing cleared or never assigned a gene symbol, but
#' chromatin evidence later reclassifies the anchor to P/G (H3K4me3+) or
#' dual (H3K4me1+ H3K4me3+), the original gene symbol must be restored.
#' This covers eP/eG ->P/G/dual (expression refinement cleared the gene),
#' E ->P (enhancers reclassified as promoters by chromatin evidence), and
#' G ->P (gene bodies reclassified as promoters). Without restoration, the
#' reclassified anchor is promoter-like by type but gene-less -- it is
#' silently dropped from promoter-centric stats and target-gene assignment.
#'
#' @param loop_df Loop annotation data frame with anchor1_gene/anchor2_gene.
#' @param reclass_map Reclassification map from \code{\link{.chromatin_reclassify}}.
#' @param annotation_res The full annotation result (may carry
#'   \code{looplook_anchor_state} attribute with original gene symbols).
#' Internal: Re-annotate TSS-proximal genes for E->P / G->P anchors
#'
#' For anchors reclassified from E or G to P, use the TxDb to find the
#' nearest TSS-proximal gene.  Used by \code{\link{.chromatin_restore_genes}}
#' when the original annotation did not record a promoter gene.
#'
#' @param anchor_ids Character vector of anchor IDs needing re-annotation.
#' @param validation Data frame with columns chr, start, end, anchor_id.
#' @param txdb_obj A TxDb object.
#' @return Named character vector of gene symbols (NA where no gene found).
#' @keywords internal
#' @noRd
.reannotate_tss_genes <- function(anchor_ids, validation, txdb_obj,
                                  org_db_pkg = NULL,
                                  tss_region = c(-2000, 2000),
                                  gene_before = NULL,
                                  gene_expr_map = NULL,
                                  min_expr = 0,
                                  conflict_strategy = c("biotype_first", "expression_first"),
                                  co_dominance_ratio = 0.1,
                                  biotype_order = c("protein", "small_ncRNA", "antisense", "lncRNA", "pseudogene"),
                                  quiet = FALSE) {
  val_sub <- validation[match(anchor_ids, validation$anchor_id), ]
  gr <- GenomicRanges::GRanges(
    seqnames = val_sub$chr,
    ranges   = IRanges::IRanges(val_sub$start, val_sub$end),
    strand   = "*"
  )
  # Explicitly store the anchor ID as a metadata column so that it
  # survives ChIPseeker's output format.  GRanges names() alone are
  # not sufficient -- ChIPseeker may or may not preserve them as the
  # first column of the data frame.
  S4Vectors::mcols(gr)$looplook_anchor_id <- val_sub$anchor_id
  # Delegate to ChIPseeker::annotatePeak -- same engine used by the main
  # annotation pipeline.  It handles strand, multiple transcripts,
  # distance-to-TSS, and OrgDb SYMBOL mapping correctly.
  anno <- suppressMessages(
    if (!is.null(org_db_pkg) && nzchar(org_db_pkg)) {
      ChIPseeker::annotatePeak(
        gr,
        TxDb = txdb_obj, tssRegion = tss_region,
        annoDb = org_db_pkg, verbose = FALSE
      )
    } else {
      ChIPseeker::annotatePeak(
        gr,
        TxDb = txdb_obj, tssRegion = tss_region,
        verbose = FALSE
      )
    }
  )
  df <- as.data.frame(anno)
  # Align ChIPseeker output to input GRanges.  Three-tier strategy:
  #   1. metadata_id:     explicit looplook_anchor_id in mcols (preferred)
  #   2. coordinate_join: match by seqnames/start/end (fallback)
  #   3. unmatched:       rows without ID or coordinate match -> NA
  input_names <- val_sub$anchor_id
  alignment_method <- rep(NA_character_, length(input_names))

  if ("looplook_anchor_id" %in% colnames(df)) {
    out_names <- as.character(df$looplook_anchor_id)
    use_metadata <- TRUE
  } else {
    # ChIPseeker dropped the metadata column (e.g. older version).
    # Fall back to the first column; coordinate join will validate.
    out_names <- as.character(df[[1]])
    use_metadata <- FALSE
  }

  # Legacy id_match test (only when metadata column is absent).
  # When the ChIPseeker first column preserves GRanges names, this
  # short-circuits coordinate join.  Defined here because the branch
  # is shared with the metadata-ID path above.
  ok_id <- !use_metadata && length(out_names) == length(input_names) &&
    all(out_names %in% input_names) &&
    !anyDuplicated(out_names)

  # Prefer explicit metadata ID matching.  Unlike coordinate join, this
  # correctly distinguishes anchors with identical coordinates but
  # different anchor_id provenance.
  if (use_metadata) {
    id_match_idx <- match(input_names, out_names)
    n_matched <- sum(!is.na(id_match_idx))
    if (n_matched == length(input_names)) {
      # All anchors matched by metadata ID -- reorder and return.
      df <- df[id_match_idx, , drop = FALSE]
      alignment_method <- rep("metadata_id", length(input_names))
      out_names <- input_names
    } else {
      # Partial ID match: align what we can, mark the rest unmatched.
      alignment_method <- ifelse(
        !is.na(id_match_idx), "metadata_id", "unmatched"
      )
      df_clean <- df[1L, , drop = FALSE]
      df_clean[1L, ] <- NA
      df_clean <- df_clean[rep(1L, length(input_names)), , drop = FALSE]
      rownames(df_clean) <- NULL
      ok_rows <- which(!is.na(id_match_idx))
      df_clean[ok_rows, ] <- df[id_match_idx[ok_rows], , drop = FALSE]
      df <- df_clean
      out_names <- input_names
      warning(
        "ChIPseeker output metadata IDs partially matched: ",
        n_matched, "/", length(input_names),
        " anchors aligned by looplook_anchor_id, ",
        sum(is.na(id_match_idx)), " are unmatched and set to NA.",
        call. = FALSE
      )
    }
  } else if (ok_id) {
    # Legacy path: ChIPseeker first column preserved GRanges names.
    # Reorder to match input order (ChIPseeker may change row order).
    id_idx <- match(input_names, out_names)
    df <- df[id_idx, , drop = FALSE]
    out_names <- input_names
    alignment_method <- rep("id_match", length(input_names))
  }

  # Coordinate join: used when metadata ID matching is unavailable or
  # the legacy id_match fails.
  if (all(is.na(alignment_method)) && all(c("seqnames", "start", "end") %in% colnames(df))) {
    # Build coordinate keys for both sides
    df_coord <- paste0(df$seqnames, ":", df$start, "-", df$end)
    gr_coord <- paste0(
      GenomicRanges::seqnames(gr), ":",
      GenomicRanges::start(gr), "-",
      GenomicRanges::end(gr)
    )
    # Detect duplicate coordinates -- ambiguous assignment, mark as unmatched.
    # match() returns the first hit, which is silently wrong for duplicates.
    dup_input <- duplicated(gr_coord) | duplicated(gr_coord, fromLast = TRUE)
    dup_output <- duplicated(df_coord) | duplicated(df_coord, fromLast = TRUE)
    n_dup_input <- sum(dup_input)
    n_dup_output <- sum(dup_output)
    if (n_dup_input > 0 || n_dup_output > 0) {
      warning(
        "Duplicate coordinates detected in TSS alignment: ",
        n_dup_input, " input anchor(s) and ", n_dup_output,
        " ChIPseeker output row(s). Rows with duplicate coordinates ",
        "will be marked as unmatched (gene_after=NA). ",
        "Ensure that all anchor coordinates are unique ",
        "before running TSS annotation.",
        call. = FALSE
      )
    }
    coord_match_raw <- match(gr_coord, df_coord)
    # Save ambiguity mask BEFORE clearing matches so the
    # alignment_method labels are accurate.
    ambiguous_from_input <- dup_input
    ambiguous_from_output <- !is.na(coord_match_raw) &
      dup_output[coord_match_raw]
    ambiguous <- ambiguous_from_input | ambiguous_from_output
    # Invalidate matches that involve duplicate coordinates
    coord_match <- coord_match_raw
    coord_match[ambiguous] <- NA_integer_
    n_matched <- sum(!is.na(coord_match))
    if (n_matched == length(input_names) &&
      !anyDuplicated(df_coord) && !anyDuplicated(gr_coord) && !any(ambiguous)) {
      # Perfect coordinate match -- reorder df rows to match input order
      df <- df[coord_match, , drop = FALSE]
      out_names <- input_names
      alignment_method <- rep("coordinate_join", length(input_names))
    } else if (n_matched > 0 || any(ambiguous)) {
      # Partial match or duplicates detected: align matched rows by
      # coordinate.  Unmatched and duplicate-ambiguous rows are
      # cleared to NA.
      out_names <- input_names
      alignment_method <- dplyr::case_when(
        !is.na(coord_match) ~ "coordinate_join",
        ambiguous ~ "duplicate_coordinate_ambiguous",
        TRUE ~ "unmatched"
      )
      ok_rows <- which(!is.na(coord_match))
      # Build a clean template: all rows NA, then fill matched ones
      df_clean <- df[1L, , drop = FALSE]
      df_clean[1L, ] <- NA
      df_clean <- df_clean[rep(1L, length(input_names)), , drop = FALSE]
      rownames(df_clean) <- NULL
      df_clean[ok_rows, ] <- df[coord_match[ok_rows], , drop = FALSE]
      df <- df_clean
      warning(
        "ChIPseeker output peak IDs do not match input GRanges names. ",
        "Partial coordinate-based alignment: ", n_matched, "/",
        length(input_names), " anchors matched by position, ",
        sum(is.na(coord_match)), " are unmatched and set to NA. ",
        "Check that the TxDb covers all anchor coordinates, ",
        "or set tss_region to a wider interval.",
        call. = FALSE
      )
    } else {
      # Zero coordinate matches: clear all rows to NA rather than
      # trusting row-order alignment, which would silently assign
      # wrong genes to anchors.
      out_names <- input_names
      alignment_method <- ifelse(ambiguous,
        "duplicate_coordinate_ambiguous", "unmatched"
      )
      df <- df[rep(1L, length(input_names)), , drop = FALSE]
      df[] <- NA
      rownames(df) <- NULL
      warning(
        "ChIPseeker output peak IDs do not match input GRanges names, ",
        "and coordinate-based alignment yielded zero matches. ",
        "All ", length(input_names), " anchors are set to NA ",
        "(gene_after=NA, TSS_supported=FALSE). ",
        "Check that the TxDb covers all anchor coordinates, ",
        "or set tss_region to a wider interval.",
        call. = FALSE
      )
    }
  } else if (all(is.na(alignment_method))) {
    # Coordinates not available in ChIPseeker output -- cannot
    # validate alignment at all.  Set all rows to NA.
    out_names <- input_names
    alignment_method <- rep("unmatched", length(input_names))
    df <- df[rep(1L, length(input_names)), , drop = FALSE]
    df[] <- NA
    rownames(df) <- NULL
    warning(
      "ChIPseeker output peak IDs do not match input GRanges names, ",
      "and coordinates are not available for fallback matching. ",
      "All ", length(input_names), " anchors are set to NA. ",
      "This may indicate a ChIPseeker version change or custom TxDb.",
      call. = FALSE
    )
  }
  # Only keep genes where ChIPseeker classifies the anchor as "Promoter"
  # (i.e. within tss_region of a TSS).  For Distal Intergenic or other
  # categories, the nearest gene may be far from the anchor and is NOT
  # a valid TSS-proximal assignment.
  is_promoter <- grepl("^Promoter", df$annotation)
  # Resolve multi-gene conflicts using the same engine as basic annotation.
  # ChIPseeker SYMBOL may contain "geneA/geneB" for overlapping TSS regions;
  # resolve_gene_conflicts applies biotype priority + expression tiebreaker.
  if (any(is_promoter, na.rm = TRUE) && !is.null(org_db_pkg) && nzchar(org_db_pkg)) {
    resolved_df <- format_annotation_columns(df)
    resolved_df <- resolve_gene_conflicts(
      resolved_df, txdb_obj, org_db_pkg, tss_region,
      gene_expr_map = gene_expr_map,
      min_expr = min_expr,
      conflict_strategy = conflict_strategy,
      co_dominance_ratio = co_dominance_ratio,
      biotype_order = biotype_order,
      unmeasured_policy = "keep"
    )
    if ("SYMBOL" %in% colnames(resolved_df)) {
      df$SYMBOL <- resolved_df$SYMBOL
    }
    # Sync annotation from resolved result and recalculate is_promoter.
    # resolve_gene_conflicts may upgrade the annotation (e.g. from
    # "Exon" to "Promoter" when the anchor spans a TSS), so is_promoter
    # must reflect the resolved state, not the raw ChIPseeker output.
    if ("annotation" %in% colnames(resolved_df)) {
      df$annotation <- resolved_df$annotation
    }
    is_promoter <- grepl("^Promoter", df$annotation)
  }
  # Per-row SYMBOL ->geneId fallback.  When OrgDb is unavailable
  # (non-model species) or missing for some rows, ChIPseeker may
  # only populate geneId.  Guard against missing columns: coalesce
  # with a zero-length vector is undefined.
  symbol_vec <- if ("SYMBOL" %in% colnames(df)) {
    as.character(df$SYMBOL)
  } else {
    rep(NA_character_, nrow(df))
  }
  geneid_vec <- if ("geneId" %in% colnames(df)) {
    as.character(df$geneId)
  } else {
    rep(NA_character_, nrow(df))
  }
  sym_raw <- dplyr::coalesce(symbol_vec, geneid_vec)
  gene_after <- sym_raw
  gene_after[!is_promoter] <- NA_character_
  gene_before_vec <- if (!is.null(gene_before)) {
    gene_before[out_names]
  } else {
    rep(NA_character_, length(out_names))
  }
  gene_before_vec[is.na(gene_before_vec)] <- NA_character_

  result <- data.frame(
    anchor_id = out_names,
    gene_before = gene_before_vec,
    gene_after = gene_after,
    annotation = df$annotation,
    distanceToTSS = if (!is.null(df$distanceToTSS)) df$distanceToTSS else NA_integer_,
    TSS_supported = is_promoter,
    assignment_method = ifelse(is_promoter, "ChIPseeker_TxDb", "no_TSS_support"),
    alignment_method = alignment_method,
    stringsAsFactors = FALSE
  )
  rownames(result) <- NULL

  n_promoter <- sum(is_promoter, na.rm = TRUE)
  if (n_promoter < length(anchor_ids) && !quiet) {
    message(sprintf(
      "    TSS re-annotation: %d / %d anchors assigned a promoter gene (others lack TSS support)",
      n_promoter, length(anchor_ids)
    ))
  }
  result
}

#' @return The \code{loop_df} with gene symbols restored where applicable.
#' @keywords internal
#' @noRd
.chromatin_restore_genes <- function(loop_df, reclass_map, annotation_res,
                                     validation = NULL, txdb_obj = NULL,
                                     org_db_pkg = NULL,
                                     tss_region = c(-2000, 2000),
                                     gene_expr_map = NULL,
                                     min_expr = 0,
                                     conflict_strategy = c("biotype_first", "expression_first"),
                                     co_dominance_ratio = 0.1,
                                     biotype_order = c("protein", "small_ncRNA", "antisense", "lncRNA", "pseudogene"),
                                     quiet = FALSE) {
  # All return paths must emit list(loop_df = ..., tss_reannotation = ...)
  # so the caller can always destructure with $loop_df / $tss_reannotation.
  tss_reann_out <- NULL

  restore <- reclass_map[
    reclass_map$changed &
      reclass_map$old_type %in% c("eP", "eG", "E", "G") &
      reclass_map$new_type %in% c("P", "G", "dual"),
  ]
  if (nrow(restore) == 0L) {
    return(list(
      loop_df = loop_df, tss_reannotation = NULL,
      tss_assignment_provenance = .empty_tss_assignment_provenance()
    ))
  }

  anchor_state <- attr(annotation_res, "looplook_anchor_state", exact = TRUE)
  if (is.null(anchor_state) || !"map_info" %in% names(anchor_state)) {
    warning(
      "Anchor state is unavailable; TSS assignment provenance ",
      "cannot be reconstructed. Chromatin types are retained, ",
      "but promoter-gene confidence is unresolved.",
      call. = FALSE
    )
    return(list(
      loop_df = loop_df, tss_reannotation = NULL,
      tss_assignment_provenance = .empty_tss_assignment_provenance()
    ))
  }
  gene_lookup <- setNames(
    anchor_state$map_info$SYMBOL,
    anchor_state$map_info$anchor_id
  )
  # Fall back to positional_gene where SYMBOL is NA
  if ("positional_gene" %in% colnames(anchor_state$map_info)) {
    pos_lookup <- setNames(
      anchor_state$map_info$positional_gene,
      anchor_state$map_info$anchor_id
    )
    na_idx <- is.na(gene_lookup) | gene_lookup == ""
    gene_lookup[na_idx] <- pos_lookup[names(gene_lookup)[na_idx]]
  }

  # Anchors that must be re-annotated against the TxDb:
  # 1. Any anchor whose gene is missing/empty
  # 2. E->P and E->dual: the positional gene of an E anchor is a distal
  #    nearest gene from ChIPseeker, NOT a TSS-validated promoter gene.
  #    Reusing it would assign an untested gene to a newly classified
  #    promoter, contaminating promoter stats and target assignment.
  requires_tss <- restore$old_type %in% c("E", "G", "eG") &
    restore$new_type %in% c("P", "dual")
  # Save original gene before clearing, so TSS provenance can report
  # the full before->after transition (e.g. distal GeneA ->TSS GeneB).
  gene_lookup_before_tss <- gene_lookup
  # Re-annotate TSS for E/G/eG->P/dual.  Original gene is kept:
  # if TxDb finds a TSS-supported gene it replaces the old one;
  # otherwise the original gene stays (warning issued below).
  # untested gene to a newly classified promoter.
  need_reannotate <- unique(c(
    restore$anchor_id[is.na(gene_lookup[restore$anchor_id]) |
      gene_lookup[restore$anchor_id] == ""],
    restore$anchor_id[requires_tss]
  ))
  if (length(need_reannotate) > 0 && !is.null(txdb_obj) &&
    !is.null(validation) && "chr" %in% colnames(validation)) {
    reann_df <- .reannotate_tss_genes(
      need_reannotate, validation, txdb_obj, org_db_pkg, tss_region,
      gene_before = gene_lookup_before_tss[need_reannotate],
      gene_expr_map = gene_expr_map,
      min_expr = min_expr,
      conflict_strategy = conflict_strategy,
      co_dominance_ratio = co_dominance_ratio,
      biotype_order = biotype_order,
      quiet = quiet
    )
    # Update gene_lookup with TSS-annotated genes (only where supported)
    reann_hits <- reann_df[reann_df$TSS_supported &
      !is.na(reann_df$gene_after) &
      nzchar(reann_df$gene_after), ]
    gene_lookup[reann_hits$anchor_id] <- reann_hits$gene_after
    n_found <- nrow(reann_hits)
    if (n_found > 0 && !quiet) {
      message(sprintf(
        "    TSS re-annotation: %d / %d gene-less E->P / G->P anchors assigned genes via TxDb",
        n_found, length(need_reannotate)
      ))
    }
    # Store for stable output (used below in out$tss_reannotation)
    tss_reann_out <- reann_df
  } else {
    tss_reann_out <- NULL
  }

  # Build TSS provenance columns for all anchors
  restore$TSS_supported <- NA
  restore$TSS_support_status <- "not_required"
  restore$Gene_Assignment_Confidence <- "structural"
  restore$Gene_Assignment_Evidence <- "positional_annotation"

  # Backfill requires_tss anchors that didn't get TxDb validation:
  # mark as TSS_validation_unavailable instead of not_required
  if (.is_null_or_empty(tss_reann_out)) {
    no_txdb_idx <- which(requires_tss)
    if (length(no_txdb_idx) > 0) {
      restore$TSS_support_status[no_txdb_idx] <-
        "TSS_validation_unavailable"
      has_gene <- !is.na(gene_lookup[restore$anchor_id[no_txdb_idx]]) &
        gene_lookup[restore$anchor_id[no_txdb_idx]] != ""
      restore$Gene_Assignment_Confidence[no_txdb_idx] <-
        ifelse(has_gene, "provisional", "unresolved")
      is_g_origin <- restore$old_type[no_txdb_idx] %in% c("G", "eG")
      restore$Gene_Assignment_Evidence[no_txdb_idx] <-
        ifelse(has_gene,
          ifelse(is_g_origin,
            "gene_body_host+promoter_like_chromatin",
            "positional_gene+promoter_like_chromatin"
          ),
          "promoter_like_chromatin_only"
        )
    }
  }

  # Mark TSS-reannotated anchors with their support status
  if (!is.null(tss_reann_out) && nrow(tss_reann_out) > 0) {
    tss_map <- setNames(
      tss_reann_out$TSS_supported,
      tss_reann_out$anchor_id
    )
    reann_genes <- setNames(
      tss_reann_out$gene_after,
      tss_reann_out$anchor_id
    )
    for (j in seq_len(nrow(restore))) {
      aid <- restore$anchor_id[j]
      if (!aid %in% names(tss_map)) next
      if (isTRUE(tss_map[aid]) &&
        !is.na(reann_genes[aid]) && nzchar(reann_genes[aid])) {
        # Annotated TSS found and gene assigned
        restore$TSS_supported[j] <- TRUE
        restore$TSS_support_status[j] <- "annotated_TSS_supported"
        restore$Gene_Assignment_Confidence[j] <- "high"
        restore$Gene_Assignment_Evidence[j] <- "annotated_TSS"
      } else if (requires_tss[j] &&
        !is.na(gene_lookup[aid]) && gene_lookup[aid] != "") {
        # No annotated TSS, but original gene retained.
        # Distinguish G/eG structural host gene from E positional nearest gene.
        restore$TSS_supported[j] <- FALSE
        restore$TSS_support_status[j] <- "no_annotated_TSS"
        restore$Gene_Assignment_Confidence[j] <- "provisional"
        is_g <- restore$old_type[j] %in% c("G", "eG")
        restore$Gene_Assignment_Evidence[j] <-
          if (is_g) "gene_body_host+promoter_like_chromatin"
          else "positional_gene+promoter_like_chromatin"
      } else if (requires_tss[j]) {
        # No TSS, no gene: unresolved
        restore$TSS_supported[j] <- FALSE
        restore$TSS_support_status[j] <- "no_annotated_TSS"
        restore$Gene_Assignment_Confidence[j] <- "unresolved"
        restore$Gene_Assignment_Evidence[j] <-
          "promoter_like_chromatin_only"
      }
    }
  }

  # Role correction for non-promoter anchors that became P/dual without TSS
  # is now handled in .chromatin_annotate_links() after TSS provenance columns
  # are available.  The gene_lookup in this function retains all genes so that
  # loop_df anchor_gene columns are preserved for downstream role correction.

  # Warn about TSS anchors that kept their original gene (no TSS found)
  tss_kept_original <- which(
    restore$TSS_support_status == "no_annotated_TSS" &
      restore$Gene_Assignment_Confidence == "provisional"
  )
  tss_not_found <- which(
    restore$TSS_support_status == "no_annotated_TSS" &
      restore$Gene_Assignment_Confidence == "unresolved"
  )
  if (length(tss_kept_original) > 0) {
    warning(
      sprintf(
        paste0(
          "%d anchor(s) reclassified to P/dual but no TSS-validated ",
          "gene was found by TxDb. The original positional gene is ",
          "for downstream role-aware interpretation. ",
          "Role-specific classification (gene_body / positional_candidate) ",
          "is reported during chromatin link annotation. ",
          "Affected: %s"
        ),
        length(tss_kept_original),
        paste(head(restore$anchor_id[tss_kept_original], 5), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  if (length(tss_not_found) > 0) {
    warning(
      sprintf(
        paste0(
          "%d promoter-like anchor(s) have neither an annotated TSS ",
          "nor a positional candidate gene. Promoter assignment ",
          "confidence is unresolved; verify with orthogonal data. ",
          "Affected: %s"
        ),
        length(tss_not_found),
        paste(head(restore$anchor_id[tss_not_found], 5), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  n_restored <- 0L
  n_skipped <- 0L
  for (i in seq_len(nrow(restore))) {
    aid <- restore$anchor_id[i]
    orig_gene <- gene_lookup[aid]
    # E/G/eG->P/dual anchors: gene was re-annotated against TxDb.
    # If TSS found a better gene, gene_lookup already has it.
    # If not, original gene is kept (warning emitted above).
    # Ordinary restoration fills empty genes; TSS anchors force
    # overwrite even when already non-empty.
    is_tss_anchor <- aid %in% restore$anchor_id[requires_tss]

    if (is.na(orig_gene) || orig_gene == "") {
      # Anchor genuinely has no gene: nothing to restore.
      # TSS anchors that couldn't get a gene are counted above
      # in tss_not_found with the consolidated warning.
      if (is_tss_anchor) {
        if ("a1_id" %in% colnames(loop_df)) {
          loop_df$anchor1_gene[loop_df$a1_id == aid] <- NA_character_
        }
        if ("a2_id" %in% colnames(loop_df)) {
          loop_df$anchor2_gene[loop_df$a2_id == aid] <- NA_character_
        }
        n_restored <- n_restored + sum(loop_df$a1_id == aid | loop_df$a2_id == aid, na.rm = TRUE)
      } else {
        n_skipped <- n_skipped + 1L
      }
      next
    }

    # anchor1
    if ("a1_id" %in% colnames(loop_df)) {
      if (is_tss_anchor) {
        loop_df$anchor1_gene[loop_df$a1_id == aid] <- orig_gene
        n_restored <- n_restored + sum(loop_df$a1_id == aid, na.rm = TRUE)
      } else {
        hits1 <- which(loop_df$a1_id == aid &
          (is.na(loop_df$anchor1_gene) | loop_df$anchor1_gene == ""))
        if (length(hits1) > 0) {
          loop_df$anchor1_gene[hits1] <- orig_gene
          n_restored <- n_restored + length(hits1)
        }
      }
    }
    # anchor2
    if ("a2_id" %in% colnames(loop_df)) {
      if (is_tss_anchor) {
        loop_df$anchor2_gene[loop_df$a2_id == aid] <- orig_gene
        n_restored <- n_restored + sum(loop_df$a2_id == aid, na.rm = TRUE)
      } else {
        hits2 <- which(loop_df$a2_id == aid &
          (is.na(loop_df$anchor2_gene) | loop_df$anchor2_gene == ""))
        if (length(hits2) > 0) {
          loop_df$anchor2_gene[hits2] <- orig_gene
          n_restored <- n_restored + length(hits2)
        }
      }
    }
  }

  if ((n_restored > 0 || n_skipped > 0) && !quiet) {
    message(sprintf(
      "    Gene restoration: %d entries restored, %d anchors skipped (no gene in anchor_state).",
      n_restored, n_skipped
    ))
    if (n_skipped > 0) {
      skipped_types <- unique(restore$old_type[
        vapply(restore$anchor_id, function(aid) {
          g <- gene_lookup[aid]
          is.na(g) || g == ""
        }, logical(1))
      ])
      message(sprintf(
        "    Skipped anchor types: %s. E->P anchors typically lack genes (enhancers do not overlap annotated TSS) and will be absent from promoter_centric_stats.",
        paste(skipped_types, collapse = ", ")
      ))
    }
  }
  # Build TSS provenance from the restore data frame
  txdb_gene_map <- if (!is.null(tss_reann_out) && nrow(tss_reann_out) > 0L) {
    setNames(tss_reann_out$gene_after, tss_reann_out$anchor_id)
  } else {
    character()
  }
  tss_provenance_out <- restore %>%
    dplyr::transmute(
      anchor_id,
      positional_type = old_type,
      final_type = new_type,
      gene_before = unname(gene_lookup_before_tss[anchor_id]),
      TxDb_gene = unname(txdb_gene_map[anchor_id]),
      final_gene = unname(gene_lookup[anchor_id]),
      TSS_supported,
      TSS_support_status,
      Gene_Assignment_Confidence,
      Gene_Assignment_Evidence
    )
  list(
    loop_df = loop_df,
    tss_reannotation = tss_reann_out,
    tss_assignment_provenance = tss_provenance_out
  )
}

#' Internal: Recompute loop_type from updated anchor types, and
#' rebuild Putative_Target_Genes and Cluster_All_Genes from
#' current anchor types and genes after chromatin refinement.
#'
#' Rebuilds loop-level Putative_Target_Genes and Promoter_Target_Genes
#' from post-chromatin anchor states and stored topology.  When expression
#' provenance is available, the rebuilt target sets are re-filtered using
#' the stored effective expression threshold.
#'
#' @param loop_df Loop annotation data frame with updated anchor columns.
#' @param expr_vals Named numeric vector of per-gene expression means.
#' @param expr_threshold Expression activity threshold.
#' @param ego_list_loop Per-anchor ego network vertex lists.
#' @param map_info Anchor-level map with type_code and SYMBOL.
#' @param graph igraph object for vertex ID resolution.
#' @return The loop_df with Putative_Target_Genes, Promoter_Target_Genes,
#'   All_Anchor_Genes, Cluster_All_Genes, and related columns rebuilt.
#' @keywords internal
#' @noRd
.rebuild_loop_gene_summaries <- function(loop_df, expr_vals = NULL,
                                         expr_threshold = NULL,
                                         ego_list_loop = NULL,
                                         map_info = NULL,
                                         graph = NULL) {
  if (!all(c("anchor1_type", "anchor2_type", "anchor1_gene", "anchor2_gene") %in%
    colnames(loop_df))) {
    return(loop_df)
  }
  is_promoter <- .is_target_promoter_like
  is_gene_body <- .is_target_gene_body_like
  has_roles <- all(c("anchor1_gene_role", "anchor2_gene_role") %in%
    colnames(loop_df))

  # All_Anchor_Genes: strictly structural -- promoter-role and gene_body-role
  # anchors only.  This reflects the anchor's own classification, not whether
  # the gene is valid.  After chromatin refinement:
  #   - P/eP/G/eG -> E downgrades produce enhancer_candidate (strict=TRUE).
  #     Their genes are structurally valid but the anchor is now enhancer-class,
  #     so they are excluded from All_Anchor_Genes.  They remain visible in
  #     Putative_Target_Genes (topology expansion, strict=TRUE),
  #     All_Loop_Connected_Genes, and compete in Assigned_Target_Genes at rank 3.
  #   - E anchors (any role): ChIPseeker nearest genes are NOT structural
  #     anchor genes and are excluded regardless of role.
  #   - Positional candidates are excluded from All_Anchor_Genes but remain
  #     in All_Loop_Connected_Genes and Candidate_Positional_Genes.
  has_strict <- all(c("anchor1_strict_eligible", "anchor2_strict_eligible") %in%
    colnames(loop_df))
  if (has_roles) {
    loop_df <- loop_df %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        .p1 = if (isTRUE(anchor1_gene_role == "promoter") &&
                  (!has_strict || isTRUE(anchor1_strict_eligible)))
              extract_genes(anchor1_gene) else NA_character_,
        .p2 = if (isTRUE(anchor2_gene_role == "promoter") &&
                  (!has_strict || isTRUE(anchor2_strict_eligible)))
              extract_genes(anchor2_gene) else NA_character_,
        .g1 = if (isTRUE(anchor1_gene_role == "gene_body") &&
                  (!has_strict || isTRUE(anchor1_strict_eligible)))
              extract_genes(anchor1_gene) else NA_character_,
        .g2 = if (isTRUE(anchor2_gene_role == "gene_body") &&
                  (!has_strict || isTRUE(anchor2_strict_eligible)))
              extract_genes(anchor2_gene) else NA_character_,
        All_Anchor_Genes = extract_genes(c(.p1, .p2, .g1, .g2))
      ) %>%
      dplyr::ungroup()
  } else {
    loop_df <- loop_df %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        All_Anchor_Genes = .loop_locus_genes(
          anchor1_type, anchor2_type,
          anchor1_gene, anchor2_gene
        ),
        .p1 = if (is_promoter(anchor1_type)) extract_genes(anchor1_gene) else NA_character_,
        .p2 = if (is_promoter(anchor2_type)) extract_genes(anchor2_gene) else NA_character_,
        .g1 = if (is_gene_body(anchor1_type)) extract_genes(anchor1_gene) else NA_character_,
        .g2 = if (is_gene_body(anchor2_type)) extract_genes(anchor2_gene) else NA_character_
      ) %>%
      dplyr::ungroup()
  }

  # Putative_Target_Genes: uses topology-expanded promoter- and
  # gene-body-associated candidates, with promoter-side restriction
  # for gene-body--promoter loops.  Direct anchor genes are used only
  # when the topology-derived set is empty.  E anchor nearest genes
  # are NOT targets.
  if (!is.null(ego_list_loop) && !is.null(map_info) &&
    all(c("a1_id", "a2_id") %in% colnames(loop_df))) {
    # Recompute neighbor_hop gene expansion from current map_info.
    # This mirrors the initial annotation's Putative_Target_Genes
    # computation, preserving neighbor_hop semantics.
    map_info$SYMBOL <- trimws(map_info$SYMBOL)
    has_roles <- "effective_gene_role" %in% colnames(map_info)
    has_map_strict <- "strict_assignment_eligible" %in% colnames(map_info)
    if (has_roles) {
      # Use strict eligibility to filter topology-expanded gene sources.
      # This excludes E->E nearest genes (strict=FALSE) and E->P/dual
      # unresolved (positional_candidate, strict=FALSE), while retaining
      # G->P/dual host genes and P/G->E structurally-supported candidates.
      if (has_map_strict) {
        valid_pg <- map_info %>%
          dplyr::filter(
            strict_assignment_eligible %in% TRUE,
            !is.na(SYMBOL), SYMBOL != ""
          )
        valid_p <- map_info %>%
          dplyr::filter(
            effective_gene_role == "promoter",
            strict_assignment_eligible %in% TRUE,
            !is.na(SYMBOL), SYMBOL != ""
          )
      } else {
        valid_pg <- map_info %>%
          dplyr::filter(
            effective_gene_role %in% c("promoter", "gene_body"),
            !is.na(SYMBOL), SYMBOL != ""
          )
        valid_p <- map_info %>%
          dplyr::filter(
            effective_gene_role == "promoter",
            !is.na(SYMBOL), SYMBOL != ""
          )
      }
    } else {
      valid_pg <- map_info %>%
        dplyr::filter(
          (.is_target_promoter_like(type_code) | .is_target_gene_body_like(type_code)) &
            !is.na(SYMBOL) & SYMBOL != ""
        )
      valid_p <- map_info %>%
        dplyr::filter(
          .is_target_promoter_like(type_code),
          !is.na(SYMBOL) & SYMBOL != ""
        )
    }
    lookup_pg <- setNames(valid_pg$SYMBOL, valid_pg$anchor_id)
    .anchor_topo <- function(aid) {
      n <- .igraph_vertex_ids(ego_list_loop[[aid]], graph)
      if (is.null(n)) n <- character(0)
      extract_genes(lookup_pg[c(aid, n)])
    }
    # Promoter-only lookup
    lookup_p <- setNames(valid_p$SYMBOL, valid_p$anchor_id)
    .anchor_topo_p <- function(aid) {
      n <- .igraph_vertex_ids(ego_list_loop[[aid]], graph)
      if (is.null(n)) n <- character(0)
      extract_genes(lookup_p[c(aid, n)])
    }
    topo_genes <- vapply(seq_len(nrow(loop_df)), function(i) {
      pg1 <- .anchor_topo(loop_df$a1_id[i])
      pg2 <- .anchor_topo(loop_df$a2_id[i])
      if (has_roles) {
        r1 <- loop_df$anchor1_gene_role[i]
        r2 <- loop_df$anchor2_gene_role[i]
      } else {
        r1 <- "other"
        r2 <- "other"
      }
      # G-P/P-G asymmetry: if one side is gene_body and the other is
      # promoter, use the promoter-side genes only.
      if ((has_roles && isTRUE(r1 == "gene_body") && isTRUE(r2 == "promoter")) ||
        (!has_roles && is_gene_body(loop_df$anchor1_type[i]) && is_promoter(loop_df$anchor2_type[i]))) {
        return(extract_genes(pg2))
      }
      if ((has_roles && isTRUE(r1 == "promoter") && isTRUE(r2 == "gene_body")) ||
        (!has_roles && is_promoter(loop_df$anchor1_type[i]) && is_gene_body(loop_df$anchor2_type[i]))) {
        return(extract_genes(pg1))
      }
      extract_genes(c(pg1, pg2))
    }, character(1))
    loop_df$.topo_genes <- topo_genes

    # Promoter-only topology, same asymmetry logic
    topo_promoter <- vapply(seq_len(nrow(loop_df)), function(i) {
      p_p1 <- .anchor_topo_p(loop_df$a1_id[i])
      p_p2 <- .anchor_topo_p(loop_df$a2_id[i])
      if (has_roles) {
        r1 <- loop_df$anchor1_gene_role[i]
        r2 <- loop_df$anchor2_gene_role[i]
      } else {
        r1 <- "other"
        r2 <- "other"
      }
      if ((has_roles && isTRUE(r1 == "gene_body") && isTRUE(r2 == "promoter")) ||
        (!has_roles && is_gene_body(loop_df$anchor1_type[i]) && is_promoter(loop_df$anchor2_type[i]))) {
        return(extract_genes(p_p2))
      }
      if ((has_roles && isTRUE(r1 == "promoter") && isTRUE(r2 == "gene_body")) ||
        (!has_roles && is_promoter(loop_df$anchor1_type[i]) && is_gene_body(loop_df$anchor2_type[i]))) {
        return(extract_genes(p_p1))
      }
      extract_genes(c(p_p1, p_p2))
    }, character(1))
    loop_df$.topo_promoter <- topo_promoter
  } else {
    loop_df$.topo_genes <- NA_character_
    loop_df$.topo_promoter <- NA_character_
  }

  # Build structural target columns.  When neighbor_hop data exists,
  # use topology-expanded genes; otherwise fall back to direct anchor
  # genes.  Promoter targets only use P/eP/dual anchors.
  # When effective_gene_role columns are present (post-chromatin),
  # they take precedence over type-based classification.
  loop_df <- loop_df %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      Putative_Target_Genes = dplyr::case_when(
        !is.na(.topo_genes) & nzchar(.topo_genes) ~ .topo_genes,
        !is.na(.p1) | !is.na(.p2) ~ extract_genes(c(.p1, .p2)),
        !is.na(.g1) | !is.na(.g2) ~ extract_genes(c(.g1, .g2)),
        TRUE ~ NA_character_
      ),
      Promoter_Target_Genes = dplyr::case_when(
        !is.na(.topo_promoter) & nzchar(.topo_promoter) ~
          extract_genes(c(.topo_promoter)),
        !is.na(.p1) | !is.na(.p2) ~ extract_genes(c(.p1, .p2)),
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-.p1, -.p2, -.g1, -.g2, -.topo_genes, -.topo_promoter)

  # Rebuild Cluster_All_Genes and Cluster_Locus_Genes from the new values
  if ("cluster_id" %in% colnames(loop_df)) {
    cluster_genes <- loop_df %>%
      dplyr::filter(!is.na(cluster_id)) %>%
      dplyr::group_by(cluster_id) %>%
      dplyr::summarise(
        Cluster_All_Genes = extract_genes(Putative_Target_Genes),
        Cluster_Locus_Genes = extract_genes(All_Anchor_Genes),
        .groups = "drop"
      )
    loop_df <- loop_df %>%
      dplyr::select(-dplyr::any_of(c("Cluster_All_Genes", "Cluster_Locus_Genes"))) %>%
      dplyr::left_join(cluster_genes, by = "cluster_id")
  }

  # Rebuild or clear expression-dependent fields.
  # When expression data is available (stored in anchor_state by prior
  # expression refinement), recompute active target genes from the new
  # Putative_Target_Genes.  Otherwise clear to NA -- chromatin alone
  # cannot assess transcriptional activity.
  has_expr <- !is.null(expr_vals) &&
    !is.null(expr_threshold)
  if (has_expr) {
    .filter_active <- function(gene_string) {
      genes <- clean_gene_names(gene_string, ";")
      if (length(genes) == 0) {
        return(list(genes = NA_character_, state = "no_target"))
      }
      keys <- toupper(genes)
      v <- expr_vals[keys]
      measured <- keys %in% names(expr_vals)
      active <- measured & .passes_expression_threshold(v, expr_threshold)
      n_measured <- sum(measured)
      n_total <- length(genes)
      if (any(active)) {
        state <- if (all(measured)) {
          "active"
        } else if (any(active)) {
          "active_partial"
        } else {
          "active"
        }
        return(list(
          genes = .collapse_genes(genes[active]),
          state = state
        ))
      }
      # None active -- distinguish silent from unmeasured
      if (n_measured == 0L) {
        return(list(genes = NA_character_, state = "unmeasured"))
      } else if (n_measured == n_total) {
        return(list(genes = NA_character_, state = "measured_silent"))
      } else {
        return(list(genes = NA_character_, state = "mixed"))
      }
    }
    loop_df <- loop_df %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        res = list(.filter_active(Putative_Target_Genes)),
        res_promoter = list(.filter_active(Promoter_Target_Genes)),
        Putative_Target_Genes = res$genes,
        Promoter_Target_Genes = res_promoter$genes,
        Target_Expression_State = res$state
      ) %>%
      dplyr::ungroup()
    loop_df$res <- NULL
    loop_df$res_promoter <- NULL
    loop_df$Has_Active_Target <- !is.na(loop_df$Putative_Target_Genes) &
      loop_df$Putative_Target_Genes != ""
    # unmeasured / mixed ->NA (not FALSE): we cannot make a negative
    # claim when some or all genes lack expression measurements
    unknown_idx <- loop_df$Target_Expression_State %in% c("unmeasured", "mixed")
    loop_df$Has_Active_Target[unknown_idx] <- NA
    loop_df$Retained_In_Functional_Network <- loop_df$Has_Active_Target
    loop_df$Refinement_Action <- dplyr::case_when(
      loop_df$Target_Expression_State == "active" ~ "chromatin_retained",
      loop_df$Target_Expression_State == "active_partial" ~ "chromatin_partial",
      loop_df$Target_Expression_State == "measured_silent" ~ "chromatin_silent",
      loop_df$Target_Expression_State == "unmeasured" ~ "chromatin_unmeasured",
      loop_df$Target_Expression_State == "mixed" ~ "chromatin_mixed",
      TRUE ~ "chromatin_no_target"
    )
    if ("cluster_id" %in% colnames(loop_df)) {
      clust_active <- loop_df %>%
        dplyr::filter(!is.na(cluster_id)) %>%
        dplyr::group_by(cluster_id) %>%
        dplyr::summarise(
          Cluster_All_Genes = extract_genes(Putative_Target_Genes),
          .groups = "drop"
        )
      loop_df <- loop_df %>%
        dplyr::select(-dplyr::any_of("Cluster_All_Genes")) %>%
        dplyr::left_join(clust_active, by = "cluster_id")
    }
  } else {
    # Chromatin-only: keep structural targets intact; unconditionally
    # create expression-status columns so all refinement objects have
    # a stable output schema.
    loop_df$Target_Expression_State <- "not_assessed"
    loop_df$Has_Active_Target <- NA
    loop_df$Retained_In_Functional_Network <- NA
    loop_df$Refinement_Action <- "not_assessed"
  }

  loop_df
}

.chromatin_recompute_loop_types <- function(loop_df) {
  if (!all(c("anchor1_type", "anchor2_type") %in% colnames(loop_df))) {
    return(loop_df)
  }
  loop_df$loop_type <- vapply(seq_len(nrow(loop_df)), function(i) {
    .loop_type_code(loop_df$anchor1_type[i], loop_df$anchor2_type[i])
  }, character(1))
  loop_df
}

#' Internal: Annotate target gene links with chromatin provenance columns
#'
#' Adds `anchor_type_before_chromatin`, `anchor_type_after_chromatin`,
#' `chromatin_target_action`, `chromatin_confidence`, and
#' `chromatin_evidence` columns to the target_gene_links data frame.
#'
#' @param new_links Data frame from `.build_target_gene_links`.
#' @param reclass_map Reclassification data frame from `.chromatin_reclassify`.
#' @param annotation_res The chromatin-refined annotation result (may contain
#'   `chromatin_validation`).
#' @return The `new_links` data frame with chromatin provenance columns appended.
#' @keywords internal
#' @noRd
.chromatin_annotate_links <- function(new_links, reclass_map, annotation_res,
                                      tss_assignment_provenance = NULL) {
  if (nrow(new_links) == 0 || !"anchor_id" %in% colnames(new_links)) {
    return(new_links)
  }
  type_before <- setNames(reclass_map$old_type, reclass_map$anchor_id)
  type_after <- setNames(reclass_map$new_type, reclass_map$anchor_id)
  changed_map <- setNames(reclass_map$changed, reclass_map$anchor_id)
  new_links$anchor_type_before_chromatin <- type_before[new_links$anchor_id]
  new_links$anchor_type_after_chromatin <- type_after[new_links$anchor_id]
  new_links$chromatin_target_action <- ifelse(
    !is.na(new_links$anchor_type_before_chromatin) &
      new_links$anchor_type_before_chromatin != new_links$anchor_type_after_chromatin,
    paste0(
      new_links$anchor_type_before_chromatin, "->",
      new_links$anchor_type_after_chromatin
    ),
    "unchanged"
  )
  if (!is.null(annotation_res$chromatin_validation)) {
    cv <- annotation_res$chromatin_validation
    conf_map <- setNames(as.character(cv$enhancer_evidence), cv$anchor_id)
    evid_map <- setNames(cv$evidence, cv$anchor_id)
    new_links$chromatin_confidence <- conf_map[new_links$anchor_id]
    new_links$chromatin_evidence <- evid_map[new_links$anchor_id]
  }
  # Join TSS provenance from chromatin restoration, keyed by anchor_id.
  # tss_assignment_provenance contains TSS_supported, TSS_support_status,
  # Gene_Assignment_Confidence, and Gene_Assignment_Evidence.
  # is_expanded_path is computed in the unified .mark_target_gene_link_membership()
  if (!is.null(tss_assignment_provenance) && nrow(tss_assignment_provenance) > 0) {
    prov_cols <- intersect(
      c(
        "anchor_id", "TSS_supported", "TSS_support_status",
        "Gene_Assignment_Confidence", "Gene_Assignment_Evidence"
      ),
      colnames(tss_assignment_provenance)
    )
    tss_map <- tss_assignment_provenance[, prov_cols, drop = FALSE] %>%
      dplyr::distinct(anchor_id, .keep_all = TRUE)
    idx <- match(new_links$anchor_id, tss_map$anchor_id)
    for (nm in setdiff(prov_cols, "anchor_id")) {
      new_links[[nm]] <- tss_map[[nm]][idx]
    }
  }

  # --- Origin-aware role correction and strict eligibility ---
  #
  # Uses the same .resolve_chromatin_gene_role() resolver as the anchor-level
  # (map_info / loop_df), guaranteeing consistent role/strict across all
  # output layers.  The resolver combines final chromatin type with origin
  # (old_type) and TSS validation to distinguish chromatin-state
  # classification from gene-identity confidence.
  if ("anchor_type_before_chromatin" %in% colnames(new_links) &&
    "anchor_type_after_chromatin" %in% colnames(new_links)) {

    has_chromatin <- !is.na(new_links$anchor_type_before_chromatin) &
      !is.na(new_links$anchor_type_after_chromatin)
    has_tss  <- !is.na(new_links$TSS_supported) & new_links$TSS_supported
    has_gene <- !is.na(new_links$gene) & nzchar(new_links$gene)

    resolved <- .resolve_chromatin_gene_role(
      old_type     = new_links$anchor_type_before_chromatin,
      final_type   = new_links$anchor_type_after_chromatin,
      tss_supported = has_tss,
      has_gene     = has_gene
    )

    idx <- which(has_chromatin)
    if (length(idx) > 0) {
      new_links$gene_role[idx] <- resolved$role[idx]
      new_links$strict_assignment_eligible[idx] <- resolved$strict[idx]
    }
  }

  # --- QC: count role transitions and eligibility changes ---
  n_promoter_gained <- length(unique(na.omit(new_links$anchor_id[
    new_links$gene_role == "promoter" &
      new_links$anchor_type_before_chromatin %in% c("G", "eG")
  ])))
  n_promoter_gained_tss <- length(unique(na.omit(new_links$anchor_id[
    new_links$gene_role == "promoter" &
      new_links$anchor_type_before_chromatin == "E"
  ])))
  n_promoter_lost <- length(unique(na.omit(new_links$anchor_id[
    new_links$gene_role != "promoter" &
      new_links$anchor_type_before_chromatin %in% c("P", "eP")
  ])))
  n_e_to_p_unresolved <- length(unique(na.omit(new_links$anchor_id[
    new_links$gene_role == "positional_candidate" &
      new_links$anchor_type_before_chromatin == "E" &
      new_links$anchor_type_after_chromatin %in% c("P", "dual") &
      !(!is.na(new_links$TSS_supported) & new_links$TSS_supported)
  ])))
  n_e_unchanged_blocked <- length(unique(na.omit(new_links$anchor_id[
    new_links$gene_role == "enhancer_candidate" &
      new_links$anchor_type_before_chromatin == "E" &
      new_links$anchor_type_after_chromatin == "E" &
      !is.na(new_links$strict_assignment_eligible) &
      !as.logical(new_links$strict_assignment_eligible)
  ])))
  n_enhancer_candidate <- length(unique(na.omit(new_links$anchor_id[
    new_links$gene_role == "enhancer_candidate"
  ])))
  n_positional <- length(unique(na.omit(new_links$anchor_id[
    !is.na(new_links$strict_assignment_eligible) &
      !as.logical(new_links$strict_assignment_eligible)
  ])))

  if (n_promoter_gained > 0) {
    message("  chromatin role: ", n_promoter_gained,
      " anchor(s) gained promoter role (G/eG -> P/dual, host-gene)")
  }
  if (n_promoter_gained_tss > 0) {
    message("  chromatin role: ", n_promoter_gained_tss,
      " anchor(s) gained promoter role (E -> P/dual, TSS-validated)")
  }
  if (n_e_to_p_unresolved > 0) {
    message("  chromatin role: ", n_e_to_p_unresolved,
      " anchor(s) E -> P/dual with unresolved gene identity ",
      "(kept as positional_candidate)")
  }
  if (n_promoter_lost > 0) {
    message("  chromatin role: ", n_promoter_lost,
      " anchor(s) lost promoter role (P/eP -> G/E)")
  }
  if (n_e_unchanged_blocked > 0) {
    message("  chromatin role: ", n_e_unchanged_blocked,
      " unchanged E anchor(s) blocked from strict target ",
      "(nearest-gene not loop-supported)")
  }
  if (n_enhancer_candidate > 0) {
    message("  chromatin role: ", n_enhancer_candidate,
      " anchor(s) assigned enhancer_candidate role")
  }
  if (n_positional > 0) {
    message("  chromatin role: ", n_positional,
      " anchor(s) excluded from strict assignment (positional_candidate)")
  }

  new_links
}

#' Internal: Clear stale 3D target assignment columns before chromatin recomputation
#'
#' Sets all strict and filled target gene columns to NA so that an empty
#' `new_links` result does not silently retain pre-chromatin values.
#'
#' @param bed_info Target annotation data frame, or NULL.
#' @return The input `bed_info` with strict columns cleared.
#' @keywords internal
#' @noRd
.clear_target_summary_columns <- function(bed_info) {
  if (is.null(bed_info)) {
    return(NULL)
  }
  strict_cols <- c(
    "All_Loop_Connected_Genes",
    "Regulated_promoter_genes",
    "Assigned_Target_Genes",
    "Expanded_Target_Genes",
    "Candidate_Positional_Genes",
    "All_Loop_Connected_Genes_Filled",
    "Regulated_promoter_genes_Filled",
    "Assigned_Target_Genes_Filled",
    "Regulated_promoter_Evidence",
    "Regulated_promoter_Fallback_Evidence"
  )
  for (col in strict_cols) {
    bed_info[[col]] <- NA_character_
  }
  bed_info
}

#' Internal: Aggregate gene assignment columns from chromatin-aware loop links
#'
#' Computes `All_Loop_Connected_Genes`, `Regulated_promoter_genes`, and
#' `Assigned_Target_Genes` by aggregating over gene_role and evidence fields.
#'
#' @param new_links Target gene links data frame.
#' @param bed_info Target annotation data frame.
#' @return The `bed_info` with aggregated columns joined in.
#' @keywords internal
#' @noRd
.aggregate_strict_targets <- function(new_links, bed_info,
                                      target_priority = c("promoter_then_distance", "distance_then_role"),
                                      max_primary_hop = 1L) {
  target_priority <- match.arg(target_priority)
  if (is.null(bed_info) || .is_null_or_empty(new_links)) {
    return(bed_info)
  }
  loop_links <- new_links %>%
    dplyr::filter(source == "loop_anchor", !is.na(gene), gene != "")

  if (nrow(loop_links) == 0) {
    return(bed_info)
  }

  all_agg <- loop_links %>%
    dplyr::group_by(input_id) %>%
    dplyr::summarise(
      All_Loop_Connected_Genes = .collapse_genes(gene),
      .groups = "drop"
    )

  # Normalise schema before any filtering: old objects may lack
  # strict_assignment_eligible.  Old objects lacking strict provenance are
  # normalised conservatively: only promoter and gene_body are presumed
  # eligible; enhancer_candidate and positional_candidate are not.
  if (!"strict_assignment_eligible" %in% colnames(loop_links)) {
    loop_links$strict_assignment_eligible <- NA
  }
  is_positional <- !is.na(loop_links$gene_role) &
    loop_links$gene_role == "positional_candidate"
  loop_links$strict_assignment_eligible[is_positional] <- FALSE
  # NA strict on old objects: conservative -- only promoter/gene_body
  # are presumed eligible; enhancer_candidate with unknown provenance
  # defaults to ineligible.
  na_strict <- is.na(loop_links$strict_assignment_eligible)
  if (any(na_strict)) {
    loop_links$strict_assignment_eligible[na_strict] <-
      loop_links$gene_role[na_strict] %in% c("promoter", "gene_body")
  }

  # For strict target aggregation, filter to eligible links only.
  strict_links <- loop_links %>%
    dplyr::filter(
      !is.na(gene_role) & gene_role != "positional_candidate",
      is.na(strict_assignment_eligible) | strict_assignment_eligible
    )

  # Compute path_rank and role_rank for strict-eligible links
  ranked_links <- strict_links %>%
    dplyr::mutate(
      role_rank = dplyr::case_when(
        gene_role == "promoter"           ~ 1L,
        gene_role == "gene_body"          ~ 2L,
        gene_role == "enhancer_candidate"  ~ 3L,
        TRUE                              ~ 4L
      ),
      path_rank = dplyr::coalesce(as.numeric(path_length), Inf)
    )

  # Split by path rank: primary (finite, <= max_primary_hop) vs expanded (finite, > max_primary_hop)
  has_path <- is.finite(ranked_links$path_rank)
  primary_links <- ranked_links %>% dplyr::filter(has_path, path_rank <= max_primary_hop)
  expanded_links <- ranked_links %>% dplyr::filter(has_path, path_rank > max_primary_hop)

  # Primary-only: regulated promoter genes (P1-1: use primary, not all loop_links)
  empty_promoter <- data.frame(
    input_id = character(),
    Regulated_promoter_genes = character(),
    stringsAsFactors = FALSE
  )
  promoter_agg <- empty_promoter
  if (nrow(primary_links) > 0) {
    promoter_agg <- primary_links %>%
      dplyr::filter(gene_role == "promoter") %>%
      dplyr::group_by(input_id) %>%
      dplyr::summarise(
        Regulated_promoter_genes = .collapse_genes(gene),
        .groups = "drop"
      )
  }

  # Expand-only: anti-join with primary gene pairs (P1-2: no duplicate genes)
  empty_expanded <- data.frame(
    input_id = character(),
    Expanded_Target_Genes = character(),
    stringsAsFactors = FALSE
  )
  expanded_agg <- empty_expanded
  if (nrow(expanded_links) > 0) {
    primary_gene_pairs <- primary_links %>%
      dplyr::distinct(input_id, gene)
    expanded_only <- expanded_links %>%
      dplyr::anti_join(primary_gene_pairs, by = c("input_id", "gene"))
    if (nrow(expanded_only) > 0) {
      expanded_agg <- expanded_only %>%
        dplyr::group_by(input_id) %>%
        dplyr::summarise(
          Expanded_Target_Genes = .collapse_genes(gene),
          .groups = "drop"
        )
    }
  }

  # --- Primary target aggregation ---
  empty_assigned <- data.frame(
    input_id = character(),
    Assigned_Target_Genes = character(),
    stringsAsFactors = FALSE
  )
  assigned_agg <- empty_assigned
  if (nrow(primary_links) > 0) {
    tmp <- primary_links %>%
      dplyr::group_by(input_id)

    if (target_priority == "promoter_then_distance") {
      tmp <- tmp %>%
        dplyr::filter(role_rank == min(role_rank, na.rm = TRUE)) %>%
        dplyr::filter(path_rank == min(path_rank, na.rm = TRUE))
    } else {
      tmp <- tmp %>%
        dplyr::filter(path_rank == min(path_rank, na.rm = TRUE)) %>%
        dplyr::filter(role_rank == min(role_rank, na.rm = TRUE))
    }

    assigned_agg <- tmp %>%
      dplyr::summarise(
        Assigned_Target_Genes = .collapse_genes(gene),
        .groups = "drop"
      )
  }

  # --- Direct P--P co-assignment (promoter_then_distance only) ---
  # For direct promoter--promoter contacts, the technical direction
  # (path0 vs path1) is determined by which endpoint the input peak
  # hits, not by biological regulatory direction.  Co-assign both
  # strict promoter endpoints when they share a single direct loop.
  has_loop_id  <- "loop_ID"  %in% colnames(primary_links)
  has_anchor_id <- "anchor_id" %in% colnames(primary_links)
  if (target_priority == "promoter_then_distance" && nrow(primary_links) > 0 &&
      has_loop_id && has_anchor_id) {
    pp_links <- primary_links %>%
      dplyr::filter(
        gene_role == "promoter",
        strict_assignment_eligible %in% TRUE,
        path_rank %in% c(0, 1),
        !is.na(loop_ID), loop_ID != "",
        !is.na(anchor_id), anchor_id != ""
      )
    if (nrow(pp_links) > 0) {
      pp_agg <- pp_links %>%
        dplyr::group_by(input_id, loop_ID) %>%
        dplyr::filter(
          any(path_rank == 0),
          any(path_rank == 1),
          dplyr::n_distinct(anchor_id) >= 2
        ) %>%
        dplyr::summarise(
          pp_genes = .collapse_genes(gene),
          .groups = "drop"
        ) %>%
        dplyr::group_by(input_id) %>%
        dplyr::summarise(
          PP_CoAssigned_Genes = .collapse_genes(
            unlist(strsplit(paste(pp_genes, collapse = ";"), ";", fixed = TRUE))
          ),
          .groups = "drop"
        )
      # Union PP co-assigned genes with existing Assigned
      assigned_agg <- assigned_agg %>%
        dplyr::full_join(pp_agg, by = "input_id") %>%
        dplyr::rowwise() %>%
        dplyr::mutate(
          Assigned_Target_Genes = {
            existing <- if (!is.na(Assigned_Target_Genes) && Assigned_Target_Genes != "") {
              strsplit(Assigned_Target_Genes, ";", fixed = TRUE)[[1]]
            } else {
              character(0)
            }
            pp <- if (!is.na(PP_CoAssigned_Genes) && PP_CoAssigned_Genes != "") {
              strsplit(PP_CoAssigned_Genes, ";", fixed = TRUE)[[1]]
            } else {
              character(0)
            }
            result <- .collapse_genes(c(existing, pp))
            if (result == "") NA_character_ else result
          }
        ) %>%
        dplyr::ungroup() %>%
        dplyr::select(-PP_CoAssigned_Genes)
    }
  }

  # --- Positional candidate aggregation ---
  empty_positional <- data.frame(
    input_id = character(),
    Candidate_Positional_Genes = character(),
    stringsAsFactors = FALSE
  )
  # Schema already normalised above (strict_assignment_eligible column
  # is guaranteed to exist and positional_candidate is always FALSE).

  positional_agg <- empty_positional
  if (nrow(loop_links) > 0) {
    pos_links <- loop_links %>%
      dplyr::filter(
        gene_role == "positional_candidate" |
          (!is.na(strict_assignment_eligible) & !strict_assignment_eligible)
      )
    # Candidate_Positional_Genes is the complement of strict assignment:
    # a gene that already has strict (promoter/gene-body) evidence for
    # the same input_id is a confirmed target, not a candidate.
    # Anti-join strict gene pairs so the two sets are mutually exclusive.
    if (nrow(pos_links) > 0 && nrow(strict_links) > 0) {
      strict_pairs <- strict_links %>%
        dplyr::distinct(input_id, gene)
      pos_links <- pos_links %>%
        dplyr::anti_join(strict_pairs, by = c("input_id", "gene"))
    }
    if (nrow(pos_links) > 0) {
      positional_agg <- pos_links %>%
        dplyr::group_by(input_id) %>%
        dplyr::summarise(
          Candidate_Positional_Genes = .collapse_genes(gene),
          .groups = "drop"
        )
    }
  }

  # --- Merge all columns back to bed_info ---
  bed_info %>%
    dplyr::select(-dplyr::any_of(c(
      "All_Loop_Connected_Genes",
      "Regulated_promoter_genes",
      "Assigned_Target_Genes",
      "Expanded_Target_Genes",
      "Candidate_Positional_Genes"
    ))) %>%
    dplyr::left_join(all_agg, by = "input_id") %>%
    dplyr::left_join(promoter_agg, by = "input_id") %>%
    dplyr::left_join(assigned_agg, by = "input_id") %>%
    dplyr::left_join(expanded_agg, by = "input_id") %>%
    dplyr::left_join(positional_agg, by = "input_id")
}

#' Internal: Rebuild Filled columns and evidence strings after chromatin recomputation
#'
#' Applies the linear nearest-gene fallback mechanism (`*_Filled` columns) and
#' recomputes `Regulated_promoter_Evidence` and
#' `Regulated_promoter_Fallback_Evidence` from the updated 3D assignments.
#'
#' @param bed_info Target annotation data frame.
#' @param current_links Expression-qualified target gene links
#'   (all structural links in chromatin-only mode; only links with
#'   \code{Passes_Expression_Filter == TRUE} in expression+chromatin mode).
#' @return The `bed_info` with `*_Filled` and evidence columns rebuilt.
#' @keywords internal
#' @noRd
.fill_fallback_targets <- function(bed_info, current_links) {
  if (is.null(bed_info)) {
    return(NULL)
  }

  # Rebuild evidence from CURRENT links only -- inactive links must not
  # contribute evidence to the current summary.
  evidence_df <- .summarise_regulated_promoter_evidence(current_links)
  bed_info <- bed_info %>%
    dplyr::select(-dplyr::any_of(c("Regulated_promoter_Evidence"))) %>%
    dplyr::left_join(evidence_df, by = "input_id")
  bed_info$Regulated_promoter_Evidence <- ifelse(
    is.na(bed_info$Regulated_promoter_Evidence) |
      bed_info$Regulated_promoter_Evidence == "",
    "none",
    bed_info$Regulated_promoter_Evidence
  )

  # Build active-only linear fallback from current_links -- this
  # guarantees Filled columns only receive genes that pass the
  # current refinement stage's expression filter.
  fallback_map <- current_links %>%
    dplyr::filter(
      source == "linear_annotation",
      !is.na(gene), gene != ""
    ) %>%
    dplyr::group_by(input_id) %>%
    dplyr::summarise(
      current_linear_fallback = .collapse_genes(gene),
      .groups = "drop"
    )
  bed_info <- bed_info %>%
    dplyr::left_join(fallback_map, by = "input_id")
  fallback_vec <- bed_info$current_linear_fallback
  ann_vec <- if ("annotation" %in% colnames(bed_info)) {
    bed_info$annotation
  } else {
    rep(NA_character_, nrow(bed_info))
  }
  fallback_evidence <- .fallback_evidence_from_annotation(ann_vec)

  bed_info %>%
    dplyr::mutate(
      All_Loop_Connected_Genes_Filled = .fill_target_gene_fallback(
        All_Loop_Connected_Genes, fallback_vec
      ),
      Regulated_promoter_genes_Filled = .fill_target_gene_fallback(
        Regulated_promoter_genes, fallback_vec
      ),
      Assigned_Target_Genes_Filled = .fill_target_gene_fallback(
        Assigned_Target_Genes, fallback_vec
      ),
      Regulated_promoter_Fallback_Evidence =
        dplyr::case_when(
          !is.na(Regulated_promoter_genes) &
            Regulated_promoter_genes != "" ~ "none",
          !is.na(Regulated_promoter_genes_Filled) &
            Regulated_promoter_genes_Filled != "" ~ fallback_evidence,
          TRUE ~ "none"
        )
    ) %>%
    dplyr::select(-dplyr::any_of("current_linear_fallback"))
}
#' Internal: Unified Target Annotation Finalizer
#'
#' Shared across basic annotation, expression refinement, and chromatin
#' refinement.  Aggregates current target gene links into strict columns
#' (policy-controlled priority), builds linear fallback columns, and marks
#' membership on the full provenance table.
#'
#' @param bed_info Target annotation data frame with \code{input_id}.
#' @param all_target_gene_links Full provenance table.
#' @param has_expression Logical.  If \code{TRUE}, only links with
#'   \code{Passes_Expression_Filter == TRUE} contribute to current summary.
#' @return A list with \code{target_annotation} and \code{target_gene_links}.
#' @keywords internal
#' @noRd
.finalize_current_target_annotation <- function(
  bed_info,
  all_target_gene_links,
  has_expression = FALSE,
  target_priority = c("promoter_then_distance", "distance_then_role"),
  max_primary_hop = 1L
) {
  # Always clear stale summary columns first: if current links are empty,
  # strict and filled columns must be NA -- never retain previous values.
  bed_info <- .clear_target_summary_columns(bed_info)

  if (.is_null_or_empty(all_target_gene_links)) {
    if (!is.null(bed_info)) {
      bed_info$Regulated_promoter_Evidence <- "none"
      bed_info$Regulated_promoter_Fallback_Evidence <- "none"
    }
    return(list(
      target_annotation = bed_info,
      target_gene_links = all_target_gene_links
    ))
  }

  current_links <- if (has_expression) {
    all_target_gene_links %>%
      dplyr::filter(Passes_Expression_Filter %in% TRUE)
  } else {
    all_target_gene_links
  }

  bed_info <- .aggregate_strict_targets(current_links, bed_info,
    target_priority = target_priority,
    max_primary_hop = max_primary_hop
  )
  bed_info <- .fill_fallback_targets(bed_info, current_links)
  all_target_gene_links <- .mark_target_gene_link_membership(
    all_target_gene_links, bed_info
  )

  if (!is.null(bed_info)) {
    all_target_gene_links$retained_after_refinement <-
      all_target_gene_links$in_regulated_promoter |
        all_target_gene_links$in_assigned_target |
        all_target_gene_links$in_all_loop_connected |
        all_target_gene_links$in_regulated_promoter_filled |
        all_target_gene_links$in_assigned_target_filled
  }

  list(
    target_annotation = bed_info,
    target_gene_links = all_target_gene_links
  )
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
                                                 reclass_map,
                                                 tss_assignment_provenance = NULL) {
  anchor_state <- attr(original_res, "looplook_anchor_state")
  if (is.null(anchor_state)) {
    stop(
      "Cannot recompute chromatin-aware target links because ",
      "`looplook_anchor_state` is missing from the input object. ",
      "Re-run `annotate_peaks_and_loops()` or set ",
      "`recompute_targets = FALSE` to preserve upstream target tables.",
      call. = FALSE
    )
  }

  map_info <- anchor_state$map_info
  gr_anchors <- anchor_state$gr_anchors
  ego_list_target <- anchor_state$ego_list_target
  g <- anchor_state$g
  neighbor_hop <- anchor_state$neighbor_hop

  if (.is_null_or_empty(map_info)) {
    return(list(target_annotation = NULL, target_gene_links = NULL))
  }

  # --- 1. Update map_info type_code for reclassified anchors ---
  type_updates <- setNames(reclass_map$new_type, reclass_map$anchor_id)
  matched <- match(map_info$anchor_id, names(type_updates))
  update_idx <- which(!is.na(matched) & reclass_map$changed[matched])
  if (length(update_idx) > 0) {
    map_info$type_code[update_idx] <-
      type_updates[map_info$anchor_id[update_idx]]
  }

  # --- 2. Update map_info type_code ---
  map_info$SYMBOL <- trimws(map_info$SYMBOL)

  # --- 3. Build new target_gene_links ---
  hit_df <- anchor_state$target_hit_df
  if (.is_null_or_empty(hit_df)) {
    warning("No stored peak-to-anchor hit_df found; ",
      "target_gene_links will be linear-only.",
      call. = FALSE
    )
    hit_df <- data.frame(
      qid = integer(0), anchor_id = character(0),
      stringsAsFactors = FALSE
    )
  }
  new_links <- .build_target_gene_links(
    hit_df = hit_df,
    bed_info = original_res$target_annotation,
    loop_annotation_final = annotation_res$loop_annotation,
    map_info = map_info,
    ego_list_target = ego_list_target,
    g = anchor_state$g,
    neighbor_hop = if (is.null(neighbor_hop)) 0L else neighbor_hop
  )

  # --- 3.5 Inherit expression provenance ---
  # Chromatin recomputation rebuilds target_gene_links from scratch,
  # which drops expression columns added by prior expression refinement.
  # Recovery is done in three steps:
  #   A. Pair-level join (input_id + gene) from old links
  #   B. Gene-level join from old links (for new E->P gene assignments)
  #   C. Direct lookup against anchor_state$expr_vals (for genes never
  #      appearing in any pre-chromatin link)
  expr_cols <- c("Mean_Expression", "Passes_Expression_Filter", "Expression_State")
  old_links <- original_res$target_gene_links
  has_old <- !is.null(old_links) && nrow(old_links) > 0
  expr_cols_present <- if (has_old) intersect(expr_cols, colnames(old_links)) else character(0)

  if (has_old && length(expr_cols_present) > 0) {
    # Step A: pair-level join
    expr_map <- old_links %>%
      dplyr::select(dplyr::all_of(c("input_id", "gene", expr_cols_present))) %>%
      dplyr::distinct()
    expr_map <- expr_map[!duplicated(expr_map[, c("input_id", "gene")]), ]
    new_links <- dplyr::left_join(new_links, expr_map, by = c("input_id", "gene"))

    # Step B: gene-level join for pairs not in old links
    gene_expr <- old_links %>%
      dplyr::select(dplyr::all_of(c("gene", expr_cols_present))) %>%
      dplyr::distinct()
    gene_expr <- gene_expr[!duplicated(gene_expr$gene), ]
    new_links <- dplyr::left_join(new_links, gene_expr,
      by = "gene",
      suffix = c("", ".gene")
    )
    for (col in expr_cols_present) {
      gene_col <- paste0(col, ".gene")
      if (gene_col %in% colnames(new_links)) {
        new_links[[col]] <- dplyr::coalesce(new_links[[col]], new_links[[gene_col]])
        new_links[[gene_col]] <- NULL
      }
    }
  }

  # Step C: direct lookup from anchor_state expression vector.
  # Runs regardless of whether old_links exists -- handles the case
  # where expression refinement ran but produced zero target links.
  if (nrow(new_links) > 0 && !is.null(anchor_state$expr_vals)) {
    # Ensure expression columns exist before checking for NAs.
    # When old_links is empty and steps A/B don't run, these
    # columns are absent; is.na(NULL) is logical(0), skipping
    # the entire block silently.
    if (!"Mean_Expression" %in% colnames(new_links)) {
      new_links$Mean_Expression <- NA_real_
    }
    if (!"Passes_Expression_Filter" %in% colnames(new_links)) {
      new_links$Passes_Expression_Filter <- NA
    }
    vals <- anchor_state$expr_vals
    thr <- anchor_state$expr_threshold
    still_missing <- is.na(new_links$Mean_Expression)
    if (any(still_missing)) {
      missing_genes <- unique(new_links$gene[still_missing])
      lookup <- data.frame(
        gene = missing_genes,
        Mean_Expression = as.numeric(vals[toupper(missing_genes)]),
        stringsAsFactors = FALSE
      )
      if (!is.null(thr)) {
        lookup$Passes_Expression_Filter <- .passes_expression_threshold(
          lookup$Mean_Expression, thr
        )
      } else {
        lookup$Passes_Expression_Filter <- NA
      }
      new_links <- dplyr::left_join(new_links, lookup,
        by = "gene",
        suffix = c("", ".direct")
      )
      for (col in expr_cols) {
        direct_col <- paste0(col, ".direct")
        if (direct_col %in% colnames(new_links)) {
          new_links[[col]] <- dplyr::coalesce(new_links[[col]], new_links[[direct_col]])
          new_links[[direct_col]] <- NULL
        }
      }
    }
  }

  # Step D: compute Expression_State for rows still missing it.
  # Distinguishes unmeasured (gene not in expression data) from
  # measured_silent (gene measured but below threshold).
  if (nrow(new_links) > 0 && !is.null(anchor_state$expr_vals)) {
    if (!"Expression_State" %in% colnames(new_links)) {
      new_links$Expression_State <- NA_character_
    }
    need_state <- is.na(new_links$Expression_State) |
      new_links$Expression_State == ""
    if (any(need_state)) {
      vals <- anchor_state$expr_vals
      has_threshold <- !is.null(anchor_state$expr_threshold)
      genes_need <- new_links$gene[need_state]
      keys <- toupper(genes_need)
      measured <- keys %in% toupper(names(vals))
      mean_expr <- as.numeric(unname(vals[keys]))
      if (has_threshold) {
        thr <- anchor_state$expr_threshold
        passes <- .passes_expression_threshold(mean_expr, thr)
        new_links$Expression_State[need_state] <- dplyr::case_when(
          !measured ~ "unmeasured",
          measured & passes ~ "active",
          measured & !passes ~ "measured_silent",
          TRUE ~ "unmeasured"
        )
      } else {
        # Threshold not set: cannot determine active vs silent.
        # All measured genes are "measured_not_assessed".
        new_links$Expression_State[need_state] <- dplyr::case_when(
          !measured ~ "unmeasured",
          measured ~ "measured_not_assessed",
          TRUE ~ "unmeasured"
        )
      }
      # Force-rewrite Passes_Expression_Filter for all need_state rows.
      # This guarantees the invariant:
      #   active ->TRUE, measured_silent ->FALSE, all others ->NA
      # Step C may have written FALSE for unmeasured genes
      # (Mean_Expression=NA), which must be corrected here.
      new_links$Passes_Expression_Filter[need_state] <- dplyr::case_when(
        new_links$Expression_State[need_state] == "active" ~ TRUE,
        new_links$Expression_State[need_state] == "measured_silent" ~ FALSE,
        TRUE ~ NA
      )
    }
  }

  # --- 4. Annotate with chromatin provenance ---
  new_links <- .chromatin_annotate_links(new_links, reclass_map, annotation_res,
    tss_assignment_provenance = tss_assignment_provenance
  )

  # --- 5. Rebuild target_annotation from chromatin-aware links ---
  # Correct order: rebuild strict ->filter by expression ->build Filled
  # ->mark membership.  This ensures Filled uses expression-filtered
  # strict columns as its base, so that when strict is filtered empty
  # the linear fallback can fill active genes.
  has_expr <- !is.null(anchor_state$expr_vals) &&
    !is.null(anchor_state$expr_threshold)
  wl <- if (has_expr) {
    names(anchor_state$expr_vals)[
      .passes_expression_threshold(
        anchor_state$expr_vals, anchor_state$expr_threshold
      )
    ]
  } else {
    NULL
  }

  new_bed_info <- .clear_target_summary_columns(
    original_res$target_annotation
  )

  # Unified target finalizer: when expression data exists, only links
  # passing the expression filter participate in current target columns;
  # a silent promoter must not block an active gene-body link.
  up_target_priority <- if (!is.null(original_res$metadata$parameters$target_priority)) {
    original_res$metadata$parameters$target_priority
  } else {
    "promoter_then_distance"
  }
  up_max_primary_hop <- if (!is.null(original_res$metadata$parameters$max_primary_hop)) {
    original_res$metadata$parameters$max_primary_hop
  } else {
    1L
  }
  final <- .finalize_current_target_annotation(
    bed_info = new_bed_info,
    all_target_gene_links = new_links,
    has_expression = has_expr,
    target_priority = up_target_priority,
    max_primary_hop = up_max_primary_hop
  )

   # --- End-to-end QC: verify chromatin reclassification propagates to output ---
  n_reclassified <- sum(reclass_map$changed, na.rm = TRUE)
  # Determine which reclassified anchors are actually connected to targets
  changed_ids <- reclass_map$anchor_id[reclass_map$changed %in% TRUE]
  target_anchor_ids <- unique(new_links$anchor_id[
    !is.na(new_links$anchor_id) & new_links$anchor_id != "" &
    !is.na(new_links$gene_role) & new_links$gene_role %in%
      c("promoter", "gene_body", "enhancer_candidate")
  ])
  target_changed_ids <- intersect(changed_ids, target_anchor_ids)
  if (n_reclassified > 0 && length(target_changed_ids) > 0) {
    old_bed <- original_res$target_annotation
    new_bed <- final$target_annotation
    if (!is.null(old_bed) && !is.null(new_bed) &&
      "input_id" %in% colnames(old_bed) && "input_id" %in% colnames(new_bed)) {
      # Canonicalise gene strings for robust set comparison
      .canon_gene <- function(x) {
        if (is.na(x) || !nzchar(trimws(x))) return(NA_character_)
        g <- trimws(unlist(strsplit(x, ";", fixed = TRUE)))
        g <- g[!is.na(g) & nzchar(g)]
        if (length(g) == 0) return(NA_character_)
        .collapse_genes(g)
      }
      old_cmp <- data.frame(
        input_id = old_bed$input_id,
        old_reg = vapply(if ("Regulated_promoter_genes" %in% colnames(old_bed))
          old_bed$Regulated_promoter_genes else rep(NA_character_, nrow(old_bed)),
          .canon_gene, character(1), USE.NAMES = FALSE),
        old_assign = vapply(if ("Assigned_Target_Genes" %in% colnames(old_bed))
          old_bed$Assigned_Target_Genes else rep(NA_character_, nrow(old_bed)),
          .canon_gene, character(1), USE.NAMES = FALSE),
        stringsAsFactors = FALSE
      )
      new_cmp <- data.frame(
        input_id = new_bed$input_id,
        new_reg = vapply(if ("Regulated_promoter_genes" %in% colnames(new_bed))
          new_bed$Regulated_promoter_genes else rep(NA_character_, nrow(new_bed)),
          .canon_gene, character(1), USE.NAMES = FALSE),
        new_assign = vapply(if ("Assigned_Target_Genes" %in% colnames(new_bed))
          new_bed$Assigned_Target_Genes else rep(NA_character_, nrow(new_bed)),
          .canon_gene, character(1), USE.NAMES = FALSE),
        stringsAsFactors = FALSE
      )
      cmp <- merge(old_cmp, new_cmp, by = "input_id", all = TRUE, sort = FALSE)
      # Treat NA and empty as equivalent null
      .eq <- function(x, y) {
        xn <- is.na(x) | !nzchar(x)
        yn <- is.na(y) | !nzchar(y)
        (xn & yn) | (!xn & !yn & x == y)
      }
      reg_changed <- sum(!.eq(cmp$old_reg, cmp$new_reg), na.rm = TRUE)
      assign_changed <- sum(!.eq(cmp$old_assign, cmp$new_assign), na.rm = TRUE)

      if (reg_changed == 0 && assign_changed == 0) {
        warning(
          "Chromatin reclassification affected ", length(target_changed_ids),
          " target-connected anchors, but no Regulated_promoter_genes or ",
          "Assigned_Target_Genes row changed. ",
          "This can be valid when roles remain equivalent, genes are retained ",
          "through alternative links, or P--P union preserves the final set. ",
          "Inspect link-level provenance if unexpected.",
          call. = FALSE
        )
      } else {
        message(
          "  chromatin target QC: ", reg_changed, " rows changed in ",
          "Regulated_promoter_genes, ", assign_changed,
          " rows changed in Assigned_Target_Genes (",
          n_reclassified, " total reclassified, ",
          length(target_changed_ids), " target-connected)"
        )
      }
    }
  }

  list(
    target_annotation = .relocate_target_annotation_columns(final$target_annotation),
    target_gene_links = final$target_gene_links
  )
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
  # Use canonical Linked_Loop_IDs produced by main annotation via
  # cascade-filtered peak--anchor overlaps. This avoids re-computing
  # overlaps with missing thresholds or mismatched seqlevels.
  if (!"Linked_Loop_IDs" %in% colnames(bed_info)) {
    return(NULL)
  }
  accepted_ids <- clean_gene_names(bed_info$Linked_Loop_IDs, "[;,]")
  if (length(accepted_ids) == 0) {
    return(NULL)
  }
  tgt_loops <- loop_df %>%
    dplyr::filter(loop_ID %in% accepted_ids)
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
    message("Target Sankey skipped: no target BED data. Provide a target_bed to annotate_peaks_and_loops.")
    return(NULL)
  }
  if (!"SANKEY_RAW_GENES" %in% colnames(bed_info)) {
    message("Target Sankey skipped: missing SANKEY_RAW_GENES column.")
    return(NULL)
  }
  if (!requireNamespace("networkD3", quietly = TRUE) || !requireNamespace("htmlwidgets", quietly = TRUE)) {
    message("Target Sankey skipped: install 'networkD3' and 'htmlwidgets' packages.")
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
    if (any(toupper(gs) %in% toupper(wl))) {
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
    message("Target Sankey skipped: no rows with valid L1/L2/L3 classification.")
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
    message("Target Sankey skipped: all rows filtered after label mapping.")
    return(NULL)
  }

  nodes <- unique(c(
    unlist(sankey_df$L1_Label, use.names = FALSE),
    unlist(sankey_df$L2_Label, use.names = FALSE),
    unlist(sankey_df$L3_Label, use.names = FALSE)
  ))
  if (length(nodes) < 2) {
    message("Target Sankey skipped: fewer than 2 unique nodes after filtering.")
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
  color_scale <- paste0(
    'd3.scaleOrdinal().range(["',
    paste(sankey_colors, collapse = '","'), '"])'
  )
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
  sn <- htmlwidgets::onRender(sn, '
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
          gradient.append("stop").attr("offset", "0%")
            .attr("stop-color", sourceColor);
          gradient.append("stop").attr("offset", "100%")
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
    ')
  sn
}

#' Internal: Build Refinement Karyotype Heatmaps
#'
#' Generates \code{Refined_Karyo_CurrentTargets} and \code{Refined_Karyo_TargetGenes}
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
  txdb_pkg <- species_txdb_pkg(species)
  org_db <- species_orgdb_pkg(species)
  if (is.null(txdb_pkg) || is.null(org_db) ||
    !requireNamespace(txdb_pkg, quietly = TRUE) ||
    !requireNamespace(org_db, quietly = TRUE)) {
    message("Refinement karyotype plots skipped: TxDb/OrgDb packages not available.")
    return(plot_list)
  }
  txdb_obj <- utils::getFromNamespace(txdb_pkg, txdb_pkg)

  # Determine which gene sets need karyotypes
  need_current <- "Putative_Target_Genes" %in% colnames(loop_df)
  need_target <- !is.null(bed_info) &&
    "Assigned_Target_Genes_Filled" %in% colnames(bed_info)
  if (!need_current && !need_target) {
    return(plot_list)
  }

  g_current <- if (need_current) {
    clean_gene_names(loop_df$Putative_Target_Genes, ";")
  } else {
    character(0)
  }
  g_target <- if (need_target) {
    clean_gene_names(bed_info$Assigned_Target_Genes_Filled, ";")
  } else {
    character(0)
  }

  g_all <- unique(c(g_current, g_target))
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

  if (length(g_current) > 0) {
    target_genes_gr <- all_genes_gr[
      S4Vectors::mcols(all_genes_gr)$SYMBOL %in% g_current
    ]
    plot_list$Refined_Karyo_CurrentTargets <- draw_karyo_heatmap_internal(
      target_genes_gr,
      "Current Refined Target Genes", karyo_bin_size, 0.99, txdb_obj, species,
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
#' @param subtitle Character or NULL. Custom subtitle. NULL uses the default
#'   expression-refinement subtitle. Default: \code{NULL}.
#' @return A \code{ggplot} object.
#' Internal: Build Loop-Type UpSet Plot
#' @keywords internal
#'
#' Replaces the rose/coxcomb chart with an UpSet-style combinatorial view.
#' Each loop type (e.g. "E-P", "dual-E") is an intersection of two anchor-
#' type sets.  The top panel shows counts on a log10 axis to compress wide
#' dynamic range; the bottom panel shows which anchor types participate in
#' each combination via a dot-and-line matrix.  This format is standard in
#' 3D-genomics papers and highlights the combinatorial nature of chromatin-
#' aware anchor classification.
#'
#' @param loop_df Refined loop annotation data frame with a \code{loop_type} column.
#' @param project_name Character. Plot title prefix.
#' @param subtitle Character or \code{NULL}. Plot subtitle.
#' @return A \code{ggplot} object assembled with \code{patchwork}.
#' Internal: Build Rose Plot (Expression Refinement)
#' @keywords internal
#' @noRd
#' @noRd
.build_rose_plot <- function(loop_df, custom_colors, project_name, subtitle = NULL) {
  rose_data <- loop_df %>%
    dplyr::group_by(loop_type) %>%
    dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(
      fraction = count / sum(count),
      legend_label = paste0(
        loop_type, " (n=", count, ", ",
        round(fraction * 100, 1), "%)"
      ),
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
    ggplot2::labs(
      title = paste0(project_name, ": Structural Loop Types After Refinement"),
      subtitle = if (is.null(subtitle)) "Full refined network. For expression-supported subset, filter on Retained_In_Functional_Network." else subtitle
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = "right",
      legend.text = ggplot2::element_text(size = 10)
    )
}

#' Internal: Build Loop-Type UpSet Plot (Chromatin Refinement)
#' @keywords internal
#' @noRd
.build_loop_type_upset <- function(loop_df, project_name, subtitle = NULL) {
  # --- per-combination counts ---
  type_counts <- table(loop_df$loop_type)
  type_counts <- sort(type_counts, decreasing = TRUE)
  comb_names <- names(type_counts)
  n_comb <- length(comb_names)
  if (n_comb == 0) {
    return(NULL)
  }

  # --- binary matrix: rows = combinations, cols = anchor types ---
  anchor_types <- c("E", "P", "G", "eP", "eG", "dual")
  # Keep only types that actually appear in the data
  all_parts <- unique(unlist(strsplit(comb_names, "-")))
  anchor_types <- intersect(anchor_types, all_parts)
  n_types <- length(anchor_types)
  if (n_types < 2) {
    return(NULL)
  }

  bin_mat <- matrix(0L,
    nrow = n_comb, ncol = n_types,
    dimnames = list(comb_names, anchor_types)
  )
  for (i in seq_len(n_comb)) {
    parts <- strsplit(comb_names[i], "-")[[1]]
    bin_mat[i, intersect(parts, anchor_types)] <- 1L
  }

  # --- colours: Wong palette, consistent with Sankey ---
  type_pal <- c(
    "E" = "#E69F00", "P" = "#0072B2", "G" = "#CC79A7",
    "eP" = "#009E73", "eG" = "#56B4E9", "dual" = "#CC0000"
  )
  type_cols <- type_pal[anchor_types]

  # --- top panel: bar chart with log10 y ---
  bar_df <- data.frame(
    combination = factor(comb_names, levels = comb_names),
    count = as.integer(type_counts),
    stringsAsFactors = FALSE
  )
  # Guard against log10(0) -- should never happen after table()
  bar_df$log_count <- log10(pmax(bar_df$count, 1))

  y_breaks <- pretty(bar_df$log_count, n = 4)
  y_labels <- vapply(y_breaks, function(b) {
    scales::comma(10^b, accuracy = 1)
  }, character(1))

  p_bar <- ggplot2::ggplot(bar_df, ggplot2::aes(x = combination, y = log_count)) +
    ggplot2::geom_col(fill = "#555555", width = 0.7) +
    ggplot2::scale_y_continuous(
      name = "Intersection size (log10)",
      breaks = y_breaks, labels = y_labels, expand = c(0, 0.05)
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(color = "gray92", linewidth = 0.25),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10),
      plot.margin = ggplot2::margin(b = 0, t = 5, l = 5, r = 5)
    ) +
    ggplot2::labs(
      title = paste0(project_name, ": Loop Type Combinations"),
      subtitle = if (is.null(subtitle)) "" else subtitle
    )

  # --- bottom panel: dot matrix ---
  dot_df <- expand.grid(
    combination = factor(comb_names, levels = comb_names),
    anchor_type = factor(anchor_types, levels = rev(anchor_types)),
    stringsAsFactors = FALSE
  )
  dot_df$present <- as.logical(bin_mat[cbind(
    as.character(dot_df$combination),
    as.character(dot_df$anchor_type)
  )])

  p_dot <- ggplot2::ggplot(dot_df, ggplot2::aes(x = combination, y = anchor_type)) +
    # Absent dots (small, grey, back-most layer)
    ggplot2::geom_point(
      data = subset(dot_df, !present),
      color = "#E8E8E8", size = 1.0
    ) +
    # Connection lines -- above grey dots, behind coloured dots
    ggplot2::geom_line(
      data = subset(dot_df, present),
      ggplot2::aes(group = combination), color = "#999999", linewidth = 0.5
    ) +
    # Present dots (coloured, on top of lines)
    ggplot2::geom_point(
      data = subset(dot_df, present),
      ggplot2::aes(color = anchor_type), size = 3.2
    ) +
    # Self-loop arc -- same colour as lines
    {
      if (requireNamespace("ggforce", quietly = TRUE)) {
        ggforce::geom_arc(
          data = subset(dot_df, present & grepl("^([A-Za-z]+)-\\1$", as.character(combination))),
          ggplot2::aes(
            x0 = as.numeric(combination), y0 = as.numeric(anchor_type),
            r = 0.22, start = pi, end = 2 * pi
          ),
          color = "#999999", linewidth = 0.5, inherit.aes = FALSE
        )
      }
    } +
    ggplot2::scale_color_manual(values = type_cols, guide = "none") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 9),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(t = 0, b = 5, l = 5, r = 5)
    )

  # --- side panel: set-size bar (loops per anchor type, not combinations) ---
  set_counts <- colSums(bin_mat * bar_df$count)
  set_df <- data.frame(
    anchor_type = factor(names(set_counts), levels = rev(anchor_types)),
    count = as.integer(set_counts),
    stringsAsFactors = FALSE
  )
  set_df$log_count <- log10(pmax(set_df$count, 1))
  side_breaks <- pretty(c(0, max(set_df$log_count)), n = 3)

  p_side <- ggplot2::ggplot(set_df, ggplot2::aes(x = log_count, y = anchor_type)) +
    ggplot2::geom_col(ggplot2::aes(fill = anchor_type), width = 0.7) +
    ggplot2::scale_fill_manual(values = type_cols, guide = "none") +
    ggplot2::scale_x_reverse(
      name = "Set size (log10)",
      breaks = side_breaks,
      labels = scales::comma(10^side_breaks, accuracy = 1),
      expand = c(0, 0.05)
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      axis.title.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(color = "gray92", linewidth = 0.25),
      panel.grid.major.y = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(l = 5, r = 5, t = 5, b = 5)
    )

  # --- assemble: three-panel UpSet via explicit design ---
  # Layout:  A B
  #          A C
  # where A = set size (left), B = intersection bar (top-right),
  # C = combination matrix (bottom-right).
  if (requireNamespace("patchwork", quietly = TRUE)) {
    upset_design <- "AB\nAC"
    p_combined <- p_side + p_bar + p_dot +
      patchwork::plot_layout(
        design = upset_design,
        widths = c(1, 6),
        heights = c(3, 2)
      )
    return(p_combined)
  }
  # fallback: bar only (no side panel or dot matrix)
  return(p_bar +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 9)))
}
#' Internal: Build Refinement Visualization Suite
#'
#' Generates diagnostic plots for expression-aware refinement: dumbbell
#' comparison, donut, Sankey, karyotype heatmaps, and rose plot.
#'
#' @param original_loop_df Loop annotation before refinement.
#' @param loop_df Loop annotation after refinement.
#' Internal: Build Chromatin Refinement Dumbbell Plot
#'
#' Compares anchor-type counts before vs. after chromatin-aware
#' reclassification. Dark green academic palette for publication use.
#'
#' @param reclass_map Reclassification data frame from .chromatin_reclassify.
#' @param project_name Character. Project prefix for the plot title.
#' @return A ggplot object.
#' @keywords internal
#' @noRd
.build_chromatin_dumbbell_plot <- function(reclass_map, project_name) {
  df_before <- reclass_map %>%
    dplyr::group_by(old_type) %>%
    dplyr::summarise(Original = dplyr::n(), .groups = "drop") %>%
    dplyr::rename(type = old_type)
  df_after <- reclass_map %>%
    dplyr::group_by(new_type) %>%
    dplyr::summarise(Refined = dplyr::n(), .groups = "drop") %>%
    dplyr::rename(type = new_type)
  df_dumbbell <- dplyr::full_join(df_before, df_after, by = "type") %>%
    dplyr::mutate(
      Original = ifelse(is.na(Original), 0L, Original),
      Refined  = ifelse(is.na(Refined), 0L, Refined)
    ) %>%
    dplyr::arrange(dplyr::desc(Refined))
  df_dumbbell$type <- factor(df_dumbbell$type, levels = rev(df_dumbbell$type))
  df_long <- df_dumbbell %>%
    tidyr::pivot_longer(
      cols = c("Original", "Refined"),
      names_to = "Source", values_to = "Count"
    )

  green_colors <- c("Original" = "#999999", "Refined" = "#2E8B57")

  ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = df_dumbbell,
      ggplot2::aes(y = type, yend = type, x = Original, xend = Refined),
      color = "#b2b2b2", linewidth = 0.8
    ) +
    ggplot2::geom_point(
      data = df_long,
      ggplot2::aes(x = Count, y = type, color = Source), size = 3.5
    ) +
    ggplot2::scale_color_manual(values = green_colors) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::labs(
      title = paste0(project_name, ": Anchor Refinement by Chromatin"),
      subtitle = paste0(
        sum(reclass_map$changed), " / ", nrow(reclass_map),
        " anchors reclassified"
      ),
      x = "Number of Anchors", y = "Anchor Type"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 11, color = "#444444"),
      legend.position = "top",
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )
}

#' Internal: Build Chromatin Refinement Sankey Flow Plot
#'
#' Shows anchor-type transitions from old to new as a left-to-right flow
#' diagram. Uses ggforce::geom_parallel_sets for publication-quality rendering.
#'
#' @param reclass_map Reclassification data frame from .chromatin_reclassify.
#' @param project_name Character. Project prefix for the plot title.
#' @return A ggplot object.
#' @keywords internal
#' @noRd
.build_chromatin_sankey_plot <- function(reclass_map, project_name, sankey_colors = NULL) {
  if (.is_null_or_empty(reclass_map)) {
    return(NULL)
  }
  if (!requireNamespace("networkD3", quietly = TRUE) ||
    !requireNamespace("htmlwidgets", quietly = TRUE)) {
    message("Chromatin Sankey skipped: install 'networkD3' and 'htmlwidgets' packages.")
    return(NULL)
  }

  # Build transition counts: old_type -> new_type
  flow_df <- reclass_map %>%
    dplyr::filter(!is.na(old_type), !is.na(new_type)) %>%
    dplyr::group_by(old_type, new_type) %>%
    dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
    dplyr::filter(count > 0)

  if (nrow(flow_df) == 0) {
    message("Chromatin Sankey skipped: no reclassification flows to display.")
    return(NULL)
  }

  # Compute per-type totals for labels
  old_totals <- flow_df %>%
    dplyr::group_by(old_type) %>%
    dplyr::summarise(total = sum(count), .groups = "drop")
  new_totals <- flow_df %>%
    dplyr::group_by(new_type) %>%
    dplyr::summarise(total = sum(count), .groups = "drop")

  n_total <- nrow(reclass_map)
  # Build nodes: type name + count + percentage (column headers via D3)
  old_nodes <- paste0(
    old_totals$old_type,
    " (n=", format(old_totals$total, big.mark = ","),
    ", ", round(old_totals$total / n_total * 100, 1), "%)"
  )
  new_nodes <- paste0(
    new_totals$new_type,
    " (n=", format(new_totals$total, big.mark = ","),
    ", ", round(new_totals$total / n_total * 100, 1), "%)"
  )
  nodes <- data.frame(
    name = c(old_nodes, new_nodes),
    stringsAsFactors = FALSE
  )

  # Build links: source -> target with value = count
  old_lookup <- setNames(
    seq_along(old_nodes) - 1L,
    old_totals$old_type
  )
  new_lookup <- setNames(
    seq_along(new_nodes) - 1L +
      length(old_nodes),
    new_totals$new_type
  )

  links <- flow_df %>%
    dplyr::transmute(
      source = old_lookup[as.character(old_type)],
      target = new_lookup[as.character(new_type)],
      value  = count
    )

  n_changed <- sum(reclass_map$changed)

  # Default: Wong colourblind-safe palette with type-to-colour mapping.
  # Users can override specific types via sankey_colors (e.g. sankey_colors = c("dual" = "#8B0000")).
  default_palette <- c(
    "eP"   = "#009E73", # bluish-green
    "eG"   = "#56B4E9", # sky blue
    "P"    = "#0072B2", # blue
    "G"    = "#CC79A7", # reddish-purple
    "E"    = "#E69F00", # orange
    "dual" = "#CC0000" # red -- clearly distinct from E orange #E69F00
  )
  if (!is.null(sankey_colors)) {
    if (!is.character(sankey_colors) || is.null(names(sankey_colors))) {
      warning("sankey_colors must be a named character vector; using default palette.",
        call. = FALSE
      )
    } else {
      valid <- intersect(names(sankey_colors), names(default_palette))
      if (length(valid) > 0) default_palette[valid] <- sankey_colors[valid]
    }
  }

  all_types <- unique(c(
    as.character(old_totals$old_type),
    as.character(new_totals$new_type)
  ))
  type_colors <- default_palette[all_types]
  # Fallback for any type not in the predefined palette
  missing_types <- all_types[is.na(type_colors)]
  if (length(missing_types) > 0) {
    fallback <- setNames(get_colors(length(missing_types), "Set2"), missing_types)
    type_colors[is.na(type_colors)] <- fallback[missing_types]
  }
  names(type_colors) <- all_types
  # Per-node colours: extract type from node name ("eP (n=6910)" -> "eP")
  # then look up the type's colour.  Each node has a unique name (count
  # in parentheses), so D3 sees N distinct domain values and assigns the
  # N colours in order.  Variable renamed from sankey_colors to avoid
  # overwriting the user-supplied parameter of the same name.
  type_json <- jsonlite::toJSON(as.list(type_colors), auto_unbox = TRUE)

  sn <- networkD3::sankeyNetwork(
    Links = links, Nodes = nodes,
    Source = "source", Target = "target",
    Value = "value", NodeID = "name",
    units = "", fontSize = 12, nodeWidth = 30,
    iterations = 0
  )

  sn$sizingPolicy$defaultWidth <- "100%"
  sn$sizingPolicy$defaultHeight <- "450px"

  # Node colours via onRender -- bypass networkD3's internal node reordering
  # by applying the type->colour map directly.
  sn <- htmlwidgets::onRender(sn, sprintf('function(el, x) {
  var typeColors = %s;
  var svg = d3.select(el).select("svg");
  function nodeType(name) {
    if (!name) return "";
    return name.replace(/ \\(n=.*$/, "");
  }
  svg.selectAll(".node").each(function(d) {
    var t = nodeType(d.name);
    if (t && typeColors[t]) {
      d3.select(this).select("rect").style("fill", typeColors[t]);
    }
  });
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
          // Before / After column labels
          var nodes = [];
          svg.selectAll(".node").each(function(d) { nodes.push(d); });
          if (nodes.length > 0) {
            var leftX  = d3.min(nodes, function(d) { return d.x; });
            var rightX = d3.max(nodes, function(d) { return d.x + d.dx; });
            var nodeW  = nodes[0].dx;
            var topY   = d3.min(nodes, function(d) { return d.y; }) - 20;
            svg.append("text")
              .attr("x", leftX + nodeW / 2)
              .attr("y", topY)
              .attr("text-anchor", "middle")
              .style("font-size", "12px")
              .style("font-weight", "bold")
              .style("fill", "#555555")
              .text("Before");
            svg.append("text")
              .attr("x", rightX - nodeW / 2)
              .attr("y", topY)
              .attr("text-anchor", "middle")
              .style("font-size", "12px")
              .style("font-weight", "bold")
              .style("fill", "#555555")
              .text("After");
          }
        }
        ', type_json))
  sn
}

#' Internal: Build Chromatin Mark Enrichment Heatmap
#'
#' Aggregated ComplexHeatmap showing the percentage of anchors in each
#' reclassification group that overlap each chromatin mark. One row per
#' reclassification type, one column per mark. Cell colour encodes %%
#' positive (white = 0%%, deep green = 100%%), with numeric labels.
#' Left annotation: horizontal bar showing the number of anchors (N)
#' per reclassification group.
#'
#' @param validation Chromatin validation data frame with mark columns.
#' @param reclass_map Reclassification data frame from .chromatin_reclassify.
#' @param project_name Character. Project prefix for the plot title.
#' @return A ComplexHeatmap grob, or NULL if unavailable.
#' @keywords internal
#' @noRd
.build_chromatin_mark_heatmap <- function(validation, reclass_map, project_name) {
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE) ||
    !requireNamespace("circlize", quietly = TRUE)) {
    message("Chromatin MarkHeatmap skipped: install 'ComplexHeatmap' and 'circlize' packages.")
    return(NULL)
  }
  if (.is_null_or_empty(validation)) {
    return(NULL)
  }
  if (.is_null_or_empty(reclass_map)) {
    return(NULL)
  }

  changed_ids <- reclass_map$anchor_id[reclass_map$changed]
  if (length(changed_ids) == 0) {
    message("Chromatin MarkHeatmap skipped: no anchors were reclassified.")
    return(NULL)
  }

  mark_cols <- c("H3K4me1", "H3K27ac", "ATAC", "H3K4me3", "H3K27me3")
  available_marks <- intersect(mark_cols, colnames(validation))
  if (length(available_marks) == 0) {
    message("Chromatin MarkHeatmap skipped: no mark columns found in validation data.")
    return(NULL)
  }

  df <- validation %>%
    dplyr::filter(anchor_id %in% changed_ids) %>%
    dplyr::left_join(
      reclass_map %>% dplyr::select(anchor_id, old_type, new_type),
      by = "anchor_id"
    ) %>%
    dplyr::mutate(Reclassification = paste0(old_type, " -> ", new_type))

  # Compute % positive per mark x reclassification group
  reclass_counts <- df %>%
    dplyr::count(Reclassification, name = "N") %>%
    dplyr::arrange(dplyr::desc(N))

  pct_mat <- matrix(NA_real_,
    nrow = nrow(reclass_counts),
    ncol = length(available_marks),
    dimnames = list(reclass_counts$Reclassification, available_marks)
  )

  for (rc in reclass_counts$Reclassification) {
    sub <- df[df$Reclassification == rc, , drop = FALSE]
    for (mk in available_marks) {
      vals <- sub[[mk]]
      pct_mat[rc, mk] <- mean(!is.na(vals) & vals, na.rm = TRUE) * 100
    }
  }

  row_labs <- reclass_counts$Reclassification

  # Green gradient: light 0%, mid 50%, deep 100%
  col_fun <- circlize::colorRamp2(
    c(0, 50, 100),
    c("#F7FCF5", "#74C476", "#00441B")
  )

  # Cell labels: compact, only show non-zero
  cell_fun <- function(j, i, x, y, width, height, fill) {
    val <- pct_mat[i, j]
    if (!is.na(val) && val > 0) {
      grid::grid.text(sprintf("%.0f%%", val), x, y,
        gp = grid::gpar(
          fontsize = 7, fontface = "plain",
          col = if (val > 60) "white" else "#333333"
        )
      )
    }
  }

  # --- Left annotation: N per reclassification type ---
  n_counts <- reclass_counts$N
  names(n_counts) <- reclass_counts$Reclassification

  ha_left <- ComplexHeatmap::rowAnnotation(
    N = ComplexHeatmap::anno_barplot(
      n_counts,
      gp = grid::gpar(fill = "#8BA398", col = NA),
      bar_width = 0.7,
      axis = TRUE,
      axis_param = list(gp = grid::gpar(fontsize = 7))
    ),
    annotation_label = "N",
    annotation_name_gp = grid::gpar(fontsize = 9, fontface = "bold"),
    show_legend = FALSE
  )

  n_changed <- sum(reclass_map$changed)

  ComplexHeatmap::Heatmap(
    pct_mat,
    name = "%",
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_names_side = "left",
    row_labels = row_labs,
    row_names_gp = grid::gpar(fontsize = 10, fontface = "bold"),
    column_names_gp = grid::gpar(fontsize = 10, fontface = "bold"),
    column_title = paste0(
      project_name, ": Mark Enrichment (",
      n_changed, " reclassified anchors)"
    ),
    column_title_gp = grid::gpar(fontsize = 12, fontface = "bold"),
    cell_fun = cell_fun,
    rect_gp = grid::gpar(col = "white", lwd = 0.5),
    border = TRUE,
    left_annotation = ha_left,
    heatmap_legend_param = list(
      title = "% Positive",
      title_gp = grid::gpar(fontsize = 9, fontface = "bold"),
      labels_gp = grid::gpar(fontsize = 8),
      at = c(0, 25, 50, 75, 100),
      legend_height = unit(2, "cm")
    ),
    width = unit(length(available_marks) * 1.8, "cm"),
    height = unit(nrow(pct_mat) * 1.2, "cm")
  )
}

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

#' Internal: Compute H3K4me1 / H3K4me3 Signal Ratio from bigWig
#'
#' Extracts mean bigWig signal per anchor and computes the H3K4me1/H3K4me3
#' Quantifies relative H3K4me1 and H3K4me3 signal among dual-mark-positive
#' anchors.  A me1-dominant signature supports enhancer-like chromatin,
#' whereas a below-threshold signature is treated as promoter-leaning under
#' the operational cutoff.  This is a chromatin signature metric, not
#' evidence of dual regulatory function.  Only computes for anchors where
#' BOTH marks are BED-positive.
#'
#' @param gr_anchors A GRanges of candidate anchors with \code{anchor_id}.
#' @param mark_matrix Data frame with logical columns H3K4me1, H3K4me3.
#' @param bw_me1 Path to H3K4me1 bigWig file.
#' @param bw_me3 Path to H3K4me3 bigWig file.
#' @return A data frame with columns: \code{anchor_id}, \code{H3K4me1_anchor_mean},
#'   \code{H3K4me3_anchor_mean}, \code{H3K4me1_covered_mean},
#'   \code{H3K4me3_covered_mean}, \code{H3K4me1_AUC}, \code{H3K4me3_AUC},
#'   \code{H3K4me1_coverage_fraction}, \code{H3K4me3_coverage_fraction},
#'   \code{log2_ratio} (primary: anchor-mean-based), \code{log2_ratio_covered}
#'   (diagnostic: covered-mean-based), and \code{pseudocount}.
#' @keywords internal
#' @noRd
.compute_bw_ratio <- function(gr_anchors, mark_matrix, bw_me1, bw_me3) {
  .require_pkg("rtracklayer", "reading bigWig files", "stop")
  if (!file.exists(bw_me1)) stop("bigWig file not found: ", bw_me1)
  if (!file.exists(bw_me3)) stop("bigWig file not found: ", bw_me3)

  # dual_idx selects anchors where BOTH H3K4me1 and H3K4me3 BED-overlap
  # positively tested. We must NOT use isTRUE() here: isTRUE(x) returns
  # FALSE for any logical vector of length != 1, so on a multi-anchor
  # mark_matrix it would silently collapse dual_idx to integer(0),
  # disabling the entire bigWig ratio feature and forcing every true
  # dual-signature anchor to be misclassified as a simple promoter in
  # .chromatin_reclassify. Use the vectorised form with explicit NA
  # handling instead.
  dual_idx <- which(!is.na(mark_matrix$H3K4me1) & mark_matrix$H3K4me1 &
    !is.na(mark_matrix$H3K4me3) & mark_matrix$H3K4me3)
  if (length(dual_idx) == 0) {
    return(data.frame(
      anchor_id = character(),
      H3K4me1_anchor_mean = numeric(), H3K4me3_anchor_mean = numeric(),
      H3K4me1_covered_mean = numeric(), H3K4me3_covered_mean = numeric(),
      H3K4me1_AUC = numeric(), H3K4me3_AUC = numeric(),
      H3K4me1_coverage_fraction = numeric(),
      H3K4me3_coverage_fraction = numeric(),
      log2_ratio = numeric(), log2_ratio_covered = numeric(),
      pseudocount = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  gr_dual <- gr_anchors[dual_idx]
  anchor_w <- GenomicRanges::width(gr_dual)

  sig1 <- tryCatch(rtracklayer::import.bw(bw_me1, which = gr_dual),
    error = function(e) stop("Failed to read ", bw_me1, ": ", e$message)
  )
  sig3 <- tryCatch(rtracklayer::import.bw(bw_me3, which = gr_dual),
    error = function(e) stop("Failed to read ", bw_me3, ": ", e$message)
  )

  # Per-anchor signal metrics:
  #   AUC          = SUM(score_i x overlap_width_i)    total signal
  #   covered_mean = AUC / SUM(overlap_width_i)        density where mark exists
  #   anchor_mean  = AUC / anchor_width              density across full anchor
  #
  # For the H3K4me1/H3K4me3 ratio, anchor_mean is the correct metric:
  # both marks share the same denominator (anchor_width), so their
  # ratio simplifies to AUC1 / AUC3.  Covered-mean uses different
  # denominators per mark (different covered widths), which can invert
  # the dominance call.
  .anchor_metrics <- function(gr_dual, sig, anchor_w) {
    hits <- GenomicRanges::findOverlaps(gr_dual, sig)
    n <- length(gr_dual)
    covered_mean <- rep(NA_real_, n)
    anchor_mean <- rep(NA_real_, n)
    covered_w <- rep(0, n)
    auc <- rep(0, n)
    if (length(hits) == 0) {
      return(data.frame(
        covered_mean = covered_mean, anchor_mean = anchor_mean,
        covered_w = covered_w, auc = auc,
        stringsAsFactors = FALSE
      ))
    }
    qh <- S4Vectors::queryHits(hits)
    sh <- S4Vectors::subjectHits(hits)
    w <- GenomicRanges::width(GenomicRanges::pintersect(gr_dual[qh], sig[sh]))
    score_w <- sig$score[sh] * w
    auc[unique(qh)] <- tapply(score_w, qh, sum)
    covered_w[unique(qh)] <- tapply(w, qh, sum)
    covered_mean[unique(qh)] <- auc[unique(qh)] / covered_w[unique(qh)]
    anchor_mean[unique(qh)] <- auc[unique(qh)] / anchor_w[unique(qh)]
    data.frame(
      covered_mean = covered_mean, anchor_mean = anchor_mean,
      covered_w = covered_w, auc = auc,
      stringsAsFactors = FALSE
    )
  }
  m1 <- .anchor_metrics(gr_dual, sig1, anchor_w)
  m3 <- .anchor_metrics(gr_dual, sig3, anchor_w)

  # Data-driven pseudocount from anchor-wide means
  pooled <- c(m1$anchor_mean, m3$anchor_mean)
  pooled_pos <- pooled[pooled > 0 & !is.na(pooled)]
  pseudocount <- if (length(pooled_pos) >= 10) {
    max(stats::quantile(pooled_pos, 0.01, na.rm = TRUE), 1e-6)
  } else {
    1e-6
  }

  # Primary log2 ratio: anchor-wide mean with symmetric pseudocount
  log2_ratio <- log2(m1$anchor_mean + pseudocount) -
    log2(m3$anchor_mean + pseudocount)

  # Diagnostic: covered-mean log2 ratio (for comparison)
  log2_ratio_covered <- log2(m1$covered_mean + pseudocount) -
    log2(m3$covered_mean + pseudocount)

  n_na <- sum(is.na(log2_ratio))
  if (n_na > 0) {
    message(sprintf(
      "    bigWig ratio: %d / %d dual-positive anchors have NA ratio (signal missing in one or both marks; quantitative dominance unresolved)",
      n_na, length(log2_ratio)
    ))
  }
  message(sprintf(
    "    bigWig pseudocount: %.3g (1st percentile of pooled non-zero anchor-wide means, n=%d)",
    pseudocount, length(pooled_pos)
  ))
  data.frame(
    anchor_id = gr_dual$anchor_id,
    H3K4me1_anchor_mean = m1$anchor_mean,
    H3K4me3_anchor_mean = m3$anchor_mean,
    H3K4me1_covered_mean = m1$covered_mean,
    H3K4me3_covered_mean = m3$covered_mean,
    H3K4me1_AUC = m1$auc,
    H3K4me3_AUC = m3$auc,
    H3K4me1_coverage_fraction = m1$covered_w / anchor_w,
    H3K4me3_coverage_fraction = m3$covered_w / anchor_w,
    log2_ratio = log2_ratio,
    log2_ratio_covered = log2_ratio_covered,
    pseudocount = pseudocount,
    stringsAsFactors = FALSE
  )
}

#' @title Orthogonal validation of eP/eG reclassification by chromatin marks
#'
#' @description
#' Validates the expression-aware eP/eG reclassification produced by
#' \code{\link{refine_loop_anchors_by_expression}}, or the raw P/G/E
#' classification from \code{\link{annotate_peaks_and_loops}}, against
#' user-supplied chromatin mark BED files.
#'
#' @param annotation_res List. The raw foundational output object returned by
#'   \code{\link{annotate_peaks_and_loops}} or the refined output from
#'   \code{\link{refine_loop_anchors_by_expression}}.
#' @param chromatin_beds Named list of BED file paths. Accepts up to five
#'   canonical histone marks: \code{H3K4me1}, \code{H3K27ac}, \code{ATAC},
#'   \code{H3K27me3}, and \code{H3K4me3}. Unknown names are dropped with a
#'   warning; unmatched case is resolved to the canonical names. An empty
#'   list classifies every anchor as \code{"uncertain"}.
#' @param anchor_gap Integer. Gap tolerance for mark overlap.  Default
#'   \code{200}.
#' @param anchor_min_overlap Integer. Minimum overlap for mark overlap.
#'   Default \code{100}.
#' @param candidate_types Character vector or \code{NULL}. Anchor types to
#'   validate. \code{NULL} (default): uses all five types
#'   \code{c("P","G","E","eP","eG")} regardless of input source.
#' @param species Character. Genome assembly string. Default: \code{"hg38"}.
#'   Accepts any assembly string; used for seqlevel harmonization.
#' @param quiet Logical. Suppress messages. Default: \code{FALSE}.
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{anchor_id, chr, start, end, anchor_type, anchor_gene, cluster_id}{Anchor identifiers.}
#'   \item{H3K4me1, H3K27ac, ATAC, H3K27me3, H3K4me3}{Logical or \code{NA}. \code{TRUE} = peak detected, \code{FALSE} = peak not called (below detection threshold in this dataset), \code{NA} = mark not provided.}
#'   \item{enhancer_evidence}{Factor (enhancer evidence level): \code{canonical} > \code{strong} > \code{supported} > \code{limited} > \code{uncertain}.}
#'   \item{evidence}{Human-readable string of the constituting marks (e.g. \code{"H3K4me1+; H3K27ac+; H3K4me3-; H3K27me3-"}). Anchors with tested-positive H3K4me3 at \code{supported} confidence are annotated with a \code{"promoter_like"} tag.}
#' }
#'
#' @details
#' Enhancer evidence levels follow ENCODE-inspired active-enhancer criteria. Each
#' anchor is classified into the highest applicable level:
#' \describe{
#'   \item{canonical}{All 5 marks tested: H3K4me1+, H3K27ac+, ATAC+, H3K27me3-, H3K4me3-. This is a chromatin signature match, not functional validation.}
#'   \item{strong}{H3K4me1+ and (H3K27ac+ or ATAC+); H3K4me3 must be absent when tested.}
#'   \item{supported}{Any of H3K4me1, H3K27ac, or ATAC positive. Anchors with H3K4me3+
#'     receive a \code{"promoter_like"} evidence tag (note: the evidence tag alone does
#'     not force reclassification; triple-positive anchors retain eP/eG via
#'     \code{conflicting_marks}).}
#'   \item{limited}{Tested negative markers (H3K27me3-, H3K4me3-) present but no active marks.}
#'   \item{uncertain}{Marks tested but none identified; or no marks tested at all.}
#' }
#' Marks not provided are recorded as \code{NA} and never treated as negative evidence.
#' Case-insensitive name matching allows flexible input (e.g., \code{h3k4me1}
#' normalised to \code{H3K4me1}).
#'
#' @export
#' @examples
#' rdata_path <- system.file("extdata", "analysis_results.RData", package = "looplook")
#' tmp <- new.env()
#' load(rdata_path, envir = tmp)
#' raw_annotation <- tmp[[ls(tmp)[1]]]
#'
#' # Create dummy chromatin BED files for demonstration
#' bed_dir <- tempdir()
#' writeLines("chr1\t109492000\t109496500", file.path(bed_dir, "H3K4me1.bed"))
#' writeLines("chr1\t116686500\t116690000", file.path(bed_dir, "H3K27ac.bed"))
#'
#' # Run validation (using raw annotation; pass refined for eP/eG only)
#' result <- validate_epeG_by_chromatin(
#'   annotation_res = raw_annotation,
#'   chromatin_beds = list(
#'     H3K4me1 = file.path(bed_dir, "H3K4me1.bed"),
#'     H3K27ac = file.path(bed_dir, "H3K27ac.bed")
#'   ),
#'   quiet = TRUE
#' )
#' table(result$enhancer_evidence)
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
  log_message <- .make_log_message(quiet)

  # ---- 0. Validate inputs ----
  .validate_chromatin_overlap_inputs(
    chromatin_beds, candidate_types,
    anchor_gap, anchor_min_overlap
  )

  # ---- 1. Identify candidate anchors ----
  epeG_anchors <- .extract_epeG_anchors(
    annotation_res, log_message,
    candidate_types
  )

  # ---- 2. Validate chromatin_beds input ----
  known_marks <- c("H3K4me1", "H3K27ac", "ATAC", "H3K27me3", "H3K4me3")
  if (length(chromatin_beds) == 0) {
    warning("No chromatin_beds provided; all anchors classified as 'uncertain'.",
      call. = FALSE
    )
  }
  # Case-insensitive matching of user-provided mark names against canonical
  # names (e.g., "h3k4me1"  ->  "H3K4me1"). Unmatched names are warned and dropped.
  mark_lookup <- setNames(known_marks, toupper(known_marks))
  bed_names <- names(chromatin_beds)
  matched_idx <- toupper(bed_names) %in% names(mark_lookup)
  if (any(!matched_idx)) {
    warning("Unrecognised chromatin_beds name(s) ignored: ",
      paste(bed_names[!matched_idx], collapse = ", "),
      ". Expected: ", paste(known_marks, collapse = ", "),
      call. = FALSE
    )
  }
  names(chromatin_beds)[matched_idx] <- unname(mark_lookup[toupper(bed_names[matched_idx])])
  if (anyDuplicated(names(chromatin_beds))) {
    stop("chromatin_beds contains duplicate mark names after case-insensitive normalisation; ",
      "ensure each mark appears only once.",
      call. = FALSE
    )
  }
  provided_marks <- intersect(names(chromatin_beds), known_marks)

  # ---- 3+4. Overlap marks with anchors ----
  overlap_res <- .overlap_chromatin_marks(
    epeG_anchors, chromatin_beds, provided_marks, known_marks,
    anchor_gap, anchor_min_overlap, log_message
  )
  mark_matrix <- overlap_res$mark_matrix
  valid_provided_marks <- overlap_res$valid_provided_marks

  # ---- 5. Assign enhancer evidence levels and evidence ----
  result <- .assign_chromatin_confidence(
    epeG_anchors, mark_matrix, valid_provided_marks, known_marks
  )

  # ---- 6. Summary ----
  log_message("--- Validation Summary ---")
  tab <- table(result$enhancer_evidence, useNA = "ifany")
  for (lvl in names(tab)) {
    log_message(sprintf("  %-16s: %d anchors", lvl, tab[lvl]))
  }
  n_promoter_like <- sum(result$is_promoter_like, na.rm = TRUE)
  if (n_promoter_like > 0) {
    log_message(sprintf("  %-16s: %d anchors (H3K4me3+ promoter signal; may reflect active or bivalent promoter)", "promoter-like", n_promoter_like))
  }
  log_message("--- End Validation ---")

  out <- result %>%
    dplyr::arrange(enhancer_evidence, anchor_id) %>%
    dplyr::select(
      anchor_id, chr, start, end,
      anchor_type, anchor_gene, cluster_id,
      H3K4me1, H3K27ac, ATAC, H3K27me3, H3K4me3,
      enhancer_evidence, evidence
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
    enhancer_evidence = factor(levels = c(
      "canonical", "strong",
      "supported", "limited", "uncertain"
    )),
    evidence = character(),
    is_promoter_like = logical(),
    stringsAsFactors = FALSE
  )
}

#' Internal: Extract eP/eG (or P/G) anchors from annotation results
#' @keywords internal
#' @noRd
.extract_epeG_anchors <- function(annotation_res, log_message, candidate_types = NULL) {
  loop_df <- annotation_res$loop_annotation
  if (is.null(loop_df)) stop("'annotation_res$loop_annotation' is missing.")

  # Use canonical merged anchor coordinates from anchor_state when
  # available.  After anchor merge, a single anchor_id can appear with
  # multiple raw loop coordinates; reading from anchor_state ensures
  # one row per anchor ID with the correct merged hull coordinates.
  anchor_state <- attr(annotation_res, "looplook_anchor_state", exact = TRUE)
  if (!is.null(anchor_state) && !is.null(anchor_state$gr_anchors) &&
    !is.null(anchor_state$map_info)) {
    # Delegate to the canonical registry helper.  This ensures the
    # same anchor_id <-> coordinate mapping used by distal stats,
    # motif, and chromatin refinement.
    reg <- .get_anchor_registry(annotation_res)
    m <- S4Vectors::mcols(reg)
    anchor_map <- data.frame(
      anchor_id = as.character(m$anchor_id),
      chr = as.character(GenomicRanges::seqnames(reg)),
      start = as.integer(GenomicRanges::start(reg)),
      end = as.integer(GenomicRanges::end(reg)),
      anchor_type = as.character(m$type_code),
      anchor_gene = as.character(m$SYMBOL),
      cluster_id = if ("cluster_id" %in% colnames(m)) {
        as.character(m$cluster_id)
      } else {
        NA_character_
      },
      stringsAsFactors = FALSE
    )
    anchor_map <- anchor_map[!is.na(anchor_map$anchor_id) &
      anchor_map$anchor_id != "", ]
  } else if (all(c("a1_id", "a2_id") %in% colnames(loop_df))) {
    anchor_map <- dplyr::bind_rows(
      loop_df %>% dplyr::select(
        anchor_id = a1_id, chr = chr1,
        start = start1, end = end1, anchor_type = anchor1_type,
        anchor_gene = anchor1_gene, cluster_id
      ),
      loop_df %>% dplyr::select(
        anchor_id = a2_id, chr = chr2,
        start = start2, end = end2, anchor_type = anchor2_type,
        anchor_gene = anchor2_gene, cluster_id
      )
    ) %>% dplyr::distinct(anchor_id, .keep_all = TRUE)
  } else {
    anchor_map <- dplyr::bind_rows(
      loop_df %>% dplyr::mutate(anchor_id = paste(chr1, start1, end1,
        sep = "_"
      )) %>% dplyr::select(anchor_id,
        chr = chr1,
        start = start1, end = end1, anchor_type = anchor1_type,
        anchor_gene = anchor1_gene, cluster_id
      ),
      loop_df %>% dplyr::mutate(anchor_id = paste(chr2, start2, end2,
        sep = "_"
      )) %>% dplyr::select(anchor_id,
        chr = chr2,
        start = start2, end = end2, anchor_type = anchor2_type,
        anchor_gene = anchor2_gene, cluster_id
      )
    ) %>% dplyr::distinct(anchor_id, .keep_all = TRUE)
  }

  type_label <- "P/G/E/eP/eG"
  type_filter <- if (!is.null(candidate_types)) {
    candidate_types
  } else {
    # Always validate all positional types (P/G/E) plus any expression-
    # silenced types (eP/eG).  Chromatin evidence can identify dual-
    # function elements, unannotated promoters, and intronic enhancers
    # regardless of transcriptional state.  Restricting to eP/eG on
    # refined input would miss biologically important reclassifications
    # such as P->dual (promoter-proximal enhancer) or E->P (unannotated
    # promoter with H3K4me3+).
    c("P", "G", "E", "eP", "eG")
  }
  epeG_anchors <- anchor_map %>% dplyr::filter(anchor_type %in% type_filter)
  if (nrow(epeG_anchors) == 0) {
    message("No ", type_label, " anchors found. Returning empty result.")
    return(data.frame(
      anchor_id = character(), chr = character(),
      start = integer(), end = integer(), anchor_type = character(),
      anchor_gene = character(), cluster_id = character(),
      stringsAsFactors = FALSE
    ))
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
  if (inherits(epeG_anchors, "data.frame") && nrow(epeG_anchors) == 0) {
    return(list(
      mark_matrix = data.frame(),
      valid_provided_marks = character(0)
    ))
  }
  gr_anchors <- .with_known_upstream_noise_suppressed(
    GenomicRanges::makeGRangesFromDataFrame(epeG_anchors, keep.extra.columns = TRUE)
  )
  mark_matrix <- as.data.frame(
    matrix(NA,
      nrow = nrow(epeG_anchors), ncol = length(known_marks),
      dimnames = list(NULL, known_marks)
    )
  )
  # Validate required marks before reading any data.  A mark that is
  # named in chromatin_beds but has an invalid path must NOT be silently
  # encoded as not_called (FALSE) -- that would contaminate reclassification
  # rules with false negative evidence.
  required_marks <- c("H3K4me1", "H3K4me3")
  path_ok <- vapply(chromatin_beds[provided_marks], function(x) {
    is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x) && file.exists(x)
  }, logical(1))
  invalid_marks <- provided_marks[!path_ok]
  bad_required <- intersect(required_marks, invalid_marks)
  if (length(bad_required) > 0) {
    stop("Required chromatin BED file(s) not found or invalid: ",
      paste(bad_required, collapse = ", "),
      call. = FALSE
    )
  }
  if (length(invalid_marks) > 0) {
    warning("Optional chromatin BED file(s) not found, skipped: ",
      paste(invalid_marks, collapse = ", "),
      call. = FALSE
    )
  }
  valid_provided_marks <- provided_marks[path_ok]
  for (mark in valid_provided_marks) {
    bed_path <- chromatin_beds[[mark]]
    log_message(sprintf("  Overlapping with %s ...", mark))
    mark_gr <- read_simple_bed(bed_path, quiet = TRUE)
    if (length(mark_gr) == 0) {
      if (mark %in% required_marks) {
        stop("Required chromatin mark `", mark,
          "` BED file contains zero peaks. An empty required-mark ",
          "file cannot safely be interpreted as genome-wide ",
          "negative evidence. Check assay QC and peak calling.",
          call. = FALSE
        )
      }
      log_message(
        "    BED file for ", mark, " contains zero peaks; ",
        "all anchors will be marked as tested-but-not-called (FALSE)."
      )
      mark_matrix[[mark]][] <- FALSE
      next
    }
    mark_gr <- .harmonize_seqlevels(mark_gr, gr_anchors, mark)
    anchor_chrs <- unique(as.character(GenomicRanges::seqnames(gr_anchors)))
    mark_chrs <- unique(as.character(GenomicRanges::seqnames(mark_gr)))
    shared_chrs <- intersect(anchor_chrs, mark_chrs)
    if (length(shared_chrs) == 0L) {
      stop("Chromatin mark `", mark,
        "` has no seqlevels in common with candidate anchors. ",
        "Check genome assembly and chromosome naming style.",
        call. = FALSE
      )
    }
    hits <- GenomicRanges::findOverlaps(gr_anchors, mark_gr, maxgap = anchor_gap)
    if (anchor_min_overlap >= 1L && length(hits) > 0) {
      q_gr <- gr_anchors[S4Vectors::queryHits(hits)]
      s_gr <- mark_gr[S4Vectors::subjectHits(hits)]
      overlap_w <- GenomicRanges::width(GenomicRanges::pintersect(q_gr, s_gr))
      hits <- hits[overlap_w >= anchor_min_overlap]
    }
    hit_anchors <- unique(S4Vectors::queryHits(hits))
    mark_matrix[[mark]][hit_anchors] <- TRUE
    n_hit <- length(hit_anchors)
    log_message(sprintf(
      "    %d / %d anchors overlap %s peaks",
      n_hit, nrow(epeG_anchors), mark
    ))
    if (n_hit == 0L && mark %in% required_marks) {
      warning("Required chromatin mark `", mark,
        "` shares seqlevels with candidate anchors but no anchor ",
        "passes the overlap criteria. All candidate anchors will be ",
        "recorded as tested-but-not-called. Users are responsible for ",
        "ensuring compatible genome assemblies, chromosome naming, ",
        "and overlap parameters.",
        call. = FALSE
      )
    }
  }
  # Mark status is three-state:
  #   TRUE       = detected (peak overlaps anchor in provided BED)
  #   FALSE      = not_called (BED provided but no peak at this anchor)
  #   NA         = not_tested (BED not provided)
  # FALSE != biological absence.  Peak not called may reflect sequencing
  # depth, peak-caller threshold, mappability, or other technical factors.
  # Downstream rules that use "not_called" as exclusionary evidence should
  # treat it as moderate/low-confidence negative evidence, not definitive.
  for (mark in valid_provided_marks) {
    na_idx <- is.na(mark_matrix[[mark]])
    mark_matrix[[mark]][na_idx] <- FALSE
  }
  list(mark_matrix = mark_matrix, valid_provided_marks = valid_provided_marks)
}

#' Internal: Assign enhancer evidence levels and evidence strings
#' @keywords internal
#' @noRd
.assign_chromatin_confidence <- function(epeG_anchors, mark_matrix,
                                         provided_marks, known_marks) {
  if (nrow(mark_matrix) == 0) {
    return(.empty_validation_result())
  }
  result <- epeG_anchors
  result$H3K4me1 <- mark_matrix$H3K4me1
  result$H3K27ac <- mark_matrix$H3K27ac
  result$ATAC <- mark_matrix$ATAC
  result$H3K27me3 <- mark_matrix$H3K27me3
  result$H3K4me3 <- mark_matrix$H3K4me3

  has_positive <- vapply(seq_len(nrow(result)), function(i) {
    isTRUE(result$H3K4me1[i]) || isTRUE(result$H3K27ac[i]) || isTRUE(result$ATAC[i])
  }, logical(1))
  all_five <- length(provided_marks) == 5
  .is_absent <- function(x) !is.na(x) && !x

  result$enhancer_evidence <- vapply(seq_len(nrow(result)), function(i) {
    h3k4me1_t <- !is.na(result$H3K4me1[i])
    h3k27ac_t <- !is.na(result$H3K27ac[i])
    atac_t <- !is.na(result$ATAC[i])
    h3k27me3_t <- !is.na(result$H3K27me3[i])
    h3k4me3_t <- !is.na(result$H3K4me3[i])
    neg <- c(
      if (.is_absent(result$H3K27me3[i])) "H3K27me3-",
      if (.is_absent(result$H3K4me3[i])) "H3K4me3-"
    )
    if (all_five && h3k4me1_t && h3k27ac_t && atac_t && h3k27me3_t && h3k4me3_t &&
      result$H3K4me1[i] && result$H3K27ac[i] && result$ATAC[i] &&
      .is_absent(result$H3K27me3[i]) && .is_absent(result$H3K4me3[i])) {
      return("canonical")
    }
    if (h3k4me1_t && result$H3K4me1[i] &&
      # H3K4me3, if tested, must be absent: an anchor with both
      # H3K4me1 and H3K4me3 is promoter-like, not a distal enhancer.
      # When H3K4me3 is not provided (NA), we do not penalise.
      (is.na(result$H3K4me3[i]) || .is_absent(result$H3K4me3[i])) &&
      ((h3k27ac_t && result$H3K27ac[i]) || (atac_t && result$ATAC[i]))) {
      return("strong")
    }
    if (has_positive[i]) {
      return("supported")
    }
    if (length(neg) > 0) {
      return("limited")
    }
    return("uncertain")
  }, character(1))

  result$enhancer_evidence <- factor(result$enhancer_evidence,
    levels = c("canonical", "strong", "supported", "limited", "uncertain")
  )

  result$evidence <- vapply(seq_len(nrow(result)), function(i) {
    parts <- c()
    if (isTRUE(result$H3K4me1[i])) parts <- c(parts, "H3K4me1+")
    if (isTRUE(result$H3K27ac[i])) parts <- c(parts, "H3K27ac+")
    if (isTRUE(result$ATAC[i])) parts <- c(parts, "ATAC+")
    if (isTRUE(result$H3K27me3[i])) parts <- c(parts, "H3K27me3+")
    if (isTRUE(result$H3K4me3[i])) parts <- c(parts, "H3K4me3+")
    if (isTRUE(!is.na(result$H3K27me3[i]) && !result$H3K27me3[i])) {
      parts <- c(parts, "H3K27me3-")
    }
    if (isTRUE(!is.na(result$H3K4me3[i]) && !result$H3K4me3[i])) {
      parts <- c(parts, "H3K4me3-")
    }
    # Flag: H3K4me3 tested & positive indicates promoter signal, regardless
    # of whether H3K27me3 is absent (active promoter) or present (bivalent/
    # poised promoter -- hallmark of silenced developmental genes). Both
    # should be reverted to P by chromatin reclassification, not kept as eP.
    # Promoter evidence is independent of enhancer evidence level; H3K4me3+
    # alone with no enhancer marks is still strong promoter evidence.
    if (isTRUE(result$H3K4me3[i])) {
      parts <- c(parts, "promoter_like")
    }
    if (length(parts) == 0) {
      return("no_data")
    }
    paste(parts, collapse = "; ")
  }, character(1))
  # H3K4me3+ alone is sufficient promoter-like evidence, independent of
  # enhancer evidence tier.  A promoter with H3K4me3 but without
  # H3K4me1/H3K27ac/ATAC is still a promoter; it should not be
  # misclassified as having "uncertain" regulatory identity.
  result$is_promoter_like <- !is.na(result$H3K4me3) & result$H3K4me3
  result
}

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "contact_id", "contact_id1", "contact_id2",
    "Distal_Type", "dual_ratio_state",
    "Gene_Assignment_Confidence", "Gene_Assignment_Evidence",
    "h3k4me1_not_called", "h3k4me3_not_called",
    "is_promoter_after",
    "loop_ID_path",
    "measured",
    "n_Linked_EnhancerLike", "n_Linked_EnhancerLike_Filtered",
    "Expanded_Target_Genes",
    "n_Unique_Contacts", "n_Unique_Contacts_Filtered",
    "pair", "path_length", "path_rank",
    "proximate_promoter_gene",
    "res", "res_promoter", "role_rank",
    "TSS_support_status", "TSS_supported",
    "TSS_supported", "Gene_Assignment_Confidence",
    "was_promoter_before"
  ))
}
