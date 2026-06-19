#' Internal Package Imports
#'
#' @name looplook_imports
#' @noRd
#' @importFrom grDevices col2rgb rgb
#' @importFrom stats fisher.test median na.omit p.adjust quantile reorder runif setNames t.test var wilcox.test
#' @importFrom utils head read.table write.csv
NULL

if (getRversion() >= "2.15.1") {
    utils::globalVariables(c(
        ".", ".data", "All_Anchor_Genes", "All_Loop_Connected_Genes", "Assigned_Target_Genes",
        "CleanLoopType", "Conn_Group", "Count", "Degree", "Description",
        "Description_unique", "Distal_Anchor_ID", "Dominant_Interaction",
        "Dominant_Interaction_Filtered", "Expression", "Expression_Status", "FDR",
        "Feature", "Filtered", "Fraction", "GENENAME", "GENETYPE", "Gene",
        "Genomic_Distribution", "Group", "High_Connectivity_Gene",
        "Is_Active_Gene", "Is_High_Connectivity_Distal_Element",
        "Is_High_Connectivity_Gene", "Is_High_Distal_Connectivity_Gene", "L1_Raw",
        "L2_Raw", "L3_Raw", "LFC", "Label", "LabelText", "Label_Text",
        "Linked_Loop_IDs", "Log10Degree", "LogFDR", "Loop_Type",
        "Mean_Expression_Temp", "MotifLabel", "ONTOLOGY", "OddsRatio", "Original",
        "PlotFamily", "Putative_Target_Genes", "Rank",
        "Regulated_promoter_genes", "SANKEY_RAW_GENES", "SYMBOL", "SampleID",
        "Simplified", "Source", "Target_Genes",
        "Target_Genes_Filtered", "Total_Loops", "Total_Loops_Filtered",
        "Unique_Gene_Count", "a1_id", "a2_id", "Active_Target_Genes", "all_cluster_loop_genes", "all_of",
        "anchor1_gene", "anchor1_type", "anchor2_gene",
        "anchor2_type", "anchor_id", "annotation", "chr",
        "cluster_id", "col2rgb", "combined_score", "count", "deg", "detail_anno",
        "elementNROWS", "everything", "expansion", "final_color", "final_fill",
        "final_label", "final_symbols", "fisher.test", "fraction",
        "geneList",
        "gene_id", "gene_level", "geom_hline", "group", "has_active", "head",
        "is_active",
        "hjust", "install.packages", "is_e_type", "is_lower_e", "label",
        "label_text", "label_x", "len", "lfc", "linked_loops",
        "logP", "loop_ID", "loop_genes", "loop_i",
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
        "Loop_Connection", "Neighbor_Gene", "Neighbor_Type", "s1", "s2", "x", "y",
        "Conn_Group_jitter", "Conn_Group_num", "Conn_Group_slab",
        "left_mid", "right_mid", ".fallback_ptg",
        "Regulated_promoter_Evidence", "Regulated_promoter_Fallback_Evidence",
        "Refined", "Retained_In_Functional_Network",
        "gene", "input_id", "evidence", "gene_role", "source", "anchor_role", "used_as_fallback",
        "in_regulated_promoter", "in_assigned_target", "in_all_loop_connected",
        "in_regulated_promoter_filled", "in_assigned_target_filled",
        "opposite_anchor_id", "local_anchor_id", "Mean_Expression",
        "Passes_Expression_Filter", "retained_after_refinement",
        "H3K4me1", "H3K27ac", "ATAC", "H3K27me3", "H3K4me3",
        "anchor_type", "anchor_gene", "confidence", "evidence",
        "priority"
    ))
}

# Only suppress specific non-actionable upstream noise from third-party packages.
# Only suppress truly non-actionable third-party noise (deprecation warnings
# from upstream packages). Warnings about data integrity -- including
# out-of-bound ranges, which often indicate genome-build mismatches or
# coordinate errors -- must remain visible.
.with_known_upstream_noise_suppressed <- function(expr) {
    withCallingHandlers(
        expr,
        warning = function(w) {
            msg <- conditionMessage(w)
            if (
                grepl("S4Vectors:::anyMissing\\(\\).*deprecated", msg) ||
                    grepl("`aes_()` was deprecated in ggplot2 3.0.0.", msg, fixed = TRUE) ||
                    grepl("`aes_string()` was deprecated in ggplot2 3.0.0.", msg, fixed = TRUE) ||
                    grepl("Using `size` aesthetic for lines was deprecated", msg, fixed = TRUE) ||
                    grepl("The `size` argument of `element_line()` is deprecated", msg, fixed = TRUE)
            ) {
                invokeRestart("muffleWarning")
            }
        },
        message = function(m) {
            msg <- conditionMessage(m)
            if (
                grepl("'select\\(\\)' returned 1:(1|many) mapping between keys and columns", msg) ||
                    grepl("genes were dropped because they have exons located on both strands", msg, fixed = TRUE) ||
                    grepl("Scale for colour is already present.", msg, fixed = TRUE) ||
                    grepl("Adding another scale for colour, which will replace the existing scale.", msg, fixed = TRUE)
            ) {
                invokeRestart("muffleMessage")
            }
        }
    )
}

# Quiet mode should only silence informational messages; warnings still surface.
.with_messages_silenced <- function(expr) {
    withCallingHandlers(
        expr,
        message = function(m) {
            invokeRestart("muffleMessage")
        }
    )
}

#' Internal: Harmonize Seqlevels Style
#'
#' Converts the seqlevels style of \code{gr} to match \code{ref_gr}.
#' Emits a message when conversion is performed and a warning if conversion fails.
#'
#' @param gr A GRanges object to potentially convert.
#' @param ref_gr A GRanges object whose seqlevels style is the target.
#' @param label Character. Label for diagnostic messages (e.g., "blacklist", "ROI").
#' @return The input \code{gr} with seqlevels style matching \code{ref_gr} (if possible).
#' @keywords internal
#' @noRd
.harmonize_seqlevels <- function(gr, ref_gr, label = "") {
    if (length(gr) == 0 || length(ref_gr) == 0) {
        return(gr)
    }

    style_gr <- tryCatch(GenomeInfoDb::seqlevelsStyle(gr), error = function(e) NULL)
    style_ref <- tryCatch(GenomeInfoDb::seqlevelsStyle(ref_gr), error = function(e) NULL)

    if (is.null(style_gr) || is.null(style_ref)) {
        return(gr)
    }

    if (length(style_gr) > 0 && length(style_ref) > 0 && style_gr[1] != style_ref[1]) {
        overlap_before <- length(GenomicRanges::intersect(
            GenomeInfoDb::seqlevels(gr), GenomeInfoDb::seqlevels(ref_gr)
        ))
        tryCatch(
            {
                GenomeInfoDb::seqlevelsStyle(gr) <- style_ref[1]
                overlap_after <- length(GenomicRanges::intersect(
                    GenomeInfoDb::seqlevels(gr), GenomeInfoDb::seqlevels(ref_gr)
                ))
                if (nzchar(label)) {
                    message(
                        "Seqlevels style harmonized for ", label, ": ",
                        style_gr[1], " -> ", style_ref[1],
                        " (overlapping seqlevels: ", overlap_before, " -> ", overlap_after, ")"
                    )
                }
            },
            error = function(e) {
                warning("Failed to harmonize seqlevels for ", label, ": ",
                    conditionMessage(e),
                    call. = FALSE
                )
            }
        )
    }
    gr
}

#' Internal: Filter Target Regions Overlapping Loop Anchors
#'
#' Returns the subset of target genomic ranges that overlap at least one
#' loop anchor. Used to focus annotation on features with 3D connectivity.
#'
#' @param target A GRanges object of target features (e.g., ChIP-seq peaks).
#' @param anchors A GRanges object of loop anchors.
#' @param min_overlap Integer. Minimum number of overlapping anchors required
#'   to retain a target. Default: \code{1L}.
#' @return The filtered \code{target} GRanges containing only regions that
#'   overlap at least \code{min_overlap} anchors.
#' @keywords internal
#' @noRd
.filter_target_anchor_overlaps <- function(target, anchors, min_overlap = 1L) {
    if (length(target) == 0 || length(anchors) == 0) {
        return(target[0])
    }
    anchors <- .harmonize_seqlevels(anchors, target, "anchors")
    hits <- GenomicRanges::findOverlaps(target, anchors)
    if (length(hits) == 0) {
        return(target[0])
    }
    hit_counts <- table(S4Vectors::queryHits(hits))
    keep <- as.integer(names(hit_counts[hit_counts >= min_overlap]))
    target[sort(keep)]
}

#' Internal: Safe FindOverlaps with Seqlevels Harmonization
#'
#' Wrapper around \code{GenomicRanges::findOverlaps} that automatically
#' harmonizes seqlevels style before computing overlaps.
#'
#' @param query A GRanges or GInteractions object.
#' @param subject A GRanges object.
#' @param label Character. Label for diagnostic messages.
#' Internal: Clean Gene Name Vector
#'
#' Removes empty strings, NA values, and duplicate entries from gene identifiers.
#' Optionally splits concatenated strings (e.g., "TP53;BRCA1") before cleaning.
#'
#' @details
#' When no valid genes remain after cleaning, returns a zero-length character
#' vector (\code{character(0)}). This differs from \code{\link{extract_genes}},
#' which returns \code{NA_character_} in that case. Callers should use
#' \code{length(x) > 0} rather than \code{!is.na(x)} to test for empty results.
#'
#' @param x Character vector of gene names, possibly containing delimiters.
#' @param split Character. If non-NULL, a regex passed to \code{\link{strsplit}}
#'   to split concatenated gene strings (e.g., \code{"[;,]"}). Set to \code{NULL}
#'   if \code{x} is already a clean character vector.
#' @return A unique, non-empty, non-NA character vector, or \code{character(0)}
#'   if no valid genes remain.
#' @keywords internal
#' @noRd
clean_gene_names <- function(x, split = NULL) {
    if (is.null(x) || length(x) == 0) {
        return(character(0))
    }
    if (!is.null(split)) x <- unlist(strsplit(as.character(x), split))
    x <- unique(trimws(as.character(x)))
    x[x != "" & !is.na(x)]
}

.get_org_db_obj <- function(org_db) {
    if (any(inherits(org_db, c("OrgDb", "AnnotationDb")))) {
        return(org_db)
    }
    if (is.character(org_db) && length(org_db) == 1L && nzchar(org_db)) {
        if (!requireNamespace(org_db, quietly = TRUE)) {
            stop("Package '", org_db, "' is required but not installed.")
        }
        return(utils::getFromNamespace(org_db, org_db))
    }
    stop("`org_db` must be an OrgDb/AnnotationDb object or an installed package name.")
}

.extract_txdb_gene_ids <- function(genes_gr) {
    if (is.null(genes_gr) || length(genes_gr) == 0) {
        return(character(0))
    }
    gene_ids <- if ("gene_id" %in% colnames(S4Vectors::mcols(genes_gr))) {
        S4Vectors::mcols(genes_gr)$gene_id
    } else {
        names(genes_gr)
    }
    gene_ids <- trimws(as.character(gene_ids))
    gene_ids[gene_ids == "" | is.na(gene_ids)] <- NA_character_
    gene_ids
}

.empty_orgdb_gene_map <- function(gene_ids, columns) {
    gene_ids <- clean_gene_names(gene_ids)
    out <- data.frame(gene_id = gene_ids, stringsAsFactors = FALSE)
    for (col in columns) out[[col]] <- NA_character_
    out
}

.detect_orgdb_keytype <- function(
  gene_ids, org_db,
  score_column = "SYMBOL",
  candidate_keytypes = c("ENTREZID", "ENSEMBL", "SYMBOL", "GENENAME", "REFSEQ")
) {
    gene_ids <- clean_gene_names(gene_ids)
    org_db_obj <- .get_org_db_obj(org_db)
    if (length(gene_ids) == 0) {
        return(list(
            keytype = NA_character_,
            hit_rate = 0,
            hits = integer(0),
            score_column = score_column,
            org_db_obj = org_db_obj
        ))
    }

    valid_keys <- AnnotationDbi::keytypes(org_db_obj)
    valid_cols <- AnnotationDbi::columns(org_db_obj)
    score_column <- if (score_column %in% valid_cols) score_column else valid_cols[1]
    candidate_keytypes <- intersect(candidate_keytypes, valid_keys)
    if (length(candidate_keytypes) == 0) {
        return(list(
            keytype = NA_character_,
            hit_rate = 0,
            hits = integer(0),
            score_column = score_column,
            org_db_obj = org_db_obj
        ))
    }

    hit_counts <- vapply(candidate_keytypes, function(keytype) {
        mapped <- tryCatch(
            withCallingHandlers(
                AnnotationDbi::mapIds(
                    org_db_obj,
                    keys = gene_ids,
                    column = score_column,
                    keytype = keytype,
                    multiVals = "first"
                ),
                warning = function(w) {
                    msg <- conditionMessage(w)
                    if (grepl("None of the keys entered are valid keys for", msg, fixed = TRUE)) {
                        invokeRestart("muffleWarning")
                    }
                }
            ),
            error = function(e) {
                setNames(rep(NA_character_, length(gene_ids)), gene_ids)
            }
        )
        sum(!is.na(mapped) & mapped != "")
    }, integer(1))

    if (length(hit_counts) == 0 || max(hit_counts) == 0) {
        return(list(
            keytype = NA_character_,
            hit_rate = 0,
            hits = hit_counts,
            score_column = score_column,
            org_db_obj = org_db_obj
        ))
    }

    best_idx <- which.max(hit_counts)
    list(
        keytype = candidate_keytypes[[best_idx]],
        hit_rate = hit_counts[[best_idx]] / length(gene_ids),
        hits = hit_counts,
        score_column = score_column,
        org_db_obj = org_db_obj
    )
}

.map_txdb_gene_ids <- function(
  gene_ids, org_db, columns = "SYMBOL",
  context = "TxDb gene_id mapping",
  candidate_keytypes = c("ENTREZID", "ENSEMBL", "SYMBOL", "GENENAME", "REFSEQ"),
  warn = TRUE,
  min_hit_rate = 0.25
) {
    gene_ids <- clean_gene_names(gene_ids)
    org_db_obj <- .get_org_db_obj(org_db)
    valid_cols <- AnnotationDbi::columns(org_db_obj)
    columns <- intersect(columns, valid_cols)
    if (length(columns) == 0) {
        out <- .empty_orgdb_gene_map(gene_ids, "SYMBOL")
        attr(out, "keytype") <- NA_character_
        attr(out, "hit_rate") <- 0
        return(out)
    }

    det <- .detect_orgdb_keytype(
        gene_ids = gene_ids,
        org_db = org_db_obj,
        score_column = if ("SYMBOL" %in% columns) "SYMBOL" else columns[1],
        candidate_keytypes = candidate_keytypes
    )
    if (is.na(det$keytype)) {
        if (warn && length(gene_ids) > 0) {
            warning(
                sprintf(
                    "%s: unable to match TxDb gene_id values against supported OrgDb keytypes; raw gene_id values will be retained where needed.",
                    context
                ),
                call. = FALSE
            )
        }
        out <- .empty_orgdb_gene_map(gene_ids, columns)
        attr(out, "keytype") <- NA_character_
        attr(out, "hit_rate") <- 0
        return(out)
    }

    if (warn && det$hit_rate < min_hit_rate) {
        warning(
            sprintf(
                "%s: low OrgDb mapping rate for TxDb gene_id values (best keytype = %s, %.1f%% matched).",
                context, det$keytype, det$hit_rate * 100
            ),
            call. = FALSE
        )
    }

    raw_map <- tryCatch(
        .with_known_upstream_noise_suppressed(AnnotationDbi::select(
            org_db_obj,
            keys = gene_ids,
            columns = columns,
            keytype = det$keytype
        )),
        error = function(e) NULL
    )
    if (is.null(raw_map) || nrow(raw_map) == 0) {
        out <- .empty_orgdb_gene_map(gene_ids, columns)
        attr(out, "keytype") <- det$keytype
        attr(out, "hit_rate") <- det$hit_rate
        return(out)
    }

    raw_map$gene_id <- as.character(raw_map[[det$keytype]])
    raw_map <- raw_map[, c("gene_id", setdiff(colnames(raw_map), c(det$keytype, "gene_id"))), drop = FALSE]
    attr(raw_map, "keytype") <- det$keytype
    attr(raw_map, "hit_rate") <- det$hit_rate
    raw_map
}


#' Internal: Update Anchor State After Expression Refinement
#'
#' Synchronises \code{anchor_state$map_info$type_code} with the expression-refined
#' anchor types in \code{loop_df} so that downstream chromatin-aware target
#' recomputation uses the updated eP/eG classifications.
#'
#' @param anchor_state The \code{looplook_anchor_state} attribute from
#'   \code{\link{annotate_peaks_and_loops}} output.
#' @param loop_df Refined loop annotation data frame with \code{a1_id},
#'   \code{a2_id}, \code{anchor1_type}, \code{anchor2_type} columns.
#' @return The \code{anchor_state} list with \code{map_info$type_code} updated.
#' @keywords internal
#' @noRd
.update_anchor_state_from_loop_df <- function(anchor_state, loop_df) {
    if (is.null(anchor_state) || !"map_info" %in% names(anchor_state)) {
        return(anchor_state)
    }
    required_cols <- c("a1_id", "a2_id", "anchor1_type", "anchor2_type")
    if (!all(required_cols %in% colnames(loop_df))) {
        return(anchor_state)
    }

    type_map <- dplyr::bind_rows(
        loop_df %>%
            dplyr::transmute(
                anchor_id = as.character(.data$a1_id),
                type_code_new = as.character(.data$anchor1_type)
            ),
        loop_df %>%
            dplyr::transmute(
                anchor_id = as.character(.data$a2_id),
                type_code_new = as.character(.data$anchor2_type)
            )
    ) %>%
        dplyr::filter(!is.na(.data$anchor_id), .data$anchor_id != "") %>%
        dplyr::distinct(.data$anchor_id, .keep_all = TRUE)

    map_info <- anchor_state$map_info
    idx <- match(map_info$anchor_id, type_map$anchor_id)
    hit <- !is.na(idx)
    map_info$type_code[hit] <- type_map$type_code_new[idx[hit]]

    anchor_state$map_info <- map_info
    anchor_state
}

#' Internal: Collapse Delimited Gene String
#'
#' Splits a semicolon-delimited gene string, removes duplicates and NAs,
#' and recollapses into a single string. Returns \code{NA_character_} if
#' no valid genes remain.
#'
#' @param genes_vec Character vector of delimited gene strings.
#' @return A single semicolon-delimited string, or \code{NA_character_}.
#' @keywords internal
#' @noRd
extract_genes <- function(genes_vec) {
    res <- unique(na.omit(trimws(unlist(strsplit(as.character(genes_vec), ";")))))
    res <- res[nzchar(res)]
    if (length(res) == 0) {
        return(NA_character_)
    }
    paste(res, collapse = ";")
}


#' Internal: Resolve TxDb package name from species
#'
#' Maps a genome assembly shorthand to the corresponding Bioconductor
#' TxDb annotation package name.
#'
#' @param species Character. One of \code{"hg38"}, \code{"hg19"},
#'   \code{"mm10"}, \code{"mm9"}.
#' @return Character. TxDb package name.
#' @keywords internal
#' @noRd
species_txdb_pkg <- function(species) {
    switch(species,
        hg38 = "TxDb.Hsapiens.UCSC.hg38.knownGene",
        hg19 = "TxDb.Hsapiens.UCSC.hg19.knownGene",
        mm10 = "TxDb.Mmusculus.UCSC.mm10.knownGene",
        mm9  = "TxDb.Mmusculus.UCSC.mm9.knownGene",
        stop("Species not supported: ", species)
    )
}

#' Internal: Resolve OrgDb package name from species
#'
#' @param species Character. One of \code{"hg38"}, \code{"hg19"},
#'   \code{"mm10"}, \code{"mm9"}.
#' @return Character. OrgDb package name.
#' @keywords internal
#' @noRd
species_orgdb_pkg <- function(species) {
    switch(species,
        hg38 = "org.Hs.eg.db",
        hg19 = "org.Hs.eg.db",
        mm10 = "org.Mm.eg.db",
        mm9 = "org.Mm.eg.db",
        stop("Species not supported: ", species)
    )
}

#' Internal: Resolve BSgenome package name from species
#'
#' @param species Character. One of \code{"hg38"}, \code{"hg19"},
#'   \code{"mm10"}, \code{"mm9"}.
#' @return Character. BSgenome package name, or \code{NULL}.
#' @keywords internal
#' @noRd
species_bsgenome_pkg <- function(species) {
    switch(species,
        hg38 = "BSgenome.Hsapiens.UCSC.hg38",
        hg19 = "BSgenome.Hsapiens.UCSC.hg19",
        mm10 = "BSgenome.Mmusculus.UCSC.mm10",
        mm9  = "BSgenome.Mmusculus.UCSC.mm9",
        NULL
    )
}

#' Internal: Resolve Gene Conflicts via Biotype Priority Then Expression
#'
#' For each genomic range, identifies all promoter-overlapping genes,
#' resolves conflicts using a two-stage strategy: (1) biotype priority
#' (protein-coding > small-ncRNA > antisense > lncRNA/ncRNA > pseudogene),
#' optionally overridden by \code{biotype_order}, and
#' then (2) expression-aware filtering within the selected biotype tier.
#' If any gene in the best tier is expressed (\code{tpm >= min_expr}), only
#' expressed candidates are retained; otherwise all candidates in that tier
#' are kept. When multiple candidates share the same biotype rank, a
#' co-dominant expression rule is applied: all genes with expression >= 10\%
#' of the group maximum are retained (collapsed with ";").
#'
#' @param current_anno_df Data frame with columns suitable for
#'   \code{\link[GenomicRanges]{makeGRangesFromDataFrame}}.
#' @param txdb_obj A \code{TxDb} object for gene coordinate lookup.
#' @param org_db_pkg Character. Organism database package name.
#' @param tss_region Numeric vector of length 2. TSS region for promoter
#'   definition, e.g. \code{c(-2000, 2000)}.
#' @param gene_expr_map Named numeric vector of per-gene expression values,
#'   or \code{NULL} if unavailable.
#' @param min_expr Numeric. Minimum expression value for a gene to be
#'   considered active during conflict resolution. Default: \code{0}.
#' @param conflict_strategy Character. Conflict resolution order.
#'   \code{"biotype_first"} (default): select the best biotype tier first,
#'   then apply expression filtering within that tier. This is the more
#'   conservative default -- a silent protein-coding gene is preferred over
#'   a highly expressed lncRNA at the same locus.
#'   \code{"expression_first"}: apply expression filtering across all
#'   biotypes first, then pick the best biotype among expressed candidates.
#' @param co_dominance_ratio Numeric (0-1). In the expression tiebreaker step,
#'   genes with expression >= \code{co_dominance_ratio * max(expression)} in the
#'   group are retained together. Default: \code{0.1} (i.e. within one order of
#'   magnitude). Lower values (e.g. \code{0.01}) retain more co-dominant
#'   candidates; higher values (e.g. \code{0.5}) are more stringent.
#' @param biotype_order Character vector. Custom ordering of biotype
#'   categories for conflict resolution. Each element must be one of five
#'   keywords (matched case-insensitively against the \code{GENETYPE} column):
#'   \code{"protein"} (protein-coding),
#'   \code{"small_ncRNA"} (miRNA, snoRNA, snRNA, rRNA, scaRNA),
#'   \code{"antisense"},
#'   \code{"lncRNA"} (lncRNA and other ncRNA),
#'   \code{"pseudogene"}.
#'   The order determines priority: rank 1 = highest. Categories not listed
#'   keep their default relative order and are appended after the listed ones.
#'   Default: \code{c("protein", "small_ncRNA", "antisense", "lncRNA",
#'   "pseudogene")}. To prioritise lncRNAs over protein-coding genes while
#'   keeping everything else as-is, set \code{c("lncRNA", "protein")}.
#' @return The input data frame with \code{SYMBOL} and \code{annotation}
#'   columns resolved.
#' @importFrom GenomicRanges makeGRangesFromDataFrame findOverlaps
#' @importFrom GenomicFeatures genes promoters
#' @importFrom S4Vectors queryHits subjectHits
#' @export
resolve_gene_conflicts <- function(
  current_anno_df, txdb_obj, org_db_pkg,
  tss_region, gene_expr_map, min_expr = 0,
  conflict_strategy = c("biotype_first", "expression_first"),
  co_dominance_ratio = 0.1,
  biotype_order = c("protein", "small_ncRNA", "antisense", "lncRNA", "pseudogene")
) {
    if (nrow(current_anno_df) == 0) {
        return(current_anno_df)
    }
    conflict_strategy <- match.arg(conflict_strategy)

    # Build candidate gene map from promoter overlaps
    candidates <- .resolve_prepare_candidates(
        current_anno_df, txdb_obj, org_db_pkg,
        tss_region, gene_expr_map, min_expr,
        biotype_order = biotype_order
    )

    if (nrow(candidates) > 0) {
        # Apply conflict resolution strategy
        gene_map <- attr(candidates, "gene_map")
        resolved_candidates <- candidates %>%
            dplyr::left_join(gene_map, by = "gene_id") %>%
            dplyr::group_by(query_idx)

        if (conflict_strategy == "biotype_first") {
            resolved_candidates <- resolved_candidates %>%
                dplyr::filter(type_rank == min(type_rank, na.rm = TRUE)) %>%
                dplyr::mutate(has_active = any(is_active)) %>%
                dplyr::filter(!has_active | is_active)
        } else {
            resolved_candidates <- resolved_candidates %>%
                dplyr::mutate(has_active = any(is_active)) %>%
                dplyr::filter(!has_active | is_active) %>%
                dplyr::filter(type_rank == min(type_rank, na.rm = TRUE))
        }

        resolved_candidates <- resolved_candidates %>%
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
                            threshold <- max(max_tpm * co_dominance_ratio, 1e-6)
                            active_genes <- genes[tpms >= threshold]
                            if (length(active_genes) == 0) {
                                paste(sort(unique(genes)), collapse = ";")
                            } else {
                                paste(sort(unique(active_genes)), collapse = ";")
                            }
                        }
                    }
                }
            ) %>%
            dplyr::ungroup() %>%
            dplyr::filter(!is.na(final_symbols))

        if (!"SYMBOL" %in% colnames(current_anno_df)) {
            current_anno_df$SYMBOL <- NA_character_
        }
        if (!"annotation" %in% colnames(current_anno_df)) {
            current_anno_df$annotation <- NA_character_
        }

        match_idx <- match(
            resolved_candidates$query_idx,
            seq_len(nrow(current_anno_df))
        )
        valid_idx <- !is.na(match_idx)

        if (any(valid_idx)) {
            safe_match <- match_idx[valid_idx]
            current_anno_df$SYMBOL[safe_match] <-
                resolved_candidates$final_symbols[valid_idx]
            current_anno_df$annotation[safe_match] <- "Promoter"
        }
    }

    return(current_anno_df)
}

#' Internal: Prepare candidate gene map for conflict resolution
#'
#' Finds promoter overlaps, maps gene IDs to symbols and biotypes,
#' and scores expression activity.
#'
#' @return A data frame of candidates with attribute \code{gene_map}.
#' @keywords internal
#' @noRd
.resolve_prepare_candidates <- function(
    current_anno_df, txdb_obj, org_db_pkg,
    tss_region, gene_expr_map, min_expr,
    biotype_order = c("protein", "small_ncRNA", "antisense", "lncRNA", "pseudogene")
) {
    gr_input <- .with_known_upstream_noise_suppressed(
        GenomicRanges::makeGRangesFromDataFrame(current_anno_df,
            keep.extra.columns = TRUE
        )
    )
    all_genes <- .with_known_upstream_noise_suppressed(GenomicFeatures::genes(txdb_obj))
    promoters_gr <- .with_known_upstream_noise_suppressed(
        GenomicFeatures::promoters(all_genes,
            upstream = abs(tss_region[1]),
            downstream = abs(tss_region[2])
        )
    )
    hits <- .with_known_upstream_noise_suppressed(GenomicRanges::findOverlaps(
        gr_input, promoters_gr
    ))

    if (length(hits) == 0) {
        return(data.frame(
            query_idx = integer(), gene_id = character(),
            stringsAsFactors = FALSE
        ))
    }

    candidates <- data.frame(
        query_idx = S4Vectors::queryHits(hits),
        gene_id = .extract_txdb_gene_ids(all_genes)[S4Vectors::subjectHits(hits)],
        stringsAsFactors = FALSE
    )

    cols_to_get <- "SYMBOL"
    valid_cols <- AnnotationDbi::columns(.get_org_db_obj(org_db_pkg))
    has_genetype <- "GENETYPE" %in% valid_cols
    if (has_genetype) cols_to_get <- c(cols_to_get, "GENETYPE")

    gene_map <- .map_txdb_gene_ids(
        gene_ids = unique(candidates$gene_id),
        org_db = org_db_pkg,
        columns = cols_to_get,
        context = "resolve_gene_conflicts",
        warn = TRUE
    )

    gene_map$tpm <- if (!is.null(gene_expr_map)) {
        expr_upper <- setNames(gene_expr_map, toupper(names(gene_expr_map)))
        val <- expr_upper[toupper(gene_map$SYMBOL)]
        ifelse(is.na(val), 0, val)
    } else {
        0
    }
    gene_map$is_active <- if (min_expr == 0) {
        gene_map$tpm > 0
    } else {
        gene_map$tpm >= min_expr
    }

    if (has_genetype) {
        # Build biotype pattern lookup
        biotype_patterns <- list(
            protein     = "protein",
            small_ncRNA = "miRNA|snoRNA|snRNA|rRNA|scaRNA",
            antisense   = "antisense",
            lncRNA      = "lncRNA|ncrna",
            pseudogene  = "pseudo"
        )
        default_order <- c("protein", "small_ncRNA", "antisense", "lncRNA", "pseudogene")
        # Listed categories: ranks 1, 2, 3, ...
        # Unlisted categories: keep default relative order, appended after listed
        unlisted <- setdiff(default_order, biotype_order)
        full_order <- c(biotype_order, unlisted)
        gene_map$type_rank <- length(full_order) + 1L  # fallback for unknown biotypes
        for (i in seq_along(full_order)) {
            key <- full_order[[i]]
            pattern <- biotype_patterns[[key]]
            if (!is.null(pattern)) {
                gene_map <- gene_map %>%
                    dplyr::mutate(
                        type_rank = dplyr::if_else(
                            grepl(pattern, GENETYPE, ignore.case = TRUE),
                            as.integer(i), type_rank
                        )
                    )
            }
        }
    } else {
        gene_map$type_rank <- 1
    }

    attr(candidates, "gene_map") <- gene_map
    candidates
}


#' Internal: Check if anchor type is promoter-like (P or dual)
#' @keywords internal
#' @noRd
.is_promoter_like <- function(t) t %in% c("P", "dual")

#' Internal: Check if anchor type is distal/reenhancer-like (E, eP, eG, dual)
#' @keywords internal
#' @noRd
.is_distal_like <- function(t) t %in% c("E", "eP", "eG", "dual")

#' Internal: Check if anchor type is gene-body-like (G or eG)
#' @keywords internal
#' @noRd
.is_gene_body_like <- function(t) t %in% c("G", "eG")

#' Internal: Reclassify Anchor by Expression
#'
#' Given an anchor's gene symbol and type, filters to active genes (present in
#' \code{allow}) and optionally reclassifies silent promoters/gene bodies.
#'
#' @param g Character. Semicolon-delimited gene string.
#' @param t Character. Anchor type code (P, E, G, eP, eG).
#' @param allow Character vector. Whitelist of active gene symbols.
#' @param down Logical. If \code{TRUE}, reclassify silent P->eP and G->eG.
#' @return A list with \code{type} and \code{gene}.
#' @keywords internal
#' @noRd
clean_anchor <- function(g, t, allow, down) {
    g_char <- as.character(g)
    t_char <- as.character(t)
    if (is.na(g_char) || g_char == "") {
        return(list(type = t_char, gene = NA_character_))
    }
    gs <- unlist(strsplit(g_char, ";"))
    active_gs <- trimws(gs[toupper(trimws(gs)) %in% toupper(allow)])
    if (length(active_gs) > 0) {
        return(list(type = t_char, gene = paste(unique(active_gs), collapse = ";")))
    }
    if (down) {
        new_type <- dplyr::case_when(
            t_char == "P" ~ "eP", t_char == "G" ~ "eG", TRUE ~ t_char
        )
        return(list(type = new_type, gene = NA_character_))
    }
    return(list(type = t_char, gene = NA_character_))
}



#' Internal: Compute Dominant Interaction Type
#'
#' Returns the most frequent non-NA value in a vector.
#'
#' @param x Character vector of interaction type codes.
#' @return The modal value, or \code{NA_character_} if empty.
#' @keywords internal
#' @noRd
.get_dom <- function(x) {
    # Returns the modal (most frequent) non-NA value. When multiple values
    # tie for the highest frequency, the first in alphabetical order is
    # returned (which.max on a named table breaks ties by position).
    # This is acceptable for dominant-interaction labelling but worth noting
    # for small datasets where ties are more likely.
    x <- x[!is.na(x)]
    if (length(x) == 0) {
        return(NA_character_)
    }
    names(which.max(table(x)))
}

#' Internal: Look Up Per-Gene Mean Expression
#'
#' @param g Character. Gene symbol.
#' @param vals Named numeric vector of per-gene mean expression.
#' @return Numeric expression value, or 0 if missing.
#' @keywords internal
#' @noRd
.get_expr <- function(g, vals) {
    e <- vals[g]
    e[is.na(e)] <- 0
    return(e)
}

#' Internal: Compute Raw Promoter-Centric Statistics
#'
#' Builds a per-gene summary from promoter-anchored loop rows.
#'
#' @param loop_df Loop annotation data frame.
#' @return A data frame with columns \code{Gene}, \code{Total_Loops_Filtered},
#'   \code{n_Linked_Promoters_Filtered}, \code{n_Linked_Distal_Filtered},
#'   \code{Dominant_Interaction_Filtered}.
#' @keywords internal
#' @noRd
.compute_raw_promoter_stats <- function(loop_df) {
    raw_stats_df <- dplyr::bind_rows(
        loop_df %>% dplyr::filter(.is_promoter_like(anchor1_type) & !is.na(anchor1_gene)) %>%
            dplyr::select(
                Gene = anchor1_gene, Neighbor_Type = anchor2_type,
                Loop_Type = loop_type
            ) %>% dplyr::mutate(Gene = as.character(Gene)),
        loop_df %>% dplyr::filter(.is_promoter_like(anchor2_type) & !is.na(anchor2_gene)) %>%
            dplyr::select(
                Gene = anchor2_gene, Neighbor_Type = anchor1_type,
                Loop_Type = loop_type
            ) %>% dplyr::mutate(Gene = as.character(Gene))
    ) %>%
        tidyr::separate_rows(Gene, sep = ";") %>%
        dplyr::mutate(Gene = trimws(Gene)) %>%
        dplyr::filter(Gene != "" & !is.na(Gene)) %>%
        dplyr::group_by(Gene) %>%
        dplyr::summarise(
            Total_Loops_Filtered = dplyr::n(),
            n_Linked_Promoters_Filtered = sum(.is_promoter_like(Neighbor_Type), na.rm = TRUE),
            n_Linked_Distal_Filtered = sum(
                (.is_distal_like(Neighbor_Type) | .is_gene_body_like(Neighbor_Type)),
                na.rm = TRUE
            ),
            Dominant_Interaction_Filtered = .get_dom(Loop_Type),
            .groups = "drop"
        )
    raw_stats_df
}

#' Internal: Build Promoter-Centric Summary Data Frame
#'
#' Merges raw refined stats with upstream promoter stats, appends
#' expression and connectivity classification columns.
#'
#' @param raw_stats_df Raw per-gene summary from
#'   \code{\link{.compute_raw_promoter_stats}}.
#' @param upstream_promoter_stats Upstream promoter stats (or NULL).
#' @param vals Named numeric vector of per-gene mean expression.
#' @param threshold Numeric. Expression threshold.
#' @param hub_percentile Numeric. Quantile for hub cutoff.
#' @return A data frame of promoter-centric statistics, or \code{NULL}.
#' @keywords internal
#' @noRd
.build_promoter_centric_df <- function(
  raw_stats_df, upstream_promoter_stats,
  vals, threshold, hub_percentile
) {
    empty_promoter_df <- data.frame(
        Gene = character(), Total_Loops = integer(), n_Linked_Promoters = integer(),
        n_Linked_Distal = integer(), Dominant_Interaction = character(),
        Is_High_Connectivity_Gene = character(), Is_High_Distal_Connectivity_Gene = character(),
        Is_Active_Gene = character(), stringsAsFactors = FALSE
    )
    if (nrow(raw_stats_df) == 0) {
        return(empty_promoter_df)
    }
    final_cutoff <- max(stats::quantile(
        raw_stats_df$Total_Loops_Filtered, hub_percentile,
        na.rm = TRUE
    ), 3)
    distal_cutoff <- max(stats::quantile(
        raw_stats_df$n_Linked_Distal_Filtered, hub_percentile,
        na.rm = TRUE
    ), 2)

    if (!is.null(upstream_promoter_stats)) {
        promoter_centric_df <- upstream_promoter_stats %>%
            dplyr::left_join(raw_stats_df, by = "Gene") %>%
            dplyr::mutate(
                Total_Loops = dplyr::coalesce(Total_Loops_Filtered, 0),
                n_Linked_Promoters = dplyr::coalesce(n_Linked_Promoters_Filtered, 0),
                n_Linked_Distal = dplyr::coalesce(n_Linked_Distal_Filtered, 0),
                Dominant_Interaction = dplyr::coalesce(
                    Dominant_Interaction_Filtered, "None"
                ),
                Mean_Expression_Temp = .get_expr(Gene, vals),
                Is_Active_Gene = dplyr::if_else(
                    Mean_Expression_Temp >= threshold, "Yes", "No"
                ),
                Is_High_Connectivity_Gene = dplyr::if_else(
                    Total_Loops >= final_cutoff, "Yes", "No"
                ),
                Is_High_Distal_Connectivity_Gene = dplyr::if_else(
                    n_Linked_Distal >= distal_cutoff, "Yes", "No"
                )
            ) %>%
            dplyr::select(
                Gene, Total_Loops, n_Linked_Promoters, n_Linked_Distal,
                Dominant_Interaction, Is_High_Connectivity_Gene,
                Is_High_Distal_Connectivity_Gene, Is_Active_Gene,
                dplyr::everything()
            ) %>%
            dplyr::select(-any_of(c(
                "Total_Loops_Filtered",
                "n_Linked_Promoters_Filtered", "n_Linked_Distal_Filtered",
                "Dominant_Interaction_Filtered", "Is_Regulatory_Hub",
                "Mean_Expression_Temp", "n_Linked_Enhancers", "n_Linked_GeneBodies"
            )))
    } else {
        promoter_centric_df <- raw_stats_df %>%
            dplyr::rename(
                Total_Loops = Total_Loops_Filtered,
                n_Linked_Promoters = n_Linked_Promoters_Filtered,
                n_Linked_Distal = n_Linked_Distal_Filtered,
                Dominant_Interaction = Dominant_Interaction_Filtered
            ) %>%
            dplyr::mutate(
                Mean_Expression_Temp = .get_expr(Gene, vals),
                Is_Active_Gene = dplyr::if_else(
                    Mean_Expression_Temp >= threshold, "Yes", "No"
                ),
                Is_High_Connectivity_Gene = dplyr::if_else(
                    Total_Loops >= final_cutoff, "Yes", "No"
                ),
                Is_High_Distal_Connectivity_Gene = dplyr::if_else(
                    n_Linked_Distal >= distal_cutoff, "Yes", "No"
                )
            ) %>%
            dplyr::select(
                Gene, Total_Loops, n_Linked_Promoters, n_Linked_Distal,
                Dominant_Interaction, Is_High_Connectivity_Gene,
                Is_High_Distal_Connectivity_Gene, Is_Active_Gene,
                dplyr::everything()
            ) %>%
            dplyr::select(-any_of("Mean_Expression_Temp"))
    }
    promoter_centric_df <- promoter_centric_df %>%
        dplyr::arrange(dplyr::desc(n_Linked_Distal))
    promoter_centric_df
}

#' Internal: Build Distal Element Connectivity Data Frame
#'
#' @param loop_df Loop annotation data frame with anchor-level columns.
#' @param hub_percentile Numeric. Quantile for hub cutoff.
#' @return A data frame of distal element statistics, or \code{NULL}.
#' @keywords internal
#' @noRd
.build_distal_element_df <- function(loop_df, hub_percentile) {
    if (!"a1_id" %in% colnames(loop_df)) {
        return(NULL)
    }
    distal_raw_df <- dplyr::bind_rows(
        loop_df %>% dplyr::filter((.is_distal_like(anchor1_type) | .is_gene_body_like(anchor1_type))) %>%
            dplyr::select(
                Distal_Anchor_ID = a1_id, Neighbor_Type = anchor2_type,
                Loop_Type = loop_type, Neighbor_Gene = anchor2_gene
            ),
        loop_df %>% dplyr::filter((.is_distal_like(anchor2_type) | .is_gene_body_like(anchor2_type))) %>%
            dplyr::select(
                Distal_Anchor_ID = a2_id, Neighbor_Type = anchor1_type,
                Loop_Type = loop_type, Neighbor_Gene = anchor1_gene
            )
    ) %>%
        dplyr::group_by(Distal_Anchor_ID) %>%
        dplyr::summarise(
            Total_Loops_Filtered = dplyr::n(),
            n_Linked_Distal_Filtered = sum(
                (.is_distal_like(Neighbor_Type) | .is_gene_body_like(Neighbor_Type)),
                na.rm = TRUE
            ),
            n_Linked_Promoters_Filtered = sum(.is_promoter_like(Neighbor_Type), na.rm = TRUE),
            Dominant_Interaction_Filtered = .get_dom(Loop_Type),
            Target_Genes_Filtered = extract_genes(
                Neighbor_Gene[.is_promoter_like(Neighbor_Type)]
            ),
            .groups = "drop"
        )

    anchor_map <- dplyr::bind_rows(
        loop_df %>% dplyr::select(
            anchor_id = a1_id, chr = chr1,
            start = start1, end = end1, cluster_id
        ),
        loop_df %>% dplyr::select(
            anchor_id = a2_id, chr = chr2,
            start = start2, end = end2, cluster_id
        )
    ) %>% dplyr::distinct()

    if (nrow(distal_raw_df) == 0) {
        return(NULL)
    }
    final_cutoff_dist <- max(stats::quantile(
        distal_raw_df$Total_Loops_Filtered, hub_percentile,
        na.rm = TRUE
    ), 3)
    temp_df <- distal_raw_df %>%
        dplyr::rename(
            Total_Loops = Total_Loops_Filtered,
            n_Linked_Distal = n_Linked_Distal_Filtered,
            n_Linked_Promoters = n_Linked_Promoters_Filtered,
            Dominant_Interaction = Dominant_Interaction_Filtered,
            Target_Genes = Target_Genes_Filtered
        ) %>%
        dplyr::mutate(
            Is_High_Connectivity_Distal_Element = dplyr::if_else(
                Total_Loops >= final_cutoff_dist, "Yes", "No"
            )
        )
    temp_df <- temp_df %>%
        dplyr::select(-any_of(c(
            "chr", "start", "end", "cluster_id",
            "Distal_Type", "Distal_Type_Filtered", "Total_Loops_Filtered",
            "Target_Genes_Filtered", "n_Linked_Distal_Filtered",
            "n_Linked_Promoters_Filtered", "Dominant_Interaction_Filtered"
        )))
    distal_element_df <- temp_df %>%
        dplyr::left_join(anchor_map,
            by = c("Distal_Anchor_ID" = "anchor_id")
        ) %>%
        dplyr::select(
            chr, start, end, cluster_id, Total_Loops,
            n_Linked_Promoters, n_Linked_Distal, Dominant_Interaction,
            any_of("Is_High_Connectivity_Distal_Element"), Target_Genes
        ) %>%
        dplyr::filter(Total_Loops > 0) %>%
        dplyr::arrange(dplyr::desc(Total_Loops))
    distal_element_df
}


#' Internal: Compute Refined Network Statistics
#'
#' Recalculates promoter-centric and distal-element connectivity statistics
#' after expression-aware filtering, merging with upstream annotation stats
#' where available.
#'
#' @param loop_df Loop annotation data frame after expression filtering.
#' @param upstream_promoter_stats Upstream promoter-centric stats (or NULL).
#' @param vals Named numeric vector of per-gene mean expression.
#' @param threshold Numeric. Expression threshold for active gene classification.
#' @param hub_percentile Numeric. Quantile for hub cutoff.
#' @return A list with \code{promoter_centric} and \code{distal_element}
#'   data frames.
#' @importFrom stats quantile
#' @keywords internal
#' @noRd
compute_refined_stats <- function(
  loop_df, upstream_promoter_stats,
  vals, threshold, hub_percentile
) {
    raw_stats_df <- .compute_raw_promoter_stats(loop_df)
    promoter_centric_df <- .build_promoter_centric_df(
        raw_stats_df, upstream_promoter_stats,
        vals, threshold, hub_percentile
    )
    distal_element_df <- .build_distal_element_df(loop_df, hub_percentile)

    list(
        promoter_centric = promoter_centric_df,
        distal_element = distal_element_df
    )
}


#' Internal: Load Expression Matrix
#'
#' Reads a normalized expression matrix (TPM/FPKM), sets gene identifiers as
#' the first column, validates the requested sample columns, and returns
#' per-gene mean expression values. Sample column names must be unique; missing
#' or duplicated selections raise an error instead of being silently dropped.
#'
#' @param expr_matrix_file Character. Path to the expression matrix file.
#' @param sample_columns Character or integer vector. Sample columns to average.
#'   Character selections must exactly match unique column names.
#' @return Named numeric vector of per-gene mean expression values.
#' @importFrom data.table fread
#' @keywords internal
#' @noRd
load_expression_matrix <- function(expr_matrix_file, sample_columns = NULL) {
    if (!file.exists(expr_matrix_file)) {
        stop("Expression matrix file not found: ", expr_matrix_file)
    }
    d <- as.data.frame(data.table::fread(
        expr_matrix_file,
        data.table = FALSE,
        showProgress = FALSE
    ))
    if (ncol(d) < 2) {
        stop(
            "Expression matrix must contain a gene identifier column and at least one sample column."
        )
    }

    gene_ids <- trimws(as.character(d[[1]]))

    # Diagnostic: check for potential case mismatch with OrgDb conventions
    n_nonempty <- sum(nzchar(gene_ids))
    if (n_nonempty > 0) {
        frac_upper <- mean(grepl("^[A-Z0-9.]+$", gene_ids[nzchar(gene_ids)]))
        if (frac_upper < 0.5) {
            warning(
                "Only ", round(frac_upper * 100), "% of expression matrix gene identifiers ",
                "are all-uppercase. Matching is case-insensitive (toupper() on both sides), ",
                "so this will NOT cause misclassification. However, if your expression ",
                "matrix uses a different identifier convention (e.g., Ensembl IDs while ",
                "OrgDb returns SYMBOL), cross-check that genes map correctly. ",
                "Human OrgDb returns UPPERCASE (e.g., 'TP53'); mouse OrgDb returns ",
                "Title Case (e.g., 'Trp53').",
                call. = FALSE
            )
        }
    }

    dup_genes <- unique(gene_ids[duplicated(gene_ids) & nzchar(gene_ids)])
    if (length(dup_genes) > 0) {
        warning(
            "Expression matrix contains ", length(dup_genes),
            " duplicated gene identifier(s). Only the first occurrence of each duplicated gene is retained. ",
            "Consider aggregating duplicates before calling looplook.",
            call. = FALSE
        )
    }
    sample_names <- colnames(d)[-1]
    if (any(is.na(sample_names)) || any(!nzchar(sample_names))) {
        stop("Expression matrix contains empty sample column names.")
    }
    dup_sample_names <- unique(sample_names[duplicated(sample_names)])
    if (length(dup_sample_names) > 0) {
        stop(
            "Expression matrix contains duplicated sample column names: ",
            paste(dup_sample_names, collapse = ", "),
            ". Rename columns uniquely before calling looplook."
        )
    }
    d <- d[, -1, drop = FALSE]
    colnames(d) <- sample_names

    if (is.null(sample_columns)) {
        sub_mat <- d
    } else if (is.character(sample_columns)) {
        dup_requested <- unique(sample_columns[duplicated(sample_columns)])
        if (length(dup_requested) > 0) {
            stop(
                "`sample_columns` contains duplicates: ",
                paste(dup_requested, collapse = ", "),
                "."
            )
        }

        missing_cols <- setdiff(sample_columns, sample_names)
        if (length(missing_cols) > 0) {
            stop(
                "Requested sample columns not found in expression matrix: ",
                paste(missing_cols, collapse = ", ")
            )
        }
        sub_mat <- d[, sample_columns, drop = FALSE]
    } else {
        sample_columns <- as.integer(sample_columns)
        dup_requested <- unique(sample_columns[duplicated(sample_columns)])
        if (length(dup_requested) > 0) {
            stop(
                "`sample_columns` contains duplicated column indices: ",
                paste(dup_requested, collapse = ", "),
                "."
            )
        }
        if (any(is.na(sample_columns)) || any(sample_columns < 1L | sample_columns > ncol(d))) {
            stop("`sample_columns` contains invalid column indices for the expression matrix.")
        }
        sub_mat <- d[, sample_columns, drop = FALSE]
    }
    if (ncol(sub_mat) == 0) stop("No valid sample columns found in expression matrix.")

    numeric_token_pattern <- "^[+-]?(?:Inf|NaN|(?:\\d+(?:\\.\\d*)?|\\.\\d+)(?:[eE][+-]?\\d+)?)$"
    sub_mat_parsed <- lapply(sub_mat, function(x) {
        raw_x <- trimws(as.character(x))
        is_missing <- is.na(raw_x) | raw_x == ""
        is_numeric <- grepl(numeric_token_pattern, raw_x)
        bad <- !is_missing & !is_numeric
        values <- rep(NA_real_, length(raw_x))
        if (any(is_numeric)) {
            values[is_numeric] <- as.numeric(raw_x[is_numeric])
        }
        list(values = values, bad = bad)
    })
    sub_mat_num <- lapply(sub_mat_parsed, `[[`, "values")
    bad_cols <- names(sub_mat)[vapply(sub_mat_parsed, function(x) any(x$bad), logical(1))]
    if (length(bad_cols) > 0) {
        stop(
            "Expression matrix contains non-numeric values in sample columns: ",
            paste(bad_cols, collapse = ", ")
        )
    }
    sub_mat_num <- as.data.frame(sub_mat_num, check.names = FALSE)

    if (is.null(sample_columns)) {
        message("sample_columns not specified; averaging all ", ncol(sub_mat_num),
                " sample columns for baseline expression")
    }
    vals <- if (ncol(sub_mat_num) > 1) {
        rowMeans(sub_mat_num, na.rm = TRUE)
    } else {
        sub_mat_num[[1]]
    }
    vals[is.nan(vals)] <- NA_real_
    names(vals) <- gene_ids
    vals
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
#' @noRd
get_colors <- function(n, palette_input) {
    if (n <= 0) {
        return(character(0))
    }
    safe_n <- max(1, n)

    if (length(palette_input) == 1 && palette_input %in% row.names(RColorBrewer::brewer.pal.info)) {
        max_avail <- RColorBrewer::brewer.pal.info[palette_input, "maxcolors"]
        pal_n <- min(max(3, safe_n), max_avail)
        pal <- RColorBrewer::brewer.pal(pal_n, palette_input)
        cols <- grDevices::colorRampPalette(pal)(safe_n)
        return(cols)
    } else if (length(palette_input) >= 1) {
        if (length(palette_input) < safe_n) {
            cols <- rep(palette_input, length.out = safe_n)
        } else {
            cols <- palette_input[seq_len(safe_n)]
        }
        return(cols)
    } else {
        cols <- scales::hue_pal()(safe_n)
        return(cols)
    }
}


#' Internal: Draw Karyotype Heatmap
#'
#' Creates a genome-wide heatmap of genomic feature density (e.g., loops) across
#' standard chromosomes, binned by a fixed window size, and rendered as a
#' karyotype plot. Only canonical chromosomes (chr1-22,X,Y for human;
#' chr1-19,X,Y for mouse) are displayed; chrM and alternate contigs are
#' silently excluded.
#'
#' @param gr_data (GRanges) Genomic ranges to visualize (e.g., loop anchors).
#' @param title_prefix (character) Subtitle descriptor (e.g., sample name).
#' @param bin_size (integer) Bin width in base pairs (e.g., 1e6 for 1 Mb).
#' @param sat_level (numeric) Quantile (0-1) for color saturation (e.g., 0.95).
#' @param ref_txdb (TxDb or similar) Reference genome annotation for chromosome lengths.
#' @param plot_species (character) Genome build/species code (e.g., "hg38", "mm10").
#' @param unit_label (character) Unit for load annotation (e.g., "loops").
#' @param custom_colors (character vector) Optional custom color palette.
#' @keywords internal
#' @importFrom GenomeInfoDb seqinfo seqlevelsStyle keepSeqlevels seqlengths seqlevels
#' @importFrom GenomicRanges GRanges tileGenome countOverlaps
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors mcols
#' @importFrom grDevices colorRampPalette
#' @importFrom karyoploteR plotKaryotype kpRect getDefaultPlotParams
#' @importFrom fields image.plot
#' @return A \code{looplook_karyo} object wrapping a rendered PNG payload. Use
#'   \code{print()} to display.
#' @noRd
draw_karyo_heatmap_internal <- function(gr_data, title_prefix, bin_size, sat_level, ref_txdb, plot_species, unit_label, custom_colors = NULL) {
    standard_chroms <- paste0("chr", c(seq_len(22), "X", "Y"))
    if (grepl("mm", plot_species, fixed = TRUE)) standard_chroms <- paste0("chr", c(seq_len(19), "X", "Y"))

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
        full_genome_gr <- .with_known_upstream_noise_suppressed(
            GenomicRanges::GRanges(seqnames = valid_chroms, ranges = IRanges::IRanges(start = 1, end = GenomeInfoDb::seqlengths(std_seqinfo)[valid_chroms]))
        )
        GenomeInfoDb::seqinfo(full_genome_gr) <- std_seqinfo[valid_chroms]
        tiles <- .with_known_upstream_noise_suppressed(
            GenomicRanges::tileGenome(GenomeInfoDb::seqinfo(full_genome_gr), tilewidth = bin_size, cut.last.tile.in.chrom = TRUE)
        )
        hits <- GenomicRanges::countOverlaps(tiles, gr_data)

        bin_size_mb <- bin_size / 1e6
        median_val <- median(hits[hits > 0], na.rm = TRUE)
        if (is.na(median_val)) median_val <- 0

        heatmap_colors <- if (is.null(custom_colors)) c("#FFFFFF", "#FFFFCC", "#FFEDA0", "#FED976", "#FEB24C", "#FD8D3C", "#FC4E2A", "#E31A1C", "#BD0026", "#800026", "#000000") else custom_colors

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

        # Render once to PNG and keep the bytes in-memory so deferred report
        # rendering does not depend on a temp file surviving until a later chunk.
        f <- tempfile(fileext = ".png")
        grDevices::png(f, width = 10, height = 8, units = "in", res = 150)
        needs_close <- TRUE
        on.exit(if (needs_close) try(grDevices::dev.off(), silent = TRUE), add = TRUE)

        graphics::par(oma = c(2, 2, 6, 2))
        pp <- karyoploteR::getDefaultPlotParams(plot.type = 1)
        pp$leftmargin <- 0.08
        pp$rightmargin <- 0.08
        pp$data1height <- 100
        kp <- .with_known_upstream_noise_suppressed(
            karyoploteR::plotKaryotype(genome = plot_species, plot.type = 1, chromosomes = valid_chroms, plot.params = pp, main = NULL)
        )
        karyoploteR::kpRect(kp, data = tiles, y0 = 0, y1 = 1, col = S4Vectors::mcols(tiles)$color, border = NA)
        main_title <- paste0("Loop Analysis: ", title_prefix, "\n(Genomic Load: Median ~", round(median_val / bin_size_mb, 1), " ", unit_label, "/MB)")
        graphics::mtext(main_title, side = 3, line = 1, outer = TRUE, cex = 1.2, font = 2)
        fields::image.plot(
            legend.only = TRUE, zlim = c(0, max_load), col = cols,
            legend.lab = paste0("Load (", unit_label, "/MB)"), legend.mar = 4.5, smallplot = c(0.88, 0.91, 0.3, 0.7)
        )

        needs_close <- FALSE
        grDevices::dev.off()
        png_raw <- NULL
        if (file.exists(f)) {
            png_raw <- tryCatch(
                readBin(f, what = "raw", n = file.info(f)$size),
                error = function(e) NULL
            )
            if (!is.null(png_raw)) {
                unlink(f)
            }
        }

        payload <- list(type = "karyo_heatmap")
        if (!is.null(png_raw)) {
            payload$png_raw <- png_raw
        } else if (file.exists(f)) {
            payload$file <- f
        }

        structure(payload, class = "looplook_karyo")
    }
}

#' Print looplook karyogram
#'
#' Displays a previously captured karyotype heatmap. Uses embedded PNG bytes
#' when available, otherwise falls back to a file-backed image.
#'
#' @param x A \code{looplook_karyo} object.
#' @param ... Additional arguments (unused).
#' @return The input \code{looplook_karyo} object, returned invisibly after
#'   drawing the image.
#' @export
print.looplook_karyo <- function(x, ...) {
    if (!is.null(x$png_raw)) {
        f <- tempfile(fileext = ".png")
        writeBin(x$png_raw, f)
        on.exit(unlink(f), add = TRUE)
        if (requireNamespace("png", quietly = TRUE)) {
            img <- png::readPNG(f)
            grid::grid.newpage()
            grid::grid.raster(img, interpolate = TRUE)
        } else {
            utils::browseURL(f)
        }
        return(invisible(x))
    }

    if (is.null(x$file) || !file.exists(x$file)) {
        stop("Karyo image is unavailable.", call. = FALSE)
    }
    if (requireNamespace("png", quietly = TRUE)) {
        img <- png::readPNG(x$file)
        grid::grid.newpage()
        grid::grid.raster(img, interpolate = TRUE)
    } else {
        utils::browseURL(x$file)
    }
    invisible(x)
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
#' @noRd
draw_circular_bar_plot <- function(data_df, project_name, filename = NULL, color_vec) {
    color_vec <- as.character(color_vec)
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
        ggplot2::geom_col(ggplot2::aes(y = Unique_Gene_Count), width = 0.8, color = "white", linewidth = 0.2) +
        ggplot2::geom_text(ggplot2::aes(y = Unique_Gene_Count + max_gene_count * 0.02, label = Label_Text), hjust = 0, size = 3.5, fontface = "bold") +
        ggplot2::coord_polar(theta = "y", start = 0, clip = "off") +
        ggplot2::scale_fill_manual(values = color_vec, name = "Loop Type") +
        ggplot2::theme_minimal(base_size = 14) +
        ggplot2::theme(axis.title = ggplot2::element_blank(), axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), panel.grid = ggplot2::element_blank(), legend.position = "right", plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"), plot.margin = ggplot2::margin(t = 20, r = 100, b = 20, l = 20, unit = "pt")) +
        ggplot2::scale_y_continuous(limits = c(-max_gene_count * 0.4, max_gene_count * 1.3)) +
        ggplot2::labs(title = paste0(project_name, ": Unique Target Genes (Ascending)"))

    return(p)
}

#' Internal: Simplify Genomic Annotation to Broad Categories
#'
#' Collapses detailed ChIPseeker annotation strings into five broad
#' categories: Promoter, Intron, Exon, Distal Intergenic, and Downstream.
#' Anything unrecognised is labelled "Others".
#'
#' @param x Character vector of annotation strings.
#' @return Character vector of simplified categories.
#' @keywords internal
#' @noRd
simplify_annotation <- function(x) {
    vapply(x, function(s) {
        if (grepl("Promoter", s, ignore.case = TRUE)) {
            return("Promoter")
        }
        if (grepl("Intron", s, ignore.case = TRUE)) {
            return("Intron")
        }
        if (grepl("Exon", s, ignore.case = TRUE)) {
            return("Exon")
        }
        if (grepl("Intergenic", s, ignore.case = TRUE)) {
            return("Distal Intergenic")
        }
        if (grepl("Downstream", s, ignore.case = TRUE)) {
            return("Downstream")
        }
        return("Others")
    }, FUN.VALUE = character(1))
}


#' Internal: Draw Pie Chart with Outside Labels
#'
#' Simplified pie chart with labels placed outside the slices, using
#' RColorBrewer palettes for genomic annotation categories.
#'
#' @param data_df Data frame with an annotation column.
#' @param group_col Character. Column name for grouping.
#' @param title Character. Plot title.
#' @param palette Character. RColorBrewer palette name.
#' @return A ggplot object, or NULL if data is empty.
#' @importFrom ggplot2 ggplot aes geom_bar geom_segment geom_text coord_polar
#'   xlim scale_fill_brewer theme_void labs theme element_text
#' @keywords internal
#' @noRd
draw_pie_with_outside_labels <- function(data_df, group_col, title, palette) {
    if (is.null(data_df) || nrow(data_df) == 0) {
        return(NULL)
    }

    plot_data <- data_df
    plot_data$Simplified <- simplify_annotation(plot_data[[group_col]])

    stats <- plot_data %>%
        dplyr::group_by(Simplified) %>%
        dplyr::summarise(Count = dplyr::n(), .groups = "drop") %>%
        dplyr::mutate(
            Fraction = Count / sum(Count),
            LabelText = ifelse(Fraction >= 0.01, paste0(Count, " (", round(Fraction * 100, 1), "%)"), "")
        ) %>%
        dplyr::arrange(dplyr::desc(Simplified))

    if (nrow(stats) == 0) {
        return(NULL)
    }

    stats <- stats %>%
        dplyr::mutate(
            ymax = cumsum(Fraction),
            ymin = c(0, head(ymax, n = -1)),
            ypos = (ymax + ymin) / 2,
            hjust = ifelse(ypos < 0.5, 0, 1)
        )

    ggplot2::ggplot(stats, ggplot2::aes(y = Fraction, fill = Simplified)) +
        ggplot2::geom_bar(ggplot2::aes(x = 1), width = 1, stat = "identity", color = "white") +
        ggplot2::geom_segment(
            data = subset(stats, LabelText != ""),
            ggplot2::aes(x = 1.51, xend = 1.62, y = ypos, yend = ypos),
            color = "grey50", linewidth = 0.5
        ) +
        ggplot2::geom_text(
            data = subset(stats, LabelText != ""),
            ggplot2::aes(x = 1.65, y = ypos, label = LabelText, hjust = hjust),
            size = 3.5, fontface = "bold", check_overlap = FALSE
        ) +
        ggplot2::coord_polar("y", start = 0) +
        ggplot2::xlim(0.5, 2.5) +
        ggplot2::scale_fill_brewer(palette = palette, name = "Genomic Feature") +
        ggplot2::theme_void() +
        ggplot2::labs(title = title) +
        ggplot2::theme(
            legend.position = "bottom",
            plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 12)
        )
}

#' Internal: Record Versions of Key External Databases
#'
#' Attempts to look up the installed versions of annotation packages
#' commonly used by looplook. Missing packages are recorded as \code{NA}.
#' @param species Character. Genome assembly (e.g. \code{"hg38"}).
#' @return A named list of version strings.
#' @keywords internal
#' @noRd
.record_database_versions <- function(species = NULL) {
    get_pkg_version <- function(pkg) {
        if (!is.null(pkg) && requireNamespace(pkg, quietly = TRUE))
            as.character(utils::packageVersion(pkg))
        else
            NA_character_
    }
    txdb_pkg <- if (!is.null(species)) tryCatch(species_txdb_pkg(species), error = function(e) NULL) else NULL
    orgdb_pkg <- if (!is.null(species)) tryCatch(species_orgdb_pkg(species), error = function(e) NULL) else NULL
    bsgenome_pkg <- if (!is.null(species)) tryCatch(species_bsgenome_pkg(species), error = function(e) NULL) else NULL

    list(
        TxDb       = get_pkg_version(txdb_pkg),
        OrgDb      = get_pkg_version(orgdb_pkg),
        BSgenome   = get_pkg_version(bsgenome_pkg),
        JASPAR     = get_pkg_version("JASPAR2020"),
        clusterProfiler = get_pkg_version("clusterProfiler"),
        STRINGdb   = get_pkg_version("STRINGdb"),
        motifmatchr = get_pkg_version("motifmatchr"),
        txdb_pkg   = txdb_pkg,
        orgdb_pkg  = orgdb_pkg,
        bsgenome_pkg = bsgenome_pkg
    )
}

#' Internal: Build Run Metadata
#'
#' Creates a standardised metadata list recording package version,
#' call timestamp, key parameters, and genome build for provenance tracking.
#'
#' @param fun Character. Function name.
#' @param params Named list of relevant parameters.
#' @param genome_build Character or NULL. Genome assembly if applicable.
#' @param score_semantics Character or NULL. Score interpretation note.
#' @return A list of metadata.
#' @keywords internal
#' @noRd
.build_looplook_metadata <- function(fun, params = list(),
                                      genome_build = NULL,
                                      score_semantics = NULL,
                                      diagnostics = NULL,
                                      database_versions = NULL) {
    m <- list(
        package = "looplook",
        version = as.character(utils::packageVersion("looplook")),
        function_name = fun,
        call_time = Sys.time(),
        parameters = params,
        genome_build = genome_build,
        score_semantics = score_semantics,
        r_version = R.version.string,
        platform = R.version$platform
    )
    if (!is.null(diagnostics)) m$diagnostics <- diagnostics
    if (!is.null(database_versions)) m$database_versions <- database_versions
    m
}
