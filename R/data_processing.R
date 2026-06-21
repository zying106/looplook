#' Internal: Validate and Normalize BEDPE Data Frame
#'
#' Shared validation logic for BEDPE input. Checks column count, converts
#' coordinates to numeric, validates start < end, warns on narrow/wide anchors,
#' normalizes anchor order, and converts 0-based BEDPE to 1-based coordinates.
#'
#' @param df A data frame (from \code{data.table::fread}).
#' @param quiet Logical. Suppress data-quality warnings.
#' @return A validated, normalized data frame with 1-based coordinates.
#' @keywords internal
#' @noRd
.validate_bedpe_df <- function(df, quiet = FALSE) {
    if (ncol(df) < 6) {
        stop("BEDPE file must have at least 6 columns.", call. = FALSE)
    }

    coord_cols <- c(2, 3, 5, 6)
    for (cc in coord_cols) {
        df[[cc]] <- suppressWarnings(as.numeric(df[[cc]]))
    }
    coord_mat <- cbind(df[[2]], df[[3]], df[[5]], df[[6]])
    if (any(is.na(coord_mat))) {
        stop("BEDPE file contains non-numeric coordinate columns.", call. = FALSE)
    }
    if (any(df[[2]] >= df[[3]], na.rm = TRUE) || any(df[[5]] >= df[[6]], na.rm = TRUE)) {
        stop("BEDPE file contains rows with start >= end (zero-width or invalid).", call. = FALSE)
    }

    anchor_widths <- c(df[[3]] - df[[2]], df[[6]] - df[[5]])
    narrow <- sum(anchor_widths < 10, na.rm = TRUE)
    wide   <- sum(anchor_widths > 50000, na.rm = TRUE)
    if (!quiet && narrow > 0) {
        warning(narrow, " anchor(s) are narrower than 10 bp. ",
                "Very narrow anchors may indicate data artefacts.", call. = FALSE)
    }
    if (!quiet && wide > 0) {
        warning(wide, " anchor(s) are wider than 50 kb. ",
                "Unusually wide anchors may represent broad domains (e.g. super-enhancers).",
                call. = FALSE)
    }

    df <- as.data.frame(df)

    # Normalise anchor order: first anchor lexicographically <= second anchor.
    # Only columns 1-6 are swapped; columns 7+ (name, score, strand, etc.)
    # stay in place because they are interaction-level metadata, not
    # per-anchor attributes. If columns 7+ represent per-anchor data
    # (e.g. strand1 / strand2), reorder them before calling this function.
    swap <- (df[, 1] > df[, 4]) | (df[, 1] == df[, 4] & df[, 2] > df[, 5])
    if (any(swap)) {
        df[swap, seq_len(6)] <- df[swap, c(4, 5, 6, 1, 2, 3)]
    }

    df[[2]] <- df[[2]] + 1L
    df[[5]] <- df[[5]] + 1L

    df
}

#' Internal: Compute Numeric Fraction of a Column
#'
#' Returns the fraction of non-empty, non-NA values that can be parsed as
#' numeric.  Used for auto-detecting score columns in BEDPE/BED input.
#'
#' @param x A vector (typically a column from a data frame).
#' @return A numeric value between 0 and 1.
#' @keywords internal
#' @noRd
.numeric_ratio <- function(x) {
    x <- as.character(x)
    x <- x[!is.na(x) & trimws(x) != ""]
    if (length(x) == 0) {
        return(0)
    }
    mean(!is.na(suppressWarnings(as.numeric(x))))
}

#' Read BEDPE File into a GInteractions Object
#'
#' Reads a standard BEDPE file and converts it into a Bioconductor
#' \code{\link[InteractionSet]{GInteractions}} object.
#'
#' @details
#' \strong{Anchor Normalization:}
#' Anchor order is automatically normalized so that the first anchor is
#' lexicographically less than or equal to the second (e.g., chr1 < chr2),
#' ensuring compatibility with \code{GInteractions(mode = "strict")}.
#'
#' \strong{Score Detection:}
#' The function attempts to automatically detect a numeric score column.
#' It checks the 8th column first (standard for many tools); if not numeric,
#' it falls back to the 7th column. Non-numeric values are treated as 0.
#' When \code{score_col} is specified explicitly, the function requires that
#' column to contain predominantly numeric values (>=50\%) and will stop
#' otherwise.
#'
#' @param bedpe_file Character. Path to a BEDPE file. Must contain at least six columns:
#'   \code{chr1, start1, end1, chr2, start2, end2}.
#' @param score_col Integer or \code{NULL}. Column index to use as interaction
#'   score (e.g. PET count, -log10(p-value)). If \code{NULL} (default), column 8
#'   is tried first, then column 7; if neither is numeric, scores default to 0.
#'   Set explicitly when the score column position differs from the standard or
#'   when auto-detection picks the wrong column. Note: downstream filtering
#'   parameters such as \code{min_score} in
#'   \code{\link{consolidate_chromatin_loops}} assume \emph{higher scores = better
#'   interactions}. If your score column contains p-values or other metrics where
#'   \emph{smaller is better}, convert it first (e.g. \code{-log10(p)}) or filter
#'   manually after consolidation.
#' @param quiet Logical. If \code{TRUE}, suppress data-quality warnings
#'   (e.g., unusually narrow/wide anchors, p-value-like scores).
#'   Errors are never suppressed. Default: \code{FALSE}.
#' @return A \code{\link[InteractionSet]{GInteractions}} object with a \code{score} metadata column
#'   (defaulting to 0 if not provided).
#' @importFrom data.table fread
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom InteractionSet GInteractions
#' @importFrom S4Vectors mcols mcols<-
#' @export
#' @examples
#' # 1. Locate the example BEDPE file included in the package
#' # system.file finds the absolute path to 'inst/extdata/example_loops_1.bedpe'
#' bedpe_path <- system.file("extdata", "example_loops_1.bedpe", package = "looplook")
#'
#' # 2. Run the function (ensure file was found)
#' gi <- bedpe_to_gi(bedpe_path)
#'
#' # 3. Inspect the result
#' print(gi)
#'
#' # Check the imported score column
#' S4Vectors::mcols(gi)$score
bedpe_to_gi <- function(bedpe_file, score_col = NULL, quiet = FALSE) {
    if (is.null(bedpe_file) || !file.exists(bedpe_file)) {
        stop("BEDPE file does not exist or path is invalid: ", bedpe_file)
    }
    df <- data.table::fread(bedpe_file, header = FALSE)

    df <- .validate_bedpe_df(df, quiet = quiet)

    # Score detection
    .is_pvalue_like <- function(x) {
        nums <- suppressWarnings(as.numeric(as.character(x)))
        nums <- nums[!is.na(nums) & trimws(as.character(x)) != ""]
        if (length(nums) == 0) return(FALSE)
        mean(nums >= 0 & nums <= 1, na.rm = TRUE) >= 0.8
    }

    final_scores <- rep(0, nrow(df))
    found <- FALSE
    score_is_pvalue <- FALSE

    if (!is.null(score_col)) {
        if (!is.numeric(score_col) || length(score_col) != 1L ||
            is.na(score_col) || score_col < 1 ||
            score_col != as.integer(score_col)) {
            stop("score_col must be a single positive integer column index.", call. = FALSE)
        }
        if (score_col > ncol(df)) {
            stop("score_col = ", score_col, " exceeds file column count (", ncol(df), ")")
        }
        if (.numeric_ratio(df[[score_col]]) < 0.5) {
            stop("score_col = ", score_col,
                 " does not contain predominantly numeric values.",
                 call. = FALSE)
        }
        final_scores <- suppressWarnings(as.numeric(df[[score_col]]))
        score_is_pvalue <- .is_pvalue_like(df[[score_col]])
        found <- TRUE
    } else {
        if (ncol(df) >= 8 && .numeric_ratio(df[[8]]) >= 0.5) {
            final_scores <- suppressWarnings(as.numeric(df[[8]]))
            score_is_pvalue <- .is_pvalue_like(df[[8]])
            found <- TRUE
        }
        if (!found && ncol(df) >= 7 && .numeric_ratio(df[[7]]) >= 0.5) {
            final_scores <- suppressWarnings(as.numeric(df[[7]]))
            score_is_pvalue <- .is_pvalue_like(df[[7]])
            found <- TRUE
        }
        if (!quiet && !found && ncol(df) >= 7) {
            warning("Scores defaulted to 0: no column beyond 6 had >50% numeric values")
        }
        if (!quiet && found && score_is_pvalue) {
            warning("The auto-detected score column contains values predominantly ",
                    "in [0, 1], which resemble p-values. Downstream ",
                    "min_score / min_raw_score filtering assumes ",
                    "higher scores = better interactions. ",
                    "If this column contains p-values, convert it first ",
                    "(e.g., -log10(p)) or set score_col explicitly to another column.",
                    call. = FALSE)
        }
    }

    final_scores[is.na(final_scores)] <- 0

    gr1 <- .with_known_upstream_noise_suppressed(
        GenomicRanges::GRanges(df[, 1], IRanges::IRanges(df[, 2], df[, 3]))
    )
    gr2 <- .with_known_upstream_noise_suppressed(
        GenomicRanges::GRanges(df[, 4], IRanges::IRanges(df[, 5], df[, 6]))
    )
    gi <- .with_known_upstream_noise_suppressed(
        InteractionSet::GInteractions(gr1, gr2, mode = "strict")
    )
    S4Vectors::mcols(gi)$score <- final_scores

    return(gi)
}


#' Spatial Clustering of GInteractions
#'
#' Merges spatially proximal chromatin loops into consensus interactions using graph-based clustering.
#' Loops are considered overlapping if both anchors are within \code{gap} bp of each other.
#' Each resulting cluster is represented by the \strong{union genomic range (min start to max end)} spanning all its members.
#'
#' @param gi A \code{\link[InteractionSet]{GInteractions}} object.
#' @param gap Numeric. Maximum distance (in base pairs) allowed between anchors to consider two loops overlapping. Default: 1000.
#' @return A list with two elements:
#' \describe{
#'   \item{\code{gi}}{Reduced \code{\link[InteractionSet]{GInteractions}} object, one per cluster.}
#'   \item{\code{membership}}{Integer vector indicating cluster assignment for each input loop.}
#' }
#' Metadata columns include \code{cluster_id}, \code{n_members}, \code{n_reps}, and averaged \code{score}.
#' @importFrom igraph make_empty_graph add_edges components
#' @importFrom S4Vectors queryHits subjectHits mcols
#' @importFrom GenomicRanges GRangesList
#' @export
#' @examples
#' # 1. Load example data (loops that are close to each other)
#' bedpe_path <- system.file("extdata", "example_loops_1.bedpe", package = "looplook")
#' # Convert BEDPE to GInteractions object
#' gi_raw <- bedpe_to_gi(bedpe_path)
#'
#' # 2. Run clustering
#' res <- reduce_ginteractions(gi_raw, gap = 1000)
#'
#' # 3. Inspect results
#' print(res$gi)
#' head(res$membership)
#' table(res$membership)
reduce_ginteractions <- function(gi, gap = 1000) {
    if (length(gi) == 0) {
        return(list(gi = gi, membership = integer(0)))
    }

    dt <- gi_to_dt(gi)
    gap_diag <- .diagnose_gap(dt, gap, message)
    dt <- cluster_loops_dt(dt, gap)
    reduced_dt <- reduce_clusters_dt(dt)
    gi_red <- dt_to_gi(reduced_dt)
    cluster_diag <- .diagnose_clusters(gi_red, reduced_dt, gap, message,
                       med_width = if (!is.null(gap_diag)) gap_diag$anchor_width_median else NULL,
                       n_input_loops = nrow(dt))

    attr(gi_red, "looplook_gap_diagnosis") <- gap_diag
    attr(gi_red, "looplook_cluster_diagnosis") <- cluster_diag

    list(gi = gi_red, membership = dt$cluster)
}

#' Read a Simple BED File into a GRanges Object
#'
#' Reads the first three columns of a BED file (chrom, start, end) and returns a
#' \code{\link[GenomicRanges]{GRanges}} object. Additional columns are ignored.
#'
#' @param bed_file Character. Path to a BED file (must be tab-delimited).
#'   ENCODE narrowPeak and broadPeak formats are also accepted (the first three
#'   columns are chr, start, end in all three formats).
#' @param quiet Logical. If \code{TRUE}, suppress data-quality warnings
#'   (e.g., unusually narrow/wide intervals).
#'   Errors are never suppressed. Default: \code{FALSE}.
#' @return A \code{\link[GenomicRanges]{GRanges}} object, or \code{NULL} if \code{bed_file} is \code{NULL}.
#' @importFrom data.table fread
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @export
#' @examples
#' # 1. Locate the example BED file included in the package
#' bed_path <- system.file("extdata", "example_peaks.bed", package = "looplook")
#'
#' # 2. Run the function
#' gr <- read_simple_bed(bed_path)
#'
#' # 3. Inspect the result
#' print(gr)
#' length(gr)
read_simple_bed <- function(bed_file, quiet = FALSE) {
    if (is.null(bed_file)) {
        return(NULL)
    }
    if (!file.exists(bed_file)) {
        stop("BED file does not exist: ", bed_file)
    }

    df <- data.table::fread(bed_file, header = FALSE, select = c(1, 2, 3))
    # Auto-detect header: only examine the first row to avoid misidentifying
    # a malformed data row as a header line.
    has_header <- nrow(df) > 0 && (
        is.na(suppressWarnings(as.numeric(df[[2]][1]))) ||
            is.na(suppressWarnings(as.numeric(df[[3]][1]))))
    if (has_header) {
        df <- data.table::fread(bed_file, header = FALSE, select = c(1, 2, 3), skip = 1)
    }
    # Validate all coordinate columns are numeric after header handling
    df[[2]] <- suppressWarnings(as.numeric(df[[2]]))
    df[[3]] <- suppressWarnings(as.numeric(df[[3]]))
    if (any(is.na(df[[2]]) | is.na(df[[3]]))) {
        stop("BED file contains non-numeric start/end coordinates after header handling.",
            call. = FALSE)
    }
    if (nrow(df) > 0 && any(df[[2]] >= df[[3]], na.rm = TRUE)) {
        stop("BED file contains intervals with start >= end (zero-width or invalid).",
             call. = FALSE)
    }

    # Check for unusually narrow or wide intervals (BED 0-based: width = end - start)
    bed_widths <- df[[3]] - df[[2]]
    narrow_bed <- sum(bed_widths < 10, na.rm = TRUE)
    wide_bed   <- sum(bed_widths > 50000, na.rm = TRUE)
    if (!quiet && narrow_bed > 0) {
        warning(narrow_bed, " interval(s) are narrower than 10 bp. ",
                "Very narrow intervals may indicate data artefacts.",
                call. = FALSE)
    }
    if (!quiet && wide_bed > 0) {
        warning(wide_bed, " interval(s) are wider than 50 kb. ",
                "Unusually wide intervals may represent broad domains ",
                "(e.g. super-enhancers) rather than focal peaks.",
                call. = FALSE)
    }

    # BED is 0-based half-open; GRanges expects 1-based closed
    .with_known_upstream_noise_suppressed(
        GenomicRanges::GRanges(df[[1]], IRanges::IRanges(df[[2]] + 1, df[[3]]))
    )
}


#' Consolidate and Integrate Chromatin Loops from Replicates or Multiple Sources
#'
#' @description
#' This function consolidates chromatin loops from multiple BEDPE files. It is designed for two main purposes:
#' \enumerate{
#'   \item \strong{Replicate Consolidation}: Merging biological or technical replicates to identify high-confidence, reproducible loops (e.g., 3 replicates of H3K27ac HiChIP).
#'   \item \strong{Multi-Omics Integration}: The framework can be used to identify multi-source consensus by integrating datasets from various experimental designs, such as HiChIP assays targeting different factors (e.g., integrating \strong{H3K27ac} and \strong{H3K4me3}, or overlapping Hi-C with ChIA-PET).
#' }
#'
#' The function supports three modes:
#' \itemize{
#'   \item \code{"consensus"}: Implements graph-based connected component analysis to cluster spatially proximal anchors across samples. Only retains clusters detected in >= min_consensus biological replicates.
#'   \item \code{"intersect"}: Reference-based filtering. Retains loops in File 1 whose anchors overlap with loops in every other file within the specified \code{gap} tolerance. Coordinates and scores are inherited from File 1 without merging.
#'   \item \code{"union"}: Retains all chromatin interactions across the entire cohort, ideal for exploratory pan-tissue analyses.
#' }
#'
#' \strong{Connected-component chaining}: Graph-based clustering may
#' transitively chain loci (A-B and B-C merge, pulling A and C into the same
#' cluster even if they are far apart). By default, a warning is emitted when
#' any cluster span exceeds the chaining threshold
#' (\code{max(3*gap, 5*med_width*gap/1000)}). Use \code{chaining_policy} to
#' control this behavior (\code{"warn"}, \code{"none"}, \code{"drop"}, or
#' \code{"error"}). \strong{Note:} the default \code{"warn"} policy does
#' \emph{not} remove affected clusters — their inflated coordinates (from
#' \code{min(start)} to \code{max(end)} across all chained members) flow
#' into downstream analyses unchanged. When chaining is a concern for
#' downstream overlapping operations, use \code{chaining_policy = "drop"}
#' to exclude wide clusters, or reduce \code{gap} to prevent chaining
#' from forming.
#'
#' It also supports a \strong{two-stage filtering strategy} to improve signal-to-noise ratio:
#' \itemize{
#'   \item \strong{Pre-filtering} (\code{min_raw_score}): Removes low-confidence noise (e.g., singleton reads) from raw files \emph{before} merging.
#'   \item \strong{Post-filtering} (\code{min_score}): Filters the final consensus loops based on their aggregated score.
#'   \item \strong{Replicate-Balanced Aggregation}: In \code{"consensus"} and \code{"union"} modes, each cluster score is computed as the mean of per-source mean scores, so one replicate with many fragmented calls cannot dominate the final representative score.
#' }
#'
#' @param files Character vector. Paths to BEDPE files (at least two).
#' @param gap Numeric. Maximum distance (in base pairs) between loop anchors to consider them as overlapping. Default: \code{1000} (1 kb). Typical range: 500-5000 bp depending on data resolution. Use \code{gap = 0} for strict overlap only.
#' @param mode Character. Choose one of the following: "consensus", "intersect", "union". Merge strategy:
#'   \itemize{
#'     \item \code{"consensus"}: Graph-based clustering to find a consensus set supported by a majority of samples (default). Formerly "reproducible".
#'     \item \code{"intersect"}: Strict reference-based filtering (keeps loops in File 1 supported by ALL other files).
#'     \item \code{"union"}: Merges all detected loops into a comprehensive map.
#'   }
#' @param min_consensus Integer. Minimum number of replicates a loop must appear in
#'   (only effective when \code{mode = "consensus"}).
#'   If \code{NULL} (default), the threshold is automatically calculated:
#'   \itemize{
#'     \item For 2-3 replicates: Requires all (100\%).
#'     \item For \eqn{N \ge 4}: Requires \eqn{\ge 75\%} of replicates
#'       (\code{ceiling(0.75 * n_reps)}; e.g., 3 for N=4, 4 for N=5, 6 for N=8).
#'   }
#' @param min_raw_score Numeric. \strong{Pre-filtering threshold}. Loops with a raw score (e.g., read count) below this value in individual files will be discarded \strong{before} any merging or intersection.
#'   \itemize{
#'     \item Recommended value: \code{2} (to remove singleton noise loops with count=1).
#'     \item Default: \code{NULL} (no pre-filtering).
#'   }
#' @param min_score Numeric. \strong{Post-filtering threshold}. Minimum score to keep a consolidated loop \strong{after} merging.
#'   \itemize{
#'     \item For \code{"consensus"} and \code{"union"} modes, this filters loops based on a replicate-balanced representative score (mean of per-source cluster means).
#'     \item For \code{"intersect"} mode, this filters the retained File 1 loops by their original score.
#'     \item Default: \code{NULL} (no post-filtering).
#'   }
#' @param score_col Integer or \code{NULL}. Column index to use as interaction
#'   score when reading BEDPE files. Passed to \code{\link{bedpe_to_gi}}.
#'   If \code{NULL} (default), auto-detection is used (see \code{\link{bedpe_to_gi}}).
#' @param chaining_policy Character. Controls behaviour when connected-component
#'   chaining produces clusters with span exceeding the chaining threshold
#'   (\code{max(3*gap, 5*med_width*gap/1000)}).
#'   \code{"warn"} (default): emit a warning and retain all clusters
#'   \emph{unchanged}, including their potentially inflated coordinates.
#'   \code{"none"}: silently accept all clusters.
#'   \code{"drop"}: remove clusters exceeding the threshold (safe for
#'   downstream overlapping analyses that may be affected by inflated spans).
#'   \code{"error"}: stop with an error.
#' @param blacklist_species Character. Species/build for built-in ENCODE
#'   blacklist (\code{"hg38"}, \code{"hg19"}, \code{"mm10"}, \code{"mm9"}),
#'   or a path to a custom BED file. When a species name is recognised, the
#'   bundled blacklist is used; otherwise the value is treated as a file path.
#' @param region_of_interest Character. Path to BED file defining regions of interest (ROI). Only loops overlapping these regions will be kept.
#' @param roi_mode Character. How loops must overlap \code{region_of_interest}.
#'   \code{"any"} (default): keep loops where \emph{either} anchor overlaps the ROI
#'   (suitable for promoter-centric or enhancer-gene queries).
#'   \code{"both"}: keep loops where \emph{both} anchors overlap the ROI
#'   (suitable for TAD confinement or domain-internal interaction queries).
#' @param out_file Character. The file name (including the file path) for saving results in the extended BEDPE format.
#' @param write_output Logical. If \code{TRUE} (default), write the consolidated BEDPE file when \code{out_file} is provided. If \code{FALSE}, return the \code{GInteractions} object without creating directories or files.
#' @param quiet Logical. If \code{TRUE}, suppress progress messages while preserving warnings. Default: \code{FALSE}.
#' @return A filtered \code{\link[InteractionSet]{GInteractions}} object with metadata columns:
#'   \describe{
#'     \item{\code{score}}{Replicate-balanced consensus score.}
#'     \item{\code{n_members}}{Number of raw loops merged into this entry
#'       (1 for intersect mode where no coordinate merging occurs).}
#'     \item{\code{n_reps}}{Number of input files that support this entry.}
#'     \item{\code{cluster_id}}{Connected-component cluster identifier.}
#'   }
#'   The returned object carries a \code{looplook_metadata} attribute (access via
#'   \code{attr(x, "looplook_metadata")}) with package version, call parameters,
#'   diagnostics, and database versions.
#'   When \code{write_output = TRUE} and \code{out_file} is provided, an extended
#'   BEDPE file is written with the additional columns \code{n_members} and
#'   \code{n_reps} appended after the standard BEDPE fields.
#' @importFrom data.table fread
#' @importFrom GenomicRanges GRanges seqnames start end
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors mcols mcols<- queryHits subjectHits
#' @importFrom igraph make_empty_graph add_edges components
#' @export
#' @examples
#' # 1. Get paths to example BEDPE files included in the package
#' f1 <- system.file("extdata", "example_loops_1.bedpe", package = "looplook")
#' f2 <- system.file("extdata", "example_loops_2.bedpe", package = "looplook")
#'
#' # 2. Run consolidation (ensure files exist)
#' # Example A: Intersect Mode
#' # Only keeps loops present in f1 that are also supported by f2
#' res_intersect <- consolidate_chromatin_loops(
#'     files = c(f1, f2),
#'     mode = "intersect",
#'     gap = 1000,
#'     out_file = tempfile(fileext = ".bedpe")
#' )
#'
#' # Example B: Consensus Mode (formerly Reproducible)
#' # Finds consensus loops supported by both replicates (default for N=2)
#' res_consensus <- consolidate_chromatin_loops(
#'     files = c(f1, f2),
#'     mode = "consensus",
#'     gap = 1000,
#'     out_file = tempfile(fileext = ".bedpe")
#' )
#'
#' # Example C: Union Mode
#' # Merges all loops into a single map
#' res_union <- consolidate_chromatin_loops(
#'     files = c(f1, f2),
#'     mode = "union",
#'     gap = 1000,
#'     out_file = tempfile(fileext = ".bedpe")
#' )
#'
#' # Example D: Dual Filtering Strategy (Recommended for HiChIP)
#' # 1. Pre-filter: Discard singletons (score < 2) to remove noise.
#' # 2. Merge: Find loops present in both replicates.
#' # 3. Post-filter: Keep only strong consensus loops (score > 5).
#' res_clean <- consolidate_chromatin_loops(
#'     files = c(f1, f2),
#'     mode = "consensus",
#'     min_raw_score = 2, # Pre-filter (remove noise)
#'     min_score = 5, # Post-filter (keep strong loops)
#'     gap = 1000,
#'     out_file = tempfile(fileext = ".bedpe")
#' )
#'
#' # Inspect results
#' length(res_intersect)
#' length(res_clean)
consolidate_chromatin_loops <- function(
  files = NULL,
  gap = 1000,
  mode = c("consensus", "intersect", "union"),
  min_consensus = NULL,
  score_col = NULL,
  min_raw_score = NULL,
  min_score = NULL,
  chaining_policy = c("warn", "none", "drop", "error"),
  blacklist_species = NULL,
  region_of_interest = NULL,
  roi_mode = c("any", "both"),
  out_file = NULL,
  write_output = TRUE,
  quiet = FALSE
) {
    log_message <- function(...) {
        if (!quiet) message(...)
    }

    if (is.null(files) || length(files) < 2) {
        stop("`files` must contain at least two BEDPE file paths.", call. = FALSE)
    }
    missing_files <- files[!file.exists(files)]
    if (length(missing_files) > 0) {
        stop("BEDPE files not found: ", paste(missing_files, collapse = ", "), call. = FALSE)
    }
    mode <- match.arg(mode)
    roi_mode <- match.arg(roi_mode)
    chaining_policy <- match.arg(chaining_policy)

    # Parameter validation
    if (!is.numeric(gap) || length(gap) != 1L || is.na(gap) || gap < 0)
        stop("`gap` must be a non-negative number", call. = FALSE)
    if (!is.null(min_consensus) &&
        (!is.numeric(min_consensus) || length(min_consensus) != 1L ||
         is.na(min_consensus) || min_consensus < 1))
        stop("`min_consensus` must be a positive integer or NULL", call. = FALSE)
    if (!is.null(min_raw_score) &&
        (!is.numeric(min_raw_score) || length(min_raw_score) != 1L ||
         is.na(min_raw_score)))
        stop("`min_raw_score` must be a single number or NULL", call. = FALSE)
    if (!is.null(min_score) &&
        (!is.numeric(min_score) || length(min_score) != 1L ||
         is.na(min_score)))
        stop("`min_score` must be a single number or NULL", call. = FALSE)
    if (!is.null(blacklist_species)) {
        if (!is.character(blacklist_species) || length(blacklist_species) != 1L ||
            is.na(blacklist_species) || !nzchar(blacklist_species))
            stop("`blacklist_species` must be a non-empty string or NULL", call. = FALSE)
    }
    if (!is.null(region_of_interest)) {
        if (!is.character(region_of_interest) || length(region_of_interest) != 1L ||
            is.na(region_of_interest) || !nzchar(region_of_interest))
            stop("`region_of_interest` must be a non-empty file path or NULL", call. = FALSE)
        if (!file.exists(region_of_interest))
            stop("`region_of_interest` file not found: ", region_of_interest, call. = FALSE)
    }

    n_reps <- length(files)

    # --- Read & pre-filter files ---
    gi_list <- .consolidate_read_files(
        files, score_col, min_raw_score, quiet, log_message
    )

    total_loops <- sum(vapply(gi_list, length, integer(1)))
    if (total_loops == 0) {
        if (!quiet) {
            message("All loops filtered out by min_raw_score = ", min_raw_score,
                    ". Returning empty result.")
        }
        empty_gi <- InteractionSet::GInteractions(
            GenomicRanges::GRanges(), GenomicRanges::GRanges()
        )
        S4Vectors::mcols(empty_gi)$score <- numeric(0)
        S4Vectors::mcols(empty_gi)$n_members <- integer(0)
        S4Vectors::mcols(empty_gi)$n_reps <- integer(0)
        return(empty_gi)
    }

    # --- Merge by mode ---
    if (mode == "intersect") {
        result_gi <- .consolidate_intersect(gi_list, gap, n_reps, log_message)
    } else {
        result_gi <- .consolidate_cluster_mode(
            gi_list, mode, gap, min_consensus, n_reps,
            chaining_policy, log_message
        )
    }

    # --- Post-filters ---
    result_gi <- .consolidate_post_filters(
        result_gi, min_score, blacklist_species,
        region_of_interest, roi_mode, quiet, log_message
    )

    # --- Export ---
    if (write_output && !is.null(out_file)) {
        .consolidate_export(result_gi, out_file, log_message)
    }

    log_message("Finished! Final loops: ", length(result_gi))

    # Collect diagnostic data attached by cluster_mode / intersect
    gap_diag <- attr(result_gi, "looplook_gap_diagnosis")
    cluster_diag <- attr(result_gi, "looplook_cluster_diagnosis")

    attr(result_gi, "looplook_metadata") <- .build_looplook_metadata(
        fun = "consolidate_chromatin_loops",
        params = list(
            files = basename(files), n_files = n_reps,
            gap = gap, mode = mode,
            min_consensus = min_consensus,
            min_raw_score = min_raw_score, min_score = min_score,
            chaining_policy = chaining_policy,
            roi_mode = roi_mode,
            blacklist_species = blacklist_species
        ),
        score_semantics = if (is.null(min_raw_score) && is.null(min_score)) {
            "raw score (higher = better)"
        } else {
            "filtered score (higher = better)"
        },
        diagnostics = list(
            gap = gap_diag,
            cluster = cluster_diag
        )
    )
    result_gi
}

# --- Helpers extracted from consolidate_chromatin_loops ---

#' Internal: Read and pre-filter BEDPE files for consolidation
#' @keywords internal
#' @noRd
.consolidate_read_files <- function(
    files, score_col, min_raw_score, quiet, log_message
) {
    log_message(">>> Reading BEDPE files")
    gi_list <- lapply(files, function(f) {
        gi <- bedpe_to_gi(f, score_col = score_col, quiet = quiet)
        if (!is.null(min_raw_score)) {
            if ("score" %in% colnames(S4Vectors::mcols(gi))) {
                keep_idx <- S4Vectors::mcols(gi)$score >= min_raw_score
                gi <- gi[keep_idx]
            }
        }
        return(gi)
    })
    for (i in seq_along(gi_list)) {
        if (length(gi_list[[i]]) > 0) {
            S4Vectors::mcols(gi_list[[i]])$source <- i
        }
        log_message("    File ", i, ": ", length(gi_list[[i]]), " loops")
    }
    gi_list
}

#' Internal: Intersect-mode consolidation (reference-based, no coordinate merging)
#' @keywords internal
#' @noRd
.consolidate_intersect <- function(gi_list, gap, n_reps, log_message) {
    log_message(">>> Intersect mode: Reference-based filtering (No Coordinate Merging)")
    log_message("    Base: File 1. Criterion: Must overlap with ALL other files.")
    current_gi <- gi_list[[1]]
    for (i in 2:n_reps) {
        if (length(current_gi) == 0) break
        log_message("    Intersecting with File ", i, "...")
        hits <- InteractionSet::findOverlaps(
            current_gi, gi_list[[i]], maxgap = gap, use.region = "both"
        )
        keep_idx <- unique(S4Vectors::queryHits(hits))
        current_gi <- current_gi[keep_idx]
    }
    S4Vectors::mcols(current_gi)$n_reps <- n_reps
    S4Vectors::mcols(current_gi)$n_members <- 1L
    S4Vectors::mcols(current_gi)$cluster_id <- as.character(seq_along(current_gi))
    current_gi
}

#' Internal: Cluster-mode consolidation (consensus/union) with chaining check
#' @keywords internal
#' @noRd
.consolidate_cluster_mode <- function(
    gi_list, mode, gap, min_consensus, n_reps,
    chaining_policy, log_message
) {
    log_message(">>> Clustering mode (Union/Consensus): Merging coordinates via Graph")

    # Seqlevels consistency check across input files.
    # chr1 vs 1, or mixed UCSC/Ensembl styles, cause silent false negatives
    # in the string-based chr matching inside cluster_loops_dt().
    sl_styles <- unique(vapply(gi_list, function(gi) {
        tryCatch(GenomeInfoDb::seqlevelsStyle(InteractionSet::anchors(gi, "first"))[1],
                 error = function(e) NA_character_)
    }, character(1)))
    sl_styles <- sl_styles[!is.na(sl_styles)]
    if (length(sl_styles) > 1L) {
        warning("Input BEDPE files have mixed seqlevel styles: ",
                paste(sl_styles, collapse = ", "), ". ",
                "Loop clustering uses string-based chromosome matching. ",
                "Mismatched styles (e.g. 'chr1' vs '1') will produce ",
                "false negatives. Harmonise your input files first.",
                call. = FALSE)
    }

    combined_dt <- data.table::rbindlist(lapply(gi_list, gi_to_dt))

    # Pre-clustering gap diagnosis; capture anchor width for post-clustering use
    gap_diag <- .diagnose_gap(combined_dt, gap, log_message)
    med_width <- if (!is.null(gap_diag)) gap_diag$anchor_width_median else NULL
    n_input <- nrow(combined_dt)

    clustered <- cluster_loops_dt(combined_dt, gap)
    reduced_dt <- reduce_clusters_dt(clustered)

    if (mode == "consensus") {
        if (is.null(min_consensus)) {
            if (n_reps <= 3) {
                min_consensus <- n_reps  # 100% for 2-3 replicates
            } else {
                min_consensus <- ceiling(0.75 * n_reps)  # >=75% for N >= 4
            }
        }
        log_message(">>> Consensus mode: Keeping clusters in >= ", min_consensus, " replicates")
        reduced_dt <- reduced_dt[n_reps >= min_consensus]
    } else {
        log_message(">>> Union mode: Keeping all clusters")
    }

    result_gi <- dt_to_gi(reduced_dt)

    # Post-clustering diagnosis with anchor-width-aware chaining threshold
    cluster_diag <- .diagnose_clusters(result_gi, reduced_dt, gap, log_message,
                       med_width = med_width, n_input_loops = n_input)

    # Attach diagnostic metrics to result for programmatic access
    attr(result_gi, "looplook_gap_diagnosis") <- gap_diag
    attr(result_gi, "looplook_cluster_diagnosis") <- cluster_diag

    # Chaining span check (uses combined threshold: anchor-width-aware)
    if (chaining_policy != "none" && length(result_gi) > 0) {
        a1 <- InteractionSet::anchors(result_gi, "first")
        a2 <- InteractionSet::anchors(result_gi, "second")
        span1 <- GenomicRanges::end(a1) - GenomicRanges::start(a1) + 1L
        span2 <- GenomicRanges::end(a2) - GenomicRanges::start(a2) + 1L
        span_threshold <- if (!is.null(med_width) && is.finite(med_width) && med_width > 0) {
            max(3 * gap, 5 * med_width * gap / 1000)
        } else {
            3 * gap
        }
        nm <- S4Vectors::mcols(result_gi)$n_members
        # Only flag clusters with n_members > 2 as potential chaining.
        # n_members == 2 typically means the same loop detected in two
        # replicates; wide span is from naturally broad anchors (super-
        # enhancers), not from chaining through intermediate loops.
        wide_idx <- which(pmax(span1, span2) > span_threshold & nm > 2)
        n_wide <- length(wide_idx)
        if (n_wide > 0) {
            if (chaining_policy == "warn") {
                warning(
                    n_wide, " cluster(s) have max_span > chaining threshold (",
                    format(round(span_threshold), big.mark = ","), " bp). ",
                    "Connected-component clustering may have chained ",
                    "through intermediate loops. Consider reducing 'gap' ",
                    "or inspecting clusters with large 'n_members'.",
                    call. = FALSE
                )
            } else if (chaining_policy == "drop") {
                log_message(
                    ">>> Dropping ", n_wide, " cluster(s) with max_span ",
                    "> chaining threshold (", format(round(span_threshold), big.mark = ","), " bp)"
                )
                result_gi <- result_gi[-wide_idx]
            } else if (chaining_policy == "error") {
                stop(
                    n_wide, " cluster(s) have max_span > chaining threshold (",
                    format(round(span_threshold), big.mark = ","), " bp). ",
                    "Connected-component clustering may have chained ",
                    "through intermediate loops. Set chaining_policy to ",
                    "'warn', 'drop', or reduce 'gap'."
                )
            }
        }
    }
    result_gi
}

#' Internal: Apply post-merge filters (min_score, blacklist, ROI)
#' @keywords internal
#' @noRd
.consolidate_post_filters <- function(
    result_gi, min_score, blacklist_species,
    region_of_interest, roi_mode, quiet, log_message
) {
    if (!is.null(min_score)) {
        keep <- S4Vectors::mcols(result_gi)$score >= min_score
        result_gi <- result_gi[keep]
    }

    if (!is.null(blacklist_species)) {
        known_lists <- list(
            "hg38" = "hg38-blacklist.v2.bed",
            "hg19" = "hg19-blacklist.v2.bed",
            "mm10" = "mm10-blacklist.v2.bed",
            "mm9"  = "mm9-blacklist.v2.bed"
        )
        blacklist_path <- NULL
        if (blacklist_species %in% names(known_lists)) {
            blacklist_path <- system.file(
                "extdata", known_lists[[blacklist_species]], package = "looplook"
            )
        }
        if (is.null(blacklist_path) || blacklist_path == "") {
            blacklist_path <- blacklist_species
        }
        if (file.exists(blacklist_path)) {
            log_message(">>> Filtering blacklist: ", basename(blacklist_path))
            bl <- read_simple_bed(blacklist_path, quiet = quiet)
            bl <- .harmonize_seqlevels(
                bl, InteractionSet::anchors(result_gi, "first"), "blacklist"
            )
            h1 <- InteractionSet::findOverlaps(
                InteractionSet::anchors(result_gi, "first"), bl
            )
            h2 <- InteractionSet::findOverlaps(
                InteractionSet::anchors(result_gi, "second"), bl
            )
            bad <- unique(c(S4Vectors::queryHits(h1), S4Vectors::queryHits(h2)))
            if (length(bad)) result_gi <- result_gi[-bad]
        } else {
            stop("Blacklist file not found: ", blacklist_species, call. = FALSE)
        }
    }

    if (!is.null(region_of_interest)) {
        log_message(">>> Filtering by region of interest (", roi_mode, "): ",
                     basename(region_of_interest))
        if (file.exists(region_of_interest)) {
            tg <- read_simple_bed(region_of_interest, quiet = quiet)
            tg <- .harmonize_seqlevels(
                tg, InteractionSet::anchors(result_gi, "first"), "ROI"
            )
            h1 <- InteractionSet::findOverlaps(
                InteractionSet::anchors(result_gi, "first"), tg
            )
            h2 <- InteractionSet::findOverlaps(
                InteractionSet::anchors(result_gi, "second"), tg
            )
            keep <- if (roi_mode == "any") {
                base::union(
                    S4Vectors::queryHits(h1), S4Vectors::queryHits(h2)
                )
            } else {
                base::intersect(
                    S4Vectors::queryHits(h1), S4Vectors::queryHits(h2)
                )
            }
            if (length(keep) > 0) {
                result_gi <- result_gi[keep]
                log_message("    Kept ", length(result_gi),
                            " loops overlapping ROI.")
            } else {
                log_message("    No loops overlapped with the ROI. Returning empty set.")
                result_gi <- result_gi[0]
            }
        } else {
            warning("Region of interest file not found: ", region_of_interest)
        }
    }
    result_gi
}

#' Internal: Export consolidated loops to BEDPE file
#' @keywords internal
#' @noRd
.consolidate_export <- function(result_gi, out_file, log_message) {
    a1 <- InteractionSet::anchors(result_gi, "first")
    a2 <- InteractionSet::anchors(result_gi, "second")
    # BEDPE export: convert back from 1-based closed to 0-based half-open
    out_df <- data.frame(
        chr1 = as.character(GenomicRanges::seqnames(a1)),
        start1 = GenomicRanges::start(a1) - 1L,
        end1 = GenomicRanges::end(a1),
        chr2 = as.character(GenomicRanges::seqnames(a2)),
        start2 = GenomicRanges::start(a2) - 1L,
        end2 = GenomicRanges::end(a2),
        name = ".",
        score = round(S4Vectors::mcols(result_gi)$score, 2),
        n_members = if (!is.null(S4Vectors::mcols(result_gi)$n_members)) {
            S4Vectors::mcols(result_gi)$n_members
        } else {
            rep(1L, length(result_gi))
        },
        n_reps = if (!is.null(S4Vectors::mcols(result_gi)$n_reps)) {
            S4Vectors::mcols(result_gi)$n_reps
        } else {
            rep(1L, length(result_gi))
        },
        stringsAsFactors = FALSE
    )
    tryCatch({
        dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
        utils::write.table(out_df, out_file,
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE
        )
        log_message("Finished! Saved to ", out_file)
    }, error = function(e)
        warning("Failed to save consolidated BEDPE: ", conditionMessage(e),
                call. = FALSE)
    )
}


gi_to_dt <- function(gi) {
    # Extracts 1-based closed coordinates from GInteractions (GRanges).
    # All downstream cluster_loops_dt gap calculations operate on this
    # 1-based closed system -- consistent with GRanges semantics and
    # Bioconductor findOverlaps(maxgap=...) conventions.
    a1 <- InteractionSet::anchors(gi, "first")
    a2 <- InteractionSet::anchors(gi, "second")

    data.table::data.table(
        chr1 = as.character(GenomicRanges::seqnames(a1)),
        start1 = GenomicRanges::start(a1),
        end1 = GenomicRanges::end(a1),
        chr2 = as.character(GenomicRanges::seqnames(a2)),
        start2 = GenomicRanges::start(a2),
        end2 = GenomicRanges::end(a2),
        score = if (!is.null(S4Vectors::mcols(gi)$score)) {
            S4Vectors::mcols(gi)$score
        } else {
            0
        },
        source = if (!is.null(S4Vectors::mcols(gi)$source)) {
            S4Vectors::mcols(gi)$source
        } else {
            1L
        }
    )
}

reduce_clusters_dt <- function(dt) {
    cluster_coords <- dt[, list(
        chr1 = chr1[1],
        start1 = min(start1),
        end1 = max(end1),
        chr2 = chr2[1],
        start2 = min(start2),
        end2 = max(end2),
        n_members = .N
    ), by = cluster]

    cluster_scores <- dt[, .(
        score = .mean_or_na(score)
    ), by = .(cluster, source)][!is.na(score), .(
        score = .mean_or_na(score),
        n_reps = .N
    ), by = cluster]

    reduced_dt <- merge(
        cluster_coords,
        cluster_scores,
        by = "cluster",
        all.x = TRUE,
        sort = FALSE
    )
    data.table::setcolorder(reduced_dt, c(
        "cluster", "chr1", "start1", "end1", "chr2", "start2", "end2",
        "score", "n_members", "n_reps"
    ))
    reduced_dt
}

.mean_or_na <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0) {
        return(NA_real_)
    }
    mean(x)
}

dt_to_gi <- function(dt) {
    gr1 <- .with_known_upstream_noise_suppressed(
        GenomicRanges::GRanges(
            dt$chr1, IRanges::IRanges(dt$start1, dt$end1)
        )
    )
    gr2 <- .with_known_upstream_noise_suppressed(
        GenomicRanges::GRanges(
            dt$chr2, IRanges::IRanges(dt$start2, dt$end2)
        )
    )

    gi <- .with_known_upstream_noise_suppressed(
        InteractionSet::GInteractions(gr1, gr2, mode = "strict")
    )
    S4Vectors::mcols(gi)$cluster_id <- dt$cluster
    S4Vectors::mcols(gi)$n_members <- dt$n_members
    S4Vectors::mcols(gi)$score <- dt$score
    S4Vectors::mcols(gi)$n_reps <- dt$n_reps
    gi
}

#' Internal: Diagnose Gap Parameter for Spatial Clustering
#'
#' Examines anchor width distribution and inter-anchor distances to
#' assess whether the current \code{gap} setting is appropriate.
#' Emits diagnostic messages with data-driven recommendations.
#' Called automatically by \code{\link{consolidate_chromatin_loops}}
#' before clustering when \code{quiet = FALSE}.
#'
#' @param dt A data.table with columns chr1, start1, end1, chr2, start2, end2.
#' @param gap Numeric. The gap parameter to diagnose.
#' @param log_message Function. Message output function.
#' @return Invisibly returns a list of diagnostic metrics, or \code{NULL}
#'   if data is too small for diagnosis.
#' @keywords internal
#' @noRd
.diagnose_gap <- function(dt, gap, log_message) {
    if (nrow(dt) < 10) {
        log_message("  [i]  Too few loops (< 10) for gap diagnosis; skipping.")
        return(invisible(NULL))
    }

    # Anchor width in bp (1-based closed: end - start + 1)
    w1 <- dt$end1 - dt$start1 + 1L
    w2 <- dt$end2 - dt$start2 + 1L
    all_widths <- c(w1, w2)

    med_w <- stats::median(all_widths, na.rm = TRUE)
    q25_w <- stats::quantile(all_widths, 0.25, na.rm = TRUE, names = FALSE)
    q75_w <- stats::quantile(all_widths, 0.75, na.rm = TRUE, names = FALSE)

    # True nearest-neighbour distance via distanceToNearest (per-chromosome).
    .true_nearest_dist <- function(chr_vec, start_vec, end_vec) {
        gr <- GenomicRanges::GRanges(chr_vec, IRanges::IRanges(start_vec, end_vec))
        hits <- GenomicRanges::distanceToNearest(gr)
        if (length(hits) == 0) return(numeric(0))
        S4Vectors::mcols(hits)$distance
    }
    dist1 <- .true_nearest_dist(dt$chr1, dt$start1, dt$end1)
    dist2 <- .true_nearest_dist(dt$chr2, dt$start2, dt$end2)
    all_dists <- c(dist1, dist2)

    med_dist <- if (length(all_dists) > 0) stats::median(all_dists, na.rm = TRUE) else NA_real_
    near_fraction <- if (length(all_dists) > 0) mean(all_dists <= gap, na.rm = TRUE) else NA_real_

    # Effective width ratio: (anchor_width + 2*gap) / anchor_width
    # Measures how much gap inflates the anchor"s reach
    effective_ratio <- (med_w + 2 * gap) / med_w

    # --- Diagnostic messages ---
    log_message("--- Gap Diagnosis ---")
    log_message(sprintf(
        "  Anchor width: median = %s bp  (IQR: %s - %s bp)",
        format(med_w, big.mark = ","),
        format(q25_w, big.mark = ","),
        format(q75_w, big.mark = ",")
    ))
    if (!is.na(med_dist)) {
        log_message(sprintf(
            "  Adjacent-anchor gap: median = %s bp  (%.0f%% within current gap = %s bp)",
            format(med_dist, big.mark = ","),
            near_fraction * 100, format(gap, big.mark = ",")
        ))
    }
    log_message(sprintf(
        "  Gap / median_anchor_width ratio: %.1fx  (effective width expansion: %.1fx)",
        gap / med_w, effective_ratio
    ))

    # Risk assessment
    # Use data-type-aware thresholds for better biological relevance
    if (effective_ratio > 5 && gap > 500 && med_w < 500) {
        log_message(sprintf(
            "  [!!] RISK HIGH: gap=%s bp inflates narrow anchors (median %s bp) by %.0fx.",
            format(gap, big.mark = ","), format(med_w, big.mark = ","), effective_ratio
        ))
        log_message("  Clustering is dominated by gap rather than anchor positions.")
        log_message(sprintf(
            "  Consider reducing gap to <= %s bp (2x median width) for these narrow peaks.",
            format(round(med_w * 2), big.mark = ",")
        ))
        log_message(sprintf(
            "  Suggested: gap = %s  (or gap = %s for conservative merging)",
            format(round(med_w), big.mark = ","),
            format(round(med_w * 0.5), big.mark = ",")
        ))
    } else if (effective_ratio > 3 && med_w < 1000) {
        log_message(sprintf(
            "  [!]  RISK MODERATE: gap=%s bp is %sx wider than median anchor (%s bp).",
            format(gap, big.mark = ","), format(round(gap / med_w, 1), big.mark = ","), format(med_w, big.mark = ",")
        ))
        log_message("  Independent loops in dense regions may be merged.")
        log_message(sprintf(
            "  If your data is TF HiChIP or ChIA-PET with narrow loop anchors, ",
            "consider reducing gap to %s-%s bp.",
            format(round(med_w * 0.5), big.mark = ","),
            format(round(med_w * 2), big.mark = ",")
        ))
    } else if (effective_ratio > 2.5 && gap > med_w * 0.5) {
        # Additional check for broad enhancer data where gap exceeds half the anchor width
        log_message(sprintf(
            "  [!]  RISK MODERATE: gap=%s bp is %.1fx the median anchor width (%s bp).",
            format(gap, big.mark = ","), gap / med_w, format(med_w, big.mark = ",")
        ))
        log_message("  For broad regulatory domains, consider whether this gap merges independent elements.")
        log_message(sprintf(
            "  If your data contains dense regulatory regions, consider reducing gap to %s-%s bp.",
            format(round(med_w * 0.3), big.mark = ","),
            format(round(med_w * 0.5), big.mark = ",")
        ))
    } else {
        log_message("  Gap appears appropriate for the anchor width distribution.")
    }

    # Adjacent-anchor distance insight
    if (!is.na(near_fraction) && near_fraction < 0.01 && med_w < 1000 && gap > 1000) {
        log_message(sprintf(
            "  [i]  Only %.1f%% of adjacent anchors fall within gap=%s bp. ",
            near_fraction * 100, format(gap, big.mark = ",")
        ))
        log_message(sprintf(
            "  Gap may be unnecessarily small; most anchors are far apart. ",
            "Increasing gap would not cause over-merging for this dataset."
        ))
    }

    # Data type inference (based on loop anchor width, not peak width)
    if (med_w < 300) {
        log_message("  [i]  TF-binding resolution (CTCF, Pol II HiChIP / ChIA-PET).")
        dt_label <- "TF_binding"
    } else if (med_w < 1000) {
        log_message("  [i]  Narrow enhancer anchors -- typical of ATAC-seq HiChIP or focal histone marks.")
        dt_label <- "narrow_enhancer"
    } else if (med_w < 3000) {
        log_message("  [i]  Broad enhancer anchors -- typical of H3K27ac / H3K4me3 HiChIP or PLAC-seq.")
        dt_label <- "broad_enhancer"
    } else if (med_w < 10000) {
        log_message("  [i]  Super-enhancer-like anchors -- very broad regulatory domains.")
        dt_label <- "super_enhancer"
    } else {
        log_message("  [i]  Hi-C bin / Capture Hi-C bait-level resolution.")
        dt_label <- "HiC_bin"
    }
    log_message("--- End Gap Diagnosis ---")

    invisible(list(
        anchor_width_median = med_w,
        anchor_width_iqr = c(q25_w, q75_w),
        adjacent_gap_median = med_dist,
        adjacent_gap_near_fraction = near_fraction,
        effective_width_ratio = effective_ratio,
        data_type = dt_label
    ))
}

#' Internal: Diagnose Clustering Results
#'
#' Reports key cluster-level statistics after graph-based clustering:
#' cluster count, membership distribution, anchor span distribution,
#' consensus-filter impact, and chaining risk assessment.
#'
#' @param result_gi A \code{GInteractions} object (one per cluster).
#' @param reduced_dt A data.table of reduced cluster data (with
#'   \code{n_members} and \code{n_reps} columns).
#' @param gap Numeric. The gap parameter used for clustering.
#' @param log_message Function. Message output function.
#' @param med_width Numeric or NULL. Median anchor width from pre-clustering
#'   diagnosis. When provided, the chaining threshold is
#'   \code{max(3*gap, 5*med_width*gap/1000)} instead of \code{3*gap} alone, preventing
#'   false alarms when gap is small relative to anchor width.
#' @param n_input_loops Integer or NULL. Total number of input loops before
#'   consensus filtering. When provided, the consensus retention rate is
#'   reported.
#' @return Invisibly returns a list of cluster statistics.
#' @keywords internal
#' @noRd
.diagnose_clusters <- function(result_gi, reduced_dt, gap, log_message,
                                med_width = NULL, n_input_loops = NULL) {
    n_clusters <- length(result_gi)
    if (n_clusters == 0) {
        log_message("--- Post-Clustering Diagnosis ---")
        log_message("  No clusters survived filtering.")
        log_message("--- End Post-Clustering Diagnosis ---")
        return(invisible(NULL))
    }

    # --- Cluster size (n_members) distribution ---
    nm <- reduced_dt$n_members
    total_survived <- sum(nm)

    log_message("--- Post-Clustering Diagnosis ---")
    log_message(sprintf(
        "  Clusters formed: %s  (from %s loops surviving consensus)",
        format(n_clusters, big.mark = ","),
        format(total_survived, big.mark = ",")
    ))
    log_message(sprintf(
        "  Members per cluster: median = %.0f, IQR = %.0f-%.0f, max = %.0f",
        stats::median(nm), stats::quantile(nm, 0.25, names = FALSE),
        stats::quantile(nm, 0.75, names = FALSE), max(nm)
    ))

    # --- Consensus retention ---
    if (!is.null(n_input_loops) && n_input_loops > 0) {
        retention <- total_survived / n_input_loops * 100
        log_message(sprintf(
            "  Consensus retention: %s / %s input loops (%.1f%%)",
            format(total_survived, big.mark = ","),
            format(n_input_loops, big.mark = ","),
            retention
        ))
        if (retention < 20) {
            log_message("  [i]  Low retention -- many loops failed consensus. Gap may be too small for reproducible calls across replicates.")
        }
    }

    # --- Anchor span distribution ---
    a1 <- InteractionSet::anchors(result_gi, "first")
    a2 <- InteractionSet::anchors(result_gi, "second")
    span1 <- GenomicRanges::end(a1) - GenomicRanges::start(a1) + 1L
    span2 <- GenomicRanges::end(a2) - GenomicRanges::start(a2) + 1L
    max_span <- pmax(span1, span2)

    # --- Chaining threshold: combine gap-based and anchor-width-based ---
    chain_threshold <- if (!is.null(med_width) && is.finite(med_width) && med_width > 0) {
        max(3 * gap, 5 * med_width * gap / 1000)
    } else {
        3 * gap
    }
    threshold_src <- if (!is.null(med_width) && is.finite(med_width) && 5 * med_width * gap / 1000 > 3 * gap) {
        sprintf("5xmed_width(%s) x gap/1000", format(med_width, big.mark = ","))
    } else {
        sprintf("3xgap(%s)", format(gap, big.mark = ","))
    }

    log_message(sprintf(
        "  Cluster span: median = %s  |  max = %s  |  threshold = %s = %s bp",
        format(round(stats::median(max_span)), big.mark = ","),
        format(max(max_span), big.mark = ","),
        threshold_src,
        format(round(chain_threshold), big.mark = ",")
    ))

    # --- Top spans: show actual values of largest clusters ---
    # n_members = 2 usually means the same loop detected in two replicates
    # -- not pathological chaining.  Only flag clusters with n_members > 2.
    n_above <- sum(max_span > chain_threshold & nm > 2)
    n_show <- min(max(n_above, 3L), 5L)
    top_idx <- head(order(max_span, decreasing = TRUE), n_show)
    log_message("  Largest cluster spans:")
    for (i in top_idx) {
        flag <- if (max_span[i] > chain_threshold && nm[i] > 2) " [!]" else ""
        log_message(sprintf(
            "    #%s: max_span = %s bp, n_members = %.0f, n_reps = %.0f%s",
            i, format(max_span[i], big.mark = ","),
            nm[i], reduced_dt$n_reps[i], flag
        ))
    }

    # --- Chaining risk summary ---
    pct_chain <- n_above / n_clusters * 100
    if (n_above == 0) {
        log_message(sprintf("  Chaining: 0/%s above threshold -- PASS.", n_clusters))
    } else if (pct_chain > 20) {
        log_message(sprintf(
            "  [!!] Chaining: %s/%s (%.0f%%) above threshold -- EXCESSIVE. Consider reducing gap.",
            n_above, n_clusters, pct_chain
        ))
    } else if (pct_chain > 5) {
        log_message(sprintf(
            "  [!]  Chaining: %s/%s (%.0f%%) above threshold -- MODERATE. Inspect flagged clusters above.",
            n_above, n_clusters, pct_chain
        ))
    } else {
        log_message(sprintf(
            "  Chaining: %s/%s (%.0f%%) above threshold -- minimal, acceptable.",
            n_above, n_clusters, pct_chain
        ))
    }

    # --- Membership skew warning ---
    if (max(nm) > 10 && max(nm) / stats::median(nm) > 5) {
        log_message(sprintf(
            "  [i]  Skewed cluster sizes: max members (%.0f) >> median (%.0f). ",
            max(nm), stats::median(nm)
        ))
        log_message("  A few super-clusters may absorb independent events. Inspect the largest clusters.")
    }

    log_message("--- End Post-Clustering Diagnosis ---")

    invisible(list(
        n_clusters = n_clusters,
        n_input_loops = n_input_loops,
        n_survived_loops = total_survived,
        members_median = stats::median(nm),
        members_max = max(nm),
        span_median = stats::median(max_span),
        span_max = max(max_span),
        n_above_chain_threshold = n_above,
        pct_above_chain_threshold = pct_chain
    ))
}

#' Internal: Core clustering for one chromosome-pair subset
#' @keywords internal
#' @noRd
.cluster_loops_inner <- function(dt, gap) {
    n_loops <- nrow(dt)
    if (n_loops == 0) return(dt)
    # Remap idx to local range (1..n_loops) so igraph vertex IDs stay in bounds.
    # The original idx values are from the global table before chr-pair splitting.
    dt[, idx := seq_len(n_loops)]

    hits <- dt[dt, on = .(
        chr1 = chr1,
        a1_l <= end1,
        a1_r >= start1,
        chr2 = chr2,
        a2_l <= end2,
        a2_r >= start2
    ), nomatch = NULL, allow.cartesian = TRUE]

    edges <- hits[idx < i.idx, .(from = idx, to = i.idx)]

    g <- igraph::make_empty_graph(n = n_loops, directed = FALSE)
    if (nrow(edges) > 0) {
        edge_vec <- as.vector(rbind(edges$from, edges$to))
        g <- igraph::add_edges(g, edge_vec)
    }

    comp <- igraph::components(g)
    dt[, cluster := comp$membership]
    dt
}

cluster_loops_dt <- function(dt, gap) {
    # All coordinates are 1-based closed (from gi_to_dt).
    # gap=0: only overlapping anchors merge (end1 >= start2 and end2 >= start1).
    #   Touching boundaries (end1 == start2) DO merge -- consistent with
    #   Bioconductor"s findOverlaps(maxgap=0) semantics.
    # gap=N: anchors within N bp merge (end1 + N >= start2).
    #
    # Clustering is batched by chromosome pair (chr1, chr2) to avoid O(N^2)
    # self-join on large datasets.  The join condition requires chr1 == chr1
    # and chr2 == chr2, so splitting is lossless.
    dt[, idx := .I]
    dt[, `:=`(
        a1_l = start1 - gap,
        a1_r = end1 + gap,
        a2_l = start2 - gap,
        a2_r = end2 + gap
    )]

    n_loops <- nrow(dt)
    if (n_loops > 50000) {
        warning(
            "Clustering ", n_loops, " loops with gap = ", gap, " bp. ",
            "Large datasets may require significant memory. ",
            "Consider pre-filtering or reducing gap.",
            call. = FALSE
        )
    }

    # Split by chromosome pair and cluster independently
    dt_list <- split(dt, by = c("chr1", "chr2"), drop = TRUE)
    n_chunks <- length(dt_list)
    if (n_chunks > 10) message("Clustering loops across ", n_chunks, " chromosome pairs...")
    clust_offset <- 0L
    result_list <- vector("list", length(dt_list))
    for (i in seq_along(dt_list)) {
        chunk <- .cluster_loops_inner(dt_list[[i]], gap)
        if (clust_offset > 0L) {
            chunk[, cluster := cluster + clust_offset]
        }
        clust_offset <- max(chunk$cluster)
        result_list[[i]] <- chunk
    }
    dt <- data.table::rbindlist(result_list)

    dt[, `:=`(idx = NULL, a1_l = NULL, a1_r = NULL, a2_l = NULL, a2_r = NULL)]
    dt
}

if (getRversion() >= "2.15.1") {
    utils::globalVariables(c(
        "V1", "V2", "V3", "V4", "V5", "V6", "V7",
        "chr1", "start1", "end1", "chr2", "start2", "end2",
        "idx", "i.idx", "cluster", "score", "source", "n_members", "n_reps",
        "a1_l", "a1_r", "a2_l", "a2_r", ".N", ".I", ".SD", ".SDcols", "..coord_cols"
    ))
}
