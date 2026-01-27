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
#'
#' @param bedpe_file Character. Path to a BEDPE file. Must contain at least six columns:
#'   \code{chr1, start1, end1, chr2, start2, end2}.
#' @return A \code{\link[InteractionSet]{GInteractions}} object with a \code{score} metadata column
#'   (defaulting to 0 if not provided).
#' @importFrom data.table fread
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom InteractionSet GInteractions
#' @importFrom S4Vectors mcols<-
#' @export
#' @examples
#' # 1. Locate the example BEDPE file included in the package
#' # system.file finds the absolute path to 'inst/extdata/example_loops_1.bedpe'
#' bedpe_path <- system.file("extdata", "example_loops_1.bedpe", package = "looplook")
#'
#' # 2. Run the function (ensure file was found)
#' if (bedpe_path != "") {
#'   gi <- bedpe_to_gi(bedpe_path)
#'
#'   # 3. Inspect the result
#'   print(gi)
#'
#'   # Check the imported score column
#'   S4Vectors::mcols(gi)$score
#' }
bedpe_to_gi <- function(bedpe_file) {
  # Robust read with fallback
  df <- tryCatch(
    {
      data.table::fread(bedpe_file, header = FALSE, sep = "\t")
    },
    error = function(e) {
      data.table::fread(bedpe_file, header = FALSE)
    }
  )

  if (ncol(df) < 6) {
    stop("BEDPE file must have at least 6 columns: ", bedpe_file)
  }

  df <- as.data.frame(df)

  # Normalize anchor order: ensure (chr1, start1) <= (chr2, start2)
  swap <- (df[, 1] > df[, 4]) | (df[, 1] == df[, 4] & df[, 2] > df[, 5])
  if (any(swap)) {
    df[swap, seq_len(6)] <- df[swap, c(4, 5, 6, 1, 2, 3)]
  }

  gr1 <- GenomicRanges::GRanges(df[, 1], IRanges::IRanges(df[, 2], df[, 3]))
  gr2 <- GenomicRanges::GRanges(df[, 4], IRanges::IRanges(df[, 5], df[, 6]))

  gi <- InteractionSet::GInteractions(gr1, gr2, mode = "strict")

  # Smart Score Detection
  final_scores <- rep(0, nrow(df))
  found <- FALSE

  if (ncol(df) >= 8) {
    v8 <- as.numeric(df[, 8])
    if (sum(!is.na(v8)) > (nrow(df) * 0.5)) {
      final_scores <- v8
      found <- TRUE
    }
  }
  if (!found && ncol(df) >= 7) {
    v7 <- as.numeric(df[, 7])
    if (sum(!is.na(v7)) > (nrow(df) * 0.5)) {
      final_scores <- v7
    }
  }

  final_scores[is.na(final_scores)] <- 0
  S4Vectors::mcols(gi)$score <- final_scores

  return(gi)
}


#' Spatial Clustering of GInteractions
#'
#' Merges spatially proximal chromatin loops into consensus interactions using graph-based clustering.
#' Loops are considered overlapping if both anchors are within \code{gap} bp of each other.
#' Each resulting cluster is represented by the **union genomic range (min start to max end)** spanning all its members.
#'
#' @param gi A \code{\link[InteractionSet]{GInteractions}} object.
#' @param gap Numeric. Maximum distance (in base pairs) allowed between anchors to consider two loops overlapping. Default: 1000.
#' @return A list with two elements:
#'   \itemize{
#'     \item{\code{gi}}{Reduced \code{GInteractions} object, one per cluster.}
#'     \item{\code{membership}}{Integer vector indicating cluster assignment for each input loop.}
#'   }
#'   Metadata columns include \code{cluster_id}, \code{n_members}, and averaged \code{score}.
#' @importFrom igraph make_empty_graph add_edges components
#' @importFrom S4Vectors queryHits subjectHits mcols
#' @importFrom GenomicRanges GRangesList
#' @export
#' @examples
#' # 1. Load example data (loops that are close to each other)
#' bedpe_path <- system.file("extdata", "example_loops_1.bedpe", package = "looplook")
#'
#' if (bedpe_path != "") {
#'   # Convert BEDPE to GInteractions object
#'   gi_raw <- bedpe_to_gi(bedpe_path)
#'
#'   # 2. Run clustering
#'   # Merge loops if their anchors are within 1000bp
#'   res <- reduce_ginteractions(gi_raw, gap = 1000)
#'
#'   # 3. Inspect results
#'   # The 'gi' element contains the merged consensus loops
#'   print(res$gi)
#'
#'   # The 'membership' vector tells which original loop belongs to which cluster
#'   head(res$membership)
#'
#'   # Check cluster sizes (how many loops were merged into each cluster)
#'   table(res$membership)
#' }
reduce_ginteractions <- function(gi, gap = 1000) {
  if (length(gi) == 0) {
    return(list(gi = gi, membership = integer(0)))
  }

  dt <- gi_to_dt(gi)
  dt <- cluster_loops_dt(dt, gap)
  reduced_dt <- reduce_clusters_dt(dt)
  gi_red <- dt_to_gi(reduced_dt)

  list(gi = gi_red, membership = dt$cluster)
}

#' Read a Simple BED File into a GRanges Object
#'
#' Reads the first three columns of a BED file (chrom, start, end) and returns a
#' \code{\link[GenomicRanges]{GRanges}} object. Additional columns are ignored.
#'
#' @param bed_file Character. Path to a BED file (must be tab-delimited).
#' @return A \code{\link[GenomicRanges]{GRanges}} object, or \code{NULL} if \code{bed_file} is \code{NULL}.
#' @importFrom data.table fread
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @export
#' @examples
#' # 1. Locate the example BED file included in the package
#' bed_path <- system.file("extdata", "example_peaks.bed", package = "looplook")
#'
#' # 2. Run the function (ensure file was found)
#' if (bed_path != "") {
#'   # Read BED file into a GRanges object
#'   gr <- read_simple_bed(bed_path)
#'
#'   # 3. Inspect the result
#'   print(gr)
#'
#'   # Check how many peaks were loaded
#'   length(gr)
#' }
read_simple_bed <- function(bed_file) {
  if (is.null(bed_file)) {
    return(NULL)
  }
  if (!file.exists(bed_file)) {
    stop("BED file does not exist: ", bed_file)
  }

  df <- data.table::fread(bed_file, header = FALSE, select = c(1, 2, 3))
  GenomicRanges::GRanges(df$V1, IRanges::IRanges(df$V2, df$V3))
}


#' Consolidate and Integrate Chromatin Loops from Multiple Sources
#'
#' @description
#' This function consolidates chromatin loops from multiple BEDPE files. It is designed for two main purposes:
#' \enumerate{
#'   \item **Replicate Consolidation**: Merging biological or technical replicates to identify high-confidence, reproducible loops (e.g., 3 replicates of H3K27ac HiChIP).
#'   \item **Multi-Omics Integration**: Intersecting or merging loops from distinct chromatin features to find shared regulatory architectures (e.g., integrating **H3K27ac** and **H3K4me3** HiChIP data, or overlapping Hi-C with ChIA-PET).
#' }
#'
#' The function supports three modes:
#' \itemize{
#'   \item \code{"reproducible"}: Graph-based clustering to find a consensus set supported by a majority of samples.
#'   \item \code{"intersect"}: Strict reference-based filtering (keeps loops in File 1 supported by ALL other files).
#'   \item \code{"union"}: Merges all detected loops into a comprehensive map.
#' }
#' Optional filtering by blacklist regions or target genomic regions (e.g., promoters) is supported.
#'
#' @param files Character vector. Paths to BEDPE files (at least two).
#' @param gap Numeric. Distance (bp) to consider loops as overlapping. Default 1000.
#' @param mode Character. Choose one of the following: intersect, union, reproducible. Merge strategy:
#'   \itemize{
#'     \item{"intersect"}{ Strict reference-based filtering (keeps loops in File 1 supported by ALL other files)}
#'     \item{"union"}{ Merges all detected loops into a comprehensive map.}
#'     \item{"reproducible"}{ Graph-based clustering to find a consensus set supported by a majority of samples.}
#'   }
#' @param min_reproducible Integer. Minimum number of replicates a loop must appear in
#'   (only effective when \code{mode = "reproducible"}).
#'   If \code{NULL} (default), the threshold is automatically calculated:
#'   \itemize{
#'     \item For 2 replicates: Requires both (2/2).
#'     \item For >2 replicates: Requires strict majority (>75% support).
#'     (e.g., 3 for N=3, 4 for N=4, 4 for N=5).
#'   }
#' @param min_score Numeric. Minimum score threshold to keep a loop.
#' @param blacklist_species Character. Species/build for built-in blacklist
#'   (e.g., "hg38", "hg19", "mm10", "mm9"), or a path to a custom BED file.
#' @param target_bed_file Character. Path to BED file. Only loops overlapping these regions will be kept.
#' @param out_file Character. The file name (including the file path) for saving results in the extended BEDPE format.
#' @return A filtered \code{\link[InteractionSet]{GInteractions}} object.
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
#' if (f1 != "" && f2 != "") {
#'   # Example A: Intersect Mode
#'   # Only keeps loops present in f1 that are also supported by f2
#'   res_intersect <- consolidate_chromatin_loops(
#'     files = c(f1, f2),
#'     mode = "intersect",
#'     gap = 1000,
#'     out_file = tempfile(fileext = ".bedpe")
#'   )
#'
#'   # Example B: Reproducible Mode
#'   # Finds consensus loops supported by both replicates (default for N=2)
#'   res_repro <- consolidate_chromatin_loops(
#'     files = c(f1, f2),
#'     mode = "reproducible",
#'     gap = 1000,
#'     out_file = tempfile(fileext = ".bedpe")
#'   )
#'
#'   # Example C: Union Mode
#'   # Merges all loops into a single map
#'   res_union <- consolidate_chromatin_loops(
#'     files = c(f1, f2),
#'     mode = "union",
#'     gap = 1000,
#'     out_file = tempfile(fileext = ".bedpe")
#'   )
#'
#'   # Inspect results
#'   length(res_intersect)
#'   length(res_union)
#' }
consolidate_chromatin_loops <- function(
  files = NULL,
  gap = 1000,
  mode = c("reproducible", "intersect", "union"),
  min_reproducible = NULL,
  min_score = NULL,
  blacklist_species = NULL,
  target_bed_file = NULL,
  out_file = NULL
) {
  stopifnot(length(files) >= 2)
  mode <- match.arg(mode)
  n_reps <- length(files)

  message(">>> Reading BEDPE files")
  gi_list <- lapply(files, bedpe_to_gi)

  for (i in seq_along(gi_list)) {
    S4Vectors::mcols(gi_list[[i]])$source <- i
    message("    File", i, ": ", length(gi_list[[i]]), " loops")
  }

  result_gi <- NULL

  # ==========================================
  # PATH A: STRICT INTERSECT (No Clustering)
  # ==========================================
  if (mode == "intersect") {
    message(">>> Intersect mode: Reference-based filtering (No Coordinate Merging)")
    message("    Base: File 1. Criterion: Must overlap with ALL other files.")

    # Start with File 1 as the reference universe
    current_gi <- gi_list[[1]]

    # Iteratively filter against File 2, 3, ... N
    for (i in 2:n_reps) {
      if (length(current_gi) == 0) break

      message("    Intersecting with File ", i, "...")

      # Use S4 findOverlaps for robust gap checking
      hits <- InteractionSet::findOverlaps(
        current_gi,
        gi_list[[i]],
        maxgap = gap,
        use.region = "both"
      )

      # Keep only loops in current_gi that have a hit
      keep_idx <- unique(S4Vectors::queryHits(hits))
      current_gi <- current_gi[keep_idx]
    }

    result_gi <- current_gi

    # In intersect mode, n_members/n_reps is logically N (since we filtered for it)
    S4Vectors::mcols(result_gi)$n_reps <- n_reps
    # Set n_members to 1 (representing the single representative from File 1) or N?
    # Let's use N to be consistent with the logic that N supports it.
    S4Vectors::mcols(result_gi)$n_members <- n_reps

    # ==========================================
    # PATH B: CLUSTERING (Reproducible / Union)
    # ==========================================
  } else {
    message(">>> Clustering mode (Union/Reproducible): Merging coordinates via Graph")

    # 1. Combine
    combined_dt <- data.table::rbindlist(lapply(gi_list, gi_to_dt))

    # 2. Cluster
    clustered <- cluster_loops_dt(combined_dt, gap)

    # 3. Reduce
    reduced_dt <- reduce_clusters_dt(clustered)

    # 4. Filter Logic
    if (mode == "reproducible") {
      if (is.null(min_reproducible)) {
        if (n_reps == 2) {
          min_reproducible <- 2
        } else {
          # > 75% logic
          min_reproducible <- floor(0.75 * n_reps) + 1
        }
      }
      message(">>> Reproducible mode: Keeping clusters in >= ", min_reproducible, " replicates")
      reduced_dt <- reduced_dt[n_reps >= min_reproducible]
    } else {
      message(">>> Union mode: Keeping all clusters")
    }

    result_gi <- dt_to_gi(reduced_dt)
  }

  # ==========================================
  # POST-PROCESSING (Filters & Output)
  # ==========================================

  # Score Filter
  if (!is.null(min_score)) {
    keep <- S4Vectors::mcols(result_gi)$score > min_score
    result_gi <- result_gi[keep]
  }

  # Blacklist Filter
  if (!is.null(blacklist_species)) {
    known_lists <- list(
      "hg38" = "hg38-blacklist.v2.bed",
      "hg19" = "hg19-blacklist.v2.bed",
      "mm10" = "mm10-blacklist.v2.bed",
      "mm9"  = "mm9-blacklist.v2.bed"
    )
    blacklist_path <- NULL
    if (blacklist_species %in% names(known_lists)) {
      blacklist_path <- system.file("extdata", known_lists[[blacklist_species]], package = "looplook")
    }
    if (is.null(blacklist_path) || blacklist_path == "") {
      blacklist_path <- blacklist_species
    }
    if (file.exists(blacklist_path)) {
      message(">>> Filtering blacklist: ", basename(blacklist_path))
      bl <- read_simple_bed(blacklist_path)
      h1 <- InteractionSet::findOverlaps(InteractionSet::anchors(result_gi, "first"), bl)
      h2 <- InteractionSet::findOverlaps(InteractionSet::anchors(result_gi, "second"), bl)
      bad <- unique(c(S4Vectors::queryHits(h1), S4Vectors::queryHits(h2)))
      if (length(bad)) result_gi <- result_gi[-bad]
    } else {
      warning("Blacklist file not found: ", blacklist_species)
    }
  }

  # Target Region Filter
  if (!is.null(target_bed_file)) {
    message(">>> Filtering targets: ", basename(target_bed_file))
    if (file.exists(target_bed_file)) {
      tg <- read_simple_bed(target_bed_file)
      h1 <- InteractionSet::findOverlaps(InteractionSet::anchors(result_gi, "first"), tg)
      h2 <- InteractionSet::findOverlaps(InteractionSet::anchors(result_gi, "second"), tg)
      keep <- unique(c(S4Vectors::queryHits(h1), S4Vectors::queryHits(h2)))
      result_gi <- result_gi[keep]
    } else {
      warning("Target BED file not found.")
    }
  }

  # Output to file
  if (!is.null(out_file)) {
    a1 <- InteractionSet::anchors(result_gi, "first")
    a2 <- InteractionSet::anchors(result_gi, "second")

    out_df <- data.frame(
      chr1 = GenomicRanges::seqnames(a1),
      start1 = GenomicRanges::start(a1),
      end1 = GenomicRanges::end(a1),
      chr2 = GenomicRanges::seqnames(a2),
      start2 = GenomicRanges::start(a2),
      end2 = GenomicRanges::end(a2),
      name = ".",
      score = round(S4Vectors::mcols(result_gi)$score, 2),
      n_members = if (!is.null(S4Vectors::mcols(result_gi)$n_reps)) S4Vectors::mcols(result_gi)$n_reps else 1,
      stringsAsFactors = FALSE
    )

    dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
    utils::write.table(out_df, out_file,
      sep = "\t",
      row.names = FALSE, col.names = FALSE, quote = FALSE
    )
    message("Finished! Saved to ", out_file)
  }

  message("inished! Final loops: ", length(result_gi))
  result_gi
}


# --- Helpers ---

gi_to_dt <- function(gi) {
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
  dt[, list(
    chr1 = chr1[1],
    start1 = min(start1),
    end1 = max(end1),
    chr2 = chr2[1],
    start2 = min(start2),
    end2 = max(end2),
    score = mean(score, na.rm = TRUE),
    n_members = .N,
    n_reps = data.table::uniqueN(source)
  ), by = cluster]
}

dt_to_gi <- function(dt) {
  gr1 <- GenomicRanges::GRanges(
    dt$chr1, IRanges::IRanges(dt$start1, dt$end1)
  )
  gr2 <- GenomicRanges::GRanges(
    dt$chr2, IRanges::IRanges(dt$start2, dt$end2)
  )

  gi <- InteractionSet::GInteractions(gr1, gr2, mode = "strict")
  S4Vectors::mcols(gi)$cluster_id <- dt$cluster
  S4Vectors::mcols(gi)$n_members <- dt$n_members
  S4Vectors::mcols(gi)$score <- dt$score
  S4Vectors::mcols(gi)$n_reps <- dt$n_reps
  gi
}

cluster_loops_dt <- function(dt, gap) {
  dt[, idx := .I]
  dt[, `:=`(
    a1_l = start1 - gap,
    a1_r = end1 + gap,
    a2_l = start2 - gap,
    a2_r = end2 + gap
  )]

  hits <- dt[dt, on = .(
    chr1 = chr1,
    a1_l <= start1,
    a1_r >= end1,
    chr2 = chr2,
    a2_l <= start2,
    a2_r >= end2
  ), nomatch = NULL, allow.cartesian = TRUE]

  edges <- hits[idx < i.idx, .(from = idx, to = i.idx)]

  g <- igraph::make_empty_graph(n = nrow(dt), directed = FALSE)
  if (nrow(edges) > 0) {
    edge_vec <- as.vector(t(as.matrix(edges)))
    g <- igraph::add_edges(g, edge_vec)
  }

  comp <- igraph::components(g)
  dt[, cluster := comp$membership]
  dt[, `:=`(idx = NULL, a1_l = NULL, a1_r = NULL, a2_l = NULL, a2_r = NULL)]
  dt
}

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "V1", "V2", "V3", "V4", "V5", "V6", "V7",
    "chr1", "start1", "end1", "chr2", "start2", "end2",
    "idx", "i.idx", "cluster", "score", "source", "n_members", "n_reps",
    "a1_l", "a1_r", "a2_l", "a2_r", ".N", ".I"
  ))
}

#' Plot Contact Matrix Heatmap
#'
#' Visualizes chromatin interactions (GInteractions) as a rotated triangular heatmap.
#' Supports custom shapes (diamond/triangle) with perfect boundary alignment.
#'
#' @param gi A \code{\link[InteractionSet]{GInteractions}} object.
#' @param region Character. Genomic region to plot (e.g., "chr1:100000-200000").
#' @param resolution Integer. Bin size in base pairs. Default 10000.
#' @param palette Character vector. Color gradient.
#' @param max_score Numeric. Cap the score at this value.
#' @param ylim_bins Numeric. Limit Y-axis height (zoom in).
#' @param shape Character. Shape of the pixels: "diamond" (standard Hi-C), "triangle", or "rect". Default "diamond".
#' @return A ggplot object.
#' @importFrom GenomicRanges GRanges start end seqnames
#' @importFrom InteractionSet findOverlaps anchors
#' @importFrom data.table data.table setkey .N :=
#' @importFrom ggplot2 ggplot geom_tile geom_polygon geom_point geom_path scale_fill_gradientn scale_color_gradientn scale_x_continuous theme_classic theme element_blank element_text coord_fixed aes margin ggtitle guides guide_colorbar
#' @importFrom S4Vectors queryHits subjectHits mcols
#' @export
#' @examples
#' # 1. Load example data
#' bedpe_path <- system.file("extdata", "example_loops_1.bedpe", package = "looplook")
#'
#' if (bedpe_path != "") {
#'   # Convert BEDPE to GInteractions object
#'   gi <- bedpe_to_gi(bedpe_path)
#'
#'   # 2. Define a region of interest
#'   region_str <- "chr1:9893482-13450001"
#'
#'   # 3. Plot the contact heatmap
#'   # We use a resolution of 10kb to match the scale of the example loops
#'   p <- plot_contact_matrix(
#'     gi = gi,
#'     region = region_str,
#'     resolution = 5000,
#'     shape = "diamond"
#'   )
#'
#'   print(p)
#' }
plot_contact_matrix <- function(gi, region, resolution = 10000,
                                palette = c("white", "orange", "red"),
                                max_score = NULL,
                                ylim_bins = NULL,
                                shape = c("diamond", "triangle", "rect")) {
  shape <- match.arg(shape)

  # 1. Parse Region
  reg_gr <- .parse_region(region)
  chr_name <- as.character(GenomicRanges::seqnames(reg_gr))
  reg_start <- GenomicRanges::start(reg_gr)
  reg_end <- GenomicRanges::end(reg_gr)

  message(">>> Plotting region: ", region, " @ ", resolution / 1000, "kb | Shape: ", shape)

  # 2. Filter Loops
  hits <- InteractionSet::findOverlaps(gi, reg_gr, type = "within", use.region = "both")
  if (length(hits) == 0) {
    warning("No loops found.")
    return(NULL)
  }
  sub_gi <- gi[unique(S4Vectors::queryHits(hits))]

  # 3. Binning
  a1 <- InteractionSet::anchors(sub_gi, "first")
  a2 <- InteractionSet::anchors(sub_gi, "second")
  scores <- if (!is.null(S4Vectors::mcols(sub_gi)$score)) S4Vectors::mcols(sub_gi)$score else 1

  dt <- data.table::data.table(
    s1 = GenomicRanges::start(a1), e1 = GenomicRanges::end(a1),
    s2 = GenomicRanges::start(a2), e2 = GenomicRanges::end(a2),
    score = scores
  )
  dt[, `:=`(bin1 = floor(((s1 + e1) / 2 - reg_start) / resolution), bin2 = floor(((s2 + e2) / 2 - reg_start) / resolution))]
  dt[bin1 > bin2, c("bin1", "bin2") := list(bin2, bin1)]

  max_bin <- floor((reg_end - reg_start) / resolution)
  dt <- dt[bin1 >= 0 & bin2 >= 0 & bin1 <= max_bin & bin2 <= max_bin]

  # Aggregate
  mat_dt <- dt[, .(value = sum(score)), by = .(bin1, bin2)]

  # 4. Coordinate Calculation
  mat_dt[, `:=`(x_rot = (bin1 + bin2) / 2, y_rot = (bin2 - bin1) / 2)]

  if (!is.null(max_score)) mat_dt[value > max_score, value := max_score]

  # === Shape Logic ===
  plot_layer <- NULL

  if (shape == "diamond") {
    mat_dt[, id := .I]
    poly_dt <- mat_dt[, .(
      x = c(x_rot, x_rot + 0.5, x_rot, x_rot - 0.5),
      y = c(y_rot + 0.5, y_rot, y_rot - 0.5, y_rot),
      value = value
    ), by = id]

    plot_layer <- ggplot2::geom_polygon(
      data = poly_dt,
      ggplot2::aes(x = x, y = y, fill = value, group = id)
    )
  } else if (shape == "triangle") {
    plot_layer <- ggplot2::geom_point(
      data = mat_dt,
      ggplot2::aes(x = x_rot, y = y_rot, color = value),
      shape = 17, size = 1.5
    )
  } else {
    plot_layer <- ggplot2::geom_tile(
      data = mat_dt,
      ggplot2::aes(x = x_rot, y = y_rot, fill = value),
      width = 1, height = 1
    )
  }

  # === 5. Border (Fix Overlap) ===
  # Expand border by 0.5 units to wrap AROUND the diamonds
  expand <- 0.5
  triangle_border <- data.frame(
    x = c(-expand, max_bin + expand, max_bin / 2, -expand),
    y = c(0, 0, max_bin / 2 + expand, 0)
  )

  max_h <- if (!is.null(ylim_bins)) ylim_bins else (max_bin / 2 + expand)
  coord_system <- ggplot2::coord_fixed(ratio = 0.5, ylim = c(0, max_h))

  # 6. Assemble Plot
  p <- ggplot2::ggplot() +
    plot_layer +
    ggplot2::geom_path(data = triangle_border, ggplot2::aes(x = x, y = y), color = "black", size = 0.5) +
    ggplot2::scale_fill_gradientn(colors = palette, name = "Score") +
    ggplot2::scale_color_gradientn(colors = palette, name = "Score") +
    ggplot2::scale_x_continuous(
      name = paste0(chr_name, " (Mb)"),
      breaks = seq(0, max_bin, length.out = 5),
      labels = function(b) sprintf("%.2f", (reg_start + b * resolution) / 1e6)
    ) +
    coord_system +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.line.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5),
      legend.position = "right"
    ) +
    ggplot2::ggtitle(paste0(chr_name, ":", prettyNum(reg_start, big.mark = ","), "-", prettyNum(reg_end, big.mark = ",")))

  return(p)
}

# Helper
.parse_region <- function(region_str) {
  parts <- strsplit(region_str, ":")[[1]]
  if (length(parts) != 2) stop("Invalid region format. Use 'chr:start-end'")
  chrom <- parts[1]
  coords <- strsplit(parts[2], "-")[[1]]
  if (length(coords) != 2) stop("Invalid coordinates. Use 'chr:start-end'")
  return(GenomicRanges::GRanges(chrom, IRanges::IRanges(as.numeric(coords[1]), as.numeric(coords[2]))))
}

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "s1", "e1", "s2", "e2", "bin1", "bin2", "value", "x_rot", "y_rot", "x", "y", "id"
  ))
}
