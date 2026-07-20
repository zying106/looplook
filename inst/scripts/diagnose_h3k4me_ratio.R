#!/usr/bin/env Rscript
#' Diagnose H3K4me1 / H3K4me3 Signal Ratio at Dual-Positive Anchors
#'
#' Computes bigWig signal ratio using the same algorithm as
#' refine_loop_anchors_by_chromatin(): AUC / anchor_width with a data-driven
#' pseudocount, matching the package internal .compute_bw_ratio().
#'
#' Usage:
#'   Rscript diagnose_h3k4me_ratio.R anchors.bed H3K4me1.bed H3K4me3.bed \\
#'       H3K4me1.bw H3K4me3.bw [--anchor-gap 200] [--anchor-min-overlap 100] \\
#'       [--threshold 3] [--out plot.pdf]
#'
#' Output:
#'   - Console: n_dual_BED_positive, n_signal_resolved, pseudocount,
#'     quantiles, and % above specified threshold
#'   - PDF (if --out): histogram + density of log2 ratios
#'
#' Required packages: rtracklayer, GenomicRanges, ggplot2
# ----------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript diagnose_h3k4me_ratio.R anchors.bed H3K4me1.bed H3K4me3.bed H3K4me1.bw H3K4me3.bw [--anchor-gap 200] [--anchor-min-overlap 100] [--threshold 3] [--out plot.pdf]")
}

bed_path    <- args[1]
me1_peaks   <- args[2]
me3_peaks   <- args[3]
bw_me1      <- args[4]
bw_me3      <- args[5]

# Parse optional named args
named_integer <- function(key, default) {
  idx <- which(args == key)
  if (!length(idx)) return(default)
  if (idx[1] >= length(args)) stop("Missing value for ", key)
  value <- suppressWarnings(as.integer(args[idx[1] + 1L]))
  if (length(value) != 1L || is.na(value) || value < 0L) stop("Invalid non-negative integer for ", key)
  value
}
named_numeric <- function(key, default) {
  idx <- which(args == key)
  if (!length(idx)) return(default)
  value <- suppressWarnings(as.numeric(args[idx[1] + 1L]))
  if (length(value) != 1L || is.na(value) || !is.finite(value))
    stop("Invalid numeric value for ", key)
  value
}
anchor_gap         <- named_integer("--anchor-gap", 200L)
anchor_min_overlap <- named_integer("--anchor-min-overlap", 100L)
threshold          <- named_numeric("--threshold", 3)
out_plot           <- if ("--out" %in% args) args[which(args == "--out") + 1] else NULL

suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
  library(ggplot2)
})

read_bed <- function(f) {
  df <- read.table(f, sep = "\t", header = FALSE, comment.char = "#")
  colnames(df)[1:3] <- c("chr", "start", "end")
  df$start <- df$start + 1L
  GRanges(df$chr, IRanges(df$start, df$end))
}

# --- Load anchors and peak files ---
message("Loading data ...")
anchors  <- read_bed(bed_path)
pk_me1   <- read_bed(me1_peaks)
pk_me3   <- read_bed(me3_peaks)

# --- Overlap logic: same cascade as refine_loop_anchors_by_chromatin ---
# Step 1: find candidates within gap tolerance
h1_raw <- findOverlaps(anchors, pk_me1, maxgap = anchor_gap)
h2_raw <- findOverlaps(anchors, pk_me3, maxgap = anchor_gap)
# Step 2: require minimum physical overlap (anchor_min_overlap)
ol1 <- pintersect(anchors[queryHits(h1_raw)], pk_me1[subjectHits(h1_raw)])
ol2 <- pintersect(anchors[queryHits(h2_raw)], pk_me3[subjectHits(h2_raw)])
keep1 <- width(ol1) >= anchor_min_overlap
keep2 <- width(ol2) >= anchor_min_overlap
h1 <- h1_raw[keep1]
h2 <- h2_raw[keep2]

me1_hits <- unique(queryHits(h1))
me3_hits <- unique(queryHits(h2))
dual_idx <- sort(intersect(me1_hits, me3_hits))
me1_only <- setdiff(me1_hits, dual_idx)
me3_only <- setdiff(me3_hits, dual_idx)

message(sprintf("Total anchors:        %d", length(anchors)))
message(sprintf("H3K4me1+ only:       %d", length(me1_only)))
message(sprintf("H3K4me3+ only:       %d", length(me3_only)))
message(sprintf("Dual-positive (BED): %d", length(dual_idx)))

if (length(dual_idx) == 0) stop("No dual-positive anchors found; nothing to diagnose.")

# --- bigWig signal: AUC / anchor_width (matching .compute_bw_ratio) ---
gr_dual <- anchors[dual_idx]
anchor_width <- GenomicRanges::width(gr_dual)

s1 <- import.bw(bw_me1, which = gr_dual)
s3 <- import.bw(bw_me3, which = gr_dual)

.anchor_mean <- function(gr_dual, sig, anchor_width) {
  hits <- findOverlaps(gr_dual, sig)
  n <- length(gr_dual)
  anchor_mean <- rep(NA_real_, n)
  if (length(hits) == 0) return(anchor_mean)
  qh <- queryHits(hits); sh <- subjectHits(hits)
  ol <- pintersect(gr_dual[qh], sig[sh])
  score_width <- sig$score[sh] * width(ol)
  hit_ids <- unique(qh)
  auc <- tapply(score_width, qh, sum)
  anchor_mean[hit_ids] <- as.numeric(auc) / anchor_width[hit_ids]
  anchor_mean
}

mean1 <- .anchor_mean(gr_dual, s1, anchor_width)
mean3 <- .anchor_mean(gr_dual, s3, anchor_width)

# Data-driven pseudocount: 1st percentile of ALL positive anchor-means
# (matching .compute_bw_ratio() exactly)
pooled <- c(mean1, mean3)
pooled_pos <- pooled[pooled > 0 & !is.na(pooled)]
pseudocount <- if (length(pooled_pos) >= 10) {
  max(quantile(pooled_pos, 0.01, na.rm = TRUE, names = FALSE), 1e-6)
} else {
  1e-6
}

log2_ratio <- log2((mean1 + pseudocount) / (mean3 + pseudocount))
resolved <- !is.na(log2_ratio) & is.finite(log2_ratio)
n_resolved <- sum(resolved)

message(sprintf("Anchors with resolved signal: %d", n_resolved))
message(sprintf("pseudocount: %g", pseudocount))
if (n_resolved == 0) stop("No dual anchors had positive signal for both marks.")

# --- Diagnostics ---
cat(sprintf("\n=== H3K4me1 / H3K4me3 (anchor_mean, pseudocount = %g, threshold = %g) ===\n",
    pseudocount, threshold))
r <- log2_ratio[resolved]
print(round(quantile(r, probs = c(0, 0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99, 1)), 3))

cat(sprintf("\n--- Threshold diagnostics ---\n"))
for (t in c(1, 1.5, 2, 3, 4, 5)) {
  pct <- 100 * mean(r >= log2(t))
  cat(sprintf("Ratio >= %4.1f:  %.1f%%%s\n", t, pct, if (t == threshold) "  (default threshold)" else ""))
}
cat(sprintf("Below threshold (< %g): %.1f%% (classified as P/not dual)\n",
    threshold, 100 * mean(r < log2(threshold))))

# --- Plot ---
if (!is.null(out_plot)) {
  lt <- log2(threshold)
  df <- data.frame(log2_ratio = r)
  p <- ggplot(df, aes(x = log2_ratio)) +
    geom_histogram(aes(y = after_stat(density)), bins = 80,
                   fill = "#5D6D7E", alpha = 0.7, boundary = 0) +
    geom_density(color = "#E41A1C", linewidth = 1) +
    geom_vline(xintercept = lt, linetype = "dashed",
               color = "#27AE60", linewidth = 0.8) +
    annotate("text", x = lt, y = Inf, vjust = 2,
             label = sprintf("threshold=%g", threshold),
             color = "#27AE60", size = 3.5) +
    scale_x_continuous(breaks = seq(-10, 10, 1)) +
    labs(x = expression(H3K4me1 / H3K4me3~(log[2])),
         y = "Density",
         title = sprintf("H3K4me1 / H3K4me3 Signal Ratio\n(%d resolved dual-positive anchors)", n_resolved),
         subtitle = sprintf("AUC/anchor_width, pseudocount = %g", pseudocount)) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  ggsave(out_plot, p, width = 8, height = 5)
  message("\nPlot saved to ", out_plot)
}

message("\nDone.")
