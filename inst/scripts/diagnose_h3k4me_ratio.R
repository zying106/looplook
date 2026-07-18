#!/usr/bin/env Rscript
#' Diagnose H3K4me1 / H3K4me3 Signal Ratio for Dual-Function Classification
#'
#' Computes overlap-width-weighted bigWig signal per anchor at dual-positive
#' loci (both H3K4me1 and H3K4me3 peaks), and at negative-control regions
#' (neither mark called) to estimate a data-driven ratio threshold.
#'
#' Usage:
#'   Rscript diagnose_h3k4me_ratio.R anchors.bed H3K4me1.bed H3K4me3.bed H3K4me1.bw H3K4me3.bw
#'   Rscript diagnose_h3k4me_ratio.R anchors.bed H3K4me1.bed H3K4me3.bed H3K4me1.bw H3K4me3.bw --k27me3 H3K27me3.bed --out plot.pdf
#'
#' With --k27me3: uses H3K27me3-only anchors (K27me3+ & K4me1- & K4me3-) as
#'   negative control.  Without: uses peak-negative anchors as background.
#'
#' Output:
#'   - Console: quantiles, threshold diagnostics, data-driven recommendation
#'   - PDF (if --out): dual-positive and background ratio distributions
#'
#' Required packages: rtracklayer, GenomicRanges, ggplot2
# ----------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript diagnose_h3k4me_ratio.R anchors.bed H3K4me1.bed H3K4me3.bed H3K4me1.bw H3K4me3.bw [--k27me3 H3K27me3.bed] [--out plot.pdf]")
}

bed_path   <- args[1]
me1_peaks  <- args[2]
me3_peaks  <- args[3]
bw_me1     <- args[4]
bw_me3     <- args[5]
k27me3_bed <- if ("--k27me3" %in% args) args[which(args == "--k27me3") + 1] else NULL
out_plot   <- if ("--out" %in% args) args[which(args == "--out") + 1] else NULL

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

# Overlap-width-weighted mean signal (matches .compute_bw_ratio in package)
weighted_mean_signal <- function(gr, sig) {
  v <- rep(NA_real_, length(gr))
  hits <- findOverlaps(gr, sig)
  if (length(hits) == 0) return(v)
  qh <- queryHits(hits); sh <- subjectHits(hits)
  w <- width(pintersect(gr[qh], sig[sh]))
  v[unique(qh)] <- tapply(sig$score[sh] * w, qh, sum) / tapply(w, qh, sum)
  v
}

# --- Load ---
message("Loading data ...")
anchors  <- read_bed(bed_path)
pk_me1   <- read_bed(me1_peaks)
pk_me3   <- read_bed(me3_peaks)

# --- Dual-positive anchors ---
h1 <- findOverlaps(anchors, pk_me1, maxgap = 200L)
h2 <- findOverlaps(anchors, pk_me3, maxgap = 200L)
dual_idx <- base::intersect(unique(queryHits(h1)), unique(queryHits(h2)))
gr_dual <- anchors[dual_idx]

message(sprintf("Total anchors:  %d", length(anchors)))
message(sprintf("H3K4me1+ only:  %d", length(setdiff(unique(queryHits(h1)), dual_idx))))
message(sprintf("H3K4me3+ only:  %d", length(setdiff(unique(queryHits(h2)), dual_idx))))
message(sprintf("Dual-positive:  %d", length(dual_idx)))

# --- Negative-control anchors ---
if (!is.null(k27me3_bed)) {
  pk_k27 <- read_bed(k27me3_bed)
  hk <- findOverlaps(anchors, pk_k27, maxgap = 200L)
  neg_idx <- setdiff(unique(queryHits(hk)), union(unique(queryHits(h1)), unique(queryHits(h2))))
  neg_label <- "H3K27me3-only"
} else {
  neg_idx <- setdiff(seq_along(anchors), union(unique(queryHits(h1)), unique(queryHits(h2))))
  neg_label <- "peak-negative"
}
gr_neg <- anchors[neg_idx]
message(sprintf("Negative control (%s): %d anchors", neg_label, length(neg_idx)))

# --- Weighted bigWig signal ---
compute_ratios <- function(gr, label) {
  if (length(gr) == 0) return(numeric(0))
  s1 <- import.bw(bw_me1, which = gr)
  s3 <- import.bw(bw_me3, which = gr)
  v1 <- weighted_mean_signal(gr, s1)
  v3 <- weighted_mean_signal(gr, s3)
  valid <- !is.na(v1) & !is.na(v3) & v1 > 0 & v3 > 0
  if (sum(valid) == 0) return(numeric(0))
  message(sprintf("  %s: %d / %d have valid signal", label, sum(valid), length(gr)))
  v1[valid] / v3[valid]
}

message("\nComputing weighted bigWig signal ratios ...")
ratio_dual   <- compute_ratios(gr_dual, "Dual-positive")
ratio_neg    <- compute_ratios(gr_neg,  neg_label)

if (length(ratio_dual) == 0) stop("No dual anchors had measurable signal for both marks.")

# --- Diagnostics ---
cat("\n=== Dual-Positive Anchor H3K4me1 / H3K4me3 Ratio ===\n")
print(round(quantile(ratio_dual, probs = c(0, 0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99, 1)), 3))

cat(sprintf("\n--- Fraction of dual-positive anchors at threshold ---\n"))
for (th in c(1, 1.5, 2, 3, 5)) {
  cat(sprintf("Ratio >= %.1f:  %.1f%%\n", th, 100 * mean(ratio_dual >= th)))
}

# --- Data-driven recommendation from negative control ---
if (length(ratio_neg) >= 10) {
  neg_p95 <- round(quantile(ratio_neg, 0.95, na.rm = TRUE), 2)
  cat(sprintf("\n=== Negative-control (%s) ratio distribution ===\n", neg_label))
  print(round(quantile(ratio_neg, probs = c(0, 0.50, 0.90, 0.95, 0.99, 1)), 3))
  cat(sprintf("\n--- Recommendation ---\n"))
  cat(sprintf("95th percentile of background ratio: %.2f\n", neg_p95))
  cat(sprintf("Suggested bw_ratio_threshold: %.1f\n",
      max(2, ceiling(neg_p95 * 10) / 10)))
  cat(sprintf("(i.e. only anchors with ratio > %.1f× background are likely true dual)\n",
      max(2, ceiling(neg_p95 * 10) / 10)))
  cat(sprintf("Current default bw_ratio_threshold = 3 would retain %.1f%% of dual-positives.\n",
      100 * mean(ratio_dual >= 3)))
} else {
  cat(sprintf("\n--- Recommendation ---\n"))
  cat(sprintf("Insufficient negative-control anchors (n=%d) for data-driven threshold.\n", length(ratio_neg)))
  cat(sprintf("Using default bw_ratio_threshold = 3.\n"))
}

# --- Plot ---
if (!is.null(out_plot)) {
  df_dual <- data.frame(log2_ratio = log2(ratio_dual + 1e-6), Group = "Dual-positive")
  df_all  <- df_dual
  if (length(ratio_neg) >= 10) {
    df_neg <- data.frame(log2_ratio = log2(ratio_neg + 1e-6), Group = neg_label)
    df_all <- rbind(df_dual, df_neg)
  }
  p <- ggplot(df_all, aes(x = log2_ratio, fill = Group)) +
    geom_histogram(aes(y = after_stat(density)), bins = 80,
                   alpha = 0.5, position = "identity", boundary = 0) +
    geom_density(aes(color = Group), linewidth = 1) +
    geom_vline(xintercept = log2(c(1, 2, 3)), linetype = "dashed",
               color = c("#999999", "#27AE60", "#E41A1C"), linewidth = 0.8) +
    annotate("text", x = log2(c(1, 2, 3)), y = Inf, vjust = 2,
             label = c("1x", "2x", "3x"),
             color = c("#999999", "#27AE60", "#E41A1C"), size = 3.5) +
    scale_x_continuous(breaks = seq(-5, 5, 1),
                       labels = c("1/32","1/16","1/8","1/4","1/2","1",
                                  "2","4","8","16","32")) +
    scale_fill_manual(values = c("Dual-positive" = "#5D6D7E", "#D95F02")) +
    scale_color_manual(values = c("Dual-positive" = "#E41A1C", "#D95F02")) +
    labs(x = expression(H3K4me1 / H3K4me3~(log[2])), y = "Density",
         title = "H3K4me1 / H3K4me3 Signal Ratio",
         subtitle = sprintf("Dual-positive (n=%d) vs %s background (n=%d)",
           length(ratio_dual), neg_label, length(ratio_neg))) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  ggsave(out_plot, p, width = 8, height = 5)
  message("\nPlot saved to ", out_plot)
}

message("\nDone.")
