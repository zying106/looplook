#!/usr/bin/env Rscript
#' Diagnose H3K4me1 / H3K4me3 Signal Ratio at Dual-Positive Anchors
#'
#' Computes bigWig signal ratio only for anchors where BOTH H3K4me1 and
#' H3K4me3 peaks overlap (i.e., dual candidates in the BED-based pipeline).
#' These are the only anchors that need bigWig to resolve dual vs P.
#'
#' Usage:
#'   Rscript diagnose_h3k4me_ratio.R anchors.bed H3K4me1.bed H3K4me3.bed H3K4me1.bw H3K4me3.bw
#'   Rscript diagnose_h3k4me_ratio.R anchors.bed H3K4me1.bed H3K4me3.bed H3K4me1.bw H3K4me3.bw --out plot.pdf
#'
#' Output:
#'   - Console: quantiles + threshold diagnostics (dual candidates only)
#'   - PDF (if --out): histogram + density of log2 ratios
#'
#' Required packages: rtracklayer, GenomicRanges, ggplot2
# ----------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript diagnose_h3k4me_ratio.R anchors.bed H3K4me1.bed H3K4me3.bed H3K4me1.bw H3K4me3.bw [--out plot.pdf]")
}

bed_path   <- args[1]
me1_peaks  <- args[2]
me3_peaks  <- args[3]
bw_me1     <- args[4]
bw_me3     <- args[5]
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

# --- Load anchors and peak files ---
message("Loading data ...")
anchors  <- read_bed(bed_path)
pk_me1   <- read_bed(me1_peaks)
pk_me3   <- read_bed(me3_peaks)

# --- Find dual-positive anchors (both marks overlap) ---
h1 <- findOverlaps(anchors, pk_me1, maxgap = 200L)
h2 <- findOverlaps(anchors, pk_me3, maxgap = 200L)
dual_idx <- base::intersect(unique(queryHits(h1)), unique(queryHits(h2)))
gr_dual <- anchors[dual_idx]

message(sprintf("Total anchors:  %d", length(anchors)))
message(sprintf("H3K4me1+ only:  %d", length(unique(queryHits(h1)))))
message(sprintf("H3K4me3+ only:  %d", length(unique(queryHits(h2)))))
message(sprintf("Dual-positive:  %d", length(dual_idx)))

if (length(dual_idx) == 0) stop("No dual-positive anchors found; nothing to diagnose.")

# --- Extract bigWig signal at dual anchors ---
s1 <- import.bw(bw_me1, which = gr_dual)
s3 <- import.bw(bw_me3, which = gr_dual)

hits1 <- findOverlaps(gr_dual, s1)
hits3 <- findOverlaps(gr_dual, s3)

sig1 <- rep(NA_real_, length(gr_dual))
sig3 <- rep(NA_real_, length(gr_dual))
if (length(hits1) > 0) sig1[unique(queryHits(hits1))] <- tapply(s1$score[subjectHits(hits1)], queryHits(hits1), mean)
if (length(hits3) > 0) sig3[unique(queryHits(hits3))] <- tapply(s3$score[subjectHits(hits3)], queryHits(hits3), mean)

valid <- !is.na(sig1) & !is.na(sig3) & sig1 > 0 & sig3 > 0
ratio_me1_me3 <- sig1[valid] / sig3[valid]
n_valid <- sum(valid)

message(sprintf("Anchors with valid signal for both marks: %d", n_valid))
if (n_valid == 0) stop("No dual anchors had measurable signal for both marks.")

# --- Diagnostics ---
cat("\n=== H3K4me1 / H3K4me3 (dual-positive anchors only) ===\n")
print(round(quantile(ratio_me1_me3, probs = c(0, 0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99, 1)), 3))

cat(sprintf("\n--- Threshold diagnostics ---\n"))
cat(sprintf("Ratio >= 1.0:  %.1f%%  (H3K4me1 >= H3K4me3)\n", 100 * mean(ratio_me1_me3 >= 1.0)))
cat(sprintf("Ratio >= 1.5:  %.1f%%\n", 100 * mean(ratio_me1_me3 >= 1.5)))
cat(sprintf("Ratio >= 2.0:  %.1f%%  (threshold 2; these stay dual)\n", 100 * mean(ratio_me1_me3 >= 2.0)))
cat(sprintf("Ratio >= 3.0:  %.1f%%  (threshold 3; more stringent)\n", 100 * mean(ratio_me1_me3 >= 3.0)))
cat(sprintf("Ratio <  2.0:  %.1f%%  (reclassified to P)\n", 100 * mean(ratio_me1_me3 < 2.0)))

# --- Plot ---
if (!is.null(out_plot)) {
  df <- data.frame(log2_ratio = log2(ratio_me1_me3 + 1e-6))
  p <- ggplot(df, aes(x = log2_ratio)) +
    geom_histogram(aes(y = after_stat(density)), bins = 80,
                   fill = "#5D6D7E", alpha = 0.7, boundary = 0) +
    geom_density(color = "#E41A1C", linewidth = 1) +
    geom_vline(xintercept = log2(c(1, 2, 3)), linetype = "dashed",
               color = c("#999999", "#27AE60", "#E41A1C"), linewidth = 0.8) +
    annotate("text", x = log2(c(1, 2, 3)), y = Inf, vjust = 2,
             label = c("1x", "2x", "3x"),
             color = c("#999999", "#27AE60", "#E41A1C"), size = 3.5) +
    scale_x_continuous(breaks = seq(-5, 5, 1),
                       labels = c("1/32","1/16","1/8","1/4","1/2","1",
                                  "2","4","8","16","32")) +
    labs(x = expression(H3K4me1 / H3K4me3~(log[2])),
         y = "Density",
         title = sprintf("H3K4me1 / H3K4me3 Signal Ratio\n(%d dual-positive anchors)", n_valid),
         subtitle = "Ratio >= 2 (green) → dual ; Ratio < 2 → P") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  ggsave(out_plot, p, width = 8, height = 5)
  message("\nPlot saved to ", out_plot)
}

message("\nDone.")
