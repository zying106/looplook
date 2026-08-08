#' looplook: A Triple-Layer Orthogonal Annotation Framework for Chromatin Interactions
#'
#' looplook classifies loop anchors -- the fundamental units connecting 3D
#' chromatin contacts to downstream peak and feature annotation -- through a
#' \strong{three-layer orthogonal annotation framework} in which each layer
#' integrates an independent category of experimental evidence:
#' \itemize{
#'   \item \strong{Genomic annotation} provides the baseline anchor assignment
#'     from genomic coordinates (TSS and gene-body overlap).
#'   \item \strong{Expression-aware refinement} adds transcriptional-activity
#'     context from transcriptomic data (CAGE-seq, RNA-seq, TT-seq, SLAM-seq).
#'   \item \strong{Chromatin-Aware Refinement} supplies orthogonal
#'     histone-mark validation (ChIP-seq, CUT&Tag, ATAC-seq).
#' }
#' These anchor-level classifications propagate through the loop network to
#' determine target-gene assignments for every auxiliary feature (ChIP-seq
#' peaks, ATAC-seq regions, GWAS variants). Because the layers function
#' independently, users may deploy any subset matched to their available data
#' modalities.
#'
#' @section Core modules:
#' \itemize{
#'   \item \strong{Consolidation & Consensus:} Graph-based clustering to harmonize
#'     replicates and build high-confidence 3D interactomes
#'     (\code{\link{consolidate_chromatin_loops}}).
#'   \item \strong{3D-Guided Annotation:} Hierarchical peak-to-gene mapping with
#'     expression-aware conflict resolution
#'     (\code{\link{annotate_peaks_and_loops}}).
#'   \item \strong{Expression-Aware Refinement:} Transcriptome-aware filtering
#'     and topological reclassification of silent regulatory elements
#'     (\code{\link{refine_loop_anchors_by_expression}}).
#'   \item \strong{Chromatin-Aware Refinement:} Chromatin mark-based anchor
#'     reclassification, dual-signature element detection, and orthogonal
#'     validation of expression-derived hypotheses
#'     (\code{\link{refine_loop_anchors_by_chromatin}},
#'     \code{\link{validate_epeG_by_chromatin}}).
#'   \item \strong{Automated Profiling:} End-to-end multi-omics analysis including
#'     GSEA, GO enrichment, motif scanning, and PPI networks
#'     (\code{\link{profile_target_genes}}).
#'   \item \strong{Visualization:} IGV-style multi-track plots and flower plots
#'     (\code{\link{plot_peaks_interactions}},
#'     \code{\link{draw_flower_simplified}}).
#' }
#'
#' @section Data I/O:
#' \itemize{
#'   \item \strong{BEDPE:} Read and convert chromatin interaction data
#'     (\code{\link{bedpe_to_gi}}).
#'   \item \strong{BED:} Read simple genomic region files
#'     (\code{\link{read_simple_bed}}).
#'   \item \strong{Spatial Clustering:} Merge proximal chromatin loops
#'     (\code{\link{reduce_ginteractions}}).
#' }
#'
#' @author
#' \strong{Maintainer}: Ying ZHANG \email{12207129@zju.edu.cn} (ORCID: 0009-0005-9644-7062)
#'
#' Contributors: Xingze HUANG \email{22407026@zju.edu.cn} (ORCID: 0009-0002-9286-1344); Ye CHEN \email{chenyephd@zju.edu.cn}
#'
#' Funding: Liang XU \email{xuliang.phd@zju.edu.cn}
#'
#' @seealso
#' \itemize{
#'   \item \url{https://github.com/zying106/looplook}
#'   \item \url{https://zying106.github.io/looplook/}
#'   \item Report bugs at \url{https://github.com/zying106/looplook/issues}
#' }
#'
#' @return \code{NULL}
#' @aliases looplook
"_PACKAGE"
