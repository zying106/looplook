#' looplook: Integrative Suite for Target Assignment and Functional Annotation of Chromatin Interactions
#'
#' looplook is a versatile R package for target assignment and functional
#' annotation of chromatin interactions. It leverages customizable genomic
#' feature integration, expression-aware refinement, and connected components
#' clustering to prioritize biologically relevant chromatin interactions. The
#' package also provides flexible tools for multi-omics data integration,
#' functional annotation, and data visualization.
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
#'     reclassification, dual-function element detection, and orthogonal
#'     validation of expression-derived hypotheses
#'     (\code{\link{refine_loop_anchors_by_chromatin}},
#'     \code{\link{validate_epeG_by_chromatin}}).
#'   \item \strong{Automated Profiling:} End-to-end multi-omics analysis including
#'     GSEA, GO enrichment, motif scanning, and PPI networks
#'     (\code{\link{profile_target_genes}}).
#'   \item \strong{Visualization:} IGV-style multi-track plots, flower plots, and
#'     UpSet intersection diagrams (\code{\link{plot_peaks_interactions}},
#'     \code{\link{draw_flower_simplified}}, \code{\link{draw_upset_intersections}}).
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
