# looplook: A Triple-Layer Orthogonal Annotation Framework for Chromatin Interactions

looplook classifies loop anchors – the fundamental units connecting 3D
chromatin contacts to downstream peak and feature annotation – through a
**three-layer orthogonal annotation framework** in which each layer
integrates an independent category of experimental evidence:

- **Genomic annotation** provides the baseline anchor assignment from
  genomic coordinates (TSS and gene-body overlap).

- **Expression-aware refinement** adds transcriptional-activity context
  from transcriptomic data (CAGE-seq, RNA-seq, TT-seq, SLAM-seq).

- **Chromatin-Aware Refinement** supplies orthogonal histone-mark
  validation (ChIP-seq, CUT&Tag, ATAC-seq).

These anchor-level classifications propagate through the loop network to
determine target-gene assignments for every auxiliary feature (ChIP-seq
peaks, ATAC-seq regions, GWAS variants). Because the layers function
independently, users may deploy any subset matched to their available
data modalities.

## Value

`NULL`

## Core modules

- **Consolidation & Consensus:** Graph-based clustering to harmonize
  replicates and build high-confidence 3D interactomes
  ([`consolidate_chromatin_loops`](https://zying106.github.io/looplook/reference/consolidate_chromatin_loops.md)).

- **3D-Guided Annotation:** Hierarchical peak-to-gene mapping with
  expression-aware conflict resolution
  ([`annotate_peaks_and_loops`](https://zying106.github.io/looplook/reference/annotate_peaks_and_loops.md)).

- **Expression-Aware Refinement:** Transcriptome-aware filtering and
  topological reclassification of silent regulatory elements
  ([`refine_loop_anchors_by_expression`](https://zying106.github.io/looplook/reference/refine_loop_anchors_by_expression.md)).

- **Chromatin-Aware Refinement:** Chromatin mark-based anchor
  reclassification, dual-signature element detection, and orthogonal
  validation of expression-derived hypotheses
  ([`refine_loop_anchors_by_chromatin`](https://zying106.github.io/looplook/reference/refine_loop_anchors_by_chromatin.md),
  [`validate_epeG_by_chromatin`](https://zying106.github.io/looplook/reference/validate_epeG_by_chromatin.md)).

- **Automated Profiling:** End-to-end multi-omics analysis including
  GSEA, GO enrichment, motif scanning, and PPI networks
  ([`profile_target_genes`](https://zying106.github.io/looplook/reference/profile_target_genes.md)).

- **Visualization:** IGV-style multi-track plots and flower plots
  ([`plot_peaks_interactions`](https://zying106.github.io/looplook/reference/plot_peaks_interactions.md),
  [`draw_flower_simplified`](https://zying106.github.io/looplook/reference/draw_flower_simplified.md)).

## Data I/O

- **BEDPE:** Read and convert chromatin interaction data
  ([`bedpe_to_gi`](https://zying106.github.io/looplook/reference/bedpe_to_gi.md)).

- **BED:** Read simple genomic region files
  ([`read_simple_bed`](https://zying106.github.io/looplook/reference/read_simple_bed.md)).

- **Spatial Clustering:** Merge proximal chromatin loops
  ([`reduce_ginteractions`](https://zying106.github.io/looplook/reference/reduce_ginteractions.md)).

## See also

- <https://github.com/zying106/looplook>

- <https://zying106.github.io/looplook/>

- Report bugs at <https://github.com/zying106/looplook/issues>

## Author

**Maintainer**: Ying ZHANG <12207129@zju.edu.cn> (ORCID:
0009-0005-9644-7062)

Contributors: Xingze HUANG <22407026@zju.edu.cn> (ORCID:
0009-0002-9286-1344); Ye CHEN <chenyephd@zju.edu.cn>

Funding: Liang XU <xuliang.phd@zju.edu.cn>
