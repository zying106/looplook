# Package index

## Module 1: Data Consolidation & Preprocessing

Foundational data-cleaning engine for merging and standardizing 3D
interactomes.

- [`consolidate_chromatin_loops()`](https://zying106.github.io/looplook/reference/consolidate_chromatin_loops.md)
  : Consolidate and Integrate Chromatin Loops from Multiple Sources
- [`bedpe_to_gi()`](https://zying106.github.io/looplook/reference/bedpe_to_gi.md)
  : Read BEDPE File into a GInteractions Object
- [`reduce_ginteractions()`](https://zying106.github.io/looplook/reference/reduce_ginteractions.md)
  : Spatial Clustering of GInteractions

## Module 2: 3D-Guided Peak Annotation

The core mapping engine for resolving 1D-to-3D spatial target
assignments.

- [`annotate_peaks_and_loops()`](https://zying106.github.io/looplook/reference/annotate_peaks_and_loops.md)
  : Advanced 3D Chromatin Loop Annotation & Multi-Omics Feature Mapping
  Engine

## Module 3: Expression-Aware Refinement

Advanced logic for eliminating transcriptionally silent contacts using
RNA-seq data.

- [`refine_loop_anchors_by_expression()`](https://zying106.github.io/looplook/reference/refine_loop_anchors_by_expression.md)
  : Transcriptome-Guided 3D Network Refinement & Anchor Reclassification

## Module 4: Automated Functional Profiling

End-to-end integration with JASPAR motifs, GO enrichment, and PPI
networks.

- [`profile_target_genes()`](https://zying106.github.io/looplook/reference/profile_target_genes.md)
  : Comprehensive Functional Profiling of Loop-Associated Target Genes

## Module 5: IGV-Style & Statistical Visualization

Publication-ready tools for rendering genomic tracks and statistical
summaries.

- [`plot_peaks_interactions()`](https://zying106.github.io/looplook/reference/plot_peaks_interactions.md)
  : Plot HiChIP/ChIA-PET Loops in IGV-Style Integrative View
- [`draw_flower_simplified()`](https://zying106.github.io/looplook/reference/draw_flower_simplified.md)
  : Draw Simplified Flower Plot for Core vs. Unique Genes
- [`draw_upset_intersections()`](https://zying106.github.io/looplook/reference/draw_upset_intersections.md)
  : Generate UpSet Plot for Gene Set Intersections

## Auxiliary Utilities

Internal helper functions for file handling and data parsing.

- [`read_simple_bed()`](https://zying106.github.io/looplook/reference/read_simple_bed.md)
  : Read a Simple BED File into a GRanges Object
