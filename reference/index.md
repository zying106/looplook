# Package index

## Module 1: Data Consolidation & Preprocessing

Foundational data-cleaning engine for merging and standardizing 3D
interactomes.

- [`consolidate_chromatin_loops()`](https://zying106.github.io/looplook/reference/consolidate_chromatin_loops.md)
  : Consolidate and Integrate Chromatin Loops from Replicates or
  Multiple Sources
- [`filter_chromatin_loops()`](https://zying106.github.io/looplook/reference/filter_chromatin_loops.md)
  : Filter Chromatin Loops by Blacklist and/or Region of Interest
- [`bedpe_to_gi()`](https://zying106.github.io/looplook/reference/bedpe_to_gi.md)
  : Read BEDPE File into a GInteractions Object
- [`reduce_ginteractions()`](https://zying106.github.io/looplook/reference/reduce_ginteractions.md)
  : Spatial Clustering of GInteractions

## Module 2: 3D-Guided Peak Annotation

The core mapping engine for resolving 1D-to-3D spatial target
assignments.

- [`annotate_peaks_and_loops()`](https://zying106.github.io/looplook/reference/annotate_peaks_and_loops.md)
  : Spatial mapping of 1D genomic features to 3D chromatin interaction
  targets
- [`resolve_gene_conflicts()`](https://zying106.github.io/looplook/reference/resolve_gene_conflicts.md)
  : Resolve Gene Conflicts via Biotype Priority Then Expression

## Module 3: Expression-Aware & Chromatin-Aware Refinement

Advanced logic for eliminating transcriptionally silent contacts and
validating regulatory elements.

- [`refine_loop_anchors_by_expression()`](https://zying106.github.io/looplook/reference/refine_loop_anchors_by_expression.md)
  : Expression-Aware refinement of loop anchors and target linkages
- [`refine_loop_anchors_by_chromatin()`](https://zying106.github.io/looplook/reference/refine_loop_anchors_by_chromatin.md)
  : Chromatin-Aware refinement of loop anchor classification
- [`validate_epeG_by_chromatin()`](https://zying106.github.io/looplook/reference/validate_epeG_by_chromatin.md)
  : Orthogonal validation of eP/eG reclassification by chromatin marks

## Module 4: Automated Functional Profiling

End-to-end integration with JASPAR motifs, GO enrichment, and PPI
networks.

- [`profile_target_genes()`](https://zying106.github.io/looplook/reference/profile_target_genes.md)
  : Integrative functional annotation and profiling of target genes

## Module 5: IGV-Style & Statistical Visualization

Publication-ready tools for rendering genomic tracks and statistical
summaries.

- [`plot_peaks_interactions()`](https://zying106.github.io/looplook/reference/plot_peaks_interactions.md)
  : Integrative visualization of 3D chromatin loops and genomic features
- [`draw_flower_simplified()`](https://zying106.github.io/looplook/reference/draw_flower_simplified.md)
  : Draw Simplified Flower Plot for Core vs. Unique Genes

## Package Overview & Report

Package-level documentation and one-click HTML report generation.

- [`looplook-package`](https://zying106.github.io/looplook/reference/looplook-package.md)
  [`looplook`](https://zying106.github.io/looplook/reference/looplook-package.md)
  : looplook: A Triple-Layer Orthogonal Annotation Framework for
  Chromatin Interactions
- [`looplook_report()`](https://zying106.github.io/looplook/reference/looplook_report.md)
  : Render a Publication-Ready 3D Annotation Report

## Auxiliary Utilities

Internal helper functions for file handling and data parsing.

- [`read_simple_bed()`](https://zying106.github.io/looplook/reference/read_simple_bed.md)
  : Read a Simple BED File into a GRanges Object
- [`print(`*`<looplook_karyo>`*`)`](https://zying106.github.io/looplook/reference/print.looplook_karyo.md)
  : Print looplook karyogram
