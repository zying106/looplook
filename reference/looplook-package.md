# looplook: Integrative Suite for Target Assignment and Functional Annotation of Chromatin Interactions

looplook is a versatile R package for target assignment and functional
annotation of chromatin interactions. It leverages customizable genomic
feature integration, expression-aware refinement, and connected
components clustering to prioritize biologically relevant chromatin
interactions. The package also provides flexible tools for multi-omics
data integration, functional annotation, and data visualization.

## Value

`NULL`

## Core modules

- **Consolidation & Consensus:** Graph-based clustering to harmonize
  replicates and build high-confidence 3D interactomes
  ([`consolidate_chromatin_loops`](https://zying106.github.io/looplook/reference/consolidate_chromatin_loops.md)).

- **3D-Guided Annotation:** Hierarchical peak-to-gene mapping with
  expression-aware conflict resolution
  ([`annotate_peaks_and_loops`](https://zying106.github.io/looplook/reference/annotate_peaks_and_loops.md)).

- **Expression-Aware Refinement:** Transcriptome-guided filtering and
  topological reclassification of silent regulatory elements
  ([`refine_loop_anchors_by_expression`](https://zying106.github.io/looplook/reference/refine_loop_anchors_by_expression.md)).

- **Automated Profiling:** End-to-end multi-omics analysis including
  GSEA, GO enrichment, motif scanning, and PPI networks
  ([`profile_target_genes`](https://zying106.github.io/looplook/reference/profile_target_genes.md)).

- **Visualization:** IGV-style multi-track plots, flower plots, and
  UpSet intersection diagrams
  ([`plot_peaks_interactions`](https://zying106.github.io/looplook/reference/plot_peaks_interactions.md),
  [`draw_flower_simplified`](https://zying106.github.io/looplook/reference/draw_flower_simplified.md),
  [`draw_upset_intersections`](https://zying106.github.io/looplook/reference/draw_upset_intersections.md)).

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
