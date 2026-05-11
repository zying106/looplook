# Looplook: An integrative suite for target assignment and functional annotation of chromatin interactions

------------------------------------------------------------------------

## Introduction

Welcome to **`looplook`**, a versatile R/Bioconductor toolkit developed
to integrate three-dimensional (3D) chromatin architecture data (e.g.,
HiChIP, ChIA-PET, Hi-C) with other tabular omics datasets, including
transcriptomics, chromatin accessibility, protein-DNA interactions
(derived from ChIP-seq, CUT&Tag, or CUT&RUN), and genetic variants
annotated by genome-wide association studies (GWAS).

Numerous studies have demonstrated that distal regulatory elements
physically interact with target gene promoters via **3D chromatin
looping**, thereby regulating the expression of genes located tens of
kilobases to megabases away along the linear genome. However,
conventional annotation strategies typically assign putative regulatory
elements (e.g., informed by ChIP-seq peaks) to their **nearest
neighboring genes in cis**, which often deviates from actual biological
mechanisms. Consequently, the accurate assignment of non-coding genetic
variants or orphan regulatory peaks to their **cognate target genes**
represents a major bottleneck in the functional annotation of genomic
regulatory elements.

To overcome this challenge, **`looplook`** systematically prioritizes
experimentally validated spatial chromatin contacts to batch-annotate
thousands of regulatory elements at a **genome-wide, high-throughput
scale**, thereby identifying their candidate target genes with high
confidence and unprecedented efficiency.

Beyond its utility as an integrative tool for target gene annotation,
**`looplook`** serves as a **standalone platform** dedicated to
chromatin loop analysis. Even in the absence of auxiliary omics data,
this tool enables systematic annotation of the 3D chromatin interactome
itself, classification of complex spatial interaction topologies (e.g.,
Enhancer-Promoter, Promoter-Promoter interactions), and quantification
of node connectivity, which facilitates the identification of **dense
regulatory hubs** and **enhancer cliques** (e.g., super-enhancers) that
drive cell-type-specific transcriptional programs.

### Key Features & Capabilities

1.  **The “3D Spatial Bridge” for Multi-Omics**: Integrates auxiliary
    genomic features (e.g., GWAS risk SNPs, ChIP-seq peaks) with 3D
    chromatin loops. It implements a rigorous “Smart Fallback” logic
    that prioritizes distal loop-mediated target genes while reverting
    to the nearest neighboring genes when no spatial interactions are
    detected, ensuring comprehensive and gapless target annotation.

2.  **Comprehensive Loop Annotation & Topological Hub Detection**: In
    addition to loop mapping, **`looplook`** enables biologically
    informed annotation of the 3D chromatin interactome. It helps to
    categorize spatial interaction types and computes spatial node
    degrees to identify candidate 3D regulatory hubs.

3.  **Expression-Aware Refinement**: Incorporates quantitative
    expressional data to filter out spurious structural contacts. As
    some chromatin regions may exhibit both enhancer and promoter
    activities, **`looplook`** provides a function to reclassify loop
    anchors which associate with transcriptionally silent reference
    promoters into **enhancer-like elements (eP)**, yielding a refined
    network of regulatory interactions.

4.  **Reproducible Consolidation & Multi-source Consensus**: Utilizes
    graph-theoretic clustering to effectively harmonize biological
    replicates and mitigate potential technical noise. The framework can
    be used to identify multi-source consensus by integrating datasets
    from various experimental designs, such as HiChIP assays targeting
    different factors (e.g., **H3K27ac and Pol II**). This enables
    **orthogonal validation** to facilitate the construction of a
    high-confidence, unified 3D chromatin interactome based on
    consistent spatial features.

5.  **Automated Multi-Omics Functional Interpretations**: Provides an
    end-to-end analytical engine for generating publication-grade
    visualizations and functional interpretations. This includes
    assessments of network connectivity and expression dynamics,
    transcription factor motif enrichment (**SeqLogos**), and pathway
    association networks (**Cnetplots**).

### Comparison with Existing Tools

`looplook` bridges the gap between conventional 1D linear annotation
strategies and specialized 3D chromatin interaction tools, providing a
unified environment for target assignment, expression-aware refinement,
consensus consolidation, and functional inference. Table
@ref(tab:comparison-table) highlights the core capabilities of
`looplook` relative to existing tools in the ecosystem.

| Feature | ChIPseeker | GREAT | GenomicInteractions | FUMA | ABC Model | looplook |
|:---|:--:|:--:|:--:|:--:|:--:|:--:|
| Primary annotation scheme | 1D | 1D | 3D | 3D | 3D | 1D & 3D |
| User-input 3D interactome | ✘ | ✘ | ✔ | ✘ | ✔ | ✔ |
| User-input omics data | ✔ | ✔ | ✘ | ✘ | ✔ | ✔ |
| Graph-based topological consensus | ✘ | ✘ | ✘ | ✘ | ✘ | ✔ |
| Expression-guided refinement | ✘ | ✘ | ✘ | Partial | ✘ | ✔ |
| Topological reclassification | ✘ | ✘ | ✘ | ✘ | ✘ | ✔ |
| Multi-hop network diffusion | ✘ | ✘ | ✘ | ✘ | ✘ | ✔ |
| Smart linear fallback | ✘ | ✘ | ✘ | ✘ | ✘ | ✔ |
| Functional inference | ✘ | ✔ | ✘ | ✘ | ✘ | ✔ |
| Multi-track visualization | ✘ | ✘ | ✘ | ✘ | ✘ | ✔ |
| Local execution | ✔ | ✘ | ✔ | ✘ | ✔ | ✔ |
| Ecosystem | R | Web | R | Web | Python | R |

Comparison of core features between looplook and existing genomic
annotation tools. {.table}

*Notes:* “1D” = linear proximity-based mapping; “3D” = mapping via
spatial chromatin topology; “✔” = native support; “✘” = not supported;
“Partial” = partial support (FUMA relies on precompiled data rather than
user-provided input).

`looplook` is the only tool in this comparison that unifies both 1D and
3D annotation paradigms within a single R-based workflow, while
additionally offering graph-theoretic consensus consolidation,
topological reclassification of silent regulatory elements, and local
execution for data privacy — capabilities not natively available in any
of the alternatives listed above.

------------------------------------------------------------------------

## Installation

To ensure full functionality, install `looplook` along with its
recommended Bioconductor annotation dependencies:

``` r

if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# Install annotation databases required for human (hg38) analysis
BiocManager::install(c("TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db"))

# Install looplook from Bioconductor (once accepted)
# BiocManager::install("looplook")

# Or install the development version from GitHub
devtools::install_github("zying106/looplook")
```

``` r

library(looplook)
# Create a temporary directory to store all output plots and Excel files
out_dir <- tempdir()
```

------------------------------------------------------------------------

## Module 1: Data Consolidation & Preprocessing

In 3D genomics analyses, individual replicates typically exhibit certain
degrees of **inconsistency or noise**. The
**`consolidate_chromatin_loops`** function serves as the foundational
data-cleaning module, merging multiple replicates into a standardized,
unified 3D chromatin interaction coordinate framework.

This module identifies **multi-source consensus** by integrating
datasets from varied experimental designs, such as HiChIP assays
targeting H3K27ac and Pol II. By capturing shared spatial interactions
across these inputs, the framework enables **orthogonal validation** to
**mitigate potential technical artifacts** and assay-specific biases.
This approach leverages consistent spatial features to provide a robust
foundation for downstream analyses supported by multiple lines of
experimental evidence.

#### Parameter Strategy

- **`mode`**: Defines the overarching merging algorithm.
  - `"consensus"` *(Recommended)*: Implements graph-based connected
    component analysis to cluster spatially proximal anchors across
    samples. Only retains clusters detected in \>= `min_consensus`
    biological replicates.
  - `"intersect"`: Enforces strict reference-based filtering, retaining
    loops that show full genomic overlap with the reference file (File
    1).
  - `"union"`: Retains all chromatin interactions across the entire
    cohort, ideal for exploratory pan-tissue analyses.
- **`min_consensus`**: When using `"consensus"` mode, this parameter
  defines the minimum number of biological replicates in which a loop
  cluster must be detected to be retained in the final dataset. If set
  to `NULL`, the algorithm dynamically computes a strict majority
  threshold (e.g., \>= 75% of replicates).
- **`gap`**: Defines the maximum allowable spatial distance (in base
  pairs) between loop anchors for them to be considered as part of the
  same physical cluster (default: `1000`).
- **`min_raw_score` vs. `min_score` (The Dual-Filter)**:
  - `min_raw_score` acts as a pre-filter applied to individual BEDPE
    files before clustering (e.g., removing singleton noise where PET
    count \< 2) to substantially reduce computational memory overhead.
  - `min_score` serves as a post-filter applied to the final merged
    chromatin interactome to ensure high-confidence interactions. In
    `"consensus"` and `"union"` modes, this representative score is
    computed in a replicate-balanced manner: loop scores are averaged
    within each replicate first, then averaged across replicates.
- **`blacklist_species`**: Automatically excludes chromatin loops
  overlapping with high-variance, artifact-prone genomic regions (e.g.,
  centromeres, telomeres) by integrating the official ENCODE blacklist
  for specified species (e.g., `"hg38"`, `"mm10"`).
- **`region_of_interest`**: Accepts an auxiliary BED file (e.g., a
  specific disease-associated locus or ChIP-seq peak set) to exclude
  global background interactions, outputting only loops with physical
  connectivity to the target genomic region.

### Example 1: Building a Global Consensus Interactome

This represents the standard workflow for generating a high-confidence,
whole-genome 3D chromatin interaction dataset via integration of
biological replicates and removal of blacklist-associated artifacts.

``` r

# Locate example BEDPE replicates
f1 <- system.file("extdata", "example_loops_1.bedpe", package = "looplook")
f2 <- system.file("extdata", "example_loops_2.bedpe", package = "looplook")
global_out <- file.path(out_dir, "consensus_loops_global.bedpe")

# Execute consensus merging with strict Quality Control
consensus_global <- consolidate_chromatin_loops(
  files = c(f1, f2),
  mode = "consensus",
  gap = 1000,
  min_raw_score = 2, # Pre-filter sequencing noise
  blacklist_species = "hg38", # Apply artifact blacklist
  out_file = global_out
)
```

### Example 2: Feature-Driven Loop Filtering

In most scenarios, researchers are interested in 3D chromatin
interactions associated with specific biological features, such as
active enhancers (marked by H3K27ac) or particular transcription
factors.

By providing a BED file to `region_of_interest`, the tool enables
efficient “on-the-fly” capture, retaining only those loops whose anchors
overlap with the specified functional regions.

``` r

# Locate a BED file containing H3K27ac ChIP-seq peaks
h3k27ac_peaks <- system.file("extdata", "example_k27ac_peaks.bed", package = "looplook")
targeted_out <- file.path(out_dir, "H3K27ac_anchored_loops.bedpe")

# Execute consensus merging, strictly capturing H3K27ac-associated loops
consensus_targeted <- consolidate_chromatin_loops(
  files = c(f1, f2),
  mode = "consensus",
  gap = 1000,
  min_raw_score = 2,
  region_of_interest = h3k27ac_peaks, # Surgical capture via 1D features
  out_file = targeted_out
)
```

### Example 3: Multi-source Consensus for Orthogonal Validation

In integrative multi-omics analyses, orthogonal validation across
independent functional assays is critical for identifying robust 3D
regulatory interactions. For instance, enhancer-promoter (E-P) loops are
typically supported by the spatial co-localization of H3K27ac (marking
active chromatin) and RNA polymerase II (indicative of transcriptional
machinery engagement). By utilizing the `"consensus"` mode to intersect
independent loop datasets, `looplook` facilitates the identification of
chromatin interactions supported by multi-source evidence. This rigorous
intersection strategy effectively mitigates assay-specific technical
biases and artifacts, yielding a refined set of topological hubs that
couple physical chromatin architecture with active transcriptional
regulation.

``` r

loop_k27ac <- system.file("extdata", "example_loops_H3K27ac.bedpe", package = "looplook")
loop_pol2 <- system.file("extdata", "example_loops_pol2.bedpe", package = "looplook")
dual_functional_out <- file.path(out_dir, "Dual_Functional_Consensus.bedpe")

consensus_dual <- consolidate_chromatin_loops(
  files = c(loop_k27ac, loop_pol2),
  mode = "consensus",
  gap = 1000,
  min_raw_score = 2,
  out_file = dual_functional_out
)
```

------------------------------------------------------------------------

#### Module 2: 3D-Guided Annotation & Mapping

This module constitutes the core mapping engine of **`looplook`**. Using
**`annotate_peaks_and_loops`**, users can classify the topological
architecture of the interactome. Optionally, users may supply a
`target_bed` (e.g., GWAS loci, eQTLs, ChIP-seq peaks) to trace
non-coding regulatory signals to their putative functional target genes.

To resolve mapping ambiguities within densely populated gene loci, where
a single loop anchor overlaps multiple candidate genes, the engine
implements a rigorous three-step hierarchical annotation pipeline:

1.  **Expression Pre-filter**: Eliminates transcriptionally silent genes
    based on the user-provided expression matrix.
2.  **Functional Biotype Prioritization**: Prioritizes the remaining
    candidates by functional class in the following order: Protein
    Coding \> lncRNA \> Pseudogene.
3.  **Dominant Expression Tiebreaker**: Designates the gene with the
    highest transcriptional abundance as the target gene to further
    resolve any remaining mapping ambiguities.

#### Parameter Strategy & Core Inputs

This multi-omic integration relies on several key parameters:

- **`bedpe_file`**: The 3D structural framework (e.g., the
  high-confidence consensus interactome generated in Module 1).
- **`target_bed`**: The auxiliary genomic features of interest (e.g., a
  BED file of GWAS SNPs, ATAC-seq peaks, or ChIP-seq peaks) awaiting
  spatial target assignment.
- **`expr_matrix_file` & `sample_columns`**: Optional but strongly
  recommended. An expression matrix (e.g., RNA-seq) enables the
  Expression Pre-filter and Tiebreaker logic, substantially reducing
  false-positive gene assignments.
- **`species`**: Specifies the genome assembly (e.g., `"hg38"`,
  `"mm10"`). The engine automatically loads the corresponding
  Bioconductor `TxDb` and `org.db` annotations.
- **`neighbor_hop`**: An advanced topological parameter for network
  traversal:
  - `0` *(Default)*: Restricts annotation to direct physical contacts
    (Anchor A ↔︎ Anchor B).
  - `1` *(Hub Mode)*: Evaluates secondary network effects within
    enhancer cliques (If Anchor A ↔︎ Anchor B ↔︎ Anchor C, this function
    enables topological linkage between A and C).
- **`tss_region`**: Defines the genomic window surrounding the
  transcription start site (TSS) used to define promoter regions
  (default: `c(-3000, 3000)` bp).

### Example A: Integrative Analysis (Loops + Genomic Features + RNA-seq)

This represents the comprehensive “3D Spatial Bridge” scenario. It links
external genomic feature (peaks) to 3D‑inferred target genes, while
simultaneously resolving ambiguous assignments in dense genomic regions
by weighing gene expression.

``` r

# Locate auxiliary genomic features and RNA-seq expression matrix
expr_path <- system.file("extdata", "example_tpm.txt", package = "looplook")
atac_path <- system.file("extdata", "example_peaks.bed", package = "looplook")
if (requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE) &&
  requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  res_integrated <- annotate_peaks_and_loops(
    bedpe_file = global_out, # Uses the global consensus network from Module 1
    target_bed = atac_path, # External genomic features
    expr_matrix_file = expr_path, # Activates the rigorous 3-step tiebreaker
    sample_columns = c("con1", "con2"),
    species = "hg38",
    neighbor_hop = 0, # Focus on direct physical contacts
    hub_percentile = 0.95, # Top 5% nodes defined as hubs
    out_dir = out_dir,
    project_name = "Example_HiChIP_Integrative"
  )
}
```

### Example B: Intrinsic Loop Profiling (Loops ONLY)

For users focused exclusively on 3D interactome architecture without
auxiliary data related to genomic features, looplook functions as a
standalone platform for loop profiling and structural hub
identification.

``` r

if (requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE) &&
  requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  res_loops_only <- annotate_peaks_and_loops(
    bedpe_file = global_out,
    target_bed = NULL, # Omit genomic features entirely
    expr_matrix_file = expr_path, # Resolves multi-gene conflicts at loop anchors
    sample_columns = c("con1", "con2"),
    species = "hg38",
    neighbor_hop = 0,
    out_dir = out_dir,
    project_name = "Example_Loops_Only"
  )
}
```

#### Deep Dive: Output Data Dictionary

The function exports a comprehensive Excel workbook
(`*_Basic_Results.xlsx`) containing multi-layered annotations.
Interpretation of spatial regulatory logic depends critically on the
hierarchical structure of its columns:

##### 1. `target_annotation` (The Genomic Feature-to-3D Mapping Result)

*(Generated only when `target_bed` is provided)*

- **`annotation`**: Genomic context of the input peaks (e.g., “Distal
  Intergenic”, “Promoter”).
- **`Linked_Loop_IDs`**: Identifiers of 3D loops that overlap the peak.
- **`All_Loop_Connected_Genes`**: All genes (including promoters and
  gene bodies) topologically linked to the peak via 3D looping.
- **`Regulated_promoter_genes`**: A stringent subset restricted to genes
  whose promoters are directly contacted by the peak via 3D looping.
- **`Assigned_Target_Genes`**: Target genes assigned *exclusively* via
  3D loops. If the peak does not overlap any structural loop, this value
  remains empty (`NA`).
- **`*_Filled` Columns (The “Smart Fallback” Logic)**: Columns suffixed
  with `_Filled` (e.g., `Assigned_Target_Genes_Filled`,
  `Regulated_promoter_genes_Filled`) provide complete, gapless
  annotations. For peaks within looped regions, they retain 3D-refined
  distal targets; for peaks in unlooped regions, they will be annotated
  to the nearest linear gene.

##### 2. `loop_annotation` (3D Network Architecture)

- **`loop_type`**: Topological classification (e.g., E-P, P-P, E-E)
  determined by the biotype-aware logic.
- **`All_Anchor_Genes`**: All genes located linearly within the left and
  right anchors of a given loop.
- **`Putative_Target_Genes`**: Biologically resolved target genes. For
  an E-P loop, this field will report only the gene at the Promoter
  anchor.
- **`Cluster_All_Genes`**: All genes within the broader 3D cluster or
  clique in which the loop resides, reflecting the extended topological
  neighborhood.

##### 3. `anchor_loci_annotation` (Non-redundant Spatial Anchor Footprints)

- **`Cluster_Locus_Genes`**: All genes from loop anchors within a
  spatially interconnected component (defined by network-connected
  components), representing the full gene set associated with the
  deduplicated anchor loci.
- **Backward compatibility**: `anchor_annotation` is retained as an
  alias of `anchor_loci_annotation` for older scripts.

##### 4. `promoter_centric_stats` & `distal_element_stats` (Topological Hub Detection)

These metrics quantify interactome connectivity to identify candidate
regulatory hubs from complementary perspectives:

**A. Gene-Centric (`promoter_centric_stats`)**

- **`Total_Loops` / `n_Linked_Promoters` / `n_Linked_Distal`**:
  Deconstructs a gene’s 3D connectome. High `n_Linked_Distal` suggests
  potential regulation by enhancer clusters, whereas high
  `n_Linked_Promoters` may indicate possible participation in multi-gene
  transcription.
- **`Dominant_Interaction`**: The predominant topological class
  governing the gene (e.g., E-P vs. P-P).
- **`Is_High_Connectivity_Gene`** *(and Distal variant)*: Binary flags
  defined by the `hub_percentile` to prioritize putative 3D regulatory
  hubs.

**B. Enhancer-Centric (`distal_element_stats`)**

- **`Target_Genes`**: Promoters physically contacted by the distal
  anchor.
- **`n_Linked_Promoters`**: Facilitates identification of candidate
  master enhancers that contact multiple distinct genes.
- **`Is_High_Connectivity_Distal_Element`**: Flags highly connected
  non-coding elements that may act as 3D hubs.

#### Deep Dive: Output Visualizations

All plots are returned as objects in `res_integrated$plots$`, enabling
direct inspection and customisation with standard **ggplot2** layers
without re-running the pipeline.

##### 1. Global 3D Network Profiling

*(Always Generated)*

| Plot key | Description |
|----|----|
| `Basic_Donut` | Loop-type frequency (donut chart) |
| `Basic_Circular` | Unique target genes per loop type (circular bar) |
| `Basic_Flower` | Shared vs. unique genes across loop types |
| `Karyo_Anchors` | Genome-wide anchor density heatmap |
| `Karyo_LoopGenes` | Genome-wide loop-associated gene density |
| `Anchor_Genomic_Distribution` | Genomic context of anchors (pie) |
| `Target_Genomic_Distribution` | Genomic context of input features (pie) |

##### 2. Genomic Feature-to-3D Target Profiling

*(Requires `target_bed`)*

| Plot key | Description |
|----|----|
| `Target_Rose` | Loop types linked to input features (coxcomb) |
| `Target_Loop_Genomic_Distribution` | Genomic context of loop-connected features (pie) |
| `Karyo_TargetGenes` | Chromosomal density of assigned target genes |

##### Customising Plots

All ggplot-based plots accept standard **ggplot2** layers for
recolouring, theme adjustment, or annotation:

``` r

# Access a stored plot
p <- res_integrated$plots$Basic_Donut

# Extract data and original label mapping
donut_data <- p$data
labels_vec <- setNames(donut_data$legend_label, donut_data$loop_type)

# Apply a custom colour palette
p + ggplot2::scale_fill_manual(
  values = RColorBrewer::brewer.pal(length(labels_vec), "Dark2"),
  labels = labels_vec
) + ggplot2::theme(legend.position = "bottom")
```

To save any plot to disk, use
[`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html):

``` r

ggplot2::ggsave(
  filename = "my_custom_donut.pdf",
  plot = res_integrated$plots$Basic_Donut,
  width = 8, height = 6
)
```

For `looplook_karyo` objects (karyotype heatmaps), use
[`print()`](https://rdrr.io/r/base/print.html) to display:

``` r

print(res_integrated$plots$Karyo_Anchors) # renders via grid or opens in browser
```

------------------------------------------------------------------------

## Module 3: Expression-Aware Refinement

Physical proximity is a structural prerequisite, but not a direct proxy
for active transcriptional regulation. This module integrates
quantitative transcriptome data to systematically filter out
transcriptionally silent physical chromatin contacts, ensuring only
functionally relevant 3D interactions are retained.

#### Parameter Strategy & Core Inputs

The `refine_loop_anchors_by_expression` function implements functional
filtration via the following key parameters:

- **`annotation_res`**: Foundational 3D annotation output generated by
  Module 2.
- **`expr_matrix_file` & `sample_columns`**: Quantitative expression
  matrix (e.g., TPM, FPKM) and the specific replicates to average,
  establishing a robust baseline expression profile.
- **`threshold` & `unit_type`**: Defines the quantitative expression
  cutoff (e.g., `threshold = 1`, `unit_type = "TPM"`) required to
  consider a gene biologically active. This parameter enables downward
  compatibility with various normalization methods.
- **`reclassify_by_expression`**: A transformative logical parameter
  (`TRUE`/`FALSE`). When enabled, transcriptionally silent Promoters
  (`P`) and Gene Bodies (`G`) are not simply discarded; instead, they
  are reclassified as enhancer-like regulatory elements (`eP`, `eG`).
  This correction refines the regulatory topology (e.g., reclassifying a
  functionally silent `P-P` loop into a curated `eP-P` loop).

### Example A: Standard Filtration (Strict Removal)

In this scenario, a strict expression threshold eliminates loops where
putative target genes are transcriptionally silent, without modifying
the original structural classifications.

``` r

rdata_path <- system.file("extdata", "analysis_results.RData", package = "looplook")
load(rdata_path) # loads res_integrated into the environment

res_basic <- refine_loop_anchors_by_expression(
  annotation_res = res_integrated,
  expr_matrix_file = expr_path,
  sample_columns = c("con1", "con2"),
  threshold = 1.0,
  unit_type = "TPM",
  reclassify_by_expression = FALSE,
  out_dir = out_dir,
  project_name = "Example_Basic_Filter"
)
```

### Example B: Expression-Aware Reclassification (Recommended)

This advanced mode is strongly recommended, as it preserves the valuable
3D structural framework while refining functional annotations. A silent
promoter anchor is reclassified as an enhancer (`eP`), acknowledging
that while its reference gene is transcriptionally inactive, this region
itself may function as an active regulatory element for a distant gene.

``` r

refined_res <- refine_loop_anchors_by_expression(
  annotation_res = res_integrated,
  expr_matrix_file = expr_path,
  sample_columns = c("con1", "con2"),
  threshold = 1.0,
  unit_type = "TPM",
  reclassify_by_expression = TRUE,
  out_dir = out_dir,
  project_name = "Example_Reclassified_Filter"
)
```

#### Deep Dive: Filtration Visualizations

This module automatically generates a specialized suite of visual
diagnostics to quantify transcriptome-guided refinement efficacy and
characterize the surviving functional network:

##### 1. Global Filtration Profiling

*(Always Generated)*

- **Filtration Effect Dumbbell** (`*_Comparison_Dumbbell.pdf`):
  Quantifies the “cleaning power” of the expression threshold by
  visualizing the numerical contrast between original structural loops
  and surviving functional loops across all topological classes,
  highlighting the importance of removing untranscribed structural
  noise.

- **Refined Loop Proportion Rose** (`*_Rose.pdf`): Coxcomb chart
  illustrating the topological distribution (by count) of the
  interactome after expression-guided reclassification and filtration.

- **Refined Active Genes Karyotype** (`*_Refined_Karyo_Active.pdf`):
  Genome-wide ideograms mapping the chromosomal density of
  transcriptionally active loop-associated genes that passed the
  expression threshold.

##### 2. Genomic Feature-to-3D Target Filtration Profiling

*(Generated only if `target_bed` was provided in Module 2)*

- **Multi-Omics Sankey Diagram** (`Target_Sankey` htmlwidget in-memory;
  `target_sankey.html` alongside the report when rendered): Core
  integration diagnostic that maps the fate of inputted genomic features
  across three logical tiers: Genomic Region (L1) → 3D Loop Connection
  (L2) → Target Expression Status (L3). This intuitively reveals the
  proportion of GWAS/ChIP-seq peaks structurally connected to
  functionally active genes. When generated through
  [`looplook_report()`](https://zying106.github.io/looplook/reference/looplook_report.md),
  keep `target_sankey.html` in the same output directory as the main
  HTML report.

- **Refined Target Loop Donut** (`*_Target_Loop_Donut.pdf`): Visualizes
  the topological distribution of functionally active loops that bridge
  user-defined genomic features.

- **Refined Target Genes Karyotype**
  (`*_Refined_Karyo_TargetGenes.pdf`): Chromosomal density map of
  strictly refined putative target genes, facilitating visualization of
  transcriptionally verified multi-omics loci.

------------------------------------------------------------------------

### Module 4: Automated Functional Profiling

Following identification of high-confidence putative target genes,
`profile_target_genes` provides a fully automated, end-to-end
multi-omics analysis pipeline. It integrates 3D genomic interactions
with transcriptomic data to systematically characterize the functional
landscape and regulatory mechanisms underlying the identified target
genes.

#### Parameter Strategy: A Highly Modular Pipeline

This function supports highly flexible experimental design through a set
of routing parameters that define the gene set to be analyzed, as well
as toggle parameters that enable or disable specific downstream
biological modules.

##### 1. Core Data Inputs

- **`annotation_res`**: The foundational annotation result, inherited
  directly from Module 2 or the refined output from Module 3.
- **`diff_file` & `lfc_col`**: Differential expression summary (e.g.,
  from DESeq2) required for fold change profiling and gene set
  enrichment analysis (GSEA).
- **`expr_matrix_file` & `metadata_file`**: Normalized expression matrix
  and sample metadata, which are essential for generating heatmaps and
  Raincloud plots.

##### 2. Target Selection & Mapping Strategies

This section contains the most critical routing parameters that define
the biological scope and stringency of downstream analyses:

- **`target_source`** (The Biological Scope):
  - `"targets"` *(Variant/Peak-Centric)*: Focuses exclusively on the
    putative genes regulated by user-defined genomic features (e.g.,
    GWAS SNPs, ATAC-seq peaks).
  - `"loops"` *(Global Network-Centric)*: Analyzes the full 3D
    interactome, including all genes associated with spatial loops,
    without considering user-defined genomic features.
  - `c("loops", "targets")`: Analyzes both scopes simultaneously, which
    is strongly recommended for comparative functional profiling.
- **`target_mapping_mode`**:
  - `"all"` *(Default)*: Retains broad 3D regulatory targets, including
    distal enhancers looping to gene bodies.
  - `"promoter"`: Applies stringent filtering, requiring loops to anchor
    explicitly at canonical promoter regions.
- **`include_Filled`** (The Stringency Toggle): Logical parameter
  controlling the purity of spatial regulation:
  - `TRUE` *(Hybrid Mode)*: Utilizes the comprehensively merged
    annotation (`_Filled` columns from Module 2). It prioritizes 3D
    loop-derived targets while rescuing unlooped genomic features
    (peaks) by assigning them the nearest genes, ensuring no data loss
    and a complete functional overview.
  - `FALSE` *(Pure Spatial Mode)*: Highly stringent. It only keeps 3D
    loop-derived targets with clear annotation.
- **`use_nearest_gene`**: Logical parameter. If `TRUE`, the pipeline
  bypasses 3D loop-based assignment and uses only the nearest
  neighboring gene, serving as a conventional baseline reference for
  Control comparisons.

##### 3. Downstream Functional Analyses

- **`run_go`**: Executes Gene Ontology (Biological Process, BP)
  enrichment and generates Divergent Concept Networks.
- **`run_motif` & `genome_id`**: Scans proximal and distal anchor
  sequences against JASPAR core motifs (requires `BSgenome`).
- **`run_ppi` & `ppi_score`**: Constructs Protein-Protein Interaction
  (PPI) networks via the STRING database.

### Example A: Comprehensive Integrative Profiling (Recommended)

Ideal for complete multi-omics integration (e.g., HiChIP + ATAC-seq).
This mode comprehensively maps and profiles target genes of interest.

``` r

diff_path <- system.file("extdata", "example_deg.txt", package = "looplook")
meta_path <- system.file("extdata", "example_coldata.txt", package = "looplook")

res_A <- profile_target_genes(
  annotation_res = res_integrated,
  diff_file = diff_path,
  lfc_col = "log2FoldChange",
  expr_matrix_file = expr_path,
  metadata_file = meta_path,
  target_source = c("loops", "targets"),
  target_mapping_mode = "all",
  include_Filled = TRUE,
  use_nearest_gene = FALSE,
  project_name = "Scenario_A_Integrative",
  run_motif = FALSE,
  run_go = FALSE,
  run_ppi = FALSE
)
```

### Example B: Peak-Driven Strict Promoter Profiling (High Stringency)

This mode applies a rigorous dual-filtering strategy to yield a
high-confidence target gene set.

First, `target_source = "targets"` restricts profiling to genes
regulated specifically by user-defined genomic features (e.g., GWAS SNPs
or TF binding sites), excluding background spatial loops. Second,
`target_mapping_mode = "promoter"` enforces topological stringency:
loops linking peaks to target genes must directly contact canonical or
reference promoters (i.e., strict `E-P` or `P-P` loops), excluding more
permissive connections such as enhancer-gene body (`E-G`) interactions.

``` r

res_B <- profile_target_genes(
  annotation_res = res_integrated,
  diff_file = diff_path,
  lfc_col = "log2FoldChange",
  expr_matrix_file = expr_path,
  metadata_file = meta_path,
  target_source = "targets", # 1. Focus exclusively on inputted genomic features
  target_mapping_mode = "promoter", # 2. Require strict loop-to-promoter collision
  include_Filled = TRUE,
  use_nearest_gene = FALSE,
  project_name = "Scenario_B_StrictPromoter",
  run_motif = FALSE, # Set TRUE in real analysis to scan JASPAR motifs
  run_go = FALSE, # Set TRUE in real analysis for pathway enrichment
  run_ppi = FALSE
)
```

### Example C: 1D Linear Annotation Baseline (The Control)

By setting `use_nearest_gene = TRUE`, this mode intentionally bypasses
the 3D spatial topology and assigns features exclusively to their
nearest neighboring gene. It serves as a conventional “baseline
control.” Comparing results against Example A or B enables rigorous
demonstration of the novel functional insights gained from 3D chromatin
looping analyses relative to conventional method.

``` r

res_C <- profile_target_genes(
  annotation_res = res_integrated,
  diff_file = diff_path,
  lfc_col = "log2FoldChange",
  expr_matrix_file = expr_path,
  metadata_file = meta_path,
  target_source = "targets", # Focus on inputted peaks
  target_mapping_mode = "all",
  include_Filled = FALSE,
  use_nearest_gene = TRUE, # Nearest neighboring gene (The Control)
  project_name = "Scenario_C_LinearControl",
  run_motif = FALSE, # Set TRUE in real analysis to scan JASPAR motifs
  run_go = FALSE, # Set TRUE in real analysis for pathway enrichment
  run_ppi = FALSE
)
```

#### Deep Dive: Functional Visualizations

Results are returned as a named list, one element per `target_source`
(e.g., `"loops"`, `"targets"`). Each element contains `plots` (ggplot
objects), `target_gene_sets` (character vectors), and `go_results` (data
frames).

``` r

# Accessing results from a profiling run
names(res_profile) # one element per target_source ("loops", "targets")
names(res_profile$loops) # "go_results" "target_gene_sets" "plots"

# --- Target genes ---
res_profile$loops$target_gene_sets # named list of character vectors
names(res_profile$loops$target_gene_sets) # gene set keys (e.g. "All", "Up", "Down")

# --- GO enrichment table ---
go_df <- res_profile$loops$go_results # data frame
head(go_df[
  order(go_df$pvalue), # sort by significance
  c("Description", "ONTOLOGY", "pvalue", "Count", "geneID")
])

# --- All plot keys ---
names(res_profile$loops$plots)
```

##### Available Plot Keys

| Key                    | Description                                |
|------------------------|--------------------------------------------|
| `LFC_Violin`           | Differential expression violin + boxplot   |
| `Heatmap`              | Expression heatmap with sample annotations |
| `Scatter`              | Connectivity vs. expression scatter        |
| `Raincloud_LFC`        | Raincloud plot of log2 fold changes        |
| `Raincloud_Expr`       | Raincloud plot of expression levels        |
| `GSEA`                 | Gene Set Enrichment Analysis plot          |
| `GO_Network`           | GO enrichment cnetplot                     |
| `GO_Lollipop`          | Summary GO lollipop facet                  |
| `PPI`                  | STRING protein-protein interaction network |
| `Proximal_Motif_Bar`   | Proximal anchor motif enrichment barplot   |
| `Proximal_Motif_Logos` | Proximal anchor motif sequence logos       |
| `Proximal_Motif_Rank`  | Proximal anchor motif rank scatter         |
| `Distal_Motif_Bar`     | Distal anchor motif enrichment barplot     |
| `Distal_Motif_Logos`   | Distal anchor motif sequence logos         |
| `Distal_Motif_Rank`    | Distal anchor motif rank scatter           |

All ggplot-based plots accept standard **ggplot2** `+` layers for
recolouring or theme adjustment (see Module 2 for examples).

------------------------------------------------------------------------

### Module 5: IGV-Style Track Visualization

To visualize the local spatial interactome with high resolution, this
module generates an IGV-style multi-tiered genomic browser view via the
`plot_peaks_interactions` function. This visualization integrates
multi-omics datasets into three distinct, non-overlapping genomic
tracks, enabling intuitive assessment of 3D structural context:

- **Loop Track (Top)**: Depicts 3D chromatin interactions as Bezier arcs
  connecting paired spatial anchors. When `score_to_alpha = TRUE` is
  enabled, the opacity (alpha channel) of each arc is scaled
  proportionally to the quantitative interaction score.
- **Target Track (Middle)**: Renders user-provided genomic features
  (e.g., ChIP-seq peaks or GWAS SNPs) from the input BED file,
  facilitating direct visualization of their spatial overlap with the
  loop anchors.
- **Gene Track (Bottom)**: Retrieves the longest canonical transcript
  for each gene using standard Bioconductor annotation databases, with
  exons and introns mapped according to strand directionality.
  Overlapping gene models are vertically separated to avoid label
  collisions and improve readability.

#### Parameter Strategy & Core Inputs

The function provides fine-grained control over genomic locus selection,
data filtering, and visual aesthetics through the following arguments:

##### 1. Data Inputs & Genomic Coordinates

- **`bedpe_file`**: Character. Specifies the path to the BEDPE-format
  spatial interaction file. The function utilizes the first six columns
  for paired anchor coordinates and the seventh column (if present) as
  the quantitative interaction score.
- **`target_bed`**: Character (Optional). Specifies the path to a
  user-supplied genomic feature file in standard BED format, rendered in
  the middle track to assess overlap with 3D anchors.
- **`chr`, `from`, `to`**: Specifies the exact genomic window for
  visualization. If `chr` is omitted (`NULL`), the function defaults to
  the chromosome with the highest loop count in the BEDPE file.
- **`min_score`**: Numeric (Optional). Applies a threshold to filter the
  displayed loops based on the interaction score column.

##### 2. Annotation Dependencies

- **`species`**: Character. Specifies the reference genome assembly
  (supported: `"hg38"`, `"hg19"`, `"mm10"`, or `"mm9"`). This parameter
  triggers automatic loading of the corresponding Bioconductor `TxDb`
  and `OrgDb` annotation packages for gene track rendering.

##### 3. Aesthetic Mapping & Track Configuration

- **`score_to_alpha`**: Logical flag. When `TRUE`, the quantitative
  interaction score is mapped to the transparency of Bezier arcs,
  enabling visual differentiation of interaction strength.
- **`max_levels`**: Integer. Defines the maximum number of visible
  stacking bands for loop arcs. Overlapping loops are separated
  vertically up to this cap; denser regions are compressed into the
  available height while preserving relative layering.
- **`base_anchor_height`**: Numeric. Sets the vertical thickness of the
  rectangular anchors plotted at the base of the loop arcs.
- **Color Controls** (`loop_color`, `anchor_color`, `overlap_color`,
  `exon_color`, `intron_color`): Character strings specifying custom
  colors for structural elements across all three tracks.
- **`save_file`**: Character (Optional). Specifies the output file path
  and format (e.g., `".pdf"`) for exporting the final ggplot object.

``` r

if (requireNamespace("ggplot2", quietly = TRUE) &&
  requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE) &&
  requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  # If 'from' and 'to' are omitted, it automatically detects the densest viewport.
  track_plot <- plot_peaks_interactions(
    bedpe_file = f1,
    target_bed = atac_path,
    chr = "chr1",
    from = 11884299,
    to = 12106581,
    species = "hg38",
    save_file = file.path(out_dir, "Locus_Track.pdf")
  )
}
```

------------------------------------------------------------------------

## One-Click Parameterised Report

For rapid sharing with collaborators, a publication-ready HTML report
can be generated with a single function call. The template automatically
detects which pipeline stages are configured and hides unavailable
sections:

``` r

library(looplook)

# ── Required ──────────────────────────────────────────────────
#   bedpe_file         Chromatin loops in BEDPE format
# ── Required for annotation ───────────────────────────────────
#   target_bed         Genomic features (ATAC-seq peaks, etc.)
# ── Required for refinement + profiling ───────────────────────
#   expr_matrix_file   Normalised expression matrix (TPM/FPKM)
#   diff_file          Differential expression table (CSV/TSV)
#   metadata_file      Sample metadata with Group column
# ── Optional (defaults shown) ─────────────────────────────────
#   species = "hg38"   project_name = "looplook Analysis"
#   threshold = 1.0    unit_type = "TPM"
#   run_go = TRUE      run_ppi = FALSE     run_motif = FALSE
#   ... (see ?looplook_report for full list)
# ── Reproducibility ────────────────────────────────────────
#   The profiling stage uses R's random number generator for GSEA gene-set
#   down-sampling and motif background sampling. For fully reproducible
#   reports, call set.seed() before looplook_report():
#     set.seed(42)
#     looplook_report(...)

# Minimal: annotation only
looplook_report(
  bedpe_file   = system.file("extdata", "example_loops_1.bedpe", package = "looplook"),
  target_bed   = system.file("extdata", "example_peaks.bed", package = "looplook"),
  project_name = "Quick Annotation"
)

# Full pipeline: annotation + refinement + profiling
looplook_report(
  bedpe_file       = "my_loops.bedpe",
  target_bed       = "my_peaks.bed",
  expr_matrix_file = "my_tpm.txt",
  diff_file        = "my_deg.csv",
  metadata_file    = "my_coldata.txt",
  project_name     = "Integrative HiChIP Analysis",
  run_go           = TRUE
)

# Reuse precomputed results: skip annotation and jump to refinement + profiling.
# precomputed_res accepts a file path (.RData) or an in-memory list object.
# bedpe_file is optional when precomputed_res is set.

# Option A — save and reload from .RData:
res <- annotate_peaks_and_loops(
  bedpe_file = "my_loops.bedpe",
  target_bed = "my_peaks.bed",
  project_name = "MyProject"
)
save(res, file = "annotation_results.RData")

looplook_report(
  precomputed_res  = "annotation_results.RData",
  expr_matrix_file = "my_tpm.txt",
  diff_file        = "my_deg.csv",
  metadata_file    = "my_coldata.txt",
  run_go           = TRUE
)

# Option B — pass the object directly:
looplook_report(
  precomputed_res  = res,
  expr_matrix_file = "my_tpm.txt",
  diff_file        = "my_deg.csv",
  metadata_file    = "my_coldata.txt",
  run_go           = TRUE
)
```

The template is also accessible from the RStudio menu: **File → New File
→ R Markdown → From Template → looplook Report**.

------------------------------------------------------------------------

## Session Info

``` r

sessionInfo()
#> R version 4.6.0 (2026-04-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.4 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] looplook_0.99.8  dplyr_1.2.1      BiocStyle_2.40.0
#> 
#> loaded via a namespace (and not attached):
#>   [1] RColorBrewer_1.1-3          rstudioapi_0.18.0          
#>   [3] jsonlite_2.0.0              magrittr_2.0.5             
#>   [5] GenomicFeatures_1.64.0      farver_2.1.2               
#>   [7] rmarkdown_2.31              fs_2.1.0                   
#>   [9] BiocIO_1.22.0               fields_17.3                
#>  [11] ragg_1.5.2                  vctrs_0.7.3                
#>  [13] memoise_2.0.1               Rsamtools_2.28.0           
#>  [15] RCurl_1.98-1.18             base64enc_0.1-6            
#>  [17] htmltools_0.5.9             S4Arrays_1.12.0            
#>  [19] curl_7.1.0                  SparseArray_1.12.2         
#>  [21] Formula_1.2-5               sass_0.4.10                
#>  [23] bslib_0.10.0                htmlwidgets_1.6.4          
#>  [25] desc_1.4.3                  plyr_1.8.9                 
#>  [27] cachem_1.1.0                GenomicAlignments_1.48.0   
#>  [29] igraph_2.3.1                lifecycle_1.0.5            
#>  [31] pkgconfig_2.0.3             Matrix_1.7-5               
#>  [33] R6_2.6.1                    fastmap_1.2.0              
#>  [35] MatrixGenerics_1.24.0       digest_0.6.39              
#>  [37] colorspace_2.1-2            AnnotationDbi_1.74.0       
#>  [39] S4Vectors_0.50.0            regioneR_1.44.0            
#>  [41] bezier_1.1.2                textshaping_1.0.5          
#>  [43] Hmisc_5.2-5                 GenomicRanges_1.64.0       
#>  [45] RSQLite_3.52.0              httr_1.4.8                 
#>  [47] polyclip_1.10-7             abind_1.4-8                
#>  [49] compiler_4.6.0              bit64_4.8.0                
#>  [51] withr_3.0.2                 htmlTable_2.5.0            
#>  [53] S7_0.2.2                    backports_1.5.1            
#>  [55] BiocParallel_1.46.0         DBI_1.3.0                  
#>  [57] UpSetR_1.4.0                ggforce_0.5.0              
#>  [59] maps_3.4.3                  MASS_7.3-65                
#>  [61] DelayedArray_0.38.1         rjson_0.2.23               
#>  [63] tools_4.6.0                 foreign_0.8-91             
#>  [65] zip_2.3.3                   karyoploteR_1.38.0         
#>  [67] nnet_7.3-20                 glue_1.8.1                 
#>  [69] restfulr_0.0.16             InteractionSet_1.40.0      
#>  [71] grid_4.6.0                  checkmate_2.3.4            
#>  [73] cluster_2.1.8.2             generics_0.1.4             
#>  [75] gtable_0.3.6                BSgenome_1.80.0            
#>  [77] tidyr_1.3.2                 ensembldb_2.36.0           
#>  [79] data.table_1.18.4           XVector_0.52.0             
#>  [81] BiocGenerics_0.58.0         ggrepel_0.9.8              
#>  [83] pillar_1.11.1               stringr_1.6.0              
#>  [85] spam_2.11-3                 tweenr_2.0.3               
#>  [87] lattice_0.22-9              rtracklayer_1.72.0         
#>  [89] bit_4.6.0                   biovizBase_1.60.0          
#>  [91] tidyselect_1.2.1            Biostrings_2.80.0          
#>  [93] knitr_1.51                  gridExtra_2.3              
#>  [95] bookdown_0.46               ProtGenerics_1.44.0        
#>  [97] IRanges_2.46.0              Seqinfo_1.2.0              
#>  [99] SummarizedExperiment_1.42.0 stats4_4.6.0               
#> [101] xfun_0.57                   Biobase_2.72.0             
#> [103] matrixStats_1.5.0           stringi_1.8.7              
#> [105] UCSC.utils_1.8.0            lazyeval_0.2.3             
#> [107] yaml_2.3.12                 evaluate_1.0.5             
#> [109] codetools_0.2-20            cigarillo_1.2.0            
#> [111] tibble_3.3.1                BiocManager_1.30.27        
#> [113] cli_3.6.6                   rpart_4.1.27               
#> [115] systemfonts_1.3.2           jquerylib_0.1.4            
#> [117] dichromat_2.0-0.1           Rcpp_1.1.1-1.1             
#> [119] GenomeInfoDb_1.48.0         png_0.1-9                  
#> [121] XML_3.99-0.23               parallel_4.6.0             
#> [123] pkgdown_2.2.0               ggplot2_4.0.3              
#> [125] blob_1.3.0                  dotCall64_1.2              
#> [127] AnnotationFilter_1.36.0     bitops_1.0-9               
#> [129] viridisLite_0.4.3           VariantAnnotation_1.58.0   
#> [131] scales_1.4.0                purrr_1.2.2                
#> [133] openxlsx_4.2.8.1            crayon_1.5.3               
#> [135] bamsignals_1.43.0           rlang_1.2.0                
#> [137] KEGGREST_1.52.0
```
