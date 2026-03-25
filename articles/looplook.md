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

------------------------------------------------------------------------

## Installation

To ensure full functionality, install `looplook` along with its
recommended Bioconductor annotation dependencies:

``` r
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# Install annotation databases required for human (hg38) analysis
BiocManager::install(c("TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db"))

# Install looplook from GitHub
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
    samples. Only retains clusters detected in ≥ `min_consensus`
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
  threshold (e.g., ≥ 75% of replicates).
- **`gap`**: Defines the maximum allowable spatial distance (in base
  pairs) between loop anchors for them to be considered as part of the
  same physical cluster (default: `1000`).
- **`min_raw_score` vs. `min_score` (The Dual-Filter)**:
  - `min_raw_score` acts as a pre-filter applied to individual BEDPE
    files before clustering (e.g., removing singleton noise where PET
    count \< 2) to substantially reduce computational memory overhead.
  - `min_score` serves as a post-filter applied to the final merged
    chromatin interactome to ensure high-confidence interactions.
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

##### 3. `anchor_annotation` (Spatial Hub Footprints)

- **`Cluster_Locus_Genes`**: All genes from loop anchors within a
  spatially interconnected hub (defined by network-connected
  components), representing the full gene set associated with the 3D
  hub.

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

The module automatically generates up to 10 publication-grade PDF plots,
enabling multi-scale interrogation of the 3D interactome.

##### 1. Global 3D Network Profiling

*(Always Generated)*

- **Loop Type Proportions** (`*_Basic_Donut.pdf` &
  `*_Basic_Circular.pdf`): The Donut chart visualizes the raw physical
  frequency of each topological interaction. The Circular plot
  quantifies unique target genes per loop type, intuitively illustrating
  how certain loop types might disproportionately regulate target genes.

- **Karyotype Density Heatmaps** (`*_Basic_Karyo_Anchors.pdf` &
  `*_Basic_Karyo_LoopGenes.pdf`): Genome-wide ideograms illustrating the
  density of 3D anchors and their associated genes, highlighting
  architectural hotspots.

- **Topological Overlap** (`*_Basic_Flower.pdf`): A simplified flower
  plot showing the intersection of target genes (Core vs. Unique) across
  different loop types.

- **Anchor Genomic Distribution**
  (`*_Basic_Anchor_Genomic_Distribution.pdf`): Pie chart summarizing the
  genomic context of 3D anchors (e.g., Distal Intergenic, Intron,
  Promoter).

##### 2. Genomic Feature-to-3D Target Profiling

*(Requires `target_bed`)*

- **Target Connected Rose Chart** (`*_Target_Rose.pdf`): Coxcomb chart
  illustrating the dominant 3D topologies associated with user-defined
  genomic features.

- **Target Genomic Distribution Pies**
  (`*_Target_Genomic_Distribution.pdf` &
  `*_Target_Loop_Genomic_Distribution.pdf`): Compares the genomic
  distribution of input genomic features (peaks) versus those
  successfully linked via 3D loops. This comparison is insightful for
  demonstrating the enrichment of spatial interactions.

- **Target Genes Karyotype** (`*_Basic_Karyo_TargetGenes.pdf`):
  Chromosomal density map of final putative target genes, identifying
  functionally active multi-omics loci.

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
res_integrated <- system.file("extdata", "analysis_results.RData", package = "looplook")

if (res_integrated != "") {
  load(res_integrated)
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
}
```

### Example B: Expression-Aware Reclassification (Recommended)

This advanced mode is strongly recommended, as it preserves the valuable
3D structural framework while refining functional annotations. A silent
promoter anchor is reclassified as an enhancer (`eP`), acknowledging
that while its reference gene is transcriptionally inactive, this region
itself may function as an active regulatory element for a distant gene.

``` r
if (exists("res_integrated")) {
  refined_res <- refine_loop_anchors_by_expression(
    annotation_res = res_integrated,
    expr_matrix_file = expr_path,
    sample_columns = c("con1", "con2"),
    threshold = 1.0,
    unit_type = "TPM",
    reclassify_by_expression = TRUE, # Activates biological reclassification (e.g., P -> eP)
    out_dir = out_dir,
    project_name = "Example_Reclassified_Filter"
  )
}
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

- **Multi-Omics Sankey Diagram** (`*_Target_Sankey.pdf` / `.html`): Core
  integration diagnostic that maps the fate of inputted genomic features
  across three logical tiers: Genomic Region (L1) → 3D Loop Connection
  (L2) → Target Expression Status (L3). This intuitively reveals the
  proportion of GWAS/ChIP-seq peaks structurally connected to
  functionally active genes.

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

if (exists("refined_res")) {
  res_A <- profile_target_genes(
    annotation_res = res_integrated,
    diff_file = diff_path,
    lfc_col = "log2FoldChange",
    expr_matrix_file = expr_path,
    metadata_file = meta_path,
    target_source = c("loops", "targets"), # Analyzes both sources
    target_mapping_mode = "all",
    include_Filled = TRUE,
    use_nearest_gene = FALSE,
    project_name = "Scenario_A_Integrative",
    out_dir = out_dir,
    run_motif = FALSE, # Set TRUE in real analysis to scan JASPAR motifs
    run_go = FALSE, # Set TRUE in real analysis for pathway enrichment
    run_ppi = FALSE
  )
}
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
if (exists("refined_res")) {
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
    out_dir = out_dir,
    run_go = FALSE
  )
}
```

### Example C: 1D Linear Annotation Baseline (The Control)

By setting `use_nearest_gene = TRUE`, this mode intentionally bypasses
the 3D spatial topology and assigns features exclusively to their
nearest neighboring gene. It serves as a conventional “baseline
control.” Comparing results against Example A or B enables rigorous
demonstration of the novel functional insights gained from 3D chromatin
looping analyses relative to conventional method.

``` r
if (exists("refined_res")) {
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
    out_dir = out_dir,
    run_go = FALSE
  )
}
```

#### Deep Dive: Functional Visualizations

This module executes an extensive downstream pipeline that generates a
complete, publication-ready figure panel. Rather than simple descriptive
summaries, these visualizations are designed to test research hypotheses
related to multi-omics integration.

##### 1. Expression & Topological Dynamics

- **Transcriptional Profiling** (`*_LFC_Violin.pdf` &
  `*_Expression_Heatmap.pdf`): Integrates quantitative transcriptome
  data (e.g., RNA-seq) and differential expression statistics to
  evaluate whether 3D-annotated target genes exhibit statistically
  significant expression shifts relative to the global genome
  background. This provides a quantitative measure of functional
  relevance and correlative evidence of regulatory modulation.
- **Topological-Transcriptional Association**
  (`*_Connectivity_Scatter.pdf` & `*_Raincloud_*.pdf`): Cross-references
  3D connectivity (e.g., number of linked enhancers) with transcript
  abundance and fold changes. This may help to test the “hub hypothesis”
  (whether increased spatial connectivity correlates with more dynamic
  transcriptional output).

##### 2. Regulatory Motif Scanning

- **Spatially Asymmetric Motif Signatures** (`*_Motif_*_RankScatter.pdf`
  & `*_Logos.pdf`): Separately analyses transcription factor motifs
  enriched in the promoter-proximal anchors and distal (enhancer)
  anchors.

##### 3. Pathway & Network Enrichment

- **Spatial Gene Set Enrichment** (`*_GSEA.pdf`): Projects 3D-derived
  target gene sets onto the globally ranked differential expression
  landscape. Rather than relying on pre-defined pathways, it treats the
  specific spatial network as a biologically meaningful gene set, and
  examines whether the 3D interactome contributes to the observed
  transcriptomic phenotype.
- **Core Effector Networks** (`*_GO_Network.pdf`, `*_GO_Lollipop_*.pdf`
  & `*_PPI_Network_*.pdf`): Derived from Gene Ontology and STRING
  analyses. Rather than simply listing enriched terms, the Divergent
  Concept Networks and PPI graphs link key pathways to high-magnitude
  differentially expressed hub genes, bridging spatial chromatin
  architecture, transcriptional regulation, and functional protein
  networks.

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
- **`max_levels`**: Integer. Defines the maximum vertical stacking limit
  for the loop arcs to manage rendering height and readability.
- **`base_anchor_height`**: Numeric. Sets the vertical thickness of the
  rectangular anchors plotted at the base of the loop arcs.
- **Color Controls** (`loop_color`, `anchor_color`, `overlap_color`,
  `exon_color`, `intron_color`): Character strings specifying custom
  colors for structural elements across all three tracks.
- **`save_file`**: Character (Optional). Specifies the output file path
  and format (e.g., `".pdf"`) for exporting the final ggplot object.

``` r
if (requireNamespace("ggplot2", quietly = TRUE)) {
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

## Session Info

``` r
sessionInfo()
#> R version 4.5.3 (2026-03-11)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.3 LTS
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
#> [1] stats4    stats     graphics  grDevices utils     datasets  methods  
#> [8] base     
#> 
#> other attached packages:
#> [1] org.Hs.eg.db_3.22.0  AnnotationDbi_1.72.0 IRanges_2.44.0      
#> [4] S4Vectors_0.48.0     Biobase_2.70.0       BiocGenerics_0.56.0 
#> [7] generics_0.1.4       looplook_0.99.4      BiocStyle_2.38.0    
#> 
#> loaded via a namespace (and not attached):
#>   [1] fs_2.0.0                                
#>   [2] ProtGenerics_1.42.0                     
#>   [3] matrixStats_1.5.0                       
#>   [4] bitops_1.0-9                            
#>   [5] enrichplot_1.30.5                       
#>   [6] httr_1.4.8                              
#>   [7] RColorBrewer_1.1-3                      
#>   [8] InteractionSet_1.38.0                   
#>   [9] tools_4.5.3                             
#>  [10] backports_1.5.0                         
#>  [11] R6_2.6.1                                
#>  [12] lazyeval_0.2.2                          
#>  [13] withr_3.0.2                             
#>  [14] gridExtra_2.3                           
#>  [15] cli_3.6.5                               
#>  [16] textshaping_1.0.5                       
#>  [17] scatterpie_0.2.6                        
#>  [18] labeling_0.4.3                          
#>  [19] sass_0.4.10                             
#>  [20] S7_0.2.1                                
#>  [21] pkgdown_2.2.0                           
#>  [22] Rsamtools_2.26.0                        
#>  [23] systemfonts_1.3.2                       
#>  [24] yulab.utils_0.2.4                       
#>  [25] foreign_0.8-91                          
#>  [26] DOSE_4.4.0                              
#>  [27] R.utils_2.13.0                          
#>  [28] dichromat_2.0-0.1                       
#>  [29] plotrix_3.8-14                          
#>  [30] BSgenome_1.78.0                         
#>  [31] maps_3.4.3                              
#>  [32] rstudioapi_0.18.0                       
#>  [33] RSQLite_2.4.6                           
#>  [34] gridGraphics_0.5-1                      
#>  [35] TxDb.Hsapiens.UCSC.hg19.knownGene_3.22.1
#>  [36] BiocIO_1.20.0                           
#>  [37] gtools_3.9.5                            
#>  [38] dplyr_1.2.0                             
#>  [39] zip_2.3.3                               
#>  [40] GO.db_3.22.0                            
#>  [41] Matrix_1.7-4                            
#>  [42] abind_1.4-8                             
#>  [43] R.methodsS3_1.8.2                       
#>  [44] lifecycle_1.0.5                         
#>  [45] yaml_2.3.12                             
#>  [46] SummarizedExperiment_1.40.0             
#>  [47] gplots_3.3.0                            
#>  [48] qvalue_2.42.0                           
#>  [49] SparseArray_1.10.9                      
#>  [50] grid_4.5.3                              
#>  [51] blob_1.3.0                              
#>  [52] promises_1.5.0                          
#>  [53] crayon_1.5.3                            
#>  [54] ggtangle_0.1.1                          
#>  [55] lattice_0.22-9                          
#>  [56] cowplot_1.2.0                           
#>  [57] GenomicFeatures_1.62.0                  
#>  [58] cigarillo_1.0.0                         
#>  [59] chromote_0.5.1                          
#>  [60] KEGGREST_1.50.0                         
#>  [61] pillar_1.11.1                           
#>  [62] knitr_1.51                              
#>  [63] fgsea_1.36.2                            
#>  [64] GenomicRanges_1.62.1                    
#>  [65] rjson_0.2.23                            
#>  [66] boot_1.3-32                             
#>  [67] codetools_0.2-20                        
#>  [68] fastmatch_1.1-8                         
#>  [69] glue_1.8.0                              
#>  [70] ggiraph_0.9.6                           
#>  [71] ggfun_0.2.0                             
#>  [72] fontLiberation_0.1.0                    
#>  [73] data.table_1.18.2.1                     
#>  [74] vctrs_0.7.2                             
#>  [75] png_0.1-9                               
#>  [76] treeio_1.34.0                           
#>  [77] spam_2.11-3                             
#>  [78] gtable_0.3.6                            
#>  [79] cachem_1.1.0                            
#>  [80] xfun_0.57                               
#>  [81] openxlsx_4.2.8.1                        
#>  [82] S4Arrays_1.10.1                         
#>  [83] Seqinfo_1.0.0                           
#>  [84] fields_17.1                             
#>  [85] nlme_3.1-168                            
#>  [86] ggtree_4.0.5                            
#>  [87] bit64_4.6.0-1                           
#>  [88] fontquiver_0.2.1                        
#>  [89] UpSetR_1.4.0                            
#>  [90] GenomeInfoDb_1.46.2                     
#>  [91] data.tree_1.2.0                         
#>  [92] bslib_0.10.0                            
#>  [93] KernSmooth_2.23-26                      
#>  [94] otel_0.2.0                              
#>  [95] rpart_4.1.24                            
#>  [96] colorspace_2.1-2                        
#>  [97] DBI_1.3.0                               
#>  [98] Hmisc_5.2-5                             
#>  [99] nnet_7.3-20                             
#> [100] tidyselect_1.2.1                        
#> [101] processx_3.8.6                          
#> [102] bit_4.6.0                               
#> [103] compiler_4.5.3                          
#> [104] curl_7.0.0                              
#> [105] htmlTable_2.4.3                         
#> [106] bezier_1.1.2                            
#> [107] desc_1.4.3                              
#> [108] fontBitstreamVera_0.1.1                 
#> [109] DelayedArray_0.36.0                     
#> [110] bookdown_0.46                           
#> [111] rtracklayer_1.70.1                      
#> [112] checkmate_2.3.4                         
#> [113] scales_1.4.0                            
#> [114] caTools_1.18.3                          
#> [115] ChIPseeker_1.46.1                       
#> [116] rappdirs_0.3.4                          
#> [117] stringr_1.6.0                           
#> [118] digest_0.6.39                           
#> [119] rmarkdown_2.30                          
#> [120] XVector_0.50.0                          
#> [121] htmltools_0.5.9                         
#> [122] pkgconfig_2.0.3                         
#> [123] base64enc_0.1-6                         
#> [124] MatrixGenerics_1.22.0                   
#> [125] regioneR_1.42.0                         
#> [126] fastmap_1.2.0                           
#> [127] ensembldb_2.34.0                        
#> [128] rlang_1.1.7                             
#> [129] htmlwidgets_1.6.4                       
#> [130] UCSC.utils_1.6.1                        
#> [131] farver_2.1.2                            
#> [132] jquerylib_0.1.4                         
#> [133] karyoploteR_1.36.0                      
#> [134] jsonlite_2.0.0                          
#> [135] BiocParallel_1.44.0                     
#> [136] GOSemSim_2.36.0                         
#> [137] R.oo_1.27.1                             
#> [138] VariantAnnotation_1.56.0                
#> [139] RCurl_1.98-1.18                         
#> [140] magrittr_2.0.4                          
#> [141] Formula_1.2-5                           
#> [142] ggplotify_0.1.3                         
#> [143] dotCall64_1.2                           
#> [144] patchwork_1.3.2                         
#> [145] Rcpp_1.1.1                              
#> [146] ape_5.8-1                               
#> [147] ggnewscale_0.5.2                        
#> [148] bamsignals_1.42.0                       
#> [149] gdtools_0.5.0                           
#> [150] stringi_1.8.7                           
#> [151] MASS_7.3-65                             
#> [152] plyr_1.8.9                              
#> [153] parallel_4.5.3                          
#> [154] ggrepel_0.9.8                           
#> [155] Biostrings_2.78.0                       
#> [156] splines_4.5.3                           
#> [157] ps_1.9.1                                
#> [158] igraph_2.2.2                            
#> [159] reshape2_1.4.5                          
#> [160] XML_3.99-0.23                           
#> [161] evaluate_1.0.5                          
#> [162] biovizBase_1.58.0                       
#> [163] BiocManager_1.30.27                     
#> [164] tweenr_2.0.3                            
#> [165] networkD3_0.4.1                         
#> [166] tidyr_1.3.2                             
#> [167] webshot2_0.1.2                          
#> [168] purrr_1.2.1                             
#> [169] polyclip_1.10-7                         
#> [170] ggplot2_4.0.2                           
#> [171] ggforce_0.5.0                           
#> [172] restfulr_0.0.16                         
#> [173] AnnotationFilter_1.34.0                 
#> [174] tidytree_0.4.7                          
#> [175] tidydr_0.0.6                            
#> [176] later_1.4.8                             
#> [177] viridisLite_0.4.3                       
#> [178] ragg_1.5.2                              
#> [179] TxDb.Hsapiens.UCSC.hg38.knownGene_3.22.0
#> [180] tibble_3.3.1                            
#> [181] websocket_1.4.4                         
#> [182] aplot_0.2.9                             
#> [183] memoise_2.0.1                           
#> [184] GenomicAlignments_1.46.0                
#> [185] cluster_2.1.8.2
```
