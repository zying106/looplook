# A core analytical engine enabling high-resolution spatial mapping of genomic features/peaks to 3D chromatin interaction targets

A highly sophisticated, dual-purpose 3D genomic engine designed to:

1.  **Annotate Chromatin Loops:** Accurately classify 3D spatial
    interactions (e.g., Enhancer-Promoter, Promoter-Promoter) using a
    strict structural hierarchy.

2.  **Map External Features (genomic features to 3D):** Act as a **"3D
    Spatial Bridge"** to link auxiliary 1D genomic features (e.g., GWAS
    risk SNPs, ATAC-seq peaks, ChIP-seq binding sites) to their genuine
    regulated target genes, fundamentally outperforming simple "nearest
    gene" heuristics.

**Core Philosophy: Bridging the genomic features-3D Gap** Traditional
annotations typically assign non-coding variants to the nearest linear
gene. This engine prioritizes **physical 3D chromatin contacts**. If a
GWAS SNP or enhancer lands in an anchor looping to a distal gene, the
distal gene is accurately assigned. If no spatial loop exists, the
engine intelligently falls back to the nearest active local gene (Smart
Fallback), ensuring comprehensive and gapless annotation coverage.

**Key Algorithmic Innovations:**

- **Biotype & Expression-Aware Conflict Resolution:** When an anchor
  overlaps multiple promoters (e.g., dense gene loci or bidirectional
  promoters), it executes a rigorous 3-step resolution:

  1.  *Expression Pre-filter:* Eliminates transcriptionally silent genes
      based on the user-provided expression matrix.

  2.  *Functional Biotype Prioritization::* Prioritizes the remaining
      candidates by functional class in the following order:
      `Protein Coding > Antisense > lncRNA > Pseudogene`.

  3.  *Dominant Expression Tiebreaker:* Designates the gene with the
      highest transcriptional abundance as the target gene to further
      resolve any remaining mapping ambiguities.

- **Dynamic Topology Control (`neighbor_hop`):** Controls network
  diffusion depth. Allows signal propagation from a SNP -\> Enhancer -\>
  Hub -\> Target Gene, which is ideal for uncovering multi-way
  super-enhancer cliques.

- **Hub Detection:** Calculates precise node degrees to identify
  high-connectivity Promoter and Distal hubs based on the
  `hub_percentile`.

## Usage

``` r
annotate_peaks_and_loops(
  bedpe_file,
  target_bed = NULL,
  species = "hg38",
  tss_region = c(-2000, 2000),
  out_dir = "./results",
  expr_matrix_file = NULL,
  sample_columns = NULL,
  project_name = "HiChIP",
  color_palette = "Set2",
  karyo_bin_size = 1e+05,
  neighbor_hop = 0,
  hub_percentile = 0.95
)
```

## Arguments

- bedpe_file:

  Character. Path to the BEDPE file containing 3D loop coordinates (The
  "3D Bridge").

- target_bed:

  Optional character. Path to an auxiliary BED file (e.g., GWAS summary
  stats, eQTLs, ChIP-seq/ATAC-seq peaks). If provided, maps these 1D
  regulatory regions to 3D target genes.

- species:

  Character. Reference genome build. Supported: `"hg38"`, `"hg19"`,
  `"mm10"`, `"mm9"`.

- tss_region:

  Numeric vector. Defines the promoter window around the TSS. Default:
  `c(-2000, 2000)`.

- out_dir:

  Character. Directory for saving generated PDF plots and Excel results.

- expr_matrix_file:

  Optional character. Path to the normalized RNA-seq matrix (e.g.,
  TPM/FPKM). **Highly recommended** as it fuels the intelligent conflict
  resolution and filters out silent elements.

- sample_columns:

  Character vector or Integer indices. Specifies which columns in
  `expr_matrix_file` to average for the baseline expression reference.

- project_name:

  Character. Prefix for all output files and plot titles.

- color_palette:

  Character. Color brewer palette name for visualizations (e.g.,
  `"Set2"`).

- karyo_bin_size:

  Numeric. Bin size for Karyotype genomic density heatmaps. Default:
  `1e5` (100kb).

- neighbor_hop:

  Integer. **Network Diffusion Depth**. `0` (Default) captures direct
  physical connections only; `1` captures indirect 1-hop hub-mediated
  connections.

- hub_percentile:

  Numeric (0-1). Top percentile threshold for defining high-connectivity
  "Regulatory Hubs" (default: `0.95`).

## Value

An invisible list of comprehensive data frames (also auto-saved as a
multi-sheet `.xlsx`):

- `target_annotation`: Detailed annotation of `target_bed`. Features the
  important column: `Assigned_Target_Genes_Filled` (Loop-prioritized
  target, gracefully falling back to the local nearest gene if
  unlooped).

- `loop_annotation`: The fully annotated 3D network. Features
  `Putative_Target_Genes` capturing precisely resolved regulatory
  targets.

- `promoter_centric_stats`: Gene-level topological statistics. Contains
  structural node degrees and identifies `Is_High_Connectivity_Gene`.

- `distal_element_stats`: Enhancer-level topological statistics,
  highlighting critical regulatory anchors orchestrating multiple
  promoters.

- `anchor_annotation`: Locus-level summarization of all 1D anchor
  footprints.

## Examples

``` r
# 1. Get paths to example files included in the package
bedpe_path <- system.file("extdata", "example_loops_1.bedpe", package = "looplook")
bed_path <- system.file("extdata", "example_peaks.bed", package = "looplook")
expr_path <- system.file("extdata", "example_tpm.txt", package = "looplook")

# 2. Check if files and required annotation databases exist
if (bedpe_path != "" && bed_path != "" && expr_path != "" &&
  requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE) &&
  requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  # =========================================================================
  # Example A: Integrative Analysis (Loops + GWAS/ChIP-seq + Expression)
  # =========================================================================
  # Here, `target_bed` can represent GWAS risk SNPs, ATAC-seq peaks, or
  # ChIP-seq peaks. The engine bridges these 1D regions to their 3D genes.
  res_integrated <- annotate_peaks_and_loops(
    bedpe_file = bedpe_path,
    target_bed = bed_path,
    expr_matrix_file = expr_path,
    sample_columns = c("con1", "con2"),
    species = "hg38",
    tss_region = c(-2000, 2000),
    out_dir = tempdir(),
    color_palette = "Set2",
    karyo_bin_size = 1e5,
    neighbor_hop = 0,
    hub_percentile = 0.95,
    project_name = "Example_HiChIP_Integrative"
  )

  # View integrated annotations linking external peaks/SNPs to target genes
  head(res_integrated$target_annotation)

  # =========================================================================
  # Example B: Deep Analysis of Loops ONLY (No auxiliary peaks)
  # =========================================================================
  res_loops_only <- annotate_peaks_and_loops(
    bedpe_file = bedpe_path,
    target_bed = NULL,
    species = "hg38",
    expr_matrix_file = expr_path,
    sample_columns = c("con1", "con2"),
    tss_region = c(-2000, 2000),
    out_dir = tempdir(),
    color_palette = "Set1",
    karyo_bin_size = 1e5,
    neighbor_hop = 0,
    hub_percentile = 0.95,
    project_name = "Example_Loops_Only"
  )

  # View standalone loop annotations
  head(res_loops_only$loop_annotation)
}
#> 
#> Step 0: Loading expression data...
#>     >>> Expression loaded for 943 genes.
#> Step 1: Reading BEDPE file...
#> Step 2: Clustering loops...
#> Step 3: Biological Classification & Topology...
#> 
#> >> Using Genome: hg38 ...
#> >> Using Genome: hg38 ...
#> >> Using Genome: hg38 ...
#> Loading required package: org.Hs.eg.db
#> Loading required package: AnnotationDbi
#> Loading required package: stats4
#> Loading required package: BiocGenerics
#> Loading required package: generics
#> 
#> Attaching package: 'generics'
#> The following objects are masked from 'package:base':
#> 
#>     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
#>     setequal, union
#> 
#> Attaching package: 'BiocGenerics'
#> The following objects are masked from 'package:stats':
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from 'package:base':
#> 
#>     Filter, Find, Map, Position, Reduce, anyDuplicated, aperm, append,
#>     as.data.frame, basename, cbind, colnames, dirname, do.call,
#>     duplicated, eval, evalq, get, grep, grepl, is.unsorted, lapply,
#>     mapply, match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
#>     rank, rbind, rownames, sapply, saveRDS, table, tapply, unique,
#>     unsplit, which.max, which.min
#> Loading required package: Biobase
#> Welcome to Bioconductor
#> 
#>     Vignettes contain introductory material; view with
#>     'browseVignettes()'. To cite Bioconductor, see
#>     'citation("Biobase")', and for packages 'citation("pkgname")'.
#> Loading required package: IRanges
#> Loading required package: S4Vectors
#> 
#> Attaching package: 'S4Vectors'
#> The following object is masked from 'package:utils':
#> 
#>     findMatches
#> The following objects are masked from 'package:base':
#> 
#>     I, expand.grid, unname
#> 'select()' returned 1:many mapping between keys and columns
#>   2169 genes were dropped because they have exons located on both strands of
#>   the same reference sequence or on more than one reference sequence, so cannot
#>   be represented by a single genomic range.
#>   Use 'single.strand.genes.only=FALSE' to get all the genes in a GRangesList
#>   object, or use suppressMessages() to suppress this message.
#> Warning: GRanges object contains 8 out-of-bound ranges located on sequences
#>   chr4_GL000257v2_alt, chr16_KI270728v1_random, chr22_KI270731v1_random,
#>   chr3_GL000221v1_random, chr1_KI270706v1_random, chr16_GL383556v1_alt, and
#>   chrUn_KI270748v1. Note that ranges located on a sequence whose length is
#>   unknown (NA) or on a circular sequence are not considered out-of-bound (use
#>   seqlengths() and isCircular() to get the lengths and circularity flags of the
#>   underlying sequences). You can use trim() to trim these ranges. See
#>   ?`trim,GenomicRanges-method` for more information.
#> 'select()' returned 1:1 mapping between keys and columns
#>     Calculating Topology (Hops)...
#> Step 4: Constructing Loop Tables...
#> >> Using Genome: hg38 ...
#> >> Using Genome: hg38 ...
#> >> Using Genome: hg38 ...
#> 'select()' returned 1:many mapping between keys and columns
#>     Generating Promoter Centric Stats...
#>     Generating Distal Element Stats...
#> Step 5: Integrating Target Annotations...
#> >> Using Genome: hg38 ...
#> >> Using Genome: hg38 ...
#> >> Using Genome: hg38 ...
#> 'select()' returned 1:many mapping between keys and columns
#>     Refining Target annotation...
#>   2169 genes were dropped because they have exons located on both strands of
#>   the same reference sequence or on more than one reference sequence, so cannot
#>   be represented by a single genomic range.
#>   Use 'single.strand.genes.only=FALSE' to get all the genes in a GRangesList
#>   object, or use suppressMessages() to suppress this message.
#> Warning: GRanges object contains 8 out-of-bound ranges located on sequences
#>   chr4_GL000257v2_alt, chr16_KI270728v1_random, chr22_KI270731v1_random,
#>   chr3_GL000221v1_random, chr1_KI270706v1_random, chr16_GL383556v1_alt, and
#>   chrUn_KI270748v1. Note that ranges located on a sequence whose length is
#>   unknown (NA) or on a circular sequence are not considered out-of-bound (use
#>   seqlengths() and isCircular() to get the lengths and circularity flags of the
#>   underlying sequences). You can use trim() to trim these ranges. See
#>   ?`trim,GenomicRanges-method` for more information.
#> 'select()' returned 1:1 mapping between keys and columns
#> Step 6: Generating Visualizations...
#> Warning: Ignoring unknown parameters: `size`
#>  Saved: /tmp/RtmpRjnvH9/Example_HiChIP_Integrative_Basic_Circular.pdf
#>   2169 genes were dropped because they have exons located on both strands of
#>   the same reference sequence or on more than one reference sequence, so cannot
#>   be represented by a single genomic range.
#>   Use 'single.strand.genes.only=FALSE' to get all the genes in a GRangesList
#>   object, or use suppressMessages() to suppress this message.
#> 'select()' returned 1:1 mapping between keys and columns
#>     Saved Heatmap: /tmp/RtmpRjnvH9/Example_HiChIP_Integrative_Basic_Karyo_LoopGenes.pdf
#>     Saved Heatmap: /tmp/RtmpRjnvH9/Example_HiChIP_Integrative_Basic_Karyo_Anchors.pdf
#> Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
#> ℹ Please use `linewidth` instead.
#> ℹ The deprecated feature was likely used in the looplook package.
#>   Please report the issue at <https://github.com/zying106/looplook/issues>.
#>     Saved (Simplified Flower Plot with inner counts): /tmp/RtmpRjnvH9/Example_HiChIP_Integrative_Basic_Flower.pdf
#>     Plotting Target Visualizations...
#>   2169 genes were dropped because they have exons located on both strands of
#>   the same reference sequence or on more than one reference sequence, so cannot
#>   be represented by a single genomic range.
#>   Use 'single.strand.genes.only=FALSE' to get all the genes in a GRangesList
#>   object, or use suppressMessages() to suppress this message.
#> 'select()' returned 1:1 mapping between keys and columns
#>     Saved Heatmap: /tmp/RtmpRjnvH9/Example_HiChIP_Integrative_Basic_Karyo_TargetGenes.pdf
#>     Generating Pie Chart 0: All Anchors Genomic Distribution...
#>     Generating Pie Chart 1: All Targets Distribution...
#>     Generating Pie Chart 2: Loop-Connected Targets Distribution...
#> Step 7: Exporting to Excel...
#>     Excel file saved.
#> Analysis Complete.
#> Step 0: Loading expression data...
#>     >>> Expression loaded for 943 genes.
#> Step 1: Reading BEDPE file...
#> Step 2: Clustering loops...
#> Step 3: Biological Classification & Topology...
#> >> Using Genome: hg38 ...
#> >> Using Genome: hg38 ...
#> >> Using Genome: hg38 ...
#> 'select()' returned 1:many mapping between keys and columns
#>   2169 genes were dropped because they have exons located on both strands of
#>   the same reference sequence or on more than one reference sequence, so cannot
#>   be represented by a single genomic range.
#>   Use 'single.strand.genes.only=FALSE' to get all the genes in a GRangesList
#>   object, or use suppressMessages() to suppress this message.
#> Warning: GRanges object contains 8 out-of-bound ranges located on sequences
#>   chr4_GL000257v2_alt, chr16_KI270728v1_random, chr22_KI270731v1_random,
#>   chr3_GL000221v1_random, chr1_KI270706v1_random, chr16_GL383556v1_alt, and
#>   chrUn_KI270748v1. Note that ranges located on a sequence whose length is
#>   unknown (NA) or on a circular sequence are not considered out-of-bound (use
#>   seqlengths() and isCircular() to get the lengths and circularity flags of the
#>   underlying sequences). You can use trim() to trim these ranges. See
#>   ?`trim,GenomicRanges-method` for more information.
#> 'select()' returned 1:1 mapping between keys and columns
#>     Calculating Topology (Hops)...
#> Step 4: Constructing Loop Tables...
#> >> Using Genome: hg38 ...
#> >> Using Genome: hg38 ...
#> >> Using Genome: hg38 ...
#> 'select()' returned 1:many mapping between keys and columns
#>     Generating Promoter Centric Stats...
#>     Generating Distal Element Stats...
#> Step 6: Generating Visualizations...
#> Warning: Ignoring unknown parameters: `size`
#>  Saved: /tmp/RtmpRjnvH9/Example_Loops_Only_Basic_Circular.pdf
#>   2169 genes were dropped because they have exons located on both strands of
#>   the same reference sequence or on more than one reference sequence, so cannot
#>   be represented by a single genomic range.
#>   Use 'single.strand.genes.only=FALSE' to get all the genes in a GRangesList
#>   object, or use suppressMessages() to suppress this message.
#> 'select()' returned 1:1 mapping between keys and columns
#>     Saved Heatmap: /tmp/RtmpRjnvH9/Example_Loops_Only_Basic_Karyo_LoopGenes.pdf
#>     Saved Heatmap: /tmp/RtmpRjnvH9/Example_Loops_Only_Basic_Karyo_Anchors.pdf
#>     Saved (Simplified Flower Plot with inner counts): /tmp/RtmpRjnvH9/Example_Loops_Only_Basic_Flower.pdf
#> Step 7: Exporting to Excel...
#>     Excel file saved.
#> Analysis Complete.
#> # A tibble: 6 × 16
#>   loop_ID chr1     start1      end1 chr2     start2    end2 cluster_id loop_type
#>   <chr>   <chr>     <int>     <int> <chr>     <int>   <int> <chr>      <chr>    
#> 1 L1      chr1  112721687 112724586 chr1  112791557  1.13e8 1          E-P      
#> 2 L2      chr1  113867048 113875597 chr1  113927096  1.14e8 1          P-P      
#> 3 L3      chr1  115060199 115071191 chr1  116162680  1.16e8 1          P-P      
#> 4 L4      chr1  118807309 118811756 chr1  119277565  1.19e8 1          E-G      
#> 5 L5      chr1  112399375 112406783 chr1  113756152  1.14e8 1          G-P      
#> 6 L6      chr1   12793660  12795281 chr1   12809516  1.28e7 2          E-P      
#> # ℹ 7 more variables: anchor1_gene <chr>, anchor1_type <chr>,
#> #   anchor2_gene <chr>, anchor2_type <chr>, Cluster_All_Genes <chr>,
#> #   All_Anchor_Genes <chr>, Putative_Target_Genes <chr>
```
