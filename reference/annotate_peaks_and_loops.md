# Spatial mapping of 1D genomic features to 3D chromatin interaction targets

A dual-purpose analytical framework designed to integrate 1D genomic
features with 3D chromatin architecture.

1.  **Loop Annotation:** Classifies 3D spatial interactions (e.g.,
    Enhancer-Promoter, Promoter-Promoter) using a defined structural
    hierarchy.

2.  **Feature-to-Target Mapping:** Links 1D genomic features (e.g., GWAS
    risk SNPs, ATAC-seq peaks, ChIP-seq binding sites) to putative
    target genes via 3D chromatin contacts, providing a spatial
    complement to linear proximity-based assignments.

## Usage

``` r
annotate_peaks_and_loops(
  bedpe_file,
  target_bed = NULL,
  txdb = NULL,
  org_db = NULL,
  species = "hg38",
  tss_region = c(-2000, 2000),
  out_dir = "./results",
  expr_matrix_file = NULL,
  sample_columns = NULL,
  project_name = "HiChIP",
  color_palette = "Set2",
  karyo_bin_size = 1e+05,
  neighbor_hop = 0,
  hub_percentile = 0.95,
  write_output = TRUE,
  quiet = FALSE
)
```

## Arguments

- bedpe_file:

  Character. Path to a BEDPE file (at least 6 columns: chr1, start1,
  end1, chr2, start2, end2).

- target_bed:

  Optional path to a BED file of genomic features (e.g., ChIP-seq peaks,
  GWAS SNPs). When provided, these 1D regions are mapped to 3D target
  genes. Default: `NULL`.

- txdb:

  A [`TxDb`](https://rdrr.io/pkg/GenomicFeatures/man/TxDb-class.html)
  object, a package name string, or `NULL` to auto-resolve from
  `species`. Default: `NULL`.

- org_db:

  An `OrgDb` object, a package name string, or `NULL` to auto-resolve
  from `species`. Default: `NULL`.

- species:

  Character. Genome assembly used when `txdb` and `org_db` are `NULL`.
  One of `"hg38"`, `"hg19"`, `"mm10"`, `"mm9"`. Default: `"hg38"`.

- tss_region:

  Numeric vector of length 2. Promoter window around the TSS in bp.
  Default: `c(-2000, 2000)`.

- out_dir:

  Character. Output directory for the Excel results file. Default:
  `"./results"`.

- expr_matrix_file:

  Optional path to a normalised expression matrix (TPM/FPKM, genes ×
  samples). Enables expression-aware conflict resolution. Default:
  `NULL`.

- sample_columns:

  Character vector or integer indices. Columns in `expr_matrix_file` to
  average for baseline expression. Default: `NULL`.

- project_name:

  Character. Prefix for output files and plot titles. Default:
  `"HiChIP"`.

- color_palette:

  Character. RColorBrewer palette name. Default: `"Set2"`.

- karyo_bin_size:

  Integer. Bin width in bp for karyotype heatmaps. Default: `1e5`.

- neighbor_hop:

  Integer. k-hop ego-network expansion order via
  [`igraph::ego()`](https://r.igraph.org/reference/ego.html). `0`
  restricts to direct contacts. Default: `0`.

- hub_percentile:

  Numeric (0–1). Node-degree quantile for hub detection. Default:
  `0.95`.

- write_output:

  Logical. If `TRUE` (default), write the Excel workbook to `out_dir`.
  If `FALSE`, return results without creating directories or files.

- quiet:

  Logical. If `TRUE`, suppress progress messages while preserving
  warnings. Default: `FALSE`.

## Value

An invisible named list:

- `target_annotation` — Peak-to-gene assignments with
  `Assigned_Target_Genes_Filled` (3D-prioritised, falling back to
  nearest gene).

- `loop_annotation` — Annotated 3D interactome with
  `Putative_Target_Genes`.

- `anchor_loci_annotation` — Non-redundant anchor-locus genomic
  classifications after within-cluster interval reduction.

- `anchor_annotation` — Backward-compatible alias of
  `anchor_loci_annotation`.

- `promoter_centric_stats` — Gene-level connectivity statistics.

- `distal_element_stats` — Distal-element connectivity statistics.

- `plots` — Named list of ggplot objects (donut, karyotype, rose,
  flower).

- `plot_list` — Backward-compatible alias of `plots`.

If `write_output = TRUE`, also writes a multi-sheet Excel workbook to
`out_dir`.

## Details

**Mapping Strategy and Fallback Mechanism** The method prioritizes
physical 3D chromatin contacts. If a genomic feature overlaps an anchor
looping to a distal gene, the distal gene is assigned as the target. In
the absence of spatial loop evidence, the function implements a linear
proximity-based fallback to assign the nearest active local gene,
ensuring continuous annotation coverage.

**Hierarchical Conflict Resolution** To address complex loci where a
single anchor overlaps multiple promoters (e.g., dense gene clusters or
bidirectional promoters), the function executes a 3-step resolution:

1.  *Expression Filter:* Excludes transcriptionally silent genes using a
    user-provided expression matrix.

2.  *Biotype Prioritization:* Ranks remaining candidates by functional
    class:
    `Protein Coding > small-ncRNA (miRNA, snoRNA, snRNA, rRNA, scaRNA) > Antisense > lncRNA/ncRNA > Pseudogene`.

3.  *Expression Tiebreaker:* Resolves remaining ambiguities by
    designating the gene with the highest transcriptional abundance as
    the primary target.

**Network Topology Analysis**

- **Ego-Network Expansion (`neighbor_hop`):** Implements k-hop
  neighborhood expansion via
  [`igraph::ego()`](https://r.igraph.org/reference/ego.html). A value of
  `0` restricts to direct contacts, while `1` includes secondary
  contacts to capture broader regulatory cliques.

- **Hub Detection:** Utilizes a node-degree quantile threshold
  (`hub_percentile`) to identify highly connected regulatory elements.

## Examples

``` r
# Minimal runnable example for package checks
if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  txdb_example <- AnnotationDbi::loadDb(
    system.file("extdata", "hg19_knownGene_sample.sqlite", package = "GenomicFeatures")
  )
  bedpe_path <- tempfile(fileext = ".bedpe")
  writeLines(
    "chr6\t10412000\t10412600\tchr6\t10415000\t10415600",
    bedpe_path
  )

  res <- annotate_peaks_and_loops(
    bedpe_file = bedpe_path,
    txdb = txdb_example,
    org_db = "org.Hs.eg.db",
    species = "hg19",
    out_dir = tempdir(),
    project_name = "Quick_Example",
    write_output = FALSE,
    quiet = TRUE
  )
  head(res$loop_annotation, 1)
}
```
