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
  min_expr = 0,
  conflict_strategy = c("biotype_first", "expression_first"),
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

- min_expr:

  Numeric. Minimum expression value for a gene to be considered active
  during anchor-level conflict resolution. Used only when
  `expr_matrix_file` is provided. Default: `0` (any detectable
  expression qualifies). Increase to `1` or higher to require stronger
  evidence. See
  [`refine_loop_anchors_by_expression`](https://zying106.github.io/looplook/reference/refine_loop_anchors_by_expression.md)
  for the separate `threshold` parameter that controls promoter
  reclassification.

- conflict_strategy:

  Character. Conflict resolution order for overlapping gene assignments.
  `"biotype_first"` (default): select the best biotype tier first, then
  apply expression filtering within that tier. `"expression_first"`:
  apply expression filtering across all biotypes first, then pick the
  best biotype among expressed candidates.

- write_output:

  Logical. If `TRUE` (default), write the Excel workbook to `out_dir`.
  If `FALSE`, return results without creating directories or files.

- quiet:

  Logical. If `TRUE`, suppress progress messages while preserving
  warnings. Default: `FALSE`.

## Value

An invisible named list:

- `target_annotation` — Target features (peaks) with gene assignments.
  Key columns include:

  - `All_Loop_Connected_Genes`: All genes from loop-connected anchors
    (P/G types).

  - `Regulated_promoter_genes`: Promoter genes supported by loop-anchor
    context.

  - `Assigned_Target_Genes`: Promoter-first 3D assignment (prioritises P
    \> G \> E).

  - `*_Filled` variants: Linear nearest-gene fallback when strict 3D
    assignments are empty.

  - `Regulated_promoter_Evidence`: Provenance of
    `Regulated_promoter_genes` (e.g., `local_promoter_overlap`,
    `direct_opposite_promoter`). **Read with**
    `Regulated_promoter_genes`; do not cross-reference with
    `Assigned_Target_Genes` or other columns.

  - `Regulated_promoter_Fallback_Evidence`: Provenance of
    `Regulated_promoter_genes_Filled`. **Read with**
    `Regulated_promoter_genes_Filled`; indicates which `*_Filled` column
    supplied the fallback gene.

- `target_gene_links` — Long-format peak-gene provenance table. Each row
  records one peak-gene linkage with full provenance. **Read**
  `evidence`, `anchor_role`, and `gene_role` **together as a group** —
  they jointly describe how each gene was assigned to each peak; do not
  interpret any one column in isolation.

  - `input_id`, `loop_ID`, `anchor_id`: Identifiers.

  - `gene`: Linked gene symbol.

  - `gene_role`: `"promoter"`, `"gene_body"`, or `"linear_annotation"`.

  - `source`: `"loop_anchor"` (3D-derived) or `"linear_annotation"`
    (nearest gene).

  - `evidence`: Provenance label — `"local_promoter_overlap"` (peak
    overlaps anchor promoter), `"direct_opposite_promoter"` (opposite
    anchor is promoter), `"gene_body_context"` (gene body linkage),
    `"expanded_promoter_loop"` (via ego-network expansion),
    `"linear_annotation"` (direct nearest gene), or `"linear_fallback"`
    (filled when 3D assignment was empty).

  - `anchor_role`: `"local_anchor"`, `"opposite_anchor"`,
    `"expanded_anchor"`, or `"linear_annotation"`.

  - `used_as_fallback`: Logical. `TRUE` when this link was added via the
    `*_Filled` linear nearest-gene fallback mechanism.

  - `in_regulated_promoter` through `in_assigned_target_filled`: Logical
    membership flags indicating which target annotation column(s) this
    gene appears in.

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
physical 3D chromatin contacts while keeping strict and fallback
semantics separate. `Regulated_promoter_genes` reports promoter genes
supported by loop-anchor context, `Assigned_Target_Genes` preserves the
historical promoter-first 3D assignment, and `*_Filled` columns add a
linear nearest-gene fallback only when strict 3D assignments are empty.
`Regulated_promoter_Evidence`, `Regulated_promoter_Fallback_Evidence`,
and `target_gene_links` provide row-level provenance for these
decisions.

**Hierarchical Conflict Resolution** To address complex loci where a
single anchor overlaps multiple promoters (e.g., dense gene clusters or
bidirectional promoters), the function executes a 3-step resolution:

1.  *Biotype Prioritization:* Selects the highest-priority candidates by
    functional class:
    `Protein Coding > small-ncRNA (miRNA, snoRNA, snRNA, rRNA, scaRNA) > Antisense > lncRNA/ncRNA > Pseudogene`.

2.  *Expression Filter:* Within the selected biotype tier, excludes
    transcriptionally silent genes using a user-provided expression
    matrix. If no gene in the tier is expressed, all candidates in that
    tier are retained.

3.  *Expression Tiebreaker:* Among remaining candidates of equal biotype
    priority, retains all genes whose expression is within one order of
    magnitude of the highest-expressing candidate (i.e., expression \>=
    10\\

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
