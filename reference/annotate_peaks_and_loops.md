# Spatial mapping of 1D genomic features to 3D chromatin interaction targets

A function for loop annotation and target mapping, designed to integrate
1D genomic features with 3D chromatin architecture.

1.  **Loop Annotation:** Classifies 3D spatial interactions (e.g.,
    distal-to-promoter, promoter-to-promoter) using positional anchor
    labels relative to gene annotations. **Important:** anchor type
    `"E"` denotes a *positional* distal/intergenic classification – it
    does **not** imply functional enhancer activity. Orthogonal
    chromatin data are required for functional interpretation.

2.  **Feature-to-Target Mapping:** Links 1D genomic features (e.g., GWAS
    risk SNPs, ATAC-seq peaks, ChIP-seq binding sites) to putative
    target genes via 3D chromatin contacts, providing a loop-based
    alternative to linear proximity-based assignments.

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
  co_dominance_ratio = 0.1,
  anchor_gap = -1L,
  anchor_min_overlap = 1L,
  anchor_min_frac = 0,
  write_output = TRUE,
  quiet = FALSE
)
```

## Arguments

- bedpe_file:

  Character. Path to a BEDPE file (at least 6 columns: chr1, start1,
  end1, chr2, start2, end2). Additional columns beyond 6 are retained in
  the output; anchor swapping only affects columns 1-6.

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
  One of `"hg38"`, `"hg19"`, `"mm10"`, `"mm9"`. Default: `"hg38"`. The
  package is currently tested on human and mouse; the architecture
  supports extension to other species by adding entries to
  `species_txdb_pkg()`, `species_orgdb_pkg()`, and
  `species_bsgenome_pkg()`.

- tss_region:

  Numeric vector of length 2. Promoter window around the TSS in bp.
  Default: `c(-2000, 2000)` (typical for mammalian protein-coding genes;
  may need widening for broad domains like HOX clusters, or narrowing
  for compact genomes).

- out_dir:

  Character. Output directory for the Excel results file. Default:
  `"./results"`.

- expr_matrix_file:

  Optional path to a normalised expression matrix (TPM/FPKM, genes x
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
  restricts to direct contacts. Default: `0`. Target gene assignment
  searches one additional hop (`neighbor_hop + 1`) to capture genes at
  the opposite anchor of directly connected loops.

- hub_percentile:

  Numeric (0-1). Loop-count quantile for hub detection. Default: `0.95`.
  Genes or distal elements with connectivity at or above this quantile
  are flagged as hubs. A minimum floor of 3 (promoter-centric) or 2
  (distal) is applied to avoid false hubs in small datasets.

- min_expr:

  Numeric. Minimum expression value for a gene to be considered active
  during anchor-level conflict resolution. Used only when
  `expr_matrix_file` is provided. Default: `0` (any non-zero expression,
  i.e. TPM \> 0). Increase to `1` or higher to require stronger
  evidence. Note: when `min_expr = 0`, the code uses a strict
  greater-than comparison (`> 0`) to exclude truly undetected genes;
  when `min_expr >= 1`, it uses `>= min_expr`. See
  [`refine_loop_anchors_by_expression`](https://zying106.github.io/looplook/reference/refine_loop_anchors_by_expression.md)
  for the separate `threshold` parameter that controls promoter
  reclassification.

- conflict_strategy:

  Character. Conflict resolution order for overlapping gene assignments.
  `"biotype_first"` (default): select the best biotype tier first, then
  apply expression filtering within that tier. `"expression_first"`:
  apply expression filtering across all biotypes first, then pick the
  best biotype among expressed candidates.

- co_dominance_ratio:

  Numeric (0-1). In the expression tiebreaker step, genes with
  expression \>= `co_dominance_ratio x max(expression)` in the group are
  retained together. Default: `0.1` (i.e. within one order of
  magnitude). Lower values (e.g. `0.01`) retain more co-dominant
  candidates; higher values (e.g. `0.5`) are more stringent.

- anchor_gap:

  Integer. Search radius: how far apart (bp) can a peak and loop anchor
  be for the peak to be considered "near" the anchor? `-1L` (default):
  strict physical overlap required (GenomicRanges default – peak and
  anchor must share at least 1 bp). `0L`: adjacent intervals (peak end
  == anchor start) also count. `>0`: explicit gap tolerance (e.g. `200`
  for cross-experiment integration). When `>= 0`, the result includes
  both physically overlapping pairs AND proximity-only pairs (within gap
  but no actual overlap). Use `anchor_min_overlap > 1` to require actual
  physical overlap among these candidates.

- anchor_min_overlap:

  Integer. After candidate pairs are found (via `anchor_gap`), how many
  base pairs of actual physical overlap are required? Default `1L`: any
  touch counts (including proximity-only hits when `anchor_gap >= 0`).
  Increase to `10-100` to filter out spurious boundary overlaps. Setting
  this `> 1` with `anchor_gap >= 0` ensures that even with gap
  tolerance, only pairs with genuine physical overlap are retained.

- anchor_min_frac:

  Numeric (0-1). After the first two filters, what fraction of the
  *peak* width must physically overlap the anchor? Default `0`: any
  fraction accepted. Set to `0.1-0.5` when peaks are broad (e.g. H3K27ac
  domains, 2-5 kb) so a 1 bp overlap does not link the entire broad
  peak. Ignored for point features (SNPs, eQTLs). Applied last, only to
  pairs that passed `anchor_gap` and `anchor_min_overlap`.

- write_output:

  Logical. If `TRUE` (default), write the Excel workbook to `out_dir`.
  If `FALSE`, return results without creating directories or files.

- quiet:

  Logical. If `TRUE`, suppress progress messages while preserving
  warnings. Default: `FALSE`.

## Value

An invisible named list:

- `target_annotation` – Target features (peaks) with gene assignments.
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

- `target_gene_links` – Long-format peak-gene provenance table. Each row
  records one peak-gene linkage with full provenance. **Read**
  `evidence`, `anchor_role`, and `gene_role` **together as a group** –
  they jointly describe how each gene was assigned to each peak; do not
  interpret any one column in isolation.

  - `input_id`, `loop_ID`, `anchor_id`: Identifiers.

  - `gene`: Linked gene symbol.

  - `gene_role`: `"promoter"`, `"gene_body"`, or `"linear_annotation"`.

  - `source`: `"loop_anchor"` (3D-derived) or `"linear_annotation"`
    (nearest gene).

  - `evidence`: Provenance label – `"local_promoter_overlap"` (peak
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

- `loop_annotation` – Annotated 3D interactome with
  `Putative_Target_Genes`.

- `anchor_loci_annotation` – Non-redundant anchor-locus genomic
  classifications after within-cluster interval reduction.

- `anchor_annotation` – Backward-compatible alias of
  `anchor_loci_annotation`.

- `promoter_centric_stats` – Gene-level connectivity statistics.

- `distal_element_stats` – Distal-element connectivity statistics.

- `plots` – Named list of ggplot/grob objects: `Basic_Donut`,
  `Basic_Circular`, `Basic_Flower`, `Karyo_LoopGenes`, `Karyo_Anchors`,
  `Anchor_Genomic_Distribution`, and (when `target_bed` is provided)
  `Karyo_TargetGenes`, `Target_Rose`, `Target_Genomic_Distribution`,
  `Target_Loop_Genomic_Distribution`.

- `plot_list` – Backward-compatible alias of `plots`.

- `metadata` – Internal metadata list (parameters, versions). Not
  intended for direct use.

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
    priority, retains all genes whose expression \>=
    `co_dominance_ratio` x the group maximum. Default `0.1` (one order
    of magnitude; ~3.3-fold in log2 space). This co-dominant rule
    preserves functionally redundant candidates such as bidirectional
    promoter pairs, where co-expressed partners typically fall within
    2-3 fold of each other. Edge case: when all candidates in the best
    biotype tier have TPM = 0, all are retained (the tiebreaker cannot
    distinguish them), which may produce multi-gene assignments at
    transcriptionally silent loci.

**Network Topology Analysis**

- **Ego-Network Expansion (`neighbor_hop`):** Implements k-hop
  neighborhood expansion via
  [`igraph::ego()`](https://r.igraph.org/reference/ego.html). A value of
  `0` restricts to direct contacts, while `1` includes secondary
  contacts to capture broader regulatory cliques.

- **Hub Detection:** Utilizes a node-degree quantile threshold
  (`hub_percentile`) to identify highly connected regulatory elements.

**Peak-Anchor Overlap Control** When `target_bed` is provided, three
parameters control how target peaks are matched to loop anchors. They
act as a cascade of increasingly stringent filters:

1.  `anchor_gap` – expands the search radius around each anchor.

2.  `anchor_min_overlap` – requires a minimum physical overlap in bp.

3.  `anchor_min_frac` – requires the overlap to cover a minimum fraction
    of the peak.

The table below lists suggested starting points for common experimental
designs.

|  |  |  |  |
|----|----|----|----|
| **Experimental design** | `anchor_gap` | `anchor_min_overlap` | `anchor_min_frac` |
| Same-experiment HiChIP / ChIA-PET (default) | `0` | `1` | `0` |
| Cross-experiment (e.g. ATAC-seq peaks x HiChIP loops) | `200-500` | `10` | `0` |
| Broad histone-mark peaks (H3K27ac, 2-5 kb) | `0` | `100` | `0.1` |
| Super-enhancers / wide domains (20-80 kb) | `0` | `500` | `0.05` |
| Point features (GWAS SNPs, eQTLs, 1 bp) | `0` | `1` | `0` |
| Stringent (high-confidence only) | `0` | `100` | `0.5` |

When `quiet = FALSE` (default), the function prints a diagnostic line
reporting how many peaks overlapped loop anchors after filtering,
helping you tune these thresholds for your data.

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
