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
  anchor_merge_gap = 0,
  out_dir = "./results",
  expr_matrix_file = NULL,
  sample_columns = NULL,
  project_name = "HiChIP",
  color_palette = "Set2",
  karyo_bin_size = 1e+05,
  neighbor_hop = 0,
  hub_percentile = 0.95,
  hub_metric = c("unique_contacts", "total_loops"),
  min_expr = 0,
  conflict_strategy = c("biotype_first", "expression_first"),
  co_dominance_ratio = 0.1,
  biotype_order = c("protein", "small_ncRNA", "antisense", "lncRNA", "pseudogene"),
  anchor_gap = -1L,
  anchor_min_overlap = 1L,
  anchor_min_frac = 0,
  write_output = TRUE,
  quiet = FALSE,
  target_priority = c("promoter_then_distance", "distance_then_role")
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

  Character. Genome assembly string (e.g. `"hg38"`, `"mm10"`,
  `"danRer11"`, `"dm6"`). Default: `"hg38"`. When `txdb` and `org_db`
  are `NULL`, auto-resolved from built-in species (hg38/hg19/mm10/mm9);
  for any other species you must pass `txdb` and `org_db` as objects or
  package name strings directly (e.g.
  `txdb = TxDb.Dmelanogaster.UCSC.dm6.ensGene`,
  `org_db = "org.Dm.eg.db"`). The `species` string is also used for
  karyotype ideograms and JASPAR motif species filtering.

- tss_region:

  Numeric vector of length 2. Promoter window around the TSS in bp.
  Default: `c(-2000, 2000)` (+/-2 kb; typical for mammalian
  protein-coding genes; may need widening for broad domains like HOX
  clusters, or narrowing for compact genomes).

- anchor_merge_gap:

  Integer. Merge loop anchors within this many bp on the same chromosome
  before building the connectivity graph. Default `0` (no merging;
  anchors must match exactly on chr, start, end). Set to `50--200` when
  input BEDPE comes directly from a loop caller without prior replicate
  consolidation – the same biological anchor may appear with slightly
  different boundaries in different loop rows, which would otherwise
  fragment the graph. Not needed when input is pre-consolidated via
  [`consolidate_chromatin_loops`](https://zying106.github.io/looplook/reference/consolidate_chromatin_loops.md).

- out_dir:

  Character. Output directory for the Excel results file. Default:
  `"./results"`.

- expr_matrix_file:

  Optional path to a normalised expression matrix (genes x samples).
  Accepts steady-state RNA-seq (TPM/FPKM), nascent transcription data
  (NET-seq, PRO-seq, GRO-seq, TT-seq), or CAGE-seq. Enables
  expression-aware conflict resolution. Default: `NULL`.

- sample_columns:

  Character vector or integer indices. Columns in `expr_matrix_file` to
  average for baseline expression. Default: `NULL`.

- project_name:

  Character. Prefix for output files and plot titles. Default:
  `"HiChIP"`.

- color_palette:

  Character. RColorBrewer palette name. Default: `"Set2"`.

- karyo_bin_size:

  Integer. Bin width in base pairs (bp) for karyotype heatmaps. Default:
  `1e5` (100 kb). Typical range: 5e4-5e5 depending on genome size and
  resolution.

- neighbor_hop:

  Integer. k-hop ego-network expansion order via
  [`igraph::ego()`](https://r.igraph.org/reference/ego.html). `0`
  (default) restricts to direct loop contacts. `1` additionally includes
  2-hop expanded targets for exploratory network analysis. Values
  greater than `1` are not supported. Target gene assignment searches
  one additional hop (`neighbor_hop + 1`) to capture genes at the
  opposite anchor of directly connected loops.

- hub_percentile:

  Numeric (0-1). Loop-count quantile for hub detection. Default: `0.95`.
  Genes or distal elements with connectivity at or above this quantile
  are flagged as hubs. A minimum floor of 3 (promoter-centric) or 2
  (distal) is applied to avoid false hubs in small datasets.

- hub_metric:

  Character. Which connectivity count to use for hub detection.
  `"unique_contacts"` (default): counts distinct neighbour anchor IDs,
  robust to duplicate/replicate loop rows. `"total_loops"`: counts all
  loop rows (backward compatible; may inflate hub calls when input
  contains unconsolidated replicates).

- min_expr:

  Numeric. Minimum expression value for a gene to be considered active
  during anchor-level conflict resolution. Used only when
  `expr_matrix_file` is provided. Default: `0` (any non-zero
  expression). Increase to `1` or higher to require stronger evidence.
  Note: when `min_expr = 0`, the code uses a strict greater-than
  comparison (`> 0`) to exclude truly undetected genes; when
  `min_expr >= 1`, it uses `>= min_expr`. For nascent transcription data
  (NET-seq, PRO-seq), gene-body aggregated signal should be used as the
  quantitative input. See
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

- biotype_order:

  Character vector. Custom ordering of biotype categories for gene
  conflict resolution (passed to
  [`resolve_gene_conflicts`](https://zying106.github.io/looplook/reference/resolve_gene_conflicts.md)).
  Five keywords: `"protein"`, `"small_ncRNA"`, `"antisense"`,
  `"lncRNA"`, `"pseudogene"`. Listed categories get top priority;
  unlisted keep their default relative order appended at the bottom.
  Default:
  `c("protein", "small_ncRNA", "antisense", "lncRNA", "pseudogene")`.

- anchor_gap:

  Integer. Candidate search radius in bp between a peak and a loop
  anchor. `-1L` (default): strict physical overlap required. When
  `>= 0`, expands the candidate search window (e.g. `200` for
  cross-experiment coordinate shifts). **Note:** final retention always
  requires at least `anchor_min_overlap` bp of physical overlap; this
  parameter only controls which candidates are evaluated, not whether
  proximity-only pairs are retained.

- anchor_min_overlap:

  Integer. Minimum required physical overlap (bp) between a peak and an
  anchor. Default `1L`: at least 1 bp of actual sequence overlap
  required – proximity-only pairs within the `anchor_gap` window but
  without physical overlap are excluded. Increase to `10-100` for broad
  peaks to avoid spurious boundary overlaps.

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

- target_priority:

  Character. How to prioritise among multiple candidate target genes per
  input feature. The policy applies only within primary target links
  (default `path_length <= 1`); longer paths are reported separately as
  `Expanded_Target_Genes` and do not participate in
  `Assigned_Target_Genes` selection. `"promoter_then_distance"`
  (default): within primary links, promoter evidence dominates – all
  promoter-linked genes beat all gene-body genes regardless of path
  length; within each tier shorter paths win. Exception: direct strict
  promoter–promoter contacts (same loop, path0 + path1) co-assign both
  endpoints via union, because the technical endpoint orientation does
  not reflect biological regulatory direction. `"distance_then_role"`:
  within primary links, path-length dominates – the closest linked gene
  wins; at equal distance promoter beats gene-body (legacy behaviour).
  The policy affects `Assigned_Target_Genes` only.
  `Regulated_promoter_genes` always reports all promoter-linked genes
  regardless of the chosen policy.

## Value

A named list:

- `target_annotation` – Target features (peaks) with gene assignments.
  Key columns include:

  - `All_Loop_Connected_Genes`: Inclusive provenance union of all
    loop-anchor gene links. May include strict assignment-eligible
    targets and non-strict positional/enhancer candidates. Not a
    confirmed target-gene set.

  - `Regulated_promoter_genes`: Promoter genes supported by loop-anchor
    context.

  - `Assigned_Target_Genes`: Policy-based 3D assignment (default:
    promoter-first, then shorter path wins; see `target_priority`).

  - `*_Filled` variants: Linear nearest-gene fallback when strict 3D
    assignments are empty.

  - `Regulated_promoter_Evidence`: Provenance of
    `Regulated_promoter_genes` (e.g., `local_promoter_overlap`,
    `distal_promoter`). **Read with** `Regulated_promoter_genes`; do not
    cross-reference with `Assigned_Target_Genes` or other columns.

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

  - `gene_role`: `"promoter"`, `"gene_body"`, `"enhancer_candidate"`,
    `"positional_candidate"`, or `"linear_annotation"`.

  - `source`: `"loop_anchor"` (3D-derived) or `"linear_annotation"`
    (nearest gene).

  - `evidence`: Provenance label – `"local_promoter_overlap"` (peak
    overlaps anchor promoter), `"distal_promoter"` (promoter on the
    distal loop anchor), `"gene_body_context"` /
    `"distal_gene_body_context"` (gene body linkage),
    `"local_enhancer_candidate"` / `"distal_enhancer_candidate"` /
    `"expanded_enhancer_candidate"` (enhancer-associated linkage),
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
  `Putative_Target_Genes` (all P/G anchor genes connected through the
  loop network) and `Promoter_Target_Genes` (promoter-only subset,
  P-side genes). G–P/P–G asymmetry is respected. Does not include linear
  nearest-gene fallback.

- `anchor_loci_annotation` – Non-redundant anchor-locus genomic
  classifications after within-cluster interval reduction.

- `anchor_annotation` – Backward-compatible alias of
  `anchor_loci_annotation`.

- `promoter_centric_stats` – Gene-level connectivity statistics.

- `distal_element_stats` – Distal anchor connectivity (E, dual, G, eG).

- `plots` – Named list of ggplot/grob objects: `Basic_Donut`,
  `Basic_Circular`, `Basic_Flower`, `Karyo_LoopGenes`, `Karyo_Anchors`,
  `Anchor_Genomic_Distribution`, and (when `target_bed` is provided)
  `Karyo_TargetGenes`, `Target_Rose`, `Target_Genomic_Distribution`,
  `Target_Loop_Genomic_Distribution`.

- `plot_list` – Backward-compatible alias of `plots`.

- `metadata` – Internal metadata list (parameters, versions). Not
  intended for direct use.

The returned object carries a `looplook_anchor_state` attribute (access
via `attr(result, "looplook_anchor_state")`) containing internal anchor
topology data required by
[`refine_loop_anchors_by_chromatin`](https://zying106.github.io/looplook/reference/refine_loop_anchors_by_chromatin.md)
with `recompute_targets = TRUE`. If `write_output = TRUE`, also writes a
multi-sheet Excel workbook to `out_dir`.

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
| Same-experiment HiChIP / ChIA-PET (default) | `-1` | `1` | `0` |
| Cross-experiment (e.g. ATAC-seq peaks x HiChIP loops) | `200-500` | `10` | `0` |
| Broad histone-mark peaks (H3K27ac, 2-5 kb) | `0` | `100` | `0.1` |
| Super-enhancers / wide domains (20-80 kb) | `0` | `500` | `0.05` |
| Point features (GWAS SNPs, eQTLs, 1 bp) | `0` | `1` | `0` |
| Stringent (high-confidence only) | `0` | `100` | `0.5` |

When `quiet = FALSE` (default), the function prints a diagnostic line
reporting how many peaks overlapped loop anchors after filtering,
helping you tune these thresholds for your data.

## See also

[`refine_loop_anchors_by_expression`](https://zying106.github.io/looplook/reference/refine_loop_anchors_by_expression.md)
for expression-aware refinement,
[`refine_loop_anchors_by_chromatin`](https://zying106.github.io/looplook/reference/refine_loop_anchors_by_chromatin.md)
for chromatin-aware reclassification,
[`profile_target_genes`](https://zying106.github.io/looplook/reference/profile_target_genes.md)
for automated functional profiling.

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
