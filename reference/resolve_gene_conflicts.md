# Internal: Resolve Gene Conflicts via Biotype Priority Then Expression

For each genomic range, identifies all promoter-overlapping genes,
resolves conflicts using a two-stage strategy: (1) biotype priority
(protein-coding \> small-ncRNA \> antisense \> lncRNA/ncRNA \>
pseudogene), optionally overridden by `biotype_order`, and then (2)
expression-aware filtering within the selected biotype tier. If any gene
in the best tier is expressed (`tpm >= min_expr`), only expressed
candidates are retained; otherwise all candidates in that tier are kept.
When multiple candidates share the same biotype rank, a co-dominant
expression rule is applied: all genes with expression \>= 10\\ of the
group maximum are retained (collapsed with ";").

## Usage

``` r
resolve_gene_conflicts(
  current_anno_df,
  txdb_obj,
  org_db_pkg,
  tss_region,
  gene_expr_map,
  min_expr = 0,
  conflict_strategy = c("biotype_first", "expression_first"),
  co_dominance_ratio = 0.1,
  biotype_order = c("protein", "small_ncRNA", "antisense", "lncRNA", "pseudogene")
)
```

## Arguments

- current_anno_df:

  Data frame with columns suitable for
  [`makeGRangesFromDataFrame`](https://rdrr.io/pkg/GenomicRanges/man/makeGRangesFromDataFrame.html).

- txdb_obj:

  A `TxDb` object for gene coordinate lookup.

- org_db_pkg:

  Character. Organism database package name.

- tss_region:

  Numeric vector of length 2. TSS region for promoter definition, e.g.
  `c(-2000, 2000)`.

- gene_expr_map:

  Named numeric vector of per-gene expression values, or `NULL` if
  unavailable.

- min_expr:

  Numeric. Minimum expression value for a gene to be considered active
  during conflict resolution. Default: `0`.

- conflict_strategy:

  Character. Conflict resolution order. `"biotype_first"` (default):
  select the best biotype tier first, then apply expression filtering
  within that tier. This is the more conservative default – a silent
  protein-coding gene is preferred over a highly expressed lncRNA at the
  same locus. `"expression_first"`: apply expression filtering across
  all biotypes first, then pick the best biotype among expressed
  candidates.

- co_dominance_ratio:

  Numeric (0-1). In the expression tiebreaker step, genes with
  expression \>= `co_dominance_ratio * max(expression)` in the group are
  retained together. Default: `0.1` (i.e. within one order of
  magnitude). Lower values (e.g. `0.01`) retain more co-dominant
  candidates; higher values (e.g. `0.5`) are more stringent.

- biotype_order:

  Character vector. Custom ordering of biotype categories for conflict
  resolution. Each element must be one of five keywords (matched
  case-insensitively against the `GENETYPE` column): `"protein"`
  (protein-coding), `"small_ncRNA"` (miRNA, snoRNA, snRNA, rRNA,
  scaRNA), `"antisense"`, `"lncRNA"` (lncRNA and other ncRNA),
  `"pseudogene"`. The order determines priority: rank 1 = highest.
  Categories not listed keep their default relative order and are
  appended after the listed ones. Default:
  `c("protein", "small_ncRNA", "antisense", "lncRNA", "pseudogene")`. To
  prioritise lncRNAs over protein-coding genes while keeping everything
  else as-is, set `c("lncRNA", "protein")`.

## Value

The input data frame with `SYMBOL` and `annotation` columns resolved.

## Examples

``` r
# This example uses a sample TxDb included with GenomicFeatures.
# It demonstrates the core conflict-resolution logic on a synthetic
# genomic region near the hg19 HLA locus (chr6:29,940,000-29,950,000).
if (requireNamespace("org.Hs.eg.db", quietly = TRUE) &&
    requireNamespace("GenomicFeatures", quietly = TRUE)) {
    txdb <- AnnotationDbi::loadDb(
        system.file("extdata", "hg19_knownGene_sample.sqlite",
                     package = "GenomicFeatures")
    )

    # Minimal genomic region: seqnames, start, end
    test_region <- data.frame(
        seqnames = "chr6",
        start    = 29940000L,
        end      = 29950000L,
        stringsAsFactors = FALSE
    )

    # Run with default settings (biotype-first, no expression data)
    result <- resolve_gene_conflicts(
        current_anno_df    = test_region,
        txdb_obj           = txdb,
        org_db_pkg         = "org.Hs.eg.db",
        tss_region         = c(-2000, 2000),
        gene_expr_map      = NULL,
        conflict_strategy  = "biotype_first"
    )

    # Inspect resolved gene symbols (columns appear when genes are found)
    if ("SYMBOL" %in% colnames(result)) print(result$SYMBOL)
}
```
