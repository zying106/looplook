# Internal: Resolve Gene Conflicts via Expression & Biotype

For each genomic range in an annotation data frame, identifies all genes
whose promoters overlap the range, then selects the best candidate using
expression level and biotype priority (protein-coding \> antisense \>
lncRNA \> pseudo). When multiple genes have similar expression, all are
retained (collapsed with ";").

## Usage

``` r
resolve_gene_conflicts(
  current_anno_df,
  txdb_obj,
  org_db_pkg,
  tss_region,
  gene_expr_map
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

## Value

The input data frame with `SYMBOL` and `annotation` columns resolved.
