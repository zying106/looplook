# Internal: Resolve Gene Conflicts via Biotype Priority Then Expression

For each genomic range, identifies all promoter-overlapping genes,
resolves conflicts using a two-stage strategy: (1) biotype priority
(protein-coding \> small-ncRNA \> antisense \> lncRNA/ncRNA \>
pseudogene), then (2) expression-aware filtering within the selected
biotype tier. If any gene in the best tier is expressed
(`tpm >= min_expr`), only expressed candidates are retained; otherwise
all candidates in that tier are kept. When multiple candidates share the
same biotype rank, a co-dominant expression rule is applied: all genes
with expression \>= 10\\ of the group maximum are retained (collapsed
with ";").

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
  co_dominance_ratio = 0.1
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

## Value

The input data frame with `SYMBOL` and `annotation` columns resolved.
