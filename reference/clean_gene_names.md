# Internal: Clean Gene Name Vector

Removes empty strings, NA values, and duplicate entries from gene
identifiers. Optionally splits concatenated strings (e.g., "TP53;BRCA1")
before cleaning.

## Usage

``` r
clean_gene_names(x, split = NULL)
```

## Arguments

- x:

  Character vector of gene names, possibly containing delimiters.

- split:

  Character. If non-NULL, a regex passed to
  [`strsplit`](https://rdrr.io/r/base/strsplit.html) to split
  concatenated gene strings (e.g., `"[;,]"`). Set to `NULL` if `x` is
  already a clean character vector.

## Value

A unique, non-empty, non-NA character vector.
