# Internal: Collapse Delimited Gene String

Splits a semicolon-delimited gene string, removes duplicates and NAs,
and recollapses into a single string. Returns `NA_character_` if no
valid genes remain.

## Usage

``` r
extract_genes(genes_vec)
```

## Arguments

- genes_vec:

  Character vector of delimited gene strings.

## Value

A single semicolon-delimited string, or `NA_character_`.
