# Internal: Reclassify Anchor by Expression

Given an anchor's gene symbol and type, filters to active genes (present
in `allow`) and optionally reclassifies silent promoters/enhancers.

## Usage

``` r
clean_anchor(g, t, allow, down)
```

## Arguments

- g:

  Character. Semicolon-delimited gene string.

- t:

  Character. Anchor type code (P, E, G, eP, eG).

- allow:

  Character vector. Whitelist of active gene symbols.

- down:

  Logical. If `TRUE`, reclassify silent P→eP and G→eG.

## Value

A list with `type` and `gene`.
