# Internal: Compute Refined Network Statistics

Recalculates promoter-centric and distal-element connectivity statistics
after expression-aware filtering, merging with upstream annotation stats
where available.

## Usage

``` r
compute_refined_stats(
  loop_df,
  upstream_promoter_stats,
  vals,
  threshold,
  hub_percentile
)
```

## Arguments

- loop_df:

  Loop annotation data frame after expression filtering.

- upstream_promoter_stats:

  Upstream promoter-centric stats (or NULL).

- vals:

  Named numeric vector of per-gene mean expression.

- threshold:

  Numeric. Expression threshold for active gene classification.

- hub_percentile:

  Numeric. Quantile for hub cutoff.

## Value

A list with `promoter_centric` and `distal_element` data frames.
