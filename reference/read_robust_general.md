# Robust Data Reader

Safely reads standard genomic formats utilizing
[`data.table::fread`](https://rdrr.io/pkg/data.table/man/fread.html) for
intelligent format inference.

## Usage

``` r
read_robust_general(
  f,
  header = FALSE,
  row_name = NULL,
  desc = "file",
  min_cols = 3
)
```

## Arguments

- f:

  Character. Path to the input file.

- header:

  Logical. Whether the file contains a header row.

- row_name:

  Integer or NULL. Column index to be used as row names.

- desc:

  Character. Short description for error logging.

- min_cols:

  Integer. Minimum number of columns required.

## Value

A data frame.
