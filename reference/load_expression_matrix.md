# Internal: Load Expression Matrix

Reads a normalized expression matrix (TPM/FPKM), sets gene identifiers
as row names, extracts the specified sample columns, and returns
per-gene mean expression values.

## Usage

``` r
load_expression_matrix(expr_matrix_file, sample_columns = NULL)
```

## Arguments

- expr_matrix_file:

  Character. Path to the expression matrix file.

- sample_columns:

  Character or integer vector. Sample columns to average.

## Value

Named numeric vector of per-gene mean expression values.
