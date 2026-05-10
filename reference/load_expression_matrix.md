# Internal: Load Expression Matrix

Reads a normalized expression matrix (TPM/FPKM), sets gene identifiers
as the first column, validates the requested sample columns, and returns
per-gene mean expression values. Sample column names must be unique;
missing or duplicated selections raise an error instead of being
silently dropped.

## Usage

``` r
load_expression_matrix(expr_matrix_file, sample_columns = NULL)
```

## Arguments

- expr_matrix_file:

  Character. Path to the expression matrix file.

- sample_columns:

  Character or integer vector. Sample columns to average. Character
  selections must exactly match unique column names.

## Value

Named numeric vector of per-gene mean expression values.
