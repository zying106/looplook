# Print looplook karyogram

Displays a previously captured karyotype heatmap. Uses embedded PNG
bytes when available, otherwise falls back to a file-backed image.

## Usage

``` r
# S3 method for class 'looplook_karyo'
print(x, ...)
```

## Arguments

- x:

  A `looplook_karyo` object.

- ...:

  Additional arguments (unused).

## Value

The input `looplook_karyo` object, returned invisibly after drawing the
image.
