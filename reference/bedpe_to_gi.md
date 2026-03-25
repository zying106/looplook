# Read BEDPE File into a GInteractions Object

Reads a standard BEDPE file and converts it into a Bioconductor
[`GInteractions`](https://rdrr.io/pkg/InteractionSet/man/GInteractions-class.html)
object.

## Usage

``` r
bedpe_to_gi(bedpe_file)
```

## Arguments

- bedpe_file:

  Character. Path to a BEDPE file. Must contain at least six columns:
  `chr1, start1, end1, chr2, start2, end2`.

## Value

A
[`GInteractions`](https://rdrr.io/pkg/InteractionSet/man/GInteractions-class.html)
object with a `score` metadata column (defaulting to 0 if not provided).

## Details

**Anchor Normalization:** Anchor order is automatically normalized so
that the first anchor is lexicographically less than or equal to the
second (e.g., chr1 \< chr2), ensuring compatibility with
`GInteractions(mode = "strict")`.

**Score Detection:** The function attempts to automatically detect a
numeric score column. It checks the 8th column first (standard for many
tools); if not numeric, it falls back to the 7th column. Non-numeric
values are treated as 0.

## Examples

``` r
# 1. Locate the example BEDPE file included in the package
# system.file finds the absolute path to 'inst/extdata/example_loops_1.bedpe'
bedpe_path <- system.file("extdata", "example_loops_1.bedpe", package = "looplook")

# 2. Run the function (ensure file was found)
if (bedpe_path != "") {
  gi <- bedpe_to_gi(bedpe_path)

  # 3. Inspect the result
  print(gi)

  # Check the imported score column
  S4Vectors::mcols(gi)$score
}
#> StrictGInteractions object with 300 interactions and 1 metadata column:
#>         seqnames1             ranges1     seqnames2             ranges2 |
#>             <Rle>           <IRanges>         <Rle>           <IRanges> |
#>     [1]      chr1 109492652-109495835 ---      chr1 110644072-110649627 |
#>     [2]      chr1 116687216-116689300 ---      chr1 116844969-116856809 |
#>     [3]      chr1 116876465-116881376 ---      chr1 117819956-117830132 |
#>     [4]      chr1 146374138-146379497 ---      chr1 146442737-146445345 |
#>     [5]      chr1     1370034-1377149 ---      chr1     1469037-1474411 |
#>     ...       ...                 ... ...       ...                 ... .
#>   [296]      chr1     1011273-1015982 ---      chr1     1018256-1022868 |
#>   [297]      chr1     1342739-1364987 ---      chr1     1397119-1403614 |
#>   [298]      chr1 147271364-147274896 ---      chr1 147290756-147294529 |
#>   [299]      chr1 108639596-108642391 ---      chr1 109098569-109102247 |
#>   [300]      chr1   11011474-11014295 ---      chr1   12745038-12747130 |
#>             score
#>         <numeric>
#>     [1]         1
#>     [2]         1
#>     [3]         2
#>     [4]         1
#>     [5]         3
#>     ...       ...
#>   [296]        12
#>   [297]         7
#>   [298]         4
#>   [299]         1
#>   [300]         1
#>   -------
#>   regions: 385 ranges and 0 metadata columns
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
#>   [1]   1   1   2   1   3   1   2   8   3   1   1   4   2   1   1   1   1   1
#>  [19]   2   8   3   1   1   1   1 219   7   2   1   1   2   1   1   1   1   1
#>  [37]   1   3   1   1   2  44   1   1   2   4   1   1   1   3   1   1   5   2
#>  [55]   1   1   3   2   1   1   3   1   1   1   1   1   1   5   1   2   1   1
#>  [73]   1   2   1   2   5   1   2   1   1   1   1   1   2   1   1   1   8   1
#>  [91]   3   1   1   2   3   1   1   2   2   1   1   1   1   1   1   1   1   1
#> [109]  11   7   1   9   1   1   1   1   1   1   1   2   1   3   3   1   2   1
#> [127]   2   1   1   2   1   1   1   1   4   1   1   1   1   1   3   1   3   2
#> [145]   7   1   1   1   2   3   1   2   2   1   2   9   1   2   5   1   3   4
#> [163]   1   8   1   1   1   3   1   1   1   7   7   1   2   1   2   1   5   1
#> [181]   1   1   1   1   1   2   5   1   1   1  12   2   1   1   1   1   2  14
#> [199]   2   1   1   1  33   1   1   1   1   4   1   1   1   2   1  36   1   3
#> [217]   2   1   1   2   2   1   1   1   1   5  10   1   1   1   1   1   3  36
#> [235]   1   2   1   3   1   1   2   1   1   3   1   1   2   5   1   1   1   3
#> [253]   1   9   2   1   1   1   1   4   1   1   9   1   2   4   1   2   1   1
#> [271]   1  19   1  10   1   2   2   1  10   1   1   2   1   1   3   1   1   1
#> [289]   1   4   1   1   2   1   2  12   7   4   1   1
```
