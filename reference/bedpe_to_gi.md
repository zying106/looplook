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
#> StrictGInteractions object with 1000 interactions and 1 metadata column:
#>          seqnames1             ranges1     seqnames2             ranges2 |
#>              <Rle>           <IRanges>         <Rle>           <IRanges> |
#>      [1]      chr1 112721687-112724586 ---      chr1 112791557-112793600 |
#>      [2]      chr1 113867048-113875597 ---      chr1 113927096-113933645 |
#>      [3]      chr1 115060199-115071191 ---      chr1 116162680-116169369 |
#>      [4]      chr1 118807309-118811756 ---      chr1 119277565-119282472 |
#>      [5]      chr1 112399375-112406783 ---      chr1 113756152-113762298 |
#>      ...       ...                 ... ...       ...                 ... .
#>    [996]      chr1   12374800-12386079 ---      chr1   12848181-12853963 |
#>    [997]      chr1 109654027-109657363 ---      chr1 110405568-110409319 |
#>    [998]      chr1 115176997-115186579 ---      chr1 115779639-115785509 |
#>    [999]      chr1   12374800-12386079 ---      chr1   12548119-12551309 |
#>   [1000]      chr1   10054549-10058953 ---      chr1   11801550-11808331 |
#>              score
#>          <numeric>
#>      [1]         1
#>      [2]         4
#>      [3]         2
#>      [4]         1
#>      [5]         2
#>      ...       ...
#>    [996]         1
#>    [997]         1
#>    [998]         1
#>    [999]         2
#>   [1000]         1
#>   -------
#>   regions: 689 ranges and 0 metadata columns
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
#>    [1]   1   4   2   1   2   1   1   1   2   1   1   2   1   3   1   2   1   7
#>   [19]   1   1   1   2   1   2   1   1   1   1   6   1   4   1   1   1   1   2
#>   [37]   2   2   1   1   7   4   2   1   1   3   1   8   1   1   2   1   1   1
#>   [55]   2   4   2   1   7   1   5   3   1   3   4   1   1   3   1   1   1   1
#>   [73]   1   1   1   2   1   1   1   1   1   1   1   1   2   4   2   2   5   4
#>   [91]   1   1   1   1   6   1   1   2   1   1   2   3   1   3   2   4   3   1
#>  [109]   1   1   1   2   2   1   1   2   1   1   2   1   1   1   1   1   1   4
#>  [127]   1   1   1   1  16   1   4   1   1   1   1   1   1   1   1   4   1   1
#>  [145]   3   2   1   3   1   1   1   6   1   7   2   1   1   1   1   4   1   2
#>  [163]   1   1   1   2   1  12   3   1   4  12   2   1  14   1   1   1   1   2
#>  [181]   1   2  10   4   1   1   1   1   1   1   1   2   1   1   6   1   1  44
#>  [199]   1   2   1   1   1   2   1   1   2   2   1   4   2   1   1   3   1   1
#>  [217]   1   9   1   4   1   1   1   1   1   5   2   1   1   2   1   1   8   1
#>  [235]   2   1   1   1   2   1   1   1   1  13   1   1   4   1   3   2   5   1
#>  [253]   1  14  22   1   1   1   2   1   2   2   1   3   4   7   1   1   2   3
#>  [271]   1   1   1   4   1   1   4   1   6   1  10   1   1   1   1   1   2   1
#>  [289]   1   1   1   1   1   5   2   1   1   1   1   3   6   1   1   1   1   1
#>  [307]   1  32   1   1   1   1   1   1   1   1   1   1   1   1   2   1   4   5
#>  [325]   3   1   2   2   1   1   1   2   1   1   1   5   1   1   1   3   1   1
#>  [343]   1  23   1   1   1   2   1   1   1   1   3   1   1   1   1   1   7   2
#>  [361]   1   1   1   4   1   1   1   1   1   1   2   1   2   2   1   1   2   1
#>  [379]   1   1   1   1   1   1   1   5   1   1   1   1  14   1   1   1   1   1
#>  [397]   1   2   1   1   1   1   1   1   1  10   1  21   1   1   5   1   1   1
#>  [415]   1   1   1   1   1   1   1   2   1   1   1   1   3   9   2   1   4   2
#>  [433]   2   1   1  19   1   1   1   8   1   1   1   1   1   3   1   1   1   2
#>  [451]  70   6   2   1   5   2   1   2   1   1   1   1   2   1   2   2   1   1
#>  [469]   2   1   2   1   2   1   1   1   1  36   1   1   1   3   7   1   1   2
#>  [487]   2   1   1   2   2   1   2   1   1   1   1  23   5   1   1   2   2   1
#>  [505]   1   1   1   1  44   1   3   1   1  12   4   1   4   1   1   3   1  10
#>  [523]   5   3   1   1   6   5   1   1   1   3   1   2   7   2   1   2   1   1
#>  [541]   1   2   1   1   1   1   1   2   1   1   1   1   1 219   9   1   1   3
#>  [559]  21   1   1   1   1 244   1   3   1   1   1   1   1   2   1   1   2   1
#>  [577]   1   2   1   3   1   1   1   1   3   1   1   4   1   6   2   1   1   2
#>  [595]   3   1   1   1   5   1   1   1   1   1   1   4   7   1   1   3   3   6
#>  [613]   1   2   1   2   1   1   2   1   1   1   1   2   1   3   1   1   1   1
#>  [631]   1  12   1   1   8   1   1   1   2   2   1   1   2   1   1   1   1   1
#>  [649]   2   1   2   9   1   8   1   4   2   2  17   2   1   1   1   1  10   1
#>  [667]   1   1   1   1   1   1   1   1   8   1   3   1   1   1   1   1   1   2
#>  [685]   1   1   1   3   1   2   2   1   2   2   1   1   1   2   1   2   1   2
#>  [703]   2   2   2   1   1   1   1   1   1   1   1   1   1   1   4   1   2   2
#>  [721]   1   1   3   1   3   2   1   1  11   1   8   3   2   1   1   1   1   1
#>  [739]   2  13   1   1   3   2   2   3   1   1   3   2   1   1   1   1   1   3
#>  [757]  13   1   1   3  12   1   1   3   1   1   1   1   1   2   1   2   2   2
#>  [775]   1   1   1   4   1   1   1   1   1   1   3   2   1   2   8   1   1   1
#>  [793]   1   1   3   5   2   7   2   1   1   5   2   2   1   1   2   1   2   1
#>  [811]   1   3   1   4   1   1   1   1   1   1   1   1   2   1   1   1   2   2
#>  [829]   1   1  12  12   1   1   3   1   1   8   1   1   2   2   1   3   5   1
#>  [847]   1   2   2   1   1   2   1   1   1   2   3   1   1   2   9   3   1  28
#>  [865]   1   1   1   2   2   3   1   1   1   1   2   2   1   3   1   1   1   1
#>  [883]  10   4   1   5  17   2   1   1   4   1   1   5   1   3   1  10   1   2
#>  [901]   2   1   2   4   7   1   8   1   1   1   2   2   4   1   1   2   1   1
#>  [919]   1  12   1   2   1   1  12   1   3  33   1   1   1   1   4   2   2   1
#>  [937]   9   1   1   1   1   1   1   3   1   4   2   1   2   1   3   1   1   1
#>  [955]   2   1   1   1   2   1   1   2   1   1   3   3   7   1   1   1   1  36
#>  [973]   1   1   1   4   2   2   1   1   2   1   9   1   1   1   1   1   1   5
#>  [991]   2   2   1   1   4   1   1   1   2   1
```
