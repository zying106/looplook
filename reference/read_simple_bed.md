# Read a Simple BED File into a GRanges Object

Reads the first three columns of a BED file (chrom, start, end) and
returns a
[`GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
object. Additional columns are ignored.

## Usage

``` r
read_simple_bed(bed_file)
```

## Arguments

- bed_file:

  Character. Path to a BED file (must be tab-delimited).

## Value

A [`GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
object, or `NULL` if `bed_file` is `NULL`.

## Examples

``` r
# 1. Locate the example BED file included in the package
bed_path <- system.file("extdata", "example_peaks.bed", package = "looplook")

# 2. Run the function (ensure file was found)
if (bed_path != "") {
  # Read BED file into a GRanges object
  gr <- read_simple_bed(bed_path)

  # 3. Inspect the result
  print(gr)

  # Check how many peaks were loaded
  length(gr)
}
#> GRanges object with 2000 ranges and 0 metadata columns:
#>          seqnames              ranges strand
#>             <Rle>           <IRanges>  <Rle>
#>      [1]     chr1   36118048-36118548      *
#>      [2]     chr1   64072284-64073108      *
#>      [3]     chr1 119630541-119633050      *
#>      [4]     chr1 209038325-209040095      *
#>      [5]     chr1   28064912-28068276      *
#>      ...      ...                 ...    ...
#>   [1996]     chr1 226777783-226779926      *
#>   [1997]     chr1   25642050-25642876      *
#>   [1998]     chr1   84602828-84605421      *
#>   [1999]     chr1     1437476-1439020      *
#>   [2000]     chr1   21320654-21327079      *
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
#> [1] 2000
```
