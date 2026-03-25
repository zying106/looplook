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
#> GRanges object with 300 ranges and 0 metadata columns:
#>         seqnames              ranges strand
#>            <Rle>           <IRanges>  <Rle>
#>     [1]     chr1 171285865-171286512      *
#>     [2]     chr1   69582109-69582609      *
#>     [3]     chr1   93669531-93670158      *
#>     [4]     chr1 230112638-230115006      *
#>     [5]     chr1   83034931-83035876      *
#>     ...      ...                 ...    ...
#>   [296]     chr1     2193690-2196269      *
#>   [297]     chr1 101024799-101027134      *
#>   [298]     chr1   81834837-81835950      *
#>   [299]     chr1 222711087-222714102      *
#>   [300]     chr1 205449123-205449901      *
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
#> [1] 300
```
