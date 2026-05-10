# Internal: Simplify Genomic Annotation to Broad Categories

Collapses detailed ChIPseeker annotation strings into five broad
categories: Promoter, Intron, Exon, Distal Intergenic, and Downstream.
Anything unrecognised is labelled "Others".

## Usage

``` r
simplify_annotation(x)
```

## Arguments

- x:

  Character vector of annotation strings.

## Value

Character vector of simplified categories.
