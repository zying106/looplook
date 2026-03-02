# Perform GO Enrichment and Generate Network Plot

Runs Biological Process GO enrichment via `clusterProfiler` and
visualizes the top pathways and their core Hub genes in a highly
readable, non-overlapping divergent concept network (Cnetplot).

## Usage

``` r
run_go_enrichment(
  genes,
  org_db,
  universe_genes,
  cnet_nSample = 50,
  project_name = "Analysis",
  out_dir = "./"
)
```

## Arguments

- genes:

  Character vector of target gene symbols.

- org_db:

  Character. Organism database (e.g., "org.Hs.eg.db").

- universe_genes:

  Named numeric vector (LFC) used as background universe and color
  scale.

- cnet_nSample:

  Integer. Number of top GO pathways to display in the network.

- project_name:

  Character. Prefix for outputs.

- out_dir:

  Character. Output directory.

## Value

Data frame of raw GO enrichment results. Saves a Cnetplot PDF.
