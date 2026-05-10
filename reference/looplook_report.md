# Render a Publication-Ready 3D Annotation Report

One-click parameterised R Markdown report that executes the full
looplook pipeline (annotation → refinement → profiling) and renders an
interpretation-ready HTML document suitable for sharing with
collaborators.

## Usage

``` r
looplook_report(
  bedpe_file,
  target_bed = NULL,
  expr_matrix_file = NULL,
  sample_columns = NULL,
  species = "hg38",
  project_name = "looplook Analysis",
  out_dir = "looplook_results",
  threshold = 1,
  unit_type = "TPM",
  reclassify_by_expression = TRUE,
  tss_region = c(-2000, 2000),
  neighbor_hop = 0,
  hub_percentile = 0.95,
  color_palette = "Set2",
  target_mapping_mode = "all",
  include_Filled = TRUE,
  use_nearest_gene = FALSE,
  target_source = "targets",
  loop_types = c("E-P", "P-P"),
  stat_test = "wilcox.test",
  run_go = TRUE,
  run_ppi = FALSE,
  run_motif = FALSE,
  genome_id = species,
  motif_p_thresh = 1e-04,
  motif_ntop = 5,
  ppi_score = 400,
  ppi_nSample = 400,
  heatmap_nSample = 99999,
  gsea_nSample = 99999,
  cnet_nSample = 50,
  karyo_bin_size = 1e+05,
  diff_file = NULL,
  lfc_col = "log2FoldChange",
  metadata_file = NULL,
  precomputed_res = NULL,
  output_file = NULL,
  quiet = FALSE,
  ...
)
```

## Arguments

- bedpe_file:

  Path to a BEDPE file of chromatin loops.

- target_bed:

  Optional path to a BED file of genomic features.

- expr_matrix_file:

  Optional path to a normalised expression matrix.

- sample_columns:

  Sample columns in the expression matrix to average.

- species:

  Genome assembly (`"hg38"`, `"hg19"`, `"mm10"`, `"mm9"`).

- project_name:

  Character. Project prefix for the report title.

- out_dir:

  Output directory. Created if missing.

- threshold:

  Numeric. Expression threshold for active gene classification.

- unit_type:

  Character. Expression unit label for plot annotations. Default
  `"TPM"`.

- reclassify_by_expression:

  Logical. Reclassify silent promoters as eP/eG.

- tss_region:

  Numeric vector of length 2. TSS flanking region in bp. Default
  `c(-2000, 2000)`.

- neighbor_hop:

  Integer. Neighbourhood extension for fragment overlap. Default `0`.

- hub_percentile:

  Numeric. Quantile threshold for hub classification. Default `0.95`.

- color_palette:

  Character. RColorBrewer qualitative palette name. Default `"Set2"`.

- target_mapping_mode:

  Character. Mapping strategy for target genes. Default `"all"`.

- include_Filled:

  Logical. Include comprehensively merged gene assignments. Default
  `TRUE`.

- use_nearest_gene:

  Logical. Bypass 3D loop-based assignment and use linear proximity.
  Default `FALSE`.

- target_source:

  Character vector. Source of target genes to profile. Default
  `"targets"`.

- loop_types:

  Character vector. Loop types to include in profiling. Default
  `c("E-P", "P-P")`.

- stat_test:

  Character. Statistical test for violin comparisons. Default
  `"wilcox.test"`.

- run_go:

  Logical. Run GO enrichment (requires clusterProfiler).

- run_ppi:

  Logical. Run protein-protein interaction network analysis. Default
  `FALSE`.

- run_motif:

  Logical. Run transcription factor motif analysis. Default `FALSE`.

- genome_id:

  Character. Reference genome for motif scanning. Defaults to `species`.

- motif_p_thresh:

  Numeric. P-value threshold for motif enrichment. Default `1e-4`.

- motif_ntop:

  Numeric. Number of top motifs to display. Default `5`.

- ppi_score:

  Numeric. Minimum STRING combined score. Default `400`.

- ppi_nSample:

  Numeric. Maximum genes to include in PPI. Default `400`.

- heatmap_nSample:

  Numeric. Maximum genes in expression heatmap. Default `99999`.

- gsea_nSample:

  Numeric. Maximum genes sampled for GSEA. Default `99999`.

- cnet_nSample:

  Numeric. Number of GO terms in cnetplot. Default `50`.

- karyo_bin_size:

  Numeric. Bin size for karyotype heatmaps. Default `100000`.

- diff_file:

  Optional differential expression result file.

- lfc_col:

  Column name for log2 fold change in `diff_file`.

- metadata_file:

  Optional sample metadata file.

- precomputed_res:

  Optional. Either a `.RData` file path or an in-memory list object
  returned by `annotate_peaks_and_loops`. When provided, annotation is
  skipped and refinement starts from this object.

- output_file:

  Character. Output HTML file name. `NULL` derives from `project_name`.

- quiet:

  Logical. Suppress rendering output. Default `FALSE`.

- ...:

  Additional arguments passed to
  [`rmarkdown::render`](https://rdrr.io/pkg/rmarkdown/man/render.html).

## Value

The path to the generated HTML report (invisibly).

## Details

The profiling stage uses R's random number generator for GSEA gene-set
down-sampling (`gsea_nSample`) and motif background sampling. Call
[`set.seed()`](https://rdrr.io/r/base/Random.html) before
`looplook_report()` for fully reproducible results.

## Examples

``` r
# \donttest{
looplook_report(
  bedpe_file   = system.file("extdata", "example_loops_1.bedpe", package = "looplook"),
  target_bed   = system.file("extdata", "example_peaks.bed", package = "looplook"),
  project_name = "My Study",
  out_dir      = tempdir()
)
#> >> Generating looplook report: My_Study_Report.html
#> 
#> 
#> processing file: skeleton.Rmd
#> 1/36                            
#> 2/36 [setup]                    
#> 3/36                            
#> 4/36 [summary-block]            
#> 5/36                            
#> 6/36 [annotation-core]          
#> Error in .resolve_db(txdb, tx_species, "TxDb"): TxDb 'TxDb.Hsapiens.UCSC.hg38.knownGene' not installed
#> Error: object 'res' not found
#> Error: object 'res' not found
#> Error: object 'ta' not found
#> Error: object 'ta' not found
#> Error: object 'n_targets' not found
#> Error: object 'la' not found
#> Error: object 'n_loop_types' not found
#> Error: object 'n_targets' not found
#> Error: object 'la' not found
#> Error: object 'n_targets' not found
#> 7/36                            
#> 8/36 [annotation-donut]         
#> Error: object 'res' not found
#> 9/36                            
#> 10/36 [annotation-table-top]     
#> Error: object 'ta' not found
#> 11/36                            
#> 12/36 [annotation-table-loops]   
#> Error: object 'res' not found
#> Error: object 'la' not found
#> 13/36                            
#> 14/36 [annotation-table-promoter]
#> Error: object 'res' not found
#> Error: object 'pc' not found
#> 15/36                            
#> 16/36 [annotation-table-distal]  
#> Error: object 'res' not found
#> Error in if (!is.null(de) && nrow(de) > 0) {    render_table(de[order(de$n_Linked_Promoters, decreasing = TRUE),         ])}: missing value where TRUE/FALSE needed
#> 17/36                            
#> 18/36 [annotation-plots]         
#> Error: object 'res' not found
#> 19/36                            
#> 20/36 [annotation-karyo]         
#> Error: object 'res' not found
#> 21/36                            
#> 22/36 [refinement-run]           
#> 23/36                            
#> 24/36 [refinement-summary]       
#> 25/36                            
#> 26/36 [refinement-results]       
#> 27/36                            
#> 28/36 [refinement-sankey]        
#> 29/36                            
#> 30/36 [functional-profiling]     
#> 31/36                            
#> 32/36 [profiling-summary]        
#> 33/36                            
#> 34/36 [profiling-results]        
#> 35/36                            
#> 36/36 [sessioninfo]              
#> output file: skeleton.knit.md
#> /opt/hostedtoolcache/pandoc/3.8.3/x64/pandoc +RTS -K512m -RTS skeleton.knit.md --to html4 --from markdown+autolink_bare_uris+tex_math_single_backslash --output /tmp/Rtmpxwk3qx/My_Study_Report.html --lua-filter /home/runner/work/_temp/R-lib/bookdown/rmarkdown/lua/custom-environment.lua --lua-filter /home/runner/work/_temp/R-lib/rmarkdown/rmarkdown/lua/pagebreak.lua --lua-filter /home/runner/work/_temp/R-lib/rmarkdown/rmarkdown/lua/latex-div.lua --lua-filter /home/runner/work/_temp/R-lib/rmarkdown/rmarkdown/lua/table-classes.lua --embed-resources --standalone --wrap preserve --variable bs3=TRUE --section-divs --table-of-contents --toc-depth 3 --variable toc_float=1 --variable toc_selectors=h1,h2,h3 --variable toc_collapsed=1 --variable toc_smooth_scroll=1 --variable toc_print=1 --template /tmp/Rtmpxwk3qx/BiocStyle/template.html --syntax-highlighting none --variable highlightjs=1 --number-sections --variable theme=bootstrap --css /home/runner/work/_temp/R-lib/BiocStyle/resources/html/bioconductor.css --mathjax --variable 'mathjax-url=https://mathjax.rstudio.com/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML' --include-in-header /tmp/Rtmpxwk3qx/rmarkdown-str9b496e71ea3b.html --variable code_folding=hide --variable code_menu=1 
#> 
#> Output created: /tmp/Rtmpxwk3qx/My_Study_Report.html
#> >> Report saved to: /tmp/Rtmpxwk3qx/My_Study_Report.html
# }
```
