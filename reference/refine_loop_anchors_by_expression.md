# Expression-Aware refinement of loop anchors and target linkages

Integrates quantitative RNA-seq data (e.g., TPM/FPKM) with 3D structural
data to filter and reclassify regulatory elements, deriving a
functionally active regulatory network from physical chromatin contacts.

## Usage

``` r
refine_loop_anchors_by_expression(
  annotation_res,
  expr_matrix_file = NULL,
  sample_columns = NULL,
  threshold = 1,
  unit_type = "TPM",
  species = "hg38",
  out_dir = "./results/filtered",
  project_name = "HiChIP",
  color_palette = "Set2",
  karyo_bin_size = 1e+05,
  reclassify_by_expression = TRUE,
  hub_percentile = 0.95
)
```

## Arguments

- annotation_res:

  List. The raw foundational output object returned by
  [`annotate_peaks_and_loops`](https://zying106.github.io/looplook/reference/annotate_peaks_and_loops.md).

- expr_matrix_file:

  Path to a normalised expression matrix (TPM/FPKM, genes × samples).
  Required for refinement. Default: `NULL`.

- sample_columns:

  Character vector or integer indices. Columns in `expr_matrix_file` to
  average. Default: `NULL`.

- threshold:

  Numeric. Minimum expression (e.g. TPM \> 1) for a gene to be
  considered active. Default: `1`.

- unit_type:

  Character. Expression unit for plot labels (e.g., `"TPM"`). Default:
  `"TPM"`.

- species:

  Character. Genome assembly. One of `"hg38"`, `"hg19"`, `"mm10"`,
  `"mm9"`. Default: `"hg38"`.

- out_dir:

  Character. Output directory for the Excel results file. Default:
  `"./results/filtered"`.

- project_name:

  Character. Prefix for output files (automatically appends
  `"_Filtered"`). Default: `"HiChIP"`.

- color_palette:

  Character. RColorBrewer palette name. Default: `"Set2"`.

- karyo_bin_size:

  Integer. Bin width in bp for karyotype heatmaps. Default: `1e5`.

- reclassify_by_expression:

  Logical. If `TRUE` (default), silent promoters (P) and gene bodies (G)
  are reclassified as eP/eG.

- hub_percentile:

  Numeric (0–1). Node-degree quantile for hub detection. Default:
  `0.95`.

## Value

An invisible named list:

- `loop_annotation` — Filtered 3D network with updated `loop_type`
  (e.g., eP-P).

- `anchor_annotation` — Cluster annotations with expressed targets.

- `promoter_centric_stats` — Gene-level connectivity statistics.

- `distal_element_stats` — Distal-element connectivity statistics.

- `target_annotation` — External features linked to active loop
  components.

- `plot_list` — Named list of ggplot objects (dumbbell, rose,
  karyotype).

Also writes `_Refined_Results.xlsx` to `out_dir`.

## Details

**Algorithmic Framework:**

- **Target Filtration:** Parses merged gene assignments (e.g.,
  `"GeneA;GeneB"`), evaluates individual genes against a defined
  expression threshold, and retains only transcriptionally active
  targets.

- **Biological Reclassification:** Reclassifies physically annotated
  promoters (`P`) and gene bodies (`G`) lacking active transcription as
  enhancer-like regulatory elements (`eP`, `eG`). This adjusts the
  regulatory syntax to reflect functional states (e.g., reannotating a
  silent `P-P` loop to an `eP-P` interaction).

- **Structural Hub Preservation:** Inherits foundational 3D Hub
  classifications (e.g., `Is_High_Connectivity_Gene`) derived from the
  raw physical network. This decouples intrinsic structural network
  topology from tissue-specific transcriptional activation states.

- **External Target Refinement:** Filters auxiliary target mapping
  columns (e.g., `Assigned_Target_Genes_Filled`) based on expression
  criteria, ensuring that mapped 1D genomic features are exclusively
  linked to active genes.

## Examples

``` r
# 1. Get paths to the required example files in the package
rdata_path <- system.file("extdata", "analysis_results.RData", package = "looplook")
expr_path <- system.file("extdata", "example_tpm.txt", package = "looplook")

# Safely load the pre-computed annotation result from RData
temp_env <- new.env()
load(rdata_path, envir = temp_env)
# Extract the first object found in the RData file
raw_annotation <- temp_env[[ls(temp_env)[1]]]

# =========================================================================
# Example : Advanced filtering WITH Transcriptome-Guided Reclassification
# =========================================================================
res_reclassified <- refine_loop_anchors_by_expression(
  annotation_res = raw_annotation,
  expr_matrix_file = expr_path,
  sample_columns = c("con1", "con2"),
  threshold = 1.0,
  unit_type = "TPM",
  species = "hg38",
  out_dir = tempdir(),
  project_name = "Example_Reclassified",
  reclassify_by_expression = TRUE
)
#> >>> [Refinement] Project Name: Example_Reclassified_Filtered
#> >>> [Step 1] Loading Data & Expression Matrix...
#>     [Info] 'a1_id'/'a2_id' columns missing. Reconstructing from coordinates...
#> Warning: Detected 4 column names but the data has 5 columns (i.e. invalid file). Added an extra default column name for the first column which is guessed to be row names or an index. Use setnames() afterwards if this guess is not correct, or fix the file write command that created the file to create a valid file.
#>     >>> Active Genes (> 1 TPM): 605
#> >>> [Step 2] Updating Anchors & Loops...
#> >>> [Step 3] Updating Stats...
#> >>> [Step 4] Refining Target Annotations...
#> >>> [Step 5] Generating Visualizations (Returning plot objects)...
#> >>> [Step 6] Exporting Refined Results...
#>     Excel saved.
#> Refinement Complete.

# View the biologically corrected loop types (e.g., transition from P-P to eP-P)
print(table(res_reclassified$loop_annotation$loop_type))
#> 
#>   E-E   E-G   E-P  E-eG  E-eP   G-G   G-P  G-eG  G-eP   P-P  P-eG  P-eP eG-eG 
#>     4     5    14    11    15     2    19    11    22    13    23    50    23 
#> eG-eP eP-eP 
#>    43    45 
```
