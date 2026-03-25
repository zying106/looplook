# =========================================================================
# Script and provenance for generating example data in `inst/extdata`
# Package: looplook
# =========================================================================
#
# OVERVIEW:
# The data provided in `inst/extdata` are strictly minimalistic, downsampled
# subsets of genomic and transcriptomic datasets. They are solely intended
# for demonstrating package functionality, running vignettes, and executing
# unit tests within an acceptable time frame.
#
# LICENSE:
# All example datasets are distributed under the same license as the
# `looplook` package and are free for open-source usage.
#
# =========================================================================
# 1. ENCODE Blacklist Regions
# Files: hg19-blacklist.v2.bed, hg38-blacklist.v2.bed,
#        mm9-blacklist.v2.bed, mm10-blacklist.v2.bed
# =========================================================================
# Source: Standard ENCODE universally blocked regions (v2).
# Downloaded directly from the ENCODE portal and Boyle Lab GitHub repository.
# URL: https://github.com/Boyle-Lab/Blacklist/
# These are kept in their original BED format as standard reference files.
#
# =========================================================================
# 2. 3D Genomic and Epigenomic Data (HiChIP / ChIP-seq)
# Files: example_peaks.bed, example_k27ac_peaks.bed,
#        example_loops_1.bedpe, example_loops_2.bedpe,
#        example_loops_H3K27ac.bedpe, example_loops_pol2.bedpe
# =========================================================================
# Source: Derived from publicly available HiChIP and ChIP-seq datasets (GEO: GSE213300,GSE111253).
#
# Raw data processing pipeline (Bash Pseudo-code):
# ```bash
# # Alignment and Peak Calling
# bowtie2 -p 25 -x ${index} -U ${fq1} | samtools sort -O bam -@ 10 -o - > ./mapping/${sample}_bowtie2.bam
# macs2 callpeak -t ../mapping/${name}_bowtie2.bam.unique.bam_rmdup.bam -c ../mapping/${input}_bowtie2.bam.unique.bam_rmdup.bam -f AUTO -g hs --nomodel --extsize ${lenth} -n ${name} -B --SPMR -q 0.01 --outdir ./ 2>./${name}_macs2Peak_summary.txt
# ```
#
# Subsetting process for example data (Bash Pseudo-code):
# ```bash
# # Randomly sample exact number of lines using `shuf` to minimize package size
# shuf -n 3wc-l00 raw_peaks_full.bed > example_peaks.bed
# shuf -n 600  raw_k27ac_peaks_full.bed > example_k27ac_peaks.bed
#
# shuf -n 300 raw_loops_1_full.bedpe > example_loops_1.bedpe
# shuf -n 300 raw_loops_2_full.bedpe > example_loops_2.bedpe
#
# shuf -n 300  raw_H3K27ac_loops_full.bedpe > example_loops_H3K27ac.bedpe
# shuf -n 300  raw_Pol2_loops_full.bedpe > example_loops_pol2.bedpe
# ```
#
# =========================================================================
# 3. Pre-computed Result Objects
# File: analysis_results.RData
# =========================================================================
# Source: Generated internally by running `looplook` core functions on the
# example data above.
# Reason: Provided to speed up unit tests and vignette compilation, avoiding
# redundant heavy computations during package checks.
#
# Generation process (R Pseudo-code):
# ```R
# library(looplook)
#
# # 1. Get paths to example files included in the package
# bedpe_path <- system.file("extdata", "example_loops_1.bedpe", package = "looplook")
# bed_path <- system.file("extdata", "example_peaks.bed", package = "looplook")
#
# # 2. Check if files and required annotation databases exist
# if (bedpe_path != "" && bed_path != "" &&
#     requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE) &&
#     requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
#
#   # =========================================================================
#   # Example A: Integrative Analysis (Loops + Target BED like ATAC-seq)
#   # =========================================================================
#   res_integrated <- annotate_peaks_and_loops(
#     bedpe_file = bedpe_path,
#     target_bed = bed_path,
#     species = "hg38",
#     tss_region = c(-2000, 2000),
#     out_dir = tempdir(),
#     color_palette = "Set2",
#     karyo_bin_size = 1e5,
#     neighbor_hop = 0,
#     hub_percentile = 0.95,
#     project_name = "Example_HiChIP_Integrative"
#   )
#
#   # Save the output for downstream testing
#   save(res_integrated, file = "inst/extdata/analysis_results.RData")
# }
# ```
#
# =========================================================================
# 4. Transcriptomic Data (RNA-seq)
# Files: example_tpm.txt, example_deg.txt, example_coldata.txt
# =========================================================================
# Source: Derived from publicly available RNA-seq datasets (GEO: GSE111253).
#
# Raw data processing pipeline (Bash Pseudo-code):
# ```bash
# # Mapping and quantification using STAR and RSEM
# STAR --runThreadN 12 --genomeDir path/to/STAR_index --readFilesCommand zcat --sjdbOverhang 99 --sjdbGTFfile path/to/gtf --readFilesIn ${fq1} ${fq2} --outFileNamePrefix ./mapping/${sample} --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 8 --quantMode TranscriptomeSAM GeneCounts
# rsem-calculate-expression --forward-prob 0.5 --paired-end --no-bam-output --alignments -p 20 -q ./mapping/${sample}Aligned.toTranscriptome.out.bam path/to/RSEM_ref ./rsem/${sample}
# ```
#
# Subsetting process for example data (R Pseudo-code):
# ```R
# # To match the downsampled 3D data, we filter the full RNA-seq matrix
# # to include only the genes assigned in `res_integrated` (generated in Section 3).
#
# # 1. Load full TPM matrix and DESeq2 output
# full_tpm <- read.delim("raw_data/Full_TPM.txt")
# full_deg <- read.csv("raw_data/Full_DESeq2_res.csv")
#
# # 2. Subset TPM matrix
# target_genes <- res_integrated$target_annotation$SYMBOL
# tpm <- full_tpm[row.names(full_tpm) %in% target_genes, 1:4, drop = FALSE]
# colnames(tpm) <- c("con1", "con2", "trt1", "trt2")
#
# write.table(tpm, "inst/extdata/example_tpm.txt",
#             quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t")
#
# # 3. Subset and clean DESeq2 results
# deg <- full_deg[full_deg$Gene_name %in% target_genes, c(5, 17:19)]
# deg[is.na(deg)] <- 0
# deg <- deg[!duplicated(deg$Gene_name), ]
#
# row.names(deg) <- deg$Gene_name
# deg <- deg[, -1]
#
# write.table(deg, "inst/extdata/example_deg.txt",
#             quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t")
# ```
# 4. Subset Sample Metadata
# full_coldata <- read.delim("raw_data/Full_coldata.txt")
# coldata <- full_coldata[full_coldata$SampleID %in% colnames(tpm), ]
# write.table(coldata, "inst/extdata/example_coldata.txt",
#             quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
