# scRNA Version 4

## Branch proj/p12068_v6 (Manuscript Version)

## Methods

### scRNAseq analysis

The raw sequence data (FASTQ files) were first processed with 10X Genomics CellRanger (v6) software to compute the cell/barcode by gene count matrices using the `cellranger count` command with the refdata-gex-GRCh38-2020-A reference database. The count matrix was then processed with a series of R scripts (R version 4) using the Seurat (version 4) R library38. The cells in each sample were filtered to remove those that failed any of the following QC measures: at least 1,000 detected genes per cell; at least 2,500 UMIs per cell; less than 20% of the cells reads mapping to mitochondrial genes. All ribosomal and mitochondrial genes were also explicitly filtered out of any downstream analysis. Each cell was scored for its cell cycle phase using Seurat’s CellCycleScoring function. All samples were then integrated into one normalized dataset using Seurat’s standard workflow, using SCTransform to normalize the data before merging. This both normalizes and regresses out the cell cycle component. We then ran the PCA projection (RunPCA) and computed the UMAP embedding (RunUMAP) using the first 20 PCA components. The data was then clustered using the FindNeighbors/FindClusters functions with resolutions of 0.1 and 0.5 in the clustering step. We then found cluster specific marker genes which were used to manually annotate cell type identity using comparison to existing relevant single cell dataset annotations. UMAP plots were then created showing the cell level intensity of selected genes. All custom script used in data processing are available here: https://github.com/soccin/scRNA/tree/proj/p12068_v6
