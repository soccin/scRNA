# Methods

## Alignment and Gene Counts

The FASTQ files were first processed using CellRanger V6 with the mouse genome and transcript reference data (refdata-gex-mm10-2020-A) supplied by 10X genomics. The output from CellRanger was then analyzed using the Seurat R package version 4. The filtered barcode/feature matrix was read using Seurat read10X command and the percentage of mitochondiral reads per cell was computed for subsequence filtering. Cells were filtered using the following criteria: percent mitochondiral read greater than 20% or the number of genes was less than 1500 or the number of distinct molecules was less than 5000. If any of these conditions were true the cell was filtered out of the downstream analysis. Next we scored the cell cycle phase of each cell usign Seurat's scoreCellCycle functions. The list of cell cycle genes for mouse was determined by taking the supplied list of human cell cycle genes and mapping them to mouse genes using the bioMart R package to do the homology mapping.

## Seurat Analysis

The filtered data was normalized and scaled using the SCTransform method from Seurat. For the transform we regressed against the cell cycle scores previously computed. After normlazation we computed the PCA coordinates and retained the first 20 coordiantes in the clustering and projection analysis. For cluster we used the Seurat FindNeighbor and FindClusters functions with several resolution values and after manual inspection fixed on a resolution value of 0.2 for subsequent work. We also computed the UMAP project using RunUMAP and 20 pca coordinates. Cluster specific marker genes were computed using FindAllMarkers with a cutoff of 0.25 in the log fold change and a minimum percentage of 25%. For two specifed gene sets: Adipocyte and Skeletal Muscle, we computed a score for each using AddModuleScore.

