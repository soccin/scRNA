# scRNA Version 4

## Branch MASTER

10X data

## Current pipeline

Multiple stages

- `doSeuratV5_01.R`: Initial QC to check filtering paramters and check cell cycle regression plots to see if cell cycle regression is needed.

- `doSeuratV5_02.R`: Stage 2, redo the QC-filter (but can reset parameters) but now move the Cell Cycle regression into the main script and follow the Seurat suggestion of only regressing out the variable features.



