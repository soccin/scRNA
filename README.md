# scRNA Version 4

## Branch dev/2023

Major refactor of workflow

Needs R>=4.x.

Major workflow refactor.

- Default is to use simple merge for multiple samples but can set `COMBINE=INTEGRATE` to get the `Seurat` integration workflow. This is now doable from the `CMD.p12` script.

- Same two pass of stage 1 `doSeuratV5_01.R` to do QC and determine filtering options. Filtering levels now computed using `MED +/- 3*MAD` criteria.

- Cell Cycle QC done after filtering in stage 1. And added a script for just filtering.

- Added gene filter to Stage 2a

- QC Filter tables

## Current pipeline

Multiple stages

- `doSeuratV5_01.R`: Initial QC to check filtering paramters and check cell cycle regression plots to see if cell cycle regression is needed.

- Stage 2: Multiple parts now:

    - `doSeuratV5_02a_MergeOnly.R` or `doSeuratV5_02a_IntegrateData.R` for merge or integration steps
    - `doSeuratV5_02b.R` PCA and clustering
