# scRNA Version 4

## Branch master

Added code to work with xenograft samples (human samples with mouse "contamination") (rebased on 2023-11-20)

## Summary

Needs R>=4.2.

This one has the input fixed to deal with cellranger arc output (Chromium Single Cell Multiome ATAC + Gene Expression) also to deal with external groups running of cellranger (SAIL)

For multiplex runs the auto finding of the data dirs does not work. Will need to create and pass an explicity `config.yaml` file with the fillowing info:

```
genome: mm10
multi:
  gex: "Gene Expression"
inputs:
- dir: /.../cellranger/SAMPID/outs/per_sample_outs/SUBSAMPID/count/sample_filtered_feature_bc_matrix
  sid: SAMID-SUBSAMPID
- dir: <NEXTPATH>
  sid: <NEXTSID>
```

The script 'scRNA/getCRConfigTmpl.sh <CELL_RANGER_DIR>' will create a possible one. Remember to edit the genome line

## Previous Updates

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
