# scRNA Version 4

## Branch Master - ver/Seurat4_v4

- Added support for multi-output cellranger runs `cellranger multi`

- Merged in version to work with output from SAIL. 

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

