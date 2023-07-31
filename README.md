# Single Cell RNAseq Analysis Code

Code for the analysis of single cell RNAseq data in: Nacev, _et. al._, _Cancer-associated Histone H3 N-terminal arginine mutations disrupt PRC2 activity and impair differentiation_

## Installation

### Requirements:

- OS: Linux CentOS-7 ver 7.9.2009 (may work on other Linux's untested)

- CellRanger Version 6 (https://support.10xgenomics.com)

- R Version >= 4.1.x (https://www.r-project.org)

- Seurat Version 4.2.0 (https://satijalab.org/seurat)

The code follows the standard Seurat workflow: (https://satijalab.org/seurat/articles/pbmc3k_tutorial.html) using simple merging and `SCTransform` for normalization with cell cycle regression.

### Setup

To run the code the following directory struction is need as it is hardwired in the code. Create a directory for the analysis, it can be called anything. Here we will refer to by the variable `$AROOT`.

Then in `$AROOT` create the following subfolders:
- `$AROOT\analysis`
- `$AROOT\cellRanger`

In `$AROOT\cellRanger` you need to run `cellranger count` on the FASTQ files. The exact commands we ran are included in the file `CMDS.CellRanger.p12553`. Adjust to the correct path to the FASTQ files on your system however the `id` names must remain the same.


```
cd $AROOT\analysis

git clone -b proj/p12553 git@github.com:soccin/scRNA.git

cp scRNA/pass_00_PARAMS.yaml .
cp scRNA/pass_02b__MetaData__Res0.5_Cluster_6.csv .
```


## Run Analysis

Once the install is complete you can run the analysis by going to the `$AROOT` folder and running the folling scripts
```
./scRNA/CMDS.p12.CellFilter
Rscript scRNA/doSeuratV5_02d.R pass_02b_PARAMS.yaml 0.2
```

Expected run time:
```
real    39m57.509s

real    12m36.144s
```

Output:
```
seuratQC_p12553_01_CellCycle.pdf
seuratQC_p12553_02_Filter.pdf
pass_01_PARAMS.yaml
pass_01_SObj_9c1ceafabf9f90c54b80e1a146236ea4_d10X.orig_.rda
seuratQC_p12553_11_PostIntegrateCC.pdf
seuratQC_p12553_12_VariableFeatures.pdf
pass_02_PARAMS.yaml
pass_02a_SObj_6197c2f6a624a6175059ae5eb2dfdca7_d10X.integrate_.rda
seuratQC_p12553_21_PCADimMetric.pdf
seuratQC_p12553_22_UMAP_20_.pdf
seuratQC_p12553_23_ClusterChart_20_.pdf
pass_02b_PARAMS.yaml
pass_02b_SObj_958bba8282164eabdbc603b9de801ec1_s1_.rda
tblClusterMarkers_SCT_snn_res.0.2_FDR_0.05_logFC_1.xlsx
seuratQC_p12553_21d_ClusterMarkers_SCT_snn_res.0.2_FDR_0.05_logFC_1.pdf
seuratQC_p12553_22c_ClusterMarkersDot_SCT_snn_res.0.2_FDR_0.05_logFC_1_.pdf
```
