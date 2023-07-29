# Single Cell RNAseq Analysis Code

Code for the analysis of single cell RNAseq data in: Nacev, _et. al._, _Cancer-associated Histone H3 N-terminal arginine mutations disrupt PRC2 activity and impair differentiation_

## Installation

### Requirements:

- OS: Linux (tested on CentOS-7 ver 7.9.2009)

- CellRanger Version 6 (https://support.10xgenomics.com)

- R Version >= 4.1.x (https://www.r-project.org)

- Seurat Version 4.2.0 (https://satijalab.org/seurat)


### Setup

To run the code the following directory struction is need as it is hardwired in the code. Create a directory for the analysis, it can be called anything. Here we will refer to by the variable `$AROOT`.

Then in `$AROOT` create the following subfolders:
- `$AROOT\analysis`
- `$AROOT\cellRanger`

In `$AROOT\cellRanger` you need to run `cellranger count` on the FASTQ files. The exact commands we ran are included in the file `CMDS.CellRanger.p12553`. Adjust to the correct path to the FASTQ files on your system however the `id` names must remain the same.


```
git clone -b proj/p12553 git@github.com:soccin/scRNA.git
cp scRNA/pass_00_PARAMS.yaml .

```

## Run Analysis




