suppressPackageStartupMessages(require(stringr))

usage="
usage: doSeuratV5_02.R PARAMS.yaml

    PARAMS.yaml     parameter file from pass1

"

cArgs=commandArgs(trailing=T)

if(len(cArgs)!=1) {
    cat(usage)
    quit()
}

library(yaml)
args=read_yaml(cArgs[1])

##############################################################################
cat("\n=========================================================\n")
cat(str(args))
cat("\n")
##############################################################################

suppressPackageStartupMessages({
    library(Seurat)
    library(patchwork)
    library(tidyverse)
})

source("seuratTools.R")

glbs=args$glbs
ap=args$algoParams

plotNo<-makeAutoIncrementor(10)

d10X=readRDS(args$PASS1.RDAFile)

if(args$MERGE & len(d10X)>1) {
    cat("\nMerging sample files...")
    s.merge=merge(d10X[[1]],d10X[-1],project=args$PROJNAME)
    d10X=list()
    d10X[[args$PROJNAME]]=s.merge
    cat("done\n\n")
}

glb.digest=digest::digest(d10X)
cat("digest=",digest::digest(d10X),"\n")

cat("\nDoQCandFilter\n")
for(ii in seq(d10X)) {
    print(ii)
    ret=doQCandFilter(d10X[[ii]], ap$MIN_NCOUNT_RNA, ap$MIN_FEATURE_RNA, ap$PCT_MITO)
    d10X[[ii]]=ret$so
}

if(args$DEBUG) {

    c("\n\nDEBUG SET; Subset data\n\n")
    set.seed(ap$SEED)
    cat("\nDEBUG::subset\n")
    for(ii in seq(d10X)) {
        print(ii)
        xx=d10X[[ii]]
        d10X[[ii]]=subset(xx,cells=Cells(xx)[runif(nrow(xx@meta.data))<args$DOWNSAMPLE])
    }
    cat("digest=",digest::digest(d10X),"\n")

}

if(len(d10X)>1) {
    cat("\n\n The rest of this workflow only works on merged datasets\n\n")
    quit()
}

stop("DDDDDDDDD")

ret=regressCellCycle(d10X[[1]])

s1=ret$so

pdf(file=cc("seuratQC",plotNo(),"PostCCRegress.pdf"),width=11,height=8.5)

plotCellCycle(s1,"Post CC Regression")

dev.off()

# PCA

# https://satijalab.org/seurat/archive/v3.0/pbmc3k_tutorial.html
# Perform linear dimensional reduction

