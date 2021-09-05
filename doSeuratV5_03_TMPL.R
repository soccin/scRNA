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

if(R.Version()$major<4) {
    cat(usage)
    cat("\n\nThis script needs version(R).major>=4\n\n")
    quit()
}

source(file.path(SDIR,"seuratTools.R"))
source(file.path(SDIR,"plotTools.R"))

library(yaml)
args=read_yaml(cArgs[1])

glbs=args$glbs
ap=args$algoParams

plotNo<-makeAutoIncrementor(20)

suppressPackageStartupMessages({
    library(Seurat)
    library(patchwork)
    library(tidyverse)
    library(openxlsx)
    library(pals)
})



##########################################################################
#
# INCLUDE BREAK
#
##########################################################################

obj=readRDS(args$PASS2.RDAFile)



# clusterRes="RNA_snn_res.0.5"
# s1=obj$s1
# maxClusters=nlevels(s1@meta.data[[clusterRes]])
# pal1=cols25(maxClusters)

# pg=DimPlot(s1, reduction = "umap", group.by="orig.ident", split.by="orig.ident", label.size=6,ncol=2) + scale_color_manual(values=c("darkviolet","red","grey")) + ggtitle(clusterRes)

# pdf(file=cc("seuratQC",args$PROJNAME,plotNo(),"UMAP","v2",".pdf"),width=11,height=8.5)
# print(pg)
# dev.off()
