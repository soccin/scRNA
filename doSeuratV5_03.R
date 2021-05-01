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

source("seuratTools.R")
source("tools.R")
source("doQCandFilter.R")

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
s1=obj$s1

genes=scan("geneList.txt","")
pg=FeaturePlot(s1,features=genes,combine=F)
pg1=VlnPlot(s1,features=genes,combine=F,pt.size=.025)

pgC=list()
for(ii in seq(len(pg))) {
    pgC[[len(pgC)+1]]=pg[[ii]]
    pgC[[len(pgC)+1]]=pg1[[ii]]
}
pp=paginatePlots(pgC,2,2)
pdf(file=cc("seuratQC",args$PROJNAME,plotNo(),"GenePlot","v2",".pdf"),width=11,height=8.5)
print(pp)
dev.off()

genes1=intersect(rownames(obj$s1@assays$RNA@counts),genes)

xx=t(as.matrix(obj$s1@assays$RNA@counts[genes1,]))
rownames(xx)=map(rownames(xx),digest::digest,algo="md5") %>% unlist
xx=xx[order(-rowSums(xx)),]



# clusterRes="RNA_snn_res.0.5"
# s1=obj$s1
# maxClusters=nlevels(s1@meta.data[[clusterRes]])
# pal1=cols25(maxClusters)

# pg=DimPlot(s1, reduction = "umap", group.by="orig.ident", split.by="orig.ident", label.size=6,ncol=2) + scale_color_manual(values=c("darkviolet","red","grey")) + ggtitle(clusterRes)

# pdf(file=cc("seuratQC",args$PROJNAME,plotNo(),"UMAP","v2",".pdf"),width=11,height=8.5)
# print(pg)
# dev.off()
