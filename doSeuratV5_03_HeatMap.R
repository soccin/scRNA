suppressPackageStartupMessages(require(stringr))

usage="
usage: doSeuratV5_02.R PARAMS.yaml

    PARAMS.yaml     parameter file from pass1

"

cArgs=commandArgs(trailing=T)

if(len(cArgs)<1) {
    cat(usage)
    quit()
}

if(R.Version()$major<4) {
    cat(usage)
    cat("\n\nThis script needs version(R).major>=4\n\n")
    quit()
}

if(Sys.getenv("SDIR")=="") {
    file.arg=grep("--file=",commandArgs(),value=T)
    if(len(file.arg)>0) {
        SDIR=dirname(gsub(".*=","",file.arg))
    } else {
        SDIR="."
    }
} else {
    SDIR=Sys.getenv("SDIR")
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

obj=readRDS(args$PASS2b.RDAFile)
genes=read_csv(cArgs[2])

sr=subset(obj,cells=Cells(obj)[runif(nrow(obj@meta.data))<0.1])
DefaultAssay(sr)="SCT"

pg=DoHeatmap(sr,assay="SCT",features=genes$Gene) + scale_fill_gradientn(colors = c("blue", "white", "red"))

png(filename="heatmap.png",type="cairo",units="in",width=14,height=10,pointsize=12,res=150)
print(pg)
dev.off()

# clusterRes="RNA_snn_res.0.5"
# s1=obj$s1
# maxClusters=nlevels(s1@meta.data[[clusterRes]])
# pal1=cols25(maxClusters)

# pg=DimPlot(s1, reduction = "umap", group.by="orig.ident", split.by="orig.ident", label.size=6,ncol=2) + scale_color_manual(values=c("darkviolet","red","grey")) + ggtitle(clusterRes)

# pdf(file=get_plot_filename(plotNo(),"UMAP","v2",".pdf"),width=11,height=8.5)
# print(pg)
# dev.off()
