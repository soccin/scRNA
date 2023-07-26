#'
#' Phase-II.c
#'
#' Do Find Cluster Markers
#'
STAGE=3

suppressPackageStartupMessages(require(stringr))

usage="
usage: doSeuratV5_02d.R PARAMS.yaml

    PARAMS.yaml     parameter file from pass2b
    ClusterRes      Resolution value of clusters to use
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

library(yaml)
args=read_yaml(cArgs[1])

SDIR="scRNA"

# if(Sys.getenv("SDIR")=="") {
#     #
#     # getSDIR defined in .Rprofile
#     #
#     SDIR=getSDIR()
# } else {
#     SDIR=Sys.getenv("SDIR")
# }

##############################################################################
cat("\n=========================================================\n")
cat(str(args))
cat("\n")
##############################################################################

suppressPackageStartupMessages({
    library(Seurat)
    library(patchwork)
    library(tidyverse)
    library(openxlsx)
    library(pals)
})

source(file.path(SDIR,"seuratTools.R"))
source(file.path(SDIR,"plotTools.R"))

glbs=args$glbs
ap=args$algoParams

plotNo<-makeAutoIncrementor(30)

##########################################################################
#
# INCLUDE BREAK
#
##########################################################################


s1=readRDS(args$PASS2b.RDAFile)

if(len(cArgs)>1) {
    cResValue=cArgs[2]
    clusterRes=grep(cResValue,grep("res\\.",colnames(s1@meta.data),value=T),value=T)
} else {
    cat("\n\nNeed to specify a cluster resolution that will match with grep\n")
    cat("possible resolutions =>\n")
    cat("  ",paste(grep("res\\.",colnames(s1@meta.data),value=T),collapse=", "),"\n\n")
    quit()
}

if(!clusterRes %in% colnames(s1@meta.data)) {
    cat(paste("\n\n", "Invalid Cluster Resolution",clusterRes, "\n"))
    cat("Valid values\n\t",paste(grep("integrated_snn_res.",colnames(s1@meta.data),value=T),collapse="; "),"\n\n")
    rlang::abort("\n\nInvalid Cluster Resolution\n\n")
}

s1=SetIdent(s1,value=clusterRes)

halt("BREAK")

DefaultAssay(s1)="RNA"

s1=FindVariableFeatures(s1)

genesOfInterest=c("Tyk2","Tnfrsf8")
residual.features=union(genesOfInterest,VariableFeatures(s1))

s2=SCTransform(s1,vars.to.regress = c('S.Score', 'G2M.Score'),residual.features=residual.features)

nCells=min(ncol(s2),5000)
sh=s2[,sample(colnames(s1),nCells)]

ph1=DoHeatmap(sh,features=genesOfInterest)
ph2=DoHeatmap(sh,features="Tyk2")

clustTag=gsub("inte.*res.|SCT.*res.","cRes_",clusterRes)

pdf(file=get_plot_filename(plotNo(),"ClusterHeatmapV2",clustTag,".pdf"),width=11,height=8.5)
print(ph1)
print(ph2)
dev.off()

