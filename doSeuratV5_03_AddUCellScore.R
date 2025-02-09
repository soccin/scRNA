suppressPackageStartupMessages(require(stringr))

usage="
usage: doSeuratV5_03_AddUCellScore.R [CLUSTER_RES=res] PARAMS_2b.yaml ModuleFile

    CLUSTER_RES     optional: resolution of clusters to use for cluster level assigments
    PARAMS_2b.yaml  parameter file from pass2b
    ModuleFile      File with Module Genes (name of module from filename)

"

cArgs=commandArgs(trailing=T)
#
# Separate out any options arguments
#
optionals=grep("=",cArgs,value=T)

oArgs=list(CLUSTER_RES=NULL)
if(len(optionals)>0) {
    require(stringr, quietly = T, warn.conflicts=F)
    parseArgs=str_match(optionals,"(.*)=(.*)")
    aa=apply(parseArgs,1,function(x){oArgs[[str_trim(x[2])]]<<-str_trim(x[3])})
}

cArgs=grep("=",cArgs,value=T,invert=T)

if(len(cArgs)!=2) {
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

plotNo<-makeAutoIncrementor(70)

moduleFile=cArgs[2]

if(!grepl("\\.(xlsx|csv)$",moduleFile)) {
    cat("\n    Note implemented: Only XLSX modules files currently working\n\n")
    quit()
}

args$MODULES=normalizePath(moduleFile)
modFileExt=tools::file_ext(moduleFile)

suppressPackageStartupMessages({
    library(Seurat)
    library(patchwork)
    library(tidyverse)
    library(readxl)
    library(openxlsx)
    library(pals)
    library(UCell)
})

##########################################################################
#
# INCLUDE BREAK
#
##########################################################################

s1=readRDS(args$PASS2b.RDAFile)

cat("Using SCT assay\n")
DefaultAssay(s1)="SCT"

if(modFileExt=="xlsx") {
    if("GeneList" %in% readxl::excel_sheets(moduleFile)) {
        moduleTbl=read_xlsx(moduleFile,sheet="GeneList")
    } else {
        moduleTbl=read_xlsx(moduleFile)
    }
} else if(modFileExt=="csv") {
    moduleTbl=read_csv(moduleFile)
}

modules=split(moduleTbl$Genes,moduleTbl$Module)
names(modules)=make.names(names(modules))

if(!is.null(oArgs$CLUSTER_RES)) {
    clusterRes=oArgs$CLUSTER_RES
} else {
    clusterRes="0.1"
}
clusterRes=grep(paste0(clusterRes,"$"),colnames(s1@meta.data),value=T)

#
# 
#

s1=AddModuleScore_UCell(s1,features=modules,name="_UC")

