suppressPackageStartupMessages(require(stringr))

usage="
usage: doSeuratV5_03_CTSingleR.R MODULE_FILE=file PARAMS_2b.yaml

    PARAMS_2b.yaml     parameter file from pass2b [post PCA]

"

cArgs=commandArgs(trailing=T)

#
# Separate out any options arguments
#
optionals=grep("=",cArgs,value=T)

oArgs=list(MODULE_FILE=NULL)
if(len(optionals)>0) {
    require(stringr, quietly = T, warn.conflicts=F)
    parseArgs=str_match(optionals,"(.*)=(.*)")
    aa=apply(parseArgs,1,function(x){oArgs[[str_trim(x[2])]]<<-str_trim(x[3])})
}

cArgs=grep("=",cArgs,value=T,invert=T)

if(len(cArgs)!=1) {
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

plotNo<-makeAutoIncrementor(00)

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

s1=readRDS(args$PASS2b.RDAFile)

halt("BREAK")

require(ggforce)

pcc=list()
pltClusterModules<-function(dd,pageI) {
    ggplot(dd,aes(Module,PCT,fill=Module)) +
        geom_bar(stat="identity") +
        facet_wrap_paginate(~Cluster,ncol=3,nrow=2,page=pageI) +
        scale_fill_manual(values=cols25()) +
        theme_light() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        ggtitle(paste("Cluster Resolution",clusterI))
}
for(clusterI in grep("integrated_snn",colnames(s1@meta.data),value=T)) {
    print(clusterI)
    sc=s1@meta.data %>%
        count(.data[[clusterI]],Module,.drop=F) %>%
        tibble %>%
        rename(Cluster=1) %>%
        group_by(Cluster) %>%
        mutate(PCT=n/sum(n)) %>%
        ungroup
    pg=list()
    pg[[1]]=pltClusterModules(sc,1)
    for(jj in 2:n_pages(pg[[1]])) {
        print(jj)
        pg[[jj]]=pltClusterModules(sc,jj)
    }
    pcc=c(pcc,pg)
}

pdf(file=cc("seuratQC",args$PROJNAME,cc("b",plotNo()),"CellTypes",args$algoParams$NDIMS,".pdf"),width=14,height=8.5)
print(pcc)
dev.off()
