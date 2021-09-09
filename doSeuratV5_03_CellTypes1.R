suppressPackageStartupMessages(require(stringr))

usage="
usage: doSeuratV5_03_CellTypes1.R PARAMS.yaml

    PARAMS.yaml     parameter file from pass2

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

plotNo<-makeAutoIncrementor(30)

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

if(is.null(args$optionals$MODULE_FILE)) {
    cat("\n\nNeed module file and module scores\n\n")
    quit()
}

moduleTbl=read_tsv(args$optionals$MODULE_FILE)
colnames(moduleTbl)=c("Module","Gene")
modules=split(moduleTbl$Gene,moduleTbl$Module)

renameModules=names(modules)
names(renameModules)=paste0("Modules",seq(names(modules)))

ss=s1@meta.data %>%
    rownames_to_column("CellID") %>%
    select(CellID,matches("^Modules")) %>%
    gather(Module,Score,-CellID) %>%
    mutate(Module=recode(Module,!!!renameModules)) %>%
    arrange(CellID,desc(Score),Module) %>%
    group_by(CellID) %>%
    mutate(Rank=row_number()) %>%
    filter(Rank==1) %>%
    select(-Rank) %>%
    ungroup() %>%
    column_to_rownames("CellID")

ss=ss[rownames(s1@meta.data),]

s1=AddMetaData(object=s1,metadata=factor(ss$Module),col.name="Module")
s1=AddMetaData(object=s2,metadata=ss$Score,col.name="Module.Score")





