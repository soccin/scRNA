suppressPackageStartupMessages(require(stringr))

usage="
usage: doSeuratV5_03_CTSingleR.R [CLUSTER_RES=res] PARAMS_2b.yaml

    PARAMS_2b.yaml     parameter file from pass2b [post PCA]
    CLUSTER_RES        optional: resolution of clusters to use for cluster level assigments

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
    library(SingleR)
    library(celldex)
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

if(args$glbs$genome=="mm10") {
    atlas=celldex::MouseRNAseqData()
} else {
    cat("\n    Genome",args$glbs$genome,"not implemented\n\n")
    halt("FATAL ERROR::CTSingleR")
}

DefaultAssay(s1)="RNA"
s2=DietSeurat(s1)
s2=NormalizeData(s1)
sce=as.SingleCellExperiment(DietSeurat(s2))

# singleR_to_long<-function(pred,method){
#     pred %>%
#         data.frame %>%
#         rownames_to_column("CellID") %>%
#         select(CellID,Labels=pruned.labels) %>%
#         mutate(Method=method) %>%
#         tibble
# }

# pred_cell_main=SingleR(test=sce,ref=atlas,assay.type.test="logcounts",labels=atlas$label.main,prune=T)
# cellTypes=singleR_to_long(pred_cell_main,"Cell.Main")
pred_cell_fine=SingleR::SingleR(test=sce,ref=atlas,assay.type.test="logcounts",labels=atlas$label.fine,prune=T)

md=s1@meta.data %>% data.frame(check.names=F) %>% rownames_to_column("CellID") %>% tibble
md=md %>% left_join(
            pred_cell_fine %>%
            data.frame %>%
            rownames_to_column("CellID") %>%
            select(CellID,CT_Fine=pruned.labels)
        )

for(cres in grep("_res",colnames(md),value=T)) {
    print(cres)
    res=gsub(".*_res.","",cres)
    cpred=SingleR::SingleR(test=sce,ref=atlas,assay.type.test="logcounts",labels=atlas$label.fine,clusters=md[[cres]],prune=T)
    ctbl=cpred %>% data.frame %>% rownames_to_column(cres) %>% select(cres,Labels=pruned.labels)
    colnames(ctbl)[2]=paste0("CTC_Fine_",res)
    md=left_join(md,ctbl)
}

s1@meta.data=md %>% column_to_rownames("CellID")

ctNames=sort(unique(md$CT_Fine))
ctCols=pals::cols25(len(ctNames))
names(ctCols)=ctNames

pg=DimPlot(s1,group.by="CT_Fine",cols=ctCols)

pdf(file=cc("seuratQC",args$PROJNAME,cc("b",plotNo()),"CellTypes","SingleR",".pdf"),width=12,height=8.5)
print(pg)
dev.off()
write_csv(md,"cellTypes_SingleR.csv.gz")


