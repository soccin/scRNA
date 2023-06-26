suppressPackageStartupMessages(require(stringr))

usage="
usage: doSeuratV5_03_CTSingleR.R [CLUSTER_RES=res] PARAMS_2b.yaml

    PARAMS_2b.yaml     parameter file from pass2b [post PCA]
    CLUSTER_RES        optional: resolution of clusters to use for cluster level assigments
    ATLAS_TAG          Atlas to use [ImmGenData,MouseRNAseqData]
"

STAGE=3

cArgs=commandArgs(trailing=T)

#
# Separate out any options arguments
#
optionals=grep("=",cArgs,value=T)

oArgs=list(CLUSTER_RES=NULL,ATLAS_TAG="MouseRNAseqData")
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

ATLAS_TAG=oArgs$ATLAS_TAG


##########################################################################
#
# INCLUDE BREAK
#
##########################################################################

s1=readRDS(args$PASS2b.RDAFile)

if(args$glbs$genome=="mm10") {
    if(ATLAS_TAG=="MouseRNAseqData") {
        atlas=celldex::MouseRNAseqData()
    } else if(ATLAS_TAG=="ImmGenData") {
        atlas=celldex::ImmGenData()
    } else {
        cat("\n\nUnknown ATLAS",ATLAS_TAG,"\n")
        quit()
    }
} else {
    cat("\n    Genome",args$glbs$genome,"not implemented\n\n")
    halt("FATAL ERROR::CTSingleR")
}



DefaultAssay(s1)="RNA"
DEBUG=FALSE
if(DEBUG) {
    set.seed(31415)
    s1.o=s1
    s1=subset(s1,downsample=1000)
}

.seurat_to_sce<-function(sobj) {

    s2=DietSeurat(sobj)
    s2=NormalizeData(s2)
    sce=as.SingleCellExperiment(s2)
    sce

}

library(memoise)

mCache=file.path("/fscratch/socci/_RCache_",digest::digest(getwd(),algo="md5"))
fs::dir_create(mCache)

cdb <- cachem::cache_disk(mCache)
seurat_to_sce=memoise(.seurat_to_sce,cache=cdb)

#halt("DB")

sce=seurat_to_sce(s1)


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

ATLAS_LEVEL="fine"
if(ATLAS_LEVEL=="fine") {
    pred_cell=SingleR::SingleR(test=sce,ref=atlas,assay.type.test="logcounts",labels=atlas$label.fine,prune=T)
} else {
    pred_cell=SingleR::SingleR(test=sce,ref=atlas,assay.type.test="logcounts",labels=atlas$label.main,prune=T)
}

md=s1@meta.data %>% data.frame(check.names=F) %>% rownames_to_column("CellID") %>% tibble
md=md %>% left_join(
            pred_cell %>%
            data.frame %>%
            rownames_to_column("CellID") %>%
            select(CellID,CT_Main=pruned.labels)
        )

for(cres in grep("_res",colnames(md),value=T)) {
    print(cres)
    res=gsub(".*_res.","",cres)
    cpred=SingleR::SingleR(test=sce,ref=atlas,assay.type.test="logcounts",labels=atlas$label.fine,clusters=md[[cres]],prune=T)
    ctbl=cpred %>% data.frame %>% rownames_to_column(cres) %>% select(cres,Labels=pruned.labels)
    colnames(ctbl)[2]=paste0("CTC_Main_",res)
    md=left_join(md,ctbl)
}

md=md %>% mutate(CT=gsub(" \\(.*","",CT_Main))

s1@meta.data=md %>% column_to_rownames("CellID")

ctNames=sort(unique(atlas$label.main))
ctCols=pals::cols25(len(ctNames))
names(ctCols)=ctNames
ctCols=ctCols[sort(unique(md$CT[!is.na(md$CT)]))]

pg=DimPlot(s1,group.by="CT",cols=ctCols)

pdf(file=get_plot_filename(cc("b",plotNo()),"CellTypes","SingleR",ATLAS_TAG,ATLAS_LEVEL,".pdf"),width=12,height=8.5)
print(pg)
dev.off()
write_csv(md,cc("cellTypes_SingleR",ATLAS_TAG,ATLAS_LEVEL,".csv.gz"))

