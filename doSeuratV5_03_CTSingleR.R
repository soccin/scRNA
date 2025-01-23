suppressPackageStartupMessages(require(stringr))

usage="
usage: doSeuratV5_03_CTSingleR.R [CLUSTER_RES=res] PARAMS_2b.yaml

    PARAMS_2b.yaml     parameter file from pass2b [post PCA]
    CLUSTER_RES        optional: resolution of clusters to use for cluster level assigments
    ATLAS_TAG          Atlas to use [ImmGenData,MouseRNAseqData]
    ATLAS_LEVEL        Level to use in atlas [main,fine]
"

STAGE=4

cArgs=commandArgs(trailing=T)

#
# Separate out any options arguments
#
optionals=grep("=",cArgs,value=T)

oArgs=list(CLUSTER_RES=NULL,ATLAS_TAG="MouseRNAseqData",ATLAS_LEVEL="main")
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

if(packageVersion("celldex")<"1.12") {
    cat(usage)
    cat("\n\nThis script needs Version(celldex)>=1.12\n")
    cat("Load conda environment here: Work/CONDA/VirtualEnvs/scRNA_V4\n\n")
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
ATLAS_LEVEL=oArgs$ATLAS_LEVEL

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
} else if(args$glbs$genome=="hg38") {
    atlas=celldex::HumanPrimaryCellAtlasData()
    ATLAS_TAG="HumanPrimaryCellAtlas"
} else {
    cat("\n    Genome",args$glbs$genome,"not implemented\n\n")
    rlang::abort("FATAL ERROR::CTSingleR")
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

md$CT=fct_lump_n(factor(md$CT),10)
md$CT[is.na(md$CT)]="Other"

s1@meta.data=md %>% column_to_rownames("CellID")

ctNames=sort(unique(md$CT[!is.na(md$CT)]))

nCols=len(ctNames)

if(nCols<=11) {
    #ctCols=c(RColorBrewer::brewer.pal(10,"Paired"),"grey35")
    ctCols=c(pals::cols25(10),"grey35")
} else {
    ctCols=pals::cols25(nCols)
}
names(ctCols)=ctNames
ctCols=ctCols[sort(unique(md$CT[!is.na(md$CT)]))]

pg=DimPlot(s1,group.by="CT",cols=ctCols)

pdf(file=get_plot_filename(cc("b",plotNo()),"CellTypes","SingleR",ATLAS_TAG,ATLAS_LEVEL,".pdf"),width=12,height=8.5)
print(pg)
dev.off()
write_csv(md,cc("cellTypes_SingleR",ATLAS_TAG,ATLAS_LEVEL,".csv.gz"))

if(ATLAS_TAG=="ImmGenData") {

annote=read_csv(file.path(SDIR,"data","celldex_ImmGenData_Annote.csv.gz")) %>% 
    select(label.main,label.fine,label.ont,onto.name) %>%
    distinct(label.fine,.keep_all=T)

mdfine=md %>% 
    select(CellID,matches("^CT_|^CTC_")) %>% 
    gather(Class,Type,-CellID) %>% 
    left_join(annote,by=c(Type="label.fine")) %>% 
    select(CellID,Class,onto.name) %>% 
    mutate(Class=gsub("_Main","_Fine",Class)) %>%
    spread(Class,onto.name) %>%
    mutate(CT_Fine=factor(CT_Fine) %>% 
                    fct_infreq %>%
                    fct_na_value_to_level("Unknown") %>%
                    fct_lump_n(24) %>%
                    fct_recode(Other="Unknown")
                    ) 

s1@meta.data=left_join(md,mdfine) %>% column_to_rownames("CellID")

pf=DimPlot(s1,group.by="CT_Fine",cols=pals::cols25()) + guides(color=guide_legend(ncol=1,override.aes=list(size=5)))

pdf(file=get_plot_filename(cc("b",plotNo()),"CellTypes","SingleR",ATLAS_TAG,"FineOnto",".pdf"),width=15,height=8.5)
print(pf)
dev.off()

write_csv(left_join(md,mdfine),cc("cellTypes_SingleR",ATLAS_TAG,ATLAS_LEVEL,".csv.gz"))
}
