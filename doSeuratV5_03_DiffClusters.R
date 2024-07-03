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

if(Sys.getenv("SDIR")=="") {
    file.arg=grep("--file=",commandArgs(),value=T)
    if(len(file.arg)>0) {
        SDIR=dirname(gsub(".*=","",file.arg))
    } else {
        if(len(fs::dir_ls(regex="scRNA"))>0) {
            SDIR=fs::dir_ls(regex="scRNA")
        } else {
            SDIR="."
        }
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

obj=readRDS(args$PASS2b.RDAFile)

so=obj

so@meta.data$Clusters=so@meta.data$SCT_snn_res.0.1
Idents(so)="Clusters"

library(fgsea)
library(msigdbr)

msigdbr_df = msigdbr(species = "Homo sapiens")
msigdbr_list = split(x = msigdbr_df$gene_symbol, f = msigdbr_df$gs_name)
pathdb=msigdbr_df %>% distinct(gs_id,.keep_all=T) %>% select(-gene_symbol,-entrez_gene,-ensembl_gene,-human_gene_symbol,-human_entrez_gene,-human_ensembl_gene)

pathways=list()
comps=tribble(~ClustA,~ClustB,1,0,1,2,1,3)

for(ci in transpose(comps)) {

    compName=paste(rev(ci),collapse="_vs_")

    fm=FindMarkers(so,ident.1=ci$ClustA,ident.2=ci$ClustB)
    gstats=fm$avg_log2FC
    names(gstats)=rownames(fm)
    fg=fgsea(msigdbr_list,gstats,minSize=15,maxSize=500)

    pt=tibble(fg) %>%
        filter(padj<0.05) %>%
        arrange(pval) %>%
        rowwise %>%
        mutate(leadingEdge=paste((leadingEdge),collapse=";")) %>%
        left_join(pathdb,by=c(pathway="gs_name"))

    pathways[[compName]]=pt

}

pt_df=map(pathways,data.frame)

openxlsx::write.xlsx(pt_df,cc(args$PROJNAME,"ClusterPathways","V1.xlsx"))
