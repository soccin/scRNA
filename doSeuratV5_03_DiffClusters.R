suppressPackageStartupMessages(require(stringr))

usage="
usage: doSeuratV5_02.R PARAMS.yaml diffParams.yaml

    PARAMS.yaml     parameter file from pass1
    diffParams.yaml parameters for differential analysis

"

cArgs=commandArgs(trailing=T)

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
diffParams=read_yaml(cArgs[2])

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

Idents(so)=diffParams$groupVar

library(fgsea)
library(msigdbr)

if(args$glbs$genome %in% c("human","xenograft")) {

    msigdb_species="Homo sapiens"

} else {

    cat("\n\nUnknown genome for msigdb",args$glbs$genome,"\n\n")
    rlang::abort("FATAL::ERROR")

}

msigdbr_df = msigdbr(species = msigdb_species)
msigdbr_list = split(x = msigdbr_df$gene_symbol, f = msigdbr_df$gs_name)
pathdb=msigdbr_df %>% distinct(gs_id,.keep_all=T) %>% select(-gene_symbol,-entrez_gene,-ensembl_gene,-human_gene_symbol,-human_entrez_gene,-human_ensembl_gene)

pathways=list()
diffTbl=list()

comps=diffParams$comps %>% map(as_tibble) %>% bind_rows

for(ci in transpose(comps)) {

    compName=paste(rev(ci),collapse="_vs_")

    fm=FindMarkers(so,ident.1=ci$GroupB,ident.2=ci$GroupA)
    fm=fm %>% rownames_to_column("Gene") %>% arrange(desc(abs(avg_log2FC)))
    colnames(fm)[4]=paste0("pct.",ci$GroupB)
    colnames(fm)[5]=paste0("pct.",ci$GroupA)
    diffTbl[[compName]]=fm %>% filter(abs(avg_log2FC)>log2(1.5))

    gstats=fm$avg_log2FC
    names(gstats)=fm$Gene
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

openxlsx::write.xlsx(diffTbl,cc(args$PROJNAME,"DiffGenesSortAbsFoldChange","V1.xlsx"))
