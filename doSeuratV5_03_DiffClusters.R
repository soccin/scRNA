suppressPackageStartupMessages(require(stringr))

usage="
usage: doSeuratV5_02.R PARAMS.yaml diffParams.yaml

    PARAMS.yaml     parameter file from pass1
    diffParams.yaml parameters for differential analysis
                    See: scRNA/diffParams.Example.yaml for details
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

if(!is.null(diffParams$metadata)) {
    cat("\n\tReplacing default metdata with",diffParams$metadata,"\n\n")
    md=read_csv(diffParams$metadata) %>% column_to_rownames("CellID") %>% as.data.frame
    so@meta.data=md
}

deTest="wilcox"
if(!is.null(diffParams$method)) {
    cat("\n\tUsing diff method",diffParams$method,"\n\n")
    deTest=diffParams$method
}


Idents(so)=diffParams$groupVar

library(fgsea)
library(msigdbr)

if(args$glbs$genome %in% c("human","xenograft","hg38")) {

    msigdb_species="Homo sapiens"

} else if(args$glbs$genome %in% c("mm10")) {

    msigdb_species="Mus musculus"

} else {

    cat("\n\nUnknown genome for msigdb",args$glbs$genome,"\n\n")
    rlang::abort("FATAL::ERROR")

}

msigdbr_df = msigdbr(species = msigdb_species)
msigdbr_list = split(x = msigdbr_df$gene_symbol, f = msigdbr_df$gs_name)
pathdb=msigdbr_df %>%
    distinct(gs_id,.keep_all=T) %>%
    select(
        -gene_symbol,-entrez_gene,-ensembl_gene,-human_gene_symbol,
        -human_entrez_gene,-human_ensembl_gene
    )

pathways=list()
diffTbl=list()

comps=diffParams$comps %>% map(as_tibble) %>% bind_rows
grpLevels=unique(so@meta.data[[diffParams$groupVar]])

countTbl=so@meta.data %>% count(Group=.[[diffParams$groupVar]],name="Count")

for(ci in transpose(comps)) {

    compName=paste(rev(ci),collapse="_vs_")

    if(!(ci$GroupA %in% grpLevels && ci$GroupB %in% grpLevels)){
        cat("\n\tMissing Group in comparison",compName,"\n")
        cat("\t   ci$GroupA",ci$GroupA,ci$GroupA %in% grpLevels,"\n")
        cat("\t   ci$GroupB",ci$GroupB,ci$GroupB %in% grpLevels,"\n\n")
        next
    }

    fm=NULL

    res=try({fm=FindMarkers(so,ident.1=ci$GroupB,ident.2=ci$GroupA,test.use=deTest)})
    if(class(res)=="try-error") {
        cat("\n\n",compName,res,"\n\n")
        next
    }

    fm=fm %>% rownames_to_column("Gene") %>% arrange(desc(abs(avg_log2FC)))
    colnames(fm)[4]=paste0("pct.",ci$GroupB)
    colnames(fm)[5]=paste0("pct.",ci$GroupA)
    diffTbl[[compName]]=fm %>%
        filter(abs(avg_log2FC)>log2(1.5) & p_val_adj<0.05) %>%
        select(-p_val,-p_val_adj)

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

counts=md %>% count(.[[diffParams$groupVar]])
colnames(counts)[1]=diffParams$groupVar

fix_names_for_excel<-function(ss) {
    gsub("-","_",ss) %>% gsub("_vs_","-",.)
}

compNames=names(diffTbl)
excelNames=
    substr(
        paste0(
            "C",
            sprintf("%02d",seq(compNames)),
            "_",
            fix_names_for_excel(compNames)
            ),
            1,31
        )

cTbl=comps %>% mutate(CompName=cc(GroupB,"vs",GroupA))

compManifest=tibble(TAG=excelNames,CompName=compNames) %>% 
    left_join(cTbl) %>% 
    select(TAG,GroupB,GroupA,CompName)

names(pt_df)=excelNames
names(diffTbl)=excelNames

openxlsx::write.xlsx(
    c(list(comps=compManifest), pt_df),
    cc(args$PROJNAME,"ClusterPathways",diffParams$groupVar,deTest,"V1.xlsx")
    )

openxlsx::write.xlsx(
    c(list(counts=counts,comps=compManifest),diffTbl),
    cc(args$PROJNAME,"DiffGenesSortAbsFoldChange",diffParams$groupVar,deTest,"V1.xlsx")
    )
