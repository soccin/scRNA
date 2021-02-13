suppressPackageStartupMessages({
    library(tidyverse)
    library(Seurat)
    library(patchwork)
    library(pals)
})

source("seuratTools.R");
source("tools.R")

pipeTag=paste0(basename(getwd()),"_v2")

d10X.integrated=readRDS("obj__d10X.integrated.rda")

md=tibble(d10X.integrated@meta.data) %>%
    mutate(sample=gsub("\\d+$","",orig.ident) %>% gsub("\\.$","",.)) %>%
    mutate(sample=factor(sample,levels=c("Sham","RT","Tmem.RT"))) %>%
    mutate(cluster=forcats::fct_relevel(seurat_clusters,rev))

tbl=md %>% count(sample,cluster) %>% group_by(sample) %>% mutate(total=sum(n)) %>% mutate(P=n/total)

numClust=nlevels(md$cluster)

source("getCellTypesGenesLists.R")
features=getCellTypeGeneLists(rownames(d10X.integrated@assays$RNA@counts))

defaultAssay="Integrated"
defaultAssay="RNA"

if(defaultAssay=="RNA") {
    so=d10X.integrated
    DefaultAssay(so)="RNA"
    so=NormalizeData(so)
    so=FindVariableFeatures(so)
    so=ScaleData(so)
    so@meta.data$sample=gsub("\\d$","",so@meta.data$orig.ident)

    xx=AddModuleScore(object=so,features=features,name="CellType.Score.")
} else {
    xx=AddModuleScore(object=d10X.integrated,features=features,name="CellType.Score.")
}

cellTypes=xx@meta.data %>%
    rownames_to_column("CellID") %>%
    tibble %>% select(CellID,matches("\\.Score\\.")) %>%
    gather(Score,Value,-CellID) %>%
    group_by(CellID) %>%
    arrange(desc(Value)) %>%
    slice(1)

newMetaData=xx@meta.data %>%
    rownames_to_column("CellID") %>%
    tibble %>%
    left_join(cellTypes) %>%
    select(-matches("CellType.Score.")) %>%
    mutate(CellType=names(features)[gsub(".*Score.","",Score) %>% as.numeric]) %>%
    mutate(CellType=ifelse(Value>0,CellType,"UNK")) %>%
    mutate(sample=gsub("[.0-9]+$","",orig.ident)) %>%
    mutate(sample=factor(sample,levels=c("Sham","RT","Tmem.RT"))) %>%
    data.frame %>%
    column_to_rownames("CellID")

if(defaultAssay=="Integrated") {
    so=d10X.integrated
}

so@meta.data=newMetaData

nCellTypes=len(unique(cellTypes$Score))
clusterTypeTbl=tibble(newMetaData) %>%
    count(seurat_clusters,CellType) %>%
    group_by(seurat_clusters) %>%
    mutate(PCT=n/sum(n))

cellTypes=clusterTypeTbl %>% ungroup %>% distinct(CellType) %>% pull %>% sort

clusterTypeTbl=clusterTypeTbl %>%
    mutate(CellType=factor(CellType,levels=rev(cellTypes)))

superMajority=clusterTypeTbl %>% arrange(desc(PCT)) %>% filter(PCT>.6)

pass1=superMajority %>% pull(seurat_clusters)

