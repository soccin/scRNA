suppressPackageStartupMessages({
    library(tidyverse)
    library(Seurat)
    library(patchwork)
})

source("seuratTools.R");
source("tools.R")

d10X.integrated=readRDS("obj__d10X.integrated.rda")

md=tibble(d10X.integrated@meta.data) %>%
    mutate(sample=gsub("\\d+$","",orig.ident) %>% gsub("\\.$","",.)) %>%
    mutate(sample=factor(sample,levels=c("Sham","RT","Tmem.RT"))) %>%
    mutate(cluster=forcats::fct_relevel(seurat_clusters,rev))

