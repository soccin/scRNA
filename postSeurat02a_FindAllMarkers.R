suppressPackageStartupMessages({
    library(tidyverse)
    library(Seurat)
    library(patchwork)
    library(pals)
    require(dplyr)
    require(openxlsx)
    library(monocle)
})

source("seuratTools.R");
source("tools.R")

pipeTag=paste0(basename(getwd()),"_v2")

d10X.integrated=readRDS("obj__d10X.integrated.rda")

so=d10X.integrated
DefaultAssay(so)="RNA"
so=NormalizeData(so)
so=FindVariableFeatures(so)
so=ScaleData(so)
so@meta.data$sample=gsub("\\d+$","",so@meta.data$orig.ident) %>% gsub("\\.$","",.)

options("mc.cores"=12)
clusterMarkers=FindAllMarkers(so,only.pos=TRUE,logfc.threshold=0.25,min.pct = 0.25)

pct=c(clusterMarkers$pct.1,clusterMarkers$pct.2)
ps=min(min(pct[pct>0]),min(1-pct[pct<1]))/2

FDR.cut=0.05
logFC.cut=1
filterCLTable<-function(clm) {
    tibble(clm) %>%
        mutate(lor.1=log((pct.1+ps)/(1-pct.1+ps))) %>%
        mutate(lor.2=log((pct.2+ps)/(1-pct.2+ps))) %>%
        mutate(lOR=lor.1-lor.2) %>%
        arrange(desc(lOR))
}

cl=filterCLTable(clusterMarkers) %>% select(cluster,gene,p_val_adj,avg_logFC,pct.1,pct.2,lOR)

top6ClusterMarkers=clusterMarkers %>%
    tibble %>%
    group_by(cluster) %>%
    arrange(cluster,p_val_adj) %>%
    mutate(N=row_number()) %>%
    filter(N<=6)

cl=cl %>%
    full_join(top6ClusterMarkers %>% select(cluster,gene,N),by=c("cluster","gene")) %>%
    rename(Top6.Orig=N) %>%
    filter(!is.na(Top6.Orig) | (p_val_adj<FDR.cut & avg_logFC>logFC.cut))

clusters=cl %>% distinct(cluster) %>% arrange(cluster) %>% pull %>% levels

clusterMarkerTbl=list()
for(ci in clusters) {
    clusterMarkerTbl[[ci]]=cl %>% filter(cluster==ci)
}

geneCounts=cl %>% count(cluster)




ll=c(list(GeneCounts=geneCounts,AllCluster=cl),clusterMarkerTbl)

xfile=cc("tblClusterMarkers",pipeTag,"FDR",FDR.cut,"logFC",logFC.cut)

write.xlsx(ll,paste0(xfile,".xlsx"))


data <- as(as.matrix(so@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = so@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

cds <- newCellDataSet(data,
                         phenoData = pd,
                         featureData = fd,
                         lowerDetectionLimit = 0.5,
                         expressionFamily = negbinomial.size())

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

cds <- detectGenes(cds, min_expr = 0.1)

expressed_genes <- row.names(subset(fData(cds), num_cells_expressed >= 10))
markerGenes=cl %>% distinct(gene) %>% arrange(gene) %>% pull


# > system.time({dd=differentialGeneTest(cds[markerGenes,],fullModelFormulaStr="~seurat_clusters",cores=4)})
#    user  system elapsed
#   8.257   3.817  87.003

#dm=differentialGeneTest(cds[markerGenes,],fullModelFormulaStr="~seurat_clusters",cores=4)

#dex=differentialGeneTest(cds[expressed_genes,],fullModelFormulaStr="~seurat_clusters",cores=8)

##
# timing

# > dim(cds)
# Features  Samples
#    32285    17894

# > len(ic)
# [1] 1000

# > system.time({dd=differentialGeneTest(cds[markerGenes,ic],fullModelFormulaStr="~seurat_clusters",cores=4)})
#    user  system elapsed
#   7.473   3.911  32.357

# > len(ic)
# [1] 2000

# > system.time({dd=differentialGeneTest(cds[markerGenes,ic],fullModelFormulaStr="~seurat_clusters",cores=4)})
#    user  system elapsed
#   6.852   3.219  34.873

# > len(ic)
# [1] 4000

# > system.time({dd=differentialGeneTest(cds[markerGenes,ic],fullModelFormulaStr="~seurat_clusters",cores=4)})
#    user  system elapsed
#   7.569   3.241  40.288
