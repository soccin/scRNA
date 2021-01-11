library(Seurat)
library(cowplot)
library(ggplot2)
library(patchwork)

source("seuratTools.R")

## QC
PCT_MITO=10
MIN_FEATURE_RNA=1500
MIN_NCOUNT_RNA=5000

cc.genes=lapply(cc.genes.updated.2019,function(x){convertGeneSymbolsHumanToMouse(x)})

doQCandFilter <- function(so) {

    prjName=as.character(so@meta.data$orig.ident[1])

    pg0=VlnPlot(so, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

    plot1 <- FeatureScatter(so, feature1 = "nCount_RNA", feature2 = "percent.mt")
    plot2 <- FeatureScatter(so, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    #CombinePlots(plots = list(plot1, plot2))

    keep=(
        so@meta.data$nFeature_RNA > MIN_FEATURE_RNA &
        so@meta.data$nCount_RNA > MIN_NCOUNT_RNA &
        so@meta.data$percent.mt < PCT_MITO
        )

    knitr::kable(table(keep))

    so <- subset(so,
            subset = nFeature_RNA > MIN_FEATURE_RNA & nCount_RNA > MIN_NCOUNT_RNA & percent.mt < PCT_MITO
            )

    pg1=VlnPlot(so, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

    pdf(file=cc("seuratQC",prjName,"01.pdf"),height=8.5,width=11)
    cat(cc("seuratQC",prjName,"01.pdf"),"\n")
    print(pg0)
    print(plot1+plot2)
    print(pg1)
    dev.off()

    so

}

preProcessSO<-function(so) {
    so=NormalizeData(so);
    so=FindVariableFeatures(so);
    so=ScaleData(so,features=rownames(so))
}


plotCellCycle<-function(sc) {
    sc=RunPCA(sc,features=c(cc.genes$s.genes,cc.genes$g2m.genes))
    pg=DimPlot(sc,group.by="Phase") + ggtitle(paste(sc@project.name,"Cell Cycle PCA Projection"))
    pg
}

#############################################################################################
#
if(interactive()) {stop("INCLUDE")}
#####

checkFiles=dir(pattern="CHECKPOINT.*.Rdata")
if(len(checkFiles)>1) {
    cat("\n\n Too many CHECKPOINT files need to set explicitly\n")
}
load(checkFiles)


projectTag="p11206"
d10X.integrated@project.name=projectTag
DefaultAssay(d10X.integrated) <- "integrated"

p0 <- DimPlot(d10X.integrated, reduction = "umap", label = TRUE)
p1 <- DimPlot(d10X.integrated, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(d10X.integrated, reduction = "umap", split.by = "orig.ident",ncol=2)
p3 <- plotCellCycle(d10X.integrated)

pdf(file=cc("seuratAnalysis_PostMerge_UMAP",projectTag,".pdf"),width=11,height=8.5)

print(p0)
print(p0+p3)
print(p0+p1)
print(p2)

dev.off()

##
## Find Cluster Markers
##

so=d10X.integrated
DefaultAssay(so)="RNA"
so=NormalizeData(so)
so=FindVariableFeatures(so)
so=ScaleData(so)

if(!interactive()) {quit()}
# stop("\n\nInit and Load\n\n")

so=d10X.integrated
DefaultAssay(so)="RNA"
so=NormalizeData(so)
so=FindVariableFeatures(so)
so=ScaleData(so)

so@meta.data$Types=gsub("\\d$","",so@meta.data$orig.ident)

clusterMarkers=FindAllMarkers(so,only.pos=TRUE,logfc.threshold=0.25,min.pct = 0.25)

top6ClusterMarkers=clusterMarkers %>%
    tibble %>%
    group_by(cluster) %>%
    arrange(cluster,p_val_adj) %>%
    mutate(N=row_number()) %>%
    filter(N<=6)

pdf(file=cc("seuratAnalysis_ClusterMarkers",projectTag,".pdf"),width=11,height=8.5)

for(cii in levels(clusterMarkers$cluster)) {
    print(cii)
    gene.ii=top6ClusterMarkers %>% filter(cluster==cii) %>% pull(gene)
    p1=VlnPlot(so,features=gene.ii,pt.size=.25)
    p2=DotPlot(so,features=gene.ii,split.by="Types",cols=RColorBrewer::brewer.pal(3,"Dark2"))
    print(p1)
    print(p2)
}

dev.off()



# pg.dot=DotPlot(so,features=rev(gg),cols = c("blue", "red"), dot.scale = 10,split.by="orig.ident",group.by="seurat_clusters")+RotatedAxis()

# stop("DMSO")
# pdf(file="diffGenes__pct0.1_IBR10_vs_DMSO__FC_2_FDR_0.05.pdf",height=8.5,width=11)
# print(pg)
# print(pg.dot)
# dev.off()

# cellTypes=paste0("c",d10X.integrated@meta.data$seurat_clusters,"_",d10X.integrated@meta.data$orig.ident)
# cellTypes=factor(cellTypes,levels=sort(unique(cellTypes)))
# d10X.integrated$cellType.Treat=cellTypes
# Idents(d10X.integrated)="cellType.Treat"

# markers=list()

# for(ci in levels(d10X.integrated@meta.data$seurat_clusters)) {
#     grp1=paste0("c",ci,"_IBR10")
# stop("DMSO")
#     grp2=paste0("c",ci,"_DMSO")
#     cat(grp1,"vs",grp2,"...")
#     markers[[paste0("c",ci)]]=FindMarkers(so,ident.1=grp1,ident.2=grp2,pct.min=.1,) %>%
#         rownames_to_column("Gene") %>%
#         as_tibble %>%
#         mutate(avg_log2FC=log2(exp(1))*avg_logFC) %>%
#         mutate(FoldChange=ifelse(avg_log2FC>0,2^avg_log2FC,-2^(-avg_log2FC))) %>%
#         select(Gene,p_val_adj,FoldChange,avg_log2FC,pct.1,pct.2,p_val) %>%
#         mutate(OR.1=pct.1/(1-pct.1),OR.2=pct.2/(1-pct.2),lOR=log10(OR.1/OR.2)) %>%
#         arrange(p_val,desc(abs(FoldChange)))

#     cat("\n")
# }



save.image(cc("CHECKPOINT",DATE(),digest::digest(d10X),".Rdata"),compress=T)

