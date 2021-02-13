suppressPackageStartupMessages(require(stringr))

cArgs=commandArgs(trailing=T)
args=list(DEBUG=FALSE)

ii=grep("=",cArgs)
if(len(ii)>0) {
    parseArgs=str_match(cArgs[ii],"(.*)=(.*)")
    aa=apply(parseArgs,1,function(x){args[[str_trim(x[2])]]<<-str_trim(x[3])})
}
args$DEBUG=as.logical(args$DEBUG)

cat("\n=========================================================\n")
cat(str(args))
cat("\n")

args$DEBUG=as.logical(args$DEBUG)

##############################################################################
##############################################################################

suppressPackageStartupMessages({
    library(Seurat)
    library(cowplot)
    library(ggplot2)
    library(patchwork)
})

source("seuratTools.R")

dataFolders=scan("dataX10InputDirs","")
names(dataFolders)=gsub("_",".",gsub(".outs.*","",gsub(".*/s_","",dataFolders)))

#dmso <- read10XIntoSeuratObj("data/HG19+MM10/s_DMSO/outs/filtered_feature_bc_matrix","DMSO")
#ibr10 <- read10XIntoSeuratObj("data/HG19+MM10/s_Ibr10/outs/filtered_feature_bc_matrix","IBR10")

d10X=list()

for(ii in seq(len(dataFolders))) {
    sampleName=names(dataFolders)[ii]
    cat("Reading Sample =",sampleName,"...")
    d10X[[sampleName]]=read10XMouseIntoSeuratObj(dataFolders[ii],sampleName)
    cat("\n")
}

cat("digest=",digest::digest(d10X),"\n")

## QC
PCT_MITO=10
MIN_FEATURE_RNA=1500
MIN_NCOUNT_RNA=5000

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

# dmso=doQCandFilter(dmso)
# ibr10=doQCandFilter(ibr10)

cat("\nDoQCandFilter\n")
for(ii in seq(d10X)) {
    print(ii)
    d10X[[ii]]=doQCandFilter(d10X[[ii]])
}

if(args$DEBUG) {

    c("\n\nDEBUG SET; Subset data\n\n")

    set.seed(101)

# xx=dmso
# dmso=subset(xx,cells=Cells(xx)[runif(nrow(xx@meta.data))<.05])
# xx=ibr10
# ibr10=subset(xx,cells=Cells(xx)[runif(nrow(xx@meta.data))<.05])

    cat("\nDEBUG::subset\n")

    for(ii in seq(d10X)) {
        print(ii)
        xx=d10X[[ii]]
        d10X[[ii]]=subset(xx,cells=Cells(xx)[runif(nrow(xx@meta.data))<.1])
    }

    cat("digest=",digest::digest(d10X),"\n")

}

# dmso=scoreCellCycle(dmso)
# ibr10=scoreCellCycle(ibr10)
cat("\nScoreCellCycle\n")
for(ii in seq(d10X)) {
    print(ii)
    d10X[[ii]]=scoreCellCycle(d10X[[ii]])
}

preProcessSO<-function(so) {
    so=NormalizeData(so);
    so=FindVariableFeatures(so);
    so=ScaleData(so,features=rownames(so))
}

cc.genes=lapply(cc.genes.updated.2019,function(x){convertGeneSymbolsHumanToMouse(x)})

plotCellCycle<-function(sc) {
    sc=RunPCA(sc,features=c(cc.genes$s.genes,cc.genes$g2m.genes))
    pg=DimPlot(sc,group.by="Phase") + ggtitle(paste(sc@project.name,"Cell Cycle PCA Projection"))
    pg
}

# pcc.dmso=plotCellCycle(preProcessSO(dmso))
# pcc.ibr10=plotCellCycle(preProcessSO(ibr10))

pcc=list()
cat("\nPlotCellCycle\n")
for(ii in seq(d10X)) {
    print(ii)
    pcc[[ii]]=plotCellCycle(preProcessSO(d10X[[ii]]))
}

save.image(cc("CHECKPOINT",DATE(),digest::digest(d10X),".Rdata"),compress=T)

pdf(file="seuratQC_CellCycle.pdf",width=8.5,height=11)
    wrap_plots(pcc[1:4],ncol=2)
    wrap_plots(pcc[5:6],ncol=2)
dev.off()

#
# SCTransform Normalizes
#
d10X.list=d10X
for (i in 1:length(d10X.list)) {
    d10X.list[[i]] <- SCTransform(d10X.list[[i]], vars.to.regress = c('S.Score', 'G2M.Score'), verbose = T)
}

d10X.features <- SelectIntegrationFeatures(object.list = d10X.list, nfeatures = 3000)
d10X.list <- PrepSCTIntegration(object.list = d10X.list, anchor.features = d10X.features, verbose=T)


d10X.anchors <- FindIntegrationAnchors(object.list = d10X.list, normalization.method = "SCT",
    anchor.features = d10X.features, verbose = T)
d10X.integrated <- IntegrateData(anchorset = d10X.anchors, normalization.method = "SCT",
    verbose = T)

DefaultAssay(d10X.integrated) <- "integrated"

# This shoudld not be needed here, data scaled prior to IntegrateData
# see (https://satijalab.org/seurat/v3.0/integration.html)
#
#d10X.integrated <- ScaleData(d10X.integrated)

d10X.integrated <- RunPCA(d10X.integrated,verbose=T)
d10X.integrated <- RunUMAP(d10X.integrated, reduction = "pca", dims = 1:20)
d10X.integrated <- FindNeighbors(d10X.integrated, reduction = "pca", dims = 1:20)
d10X.integrated <- FindClusters(d10X.integrated, resolution = 0.2)

p1 <- DimPlot(d10X.integrated, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(d10X.integrated, reduction = "umap", label = TRUE)

pdf(file="seuratQC_PostMerge.pdf",width=11,height=8.5)

plot_grid(p1, p2)

d10X.integrated@project.name="p11206"
plotCellCycle(d10X.integrated)

p3=DimPlot(d10X.integrated, reduction = "umap", split.by = "orig.ident")
print(p3)

dev.off()

save.image(cc("CHECKPOINT",DATE(),digest::digest(d10X),".Rdata"),compress=T)


DefaultAssay(d10X.integrated) <- "RNA"
clusters=levels(d10X.integrated@meta.data$seurat_clusters)

clusterMarkers=list()
for(ci in clusters) {
    clusterMarkers[[ci]]=NULL
    try({
        clusterMarkers[[ci]]=FindConservedMarkers(d10X.integrated, ident.1 = ci, grouping.var = "orig.ident")
    })
}

save.image(cc("CHECKPOINT",DATE(),digest::digest(d10X),".Rdata"),compress=T)

Idents(d10X.integrated)="orig.ident"

so=d10X.integrated
DefaultAssay(so)="RNA"
so=NormalizeData(so)
so=FindVariableFeatures(so)
so=ScaleData(so)

library(tidyverse)
library(openxlsx)

#markers.all=FindMarkers(so,ident.1="IBR10",ident.2="DMSO",pct=.1)

#save.image(cc("CHECKPOINT",DATE(),digest::digest(d10X),".Rdata"),compress=T)

#annote=readRDS("Homo_sapiens.GRCh37.75__GeneIDNameClass.rda")
#codingGenes=annote %>% filter(class=="protein_coding") %>% distinct(gene_name) %>% pull(gene_name)

# tbl.1=markers.all %>%
#     rownames_to_column("Gene") %>%
#     as_tibble %>%
#     mutate(avg_log2FC=log2(exp(1))*avg_logFC) %>%
#     mutate(FoldChange=ifelse(avg_log2FC>0,2^avg_log2FC,-2^(-avg_log2FC))) %>%
#     select(Gene,p_val_adj,FoldChange,avg_log2FC,pct.1,pct.2,p_val) %>%
#     mutate(OR.1=pct.1/(1-pct.1),OR.2=pct.2/(1-pct.2),lOR=log10(OR.1/OR.2))


# tbl.f=tbl.1 %>% arrange(desc(abs(FoldChange))) %>% filter(abs(FoldChange)>2 & p_val_adj<0.05)

    # filter(Gene %in% codingGenes) %>%
    # filter(abs(avg_log2FC)>log2(2) & p_val_adj<0.05)

# pct.clusters=so@meta.data %>% rownames_to_column("Cell") %>% as_tibble %>% count(orig.ident,seurat_clusters) %>% group_by(orig.ident) %>% mutate(PCT=n/sum(n)) %>% select(-n) %>% spread(seurat_clusters,PCT)

# write.xlsx(list(diffGenes=tbl.f,PCT.Cluster=pct.clusters),
#     "findMarkersDefault_pct0.1_IBR10_vs_DMSO__FC_2_FDR_0.05.xlsx")

# tbl.f %>% pull(Gene) -> gg
# plts=VlnPlot(so,features=gg,combine=F,pt.size=0,split.by="orig.ident",group.by="seurat_clusters")

# pg=list()
# for(ii in seq(ceiling(len(gg)/6))-1) {
#     if(ii<len(gg)) {
#         pg[[paste0("p",ii)]]=CombinePlots(plots=plts[((ii*6)+1):len(plts)],ncol=3,nrow=2,legend="right")
#     }
# }

if(!interactive()) {quit()}
stop("DDDDD")

pg.dot=DotPlot(so,features=rev(gg),cols = c("blue", "red"), dot.scale = 10,split.by="orig.ident",group.by="seurat_clusters")+RotatedAxis()

stop("DMSO")
pdf(file="diffGenes__pct0.1_IBR10_vs_DMSO__FC_2_FDR_0.05.pdf",height=8.5,width=11)
print(pg)
print(pg.dot)
dev.off()

cellTypes=paste0("c",d10X.integrated@meta.data$seurat_clusters,"_",d10X.integrated@meta.data$orig.ident)
cellTypes=factor(cellTypes,levels=sort(unique(cellTypes)))
d10X.integrated$cellType.Treat=cellTypes
Idents(d10X.integrated)="cellType.Treat"


so=d10X.integrated
DefaultAssay(so)="RNA"
so=NormalizeData(so)
so=FindVariableFeatures(so)
so=ScaleData(so)

markers=list()

for(ci in levels(d10X.integrated@meta.data$seurat_clusters)) {
    grp1=paste0("c",ci,"_IBR10")
stop("DMSO")
    grp2=paste0("c",ci,"_DMSO")
    cat(grp1,"vs",grp2,"...")
    markers[[paste0("c",ci)]]=FindMarkers(so,ident.1=grp1,ident.2=grp2,pct.min=.1,) %>%
        rownames_to_column("Gene") %>%
        as_tibble %>%
        mutate(avg_log2FC=log2(exp(1))*avg_logFC) %>%
        mutate(FoldChange=ifelse(avg_log2FC>0,2^avg_log2FC,-2^(-avg_log2FC))) %>%
        select(Gene,p_val_adj,FoldChange,avg_log2FC,pct.1,pct.2,p_val) %>%
        mutate(OR.1=pct.1/(1-pct.1),OR.2=pct.2/(1-pct.2),lOR=log10(OR.1/OR.2)) %>%
        arrange(p_val,desc(abs(FoldChange)))

    cat("\n")
}


tbl.2=bind_rows(markers,.id="Comp") %>%
    gather(Stat,Value,-Comp,-Gene) %>%
    unite(Comp.Stat,Comp,Stat,sep="__") %>%
    spread(Comp.Stat,Value) %>%
    mutate_at(vars(matches("p_val")),~replace(., is.na(.), 1))

fs=tbl.2 %>%
    select(Gene,matches("p_val$")) %>%
    mutate(FS=rowSums(-2*log(.[-1]+.Machine$double.xmin))) %>%
    arrange(desc(FS))

fs %>%
    head %>%
    pull(Gene) -> g1


pg=VlnPlot(so,features=g1)

print(pg)

save.image(cc("CHECKPOINT",DATE(),digest::digest(d10X),".Rdata"),compress=T)

stop("END OF SCRIPT")

stop("DMSO")
write.xlsx(tbl.1,"findMarkersDefault_pct0.1_IBR10_vs_DMSO__logFC_2_FDR_0.05.xlsx")

tbl.1 %>% head %>% pull(Gene) -> g1
head(tbl.1[tbl.1$FoldChange<0,]) %>% pull(Gene) -> gneg

pg=VlnPlot(d10X.integrated,features=g1)
pneg=VlnPlot(d10X.integrated,features=gneg)
stop("DMSO")
pdf(file="top6Genes___pct0.1_IBR10_vs_DMSO__logFC_2_FDR_0.05.pdf",width=11,height=8.5)
print(pg)
dev.off()

save.image(cc("CHECKPOINT.P2",DATE(),digest::digest(d10X),".Rdata"),compress=T)

if(0) {
    library(Seurat)
    library(cowplot)
    load("checkpointP2.Rdata")
}

so=d10X.integrated
so=CellCycleScoring(so,s.features=cc.genes$s.genes,g2m.features=cc.genes$g2m.genes,set.ident=T)
so=RunPCA(so, features=c(cc.genes$s.genes,cc.genes$g2m.genes))
PCAPlot(so,split="orig.ident")
DimPlot(so,split="orig.ident")


cellTypes=paste0(so@meta.data$Phase,"_",so@meta.data$orig.ident)
cellTypes=factor(cellTypes,levels=sort(unique(cellTypes)))
so$cellType.Phase=cellTypes
Idents(so)="cellType.Phase"

markers.0=markers
markers=list()


for(ci in sort(unique(so@meta.data$Phase))) {
    grp1=paste0(ci,"_IBR10")
stop("DMSO")
    grp2=paste0(ci,"_DMSO")
    cat(grp1,"vs",grp2,"...")
    markers[[ci]]=FindMarkers(so,ident.1=grp1,ident.2=grp2)
    cat("\n")
}


tbl.cc=map(markers,rownames_to_column,"Gene") %>%
    bind_rows(.id="Comp") %>%
    as_tibble %>%
    gather(Stat,Value,-Comp,-Gene) %>%
    unite(Comp.Stat,Comp,Stat,sep="__") %>%
    spread(Comp.Stat,Value) %>%
    mutate_at(vars(matches("p_val")),~replace(., is.na(.), 1)) %>%
    mutate(FS=-2*((log(G1__p_val)+log(G2M__p_val)+log(S__p_val)))) %>%
    arrange(desc(FS))

save.image(cc("CHECKPOINT.P2",DATE(),digest::digest(d10X),".Rdata"),compress=T)
