
qcPlotNum = 0
nextQCPlotFileName <- function(...) {
    qcPlotNum <<- qcPlotNum + 1
    paste0("qc",sprintf("%02d",qcPlotNum),cc("Seurat",TAGS,...),".pdf")
}


if(interactive()) {
    pdf.orig<-pdf
    pdf<-function(...) {x11()}
    dev.off.orig<-dev.off
    dev.off<-function(...){}
}

PDIR="scRNA"
if(PDIR=="") stop("Set SDIR")

#########################################################################################
#########################################################################################
library(tidyverse)
library(Seurat)
#########################################################################################
#stop("INCLUDE")
#########################################################################################

TAGS="PlateReg"
PROJNO="10223_B"
countsFile="pipeline/gene/counts_gene/Proj_10223_B_htseq_all_samples.txt"
metaDataFile="metadata.csv"
genesOfInterest=c("ID1","ID3")

#
# Primary QC filtering paramters
#

PCT_MITO=25
MIN_FEATURE_RNA=1000
MIN_NCOUNT_RNA=0.25e6

dx=read_tsv(countsFile, col_types = cols())

ds <- dx %>%
    filter(!is.na(GeneSymbol)) %>%
    select(-GeneSymbol) %>%
    data.frame %>%
    column_to_rownames("GeneID") %>%
    as.matrix

annote =  dx %>% filter(!is.na(GeneSymbol)) %>% select(1,2)
mitoGenes=scan(file.path(PDIR,"mitoGenes.id.txt"),"")
annote=annote %>% mutate(MitoGene=gsub("\\.\\d+$","",GeneID) %in% mitoGenes)

rownames(ds)=scater::uniquifyFeatureNames(annote$GeneID,annote$GeneSymbol)

metadata=read_csv(metaDataFile)

gPos=data.frame(Pos=t(ifelse(ds[genesOfInterest,]>0,"Pos","Neg")))
gPos=gPos %>% rownames_to_column("Sample") %>% as_tibble
metadata=full_join(gPos,metadata)

#
# Seurat needs metadata as data from with Samples in rownames
#
metadata=data.frame(metadata) %>% column_to_rownames("Sample")

if(!all(colnames(ds)==rownames(metadata))) {
    cat("\n    ERROR cols/rows of data & metadata do not match up\n\n")
    stop("FATAL:ERROR")
}

so=CreateSeuratObject(counts=ds,project=PROJNO,min.cells=3,min.features=200,meta.data=metadata)

so[["percent.mt"]]=PercentageFeatureSet(so, pattern = "^MT-")

pdf(nextQCPlotFileName(),width=11,height=8.5)

VlnPlot(so, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, , group.by="Plate")

plot1 <- FeatureScatter(so, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by="Plate")
plot2 <- FeatureScatter(so, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by="Plate")
CombinePlots(plots = list(plot1, plot2))

dev.off()

keep=(
        so@meta.data$nFeature_RNA > MIN_FEATURE_RNA &
        so@meta.data$nCount_RNA > MIN_NCOUNT_RNA &
        so@meta.data$percent.mt < PCT_MITO
        )

knitr::kable(table(keep,so@meta.data$Plate))

so.orig=so
so <- subset(so,
        subset = nFeature_RNA > MIN_FEATURE_RNA & nCount_RNA > MIN_NCOUNT_RNA & percent.mt < PCT_MITO
        )
#so@meta.data$Plate=droplevels(so@meta.data$Plate)

so <- NormalizeData(so)

# ## Identification of highly variable features (feature selection)
# We next calculate a subset of features that exhibit high cell-to-cell
# variation in the dataset (i.e, they are highly expressed in some cells,
# and lowly expressed in others). We and others have found that focusing
# on these genes in downstream analysis helps to highlight biological
# signal in single-cell datasets.

so <- FindVariableFeatures(so, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(so), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(so)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

pdf(nextQCPlotFileName("VarianceTopGenes"),width=11,height=8.5)
print(plot2)
dev.off()

#CombinePlots(plots = list(plot1, plot2))

# ## Scaling the data
# Next, we apply a linear transformation (‘scaling’) that is a standard
# pre-processing step prior to dimensional reduction techniques like PCA.
# The ScaleData function:
#
# - Shifts the expression of each gene, so that the mean expression across
#   cells is 0
# - Scales the expression of each gene, so that the variance across cells is 1
#    - This step gives equal weight in downstream analyses, so that
#      highly-expressed genes do not dominate
# - The results of this are stored in `so[["RNA"]]@scale.data`

all.genes <- rownames(so)
so.orig=so

#so <- ScaleData(so, features = all.genes)
so <- ScaleData(so, features = all.genes, vars.to.regress = c("Plate"))

#
## Cell Cycle Analysis
#


cc.genes=cc.genes.updated.2019

so=CellCycleScoring(so,s.features=cc.genes$s.genes,g2m.features=cc.genes$g2m.genes,set.ident=T)
pdf(nextQCPlotFileName("CellCycle"),width=11,height=8.5)
RidgePlot(so, features=c("PCNA","CDC20","AURKA"),ncol=2)
so=RunPCA(so, features=c(cc.genes$s.genes,cc.genes$g2m.genes))
DimPlot(so) + ggtitle("Pre Cell Cycle Regression")

#DO_CELL_CYCLE=FALSE
if(!exists("DO_CELL_CYCLE")) {
    stop("NEED TO SET DO_CELL_CYCLE")
}

if(DO_CELL_CYCLE) {
    #
    # Regress out cell cycle scores during data scaling also remove Plate covariate
    #

    #so <- ScaleData(so, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(so))
    so <- ScaleData(so, vars.to.regress = c("S.Score", "G2M.Score", "Plate"), features = rownames(so))
    so <- RunPCA(so, features=c(cc.genes$s.genes,cc.genes$g2m.genes))
    DimPlot(so) + ggtitle("Post Cell Cycle Regression")

}

dev.off()

so <- RunPCA(so, features = VariableFeatures(so), nfeatures.print = 10)
VizDimLoadings(so, dims = 1:2, reduction = "pca")

pdf(nextQCPlotFileName("PCA"),width=11,height=8.5)
DimPlot(so, reduction = "pca", group.by="Plate",pt.size = 3)

DimHeatmap(so, dims = 1, cells = 500, balanced = TRUE)
dev.off()

so <- JackStraw(so, num.replicate = 100)
so <- ScoreJackStraw(so, dims = 1:20)

pdf(nextQCPlotFileName("Dim"),width=11,height=8.5)
JackStrawPlot(so, dims = 1:10)
ElbowPlot(so)
dev.off()

JS(so[["pca"]])$overall.p.values
which(p.adjust(JS(so[["pca"]])$overall.p.values[,2])<0.05)

so.orig=so

#dimsToUse=1:7

if(!exists("dimsToUse")) {
    cat("\n\n    NEED TO SET dimsToUse\n\n\n")
}


so=FindNeighbors(so,dims=dimsToUse)
so=FindClusters(so,resolution=0.5)
so=RunUMAP(so,dims=dimsToUse)
DimPlot(so,reduction="umap",ncol=2,pt.size=3,group.by="seurat_clusters")


pdf(nextQCPlotFileName("UMAP"),width=11,height=8.5)
DimPlot(so,reduction="umap",ncol=2,pt.size=3,group.by="seurat_clusters") +
    ggtitle(paste("By Clusters | Dimensions =",dimsToUse[1],"...",dimsToUse[len(dimsToUse)]))

DimPlot(so,reduction="umap",ncol=2,pt.size=3,group.by="Phase") +
    ggtitle(paste("By Cell Cycle | Dimensions =",dimsToUse[1],"...",dimsToUse[len(dimsToUse)]))

DimPlot(so,reduction="umap",ncol=2,pt.size=3,group.by="Plate") +
    ggtitle(paste("By Plate | Dimensions =",dimsToUse[1],"...",dimsToUse[len(dimsToUse)]))

DimPlot(so,reduction="umap",ncol=2,pt.size=3,group.by="Pos.ID1") +
    ggtitle("ID1")

DimPlot(so,reduction="umap",ncol=2,pt.size=3,group.by="Pos.ID3") +
    ggtitle("ID3")

for(gi in genesOfInterest) {
    FeaturePlot(so,c(gi),pt.size=3)
}

FeaturePlot(so,genesOfInterest,pt.size=3)
VlnPlot(so,genesOfInterest)

dev.off()

so.markers <- FindAllMarkers(so, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
so.markers %>% group_by(cluster) %>% top_n(n = 4, wt = avg_logFC)

# xx=so.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
# so.markers=xx[xx$p_val_adj < 0.05,]

marker.genes=so.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_logFC) %>% pull(gene)

pdf(nextQCPlotFileName("Markers"),width=8.5,height=11)
VlnPlot(so,marker.genes,ncol=2)
FeaturePlot(so,marker.genes,ncol=2)
dev.off()

