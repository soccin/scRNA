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

