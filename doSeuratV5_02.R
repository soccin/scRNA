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

library(yaml)
args=read_yaml(cArgs[1])

##############################################################################
cat("\n=========================================================\n")
cat(str(args))
cat("\n")
##############################################################################

suppressPackageStartupMessages({
    library(Seurat)
    library(patchwork)
    library(tidyverse)
})

source("seuratTools.R")
source("doQCandFilter.R")

glbs=args$glbs
ap=args$algoParams

plotNo<-makeAutoIncrementor(10)

d10X=readRDS(args$PASS1.RDAFile)

if(args$MERGE & len(d10X)>1) {
    cat("\nMerging sample files...")
    s.merge=merge(d10X[[1]],d10X[-1],project=args$PROJNAME)
    d10X=list()
    d10X[[args$PROJNAME]]=s.merge
    cat("done\n\n")
}

glb.digest=digest::digest(d10X)
cat("digest=",digest::digest(d10X),"\n")

cat("\nDoQCandFilter\n")
for(ii in seq(d10X)) {
    print(ii)
    ret=doQCandFilter(d10X[[ii]], ap$MIN_NCOUNT_RNA, ap$MIN_FEATURE_RNA, ap$PCT_MITO)
    d10X[[ii]]=ret$so
}

if(args$DEBUG) {

    c("\n\nDEBUG SET; Subset data\n\n")
    set.seed(ap$SEED)
    cat("\nDEBUG::subset\n")
    for(ii in seq(d10X)) {
        print(ii)
        xx=d10X[[ii]]
        d10X[[ii]]=subset(xx,cells=Cells(xx)[runif(nrow(xx@meta.data))<args$DOWNSAMPLE])
    }
    cat("digest=",digest::digest(d10X),"\n")

}

if(len(d10X)>1) {
    cat("\n\n The rest of this workflow only works on merged datasets\n\n")
    quit()
}

##
# Cell-Cycle Scoring and Regression
# https://satijalab.org/seurat/archive/v3.1/cell_cycle_vignette.html
#

ap$NFEATURES=2000

so=d10X[[1]]

cellCycle.genes = getCellCycleGenes(glbs$genome)

so <- NormalizeData(so)
so <- FindVariableFeatures(so, selection.method="vst", nfeatures = ap$NFEATURES)


pv1=VariableFeaturePlot(so)
top10 <- head(VariableFeatures(so), 10)
pv2 <- LabelPoints(plot = pv1, points = top10, repel = TRUE, xnudge=0, ynudge=0)
pdf(file=cc("seuratQC",plotNo(),"VariableFeatures.pdf"),width=11,height=8.5)
print(pv2)
dev.off()

so <- ScaleData(so, features=rownames(so))

cellCycle.genes = getCellCycleGenes(glbs$genome)

so=CellCycleScoring(so,
                    s.features=cellCycle.genes$s.genes,
                    g2m.features=cellCycle.genes$g2m.genes,
                    set.ident=T
                    )

if(interactive()) {
    #top("DDD");
    #s1=readRDS("ccRegression_nFeat_5000_.rda");ap$NFEATURES=len(VariableFeatures(so));stop("BREAK")
}

s1=ScaleData(so, vars.to.regress = c("S.Score", "G2M.Score"), features = VariableFeatures(so))

#
# After Cell Cycle regression Ident -> Phase reset to orig.ident
#
s1 <- SetIdent(s1,value="orig.ident")


saveRDS(s1,cc("ccRegression","nFeat",len(VariableFeatures(so)),".rda"),compress=T)

pdf(file=cc("seuratQC",plotNo(),"PostCCRegress.pdf"),width=11,height=8.5)

plotCellCycle(s1,"Post CC Regression")

dev.off()

# PCA

# https://satijalab.org/seurat/archive/v3.0/s13k_tutorial.html
# Perform linear dimensional reduction

s1=RunPCA(s1,features=VariableFeatures(s1))

# # Determine the ‘dimensionality’ of the dataset

# To overcome the extensive technical noise in any single feature for scRNA-seq data,
# Seurat clusters cells based on their PCA scores, with each PC essentially representing
# a ‘metafeature’ that combines information across a correlated feature set. The top
# 3 principal components therefore represent a robust compression of the dataset. However,
# how many componenets should we choose to include? 10? 20? 100?

# In Macosko et al, (http://www.cell.com/abstract/S0092-8674(15)00549-8)
# we implemented a resampling test inspired by the JackStraw procedure. We randomly permute
# a subset of the data (1% by default) and rerun PCA, constructing a ‘null distribution’ of
# feature scores, and repeat this procedure. We identify ‘significant’ PCs as those who
# have a strong enrichment of low p-value features.

nReps=100
nDims=50 # Default
s1 <- JackStraw(s1, num.replicate = nReps, dims=nDims)

saveRDS(s1,cc("jsSampling","nReps",nReps,"Dims",nDims,".rda"),compress=T)

s1 <- ScoreJackStraw(s1, dims = 1:nDims)

p.js1=JackStrawPlot(s1, dims = 1:nDims)
p.elbow=ElbowPlot(s1,ndims=nDims) + geom_hline(yintercept=0,color="grey",size=2)
pdf(file=cc("seuratQC",plotNo(),"PCADimMetric.pdf"),width=11,height=8.5)
print(p.js1)
print(p.elbow)
dev.off()

stop("CONTINUE")
nDims=40

s1 <- FindNeighbors(s1, dims = 1:nDims)
s1 <- FindClusters(s1, resolution = c(0.1,0.2,0.5,0.8))
s1 <- RunUMAP(s1, dims = 1:nDims)


pum.1=DimPlot(s1, reduction = "umap",label=T)
pum.2=DimPlot(s1, reduction = "umap", group.by="orig.ident")
pum.3=DimPlot(s1, reduction = "umap", group.by="Phase")
pnull=ggplot()+theme_void()

pdf(file=cc("seuratQC",plotNo(),"UMAP",nDims,".pdf"),width=11,height=8.5)
print((pum.1+pum.2)/(pum.3+pnull))
dev.off()

pum.1=DimPlot(fcs1, reduction = "umap", label=T)
pum.2=DimPlot(fcs1, reduction = "umap", group.by="orig.ident")
pum.3=DimPlot(fcs1, reduction = "umap", group.by="Phase")
print((pum.1+pum.2)/(pum.3+pnull))

#s1@meta.data %>% tibble %>% count(seurat_clusters,Phase) %>% group_by(Phase) %>% mutate(PCT=n/sum(n)) %>% select(-n) %>% spread(Phase,PCT)

s1 <- RunUMAP(s1, dims = 1:nDims, spread=1, min.dist=0.6);
p1=DimPlot(s1, reduction = "umap", label=T) + ggtitle(paste("Spread=1 min.dist=0.6"));
s1 <- RunUMAP(s1, dims = 1:nDims, spread=1, min.dist=0.3);
p2=DimPlot(s1, reduction = "umap", label=T) + ggtitle(paste("Spread=1 min.dist=0.3"));
s1 <- RunUMAP(s1, dims = 1:nDims, spread=1, min.dist=0.1);
p3=DimPlot(s1, reduction = "umap", label=T) + ggtitle(paste("Spread=1 min.dist=0.1"));
s1 <- RunUMAP(s1, dims = 1:nDims, spread=1, min.dist=0.05);
p4=DimPlot(s1, reduction = "umap", label=T) + ggtitle(paste("Spread=1 min.dist=0.05"));
