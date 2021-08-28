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

if(Sys.getenv("SDIR")=="") {
    SDIR="."
} else {
    SDIR=Sys.getenv("SDIR")
}

##############################################################################
cat("\n=========================================================\n")
cat(str(args))
cat("\n")
##############################################################################

suppressPackageStartupMessages({
    library(Seurat)
    library(patchwork)
    library(tidyverse)
    library(openxlsx)
})

source(file.path(SDIR,"seuratTools.R"))
source(file.path(SDIR,"plotTools.R"))
source(file.path(SDIR,"doQCandFilter.R"))

glbs=args$glbs
ap=args$algoParams

plotNo<-makeAutoIncrementor(10)

##########################################################################
#
# INCLUDE BREAK
#
##########################################################################


d10X=readRDS(args$PASS1.RDAFile)

glb.digest=digest::digest(d10X)
cat("digest=",digest::digest(d10X),"\n")

cat("\nDoQCandFilter\n")
for(ii in seq(d10X)) {
    print(ii)
    ret=doQCandFilter(d10X[[ii]], ap$MIN_NCOUNT_RNA, ap$MIN_FEATURE_RNA, ap$PCT_MITO)
    d10X[[ii]]=ret$so
}


##############################################################################
# Check Cell Cycle
#
cat("\nScoreCellCycle\n")

for(ii in seq(d10X)) {
    print(ii)
    d10X[[ii]]=scoreCellCycle(d10X[[ii]])
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



##############################################################################
# Do SCTransform
#
cat("\nSCTransform\n")

d10X.int=list()
for(ii in seq(d10X)) {
    print(ii)
    d10X.int[[ii]]=SCTransform(d10X[[ii]],vars.to.regress = c('S.Score', 'G2M.Score'))
}

##############################################################################
# Do Integration0
#
cat("\nSCTransform\n")

features <- SelectIntegrationFeatures(object.list = d10X.int, nfeatures = 3000)
d10X.int <- PrepSCTIntegration(object.list = d10X.int, anchor.features = features)

anchors <- FindIntegrationAnchors(
        object.list = d10X.int,
        normalization.method = "SCT",
        anchor.features = features
        )

d10X.integrate <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

d10X.integrate <- RunPCA(d10X.integrate, verbose = FALSE)
d10X.integrate <- RunUMAP(d10X.integrate, reduction = "pca", dims = 1:30)

DimPlot(d10X.integrate,reduction="umap") + scale_color_brewer(palette="Paired")

cellCycle.genes = getCellCycleGenes(glbs$genome)

so=CellCycleScoring(so,
                    s.features=cellCycle.genes$s.genes,
                    g2m.features=cellCycle.genes$g2m.genes,
                    set.ident=T
                    )

cc.meta.data=so@meta.data[,c("S.Score","G2M.Score","Phase")]

d10X.integrate=AddMetaData(d10X.integrate,cc.meta.data$Phase,"Phase")
d10X.integrate=AddMetaData(d10X.integrate,cc.meta.data$G2M.Score,"G2M.Score")
d10X.integrate=AddMetaData(d10X.integrate,cc.meta.data$S.Score,"S.Score")
d10X.integrate=AddMetaData(d10X.integrate,cc.meta.data$S.Score-cc.meta.data$G2M.Score,"CC.Difference")

Idents(d10X.integrate)="Phase"

pdf(file=cc("seuratQC",args$PROJNAME,plotNo(),"PostIntegrateCC.pdf"),width=11,height=8.5)
plotCellCycle(d10X.integrate,"Post Integration CC Regression")
DimPlot(d10X.integrate,reduction="umap",group.by="Phase")
dev.off()

halt("BREAK-POINT 149")

so=d10X.integrate

pv1=VariableFeaturePlot(so)
top10 <- head(VariableFeatures(so), 10)
pv2 <- LabelPoints(plot = pv1, points = top10, repel = TRUE, xnudge=0, ynudge=0)
pdf(file=cc("seuratQC",args$PROJNAME,plotNo(),"VariableFeatures.pdf"),width=11,height=8.5)
print(pv2)
dev.off()

# PCA

# https://satijalab.org/seurat/archive/v3.0/s13k_tutorial.html
# Perform linear dimensional reduction

s1=so
s1=RunPCA(s1,features=VariableFeatures(s1),approx=FALSE)

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
pdf(file=cc("seuratQC",args$PROJNAME,plotNo(),"PCADimMetric.pdf"),width=11,height=8.5)
print(p.js1)
print(p.elbow)
dev.off()

# stop("CONTINUE")
# if(interactive()) {s1=readRDS("jsSampling_nReps_100_Dims_50_.rda")}

#
# Look at PCADimMetric.pdf to decided
#

#stop("\n\n CHECK PCA AND CONTINUE\n\n")

ap$NDIMS=40
nDims=40

ap$ClusterResolutions=c(0.1,0.2,0.5,0.8)

s1 <- FindNeighbors(s1, dims = 1:nDims)
s1 <- FindClusters(s1, resolution = ap$ClusterResolutions)
s1 <- RunUMAP(s1, dims = 1:nDims)

library(pals)
maxClusters=s1@meta.data %>% tibble %>% distinct(RNA_snn_res.0.8) %>% pull %>% len
if(maxClusters>33) {
    save.image(cc("CHECKPOINT",DATE(),".RData"),compress=T)
    stop("\n\nTOO MANY CLUSTERS\n\n")
}

pal1=c(cols25(maxClusters),brewer.dark2(8))
pu=list()
pu[[1]]=DimPlot(s1, reduction = "umap", label=T, group.by="RNA_snn_res.0.1", label.size=6) + scale_color_manual(values=pal1) + ggtitle("RNA_snn_res.0.1")
pu[[2]]=DimPlot(s1, reduction = "umap", label=T, group.by="RNA_snn_res.0.2", label.size=6) + scale_color_manual(values=pal1) + ggtitle("RNA_snn_res.0.2")
pu[[3]]=DimPlot(s1, reduction = "umap", label=T, group.by="RNA_snn_res.0.5", label.size=6) + scale_color_manual(values=pal1) + ggtitle("RNA_snn_res.0.5")
pu[[4]]=DimPlot(s1, reduction = "umap", group.by="orig.ident") + scale_color_brewer(palette="Paired")
pu[[5]]=DimPlot(s1, reduction = "umap", group.by="Phase")

pdf(file=cc("seuratQC",args$PROJNAME,plotNo(),"UMAP",nDims,".pdf"),width=11,height=8.5)
print(pu)
dev.off()

clusterRes="RNA_snn_res.0.5"
#s1=SetIdent(s1,value="RNA_snn_res.0.2")
s1=SetIdent(s1,value=clusterRes)

clusterMarkers=FindAllMarkers(s1,only.pos=TRUE,logfc.threshold=0.25,min.pct = 0.25)

pct=c(clusterMarkers$pct.1,clusterMarkers$pct.2)
ps=min(min(pct[pct>0]),min(1-pct[pct<1]))/2

FDR.cut=0.05
logFC.cut=1
filterCLTable<-function(clm) {
    tibble(clm) %>%
        mutate(lor.1=log((pct.1+ps)/(1-pct.1+ps))) %>%
        mutate(lor.2=log((pct.2+ps)/(1-pct.2+ps))) %>%
        mutate(lOR=lor.1-lor.2) %>%
        arrange(desc(avg_logFC))
}

#

cl=filterCLTable(clusterMarkers) %>%
    dplyr::select(cluster,gene,p_val_adj,avg_logFC,pct.1,pct.2,lOR) %>%
    filter(p_val_adj<FDR.cut & avg_logFC>logFC.cut)

clusters=cl %>% distinct(cluster) %>% arrange(cluster) %>% pull %>% levels


clusterMarkerTbl=list()
for(ci in clusters) {
    clusterMarkerTbl[[ci]]=cl %>% filter(cluster==ci)
}

geneCounts=cl %>% count(cluster)

ll=c(list(GeneCounts=geneCounts,AllCluster=cl),clusterMarkerTbl)

xfile=cc("tblClusterMarkers","",clusterRes,"FDR",FDR.cut,"logFC",logFC.cut)

write.xlsx(ll,paste0(xfile,".xlsx"))

args$algoParams=ap
obj=list(
    args=args,
    s1=s1,
    clusterMarkers=clusterMarkers,
    cl=cl,
    clusterRes=clusterRes
    )

args.digest.orig=digest::digest(obj)
args$PASS2.RDAFile=cc("pass_02",args.digest.orig,"OBJ",".rda")
obj$args=args

saveRDS(obj,args$PASS2.RDAFile,compress=T)
write_yaml(args,cc("pass_02","PARAMS.yaml"))


plt.cmark=cl %>% group_split(cluster) %>% map(plotClusterMarkers,s1,pal1)

cmFile=cc("seuratQC",args$PROJNAME,plotNo(),"ClusterMarkers","%03d",".png")

png(filename=cmFile,
    type="cairo",
    units="in",
    width=14,
    height=8.5,
    pointsize=12,
    res=96)

print(plt.cmark)

dev.off()

g1=cl %>% distinct(cluster,.keep_all=T) %>% pull(gene)
g2=cl %>% filter(!(gene %in% g1)) %>% distinct(cluster,.keep_all=T) %>% pull(gene)

dot.gene.lOR=cl %>%
    filter(gene %in% c(g1,g2)) %>%
    arrange(desc(lOR)) %>%
    distinct(gene,.keep_all=T) %>%
    head(12) %>%
    pull(gene)

dot.gene.lFC=cl %>%
    filter(gene %in% c(g1,g2)) %>%
    arrange(desc(avg_logFC)) %>%
    distinct(gene,.keep_all=T) %>%
    head(12) %>%
    pull(gene)

pc1=DotPlot(s1,features=dot.gene.lOR)
pc2=DotPlot(s1,features=dot.gene.lFC)

pdf(file=cc("seuratQC",args$PROJNAME,plotNo(),"ClusterMarkersDot",".pdf"),width=11,height=8.5)
print(pc1)
print(pc2)
dev.off()
