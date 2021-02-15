suppressPackageStartupMessages(require(stringr))

usage="
usage: doSeuratV5.R [DEBUG=${DEBUG}] [MERGE=${MERGE}] [PROJNAME=${PROJNAME}] 10XDir1 [10XDir2 ... ]

    DEBUG        Set DEBUG mode (downsample to 10%)
    MERGE        If True then merge samples with simple merge
    PROJNAME     Set name of project. Must either be used or there
                 must be a file PROJNAME in folder with name

"

cArgs=commandArgs(trailing=T)
args=list(DEBUG=FALSE,MERGE=TRUE,PROJNAME="scRNA")
usage=str_interp(usage,args)

ii=grep("=",cArgs)
if(len(ii)>0) {
    parseArgs=str_match(cArgs[ii],"(.*)=(.*)")
    aa=apply(parseArgs,1,function(x){args[[str_trim(x[2])]]<<-str_trim(x[3])})
}

args$DEBUG=as.logical(args$DEBUG)
args$MERGE=as.logical(args$MERGE)

if(args$PROJNAME=="scRNA") {
    if(file.exists("PROJNAME")) {
        args$PROJNAME=scan("PROJNAME","",quiet=T)
    } else {
        cat(usage)
        quit()
    }
}

argv=grep("=",cArgs,value=T,invert=T)

##############################################################################
if(len(argv)<1) {
    cat(usage)
    quit()
}
cat("\n=========================================================\n")
args[["10XDirs"]]=paste0(argv,collapse=", ")
cat(str(args))
cat("\n")
##############################################################################

suppressPackageStartupMessages({
    library(Seurat)
    library(cowplot)
    library(patchwork)
    library(tidyverse)
    library(gridExtra)
})

source("seuratTools.R")

dataFolders=argv
sampleIDs=gsub("_",".",gsub(".outs.*","",gsub(".*/s_","",dataFolders)))
names(dataFolders)=sampleIDs

d10X=list()

projectName="p11533"

for(ii in seq(len(dataFolders))) {
    sampleName=sampleIDs[ii]
    cat("Reading Sample =",sampleName,"...")
    d10X[[sampleName]] <- read10XDataFolderAsSeuratObj(dataFolders[ii],projectName)
    cat("\n")
}
d10X.orig=d10X

if(args$MERGE & len(d10X)>1) {
    cat("\nMerging sample files...")
    s.merge=merge(d10X[[1]],d10X[-1],project=args$PROJNAME)
    d10X=list()
    d10X[[args$PROJNAME]]=s.merge
    cat("done\n\n")
}

glb.digest=digest::digest(d10X)
cat("digest=",digest::digest(d10X),"\n")

## QC
PCT_MITO=10
MIN_FEATURE_RNA=1000
MIN_NCOUNT_RNA=2500

algoParams=list()
algoParams$PCT_MITO=PCT_MITO
algoParams$MIN_FEATURE_RNA=MIN_FEATURE_RNA
algoParams$MIN_NCOUNT_RNA=MIN_NCOUNT_RNA

# as_tibble(algoParams) %>% gather(param,Value)

stats=list()

doQCandFilter <- function(so) {


    sampleId=unique(so@meta.data$orig.ident)
    if(len(sampleId)>1) {
        sampleId=cc("MERGE",so@project.name)
    }

    pg0=VlnPlot(so, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

    max.nCount_RNA=max(so@meta.data$nCount_RNA)
    max.nFeature_RNA=max(so@meta.data$nFeature_RNA)
    max.pct.mito=max(so@meta.data$percent.mt)

    tbl=so@meta.data %>% tibble %>%
        count(
            Count.RNA=nCount_RNA>MIN_NCOUNT_RNA,
            Num.Features=nFeature_RNA>MIN_FEATURE_RNA,
            PCT.MT=percent.mt<PCT_MITO
            ) %>%
        mutate(PCT=round(100*n/sum(n),1))

    tbl2=so@meta.data %>% tibble %>%
        mutate(
            nCount_RNA=nCount_RNA>MIN_NCOUNT_RNA,
            nFeature_RNA=nFeature_RNA>MIN_FEATURE_RNA,
            percent.mt=percent.mt<PCT_MITO
            ) %>%
        summarize_if(is.logical,~round(100*(1-sum(.)/n()),1)) %>%
        gather(Metric,PCT.Fail)

    plot1 <- FeatureScatter(so, feature1 = "nCount_RNA", feature2 = "percent.mt") +
                geom_hline(yintercept=PCT_MITO,col="grey",alpha=0.75) +
                geom_vline(xintercept=MIN_NCOUNT_RNA,col="grey",alpha=0.75) +
                annotation_custom(tableGrob(tbl,rows=NULL),
                    xmin=max.nCount_RNA/2,ymin=50,xmax=max.nCount_RNA/2,ymax=50) +
                scale_x_continuous(breaks=sort(c(MIN_NCOUNT_RNA,seq(0,1e6,by=20000)))) +
                scale_y_continuous(breaks=sort(c(PCT_MITO,seq(0,100,by=20))))

    plot2 <- FeatureScatter(so, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
                geom_hline(yintercept=MIN_FEATURE_RNA,col="grey",alpha=0.75) +
                geom_vline(xintercept=MIN_NCOUNT_RNA,col="grey",alpha=0.75) +
                scale_x_continuous(breaks=sort(c(MIN_NCOUNT_RNA,seq(0,1e6,by=20000)))) +
                scale_y_continuous(breaks=sort(c(MIN_FEATURE_RNA,seq(0,1e5,by=2000)))) +
                annotation_custom(tableGrob(tbl2,rows=NULL),xmin=max.nCount_RNA/3,ymin=max.nFeature_RNA/10)
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

    pdf(file=cc("seuratQC",sampleId,"01.pdf"),height=8.5,width=11)
    cat(cc("seuratQC",sampleId,"01.pdf"),"\n")
    print(pg0)
    print(plot1)
    print(plot2)
    dev.off()

    stats[[sampleId]]=keep

    so

}

cat("\nDoQCandFilter\n")
for(ii in seq(d10X)) {
    print(ii)
    d10X[[ii]]=doQCandFilter(d10X[[ii]])
}

if(args$DEBUG) {

    c("\n\nDEBUG SET; Subset data\n\n")

    set.seed(101)

    cat("\nDEBUG::subset\n")

    for(ii in seq(d10X)) {
        print(ii)
        xx=d10X[[ii]]
        d10X[[ii]]=subset(xx,cells=Cells(xx)[runif(nrow(xx@meta.data))<.1])
    }

    cat("digest=",digest::digest(d10X),"\n")

}

cat("\nScoreCellCycle\n")
for(ii in seq(d10X)) {
    print(ii)
    d10X[[ii]]=scoreCellCycle(d10X[[ii]])
}

pcc=list()
cat("\nPlotCellCycle\n")
for(ii in seq(d10X)) {
    pcc[[ii]]=plotCellCycle(preProcessSO(d10X[[ii]]))
}

pdf(file="seuratQC_CellCycle.pdf",width=11,height=11)

if(len(pcc)>1) {
    nPages=ceiling(len(pcc)/4)
    for(ii in seq(nPages)) {
        jj=(1:4)+4*(ii-1)
        jj=intersect(seq(pcc),jj)
        pOut=pcc[jj]
        if(len(pOut)<4){
            pOut=c(pOut,rep(list(ggplot()+theme_void()),4-len(jj)))
        }
        print(wrap_plots(pOut,ncol=2))
    }
} else {
    print(pcc[[1]])
}

dev.off()

save.image(cc("CHECKPOINT",DATE(),glb.digest,".Rdata"),compress=T)


if(len(d10X)>1) {
    cat("\n\n The rest of this workflow only works on merged datasets\n\n")
    quit()
}


stop("Continue workign")

#
# SCTransform Normalizes
#

cat("\n\n  SCTransform and normalize\n")

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

save.image(cc("CHECKPOINT",DATE(),glb.digest,".Rdata"),compress=T)

##
## Find Cluster Markers
##

cat("\n\n  Find Cluster Markers\n")

so=d10X.integrated
DefaultAssay(so)="RNA"
so=NormalizeData(so)
so=FindVariableFeatures(so)
so=ScaleData(so)

so@meta.data$Types=gsub("\\d$","",so@meta.data$orig.ident)

clusterMarkers=FindAllMarkers(so,only.pos=TRUE,logfc.threshold=0.25,min.pct = 0.25)

require(dplyr)

top6ClusterMarkers=clusterMarkers %>%
    tibble %>%
    group_by(cluster) %>%
    arrange(cluster,p_val_adj) %>%
    mutate(N=row_number()) %>%
    filter(N<=6)

png(filename=cc("seuratAnalysis_ClusterMarkers",projectTag,"%03d.png"),
    type="cairo",
    units="in",
    width=11,
    height=8.5,
    pointsize=12,
    res=96)

#pdf(file=cc("seuratAnalysis_ClusterMarkers",projectTag,".pdf"),width=11,height=8.5)

for(cii in levels(clusterMarkers$cluster)) {
    print(cii)
    gene.ii=top6ClusterMarkers %>% filter(cluster==cii) %>% pull(gene)
    p1=VlnPlot(so,features=gene.ii,pt.size=.25)
    p2=DotPlot(so,features=gene.ii,split.by="Types",cols=RColorBrewer::brewer.pal(3,"Dark2"))
    print(p1)
    print(p2)
}

dev.off()

save.image(cc("CHECKPOINT",DATE(),digest::digest(d10X),".Rdata"),compress=T)

saveRDS(d10X.integrated,"obj__d10X.integrated.rda",compress=T)
saveRDS(d10X.integrated@meta.data,"obj__d10X.integrated_meta.data.rda",compress=T)
