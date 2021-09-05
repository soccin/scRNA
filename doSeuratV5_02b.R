suppressPackageStartupMessages(require(stringr))

usage="
usage: doSeuratV5_02b.R PARAMS.yaml

    PARAMS.yaml     parameter file from pass1

"

cArgs=commandArgs(trailing=T)

if(len(cArgs)!=1) {
    cat(usage)
    quit()
}

if(R.Version()$major<4) {
    cat(usage)
    cat("\n\nThis script needs version(R).major>=4\n\n")
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
})

source(file.path(SDIR,"seuratTools.R"))
source(file.path(SDIR,"plotTools.R"))

glbs=args$glbs
ap=args$algoParams

plotNo<-makeAutoIncrementor(20)

##########################################################################
#
# INCLUDE BREAK
#
##########################################################################


###d10X=readRDS(args$PASS1.RDAFile)

glb.digest=digest::digest(d10X)
cat("digest=",digest::digest(d10X),"\n")

##########################################################################
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

nDims=50 # Default

p.elbow=ElbowPlot(s1,ndims=nDims) + geom_hline(yintercept=0,color="grey",size=2)
pdf(file=cc("seuratQC",args$PROJNAME,plotNo(),"PCADimMetric.pdf"),width=11,height=8.5)
print(p.elbow)
dev.off()

#stop("\n\n CHECK PCA AND CONTINUE\n\n")

ap$NDIMS=20
nDims=20

ap$ClusterResolutions=c(0.1,0.2,0.5,0.8)

s1 <- FindNeighbors(s1, dims = 1:nDims)
s1 <- FindClusters(s1, resolution = ap$ClusterResolutions)
s1 <- RunUMAP(s1, dims = 1:nDims)

library(pals)
maxClusters=s1@meta.data %>% tibble %>% distinct(integrated_snn_res.0.8) %>% pull %>% len
if(maxClusters>33) {
    save.image(cc("CHECKPOINT",DATE(),".RData"),compress=T)
    stop("\n\nTOO MANY CLUSTERS\n\n")
}

pal1=c(cols25(maxClusters),brewer.dark2(8))
pu=list()
pu[[1]]=DimPlot(s1, reduction = "umap", label=T, group.by="integrated_snn_res.0.1", label.size=6) + scale_color_manual(values=pal1) + ggtitle("integrated_snn_res.0.1")
pu[[2]]=DimPlot(s1, reduction = "umap", label=T, group.by="integrated_snn_res.0.2", label.size=6) + scale_color_manual(values=pal1) + ggtitle("integrated_snn_res.0.2")
pu[[3]]=DimPlot(s1, reduction = "umap", label=T, group.by="integrated_snn_res.0.5", label.size=6) + scale_color_manual(values=pal1) + ggtitle("integrated_snn_res.0.5")
pu[[4]]=DimPlot(s1, reduction = "umap", group.by="orig.ident") + scale_color_brewer(palette="Paired")
pu[[5]]=DimPlot(s1, reduction = "umap", group.by="Phase")

pdf(file=cc("seuratQC",args$PROJNAME,plotNo(),"UMAP",nDims,".pdf"),width=11,height=8.5)
print(pu)
dev.off()

clusterRes="integrated_snn_res.0.5"
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
        arrange(desc(avg_log2FC))
}

#

halt("BUG")

cl=filterCLTable(clusterMarkers) %>%
    dplyr::select(cluster,gene,p_val_adj,avg_log2FC,pct.1,pct.2,lOR) %>%
    filter(p_val_adj<FDR.cut & avg_log2FC>logFC.cut)

clusters=cl %>% distinct(cluster) %>% arrange(cluster) %>% pull %>% levels


clusterMarkerTbl=list()
for(ci in clusters) {
    clusterMarkerTbl[[ci]]=cl %>% filter(cluster==ci)
}

geneCounts=cl %>% count(cluster)

ll=c(list(GeneCounts=geneCounts,AllCluster=cl),clusterMarkerTbl)

xfile=cc("tblClusterMarkers","",clusterRes,"FDR",FDR.cut,"logFC",logFC.cut)

write.xlsx(ll,paste0(xfile,".xlsx"),overwrite=T)

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
mergePNGs(cmFile)

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
    arrange(desc(avg_log2FC)) %>%
    distinct(gene,.keep_all=T) %>%
    head(12) %>%
    pull(gene)

pc1=DotPlot(s1,features=dot.gene.lOR)
pc2=DotPlot(s1,features=dot.gene.lFC)

pdf(file=cc("seuratQC",args$PROJNAME,plotNo(),"ClusterMarkersDot",".pdf"),width=11,height=8.5)
print(pc1)
print(pc2)
dev.off()
