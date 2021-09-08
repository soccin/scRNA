#'
#' Phase-II.b
#'
#' Do PCA, UMAP, Module Scores
#'

suppressPackageStartupMessages(require(stringr))

usage="
usage: doSeuratV5_02b.R [MODULE_FILE=file] PARAMS.yaml

    PARAMS.yaml     parameter file from pass2[integration]

    OPTIONAL:
        MODULE_FILE     TSV file with list of genes for modules

"

cArgs=commandArgs(trailing=T)

#
# Separate out any options arguments
#
optionals=grep("=",cArgs,value=T)

oArgs=list(MODULE_FILE=NULL)
if(len(optionals)>0) {
    require(stringr, quietly = T, warn.conflicts=F)
    parseArgs=str_match(optionals,"(.*)=(.*)")
    aa=apply(parseArgs,1,function(x){oArgs[[str_trim(x[2])]]<<-str_trim(x[3])})
}

cArgs=grep("=",cArgs,value=T,invert=T)

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
    file.arg=grep("--file=",commandArgs(),value=T)
    if(len(file.arg)>0) {
        SDIR=dirname(gsub(".*=","",file.arg))
    } else {
        SDIR="."
    }
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

cat("\nLoading pass2 RDA file ...")
d10X.integrate=readRDS(args$PASS2.RDAFile)
cat(" done\n\n")
cat("digest=",digest::digest(d10X.integrate),"\n")

##########################################################################
# PCA

# https://satijalab.org/seurat/archive/v3.0/s13k_tutorial.html
# Perform linear dimensional reduction

cat("\nComputing PCA ...")
s1=d10X.integrate
s1=RunPCA(s1,features=VariableFeatures(s1),approx=FALSE)
cat(" done\n\n")

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

cat("\nClustering ...")
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
pal2=c(brewer.paired(20))
pu=list()
pu[[1]]=DimPlot(s1, reduction = "umap", label=T, group.by="integrated_snn_res.0.1", label.size=6) + scale_color_manual(values=pal1) + ggtitle("integrated_snn_res.0.1")
pu[[2]]=DimPlot(s1, reduction = "umap", label=T, group.by="integrated_snn_res.0.2", label.size=6) + scale_color_manual(values=pal1) + ggtitle("integrated_snn_res.0.2")
pu[[3]]=DimPlot(s1, reduction = "umap", label=T, group.by="integrated_snn_res.0.5", label.size=6) + scale_color_manual(values=pal1) + ggtitle("integrated_snn_res.0.5")
pu[[4]]=DimPlot(s1, reduction = "umap", group.by="orig.ident") + scale_color_manual(values=pal2)
pu[[5]]=DimPlot(s1, reduction = "umap", group.by="Phase")

pdf(file=cc("seuratQC",args$PROJNAME,plotNo(),"UMAP",nDims,".pdf"),width=11,height=8.5)
print(pu)
dev.off()
cat(" done\n\n")

if(!is.null(oArgs$MODULE_FILE)) {
    oArgs$MODULE_FILE=normalizePath(oArgs$MODULE_FILE)
    DefaultAssay(s1)="SCT"
    moduleTbl=read_tsv(oArgs$MODULE_FILE)
    colnames(moduleTbl)=c("Module","Gene")
    modules=split(moduleTbl$Gene,moduleTbl$Module)
    s1=AddModuleScore(s1,features=modules,name="Modules")

    cat("\nPlot modules ...")
    pm=list()
    for(ii in seq(len(modules))) {
        print(ii)
        pp=FeaturePlot(s1,
            features=paste0("Modules",ii),
            max.cutoff="q95",min.cutoff="q05",
            combine=F)

        pm[[ii]]=pp[[1]] + ggtitle(names(modules)[ii])
    }
    cat(" done\n\n")
    pfile=cc("seuratQC",args$PROJNAME,plotNo(),"ModuleScores_%03d.png")
    pngCairo(pfile,width=11,height=8.5)
    print(paginatePlots(pm,2,2,FALSE))
    dev.off()
    mergePNGs(pfile)

}

#
# Collect any changes in globals and parameters
#

DefaultAssay(s1)="integrated"

args$algoParams=ap
args$glbs=glbs
args$optionals=oArgs
args$GIT.Describe=git.describe(SDIR)
args.digest.orig=digest::digest(args)

args$PASS2b.RDAFile=cc("pass_02b",args.digest.orig,"s1",".rda")
write_yaml(args,cc("pass_02b","PARAMS.yaml"))

cat("\nSaving rda object ...")
saveRDS(s1,args$PASS2b.RDAFile,compress=T)
cat(" done\n\n")

