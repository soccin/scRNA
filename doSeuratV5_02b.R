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
    #
    # getSDIR defined in .Rprofile
    #
    SDIR=getSDIR()
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
for(ci in grep("integrated_snn_",colnames(s1@meta.data),value=T)) {

    clusterLevels=s1@meta.data[[ci]]

    s1@meta.data[[ci]]=factor(as.numeric(as.character(clusterLevels))+1,levels=sort(as.numeric(levels(clusterLevels)))+1)

    pu[[len(pu)+1]] <- DimPlot(s1, reduction = "umap", label=T, group.by=ci, label.size=6) + scale_color_manual(values=pal1) + ggtitle(ci)

}

pu[[len(pu)+1]] <- DimPlot(s1, reduction = "umap", group.by="SampleID") + scale_color_manual(values=cols25())
pu[[len(pu)+1]] <- DimPlot(s1, reduction = "umap", group.by="Phase")

pdf(file=cc("seuratQC",args$PROJNAME,plotNo(),"UMAP",nDims,".pdf"),width=11,height=8.5)
print(pu)
dev.off()
cat(" done\n\n")

md=s1@meta.data %>% rownames_to_column("CellID") %>% tibble
pc=list()
for(clusterI in grep("integrated_snn_res",colnames(md),value=T)) {

    cLevels=sort(as.numeric(levels(md[[clusterI]])))

    nClusters=len(cLevels)
    nSamples=md %>% distinct(SampleID) %>% nrow

    cTbl=md %>%
        count(SampleID,.data[[clusterI]]) %>%
        rename(Clusters=all_of(clusterI)) %>%
        mutate(Clusters=factor(Clusters,levels=cLevels)) %>%
        group_by(SampleID) %>%
        mutate(Total=sum(n)) %>%
        ungroup %>%
        mutate(nSampleNorm=n/Total) %>%
        arrange(Clusters,SampleID)

    S_sample=cTbl %>%
        select(SampleID,Clusters,n) %>%
        group_by(SampleID) %>%
        mutate(Total=sum(n)) %>%
        mutate(PCT=n/Total) %>%
        arrange(SampleID) %>%
        summarize(S=-sum(PCT*log(PCT))/log(nClusters)) %>%
        mutate(Clusters=0)

    S_cluster=cTbl %>%
        select(SampleID,Clusters,n) %>%
        group_by(SampleID) %>%
        mutate(Norm=sum(n)) %>%
        mutate(n=n/Norm) %>%
        group_by(Clusters) %>%
        mutate(Total=sum(n)) %>%
        mutate(PCT=n/Total) %>%
        arrange(Clusters) %>%
        summarize(S=-sum(PCT*log(PCT))/log(nSamples)) %>%
        mutate(SampleID="")

    pc[[len(pc)+1]]=ggplot(cTbl,aes(y=Clusters,x=nSampleNorm,fill=SampleID)) +
        geom_bar(position="fill", stat="identity") +
        scale_fill_manual(values=cols25()) +
        theme_light(base_size=18) +
        ggtitle(paste(clusterI,"Sample Cell Count Normalized")) +
        geom_label(data=S_cluster,aes(x=1,y=Clusters,label=sprintf("%.3f",S)),fill="white",hjust=1.1,size=4,label.size=.3)

    pc[[len(pc)+1]]=ggplot(cTbl,aes(fill=Clusters,x=n,y=SampleID)) +
        geom_bar(position="fill", stat="identity") +
        scale_fill_manual(values=pal1) +
        theme_light(base_size=18) +
        ggtitle(clusterI) +
        geom_label(data=S_sample,aes(x=1,y=SampleID,label=sprintf("%.3f",S)),fill="white",hjust=1.1,size=4,label.size=.3)

}

pdf(file=cc("seuratQC",args$PROJNAME,plotNo(),"ClusterChart",nDims,".pdf"),width=14,height=8.5)
print(pc)
dev.off()

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

args$PASS2b.RDAFile=cc("pass_02b","SObj",args.digest.orig,"s1",".rda")
write_yaml(args,cc("pass_02b","PARAMS.yaml"))

cat("\nSaving rda object ...")
saveRDS(s1,args$PASS2b.RDAFile,compress=T)
cat(" done\n\n")

