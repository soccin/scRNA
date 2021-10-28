suppressPackageStartupMessages(require(stringr))

usage="
usage: doSeuratV5_02.R pass_02b_PARAMS.yaml PCA_DIMS

    pass_02b_PARAMS.yaml        Parameter file from pass2 post PCA
    PCA_DIMS                    Number of PCA Dimensions to use
"

cArgs=commandArgs(trailing=T)

if(len(cArgs)!=2) {
    cat(usage)
    quit()
}

if(R.Version()$major<4) {
    cat(usage)
    cat("\n\nThis script needs version(R).major>=4\n\n")
    quit()
}

if(Sys.getenv("SDIR")=="") {
    #
    # getSDIR defined in .Rprofile
    #
    SDIR=getSDIR()
} else {
    SDIR=Sys.getenv("SDIR")
}

source(file.path(SDIR,"seuratTools.R"))
source(file.path(SDIR,"plotTools.R"))

library(yaml)
args=read_yaml(cArgs[1])
nDims=as.numeric(cArgs[2])

glbs=args$glbs
ap=args$algoParams
ap$NDIMS=nDims
ap$ClusterResolutions=c(0.1,0.2,0.5,0.8)

plotNo<-makeAutoIncrementor(20)

suppressPackageStartupMessages({
    library(Seurat)
    library(patchwork)
    library(tidyverse)
    library(openxlsx)
    library(pals)
})



##########################################################################
#
# INCLUDE BREAK
#
##########################################################################

s1=readRDS(args$PASS2b.RDAFile)


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
    cTbl=md %>% count(SampleID,.data[[clusterI]]) %>% rename(Clusters=all_of(clusterI)) %>% mutate(Clusters=factor(Clusters,levels=cLevels))
    pc[[len(pc)+1]]=ggplot(cTbl,aes(y=Clusters,x=n,fill=SampleID)) + geom_bar(position="fill", stat="identity") + scale_fill_manual(values=cols25()) + theme_light(base_size=18) + ggtitle(clusterI)
    pc[[len(pc)+1]]=ggplot(cTbl,aes(fill=Clusters,x=n,y=SampleID)) + geom_bar(position="fill", stat="identity") + scale_fill_manual(values=pal1) + theme_light(base_size=18) + ggtitle(clusterI)
}

pdf(file=cc("seuratQC",args$PROJNAME,plotNo(),"ClusterChart",nDims,".pdf"),width=14,height=8.5)
print(pc)
dev.off()
