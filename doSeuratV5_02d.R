#'
#' Phase-II.c
#'
#' Do Find Cluster Markers
#'
STAGE=3

suppressPackageStartupMessages(require(stringr))

usage="
usage: doSeuratV5_02d.R PARAMS.yaml

    PARAMS.yaml     parameter file from pass2b
    ClusterRes      Resolution value of clusters to use
"

cArgs=commandArgs(trailing=T)

if(len(cArgs)<1) {
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
    library(openxlsx)
    library(pals)
})

source(file.path(SDIR,"seuratTools.R"))
source(file.path(SDIR,"plotTools.R"))

glbs=args$glbs
ap=args$algoParams

plotNo<-makeAutoIncrementor(30)

##########################################################################
#
# INCLUDE BREAK
#
##########################################################################


s1=readRDS(args$PASS2b.RDAFile)

if(len(cArgs)>1) {
    cResValue=cArgs[2]
    clusterRes=grep(cResValue,grep("res\\.",colnames(s1@meta.data),value=T),value=T)
} else {
    cat("\n\nNeed to specify a cluster resolution that will match with grep\n")
    cat("possible resolutions =>\n")
    cat("  ",paste(grep("res\\.",colnames(s1@meta.data),value=T),collapse=", "),"\n\n")
    quit()
}

if(!clusterRes %in% colnames(s1@meta.data)) {
    cat(paste("\n\n", "Invalid Cluster Resolution",clusterRes, "\n"))
    cat("Valid values\n\t",paste(grep("integrated_snn_res.",colnames(s1@meta.data),value=T),collapse="; "),"\n\n")
    rlang::abort("\n\nInvalid Cluster Resolution\n\n")
}

s1=SetIdent(s1,value=clusterRes)

clusterMarkers=FindAllMarkers(s1,only.pos=TRUE,logfc.threshold=0.25,min.pct = 0.25)

if(nrow(clusterMarkers)==0) {
    cat("\n\tNo cluster Markers Found\n\n")
    quit()
}

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

clustTag=gsub("inte.*res.","cRes_",clusterRes)

xfile=cc("tblClusterMarkers",clustTag,"FDR",FDR.cut,"logFC",logFC.cut)

write.xlsx(ll,paste0(xfile,".xlsx"),overwrite=T)

pal1=c(cols25(25),brewer.dark2(8))

plt.cmark=cl %>% group_split(cluster) %>% map(plotClusterMarkers,s1,pal1)

cmFile=get_plot_filename(plotNo(),"ClusterMarkers",clustTag,"FDR",FDR.cut,"logFC",logFC.cut,"%03d",".png")

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

pc1=DotPlot(s1,features=dot.gene.lOR) + scale_x_discrete(guide = guide_axis(n.dodge = 2))
pc2=DotPlot(s1,features=dot.gene.lFC) + scale_x_discrete(guide = guide_axis(n.dodge = 2))

pdf(file=get_plot_filename(plotNo(),"ClusterMarkersDot",clustTag,"FDR",FDR.cut,"logFC",logFC.cut,".pdf"),width=11,height=8.5)
print(pc1)
print(pc2)
dev.off()

nClust=nlevels(cl$cluster)
nGenes=floor(75/nClust)
genesHeat=cl %>% group_by(cluster) %>% top_n(n=nGenes,wt=avg_log2FC) %>% pull(gene) %>% unique
nCells=min(ncol(s1),5000)
sh=s1[,sample(colnames(s1),nCells)]
ph=DoHeatmap(sh,features=genesHeat)

pdf(file=get_plot_filename(plotNo(),"ClusterHeatmap",clustTag,"FDR",FDR.cut,"logFC",logFC.cut,".pdf"),width=11,height=8.5)
print(ph)
dev.off()

umapGenes=cl %>%
    filter(pct.1>.75 & pct.2<.25) %>%
    group_by(cluster) %>%
    top_n(n=3,wt=lOR) %>%
    ungroup %>%
    arrange(cluster,desc(lOR)) %>%
    distinct(gene) %>%
    pull

gumap=FeaturePlot(s1,features=umapGenes,combine=F,max.cutoff="q95")
pgu=paginatePlots(gumap,2,3,oneLegend=F)

umFile=get_plot_filename(plotNo(),"ClusterUMAP",clustTag,"FDR",FDR.cut,"logFC",logFC.cut,"%03d",".png")

png(filename=umFile,
    type="cairo",
    units="in",
    width=14,
    height=8.5,
    pointsize=12,
    res=96)

print(pgu)

dev.off()
mergePNGs(umFile)
