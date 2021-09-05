#'
#' Phase-II.b
#'
#' Do ???
#'

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


d10X.????=readRDS(args$????)

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
