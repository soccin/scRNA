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
    library(cowplot)
    library(patchwork)
    library(tidyverse)
    library(gridExtra)
})

source("seuratTools.R")

plotNo<-makeAutoIncrementor(10)

d10X.orig=readRDS(args$PASS1.RDAFile)

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

    pdf(file=cc("seuratQC",plotNo(),sampleId,"01.pdf"),height=8.5,width=11)
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
        d10X[[ii]]=subset(xx,cells=Cells(xx)[runif(nrow(xx@meta.data))<args$DOWNSAMPLE])
    }

    cat("digest=",digest::digest(d10X),"\n")

}

stop("DDDDD")

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

pdf(file=cc("seuratQC",plotNo(),"CellCycle.pdf"),width=11,height=11)

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



if(len(d10X)>1) {
    cat("\n\n The rest of this workflow only works on merged datasets\n\n")
    quit()
}


s0=d10X[[1]]

ret=regressCellCycle(d10X[[1]])

s1=ret$so

save.image(cc("CHECKPOINT",DATE(),glb.digest,".Rdata"),compress=T)


stop("Continue working")


pdf(file=cc("seuratQC",plotNo(),"PostCCRegress.pdf"),width=11,height=8.5)

plotCellCycle(s1,"Post CC Regression")

p3=DimPlot(s1, reduction = "umap", split.by = "orig.ident")
print(p3)

dev.off()

