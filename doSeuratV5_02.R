suppressPackageStartupMessages(require(stringr))

usage="
usage: doSeuratV5.R [DEBUG=${DEBUG}] [MERGE=${MERGE}] [PROJNAME=${PROJNAME}] 10XDir1 [10XDir2 ... ]

    DEBUG        Set DEBUG mode (downsample to 10%) [${DEBUG}]
    DOWNSAMPLE   Set amount of DEBUG downsample [${DOWNSAMPLE}]
    MERGE        If True then merge samples with simple merge [${MERGE}]
    PROJNAME     Set name of project. Must either be used or there
                 must be a file PROJNAME in folder with name

"

cArgs=commandArgs(trailing=T)
args=list(DEBUG=FALSE,MERGE=TRUE,PROJNAME="scRNA",DOWNSAMPLE=0.1)
usage=str_interp(usage,args)

ii=grep("=",cArgs)
if(len(ii)>0) {
    parseArgs=str_match(cArgs[ii],"(.*)=(.*)")
    aa=apply(parseArgs,1,function(x){args[[str_trim(x[2])]]<<-str_trim(x[3])})
}

args$DEBUG=as.logical(args$DEBUG)
args$MERGE=as.logical(args$MERGE)
args$DOWNSAMPLE=as.numeric(args$DOWNSAMPLE)

if(args$PROJNAME=="scRNA") {
    if(file.exists("PROJNAME")) {
        args$PROJNAME=scan("PROJNAME","",quiet=T)
    } else {
        cat(paste0(usage,"  You have not set a project name\n\n"))
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

plotNo<-makeAutoIncrementor()

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

