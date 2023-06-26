#'
#' Phase-II
#'
#' Do SCTransform but simple MERGE instead of Integrate
#' Also Switch to not regress Cell Cycle Out
#'
STAGE=1

suppressPackageStartupMessages(require(stringr))

usage="
usage: doSeuratV5_02a_MergeOnly.R [CC_REGRESS=TRUE] [CELL_FILTER=filterFile.csv] [GENE_FILTER=geneFile] PARAMS.yaml

    PARAMS.yaml     parameter file from pass1
    CC_REGRESS      Flag to control cell cycle regression [Def: TRUE]
    CELL_FILTER     File of cells to filter out
    GENE_FILTER     File of genes to filter out

"
cArgs=commandArgs(trailing=T)
args=list(CC_REGRESS=TRUE,CELL_FILTER=NULL,GENE_FILTER=NULL)
usage=str_interp(usage,args)

ii=grep("=",cArgs)
if(len(ii)>0) {
    parseArgs=str_match(cArgs[ii],"(.*)=(.*)")
    aa=apply(parseArgs,1,function(x){args[[str_trim(x[2])]]<<-str_trim(x[3])})
}

args$CC_REGRESS=as.logical(args$CC_REGRESS)

argv=grep("=",cArgs,value=T,invert=T)

if(len(argv)<1) {
    cat(usage)
    quit()
}

if(R.Version()$major<4) {
    cat(usage)
    cat("\n\nThis script needs version(R).major>=4\n\n")
    quit()
}

library(yaml)
args1=read_yaml(argv[1])

args=modifyList(args,args1)

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
source(file.path(SDIR,"qcAndFilter.R"))

glbs=args$glbs
ap=args$algoParams

plotNo<-makeAutoIncrementor(10)

##########################################################################
#
# INCLUDE BREAK
#
##########################################################################


d10X=readRDS(args$PASS1.RDAFile)

#
# Add gene filters also
# 

if(!is.null(args$GENE_FILTER)) {
    cat("\n   Running gene file with file ",args$GENE_FILTER,"\n\n")
    rna=d10X[[1]]@assays$RNA
    allGenes=rownames(rna@counts)
    genesToFilter=scan(args$GENE_FILTER,"")
    genesToKeep=setdiff(allGenes,genesToFilter)
}

cat("digest=",digest::digest(d10X),"\n")

##############################################################################
# Do QC and Filter Cells
#
cat("\nDoQCandFilter\n")

qTbls=list()
for(ii in seq(d10X)) {
    print(ii)

    so = d10X[[ii]]
    md=so@meta.data
    if(is.null(args$SAMPLE_MANIFEST)) {
        md$SampleID=md$orig.ident
    } else {
        manifest=read_csv(args$SAMPLE_MANIFEST)
        md=md %>% rownames_to_column("CELLID") %>% left_join(manifest,by="orig.ident") %>% column_to_rownames("CELLID")
    }
    so@meta.data=md

    qTbls[[len(qTbls)+1]]=get_qc_tables(so@meta.data,ap)

    so = apply_filter01(so,ap)

    if(!is.null(args$GENE_FILTER)) {
        so = so[genesToKeep,]
    }

    d10X[[ii]]=so

}

######################################################################
# Post filter QC report
#

names(qTbls)=names(d10X)

tblCutOff=map(qTbls,1) %>% bind_rows(.id="SampleID") %>% spread(SampleID,Cutoff)
tblFailN=map(qTbls,5) %>% bind_rows
tblFailPCT=map(qTbls,4) %>% bind_rows %>% mutate_if(is.numeric,\(x) sprintf("%.2f%%",round(100*x,2)))

tblTotals=map(qTbls,2) %>%
    bind_rows(.id="SampleID") %>%
    mutate(Pass=Count.RNA&Num.Features&PCT.MT) %>%
    group_by(SampleID,Pass) %>%
    summarize(N=sum(n)) %>%
    ungroup %>%
    mutate(Pass=ifelse(Pass,"Keep","Fail")) %>%
    spread(Pass,N) %>%
    mutate(N=Fail+Keep) %>%
    mutate(PCT.Keep=sprintf("%.2f%%",100*Keep/N)) %>%
    select(SampleID,N,Fail,Keep,PCT.Keep)

ptb=list()
ptb[[1]]=ggplot()+theme_void()+annotation_custom(tableGrob(tblCutOff,rows=NULL))
ptb[[2]]=ggplot()+theme_void()+annotation_custom(tableGrob(tblFailN,rows=NULL))
ptb[[3]]=ggplot()+theme_void()+annotation_custom(tableGrob(tblTotals,rows=NULL))
ptb[[4]]=ggplot()+theme_void()+annotation_custom(tableGrob(tblFailPCT,rows=NULL))

pdf(file=get_plot_filename(plotNo(),"PostFilterQCTbls.pdf"),width=11,height=8.5)
print(ptb[[1]]/(ptb[[2]]+ptb[[3]])/ptb[[4]])
dev.off()

######################################################################
######################################################################

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
# Filter cells if file provided
#
# The file is a list of cells to _REMOVE_
# It must contain valid cell barcodes in a column
# called $CellID
#

if(!is.null(args$CELL_FILTER)) {

    c("\n\nFiltering Cells...\n\n")
    cellFilter=read_csv(args$CELL_FILTER)
    for(ii in seq(d10X)) {
        print(ii)
        so=d10X[[ii]]
        d10X[[ii]]=subset(so,cells=setdiff(Cells(so),cellFilter$CellID))
    }
}

##############################################################################
# Check Cell Cycle
#
cat("\nScoreCellCycle\n")

for(ii in seq(d10X)) {
    print(ii)
    d10X[[ii]]=scoreCellCycle(d10X[[ii]])
}

##############################################################################
# Merge samples if MERGE set
#
cat("\nMERGE Samples\n")

cat("\nMerging sample files...")
merge=merge(d10X[[1]],d10X[-1],project=args$PROJNAME)
cat("done\n\n")

##############################################################################
# Do SCTransform
# - This is where the cell cycle gets regressed out
#
cat("\nSCTransform\n")

if(args$CC_REGRESS) {
    cat("\n\nRegress out Cell Cycle\n\n")
    d10X.integrate=SCTransform(merge,vars.to.regress = c('S.Score', 'G2M.Score'))
} else {
    cat("\n\nNO Cell Cycle regression done\n\n")
    d10X.integrate=SCTransform(merge)
}

##############################################################################
# PCA + UMAP on merge/integrated data
#

d10X.integrate <- RunPCA(d10X.integrate, verbose = FALSE)
d10X.integrate <- RunUMAP(d10X.integrate, reduction = "pca", dims = 1:30)

cellCycle.genes = getCellCycleGenes(glbs$genome)

so=CellCycleScoring(d10X.integrate,
                    s.features=cellCycle.genes$s.genes,
                    g2m.features=cellCycle.genes$g2m.genes,
                    set.ident=T
                    )

cc.meta.data=so@meta.data[,c("S.Score","G2M.Score","Phase")]

d10X.integrate=AddMetaData(d10X.integrate,cc.meta.data$Phase,"Phase")
d10X.integrate=AddMetaData(d10X.integrate,cc.meta.data$G2M.Score,"G2M.Score")
d10X.integrate=AddMetaData(d10X.integrate,cc.meta.data$S.Score,"S.Score")
d10X.integrate=AddMetaData(d10X.integrate,cc.meta.data$S.Score-cc.meta.data$G2M.Score,"CC.Difference")

#
# Add SampleID metadata, if there is a manifest
# use that for the id's otherwise make them orig.ident
#
md=d10X.integrate@meta.data
if(is.null(args$SAMPLE_MANIFEST)) {
    md$SampleID=md$orig.ident
} else {
    cat("\n\n  Fixing sample names\n\n")
    manifest=read_csv(args$SAMPLE_MANIFEST)
    md=md %>% rownames_to_column("CELLID") %>% left_join(manifest,by="orig.ident") %>% column_to_rownames("CELLID")
}
d10X.integrate@meta.data=md

Idents(d10X.integrate)="Phase"

pdf(file=get_plot_filename(plotNo(),"PostIntegrateCC.pdf"),width=11,height=8.5)
plotCellCycle(d10X.integrate,"Post Integration CC Regression")
DimPlot(d10X.integrate,reduction="umap",group.by="Phase")
dev.off()

Idents(d10X.integrate)<-"SampleID"

so=FindVariableFeatures(d10X.integrate)
pv1=VariableFeaturePlot(so)
top10 <- head(VariableFeatures(so), 10)
pv2 <- LabelPoints(plot = pv1, points = top10, repel = TRUE, xnudge=0, ynudge=0)
pdf(file=get_plot_filename(plotNo(),"VariableFeatures.pdf"),width=11,height=8.5)
print(pv2)
dev.off()

args$glbs=glbs
args$CombineMethod="merge"
args$CombineMethodARGS=list()
args$algoParams=ap
args$GIT.Describe=git.describe(SDIR)
args.digest.orig=digest::digest(args)

args$PASS2.RDAFile=cc("pass_02a","SObj",args.digest.orig,"d10X.merge",".rda")
write_yaml(args,cc("pass_02","PARAMS.yaml"))

saveRDS(d10X.integrate,args$PASS2.RDAFile,compress=T)
