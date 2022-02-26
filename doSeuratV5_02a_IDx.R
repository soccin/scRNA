#'
#' Phase-II
#'
#' Do SCTransform but simple MERGE instead of Integrate
#' Also Swith to not regress Cell Cycle Out
#'

suppressPackageStartupMessages(require(stringr))

usage="
usage: doSeuratV5_02.R PARAMS.yaml

    PARAMS.yaml     parameter file from pass1
    CC_REGRESS      Flag to control cell cycle regression [Def: FALSE]

"
cArgs=commandArgs(trailing=T)
args=list(CC_REGRESS=FALSE)
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

args=c(args,args1)

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
source(file.path(SDIR,"doQCandFilter.R"))

glbs=args$glbs
ap=args$algoParams

plotNo<-makeAutoIncrementor(10)

##########################################################################
#
# INCLUDE BREAK
#
##########################################################################


d10X=readRDS(args$PASS1.RDAFile)

cat("digest=",digest::digest(d10X),"\n")

cat("\nDoQCandFilter\n")
for(ii in names(d10X)) {
    print(ii)
    if(grepl("\\.ID\\d",ii)) {
        cat("Skipping")
    } else {
        ret=doQCandFilter(d10X[[ii]], ap$MIN_NCOUNT_RNA, ap$MIN_FEATURE_RNA, ap$PCT_MITO)
        d10X[[ii]]=ret$so
    }
    cat("\n")
}

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
cat("\nMERGE Samples - IDx join\n")

merge=d10X[[which(!grepl("\\.ID\\d",names(d10X)))[1]]]
merge=AddMetaData(merge,gsub(".*bc:","",rownames(merge@meta.data)),"BC")

md=merge@meta.data
expr=list()

for(ij in grep("\\.ID\\d",names(d10X),value=T)) {

    gene=gsub(".*\\.ID","ID",ij)
    geneTag=paste0("pcr.",gene)
    cat(ij,gene,"\n")
    rna=d10X[[ij]]@assays$RNA
    g.exp=as.numeric(rna[gene,])
    g.exp=ifelse(is.na(g.exp),0,g.exp)
    xx=tibble(BC=gsub(".*bc:","",colnames(rna)),Gene=g.exp)
    colnames(xx)[2]=geneTag
    expr[[geneTag]]=xx
    md=left_join(md,xx)
    md[[geneTag]]=ifelse(is.na(md[[geneTag]]),0,md[[geneTag]])
    merge=AddMetaData(merge,md[[geneTag]],geneTag)

}

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

pdf(file=cc("seuratQC",args$PROJNAME,plotNo(),"PostIntegrateCC.pdf"),width=11,height=8.5)
plotCellCycle(d10X.integrate,"Post Integration CC Regression")
DimPlot(d10X.integrate,reduction="umap",group.by="Phase")
dev.off()

Idents(d10X.integrate)<-"SampleID"

so=FindVariableFeatures(d10X.integrate)
pv1=VariableFeaturePlot(so)
top10 <- head(VariableFeatures(so), 10)
pv2 <- LabelPoints(plot = pv1, points = top10, repel = TRUE, xnudge=0, ynudge=0)
pdf(file=cc("seuratQC",args$PROJNAME,plotNo(),"VariableFeatures.pdf"),width=11,height=8.5)
print(pv2)
dev.off()

args$glbs=glbs
args$algoParams=ap
args$GIT.Describe=git.describe(SDIR)
args.digest.orig=digest::digest(args)

args$PASS2.RDAFile=cc("pass_02a","SObj",args.digest.orig,"d10X.integrate",".rda")
write_yaml(args,cc("pass_02","PARAMS.yaml"))

saveRDS(d10X.integrate,args$PASS2.RDAFile,compress=T)


