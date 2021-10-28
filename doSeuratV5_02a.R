#'
#' Phase-II
#'
#' Do SCTransform and Integration
#'

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
for(ii in seq(d10X)) {
    print(ii)
    ret=doQCandFilter(d10X[[ii]], ap$MIN_NCOUNT_RNA, ap$MIN_FEATURE_RNA, ap$PCT_MITO)
    d10X[[ii]]=ret$so
}


##############################################################################
# Check Cell Cycle
#
cat("\nScoreCellCycle\n")

for(ii in seq(d10X)) {
    print(ii)
    d10X[[ii]]=scoreCellCycle(d10X[[ii]])
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
# Do SCTransform
# - This is where the cell cycle gets regressed out
#
cat("\nSCTransform\n")

d10X.int=list()
for(ii in seq(d10X)) {
    print(ii)
    d10X.int[[ii]]=SCTransform(d10X[[ii]],vars.to.regress = c('S.Score', 'G2M.Score'))
}

##############################################################################
# Do Integration
#
cat("\nSCTransform\n")

features <- SelectIntegrationFeatures(object.list = d10X.int, nfeatures = 3000)
d10X.int <- PrepSCTIntegration(object.list = d10X.int, anchor.features = features)

anchors <- FindIntegrationAnchors(
        object.list = d10X.int,
        normalization.method = "SCT",
        anchor.features = features
        )

d10X.integrate <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

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

args$PASS2.RDAFile=cc("pass_02a",args.digest.orig,"d10X.integrate",".rda")
write_yaml(args,cc("pass_02","PARAMS.yaml"))

saveRDS(d10X.integrate,args$PASS2.RDAFile,compress=T)
