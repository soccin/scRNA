suppressPackageStartupMessages(require(stringr))

usage="
usage:
        doSeuratV5_01.R [DEBUG=${DEBUG}] [MERGE=${MERGE}] [PROJNAME=${PROJNAME}] 10XDir1 [10XDir2 ... ]

    DEBUG        Set DEBUG mode (downsample to 10%) [${DEBUG}]
    DOWNSAMPLE   Set amount of DEBUG downsample [${DOWNSAMPLE}]
    MERGE        If True then merge samples with simple merge [${MERGE}]
    PROJNAME     Set name of project. Must either be used or there
                 must be a file PROJNAME in folder with name

  Pass 1, check QC and filtering levels and also generate Cell Cycle
  Regression plots to see if C.C. needs to be regressed out.

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

if(Sys.getenv("SDIR")=="") {
    #
    # getSDIR defined in .Rprofile
    #
    SDIR=getSDIR()
} else {
    SDIR=Sys.getenv("SDIR")
}

##############################################################################
if(len(argv)<1 && !file.exists("pass_00_PARAMS.yaml")) {
    cat(usage)
    quit()
}

if(R.Version()$major<4) {
    cat(usage)
    cat("\n\nThis script needs version(R).major>=4\n\n")
    quit()
}


cat("\n=========================================================\n")
args[["10XDirs"]]=paste0(argv,collapse=", ")
cat(str(args))
cat("\n")
##############################################################################

suppressPackageStartupMessages({
    library(Seurat)
    library(patchwork)
    library(tidyverse)
    library(yaml)
})

source(file.path(SDIR,"seuratTools.R"))
source(file.path(SDIR,"plotTools.R"))
source(file.path(SDIR,"qcAndFilter.R"))

##############################################################################
# Set up global variables, parameters and defaults
#

plotNo<-makeAutoIncrementor()

if(file.exists("pass_00_PARAMS.yaml")) {
    cat("\n   Loading defaults from pass_00_PARAMS.yaml\n\n\n")
    args00=read_yaml("pass_00_PARAMS.yaml")
} else {
    args00=list()
    cat("\n   No pass_00_PARAMS file; using default PARAMS\n\n\n")
}

#halt("FIX SAMPLE IDS")

dataFolders=argv
sampleIDs=gsub("_",".",gsub(".outs.*","",gsub("^s_","",basename(dataFolders))))
names(dataFolders)=sampleIDs


#
# QC Defaults
#
MIN_FEATURE_RNA=1500
MIN_NCOUNT_RNA=5000
PCT_MITO=10

getDefault<-function(ll,key) ifelse(is.null(ll[[key]]),get(key),ll[[key]])

if(exists("args00")) {
    MIN_FEATURE_RNA=getDefault(args00$algoParams,"MIN_FEATURE_RNA")
    MIN_NCOUNT_RNA=getDefault(args00$algoParams,"MIN_NCOUNT_RNA")
    PCT_MITO=getDefault(args00$algoParams,"PCT_MITO")
}

algoParams=list()
algoParams$MIN_FEATURE_RNA=MIN_FEATURE_RNA
algoParams$MIN_NCOUNT_RNA=MIN_NCOUNT_RNA
algoParams$PCT_MITO=PCT_MITO
algoParams$SEED=101


##############################################################################
# Read 10X data
#

set.seed(algoParams$SEED)

d10X=list()
for(ii in seq(len(dataFolders))) {
    sampleName=sampleIDs[ii]
    cat("Reading Sample =",sampleName,"...")
    d10X[[sampleName]] <- read10XDataFolderAsSeuratObj(dataFolders[ii],args$PROJNAME)
    cat("\n")
}


d10X.orig=d10X

glb.digest=digest::digest(d10X.orig)
cat("digest=",digest::digest(d10X.orig),"\n")

##############################################################################
# Downsample if DEBUG set
#

if(args$DEBUG) {

    c("\n\nDEBUG SET; Subset data\n\n")

    cat("\nDEBUG::subset\n")

    for(ii in seq(d10X)) {
        print(ii)
        xx=d10X[[ii]]
        d10X[[ii]]=subset(xx,cells=Cells(xx)[runif(nrow(xx@meta.data))<args$DOWNSAMPLE])
    }

    cat("digest=",digest::digest(d10X),"\n")

}

##############################################################################
# Merge samples if MERGE set
#
cat("\nMERGE Samples\n")

if(args$MERGE & len(d10X)>1) {
    cat("\nMerging sample files...")
    s.merge=merge(d10X[[1]],d10X[-1],project=args$PROJNAME)
    d10X=list()
    d10X[[args$PROJNAME]]=s.merge
    cat("done\n\n")
}

#
# Add SampleID metadata, if there is a manifest
# use that for the id's otherwise make them orig.ident
#

md=d10X[[1]]@meta.data
if(is.null(args00$SAMPLE_MANIFEST)) {
    md$SampleID=md$orig.ident
} else {
    args$SAMPLE_MANIFEST=args00$SAMPLE_MANIFEST
    manifest=read_csv(args00$SAMPLE_MANIFEST)
    md=md %>% rownames_to_column("CELLID") %>% left_join(manifest,by="orig.ident") %>% column_to_rownames("CELLID")
}
d10X[[1]]@meta.data=md


##############################################################################
# Save load cache file for faster reloading by filter test script
#   doSeuratV5_01_Filter.R
#

args$glbs=glbs
args$algoParams=algoParams

args$GIT.Describe=git.describe(SDIR)
args.digest.orig=digest::digest(args)

cobj=list(args=args,cArgs=cArgs,glbs=glbs,d10X.orig=d10X.orig,d10X=d10X)
args$PASS0.CACHE=cc("preFilter","CACHE",substr(digest::digest(cArgs),1,7),".rda")
write_yaml(args,cc("preFilter","PARAMS.yaml"))
saveRDS(cobj,args$PASS0.CACHE,compress=T)

##############################################################################
# Do Stage-I QC
#
cat("\nDoQCandFilter\n")

Idents(d10X[[1]])<-"SampleID"

cat("md5(dX10) =",digest::digest(d10X),"\n")

for(ii in seq(d10X)) {
    print(ii)
    ret=qcSamples(d10X[[ii]])
    pdf(file=cc("seuratQC",args$PROJNAME,plotNo(),"QC.pdf"),width=11,height=8.5)
    print(ret$plts)
    dev.off
}

args.digest.orig=digest::digest(args)
args$PASS1.RDAFile=cc("pass_01","SObj",args.digest.orig,"d10X.orig",".rda")
write_yaml(args,cc("pass_01","PARAMS.yaml"))

saveRDS(d10X.orig,args$PASS1.RDAFile,compress=T)

##############################################################################
##############################################################################
# Move cell cycle here and run if filtering is set by pass_00_PARAMS.yaml
#     if(file.exists("pass_00_PARAMS.yaml"))
##############################################################################
##############################################################################

##############################################################################
# Check Cell Cycle
#
cat("\nScoreCellCycle\n")

for(ii in seq(d10X)) {
    print(ii)
    d10X[[ii]]=scoreCellCycle(d10X[[ii]])
}

pcc=list()
cat("\nPlotCellCycle\n")

for(ii in seq(d10X)) {
    pcc[[ii]]=plotCellCycle(preProcessSO(d10X[[ii]]),names(d10X)[ii])
}

pdf(file=cc("seuratQC",args$PROJNAME,plotNo(),"CellCycle.pdf"),width=11,height=8.5)

# if(len(pcc)>1) {
#     nPages=ceiling(len(pcc)/4)
#     for(ii in seq(nPages)) {
#         jj=(1:4)+4*(ii-1)
#         jj=intersect(seq(pcc),jj)
#         pOut=pcc[jj]
#         if(len(pOut)<4){
#             pOut=c(pOut,rep(list(ggplot()+theme_void()),4-len(jj)))
#         }
#         print(wrap_plots(pOut,ncol=2))
#     }
# } else {
   print(pcc)
#}

dev.off()

