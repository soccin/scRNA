suppressPackageStartupMessages(require(stringr))
STAGE=1

usage="
usage:
        doSeuratV5_01.R [CONFIG=${CONFIG}] [DEBUG=${DEBUG}] [MERGE=${MERGE}] [PROJNAME=${PROJNAME}] 10XDir1 [10XDir2 ... ]

    DEBUG        Set DEBUG mode (downsample to 10%) [${DEBUG}]
    DOWNSAMPLE   Set amount of DEBUG downsample [${DOWNSAMPLE}]
    MERGE        If True then merge samples with simple merge [${MERGE}]
    PROJNAME     Set name of project. Must either be used or there
                 must be a file PROJNAME in folder with name
    CONFIG       YAML file used to specify explicit paths to cellRanger data folder and genome
                 if used then directories not read from command line

  Pass 1, check QC and filtering levels and also generate Cell Cycle
  Regression plots to see if C.C. needs to be regressed out.

"

cArgs=commandArgs(trailing=T)
args=list(DEBUG=FALSE,MERGE=TRUE,PROJNAME="scRNA",DOWNSAMPLE=0.1,CONFIG=".none")
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

if(args$CONFIG!=".none") {
    cat("\nUsing config file for data paths\n\n")
    library(yaml)
    require(dplyr)

    if(!file.exists(args$CONFIG)) {
        cat("\n\n\tConfig file",args$CONFIG,"does not exist\n\n")
        stop("FATAL:ERROR")
    }

    config=read_yaml(args$CONFIG)
    inputs=purrr::map(config$inputs,as_tibble) %>% bind_rows
    dataFolders=inputs$dir
    sampleIDs=inputs$sid
    names(dataFolders)=sampleIDs

} else {

    #halt("FIX SAMPLE IDS")
    #
    dataFolders=argv
    sampleIDs=gsub("_",".",gsub(".outs.*","",gsub("^s_","",basename(dataFolders))))
    names(dataFolders)=sampleIDs
    config=NULL

}

if(len(dataFolders)<1) {
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



#
# QC Defaults
#
MIN_FEATURE_RNA=1500
MIN_NCOUNT_RNA=5000
PCT_MITO=10
METHOD="_auto"
XENO_PCTCUT=80

getDefault<-function(ll,key) ifelse(is.null(ll[[key]]),get(key),ll[[key]])

if(exists("args00")) {
    MIN_FEATURE_RNA=getDefault(args00$algoParams,"MIN_FEATURE_RNA")
    MIN_NCOUNT_RNA=getDefault(args00$algoParams,"MIN_NCOUNT_RNA")
    PCT_MITO=getDefault(args00$algoParams,"PCT_MITO")
    METHOD=getDefault(args00$algoParams,"METHOD")
    XENO_PCTCUT=getDefault(args00$algoParams,"XENO_PCTCUT")
}

algoParams=list()
algoParams$MIN_FEATURE_RNA=MIN_FEATURE_RNA
algoParams$MIN_NCOUNT_RNA=MIN_NCOUNT_RNA
algoParams$PCT_MITO=PCT_MITO
algoParams$METHOD=METHOD
algoParams$XENO_PCTCUT=XENO_PCTCUT
algoParams$SEED=101


##############################################################################
# Read 10X data
#

set.seed(algoParams$SEED)

d10X=list()
for(ii in seq(len(dataFolders))) {
    sampleName=sampleIDs[ii]
    cat("Reading Sample =",sampleName,"...")

    if(!is.null(config$multi)) {

        #
        # Need to explicity set genome also from config
        #
        glbs$genome <- genomes[config$genome]

        #read_10X_multi_as_SeurateObject(bc_matrix_dir,gexSlot,genome,projName)

        d10X[[sampleName]] <- read_10X_multi_as_SeurateObject(
                                        dataFolders[ii],
                                        sampleName,
                                        config$multi$gex,
                                        glbs$genome,
                                        args$PROJNAME
                                    )

    } else {

        # Orig default path

        #
        # If we are using a config file (in non-multi mode)
        # check if genome was set and set the global
        #

        if(!is.null(config$genome)) {
           glbs$genome <- genomes[config$genome]
        }

        d10X[[sampleName]] <- read10XDataFolderAsSeuratObj(dataFolders[ii],args$PROJNAME)

    }

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
    pdf(file=get_plot_filename(plotNo(),"QC.pdf"),width=11,height=8.5)
    print(ret$plts)
    dev.off
}

args.digest.orig=digest::digest(args)
args$PASS1.RDAFile=cc("pass_01","SObj",args.digest.orig,"d10X.orig",".rda")
write_yaml(args,cc("pass_01","PARAMS.yaml"))

saveRDS(d10X.orig,args$PASS1.RDAFile,compress=T)

##############################################################################
##############################################################################
# Dump MAD3 levels for possible pass_00_PARAMS
##############################################################################
##############################################################################

cutOffsMAD3=ret$stats[[1]]$Cutoff
names(cutOffsMAD3)=ret$stats[[1]]$Feature
cutOffsMAD3=as.list(cutOffsMAD3)

args00Mad3=list(
    algoParams=list(
        MIN_FEATURE_RNA=cutOffsMAD3$nFeature_RNA,
        MIN_NCOUNT_RNA=cutOffsMAD3$nCount_RNA,
        PCT_MITO=cutOffsMAD3$percent.mt,
        METHOD="MAD3"
    )
)
if(glbs$genome=="xenograft") { # Set xeno percent cut for cell filter later
    args00Mad3$algoParams$XENO_PCTCUT=80
}

write_yaml(args00Mad3,"pass_00_PARAMS.yaml.mad3")

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

pdf(file=get_plot_filename(plotNo(),"CellCycle.pdf"),width=11,height=8.5)

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


if(glbs$genome=="xenograft") { # Plot xenoqc graphs

md=d10X[[1]]@meta.data %>%
    tibble %>%
    mutate(Species=case_when(percent.Hs>95 ~ "Human", percent.Mm>95 ~ "Mouse", T ~ "Mixed")) %>%
    mutate(Species=factor(Species))

px1=md %>%
    count(SampleID,Species,.drop=F) %>%
    ggplot(aes(SampleID,n,fill=Species)) +
        theme_light(16) +
        geom_col(position='dodge') +
        scale_fill_brewer(palette="Dark2",drop=F) +
        coord_flip() +
        theme(panel.grid.minor=element_blank())

px2=px1 + scale_y_continuous(trans="log1p",breaks=c(0,1,10,100,1000,10000))

px3=md %>% gather(Species2,PCT,percent.Hs,percent.Mm) %>%
    ggplot(aes(nCount_RNA,PCT,color=Species)) +
    geom_point(alpha=.4) +
    scale_x_log10() +
    theme_light(12) +
    facet_wrap(~SampleID)

xPos=10^floor(log10(max(md$nCount_RNA)))
adf=md %>% group_by(SampleID) %>%
    summarize(PCT.Human=sprintf("Percent Human = %.1f%%",mean(Species=="Human")*100)) %>%
    mutate(x=xPos,y=50)
px3=px3 + geom_label(data=adf,aes(x=x,y=y,label=PCT.Human),color="black",hjust=1,size=4)


pdf(file=get_plot_filename(plotNo(),"XenoStats.pdf"),width=11,height=8.5)
print(px1)
print(px2)
print(px3)
dev.off()
}
#d10X[[1]]@meta.data %>% tibble %>% ggplot(aes(percent.Hs,percent.Mm)) + geom_point()

