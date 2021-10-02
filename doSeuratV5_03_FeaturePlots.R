suppressPackageStartupMessages(require(stringr))

usage="
usage: doSeuratV5_02.R PARAMS.yaml GeneListFile

    PARAMS.yaml     parameter file from pass1
    GeneListFile    List of Genes to Plot
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
    file.arg=grep("--file=",commandArgs(),value=T)
    if(len(file.arg)>0) {
        SDIR=dirname(gsub(".*=","",file.arg))
    } else {
        SDIR="."
    }
} else {
    SDIR=Sys.getenv("SDIR")
}

source(file.path(SDIR,"seuratTools.R"))
source(file.path(SDIR,"plotTools.R"))

library(yaml)
args=read_yaml(cArgs[1])

glbs=args$glbs
ap=args$algoParams

plotNo<-makeAutoIncrementor(40)

suppressPackageStartupMessages({
    library(Seurat)
    library(patchwork)
    library(tidyverse)
    library(openxlsx)
    library(pals)
})

geneListFile=cArgs[2]
genes=scan(geneListFile,"")
if(args$glbs$genome=="mm10") {
    genes=str_to_title(genes)
}



##########################################################################
#
# INCLUDE BREAK
#
##########################################################################

s1=readRDS(args$PASS2b.RDAFile)

pp=FeaturePlot(s1,features=genes,combine=F,order=T,max.cutoff="q95",min.cutoff=0)

pgL=list()
for(ii in seq(genes)) {
    pg=NULL
    print(genes[ii])
    try({pg=FeaturePlot(s1,features=genes[ii],combine=F,order=T,max.cutoff="q95",min.cutoff=0,split.by="orig.ident")})
    if(!is.null(pg)) {
        pgL[[genes[ii]]]=pg
    }
}

pfile=cc("seuratQC",args$PROJNAME,plotNo(),"GeneUMAPs_%03d.png")
pngCairo(pfile,width=11,height=8.5)
pp1=paginatePlots(pp,2,3,oneLegend=F)
pp2=map(pgL,paginatePlots,2,2,oneLegend=F)
print(pp1)
print(pp2)
dev.off()

mergePNGs(pfile)
