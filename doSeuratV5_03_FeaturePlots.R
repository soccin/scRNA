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

suppressPackageStartupMessages(require(stringr))

usage="
usage: doSeuratV5_03_FeaturePlots.R [CRES=clusterResolution] PARAMS_2b.yaml GeneListFile

    PARAMS_2b.yaml  parameter file from pass2b
    GeneListFile    List of Genes to Plot

  Optional:
    CRES=resNumber  Cluster Resolution to use (eg: CRES=0.2)

"

cArgs=commandArgs(trailing=T)
args=list(CRES=NULL)
usage=str_interp(usage,args)

ii=grep("=",cArgs)
if(len(ii)>0) {
    parseArgs=str_match(cArgs[ii],"(.*)=(.*)")
    aa=apply(parseArgs,1,function(x){args[[str_trim(x[2])]]<<-str_trim(x[3])})
}

argv=grep("=",cArgs,value=T,invert=T)

suppressPackageStartupMessages(library(yaml))
args0=read_yaml(argv[1])
args=c(args0,args)

source(file.path(SDIR,"seuratTools.R"))
source(file.path(SDIR,"plotTools.R"))

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

geneListFile=argv[2]
genes=scan(geneListFile,"")
if(args$glbs$genome=="mm10") {
    genes=str_to_title(genes)
}


##########################################################################
#
# INCLUDE BREAK
#halt("INCLUDE")
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


pfile=cc("seuratQC",args$PROJNAME,plotNo(),"GeneUMAPs",basename(geneListFile),"%03d.png")
pngCairo(pfile,width=11,height=8.5)
pp1=paginatePlots(pp,2,3,oneLegend=F)
pp2=map(pgL,paginatePlots,2,2,oneLegend=F)
print(pp1)
print(pp2)
dev.off()

mergePNGs(pfile)
