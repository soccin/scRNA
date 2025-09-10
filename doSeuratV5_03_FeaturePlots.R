usage="
usage: doSeuratV5_03_FeaturePlots.R [CRES=clusterResolution] PARAMS_2b.yaml GeneListFile

    PARAMS_2b.yaml  parameter file from pass2b
    GeneListFile    List of Genes to Plot

  Optional:
    CRES=resNumber  Cluster Resolution to use (eg: CRES=0.2)

"
STAGE=3

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

cArgs=commandArgs(trailing=T)
args=list(CRES=NULL)
usage=str_interp(usage,args)

if(len(cArgs)<2) {
    cat(usage)
    quit()
}

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
DefaultAssay(s1)="SCT"

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

nSamples=len(unique(s1@meta.data$orig.ident))
if(nSamples<=4) {
    splot.nrow=2
    splot.ncol=2
} else if(nSamples <=6) {
    splot.nrow=2
    splot.ncol=3
} else {
    splot.nrow=3
    splot.ncol=3
}

pfile=get_plot_filename(plotNo(),"GeneUMAPs",basename(geneListFile),"%03d.png")
pngCairo(pfile,width=11,height=8.5)
pp1=paginatePlots(pp,2,2,oneLegend=F)
pp2=map(pgL,paginatePlots,splot.nrow,splot.ncol,oneLegend=T)
print(pp1)
print(pp2)
dev.off()

mergePNGs(pfile)
fs::dir_ls(dirname(pfile),regex=basename(gsub("_%03d.*",".*.png$",pfile))) %>% map(fs::file_delete)

if(!is.null(args$CRES)) {

    library(pals)
    pal1=c(cols25(),brewer.dark2(8))

    clusterRes=grep(args$CRES,grep("res",colnames(s1@meta.data),value=T),value=T)

    pv=VlnPlot(s1,features=genes,combine=F,group.by=clusterRes,col=pal1)
    pp3=paginatePlots(pv,2,3)

    pfile=get_plot_filename(plotNo(),"GeneVlnPlt",
                basename(geneListFile),"res",args$CRES,"%03d.png")
    pngCairo(pfile,width=11,height=8.5)
    print(pp3)
    dev.off()
    mergePNGs(pfile)

}
