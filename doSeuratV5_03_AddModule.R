suppressPackageStartupMessages(require(stringr))

usage="
usage: doSeuratV5_03_FeaturePlots.R PARAMS_2b.yaml ModuleFile

    PARAMS_2b.yaml  parameter file from pass2b
    ModuleFile      File with Module Genes (name of module from filename)

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

plotNo<-makeAutoIncrementor(50)

moduleFile=cArgs[2]

if(!grepl("\\.(xlsx|csv)$",moduleFile)) {
    cat("\n    Note implemented: Only XLSX modules files currently working\n\n")
    quit()
}

args$MODULES=normalizePath(moduleFile)
modFileExt=tools::file_ext(moduleFile)

suppressPackageStartupMessages({
    library(Seurat)
    library(patchwork)
    library(tidyverse)
    library(readxl)
    library(openxlsx)
    library(pals)
})

##########################################################################
#
# INCLUDE BREAK
#
##########################################################################

s1=readRDS(args$PASS2b.RDAFile)
DefaultAssay(s1)="SCT"

if(modFileExt=="xlsx") {
    moduleTbl=read_xlsx(moduleFile,sheet="GeneList")
} else if(modFileExt=="csv") {
    moduleTbl=read_csv(moduleFile)
}

modules=split(moduleTbl$Genes,moduleTbl$Module)

s1=AddModuleScore(s1,features=modules,name="Modules")

clusterRes="SCT_snn_res.0.1"

cat("\nPlot modules ...")
pm=list()
pn=list()
for(ii in seq(len(modules))) {
    print(ii)
    modTag=paste0("Modules",ii)
    pp=FeaturePlot(s1,
        features=modTag,
        max.cutoff="q95",min.cutoff="q05",
        combine=F)

    pm[[ii]]=pp[[1]] + ggtitle(names(modules)[ii])

    pn[[ii]] = ggplot(s1@meta.data,aes_string(clusterRes,modTag,fill=clusterRes)) +
        geom_violin() +
        theme_light() +
        geom_jitter(alpha=.1,size=.7,width=.2) +
        ylab("Module Score") +
        ggtitle(names(modules)[ii]) +
        theme(legend.position = "none")


}

cat(" done\n\n")
pfile=cc("seuratQC",args$PROJNAME,plotNo(),"ModuleScores_%03d.png")
pngCairo(pfile,width=11,height=8.5)
print(paginatePlots(pm,2,2,FALSE))
dev.off()
mergePNGs(pfile)

pfile=cc("seuratQC",args$PROJNAME,plotNo(),"ModuleDistribution_%03d.png")
pngCairo(pfile,width=8.5,height=11)
print(paginatePlots(pn,3,1,FALSE))
dev.off()
mergePNGs(pfile)














