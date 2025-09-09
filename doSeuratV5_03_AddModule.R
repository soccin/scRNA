suppressPackageStartupMessages(require(stringr))

usage="
usage: doSeuratV5_03_AddModule.R [CLUSTER_RES=res] PARAMS_2b.yaml ModuleFile

    CLUSTER_RES     optional: resolution of clusters to use for cluster level assigments
    PARAMS_2b.yaml  parameter file from pass2b
    ModuleFile      File with Module Genes (name of module from filename)

"

cArgs=commandArgs(trailing=T)
#
# Separate out any options arguments
#
optionals=grep("=",cArgs,value=T)

oArgs=list(CLUSTER_RES=NULL)
if(len(optionals)>0) {
    require(stringr, quietly = T, warn.conflicts=F)
    parseArgs=str_match(optionals,"(.*)=(.*)")
    aa=apply(parseArgs,1,function(x){oArgs[[str_trim(x[2])]]<<-str_trim(x[3])})
}

cArgs=grep("=",cArgs,value=T,invert=T)

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
moduleTag=basename(moduleFile) %>% tools::file_path_sans_ext()

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
    if("GeneList" %in% readxl::excel_sheets(moduleFile)) {
        moduleTbl=read_xlsx(moduleFile,sheet="GeneList")
    } else {
        moduleTbl=read_xlsx(moduleFile)
    }
} else if(modFileExt=="csv") {
    moduleTbl=read_csv(moduleFile)
}

modules=split(moduleTbl$Genes,moduleTbl$Module)

s1=AddModuleScore(s1,features=modules,name="Modules")

if(!is.null(oArgs$CLUSTER_RES)) {
    clusterRes=oArgs$CLUSTER_RES
} else {
    clusterRes="0.1"
}
clusterRes=grep(paste0(clusterRes,"$"),colnames(s1@meta.data),value=T)

cat("\nPlot modules ...")
pm=list()
pn=list()
po=list()
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

    po[[ii]] = ggplot(s1@meta.data,aes_string("SampleID",modTag,fill="SampleID")) +
        geom_violin() +
        theme_light() +
        geom_jitter(alpha=.1,size=.7,width=.2) +
        ylab("Module Score") +
        ggtitle(names(modules)[ii]) +
        theme(legend.position = "none")


}

cat(" done\n\n")
pfile=get_plot_filename(plotNo(),moduleTag,"ModuleScores_%03d.png")
pngCairo(pfile,width=11,height=8.5)
print(paginatePlots(pm,2,2,FALSE))
dev.off()
mergePNGs(pfile)

pfile=get_plot_filename(plotNo(),moduleTag,"ModuleDistribution_%03d.png")
pngCairo(pfile,width=8.5,height=11)
print(paginatePlots(pn,3,1,FALSE))
dev.off()
mergePNGs(pfile)

#
# Dump metadata
#

md=s1@meta.data %>% data.frame %>% rownames_to_column("CellID") %>% tibble
moduleCols=grep("^Modules\\d+",names(md))
colnames(md)[moduleCols]=paste0("mod.",names(modules))
write_csv(md,cc("metaData",moduleTag,"AddModules.csv"))










