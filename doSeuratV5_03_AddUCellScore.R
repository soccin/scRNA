suppressPackageStartupMessages(require(stringr))

usage="
usage: doSeuratV5_03_AddUCellScore.R [CLUSTER_RES=res] PARAMS_2b.yaml ModuleFile

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

plotNo<-makeAutoIncrementor(70)

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
    library(UCell)
})

##########################################################################
#
# INCLUDE BREAK
#
##########################################################################

s1=readRDS(args$PASS2b.RDAFile)
genes=rownames(s1@assays[[DefaultAssay(s1)]])

cat("Using SCT assay\n")
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
names(modules)=make.names(names(modules))

TAG="Intersect"
if(TAG=="Intersect") {
    modules0=modules
    modules=map(modules,\(x){intersect(x,genes)})
}

if(!is.null(oArgs$CLUSTER_RES)) {
    clusterRes=oArgs$CLUSTER_RES
} else {
    clusterRes="0.1"
}
clusterRes=grep(paste0(clusterRes,"$"),colnames(s1@meta.data),value=T)

s1@meta.data[[clusterRes]]=factor(as.numeric(s1@meta.data[[clusterRes]]))

#
# 
#

s1=AddModuleScore_UCell(s1,features=modules,name="_UC")
uModules=grep("_UC$",colnames(s1@meta.data),value=T)
s1=SmoothKNN(s1,signature.names=uModules,reduction="umap")

UCellCutoff=0.2

cat("\nPlot modules ...")
pm=list()
pn=list()
for(ii in seq(len(uModules))) {
    print(ii)
    modTag=gsub("_UC$","",uModules[ii])
    pp=FeaturePlot(s1,
        features=cc(uModules[ii],"kNN"),
        max.cutoff=UCellCutoff,min.cutoff=0,
        cols=c("grey90","darkred"),
        combine=F)

    pm[[ii]]=pp[[1]] + ggtitle(modTag)

    pn[[ii]] = ggplot(s1@meta.data,aes_string(clusterRes,uModules[ii],fill=clusterRes)) +
        geom_violin() +
        theme_light() +
        geom_jitter(alpha=.1,size=.7,width=.2) +
        ylab("Module Score") +
        ggtitle(names(modules)[ii]) +
        theme(legend.position = "none")


}

STAGE=7
cat(" done\n\n")
pfile=get_plot_filename(plotNo(),cc("UCellModuleScores",TAG,"%03d.png"))
pngCairo(pfile,width=11,height=8.5)
print(paginatePlots(pm,2,2,FALSE))
dev.off()
mergePNGs(pfile)

pfile=get_plot_filename(plotNo(),cc("UCellModuleDistribution",TAG,"%03d.png"))
pngCairo(pfile,width=8.5,height=11)
print(paginatePlots(pn,3,1,FALSE))
dev.off()
mergePNGs(pfile)

fs::file_delete(fs::dir_ls("results/stage7",regex=".png"))

mTbl=s1@meta.data %>%
    rownames_to_column("CellID") %>%
    tibble %>%
    select(CellID,SampleID,all_of(uModules)) %>%
    gather(Module,Score,all_of(uModules))
mS=mTbl %>% group_by(CellID) %>% slice_max(Score) %>% slice(1)

pTbl=list()
thetas=c(0.05,0.1,0.15,0.2)
for(theta in thetas) {
    pTbl[[len(pTbl)+1]]=mS %>%
        group_by(SampleID) %>%
        summarize(PCT.Unknown=mean(Score<=theta)) %>%
        mutate(Threshold=theta)
}
pTbl=bind_rows(pTbl) %>% spread(SampleID,PCT.Unknown)

tfile=get_plot_filename(plotNo(),cc("UCellPCTUnknown",TAG,".xlsx"))
openxlsx::write.xlsx(pTbl,tfile)

