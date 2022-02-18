suppressPackageStartupMessages(require(stringr))

usage="
usage: doSeuratV5_03_FeaturePlots.R [CRES=clusterResolution] PARAMS_2b.yaml ModuleFile

    PARAMS_2b.yaml   parameter file from pass2b
    ModuleFile       File with Module Genes (name of module from filename)

    Optional:
      CRES=resNumber Cluster Resolution to use (eg: CRES=0.2)

"

# cArgs=commandArgs(trailing=T)

# if(len(cArgs)!=2) {
#     cat(usage)
#     quit()
# }

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

cArgs=commandArgs(trailing=T)
args=list(CRES=NULL)
usage=str_interp(usage,args)

ii=grep("=",cArgs)
if(len(ii)>0) {
    parseArgs=str_match(cArgs[ii],"(.*)=(.*)")
    aa=apply(parseArgs,1,function(x){args[[str_trim(x[2])]]<<-str_trim(x[3])})
}

argv=grep("=",cArgs,value=T,invert=T)

if(len(argv)!=2) {
    cat(usage)
    quit()
}

suppressPackageStartupMessages(library(yaml))
args0=read_yaml(argv[1])
args=c(args0,args)

glbs=args$glbs
ap=args$algoParams

plotNo<-makeAutoIncrementor(50)

moduleFile=argv[2]

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

if(is.null(args$CRES)) {
    clusterResolutions=c("0.1")
} else {
    clusterResolutions=strsplit(args$CRES,",")[[1]]
}
cTag=gsub("_res.*","_res",grep("_res",colnames(s1@meta.data),value=T)[1])

halt()

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

    for(ci in clusterResolutions) {

        clusterRes=paste0(cTag,".",ci)

        pn[[len(pn)+1]] = ggplot(s1@meta.data,aes_string(clusterRes,modTag,fill=clusterRes)) +
            geom_violin() +
            theme_light() +
            geom_jitter(alpha=.1,size=.7,width=.2) +
            ylab("Module Score") +
            ggtitle(names(modules)[ii]) +
            theme(legend.position = "none")

    }

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

#
# Dump metadata
#

md=s1@meta.data %>% data.frame %>% rownames_to_column("CellID") %>% tibble
moduleCols=grep("^Modules\\d+",names(md))
colnames(md)[moduleCols]=paste0("mod.",names(modules))
write_csv(md,"metaData_AddModules.csv")










