suppressPackageStartupMessages({
    library(Seurat)
    library(cowplot)
    library(ggplot2)
    library(patchwork)
    library(tidyverse)
})
source("seuratTools.R");
source("tools.R")

pipeTag=paste0(basename(getwd()),"_v2")

d10X.integrated=readRDS("obj__d10X.integrated.rda")

so=d10X.integrated
DefaultAssay(so)="RNA"
so=NormalizeData(so)
so=FindVariableFeatures(so)
so=ScaleData(so)
so@meta.data$sample=factor(
    gsub("\\d+$","",so@meta.data$orig.ident) %>% gsub("\\.$","",.),
    levels=c("Sham","RT","Tmem.RT")
    )

pltTitle=ggtitle(paste("SeuratV4",so@project.name,pipeTag))

cArgs=commandArgs(trailing=T)

genes=scan(cArgs[1],"") %>% str_to_title(.)
gTag=basename(cArgs)

if(len(intersect(genes,rownames(so@assays$RNA@counts)))!=len(genes)) {
    cat("\n\nMissing genes:\n")
    cat(paste(setdiff(genes,rownames(so@assays$RNA@counts))),"\n")
    cat("\n")
    quit()
}



pg.fp=FeaturePlot(so,features=genes,max.cutoff="q98",min.cutoff="q9",combine=F)
pNull=list(ggplot()+theme_void())

for(gi in genes) {
    pgi=FeaturePlot(so,features=gi,max.cutoff="q98",min.cutoff="q9",combine=F,split.by="sample")
    pg.fp=c(pg.fp,pgi,pNull)
}

pp=paginatePlots(pg.fp,2,2)

pfile2=cc("pltSeuratV1",gTag,"FeaturePlot",d10X.integrated@project.name,"Pipe",pipeTag,"%02d.png")
png(filename=pfile2,
    type="cairo",
    units="in",
    width=11,
    height=8.5,
    pointsize=12,
    res=300)

xx=map(pp,print)

dev.off()


gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


pfile3=cc("pltSeuratV1",gTag,"DotVlnPlot",d10X.integrated@project.name,"Pipe",pipeTag,".pdf")

dp1=DotPlot(so,features=genes,dot.scale=10)
dp2=DotPlot(so,features=genes,dot.scale=6,split.by="sample",cols=gg_color_hue(3))

vv=VlnPlot(so,features=genes,group.by="seurat_clusters",pt.size=0,combine=F)
pp1=paginatePlots(vv,3,2)

vv=VlnPlot(so,features=genes,split.by="seurat_clusters",group.by="sample",pt.size=0,combine=F)
pp2=paginatePlots(vv,3,1)

vv=VlnPlot(so,features=genes,group.by="seurat_clusters",split.by="sample",pt.size=0,combine=F)
pp3=paginatePlots(vv,3,1)

pdf(file=pfile3,width=10,height=11)
print(dp1)
print(dp2)
xx=map(c(pp1,pp2,pp3),print)
dev.off()



