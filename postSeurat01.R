suppressPackageStartupMessages({
    library(Seurat)
    library(cowplot)
    library(ggplot2)
    library(patchwork)
    library(tidyverse)
})
source("seuratTools.R");
source("tools.R")

load(tail(dir(pattern="CHECK.*Rdata"),1))

pipeTag=basename(getwd())

pltTitle=ggtitle(paste("SeuratV4",d10X.integrated@project.name,pipeTag))

s1=d10X.integrated
s1@meta.data$sample=gsub("\\d$","",s1@meta.data$orig.ident)

p1=DimPlot(s1, reduction = "umap", label=T)
p2=DimPlot(s1, reduction = "umap", group.by = "sample")

p3=ggplot(s1@meta.data,aes(sample,group=seurat_clusters)) +
    geom_bar(aes(y=..prop..,fill=factor(..x..))) +
    facet_wrap(~seurat_clusters) +
    scale_y_continuous(labels=scales::percent) +
    labs(y="Percent",fill="Sample") +
    scale_fill_discrete(labels=sort(unique(s1@meta.data$sample))) +
    geom_text(aes( label = scales::percent(..prop..,accuracy=2), y= ..prop.. ), stat= "count", vjust = 1.2, color="white")

p4=ggplot(s1@meta.data,aes(seurat_clusters,group=sample)) + geom_bar(aes(y=..prop..,fill=factor(..x..))) + facet_wrap(~sample) + scale_y_continuous(labels=scales::percent) + labs(y="Percent",fill="Cluster") + scale_fill_discrete(labels=sort(unique(s1@meta.data$seurat_clusters)))

pfile=cc("pltSeuratV1",d10X.integrated@project.name,"Pipe",pipeTag,"%02d.png")
png(filename=pfile,
    type="cairo",
    units="in",
    width=11,
    height=8.5,
    pointsize=12,
    res=150)

print(p1+pltTitle)
print(p2+pltTitle)
print(p3+pltTitle)
print(p4+pltTitle)

dev.off()

so=d10X.integrated
DefaultAssay(so)="RNA"
so=NormalizeData(so)
so=FindVariableFeatures(so)
so=ScaleData(so)
so@meta.data$sample=gsub("\\d$","",so@meta.data$orig.ident)

genes=scan("markerGenes01","")

cat("\n  Missing Genes =",paste(setdiff(genes,rownames(so@assays$RNA@data))),"\n\n")

genes=intersect(genes,rownames(so@assays$RNA@data))

fp=FeaturePlot(so,features=genes,max.cutoff="q98",min.cutoff="q9",combine=F)


pfile2=cc("pltSeuratV1","FeaturePlot",d10X.integrated@project.name,"Pipe",pipeTag,"%02d.png")

pp=paginatePlots(fp,2,2)

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


pfile3=cc("pltSeuratV1","DotVlnPlot",d10X.integrated@project.name,"Pipe",pipeTag,".pdf")

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
