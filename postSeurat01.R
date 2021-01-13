suppressPackageStartupMessages({
    library(Seurat)
    library(cowplot)
    library(ggplot2)
    library(patchwork)
    library(tidyverse)
})
source("seuratTools.R");

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

pfile=cc("pltSeuratV1",d10X.integrated@project.name,"Pipe:",pipeTag,"%02d.pdf")
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
