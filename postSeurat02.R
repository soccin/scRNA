suppressPackageStartupMessages({
    library(tidyverse)
    library(Seurat)
    library(patchwork)
})

source("seuratTools.R");
source("tools.R")

pipeTag=basename(getwd())

d10X.integrated=readRDS("obj__d10X.integrated.rda")

md=tibble(d10X.integrated@meta.data) %>%
    mutate(sample=gsub("\\d+$","",orig.ident) %>% gsub("\\.$","",.)) %>%
    mutate(sample=factor(sample,levels=c("Sham","RT","Tmem.RT"))) %>%
    mutate(cluster=forcats::fct_relevel(seurat_clusters,rev))

tbl=md %>% count(sample,cluster) %>% group_by(sample) %>% mutate(total=sum(n)) %>% mutate(P=n/total)


numClust=nlevels(md$cluster)
colors=rev(rep(RColorBrewer::brewer.pal(10,"Paired"),ceiling(numClust/10))[1:numClust])

# pg=ggplot(tbl,aes(sample,P,fill=sample))+geom_bar(stat="identity")

pg=md %>%
    count(sample,cluster) %>%
    ggplot(aes(sample,n,fill=cluster)) +
    geom_bar(stat="identity",position="fill") +
    scale_fill_manual(values=colors) +
    labs(x="Sample",y="Percentage",fill="Cluster") +
    scale_y_continuous(labels=scales::percent)


pg0=ggplot(tbl,aes(sample,P,fill=sample)) +
    geom_bar(stat="identity") +
    scale_y_continuous(labels=scales::percent) +
    labs(x="Sample",y="Percent",fill="Sample")

pg1=pg0 + facet_wrap(~cluster) +
    geom_text(aes( sample, P, label = scales::percent(P,accuracy=0.1)), stat= "identity", vjust = 0, color="black", size=3)

pg2=pg0 + facet_wrap(~cluster,scales="free") +
    geom_text(aes( sample, P, label = scales::percent(P,accuracy=0.1)), stat= "identity", vjust = 1.2, color="white", size=3)

pfile3=cc("pltSeuratV1","percentClusterBySample",d10X.integrated@project.name,"Pipe",pipeTag,".pdf")
pdf(file=pfile3,width=11,height=8.5)
print(pg)
print(pg1)
print(pg2)
dev.off()

