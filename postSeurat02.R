suppressPackageStartupMessages({
    library(tidyverse)
    library(Seurat)
    library(patchwork)
    library(pals)
})

source("seuratTools.R");
source("tools.R")

pipeTag=paste0(basename(getwd()),"_v2")

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

source("getCellTypesGenesLists.R")
features=getCellTypeGeneLists(rownames(d10X.integrated@assays$RNA@counts))

defaultAssay="Integrated"
defaultAssay="RNA"

if(defaultAssay=="RNA") {
    so=d10X.integrated
    DefaultAssay(so)="RNA"
    so=NormalizeData(so)
    so=FindVariableFeatures(so)
    so=ScaleData(so)
    so@meta.data$sample=gsub("\\d$","",so@meta.data$orig.ident)

    xx=AddModuleScore(object=so,features=features,name="CellType.Score.")
} else {
    xx=AddModuleScore(object=d10X.integrated,features=features,name="CellType.Score.")
}

# ctScores=xx@meta.data %>% rownames_to_column("CellID") %>% tibble %>% select(CellID,matches("\\.Score\\."))
# scrCols=grep("\\.Score\\.",colnames(ctScores),value=T)

cellTypes=xx@meta.data %>%
    rownames_to_column("CellID") %>%
    tibble %>% select(CellID,matches("\\.Score\\.")) %>%
    gather(Score,Value,-CellID) %>%
    group_by(CellID) %>%
    arrange(desc(Value)) %>%
    slice(1)

newMetaData=xx@meta.data %>%
    rownames_to_column("CellID") %>%
    tibble %>%
    left_join(cellTypes) %>%
    select(-matches("CellType.Score.")) %>%
    mutate(CellType=names(features)[gsub(".*Score.","",Score) %>% as.numeric]) %>%
    mutate(CellType=ifelse(Value>0,CellType,"UNK")) %>%
    mutate(sample=gsub("[.0-9]+$","",orig.ident)) %>%
    mutate(sample=factor(sample,levels=c("Sham","RT","Tmem.RT"))) %>%
    data.frame %>%
    column_to_rownames("CellID")

if(defaultAssay=="Integrated") {
    so=d10X.integrated
}

so@meta.data=newMetaData

nCellTypes=len(unique(cellTypes$Score))

clsPal01=rev(kelly(22))
pal=c("#BEBEBE66",clsPal01[1:nCellTypes])

#ctCls=RColorBrewer::brewer.pal(10,"Paired")
#ctCls[len(ctCls)]="#EBEBEB"

ctCls=rev(pal)

pg3=DimPlot(so,label=T)
pg3b=DimPlot(so,label=T,split.by="sample",ncol=2)
pg4=DimPlot(so,group.by=c("CellType"),cols=ctCls)
pg4b=DimPlot(so,group.by=c("CellType"),cols=ctCls,split.by="sample",ncol=2)

samples=distinct(so@meta.data,sample) %>% pull(sample)
pgN=list()
for(si in samples) {
    print(si)
    so.ss=subset(so,subset=sample==si)
    pp=DimPlot(so.ss,label=T)+ggtitle(paste("Sample:",si))
    pgN[[len(pgN)+1]]=pp
}
for(si in samples) {
    print(si)
    so.ss=subset(so,subset=sample==si)
    pp=DimPlot(so.ss,group.by=c("CellType"),cols=ctCls)+ggtitle(paste("Sample:",si))
    pgN[[len(pgN)+1]]=pp
}


clusterTypeTbl=tibble(newMetaData) %>%
    count(seurat_clusters,CellType) %>%
    group_by(seurat_clusters) %>%
    mutate(PCT=n/sum(n))

cellTypes=clusterTypeTbl %>% ungroup %>% distinct(CellType) %>% pull %>% sort

clusterTypeTbl=clusterTypeTbl %>%
    mutate(CellType=factor(CellType,levels=rev(cellTypes)))

superMajority=clusterTypeTbl %>% arrange(desc(PCT)) %>% filter(PCT>.6)

pass1=superMajority %>% pull(seurat_clusters)

names(pal)=levels(cellTypes)

pgA=clusterTypeTbl %>% ggplot(aes(seurat_clusters,PCT,fill=CellType))+geom_bar(stat="identity") + scale_fill_manual(values=pal,name="Cell Types",drop=F) + theme_bw(base_size=18)+ scale_y_continuous(labels=scales::percent) + ggtitle("All Clusters")
pgB=clusterTypeTbl %>% filter(!(seurat_clusters %in% pass1)) %>% ggplot(aes(seurat_clusters,PCT,fill=CellType))+geom_bar(stat="identity") + scale_fill_manual(values=pal,name="Cell Types",drop=F) + theme_bw(base_size=18)+ scale_y_continuous(labels=scales::percent) + ggtitle("Mixed Clusters")
pgC=clusterTypeTbl %>% filter(seurat_clusters %in% pass1) %>%ggplot(aes(seurat_clusters,PCT,fill=CellType))+geom_bar(stat="identity") + scale_fill_manual(values=pal,name="Cell Types",drop=F) + theme_bw(base_size=18)+ scale_y_continuous(labels=scales::percent) + ggtitle("Super Majority Clusters")

pfile4=cc("pltSeuratV1","CellTypeAutoId","DefAssay",defaultAssay,d10X.integrated@project.name,"Pipe",pipeTag,".pdf")
pdf(file=pfile4,width=14,height=8.5)
print(pg3+pg4)
print(pg3b)
print(pg4b)
print(pgN)
print(pgA)
print(pgB)
print(pgC)
dev.off()

