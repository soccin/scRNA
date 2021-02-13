source("postSeurat02.R")
features=list(NKCells=c("Nkg7", "Gzma", "Gzmb"))
features
d10X.integrated@assays$RNA
d10X.integrated@assays$RNA@meta.features
d10X.integrated@assays$RNA@counts
rownames(d10X.integrated@assays$RNA@counts)
intersect(features$NKCells)
intersect(features$NKCells,rownames(d10X.integrated@assays$RNA@counts))
AddModuleScore(object=d10X.integrated,features=features,ctrl=5,name="Features")
xx=AddModuleScore(object=d10X.integrated,features=features,ctrl=5,name="Features")
xx@meta.data
head(xx@meta.data)
ggplot(xx@meta.data,aes(seurat_clusters,Features1))
ggplot(xx@meta.data,aes(seurat_clusters,Features1)) + geom_violin()
xx=AddModuleScore(object=d10X.integrated,features=features,name="Features")
ggplot(xx@meta.data,aes(seurat_clusters,Features1)) + geom_violin()
ggplot(xx@meta.data,aes(seurat_clusters,Features1)) + geom_boxplot()
features=list(NKCells=c("Nkg7", "Gzma", "Gzmb"),BCells=c("Cd79a","Cd79b","Cd19"))

head(xx@meta.data)
ggplot(xx@meta.data,aes(seurat_clusters,Features1)) + geom_boxplot()
ggplot(xx@meta.data,aes(seurat_clusters,Features2)) + geom_boxplot()
ggplot(xx@meta.data,aes(seurat_clusters,Features2)) + geom_violin()
tibble(xx@meta.data)
tibble(xx@meta.data) %>% group_by(seurat_clusters) %>% summarize(scr1=median(Features1),scr2=median(Features2))
