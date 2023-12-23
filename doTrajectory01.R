library(dplyr)
library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(patchwork)

params=yaml::read_yaml("pass_02b_PARAMS.yaml")
so=readRDS(params$PASS2b.RDAFile)

cds=as.cell_data_set(so)

#
# cluster_method = c("leiden", "louvain"),
#

#cds=cluster_cells(cds)
cds=cluster_cells(cds,resolution=0.0001)

s1=as.Seurat(cds,assay=NULL)

cds=learn_graph(cds,use_partition = FALSE, verbose = TRUE)

#cdsA=order_cells(cds,root_cells=colnames(cds[,colData(cds)$orig.ident=="WT_DSMO"]))
cdsA=order_cells(cds,root_cells=colnames(cds[,colData(cds)$SCT_snn_res.0.2==10]))

pg0=plot_cells(cdsA,color_cells_by = "cluster",label_cell_groups = FALSE,label_groups_by_cluster=FALSE,label_leaves=FALSE,label_branch_points=FALSE,label_roots = FALSE,trajectory_graph_color = "grey25", cell_size=1)
pg1=plot_cells(cdsA,color_cells_by = "pseudotime",group_cells_by = "cluster",label_cell_groups = FALSE,label_groups_by_cluster=FALSE,label_leaves=FALSE,label_branch_points=FALSE,label_roots = FALSE,trajectory_graph_color = "grey60", cell_size=1)

#s1=as.Seurat(cdsA,assay=NULL)

s2=AddMetaData(s1,pseudotime(cdsA),"pseudotime")

ps1=FeaturePlot(s2,"pseudotime",pt.size=2)

pl=FeaturePlot(s2,"pseudotime",split.by="SampleID",pt.size=2,combine=F)

ps2=Reduce(`+`,pl) + plot_layout(guide="collect")

pfile=cc("monocle3","pseudotime","cres_Default_RootCls10",".pdf")
#pfile=cc("monocle3","pseudotime","cres_Default",".pdf")
#pfile=cc("monocle3","pseudotime","cres_0.0004",".pdf")
pdf(file=pfile,height=8,width=8)
print(pg0)
print(pg1)
print(ps1)
print(ps2)
dev.off()

