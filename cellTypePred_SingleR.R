require(celldex)
library(SingleR)
library(Seurat)

so=readRDS("pass_02b_SObj_d0577b7ca8d3271c29bdc6164869267a_s1_.rda")
so.sce=as.SingleCellExperiment(so)
mcells=celldex::MouseRNAseqData()

so.test=so.sce[,1:100]
system.time({pred=SingleR(test=so.test,ref=mcells,assay.type.test=1,labels=mcells$label.main,prune=T)})

library(tidyverse)
predT=data.frame(pred) %>% rownames_to_column("CellID") %>% tibble %>% select(-matches("scores"),everything())

write_csv(predT,"pred_Mouse.csv.gz")