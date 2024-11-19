#d10X=readRDS(args$PASS1.RDAFile)
#halt("CHECK GENE FILTER")

if(!is.null(args$GENE_FILTER)) {
    rna=d10X[[1]]@assays$RNA
    allGenes=rownames(rna@counts)
    genesToFilter=scan(args$GENE_FILTER,"")
    genesToKeep=setdiff(allGenes,genesToFilter)
}

for(ii in seq(d10X)) {
    print(ii)

    ret=doQCandFilter(d10X[[ii]], ap$MIN_NCOUNT_RNA, ap$MIN_FEATURE_RNA, ap$PCT_MITO)

    if(!is.null(args$GENE_FILTER)) {
        dn=ret$so
        dn@assays$RNA@counts=dn@assays$RNA@counts[genesToKeep,]
        dn@assays$RNA@data=dn@assays$RNA@data[genesToKeep,]
        dn@assays$RNA@meta.features=dn@assays$RNA@meta.features[genesToKeep,]
        d10X[[ii]]=dn
    } else {
        d10X[[ii]]=ret$so
    }

}
