library(ggrastr)
library(ggforce)

plot_nCount_vs_nFeature<-function(obj,ncol=3,nrow=2) {

    if(class(obj)[1]=="Seurat") {
        md=tibble(obj@meta.data)
    } else {
        md=obj
    }

    pg=md %>%
        ggplot(aes(nCount_RNA,nFeature_RNA,color=percent.mt)) +
        rasterize(geom_point()) +
        scale_x_log10() +
        scale_y_log10() +
        scale_colour_gradient(low = "#deebf7", high = "#08306b") +
        theme_classic() +
        stat_smooth(method=lm)

    nPages=n_pages(pg+facet_wrap_paginate(~SampleID,ncol=ncol,nrow=nrow,page=1))

    pp=list()
    for(page in seq(nPages)) {
        pp[[page]]=pg+facet_wrap_paginate(~SampleID,ncol=ncol,nrow=nrow,page=page)
    }

    pp

}