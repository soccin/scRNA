plot_cell_qc <- function(obj,feature,maxVal=NA,lower=T,numCellsPerSample=5000) {

    if(class(obj)[1]=="Seurat") {
        md=tibble(so@meta.data)
    } else {
        md=obj
    }

    nSamples=len(unique(md$SampleID))
    if(is.na(maxVal)) {
        maxVal=unname(quantile(md[[feature]],.95))
    }

    numCells=numCellsPerSample*nSamples

    #
    # from:
    #    http://bioconductor.org/books/3.16/OSCA.basic/quality-control.html#identifying-low-quality-cells
    #

    medFeat=median(md[[feature]])
    if(lower) {
        cut3mad=exp(log(medFeat)-3*mad(log(md[[feature]])))
    } else {
        cut3mad=medFeat+3*mad(md[[feature]])
        maxVal=max(maxVal,cut3mad)
    }

    smd=md %>% filter(.data[[feature]]<maxVal) %>% sample_n(min(nrow(md),numCells))
    yBreaks=sort(c(scales::breaks_pretty(3)(smd[[feature]]),medFeat,cut3mad))
    yBreaksLabels=sort(c(scales::breaks_pretty(3)(smd[[feature]]),
        round(medFeat,1),round(cut3mad,1)))

    cat(yBreaks,"\n")

    md %>%
        ggplot(aes_string("SampleID",feature,fill="SampleID")) +
        cowplot::theme_half_open() +
        NoLegend() +
        theme(axis.text.x = element_text(angle = 45, hjust=1)) +
        coord_flip() +
        rasterize(geom_jitter(alpha=0.1,size=0.25,data=smd)) +
        scale_y_continuous(limit=c(0,maxVal),breaks=yBreaks,labels=yBreaksLabels) +
        geom_violin(alpha=.85) +
        geom_hline(yintercept=c(medFeat,cut3mad),color=c("black","grey75"),size=1,alpha=.5) +
        geom_boxplot(fill="grey50",alpha=.2,outlier.shape=NA,coef=0,width=1/10) +
        xlab("") + ylab("") +
        ggtitle(paste(feature))

#        scale_y_continuous(trans="log1p",limit=c(0,max(md[[feature]])),breaks=yBreaks,labels=yBreaksLabels) +


}
