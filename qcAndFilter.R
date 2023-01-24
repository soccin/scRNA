suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(gridExtra)
})

doQCandFilter <- function(so,MIN_NCOUNT_RNA,MIN_FEATURE_RNA,PCT_MITO) {

    sampleId=unique(so@meta.data$orig.ident)
    n.samples=len(sampleId)
    if(n.samples>1) {
        sampleId=cc("MERGE",so@project.name)
    }

    pg0=VlnPlot(so, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size=0,combine=F)

    pg0[[1]]=transformVlnPlot(pg0[[1]],MIN_FEATURE_RNA)
    pg0[[2]]=transformVlnPlot(pg0[[2]],MIN_NCOUNT_RNA)
    pg0[[3]]=transformVlnPlot(pg0[[3]],PCT_MITO,50)

    max.nCount_RNA=max(so@meta.data$nCount_RNA)
    max.nFeature_RNA=max(so@meta.data$nFeature_RNA)
    max.pct.mito=max(so@meta.data$percent.mt)

    tbl=so@meta.data %>% tibble %>%
        count(
            Count.RNA=nCount_RNA>MIN_NCOUNT_RNA,
            Num.Features=nFeature_RNA>MIN_FEATURE_RNA,
            PCT.MT=percent.mt<PCT_MITO
            ) %>%
        mutate(PCT=round(100*n/sum(n),1))

    tbl2=so@meta.data %>% tibble %>%
        mutate(
            nCount_RNA.FAIL=nCount_RNA>MIN_NCOUNT_RNA,
            nFeature_RNA.FAIL=nFeature_RNA>MIN_FEATURE_RNA,
            percent.mt.FAIL=percent.mt<PCT_MITO
            ) %>%
        group_by(SampleID) %>%
        summarize_if(is.logical,~round(100*(1-sum(.)/n()),1))

    plot1 <- FeatureScatter(so, feature1 = "nCount_RNA", feature2 = "percent.mt") +
                geom_hline(yintercept=PCT_MITO,col="grey",alpha=0.75) +
                geom_vline(xintercept=MIN_NCOUNT_RNA,col="grey",alpha=0.75) +
                annotation_custom(tableGrob(tbl,rows=NULL),
                    xmin=max.nCount_RNA/2,ymin=50,xmax=max.nCount_RNA/2,ymax=50) +
                scale_x_continuous(breaks=sort(c(MIN_NCOUNT_RNA,seq(0,1e6,by=20000)))) +
                scale_y_continuous(breaks=sort(c(PCT_MITO,seq(0,100,by=20))))

    plot2 <- FeatureScatter(so, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
                geom_hline(yintercept=MIN_FEATURE_RNA,col="grey",alpha=0.75) +
                geom_vline(xintercept=MIN_NCOUNT_RNA,col="grey",alpha=0.75) +
                scale_x_continuous(breaks=sort(c(MIN_NCOUNT_RNA,seq(0,1e6,by=20000)))) +
                scale_y_continuous(breaks=sort(c(MIN_FEATURE_RNA,seq(0,1e5,by=2000)))) +
                annotation_custom(tableGrob(tbl2,rows=NULL),xmin=max.nCount_RNA/3,ymin=max.nFeature_RNA/10)

    keep=(
        so@meta.data$nFeature_RNA > MIN_FEATURE_RNA &
        so@meta.data$nCount_RNA > MIN_NCOUNT_RNA &
        so@meta.data$percent.mt < PCT_MITO
        )

    knitr::kable(table(keep))

    so <- subset(so,
            subset = nFeature_RNA > MIN_FEATURE_RNA & nCount_RNA > MIN_NCOUNT_RNA & percent.mt < PCT_MITO
            )

    #pg1=VlnPlot(so, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

    list(so=so, plts=list(pg0,plot1,plot2), stats=keep)

}
