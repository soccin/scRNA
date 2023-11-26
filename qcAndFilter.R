suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(gridExtra)
})

source(file.path(SDIR,"rsrc","plot_cell_qc.R"))
source(file.path(SDIR,"rsrc","plot_nCount_vs_nFeature.R"))
source(file.path(SDIR,"rsrc","plot_highest_exprs_genes.R"))

get_mad3Cutoff<-function(xx,lower=T) {
    if(lower) {
        exp(median(log(xx))-3*mad(log(xx)))
    } else {
        median(xx)+3*mad(xx)
    }
}

apply_filter01 <- function(so,filters) {

    subset(so,
           subset = nFeature_RNA > filters$MIN_FEATURE_RNA
                    & nCount_RNA > filters$MIN_NCOUNT_RNA
                    & percent.mt < filters$PCT_MITO,
           return.null=TRUE
        )

}

get_qc_tables <- function(md,filters) {

    MIN_NCOUNT_RNA=filters$MIN_NCOUNT_RNA
    MIN_FEATURE_RNA=filters$MIN_FEATURE_RNA
    PCT_MITO=filters$PCT_MITO
    METHOD=filters$METHOD

    tbl0=tibble(
        Feature=c("nCount_RNA","nFeature_RNA","percent.mt"),
        Cutoff=c(MIN_NCOUNT_RNA,MIN_FEATURE_RNA,PCT_MITO),
        Method=c(METHOD)
    )

    tbl1=md %>%
        count(
            Count.RNA=nCount_RNA>MIN_NCOUNT_RNA,
            Num.Features=nFeature_RNA>MIN_FEATURE_RNA,
            PCT.MT=percent.mt<PCT_MITO
            ) %>%
        mutate(PCT=round(100*n/sum(n),1))

    nFail=c(
        sum(md$nCount_RNA<=MIN_NCOUNT_RNA),
        sum(md$nFeature_RNA<=MIN_FEATURE_RNA),
        sum(md$percent.mt>=PCT_MITO)
    )

    pctFail=c(
        mean(md$nCount_RNA<=MIN_NCOUNT_RNA),
        mean(md$nFeature_RNA<=MIN_FEATURE_RNA),
        mean(md$percent.mt>=PCT_MITO)
    )

    tbl2=tibble(
        Feature=c("nCount_RNA","nFeature_RNA","percent.mt"),
        N.Fail=nFail,PCT.Fail=pctFail
    )


    tbl3a=md %>% group_by(SampleID) %>% summarize(Fail.Feature_RNA=mean(nFeature_RNA<=MIN_FEATURE_RNA),Fail.Count_RNA=mean(nCount_RNA<=MIN_NCOUNT_RNA),Fail.Mito=mean(percent.mt>=PCT_MITO))
    tbl3b=md %>% group_by(SampleID) %>% summarize(Fail.Feature_RNA=sum(nFeature_RNA<=MIN_FEATURE_RNA),Fail.Count_RNA=sum(nCount_RNA<=MIN_NCOUNT_RNA),Fail.Mito=sum(percent.mt>=PCT_MITO))

    list(tbl0,tbl1,tbl2,tbl3a,tbl3b)

}

plot_qTables <- function(qTbls) {

    ptb=list()
    ptb[[len(ptb)+1]]=ggplot()+theme_void()+annotation_custom(tableGrob(mutate_if(qTbls[[1]],is.numeric,~round(.,2)),rows=NULL))
    ptb[[len(ptb)+1]]=ggplot()+theme_void()+annotation_custom(tableGrob(mutate_if(qTbls[[2]],is.numeric,~round(.,2)),rows=NULL))
    tt=qTbls[[3]] %>% mutate(PCT.Fail=round(PCT.Fail*100,2))
    ptb[[len(ptb)+1]]=ggplot()+theme_void()+annotation_custom(tableGrob(mutate_if(tt,is.numeric,~round(.,2)),rows=NULL))
    tt=qTbls[[4]] %>% mutate_if(is.numeric,\(x)x*100)
    ptb[[len(ptb)+1]]=ggplot()+theme_void()+annotation_custom(tableGrob(mutate_if(tt,is.numeric,~round(.,2)),rows=NULL))
    ptb[[len(ptb)+1]]=ggplot()+theme_void()+annotation_custom(tableGrob(mutate_if(qTbls[[5]],is.numeric,~round(.,2)),rows=NULL))

    ptb

}

qcSamples <- function(so) {

    sampleId=unique(so@meta.data$orig.ident)
    n.samples=len(sampleId)
    if(n.samples>1) {
        sampleId=cc("MERGE",so@project.name)
    }

    md=tibble(so@meta.data)

    pq=list()
    pq[[len(pq)+1]]=plot_cell_qc(md,"nFeature_RNA")
    pq[[len(pq)+1]]=plot_cell_qc(md,"nCount_RNA")
    pq[[len(pq)+1]]=plot_cell_qc(md,"percent.mt",lower=F,maxVal=50)

    p2=plot_nCount_vs_nFeature(so,ncol=3,nrow=2)

    pq=c(pq,p2)

    pq[[len(pq)+1]]=plot_highest_exprs_genes(so)
    pq[[len(pq)+1]]=plot_highest_exprs_genes(so,filterRibo=T)

    MIN_NCOUNT_RNA=get_mad3Cutoff(md$nCount_RNA)
    MIN_FEATURE_RNA=get_mad3Cutoff(md$nFeature_RNA)
    PCT_MITO=get_mad3Cutoff(md$percent.mt,lower=F)

    filters=list(
        MIN_NCOUNT_RNA=MIN_NCOUNT_RNA,
        MIN_FEATURE_RNA=MIN_FEATURE_RNA,
        PCT_MITO=PCT_MITO
    )

    qTbls=get_qc_tables(md,filters)

    ptb=plot_qTables(qTbls)

    pq[[len(pq)+1]]=(ptb[[1]]+ptb[[2]]+ptb[[3]])/(ptb[[4]]+ptb[[5]])

    list(plts=pq, stats=qTbls)

}

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
