#
# Seurat Stuff
#
       # if (is.null(x = cols)) {
       #      cols <- hue_pal()(length(x = levels(x = idents)))
       #      cols <- Interleave(cols, InvertHex(hexadecimal = cols))
       #  }

get_plot_filename()

plotFilePrefix="seuratPlot"

get_plot_filename<-function(...) {
    argv=list(...)
    odir="results/stage1"
    if(!is.null(argv$stage)) {
        odir=file.path("results",paste0("stage",argv$stage))
        argv$stage=NULL
    } else if(exists("STAGE")) {
        odir=file.path("results",paste0("stage",STAGE))  
    }
    fs::dir_create(odir)
    base=cc(plotFilePrefix,paste0(argv,collapse="_"))
    fname=file.path(odir,base)
    cat("fname =",fname,"\n")
    #fname
}

pngCairo<-function(filename,width=14,height=8.5,pointsize=12,res=150) {

    png(filename,type="cairo",units="in",
        width=width,height=height,pointsize=pointsize,res=res)

}

mergePNGs<-function(fileSpec) {
    fileRe=gsub("_%\\d+d",".*",fileSpec)
    pdfFile=gsub("_%\\d+d.*",".pdf",fileSpec)
    system2("convert",c(sort(dir_ls(regex=fileRe)),pdfFile),stderr=cc("stderr","mergePNGs","convert",DATE()))
}

paginatePlots<-function(plts.o,pRows,pCols,oneLegend=T) {

    nPlots=pRows*pCols

    pp=list()
    page=1
    currPlot=1
    pa=NULL

    if(oneLegend) {
        plts=map(plts.o,~.+theme(legend.position = "none"))
        plts[[pCols]]=plts.o[[min(pCols,len(plts.o))]]
    } else {
        plts=plts.o
    }


    for(ii in seq(len(plts))) {
        if(is.null(pa)) {
            pa=plts[[ii]]
        } else {
            pa=pa+plts[[ii]]
            if(len(pa$patches$plots)==(nPlots-1)) {
                pp[[len(pp)+1]]=pa + plot_layout(ncol=pCols,nrow=pRows)
                pa=NULL
            }
        }
    }

    if(!is.null(pa))
        pp[[len(pp)+1]]=pa + plot_layout(ncol=pCols,nrow=pRows)

    pp

}

minor.breaks.log10 <- function(major) {
    map(major,function(x){(0:9)*(10^x)}) %>% unlist
}

transformVlnPlot<-function(gg,qcFilterLevels,maxValue=NA) {

    if(is.na(maxValue)) {
        maxValue=quantile(gg$data[[1]],.95)
        cat("maxValue =",maxValue,"\n")
    }

    minValue=0

    newBreaks=scales::breaks_extended(5)(c(minValue,maxValue))
    newBreaks=unique(sort(c(newBreaks,qcFilterLevels*c(1/2,1,2))))

    xLabels=str_wrap(gsub("[-_.]"," ",levels(gg$data[[2]])),15)

    gret=gg +
        geom_jitter(alpha=0.1,size=0.25) +
        geom_hline(yintercept=c(2,1,.5)*qcFilterLevels,col="#BEBEBEBE",size=c(1,1.4,1)) +
        scale_y_continuous(limit=c(minValue,maxValue), breaks=newBreaks) +
        scale_x_discrete(labels=xLabels) +
        theme(axis.title.y=element_blank()) +
        coord_flip() +
        NoLegend()

    pct.pass=round(100*tapply(gg$data[[1]]>qcFilterLevels,gg$data[[2]],mean) %>% unname,0)
    pct.pass=paste0(pct.pass,"%")
    gret = gret + annotate("text",x=seq(pct.pass),y=0,label=pct.pass,hjust="right")

    gret

}

transformVlnPlot0<-function(gg,maxValue=NA) {

    if(is.na(maxValue)) {
        maxValue=quantile(gg$data[[1]],.95)
        cat("maxValue =",maxValue,"\n")
    }

    minValue=0

    newBreaks=scales::breaks_extended(5)(c(minValue,maxValue))

    xLabels=str_wrap(gsub("[-_.]"," ",levels(gg$data[[2]])),15)

    gret=gg +
        geom_jitter(alpha=0.1,size=0.25) +
        scale_y_continuous(limit=c(minValue,maxValue), breaks=newBreaks) +
        scale_x_discrete(labels=xLabels) +
        theme(axis.title.y=element_blank()) +
        coord_flip() +
        NoLegend()

    gret

}
