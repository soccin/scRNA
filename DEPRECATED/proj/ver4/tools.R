paginatePlots<-function(plts,pRows,pCols) {

    nPlots=pRows*pCols

    pp=list()
    page=1
    currPlot=1
    pa=NULL

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
