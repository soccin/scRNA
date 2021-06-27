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
