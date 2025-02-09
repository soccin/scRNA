cat("\nPlot modules ...")
pm=list()
pn=list()
for(ii in seq(len(modules))) {
    print(ii)
    modTag=paste0("Modules",ii)
    pp=FeaturePlot(s1,
        features=modTag,
        max.cutoff="q95",min.cutoff="q05",
        combine=F)

    pm[[ii]]=pp[[1]] + ggtitle(names(modules)[ii])

    pn[[ii]] = ggplot(s1@meta.data,aes_string(clusterRes,modTag,fill=clusterRes)) +
        geom_violin() +
        theme_light() +
        geom_jitter(alpha=.1,size=.7,width=.2) +
        ylab("Module Score") +
        ggtitle(names(modules)[ii]) +
        theme(legend.position = "none")


}

cat(" done\n\n")
pfile=get_plot_filename(plotNo(),"ModuleScores_%03d.png")
pngCairo(pfile,width=11,height=8.5)
print(paginatePlots(pm,2,2,FALSE))
dev.off()
mergePNGs(pfile)

pfile=get_plot_filename(plotNo(),"ModuleDistribution_%03d.png")
pngCairo(pfile,width=8.5,height=11)
print(paginatePlots(pn,3,1,FALSE))
dev.off()
mergePNGs(pfile)

#
# Dump metadata
#

md=s1@meta.data %>% data.frame %>% rownames_to_column("CellID") %>% tibble
moduleCols=grep("^Modules\\d+",names(md))
colnames(md)[moduleCols]=paste0("mod.",names(modules))
write_csv(md,"metaData_AddModules.csv")










