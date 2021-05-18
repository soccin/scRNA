regressCellCycle <- function(so,saveVar=T) {

    checksum=digest::digest(so)
    projName=so@project.name
    cacheFile=cc("regressCellCycle",projName,checksum,".rda")

    if(file.exists(cacheFile)) {
        cat("Reading cache",cacheFile,"...")
        xx=readRDS(cacheFile)
        cat("done\n")
        return(xx)
    }

    cat("\n\n=========================================\n\nFile =",cacheFile,"\n")

    all.genes <- rownames(so)

    so <- NormalizeData(so)
    so <- FindVariableFeatures(so, selection.method="vst")
    so <- ScaleData(so, features = all.genes)

    #
    ## Cell Cycle Analysis
    #

    if(glbs$genome=="mm10") {

        stop("\n\nIs this code still used?]\nProbably not a good idea-custom genomes\n\n")

        cellCycle.genes=genes.cellCycle.mm10
        rPlot.features=c("Pcna","Cdc20","Aurka")

    } else {

       stop(paste("seuratTools::scoreCellCycle::Unknown genome",gbls$genome,"Need to implement"))

    }

    cc.features=c(cellCycle.genes$s.genes,cellCycle.genes$g2m.genes)

    so=CellCycleScoring(so,
                    s.features=cellCycle.genes$s.genes,
                    g2m.features=cellCycle.genes$g2m.genes,
                    set.ident=T)

    cc.plts=list()

    if("Phase" %in% colnames(so@meta.data)) {

        cc.plts[[1]]=RidgePlot(so, features=rPlot.features,ncol=2)
        so=RunPCA(so, features=cc.features)
        cc.plts[[2]]=DimPlot(so) + ggtitle("Pre Cell Cycle Regression")

        so <- ScaleData(so, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(so))
        so <- RunPCA(so, features=cc.features)
        cc.plts[[3]]=DimPlot(so) + ggtitle("Post Cell Cycle Regression")


    } else {
        cat("\n\nCan not fit/score cell cycle phase\n\n")
    }

    if(saveVar) {
        saveRDS(list(so=so,cc.plts=cc.plts),cacheFile,compress=T)
    }

    return(list(so=so,cc.plts=cc.plts))

}

