read10XIntoSeuratObj <- function(fname,projName) {

    #
    # This is for Human+Mouse data
    #

    xx <- Read10X(fname)
    xx <- xx[grep("^hg19",rownames(xx)),]
    rownames(xx) <- gsub("^hg19_","",rownames(xx))
    so <- CreateSeuratObject(counts = xx, project=projName)
    so <- RenameCells(so,paste0(projName,"_",colnames(x=so)))
    so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")

    so

}

read10XMouseIntoSeuratObj <- function(fname,projName) {

    #
    # This is for Mouse
    #

    xx <- Read10X(fname)
    so <- CreateSeuratObject(counts = xx, project=projName)
    so <- RenameCells(so,paste0(projName,"_",colnames(x=so)))
    so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^mt-")

    so

}

preProcessSO<-function(so) {
    so=NormalizeData(so);
    so=FindVariableFeatures(so);
    so=ScaleData(so,features=rownames(so))
}

plotCellCycle<-function(sc) {
    sc=RunPCA(sc,features=c(cc.genes$s.genes,cc.genes$g2m.genes))
    pg=DimPlot(sc,group.by="Phase") + ggtitle(paste(sc@project.name,"Cell Cycle PCA Projection"))
    pg
}

convertGeneSymbolsHumanToMouse <- function(hgg) {

    require("biomaRt")
    human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

    genesV2 = getLDS(
        attributes = c("hgnc_symbol"),
        filters = "hgnc_symbol",
        values = hgg,
        mart = human,
        attributesL = c("mgi_symbol"),
        martL = mouse,
        uniqueRows=T)

    mgg <- unique(genesV2[, 2])

    mgg

}

scoreCellCycle <- function(dorig) {

    so <- dorig
    so <- NormalizeData(so)
    so <- FindVariableFeatures(so, selection.method="vst")
    so <- ScaleData(so, features=rownames(so))

    cc.genes=lapply(cc.genes.updated.2019,function(x){convertGeneSymbolsHumanToMouse(x)})

    so=CellCycleScoring(so,s.features=cc.genes$s.genes,g2m.features=cc.genes$g2m.genes,set.ident=T)

    cc.meta.data=so@meta.data[,c("S.Score","G2M.Score","Phase")]

    dorig=AddMetaData(dorig,cc.meta.data$Phase,"Phase")
    dorig=AddMetaData(dorig,cc.meta.data$G2M.Score,"G2M.Score")
    dorig=AddMetaData(dorig,cc.meta.data$S.Score,"S.Score")

    Idents(dorig)="Phase"
    return(dorig)

}

regressCellCycle <- function(so,saveVar=T) {

    stop("Need to find MOUSE Cell Cycle Genes")

    checksum=digest::digest(so)
    projName=as.character(so@meta.data$orig.ident[1])
    cacheFile=cc("regressCellCycle",projName,checksum,".rda")

    if(file.exists(cacheFile)) {
        cat("Reading cache",cacheFile,"...")
        xx=readRDS(cacheFile)
        cat("done\n")
        return(xx)
    }

    cat("\n\n=========================================\n\nFile =",cacheFile,"\n")

    all.genes <- rownames(so)

    so <- ScaleData(so, features = all.genes)

    #
    ## Cell Cycle Analysis
    #

    cc.genes=cc.genes.updated.2019

    so=CellCycleScoring(so,s.features=cc.genes$s.genes,g2m.features=cc.genes$g2m.genes,set.ident=T)

    if("Phase" %in% colnames(so@meta.data)) {

        RidgePlot(so, features=c("PCNA","CDC20","AURKA"),ncol=2)
        so=RunPCA(so, features=c(cc.genes$s.genes,cc.genes$g2m.genes))
        DimPlot(so) + ggtitle("Pre Cell Cycle Regression")

        so <- ScaleData(so, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(so))
        so <- RunPCA(so, features=c(cc.genes$s.genes,cc.genes$g2m.genes))
        DimPlot(so) + ggtitle("Post Cell Cycle Regression")


    } else {
        cat("\n\nCan not fit/score cell cycle phase\n\n")
    }

    if(saveVar) {
        saveRDS(so,cacheFile,compress=T)
    }
    return(so)

}

