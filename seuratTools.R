library(Seurat)
library(fs)

glbs=list()

genomes=c("refdata-gex-mm10-2020-A"="mm10")

genes.cellCycle.hg19=cc.genes.updated.2019
delayedAssign("genes.cellCycle.mm10",loadCellCycleGenes())

loadCellCycleGenes <- function() {

    lapply(cc.genes.updated.2019,function(x){convertGeneSymbolsHumanToMouse(x)})

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

read10XDataFolderAsSeuratObj<-function(cellRangerDir,projName) {

    cmdlineFile=dir_ls(cellRangerDir,regex="_cmd")
    if(len(cmdlineFile)!=1) {
        cat("\n\n  cellRangerDir =",cellRangerDir,"not a valid CellRanger Directory\n\n")
        stop("FATAL ERROR")
    }
    cmdline=scan(cmdlineFile,"",quiet = TRUE)

    genomeFile=grep("--transcriptome=",cmdline,value=T) %>% basename
    genome=genomes[genomeFile]
    glbs$genome <<- union(glbs$genome,genome)
    if(is.na(genome)) {
        cat("\n\n  Unknown genomeFile =",genomeFile,"\n\n")
        stop("FATAL ERROR")
    }

    dataDir=dir_ls(cellRangerDir,regex="outs/filtered_feature_bc_matrix$",recurs=T)
    if(len(dataDir)!=1) {
        cat("\n\n  Can not find bc_matrix folder\n\n")
        stop("FATAL ERROR")
    }

    xx <- Read10X(dataDir)

    so <- CreateSeuratObject(counts = xx, project=projName,)
    so <- RenameCells(so,paste0(projName,"_",colnames(x=so)))
    if(genome=="mm10") {
        so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^mt-")
    } else {
        stop(paste("Unknown genome",genome,"Should not get here"))
    }

    sampleId=grep("--id=",cmdline,value=T) %>% gsub(".*id=","",.) %>% gsub("^s_","",.)

    so@meta.data$orig.ident=sampleId
    Idents(so)<-"orig.ident"

    so

}


scoreCellCycle <- function(dorig) {

    so <- dorig
    so <- NormalizeData(so)
    so <- FindVariableFeatures(so, selection.method="vst")
    so <- ScaleData(so, features=rownames(so))

    if(glbs$genome=="mm10") {

        cellCycle.genes=genes.cellCycle.mm10

    } else {

       stop(paste("seuratTools::scoreCellCycle::Unknown genome",gbls$genome,"Need to implement"))

    }

    so=CellCycleScoring(so,
                        s.features=cellCycle.genes$s.genes,
                        g2m.features=cellCycle.genes$g2m.genes,
                        set.ident=T
                        )

    cc.meta.data=so@meta.data[,c("S.Score","G2M.Score","Phase")]

    dorig=AddMetaData(dorig,cc.meta.data$Phase,"Phase")
    dorig=AddMetaData(dorig,cc.meta.data$G2M.Score,"G2M.Score")
    dorig=AddMetaData(dorig,cc.meta.data$S.Score,"S.Score")

    Idents(dorig)="Phase"
    return(dorig)

}

preProcessSO<-function(so) {
    so=NormalizeData(so);
    so=FindVariableFeatures(so);
    so=ScaleData(so,features=rownames(so))
}

plotCellCycle<-function(sc) {

    if(glbs$genome=="mm10") {

        cellCycle.genes=genes.cellCycle.mm10

    } else {

       stop(paste("seuratTools::scoreCellCycle::Unknown genome",gbls$genome,"Need to implement"))

    }


    sc=RunPCA(sc,features=c(cellCycle.genes$s.genes,cellCycle.genes$g2m.genes))
    pg=DimPlot(sc,group.by="Phase") + ggtitle(paste(sc@project.name,unname(sc$orig.ident[1]),"Cell Cycle PCA Projection"))
    pg
}



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

    so <- ScaleData(so, features = all.genes)

    #
    ## Cell Cycle Analysis
    #

    if(glbs$genome=="mm10") {

        cellCycle.genes=genes.cellCycle.mm10
        rPlot.features=c("Pcna","Cdc20","Aurka")

    } else {

       stop(paste("seuratTools::scoreCellCycle::Unknown genome",gbls$genome,"Need to implement"))

    }

    cc.features=c(cellCycle.genes$s.genes,cellCycle.genes$g2m.genes)

    so=CellCycleScoring(so,s.features=cellCycle.genes$s.genes,g2m.features=cellCycle.genes$g2m.genes,set.ident=T)

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
        saveRDS(so,cacheFile,compress=T)
    }

    return(list(so=so,cc.plts=cc.plts))

}

