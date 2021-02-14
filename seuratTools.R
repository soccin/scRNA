require(fs)

genomes=c("refdata-gex-mm10-2020-A"="mm10")

read10XDataFolderAsSeuratObj<-function(cellRangerDir,projName) {

    cmdlineFile=dir_ls(cellRangerDir,regex="_cmd")
    if(len(cmdlineFile)!=1) {
        cat("\n\n  cellRangerDir =",cellRangerDir,"not a valid CellRanger Directory\n\n")
        stop("FATAL ERROR")
    }
    cmdline=scan(cmdlineFile,"",quiet = TRUE)

    genomeFile=grep("--transcriptome=",cmdline,value=T) %>% basename
    genome=genomes[genomeFile]
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
        attr(so,"genome")<-"mm10"
    } else {
        stop(paste("Unknown genome",genome,"Should not get here"))
    }

    sampleId=grep("--id=",cmdline,value=T) %>% gsub(".*id=","",.) %>% gsub("^s_","",.)

    so@meta.data$orig.ident=sampleId
    Idents(so)<-"orig.ident"

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

    cat("\n\n   Adding Pclaf to cc.genes\n")
    cc.genes=c(cc.genes,"Pclaf")

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

