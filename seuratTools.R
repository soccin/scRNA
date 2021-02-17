library(Seurat)
library(fs)

if(!exists("glbs")) {
    glbs=list()
}

genomes=c("refdata-gex-mm10-2020-A"="mm10")

genes.cellCycle.hg19=cc.genes.updated.2019
delayedAssign("genes.cellCycle.mm10",loadCellCycleGenes())

makeAutoIncrementor <- function(init=0) {
    count <- init
    function() {
        count <<- count + 1
        sprintf("%02d",count)
    }
}


loadCellCycleGenes <- function() {

    lapply(cc.genes.updated.2019,function(x){convertGeneSymbolsHumanToMouse(x)})

}

library(memoise)
m.cache=cache_filesystem("__R_CACHE")

.convertGeneSymbolsHumanToMouse.RAW <- function(hgg) {

    require("biomaRt")
    human = useMart("ensembl", dataset = "hsapiens_gene_ensembl", host="useast.ensembl.org")
    mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl", host="useast.ensembl.org")

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

convertGeneSymbolsHumanToMouse <- memoise(.convertGeneSymbolsHumanToMouse.RAW,cache=m.cache)

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

getCellCycleGenes<-function(genome) {

    if(genome=="mm10") {

        cellCycle.genes=genes.cellCycle.mm10

    } else {

       stop(paste("seuratTools::scoreCellCycle::Unknown genome",gbls$genome,"Need to implement"))

    }

    cellCycle.genes

}


scoreCellCycle <- function(dorig) {

    so <- dorig
    so <- NormalizeData(so)
    so <- FindVariableFeatures(so, selection.method="vst")
    so <- ScaleData(so, features=rownames(so))

    cellCycle.genes = getCellCycleGenes(glbs$genome)

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

plotCellCycle<-function(sc,title="") {

    if(glbs$genome=="mm10") {

        cellCycle.genes=genes.cellCycle.mm10

    } else {

       stop(paste("seuratTools::scoreCellCycle::Unknown genome",gbls$genome,"Need to implement"))

    }

    sc=RunPCA(sc,features=c(cellCycle.genes$s.genes,cellCycle.genes$g2m.genes))
    pg=DimPlot(sc,group.by="Phase") + ggtitle(paste("Cell Cycle PCA Projection",sc@project.name,title))
    pg
}

