library(Seurat)
library(fs)
library(ggpubr)

git.describe<-function(REPO="."){

    system2(
        "git",
        c(
            paste0("--git-dir=",file.path(REPO,".git")),
            paste0("--work-tree=",REPO),
            "describe","--tags", "--always", "--long", "--dirty='-UNCOMMITED'"),
        stdout=T
        )

}

if(!exists("glbs")) {
    glbs=list()
}

makeAutoIncrementor <- function(init=0) {
    count <- init
    function() {
        count <<- count + 1
        sprintf("%02d",count)
    }
}

#genomes=c("refdata-gex-mm10-2020-A"="mm10","refdata-gex-GRCh38-2020-A"="hg38")
genomes=unlist(yaml::read_yaml(file.path(SDIR,"data/genomes.yaml"))$genomes)

extractProjNoFromPath<-function(pp) {

    pnum=grep("Project_",strsplit(pp,"/")[[1]],value=T) %>%
        unique %>%
        gsub("Project_","",.) %>% paste0(.,collapse=",")
    if(len(pnum)>0) {
        return(paste0("p",pnum))
    } else {
        return("")
    }

}

get_genome_from_cellranger<-function(cellRangerDir) {

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

    genome

}

read10XDataFolderAsSeuratObj<-function(cellRangerDir,projName) {

    #
    # This version uses the full cellRanger output to infer the Genome
    #

    genome=get_genome_from_cellranger(cellRangerDir)

    dataDir=dir_ls(cellRangerDir,regex="outs/filtered_feature_bc_matrix$",recurs=T)
    if(len(dataDir)!=1) {
        cat("\n\n  Can not find bc_matrix folder\n\n")
        stop("FATAL ERROR")
    }

    xx <- Read10X(dataDir)

    so <- CreateSeuratObject(counts = xx, project=projName)

    cell.barcode=basename(colnames(so)) %>% gsub("filtered_feature_bc_matrix_","",.)
    sampleId=grep("--id=",cmdline,value=T) %>% gsub(".*id=","",.) %>% gsub("^s_","",.)
    pNum=extractProjNoFromPath(grep("--fastq",cmdline,value=T))

    if(pNum!="") {projName=pNum}

    so <- RenameCells(so,new.names=paste0(projName,":s:",sampleId,":bc:",cell.barcode))

    if(genome=="mm10") {
        so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^mt-")
    } else if(genome=="hg38") {
        so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
    } else {
        stop(paste("Unknown genome",genome,"Should not get here"))
    }


    so@meta.data$orig.ident=sampleId
    Idents(so)<-"orig.ident"

    so

}

read_10X_multi_as_SeurateObject<-function(bc_matrix_dir,sampleName,gexSlot,genome,projName) {

    #
    # This require the genome to be passed in
    #
    # It will also handle multi objects if the Gene Expression slot
    # is properly named
    #

    xx <- Read10X(bc_matrix_dir)

    so <- CreateSeuratObject(counts = xx[[gexSlot]], project=projName)

    cell.barcode=basename(colnames(so)) %>% gsub("filtered_feature_bc_matrix_","",.)

    so <- RenameCells(so,new.names=paste0(projName,":bc:",cell.barcode))

    if(genome=="mm10") {
        so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^mt-")
    } else if(genome=="hg38") {
        so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
    } else {
        stop(paste("Unknown genome",genome,"Should not get here"))
    }

    so@meta.data$orig.ident=sampleName

    Idents(so)<-"orig.ident"

    so

}

##############################################################################
# Cell Cycle Functions
#

genes.cellCycle.hg19=cc.genes.updated.2019
genes.cellCycle.hg38=cc.genes.updated.2019
load(file.path(SDIR,"data/cc.genes.mouse.v2.rda"))
genes.cellCycle.mm10=cc.genes.mouse.v2

getCellCycleGenes<-function(genome) {

    if(genome=="mm10") {

        cellCycle.genes=genes.cellCycle.mm10

    } else if(genome=="hg38") {

        cellCycle.genes=genes.cellCycle.hg38

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

    } else if(glbs$genome=="hg38") {

        cellCycle.genes=genes.cellCycle.hg38

    } else {

       stop(paste("seuratTools::scoreCellCycle::Unknown genome",gbls$genome,"Need to implement"))

    }

    cc.features=intersect(unlist(cellCycle.genes),VariableFeatures(sc))

    sc=RunPCA(sc,features=cc.features,approx=FALSE)
    pg=DimPlot(sc,group.by="Phase",reduction="pca") + ggtitle(paste("Cell Cycle PCA Projection",sc@project.name,title))
    pg
}

plotClusterMarkers<-function(cli,s0,clusterColors) {
    # High Quality Genes First
    cli=cli %>% arrange(desc(avg_log2FC+lOR))
    gHQ=cli %>% filter(lOR>2 & avg_log2FC>2) %>% pull(gene)
    gMQ=cli %>% filter(lOR>2) %>% pull(gene)
    gLQ=pull(cli,gene)

    gg=union(gHQ,gMQ) %>% union(.,gLQ) %>% unique(.)
    gg=head(gg,10)
    tbl=cli %>%
        filter(gene %in% gg) %>%
        arrange(factor(gene,levels=gg)) %>%
        select(-lOR) %>%
        mutate_if(is.numeric,~round(.,2)) %>%
        rename(No=1,p.val=3,logFC=4)

    p0=ggtexttable(tbl,rows=NULL)
    vv=VlnPlot(s0,features=gg,pt.size=.025, cols=clusterColors, combine=F)

    pg=paginatePlots(c(list(p0),vv),2,3)

    return(pg)

}
