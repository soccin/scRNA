plot_highest_exprs_genes<-function(obj,n=25,filterRibo=F,filterMito=F) {

    cat("class(obj) =",class(obj),"\n")

    # Compute the relative expression of each gene per cell

    # Use sparse matrix operations, if your dataset is large, doing matrix
    # devisions the regular way will take a very long time.
    #
    # ref: https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/scater/scater_01_qc.html

    if(class(obj)[1]=="Seurat") {

        C=obj@assays$RNA@counts
        C@x=C@x/rep.int(colSums(C), diff(C@p))


    } else if(class(obj)[1]=="dgCMatrix") {
        #
        # Assume above was done of seurat object
        #
        C=obj
    } else {
        cat("\n   Unknown Object ",paste0("[",class(obj),"]"),"\n\n")
        rlang::abort("plot_highest_exprs_genes::Unknown Object")
    }

    titleB=c()

    if(filterRibo) {
        if(glbs$genome=="hg38") {
           C=C[grep("^RP[LSP]",rownames(C),invert=T),]
        } else if(glbs$genome=="mm10") {
            C=C[grep("^Rp[lsp]",rownames(C),invert=T),]
        } else if(glbs$genome=="xenograft") {
           C=C[grep("^RP[LSP]",rownames(C),invert=T),]
        } else {
           rlang::abort(paste("plot_highest_exprs_genes::38::UNKNOWN Genome",glbs$genome))
        }
        titleB=c(titleB,"Filter Ribo")
    }

    if(filterMito) {
        if(glbs$genome=="hg38") {
            C=C[grep("^MT-",rownames(C),invert=T),]
        } else if(glbs$genome=="xenograft") {
            C=C[grep("^MT-",rownames(C),invert=T),]
        } else if(glbs$genome=="mm10") {
            C=C[grep("^mt-",rownames(C),invert=T),]
        } else {
           rlang::abort(paste("plot_highest_exprs_genes::51::UNKNOWN Genome",glbs$genome))
        }
        titleB=c(titleB,"Filter Mito")
    }

    most_expressed=order(Matrix::rowSums(C), decreasing = T)[n:1]
    Ct=t(as.matrix(C[most_expressed,]))

    dc=data.frame(Ct,check.names=F) %>%
        rownames_to_column("CellId") %>%
        gather(Gene,PCT,-CellId) %>%
        tibble %>%
        mutate(Gene=factor(Gene,levels=colnames(Ct)))

    titleA="Percent Total Count per Cell"

    ggplot(dc,aes(Gene,PCT,fill=Gene)) +
        rasterize(geom_boxplot(outlier.shape=95,outlier.size=3,outlier.alpha=.6)) +
        scale_y_continuous(labels=scales::percent) +
        theme_light() +
        NoLegend() +
        theme(axis.text.x = element_text(angle = 45, hjust=1)) +
        xlab("") +
        ylab("% per cell") +
        ggtitle(paste(titleA,titleB,sep=" / "))

}
