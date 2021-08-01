
cat("Merging with new workflow based on SCTransform")


d10X = lapply(d10X,FUN=SCTransform)
features <- SelectIntegrationFeatures(object.list = d10X, nfeatures = 3000)
d10X <- PrepSCTIntegration(object.list = d10X, anchor.features = features)

anchors <- FindIntegrationAnchors(
        object.list = d10X,
        normalization.method = "SCT",
        anchor.features = features
        )

d10X.integrate <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

d10X.integrate <- RunPCA(d10X.integrate, verbose = FALSE)
d10X.integrate <- RunUMAP(d10X.integrate, reduction = "pca", dims = 1:30)

DimPlot(d10X.integrate,reduction="umap")


