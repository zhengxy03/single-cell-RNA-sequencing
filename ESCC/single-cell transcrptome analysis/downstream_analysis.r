#GSE203115
merged_seurat_obj <- NormalizeData(merged_seurat_obj)
merged_seurat_obj <- FindVariableFeatures(merged_seurat_obj, selection.method = "vst", nfeatures = 2000)
merged_seurat_obj <- ScaleData(merged_seurat_obj, features = rownames(merged_seurat_obj))
merged_seurat_obj$response_status <- ifelse(
  merged_seurat_obj$orig.ident %in% c("Sample3"), "Non-Responder", "Responder"
)


#PCA
merged_seurat_obj <- RunPCA(merged_seurat_obj, features = VariableFeatures(object = merged_seurat_obj))

#tSNE
merged_seurat_obj <- RunTSNE(merged_seurat_obj, dims = 1:10)

#UMAP
merged_seurat_obj <- RunUMAP(merged_seurat_obj, dims = 1:10)

#unsupervised clusting
merged_seurat_obj <- FindNeighbors(merged_seurat_obj, dims = 1:10)
merged_seurat_obj <- FindClusters(merged_seurat_obj, resolution = 0.2)

DimPlot(merged_seurat_obj, reduction = "umap", label = TRUE)