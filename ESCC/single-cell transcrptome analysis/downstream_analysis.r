#GSE203115
merged_seurat_obj <- NormalizeData(merged_seurat_obj)
merged_seurat_obj <- FindVariableFeatures(merged_seurat_obj, selection.method = "vst", nfeatures = 2000)
merged_seurat_obj <- ScaleData(merged_seurat_obj, features = rownames(merged_seurat_obj))

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

#differentiate non-responder and responder
DimPlot(merged_seurat_obj, 
        reduction = "umap", 
        group.by = "response_status",  
        cols = c("Non-Responder" = "red", "Responder" = "blue"),
        pt.size = 1)

#subsets of non-responder and responder
non_responder_obj <- subset(merged_seurat_obj, subset = response_status == "Non-Responder")
responder_obj <- subset(merged_seurat_obj, subset = response_status == "Responder")

num_clusters <- length(unique(non_responder_obj$seurat_clusters))
grey_colors <- rep("grey", num_clusters)
DimPlot(non_responder_obj, 
        reduction = "umap", 
        group.by = "seurat_clusters",  
        cols = grey_colors,  
        pt.size = 1.5,  
        label = TRUE) +  
  ggtitle("Non-Responder") +
  theme(legend.position = "none")

#DEGs
markers <- FindAllMarkers(merged_seurat_obj, 
                          only.pos = TRUE, 
                          min.pct = 0.1, 
                          logfc.threshold = 0.25, 
                          test.use = "wilcox")
#significant markers
significant_markers <- subset(markers, p_val_adj < 0.05)


write.csv(df_marker,"marker.csv")

#cell annotation
new.cluster.ids <- c("T cell", "Fibroblast", "Myeloid", "Endothelia cell", "Epithelia cell","Epithelia cell", "B cell & mast cell", "SMC & pericyte cell", "plasma", "Double cell")
names(new.cluster.ids) <- levels(merged_seurat_obj)
merged_seurat_obj <- RenameIdents(merged_seurat_obj, new.cluster.ids)
DimPlot(merged_seurat_obj, reduction = "umap", label = TRUE, pt.size = 0.5)

#get T cells subset
t_cells <- subset(merged_seurat_obj, subset = seurat_clusters == "0")
t_cells <- NormalizeData(t_cells)
t_cells <- FindVariableFeatures(t_cells, selection.method = "vst", nfeatures = 2000)
t_cells <- ScaleData(t_cells, features = rownames(t_cells))

#PCA
t_cells <- RunPCA(t_cells, features = VariableFeatures(object = t_cells))

#UMAP
t_cells <- RunUMAP(t_cells, dims = 1:10)

#annotation
t_cells <- FindNeighbors(t_cells, dims = 1:10)
t_cells <- FindClusters(t_cells, resolution = 0.5)

DimPlot(t_cells, reduction = "umap", label = TRUE)

t_cell_markers <- FindAllMarkers(t_cells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "wilcox")
t_significant_markers <- subset(t_cell_markers, p_val_adj < 0.05)
write.csv(t_significant_markers, "t_all_marker.csv")

#t cells annotation
new.cluster.ids <- c("CD8+Tem", "CD8+Tex", "CD4+Tcm", "CD4+Treg", "CD8+Tex", "CD8+Tex", "Cycling T", "Cycling T", "Cycling T", "NK")
names(new.cluster.ids) <- levels(t_cells)
t_cells <- RenameIdents(t_cells, new.cluster.ids)
DimPlot(t_cells, reduction = "umap", label = TRUE, pt.size = 0.5)

#VlnPlot
genes_to_plot <- c("GZMK", "CXCL13", "CTLA4", "FOXP3", "ICOS", "SELL", "HAVCR2", "PRDM1")
VlnPlot(t_cells, features = genes_to_plot, group.by = "seurat_clusters", pt.size = 0)

#bubble heat map
DotPlot(t_cells, features = genes_to_plot, group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Monocle2
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("monocle")
library(monocle)

CD8 <- subset(t_cells, subset = seurat_clusters %in% c(0, 1, 4, 5))
CD4 <- subset(t_cells, subset = seurat_clusters %in% c(2,3))

