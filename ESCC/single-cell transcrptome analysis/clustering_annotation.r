ElbowPlot(merged_seurat_obj)
merged_seurat_obj <- FindNeighbors(merged_seurat_obj, reduction = "harmony", dims = 1:20)
seq <- seq(0.3, 1.5, by = 0.1)
for (res in seq){
    merged_seurat_obj <- FindClusters(merged_seurat_obj, resolution = res)
}

library(clustree)
library(patchwork)
p3 <- clustree(merged_seurat_obj, prefix = 'RNA_snn_res.') + coord_flip()
p4 <- DimPlot(merged_seurat_obj, group.by = 'RNA_snn_res.0.5', label = T) +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  ggtitle(NULL) +
  guides(color = guide_legend(title = NULL))
#merged_seurat_obj <- FindClusters(merged_seurat_obj, resolution = 0.5)
#merged_seurat_obj <- RunUMAP(merged_seurat_obj, reduction = "harmony", dims = 1:20)
#DimPlot(merged_seurat_obj, reduction = "umap", label = TRUE)

#annotation
markers <- FindAllMarkers(object = merged_seurat_obj, 
                                  test.use = "roc", 
                                  only.pos = TRUE, 
                                  min.pct = 0.25, 
                                  thresh.use = 0.25)

significant_markers <- subset(markers, myAUC > 0.7)
write.csv(significant_markers,"marker.csv")


markers2 <- FindAllMarkers(merged_seurat_obj, 
                          only.pos = TRUE, 
                          min.pct = 0.1, 
                          logfc.threshold = 0.25, 
                          test.use = "wilcox")
#significant markers
significant_markers2 <- subset(markers2, p_val_adj < 0.05)
write.csv(significant_markers2,"marker2.csv")


