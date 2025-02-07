merged_seurat_obj <- FindNeighbors(merged_seurat_obj, reduction = "harmony", dims = 1:20)
merged_seurat_obj <- FindClusters(merged_seurat_obj, resolution = 1.5)
merged_seurat_obj <- RunUMAP(merged_seurat_obj, reduction = "harmony", dims = 1:20)
DimPlot(merged_seurat_obj, reduction = "umap", label = TRUE)

cluster_markers <- FindAllMarkers(object = merged_seurat_obj, 
                                  test.use = "roc", 
                                  only.pos = TRUE, 
                                  min.pct = 0.25, 
                                  thresh.use = 0.25)

#plot according to sample types
DimPlot(merged_seurat_obj, reduction = "umap", group.by = "sample_type", pt.size = 0.5)

p1 <- DimPlot(merged_seurat_obj, reduction = "umap", group.by = "sample_type", pt.size = 0.5) +
  # 设置横纵坐标标签
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  # 移除标题
  ggtitle(NULL) +
  # 移除图例标题（如果需要）
  guides(color = guide_legend(title = NULL))

print(p1)