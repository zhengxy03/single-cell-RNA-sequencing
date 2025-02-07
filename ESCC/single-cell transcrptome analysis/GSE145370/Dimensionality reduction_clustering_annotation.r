merged_seurat_obj <- FindNeighbors(merged_seurat_obj, reduction = "harmony", dims = 1:20)
merged_seurat_obj <- FindClusters(merged_seurat_obj, resolution = 1.5)
seq <- seq(0.3, 1.5, by = 0.1)
for (res in seq){
    merged_seurat_obj <- FindClusters(merged_seurat_obj, resolution = res)
}
merged_seurat_obj %>% 
  SetIdent(value = 'RNA_snn_res.0.5') %>%
  FindAllMarkers(test.use = 'roc') %>%
  filter(myAUC > 0.6) %>%
  count(cluster, name = 'number')

merged_seurat_obj <- RunUMAP(merged_seurat_obj, reduction = "harmony", dims = 1:20)
DimPlot(merged_seurat_obj, reduction = "umap", label = TRUE)

cluster_markers <- FindAllMarkers(object = merged_seurat_obj, 
                                  test.use = "roc", 
                                  only.pos = TRUE, 
                                  min.pct = 0.25, 
                                  thresh.use = 0.25)

#plot according to sample types
p1 <- DimPlot(merged_seurat_obj, reduction = "umap", group.by = "sample_type", pt.size = 0.5) +
  # 设置横纵坐标标签
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  # 移除标题
  ggtitle(NULL) +
  # 移除图例标题（如果需要）
  guides(color = guide_legend(title = NULL))

print(p1)

#plot according to sample sources
p2 <- DimPlot(merged_seurat_obj, reduction = "umap", group.by = "sample_sources", pt.size = 0.5) +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  ggtitle(NULL) +
  guides(color = guide_legend(title = NULL))

print(p2)

