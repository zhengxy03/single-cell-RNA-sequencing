ElbowPlot(merged_seurat_obj)
merged_seurat_obj <- FindNeighbors(merged_seurat_obj, reduction = "harmony", dims = 1:20)
#seq <- seq(0.3, 1.5, by = 0.1)
#for (res in seq){
    merged_seurat_obj <- FindClusters(merged_seurat_obj, resolution = res)
}

#library(clustree)
#library(patchwork)
#p3 <- clustree(merged_seurat_obj, prefix = 'RNA_snn_res.') + coord_flip()
#p4 <- DimPlot(merged_seurat_obj, group.by = 'RNA_snn_res.0.5', label = T) +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  ggtitle(NULL) +
  guides(color = guide_legend(title = NULL))

library(ggplot2)
merged_seurat_obj <- FindClusters(merged_seurat_obj, resolution = 0.5)
merged_seurat_obj <- RunUMAP(merged_seurat_obj, reduction = "harmony", dims = 1:20)

DimPlot(merged_seurat_obj, reduction = "umap", label = TRUE, pt.size = 1) +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  ggtitle(NULL) +
  guides(color = guide_legend(title = NULL, override.aes = list(size = 5))) +
  theme(
    text = element_text(size = 8, face = "bold"),
    axis.text.x = element_text(size = 12), 
    axis.text.y = element_text(size = 12), 
    axis.title.x = element_text(size = 12), 
    axis.title.y = element_text(size = 12), 
    plot.title = element_text(size = 12),
    legend.text = element_text(size = 8), # 图例文本字体大小
    legend.title = element_text(size = 8) # 图例标题字体大小（如果有标题的话，这里原代码设置为NULL）
  )
#annotation
#markers <- FindAllMarkers(object = merged_seurat_obj, 
                                  test.use = "roc", 
                                  only.pos = TRUE, 
                                  min.pct = 0.1, 
                                  thresh.use = 0.25)

#significant_markers <- subset(markers, myAUC > 0.7)
#write.csv(significant_markers,"marker.csv")

markers2 <- FindAllMarkers(merged_seurat_obj, 
                          only.pos = TRUE, 
                          min.pct = 0.25, 
                          logfc.threshold = 0.25, 
                          test.use = "wilcox")
#significant markers
significant_markers2 <- subset(markers2, p_val_adj < 0.05)
significant_markers2 <- significant_markers2 %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.csv(significant_markers2,"marker_top.csv")

identity_mapping <- c(
  "0" = "T cell",
  "1" = "B cell",
  "2" = "T cell",
  "3" = "Fibroblast",
  "4" = "Plasma",
  "5" = "Macrophage",
  "6" = "Dendritic cell",
  "7" = "Endothelial cell",
  "8" = "Monocyte",
  "9" = "Epithelial cell",
  "10" = "Fibroblast",
  "11" = "Mast cell",
  "12" = "Epithelial cell",
  "13" = "Proliferating cell",
  "14" = "Pericyte",
  "15" = "T cell",
  "16" = "T cell",
  "17" = "Fibroblast",
  "18" = "T cell",
  "19" = "Plasma",
  "20" = "Plasma",
  "21" = "Proliferating cell",
  "22" = "Epithelial cell"
)
cell_type <- identity_mapping[merged_seurat_obj@meta.data$seurat_clusters]
merged_seurat_obj@meta.data$cell_type <- cell_type

DimPlot(merged_seurat_obj, reduction = "umap", label = TRUE, pt.size = 1, group.by = "cell_type") +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  ggtitle(NULL) +
  guides(color = guide_legend(title = NULL, override.aes = list(size = 5))) +
  theme(
    text = element_text(size = 8, face = "bold"),
    axis.text.x = element_text(size = 12), 
    axis.text.y = element_text(size = 12), 
    axis.title.x = element_text(size = 12), 
    axis.title.y = element_text(size = 12), 
    plot.title = element_text(size = 12),
    legend.text = element_text(size = 8), # 图例文本字体大小
    legend.title = element_text(size = 8) # 图例标题字体大小（如果有标题的话，这里原代码设置为NULL）
  )