# 提取成纤维细胞子集
fibro <- subset(merged_seurat_obj, subset = cell_type == "Fibroblast")
fibro <- NormalizeData(fibro)
fibro <- FindVariableFeatures(fibro, nfeatures = 2000)
hvgs <- VariableFeatures(fibro)
fibro <- ScaleData(fibro, features = hvgs)
fibro <- RunPCA(fibro, features = hvgs, npcs = 20)

# UMAP降维
fibro <- RunUMAP(fibro, dims = 1:15)

# 细胞聚类分析
fibro <- FindNeighbors(fibro, dims = 1:15)
fibro <- FindClusters(fibro, resolution = 0.3)

# 可视化UMAP
DimPlot(fibro, reduction = "umap", label = TRUE) +
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
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8)
  )

# 标记基因鉴定
fibro_markers <- FindAllMarkers(fibro, only.pos = TRUE, min.pct = 0.25, 
                               logfc.threshold = 0.25, test.use = "wilcox")
fibro_significant_markers <- subset(fibro_markers, p_val_adj < 0.05)
write.csv(fibro_significant_markers, "fibro_all_marker.csv")

# 筛选Top标记基因
fibro_significant_markers <- fibro_significant_markers %>% 
  group_by(cluster) %>% 
  top_n(n = 50, wt = avg_log2FC)
write.csv(fibro_significant_markers, "fibro_top_marker.csv")