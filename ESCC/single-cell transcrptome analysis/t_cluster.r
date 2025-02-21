t_cells <- subset(merged_seurat_obj, subset = cell_type == "T cell")
t_cells <- NormalizeData(t_cells)
t_cells <- FindVariableFeatures(t_cells, nfeatures = 2000)
hvgs <- VariableFeatures(t_cells)
t_cells <- ScaleData(t_cells, features = hvgs)
t_cells <- RunPCA(t_cells, features = hvgs, npcs = 20)
t_cells <- RunHarmony(t_cells, "sample_sources")
t_cells <- RunUMAP(t_cells, dims = 1:15, reduction = "harmony")
t_cells <- FindNeighbors(t_cells, dims = 1:15, reduction = "harmony")
t_cells <- FindClusters(t_cells, resolution = 0.3)

DimPlot(t_cells, reduction = "umap", label = TRUE) +
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
t_cell_markers <- FindAllMarkers(t_cells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "wilcox")
t_significant_markers <- subset(t_cell_markers, p_val_adj < 0.05)
write.csv(t_significant_markers, "t_all_marker.csv")
t_significant_markers <- t_significant_markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
write.csv(t_significant_markers, "t_top_marker.csv")

b_cluster_ids <- c(6, 9, 10)
t_cells <- subset(t_cells, idents = setdiff(levels(Idents(t_cells)), b_cluster_ids))

identity_mapping <- c(
  "0" = "Tn",
  "1" = "CD8+Tem",
  "2" = "Treg", 
  "3" = "CD8+Teff",
  "4" = "CD8+Teff",
  "5" = "CD8+Tex",
  "6" = "TNK",
  "7" = "TNK",
  "8" = "CD4+Tex",
  "9" = "Th17"
)
cell_type <- identity_mapping[t_cells@meta.data$seurat_clusters]
t_cells@meta.data$cell_type <- cell_type

DimPlot(t_cells, reduction = "umap", label = TRUE, group.by = "cell_type") +
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

CD8 <- subset(t_cells, subset = seurat_clusters %in% c(1, 3, 4, 5))
#convert Seurat_obj to SingleCellExperiment_obj
cd8_sce <- as.SingleCellExperiment(CD8)

#Run slingshot
cd8_sce <- slingshot(cd8_sce, clusterLabels = "seurat_clusters", reducedDim = "UMAP")
trajectories <- slingCurves(cd8_sce)
length(trajectories)

#plot
umap_data <- as.data.frame(Embeddings(CD8, reduction = "umap"))
umap_data$cell_type <- CD8$cell_type
curve_data <- as.data.frame(trajectories[[1]]$s[trajectories[[1]]$ord, ])

ggplot(umap_data, aes(x = umap_1, y = umap_2, color = cell_type)) +
  geom_point(size = 2, alpha = 0.8, shape = 16) +
  geom_path(data = curve_data, aes(x = umap_1, y = umap_2), 
            color = "black", size = 1.5, linetype = "solid") +
  scale_color_npg() +  # 使用 Nature Publishing Group 配色
  theme_classic() +
  ggtitle("Trajectory Analysis of CD8 T Cells") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), 
    legend.position = "right",
    legend.title = element_blank()  # 隐藏图例标题
  ) +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  ggtitle(NULL) +
  guides(color = guide_legend(title = NULL))

#all curves
curve_data_list <- lapply(trajectories, function(traj) {
  as.data.frame(traj$s[traj$ord, ])
})
curve_data <- do.call(rbind, curve_data_list)
curve_data$trajectory <- rep(seq_along(trajectories), sapply(curve_data_list, nrow))

ggplot() +
  # 绘制单细胞散点图
  geom_point(data = umap_data, aes(x = umap_1, y = umap_2, color = cell_type), 
             size = 2, alpha = 0.8, shape = 16) +
  # 绘制六条轨迹，指定轨迹颜色为黑色
  geom_path(data = curve_data, aes(x = umap_1, y = umap_2, group = trajectory), 
            color = "black", size = 1.5, linetype = "solid") +
  # 设置细胞类型的颜色比例尺（可根据需要调整）
  scale_color_npg() +
  # 使用经典主题
  theme_classic() +
  # 设置主题相关属性
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), 
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) +
  # 设置坐标轴标签
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  # 设置图形标题
  ggtitle("Trajectory Analysis of CD8 T Cells")