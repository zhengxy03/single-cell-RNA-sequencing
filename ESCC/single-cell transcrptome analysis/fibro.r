fibroblasts <- subset(merged_seurat_obj, subset = cell_type == "Fibroblast")
fibroblasts <- NormalizeData(fibroblasts)
fibroblasts <- FindVariableFeatures(fibroblasts, nfeatures = 2000)
hvgs <- VariableFeatures(fibroblasts)
fibroblasts <- ScaleData(fibroblasts, features = hvgs)
fibroblasts <- RunPCA(fibroblasts, features = hvgs, npcs = 20)

fibroblasts <- RunHarmony(fibroblasts, "sample_sources")
fibroblasts <- RunUMAP(fibroblasts, dims = 1:15, reduction = "harmony")
fibroblasts <- FindNeighbors(fibroblasts, dims = 1:15, reduction = "harmony")
fibroblasts <- FindClusters(fibroblasts, resolution = 0.3)

DimPlot(fibroblasts, reduction = "umap", label = TRUE) +
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

fibroblast_markers <- FindAllMarkers(fibroblasts, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "wilcox")
fibroblast_significant_markers <- subset(fibroblast_markers, p_val_adj < 0.05)
#write.csv(fibroblast_significant_markers, "fibroblast_all_marker.csv")
fibroblast_significant_markers <- fibroblast_significant_markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
write.csv(fibroblast_significant_markers, "fibroblast_top_marker.csv")

neg_fibroblast_significant_markers <- fibroblast_significant_markers %>% 
  group_by(cluster) %>% 
  top_n(n = -50, wt = avg_log2FC)
write.csv(neg_fibroblast_significant_markers, "neg_fibroblast_top_marker.csv")

immune_cluster_ids <- c(9)
fibroblasts <- subset(fibroblasts, idents = setdiff(levels(Idents(fibroblasts)), immune_cluster_ids))

FeaturePlot(merged_seurat_obj, features = "CXCL5")
fibro_T <- subset(fibroblasts, subset = sample_type == "tumor")
fibro_N <- subset(fibroblasts, subset = sample_type == "normal")

identity_mapping <- c(
    "0" = "CAF1",
    "1" = "CAF2",
    "2" = "Pericyte",
    "3" = "Pericyte",
    "4" = "CAF3",
    "5" = "CAF4",
    "6" = "Pericyte",
    "7" = "NAF1",
    "8" = "myCAF5",
    "9" = "apCAF6",
    "10" = "mCAF7",
    "11" = "NAF2",
    "12" = "NAF3",
    "13" = "NAF4"
)

cell_type <- identity_mapping[fibroblasts@meta.data$seurat_clusters]
fibroblasts@meta.data$cell_type <- cell_type

DimPlot(fibroblasts, reduction = "umap", label = TRUE, group.by = "cell_type") +
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

#trajectory
fibroblast_sce <- as.SingleCellExperiment(fibroblasts)

# 运行 slingshot 进行轨迹分析
fibroblast_sce <- slingshot(fibroblast_sce, clusterLabels = "seurat_clusters", reducedDim = "UMAP")
trajectories <- slingCurves(fibroblast_sce)
length(trajectories)

# 提取 UMAP 数据
umap_data <- as.data.frame(Embeddings(fibroblasts, reduction = "umap"))
umap_data$types <- fibroblasts$cell_type

# 提取轨迹曲线数据
curve_data <- as.data.frame(trajectories[[1]]$s[trajectories[[1]]$ord, ])

# 绘制轨迹分析图
ggplot(umap_data, aes(x = umap_1, y = umap_2, color = cell_type)) +
  geom_point(size = 2, alpha = 0.8, shape = 16) +
  geom_path(data = curve_data, aes(x = umap_1, y = umap_2), 
            color = "black", size = 1.5, linetype = "solid") +
  scale_color_npg() +  # 使用 Nature Publishing Group 配色
  theme_classic() +
  ggtitle("Trajectory Analysis of Fibroblast Cells") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), 
    legend.position = "right",
    legend.title = element_blank()  # 隐藏图例标题
  ) +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  guides(color = guide_legend(title = NULL))


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


ggplot(umap_data, aes(x = umap_1, y = umap_2, color = cell_type)) +
  geom_point(size = 2, alpha = 0.8, shape = 16) +
  geom_path(data = curve_data, aes(x = umap_1, y = umap_2), 
            color = "black", linewidth = 1.5, linetype = "solid") +
  scale_color_brewer(palette = "Set3") +  # 使用能处理更多颜色的调色板
  theme_classic() +
  ggtitle("Trajectory Analysis of Fibroblast Cells") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), 
    legend.position = "right",
    legend.title = element_blank()  # 隐藏图例标题
  ) +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  guides(color = guide_legend(title = NULL))