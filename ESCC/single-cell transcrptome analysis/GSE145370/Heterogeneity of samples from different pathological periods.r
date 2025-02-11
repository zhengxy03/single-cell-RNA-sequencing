tumor <- subset(merged_seurat_obj, subset = sample_type == "Tumor")
identity_mapping <- c(
  "Sample1" = "ⅢB",
  "Sample3" = "ⅢA",
  "Sample5" = "ⅡA",
  "Sample7" = "ⅡA",
  "Sample9" = "ⅢB",
  "Sample11" = "ⅡB",
  "Sample13" = "ⅢB"
)
period <- identity_mapping[tumor@meta.data$orig.ident]
tumor@meta.data$period <- period

#proportion
cell_counts2 <- tumor@meta.data %>%
  group_by(sample, ident) %>%
  summarise(count = n()) %>%
  ungroup()

sample_totals <- cell_counts2 %>%
  group_by(sample) %>%
  summarise(total = sum(count)) %>%
  ungroup()

cell_proportion2 <- cell_counts2 %>%
  left_join(sample_totals, by = "sample") %>%
  mutate(proportion = count / total)


ggplot(cell_proportion2, aes(x = ident, y = proportion, fill = sample)) +
  geom_bar(stat = "identity", position = "dodge") +  # 绘制条形图
  labs(x = "Cell Type", y = "Proportion", fill = "Sample Type") +  # 设置坐标轴和图例标题
  theme_classic() +  # 使用经典主题
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # 调整横轴标签角度
    axis.line = element_line(color = "black"),  # 设置坐标轴颜色
    panel.grid.major = element_blank(),  # 移除主要网格线
    panel.grid.minor = element_blank(),  # 移除次要网格线
    legend.position = "right"  # 设置图例位置
  )


#t cells
t_cells <- subset(tumor, subset = ident == "T cell")
t_cells <- NormalizeData(t_cells)
t_cells <- FindVariableFeatures(t_cells, selection.method = "vst", nfeatures = 2000)
t_cells <- ScaleData(t_cells, features = rownames(t_cells))

t_cells <- RunPCA(t_cells, features = VariableFeatures(object = t_cells))
ElbowPlot(t_cells)
t_cells <- RunUMAP(t_cells, dims = 1:15)

t_cells <- FindNeighbors(t_cells, dims = 1:15)
seq <- seq(0.3, 1.5, by = 0.1)
for (res in seq){
    merged_seurat_obj <- FindClusters(merged_seurat_obj, resolution = res)
}
clustree(t_cells, prefix = 'RNA_snn_res.') + coord_flip()
DimPlot(t_cells, group.by = 'RNA_snn_res.0.5', label = T) +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  ggtitle(NULL) +
  guides(color = guide_legend(title = NULL))
#t_cells <- FindClusters(t_cells, resolution = 0.5)
DimPlot(t_cells, reduction = "umap", group.by = "seurat_clusters", label = TRUE, pt.size = 1.5) +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  ggtitle(NULL) +
  guides(color = guide_legend(title = NULL))

t_cell_markers <- FindAllMarkers(t_cells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "wilcox")
t_significant_markers <- subset(t_cell_markers, p_val_adj < 0.05)
write.csv(t_significant_markers, "t_all_marker.csv")

#t cells annotation
identity_mapping <- c("0" = "CD8+Tex", "1" = "CD8+Tem", "2" = "CD8+Tex", "3" = "CD8+Tcm", "4" = "CD4+Treg", "5" = "CD4+Treg", "6" = "CD8+Tex", "7" = "CD8+Trm", "8" = "CD4+Treg", "9" = "CD8+Teff", "10" = "TNK", "11" = "CD8+Teff", "12" = "CD4+Treg", "13" = "CD8+Tcm")

types <- identity_mapping[t_cells@meta.data$seurat_clusters]
t_cells@meta.data$types <- types
DimPlot(t_cells, reduction = "umap", group.by = "types", label = TRUE, pt.size = 1) +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  ggtitle(NULL) +
  guides(color = guide_legend(title = NULL))


#trajectory

library(slingshot)
library(ggsci)

CD8 <- subset(t_cells, subset = seurat_clusters %in% c(0, 1, 2, 3, 6, 7, 9, 11, 13))
#convert Seurat_obj to SingleCellExperiment_obj
cd8_sce <- as.SingleCellExperiment(CD8)

#Run slingshot
cd8_sce <- slingshot(cd8_sce, clusterLabels = "seurat_clusters", reducedDim = "UMAP")
trajectories <- slingCurves(cd8_sce)
length(trajectories)

#plot
umap_data <- as.data.frame(Embeddings(CD8, reduction = "umap"))
umap_data$types <- CD8$types
curve_data <- as.data.frame(trajectories[[1]]$s[trajectories[[1]]$ord, ])

ggplot(umap_data, aes(x = umap_1, y = umap_2, color = types)) +
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
  geom_point(data = umap_data, aes(x = umap_1, y = umap_2, color = types), 
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

#types
cd8_sce <- slingshot(cd8_sce, clusterLabels = "types", reducedDim = "UMAP")
trajectories <- slingCurves(cd8_sce)
length(trajectories)

umap_data <- as.data.frame(Embeddings(CD8, reduction = "umap"))
umap_data$types <- CD8$types
curve_data_list <- lapply(trajectories, function(traj) {
  as.data.frame(traj$s[traj$ord, ])
})
curve_data <- do.call(rbind, curve_data_list)
curve_data$trajectory <- rep(seq_along(trajectories), sapply(curve_data_list, nrow))

ggplot() +
  # 绘制单细胞散点图
  geom_point(data = umap_data, aes(x = umap_1, y = umap_2, color = types), 
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



#monocle3
cds_data <- GetAssayData(CD8, assay = "RNA", layer = "counts")
meta_data <- CD8@meta.data
gene_ann <- data.frame(gene_short_name = rownames(cds_data), row.names = rownames(cds_data))

cds <- new_cell_data_set(
  cds_data,
  cell_metadata = meta_data,
  gene_metadata = gene_ann
)
cds <- preprocess_cds(cds, num_dim = 50)  # 使用前 50 个主成分
cds <- reduce_dimension(cds)

cds <- cluster_cells(cds)
cds <- learn_graph(cds)
plot_cells(cds, color_cells_by = "seurat_clusters", label_groups_by_cluster = FALSE)

#View Pseudo Time Distribution
cds <- order_cells(cds)
plot_cells(cds, color_cells_by = "pseudotime")

pseudotime <- pseudotime(cds)
cell_order <- order(pseudotime)
#extract cluster info
cluster_info <- CD8@meta.data$seurat_clusters
#create dataframe
plot_data <- data.frame(
  pseudotime = pseudotime[cell_order],  # 按伪时间排序
  cluster = cluster_info[cell_order]    # 按伪时间排序
)

ggplot(plot_data, aes(x = pseudotime, fill = cluster)) +
  geom_density(alpha = 0.5, bw = 1) +  # 调整带宽参数
  facet_wrap(~cluster, ncol = 1) +  # 按簇分面显示，每列一个簇
  theme_minimal() +
  theme(
    panel.grid = element_blank(),  # 隐藏网格线
    axis.text.y = element_blank(),  # 隐藏纵坐标刻度
    axis.ticks.y = element_blank(),  # 隐藏纵坐标刻度线
    strip.text = element_blank()  # 隐藏分面标题
  ) +
  labs(
    x = "Pseudotime",  # 横坐标标题
    fill = "Cluster"   # 图例标题
  ) +
  scale_fill_viridis_d()  # 使用 viridis 调色板