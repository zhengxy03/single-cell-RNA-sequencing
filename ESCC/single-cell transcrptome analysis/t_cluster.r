t_cells <- subset(merged_seurat_obj, subset = cell_type == "T cell")
t_cells <- NormalizeData(t_cells)
t_cells <- FindVariableFeatures(t_cells, nfeatures = 2000)
hvgs <- VariableFeatures(t_cells)
t_cells <- ScaleData(t_cells, features = hvgs)
t_cells <- RunPCA(t_cells, features = hvgs, npcs = 20)
t_cells <- RunHarmony(t_cells, "orig.ident")
t_cells <- RunUMAP(t_cells, dims = 1:15, reduction = "harmony")
t_cells <- FindNeighbors(t_cells, dims = 1:15, reduction = "harmony")
t_cells <- FindClusters(t_cells, resolution = 0.3)




# 获取图例的个数和名称长度
seurat_clusters <- as.character(unique(t_cells@meta.data$seurat_clusters))  # 转换为字符向量
num_legend_items <- length(seurat_clusters)  # 图例的个数
max_label_length <- max(nchar(seurat_clusters))  # 图例名称的最大长度

# 动态计算图片尺寸
base_width <- 3000  # 基础宽度
base_height <- 3000  # 基础高度
legend_width_factor <- 100  # 每个图例项增加的宽度
label_length_factor <- 10  # 每个字符增加的宽度

# 计算动态宽度
dynamic_width <- base_width + (num_legend_items * legend_width_factor) + (max_label_length * label_length_factor)

png("t_clusters.png", width = dynamic_width, height = base_height, res = 300)
DimPlot(t_cells, reduction = "umap", label = TRUE, pt.size = 2, label.size = 8) +
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    ggtitle(NULL) +
    scale_color_npg() +
    coord_fixed(ratio = 1) +
    guides(color = guide_legend(title = NULL, override.aes = list(size = 5))) +
    theme(
        text = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 28, color = "black"),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title.x = element_text(size = 36, face = "bold", color = "black"),
        axis.title.y = element_text(size = 36, face = "bold", color = "black", margin = margin(r = 20)),  # 增加右侧间距
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 28, face = "bold", color = "black"),
        legend.title = element_text(size = 28, face = "bold", color = "black"),
        legend.position = "right",
        legend.box.margin = margin(0, 0, 0, 0),
        legend.key = element_blank(),
        legend.background = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5),
        aspect.ratio = 1,
        plot.margin = margin(10, 50, 10, 10)
    )
dev.off()



#annotation
t_cell_markers <- FindAllMarkers(t_cells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "wilcox")
t_significant_markers <- subset(t_cell_markers, p_val_adj < 0.05)
write.csv(t_significant_markers, "t_all_marker.csv")
t_significant_markers <- t_significant_markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
write.csv(t_significant_markers, "t_top_marker.csv")

b_cluster_ids <- c(6, 8, 9)
t_cells <- subset(t_cells, idents = setdiff(levels(Idents(t_cells)), b_cluster_ids))

png("t_fature.png", width = dynamic_width, height = base_height, res = 300)
FeaturePlot(t_cells, features = c("CD4", "CD8A", "KLRB1", "FOXP3", 
                  "CRTAM", "IFNG", "FOSB","ITGA1", 
                  "KLRF1", "HAVCR2", "CXCL13")) 
dev.off()

identity_mapping <- c(
  "0" = "KLRB1+ MAIT",
  "1" = "FOXP3+ Treg",
  "2" = "CRTAM+ CD8+Trm", 
  "3" = "IFNG+ CD8+Teff",
  "4" = "Early Activated T",
  "5" = "ITGA1+ CD8+Trm",
  "6" = "TNK",
  "7" = "HAVCR2+ CD8+Tex",
  "8" = "CXCL13+ Tfh"
)
cell_type <- identity_mapping[t_cells@meta.data$seurat_clusters]
t_cells@meta.data$cell_type <- cell_type

# 将 cell_type 列转换为因子，并指定顺序
t_cells@meta.data$cell_type <- factor(t_cells@meta.data$cell_type, levels = identity_mapping)

# 获取唯一的细胞类型并转换为字符向量
cell_types <- as.character(unique(t_cells@meta.data$cell_type))
num_legend_items <- length(cell_types)  # 图例的个数
max_label_length <- max(nchar(cell_types))  # 图例名称的最大长度

# 动态计算图片尺寸
base_width <- 3000  # 基础宽度
base_height <- 3000  # 基础高度
legend_width_factor <- 100  # 每个图例项增加的宽度
label_length_factor <- 10  # 每个字符增加的宽度

# 计算动态宽度
dynamic_width <- base_width + (num_legend_items * legend_width_factor) + (max_label_length * label_length_factor)

# 保存图片
png("t_annotation.png", width = dynamic_width, height = base_height, res = 300)
DimPlot(t_cells, reduction = "umap", label = TRUE, pt.size = 2, label.size = 8, group.by = "cell_type") +
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    ggtitle(NULL) +
    scale_color_npg() +
    coord_fixed(ratio = 1) +
    guides(color = guide_legend(title = NULL, override.aes = list(size = 5))) +
    theme(
        text = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 28, color = "black"),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title.x = element_text(size = 36, face = "bold", color = "black"),
        axis.title.y = element_text(size = 36, face = "bold", color = "black", margin = margin(r = 20)),  # 增加右侧间距
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 28, face = "bold", color = "black"),
        legend.title = element_text(size = 28, face = "bold", color = "black"),
        legend.position = "right",
        legend.box.margin = margin(0, 0, 0, 0),
        legend.key = element_blank(),
        legend.background = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5),
        aspect.ratio = 1,
        plot.margin = margin(10, 50, 10, 10)
    )
dev.off()

CD8 <- subset(t_cells, subset = seurat_clusters %in% c(2, 3, 5, 7))
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