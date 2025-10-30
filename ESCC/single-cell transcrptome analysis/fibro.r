fibroblasts <- subset(merged_seurat_obj, subset = cell_type == "Fibroblast")
fibroblasts <- NormalizeData(fibroblasts)
fibroblasts <- FindVariableFeatures(fibroblasts, nfeatures = 2000)
hvgs <- VariableFeatures(fibroblasts)
fibroblasts <- ScaleData(fibroblasts, features = hvgs)
fibroblasts <- RunPCA(fibroblasts, features = hvgs, npcs = 20)
library(harmony)
fibroblasts <- RunHarmony(fibroblasts, "sample_sources")
fibroblasts <- RunUMAP(fibroblasts, dims = 1:20, reduction = "harmony")
fibroblasts <- FindNeighbors(fibroblasts, dims = 1:20, reduction = "harmony")
fibroblasts <- FindClusters(fibroblasts, resolution = 0.3)

# 获取图例的个数和名称长度
seurat_clusters <- as.character(unique(fibroblasts@meta.data$seurat_clusters))  # 转换为字符向量
num_legend_items <- length(seurat_clusters)  # 图例的个数
max_label_length <- max(nchar(seurat_clusters))  # 图例名称的最大长度

# 动态计算图片尺寸
base_width <- 3000  # 基础宽度
base_height <- 3000  # 基础高度
legend_width_factor <- 100  # 每个图例项增加的宽度
label_length_factor <- 10  # 每个字符增加的宽度

# 计算动态宽度
dynamic_width <- base_width + (num_legend_items * legend_width_factor) + (max_label_length * label_length_factor)

npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(15)
png("fibro_clusters.png", width = dynamic_width, height = base_height, res = 300)
DimPlot(fibroblasts, reduction = "umap", label = TRUE, pt.size = 2, label.size = 8) +
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    ggtitle(NULL) +
    scale_color_manual(values = npg_extended) +
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


fibroblast_markers <- FindAllMarkers(fibroblasts, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "wilcox")
fibroblast_significant_markers <- subset(fibroblast_markers, p_val_adj < 0.05)
write.csv(fibroblast_significant_markers, "fibroblast_all_marker.csv")
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


# 获取图例的个数和名称长度
seurat_clusters <- as.character(unique(fibro_T@meta.data$seurat_clusters))  # 转换为字符向量
num_legend_items <- length(seurat_clusters)  # 图例的个数
max_label_length <- max(nchar(seurat_clusters))  # 图例名称的最大长度

# 动态计算图片尺寸
base_width <- 3000  # 基础宽度
base_height <- 3000  # 基础高度
legend_width_factor <- 100  # 每个图例项增加的宽度
label_length_factor <- 10  # 每个字符增加的宽度

# 计算动态宽度
dynamic_width <- base_width + (num_legend_items * legend_width_factor) + (max_label_length * label_length_factor)

npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(15)
png("fibro_T_clusters.png", width = dynamic_width, height = base_height, res = 300)
DimPlot(fibro_T, reduction = "umap", label = TRUE, pt.size = 2, label.size = 8) +
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    ggtitle(NULL) +
    scale_color_manual(values = npg_extended) +
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





unique_periods <- sort(unique(fibroblasts@meta.data$period))  # 按字母顺序排序
# 如果需要手动指定顺序，可以这样做：
# unique_periods <- c("period_A", "period_B", "period_C", "period_D", "period_E", "period_F")

# 创建一个空列表，用于存储每个 period 的图
plot_list <- list()

# 遍历每个 period，生成单独的图
for (period in unique_periods) {
    # 创建颜色映射：当前 period 为 npg 红色，其他为灰色
    color_mapping <- setNames(
        ifelse(unique_periods == period, pal_npg()(1), "gray"),  # npg 红色
        unique_periods
    )
    
    # 绘制 DimPlot
    p <- DimPlot(fibroblasts, reduction = "umap", label = FALSE, pt.size = 1, group.by = "period") +
        scale_color_manual(values = color_mapping) +  # 使用自定义颜色映射
        ggtitle(period) +  # 设置标题为当前 period
        theme(
            legend.position = "none",  # 隐藏图例
            plot.title = element_text(size = 24, face = "bold", hjust = 0.5),  # 标题样式
            axis.title.x = element_text(size = 20, face = "bold", color = "black"),  # X 轴标题
            axis.title.y = element_text(size = 20, face = "bold", color = "black"),  # Y 轴标题
            axis.text.x = element_text(size = 16, color = "black"),  # X 轴刻度
            axis.text.y = element_text(size = 16, color = "black")   # Y 轴刻度
        )
    
    # 将图添加到列表中
    plot_list[[period]] <- p
}

# 使用 patchwork 将图形排列在一起（每行三个）
combined_plot <- wrap_plots(plot_list, ncol = 4)  # 每行三个图
png("fibro_period_umap.png", width = 6000, height =3000, res = 300)
print(combined_plot)
dev.off()



library(ggrepel)
identity_mapping <- c(
    "0" = "ACTG2+ myCAFs",
    "1" = "CD34+ Fib progenitors",
    "2" = "SLPI+ NMFs",
    "3" = "IGF1+ inflammatory NAFs",
    "4" = "MMP1+ myCAFs",
    "5" = "CXCL1+ iCAFs",
    "6" = "CCL20+ iCAFs",
    "7" = "COL6A5+ myCAFs",
    "8" = "vCAFs",
    "9" = "COL27A1+ myCAFs",
    "10" = "TIMP1+ myCAFs",
    "11" = "Immune-associated NAFs",
    "12" = "Proliferative CAFs",
    "13" = "EndMT CAFs"
)

# 假设 identity_mapping 已经定义
cell_type <- identity_mapping[fibroblasts@meta.data$seurat_clusters]
fibroblasts@meta.data$cell_type <- cell_type

# 将 cell_type 列转换为因子，并指定顺序
fibroblasts@meta.data$cell_type <- factor(fibroblasts@meta.data$cell_type, levels = identity_mapping)

# 获取唯一的细胞类型并转换为字符向量
cell_types <- as.character(unique(fibroblasts@meta.data$cell_type))
num_legend_items <- length(cell_types)  # 图例的个数
max_label_length <- max(nchar(cell_types))  # 图例名称的最大长度

# 动态计算图片尺寸
base_width <- 5000  # 基础宽度
base_height <- 5000  # 基础高度
legend_width_factor <- 100  # 每个图例项增加的宽度
label_length_factor <- 10  # 每个字符增加的宽度

# 计算动态宽度
dynamic_width <- base_width + (num_legend_items * legend_width_factor) + (max_label_length * label_length_factor)

# 保存图片

umap_data <- as.data.frame(fibroblasts@reductions$umap@cell.embeddings)
umap_data$cell_type <- fibroblasts@meta.data$cell_type

# 检查列名
print(head(umap_data))

# 确保列名正确
colnames(umap_data) <- c("umap_1", "umap_2", "cell_type")

# 计算每个细胞类型的中心点
centroids <- umap_data %>%
    group_by(cell_type) %>%
    summarise(
        umap_1 = median(umap_1),
        umap_2 = median(umap_2)
    )

# 绘制 UMAP 图
p <- DimPlot(fibroblasts, reduction = "umap", label = FALSE, pt.size = 2, group.by = "cell_type") +
    geom_text_repel(
        data = centroids,  # 使用中心点数据
        aes(x = umap_1, y = umap_2, label = cell_type),  # 指定 x 和 y 的美学映射
        size = 8,  # 标签字体大小
        box.padding = 0.5,  # 标签与点之间的间距
        point.padding = 0.5,  # 标签之间的间距
        max.overlaps = Inf,  # 允许的最大重叠次数
        force = 1,  # 调整标签的排斥力
        min.segment.length = 0  # 强制显示所有标签的连接线
    ) +
    xlab("UMAP_1") +  # 添加 x 轴标签
    ylab("UMAP_2") +  # 添加 y 轴标签
    ggtitle(NULL) +  # 移除标题
    scale_color_manual(values = npg_extended) +  # 使用自定义颜色
    coord_fixed(ratio = 1) +  # 固定坐标轴比例
    guides(color = guide_legend(title = NULL, override.aes = list(size = 5))) +  # 调整图例
    theme(
        text = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 40, color = "black"),
        axis.text.y = element_text(size = 40, color = "black"),
        axis.title.x = element_text(size = 56, face = "bold", color = "black"),
        axis.title.y = element_text(size = 56, face = "bold", color = "black", margin = margin(r = 20)),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 40, face = "bold", color = "black"),
        legend.title = element_text(size = 40, face = "bold", color = "black"),
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
    ) + scale_x_continuous(limits = c(-10, 10.5))

# 保存图片
png("fibroblasts_annotation.png", width = dynamic_width, height = base_height, res = 300)
print(p)
dev.off()






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
  scale_color_brewer(palette = "Set3") +
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


