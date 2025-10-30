fibroblasts <- subset(merged_seurat_obj, subset = cell_type == "Fibroblast")
library(Seurat)
fibroblasts <- subset(fibroblasts,subset = seurat_clusters %in% c(0,1,2,3,4,5,6,7,8,9,11,12))
fibroblasts <- NormalizeData(fibroblasts)
fibroblasts <- FindVariableFeatures(fibroblasts, nfeatures = 2000)
hvgs <- VariableFeatures(fibroblasts)
fibroblasts <- ScaleData(fibroblasts, features = hvgs)
fibroblasts <- RunPCA(fibroblasts, features = hvgs, npcs = 20)
library(harmony)
fibroblasts <- RunHarmony(fibroblasts, "sample_sources")
fibroblasts <- RunUMAP(fibroblasts, dims = 1:20, reduction = "harmony")
fibroblasts <- FindNeighbors(fibroblasts, dims = 1:20, reduction = "harmony")
fibroblasts <- FindClusters(fibroblasts, resolution = 0.4)
saveRDS(fibroblasts,file="fibro_remove.rds")
library(ggplot2)
library(dplyr)
library(ggsci)
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
pdf("fibro_clusters.pdf", width = dynamic_width/300, height = base_height/300)
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

fibroblasts <- subset(fibroblasts, subset = seurat_clusters %in% c(0,1,2,3,4,5,6,7,8,9,10,11))

library(ggrepel)
identity_mapping <- c(
    "0" = "ACTG2+ myCAF",
    "1" = "CD14+ iCAF",
    "2" = "CXCL12+ perivascular Fib",#血管周成纤维
    "3" = "MMP1+ myCAFs",
    "4" = "CCL20+ iCAFs",
    "5" = "IGF1+ iNAF",
    "6" = "RGS+ perivascular Fib",
    "7" = "CLEC14A+ perivascular Fib",
    "8" = "CXCL1+ iCAF",
    "9" = "VEGFD+ perivascular Fib",
    "10" ="MYL9+ myCAF",
    "11" = "CCL13+ iCAF"
)

npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(12)
# 假设 identity_mapping 已经定义
cell_type <- identity_mapping[fibroblasts@meta.data$seurat_clusters]
fibroblasts@meta.data$cell_type <- cell_type

# 将 cell_type 列转换为因子，并指定顺序
#fibroblasts@meta.data$cell_type <- factor(fibroblasts@meta.data$cell_type, levels = identity_mapping)

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
pdf("fibroblasts_annotation.pdf", width = dynamic_width/300, height = base_height/300)
print(p)
dev.off()


#不同组织中的分布
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(12)

fibroblasts$sample_type <- factor(fibroblasts$sample_type, 
                              levels = c("Normal", "Tumor", "mLN","nLN"))

pdf("fibroblasts_cluster_by_sampletype.pdf", width = 9000/300, height = 3000/300)

p <- DimPlot(fibroblasts, 
             reduction = "umap", 
             label = FALSE, 
             pt.size = 1, 
             group.by = "cell_type", 
             label.size = 8,
             split.by = "sample_type",
             ncol = 3) + 
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    ggtitle(NULL) +
    scale_color_manual(values = npg_extended) +
    coord_fixed(ratio = 1) +
    guides(color = guide_legend(title = NULL, override.aes = list(size = 5))) +
    theme(
        text = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 20, color = "black"),
        axis.text.y = element_text(size = 20, color = "black"),
        axis.title.x = element_text(size = 24, face = "bold", color = "black"),
        axis.title.y = element_text(size = 24, face = "bold", color = "black", margin = margin(r = 20)),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 20, face = "bold", color = "black"),
        legend.title = element_text(size = 20, face = "bold", color = "black"),
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
        plot.margin = margin(10, 50, 10, 10),
        strip.text = element_text(size = 16, face = "bold", 
                                 margin = margin(b = 15)),
        panel.spacing = unit(1.5, "lines")
    )

print(p)
dev.off()

#比例气泡图
# 计算每个sample_type和period1组合中的cell_type比例
cell_proportions <- fibroblasts@meta.data %>%
  group_by(sample_type, period1, cell_type) %>%
  summarise(count = n()) %>%
  group_by(sample_type, period1) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()

# 转换罗马数字为阿拉伯数字
cell_proportions <- cell_proportions %>%
  mutate(period1 = case_when(
    period1 == "Ⅰ" ~ "1",
    period1 == "Ⅱ" ~ "2", 
    period1 == "Ⅲ" ~ "3",
    TRUE ~ as.character(period1)
  ))

# 创建组合标签并设置正确的顺序
cell_proportions <- cell_proportions %>%
  mutate(
    sample_period = case_when(
      sample_type == "Tumor" ~ paste0("Tumor_", period1),
      sample_type == "mLN" ~ "mLN",
      sample_type == "Normal" ~ paste0("Normal_", period1),
      TRUE ~ paste0(sample_type, "_", period1)
    )
  )

# 设置样本顺序
sample_period_order <- c(
  "Tumor_1", "Tumor_2", "Tumor_3",
  "mLN",
  "Normal_1", "Normal_2", "Normal_3"
)

cell_proportions$sample_period <- factor(cell_proportions$sample_period, 
                                        levels = sample_period_order)

# 设置颜色
npg_pal <- pal_npg()(10)
cell_type_colors <- colorRampPalette(npg_pal)(length(unique(cell_proportions$cell_type)))

# 绘制气泡图
p <- ggplot(cell_proportions, aes(x = cell_type, y = sample_period, 
                                 size = proportion, color = cell_type)) +
  geom_point(alpha = 0.8) +  # 使用geom_point绘制气泡
  scale_size_continuous(name = "Proportion", 
                       range = c(1, 10),  # 调整气泡大小范围
                       limits = c(0, 1),
                       breaks = c(0.2, 0.4, 0.6, 0.8)) +  # 设置图例刻度
  scale_color_manual(values = cell_type_colors) +
  scale_y_discrete(limits = rev(levels(cell_proportions$sample_period))) +  # 反转y轴顺序
  labs(x = "Cell Type", y = "Sample Type & Period", 
       title = "Fibroblast Proportions Across Samples") +
  theme_classic() +  # 使用theme_classic，它自带坐标轴线
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
    axis.text.y = element_text(size = 11, face = "bold"),
    axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 12, face = "bold", margin = margin(r = 10)),
    legend.position = "right",  # 显示图例在右侧
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    # 移除网格线
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # 确保坐标轴框线显示
    axis.line = element_line(color = "black", linewidth = 0.5),
    panel.background = element_rect(fill = "white", color = "black", linewidth = 0.5),  # 添加面板边框
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  ) +
  guides(color = "none")  # 移除颜色图例，只保留大小图例

# 显示图形
print(p)

# 保存图形
ggsave("fibroblasts_bubble_plot.pdf", p, width = 25, height = 10, dpi = 300)