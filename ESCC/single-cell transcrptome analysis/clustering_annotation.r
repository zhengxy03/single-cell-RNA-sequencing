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
library(ggsci)
merged_seurat_obj <- FindClusters(merged_seurat_obj, resolution = 0.5)
merged_seurat_obj <- RunUMAP(merged_seurat_obj, reduction = "harmony", dims = 1:20)

npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(23)

# 获取图例的个数和名称长度
seurat_clusters <- as.character(unique(merged_seurat_obj@meta.data$seurat_clusters))  # 转换为字符向量
num_legend_items <- length(seurat_clusters)  # 图例的个数
max_label_length <- max(nchar(seurat_clusters))  # 图例名称的最大长度

# 动态计算图片尺寸
base_width <- 3000  # 基础宽度
base_height <- 3000  # 基础高度
legend_width_factor <- 100  # 每个图例项增加的宽度
label_length_factor <- 10  # 每个字符增加的宽度

# 计算动态宽度
dynamic_width <- base_width + (num_legend_items * legend_width_factor) + (max_label_length * label_length_factor)

# 导出图片
png("clusters.png", width = dynamic_width, height = base_height, res = 300)
DimPlot(merged_seurat_obj, reduction = "umap", label = TRUE, pt.size = 1, label.size = 8) +
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

#plot on samples
npg_extended <- colorRampPalette(npg_pal)(28)
png("sample.png", width = 3000, height = 3000, res = 300)  # 设置高分辨率和尺寸
DimPlot(merged_seurat_obj, reduction = "umap", label = FALSE, pt.size = 1, group.by = "orig.ident") +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  ggtitle(NULL) +
  scale_color_manual(values = npg_extended) +
  coord_fixed(ratio = 1) +
  guides(color = guide_legend(title = NULL, override.aes = list(size = 5))) +
  theme(
    text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 14, face = "bold", color = "black"),
    axis.title.y = element_text(size = 14, face = "bold", color = "black"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 10, face = "bold", color = "black"),
    legend.title = element_text(size = 10, face = "bold", color = "black"),
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
#markers <- FindAllMarkers(object = merged_seurat_obj, 
                                  test.use = "roc", 
                                  only.pos = TRUE, 
                                  min.pct = 0.1, 
                                  thresh.use = 0.25)

#significant_markers <- subset(markers, myAUC > 0.7)
#write.csv(significant_markers,"marker.csv")

markers <- FindAllMarkers(merged_seurat_obj, 
                          only.pos = TRUE, 
                          min.pct = 0.25, 
                          logfc.threshold = 0.25, 
                          test.use = "wilcox")
#significant markers
library(dplyr)
significant_markers <- subset(markers, p_val_adj < 0.05)
significant_markers <- significant_markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.csv(significant_markers,"marker_top.csv")

identity_mapping <- c(
    "0" = "T cell",
    "1" = "B cell",
    "2" = "T cell",
    "3" = "Fibroblast",
    "4" = "Plasma",
    "5" = "Monocyte",
    "6" = "Endothelial cell",
    "7" = "Macrophage",
    "8" = "Dendritic cell",
    "9" = "Fibroblast",
    "10" = "Epithelial cell",
    "11" = "Mast cell",
    "12" = "Epithelial cell",
    "13" = "Pericyte",
    "14" = "Proliferating cell",
    "15" = "T cell",
    "16" = "NK cell",
    "17" = "Fibroblast",
    "18" = "T cell",
    "19" = "T cell",
    "20" = "Plasma",
    "21" = "Plasma",
    "22" = "Epithelial cell"
)

cell_type <- identity_mapping[merged_seurat_obj@meta.data$seurat_clusters]
merged_seurat_obj@meta.data$cell_type <- cell_type

npg_extended <- colorRampPalette(npg_pal)(23)

# 获取图例的个数和名称长度
cell_types <- unique(merged_seurat_obj@meta.data$cell_type)
num_legend_items <- length(cell_types)  # 图例的个数
max_label_length <- max(nchar(cell_types))  # 图例名称的最大长度

# 动态计算图片尺寸
base_width <- 3000  # 基础宽度
base_height <- 3000  # 基础高度
legend_width_factor <- 100  # 每个图例项增加的宽度
label_length_factor <- 10  # 每个字符增加的宽度

# 计算动态宽度
dynamic_width <- base_width + (num_legend_items * legend_width_factor) + (max_label_length * label_length_factor)

# 导出图片
png("annotation.png", width = dynamic_width, height = base_height, res = 300)
DimPlot(merged_seurat_obj, reduction = "umap", label = TRUE, pt.size = 1, group.by = "cell_type", label.size = 8) +
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

#plot according to sample types
png("sampletype.png", width = 3000, height = 3000, res = 300)  # 设置高分辨率和尺寸
DimPlot(merged_seurat_obj, reduction = "umap", label = FALSE, pt.size = 1, group.by = "sample_type") +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  ggtitle(NULL) +
  scale_color_npg() +
  coord_fixed(ratio = 1) +
  guides(color = guide_legend(title = NULL, override.aes = list(size = 5))) +
  theme(
    text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 14, face = "bold", color = "black"),
    axis.title.y = element_text(size = 14, face = "bold", color = "black"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 10, face = "bold", color = "black"),
    legend.title = element_text(size = 10, face = "bold", color = "black"),
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

#plot according to sample sources
DimPlot(merged_seurat_obj, reduction = "umap", group.by = "sample_sources", pt.size = 1) +
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

#plot according to sample
DimPlot(merged_seurat_obj, reduction = "umap", group.by = "orig.ident", pt.size = 1) +
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

#cell proportion based on sample
proportion_data <- merged_seurat_obj@meta.data %>%
  group_by(orig.ident, cell_type) %>% summarise(count = n()) %>% mutate(proportion = count / sum(count))

npg_extended <- colorRampPalette(npg_pal)(13)
png("sample_prop1.png", width = 6000, height = 3000, res = 300) 

ggplot(proportion_data, aes(x = orig.ident, y = proportion, fill = cell_type)) +
  scale_fill_manual(values = npg_extended) +
  geom_bar(stat = "identity", position = "stack") +  # 堆叠柱状图
  labs(x = "", y = "Proportion", fill = "Cell Type") +  # 设置坐标轴和图例标题
  theme_classic() +  # 使用经典主题（独立坐标轴）
  theme(
    axis.text.x = element_text(size = 18, angle = 90, hjust = 1, vjust = 0.5),  # 调整横轴字体大小、角度和对齐方式
    axis.line = element_line(color = "black"),  # 设置坐标轴颜色
    panel.grid.major = element_blank(),  # 移除主要网格线
    panel.grid.minor = element_blank(),   # 移除次要网格线
    legend.text = element_text(size = 28), # 图例文本字体大小
    legend.title = element_text(size = 28) # 图例标题字体大小（如果有标题的话，这里原代码设置为NULL）
  )

dev.off()

#cell proportion based on sampletype
proportion_data <- merged_seurat_obj@meta.data %>%
    group_by(sample_type, cell_type) %>% summarise(count = n()) %>% mutate(proportion = count / sum(count))

npg_extended <- colorRampPalette(npg_pal)(13)
png("sampletype_prop1.png", width = 6000, height = 3000, res = 300)  # 设置高分辨率和尺寸
ggplot(proportion_data, aes(x = "", y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", width = 1) +         # 堆叠柱状图
  coord_polar(theta = "y") +                       # 转换为饼图
  scale_fill_manual(values = npg_extended) +       # 使用自定义颜色
  theme_void() +                                   # 空白背景
  labs(fill = "Cell Type") +
  theme(
    legend.position = "right",                     # 图例放在右侧
    plot.title = element_blank(),                  # 移除标题
    # 分面标签设置
    strip.placement = "outside",                   # 标签放在绘图区域外
    strip.text = element_text(                     # 标签样式
      size = 40,                                   # 字体大小
      face = "bold",                               # 加粗
      margin = margin(b = 10)                      # 下方留白（避免与饼图重叠）
    ),
    legend.text = element_text(size = 36),         # 图例文本字体大小
    legend.title = element_text(size = 36)                 # 移除图例标题
  ) +
  facet_wrap(
    ~ sample_type,
    ncol = 2,
    strip.position = "bottom"                      # 标签放在下方
  )
dev.off()

png("sampletype_prop2.png", width = 3000, height = 3000, res = 300)  # 设置高分辨率和尺寸
ggplot(proportion_data, aes(x = sample_type, y = proportion, fill = cell_type)) +
    scale_fill_manual(values = npg_extended) +
    geom_bar(stat = "identity", position = "stack") +  # 堆叠柱状图
    labs(x = "Sample Type", y = "Proportion", fill = "Cell Type") +  # 设置坐标轴和图例标题
    theme_classic() +  # 使用经典主题（独立坐标轴）
    theme(
        axis.text.x = element_text(hjust = 0.5, vjust = 0.5, size = 28),  # 调整横轴标签角度和对齐方式
        axis.line = element_line(color = "black"),  # 设置坐标轴颜色
        panel.grid.major = element_blank(),  # 移除主要网格线
        panel.grid.minor = element_blank(),   # 移除次要网格线
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 40),

        axis.text.y = element_text(size = 28),
        legend.text = element_text(size = 36),         # 图例文本字体大小
        legend.title = element_text(size = 40)             
    )
dev.off()

#cell proportion based on samplesources
proportion_data <- merged_seurat_obj@meta.data %>%
  group_by(sample_sources, cell_type) %>% summarise(count = n()) %>% mutate(proportion = count / sum(count))

ggplot(proportion_data, aes(x = sample_sources, y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack") +  # 堆叠柱状图
  labs(x = "Sample Sources", y = "Proportion", fill = "Cell Type") +  # 设置坐标轴和图例标题
  theme_classic() +  # 使用经典主题（独立坐标轴）
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),  # 调整横轴标签角度和对齐方式
    axis.line = element_line(color = "black"),  # 设置坐标轴颜色
    panel.grid.major = element_blank(),  # 移除主要网格线
    panel.grid.minor = element_blank()   # 移除次要网格线
  )

#smaple source proportion
proportion_data <- merged_seurat_obj@meta.data %>%
  group_by(sample_sources) %>% summarise(count = n()) %>% mutate(proportion = count / sum(count))

ggplot(proportion_data, aes(x = "", y = proportion, fill = sample_sources)) +
  geom_bar(stat = "identity", width = 1) +  # 堆叠柱状图
  coord_polar(theta = "y") +                # 转换为饼图
  labs(
    title = "Sample Source Proportions",
    fill = "Sample Source",
    x = NULL,
    y = NULL
  ) +
  theme_void() +                           # 移除坐标轴和背景
  theme(
    plot.title = element_text(hjust = 0.5),  # 标题居中
    legend.position = "right"               # 图例在右侧
  )

#箱线图
library(ggpubr)
cell_counts <- merged_seurat_obj@meta.data %>%
  group_by(orig.ident, cell_type, sample_type) %>%
  summarise(count = n()) %>%
  ungroup()


sample_totals <- cell_counts %>%
  group_by(orig.ident) %>%
  summarise(total = sum(count)) %>%
  ungroup()

cell_proportion <- cell_counts %>%
  left_join(sample_totals, by = "orig.ident") %>%
  mutate(proportion = count / total)

png("箱图.png", width = 6000, height = 3000, res = 300)
ggplot(cell_proportion, aes(x = cell_type, y = proportion, fill = sample_type)) +
  geom_boxplot() +  # 绘制箱线图
  stat_compare_means(  # 添加显著性标记
    aes(group = sample_type),  # 按 sample_type 分组比较
    method = "wilcox.test",    # 使用 Wilcoxon 检验
    label = "p.signif",        # 显示显著性符号（如 *）
    label.y = 0.6,             # 调整标签的垂直位置
    size = 12                   # 标签字体大小
  ) +
  labs(x = "", y = "Proportion", fill = "Sample Type") +  # 设置坐标轴和图例标题
  theme_classic() +  # 使用经典主题
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 28),  # 调整横轴标签角度
    axis.text.y = element_text(size = 28),
    axis.title.y = element_text(size = 40),
    axis.line = element_line(color = "black"),  # 设置坐标轴颜色
    panel.grid.major = element_blank(),  # 移除主要网格线
    panel.grid.minor = element_blank(),  # 移除次要网格线
    legend.text = element_text(size = 36),         # 图例文本字体大小
    legend.title = element_text(size = 40),
    legend.position = "right"  # 设置图例位置在右侧
  ) +
  scale_fill_npg() +
  scale_y_continuous(limits = c(0, 0.6)) 
dev.off()


png("count_箱图.png", width = 6000, height = 3000, res = 300)
ggplot(cell_proportion, aes(x = cell_type, y = count, fill = sample_type)) +
    geom_boxplot() +  # 绘制箱线图
    stat_compare_means(  # 添加显著性标记
        aes(group = sample_type),  # 按 sample_type 分组比较
        method = "wilcox.test",    # 使用 Wilcoxon 检验
        label = "p.signif",        # 显示显著性符号（如 *）
    # 调整标签的垂直位置
        size = 12                   # 标签字体大小
    ) +
    labs(x = "", y = "Count", fill = "Sample Type") +  # 设置坐标轴和图例标题
    theme_classic() +  # 使用经典主题
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 28),  # 调整横轴标签角度
        axis.text.y = element_text(size = 28),
        axis.title.y = element_text(size = 40),
        axis.line = element_line(color = "black"),  # 设置坐标轴颜色
        panel.grid.major = element_blank(),  # 移除主要网格线
        panel.grid.minor = element_blank(),  # 移除次要网格线
        legend.text = element_text(size = 36),         # 图例文本字体大小
        legend.title = element_text(size = 40),
        legend.position = "right"  # 设置图例位置在右侧
    ) +
    scale_fill_npg() +
    scale_y_continuous(limits = c(0, 4000)) 
dev.off()

identity_mapping <- c(
  "Sample1" = "Ⅲb",
  "Sample2" = "Ⅱa",
  "Sample3" = "Ⅲb",
  "Sample4" = "Ⅱa",
  "Sample5" = "Ⅰa",
  "Sample6" = "Ⅱa",
  "Sample7" = "Ⅲb",
  "Sample8" = "Ⅱa",
  "Sample9" = "Ⅲb",
  "Sample10" = "Ⅱb",
  "Sample11" = "Ⅱb",
  "Sample12" = "Ⅱb",
  "Sample13" = "Ⅱb",
  "Sample14" = "Ⅲa",
  "Sample15" = "Ⅰb",
  "Sample16" = "Ⅰb",
  "Sample17" = "Ⅰb",
  "Sample18" = "Ⅰb",
  "Sample19" = "Ⅲc",
  "Sample20" = "Ⅲc",
  "Sample21" = "Ⅲa",
  "Sample22" = "Ⅲc",
  "Sample23" = "Ⅲc",
  "Sample24" = "Ⅲa",
  "Sample25" = "Ⅲa",
  "Sample26" = "Ⅰa",
  "Sample27" = "0",
  "Sample28" = "0"
)
period <- identity_mapping[merged_seurat_obj@meta.data$orig.ident]
merged_seurat_obj@meta.data$period <- period


png("period.png", width = 3000, height = 3000, res = 300)  # 设置高分辨率和尺寸
DimPlot(merged_seurat_obj, reduction = "umap", label = FALSE, pt.size = 1, group.by = "period") +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  ggtitle(NULL) +
  scale_color_npg() +
  coord_fixed(ratio = 1) +
  guides(color = guide_legend(title = NULL, override.aes = list(size = 5))) +
  theme(
    text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 14, face = "bold", color = "black"),
    axis.title.y = element_text(size = 14, face = "bold", color = "black"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 10, face = "bold", color = "black"),
    legend.title = element_text(size = 10, face = "bold", color = "black"),
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


proportion_data <- merged_seurat_obj@meta.data %>%
    group_by(period, cell_type) %>% summarise(count = n()) %>% mutate(proportion = count / sum(count))

npg_extended <- colorRampPalette(npg_pal)(13)
png("period_prop1.png", width = 6000, height = 3000, res = 300)  # 设置高分辨率和尺寸
ggplot(proportion_data, aes(x = "", y = proportion, fill = cell_type)) +
    geom_bar(stat = "identity", width = 1) +         # 堆叠柱状图
    coord_polar(theta = "y") +                       # 转换为饼图
    scale_fill_manual(values = npg_extended) +
    labs(fill = "Cell Type") +
    theme_void() +                                   # 空白背景
    theme(
        legend.position = "right",                     # 图例放在右侧
        plot.title = element_blank(),                  # 移除标题
        # 分面标签设置
        strip.placement = "outside",                   # 标签放在绘图区域外
        strip.text = element_text(                     # 标签样式
            size = 40,                                   # 字体大小
            face = "bold",                               # 加粗
            margin = margin(b = 10)                      # 下方留白（避免与饼图重叠）
        ),
        legend.text = element_text(size = 36),         # 图例文本字体大小
        legend.title = element_text(size = 36)  
    ) +
    facet_wrap(
        ~ period,
        ncol = 4,
        strip.position = "bottom"                      # 标签放在下方
    )
dev.off()

png("period_prop2.png", width = 6000, height = 3000, res = 300)
cell_counts <- merged_seurat_obj@meta.data %>%
  group_by(orig.ident, cell_type, period) %>%
  summarise(count = n()) %>%
  ungroup()

ggplot(cell_counts, aes(x = cell_type, y = count, fill = period)) +
  geom_bar(stat = "identity", position = "stack") +  # 堆叠柱状图
  labs(x = "Cell Type", y = "Count", fill = "Stage") +  # 设置坐标轴和图例标题
  scale_fill_npg() +  # 使用ggsci的npg配色
  theme_classic() +  # 使用经典主题
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 28),  # 调整横轴标签角度
    axis.text.y = element_text(size = 28),
    axis.title.x = element_blank(),  # 横坐标标题字体大小
    axis.title.y = element_text(size = 40),  # 纵坐标标题字体大小
    legend.text = element_text(size = 36),  # 图例文本字体大小
    legend.title = element_text(size = 40),  # 图例标题字体大小
    legend.position = "right"  # 将图例放在顶部
  )
dev.off()

tumor <- subset(merged_seurat_obj, subset = sample_type == "tumor")
png("t_period_prop2.png", width = 6000, height = 3000, res = 300)
cell_counts <- tumor@meta.data %>%
  group_by(orig.ident, cell_type, period) %>%
  summarise(count = n()) %>%
  ungroup()

ggplot(cell_counts, aes(x = cell_type, y = count, fill = period)) +
  geom_bar(stat = "identity", position = "stack") +  # 堆叠柱状图
  labs(x = "Cell Type", y = "Count", fill = "Stage") +  # 设置坐标轴和图例标题
  scale_fill_npg() +  # 使用ggsci的npg配色
  theme_classic() +  # 使用经典主题
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 28),  # 调整横轴标签角度
    axis.text.y = element_text(size = 28),
    axis.title.x = element_blank(),  # 横坐标标题字体大小
    axis.title.y = element_text(size = 40),  # 纵坐标标题字体大小
    legend.text = element_text(size = 36),  # 图例文本字体大小
    legend.title = element_text(size = 40),  # 图例标题字体大小
    legend.position = "right"  # 将图例放在顶部
  )
dev.off()

fibro_data <- subset(cell_proportion, cell_type == "Fibroblast")

#DEGs

genes_to_plot <- c(
    # 淋巴细胞系
    "CD3D", "CD3E", "NKG7",          # T细胞
    "BANK1", "CD79A", "MS4A1",        # B细胞
    "GNLY", "KLRC1", "KLRB1",       # NK细胞
    "IGHG1", "IGHG3", "JCHAIN",      # 浆细胞
    
    # 髓系免疫细胞
    "C1QA", "C1QB", "APOC1",         # 巨噬细胞
    "FCN1", "CD300E", "SLC11A1",          # 单核细胞
    "CD86", "PKIB", "CLEC10A",      # 树突细胞（修正HAL-DRA为HLA-DRA）
    
    # 结构细胞
    "SFN", "KRT19", "KRT17",       # 上皮细胞
    "SFRP2", "COL1A2", "FGF7",     # 成纤维细胞
    "CLDN5", "PECAM1", "RAMP2",      # 内皮细胞
    "RGS5", "MCAM", "ACTA2",         # 周细胞
    
    # 特殊功能细胞
    "CPA3", "TPSAB1", "TPSB2",       # 肥大细胞
    "TOP2A", "STMN1", "MKI67"         # 增殖细胞
)
png("dotplot.png", width = 8000, height = 3000, res = 300)  # 设置高分辨率和尺寸

DotPlot(merged_seurat_obj, 
        features = genes_to_plot, 
        group.by = "cell_type",
        dot.scale = 15) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 28),  # 横坐标基因名字体大小
    axis.text.y = element_text(size = 28),  # 纵坐标细胞类型字体大小
    axis.title.x = element_blank(),  # 去除横坐标标题
    axis.title.y = element_text(size = 40),
    panel.grid.major = element_line(color = "grey70", linewidth = 0.5),  # 添加网格线
    panel.grid.minor = element_line(color = "grey80", linewidth = 0.3),  # 细网格线 
    legend.text = element_text(size = 28, face = "bold", color = "black"),
    legend.title = element_text(size = 28, face = "bold", color = "black"),
    legend.position = "right",
    legend.box.margin = margin(0, 0, 0, 0),
    legend.key = element_blank(),
    legend.background = element_blank()
  ) +
  scale_color_gsea()  # 使用ggsci的GSEA配色
dev.off()

FeaturePlot(merged_seurat_obj, features = genes_to_plot)