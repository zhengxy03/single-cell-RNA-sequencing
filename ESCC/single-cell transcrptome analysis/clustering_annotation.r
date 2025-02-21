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
png("clusters.png", width = 3000, height = 3000, res = 300)  # 设置高分辨率和尺寸
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
    axis.title.y = element_text(size = 36, face = "bold", color = "black"),
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
png("annotation.png", width = 3000, height = 3000, res = 300)  # 设置高分辨率和尺寸
DimPlot(merged_seurat_obj, reduction = "umap", label = TRUE, pt.size = 1, group.by = "cell_type") +
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
    axis.text.x = element_text(size = 14, angle = 90, hjust = 1, vjust = 0.5),  # 调整横轴字体大小、角度和对齐方式
    axis.line = element_line(color = "black"),  # 设置坐标轴颜色
    panel.grid.major = element_blank(),  # 移除主要网格线
    panel.grid.minor = element_blank()   # 移除次要网格线
  )

dev.off()

#cell proportion based on sampletype
proportion_data <- merged_seurat_obj@meta.data %>%
    group_by(sample_type, cell_type) %>% summarise(count = n()) %>% mutate(proportion = count / sum(count))

npg_extended <- colorRampPalette(npg_pal)(13)
png("sampletype_prop1.png", width = 3000, height = 3000, res = 300)  # 设置高分辨率和尺寸
ggplot(proportion_data, aes(x = "", y = proportion, fill = cell_type)) +
    geom_bar(stat = "identity", width = 1) +         # 堆叠柱状图
    coord_polar(theta = "y") +                       # 转换为饼图
    scale_fill_manual(values = npg_extended) +
    theme_void() +                                   # 空白背景
    theme(
        legend.position = "right",                     # 图例放在右侧
        plot.title = element_blank(),                  # 移除标题
        # 分面标签设置
        strip.placement = "outside",                   # 标签放在绘图区域外
        strip.text = element_text(                     # 标签样式
            size = 12,                                   # 字体大小
            face = "bold",                               # 加粗
            margin = margin(b = 10)                      # 下方留白（避免与饼图重叠）
        )
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
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5),  # 调整横轴标签角度和对齐方式
        axis.line = element_line(color = "black"),  # 设置坐标轴颜色
        panel.grid.major = element_blank(),  # 移除主要网格线
        panel.grid.minor = element_blank()   # 移除次要网格线
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
    label.y = 1.1,             # 调整标签的垂直位置
    size = 4                   # 标签字体大小
  ) +
  labs(x = "", y = "Proportion", fill = "Sample Type") +  # 设置坐标轴和图例标题
  theme_classic() +  # 使用经典主题
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),  # 调整横轴标签角度
    axis.line = element_line(color = "black"),  # 设置坐标轴颜色
    panel.grid.major = element_blank(),  # 移除主要网格线
    panel.grid.minor = element_blank(),  # 移除次要网格线
    legend.position = "right"  # 设置图例位置在右侧
  ) +
  scale_fill_npg()
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
    theme_void() +                                   # 空白背景
    theme(
        legend.position = "right",                     # 图例放在右侧
        plot.title = element_blank(),                  # 移除标题
        # 分面标签设置
        strip.placement = "outside",                   # 标签放在绘图区域外
        strip.text = element_text(                     # 标签样式
            size = 12,                                   # 字体大小
            face = "bold",                               # 加粗
            margin = margin(b = 10)                      # 下方留白（避免与饼图重叠）
        )
    ) +
    facet_wrap(
        ~ period,
        ncol = 4,
        strip.position = "bottom"                      # 标签放在下方
    )
dev.off()
