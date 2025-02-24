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




immune_cell_types <- c("T cell", "B cell", "Macrophages", "NK cell", "Dendritic cell", "Monocyte", "Plasma")
immune_seurat <- subset(merged_seurat_obj, cell_type %in% immune_cell_types)

samples <- unique(immune_seurat$orig.ident)

# 计算每个样本中T细胞的占比
tcell_percent <- sapply(samples, function(s) {
  sample_cells <- subset(immune_seurat, orig.ident == s)
  tcell_count <- sum(sample_cells$cell_type == "T cell")
  total_cells <- ncol(sample_cells)
  return((tcell_count / total_cells) * 100)
})

# 创建数据框
df_tcell <- data.frame(
  Sample = samples,
  Tcell_Percent = tcell_percent,
  Group = immune_seurat@meta.data[match(samples, immune_seurat$orig.ident), "sample_type"]  # 添加分组信息
)

tcell_matrix <- matrix(df_tcell$Tcell_Percent, nrow = 1)
rownames(tcell_matrix) <- "T cells"
colnames(tcell_matrix) <- df_tcell$Sample