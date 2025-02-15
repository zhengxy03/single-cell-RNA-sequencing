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
merged_seurat_obj <- FindClusters(merged_seurat_obj, resolution = 0.5)
merged_seurat_obj <- RunUMAP(merged_seurat_obj, reduction = "harmony", dims = 1:20)

DimPlot(merged_seurat_obj, reduction = "umap", label = TRUE, pt.size = 1) +
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
significant_markers <- subset(markers, p_val_adj < 0.05)
significant_markers <- significant_markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.csv(significant_markers,"marker_top.csv")

identity_mapping <- c(
  "0" = "T cell",
  "1" = "B cell",
  "2" = "T cell",
  "3" = "Fibroblast",
  "4" = "Plasma",
  "5" = "Macrophage",
  "6" = "Dendritic cell",
  "7" = "Endothelial cell",
  "8" = "Monocyte",
  "9" = "Epithelial cell",
  "10" = "Fibroblast",
  "11" = "Mast cell",
  "12" = "Epithelial cell",
  "13" = "Proliferating cell",
  "14" = "Pericyte",
  "15" = "T cell",
  "16" = "T cell",
  "17" = "Fibroblast",
  "18" = "T cell",
  "19" = "Plasma",
  "20" = "Plasma",
  "21" = "Proliferating cell",
  "22" = "Epithelial cell"
)

identity_mapping <- c(
  "0" = "T cell",
  "1" = "B cell",
  "2" = "T cell",
  "3" = "Fibroblast",
  "4" = "Plasma",
  "5" = "Macrophage",
  "6" = "Dendritic cell",
  "7" = "Endothelial cell",
  "8" = "Monocyte",
  "9" = "Epithelial cell",
  "10" = "Fibroblast",
  "11" = "Mast cell",
  "12" = "Epithelial cell",
  "13" = "Proliferating cell",
  "14" = "Pericyte",
  "15" = "T cell",
  "16" = "T cell",
  "17" = "Fibroblast",
  "18" = "T cell",
  "19" = "Plasma",
  "20" = "Plasma",
  "21" = "Proliferating cell",
  "22" = "Epithelial cell"
)

cell_type <- identity_mapping[merged_seurat_obj@meta.data$seurat_clusters]
merged_seurat_obj@meta.data$cell_type <- cell_type

DimPlot(merged_seurat_obj, reduction = "umap", label = TRUE, pt.size = 1, group.by = "cell_type") +
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

#plot according to sample types
DimPlot(merged_seurat_obj, reduction = "umap", group.by = "sample_type", pt.size = 1) +
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

ggplot(proportion_data, aes(x = orig.ident, y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack") +  # 堆叠柱状图
  labs(x = "Sample", y = "Proportion", fill = "Cell Type") +  # 设置坐标轴和图例标题
  theme_classic() +  # 使用经典主题（独立坐标轴）
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),  # 调整横轴标签角度和对齐方式
    axis.line = element_line(color = "black"),  # 设置坐标轴颜色
    panel.grid.major = element_blank(),  # 移除主要网格线
    panel.grid.minor = element_blank()   # 移除次要网格线
  )

#cell proportion based on sampletype
proportion_data <- merged_seurat_obj@meta.data %>%
  group_by(sample_type, cell_type) %>% summarise(count = n()) %>% mutate(proportion = count / sum(count))

ggplot(proportion_data, aes(x = sample_type, y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack") +  # 堆叠柱状图
  labs(x = "Sample Type", y = "Proportion", fill = "Cell Type") +  # 设置坐标轴和图例标题
  theme_classic() +  # 使用经典主题（独立坐标轴）
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),  # 调整横轴标签角度和对齐方式
    axis.line = element_line(color = "black"),  # 设置坐标轴颜色
    panel.grid.major = element_blank(),  # 移除主要网格线
    panel.grid.minor = element_blank()   # 移除次要网格线
  )

ggplot(proportion_data, aes(x = "", y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", width = 1) +         # 堆叠柱状图
  coord_polar(theta = "y") +                       # 转换为饼图
  labs(title = "Proportion of Cell Types", fill = "Cell Type") +
  theme_void() +                                   # 空白背景
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5),        # 标题居中
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

ggplot(cell_proportion, aes(x = cell_type, y = proportion, fill = sample_type)) +
  geom_boxplot() +  # 绘制箱线图
  stat_compare_means(  # 添加显著性标记
    aes(group = sample_type),  # 按 sample_type 分组比较
    method = "wilcox.test",    # 使用 Wilcoxon 检验
    label = "p.signif",        # 显示显著性符号（如 *）
    label.y = 1.1,             # 调整标签的垂直位置
    size = 4                   # 标签字体大小
  ) +
  labs(x = "Cell Type", y = "Proportion", fill = "Sample Type") +  # 设置坐标轴和图例标题
  theme_classic() +  # 使用经典主题
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # 调整横轴标签角度
    axis.line = element_line(color = "black"),  # 设置坐标轴颜色
    panel.grid.major = element_blank(),  # 移除主要网格线
    panel.grid.minor = element_blank(),  # 移除次要网格线
    legend.position = "right"  # 设置图例位置在右侧
  ) +
  scale_fill_manual(values = c("tumor" = "red", "normal" = "blue"))