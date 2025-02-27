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
count_data <- merged_seurat_obj@meta.data %>%
  group_by(orig.ident, cell_type) %>% 
  summarise(count = n())

# 从原始元数据中获取 orig.ident 与 sample_type 的对应关系
sample_type_info <- merged_seurat_obj@meta.data %>%
  select(orig.ident, sample_type) %>%
  distinct()

# 合并 sample_type 信息到 count_data
count_data <- left_join(count_data, sample_type_info, by = "orig.ident")

# 分离肿瘤和正常样本的 orig.ident 并确定顺序
tumor_orig_ident <- count_data %>% 
  filter(sample_type == "tumor") %>% 
  pull(orig.ident) %>% 
  unique()
normal_orig_ident <- count_data %>% 
  filter(sample_type == "normal") %>% 
  pull(orig.ident) %>% 
  unique()
new_order <- c(tumor_orig_ident, normal_orig_ident)

npg_extended <- colorRampPalette(pal_npg("nrc")(10))(length(unique(count_data$cell_type)))

png("sample_count_stack.png", width = 6000, height = 3000, res = 300) 

ggplot(count_data, aes(x = orig.ident, y = count, fill = cell_type)) +
  scale_x_discrete(limits = new_order) +
  scale_fill_manual(values = npg_extended) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "", y = "Cell Count", fill = "Cell Type") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 18, angle = 90, hjust = 1, vjust = 0.5),
    axis.line = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 28),
    legend.title = element_text(size = 28)
  )

dev.off()

#cell proportion based on sample
proportion_data <- merged_seurat_obj@meta.data %>%
  group_by(orig.ident, cell_type) %>% 
  summarise(count = n()) %>% 
  mutate(proportion = count / sum(count))

# 从原始元数据中获取 orig.ident 与 sample_type 的对应关系
sample_type_info <- merged_seurat_obj@meta.data %>%
  select(orig.ident, sample_type) %>%
  distinct()

# 合并 sample_type 信息到 proportion_data
proportion_data <- left_join(proportion_data, sample_type_info, by = "orig.ident")

# 分离肿瘤和正常样本的 orig.ident 并确定顺序
tumor_orig_ident <- proportion_data %>% 
  filter(sample_type == "tumor") %>% 
  pull(orig.ident) %>% 
  unique()
normal_orig_ident <- proportion_data %>% 
  filter(sample_type == "normal") %>% 
  pull(orig.ident) %>% 
  unique()
new_order <- c(tumor_orig_ident, normal_orig_ident)

npg_extended <- colorRampPalette(pal_npg("nrc")(10))(length(unique(proportion_data$cell_type)))

png("sample_prop1.png", width = 6000, height = 3000, res = 300) 

ggplot(proportion_data, aes(x = orig.ident, y = proportion, fill = cell_type)) +
  scale_x_discrete(limits = new_order) +
  scale_fill_manual(values = npg_extended) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "", y = "Proportion", fill = "Cell Type") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 18, angle = 90, hjust = 1, vjust = 0.5),
    axis.line = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 28),
    legend.title = element_text(size = 28)
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

# 将样本类型（sample_type）和细胞类型（cell_type）分开
cell_counts <- cell_counts %>%
  left_join(sample_totals, by = "orig.ident")

# 对每个细胞类型构建列联表并进行卡方检验
chisq_results <- cell_counts %>%
  group_by(cell_type) %>%
  summarise(
    p_value = {
      # 提取正常样本和肿瘤样本的细胞数量和总细胞数
      normal_count = sum(count[sample_type == "normal"])
      tumor_count = sum(count[sample_type == "tumor"])
      normal_total = sum(total[sample_type == "normal"])
      tumor_total = sum(total[sample_type == "tumor"])
      
      # 构建 2x2 列联表
      cont_table <- matrix(
        c(normal_count, tumor_count, normal_total - normal_count, tumor_total - tumor_count),
        nrow = 2
      )
      
      # 进行卡方检验
      chisq.test(cont_table)$p.value
    }
  ) %>%
  ungroup()

# 添加显著性标记
chisq_results <- chisq_results %>%
  mutate(
    significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

# 绘制基于 count 的箱线图并添加显著性标记
png("箱图_count.png", width = 6000, height = 3000, res = 300)
ggplot(cell_counts, aes(x = cell_type, y = count, fill = sample_type)) +
  geom_boxplot() +  # 绘制箱线图
  geom_text(
    data = chisq_results, 
    aes(x = cell_type, y = max(cell_counts$count) * 0.9, label = significance),  # 调整 y 值
    size = 12, 
    vjust = 0.5,  # 调整 vjust 参数
    inherit.aes = FALSE  # 忽略父图层的 aes 映射
  ) +  # 添加显著性标记
  labs(x = "", y = "Cell Count", fill = "Sample Type") +  # 设置坐标轴和图例标题
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
  scale_y_continuous(limits = c(0, 6500)) 
dev.off()


identity_mapping <- c(
  "Sample1" = "Ⅱa",
  "Sample2" = "Ⅰa",
  "Sample3" = "Ⅲb",
  "Sample4" = "Ⅱa",
  "Sample5" = "Ⅲb",
  "Sample6" = "Ⅱa",
  "Sample7" = "Ⅲb",
  "Sample8" = "Ⅱb",
  "Sample9" = "Ⅱb",
  "Sample10" = "Ⅱb",
  "Sample11" = "Ⅱb",
  "Sample12" = "Ⅲa",
  "Sample13" = "Ⅰb",
  "Sample14" = "Ⅰb",
  "Sample15" = "Ⅰb",
  "Sample16" = "Ⅰb",
  "Sample17" = "Ⅲc",
  "Sample18" = "Ⅲc",
  "Sample19" = "Ⅲa",
  "Sample20" = "Ⅲc",
  "Sample21" = "Ⅲc",
  "Sample22" = "Ⅲa",
  "Sample23" = "Ⅲa",
  "Sample24" = "Ⅰa",
  "Sample25" = "0",
  "Sample26" = "0"
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

#patients
identity_mapping <- c(
  "Sample1" = "P14",
  "Sample2" = "P3",
  "Sample3" = "P15",
  "Sample4" = "P6",
  "Sample5" = "P11",
  "Sample6" = "P6",
  "Sample7" = "P11",
  "Sample8" = "P7",
  "Sample9" = "P7",
  "Sample10" = "P8",
  "Sample11" = "P8",
  "Sample12" = "P9",
  "Sample13" = "P4",
  "Sample14" = "P4",
  "Sample15" = "P5",
  "Sample16" = "P5",
  "Sample17" = "P12",
  "Sample18" = "P12",
  "Sample19" = "P9",
  "Sample20" = "P13",
  "Sample21" = "P13",
  "Sample22" = "P10",
  "Sample23" = "P10",
  "Sample24" = "P2",
  "Sample25" = "P1",
  "Sample26" = "P1"
)

patient <- identity_mapping[merged_seurat_obj@meta.data$orig.ident]
merged_seurat_obj@meta.data$patient <- patient

count_data <- merged_seurat_obj@meta.data %>%
  group_by(orig.ident, cell_type) %>% 
  summarise(count = n())

# 从原始元数据中获取 orig.ident、patient 与 sample_type 的对应关系
patient_info <- merged_seurat_obj@meta.data %>%
  select(orig.ident, patient, sample_type) %>%
  distinct()

# 合并 patient 信息到 count_data
count_data <- left_join(count_data, patient_info, by = "orig.ident")

# 定义新的顺序
# 这里将 patient 和 sample_type 组合成一个新的因子
count_data$patient_sample_type <- paste(count_data$patient, count_data$sample_type, sep = "_")
new_order <- count_data$patient_sample_type %>% unique()

# 按照肿瘤在前、正常在后的顺序对新顺序进行排序
tumor_order <- new_order[grepl("tumor", new_order)]
normal_order <- new_order[grepl("normal", new_order)]
new_order <- c(tumor_order, normal_order)

npg_extended <- colorRampPalette(pal_npg("nrc")(10))(length(unique(count_data$cell_type)))

png("sample_count_stack.png", width = 6000, height = 3000, res = 300) 

ggplot(count_data, aes(x = patient_sample_type, y = count, fill = cell_type)) +
  scale_x_discrete(limits = new_order) +
  scale_fill_manual(values = npg_extended) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "", y = "Cell Count", fill = "Cell Type") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 18, angle = 90, hjust = 1, vjust = 0.5),
    axis.line = element_line(color = "black"),
    axis.title.y = element_text(size = 28),
    axis.text.y = element_text(size = 18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 28),
    legend.title = element_text(size = 28)
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