library(ggplot2)
library(ggsci)
library(patchhwork)
library(Seurat)
library(dplyr)

npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(17)
#sample_type
# 获取唯一的 sample_type 值，并按字母顺序排序
unique_sample_types <- sort(unique(epi@meta.data$sample_type))  # 按字母顺序排序

# 如果需要手动指定顺序，可以这样做：
# unique_sample_types <- c("Tumor", "Normal", "Blood")

# 创建一个空列表，用于存储每个 sample_type 的图
plot_list <- list()

# 遍历每个 sample_type，生成单独的图
for (sample in unique_sample_types) {
    # 创建颜色映射：当前 sample_type 为 npg 红色，其他为灰色
    color_mapping <- setNames(
        ifelse(unique_sample_types == sample, pal_npg()(1), "gray"),  # npg 红色
        unique_sample_types
    )
    
    # 绘制 DimPlot
    p <- DimPlot(epi, reduction = "umap", label = FALSE, pt.size = 1, group.by = "sample_type") +
        scale_color_manual(values = color_mapping) +  # 使用自定义颜色映射
        ggtitle(sample) +  # 设置标题为当前 sample_type
        theme(
            legend.position = "none",  # 隐藏图例
            plot.title = element_text(size = 24, face = "bold", hjust = 0.5),  # 标题样式
            axis.title.x = element_text(size = 20, face = "bold", color = "black"),  # X 轴标题
            axis.title.y = element_text(size = 20, face = "bold", color = "black"),  # Y 轴标题
            axis.text.x = element_text(size = 16, color = "black"),  # X 轴刻度
            axis.text.y = element_text(size = 16, color = "black")   # Y 轴刻度
        )
    
    # 将图添加到列表中
    plot_list[[sample]] <- p
}

# 使用 patchwork 将图形排列在一起（每行两个）
combined_plot <- wrap_plots(plot_list, ncol = 2)  

# 保存为 PNG 文件
png("epi_sample_type_umap.png", width = 6000, height = 3000, res = 300)
print(combined_plot)
dev.off()

#period1 and sample_type
unique_periods <- sort(unique(epi@meta.data$period1))

# 获取唯一的 sample_type 值，并按字母顺序排序
unique_sample_types <- sort(unique(epi@meta.data$sample_type))

# 创建一个空列表，用于存储所有图
all_plot_list <- list()

# 遍历每个 period1
for (period in unique_periods) {
    # 遍历每个 sample_type
    for (sample in unique_sample_types) {
        # 创建一个辅助列，将 period1 和 sample_type 组合起来
        epi$combined_group <- paste(epi@meta.data$period1, epi@meta.data$sample_type, sep = "_")
        current_group <- paste(period, sample, sep = "_")
        
        # 创建颜色映射：当前 period1 和 sample_type 组合为 npg 红色，其他为灰色
        color_mapping <- setNames(
            ifelse(epi$combined_group == current_group, pal_npg()(1), "gray"),
            epi$combined_group
        )

        # 绘制 DimPlot
        p <- DimPlot(epi, reduction = "umap", label = FALSE, pt.size = 1, group.by = "combined_group") +
            scale_color_manual(values = color_mapping) +
            ggtitle(paste(period, sample, sep = " - ")) +
            theme(
                legend.position = "none",
                plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
                axis.title.x = element_text(size = 20, face = "bold", color = "black"),
                axis.title.y = element_text(size = 20, face = "bold", color = "black"),
                axis.text.x = element_text(size = 16, color = "black"),
                axis.text.y = element_text(size = 16, color = "black")
            )

        # 将图添加到总列表中
        all_plot_list[[paste(period, sample, sep = "_")]] <- p
    }
}

# 使用 patchwork 将图形排列在一起
# 这里根据实际图形数量动态调整每行的图形数量


n_periods <- length(unique_periods)
n_samples <- length(unique_sample_types)

# 创建一个新的列表来存储重新排列后的图形
reordered_plot_list <- list()

# 重新排列顺序
for (i in 1:n_samples) {
    for (j in 1:n_periods) {
        # 计算原始列表中的索引
        original_index <- (j - 1) * n_samples + i
        # 将图形添加到新的列表中
        reordered_plot_list[[length(reordered_plot_list) + 1]] <- all_plot_list[[original_index]]
    }
}

# 使用 wrap_plots 将图形排列成两行四列
combined_plot <- wrap_plots(reordered_plot_list, ncol = 3)

# 保存为 PNG 文件
png("epi_period1_sample_type_umap.png", width = 6000, height = 3000, res = 300)
print(combined_plot)
dev.off()  

#period2
unique_periods <- sort(unique(epi@meta.data$period2))

# 获取唯一的 sample_type 值，并按字母顺序排序
unique_sample_types <- sort(unique(epi@meta.data$sample_type))

# 创建一个空列表，用于存储所有图
all_plot_list <- list()

# 遍历每个 period2
for (period in unique_periods) {
    # 遍历每个 sample_type
    for (sample in unique_sample_types) {
        # 创建一个辅助列，将 period2 和 sample_type 组合起来
        epi$combined_group <- paste(epi@meta.data$period2, epi@meta.data$sample_type, sep = "_")
        current_group <- paste(period, sample, sep = "_")
        
        # 创建颜色映射：当前 period2 和 sample_type 组合为 npg 红色，其他为灰色
        color_mapping <- setNames(
            ifelse(epi$combined_group == current_group, pal_npg()(1), "gray"),
            epi$combined_group
        )

        # 绘制 DimPlot
        p <- DimPlot(epi, reduction = "umap", label = FALSE, pt.size = 1, group.by = "combined_group") +
            scale_color_manual(values = color_mapping) +
            ggtitle(paste(period, sample, sep = " - ")) +
            theme(
                legend.position = "none",
                plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
                axis.title.x = element_text(size = 20, face = "bold", color = "black"),
                axis.title.y = element_text(size = 20, face = "bold", color = "black"),
                axis.text.x = element_text(size = 16, color = "black"),
                axis.text.y = element_text(size = 16, color = "black")
            )

        # 将图添加到总列表中
        all_plot_list[[paste(period, sample, sep = "_")]] <- p
    }
}

# 使用 patchwork 将图形排列在一起
# 这里根据实际图形数量动态调整每行的图形数量


n_periods <- length(unique_periods)
n_samples <- length(unique_sample_types)

# 创建一个新的列表来存储重新排列后的图形
reordered_plot_list <- list()

# 重新排列顺序
for (i in 1:n_samples) {
    for (j in 1:n_periods) {
        # 计算原始列表中的索引
        original_index <- (j - 1) * n_samples + i
        # 将图形添加到新的列表中
        reordered_plot_list[[length(reordered_plot_list) + 1]] <- all_plot_list[[original_index]]
    }
}

# 使用 wrap_plots 将图形排列成两行四列
combined_plot <- wrap_plots(reordered_plot_list, ncol = 4)

# 保存为 PNG 文件
png("epi_period2_sample_type_umap.png", width = 6000, height = 6000, res = 300)
print(combined_plot)
dev.off()

#patients
#patients



# 假设 identity_mapping 已经定义
cell_type2 <- identity_mapping[epi@meta.data$seurat_clusters]
epi@meta.data$cell_type2 <- cell_type2
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(17)

patients <- unique(epi@meta.data$patients)
patient_numbers <- as.numeric(gsub("[^0-9]", "", patients))
patients_ordered <- patients[order(patient_numbers)]
non_numeric_patients <- patients[is.na(patient_numbers)]  # 提取没有数字部分的患者名称
patients_ordered <- c(patients_ordered, non_numeric_patients)

epi@meta.data$patients <- factor(
  epi@meta.data$patients,
  levels = patients_ordered  # 按数字顺序设置水平
)

npg_extended <- colorRampPalette(npg_pal)(70)
proportion_data <- epi@meta.data %>%
    group_by(cell_type2, patients) %>% summarise(count = n()) %>% mutate(proportion = count / sum(count))


png("epi_celltype_patients_prop1.png", width = 6000, height = 9000, res = 300)  # 设置高分辨率和尺寸
ggplot(proportion_data, aes(x = "", y = proportion, fill = patients)) +
  geom_bar(stat = "identity", width = 1) +         # 堆叠柱状图
  coord_polar(theta = "y") +                       # 转换为饼图
  scale_fill_manual(values = npg_extended) +       # 使用自定义颜色
  theme_void() +                                   # 空白背景
  labs(fill = "Patients") +
  theme(
    legend.position = "right",                     # 图例放在右侧
    plot.title = element_blank(),                  # 移除标题
    # 分面标签设置
    strip.placement = "outside",                   # 标签放在绘图区域外
    strip.text = element_text(                     # 标签样式
      size = 20,                                   # 字体大小
      face = "bold",                               # 加粗
      margin = margin(b = 10)                      # 下方留白（避免与饼图重叠）
    ),
    legend.text = element_text(size = 36),         # 图例文本字体大小
    legend.title = element_text(size = 36)         # 移除图例标题
  ) +
  facet_wrap(
    ~ cell_type2,
    ncol = 3,
    strip.position = "bottom"                      # 标签放在下方
  )
dev.off()


npg_extended <- colorRampPalette(npg_pal)(17)

# 假设 sample_type 有两个水平，例如 "Normal" 和 "Tumor"
# 查看 sample_type 的水平
sample_type_levels <- levels(epi@meta.data$sample_type)

# 假设 sample_type_levels 为 c("Normal", "Tumor")，根据实际数据调整
# 创建颜色映射表，将 "Normal" 和 "Tumor" 的颜色互换
color_mapping <- setNames(npg_extended, sample_type_levels)
# 互换颜色
color_mapping["Normal"] <- npg_extended[2]
color_mapping["Tumor"] <- npg_extended[1]

proportion_data <- epi@meta.data %>%
  group_by(period1, cell_type2, sample_type) %>%
  summarise(count = n()) %>%
  group_by(period1, cell_type2) %>%
  mutate(proportion = count / sum(count))

# 打开 PNG 文件设备
png("epi_celltype_prop1.png", width = 4000, height = 4000, res = 300)

# 绘制饼图
ggplot(proportion_data, aes(x = "", y = proportion, fill = sample_type)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = color_mapping) +  # 使用自定义颜色映射
  labs(fill = "Sample Type") +
  theme_void() +
  theme(
    legend.position = "right",
    plot.title = element_blank(),
    strip.placement = "outside",
    strip.text = element_text(
      size = 28,
      face = "bold",
      margin = margin(b = 10)
    ),
    legend.text = element_text(size = 36),
    legend.title = element_text(size = 36)
  ) +
  facet_grid(
    cell_type2 ~ period1,
    switch = "y",
    labeller = label_value
  )

# 关闭文件设备
dev.off()

npg_extended <- colorRampPalette(npg_pal)(17)

# 假设 sample_type 有两个水平，例如 "Normal" 和 "Tumor"
# 查看 sample_type 的水平
sample_type_levels <- levels(epi@meta.data$sample_type)

# 假设 sample_type_levels 为 c("Normal", "Tumor")，根据实际数据调整
# 创建颜色映射表，将 "Normal" 和 "Tumor" 的颜色互换
color_mapping <- setNames(npg_extended, sample_type_levels)
# 互换颜色
color_mapping["Normal"] <- npg_extended[2]
color_mapping["Tumor"] <- npg_extended[1]

proportion_data <- epi@meta.data %>%
  group_by(period2, cell_type2, sample_type) %>%
  summarise(count = n()) %>%
  group_by(period2, cell_type2) %>%
  mutate(proportion = count / sum(count))

# 打开 PNG 文件设备
png("epi_celltype_prop2.png", width = 4000, height = 4000, res = 300)

# 绘制饼图
ggplot(proportion_data, aes(x = "", y = proportion, fill = sample_type)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = color_mapping) +  # 使用自定义颜色映射
  labs(fill = "Sample Type") +
  theme_void() +
  theme(
    legend.position = "right",
    plot.title = element_blank(),
    strip.placement = "outside",
    strip.text = element_text(
      size = 28,
      face = "bold",
      margin = margin(b = 10)
    ),
    legend.text = element_text(size = 36),
    legend.title = element_text(size = 36)
  ) +
  facet_grid(
    cell_type2 ~ period2,
    switch = "y",
    labeller = label_value
  )

# 关闭文件设备
dev.off()    