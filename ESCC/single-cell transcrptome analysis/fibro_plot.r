#sample_type
# 获取唯一的 sample_type 值，并按字母顺序排序
unique_sample_types <- sort(unique(fibroblasts@meta.data$sample_type))  # 按字母顺序排序

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
    p <- DimPlot(fibroblasts, reduction = "umap", label = FALSE, pt.size = 1, group.by = "sample_type") +
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

# 使用 patchwork 将图形排列在一起（每行四个）
combined_plot <- wrap_plots(plot_list, ncol = 2)  # 每行四个图

# 保存为 PNG 文件
png("fibro_sample_type_umap.png", width = 6000, height = 3000, res = 300)
print(combined_plot)
dev.off()



#period2
# 获取唯一的 period2 值，并按字母顺序排序
unique_periods <- sort(unique(fibroblasts@meta.data$period2))  # 按字母顺序排序

# 如果需要手动指定顺序，可以这样做：
# unique_periods <- c("period_A", "period_B", "period_C", "period_D", "period_E", "period_F")

# 创建一个空列表，用于存储每个 period2 的图
plot_list <- list()

# 遍历每个 period2，生成单独的图
for (period in unique_periods) {
    # 创建颜色映射：当前 period2 为 npg 红色，其他为灰色
    color_mapping <- setNames(
        ifelse(unique_periods == period, pal_npg()(1), "gray"),  # npg 红色
        unique_periods
    )
    
    # 绘制 DimPlot
    p <- DimPlot(fibroblasts, reduction = "umap", label = FALSE, pt.size = 1, group.by = "period2") +
        scale_color_manual(values = color_mapping) +  # 使用自定义颜色映射
        ggtitle(period) +  # 设置标题为当前 period2
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

# 使用 patchwork 将图形排列在一起（每行四个）
combined_plot <- wrap_plots(plot_list, ncol = 4)  # 每行四个图

# 保存为 PNG 文件
png("fibro_period2_umap.png", width = 6000, height = 3000, res = 300)
print(combined_plot)
dev.off()


#period1
# 获取唯一的 period1 值，并按字母顺序排序
unique_periods <- sort(unique(fibroblasts@meta.data$period1))  # 按字母顺序排序

# 如果需要手动指定顺序，可以这样做：
# unique_periods <- c("period_A", "period_B", "period_C", "period_D", "period_E", "period_F")

# 创建一个空列表，用于存储每个 period1 的图
plot_list <- list()

# 遍历每个 period1，生成单独的图
for (period in unique_periods) {
    # 创建颜色映射：当前 period1 为 npg 红色，其他为灰色
    color_mapping <- setNames(
        ifelse(unique_periods == period, pal_npg()(1), "gray"),  # npg 红色
        unique_periods
    )
    
    # 绘制 DimPlot
    p <- DimPlot(fibroblasts, reduction = "umap", label = FALSE, pt.size = 1, group.by = "period1") +
        scale_color_manual(values = color_mapping) +  # 使用自定义颜色映射
        ggtitle(period) +  # 设置标题为当前 period1
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

# 使用 patchwork 将图形排列在一起（每行四个）
combined_plot <- wrap_plots(plot_list, ncol = 3)  # 每行四个图

# 保存为 PNG 文件
png("fibro_period1_umap.png", width = 9000, height = 3000, res = 300)
print(combined_plot)
dev.off()

#Ⅰ
library(Seurat)
library(patchwork)
library(ggsci)

# 获取唯一的 period1 值，并按字母顺序排序
unique_periods <- sort(unique(fibroblasts@meta.data$period1))

# 获取唯一的 sample_type 值，并按字母顺序排序
unique_sample_types <- sort(unique(fibroblasts@meta.data$sample_type))

# 创建一个空列表，用于存储所有图
all_plot_list <- list()

# 遍历每个 period1
for (period in unique_periods) {
    # 遍历每个 sample_type
    for (sample in unique_sample_types) {
        # 创建一个辅助列，将 period1 和 sample_type 组合起来
        fibroblasts$combined_group <- paste(fibroblasts@meta.data$period1, fibroblasts@meta.data$sample_type, sep = "_")
        current_group <- paste(period, sample, sep = "_")
        
        # 创建颜色映射：当前 period1 和 sample_type 组合为 npg 红色，其他为灰色
        color_mapping <- setNames(
            ifelse(fibroblasts$combined_group == current_group, pal_npg()(1), "gray"),
            fibroblasts$combined_group
        )

        # 绘制 DimPlot
        p <- DimPlot(fibroblasts, reduction = "umap", label = FALSE, pt.size = 1, group.by = "combined_group") +
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
png("fibro_period1_sample_type_umap.png", width = 6000, height = 3000, res = 300)
print(combined_plot)
dev.off()

#
# 获取唯一的 period2 值，并按字母顺序排序
unique_periods <- sort(unique(fibroblasts@meta.data$period2))

# 获取唯一的 sample_type 值，并按字母顺序排序
unique_sample_types <- sort(unique(fibroblasts@meta.data$sample_type))

# 创建一个空列表，用于存储所有图
all_plot_list <- list()

# 遍历每个 period2
for (period in unique_periods) {
    # 遍历每个 sample_type
    for (sample in unique_sample_types) {
        # 创建一个辅助列，将 period2 和 sample_type 组合起来
        fibroblasts$combined_group <- paste(fibroblasts@meta.data$period2, fibroblasts@meta.data$sample_type, sep = "_")
        current_group <- paste(period, sample, sep = "_")
        
        # 创建颜色映射：当前 period2 和 sample_type 组合为 npg 红色，其他为灰色
        color_mapping <- setNames(
            ifelse(fibroblasts$combined_group == current_group, pal_npg()(1), "gray"),
            fibroblasts$combined_group
        )

        # 绘制 DimPlot
        p <- DimPlot(fibroblasts, reduction = "umap", label = FALSE, pt.size = 1, group.by = "combined_group") +
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
combined_plot <- wrap_plots(all_plot_list, ncol = 4)

# 保存为 PNG 文件
png("fibro_period2_sample_type_umap.png", width = 6000, height = 6000, res = 300)
print(combined_plot)
dev.off()


npg_extended <- colorRampPalette(npg_pal)(14)

# 假设 sample_type 有两个水平，例如 "Normal" 和 "Tumor"
# 查看 sample_type 的水平
sample_type_levels <- levels(fibroblasts@meta.data$sample_type)

# 假设 sample_type_levels 为 c("Normal", "Tumor")，根据实际数据调整
# 创建颜色映射表，将 "Normal" 和 "Tumor" 的颜色互换
color_mapping <- setNames(npg_extended, sample_type_levels)
# 互换颜色
color_mapping["Normal"] <- npg_extended[2]
color_mapping["Tumor"] <- npg_extended[1]

proportion_data <- fibroblasts@meta.data %>%
  group_by(period1, cell_type, sample_type) %>%
  summarise(count = n()) %>%
  group_by(period1, cell_type) %>%
  mutate(proportion = count / sum(count))

# 打开 PNG 文件设备
png("fibro_celltype_prop1.png", width = 4000, height = 4000, res = 300)

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
    cell_type ~ period1,
    switch = "y",
    labeller = label_value
  )

# 关闭文件设备
dev.off()

npg_extended <- colorRampPalette(npg_pal)(14)

# 假设 sample_type 有两个水平，例如 "Normal" 和 "Tumor"
# 查看 sample_type 的水平
sample_type_levels <- levels(fibroblasts@meta.data$sample_type)

# 假设 sample_type_levels 为 c("Normal", "Tumor")，根据实际数据调整
# 创建颜色映射表，将 "Normal" 和 "Tumor" 的颜色互换
color_mapping <- setNames(npg_extended, sample_type_levels)
# 互换颜色
color_mapping["Normal"] <- npg_extended[2]
color_mapping["Tumor"] <- npg_extended[1]

proportion_data <- fibroblasts@meta.data %>%
  group_by(period2, cell_type, sample_type) %>%
  summarise(count = n()) %>%
  group_by(period2, cell_type) %>%
  mutate(proportion = count / sum(count))

# 打开 PNG 文件设备
png("fibro_celltype_prop2.png", width = 6000, height = 4000, res = 300)

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
    cell_type ~ period2,
    switch = "y",
    labeller = label_value
  )

# 关闭文件设备
dev.off()

#patients
patients <- unique(fibroblasts@meta.data$patients)
patient_numbers <- as.numeric(gsub("[^0-9]", "", patients))
patients_ordered <- patients[order(patient_numbers)]
non_numeric_patients <- patients[is.na(patient_numbers)]  # 提取没有数字部分的患者名称
patients_ordered <- c(patients_ordered, non_numeric_patients)

fibroblasts@meta.data$patients <- factor(
  fibroblasts@meta.data$patients,
  levels = patients_ordered  # 按数字顺序设置水平
)

npg_extended <- colorRampPalette(npg_pal)(70)
proportion_data <- fibroblasts@meta.data %>%
    group_by(cell_type, patients) %>% summarise(count = n()) %>% mutate(proportion = count / sum(count))


png("fibro_celltype_patients_prop1.png", width = 6000, height = 6000, res = 300)  # 设置高分辨率和尺寸
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
    legend.title = element_text(size = 36)                 # 移除图例标题
  ) +
  facet_wrap(
    ~ cell_type,
    ncol = 3,
    strip.position = "bottom"                      # 标签放在下方
  )
dev.off()



library(ggpubr)

cell_counts <- fibroblasts@meta.data %>%  # 将 merged_seurat_obj 替换为 fibroblasts
  group_by(patients, cell_type, sample_type) %>%
  summarise(count = n()) %>%
  ungroup()

# 计算每个样本的总细胞数
sample_totals <- cell_counts %>%
  group_by(patients) %>%
  summarise(total = sum(count)) %>%
  ungroup()

# 将样本总细胞数合并到 cell_counts 中
cell_counts <- cell_counts %>%
  left_join(sample_totals, by = "patients")

# 计算细胞类型比例
cell_proportions <- fibroblasts@meta.data %>%  # 将 merged_seurat_obj 替换为 fibroblasts
  group_by(patients, cell_type, sample_type) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  left_join(sample_totals, by = "patients") %>%
  mutate(proportion = count / total)

# 对每个细胞类型构建列联表并进行卡方检验
chisq_results <- cell_proportions %>%
  group_by(cell_type) %>%
  summarise(
    p_value = {
      normal_count <- sum(count[sample_type == "Normal"])
      tumor_count <- sum(count[sample_type == "Tumor"])
      normal_total <- sum(total[sample_type == "Normal"])
      tumor_total <- sum(total[sample_type == "Tumor"])
      
      cont_table <- matrix(
        c(normal_count, tumor_count, normal_total - normal_count, tumor_total - tumor_count),
        nrow = 2
      )
      
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

# 绘制基于比例的箱线图并添加显著性标记
png("箱图_fibro_proportion.png", width = 6000, height = 3000, res = 300)
ggplot(cell_proportions, aes(x = cell_type, y = proportion, fill = sample_type)) +
  geom_boxplot() +
  geom_text(
    data = chisq_results, 
    aes(x = cell_type, y = max(cell_proportions$proportion) * 1.05, label = significance),
    size = 12, 
    vjust = 0.5, 
    inherit.aes = FALSE 
  ) +
  labs(x = "", y = "Cell Proportion", fill = "Sample Type") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 28),
    axis.text.y = element_text(size = 28),
    axis.title.y = element_text(size = 40),
    axis.line = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 36),
    legend.title = element_text(size = 40),
    legend.position = "right" 
  ) +
  scale_fill_npg() +
  scale_y_continuous(limits = c(0, 1))
dev.off()

library(ggpubr)
cell_counts <- fibroblasts@meta.data %>%
  group_by(sample, period1, sample_type, cell_type) %>%  # ✅ 新增sample分组
  summarise(count = n()) %>%
  ungroup()
total_counts <- cell_counts %>%
  group_by(sample, period1, sample_type) %>%
  summarise(total_count = sum(count)) %>%
  ungroup()

# 合并数据并计算比例
cell_proportions <- cell_counts %>%
  left_join(total_counts, by = c("sample", "period1", "sample_type")) %>%
  mutate(proportion = count / total_count)


library(tidyr)
cell_proportions <- cell_proportions %>%
    mutate(
        cell_type = gsub("\\+", "_", cell_type),
        cell_type = ifelse(cell_type == "", "Unknown", cell_type)  # 直接填充为 "Unknown"
    )



# 对每个 cell_type 构建列联表并进行卡方检验
all_sample_types <- unique(cell_proportions$sample_type)
all_periods <- unique(cell_proportions$period1)

# 生成所有 period1 的两两组合（确保顺序正确）
period_pairs <- combn(all_periods, 2, simplify = FALSE) %>%
  map(~ sort(.x)) %>%  # 确保组合顺序一致（如 I vs II，而不是 II vs I）
  unique()  # 去重

# 初始化结果列表
all_results <- list()

for (sample_type in all_sample_types) {
  for (pair in period_pairs) {
    # 提取当前 sample_type 和 period pair 的数据
    subset_data <- cell_proportions %>%
      filter(sample_type == !!sample_type, period1 %in% pair)
    
    # 检查每个 cell_type 是否在两个 period1 中都有数据
    valid_cell_types <- subset_data %>%
      group_by(cell_type) %>%
      summarise(
        n_periods = n_distinct(period1),
        has_valid_counts = all(!is.na(count)) & all(!is.na(total_count))
      ) %>%
      filter(n_periods == 2 & has_valid_counts) %>%
      pull(cell_type)
    
    # 只对有效的 cell_type 进行计算
    if (length(valid_cell_types) > 0) {
      chisq_results <- subset_data %>%
        filter(cell_type %in% valid_cell_types) %>%
        group_by(cell_type) %>%
        summarise(
          p_value = {
            # 提取两个时期的细胞数量和总细胞数
            period1_data <- filter(cur_data(), period1 == pair[1])
            period2_data <- filter(cur_data(), period1 == pair[2])
            
            # 如果数据缺失，返回 NA
            if (nrow(period1_data) == 0 || nrow(period2_data) == 0) return(NA)
            
            period1_count <- sum(period1_data$count, na.rm = TRUE)
            period2_count <- sum(period2_data$count, na.rm = TRUE)
            period1_total <- sum(period1_data$total_count, na.rm = TRUE)
            period2_total <- sum(period2_data$total_count, na.rm = TRUE)
            
            # 如果计数为 0，返回 NA
            if (period1_count == 0 || period2_count == 0 || 
                period1_total == 0 || period2_total == 0) return(NA)
            
            # 构建 2x2 列联表
            cont_table <- matrix(
              c(period1_count, period2_count,
                period1_total - period1_count,
                period2_total - period2_count),
              nrow = 2
            )
            
            # 进行卡方检验（如果失败则返回 NA）
            tryCatch(
              chisq.test(cont_table)$p.value,
              error = function(e) NA
            )
          },
          period1_comparison = paste(pair, collapse = " vs "),
          sample_type = sample_type
        ) %>%
        ungroup() %>%
        filter(!is.na(p_value))  # 移除无效计算结果
      
      # 添加显著性标记
      if (nrow(chisq_results) > 0) {
        chisq_results <- chisq_results %>%
          mutate(
            significance = case_when(
              p_value < 0.001 ~ "***",
              p_value < 0.01 ~ "**",
              p_value < 0.05 ~ "*",
              TRUE ~ "ns"
            )
          )
        
        # 将结果添加到列表中
        all_results <- append(all_results, list(chisq_results))
      }
    }
  }
}

# 合并所有结果
final_results <- bind_rows(all_results)

# 提取每个细胞类型在不同 period1_comparison 之下的 significance 的第一条记录
# （前面的代码保持不变，直到 first_significance 部分）

# 提取每个细胞类型在不同 period1_comparison 和 sample_type 之下的第一条记录
first_significance <- final_results %>%
    group_by(cell_type, period1_comparison, sample_type) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    select(cell_type, period1_comparison, significance, sample_type)

# （后面的绘图代码修改如下）

plots <- list()

for (cell_type in unique(cell_proportions$cell_type)) {
    for (sample_type in unique(cell_proportions$sample_type)) {
        # 筛选出当前 cell_type 和 sample_type 的数据
        subset_data <- cell_proportions %>%
            filter(cell_type == !!cell_type & sample_type == !!sample_type)
        
        # 获取当前 cell_type 和 sample_type 的显著性结果（关键修改：添加 sample_type 筛选）
        sig_data <- first_significance %>%
            filter(cell_type == !!cell_type & sample_type == !!sample_type) %>%
            separate(period1_comparison, into = c("group1", "group2"), sep = " vs ")
        
        # 创建显著性标记的列表
        comparisons <- list()
        annotations <- list()
        
        if (nrow(sig_data) > 0) {
            for (i in 1:nrow(sig_data)) {
                comparisons[[i]] <- c(sig_data$group1[i], sig_data$group2[i])
                annotations[[i]] <- sig_data$significance[i]
            }
        }
        
        # 创建箱线图（保持不变）
        p <- ggplot(subset_data, aes(x = period1, y = proportion, fill = period1)) +
            geom_boxplot() +
            labs(x = "", y = "", fill = "Period") +
            theme_classic() +
            theme(
                axis.text.x = element_text(angle = 0, hjust = 1, size = 36),
                axis.text.y = element_text(size = 28),
                axis.title.y = element_blank(),
                axis.line = element_line(color = "black"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),

                legend.position = "none",
                plot.title = element_text(size = 28)
            ) +
            scale_fill_manual(values = c("I" = "#F8766D", "II" = "#00BFC4", "III" = "#619CFF")) +
            scale_y_continuous(limits = c(0, 1)) +
            ggtitle(paste(cell_type, sample_type))
        
        # 添加显著性标记（仅当前 sample_type 的比较）
        if (length(comparisons) > 0) {
            for (i in 1:length(comparisons)) {
                p <- p + geom_signif(
                    comparisons = list(comparisons[[i]]),
                    annotations = annotations[[i]],
                    y_position = 0.9 - (i-1)*0.1,  # 调整垂直位置
                    textsize = 8,
                    vjust = 0.5,
                    tip_length = 0.01
                )
            }
        }
        
        plots <- append(plots, list(p))
    }
}

# （合并和保存代码保持不变）

# 计算所需的行数，确保所有图都能合理显示
n_plots <- length(plots)
ncol <- 4
nrow <- ceiling(n_plots / ncol)

# 使用 patchwork 合并所有图
combined_plot <- wrap_plots(plots, ncol = ncol, nrow = nrow) + 
  plot_layout(guides = 'collect')

# 保存合并后的图
ggsave("fibro_combined_boxplots_with_signif.png", combined_plot, 
       width = 24, height = 6 * nrow, units = "in", dpi = 300)