# DC二合一
```
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(patchwork)

# ===== 第一部分：准备柱状图数据 =====
cell_proportions_bar <- DC@meta.data %>%
  group_by(Group, Response, cell_type) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(Group, Response) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()

# 设置Response的顺序为YES, NO
cell_proportions_bar$Response <- factor(cell_proportions_bar$Response, levels = c("YES", "NO"))

all_cell_types <- unique(DC@meta.data$cell_type)
cell_proportions_complete <- cell_proportions_bar %>%
  complete(Group, Response, cell_type = all_cell_types, fill = list(count = 0, proportion = 0))

# 保持Response的顺序
cell_proportions_complete$Response <- factor(cell_proportions_complete$Response, levels = c("YES", "NO"))

celltype_alphabetical <- sort(all_cell_types)
cell_proportions_complete$cell_type <- factor(cell_proportions_complete$cell_type, levels = celltype_alphabetical)

# ===== 第二部分：准备箱图数据（包含显著性检验） =====
process_group_data <- function(current_group) {
  cell_counts <- DC@meta.data %>%
    filter(Group == current_group) %>%
    group_by(patients, cell_type, Response) %>%
    summarise(count = n(), .groups = "drop")
  
  sample_totals <- cell_counts %>%
    group_by(patients) %>%
    summarise(total = sum(count), .groups = "drop")
  
  cell_proportions <- cell_counts %>%
    left_join(sample_totals, by = "patients") %>%
    mutate(proportion = count / total)
  
  # 设置Response的顺序为YES, NO
  cell_proportions$Response <- factor(cell_proportions$Response, levels = c("YES", "NO"))
  
  # 找出缺失的细胞类型并添加0值（用于创建横线占位）
  celltype_response_presence <- cell_proportions %>%
    group_by(cell_type, Response) %>%
    summarise(has_data = n() > 0, .groups = "drop") %>%
    pivot_wider(names_from = Response, values_from = has_data, values_fill = FALSE)
  
  missing_in_NO <- celltype_response_presence %>%
    filter(YES == TRUE & NO == FALSE) %>%
    pull(cell_type)
  
  missing_in_YES <- celltype_response_presence %>%
    filter(YES == FALSE & NO == TRUE) %>%
    pull(cell_type)
  
  plot_data <- cell_proportions
  
  # 为缺失的细胞类型添加0值（用于创建横线占位）
  if (length(missing_in_NO) > 0) {
    NO_patients <- unique(cell_proportions$patients[cell_proportions$Response == "NO"])
    missing_records_NO <- expand_grid(
      patients = NO_patients,
      cell_type = missing_in_NO,
      Response = "NO"
    ) %>%
    left_join(sample_totals, by = "patients") %>%
    mutate(count = 0, proportion = 0)
    
    plot_data <- plot_data %>% bind_rows(missing_records_NO)
  }
  
  if (length(missing_in_YES) > 0) {
    YES_patients <- unique(cell_proportions$patients[cell_proportions$Response == "YES"])
    missing_records_YES <- expand_grid(
      patients = YES_patients,
      cell_type = missing_in_YES,
      Response = "YES"
    ) %>%
    left_join(sample_totals, by = "patients") %>%
    mutate(count = 0, proportion = 0)
    
    plot_data <- plot_data %>% bind_rows(missing_records_YES)
  }
  
  # 统计检验
  all_celltypes <- unique(plot_data$cell_type)
  full_combinations <- expand_grid(
    cell_type = all_celltypes,
    Response = factor(c("YES", "NO"), levels = c("YES", "NO"))
  )
  
  cell_response_counts <- cell_counts %>%
    group_by(cell_type, Response) %>%
    summarise(total_count = sum(count), .groups = "drop") %>%
    right_join(full_combinations, by = c("cell_type", "Response")) %>%
    mutate(total_count = replace(total_count, is.na(total_count), 0)) %>%
    arrange(cell_type, Response)
  
  global_response_totals <- cell_counts %>%
    group_by(Response) %>%
    summarise(response_total = sum(count), .groups = "drop") %>%
    right_join(data.frame(Response = factor(c("YES", "NO"), levels = c("YES", "NO"))), by = "Response") %>%
    mutate(response_total = replace(response_total, is.na(response_total), 0)) %>%
    arrange(Response)
  
  chisq_results <- cell_response_counts %>%
    group_by(cell_type) %>%
    summarise(
      p_value = {
        yes_count <- total_count[Response == "YES"]
        no_count <- total_count[Response == "NO"]
        yes_total <- global_response_totals$response_total[global_response_totals$Response == "YES"]
        no_total <- global_response_totals$response_total[global_response_totals$Response == "NO"]
        
        if (sum(c(yes_count, no_count)) < 5) {
          NA_real_
        } else {
          cont_table <- matrix(c(yes_count, yes_total - yes_count,
                                 no_count, no_total - no_count),
                               nrow = 2, byrow = FALSE)
          
          if (any(cont_table < 0) || sum(cont_table) == 0) {
            NA_real_
          } else if (any(cont_table < 5)) {
            fisher.test(cont_table)$p.value
          } else {
            chisq.test(cont_table)$p.value
          }
        }
      },
      .groups = "drop"
    ) %>%
    mutate(
      significance = case_when(
        is.na(p_value) ~ "ns",
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**", 
        p_value < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    )
  
  return(list(plot_data = plot_data, chisq_results = chisq_results))
}

# 处理两组数据
group1_data <- process_group_data("Group1")
group2_data <- process_group_data("Group2")

# 合并箱图数据
plot_data_box_combined <- bind_rows(
  group1_data$plot_data %>% mutate(Group = "Group1"),
  group2_data$plot_data %>% mutate(Group = "Group2")
)

chisq_results_combined <- bind_rows(
  group1_data$chisq_results %>% mutate(Group = "Group1"),
  group2_data$chisq_results %>% mutate(Group = "Group2")
)

# 设置细胞类型顺序
plot_data_box_combined$cell_type <- factor(plot_data_box_combined$cell_type, levels = celltype_alphabetical)
chisq_results_combined$cell_type <- factor(chisq_results_combined$cell_type, levels = celltype_alphabetical)

# 计算y轴最大值
y_max <- max(plot_data_box_combined$proportion, na.rm = TRUE) * 1.1

# ===== 第三部分：创建图形 =====
# 创建柱状图（左侧）- 两行布局
bar_plot <- ggplot(cell_proportions_complete, aes(x = cell_type, y = proportion, fill = Response)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(
    values = c("YES" = "#E64B35FF", "NO" = "#4DBBD5FF"), 
    breaks = c("YES", "NO")  # 确保图例顺序
  ) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "", y = "Proportion", title = "Bar Plot") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black", face = "bold"),
    axis.text.y = element_text(size = 12, color = "black", face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold", color = "black"),
    axis.line = element_line(color = "black"),
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  ) +
  facet_grid(Group ~ ., scales = "free_x")

# 创建箱图（右侧）- 两行布局，添加显著性标记
box_plot <- ggplot(plot_data_box_combined, aes(x = cell_type, y = proportion, fill = Response)) +
  geom_boxplot(width = 0.7) +
  # 为0值数据添加横线
  geom_segment(
    data = plot_data_box_combined %>% filter(proportion == 0),
    aes(x = as.numeric(cell_type) - 0.35, xend = as.numeric(cell_type) + 0.35,
        y = 0, yend = 0),
    color = "black", linewidth = 1, inherit.aes = FALSE
  ) +
  # 添加显著性标记
  geom_text(
    data = chisq_results_combined, 
    aes(x = cell_type, y = y_max, label = significance),
    size = 6, 
    inherit.aes = FALSE
  ) +
  scale_fill_manual(
    values = c("YES" = "#E64B35FF", "NO" = "#4DBBD5FF"),
    breaks = c("YES", "NO")  # 确保图例顺序
  ) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, y_max)) +
  labs(x = "Cell Type", y = "Cell Proportion", title = "Box Plot") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black", face = "bold"),
    axis.text.y = element_text(size = 12, color = "black", face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold", color = "black"),
    axis.line = element_line(color = "black"),
    panel.grid = element_blank(),
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  ) +
  facet_grid(Group ~ ., scales = "free_x")

# ===== 第四部分：合并图形 =====
# 创建图例
legend_plot <- ggplot(cell_proportions_complete, aes(x = cell_type, y = proportion, fill = Response)) +
  geom_col(position = position_dodge(width = 0.8)) +
  scale_fill_manual(
    values = c("YES" = "#E64B35FF", "NO" = "#4DBBD5FF"),
    breaks = c("YES", "NO")  # 确保图例顺序
  ) +
  theme(legend.position = "top",
        legend.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"))

legend <- get_legend(legend_plot)

# 组合图形
combined_plot <- (bar_plot | box_plot) +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")

# 添加图例
final_plot <- wrap_plots(legend, combined_plot, ncol = 1, heights = c(0.05, 0.95))

# 保存图形
pdf("DC_combined_bar_box_plots.pdf", width = 6000/300, height = 3000/300)
print(final_plot)
dev.off()

cat("已完成合并图形\n")
```
# DC 
```
library(Seurat)
library(ggplot2)
library(ggsci)
library(dplyr)
plot_data <- DC


npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(14)

pdf("DC_combined_annotation_by_group.pdf", width = 4000/300, height = 2000/300)

p <- DimPlot(plot_data, 
             reduction = "umap", 
             label = TRUE, 
             pt.size = 1, 
             group.by = "cell_type", 
             label.size = 4,
             split.by = "Group",  
             ncol = 2) + 
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


#箱图
library(ggpubr)
library(dplyr)
library(ggsci)
library(tidyr)

# 统计每个样本-细胞类型-分组的细胞数
cell_counts <- DC@meta.data %>%
  group_by(patients, cell_type, Group) %>%
  summarise(count = n(), .groups = "drop")

# 计算每个样本的总细胞数
sample_totals <- cell_counts %>%
  group_by(patients) %>%
  summarise(total = sum(count), .groups = "drop")

# 计算细胞类型比例
cell_proportions <- cell_counts %>%
  left_join(sample_totals, by = "patients") %>%
  mutate(proportion = count / total)

# === 找出哪个细胞类型在Group2中缺失 ===
# 检查每个细胞类型在两个Group中的存在情况
celltype_group_presence <- cell_proportions %>%
  group_by(cell_type, Group) %>%
  summarise(has_data = n() > 0, .groups = "drop") %>%
  pivot_wider(names_from = Group, values_from = has_data, values_fill = FALSE)

# 找出在Group1中存在但在Group2中缺失的细胞类型
missing_in_group2 <- celltype_group_presence %>%
  filter(Group1 == TRUE & Group2 == FALSE) %>%
  pull(cell_type)

cat("在Group2中缺失的细胞类型:", paste(missing_in_group2, collapse = ", "), "\n")

# === 只为缺失的细胞类型在Group2中添加0值 ===
if (length(missing_in_group2) > 0) {
  # 获取Group2的所有患者
  group2_patients <- unique(cell_proportions$patients[cell_proportions$Group == "Group2"])
  
  # 只为缺失的细胞类型创建0值记录
  missing_records <- expand_grid(
    patients = group2_patients,
    cell_type = missing_in_group2,
    Group = "Group2"
  ) %>%
  left_join(sample_totals, by = "patients") %>%
  mutate(count = 0, proportion = 0)
  
  # 合并到绘图数据
  plot_data <- cell_proportions %>% bind_rows(missing_records)
} else {
  plot_data <- cell_proportions
}

# === 关键修复：确保分组数量匹配 ===
# 获取所有细胞类型和分组
all_celltypes <- unique(plot_data$cell_type)
all_groups <- unique(plot_data$Group)
n_groups <- length(all_groups)

# 构建完整的细胞类型-分组组合并填充0
full_combinations <- expand.grid(
  cell_type = all_celltypes,
  Group = all_groups,
  stringsAsFactors = FALSE
)

# 按细胞类型和分组汇总计数
cell_group_counts <- cell_counts %>%
  group_by(cell_type, Group) %>%
  summarise(total_count = sum(count), .groups = "drop") %>%
  right_join(full_combinations, by = c("cell_type", "Group")) %>%
  mutate(total_count = replace(total_count, is.na(total_count), 0)) %>%
  arrange(cell_type, Group)  # 确保分组顺序一致

# 计算每个分组的总细胞数（全局）
global_group_totals <- cell_counts %>%
  group_by(Group) %>%
  summarise(group_total = sum(count), .groups = "drop") %>%
  right_join(data.frame(Group = all_groups), by = "Group") %>%
  mutate(group_total = replace(group_total, is.na(group_total), 0)) %>%
  arrange(Group)  # 保持分组顺序一致

# === 修复卡方检验逻辑 ===
chisq_results <- cell_group_counts %>%
  group_by(cell_type) %>%
  summarise(
    p_value = {
      # 获取当前细胞类型在所有分组的计数（顺序与global_group_totals一致）
      obs <- total_count[match(all_groups, Group)]
      
      # 获取对应分组的总细胞数
      group_tot <- global_group_totals$group_total
      
      # 检查数据有效性
      if (sum(obs) < 5 || length(obs) != length(group_tot)) {
        NA_real_
      } else {
        # 计算期望比例（基于全局分组总细胞数）
        expected_prop <- group_tot / sum(group_tot)
        
        # 处理两组情况（Fisher精确检验优先）
        if (length(obs) == 2) {
          # 构建2x2列联表（细胞类型计数 vs 其他细胞计数）
          other_counts <- group_tot - obs
          cont_table <- matrix(c(obs[1], other_counts[1],
                                 obs[2], other_counts[2]),
                               nrow = 2, byrow = TRUE)
          
          # 检查列联表有效性
          if (any(cont_table < 0) || sum(cont_table) == 0) {
            NA_real_
          } else if (any(cont_table < 5)) {
            fisher.test(cont_table)$p.value
          } else {
            chisq.test(cont_table)$p.value
          }
        } else {
          # 多组情况使用卡方检验
          if (any(obs * expected_prop < 5)) {
            NA_real_  # 期望频数不足时返回NA
          } else {
            chisq.test(obs, p = expected_prop)$p.value
          }
        }
      }
    },
    .groups = "drop"
  ) %>%
  # 添加显著性标记
  mutate(
    significance = case_when(
      is.na(p_value) ~ "ns",
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**", 
      p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

# === 按首字母排序 ===
celltype_alphabetical <- sort(all_celltypes)
plot_data$cell_type <- factor(plot_data$cell_type, levels = celltype_alphabetical)
chisq_results$cell_type <- factor(chisq_results$cell_type, levels = celltype_alphabetical)

# === 绘制图形 ===
y_max <- max(plot_data$proportion, na.rm = TRUE) * 1.1

pdf("DC_box_proportion_by_group.pdf", width = 6000/300, height = 3000/300)
ggplot(plot_data, aes(x = cell_type, y = proportion, fill = Group)) +
  geom_boxplot(width = 0.7) +
  geom_text(
    data = chisq_results, 
    aes(x = cell_type, y = y_max, label = significance),
    size = 12, 
    inherit.aes = FALSE
  ) +
  labs(x = "", y = "Cell Proportion", fill = "Group") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.y = element_text(size = 32),
    axis.line = element_line(color = "black"),
    panel.grid = element_blank(),
    legend.text = element_text(size = 36),
    legend.title = element_text(size = 40),
    legend.position = "right"
  ) +
  scale_fill_npg() +
  scale_y_continuous(limits = c(0, y_max))
dev.off()
```