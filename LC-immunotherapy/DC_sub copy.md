

# 饼图
```
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(16)
library(dplyr)

LUAD <- subset(myeloid, subset = cancer_type == "LUAD")
LUSC <- subset(myeloid, subset = cancer_type == "LUSC")


DC_data <- LUSC@meta.data %>%
  filter(major_cell == "DC") %>%
  group_by(LN_group, major_cell, cell_type) %>%
  summarise(count = n()) %>%
  group_by(LN_group, major_cell) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup() %>%
  mutate(group_major = paste0(LN_group, "_", major_cell)) %>%
  mutate(group_major = factor(group_major,
                             levels = c("LN+_DC", "LN-DC")))


pdf("LUSC_Group_DC_cell_type_pie.pdf", width = 4000/300, height = 3000/300) 

ggplot(DC_data, aes(x = "", y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack", width = 1) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = npg_extended) +
  facet_wrap(~ group_major, ncol = 2) +  # 改为2列
  labs(x = "", y = "", fill = "Cell Type", 
       title = "DC Cell Type Composition") +
  theme_void() +
  theme(
    legend.text = element_text(size = 28),
    legend.title = element_text(size = 28),
    strip.text = element_text(
      size = 20,
      face = "bold",
      margin = margin(b = 15)
    ),
    legend.position = "right",
    plot.title = element_text(size = 32, face = "bold", hjust = 0.5, margin = margin(b = 20))
  )

dev.off()
#tiff
ggsave(
  filename = "LUSC_Group_DC_cell_type_pie.tiff",
  plot = last_plot(),
  device = "tiff",
  width = 4000/300,  # 单位：英寸
  height = 3000/300,
  units = "in",
  dpi = 300,
  compression = "lzw",
  bg = "white"
)




#LUAD
DC_data <- LUAD@meta.data %>%
  filter(major_cell == "DC") %>% 
  group_by(LN_group, major_cell, cell_type) %>%
  summarise(count = n()) %>%
  group_by(LN_group, major_cell) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup() %>%
  mutate(group_major = paste0(LN_group, "_", major_cell)) %>%
  mutate(group_major = factor(group_major,
                             levels = c("LN+_DC", "LN-_DC")))


pdf("LUAD_Group_DC_cell_type_pie.pdf", width = 4000/300, height = 3000/300) 

ggplot(DC_data, aes(x = "", y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack", width = 1) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = npg_extended) +
  facet_wrap(~ group_major, ncol = 2) +  # 改为2列
  labs(x = "", y = "", fill = "Cell Type", 
       title = "DC Cell Type Composition") +
  theme_void() +
  theme(
    legend.text = element_text(size = 28),
    legend.title = element_text(size = 28),
    strip.text = element_text(
      size = 20,
      face = "bold",
      margin = margin(b = 15)
    ),
    legend.position = "right",
    plot.title = element_text(size = 32, face = "bold", hjust = 0.5, margin = margin(b = 20))
  )

dev.off()


ggsave(
  filename = "LUAD_Group_DC_cell_type_pie.tiff",
  plot = last_plot(),
  device = "tiff",
  width = 4000/300,  # 单位：英寸
  height = 3000/300,
  units = "in",
  dpi = 300,
  compression = "lzw",
  bg = "white"
)


```
# 箱图(LUSC)
```
library(ggpubr)
library(dplyr)
library(ggsci)
library(tidyr)
LUSC <- subset(myeloid,subset=cancer_type=="LUSC")

macro <- subset(LUSC,subset=major_cell=="Macro")
DC <- subset(LUSC,subset=major_cell=="DC")

# 统计每个样本-细胞类型-LN分组的细胞数
cell_counts <- DC@meta.data %>%
  group_by(patients, cell_type, LN_group) %>%  # Group改为LN_group
  summarise(count = n(), .groups = "drop")

# 计算每个样本的总细胞数
sample_totals <- cell_counts %>%
  group_by(patients) %>%
  summarise(total = sum(count), .groups = "drop")

# 计算细胞类型比例
cell_proportions <- cell_counts %>%
  left_join(sample_totals, by = "patients") %>%
  mutate(proportion = count / total)

# === 找出哪个细胞类型在LN-中缺失 ===
# 检查每个细胞类型在两个LN_group中的存在情况
celltype_group_presence <- cell_proportions %>%
  group_by(cell_type, LN_group) %>%  # Group改为LN_group
  summarise(has_data = n() > 0, .groups = "drop") %>%
  pivot_wider(names_from = LN_group, values_from = has_data, values_fill = FALSE)

# 找出在LN+中存在但在LN-中缺失的细胞类型
missing_in_ln_neg <- celltype_group_presence %>%
  filter(`LN+` == TRUE & `LN-` == FALSE) %>%
  pull(cell_type)

cat("在LN-中缺失的细胞类型:", paste(missing_in_ln_neg, collapse = ", "), "\n")

# === 只为缺失的细胞类型在LN-中添加0值 ===
if (length(missing_in_ln_neg) > 0) {
  # 获取LN-的所有患者
  ln_neg_patients <- unique(cell_proportions$patients[cell_proportions$LN_group == "LN-"])
  
  # 只为缺失的细胞类型创建0值记录
  missing_records <- expand_grid(
    patients = ln_neg_patients,
    cell_type = missing_in_ln_neg,
    LN_group = "LN-"
  ) %>%
  left_join(sample_totals, by = "patients") %>%
  mutate(count = 0, proportion = 0)
  
  # 合并到绘图数据
  plot_data <- cell_proportions %>% bind_rows(missing_records)
} else {
  plot_data <- cell_proportions
}

# === 关键修复：确保分组数量匹配 ===
# 获取所有细胞类型和LN_group
all_celltypes <- unique(plot_data$cell_type)
all_ln_groups <- c("LN+", "LN-")  # 明确指定两个分组
n_groups <- length(all_ln_groups)

# 构建完整的细胞类型-LN_group组合并填充0
full_combinations <- expand.grid(
  cell_type = all_celltypes,
  LN_group = all_ln_groups,
  stringsAsFactors = FALSE
)

# 按细胞类型和LN_group汇总计数
cell_group_counts <- cell_counts %>%
  group_by(cell_type, LN_group) %>%  # Group改为LN_group
  summarise(total_count = sum(count), .groups = "drop") %>%
  right_join(full_combinations, by = c("cell_type", "LN_group")) %>%  # Group改为LN_group
  mutate(total_count = replace(total_count, is.na(total_count), 0)) %>%
  arrange(cell_type, LN_group)  # 确保分组顺序一致

# 计算每个LN_group的总细胞数（全局）
global_group_totals <- cell_counts %>%
  group_by(LN_group) %>%  # Group改为LN_group
  summarise(group_total = sum(count), .groups = "drop") %>%
  right_join(data.frame(LN_group = all_ln_groups), by = "LN_group") %>%  # Group改为LN_group
  mutate(group_total = replace(group_total, is.na(group_total), 0)) %>%
  arrange(LN_group)  # 保持分组顺序一致

# === 修复卡方检验逻辑 ===
chisq_results <- cell_group_counts %>%
  group_by(cell_type) %>%
  summarise(
    p_value = {
      # 获取当前细胞类型在两个LN_group的计数
      obs <- c(
        total_count[LN_group == "LN+"],
        total_count[LN_group == "LN-"]
      )
      
      # 获取对应LN_group的总细胞数
      group_tot <- c(
        global_group_totals$group_total[global_group_totals$LN_group == "LN+"],
        global_group_totals$group_total[global_group_totals$LN_group == "LN-"]
      )
      
      # 检查数据有效性
      if (sum(obs) < 5 || length(obs) != 2) {
        NA_real_
      } else {
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

# 设置LN_group的顺序
plot_data$LN_group <- factor(plot_data$LN_group, levels = c("LN+", "LN-"))

# === 绘制图形 ===
y_max <- max(plot_data$proportion, na.rm = TRUE) * 1.1

pdf("LUSC_DC_box_proportion_by_ln_group.pdf", width = 6000/300, height = 3000/300)
ggplot(plot_data, aes(x = cell_type, y = proportion, fill = LN_group)) +  # Group改为LN_group
  geom_boxplot(width = 0.7) +
  geom_text(
    data = chisq_results, 
    aes(x = cell_type, y = y_max, label = significance),
    size = 12, 
    inherit.aes = FALSE
  ) +
  labs(x = "", y = "Cell Proportion", fill = "LN Group") +  # Group改为LN Group
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

p <- ggplot(plot_data, aes(x = cell_type, y = proportion, fill = LN_group)) +
  geom_boxplot(width = 0.7) +
  geom_text(
    data = chisq_results, 
    aes(x = cell_type, y = y_max, label = significance),
    size = 12, 
    inherit.aes = FALSE
  ) +
  labs(x = "", y = "Cell Proportion", fill = "LN Group") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.y = element_text(size = 32),
    axis.line = element_line(color = "black"),
    panel.grid = element_blank(),
    legend.text = element_text(size = 36),
    legend.title = element_text(size = 40),
    legend.position = "right",
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  scale_fill_npg() +
  scale_y_continuous(limits = c(0, y_max))

# 使用ggsave保存为TIFF
ggsave("LUSC_DC_box_proportion_by_ln_group.tiff",
       plot = p,
       width = 6000/300,  # 20英寸
       height = 3000/300, # 10英寸
       units = "in",
       dpi = 300,
       device = "tiff",
       compression = "lzw",
       bg = "white")
```
# DC二合一（LUSC）
```
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(patchwork)

# ===== 第一部分：准备柱状图数据 =====
prepare_bar_plot_data <- function(DC_obj) {
  cell_proportions_bar <- DC_obj@meta.data %>%
    group_by(LN_group, pathological_response, cell_type) %>%
    summarise(count = n(), .groups = 'drop') %>%
    group_by(LN_group, pathological_response) %>%
    mutate(proportion = count / sum(count)) %>%
    ungroup()
  
  # 设置LN组顺序：LN+在前，LN-在后
  cell_proportions_bar$LN_group <- factor(
    cell_proportions_bar$LN_group,
    levels = c("LN+", "LN-")
  )
  
  # 设置病理反应顺序
  cell_proportions_bar$pathological_response <- factor(
    cell_proportions_bar$pathological_response, 
    levels = c("pCR", "MPR", "non-MPR")
  )
  
  all_cell_types <- unique(DC_obj@meta.data$cell_type)
  cell_proportions_complete <- cell_proportions_bar %>%
    complete(LN_group, pathological_response, cell_type = all_cell_types, 
             fill = list(count = 0, proportion = 0))
  
  # 保持LN组和病理反应顺序
  cell_proportions_complete$LN_group <- factor(
    cell_proportions_complete$LN_group,
    levels = c("LN+", "LN-")
  )
  
  cell_proportions_complete$pathological_response <- factor(
    cell_proportions_complete$pathological_response, 
    levels = c("pCR", "MPR", "non-MPR")
  )
  
  celltype_alphabetical <- sort(all_cell_types)
  cell_proportions_complete$cell_type <- factor(
    cell_proportions_complete$cell_type, 
    levels = celltype_alphabetical
  )
  
  return(list(
    bar_data = cell_proportions_complete,
    celltypes = celltype_alphabetical
  ))
}

# ===== 第二部分：数据处理模块 =====

# 模块1: 提取和计算细胞计数
extract_cell_counts <- function(DC_obj, ln_group) {
  cell_counts <- DC_obj@meta.data %>%
    filter(LN_group == ln_group) %>%
    group_by(patients, cell_type, pathological_response) %>%
    summarise(count = n(), .groups = "drop")
  
  return(cell_counts)
}

# 模块2: 计算样本总数
calculate_sample_totals <- function(cell_counts) {
  sample_totals <- cell_counts %>%
    group_by(patients) %>%
    summarise(total = sum(count), .groups = "drop")
  
  return(sample_totals)
}

# 模块3: 计算细胞比例
calculate_cell_proportions <- function(cell_counts, sample_totals) {
  cell_proportions <- cell_counts %>%
    left_join(sample_totals, by = "patients") %>%
    mutate(proportion = count / total)
  
  # 设置病理反应顺序
  cell_proportions$pathological_response <- factor(
    cell_proportions$pathological_response, 
    levels = c("pCR", "MPR", "non-MPR")
  )
  
  return(cell_proportions)
}

# 模块4: 补充缺失数据
supplement_missing_data <- function(cell_proportions, sample_totals) {
  # 找出缺失的病理反应并添加0值
  celltype_response_presence <- cell_proportions %>%
    group_by(cell_type, pathological_response) %>%
    summarise(has_data = n() > 0, .groups = "drop") %>%
    pivot_wider(names_from = pathological_response, values_from = has_data, values_fill = FALSE)
  
  # 检查缺失的组合
  all_responses <- c("pCR", "MPR", "non-MPR")
  missing_combinations <- list()
  
  for (resp in all_responses) {
    missing_celltypes <- celltype_response_presence %>%
      filter(.data[[resp]] == FALSE) %>%
      pull(cell_type)
    
    if (length(missing_celltypes) > 0) {
      missing_combinations[[resp]] <- missing_celltypes
    }
  }
  
  plot_data <- cell_proportions
  
  # 为缺失的组合添加0值
  for (resp in names(missing_combinations)) {
    resp_patients <- unique(cell_proportions$patients[cell_proportions$pathological_response == resp])
    
    if (length(resp_patients) > 0) {
      missing_records <- expand_grid(
        patients = resp_patients,
        cell_type = missing_combinations[[resp]],
        pathological_response = resp
      ) %>%
      left_join(sample_totals, by = "patients") %>%
      mutate(count = 0, proportion = 0)
      
      plot_data <- plot_data %>% bind_rows(missing_records)
    }
  }
  
  return(plot_data)
}

# 模块5: 执行统计检验
perform_statistical_tests <- function(cell_counts, plot_data) {
  # 统计检验 - 三组比较
  all_celltypes <- unique(plot_data$cell_type)
  full_combinations <- expand_grid(
    cell_type = all_celltypes,
    pathological_response = factor(c("pCR", "MPR", "non-MPR"), 
                                  levels = c("pCR", "MPR", "non-MPR"))
  )
  
  cell_response_counts <- cell_counts %>%
    group_by(cell_type, pathological_response) %>%
    summarise(total_count = sum(count), .groups = "drop") %>%
    right_join(full_combinations, by = c("cell_type", "pathological_response")) %>%
    mutate(total_count = replace(total_count, is.na(total_count), 0))
  
  # 计算全局病理反应组总数
  global_response_totals <- cell_counts %>%
    group_by(pathological_response) %>%
    summarise(response_total = sum(count), .groups = "drop") %>%
    right_join(data.frame(pathological_response = factor(c("pCR", "MPR", "non-MPR"), 
                                                        levels = c("pCR", "MPR", "non-MPR"))), 
              by = "pathological_response") %>%
    mutate(response_total = replace(response_total, is.na(response_total), 0))
  
  # 卡方检验（三组比较）
  chisq_results <- cell_response_counts %>%
    group_by(cell_type) %>%
    summarise(
      p_value = calculate_p_value_for_celltype(cur_data(), global_response_totals),
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
  
  return(chisq_results)
}

# 模块6: 为单个细胞类型计算p值
calculate_p_value_for_celltype <- function(cell_data, global_response_totals) {
  # 获取三组的计数
  pcr_count <- cell_data$total_count[cell_data$pathological_response == "pCR"]
  mpr_count <- cell_data$total_count[cell_data$pathological_response == "MPR"]
  non_mpr_count <- cell_data$total_count[cell_data$pathological_response == "non-MPR"]
  
  pcr_total <- global_response_totals$response_total[global_response_totals$pathological_response == "pCR"]
  mpr_total <- global_response_totals$response_total[global_response_totals$pathological_response == "MPR"]
  non_mpr_total <- global_response_totals$response_total[global_response_totals$pathological_response == "non-MPR"]
  
  total_counts <- sum(c(pcr_count, mpr_count, non_mpr_count))
  
  if (total_counts < 15 || any(c(pcr_count, mpr_count, non_mpr_count) < 5)) {
    return(NA_real_)
  }
  
  # 构建3x2列联表
  other_pcr <- pcr_total - pcr_count
  other_mpr <- mpr_total - mpr_count
  other_non_mpr <- non_mpr_total - non_mpr_count
  
  cont_table <- matrix(c(pcr_count, other_pcr,
                        mpr_count, other_mpr,
                        non_mpr_count, other_non_mpr),
                      nrow = 3, ncol = 2, byrow = TRUE)
  
  if (any(cont_table < 0) || sum(cont_table) == 0) {
    return(NA_real_)
  } else if (any(cont_table < 5)) {
    return(fisher.test(cont_table)$p.value)
  } else {
    return(chisq.test(cont_table)$p.value)
  }
}

# 模块7: 处理单个LN组数据
process_ln_group_data <- function(DC_obj, ln_group) {
  # 提取细胞计数
  cell_counts <- extract_cell_counts(DC_obj, ln_group)
  
  # 计算样本总数
  sample_totals <- calculate_sample_totals(cell_counts)
  
  # 计算细胞比例
  cell_proportions <- calculate_cell_proportions(cell_counts, sample_totals)
  
  # 补充缺失数据
  plot_data <- supplement_missing_data(cell_proportions, sample_totals)
  
  # 执行统计检验
  chisq_results <- perform_statistical_tests(cell_counts, plot_data)
  
  return(list(
    plot_data = plot_data,
    chisq_results = chisq_results
  ))
}

# ===== 第三部分：主数据处理流程（确保LN+在前） =====
process_box_plot_data <- function(DC_obj, celltype_alphabetical) {
  # 按顺序处理：先LN+，后LN-
  ln_groups_ordered <- c("LN+", "LN-")
  
  # 初始化列表存储结果
  plot_data_list <- list()
  stats_list <- list()
  
  for (ln_group in ln_groups_ordered) {
    group_data <- process_ln_group_data(DC_obj, ln_group)
    
    plot_data_list[[ln_group]] <- group_data$plot_data %>% mutate(LN_group = ln_group)
    stats_list[[ln_group]] <- group_data$chisq_results %>% mutate(LN_group = ln_group)
  }
  
  # 合并数据，保持LN+在前
  plot_data_box_combined <- bind_rows(plot_data_list)
  chisq_results_combined <- bind_rows(stats_list)
  
  # 设置LN组因子顺序：LN+在前
  plot_data_box_combined$LN_group <- factor(
    plot_data_box_combined$LN_group,
    levels = c("LN+", "LN-")
  )
  
  chisq_results_combined$LN_group <- factor(
    chisq_results_combined$LN_group,
    levels = c("LN+", "LN-")
  )
  
  # 确保数据按病理反应顺序排序
  plot_data_box_combined <- plot_data_box_combined %>%
    arrange(LN_group, cell_type, pathological_response) %>%
    mutate(
      pathological_response = factor(pathological_response, 
                                    levels = c("pCR", "MPR", "non-MPR")),
      order_factor = interaction(cell_type, pathological_response, sep = "_")
    )
  
  # 设置细胞类型顺序
  plot_data_box_combined$cell_type <- factor(plot_data_box_combined$cell_type, 
                                            levels = celltype_alphabetical)
  chisq_results_combined$cell_type <- factor(chisq_results_combined$cell_type, 
                                            levels = celltype_alphabetical)
  
  # 计算y轴最大值
  y_max <- max(plot_data_box_combined$proportion, na.rm = TRUE) * 1.1
  
  return(list(
    box_data = plot_data_box_combined,
    stats_results = chisq_results_combined,
    y_max = y_max
  ))
}

# ===== 第四部分：绘图模块 =====

# 模块8: 创建柱状图
create_bar_plot <- function(bar_data, celltypes) {
  # 设置颜色
  response_colors <- c("pCR" = "#00A087FF", "MPR" = "#F39B7FFF", "non-MPR" = "#E64B35FF")
  
  p <- ggplot(bar_data, aes(x = cell_type, y = proportion, fill = pathological_response)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    scale_fill_manual(
      values = response_colors,
      breaks = c("pCR", "MPR", "non-MPR"),
      labels = c("pCR", "MPR", "non-MPR")
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
    facet_grid(LN_group ~ ., scales = "free_x")
  
  return(p)
}

# 模块9: 创建箱图
create_box_plot <- function(box_data, stats_results, y_max, celltypes) {
  # 设置颜色
  response_colors <- c("pCR" = "#00A087FF", "MPR" = "#F39B7FFF", "non-MPR" = "#E64B35FF")
  
  p <- ggplot(box_data, aes(x = cell_type, y = proportion, fill = pathological_response)) +
    geom_boxplot(width = 0.7, aes(group = interaction(cell_type, pathological_response))) +
    # 为0值数据添加横线
    geom_segment(
      data = box_data %>% filter(proportion == 0),
      aes(x = as.numeric(cell_type) - 0.25, 
          xend = as.numeric(cell_type) + 0.25,
          y = 0, yend = 0),
      color = "black", linewidth = 1, inherit.aes = FALSE
    ) +
    # 添加显著性标记
    geom_text(
      data = stats_results, 
      aes(x = cell_type, y = y_max, label = significance),
      size = 6, 
      inherit.aes = FALSE
    ) +
    scale_fill_manual(
      values = response_colors,
      breaks = c("pCR", "MPR", "non-MPR"),
      labels = c("pCR", "MPR", "non-MPR")
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
    facet_grid(LN_group ~ ., scales = "free_x")
  
  return(p)
}

# 模块10: 创建图例
create_legend <- function(bar_data) {
  response_colors <- c("pCR" = "#00A087FF", "MPR" = "#F39B7FFF", "non-MPR" = "#E64B35FF")
  
  legend_plot <- ggplot(bar_data, 
                        aes(x = cell_type, y = proportion, fill = pathological_response)) +
    geom_col(position = position_dodge(width = 0.8)) +
    scale_fill_manual(
      values = response_colors,
      breaks = c("pCR", "MPR", "non-MPR"),
      labels = c("pCR", "MPR", "non-MPR"),
      name = "Pathological Response"
    ) +
    theme(legend.position = "top",
          legend.text = element_text(size = 14, face = "bold"),
          legend.title = element_text(size = 16, face = "bold"))
  
  legend <- get_legend(legend_plot)
  return(legend)
}

# ===== 第五部分：主程序 =====

# 1. 准备柱状图数据
bar_data_result <- prepare_bar_plot_data(DC)
bar_data <- bar_data_result$bar_data
celltypes <- bar_data_result$celltypes

# 2. 处理箱图数据（确保LN+在前）
box_data_result <- process_box_plot_data(DC, celltypes)
box_data <- box_data_result$box_data
stats_results <- box_data_result$stats_results
y_max <- box_data_result$y_max

# 3. 创建图形
bar_plot <- create_bar_plot(bar_data, celltypes)
box_plot <- create_box_plot(box_data, stats_results, y_max, celltypes)
legend <- create_legend(bar_data)

# 4. 组合图形
combined_plot <- (bar_plot | box_plot) +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")

# 添加图例
final_plot <- wrap_plots(legend, combined_plot, ncol = 1, heights = c(0.05, 0.95))

# 5. 保存图形
pdf("LUSC_DC_combined_bar_box_plots.pdf", width = 6000/300, height = 3000/300)
print(final_plot)
dev.off()

# 同时保存TIFF格式
ggsave("LUSC_DC_combined_bar_box_plots.tiff",
       plot = final_plot,
       width = 6000/300,  # 20英寸 (6000/300 = 20)
       height = 3000/300, # 10英寸 (3000/300 = 10)
       units = "in",
       dpi = 300,
       device = "tiff",
       compression = "lzw",
       bg = "white")

cat("完成！已生成PDF和TIFF格式的图形。\n")
cat("图形顺序：LN+在上方，LN-在下方\n")
cat("输出文件：\n")
cat("  - DC_combined_bar_box_plots.pdf\n")
cat("  - DC_combined_bar_box_plots.tiff\n")
```
# 二合一（LUAD）