#significance
count_table <- fibroblasts@meta.data %>%
    dplyr::filter(sample_type %in% c("Tumor", "Normal")) %>%
    dplyr::count(cell_type, sample_type) %>%  # 直接计算各组合的细胞数
    tidyr::pivot_wider(
        names_from = sample_type,
        values_from = n,
        values_fill = 0  # 缺失组合补0（如某细胞类型仅存在于Tumor）
    ) %>%
    tibble::column_to_rownames("cell_type") %>%
    as.matrix()  # 转换为检验所需的矩阵格式
get_significance <- function(count_table) {
    results <- NULL
    for (i in 1:nrow(count_table)) {
        a <- count_table[i, "Tumor"]  # 当前细胞Tumor计数
        b <- count_table[i, "Normal"]# 当前细胞Normal计数
        c <- sum(count_table[-i, "Tumor"]) # 其他细胞Tumor合计
        d <- sum(count_table[-i, "Normal"])# 其他细胞Normal合计
        
        mat <- matrix(c(a, c, b, d), 2, dimnames=list(c("CellType", "Other"), c("Tumor", "Normal")))
        method <- if (min(mat) < 5) "fisher" else "chisq"
        p_val <- switch(method,
                        fisher = fisher.test(mat)$p.value,
                        chisq = chisq.test(mat)$p.value
        )
        
        results <- rbind(results, data.frame(
            cell_type = rownames(count_table)[i],
            p_value = p_val,
            method = method,
            stringsAsFactors = FALSE
        ))
    }
    return(results)
}

# 2. 运行检验+FDR校正
signif_results <- get_significance(count_table) %>%
    dplyr::mutate(
        adj_p = p.adjust(p_value, "fdr"),
        significant = adj_p < 0.05
    ) %>%
    dplyr::arrange(adj_p) %>%
    dplyr::select(cell_type, method, p_value, adj_p, significant)

#plot
library(ggpubr)
cell_counts <- fibroblasts@meta.data %>%
  group_by(orig.ident, cell_type, sample_type) %>%
  summarise(count = n()) %>%
  ungroup()

# 计算每个样本的总细胞数
sample_totals <- cell_counts %>%
  group_by(orig.ident) %>%
  summarise(total = sum(count)) %>%
  ungroup()

# 将样本总细胞数合并到 cell_counts 中
cell_counts <- cell_counts %>%
  left_join(sample_totals, by = "orig.ident")

# 计算细胞类型比例
cell_proportions <- fibroblasts@meta.data %>%
  group_by(orig.ident, cell_type, sample_type) %>%
  summarise(count = n()) %>%  # 计算每种细胞类型的数量
  ungroup() %>%
  left_join(sample_totals, by = "orig.ident") %>%  # 合并样本总细胞数
  mutate(proportion = count / total)  # 计算细胞类型比例


# 添加显著性标记
chisq_results <- signif_results %>%
  mutate(
    significance = case_when(
      adj_p < 0.001 ~ "***",
      adj_p < 0.01 ~ "**",
      adj_p < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )
color_mapping <- c("Tumor" = "#F8766D", "Normal" = "#00BFC4")  
# 绘制基于比例的箱线图并添加显著性标记
png("fibro_箱图_proportion.png", width = 6000, height = 3000, res = 300)
ggplot(cell_proportions, aes(x = cell_type, y = proportion, fill = sample_type)) +
  geom_boxplot() +  # 绘制箱线图
  geom_text(
    data = chisq_results, 
    aes(x = cell_type, y = max(cell_proportions$proportion) * 1.05, label = significance),  # 调整 y 值
    size = 12, 
    vjust = 0.5,  # 调整 vjust 参数
    inherit.aes = FALSE  # 忽略父图层的 aes 映射
  ) +  # 添加显著性标记
  labs(x = "", y = "Cell Proportion", fill = "Sample Type") +  # 设置坐标轴和图例标题
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
  scale_fill_manual(values = color_mapping)
dev.off()



#period1
period1 <- subset(fibroblasts, subset = period1 == "I")
# 1. 数据聚合（保留所有细胞类型，含零值）
count_table <- period1@meta.data %>%
  dplyr::filter(sample_type %in% c("Tumor", "Normal")) %>%
  dplyr::count(cell_type, sample_type) %>%
  tidyr::pivot_wider(
    names_from = sample_type,
    values_from = n,
    values_fill = 0  # 保留零值
  ) %>%
  tibble::column_to_rownames("cell_type") %>%
  as.matrix()

# 2. 零值检验标记（记空而非过滤）
get_significance <- function(count_table) {
  results <- NULL
  
  for (i in 1:nrow(count_table)) {
    ct <- rownames(count_table)[i]
    a <- count_table[i, "Tumor"]
    b <- count_table[i, "Normal"]
    
    # 零值标记逻辑（记空不检验）
    if (a == 0 || b == 0) {
      results <- rbind(results, data.frame(
        cell_type = ct,
        p_value = NA,
        method = "ZeroExcluded",
        a = a, b = b,
        stringsAsFactors = FALSE
      ))
      next  # 跳过检验
    }
    
    # 非零值检验（保留原逻辑）
    c <- sum(count_table[-i, "Tumor"])
    d <- sum(count_table[-i, "Normal"])
    mat <- matrix(c(a, c, b, d), 2)
    method <- if (min(mat) < 5) "fisher" else "chisq"
    p_val <- switch(method,
      fisher = fisher.test(mat)$p.value,
      chisq = chisq.test(mat)$p.value
    )
    
    results <- rbind(results, data.frame(
      cell_type = ct,
      p_value = p_val,
      method = method,
      a = a, b = b,
      stringsAsFactors = FALSE
    ))
  }
  return(results)
}

# 3. 结果整理（空值记为""）
# 修复后完整代码（关键修改已标注）
signif_results <- get_significance(count_table) %>%
    dplyr::mutate(
        adj_p = ifelse(is.na(p_value), NA, p.adjust(p_value, "fdr")),
        significance = case_when(
            is.na(p_value) ~ "",
            adj_p < 0.001 ~ "***",
            adj_p < 0.01 ~ "**",
            adj_p < 0.05 ~ "*",
            TRUE ~ "ns"
        )
    ) %>%
    dplyr::arrange(cell_type)

cell_counts <- period1@meta.data %>%
  group_by(orig.ident, cell_type, sample_type) %>%
  summarise(count = n()) %>%
  ungroup()

# 计算每个样本的总细胞数
sample_totals <- cell_counts %>%
  group_by(orig.ident) %>%
  summarise(total = sum(count)) %>%
  ungroup()

# 将样本总细胞数合并到 cell_counts 中
cell_counts <- cell_counts %>%
  left_join(sample_totals, by = "orig.ident")

# 计算细胞类型比例
cell_proportions <- period1@meta.data %>%
  group_by(orig.ident, cell_type, sample_type) %>%
  summarise(count = n()) %>%  # 计算每种细胞类型的数量
  ungroup() %>%
  left_join(sample_totals, by = "orig.ident") %>%  # 合并样本总细胞数
  mutate(proportion = count / total)  # 计算细胞类型比例


# 核心修复：确保绘图数据含proportion列（新增校验）
if (!"proportion" %in% colnames(cell_proportions)) {
    stop("错误：cell_proportions缺少'proportion'列，请检查数据聚合步骤")
}


png("fibro_period1_箱图_proportion.png", width = 6000, height = 3000, res = 300)
ggplot(cell_proportions, aes(x = cell_type, y = proportion, fill = sample_type)) +
    geom_boxplot() +
    geom_text(
        data = signif_results %>% filter(significance != ""),
        # 关键修复：使用绘图数据的最大值（来自cell_proportions）
        aes(x = cell_type, y = max(cell_proportions$proportion, na.rm = TRUE) * 1.05, label = significance),
        size = 12, vjust = 0.5, inherit.aes = FALSE,
        color = "#222222"  # 深灰色标记
    ) +
    # 新增：零值细胞底部注释（可选）
    geom_text(
        data = cell_proportions %>% filter(proportion == 0),  # 直接基于绘图数据
        aes(x = cell_type, y = 0, label = "0", vjust = 1.5),
        size = 10, color = "#666666", inherit.aes = FALSE
    ) +
    labs(title = "Period I", x = "", y = "Cell Proportion", fill = "Sample Type") +
    theme_classic() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 28),
        axis.text.y = element_text(size = 28),
        axis.title = element_text(size = 28),
        legend.text = element_text(size = 36),        
        legend.title = element_text(size = 40)
    ) +
    scale_fill_manual(values = color_mapping)
dev.off()