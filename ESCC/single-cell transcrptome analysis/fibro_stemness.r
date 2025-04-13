library(Seurat)
library(GSVA)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(ggsci)
fibroblasts <- readRDS("fibro_modify.rds")


stem_pathways <- read.gmt("stem_genesets.v2024.1.Hs.gmt")
stemness_genes <- unique(stem_pathways$gene)

fibroblasts_genes <- rownames(fibroblasts)
matched_genes <- stemness_genes[stemness_genes %in% fibroblasts_genes]
length(matched_genes)
stemness_genes <- matched_genes

fibroblasts <- AddModuleScore(fibroblasts, features = list(stemness_genes), name = "Stemness_Score")
head(fibroblasts@meta.data)
summary(fibroblasts@meta.data$Stemness_Score1)
fibroblasts@meta.data$Stemness_Score_Z <- scale(fibroblasts@meta.data$Stemness_Score1)
summary(fibroblasts@meta.data$Stemness_Score_Z)

npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(4)



p <- FeaturePlot(
  object = fibroblasts,
  features = "Stemness_Score1",
  label = TRUE,
  repel = TRUE,
  order = TRUE,
  pt.size = 1
) +
  scale_color_gradient2(
    low = "#4DBBD5FF",  # 低值颜色
    mid = "white",      # 中值颜色
    high = "#E64B35FF", # 高值颜色
    midpoint = 0,       # 中值点（对应Stemness_Score1的中值）
    limits = c(-0.6, 0.8)  # 固定图例范围
  )
ggsave("fibro_stem_umap.png", plot = p)

npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(15)
p <- VlnPlot(
  object = fibroblasts,
  features = "Stemness_Score_Z",
  group.by = "cell_type",
  pt.size = 0,
  cols = npg_extended
) +
  theme(legend.position = "none",
  axis.title.x = element_blank(), 
  axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("fibro_stem_vln.png", plot = p)


stemness_genes <- c(
    "SOX2",    # 维持干性和多能性的关键转录因子
    "NANOG",   # 多能性标志基因

    "MYC",     # 调控细胞增殖和干性
    "KLF4",    # 多能性相关基因
    "CD44",    # 肿瘤干细胞标志物
    "ALDH1A1", # 醛脱氢酶，与干性相关
    "BMI1",    # 调控干细胞自我更新
    "NOTCH1",  # Notch信号通路，与干性维持相关
    "ESR1",     # 雌激素受体，与某些癌症的干性相关
    "KDM5B",
    "EPAS1",
    "CTNNB1",
    "HIF1A",
    "EZH2"
)


proliferation_genes <- c("CDK1", "CCNB1", "CDC20", "PLK1", "AURKA", "BUB1", "MKI67", "PCNA", "TOP2A", "ESPL1")
fibroblasts_genes <- rownames(fibroblasts)
matched_genes <- proliferation_genes[proliferation_genes %in% fibroblasts_genes]
length(matched_genes)
proliferation_genes <- matched_genes

fibroblasts <- AddModuleScore(fibroblasts, features = list(proliferation_genes), name = "Proliferation_Score")

summary(fibroblasts@meta.data$Proliferation_Score1)
fibroblasts@meta.data$Proliferation_Score_Z <- scale(fibroblasts@meta.data$Proliferation_Score1)
summary(fibroblasts@meta.data$Proliferation_Score_Z)

npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(4)



p <- FeaturePlot(
  object = fibroblasts,
  features = "Proliferation_Score1",
  label = TRUE,
  repel = TRUE,
  order = TRUE,
  pt.size = 1
) +
  scale_color_gradient2(
    low = "#4DBBD5FF",  # 低值颜色
    mid = "white",      # 中值颜色
    high = "#E64B35FF", # 高值颜色
    midpoint = 0       # 中值点（对应Stemness_Score1的中值）

  )
ggsave("fibro_proliferation_umap.png", plot = p)

npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(15)
p <- VlnPlot(
  object = fibroblasts,
  features = "Proliferation_Score_Z",
  group.by = "cell_type",
  pt.size = 0,
  cols = npg_extended
) +
  theme(legend.position = "none",
  axis.title.x = element_blank(), 
  axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("fibro_proliferation_vln.png", plot = p)

#period1
library(patchwork)

# 1. 准备数据
all_sample_types <- as.character(unique(fibroblasts@meta.data$sample_type))
all_cell_types <- as.character(unique(fibroblasts@meta.data$cell_type))

# 2. 创建绘图函数
create_violin <- function(sample, cell) {
    cells_to_keep <- fibroblasts@meta.data$sample_type == sample & 
        fibroblasts@meta.data$cell_type == cell
    
    if (sum(cells_to_keep) == 0) return(NULL)
    
    subset_data <- subset(fibroblasts, cells = colnames(fibroblasts)[cells_to_keep])
    
    npg_pal <- pal_npg()(10)
    unique_periods <- unique(subset_data@meta.data$period1)

    npg_extended <- c("#F8766D","#00BFC4", "#619CFF")
    VlnPlot(
        object = subset_data,
        features = "Proliferation_Score_Z",
        group.by = "period1",
        pt.size = 0,
        cols = npg_extended
    ) +
        theme(
            legend.position = "none",
            axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
            plot.title = element_text(size = 10, face = "bold"),
            plot.margin = unit(c(1, 1, 1, 1), "mm")
        ) +
        ggtitle(paste0(sample, "\n", cell))
}

# 3. 生成所有子图
plot_list <- list()
for (sample in all_sample_types) {
    for (cell in all_cell_types) {
        p <- create_violin(sample, cell)
        if (!is.null(p)) plot_list[[paste(sample, cell, sep = "_")]] <- p
    }
}

# 4. 计算布局行列数（假设每行4个子图）
ncol <- 7
nrow <- ceiling(length(plot_list) / ncol)

# 5. 合并所有子图
combined_plot <- wrap_plots(plot_list, ncol = ncol, nrow = nrow) +
    plot_annotation(
        title = "Proliferation Score Across Conditions",
        theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
    )

# 6. 保存合并后的图（调整尺寸适应子图数量）
ggsave("combined_proliferation_violins.png", 
       plot = combined_plot,
       width = ncol * 3, 
       height = nrow * 3,
       limitsize = FALSE)









# --------------- 显著性计算（基于增殖评分，复用你的卡方检验逻辑） ---------------
all_signif <- list()

for (sample in all_sample_types) {
    for (cell in all_cell_types) {
        # 提取当前样本-细胞的增殖评分数据
        subset_meta <- fibroblasts@meta.data %>%
            filter(sample_type == sample, cell_type == cell)
        
        if (nrow(subset_meta) == 0) next  # 无细胞，跳过
        
        # 获取所有非空的period1（至少1个细胞）
        valid_periods <- subset_meta %>%
            group_by(period1) %>%
            filter(n() > 0) %>%  # 过滤空period1
            pull(period1) %>% unique()
        
        if (length(valid_periods) < 2) next  # 不足2组，无法比较
        
        # 生成period1两两组合（与你的逻辑一致）
        period_pairs <- combn(valid_periods, 2, simplify = FALSE) %>%
            map(sort) %>% unique()
        
        # 对每个period pair进行卡方检验
        for (pair in period_pairs) {
            # 提取两个period的数据
            period1_data <- subset_meta %>% filter(period1 == pair[1])
            period2_data <- subset_meta %>% filter(period1 == pair[2])
            
            # 检查数据有效性（复用你的验证逻辑）
            if (nrow(period1_data) == 0 || nrow(period2_data) == 0) next
            
            # 构建列联表（基于细胞数量）
            cont_table <- matrix(
                c(nrow(period1_data), nrow(period2_data),
                  nrow(subset_meta) - nrow(period1_data),
                  nrow(subset_meta) - nrow(period2_data)),
                nrow = 2
            )
            
            # 卡方检验（复用你的代码逻辑）
            chisq_res <- tryCatch(
                chisq.test(cont_table)$p.value,
                error = function(e) NA
            )
            
            if (is.na(chisq_res)) next
            
            # 记录结果（添加 period1_comparison）
            all_signif <- append(all_signif, list(
                data.frame(
                    cell_type = cell,
                    sample_type = sample,
                    period1_comparison = paste(pair, collapse = " vs "),  # 关键修复
                    p_value = chisq_res,
                    stringsAsFactors = FALSE
                )
            ))
        }
    }
}

# 合并结果并添加显著性标记（与你的逻辑一致）
final_signif <- bind_rows(all_signif) %>%
    mutate(
        significance = case_when(
            p_value < 0.001 ~ "***",
            p_value < 0.01 ~ "**",
            p_value < 0.05 ~ "*",
            TRUE ~ "ns"
        )
    ) %>%
    filter(significance != "ns")  # 仅保留显著结果

# --------------- 绘图函数修改（添加显著性标记） ---------------
create_violin <- function(sample, cell) {
    cells_to_keep <- fibroblasts@meta.data$sample_type == sample & 
        fibroblasts@meta.data$cell_type == cell
    
    if (sum(cells_to_keep) == 0) return(NULL)
    
    subset_data <- subset(fibroblasts, cells = colnames(fibroblasts)[cells_to_keep])
    meta <- subset_data@meta.data  # 提取元数据
    
    # 你的原始绘图代码（完全保留）
    p <- VlnPlot(
        object = subset_data,
        features = "Proliferation_Score_Z",
        group.by = "period1",
        pt.size = 0,
        cols = npg_extended  # 保留你的配色
    ) +
        theme(
            legend.position = "none",
            axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, size = 8),  # 保留你的字体设置
            plot.title = element_text(size = 10, face = "bold"),
            plot.margin = unit(c(1, 1, 1, 1), "mm")  # 保留边距
        ) +
        ggtitle(paste0(sample, "\n", cell))  # 保留标题格式
    
    # 仅添加显著性标记（复用你的代码逻辑）
    sig_data <- final_signif %>%
        filter(sample_type == sample, cell_type == cell) %>%
        separate(period1_comparison, into = c("group1", "group2"), sep = " vs ") %>%  # 关键修复
        arrange(desc(p_value))  # 按显著性排序，避免标记重叠
    
    if (nrow(sig_data) > 0) {
        # 动态计算y位置（基于当前图的最大评分+0.1，与你的历史代码一致）
        y_pos <- max(meta$Proliferation_Score_Z, na.rm = TRUE) + 0.1
        
        # 绘制显著性标记（完全复用你的geom_signif参数）
        p <- p + geom_signif(
            data = sig_data,
            aes(xmin = group1, xmax = group2, annotations = significance),
            comparisons = sig_data %>% select(group1, group2) %>% as.list(),
            y_position = y_pos,
            textsize = 3,  # 保留你的字体大小
            tip_length = 0.01,  # 保留线尖长度
            vjust = 0.5,  # 垂直居中
            manual = TRUE  # 手动指定位置
        )
    }
    
    return(p)
}

# --------------- 后续代码与你的原始逻辑完全一致 ---------------
# 生成所有子图（循环结构不变）
plot_list <- list()
for (sample in all_sample_types) {
    for (cell in all_cell_types) {
        p <- create_violin(sample, cell)
        if (!is.null(p)) plot_list[[paste(sample, cell, sep = "_")]] <- p
    }
}

# 合并和保存（参数完全不变）
ncol <- 7
nrow <- ceiling(length(plot_list) / ncol)
combined_plot <- wrap_plots(plot_list, ncol = ncol, nrow = nrow) +
    plot_annotation(title = "Proliferation Score Across Conditions")
ggsave("combined_proliferation_violins.png", combined_plot, 
       width = ncol * 3, height = nrow * 3, limitsize = FALSE)