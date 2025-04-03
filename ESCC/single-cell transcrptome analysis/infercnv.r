epi_T <- readRDS("epi_T.rds")
#expr_matrix <- GetAssayData(epi_T, layer = "counts")
T <- subset(epi_T, subset = seurat_clusters %in% c(0,1))
epi <- subset(epi_T, subset = seurat_clusters %in% c(2,9,13))
# 获取表达矩阵

T_expr_matrix <- GetAssayData(T, layer = "counts")
epi_expr_matrix <- GetAssayData(epi, layer = "counts")
T_matrix <- T_expr_matrix  # 所有 T 细胞
Epi_matrix <- epi_expr_matrix 
# 找到共同的基因
common_genes <- intersect(rownames(T_expr_matrix), rownames(epi_expr_matrix))

# 过滤表达矩阵，只保留共同基因
T_expr_matrix <- T_expr_matrix[common_genes, ]
epi_expr_matrix <- epi_expr_matrix[common_genes, ]

# 合并表达矩阵
expression_matrix <- cbind(T_expr_matrix, epi_expr_matrix)

epi_clusters <- epi$seurat_clusters
# 生成注释文件
annotations_file <- data.frame(
    cells = colnames(expression_matrix),  # 细胞名称
    cell_type = ifelse(colnames(expression_matrix) %in% colnames(T_expr_matrix), "T cells", 
                      paste0("Epi_Cluster_", epi_clusters[match(colnames(expression_matrix), colnames(epi_expr_matrix))]))  # 细胞类型
)

# 加载必要的库
library(infercnv)
library(AnnoProbe)
library(biomaRt)

# 获取基因名称
gene_names <- rownames(expression_matrix)

# 基因注释
gene_annotation <- annoGene(
    IDs = gene_names,  # 基因名称
    ID_type = "SYMBOL",  # 基因名称类型（如 SYMBOL 或 ENSEMBL）
    species = "human"  # 物种
)

# 生成基因位置参考文件
gene_location_ref <- data.frame(
    gene = gene_annotation$SYMBOL,  # 基因符号
    chr = gene_annotation$chr,  # 染色体名称
    start = gene_annotation$start,  # 基因起始位置
    stop = gene_annotation$end
)

# 保存基因位置参考文件
write.table(gene_location_ref, file = "gene_location_ref.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# 读取基因位置参考文件
gene_location_ref <- read.table("gene_location_ref.txt", header = FALSE, sep = "\t")

# 检查重复的基因名称
duplicated_genes <- gene_location_ref$V1[duplicated(gene_location_ref$V1)]
print(duplicated_genes)
# 去除重复的基因名称
gene_location_ref <- gene_location_ref[!duplicated(gene_location_ref$V1), ]

# 保存为文件
write.table(gene_location_ref, file = "gene_location_ref.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
groupFiles <- 'groupFiles.txt'
write.table(annotations_file, file = groupFiles, sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE)

# 3. 保存基因位置参考文件到文件
geneFile <- 'geneFile.txt'
write.table(gene_location_ref, file = geneFile, sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE)

# 4. 创建 infercnv 对象
infercnv_obj <- CreateInfercnvObject(
    raw_counts_matrix = expression_matrix,  # 表达矩阵文件路径
    annotations_file = groupFiles,  # 注释文件路径
    delim = "\t",  # 文件分隔符
    gene_order_file = geneFile,  # 基因位置参考文件路径
    ref_group_names = c("T cells")  # 以 T 细胞作为参考
)

infercnv_obj <- infercnv::run(
    infercnv_obj,  # infercnv 对象
    cutoff = 0.2,  # 设置 cutoff 值（默认值为 0.1）
    out_dir = "./infercnv",  # 输出目录
    cluster_by_groups = TRUE,  # 按组聚类
    denoise = TRUE,  # 去噪
    HMM = FALSE  # 使用 HMM 模型
)



#抽样
set.seed(123)

# 设置抽样数量（例如抽取 10,000 个细胞）
n_samples <- 10000
sampled_epi_cells <- sample(colnames(epi_expr_matrix), size = n_samples, replace = FALSE)
sampled_expression_matrix <- epi_expr_matrix[, sampled_epi_cells]
sampled_annotations_file <- data.frame(
    cells = sampled_epi_cells,
    cell_type = paste0("Epi_Cluster_", epi$seurat_clusters[match(sampled_epi_cells, colnames(epi))])
)
groupFiles <- "sampled_groupFiles.txt"
write.table(sampled_annotations_file, file = groupFiles, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

geneFile <- "geneFile.txt"
write.table(gene_location_ref, file = geneFile, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

T_ref_matrix <- T_expr_matrix[, ]

infercnv_obj <- CreateInfercnvObject(
    raw_counts_matrix = sampled_expression_matrix,
    annotations_file = groupFiles,
    delim = "\t",
    gene_order_file = geneFile,
    reference_group_matrix = T_ref_matrix  # 参考组仍为 T 细胞（需单独准备 T 细胞矩阵）
)

sampled_cells <- sample(colnames(expression_matrix), size = n_samples, replace = FALSE)

sampled_expression_matrix <- expression_matrix[, sampled_cells]
sampled_annotations_file <- annotations_file[annotations_file$cells %in% sampled_cells, ]

groupFiles <- "sampled_groupFiles.txt"
write.table(sampled_annotations_file, file = groupFiles, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

geneFile <- "geneFile.txt"
write.table(gene_location_ref, file = geneFile, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

infercnv_obj <- CreateInfercnvObject(
    raw_counts_matrix = sampled_expression_matrix,  # 直接使用表达矩阵
    annotations_file = groupFiles,  # 注释文件路径
    delim = "\t",  # 文件分隔符
    gene_order_file = geneFile,  # 基因位置参考文件路径
    ref_group_names = c("T cells")  # 以 T 细胞作为参考
)

infercnv_obj <- infercnv::run(
    infercnv_obj,
    cutoff = 0.2,  # 设置 cutoff 值
    out_dir = "infercnv_output",  # 输出目录
    cluster_by_groups = TRUE,  # 按组聚类
    denoise = TRUE,  # 去噪
    HMM = FALSE  # 是否使用 HMM 模型
)

#cnv score
epi$is_sampled <- colnames(epi) %in% sampled_cells
epi_sampled <- subset(epi, cells = sampled_cells)

infercnv_obj <- readRDS("infercnv_output/run.final.infercnv_obj")
cnv_matrix <- infercnv_obj@expr.data 
cnv_cells <- colnames(cnv_matrix)
sampled_cells_infercnv <- cnv_cells[cnv_cells %in% sampled_cells]
cnv_sampled <- cnv_matrix[, sampled_cells_infercnv, drop = FALSE]

normalize_to_minus1_1 <- function(x) {
  2 * (x - min(x)) / (max(x) - min(x)) - 1
}
cnv_norm <- t(apply(cnv_sampled, 1, normalize_to_minus1_1))
cnv_score <- data.frame(
  cell = colnames(cnv_norm),
  cnv_score = colSums(cnv_norm^2)
)
epi_sampled$cnv_score <- cnv_score$cnv_score[match(colnames(epi_sampled), cnv_score$cell)]


cell_cnv_score <- colSums(abs(infercnv_obj@expr.data - 1))




#正确版本 
epi_T <- readRDS("epi_T.rds")
T <- subset(epi_T, seurat_clusters %in% c(0,1))
epi <- readRDS("epi_gsva.rds")

# 获取表达矩阵（不合并）
T_expr_matrix <- GetAssayData(T, layer = "counts")
epi_expr_matrix <- GetAssayData(epi, layer = "counts")

# 找到共同基因
common_genes <- intersect(rownames(T_expr_matrix), rownames(epi_expr_matrix))
T_expr_matrix <- T_expr_matrix[common_genes, ]
epi_expr_matrix <- epi_expr_matrix[common_genes, ]

# --------------- 关键修改开始 ---------------
# 不合并矩阵，分别处理
T_matrix <- T_expr_matrix  # 所有 T 细胞
Epi_matrix <- epi_expr_matrix  # 所有 Epi 细胞

# 仅对 Epi 细胞抽样
set.seed(123)
n_samples <- 10000
sampled_Epi_cells <- sample(colnames(Epi_matrix), size = n_samples, replace = FALSE)
sampled_T_cells <- sample(colnames(T_matrix), size = n_samples, replace = FALSE)
# 构建 infercnv 输入（参考组全量 + 目标组抽样）
infercnv_input <- cbind(
  T_matrix[, sampled_T_cells],  # 所有 T 细胞（参考组）
  Epi_matrix[, sampled_Epi_cells]  # 抽样的 Epi 细胞
)

# 生成注释文件
annotations_file <- data.frame(
  cells = colnames(infercnv_input),
  cell_type = ifelse(
    colnames(infercnv_input) %in% colnames(T_matrix),
    "T cells",
    "Epi"
  )
)
groupFiles <- 'groupFiles.txt'
write.table(annotations_file, file = groupFiles, sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE)

# 基因注释（保持原代码）
gene_names <- rownames(infercnv_input)
gene_annotation <- annoGene(IDs = gene_names, ID_type = "SYMBOL", species = "human")
gene_location_ref <- data.frame(
  gene = gene_annotation$SYMBOL,
  chr = gene_annotation$chr,
  start = gene_annotation$start,
  stop = gene_annotation$end
)
write.table(gene_location_ref[!duplicated(gene_location_ref$gene), ], "geneFile.txt", sep = "\t", row.names = FALSE, col.names = FALSE)

# 运行 infercnv（禁用过滤）
infercnv_obj <- CreateInfercnvObject(
    raw_counts_matrix = infercnv_input,  # 表达矩阵文件路径
    annotations_file = groupFiles,  # 注释文件路径
    delim = "\t",  # 文件分隔符
    gene_order_file = geneFile,  # 基因位置参考文件路径
    ref_group_names = c("T cells")  # 以 T 细胞作为参考
)

infercnv_obj <- infercnv::run(
    infercnv_obj,
    cutoff = 0.2,  # 设置 cutoff 值
    out_dir = "infercnv_output",  # 输出目录
    cluster_by_groups = TRUE,  # 按组聚类
    denoise = TRUE,  # 去噪
    HMM = FALSE  # 是否使用 HMM 模型
)

# 提取 Epi 细胞的 CNV 矩阵（顺序一致）
Epi_cnv_matrix <- infercnv_obj@expr.data[, colnames(Epi_matrix) %in% sampled_Epi_cells]

# 赋值到 Seurat 对象
epi_sampled <- epi[, sampled_Epi_cells]
epi_sampled$cnv_score <- colSums(t(apply(Epi_cnv_matrix, 1, function(x) {
  norm_x <- 2 * ((x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))) - 1
  norm_x^2
})))

epi_sampled <- epi[, sampled_Epi_cells]
epi_in_expr <- intersect(colnames(epi_sampled), colnames(infercnv_obj@expr.data))
epi_sampled$cnv_score <- colSums(t(apply(infercnv_obj@expr.data[, epi_in_expr], 1, function(x) {
    ((x - mean(x)) / sd(x))^2  # Z-score平方和
})))

summary(as.vector(epi_sampled@meta.data$cnv_score))

library(ggplot2)
library(ggpubr)
library(rstatix)
library(ggsci)
data <- epi_sampled@meta.data %>%
    mutate(seurat_clusters = as.factor(seurat_clusters))

# 2. 扩展NPG调色板（17色，Nature风格）
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(17)

# 3. 绘制无点箱图（无网格/无标题/大字体）
ggplot(data, aes(x = seurat_clusters, y = cnv_score, fill = seurat_clusters)) +
    geom_boxplot(
        outlier.shape = NA,       # 隐藏极值点
        width = 0.7,              # 加宽箱体
        color = "#222222",        # 黑色边框
        alpha = 0.8               # 半透明填充
    ) +
    scale_y_continuous(
        limits = c(1000, 10000),  # 固定Y轴范围
        expand = expansion(mult = 0)
    ) +
    scale_fill_manual(values = npg_extended) +  # 扩展Nature配色
    theme_minimal(base_size = 14) +            # 基础字体14pt
    theme(
        panel.grid = element_blank(),            # 移除所有网格线
        plot.title = element_blank(),            # 无标题
        axis.title = element_text(face = "bold"),# 加粗轴标题
        axis.text.x = element_text(angle = 45, hjust = 1), # 倾斜X轴标签
        legend.position = "none",                # 隐藏图例
        axis.line = element_line(color = "black", size = 0.8) # 强化坐标轴
    ) +
    labs(x = "Epi亚群", y = "CNV分数")  # 保留大字体坐标轴

#signf
pairwise_stats <- data %>%
  mutate(seurat_clusters = as.factor(seurat_clusters)) %>%  # 确保是因子
  pairwise_wilcox_test(
    formula = cnv_score ~ seurat_clusters,  # 公式正确
    p.adjust.method = "fdr"
  ) %>%
  # 格式化表格（此时group1/group2已生成）
  mutate(
    Comparison = paste(group1, group2, sep = " vs "),
    significance = case_when(
      p.adj < 0.001 ~ "***",
      p.adj < 0.01 ~ "**",
      p.adj < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  select(Comparison, p, p.adj, significance) %>%
  arrange(p.adj)

library(Seurat)

FeaturePlot(
    object = epi_sampled,
    features = "cnv_score",
    reduction = "umap",  # 替换为你的降维方法
    cols = c("#F5F5F5", "#FF6B6B", "#DC143C"),  # 浅灰→红→深红（>10000固定深红）
    min.cutoff = 1000,  # 1000以下统一浅灰
    max.cutoff = 10000,  # 10000以上固定深红
    pt.size = 0.8,
    blend = TRUE,
    order = TRUE
) +
    theme_void() +
    theme(
        legend.position = "right",
        legend.direction = "vertical",
        legend.key.width = unit(0.8, "cm"),
        plot.margin = margin(10, 20, 10, 10)  # 扩大右侧空间放图例
    )


#T
T_sampled <- T[, sampled_T_cells]
T_in_expr <- intersect(colnames(T_sampled), colnames(infercnv_obj@expr.data))
T_sampled$cnv_score <- colSums(t(apply(infercnv_obj@expr.data[, T_in_expr], 1, function(x) {
    ((x - mean(x)) / sd(x))^2  # Z-score平方和
})))
summary(as.vector(T_sampled@meta.data$cnv_score))
t_mean <- 4501.0  # T细胞第三分位数（标量）

epi_nonparametric <- epi_sampled@meta.data %>%
  as_tibble() %>%
  mutate(seurat_clusters = as.character(seurat_clusters)) %>%  # 确保字符型分组
  group_by(seurat_clusters) %>%
  group_modify(~ {  # 对每个亚群独立操作
    cluster_data <- .x$cnv_score  # 提取原始CNV分数（非均值）
    
    # 单样本 Wilcoxon 检验（中位数 vs t_q3）
    test_result <- wilcox.test(
      cluster_data,
      mu = t_mean,
      alternative = "greater", 
      exact = FALSE  # 大样本时使用正态近似
    )
    
    # 整理结果
    tibble(
      n_cells = nrow(.x),
      median_cnv = median(cluster_data),
      p_value = test_result$p.value,
      statistic = test_result$statistic
    )
  }) %>%
  ungroup() %>%
  mutate(
    p_adj = p.adjust(p_value, "fdr"),  # FDR校正
    significance = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE          ~ "ns"
    )
  )


#文献：
# 0. 依赖包
library(Seurat)
library(dplyr)

# 1. 数据对齐
epi_sampled <- epi_sampled[, colnames(epi_sampled) %in% colnames(infercnv_obj@expr.data)]

# 2. 定义高置信度恶性基础集
cutoff <- quantile(epi_sampled$cnv_score, 0.95)
high_cnv_cells <- WhichCells(epi_sampled, cnv_score >= cutoff)

# 3. 计算参考谱（Top5%平均）
avg_high_cnv <- rowMeans(infercnv_obj@expr.data[, high_cnv_cells])

# 4. 单细胞相关性（核心步骤）
epi_sampled$cnv_cor <- apply(
  infercnv_obj@expr.data[, colnames(epi_sampled)],
  2,
  cor,
  y = avg_high_cnv,
  method = "spearman"  # 或"spearman"
)

# 5. 分类（双重筛选）
epi_sampled$malignant <- ifelse(
  epi_sampled$cnv_cor > 0.15 & epi_sampled$cnv_score >= cutoff*0.9,  # 放宽CNV阈值
  "Malignant",
  "Normal"
)

# 6. 保存结果
saveRDS(epi_sampled, "epi_malignant_annotated.rds")

# 原代码第1行改为聚焦恶性标签（删除sample_type相关）
unique_malignant <- c("Malignant", "Normal")  # 固定顺序

# 原循环删除，直接绘制单图（颜色逻辑完全保留！）
p <- DimPlot(epi_sampled, reduction = "umap", label = FALSE, pt.size = 1, group.by = "malignant") +  # 改1：group.by
  scale_color_manual(
    values = setNames(  # 改2：仅替换分组名，颜色函数不变
      ifelse(unique_malignant == "Malignant", pal_npg()(1), "gray"),  # 复用原npg红逻辑
      unique_malignant
    )
  ) +
  ggtitle("恶性细胞分布") +  # 标题更名
  theme(
    legend.position = "right",  # 恢复图例（原代码隐藏）
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.text = element_text(size = 16, color = "black")
  )

# 保存代码完全保留
png("malignant_umap.png", width = 6000, height = 3000, res = 300)
print(p)
dev.off()