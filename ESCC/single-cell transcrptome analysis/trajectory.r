library(Seurat)
library(monocle)
library(ggplot2)
trace('project2MST', edit = T, where = asNamespace("monocle"))

setwd("/share/home/wangq/zxy/ESCC")
t_cells <- readRDS("t_cells.rds")
CD8 <- subset(t_cells, subset = seurat_clusters %in% c(0, 5, 8))
expression_matrix <- LayerData(CD8, assay = "RNA", layer = "data")
# 提取细胞元数据
cell_metadata <- CD8@meta.data

# 提取基因元数据
gene_annotation <- data.frame(gene_short_name = rownames(expression_matrix))
rownames(gene_annotation) <- rownames(expression_matrix)
cds <- newCellDataSet(expression_matrix,
                      phenoData = new("AnnotatedDataFrame", data = cell_metadata),
                      featureData = new("AnnotatedDataFrame", data = gene_annotation),
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

# 选择高变基因
cds <- detectGenes(cds, min_expr = 0.1)
print(head(fData(cds)))

disp_table <- dispersionTable(cds)
ordering_genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
length(ordering_genes)
cds <- setOrderingFilter(cds, ordering_genes)

# 降维
cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')

# 轨迹推断
cds <- orderCells(cds)
cds <- orderCells(cds, root_state  = 4)




#抽样
expression_matrix <- LayerData(fibroblasts, assay = "RNA", layer = "data")
cell_metadata <- fibroblasts@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(expression_matrix))
rownames(gene_annotation) <- rownames(expression_matrix)

# 随机选择部分细胞
set.seed(123)
subset_cells <- sample(colnames(expression_matrix), size = 20000)  # 选择 10,000 个细胞
expression_matrix <- expression_matrix[, subset_cells]
cell_metadata <- cell_metadata[subset_cells, ]

# 创建 CellDataSet 对象
cds <- newCellDataSet(expression_matrix,
                      phenoData = new("AnnotatedDataFrame", data = cell_metadata),
                      featureData = new("AnnotatedDataFrame", data = gene_annotation),
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())

# 归一化并估计离散度
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

# 检测基因并选择高变基因
cds <- detectGenes(cds, min_expr = 0.1)
disp_table <- dispersionTable(cds)
ordering_genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
cds <- setOrderingFilter(cds, ordering_genes)

cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')

# 轨迹推断
cds <- orderCells(cds)



# 绘制轨迹图
p <- plot_cell_trajectory(cds, color_by = "seurat_clusters")
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(25)
p <- p + scale_color_manual(values = npg_extended)

p2 <- plot_cell_trajectory(cds, color_by = "seurat_clusters") + facet_wrap("~seurat_clusters", nrow = 2) + scale_color_manual(values = npg_extended)

p3 <- plot_cell_trajectory(cds, color_by = "cell_type") + scale_color_manual(values = npg_extended)

p4 <- plot_cell_trajectory(cds, color_by = "cell_type") + facet_wrap("~cell_type", nrow = 2) + scale_color_manual(values = npg_extended)

p <- plot_cell_trajectory(cds, color_by = "sample_type") +
    scale_color_manual(values = npg_extended)

# 保存为PNG文件
ggsave("trajectory_plot_epi.png", plot = p, width = 8, height = 6, dpi = 300)
ggsave("trajectory_plot_epi_split.png", plot = p2, width = 16, height = 6, dpi = 300)


start_color <- npg_extended[1]
end_color <- npg_extended[3]

# 绘制轨迹图并设置连续型颜色比例尺
p5 <- plot_cell_trajectory(cds, color_by = "Pseudotime") +
  scale_color_gradient(low = start_color, high = end_color)

# 保存图片
ggsave("trajectory_plot_fibro_pseudo.png", plot = p5, width = 16, height = 6, dpi = 300)


p <- plot_cell_trajectory(cds, color_by = "period1") +
    scale_color_manual(values = npg_extended)
ggsave("trajectory_plot_fibro_period1.png", plot = p, width = 16, height = 6, dpi = 300)

p <- plot_cell_trajectory(cds, color_by = "period1") + facet_wrap("~period1", nrow = 2) +
    scale_color_manual(values = npg_extended)
ggsave("trajectory_plot_fibro_period1_split.png", plot = p, width = 16, height = 6, dpi = 300)



# 创建一个新的组合分组变量
cds$combined_group <- paste(cds$period1, cds$sample_type, sep = "_")

# 获取组合分组的唯一值
unique_combined_groups <- unique(cds$combined_group)

# 生成颜色映射
n_colors <- length(unique_combined_groups)
npg_extended <- pal_npg("nrc")(n_colors)
color_mapping <- setNames(npg_extended, unique_combined_groups)

# 绘制轨迹图
p <- plot_cell_trajectory(cds, color_by = "combined_group") +
  scale_color_manual(values = color_mapping) + facet_wrap("~combined_group", nrow = 3) +
  ggtitle("Trajectory Plot by Period1 and Sample Type")

# 保存图形
ggsave("trajectory_plot_fibro_period1_sample_type_split.png", plot = p, width = 10, height = 6, dpi = 300)


period_info <- pData(cds)$period1  # 假设分期信息在 "period2" 列，可按需修改
unique_periods <- unique(period_info)

# 创建一个空列表来存储每个分期的绘图
plot_list <- list()

for (period in unique_periods) {
  # 筛选出当前分期的细胞
  period_cds <- cds[, pData(cds)$period1 == period]
  
  # 打印筛选后的数据维度
  print(paste("Period:", period))
  print(dim(period_cds))
  
  # 重新进行降维
  cat("Before reduceDimension:\n")
  print(dim(period_cds))
  period_cds <- reduceDimension(period_cds, max_components = 2, method = 'DDRTree')
  cat("After reduceDimension:\n")
  print(dim(period_cds))
  
  # 重新进行轨迹推断
  cat("Before orderCells:\n")
  print(dim(period_cds))
  period_cds <- orderCells(period_cds)
  cat("After orderCells:\n")
  print(dim(period_cds))
  
  # 绘制轨迹图，按 cell_type 着色
  p <- plot_cell_trajectory(period_cds, color_by = "cell_type") +
    ggtitle(paste("Trajectory for Period", period)) +
    scale_color_manual(values = npg_extended)
  
  # 将图添加到列表中
  plot_list[[as.character(period)]] <- p
}
plot_list <- list(plot_list[[2]], plot_list[[1]], plot_list[[3]])
# 使用 patchwork 包将所有图组合在一起
library(patchwork)
combined_plot <- wrap_plots(plot_list, ncol = 3)  # 每行 2 个图，可按需调整

# 保存组合后的图
png("trajectory_by_fibro_period1_celltype.png", width = 6000, height = 3000, res = 300)
print(combined_plot)
dev.off()

combined_plot <- wrap_plots(lapply(plot_list, function(p) {
  p + facet_wrap(~sample_type, nrow = 2)
}), ncol = 2)

# 保存组合后的图
png("trajectory_by_fibro_period1_celltype_sampletype.png", width = 6000, height = 3000, res = 300)
print(combined_plot)
dev.off()

#密度图
df <- pData(cds)
density_plot <- ggplot(df, aes(Pseudotime, colour = cell_type, fill = cell_type)) +
  geom_density(bw = 0.5, size = 1, alpha = 0.5) +  # 密度图
  scale_color_manual(values = npg_extended) +      # 设置线条颜色
  scale_fill_manual(values = npg_extended) +       # 设置填充颜色（与线条颜色一致）
  theme(
    panel.background = element_blank(),           # 去除背景颜色
    panel.grid = element_blank(),                 # 去除网格线
    axis.line = element_line(color = "black")     # 保留坐标轴线
  )
ggsave("cell_density_along_pseudotime.png", plot = density_plot, width = 16, height = 6, dpi = 300)

density_plot_split <- ggplot(df, aes(Pseudotime, colour = cell_type, fill = cell_type)) +
  geom_density(bw = 0.5, size = 1, alpha = 0.5) +  # 密度图
  scale_color_manual(values = npg_extended) +      # 设置线条颜色
  scale_fill_manual(values = npg_extended) +       # 设置填充颜色（与线条颜色一致）
  theme(
    panel.background = element_blank(),           # 去除背景颜色
    panel.grid = element_blank(),                 # 去除网格线
    axis.line = element_line(color = "black")     # 保留坐标轴线
  ) + 
    facet_wrap("~cell_type", nrow = 2)
ggsave("cell_density_along_pseudotim_split.png", plot = density_plot_split, width = 16, height = 6, dpi = 300)

#差异基因
expressed_genes=row.names(subset(fData(cds),num_cells_expressed>=10)) #在部分基因里面找
pseudotime_de <- differentialGeneTest(cds[expressed_genes,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
pseudotime_de <- pseudotime_de[order(pseudotime_de$qval), ]
states_de <- differentialGeneTest(cds[expressed_genes,],
                                  fullModelFormulaStr = "~State")
states_de <- states_de[order(states_de$qval), ]

saveRDS(cds, file = "test_monocle.rds")
write.table(pseudotime_de, file = "pseudotime_de.rds", quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
write.table(states_de, file = "states_de.rds", quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)


#驱动基因
trace('buildBranchCellDataSet', edit = T, where = asNamespace("monocle"))
beam_res <- BEAM(
  cds, 
  branch_point = 1,          # 分支点编号（通过plot_cell_trajectory确定）
  cores = 4,                  # 多线程加速
  progenitor_method = "duplicate"
)