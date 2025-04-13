library(Seurat)
library(monocle)
library(ggplot2)
library(ggsci)
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
cds <- orderCells(cds, root_state  = 4)


# 绘制轨迹图
p <- plot_cell_trajectory(cds, color_by = "State") +
  facet_wrap("~State", nrow = 2) + scale_color_manual(values = npg_extended)
p <- plot_cell_trajectory(cds, color_by = "seurat_clusters")
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(15)
p <- p + scale_color_manual(values = npg_extended)

p2 <- plot_cell_trajectory(cds, color_by = "seurat_clusters") + facet_wrap("~seurat_clusters", nrow = 2) + scale_color_manual(values = npg_extended)

p3 <- plot_cell_trajectory(cds, color_by = "cell_type") + scale_color_manual(values = npg_extended) +
  theme(legend.text = element_text(size = 18),
        legend.title = element_blank())

png("traj_fib_celltype.png", width = 6000, height = 3000, res = 300)
print(p3)
dev.off()

p4 <- plot_cell_trajectory(cds, color_by = "cell_type") + facet_wrap("~cell_type", nrow = 2) + scale_color_manual(values = npg_extended) +
  theme(legend.text = element_text(size = 18),
        legend.title = element_blank())
png("traj_fib_celltype_split.png", width = 6000, height = 3000, res = 300)
print(p4)
dev.off()

p4 <- plot_cell_trajectory(cds, color_by = "period1") + facet_wrap("~period1", nrow = 3) + scale_color_manual(values = npg_extended) +
  theme(legend.text = element_text(size = 18),
        legend.title = element_blank())
png("traj_fib_period1_split.png", width = 6000, height = 3000, res = 300)
print(p4)
dev.off()


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


#period and sampletype
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
saveRDS(plot_list, file = "fibro_plot_list.rds")
plot_list <- list(plot_list[[2]], plot_list[[1]], plot_list[[3]])
# 使用 patchwork 包将所有图组合在一起
for (i in 2:length(plot_list)) {
  plot_list[[i]] <- plot_list[[i]] + guides(color = "none")
}
# 组合图形
combined_plot <- wrap_plots(plot_list, ncol = 3) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 18)) # 修改这里来增大图例字体大小

# 保存组合后的图
png("trajectory_by_fibro_period1_celltype.png", width = 6000, height = 3000, res = 300)
print(combined_plot)
dev.off()

#sampletype
faceted_plot_list <- lapply(plot_list, function(p) {
  p + facet_wrap(~sample_type, nrow = 2)
})

# 隐藏除第一个图之外的所有图的图例
for (i in 2:length(faceted_plot_list)) {
  faceted_plot_list[[i]] <- faceted_plot_list[[i]] + guides(color = "none")
}

# 使用 patchwork 组合图形，并设置图例位置
library(patchwork)
combined_plot <- wrap_plots(faceted_plot_list, ncol = 3) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 18))

# 保存组合后的图
png("trajectory_by_fibro_period1_celltype_sampletype.png", width = 6000, height = 3000, res = 300)
print(combined_plot)
dev.off()  


#sampletype and celltype
faceted_plot_list <- lapply(plot_list, function(p) {
  p + facet_grid(cell_type ~ sample_type)
})

# 隐藏除第一个图之外的所有图的图例
for (i in 2:length(faceted_plot_list)) {
  faceted_plot_list[[i]] <- faceted_plot_list[[i]] + guides(color = "none")
}

# 使用 patchwork 组合图形，并设置图例位置
library(patchwork)
combined_plot <- wrap_plots(faceted_plot_list, ncol = 3) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 18))

# 保存组合后的图
png("trajectory_by_fibro_period1_celltype_sampletype_split.png", width = 6000, height = 4000, res = 300)
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

df <- pData(cds)

# 创建自定义的标签函数，移除 cell_type 的显示
my_labeller <- function(variable, value) {
  if (variable == "cell_type") {
    return("")
  }
  return(value)
}

# 创建密度图并添加分面
density_plot <- ggplot(df, aes(Pseudotime, colour = cell_type, fill = cell_type)) +
  geom_density(bw = 0.5, size = 1, alpha = 0.5) +  # 密度图
  scale_color_manual(values = npg_extended) +      # 设置线条颜色
  scale_fill_manual(values = npg_extended) +       # 设置填充颜色（与线条颜色一致）
  theme(
    panel.background = element_blank(),           # 去除背景颜色
    panel.grid = element_blank(),                 # 去除网格线
    axis.line = element_line(color = "black"),     # 保留坐标轴线
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 18),
    strip.text = element_text(size = 22)  # 增大分面标题字体大小
  ) +
  # 使用 facet_grid 进行三维分面，并使用自定义标签
  facet_grid(period1 ~ cell_type + sample_type, labeller = my_labeller)

# 保存图形
png("cell_density_along_pseudotime.png", width = 9000, height = 3000, res = 300)
print(density_plot)
dev.off()

#差异基因
#expressed_genes=row.names(subset(fData(cds),num_cells_expressed>=10)) #在部分基因里面找
ordering_genes <- row.names(subset(fData(cds), use_for_ordering == TRUE))
high_quality_ordering_genes <- row.names(subset(fData(cds), 
                                              use_for_ordering == TRUE & 
                                              num_cells_expressed >= 10 &
                                              mean_expression > 0.5))
stage_markers <- differentialGeneTest(cds[high_quality_ordering_genes,], 
                                    fullModelFormulaStr = "~period1", 
                                    cores = 24)
#stage_markers <- differentialGeneTest(cds[expressed_genes,], 
                                    fullModelFormulaStr = "~period1", 
                                    cores = 4)
stage_markers <- stage_markers[order(stage_markers$qval), ]
write.csv(marker_fib, "stage_markers_fib.csv")
sig_genes <- subset(stage_markers, qval < 0.01)
ordering_genes <- row.names(subset(fData(cds), use_for_ordering == TRUE))
sig_ordering_genes <- sig_genes[row.names(sig_genes) %in% ordering_genes, ]
sig_ordering_genes <- sig_ordering_genes[order(sig_ordering_genes$qval), ]
write.csv(sig_ordering_genes, "stage_markers_fib_sig.csv")

top_genes <- row.names(head(sig_ordering_genes, 52))
genes_to_remove <- c("ENSG00000286533", "LINC01705")
top_genes <- setdiff(top_genes, genes_to_remove)
p <- plot_pseudotime_heatmap(
  cds[top_genes, ],
  show_rownames = TRUE,
  hmcols = viridis::inferno(100),
  cluster_rows = TRUE,
  use_gene_short_name = TRUE,
  return_heatmap = TRUE  # 确保返回ggplot对象
)
ggsave("topgene_pseudotime_heatmap.png", plot = p)

pseudotime_de <- differentialGeneTest(cds[expressed_genes,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
pseudotime_de <- pseudotime_de[order(pseudotime_de$qval), ]
states_de <- differentialGeneTest(cds[expressed_genes,],
                                  fullModelFormulaStr = "~State")
states_de <- states_de[order(states_de$qval), ]

saveRDS(cds, file = "test_monocle.rds")
write.table(pseudotime_de, file = "pseudotime_de.rds", quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
write.table(states_de, file = "states_de.rds", quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)

#pheatmap
pData(cds)$period1 <- factor(pData(cds)$period1, 
                            levels = c("I", "II", "III"))
top_genes <- row.names(head(sig_ordering_genes, 52))

genes_to_remove <- c("ENSG00000286533", "LINC01705")
top_genes <- setdiff(top_genes, genes_to_remove)

# 验证结果
length(top_genes) 
library(pheatmap)

# 按period1排序样本
sample_order <- order(pData(cds)$period1)
expr_ordered <- exprs(cds[top_genes, sample_order])

# 创建分组注释
annotation_col <- data.frame(
  Period = pData(cds)$period1[sample_order],
  row.names = colnames(expr_ordered)
)

# 自定义颜色
period_colors <- c("I" = "#1f77b4", "II" = "#ff7f0e", "III" = "#2ca02c")

mean_expr <- do.call(cbind, lapply(levels(pData(cds)$period1), function(grp) {
  cells_in_group <- row.names(subset(pData(cds), period1 == grp))
  Matrix::rowMeans(exprs(cds[top_genes, cells_in_group]))
}))
colnames(mean_expr) <- levels(pData(cds)$period1)


# 检查结果
head(mean_expr)
annotation_col <- data.frame(
  Period = colnames(mean_expr),  # 列名即为分组
  row.names = colnames(mean_expr)  # 行名与矩阵列名一致
)

# 3. 定义分组颜色
period_colors <- c("I" = "#1f77b4", "II" = "#ff7f0e", "III" = "#2ca02c")



p <- pheatmap(
  mean_expr,
  scale = "row",                 # 按基因标准化（行）
  cluster_cols = FALSE,          # 禁用列聚类 → 保持I-II-III顺序
  show_colnames = TRUE,          # 显示列名（I/II/III）
  show_rownames = TRUE,          # 显示基因名
  color = colorRampPalette(c("#13007D", "white", "#9B0000"))(100),  # 高对比度
  annotation_col = annotation_col,  # 添加分组条
  annotation_colors = list(Period = period_colors),  # 分组颜色
  gaps_col = 1:2,                # 在I-II和II-III之间画分隔线（因只有3列）
  main = "Average Expression by Period",
  border_color = NA,             # 无单元格边框
  cellwidth = 100,                # 固定列宽
  cellheight = 12                # 固定行高
)
ggsave("topgene_period1_heatmap.png", plot = p, width = 8, height = 10)



#驱动基因
trace('buildBranchCellDataSet', edit = T, where = asNamespace("monocle"))
expressed_genes <- row.names(subset(fData(cds), 
                                  num_cells_expressed >= 10))  # 至少在10个细胞中表达
cds_filtered <- cds[expressed_genes, ]

var_genes <- dispersionTable(cds)
var_genes <- var_genes[order(-var_genes$dispersion_empirical), ]
top_genes <- var_genes$gene_id[1:2000]  # 取前2000个高变基因

beam_res <- BEAM(
  cds[top_genes, ],  # 只分析高变基因
  branch_point = 1,
  cores = 22,
  progenitor_method = "sequential_split"
)
sig_genes <- subset(beam_res, qval < 0.05)
sig_genes <- sig_genes[order(sig_genes$qval), ]

gene_list <- c(MSTN, ZNF750, TAS2R19, POTEF, HOXD11, IGLV1-47)
plot_genes_in_pseudotime(
  cds[top_genes[1:5], ],
  color_by = "Branch",
  trend_formula = "~ sm.ns(Pseudotime, df=3)",
  panel_order = top_genes[1:5],  # 控制基因顺序
  ncol = 2  # 分两列显示
) +
  theme(legend.position = "right")