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

# 绘制轨迹图
p <- plot_cell_trajectory(cds, color_by = "seurat_clusters")
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(15)
p <- p + scale_color_manual(values = npg_extended)

p2 <- plot_cell_trajectory(cds, color_by = "seurat_clusters") + facet_wrap("~seurat_clusters", nrow = 2) + scale_color_manual(values = npg_extended)

p3 <- plot_cell_trajectory(cds, color_by = "cell_type") + scale_color_manual(values = npg_extended)

p4 <- plot_cell_trajectory(cds, color_by = "cell_type") + facet_wrap("~cell_type", nrow = 2) + scale_color_manual(values = npg_extended)



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