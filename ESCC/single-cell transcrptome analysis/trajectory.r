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

disp_table <- dispersionTable(cds)
ordering_genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
cds <- setOrderingFilter(cds, ordering_genes)

# 降维
cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')

# 轨迹推断
cds <- orderCells(cds)

# 绘制轨迹图
p <- plot_cell_trajectory(cds, color_by = "State")

# 保存为PNG文件
ggsave("trajectory_plot.png", plot = p, width = 8, height = 6, dpi = 300)