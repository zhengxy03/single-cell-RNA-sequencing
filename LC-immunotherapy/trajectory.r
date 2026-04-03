library(Seurat)
library(monocle)
library(ggplot2)
library(ggsci)
library(igraph)  # 添加这一行

# 加载完整修复脚本
source("monocle_fix_complete.R")

Malignant <- readRDS("malignant_unknown_anno_0.6.rds")

# 随机抽样20000个细胞
random_seed <- sample(1:10000, 1)
set.seed(random_seed)
cat("随机种子:", random_seed, "\n")
sample_cells <- sample(colnames(Malignant), size = 30000)
Malignant_sub <- Malignant[, sample_cells]
cat("抽样后细胞数:", ncol(Malignant_sub), "\n")

# 使用 RNA assay
Malignant_sub <- NormalizeData(Malignant_sub, assay = "RNA")
Malignant_sub <- FindVariableFeatures(Malignant_sub, nfeatures = 2000, assay = "RNA")
hvgs <- VariableFeatures(Malignant_sub, assay = "RNA")

# 提取表达矩阵 - 保持稀疏格式（不转换）
expression_matrix <- LayerData(Malignant_sub, assay = "RNA", layer = "counts")
# 不要转成密集矩阵！保持 sparse matrix

cell_metadata <- Malignant_sub@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(expression_matrix))
rownames(gene_annotation) <- rownames(expression_matrix)

# 只保留存在的基因
hvgs_valid <- hvgs[hvgs %in% rownames(expression_matrix)]
cat("有效高变基因数:", length(hvgs_valid), "\n")

# 提取高变基因表达矩阵（保持稀疏）
expression_matrix_hvgs <- expression_matrix[hvgs_valid, ]
gene_annotation_hvgs <- gene_annotation[hvgs_valid, , drop = FALSE]

cat("高变基因数:", nrow(expression_matrix_hvgs), "\n")
cat("细胞数:", ncol(expression_matrix_hvgs), "\n")
cat("矩阵类型:", class(expression_matrix_hvgs)[1], "\n")

# 创建 CellDataSet 对象
cds <- newCellDataSet(expression_matrix_hvgs,
                      phenoData = new("AnnotatedDataFrame", data = cell_metadata),
                      featureData = new("AnnotatedDataFrame", data = gene_annotation_hvgs),
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())

# 归一化
cds <- estimateSizeFactors(cds)

# 估计离散度（如果报错，用 pooled 方法）
cds <- estimateDispersions(cds, method = "pooled")
cds <- setOrderingFilter(cds, hvgs_valid)

# 降维
cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')
gc()
cds <- orderCells(cds)
# 保存
saveRDS(cds, "malignant_unknown_cds_30000cells_hvg_6.rds")