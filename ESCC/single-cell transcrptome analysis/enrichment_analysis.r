library(Seurat)
library(GSVA)
library(clusterProfiler)
library(org.Hs.eg.db)

hallmark_pathways <- read.gmt("h.all.v2024.1.Hs.symbols.gmt")
hallmark_list <- split(hallmark_pathways$gene, hallmark_pathways$term)

expression_matrix <- GetAssayData(epi, layer = "data")
param <- GSVA::gsvaParam(
  exprData = expression_matrix,  # 表达矩阵
  geneSets = hallmark_list,      # 基因集列表
  kcdf = "Gaussian"              # 核函数
)
gsva_results <- GSVA::gsva(param)
gsva_results <- t(gsva_results)
colnames(gsva_results) <- gsub("HALLMARK_", "", colnames(gsva_results))
epi@meta.data <- cbind(epi@meta.data, gsva_results)

epi@meta.data$Epi_cluster <- paste0("Epi", epi@meta.data$seurat_clusters)

#top 50
subtypes <- unique(epi@meta.data$seurat_clusters)

p_values <- numeric(length = ncol(gsva_results))
names(p_values) <- colnames(gsva_results)

for (pathway in colnames(gsva_results)) {
  # 提取当前通路的 GSVA 得分
  scores <- gsva_results[, pathway]
  
  # 分组比较（例如亚型1 vs 亚型2）
  group1 <- scores[epi@meta.data$seurat_clusters == subtypes[1]]
  group2 <- scores[epi@meta.data$seurat_clusters == subtypes[2]]
  
  # 进行 t 检验
  test_result <- t.test(group1, group2)
  
  # 保存 p 值
  p_values[pathway] <- test_result$p.value
}
top_50_pathways <- names(sort(p_values))[1:50]
print(top_50_pathways)

top_50_gsva <- gsva_results[, top_50_pathways]
mean_gsva_by_epi_cluster <- aggregate(top_50_gsva, by = list(epi@meta.data$Epi_cluster), FUN = mean)
rownames(mean_gsva_by_epi_cluster) <- mean_gsva_by_epi_cluster$Group.1
mean_gsva_by_epi_cluster <- mean_gsva_by_epi_cluster[, -1]  # 去掉分组列

# 转置矩阵，使行为通路，列为 Epi_cluster
mean_gsva_by_epi_cluster <- t(mean_gsva_by_epi_cluster)

library(pheatmap)
pheatmap(
  mean_gsva_by_epi_cluster,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_colnames = TRUE,
  show_rownames = TRUE,
  fontsize_row = 8,
  fontsize_col = 10
)

#fibroblasts
#多核运行，limma比较差异性
library(Seurat)
library(GSVA)
library(clusterProfiler)
library(org.Hs.eg.db)
library(BiocParallel)
library(parallel)
library(doParallel)
library(limma)
library(pheatmap)
fibroblasts <- readRDS("fibro_modify.rds")
hallmark_pathways <- read.gmt("h.all.v2024.1.Hs.symbols.gmt")
hallmark_list <- split(hallmark_pathways$gene, hallmark_pathways$term)

expression_matrix <- GetAssayData(fibroblasts, layer = "data")
param <- GSVA::gsvaParam(
  exprData = expression_matrix,
  geneSets = hallmark_list,
  kcdf = "Gaussian"
)


bp <- MulticoreParam(workers = detectCores() - 1)

gsva_results <- GSVA::gsva(param, BPPARAM = bp)

saveRDS(gsva_results, file = "fibro_gsva_results.rds")
rownames(gsva_results) <- gsub("HALLMARK_", "", rownames(gsva_results))
gsva_matrix <- gsva_results  

cluster_info <- fibroblasts@meta.data$seurat_clusters

design <- model.matrix(~ 0 + factor(cluster_info))
colnames(design) <- levels(factor(cluster_info))  # 设置列名为 cluster 名称


fit <- lmFit(gsva_matrix, design)

colnames(design) <- paste0("Fib_", levels(factor(cluster_info)))


contrast_matrix <- makeContrasts(
  contrasts = paste(colnames(design), collapse = "-"),
  levels = design
)
# 计算对比
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)



diff_pathways <- topTable(fit2, number = 50, adjust.method = "BH", sort.by = "p")

# 提取前50条差异通路的名称
top50_pathways <- rownames(diff_pathways)[1:50]

# 提取前50条差异通路的GSVA分数
top50_gsva_matrix <- gsva_matrix[top50_pathways, ]

# 转置矩阵，使行是样本，列是通路
top50_gsva_matrix_transposed <- t(top50_gsva_matrix)
library(tibble)

gsva_df <- as.data.frame(top50_gsva_matrix_transposed)
gsva_df$seurat_clusters <- epi@meta.data[rownames(gsva_df), "seurat_clusters"]
average_pathway_scores <- gsva_df %>%
    group_by(seurat_clusters) %>%
    summarise(across(where(is.numeric), mean, na.rm = TRUE))
average_pathway_scores_matrix <- as.matrix(average_pathway_scores[, -1])  # 移除第一列（seurat_clusters）
rownames(average_pathway_scores_matrix) <- average_pathway_scores$seurat_clusters

# 转置矩阵（行为通路，列为簇）
average_pathway_scores_matrix <- t(average_pathway_scores_matrix)

png("epi_top50_average_pathway_scores_heatmap.png", width = 2000, height = 2000, res = 300)
pheatmap(
    average_pathway_scores_matrix,
    scale = "row",  # 按行标准化
    clustering_method = "complete",
    show_rownames = TRUE,
    show_colnames = TRUE,
    color = colorRampPalette(c("blue", "white", "red"))(100),
    main = "Average GSVA Scores by Cluster"
)
dev.off()




# 保存结果
write.csv(diff_pathways, file = "diff_pathways.csv")

