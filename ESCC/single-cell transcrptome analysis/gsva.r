library(Seurat)
library(GSVA)
library(clusterProfiler)
library(org.Hs.eg.db)

hallmark_pathways <- read.gmt("h.all.v2024.1.Hs.symbols.gmt")
hallmark_list <- split(hallmark_pathways$gene, hallmark_pathways$term)
epi <- readRDS("epi_filtered.rds")
expression_matrix <- GetAssayData(epi, layer = "data")
param <- GSVA::gsvaParam(
  exprData = expression_matrix,  # 表达矩阵
  geneSets = hallmark_list,      # 基因集列表
  kcdf = "Gaussian"              # 核函数
)
gsva_results <- GSVA::gsva(param)
saveRDS(gsva_results, file = "epi_gsva_results.rds")
gsva_results <- t(gsva_results)

colnames(gsva_results) <- gsub("HALLMARK_", "", colnames(gsva_results))
epi@meta.data <- cbind(epi@meta.data, gsva_results)

epi@meta.data$Epi_cluster <- paste0("C", epi@meta.data$seurat_clusters)

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
saveRDS(top_50_gsva,file="epi_top50.rds")
mean_gsva_by_epi_cluster <- aggregate(top_50_gsva, by = list(epi@meta.data$Epi_cluster), FUN = mean)
rownames(mean_gsva_by_epi_cluster) <- mean_gsva_by_epi_cluster$Group.1
mean_gsva_by_epi_cluster <- mean_gsva_by_epi_cluster[, -1]  # 去掉分组列

# 转置矩阵，使行为通路，列为 Epi_cluster
mean_gsva_by_epi_cluster <- t(mean_gsva_by_epi_cluster)

library(pheatmap)
pdf("gsva_heatmap.pdf", width = 10, height = 8)
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
dev.off()