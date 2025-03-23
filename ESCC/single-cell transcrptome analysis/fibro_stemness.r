library(Seurat)
library(GSVA)
library(clusterProfiler)
library(org.Hs.eg.db)
fibroblasts <- readRDS("fibro_modify.rds")
stem_pathways <- read.gmt("stem_genesets.v2024.1.Hs.gmt")
stemness_genes <- unique(stem_pathways$gene)

fibroblasts_genes <- rownames(fibroblasts)
matched_genes <- stemness_genes[stemness_genes %in% fibroblasts_genes]
length(matched_genes)
stemness_genes <- matched_genes

library(parallel)

num_cores <- detectCores()
cluster <- makeCluster(num_cores - 1)

clusterExport(cluster, c("fibroblasts", "stemness_genes", "AddModuleScore"))

n_splits <- num_cores - 1 # 分割的子集数量，与使用的核心数相同
cell_splits <- split(Cells(fibroblasts), cut(seq_along(Cells(fibroblasts)), n_splits))

compute_module_score <- function(cells_subset) {
  subset_fibroblasts <- subset(fibroblasts, cells = cells_subset)
  subset_fibroblasts <- AddModuleScore(subset_fibroblasts, features = stemness_genes, name = "Stemness_Score")
  return(subset_fibroblasts)
}

subset_results <- parLapply(cluster, cell_splits, compute_module_score)


merged_fibroblasts <- merge(subset_results[[1]], y = subset_results[-1])

saveRDS(merged_fibroblasts, file = "fibro_stemness.rds")
# 停止集群
stopCluster(cluster)










fibroblasts <- AddModuleScore(fibroblasts, features = stemness_genes, name = "Stemness_Score")
head(fibroblasts@meta.data)
