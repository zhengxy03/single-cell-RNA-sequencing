# 清空环境
rm(list = ls()) 

# 设置内存上限
options(future.globals.maxSize = 50 * 1024^3)  # 50GB

# 加载必要的包
library(Seurat)
library(Matrix)
library(mclust)
library(ggplot2)
library(pheatmap)
library(parallel)  # 用于多线程计算

# 加载修改后的脚本（获取函数定义）
source("cellcharter_R_parallel.R")

# 加载第三个中间文件
intermediate <- readRDS("./cellcharter_results_parallel/intermediates/seurat_before_colnames.rds")
mPT <- intermediate$seurat_obj
features <- intermediate$features

# 从这里开始继续执行后续步骤
cat("Loaded intermediate file successfully!\n")
cat("Features dimensions:", dim(features), "\n")
cat("Seurat object dimensions:", dim(mPT), "\n")

# 转置特征矩阵（如果需要）
if (nrow(features) == ncol(mPT)) {
  features <- t(features)
  cat("Transposed features dimensions:", dim(features), "\n")
}

# 设置列名
colnames(features) <- colnames(mPT)
cat("Set cell names successfully!\n")

# 继续执行聚合过程
n_layers <- 3
aggregations <- "mean"
use_rep <- "pca"
out_key <- "cellcharter"
sample_key <- NULL
x_coord_col <- "CenterX_global_px"
y_coord_col <- "CenterY_global_px"
n_cores <- 20

# 处理 n_layers 参数
if (is.numeric(n_layers) && length(n_layers) == 1) {
  n_layers <- 0:n_layers
}

# 初始化聚合结果
aggregated_features <- list()

# 对每个层进行聚合
for (layer in n_layers) {
  message(paste("Processing layer", layer, "..."))
  
  # 创建邻接矩阵的幂
  if (layer == 0) {
    # 0层表示自身
    adj_matrix_layer <- Matrix::Diagonal(ncol(mPT))
  } else {
    # 对于多层，使用邻接矩阵的幂
    adj_matrix_layer <- Matrix::Matrix(0, nrow = ncol(mPT), ncol = ncol(mPT), sparse = TRUE)
    
    # 计算邻接矩阵的幂
    temp_adj <- mPT@graphs$spatial
    for (i in 1:layer) {
      adj_matrix_layer <- adj_matrix_layer + temp_adj
      temp_adj <- temp_adj %*% mPT@graphs$spatial
    }
    
    # 二值化邻接矩阵
    adj_matrix_layer@x[adj_matrix_layer@x > 0] <- 1
  }
  
  # 移除对角线（避免自循环）
  if (is(adj_matrix_layer, "dgCMatrix")) {
    # 对于稀疏矩阵，设置对角线为0
    for (i in 1:nrow(adj_matrix_layer)) {
      idx <- which(adj_matrix_layer@i + 1 == i & adj_matrix_layer@p[i] < adj_matrix_layer@p[i+1])
      if (length(idx) > 0) {
        pos <- adj_matrix_layer@p[i] + idx
        adj_matrix_layer@x[pos] <- 0
      }
    }
  } else {
    # 对于密集矩阵，设置对角线为0
    diag(adj_matrix_layer) <- 0
  }
  
  # 归一化邻接矩阵
  if (is(adj_matrix_layer, "dgCMatrix")) {
    # 对于稀疏矩阵，计算每行的和
    row_sums <- Matrix::rowSums(adj_matrix_layer)
    row_sums[row_sums == 0] <- 1  # 避免除以0
    
    # 归一化
    for (i in 1:nrow(adj_matrix_layer)) {
      start <- adj_matrix_layer@p[i]
      end <- adj_matrix_layer@p[i+1] - 1
      if (start <= end) {
        adj_matrix_layer@x[start:end] <- adj_matrix_layer@x[start:end] / row_sums[i]
      }
    }
  } else {
    # 对于密集矩阵，归一化
    row_sums <- rowSums(adj_matrix_layer)
    row_sums[row_sums == 0] <- 1  # 避免除以0
    adj_matrix_layer <- adj_matrix_layer / row_sums
  }
  
  # 执行聚合
  if (is(features, "dgCMatrix")) {
    aggregated <- features %*% adj_matrix_layer
  } else {
    aggregated <- features %*% adj_matrix_layer
  }
  
  # 存储聚合结果
  aggregated_features[[as.character(layer)]] <- aggregated
}

# 合并聚合结果
combined_features <- do.call(rbind, aggregated_features)
cat("Combined features dimensions:", dim(combined_features), "\n")

# 创建新的 assay
if (!out_key %in% names(mPT@assays)) {
  # 添加新的 assay
  mPT <- AddAssay(mPT, object = as.matrix(combined_features), assay = out_key)
  cat("Added new assay:", out_key, "\n")
} else {
  # 更新现有 assay
  mPT@assays[[out_key]] <- CreateAssayObject(counts = as.matrix(combined_features))
  cat("Updated existing assay:", out_key, "\n")
}

# 保存聚合后的结果
#saveRDS(mPT, "./cellcharter_results_parallel/intermediates/seurat_after_aggregation.rds")
#cat("Saved seurat_after_aggregation.rds\n")

# 继续执行聚类
message("Clustering cells...")
mPT <- cluster_cells(mPT, n_clusters = 10)

# 保存聚类后的结果
saveRDS(mPT, "./cellcharter_results_parallel/intermediates/seurat_after_clustering.rds")
cat("Saved seurat_after_clustering.rds\n")

# 继续执行空间域分析
message("Analyzing spatial domains...")
mPT <- analyze_spatial_domains(
  mPT, 
  out_dir = "./cellcharter_results_parallel",
  cell_type_col = "detailed",
  n_cores = n_cores
)

# 保存最终结果
#saveRDS(mPT, "./cellcharter_results_parallel/seurat_final.rds")
#cat("Saved final results to seurat_final.rds\n")

# 保存到主对象
saveRDS(mPT, "mPT_cellcharter_results.rds")
cat("Saved final results to mPT_cellcharter_results.rds\n")