# 改进的 CellCharter 风格聚类（R 实现）
# 适用于 mPT 对象，坐标在 metadata，细胞类型在 detailed
# 不强制空间连续性，基于局部邻居组成聚类

library(Seurat)
library(Matrix)
library(mclust)
library(RANN)
library(ggplot2)
library(pheatmap)
mPT <- readRDS("mPT_detailed.rds")


# ==================== 参数设置 ====================
k_neighbors <- 20          # 最近邻数
n_layers <- 3              # 聚合邻居层数（0=自身，1=1阶邻居，2=2阶邻居…）
use_pca_dims <- 1:20       # 使用的 PCA 维度
candidate_k <- 5:20        # 待评估的聚类数范围
n_repeats <- 10            # 每个 k 的重复次数（抽样数据上的重复次数）
model_name <- "VVV"        # 固定协方差模型，提高稳定性
x_coord <- "CenterX_global_px"
y_coord <- "CenterY_global_px"
sample_col <- "sample"

# 抽样大小（用于稳定性分析）
sample_size <- 10000       # 可根据数据量调整，如 5000-20000

# -------------------- 1. 提取坐标和 PCA 嵌入（全量）--------------------
coords <- mPT@meta.data[, c(x_coord, y_coord)]
coords <- as.matrix(coords)

# 确保已运行 PCA
if (!"pca" %in% names(mPT@reductions)) {
  mPT <- RunPCA(mPT, npcs = 30)
}
pca_emb <- Embeddings(mPT, reduction = "pca")[, use_pca_dims]

n_cells <- nrow(coords)

# -------------------- 2. 构建 k 近邻图（全量）--------------------
nn <- nn2(coords, coords, k = k_neighbors + 1)
adj_list <- lapply(1:n_cells, function(i) nn$nn.idx[i, 2:(k_neighbors+1)])  # 排除自身

# -------------------- 3. 计算多圈邻居的平均 PCA 特征（全量）--------------------
n_pca <- length(use_pca_dims)
aggregated_features <- matrix(0, n_cells, (n_layers + 1) * n_pca)

# 第 0 层：自身
aggregated_features[, 1:n_pca] <- pca_emb

# 第 1 至 n_layers 层
current_adj <- adj_list
for (layer in 1:n_layers) {
  layer_means <- t(sapply(1:n_cells, function(i) {
    neigh <- current_adj[[i]]
    if (length(neigh) == 0) return(rep(0, n_pca))
    colMeans(pca_emb[neigh, , drop = FALSE])
  }))
  
  start_col <- layer * n_pca + 1
  end_col <- (layer + 1) * n_pca
  aggregated_features[, start_col:end_col] <- layer_means
  
  if (layer < n_layers) {
    new_adj <- vector("list", n_cells)
    for (i in 1:n_cells) {
      expanded <- unique(c(current_adj[[i]], unlist(current_adj[current_adj[[i]]])))
      expanded <- expanded[expanded != i]
      new_adj[[i]] <- expanded
    }
    current_adj <- new_adj
  }
}

# 标准化特征（有助于 GMM）
aggregated_features <- scale(aggregated_features)

# -------------------- 4. 抽样进行稳定性分析 --------------------
set.seed(123)
sample_idx <- sample(1:n_cells, size = min(sample_size, n_cells))
features_sample <- aggregated_features[sample_idx, ]

cat("抽样细胞数:", length(sample_idx), "\n")
cat("开始稳定性分析，候选 k 范围:", paste(candidate_k, collapse = ", "), "\n")

all_labels <- list()
for (k in candidate_k) {
  cat("\n当前 k =", k, "\n")
  for (rep in 1:n_repeats) {
    set.seed(rep * 12345)
    fit <- Mclust(features_sample, G = k, modelNames = model_name, verbose = FALSE)
    all_labels[[paste(k, rep, sep = "_")]] <- fit$classification
  }
}

# 计算每个 k 的平均稳定性（ARI）
stability_scores <- sapply(candidate_k, function(k) {
  idx <- grep(paste0("^", k, "_"), names(all_labels))
  labels_list <- all_labels[idx]
  n_reps <- length(labels_list)
  if (n_reps < 2) return(NA)
  
  ari_vals <- c()
  for (i in 1:(n_reps-1)) {
    for (j in (i+1):n_reps) {
      ari_vals <- c(ari_vals, adjustedRandIndex(labels_list[[i]], labels_list[[j]]))
    }
  }
  mean(ari_vals, na.rm = TRUE)
})

# 绘制稳定性曲线
pdf("cn_stability_curve.pdf", width = 8, height = 6)
plot(candidate_k, stability_scores, type = "b", 
     xlab = "Number of clusters (k)", ylab = "Stability (mean ARI)",
     main = "Cluster Stability (sampling)")
abline(v = candidate_k[which.max(stability_scores)], col = "red", lty = 2)
text(candidate_k[which.max(stability_scores)], 
     max(stability_scores, na.rm = TRUE), 
     labels = paste("k =", candidate_k[which.max(stability_scores)]),
     pos = 3, col = "red")
dev.off()

# 选择最佳 k
best_k <- candidate_k[which.max(stability_scores)]
cat("\n最佳聚类数 (最高稳定性):", best_k, "\n")

# -------------------- 5. 用最佳 k 对全量数据进行最终聚类 --------------------
set.seed(123)
final_fit <- Mclust(aggregated_features, G = best_k, modelNames = model_name, verbose = FALSE)
mPT$cn_label <- factor(final_fit$classification)

# -------------------- 6. 可视化与保存 --------------------
# 空间分布图（按样本分面）
plot_df <- cbind(coords, cn = mPT$cn_label, sample = mPT@meta.data[[sample_col]])

# 使用 aes() 替代 aes_string
p <- ggplot(plot_df, aes(x = !!sym(x_coord), y = !!sym(y_coord), color = cn)) +
  geom_point(size = 0.3, alpha = 0.6) +
  facet_wrap(~ sample, scales = "free") +
  # coord_fixed()  # 删除此行，因为与自由刻度冲突
  theme_void() +
  labs(title = paste("Cell Neighborhoods (k =", best_k, ")")) +
  theme(legend.position = "bottom")

ggsave("cn_spatial_distribution.pdf", plot = p, width = 12, height = 8)
# 细胞组成热图
cn_comp <- table(mPT$cn_label, mPT$detailed)
cn_prop <- prop.table(cn_comp, margin = 1)
p2 <- pheatmap(cn_prop, 
         main = paste("Cell composition by CN (k =", best_k, ")"),
         fontsize = 8,
         filename = "cn_composition_heatmap.pdf",
         width = 10, height = 8)
ggsave("cn_cellgroup_heatmap.pdf", plot = p2, width = 12, height = 8)

# 保存带有 CN 标签的 Seurat 对象
saveRDS(mPT, "mPT_with_cn.rds")

cat("\n分析完成！结果已保存。\n")