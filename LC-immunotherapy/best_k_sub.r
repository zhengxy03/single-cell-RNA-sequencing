# ============================
# CellCharter 最佳 k 值选择
# 只运行 k 值选择，不进行后续分析
# ============================

rm(list = ls())
options(future.glob.maxSize = 100 * 1024^3)

# 加载包
library(Seurat)
library(Matrix)
library(mclust)
library(parallel)
library(ggplot2)

# 加载修复版函数
source("cellcharter_R_parallel.R")

# =====================
# 1. 加载你的原始空间数据
# =====================
message("=== 1. 加载原始 mPT 对象 ===")
mPT <- readRDS("mPT_detailed_new.rds")  # 改成你的路径


# =====================
# 2. 运行 CellCharter 聚合邻居
# =====================
message("=== 2. 运行 CellCharter 聚合邻居 ===")
out_dir <- "./cellcharter_results_autok"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
idir <- file.path(out_dir, "intermediates")
if (!dir.exists(idir)) dir.create(idir, recursive = TRUE)

mPT <- aggregate_neighbors(mPT, n_layers = 3, use_rep = "pca",
                          x_coord_col = "CenterX_global_px", y_coord_col = "CenterY_global_px",
                          n_cores = 20, save_intermediates = FALSE,
                          intermediates_dir = idir)


# =====================
# 3. 运行自动选择 k 值（多个范围）
# =====================
message("=== 3. 运行自动选择 k 值 ===")

# 定义单一 k 值范围
k_range <- 2:25
range_name <- "2-25"

# 存储结果
all_results <- list()

message(paste("=== 分析 k 范围:", range_name, "==="))
range_out_dir <- file.path(out_dir, paste0("autok_", range_name))

# 检查是否已经存在结果
result_file <- file.path(range_out_dir, "autok_results.rds")
if (file.exists(result_file)) {
  message(paste("发现已存在的结果:", result_file))
  message("加载现有结果，跳过重新计算...")
  autok_results <- readRDS(result_file)
} else {
  message("未发现现有结果，开始计算...")
  # 运行自动选择 k 值
  autok_results <- cluster_autok(mPT, n_clusters_range = k_range, max_runs = 10, 
                               n_cores = 20, out_dir = range_out_dir,
                               sample_size = 0.1)
}

# 存储结果，包括 k 范围
all_results[[range_name]] <- list(
  best_k = autok_results$best_k,
  best_k_candidates = autok_results$best_k_candidates,
  stability_mean = autok_results$stability_mean,
  k_range = k_range
)


# =====================
# 4. 合并结果并绘制综合稳定性曲线
# =====================
message("=== 4. 绘制综合稳定性曲线 ===")

# 准备数据
all_stability <- data.frame()
for (range_name in names(all_results)) {
  result <- all_results[[range_name]]
  k_range <- result$k_range[-1]  # 去掉第一个值，因为稳定性是 k 和 k+1 之间的
  stability <- result$stability_mean
  
  # 确保长度匹配
  if (length(k_range) == length(stability)) {
    temp_df <- data.frame(
      k = k_range,
      stability = stability,
      range = range_name
    )
    all_stability <- rbind(all_stability, temp_df)
  }
}

# 绘制综合稳定性曲线
pdf(file.path(out_dir, "stability_curve_combined.pdf"), width = 12, height = 6)
ggplot(all_stability, aes(x = k, y = stability, color = range)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  labs(
    title = "Cluster stability across k values",
    x = "Number of clusters (k)",
    y = "Stability score",
    color = "k range"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )
dev.off()

# 找到全局最佳 k 值
global_best_k <- all_stability[which.max(all_stability$stability), "k"]
message(paste("全局最佳 k 值:", global_best_k))

# 保存综合结果
saveRDS(list(
  all_results = all_results,
  all_stability = all_stability,
  global_best_k = global_best_k
), file.path(out_dir, "autok_combined_results.rds"))

message("✅ 最佳 k 值选择完成！")
message("综合稳定性曲线已保存到: stability_curve_combined.pdf")
message("最佳 k 值: ", global_best_k)