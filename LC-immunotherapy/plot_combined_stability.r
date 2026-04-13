# ============================
# 绘制合并的稳定性曲线
# 从已生成的rds文件中读取数据
# ============================

library(ggplot2)

# 设置输出目录
out_dir <- "./cellcharter_results_autok"

# 读取综合结果
combined_results <- readRDS(file.path(out_dir, "autok_combined_results.rds"))
all_stability <- combined_results$all_stability

# 或者如果综合结果不存在，可以从单独的结果文件中读取
# 读取2-15范围的结果
result_2_15 <- readRDS(file.path(out_dir, "autok_2-15/autok_results.rds"))
stability_2_15 <- data.frame(
  k = 3:15,  # 2-15范围的稳定性对应k=3到15
  stability = result_2_15$stability_mean,
  range = "2-15"
)

# 读取16-25范围的结果
result_16_25 <- readRDS(file.path(out_dir, "autok_16-25/autok_results.rds"))
stability_16_25 <- data.frame(
  k = 17:25,  # 16-25范围的稳定性对应k=17到25
  stability = result_16_25$stability_mean,
  range = "16-25"
)

# 合并数据
all_stability <- rbind(stability_2_15, stability_16_25)

# 绘制合并的稳定性曲线（与单独的画法一致）
pdf(file.path(out_dir, "stability_curve_combined.pdf"), width = 10, height = 6)

plot(all_stability$k, all_stability$stability, type = "b", 
     xlab = "Number of clusters (k)", 
     ylab = "Stability score", 
     main = "Cluster stability across k values",
     pch = 16, col = "black")

# 找到局部最大值
find_peaks <- function(x) {
  peaks <- c()
  for (i in 2:(length(x)-1)) {
    if (x[i] > x[i-1] && x[i] > x[i+1]) {
      peaks <- c(peaks, i)
    }
  }
  return(peaks)
}

# 只标记全局最大值（最佳k值）为一个红点
best_k_idx <- which.max(all_stability$stability)
points(all_stability$k[best_k_idx], all_stability$stability[best_k_idx], 
       col = "red", pch = 16, cex = 1.5)

# 标记最佳k值
best_k <- all_stability$k[which.max(all_stability$stability)]
abline(v = best_k, col = "blue", lty = 2)

# 添加图例
legend("topright", 
       c("Stability", "Local maxima", "Best k"), 
       col = c("black", "red", "blue"), 
       pch = c(16, 16, NA), 
       lty = c(1, NA, 2))

dev.off()

message("✅ 合并稳定性曲线已保存到: stability_curve_combined.pdf")
message("最佳 k 值: ", best_k)