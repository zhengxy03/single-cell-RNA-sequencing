# ============================
# CellCharter 固定 k 值分析
# 使用设定的 k 值进行完整分析
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
# 2. 运行 CellCharter（使用固定 k 值）
# =====================
message("=== 2. 开始运行 CellCharter 固定 k 值分析 ===")

# 🔥 核心：这里填你设定的 k 值
# 软件内部会自动分配：该快的快，该稳的稳
mPT <- run_cellcharter_mPT(
  n_layers = 3,
  n_clusters = 19,        # 🔥 最佳k值 = 19
  n_cores = 20,           # 你想给多少就给多少！
  save_intermediates = FALSE
)


message("✅ 全部完成！")

# 保存最终结果
saveRDS(mPT, "mPT_cellcharter_FINAL.rds")
message("结果已保存到: mPT_cellcharter_FINAL.rds")
