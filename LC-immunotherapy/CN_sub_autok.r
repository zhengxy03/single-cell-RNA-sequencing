# ============================
# CellCharter 从头跑 → 智能多核版（自动选择 k 值）
# 最快速度 + 绝不报错
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
# 2. 运行 CellCharter（智能多核 + 自动选择 k 值）
# =====================
message("=== 2. 开始运行 CellCharter 自动选择 k 值 ===")

# 🔥 核心：这里填你最大核心数（比如 20）
# 软件内部会自动分配：该快的快，该稳的稳
mPT <- run_cellcharter_autok(
  n_layers = 3,
  n_clusters_range = 2:15,  # 搜索 k 值的范围
  n_cores = 20,        # 你想给多少就给多少！
  save_intermediates = FALSE,
  sample_size = 0.1    # 抽样比例，默认10%，可以根据需要调整
)


message("✅ 全部完成！")
