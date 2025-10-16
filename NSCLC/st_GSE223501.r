library(Seurat)
library(readr)  # 用于读取csv（比read.csv更快）

# 1. 导入表达矩阵（counts.csv）
setwd("E:/zxy/GSE223501_RAW/PA019")
counts <- read.csv("counts.csv.gz", row.names = 1)  # 假设第一列是基因名，设为行名
# 确保矩阵是 numeric 类型（非字符）
counts <- as.matrix(counts)
# 检查格式：行=基因，列=细胞
cat("表达矩阵维度（基因数×细胞数）：", dim(counts), "\n")
head(rownames(counts))  # 查看前几个基因名
head(colnames(counts))  # 查看前几个细胞ID

# 2. 导入坐标文件（coords.csv）
coords <- read.csv("coords.csv.gz", row.names = "cell_id")  # 以cell_id列为行名
# 检查坐标列名（通常是x、y，或row、col等）
cat("坐标列名：", colnames(coords), "\n")
# 保留坐标列（假设列名为x和y，根据实际修改）
coords <- coords[, c("x", "y")]  # 只保留x和y列
# 检查坐标行名（细胞ID）是否与表达矩阵列名匹配
common_cells <- intersect(colnames(counts), rownames(coords))
cat("共有的细胞数：", length(common_cells), "\n")  # 必须>0，否则无法匹配

# 3. 筛选表达矩阵和坐标，只保留共有的细胞
counts <- counts[, common_cells, drop = FALSE]
coords <- coords[common_cells, , drop = FALSE]

# 4. 创建Seurat对象
seu <- CreateSeuratObject(
  counts = counts,
  project = "SpatialData",
  assay = "Spatial"  # 空间数据通常用"Spatial" assay
)

# 5. 添加空间坐标到Seurat对象
# 方法1：适用于Visium-like数据（推荐）
seu[["image"]] <- CreateSpatialImage(
  image = matrix(0, nrow = 1, ncol = 1),  # 无需原始图像时用占位矩阵
  coordinates = coords  # 传入坐标数据框
)

# 方法2：若方法1报错（旧版本Seurat），直接添加到meta.data并手动关联
seu@meta.data <- cbind(seu@meta.data, coords)  # 坐标添加到元数据
# 后续绘图时可直接用ggplot2基于meta.data中的x和y绘图

# 6. 验证对象是否正确
print(seu)  # 查看对象信息（应包含细胞数、基因数）
head(seu@meta.data)  # 查看元数据（应包含x、y坐标）