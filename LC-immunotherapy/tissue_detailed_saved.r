library(Seurat)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(deldir)
library(doParallel)
library(foreach)

# ===================== 数据加载 =====================
obj <- readRDS("YA2025263-1_fin.rds")
Malignant <- readRDS("malignant_anno.rds")
fib <- readRDS("fib_anno_new.rds")
macro <- readRDS("macro_anno_new.rds")
t_cells <- readRDS("t_anno.rds")
DC <- readRDS("DC_anno.rds")
B <- readRDS("B_anno.rds")
# ===================== 定义样本与组织的对应关系 =====================
sample_tissue <- data.frame(
  sample = c("A1", "A2", "A3", "A4", "A5", 
             "B1", "B2", "B3", "B4", "B5",
             "C1", "C2", "C3", "C4", "C5",
             "D1", "D2", "D3", "D4", "D5"),
  tissue = c("nmPT", "negLN", "nmPT", "negLN", "nmPT",
             "negLN", "nmPT", "negLN", "mPT", "negLN",
             "metLN", "mPT", "negLN", "metLN", "mPT",
             "negLN", "metLN", "mPT", "negLN", "metLN")
)
obj$tissue <- sample_tissue$tissue[match(obj$sample, sample_tissue$sample)]

# ===================== 细胞类型注释 =====================
obj$CellType <- recode(obj$CellType,
                       #"Malignant cells" = "unknown",
                       "Basal cells" = "Malignant cells"
)

# 检查
table(obj$CellType)

obj@meta.data$sub_cell_type <- NA
common_cells <- intersect(rownames(Malignant@meta.data), rownames(obj@meta.data))
obj@meta.data[common_cells, "sub_cell_type"] <- Malignant@meta.data[common_cells, "sub_cell_type"]

common_cells <- intersect(rownames(fib@meta.data), rownames(obj@meta.data))
obj@meta.data[common_cells, "sub_cell_type"] <- fib@meta.data[common_cells, "sub_cell_type"]

common_cells <- intersect(rownames(macro@meta.data), rownames(obj@meta.data))
obj@meta.data[common_cells, "sub_cell_type"] <- macro@meta.data[common_cells, "sub_cell_type"]

common_cells <- intersect(rownames(t_cells@meta.data), rownames(obj@meta.data))
obj@meta.data[common_cells, "sub_cell_type"] <- t_cells@meta.data[common_cells, "sub_cell_type"]

common_cells <- intersect(rownames(DC@meta.data), rownames(obj@meta.data))
obj@meta.data[common_cells, "sub_cell_type"] <- DC@meta.data[common_cells, "sub_cell_type"]

common_cells <- intersect(rownames(B@meta.data), rownames(obj@meta.data))
obj@meta.data[common_cells, "sub_cell_type"] <- B@meta.data[common_cells, "sub_cell_type"]

cat("成功注释的细胞数：", sum(!is.na(obj@meta.data$sub_cell_type)), "\n")
cat("obj中sub_cell_type列的唯一值：", paste(unique(obj@meta.data$sub_cell_type), collapse = ", "), "\n")


print(class(obj@meta.data$CellType))  # 查看类型（应该是factor）
cell_type_levels <- levels(obj@meta.data$CellType)  # 提取因子水平（数字→名称映射）
cat("CellType因子水平（数字→名称）：\n")
print(cell_type_levels)

# 1. 重置detailed列（先清空错误赋值）
obj@meta.data$detailed <- NA

# 2. 优先级1：填充sub_cell_type（非NA的细胞，保留原名称）
sub_cell_idx <- !is.na(obj@meta.data$sub_cell_type)
obj@meta.data$detailed[sub_cell_idx] <- as.character(obj@meta.data$sub_cell_type[sub_cell_idx])
cat("优先级1（sub_cell_type）填充：", sum(sub_cell_idx), "个细胞\n")



# 4. 优先级3：填充CellType（剩余细胞，关键：转字符+用因子水平映射真实名称）
cell_type_idx <- is.na(obj@meta.data$detailed) & !is.na(obj@meta.data$CellType)
# 方法：如果CellType是因子，用levels映射；否则直接转字符
if (is.factor(obj@meta.data$CellType)) {
  # 因子类型：用因子水平把数字编码转成真实名称
  obj@meta.data$detailed[cell_type_idx] <- cell_type_levels[as.integer(obj@meta.data$CellType[cell_type_idx])]
} else {
  # 字符类型：直接转字符
  obj@meta.data$detailed[cell_type_idx] <- as.character(obj@meta.data$CellType[cell_type_idx])
}
cat("优先级3（CellType）填充：", sum(cell_type_idx), "个细胞\n")

# 验证修复结果
cat("\n=== 修复后detailed列前10行 ===")
print(head(obj@meta.data[, c("sub_cell_type", "CellType", "detailed")], 10))
cat("\n=== 修复后detailed列前15个类型 ===")
print(head(sort(table(obj@meta.data$detailed), decreasing = TRUE), 15))

remove_types <- c("B cells", "Dendritic cells", "Germinal center B cells", 
                  "immature B cells", "Malignant cells", "T cells","Plasmacytoid dendritic cells")

# 方法：保留不在移除列表中的细胞
cells_keep <- which(!(obj@meta.data$detailed %in% remove_types))
obj <- obj[, cells_keep]

# 检查结果
ncol(obj)
table(obj@meta.data$detailed)

mPT <- subset(obj,subset=tissue=="mPT")
saveRDS(mPT,file="mPT_detailed.rds")
saveRDS(obj,file="obj_filtered_detailed.rds")
nmPT <- subset(obj,subset=tissue=="nmPT")
saveRDS(nmPT,file="nmPT_detailed.rds")
