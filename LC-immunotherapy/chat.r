obj <- readRDS("YA2025263-1_fin.rds")
library(Seurat)

library(ggplot2)
library(dplyr)
Malignant <- readRDS("malignant_anno.rds")

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

obj$CellType <- recode(obj$CellType,
                       "Malignant cells" = "unknown",
                       "Basal cells" = "Malignant cells"
)

# 检查
table(obj$CellType)

obj@meta.data$sub_cell_type <- NA
common_cells <- intersect(rownames(Malignant@meta.data), rownames(obj@meta.data))
obj@meta.data[common_cells, "sub_cell_type"] <- Malignant@meta.data[common_cells, "sub_cell_type"]

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



coords <- obj@meta.data[, c("CenterX_global_px", "CenterY_global_px")]
head(coords)

# 获取细胞barcodes（就是行名）
barcodes <- rownames(coords)

labels <- as.character(obj@meta.data$detailed)  # 或者用 obj@meta.data$sub_cell_type
names(labels) <- barcodes

# 进行匹配（坐标和标签的行名应该一致，这里做一下检查）
common <- intersect(rownames(coords), names(labels))
print(paste("共同barcodes数量:", length(common)))




# 检查细胞类型分布
print("细胞类型分布:")
print(table(labels))

# 处理数据
coords <- coords[common, , drop = FALSE]
labels <- labels[common]
barcodes <- common

# 移除NA值的细胞
valid_cells <- !is.na(labels)
coords <- coords[valid_cells, , drop = FALSE]
labels <- labels[valid_cells]
barcodes <- barcodes[valid_cells]

print(paste("有效细胞数量:", length(labels)))
print("最终细胞类型分布:")
print(table(labels))

# 继续处理...
xy <- coords
colnames(xy) <- c("x","y")

print("开始计算Delaunay三角剖分...")
library(deldir)

x <- xy[,1]
y <- xy[,2]
#rw <- c(min(x), max(x), min(y), max(y))

#deld <- deldir(x, y, rw = rw)
set.seed(42)
sample_size <- min(700000, length(x))
sample_idx <- sample(length(x), sample_size)

x_sub <- x[sample_idx]
y_sub <- y[sample_idx]
labels_sub <- labels[sample_idx]

deld <- deldir(x_sub, y_sub, rw = c(range(x_sub), range(y_sub)))



segs <- deld$delsgs

print(paste("生成三角边数量:", nrow(segs)))

# 使用索引
edges <- cbind(segs$ind1, segs$ind2)

# 去除自环与重复
edges <- edges[edges[,1] != edges[,2], , drop=FALSE]
edges <- t(apply(edges, 1, function(x) sort(x)))
edges <- unique(edges)
edges_df <- data.frame(from = edges[,1], to = edges[,2])

print(paste("最终边数量:", nrow(edges_df)))

# 获取所有细胞类型
types <- sort(unique(labels))
K <- length(types)

print(paste("细胞类型数量:", K))
print("细胞类型:")
print(types)

# 将索引映射到类型
type_by_index <- labels
index_to_type <- type_by_index

# For edges_df, get types:
t1 <- index_to_type[edges_df$from]
t2 <- index_to_type[edges_df$to]

# 构建矩阵
mat_obs <- matrix(0, nrow = K, ncol = K, dimnames = list(types, types))

for(i in seq_along(t1)){
    a <- t1[i]
    b <- t2[i]
    mat_obs[a,b] <- mat_obs[a,b] + 1
    mat_obs[b,a] <- mat_obs[b,a] + 1
}

print("观察到的细胞-细胞接触矩阵：")
print(mat_obs)

# 排列检验
library(doParallel)
library(foreach)

nperm <- 1000
ncores <- parallel::detectCores() - 1
ncores <- max(1, ncores)

print(paste("使用", ncores, "个核心进行排列检验"))

cl <- makeCluster(ncores)
registerDoParallel(cl)

from_idx <- edges_df$from
to_idx <- edges_df$to

perm_counts <- foreach(p = 1:nperm, .packages = c(), .combine = rbind) %dopar% {
    set.seed(p + 12345)
    perm_labels <- sample(index_to_type, length(index_to_type), replace = FALSE)
    
    pt1 <- perm_labels[from_idx]
    pt2 <- perm_labels[to_idx]
    
    mat <- matrix(0, nrow = K, ncol = K, dimnames = list(types, types))
    
    for(i in seq_along(pt1)){
        a <- pt1[i]
        b <- pt2[i]
        mat[a,b] <- mat[a,b] + 1
        mat[b,a] <- mat[b,a] + 1
    }
    as.vector(mat)
}

stopCluster(cl)

print("排列检验完成")

# 继续剩余的分析代码...
obs_vec <- as.vector(mat_obs)
mu_rand <- colMeans(perm_counts)
sd_rand <- apply(perm_counts, 2, sd)

# z score
z_vec <- (obs_vec - mu_rand) / (sd_rand + 1e-8)

# 经验 p 值
p_emp <- sapply(seq_along(obs_vec), function(i){
  perm_i <- perm_counts[, i]
  obs_i <- obs_vec[i]
  mu_i <- mu_rand[i]
  p_val = (sum(abs(perm_i - mu_i) >= abs(obs_i - mu_i)) + 1) / (nperm + 1)
  p_val
})

# 转回矩阵形式
mat_mu <- matrix(mu_rand, nrow = K, ncol = K, dimnames = list(types, types))
mat_sd <- matrix(sd_rand, nrow = K, ncol = K, dimnames = list(types, types))
mat_z <- matrix(z_vec, nrow = K, ncol = K, dimnames = list(types, types))
mat_p <- matrix(p_emp, nrow = K, ncol = K, dimnames = list(types, types))

# z-score 热图
library(pheatmap)
p <- pheatmap(mat_z,
         cluster_rows = TRUE, 
         cluster_cols = TRUE,
         main = "Z-score of contact enrichment (Obs vs Random)",
         fontsize = 10)
ggsave("chat_heatmap.png",plot=p)
print("分析完成！")

# 保存结果
save(mat_obs, mat_z, mat_p, mat_mu, file = "spatial_contact_analysis_results.RData")

