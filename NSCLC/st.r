#10× Visium
setwd("E:/project/nsclc/GSE189487_RAW")
list.files()

library(Seurat)
library(Matrix)
library(dplyr)

data_dir <- "TD1/"
expr <- ReadMtx(
    mtx = file.path(data_dir, "matrix.mtx.gz"),
    features = file.path(data_dir, "features.tsv.gz"),
    cells = file.path(data_dir, "barcodes.tsv.gz")
)
coords <- read.csv(file.path(data_dir, "tissue_positions_list.csv.gz"), header=FALSE)
colnames(coords) <- c("barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres")
rownames(coords) <- coords$barcode
coords <- coords[, c("pxl_row_in_fullres", "pxl_col_in_fullres")]

seu <- CreateSeuratObject(counts = expr, assay = "Spatial")
common_cells <- intersect(colnames(seu), rownames(coords))
seu <- subset(seu, cells = common_cells)
coords <- coords[common_cells, , drop = FALSE]
seu[["image"]] <- new(Class = "SlideSeq", coordinates = coords, assay = "Spatial")

seu <- SCTransform(seu, assay = "Spatial", verbose = FALSE)

seu <- RunPCA(seu, verbose = FALSE)

seu <- FindNeighbors(seu, dims = 1:30)
seu <- FindClusters(seu, resolution = 0.5)

seu <- RunUMAP(seu, dims = 1:30)
SpatialPlot(seu, group.by = "seurat_clusters", image.alpha = 0)

markers <- FindAllMarkers(seu,
                          only.pos = TRUE,
                          test.use = "wilcox",
                          min.pct = 0.1,
                          logfc.threshold = 0.25)

marker_genes <- markers %>%
  filter(p_val_adj < 0.05, avg_log2FC > 0.25)
significant_markers <- subset(markers, p_val_adj < 0.05)
significant_markers <- significant_markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.csv(significant_markers,"marker_top.csv")

identity_mapping <- c(
    "0" = "Cilited cell",
    "1" = "Alveolar cell",
    "2" = "Fibroblast",
    "3" = "Alveolar Macrophage",
    "4" = "B cell",
    "5" = "Epithelial cell",
    "6" = "Cilited cell",
    "7" = "Secretory club cell",
    "8" = "Plasma"
)

cell_type <- identity_mapping[seu@meta.data$seurat_clusters]
seu@meta.data$cell_type <- cell_type

SpatialPlot(seu, group.by = "cell_type", image.alpha = 0, pt.size = 2)

#互作 
library(deldir)
library(RANN)      # 最近邻查询（用于精确匹配端点）
library(igraph)
library(pheatmap)
library(foreach)
library(doParallel)

coords <- seu[["image"]]@coordinates
# 确认行名和 meta 的顺序一致（若不一致，我们会用索引对应）
head(coords)
n_cells <- nrow(coords)
barcodes <- rownames(coords)

# 取细胞标签（factor）
labels <- as.character(seu@meta.data$cell_type)
names(labels) <- rownames(seu@meta.data)  # 保证名字为 barcode

# 对齐 coords 与 labels（取交集）
common <- intersect(rownames(coords), names(labels))
coords <- coords[common, , drop = FALSE]
labels <- labels[common]
barcodes <- common

xy <- coords[,1:2]
colnames(xy) <- c("x","y")
deld <- deldir(xy$x, xy$y)
segs <- deld$delsgs  # 三角边段，包含 x1,y1,x2,y2

# 为每个端点找到对应的原始点索引（用 RANN 的最近邻，避免小数误差）
pts_matrix <- as.matrix(xy)
# 构建查询表：把 segs 的端点连成矩阵
endpts <- rbind(
  data.frame(x=segs$x1, y=segs$y1),
  data.frame(x=segs$x2, y=segs$y2)
)
nn <- RANN::nn2(pts_matrix, query = as.matrix(endpts), k = 1)
idxs <- nn$nn.idx[,1]

# 组成边列表（每一对端点）
nseg <- nrow(segs)
edge_idx1 <- idxs[1:nseg]
edge_idx2 <- idxs[(nseg+1):(2*nseg)]

# 转换为条目对（行索引）
edges <- cbind(edge_idx1, edge_idx2)
# 去除自环与重复（保持 i < j）
edges <- edges[edges[,1] != edges[,2], , drop=FALSE]
edges <- t(apply(edges, 1, function(x) sort(x)))
edges <- unique(edges)
edges_df <- data.frame(from = edges[,1], to = edges[,2])

types <- sort(unique(labels))
K <- length(types)
# 将索引映射到类型
type_by_index <- labels  # 名字为 barcode，各元素为类型
index_to_type <- type_by_index[barcodes] # 确保顺序与 coords/barcodes 一致

# For edges_df, get types:
t1 <- index_to_type[edges_df$from]
t2 <- index_to_type[edges_df$to]

# 构建矩阵
mat_obs <- matrix(0, nrow = K, ncol = K, dimnames = list(types, types))
for(i in seq_along(t1)){
  a <- t1[i]; b <- t2[i]
  mat_obs[a,b] <- mat_obs[a,b] + 1
  mat_obs[b,a] <- mat_obs[b,a] + 1  # 对称计数（若你想只计一次可去掉这一行）
}
# 如果希望对角线表示 A-A 边数量（现在对称计两次），可以除以 2 或用不同计法
# 这里 mat_obs 为对称矩阵，A-B 在两个方向上均计数（常见做法之一）
mat_obs

nperm <- 1000   # 可按需修改
ncores <- parallel::detectCores() - 1
ncores <- max(1, ncores)
cl <- makeCluster(ncores)
registerDoParallel(cl)

# 事先准备 edges 的索引，避免在循环中重复计算
from_idx <- edges_df$from
to_idx <- edges_df$to

perm_counts <- foreach(p = 1:nperm, .packages = c(), .combine = rbind) %dopar% {
    set.seed(p + 12345)
    perm_labels <- sample(index_to_type, length(index_to_type), replace = FALSE)

    pt1 <- perm_labels[from_idx]
    pt2 <- perm_labels[to_idx]

    # ✅ 修复点：加 dimnames
    mat <- matrix(0, nrow = K, ncol = K, dimnames = list(types, types))

    for(i in seq_along(pt1)){
        a <- pt1[i]; b <- pt2[i]
        mat[a,b] <- mat[a,b] + 1
        mat[b,a] <- mat[b,a] + 1
    }
    as.vector(mat)
}

stopCluster(cl)
dim(perm_counts)  # nperm x (K*K)

# 观测向量（与 perm_counts 展开方式一致）
obs_vec <- as.vector(mat_obs)

mu_rand <- colMeans(perm_counts)
sd_rand <- apply(perm_counts, 2, sd)

# z score
z_vec <- (obs_vec - mu_rand) / (sd_rand + 1e-8)
# 经验 p 值（双边）: 计算 perm 中有多少次 >= obs 或 <= obs
p_emp <- sapply(seq_along(obs_vec), function(i){
  perm_i <- perm_counts[, i]
  obs_i <- obs_vec[i]
  # 双侧 p 值（包括等于）
  p = (sum(perm_i >= obs_i) + sum(perm_i <= obs_i)) / (2 * nperm)  # 这是不太常规的，改成两边更严谨:
  # 更常用：双侧 p = (count(|perm - mu| >= |obs - mu|) + 1) / (nperm + 1)
  mu_i <- mu_rand[i]
  p2 = (sum(abs(perm_i - mu_i) >= abs(obs_i - mu_i)) + 1) / (nperm + 1)
  p2
})

# 转回矩阵形式
mat_mu <- matrix(mu_rand, nrow = K, ncol = K, dimnames = list(types, types))
mat_sd <- matrix(sd_rand, nrow = K, ncol = K, dimnames = list(types, types))
mat_z <- matrix(z_vec, nrow = K, ncol = K, dimnames = list(types, types))
mat_p <- matrix(p_emp, nrow = K, ncol = K, dimnames = list(types, types))

# z-score 热图（超出 0 为富集）
pheatmap(mat_z,
         cluster_rows = TRUE, cluster_cols = TRUE,
         main = "Z-score of contact enrichment (Obs vs Random)",
         fontsize = 10)

# p-value 热图（可用 -log10）
neglogp <- -log10(mat_p)
neglogp[is.infinite(neglogp)] <- max(neglogp[is.finite(neglogp)], na.rm=TRUE)
pheatmap(neglogp,
         cluster_rows = TRUE, cluster_cols = TRUE,
         main = "-log10(empirical p) of contact enrichment",
         fontsize = 10)

#assortativity
types <- sort(unique(labels))
K <- length(types)

# 初始化混合矩阵
mix_mat <- matrix(0, nrow = K, ncol = K, dimnames = list(types, types))

# edges_df 中 from/to 为行号或 barcode，对应 labels
t1 <- labels[edges_df$from]
t2 <- labels[edges_df$to]

for(i in seq_along(t1)){
  a <- t1[i]; b <- t2[i]
  mix_mat[a,b] <- mix_mat[a,b] + 1
  mix_mat[b,a] <- mix_mat[b,a] + 1  # 对称
}

# 可选：除以 2 保持对称，只计一次
mix_mat <- mix_mat / 2

# 归一化
e <- mix_mat / sum(mix_mat)
a <- rowSums(e)
b <- colSums(e)

denom <- 1 - sum(a*b)
r_manual <- ifelse(denom == 0, NA, (sum(diag(e)) - sum(a*b)) / denom)
r_manual

#each cell type
# mix_mat 已经是对称矩阵，元素 = 边数
e <- mix_mat / sum(mix_mat)   # 归一化成概率矩阵

a <- rowSums(e)
b <- colSums(e)

r_local <- numeric(length = nrow(e))
names(r_local) <- rownames(e)

for(i in seq_along(r_local)){
  denom <- 1 - a[i]*b[i]
  r_local[i] <- ifelse(denom==0, NA, (e[i,i] - a[i]*b[i]) / denom)
}

r_local

library(ggplot2)

df_r <- data.frame(
  cell_type = names(r_local),
  assortativity = r_local
)

ggplot(df_r, aes(x = cell_type, y = assortativity, fill = cell_type)) +
  geom_bar(stat="identity") +
  ylim(-1,1) +
  geom_hline(yintercept = 0, linetype="dashed") +
  ylab("Local assortativity") +
  ggtitle("Assortativity per cell type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=45, hjust=1))

#cell type communication
types <- rownames(mix_mat)
K <- length(types)
e <- mix_mat / sum(mix_mat)   # 归一化
a <- rowSums(e)
b <- colSums(e)

r_pair <- matrix(0, nrow = K, ncol = K, dimnames = list(types, types))

for(i in 1:K){
  for(j in 1:K){
    denom <- 1 - a[i]*b[j]
    r_pair[i,j] <- ifelse(denom==0, NA, (e[i,j] - a[i]*b[j]) / denom)
  }
}

r_pair

library(pheatmap)

pheatmap(r_pair,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Pairwise assortativity (r_AB)",
         color = colorRampPalette(c("blue","white","red"))(50),
         fontsize = 10)


#GeoMx DSP
BiocManager::install(c("NanoStringNCTools", "limma", "edgeR"))
remotes::install_github("Nanostring-Biostats/GeoMxTools", dependencies = TRUE)

library(GeomxTools)
library(NanoStringNCTools)
library(limma)
library(edgeR)

dcc_gz <- list.files(pattern = "\\.dcc\\.gz$", full.names = TRUE)
sapply(dcc_gz, R.utils::gunzip, overwrite = TRUE)
dcc_files <- list.files(pattern = "\\.dcc$", full.names = TRUE)
raw_data <- readNanoStringGeoMxSet(dccFiles = dcc_files)
