# download space ranger
```
curl -o spaceranger-4.0.1.tar.gz "https://cf.10xgenomics.com/releases/spatial-exp/spaceranger-4.0.1.tar.gz?Expires=1760836993&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=TZ-xX2mvIbCK768RyFYO93ef5QLL4UiZD-SZyYrPVLJyXo4I1IZhiJEb9kzgzzgx4mwSYZaptAmvus28pmkjfWRloTY1k8fO-TalrO7erewbnwfAMk2vMFi67PDdyJy0XJbNQHdOtOdFym-k3PdfKX9s4-hYr-M-VMow6Wj5Jhk0VMSnTLkr~SgLk7pC2u1C~jCI43wXaGCqMvI2u~p~JOWHFYmaLIZ0HSsk4CFnjgzqfOvwD46ZRujUApCSruZ6XGFbHWclfEq3p5gw-kj-zL-5LjDt0kEgKl-Jl~Det-TrKkpJB7mDINLDZAPxKlKUZFzvHK7~PvJ~yPgDHDWrfg__"
tar xvfz spaceranger-4.0.1.tar.gz
cd spaceranger-4.0.1/bin
export PATH="$(pwd):$PATH"
```
# download mouse reference genome
```
curl -O "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz"
tar xvfz refdata-gex-GRCh38-2024-A.tar.gz
```
# rename fastq file
e.g SampleName1_S00_L001_R1_001.fastq.gz<br>
    SampleName1_S00_L001_R2_001.fastq.gz<br>
    SampleName2_S01_L001_R1_001.fastq.gz<br>
    SampleName2_S01_L001_R2_001.fastq.gz<br>
# download Loupe Browser and modify the H&E images
# run spaceranger
```
spaceranger count --id=CRR857759 \
                   --transcriptome=/share/home/wangq/zxy/CRA012359/refgenome/refdata-gex-GRCh38-2024-A \
                   --fastqs=/share/home/wangq/zxy/CRA012359_1 \
                   --sample=CRR857759 \
                   --image=/share/home/wangq/zxy/CRA012359_1/lgin_escc.tif \
                   --localcores=24 \
                   --localmem=112 \
                   --create-bam=false \
                   --slide=V19L01-041 \
                   --slidefile=/share/home/wangq/zxy/CRA012359/V19L01-041.gpr \
                   --loupe-alignment=/share/home/wangq/zxy/CRA012359_1/V19L01-041-A1_2.json \
                   --area=A1
```
# downstream analysis
```R
#导入
library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)

setwd("/share/home/wangq/zxy/CRA012359_1")
counts_matrix1 <- Read10X(data.dir = "./CRR857758/outs/filtered_feature_bc_matrix")
seu1 <- CreateSeuratObject(counts = counts_matrix1, assay = "Spatial",  project = "Sample1")
coords1 <- read.csv(file.path("./CRR857758/outs/spatial/tissue_positions.csv"), header=FALSE)

counts_matrix2 <- Read10X(data.dir = "./CRR857759/outs/filtered_feature_bc_matrix")
seu2 <- CreateSeuratObject(counts = counts_matrix2, assay = "Spatial",  project = "Sample2")
coords2 <- read.csv(file.path("./CRR857759/outs/spatial/tissue_positions.csv"), header=FALSE)
#修改原始文件
coords1 <- coords1[-1, ]
colnames(coords1) <- c("barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col")

common_cells <- intersect(colnames(seu1), coords1$barcode)
seu1 <- subset(seu1, cells = common_cells)
coords1 <- coords1[match(common_cells, coords1$barcode), ]

coords_matrix1 <- as.matrix(coords1[, c("pxl_row", "pxl_col")])
rownames(coords_matrix1) <- coords1$barcode
coords_df1 <- as.data.frame(coords_matrix1)

seu1[["image"]] <- new(Class = "SlideSeq", coordinates = coords_df1, assay = "Spatial")

coords2 <- coords2[-1, ]
colnames(coords2) <- c("barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col")

common_cells <- intersect(colnames(seu2), coords2$barcode)
seu2 <- subset(seu2, cells = common_cells)
coords2 <- coords2[match(common_cells, coords2$barcode), ]

coords_matrix2 <- as.matrix(coords2[, c("pxl_row", "pxl_col")])
rownames(coords_matrix2) <- coords2$barcode
coords_df2 <- as.data.frame(coords_matrix2)


seu2[["image"]] <- new(Class = "SlideSeq", coordinates = coords_df2, assay = "Spatial")







seu1 <- subset(seu1, subset = nCount_Spatial > 0)
seu1 <- SCTransform(seu1, assay = "Spatial", verbose = FALSE)
seu1 <- RunPCA(seu1, verbose = FALSE)

seu1 <- FindNeighbors(seu1, dims = 1:30)
seu1 <- FindClusters(seu1, resolution = 0.5)

seu1 <- RunUMAP(seu1, dims = 1:30)
SpatialPlot(seu1, group.by = "seurat_clusters", image.alpha = 0)


```
```
library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
#import data
setwd("F:/CRA012359/")
# Load the expression data
expr.data1 <- Seurat::Read10X_h5(filename =  "./CRR857758/filtered_feature_bc_matrix.h5" )
seu1 <- Seurat::CreateSeuratObject(counts = expr.data1, project = 'CRR857758', assay = 'Spatial')
seu1$slice <- 1
seu1$region <- 'nor_inf'
# Load the image data
img1 <- Seurat::Read10X_Image(image.dir = './CRR857758/spatial')
Seurat::DefaultAssay(object = img1) <- 'Spatial'
img1 <- img1[colnames(x = seu1)]
seu1[['image']] <- img1
SpatialFeaturePlot(seu1, features = "nCount_Spatial")

expr.data2 <- Seurat::Read10X_h5(filename =  "./CRR857759/filtered_feature_bc_matrix.h5" )
seu2 <- Seurat::CreateSeuratObject(counts = expr.data2, project = 'CRR8577599', assay = 'Spatial')
seu2$slice <- 1
seu2$region <- 'lgin_escc'
# Load the image data
img2 <- Seurat::Read10X_Image(image.dir = './CRR857759/spatial')
Seurat::DefaultAssay(object = img2) <- 'Spatial'
img2 <- img2[colnames(x = seu2)]
seu2[['image']] <- img2
SpatialFeaturePlot(seu2, features = "nCount_Spatial")
```
# spots distribution
```
#seu1
seu1 <- SCTransform(seu1, assay = "Spatial", verbose = FALSE)
seu1 <- RunPCA(seu1, verbose = FALSE)

seu1 <- FindNeighbors(seu1, dims = 1:20)
seu1 <- FindClusters(seu1, resolution = 0.5)

seu1 <- RunUMAP(seu1, dims = 1:20)
SpatialPlot(seu1, group.by = "seurat_clusters", image.alpha = 1, pt.size.factor = 3)

#seu2
seu2 <- SCTransform(seu2, assay = "Spatial", verbose = FALSE)
seu2 <- RunPCA(seu2, verbose = FALSE)

seu2 <- FindNeighbors(seu2, dims = 1:20)
seu2 <- FindClusters(seu2, resolution = 0.2)

seu2 <- RunUMAP(seu2, dims = 1:20)
SpatialPlot(seu2, group.by = "seurat_clusters", image.alpha = 1)


merged_seu <- merge(seu1, seu2, add.cell.ids = c("S1", "S2"))
DefaultAssay(merged_seu) <- "Spatial" 
merged_seu <- JoinLayers(merged_seu)

merged_seu <- SCTransform(merged_seu, assay = "Spatial", verbose = FALSE)
merged_seu <- RunPCA(merged_seu, verbose = FALSE)

merged_seu <- FindNeighbors(merged_seu, dims = 1:20)
merged_seu <- FindClusters(merged_seu, resolution = 0.5)

merged_seu <- RunUMAP(merged_seu, dims = 1:20)
SpatialPlot(merged_seu, group.by = "seurat_clusters", image.alpha = 1, pt.size.factor = 3)
```
# cell annotation
```
markers_merged <- FindAllMarkers(merged_seu,
                          only.pos = TRUE,
                          test.use = "wilcox",
                          min.pct = 0.1,
                          logfc.threshold = 0.25)

significant_markers_merged <- subset(markers_merged, p_val_adj < 0.05)
significant_markers_merged <- significant_markers_merged %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.csv(significant_markers_merged,"marker_top_merged.csv")

markers_seu1 <- FindAllMarkers(seu1,
                          only.pos = TRUE,
                          test.use = "wilcox",
                          min.pct = 0.1,
                          logfc.threshold = 0.25)

significant_markers_seu1 <- subset(markers_seu1, p_val_adj < 0.05)
significant_markers_seu1 <- significant_markers_seu1 %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.csv(significant_markers_seu1,"marker_top_seu1.csv")


markers_seu2 <- FindAllMarkers(seu2,
                          only.pos = TRUE,
                          test.use = "wilcox",
                          min.pct = 0.1,
                          logfc.threshold = 0.25)

significant_markers_seu2 <- subset(markers_seu2, p_val_adj < 0.05)
significant_markers_seu2 <- significant_markers_seu2 %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.csv(significant_markers_seu2,"marker_top_seu2.csv")


identity_mapping <- c(
    "0" = "Epithelial-stromal crosstalk region",
    "1" = "Muscle cell region",
    "2" = "ESCC development region",
    "3" = "Smooth muscle cell region",
    "4" = "Transitional epithelial region"
)
cell_type <- identity_mapping[merged_seu@meta.data$seurat_clusters]
merged_seu@meta.data$cell_type <- cell_type


SpatialPlot(merged_seu, group.by = "cell_type", image.alpha =1, pt.size = 3, label.size = 8) 

DimPlot(merged_seu, reduction = "umap", group.by = "cell_type", pt.size = 1, label.size = 8)

prop <- read.csv("st_sum.csv")
npg_pal <- pal_npg()(6)
npg_extended <- colorRampPalette(npg_pal)(6)

p1 <- ggplot(prop, aes(x = factor(sample, levels = c("NOR","INF", "LGIN","HGIN","ESCC")), 
                       y = proportion, 
                       fill = cell_type)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = npg_extended) +
    labs(title = "Cell Type Composition by Sample",
         x = "Sample", 
         y = "Percentage (%)",
         fill = "Cell Type") +
    theme_classic() +  # 使用classic主题，自带坐标轴且无网格线
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p1)
```
# 互作 
```
library(deldir)
library(RANN)      # 最近邻查询（用于精确匹配端点）
library(igraph)
library(pheatmap)
library(foreach)
library(doParallel)


coords_list <- list()
for (img_name in image_names) {
    coords_list[[img_name]] <- GetTissueCoordinates(merged_seu, image = img_name)
}

# 访问特定图像的坐标
coords_image1 <- coords_list[[image_names[1]]]
coords_image2 <- coords_list[[image_names[2]]]

# 首先分析第一张图像
image_name1 <- image_names[1]
cat("=== 分析第一张图像:", image_name1, "===\n")

# 获取第一张图像的坐标
coords1 <- GetTissueCoordinates(merged_seu, image = image_name1)

# 确认行名和 meta 的顺序一致
head(coords1)
n_cells1 <- nrow(coords1)
barcodes1 <- rownames(coords1)

# 取细胞标签（factor）
labels1 <- as.character(merged_seu@meta.data$cell_type)
names(labels1) <- rownames(merged_seu@meta.data)

# 对齐 coords 与 labels（取交集）
common1 <- intersect(rownames(coords1), names(labels1))
coords1 <- coords1[common1, , drop = FALSE]
labels1 <- labels1[common1]
barcodes1 <- common1

xy1 <- coords1[,1:2]
colnames(xy1) <- c("x","y")
deld1 <- deldir(xy1$x, xy1$y)
segs1 <- deld1$delsgs

# 为每个端点找到对应的原始点索引
pts_matrix1 <- as.matrix(xy1)
endpts1 <- rbind(
  data.frame(x=segs1$x1, y=segs1$y1),
  data.frame(x=segs1$x2, y=segs1$y2)
)
nn1 <- RANN::nn2(pts_matrix1, query = as.matrix(endpts1), k = 1)
idxs1 <- nn1$nn.idx[,1]

# 组成边列表
nseg1 <- nrow(segs1)
edge_idx1_1 <- idxs1[1:nseg1]
edge_idx2_1 <- idxs1[(nseg1+1):(2*nseg1)]

edges1 <- cbind(edge_idx1_1, edge_idx2_1)
edges1 <- edges1[edges1[,1] != edges1[,2], , drop=FALSE]
edges1 <- t(apply(edges1, 1, function(x) sort(x)))
edges1 <- unique(edges1)
edges_df1 <- data.frame(from = edges1[,1], to = edges1[,2])

types1 <- sort(unique(labels1))
K1 <- length(types1)
type_by_index1 <- labels1
index_to_type1 <- type_by_index1[barcodes1]

# For edges_df, get types:
t1_1 <- index_to_type1[edges_df1$from]
t2_1 <- index_to_type1[edges_df1$to]

# 构建矩阵
mat_obs1 <- matrix(0, nrow = K1, ncol = K1, dimnames = list(types1, types1))
for(i in seq_along(t1_1)){
  a <- t1_1[i]; b <- t2_1[i]
  mat_obs1[a,b] = mat_obs1[a,b] + 1
  mat_obs1[b,a] = mat_obs1[b,a] + 1
}

# 置换检验
nperm <- 1000
ncores <- parallel::detectCores() - 1
ncores <- max(1, ncores)
cl <- makeCluster(ncores)
registerDoParallel(cl)

from_idx1 <- edges_df1$from
to_idx1 <- edges_df1$to

perm_counts1 <- foreach(p = 1:nperm, .packages = c(), .combine = rbind) %dopar% {
    set.seed(p + 12345)
    perm_labels <- sample(index_to_type1, length(index_to_type1), replace = FALSE)
    
    pt1 <- perm_labels[from_idx1]
    pt2 <- perm_labels[to_idx1]
    
    mat <- matrix(0, nrow = K1, ncol = K1, dimnames = list(types1, types1))
    
    for(i in seq_along(pt1)){
        a <- pt1[i]; b <- pt2[i]
        mat[a,b] = mat[a,b] + 1
        mat[b,a] = mat[b,a] + 1
    }
    as.vector(mat)
}

stopCluster(cl)

# 计算统计量
obs_vec1 <- as.vector(mat_obs1)
mu_rand1 <- colMeans(perm_counts1)
sd_rand1 <- apply(perm_counts1, 2, sd)

z_vec1 <- (obs_vec1 - mu_rand1) / (sd_rand1 + 1e-8)

p_emp1 <- sapply(seq_along(obs_vec1), function(i){
  perm_i <- perm_counts1[, i]
  obs_i <- obs_vec1[i]
  mu_i <- mu_rand1[i]
  p_val <- (sum(abs(perm_i - mu_i) >= abs(obs_i - mu_i)) + 1) / (nperm + 1)
  p_val
})

mat_mu1 <- matrix(mu_rand1, nrow = K1, ncol = K1, dimnames = list(types1, types1))
mat_sd1 <- matrix(sd_rand1, nrow = K1, ncol = K1, dimnames = list(types1, types1))
mat_z1 <- matrix(z_vec1, nrow = K1, ncol = K1, dimnames = list(types1, types1))
mat_p1 <- matrix(p_emp1, nrow = K1, ncol = K1, dimnames = list(types1, types1))

# 绘制第一张图像的热图
pheatmap(mat_z1,
         cluster_rows = TRUE, cluster_cols = TRUE,
         main = paste("Z-score of contact enrichment -", image_name1),
         fontsize = 10)

cat("=== 第一张图像分析完成 ===\n")
cat("细胞类型:", types1, "\n")
cat("细胞数量:", length(common1), "\n")
cat("边数量:", nrow(edges_df1), "\n")

# 现在分析第二张图像
image_name2 <- image_names[2]
cat("\n=== 分析第二张图像:", image_name2, "===\n")

# 获取第二张图像的坐标
coords2 <- GetTissueCoordinates(merged_seu, image = image_name2)

# 确认行名和 meta 的顺序一致
head(coords2)
n_cells2 <- nrow(coords2)
barcodes2 <- rownames(coords2)

# 取细胞标签（factor）
labels2 <- as.character(merged_seu@meta.data$cell_type)
names(labels2) <- rownames(merged_seu@meta.data)

# 对齐 coords 与 labels（取交集）
common2 <- intersect(rownames(coords2), names(labels2))
coords2 <- coords2[common2, , drop = FALSE]
labels2 <- labels2[common2]
barcodes2 <- common2

xy2 <- coords2[,1:2]
colnames(xy2) <- c("x","y")
deld2 <- deldir(xy2$x, xy2$y)
segs2 <- deld2$delsgs

# 为每个端点找到对应的原始点索引
pts_matrix2 <- as.matrix(xy2)
endpts2 <- rbind(
  data.frame(x=segs2$x1, y=segs2$y1),
  data.frame(x=segs2$x2, y=segs2$y2)
)
nn2 <- RANN::nn2(pts_matrix2, query = as.matrix(endpts2), k = 1)
idxs2 <- nn2$nn.idx[,1]

# 组成边列表
nseg2 <- nrow(segs2)
edge_idx1_2 <- idxs2[1:nseg2]
edge_idx2_2 <- idxs2[(nseg2+1):(2*nseg2)]

edges2 <- cbind(edge_idx1_2, edge_idx2_2)
edges2 <- edges2[edges2[,1] != edges2[,2], , drop=FALSE]
edges2 <- t(apply(edges2, 1, function(x) sort(x)))
edges2 <- unique(edges2)
edges_df2 <- data.frame(from = edges2[,1], to = edges2[,2])

types2 <- sort(unique(labels2))
K2 <- length(types2)
type_by_index2 <- labels2
index_to_type2 <- type_by_index2[barcodes2]

# For edges_df, get types:
t1_2 <- index_to_type2[edges_df2$from]
t2_2 <- index_to_type2[edges_df2$to]

# 构建矩阵
mat_obs2 <- matrix(0, nrow = K2, ncol = K2, dimnames = list(types2, types2))
for(i in seq_along(t1_2)){
  a <- t1_2[i]; b <- t2_2[i]
  mat_obs2[a,b] = mat_obs2[a,b] + 1
  mat_obs2[b,a] = mat_obs2[b,a] + 1
}

# 置换检验
cl <- makeCluster(ncores)
registerDoParallel(cl)

from_idx2 <- edges_df2$from
to_idx2 <- edges_df2$to

perm_counts2 <- foreach(p = 1:nperm, .packages = c(), .combine = rbind) %dopar% {
    set.seed(p + 12345)
    perm_labels <- sample(index_to_type2, length(index_to_type2), replace = FALSE)
    
    pt1 <- perm_labels[from_idx2]
    pt2 <- perm_labels[to_idx2]
    
    mat <- matrix(0, nrow = K2, ncol = K2, dimnames = list(types2, types2))
    
    for(i in seq_along(pt1)){
        a <- pt1[i]; b <- pt2[i]
        mat[a,b] = mat[a,b] + 1
        mat[b,a] = mat[b,a] + 1
    }
    as.vector(mat)
}

stopCluster(cl)

# 计算统计量
obs_vec2 <- as.vector(mat_obs2)
mu_rand2 <- colMeans(perm_counts2)
sd_rand2 <- apply(perm_counts2, 2, sd)

z_vec2 <- (obs_vec2 - mu_rand2) / (sd_rand2 + 1e-8)

p_emp2 <- sapply(seq_along(obs_vec2), function(i){
  perm_i <- perm_counts2[, i]
  obs_i <- obs_vec2[i]
  mu_i <- mu_rand2[i]
  p_val <- (sum(abs(perm_i - mu_i) >= abs(obs_i - mu_i)) + 1) / (nperm + 1)
  p_val
})

mat_mu2 <- matrix(mu_rand2, nrow = K2, ncol = K2, dimnames = list(types2, types2))
mat_sd2 <- matrix(sd_rand2, nrow = K2, ncol = K2, dimnames = list(types2, types2))
mat_z2 <- matrix(z_vec2, nrow = K2, ncol = K2, dimnames = list(types2, types2))
mat_p2 <- matrix(p_emp2, nrow = K2, ncol = K2, dimnames = list(types2, types2))

# 绘制第二张图像的热图
pheatmap(mat_z2,
         cluster_rows = TRUE, cluster_cols = TRUE,
         main = paste("Z-score of contact enrichment -", image_name2),
         fontsize = 10)

cat("=== 第二张图像分析完成 ===\n")
cat("细胞类型:", types2, "\n")
cat("细胞数量:", length(common2), "\n")
cat("边数量:", nrow(edges_df2), "\n")

# 比较两张图像的结果
cat("\n=== 两张图像比较 ===\n")
cat("图像1 -", image_name1, ": 细胞数 =", length(common1), ", 边数 =", nrow(edges_df1), "\n")
cat("图像2 -", image_name2, ": 细胞数 =", length(common2), ", 边数 =", nrow(edges_df2), "\n")

# 如果需要比较 z-score 差异
if (all(dim(mat_z1) == dim(mat_z2)) && all(rownames(mat_z1) == rownames(mat_z2))) {
    z_diff <- mat_z1 - mat_z2
    pheatmap(z_diff,
             cluster_rows = TRUE, cluster_cols = TRUE,
             main = "Z-score difference (Image1 - Image2)",
             fontsize = 10)
}


#两张图混合
# 获取两张图像的坐标并合并
coords1 <- GetTissueCoordinates(merged_seu, image = image_names[1])
coords2 <- GetTissueCoordinates(merged_seu, image = image_names[2])

# 为每个坐标添加样本标识
coords1$sample <- image_names[1]
coords2$sample <- image_names[2]

# 合并坐标
all_coords <- rbind(coords1, coords2)

# 确认行名和 meta 的顺序一致
head(all_coords)
n_cells <- nrow(all_coords)
barcodes <- rownames(all_coords)

# 取细胞标签（factor）
labels <- as.character(merged_seu@meta.data$cell_type)
names(labels) <- rownames(merged_seu@meta.data)  # 保证名字为 barcode

# 对齐 coords 与 labels（取交集）
common <- intersect(rownames(all_coords), names(labels))
all_coords <- all_coords[common, , drop = FALSE]
labels <- labels[common]
barcodes <- common

xy <- all_coords[,1:2]
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

# z-score 热图（合并两张图的数据）
pheatmap(mat_z,
         cluster_rows = TRUE, 
         cluster_cols = TRUE,
         main = "Z-score of contact enrichment (Combined two images)",
         fontsize = 10)



# 输出统计信息
cat("=== 合并两张图像的分析结果 ===\n")
cat("总细胞数量:", length(common), "\n")
cat("细胞类型:", types, "\n")
cat("边数量:", nrow(edges_df), "\n")
cat("图像1细胞数:", sum(all_coords$sample == image_names[1]), "\n")
cat("图像2细胞数:", sum(all_coords$sample == image_names[2]), "\n")
cat("Z-score 范围:", range(mat_z), "\n")

# 查看观测矩阵
print("观测矩阵:")
print(mat_obs)

# 查看 Z-score 矩阵
print("Z-score 矩阵:")
print(round(mat_z, 2))













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
