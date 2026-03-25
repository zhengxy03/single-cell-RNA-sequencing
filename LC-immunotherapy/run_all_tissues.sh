#!/bin/bash
# 文件：run_all_tissues.sh
# 功能：为每个组织生成独立的R脚本并提交作业

# 定义组织列表（根据你的数据）
TISSUES=("mPT" "nmPT" "metLN" "negLN")

# 遍历每个组织
for TISSUE in "${TISSUES[@]}"; do
    echo "准备提交作业：$TISSUE"
    
    # 为每个组织生成独立的R脚本
    cat > spatial_analysis_${TISSUE}.R << 'EOF'
# ===================== 加载必要的包 =====================
library(Seurat)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(deldir)
library(doParallel)
library(foreach)

# ===================== 获取命令行参数 =====================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    stop("请指定组织名称！")
}
TISSUE_NAME <- args[1]
cat("分析组织：", TISSUE_NAME, "\n")

# ===================== 数据加载 =====================
obj <- readRDS("YA2025263-1_fin.rds")
Malignant <- readRDS("malignant_anno.rds")

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
                       "Malignant cells" = "unknown",
                       "Basal cells" = "Malignant cells")

# 初始化detailed列
obj@meta.data$detailed <- NA

# 优先级1：填充sub_cell_type
common_cells <- intersect(rownames(Malignant@meta.data), rownames(obj@meta.data))
obj@meta.data[common_cells, "sub_cell_type"] <- Malignant@meta.data[common_cells, "sub_cell_type"]
sub_cell_idx <- !is.na(obj@meta.data$sub_cell_type)
obj@meta.data$detailed[sub_cell_idx] <- as.character(obj@meta.data$sub_cell_type[sub_cell_idx])

# 优先级2：填充CellType
cell_type_levels <- levels(obj@meta.data$CellType)
cell_type_idx <- is.na(obj@meta.data$detailed) & !is.na(obj@meta.data$CellType)
if (is.factor(obj@meta.data$CellType)) {
  obj@meta.data$detailed[cell_type_idx] <- cell_type_levels[as.integer(obj@meta.data$CellType[cell_type_idx])]
} else {
  obj@meta.data$detailed[cell_type_idx] <- as.character(obj@meta.data$CellType[cell_type_idx])
}

# ===================== 定义空间邻域分析函数 =====================
run_spatial_analysis <- function(coords_sub, labels_sub, tissue_name, nperm = 1000) {
  
  cat("\n========================================\n")
  cat("开始分析：", tissue_name, "\n")
  cat("细胞数量：", length(labels_sub), "\n")
  cat("========================================\n")
  
  types <- sort(unique(labels_sub))
  K <- length(types)
  
  if (K < 2) {
    cat("警告：只有", K, "种细胞类型，跳过分析\n")
    return(NULL)
  }
  
  cat("细胞类型数量：", K, "\n")
  
  # 提取坐标
  x <- coords_sub[,1]
  y <- coords_sub[,2]
  
  # 计算Delaunay三角剖分
  cat("计算Delaunay三角剖分...\n")
  deld <- deldir(x, y, rw = c(range(x), range(y)))
  segs <- deld$delsgs
  cat("生成三角边数量：", nrow(segs), "\n")
  
  # 构建边
  edges <- cbind(segs$ind1, segs$ind2)
  edges <- edges[edges[,1] != edges[,2], , drop=FALSE]
  edges <- t(apply(edges, 1, function(x) sort(x)))
  edges <- unique(edges)
  edges_df <- data.frame(from = edges[,1], to = edges[,2])
  cat("最终边数量：", nrow(edges_df), "\n")
  
  type_by_index <- labels_sub
  
  # 构建观察矩阵
  t1 <- type_by_index[edges_df$from]
  t2 <- type_by_index[edges_df$to]
  
  mat_obs <- matrix(0, nrow = K, ncol = K, dimnames = list(types, types))
  for(i in seq_along(t1)) {
    a <- t1[i]
    b <- t2[i]
    mat_obs[a,b] <- mat_obs[a,b] + 1
    mat_obs[b,a] <- mat_obs[b,a] + 1
  }
  
  # ===== 排列检验 =====
  cat("\n开始排列检验（", nperm, "次）...\n")
  
  from_idx <- edges_df$from
  to_idx <- edges_df$to
  n_cells <- length(type_by_index)
  
  ncores <- 24
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  
  perm_counts <- foreach(p = 1:nperm, .packages = c(), .combine = rbind) %dopar% {
    set.seed(p + 12345)
    perm_labels <- sample(type_by_index, n_cells, replace = FALSE)
    
    pt1 <- perm_labels[from_idx]
    pt2 <- perm_labels[to_idx]
    
    mat <- matrix(0, nrow = K, ncol = K, dimnames = list(types, types))
    for(i in seq_along(pt1)) {
      a <- pt1[i]
      b <- pt2[i]
      mat[a,b] <- mat[a,b] + 1
      mat[b,a] <- mat[b,a] + 1
    }
    as.vector(mat)
  }
  
  stopCluster(cl)
  cat("排列检验完成\n")
  
  # 计算结果
  obs_vec <- as.vector(mat_obs)
  mu_rand <- colMeans(perm_counts)
  sd_rand <- apply(perm_counts, 2, sd)
  
  z_vec <- (obs_vec - mu_rand) / (sd_rand + 1e-8)
  
  p_emp <- sapply(seq_along(obs_vec), function(i) {
    perm_i <- perm_counts[, i]
    obs_i <- obs_vec[i]
    mu_i <- mu_rand[i]
    p_val <- (sum(abs(perm_i - mu_i) >= abs(obs_i - mu_i)) + 1) / (nperm + 1)
    return(p_val)
  })
  
  mat_mu <- matrix(mu_rand, nrow = K, ncol = K, dimnames = list(types, types))
  mat_sd <- matrix(sd_rand, nrow = K, ncol = K, dimnames = list(types, types))
  mat_z <- matrix(z_vec, nrow = K, ncol = K, dimnames = list(types, types))
  mat_p <- matrix(p_emp, nrow = K, ncol = K, dimnames = list(types, types))
  
  results <- list(
    tissue = tissue_name,
    mat_obs = mat_obs,
    mat_z = mat_z,
    mat_p = mat_p,
    mat_mu = mat_mu,
    mat_sd = mat_sd,
    cell_types = types,
    n_cells = length(labels_sub),
    n_edges = nrow(edges_df)
  )
  
  return(results)
}

# ===================== 分析指定的组织 =====================
cat("\n########## 处理组织：", TISSUE_NAME, " ##########\n")

# 筛选该组织的细胞
tissue_cells <- obj@meta.data[obj@meta.data$tissue == TISSUE_NAME, ]

if (nrow(tissue_cells) == 0) {
  stop("组织 ", TISSUE_NAME, " 没有细胞！")
}

# 获取坐标和标签
coords_tissue <- tissue_cells[, c("CenterX_global_px", "CenterY_global_px")]
coords_tissue <- coords_tissue[!is.na(coords_tissue[,1]) & !is.na(coords_tissue[,2]), , drop = FALSE]
labels_tissue <- tissue_cells$detailed[rownames(coords_tissue)]

# 移除NA标签
valid_idx <- !is.na(labels_tissue)
coords_tissue <- coords_tissue[valid_idx, , drop = FALSE]
labels_tissue <- labels_tissue[valid_idx]

cat("有效细胞数量：", length(labels_tissue), "\n")
cat("细胞类型分布：\n")
print(table(labels_tissue))

# 确保坐标是数值型
coords_tissue <- as.matrix(coords_tissue)

# 运行分析
result <- run_spatial_analysis(coords_tissue, labels_tissue, TISSUE_NAME, nperm = 1000)

if (!is.null(result)) {
  # 绘制并保存热图
  if (nrow(result$mat_z) > 1 && ncol(result$mat_z) > 1) {
    pdf(paste0("spatial_contact_", TISSUE_NAME, ".pdf"), width = 12, height = 10)
    pheatmap(result$mat_z,
             cluster_rows = TRUE,
             cluster_cols = TRUE,
             main = paste("Z-score of contact enrichment -", TISSUE_NAME),
             fontsize = 8)
    dev.off()
    cat("保存热图：", paste0("spatial_contact_", TISSUE_NAME, ".pdf"), "\n")
  }
  
  # 保存结果文件
  saveRDS(result, file = paste0("spatial_results_", TISSUE_NAME, ".rds"))
  cat("保存结果：", paste0("spatial_results_", TISSUE_NAME, ".rds"), "\n")
}

cat("\n===== 分析完成！=====\n")
EOF

    # 为每个组织创建独立的作业提交脚本
    cat > submit_${TISSUE}.sh << EOF
#!/bin/bash
#BSUB -q mpi
#BSUB -n 24
#BSUB -J spatial_${TISSUE}
#BSUB -o spatial_${TISSUE}_%J.out
#BSUB -e spatial_${TISSUE}_%J.err

# 加载R模块
module load R/4.2.0

# 运行R脚本，传入组织名称作为参数
Rscript spatial_analysis_${TISSUE}.R ${TISSUE}

echo "Spatial analysis for ${TISSUE} completed at \$(date)"
EOF

    # 提交作业
    bsub < submit_${TISSUE}.sh
    
    echo "已提交作业：spatial_${TISSUE}"
    echo "----------------------------------------"
done

echo "所有作业已提交完成！"