library(Seurat)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(deldir)
library(doParallel)
library(foreach)

# ===================== 数据加载 =====================
obj <- readRDS("YA2025263-1_fin.rds")

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
                       "Basal cells" = "Malignant cells")

# ===================== 提取 nmPT 样本 =====================
nmPT <- subset(obj, subset = tissue == "nmPT")
samples <- unique(nmPT$sample)
cat("nmPT 样本:", paste(samples, collapse = ", "), "\n")

# ===================== 定义空间分析函数 =====================
run_spatial_analysis <- function(coords_sub, labels_sub, sample_name, nperm = 1000) {
  
  cat("\n========================================\n")
  cat("开始分析样本：", sample_name, "\n")
  cat("细胞数量：", length(labels_sub), "\n")
  cat("========================================\n")
  
  # 检查细胞类型数量
  types <- sort(unique(labels_sub))
  K <- length(types)
  
  if (K < 2) {
    cat("警告：", sample_name, "只有", K, "种细胞类型，跳过分析\n")
    return(NULL)
  }
  
  cat("细胞类型数量：", K, "\n")
  cat("细胞类型：", paste(types, collapse = ", "), "\n")
  
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
  
  # 构建观察矩阵
  type_by_index <- labels_sub
  t1 <- type_by_index[edges_df$from]
  t2 <- type_by_index[edges_df$to]
  
  mat_obs <- matrix(0, nrow = K, ncol = K, dimnames = list(types, types))
  for(i in seq_along(t1)) {
    a <- t1[i]
    b <- t2[i]
    mat_obs[a,b] <- mat_obs[a,b] + 1
    mat_obs[b,a] <- mat_obs[b,a] + 1
  }
  
  # 排列检验
  cat("\n开始排列检验（", nperm, "次）...\n")
  
  from_idx <- edges_df$from
  to_idx <- edges_df$to
  n_cells <- length(type_by_index)
  
  ncores <- parallel::detectCores() - 1
  ncores <- max(1, min(ncores, 24))
  
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
  
  # z-score
  z_vec <- (obs_vec - mu_rand) / (sd_rand + 1e-8)
  
  # 经验p值
  p_emp <- sapply(seq_along(obs_vec), function(i) {
    perm_i <- perm_counts[, i]
    obs_i <- obs_vec[i]
    mu_i <- mu_rand[i]
    p_val <- (sum(abs(perm_i - mu_i) >= abs(obs_i - mu_i)) + 1) / (nperm + 1)
    return(p_val)
  })
  
  # 转为矩阵
  mat_mu <- matrix(mu_rand, nrow = K, ncol = K, dimnames = list(types, types))
  mat_sd <- matrix(sd_rand, nrow = K, ncol = K, dimnames = list(types, types))
  mat_z <- matrix(z_vec, nrow = K, ncol = K, dimnames = list(types, types))
  mat_p <- matrix(p_emp, nrow = K, ncol = K, dimnames = list(types, types))
  
  # 计算 log2FC
  mat_fc <- mat_obs / (mat_mu + 1e-8)
  mat_log2fc <- log2(mat_fc)
  
  # 返回结果
  results <- list(
    sample = sample_name,
    mat_obs = mat_obs,
    mat_z = mat_z,
    mat_p = mat_p,
    mat_mu = mat_mu,
    mat_sd = mat_sd,
    mat_log2fc = mat_log2fc,
    cell_types = types,
    n_cells = length(labels_sub),
    n_edges = nrow(edges_df)
  )
  
  return(results)
}

# ===================== 按样本分别分析 =====================
all_results <- list()

for (sample_id in samples) {
  cat("\n\n########## 处理样本：", sample_id, " ##########\n")
  
  # 提取该样本的细胞
  sample_cells <- subset(nmPT, subset = sample == sample_id)
  
  if (nrow(sample_cells@meta.data) == 0) {
    cat("警告：", sample_id, "没有细胞，跳过\n")
    next
  }
  
  # 获取坐标和标签
  coords <- sample_cells@meta.data[, c("CenterX_global_px", "CenterY_global_px")]
  coords <- coords[!is.na(coords[,1]) & !is.na(coords[,2]), , drop = FALSE]
  labels <- as.character(sample_cells@meta.data$detailed)
  names(labels) <- rownames(coords)
  
  # 移除NA标签
  valid_idx <- !is.na(labels)
  coords <- coords[valid_idx, , drop = FALSE]
  labels <- labels[valid_idx]
  
  cat("有效细胞数量：", length(labels), "\n")
  cat("细胞类型分布：\n")
  print(table(labels))
  
  # 过滤稀有细胞类型（可选，数量 < 50 的过滤）
  cell_counts <- table(labels)
  rare_threshold <- 50
  keep_types <- names(cell_counts[cell_counts >= rare_threshold])
  keep_cells <- labels %in% keep_types
  coords <- coords[keep_cells, , drop = FALSE]
  labels <- labels[keep_cells]
  
  cat("过滤稀有细胞后数量：", length(labels), "\n")
  
  if (length(labels) < 100) {
    cat("警告：", sample_id, "过滤后细胞太少，跳过\n")
    next
  }
  
  # 确保坐标是数值型
  coords <- as.matrix(coords)
  colnames(coords) <- c("x", "y")
  
  # 运行分析
  result <- run_spatial_analysis(coords, labels, sample_id, nperm = 1000)
  
  if (!is.null(result)) {
    all_results[[sample_id]] <- result
    
    # 绘制并保存热图
    if (nrow(result$mat_log2fc) > 1 && ncol(result$mat_log2fc) > 1) {
      # ===== 修改：处理 -Inf 和 +Inf，然后裁剪 =====
      mat_plot <- result$mat_log2fc
      
      # 替换 -Inf 为 -10
      mat_plot[is.infinite(mat_plot) & mat_plot < 0] <- -10
      # 替换 +Inf 为 10
      mat_plot[is.infinite(mat_plot) & mat_plot > 0] <- 10
      
      # 裁剪超出 [-10, 10] 的范围
      mat_plot[mat_plot > 10] <- 10
      mat_plot[mat_plot < -10] <- -10
      
      # 显著性标注
      signif_symbols <- matrix("", nrow = length(result$cell_types), 
                                ncol = length(result$cell_types),
                                dimnames = list(result$cell_types, result$cell_types))
      signif_symbols[result$mat_p < 0.001] <- "***"
      signif_symbols[result$mat_p < 0.01 & result$mat_p >= 0.001] <- "**"
      signif_symbols[result$mat_p < 0.05 & result$mat_p >= 0.01] <- "*"
      
      # 热图
      p <- pheatmap(mat_plot,
                    cluster_rows = TRUE,
                    cluster_cols = TRUE,
                    main = paste("Log2(Fold Change) -", sample_id),
                    fontsize = 10,
                    display_numbers = signif_symbols,
                    number_color = "black",
                    fontsize_number = 8,
                    #filename = paste0("nmPT_", sample_id, "_contact_log2fc.pdf"),
                    width = 12, height = 10)
      
      cat("保存热图：nmPT_", sample_id, "_contact_log2fc.pdf\n", sep = "")
    }
    
    # 保存结果文件
    #saveRDS(result, file = paste0("nmPT_", sample_id, "_results.rds"))
    cat("保存结果：nmPT_", sample_id, "_results.rds\n", sep = "")
  }
}
saveRDS(all_results, file = "nmPT_all_results_subtype_interaction.rds")
# ===================== 整合多样本结果 =====================
all_cell_types <- unique(unlist(lapply(all_results, function(x) x$cell_types)))
cat("总细胞类型数量：", length(all_cell_types), "\n")
cat("所有细胞类型：", paste(all_cell_types, collapse = ", "), "\n")

if (length(all_cell_types) >= 2) {
  # 计算平均 log2FC（初始化为0）
  avg_log2fc <- matrix(0, nrow = length(all_cell_types), ncol = length(all_cell_types),
                       dimnames = list(all_cell_types, all_cell_types))
  
  # 记录每个单元格有多少个有效样本
  n_effective <- matrix(0, nrow = length(all_cell_types), ncol = length(all_cell_types),
                         dimnames = list(all_cell_types, all_cell_types))
  
  # 存储合并 chi2 和自由度（也要设置 dimnames）
  combined_chi2 <- matrix(0, nrow = length(all_cell_types), ncol = length(all_cell_types),
                          dimnames = list(all_cell_types, all_cell_types))
  combined_df <- matrix(0, nrow = length(all_cell_types), ncol = length(all_cell_types),
                        dimnames = list(all_cell_types, all_cell_types))
  
  for (sample_id in names(all_results)) {
    res <- all_results[[sample_id]]
    
    # 当前样本的细胞类型
    current_types <- res$cell_types
    
    # 只处理当前样本中存在的细胞类型
    for (i in seq_along(current_types)) {
      for (j in seq_along(current_types)) {
        ct_i <- current_types[i]
        ct_j <- current_types[j]
        
        avg_log2fc[ct_i, ct_j] <- avg_log2fc[ct_i, ct_j] + res$mat_log2fc[ct_i, ct_j]
        n_effective[ct_i, ct_j] <- n_effective[ct_i, ct_j] + 1
        
        combined_chi2[ct_i, ct_j] <- combined_chi2[ct_i, ct_j] + (-2 * log(res$mat_p[ct_i, ct_j] + 1e-10))
        combined_df[ct_i, ct_j] <- combined_df[ct_i, ct_j] + 2
      }
    }
  }
  
  # 计算平均 log2FC
  avg_log2fc <- avg_log2fc / n_effective
  
  # 计算合并 p 值
  combined_p <- matrix(1, nrow = length(all_cell_types), ncol = length(all_cell_types),
                       dimnames = list(all_cell_types, all_cell_types))
  for (i in seq_along(all_cell_types)) {
    for (j in seq_along(all_cell_types)) {
      if (combined_df[i, j] > 0) {
        combined_p[i, j] <- pchisq(combined_chi2[i, j], df = combined_df[i, j], lower.tail = FALSE)
      }
    }
  }
  
  # 处理 NaN 和 Inf
  avg_log2fc[is.nan(avg_log2fc)] <- 0
  avg_log2fc[is.infinite(avg_log2fc) & avg_log2fc < 0] <- -10
  avg_log2fc[is.infinite(avg_log2fc) & avg_log2fc > 0] <- 10
  avg_log2fc[avg_log2fc > 10] <- 10
  avg_log2fc[avg_log2fc < -10] <- -10
  
  # 显著性标注（只标星号）
  signif_symbols <- matrix("", nrow = length(all_cell_types), ncol = length(all_cell_types),
                           dimnames = list(all_cell_types, all_cell_types))
  signif_symbols[combined_p < 0.001] <- "***"
  signif_symbols[combined_p < 0.01 & combined_p >= 0.001] <- "**"
  signif_symbols[combined_p < 0.05 & combined_p >= 0.01] <- "*"
  
  # 绘制平均热图
  p_avg <- pheatmap(avg_log2fc,
                    cluster_rows = TRUE,
                    cluster_cols = TRUE,
                    main = "Average Log2(Fold Change) across nmPT samples",
                    fontsize = 10,
                    display_numbers = signif_symbols,
                    number_color = "black",
                    fontsize_number = 8,
                    filename = "nmPT_all_samples_average_contact_log2fc_detailed.pdf",
                    width = 14, height = 12)
  
  cat("保存平均热图：nmPT_all_samples_average_contact_log2fc.pdf\n")
  
  # 保存整合结果
  integrated_results <- list(
    all_cell_types = all_cell_types,
    avg_log2fc = avg_log2fc,
    combined_p = combined_p,
    n_effective = n_effective,
    n_samples = length(all_results),
    individual_results = all_results
  )
  
  saveRDS(integrated_results, file = "nmPT_integrated_results_detailed.rds")
  cat("保存整合结果：nmPT_integrated_results.rds\n")
}

# ===================== 打印汇总信息 =====================
cat("\n\n========== 分析汇总 ==========\n")
for (sample_id in names(all_results)) {
  res <- all_results[[sample_id]]
  cat(sample_id, ": ", res$n_cells, "个细胞, ", 
      length(res$cell_types), "种细胞类型, ", 
      res$n_edges, "条边\n", sep = "")
}

cat("\n===== 所有分析完成！=====\n")