# CellCharter R implementation (Parallel version)
# STABLE + INTELLIGENT MULTICORE
library(Seurat)
library(Matrix)
library(mclust)
library(parallel)

aggregate_neighbors <- function(seurat_obj, n_layers = 3, aggregations = "mean", 
                               use_rep = "pca", out_key = "cellcharter", 
                               sample_key = NULL, x_coord_col = "CenterX_global_px", 
                               y_coord_col = "CenterY_global_px", n_cores = 1, 
                               save_intermediates = FALSE, intermediates_dir = NULL) {
  n_cells <- ncol(seurat_obj)
  
  if (!x_coord_col %in% colnames(seurat_obj@meta.data) || !y_coord_col %in% colnames(seurat_obj@meta.data)) {
    stop("Spatial coordinates not found")
  }
  
  spatial_coord <- seurat_obj@meta.data[, c(x_coord_col, y_coord_col)]
  colnames(spatial_coord) <- c("x", "y")
  
  if (!"spatial" %in% names(seurat_obj@graphs)) {
    message("Creating spatial graph...")
    spatial_embeddings <- as.matrix(spatial_coord)
    colnames(spatial_embeddings) <- c("CCSPATIAL_1", "CCSPATIAL_2")
    seurat_obj@reductions[["cellcharter_spatial"]] <- CreateDimReducObject(
      embeddings = spatial_embeddings, key = "CCSPATIAL_"
    )
    seurat_obj <- FindNeighbors(seurat_obj, reduction = "cellcharter_spatial", dims = 1:2, k = 20, graph.name = "spatial")
  }
  
  adj_matrix <- seurat_obj@graphs$spatial
  if (inherits(adj_matrix, "dgCMatrix")) {
    diag(adj_matrix) <- 0
    adj_matrix <- drop0(adj_matrix)
  } else {
    diag(adj_matrix) <- 0
  }
  
  if (!is.null(use_rep)) {
    features <- t(Embeddings(seurat_obj, reduction = use_rep))
  } else {
    features <- as.matrix(GetAssayData(seurat_obj, layer = "data"))
  }
  colnames(features) <- colnames(seurat_obj)
  
  if (is.numeric(n_layers) && length(n_layers) == 1) n_layers <- 0:n_layers
  if (is.character(aggregations)) aggregations <- list(aggregations)
  
  aggregated_list <- list()
  
  for (layer in n_layers) {
    message(paste("Layer", layer))
    if (layer == 0) {
      for (agg in aggregations) {
        aggregated_list[[paste("layer", layer, agg, sep = "_")]] <- features
      }
    } else {
      adj_power <- adj_matrix
      if (layer > 1) for (i in 2:layer) adj_power <- adj_power %*% adj_power
      rs <- rowSums(adj_power)
      rs[rs == 0] <- 1
      adj_norm <- adj_power / rs
      
      for (agg in aggregations) {
        if (agg == "mean") {
          aggregated_list[[paste("layer", layer, agg, sep = "_")]] <- features %*% t(adj_norm)
        } else if (agg == "var") {
          mean_agg <- features %*% t(adj_norm)
          mean_squared_agg <- (features^2) %*% t(adj_norm)
          var_agg <- mean_squared_agg - mean_agg^2
          aggregated_list[[paste("layer", layer, agg, sep = "_")]] <- var_agg
        }
      }
    }
  }
  
  aggregated_features <- do.call(rbind, aggregated_list)
  seurat_obj[[out_key]] <- CreateAssayObject(counts = as.matrix(aggregated_features))
  DefaultAssay(seurat_obj) <- out_key
  
  if (save_intermediates) {
    saveRDS(list(seurat_obj = seurat_obj, combined_features = aggregated_features),
            file.path(intermediates_dir, "seurat_after_combining.rds"))
  }
  return(seurat_obj)
}

# --------------------------
# GMM 聚类 → 自动限制 6 核
# --------------------------
cluster_cells <- function(seurat_obj, n_clusters = 5, cluster_key = "cellcharter_cluster", n_cores = 6) {
  message("Clustering cells (GMM) - AUTO LIMIT TO 6 CORES")
  features <- t(as.matrix(GetAssayData(seurat_obj, layer = "data")))
  
  if (ncol(features) > 50) {
    message("Reducing to 50 dims for stability")
    features <- prcomp(features, rank. = 50)$x
  }
  
  safe_cores <- min(n_cores, 6)
  message(paste("Using", safe_cores, "cores for GMM"))
  gmm <- Mclust(features, G = n_clusters, modelNames = "EEE", ncores = safe_cores)
  seurat_obj@meta.data[[cluster_key]] <- as.factor(gmm$classification)
  return(seurat_obj)
}

# --------------------------
# 边界细胞 → 必须单线程（不丢细胞）
# --------------------------
identify_boundary_cells <- function(seurat_obj, cluster_col = "cellcharter_cluster", n_cores = 1) {
  message("Identifying boundary cells (SINGLE THREAD - NO BUG)")
  clu <- seurat_obj@meta.data[[cluster_col]]
  adj <- seurat_obj@graphs$spatial
  n <- length(clu)
  res <- logical(n)
  
  for (i in 1:n) {
    nbs <- which(adj[i,] > 0)
    res[i] <- length(nbs) > 0 && any(clu[nbs] != clu[i])
  }
  
  seurat_obj@meta.data$boundary_cell <- res
  return(seurat_obj)
}

# --------------------------
# 形状计算 → 多核拉满
# --------------------------
calculate_shape_metrics <- function(seurat_obj, cluster_col = "cellcharter_cluster", n_cores = 1) {
  message("Calculating shape metrics (MULTICORE)")
  clu <- seurat_obj@meta.data[[cluster_col]]
  coords <- seurat_obj@meta.data[, c("CenterX_global_px", "CenterY_global_px")]
  colnames(coords) <- c("x", "y")
  uq <- unique(clu)
  
  fun <- function(cid) {
    idx <- which(clu == cid)
    co <- coords[idx,]
    ar <- 0
    if (nrow(co) > 2) {
      hu <- chull(co$x, co$y)
      hc <- co[hu,]
      ar <- 0.5 * abs(sum(hc$x[-1] * hc$y[-nrow(hc)] - hc$x[-nrow(hc)] * hc$y[-1]))
    }
    cx <- mean(co$x)
    cy <- mean(co$y)
    el <- 1
    if (nrow(co) > 1) {
      cv <- cov(co)
      eg <- eigen(cv, only.values = TRUE)$values
      if (min(eg) > 0) el <- sqrt(max(eg) / min(eg))
    }
    data.frame(cluster = cid, area = ar, centroid_x = cx, centroid_y = cy, elongation = el, n_cells = length(idx))
  }
  
  if (n_cores > 1) {
    out <- mclapply(uq, fun, mc.cores = n_cores)
  } else {
    out <- lapply(uq, fun)
  }
  do.call(rbind, out)
}

# --------------------------
# 空间域分析
# --------------------------
analyze_spatial_domains <- function(seurat_obj, out_dir = "./cellcharter_results_parallel", 
                                   cell_type_col = "detailed", n_cores = 1) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  seurat_obj <- identify_boundary_cells(seurat_obj, n_cores = 1)
  shapes <- calculate_shape_metrics(seurat_obj, n_cores = n_cores)
  write.csv(shapes, file.path(out_dir, "shape_metrics.csv"), row.names = FALSE)
  ct <- table(seurat_obj@meta.data[[cell_type_col]], seurat_obj@meta.data$cellcharter_cluster)
  write.csv(as.data.frame.matrix(ct), file.path(out_dir, "domain_cell_type_composition.csv"))
  return(seurat_obj)
}

# --------------------------
# 聚类稳定性分析
# --------------------------
cluster_autok <- function(seurat_obj, n_clusters_range = 2:10, max_runs = 10, convergence_tol = 0.01, 
                         use_rep = "cellcharter", n_cores = 6, out_dir = "./cellcharter_autok",
                         sample_size = 0.1) {  # 添加抽样比例参数，默认10%
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  if (!use_rep %in% names(seurat_obj@assays)) {
    stop(paste(use_rep, "not found in seurat_obj@assays"))
  }
  
  # 获取所有细胞的特征
  all_features <- t(as.matrix(GetAssayData(seurat_obj, assay = use_rep, layer = "data")))
  
  # 抽样
  n_cells <- nrow(all_features)
  if (sample_size < 1) {
    sample_size <- max(1, round(n_cells * sample_size))
  }
  sample_indices <- sample(1:n_cells, size = sample_size, replace = FALSE)
  features <- all_features[sample_indices, ]
  
  message(paste("Sampling", length(sample_indices), "cells for stability analysis"))
  
  if (ncol(features) > 50) {
    message("Reducing to 50 dims for stability")
    pca <- prcomp(features, rank. = 50)
    features <- pca$x
  }
  
  labels <- list()
  best_models <- list()
  stability <- list()
  
  previous_stability <- NULL
  
  for (i in 1:max_runs) {
    message(paste("Iteration", i, "/", max_runs))
    new_labels <- list()
    
    for (k in n_clusters_range) {
      message(paste("Clustering with k =", k))
      gmm <- Mclust(features, G = k, modelNames = "EEE", ncores = min(n_cores, 6))
      new_labels[[as.character(k)]] <- gmm$classification
      
      if (!as.character(k) %in% names(best_models) || gmm$loglik > best_models[[as.character(k)]]$loglik) {
        best_models[[as.character(k)]] <- gmm
      }
    }
    
    if (i > 1) {
      current_stability <- list()
      for (k_idx in 1:(length(n_clusters_range)-1)) {
        k_val <- n_clusters_range[k_idx]
        k_plus_1_val <- n_clusters_range[k_idx+1]
        
        for (j in 1:(i-1)) {
          # 计算 k 和 k+1 之间的稳定性
          fm_score <- fowlkes_mallows_score(new_labels[[as.character(k_val)]], labels[[as.character(k_plus_1_val)]][[j]])
          current_stability[[length(current_stability) + 1]] <- fm_score
        }
      }
      
      stability[[i-1]] <- current_stability
      
      if (!is.null(previous_stability)) {
        # 计算稳定性变化
        stability_change <- mean(abs(unlist(stability[[i-1]]) - unlist(previous_stability)) / abs(unlist(previous_stability) + 1e-9))
        if (stability_change < convergence_tol) {
          message(paste("Convergence with a change in stability of", stability_change, "reached after", i, "iterations"))
          
          # 保存新标签
          for (k in names(new_labels)) {
            if (!k %in% names(labels)) {
              labels[[k]] <- list()
            }
            labels[[k]][[i]] <- new_labels[[k]]
          }
          break
        }
      }
      
      previous_stability <- current_stability
    }
    
    # 保存新标签
    for (k in names(new_labels)) {
      if (!k %in% names(labels)) {
        labels[[k]] <- list()
      }
      labels[[k]][[i]] <- new_labels[[k]]
    }
  }
  
  # 计算稳定性矩阵
  if (length(stability) > 0) {
    # 计算实际的稳定性矩阵维度
    n_rows <- length(n_clusters_range) - 1
    n_cols <- length(stability)
    
    # 创建稳定性矩阵
    stability_matrix <- matrix(0, nrow = n_rows, ncol = n_cols)
    
    # 填充稳定性矩阵
    for (i in 1:n_cols) {
      # 确保不超出范围
      n_vals <- length(stability[[i]])
      if (n_vals > 0) {
        # 只填充有效的值
        for (j in 1:min(n_rows, n_vals)) {
          stability_matrix[j, i] <- stability[[i]][[j]]
        }
      }
    }
    
    # 计算平均稳定性
    stability_mean <- rowMeans(stability_matrix)
    
    # 寻找稳定性曲线的局部最大值
    peaks <- find_peaks(stability_mean)
    best_k_candidates <- n_clusters_range[peaks + 1]
    
    # 选择最佳 k 值
    best_k <- n_clusters_range[which.max(stability_mean) + 1]
    
    # 保存结果
    saveRDS(list(
      labels = labels,
      best_models = best_models,
      stability_matrix = stability_matrix,
      stability_mean = stability_mean,
      best_k = best_k,
      best_k_candidates = best_k_candidates,
      sample_indices = sample_indices
    ), file.path(out_dir, "autok_results.rds"))
    
    # 可视化稳定性曲线
    pdf(file.path(out_dir, "stability_curve.pdf"))
    plot(n_clusters_range[-1], stability_mean, type = "b", xlab = "Number of clusters (k)", 
         ylab = "Stability score", main = "Cluster stability across k values")
    points(n_clusters_range[peaks + 1], stability_mean[peaks], col = "red", pch = 16, cex = 1.5)
    abline(v = best_k, col = "blue", lty = 2)
    legend("topright", c("Stability", "Local maxima", "Best k"), 
           col = c("black", "red", "blue"), pch = c(1, 16, NA), lty = c(1, NA, 2))
    dev.off()
    
    message(paste("Best k found:", best_k))
    message(paste("Local maxima k candidates:", paste(best_k_candidates, collapse = ", ")))
    
    return(list(
      best_k = best_k,
      best_k_candidates = best_k_candidates,
      stability_mean = stability_mean,
      best_models = best_models
    ))
  } else {
    stop("Could not compute stability")
  }
}

# Fowlkes-Mallows 评分函数
fowlkes_mallows_score <- function(labels1, labels2) {
  contingency <- table(labels1, labels2)
  tp <- sum(contingency^2) - sum(rowSums(contingency)^2) - sum(colSums(contingency)^2) + sum(contingency)^2
  tp <- tp / 2
  fp <- (sum(rowSums(contingency)^2) - sum(contingency)) / 2
  fn <- (sum(colSums(contingency)^2) - sum(contingency)) / 2
  precision <- tp / (tp + fp)
  recall <- tp / (tp + fn)
  return(sqrt(precision * recall))
}

# 寻找局部最大值
find_peaks <- function(x) {
  peaks <- c()
  for (i in 2:(length(x)-1)) {
    if (x[i] > x[i-1] && x[i] > x[i+1]) {
      peaks <- c(peaks, i)
    }
  }
  return(peaks)
}

# --------------------------
# 主流程
# --------------------------
run_cellcharter <- function(seurat_obj, n_layers = 3, n_clusters = 5, use_rep = "pca", cell_type_col = "detailed",
                           out_dir = "./cellcharter_results_parallel", x_coord_col = "CenterX_global_px", 
                           y_coord_col = "CenterY_global_px", n_cores = 20, save_intermediates = FALSE) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  idir <- file.path(out_dir, "intermediates")
  if (save_intermediates && !dir.exists(idir)) dir.create(idir, recursive = TRUE)
  
  seurat_obj <- aggregate_neighbors(seurat_obj, n_layers = n_layers, use_rep = use_rep,
                                   x_coord_col = x_coord_col, y_coord_col = y_coord_col,
                                   n_cores = n_cores, save_intermediates = save_intermediates,
                                   intermediates_dir = idir)
  
  seurat_obj <- cluster_cells(seurat_obj, n_clusters = n_clusters, n_cores = n_cores)
  seurat_obj <- analyze_spatial_domains(seurat_obj, out_dir = out_dir, cell_type_col = cell_type_col, n_cores = n_cores)
  return(seurat_obj)
}

run_cellcharter_autok <- function(n_layers = 3, n_clusters_range = 2:10, n_cores = 20, save_intermediates = TRUE, sample_size = 0.1) {
  if (!exists("mPT")) stop("mPT not found")
  mPT <- get("mPT", envir = .GlobalEnv)
  
  # 运行 CellCharter 聚合邻居
  out_dir <- "./cellcharter_results_autok"
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  idir <- file.path(out_dir, "intermediates")
  if (save_intermediates && !dir.exists(idir)) dir.create(idir, recursive = TRUE)
  
  mPT <- aggregate_neighbors(mPT, n_layers = n_layers, use_rep = "pca",
                            x_coord_col = "CenterX_global_px", y_coord_col = "CenterY_global_px",
                            n_cores = n_cores, save_intermediates = save_intermediates,
                            intermediates_dir = idir)
  
  # 运行自动选择 k 值（使用抽样）
  autok_results <- cluster_autok(mPT, n_clusters_range = n_clusters_range, max_runs = 10, 
                               n_cores = n_cores, out_dir = file.path(out_dir, "autok"),
                               sample_size = sample_size)
  
  # 使用最佳 k 值进行聚类（使用全部细胞）
  best_k <- autok_results$best_k
  message(paste("Using best k:", best_k))
  
  mPT <- cluster_cells(mPT, n_clusters = best_k, n_cores = n_cores)
  mPT <- analyze_spatial_domains(mPT, out_dir = out_dir, cell_type_col = "detailed", n_cores = n_cores)
  
  # 保存结果
  saveRDS(mPT, file.path(out_dir, "mPT_cellcharter_AUTO_K_FINAL.rds"))
  assign("mPT", mPT, envir = .GlobalEnv)
  return(mPT)
}

run_cellcharter_mPT <- function(n_layers = 3, n_clusters = 5, n_cores = 20, save_intermediates = TRUE) {
  if (!exists("mPT")) stop("mPT not found")
  mPT <- get("mPT", envir = .GlobalEnv)
  mPT <- run_cellcharter(mPT, n_layers = n_layers, n_clusters = n_clusters, n_cores = n_cores,
                        use_rep = "pca", cell_type_col = "detailed",
                        x_coord_col = "CenterX_global_px", y_coord_col = "CenterY_global_px",
                        save_intermediates = save_intermediates)
  assign("mPT", mPT, envir = .GlobalEnv)
  return(mPT)
}