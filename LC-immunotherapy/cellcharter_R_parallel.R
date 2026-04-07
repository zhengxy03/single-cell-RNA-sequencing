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