# CellCharter R implementation
# Based on the Python package cellcharter-0.3.7

library(Seurat)
library(Matrix)
library(mclust)

#' Aggregate spatial neighbor features
#' 
#' Similar to cellcharter's aggregate_neighbors function, this function aggregates
#' features from spatial neighbors at different layers and concatenates them with
#' the original cell features.
#' 
#' @param seurat_obj Seurat object with spatial coordinates
#' @param n_layers Integer or list specifying which neighborhood layers to aggregate
#' @param aggregations Character or list of aggregation methods (default: "mean")
#' @param use_rep Name of the dimensionality reduction to use (default: "pca")
#' @param out_key Name of the output key to store the aggregated features
#' @param sample_key Name of the sample identifier in meta.data (default: NULL)
#' @param x_coord_col Name of the metadata column containing x-coordinates (default: "CenterX_global_px")
#' @param y_coord_col Name of the metadata column containing y-coordinates (default: "CenterY_global_px")
#' 
#' @return Seurat object with aggregated features stored in assay
#' 
aggregate_neighbors <- function(seurat_obj, n_layers = 3, aggregations = "mean", 
                               use_rep = "pca", out_key = "cellcharter", 
                               sample_key = NULL, x_coord_col = "CenterX_global_px", 
                               y_coord_col = "CenterY_global_px") {
  # Check if spatial coordinates are present in metadata
  if (!x_coord_col %in% colnames(seurat_obj@meta.data) || !y_coord_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Spatial coordinates not found in metadata. Please check", x_coord_col, "and", y_coord_col))
  }
  
  # Get spatial coordinates
  spatial_coord <- seurat_obj@meta.data[, c(x_coord_col, y_coord_col)]
  colnames(spatial_coord) <- c("x", "y")
  
  # Create spatial neighbor graph if not present
  if (!"spatial" %in% names(seurat_obj@graphs)) {
    message("Creating spatial neighbor graph...")
    # Use a unique name for the spatial coordinates reduction to avoid conflicts
    spatial_reduction_name <- "cellcharter_spatial"
    
    # Add spatial coordinates to the Seurat object for FindNeighbors
    spatial_embeddings <- as.matrix(spatial_coord)
    colnames(spatial_embeddings) <- c("CCSPATIAL_1", "CCSPATIAL_2")
    seurat_obj@reductions[[spatial_reduction_name]] <- CreateDimReducObject(
      embeddings = spatial_embeddings,
      key = "CCSPATIAL_"
    )
    
    seurat_obj <- FindNeighbors(seurat_obj, 
                               reduction = spatial_reduction_name,
                               dims = 1:2,
                               k = 20,
                               graph.name = "spatial")
  }
  
  # Get adjacency matrix
  adj_matrix <- seurat_obj@graphs$spatial
  
  # Set diagonal to 0 (remove self-loops) without converting to dense
  if (is(adj_matrix, "dgCMatrix")) {
    # For sparse matrix, set diagonal to 0
    # Create a copy to avoid modifying the original graph
    adj_matrix <- as(adj_matrix, "dgCMatrix")
    # Get the diagonal elements and set them to 0
    for (i in 1:nrow(adj_matrix)) {
      # Find if there's a diagonal element
      idx <- which(adj_matrix@i + 1 == i & adj_matrix@p[i] < adj_matrix@p[i+1])
      if (length(idx) > 0) {
        pos <- adj_matrix@p[i] + idx
        adj_matrix@x[pos] <- 0
      }
    }
    adj_matrix <- drop0(adj_matrix)
  } else {
    # For dense matrix, set diagonal to 0
    diag(adj_matrix) <- 0
  }
  
  # Get features to aggregate
  if (!is.null(use_rep)) {
    features <- Embeddings(seurat_obj, reduction = use_rep)
  } else {
    features <- GetAssayData(seurat_obj, slot = "data")
    if (is(features, "dgCMatrix")) {
      features <- as.matrix(features)
    }
  }
  
  # Process n_layers parameter
  if (is.numeric(n_layers) && length(n_layers) == 1) {
    n_layers <- 0:n_layers
  }
  
  # Process aggregations parameter
  if (is.character(aggregations) && length(aggregations) == 1) {
    aggregations <- list(aggregations)
  }
  
  # Initialize aggregated features list
  aggregated_list <- list()
  
  # Add layer 0 (original features)
  if (0 %in% n_layers) {
    aggregated_list[["layer_0"]] <- features
  }
  
  # Process each layer
  for (layer in n_layers[n_layers > 0]) {
    # Compute adjacency matrix for current layer
    adj_layer <- adj_matrix
    if (layer > 1) {
      for (i in 2:layer) {
        # Use sparse matrix multiplication if possible
        adj_layer <- adj_layer %*% adj_matrix
      }
    }
    
    # Remove self-loops
    if (is(adj_layer, "dgCMatrix")) {
      # For sparse matrix, set diagonal to 0
      for (i in 1:nrow(adj_layer)) {
        # Find if there's a diagonal element
        idx <- which(adj_layer@i + 1 == i & adj_layer@p[i] < adj_layer@p[i+1])
        if (length(idx) > 0) {
          pos <- adj_layer@p[i] + idx
          adj_layer@x[pos] <- 0
        }
      }
      adj_layer <- drop0(adj_layer)
    } else {
      # For dense matrix, set diagonal to 0
      diag(adj_layer) <- 0
    }
    
    # Normalize (handle sparse and dense matrices)
    row_sums <- rowSums(adj_layer)
    row_sums[row_sums == 0] <- 1  # Avoid division by zero
    
    if (is(adj_layer, "dgCMatrix")) {
      # For sparse matrix, normalize by row sums
      adj_layer_normalized <- adj_layer
      # Iterate through each row
      for (i in 1:nrow(adj_layer)) {
        if (row_sums[i] > 0) {
          # Get the range of non-zero elements in this row
          start_idx <- adj_layer_normalized@p[i] + 1
          end_idx <- adj_layer_normalized@p[i + 1]
          # Check if there are any non-zero elements in this row
          if (start_idx <= end_idx) {
            # Normalize the non-zero elements
            adj_layer_normalized@x[start_idx:end_idx] <- adj_layer_normalized@x[start_idx:end_idx] / row_sums[i]
          }
        }
      }
    } else {
      # For dense matrix
      adj_layer_normalized <- adj_layer / row_sums
    }
    
    # Apply aggregations
    for (agg in aggregations) {
      if (agg == "mean") {
        aggregated <- adj_layer_normalized %*% features
        aggregated_list[[paste0("layer_", layer, "_mean")]] <- aggregated
      } else if (agg == "var") {
        # Calculate variance
        mean_agg <- adj_layer_normalized %*% features
        mean_sq_agg <- adj_layer_normalized %*% (features^2)
        var_agg <- mean_sq_agg - (mean_agg^2)
        aggregated_list[[paste0("layer_", layer, "_var")]] <- var_agg
      }
    }
  }
  
  # Concatenate all aggregated features
  if (length(aggregated_list) > 0) {
    aggregated_features <- do.call(cbind, aggregated_list)
    
    # Store in Seurat object
    seurat_obj[[out_key]] <- CreateAssayObject(counts = aggregated_features)
    DefaultAssay(seurat_obj) <- out_key
  }
  
  return(seurat_obj)
}

#' Cluster cells based on aggregated neighbor features
#' 
#' Similar to cellcharter's Cluster function, this function uses Gaussian Mixture
#' Models to cluster cells based on the aggregated neighbor features.
#' 
#' @param seurat_obj Seurat object with aggregated features
#' @param n_clusters Number of clusters to identify
#' @param use_assay Name of the assay containing aggregated features
#' @param out_col Name of the metadata column to store cluster labels
#' @param ... Additional parameters passed to Mclust
#' 
#' @return Seurat object with cluster labels
#' 
cluster_cells <- function(seurat_obj, n_clusters = 10, use_assay = "cellcharter", 
                         out_col = "cellcharter_cluster", ...) {
  # Get aggregated features
  if (!use_assay %in% names(seurat_obj@assays)) {
    stop(paste("Assay", use_assay, "not found. Run aggregate_neighbors first."))
  }
  
  features <- GetAssayData(seurat_obj, assay = use_assay, slot = "counts")
  
  # Run GMM clustering
  message(paste("Running Gaussian Mixture Model with", n_clusters, "clusters..."))
  gmm_fit <- Mclust(features, G = n_clusters, ...)
  
  # Add cluster labels to metadata
  seurat_obj@meta.data[[out_col]] <- factor(gmm_fit$classification)
  
  return(seurat_obj)
}

#' Analyze spatial domains
#' 
#' Analyzes spatial domains by identifying boundaries and computing shape metrics
#' 
#' @param seurat_obj Seurat object with cluster labels
#' @param cluster_col Name of the metadata column containing cluster labels
#' @param out_dir Directory to save output plots
#' 
#' @return Seurat object with additional spatial domain metrics
#' 
analyze_spatial_domains <- function(seurat_obj, cluster_col = "cellcharter_cluster", 
                                   out_dir = ".") {
  # Ensure output directory exists
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  # Plot spatial domains
  p <- SpatialDimPlot(seurat_obj, group.by = cluster_col, pt.size.factor = 1.5)
  ggsave(file.path(out_dir, "spatial_domains.png"), p, width = 10, height = 8)
  
  # Compute domain metrics
  clusters <- seurat_obj@meta.data[[cluster_col]]
  spatial_coord <- GetTissueCoordinates(seurat_obj)
  
  # Calculate domain sizes
  domain_sizes <- table(clusters)
  print("Domain sizes:")
  print(domain_sizes)
  
  # Add domain size to metadata
  seurat_obj@meta.data$domain_size <- as.numeric(domain_sizes[clusters])
  
  # Identify boundary cells
  message("Identifying boundary cells...")
  boundary_cells <- identify_boundary_cells(seurat_obj, cluster_col)
  seurat_obj@meta.data$boundary_cell <- boundary_cells
  
  # Plot boundary cells
  p_boundary <- SpatialDimPlot(seurat_obj, group.by = "boundary_cell", pt.size.factor = 1.5)
  ggsave(file.path(out_dir, "boundary_cells.png"), p_boundary, width = 10, height = 8)
  
  # Calculate domain shape metrics
  message("Calculating domain shape metrics...")
  shape_metrics <- calculate_shape_metrics(seurat_obj, cluster_col, spatial_coord)
  write.csv(shape_metrics, file.path(out_dir, "domain_shape_metrics.csv"))
  
  return(seurat_obj)
}

#' Identify boundary cells
#' 
#' Identifies cells at the boundary of spatial domains
#' 
#' @param seurat_obj Seurat object with cluster labels
#' @param cluster_col Name of the metadata column containing cluster labels
#' 
#' @return Logical vector indicating boundary cells
#' 
identify_boundary_cells <- function(seurat_obj, cluster_col) {
  # Get spatial neighbor graph
  if (!"spatial" %in% names(seurat_obj@graphs)) {
    seurat_obj <- FindSpatialNeighbors(seurat_obj, k = 20)
  }
  
  adj_matrix <- seurat_obj@graphs$spatial
  clusters <- seurat_obj@meta.data[[cluster_col]]
  
  # Convert to dense matrix if sparse
  if (is(adj_matrix, "dgCMatrix")) {
    adj_matrix <- as.matrix(adj_matrix)
  }
  
  # Identify boundary cells
  boundary <- logical(length(clusters))
  for (i in 1:length(clusters)) {
    neighbors <- which(adj_matrix[i, ] > 0)
    if (length(neighbors) > 0) {
      neighbor_clusters <- clusters[neighbors]
      boundary[i] <- any(neighbor_clusters != clusters[i])
    }
  }
  
  return(boundary)
}

#' Calculate domain shape metrics
#' 
#' Calculates shape metrics for each spatial domain
#' 
#' @param seurat_obj Seurat object with cluster labels
#' @param cluster_col Name of the metadata column containing cluster labels
#' @param spatial_coord Spatial coordinates of cells
#' 
#' @return Data frame with shape metrics for each domain
#' 
calculate_shape_metrics <- function(seurat_obj, cluster_col, spatial_coord) {
  clusters <- seurat_obj@meta.data[[cluster_col]]
  unique_clusters <- unique(clusters)
  
  metrics <- data.frame()
  
  for (cluster in unique_clusters) {
    # Get cells in this cluster
    cluster_cells <- which(clusters == cluster)
    cluster_coord <- spatial_coord[cluster_cells, ]
    
    # Calculate area (convex hull area)
    if (nrow(cluster_coord) > 2) {
      hull <- chull(cluster_coord)
      hull_coord <- cluster_coord[hull, ]
      area <- polygon_area(hull_coord[, 1], hull_coord[, 2])
    } else {
      area <- 0
    }
    
    # Calculate centroid
    centroid <- colMeans(cluster_coord)
    
    # Calculate elongation (eccentricity)
    if (nrow(cluster_coord) > 1) {
      cov_mat <- cov(cluster_coord)
      eig <- eigen(cov_mat)
      eigenvalues <- eig$values
      if (sum(eigenvalues) > 0) {
        elongation <- max(eigenvalues) / sum(eigenvalues)
      } else {
        elongation <- 0
      }
    } else {
      elongation <- 0
    }
    
    # Add to metrics
    metrics <- rbind(metrics, data.frame(
      cluster = cluster,
      size = length(cluster_cells),
      area = area,
      centroid_x = centroid[1],
      centroid_y = centroid[2],
      elongation = elongation
    ))
  }
  
  return(metrics)
}

#' Calculate polygon area
#' 
#' Calculates the area of a polygon given its coordinates
#' 
#' @param x X-coordinates
#' @param y Y-coordinates
#' 
#' @return Area of the polygon
#' 
polygon_area <- function(x, y) {
  n <- length(x)
  if (n < 3) return(0)
  
  area <- 0
  for (i in 1:n) {
    j <- ifelse(i == n, 1, i + 1)
    area <- area + (x[i] * y[j]) - (x[j] * y[i])
  }
  
  return(abs(area) / 2)
}

#' Main CellCharter workflow for Seurat objects
#' 
#' Runs the complete CellCharter workflow on a Seurat object
#' 
#' @param seurat_obj Seurat object with spatial data
#' @param n_layers Number of neighborhood layers to aggregate
#' @param n_clusters Number of spatial domains to identify
#' @param use_rep Dimensionality reduction to use
#' @param cell_type_col Name of the metadata column containing cell types
#' @param out_dir Output directory for results
#' @param x_coord_col Name of the metadata column containing x-coordinates
#' @param y_coord_col Name of the metadata column containing y-coordinates
#' 
#' @return Seurat object with spatial domain analysis
#' 
run_cellcharter <- function(seurat_obj, n_layers = 3, n_clusters = 10, 
                           use_rep = "pca", cell_type_col = "sub_cell_type",
                           out_dir = "./cellcharter_results", 
                           x_coord_col = "CenterX_global_px", 
                           y_coord_col = "CenterY_global_px") {
  # Create output directory
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  # Step 1: Aggregate neighbor features
  message("Aggregating neighbor features...")
  seurat_obj <- aggregate_neighbors(seurat_obj, 
                                   n_layers = n_layers, 
                                   use_rep = use_rep,
                                   x_coord_col = x_coord_col,
                                   y_coord_col = y_coord_col)
  
  # Step 2: Cluster cells
  message("Clustering cells...")
  seurat_obj <- cluster_cells(seurat_obj, 
                             n_clusters = n_clusters)
  
  # Step 3: Analyze spatial domains
  message("Analyzing spatial domains...")
  seurat_obj <- analyze_spatial_domains(seurat_obj, 
                                       out_dir = out_dir)
  
  # Step 4: Cell type enrichment in domains
  if (!is.null(cell_type_col) && cell_type_col %in% colnames(seurat_obj@meta.data)) {
    message("Analyzing cell type enrichment in domains...")
    
    # Create contingency table
    contingency_table <- table(seurat_obj@meta.data[[cell_type_col]], 
                              seurat_obj@meta.data$cellcharter_cluster)
    
    # Save contingency table
    write.csv(contingency_table, file.path(out_dir, "cell_type_enrichment.csv"))
    
    # Plot cell type distribution per domain
    p <- DimPlot(seurat_obj, group.by = cell_type_col, split.by = "cellcharter_cluster")
    ggsave(file.path(out_dir, "cell_type_by_domain.png"), p, width = 15, height = 10)
  }
  
  message("CellCharter workflow completed!")
  return(seurat_obj)
}

# Example usage for mPT Seurat object
#' Run CellCharter on mPT Seurat object
#' 
#' This function runs the complete CellCharter workflow on the mPT Seurat object
#' with detailed in metadata.
#' 
#' @param n_layers Number of neighborhood layers to aggregate (default: 3)
#' @param n_clusters Number of spatial domains to identify (default: 10)
#' @param use_rep Dimensionality reduction to use (default: "pca")
#' @param out_dir Output directory for results (default: "./cellcharter_results")
#' @param x_coord_col Name of the metadata column containing x-coordinates (default: "CenterX_global_px")
#' @param y_coord_col Name of the metadata column containing y-coordinates (default: "CenterY_global_px")
#' 
#' @return Seurat object with spatial domain analysis
#' 
run_cellcharter_mPT <- function(n_layers = 3, n_clusters = 10, 
                               use_rep = "pca", 
                               out_dir = "./cellcharter_results",
                               x_coord_col = "CenterX_global_px", 
                               y_coord_col = "CenterY_global_px") {
  # Check if mPT exists
  if (!exists("mPT")) {
    stop("mPT Seurat object not found. Please load it first.")
  }
  
  # Check if detailed exists in metadata
  if (!"detailed" %in% colnames(mPT@meta.data)) {
    stop("detailed not found in mPT metadata.")
  }
  
  # Check if spatial coordinates exist in metadata
  if (!x_coord_col %in% colnames(mPT@meta.data) || !y_coord_col %in% colnames(mPT@meta.data)) {
    stop(paste("Spatial coordinates not found in metadata. Please check", x_coord_col, "and", y_coord_col))
  }
  
  # Run CellCharter workflow
  mPT <- run_cellcharter(mPT, 
                        n_layers = n_layers, 
                        n_clusters = n_clusters, 
                        use_rep = use_rep, 
                        cell_type_col = "detailed",
                        out_dir = out_dir,
                        x_coord_col = x_coord_col,
                        y_coord_col = y_coord_col)
  
  return(mPT)
}

# Example usage
# Load required packages
# library(Seurat)
# library(Matrix)
# library(mclust)
# 
# Load your mPT Seurat object
# mPT <- readRDS("path/to/mPT.rds")
# 
# Run CellCharter
# mPT <- run_cellcharter_mPT()
# 
# Save the results
# saveRDS(mPT, "path/to/mPT_cellcharter.rds")

# Detailed usage example
#' Detailed usage example for CellCharter R
#' 
#' 1. Load required packages:
#'    library(Seurat)
#'    library(Matrix)
#'    library(mclust)
#'    source("cellcharter_R.R")
#' 
#' 2. Load your mPT Seurat object:
#'    mPT <- readRDS("path/to/mPT.rds")
#' 
#' 3. (Optional) Run dimensionality reduction if not already done:
#'    mPT <- NormalizeData(mPT)
#'    mPT <- FindVariableFeatures(mPT)
#'    mPT <- ScaleData(mPT)
#'    mPT <- RunPCA(mPT)
#' 
#' 4. Run CellCharter workflow:
#'    mPT <- run_cellcharter_mPT(n_layers = 3, n_clusters = 10)
#' 
#' 5. Explore the results:
#'    # View spatial domains
#'    SpatialDimPlot(mPT, group.by = "cellcharter_cluster")
#'    
#'    # View boundary cells
#'    SpatialDimPlot(mPT, group.by = "boundary_cell")
#'    
#'    # View cell type distribution in domains
#'    DimPlot(mPT, group.by = "detailed", split.by = "cellcharter_cluster")
#' 
#' 6. Save the results:
#'    saveRDS(mPT, "path/to/mPT_cellcharter.rds")
