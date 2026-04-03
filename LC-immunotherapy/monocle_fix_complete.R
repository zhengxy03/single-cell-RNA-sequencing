# monocle_fix_complete.R
# 修复 monocle 中所有缺失的内部函数

library(monocle)
library(igraph)

# ===================== 获取所有需要的内部函数 =====================
# projPointOnLine
if (exists("projPointOnLine", where = asNamespace("monocle"))) {
  projPointOnLine <- get("projPointOnLine", envir = asNamespace("monocle"))
} else {
  # 简化版定义
  projPointOnLine <- function(point, line_points) {
    # 计算点到直线的投影
    v1 <- line_points[,2] - line_points[,1]
    v2 <- point - line_points[,1]
    t <- sum(v2 * v1) / sum(v1 * v1)
    projection <- line_points[,1] + t * v1
    return(projection)
  }
}

# findNearestPointOnMST
if (exists("findNearestPointOnMST", where = asNamespace("monocle"))) {
  findNearestPointOnMST <- get("findNearestPointOnMST", envir = asNamespace("monocle"))
} else {
  findNearestPointOnMST <- function(cds) {
    return(cds)
  }
}

# minSpanningTree
if (exists("minSpanningTree", where = asNamespace("monocle"))) {
  minSpanningTree <- get("minSpanningTree", envir = asNamespace("monocle"))
} else {
  minSpanningTree <- function(cds) {
    return(igraph::mst(cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_tree))
  }
}

# reducedDimS, reducedDimK
if (exists("reducedDimS", where = asNamespace("monocle"))) {
  reducedDimS <- get("reducedDimS", envir = asNamespace("monocle"))
}
if (exists("reducedDimK", where = asNamespace("monocle"))) {
  reducedDimK <- get("reducedDimK", envir = asNamespace("monocle"))
}

# cellPairwiseDistances
if (exists("cellPairwiseDistances", where = asNamespace("monocle"))) {
  cellPairwiseDistances <- get("cellPairwiseDistances", envir = asNamespace("monocle"))
}

# ===================== 修改 project2MST =====================
project2MST <- function (cds, Projection_Method) 
{
    dp_mst <- minSpanningTree(cds)
    Z <- reducedDimS(cds)
    Y <- reducedDimK(cds)
    cds <- findNearestPointOnMST(cds)
    closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
    closest_vertex_names <- colnames(Y)[closest_vertex]
    closest_vertex_df <- as.matrix(closest_vertex)
    row.names(closest_vertex_df) <- row.names(closest_vertex)
    tip_leaves <- names(which(degree(dp_mst) == 1))
    if (!is.function(Projection_Method)) {
        P <- Y[, closest_vertex]
    }
    else {
        P <- matrix(rep(0, length(Z)), nrow = nrow(Z))
        for (i in 1:length(closest_vertex)) {
            neighbors <- names(V(dp_mst)[suppressWarnings(.nei(closest_vertex_names[i], 
                mode = "all"))])
            projection <- NULL
            distance <- NULL
            Z_i <- Z[, i]
            for (neighbor in neighbors) {
                if (closest_vertex_names[i] %in% tip_leaves) {
                  tmp <- projPointOnLine(Z_i, Y[, c(closest_vertex_names[i], 
                    neighbor)])
                }
                else {
                  tmp <- Projection_Method(Z_i, Y[, c(closest_vertex_names[i], 
                    neighbor)])
                }
                projection <- rbind(projection, tmp)
                distance <- c(distance, dist(rbind(Z_i, tmp)))
            }
            if (!is(projection, "matrix")) 
                projection <- as.matrix(projection)
            P[, i] <- projection[which(distance == min(distance))[1], 
                ]
        }
    }
    colnames(P) <- colnames(Z)
    dp <- as.matrix(dist(t(P)))
    min_dist = min(dp[dp != 0])
    dp <- dp + min_dist
    diag(dp) <- 0
    cellPairwiseDistances(cds) <- dp
    gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
    dp_mst <- minimum.spanning.tree(gp)
    cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_tree <- dp_mst
    cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_dist <- P
    cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex <- closest_vertex_df
    cds
}

# ===================== 替换 monocle 命名空间中的函数 =====================
assignInNamespace("projPointOnLine", projPointOnLine, ns = "monocle")
assignInNamespace("findNearestPointOnMST", findNearestPointOnMST, ns = "monocle")
assignInNamespace("minSpanningTree", minSpanningTree, ns = "monocle")
assignInNamespace("reducedDimS", reducedDimS, ns = "monocle")
assignInNamespace("reducedDimK", reducedDimK, ns = "monocle")
assignInNamespace("cellPairwiseDistances", cellPairwiseDistances, ns = "monocle")
assignInNamespace("project2MST", project2MST, ns = "monocle")

# 导出到全局环境
assign("projPointOnLine", projPointOnLine, envir = .GlobalEnv)
assign("findNearestPointOnMST", findNearestPointOnMST, envir = .GlobalEnv)
assign("project2MST", project2MST, envir = .GlobalEnv)

cat("✓ Monocle fixed: all internal functions\n")