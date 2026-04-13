# ============================
# CellCharter 结果全面可视化
# ============================

rm(list = ls())
options(future.globals.maxSize = 100 * 1024^3)

# 加载包
library(Seurat)
library(ggplot2)
library(patchwork)
library(pheatmap)
library(vegan)
library(dplyr)
library(tidyr)
library(RColorBrewer)

# 加载结果
message("=== 加载 CellCharter 结果 ===")
mPT <- readRDS("mPT_cellcharter_FINAL.rds")

# 确保邻域标签是因子
mPT$cellcharter_cluster <- as.factor(mPT$cellcharter_cluster)

# 获取邻域数量
n_clusters <- length(unique(mPT$cellcharter_cluster))
message(paste("邻域数量:", n_clusters))

# 检查是否有sample列
has_sample <- "sample" %in% colnames(mPT@meta.data)
if (!has_sample) {
  message("未找到sample列，将创建默认sample")
  mPT@meta.data$sample <- "all_samples"
}

# 创建19种颜色方案
if (n_clusters <= 12) {
  cluster_colors <- brewer.pal(n_clusters, "Set3")
} else {
  colors1 <- brewer.pal(12, "Set3")
  colors2 <- brewer.pal(8, "Pastel2")
  colors3 <- brewer.pal(9, "Pastel1")
  cluster_colors <- c(colors1, colors2, colors3)[1:n_clusters]
}

# 检查空间坐标
x_col <- "CenterX_global_px"
y_col <- "CenterY_global_px"
has_spatial <- x_col %in% colnames(mPT@meta.data) && y_col %in% colnames(mPT@meta.data)

if (has_spatial) {
  message("使用空间坐标进行可视化")
  coords <- mPT@meta.data[, c(x_col, y_col)]
  colnames(coords) <- c("COORD_1", "COORD_2")
  mPT[['spatial_coords']] <- CreateDimReducObject(
    embeddings = as.matrix(coords),
    key = "COORD_"
  )
  reduction_name <- "spatial_coords"
} else {
  message("未找到空间坐标，计算UMAP")
  if (!"umap" %in% names(mPT@reductions)) {
    mPT <- RunUMAP(mPT, dims = 1:50)
  }
  reduction_name <- "umap"
}

samples <- unique(mPT@meta.data$sample)

# =====================
# 1. 每个样本的邻域空间点图
# =====================
message("=== 1. 生成每个样本的邻域空间点图 ===")

for (sample in samples) {
  message(paste("处理样本:", sample))
  
  sample_data <- mPT@meta.data[mPT@meta.data$sample == sample, ]
  sample_cells <- rownames(sample_data)
  temp_seurat <- subset(mPT, cells = sample_cells)
  
  p <- DimPlot(temp_seurat, 
             reduction = reduction_name,
             group.by = "cellcharter_cluster", 
             pt.size = 0.3,
             label = FALSE,
             cols = cluster_colors) +
    ggtitle(paste("Sample:", sample)) +
    theme(legend.position = "bottom")
  
  sample_filename <- paste0("spatial_domain_sample_", gsub(" ", "_", sample), ".pdf")
  ggsave(sample_filename, p, width = 10, height = 8, dpi = 300)
  message(paste("保存样本图:", sample_filename))
}

# =====================
# 2. 所有样本合并的空间点图
# =====================
message("=== 2. 生成所有样本合并的空间点图 ===")

if (has_spatial) {
  p_combined <- DimPlot(mPT, 
                       reduction = reduction_name,
                       group.by = "cellcharter_cluster", 
                       pt.size = 0.3,
                       label = FALSE,
                       cols = cluster_colors) +
    ggtitle("All Samples - Spatial Domains") +
    theme(legend.position = "bottom")
  
  ggsave("spatial_domain_all_samples.pdf", p_combined, width = 14, height = 10, dpi = 300)
  message("保存合并空间图: spatial_domain_all_samples.pdf")
}

# =====================
# 3. 每个邻域细胞组成热图
# =====================
message("=== 3. 生成邻域细胞组成热图 ===")

cluster_composition <- table(mPT@meta.data$cellcharter_cluster, 
                           mPT@meta.data$detailed)

cluster_composition_df <- as.data.frame(cluster_composition)
colnames(cluster_composition_df) <- c("Cluster", "CellType", "Count")

cluster_totals <- aggregate(Count ~ Cluster, data = cluster_composition_df, sum)

cluster_composition_df <- merge(cluster_composition_df, cluster_totals, by = "Cluster")
cluster_composition_df$Proportion <- cluster_composition_df$Count.x / cluster_composition_df$Count.y

heatmap_data <- cluster_composition_df %>%
  select(Cluster, CellType, Proportion) %>%
  pivot_wider(names_from = CellType, values_from = Proportion, values_fill = 0)

cluster_numeric <- as.numeric(as.character(heatmap_data$Cluster))
heatmap_data <- heatmap_data[order(cluster_numeric), ]
row_labels <- paste0("CN", cluster_numeric[order(cluster_numeric)])
rownames(heatmap_data) <- row_labels
heatmap_matrix <- as.matrix(heatmap_data[, -1])

p_heatmap <- pheatmap(heatmap_matrix,
                    cluster_rows = TRUE,
                    cluster_cols = TRUE,
                    display_numbers = FALSE,
                    main = "Cell Type Composition by Spatial Domain",
                    fontsize = 10,
                    fontsize_row = 10,
                    labels_row = row_labels,
                    color = colorRampPalette(c("#f0f0f0", "#3182bd", "#08519c"))(100))

pdf("cell_type_composition_heatmap.pdf", width = 14, height = 10)
print(p_heatmap)
dev.off()
message("保存热图: cell_type_composition_heatmap.pdf")

# =====================
# 4. Shannon多样性指数小提琴图
# =====================
message("=== 4. 生成Shannon多样性指数小提琴图 ===")

shannon_diversity <- data.frame()

for (sample in samples) {
  sample_cells <- mPT@meta.data[mPT@meta.data$sample == sample, ]
  for (cluster in unique(sample_cells$cellcharter_cluster)) {
    cluster_cells <- sample_cells[sample_cells$cellcharter_cluster == cluster, ]
    cell_type_counts <- table(cluster_cells$detailed)
    
    if (length(cell_type_counts) > 1) {
      proportions <- cell_type_counts / sum(cell_type_counts)
      shannon <- -sum(proportions * log(proportions))
    } else {
      shannon <- 0
    }
    
    shannon_diversity <- rbind(shannon_diversity, 
                              data.frame(sample = sample, 
                                        cluster = as.numeric(as.character(cluster)), 
                                        shannon = shannon))
  }
}

shannon_diversity$cluster <- factor(shannon_diversity$cluster, 
                                     levels = sort(unique(shannon_diversity$cluster)))

n_clusters_shannon <- length(unique(shannon_diversity$cluster))
light_colors <- cluster_colors[1:n_clusters_shannon]
names(light_colors) <- levels(shannon_diversity$cluster)

p_shannon <- ggplot(shannon_diversity, aes(x = cluster, y = shannon, fill = cluster)) +
  geom_violin(alpha = 0.7, trim = TRUE, linewidth = 0.5, color = "black") +
  geom_jitter(width = 0.2, height = 0, size = 0.8, alpha = 0.6, color = "black") +
  scale_fill_manual(values = light_colors, name = "Spatial Domain") +
  labs(title = "Shannon Diversity Index by Spatial Domain",
       x = "Spatial Domain (CN)",
       y = "Shannon Diversity Index") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.position = "none",  # 移除图例
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_line(color = "gray95")
  )

# 保存Shannon多样性图
ggsave("shannon_diversity.pdf", p_shannon, width = 14, height = 6, dpi = 300)  # 调整高度使其更扁
message("保存Shannon多样性图: shannon_diversity.pdf")

# =====================
# 5. 每个样本中邻域组成饼图（合并，每行2个样本）
# =====================
message("=== 5. 生成合并的邻域组成饼图 ===")

all_clusters <- sort(as.numeric(unique(mPT@meta.data$cellcharter_cluster)))
all_pie_data <- data.frame()

for (sample in samples) {
  sample_data <- mPT@meta.data[mPT@meta.data$sample == sample, ]
  cluster_counts <- table(sample_data$cellcharter_cluster)
  
  for (cl in all_clusters) {
    count <- ifelse(as.character(cl) %in% names(cluster_counts), 
                    cluster_counts[as.character(cl)], 0)
    all_pie_data <- rbind(all_pie_data, data.frame(
      sample = sample,
      cluster = cl,
      count = count
    ))
  }
}

all_pie_data <- all_pie_data %>%
  group_by(sample) %>%
  mutate(proportion = count / sum(count))

all_pie_data$cluster <- factor(all_pie_data$cluster, levels = all_clusters)

pie_colors <- cluster_colors
names(pie_colors) <- all_clusters

ncol_pie <- 2
nrow_pie <- ceiling(length(samples) / ncol_pie)
pie_height <- max(6, nrow_pie * 4)

p_pie <- ggplot(all_pie_data, aes(x = "", y = proportion, fill = cluster)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  facet_wrap(~ sample, ncol = ncol_pie) +
  scale_fill_manual(values = pie_colors, name = "Spatial Domain (CN)") +
  theme_void() +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  ) +
  labs(title = "Spatial Domain Composition by Sample")

ggsave("sample_composition_pie_combined.pdf", p_pie, 
       width = 10, height = pie_height, dpi = 300)
message("保存合并饼图: sample_composition_pie_combined.pdf")

# =====================
# 6. 全部CN的UMAP图（参考）
# =====================
message("=== 6. 生成全部CN的UMAP图 ===")

if (!"umap" %in% names(mPT@reductions)) {
  message("计算UMAP...")
  mPT <- RunUMAP(mPT, dims = 1:50)
}

p_umap_all <- DimPlot(mPT, 
                     reduction = "umap",
                     group.by = "cellcharter_cluster", 
                     pt.size = 0.8, 
                     label = TRUE, 
                     label.size = 3,
                     cols = cluster_colors) +
  ggtitle("UMAP of Spatial Domains - All CNs") +
  theme(legend.position = "right")

ggsave("domain_umap_all.pdf", p_umap_all, width = 12, height = 10, dpi = 300)
message("保存全部CN的UMAP图: domain_umap_all.pdf")

# =====================
# 6. 每个CN单独高亮的UMAP图（合并）
# =====================
message("=== 6. 生成每个CN单独高亮的UMAP图 ===")

if (!"umap" %in% names(mPT@reductions)) {
  message("计算UMAP...")
  mPT <- RunUMAP(mPT, dims = 1:50)
}

# 获取UMAP坐标
umap_coords <- Embeddings(mPT, reduction = "umap")
umap_df <- data.frame(
  UMAP_1 = umap_coords[, 1],
  UMAP_2 = umap_coords[, 2],
  CN = mPT$cellcharter_cluster
)

# 获取所有CN编号
cn_numbers <- sort(as.numeric(unique(mPT$cellcharter_cluster)))
n_cn <- length(cn_numbers)

# 创建存储所有子图的列表
umap_plots <- list()

for (cn in cn_numbers) {
  # 创建高亮特定CN的标签
  umap_df$highlight <- ifelse(umap_df$CN == cn, as.character(umap_df$CN), "Other")
  
  # 颜色：该CN用对应颜色，其他用浅灰色
  highlight_colors <- c(cluster_colors[which(cn_numbers == cn)], "#E8E8E8")
  names(highlight_colors) <- c(as.character(cn), "Other")
  
  p <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = highlight)) +
    geom_point(size = 0.5, alpha = 0.7) +
    scale_color_manual(values = highlight_colors, name = "CN") +
    theme_minimal() +
    labs(title = paste("CN", cn, "highlighted"),
         x = "UMAP_1", y = "UMAP_2") +
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
      legend.position = "right",
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", color = NA)
    )
  
  umap_plots[[as.character(cn)]] <- p
}

# 将所有UMAP图合并到一张图上（4列5行，共20个位置，19个图）
umap_combined <- wrap_plots(umap_plots, ncol = 4, nrow = 5)

# 添加总标题
umap_combined <- umap_combined + 
  plot_annotation(
    title = "UMAP of Spatial Domains - Each CN Highlighted Separately",
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  )

# 保存合并的UMAP图
ggsave("domain_umap_combined.pdf", umap_combined, width = 16, height = 20, dpi = 300)
message("保存合并UMAP图: domain_umap_combined.pdf")

# =====================
# 7. 统计摘要
# =====================
message("=== 7. 生成统计摘要 ===")

cluster_stats <- mPT@meta.data %>%
  group_by(cellcharter_cluster) %>%
  summarise(
    n_cells = n(),
    n_boundary = sum(boundary_cell),
    boundary_ratio = sum(boundary_cell) / n()
  )

shannon_stats <- shannon_diversity %>%
  group_by(cluster) %>%
  summarise(
    mean_shannon = mean(shannon),
    sd_shannon = sd(shannon),
    min_shannon = min(shannon),
    max_shannon = max(shannon)
  )

write.csv(cluster_stats, "cluster_statistics.csv", row.names = FALSE)
write.csv(shannon_stats, "shannon_diversity_statistics.csv", row.names = FALSE)

message("=== 邻域统计摘要 ===")
print(cluster_stats)
message("=== Shannon多样性统计摘要 ===")
print(shannon_stats)

message("✅ 所有分析完成！")
message("生成的图表：")
message("1. spatial_domain_sample_*.pdf - 每个样本单独的邻域空间分布")
message("2. spatial_domain_all_samples.pdf - 所有样本合并的空间点图")
message("3. cell_type_composition_heatmap.pdf - 邻域细胞组成热图")
message("4. shannon_diversity.pdf - Shannon多样性指数小提琴图")
message("5. sample_composition_pie_combined.pdf - 合并的邻域组成饼图（每行2个）")
message("6. domain_umap_all.pdf - 全部CN的UMAP图")
message("7. cluster_statistics.csv - 邻域统计摘要")
message("8. shannon_diversity_statistics.csv - Shannon多样性统计摘要")