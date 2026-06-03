# metadata pretreatment
```R
obj <- readRDS("YA2025263-1_fin.rds")
# 查看每个 assay 的特征数
lapply(Assays(obj), function(assay) {
  length(rownames(obj[[assay]]))
})
library(Seurat)

library(ggplot2)
library(dplyr)
sample_tissue <- data.frame(
  sample = c("A1", "A2", "A3", "A4", "A5", 
             "B1", "B2", "B3", "B4", "B5",
             "C1", "C2", "C3", "C4", "C5",
             "D1", "D2", "D3", "D4", "D5"),
  tissue = c("NRT", "NRLN", "NRT", "NRLN", "NRT",
             "NRLN", "NRT", "NRLN", "RT", "RLN_N",
             "RLN_P", "RT", "RLN_N", "RLN_P", "RT",
             "RLN_N", "RLN_P", "RT", "RLN_N", "RLN_P")
)
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
                       #"Malignant cells" = "unknown",
                       "Basal cells" = "Malignant cells")

# 检查
table(obj$CellType)


#按照fov显示细胞分布
spatial_data <- obj@meta.data %>%
  select(CenterX_local_px, CenterY_local_px, CellType, fov, sample) %>%
  filter(!is.na(CenterX_local_px), !is.na(CenterY_local_px))

fov_list <- unique(spatial_data$fov)
for(f in fov_list) {
  data_fov <- spatial_data %>% filter(fov == f)
  
  p <- ggplot(data_fov, aes(x = CenterX_local_px, y = CenterY_local_px, color = CellType)) +
    geom_point(size = 1, alpha = 0.7) +
    coord_fixed() +
    theme_minimal() +
    labs(title = paste("FOV", f, "- Cell Types Distribution"),
         x = "X Coordinate", y = "Y Coordinate") +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5))
  
  print(p)

  ggsave(paste0("FOV_", f, "_celltypes.pdf"), p, width = 8, height = 6)
}

#按照样本显示细胞分布
library(ggplot2)
library(dplyr)
library(RColorBrewer)

spatial_data <- obj@meta.data %>%
  select(
    CenterX_global_px, CenterY_global_px,
    CellType, sample, fov,
    Area.um2, Circularity, AspectRatio, NucArea,
    Mean.PanCK, Mean.CD45, Mean.DAPI
  ) %>%
  filter(!is.na(CenterX_global_px), !is.na(CenterY_global_px))

all_cell_types <- unique(spatial_data$CellType)

cell_colors <- setNames(
  brewer.pal(min(length(all_cell_types), 12), "Set3"),
  all_cell_types[1:min(length(all_cell_types), 12)]
)

if(length(all_cell_types) > 12) {
  more_colors <- colorRampPalette(brewer.pal(8, "Pastel2"))(length(all_cell_types) - 12)
  names(more_colors) <- all_cell_types[13:length(all_cell_types)]
  cell_colors <- c(cell_colors, more_colors)
}

unique_samples <- unique(spatial_data$sample)
for(s in unique_samples) {
  sample_data <- spatial_data %>% filter(sample == s)
  
  # 使用"um"代替"μm"
  p <- ggplot(sample_data, 
              aes(x = CenterX_global_px, 
                  y = CenterY_global_px)) +
    geom_point(aes(color = CellType, size = Area.um2), 
               alpha = 0.8) +  # 稍微提高透明度让浅色更明显
    scale_size_continuous(range = c(0.5, 3)) +
    scale_color_manual(values = cell_colors) +  # 使用浅色系颜色
    coord_fixed() +
    theme_minimal() +
    labs(title = paste("Sample", s, "- Cell Types Distribution"),
         subtitle = paste("Total cells:", nrow(sample_data),
                         "| FOVs:", length(unique(sample_data$fov))),
         x = "X Coordinate (um)", 
         y = "Y Coordinate (um)",
         color = "Cell Type",
         size = "Cell Area (um²)") +
    theme(
      # 添加方框
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      # 白色背景
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      legend.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = "white", color = NA),
      # 黑色文字
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold", color = "black"),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "black"),
      legend.text = element_text(color = "black", size = 8),
      legend.title = element_text(color = "black", size = 9, face = "bold"),
      axis.title = element_text(color = "black", size = 10),
      axis.text = element_text(color = "black", size = 8),
      axis.line = element_line(color = "black", size = 0.3),
      axis.ticks = element_line(color = "black", size = 0.3),
      panel.grid = element_line(color = "gray90", size = 0.2),  # 浅灰色网格线
      legend.position = "right"
    )

  ggsave(paste0("sample_", s, "_distribution.pdf"), 
         p, width = 14, height = 10, bg = "white")
  
  print(p)
}

#只显示T 巨噬和树突
library(ggplot2)
library(dplyr)
library(RColorBrewer)

# 准备完整的空间数据
spatial_data <- obj@meta.data %>%
  select(
    CenterX_global_px, CenterY_global_px,
    CellType, sample, fov,
    Area.um2, Circularity, AspectRatio, NucArea,
    Mean.PanCK, Mean.CD45, Mean.DAPI
  ) %>%
  filter(!is.na(CenterX_global_px), !is.na(CenterY_global_px))

# 定义要显示的细胞类型（根据您的实际细胞类型名称调整）
target_cell_types <- c("T cells", "Macrophages", "Dendritic cells")

# 检查数据中实际存在的细胞类型
available_types <- intersect(target_cell_types, unique(spatial_data$CellType))
print(paste("Available cell types:", paste(available_types, collapse = ", ")))

# 为这三种细胞类型分配浅色系颜色
cell_colors <- c(
  "T cells" = "#c96a61",        # 浅粉色
  "Macrophages" = "#B3CDE3",    # 浅蓝色
  "Dendritic cells" = "#CCEBC5" # 浅绿色
)

# 只保留这三种细胞类型的数据
spatial_filtered <- spatial_data %>%
  filter(CellType %in% available_types)

# 获取所有样本
unique_samples <- unique(spatial_filtered$sample)

# 为每个样本创建PDF
for(s in unique_samples) {
  sample_data <- spatial_filtered %>% filter(sample == s)
  
  # 如果该样本没有目标细胞类型，跳过
  if(nrow(sample_data) == 0) {
    print(paste("Sample", s, "has no target cell types, skipping..."))
    next
  }
  
  # 打开PDF设备
  pdf(paste0("sample_", s, "_immune_cells.pdf"), width = 10, height = 8)
  
  # 创建图形
  p <- ggplot(sample_data, 
              aes(x = CenterX_global_px, 
                  y = CenterY_global_px)) +
    geom_point(aes(color = CellType, size = Area.um2), 
               alpha = 0.8) +
    scale_size_continuous(range = c(0.5, 4)) +
    scale_color_manual(values = cell_colors) +
    coord_fixed() +
    theme_minimal() +
    labs(title = paste("Sample", s, "- Immune Cells Distribution"),
         subtitle = paste("T cells:", sum(sample_data$CellType == "T cells"),
                         "| Macrophages:", sum(sample_data$CellType == "Macrophages"),
                         "| Dendritic cells:", sum(sample_data$CellType == "Dendritic cells")),
         x = "X Coordinate (μm)", 
         y = "Y Coordinate (μm)",
         color = "Cell Type",
         size = "Cell Area (μm²)") +
    theme(
      # 添加方框
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      # 白色背景
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      legend.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = "white", color = NA),
      # 黑色文字
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold", color = "black"),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "black"),
      legend.text = element_text(color = "black", size = 9),
      legend.title = element_text(color = "black", size = 10, face = "bold"),
      axis.title = element_text(color = "black", size = 10),
      axis.text = element_text(color = "black", size = 8),
      axis.line = element_line(color = "black", size = 0.3),
      axis.ticks = element_line(color = "black", size = 0.3),
      panel.grid = element_line(color = "gray90", size = 0.2),
      legend.position = "right"
    )
  
  print(p)
  dev.off()
  
  print(paste("PDF saved: sample_", s, "_immune_cells.pdf", sep = ""))
}

# 可选：创建一个包含所有样本的汇总PDF
p <- ggplot(spatial_filtered, 
            aes(x = CenterX_global_px, 
                y = CenterY_global_px)) +
  geom_point(aes(color = CellType, size = Area.um2), 
             alpha = 0.7) +
  scale_size_continuous(range = c(0.5, 3)) +
  scale_color_manual(values = cell_colors) +
  coord_fixed() +
  facet_wrap(~sample, ncol = 2) +  # 移除 scales = "free"
  theme_minimal() +
  labs(title = "Immune Cells Distribution Across Samples",
       x = "X Coordinate (μm)", 
       y = "Y Coordinate (μm)",
       color = "Cell Type",
       size = "Area (μm²)") +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    strip.text = element_text(size = 12, face = "bold"),
    strip.background = element_rect(fill = "gray95", color = "black"),
    legend.text = element_text(color = "black", size = 9),
    legend.title = element_text(color = "black", size = 10, face = "bold"),
    axis.title = element_text(color = "black", size = 9),
    axis.text = element_text(color = "black", size = 7),
    panel.grid = element_line(color = "gray95", size = 0.1),
    legend.position = "bottom",
    legend.box = "horizontal"
  )

# 保存为PDF
pdf("all_samples_immune_cells_combined.pdf", width = 14, height = 12)
print(p)
dev.off()

#按组织
library(ggplot2)
library(dplyr)
library(RColorBrewer)

# 准备完整的空间数据
spatial_data <- obj@meta.data %>%
  select(
    CenterX_global_px, CenterY_global_px,
    CellType, sample, fov, tissue,
    Area.um2, Circularity, AspectRatio, NucArea,
    Mean.PanCK, Mean.CD45, Mean.DAPI
  ) %>%
  filter(!is.na(CenterX_global_px), !is.na(CenterY_global_px))

# 定义要显示的细胞类型
target_cell_types <- c("T cells", "Macrophages", "Dendritic cells")
available_types <- intersect(target_cell_types, unique(spatial_data$CellType))
print(paste("Available cell types:", paste(available_types, collapse = ", ")))

# 为这三种细胞类型分配浅色系颜色
cell_colors <- c(
  "T cells" = "#FBB4AE",        # 浅粉色
  "Macrophages" = "#B3CDE3",    # 浅蓝色
  "Dendritic cells" = "#CCEBC5" # 浅绿色
)

# 只保留这三种细胞类型的数据
spatial_filtered <- spatial_data %>%
  filter(CellType %in% available_types) %>%
  filter(!is.na(tissue))

# 检查tissue类型
tissue_types <- unique(spatial_filtered$tissue)
print(paste("Tissue types:", paste(tissue_types, collapse = ", ")))

for(t in tissue_types) {
  tissue_data <- spatial_filtered %>% filter(tissue == t)
  
  if(nrow(tissue_data) == 0) next
  
  # 计算该tissue中的样本数
  samples_in_tissue <- unique(tissue_data$sample)
  n_samples <- length(samples_in_tissue)
  
  p <- ggplot(tissue_data,
              aes(x = CenterX_global_px,
                  y = CenterY_global_px)) +
    geom_point(aes(color = CellType, size = Area.um2),
               alpha = 0.7) +
    scale_size_continuous(range = c(0.5, 3)) +
    scale_color_manual(values = cell_colors) +
    facet_wrap(~sample, scales = "free", ncol = 3) +  # 添加 scales = "free"，每个子图有自己的坐标轴
    theme_minimal() +
    labs(title = paste("Tissue:", t, "- Immune Cells Distribution"),
         x = "X Coordinate (um)",
         y = "Y Coordinate (um)",
         color = "Cell Type",
         size = "Area (um²)") +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      strip.text = element_text(size = 10, face = "bold"),
      strip.background = element_rect(fill = "gray90"),
      legend.position = "bottom",
      legend.box = "horizontal",
      axis.text = element_text(size = 6),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_line(color = "gray95")
    ) +
    # 移除 coord_fixed()，因为它与 scales = "free" 冲突
    NULL
  
  # 保存每个tissue的PDF
  pdf(paste0("immune_cells_tissue_", gsub("/", "_", t), ".pdf"),
      width = max(20, n_samples * 3),
      height = 8)
  print(p)
  dev.off()
  
  print(paste("Saved: immune_cells_tissue_", t, ".pdf"))
}

library(tidyverse)

# 计算各亚群在不同组织中的占比
prop_data <- obj@meta.data %>%
  group_by(tissue, CellType) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(tissue) %>%
  mutate(proportion = count / sum(count) * 100) %>%
  # 筛选出RT和NRT，并指定顺序
  filter(tissue %in% c("RT", "NRT")) %>%
  # 将tissue转为因子，固定显示顺序
  mutate(tissue = factor(tissue, levels = c("RT", "NRT")))

# 获取CellType的颜色映射 - 修复版本
# 直接从meta.data中提取唯一的CellType和对应的颜色
celltype_color_df <- obj@meta.data %>%
  select(CellType, CellType_Colors) %>%
  distinct() %>%
  # 确保CellType_Colors是字符类型
  mutate(CellType_Colors = as.character(CellType_Colors))

# 创建命名的颜色向量
celltype_colors <- setNames(celltype_color_df$CellType_Colors, 
                           celltype_color_df$CellType)

# 确保prop_data中的CellType是字符类型
prop_data <- prop_data %>%
  mutate(CellType = as.character(CellType))

# 绘制堆叠条形图
p <- ggplot(prop_data, aes(x = tissue, y = proportion, fill = CellType)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = celltype_colors) +  # 使用CellType_Colors中的颜色
  labs(title = "Cell Type Distribution: RT vs NRT",
       x = "Tissue", y = "Proportion (%)", fill = "Cell Type") +
  theme_minimal() +
  theme(
    # 坐标轴线条和刻度
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.15, "cm"),
    # 坐标轴文本
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 10),
    axis.text.y = element_text(color = "black", size = 10),
    # 坐标轴标题
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    # 面板背景
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    # 图例
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )

# 查看图形
print(p)

# 保存为PDF
pdf("obj_celltype_distribution_RT_NRT.pdf", width = 4, height = 6)
print(p)
dev.off()

prop_data <- obj@meta.data %>%
  group_by(tissue, CellType) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(tissue) %>%
  mutate(proportion = count / sum(count) * 100) %>%
  # 筛选出RLN_P、RLN_N、NRLN，并指定顺序
  filter(tissue %in% c("RLN_P", "RLN_N", "NRLN")) %>%
  # 将tissue转为因子，固定显示顺序
  mutate(tissue = factor(tissue, levels = c("RLN_P", "RLN_N", "NRLN")))

# 获取CellType的颜色映射
celltype_color_df <- obj@meta.data %>%
  select(CellType, CellType_Colors) %>%
  distinct() %>%
  mutate(CellType_Colors = as.character(CellType_Colors))

# 创建命名的颜色向量
celltype_colors <- setNames(celltype_color_df$CellType_Colors, 
                           celltype_color_df$CellType)

# 确保prop_data中的CellType是字符类型
prop_data <- prop_data %>%
  mutate(CellType = as.character(CellType))

# 绘制堆叠条形图
p <- ggplot(prop_data, aes(x = tissue, y = proportion, fill = CellType)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = celltype_colors) +
  labs(title = "Cell Type Distribution by Tissue",
       x = "Tissue", y = "Proportion (%)", fill = "Cell Type") +
  theme_minimal() +
  theme(
    # 坐标轴线条和刻度
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.15, "cm"),
    # 坐标轴文本
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 10),
    axis.text.y = element_text(color = "black", size = 10),
    # 坐标轴标题
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    # 面板背景
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    # 图例
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )

# 查看图形
print(p)

# 保存为PDF
pdf("obj_celltype_distribution_RLN_P_RLN_N_NRLN.pdf", width = 6, height = 6)
print(p)
dev.off()

```
# t cells
```R
library(Seurat)
library(ggplot2)
library(dplyr)
obj <- readRDS("YA2025263-1_fin.rds")
#t_cells <- subset(obj,subset=CellType %in% c("T cells","NK cells"))
t_cells <- subset(obj,subset=CellType %in% c("T cells"))
#saveRDS(t_cells,file="t_NK.rds")
#t_cells <- subset(t_cells, subset = tissue %in% c("NRT","RT","NRLN","RLN_N","RLN_P"))
t_cells <- NormalizeData(t_cells)
t_cells <- FindVariableFeatures(t_cells, nfeatures = 3000)
hvgs <- VariableFeatures(t_cells)
t_cells <- ScaleData(t_cells, features = hvgs)
t_cells <- RunPCA(t_cells, features = hvgs, npcs = 50)
p <- ElbowPlot(t_cells, ndims = 30)
ggsave("p3.png",plot=p)

pca_var <- t_cells[["pca"]]@stdev^2 / sum(t_cells[["pca"]]@stdev^2)
cumsum_var <- cumsum(pca_var)

# 看达到多少 PC 时累计方差 > 70% / 80% / 90%
which(cumsum_var > 0.7)[1]   # 70%
which(cumsum_var > 0.8)[1]   # 80%
which(cumsum_var > 0.9)[1]   # 90%
#library(harmony)
#Malignant <- RunHarmony(Malignant,"tissue")


t_cells <- FindNeighbors(t_cells, dims = 1:32)
t_cells <- FindClusters(t_cells, resolution = 0.5)
t_cells <- RunUMAP(t_cells, dims = 1:32)
#saveRDS(t_cells,file="t_cells_new.rds")


t_cell_markers <- FindAllMarkers(t_cells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "wilcox")
t_significant_markers <- subset(t_cell_markers, p_val_adj < 0.05)
#write.csv(t_significant_markers, "t_all_marker.csv")
t_significant_markers <- t_significant_markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
write.csv(t_significant_markers, "t_top_marker_50.csv")


t_cell_markers <- FindAllMarkers(t_cells, only.pos = TRUE, min.pct = 0.01, logfc.threshold = 0.1, test.use = "wilcox")
t_significant_markers <- subset(t_cell_markers, p_val_adj < 0.05)
#write.csv(t_significant_markers, "t_all_marker.csv")
t_significant_markers <- t_significant_markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
write.csv(t_significant_markers, "t_top_marker_50.csv")



CD8 <- subset(t_cells, subset = CD4 < 0.5)
#CD8 <- subset(t_cells,subset=seurat_clusters %in% c(1,2,8,10))





genes_to_plot <- c("CD4","CD8A","CD8B","GZMK", "CTLA4", "FOXP3", "SELL", "HAVCR2","CD68")
p <- VlnPlot(t_cells, features = genes_to_plot, group.by = "seurat_clusters", pt.size = 0)
ggsave("p4.png",plot=p)
p <- FeaturePlot(t_cells, 
                  features = genes_to_plot, 
                  order = TRUE)
ggsave("p5.png",plot=p)
p <- DotPlot(t_cells, 
        features = genes_to_plot,
        group.by = "seurat_clusters") + 
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# ==================== 扩展版 T 细胞亚型基因集 ====================

# 1. Naive T cell（初始 T 细胞）
naive_genes <- c(
  "LEF1", "TCF7", "SELL", "CCR7", "IL7R", 
  "CD27", "CD28", "MAL", "S1PR1", "FOXP1",
  "BACH2", "ID3", "CD44", "LTB", "GPR183",
  "TXNIP", "NOSIP", "P2RY10", "RGS1", "EBF1",
  "TCF7L2", "MYB", "BCL6", "IL6ST", "TNS1"
)

# 2. Central memory T cell (Tcm，中心记忆)
tcm_genes <- c(
  "SELL", "CCR7", "IL7R", "TCF7", "LEF1",
  "CD27", "CD28", "LTB", "GPR183", "BACH2",
  "IL6ST", "MYB", "BCL2", "TBKBP1", "KLF2",
  "IL2RG", "CD3E", "TRAF3", "PPP2R2C", "GPR33"
)

# 3. Effector memory T cell (Tem，效应记忆)
tem_genes <- c(
  "GZMA", "GZMB", "GZMK", "PRF1", "NKG7",
  "IFNG", "CCL5", "CCL4", "CCL3", "KLRG1",
  "CX3CR1", "CXCR4", "CD44", "IL2RB", "EOMES",
  "TBX21", "ZEB2", "ID2", "FGFBP2", "GZMH",
  "GZMM", "GNLY", "FASLG", "TNF", "LTB"
)

# 4. Effector T cell (Teff，效应 T 细胞)
effector_genes <- c(
  "GZMA", "GZMB", "PRF1", "NKG7", "IFNG",
  "TNF", "CCL5", "CCL4", "CCL3", "XCL1", "XCL2",
  "TBX21", "EOMES", "ID2", "KLRG1", "CX3CR1",
  "FASLG", "GZMK", "GZMM", "GZMH", "GNLY",
  "CD44", "IL2RA", "IL2RB", "CD69", "HSPA1A"
)

# 5. Treg（调节 T 细胞）
treg_genes <- c(
  "FOXP3", "CTLA4", "IL2RA", "TIGIT", "IKZF2",
  "IKZF4", "LRRC32", "TNFRSF18", "TNFRSF4", "ENTPD1",
  "LAYN", "GITR", "CCR8", "CCR4", "CTLA4",
  "FOXP1", "BATF", "IL10", "TGFB1", "LAG3",
  "CD69", "CD25", "RTKN2", "FCRL3", "S100A4"
)

# 6. Exhausted T cell（耗竭 T 细胞）
exhausted_genes <- c(
  "PDCD1", "LAG3", "TIGIT", "TOX", "HAVCR2",
  "CTLA4", "ENTPD1", "CD160", "BTLA", "CD244",
  "LAYN", "CXCL13", "GZMB", "PRF1", "IFNG",
  "TBX21", "EOMES", "ID2", "ITGAE", "LILRB4",
  "RGS16", "PLAAT4", "ISG15", "IFI6", "BATF"
)

# 7. Th1（辅助 T 细胞 1 型）
th1_genes <- c(
  "TBX21", "IFNG", "STAT1", "STAT4", "IL12RB1",
  "IL12RB2", "CXCR3", "CXCL10", "CXCL11", "CCL3",
  "CCL4", "GZMA", "GZMB", "PRF1", "TNF",
  "IL18R1", "IL18RAP", "HLA-DRA", "HLA-DRB1", "CD38"
)

# 8. Th2（辅助 T 细胞 2 型）
th2_genes <- c(
  "GATA3", "IL4", "IL5", "IL13", "STAT6",
  "IL9", "IL10", "CCL11", "CCL17", "CCL22",
  "CCR3", "CCR4", "CCR8", "CRTH2", "IL1RL1",
  "IL17RB", "IL25", "IL33", "HLA-DR", "CD200R"
)

# 9. Th17（辅助 T 细胞 17 型）
th17_genes <- c(
  "RORC", "RORA", "IL17A", "IL17F", "IL21",
  "IL22", "IL23R", "CCR6", "CCL20", "CARD9",
  "STAT3", "IRF4", "BATF", "RBPJ", "AHR",
  "CSF2", "TGFB1", "IL6ST", "IL1R1", "IL12RB1"
)

# 10. Tfh（滤泡辅助 T 细胞）
tfh_genes <- c(
  "CXCR5", "BCL6", "ICOS", "PDCD1", "IL21",
  "IL6ST", "STAT3", "MAF", "BATF", "TOX",
  "CXCL13", "CCL19", "CD200", "CD40LG", "CD28",
  "SAP", "SH2D1A", "IL4", "IL10", "TNFRSF4"
)

# 11. Proliferating T cell（增殖 T 细胞）
proliferating_genes <- c(
  "MKI67", "TOP2A", "PCNA", "CCNB1", "CCNA2",
  "CDK1", "CDK2", "CDC20", "BIRC5", "AURKA",
  "AURKB", "PLK1", "TYMS", "RRM2", "ASPM",
  "CENPF", "CENPA", "KIF11", "KIF23", "KIF2C"
)

# 12. Cytotoxic T cell（细胞毒性评分，综合 GZMA/B/PRF1 等）
cytotoxic_genes <- c(
  "GZMA", "GZMB", "GZMH", "GZMK", "GZMM",
  "PRF1", "NKG7", "GNLY", "FASLG", "IFNG",
  "TNF", "CCL5", "CCL4", "CCL3", "XCL1", "XCL2"
)
activated_genes <- c(
  "CD69", "TNF", "IL23A", "CD5", "FOSL1", "JUNB", "JUND", 
  "IRF3", "BANK1", "S1PR4", "AKT1S1", "CCL19", "TRBV12-3",
  "GADD45B", "HSPA1B", "DNAJB1", "NLRP1", "DGKA", "AREL1",
  "CXCL12", "MADCAM1"
)

# 过滤存在的基因
activated_genes <- activated_genes[activated_genes %in% rownames(t_cells)]


# 过滤存在的基因（确保基因在数据中）
naive_genes <- naive_genes[naive_genes %in% rownames(t_cells)]
tcm_genes <- tcm_genes[tcm_genes %in% rownames(t_cells)]
tem_genes <- tem_genes[tem_genes %in% rownames(t_cells)]
treg_genes <- treg_genes[treg_genes %in% rownames(t_cells)]
exhausted_genes <- exhausted_genes[exhausted_genes %in% rownames(t_cells)]
effector_genes <- effector_genes[effector_genes %in% rownames(t_cells)]
proliferating_genes <- proliferating_genes[proliferating_genes %in% rownames(t_cells)]
th1_genes <- th1_genes[th1_genes %in% rownames(t_cells)]
th2_genes <- th2_genes[th2_genes %in% rownames(t_cells)]
th17_genes <- th17_genes[th17_genes %in% rownames(t_cells)]
tfh_genes <- tfh_genes[tfh_genes %in% rownames(t_cells)]
cytotoxic_genes <- cytotoxic_genes[cytotoxic_genes %in% rownames(t_cells)]

# 打印实际使用的基因数量
cat("Activated:", length(activated_genes), "\n")
cat("Naive:", length(naive_genes), "\n")
cat("Tcm:", length(tcm_genes), "\n")
cat("Tem:", length(tem_genes), "\n")
cat("Treg:", length(treg_genes), "\n")
cat("Exhausted:", length(exhausted_genes), "\n")
cat("Effector:", length(effector_genes), "\n")
cat("Proliferating:", length(proliferating_genes), "\n")
cat("Th1:", length(th1_genes), "\n")
cat("Th2:", length(th2_genes), "\n")
cat("Th17:", length(th17_genes), "\n")
cat("Tfh:", length(tfh_genes), "\n")
cat("Cytotoxic:", length(cytotoxic_genes), "\n")

# 打分
t_cells <- AddModuleScore(t_cells, features = list(naive_genes), name = "Naive_Score")
t_cells <- AddModuleScore(t_cells, features = list(tcm_genes), name = "Tcm_Score")
t_cells <- AddModuleScore(t_cells, features = list(tem_genes), name = "Tem_Score")
t_cells <- AddModuleScore(t_cells, features = list(treg_genes), name = "Treg_Score")
t_cells <- AddModuleScore(t_cells, features = list(exhausted_genes), name = "Exhausted_Score")
t_cells <- AddModuleScore(t_cells, features = list(effector_genes), name = "Effector_Score")
t_cells <- AddModuleScore(t_cells, features = list(proliferating_genes), name = "Proliferating_Score")
t_cells <- AddModuleScore(t_cells, features = list(th1_genes), name = "Th1_Score")
t_cells <- AddModuleScore(t_cells, features = list(th2_genes), name = "Th2_Score")
t_cells <- AddModuleScore(t_cells, features = list(th17_genes), name = "Th17_Score")
t_cells <- AddModuleScore(t_cells, features = list(tfh_genes), name = "Tfh_Score")
t_cells <- AddModuleScore(t_cells, features = list(cytotoxic_genes), name = "Cytotoxic_Score")
t_cells <- AddModuleScore(t_cells, features = list(activated_genes), name = "Activated_Score")

# 计算每个 cluster 的平均评分
score_cols <- c("Naive_Score1","Tcm_Score1","Tem_Score1", "Treg_Score1", "Exhausted_Score1", 
                "Effector_Score1", "Proliferating_Score1", "Th1_Score1",
                "Th2_Score1", "Th17_Score1", "Tfh_Score1", "Activated_Score1")

score_avg <- t_cells@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(across(all_of(score_cols), mean))

# 转换为矩阵用于热图
score_mat <- as.matrix(score_avg[, -1])
rownames(score_mat) <- score_avg$seurat_clusters

library(pheatmap)
# 热图
p <- pheatmap(score_mat, 
         scale = "row", 
         cluster_rows = TRUE, 
         cluster_cols = FALSE,
         main = "T cell subset scores by cluster",
         fontsize = 10,
         color = colorRampPalette(c("blue", "white", "red"))(100))
ggsave("t_score.png",plot=p)

t_cells<- subset(t_cells,subset=seurat_clusters %in% c(0,1,2,3,4,5,6,7))
t_cells <- RunUMAP(t_cells, dims = 1:32)

seurat_clusters <- as.character(unique(t_cells@meta.data$seurat_clusters))  # 转换为字符向量
num_legend_items <- length(seurat_clusters)  # 图例的个数
max_label_length <- max(nchar(seurat_clusters))  # 图例名称的最大长度

base_width <- 3000  # 基础宽度
base_height <- 3000  # 基础高度
legend_width_factor <- 100  # 每个图例项增加的宽度
label_length_factor <- 10  # 每个字符增加的宽度

dynamic_width <- base_width + (num_legend_items * legend_width_factor) + (max_label_length * label_length_factor)

library(RColorBrewer)
library(ggplot2)

if(num_legend_items <= 12) {
  light_colors <- brewer.pal(num_legend_items, "Set3")
} else {
  # 如果超过12个，用colorRampPalette扩展
  light_colors <- colorRampPalette(brewer.pal(12, "Set3"))(num_legend_items)
}


pdf("t_cluster.pdf", width = dynamic_width/300, height = base_height/300)

DimPlot(t_cells, reduction = "umap", label = TRUE, pt.size = 1, label.size = 8) +
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    ggtitle(NULL) +
    scale_color_manual(values = light_colors) +  # 使用浅色系颜色
    coord_fixed(ratio = 1) +
    guides(color = guide_legend(title = NULL, override.aes = list(size = 5))) +
    theme(
        text = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 28, color = "black"),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title.x = element_text(size = 36, face = "bold", color = "black"),
        axis.title.y = element_text(size = 36, face = "bold", color = "black", margin = margin(r = 20)),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 28, face = "bold", color = "black"),
        legend.title = element_text(size = 28, face = "bold", color = "black"),
        legend.position = "right",
        legend.box.margin = margin(0, 0, 0, 0),
        legend.key = element_blank(),
        legend.background = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5),
        aspect.ratio = 1,
        plot.margin = margin(10, 50, 10, 10)
    )

dev.off()



identity_mapping <- c(
  "0" = "T_Central-memory",
  "1" = "T_Activated",
  "2" = "T_Naive",
  "3" = "T_Effector-memory",
  "4" = "T_Naive/Memory",
  "5" = "T_Treg",
  "6" = "T_Naive",
  "7" = "T_Th1"
)
sub_cell_type <- identity_mapping[t_cells@meta.data$seurat_clusters]
t_cells@meta.data$sub_cell_type <- sub_cell_type

library(ggplot2)
library(ggsci)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(14)
cluster_colors <- c(
  # 浅色系在前（浅蓝第一个）
  "#A0CBE8",   # 1 浅蓝（你要的放第一个）
  "#FFBE7D",   # 2 浅橙
  "#8CD17D",   # 3 浅绿
  "#86BCB6",   # 4 浅青
  "#FF9D9A",   # 5 浅红
  "#FABFD2",   # 6 浅粉
  "#D4A6C8",   # 7 浅紫
  "#F1CE63",   # 8 淡黄
  "#D7B5A6",   # 9 浅棕
  "#B07AA1",   # 10 紫灰
  "#BAB0AC",   # 11 浅灰
  
  # 深色系在后
  "#4E79A7",   # 12 深蓝
  "#F28E2B",   # 13 橙
  "#59A14F",   # 14 深绿
  "#E15759",   # 15 红
  "#499894",   # 16 青绿
  "#D37295",   # 17 玫红
  "#B6992D",   # 18 土黄
  "#9D7660",   # 19 棕
  "#79706E"    # 20 深灰
)

cell_types <- as.character(unique(t_cells@meta.data$sub_cell_type))
num_legend_items <- length(cell_types)  # 图例的个数
max_label_length <- max(nchar(cell_types))  # 图例名称的最大长度

base_width <- 3000  # 基础宽度
base_height <- 3000  # 基础高度
legend_width_factor <- 100  # 每个图例项增加的宽度
label_length_factor <- 10  # 每个字符增加的宽度
dynamic_width <- base_width + (num_legend_items * legend_width_factor) + (max_label_length * label_length_factor)
pdf("t_annotation.pdf", width = dynamic_width/300, height = base_height/300)
DimPlot(t_cells, reduction = "umap", label = TRUE, pt.size = 2, group.by = "sub_cell_type", label.size = 4) +
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    ggtitle(NULL) +
    scale_color_manual(values = cluster_colors) +
    coord_fixed(ratio = 1) +
    guides(color = guide_legend(title = NULL, override.aes = list(size = 5))) +
    theme(
        text = element_text(size = 8, face = "bold"),
        axis.text.x = element_text(size = 28, color = "black"),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title.x = element_text(size = 36, face = "bold", color = "black"),
        axis.title.y = element_text(size = 36, face = "bold", color = "black", margin = margin(r = 20)),  # 增加右侧间距
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 16, face = "bold", color = "black"),
        legend.title = element_text(size = 16, face = "bold", color = "black"),
        legend.position = "right",
        legend.box.margin = margin(0, 0, 0, 0),
        legend.key = element_blank(),
        legend.background = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5),
        aspect.ratio = 1,
        plot.margin = margin(10, 50, 10, 10)
    )
dev.off()
saveRDS(t_cells,file="t_anno.rds")

#比例
library(ggplot2)
library(dplyr)

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
t_cells$tissue <- sample_tissue$tissue[match(t_cells$sample, sample_tissue$sample)]


library(tidyverse)

# 计算各亚群在不同组织中的占比
prop_data <- t_cells@meta.data %>%
  group_by(tissue, sub_cell_type) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(tissue) %>%
  mutate(proportion = count / sum(count) * 100) %>%
  filter(tissue %in% c("mPT", "nmPT")) %>%
  mutate(tissue = factor(tissue, levels = c("mPT", "nmPT")))

library(ggplot2)
library(ggsci)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(7)
cluster_colors <- c(
  # 浅色系在前（浅蓝第一个）
  "#A0CBE8",   # 1 浅蓝（你要的放第一个）
  "#FFBE7D",   # 2 浅橙
  "#8CD17D",   # 3 浅绿
  "#86BCB6",   # 4 浅青
  "#FF9D9A",   # 5 浅红
  "#FABFD2",   # 6 浅粉
  "#D4A6C8",   # 7 浅紫
  "#F1CE63",   # 8 淡黄
  "#D7B5A6",   # 9 浅棕
  "#B07AA1",   # 10 紫灰
  "#BAB0AC",   # 11 浅灰
  
  # 深色系在后
  "#4E79A7",   # 12 深蓝
  "#F28E2B",   # 13 橙
  "#59A14F",   # 14 深绿
  "#E15759",   # 15 红
  "#499894",   # 16 青绿
  "#D37295",   # 17 玫红
  "#B6992D",   # 18 土黄
  "#9D7660",   # 19 棕
  "#79706E"    # 20 深灰
)


prop_data <- prop_data %>%
  mutate(sub_cell_type = as.character(sub_cell_type))

# 绘制堆叠条形图
p <- ggplot(prop_data, aes(x = tissue, y = proportion, fill = sub_cell_type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = cluster_colors) +
  labs(title = "Cell Type Distribution: mPT vs nmPT",
       x = "Tissue", y = "Proportion (%)", fill = "Cell Type") +
  theme_minimal() +
  theme(
    # 坐标轴线条和刻度
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.15, "cm"),
    # 坐标轴文本
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 10),
    axis.text.y = element_text(color = "black", size = 10),
    # 坐标轴标题
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    # 面板背景
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    # 图例
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )


# 保存为PDF
pdf("t_celltype_distribution_tumor.pdf", width = 4, height = 6)
print(p)
dev.off()

prop_data <- t_cells@meta.data %>%
  group_by(tissue, sub_cell_type) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(tissue) %>%
  mutate(proportion = count / sum(count) * 100) %>%
  filter(tissue %in% c("metLN", "negLN")) %>%
  mutate(tissue = factor(tissue, levels = c("metLN", "negLN")))

library(ggplot2)
library(ggsci)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(7)


prop_data <- prop_data %>%
  mutate(sub_cell_type = as.character(sub_cell_type))

# 绘制堆叠条形图
p <- ggplot(prop_data, aes(x = tissue, y = proportion, fill = sub_cell_type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = cluster_colors) +
  labs(title = "Cell Type Distribution: metLN vs negLN",
       x = "Tissue", y = "Proportion (%)", fill = "Cell Type") +
  theme_minimal() +
  theme(
    # 坐标轴线条和刻度
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.15, "cm"),
    # 坐标轴文本
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 10),
    axis.text.y = element_text(color = "black", size = 10),
    # 坐标轴标题
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    # 面板背景
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    # 图例
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )


# 保存为PDF
pdf("t_celltype_distribution_LN.pdf", width = 4, height = 6)
print(p)
dev.off()

#negLN
sample_tissue <- data.frame(
  sample = c("A1", "A2", "A3", "A4", "A5", 
             "B1", "B2", "B3", "B4", "B5",
             "C1", "C2", "C3", "C4", "C5",
             "D1", "D2", "D3", "D4", "D5"),
  tissue2 = c("NRT", "NRLN", "NRT", "NRLN", "NRT",
             "NRLN", "NRT", "NRLN", "RT", "RLN_N",
             "RLN_P", "RT", "RLN_N", "RLN_P", "RT",
             "RLN_N", "RLN_P", "RT", "RLN_N", "RLN_P")
)
t_cells$tissue2 <- sample_tissue$tissue2[match(t_cells$sample, sample_tissue$sample)]

```
# t sc
```R
obj <- readRDS("sc_obj_anno.rds")
t_cells <- subset(obj,subset=cell_type %in% c("T cells"))
library(Seurat)
library(ggplot2)
library(dplyr)
t_cells <- NormalizeData(t_cells)
t_cells <- FindVariableFeatures(t_cells, nfeatures = 3000)
hvgs <- VariableFeatures(t_cells)
t_cells <- ScaleData(t_cells, features = hvgs)
t_cells <- RunPCA(t_cells, features = hvgs, npcs = 50)
p <- ElbowPlot(t_cells, ndims = 30)
ggsave("p3.png",plot=p)

pca_var <- t_cells[["pca"]]@stdev^2 / sum(t_cells[["pca"]]@stdev^2)
cumsum_var <- cumsum(pca_var)

# 看达到多少 PC 时累计方差 > 70% / 80% / 90%
which(cumsum_var > 0.7)[1]   # 70%
which(cumsum_var > 0.8)[1]   # 80%
which(cumsum_var > 0.9)[1]   # 90%
library(harmony)
t_cells <- RunHarmony(t_cells,"dataset")


t_cells <- FindNeighbors(t_cells, reduction = "harmony",dims = 1:22)
t_cells <- FindClusters(t_cells, resolution = 0.3)
t_cells <- RunUMAP(t_cells, reduction = "harmony",dims = 1:22)

library(ggplot2)
library(ggsci)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(15)

seurat_clusters <- as.character(unique(t_cells@meta.data$seurat_clusters))  # 转换为字符向量
num_legend_items <- length(seurat_clusters)  # 图例的个数
max_label_length <- max(nchar(seurat_clusters))  # 图例名称的最大长度

base_width <- 3000  # 基础宽度
base_height <- 3000  # 基础高度
legend_width_factor <- 100  # 每个图例项增加的宽度
label_length_factor <- 10  # 每个字符增加的宽度
dynamic_width <- base_width + (num_legend_items * legend_width_factor) + (max_label_length * label_length_factor)
pdf("t_sc_clusters.pdf", width = dynamic_width/300, height = base_height/300)  # 转换为英寸

DimPlot(t_cells, reduction = "umap", label = TRUE, pt.size = 1, label.size = 8) +
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    ggtitle(NULL) +
    scale_color_manual(values = npg_extended) +
    coord_fixed(ratio = 1) +
    guides(color = guide_legend(title = NULL, override.aes = list(size = 5))) +
    theme(
        text = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 28, color = "black"),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title.x = element_text(size = 36, face = "bold", color = "black"),
        axis.title.y = element_text(size = 36, face = "bold", color = "black", margin = margin(r = 20)),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 28, face = "bold", color = "black"),
        legend.title = element_text(size = 28, face = "bold", color = "black"),
        legend.position = "right",
        legend.box.margin = margin(0, 0, 0, 0),
        legend.key = element_blank(),
        legend.background = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5),
        aspect.ratio = 1,
        plot.margin = margin(10, 50, 10, 10)
    )

dev.off()

p <- DimPlot(t_cells, reduction = "umap", group.by = "dataset", label = FALSE)
ggsave("t_orig.png",plot=p)


t_cell_markers <- FindAllMarkers(t_cells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "wilcox")
t_significant_markers <- subset(t_cell_markers, p_val_adj < 0.05)
#write.csv(t_significant_markers, "t_all_marker.csv")
t_significant_markers <- t_significant_markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
write.csv(t_significant_markers, "t_sc_top_marker_50.csv")


t_cells <- subset(t_cells, idents = c(1,2,3,4,7,8,11))
t_cells <- RunUMAP(t_cells, reduction = "harmony",dims = 1:22)

identity_mapping <- c(
  "1" = "T_Naive",
  "2" = "T_Effector-memory",
  "3" = "T_Treg",
  "4" = "T_Exhausted",
  "7" = "T_Tfh",
  "8" = "T_NKT",
  "11" = "T_Activated"
)
t_cells$seurat_clusters <- factor(t_cells$seurat_clusters)
sub_cell_type <- identity_mapping[t_cells@meta.data$seurat_clusters]
t_cells@meta.data$sub_cell_type <- sub_cell_type

library(ggplot2)
library(ggsci)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(14)
cluster_colors <- c(
  # 浅色系在前（浅蓝第一个）
  "#A0CBE8",   # 1 浅蓝（你要的放第一个）
  "#FFBE7D",   # 2 浅橙
  "#8CD17D",   # 3 浅绿
  "#86BCB6",   # 4 浅青
  "#FF9D9A",   # 5 浅红
  "#FABFD2",   # 6 浅粉
  "#D4A6C8",   # 7 浅紫
  "#F1CE63",   # 8 淡黄
  "#D7B5A6",   # 9 浅棕
  "#B07AA1",   # 10 紫灰
  "#BAB0AC",   # 11 浅灰
  
  # 深色系在后
  "#4E79A7",   # 12 深蓝
  "#F28E2B",   # 13 橙
  "#59A14F",   # 14 深绿
  "#E15759",   # 15 红
  "#499894",   # 16 青绿
  "#D37295",   # 17 玫红
  "#B6992D",   # 18 土黄
  "#9D7660",   # 19 棕
  "#79706E"    # 20 深灰
)

cell_types <- as.character(unique(t_cells@meta.data$sub_cell_type))
num_legend_items <- length(cell_types)  # 图例的个数
max_label_length <- max(nchar(cell_types))  # 图例名称的最大长度

base_width <- 3000  # 基础宽度
base_height <- 3000  # 基础高度
legend_width_factor <- 100  # 每个图例项增加的宽度
label_length_factor <- 10  # 每个字符增加的宽度
dynamic_width <- base_width + (num_legend_items * legend_width_factor) + (max_label_length * label_length_factor)
pdf("t_sc_annotation.pdf", width = dynamic_width/300, height = base_height/300)
DimPlot(t_cells, reduction = "umap", label = TRUE, pt.size = 1, group.by = "sub_cell_type", label.size = 4) +
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    ggtitle(NULL) +
    scale_color_manual(values = cluster_colors) +
    coord_fixed(ratio = 1) +
    guides(color = guide_legend(title = NULL, override.aes = list(size = 5))) +
    theme(
        text = element_text(size = 8, face = "bold"),
        axis.text.x = element_text(size = 28, color = "black"),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title.x = element_text(size = 36, face = "bold", color = "black"),
        axis.title.y = element_text(size = 36, face = "bold", color = "black", margin = margin(r = 20)),  # 增加右侧间距
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 16, face = "bold", color = "black"),
        legend.title = element_text(size = 16, face = "bold", color = "black"),
        legend.position = "right",
        legend.box.margin = margin(0, 0, 0, 0),
        legend.key = element_blank(),
        legend.background = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5),
        aspect.ratio = 1,
        plot.margin = margin(10, 50, 10, 10)
    )
dev.off()
saveRDS(t_cells,file="t_sc_anno.rds")
```
# t mapping
```R
#heatmap
library(pheatmap)
library(ggplot2)
library(dplyr)

# ==================== 1. 获取单细胞 t_cells_sc 的 top 200 基因 ====================
cat("正在计算 t_cells_sc 的标记基因...\n")
t_cells_sc_markers <- FindAllMarkers(t_cells_sc, 
                                        only.pos = TRUE, 
                                        min.pct = 0.25, 
                                        logfc.threshold = 0.25, 
                                        test.use = "wilcox")
t_cells_sc_significant <- subset(t_cells_sc_markers, p_val_adj < 0.05)
t_cells_sc_top <- t_cells_sc_significant %>% 
                    group_by(cluster) %>% 
                    top_n(n = 500, wt = avg_log2FC)
sc_genes <- unique(t_cells_sc_top$gene)
cat("单细胞 top 200 基因数:", length(sc_genes), "\n")

# ==================== 2. 获取空间 t_cells 的 top 200 基因 ====================
cat("正在计算 t_cells 的标记基因...\n")
t_cells_markers <- FindAllMarkers(t_cells, 
                                     only.pos = TRUE, 
                                     min.pct = 0.25, 
                                     logfc.threshold = 0.25, 
                                     test.use = "wilcox")
t_cells_significant <- subset(t_cells_markers, p_val_adj < 0.05)
t_cells_top <- t_cells_significant %>% 
                 group_by(cluster) %>% 
                 top_n(n = 500, wt = avg_log2FC)
sp_genes <- unique(t_cells_top$gene)
cat("空间 top 200 基因数:", length(sp_genes), "\n")

# ==================== 3. 计算重叠基因 ====================
overlap.genes <- intersect(sc_genes, sp_genes)
cat("重叠基因数:", length(overlap.genes), "\n")

# 如果重叠基因太少，给出警告
if (length(overlap.genes) < 50) {
  cat("警告：重叠基因较少，建议增大 top_n 或使用高变基因\n")
}

# ==================== 4. 获取细胞类型 ====================
sc_types <- unique(t_cells_sc$sub_cell_type)
sp_types <- unique(t_cells$sub_cell_type)
cat("单细胞类型数:", length(sc_types), "\n")
cat("空间类型数:", length(sp_types), "\n")

# ==================== 5. 获取表达数据 ====================
# 使用 counts layer 并做 log1p 转换
expr_sc <- GetAssayData(t_cells_sc, assay = "RNA", layer = "counts")[overlap.genes, , drop = FALSE]
expr_sp <- GetAssayData(t_cells, assay = "RNA", layer = "counts")[overlap.genes, , drop = FALSE]

expr_sc <- log1p(expr_sc)
expr_sp <- log1p(expr_sp)

min_cells <- 5

# ==================== 6. 计算平均表达谱 ====================
# 单细胞
sc_avg <- matrix(NA, nrow = length(sc_types), ncol = length(overlap.genes))
rownames(sc_avg) <- sc_types
colnames(sc_avg) <- overlap.genes
for (i in seq_along(sc_types)) {
  ct <- sc_types[i]
  cells <- colnames(t_cells_sc)[t_cells_sc$sub_cell_type == ct]
  if (length(cells) >= min_cells) {
    sc_avg[i, ] <- rowMeans(expr_sc[, cells, drop = FALSE])
  } else {
    cat("单细胞类型", ct, "细胞数不足", min_cells, "，跳过\n")
  }
}

# 空间
sp_avg <- matrix(NA, nrow = length(sp_types), ncol = length(overlap.genes))
rownames(sp_avg) <- sp_types
colnames(sp_avg) <- overlap.genes
for (i in seq_along(sp_types)) {
  ct <- sp_types[i]
  cells <- colnames(t_cells)[t_cells$sub_cell_type == ct]
  if (length(cells) >= min_cells) {
    sp_avg[i, ] <- rowMeans(expr_sp[, cells, drop = FALSE])
  } else {
    cat("空间类型", ct, "细胞数不足", min_cells, "，跳过\n")
  }
}

# 移除全NA的行
sc_avg <- sc_avg[complete.cases(sc_avg), , drop = FALSE]
sp_avg <- sp_avg[complete.cases(sp_avg), , drop = FALSE]
cat("有效单细胞类型数:", nrow(sc_avg), "\n")
cat("有效空间类型数:", nrow(sp_avg), "\n")

# ==================== 7. 计算相关性矩阵 ====================
cor_matrix <- matrix(NA, nrow = nrow(sc_avg), ncol = nrow(sp_avg),
                     dimnames = list(rownames(sc_avg), rownames(sp_avg)))
for (i in 1:nrow(sc_avg)) {
  for (j in 1:nrow(sp_avg)) {
    cor_matrix[i, j] <- cor(sc_avg[i, ], sp_avg[j, ], method = "spearman", use = "complete.obs")
  }
}

cat("相关性范围:", range(cor_matrix, na.rm = TRUE), "\n")
plot_width <- 14   # 增加宽度
plot_height <- 12  # 增加高度
# ==================== 8. 绘制热图 ====================
pdf("t_cells_sc_subcelltype_correlation.pdf", width = plot_width, height = plot_height)

# 计算需要的边距：左边距（行名）+ 右边距（图例）+ 底部边距（列名）+ 顶部边距（标题）
# 使用 cellwidth 和 cellheight 控制每个格子的大小，使热图变小
pheatmap(cor_matrix,
         main = paste("Spearman correlation (", length(overlap.genes), " genes)", sep = ""),
         fontsize = 12,                    # 增大基础字体（原 8 → 12）
         fontsize_row = 10,                # 行名大小（原 8 → 10）
         fontsize_col = 10,                # 列名大小（原 8 → 10）
         fontsize_number = 6,              # 数字大小（如果显示）
         color = colorRampPalette(c("blue", "white", "red"))(100),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         display_numbers = FALSE,
         # 关键：控制格子大小，让热图变小
         cellwidth = 18,                   # 每个格子宽度（像素）
         cellheight = 18,                  # 每个格子高度（像素）
         # 调整边距（增加底部和左边空间给长标签）
         margins = c(12, 12))              # 底部边距12，左边边距12

dev.off()


```
# exhausted
```R
# 耗竭 T 细胞核心基因集（Tex，Exhausted T cells）
exhausted_genes <- c(
  # 免疫检查点分子
  "PDCD1", "LAG3", "TIGIT", "HAVCR2", "CTLA4",
  # 耗竭相关转录因子
  "TOX", "TOX2", "EOMES", "PRDM1", "IRF4", "BATF",
  # 共抑制分子
  "CD160", "BTLA", "CD244", "ENTPD1", "LAYN",
  # 效应分子（在耗竭中中高水平）
  "GZMB", "GZMK", "PRF1", "IFNG", "TNF",
  # 趋化因子
  "CXCL13", "CCL3", "CCL4", "CCL5",
  # 其他耗竭相关
  "RGS16", "PLAAT4", "ISG15", "IFI6"
)

# 可选：更精简的核心耗竭基因集
exhausted_core_genes <- c("PDCD1", "LAG3", "TIGIT", "HAVCR2", "TOX", "CXCL13")

# 过滤存在的基因
exhausted_genes <- exhausted_genes[exhausted_genes %in% rownames(t_cells)]
cat("耗竭基因集中实际存在的基因数:", length(exhausted_genes), "\n")

# 计算耗竭评分
t_cells <- AddModuleScore(t_cells, 
                          features = list(exhausted_genes), 
                          name = "Exhausted_Score")

# 查看评分列名（通常为 "Exhausted_Score1"）
head(t_cells$Exhausted_Score1)

library(ggplot2)
library(ggpubr)
library(dplyr)

# 准备数据
plot_data <- t_cells@meta.data %>%
  filter(!is.na(tissue), !is.na(Exhausted_Score1)) %>%
  mutate(tissue = factor(tissue, levels = c("mPT", "nmPT", "metLN", "negLN")))

# 计算 y 轴上限
y_max <- max(plot_data$Exhausted_Score1, na.rm = TRUE)

# 只显示 4 个关键比较的横线位置（合理分布）
p1 <- y_max * 1.00   # mPT vs nmPT（较低）
p2 <- y_max * 1.08   # mPT vs metLN
p3 <- y_max * 1.16   # mPT vs negLN（最高）
p4 <- y_max * 1.00  # metLN vs negLN（较低）

# 绘图
p <- ggplot(plot_data, aes(x = tissue, y = Exhausted_Score1, fill = tissue)) +
  geom_violin(alpha = 0.7, trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.15, fill = "white", alpha = 0.8, outlier.shape = NA) +
  scale_fill_manual(values = c("mPT" = "#d25c5e", "nmPT" = "#377EB8", 
                                "metLN" = "#57a554", "negLN" = "#984EA3")) +
  labs(title = "T Cell Exhaustion Score Across Tissues",
       subtitle = "Kruskal-Wallis test, p < 0.0001; Dunn's post-hoc test with FDR correction",
       x = NULL, y = "Exhaustion Score") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 9, color = "gray30"),
    legend.position = "none"
  ) +
  # 只保留 4 个关键比较
  geom_signif(comparisons = list(c("mPT", "nmPT")),
              annotations = "***", y_position = p1,
              tip_length = 0.01, size = 0.4, textsize = 3.5) +
  geom_signif(comparisons = list(c("mPT", "metLN")),
              annotations = "****", y_position = p2,
              tip_length = 0.01, size = 0.4, textsize = 3.5) +
  geom_signif(comparisons = list(c("mPT", "negLN")),
              annotations = "****", y_position = p3,
              tip_length = 0.01, size = 0.4, textsize = 3.5) +
  geom_signif(comparisons = list(c("metLN", "negLN")),
              annotations = "****", y_position = p4,
              tip_length = 0.01, size = 0.4, textsize = 3.5)

# 保存，高度设为 5（比之前的 4 稍高一点，给横线留空间）
ggsave("exhaustion_score_violin.pdf", p, width = 6.5, height = 4, dpi = 300)

```
# CD8
```R
#CD8 <- subset(t_cells, subset = CD4 < 0.5)
#t_cells <- subset(t_cells, subset = tissue %in% c("NRT","RT","NRLN","RLN_N","RLN_P"))


CD8 <- NormalizeData(CD8)
CD8 <- FindVariableFeatures(CD8, nfeatures = 2000)
hvgs <- VariableFeatures(CD8)
CD8 <- ScaleData(CD8, features = hvgs)
CD8 <- RunPCA(CD8, features = hvgs, npcs = 20)
CD8 <- RunUMAP(CD8, dims = 1:20)
CD8 <- FindNeighbors(CD8, dims = 1:20)
CD8 <- FindClusters(CD8, resolution = 0.8)

prop_data <- CD8@meta.data %>%
  group_by(tissue, seurat_clusters) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(tissue) %>%
  mutate(proportion = count / sum(count) * 100)

# 堆叠条形图
p <- ggplot(prop_data, aes(x = tissue, y = proportion, fill = as.factor(seurat_clusters))) +
  geom_bar(stat = "identity", position = "stack") 
  labs(title = "Cluster Distribution by Tissue",
       x = "Tissue", y = "Proportion (%)", fill = "Cluster") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("p7.png",plot=p,bg="white")





CD8_cell_markers <- FindAllMarkers(CD8, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "wilcox")
CD8_significant_markers <- subset(CD8_cell_markers, p_val_adj < 0.05)
#write.csv(CD8_significant_markers, "CD8_all_marker.csv")
CD8_significant_markers <- CD8_significant_markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.csv(CD8_significant_markers, "CD8_top_marker_20.csv")

CD8 <- subset(CD8,subset=seurat_clusters %in% c(2,5,6,7))

CD8 <- subset(CD8,subset=seurat_clusters %in% c(0,1,2,3,4,5,7,8,9,10))
CD8 <- subset(CD8,subset=seurat_clusters %in% c(1,3,4,8,9,10))
CD8 <- subset(CD8,subset=seurat_clusters %in% c(1,2,3,4,8,9,10))

CD8_cell_markers <- FindAllMarkers(CD8, only.pos = TRUE, min.pct = 0.01, logfc.threshold = 0.1, test.use = "wilcox")
CD8_significant_markers <- subset(CD8_cell_markers, p_val_adj < 1)
write.csv(CD8_significant_markers, "CD8_all_marker_3.csv")
CD8_significant_markers <- CD8_significant_markers %>% group_by(cluster) %>% top_n(n = 200, wt = avg_log2FC)
write.csv(CD8_significant_markers, "CD8_top_marker_200_3.csv")

#identity_mapping <- c(
  "0" = "Activited_T",
  "1" = "Quiescent_T",
  "2" = "Activited_T",
  "3" = "Memory_T",
  "4" = "Effector_T",
  "5" = "Exhausted_T",
  "6" = "FRC",
  "7" = "Memory_T",
  "8" = "Exhausted_T",
  "9" = "Effector_T",
  "10" = "γδ_T"
)

identity_mapping <- c(
  "0" = "Quiescent_T",
  "1" = "Activited_T",
  "2" = "Memory_T_1",
  "3" = "Native_T",
  "4" = "Memory_T_2",
  "5" = "Exhausted_T_1",  
  "6" = "Effector_T_1",
  "7" = "Exhausted_T_2",
  "8" = "Effector_T_2",
  "9" = "γδ_T"
)


identity_mapping <- c(
  "0" = "Quiescent_T",
  "1" = "Central_Memory_T_1",
  "2" = "Effector_Memory_T",
  "3" = "Exhausted_T",
  "4" = "Terminal_Effector_T_1",
  "5" = "Terminal_Effector_T_2",  
  "6" = "Central_Memory_T_2",
  "7" = "NKT/γδ_T"
)


identity_mapping <- c(
  "0" = "Quiescent_T",
  "1" = "Memory_T",
  "2" = "Other_T",
  "3" = "Effector_T",
  "4" = "Effector_T",
  "5" = "Memory_T",  
  "6" = "Memory_T",
  "7" = "Memory_T",
  "8" = "Effector_T",
  "9" = "Memory_T",
  "10" = "Exhausted_T",
  "11" = "Effector_T",
  "12" = "Effector_T",
  "13" = "Effector_T"
)

identity_mapping <- c(
  "0" = "Effector_T",
  "1" = "Naive_T",
  "2" = "Memory_T",
  "3" = "Memory_T",
  "4" = "Effector_T",
  "5" = "Exhausted_T",  
  "6" = "Memory_T",
  "7" = "Effector_T",
  "8" = "Effector_T",
  "9" = "Memory_T",
  "10" = "Effector_T",
  "11" = "Other_T",
  "12" = "Other_T"
)
major_cell_type <- identity_mapping[CD8@meta.data$seurat_clusters]
CD8@meta.data$major_cell_type <- major_cell_type

prop_data <- CD8@meta.data %>%
  group_by(tissue, major_cell_type) %>%
  summarise(count = n()) %>%
  group_by(tissue) %>%
  mutate(proportion = count / sum(count) * 100)

p <- ggplot(prop_data, aes(x = tissue, y = proportion, fill = major_cell_type)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  labs(title = "T Cell Subset Proportions by Tissue",
       x = "Tissue", y = "Proportion (%)", fill = "Cell Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("p8.png",plot=p,bg="white")
# 1.2 分组条形图（更清楚看每个亚群）
p <- ggplot(prop_data, aes(x = tissue, y = proportion, fill = tissue)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~major_cell_type, scales = "free_y") +
  theme_minimal() +
  labs(title = "T Cell Subset Proportions by Tissue",
       x = "Tissue", y = "Proportion (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

#pie
library(tidyverse)
library(patchwork)
library(cowplot)
library(RColorBrewer)

# 计算各组织中不同细胞类型的比例
prop_data <- CD8@meta.data %>%
  group_by(tissue, major_cell_type) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(tissue) %>%
  mutate(proportion = count / sum(count) * 100) %>%
  # 筛选并排序组织
  filter(tissue %in% c("RLN_P", "RLN_N", "NRLN", "RT", "NRT")) %>%
  mutate(tissue = factor(tissue, levels = c("RLN_P", "RLN_N", "NRLN", "RT", "NRT")))

# 获取所有细胞类型
cell_types <- unique(prop_data$major_cell_type)
# 使用RColorBrewer的浅色系调色板
if(length(cell_types) <= 12) {
  celltype_colors <- setNames(brewer.pal(length(cell_types), "Set3"), cell_types)
} else {
  set.seed(42)
  library(colorspace)
  random_colors <- qualitative_hcl(length(cell_types), palette = "Pastel 3")
  celltype_colors <- setNames(random_colors, cell_types)
}

# 创建饼图函数（带图例）
create_pie_chart <- function(tissue_name, data, colors) {
  tissue_data <- data %>% filter(tissue == tissue_name)
  
  p <- ggplot(tissue_data, aes(x = "", y = proportion, fill = major_cell_type)) +
    geom_bar(stat = "identity", width = 1, color = "white", linewidth = 0.2) +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = colors, name = "Cell Type") +
    theme_void() +
    labs(title = tissue_name) +
    theme(
      plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
      legend.position = "right",  # 显示图例
      legend.title = element_text(size = 8, face = "bold"),
      legend.text = element_text(size = 7),
      legend.key.size = unit(0.4, "cm")
    )
  return(p)
}

# 第一行：RLN_P, RLN_N, NRLN
row1_tissues <- c("RLN_P", "RLN_N", "NRLN")
row1_plots <- map(row1_tissues, ~create_pie_chart(.x, prop_data, celltype_colors))

# 第二行：RT, NRT
row2_tissues <- c("RT", "NRT")
row2_plots <- map(row2_tissues, ~create_pie_chart(.x, prop_data, celltype_colors))

# 拼图：第一行3个，第二行2个
row1 <- wrap_plots(row1_plots, nrow = 1)
row2 <- wrap_plots(row2_plots, nrow = 1)

# 组合两行
combined_plot <- plot_grid(
  row1,
  row2,
  ncol = 1,
  rel_heights = c(1, 1)
) +
  plot_annotation(title = "T Cell Subset Distribution by Tissue") &
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

# 查看图形
print(combined_plot)

# 保存为PDF
ggsave("pie_charts_by_tissue_with_legend.pdf", combined_plot, width = 16, height = 8, bg = "white")

seurat_clusters <- as.character(unique(CD8@meta.data$seurat_clusters))  # 转换为字符向量
num_legend_items <- length(seurat_clusters)  # 图例的个数
max_label_length <- max(nchar(seurat_clusters))  # 图例名称的最大长度

base_width <- 3000  # 基础宽度
base_height <- 3000  # 基础高度
legend_width_factor <- 100  # 每个图例项增加的宽度
label_length_factor <- 10  # 每个字符增加的宽度

dynamic_width <- base_width + (num_legend_items * legend_width_factor) + (max_label_length * label_length_factor)

library(RColorBrewer)
library(ggplot2)

if(num_legend_items <= 12) {
  light_colors <- brewer.pal(num_legend_items, "Set3")
} else {
  # 如果超过12个，用colorRampPalette扩展
  light_colors <- colorRampPalette(brewer.pal(12, "Set3"))(num_legend_items)
}


pdf("CD8_cluster.pdf", width = dynamic_width/300, height = base_height/300)

DimPlot(CD8, reduction = "umap", label = TRUE, pt.size = 1, label.size = 8) +
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    ggtitle(NULL) +
    scale_color_manual(values = light_colors) +  # 使用浅色系颜色
    coord_fixed(ratio = 1) +
    guides(color = guide_legend(title = NULL, override.aes = list(size = 5))) +
    theme(
        text = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 28, color = "black"),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title.x = element_text(size = 36, face = "bold", color = "black"),
        axis.title.y = element_text(size = 36, face = "bold", color = "black", margin = margin(r = 20)),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 28, face = "bold", color = "black"),
        legend.title = element_text(size = 28, face = "bold", color = "black"),
        legend.position = "right",
        legend.box.margin = margin(0, 0, 0, 0),
        legend.key = element_blank(),
        legend.background = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5),
        aspect.ratio = 1,
        plot.margin = margin(10, 50, 10, 10)
    )

dev.off()
```
# memory
```R
#memory <- subset(CD8,subset=major_cell_type =="Memory_T")
memory <- subset(CD8,subset=seurat_clusters %in% c(1,2,6))
memory <- NormalizeData(memory)
memory <- FindVariableFeatures(memory, nfeatures = 2000)
hvgs <- VariableFeatures(memory)
memory <- ScaleData(memory, features = hvgs)
memory <- RunPCA(memory, features = hvgs, npcs = 20)
memory <- RunUMAP(memory, dims = 1:20)
memory <- FindNeighbors(memory, dims = 1:20)
memory <- FindClusters(memory, resolution = 0.8)

memory_cell_markers <- FindAllMarkers(memory, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "wilcox")
memory_significant_markers <- subset(memory_cell_markers, p_val_adj < 0.05)
# write.csv(memory_significant_markers, "memory_all_marker.csv")
memory_significant_markers <- memory_significant_markers %>% group_by(cluster) %>% top_n(n =50, wt = avg_log2FC)
write.csv(memory_significant_markers, "memory_top_marker_50.csv")

memory_cell_markers <- FindAllMarkers(memory, only.pos = TRUE, min.pct = 0.01, logfc.threshold = 0.1, test.use = "wilcox")
memory_significant_markers <- subset(memory_cell_markers, p_val_adj < 1)
write.csv(memory_significant_markers, "memory_all_marker.csv")
memory_significant_markers <- memory_significant_markers %>% group_by(cluster) %>% top_n(n = 200, wt = avg_log2FC)
write.csv(memory_significant_markers, "memory_top_marker_200.csv")

prop_data <- memory@meta.data %>%
  group_by(tissue, seurat_clusters) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(tissue) %>%
  mutate(proportion = count / sum(count) * 100)

# 堆叠条形图
p <- ggplot(prop_data, aes(x = tissue, y = proportion, fill = as.factor(seurat_clusters))) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_brewer(palette = "Set3") +  # 浅色系
  labs(title = "Cluster Distribution by Tissue",
       x = "Tissue", y = "Proportion (%)", fill = "Cluster") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("p9.png",plot=p,bg="white")

seurat_clusters <- as.character(unique(memory@meta.data$seurat_clusters))  # 转换为字符向量
num_legend_items <- length(seurat_clusters)  # 图例的个数
max_label_length <- max(nchar(seurat_clusters))  # 图例名称的最大长度

base_width <- 3000  # 基础宽度
base_height <- 3000  # 基础高度
legend_width_factor <- 100  # 每个图例项增加的宽度
label_length_factor <- 10  # 每个字符增加的宽度

dynamic_width <- base_width + (num_legend_items * legend_width_factor) + (max_label_length * label_length_factor)

library(RColorBrewer)
library(ggplot2)

if(num_legend_items <= 12) {
  light_colors <- brewer.pal(num_legend_items, "Set3")
} else {
  # 如果超过12个，用colorRampPalette扩展
  light_colors <- colorRampPalette(brewer.pal(12, "Set3"))(num_legend_items)
}


pdf("memory_cluster.pdf", width = dynamic_width/300, height = base_height/300)

DimPlot(memory, reduction = "umap", label = TRUE, pt.size = 1, label.size = 8) +
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    ggtitle(NULL) +
    scale_color_manual(values = light_colors) +  # 使用浅色系颜色
    coord_fixed(ratio = 1) +
    guides(color = guide_legend(title = NULL, override.aes = list(size = 5))) +
    theme(
        text = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 28, color = "black"),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title.x = element_text(size = 36, face = "bold", color = "black"),
        axis.title.y = element_text(size = 36, face = "bold", color = "black", margin = margin(r = 20)),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 28, face = "bold", color = "black"),
        legend.title = element_text(size = 28, face = "bold", color = "black"),
        legend.position = "right",
        legend.box.margin = margin(0, 0, 0, 0),
        legend.key = element_blank(),
        legend.background = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5),
        aspect.ratio = 1,
        plot.margin = margin(10, 50, 10, 10)
    )

dev.off()


identity_mapping <- c(
  "0" = "TCF7+CCR7+Tmem",
  "1" = "GZMK+CXCR6+Tmem",
  "2" = "LTB+ARHGDIB+Tmem",
  "3" = "CD6+DUSP2+Tmem"
)

sub_cell_type <- identity_mapping[memory@meta.data$seurat_clusters]
memory@meta.data$sub_cell_type <- sub_cell_type

prop_data <- memory@meta.data %>%
  group_by(tissue, sub_cell_type) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(tissue) %>%
  mutate(proportion = count / sum(count) * 100) %>%
  # 筛选出需要的3个组织，并指定顺序
  filter(tissue %in% c("RLN_P", "RLN_N", "NRLN")) %>%
  # 将tissue转为因子，固定显示顺序
  mutate(tissue = factor(tissue, levels = c("RLN_P", "RLN_N", "NRLN")))

# 绘制水平堆叠条形图
p <- ggplot(prop_data, aes(x = proportion, y = tissue, fill = as.factor(sub_cell_type))) +  # 交换x和y
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_brewer(palette = "Set3") +  # 浅色系
  labs(title = "Cluster Distribution by Tissue",
       x = "Proportion (%)", y = "Tissue", fill = "Cluster") +  # 交换x和y轴标签
  theme_minimal() +
  theme(
    # 坐标轴线条和刻度
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.15, "cm"),
    # 坐标轴文本 - 注意y轴文本现在不需要旋转了
    axis.text.x = element_text(color = "black", size = 10),
    axis.text.y = element_text(color = "black", size = 10),
    # 坐标轴标题
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    # 面板背景
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    # 图例
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )

# 查看图形
print(p)

# 保存为PDF
pdf("cluster_distribution_by_tissue_LN_horizontal.pdf", width = 8, height = 3)
print(p)
dev.off()

genes_to_plot <- c("CCL5","NKG7",  "TCF7", "CCR7", "CD6", "DUSP2", "LTB", "ARHGDIB")
VlnPlot(CD8, features = genes_to_plot, group.by = "seurat_clusters", pt.size = 0)

#bubble heat map
p <- DotPlot(memory, features = genes_to_plot, group.by = "sub_cell_type") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_color_gradient(low = "#00CED1", high = "#FF4500") 
ggsave("dotplot_memory_cell_type.pdf",plot=p,width=8,height=3,bg="white")
```
# 效应T
```R
CD8 <- readRDS("CD8_final.rds")

identity_mapping <- c(
  "0" = "Quiescent_T",
  "1" = "Memory_T",
  "2" = "Effector_T",
  "3" = "Exhausted_T",
  "4" = "Effector_T",
  "5" = "Effector_T",  
  "6" = "Memory_T",
  "7" = "Effector_T"
)

major_cell_type <- identity_mapping[CD8@meta.data$seurat_clusters]
CD8@meta.data$major_cell_type <- major_cell_type

prop_data <- CD8@meta.data %>%
  group_by(tissue, major_cell_type) %>%
  summarise(count = n()) %>%
  group_by(tissue) %>%
  mutate(proportion = count / sum(count) * 100)

ggplot(prop_data, aes(x = tissue, y = proportion, fill = major_cell_type)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  labs(title = "T Cell Subset Proportions by Tissue",
       x = "Tissue", y = "Proportion (%)", fill = "Cell Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

effector <- subset(CD8,subset=seurat_clusters %in% c(2,4,5,7))
effector <- NormalizeData(effector)
effector <- FindVariableFeatures(effector, nfeatures = 2000)

hvgs <- VariableFeatures(effector)

effector <- ScaleData(effector, features = hvgs)
effector <- RunPCA(effector, features = hvgs, npcs = 20)
effector <- RunUMAP(effector, dims = 1:20)
effector <- FindNeighbors(effector, dims = 1:20)
effector <- FindClusters(effector, resolution = 0.8)
DimPlot(effector, reduction = "umap", label = TRUE, pt.size = 1, label.size = 8)

prop_data <- effector@meta.data %>%
  group_by(tissue, seurat_clusters) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(tissue) %>%
  mutate(proportion = count / sum(count) * 100)

# 堆叠条形图
ggplot(prop_data, aes(x = tissue, y = proportion, fill = as.factor(seurat_clusters))) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_brewer(palette = "Set3") +  # 浅色系
  labs(title = "Cluster Distribution by Tissue",
       x = "Tissue", y = "Proportion (%)", fill = "Cluster") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

effector_cell_markers <- FindAllMarkers(effector, only.pos = TRUE, min.pct = 0.01, logfc.threshold = 0.1, test.use = "wilcox")
effector_significant_markers <- subset(effector_cell_markers, p_val_adj < 1)
write.csv(effector_significant_markers, "effector_all_marker.csv")
effector_significant_markers <- effector_significant_markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
write.csv(effector_significant_markers, "effector_top_marker_50.csv")

identity_mapping <- c(
  "0" = "COX1+Teff",
  "1" = "GZMK+Teff",
  "2" = "ITGAL+Teff",
  "3" = "CCL4+Teff",
  "4" = "GNLY+Teff",
  "5" = "TYROBP+Teff"
)

sub_cell_type <- identity_mapping[effector@meta.data$seurat_clusters]
effector@meta.data$sub_cell_type <- sub_cell_type


sub_cell_types <- as.character(unique(effector@meta.data$sub_cell_type))  # 转换为字符向量
num_legend_items <- length(sub_cell_types)  # 图例的个数
max_label_length <- max(nchar(sub_cell_types))  # 图例名称的最大长度

base_width <- 3000  # 基础宽度
base_height <- 3000  # 基础高度
legend_width_factor <- 100  # 每个图例项增加的宽度
label_length_factor <- 10  # 每个字符增加的宽度

dynamic_width <- base_width + (num_legend_items * legend_width_factor) + (max_label_length * label_length_factor)

library(RColorBrewer)
library(ggplot2)

if(num_legend_items <= 12) {
  light_colors <- brewer.pal(num_legend_items, "Set3")
} else {
  # 如果超过12个，用colorRampPalette扩展
  light_colors <- colorRampPalette(brewer.pal(12, "Set3"))(num_legend_items)
}


pdf("effector_anno.pdf", width = dynamic_width/300, height = base_height/300)

DimPlot(effector, reduction = "umap", label = TRUE, pt.size = 1, label.size = 8,group.by="sub_cell_type") +
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    ggtitle(NULL) +
    scale_color_manual(values = light_colors) +  # 使用浅色系颜色
    coord_fixed(ratio = 1) +
    guides(color = guide_legend(title = NULL, override.aes = list(size = 5))) +
    theme(
        text = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 28, color = "black"),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title.x = element_text(size = 36, face = "bold", color = "black"),
        axis.title.y = element_text(size = 36, face = "bold", color = "black", margin = margin(r = 20)),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 28, face = "bold", color = "black"),
        legend.title = element_text(size = 28, face = "bold", color = "black"),
        legend.position = "right",
        legend.box.margin = margin(0, 0, 0, 0),
        legend.key = element_blank(),
        legend.background = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5),
        aspect.ratio = 1,
        plot.margin = margin(10, 50, 10, 10)
    )

dev.off()
genes_to_plot <- c("COX1")
VlnPlot(effector, features = genes_to_plot, group.by = "seurat_clusters", pt.size = 0)

prop_data <- effector@meta.data %>%
  group_by(tissue, sub_cell_type) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(tissue) %>%
  mutate(proportion = count / sum(count) * 100) %>%
  # 筛选出需要的3个组织，并指定顺序
  filter(tissue %in% c("RLN_P", "RLN_N", "NRLN")) %>%
  # 将tissue转为因子，固定显示顺序
  mutate(tissue = factor(tissue, levels = c("RLN_P", "RLN_N", "NRLN")))

# 绘制水平堆叠条形图
p <- ggplot(prop_data, aes(x = proportion, y = tissue, fill = as.factor(sub_cell_type))) +  # 交换x和y
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_brewer(palette = "Set3") +  # 浅色系
  labs(title = "Cluster Distribution by Tissue",
       x = "Proportion (%)", y = "Tissue", fill = "Cluster") +  # 交换x和y轴标签
  theme_minimal() +
  theme(
    # 坐标轴线条和刻度
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.15, "cm"),
    # 坐标轴文本 - 注意y轴文本现在不需要旋转了
    axis.text.x = element_text(color = "black", size = 10),
    axis.text.y = element_text(color = "black", size = 10),
    # 坐标轴标题
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    # 面板背景
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    # 图例
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )

# 查看图形
print(p)

# 保存为PDF
pdf("effector_cluster_distribution_by_tissue_LN.pdf", width = 8, height = 3)
print(p)
dev.off()
```
# 空间
```R
library(ggplot2)
library(dplyr)
library(RColorBrewer)

# ===================== 1. 数据准备（完全复用你的列名） =====================
spatial_data <- memory@meta.data %>%
  select(
    CenterX_global_px, CenterY_global_px,
    sub_cell_type, sample, fov, tissue,  # 替换CellType为sub_cell_type
    Area.um2, Circularity, AspectRatio, NucArea,
    Mean.PanCK, Mean.CD45, Mean.DAPI
  ) %>%
  filter(!is.na(CenterX_global_px), !is.na(CenterY_global_px))

# ===================== 2. 定义目标亚群（替换为GZMK+CXCR6+Tmem） =====================
# 定义要显示的亚群（仅GZMK+CXCR6+Tmem）
target_subtypes <- c("GZMK+CXCR6+Tmem")
available_subtypes <- intersect(target_subtypes, unique(spatial_data$sub_cell_type))
print(paste("Available subtypes:", paste(available_subtypes, collapse = ", ")))

# 为目标亚群分配你指定的浅粉色（和你T cells颜色一致）
cell_colors <- c(
  "GZMK+CXCR6+Tmem" = "#FBB4AE"  # 浅粉色
)

# ===================== 3. 筛选数据（仅保留目标亚群+非空tissue） =====================
spatial_filtered <- spatial_data %>%
  filter(sub_cell_type %in% available_subtypes) %>%
  filter(!is.na(tissue))

# 检查tissue类型（和你代码逻辑一致）
tissue_types <- unique(spatial_filtered$tissue)
print(paste("Tissue types:", paste(tissue_types, collapse = ", ")))

# ===================== 4. 按tissue循环绘图（完全复刻你的逻辑） =====================
for(t in tissue_types) {
  tissue_data <- spatial_filtered %>% filter(tissue == t)
  
  # 跳过无数据的tissue
  if(nrow(tissue_data) == 0) {
    print(paste("Tissue", t, "has no GZMK+CXCR6+Tmem cells, skipping..."))
    next
  }
  
  # 计算该tissue中的样本数（用于适配PDF宽度）
  samples_in_tissue <- unique(tissue_data$sample)
  n_samples <- length(samples_in_tissue)
  
  # 绘图：100%复用你的参数（scales="free"、facet_wrap、主题等）
  p <- ggplot(tissue_data,
              aes(x = CenterX_global_px,
                  y = CenterY_global_px)) +
    geom_point(aes(color = sub_cell_type, size = Area.um2),  # 替换CellType为sub_cell_type
               alpha = 0.7) +
    scale_size_continuous(range = c(0.5, 3)) +  # 和你一致的大小范围
    scale_color_manual(values = cell_colors) +  # 专属浅粉色
    facet_wrap(~sample, scales = "free", ncol = 3) +  # 保留scales="free"，和你一致
    theme_minimal() +
    labs(title = paste("Tissue:", t, "- GZMK+CXCR6+Tmem Distribution"),  # 修改标题
         x = "X Coordinate (um)",
         y = "Y Coordinate (um)",
         color = "Subtype",  # 替换Cell Type为Subtype
         size = "Area (um²)") +
    theme(
      # 完全复用你的主题样式
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      strip.text = element_text(size = 10, face = "bold"),
      strip.background = element_rect(fill = "gray90"),
      legend.position = "bottom",
      legend.box = "horizontal",
      axis.text = element_text(size = 6),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_line(color = "gray95")
    ) +
    # 移除coord_fixed()（和你一样，避免与scales="free"冲突）
    NULL
  
  # 保存每个tissue的PDF（复用你的命名+宽度适配逻辑）
  pdf(paste0("GZMK+CXCR6+Tmem_tissue_", gsub("/", "_", t), ".pdf"),
      width = max(20, n_samples * 3),  # 适配样本数的宽度
      height = 8)
  print(p)
  dev.off()
  
  print(paste("Saved: GZMK+CXCR6+Tmem_tissue_", t, ".pdf"))
}

#全部细胞的空间
obj@meta.data$sub_cell_type <- NA
common_cells <- intersect(rownames(memory@meta.data), rownames(obj@meta.data))
obj@meta.data[common_cells, "sub_cell_type"] <- memory@meta.data[common_cells, "sub_cell_type"]

cat("成功注释的细胞数：", sum(!is.na(obj@meta.data$sub_cell_type)), "\n")
cat("obj中sub_cell_type列的唯一值：", paste(unique(obj@meta.data$sub_cell_type), collapse = ", "), "\n")

obj@meta.data$major_cell_type <- NA
common_cells <- intersect(rownames(CD8@meta.data), rownames(obj@meta.data))
obj@meta.data[common_cells, "major_cell_type"] <- CD8@meta.data[common_cells, "major_cell_type"]

cat("成功注释的细胞数：", sum(!is.na(obj@meta.data$major_cell_type)), "\n")
cat("obj中major_cell_type列的唯一值：", paste(unique(obj@meta.data$major_cell_type), collapse = ", "), "\n")

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

# 3. 优先级2：填充major_cell_type（sub_cell_type为NA的细胞，转字符避免数字）
major_cell_idx <- is.na(obj@meta.data$detailed) & !is.na(obj@meta.data$major_cell_type)
obj@meta.data$detailed[major_cell_idx] <- as.character(obj@meta.data$major_cell_type[major_cell_idx])
cat("优先级2（major_cell_type）填充：", sum(major_cell_idx), "个细胞\n")

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
print(head(obj@meta.data[, c("sub_cell_type", "major_cell_type", "CellType", "detailed")], 10))
cat("\n=== 修复后detailed列前15个类型 ===")
print(head(sort(table(obj@meta.data$detailed), decreasing = TRUE), 15))

#


library(ggplot2)
library(dplyr)
library(RColorBrewer)

# ===================== 第一步：数据准备 + 统一添加分面变量 =====================
spatial_data <- obj@meta.data %>%
  select(
    CenterX_global_px, CenterY_global_px,
    detailed, sample, fov, tissue,
    Area.um2, Circularity, AspectRatio, NucArea,
    Mean.PanCK, Mean.CD45, Mean.DAPI
  ) %>%
  filter(!is.na(CenterX_global_px), !is.na(CenterY_global_px)) %>%
  # 1. 筛选目标tissue
  filter(tissue %in% c("RLN_P", "RLN_N", "NRLN", "RT", "NRT")) %>%
  filter(!is.na(tissue)) %>%
  # 2. 标注plot_type（图层顺序：从底层到顶层）
  mutate(
    plot_type = case_when(
      detailed == "GZMK+CXCR6+Tmem" ~ "GZMK+CXCR6+Tmem",  # 最顶层
      detailed %in% c("CD6+DUSP2+Tmem", "LTB+ARHGDIB+Tmem", "TCF7+CCR7+Tmem") |
        grepl("Memory_T|Exhausted_T|Effector_T|Quiescent_T|NKT", detailed) ~ "Other T cells",  # 第四层
      grepl("Malignant cells", detailed, ignore.case = TRUE) ~ "Malignant cells",  # 第三层
      grepl("Dendritic|DC", detailed, ignore.case = TRUE) ~ "Dendritic cells",  # 第二层
      TRUE ~ "Other cells"  # 最底层
    ),
    # 3. 提前创建tissue_sample分面变量（所有数据都包含）
    tissue = factor(tissue, levels = c("RLN_P", "RLN_N", "NRLN", "RT", "NRT")),
    tissue_sample = paste(tissue, sample, sep = "_")  # 分面核心变量
  )

# 定义颜色（GZMK+CXCR6+Tmem改为深蓝色，其他不变）
plot_colors <- c(
  "GZMK+CXCR6+Tmem" = "#1E3A8A",    # 深蓝色（最顶层）
  "Other T cells" = "#B3CDE3",      # 浅蓝色（第四层）
  "Malignant cells" = "#8B3A3A",    # 深红褐色（第三层）
  "Dendritic cells" = "#CCEBC5",    # 浅绿色（第二层）
  "Other cells" = "#DCDCDC"         # 浅灰色（最底层）
)

# 调整图例顺序（从上到下 = 从顶层到底层）
legend_order <- c("GZMK+CXCR6+Tmem", "Other T cells", "Malignant cells", 
                  "Dendritic cells", "Other cells")

# ===================== 第二步：分组1绘图（RLN_P,RLN_N,NRLN → 一张图） =====================
group1_tissues <- c("RLN_P", "RLN_N", "NRLN")
group1_data <- spatial_data %>% filter(tissue %in% group1_tissues)

# 检查分组1数据
if(nrow(group1_data) == 0) {
  stop("分组1（RLN_P/RLN_N/NRLN）无数据！")
}

# 按tissue_sample排序（确保分面顺序正确）
group1_data$tissue_sample <- factor(
  group1_data$tissue_sample,
  levels = unique(group1_data$tissue_sample[order(group1_data$tissue, group1_data$sample)])
)

# 计算分面行数（每行3个sample）
total_facets1 <- length(unique(group1_data$tissue_sample))
n_rows1 <- ceiling(total_facets1 / 3)

# 绘图：图层顺序从底层到顶层
p_group1 <- ggplot(group1_data,
                   aes(x = CenterX_global_px, y = CenterY_global_px, 
                       color = plot_type, size = Area.um2)) +
  # 第1层（最底层）：Other cells
  geom_point(data = filter(group1_data, plot_type == "Other cells"),
             alpha = 0.7) +
  # 第2层：Dendritic cells
  geom_point(data = filter(group1_data, plot_type == "Dendritic cells"),
             alpha = 0.7) +
  # 第3层：Malignant cells
  geom_point(data = filter(group1_data, plot_type == "Malignant cells"),
             alpha = 0.8, size = 1.2) +
  # 第4层：Other T cells
  geom_point(data = filter(group1_data, plot_type == "Other T cells"),
             alpha = 0.7) +
  # 第5层（最顶层）：GZMK+CXCR6+Tmem（深蓝色）
  geom_point(data = filter(group1_data, plot_type == "GZMK+CXCR6+Tmem"),
             alpha = 0.9, size = 1.5) +  # 最亮最大
  # 样式设置
  scale_size_continuous(range = c(0.1, 0.5), guide = "none") +  # 隐藏size图例
  scale_color_manual(values = plot_colors, breaks = legend_order) +
  facet_wrap(~tissue_sample, scales = "free", ncol = 3, nrow = n_rows1) +
  theme_minimal() +
  labs(title = "Cell Distribution - RLN_P / RLN_N / NRLN",
       x = "X Coordinate (um)",
       y = "Y Coordinate (um)",
       color = "Cell Type") +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    strip.text = element_text(size = 10, face = "bold"),
    strip.background = element_rect(fill = "gray90"),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9, face = "bold"),
    axis.text = element_text(size = 6),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_line(color = "gray95")
  )

# 保存分组1PDF
pdf("cell_distribution_RLN_P_RLN_N_NRLN.pdf",
    width = 15, height = 15)
print(p_group1)
dev.off()
print("Saved: cell_distribution_RLN_P_RLN_N_NRLN.pdf")

group2_tissues <- c("RT", "NRT")
group2_data <- spatial_data %>% filter(tissue %in% group2_tissues)

# 检查分组2数据
if(nrow(group2_data) == 0) {
  stop("分组2（RT/NRT）无数据！")
}

# 按tissue_sample排序
group2_data$tissue_sample <- factor(
  group2_data$tissue_sample,
  levels = unique(group2_data$tissue_sample[order(group2_data$tissue, group2_data$sample)])
)

# 计算分面行数
total_facets2 <- length(unique(group2_data$tissue_sample))
n_rows2 <- ceiling(total_facets2 / 3)

# 绘图：图层顺序从底层到顶层
p_group2 <- ggplot(group2_data,
                   aes(x = CenterX_global_px, y = CenterY_global_px, 
                       color = plot_type, size = Area.um2)) +
  # 第1层（最底层）：Other cells
  geom_point(data = filter(group2_data, plot_type == "Other cells"),
             alpha = 0.7) +
  # 第2层：Dendritic cells
  geom_point(data = filter(group2_data, plot_type == "Dendritic cells"),
             alpha = 0.7) +
  # 第3层：Malignant cells
  geom_point(data = filter(group2_data, plot_type == "Malignant cells"),
             alpha = 0.8, size = 1.2) +
  # 第4层：Other T cells
  geom_point(data = filter(group2_data, plot_type == "Other T cells"),
             alpha = 0.7) +
  # 第5层（最顶层）：GZMK+CXCR6+Tmem（深蓝色）
  geom_point(data = filter(group2_data, plot_type == "GZMK+CXCR6+Tmem"),
             alpha = 0.9, size = 1.5) +
  # 样式设置
  scale_size_continuous(range = c(0.1, 0.5), guide = "none") +  # 隐藏size图例
  scale_color_manual(values = plot_colors, breaks = legend_order) +
  facet_wrap(~tissue_sample, scales = "free", ncol = 3, nrow = n_rows2) +
  theme_minimal() +
  labs(title = "Cell Distribution - RT / NRT",
       x = "X Coordinate (um)",
       y = "Y Coordinate (um)",
       color = "Cell Type") +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    strip.text = element_text(size = 10, face = "bold"),
    strip.background = element_rect(fill = "gray90"),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9, face = "bold"),
    axis.text = element_text(size = 6),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_line(color = "gray95")
  )

# 保存分组2PDF
pdf("eff_cell_distribution_RT_NRT.pdf", width = 15, height = 10)
print(p_group2)
dev.off()
print("Saved: cell_distribution_RT_NRT.pdf")


#表达水平
gene_name <- "CXCR6"
gene_exp <- FetchData(obj, vars = gene_name)
spatial_data$gene_exp <- gene_exp[rownames(spatial_data), gene_name]
spatial_data_sorted <- spatial_data %>%
  arrange(gene_exp)   # 从低到高排序

p_gene <- ggplot(spatial_data_sorted,
                 aes(x = CenterX_global_px,
                     y = CenterY_global_px)) +
  geom_point(aes(color = gene_exp),
             size = 1.2,
             alpha = 0.9) +
  scale_color_gradientn(
    colors = c("lightgrey", "yellow", "red"),
    name = paste0(gene_name, " expression")
  ) +
  coord_fixed() +
  theme_minimal()
ggsave("p11.png",plot=p_gene,bg="white")
```
# 统计比例
```R
#分布
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(17)

# 修改后的颜色（替换了第5个橙色、第8个粉色、第9个灰色）
cluster_colors <- c(
  "#76aef3", "#377EB8", "#4DAF4A", "#984EA3", "#1B9E77",  # 第5个：橙色 → 蓝绿色
  "#f2f286", "#A65628", "#D95F02", "#7570B3",            # 第8个：粉色 → 橙色；第9个：灰色 → 紫色
  "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", 
  "#FFD92F", "#E5C494", "#f18686"
)
cluster_colors <- c(
  # 浅色系在前（浅蓝第一个）
  "#A0CBE8",   # 1 浅蓝（你要的放第一个）
  "#FFBE7D",   # 2 浅橙
  "#8CD17D",   # 3 浅绿
  "#86BCB6",   # 4 浅青
  "#FF9D9A",   # 5 浅红
  "#FABFD2",   # 6 浅粉
  "#D4A6C8",   # 7 浅紫
  "#F1CE63",   # 8 淡黄
  "#D7B5A6",   # 9 浅棕
  "#B07AA1",   # 10 紫灰
  "#BAB0AC",   # 11 浅灰
  
  # 深色系在后
  "#4E79A7",   # 12 深蓝
  "#F28E2B",   # 13 橙
  "#59A14F",   # 14 深绿
  "#E15759",   # 15 红
  "#499894",   # 16 青绿
  "#D37295",   # 17 玫红
  "#B6992D",   # 18 土黄
  "#9D7660",   # 19 棕
  "#79706E"    # 20 深灰
)

cell_types <- as.character(unique(obj@meta.data$CellType))
num_legend_items <- length(cell_types)
max_label_length <- max(nchar(cell_types))

cell_types_sorted <- sort(cell_types) 
names(cluster_colors) <- cell_types_sorted
# 将 CellType 列转换为因子，顺序按字母排序
obj@meta.data$CellType <- factor(obj@meta.data$CellType, levels = cell_types_sorted)

base_width <- 3000
base_height <- 3000
legend_width_factor <- 100
label_length_factor <- 10
dynamic_width <- base_width + (num_legend_items * legend_width_factor) + (max_label_length * label_length_factor)

pdf("annotation.pdf", width = dynamic_width/300, height = base_height/300)
DimPlot(obj, reduction = "umap", label = TRUE, pt.size = 1, 
        group.by = "CellType", label.size = 4, alpha = 0.35) +
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    ggtitle(NULL) +
    scale_color_manual(values = cluster_colors) +
    coord_fixed(ratio = 1) +
    guides(color = guide_legend(
        title = NULL, 
        override.aes = list(size = 5, alpha = 1)  # 只改透明度，不改边框
    )) +
    theme(
        text = element_text(size = 8, face = "bold"),
        axis.text.x = element_text(size = 28, color = "black"),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title.x = element_text(size = 36, face = "bold", color = "black"),
        axis.title.y = element_text(size = 36, face = "bold", color = "black", margin = margin(r = 20)),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 16, face = "bold", color = "black"),
        legend.title = element_text(size = 16, face = "bold", color = "black"),
        legend.position = "right",
        legend.box.margin = margin(0, 0, 0, 0),
        legend.key = element_blank(),
        legend.background = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5),
        aspect.ratio = 1,
        plot.margin = margin(10, 50, 10, 10)
    )
dev.off()

library(RColorBrewer)
cluster_colors <- c(brewer.pal(12, "Set3"), brewer.pal(5, "Paired"))
cell_types <- as.character(unique(obj@meta.data$CellType))
cell_types_sorted <- sort(cell_types)  # 按字母顺序排序

# 确保颜色数量与细胞类型数量匹配
if (length(cluster_colors) < length(cell_types_sorted)) {
  cluster_colors <- rep(cluster_colors, length.out = length(cell_types_sorted))
} else {
  cluster_colors <- cluster_colors[1:length(cell_types_sorted)]
}

# 将颜色按排序后的细胞类型命名
names(cluster_colors) <- cell_types_sorted

# 将 CellType 列转换为因子，顺序按字母排序
obj@meta.data$CellType <- factor(obj@meta.data$CellType, levels = cell_types_sorted)

num_legend_items <- length(cell_types_sorted)
max_label_length <- max(nchar(cell_types_sorted))

base_width <- 3000
base_height <- 3000
legend_width_factor <- 100
label_length_factor <- 10
dynamic_width <- base_width + (num_legend_items * legend_width_factor) + (max_label_length * label_length_factor)

pdf("annotation.pdf", width = dynamic_width/300, height = base_height/300)
DimPlot(obj, reduction = "umap", label = TRUE, pt.size = 1, 
        group.by = "CellType", label.size = 4, alpha = 0.35,
        cols = cluster_colors) +  # 使用 cols 参数直接传入颜色
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    ggtitle(NULL) +
    coord_fixed(ratio = 1) +
    guides(color = guide_legend(
        title = NULL, 
        override.aes = list(size = 5, alpha = 1)
    )) +
    theme(
        text = element_text(size = 8, face = "bold"),
        axis.text.x = element_text(size = 28, color = "black"),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title.x = element_text(size = 36, face = "bold", color = "black"),
        axis.title.y = element_text(size = 36, face = "bold", color = "black", margin = margin(r = 20)),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 16, face = "bold", color = "black"),
        legend.title = element_text(size = 16, face = "bold", color = "black"),
        legend.position = "right",
        legend.box.margin = margin(0, 0, 0, 0),
        legend.key = element_blank(),
        legend.background = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5),
        aspect.ratio = 1,
        plot.margin = margin(10, 50, 10, 10)
    )
dev.off()


library(tidyverse)

# 计算各亚群在不同组织中的占比
prop_data <- obj@meta.data %>%
  group_by(tissue, CellType) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(tissue) %>%
  mutate(proportion = count / sum(count) * 100) %>%
  filter(tissue %in% c("mPT", "nmPT")) %>%
  mutate(tissue = factor(tissue, levels = c("mPT", "nmPT")))

library(ggplot2)
library(ggsci)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(17)


prop_data <- prop_data %>%
  mutate(CellType = as.character(CellType))

# 绘制堆叠条形图
p <- ggplot(prop_data, aes(x = tissue, y = proportion, fill = CellType)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = cluster_colors) +
  labs(title = "Cell Type Distribution: mPT vs nmPT",
       x = "Tissue", y = "Proportion (%)", fill = "Cell Type") +
  theme_minimal() +
  theme(
    # 坐标轴线条和刻度
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.15, "cm"),
    # 坐标轴文本
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 10),
    axis.text.y = element_text(color = "black", size = 10),
    # 坐标轴标题
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    # 面板背景
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    # 图例
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )


# 保存为PDF
pdf("obj_celltype_distribution_tumor.pdf", width = 4, height = 6)
print(p)
dev.off()

filtered_meta <- obj@meta.data %>%
  filter(!(tissue == "negLN" & CellType %in% c("Malignant cells", "Alveolar type II cells")))
prop_data <- filtered_meta %>%
  group_by(tissue, CellType) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(tissue) %>%
  mutate(proportion = count / sum(count) * 100) %>%
  filter(tissue %in% c("metLN", "negLN")) %>%
  mutate(tissue = factor(tissue, levels = c("metLN", "negLN")))

library(ggplot2)
library(ggsci)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(17)


prop_data <- prop_data %>%
  mutate(CellType = as.character(CellType))

# 绘制堆叠条形图
p <- ggplot(prop_data, aes(x = tissue, y = proportion, fill = CellType)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = cluster_colors) +
  labs(title = "Cell Type Distribution: metLN vs negLN",
       x = "Tissue", y = "Proportion (%)", fill = "Cell Type") +
  theme_minimal() +
  theme(
    # 坐标轴线条和刻度
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.15, "cm"),
    # 坐标轴文本
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 10),
    axis.text.y = element_text(color = "black", size = 10),
    # 坐标轴标题
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    # 面板背景
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    # 图例
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )


# 保存为PDF
pdf("obj_celltype_distribution_RLN_P_RLN_N_NRLN.pdf", width = 4, height = 6)
print(p)
dev.off()

library(dplyr)
library(ggplot2)
library(patchwork)

# 过滤 negLN 中的恶性细胞和肺泡Ⅱ型细胞
filtered_meta <- obj@meta.data %>%
  filter(!(tissue == "negLN" & CellType %in% c("Malignant cells", "Alveolar type II cells")))

# 计算每个样本中每种细胞类型的比例
prop_by_sample <- filtered_meta %>%
  group_by(tissue, sample, CellType) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(tissue, sample) %>%
  mutate(proportion = count / sum(count) * 100) %>%
  ungroup()

# 获取所有细胞类型（按总体中位数排序）
all_cell_types <- prop_by_sample %>%
  group_by(CellType) %>%
  summarise(median_prop = median(proportion), .groups = "drop") %>%
  arrange(desc(median_prop)) %>%
  pull(CellType)

# 组织顺序和颜色
tissue_order <- c("mPT", "nmPT", "metLN", "negLN")
tissue_colors <- c("mPT" = "#E41A1C", "nmPT" = "#377EB8", 
                   "metLN" = "#4DAF4A", "negLN" = "#984EA3")

# 为每个细胞类型单独绘图
plot_list <- list()

for (i in seq_along(all_cell_types)) {
  cell_name <- all_cell_types[i]
  
  # 筛选当前细胞类型的数据
  cell_data <- prop_by_sample %>% filter(CellType == cell_name)
  
  # 确保组织顺序
  cell_data$tissue <- factor(cell_data$tissue, levels = tissue_order)
  
  p <- ggplot(cell_data, aes(x = tissue, y = proportion)) +
    # 箱线图：透明/白色，只有边框
    geom_boxplot(alpha = 0, fill = "white", color = "gray50", 
                 outlier.shape = NA, width = 0.6) +
    # 点：按组织着色，半透明
    geom_jitter(aes(color = tissue), width = 0.2, size = 1.5, alpha = 0.6) +
    scale_color_manual(values = tissue_colors) +
    labs(title = cell_name,
         x = NULL,
         y = "Proportion (%)") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y = element_text(size = 8),
      axis.title.y = element_text(size = 9, face = "bold"),
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
      legend.position = "none",  # 图例统一加，这里去掉
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_blank()
    )
  
  plot_list[[i]] <- p
}

# 拼图：每行4个
n_plots <- length(plot_list)
ncol <- 4
nrow <- ceiling(n_plots / ncol)

combined_plot <- wrap_plots(plot_list, ncol = ncol, nrow = nrow)

# 添加标题
combined_plot <- combined_plot + 
  plot_annotation(title = "Cell Type Proportion Across Tissues",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14)))

# 保存
ggsave("celltype_boxplot_transparent.pdf", combined_plot, 
       width = 14, height = 3.5 * nrow, dpi = 300, limitsize = FALSE)

library(dplyr)
library(ggplot2)
library(patchwork)

# 你指定的6种细胞类型
target_cell_types <- c("Malignant cells", "Fibroblasts", "Endothelial cells", 
                       "T cells", "Macrophages", "Dendritic cells")

# 过滤 negLN 中的恶性细胞和肺泡Ⅱ型细胞（如果需要）
filtered_meta <- obj@meta.data %>%
  filter(!(tissue == "negLN" & CellType %in% c("Malignant cells", "Alveolar type II cells")))

# 计算每个样本中每种细胞类型的比例
prop_by_sample <- filtered_meta %>%
  group_by(tissue, sample, CellType) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(tissue, sample) %>%
  mutate(proportion = count / sum(count) * 100) %>%
  ungroup()

# 筛选目标细胞类型
prop_filtered <- prop_by_sample %>%
  filter(CellType %in% target_cell_types)

# 确保 CellType 顺序（按你指定的顺序）
prop_filtered$CellType <- factor(prop_filtered$CellType, levels = target_cell_types)

# 组织顺序
tissue_order <- c("mPT", "nmPT", "metLN", "negLN")
tissue_colors <- c("mPT" = "#E41A1C", "nmPT" = "#377EB8", 
                   "metLN" = "#4DAF4A", "negLN" = "#984EA3")

# 为每个细胞类型单独绘图
plot_list <- list()

for (i in seq_along(target_cell_types)) {
  cell_name <- target_cell_types[i]
  
  # 筛选当前细胞类型的数据
  cell_data <- prop_filtered %>% filter(CellType == cell_name)
  
  # 确保组织顺序
  cell_data$tissue <- factor(cell_data$tissue, levels = tissue_order)
  
  # 简化细胞类型名称用于标题
  title_name <- case_when(
    cell_name == "Malignant cells" ~ "Malignant",
    cell_name == "Fibroblasts" ~ "Fibroblast",
    cell_name == "Endothelial cells" ~ "Endothelial",
    cell_name == "T cells" ~ "T cell",
    cell_name == "Macrophages" ~ "Macrophage",
    cell_name == "Dendritic cells" ~ "Dendritic cell",
    TRUE ~ cell_name
  )
  
  p <- ggplot(cell_data, aes(x = tissue, y = proportion)) +
    # 箱线图：透明，灰色边框（与密度图一致）
    geom_boxplot(alpha = 0, fill = "white", color = "gray50", 
                 outlier.shape = NA, width = 0.6) +
    # 点：按组织着色，半透明
    geom_jitter(aes(color = tissue), width = 0.2, size = 1.5, alpha = 0.6) +
    scale_color_manual(values = tissue_colors) +
    labs(title = title_name,
         x = NULL,
         y = NULL) +  # 去掉 y 轴标签，拼图时统一加
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 7),
      plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
      legend.position = "none",
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_blank()
    )
  
  plot_list[[i]] <- p
}

# 拼图：一行3个，共2行
combined_plot <- wrap_plots(plot_list, ncol = 3, nrow = 2)

# 添加统一的 y 轴标签
combined_plot <- combined_plot + 
  plot_annotation(
    title = "Cell Type Proportion Across Tissues",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
  ) &
  ylab("Proportion (%)")

# 提取一个图例
legend_data <- data.frame(tissue = names(tissue_colors))
legend_p <- ggplot(legend_data, aes(x = 1, y = 1, color = tissue)) +
  geom_point(size = 3) +
  scale_color_manual(values = tissue_colors, name = "Tissue") +
  theme_void() +
  theme(legend.position = "bottom", legend.title = element_text(size = 10))

legend <- cowplot::get_legend(legend_p)

# 将图例添加到拼图底部
final_plot <- combined_plot + 
  patchwork::inset_element(legend, 
                           left = 0.3, bottom = 0, right = 0.7, top = 0.05,
                           align_to = "full")

# 调整底部边距
final_plot <- final_plot & theme(plot.margin = margin(5, 5, 30, 5))

# 保存
ggsave("celltype_proportion_boxplot_6types.pdf", final_plot, 
       width = 10, height = 6, dpi = 300)

#显著性
library(dplyr)
library(ggplot2)
library(patchwork)
library(rstatix)
library(ggpubr)
library(cowplot)

target_cell_types <- c("Malignant cells", "Fibroblasts", "Endothelial cells", 
                       "T cells", "Macrophages", "Dendritic cells")

filtered_meta <- obj@meta.data %>%
  filter(!(tissue == "negLN" & CellType %in% c("Malignant cells", "Alveolar type II cells")))

prop_by_sample <- filtered_meta %>%
  group_by(tissue, sample, CellType) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(tissue, sample) %>%
  mutate(proportion = count / sum(count) * 100) %>%
  ungroup()

prop_filtered <- prop_by_sample %>%
  filter(CellType %in% target_cell_types)

prop_filtered$CellType <- factor(prop_filtered$CellType, levels = target_cell_types)

tissue_order <- c("mPT", "nmPT", "metLN", "negLN")
tissue_colors <- c("mPT" = "#E41A1C", "nmPT" = "#377EB8", 
                   "metLN" = "#4DAF4A", "negLN" = "#984EA3")

# ==========================================
# 🔥 小样本专用统计（最推荐）
# ==========================================
# 整体 KW 检验：BH 校正（应付审稿人）
kw_results <- prop_filtered %>%
  group_by(CellType) %>%
  kruskal_test(proportion ~ tissue) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")

# 两两 Dunn：原始 p，不校正（小样本必须这样）
dunn_results <- prop_filtered %>%
  group_by(CellType) %>%
  dunn_test(proportion ~ tissue, p.adjust.method = "none") %>%  # 不校正
  add_significance("p")                                         # 用原始 p

sig_dunn <- dunn_results %>%
  filter(p < 0.05) %>%
  mutate(group1 = as.character(group1),
         group2 = as.character(group2))

# ==========================================
# 绘图：和你密度图完全一样的参数
# ==========================================
plot_list <- list()

for (i in seq_along(target_cell_types)) {
  cell_name <- target_cell_types[i]
  
  cell_data <- prop_filtered %>% filter(CellType == cell_name)
  cell_data$tissue <- factor(cell_data$tissue, levels = tissue_order)
  ct_sig <- sig_dunn %>% filter(CellType == cell_name)
  
  title_name <- case_when(
    cell_name == "Malignant cells" ~ "Malignant",
    cell_name == "Fibroblasts" ~ "Fibroblast",
    cell_name == "Endothelial cells" ~ "Endothelial",
    cell_name == "T cells" ~ "T cell",
    cell_name == "Macrophages" ~ "Macrophage",
    cell_name == "Dendritic cells" ~ "Dendritic cell",
    TRUE ~ cell_name
  )
  
  p <- ggplot(cell_data, aes(x = tissue, y = proportion, color = tissue)) +
    geom_boxplot(alpha = 0, fill = "white", color = "gray50", 
                 outlier.shape = NA, width = 0.6) +
    geom_jitter(aes(color = tissue), width = 0.2, size = 1.5, alpha = 0.6) +
    scale_color_manual(values = tissue_colors) +
    scale_y_continuous(labels = scales::comma_format(accuracy = 1)) +
    labs(title = title_name, x = NULL, y = "Proportion (%)") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 7),
      plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
      legend.position = "none",
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_blank()
    )
  
  if (nrow(ct_sig) > 0) {
    comps <- lapply(1:nrow(ct_sig), function(j) {
      c(ct_sig$group1[j], ct_sig$group2[j])
    })
    annots <- ct_sig$p.signif
    
    p <- p + geom_signif(
      comparisons = comps,
      annotations = annots,
      map_signif_level = FALSE,
      tip_length = 0.008,
      size = 0.3,
      textsize = 2.5,
      vjust = 0.4,
      step_increase = 0.10,    # 和密度图一样
      margin_text = 0.08
    )
  }
  
  plot_list[[i]] <- p
}

combined_plot <- wrap_plots(plot_list, ncol = 3, nrow = 2, axes = "keep")

combined_plot <- combined_plot + 
  plot_annotation(
    title = "Cell Type Proportion Across Tissues",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12))
  )

legend_data <- data.frame(tissue = names(tissue_colors))
legend_p <- ggplot(legend_data, aes(x = 1, y = 1, color = tissue)) +
  geom_point(size = 3) +
  scale_color_manual(values = tissue_colors, name = "Tissue") +
  theme_void() +
  theme(legend.position = "bottom", legend.title = element_text(size = 10))

legend <- get_legend(legend_p)

final_plot <- combined_plot + 
  inset_element(legend, left = 0.3, bottom = 0, right = 0.7, top = 0.05, align_to = "full")

final_plot <- final_plot & theme(plot.margin = margin(5, 5, 30, 5))

ggsave("celltype_proportion_boxplot_small_sample.pdf", final_plot, width = 12, height = 7, dpi = 300)

```
# 密度
```R
# 方法1：统计每个样本的 FOV 数量
fov_count_by_sample <- obj@meta.data %>%
  group_by(sample, tissue) %>%
  summarise(n_fov = n_distinct(fov), .groups = "drop")

# 查看结果
print(fov_count_by_sample)

# 方法2：按组织类型汇总
fov_summary <- obj@meta.data %>%
  group_by(tissue, sample) %>%
  summarise(n_fov = n_distinct(fov), .groups = "drop") %>%
  group_by(tissue) %>%
  summarise(
    total_fov = sum(n_fov),
    mean_fov_per_sample = mean(n_fov),
    min_fov = min(n_fov),
    max_fov = max(n_fov),
    .groups = "drop"
  )

print(fov_summary)

library(dplyr)
library(ggplot2)
library(patchwork)

# 过滤 negLN 中的恶性细胞和肺泡Ⅱ型细胞
filtered_meta <- obj@meta.data %>%
  filter(!(tissue == "negLN" & CellType %in% c("Malignant cells", "Alveolar type II cells")))

# 计算每个样本的总组织面积（所有细胞面积之和）
tissue_area <- filtered_meta %>%
  group_by(sample, tissue) %>%
  summarise(total_area_mm2 = sum(Area.um2, na.rm = TRUE) / 1e6, .groups = "drop")

# 计算每个样本中每种细胞类型的数量和密度
density_data <- filtered_meta %>%
  group_by(sample, tissue, CellType) %>%
  summarise(count = n(), .groups = "drop") %>%
  left_join(tissue_area, by = c("sample", "tissue")) %>%
  mutate(density = count / total_area_mm2)  # cells/mm²

# 获取所有细胞类型（按中位数密度排序）
celltype_order <- density_data %>%
  group_by(CellType) %>%
  summarise(median_density = median(density, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(median_density)) %>%
  pull(CellType)

# 组织顺序
tissue_order <- c("mPT", "nmPT", "metLN", "negLN")

# 颜色方案（按组织）
tissue_colors <- c("mPT" = "#E41A1C", "nmPT" = "#377EB8", 
                   "metLN" = "#4DAF4A", "negLN" = "#984EA3")

# 为每个细胞类型单独绘图
plot_list <- list()

for (i in seq_along(celltype_order)) {
  ct <- celltype_order[i]
  
  # 筛选当前细胞类型的数据
  ct_data <- density_data %>% filter(CellType == ct)
  
  # 确保组织顺序
  ct_data$tissue <- factor(ct_data$tissue, levels = tissue_order)
  
  # 简化细胞类型名称用于标题（可选）
  title_name <- ct
  
  p <- ggplot(ct_data, aes(x = tissue, y = density)) +
    # 箱线图：透明，灰色边框
    geom_boxplot(alpha = 0, fill = "white", color = "gray50", 
                 outlier.shape = NA, width = 0.6) +
    # 点：按组织着色，半透明
    geom_jitter(aes(color = tissue), width = 0.2, size = 1.5, alpha = 0.6) +
    scale_color_manual(values = tissue_colors) +
    scale_y_log10() +
    labs(title = title_name,
         x = NULL,
         y = NULL) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 7),
      plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
      legend.position = "none",
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_blank()
    )
  
  plot_list[[i]] <- p
}

# 拼图：每行4个细胞类型
ncol <- 4
nrow <- ceiling(length(plot_list) / ncol)

# 使用 wrap_plots 拼图
combined_plot <- wrap_plots(plot_list, ncol = ncol, nrow = nrow)

# 添加统一的 y 轴标签
combined_plot <- combined_plot + 
  plot_annotation(
    title = "Cell Type Density Across Tissues",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
  ) &
  ylab(expression("Density (cells/mm"^2~")"))

# 提取一个图例（用于统一图例）
legend_data <- data.frame(tissue = names(tissue_colors))
legend_p <- ggplot(legend_data, aes(x = 1, y = 1, color = tissue)) +
  geom_point(size = 3) +
  scale_color_manual(values = tissue_colors, name = "Tissue") +
  theme_void() +
  theme(legend.position = "bottom", legend.title = element_text(size = 10))

legend <- cowplot::get_legend(legend_p)

# 将图例添加到拼图底部
final_plot <- combined_plot + 
  patchwork::inset_element(legend, 
                           left = 0.3, bottom = 0, right = 0.7, top = 0.05,
                           align_to = "full")

# 调整底部边距，为图例留空间
final_plot <- final_plot & theme(plot.margin = margin(5, 5, 30, 5))

# 保存
ggsave("celltype_density_boxplot_transparent.pdf", final_plot, 
       width = 14, height = 3 * nrow, dpi = 300, limitsize = FALSE)

#每个fov的密度
library(dplyr)
library(ggplot2)
library(patchwork)

# 查看 metadata 中是否有 fov 信息
# 从你的 head(obj@meta.data) 输出中看到有 "fov" 列
# 如果没有，可能需要从其他列提取

# 过滤 negLN 中的恶性细胞和肺泡Ⅱ型细胞
filtered_meta <- obj@meta.data %>%
  filter(!(tissue == "negLN" & CellType %in% c("Malignant cells", "Alveolar type II cells")))

# 计算每个 FOV 的总组织面积（所有细胞面积之和）
# 注意：这里按 fov 分组，而不是 sample
fov_area <- filtered_meta %>%
  group_by(fov, tissue) %>%  # 按 fov 和 tissue 分组
  summarise(total_area_mm2 = sum(Area.um2, na.rm = TRUE) / 1e6, .groups = "drop")

# 计算每个 FOV 中每种细胞类型的数量和密度
density_data <- filtered_meta %>%
  group_by(fov, tissue, CellType) %>%
  summarise(count = n(), .groups = "drop") %>%
  left_join(fov_area, by = c("fov", "tissue")) %>%
  mutate(density = count / total_area_mm2)  # cells/mm²

# 获取所有细胞类型（按中位数密度排序）
celltype_order <- density_data %>%
  group_by(CellType) %>%
  summarise(median_density = median(density, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(median_density)) %>%
  pull(CellType)

# 组织顺序
tissue_order <- c("mPT", "nmPT", "metLN", "negLN")

# 颜色方案（按组织）
tissue_colors <- c("mPT" = "#E41A1C", "nmPT" = "#377EB8", 
                   "metLN" = "#4DAF4A", "negLN" = "#984EA3")

# 为每个细胞类型单独绘图
plot_list <- list()

for (i in seq_along(celltype_order)) {
  ct <- celltype_order[i]
  
  # 筛选当前细胞类型的数据
  ct_data <- density_data %>% filter(CellType == ct)
  
  # 确保组织顺序
  ct_data$tissue <- factor(ct_data$tissue, levels = tissue_order)
  
  # 简化细胞类型名称用于标题
  title_name <- ct
  
  p <- ggplot(ct_data, aes(x = tissue, y = density)) +
    # 箱线图：透明，灰色边框
    geom_boxplot(alpha = 0, fill = "white", color = "gray50", 
                 outlier.shape = NA, width = 0.6) +
    # 点：按组织着色，每个点代表一个 FOV
    geom_jitter(aes(color = tissue), width = 0.2, size = 1.5, alpha = 0.6) +
    scale_color_manual(values = tissue_colors) +
    scale_y_log10() +
    labs(title = title_name,
         x = NULL,
         y = NULL) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 7),
      plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
      legend.position = "none",
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_blank()
    )
  
  plot_list[[i]] <- p
}

# 拼图：每行4个细胞类型
ncol <- 4
nrow <- ceiling(length(plot_list) / ncol)

# 使用 wrap_plots 拼图
combined_plot <- wrap_plots(plot_list, ncol = ncol, nrow = nrow)

# 添加统一的 y 轴标签
combined_plot <- combined_plot + 
  plot_annotation(
    title = "Cell Type Density Across Tissues (per FOV)",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
  ) &
  ylab(expression("Density (cells/mm"^2~")"))

# 提取一个图例
legend_data <- data.frame(tissue = names(tissue_colors))
legend_p <- ggplot(legend_data, aes(x = 1, y = 1, color = tissue)) +
  geom_point(size = 3) +
  scale_color_manual(values = tissue_colors, name = "Tissue") +
  theme_void() +
  theme(legend.position = "bottom", legend.title = element_text(size = 10))

legend <- cowplot::get_legend(legend_p)

# 将图例添加到拼图底部
final_plot <- combined_plot + 
  patchwork::inset_element(legend, 
                           left = 0.3, bottom = 0, right = 0.7, top = 0.05,
                           align_to = "full")

# 调整底部边距
final_plot <- final_plot & theme(plot.margin = margin(5, 5, 30, 5))

# 保存
ggsave("celltype_density_boxplot_per_fov.pdf", final_plot, 
       width = 14, height = 3 * nrow, dpi = 300, limitsize = FALSE)

for (i in seq_along(target_cell_types)) {
  ct <- target_cell_types[i]
  
  # 筛选当前细胞类型的数据
  ct_data <- density_data %>% filter(CellType == ct)
  
  # 确保组织顺序
  ct_data$tissue <- factor(ct_data$tissue, levels = tissue_order)
  
  # 简化细胞类型名称用于标题
  title_name <- case_when(
    ct == "Malignant cells" ~ "Malignant",
    ct == "Fibroblasts" ~ "Fibroblast",
    ct == "Endothelial cells" ~ "Endothelial",
    ct == "T cells" ~ "T cell",
    ct == "Macrophages" ~ "Macrophage",
    ct == "Dendritic cells" ~ "Dendritic cell",
    TRUE ~ ct
  )
  
  p <- ggplot(ct_data, aes(x = tissue, y = density)) +
    # 箱线图：透明，灰色边框
    geom_boxplot(alpha = 0, fill = "white", color = "gray50", 
                 outlier.shape = NA, width = 0.6) +
    # 点：按组织着色，每个点代表一个 FOV
    geom_jitter(aes(color = tissue), width = 0.2, size = 1.5, alpha = 0.6) +
    scale_color_manual(values = tissue_colors) +
    scale_y_log10() +
    labs(title = title_name,
         x = NULL,
         y = NULL) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 7),
      plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
      legend.position = "none",
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_blank()
    )
  
  plot_list[[i]] <- p
}

# 拼图：一行3个，共2行（因为只有6个）
ncol <- 3
nrow <- 2

# 使用 wrap_plots 拼图
combined_plot <- wrap_plots(plot_list, ncol = ncol, nrow = nrow)

# 添加统一的 y 轴标签
combined_plot <- combined_plot + 
  plot_annotation(
    title = "Cell Type Density Across Tissues (per FOV)",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
  ) &
  ylab(expression("Density (cells/mm"^2~")"))

# 提取一个图例
legend_data <- data.frame(tissue = names(tissue_colors))
legend_p <- ggplot(legend_data, aes(x = 1, y = 1, color = tissue)) +
  geom_point(size = 3) +
  scale_color_manual(values = tissue_colors, name = "Tissue") +
  theme_void() +
  theme(legend.position = "bottom", legend.title = element_text(size = 10))

legend <- cowplot::get_legend(legend_p)

# 将图例添加到拼图底部
final_plot <- combined_plot + 
  patchwork::inset_element(legend, 
                           left = 0.3, bottom = 0, right = 0.7, top = 0.05,
                           align_to = "full")

# 调整底部边距
final_plot <- final_plot & theme(plot.margin = margin(5, 5, 30, 5))

# 保存
ggsave("celltype_density_boxplot_6types_per_fov.pdf", final_plot, 
       width = 10, height = 6, dpi = 300)

#显著性
library(dplyr)
library(ggplot2)
library(patchwork)
library(rstatix)
library(ggpubr)

# 指定的6种细胞类型
target_cell_types <- c("Malignant cells", "Fibroblasts", "Endothelial cells", 
                       "T cells", "Macrophages", "Dendritic cells")

# 过滤 negLN 中的恶性细胞和肺泡Ⅱ型细胞
filtered_meta <- obj@meta.data %>%
  filter(!(tissue == "negLN" & CellType %in% c("Malignant cells", "Alveolar type II cells")))

# 计算每个 FOV 的总组织面积
fov_area <- filtered_meta %>%
  group_by(fov, tissue) %>%
  summarise(total_area_mm2 = sum(Area.um2, na.rm = TRUE) / 1e6, .groups = "drop")

# 计算每个 FOV 中每种细胞类型的密度
density_data <- filtered_meta %>%
  group_by(fov, tissue, CellType) %>%
  summarise(count = n(), .groups = "drop") %>%
  left_join(fov_area, by = c("fov", "tissue")) %>%
  mutate(density = count / total_area_mm2) %>%
  filter(CellType %in% target_cell_types)

# 组织顺序
tissue_order <- c("mPT", "nmPT", "metLN", "negLN")
density_data$tissue <- factor(density_data$tissue, levels = tissue_order)

# 颜色方案
tissue_colors <- c("mPT" = "#E41A1C", "nmPT" = "#377EB8", 
                   "metLN" = "#4DAF4A", "negLN" = "#984EA3")

# ============================================
# 1. 统计检验
# ============================================

density_data$log_density <- log10(density_data$density + 0.01)

kw_results <- density_data %>%
  group_by(CellType) %>%
  kruskal_test(log_density ~ tissue) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")

dunn_results <- density_data %>%
  group_by(CellType) %>%
  dunn_test(log_density ~ tissue, p.adjust.method = "BH") %>%
  add_significance("p.adj")

sig_dunn <- dunn_results %>%
  filter(p.adj < 0.05) %>%
  mutate(
    group1 = as.character(group1),
    group2 = as.character(group2)
  )

# ============================================
# 2. 可视化
# ============================================

plot_list <- list()

for (i in seq_along(target_cell_types)) {
  ct <- target_cell_types[i]
  
  ct_data <- density_data %>% filter(CellType == ct)
  ct_sig <- sig_dunn %>% filter(CellType == ct)
  
  title_name <- case_when(
    ct == "Malignant cells" ~ "Malignant",
    ct == "Fibroblasts" ~ "Fibroblast",
    ct == "Endothelial cells" ~ "Endothelial",
    ct == "T cells" ~ "T cell",
    ct == "Macrophages" ~ "Macrophage",
    ct == "Dendritic cells" ~ "Dendritic cell",
    TRUE ~ ct
  )
  
  p <- ggplot(ct_data, aes(x = tissue, y = density, color = tissue)) +
    geom_boxplot(alpha = 0, fill = "white", color = "gray50", 
                 outlier.shape = NA, width = 0.6) +
    geom_jitter(width = 0.2, size = 1.5, alpha = 0.6) +
    scale_color_manual(values = tissue_colors) +
    scale_y_log10(labels = scales::comma_format(accuracy = 1)) + 
    labs(title = title_name, x = NULL, y = "Density (cells/mm²)") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 7),
      plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
      legend.position = "none",
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_blank()
    )
  
  if (nrow(ct_sig) > 0) {
    comps <- lapply(1:nrow(ct_sig), function(j) {
      as.character(c(ct_sig$group1[j], ct_sig$group2[j]))
    })
    annots <- ct_sig$p.adj.signif
    
    p <- p + geom_signif(
      comparisons = comps,
      annotations = annots,
      map_signif_level = FALSE,
      tip_length = 0.008,
      size = 0.3,
      textsize = 2.5,
      vjust = 0.4,
      step_increase = 0.05,  # 👈 行间距已调小！
      margin_text = 0.08
    )
  }
  plot_list[[i]] <- p
}

# 拼图：保留所有 y 轴
combined_plot <- wrap_plots(plot_list, ncol = 3, nrow = 2, axes = "keep")

# 添加标题
combined_plot <- combined_plot + 
  plot_annotation(
    title = "Cell Type Density Across Tissues (Kruskal-Wallis + Dunn's test, FDR < 0.05)",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12))
  )

# 图例
legend_data <- data.frame(tissue = names(tissue_colors))
legend_p <- ggplot(legend_data, aes(x = 1, y = 1, color = tissue)) +
  geom_point(size = 3) +
  scale_color_manual(values = tissue_colors, name = "Tissue") +
  theme_void() +
  theme(legend.position = "bottom", legend.title = element_text(size = 10))

legend <- cowplot::get_legend(legend_p)

final_plot <- combined_plot + 
  patchwork::inset_element(legend, 
                           left = 0.3, bottom = 0, right = 0.7, top = 0.05,
                           align_to = "full")

final_plot <- final_plot & theme(plot.margin = margin(5, 5, 30, 5))

# 保存
ggsave("celltype_density_boxplot_KW_Dunn.pdf", final_plot, 
       width = 12, height = 7, dpi = 300)

cat("\n========== Done! ==========\n")


```
# malignant
```R
Malignant <- subset(obj,subset=CellType=="Malignant cells")
Malignant <- subset(Malignant,subset=seurat_clusters %in% c(0,1,2,3,5,7,9,10))
Malignant <- subset(Malignant,subset=tissue %in% c("metLN","mPT","nmPT"))
Malignant <- NormalizeData(Malignant)
Malignant <- FindVariableFeatures(Malignant, nfeatures = 3000)
hvgs <- VariableFeatures(Malignant)
Malignant <- ScaleData(Malignant, features = hvgs)
Malignant <- RunPCA(Malignant, features = hvgs, npcs = 50)
#p <- ElbowPlot(Malignant, ndims = 30)
#ggsave("p3.png",plot=p)
pca_var <- Malignant[["pca"]]@stdev^2 / sum(Malignant[["pca"]]@stdev^2)
cumsum_var <- cumsum(pca_var)

# 看达到多少 PC 时累计方差 > 70% / 80% / 90%
which(cumsum_var > 0.7)[1]   # 70%
which(cumsum_var > 0.8)[1]   # 80%
which(cumsum_var > 0.9)[1]   # 90%
#library(harmony)
#Malignant <- RunHarmony(Malignant,"tissue")

Malignant <- FindNeighbors(Malignant, dims = 1:20)
Malignant <- FindClusters(Malignant, resolution = 0.4)
Malignant <- RunUMAP(Malignant, dims = 1:20)

p <- DimPlot(Malignant, reduction = "umap", group.by = "sample", label = FALSE)
ggsave("malignant_orig.png",plot=p)

library(ggplot2)
library(ggsci)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(15)#0.6

seurat_clusters <- as.character(unique(Malignant@meta.data$seurat_clusters))  # 转换为字符向量
num_legend_items <- length(seurat_clusters)  # 图例的个数
max_label_length <- max(nchar(seurat_clusters))  # 图例名称的最大长度

base_width <- 3000  # 基础宽度
base_height <- 3000  # 基础高度
legend_width_factor <- 100  # 每个图例项增加的宽度
label_length_factor <- 10  # 每个字符增加的宽度
dynamic_width <- base_width + (num_legend_items * legend_width_factor) + (max_label_length * label_length_factor)
pdf("Malignant_clusters.pdf", width = dynamic_width/300, height = base_height/300)  # 转换为英寸

DimPlot(Malignant, reduction = "umap", label = TRUE, pt.size = 1, label.size = 8) +
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    ggtitle(NULL) +
    scale_color_manual(values = npg_extended) +
    coord_fixed(ratio = 1) +
    guides(color = guide_legend(title = NULL, override.aes = list(size = 5))) +
    theme(
        text = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 28, color = "black"),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title.x = element_text(size = 36, face = "bold", color = "black"),
        axis.title.y = element_text(size = 36, face = "bold", color = "black", margin = margin(r = 20)),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 28, face = "bold", color = "black"),
        legend.title = element_text(size = 28, face = "bold", color = "black"),
        legend.position = "right",
        legend.box.margin = margin(0, 0, 0, 0),
        legend.key = element_blank(),
        legend.background = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5),
        aspect.ratio = 1,
        plot.margin = margin(10, 50, 10, 10)
    )

dev.off()

library(dplyr)
Malignant_markers <- FindAllMarkers(Malignant, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "wilcox")
Malignant_significant_markers <- subset(Malignant_markers, p_val_adj < 0.05)
#write.csv(Malignant_significant_markers, "Malignant_all_marker.csv")
Malignant_significant_markers <- Malignant_significant_markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
write.csv(Malignant_significant_markers, "Malignant_top_marker_sct.csv")
#无 unknown
identity_mapping <- c(
  "0" = "Inflamed_1",
  "1" = "Stem-like",
  "2" = "Metabolic_1",
  "3" = "Metabolic_2",
  "4" = "Inflamed_2",
  "5" = "Inflamed_3",  
  "6" = "Inflamed_4",
  "7" = "Proliferative_1",
  "8" = "Club cell",
  "9" = "Proliferative_2",
  "10" = "Differentiated",
  "11" = "EMT-like"
)
#remove
identity_mapping <- c(
  "0" = "Metabolic_1",
  "1" = "Inflamed_1",
  "2" = "Stem-like",
  "3" = "Inflamed_2",
  "4" = "Proliferative_1",
  "5" = "Metabolic_2",  
  "6" = "Proliferative_2",
  "7" = "Differentiated"
)
#0.45
identity_mapping <- c(
  "0" = "Tumor_Metabolic",
  "1" = "Tumor_Inflamed_1",
  "2" = "Tumor_Inflamed_2",
  "3" = "Tumor_Inflamed_3",
  "4" = "Tumor_Inflamed_4",
  "5" = "Tumor_Stem-like",  
  "6" = "Tumor_Proliferative",
  "7" = "Tumor_Inflamed_5",
  "8" = "Tumor_Club-cell",
  "9" = "Tumor_Differentiated"
)
#0.6
identity_mapping <- c(
  "0" = "Tumor_Inflamed_1",
  "1" = "Tumor_Inflamed_2",
  "2" = "Tumor_Squamous_1",
  "3" = "Tumor_Metabolic_1",
  "4" = "Tumor_Inflamed_3",
  "5" = "Tumor_Stem-like",
  "6" = "Tumor_Metabolic_2",
  "7" = "Tumor_Squamous_2",
  "8" = "Tumor_Proliferative_1",
  "9" = "Tumor_Inflamed_4",
  "10" = "Tumor_Secretory",
  "11" = "Tumor_Proliferative_2",
  "12" = "Tumor_Squamous_3",
  "13" = "Tumor_EMT-like"
)

#0.35
identity_mapping <- c(
  "0" = "Tumor_Metabolic_1",
  "1" = "Tumor_Immune-inflamed_1",
  "2" = "Tumor_Immunogenic",
  "3" = "Tumor_Immune-inflamed_2",
  "4" = "Tumor_Stem-like",
  "5" = "Tumor_Immunogenic",
  "6" = "Tumor_Secretory",
  "7" = "Tumor_Proliferative",
  "8" = "Tumor_EMT-like",
  "9" = "Tumor_Metabolic_2"
)
sub_cell_type <- identity_mapping[Malignant@meta.data$seurat_clusters]
Malignant@meta.data$sub_cell_type <- sub_cell_type

library(ggplot2)
library(ggsci)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(14)
cluster_colors <- c(
  # 浅色系在前（浅蓝第一个）
  "#A0CBE8",   # 1 浅蓝（你要的放第一个）
  "#FFBE7D",   # 2 浅橙
  "#8CD17D",   # 3 浅绿
  "#86BCB6",   # 4 浅青
  "#FF9D9A",   # 5 浅红
  "#FABFD2",   # 6 浅粉
  "#D4A6C8",   # 7 浅紫
  "#F1CE63",   # 8 淡黄
  "#D7B5A6",   # 9 浅棕
  "#B07AA1",   # 10 紫灰
  "#BAB0AC",   # 11 浅灰
  
  # 深色系在后
  "#4E79A7",   # 12 深蓝
  "#F28E2B",   # 13 橙
  "#59A14F",   # 14 深绿
  "#E15759",   # 15 红
  "#499894",   # 16 青绿
  "#D37295",   # 17 玫红
  "#B6992D",   # 18 土黄
  "#9D7660",   # 19 棕
  "#79706E"    # 20 深灰
)

cell_types <- as.character(unique(Malignant@meta.data$sub_cell_type))
num_legend_items <- length(cell_types)  # 图例的个数
max_label_length <- max(nchar(cell_types))  # 图例名称的最大长度

base_width <- 3000  # 基础宽度
base_height <- 3000  # 基础高度
legend_width_factor <- 100  # 每个图例项增加的宽度
label_length_factor <- 10  # 每个字符增加的宽度
dynamic_width <- base_width + (num_legend_items * legend_width_factor) + (max_label_length * label_length_factor)
pdf("malignant_annotation.pdf", width = dynamic_width/300, height = base_height/300)
DimPlot(Malignant, reduction = "umap", label = TRUE, pt.size = 1, group.by = "sub_cell_type", label.size = 4) +
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    ggtitle(NULL) +
    scale_color_manual(values = cluster_colors) +
    coord_fixed(ratio = 1) +
    guides(color = guide_legend(title = NULL, override.aes = list(size = 5))) +
    theme(
        text = element_text(size = 8, face = "bold"),
        axis.text.x = element_text(size = 28, color = "black"),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title.x = element_text(size = 36, face = "bold", color = "black"),
        axis.title.y = element_text(size = 36, face = "bold", color = "black", margin = margin(r = 20)),  # 增加右侧间距
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 16, face = "bold", color = "black"),
        legend.title = element_text(size = 16, face = "bold", color = "black"),
        legend.position = "right",
        legend.box.margin = margin(0, 0, 0, 0),
        legend.key = element_blank(),
        legend.background = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5),
        aspect.ratio = 1,
        plot.margin = margin(10, 50, 10, 10)
    )
dev.off()


#比例
library(tidyverse)

# 计算各亚群在不同组织中的占比
prop_data <- Malignant@meta.data %>%
  group_by(tissue, sub_cell_type) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(tissue) %>%
  mutate(proportion = count / sum(count) * 100) %>%
  mutate(tissue = factor(tissue, levels = c("mPT", "nmPT","metLN")))

library(ggplot2)
library(ggsci)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(14)


prop_data <- prop_data %>%
  mutate(sub_cell_type = as.character(sub_cell_type))

# 绘制堆叠条形图
p <- ggplot(prop_data, aes(x = tissue, y = proportion, fill = sub_cell_type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = cluster_colors) +
  labs(title = "Cell Type Distribution: mPT vs nmPT",
       x = "Tissue", y = "Proportion (%)", fill = "Cell Type") +
  theme_minimal() +
  theme(
    # 坐标轴线条和刻度
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.15, "cm"),
    # 坐标轴文本
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 10),
    axis.text.y = element_text(color = "black", size = 10),
    # 坐标轴标题
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    # 面板背景
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    # 图例
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )


# 保存为PDF
pdf("Malignant_celltype_distribution_tissue.pdf", width = 4, height = 6)
print(p)
dev.off()


p <- FeaturePlot(Malignant,feature="VIM")
ggsave("feature_vim.png",plot=p)
p <- VlnPlot(Malignant,feature="VIM")
ggsave("vln_vim.png",plot=p)

#patient
sample_patient <- data.frame(
  sample = c("A1", "A2", "A3", "A4", "A5", 
             "B1", "B2", "B3", "B4", "B5",
             "C1", "C2", "C3", "C4", "C5",
             "D1", "D2", "D3", "D4", "D5"),
  patient = c("P1_nmPT", "P1_negLN", "P2_nmPT", "P2_negLN", "P3_nmPT",
             "P3_negLN", "P4_nmPT", "P4_negLN", "P5_mPT", "P5_negLN",
             "P5_metLN", "P6_mPT", "P6_negLN", "P6_metLN", "P7_mPT",
             "P7_negLN", "P7_metLN", "P8_mPT", "P8_negLN", "P8_metLN")
)
Malignant$patient <- sample_patient$patient[match(Malignant$sample, sample_patient$sample)]

library(Seurat)
library(pheatmap)

# ==================== 1. 保留所有组织（mPT, nmPT, metLN）====================
Malignant_all <- Malignant  # 直接使用原对象

# ==================== 2. 计算每个患者各亚型的比例 ====================
patient_subtype_table <- table(Malignant_all$patient, Malignant_all$sub_cell_type)
patient_subtype_prop <- prop.table(patient_subtype_table, margin = 1) * 100

# ==================== 3. 筛选主要亚型（平均占比 > 1%）====================
avg_prop <- colMeans(patient_subtype_prop)
keep_types <- names(avg_prop[avg_prop > 1])
patient_subtype_prop_filtered <- patient_subtype_prop[, keep_types]

# ==================== 4. 按患者名称排序 ====================
patient_subtype_prop_filtered <- patient_subtype_prop_filtered[order(rownames(patient_subtype_prop_filtered)), ]

# ==================== 5. 添加组织类型注释 ====================
# 从 metadata 中提取每个患者对应的组织类型
patient_tissue <- Malignant_all@meta.data %>%
  group_by(patient) %>%
  summarise(tissue = first(tissue), .groups = "drop") %>%
  arrange(patient)

# 确保顺序与热图行名一致
patient_tissue <- patient_tissue[match(rownames(patient_subtype_prop_filtered), patient_tissue$patient), ]

annotation_row <- data.frame(Tissue = patient_tissue$tissue)
rownames(annotation_row) <- patient_tissue$patient

# 组织颜色
annotation_colors <- list(
  Tissue = c("mPT" = "#E41A1C", "nmPT" = "#4DAF4A", "metLN" = "#377EB8")
)

# ==================== 6. 绘制热图 ====================
pdf("patient_subtype_heatmap_all_tissues.pdf", width = 12, height = 8)
pheatmap(patient_subtype_prop_filtered,
         scale = "none",
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         annotation_row = annotation_row,
         annotation_colors = annotation_colors,
         main = "Tumor subtype composition: mPT vs nmPT vs metLN",
         xlab = "Tumor subtype", 
         ylab = "Patient",
         fontsize = 8,
         fontsize_row = 8,
         fontsize_col = 8,
         color = colorRampPalette(c("white", "lightblue", "steelblue", "darkred"))(100),
         breaks = seq(0, 100, length.out = 101))
dev.off()

```
# trajectory
```R
library(Seurat)
library(monocle)
library(ggplot2)
library(ggsci)
library(igraph)  # 添加这一行
#trace('project2MST', edit = T, where = asNamespace("monocle"))


# 加载完整修复脚本
source("monocle_fix_complete.R")

Malignant <- readRDS("malignant_anno.rds")

# 随机抽样20000个细胞
random_seed <- sample(1:10000, 1)
set.seed(random_seed)
cat("随机种子:", random_seed, "\n")
sample_cells <- sample(colnames(Malignant), size = 20000)
Malignant_sub <- Malignant[, sample_cells]
cat("抽样后细胞数:", ncol(Malignant_sub), "\n")

# 使用 RNA assay
Malignant_sub <- NormalizeData(Malignant_sub, assay = "RNA")
Malignant_sub <- FindVariableFeatures(Malignant_sub, nfeatures = 2000, assay = "RNA")
hvgs <- VariableFeatures(Malignant_sub, assay = "RNA")

# 提取表达矩阵 - 保持稀疏格式（不转换）
expression_matrix <- LayerData(Malignant_sub, assay = "RNA", layer = "counts")
# 不要转成密集矩阵！保持 sparse matrix

cell_metadata <- Malignant_sub@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(expression_matrix))
rownames(gene_annotation) <- rownames(expression_matrix)

# 只保留存在的基因
hvgs_valid <- hvgs[hvgs %in% rownames(expression_matrix)]
cat("有效高变基因数:", length(hvgs_valid), "\n")

# 提取高变基因表达矩阵（保持稀疏）
expression_matrix_hvgs <- expression_matrix[hvgs_valid, ]
gene_annotation_hvgs <- gene_annotation[hvgs_valid, , drop = FALSE]

cat("高变基因数:", nrow(expression_matrix_hvgs), "\n")
cat("细胞数:", ncol(expression_matrix_hvgs), "\n")
cat("矩阵类型:", class(expression_matrix_hvgs)[1], "\n")

# 创建 CellDataSet 对象
cds <- newCellDataSet(expression_matrix_hvgs,
                      phenoData = new("AnnotatedDataFrame", data = cell_metadata),
                      featureData = new("AnnotatedDataFrame", data = gene_annotation_hvgs),
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())

# 归一化
cds <- estimateSizeFactors(cds)

# 估计离散度（如果报错，用 pooled 方法）
cds <- estimateDispersions(cds, method = "pooled")
cds <- setOrderingFilter(cds, hvgs_valid)

# 降维
cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')
gc()
cds <- orderCells(cds)
# 保存
saveRDS(cds, "malignant_cds.rds")


library(ggplot2)
library(ggsci)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(14)

cluster_colors <- c(
  # 浅色系在前（浅蓝第一个）
  "#A0CBE8",   # 1 浅蓝（你要的放第一个）
  "#FFBE7D",   # 2 浅橙
  "#8CD17D",   # 3 浅绿
  "#86BCB6",   # 4 浅青
  "#FF9D9A",   # 5 浅红
  "#FABFD2",   # 6 浅粉
  "#D4A6C8",   # 7 浅紫
  "#F1CE63",   # 8 淡黄
  "#D7B5A6",   # 9 浅棕
  "#B07AA1",   # 10 紫灰
  "#BAB0AC",   # 11 浅灰
  
  # 深色系在后
  "#4E79A7",   # 12 深蓝
  "#F28E2B",   # 13 橙
  "#59A14F",   # 14 深绿
  "#E15759",   # 15 红
  "#499894",   # 16 青绿
  "#D37295",   # 17 玫红
  "#B6992D",   # 18 土黄
  "#9D7660",   # 19 棕
  "#79706E"    # 20 深灰
)


p3 <- plot_cell_trajectory(cds, color_by = "sub_cell_type") + scale_color_manual(values = cluster_colors) +
  theme(legend.text = element_text(size = 18),
        legend.title = element_blank())

png("traj_malignant_celltype.png", width = 6000, height = 3000, res = 300)
print(p3)
dev.off()

p4 <- plot_cell_trajectory(cds, color_by = "sub_cell_type") + facet_wrap("~sub_cell_type", nrow = 2) + scale_color_manual(values = cluster_colors) +
  theme(legend.text = element_text(size = 18),
        legend.title = element_blank())
png("traj_malignant_celltype_split.png", width = 6000, height = 3000, res = 300)
print(p4)
dev.off()

start_color <- npg_extended[1]
end_color <- npg_extended[3]
p5 <- plot_cell_trajectory(cds, color_by = "Pseudotime") +
  scale_color_gradient(low = start_color, high = end_color)

# 保存图片
ggsave("traj_malignant_pseudo.png", plot = p5, width = 16, height = 6, dpi = 300)
```
# cellchat
```R
#全部细胞的空间
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
sample_size <- min(400000, length(x))
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
```
# CAF
```R
fib <- subset(obj,subset=CellType=="Fibroblasts")

fib <- NormalizeData(fib)
fib <- FindVariableFeatures(fib, nfeatures = 2000)
hvgs <- VariableFeatures(fib)
fib <- ScaleData(fib, features = hvgs)
fib <- RunPCA(fib, features = hvgs, npcs = 50)
p <- ElbowPlot(fib, ndims = 30)
ggsave("p3.png",plot=p)

pca_var <- fib[["pca"]]@stdev^2 / sum(fib[["pca"]]@stdev^2)
cumsum_var <- cumsum(pca_var)

# 看达到多少 PC 时累计方差 > 70% / 80% / 90%
which(cumsum_var > 0.7)[1]   # 70%
which(cumsum_var > 0.8)[1]   # 80%
which(cumsum_var > 0.9)[1]   # 90%

fib <- FindNeighbors(fib, dims = 1:30)
fib <- FindClusters(fib, resolution = 0.3)
fib <- RunUMAP(fib, dims = 1:30)

library(ggplot2)
library(ggsci)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(12)

seurat_clusters <- as.character(unique(fib@meta.data$seurat_clusters))  # 转换为字符向量
num_legend_items <- length(seurat_clusters)  # 图例的个数
max_label_length <- max(nchar(seurat_clusters))  # 图例名称的最大长度

base_width <- 3000  # 基础宽度
base_height <- 3000  # 基础高度
legend_width_factor <- 100  # 每个图例项增加的宽度
label_length_factor <- 10  # 每个字符增加的宽度
dynamic_width <- base_width + (num_legend_items * legend_width_factor) + (max_label_length * label_length_factor)
pdf("fib_clusters.pdf", width = dynamic_width/300, height = base_height/300)  # 转换为英寸

DimPlot(fib, reduction = "umap", label = TRUE, pt.size = 1, label.size = 8) +
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    ggtitle(NULL) +
    scale_color_manual(values = npg_pal) +
    coord_fixed(ratio = 1) +
    guides(color = guide_legend(title = NULL, override.aes = list(size = 5))) +
    theme(
        text = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 28, color = "black"),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title.x = element_text(size = 36, face = "bold", color = "black"),
        axis.title.y = element_text(size = 36, face = "bold", color = "black", margin = margin(r = 20)),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 28, face = "bold", color = "black"),
        legend.title = element_text(size = 28, face = "bold", color = "black"),
        legend.position = "right",
        legend.box.margin = margin(0, 0, 0, 0),
        legend.key = element_blank(),
        legend.background = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5),
        aspect.ratio = 1,
        plot.margin = margin(10, 50, 10, 10)
    )

dev.off()


library(dplyr)
fib_markers <- FindAllMarkers(fib, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "wilcox")
fib_significant_markers <- subset(fib_markers, p_val_adj < 0.05)
#write.csv(fib_significant_markers, "fib_all_marker.csv")
fib_significant_markers <- fib_significant_markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
write.csv(fib_significant_markers, "fib_top_marker_50.csv")


identity_mapping <- c(
  "0" = "Fib_myCAF",
  "1" = "Fib_apCAF_1",
  "2" = "Fib_iCAF_1",
  "3" = "Fib_apCAF_2",
  "4" = "Fib_iCAF_2",
  "5" = "Fib_matCAF"
)

sub_cell_type <- identity_mapping[fib@meta.data$seurat_clusters]
fib@meta.data$sub_cell_type <- sub_cell_type

library(ggplot2)
library(ggsci)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(5)
cluster_colors <- c(
  # 浅色系在前（浅蓝第一个）
  "#A0CBE8",   # 1 浅蓝（你要的放第一个）
  "#FFBE7D",   # 2 浅橙
  "#8CD17D",   # 3 浅绿
  "#86BCB6",   # 4 浅青
  "#FF9D9A",   # 5 浅红
  "#FABFD2",   # 6 浅粉
  "#D4A6C8",   # 7 浅紫
  "#F1CE63",   # 8 淡黄
  "#D7B5A6",   # 9 浅棕
  "#B07AA1",   # 10 紫灰
  "#BAB0AC",   # 11 浅灰
  
  # 深色系在后
  "#4E79A7",   # 12 深蓝
  "#F28E2B",   # 13 橙
  "#59A14F",   # 14 深绿
  "#E15759",   # 15 红
  "#499894",   # 16 青绿
  "#D37295",   # 17 玫红
  "#B6992D",   # 18 土黄
  "#9D7660",   # 19 棕
  "#79706E"    # 20 深灰
)
cell_types <- as.character(unique(fib@meta.data$sub_cell_type))
num_legend_items <- length(cell_types)  # 图例的个数
max_label_length <- max(nchar(cell_types))  # 图例名称的最大长度

base_width <- 3000  # 基础宽度
base_height <- 3000  # 基础高度
legend_width_factor <- 100  # 每个图例项增加的宽度
label_length_factor <- 10  # 每个字符增加的宽度
dynamic_width <- base_width + (num_legend_items * legend_width_factor) + (max_label_length * label_length_factor)
pdf("fib_annotation.pdf", width = dynamic_width/300, height = base_height/300)
DimPlot(fib, reduction = "umap", label = TRUE, pt.size = 1, group.by = "sub_cell_type", label.size = 4) +
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    ggtitle(NULL) +
    scale_color_manual(values = ) +
    coord_fixed(ratio = 1) +
    guides(color = guide_legend(title = NULL, override.aes = list(size = 5))) +
    theme(
        text = element_text(size = 8, face = "bold"),
        axis.text.x = element_text(size = 28, color = "black"),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title.x = element_text(size = 36, face = "bold", color = "black"),
        axis.title.y = element_text(size = 36, face = "bold", color = "black", margin = margin(r = 20)),  # 增加右侧间距
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 16, face = "bold", color = "black"),
        legend.title = element_text(size = 16, face = "bold", color = "black"),
        legend.position = "right",
        legend.box.margin = margin(0, 0, 0, 0),
        legend.key = element_blank(),
        legend.background = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5),
        aspect.ratio = 1,
        plot.margin = margin(10, 50, 10, 10)
    )
dev.off()

#比例
fib <- subset(fib,subset=tissue %in% c("mPT", "nmPT","metLN"))
library(tidyverse)

# 计算各亚群在不同组织中的占比
prop_data <- fib@meta.data %>%
  group_by(tissue, sub_cell_type) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(tissue) %>%
  mutate(proportion = count / sum(count) * 100) %>%
  mutate(tissue = factor(tissue, levels = c("mPT", "nmPT", "metLN")))

library(ggplot2)
library(ggsci)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(14)

prop_data <- prop_data %>%
  mutate(sub_cell_type = as.character(sub_cell_type))

# 绘制堆叠条形图
p <- ggplot(prop_data, aes(x = tissue, y = proportion, fill = sub_cell_type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = cluster_colors) +
  labs(title = "Fibroblast Subtype Distribution: mPT vs nmPT vs metLN",
       x = "Tissue", y = "Proportion (%)", fill = "Cell Type") +
  theme_minimal() +
  theme(
    # 坐标轴线条和刻度
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.15, "cm"),
    # 坐标轴文本
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 10),
    axis.text.y = element_text(color = "black", size = 10),
    # 坐标轴标题
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    # 面板背景
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    # 图例
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )

# 保存为PDF
pdf("Fibroblast_subtype_distribution_tissue.pdf", width = 4, height = 6)
print(p)
dev.off()

#sample
library(tidyverse)
library(ggplot2)

# 计算各样本中各成纤维亚型的比例
prop_data_sample <- fib@meta.data %>%
  filter(tissue %in% c("mPT", "nmPT", "metLN")) %>%
  group_by(tissue, sample, sub_cell_type) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(tissue, sample) %>%
  mutate(proportion = count / sum(count) * 100) %>%
  ungroup()

# 按组织类型排序样本（可根据需要调整）
# 先获取各组织中的样本列表
samples_mPT <- prop_data_sample %>% filter(tissue == "mPT") %>% pull(sample) %>% unique() %>% sort()
samples_nmPT <- prop_data_sample %>% filter(tissue == "nmPT") %>% pull(sample) %>% unique() %>% sort()
samples_metLN <- prop_data_sample %>% filter(tissue == "metLN") %>% pull(sample) %>% unique() %>% sort()

# 合并为统一的样本顺序：mPT → nmPT → metLN
sample_order <- c(samples_mPT, samples_nmPT, samples_metLN)

# 将样本转换为因子，按指定顺序排列
prop_data_sample$sample <- factor(prop_data_sample$sample, levels = sample_order)

# 将 sub_cell_type 转换为因子（可选，控制图例顺序）
prop_data_sample$sub_cell_type <- as.character(prop_data_sample$sub_cell_type)

# 绘图：按样本展示成纤维亚型组成
p_sample <- ggplot(prop_data_sample, aes(x = sample, y = proportion, fill = sub_cell_type)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  scale_fill_manual(values = cluster_colors, name = "Subtype") +
  labs(title = "Fibroblast subtype composition by sample",
       x = "Sample", y = "Proportion (%)") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank()
  )

# 保存图片
ggsave("Fibroblast_subtype_composition_by_sample.pdf", plot = p_sample, width = 12, height = 6, dpi = 300)
```
# CAF亚型和malignant空间分布
```R
#全部细胞的空间
obj@meta.data$sub_cell_type <- NA
common_cells <- intersect(rownames(Malignant@meta.data), rownames(obj@meta.data))
obj@meta.data[common_cells, "sub_cell_type"] <- Malignant@meta.data[common_cells, "sub_cell_type"]

cat("成功注释的细胞数：", sum(!is.na(obj@meta.data$sub_cell_type)), "\n")
cat("obj中sub_cell_type列的唯一值：", paste(unique(obj@meta.data$sub_cell_type), collapse = ", "), "\n")

common_cells <- intersect(rownames(fib@meta.data), rownames(obj@meta.data))
obj@meta.data[common_cells, "sub_cell_type"] <- fib@meta.data[common_cells, "sub_cell_type"]

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


library(ggplot2)
library(dplyr)
library(RColorBrewer)

# ===================== 第一步：数据准备 + 统一添加分面变量 =====================
spatial_data <- obj@meta.data %>%
  select(
    CenterX_global_px, CenterY_global_px,
    detailed, sample, fov, tissue,
    Area.um2, Circularity, AspectRatio, NucArea,
    Mean.PanCK, Mean.CD45, Mean.DAPI
  ) %>%
  filter(!is.na(CenterX_global_px), !is.na(CenterY_global_px)) %>%
  # 1. 筛选目标tissue
  filter(tissue %in% c("metLN", "negLN",  "mPT", "nmPT")) %>%
  filter(!is.na(tissue)) %>%
  # 2. 标注plot_type
  mutate(
    # 先分出大类
    cell_category = case_when(
      # 肿瘤亚型（全部归为"Tumor cells"）
      detailed %in% c("Inflamed_2", "Inflamed_4", "Metabolic_1", "Stem-like",
                      "Differentiated", "Inflamed_3", "Club cell", "EMT-like",
                      "Inflamed_1", "Proliferative_1", "Proliferative_2", "Metabolic_2") ~ "Tumor cells",
      
      # 成纤维亚型（保留具体亚型名称）
      detailed %in% c("apCAF", "eCAF", "myCAF", "matCAF", "iCAF") ~ detailed,
      
      # 其他细胞
      TRUE ~ "Other cells"
    ),
    # 3. 创建tissue_sample分面变量
    tissue = factor(tissue, levels = c("metLN", "negLN",  "mPT", "nmPT")),
    tissue_sample = paste(tissue, sample, sep = "_")
  )

# 定义颜色
# 肿瘤细胞统一为红色
tumor_color <- c("Tumor cells" = "#FF0000")

# 成纤维亚型颜色（使用Set3调色板，每种不同）
fib_types <- c("apCAF", "eCAF", "myCAF", "matCAF", "iCAF")
n_fib <- length(fib_types)

# 为成纤维细胞生成不同颜色
fib_colors <- brewer.pal(min(n_fib, 8), "Set2")
if(n_fib > 8) {
  fib_colors <- colorRampPalette(brewer.pal(8, "Set2"))(n_fib)
}
names(fib_colors) <- fib_types

# 其他细胞颜色
other_color <- c("Other cells" = "#DCDCDC")

# 合并所有颜色
all_colors <- c(tumor_color, fib_colors, other_color)

# 调整图例顺序
legend_order <- c("Tumor cells", fib_types, "Other cells")

# ===================== 第二步：分组1绘图（metLN, negLN） =====================
group1_tissues <- c("metLN", "negLN")
group1_data <- spatial_data %>% filter(tissue %in% group1_tissues)

# 检查分组1数据
if(nrow(group1_data) == 0) {
  stop("分组1（metLN/negLN）无数据！")
}

# 按tissue_sample排序
group1_data$tissue_sample <- factor(
  group1_data$tissue_sample,
  levels = unique(group1_data$tissue_sample[order(group1_data$tissue, group1_data$sample)])
)

# 计算分面行数（每行3个sample）
total_facets1 <- length(unique(group1_data$tissue_sample))
n_rows1 <- ceiling(total_facets1 / 3)

# 绘图
p_group1 <- ggplot(group1_data,
                   aes(x = CenterX_global_px, y = CenterY_global_px, 
                       color = cell_category, size = Area.um2)) +
  # 第1层：Other cells（灰色）
  geom_point(data = filter(group1_data, cell_category == "Other cells"),
             alpha = 0.3) +
  # 第2层：Tumor cells（红色）
  geom_point(data = filter(group1_data, cell_category == "Tumor cells"),
             alpha = 0.8, size = 0.5) +
  # 第3层：成纤维细胞（各亚型不同颜色）
  geom_point(data = filter(group1_data, cell_category %in% fib_types),
             alpha = 0.9, size = 0.5) +
  # 样式设置
  scale_size_continuous(range = c(0.1, 0.5), guide = "none") +
  scale_color_manual(values = all_colors, breaks = legend_order) +
  facet_wrap(~tissue_sample, scales = "free", ncol = 3, nrow = n_rows1) +
  theme_minimal() +
  labs(title = "Cell Distribution - Lymph Nodes",
       x = "X Coordinate (um)",
       y = "Y Coordinate (um)",
       color = "Cell Type") +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    strip.text = element_text(size = 10, face = "bold"),
    strip.background = element_rect(fill = "gray90"),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9, face = "bold"),
    axis.text = element_text(size = 6),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_line(color = "gray95")
  )

# 保存分组1PDF
pdf("fib_cell_distribution_LN.pdf", width = 15, height = 15)
print(p_group1)
dev.off()


# ===================== 第三步：分组2绘图（mPT, nmPT） =====================
group2_tissues <- c("mPT", "nmPT")
group2_data <- spatial_data %>% filter(tissue %in% group2_tissues)

if(nrow(group2_data) > 0) {
  group2_data$tissue_sample <- factor(
    group2_data$tissue_sample,
    levels = unique(group2_data$tissue_sample[order(group2_data$tissue, group2_data$sample)])
  )
  
  total_facets2 <- length(unique(group2_data$tissue_sample))
  n_rows2 <- ceiling(total_facets2 / 3)
  
  p_group2 <- ggplot(group2_data,
                     aes(x = CenterX_global_px, y = CenterY_global_px, 
                         color = cell_category, size = Area.um2)) +
    geom_point(data = filter(group2_data, cell_category == "Other cells"),
               alpha = 0.3) +
    geom_point(data = filter(group2_data, cell_category == "Tumor cells"),
               alpha = 0.8, size = 0.5) +
    geom_point(data = filter(group2_data, cell_category %in% fib_types),
               alpha = 0.9, size = 0.5) +
    scale_size_continuous(range = c(0.1, 0.5), guide = "none") +
    scale_color_manual(values = all_colors, breaks = legend_order) +
    facet_wrap(~tissue_sample, scales = "free", ncol = 3, nrow = n_rows2) +
    theme_minimal() +
    labs(title = "Cell Distribution - Primary Tumor",
         x = "X Coordinate (um)",
         y = "Y Coordinate (um)",
         color = "Cell Type") +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      strip.text = element_text(size = 10, face = "bold"),
      strip.background = element_rect(fill = "gray90"),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 9, face = "bold"),
      axis.text = element_text(size = 6),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_line(color = "gray95")
    )
  
  pdf("fib_cell_distribution_PT.pdf", width = 15, height = 15)
  print(p_group2)
  dev.off()
}
```
# 肿瘤亚型与成纤维的空间分布
```R
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ggsci)

# ===================== 第一步：数据准备 + 统一添加分面变量 =====================
spatial_data <- obj@meta.data %>%
  select(
    CenterX_global_px, CenterY_global_px,
    detailed, sample, fov, tissue,
    Area.um2, Circularity, AspectRatio, NucArea,
    Mean.PanCK, Mean.CD45, Mean.DAPI
  ) %>%
  filter(!is.na(CenterX_global_px), !is.na(CenterY_global_px)) %>%
  # 1. 筛选目标tissue
  filter(tissue %in% c("metLN", "negLN",  "mPT", "nmPT")) %>%
  filter(!is.na(tissue)) %>%
  # 2. 标注plot_type
  mutate(
    # 先分出大类
    cell_category = case_when(
      # 肿瘤亚型（保留具体亚型名称）
      detailed %in% c("Inflamed_2", "Inflamed_4", "Metabolic_1", "Stem-like",
                      "Differentiated", "Inflamed_3", "Club cell", "EMT-like",
                      "Inflamed_1", "Proliferative_1", "Proliferative_2", "Metabolic_2") ~ detailed,
      
      # 成纤维亚型（统一为"Fibroblasts"）
      detailed %in% c("apCAF", "eCAF", "myCAF", "matCAF", "iCAF") ~ "Fibroblasts",
      
      # 其他细胞
      TRUE ~ "Other cells"
    ),
    # 3. 创建tissue_sample分面变量
    tissue = factor(tissue, levels = c("metLN", "negLN",  "mPT", "nmPT")),
    tissue_sample = paste(tissue, sample, sep = "_")
  )

# 定义颜色
# 肿瘤亚型（使用 NPG 调色板）
tumor_types <- c("Inflamed_2", "Inflamed_4", "Metabolic_1", "Stem-like",
                 "Differentiated", "Inflamed_3", "Club cell", "EMT-like",
                 "Inflamed_1", "Proliferative_1", "Proliferative_2", "Metabolic_2")
n_tumor <- length(tumor_types)

# 使用 NPG 调色板生成肿瘤颜色
npg_pal <- pal_npg()(10)  # NPG 默认有10种颜色
tumor_colors <- colorRampPalette(npg_pal)(n_tumor)  # 扩展到需要的数量
names(tumor_colors) <- tumor_types

# 成纤维细胞统一为橙色/棕色系（与肿瘤的蓝紫色系形成鲜明对比）
fib_color <- c("Fibroblasts" = "#E69F00")  # 橙黄色（NPG中的标准橙色）

# 其他细胞颜色
other_color <- c("Other cells" = "#DCDCDC")

# 合并所有颜色
all_colors <- c(tumor_colors, fib_color, other_color)

# 调整图例顺序
legend_order <- c(tumor_types, "Fibroblasts", "Other cells")

# ===================== 第二步：分组1绘图（metLN, negLN） =====================
group1_tissues <- c("metLN", "negLN")
group1_data <- spatial_data %>% filter(tissue %in% group1_tissues)

# 检查分组1数据
if(nrow(group1_data) == 0) {
  stop("分组1（metLN/negLN）无数据！")
}

# 按tissue_sample排序
group1_data$tissue_sample <- factor(
  group1_data$tissue_sample,
  levels = unique(group1_data$tissue_sample[order(group1_data$tissue, group1_data$sample)])
)

# 计算分面行数（每行3个sample）
total_facets1 <- length(unique(group1_data$tissue_sample))
n_rows1 <- ceiling(total_facets1 / 3)

# 绘图
p_group1 <- ggplot(group1_data,
                   aes(x = CenterX_global_px, y = CenterY_global_px, 
                       color = cell_category, size = Area.um2)) +
  # 第1层：Other cells（灰色）
  geom_point(data = filter(group1_data, cell_category == "Other cells"),
             alpha = 0.3) +
  # 第2层：肿瘤细胞（NPG颜色）
  geom_point(data = filter(group1_data, cell_category %in% tumor_types),
             alpha = 0.8, size = 0.5) +
  # 第3层：成纤维细胞（橙色）
  geom_point(data = filter(group1_data, cell_category == "Fibroblasts"),
             alpha = 0.9, size = 0.5) +
  # 样式设置
  scale_size_continuous(range = c(0.1, 0.5), guide = "none") +
  scale_color_manual(values = all_colors, breaks = legend_order) +
  facet_wrap(~tissue_sample, scales = "free", ncol = 3, nrow = n_rows1) +
  theme_minimal() +
  labs(title = "Cell Distribution - Lymph Nodes",
       x = "X Coordinate (um)",
       y = "Y Coordinate (um)",
       color = "Cell Type") +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    strip.text = element_text(size = 10, face = "bold"),
    strip.background = element_rect(fill = "gray90"),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9, face = "bold"),
    axis.text = element_text(size = 6),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_line(color = "gray95")
  )

# 保存分组1PDF
pdf("tumor_fib_orange_LN.pdf", width = 15, height = 15)
print(p_group1)
dev.off()

# ===================== 第三步：分组2绘图（mPT, nmPT） =====================
group2_tissues <- c("mPT", "nmPT")
group2_data <- spatial_data %>% filter(tissue %in% group2_tissues)

if(nrow(group2_data) > 0) {
  group2_data$tissue_sample <- factor(
    group2_data$tissue_sample,
    levels = unique(group2_data$tissue_sample[order(group2_data$tissue, group2_data$sample)])
  )
  
  total_facets2 <- length(unique(group2_data$tissue_sample))
  n_rows2 <- ceiling(total_facets2 / 3)
  
  p_group2 <- ggplot(group2_data,
                     aes(x = CenterX_global_px, y = CenterY_global_px, 
                         color = cell_category, size = Area.um2)) +
    geom_point(data = filter(group2_data, cell_category == "Other cells"),
               alpha = 0.3) +
    geom_point(data = filter(group2_data, cell_category %in% tumor_types),
               alpha = 0.8, size = 0.5) +
    geom_point(data = filter(group2_data, cell_category == "Fibroblasts"),
               alpha = 0.9, size = 0.5) +
    scale_size_continuous(range = c(0.1, 0.5), guide = "none") +
    scale_color_manual(values = all_colors, breaks = legend_order) +
    facet_wrap(~tissue_sample, scales = "free", ncol = 3, nrow = n_rows2) +
    theme_minimal() +
    labs(title = "Cell Distribution - Primary Tumor",
         x = "X Coordinate (um)",
         y = "Y Coordinate (um)",
         color = "Cell Type") +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      strip.text = element_text(size = 10, face = "bold"),
      strip.background = element_rect(fill = "gray90"),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 9, face = "bold"),
      axis.text = element_text(size = 6),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_line(color = "gray95")
    )
  
  pdf("tumor_fib_orange_PT.pdf", width = 15, height = 15)
  print(p_group2)
  dev.off()
}
```
# 单个肿瘤样本中代谢型的分布
```R
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ggsci)

# ===================== 数据准备 =====================
spatial_data <- obj@meta.data %>%
  select(
    CenterX_global_px, CenterY_global_px,
    detailed, sample, fov, tissue,
    Area.um2, Circularity, AspectRatio, NucArea,
    Mean.PanCK, Mean.CD45, Mean.DAPI
  ) %>%
  filter(!is.na(CenterX_global_px), !is.na(CenterY_global_px)) %>%
  # 筛选目标tissue_sample
  mutate(tissue_sample = paste(tissue, sample, sep = "_")) %>%
  filter(tissue_sample == "mPT_D3") %>%
  # 标注cell_category
  mutate(
    cell_category = case_when(
      # 代谢型肿瘤（分别标注）
      detailed == "Metabolic_1" ~ "Metabolic_1",
      detailed == "Metabolic_2" ~ "Metabolic_2",
      
      # 其他肿瘤亚型（统一为一种）
      detailed %in% c("Inflamed_2", "Inflamed_4", "Stem-like",
                      "Differentiated", "Inflamed_3", "Club cell", "EMT-like",
                      "Inflamed_1", "Proliferative_1", "Proliferative_2") ~ "Other_tumor",
      
      # 成纤维细胞（统一）
      detailed %in% c("apCAF", "eCAF", "myCAF", "matCAF", "iCAF") ~ "Fibroblasts",
      
      # 巨噬细胞（精确匹配）
      detailed == "Macrophages" ~ "Macrophages",
      
      # T细胞（精确匹配）
      detailed == "T cells" ~ "T_cells",
      
      # 内皮细胞（精确匹配）
      detailed == "Endothelial cells" ~ "Endothelial",
      
      # 其他细胞
      TRUE ~ "Other cells"
    )
  )

# 定义颜色
# 代谢型肿瘤（用鲜明颜色区分）
tumor_colors <- c(
  "Metabolic_1" = "#E41A1C",  # 红色
  "Metabolic_2" = "#FF7F00"   # 橙色
)

# 其他肿瘤（用深灰色，区别于其他细胞）
other_tumor_color <- c("Other_tumor" = "#8B3A3A")  # 深红褐色

# 成纤维细胞（橙色系）
fib_color <- c("Fibroblasts" = "#f4e1b6")

# 巨噬细胞（绿色系）
macro_color <- c("Macrophages" = "#33A02C")  # 绿色

# T细胞（蓝色系）
tcell_color <- c("T_cells" = "#1F78B4")  # 蓝色

# 内皮细胞（紫色系）
endo_color <- c("Endothelial" = "#6A3D9A")  # 紫色

# 其他细胞
other_color <- c("Other cells" = "#DCDCDC")

# 合并所有颜色
all_colors <- c(tumor_colors, other_tumor_color, fib_color, 
                macro_color, tcell_color, endo_color, other_color)

# 调整图例顺序
legend_order <- c("Metabolic_1", "Metabolic_2", "Other_tumor", 
                  "Fibroblasts", "Macrophages", "T_cells", "Endothelial", 
                  "Other cells")

# 绘图
p <- ggplot(spatial_data,
            aes(x = CenterX_global_px, y = CenterY_global_px, 
                color = cell_category, size = Area.um2)) +
  # 第1层：Other cells（灰色）
  geom_point(data = filter(spatial_data, cell_category == "Other cells"),
             alpha = 0.3) +
  # 第2层：其他肿瘤细胞
  geom_point(data = filter(spatial_data, cell_category == "Other_tumor"),
             alpha = 0.7, size = 0.5) +
  # 第3层：巨噬细胞
  geom_point(data = filter(spatial_data, cell_category == "Macrophages"),
             alpha = 0.7, size = 0.5) +
  # 第4层：T细胞
  geom_point(data = filter(spatial_data, cell_category == "T_cells"),
             alpha = 0.7, size = 0.5) +
  # 第5层：内皮细胞
  geom_point(data = filter(spatial_data, cell_category == "Endothelial"),
             alpha = 0.7, size = 0.5) +
  # 第6层：成纤维细胞
  geom_point(data = filter(spatial_data, cell_category == "Fibroblasts"),
             alpha = 0.9, size = 0.6) +
  # 第7层（最顶层）：代谢型肿瘤（突出显示）
  geom_point(data = filter(spatial_data, cell_category %in% c("Metabolic_1", "Metabolic_2")),
             alpha = 0.9, size = 0.7) +
  # 样式设置
  scale_size_continuous(range = c(0.1, 0.5), guide = "none") +
  scale_color_manual(values = all_colors, breaks = legend_order) +
  theme_minimal() +
  labs(title = "Cell Distribution - mPT_D3",
       x = "X Coordinate (um)",
       y = "Y Coordinate (um)",
       color = "Cell Type") +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9, face = "bold"),
    axis.text = element_text(size = 6),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_line(color = "gray95")
  )

# 保存PDF
pdf("tumor_mPT_D3_metabolic_focus.pdf", width = 12, height = 10)
print(p)
dev.off()
print("Saved: tumor_mPT_D3_metabolic_focus.pdf")

#fov
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ggsci)

# ===================== 数据准备 =====================
spatial_data <- obj@meta.data %>%
  select(
    CenterX_global_px, CenterY_global_px,
    detailed, sample, fov, tissue,
    Area.um2, Circularity, AspectRatio, NucArea,
    Mean.PanCK, Mean.CD45, Mean.DAPI
  ) %>%
  filter(!is.na(CenterX_global_px), !is.na(CenterY_global_px)) %>%
  # 筛选目标tissue_sample和fov
  mutate(tissue_sample = paste(tissue, sample, sep = "_")) %>%
  filter(tissue_sample == "mPT_D3",
         fov %in% c(155, 159)) %>%
  # 标注cell_category
  mutate(
    cell_category = case_when(
      # 代谢型肿瘤（分别标注）
      detailed == "Metabolic_1" ~ "Metabolic_1",
      detailed == "Metabolic_2" ~ "Metabolic_2",
      
      # 其他肿瘤亚型（统一为一种）
      detailed %in% c("Inflamed_2", "Inflamed_4", "Stem-like",
                      "Differentiated", "Inflamed_3", "Club cell", "EMT-like",
                      "Inflamed_1", "Proliferative_1", "Proliferative_2") ~ "Other_tumor",
      
      # 成纤维细胞（统一）
      detailed %in% c("apCAF", "eCAF", "myCAF", "matCAF", "iCAF") ~ "Fibroblasts",
      
      # 巨噬细胞（精确匹配）
      detailed == "Macrophages" ~ "Macrophages",
      
      # T细胞（精确匹配）
      detailed == "T cells" ~ "T_cells",
      
      # 内皮细胞（精确匹配）
      detailed == "Endothelial cells" ~ "Endothelial",
      
      # 其他细胞
      TRUE ~ "Other cells"
    )
  )

# 定义颜色
all_colors <- c(
  "Metabolic_1" = "#E41A1C",      # 红色
  "Metabolic_2" = "#FF7F00",      # 橙色
  "Other_tumor" = "#8B3A3A",      # 深红褐色
  "Fibroblasts" = "#fcdf9f",      # 橙黄
  "Macrophages" = "#33A02C",      # 绿色
  "T_cells" = "#1F78B4",          # 蓝色
  "Endothelial" = "#6A3D9A",      # 紫色
  "Other cells" = "#DCDCDC"       # 浅灰
)

# 图例顺序
legend_order <- c("Metabolic_1", "Metabolic_2", "Other_tumor", 
                  "Fibroblasts", "Macrophages", "T_cells", "Endothelial", 
                  "Other cells")

# 检查数据
print(paste("筛选后细胞数量:", nrow(spatial_data)))
print("fov分布:")
print(table(spatial_data$fov))

# 计算整体边界
global_boundary <- data.frame(
  xmin = min(spatial_data$CenterX_global_px, na.rm = TRUE),
  xmax = max(spatial_data$CenterX_global_px, na.rm = TRUE),
  ymin = min(spatial_data$CenterY_global_px, na.rm = TRUE),
  ymax = max(spatial_data$CenterY_global_px, na.rm = TRUE)
)

print("整体边界:")
print(global_boundary)

# 计算右上角位置（用于放置FOV文本标签）
x_right <- global_boundary$xmax
y_top <- global_boundary$ymax
x_offset <- (global_boundary$xmax - global_boundary$xmin) * 0.02
y_offset <- (global_boundary$ymax - global_boundary$ymin) * 0.02

# 创建FOV标签数据框（放在右上角）
fov_labels <- data.frame(
  x = x_right - x_offset,
  y = y_top - y_offset,
  label = paste("FOV 155 & 159")
)

# 绘图 - 交换x和y坐标实现水平旋转，使用aes(color)自动生成图例
p <- ggplot() +
  # 第1层：所有点（使用color映射自动生成图例）
  geom_point(data = spatial_data,
             aes(x = CenterY_global_px, y = CenterX_global_px, 
                 color = cell_category),
             alpha = 0.7, size = 2.5) +
  # 添加整体边界框
  geom_rect(data = global_boundary,
            aes(xmin = ymin, xmax = ymax, ymin = xmin, ymax = xmax),
            fill = NA, color = "black", linewidth = 0.8, linetype = "solid") +
  # 添加FOV文本标签
  geom_text(data = fov_labels,
            aes(x = y, y = x, label = label),
            size = 5, color = "black", fontface = "bold", hjust = 1, vjust = 1) +
  # 颜色映射
  scale_color_manual(values = all_colors, breaks = legend_order) +
  # 样式设置
  coord_fixed(ratio = 1) +
  theme_minimal() +
  labs(title = "Cell Distribution - mPT_D3 (FOV 155 & 159 Combined)",
       x = "Y Coordinate (um) - Rotated",
       y = "X Coordinate (um) - Rotated",
       color = "Cell Type") +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_line(color = "gray95"),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10, face = "bold")
  )

# 保存PDF
pdf("tumor_mPT_D3_FOV155_159_combined_rotated.pdf", width = 14, height = 10)
print(p)
dev.off()
print("Saved: tumor_mPT_D3_FOV155_159_combined_rotated.pdf")

```
# fibro单个样本nmPT_A1
```R
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ggsci)

# ===================== 数据准备 =====================
spatial_data <- obj@meta.data %>%
  select(
    CenterX_global_px, CenterY_global_px,
    detailed, sample, fov, tissue,
    Area.um2, Circularity, AspectRatio, NucArea,
    Mean.PanCK, Mean.CD45, Mean.DAPI
  ) %>%
  filter(!is.na(CenterX_global_px), !is.na(CenterY_global_px)) %>%
  # 筛选目标tissue_sample
  mutate(tissue_sample = paste(tissue, sample, sep = "_")) %>%
  filter(tissue_sample == "nmPT_A1") %>%
  # 标注cell_category
  mutate(
    cell_category = case_when(
      # 所有肿瘤亚型（统一为一种）
      detailed %in% c("Metabolic_1", "Metabolic_2", "Inflamed_2", "Inflamed_4", 
                      "Stem-like", "Differentiated", "Inflamed_3", "Club cell", 
                      "EMT-like", "Inflamed_1", "Proliferative_1", "Proliferative_2") ~ "Tumor cells",
      
      # myCAF单独一种颜色
      detailed == "myCAF" ~ "myCAF",
      
      # 其他成纤维亚型统一为一种
      detailed %in% c("apCAF", "eCAF", "matCAF", "iCAF") ~ "Other CAFs",
      
      # 巨噬细胞（精确匹配）
      detailed == "Macrophages" ~ "Macrophages",
      
      # T细胞（精确匹配）
      detailed == "T cells" ~ "T_cells",
      
      # 内皮细胞（精确匹配）
      detailed == "Endothelial cells" ~ "Endothelial",
      
      # 其他细胞
      TRUE ~ "Other cells"
    )
  )

# 定义颜色
all_colors <- c(
  # 肿瘤统一为一种颜色（深红色）
  "Tumor cells" = "#B22234",        # 深红色
  
  # myCAF单独颜色（橙红色，突出显示）
  "myCAF" = "#D55E00",              # 橙红色
  
  # 其他CAF统一颜色（浅橙色）
  "Other CAFs" = "#F0E442",         # 黄色/浅橙色
  
  # 其他免疫细胞
  "Macrophages" = "#33A02C",       # 绿色
  "T_cells" = "#1F78B4",           # 蓝色
  "Endothelial" = "#6A3D9A",       # 紫色
  "Other cells" = "#DCDCDC"        # 浅灰
)

# 图例顺序
legend_order <- c("Tumor cells", "myCAF", "Other CAFs",
                  "Macrophages", "T_cells", "Endothelial", "Other cells")

# 检查数据
print(paste("筛选后细胞数量:", nrow(spatial_data)))
print("样本信息:")
print(paste("样本:", unique(spatial_data$tissue_sample)))

# 计算整体边界
global_boundary <- data.frame(
  xmin = min(spatial_data$CenterX_global_px, na.rm = TRUE),
  xmax = max(spatial_data$CenterX_global_px, na.rm = TRUE),
  ymin = min(spatial_data$CenterY_global_px, na.rm = TRUE),
  ymax = max(spatial_data$CenterY_global_px, na.rm = TRUE)
)

print("整体边界:")
print(global_boundary)

# 计算右上角位置（用于放置样本标签）
x_right <- global_boundary$xmax
y_top <- global_boundary$ymax
x_offset <- (global_boundary$xmax - global_boundary$xmin) * 0.02
y_offset <- (global_boundary$ymax - global_boundary$ymin) * 0.02

# 创建样本标签数据框（放在右上角）
sample_label <- data.frame(
  x = x_right - x_offset,
  y = y_top - y_offset,
  label = "nmPT_A1 (Whole Section)"
)

# 绘图 - 不旋转，点调小
p <- ggplot() +
  # 所有点（使用color映射自动生成图例）
  geom_point(data = spatial_data,
             aes(x = CenterX_global_px, y = CenterY_global_px, 
                 color = cell_category),
             alpha = 0.6, size = 0.8) +
  # 添加整体边界框
  geom_rect(data = global_boundary,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = NA, color = "black", linewidth = 0.8, linetype = "solid") +
  # 添加样本标签
  geom_text(data = sample_label,
            aes(x = x, y = y, label = label),
            size = 4, color = "black", fontface = "bold", hjust = 1, vjust = 1) +
  # 颜色映射
  scale_color_manual(values = all_colors, breaks = legend_order) +
  # 样式设置
  coord_fixed(ratio = 1) +
  theme_minimal() +
  labs(title = "Cell Distribution - nmPT_A1 (Whole Section)",
       x = "X Coordinate (um)",
       y = "Y Coordinate (um)",
       color = "Cell Type") +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_line(color = "gray95"),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10, face = "bold")
  )

# 保存PDF
pdf("nmPT_A1_whole_section.pdf", width = 12, height = 10)
print(p)
dev.off()
print("Saved: nmPT_A1_whole_section.pdf")
```
# fibro放大fov
```R
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ggsci)

# ===================== 数据准备 =====================
spatial_data <- obj@meta.data %>%
  select(
    CenterX_global_px, CenterY_global_px,
    detailed, sample, fov, tissue,
    Area.um2, Circularity, AspectRatio, NucArea,
    Mean.PanCK, Mean.CD45, Mean.DAPI
  ) %>%
  filter(!is.na(CenterX_global_px), !is.na(CenterY_global_px)) %>%
  # 筛选目标tissue_sample和fov
  mutate(tissue_sample = paste(tissue, sample, sep = "_")) %>%
  filter(tissue_sample == "mPT_D3",
         fov %in% c(155, 159)) %>%   # 改为 FOV 155 和 159
  # 标注cell_category
  mutate(
    cell_category = case_when(
      # 所有肿瘤亚型（统一为一种）
      detailed %in% c("Metabolic_1", "Metabolic_2", "Inflamed_2", "Inflamed_4", 
                      "Stem-like", "Differentiated", "Inflamed_3", "Club cell", 
                      "EMT-like", "Inflamed_1", "Proliferative_1", "Proliferative_2") ~ "Tumor cells",
      
      # eCAF单独一种颜色
      detailed == "eCAF" ~ "eCAF",
      
      # apCAF单独一种颜色
      detailed == "apCAF" ~ "apCAF",
      
      # 其他成纤维亚型统一为一种
      detailed %in% c("myCAF", "matCAF", "iCAF") ~ "Other CAFs",
      
      # 巨噬细胞（精确匹配）
      detailed == "Macrophages" ~ "Macrophages",
      
      # T细胞（精确匹配）
      detailed == "T cells" ~ "T_cells",
      
      # 内皮细胞（精确匹配）
      detailed == "Endothelial cells" ~ "Endothelial",
      
      # 其他细胞
      TRUE ~ "Other cells"
    )
  )

# 定义颜色
all_colors <- c(
  # 肿瘤统一为一种颜色（深红色）
  "Tumor cells" = "#B22234",        # 深红色
  
  # eCAF单独颜色（橙色）
  "eCAF" = "#E69F00",               # 橙色
  
  # apCAF单独颜色（橙红色）
  "apCAF" = "#D55E00",              # 橙红色/赭色
  
  # 其他CAF统一颜色（浅黄色）
  "Other CAFs" = "#F0E442",         # 黄色/浅橙色
  
  # 其他免疫细胞
  "Macrophages" = "#33A02C",       # 绿色
  "T_cells" = "#1F78B4",           # 蓝色
  "Endothelial" = "#6A3D9A",       # 紫色
  "Other cells" = "#DCDCDC"        # 浅灰
)

# 图例顺序
legend_order <- c("Tumor cells", "eCAF", "apCAF", "Other CAFs",
                  "Macrophages", "T_cells", "Endothelial", "Other cells")

# 检查数据
print(paste("筛选后细胞数量:", nrow(spatial_data)))
print("fov分布:")
print(table(spatial_data$fov))
print("cell_category分布:")
print(table(spatial_data$cell_category))

# 计算整体边界
global_boundary <- data.frame(
  xmin = min(spatial_data$CenterX_global_px, na.rm = TRUE),
  xmax = max(spatial_data$CenterX_global_px, na.rm = TRUE),
  ymin = min(spatial_data$CenterY_global_px, na.rm = TRUE),
  ymax = max(spatial_data$CenterY_global_px, na.rm = TRUE)
)

print("整体边界:")
print(global_boundary)

# 计算右上角位置（用于放置FOV文本标签）
x_right <- global_boundary$xmax
y_top <- global_boundary$ymax
x_offset <- (global_boundary$xmax - global_boundary$xmin) * 0.02
y_offset <- (global_boundary$ymax - global_boundary$ymin) * 0.02

# 创建FOV标签数据框（放在右上角）
fov_labels <- data.frame(
  x = x_right - x_offset,
  y = y_top - y_offset,
  label = paste("FOV 155 & 159")
)

# 绘图 - 交换x和y坐标实现水平旋转，使用aes(color)自动生成图例
p <- ggplot() +
  # 所有点（使用color映射自动生成图例）
  geom_point(data = spatial_data,
             aes(x = CenterY_global_px, y = CenterX_global_px, 
                 color = cell_category),
             alpha = 0.7, size = 2.5) +
  # 添加整体边界框
  geom_rect(data = global_boundary,
            aes(xmin = ymin, xmax = ymax, ymin = xmin, ymax = xmax),
            fill = NA, color = "black", linewidth = 0.8, linetype = "solid") +
  # 添加FOV文本标签
  geom_text(data = fov_labels,
            aes(x = y, y = x, label = label),
            size = 5, color = "black", fontface = "bold", hjust = 1, vjust = 1) +
  # 颜色映射
  scale_color_manual(values = all_colors, breaks = legend_order) +
  # 样式设置
  coord_fixed(ratio = 1) +
  theme_minimal() +
  labs(title = "Cell Distribution - mPT_D3 (FOV 155 & 159 Combined)",
       x = "Y Coordinate (um) - Rotated",
       y = "X Coordinate (um) - Rotated",
       color = "Cell Type") +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_line(color = "gray95"),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10, face = "bold")
  )

# 保存PDF
pdf("mPT_D3_FOV155_159_combined_rotated.pdf", width = 14, height = 10)
print(p)
dev.off()
print("Saved: mPT_D3_FOV155_159_combined_rotated.pdf")
```
# Macrophage
```R
macro <- subset(obj,subset=CellType=="Macrophages")
macro <- subset(macro,subset=seurat_clusters %in% c(0,1,2,3,4))
macro <- NormalizeData(macro)
macro <- FindVariableFeatures(macro, nfeatures = 3000)
hvgs <- VariableFeatures(macro)
macro <- ScaleData(macro, features = hvgs)
macro <- RunPCA(macro, features = hvgs, npcs = 50)
p <- ElbowPlot(macro, ndims = 30)
ggsave("elbow_macro.png",plot=p)
pca_var <- macro[["pca"]]@stdev^2 / sum(macro[["pca"]]@stdev^2)
cumsum_var <- cumsum(pca_var)

# 看达到多少 PC 时累计方差 > 70% / 80% / 90%
which(cumsum_var > 0.7)[1]   # 70%
which(cumsum_var > 0.8)[1]   # 80%
which(cumsum_var > 0.9)[1]   # 90%

macro <- FindNeighbors(macro, dims = 1:20)
macro <- FindClusters(macro, resolution = 0.4)
macro <- RunUMAP(macro, dims = 1:20)

library(ggplot2)
library(ggsci)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(11)

seurat_clusters <- as.character(unique(macro@meta.data$seurat_clusters))  # 转换为字符向量
num_legend_items <- length(seurat_clusters)  # 图例的个数
max_label_length <- max(nchar(seurat_clusters))  # 图例名称的最大长度

base_width <- 3000  # 基础宽度
base_height <- 3000  # 基础高度
legend_width_factor <- 100  # 每个图例项增加的宽度
label_length_factor <- 10  # 每个字符增加的宽度
dynamic_width <- base_width + (num_legend_items * legend_width_factor) + (max_label_length * label_length_factor)
pdf("macro_clusters.pdf", width = dynamic_width/300, height = base_height/300)  # 转换为英寸

DimPlot(macro, reduction = "umap", label = TRUE, pt.size = 1, label.size = 8) +
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    ggtitle(NULL) +
    scale_color_manual(values = npg_pal) +
    coord_fixed(ratio = 1) +
    guides(color = guide_legend(title = NULL, override.aes = list(size = 5))) +
    theme(
        text = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 28, color = "black"),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title.x = element_text(size = 36, face = "bold", color = "black"),
        axis.title.y = element_text(size = 36, face = "bold", color = "black", margin = margin(r = 20)),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 28, face = "bold", color = "black"),
        legend.title = element_text(size = 28, face = "bold", color = "black"),
        legend.position = "right",
        legend.box.margin = margin(0, 0, 0, 0),
        legend.key = element_blank(),
        legend.background = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5),
        aspect.ratio = 1,
        plot.margin = margin(10, 50, 10, 10)
    )

dev.off()

library(dplyr)
macro_markers <- FindAllMarkers(macro, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "wilcox")
macro_significant_markers <- subset(macro_markers, p_val_adj < 0.05)
#write.csv(macro_significant_markers, "macro_all_marker.csv")
macro_significant_markers <- macro_significant_markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
write.csv(macro_significant_markers, "macro_top_marker_50.csv")

identity_mapping <- c(
  "0" = "Lipid_Macro",
  "1" = "Stress-responsive_Macro",
  "2" = "Homeostatic_Macro",
  "3" = "Regulatory_Macro",
  "4" = "Matrix-remodeling_Macro"
)

identity_mapping <- c(
  "0" = "CHIT1_TAM",
  "1" = "APOE_TAM",
  "2" = "CCL19_TAM",
  "3" = "CSF1R_TAM",
  "4" = "SIGMAR1_TAM",
  "5" = "SPP1_TAM"
)

identity_mapping <- c(
  "0" = "Lipid-associated_Macro_1",
  "1" = "Lipid-associated_Macro_2",
  "2" = "Inflammatory_Macro",
  "3" = "Matrix-remodeling_Macro",
  "4" = "Stress-response_Macro",
  "5" = "Proliferative-Macro"
)
identity_mapping <- c(
  "0" = "Macro_Lipid-associated",
  "1" = "Macro_Antigen-presenting",
  "2" = "Macro_Immunoregulatory",
  "3" = "Macro_Stress-response",
  "4" = "Macro_Matrix-remodeling",
  "5" = "Macro_Homeostatic",
  "6" = "Macro_Phagocytic" 
)

sub_cell_type <- identity_mapping[macro@meta.data$seurat_clusters]
macro@meta.data$sub_cell_type <- sub_cell_type

library(ggplot2)
library(ggsci)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(7)

cluster_colors <- c(
  # 浅色系在前（浅蓝第一个）
  "#A0CBE8",   # 1 浅蓝（你要的放第一个）
  "#FFBE7D",   # 2 浅橙
  "#8CD17D",   # 3 浅绿
  "#86BCB6",   # 4 浅青
  "#FF9D9A",   # 5 浅红
  "#FABFD2",   # 6 浅粉
  "#D4A6C8",   # 7 浅紫
  "#F1CE63",   # 8 淡黄
  "#D7B5A6",   # 9 浅棕
  "#B07AA1",   # 10 紫灰
  "#BAB0AC",   # 11 浅灰
  
  # 深色系在后
  "#4E79A7",   # 12 深蓝
  "#F28E2B",   # 13 橙
  "#59A14F",   # 14 深绿
  "#E15759",   # 15 红
  "#499894",   # 16 青绿
  "#D37295",   # 17 玫红
  "#B6992D",   # 18 土黄
  "#9D7660",   # 19 棕
  "#79706E"    # 20 深灰
)


cell_types <- as.character(unique(macro@meta.data$sub_cell_type))
num_legend_items <- length(cell_types)  # 图例的个数
max_label_length <- max(nchar(cell_types))  # 图例名称的最大长度

base_width <- 3000  # 基础宽度
base_height <- 3000  # 基础高度
legend_width_factor <- 100  # 每个图例项增加的宽度
label_length_factor <- 10  # 每个字符增加的宽度
dynamic_width <- base_width + (num_legend_items * legend_width_factor) + (max_label_length * label_length_factor)
pdf("macro_annotation.pdf", width = dynamic_width/300, height = base_height/300)
DimPlot(macro, reduction = "umap", label = TRUE, pt.size = 1.5, group.by = "sub_cell_type", label.size = 4) +
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    ggtitle(NULL) +
    scale_color_manual(values = cluster_colors) +
    coord_fixed(ratio = 1) +
    guides(color = guide_legend(title = NULL, override.aes = list(size = 5))) +
    theme(
        text = element_text(size = 8, face = "bold"),
        axis.text.x = element_text(size = 28, color = "black"),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title.x = element_text(size = 36, face = "bold", color = "black"),
        axis.title.y = element_text(size = 36, face = "bold", color = "black", margin = margin(r = 20)),  # 增加右侧间距
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 16, face = "bold", color = "black"),
        legend.title = element_text(size = 16, face = "bold", color = "black"),
        legend.position = "right",
        legend.box.margin = margin(0, 0, 0, 0),
        legend.key = element_blank(),
        legend.background = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5),
        aspect.ratio = 1,
        plot.margin = margin(10, 50, 10, 10)
    )
dev.off()

#比例
library(tidyverse)

# 计算各亚群在不同组织中的占比
prop_data <- macro@meta.data %>%
  group_by(tissue, sub_cell_type) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(tissue) %>%
  mutate(proportion = count / sum(count) * 100) %>%
  filter(tissue %in% c("mPT", "nmPT")) %>%
  mutate(tissue = factor(tissue, levels = c("mPT", "nmPT")))

library(ggplot2)
library(ggsci)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(7)


prop_data <- prop_data %>%
  mutate(sub_cell_type = as.character(sub_cell_type))

# 绘制堆叠条形图
p <- ggplot(prop_data, aes(x = tissue, y = proportion, fill = sub_cell_type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = cluster_colors) +
  labs(title = "Cell Type Distribution: mPT vs nmPT",
       x = "Tissue", y = "Proportion (%)", fill = "Cell Type") +
  theme_minimal() +
  theme(
    # 坐标轴线条和刻度
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.15, "cm"),
    # 坐标轴文本
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 10),
    axis.text.y = element_text(color = "black", size = 10),
    # 坐标轴标题
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    # 面板背景
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    # 图例
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )


# 保存为PDF
pdf("macro_celltype_distribution_tumor.pdf", width = 4, height = 6)
print(p)
dev.off()

prop_data <- macro@meta.data %>%
  group_by(tissue, sub_cell_type) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(tissue) %>%
  mutate(proportion = count / sum(count) * 100) %>%
  filter(tissue %in% c("metLN", "negLN")) %>%
  mutate(tissue = factor(tissue, levels = c("metLN", "negLN")))

library(ggplot2)
library(ggsci)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(7)


prop_data <- prop_data %>%
  mutate(sub_cell_type = as.character(sub_cell_type))

# 绘制堆叠条形图
p <- ggplot(prop_data, aes(x = tissue, y = proportion, fill = sub_cell_type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = cluster_colors) +
  labs(title = "Cell Type Distribution: metLN vs negLN",
       x = "Tissue", y = "Proportion (%)", fill = "Cell Type") +
  theme_minimal() +
  theme(
    # 坐标轴线条和刻度
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.15, "cm"),
    # 坐标轴文本
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 10),
    axis.text.y = element_text(color = "black", size = 10),
    # 坐标轴标题
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    # 面板背景
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    # 图例
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )


# 保存为PDF
pdf("macro_celltype_distribution_LN.pdf", width = 4, height = 6)
print(p)
dev.off()


```
# mPT的主要细胞类型互作（肿瘤亚型与其他主要类型）
```R
# ===================== 加载必要的包 =====================
library(Seurat)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(deldir)
library(doParallel)
library(foreach)

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
                       #"Malignant cells" = "unknown",
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

mPT <- subset(obj,subset=tissue=="mPT")
#不分样本看这里，分样本去看chat_mPT.r
coords <- mPT@meta.data[, c("CenterX_global_px", "CenterY_global_px")]
head(coords)

# 获取细胞barcodes（就是行名）
barcodes <- rownames(coords)

labels <- as.character(mPT@meta.data$detailed)
names(labels) <- barcodes

label_names <- names(labels)

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
rw <- c(min(x), max(x), min(y), max(y))

deld <- deldir(x, y, rw = rw)

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

mat_fc <- mat_obs / (mat_mu + 1e-8)
mat_log2fc <- log2(mat_fc)
#diag(mat_log2fc) <- NA

# 查看范围
range(mat_log2fc, na.rm = TRUE)

# 裁剪到合理范围
#mat_log2fc[mat_log2fc > 3] <- 3
#mat_log2fc[mat_log2fc < -3] <- -3
mat_log2fc_clean <- mat_log2fc

# 将 -Inf 替换为 -10（或根据你的数据范围调整）
mat_log2fc_clean[is.infinite(mat_log2fc_clean) & mat_log2fc_clean < 0] <- -10
# 绘制热图
p <- pheatmap(mat_log2fc_clean,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Log2(Fold Change) of cell-cell contacts",
         fontsize = 10,
         na_col = "grey90")
#ggsave("mPT_spatial_contact_log2fc_maligvsothers.pdf", plot = p, width = 12, height = 10)

# 保存结果
save(mat_obs, mat_z, mat_p, mat_mu, file = "mPT_spatial_contact_analysis_results_maligvsothers.RData")

#显著性
# 创建显著性标注矩阵
signif_symbols <- matrix("", nrow = K, ncol = K, dimnames = list(types, types))

# 设置显著性阈值
signif_symbols[mat_p < 0.001] <- "***"
signif_symbols[mat_p < 0.01 & mat_p >= 0.001] <- "**"
signif_symbols[mat_p < 0.05 & mat_p >= 0.01] <- "*"

# 绘制带显著性标注的热图
p4 <- pheatmap(mat_log2fc_clean,
         cluster_rows = TRUE, 
         cluster_cols = TRUE,
         main = "Log2(Fold Change) of contact enrichment - mPT",
         fontsize = 10,
         display_numbers = signif_symbols,  # 显示显著性星号
         number_color = "black",           # 星号颜色
         fontsize_number = 12)             # 星号大小
ggsave("mPT_spatial_contact_signif_maligvsothers_log2.pdf",plot=p4,width=10,height=8)

metabolic1 <- rownames(mat_log2fc_clean)[grepl("Metabolic_1", rownames(mat_log2fc_clean))]
metabolic2 <- rownames(mat_log2fc_clean)[grepl("Metabolic_2", rownames(mat_log2fc_clean))]

# 提取互作
inter_met1 <- mat_log2fc_clean[metabolic1, ]
inter_met2 <- mat_log2fc_clean[metabolic2, ]

# 排除自身和彼此
inter_met1 <- inter_met1[!names(inter_met1) %in% c(metabolic1, metabolic2)]
inter_met2 <- inter_met2[!names(inter_met2) %in% c(metabolic1, metabolic2)]

# 计算差异
common_cells <- intersect(names(inter_met1), names(inter_met2))
diff <- inter_met2[common_cells] - inter_met1[common_cells]

# 排序
top_met2 <- sort(diff, decreasing = TRUE)[1:16]   # 代谢2更强的
top_met1 <- sort(diff, decreasing = FALSE)[1:4]  # 代谢1更强的

cat("\n========== 代谢2更强的互作 ==========\n")
print(top_met2)

cat("\n========== 代谢1更强的互作 ==========\n")
print(top_met1)

library(ggplot2)
library(patchwork)

# 创建数据框
df_met2 <- data.frame(
  cell_type = names(top_met2),
  diff = as.numeric(top_met2),
  stronger = "Metabolic2"
)

df_met1 <- data.frame(
  cell_type = names(top_met1),
  diff = as.numeric(top_met1),
  stronger = "Metabolic1"
)

df_plot <- rbind(df_met2, df_met1)

# 条形图 + 方框
p <- ggplot(df_plot, aes(x = reorder(cell_type, diff), y = diff, fill = stronger)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("Metabolic1" = "#E69F00", "Metabolic2" = "#D55E00")) +
  labs(title = "Interaction differences: Metabolic1 vs Metabolic2",
       x = "Cell type", y = "Log2FC difference (Metabolic2 - Metabolic1)",
       fill = "Stronger with") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),  # 添加方框
    plot.title = element_text(hjust = 0.5)
  )

ggsave("metabolic1_vs_metabolic2_bar.pdf", plot = p, width = 10, height = 8)


#炎症型在mPT和nmPT
library(ggplot2)

# 读取结果
mPT_res <- readRDS("mPT_integrated_results_detailed.rds")
nmPT_res <- readRDS("nmPT_integrated_results_detailed.rds")

# 免疫炎症亚型名称（根据实际行名）
I1 <- "Tumor_Immune-inflamed_1"
I2 <- "Tumor_Immune-inflamed_2"

# 获取两个矩阵中共同的细胞类型（用于比较）
common_cells <- intersect(rownames(mPT_res$avg_log2fc), rownames(nmPT_res$avg_log2fc))
# 排除 I1 和 I2 自身
other_common <- setdiff(common_cells, c(I1, I2))

cat("共有的其他细胞类型数量：", length(other_common), "\n")

# 提取互作强度（只针对共有细胞类型）
I1_mPT <- mPT_res$avg_log2fc[I1, other_common]
I2_mPT <- mPT_res$avg_log2fc[I2, other_common]
delta_mPT <- I2_mPT - I1_mPT

I1_nmPT <- nmPT_res$avg_log2fc[I1, other_common]
I2_nmPT <- nmPT_res$avg_log2fc[I2, other_common]
delta_nmPT <- I2_nmPT - I1_nmPT

# 构建数据框
df_mPT <- data.frame(cell_type = other_common, delta = delta_mPT)
df_nmPT <- data.frame(cell_type = other_common, delta = delta_nmPT)

# 去除 NA
df_mPT <- df_mPT[!is.na(df_mPT$delta), ]
df_nmPT <- df_nmPT[!is.na(df_nmPT$delta), ]

# 选择正负两端显著差异的细胞类型（可自行调整数量）
N_pos <- 10   # I2 更强的数量
N_neg <- 5    # I1 更强的数量

# mPT 图数据
pos_mPT <- head(df_mPT[order(df_mPT$delta, decreasing = TRUE), ], N_pos)
neg_mPT <- head(df_mPT[order(df_mPT$delta, decreasing = FALSE), ], N_neg)
plot_mPT <- rbind(pos_mPT, neg_mPT)
plot_mPT$direction <- ifelse(plot_mPT$delta > 0, I2, I1)

# nmPT 图数据
pos_nmPT <- head(df_nmPT[order(df_nmPT$delta, decreasing = TRUE), ], N_pos)
neg_nmPT <- head(df_nmPT[order(df_nmPT$delta, decreasing = FALSE), ], N_neg)
plot_nmPT <- rbind(pos_nmPT, neg_nmPT)
plot_nmPT$direction <- ifelse(plot_nmPT$delta > 0, I2, I1)

# 绘图函数
plot_delta <- function(df, title) {
  df$direction <- factor(df$direction, levels = c(I1, I2))
  ggplot(df, aes(x = reorder(cell_type, delta), y = delta, fill = direction)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_manual(values = setNames(c("#E69F00", "#D55E00"), c(I1, I2)),
                      name = "Stronger in this context") +
    labs(title = title, x = "Cell type", y = paste0(I2, " - ", I1, " (log2FC)")) +
    theme_minimal() +
    theme(legend.position = "bottom",
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
          plot.title = element_text(hjust = 0.5, face = "bold"))
}
# 保存图片
p_mPT <- plot_delta(plot_mPT, "Interaction strength comparison in mPT")
ggsave("mPT_I2_vs_I1_delta.pdf", p_mPT, width = 10, height = 8)

p_nmPT <- plot_delta(plot_nmPT, "Interaction strength comparison in nmPT")
ggsave("nmPT_I2_vs_I1_delta.pdf", p_nmPT, width = 10, height = 8)


# ===================== 比较 mPT vs nmPT =====================

# 1. 加载两个整合结果文件
mPT_results <- readRDS("mPT_integrated_results_detailed.rds")
nmPT_results <- readRDS("nmPT_integrated_results_detailed.rds")

# 2. 提取平均 log2FC 矩阵和 p 值矩阵
mat_log2fc_mPT <- mPT_results$avg_log2fc
mat_p_mPT <- mPT_results$combined_p

mat_log2fc_nmPT <- nmPT_results$avg_log2fc
mat_p_nmPT <- nmPT_results$combined_p

# 3. 找出两种组织共有的细胞类型
common_types <- intersect(rownames(mat_log2fc_mPT), rownames(mat_log2fc_nmPT))
cat("共有细胞类型数量：", length(common_types), "\n")

if (length(common_types) >= 2) {
  
  # 4. 计算差异矩阵（mPT - nmPT）
  mat_diff <- mat_log2fc_mPT[common_types, common_types] - 
              mat_log2fc_nmPT[common_types, common_types]
  
  # 5. 处理缺失值
  mat_diff[is.na(mat_diff)] <- 0
  
  # 6. 处理 -Inf 和 +Inf，然后裁剪
  # 替换 -Inf 为 -10
  mat_diff[is.infinite(mat_diff) & mat_diff < 0] <- -10
  # 替换 +Inf 为 10
  mat_diff[is.infinite(mat_diff) & mat_diff > 0] <- 10
  
  # 裁剪超出 [-10, 10] 的范围
  mat_diff[mat_diff > 10] <- 10
  mat_diff[mat_diff < -10] <- -10
  
  # 7. 计算差异的显著性（Fisher's method 合并 p 值）
  combined_chi2 <- matrix(0, nrow = length(common_types), ncol = length(common_types),
                          dimnames = list(common_types, common_types))
  
  for (i in seq_along(common_types)) {
    for (j in seq_along(common_types)) {
      p1 <- mat_p_mPT[common_types[i], common_types[j]]
      p2 <- mat_p_nmPT[common_types[i], common_types[j]]
      # 处理 NA 或无效 p 值
      if (is.na(p1) || is.na(p2) || p1 <= 0 || p2 <= 0) {
        combined_chi2[i, j] <- 0
      } else {
        combined_chi2[i, j] <- -2 * (log(p1 + 1e-10) + log(p2 + 1e-10))
      }
    }
  }
  
  # 8. 计算合并 p 值（自由度 = 4）
  combined_p <- pchisq(combined_chi2, df = 4, lower.tail = FALSE)
  combined_p[is.na(combined_p)] <- 1
  
  # 9. 显著性标注
  signif_detail <- matrix("", nrow = length(common_types), ncol = length(common_types),
                          dimnames = list(common_types, common_types))
  
  # mPT 中更强的交互（正差异）
  signif_detail[combined_p < 0.001 & mat_diff > 0] <- "***"
  signif_detail[combined_p < 0.01 & combined_p >= 0.001 & mat_diff > 0] <- "**"
  signif_detail[combined_p < 0.05 & combined_p >= 0.01 & mat_diff > 0] <- "*"
  
  # nmPT 中更强的交互（负差异）
  signif_detail[combined_p < 0.001 & mat_diff < 0] <- "***"
  signif_detail[combined_p < 0.01 & combined_p >= 0.001 & mat_diff < 0] <- "**"
  signif_detail[combined_p < 0.05 & combined_p >= 0.01 & mat_diff < 0] <- "*"
  
  # 10. 绘制热图
  library(pheatmap)
  
  p <- pheatmap(mat_diff,
                cluster_rows = TRUE,
                cluster_cols = TRUE,
                display_numbers = signif_detail,
                number_color = "black",
                fontsize_number = 8,
                main = "Difference in cell-cell interactions (mPT - nmPT)\nRed: stronger in mPT, Blue: stronger in nmPT\n* p < 0.05, ** p < 0.01, *** p < 0.001",
                fontsize = 10,
                breaks = seq(-10, 10, length.out = 100),
                color = colorRampPalette(c("blue", "white", "red"))(100),
                filename = "mPT_vs_nmPT_interaction_diff.pdf",
                width = 14, height = 12)
  
  cat("保存差异热图：mPT_vs_nmPT_interaction_diff.pdf\n")
  
  # 保存差异矩阵和 p 值
  write.csv(mat_diff, file = "mPT_vs_nmPT_interaction_diff_matrix.csv")
  write.csv(combined_p, file = "mPT_vs_nmPT_interaction_diff_pvalues.csv")
  
} else {
  cat("共有细胞类型不足 2 个，无法绘制热图\n")
}

# 7. 找出差异最大的细胞对
diff_vec <- as.vector(mat_diff)
names(diff_vec) <- paste(rep(common_types, each = length(common_types)), 
                         rep(common_types, times = length(common_types)), sep = "_")

top_up <- sort(diff_vec, decreasing = TRUE)[1:10]  # mPT 更强的
top_down <- sort(diff_vec, decreasing = FALSE)[1:10]  # nmPT 更强的

cat("\n========== mPT 中互作更强的细胞对 ==========\n")
print(top_up)

cat("\n========== nmPT 中互作更强的细胞对 ==========\n")
print(top_down)

# ===================== 找出炎症型 =====================
cell_types_mPT <- rownames(mat_log2fc_mPT)
cell_types_nmPT <- rownames(mat_log2fc_nmPT)

immuneinflamed1 <- cell_types_mPT[grepl("Immune-inflamed_1", cell_types_mPT)]
immuneinflamed2 <- cell_types_mPT[grepl("Immune-inflamed_2", cell_types_mPT)]

cat("代谢1:", immuneinflamed1, "\n")
cat("代谢2:", immuneinflamed2, "\n")

# ===================== 提取互作 =====================
# mPT
inter_mPT_met1 <- mat_log2fc_mPT[immuneinflamed1, ]
inter_mPT_met2 <- mat_log2fc_mPT[immuneinflamed2, ]

# nmPT
inter_nmPT_met1 <- mat_log2fc_nmPT[immuneinflamed1, ]
inter_nmPT_met2 <- mat_log2fc_nmPT[immuneinflamed2, ]

# 排除自身
inter_mPT_met1 <- inter_mPT_met1[!names(inter_mPT_met1) %in% c(immuneinflamed1, immuneinflamed2)]
inter_mPT_met2 <- inter_mPT_met2[!names(inter_mPT_met2) %in% c(immuneinflamed1, immuneinflamed2)]
inter_nmPT_met1 <- inter_nmPT_met1[!names(inter_nmPT_met1) %in% c(immuneinflamed1, immuneinflamed2)]
inter_nmPT_met2 <- inter_nmPT_met2[!names(inter_nmPT_met2) %in% c(immuneinflamed1, immuneinflamed2)]

# ===================== 找共同细胞类型 =====================
common_cells <- intersect(names(inter_mPT_met1), names(inter_nmPT_met1))

# ===================== 计算差异 =====================
diff_met1 <- inter_mPT_met1[common_cells] - inter_nmPT_met1[common_cells]
diff_met2 <- inter_mPT_met2[common_cells] - inter_nmPT_met2[common_cells]

# 排序
top_met1_mPT <- sort(diff_met1, decreasing = TRUE)[1:15]
top_met2_mPT <- sort(diff_met2, decreasing = TRUE)[1:15]

cat("\n========== 免疫炎症1在 mPT 中更强的互作 ==========\n")
print(top_met1_mPT)

cat("\n========== 免疫炎症2在 mPT 中更强的互作 ==========\n")
print(top_met2_mPT)

final_diff <- diff_met2[common_cells] - diff_met1[common_cells]

#mPT样本中
# 排序
top_met2_stronger <- sort(final_diff, decreasing = TRUE)[1:12]
top_met1_stronger <- sort(final_diff, decreasing = FALSE)[1:3]


cat("\n========== 免疫炎症2更强的互作（基于转移特异性差异）==========\n")
print(top_met2_stronger)

cat("\n========== 免疫炎症1更强的互作（基于转移特异性差异）==========\n")
print(top_met1_stronger)

# ===================== 条形图 =====================
library(ggplot2)

# 创建数据框
df_met2 <- data.frame(
  cell_type = names(top_met2_stronger),
  diff = as.numeric(top_met2_stronger),
  stronger = "Immuneinflamed2"
)

df_met1 <- data.frame(
  cell_type = names(top_met1_stronger),
  diff = as.numeric(top_met1_stronger),
  stronger = "Immuneinflamed1"
)

df_plot <- rbind(df_met2, df_met1)

# 条形
p <- ggplot(df_plot, aes(x = reorder(cell_type, diff), y = diff, fill = stronger)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("Immuneinflamed1" = "#E69F00", "Immuneinflamed2" = "#D55E00")) +
  labs(title = "Immune-inflamed subtypes show divergent interaction strength in mPT",
       x = "Cell type", 
       y = "(Immune-inflamed2_mPT - Immune-inflamed2_nmPT) - (Immune-inflamed1_mPT - Immune-inflamed1_nmPT)",
       fill = "Stronger with") +
  theme_minimal() +
  theme(legend.position = "bottom",
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        plot.title = element_text(hjust = 0.5))

ggsave("Immuneinflamed_divergence_in_mPT.pdf", plot = p, width = 10, height = 8)

#nmPT
# ===================== nmPT 样本中 =====================
# 核心修改：取负号（镜像指标）
# ===================== nmPT 样本中 =====================
# 当前 final_diff_nmPT = -final_diff（即 I1 - I2 方向）
# 要改成 I2 - I1 方向，再取一次负号
final_diff_nmPT_I2_minus_I1 <- -final_diff_nmPT

# 排序：正值 = I2 更强（放上面），负值 = I1 更强（放下面的负方向）
top_I2_in_nmPT <- sort(final_diff_nmPT_I2_minus_I1, decreasing = TRUE)[1:10]
top_I1_in_nmPT <- sort(final_diff_nmPT_I2_minus_I1, decreasing = FALSE)[1:5]

cat("\n========== 在 nmPT 中免疫炎症2更强的互作（正条）==========\n")
print(top_I2_in_nmPT)

cat("\n========== 在 nmPT 中免疫炎症1更强的互作（负条）==========\n")
print(top_I1_in_nmPT)

# ===================== 条形图 =====================
library(ggplot2)

df_I2 <- data.frame(
  cell_type = names(top_I2_in_nmPT),
  diff = as.numeric(top_I2_in_nmPT),
  stronger = "Immuneinflamed2"
)

df_I1 <- data.frame(
  cell_type = names(top_I1_in_nmPT),
  diff = as.numeric(top_I1_in_nmPT),
  stronger = "Immuneinflamed1"
)

df_plot_nmPT <- rbind(df_I2, df_I1)

p_nmPT <- ggplot(df_plot_nmPT, aes(x = reorder(cell_type, diff), y = diff, fill = stronger)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("Immuneinflamed1" = "#E69F00", "Immuneinflamed2" = "#D55E00"),
                    name = "Stronger in nmPT") +
  labs(title = "Immune-inflamed subtypes show divergent interaction strength in nmPT",
       x = "Cell type", 
       y = "(Immune-inflamed2_nmPT - Immune-inflamed2_mPT) - (Immune-inflamed1_nmPT - Immune-inflamed1_mPT)",
       fill = "Stronger in nmPT") +
  theme_minimal() +
  theme(legend.position = "bottom",
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        plot.title = element_text(hjust = 0.5))

ggsave("Immuneinflamed_divergence_in_nmPT.pdf", plot = p_nmPT, width = 10, height = 8)

#代谢
mPT_res <- readRDS("mPT_integrated_results_detailed.rds")

# 确认代谢型名称（根据实际行名修改）
metabolic1 <- "Tumor_Metabolic_1"
metabolic2 <- "Tumor_Metabolic_2"

# 其他细胞类型（排除自身）
others <- setdiff(rownames(mPT_res$avg_log2fc), c(metabolic1, metabolic2))

# 提取互作强度向量并降序排序
met1_vec <- mPT_res$avg_log2fc[metabolic1, others]
met2_vec <- mPT_res$avg_log2fc[metabolic2, others]

met1_sorted <- sort(met1_vec, decreasing = TRUE)
met2_sorted <- sort(met2_vec, decreasing = TRUE)

# 打印结果
cat("===== Metabolic_1 与其他细胞互作强度排序 (log2FC) =====\n")
print(met1_sorted)

cat("\n===== Metabolic_2 与其他细胞互作强度排序 (log2FC) =====\n")
print(met2_sorted)
```
# mPT的subcelltype分布
```R
# ===================== 加载必要的包 =====================
library(Seurat)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(deldir)
library(doParallel)
library(foreach)

# ===================== 数据加载 =====================
obj <- readRDS("YA2025263-1_fin.rds")
Malignant <- readRDS("malignant_unknown_anno_0.6.rds")
fib <- readRDS("fib_anno_new.rds")
macro <- readRDS("macro_anno_new.rds")
mast <- readRDS("mast_anno.rds")
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
                       "Basal cells" = "Malignant cells"
)

# 检查
table(obj$CellType)

obj@meta.data$sub_cell_type <- NA
common_cells <- intersect(rownames(Malignant@meta.data), rownames(obj@meta.data))
obj@meta.data[common_cells, "sub_cell_type"] <- Malignant@meta.data[common_cells, "sub_cell_type"]

common_cells <- intersect(rownames(fib@meta.data), rownames(obj@meta.data))
obj@meta.data[common_cells, "sub_cell_type"] <- fib@meta.data[common_cells, "sub_cell_type"]

common_cells <- intersect(rownames(macro@meta.data), rownames(obj@meta.data))
obj@meta.data[common_cells, "sub_cell_type"] <- macro@meta.data[common_cells, "sub_cell_type"]

common_cells <- intersect(rownames(mast@meta.data), rownames(obj@meta.data))
obj@meta.data[common_cells, "sub_cell_type"] <- mast@meta.data[common_cells, "sub_cell_type"]

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



mPT <- subset(obj,subset=tissue=="mPT")
coords <- mPT@meta.data[, c("CenterX_global_px", "CenterY_global_px")]
head(coords)

# 获取细胞barcodes（就是行名）
barcodes <- rownames(coords)

labels <- as.character(mPT@meta.data$detailed)  # 或者用 obj@meta.data$sub_cell_type
names(labels) <- barcodes

label_names <- names(labels)

# 合并
#labels <- ifelse(labels %in% c("Metabolic_1", "Metabolic_2"), "Metabolic", labels)

# 恢复名称
#names(labels) <- label_names
# 查看合并后的分布
#cat("\n=== 合并后的细胞类型分布 ===\n")
#print(table(labels))

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

#rare_threshold <- 50  # 设置阈值，建议 50

# 计算每种细胞类型的数量
#cell_counts <- table(labels)
#cat("\n=== 细胞类型数量 ===\n")
#print(cell_counts)

# 找出需要保留的细胞类型（数量 >= 阈值）
#keep_types <- names(cell_counts[cell_counts >= rare_threshold])
#cat("\n=== 保留的细胞类型（>= 50个细胞）===\n")
#print(keep_types)

# 过滤细胞
#keep_cells <- labels %in% keep_types
#coords <- coords[keep_cells, , drop = FALSE]
#labels <- labels[keep_cells]

#cat("\n=== 过滤后 ===\n")
#cat("剩余细胞数量:", nrow(coords), "\n")
#cat("细胞类型数量:", length(unique(labels)), "\n")
#print(table(labels))

# 继续处理...
xy <- coords
colnames(xy) <- c("x","y")

print("开始计算Delaunay三角剖分...")
library(deldir)

x <- xy[,1]
y <- xy[,2]
rw <- c(min(x), max(x), min(y), max(y))

deld <- deldir(x, y, rw = rw)

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

mat_fc <- mat_obs / (mat_mu + 1e-8)
mat_log2fc <- log2(mat_fc)
#diag(mat_log2fc) <- NA

# 查看范围
range(mat_log2fc, na.rm = TRUE)

# 裁剪到合理范围
mat_log2fc[mat_log2fc > 3] <- 3
mat_log2fc[mat_log2fc < -3] <- -3

# 绘制热图
p <- pheatmap(mat_log2fc,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Log2(Fold Change) of cell-cell contacts",
         fontsize = 10,
         na_col = "grey90")
ggsave("mPT_spatial_contact_log2fc_majorcelltype.pdf", plot = p, width = 12, height = 10)

# 去除对角线
#mat_z_no_diag <- mat_z
#diag(mat_z_no_diag) <- NA

# z-score 热图
#library(pheatmap)
#p <- pheatmap(mat_z_no_diag,
#         cluster_rows = TRUE, 
#         cluster_cols = TRUE,
#         main = "Z-score of contact enrichment (Obs vs Random) - excluding diagonal",
#         fontsize = 10,
#         na_col = "grey90")
#ggsave("chat_heatmap_mPT_subcelltype_no_diag_range.png",plot=p)
#print("分析完成！")

# 保存结果
save(mat_obs, mat_z, mat_p, mat_mu, file = "mPT_spatial_contact_analysis_results_major_celltype.RData")

#显著性
# 创建显著性标注矩阵
signif_symbols <- matrix("", nrow = K, ncol = K, dimnames = list(types, types))

# 设置显著性阈值
signif_symbols[mat_p < 0.001] <- "***"
signif_symbols[mat_p < 0.01 & mat_p >= 0.001] <- "**"
signif_symbols[mat_p < 0.05 & mat_p >= 0.01] <- "*"

# 绘制带显著性标注的热图
p4 <- pheatmap(mat_log2fc,
         cluster_rows = TRUE, 
         cluster_cols = TRUE,
         main = "Log2(Fold Change) of contact enrichment - mPT",
         fontsize = 10,
         display_numbers = signif_symbols,  # 显示显著性星号
         number_color = "black",           # 星号颜色
         fontsize_number = 12)             # 星号大小
ggsave("mPT_spatial_contact_signif_majorcelltype_log2.pdf",plot=p4,width=10,height=8)
```
# mPT 区域划分与空间分析
```R
library(Seurat)
library(ggplot2)
library(dplyr)
library(pheatmap)

# ===================== 数据加载 =====================
obj <- readRDS("YA2025263-1_fin.rds")
Malignant <- readRDS("malignant_unknown_anno_0.6.rds")
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

# ===================== 提取 mPT =====================
mPT <- subset(obj, subset = tissue == "mPT")
samples <- unique(mPT$sample)
cat("mPT 样本:", paste(samples, collapse = ", "), "\n")

# ===================== 循环1：提取每个样本的数据 =====================
data_list <- list()

for (sample_id in samples) {
  cat("样本:", sample_id, "\n")
  
  sample_obj <- subset(mPT, subset = sample == sample_id)
  
  coords <- sample_obj@meta.data[, c("CenterX_global_px", "CenterY_global_px")]
  coords <- coords[!is.na(coords[,1]) & !is.na(coords[,2]), , drop = FALSE]
  cell_types <- sample_obj$CellType[rownames(coords)]
  
  valid <- !is.na(cell_types)
  coords <- coords[valid, , drop = FALSE]
  cell_types <- cell_types[valid]
  
  data_list[[sample_id]] <- list(
    coords = coords,
    cell_types = cell_types
  )
  
  cat("  细胞数:", nrow(coords), "\n")
}

# ===================== 循环2a：初始化 =====================
region_list <- list()
n_iter <- 20

for (sample_id in samples) {
  cat("\n========== 样本:", sample_id, "==========\n")
  
  coords <- data_list[[sample_id]]$coords
  cell_types <- data_list[[sample_id]]$cell_types
  n_cells <- nrow(coords)
  
  cat("细胞数:", n_cells, "\n")
  
  # 存储每次迭代的区域分配
  all_assignments <- matrix(NA, nrow = n_cells, ncol = n_iter)
  rownames(all_assignments) <- rownames(coords)
  
  # 计算整体坐标范围
  x_min <- min(coords[,1])
  x_max <- max(coords[,1])
  y_min <- min(coords[,2])
  y_max <- max(coords[,2])
  
  region_list[[sample_id]] <- list(
    all_assignments = all_assignments,
    n_cells = n_cells,
    coords = coords,
    cell_types = cell_types,
    x_min = x_min, x_max = x_max,
    y_min = y_min, y_max = y_max
  )
  
  cat("  初始化完成\n")
}

# ===================== 循环2b：多次迭代聚类 =====================
target_cells_per_grid <- 50
n_regions <- 3

for (sample_id in samples) {
  cat("\n样本:", sample_id, "\n")
  
  # 获取数据
  all_assignments <- region_list[[sample_id]]$all_assignments
  coords <- region_list[[sample_id]]$coords
  cell_types <- region_list[[sample_id]]$cell_types
  n_cells <- region_list[[sample_id]]$n_cells
  x_min <- region_list[[sample_id]]$x_min
  x_max <- region_list[[sample_id]]$x_max
  y_min <- region_list[[sample_id]]$y_min
  y_max <- region_list[[sample_id]]$y_max
  
  n_grid <- max(8, round(sqrt(n_cells / target_cells_per_grid)))
  grid_width_x <- (x_max - x_min) / n_grid
  grid_width_y <- (y_max - y_min) / n_grid
  
  for (iter in 1:n_iter) {
    # 随机偏移
    x_offset <- runif(1, 0, grid_width_x)
    y_offset <- runif(1, 0, grid_width_y)
    
    # 自定义区间边界
    breaks_x <- seq(x_min + x_offset, x_max + x_offset, length.out = n_grid + 1)
    breaks_y <- seq(y_min + y_offset, y_max + y_offset, length.out = n_grid + 1)
    
    # 分配网格
    grid_x <- cut(coords[,1], breaks = breaks_x, labels = FALSE, include.lowest = TRUE)
    grid_y <- cut(coords[,2], breaks = breaks_y, labels = FALSE, include.lowest = TRUE)
    grid_id <- paste(grid_x, grid_y, sep = "_")
    
    # 计算网格组成
    all_types <- sort(unique(cell_types))
    unique_grids <- unique(grid_id)
    
    grid_composition <- matrix(0, nrow = length(unique_grids), ncol = length(all_types))
    rownames(grid_composition) <- unique_grids
    colnames(grid_composition) <- all_types
    
    for (i in seq_along(unique_grids)) {
      g <- unique_grids[i]
      cells_in_grid <- which(grid_id == g)
      type_counts <- table(cell_types[cells_in_grid])
      grid_composition[i, names(type_counts)] <- type_counts
    }
    
    grid_composition <- grid_composition / rowSums(grid_composition)
    
    # 聚类
    if (nrow(grid_composition) >= 3) {
      dist_mat <- dist(grid_composition, method = "euclidean")
      hc <- hclust(dist_mat, method = "ward.D2")
      region_labels <- cutree(hc, k = n_regions)
      
      grid_to_region <- setNames(paste0("R", region_labels), rownames(grid_composition))
      
      # 分配区域
      for (i in 1:n_cells) {
        g <- grid_id[i]
        if (g %in% names(grid_to_region)) {
          all_assignments[i, iter] <- grid_to_region[g]
        }
      }
    }
    
    if (iter %% 5 == 0) cat("  迭代", iter, "/", n_iter, "完成\n")
  }
  
  # 保存回列表
  region_list[[sample_id]]$all_assignments <- all_assignments
  cat("  完成\n")
}

# ===================== 循环2c：投票和置信度计算 =====================
final_region_list <- list()

for (sample_id in samples) {
  cat("\n样本:", sample_id, "\n")
  
  all_assignments <- region_list[[sample_id]]$all_assignments
  coords <- region_list[[sample_id]]$coords
  cell_types <- region_list[[sample_id]]$cell_types
  n_cells <- region_list[[sample_id]]$n_cells
  
  # 投票
  final_region <- rep(NA, n_cells)
  names(final_region) <- rownames(coords)
  confidence <- rep(NA, n_cells)
  
  for (i in 1:n_cells) {
    valid_assignments <- all_assignments[i, !is.na(all_assignments[i, ])]
    if (length(valid_assignments) > 0) {
      tab <- table(valid_assignments)
      final_region[i] <- names(tab)[which.max(tab)]
      confidence[i] <- max(tab) / sum(tab)
    }
  }
  
  # 检查是否有差异
  first_iter <- all_assignments[, 1]
  all_same <- TRUE
  for (iter in 2:ncol(all_assignments)) {
    if (!identical(first_iter, all_assignments[, iter])) {
      all_same <- FALSE
      break
    }
  }
  
  if (all_same) {
    cat("  警告：所有迭代结果完全相同\n")
  } else {
    cat("  迭代结果有差异，多轮投票有效\n")
  }
  
  # 置信度分布
  conf_table <- table(round(confidence[!is.na(confidence)], 2))
  cat("  置信度分布:\n")
  print(head(conf_table, 15))
  
  final_region_list[[sample_id]] <- data.frame(
    cell = rownames(coords),
    region = final_region,
    confidence = confidence,
    stringsAsFactors = FALSE
  )
  
  cat("  平均置信度:", round(mean(confidence[!is.na(confidence)]), 3), "\n")
  cat("  完成\n")
}

# 替换原 region_list
region_list <- final_region_list

# ===================== 循环3：可视化 =====================
for (sample_id in samples) {
  cat("\n样本:", sample_id, "\n")
  
  coords <- data_list[[sample_id]]$coords
  region_df <- region_list[[sample_id]]
  
  # 合并区域信息
  region_labels <- setNames(region_df$region, region_df$cell)
  confidence_labels <- setNames(region_df$confidence, region_df$cell)
  
  # 获取该样本的所有细胞
  sample_obj <- subset(mPT, subset = sample == sample_id)
  sample_meta <- sample_obj@meta.data
  
  # 只处理有坐标的细胞
  valid_cells <- rownames(coords)
  
  sample_meta$region <- NA
  sample_meta$confidence <- NA
  sample_meta[valid_cells, "region"] <- region_labels[valid_cells]
  sample_meta[valid_cells, "confidence"] <- confidence_labels[valid_cells]
  
  plot_data <- sample_meta[!is.na(sample_meta$region), ]
  
  if (nrow(plot_data) > 0) {
    # 可视化1：区域组成热图
    region_comp <- table(plot_data$region, plot_data$CellType)
    if (nrow(region_comp) > 1) {
      pheatmap(region_comp,
               main = paste("Cell composition -", sample_id),
               fontsize = 8,
               filename = paste0("mPT_", sample_id, "_region_composition.pdf"),
               width = 10, height = 6)
      cat("  保存热图\n")
    }
    
    # 可视化2：空间区域图（细胞点图）
    p1 <- ggplot(plot_data, aes(x = CenterX_global_px, y = CenterY_global_px, color = region)) +
      geom_point(size = 0.3, alpha = 0.7) +
      coord_fixed() +
      theme_void() +
      labs(title = paste("Spatial regions -", sample_id)) +
      theme(legend.position = "bottom")
    
    ggsave(paste0("mPT_", sample_id, "_spatial_regions.pdf"), 
           plot = p1, width = 8, height = 6)
    cat("  保存空间图\n")
    
    # 可视化3：置信度图
    p2 <- ggplot(plot_data, aes(x = CenterX_global_px, y = CenterY_global_px, color = confidence)) +
      geom_point(size = 0.3, alpha = 0.7) +
      coord_fixed() +
      theme_void() +
      scale_color_gradient(low = "red", high = "blue", limits = c(0, 1)) +
      labs(title = paste("Region confidence -", sample_id)) +
      theme(legend.position = "bottom")
    
    ggsave(paste0("mPT_", sample_id, "_confidence.pdf"), 
           plot = p2, width = 8, height = 6)
    cat("  保存置信度图\n")
  }
}

# ===================== 循环4：整合所有样本的区域热图 =====================
cat("\n========== 整合所有样本 ==========\n")

all_region_compositions <- list()

for (sample_id in samples) {
  region_df <- region_list[[sample_id]]
  cell_types <- data_list[[sample_id]]$cell_types
  
  valid <- !is.na(region_df$region)
  region_comp <- table(region_df$region[valid], cell_types[valid])
  
  if (nrow(region_comp) == 0) next
  
  # 转换为比例
  region_pct <- prop.table(region_comp, margin = 1) * 100
  
  comp_df <- as.data.frame.matrix(region_pct)
  comp_df$sample <- sample_id
  comp_df$region <- rownames(region_pct)
  
  all_region_compositions[[sample_id]] <- comp_df
}

# 合并所有样本
all_comp <- do.call(rbind, all_region_compositions)

# 只保留主要细胞类型（前10种）
cell_cols <- setdiff(colnames(all_comp), c("sample", "region"))
avg_pct <- colMeans(all_comp[, cell_cols], na.rm = TRUE)
top10_types <- names(sort(avg_pct, decreasing = TRUE)[1:10])

# 创建矩阵
heatmap_mat <- as.matrix(all_comp[, top10_types])
rownames(heatmap_mat) <- paste(all_comp$sample, all_comp$region, sep = "_")

# 绘制整合热图
pheatmap(heatmap_mat,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Cell type composition by region across mPT samples",
         fontsize = 8,
         filename = "mPT_all_samples_region_composition.pdf",
         width = 14, height = 12)

cat("保存整合热图: mPT_all_samples_region_composition.pdf\n")

# ===================== 查看每个样本的区域组成 =====================
for (sample_id in samples) {
  cat("\n========== 样本:", sample_id, "==========\n")
  
  region_df <- region_list[[sample_id]]
  cell_types <- data_list[[sample_id]]$cell_types
  
  # 合并数据
  cell_region <- setNames(region_df$region, region_df$cell)
  cell_types_named <- setNames(cell_types, names(cell_types))
  
  # 只保留有区域的细胞
  valid <- !is.na(cell_region)
  region_vec <- cell_region[valid]
  type_vec <- cell_types_named[names(region_vec)]
  
  # 计算每个区域的细胞类型比例
  region_comp <- table(region_vec, type_vec)
  region_pct <- prop.table(region_comp, margin = 1) * 100
  
  # 打印每个区域的主要细胞类型
  for (reg in rownames(region_pct)) {
    pct_sorted <- sort(region_pct[reg, ], decreasing = TRUE)
    top3 <- names(pct_sorted)[1:10]
    top3_pct <- round(pct_sorted[1:10], 1)
    
    cat("\n区域", reg, ":\n")
    cat("  主要细胞类型: ", paste(top3, "(", top3_pct, "%)", collapse = ", "), "\n")
    cat("  细胞数:", sum(region_comp[reg, ]), "\n")
  }
}

# ===================== k=3 统一命名 =====================
# ===================== 先应用区域命名 =====================
region_naming <- list()

# 样本 C5
region_naming[["C5"]] <- c(
  "R1" = "Immune_zone",
  "R2" = "Mixed_zone",
  "R3" = "Immune_zone"
)

# 样本 B4
region_naming[["B4"]] <- c(
  "R1" = "Mixed_zone",
  "R2" = "Tumor_zone",
  "R3" = "Vascular_zone"
)

# 样本 D3
region_naming[["D3"]] <- c(
  "R1" = "Tumor_zone",
  "R2" = "Tumor_zone",
  "R3" = "Stromal_zone"
)

# 样本 C2
region_naming[["C2"]] <- c(
  "R1" = "Stromal_zone",
  "R2" = "Tumor_zone",
  "R3" = "Vascular_zone"
)

# 应用命名
for (sample_id in names(region_list)) {
  if (sample_id %in% names(region_naming)) {
    region_df <- region_list[[sample_id]]
    mapping <- region_naming[[sample_id]]
    
    region_df$region_type <- mapping[region_df$region]
    region_list[[sample_id]] <- region_df
    
    cat("样本", sample_id, "命名完成\n")
    print(table(region_df$region_type))
  }
}

#蛋白表达验证
region_marker_summary <- list()

for (sample_id in names(region_list)) {
  cat("\n样本:", sample_id, "\n")
  
  region_df <- region_list[[sample_id]]
  sample_obj <- subset(mPT, subset = sample == sample_id)
  
  # 获取每个细胞的 PanCK 和 CD45 表达
  panck <- sample_obj$Mean.PanCK
  cd45 <- sample_obj$Mean.CD45
  names(panck) <- rownames(sample_obj@meta.data)
  names(cd45) <- rownames(sample_obj@meta.data)
  
  # 为每个区域计算平均表达
  region_summary <- data.frame()
  
  for (reg in unique(region_df$region_type[!is.na(region_df$region_type)])) {
    cells_in_region <- region_df$cell[region_df$region_type == reg]
    
    avg_panck <- mean(panck[cells_in_region], na.rm = TRUE)
    avg_cd45 <- mean(cd45[cells_in_region], na.rm = TRUE)
    
    region_summary <- rbind(region_summary, data.frame(
      region = reg,
      avg_PanCK = avg_panck,
      avg_CD45 = avg_cd45,
      n_cells = length(cells_in_region)
    ))
  }
  
  print(region_summary)
  region_marker_summary[[sample_id]] <- region_summary
}

mPT$region <- NA
mPT$region_type <- NA

# 逐个样本映射
for (sample_id in names(region_list)) {
  cat("处理样本:", sample_id, "\n")
  
  # 获取该样本的区域信息
  region_df <- region_list[[sample_id]]
  
  # 获取该样本在 mPT 中的细胞
  sample_cells <- rownames(mPT@meta.data[mPT$sample == sample_id, ])
  
  # 只处理有区域信息的细胞
  cells_with_region <- intersect(sample_cells, region_df$cell)
  
  # 映射区域
  for (cell in cells_with_region) {
    cell_region <- region_df$region[region_df$cell == cell]
    cell_region_type <- region_df$region_type[region_df$cell == cell]
    
    mPT@meta.data[cell, "region"] <- cell_region
    mPT@meta.data[cell, "region_type"] <- cell_region_type
  }
  
  cat("  映射细胞数:", length(cells_with_region), "\n")
}

print(table(mPT$region_type))

# ===================== 保存 mPT 对象 =====================
saveRDS(mPT, file = "mPT_with_regions.rds")


library(Seurat)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(deldir)
library(doParallel)
library(foreach)

# ===================== 定义细胞互作分析函数 =====================
run_spatial_analysis <- function(coords_sub, labels_sub, sample_name, region_name, nperm = 500) {
  
  cat("    细胞数:", length(labels_sub), "\n")
  
  # 检查细胞类型数量
  types <- sort(unique(labels_sub))
  K <- length(types)
  
  if (K < 2) {
    cat("    警告：只有", K, "种细胞类型，跳过\n")
    return(NULL)
  }
  
  # 确保坐标是数值矩阵
  coords_sub <- as.matrix(coords_sub)
  x <- coords_sub[,1]
  y <- coords_sub[,2]
  
  # 计算Delaunay三角剖分
  deld <- deldir(x, y, rw = c(range(x), range(y)))
  segs <- deld$delsgs
  
  if (nrow(segs) == 0) {
    cat("    警告：没有三角边，跳过\n")
    return(NULL)
  }
  
  # 构建边
  edges <- cbind(segs$ind1, segs$ind2)
  edges <- edges[edges[,1] != edges[,2], , drop=FALSE]
  edges <- t(apply(edges, 1, function(x) sort(x)))
  edges <- unique(edges)
  edges_df <- data.frame(from = edges[,1], to = edges[,2])
  
  cat("    边数:", nrow(edges_df), "\n")
  
  if (nrow(edges_df) == 0) {
    cat("    警告：没有边，跳过\n")
    return(NULL)
  }
  
  # 获取每条边的细胞类型
  t1 <- labels_sub[edges_df$from]
  t2 <- labels_sub[edges_df$to]
  
  # 构建观察矩阵
  mat_obs <- matrix(0, nrow = K, ncol = K, dimnames = list(types, types))
  
  for(i in seq_along(t1)) {
    a <- as.character(t1[i])
    b <- as.character(t2[i])
    mat_obs[a, b] <- mat_obs[a, b] + 1
    mat_obs[b, a] <- mat_obs[b, a] + 1
  }
  
  # 排列检验（简化版，减少时间）
  nperm_use <- min(nperm, 200)  # 先跑200次测试
  from_idx <- edges_df$from
  to_idx <- edges_df$to
  n_cells <- length(labels_sub)
  
  ncores <- parallel::detectCores() - 1
  ncores <- max(1, min(ncores, 8))
  
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  
  perm_counts <- foreach(p = 1:nperm_use, .packages = c(), .combine = rbind) %dopar% {
    set.seed(p + 12345)
    perm_labels <- sample(labels_sub, n_cells, replace = FALSE)
    
    pt1 <- perm_labels[from_idx]
    pt2 <- perm_labels[to_idx]
    
    mat <- matrix(0, nrow = K, ncol = K, dimnames = list(types, types))
    for(i in seq_along(pt1)) {
      a <- as.character(pt1[i])
      b <- as.character(pt2[i])
      mat[a, b] <- mat[a, b] + 1
      mat[b, a] <- mat[b, a] + 1
    }
    as.vector(mat)
  }
  
  stopCluster(cl)
  
  # 计算结果
  obs_vec <- as.vector(mat_obs)
  mu_rand <- colMeans(perm_counts)
  sd_rand <- apply(perm_counts, 2, sd)
  
  z_vec <- (obs_vec - mu_rand) / (sd_rand + 1e-8)
  
  p_emp <- sapply(seq_along(obs_vec), function(i) {
    perm_i <- perm_counts[, i]
    obs_i <- obs_vec[i]
    mu_i <- mu_rand[i]
    p_val <- (sum(abs(perm_i - mu_i) >= abs(obs_i - mu_i)) + 1) / (nperm_use + 1)
    return(p_val)
  })
  
  # 转为矩阵
  mat_z <- matrix(z_vec, nrow = K, ncol = K, dimnames = list(types, types))
  mat_p <- matrix(p_emp, nrow = K, ncol = K, dimnames = list(types, types))
  
  # 计算 log2FC
  mat_fc <- mat_obs / (mu_rand + 1e-8)
  mat_log2fc <- log2(mat_fc)
  
  return(list(
    sample = sample_name,
    region = region_name,
    mat_obs = mat_obs,
    mat_z = mat_z,
    mat_p = mat_p,
    mat_log2fc = mat_log2fc,
    cell_types = types,
    n_cells = length(labels_sub),
    n_edges = nrow(edges_df)
  ))
}

# ===================== 按区域进行细胞互作分析 =====================
all_results <- list()

for (sample_id in samples) {
  cat("\n========== 样本:", sample_id, "==========\n")
  
  region_df <- region_list[[sample_id]]
  coords <- data_list[[sample_id]]$coords
  cell_types <- data_list[[sample_id]]$cell_types
  
  # 获取该样本的所有区域类型
  region_types <- unique(region_df$region_type[!is.na(region_df$region_type)])
  cat("  区域类型:", paste(region_types, collapse = ", "), "\n")
  
  for (reg in region_types) {
    cat("\n  区域:", reg, "\n")
    
    # 获取该区域的细胞
    cells_in_region <- region_df$cell[region_df$region_type == reg]
    
    if (length(cells_in_region) < 50) {
      cat("    细胞数太少 (", length(cells_in_region), ")，跳过\n", sep = "")
      next
    }
    
    # 提取坐标和细胞类型
    coords_reg <- coords[cells_in_region, , drop = FALSE]
    cell_types_reg <- cell_types[cells_in_region]
    
    # 过滤稀有细胞类型（数量 < 10）
    type_counts <- table(cell_types_reg)
    keep_types <- names(type_counts[type_counts >= 10])
    keep_cells <- cell_types_reg %in% keep_types
    coords_reg <- coords_reg[keep_cells, , drop = FALSE]
    cell_types_reg <- cell_types_reg[keep_cells]
    
    cat("    过滤后细胞数:", nrow(coords_reg), "\n")
    cat("    细胞类型数:", length(unique(cell_types_reg)), "\n")
    
    if (length(unique(cell_types_reg)) < 2) {
      cat("    过滤后细胞类型太少，跳过\n")
      next
    }
    
    # 运行互作分析
    result <- tryCatch({
      run_spatial_analysis(coords_reg, cell_types_reg, sample_id, reg, nperm = 200)
    }, error = function(e) {
      cat("    错误:", e$message, "\n")
      return(NULL)
    })
    
    if (!is.null(result)) {
      all_results[[sample_id]][[reg]] <- result
      cat("    完成\n")
    }
  }
}

# ===================== 保存结果 =====================
saveRDS(all_results, file = "mPT_region_interaction_results.rds")
cat("\n保存结果: mPT_region_interaction_results.rds\n")

# ===================== 处理有多样本的区域（带显著性）=====================

all_region_types <- unique(unlist(lapply(all_results, function(x) names(x))))
cat("发现区域类型:", paste(all_region_types, collapse = ", "), "\n")

merged_results <- list()

# ===================== 合并多样本区域（修复版）=====================
for (reg_type in all_region_types) {
  cat("\n合并区域:", reg_type, "\n")
  
  # 收集结果
  region_results <- list()
  for (sample_id in names(all_results_detailed)) {
    if (reg_type %in% names(all_results_detailed[[sample_id]])) {
      region_results[[sample_id]] <- all_results_detailed[[sample_id]][[reg_type]]
    }
  }
  
  n_samples <- length(region_results)
  if (n_samples < 2) {
    cat("  跳过（只有1个样本）\n")
    next
  }
  
  # 所有细胞类型（并集）
  all_cell_types <- unique(unlist(lapply(region_results, function(x) x$cell_types)))
  K <- length(all_cell_types)
  cat("  细胞类型数:", K, "\n")
  
  # 初始化矩阵
  avg_log2fc <- matrix(0, nrow = K, ncol = K)
  combined_chi2 <- matrix(0, nrow = K, ncol = K)
  count <- matrix(0, nrow = K, ncol = K)
  rownames(avg_log2fc) <- colnames(avg_log2fc) <- all_cell_types
  rownames(combined_chi2) <- colnames(combined_chi2) <- all_cell_types
  rownames(count) <- colnames(count) <- all_cell_types
  
  # 累加每个样本
  for (sample_id in names(region_results)) {
    res <- region_results[[sample_id]]
    types <- res$cell_types
    
    for (i in seq_along(types)) {
      for (j in seq_along(types)) {
        ci <- types[i]
        cj <- types[j]
        
        val <- res$mat_log2fc[ci, cj]
        pv <- res$mat_p[ci, cj]
        
        if (!is.na(val) && !is.infinite(val)) {
          avg_log2fc[ci, cj] <- avg_log2fc[ci, cj] + val
          combined_chi2[ci, cj] <- combined_chi2[ci, cj] + (-2 * log(pv + 1e-10))
          count[ci, cj] <- count[ci, cj] + 1
        }
      }
    }
  }
  
  # 平均
  for (i in 1:K) {
    for (j in 1:K) {
      if (count[i, j] > 0) {
        avg_log2fc[i, j] <- avg_log2fc[i, j] / count[i, j]
      }
    }
  }
  
  # 合并 p 值
  combined_p <- matrix(1, nrow = K, ncol = K)
  for (i in 1:K) {
    for (j in 1:K) {
      if (count[i, j] > 0) {
        combined_p[i, j] <- pchisq(combined_chi2[i, j], df = 2 * count[i, j], lower.tail = FALSE)
      }
    }
  }
  rownames(combined_p) <- colnames(combined_p) <- all_cell_types
  
  # 显著性标注
  signif <- matrix("", nrow = K, ncol = K)
  rownames(signif) <- colnames(signif) <- all_cell_types
  signif[combined_p < 0.001] <- "***"
  signif[combined_p < 0.01 & combined_p >= 0.001] <- "**"
  signif[combined_p < 0.05 & combined_p >= 0.01] <- "*"
  
  # 绘制热图
  pheatmap(avg_log2fc,
           display_numbers = signif,
           number_color = "black",
           fontsize_number = 6,
           main = paste(reg_type, "region (", n_samples, "samples)"),
           filename = paste0("mPT_region_", gsub("/", "_", reg_type), "_merged.pdf"),
           width = max(10, K * 0.4),
           height = max(8, K * 0.4))
  
  cat("  保存热图: mPT_region_", reg_type, "_merged.pdf\n", sep = "")
}

# ===================== 处理单样本区域 =====================
single_results <- list()

for (reg_type in all_region_types) {
  cat("\n处理区域:", reg_type, "\n")
  
  # 收集该区域类型的所有结果
  region_results <- list()
  
  for (sample_id in names(all_results)) {
    if (reg_type %in% names(all_results[[sample_id]])) {
      region_results[[sample_id]] <- all_results[[sample_id]][[reg_type]]
    }
  }
  
  n_samples <- length(region_results)
  cat("  样本数:", n_samples, "\n")
  
  # 只处理单个样本的区域
  if (n_samples != 1) {
    cat("  跳过（不是单个样本）\n")
    next
  }
  
  sample_id <- names(region_results)[1]
  res <- region_results[[sample_id]]
  
  mat_log2fc <- res$mat_log2fc
  mat_p <- res$mat_p
  
  # 处理异常值
  mat_log2fc[is.infinite(mat_log2fc) & mat_log2fc > 0] <- 10
  mat_log2fc[is.infinite(mat_log2fc) & mat_log2fc < 0] <- -10
  mat_log2fc[is.na(mat_log2fc)] <- 0
  mat_p[is.na(mat_p)] <- 1
  
  # 显著性标注
  signif_symbols <- matrix("", nrow = nrow(mat_p), ncol = ncol(mat_p),
                           dimnames = dimnames(mat_p))
  signif_symbols[mat_p < 0.001] <- "***"
  signif_symbols[mat_p < 0.01 & mat_p >= 0.001] <- "**"
  signif_symbols[mat_p < 0.05 & mat_p >= 0.01] <- "*"
  
  # 绘制热图
  pheatmap(mat_log2fc,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           display_numbers = signif_symbols,
           number_color = "black",
           fontsize_number = 6,
           main = paste("Cell-cell interactions in", reg_type, "region (", sample_id, ")", sep = ""),
           filename = paste0("mPT_region_", gsub("/", "_", reg_type), "_", sample_id, ".pdf"),
           width = 12, height = 10)
  
  cat("  保存热图: mPT_region_", reg_type, "_", sample_id, ".pdf\n", sep = "")
  
  single_results[[reg_type]] <- list(
    mat_log2fc = mat_log2fc,
    mat_p = mat_p,
    sample = sample_id,
    cell_types = res$cell_types
  )
}

cat("\n===== 完成！=====\n")


# ===================== 只保存必要的分析结果 =====================

# 1. 保存区域划分结果（每个细胞的区域标签）
region_assignment <- list()
for (sample_id in names(region_list)) {
  region_df <- region_list[[sample_id]]
  region_assignment[[sample_id]] <- data.frame(
    cell = region_df$cell,
    region = region_df$region,
    region_type = region_df$region_type,
    confidence = region_df$confidence
  )
}
saveRDS(region_assignment, file = "mPT_region_assignment.rds")

# 2. 保存区域互作结果（核心结果）
saveRDS(all_results_detailed, file = "mPT_all_results_detailed.rds")

# 3. 保存合并后的区域互作结果
saveRDS(merged_results, file = "mPT_merged_results.rds")

# 4. 保存单样本区域互作结果
saveRDS(single_results, file = "mPT_single_results.rds")

# 5. 保存区域组成统计（轻量级）
region_composition_stats <- list()
for (sample_id in names(region_list)) {
  region_df <- region_list[[sample_id]]
  region_composition_stats[[sample_id]] <- list(
    region_types = table(region_df$region_type),
    region_confidence = aggregate(confidence ~ region_type, data = region_df, FUN = mean)
  )
}
saveRDS(region_composition_stats, file = "mPT_region_stats.rds")

cat("\n===== mPT 结果已保存（轻量级）=====\n")

```
# 不分样本的区域划分
```R
library(Seurat)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(deldir)

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
                       "Malignant cells" = "unknown",
                       "Basal cells" = "Malignant cells")

# ===================== 提取所有 mPT 细胞（不分样本）=====================
mPT <- subset(obj, subset = tissue == "mPT")
```
# 单细胞辅助注释空间组
```R
#GAE117570
library(Seurat)
library(Matrix)

# 读取数据
counts_matrix3 <- read.table("GSM3304009_P2_Tumor_processed_data.txt.gz", 
                              header = TRUE, 
                              row.names = 1,      # 第一列作为行名（基因名）
                              sep = "\t")         # 根据实际分隔符调整

# 转换为矩阵
counts_matrix3 <- as.matrix(counts_matrix3)

# 确保是数值型
storage.mode(counts_matrix3) <- "numeric"

# 转换为稀疏矩阵（节省内存）
counts_sparse <- as(counts_matrix3, "dgCMatrix")

# 创建Seurat对象
seurat_obj1 <- CreateSeuratObject(counts = counts_sparse)

#GSE127465
library(Seurat)
library(Matrix)

# 定义文件路径
files <- c(
  "GSM3635278_human_p1t1_raw_counts.tsv.gz",
  "GSM3635279_human_p1t2_raw_counts.tsv.gz", 
  "GSM3635280_human_p1t3_raw_counts.tsv.gz",
  "GSM3635281_human_p1t4_raw_counts.tsv.gz",
  "GSM3635285_human_p2t1_raw_counts.tsv.gz",
  "GSM3635286_human_p2t2_raw_counts.tsv.gz"
)

# 提取样本名（去掉前缀和后缀）
sample_names <- gsub("GSM[0-9]+_human_", "", files)  # p1t1, p1t2等
sample_names <- gsub("_raw_counts.tsv.gz", "", sample_names)

# 批量读取并创建Seurat对象列表
seurat_list <- list()

for (i in 1:length(files)) {
  
  # 读取数据
  counts <- read.table(files[i], 
                       header = TRUE, 
                       row.names = 1,
                       sep = "\t",
                       check.names = FALSE)
  
  # 转换为稀疏矩阵
  counts_sparse <- as(as.matrix(counts), "dgCMatrix")
  
  # 创建Seurat对象
  seurat_obj <- CreateSeuratObject(counts = counts_sparse, 
                                    project = sample_names[i],
                                    min.cells = 3,      # 至少在3个细胞表达的基因
                                    min.features = 200) # 至少表达200个基因的细胞
  
  # 添加样本信息
  seurat_obj$sample <- sample_names[i]
  
  # 解析患者和时间点（p1t1 = patient1, time1）
  seurat_obj$patient <- gsub("p[0-9]+", "", sample_names[i])  # 提取p1, p2
  seurat_obj$time <- gsub("t[0-9]+", "", sample_names[i])      # 提取t1, t2等
  
  seurat_list[[i]] <- seurat_obj
  
  cat("Loaded:", sample_names[i], "-", ncol(seurat_obj), "cells\n")
}

# 合并所有样本
seurat_merged <- merge(seurat_list[[1]], 
                       y = seurat_list[-1],
                       add.cell.ids = sample_names)

cat("\nTotal cells:", ncol(seurat_merged))
table(seurat_merged$sample)

#GSE153935

counts_matrix3 <- read.table("GSE153935_TLDS_AllCells.txt.gz", 
                              header = TRUE, 
                              row.names = 1,      # 第一列作为行名（基因名）
                              sep = "\t")         # 根据实际分隔符调整

# 转换为矩阵
counts_matrix3 <- as.matrix(counts_matrix3)

# 确保是数值型
storage.mode(counts_matrix3) <- "numeric"

# 转换为稀疏矩阵（节省内存）
counts_sparse <- as(counts_matrix3, "dgCMatrix")

# 创建Seurat对象
seurat_obj1 <- CreateSeuratObject(counts = counts_sparse)

#GSE148071
library(Seurat)
library(dplyr)
library(data.table)

# 设置工作目录
setwd("E:/project/LUSC-spatial/GSE148071_RAW")

# 获取所有 .txt.gz 文件列表
files <- list.files(pattern = "P\\d+_exp.txt.gz")

# 初始化列表，存储每个样本的 Seurat 对象
seurat_list <- list()

for (file in files) {
  # 提取样本名（如 P1, P2, ...）
  sample_name <- gsub("_exp.txt.gz", "", file)
  
  cat("正在读取", sample_name, "...\n")
  
  # 读取表达矩阵（fread 速度快）
  expr <- fread(file, header = TRUE, data.table = FALSE)
  
  # 第一列通常是基因名，设为行名
  rownames(expr) <- expr[, 1]
  expr <- expr[, -1]
  
  # 转换为矩阵（Seurat 要求）
  expr <- as.matrix(expr)
  
  # 创建 Seurat 对象
  seurat_obj <- CreateSeuratObject(counts = expr, project = sample_name, min.cells = 3, min.features = 200)
  
  # 添加样本信息到 metadata
  seurat_obj$sample <- sample_name
  
  seurat_list[[sample_name]] <- seurat_obj
}

# 合并所有样本（注意：不同样本的基因可能不完全一致，merge 会自动取并集，缺失填0）
combined <- merge(seurat_list[[1]], seurat_list[-1], add.cell.ids = names(seurat_list))

# 从细胞名（行名）中提取 P 数字
combined$sample_id <- gsub(".*_(P\\d+)_.*", "\\1", rownames(combined@meta.data))

# 您要保留的样本列表
selected_samples <- c("P1", "P3", "P4", "P6", "P7", "P10", "P14", "P15", "P17", "P18", "P19", 
                      "P22", "P23", "P25", "P26", "P27", "P30", "P31", "P36", "P37", "P40", "P41")

# 检查哪些样本存在
present <- intersect(selected_samples, unique(combined$sample_id))
cat("存在的样本:", paste(present, collapse = ", "), "\n")

# 提取细胞
cells_keep <- WhichCells(combined, expression = sample_id %in% present)
combined_subset <- subset(combined, cells = cells_keep)

# 查看子集后各样本的细胞数
table(combined_subset$sample_id)
combined_subset <- JoinLayers(combined_subset)
# 保存结果
saveRDS(combined_subset, file = "GSE148071_LUSC.rds")

#sciencedb02028
library(Seurat)

samples <- c("2018jz4", "2018jz11", "2018jz12",
             "2018jz38", "2018jz46", "AK658", "AK2834", "AK3274", "AK4295")

target_base <- "/share/home/wangq/zxy/NSCLC/spatial/sc/outs"

# 读取第一个样本（使用目录）
seurat_obj <- CreateSeuratObject(
  counts = Read10X(file.path(target_base, samples[1])),
  project = samples[1]
)

# 依次合并其余样本
for (i in 2:length(samples)) {
  cat("Processing:", samples[i], "\n")
  seurat_obj <- merge(
    seurat_obj,
    CreateSeuratObject(
      counts = Read10X(file.path(target_base, samples[i])),
      project = samples[i]
    )
  )
}
seurat_obj <- JoinLayers(seurat_obj)
# 保存
saveRDS(seurat_obj, file = file.path(target_base, "sciencedb.rds"))

#emtab
library(Seurat)
library(dplyr)

# 设置路径
path <- "./EMTAB"

# 获取所有样本的ID（Pxx_Tx）
files <- list.files(path, pattern = "*.tsv.gz")
sample_ids <- unique(gsub("-(barcodes|features|matrix).tsv.gz", "", files))

sample_files <- list()
for (sid in sample_ids) {
  sample_files[[sid]] <- list(
    barcodes = file.path(path, paste0(sid, "-barcodes.tsv.gz")),
    features = file.path(path, paste0(sid, "-features.tsv.gz")),
    matrix = file.path(path, paste0(sid, "-matrix.mtx.gz"))
  )
}

# 读取所有样本
seurat_list <- list()
for (sid in names(sample_files)) {
  cat("Loading:", sid, "\n")
  
  # 读取矩阵
  counts <- ReadMtx(
    mtx = sample_files[[sid]]$matrix,
    features = sample_files[[sid]]$features,
    cells = sample_files[[sid]]$barcodes
  )
  
  # 创建 Seurat 对象
  obj <- CreateSeuratObject(counts = counts, project = sid)
  obj$sample <- sid
  obj$patient <- gsub("_T.*", "", sid)  # 提取患者ID
  obj$region <- gsub("P[0-9]+_", "", sid)  # 提取区域 T1/T2/T3
  
  seurat_list[[sid]] <- obj
}

seurat_merged <- merge(seurat_list[[1]], y = seurat_list[-1], 
                       add.cell.ids = names(seurat_list))

saveRDS(seurat_merged,file="emtab.rds")


obj <- readRDS("GSE207422.rds")
obj <- subset(obj, subset = cancer_type == "LUSC")

#obj2 <- readRDS("GSE148071_LUSC.rds")
obj3 <- readRDS("sciencedb_LN.rds")
#obj3 <- subset(obj3,subset=orig.ident %in% c("AK3274","AK4295","AK2834","2018jz4","2018jz37"))
#obj3$LN_status <- ifelse(obj3$orig.ident == "2018jz37", "LN+", "LN-")
#saveRDS(obj3,file="sciencedb_LN.rds")
obj4 <- readRDS("emtab_LN.rds")
#obj4 <- subset(obj4,subset=patient %in% c("P2","P19","P4","P8","P11","P18"))
#obj4$LN_status <- ifelse(obj4$patient %in% c("P2","P8","P18"), "LN+", "LN-")
#saveRDS(obj4,file="emtab_LN.rds")
merged_seurat_obj <- merge(obj, y = list(obj3, obj4), 
                    add.cell.ids = c("GSE207422", "sciencedb", "emtab"))
merged_seurat_obj <- JoinLayers(merged_seurat_obj)

merged_seurat_obj <- merged_seurat_obj[!grepl("\\.", rownames(merged_seurat_obj)), ]
merged_seurat_obj <- merged_seurat_obj[!grepl("^MT-", rownames(merged_seurat_obj)), ]
merged_seurat_obj <- merged_seurat_obj[!grepl("^RP[SL]", rownames(merged_seurat_obj)), ]
merged_seurat_obj <- merged_seurat_obj[!grepl("AP[0-9]+$", rownames(merged_seurat_obj)), ]
nrow(merged_seurat_obj)

obj <- merged_seurat_obj
obj <- obj[!grepl("^ENSG", rownames(obj)), ]

library(dplyr)
obj$dataset <- case_when(
  obj$orig.ident == "BD" ~ "BD",
  grepl("^P", obj$orig.ident) ~ "emtab",
  TRUE ~ "sciencedb"
)
saveRDS(obj,file="sc_raw.rds")

library(harmony)
library(Seurat)
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, nfeatures = 2000)

hvgs <- VariableFeatures(obj)
obj <- ScaleData(obj, features = hvgs)
obj <- RunPCA(obj, features = hvgs, npcs = 50)


obj <- RunHarmony(obj, "dataset")

png("elbowplot_sc.png", width = 800, height = 600)
elbowplot <- ElbowPlot(obj,ndims = 30)
print(elbowplot)
dev.off()

library(ggplot2)
library(ggsci)
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:15)
obj <- FindClusters(obj, resolution = 0.2)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:15)

npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(20)
seurat_clusters <- as.character(unique(obj@meta.data$seurat_clusters))  # 转换为字符向量
num_legend_items <- length(seurat_clusters)  # 图例的个数
max_label_length <- max(nchar(seurat_clusters))  # 图例名称的最大长度

base_width <- 3000  # 基础宽度
base_height <- 3000  # 基础高度
legend_width_factor <- 100  # 每个图例项增加的宽度
label_length_factor <- 10  # 每个字符增加的宽度
dynamic_width <- base_width + (num_legend_items * legend_width_factor) + (max_label_length * label_length_factor)
pdf("sc_clusters.pdf", width = dynamic_width/300, height = base_height/300)  # 转换为英寸

DimPlot(obj, reduction = "umap", label = TRUE, pt.size = 1, label.size = 8) +
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    ggtitle(NULL) +
    scale_color_manual(values = npg_extended) +
    coord_fixed(ratio = 1) +
    guides(color = guide_legend(title = NULL, override.aes = list(size = 5))) +
    theme(
        text = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 28, color = "black"),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title.x = element_text(size = 36, face = "bold", color = "black"),
        axis.title.y = element_text(size = 36, face = "bold", color = "black", margin = margin(r = 20)),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 28, face = "bold", color = "black"),
        legend.title = element_text(size = 28, face = "bold", color = "black"),
        legend.position = "right",
        legend.box.margin = margin(0, 0, 0, 0),
        legend.key = element_blank(),
        legend.background = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5),
        aspect.ratio = 1,
        plot.margin = margin(10, 50, 10, 10)
    )

dev.off()

markers <- FindAllMarkers(obj, 
                          only.pos = TRUE, 
                          min.pct = 0.25, 
                          logfc.threshold = 0.25, 
                          test.use = "wilcox")
#significant markers
library(dplyr)
significant_markers <- subset(markers, p_val_adj < 0.05)
significant_markers <- significant_markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
write.csv(significant_markers,"sc_obj_marker_top_50.csv")
identity_mapping <- c(
    "0" = "T cells",
    "1" = "T cells",
    "2" = "B cells",
    "3" = "Macrophages",
    "4" = "Macrophages",
    "5" = "Macrophages",
    "6" = "Malignant cells",
    "7" = "DC",
    "8" = "Plasma",
    "9" = "Neutrophil",
    "10" = "Proliferative",
    "11" = "Fibroblasts",
    "12" = "Malignant cells",  
    "13" = "Mast cells",
    "14" = "T cells",
    "15" = "Ciliated cell",
    "16" = "Malignant cells", 
    "17" = "Alveolar type II cells",
    "18" = "T cells",
    "19" = "Macrophages"
)

cell_type <- identity_mapping[obj@meta.data$seurat_clusters]
obj@meta.data$cell_type <- cell_type

library(ggplot2)
library(ggsci)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(17)
cluster_colors <- c(
  # 浅色系在前（浅蓝第一个）
  "#A0CBE8",   # 1 浅蓝（你要的放第一个）
  "#FFBE7D",   # 2 浅橙
  "#8CD17D",   # 3 浅绿
  "#86BCB6",   # 4 浅青
  "#FF9D9A",   # 5 浅红
  "#FABFD2",   # 6 浅粉
  "#D4A6C8",   # 7 浅紫
  "#F1CE63",   # 8 淡黄
  "#D7B5A6",   # 9 浅棕
  "#B07AA1",   # 10 紫灰
  "#BAB0AC",   # 11 浅灰
  
  # 深色系在后
  "#4E79A7",   # 12 深蓝
  "#F28E2B",   # 13 橙
  "#59A14F",   # 14 深绿
  "#E15759",   # 15 红
  "#499894",   # 16 青绿
  "#D37295",   # 17 玫红
  "#B6992D",   # 18 土黄
  "#9D7660",   # 19 棕
  "#79706E"    # 20 深灰
)

cell_types <- as.character(unique(obj@meta.data$cell_type))
num_legend_items <- length(cell_types)  # 图例的个数
max_label_length <- max(nchar(cell_types))  # 图例名称的最大长度

base_width <- 3000  # 基础宽度
base_height <- 3000  # 基础高度
legend_width_factor <- 100  # 每个图例项增加的宽度
label_length_factor <- 10  # 每个字符增加的宽度
dynamic_width <- base_width + (num_legend_items * legend_width_factor) + (max_label_length * label_length_factor)
pdf("sc_obj_annotation.pdf", width = dynamic_width/300, height = base_height/300)
DimPlot(obj, reduction = "umap", label = TRUE, pt.size = 1, group.by = "cell_type", label.size = 4) +
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    ggtitle(NULL) +
    scale_color_manual(values = cluster_colors) +
    coord_fixed(ratio = 1) +
    guides(color = guide_legend(title = NULL, override.aes = list(size = 5))) +
    theme(
        text = element_text(size = 8, face = "bold"),
        axis.text.x = element_text(size = 28, color = "black"),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title.x = element_text(size = 36, face = "bold", color = "black"),
        axis.title.y = element_text(size = 36, face = "bold", color = "black", margin = margin(r = 20)),  # 增加右侧间距
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 16, face = "bold", color = "black"),
        legend.title = element_text(size = 16, face = "bold", color = "black"),
        legend.position = "right",
        legend.box.margin = margin(0, 0, 0, 0),
        legend.key = element_blank(),
        legend.background = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5),
        aspect.ratio = 1,
        plot.margin = margin(10, 50, 10, 10)
    )
dev.off()


saveRDS(obj,file="sc_obj_anno.rds")
obj_sc <- obj
library(pheatmap)
library(ggplot2)

# 1. 获取所有细胞类型（单细胞和空间，使用大类注释）
sc_types <- unique(obj_sc$cell_type)          # 单细胞注释列：cell_type
sp_types <- unique(obj$CellType)              # 空间注释列：CellType
cat("单细胞类型数:", length(sc_types), "\n")
cat("空间类型数:", length(sp_types), "\n")

# 2. 使用重叠基因（需事先计算 overlap.genes，或直接用共有基因）
# 假设 overlap.genes 已经存在，否则可以重新计算：
overlap.genes <- intersect(rownames(obj_sc), rownames(obj))

expr_sc <- GetAssayData(obj_sc, assay = "RNA", layer = "data")[overlap.genes, ]
expr_sp <- GetAssayData(obj, assay = "SCT", layer = "data")[overlap.genes, ]

min_cells <- 5

# 单细胞平均表达谱
sc_avg <- matrix(NA, nrow = length(sc_types), ncol = length(overlap.genes))
rownames(sc_avg) <- sc_types
colnames(sc_avg) <- overlap.genes
for (i in seq_along(sc_types)) {
  ct <- sc_types[i]
  cells <- colnames(obj_sc)[obj_sc$cell_type == ct]
  if (length(cells) >= min_cells) {
    sc_avg[i, ] <- rowMeans(expr_sc[, cells, drop = FALSE])
  } else {
    cat("单细胞类型", ct, "细胞数不足", min_cells, "，跳过\n")
  }
}

# 空间平均表达谱
sp_avg <- matrix(NA, nrow = length(sp_types), ncol = length(overlap.genes))
rownames(sp_avg) <- sp_types
colnames(sp_avg) <- overlap.genes
for (i in seq_along(sp_types)) {
  ct <- sp_types[i]
  cells <- colnames(obj)[obj$CellType == ct]
  if (length(cells) >= min_cells) {
    sp_avg[i, ] <- rowMeans(expr_sp[, cells, drop = FALSE])
  } else {
    cat("空间类型", ct, "细胞数不足", min_cells, "，跳过\n")
  }
}

# 移除全NA的行
sc_avg <- sc_avg[complete.cases(sc_avg), , drop = FALSE]
sp_avg <- sp_avg[complete.cases(sp_avg), , drop = FALSE]
cat("有效单细胞类型数:", nrow(sc_avg), "\n")
cat("有效空间类型数:", nrow(sp_avg), "\n")

cor_matrix <- matrix(NA, nrow = nrow(sc_avg), ncol = nrow(sp_avg),
                     dimnames = list(rownames(sc_avg), rownames(sp_avg)))
for (i in 1:nrow(sc_avg)) {
  for (j in 1:nrow(sp_avg)) {
    cor_matrix[i, j] <- cor(sc_avg[i, ], sp_avg[j, ], method = "spearman", use = "complete.obs")
  }
}

# 5. 绘制热图（不限制颜色范围）
pheatmap(cor_matrix,
         main = "Spearman correlation of average gene expression\nbetween scRNA-seq and spatial (major cell types)",
         fontsize = 8,
         filename = "major_celltype_correlation.pdf",
         width = 12, height = 10,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         cluster_rows = TRUE,
         cluster_cols = TRUE)




```
# Malignant sc
```R
library(Seurat)
obj_sc <- readRDS("sc_obj_anno.rds")
Malignant_sc <- subset(obj_sc,subset=cell_type=="Malignant cells")
Malignant_sc <- subset(Malignant_sc,subset=orig.ident %in% c("2018jz4","AK2834","AK3274","BD","P18_T2","P19_T2","P4_T2","P4_T3","P8_T2"))

Malignant_sc <- NormalizeData(Malignant_sc)
Malignant_sc <- FindVariableFeatures(Malignant_sc, nfeatures = 3000)
hvgs <- VariableFeatures(Malignant_sc)
Malignant_sc <- ScaleData(Malignant_sc, features = hvgs)
Malignant_sc <- RunPCA(Malignant_sc, features = hvgs, npcs = 50)
p <- ElbowPlot(Malignant_sc, ndims = 30)
ggsave("p3.png",plot=p)
# 计算累计方差比例
pca_var <- Malignant_sc[["pca"]]@stdev^2 / sum(Malignant_sc[["pca"]]@stdev^2)
cumsum_var <- cumsum(pca_var)

# 看达到多少 PC 时累计方差 > 70% / 80% / 90%
which(cumsum_var > 0.7)[1]   # 70%
which(cumsum_var > 0.8)[1]   # 80%
which(cumsum_var > 0.9)[1]   # 90%

library(harmony)
Malignant_sc <- RunHarmony(Malignant_sc, "dataset")
Malignant_sc <- FindNeighbors(Malignant_sc, reduction = "harmony", dims = 1:17)
Malignant_sc <- FindClusters(Malignant_sc, resolution = 0.3)
Malignant_sc <- RunUMAP(Malignant_sc, reduction = "harmony", dims = 1:17)


library(ggplot2)
library(ggsci)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(23)#0.6

seurat_clusters <- as.character(unique(Malignant_sc@meta.data$seurat_clusters))  # 转换为字符向量
num_legend_items <- length(seurat_clusters)  # 图例的个数
max_label_length <- max(nchar(seurat_clusters))  # 图例名称的最大长度

base_width <- 3000  # 基础宽度
base_height <- 3000  # 基础高度
legend_width_factor <- 100  # 每个图例项增加的宽度
label_length_factor <- 10  # 每个字符增加的宽度
dynamic_width <- base_width + (num_legend_items * legend_width_factor) + (max_label_length * label_length_factor)
pdf("Malignant_sc_clusters.pdf", width = dynamic_width/300, height = base_height/300)  # 转换为英寸

DimPlot(Malignant_sc, reduction = "umap", label = TRUE, pt.size = 1, label.size = 8) +
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    ggtitle(NULL) +
    scale_color_manual(values = npg_extended) +
    coord_fixed(ratio = 1) +
    guides(color = guide_legend(title = NULL, override.aes = list(size = 5))) +
    theme(
        text = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 28, color = "black"),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title.x = element_text(size = 36, face = "bold", color = "black"),
        axis.title.y = element_text(size = 36, face = "bold", color = "black", margin = margin(r = 20)),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 28, face = "bold", color = "black"),
        legend.title = element_text(size = 28, face = "bold", color = "black"),
        legend.position = "right",
        legend.box.margin = margin(0, 0, 0, 0),
        legend.key = element_blank(),
        legend.background = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5),
        aspect.ratio = 1,
        plot.margin = margin(10, 50, 10, 10)
    )

dev.off()

p <- DimPlot(Malignant_sc, reduction = "umap", group.by = "orig.ident", label = FALSE)
ggsave("malignant_orig.png",plot=p)

library(dplyr)
Malignant_markers <- FindAllMarkers(Malignant_sc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "wilcox")
Malignant_significant_markers <- subset(Malignant_markers, p_val_adj < 0.05)
#write.csv(Malignant_significant_markers, "Malignant_all_marker.csv")
Malignant_significant_markers <- Malignant_significant_markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
write.csv(Malignant_significant_markers, "Malignant_sc_top_marker_50.csv")
#Malignant_sc <- subset(Malignant_sc,subset=seurat_clusters %in% c(0,1,2,3,4,5,6,7,8,9,10,11))

#0.6
identity_mapping <- c(
  "0" = "Tumor_Metabolic",      # 代谢型
  "1" = "Tumor_Hypoxia",       # 缺氧型（最接近炎症相关）
  "2" = "Tumor_Proliferative_1",  # S期增殖
  "3" = "Tumor_Squamous",       # 鳞状分化（最明确）
  "4" = "Tumor_Proliferative_2",  # M期增殖
  "5" = "Tumor_Inflamed_1",         # 间质/EMT
  "6" = "Tumor_Inflamed_2"        # 免疫调节型
)

identity_mapping <- c(
  "0" = "Tumor_Metabolic",      # 代谢型
  "1" = "Tumor_Proliferative_1",  # S期增殖
  "2" = "Tumor_Hypoxia",       # 缺氧型（最接近炎症相关）
  "3" = "Tumor_Mesenchymal-like",
  "4" = "Tumor_Squamous",       # 鳞状分化（最明确）
  "5" = "Tumor_Proliferative_2",  # M期增殖
  "6" = "Tumor_Inflamed_1",         # 间质/EMT
  "7" = "Tumor_Inflamed_2",        # 免疫调节型
  "8" = "Tumor_EMT-like"
)

identity_mapping <- c(
  "0" = "Tumor_Homeostatic-metabolic",
  "1" = "Tumor_Differentiated",
  "2" = "Tumor_Inflamed_1",
  "3" = "Tumor_Proliferative",
  "4" = "Tumor_Immune-inflamed_1",
  "5" = "Tumor_Immune-inflamed_2",
  "6" = "Tumor_Immune-inflamed_3",
  "7" = "Tumor_Inflamed_2",
  "8" = "Tumor_Secretory",
  "9" = "Tumor_Matrix-remodeling",
  "10" = "Tumor_EMT-like"
)
sub_cell_type <- identity_mapping[Malignant_sc@meta.data$seurat_clusters]
Malignant_sc@meta.data$sub_cell_type <- sub_cell_type

library(ggplot2)
library(ggsci)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(9)
cluster_colors <- c(
  # 浅色系在前（浅蓝第一个）
  "#A0CBE8",   # 1 浅蓝（你要的放第一个）
  "#FFBE7D",   # 2 浅橙
  "#8CD17D",   # 3 浅绿
  "#86BCB6",   # 4 浅青
  "#FF9D9A",   # 5 浅红
  "#FABFD2",   # 6 浅粉
  "#D4A6C8",   # 7 浅紫
  "#F1CE63",   # 8 淡黄
  "#D7B5A6",   # 9 浅棕
  "#B07AA1",   # 10 紫灰
  "#BAB0AC",   # 11 浅灰
  
  # 深色系在后
  "#4E79A7",   # 12 深蓝
  "#F28E2B",   # 13 橙
  "#59A14F",   # 14 深绿
  "#E15759",   # 15 红
  "#499894",   # 16 青绿
  "#D37295",   # 17 玫红
  "#B6992D",   # 18 土黄
  "#9D7660",   # 19 棕
  "#79706E"    # 20 深灰
)
cell_types <- as.character(unique(Malignant_sc@meta.data$sub_cell_type))
num_legend_items <- length(cell_types)  # 图例的个数
max_label_length <- max(nchar(cell_types))  # 图例名称的最大长度

base_width <- 3000  # 基础宽度
base_height <- 3000  # 基础高度
legend_width_factor <- 100  # 每个图例项增加的宽度
label_length_factor <- 10  # 每个字符增加的宽度
dynamic_width <- base_width + (num_legend_items * legend_width_factor) + (max_label_length * label_length_factor)
pdf("malignant_sc_annotation.pdf", width = dynamic_width/300, height = base_height/300)
DimPlot(Malignant_sc, reduction = "umap", label = TRUE, pt.size = 1, group.by = "sub_cell_type", label.size = 4) +
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    ggtitle(NULL) +
    scale_color_manual(values =cluster_colors) +
    coord_fixed(ratio = 1) +
    guides(color = guide_legend(title = NULL, override.aes = list(size = 5))) +
    theme(
        text = element_text(size = 8, face = "bold"),
        axis.text.x = element_text(size = 28, color = "black"),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title.x = element_text(size = 36, face = "bold", color = "black"),
        axis.title.y = element_text(size = 36, face = "bold", color = "black", margin = margin(r = 20)),  # 增加右侧间距
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 16, face = "bold", color = "black"),
        legend.title = element_text(size = 16, face = "bold", color = "black"),
        legend.position = "right",
        legend.box.margin = margin(0, 0, 0, 0),
        legend.key = element_blank(),
        legend.background = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5),
        aspect.ratio = 1,
        plot.margin = margin(10, 50, 10, 10)
    )
dev.off()
saveRDS(Malignant_sc,file="malignant_sc_anno.rds")

#比例
library(tidyverse)

# 计算各亚群在不同 LN_status 中的占比（过滤掉 NA）
prop_data <- Malignant_sc@meta.data %>%
  filter(!is.na(LN_status)) %>%
  group_by(LN_status, sub_cell_type) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(LN_status) %>%
  mutate(proportion = count / sum(count) * 100)

library(ggplot2)
library(ggsci)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(14)

prop_data <- prop_data %>%
  mutate(sub_cell_type = as.character(sub_cell_type))

# 绘制堆叠条形图
p <- ggplot(prop_data, aes(x = LN_status, y = proportion, fill = sub_cell_type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = cluster_colors) +
  labs(title = "Cell Type Distribution: LN+ vs LN-",
       x = "LN Status", y = "Proportion (%)", fill = "Cell Type") +
  theme_minimal() +
  theme(
    # 坐标轴线条和刻度
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.15, "cm"),
    # 坐标轴文本
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 10),
    axis.text.y = element_text(color = "black", size = 10),
    # 坐标轴标题
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    # 面板背景
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    # 图例
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )

# 保存为PDF
pdf("Malignant_sc_celltype_distribution_LN_status.pdf", width = 4, height = 6)
print(p)
dev.off()

library(tidyverse)
library(ggplot2)

# ==================== 1. 计算各样本中各亚群的占比 ====================
prop_data <- Malignant_sc@meta.data %>%
  filter(!is.na(LN_status)) %>%
  group_by(LN_status, orig.ident, sub_cell_type) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(LN_status, orig.ident) %>%
  mutate(proportion = count / sum(count) * 100)

# ==================== 2. 按 LN_status 排序：LN+ 在前，LN- 在后 ====================
samples_LNpos <- prop_data %>% 
  filter(LN_status == "LN+") %>% 
  distinct(orig.ident) %>% 
  pull(orig.ident) %>% 
  sort()

samples_LNneg <- prop_data %>% 
  filter(LN_status == "LN-") %>% 
  distinct(orig.ident) %>% 
  pull(orig.ident) %>% 
  sort()

sample_order <- c(samples_LNpos, samples_LNneg)
prop_data$orig.ident <- factor(prop_data$orig.ident, levels = sample_order)

# ==================== 3. 处理亚型顺序 ====================
prop_data$sub_cell_type <- as.character(prop_data$sub_cell_type)
type_order <- sort(unique(prop_data$sub_cell_type))
prop_data$sub_cell_type <- factor(prop_data$sub_cell_type, levels = type_order)

# ==================== 4. 绘制水平堆叠条形图 ====================
p <- ggplot(prop_data, aes(x = proportion, y = orig.ident, fill = sub_cell_type)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  scale_fill_manual(values = cluster_colors, name = "Cell Type") +
  labs(title = "Cell Type Proportion by Sample (LN+ / LN-)",
       x = "Proportion (%)", y = "Sample") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 12, face = "bold"),
    panel.grid.major.x = element_line(color = "gray90"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 6),
    legend.key.size = unit(0.4, "cm")
  ) +
  guides(fill = guide_legend(nrow = 6, byrow = TRUE))   # 放在 theme() 外面

# ==================== 5. 添加 LN+ / LN- 分界线和标签 ====================
n_LNpos <- length(samples_LNpos)
if (n_LNpos > 0) {
  p <- p + geom_hline(yintercept = n_LNpos + 0.5, linetype = "dashed", color = "red", linewidth = 0.5)
}

# 添加分组标签（在 y 轴右侧或图内）
p <- p + annotate("text", x = 95, y = n_LNpos/2 + 0.5, label = "LN+", size = 4, fontface = "bold", hjust = 1) +
         annotate("text", x = 95, y = n_LNpos + (length(samples_LNneg)/2) + 0.5, label = "LN-", size = 4, fontface = "bold", hjust = 1)

# ==================== 6. 保存 ====================
pdf("Malignant_sc_celltype_proportion_by_sample_LN_status_horizontal.pdf", width = 4, height = 8)
print(p)
dev.off()

```
# malignant mapping
```R
library(Seurat)
library(dplyr)

Malignant_sc <- readRDS("malignant_sc_anno.rds")
Malignant <- readRDS("malignant_anno.rds")

# ===================== 【1】内存 + 并行关闭 =====================
options(future.globals.maxSize = 100 * 1024^3)  # 100G
future::plan("sequential")                     # 关闭并行，最关键！

# ===================== 【2】强制清理 =====================
DefaultAssay(Malignant) <- "RNA"
DefaultAssay(Malignant_sc) <- "RNA"
Malignant@assays$SCT <- NULL
Malignant_sc@assays$SCT <- NULL

# ===================== 【3】合并图层（唯一能解决你问题的步骤） =====================
Malignant <- JoinLayers(Malignant)       # 必加！
Malignant_sc <- JoinLayers(Malignant_sc) # 必加！

# ===================== 【4】重新完整跑一遍流程 =====================
# 标准化
Malignant <- NormalizeData(Malignant)
Malignant_sc <- NormalizeData(Malignant_sc)

# 高变基因
Malignant <- FindVariableFeatures(Malignant, nfeatures = 3000)
Malignant_sc <- FindVariableFeatures(Malignant_sc, nfeatures = 3000)

# 共同基因
overlap.genes <- intersect(rownames(Malignant), rownames(Malignant_sc))

# 缩放！！！（必须在 JoinLayers 之后）
Malignant <- ScaleData(Malignant, features = overlap.genes)
Malignant_sc <- ScaleData(Malignant_sc, features = overlap.genes)

# PCA
Malignant <- RunPCA(Malignant, features = overlap.genes, npcs = 50)
Malignant_sc <- RunPCA(Malignant_sc, features = overlap.genes, npcs = 50)

# ===================== 【5】锚点 + 注释 =====================
transfer.anchors <- FindTransferAnchors(
  reference = Malignant_sc,
  query = Malignant,
  dims = 1:15,                # 从20 → 15，大幅降内存！
  features = overlap.genes,
  
  # 三个内存杀手参数，我帮你全关！
  reference.assay = "RNA",
  query.assay = "RNA",
  normalization.method = "LogNormalize",
  
  # 【爆内存救星】
  k.anchor = 5,               # 越小占内存越少
  approx.pca = TRUE           # 快速+低内存
)

predictions <- TransferData(
  anchorset = transfer.anchors,
  refdata = Malignant_sc$sub_cell_type,
  dims = 1:15,
  k.weight = 10
)
Malignant$cell_type_pred <-NA
Malignant <- AddMetaData(Malignant, predictions)
Malignant$cell_type_pred <- Malignant$predicted.id
table(Malignant@meta.data$cell_type_pred,Malignant@meta.data$sub_cell_type)

library(ggplot2)
library(ggsci)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(6)

cell_types <- as.character(unique(Malignant@meta.data$cell_type_pred))
num_legend_items <- length(cell_types)  # 图例的个数
max_label_length <- max(nchar(cell_types))  # 图例名称的最大长度

base_width <- 3000  # 基础宽度
base_height <- 3000  # 基础高度
legend_width_factor <- 100  # 每个图例项增加的宽度
label_length_factor <- 10  # 每个字符增加的宽度
dynamic_width <- base_width + (num_legend_items * legend_width_factor) + (max_label_length * label_length_factor)
pdf("malignant_predict_annotation.pdf", width = dynamic_width/300, height = base_height/300)
DimPlot(Malignant, reduction = "umap", label = TRUE, pt.size = 1, group.by = "cell_type_pred", label.size = 4) +
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    ggtitle(NULL) +
    scale_color_manual(values = npg_pal) +
    coord_fixed(ratio = 1) +
    guides(color = guide_legend(title = NULL, override.aes = list(size = 5))) +
    theme(
        text = element_text(size = 8, face = "bold"),
        axis.text.x = element_text(size = 28, color = "black"),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title.x = element_text(size = 36, face = "bold", color = "black"),
        axis.title.y = element_text(size = 36, face = "bold", color = "black", margin = margin(r = 20)),  # 增加右侧间距
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 16, face = "bold", color = "black"),
        legend.title = element_text(size = 16, face = "bold", color = "black"),
        legend.position = "right",
        legend.box.margin = margin(0, 0, 0, 0),
        legend.key = element_blank(),
        legend.background = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5),
        aspect.ratio = 1,
        plot.margin = margin(10, 50, 10, 10)
    )
dev.off()

#dotplot
all_genes <- c(
  # ===== E: EMT-like (Cluster 8) =====
  "SPARCL1", "COL4A2", "FASN",
  
  # ===== I: Immune-inflamed (Cluster 1, 3) =====
  "IGKC", "IGLL5", "APOE", 
  
  # ===== I: Immunogenic (Cluster 2, 5) =====
   "TAPBP", "PSMB9",
  
  # ===== M: Metabolic_1 (Cluster 0) =====
  "AKR1C1", "NQO1", "CPT1A", "TXNRD1",
  
  # ===== M: Metabolic_2 (Cluster 9) =====
  "ALDH3A1", 
  
  # ===== P: Proliferative (Cluster 4, 7) =====
 "MCM3", "HMGCS1", "LAMC2", 
  
  # ===== S: Secretory (Cluster 6) =====
  "SCGB1A1", "SCGB3A1", "MUC4", 
  # ===== S: Stem-like (Cluster 4) =====
  "SOX2", "PTGS2", "MCM3"
)
all_genes <- unique(all_genes)  # 去重（如 CD74, G6PD, TOP2A, KRT16 等重复）

library(Seurat)
library(ggplot2)
library(patchwork)

# 修复因子水平
Malignant$sub_cell_type <- as.character(Malignant$sub_cell_type)
Malignant_sc$sub_cell_type <- as.character(Malignant_sc$sub_cell_type)

# 空间组气泡图
p1 <- DotPlot(Malignant, features = all_genes, group.by = "sub_cell_type", scale = TRUE) +
  scale_color_gradient(low = "#D6EAF8", high = "#1F618D", name = "Average Expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
  labs(title = "Spatial - Tumor Subtypes", x = "", y = "") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# 单细胞气泡图
p2 <- DotPlot(Malignant_sc, features = all_genes, group.by = "sub_cell_type", scale = TRUE) +
  scale_color_gradient(low = "#D6EAF8", high = "#1F618D", name = "Average Expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
  labs(title = "scRNA-seq - Tumor Subtypes", x = "", y = "") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# 拼图
p_combined <- p1 / p2 + 
  plot_annotation(title = "Tumor Marker Genes: Spatial vs scRNA-seq") &
  theme(plot.title = element_text(hjust = 0.5, size = 14))

ggsave("bubble_plot_tumor_final.pdf", plot = p_combined, width = 18, height = 14)

#heatmap
library(pheatmap)
library(ggplot2)
library(dplyr)

# ==================== 1. 获取单细胞 Malignant_sc 的 top 200 基因 ====================
cat("正在计算 Malignant_sc 的标记基因...\n")
Malignant_sc_markers <- FindAllMarkers(Malignant_sc, 
                                        only.pos = TRUE, 
                                        min.pct = 0.25, 
                                        logfc.threshold = 0.25, 
                                        test.use = "wilcox")
Malignant_sc_significant <- subset(Malignant_sc_markers, p_val_adj < 0.05)
Malignant_sc_top <- Malignant_sc_significant %>% 
                    group_by(cluster) %>% 
                    top_n(n = 200, wt = avg_log2FC)
sc_genes <- unique(Malignant_sc_top$gene)
cat("单细胞 top 200 基因数:", length(sc_genes), "\n")

# ==================== 2. 获取空间 Malignant 的 top 200 基因 ====================
cat("正在计算 Malignant 的标记基因...\n")
Malignant_markers <- FindAllMarkers(Malignant, 
                                     only.pos = TRUE, 
                                     min.pct = 0.25, 
                                     logfc.threshold = 0.25, 
                                     test.use = "wilcox")
Malignant_significant <- subset(Malignant_markers, p_val_adj < 0.05)
Malignant_top <- Malignant_significant %>% 
                 group_by(cluster) %>% 
                 top_n(n = 200, wt = avg_log2FC)
sp_genes <- unique(Malignant_top$gene)
cat("空间 top 200 基因数:", length(sp_genes), "\n")

# ==================== 3. 计算重叠基因 ====================
overlap.genes <- intersect(sc_genes, sp_genes)
cat("重叠基因数:", length(overlap.genes), "\n")

# 如果重叠基因太少，给出警告
if (length(overlap.genes) < 50) {
  cat("警告：重叠基因较少，建议增大 top_n 或使用高变基因\n")
}

# ==================== 4. 获取细胞类型 ====================
sc_types <- unique(Malignant_sc$sub_cell_type)
sp_types <- unique(Malignant$sub_cell_type)
cat("单细胞类型数:", length(sc_types), "\n")
cat("空间类型数:", length(sp_types), "\n")

# ==================== 5. 获取表达数据 ====================
# 使用 counts layer 并做 log1p 转换
expr_sc <- GetAssayData(Malignant_sc, assay = "RNA", layer = "counts")[overlap.genes, , drop = FALSE]
expr_sp <- GetAssayData(Malignant, assay = "RNA", layer = "counts")[overlap.genes, , drop = FALSE]

expr_sc <- log1p(expr_sc)
expr_sp <- log1p(expr_sp)

min_cells <- 5

# ==================== 6. 计算平均表达谱 ====================
# 单细胞
sc_avg <- matrix(NA, nrow = length(sc_types), ncol = length(overlap.genes))
rownames(sc_avg) <- sc_types
colnames(sc_avg) <- overlap.genes
for (i in seq_along(sc_types)) {
  ct <- sc_types[i]
  cells <- colnames(Malignant_sc)[Malignant_sc$sub_cell_type == ct]
  if (length(cells) >= min_cells) {
    sc_avg[i, ] <- rowMeans(expr_sc[, cells, drop = FALSE])
  } else {
    cat("单细胞类型", ct, "细胞数不足", min_cells, "，跳过\n")
  }
}

# 空间
sp_avg <- matrix(NA, nrow = length(sp_types), ncol = length(overlap.genes))
rownames(sp_avg) <- sp_types
colnames(sp_avg) <- overlap.genes
for (i in seq_along(sp_types)) {
  ct <- sp_types[i]
  cells <- colnames(Malignant)[Malignant$sub_cell_type == ct]
  if (length(cells) >= min_cells) {
    sp_avg[i, ] <- rowMeans(expr_sp[, cells, drop = FALSE])
  } else {
    cat("空间类型", ct, "细胞数不足", min_cells, "，跳过\n")
  }
}

# 移除全NA的行
sc_avg <- sc_avg[complete.cases(sc_avg), , drop = FALSE]
sp_avg <- sp_avg[complete.cases(sp_avg), , drop = FALSE]
cat("有效单细胞类型数:", nrow(sc_avg), "\n")
cat("有效空间类型数:", nrow(sp_avg), "\n")

# ==================== 7. 计算相关性矩阵 ====================
cor_matrix <- matrix(NA, nrow = nrow(sc_avg), ncol = nrow(sp_avg),
                     dimnames = list(rownames(sc_avg), rownames(sp_avg)))
for (i in 1:nrow(sc_avg)) {
  for (j in 1:nrow(sp_avg)) {
    cor_matrix[i, j] <- cor(sc_avg[i, ], sp_avg[j, ], method = "spearman", use = "complete.obs")
  }
}

cat("相关性范围:", range(cor_matrix, na.rm = TRUE), "\n")
plot_width <- 14   # 增加宽度
plot_height <- 12  # 增加高度
# ==================== 8. 绘制热图 ====================
pdf("malignant_sc_subcelltype_correlation.pdf", width = plot_width, height = plot_height)

# 计算需要的边距：左边距（行名）+ 右边距（图例）+ 底部边距（列名）+ 顶部边距（标题）
# 使用 cellwidth 和 cellheight 控制每个格子的大小，使热图变小
pheatmap(cor_matrix,
         main = paste("Spearman correlation (", length(overlap.genes), " genes)", sep = ""),
         fontsize = 12,                    # 增大基础字体（原 8 → 12）
         fontsize_row = 10,                # 行名大小（原 8 → 10）
         fontsize_col = 10,                # 列名大小（原 8 → 10）
         fontsize_number = 6,              # 数字大小（如果显示）
         color = colorRampPalette(c("blue", "white", "red"))(100),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         display_numbers = FALSE,
         # 关键：控制格子大小，让热图变小
         cellwidth = 18,                   # 每个格子宽度（像素）
         cellheight = 18,                  # 每个格子高度（像素）
         # 调整边距（增加底部和左边空间给长标签）
         margins = c(12, 12))              # 底部边距12，左边边距12

dev.off()
```
# different expr genes(malignant)
```R
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggrepel)

deg_met1_vs_met2 <- FindMarkers(
  Malignant,
  ident.1 = "Tumor_Metabolic_1",
  ident.2 = "Tumor_Metabolic_2",
  test.use = "wilcox",
  only.pos = FALSE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  min.diff.pct = 0.1
)

deg <- deg_met1_vs_met2
deg$gene <- rownames(deg)

deg$group <- case_when(
  deg$avg_log2FC > 0.5 & deg$p_val_adj < 0.05 ~ "Met1_high",
  deg$avg_log2FC < -0.5 & deg$p_val_adj < 0.05 ~ "Met2_high",
  TRUE ~ "NS"
)

top_deg <- deg %>% filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.5) %>% arrange(desc(abs(avg_log2FC))) %>% head(40)

# 绘图（强制所有基因都有指引线）
p <- ggplot(deg, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = group), size = 1.5, alpha = 0.8) +
  scale_color_manual(values = c("Met1_high"="red","Met2_high"="blue","NS"="gray80")) +
  geom_vline(xintercept=c(-0.5,0.5),linetype="dashed") +
  geom_hline(yintercept=-log10(0.05),linetype="dashed") +
  
  geom_text_repel(
    data = top_deg,
    aes(label=gene),
    size=3,
    
    # 🔥 强制所有标注基因都显示短线！
    min.segment.length = 0,  # 关键：再短也画线
    segment.color = "black",
    segment.size = 0.3,
    force = 5,                # 把文字推开一点
    max.overlaps = 100
  ) +
  
  labs(title="Met1 vs Met2 火山图",x="log2FC",y="-log10(adjusted pvalue)") +
  theme_bw()

ggsave("Volcano_Met1_vs_Met2_ALL_LINE.pdf", p, width=10, height=8, device="pdf")
print(p)

#代谢和其他亚型
Malignant$major_subtype <- Malignant$sub_cell_type

# 将代谢1和代谢2合并为 "Metabolic"
Malignant$major_subtype[Malignant$major_subtype %in% c("Tumor_Metabolic_1", "Tumor_Metabolic_2")] <- "Metabolic"

# 查看合并后的分布
table(Malignant$major_subtype)

# 方法1：使用 is_metabolic 列（推荐）
Malignant$is_metabolic <- ifelse(Malignant$major_subtype == "Metabolic", "Metabolic", "Other")

# 设置身份为 is_metabolic
Idents(Malignant) <- Malignant$is_metabolic

# 比较代谢型 vs 其他
deg_metabolic_vs_other <- FindMarkers(
  Malignant,
  ident.1 = "Metabolic",
  ident.2 = "Other",
  test.use = "wilcox",
  logfc.threshold = 0.25,
  min.pct = 0.25,
  only.pos = FALSE,
  min.diff.pct = 0.1
)

# 查看结果
head(deg_metabolic_vs_other, 20)
```
# GSVA
```R
library(Seurat)
library(GSVA)
library(clusterProfiler)
library(org.Hs.eg.db)

hallmark_pathways <- read.gmt("h.all.v2024.1.Hs.symbols.gmt")
hallmark_list <- split(hallmark_pathways$gene, hallmark_pathways$term)
Malignant <- readRDS("malignant_unknown_anno_0.6.rds")
expression_matrix <- GetAssayData(Malignant, layer = "data")
param <- GSVA::gsvaParam(
  exprData = expression_matrix,  # 表达矩阵
  geneSets = hallmark_list,      # 基因集列表
  kcdf = "Gaussian"              # 核函数
)
gsva_results <- GSVA::gsva(param)
saveRDS(gsva_results, file = "Malignant_gsva_results.rds")
gsva_results <- t(gsva_results)

colnames(gsva_results) <- gsub("HALLMARK_", "", colnames(gsva_results))
Malignant@meta.data <- cbind(Malignant@meta.data, gsva_results)



#top 50
subtypes <- unique(Malignant@meta.data$sub_cell_type)

p_values <- numeric(length = ncol(gsva_results))
names(p_values) <- colnames(gsva_results)

for (pathway in colnames(gsva_results)) {
  # 提取当前通路的 GSVA 得分
  scores <- gsva_results[, pathway]
  
  # 分组比较（例如亚型1 vs 亚型2）
  group1 <- scores[Malignant@meta.data$sub_cell_type == subtypes[1]]
  group2 <- scores[Malignant@meta.data$sub_cell_type == subtypes[2]]
  
  # 进行 t 检验
  test_result <- t.test(group1, group2)
  
  # 保存 p 值
  p_values[pathway] <- test_result$p.value
}
top_50_pathways <- names(sort(p_values))[1:50]
print(top_50_pathways)

top_50_gsva <- gsva_results[, top_50_pathways]
saveRDS(top_50_gsva,file="Malignant_top50.rds")
mean_gsva_by_celltype <- aggregate(top_50_gsva, by = list(Malignant@meta.data$sub_cell_type), FUN = mean)
rownames(mean_gsva_by_celltype) <- mean_gsva_by_celltype$Group.1
mean_gsva_by_celltype <- mean_gsva_by_celltype[, -1]  # 去掉分组列

# 转置矩阵，使行为通路，列为 Epi_cluster
mean_gsva_by_celltype <- t(mean_gsva_by_celltype)

library(pheatmap)
pdf("Malignant_gsva_heatmap.pdf", width = 10, height = 8)
pheatmap(
  mean_gsva_by_celltype,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_colnames = TRUE,
  show_rownames = TRUE,
  fontsize_row = 8,
  fontsize_col = 10
)
dev.off()


# ==================== 加载数据 ====================
gsva_results <- readRDS("Malignant_gsva_results.rds")
gsva_results <- t(gsva_results)

Malignant <- readRDS("malignant_unknown_anno_0.6.rds")

# ==================== 提取代谢1和代谢2的GSVA评分 ====================
gsva_df <- as.data.frame(gsva_results)
gsva_df$subtype <- Malignant$sub_cell_type[rownames(gsva_df)]

# 只保留代谢1和代谢2
gsva_met <- gsva_df[gsva_df$subtype %in% c("Tumor_Metabolic_1", "Tumor_Metabolic_2"), ]

# 计算每个通路的平均活性
mean_gsva <- aggregate(. ~ subtype, data = gsva_met, FUN = mean)
rownames(mean_gsva) <- mean_gsva$subtype
mean_gsva <- mean_gsva[, -1]

# 转置，使行为通路，列为亚型
mean_gsva_t <- t(mean_gsva)

# 计算差值（代谢1 - 代谢2）
diff <- mean_gsva_t[, "Tumor_Metabolic_1"] - mean_gsva_t[, "Tumor_Metabolic_2"]

# ==================== 选择差异最大的前10个通路 ====================
# 代谢1更强的通路（正差值）
top_met1 <- names(sort(diff, decreasing = TRUE)[1:10])
# 代谢2更强的通路（负差值）
top_met2 <- names(sort(diff, decreasing = FALSE)[1:10])

# 合并，并标注方向
selected_pathways <- c(top_met1, top_met2)
direction <- c(rep("Higher in Metabolic1", 10), rep("Higher in Metabolic2",10))
diff_values <- c(diff[top_met1], diff[top_met2])

# 简化通路名称
pathway_labels <- gsub("HALLMARK_", "", selected_pathways)

# ==================== 创建绘图数据框 ====================
plot_df <- data.frame(
  pathway = pathway_labels,
  diff = diff_values,
  direction = direction
)

# 按差值绝对值排序
plot_df <- plot_df[order(plot_df$diff, decreasing = TRUE), ]
plot_df$pathway <- factor(plot_df$pathway, levels = rev(plot_df$pathway))

# ==================== 绘制条形图 ====================
library(ggplot2)

p <- ggplot(plot_df, aes(x = pathway, y = diff, fill = direction)) +
  geom_bar(stat = "identity", width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c("Higher in Metabolic1" = "#E41A1C", 
                                "Higher in Metabolic2" = "#377EB8")) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  labs(
    title = "Pathway activity difference: Metabolic1 vs Metabolic2",
    x = "Pathway",
    y = "GSVA Score Difference (Metabolic1 - Metabolic2)",
    fill = ""
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 10),
    legend.position = "bottom",
    panel.grid.major.y = element_blank()
  )

# 保存为 PDF
ggsave("metabolic1_vs_metabolic2_pathway_diff.pdf", plot = p, width = 10, height = 8)

```
# fibro sc
```R
fib <- subset(obj,subset=cell_type=="Fibroblasts")

fib <- NormalizeData(fib)
fib <- FindVariableFeatures(fib, nfeatures = 3000)
hvgs <- VariableFeatures(fib)
fib <- ScaleData(fib, features = hvgs)
fib <- RunPCA(fib, features = hvgs, npcs = 50)
p <- ElbowPlot(fib, ndims = 30)
ggsave("p3.png",plot=p)

pca_var <- fib[["pca"]]@stdev^2 / sum(fib[["pca"]]@stdev^2)
cumsum_var <- cumsum(pca_var)

# 看达到多少 PC 时累计方差 > 70% / 80% / 90%
which(cumsum_var > 0.7)[1]   # 70%
which(cumsum_var > 0.8)[1]   # 80%
which(cumsum_var > 0.9)[1]   # 90%


library(harmony)
fib <- RunHarmony(fib, "dataset")

fib <- FindNeighbors(fib, reduction = "harmony",dims = 1:23)
fib <- FindClusters(fib, resolution = 0.3)
fib <- RunUMAP(fib,reduction = "harmony", dims = 1:23)

fib <- subset(fib,subset=seurat_clusters %in% c(0,1,2,6,7,9))
fib <- RunUMAP(fib, reduction = "harmony", dims = 1:23, return.model = TRUE)

library(ggplot2)
library(ggsci)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(12)

seurat_clusters <- as.character(unique(fib@meta.data$seurat_clusters))  # 转换为字符向量
num_legend_items <- length(seurat_clusters)  # 图例的个数
max_label_length <- max(nchar(seurat_clusters))  # 图例名称的最大长度

base_width <- 3000  # 基础宽度
base_height <- 3000  # 基础高度
legend_width_factor <- 100  # 每个图例项增加的宽度
label_length_factor <- 10  # 每个字符增加的宽度
dynamic_width <- base_width + (num_legend_items * legend_width_factor) + (max_label_length * label_length_factor)
pdf("fib_sc_clusters.pdf", width = dynamic_width/300, height = base_height/300)  # 转换为英寸

DimPlot(fib, reduction = "umap", label = TRUE, pt.size = 1, label.size = 8) +
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    ggtitle(NULL) +
    scale_color_manual(values = npg_extended) +
    coord_fixed(ratio = 1) +
    guides(color = guide_legend(title = NULL, override.aes = list(size = 5))) +
    theme(
        text = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 28, color = "black"),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title.x = element_text(size = 36, face = "bold", color = "black"),
        axis.title.y = element_text(size = 36, face = "bold", color = "black", margin = margin(r = 20)),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 28, face = "bold", color = "black"),
        legend.title = element_text(size = 28, face = "bold", color = "black"),
        legend.position = "right",
        legend.box.margin = margin(0, 0, 0, 0),
        legend.key = element_blank(),
        legend.background = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5),
        aspect.ratio = 1,
        plot.margin = margin(10, 50, 10, 10)
    )

dev.off()

p <- DimPlot(fib, reduction = "umap", group.by = "orig.ident", label = FALSE)
ggsave("fib_orig.png",plot=p)

library(dplyr)
fib_markers <- FindAllMarkers(fib, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "wilcox")
fib_significant_markers <- subset(fib_markers, p_val_adj < 0.05)
#write.csv(fib_significant_markers, "fib_all_marker.csv")
fib_significant_markers <- fib_significant_markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
write.csv(fib_significant_markers, "fib_sc_top_marker_50.csv")


identity_mapping <- c(
  "0" = "Fib_matCAF",
  "1" = "Fib_iCAF",
  "2" = "Fib_myCAF",
  "3" = "Fib_pCAF"
)

identity_mapping <- c(
  "0" = "Fib_myCAF",
  "1" = "Fib_matCAF",
  "2" = "Fib_iCAF",
  "6" = "Fib_Pericyte",
  "7" = "Fib_pCAF"
)
fib$seurat_clusters <- factor(fib$seurat_clusters)

table(fib$seurat_clusters)
sub_cell_type <- identity_mapping[fib@meta.data$seurat_clusters]
fib@meta.data$sub_cell_type <- sub_cell_type


npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(5)

cluster_colors <- c(
  # 浅色系在前（浅蓝第一个）
  "#A0CBE8",   # 1 浅蓝（你要的放第一个）
  "#FFBE7D",   # 2 浅橙
  "#8CD17D",   # 3 浅绿
  "#86BCB6",   # 4 浅青
  "#FF9D9A",   # 5 浅红
  "#FABFD2",   # 6 浅粉
  "#D4A6C8",   # 7 浅紫
  "#F1CE63",   # 8 淡黄
  "#D7B5A6",   # 9 浅棕
  "#B07AA1",   # 10 紫灰
  "#BAB0AC",   # 11 浅灰
  
  # 深色系在后
  "#4E79A7",   # 12 深蓝
  "#F28E2B",   # 13 橙
  "#59A14F",   # 14 深绿
  "#E15759",   # 15 红
  "#499894",   # 16 青绿
  "#D37295",   # 17 玫红
  "#B6992D",   # 18 土黄
  "#9D7660",   # 19 棕
  "#79706E"    # 20 深灰
)

cell_types <- as.character(unique(fib@meta.data$sub_cell_type))
num_legend_items <- length(cell_types)  # 图例的个数
max_label_length <- max(nchar(cell_types))  # 图例名称的最大长度

base_width <- 3000  # 基础宽度
base_height <- 3000  # 基础高度
legend_width_factor <- 100  # 每个图例项增加的宽度
label_length_factor <- 10  # 每个字符增加的宽度
dynamic_width <- base_width + (num_legend_items * legend_width_factor) + (max_label_length * label_length_factor)
pdf("fib_sc_annotation.pdf", width = dynamic_width/300, height = base_height/300)
DimPlot(fib, reduction = "umap", label = TRUE, pt.size = 1, group.by = "sub_cell_type", label.size = 4) +
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    ggtitle(NULL) +
    scale_color_manual(values = cluster_colors) +
    coord_fixed(ratio = 1) +
    guides(color = guide_legend(title = NULL, override.aes = list(size = 5))) +
    theme(
        text = element_text(size = 8, face = "bold"),
        axis.text.x = element_text(size = 28, color = "black"),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title.x = element_text(size = 36, face = "bold", color = "black"),
        axis.title.y = element_text(size = 36, face = "bold", color = "black", margin = margin(r = 20)),  # 增加右侧间距
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 16, face = "bold", color = "black"),
        legend.title = element_text(size = 16, face = "bold", color = "black"),
        legend.position = "right",
        legend.box.margin = margin(0, 0, 0, 0),
        legend.key = element_blank(),
        legend.background = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5),
        aspect.ratio = 1,
        plot.margin = margin(10, 50, 10, 10)
    )
dev.off()
saveRDS(fib,file="fib_sc_anno.rds")
fib_sc <- fib
```
# fib mapping
```R
library(Seurat)
library(dplyr)

fib_sc <- readRDS("fib_sc_anno.rds")
fib <- readRDS("fib_anno_new.rds")

# ===================== 【1】内存 + 并行关闭 =====================
options(future.globals.maxSize = 100 * 1024^3)  # 100G
future::plan("sequential")                     # 关闭并行，最关键！

# ===================== 【2】强制清理 =====================
DefaultAssay(fib) <- "RNA"
DefaultAssay(fib_sc) <- "RNA"
fib@assays$SCT <- NULL
fib_sc@assays$SCT <- NULL

# ===================== 【3】合并图层（唯一能解决你问题的步骤） =====================
fib <- JoinLayers(fib)       # 必加！
fib_sc <- JoinLayers(fib_sc) # 必加！

# ===================== 【4】重新完整跑一遍流程 =====================
# 标准化
fib <- NormalizeData(fib)
fib_sc <- NormalizeData(fib_sc)

# 高变基因
fib <- FindVariableFeatures(fib, nfeatures = 2000)
fib_sc <- FindVariableFeatures(fib_sc, nfeatures = 2000)

# 共同基因
overlap.genes <- intersect(rownames(fib), rownames(fib_sc))

# 缩放！！！（必须在 JoinLayers 之后）
fib <- ScaleData(fib, features = overlap.genes)
fib_sc <- ScaleData(fib_sc, features = overlap.genes)

# PCA
fib <- RunPCA(fib, features = overlap.genes, npcs = 50)
fib_sc <- RunPCA(fib_sc, features = overlap.genes, npcs = 50)

# ===================== 【5】锚点 + 注释 =====================
transfer.anchors <- FindTransferAnchors(
  reference = fib_sc,
  query = fib,
  dims = 1:15,                # 从20 → 15，大幅降内存！
  features = overlap.genes,
  
  # 三个内存杀手参数，我帮你全关！
  reference.assay = "RNA",
  query.assay = "RNA",
  normalization.method = "LogNormalize",
  
  # 【爆内存救星】
  k.anchor = 5,               # 越小占内存越少
  approx.pca = TRUE           # 快速+低内存
)

predictions <- TransferData(
  anchorset = transfer.anchors,
  refdata = fib_sc$sub_cell_type,
  dims = 1:15,
  k.weight = 4
)
fib$cell_type_pred <-NA
fib <- AddMetaData(fib, predictions)
fib$cell_type_pred <- fib$predicted.id
table(fib@meta.data$cell_type_pred,fib@meta.data$sub_cell_type)

library(ggplot2)
library(ggsci)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(6)

cell_types <- as.character(unique(Malignant@meta.data$cell_type_pred))
num_legend_items <- length(cell_types)  # 图例的个数
max_label_length <- max(nchar(cell_types))  # 图例名称的最大长度

base_width <- 3000  # 基础宽度
base_height <- 3000  # 基础高度
legend_width_factor <- 100  # 每个图例项增加的宽度
label_length_factor <- 10  # 每个字符增加的宽度
dynamic_width <- base_width + (num_legend_items * legend_width_factor) + (max_label_length * label_length_factor)
pdf("malignant_predict_annotation.pdf", width = dynamic_width/300, height = base_height/300)
DimPlot(Malignant, reduction = "umap", label = TRUE, pt.size = 1, group.by = "cell_type_pred", label.size = 4) +
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    ggtitle(NULL) +
    scale_color_manual(values = npg_pal) +
    coord_fixed(ratio = 1) +
    guides(color = guide_legend(title = NULL, override.aes = list(size = 5))) +
    theme(
        text = element_text(size = 8, face = "bold"),
        axis.text.x = element_text(size = 28, color = "black"),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title.x = element_text(size = 36, face = "bold", color = "black"),
        axis.title.y = element_text(size = 36, face = "bold", color = "black", margin = margin(r = 20)),  # 增加右侧间距
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 16, face = "bold", color = "black"),
        legend.title = element_text(size = 16, face = "bold", color = "black"),
        legend.position = "right",
        legend.box.margin = margin(0, 0, 0, 0),
        legend.key = element_blank(),
        legend.background = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5),
        aspect.ratio = 1,
        plot.margin = margin(10, 50, 10, 10)
    )
dev.off()

library(Seurat)
library(ggplot2)
library(patchwork)

# 修复因子水平
fib$sub_cell_type <- as.character(fib$sub_cell_type)
fib_sc$sub_cell_type <- as.character(fib_sc$sub_cell_type)

# 定义基因顺序（基于新注释）
all_genes <- c(
  # apCAF_1 (cluster 1)
  "CD74", "IGKC", "IGHG1/2",
  # apCAF_2 (cluster 3)
  "CD74", "HLA-DRB", "IGKC",
  # iCAF_1 (cluster 2)
  "PTEN", "SIGMAR1", "GADD45B",
  # iCAF_2 (cluster 4)
  "CFD", "C1S", "SERPINE1",
  # matCAF (cluster 5)
  "COL1A1", "COL3A1", "SPARC",
  # myCAF (cluster 0)
  "ACTA2", "TAGLN", "MYL9",
  # pCAF (单细胞特有)
  "TOP2A", "BIRC5", "RRM2"
)
all_genes <- unique(all_genes)

# 空间组气泡图
p1 <- DotPlot(fib, features = all_genes, group.by = "sub_cell_type", scale = TRUE) +
  scale_color_gradient(low = "#D6EAF8", high = "#1F618D", name = "Average Expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Spatial - Fibroblast Subtypes", x = "", y = "") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# 单细胞气泡图
p2 <- DotPlot(fib_sc, features = all_genes, group.by = "sub_cell_type", scale = TRUE) +
  scale_color_gradient(low = "#D6EAF8", high = "#1F618D", name = "Average Expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "scRNA-seq - Fibroblast Subtypes", x = "", y = "") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# 拼图
p_combined <- p1 / p2 + 
  plot_annotation(title = "Fibroblast Marker Genes: Spatial vs scRNA-seq") &
  theme(plot.title = element_text(hjust = 0.5, size = 14))

ggsave("bubble_plot_fibroblast_final.pdf", plot = p_combined, width = 14, height = 10)
```
# CN 
```R
#mPT <- readRDS("mPT_with_regions.rds")
Malignant <- readRDS("malignant_unknown_anno_0.6.rds")
fib <- readRDS("fib_anno.rds")
macro <- readRDS("macro_anno_no_rm.rds")
obj <- readRDS("YA2025263-1_fin.rds")

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
                       "Basal cells" = "Malignant cells"
)

# 检查
table(obj$CellType)

obj@meta.data$sub_cell_type <- NA
common_cells <- intersect(rownames(Malignant@meta.data), rownames(obj@meta.data))
obj@meta.data[common_cells, "sub_cell_type"] <- Malignant@meta.data[common_cells, "sub_cell_type"]

common_cells <- intersect(rownames(fib@meta.data), rownames(obj@meta.data))
obj@meta.data[common_cells, "sub_cell_type"] <- fib@meta.data[common_cells, "sub_cell_type"]

common_cells <- intersect(rownames(macro@meta.data), rownames(obj@meta.data))
obj@meta.data[common_cells, "sub_cell_type"] <- macro@meta.data[common_cells, "sub_cell_type"]

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

mPT <- subset(obj,subset=tissue=="mPT")
```
# macro sc
```R
obj_sc <- readRDS("sc_obj_anno.rds")
macro <- subset(obj_sc,subset=cell_type=="Macrophages")
library(Seurat)
macro <- NormalizeData(macro)
macro <- FindVariableFeatures(macro, nfeatures = 2000)
hvgs <- VariableFeatures(macro)
macro <- ScaleData(macro, features = hvgs)
macro <- RunPCA(macro, features = hvgs, npcs = 50)
library(ggplot2)
p <- ElbowPlot(macro, ndims = 30)
ggsave("p3.png",plot=p)

pca_var <- macro[["pca"]]@stdev^2 / sum(macro[["pca"]]@stdev^2)
cumsum_var <- cumsum(pca_var)

# 看达到多少 PC 时累计方差 > 70% / 80% / 90%
which(cumsum_var > 0.7)[1]   # 70%
which(cumsum_var > 0.8)[1]   # 80%
which(cumsum_var > 0.9)[1]   # 90%

library(harmony)
macro <- RunHarmony(macro, "dataset")
macro <- FindNeighbors(macro, reduction = "harmony",dims = 1:21)
macro <- FindClusters(macro, resolution = 0.3)
macro <- RunUMAP(macro, reduction = "harmony", dims = 1:21)

macro <- subset(macro,subset=seurat_clusters %in% c(1,2,3,5,6,8,9,10,11))
macro <- RunUMAP(macro, reduction = "harmony", dims = 1:21)

library(ggplot2)
library(ggsci)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(17)

seurat_clusters <- as.character(unique(macro@meta.data$seurat_clusters))  # 转换为字符向量
num_legend_items <- length(seurat_clusters)  # 图例的个数
max_label_length <- max(nchar(seurat_clusters))  # 图例名称的最大长度

base_width <- 3000  # 基础宽度
base_height <- 3000  # 基础高度
legend_width_factor <- 100  # 每个图例项增加的宽度
label_length_factor <- 10  # 每个字符增加的宽度
dynamic_width <- base_width + (num_legend_items * legend_width_factor) + (max_label_length * label_length_factor)
pdf("macro_sc_clusters.pdf", width = dynamic_width/300, height = base_height/300)  # 转换为英寸

DimPlot(macro, reduction = "umap", label = TRUE, pt.size = 1, label.size = 8) +
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    ggtitle(NULL) +
    scale_color_manual(values = npg_extended) +
    coord_fixed(ratio = 1) +
    guides(color = guide_legend(title = NULL, override.aes = list(size = 5))) +
    theme(
        text = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 28, color = "black"),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title.x = element_text(size = 36, face = "bold", color = "black"),
        axis.title.y = element_text(size = 36, face = "bold", color = "black", margin = margin(r = 20)),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 28, face = "bold", color = "black"),
        legend.title = element_text(size = 28, face = "bold", color = "black"),
        legend.position = "right",
        legend.box.margin = margin(0, 0, 0, 0),
        legend.key = element_blank(),
        legend.background = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5),
        aspect.ratio = 1,
        plot.margin = margin(10, 50, 10, 10)
    )

dev.off()

p <- DimPlot(macro, reduction = "umap", group.by = "dataset", label = FALSE)
ggsave("macro_orig.png",plot=p)

library(dplyr)
macro_markers <- FindAllMarkers(macro, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "wilcox")
macro_significant_markers <- subset(macro_markers, p_val_adj < 0.05)
#write.csv(macro_significant_markers, "macro_all_marker.csv")
macro_significant_markers <- macro_significant_markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
write.csv(macro_significant_markers, "macro_sc_top_marker_50.csv")

identity_mapping <- c(
  "0" = "Macro_Lipid-associated_1",
  "1" = "Macro_Inflammatory",
  "2" = "Macro_Immunoregulatory",
  "3" = "Macro_Matrix-remodeling",
  "4" = "Macro_Stress-response",
  "5" = "Macro_Lipid-associated_2",
  "6" = "Macro_Lipid-associated_3",
  "7" = "Macro_Proliferative"
)

identity_mapping <- c(
  "1"  = "Macro_Lipid-associated_1",
  "2"  = "Macro_Tissue-resident",
  "3"  = "Macro_Pro-inflammatory_1",
  "5"  = "Macro_Lipid-associated_2",
  "6"  = "Macro_Immunosuppressive",
  "8"  = "Macro_Pro-inflammatory_2",
  "9"  = "Macro_Stress-response",
  "10" = "Macro_Proliferative",
  "11" = "Macro_Complement-enriched"
)

macro$seurat_clusters <- factor(macro$seurat_clusters)

table(macro$seurat_clusters)

sub_cell_type <- identity_mapping[macro@meta.data$seurat_clusters]
macro@meta.data$sub_cell_type <- sub_cell_type


npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(5)
cluster_colors <- c(
  # 浅色系在前（浅蓝第一个）
  "#A0CBE8",   # 1 浅蓝（你要的放第一个）
  "#FFBE7D",   # 2 浅橙
  "#8CD17D",   # 3 浅绿
  "#86BCB6",   # 4 浅青
  "#FF9D9A",   # 5 浅红
  "#FABFD2",   # 6 浅粉
  "#D4A6C8",   # 7 浅紫
  "#F1CE63",   # 8 淡黄
  "#D7B5A6",   # 9 浅棕
  "#B07AA1",   # 10 紫灰
  "#BAB0AC",   # 11 浅灰
  
  # 深色系在后
  "#4E79A7",   # 12 深蓝
  "#F28E2B",   # 13 橙
  "#59A14F",   # 14 深绿
  "#E15759",   # 15 红
  "#499894",   # 16 青绿
  "#D37295",   # 17 玫红
  "#B6992D",   # 18 土黄
  "#9D7660",   # 19 棕
  "#79706E"    # 20 深灰
)
cell_types <- as.character(unique(macro@meta.data$sub_cell_type))
num_legend_items <- length(cell_types)  # 图例的个数
max_label_length <- max(nchar(cell_types))  # 图例名称的最大长度

base_width <- 3000  # 基础宽度
base_height <- 3000  # 基础高度
legend_width_factor <- 100  # 每个图例项增加的宽度
label_length_factor <- 10  # 每个字符增加的宽度
dynamic_width <- base_width + (num_legend_items * legend_width_factor) + (max_label_length * label_length_factor)
pdf("macro_sc_annotation.pdf", width = dynamic_width/300, height = base_height/300)
DimPlot(macro, reduction = "umap", label = TRUE, pt.size = 1, group.by = "sub_cell_type", label.size = 4) +
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    ggtitle(NULL) +
    scale_color_manual(values = cluster_colors) +
    coord_fixed(ratio = 1) +
    guides(color = guide_legend(title = NULL, override.aes = list(size = 5))) +
    theme(
        text = element_text(size = 8, face = "bold"),
        axis.text.x = element_text(size = 28, color = "black"),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title.x = element_text(size = 36, face = "bold", color = "black"),
        axis.title.y = element_text(size = 36, face = "bold", color = "black", margin = margin(r = 20)),  # 增加右侧间距
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 16, face = "bold", color = "black"),
        legend.title = element_text(size = 16, face = "bold", color = "black"),
        legend.position = "right",
        legend.box.margin = margin(0, 0, 0, 0),
        legend.key = element_blank(),
        legend.background = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5),
        aspect.ratio = 1,
        plot.margin = margin(10, 50, 10, 10)
    )
dev.off()
saveRDS(macro,file="macro_sc_anno.rds")
```
# macro mapping
```R
macro <- readRDS("macro_anno_no_rm.rds")
library(Seurat)
library(dplyr)

macro_sc <- readRDS("macro_sc_anno.rds")


# ===================== 【1】内存 + 并行关闭 =====================
options(future.globals.maxSize = 100 * 1024^3)  # 100G
future::plan("sequential")                     # 关闭并行，最关键！

# ===================== 【2】强制清理 =====================
DefaultAssay(macro) <- "RNA"
DefaultAssay(macro_sc) <- "RNA"
macro@assays$SCT <- NULL
macro_sc@assays$SCT <- NULL

# ===================== 【3】合并图层（唯一能解决你问题的步骤） =====================
macro <- JoinLayers(macro)       # 必加！
macro_sc <- JoinLayers(macro_sc) # 必加！

# ===================== 【4】重新完整跑一遍流程 =====================
# 标准化
macro <- NormalizeData(macro)
macro_sc <- NormalizeData(macro_sc)

# 高变基因
macro <- FindVariableFeatures(macro, nfeatures = 2000)
macro_sc <- FindVariableFeatures(macro_sc, nfeatures = 2000)

# 共同基因
overlap.genes <- intersect(rownames(macro), rownames(macro_sc))

# 缩放！！！（必须在 JoinLayers 之后）
macro <- ScaleData(macro, features = overlap.genes)
macro_sc <- ScaleData(macro_sc, features = overlap.genes)

# PCA
macro <- RunPCA(macro, features = overlap.genes, npcs = 50)
macro_sc <- RunPCA(macro_sc, features = overlap.genes, npcs = 50)

# ===================== 【5】锚点 + 注释 =====================
transfer.anchors <- FindTransferAnchors(
  reference = macro_sc,
  query = macro,
  dims = 1:15,                # 从20 → 15，大幅降内存！
  features = overlap.genes,
  
  # 三个内存杀手参数，我帮你全关！
  reference.assay = "RNA",
  query.assay = "RNA",
  normalization.method = "LogNormalize",
  
  # 【爆内存救星】
  k.anchor = 5,               # 越小占内存越少
  approx.pca = TRUE           # 快速+低内存
)

predictions <- TransferData(
  anchorset = transfer.anchors,
  refdata = macro_sc$sub_cell_type,
  dims = 1:15,
  k.weight = 10
)
macro$cell_type_pred <- NA
macro <- AddMetaData(macro, predictions)
macro$cell_type_pred <- macro$predicted.id
table(macro@meta.data$cell_type_pred, macro@meta.data$sub_cell_type)

library(ggplot2)
library(ggsci)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(6)

cell_types <- as.character(unique(macro@meta.data$cell_type_pred))
num_legend_items <- length(cell_types)  # 图例的个数
max_label_length <- max(nchar(cell_types))  # 图例名称的最大长度

base_width <- 3000  # 基础宽度
base_height <- 3000  # 基础高度
legend_width_factor <- 100  # 每个图例项增加的宽度
label_length_factor <- 10  # 每个字符增加的宽度
dynamic_width <- base_width + (num_legend_items * legend_width_factor) + (max_label_length * label_length_factor)
pdf("macro_predict_annotation.pdf", width = dynamic_width/300, height = base_height/300)
DimPlot(macro, reduction = "umap", label = TRUE, pt.size = 1, group.by = "cell_type_pred", label.size = 4) +
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    ggtitle(NULL) +
    scale_color_manual(values = npg_pal) +
    coord_fixed(ratio = 1) +
    guides(color = guide_legend(title = NULL, override.aes = list(size = 5))) +
    theme(
        text = element_text(size = 8, face = "bold"),
        axis.text.x = element_text(size = 28, color = "black"),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title.x = element_text(size = 36, face = "bold", color = "black"),
        axis.title.y = element_text(size = 36, face = "bold", color = "black", margin = margin(r = 20)),  # 增加右侧间距
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 16, face = "bold", color = "black"),
        legend.title = element_text(size = 16, face = "bold", color = "black"),
        legend.position = "right",
        legend.box.margin = margin(0, 0, 0, 0),
        legend.key = element_blank(),
        legend.background = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5),
        aspect.ratio = 1,
        plot.margin = margin(10, 50, 10, 10)
    )
dev.off()

all_genes <- c(
  # Antigen-presenting_Macro
  "IGKC", "CD74", "IGHG4",        # 用 IGHG4 替代 IGHG1/2
  # Homeostatic_Macro
  "APOE", "CD68", "TXNIP",
  # Immunoregulatory_Macro
  "CCL21", "PTGDS", "CCL19",
  # Lipid-associated_Macro
  "CHIT1", "ACP5", "CYP27A1",
  # Matrix-remodeling_Macro
  "POSTN", "COL1A1", "MMP2",
  # Phagocytic_Macro
  "KRT5", "SPP1", "KRT19",
  # Stress-response_Macro
  "PTEN", "SIGMAR1", "VPS29",
  # 单细胞特有亚型
  "S100A8", "IL1B", "VEGFA",      # Inflammatory_Macro
  "TOP2A", "MKI67", "BIRC5"       # Proliferative_Macro
)
all_genes <- unique(all_genes)

library(Seurat)
library(ggplot2)
library(patchwork)

macro$sub_cell_type <- as.character(macro$sub_cell_type)
macro_sc$sub_cell_type <- as.character(macro_sc$sub_cell_type)

# 空间组
p1 <- DotPlot(macro, features = all_genes, group.by = "sub_cell_type", scale = TRUE) +
  scale_color_gradient(low = "#D6EAF8", high = "#1F618D", name = "Average Expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
  labs(title = "Spatial - Macrophage Subtypes", x = "", y = "") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# 单细胞
p2 <- DotPlot(macro_sc, features = all_genes, group.by = "sub_cell_type", scale = TRUE) +
  scale_color_gradient(low = "#D6EAF8", high = "#1F618D", name = "Average Expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
  labs(title = "scRNA-seq - Macrophage Subtypes", x = "", y = "") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# 拼图
p_combined <- p1 / p2 + 
  plot_annotation(title = "Macrophage Marker Genes: Spatial vs scRNA-seq") &
  theme(plot.title = element_text(hjust = 0.5, size = 14))

ggsave("bubble_plot_macrophage_optimized.pdf", plot = p_combined, width = 16, height = 12)
```
# mast
```R
mast <- subset(obj, subset = CellType == "Mast cells")
library(Seurat)
mast <- NormalizeData(mast)
mast <- FindVariableFeatures(mast, nfeatures = 2000)
hvgs <- VariableFeatures(mast)
mast <- ScaleData(mast, features = hvgs)
mast <- RunPCA(mast, features = hvgs, npcs = 50)
library(ggplot2)
p <- ElbowPlot(mast, ndims = 30)
ggsave("p3.png",plot=p)

mast <- FindNeighbors(mast, dims = 1:15)
mast <- FindClusters(mast, resolution = 0.6)
mast <- RunUMAP(mast, dims = 1:15)

library(ggplot2)
library(ggsci)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(11)#0.6

seurat_clusters <- as.character(unique(mast@meta.data$seurat_clusters))  # 转换为字符向量
num_legend_items <- length(seurat_clusters)  # 图例的个数
max_label_length <- max(nchar(seurat_clusters))  # 图例名称的最大长度

base_width <- 3000  # 基础宽度
base_height <- 3000  # 基础高度
legend_width_factor <- 100  # 每个图例项增加的宽度
label_length_factor <- 10  # 每个字符增加的宽度
dynamic_width <- base_width + (num_legend_items * legend_width_factor) + (max_label_length * label_length_factor)
pdf("mast_clusters.pdf", width = dynamic_width/300, height = base_height/300)  # 转换为英寸

DimPlot(mast, reduction = "umap", label = TRUE, pt.size = 1, label.size = 8) +
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    ggtitle(NULL) +
    scale_color_manual(values = npg_extended) +
    coord_fixed(ratio = 1) +
    guides(color = guide_legend(title = NULL, override.aes = list(size = 5))) +
    theme(
        text = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 28, color = "black"),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title.x = element_text(size = 36, face = "bold", color = "black"),
        axis.title.y = element_text(size = 36, face = "bold", color = "black", margin = margin(r = 20)),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 28, face = "bold", color = "black"),
        legend.title = element_text(size = 28, face = "bold", color = "black"),
        legend.position = "right",
        legend.box.margin = margin(0, 0, 0, 0),
        legend.key = element_blank(),
        legend.background = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5),
        aspect.ratio = 1,
        plot.margin = margin(10, 50, 10, 10)
    )

dev.off()

library(dplyr)
mast_markers <- FindAllMarkers(mast, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "wilcox")
mast_significant_markers <- subset(mast_markers, p_val_adj < 0.05)
#write.csv(mast_significant_markers, "mast_all_marker.csv")
mast_significant_markers <- mast_significant_markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
write.csv(mast_significant_markers, "mast_top_marker_50.csv")


#0.6
identity_mapping <- c(
  "0" = "Mast_Stress-response",
  "1" = "Mast_Metabolic_1",
  "2" = "Mast_Homeostatic",
  "3" = "Mast_Matrix-remodeling",
  "4" = "Mast_Antigen-presenting_1",
  "5" = "Mast_Immunoregulatory",
  "6" = "Mast_Antigen-presenting_2",
  "7" = "Mast_Epithelial-associated",
  "8" = "Mast_Metabolic_2",
  "9" = "Mast_Inflammatory",
  "10" = "Mast_MCTC"
)
sub_cell_type <- identity_mapping[mast@meta.data$seurat_clusters]
mast@meta.data$sub_cell_type <- sub_cell_type

library(ggplot2)
library(ggsci)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(11)

cell_types <- as.character(unique(mast@meta.data$sub_cell_type))
num_legend_items <- length(cell_types)  # 图例的个数
max_label_length <- max(nchar(cell_types))  # 图例名称的最大长度

base_width <- 3000  # 基础宽度
base_height <- 3000  # 基础高度
legend_width_factor <- 100  # 每个图例项增加的宽度
label_length_factor <- 10  # 每个字符增加的宽度
dynamic_width <- base_width + (num_legend_items * legend_width_factor) + (max_label_length * label_length_factor)
pdf("mast_annotation-0.6.pdf", width = dynamic_width/300, height = base_height/300)
DimPlot(mast, reduction = "umap", label = TRUE, pt.size = 1, group.by = "sub_cell_type", label.size = 4) +
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    ggtitle(NULL) +
    scale_color_manual(values = npg_extended) +
    coord_fixed(ratio = 1) +
    guides(color = guide_legend(title = NULL, override.aes = list(size = 5))) +
    theme(
        text = element_text(size = 8, face = "bold"),
        axis.text.x = element_text(size = 28, color = "black"),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title.x = element_text(size = 36, face = "bold", color = "black"),
        axis.title.y = element_text(size = 36, face = "bold", color = "black", margin = margin(r = 20)),  # 增加右侧间距
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 16, face = "bold", color = "black"),
        legend.title = element_text(size = 16, face = "bold", color = "black"),
        legend.position = "right",
        legend.box.margin = margin(0, 0, 0, 0),
        legend.key = element_blank(),
        legend.background = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5),
        aspect.ratio = 1,
        plot.margin = margin(10, 50, 10, 10)
    )
dev.off()

#比例
library(tidyverse)

# 计算各亚群在不同组织中的占比
prop_data <-mast@meta.data %>%
  group_by(tissue, sub_cell_type) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(tissue) %>%
  mutate(proportion = count / sum(count) * 100) %>%
  filter(tissue %in% c("mPT", "nmPT")) %>%
  mutate(tissue = factor(tissue, levels = c("mPT", "nmPT")))

library(ggplot2)
library(ggsci)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(11)


prop_data <- prop_data %>%
  mutate(sub_cell_type = as.character(sub_cell_type))

# 绘制堆叠条形图
p <- ggplot(prop_data, aes(x = tissue, y = proportion, fill = sub_cell_type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = npg_extended) +
  labs(title = "Cell Type Distribution: mPT vs nmPT",
       x = "Tissue", y = "Proportion (%)", fill = "Cell Type") +
  theme_minimal() +
  theme(
    # 坐标轴线条和刻度
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.15, "cm"),
    # 坐标轴文本
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 10),
    axis.text.y = element_text(color = "black", size = 10),
    # 坐标轴标题
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    # 面板背景
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    # 图例
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )


# 保存为PDF
pdf("mast_celltype_distribution_tumor0.6.pdf", width = 4, height = 6)
print(p)
dev.off()

prop_data <- mast@meta.data %>%
  group_by(tissue, sub_cell_type) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(tissue) %>%
  mutate(proportion = count / sum(count) * 100) %>%
  filter(tissue %in% c("metLN", "negLN")) %>%
  mutate(tissue = factor(tissue, levels = c("metLN", "negLN")))

library(ggplot2)
library(ggsci)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(11)


prop_data <- prop_data %>%
  mutate(sub_cell_type = as.character(sub_cell_type))

# 绘制堆叠条形图
p <- ggplot(prop_data, aes(x = tissue, y = proportion, fill = sub_cell_type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = npg_extended) +
  labs(title = "Cell Type Distribution: metLN vs negLN",
       x = "Tissue", y = "Proportion (%)", fill = "Cell Type") +
  theme_minimal() +
  theme(
    # 坐标轴线条和刻度
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.15, "cm"),
    # 坐标轴文本
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 10),
    axis.text.y = element_text(color = "black", size = 10),
    # 坐标轴标题
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    # 面板背景
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    # 图例
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )


# 保存为PDF
pdf("mast_celltype_distribution_RLN_P_RLN_N_NRLN-0.6.pdf", width = 4, height = 6)
print(p)
dev.off()


```
# mast sc
```R
library(Seurat)
obj_sc <- readRDS("sc_obj_anno.rds")
mast_sc <- subset(obj, subset = cell_type == "Mast cells")
mast_sc <- NormalizeData(mast_sc)
mast_sc <- FindVariableFeatures(mast_sc, nfeatures = 2000)
hvgs <- VariableFeatures(mast_sc)
mast_sc <- ScaleData(mast_sc, features = hvgs)
mast_sc <- RunPCA(mast_sc, features = hvgs, npcs = 50)
p <- ElbowPlot(mast_sc, ndims = 30)
ggsave("p3.png",plot=p)

mast_sc <- FindNeighbors(mast_sc, dims = 1:20)
mast_sc <- FindClusters(mast_sc, resolution = 0.6)
mast_sc <- RunUMAP(mast_sc, dims = 1:20)

library(ggplot2)
library(ggsci)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(11)#0.6

seurat_clusters <- as.character(unique(mast_sc@meta.data$seurat_clusters))  # 转换为字符向量
num_legend_items <- length(seurat_clusters)  # 图例的个数
max_label_length <- max(nchar(seurat_clusters))  # 图例名称的最大长度

base_width <- 3000  # 基础宽度
base_height <- 3000  # 基础高度
legend_width_factor <- 100  # 每个图例项增加的宽度
label_length_factor <- 10  # 每个字符增加的宽度
dynamic_width <- base_width + (num_legend_items * legend_width_factor) + (max_label_length * label_length_factor)
pdf("mast_sc_clusters.pdf", width = dynamic_width/300, height = base_height/300)  # 转换为英寸

DimPlot(mast_sc, reduction = "umap", label = TRUE, pt.size = 1, label.size = 8) +
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    ggtitle(NULL) +
    scale_color_manual(values = npg_extended) +
    coord_fixed(ratio = 1) +
    guides(color = guide_legend(title = NULL, override.aes = list(size = 5))) +
    theme(
        text = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 28, color = "black"),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title.x = element_text(size = 36, face = "bold", color = "black"),
        axis.title.y = element_text(size = 36, face = "bold", color = "black", margin = margin(r = 20)),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 28, face = "bold", color = "black"),
        legend.title = element_text(size = 28, face = "bold", color = "black"),
        legend.position = "right",
        legend.box.margin = margin(0, 0, 0, 0),
        legend.key = element_blank(),
        legend.background = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5),
        aspect.ratio = 1,
        plot.margin = margin(10, 50, 10, 10)
    )

dev.off()

library(dplyr)
mast_markers <- FindAllMarkers(mast_sc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "wilcox")
mast_significant_markers <- subset(mast_markers, p_val_adj < 0.05)
#write.csv(mast_significant_markers, "mast_all_marker.csv")
mast_significant_markers <- mast_significant_markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
write.csv(mast_significant_markers, "mast_sc_top_marker_50.csv")


#0.6
identity_mapping <- c(
  "0" = "Inflammatory_Mast",
  "1" = "Classical_MCTC",
  "2" = "Tissue-resident_Mast",
  "3" = "Proliferative_Mast"
)
sub_cell_type <- identity_mapping[mast_sc@meta.data$seurat_clusters]
mast_sc@meta.data$sub_cell_type <- sub_cell_type

library(ggplot2)
library(ggsci)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(11)

cell_types <- as.character(unique(mast_sc@meta.data$sub_cell_type))
num_legend_items <- length(cell_types)  # 图例的个数
max_label_length <- max(nchar(cell_types))  # 图例名称的最大长度

base_width <- 3000  # 基础宽度
base_height <- 3000  # 基础高度
legend_width_factor <- 100  # 每个图例项增加的宽度
label_length_factor <- 10  # 每个字符增加的宽度
dynamic_width <- base_width + (num_legend_items * legend_width_factor) + (max_label_length * label_length_factor)
pdf("mast_sc_annotation-0.6.pdf", width = dynamic_width/300, height = base_height/300)
DimPlot(mast_sc, reduction = "umap", label = TRUE, pt.size = 1, group.by = "sub_cell_type", label.size = 4) +
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    ggtitle(NULL) +
    scale_color_manual(values = npg_pal) +
    coord_fixed(ratio = 1) +
    guides(color = guide_legend(title = NULL, override.aes = list(size = 5))) +
    theme(
        text = element_text(size = 8, face = "bold"),
        axis.text.x = element_text(size = 28, color = "black"),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title.x = element_text(size = 36, face = "bold", color = "black"),
        axis.title.y = element_text(size = 36, face = "bold", color = "black", margin = margin(r = 20)),  # 增加右侧间距
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 16, face = "bold", color = "black"),
        legend.title = element_text(size = 16, face = "bold", color = "black"),
        legend.position = "right",
        legend.box.margin = margin(0, 0, 0, 0),
        legend.key = element_blank(),
        legend.background = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5),
        aspect.ratio = 1,
        plot.margin = margin(10, 50, 10, 10)
    )
dev.off()

all_genes <- c(
  # 1 Antigen-presenting_Mast_1
  "HLA-DRA", "CD74", "C1QA",
  # 2 Antigen-presenting_Mast_2
  "IGKC", "B2M", "HLA-DRA",
  # 3 Classical_MCTC
  "CPA3", "KIT", "HDC",
  # 4 Epithelial-associated_Mast
  "KRT5", "KRT17", "KRT18",
  # 5 Homeostatic_Mast
  "PFN1", "B2M", "HLA-E",
  # 6 Immunoregulatory_Mast
  "IGKC", "CD74", "APOE",
  # 7 Inflammatory_Mast
  "S100A8", "S100A9", "CXCL8",
  # 8 Matrix_Remodeling_Mast
  "COL1A1", "MMP2", "SPARC",
  # 9 Metabolic_Mast_1
  "NDUFA3", "SNRPD1", "VPS29",
  # 10 Metabolic_Mast_2
  "G6PD", "PGD", "TALDO1",
  # 11 Stress-response_Mast
  "PTEN", "SIGMAR1", "KRT18",
  # 12 Proliferative_Mast (单细胞特有)
  "TOP2A", "MKI67", "AURKB",
  # 13 Tissue-resident_Mast (单细胞特有)
  "SCGB1A1", "APOE", "MBP"
)
all_genes <- unique(all_genes)

library(Seurat)
library(ggplot2)
library(patchwork)

mast$sub_cell_type <- as.character(mast$sub_cell_type)
mast_sc$sub_cell_type <- as.character(mast_sc$sub_cell_type)

p1 <- DotPlot(mast, features = all_genes, group.by = "sub_cell_type", scale = TRUE) +
  scale_color_gradient(low = "#D6EAF8", high = "#1F618D", name = "Average Expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
  labs(title = "Spatial - Mast Cell Subtypes", x = "", y = "") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

p2 <- DotPlot(mast_sc, features = all_genes, group.by = "sub_cell_type", scale = TRUE) +
  scale_color_gradient(low = "#D6EAF8", high = "#1F618D", name = "Average Expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
  labs(title = "scRNA-seq - Mast Cell Subtypes", x = "", y = "") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

p_combined <- p1 / p2 + 
  plot_annotation(title = "Mast Cell Marker Genes: Spatial vs scRNA-seq") &
  theme(plot.title = element_text(hjust = 0.5, size = 14))

ggsave("bubble_plot_mast_optimized.pdf", plot = p_combined, width = 16, height = 14)
```
# CN 
```R
#找到最佳k
bsub -q mpi -n 24 -J Rscript best_k_sub.r
#重新画best-k curve

out_dir <- "./cellcharter_results_autok"
all_stability <- readRDS("./cellcharter_results_autok/autok_combined_results.rds")
pdf(file.path(out_dir, "stability_curve_combined.pdf"), width = 12, height = 6)

# 提取数据
k_vals <- all_stability$all_stability$k
stability_vals <- all_stability$all_stability$stability
best_k <- all_stability$global_best_k

# 计算局部最大值（只针对稳定值 > 0 的区域，避免尾部零值干扰）
# 找到稳定值 > 0 的索引
pos_idx <- which(stability_vals > 0)
if (length(pos_idx) >= 3) {
  stability_pos <- stability_vals[pos_idx]
  k_pos <- k_vals[pos_idx]
  # 计算局部最大值
  peaks_idx <- which(diff(sign(diff(stability_pos))) == -2) + 1
  peaks_k <- k_pos[peaks_idx]
  peaks_stability <- stability_pos[peaks_idx]
} else {
  peaks_k <- NULL
  peaks_stability <- NULL
}

# 绘图
plot(k_vals, stability_vals,
     type = "b",           # 同时画点和线
     xlab = "Number of clusters (k)",
     ylab = "Stability score",
     main = "Cluster stability across k values",
     pch = 1,              # 空心圆点
     col = "black",
     lty = 1,
     lwd = 1.5,
     ylim = c(0, max(stability_vals) * 1.05))

# 标注局部最大值点（实心红点）
if (length(peaks_k) > 0) {
  points(peaks_k, peaks_stability,
         col = "red", pch = 16, cex = 1.5)
}

# 最佳 k 值竖线
abline(v = best_k, col = "blue", lty = 2, lwd = 1.5)



# 图例
legend("topright",
       legend = c("Stability", "Local maxima", "Best k"),
       col = c("black", "red", "blue"),
       pch = c(1, 16, NA),
       lty = c(1, NA, 2),
       bty = "n")

dev.off()
#最佳k聚类划分CN
bsub -q mpi -n 24 -J Rscript run_with_fixed_k.r
#CN画图
#visualize_cellcharter_results.r

#CN组成
# 查看每个CN的细胞组成
cn_composition <- table(mPT@meta.data$cellcharter_cluster, 
                        mPT@meta.data$CellType)

# 计算每个CN的比例
cn_prop <- prop.table(cn_composition, margin = 1) * 100

# 打印每个CN的主要细胞类型（前5种）
for (cn in rownames(cn_prop)) {
  cat("\n=== CN", cn, "===\n")
  top10 <- sort(cn_prop[cn, ], decreasing = TRUE)[1:10]
  for (i in 1:length(top10)) {
    cat(names(top10)[i], ": ", round(top10[i], 1), "%\n", sep = "")
  }
}

#CN注释

cn_annotation_detailed <- c(
  "1" = "Stroma_Immune_Zone",
  "2" = "Metabolic_Tumor_Core_A",
  "3" = "Perivascular_Stroma",
  "4" = "Tumor_Stroma_Interface",
  "5" = "Macrophage_Dominant_Zone",
  "6" = "B_Cell_Development_Zone",
  "7" = "Metabolic_Tumor_Core_B",
  "8" = "GC_B_Tumor_Interface",
  "9" = "Alveolar_Epithelium_Zone",
  "10" = "Fibroblast_Dominant_Stroma",
  "11" = "Secretory_Tumor_Zone",
  "12" = "Metabolic_Tumor_Core_C",
  "13" = "Lymphoid_Aggregate_TLS_like",
  "14" = "B_Plasma_Transition_Zone",
  "15" = "Metabolic_Tumor_Core_D",
  "16" = "Metabolic_GC_Interface",
  "17" = "Metabolic_Squamous_Mix",
  "18" = "Plasma_Cell_Rich_Zone",
  "19" = "Plasma_B_Mixed_Zone"
)
# 大类注释（用于分组统计/可视化）
cn_annotation_major <- c(
  "1" = "Stroma",
  "2" = "Tumor",
  "3" = "Stroma",
  "4" = "Tumor_Stroma",
  "5" = "Myeloid",
  "6" = "Lymphoid",
  "7" = "Tumor",
  "8" = "Lymphoid",
  "9" = "Epithelium",
  "10" = "Stroma",
  "11" = "Tumor",
  "12" = "Tumor",
  "13" = "Lymphoid",
  "14" = "Lymphoid",
  "15" = "Tumor",
  "16" = "Tumor",
  "17" = "Tumor",
  "18" = "Lymphoid",
  "19" = "Lymphoid"
)

# 应用到 Seurat 对象
mPT@meta.data$cn_detailed <- cn_annotation_detailed[as.character(mPT@meta.data$cellcharter_cluster)]
mPT@meta.data$cn_major <- cn_annotation_major[as.character(mPT@meta.data$cellcharter_cluster)]

#shannon
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
    
    # 获取该cluster对应的cn_detailed注释
    cn_name <- unique(cluster_cells$cn_detailed)
    if (length(cn_name) == 0) cn_name <- as.character(cluster)
    
    shannon_diversity <- rbind(shannon_diversity, 
                              data.frame(sample = sample, 
                                        cluster = as.numeric(as.character(cluster)),
                                        cn_detailed = cn_name,
                                        shannon = shannon))
  }
}

# 按照cn_detailed的名称排序（可按需调整顺序）
shannon_diversity$cn_detailed <- factor(shannon_diversity$cn_detailed, 
                                        levels = sort(unique(shannon_diversity$cn_detailed)))

# 获取每个cn_detailed对应的颜色（取第一个cluster的颜色作为代表）
cn_detailed_colors <- sapply(levels(shannon_diversity$cn_detailed), function(cn) {
  # 找到该cn对应的第一个cluster编号
  cluster_num <- unique(shannon_diversity$cluster[shannon_diversity$cn_detailed == cn])[1]
  if (!is.na(cluster_num) && cluster_num <= length(cluster_colors)) {
    return(cluster_colors[cluster_num])
  } else {
    return("gray60")
  }
})

p_shannon <- ggplot(shannon_diversity, aes(x = cn_detailed, y = shannon, fill = cn_detailed)) +
  geom_violin(alpha = 0.7, trim = TRUE, linewidth = 0.5, color = "black") +
  geom_jitter(width = 0.2, height = 0, size = 0.8, alpha = 0.6, color = "black") +
  scale_fill_manual(values = cn_detailed_colors, name = "Spatial Domain") +
  labs(title = "Shannon Diversity Index by Spatial Domain",
       x = "Spatial Domain (CN)",
       y = "Shannon Diversity Index") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),  # 文字较多，稍微减小字号
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.position = "none",
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_line(color = "gray95")
  )

# 保存Shannon多样性图
ggsave("shannon_diversity.pdf", p_shannon, width = 20, height = 6, dpi = 300)  # 宽度适当增加以适应更长的标签
message("保存Shannon多样性图: shannon_diversity.pdf")

#
# 指定样本和要高亮的CN注释
target_sample <- "D3"  # 替换为你的实际样本名称
target_cn <- "Metabolic_GC_Interface"

# 筛选目标样本的细胞
sample_data <- mPT@meta.data[mPT@meta.data$sample == target_sample, ]
sample_cells <- rownames(sample_data)
temp_seurat <- subset(mPT, cells = sample_cells)

# 创建一个高亮分组变量
temp_seurat@meta.data$highlight <- ifelse(
  temp_seurat@meta.data$cn_detailed == target_cn, 
  target_cn, 
  "Other"
)

# 定义颜色：目标CN使用原来的颜色，Other使用灰色
highlight_colors <- c(
  setNames(cluster_colors[which(unique(mPT@meta.data$cn_detailed) == target_cn)[1]], target_cn),
  "Other" = "lightgray"
)

# 绘图
p <- DimPlot(temp_seurat, 
             reduction = reduction_name,
             group.by = "highlight", 
             pt.size = 0.3,
             label = FALSE,
             cols = highlight_colors,
             order = c("Other", target_cn)) +  # 确保目标CN绘制在上层
  ggtitle(paste("Sample:", target_sample, "- Highlight:", target_cn)) +
  theme(legend.position = "bottom")

# 保存
sample_filename <- paste0("spatial_domain_highlight_", target_sample, "_", target_cn, ".pdf")
ggsave(sample_filename, p, width = 10, height = 8, dpi = 300)
message(paste("保存高亮图:", sample_filename))



library(ggplot2)
library(dplyr)

# 计算每个 CN 中每种细胞类型的数量
library(ggplot2)
library(dplyr)
library(scales)

stack_data_count <- mPT@meta.data %>%
  group_by(cn_detailed, CellType) %>%
  summarise(count = n(), .groups = "drop")

# 按首字母顺序排序 CN
cn_order <- sort(unique(stack_data_count$cn_detailed))
stack_data_count$cn_detailed <- factor(stack_data_count$cn_detailed, levels = cn_order)

# 获取所有细胞类型
cell_types_unique <- unique(stack_data_count$CellType)
n_types <- length(cell_types_unique)  # 应该是17

# ============================================
# 浅色系配色方案（17种颜色，足够区分）
# ============================================

# 方案1：使用 Brewer 调色板 Set3（12色）+ Paired 补充（5色）
library(RColorBrewer)
colors_set3 <- brewer.pal(12, "Set3")      # 浅色系，12种
colors_paired <- brewer.pal(5, "Paired")   # 补充5种
cell_colors_light <- c(colors_set3, colors_paired)
names(cell_colors_light) <- cell_types_unique


# 绘制数量堆叠图
p_count <- ggplot(stack_data_count, aes(x = cn_detailed, y = count, fill = CellType)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  scale_fill_manual(values = cell_colors_light, name = "Cell Type") +
  labs(title = "Cell Type Count by CN (Ordered by Name)",
       x = "Cell Neighborhood (CN)",
       y = "Cell Count") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "right",
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.4, "cm")
  )

# 保存
ggsave("cn_celltype_count_alphabetical_light.pdf", p_count, width = 16, height = 8, dpi = 300)



target_sample <- "D3"
target_cn <- "Metabolic_Tumor_Core_A"  # 改这里

# 指定要高亮的细胞类型
highlight_cell_types <- c(
  "Tumor_Metabolic_1",
  "Tumor_Metabolic_2", 
  "Mast_Metabolic_2",
  "Macro_Phagocytic"
  # 去掉 Fib_myCAF
)

sample_data <- mPT@meta.data[mPT@meta.data$sample == target_sample, ]
sample_cells <- rownames(sample_data)
temp_seurat <- subset(mPT, cells = sample_cells)

# 创建着色变量：
# - 目标CN内且属于指定细胞类型：显示具体类型名称
# - 其他所有情况：归为"Other"
temp_seurat@meta.data$color_by <- ifelse(
  temp_seurat@meta.data$cn_detailed == target_cn & 
    as.character(temp_seurat@meta.data$detailed) %in% highlight_cell_types,
  as.character(temp_seurat@meta.data$detailed),  # 显示具体细胞类型
  "Other"  # 其他所有细胞
)

# 为指定的5种细胞类型分配颜色
# 方案1：使用Seurat默认颜色
library(scales)
cell_type_colors <- setNames(
  hue_pal()(length(highlight_cell_types)),
  highlight_cell_types
)


# 构建完整颜色向量
all_colors <- c(cell_type_colors, "Other" = "lightgray")

# 设置绘图顺序：Other在底层，指定细胞类型在上层
order_categories <- c("Other", highlight_cell_types)

# 绘图
p <- DimPlot(temp_seurat, 
             reduction = reduction_name,
             group.by = "color_by", 
             pt.size = 0.3,
             label = FALSE,
             cols = all_colors,
             order = order_categories) +
  ggtitle(paste("Sample:", target_sample, "- CN:", target_cn, 
                "\n(Highlighting specific cell types within CN, others gray)")) +
  theme(legend.position = "right") +
  guides(color = guide_legend(title = "Cell Type", ncol = 1, 
                              override.aes = list(size = 3)))

# 保存
sample_filename <- paste0("spatial_domain_highlight_", target_sample, "_", target_cn, "_specific_types.pdf")
ggsave(sample_filename, p, width = 14, height = 10, dpi = 300)
message(paste("保存高亮图:", sample_filename))
```
# fov 点图
```R
library(ggplot2)
library(dplyr)

# 指定的6种细胞类型
target_cell_types <- c("Malignant cells", "Fibroblasts", "Endothelial cells", 
                       "T cells", "Macrophages", "Dendritic cells")

# 指定的6个 FOV
selected_fovs <- c(152, 153, 154, 156, 157, 158)
fov_data <- obj@meta.data %>% filter(fov %in% selected_fovs)

# 创建分组变量
fov_data$cell_group <- "Other"
for (ct in target_cell_types) {
  fov_data$cell_group[fov_data$CellType == ct] <- ct
}

# 颜色方案（RColorBrewer Set1 经典配色）
cell_colors <- c(
  "Malignant cells" = "#df928e",
  "Fibroblasts" = "#1f78b4",
  "Endothelial cells" = "#B5EFB5",
  "T cells" = "#ffc089",
  "Macrophages" = "#6a3d9a",
  "Dendritic cells" = "#dbc43f",
  "Other" = "#E0E0E0"
)
primcol2 = c('#1f78b4','#ffc089','#B5EFB5','#793c1b',
             '#6a3d9a','#333333','#ffff33',
             '#df928e')
# 只保留出现在数据中的颜色
used_colors <- cell_colors[names(cell_colors) %in% unique(fov_data$cell_group)]
if ("Other" %in% unique(fov_data$cell_group) && !"Other" %in% names(used_colors)) {
  used_colors <- c(used_colors, "Other" = "#E0E0E0")
}


p <- ggplot(fov_data, aes(x = CenterX_global_px, y = CenterY_global_px, color = cell_group)) +
  geom_point(size = 0.01, alpha = 1) +  
  scale_color_manual(values = used_colors, name = "Cell Type") +
  labs(title = "mPT FOVs",
       x = "X Coordinate (px)",
       y = "Y Coordinate (px)") +
  theme_bw() +
  theme(
    aspect.ratio = 1/2,  # 👈 就加这一行！
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    legend.position = "right",
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    # 关键：图例点大小
    legend.key.size = unit(0.8, "cm"),      # 图例键大小
    legend.text = element_text(size = 10)   # 图例文字大小
  ) +
  guides(color = guide_legend(override.aes = list(size = 5)))  # 图例点大小

ggsave("fovs_combined_real_coords_mPT_fov152-158.pdf", p, width = 7, height = 3.5, dpi = 300)

#nmpt
selected_fovs <- c(365,366,367,368,371,372,373,374)

# 筛选数据
fov_data <- obj@meta.data %>% filter(fov %in% selected_fovs)

# 创建分组变量
fov_data$cell_group <- "Other"
for (ct in target_cell_types) {
  fov_data$cell_group[fov_data$CellType == ct] <- ct
}

# 颜色方案（RColorBrewer Set1 经典配色）
cell_colors <- c(
  "Malignant cells" = "#df928e",
  "Fibroblasts" = "#1f78b4",
  "Endothelial cells" = "#B5EFB5",
  "T cells" = "#ffc089",
  "Macrophages" = "#6a3d9a",
  "Dendritic cells" = "#dbc43f",
  "Other" = "#E0E0E0"
)

# 只保留出现在数据中的颜色
used_colors <- cell_colors[names(cell_colors) %in% unique(fov_data$cell_group)]
if ("Other" %in% unique(fov_data$cell_group) && !"Other" %in% names(used_colors)) {
  used_colors <- c(used_colors, "Other" = "#E0E0E0")
}


p <- ggplot(fov_data, aes(x = CenterX_global_px, y = CenterY_global_px, color = cell_group)) +
  geom_point(size = 0.01, alpha = 1) +  
  scale_color_manual(values = used_colors, name = "Cell Type") +
  labs(title = "nmPT FOVs",
       x = "X Coordinate (px)",
       y = "Y Coordinate (px)") +
  theme_bw() +
  theme(
    aspect.ratio = 1/2,  # 👈 就加这一行！
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    legend.position = "right",
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    # 关键：图例点大小
    legend.key.size = unit(0.8, "cm"),      # 图例键大小
    legend.text = element_text(size = 10)   # 图例文字大小
  ) +
  guides(color = guide_legend(override.aes = list(size = 5)))  # 图例点大小

ggsave("fovs_combined_real_coords_nmPT_fov365-374.pdf", p, width = 7, height = 3.5, dpi = 300)

selected_fovs <- c(228,229,233,234)
fov_data <- obj@meta.data %>% filter(fov %in% selected_fovs)

# 创建分组变量
fov_data$cell_group <- "Other"
for (ct in target_cell_types) {
  fov_data$cell_group[fov_data$CellType == ct] <- ct
}

# 颜色方案（RColorBrewer Set1 经典配色）
cell_colors <- c(
  "Malignant cells" = "#df928e",
  "Fibroblasts" = "#1f78b4",
  "Endothelial cells" = "#B5EFB5",
  "T cells" = "#ffc089",
  "Macrophages" = "#6a3d9a",
  "Dendritic cells" = "#dbc43f",
  "Other" = "#E0E0E0"
)
primcol2 = c('#1f78b4','#ffc089','#B5EFB5','#793c1b',
             '#6a3d9a','#333333','#ffff33',
             '#df928e')
# 只保留出现在数据中的颜色
used_colors <- cell_colors[names(cell_colors) %in% unique(fov_data$cell_group)]
if ("Other" %in% unique(fov_data$cell_group) && !"Other" %in% names(used_colors)) {
  used_colors <- c(used_colors, "Other" = "#E0E0E0")
}


p <- ggplot(fov_data, aes(x = CenterX_global_px, y = CenterY_global_px, color = cell_group)) +
  geom_point(size = 0.04, alpha = 1) +  
  scale_color_manual(values = used_colors, name = "Cell Type") +
  labs(title = "metLN FOVs",
       x = "X Coordinate (px)",
       y = "Y Coordinate (px)") +
  theme_bw() +
  theme(
    aspect.ratio = 1/2,  # 👈 就加这一行！
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    legend.position = "right",
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    # 关键：图例点大小
    legend.key.size = unit(0.8, "cm"),      # 图例键大小
    legend.text = element_text(size = 10)   # 图例文字大小
  ) +
  guides(color = guide_legend(override.aes = list(size = 5)))  # 图例点大小

ggsave("fovs_combined_real_coords_metLN_fov228-234.pdf", p, width = 7, height = 3.5, dpi = 300)

#negLN
selected_fovs <- c(303, 304, 305, 306, 307, 308)

# 直接从 Seurat 对象的 metadata 中筛选
fov_data <- obj@meta.data %>%
  filter(!(tissue == "negLN" & CellType %in% c("Malignant cells", "Alveolar type II cells"))) %>%
  filter(fov %in% selected_fovs)

# 创建分组变量
target_cell_types <- c("Malignant cells", "Fibroblasts", "Endothelial cells", 
                       "T cells", "Macrophages", "Dendritic cells")

fov_data$cell_group <- "Other"
for (ct in target_cell_types) {
  fov_data$cell_group[fov_data$CellType == ct] <- ct
}

# 查看结果
table(fov_data$cell_group)

# 颜色方案（RColorBrewer Set1 经典配色）
cell_colors <- c(
  "Malignant cells" = "#df928e",
  "Fibroblasts" = "#1f78b4",
  "Endothelial cells" = "#B5EFB5",
  "T cells" = "#ffc089",
  "Macrophages" = "#6a3d9a",
  "Dendritic cells" = "#dbc43f",
  "Other" = "#E0E0E0"
)

# 只保留出现在数据中的颜色
used_colors <- cell_colors[names(cell_colors) %in% unique(fov_data$cell_group)]
if ("Other" %in% unique(fov_data$cell_group) && !"Other" %in% names(used_colors)) {
  used_colors <- c(used_colors, "Other" = "#E0E0E0")
}


p <- ggplot(fov_data, aes(x = CenterX_global_px, y = CenterY_global_px, color = cell_group)) +
  geom_point(size = 0.04, alpha = 1) +  
  scale_color_manual(values = used_colors, name = "Cell Type") +
  labs(title = "negLN FOVs",
       x = "X Coordinate (px)",
       y = "Y Coordinate (px)") +
  theme_bw() +
  theme(
    aspect.ratio = 1/2,  # 👈 就加这一行！
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    legend.position = "right",
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    # 关键：图例点大小
    legend.key.size = unit(0.8, "cm"),      # 图例键大小
    legend.text = element_text(size = 10)   # 图例文字大小
  ) +
  guides(color = guide_legend(override.aes = list(size = 5)))  # 图例点大小

ggsave("fovs_combined_real_coords_negLN_fov303-308.pdf", p, width = 7, height = 3.5, dpi = 300)
```
# 肿瘤亚型的空间分布
```R
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ggsai)

# ===================== 提前定义恶性细胞亚型列表 =====================
malignant_subtypes <- c(
  "Tumor_Immune-inflamed_1", "Tumor_Immune-inflamed_2",
  "Tumor_Metabolic_2", "Tumor_Metabolic_1",
  "Tumor_Immunogenic", "Tumor_Stem-like",
  "Tumor_EMT-like", "Tumor_Secretory", "Tumor_Proliferative"
)

# ===================== 第一步：数据准备 + 统一添加分面变量 =====================
spatial_data <- obj@meta.data %>%
  select(
    CenterX_global_px, CenterY_global_px,
    detailed, sample, fov, tissue,
    Area.um2, Circularity, AspectRatio, NucArea,
    Mean.PanCK, Mean.CD45, Mean.DAPI
  ) %>%
  filter(!is.na(CenterX_global_px), !is.na(CenterY_global_px)) %>%
  # 1. 筛选目标tissue
  filter(tissue %in% c("metLN", "negLN", "mPT", "nmPT")) %>%
  filter(!is.na(tissue)) %>%
  # 2. 标注cell_category（只保留恶性细胞亚型，其他归为背景）
  mutate(
    cell_category = case_when(
      detailed %in% malignant_subtypes ~ detailed,
      TRUE ~ "Other cells"
    ),
    # 3. 创建tissue_sample分面变量
    tissue = factor(tissue, levels = c("metLN", "negLN", "mPT", "nmPT")),
    tissue_sample = paste(tissue, sample, sep = "_")
  )

# 定义颜色（使用 NPG 调色板）
n_tumor <- length(malignant_subtypes)
npg_pal <- pal_npg()(10)  # NPG 默认有10种颜色
tumor_colors <- colorRampPalette(npg_pal)(n_tumor)
names(tumor_colors) <- malignant_subtypes

# 其他细胞颜色（背景）
other_color <- c("Other cells" = "#DCDCDC")

# 合并所有颜色
all_colors <- c(tumor_colors, other_color)

# 调整图例顺序
legend_order <- c(malignant_subtypes, "Other cells")

# ===================== 第二步：分组1绘图（metLN, negLN） =====================
group1_tissues <- c("metLN", "negLN")
group1_data <- spatial_data %>% filter(tissue %in% group1_tissues)

# 检查分组1数据
if(nrow(group1_data) == 0) {
  stop("分组1（metLN/negLN）无数据！")
}

# 按tissue_sample排序
group1_data$tissue_sample <- factor(
  group1_data$tissue_sample,
  levels = unique(group1_data$tissue_sample[order(group1_data$tissue, group1_data$sample)])
)

# 计算分面行数（每行3个sample）
total_facets1 <- length(unique(group1_data$tissue_sample))
n_rows1 <- ceiling(total_facets1 / 3)

# 绘图
p_group1 <- ggplot(group1_data,
                   aes(x = CenterX_global_px, y = CenterY_global_px, 
                       color = cell_category)) +
  # 第1层：Other cells（灰色背景）
  geom_point(data = filter(group1_data, cell_category == "Other cells"),
             alpha = 0.2, size = 0.3) +
  # 第2层：恶性细胞亚型（彩色）
  geom_point(data = filter(group1_data, cell_category %in% malignant_subtypes),
             alpha = 0.9, size = 0.6) +
  scale_color_manual(values = all_colors, breaks = legend_order) +
  facet_wrap(~tissue_sample, scales = "free", ncol = 3, nrow = n_rows1) +
  theme_minimal() +
  labs(title = "Malignant Cell Subtypes Distribution - Lymph Nodes",
       x = "X Coordinate (um)",
       y = "Y Coordinate (um)",
       color = "Cell Type") +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    strip.text = element_text(size = 10, face = "bold"),
    strip.background = element_rect(fill = "gray90"),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9, face = "bold"),
    axis.text = element_text(size = 6),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_line(color = "gray95")
  )

# 保存分组1PDF
pdf("malignant_subtypes_LN.pdf", width = 15, height = 15)
print(p_group1)
dev.off()

# ===================== 第三步：分组2绘图（mPT, nmPT） =====================
group2_tissues <- c("mPT", "nmPT")
group2_data <- spatial_data %>% filter(tissue %in% group2_tissues)

if(nrow(group2_data) > 0) {
  group2_data$tissue_sample <- factor(
    group2_data$tissue_sample,
    levels = unique(group2_data$tissue_sample[order(group2_data$tissue, group2_data$sample)])
  )
  
  total_facets2 <- length(unique(group2_data$tissue_sample))
  n_rows2 <- ceiling(total_facets2 / 3)
  
  p_group2 <- ggplot(group2_data,
                     aes(x = CenterX_global_px, y = CenterY_global_px, 
                         color = cell_category)) +
    geom_point(data = filter(group2_data, cell_category == "Other cells"),
               alpha = 0.2, size = 0.3) +
    geom_point(data = filter(group2_data, cell_category %in% malignant_subtypes),
               alpha = 0.9, size = 0.6) +
    scale_color_manual(values = all_colors, breaks = legend_order) +
    facet_wrap(~tissue_sample, scales = "free", ncol = 3, nrow = n_rows2) +
    theme_minimal() +
    labs(title = "Malignant Cell Subtypes Distribution - Primary Tumor",
         x = "X Coordinate (um)",
         y = "Y Coordinate (um)",
         color = "Cell Type") +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      strip.text = element_text(size = 10, face = "bold"),
      strip.background = element_rect(fill = "gray90"),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 9, face = "bold"),
      axis.text = element_text(size = 6),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_line(color = "gray95")
    )
  
  pdf("malignant_subtypes_PT.pdf", width = 15, height = 15)
  print(p_group2)
  dev.off()
}
```
# 丰度
```R
library(ggplot2)

# 筛选目标细胞类型，例如 Tumor_Metabolic_1
target_cell <- "Tumor_Immune-inflamed_1"
data_subset <- obj@meta.data %>%
  filter(sample == "C2", detailed == target_cell)

# 使用 stat_density2d 绘制密度热图
p <- ggplot(obj@meta.data %>% filter(sample == "C2"), 
            aes(x = CenterX_global_px, y = CenterY_global_px)) +
  # 背景点（所有细胞浅灰色）
  geom_point(alpha = 0.1, size = 0.1, color = "lightgray") +
  # 目标细胞密度等高线 + 填充
  stat_density2d(data = data_subset,
                 aes(fill = after_stat(density)), 
                 geom = "raster", 
                 contour = FALSE,
                 alpha = 0.7) +
  scale_fill_viridis_c(option = "plasma", name = "Density") +
  coord_fixed() +
  labs(title = paste("Spatial density of", target_cell),
       x = "X Coordinate", y = "Y Coordinate") +
  theme_minimal()

ggsave(paste0("density_", target_cell, ".pdf"), p, width = 8, height = 6)
```
# fov 肿瘤亚型
```R

library(ggplot2)
library(dplyr)

# 指定的肿瘤亚型（根据你的 detailed 列中的实际名称）
target_tumor_subtypes <- c("Tumor_Metabolic_1", "Tumor_Metabolic_2", 
                           "Tumor_Immune-inflamed_1", "Tumor_Immune-inflamed_2",
                           "Tumor_Immunogenic", "Tumor_Stem-like",
                           "Tumor_EMT-like", "Tumor_Secretory", "Tumor_Proliferative")

# 指定的 FOV
selected_fovs <- c(257, 258)

# 筛选数据：C2 样本 + 指定 FOV
fov_data <- obj@meta.data %>%
  filter(sample == "C2", fov %in% selected_fovs)

# 创建分组变量：肿瘤亚型保留具体名称，其他细胞归为 "Other"
fov_data$cell_group <- "Other"
for (ct in target_tumor_subtypes) {
  fov_data$cell_group[fov_data$detailed == ct] <- ct
}

# 颜色方案
tumor_colors <- c(
  "Tumor_Metabolic_1" = "#E41A1C",
  "Tumor_Metabolic_2" = "#B22222",
  "Tumor_Immune-inflamed_1" = "#E69F00",
  "Tumor_Immune-inflamed_2" = "#D55E00",
  "Tumor_Immunogenic" = "#4DAF4A",
  "Tumor_Stem-like" = "#377EB8",
  "Tumor_EMT-like" = "#984EA3",
  "Tumor_Secretory" = "#FFD700",
  "Tumor_Proliferative" = "#A65628"
)

other_color <- c("Other" = "#E0E0E0")
cell_colors <- c(tumor_colors, other_color)

# 只保留出现在数据中的颜色
used_colors <- cell_colors[names(cell_colors) %in% unique(fov_data$cell_group)]
if ("Other" %in% unique(fov_data$cell_group) && !"Other" %in% names(used_colors)) {
  used_colors <- c(used_colors, "Other" = "#E0E0E0")
}

# 设置因子水平，确保 Other 在底层（先绘制），肿瘤细胞在上层
fov_data$cell_group <- factor(fov_data$cell_group, 
                               levels = c("Other", target_tumor_subtypes))

# 绘图（分层绘制：先画 Other，再画肿瘤细胞）
p <- ggplot() +
  # 第1层：Other 细胞（底层）
  geom_point(data = fov_data %>% filter(cell_group == "Other"),
             aes(x = CenterX_global_px, y = CenterY_global_px),
             color = "#E0E0E0", size = 0.01, alpha = 0.5) +
  # 第2层：肿瘤亚型（上层）
  geom_point(data = fov_data %>% filter(cell_group != "Other"),
             aes(x = CenterX_global_px, y = CenterY_global_px, color = cell_group),
             size = 0.01, alpha = 0.9) +
  scale_color_manual(values = used_colors, name = "Cell Type") +
  labs(title = "mPT (C2) - FOVs 257,258: Tumor Subtypes",
       x = "X Coordinate (px)",
       y = "Y Coordinate (px)") +
  theme_bw() +
  theme(
    aspect.ratio = 1/2,
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    legend.position = "right",
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    legend.key.size = unit(0.6, "cm"),
    legend.text = element_text(size = 8)
  ) +
  guides(color = guide_legend(override.aes = list(size = 3), ncol = 1))

# 保存
ggsave("mPT_C2_fov257_258_tumor_subtypes.pdf", p, width = 7, height = 3.5, dpi = 300)
```
# DC
```R
options(future.globals.maxSize = 8000 * 1024^2)
library(Seurat)
library(ggplot2)
library(dplyr)
obj <- readRDS("YA2025263-1_fin.rds")
DC <- subset(obj,subset=CellType %in% c("Dendritic cells", "Plasmacytoid dendritic cells"))
#DC <- SCTransform(DC,
                  #vars.to.regress = "sample",
                  #variable.features.n = 3000,
                  #conserve.memory = TRUE,
                  #verbose = TRUE)
DC <- NormalizeData(DC)
DC <- FindVariableFeatures(DC, nfeatures = 3000)
hvgs <- VariableFeatures(DC)
DC <- ScaleData(DC, features = hvgs)
DC <- RunPCA(DC, features = hvgs, npcs = 50)

#DC <- FindVariableFeatures(DC, selection.method = "vst", nfeatures = 4000)
#hvgs <- VariableFeatures(DC)

# 3. 关键：Scale 时 去除 线粒体、核糖体、细胞周期 干扰
# 这是 DC 能分开的核心！
#DC <- ScaleData(
  DC,
  features = hvgs,
  vars.to.regress = c("nCount_RNA"), # 去批次/去污染
  do.center = T, do.scale = T
)



# 5. 用 DC 特征基因跑 PCA，而不是全部基因
#DC <- RunPCA(
  DC,
  npcs = 20,              # DC 亚型少，不需要50
  verbose = F
)


p <- ElbowPlot(DC, ndims = 30)
ggsave("p3.png",plot=p)
pca_var <- DC[["pca"]]@stdev^2 / sum(DC[["pca"]]@stdev^2)
cumsum_var <- cumsum(pca_var)

# 看达到多少 PC 时累计方差 > 70% / 80% / 90%
which(cumsum_var > 0.7)[1]   # 70%
which(cumsum_var > 0.8)[1]   # 80%
which(cumsum_var > 0.9)[1]   # 90%
#library(harmony)
#DC <- RunHarmony(DC,"tissue")

DC <- FindNeighbors(DC, dims = 1:31)
DC <- FindClusters(DC, resolution = 0.4)
DC <- RunUMAP(DC, dims = 1:31)

p <- DimPlot(DC, reduction = "umap", group.by = "sample", label = FALSE)
ggsave("DC_orig.png", plot = p)

library(ggplot2)
library(ggsci)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(18)  # 0.6

seurat_clusters <- as.character(unique(DC@meta.data$seurat_clusters))  # 转换为字符向量
num_legend_items <- length(seurat_clusters)  # 图例的个数
max_label_length <- max(nchar(seurat_clusters))  # 图例名称的最大长度

base_width <- 3000  # 基础宽度
base_height <- 3000  # 基础高度
legend_width_factor <- 100  # 每个图例项增加的宽度
label_length_factor <- 10  # 每个字符增加的宽度
dynamic_width <- base_width + (num_legend_items * legend_width_factor) + (max_label_length * label_length_factor)
pdf("DC_clusters.pdf", width = dynamic_width / 300, height = base_height / 300)  # 转换为英寸

DimPlot(DC, reduction = "umap", label = TRUE, pt.size = 1, label.size = 8) +
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    ggtitle(NULL) +
    scale_color_manual(values = npg_extended) +
    coord_fixed(ratio = 1) +
    guides(color = guide_legend(title = NULL, override.aes = list(size = 5))) +
    theme(
        text = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 28, color = "black"),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title.x = element_text(size = 36, face = "bold", color = "black"),
        axis.title.y = element_text(size = 36, face = "bold", color = "black", margin = margin(r = 20)),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 28, face = "bold", color = "black"),
        legend.title = element_text(size = 28, face = "bold", color = "black"),
        legend.position = "right",
        legend.box.margin = margin(0, 0, 0, 0),
        legend.key = element_blank(),
        legend.background = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5),
        aspect.ratio = 1,
        plot.margin = margin(10, 50, 10, 10)
    )

dev.off()

library(dplyr)
DC_markers <- FindAllMarkers(DC, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "wilcox")
DC_significant_markers <- subset(DC_markers, p_val_adj < 0.05)
#write.csv(DC_significant_markers, "DC_all_marker.csv")
DC_significant_markers <- DC_significant_markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
write.csv(DC_significant_markers, "DC_top_marker_50.csv")

DC <- subset(DC,subset=seurat_clusters %in% c(0,1,2,3,4,5,6))
DC <- RunUMAP(DC, dims = 1:31)
library(ggplot2)
library(pheatmap)
library(dplyr)

# ==================== 定义 DC 亚型基因集 ====================

# 1. cDC1（1型常规树突状细胞）
cdc1_genes <- c(
  "CLEC9A", "XCR1", "BATF3", "IRF8", "IDO1",
  "CADM1", "DNASE1L3", "WDFY4", "ITGAE", "CD103",
  "TLR3", "TLR8", "IFNAR1", "IFNAR2", "STAT1",
  "STAT2", "IRF1", "IRF7", "BATF2", "ZBTB46"
)

# 2. cDC2（2型常规树突状细胞）
cdc2_genes <- c(
  "CD1C", "CD1B", "CLEC10A", "FCER1A", "ITGAM",
  "CD11B", "SIRPA", "CD172A", "FCGR2A", "FCGR2B",
  "CLEC6A", "CLEC7A", "CLEC4A", "CLEC4C", "CD209",
  "IL1R2", "IL1R1", "CD2", "CD5", "CD40"
)

# 3. pDC（浆细胞样树突状细胞）
pdc_genes <- c(
  "LILRA4", "IRF7", "TCF4", "CLEC4C", "IL3RA",
  "CD123", "NRP1", "CD303", "CD304", "TLR7",
  "TLR9", "BCL11A", "SPIB", "SPI1", "IRF8",
  "PDCD1LG2", "CD300C", "CD86", "CD80", "CD40"
)

# 4. 成熟 DC（mDC）
mature_dc_genes <- c(
  "CCR7", "LAMP3", "CD83", "CD80", "CD86",
  "CD40", "CD70", "CD40LG", "TNFRSF4", "TNFRSF9",
  "CD274", "PDCD1LG2", "IDO1", "CCL19", "CCL22",
  "CCL17", "CXCL10", "CXCL11", "IL6", "IL12A"
)

# 5. 迁移性 DC（Migratory DC）
migratory_dc_genes <- c(
  "CCR7", "CCL19", "CCL21", "CCL22", "CCL17",
  "CD80", "CD86", "CD83", "LAMP3", "CD40",
  "CD44", "MMP2", "MMP9", "S1PR1", "S1PR3",
  "CXCR4", "CXCR5", "GPR183", "LTB", "LTBR"
)

# 6. 朗格汉斯细胞（Langerhans cells, LC）
langerhans_genes <- c(
  "CD207", "CD1A", "CD1C", "EPCAM", "KRT1",
  "KRT10", "LANGERIN", "CD11B", "ITGAM", "CD11C",
  "CD86", "CD80", "CD40", "CD274", "CXCR4",
  "CCR6", "CCR7", "CD24", "CD44", "CD68"
)

# 7. 炎症性 DC（Inflammatory DC, infDC）
inflammatory_dc_genes <- c(
  "S100A8", "S100A9", "S100A12", "TNF", "IL1B",
  "IL6", "CXCL8", "CCL2", "CCL3", "CCL4",
  "CXCL1", "CXCL2", "CXCL3", "CCL20", "CCL5",
  "NFKBIA", "NFKB1", "REL", "RELA", "STAT3"
)

# 8. 泛 DC 标志（pan-DC）
pan_dc_genes <- c(
  "FLT3", "ZBTB46", "IRF8", "ITGAX", "CD11C",
  "HLA-DRA", "HLA-DRB", "CD86", "CD83", "CD40"
)

# ==================== 过滤存在的基因 ====================
cdc1_genes <- cdc1_genes[cdc1_genes %in% rownames(DC)]
cdc2_genes <- cdc2_genes[cdc2_genes %in% rownames(DC)]
pdc_genes <- pdc_genes[pdc_genes %in% rownames(DC)]
mature_dc_genes <- mature_dc_genes[mature_dc_genes %in% rownames(DC)]
migratory_dc_genes <- migratory_dc_genes[migratory_dc_genes %in% rownames(DC)]
langerhans_genes <- langerhans_genes[langerhans_genes %in% rownames(DC)]
inflammatory_dc_genes <- inflammatory_dc_genes[inflammatory_dc_genes %in% rownames(DC)]

cat("cDC1 基因数:", length(cdc1_genes), "\n")
cat("cDC2 基因数:", length(cdc2_genes), "\n")
cat("pDC 基因数:", length(pdc_genes), "\n")
cat("成熟 DC 基因数:", length(mature_dc_genes), "\n")
cat("迁移性 DC 基因数:", length(migratory_dc_genes), "\n")
cat("朗格汉斯细胞 基因数:", length(langerhans_genes), "\n")
cat("炎症性 DC 基因数:", length(inflammatory_dc_genes), "\n")

# ==================== 计算评分 ====================
DC <- AddModuleScore(DC, features = list(cdc1_genes), name = "cDC1_Score")
DC <- AddModuleScore(DC, features = list(cdc2_genes), name = "cDC2_Score")
DC <- AddModuleScore(DC, features = list(pdc_genes), name = "pDC_Score")
DC <- AddModuleScore(DC, features = list(mature_dc_genes), name = "Mature_Score")
DC <- AddModuleScore(DC, features = list(migratory_dc_genes), name = "Migratory_Score")
DC <- AddModuleScore(DC, features = list(langerhans_genes), name = "LC_Score")
DC <- AddModuleScore(DC, features = list(inflammatory_dc_genes), name = "InfDC_Score")

# ==================== 计算每个 cluster 的平均评分 ====================
score_cols <- c("cDC1_Score1", "cDC2_Score1", "pDC_Score1", 
                "Mature_Score1", "LC_Score1", "InfDC_Score1")

score_avg <- DC@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(across(all_of(score_cols), mean, na.rm = TRUE))

# 转换为矩阵
score_mat <- as.matrix(score_avg[, -1])
rownames(score_mat) <- score_avg$seurat_clusters

# ==================== 绘制热图 ====================
# 按行缩放（Z-score）
pheatmap(score_mat,
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         main = "DC subset scores by cluster",
         fontsize = 10,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         filename = "DC_subset_scores_heatmap.pdf",
         width = 10, height = 8)

# 不缩放（原始值）
pheatmap(score_mat,
         scale = "none",
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         main = "DC subset scores by cluster (raw values)",
         fontsize = 10,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         filename = "DC_subset_scores_raw.pdf",
         width = 10, height = 8)



# ==================== 打印评分表 ====================
print(score_avg)

# ==================== 根据最高分自动命名建议 ====================
score_cols_short <- c("cDC1", "cDC2", "pDC", "Mature",  "LC", "InfDC")
colnames(score_mat) <- score_cols_short

# 找出每个 cluster 的最高分亚型
max_score <- apply(score_mat, 1, function(x) colnames(score_mat)[which.max(x)])
max_value <- apply(score_mat, 1, max)

cat("\n===== 命名建议 =====\n")
for (i in seq_along(max_score)) {
  cat("Cluster", rownames(score_mat)[i], ":", max_score[i], "(评分:", round(max_value[i], 3), ")\n")
}

identity_mapping <- c(
  "0"= "DC_cDC2_1",
  "1" = "DC_Mature",
  "2"= "DC_cDC2_2",
  "3" = "DC_pDC",
  "4" = "DC_Mature",    
  "5" = "DC_Inflammatory",
  "6" = "DC_cDC1"
)
DC$seurat_clusters <- factor(DC$seurat_clusters)
sub_cell_type <- identity_mapping[DC@meta.data$seurat_clusters]
DC@meta.data$sub_cell_type <- sub_cell_type

library(ggplot2)
library(ggsci)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(5)
cluster_colors <- c(
  # 浅色系在前（浅蓝第一个）
  "#A0CBE8",   # 1 浅蓝（你要的放第一个）
  "#FFBE7D",   # 2 浅橙
  "#8CD17D",   # 3 浅绿
  "#86BCB6",   # 4 浅青
  "#FF9D9A",   # 5 浅红
  "#FABFD2",   # 6 浅粉
  "#D4A6C8",   # 7 浅紫
  "#F1CE63",   # 8 淡黄
  "#D7B5A6",   # 9 浅棕
  "#B07AA1",   # 10 紫灰
  "#BAB0AC",   # 11 浅灰
  
  # 深色系在后
  "#4E79A7",   # 12 深蓝
  "#F28E2B",   # 13 橙
  "#59A14F",   # 14 深绿
  "#E15759",   # 15 红
  "#499894",   # 16 青绿
  "#D37295",   # 17 玫红
  "#B6992D",   # 18 土黄
  "#9D7660",   # 19 棕
  "#79706E"    # 20 深灰
)
cell_types <- as.character(unique(DC@meta.data$sub_cell_type))
num_legend_items <- length(cell_types)  # 图例的个数
max_label_length <- max(nchar(cell_types))  # 图例名称的最大长度

base_width <- 3000  # 基础宽度
base_height <- 3000  # 基础高度
legend_width_factor <- 100  # 每个图例项增加的宽度
label_length_factor <- 10  # 每个字符增加的宽度
dynamic_width <- base_width + (num_legend_items * legend_width_factor) + (max_label_length * label_length_factor)
pdf("DC_annotation.pdf", width = dynamic_width/300, height = base_height/300)
DimPlot(DC, reduction = "umap", label = TRUE, pt.size = 1, group.by = "sub_cell_type", label.size = 4) +
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    ggtitle(NULL) +
    scale_color_manual(values = cluster_colors ) +
    coord_fixed(ratio = 1) +
    guides(color = guide_legend(title = NULL, override.aes = list(size = 5))) +
    theme(
        text = element_text(size = 8, face = "bold"),
        axis.text.x = element_text(size = 28, color = "black"),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title.x = element_text(size = 36, face = "bold", color = "black"),
        axis.title.y = element_text(size = 36, face = "bold", color = "black", margin = margin(r = 20)),  # 增加右侧间距
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 16, face = "bold", color = "black"),
        legend.title = element_text(size = 16, face = "bold", color = "black"),
        legend.position = "right",
        legend.box.margin = margin(0, 0, 0, 0),
        legend.key = element_blank(),
        legend.background = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5),
        aspect.ratio = 1,
        plot.margin = margin(10, 50, 10, 10)
    )
dev.off()
```
# DC sc
```R
obj <- readRDS("sc_obj_anno.rds")

library(Seurat)
library(ggplot2)
library(dplyr)
DC <- subset(obj,subset=cell_type %in% c("DC"))
DC <- NormalizeData(DC)
DC <- FindVariableFeatures(DC, nfeatures = 2000)
hvgs <- VariableFeatures(DC)
DC <- ScaleData(DC, features = hvgs)
DC <- RunPCA(DC, features = hvgs, npcs = 50)

library(ggplot2)
p <- ElbowPlot(DC, ndims = 30)
ggsave("p3.png", plot = p)

pca_var <- DC[["pca"]]@stdev^2 / sum(DC[["pca"]]@stdev^2)
cumsum_var <- cumsum(pca_var)

# 看达到多少 PC 时累计方差 > 70% / 80% / 90%
which(cumsum_var > 0.7)[1]   # 70%
which(cumsum_var > 0.8)[1]   # 80%
which(cumsum_var > 0.9)[1]   # 90%

library(harmony)
DC <- RunHarmony(DC, "dataset")
DC <- FindNeighbors(DC, reduction = "harmony", dims = 1:20)
DC <- FindClusters(DC, resolution = 0.2)
DC <- RunUMAP(DC, reduction = "harmony", dims = 1:20)

library(ggplot2)
library(ggsci)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(17)

seurat_clusters <- as.character(unique(DC@meta.data$seurat_clusters))  # 转换为字符向量
num_legend_items <- length(seurat_clusters)  # 图例的个数
max_label_length <- max(nchar(seurat_clusters))  # 图例名称的最大长度

base_width <- 3000  # 基础宽度
base_height <- 3000  # 基础高度
legend_width_factor <- 100  # 每个图例项增加的宽度
label_length_factor <- 10  # 每个字符增加的宽度
dynamic_width <- base_width + (num_legend_items * legend_width_factor) + (max_label_length * label_length_factor)
pdf("DC_sc_clusters.pdf", width = dynamic_width / 300, height = base_height / 300)  # 转换为英寸

DimPlot(DC, reduction = "umap", label = TRUE, pt.size = 1, label.size = 8) +
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    ggtitle(NULL) +
    scale_color_manual(values = npg_extended) +
    coord_fixed(ratio = 1) +
    guides(color = guide_legend(title = NULL, override.aes = list(size = 5))) +
    theme(
        text = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 28, color = "black"),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title.x = element_text(size = 36, face = "bold", color = "black"),
        axis.title.y = element_text(size = 36, face = "bold", color = "black", margin = margin(r = 20)),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 28, face = "bold", color = "black"),
        legend.title = element_text(size = 28, face = "bold", color = "black"),
        legend.position = "right",
        legend.box.margin = margin(0, 0, 0, 0),
        legend.key = element_blank(),
        legend.background = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5),
        aspect.ratio = 1,
        plot.margin = margin(10, 50, 10, 10)
    )

dev.off()

p <- DimPlot(DC, reduction = "umap", group.by = "dataset", label = FALSE)
ggsave("DC_orig.png", plot = p)

library(dplyr)
DC_markers <- FindAllMarkers(DC, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "wilcox")
DC_significant_markers <- subset(DC_markers, p_val_adj < 0.05)
#write.csv(DC_significant_markers, "DC_all_marker.csv")
DC_significant_markers <- DC_significant_markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
write.csv(DC_significant_markers, "DC_sc_top_marker_50.csv")

pheatmap(score_mat,
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         main = "DC subset scores by cluster",
         fontsize = 10,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         filename = "DC_sc_subset_scores_heatmap.pdf",
         width = 10, height = 8)

# 不缩放（原始值）
pheatmap(score_mat,
         scale = "none",
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         main = "DC subset scores by cluster (raw values)",
         fontsize = 10,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         filename = "DC_sc_subset_scores_raw.pdf",
         width = 10, height = 8)

DC <- subset(DC,subset=seurat_clusters %in% c(0,1,3,5,7,9))
DC <- subset(DC,subset=seurat_clusters %in% c(0,1,3,5,7,9,10,11,12,13))
DC <- RunUMAP(DC, reduction = "harmony",dims = 1:20)

identity_mapping <- c(
  "0" = "DC_Inflammatory",
  "1" = "DC_cDC2",
  "3" = "DC_pDC",
  "5" = "DC_Langerhans",
  "7" = "DC_Mature",
  "9" = "DC_cDC1",
  "10" = "DC_Inflammatory",
  "11" = "DC_Inflammatory",
  "12" = "DC_pDC",
  "13" = "DC_Inflammatory"
)

identity_mapping <- c(
  "0" = "DC_Inflammatory",
  "1" = "DC_cDC2",
  "3" = "DC_pDC",
  "5" = "DC_Langerhans",
  "7" = "DC_Mature",
  "9" = "DC_cDC1"
)
DC$seurat_clusters <- factor(DC$seurat_clusters)
sub_cell_type <- identity_mapping[DC@meta.data$seurat_clusters]
DC@meta.data$sub_cell_type <- sub_cell_type

library(ggplot2)
library(ggsci)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(5)
cluster_colors <- c(
  # 浅色系在前（浅蓝第一个）
  "#A0CBE8",   # 1 浅蓝（你要的放第一个）
  "#FFBE7D",   # 2 浅橙
  "#8CD17D",   # 3 浅绿
  "#86BCB6",   # 4 浅青
  "#FF9D9A",   # 5 浅红
  "#FABFD2",   # 6 浅粉
  "#D4A6C8",   # 7 浅紫
  "#F1CE63",   # 8 淡黄
  "#D7B5A6",   # 9 浅棕
  "#B07AA1",   # 10 紫灰
  "#BAB0AC",   # 11 浅灰
  
  # 深色系在后
  "#4E79A7",   # 12 深蓝
  "#F28E2B",   # 13 橙
  "#59A14F",   # 14 深绿
  "#E15759",   # 15 红
  "#499894",   # 16 青绿
  "#D37295",   # 17 玫红
  "#B6992D",   # 18 土黄
  "#9D7660",   # 19 棕
  "#79706E"    # 20 深灰
)
cell_types <- as.character(unique(DC@meta.data$sub_cell_type))
num_legend_items <- length(cell_types)  # 图例的个数
max_label_length <- max(nchar(cell_types))  # 图例名称的最大长度

base_width <- 3000  # 基础宽度
base_height <- 3000  # 基础高度
legend_width_factor <- 100  # 每个图例项增加的宽度
label_length_factor <- 10  # 每个字符增加的宽度
dynamic_width <- base_width + (num_legend_items * legend_width_factor) + (max_label_length * label_length_factor)
pdf("DC_sc_annotation.pdf", width = dynamic_width/300, height = base_height/300)
DimPlot(DC, reduction = "umap", label = TRUE, pt.size = 1, group.by = "sub_cell_type", label.size = 4) +
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    ggtitle(NULL) +
    scale_color_manual(values = cluster_colors ) +
    coord_fixed(ratio = 1) +
    guides(color = guide_legend(title = NULL, override.aes = list(size = 5))) +
    theme(
        text = element_text(size = 8, face = "bold"),
        axis.text.x = element_text(size = 28, color = "black"),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title.x = element_text(size = 36, face = "bold", color = "black"),
        axis.title.y = element_text(size = 36, face = "bold", color = "black", margin = margin(r = 20)),  # 增加右侧间距
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 16, face = "bold", color = "black"),
        legend.title = element_text(size = 16, face = "bold", color = "black"),
        legend.position = "right",
        legend.box.margin = margin(0, 0, 0, 0),
        legend.key = element_blank(),
        legend.background = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5),
        aspect.ratio = 1,
        plot.margin = margin(10, 50, 10, 10)
    )
dev.off()
```
# B 
```R
options(future.globals.maxSize = 8000 * 1024^2)
library(Seurat)
library(ggplot2)
library(dplyr)
obj <- readRDS("YA2025263-1_fin.rds")
B <- subset(obj,subset=CellType %in% c("B cells", "Germinal center B cells", "immature B cells"))
B <- NormalizeData(B)
B <- FindVariableFeatures(B, nfeatures = 2000)
hvgs <- VariableFeatures(B)
B <- ScaleData(B, features = hvgs)
B <- RunPCA(B, features = hvgs, npcs = 50)

p <- ElbowPlot(B, ndims = 30)
ggsave("p3.png", plot = p)
pca_var <- B[["pca"]]@stdev^2 / sum(B[["pca"]]@stdev^2)
cumsum_var <- cumsum(pca_var)

# 看达到多少 PC 时累计方差 > 70% / 80% / 90%
which(cumsum_var > 0.7)[1]   # 70%
which(cumsum_var > 0.8)[1]   # 80%
which(cumsum_var > 0.9)[1]   # 90%

B <- FindNeighbors(B, dims = 1:30)
B <- FindClusters(B, resolution = 0.3)
B <- RunUMAP(B, dims = 1:30)

library(ggplot2)
library(ggsci)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(17)

seurat_clusters <- as.character(unique(B@meta.data$seurat_clusters))  # 转换为字符向量
num_legend_items <- length(seurat_clusters)  # 图例的个数
max_label_length <- max(nchar(seurat_clusters))  # 图例名称的最大长度

base_width <- 3000  # 基础宽度
base_height <- 3000  # 基础高度
legend_width_factor <- 100  # 每个图例项增加的宽度
label_length_factor <- 10  # 每个字符增加的宽度
dynamic_width <- base_width + (num_legend_items * legend_width_factor) + (max_label_length * label_length_factor)
pdf("B_clusters.pdf", width = dynamic_width / 300, height = base_height / 300)  # 转换为英寸

DimPlot(B, reduction = "umap", label = TRUE, pt.size = 1, label.size = 8) +
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    ggtitle(NULL) +
    scale_color_manual(values = npg_extended) +
    coord_fixed(ratio = 1) +
    guides(color = guide_legend(title = NULL, override.aes = list(size = 5))) +
    theme(
        text = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 28, color = "black"),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title.x = element_text(size = 36, face = "bold", color = "black"),
        axis.title.y = element_text(size = 36, face = "bold", color = "black", margin = margin(r = 20)),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 28, face = "bold", color = "black"),
        legend.title = element_text(size = 28, face = "bold", color = "black"),
        legend.position = "right",
        legend.box.margin = margin(0, 0, 0, 0),
        legend.key = element_blank(),
        legend.background = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5),
        aspect.ratio = 1,
        plot.margin = margin(10, 50, 10, 10)
    )

dev.off()

p <- DimPlot(B, reduction = "umap", group.by = "dataset", label = FALSE)
ggsave("B_orig.png", plot = p)

library(dplyr)
B_markers <- FindAllMarkers(B, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "wilcox")
B_significant_markers <- subset(B_markers, p_val_adj < 0.05)
#write.csv(B_significant_markers, "B_all_marker.csv")
B_significant_markers <- B_significant_markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
write.csv(B_significant_markers, "B_top_marker_50.csv")

# ==================== B细胞亚型基因集 ====================
library(Seurat)
library(dplyr)
library(pheatmap)

# ==================== B 细胞完整专属基因集（无交叉污染）====================
# 1. Immature B (未成熟B，骨髓未成熟，专属强)
immature_b_genes <- c(
  "VPREB1", "VPREB3", "IGLL1", "DNTT", "RAG1", "RAG2", "MME", "CD10"
)

transitional_b_genes <- c(
  "CD24",      # 过渡B早期标志
  "CD38",      # 过渡B早期标志
  "BTLA",      # B、T淋巴细胞衰减因子
  "MME",       # 中性内肽酶，在未成熟B表达
  "VPREB3",    # 前B细胞标志，过渡期残留
  "IGLL1",     # 与VPREB3类似，未成熟B特征
  "PTPRC",     # CD45，过渡B高表达
  "TCL1A"      # 辅助区分（TCL1A在过渡B中高，成熟B中也中高，但结合负选）
)

naive_b_genes <- c(
  "MS4A1", "CD19", "CD37", "CD79A", "CD79B", 
  "TCL1A", "LTB", "TXNIP", "CD52", "BTG1"
)

gc_b_genes <- c(
  "BCL6", "AICDA", "MEF2C", "POU2AF1", "TCL1A", 
  "LMO2", "BCL7A", "MYC", "CD38"
)

memory_b_genes <- c(
  "CD27", "CD44", "CD80", "CD86", "CXCR3", 
  "CD200", "ITGB1", "FCRL3", "FCRL4"
)

activated_b_genes <- c(
  "CD69", "EGR1", "NFKBIZ", "DUSP2", "FOSB", 
  "JUNB", "CD40", "CD70", "IL4R"
)

plasmablast_genes <- c(
  "IRF4", "PRDM1", "XBP1", "SDC1", "CD38", 
  "MKI67", "IGHM", "DERL3"
)


plasma_genes <- c(
  "JCHAIN", "MZB1", "XBP1", "IGKC", "IGLC1", 
  "IGHG1", "IGHA1", "SDC1", "FKBP11", "HERPUD1"
)
# Regulatory B（Breg，免疫抑制）

breg_genes <- c(
  "IL10", "TGFB1", "CD5", "CD1D", "PDCD1LG2", 
  "LAG3", "CD200", "CD14"
)
# Proliferating B（增殖B，不分亚型）

proliferating_b_genes <- c(
  "MKI67", "TOP2A", "PCNA", "CCNB1", "CCNA2", "CDK1",
  "BIRC5", "AURKB", "PLK1", "TYMS", "RRM2", "ASPM"
)
mature_b_genes <- c(
  "CD19", "CD22", "BLK", "BANK1", "PAX5", "POU2F2", "CIITA", "FCER2","MS4A1","CD79A", "CD79B"
)

# ==================== 过滤在B对象中存在的基因 ====================
immature_b_genes     <- immature_b_genes[immature_b_genes %in% rownames(B)]
mature_b_genes     <- mature_b_genes[mature_b_genes %in% rownames(B)]
transitional_b_genes <- transitional_b_genes[transitional_b_genes %in% rownames(B)]
naive_b_genes        <- naive_b_genes[naive_b_genes %in% rownames(B)]

gc_b_genes           <- gc_b_genes[gc_b_genes %in% rownames(B)]
memory_b_genes       <- memory_b_genes[memory_b_genes %in% rownames(B)]
activated_b_genes    <- activated_b_genes[activated_b_genes %in% rownames(B)]
plasmablast_genes    <- plasmablast_genes[plasmablast_genes %in% rownames(B)]
plasma_genes         <- plasma_genes[plasma_genes %in% rownames(B)]
breg_genes           <- breg_genes[breg_genes %in% rownames(B)]
proliferating_b_genes <- proliferating_b_genes[proliferating_b_genes %in% rownames(B)]

# 输出数量
cat("Immature B:    ", length(immature_b_genes), "\n")
cat("mature B:    ", length(mature_b_genes), "\n")
cat("Transitional B:", length(transitional_b_genes), "\n")
cat("Naive B:       ", length(naive_b_genes), "\n")
cat("GC B:          ", length(gc_b_genes), "\n")
cat("Memory B:      ", length(memory_b_genes), "\n")
cat("Activated B:   ", length(activated_b_genes), "\n")
cat("Plasmablast:   ", length(plasmablast_genes), "\n")
cat("Plasma:        ", length(plasma_genes), "\n")
cat("Breg:          ", length(breg_genes), "\n")
cat("增殖 B 基因数:", length(proliferating_b_genes), "\n")

# ==================== 模块评分 ====================
B <- AddModuleScore(B, features = list(immature_b_genes),     name = "ImmatureB_Score")
B <- AddModuleScore(B, features = list(mature_b_genes),     name = "MatureB_Score")
B <- AddModuleScore(B, features = list(transitional_b_genes), name = "TransitionalB_Score")
B <- AddModuleScore(B, features = list(naive_b_genes),        name = "NaiveB_Score")
B <- AddModuleScore(B, features = list(gc_b_genes),           name = "GCB_Score")
B <- AddModuleScore(B, features = list(memory_b_genes),       name = "MemoryB_Score")
B <- AddModuleScore(B, features = list(activated_b_genes),    name = "ActivatedB_Score")
B <- AddModuleScore(B, features = list(plasmablast_genes),    name = "Plasmablast_Score")
B <- AddModuleScore(B, features = list(plasma_genes),         name = "Plasma_Score")
B <- AddModuleScore(B, features = list(breg_genes),           name = "Breg_Score")
B <- AddModuleScore(B, features = list(proliferating_b_genes), name = "Proliferating_Score")

# ==================== 按cluster计算平均评分 ====================
score_cols <- c(
  "ImmatureB_Score1", "TransitionalB_Score1", "NaiveB_Score1",
   "GCB_Score1", "MemoryB_Score1",
  "ActivatedB_Score1", "Plasmablast_Score1", "Plasma_Score1", "Breg_Score1", "Proliferating_Score1"
)

score_avg <- B@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(across(all_of(score_cols), mean, na.rm = TRUE))

score_mat <- as.matrix(score_avg[, -1])
rownames(score_mat) <- score_avg$seurat_clusters

# ==================== 绘制热图 ====================
my_color <- colorRampPalette(c("royalblue", "white", "red"))(100)
my_breaks <- seq(-max(abs(score_mat)), max(abs(score_mat)), length.out = 101)

# 按行缩放（推荐用于文章）
pheatmap(score_mat,
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         main = "B cell subset scores by cluster",
         fontsize = 10,
         color = my_color,
         filename = "B_subset_scores_heatmap.pdf",
         width = 10, height = 8)
# 原始值热图
max_raw <- max(abs(score_mat))
my_breaks_raw <- seq(-max_raw, max_raw, length.out = 101)

pheatmap(score_mat,
         scale = "none",
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         main = "B cell subset scores (raw)",
         fontsize = 10,
         color = my_color,
         breaks = my_breaks_raw,
         filename = "B_subset_scores_raw.pdf",
         width = 10, height = 8)


# ==================== 输出表格 ====================
print(score_avg)

# ==================== 自动命名 ====================
score_cols_short <- c(
  "Immature B", "Transitional B", "Naive B", 
  "GC B", "Memory B", "Activated B", "Plasmablast", "Plasma", "Breg","Proliferating B"
)
colnames(score_mat) <- score_cols_short

max_score <- apply(score_mat, 1, function(x) colnames(score_mat)[which.max(x)])
max_value <- apply(score_mat, 1, max)

cat("\n===== B细胞最终命名建议 =====\n")
for (i in seq_along(max_score)) {
  cat("Cluster", rownames(score_mat)[i], ":", max_score[i],
      " (score =", round(max_value[i], 3), ")\n")
}



B <- subset(B,subset=seurat_clusters %in% c(0,1,2,3,4,5,6,7))
B <- RunUMAP(B, dims = 1:30)

identity_mapping <- c(
  "0" = "B_Naive_1",
  "1" = "B_Naive_2",
  "2" = "B_Plasma_1",
  "3" = "B_Plasma_2",
  "4" = "B_Naive_3",
  "5" = "B_Germinal-center",
  "6" = "B_Proliferative",
  "7" = "B_Activated"
)
B$seurat_clusters <- factor(B$seurat_clusters)
sub_cell_type <- identity_mapping[B@meta.data$seurat_clusters]
B@meta.data$sub_cell_type <- sub_cell_type

library(ggplot2)
library(ggsci)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(5)
cluster_colors <- c(
  # 浅色系在前（浅蓝第一个）
  "#A0CBE8",   # 1 浅蓝（你要的放第一个）
  "#FFBE7D",   # 2 浅橙
  "#8CD17D",   # 3 浅绿
  "#86BCB6",   # 4 浅青
  "#FF9D9A",   # 5 浅红
  "#FABFD2",   # 6 浅粉
  "#D4A6C8",   # 7 浅紫
  "#F1CE63",   # 8 淡黄
  "#D7B5A6",   # 9 浅棕
  "#B07AA1",   # 10 紫灰
  "#BAB0AC",   # 11 浅灰
  
  # 深色系在后
  "#4E79A7",   # 12 深蓝
  "#F28E2B",   # 13 橙
  "#59A14F",   # 14 深绿
  "#E15759",   # 15 红
  "#499894",   # 16 青绿
  "#D37295",   # 17 玫红
  "#B6992D",   # 18 土黄
  "#9D7660",   # 19 棕
  "#79706E"    # 20 深灰
)
cell_types <- as.character(unique(B@meta.data$sub_cell_type))
num_legend_items <- length(cell_types)  # 图例的个数
max_label_length <- max(nchar(cell_types))  # 图例名称的最大长度

base_width <- 3000  # 基础宽度
base_height <- 3000  # 基础高度
legend_width_factor <- 100  # 每个图例项增加的宽度
label_length_factor <- 10  # 每个字符增加的宽度
dynamic_width <- base_width + (num_legend_items * legend_width_factor) + (max_label_length * label_length_factor)
pdf("B_annotation.pdf", width = dynamic_width/300, height = base_height/300)
DimPlot(B, reduction = "umap", label = TRUE, pt.size = 1, group.by = "sub_cell_type", label.size = 4) +
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    ggtitle(NULL) +
    scale_color_manual(values = cluster_colors ) +
    coord_fixed(ratio = 1) +
    guides(color = guide_legend(title = NULL, override.aes = list(size = 5))) +
    theme(
        text = element_text(size = 8, face = "bold"),
        axis.text.x = element_text(size = 28, color = "black"),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title.x = element_text(size = 36, face = "bold", color = "black"),
        axis.title.y = element_text(size = 36, face = "bold", color = "black", margin = margin(r = 20)),  # 增加右侧间距
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 16, face = "bold", color = "black"),
        legend.title = element_text(size = 16, face = "bold", color = "black"),
        legend.position = "right",
        legend.box.margin = margin(0, 0, 0, 0),
        legend.key = element_blank(),
        legend.background = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5),
        aspect.ratio = 1,
        plot.margin = margin(10, 50, 10, 10)
    )
dev.off()
```
# B sc
```R
obj <- readRDS("sc_obj_anno.rds")
library(Seurat)
library(ggplot2)
library(dplyr)

B <- subset(obj,subset=cell_type %in% c("B cells"))
B <- NormalizeData(B)
B <- FindVariableFeatures(B, nfeatures = 2000)
hvgs <- VariableFeatures(B)
B <- ScaleData(B, features = hvgs)
B <- RunPCA(B, features = hvgs, npcs = 50)

p <- ElbowPlot(B, ndims = 30)
ggsave("p3.png", plot = p)
pca_var <- B[["pca"]]@stdev^2 / sum(B[["pca"]]@stdev^2)
cumsum_var <- cumsum(pca_var)

# 看达到多少 PC 时累计方差 > 70% / 80% / 90%
which(cumsum_var > 0.7)[1]   # 70%
which(cumsum_var > 0.8)[1]   # 80%
which(cumsum_var > 0.9)[1]   # 90%

library(harmony)
B <- RunHarmony(B, "dataset")
B <- FindNeighbors(B, reduction = "harmony", dims = 1:27)
B <- FindClusters(B, resolution = 0.2)
B <- RunUMAP(B, reduction = "harmony", dims = 1:27)

library(ggplot2)
library(ggsci)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(17)

seurat_clusters <- as.character(unique(B@meta.data$seurat_clusters))  # 转换为字符向量
num_legend_items <- length(seurat_clusters)  # 图例的个数
max_label_length <- max(nchar(seurat_clusters))  # 图例名称的最大长度

base_width <- 3000  # 基础宽度
base_height <- 3000  # 基础高度
legend_width_factor <- 100  # 每个图例项增加的宽度
label_length_factor <- 10  # 每个字符增加的宽度
dynamic_width <- base_width + (num_legend_items * legend_width_factor) + (max_label_length * label_length_factor)
pdf("B_sc_clusters.pdf", width = dynamic_width / 300, height = base_height / 300)  # 转换为英寸

DimPlot(B, reduction = "umap", label = TRUE, pt.size = 1, label.size = 8) +
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    ggtitle(NULL) +
    scale_color_manual(values = npg_extended) +
    coord_fixed(ratio = 1) +
    guides(color = guide_legend(title = NULL, override.aes = list(size = 5))) +
    theme(
        text = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 28, color = "black"),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title.x = element_text(size = 36, face = "bold", color = "black"),
        axis.title.y = element_text(size = 36, face = "bold", color = "black", margin = margin(r = 20)),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 28, face = "bold", color = "black"),
        legend.title = element_text(size = 28, face = "bold", color = "black"),
        legend.position = "right",
        legend.box.margin = margin(0, 0, 0, 0),
        legend.key = element_blank(),
        legend.background = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5),
        aspect.ratio = 1,
        plot.margin = margin(10, 50, 10, 10)
    )

dev.off()

p <- DimPlot(B, reduction = "umap", group.by = "dataset", label = FALSE)
ggsave("B_orig.png", plot = p)

library(dplyr)
B_markers <- FindAllMarkers(B, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "wilcox")
B_significant_markers <- subset(B_markers, p_val_adj < 0.05)
#write.csv(B_significant_markers, "B_all_marker.csv")
B_significant_markers <- B_significant_markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
write.csv(B_significant_markers, "B_sc_top_marker_50.csv")

# ==================== 绘制热图 ====================
my_color <- colorRampPalette(c("royalblue", "white", "red"))(100)
my_breaks <- seq(-max(abs(score_mat)), max(abs(score_mat)), length.out = 101)

# 按行缩放（推荐用于文章）
pheatmap(score_mat,
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         main = "B cell subset scores by cluster",
         fontsize = 10,
         color = my_color,
         filename = "B_sc_subset_scores_heatmap.pdf",
         width = 10, height = 8)
# 原始值热图
max_raw <- max(abs(score_mat))
my_breaks_raw <- seq(-max_raw, max_raw, length.out = 101)

pheatmap(score_mat,
         scale = "none",
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         main = "B cell subset scores (raw)",
         fontsize = 10,
         color = my_color,
         breaks = my_breaks_raw,
         filename = "B_sc_subset_scores_raw.pdf",
         width = 10, height = 8)


# ==================== 输出表格 ====================
print(score_avg)

# ==================== 自动命名 ====================
score_cols_short <- c(
  "Immature B", "Transitional B", "Naive B", 
  "GC B", "Memory B", "Activated B", "Plasmablast", "Plasma", "Breg","Proliferating B"
)
colnames(score_mat) <- score_cols_short

max_score <- apply(score_mat, 1, function(x) colnames(score_mat)[which.max(x)])
max_value <- apply(score_mat, 1, max)

cat("\n===== B细胞最终命名建议 =====\n")
for (i in seq_along(max_score)) {
  cat("Cluster", rownames(score_mat)[i], ":", max_score[i],
      " (score =", round(max_value[i], 3), ")\n")
}

B <- subset(B,subset=seurat_clusters %in% c(0,1,2,3,5,6,9))
B <- RunUMAP(B, reduction = "harmony",dims = 1:27)

identity_mapping <- c(
  "0" = "B_Memory_1",
  "1" = "B_Naive",
  "2" = "B_Memory_2",
  "3" = "B_Activated",
  "5" = "B_Germinal-center",
  "6" = "B_Plasma_1",
  "9" = "B_Plasma_2"
)
B$seurat_clusters <- factor(B$seurat_clusters)
sub_cell_type <- identity_mapping[B@meta.data$seurat_clusters]
B@meta.data$sub_cell_type <- sub_cell_type

library(ggplot2)
library(ggsci)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(5)
cluster_colors <- c(
  # 浅色系在前（浅蓝第一个）
  "#A0CBE8",   # 1 浅蓝（你要的放第一个）
  "#FFBE7D",   # 2 浅橙
  "#8CD17D",   # 3 浅绿
  "#86BCB6",   # 4 浅青
  "#FF9D9A",   # 5 浅红
  "#FABFD2",   # 6 浅粉
  "#D4A6C8",   # 7 浅紫
  "#F1CE63",   # 8 淡黄
  "#D7B5A6",   # 9 浅棕
  "#B07AA1",   # 10 紫灰
  "#BAB0AC",   # 11 浅灰
  
  # 深色系在后
  "#4E79A7",   # 12 深蓝
  "#F28E2B",   # 13 橙
  "#59A14F",   # 14 深绿
  "#E15759",   # 15 红
  "#499894",   # 16 青绿
  "#D37295",   # 17 玫红
  "#B6992D",   # 18 土黄
  "#9D7660",   # 19 棕
  "#79706E"    # 20 深灰
)
cell_types <- as.character(unique(B@meta.data$sub_cell_type))
num_legend_items <- length(cell_types)  # 图例的个数
max_label_length <- max(nchar(cell_types))  # 图例名称的最大长度

base_width <- 3000  # 基础宽度
base_height <- 3000  # 基础高度
legend_width_factor <- 100  # 每个图例项增加的宽度
label_length_factor <- 10  # 每个字符增加的宽度
dynamic_width <- base_width + (num_legend_items * legend_width_factor) + (max_label_length * label_length_factor)
pdf("B_sc_annotation.pdf", width = dynamic_width/300, height = base_height/300)
DimPlot(B, reduction = "umap", label = TRUE, pt.size = 1, group.by = "sub_cell_type", label.size = 4) +
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    ggtitle(NULL) +
    scale_color_manual(values = cluster_colors ) +
    coord_fixed(ratio = 1) +
    guides(color = guide_legend(title = NULL, override.aes = list(size = 5))) +
    theme(
        text = element_text(size = 8, face = "bold"),
        axis.text.x = element_text(size = 28, color = "black"),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title.x = element_text(size = 36, face = "bold", color = "black"),
        axis.title.y = element_text(size = 36, face = "bold", color = "black", margin = margin(r = 20)),  # 增加右侧间距
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 16, face = "bold", color = "black"),
        legend.title = element_text(size = 16, face = "bold", color = "black"),
        legend.position = "right",
        legend.box.margin = margin(0, 0, 0, 0),
        legend.key = element_blank(),
        legend.background = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5),
        aspect.ratio = 1,
        plot.margin = margin(10, 50, 10, 10)
    )
dev.off()
```
# mPT all subtype chat
```R
library(Seurat)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(deldir)
library(doParallel)
library(foreach)

# ===================== 数据加载 =====================
obj <- readRDS("YA2025263-1_fin.rds")
Malignant <- readRDS("malignant_anno.rds")
fib <- readRDS("fib_anno_new.rds")
macro <- readRDS("macro_anno_new.rds")
t_cells <- readRDS("t_anno.rds")
DC <- readRDS("DC_anno.rds")
B <- readRDS("B_anno.rds")
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
                       "Basal cells" = "Malignant cells"
)

# 检查
table(obj$CellType)

obj@meta.data$sub_cell_type <- NA
common_cells <- intersect(rownames(Malignant@meta.data), rownames(obj@meta.data))
obj@meta.data[common_cells, "sub_cell_type"] <- Malignant@meta.data[common_cells, "sub_cell_type"]

common_cells <- intersect(rownames(fib@meta.data), rownames(obj@meta.data))
obj@meta.data[common_cells, "sub_cell_type"] <- fib@meta.data[common_cells, "sub_cell_type"]

common_cells <- intersect(rownames(macro@meta.data), rownames(obj@meta.data))
obj@meta.data[common_cells, "sub_cell_type"] <- macro@meta.data[common_cells, "sub_cell_type"]

common_cells <- intersect(rownames(t_cells@meta.data), rownames(obj@meta.data))
obj@meta.data[common_cells, "sub_cell_type"] <- t_cells@meta.data[common_cells, "sub_cell_type"]

common_cells <- intersect(rownames(DC@meta.data), rownames(obj@meta.data))
obj@meta.data[common_cells, "sub_cell_type"] <- DC@meta.data[common_cells, "sub_cell_type"]

common_cells <- intersect(rownames(B@meta.data), rownames(obj@meta.data))
obj@meta.data[common_cells, "sub_cell_type"] <- B@meta.data[common_cells, "sub_cell_type"]

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

remove_types <- c("B cells", "Dendritic cells", "Germinal center B cells", 
                  "immature B cells", "Malignant cells", "T cells")

# 方法：保留不在移除列表中的细胞
cells_keep <- which(!(obj@meta.data$detailed %in% remove_types))
obj_filtered <- obj[, cells_keep]

# 或者用逻辑向量
keep <- !(obj@meta.data$detailed %in% remove_types)
obj_filtered <- obj[, keep]

# 检查结果
ncol(obj_filtered)
table(obj_filtered@meta.data$detailed)

#chat_mPT

#画连通图
library(ggplot2)
library(deldir)
library(igraph)

# ==================== 1. 提取 C5 样本 FOV 50 的数据 ====================
sample_cells <- subset(mPT, subset = sample == "C5" & fov == 50)

# 获取坐标和标签
coords <- sample_cells@meta.data[, c("CenterX_global_px", "CenterY_global_px")]
coords <- coords[!is.na(coords[,1]) & !is.na(coords[,2]), , drop = FALSE]
labels <- as.character(sample_cells@meta.data$detailed)
names(labels) <- rownames(coords)

# 移除NA标签
valid_idx <- !is.na(labels)
coords <- coords[valid_idx, , drop = FALSE]
labels <- labels[valid_idx]

cat("FOV 50 有效细胞数量：", length(labels), "\n")

# ==================== 2. 计算 Delaunay 三角剖分和连通图 ====================
x <- coords[,1]
y <- coords[,2]

deld <- deldir(x, y, rw = c(range(x), range(y)))
segs <- deld$delsgs

# 构建边
edges <- cbind(segs$ind1, segs$ind2)
edges <- edges[edges[,1] != edges[,2], , drop = FALSE]
edges <- t(apply(edges, 1, function(x) sort(x)))
edges <- unique(edges)
edges_df <- data.frame(from = edges[,1], to = edges[,2])

# 计算连通分量
g <- graph_from_data_frame(edges_df, directed = FALSE)
comp <- components(g)

# ==================== 3. 准备绘图数据 ====================
spatial_df <- data.frame(
  x = x,
  y = y,
  cell_type = labels,
  component = as.factor(comp$membership)
)

# 获取边坐标
edge_coords <- do.call(rbind, lapply(1:nrow(edges_df), function(i) {
  from_idx <- edges_df$from[i]
  to_idx <- edges_df$to[i]
  data.frame(
    x = spatial_df$x[from_idx],
    y = spatial_df$y[from_idx],
    xend = spatial_df$x[to_idx],
    yend = spatial_df$y[to_idx],
    component = spatial_df$component[from_idx]
  )
}))

# ==================== 4. 绘图 ====================
p <- ggplot() +
  geom_point(data = spatial_df, aes(x = x, y = y, color = cell_type), 
             size = 1, alpha = 0.7) +
  geom_segment(data = edge_coords, 
               aes(x = x, y = y, xend = xend, yend = yend, group = component),
               color = "gray40", alpha = 0.3, linewidth = 0.2) +
  theme_minimal() +
  coord_fixed() +
  labs(title = "C5 Sample - FOV 50: Spatial Cell Connectivity Network",
       x = "X Coordinate (px)", y = "Y Coordinate (px)") +
  theme(legend.position = "right")

# 保存
ggsave("C5_FOV50_spatial_network.pdf", p, width = 14, height = 10)

#剩下的，要保存各自的文件 
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
                    main = "Average Log2(Fold Change) across mPT samples",
                    fontsize = 10,
                    display_numbers = signif_symbols,
                    number_color = "black",
                    fontsize_number = 8,
                    filename = "mPT_all_samples_average_contact_log2fc_detailed_subtype.pdf",
                    width = 14, height = 12)
  
  
  # 保存整合结果
  integrated_results <- list(
    all_cell_types = all_cell_types,
    avg_log2fc = avg_log2fc,
    combined_p = combined_p,
    n_effective = n_effective,
    n_samples = length(all_results),
    individual_results = all_results
  )
  
  saveRDS(integrated_results, file = "mPT_integrated_results_detailed_subtype.rds")
}

# ===================== 全部细胞类型之间两两互作强度排序（仅显著，排除自身）=====================

# ===================== 加载数据 =====================
integrated_results <- readRDS("mPT_integrated_results_detailed_subtype.rds")
avg_log2fc <- integrated_results$avg_log2fc
combined_p <- integrated_results$combined_p
all_cell_types <- integrated_results$all_cell_types

# 构建映射 (detailed 亚群 → CellType 大类)
library(dplyr)
mapping <- mPT@meta.data %>%
  select(detailed, CellType) %>%
  distinct()

# ===================== 主循环 =====================
interactions <- data.frame()

for (ct1 in all_cell_types) {
  for (ct2 in all_cell_types) {
    
    # 1. 排除自身
    if (ct1 == ct2) next
    
    # 2. 只保留显著
    pval <- combined_p[ct1, ct2]
    log2fc <- avg_log2fc[ct1, ct2]
    if (is.na(pval) || pval >= 0.05) next
    
    # 3. 获取 CellType 大类
    ct1_big <- mapping$CellType[match(ct1, mapping$detailed)]
    ct2_big <- mapping$CellType[match(ct2, mapping$detailed)]
    
    # ===================== 【核心：严格排除】 =====================
    # 规则 1：只要 CellType 相同 → 一律排除
    if (ct1_big == ct2_big) {
      next
    }
    
    # 规则 2：额外排除 Tumor / T / B / DC / Macro / Fib 内部
    g1 <- case_when(
      startsWith(ct1, "Tumor_") ~ "Tumor",
      startsWith(ct1, "T_") ~ "T",
      startsWith(ct1, "B_") ~ "B",
      startsWith(ct1, "DC_") ~ "DC",
      startsWith(ct1, "Macro_") ~ "Macro",
      startsWith(ct1, "Fib_") ~ "Fib",
      TRUE ~ "Other"
    )
    
    g2 <- case_when(
      startsWith(ct2, "Tumor_") ~ "Tumor",
      startsWith(ct2, "T_") ~ "T",
      startsWith(ct2, "B_") ~ "B",
      startsWith(ct2, "DC_") ~ "DC",
      startsWith(ct2, "Macro_") ~ "Macro",
      startsWith(ct2, "Fib_") ~ "Fib",
      TRUE ~ "Other"
    )
    
    # 同类主群内部 → 排除
    if (g1 == g2 & g1 != "Other") {
      next
    }
    
    # ===================== 能走到这里 = 全部保留 =====================
    interactions <- rbind(interactions, data.frame(
      detailed1 = ct1,
      detailed2 = ct2,
      CellType1 = ct1_big,
      CellType2 = ct2_big,
      group1 = g1,
      group2 = g2,
      log2FC = log2fc,
      p_value = pval,
      stringsAsFactors = FALSE
    ))
  }
}

# ===================== 排序：从正到负，由大到小 =====================
if (nrow(interactions) > 0) {
  interactions_sorted <- interactions[order(-interactions$log2FC), ]
  
  interactions_sorted$signif <- ""
  interactions_sorted$signif[interactions_sorted$p_value < 0.05]  <- "*"
  interactions_sorted$signif[interactions_sorted$p_value < 0.01]  <- "**"
  interactions_sorted$signif[interactions_sorted$p_value < 0.001] <- "***"
  
  cat("\n===== 最终显著互作（正→负 大→小）=====\n")
  print(head(interactions_sorted[, c("detailed1","detailed2","CellType1","CellType2","log2FC","signif")], 30))
  
  write.csv(interactions_sorted, "mPT_subtype_significant_interactions.csv", row.names = FALSE)
  cat("\n文件已保存：final_significant_interactions.csv\n")
  
} else {
  cat("没有符合条件的显著互作\n")
}

# ===================== 去重 + 绘图 =====================
library(ggplot2)
library(dplyr)

# 1. 去重：把 A↔B 和 B↔A 视为同一对
unique_inter <- interactions_sorted %>%
  rowwise() %>%
  mutate(pair = paste(sort(c(detailed1, detailed2)), collapse = " ↔ ")) %>%
  ungroup() %>%
  distinct(pair, .keep_all = TRUE) %>%
  arrange(-log2FC)

# 2. 取去重后的 top30
top30 <- head(unique_inter, 30)

# 3. 保证顺序
top30 <- top30 %>% arrange(log2FC)
top30$pair <- factor(top30$pair, levels = top30$pair)

# 4. 绘图（正常出条形）
p <- ggplot(top30, aes(x = pair, y = log2FC, fill = log2FC > 0)) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#E74C3C", "FALSE" = "#6699ff")) +
  geom_text(aes(label = signif), hjust = -0.2, size = 3.5) +
  labs(title = "Top30 Unique Cell-Cell Interactions",
       x = "Cell Pair", y = "log2FC") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 15),
        plot.title = element_text(hjust = 0.5))

ggsave("mPT_subtype_Top30_interactions_no_duplicate.pdf", p, width = 16, height = 10)
cat("✅ 去重完成，条形图已正常生成\n")

# ===================== 最终微调版：字更大 + 图更小 + TOP10 + 拼图 =====================
library(ggplot2)
library(dplyr)
library(patchwork)

# 1. 总结果去重
unique_total <- interactions_sorted %>%
  rowwise() %>%
  mutate(pair_key = paste(sort(c(detailed1, detailed2)), collapse = "_")) %>%
  ungroup() %>%
  distinct(pair_key, .keep_all = TRUE) %>%
  arrange(desc(log2FC))

# 2. 提取肿瘤亚型
tumor_types <- grep("^Tumor_", all_cell_types, value = TRUE)

if (length(tumor_types) == 0) {
  cat("未找到 Tumor 亚型\n")
} else {
  plot_list <- list()
  
  for (tumor in tumor_types) {
    tumor_dat <- unique_total %>%
      filter(detailed1 == !!tumor | detailed2 == !!tumor)
    if (nrow(tumor_dat) == 0) next
    
    # 取 TOP10
    top10 <- head(tumor_dat, 10)
    top10 <- top10 %>% arrange(log2FC)
    top10$label <- paste(top10$detailed1, "↔", top10$detailed2)
    top10$label <- factor(top10$label, levels = unique(top10$label))
    
    # 画图：字更大 + 更清晰
    p <- ggplot(top10, aes(x = label, y = log2FC, fill = log2FC > 0)) +
      geom_col(width = 0.75) +
      coord_flip() +
      scale_fill_manual(values = c("#4A90E2", "#E74C3C")) +
      geom_text(aes(label = signif), hjust = -0.2, size = 4.5, fontface="bold") +
      labs(title = tumor, y = "log2FC", x = "") +
      theme_bw(base_size = 13) +  # 全局字更大
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face="bold"),
        axis.text.y = element_text(size = 12),  # Y轴字明显更大
        axis.text.x = element_text(size = 11),
        axis.title = element_text(size = 12),
        legend.position = "none"
      )
    
    # 保存单图：尺寸更小（宽度8 高度6）
    ggsave(paste0("mPT_subtype_", tumor, "_Top10.pdf"), p, width=8, height=6, dpi=300)
    plot_list[[tumor]] <- p
  }
  
  # 拼图：整体更小
  if (length(plot_list) > 0) {
    combined <- wrap_plots(plot_list, ncol=2) + plot_layout(guides="collect")
    ggsave("mPT_subtype_all_tumors_top10_combined.pdf", combined, width=20, height=5*ceiling(length(plot_list)/2), dpi=300)
    cat("\n✅ 画图完成：字更大 + 图更小！\n")
  }
}
```
# nmPT all subtype chat
```R
library(Seurat)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(deldir)
library(doParallel)
library(foreach)

nmPT <- readRDS("nmPT_detailed.rds")
#chat_nmPT

#后续
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
                    filename = "nmPT_all_samples_average_contact_log2fc_detailed_subtype.pdf",
                    width = 14, height = 12)
  
  # 保存整合结果
  integrated_results <- list(
    all_cell_types = all_cell_types,
    avg_log2fc = avg_log2fc,
    combined_p = combined_p,
    n_effective = n_effective,
    n_samples = length(all_results),
    individual_results = all_results
  )
  
  saveRDS(integrated_results, file = "nmPT_integrated_results_detailed_subtype.rds")
  cat("保存整合结果：nmPT_integrated_results.rds\n")
}

# ===================== 全部细胞类型之间两两互作强度排序（仅显著，排除自身）=====================

# ===================== 加载数据 =====================
integrated_results <- readRDS("nmPT_integrated_results_detailed_subtype.rds")
avg_log2fc <- integrated_results$avg_log2fc
combined_p <- integrated_results$combined_p
all_cell_types <- integrated_results$all_cell_types

# 构建映射 (detailed 亚群 → CellType 大类)
library(dplyr)
mapping <- nmPT@meta.data %>%
  select(detailed, CellType) %>%
  distinct()

# ===================== 主循环 =====================
interactions <- data.frame()

for (ct1 in all_cell_types) {
  for (ct2 in all_cell_types) {
    
    # 1. 排除自身
    if (ct1 == ct2) next
    
    # 2. 只保留显著
    pval <- combined_p[ct1, ct2]
    log2fc <- avg_log2fc[ct1, ct2]
    if (is.na(pval) || pval >= 0.05) next
    
    # 3. 获取 CellType 大类
    ct1_big <- mapping$CellType[match(ct1, mapping$detailed)]
    ct2_big <- mapping$CellType[match(ct2, mapping$detailed)]
    
    # ===================== 【核心：严格排除】 =====================
    # 规则 1：只要 CellType 相同 → 一律排除
    if (ct1_big == ct2_big) {
      next
    }
    
    # 规则 2：额外排除 Tumor / T / B / DC / Macro / Fib 内部
    g1 <- case_when(
      startsWith(ct1, "Tumor_") ~ "Tumor",
      startsWith(ct1, "T_") ~ "T",
      startsWith(ct1, "B_") ~ "B",
      startsWith(ct1, "DC_") ~ "DC",
      startsWith(ct1, "Macro_") ~ "Macro",
      startsWith(ct1, "Fib_") ~ "Fib",
      TRUE ~ "Other"
    )
    
    g2 <- case_when(
      startsWith(ct2, "Tumor_") ~ "Tumor",
      startsWith(ct2, "T_") ~ "T",
      startsWith(ct2, "B_") ~ "B",
      startsWith(ct2, "DC_") ~ "DC",
      startsWith(ct2, "Macro_") ~ "Macro",
      startsWith(ct2, "Fib_") ~ "Fib",
      TRUE ~ "Other"
    )
    
    # 同类主群内部 → 排除
    if (g1 == g2 & g1 != "Other") {
      next
    }
    
    # ===================== 能走到这里 = 全部保留 =====================
    interactions <- rbind(interactions, data.frame(
      detailed1 = ct1,
      detailed2 = ct2,
      CellType1 = ct1_big,
      CellType2 = ct2_big,
      group1 = g1,
      group2 = g2,
      log2FC = log2fc,
      p_value = pval,
      stringsAsFactors = FALSE
    ))
  }
}

# ===================== 排序：从正到负，由大到小 =====================
if (nrow(interactions) > 0) {
  interactions_sorted <- interactions[order(-interactions$log2FC), ]
  
  interactions_sorted$signif <- ""
  interactions_sorted$signif[interactions_sorted$p_value < 0.05]  <- "*"
  interactions_sorted$signif[interactions_sorted$p_value < 0.01]  <- "**"
  interactions_sorted$signif[interactions_sorted$p_value < 0.001] <- "***"
  
  cat("\n===== 最终显著互作（正→负 大→小）=====\n")
  print(head(interactions_sorted[, c("detailed1","detailed2","CellType1","CellType2","log2FC","signif")], 30))
  
  write.csv(interactions_sorted, "nmPT_subtype_significant_interactions.csv", row.names = FALSE)
  cat("\n文件已保存：final_significant_interactions.csv\n")
  
} else {
  cat("没有符合条件的显著互作\n")
}

# ===================== 去重 + 绘图 =====================
library(ggplot2)
library(dplyr)

# 1. 去重：把 A↔B 和 B↔A 视为同一对
unique_inter <- interactions_sorted %>%
  rowwise() %>%
  mutate(pair = paste(sort(c(detailed1, detailed2)), collapse = " ↔ ")) %>%
  ungroup() %>%
  distinct(pair, .keep_all = TRUE) %>%
  arrange(-log2FC)

# 2. 取去重后的 top30
top30 <- head(unique_inter, 30)

# 3. 保证顺序
top30 <- top30 %>% arrange(log2FC)
top30$pair <- factor(top30$pair, levels = top30$pair)

# 4. 绘图（正常出条形）
p <- ggplot(top30, aes(x = pair, y = log2FC, fill = log2FC > 0)) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#E74C3C", "FALSE" = "#6699ff")) +
  geom_text(aes(label = signif), hjust = -0.2, size = 3.5) +
  labs(title = "Top30 Unique Cell-Cell Interactions",
       x = "Cell Pair", y = "log2FC") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 15),
        plot.title = element_text(hjust = 0.5))

ggsave("nmPT_subtype_Top30_interactions_no_duplicate.pdf", p, width = 16, height = 10)
cat("✅ 去重完成，条形图已正常生成\n")

# ===================== 最终微调版：字更大 + 图更小 + TOP10 + 拼图 =====================
library(ggplot2)
library(dplyr)
library(patchwork)

# 1. 总结果去重
unique_total <- interactions_sorted %>%
  rowwise() %>%
  mutate(pair_key = paste(sort(c(detailed1, detailed2)), collapse = "_")) %>%
  ungroup() %>%
  distinct(pair_key, .keep_all = TRUE) %>%
  arrange(desc(log2FC))

# 2. 提取肿瘤亚型
tumor_types <- grep("^Tumor_", all_cell_types, value = TRUE)

if (length(tumor_types) == 0) {
  cat("未找到 Tumor 亚型\n")
} else {
  plot_list <- list()
  
  for (tumor in tumor_types) {
    tumor_dat <- unique_total %>%
      filter(detailed1 == !!tumor | detailed2 == !!tumor)
    if (nrow(tumor_dat) == 0) next
    
    # 取 TOP10
    top10 <- head(tumor_dat, 10)
    top10 <- top10 %>% arrange(log2FC)
    top10$label <- paste(top10$detailed1, "↔", top10$detailed2)
    top10$label <- factor(top10$label, levels = unique(top10$label))
    
    # 画图：字更大 + 更清晰
    p <- ggplot(top10, aes(x = label, y = log2FC, fill = log2FC > 0)) +
      geom_col(width = 0.75) +
      coord_flip() +
      scale_fill_manual(values = c("#4A90E2", "#E74C3C")) +
      geom_text(aes(label = signif), hjust = -0.2, size = 4.5, fontface="bold") +
      labs(title = tumor, y = "log2FC", x = "") +
      theme_bw(base_size = 13) +  # 全局字更大
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face="bold"),
        axis.text.y = element_text(size = 12),  # Y轴字明显更大
        axis.text.x = element_text(size = 11),
        axis.title = element_text(size = 12),
        legend.position = "none"
      )
    
    # 保存单图：尺寸更小（宽度8 高度6）
    ggsave(paste0("nmPT_subtype_", tumor, "_Top10.pdf"), p, width=8, height=6, dpi=300)
    plot_list[[tumor]] <- p
  }
  
  # 拼图：整体更小
  if (length(plot_list) > 0) {
    combined <- wrap_plots(plot_list, ncol=2) + plot_layout(guides="collect")
    ggsave("nmPT_subtype_all_tumors_top10_combined.pdf", combined, width=20, height=5*ceiling(length(plot_list)/2), dpi=300)
    cat("\n✅ 画图完成：字更大 + 图更小！\n")
  }
}
```
# metLN all subtype chat
```R
library(Seurat)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(deldir)
library(doParallel)
library(foreach)

metLN <- readRDS("metLN_detailed.rds")
samples <- unique(metLN$sample)

all_results <- list()

for (sample_id in samples) {
  cat("\n\n########## 处理样本：", sample_id, " ##########\n")
  
  # 提取该样本的细胞
  sample_cells <- subset(metLN, subset = sample == sample_id)
  
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
                    #filename = paste0("mPT_", sample_id, "_contact_log2fc.pdf"),
                    width = 12, height = 10)
      
      cat("保存热图：mPT_", sample_id, "_contact_log2fc.pdf\n", sep = "")
    }
    
    # 保存结果文件
    #saveRDS(result, file = paste0("mPT_", sample_id, "_results.rds"))
    cat("保存结果：mPT_", sample_id, "_results.rds\n", sep = "")
  }
}
saveRDS(all_results, file = "metLN_all_results_subtype_interaction.rds")

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
                    main = "Average Log2(Fold Change) across metLN samples",
                    fontsize = 10,
                    display_numbers = signif_symbols,
                    number_color = "black",
                    fontsize_number = 8,
                    filename = "metLN_all_samples_average_contact_log2fc_detailed_subtype.pdf",
                    width = 14, height = 12)

  
  # 保存整合结果
  integrated_results <- list(
    all_cell_types = all_cell_types,
    avg_log2fc = avg_log2fc,
    combined_p = combined_p,
    n_effective = n_effective,
    n_samples = length(all_results),
    individual_results = all_results
  )
  
  saveRDS(integrated_results, file = "metLN_integrated_results_detailed_subtype.rds")
}

# ===================== 全部细胞类型之间两两互作强度排序（仅显著，排除自身）=====================

# ===================== 加载数据 =====================
integrated_results <- readRDS("metLN_integrated_results_detailed_subtype.rds")
avg_log2fc <- integrated_results$avg_log2fc
combined_p <- integrated_results$combined_p
all_cell_types <- integrated_results$all_cell_types

# 构建映射 (detailed 亚群 → CellType 大类)
library(dplyr)
mapping <- metLN@meta.data %>%
  select(detailed, CellType) %>%
  distinct()

# ===================== 主循环 =====================
interactions <- data.frame()

for (ct1 in all_cell_types) {
  for (ct2 in all_cell_types) {
    
    # 1. 排除自身
    if (ct1 == ct2) next
    
    # 2. 只保留显著
    pval <- combined_p[ct1, ct2]
    log2fc <- avg_log2fc[ct1, ct2]
    if (is.na(pval) || pval >= 0.05) next
    
    # 3. 获取 CellType 大类
    ct1_big <- mapping$CellType[match(ct1, mapping$detailed)]
    ct2_big <- mapping$CellType[match(ct2, mapping$detailed)]
    
    # ===================== 【核心：严格排除】 =====================
    # 规则 1：只要 CellType 相同 → 一律排除
    if (ct1_big == ct2_big) {
      next
    }
    
    # 规则 2：额外排除 Tumor / T / B / DC / Macro / Fib 内部
    g1 <- case_when(
      startsWith(ct1, "Tumor_") ~ "Tumor",
      startsWith(ct1, "T_") ~ "T",
      startsWith(ct1, "B_") ~ "B",
      startsWith(ct1, "DC_") ~ "DC",
      startsWith(ct1, "Macro_") ~ "Macro",
      startsWith(ct1, "Fib_") ~ "Fib",
      TRUE ~ "Other"
    )
    
    g2 <- case_when(
      startsWith(ct2, "Tumor_") ~ "Tumor",
      startsWith(ct2, "T_") ~ "T",
      startsWith(ct2, "B_") ~ "B",
      startsWith(ct2, "DC_") ~ "DC",
      startsWith(ct2, "Macro_") ~ "Macro",
      startsWith(ct2, "Fib_") ~ "Fib",
      TRUE ~ "Other"
    )
    
    # 同类主群内部 → 排除
    if (g1 == g2 & g1 != "Other") {
      next
    }
    
    # ===================== 能走到这里 = 全部保留 =====================
    interactions <- rbind(interactions, data.frame(
      detailed1 = ct1,
      detailed2 = ct2,
      CellType1 = ct1_big,
      CellType2 = ct2_big,
      group1 = g1,
      group2 = g2,
      log2FC = log2fc,
      p_value = pval,
      stringsAsFactors = FALSE
    ))
  }
}

# ===================== 排序：从正到负，由大到小 =====================
if (nrow(interactions) > 0) {
  interactions_sorted <- interactions[order(-interactions$log2FC), ]
  
  interactions_sorted$signif <- ""
  interactions_sorted$signif[interactions_sorted$p_value < 0.05]  <- "*"
  interactions_sorted$signif[interactions_sorted$p_value < 0.01]  <- "**"
  interactions_sorted$signif[interactions_sorted$p_value < 0.001] <- "***"
  
  cat("\n===== 最终显著互作（正→负 大→小）=====\n")
  print(head(interactions_sorted[, c("detailed1","detailed2","CellType1","CellType2","log2FC","signif")], 30))
  
  write.csv(interactions_sorted, "metLN_subtype_significant_interactions.csv", row.names = FALSE)
  cat("\n文件已保存：final_significant_interactions.csv\n")
  
} else {
  cat("没有符合条件的显著互作\n")
}

# ===================== 去重 + 绘图 =====================
library(ggplot2)
library(dplyr)

# 1. 去重：把 A↔B 和 B↔A 视为同一对
unique_inter <- interactions_sorted %>%
  rowwise() %>%
  mutate(pair = paste(sort(c(detailed1, detailed2)), collapse = " ↔ ")) %>%
  ungroup() %>%
  distinct(pair, .keep_all = TRUE) %>%
  arrange(-log2FC)

# 2. 取去重后的 top30
top30 <- head(unique_inter, 30)

# 3. 保证顺序
top30 <- top30 %>% arrange(log2FC)
top30$pair <- factor(top30$pair, levels = top30$pair)

# 4. 绘图（正常出条形）
p <- ggplot(top30, aes(x = pair, y = log2FC, fill = log2FC > 0)) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#E74C3C", "FALSE" = "#6699ff")) +
  geom_text(aes(label = signif), hjust = -0.2, size = 3.5) +
  labs(title = "Top30 Unique Cell-Cell Interactions",
       x = "Cell Pair", y = "log2FC") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 15),
        plot.title = element_text(hjust = 0.5))

ggsave("metLN_subtype_Top30_interactions_no_duplicate.pdf", p, width = 16, height = 10)
cat("✅ 去重完成，条形图已正常生成\n")

# ===================== 最终微调版：字更大 + 图更小 + TOP10 + 拼图 =====================
library(ggplot2)
library(dplyr)
library(patchwork)

# 1. 总结果去重
unique_total <- interactions_sorted %>%
  rowwise() %>%
  mutate(pair_key = paste(sort(c(detailed1, detailed2)), collapse = "_")) %>%
  ungroup() %>%
  distinct(pair_key, .keep_all = TRUE) %>%
  arrange(desc(log2FC))

# 2. 提取肿瘤亚型
tumor_types <- grep("^Tumor_", all_cell_types, value = TRUE)

if (length(tumor_types) == 0) {
  cat("未找到 Tumor 亚型\n")
} else {
  plot_list <- list()
  
  for (tumor in tumor_types) {
    tumor_dat <- unique_total %>%
      filter(detailed1 == !!tumor | detailed2 == !!tumor)
    if (nrow(tumor_dat) == 0) next
    
    # 取 TOP10
    top10 <- head(tumor_dat, 10)
    top10 <- top10 %>% arrange(log2FC)
    top10$label <- paste(top10$detailed1, "↔", top10$detailed2)
    top10$label <- factor(top10$label, levels = unique(top10$label))
    
    # 画图：字更大 + 更清晰
    p <- ggplot(top10, aes(x = label, y = log2FC, fill = log2FC > 0)) +
      geom_col(width = 0.75) +
      coord_flip() +
      scale_fill_manual(values = c("#4A90E2", "#E74C3C")) +
      geom_text(aes(label = signif), hjust = -0.2, size = 4.5, fontface="bold") +
      labs(title = tumor, y = "log2FC", x = "") +
      theme_bw(base_size = 13) +  # 全局字更大
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face="bold"),
        axis.text.y = element_text(size = 12),  # Y轴字明显更大
        axis.text.x = element_text(size = 11),
        axis.title = element_text(size = 12),
        legend.position = "none"
      )
    
    # 保存单图：尺寸更小（宽度8 高度6）
    ggsave(paste0("metLN_subtype_", tumor, "_Top10.pdf"), p, width=8, height=6, dpi=300)
    plot_list[[tumor]] <- p
  }
  
  # 拼图：整体更小
  if (length(plot_list) > 0) {
    combined <- wrap_plots(plot_list, ncol=2) + plot_layout(guides="collect")
    ggsave("metLN_subtype_all_tumors_top10_combined.pdf", combined, width=20, height=5*ceiling(length(plot_list)/2), dpi=300)
    cat("\n✅ 画图完成：字更大 + 图更小！\n")
  }
}
```
# mPT subtype CN
```R
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

#CN组成
# 查看每个CN的细胞组成
cn_composition <- table(mPT@meta.data$cellcharter_cluster, 
                        mPT@meta.data$detailed)

# 计算每个CN的比例
cn_prop <- prop.table(cn_composition, margin = 1) * 100

# 打印每个CN的主要细胞类型（前5种）
for (cn in rownames(cn_prop)) {
  cat("\n=== CN", cn, "===\n")
  top10 <- sort(cn_prop[cn, ], decreasing = TRUE)[1:10]
  for (i in 1:length(top10)) {
    cat(names(top10)[i], ": ", round(top10[i], 1), "%\n", sep = "")
  }
}

#CN注释

cn_annotation_detailed <- c(
  "1"  = "Alveolar_Homeostatic_Zone",
  "2"  = "Immune_Zone",
  "3"  = "Stromal_Remodeling_Zone",
  "4"  = "TLS_Zone",
  "5"  = "Inflammatory_Vascular_Zone",
  "6"  = "Macrophage_Dominant_Zone_A",
  "7"  = "Vascular_Stromal_Zone",
  "8"  = "Pericyte_Dominant_Zone",
  "9"  = "Plasma_Enriched_Zone_A",
  "10" = "Plasma_Enriched_Zone_B",
  "11" = "Metabolic_Tumor_Core_A",
  "12" = "Metabolic_Tumor_Mast_Interface",
  "13" = "Metabolic_Tumor_Core_B",
  "14" = "Metabolic_Tumor_Stromal_Interface",
  "15" = "Plasma_Enriched_Zone_C"
)
# 应用到 Seurat 对象
mPT@meta.data$cn_detailed <- cn_annotation_detailed[as.character(mPT@meta.data$cellcharter_cluster)]

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# 定义颜色（你的20种颜色）
cluster_colors <- c(
  "#A0CBE8", "#FFBE7D", "#8CD17D", "#86BCB6", "#FF9D9A",
  "#FABFD2", "#D4A6C8", "#F1CE63", "#D7B5A6", "#B07AA1",
  "#BAB0AC", "#4E79A7", "#F28E2B", "#59A14F", "#E15759",
  "#499894", "#D37295", "#B6992D", "#9D7660", "#79706E"
)

# =====================
# 1. 准备数据
# =====================

# 假设你的 mPT 对象中有：
# - cn_detailed: CN编号 (1-15)
# - sample: 样本名称

# 计算总体比例
total_prop <- mPT@meta.data %>%
  group_by(cn_detailed) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(proportion = count / sum(count) * 100,
         type = "Overall")

# 计算各样本比例
sample_prop <- mPT@meta.data %>%
  group_by(sample, cn_detailed) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(sample) %>%
  mutate(proportion = count / sum(count) * 100,
         type = sample) %>%
  ungroup()

# 合并数据
combined_data <- bind_rows(total_prop, sample_prop)

# 确保 cn_detailed 是因子（按数字顺序）
combined_data$cn_detailed <- factor(combined_data$cn_detailed, 
                                     levels = sort(unique(combined_data$cn_detailed)))

# =====================
# 方法1：分面图（推荐）
# =====================

# 总体 + 各样本分面显示
p_facet <- ggplot(combined_data, aes(x = type, y = proportion, fill = cn_detailed)) +
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  scale_fill_manual(values = cluster_colors, name = "Cell Neighborhood") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "CN Proportion by Sample",
       x = "",
       y = "Proportion") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "right",
    legend.text = element_text(size = 8),
    strip.background = element_rect(fill = "lightgray"),
    strip.text = element_text(size = 10, face = "bold")
  )

ggsave("mPT_cn_proportion_by_sample_facet.pdf", p_facet, width = 6, height = 6, dpi = 300)

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

# 转换为宽格式
heatmap_data <- cluster_composition_df %>%
  select(Cluster, CellType, Proportion) %>%
  pivot_wider(names_from = CellType, values_from = Proportion, values_fill = 0)

# 按Cluster数字排序
cluster_numeric <- as.numeric(as.character(heatmap_data$Cluster))
heatmap_data <- heatmap_data[order(cluster_numeric), ]

# 创建行标签（使用你的注释名称）
row_labels_annotated <- cn_annotation_detailed[as.character(cluster_numeric[order(cluster_numeric)])]

# 设置行名
rownames(heatmap_data) <- row_labels_annotated
heatmap_matrix <- as.matrix(heatmap_data[, -1])

# 绘制热图（字体增大的版本）
p_heatmap <- pheatmap(heatmap_matrix,
                    cluster_rows = TRUE,
                    cluster_cols = TRUE,
                    display_numbers = FALSE,
                    main = "Cell Type Composition by Spatial Domain",
                    fontsize = 14,           # 全局字体增大
                    fontsize_row = 10,       # 行标签字体（注释名称较长，用10合适）
                    fontsize_col = 9,        # 列标签字体
                    labels_row = row_labels_annotated,
                    color = colorRampPalette(c("#f0f0f0", "#3182bd", "#08519c"))(100),
                    cellwidth = 11,          # 单元格宽度
                    cellheight = 11)         # 单元格高度

# 保存
pdf("mPT_CN_cell_type_composition_heatmap_annotated.pdf", width = 16, height = 12)
print(p_heatmap)
dev.off()

library(ggplot2)
library(dplyr)
library(scales)
library(RColorBrewer)

# 计算每个 CN 中每种细胞类型的数量
stack_data_count <- mPT@meta.data %>%
  group_by(cn_detailed, detailed) %>%
  summarise(count = n(), .groups = "drop")

# 为每个 CN 找出 top 10 细胞类型
top10_per_cn <- stack_data_count %>%
  group_by(cn_detailed) %>%
  arrange(desc(count)) %>%
  slice_head(n = 10) %>%
  ungroup()

# 将非 top 10 的细胞类型合并为 "Other"
stack_data_top10 <- stack_data_count %>%
  left_join(top10_per_cn %>% select(cn_detailed, detailed) %>% mutate(is_top10 = TRUE),
            by = c("cn_detailed", "detailed")) %>%
  mutate(detailed_grouped = ifelse(is_top10, detailed, "Other")) %>%
  group_by(cn_detailed, detailed_grouped) %>%
  summarise(count = sum(count), .groups = "drop")

# 按首字母顺序排序 CN
cn_order <- sort(unique(stack_data_top10$cn_detailed))
stack_data_top10$cn_detailed <- factor(stack_data_top10$cn_detailed, levels = cn_order)

# 获取所有需要着色的类别
all_cell_types <- unique(stack_data_top10$detailed_grouped)
n_colors_needed <- length(all_cell_types)

# 方案1：使用彩虹色（qualitative调色板）
# 提取足够多的颜色
if(n_colors_needed <= 12) {
  cell_colors <- brewer.pal(n_colors_needed, "Set3")
} else {
  # 使用 colorRampPalette 生成更多颜色
  get_palette <- colorRampPalette(brewer.pal(12, "Set3"))
  cell_colors <- get_palette(n_colors_needed)
}
names(cell_colors) <- all_cell_types

# 确保 "Other" 是灰色
if("Other" %in% names(cell_colors)) {
  cell_colors["Other"] <- "#D3D3D3"
}

# 绘制
p_count_top10 <- ggplot(stack_data_top10, aes(x = cn_detailed, y = count, fill = detailed_grouped)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  scale_fill_manual(values = cell_colors, name = "Cell Type") +
  labs(title = "Cell Type Count by CN (Top 10 per CN + Others)",
       x = "Cell Neighborhood (CN)",
       y = "Cell Count") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "right",
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.4, "cm")
  )

ggsave("mPT_cn_celltype_count_top10_per_CN.pdf", p_count_top10, width = 16, height = 8, dpi = 300)

#Tumor 的CN组成
tumor_cell_types <- grep("Tumor_", unique(mPT@meta.data$detailed), value = TRUE)

if(length(tumor_cell_types) == 0) {
  tumor_cell_types <- grep("tumor|malignant|cancer", 
                           unique(mPT@meta.data$detailed), 
                           value = TRUE, ignore.case = TRUE)
}

print(paste("找到", length(tumor_cell_types), "种肿瘤亚型："))
print(tumor_cell_types)

# =====================
# 计算每种肿瘤亚型在各CN中的占比
# =====================

tumor_cn_composition <- mPT@meta.data %>%
  filter(detailed %in% tumor_cell_types) %>%
  group_by(detailed, cellcharter_cluster) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(detailed) %>%
  mutate(proportion = count / sum(count) * 100) %>%
  ungroup() %>%
  filter(proportion > 0)

# =====================
# 确保 CN 按数字顺序排序 (CN1 到 CN15)
# =====================

# 获取所有出现的CN编号
all_cn_numbers <- sort(unique(tumor_cn_composition$cellcharter_cluster))

# 创建 CN 标签（CN1, CN2, ...）
tumor_cn_composition$cn_label <- factor(
  paste0("CN", tumor_cn_composition$cellcharter_cluster),
  levels = paste0("CN", all_cn_numbers)  # 按数字顺序设置级别
)

# =====================
# 定义颜色（CN1-CN15）
# =====================

cluster_colors <- c(
  "#A0CBE8", "#FFBE7D", "#8CD17D", "#86BCB6", "#FF9D9A",
  "#FABFD2", "#D4A6C8", "#F1CE63", "#D7B5A6", "#B07AA1",
  "#BAB0AC", "#4E79A7", "#F28E2B", "#59A14F", "#E15759"
)

# 只取前15个颜色（对应CN1-CN15）
cluster_colors_15 <- cluster_colors[1:15]
names(cluster_colors_15) <- paste0("CN", 1:15)

# 只保留实际出现的CN的颜色
cn_colors <- cluster_colors_15[paste0("CN", all_cn_numbers)]

# =====================
# 绘制分面饼图
# =====================

p_facet_pie <- ggplot(tumor_cn_composition, aes(x = "", y = proportion, fill = cn_label)) +
  geom_bar(stat = "identity", width = 1, color = "black", size = 0.5) +
  coord_polar("y", start = 0) +
  facet_wrap(~ detailed, ncol =5) +
  scale_fill_manual(values = cn_colors, name = "CN", breaks = paste0("CN", 1:15)) +
  labs(title = "mPT: CN Composition for Each Tumor Subtype") +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.4, "cm")
  )

print(p_facet_pie)
ggsave("mPT_tumor_subtype_facet_pie.pdf", p_facet_pie, 
       width = 14, height = ceiling(length(tumor_cell_types)/3) * 4, 
       dpi = 300)

#D3
library(ggplot2)
library(dplyr)

# 提取D3样本中FOV 164和165的数据
selected_fovs <- c(155, 159)
fov_data <- mPT@meta.data %>% 
  filter(sample == "D3", fov %in% selected_fovs)



# 获取实际出现的CN注释，按首字母排序
actual_cn <- sort(unique(fov_data$cn_detailed))
print("按首字母排序的CN：")
print(actual_cn)

# 颜色向量（按顺序分配）
cluster_colors <- c(
  "#A0CBE8", "#FFBE7D", "#8CD17D", "#86BCB6", "#FF9D9A",
  "#FABFD2", "#D4A6C8", "#F1CE63", "#D7B5A6", "#B07AA1",
  "#BAB0AC", "#4E79A7", "#F28E2B", "#59A14F", "#E15759",
  "#499894", "#D37295", "#B6992D", "#9D7660", "#79706E"
)

# 按顺序分配颜色
used_colors <- cluster_colors[1:length(actual_cn)]
names(used_colors) <- actual_cn


# 绘制
p <- ggplot(fov_data, aes(x = CenterX_global_px, y = CenterY_global_px, color = cn_detailed)) +
  geom_point(size = 0.5, alpha = 0.8) +  
  scale_color_manual(values = used_colors, name = "Cell Neighborhood") +
  labs(title = "mPT D3 Sample: FOV 153 & 159 CN Distribution",
       x = "X Coordinate (px)",
       y = "Y Coordinate (px)") +
  theme_bw() +
  theme(
    aspect.ratio = 1.5,
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    legend.position = "right",
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    legend.key.size = unit(0.6, "cm"),
    legend.text = element_text(size = 7)
  ) +
  guides(color = guide_legend(override.aes = list(size = 3), ncol = 1))

ggsave("mPT_D3_FOV153_159_CN.pdf", p, width = 6, height = 8, dpi = 300)


#fov 范围
fov255_data <- mPT@meta.data %>% 
  filter(sample == "C2", fov == 255)

# 查看X坐标范围
range_x <- range(fov255_data$CenterX_global_px)
print(paste("X轴范围:", range_x[1], "~", range_x[2]))

# 查看完整范围信息
summary(fov255_data$CenterX_global_px)

library(ggplot2)
library(dplyr)

# =====================
# 1. 自动提取所有CN注释并分配固定颜色
# =====================

# 从mPT对象中自动提取所有唯一的CN注释
all_cn_annotations <- sort(unique(mPT@meta.data$cn_detailed))
print("所有CN注释：")
print(all_cn_annotations)

# 颜色向量（自动扩展）
cluster_colors <- c(
  "#A0CBE8", "#FFBE7D", "#8CD17D", "#86BCB6", "#FF9D9A",
  "#FABFD2", "#D4A6C8", "#F1CE63", "#D7B5A6", "#B07AA1",
  "#BAB0AC", "#4E79A7", "#F28E2B", "#59A14F", "#E15759",
  "#499894", "#D37295", "#B6992D", "#9D7660", "#79706E"
)

# 自动分配颜色（循环使用颜色向量直到足够）
fixed_colors <- rep(cluster_colors, length.out = length(all_cn_annotations))
names(fixed_colors) <- all_cn_annotations

print("颜色映射：")
print(fixed_colors)

# =====================
# 2. 提取C2样本FOV 255和261的数据，并筛选X轴范围
# =====================

selected_fovs <- c(255, 254)
fov_data <- mPT@meta.data %>% 
  filter(sample == "C2", fov %in% selected_fovs)
         #CenterX_global_px >= 49805, CenterX_global_px <= 54031)

# =====================
# 3. 绘图（颜色自动匹配）
# =====================

p <- ggplot(fov_data, aes(x = CenterX_global_px, y = CenterY_global_px, color = cn_detailed)) +
  geom_point(size = 0.5, alpha = 0.8) +  
  scale_color_manual(values = fixed_colors, name = "Cell Neighborhood") +
  labs(title = "mPT C2 Sample: FOV 254 & 255 CN Distribution",
       x = "X Coordinate (px)",
       y = "Y Coordinate (px)") +
  theme_bw() +
  theme(
    aspect.ratio = 0.5,
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    legend.position = "right",
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    legend.key.size = unit(0.6, "cm"),
    legend.text = element_text(size = 7)
  ) +
  guides(color = guide_legend(override.aes = list(size = 3), ncol = 1))

ggsave("mPT_C2_FOV254_255_CN.pdf", p, width = 8, height = 4, dpi = 300)

#CN4 B and T
library(ggplot2)
library(dplyr)
library(patchwork)

# ===================== 提取CN4中的B细胞和T细胞亚群 =====================

# 筛选CN4的细胞
cn4_cells <- mPT@meta.data %>% filter(cn_detailed == "TLS_Zone")

# 提取B细胞亚群（以B_开头）
b_cell_types <- grep("^B_", unique(cn4_cells$detailed), value = TRUE)

# 提取T细胞亚群（以T_开头，排除Tumor）
t_cell_types <- grep("^T_", unique(cn4_cells$detailed), value = TRUE)
t_cell_types <- t_cell_types[!grepl("Tumor", t_cell_types)]

cat("B细胞亚群：", paste(b_cell_types, collapse = ", "), "\n")
cat("T细胞亚群：", paste(t_cell_types, collapse = ", "), "\n")

# ===================== 颜色方案（都从第一个颜色开始） =====================

cluster_colors <- c(
  "#A0CBE8", "#FFBE7D", "#8CD17D", "#86BCB6", "#FF9D9A",
  "#FABFD2", "#D4A6C8", "#F1CE63", "#D7B5A6", "#B07AA1",
  "#BAB0AC", "#4E79A7", "#F28E2B", "#59A14F", "#E15759",
  "#499894", "#D37295", "#B6992D", "#9D7660", "#79706E"
)

# ===================== 计算B细胞亚群占比 =====================

b_data <- cn4_cells %>%
  filter(detailed %in% b_cell_types) %>%
  group_by(detailed) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(percentage = count / sum(count) * 100)

# 按百分比排序
b_data <- b_data %>% arrange(desc(percentage))
b_data$detailed <- factor(b_data$detailed, levels = b_data$detailed)

# 分配颜色（从第一个开始）
b_colors <- cluster_colors[1:nrow(b_data)]
names(b_colors) <- b_data$detailed

# ===================== 计算T细胞亚群占比 =====================

t_data <- cn4_cells %>%
  filter(detailed %in% t_cell_types) %>%
  group_by(detailed) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(percentage = count / sum(count) * 100)

# 按百分比排序
t_data <- t_data %>% arrange(desc(percentage))
t_data$detailed <- factor(t_data$detailed, levels = t_data$detailed)

# 分配颜色（也从第一个开始）
t_colors <- cluster_colors[1:nrow(t_data)]
names(t_colors) <- t_data$detailed

# B细胞饼图（不显示比例数字）
if (nrow(b_data) > 0) {
  p_b <- ggplot(b_data, aes(x = "", y = percentage, fill = detailed)) +
    geom_bar(stat = "identity", width = 1, color = "black", size = 0.5) +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = b_colors, name = "B Cell Subtype") +
    labs(title = paste0("B Cell Composition in CN4")) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      legend.position = "right",
      legend.text = element_text(size = 8)
    )
  
  ggsave("CN4_B_cell_pie.pdf", p_b, width = 7, height = 5, dpi = 300)
  print(p_b)
}

# T细胞饼图（不显示比例数字）
if (nrow(t_data) > 0) {
  p_t <- ggplot(t_data, aes(x = "", y = percentage, fill = detailed)) +
    geom_bar(stat = "identity", width = 1, color = "black", size = 0.5) +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = t_colors, name = "T Cell Subtype") +
    labs(title = paste0("T Cell Composition in CN4")) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      legend.position = "right",
      legend.text = element_text(size = 8)
    )
  
  ggsave("CN4_T_cell_pie.pdf", p_t, width = 7, height = 5, dpi = 300)
  print(p_t)
}

# ===================== 两个饼图拼在一起 =====================

if (nrow(b_data) > 0 && nrow(t_data) > 0) {
  combined <- p_b + p_t +
    plot_annotation(
      title = "CN4 (TLS Zone) Immune Cell Composition",
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
    )
  
  ggsave("CN4_B_T_cell_pies_combined.pdf", combined, width = 12, height = 5, dpi = 300)
  print(combined)
}

```
# nmPT subtype CN
```R
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
nmPT <- readRDS("nmPT_cellcharter_FINAL.rds")

#CN组成
# 查看每个CN的细胞组成
cn_composition <- table(nmPT@meta.data$cellcharter_cluster, 
                        nmPT@meta.data$detailed)

# 计算每个CN的比例
cn_prop <- prop.table(cn_composition, margin = 1) * 100

# 打印每个CN的主要细胞类型（前5种）
for (cn in rownames(cn_prop)) {
  cat("\n=== CN", cn, "===\n")
  top10 <- sort(cn_prop[cn, ], decreasing = TRUE)[1:10]
  for (i in 1:length(top10)) {
    cat(names(top10)[i], ": ", round(top10[i], 1), "%\n", sep = "")
  }
}

#CN注释
cn_annotation_detailed <- c(
  "1"  = "Vascular_Stromal_Zone",
  "2"  = "Immune_Zone",
  "3"  = "Macrophage_Dominant_Zone_B",
  "4"  = "Macrophage_Dominant_Zone_C",
  "5"  = "TLS_Zone",
  "6"  = "Stromal_Remodeling_Zone",
  "7"  = "Stem_Proliferative_Tumor_Zone_A",
  "8"  = "Alveolar_Homeostatic_Zone",
  "9"  = "Stem_Proliferative_Tumor_Zone_B",
  "10" = "Mast_Cell_Dominant_Zone",
  "11" = "Plasma_Enriched_Zone_B",
  "12" = "Immunogenic_Tumor_Core_A",
  "13" = "Immunogenic_Tumor_Core_B",
  "14" = "Immunogenic_Tumor_Mast_Interface",
  "15" = "Plasma_Enriched_Zone_D",
  "16" = "Immunogenic_Tumor_Immune_Interface",
  "17" = "EMT_Tumor_Core",
  "18" = "Immunogenic_Tumor_Core_C",
  "19" = "Plasma_Enriched_Zone_E",
  "20" = "Stem_Proliferative_Tumor_Zone_C"
)

# 应用到 Seurat 对象
nmPT@meta.data$cn_detailed <- cn_annotation_detailed[as.character(nmPT@meta.data$cellcharter_cluster)]

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# 定义颜色（你的20种颜色）
cluster_colors <- c(
  "#A0CBE8", "#FFBE7D", "#8CD17D", "#86BCB6", "#FF9D9A",
  "#FABFD2", "#D4A6C8", "#F1CE63", "#D7B5A6", "#B07AA1",
  "#BAB0AC", "#4E79A7", "#F28E2B", "#59A14F", "#E15759",
  "#499894", "#D37295", "#B6992D", "#9D7660", "#79706E"
)

# =====================
# 1. 准备数据
# =====================

# 假设你的 nmPT 对象中有：
# - cn_detailed: CN编号 (1-15)
# - sample: 样本名称

# 计算总体比例
total_prop <- nmPT@meta.data %>%
  group_by(cn_detailed) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(proportion = count / sum(count) * 100,
         type = "Overall")

# 计算各样本比例
sample_prop <- nmPT@meta.data %>%
  group_by(sample, cn_detailed) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(sample) %>%
  mutate(proportion = count / sum(count) * 100,
         type = sample) %>%
  ungroup()

# 合并数据
combined_data <- bind_rows(total_prop, sample_prop)

# 确保 cn_detailed 是因子（按数字顺序）
combined_data$cn_detailed <- factor(combined_data$cn_detailed, 
                                     levels = sort(unique(combined_data$cn_detailed)))

# =====================
# 方法1：分面图（推荐）
# =====================

# 总体 + 各样本分面显示
p_facet <- ggplot(combined_data, aes(x = type, y = proportion, fill = cn_detailed)) +
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  scale_fill_manual(values = cluster_colors, name = "Cell Neighborhood") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "CN Proportion by Sample",
       x = "",
       y = "Proportion") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "right",
    legend.text = element_text(size = 8),
    strip.background = element_rect(fill = "lightgray"),
    strip.text = element_text(size = 10, face = "bold")
  )

ggsave("nmPT_cn_proportion_by_sample_facet.pdf", p_facet, width = 6, height = 6, dpi = 300)

# 3. 每个邻域细胞组成热图
# =====================
message("=== 3. 生成邻域细胞组成热图 ===")

cluster_composition <- table(nmPT@meta.data$cellcharter_cluster, 
                           nmPT@meta.data$detailed)

cluster_composition_df <- as.data.frame(cluster_composition)
colnames(cluster_composition_df) <- c("Cluster", "CellType", "Count")

cluster_totals <- aggregate(Count ~ Cluster, data = cluster_composition_df, sum)
cluster_composition_df <- merge(cluster_composition_df, cluster_totals, by = "Cluster")
cluster_composition_df$Proportion <- cluster_composition_df$Count.x / cluster_composition_df$Count.y

# 转换为宽格式
heatmap_data <- cluster_composition_df %>%
  select(Cluster, CellType, Proportion) %>%
  pivot_wider(names_from = CellType, values_from = Proportion, values_fill = 0)

# 按Cluster数字排序
cluster_numeric <- as.numeric(as.character(heatmap_data$Cluster))
heatmap_data <- heatmap_data[order(cluster_numeric), ]

# 创建行标签（使用你的注释名称）
row_labels_annotated <- cn_annotation_detailed[as.character(cluster_numeric[order(cluster_numeric)])]

# 设置行名
rownames(heatmap_data) <- row_labels_annotated
heatmap_matrix <- as.matrix(heatmap_data[, -1])

# 绘制热图（字体增大的版本）
p_heatmap <- pheatmap(heatmap_matrix,
                    cluster_rows = TRUE,
                    cluster_cols = TRUE,
                    display_numbers = FALSE,
                    main = "Cell Type Composition by Spatial Domain",
                    fontsize = 14,           # 全局字体增大
                    fontsize_row = 10,       # 行标签字体（注释名称较长，用10合适）
                    fontsize_col = 9,        # 列标签字体
                    labels_row = row_labels_annotated,
                    color = colorRampPalette(c("#f0f0f0", "#3182bd", "#08519c"))(100),
                    cellwidth = 11,          # 单元格宽度
                    cellheight = 11)         # 单元格高度

# 保存
pdf("nmPT_CN_cell_type_composition_heatmap_annotated.pdf", width = 16, height = 12)
print(p_heatmap)
dev.off()

library(ggplot2)
library(dplyr)
library(scales)
library(RColorBrewer)

# 计算每个 CN 中每种细胞类型的数量
stack_data_count <- nmPT@meta.data %>%
  group_by(cn_detailed, detailed) %>%
  summarise(count = n(), .groups = "drop")

# 为每个 CN 找出 top 10 细胞类型
top10_per_cn <- stack_data_count %>%
  group_by(cn_detailed) %>%
  arrange(desc(count)) %>%
  slice_head(n = 10) %>%
  ungroup()

# 将非 top 10 的细胞类型合并为 "Other"
stack_data_top10 <- stack_data_count %>%
  left_join(top10_per_cn %>% select(cn_detailed, detailed) %>% mutate(is_top10 = TRUE),
            by = c("cn_detailed", "detailed")) %>%
  mutate(detailed_grouped = ifelse(is_top10, detailed, "Other")) %>%
  group_by(cn_detailed, detailed_grouped) %>%
  summarise(count = sum(count), .groups = "drop")

# 按首字母顺序排序 CN
cn_order <- sort(unique(stack_data_top10$cn_detailed))
stack_data_top10$cn_detailed <- factor(stack_data_top10$cn_detailed, levels = cn_order)

# 获取所有需要着色的类别
all_cell_types <- unique(stack_data_top10$detailed_grouped)
n_colors_needed <- length(all_cell_types)

# 方案1：使用彩虹色（qualitative调色板）
# 提取足够多的颜色
if(n_colors_needed <= 12) {
  cell_colors <- brewer.pal(n_colors_needed, "Set3")
} else {
  # 使用 colorRampPalette 生成更多颜色
  get_palette <- colorRampPalette(brewer.pal(12, "Set3"))
  cell_colors <- get_palette(n_colors_needed)
}
names(cell_colors) <- all_cell_types

# 确保 "Other" 是灰色
if("Other" %in% names(cell_colors)) {
  cell_colors["Other"] <- "#D3D3D3"
}

# 绘制
p_count_top10 <- ggplot(stack_data_top10, aes(x = cn_detailed, y = count, fill = detailed_grouped)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  scale_fill_manual(values = cell_colors, name = "Cell Type") +
  labs(title = "Cell Type Count by CN (Top 10 per CN + Others)",
       x = "Cell Neighborhood (CN)",
       y = "Cell Count") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "right",
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.4, "cm")
  )

ggsave("nmPT_cn_celltype_count_top10_per_CN.pdf", p_count_top10, width = 16, height = 8, dpi = 300)

#比较mPT和nmPT的immune zone
library(ggplot2)
library(dplyr)
library(patchwork)

# ===================== 提取两个数据集中 Immune_Zone 的细胞 =====================

# mPT中的Immune_Zone
mPT_immune <- mPT@meta.data %>% filter(cn_detailed == "Immune_Zone")
mPT_immune$dataset <- "mPT"

# nmPT中的Immune_Zone
nmPT_immune <- nmPT@meta.data %>% filter(cn_detailed == "Immune_Zone")
nmPT_immune$dataset <- "nmPT"

# 合并数据
combined_immune <- bind_rows(mPT_immune, nmPT_immune)

cat("mPT Immune_Zone 细胞数：", nrow(mPT_immune), "\n")
cat("nmPT Immune_Zone 细胞数：", nrow(nmPT_immune), "\n")

# ===================== 1. 总细胞组成占比差异堆叠图 =====================

# 计算每个数据集中各细胞类型的占比
total_composition <- combined_immune %>%
  group_by(dataset, detailed) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(dataset) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ungroup()

# 获取top 15细胞类型（用于清晰显示）
top_cell_types <- total_composition %>%
  group_by(detailed) %>%
  summarise(total_count = sum(count)) %>%
  arrange(desc(total_count)) %>%
  slice_head(n = 15) %>%
  pull(detailed)

# 其他细胞类型合并为 "Other"
total_composition_plot <- total_composition %>%
  mutate(cell_group = ifelse(detailed %in% top_cell_types, detailed, "Other")) %>%
  group_by(dataset, cell_group) %>%
  summarise(percentage = sum(percentage), .groups = "drop")

# 颜色方案
cluster_colors <- c(
  "#A0CBE8", "#FFBE7D", "#8CD17D", "#86BCB6", "#FF9D9A",
  "#FABFD2", "#D4A6C8", "#F1CE63", "#D7B5A6", "#B07AA1",
  "#BAB0AC", "#4E79A7", "#F28E2B", "#59A14F", "#E15759",
  "#499894", "#D37295", "#B6992D", "#9D7660", "#79706E", "#D3D3D3"
)

# 绘制堆叠图
p_total <- ggplot(total_composition_plot, aes(x = dataset, y = percentage, fill = cell_group)) +
  geom_bar(stat = "identity", position = "stack", width = 0.6) +
  scale_fill_manual(values = cluster_colors[1:length(unique(total_composition_plot$cell_group))], 
                    name = "Cell Type") +
  labs(title = "Immune_Zone: Cell Composition Comparison",
       subtitle = paste("mPT (n=", nrow(mPT_immune), ") vs nmPT (n=", nrow(nmPT_immune), ")"),
       x = "",
       y = "Percentage (%)") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 11),
    legend.position = "right",
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.4, "cm")
  )

ggsave("Immune_Zone_composition_comparison_total.pdf", p_total, width = 12, height = 6, dpi = 300)
print(p_total)

# ===================== 2. T细胞亚群占比差异堆叠图 =====================

# 提取T细胞（以 T_ 开头，排除 Tumor）
t_cell_types <- grep("^T_", unique(combined_immune$detailed), value = TRUE)
t_cell_types <- t_cell_types[!grepl("Tumor", t_cell_types)]

t_composition <- combined_immune %>%
  filter(detailed %in% t_cell_types) %>%
  group_by(dataset, detailed) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(dataset) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ungroup()

if (nrow(t_composition) > 0) {
  p_t <- ggplot(t_composition, aes(x = dataset, y = percentage, fill = detailed)) +
    geom_bar(stat = "identity", position = "stack", width = 0.6) +
    scale_fill_manual(values = cluster_colors[1:length(unique(t_composition$detailed))], 
                      name = "T Cell Subtype") +
    labs(title = "Immune_Zone: T Cell Subtype Comparison",
         x = "",
         y = "Percentage (%)") +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text.x = element_text(size = 12, face = "bold"),
      axis.text.y = element_text(size = 10),
      axis.title.y = element_text(size = 11),
      legend.position = "right",
      legend.text = element_text(size = 8),
      legend.key.size = unit(0.4, "cm")
    )
  
  ggsave("Immune_Zone_composition_comparison_T_cells.pdf", p_t, width = 5, height = 5, dpi = 300)
  print(p_t)
} else {
  cat("没有找到T细胞\n")
}

#比较mPt和nmPT的TLS区
library(ggplot2)
library(dplyr)
library(patchwork)

# ===================== 提取两个数据集中 TLS_Zone 的细胞 =====================

# mPT中的TLS_Zone
mPT_tls <- mPT@meta.data %>% filter(cn_detailed == "TLS_Zone")
mPT_tls$dataset <- "mPT"

# nmPT中的TLS_Zone
nmPT_tls <- nmPT@meta.data %>% filter(cn_detailed == "TLS_Zone")
nmPT_tls$dataset <- "nmPT"

# 合并数据
combined_tls <- bind_rows(mPT_tls, nmPT_tls)

cat("mPT TLS_Zone 细胞数：", nrow(mPT_tls), "\n")
cat("nmPT TLS_Zone 细胞数：", nrow(nmPT_tls), "\n")

# ===================== 颜色方案 =====================
cluster_colors <- c(
  "#A0CBE8", "#FFBE7D", "#8CD17D", "#86BCB6", "#FF9D9A",
  "#FABFD2", "#D4A6C8", "#F1CE63", "#D7B5A6", "#B07AA1",
  "#BAB0AC", "#4E79A7", "#F28E2B", "#59A14F", "#E15759",
  "#499894", "#D37295", "#B6992D", "#9D7660", "#79706E", "#D3D3D3"
)

# ===================== 1. TLS_Zone 总细胞组成占比差异堆叠图 =====================

# 计算每个数据集中各细胞类型的占比
total_composition <- combined_tls %>%
  group_by(dataset, detailed) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(dataset) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ungroup()

# 获取top 15细胞类型
top_cell_types <- total_composition %>%
  group_by(detailed) %>%
  summarise(total_count = sum(count)) %>%
  arrange(desc(total_count)) %>%
  slice_head(n = 15) %>%
  pull(detailed)

# 其他细胞类型合并为 "Other"
total_composition_plot <- total_composition %>%
  mutate(cell_group = ifelse(detailed %in% top_cell_types, detailed, "Other")) %>%
  group_by(dataset, cell_group) %>%
  summarise(percentage = sum(percentage), .groups = "drop")

p_total <- ggplot(total_composition_plot, aes(x = dataset, y = percentage, fill = cell_group)) +
  geom_bar(stat = "identity", position = "stack", width = 0.6) +
  scale_fill_manual(values = cluster_colors[1:length(unique(total_composition_plot$cell_group))], 
                    name = "Cell Type") +
  labs(title = "TLS_Zone: Cell Composition Comparison",
       subtitle = paste("mPT (n=", nrow(mPT_tls), ") vs nmPT (n=", nrow(nmPT_tls), ")"),
       x = "",
       y = "Percentage (%)") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 11),
    legend.position = "right",
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.4, "cm")
  )

ggsave("TLS_Zone_composition_comparison_total.pdf", p_total, width = 12, height = 6, dpi = 300)
print(p_total)

# ===================== 2. TLS_Zone T细胞亚群占比差异堆叠图 =====================

# 提取T细胞（以 T_ 开头，排除 Tumor）
t_cell_types <- grep("^T_", unique(combined_tls$detailed), value = TRUE)
t_cell_types <- t_cell_types[!grepl("Tumor", t_cell_types)]

t_composition <- combined_tls %>%
  filter(detailed %in% t_cell_types) %>%
  group_by(dataset, detailed) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(dataset) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ungroup()

if (nrow(t_composition) > 0) {
  p_t <- ggplot(t_composition, aes(x = dataset, y = percentage, fill = detailed)) +
    geom_bar(stat = "identity", position = "stack", width = 0.6) +
    scale_fill_manual(values = cluster_colors[1:length(unique(t_composition$detailed))], 
                      name = "T Cell Subtype") +
    labs(title = "TLS_Zone: T Cell Subtype Comparison",
         x = "",
         y = "Percentage (%)") +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text.x = element_text(size = 12, face = "bold"),
      axis.text.y = element_text(size = 10),
      axis.title.y = element_text(size = 11),
      legend.position = "right",
      legend.text = element_text(size = 8),
      legend.key.size = unit(0.4, "cm")
    )
  
  ggsave("TLS_Zone_composition_comparison_T_cells.pdf", p_t, width = 5, height = 5, dpi = 300)
  print(p_t)
} else {
  cat("没有找到T细胞\n")
}

# ===================== 3. TLS_Zone B细胞亚群占比差异堆叠图 =====================

# 提取B细胞（以 B_ 开头）
b_cell_types <- grep("^B_", unique(combined_tls$detailed), value = TRUE)

b_composition <- combined_tls %>%
  filter(detailed %in% b_cell_types) %>%
  group_by(dataset, detailed) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(dataset) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ungroup()

if (nrow(b_composition) > 0) {
  p_b <- ggplot(b_composition, aes(x = dataset, y = percentage, fill = detailed)) +
    geom_bar(stat = "identity", position = "stack", width = 0.6) +
    scale_fill_manual(values = cluster_colors[1:length(unique(b_composition$detailed))], 
                      name = "B Cell Subtype") +
    labs(title = "TLS_Zone: B Cell Subtype Comparison",
         x = "",
         y = "Percentage (%)") +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text.x = element_text(size = 12, face = "bold"),
      axis.text.y = element_text(size = 10),
      axis.title.y = element_text(size = 11),
      legend.position = "right",
      legend.text = element_text(size = 8),
      legend.key.size = unit(0.4, "cm")
    )
  
  ggsave("TLS_Zone_composition_comparison_B_cells.pdf", p_b, width = 5, height = 5, dpi = 300)
  print(p_b)
} else {
  cat("没有找到B细胞\n")
}

#A1肿瘤
library(ggplot2)
library(dplyr)

# ===================== 提取 nmPT 中 A1 样本的数据 =====================

a1_data <- nmPT@meta.data %>% filter(sample == "A1")

# 筛选指定坐标范围
a1_subset <- a1_data %>%
  filter(
    CenterX_global_px >= 94000 & CenterX_global_px <= 108000,
    CenterY_global_px >= 5000 & CenterY_global_px <= 18000
  )

cat("筛选后细胞数：", nrow(a1_subset), "\n")

# ===================== 获取实际出现的CN =====================
actual_cn <- sort(unique(a1_subset$cn_detailed))
cat("实际出现的CN：\n")
print(actual_cn)

# ===================== 定义颜色 =====================
cluster_colors <- c(
  "#A0CBE8", "#FFBE7D", "#8CD17D", "#86BCB6", "#FF9D9A",
  "#FABFD2", "#D4A6C8", "#F1CE63", "#D7B5A6", "#B07AA1",
  "#BAB0AC", "#4E79A7", "#F28E2B", "#59A14F", "#E15759",
  "#499894", "#D37295", "#B6992D", "#9D7660", "#79706E"
)

# 将颜色直接赋给CN注释名称
used_colors <- cluster_colors[1:length(actual_cn)]
names(used_colors) <- actual_cn

# ===================== 绘制点图 =====================

p <- ggplot(a1_subset, aes(x = CenterX_global_px, y = CenterY_global_px, color = cn_detailed)) +
  geom_point(size = 0.5, alpha = 0.8) +
  scale_color_manual(values = used_colors, name = "Cell Neighborhood (CN)") +
  coord_fixed() +
  labs(
    title = "nmPT A1 Sample: Spatial Distribution of Cell Neighborhoods",
    subtitle = paste("X: 94000-108000, Y: 5000-18000 | Total cells:", nrow(a1_subset)),
    x = "X Coordinate (px)",
    y = "Y Coordinate (px)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.8, "cm"),
    aspect.ratio = 1
  ) +
   guides(color = guide_legend(override.aes = list(size = 4)))

ggsave("nmPT_A1_CN_spatial_zoom.pdf", p, width = 12, height = 8, dpi = 300)
print(p)

selected_fovs <- c(280, 281, 282, 286, 287, 288, 292, 293, 294)

b2_data <- nmPT@meta.data %>% 
  filter(sample == "B2", fov %in% selected_fovs)

# 筛选X轴范围
b2_data <- b2_data %>%
  filter(CenterX_global_px >= 73060 & CenterX_global_px <= 85024)

cat("筛选后细胞数：", nrow(b2_data), "\n")
cat("X轴范围：", range(b2_data$CenterX_global_px), "\n")
cat("Y轴范围：", range(b2_data$CenterY_global_px), "\n")

# ===================== 获取实际出现的CN =====================
actual_cn <- sort(unique(b2_data$cn_detailed))
cat("实际出现的CN：\n")
print(actual_cn)

# ===================== 定义颜色 =====================
cluster_colors <- c(
  "#A0CBE8", "#FFBE7D", "#8CD17D", "#86BCB6", "#FF9D9A",
  "#FABFD2", "#D4A6C8", "#F1CE63", "#D7B5A6", "#B07AA1",
  "#BAB0AC", "#4E79A7", "#F28E2B", "#59A14F", "#E15759",
  "#499894", "#D37295", "#B6992D", "#9D7660", "#79706E"
)

# 将颜色直接赋给CN注释名称
used_colors <- cluster_colors[1:length(actual_cn)]
names(used_colors) <- actual_cn

# ===================== 绘制点图（所有FOV在一张图上） =====================

p <- ggplot(b2_data, aes(x = CenterX_global_px, y = CenterY_global_px, color = cn_detailed)) +
  geom_point(size = 0.3, alpha = 0.8) +
  scale_color_manual(values = used_colors, name = "Cell Neighborhood (CN)") +
  coord_fixed(ratio = 1) +
  labs(
    title = "nmPT B2 Sample: Spatial Distribution of Cell Neighborhoods",
    subtitle = paste("FOV:", paste(selected_fovs, collapse = ", "), 
                     "| X: 73060-85024 | Total cells:", nrow(b2_data)),
    x = "X Coordinate (px)",
    y = "Y Coordinate (px)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
    plot.subtitle = element_text(hjust = 0.5, size = 8),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.8, "cm"),
    aspect.ratio = 1  # 强制绘图区域为正方形
  ) +
  guides(color = guide_legend(override.aes = list(size = 4)))

ggsave("nmPT_B2_selected_FOVs_CN_spatial.pdf", p, width = 9, height = 6, dpi = 300)

#Fib 的CN组成
fib_cell_types <- grep("Fib_", unique(nmPT@meta.data$detailed), value = TRUE)


fib_cell_types <- grep("Fib", 
                           unique(nmPT@meta.data$detailed), 
                           value = TRUE, ignore.case = TRUE)

print(paste("找到", length(fib_cell_types), "种肿瘤亚型："))
print(fib_cell_types)

# =====================
# 计算每种肿瘤亚型在各CN中的占比
# =====================

fib_cn_composition <- nmPT@meta.data %>%
  filter(detailed %in% fib_cell_types) %>%
  group_by(detailed, cellcharter_cluster) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(detailed) %>%
  mutate(proportion = count / sum(count) * 100) %>%
  ungroup() %>%
  filter(proportion > 0)

# =====================
# 确保 CN 按数字顺序排序 (CN1 到 CN15)
# =====================

# 获取所有出现的CN编号
all_cn_numbers <- sort(unique(fib_cn_composition$cellcharter_cluster))

# 创建 CN 标签（CN1, CN2, ...）
fib_cn_composition$cn_label <- factor(
  paste0("CN", fib_cn_composition$cellcharter_cluster),
  levels = paste0("CN", all_cn_numbers)  # 按数字顺序设置级别
)

# =====================
# 定义颜色（CN1-CN15）
# =====================

cluster_colors <- c(
  "#A0CBE8", "#FFBE7D", "#8CD17D", "#86BCB6", "#FF9D9A",
  "#FABFD2", "#D4A6C8", "#F1CE63", "#D7B5A6", "#B07AA1",
  "#BAB0AC", "#4E79A7", "#F28E2B", "#59A14F", "#E15759"
)

# 只取前15个颜色（对应CN1-CN15）
cluster_colors_20 <- cluster_colors[1:20]
names(cluster_colors_20) <- paste0("CN", 1:20)

# 只保留实际出现的CN的颜色
cn_colors <- cluster_colors_20[paste0("CN", all_cn_numbers)]

# =====================
# 绘制分面饼图
# =====================

p_facet_pie <- ggplot(fib_cn_composition, aes(x = "", y = proportion, fill = cn_label)) +
  geom_bar(stat = "identity", width = 1, color = "black", size = 0.5) +
  coord_polar("y", start = 0) +
  facet_wrap(~ detailed, ncol =4) +
  scale_fill_manual(values = cn_colors, name = "CN", breaks = paste0("CN", 1:15)) +
  labs(title = "nmPT: CN Composition for Each Fibroblasts Subtype") +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.4, "cm")
  )

print(p_facet_pie)
ggsave("nmPT_fib_subtype_facet_pie.pdf", p_facet_pie, 
       width = 14, height = ceiling(length(fib_cell_types)/3) * 4, 
       dpi = 300)
```
# metLN subtype CN
```R
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
metLN <- readRDS("metLN_cellcharter_FINAL.rds")

#CN组成
# 查看每个CN的细胞组成
cn_composition <- table(metLN@meta.data$cellcharter_cluster, 
                        metLN@meta.data$detailed)

# 计算每个CN的比例
cn_prop <- prop.table(cn_composition, margin = 1) * 100

# 打印每个CN的主要细胞类型（前5种）
for (cn in rownames(cn_prop)) {
  cat("\n=== CN", cn, "===\n")
  top10 <- sort(cn_prop[cn, ], decreasing = TRUE)[1:10]
  for (i in 1:length(top10)) {
    cat(names(top10)[i], ": ", round(top10[i], 1), "%\n", sep = "")
  }
}

#CN注释
cn_annotation_detailed <- c(
  "1"  = "B_Cell_Follicle_Zone",
  "2"  = "Plasma_Macrophage_Zone",
  "3"  = "Plasma_Enriched_Zone_F",
  "4"  = "Vascular_Stromal_Zone",
  "5"  = "FRC-rich_Lymphoid_Zone",
  "6"  = "Mast_Cell_Immune_Zone",
  "7"  = "Immune-inflamed_Tumor_Core_A",
  "8"  = "B_T_Cell_Zone",
  "9"  = "Immune-inflamed_Tumor_Core_B",
  "10" = "Macrophage_Dominant_Zone_D",
  "11" = "Immune-inflamed_Tumor_Macrophage_Interface",
  "12" = "Plasma_Enriched_Zone_G",
  "13" = "Macrophage_Dominant_Zone_E",
  "14" = "B_Cell_Immune-inflamed_Tumor_Interface"
)
# 应用到 Seurat 对象
metLN@meta.data$cn_detailed <- cn_annotation_detailed[as.character(metLN@meta.data$cellcharter_cluster)]

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# 定义颜色（你的20种颜色）
cluster_colors <- c(
  "#A0CBE8", "#FFBE7D", "#8CD17D", "#86BCB6", "#FF9D9A",
  "#FABFD2", "#D4A6C8", "#F1CE63", "#D7B5A6", "#B07AA1",
  "#BAB0AC", "#4E79A7", "#F28E2B", "#59A14F", "#E15759",
  "#499894", "#D37295", "#B6992D", "#9D7660", "#79706E"
)

# =====================
# 1. 准备数据
# =====================

# 假设你的 metLN 对象中有：
# - cn_detailed: CN编号 (1-15)
# - sample: 样本名称

# 计算总体比例
total_prop <- metLN@meta.data %>%
  group_by(cn_detailed) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(proportion = count / sum(count) * 100,
         type = "Overall")

# 计算各样本比例
sample_prop <- metLN@meta.data %>%
  group_by(sample, cn_detailed) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(sample) %>%
  mutate(proportion = count / sum(count) * 100,
         type = sample) %>%
  ungroup()

# 合并数据
combined_data <- bind_rows(total_prop, sample_prop)

# 确保 cn_detailed 是因子（按数字顺序）
combined_data$cn_detailed <- factor(combined_data$cn_detailed, 
                                     levels = sort(unique(combined_data$cn_detailed)))

# =====================
# 方法1：分面图（推荐）
# =====================

# 总体 + 各样本分面显示
p_facet <- ggplot(combined_data, aes(x = type, y = proportion, fill = cn_detailed)) +
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  scale_fill_manual(values = cluster_colors, name = "Cell Neighborhood") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "CN Proportion by Sample",
       x = "",
       y = "Proportion") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "right",
    legend.text = element_text(size = 8),
    strip.background = element_rect(fill = "lightgray"),
    strip.text = element_text(size = 10, face = "bold")
  )

ggsave("metLN_cn_proportion_by_sample_facet.pdf", p_facet, width = 6, height = 6, dpi = 300)

# 3. 每个邻域细胞组成热图
# =====================
message("=== 3. 生成邻域细胞组成热图 ===")

cluster_composition <- table(metLN@meta.data$cellcharter_cluster, 
                           metLN@meta.data$detailed)

cluster_composition_df <- as.data.frame(cluster_composition)
colnames(cluster_composition_df) <- c("Cluster", "CellType", "Count")

cluster_totals <- aggregate(Count ~ Cluster, data = cluster_composition_df, sum)
cluster_composition_df <- merge(cluster_composition_df, cluster_totals, by = "Cluster")
cluster_composition_df$Proportion <- cluster_composition_df$Count.x / cluster_composition_df$Count.y

# 转换为宽格式
heatmap_data <- cluster_composition_df %>%
  select(Cluster, CellType, Proportion) %>%
  pivot_wider(names_from = CellType, values_from = Proportion, values_fill = 0)

# 按Cluster数字排序
cluster_numeric <- as.numeric(as.character(heatmap_data$Cluster))
heatmap_data <- heatmap_data[order(cluster_numeric), ]

# 创建行标签（使用你的注释名称）
row_labels_annotated <- cn_annotation_detailed[as.character(cluster_numeric[order(cluster_numeric)])]

# 设置行名
rownames(heatmap_data) <- row_labels_annotated
heatmap_matrix <- as.matrix(heatmap_data[, -1])

# 绘制热图（字体增大的版本）
p_heatmap <- pheatmap(heatmap_matrix,
                    cluster_rows = TRUE,
                    cluster_cols = TRUE,
                    display_numbers = FALSE,
                    main = "Cell Type Composition by Spatial Domain",
                    fontsize = 14,           # 全局字体增大
                    fontsize_row = 10,       # 行标签字体（注释名称较长，用10合适）
                    fontsize_col = 9,        # 列标签字体
                    labels_row = row_labels_annotated,
                    color = colorRampPalette(c("#f0f0f0", "#3182bd", "#08519c"))(100),
                    cellwidth = 11,          # 单元格宽度
                    cellheight = 11)         # 单元格高度

# 保存
pdf("metLN_CN_cell_type_composition_heatmap_annotated.pdf", width = 16, height = 12)
print(p_heatmap)
dev.off()

library(ggplot2)
library(dplyr)
library(scales)
library(RColorBrewer)

# 计算每个 CN 中每种细胞类型的数量
stack_data_count <- metLN@meta.data %>%
  group_by(cn_detailed, detailed) %>%
  summarise(count = n(), .groups = "drop")

# 为每个 CN 找出 top 10 细胞类型
top10_per_cn <- stack_data_count %>%
  group_by(cn_detailed) %>%
  arrange(desc(count)) %>%
  slice_head(n = 10) %>%
  ungroup()

# 将非 top 10 的细胞类型合并为 "Other"
stack_data_top10 <- stack_data_count %>%
  left_join(top10_per_cn %>% select(cn_detailed, detailed) %>% mutate(is_top10 = TRUE),
            by = c("cn_detailed", "detailed")) %>%
  mutate(detailed_grouped = ifelse(is_top10, detailed, "Other")) %>%
  group_by(cn_detailed, detailed_grouped) %>%
  summarise(count = sum(count), .groups = "drop")

# 按首字母顺序排序 CN
cn_order <- sort(unique(stack_data_top10$cn_detailed))
stack_data_top10$cn_detailed <- factor(stack_data_top10$cn_detailed, levels = cn_order)

# 获取所有需要着色的类别
all_cell_types <- unique(stack_data_top10$detailed_grouped)
n_colors_needed <- length(all_cell_types)

# 方案1：使用彩虹色（qualitative调色板）
# 提取足够多的颜色
if(n_colors_needed <= 12) {
  cell_colors <- brewer.pal(n_colors_needed, "Set3")
} else {
  # 使用 colorRampPalette 生成更多颜色
  get_palette <- colorRampPalette(brewer.pal(12, "Set3"))
  cell_colors <- get_palette(n_colors_needed)
}
names(cell_colors) <- all_cell_types

# 确保 "Other" 是灰色
if("Other" %in% names(cell_colors)) {
  cell_colors["Other"] <- "#D3D3D3"
}

# 绘制
p_count_top10 <- ggplot(stack_data_top10, aes(x = cn_detailed, y = count, fill = detailed_grouped)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  scale_fill_manual(values = cell_colors, name = "Cell Type") +
  labs(title = "Cell Type Count by CN (Top 10 per CN + Others)",
       x = "Cell Neighborhood (CN)",
       y = "Cell Count") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "right",
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.4, "cm")
  )

ggsave("metLN_cn_celltype_count_top10_per_CN.pdf", p_count_top10, width = 16, height = 8, dpi = 300)

#TLS区的差异
#比较mPT、nmPT和metLN的TLS区
library(ggplot2)
library(dplyr)
library(patchwork)

# ===================== 提取三个数据集中 TLS_Zone 的细胞 =====================

# mPT中的TLS_Zone
mPT_tls <- mPT@meta.data %>% filter(cn_detailed == "TLS_Zone")
mPT_tls$dataset <- "mPT"

# nmPT中的TLS_Zone
nmPT_tls <- nmPT@meta.data %>% filter(cn_detailed == "TLS_Zone")
nmPT_tls$dataset <- "nmPT"

# metLN中的TLS_Zone（注意metLN中可能叫其他名字，根据实际注释调整）
# 根据之前的数据，metLN中CN8是B_T_Cell_Zone，CN1是B_Cell_Follicle_Zone，这些都接近TLS
# 这里使用 "B_T_Cell_Zone" 作为metLN中接近TLS的区域
metLN_tls <- metLN@meta.data %>% filter(cn_detailed == "B_T_Cell_Zone")
metLN_tls$dataset <- "metLN"

# 如果metLN有专门的TLS_Zone，用下面这行替换上面的
# metLN_tls <- metLN@meta.data %>% filter(cn_detailed == "TLS_Zone")

# 合并数据
combined_tls <- bind_rows(mPT_tls, nmPT_tls, metLN_tls)

cat("mPT TLS_Zone 细胞数：", nrow(mPT_tls), "\n")
cat("nmPT TLS_Zone 细胞数：", nrow(nmPT_tls), "\n")
cat("metLN B_T_Cell_Zone 细胞数：", nrow(metLN_tls), "\n")

# ===================== 颜色方案 =====================
cluster_colors <- c(
  "#A0CBE8", "#FFBE7D", "#8CD17D", "#86BCB6", "#FF9D9A",
  "#FABFD2", "#D4A6C8", "#F1CE63", "#D7B5A6", "#B07AA1",
  "#BAB0AC", "#4E79A7", "#F28E2B", "#59A14F", "#E15759",
  "#499894", "#D37295", "#B6992D", "#9D7660", "#79706E"
)

# ===================== 1. TLS_Zone 总细胞组成占比差异堆叠图 =====================

# 计算每个数据集中各细胞类型的占比
total_composition <- combined_tls %>%
  group_by(dataset, detailed) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(dataset) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ungroup()

# 获取top 15细胞类型
top_cell_types <- total_composition %>%
  group_by(detailed) %>%
  summarise(total_count = sum(count)) %>%
  arrange(desc(total_count)) %>%
  slice_head(n = 15) %>%
  pull(detailed)

# 其他细胞类型合并为 "Other"
total_composition_plot <- total_composition %>%
  mutate(cell_group = ifelse(detailed %in% top_cell_types, detailed, "Other")) %>%
  group_by(dataset, cell_group) %>%
  summarise(percentage = sum(percentage), .groups = "drop")

p_total <- ggplot(total_composition_plot, aes(x = dataset, y = percentage, fill = cell_group)) +
  geom_bar(stat = "identity", position = "stack", width = 0.6) +
  scale_fill_manual(values = cluster_colors[1:length(unique(total_composition_plot$cell_group))], 
                    name = "Cell Type") +
  labs(title = "TLS_Like_Zone: Cell Composition Comparison",
       subtitle = paste("mPT (n=", nrow(mPT_tls), ") | nmPT (n=", nrow(nmPT_tls), ") | metLN (n=", nrow(metLN_tls), ")"),
       x = "",
       y = "Percentage (%)") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 11),
    legend.position = "right",
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.4, "cm")
  )

ggsave("TLS_Zone_composition_comparison_three_groups.pdf", p_total, width = 14, height = 6, dpi = 300)
print(p_total)

# ===================== 2. T细胞亚群占比差异堆叠图 =====================

# 提取T细胞（以 T_ 开头，排除 Tumor）
t_cell_types <- grep("^T_", unique(combined_tls$detailed), value = TRUE)
t_cell_types <- t_cell_types[!grepl("Tumor", t_cell_types)]

t_composition <- combined_tls %>%
  filter(detailed %in% t_cell_types) %>%
  group_by(dataset, detailed) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(dataset) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ungroup()

if (nrow(t_composition) > 0) {
  p_t <- ggplot(t_composition, aes(x = dataset, y = percentage, fill = detailed)) +
    geom_bar(stat = "identity", position = "stack", width = 0.6) +
    scale_fill_manual(values = cluster_colors[1:length(unique(t_composition$detailed))], 
                      name = "T Cell Subtype") +
    labs(title = "TLS_Like_Zone: T Cell Subtype Comparison",
         x = "",
         y = "Percentage (%)") +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text.x = element_text(size = 12, face = "bold"),
      axis.text.y = element_text(size = 10),
      axis.title.y = element_text(size = 11),
      legend.position = "right",
      legend.text = element_text(size = 8),
      legend.key.size = unit(0.4, "cm")
    )
  
  ggsave("TLS_Zone_composition_comparison_T_cells_three_groups.pdf", p_t, width = 8, height = 5, dpi = 300)
  print(p_t)
} else {
  cat("没有找到T细胞\n")
}

# ===================== 3. B细胞亚群占比差异堆叠图 =====================

# 提取B细胞（以 B_ 开头）
b_cell_types <- grep("^B_", unique(combined_tls$detailed), value = TRUE)

b_composition <- combined_tls %>%
  filter(detailed %in% b_cell_types) %>%
  group_by(dataset, detailed) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(dataset) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ungroup()

if (nrow(b_composition) > 0) {
  p_b <- ggplot(b_composition, aes(x = dataset, y = percentage, fill = detailed)) +
    geom_bar(stat = "identity", position = "stack", width = 0.6) +
    scale_fill_manual(values = cluster_colors[1:length(unique(b_composition$detailed))], 
                      name = "B Cell Subtype") +
    labs(title = "TLS_Like_Zone: B Cell Subtype Comparison",
         x = "",
         y = "Percentage (%)") +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text.x = element_text(size = 12, face = "bold"),
      axis.text.y = element_text(size = 10),
      axis.title.y = element_text(size = 11),
      legend.position = "right",
      legend.text = element_text(size = 8),
      legend.key.size = unit(0.4, "cm")
    )
  
  ggsave("TLS_Zone_composition_comparison_B_cells_three_groups.pdf", p_b, width = 8, height = 5, dpi = 300)
  print(p_b)
} else {
  cat("没有找到B细胞\n")
}


#肿瘤分布于哪些CN
tumor_cell_types <- grep("Tumor_", unique(metLN@meta.data$detailed), value = TRUE)

if(length(tumor_cell_types) == 0) {
  tumor_cell_types <- grep("tumor|malignant|cancer", 
                           unique(metLN@meta.data$detailed), 
                           value = TRUE, ignore.case = TRUE)
}

print(paste("找到", length(tumor_cell_types), "种肿瘤亚型："))
print(tumor_cell_types)

# =====================
# 计算每种肿瘤亚型在各CN中的占比
# =====================

tumor_cn_composition <- metLN@meta.data %>%
  filter(detailed %in% tumor_cell_types) %>%
  group_by(detailed, cellcharter_cluster) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(detailed) %>%
  mutate(proportion = count / sum(count) * 100) %>%
  ungroup() %>%
  filter(proportion > 0)

# =====================
# 确保 CN 按数字顺序排序 (CN1 到 CNn)
# =====================

# 获取所有出现的CN编号
all_cn_numbers <- sort(unique(tumor_cn_composition$cellcharter_cluster))

# 创建 CN 标签（CN1, CN2, ...）
tumor_cn_composition$cn_label <- factor(
  paste0("CN", tumor_cn_composition$cellcharter_cluster),
  levels = paste0("CN", all_cn_numbers)  # 按数字顺序设置级别
)

# =====================
# 定义颜色（根据实际CN数量）
# =====================

cluster_colors <- c(
  "#A0CBE8", "#FFBE7D", "#8CD17D", "#86BCB6", "#FF9D9A",
  "#FABFD2", "#D4A6C8", "#F1CE63", "#D7B5A6", "#B07AA1",
  "#BAB0AC", "#4E79A7", "#F28E2B", "#59A14F", "#E15759",
  "#499894", "#D37295", "#B6992D", "#9D7660", "#79706E"
)

# 根据实际CN数量取颜色
n_cn <- length(all_cn_numbers)
cluster_colors_used <- cluster_colors[1:n_cn]
names(cluster_colors_used) <- paste0("CN", all_cn_numbers)

# 只保留实际出现的CN的颜色
cn_colors <- cluster_colors_used[paste0("CN", all_cn_numbers)]

# =====================
# 绘制分面饼图
# =====================

p_facet_pie <- ggplot(tumor_cn_composition, aes(x = "", y = proportion, fill = cn_label)) +
  geom_bar(stat = "identity", width = 1, color = "black", size = 0.5) +
  coord_polar("y", start = 0) +
  facet_wrap(~ detailed, ncol = 5) +
  scale_fill_manual(values = cn_colors, name = "CN", breaks = paste0("CN", all_cn_numbers)) +
  labs(title = "metLN: CN Composition for Each Tumor Subtype") +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.4, "cm")
  )

print(p_facet_pie)
ggsave("metLN_tumor_subtype_facet_pie.pdf", p_facet_pie, 
       width = 14, height = ceiling(length(tumor_cell_types)/3) * 4, 
       dpi = 300)

#比较mPT中肿瘤区域于metLN中的肿瘤相关区域
# 比较 mPT 的 Inflammatory_Vascular_Zone 和 metLN 中特定区域
library(ggplot2)
library(dplyr)

# ===================== 提取数据 =====================

# mPT中的 Inflammatory_Vascular_Zone
mPT_zone <- mPT@meta.data %>% filter(cn_detailed == "Inflammatory_Vascular_Zone")
mPT_zone$dataset <- "mPT_Inflam_Vasc"

# metLN中包含 "Tumor" 的CN（分开提取）
tumor_related_cn <- grep("Tumor", unique(metLN@meta.data$cn_detailed), value = TRUE)

# 分别提取每个Tumor相关的CN
metLN_tumor_list <- list()
for (cn in tumor_related_cn) {
  temp <- metLN@meta.data %>% filter(cn_detailed == cn)
  simple_name <- gsub("[-_]", "_", cn)
  simple_name <- gsub(" ", "_", simple_name)
  temp$dataset <- paste0("metLN_", simple_name)
  metLN_tumor_list[[cn]] <- temp
}

# metLN中的 Mast_Cell_Immune_Zone
metLN_mast <- metLN@meta.data %>% filter(cn_detailed == "Mast_Cell_Immune_Zone")
metLN_mast$dataset <- "metLN_Mast_Cell"

# 合并数据
combined_data <- bind_rows(mPT_zone, metLN_mast)
for (temp in metLN_tumor_list) {
  combined_data <- bind_rows(combined_data, temp)
}

# ===================== 删除指定的dataset =====================
combined_data <- combined_data %>%
  filter(dataset != "metLN_B_Cell_Immune_inflamed_Tumor_Interface")

# 验证
cat("删除后剩余的dataset：\n")
print(unique(combined_data$dataset))

# ===================== 颜色方案 =====================
cluster_colors <- c(
  "#A0CBE8", "#FFBE7D", "#8CD17D", "#86BCB6", "#FF9D9A",
  "#FABFD2", "#D4A6C8", "#F1CE63", "#D7B5A6", "#B07AA1",
  "#BAB0AC", "#4E79A7", "#F28E2B", "#59A14F", "#E15759",
  "#499894", "#D37295", "#B6992D", "#9D7660", "#79706E", "#D3D3D3"
)

# ===================== 总细胞组成占比差异堆叠图 =====================

total_composition <- combined_data %>%
  group_by(dataset, detailed) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(dataset) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ungroup()

# 获取top 15细胞类型
top_cell_types <- total_composition %>%
  group_by(detailed) %>%
  summarise(total_count = sum(count)) %>%
  arrange(desc(total_count)) %>%
  slice_head(n = 15) %>%
  pull(detailed)

# 将非top细胞类型标记为 "Other"
total_composition_plot <- total_composition %>%
  mutate(cell_group = ifelse(detailed %in% top_cell_types, detailed, "Other")) %>%
  group_by(dataset, cell_group) %>%
  summarise(percentage = sum(percentage), .groups = "drop") %>%
  mutate(cell_group = factor(cell_group, 
                              levels = c(sort(setdiff(unique(cell_group), "Other")), "Other")))

# 设置dataset因子水平
dataset_order <- unique(total_composition_plot$dataset)
total_composition_plot$dataset <- factor(total_composition_plot$dataset, levels = dataset_order)

# 计算颜色
n_colors <- length(unique(total_composition_plot$cell_group))
used_colors <- cluster_colors[1:n_colors]
names(used_colors) <- levels(total_composition_plot$cell_group)

p_total <- ggplot(total_composition_plot, aes(x = dataset, y = percentage, fill = cell_group)) +
  geom_bar(stat = "identity", position = "stack", width = 0.6) +
  scale_fill_manual(values = used_colors, name = "Cell Type") +
  labs(title = "Zone Comparison: mPT vs metLN",
       x = "",
       y = "Percentage (%)") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9, face = "bold"),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 11),
    legend.position = "right",
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.4, "cm")
  )

ggsave("Zone_comparison_mPT_metLN_total.pdf", p_total, width = 9, height = 6, dpi = 300)
print(p_total)

#CN相关性
# 比较 metLN 中四个 CN 与 mPT 炎症血管区细胞组成的 Pearson 相关性
library(ggplot2)
library(dplyr)
library(pheatmap)

# ===================== 提取各区域细胞组成 =====================

# mPT 炎症血管区
mPT_composition <- mPT@meta.data %>%
  filter(cn_detailed == "Inflammatory_Vascular_Zone") %>%
  group_by(detailed) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(percentage = count / sum(count) * 100,
         region = "mPT_Inflammatory_Vascular")

# metLN 四个CN（使用正确名称）
metLN_cn_list <- c("Immune-inflamed_Tumor_Core_A", 
                   "Immune-inflamed_Tumor_Core_B",
                   "Immune-inflamed_Tumor_Macrophage_Interface",
                   "Mast_Cell_Immune_Zone")

# 提取 metLN 各 CN 数据
metLN_compositions <- list()
for (cn in metLN_cn_list) {
  temp <- metLN@meta.data %>%
    filter(cn_detailed == cn) %>%
    group_by(detailed) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(percentage = count / sum(count) * 100,
           region = cn)
  metLN_compositions[[cn]] <- temp
}

# 合并所有数据
all_compositions <- bind_rows(mPT_composition)
for (temp in metLN_compositions) {
  all_compositions <- bind_rows(all_compositions, temp)
}

# 转换为宽格式
wide_composition <- all_compositions %>%
  select(region, detailed, percentage) %>%
  tidyr::pivot_wider(names_from = region, values_from = percentage, values_fill = 0)

wide_matrix <- as.matrix(wide_composition[, -1])
rownames(wide_matrix) <- wide_composition$detailed

# ===================== 计算 Pearson 相关性 =====================

cor_matrix_pearson <- cor(wide_matrix, method = "pearson")

print("=== Pearson 相关性矩阵 ===")
print(round(cor_matrix_pearson, 3))

# ===================== 绘制 Pearson 相关性热图 =====================

pdf("zone_composition_pearson_correlation_heatmap.pdf", width = 8, height = 8)
pheatmap(cor_matrix_pearson,
         main = "Pearson Correlation of Cell Composition",
         display_numbers = TRUE,
         number_format = "%.2f",
         color = colorRampPalette(c("#4E79A7", "white", "#E15759"))(100),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         fontsize = 10)
dev.off()

# ===================== mPT 炎症血管区 vs metLN Mast_Cell_Immune_Zone 相关性分析 =====================
library(ggplot2)
library(dplyr)

# 提取 mPT 炎症血管区的细胞组成
mPT_composition <- mPT@meta.data %>%
  filter(cn_detailed == "Inflammatory_Vascular_Zone") %>%
  group_by(detailed) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(percentage = count / sum(count) * 100,
         region = "mPT_Inflammatory_Vascular")

# 提取 metLN Mast_Cell_Immune_Zone 的细胞组成
metLN_composition <- metLN@meta.data %>%
  filter(cn_detailed == "Mast_Cell_Immune_Zone") %>%
  group_by(detailed) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(percentage = count / sum(count) * 100,
         region = "metLN_Mast_Cell_Immune_Zone")

# 合并两个区域的数据
comparison_data <- bind_rows(mPT_composition, metLN_composition)

# 转换为宽格式（每个细胞类型一行，两个区域各一列）
wide_data <- comparison_data %>%
  select(region, detailed, percentage) %>%
  tidyr::pivot_wider(names_from = region, values_from = percentage, values_fill = 0)

# 计算 Pearson 相关系数
cor_pearson <- cor(wide_data$mPT_Inflammatory_Vascular, 
                   wide_data$metLN_Mast_Cell_Immune_Zone, 
                   method = "pearson")

# 计算 Spearman 相关系数
cor_spearman <- cor(wide_data$mPT_Inflammatory_Vascular, 
                    wide_data$metLN_Mast_Cell_Immune_Zone, 
                    method = "spearman")

cat("=== 相关性分析结果 ===\n")
cat("Pearson 相关系数: ", round(cor_pearson, 4), "\n")
cat("Spearman 相关系数: ", round(cor_spearman, 4), "\n")

# ===================== 散点图 =====================

# 计算线性回归
lm_fit <- lm(metLN_Mast_Cell_Immune_Zone ~ mPT_Inflammatory_Vascular, data = wide_data)

# 获取R²和p值
r_squared <- summary(lm_fit)$r.squared
p_value <- summary(lm_fit)$coefficients[2, 4]

# 散点图
p_scatter <- ggplot(wide_data, aes(x = mPT_Inflammatory_Vascular, y = metLN_Mast_Cell_Immune_Zone)) +
  geom_point(size = 3, alpha = 0.7, color = "#4E79A7") +
  geom_smooth(method = "lm", se = TRUE, color = "#E15759", fill = "#FABFD2") +
  annotate("text", x = max(wide_data$mPT_Inflammatory_Vascular) * 0.7, 
           y = max(wide_data$metLN_Mast_Cell_Immune_Zone) * 0.95,
           label = paste("Pearson r =", round(cor_pearson, 3),
                        #"\nSpearman ρ =", round(cor_spearman, 3),
                        "\nR² =", round(r_squared, 3)),
           size = 4, hjust = 0) +
  labs(title = "mPT Inflammatory_Vascular_Zone vs metLN Mast_Cell_Immune_Zone",
       x = "mPT Inflammatory_Vascular Zone (%)",
       y = "metLN Mast_Cell_Immune_Zone (%)") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10)
  )

ggsave("mPT_vs_metLN_Mast_Cell_correlation.pdf", p_scatter, width = 6, height = 5, dpi = 300)
print(p_scatter)

#点图
selected_fovs <- c(239, 240)

d2_data <- metLN@meta.data %>% 
  filter(sample == "D2", fov %in% selected_fovs)

# 筛选X轴范围（如果需要，可以调整或删除）
# d2_data <- d2_data %>%
#   filter(CenterX_global_px >= 73060 & CenterX_global_px <= 85024)

cat("筛选后细胞数：", nrow(d2_data), "\n")
cat("X轴范围：", range(d2_data$CenterX_global_px), "\n")
cat("Y轴范围：", range(d2_data$CenterY_global_px), "\n")

# ===================== 获取实际出现的CN =====================
actual_cn <- sort(unique(d2_data$cn_detailed))
cat("实际出现的CN：\n")
print(actual_cn)

# ===================== 定义颜色 =====================
cluster_colors <- c(
  "#A0CBE8", "#FFBE7D", "#8CD17D", "#86BCB6", "#FF9D9A",
  "#FABFD2", "#D4A6C8", "#F1CE63", "#D7B5A6", "#B07AA1",
  "#BAB0AC", "#4E79A7", "#F28E2B", "#59A14F", "#E15759",
  "#499894", "#D37295", "#B6992D", "#9D7660", "#79706E"
)

# 将颜色直接赋给CN注释名称
used_colors <- cluster_colors[1:length(actual_cn)]
names(used_colors) <- actual_cn

# ===================== 绘制点图（所有FOV在一张图上） =====================

p <- ggplot(d2_data, aes(x = CenterX_global_px, y = CenterY_global_px, color = cn_detailed)) +
  geom_point(size = 0.5, alpha = 0.8) +
  scale_color_manual(values = used_colors, name = "Cell Neighborhood (CN)") +
  coord_fixed(ratio = 1) +
  labs(
    title = "metLN D2 Sample: Spatial Distribution of Cell Neighborhoods",
    subtitle = paste("FOV:", paste(selected_fovs, collapse = ", "), 
                     "| Total cells:", nrow(d2_data)),
    x = "X Coordinate (px)",
    y = "Y Coordinate (px)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    plot.subtitle = element_text(hjust = 0.5, size = 8),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.6, "cm"),
    aspect.ratio = 0.5
  ) +
  guides(color = guide_legend(override.aes = list(size = 3)))

ggsave("metLN_D2_FOV239_240_CN_spatial.pdf", p, width = 8, height = 4, dpi = 300)

```
# negLN subtype CN
```R
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
negLN <- readRDS("negLN_cellcharter_FINAL.rds")

#CN组成
# 查看每个CN的细胞组成
cn_composition <- table(negLN@meta.data$cellcharter_cluster, 
                        negLN@meta.data$detailed)

# 计算每个CN的比例
cn_prop <- prop.table(cn_composition, margin = 1) * 100

# 打印每个CN的主要细胞类型（前5种）
for (cn in rownames(cn_prop)) {
  cat("\n=== CN", cn, "===\n")
  top10 <- sort(cn_prop[cn, ], decreasing = TRUE)[1:10]
  for (i in 1:length(top10)) {
    cat(names(top10)[i], ": ", round(top10[i], 1), "%\n", sep = "")
  }
}

#CN注释
cn_annotation_detailed <- c(
  "1"  = "T_Cell_Dominant_Zone",
  "2"  = "Mast_Cell_Immune_Zone",
  "3"  = "FRC-rich_Lymphoid_Zone",
  "4"  = "Immune_Zone",
  "5"  = "Macrophage_Dominant_Zone_F",
  "6"  = "B_Cell_Follicle_Zone",
  "7"  = "Plasma_Enriched_Zone_H",
  "8"  = "Macrophage_Dominant_Zone_D",
  "9"  = "B_Germinal_Center_Zone",
  "10" = "Plasma_Enriched_Zone_I",
  "11" = "Vascular_Stromal_Zone",
  "12" = "Plasma_Enriched_Zone_G"
)
# 应用到 Seurat 对象
negLN@meta.data$cn_detailed <- cn_annotation_detailed[as.character(negLN@meta.data$cellcharter_cluster)]

library(ggplot2)
library(dplyr)
library(patchwork)

# ===================== 1. 准备数据 =====================

# 样本顺序（按首字母排序）
samples_order <- c("A2", "A4", "B1", "B3", "B5", "C3", "D1", "D4")
non_metastatic_samples <- samples_order[1:4]  # "A2", "A4", "B1", "B3"
metastatic_samples <- samples_order[5:8]      # "B5", "C3", "D1", "D4"

# 计算总体比例
total_prop <- negLN@meta.data %>%
  group_by(cn_detailed) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(proportion = count / sum(count) * 100,
         group = "Overall")

# 计算未转移样本组比例
non_metastatic_prop <- negLN@meta.data %>%
  filter(sample %in% non_metastatic_samples) %>%
  group_by(cn_detailed) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(proportion = count / sum(count) * 100,
         group = "Non-Metastatic")

# 计算转移样本组比例
metastatic_prop <- negLN@meta.data %>%
  filter(sample %in% metastatic_samples) %>%
  group_by(cn_detailed) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(proportion = count / sum(count) * 100,
         group = "Metastatic")

# 合并分组数据
group_data <- bind_rows(total_prop, non_metastatic_prop, metastatic_prop)

# 计算各样本比例
sample_prop <- negLN@meta.data %>%
  group_by(sample, cn_detailed) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(sample) %>%
  mutate(proportion = count / sum(count) * 100) %>%
  ungroup()

# 设置因子水平
group_data$cn_detailed <- factor(group_data$cn_detailed, 
                                  levels = sort(unique(group_data$cn_detailed)))
group_data$group <- factor(group_data$group, 
                            levels = c("Overall", "Non-Metastatic", "Metastatic"))

sample_prop$cn_detailed <- factor(sample_prop$cn_detailed, 
                                   levels = sort(unique(sample_prop$cn_detailed)))
sample_prop$sample <- factor(sample_prop$sample, levels = samples_order)

# ===================== 2. 颜色方案 =====================

cluster_colors <- c(
  "#A0CBE8", "#FFBE7D", "#8CD17D", "#86BCB6", "#FF9D9A",
  "#FABFD2", "#D4A6C8", "#F1CE63", "#D7B5A6", "#B07AA1",
  "#BAB0AC", "#4E79A7", "#F28E2B", "#59A14F", "#E15759",
  "#499894", "#D37295", "#B6992D", "#9D7660", "#79706E"
)

# ===================== 3. 图1：分组统计堆叠图 =====================

p_group <- ggplot(group_data, aes(x = group, y = proportion, fill = cn_detailed)) +
  geom_bar(stat = "identity", position = "fill", width = 0.6) +
  scale_fill_manual(values = cluster_colors, name = "CN") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "A: CN Composition by Group",
       x = "",
       y = "Proportion") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10, face = "bold"),
    axis.text.y = element_text(size = 9),
    axis.title.y = element_text(size = 10),
    legend.position = "none",
    panel.grid = element_blank()
  )

# ===================== 4. 图2：各样本堆叠图 =====================

p_samples <- ggplot(sample_prop, aes(x = sample, y = proportion, fill = cn_detailed)) +
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  scale_fill_manual(values = cluster_colors, name = "CN") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "B: CN Composition by Individual Sample",
       x = "",
       y = "Proportion") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    axis.title.y = element_text(size = 10),
    legend.position = "right",
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.4, "cm"),
    panel.grid = element_blank()
  )

# ===================== 5. 提取图例 =====================

library(cowplot)
legend <- get_legend(p_samples)

# 移除图2的图例
p_samples_no_legend <- p_samples + theme(legend.position = "none")

# ===================== 6. 合并两张图 =====================

combined_plot <- (p_group | p_samples_no_legend) +
  plot_annotation(
    title = "negLN: Cell Neighborhood Composition",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
  ) +
  plot_layout(widths = c(0.6, 1.4))

# 添加图例到右侧
final_plot <- plot_grid(combined_plot, legend, ncol = 2, rel_widths = c(3, 0.4))

ggsave("negLN_cn_proportion_combined.pdf", final_plot, width = 16, height = 6, dpi = 300)


# 3. 每个邻域细胞组成热图
# =====================
message("=== 3. 生成邻域细胞组成热图 ===")

cluster_composition <- table(negLN@meta.data$cellcharter_cluster, 
                           negLN@meta.data$detailed)

cluster_composition_df <- as.data.frame(cluster_composition)
colnames(cluster_composition_df) <- c("Cluster", "CellType", "Count")

cluster_totals <- aggregate(Count ~ Cluster, data = cluster_composition_df, sum)
cluster_composition_df <- merge(cluster_composition_df, cluster_totals, by = "Cluster")
cluster_composition_df$Proportion <- cluster_composition_df$Count.x / cluster_composition_df$Count.y

# 转换为宽格式
heatmap_data <- cluster_composition_df %>%
  select(Cluster, CellType, Proportion) %>%
  pivot_wider(names_from = CellType, values_from = Proportion, values_fill = 0)

# 按Cluster数字排序
cluster_numeric <- as.numeric(as.character(heatmap_data$Cluster))
heatmap_data <- heatmap_data[order(cluster_numeric), ]

# 创建行标签（使用你的注释名称）
row_labels_annotated <- cn_annotation_detailed[as.character(cluster_numeric[order(cluster_numeric)])]

# 设置行名
rownames(heatmap_data) <- row_labels_annotated
heatmap_matrix <- as.matrix(heatmap_data[, -1])

# 绘制热图（字体增大的版本）
p_heatmap <- pheatmap(heatmap_matrix,
                    cluster_rows = TRUE,
                    cluster_cols = TRUE,
                    display_numbers = FALSE,
                    main = "Cell Type Composition by Spatial Domain",
                    fontsize = 14,           # 全局字体增大
                    fontsize_row = 10,       # 行标签字体（注释名称较长，用10合适）
                    fontsize_col = 9,        # 列标签字体
                    labels_row = row_labels_annotated,
                    color = colorRampPalette(c("#f0f0f0", "#3182bd", "#08519c"))(100),
                    cellwidth = 11,          # 单元格宽度
                    cellheight = 11)         # 单元格高度

# 保存
pdf("negLN_CN_cell_type_composition_heatmap_annotated.pdf", width = 16, height = 12)
print(p_heatmap)
dev.off()


#空间点图
selected_fovs <- c(380, 382)

b1_data <- negLN@meta.data %>% 
  filter(sample == "B1", fov %in% selected_fovs)

# 筛选X轴范围（如果需要，可以调整或删除）
# b1_data <- b1_data %>%
#   filter(CenterX_global_px >= 73060 & CenterX_global_px <= 85024)

cat("筛选后细胞数：", nrow(b1_data), "\n")
cat("X轴范围：", range(b1_data$CenterX_global_px), "\n")
cat("Y轴范围：", range(b1_data$CenterY_global_px), "\n")

# ===================== 获取实际出现的CN =====================
actual_cn <- sort(unique(b1_data$cn_detailed))
cat("实际出现的CN：\n")
print(actual_cn)

# ===================== 定义颜色 =====================
cluster_colors <- c(
  "#A0CBE8", "#FFBE7D", "#8CD17D", "#86BCB6", "#FF9D9A",
  "#FABFD2", "#D4A6C8", "#F1CE63", "#D7B5A6", "#B07AA1",
  "#BAB0AC", "#4E79A7", "#F28E2B", "#59A14F", "#E15759",
  "#499894", "#D37295", "#B6992D", "#9D7660", "#79706E"
)

# 将颜色直接赋给CN注释名称
used_colors <- cluster_colors[1:length(actual_cn)]
names(used_colors) <- actual_cn

# ===================== 绘制点图（所有FOV在一张图上） =====================

p <- ggplot(b1_data, aes(x = CenterX_global_px, y = CenterY_global_px, color = cn_detailed)) +
  geom_point(size = 0.5, alpha = 0.8) +
  scale_color_manual(values = used_colors, name = "Cell Neighborhood (CN)") +
  coord_fixed(ratio = 1) +
  labs(
    title = "negLN B1 Sample: Spatial Distribution of Cell Neighborhoods",
    subtitle = paste("FOV:", paste(selected_fovs, collapse = ", "), 
                     "| Total cells:", nrow(b1_data)),
    x = "X Coordinate (px)",
    y = "Y Coordinate (px)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    plot.subtitle = element_text(hjust = 0.5, size = 8),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.6, "cm"),
    aspect.ratio = 0.5
  ) +
  guides(color = guide_legend(override.aes = list(size = 3)))

ggsave("negLN_B1_FOV380_382_CN_spatial.pdf", p, width = 8, height = 4, dpi = 300)


selected_fovs <- c(193, 195, 198, 199)

b3_data <- negLN@meta.data %>% 
  filter(sample == "B3", fov %in% selected_fovs)

# 筛选X轴范围（如果需要，可以调整或删除）
# b3_data <- b3_data %>%
#   filter(CenterX_global_px >= 73060 & CenterX_global_px <= 85024)

cat("筛选后细胞数：", nrow(b3_data), "\n")
cat("X轴范围：", range(b3_data$CenterX_global_px), "\n")
cat("Y轴范围：", range(b3_data$CenterY_global_px), "\n")

# ===================== 获取实际出现的CN =====================
actual_cn <- sort(unique(b3_data$cn_detailed))
cat("实际出现的CN：\n")
print(actual_cn)

# ===================== 定义颜色 =====================
cluster_colors <- c(
  "#A0CBE8", "#FFBE7D", "#8CD17D", "#86BCB6", "#FF9D9A",
  "#FABFD2", "#D4A6C8", "#F1CE63", "#D7B5A6", "#B07AA1",
  "#BAB0AC", "#4E79A7", "#F28E2B", "#59A14F", "#E15759",
  "#499894", "#D37295", "#B6992D", "#9D7660", "#79706E"
)

# 将颜色直接赋给CN注释名称
used_colors <- cluster_colors[1:length(actual_cn)]
names(used_colors) <- actual_cn

# ===================== 绘制点图（所有FOV在一张图上） =====================

p <- ggplot(b3_data, aes(x = CenterX_global_px, y = CenterY_global_px, color = cn_detailed)) +
  geom_point(size = 0.05, alpha = 0.8) +
  scale_color_manual(values = used_colors, name = "Cell Neighborhood (CN)") +
  coord_fixed(ratio = 1) +
  labs(
    title = "negLN B3 Sample: Spatial Distribution of Cell Neighborhoods",
    subtitle = paste("FOV:", paste(selected_fovs, collapse = ", "), 
                     "| Total cells:", nrow(b3_data)),
    x = "X Coordinate (px)",
    y = "Y Coordinate (px)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
    plot.subtitle = element_text(hjust = 0.5, size = 8),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.6, "cm"),
    aspect.ratio = 1
  ) +
  guides(color = guide_legend(override.aes = list(size = 3)))

ggsave("negLN_B3_FOV193_195_198_199_CN_spatial.pdf", p, width = 6, height = 4, dpi = 300)

#C3
selected_fovs <- c(167, 168)

c3_data <- negLN@meta.data %>% 
  filter(sample == "C3", fov %in% selected_fovs)

# 筛选X轴范围（如果需要，可以调整或删除）
# c3_data <- c3_data %>%
#   filter(CenterX_global_px >= 73060 & CenterX_global_px <= 85024)

cat("筛选后细胞数：", nrow(c3_data), "\n")
cat("X轴范围：", range(c3_data$CenterX_global_px), "\n")
cat("Y轴范围：", range(c3_data$CenterY_global_px), "\n")

# ===================== 获取实际出现的CN =====================
actual_cn <- sort(unique(c3_data$cn_detailed))
cat("实际出现的CN：\n")
print(actual_cn)

# ===================== 定义颜色 =====================
cluster_colors <- c(
  "#A0CBE8", "#FFBE7D", "#8CD17D", "#86BCB6", "#FF9D9A",
  "#FABFD2", "#D4A6C8", "#F1CE63", "#D7B5A6", "#B07AA1",
  "#BAB0AC", "#4E79A7", "#F28E2B", "#59A14F", "#E15759",
  "#499894", "#D37295", "#B6992D", "#9D7660", "#79706E"
)

# 将颜色直接赋给CN注释名称
used_colors <- cluster_colors[1:length(actual_cn)]
names(used_colors) <- actual_cn

# ===================== 绘制点图（所有FOV在一张图上） =====================

p <- ggplot(c3_data, aes(x = CenterX_global_px, y = CenterY_global_px, color = cn_detailed)) +
  geom_point(size = 0.5, alpha = 0.8) +
  scale_color_manual(values = used_colors, name = "Cell Neighborhood (CN)") +
  coord_fixed(ratio = 1) +
  labs(
    title = "negLN C3 Sample: Spatial Distribution of Cell Neighborhoods",
    subtitle = paste("FOV:", paste(selected_fovs, collapse = ", "), 
                     "| Total cells:", nrow(c3_data)),
    x = "X Coordinate (px)",
    y = "Y Coordinate (px)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    plot.subtitle = element_text(hjust = 0.5, size = 8),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.6, "cm"),
    aspect.ratio = 0.5
  ) +
  guides(color = guide_legend(override.aes = list(size = 3)))

ggsave("negLN_C3_FOV167_168_CN_spatial.pdf", p, width = 8, height = 4, dpi = 300)

selected_fovs <- c(315, 316)

d1_data <- negLN@meta.data %>% 
  filter(sample == "D1", fov %in% selected_fovs)

# 筛选X轴范围（如果需要，可以调整或删除）
# d1_data <- d1_data %>%
#   filter(CenterX_global_px >= 73060 & CenterX_global_px <= 85024)

cat("筛选后细胞数：", nrow(d1_data), "\n")
cat("X轴范围：", range(d1_data$CenterX_global_px), "\n")
cat("Y轴范围：", range(d1_data$CenterY_global_px), "\n")

# ===================== 获取实际出现的CN =====================
actual_cn <- sort(unique(d1_data$cn_detailed))
cat("实际出现的CN：\n")
print(actual_cn)

# ===================== 定义颜色 =====================
cluster_colors <- c(
  "#A0CBE8", "#FFBE7D", "#8CD17D", "#86BCB6", "#FF9D9A",
  "#FABFD2", "#D4A6C8", "#F1CE63", "#D7B5A6", "#B07AA1",
  "#BAB0AC", "#4E79A7", "#F28E2B", "#59A14F", "#E15759",
  "#499894", "#D37295", "#B6992D", "#9D7660", "#79706E"
)

# 将颜色直接赋给CN注释名称
used_colors <- cluster_colors[1:length(actual_cn)]
names(used_colors) <- actual_cn

# ===================== 绘制点图（所有FOV在一张图上） =====================

p <- ggplot(d1_data, aes(x = CenterX_global_px, y = CenterY_global_px, color = cn_detailed)) +
  geom_point(size = 0.5, alpha = 0.8) +
  scale_color_manual(values = used_colors, name = "Cell Neighborhood (CN)") +
  coord_fixed(ratio = 1) +
  labs(
    title = "negLN D1 Sample: Spatial Distribution of Cell Neighborhoods",
    subtitle = paste("FOV:", paste(selected_fovs, collapse = ", "), 
                     "| Total cells:", nrow(d1_data)),
    x = "X Coordinate (px)",
    y = "Y Coordinate (px)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    plot.subtitle = element_text(hjust = 0.5, size = 8),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.6, "cm"),
    aspect.ratio = 0.5
  ) +
  guides(color = guide_legend(override.aes = list(size = 3)))

ggsave("negLN_D1_FOV315_316_CN_spatial.pdf", p, width = 8, height = 4, dpi = 300)

```
# fib fov
```R

fib <- readRDS("fib_anno_new.rds")
obj <- readRDS("YA2025263-1_fin.rds")

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
                       "Basal cells" = "Malignant cells"
)

# 检查
table(obj$CellType)


obj@meta.data$sub_cell_type <- NA
common_cells <- intersect(rownames(fib@meta.data), rownames(obj@meta.data))
obj@meta.data[common_cells, "sub_cell_type"] <- fib@meta.data[common_cells, "sub_cell_type"]

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

library(ggplot2)
library(dplyr)

# ==================== 1. 定义要显示的细胞类型 ====================
# 成纤维亚型
target_fib_subtypes <- c("Fib_apCAF_1", "Fib_apCAF_2", "Fib_iCAF_1", 
                         "Fib_iCAF_2", "Fib_matCAF", "Fib_myCAF")

# 恶性细胞（单独作为一类）
malignant_cell <- "Malignant cells"

# 合并所有要显示的细胞类型
target_cell_types <- c(target_fib_subtypes, malignant_cell)

# ==================== 2. 筛选数据 ====================
selected_fovs <- c(272,273,274,275,276,277,278,279,280,281,282)

fov_data <- obj@meta.data %>%
  filter(sample == "B2", fov %in% selected_fovs)

# 创建分组变量：目标细胞类型保留具体名称，其他细胞归为 "Other"
fov_data$cell_group <- "Other"
for (ct in target_cell_types) {
  fov_data$cell_group[fov_data$detailed == ct] <- ct
}

# ==================== 3. 颜色方案 ====================
fib_colors <- c(
  "Fib_apCAF_1" = "#793c1b",
  "Fib_apCAF_2" = "#B5EFB5",
  "Fib_iCAF_1"  = "#ffc089",
  "Fib_iCAF_2"  = "#ffff33",
  "Fib_matCAF"  = "#6a3d9a",
  "Fib_myCAF"   = "#1f78b4"
)

malignant_color <- c("Malignant cells" = "#df928e")
other_color <- c("Other" = "#E0E0E0")

# 合并颜色
cell_colors <- c(fib_colors, malignant_color, other_color)

# 只保留出现在数据中的颜色
used_colors <- cell_colors[names(cell_colors) %in% unique(fov_data$cell_group)]
if ("Other" %in% unique(fov_data$cell_group) && !"Other" %in% names(used_colors)) {
  used_colors <- c(used_colors, "Other" = "#E0E0E0")
}

# ==================== 4. 设置绘图顺序 ====================
plot_order <- c("Other", target_fib_subtypes, malignant_cell)
fov_data$cell_group <- factor(fov_data$cell_group, levels = plot_order)

# ==================== 5. 绘图（限制 X 轴范围）====================
p <- ggplot(fov_data, aes(x = CenterX_global_px, y = CenterY_global_px, color = cell_group)) +
  geom_point(size = 0.01, alpha = 0.9) +
  scale_color_manual(values = used_colors, name = "Cell Type") +
  # 关键修改：限制 X 轴显示范围
  xlim(61565, 82820) +  # 只显示这个 X 范围内的点
  labs(title = "nmPT (B2) - FOVs 272-282: Fibroblasts Subtypes and Malignant Cells",
       x = "X Coordinate (px)",
       y = "Y Coordinate (px)") +
  theme_bw() +
  theme(
    aspect.ratio = 1/2,
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    legend.position = "right",
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    legend.key.size = unit(0.6, "cm"),
    legend.text = element_text(size = 8)
  ) +
  guides(color = guide_legend(override.aes = list(size = 3), ncol = 1))

# ==================== 6. 保存 ====================
ggsave("nmPT_B2_fov272_282_fib_malignant_xlim.pdf", p, width = 7, height = 3.5, dpi = 300)
```
# mPT CN5 chat
```R
# ===================== 
# 分析 CN5 内部的细胞互作
# =====================

library(ggplot2)
library(dplyr)
library(pheatmap)
library(doParallel)
library(foreach)
library(deldir)

mPT <- readRDS("mPT_cellcharter_FINAL.rds")
cn_composition <- table(mPT@meta.data$cellcharter_cluster, 
                        mPT@meta.data$detailed)

# 计算每个CN的比例
cn_prop <- prop.table(cn_composition, margin = 1) * 100

# 打印每个CN的主要细胞类型（前5种）
for (cn in rownames(cn_prop)) {
  cat("\n=== CN", cn, "===\n")
  top10 <- sort(cn_prop[cn, ], decreasing = TRUE)[1:10]
  for (i in 1:length(top10)) {
    cat(names(top10)[i], ": ", round(top10[i], 1), "%\n", sep = "")
  }
}

#CN注释

cn_annotation_detailed <- c(
  "1"  = "Alveolar_Homeostatic_Zone",
  "2"  = "Immune_Zone",
  "3"  = "Stromal_Remodeling_Zone",
  "4"  = "TLS_Zone",
  "5"  = "Inflammatory_Vascular_Zone",
  "6"  = "Macrophage_Dominant_Zone_A",
  "7"  = "Vascular_Stromal_Zone",
  "8"  = "Pericyte_Dominant_Zone",
  "9"  = "Plasma_Enriched_Zone_A",
  "10" = "Plasma_Enriched_Zone_B",
  "11" = "Metabolic_Tumor_Core_A",
  "12" = "Metabolic_Tumor_Mast_Interface",
  "13" = "Metabolic_Tumor_Core_B",
  "14" = "Metabolic_Tumor_Stromal_Interface",
  "15" = "Plasma_Enriched_Zone_C"
)
# 应用到 Seurat 对象
mPT@meta.data$cn_detailed <- cn_annotation_detailed[as.character(mPT@meta.data$cellcharter_cluster)]

# ===================== 定义空间分析函数 =====================
run_spatial_analysis_cn5 <- function(coords_sub, labels_sub, sample_name, nperm = 1000) {
  
  cat("\n========================================\n")
  cat("开始分析样本：", sample_name, "\n")
  cat("CN5内细胞数量：", length(labels_sub), "\n")
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

# ===================== 提取 CN5 内的细胞 =====================
samples <- unique(mPT$sample)
cat("mPT 样本:", paste(samples, collapse = ", "), "\n")

# 存储结果
cn5_results <- list()

for (sample_id in samples) {
  cat("\n\n########## 处理样本：", sample_id, " ##########\n")
  
  # 提取该样本中 CN5 的细胞
  sample_cells <- subset(mPT, subset = sample == sample_id & cn_detailed == "Inflammatory_Vascular_Zone")
  
  if (nrow(sample_cells@meta.data) == 0) {
    cat("警告：", sample_id, "CN5中没有细胞，跳过\n")
    next
  }
  
  cat("CN5细胞数量：", nrow(sample_cells@meta.data), "\n")
  
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
  
  # 过滤稀有细胞类型（数量 < 5 的过滤，因为CN5内细胞数较少）
  cell_counts <- table(labels)
  rare_threshold <- 5
  keep_types <- names(cell_counts[cell_counts >= rare_threshold])
  keep_cells <- labels %in% keep_types
  coords <- coords[keep_cells, , drop = FALSE]
  labels <- labels[keep_cells]
  
  cat("过滤稀有细胞后数量：", length(labels), "\n")
  
  if (length(labels) < 30) {
    cat("警告：", sample_id, "过滤后细胞太少，跳过\n")
    next
  }
  
  # 确保坐标是数值型
  coords <- as.matrix(coords)
  colnames(coords) <- c("x", "y")
  
  # 运行分析
  result <- run_spatial_analysis_cn5(coords, labels, sample_id, nperm = 500)
  
  if (!is.null(result)) {
    cn5_results[[sample_id]] <- result
    
    # 绘制并保存热图
    if (nrow(result$mat_log2fc) > 1 && ncol(result$mat_log2fc) > 1) {
      mat_plot <- result$mat_log2fc
      
      # 处理 Inf
      mat_plot[is.infinite(mat_plot) & mat_plot < 0] <- -10
      mat_plot[is.infinite(mat_plot) & mat_plot > 0] <- 10
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
      pdf(paste0("mPT_CN5_", sample_id, "_contact_log2fc.pdf"), width = 12, height = 10)
      pheatmap(mat_plot,
               cluster_rows = TRUE,
               cluster_cols = TRUE,
               main = paste("CN5 - Log2(FC) -", sample_id),
               fontsize = 10,
               display_numbers = signif_symbols,
               number_color = "black",
               fontsize_number = 8)
      dev.off()
      
      cat("保存热图：mPT_CN5_", sample_id, "_contact_log2fc.pdf\n", sep = "")
    }
  }
}

# ===================== 保存结果 =====================
saveRDS(cn5_results, file = "mPT_CN5_cell_interaction_results.rds")
cat("\n\n========== 分析完成 ==========\n")
cat("成功分析样本数：", length(cn5_results), "\n")

# ===================== 整合多样本CN5结果 =====================
all_cell_types <- unique(unlist(lapply(cn5_results, function(x) x$cell_types)))
cat("总细胞类型数量：", length(all_cell_types), "\n")
cat("所有细胞类型：", paste(all_cell_types, collapse = ", "), "\n")

if (length(all_cell_types) >= 2) {
  # 计算平均 log2FC（初始化为0）
  avg_log2fc <- matrix(0, nrow = length(all_cell_types), ncol = length(all_cell_types),
                       dimnames = list(all_cell_types, all_cell_types))
  
  # 记录每个单元格有多少个有效样本
  n_effective <- matrix(0, nrow = length(all_cell_types), ncol = length(all_cell_types),
                         dimnames = list(all_cell_types, all_cell_types))
  
  # 存储合并 chi2 和自由度
  combined_chi2 <- matrix(0, nrow = length(all_cell_types), ncol = length(all_cell_types),
                          dimnames = list(all_cell_types, all_cell_types))
  combined_df <- matrix(0, nrow = length(all_cell_types), ncol = length(all_cell_types),
                        dimnames = list(all_cell_types, all_cell_types))
  
  for (sample_id in names(cn5_results)) {
    res <- cn5_results[[sample_id]]
    
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
  
  # 显著性标注
  signif_symbols <- matrix("", nrow = length(all_cell_types), ncol = length(all_cell_types),
                           dimnames = list(all_cell_types, all_cell_types))
  signif_symbols[combined_p < 0.001] <- "***"
  signif_symbols[combined_p < 0.01 & combined_p >= 0.001] <- "**"
  signif_symbols[combined_p < 0.05 & combined_p >= 0.01] <- "*"
  
  # 绘制平均热图
  p_avg <- pheatmap(avg_log2fc,
                    cluster_rows = TRUE,
                    cluster_cols = TRUE,
                    main = "CN5: Average Log2(Fold Change) across mPT samples",
                    fontsize = 10,
                    display_numbers = signif_symbols,
                    number_color = "black",
                    fontsize_number = 8,
                    filename = "mPT_CN5_all_samples_average_contact_log2fc.pdf",
                    width = 14, height = 12)
  
  cat("保存平均热图：mPT_CN5_all_samples_average_contact_log2fc.pdf\n")
  
  # 保存整合结果
  integrated_results <- list(
    all_cell_types = all_cell_types,
    avg_log2fc = avg_log2fc,
    combined_p = combined_p,
    n_effective = n_effective,
    n_samples = length(cn5_results),
    individual_results = cn5_results
  )
  
  saveRDS(integrated_results, file = "mPT_CN5_integrated_results.rds")
  cat("保存整合结果：mPT_CN5_integrated_results.rds\n")
}

#炎症型肿瘤正互作
library(ggplot2)
library(dplyr)
library(patchwork)

integrated_results <- readRDS("mPT_CN5_integrated_results.rds")

# 目标细胞类型
target_tumor_types <- c("Tumor_Immune-inflamed_1", "Tumor_Immune-inflamed_2")

# 获取所有细胞类型
all_types <- integrated_results$all_cell_types

# ===================== 定义绘图函数 =====================
plot_positive_interactions <- function(tumor_type, all_types, integrated_results) {
  
  if (!tumor_type %in% all_types) {
    return(NULL)
  }
  
  other_types <- all_types[all_types != tumor_type]
  
  df <- data.frame(
    CellType = other_types,
    Log2FC = integrated_results$avg_log2fc[tumor_type, other_types],
    P_value = integrated_results$combined_p[tumor_type, other_types],
    stringsAsFactors = FALSE
  )
  
  # 只保留正值且显著
  df_positive <- df[df$Log2FC > 0 & df$P_value < 0.05, ]
  
  if (nrow(df_positive) == 0) {
    return(NULL)
  }
  
  # 按 Log2FC 排序
  df_positive <- df_positive[order(df_positive$Log2FC, decreasing = TRUE), ]
  df_positive$CellType <- factor(df_positive$CellType, levels = rev(df_positive$CellType))
  
  # 提取肿瘤类型短名称
  short_name <- ifelse(tumor_type == "Tumor_Immune-inflamed_1", "Immune-inflamed_1", "Immune-inflamed_2")
  
  p <- ggplot(df_positive, aes(x = Log2FC, y = CellType)) +
    geom_bar(stat = "identity", fill = ifelse(tumor_type == "Tumor_Immune-inflamed_1", "#E15759", "#F28E2B")) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    labs(title = short_name,
         x = "Log2(Fold Change)",
         y = "") +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      axis.title.x = element_text(size = 11),
      plot.margin = margin(5, 5, 5, 5)
    )
  
  return(p)
}
p1 <- plot_positive_interactions("Tumor_Immune-inflamed_1", all_types, integrated_results)
p2 <- plot_positive_interactions("Tumor_Immune-inflamed_2", all_types, integrated_results)

# ===================== 左右拼接 =====================
if (!is.null(p1) && !is.null(p2)) {
  combined <- p1 + p2 +
    plot_annotation(
      title = "CN5: Positive interactions of immune-inflamed tumors",
      subtitle = "Log2FC > 0, P < 0.05",
      theme = theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        plot.subtitle = element_text(hjust = 0.5, size = 12)
      )
    ) +
    plot_layout(widths = c(0.5, 0.5))
  
  ggsave("mPT_CN5_Immune_inflamed_positive_interactions_combined.pdf", 
         combined, width = 10, height = 5, dpi = 300)
  cat("保存拼接图：mPT_CN5_Immune_inflamed_positive_interactions_combined.pdf\n")
  print(combined)
} else if (!is.null(p1)) {
  ggsave("mPT_CN5_Immune_inflamed_1_positive_interactions.pdf", p1, width = 5, height = 4, dpi = 300)
  print(p1)
} else if (!is.null(p2)) {
  ggsave("mPT_CN5_Immune_inflamed_2_positive_interactions.pdf", p2, width = 5, height = 4, dpi = 300)
  print(p2)
} else {
  cat("没有找到显著的吸引互作\n")
}
```
# mPT CN4 chat
```R
# ===================== 定义空间分析函数 =====================
run_spatial_analysis_cn4 <- function(coords_sub, labels_sub, sample_name, nperm = 1000) {
  
  cat("\n========================================\n")
  cat("开始分析样本：", sample_name, "\n")
  cat("CN4内细胞数量：", length(labels_sub), "\n")
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

samples <- unique(mPT$sample)
cat("mPT 样本:", paste(samples, collapse = ", "), "\n")

# 存储结果
cn4_results <- list()

for (sample_id in samples) {
  cat("\n\n########## 处理样本：", sample_id, " ##########\n")
  
  # 提取该样本中 CN4 的细胞
  sample_cells <- subset(mPT, subset = sample == sample_id & cn_detailed == "TLS_Zone")
  
  if (nrow(sample_cells@meta.data) == 0) {
    cat("警告：", sample_id, "CN4中没有细胞，跳过\n")
    next
  }
  
  cat("CN4细胞数量：", nrow(sample_cells@meta.data), "\n")
  
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
  
  # 过滤稀有细胞类型（数量 < 5 的过滤，因为CN4内细胞数较少）
  cell_counts <- table(labels)
  rare_threshold <- 5
  keep_types <- names(cell_counts[cell_counts >= rare_threshold])
  keep_cells <- labels %in% keep_types
  coords <- coords[keep_cells, , drop = FALSE]
  labels <- labels[keep_cells]
  
  cat("过滤稀有细胞后数量：", length(labels), "\n")
  
  if (length(labels) < 30) {
    cat("警告：", sample_id, "过滤后细胞太少，跳过\n")
    next
  }
  
  # 确保坐标是数值型
  coords <- as.matrix(coords)
  colnames(coords) <- c("x", "y")
  
  # 运行分析
  result <- run_spatial_analysis_cn4(coords, labels, sample_id, nperm = 500)
  
  if (!is.null(result)) {
    cn4_results[[sample_id]] <- result
    
    # 绘制并保存热图
    if (nrow(result$mat_log2fc) > 1 && ncol(result$mat_log2fc) > 1) {
      mat_plot <- result$mat_log2fc
      
      # 处理 Inf
      mat_plot[is.infinite(mat_plot) & mat_plot < 0] <- -10
      mat_plot[is.infinite(mat_plot) & mat_plot > 0] <- 10
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
      pdf(paste0("mPT_CN4_", sample_id, "_contact_log2fc.pdf"), width = 12, height = 10)
      pheatmap(mat_plot,
               cluster_rows = TRUE,
               cluster_cols = TRUE,
               main = paste("CN4 - Log2(FC) -", sample_id),
               fontsize = 10,
               display_numbers = signif_symbols,
               number_color = "black",
               fontsize_number = 8)
      dev.off()
      
      cat("保存热图：mPT_CN4_", sample_id, "_contact_log2fc.pdf\n", sep = "")
    }
  }
}

saveRDS(cn4_results, file = "mPT_CN4_cell_interaction_results.rds")
cat("\n\n========== 分析完成 ==========\n")
cat("成功分析样本数：", length(cn4_results), "\n")

# ===================== 整合多样本CN4结果 =====================
all_cell_types <- unique(unlist(lapply(cn4_results, function(x) x$cell_types)))
cat("总细胞类型数量：", length(all_cell_types), "\n")
cat("所有细胞类型：", paste(all_cell_types, collapse = ", "), "\n")

if (length(all_cell_types) >= 2) {
  # 计算平均 log2FC（初始化为0）
  avg_log2fc <- matrix(0, nrow = length(all_cell_types), ncol = length(all_cell_types),
                       dimnames = list(all_cell_types, all_cell_types))
  
  # 记录每个单元格有多少个有效样本
  n_effective <- matrix(0, nrow = length(all_cell_types), ncol = length(all_cell_types),
                         dimnames = list(all_cell_types, all_cell_types))
  
  # 存储合并 chi2 和自由度
  combined_chi2 <- matrix(0, nrow = length(all_cell_types), ncol = length(all_cell_types),
                          dimnames = list(all_cell_types, all_cell_types))
  combined_df <- matrix(0, nrow = length(all_cell_types), ncol = length(all_cell_types),
                        dimnames = list(all_cell_types, all_cell_types))
  
  for (sample_id in names(cn4_results)) {
    res <- cn4_results[[sample_id]]
    
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
  
  # 显著性标注
  signif_symbols <- matrix("", nrow = length(all_cell_types), ncol = length(all_cell_types),
                           dimnames = list(all_cell_types, all_cell_types))
  signif_symbols[combined_p < 0.001] <- "***"
  signif_symbols[combined_p < 0.01 & combined_p >= 0.001] <- "**"
  signif_symbols[combined_p < 0.05 & combined_p >= 0.01] <- "*"
  
  # 绘制平均热图
  p_avg <- pheatmap(avg_log2fc,
                    cluster_rows = TRUE,
                    cluster_cols = TRUE,
                    main = "CN4: Average Log2(Fold Change) across mPT samples",
                    fontsize = 10,
                    display_numbers = signif_symbols,
                    number_color = "black",
                    fontsize_number = 8,
                    filename = "mPT_CN4_all_samples_average_contact_log2fc.pdf",
                    width = 14, height = 12)
  
  cat("保存平均热图：mPT_CN4_all_samples_average_contact_log2fc.pdf\n")
  
  # 保存整合结果
  integrated_results <- list(
    all_cell_types = all_cell_types,
    avg_log2fc = avg_log2fc,
    combined_p = combined_p,
    n_effective = n_effective,
    n_samples = length(cn4_results),
    individual_results = cn4_results
  )
  
  saveRDS(integrated_results, file = "mPT_CN4_integrated_results.rds")
  cat("保存整合结果：mPT_CN4_integrated_results.rds\n")
}