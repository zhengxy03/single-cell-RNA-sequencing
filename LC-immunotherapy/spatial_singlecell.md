# metadata pretreatment
```R
obj <- readRDS("YA2025263-1_fin.rds")
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
                       "Malignant cells" = "unknown",
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
t_cells <- subset(obj,subset=CellType %in% c("T cells","NK cells"))
#saveRDS(t_cells,file="t_NK.rds")
t_cells <- subset(t_cells, subset = tissue %in% c("NRT","RT","NRLN","RLN_N","RLN_P"))
t_cells <- NormalizeData(t_cells)
t_cells <- FindVariableFeatures(t_cells, nfeatures = 2000)
hvgs <- VariableFeatures(t_cells)
t_cells <- ScaleData(t_cells, features = hvgs)
t_cells <- RunPCA(t_cells, features = hvgs, npcs = 20)
#p <- ElbowPlot(t_cells, ndims = 30)
#ggsave("p3.png",plot=p)
t_cells <- RunUMAP(t_cells, dims = 1:20)
t_cells <- FindNeighbors(t_cells, dims = 1:20)
t_cells <- FindClusters(t_cells, resolution = 0.8)
saveRDS(t_cells,file="t_cells_new.rds")

t_cell_markers <- FindAllMarkers(t_cells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "wilcox")
t_significant_markers <- subset(t_cell_markers, p_val_adj < 0.05)
#write.csv(t_significant_markers, "t_all_marker.csv")
t_significant_markers <- t_significant_markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.csv(t_significant_markers, "t_top_marker_20.csv")


t_cell_markers <- FindAllMarkers(t_cells, only.pos = TRUE, min.pct = 0.01, logfc.threshold = 0.1, test.use = "wilcox")
t_significant_markers <- subset(t_cell_markers, p_val_adj < 1)
write.csv(t_significant_markers, "t_all_marker.csv")
t_significant_markers <- t_significant_markers %>% group_by(cluster) %>% top_n(n = 200, wt = avg_log2FC)
write.csv(t_significant_markers, "t_top_marker_200.csv")


t_cells <- subset(t_cells,subset=seurat_clusters %in% c(0,2,7,9,11))
#运行预处理
#saveRDS(t_cells,file="t_cells_rm.rds")
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

t_cells <- subset(t_cells,subset=seurat_clusters %in% c(0,1,5,7,9,10,11,13,14,15,16,17,18,20,21,22,23,25,26,27,28,29,30,32))
t_cells <- subset(t_cells,subset=seurat_clusters %in% c(0,1,5,7,13,22,23,25,26,29,30,32))
t_cells <- subset(t_cells,subset=seurat_clusters %in% c(0,2,3,11,12))
t_cells <- FindClusters(t_cells, resolution = 0.9)
saveRDS(t_cells,file="t_NK_rem.rds")
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
npg_extended <- colorRampPalette(npg_pal)(18)

cell_types <- as.character(unique(obj@meta.data$CellType))
num_legend_items <- length(cell_types)  # 图例的个数
max_label_length <- max(nchar(cell_types))  # 图例名称的最大长度

base_width <- 3000  # 基础宽度
base_height <- 3000  # 基础高度
legend_width_factor <- 100  # 每个图例项增加的宽度
label_length_factor <- 10  # 每个字符增加的宽度
dynamic_width <- base_width + (num_legend_items * legend_width_factor) + (max_label_length * label_length_factor)
pdf("annotation.pdf", width = dynamic_width/300, height = base_height/300)
DimPlot(obj, reduction = "umap", label = TRUE, pt.size = 1, group.by = "CellType", label.size = 4) +
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
npg_extended <- colorRampPalette(npg_pal)(18)


prop_data <- prop_data %>%
  mutate(CellType = as.character(CellType))

# 绘制堆叠条形图
p <- ggplot(prop_data, aes(x = tissue, y = proportion, fill = CellType)) +
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
pdf("obj_celltype_distribution_tumor.pdf", width = 4, height = 6)
print(p)
dev.off()

prop_data <- obj@meta.data %>%
  group_by(tissue, CellType) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(tissue) %>%
  mutate(proportion = count / sum(count) * 100) %>%
  filter(tissue %in% c("metLN", "negLN")) %>%
  mutate(tissue = factor(tissue, levels = c("metLN", "negLN")))

library(ggplot2)
library(ggsci)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(18)


prop_data <- prop_data %>%
  mutate(CellType = as.character(CellType))

# 绘制堆叠条形图
p <- ggplot(prop_data, aes(x = tissue, y = proportion, fill = CellType)) +
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
pdf("obj_celltype_distribution_RLN_P_RLN_N_NRLN.pdf", width = 4, height = 6)
print(p)
dev.off()
```
# malignant
```R
Malignant <- subset(obj,subset=CellType=="Malignant cells")
Malignant <- subset(Malignant,subset=seurat_clusters %in% c(0,1,2,3,5,7,9,10))
Malignant <- NormalizeData(Malignant)
Malignant <- FindVariableFeatures(Malignant, nfeatures = 2000)
hvgs <- VariableFeatures(Malignant)
Malignant <- ScaleData(Malignant, features = hvgs)
Malignant <- RunPCA(Malignant, features = hvgs, npcs = 50)
#p <- ElbowPlot(Malignant, ndims = 30)
#ggsave("p3.png",plot=p)

Malignant <- FindNeighbors(Malignant, dims = 1:30)
Malignant <- FindClusters(Malignant, resolution = 0.35)
Malignant <- RunUMAP(Malignant, dims = 1:30)

library(ggplot2)
library(ggsci)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(12)

seurat_clusters <- as.character(unique(Malignant@meta.data$seurat_clusters))  # 转换为字符向量
num_legend_items <- length(seurat_clusters)  # 图例的个数
max_label_length <- max(nchar(seurat_clusters))  # 图例名称的最大长度

base_width <- 3000  # 基础宽度
base_height <- 3000  # 基础高度
legend_width_factor <- 100  # 每个图例项增加的宽度
label_length_factor <- 10  # 每个字符增加的宽度
dynamic_width <- base_width + (num_legend_items * legend_width_factor) + (max_label_length * label_length_factor)
pdf("Malignant_unknown_clusters.pdf", width = dynamic_width/300, height = base_height/300)  # 转换为英寸

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
write.csv(Malignant_significant_markers, "Malignant_unknown_top_marker_50_remove.csv")

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

identity_mapping <- c(
  "0" = "Metabolic",
  "1" = "Inflamed_1",
  "2" = "Inflamed_2",
  "3" = "Inflamed_3",
  "4" = "Inflamed_4",
  "5" = "Stem-like",  
  "6" = "Proliferative",
  "7" = "Inflamed_5",
  "8" = "Club cell",
  "9" = "Differentiated"
)
sub_cell_type <- identity_mapping[Malignant@meta.data$seurat_clusters]
Malignant@meta.data$sub_cell_type <- sub_cell_type


npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(10)

cell_types <- as.character(unique(Malignant@meta.data$sub_cell_type))
num_legend_items <- length(cell_types)  # 图例的个数
max_label_length <- max(nchar(cell_types))  # 图例名称的最大长度

base_width <- 3000  # 基础宽度
base_height <- 3000  # 基础高度
legend_width_factor <- 100  # 每个图例项增加的宽度
label_length_factor <- 10  # 每个字符增加的宽度
dynamic_width <- base_width + (num_legend_items * legend_width_factor) + (max_label_length * label_length_factor)
pdf("malignant_unknown_annotation.pdf", width = dynamic_width/300, height = base_height/300)
DimPlot(Malignant, reduction = "umap", label = TRUE, pt.size = 1, group.by = "sub_cell_type", label.size = 4) +
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
prop_data <- Malignant@meta.data %>%
  group_by(tissue, sub_cell_type) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(tissue) %>%
  mutate(proportion = count / sum(count) * 100) %>%
  filter(tissue %in% c("mPT", "nmPT")) %>%
  mutate(tissue = factor(tissue, levels = c("mPT", "nmPT")))

library(ggplot2)
library(ggsci)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(10)


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
pdf("Malignant_unknown_celltype_distribution_tumor.pdf", width = 4, height = 6)
print(p)
dev.off()

prop_data <- Malignant@meta.data %>%
  group_by(tissue, sub_cell_type) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(tissue) %>%
  mutate(proportion = count / sum(count) * 100) %>%
  filter(tissue %in% c("metLN", "negLN")) %>%
  mutate(tissue = factor(tissue, levels = c("metLN", "negLN")))

library(ggplot2)
library(ggsci)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(12)


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
pdf("Malignant_celltype_distribution_RLN_P_RLN_N_NRLN.pdf", width = 4, height = 6)
print(p)
dev.off()

prop_data <- Malignant@meta.data %>%
  group_by(tissue, sub_cell_type) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(tissue) %>%
  mutate(proportion = count / sum(count) * 100) %>%
  filter(tissue %in% c("mPT", "metLN")) %>%
  mutate(tissue = factor(tissue, levels = c("mPT", "metLN")))

library(ggplot2)
library(ggsci)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(12)


prop_data <- prop_data %>%
  mutate(sub_cell_type = as.character(sub_cell_type))

# 绘制堆叠条形图
p <- ggplot(prop_data, aes(x = tissue, y = proportion, fill = sub_cell_type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = npg_extended) +
  labs(title = "Cell Type Distribution: mPT vs metLN",
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
pdf("Malignant_celltype_distribution_mPT_metLN.pdf", width = 4, height = 6)
print(p)
dev.off()

p <- FeaturePlot(Malignant,feature="VIM")
ggsave("feature_vim.png",plot=p)
p <- VlnPlot(Malignant,feature="VIM")
ggsave("vln_vim.png",plot=p)
```
# trajectory
```R
library(Seurat)
library(monocle)
library(ggplot2)
library(ggsci)
trace('project2MST', edit = T, where = asNamespace("monocle"))



#抽样
expression_matrix <- LayerData(Malignant, assay = "SCT", layer = "data")
cell_metadata <- Malignant@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(expression_matrix))
rownames(gene_annotation) <- rownames(expression_matrix)

# 随机选择部分细胞
set.seed(123)
subset_cells <- sample(colnames(expression_matrix), size = 20000)  # 选择 10,000 个细胞
expression_matrix <- expression_matrix[, subset_cells]
cell_metadata <- cell_metadata[subset_cells, ]

# 创建 CellDataSet 对象
cds <- newCellDataSet(expression_matrix,
                      phenoData = new("AnnotatedDataFrame", data = cell_metadata),
                      featureData = new("AnnotatedDataFrame", data = gene_annotation),
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())

# 归一化并估计离散度
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

# 检测基因并选择高变基因
cds <- detectGenes(cds, min_expr = 0.1)
disp_table <- dispersionTable(cds)
ordering_genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
cds <- setOrderingFilter(cds, ordering_genes)

cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')

# 轨迹推断
cds <- orderCells(cds)
cds <- orderCells(cds, root_state  = 4)

library(ggplot2)
library(ggsci)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(12)

p3 <- plot_cell_trajectory(cds, color_by = "sub_cell_type") + scale_color_manual(values = npg_extended) +
  theme(legend.text = element_text(size = 18),
        legend.title = element_blank())

png("traj_malignant_celltype.png", width = 6000, height = 3000, res = 300)
print(p3)
dev.off()

p4 <- plot_cell_trajectory(cds, color_by = "sub_cell_type") + facet_wrap("~sub_cell_type", nrow = 2) + scale_color_manual(values = npg_extended) +
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
ggsave("malignant_pseudo.png", plot = p5, width = 16, height = 6, dpi = 300)
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
Malignant <- subset(Malignant,subset=seurat_clusters %in% c(0,1,2,3,5,7,9,10))
fib <- NormalizeData(fib)
fib <- FindVariableFeatures(fib, nfeatures = 2000)
hvgs <- VariableFeatures(fib)
fib <- ScaleData(fib, features = hvgs)
fib <- RunPCA(fib, features = hvgs, npcs = 50)
fib <- FindNeighbors(fib, dims = 1:30)
fib <- FindClusters(fib, resolution = 0.3)
fib <- RunUMAP(fib, dims = 1:30)

library(ggplot2)
library(ggsci)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(10)

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
fib_markers <- FindAllMarkers(fib, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "wilcox")
fib_significant_markers <- subset(fib_markers, p_val_adj < 0.05)
#write.csv(fib_significant_markers, "fib_all_marker.csv")
fib_significant_markers <- fib_significant_markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
write.csv(fib_significant_markers, "fib_top_marker_50.csv")


identity_mapping <- c(
  "0" = "myCAF",
  "1" = "apCAF",
  "2" = "eCAF",
  "3" = "apCAF",
  "4" = "iCAF",
  "5" = "matCAF"
)

sub_cell_type <- identity_mapping[fib@meta.data$seurat_clusters]
fib@meta.data$sub_cell_type <- sub_cell_type


npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(5)

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

#比例
library(tidyverse)

# 计算各亚群在不同组织中的占比
prop_data <- fib@meta.data %>%
  group_by(tissue, sub_cell_type) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(tissue) %>%
  mutate(proportion = count / sum(count) * 100) %>%
  filter(tissue %in% c("mPT", "nmPT")) %>%
  mutate(tissue = factor(tissue, levels = c("mPT", "nmPT")))

library(ggplot2)
library(ggsci)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(5)


prop_data <- prop_data %>%
  mutate(sub_cell_type = as.character(sub_cell_type))

# 绘制堆叠条形图
p <- ggplot(prop_data, aes(x = tissue, y = proportion, fill = sub_cell_type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = npg_pal) +
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
pdf("fib_celltype_distribution_tumor.pdf", width = 4, height = 6)
print(p)
dev.off()

prop_data <- fib@meta.data %>%
  group_by(tissue, sub_cell_type) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(tissue) %>%
  mutate(proportion = count / sum(count) * 100) %>%
  filter(tissue %in% c("metLN", "negLN")) %>%
  mutate(tissue = factor(tissue, levels = c("metLN", "negLN")))

library(ggplot2)
library(ggsci)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(5)


prop_data <- prop_data %>%
  mutate(sub_cell_type = as.character(sub_cell_type))

# 绘制堆叠条形图
p <- ggplot(prop_data, aes(x = tissue, y = proportion, fill = sub_cell_type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = npg_pal) +
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
pdf("fib_celltype_distribution_RLN_P_RLN_N_NRLN.pdf", width = 4, height = 6)
print(p)
dev.off()

prop_data <- fib@meta.data %>%
  group_by(tissue, sub_cell_type) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(tissue) %>%
  mutate(proportion = count / sum(count) * 100) %>%
  filter(tissue %in% c("mPT", "metLN")) %>%
  mutate(tissue = factor(tissue, levels = c("mPT", "metLN")))

library(ggplot2)
library(ggsci)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(5)


prop_data <- prop_data %>%
  mutate(sub_cell_type = as.character(sub_cell_type))

# 绘制堆叠条形图
p <- ggplot(prop_data, aes(x = tissue, y = proportion, fill = sub_cell_type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = npg_pal) +
  labs(title = "Cell Type Distribution: mPT vs metLN",
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
pdf("fib_celltype_distribution_mPT_metLN.pdf", width = 4, height = 6)
print(p)
dev.off()
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
macro <- FindVariableFeatures(macro, nfeatures = 2000)
hvgs <- VariableFeatures(macro)
macro <- ScaleData(macro, features = hvgs)
macro <- RunPCA(macro, features = hvgs, npcs = 50)
p <- ElbowPlot(macro, ndims = 30)
ggsave("elbow_macro.png",plot=p)
macro <- FindNeighbors(macro, dims = 1:20)
macro <- FindClusters(macro, resolution = 0.3)
macro <- RunUMAP(macro, dims = 1:20)

library(ggplot2)
library(ggsci)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(6)

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

sub_cell_type <- identity_mapping[macro@meta.data$seurat_clusters]
macro@meta.data$sub_cell_type <- sub_cell_type


npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(5)

cell_types <- as.character(unique(macro@meta.data$sub_cell_type))
num_legend_items <- length(cell_types)  # 图例的个数
max_label_length <- max(nchar(cell_types))  # 图例名称的最大长度

base_width <- 3000  # 基础宽度
base_height <- 3000  # 基础高度
legend_width_factor <- 100  # 每个图例项增加的宽度
label_length_factor <- 10  # 每个字符增加的宽度
dynamic_width <- base_width + (num_legend_items * legend_width_factor) + (max_label_length * label_length_factor)
pdf("macro_annotation.pdf", width = dynamic_width/300, height = base_height/300)
DimPlot(macro, reduction = "umap", label = TRUE, pt.size = 1, group.by = "sub_cell_type", label.size = 4) +
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
npg_extended <- colorRampPalette(npg_pal)(5)


prop_data <- prop_data %>%
  mutate(sub_cell_type = as.character(sub_cell_type))

# 绘制堆叠条形图
p <- ggplot(prop_data, aes(x = tissue, y = proportion, fill = sub_cell_type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = npg_pal) +
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
npg_extended <- colorRampPalette(npg_pal)(5)


prop_data <- prop_data %>%
  mutate(sub_cell_type = as.character(sub_cell_type))

# 绘制堆叠条形图
p <- ggplot(prop_data, aes(x = tissue, y = proportion, fill = sub_cell_type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = npg_pal) +
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
pdf("macro_celltype_distribution_RLN_P_RLN_N_NRLN.pdf", width = 4, height = 6)
print(p)
dev.off()

prop_data <- macro@meta.data %>%
  group_by(tissue, sub_cell_type) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(tissue) %>%
  mutate(proportion = count / sum(count) * 100) %>%
  filter(tissue %in% c("mPT", "metLN")) %>%
  mutate(tissue = factor(tissue, levels = c("mPT", "metLN")))

library(ggplot2)
library(ggsci)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(5)


prop_data <- prop_data %>%
  mutate(sub_cell_type = as.character(sub_cell_type))

# 绘制堆叠条形图
p <- ggplot(prop_data, aes(x = tissue, y = proportion, fill = sub_cell_type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = npg_pal) +
  labs(title = "Cell Type Distribution: mPT vs metLN",
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
pdf("macro_celltype_distribution_mPT_metLN.pdf", width = 4, height = 6)
print(p)
dev.off()
```
# mPT的主要细胞类型互作
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
                       "Basal cells" = "Malignant cells"
)

# 检查
table(obj$CellType)



mPT <- subset(obj,subset=tissue=="mPT")
coords <- mPT@meta.data[, c("CenterX_global_px", "CenterY_global_px")]
head(coords)

# 获取细胞barcodes（就是行名）
barcodes <- rownames(coords)

labels <- as.character(mPT@meta.data$CellType)
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
Malignant <- readRDS("malignant_anno.rds")
fib <- readRDS("fib_anno.rds")
macro <- readRDS("macro_anno_no_rm.rds")

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
    top3 <- names(pct_sorted)[1:3]
    top3_pct <- round(pct_sorted[1:3], 1)
    
    cat("\n区域", reg, ":\n")
    cat("  主要细胞类型: ", paste(top3, "(", top3_pct, "%)", collapse = ", "), "\n")
    cat("  细胞数:", sum(region_comp[reg, ]), "\n")
  }
}

# ===================== k=3 统一命名 =====================
# ===================== 先应用区域命名 =====================
region_naming_k3 <- list()

region_naming_k3[["C5"]] <- c(
  "R1" = "Mast_cell_zone",
  "R2" = "Epithelial_zone", 
  "R3" = "Macrophage_zone"
)

region_naming_k3[["B4"]] <- c(
  "R1" = "Mixed_epithelial",
  "R2" = "Tumor_core",
  "R3" = "Vascular_zone"
)

region_naming_k3[["D3"]] <- c(
  "R1" = "Tumor_core",
  "R2" = "Tumor_edge",
  "R3" = "Stromal_zone"
)

region_naming_k3[["C2"]] <- c(
  "R1" = "Mixed_zone",
  "R2" = "Tumor_zone",
  "R3" = "Vascular_zone"
)

# 应用命名
for (sample_id in samples) {
  if (sample_id %in% names(region_naming_k3)) {
    region_df <- region_list[[sample_id]]
    mapping <- region_naming_k3[[sample_id]]
    
    region_df$region_type <- mapping[region_df$region]
    region_list[[sample_id]] <- region_df
    
    cat("样本", sample_id, "命名完成\n")
    print(table(region_df$region_type))
  }
}

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
merged_results <- list()

for (reg_type in all_region_types) {
  cat("\n合并区域:", reg_type, "\n")
  
  # 收集该区域类型的所有结果
  region_results <- list()
  
  for (sample_id in names(all_results)) {
    if (reg_type %in% names(all_results[[sample_id]])) {
      region_results[[sample_id]] <- all_results[[sample_id]][[reg_type]]
    }
  }
  
  n_samples <- length(region_results)
  cat("  样本数:", n_samples, "\n")
  
  # 只处理有多样本的区域
  if (n_samples < 2) {
    cat("  跳过（只有1个样本）\n")
    next
  }
  
  # 找出所有样本共有的细胞类型
  common_types <- Reduce(intersect, lapply(region_results, function(x) x$cell_types))
  
  if (length(common_types) < 2) {
    cat("  共有细胞类型太少，跳过\n")
    next
  }
  
  cat("  共有细胞类型数:", length(common_types), "\n")
  
  # 计算平均 log2FC 和合并 p 值
  avg_log2fc <- matrix(0, nrow = length(common_types), ncol = length(common_types),
                       dimnames = list(common_types, common_types))
  combined_chi2 <- matrix(0, nrow = length(common_types), ncol = length(common_types))
  
  for (sample_id in names(region_results)) {
    res <- region_results[[sample_id]]
    
    mat_log2fc <- res$mat_log2fc[common_types, common_types]
    mat_p <- res$mat_p[common_types, common_types]
    
    # 处理无穷值和NA
    mat_log2fc[is.infinite(mat_log2fc) & mat_log2fc > 0] <- 10
    mat_log2fc[is.infinite(mat_log2fc) & mat_log2fc < 0] <- -10
    mat_log2fc[is.na(mat_log2fc)] <- 0
    mat_p[is.na(mat_p)] <- 1
    
    avg_log2fc <- avg_log2fc + mat_log2fc
    combined_chi2 <- combined_chi2 + (-2 * log(mat_p + 1e-10))
  }
  
  avg_log2fc <- avg_log2fc / n_samples
  combined_p <- pchisq(combined_chi2, df = 2 * n_samples, lower.tail = FALSE)
  
  # 最终处理一次
  avg_log2fc[is.infinite(avg_log2fc) & avg_log2fc > 0] <- 10
  avg_log2fc[is.infinite(avg_log2fc) & avg_log2fc < 0] <- -10
  avg_log2fc[is.na(avg_log2fc)] <- 0
  combined_p[is.na(combined_p)] <- 1
  
  # 显著性标注
  signif_symbols <- matrix("", nrow = length(common_types), ncol = length(common_types),
                           dimnames = list(common_types, common_types))
  signif_symbols[combined_p < 0.001] <- "***"
  signif_symbols[combined_p < 0.01 & combined_p >= 0.001] <- "**"
  signif_symbols[combined_p < 0.05 & combined_p >= 0.01] <- "*"
  
  # 绘制热图
  pheatmap(avg_log2fc,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           display_numbers = signif_symbols,
           number_color = "black",
           fontsize_number = 6,
           main = paste("Cell-cell interactions in", reg_type, "region"),
           filename = paste0("mPT_region_", gsub("/", "_", reg_type), "_merged.pdf"),
           width = 12, height = 10)
  
  cat("  保存热图: mPT_region_", reg_type, "_merged.pdf\n", sep = "")
}

cat("\n===== 完成！=====\n")

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
