# metadata pretreatment
```R
obj <- readRDS("YA2025263-1_fin.rds")
library(Seurat)

sample_tissue <- data.frame(
  sample = c("A1", "A2", "A3", "A4", "A5", 
             "B1", "B2", "B3", "B4", "B5",
             "C1", "C2", "C3", "C4", "C5",
             "D1", "D2", "D3", "D4", "D5"),
  tissue = c("NRT", "NRLN", "", "", "NRT",
             "NRLN", "NRT", "NRLN", "RT", "RLN_N",
             "RLN_P", "RT", "RLN_N", "RLN_P", "RT",
             "RLN_N", "RLN_P", "", "", "")
)

obj$tissue <- sample_tissue$tissue[match(obj$sample, sample_tissue$sample)]

library(ggplot2)
library(dplyr)

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