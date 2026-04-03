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
npg_extended <- colorRampPalette(npg_pal)(17)

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
npg_extended <- colorRampPalette(npg_pal)(17)


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
npg_extended <- colorRampPalette(npg_pal)(17)


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
Malignant <- FindClusters(Malignant, resolution = 0.6)
Malignant <- RunUMAP(Malignant, dims = 1:30)

library(ggplot2)
library(ggsci)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(14)#0.6

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
write.csv(Malignant_significant_markers, "Malignant_unknown_top_marker_50_0.6.csv")
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
sub_cell_type <- identity_mapping[Malignant@meta.data$seurat_clusters]
Malignant@meta.data$sub_cell_type <- sub_cell_type

library(ggplot2)
library(ggsci)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(14)

cell_types <- as.character(unique(Malignant@meta.data$sub_cell_type))
num_legend_items <- length(cell_types)  # 图例的个数
max_label_length <- max(nchar(cell_types))  # 图例名称的最大长度

base_width <- 3000  # 基础宽度
base_height <- 3000  # 基础高度
legend_width_factor <- 100  # 每个图例项增加的宽度
label_length_factor <- 10  # 每个字符增加的宽度
dynamic_width <- base_width + (num_legend_items * legend_width_factor) + (max_label_length * label_length_factor)
pdf("mPT_malignant_unknown_annotation-0.6.pdf", width = dynamic_width/300, height = base_height/300)
DimPlot(mPT, reduction = "umap", label = TRUE, pt.size = 1, group.by = "sub_cell_type", label.size = 4) +
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
npg_extended <- colorRampPalette(npg_pal)(14)


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
pdf("Malignant_unknown_celltype_distribution_tumor0.6.pdf", width = 4, height = 6)
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
npg_extended <- colorRampPalette(npg_pal)(14)


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
pdf("Malignant_unknown_celltype_distribution_RLN_P_RLN_N_NRLN-0.6.pdf", width = 4, height = 6)
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
npg_extended <- colorRampPalette(npg_pal)(14)


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
pdf("Malignant_unknown_celltype_distribution_mPT_metLN-0.6.pdf", width = 4, height = 6)
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
#trace('project2MST', edit = T, where = asNamespace("monocle"))

Malignant <- readRDS("malignant_unknown_anno_0.6.rds")

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

# 保存
saveRDS(cds, "malignant_unknown_cds_0.6_20000cells_hvg.rds")
# 轨迹推断
cds <- orderCells(cds)
cds <- orderCells(cds, root_state  = 4)

library(ggplot2)
library(ggsci)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(14)

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
  "1" = "apCAF_1",
  "2" = "eCAF",
  "3" = "apCAF_2",
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

identity_mapping <- c(
  "0" = "Lipid-associated_Macro_1",
  "1" = "Lipid-associated_Macro_2",
  "2" = "Inflammatory_Macro",
  "3" = "Matrix-remodeling_Macro",
  "4" = "Stress-response_Macro",
  "5" = "Proliferative-Macro"
)


sub_cell_type <- identity_mapping[macro@meta.data$seurat_clusters]
macro@meta.data$sub_cell_type <- sub_cell_type

library(ggplot2)
library(ggsci)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(6)

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
ggsave("mPT_spatial_contact_log2fc_maligvsothers.pdf", plot = p, width = 12, height = 10)

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


# ... 运行完整的互作分析，得到 mat_log2fc_nmPT, mat_p_nmPT ...

# ===================== 比较 mPT vs nmPT =====================
mat_log2fc_mPT <- mat_log2fc_clean
mat_p_mPT <- mat_p
mat_log2fc_nmPT <- mat_log2fc_clean
mat_p_nmPT <- mat_p
# 1. 找出两种组织共有的细胞类型
common_types <- intersect(rownames(mat_log2fc_mPT), rownames(mat_log2fc_nmPT))

# 2. 计算差异矩阵（mPT - nmPT）
mat_diff <- mat_log2fc_mPT[common_types, common_types] - 
            mat_log2fc_nmPT[common_types, common_types]

# 3. 裁剪到合理范围
#mat_diff[mat_diff > 3] <- 3
#mat_diff[mat_diff < -3] <- -3

# 4. 计算差异的显著性（Fisher's method 合并 p 值）
combined_chi2 <- matrix(0, nrow = length(common_types), ncol = length(common_types),
                        dimnames = list(common_types, common_types))

for (i in 1:length(common_types)) {
  for (j in 1:length(common_types)) {
    p1 <- mat_p_mPT[common_types[i], common_types[j]]
    p2 <- mat_p_nmPT[common_types[i], common_types[j]]
    combined_chi2[i, j] <- -2 * (log(p1 + 1e-10) + log(p2 + 1e-10))
  }
}
combined_p <- pchisq(combined_chi2, df = 4, lower.tail = FALSE)

# 5. 显著性标注
signif_detail <- matrix("", nrow = length(common_types), ncol = length(common_types),
                        dimnames = list(common_types, common_types))
signif_detail[combined_p < 0.001 & mat_diff > 0] <- "***"
signif_detail[combined_p < 0.01 & combined_p >= 0.001 & mat_diff > 0] <- "**"
signif_detail[combined_p < 0.05 & combined_p >= 0.01 & mat_diff > 0] <- "*"
signif_detail[combined_p < 0.05 & mat_diff < 0] <- "*"  # 负值也标注

# 绘制热图（带显著性标注）
p <- pheatmap(mat_diff,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         display_numbers = signif_detail,  # 关键：显示显著性标注
         number_color = "black",
         fontsize_number = 8,
         main = "Difference in cell-cell interactions (mPT - nmPT)\nRed: stronger in mPT, Blue: stronger in nmPT\n* p < 0.05, ** p < 0.01, *** p < 0.001",
         fontsize = 10,
         breaks = seq(-3, 3, length.out = 100),
         color = colorRampPalette(c("blue", "white", "red"))(100))

ggsave("mPT_vs_nmPT_interaction_diff.pdf", plot = p, width = 14, height = 12)
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

# ===================== 找出代谢型 =====================
cell_types_mPT <- rownames(mat_log2fc_mPT)
cell_types_nmPT <- rownames(mat_log2fc_nmPT)

metabolic1 <- cell_types_mPT[grepl("Metabolic_1", cell_types_mPT)]
metabolic2 <- cell_types_mPT[grepl("Metabolic_2", cell_types_mPT)]

cat("代谢1:", metabolic1, "\n")
cat("代谢2:", metabolic2, "\n")

# ===================== 提取互作 =====================
# mPT
inter_mPT_met1 <- mat_log2fc_mPT[metabolic1, ]
inter_mPT_met2 <- mat_log2fc_mPT[metabolic2, ]

# nmPT
inter_nmPT_met1 <- mat_log2fc_nmPT[metabolic1, ]
inter_nmPT_met2 <- mat_log2fc_nmPT[metabolic2, ]

# 排除自身
inter_mPT_met1 <- inter_mPT_met1[!names(inter_mPT_met1) %in% c(metabolic1, metabolic2)]
inter_mPT_met2 <- inter_mPT_met2[!names(inter_mPT_met2) %in% c(metabolic1, metabolic2)]
inter_nmPT_met1 <- inter_nmPT_met1[!names(inter_nmPT_met1) %in% c(metabolic1, metabolic2)]
inter_nmPT_met2 <- inter_nmPT_met2[!names(inter_nmPT_met2) %in% c(metabolic1, metabolic2)]

# ===================== 找共同细胞类型 =====================
common_cells <- intersect(names(inter_mPT_met1), names(inter_nmPT_met1))

# ===================== 计算差异 =====================
diff_met1 <- inter_mPT_met1[common_cells] - inter_nmPT_met1[common_cells]
diff_met2 <- inter_mPT_met2[common_cells] - inter_nmPT_met2[common_cells]

# 排序
top_met1_mPT <- sort(diff_met1, decreasing = TRUE)[1:15]
top_met2_mPT <- sort(diff_met2, decreasing = TRUE)[1:15]

cat("\n========== 代谢1在 mPT 中更强的互作 ==========\n")
print(top_met1_mPT)

cat("\n========== 代谢2在 mPT 中更强的互作 ==========\n")
print(top_met2_mPT)

final_diff <- diff_met2[common_cells] - diff_met1[common_cells]

# 排序
top_met2_stronger <- sort(final_diff, decreasing = TRUE)[1:15]
top_met1_stronger <- sort(final_diff, decreasing = FALSE)[1:15]

cat("\n========== 代谢2更强的互作（基于转移特异性差异）==========\n")
print(top_met2_stronger)

cat("\n========== 代谢1更强的互作（基于转移特异性差异）==========\n")
print(top_met1_stronger)

# ===================== 条形图 =====================
library(ggplot2)

# 创建数据框
df_met2 <- data.frame(
  cell_type = names(top_met2_stronger),
  diff = as.numeric(top_met2_stronger),
  stronger = "Metabolic2"
)

df_met1 <- data.frame(
  cell_type = names(top_met1_stronger),
  diff = as.numeric(top_met1_stronger),
  stronger = "Metabolic1"
)

df_plot <- rbind(df_met2, df_met1)

# 条形图
p <- ggplot(df_plot, aes(x = reorder(cell_type, diff), y = diff, fill = stronger)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("Metabolic1" = "#E69F00", "Metabolic2" = "#D55E00")) +
  labs(title = "Functional divergence of metabolic subtypes (based on mPT vs nmPT difference)",
       x = "Cell type", 
       y = "(Metabolic2_mPT - Metabolic2_nmPT) - (Metabolic1_mPT - Metabolic1_nmPT)",
       fill = "Stronger with") +
  theme_minimal() +
  theme(legend.position = "bottom",
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        plot.title = element_text(hjust = 0.5))

ggsave("metabolic_divergence_based_on_transfer_specificity.pdf", plot = p, width = 10, height = 8)

cat("\n===== 完成！=====\n")
完整代码（你的 + 添加的）
r
# ===================== 加载数据 =====================
load("mPT_spatial_contact_analysis_results_maligvsothers.RData")
mat_log2fc_mPT <- mat_log2fc_clean
mat_p_mPT <- mat_p

load("nmPT_spatial_contact_analysis_results_maligvsothers.RData")
mat_log2fc_nmPT <- mat_log2fc_clean
mat_p_nmPT <- mat_p

# ===================== 找出代谢型 =====================
cell_types_mPT <- rownames(mat_log2fc_mPT)
cell_types_nmPT <- rownames(mat_log2fc_nmPT)

metabolic1 <- cell_types_mPT[grepl("Metabolic_1", cell_types_mPT)]
metabolic2 <- cell_types_mPT[grepl("Metabolic_2", cell_types_mPT)]

cat("代谢1:", metabolic1, "\n")
cat("代谢2:", metabolic2, "\n")

# ===================== 提取互作 =====================
# mPT
inter_mPT_met1 <- mat_log2fc_mPT[metabolic1, ]
inter_mPT_met2 <- mat_log2fc_mPT[metabolic2, ]

# nmPT
inter_nmPT_met1 <- mat_log2fc_nmPT[metabolic1, ]
inter_nmPT_met2 <- mat_log2fc_nmPT[metabolic2, ]

# 排除自身
inter_mPT_met1 <- inter_mPT_met1[!names(inter_mPT_met1) %in% c(metabolic1, metabolic2)]
inter_mPT_met2 <- inter_mPT_met2[!names(inter_mPT_met2) %in% c(metabolic1, metabolic2)]
inter_nmPT_met1 <- inter_nmPT_met1[!names(inter_nmPT_met1) %in% c(metabolic1, metabolic2)]
inter_nmPT_met2 <- inter_nmPT_met2[!names(inter_nmPT_met2) %in% c(metabolic1, metabolic2)]

# ===================== 找共同细胞类型 =====================
common_cells <- intersect(names(inter_mPT_met1), names(inter_nmPT_met1))

# ===================== 计算差异 =====================
diff_met1 <- inter_mPT_met1[common_cells] - inter_nmPT_met1[common_cells]
diff_met2 <- inter_mPT_met2[common_cells] - inter_nmPT_met2[common_cells]

# 排序
top_met1_mPT <- sort(diff_met1, decreasing = TRUE)[1:15]
top_met2_mPT <- sort(diff_met2, decreasing = TRUE)[1:15]

cat("\n========== 代谢1在 mPT 中更强的互作 ==========\n")
print(top_met1_mPT)

cat("\n========== 代谢2在 mPT 中更强的互作 ==========\n")
print(top_met2_mPT)

# ===================== 代谢2 vs 代谢1 差异分析（基于 mPT vs nmPT 差异）=====================
final_diff <- diff_met2[common_cells] - diff_met1[common_cells]

top_met2_stronger <- sort(final_diff, decreasing = TRUE)[1:10]
top_met1_stronger <- sort(final_diff, decreasing = FALSE)[1:18]

cat("\n========== 代谢2更强的互作（基于转移特异性差异）==========\n")
print(top_met2_stronger)

cat("\n========== 代谢1更强的互作（基于转移特异性差异）==========\n")
print(top_met1_stronger)

# ===================== 条形图 =====================
library(ggplot2)

df_met2 <- data.frame(
  cell_type = names(top_met2_stronger),
  diff = as.numeric(top_met2_stronger),
  stronger = "Metabolic2"
)

df_met1 <- data.frame(
  cell_type = names(top_met1_stronger),
  diff = as.numeric(top_met1_stronger),
  stronger = "Metabolic1"
)

df_plot <- rbind(df_met2, df_met1)

p <- ggplot(df_plot, aes(x = reorder(cell_type, diff), y = diff, fill = stronger)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("Metabolic1" = "#E69F00", "Metabolic2" = "#D55E00")) +
  labs(title = "Functional divergence of metabolic subtypes (based on mPT vs nmPT difference)",
       x = "Cell type", 
       y = "(Metabolic2_mPT - Metabolic2_nmPT) - (Metabolic1_mPT - Metabolic1_nmPT)",
       fill = "Stronger with") +
  theme_minimal() +
  theme(legend.position = "bottom",
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        plot.title = element_text(hjust = 0.5))

ggsave("metabolic_divergence_based_on_transfer_specificity.pdf", plot = p, width = 10, height = 8)

cat("\n===== 完成！=====\n")



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

obj <- readRDS("GSE207422.rds")
obj <- subset(obj, subset = cancer_type == "LUSC")

#obj2 <- readRDS("GSE127465.rds")
#obj3 <- readRDS("GSE153935.rds")
#obj4 <- readRDS("GSE117570_LUSC.rds")
#merged_seurat_obj <- merge(obj1, y = list(obj2, obj3, obj4), 
#                    add.cell.ids = c("GSE207422", "GSE127465", "GSE153935", "GSE117570"))
#merged_seurat_obj <- JoinLayers(merged_seurat_obj)
#library(harmony)
library(Seurat)
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, nfeatures = 2000)

hvgs <- VariableFeatures(obj)
obj <- ScaleData(obj, features = hvgs)
obj <- RunPCA(obj, features = hvgs, npcs = 20)
#obj <- RunHarmony(obj, "orig.ident")

png("elbowplot_sc.png", width = 800, height = 600)
elbowplot <- ElbowPlot(obj,ndims = 30)
print(elbowplot)
dev.off()

library(ggplot2)
library(ggsci)
obj <- FindNeighbors(obj,  dims = 1:20)
obj <- FindClusters(obj, resolution = 0.3)
obj <- RunUMAP(obj, dims = 1:20)

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
    "2" = "Macrophages",
    "3" = "Malignant cells",
    "4" = "Macrophages",
    "5" = "B cells",
    "6" = "Neutrophil",
    "7" = "Plasma",
    "8" = "Proliferative",
    "9" = "DC",
    "10" = "Malignant cells",
    "11" = "Malignant cells",
    "12" = "T cells",
    "13" = "Alveolar type II cells", 
    "14" = "Malignant cells",
    "15" = "Fibroblasts",
    "16" = "Mast cells",
    "17" = "Malignant cells",
    "18" = "Endothelial cells",
    "19" = "Ciliated cells"
)

cell_type <- identity_mapping[obj@meta.data$seurat_clusters]
obj@meta.data$cell_type <- cell_type

library(ggplot2)
library(ggsci)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(13)

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
Malignant_sc <- NormalizeData(Malignant_sc)
Malignant_sc <- FindVariableFeatures(Malignant_sc, nfeatures = 2000)
hvgs <- VariableFeatures(Malignant_sc)
Malignant_sc <- ScaleData(Malignant_sc, features = hvgs)
Malignant_sc <- RunPCA(Malignant_sc, features = hvgs, npcs = 50)
#p <- ElbowPlot(Malignant_sc, ndims = 30)
#ggsave("p3.png",plot=p)

Malignant_sc <- FindNeighbors(Malignant_sc, dims = 1:20)
Malignant_sc <- FindClusters(Malignant_sc, resolution = 0.6)
Malignant_sc <- RunUMAP(Malignant_sc, dims = 1:20)

library(ggplot2)
library(ggsci)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(9)#0.6

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

library(dplyr)
Malignant_markers <- FindAllMarkers(Malignant_sc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "wilcox")
Malignant_significant_markers <- subset(Malignant_markers, p_val_adj < 0.05)
#write.csv(Malignant_significant_markers, "Malignant_all_marker.csv")
Malignant_significant_markers <- Malignant_significant_markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
write.csv(Malignant_significant_markers, "Malignant_sc_top_marker_50_0.6_new.csv")

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

sub_cell_type <- identity_mapping[Malignant_sc@meta.data$seurat_clusters]
Malignant_sc@meta.data$sub_cell_type <- sub_cell_type

library(ggplot2)
library(ggsci)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(9)

cell_types <- as.character(unique(Malignant_sc@meta.data$sub_cell_type))
num_legend_items <- length(cell_types)  # 图例的个数
max_label_length <- max(nchar(cell_types))  # 图例名称的最大长度

base_width <- 3000  # 基础宽度
base_height <- 3000  # 基础高度
legend_width_factor <- 100  # 每个图例项增加的宽度
label_length_factor <- 10  # 每个字符增加的宽度
dynamic_width <- base_width + (num_legend_items * legend_width_factor) + (max_label_length * label_length_factor)
pdf("malignant_sc_annotation-0.6.pdf", width = dynamic_width/300, height = base_height/300)
DimPlot(Malignant_sc, reduction = "umap", label = TRUE, pt.size = 1, group.by = "sub_cell_type", label.size = 4) +
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
saveRDS(Malignant_sc,file="malignant_sc_anno.rds")
```
# malignant mapping
```R
library(Seurat)
library(dplyr)

Malignant_sc <- readRDS("malignant_sc_anno.rds")
Malignant <- readRDS("malignant_unknown_anno_0.6.rds")

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
Malignant <- FindVariableFeatures(Malignant, nfeatures = 2000)
Malignant_sc <- FindVariableFeatures(Malignant_sc, nfeatures = 2000)

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
# ==================== 1. 定义每个亚型的标记基因 ====================
subtype_genes <- list(
  Tumor_EMT_like = c("FN1", "COL1A1", "TAGLN"),
  Tumor_Inflamed_1 = c("CD74", "HLA-DRA"),
  Tumor_Inflamed_2 = c("HLA-E", "CD74"),
  Tumor_Inflamed_3 = c("APOE", "CD68", "CTSD"),
  Tumor_Inflamed_4 = c("CD74", "HLA-DRA"),
  Tumor_Metabolic_1 = c("G6PD", "PGD", "TALDO1"),
  Tumor_Metabolic_2 = c("G6PD", "SPP1", "NAMPT"),
  Tumor_Proliferative_1 = c("TOP2A", "MKI67", "TYMS"),
  Tumor_Proliferative_2 = c("TOP2A", "MKI67", "CCNB1"),
  Tumor_Secretory = c("SCGB1A1", "SCGB3A1", "TFF3"),
  Tumor_Squamous_1 = c("KRT17", "KRT14", "S100A9"),
  Tumor_Squamous_2 = c("KRT16", "KRT6A", "S100A8"),
  Tumor_Squamous_3 = c("KRT16", "S100A9", "SPRR3"),
  Tumor_Stem_like = c("SOX2", "EPCAM", "CDH1")
)

# 按细胞类型首字母排序
subtype_names <- names(subtype_genes)
sorted_subtypes <- subtype_names[order(subtype_names)]

# 按排序后的顺序提取基因
marker_genes <- unique(unlist(subtype_genes[sorted_subtypes]))

# ==================== 2. 检查基因是否存在 ====================
spatial_existing <- intersect(marker_genes, rownames(Malignant))
sc_existing <- intersect(marker_genes, rownames(Malignant_sc))
final_genes <- intersect(spatial_existing, sc_existing)

# 保持基因顺序（按细胞类型首字母排序）
final_genes <- final_genes[order(match(final_genes, marker_genes))]

# ==================== 3. 绘制气泡图 ====================
p1 <- DotPlot(Malignant, features = final_genes, group.by = "sub_cell_type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Spatial Data") +
  xlab("") + ylab("")

p2 <- DotPlot(Malignant_sc, features = final_genes, group.by = "sub_cell_type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "scRNA-seq Data") +
  xlab("") + ylab("")

# 拼图
library(patchwork)
p_combined <- p1 / p2 + plot_annotation(title = "Marker Gene Expression: Spatial vs scRNA-seq")
ggsave("bubble_plot_combined.pdf", plot = p_combined, width = 14, height = 12)
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
fib <- FindVariableFeatures(fib, nfeatures = 2000)
hvgs <- VariableFeatures(fib)
fib <- ScaleData(fib, features = hvgs)
fib <- RunPCA(fib, features = hvgs, npcs = 50)
p <- ElbowPlot(fib, ndims = 30)
ggsave("p3.png",plot=p)

fib <- FindNeighbors(fib, dims = 1:20)
fib <- FindClusters(fib, resolution = 0.3)
fib <- RunUMAP(fib, dims = 1:20)

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


library(dplyr)
fib_markers <- FindAllMarkers(fib, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "wilcox")
fib_significant_markers <- subset(fib_markers, p_val_adj < 0.05)
#write.csv(fib_significant_markers, "fib_all_marker.csv")
fib_significant_markers <- fib_significant_markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
write.csv(fib_significant_markers, "fib_sc_top_marker_50.csv")


identity_mapping <- c(
  "0" = "matCAF",
  "1" = "iCAF",
  "2" = "myCAF",
  "3" = "pCAF"
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
pdf("fib_sc_annotation.pdf", width = dynamic_width/300, height = base_height/300)
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
saveRDS(fib,file="fib_sc_anno.rds")
fib_sc <- fib
```
# fib mapping
```R
library(Seurat)
library(dplyr)

fib_sc <- readRDS("fib_sc_anno.rds")
fib <- readRDS("fib_anno.rds")

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
macro <- subset(obj,subset=cell_type=="Macrophages")
library(Seurat)
macro <- NormalizeData(macro)
macro <- FindVariableFeatures(macro, nfeatures = 2000)
hvgs <- VariableFeatures(macro)
macro <- ScaleData(macro, features = hvgs)
macro <- RunPCA(macro, features = hvgs, npcs = 50)
library(ggplot2)
p <- ElbowPlot(macro, ndims = 30)
ggsave("p3.png",plot=p)

macro <- FindNeighbors(macro, dims = 1:20)
macro <- FindClusters(macro, resolution = 0.3)
macro <- RunUMAP(macro, dims = 1:20)

library(ggplot2)
library(ggsci)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(10)

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


library(dplyr)
macro_markers <- FindAllMarkers(macro, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "wilcox")
macro_significant_markers <- subset(macro_markers, p_val_adj < 0.05)
#write.csv(macro_significant_markers, "macro_all_marker.csv")
macro_significant_markers <- macro_significant_markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
write.csv(macro_significant_markers, "macro_sc_top_marker_50.csv")

identity_mapping <- c(
  "0" = "Lipid-associated_Macro_1",
  "1" = "Inflammatory_Macro",
  "2" = "IFN_Macro",
  "3" = "Tissue-resident_Macro",
  "4" = "Stress-response_Macro",
  "5" = "Lipid-associated_Macro_2",
  "6" = "Lipid-associated_Macro_3",
  "7" = "Proliferative-Macro"
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
pdf("macro_sc_annotation.pdf", width = dynamic_width/300, height = base_height/300)
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
