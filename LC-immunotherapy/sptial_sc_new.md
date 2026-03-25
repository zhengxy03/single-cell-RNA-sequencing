# pretreatmemt
```R
obj <- readRDS("YA2025263-1_fin.rds")
library(Seurat)
library(ggplot2)
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, nfeatures = 2000)
hvgs <- VariableFeatures(obj)
obj <- ScaleData(obj, features = hvgs)
obj <- RunPCA(obj, features = hvgs, npcs = 50)
p <- ElbowPlot(obj, ndims = 50)
ggsave("elbow.png",plot=p)
obj <- RunUMAP(obj, dims = 1:20)
obj <- FindNeighbors(obj, dims = 1:20)
obj <- FindClusters(obj, resolution = 0.3, graph.name = "RNA_snn")
obj_markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1, test.use = "wilcox")
library(dplyr)
obj_significant_markers <- subset(obj_markers, p_val_adj < 0.1)
#write.csv(obj_significant_markers, "obj_all_marker.csv")
obj_significant_markers <- obj_significant_markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
write.csv(obj_significant_markers, "obj_top_marker_50.csv")
saveRDS(obj,file="obj.rds")

library(ggplot2)
library(ggsci)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(26)

seurat_clusters <- as.character(unique(obj@meta.data$seurat_clusters))  # 转换为字符向量
num_legend_items <- length(seurat_clusters)  # 图例的个数
max_label_length <- max(nchar(seurat_clusters))  # 图例名称的最大长度

base_width <- 3000  # 基础宽度
base_height <- 3000  # 基础高度
legend_width_factor <- 100  # 每个图例项增加的宽度
label_length_factor <- 10  # 每个字符增加的宽度
dynamic_width <- base_width + (num_legend_items * legend_width_factor) + (max_label_length * label_length_factor)
pdf("clusters.pdf", width = dynamic_width/300, height = base_height/300)  # 转换为英寸

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

identity_mapping <- c(
    "0" = "B cell",
    "1" = "Macrophage",
    "2" = "T cell",
    "3" = "Fibroblast",
    "4" = "Plasma",
    "5" = "Endothelial cell",
    "6" = "Plasma",
    "7" = "Plasma",
    "8" = "Lymphatic Endothelial cell",
    "9" = "unknown_1",
    "10" = "Malignant cell_1",
    "11" = "Epithelial cell",
    "12" = "Malignant cell_2",
    "13" = "Malignant cell_3", 
    "14" = "B cell",
    "15" = "Macrophage",
    "16" = "unknown_2",
    "17" = "Fibroblast",
    "18" = "unknown_3",
    "19" = "Macrophage",
    "20" = "T cell",
    "21" = "B cell",
    "22" = "B cell",
    "23" = "unknown_4",
    "24" = "Plasma"
)
major_cell_type <- identity_mapping[obj@meta.data$seurat_clusters]
obj@meta.data$major_cell_type <- major_cell_type

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

npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(24)

cell_types <- unique(obj@meta.data$major_cell_type)
num_legend_items <- length(cell_types)  # 图例的个数
max_label_length <- max(nchar(cell_types))  # 图例名称的最大长度

base_width <- 3000  # 基础宽度
base_height <- 3000  # 基础高度
legend_width_factor <- 100  # 每个图例项增加的宽度
label_length_factor <- 10  # 每个字符增加的宽度
dynamic_width <- base_width + (num_legend_items * legend_width_factor) + (max_label_length * label_length_factor)
pdf("annotation.pdf", width = dynamic_width/300, height = base_height/300)
DimPlot(obj, reduction = "umap", label = TRUE, pt.size = 1, group.by = "major_cell_type", label.size = 4) +
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
```
# 统计比例
```R
library(tidyverse)

# 计算各亚群在不同组织中的占比
prop_data <- obj@meta.data %>%
  group_by(tissue, major_cell_type) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(tissue) %>%
  mutate(proportion = count / sum(count) * 100) %>%
  filter(tissue %in% c("mPT", "nmPT")) %>%
  mutate(tissue = factor(tissue, levels = c("mPT", "nmPT")))

library(ggplot2)
library(ggsci)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(15)


prop_data <- prop_data %>%
  mutate(major_cell_type = as.character(major_cell_type))

# 绘制堆叠条形图
p <- ggplot(prop_data, aes(x = tissue, y = proportion, fill = major_cell_type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = npg_extended) +
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
pdf("obj_celltype_distribution_tumor.pdf", width = 4, height = 6)
print(p)
dev.off()

prop_data <- obj@meta.data %>%
  group_by(tissue, major_cell_type) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(tissue) %>%
  mutate(proportion = count / sum(count) * 100) %>%
  filter(tissue %in% c("metLN", "negLN")) %>%
  mutate(tissue = factor(tissue, levels = c("metLN", "negLN")))

library(ggplot2)
library(ggsci)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(15)


prop_data <- prop_data %>%
  mutate(major_cell_type = as.character(major_cell_type))

# 绘制堆叠条形图
p <- ggplot(prop_data, aes(x = tissue, y = proportion, fill = major_cell_type)) +
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
# Malignant
```R
Malignant <- subset(obj,subset=CellType=="Malignant cells")
Malignant <- NormalizeData(Malignant)
Malignant <- FindVariableFeatures(Malignant, nfeatures = 2000)
hvgs <- VariableFeatures(Malignant)
Malignant <- ScaleData(Malignant, features = hvgs)
Malignant <- RunPCA(Malignant, features = hvgs, npcs = 20)
#p <- ElbowPlot(Malignant, ndims = 30)
#ggsave("p3.png",plot=p)
Malignant <- RunUMAP(Malignant, dims = 1:20)
Malignant <- FindNeighbors(Malignant, dims = 1:20)
Malignant <- FindClusters(Malignant, resolution = 0.6)

library(ggplot2)
library(ggsci)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(9)

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
write.csv(Malignant_significant_markers, "Malignant_top_marker_50.csv")

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
sub_cell_type <- identity_mapping[Malignant@meta.data$seurat_clusters]
Malignant@meta.data$sub_cell_type <- sub_cell_type

```