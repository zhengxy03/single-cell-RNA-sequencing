epi <- subset(merged_seurat_obj, subset = seurat_clusters %in% c(2,9,13))
epi <- NormalizeData(epi)
epi <- FindVariableFeatures(epi, nfeatures = 2000)
hvgs <- VariableFeatures(epi)
epi <- ScaleData(epi, features = hvgs)
epi <- RunPCA(epi, features = hvgs, npcs = 20)
epi <- RunHarmony(epi, "sample_sources")
epi <- RunUMAP(epi, dims = 1:20, reduction = "harmony")
epi <- FindNeighbors(epi, dims = 1:20, reduction = "harmony")
epi <- FindClusters(epi, resolution = 0.1)

# 获取唯一的聚类标签并转换为字符向量
seurat_clusters <- as.character(unique(epi@meta.data$seurat_clusters))
# 计算图例的个数
num_legend_items <- length(seurat_clusters)
# 计算图例名称的最大长度
max_label_length <- max(nchar(seurat_clusters))

# 动态计算图片尺寸
base_width <- 3000  # 基础宽度
base_height <- 3000  # 基础高度
legend_width_factor <- 100  # 每个图例项增加的宽度
label_length_factor <- 10  # 每个字符增加的宽度

# 计算动态宽度
dynamic_width <- base_width + (num_legend_items * legend_width_factor) + (max_label_length * label_length_factor)

# 创建 NPG 颜色调色板
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(17)

# 打开 PNG 图形设备，设置图片尺寸和分辨率
png("epi_clusters.png", width = dynamic_width, height = base_height, res = 300)

# 绘制 UMAP 图
DimPlot(epi, reduction = "umap", label = TRUE, pt.size = 1, label.size = 8) +
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
        axis.title.y = element_text(size = 36, face = "bold", color = "black", margin = margin(r = 20)),  # 增加右侧间距
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

# 关闭图形设备
dev.off()


library(dplyr)


epi_markers <- FindAllMarkers(epi, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "wilcox")
epi_significant_markers <- subset(epi_markers, p_val_adj < 0.05)
write.csv(epi_significant_markers, "epi_all_marker.csv")
epi_significant_markers <- epi_significant_markers %>% 
  group_by(cluster) %>% 
  top_n(n = 50, wt = avg_log2FC)
write.csv(epi_significant_markers, "epi_top_marker.csv")

Proliferative 
identity_mapping <- c(
    "0" = "Proliferative tumor BC",
    "1" = "tumor differentiated SC",
    "2" = "normal differentiated GC",
    "3" = "tumor differentiated SC",
    "4" = "Proliferative tumor Basal cells",
    "5" = "tumor",
    "6" = "tumor differentiated GC",
    "7" = "invasive tumor KC",
    "8" = "tumor",
    "9" = "tumor differentiated SC",
    "10" = "tumor differentiated GC",
    "11" = "tumor differentiated",
    "12" = "免疫细胞",
    "13" = "normal BC",
    "14" = "EMT-like Epi",
    "15" = "分泌型腺上皮细胞 ",
    "16" = "神经内分泌细胞" 
)

identity_mapping <- c(
    "0" = "Proliferative cells",
    "1" = "Invasive cells",
    "2" = "Differentiated cells",
    "3" = "Invasive cells",
    "4" = "Invasive cells",
    "5" = "Invasive cells",
    "6" = "Invasive cells",
    "7" = "Keratinocytes",
    "8" = "Differentiated cells",
    "9" = "Urothelial cells",
    "10" = "Basal cells",
    "11" = "Proliferative cells",
    "12" = "Immune-associated invasive cells",
    "13" = "Hepatocytes",
    "14" = "EMT-like Epis",
    "15" = "Invasive cells",
    "16" = "Neuroendocrine cells" 
)

identity_mapping <- c(
    "0" = "Proliferative cells",
    "1" = "Invasive cells",
    "2" = "Superficial cells",
    "3" = "Immune-associated Epis",
    "4" = "Invasive cells",
    "5" = "Secretory cells",
    "6" = "Invasive cells",
    "7" = "Superficial cells",
    "8" = "Parabasal cells",
    "9" = "Parabasal cells",
    "10" = "Parabasal cells",
    "11" = "Proliferative cells",
    "12" = "Immune-associated Epis",
    "13" = "Invasive cells",
    "14" = "Invasive cells",
    "15" = "Invasive cells",
    "16" = "Neuroendocrine cells" 
)

# 假设 identity_mapping 已经定义
cell_type <- identity_mapping[epi@meta.data$seurat_clusters]
epi@meta.data$cell_type <- cell_type
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(8)

# 获取图例的个数和名称长度
cell_types <- unique(epi@meta.data$cell_type)
num_legend_items <- length(cell_types)  # 图例的个数
max_label_length <- max(nchar(cell_types))  # 图例名称的最大长度

# 动态计算图片尺寸
base_width <- 3000  # 基础宽度
base_height <- 3000  # 基础高度
legend_width_factor <- 100  # 每个图例项增加的宽度
label_length_factor <- 10  # 每个字符增加的宽度

# 计算动态宽度
dynamic_width <- base_width + (num_legend_items * legend_width_factor) + (max_label_length * label_length_factor)

# 导出图片
png("epi_annotation.png", width = dynamic_width, height = base_height, res = 300)
DimPlot(epi, reduction = "umap", label = FALSE, pt.size = 1, group.by = "cell_type", label.size = 8) +
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
        axis.title.y = element_text(size = 36, face = "bold", color = "black", margin = margin(r = 20)),  # 增加右侧间距
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 16, face = "bold", color = "black"),
        legend.title = element_text(size = 28, face = "bold", color = "black"),
        legend.position = "bottom",
        legend.justification = "center",
        legend.box.just = "center",
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

genes_to_plot <- c(

    "COL17A1", "NOTCH1","KRT5", "KRT14", "TP63", "ITGA6", "NGFR", "SOX2",

    "JAG1", "MKI67","PCNA", "TOP2A", "CCNB1", "CDK1",
    "ANXA1", "ECM1","KRT1", "KRT10", "IVL", "DSP", "FLG",
    "NFE2L2", "SPP1", "MMP11","VIM", "FN1", "SNAI1", "SNAI2", "ZEB1", "ZEB2", "CDH2", "MMP2", "MMP9"
)

png("epi_dotplot.png", width = 8000, height = 3000, res = 300)  # 设置高分辨率和尺寸
DotPlot(epi, 
        features = genes_to_plot, 
        group.by = "seurat_clusters",
        dot.scale = 15) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 28),  # 横坐标基因名字体大小
    axis.text.y = element_text(size = 28),  # 纵坐标细胞类型字体大小
    axis.title.x = element_blank(),  # 去除横坐标标题
    axis.title.y = element_text(size = 40),
    panel.grid.major = element_line(color = "grey70", linewidth = 0.5),  # 添加网格线
    panel.grid.minor = element_line(color = "grey80", linewidth = 0.3),  # 细网格线 
    legend.text = element_text(size = 28, face = "bold", color = "black"),
    legend.title = element_text(size = 28, face = "bold", color = "black"),
    legend.position = "right",
    legend.box.margin = margin(0, 0, 0, 0),
    legend.key = element_blank(),
    legend.background = element_blank()
  ) +
  scale_color_gsea()  # 使用ggsci的GSEA配色
dev.off()

identity_mapping <- c(
    "0" = "Proliferative cells_1",
    "1" = "Invasive cells_1",
    "2" = "Superficial cells_1",
    "3" = "Immune-associated Epis_1",
    "4" = "Invasive cells_2",
    "5" = "Secretory cells",
    "6" = "Invasive cells_3",
    "7" = "Superficial cells_2",
    "8" = "Parabasal cells_1",
    "9" = "Parabasal cells_2",
    "10" = "Parabasal cells_3",
    "11" = "Proliferative cells_2",
    "12" = "Immune-associated Epis_2",
    "13" = "Invasive cells_4",
    "14" = "Invasive cells_5",
    "15" = "Invasive cells_6",
    "16" = "Neuroendocrine cells" 
)

cell_type2 <- identity_mapping[epi@meta.data$seurat_clusters]
epi@meta.data$cell_type2 <- cell_type2
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(17)
epi@meta.data$cell_type2 <- factor(epi@meta.data$cell_type2, levels = identity_mapping)

# 获取图例的个数和名称长度
cell_types <- as.character(unique(epi@meta.data$cell_type2))
num_legend_items <- length(cell_types)  # 图例的个数
max_label_length <- max(nchar(cell_types))  # 图例名称的最大长度

# 动态计算图片尺寸
base_width <- 3000  # 基础宽度
base_height <- 3000  # 基础高度
legend_width_factor <- 100  # 每个图例项增加的宽度
label_length_factor <- 10  # 每个字符增加的宽度

# 计算动态宽度
dynamic_width <- base_width + (num_legend_items * legend_width_factor) + (max_label_length * label_length_factor)

# 导出图片
png("epi_annotation2.png", width = dynamic_width, height = base_height, res = 300)
DimPlot(epi, reduction = "umap", label = FALSE, pt.size = 1, group.by = "cell_type2", label.size = 8) +
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
        axis.title.y = element_text(size = 36, face = "bold", color = "black", margin = margin(r = 20)),  # 增加右侧间距
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 16, face = "bold", color = "black"),
        legend.title = element_text(size = 28, face = "bold", color = "black"),
        legend.position = "bottom",
        legend.justification = "center",
        legend.box.just = "center",
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