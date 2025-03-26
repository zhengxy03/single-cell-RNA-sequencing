epi_fib <- NormalizeData(epi_fib)
epi_fib <- FindVariableFeatures(epi_fib, nfeatures = 2000)
hvgs <- VariableFeatures(epi_fib)
epi_fib <- ScaleData(epi_fib, features = hvgs)
epi_fib <- RunPCA(epi_fib, features = hvgs, npcs = 20)
epi_fib <- RunHarmony(epi_fib, "sample_sources")
epi_fib <- RunUMAP(epi_fib, dims = 1:20, reduction = "harmony")
epi_fib <- FindNeighbors(epi_fib, dims = 1:20, reduction = "harmony")
epi_fib <- FindClusters(epi_fib, resolution = 0.3)

seurat_clusters <- as.character(unique(epi_fib@meta.data$seurat_clusters))
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
npg_extended <- colorRampPalette(npg_pal)(24)

# 打开 PNG 图形设备，设置图片尺寸和分辨率
png("epi_fib_clusters.png", width = dynamic_width, height = base_height, res = 300)

# 绘制 UMAP 图
DimPlot(epi_fib, reduction = "umap", label = TRUE, pt.size = 1, label.size = 8) +
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
epi_fib_markers <- FindAllMarkers(epi_fib, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "wilcox")
epi_fib_significant_markers <- subset(epi_fib_markers, p_val_adj < 0.05)
write.csv(epi_fib_significant_markers, "epi_fib_all_marker.csv")
epi_fib_significant_markers <- epi_fib_significant_markers %>% 
  group_by(cluster) %>% 
  top_n(n = 50, wt = avg_log2FC)
write.csv(epi_fib_significant_markers, "epi_fib_top_marker.csv")



identity_mapping <- c(
    "0" = "COL11A1+ myCAFs",
    "1" = "iCAF/myCAF transition",
    "2" = "IGF1+ iCAFs",
    "3" = "Proliferative cells",
    "4" = "Basal cells",
    "5" = "PI16+ Fib progenitors",
    "6" = "Differentiated cells",
    "7" = "Proliferative cells",
    "8" = "Keratinocytes",
    "9" = "Proliferative cells",
    "10" = "vCAFs",
    "11" = "Invasive cells",
    "12" = "Invasive cells",
    "13" = "Differentiated cells",
    "14" = "Proliferative cells",
    "15" = "Invasive cells",
    "16" = "COL7A1+ myCAFs",
    "17" = "Immune-associated NAFs",
    "18" = "Neuroendocrine cells",
    "19" = "Urothelial cells",
    "20" = "Hepatocytes",
    "21" = "Proliferative cells",
    "22" = "Proliferative myCAF",
    "23" = "EMT-like Epi"
)
cell_type <- identity_mapping[epi_fib@meta.data$seurat_clusters]
epi_fib@meta.data$cell_type <- cell_type
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(24)

# 获取图例的个数和名称长度
cell_types <- unique(epi_fib@meta.data$cell_type)
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
png("epi_fib_annotation.png", width = dynamic_width, height = base_height, res = 300)
DimPlot(epi_fib, reduction = "umap", label = FALSE, pt.size = 1, group.by = "cell_type", label.size = 8) +
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

png("epi_fib_dotplot.png", width = 8000, height = 3000, res = 300)  # 设置高分辨率和尺寸
DotPlot(epi_fib, 
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
    "0" = "COL11A1+ myCAFs",
    "1" = "iCAF/myCAF transition",
    "2" = "IGF1+ iCAFs",
    "3" = "Proliferative Epis_1",
    "4" = "Basal Epis",
    "5" = "PI16+ Fib progenitors",
    "6" = "Differentiated Epis_1",
    "7" = "Proliferative Epis_2",
    "8" = "Keratinocytes",
    "9" = "Proliferative Epis_3",
    "10" = "vCAFs",
    "11" = "Invasive Epis_1",
    "12" = "Invasive Epis_2",
    "13" = "Differentiated Epis_2",
    "14" = "Proliferative Epis_4",
    "15" = "Invasive Epis_3",
    "16" = "COL7A1+ myCAFs",
    "17" = "Immune-associated NAFs",
    "18" = "Neuroendocrine cells",
    "19" = "Urothelial cells",
    "20" = "Hepatocytes",
    "21" = "Proliferative Epi_5",
    "22" = "Proliferative myCAF",
    "23" = "EMT-like Epi"
)
cell_type <- identity_mapping[epi_fib@meta.data$seurat_clusters]
epi_fib@meta.data$cell_type <- cell_type
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(24)

# 获取图例的个数和名称长度
cell_types <- unique(epi_fib@meta.data$cell_type)
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
png("epi_fib_annotation.png", width = dynamic_width, height = base_height, res = 300)
DimPlot(epi_fib, reduction = "umap", label = FALSE, pt.size = 1, group.by = "cell_type", label.size = 8) +
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