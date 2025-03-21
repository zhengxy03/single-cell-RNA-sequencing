epi <- subset(merged_seurat_obj, subset = seurat_clusters %in% c(2,9,13))
epi <- NormalizeData(epi)
epi <- FindVariableFeatures(epi, nfeatures = 2000)
hvgs <- VariableFeatures(epi)
epi <- ScaleData(epi, features = hvgs)
epi <- RunPCA(epi, features = hvgs, npcs = 20)
epi <- RunHarmony(epi, "sample_sources")
epi <- RunUMAP(epi, dims = 1:20, reduction = "harmony")
epi <- FindNeighbors(epi, dims = 1:20, reduction = "harmony")
epi <- FindClusters(epi, resolution = 0.3)

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

# 假设 identity_mapping 已经定义
cell_type <- identity_mapping[epi@meta.data$seurat_clusters]
epi@meta.data$cell_type <- cell_type