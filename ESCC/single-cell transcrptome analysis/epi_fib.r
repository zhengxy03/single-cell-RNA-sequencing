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
    "3" = "Proliferative Basal cells",
    "4" = "Basal cells",
    "5" = "PI16+ Fib progenitors",
    "6" = "Differentiated cells",
    "7" = "Proliferative cells",
    "8" = "Keratinocytes",
    "9" = "Proliferative cells(zand)",
    "10" = "Basal cells_1",
    "11" = "Proliferative cells_2",
    "12" = "Immune-associated invasive cells",
    "13" = "Basal cells_2",
    "14" = "EMT-like Epi",
    "15" = "Invasive cells_6",
    "16" = "Differentiated cells_5" 
)
