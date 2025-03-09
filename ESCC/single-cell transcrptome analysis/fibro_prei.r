fibro_peri <- subset(merged_seurat_obj, subset = seurat_clusters %in% c(5,7,11))
fibro_peri <- NormalizeData(fibro_peri)
fibro_peri <- FindVariableFeatures(fibro_peri, nfeatures = 2000)
hvgs <- VariableFeatures(fibro_peri)
fibro_peri <- ScaleData(fibro_peri, features = hvgs)
fibro_peri <- RunPCA(fibro_peri, features = hvgs, npcs = 20)

fibro_peri <- RunHarmony(fibro_peri, "sample_sources")
fibro_peri <- RunUMAP(fibro_peri, dims = 1:20, reduction = "harmony")
fibro_peri <- FindNeighbors(fibro_peri, dims = 1:20, reduction = "harmony")
fibro_peri <- FindClusters(fibro_peri, resolution = 0.4)

seurat_clusters <- as.character(unique(fibro_peri@meta.data$seurat_clusters))  # 转换为字符向量
num_legend_items <- length(seurat_clusters)  # 图例的个数
max_label_length <- max(nchar(seurat_clusters))  # 图例名称的最大长度

# 动态计算图片尺寸
base_width <- 3000  # 基础宽度
base_height <- 3000  # 基础高度
legend_width_factor <- 100  # 每个图例项增加的宽度
label_length_factor <- 10  # 每个字符增加的宽度

# 计算动态宽度
dynamic_width <- base_width + (num_legend_items * legend_width_factor) + (max_label_length * label_length_factor)

npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(18)
png("fibro_peri_clusters.png", width = dynamic_width, height = base_height, res = 300)
DimPlot(fibro_peri, reduction = "umap", label = TRUE, pt.size = 1, label.size = 8) +
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
dev.off()
