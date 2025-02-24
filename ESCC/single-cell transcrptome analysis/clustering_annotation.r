ElbowPlot(merged_seurat_obj)
merged_seurat_obj <- FindNeighbors(merged_seurat_obj, reduction = "harmony", dims = 1:20)
#seq <- seq(0.3, 1.5, by = 0.1)
#for (res in seq){
    merged_seurat_obj <- FindClusters(merged_seurat_obj, resolution = res)
}

#library(clustree)
#library(patchwork)
#p3 <- clustree(merged_seurat_obj, prefix = 'RNA_snn_res.') + coord_flip()
#p4 <- DimPlot(merged_seurat_obj, group.by = 'RNA_snn_res.0.5', label = T) +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  ggtitle(NULL) +
  guides(color = guide_legend(title = NULL))

library(ggplot2)
library(ggsci)
merged_seurat_obj <- FindClusters(merged_seurat_obj, resolution = 0.5)
merged_seurat_obj <- RunUMAP(merged_seurat_obj, reduction = "harmony", dims = 1:20)

npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(23)

# 获取图例的个数和名称长度
seurat_clusters <- as.character(unique(merged_seurat_obj@meta.data$seurat_clusters))  # 转换为字符向量
num_legend_items <- length(seurat_clusters)  # 图例的个数
max_label_length <- max(nchar(seurat_clusters))  # 图例名称的最大长度

# 动态计算图片尺寸
base_width <- 3000  # 基础宽度
base_height <- 3000  # 基础高度
legend_width_factor <- 100  # 每个图例项增加的宽度
label_length_factor <- 10  # 每个字符增加的宽度

# 计算动态宽度
dynamic_width <- base_width + (num_legend_items * legend_width_factor) + (max_label_length * label_length_factor)

# 导出图片
png("clusters.png", width = dynamic_width, height = base_height, res = 300)
DimPlot(merged_seurat_obj, reduction = "umap", label = TRUE, pt.size = 1, label.size = 8) +
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

#plot on samples
npg_extended <- colorRampPalette(npg_pal)(28)
png("sample.png", width = 3000, height = 3000, res = 300)  # 设置高分辨率和尺寸
DimPlot(merged_seurat_obj, reduction = "umap", label = FALSE, pt.size = 1, group.by = "orig.ident") +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  ggtitle(NULL) +
  scale_color_manual(values = npg_extended) +
  coord_fixed(ratio = 1) +
  guides(color = guide_legend(title = NULL, override.aes = list(size = 5))) +
  theme(
    text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 14, face = "bold", color = "black"),
    axis.title.y = element_text(size = 14, face = "bold", color = "black"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 10, face = "bold", color = "black"),
    legend.title = element_text(size = 10, face = "bold", color = "black"),
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


#annotation
#markers <- FindAllMarkers(object = merged_seurat_obj, 
                                  test.use = "roc", 
                                  only.pos = TRUE, 
                                  min.pct = 0.1, 
                                  thresh.use = 0.25)

#significant_markers <- subset(markers, myAUC > 0.7)
#write.csv(significant_markers,"marker.csv")

markers <- FindAllMarkers(merged_seurat_obj, 
                          only.pos = TRUE, 
                          min.pct = 0.25, 
                          logfc.threshold = 0.25, 
                          test.use = "wilcox")
#significant markers
library(dplyr)
significant_markers <- subset(markers, p_val_adj < 0.05)
significant_markers <- significant_markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.csv(significant_markers,"marker_top.csv")

identity_mapping <- c(
    "0" = "T cell",
    "1" = "B cell",
    "2" = "T cell",
    "3" = "Fibroblast",
    "4" = "Plasma",
    "5" = "Epithelial cell",
    "6" = "Monocyte",
    "7" = "Macrophage",
    "8" = "Endothelial cell",
    "9" = "Dendritic cell",
    "10" = "Fibroblast",
    "11" = "Mast cell",
    "12" = "Pericyte",
    "13" = "Proliferating cell",
    "14" = "T cell",
    "15" = "NK cell",
    "16" = "Epithelial cell",
    "17" = "Fibroblast",
    "18" = "Plasma",
    "19" = "Plasma",
    "20" = "Epithelial cell"
)

cell_type <- identity_mapping[merged_seurat_obj@meta.data$seurat_clusters]
merged_seurat_obj@meta.data$cell_type <- cell_type

npg_extended <- colorRampPalette(npg_pal)(23)

# 获取图例的个数和名称长度
cell_types <- unique(merged_seurat_obj@meta.data$cell_type)
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
png("annotation.png", width = dynamic_width, height = base_height, res = 300)
DimPlot(merged_seurat_obj, reduction = "umap", label = TRUE, pt.size = 1, group.by = "cell_type", label.size = 8) +
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



#DEGs

genes_to_plot <- c(
    # 淋巴细胞系
    "CD3D", "CD3E", "NKG7",          # T细胞
    "BANK1", "CD79A", "MS4A1",        # B细胞
    "GNLY", "KLRC1", "KLRB1",       # NK细胞
    "IGHG1", "IGHG3", "JCHAIN",      # 浆细胞
    
    # 髓系免疫细胞
    "C1QA", "C1QB", "APOC1",         # 巨噬细胞
    "FCN1", "CD300E", "SLC11A1",          # 单核细胞
    "CD86", "PKIB", "CLEC10A",      # 树突细胞（修正HAL-DRA为HLA-DRA）
    
    # 结构细胞
    "SFN", "KRT19", "KRT17",       # 上皮细胞
    "SFRP2", "COL1A2", "FGF7",     # 成纤维细胞
    "CLDN5", "PECAM1", "RAMP2",      # 内皮细胞
    "RGS5", "MCAM", "ACTA2",         # 周细胞
    
    # 特殊功能细胞
    "CPA3", "TPSAB1", "TPSB2",       # 肥大细胞
    "TOP2A", "STMN1", "MKI67"         # 增殖细胞
)
png("dotplot.png", width = 8000, height = 3000, res = 300)  # 设置高分辨率和尺寸

DotPlot(merged_seurat_obj, 
        features = genes_to_plot, 
        group.by = "cell_type",
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

FeaturePlot(merged_seurat_obj, features = genes_to_plot)