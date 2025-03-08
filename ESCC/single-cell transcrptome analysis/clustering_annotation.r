png("elbowplot.png", width = 800, height = 600)
elbowplot <- ElbowPlot(merged_seurat_obj)
print(elbowplot)
dev.off()

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
npg_extended <- colorRampPalette(npg_pal)(15)

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
merged_seurat_obj <- FindClusters(merged_seurat_obj, resolution = 0.3)
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
    "1" = "T cell",
    "2" = "Epithelial cell",
    "3" = "B cell",
    "4" = "Myeloid cell",
    "5" = "Fibroblast",
    "6" = "Endothelial cell",
    "7" = "Fibroblast",
    "8" = "Plasma",
    "9" = "Epithelial cell",
    "10" = "Mast cell",
    "11" = "Pericyte",
    "12" = "Proliferating cell",
    "13" = "Tumor cell",
    "14" = "B cell"
)

cell_type <- identity_mapping[merged_seurat_obj@meta.data$seurat_clusters]
merged_seurat_obj@meta.data$cell_type <- cell_type

npg_extended <- colorRampPalette(npg_pal)(15)

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

#sample_sources
unique_sources <- sort(unique(merged_seurat_obj@meta.data$sample_sources))  # 按字母顺序排序
# 如果需要手动指定顺序，可以这样做：
# unique_periods <- c("period_A", "period_B", "period_C", "period_D", "period_E", "period_F")

# 创建一个空列表，用于存储每个 period 的图
plot_list <- list()

# 遍历每个 period，生成单独的图
for (source in unique_sources) {
    # 创建颜色映射：当前 period 为 npg 红色，其他为灰色
    color_mapping <- setNames(
        ifelse(unique_sources == source, pal_npg()(1), "gray"),  # npg 红色
        unique_sources
    )
    
    # 绘制 DimPlot
    p <- DimPlot(merged_seurat_obj, reduction = "umap", label = FALSE, pt.size = 1, group.by = "sample_sources") +
        scale_color_manual(values = color_mapping) +  # 使用自定义颜色映射
        ggtitle(source) +  # 设置标题为当前 period
        theme(
            legend.position = "none",  # 隐藏图例
            plot.title = element_text(size = 24, face = "bold", hjust = 0.5),  # 标题样式
            axis.title.x = element_text(size = 20, face = "bold", color = "black"),  # X 轴标题
            axis.title.y = element_text(size = 20, face = "bold", color = "black"),  # Y 轴标题
            axis.text.x = element_text(size = 16, color = "black"),  # X 轴刻度
            axis.text.y = element_text(size = 16, color = "black")   # Y 轴刻度
        )
    
    # 将图添加到列表中
    plot_list[[source]] <- p
}

# 使用 patchwork 将图形排列在一起（每行三个）
combined_plot <- wrap_plots(plot_list, ncol = 2)  # 每行三个图
png("combined_source_umap.png", width = 6000, height =3000, res = 300)
print(combined_plot)
dev.off()

#DEGs

genes_to_plot <- c(
    # 淋巴细胞系
    "CD3D", "CD3E", "NKG7",          # T细胞
    "BANK1", "CD79A", "MS4A1",        # B细胞

    "IGHG1", "IGHG3", "JCHAIN",      # 浆细胞
    
    # 髓系免疫细胞
    "C1QA", "C1QB", "APOC1",         # 巨噬细胞

    "CD86", "LYZ", "PKIB",      # 树突细胞（修正HAL-DRA为HLA-DRA）
    
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