# epi

epi <- NormalizeData(epi)
epi <- FindVariableFeatures(epi, nfeatures = 2000)
hvgs <- VariableFeatures(epi)
epi <- ScaleData(epi, features = hvgs)
epi <- RunPCA(epi, features = hvgs, npcs = 20)
library(harmony)
epi <- RunHarmony(epi, "sample_sources")
epi <- RunUMAP(epi, dims = 1:20, reduction = "harmony")
epi <- FindNeighbors(epi, dims = 1:20, reduction = "harmony")
epi <- FindClusters(epi, resolution = 0.1)


seurat_clusters <- as.character(unique(epi@meta.data$seurat_clusters))

num_legend_items <- length(seurat_clusters)

max_label_length <- max(nchar(seurat_clusters))


base_width <- 3000  # 基础宽度
base_height <- 3000  # 基础高度
legend_width_factor <- 100  # 每个图例项增加的宽度
label_length_factor <- 10  # 每个字符增加的宽度


dynamic_width <- base_width + (num_legend_items * legend_width_factor) + (max_label_length * label_length_factor)


npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(13)


pdf("epi_clusters.pdf", width = dynamic_width/300, height = base_height/300)

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


dev.off()

library(dplyr)


epi_markers <- FindAllMarkers(epi, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "wilcox")
epi_significant_markers <- subset(epi_markers, p_val_adj < 0.05)
write.csv(epi_significant_markers, "epi_all_marker.csv")
epi_significant_markers <- epi_significant_markers %>% 
  group_by(cluster) %>% 
  top_n(n = 50, wt = avg_log2FC)
write.csv(epi_significant_markers, "epi_top_marker.csv")


saveRDS(epi, file = "epi_cluster13.rds")


#不同组织中的分布
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(13)


epi_filtered <- subset(epi, subset = sample_type %in% c("Normal", "Tumor", "mLN"))
epi_filtered$sample_type <- factor(epi_filtered$sample_type, 
                                  levels = c("Normal", "Tumor", "mLN"))

pdf("epi_cluster_by_sampletype.pdf", width = 6000/300, height = 3000/300)

p <- DimPlot(epi_filtered, 
             reduction = "umap", 
             label = TRUE, 
             pt.size = 1, 
             group.by = "seurat_clusters", 
             label.size = 8,
             split.by = "sample_type",
             ncol = 3) +  # 改为3列，因为现在只有3个样本类型
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    ggtitle(NULL) +
    scale_color_manual(values = npg_extended) +
    coord_fixed(ratio = 1) +
    guides(color = guide_legend(title = NULL, override.aes = list(size = 5))) +
    theme(
        text = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 20, color = "black"),
        axis.text.y = element_text(size = 20, color = "black"),
        axis.title.x = element_text(size = 24, face = "bold", color = "black"),
        axis.title.y = element_text(size = 24, face = "bold", color = "black", margin = margin(r = 20)),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 20, face = "bold", color = "black"),
        legend.title = element_text(size = 20, face = "bold", color = "black"),
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
        plot.margin = margin(10, 50, 10, 10),
        strip.text = element_text(size = 16, face = "bold", 
                                 margin = margin(b = 15)),
        panel.spacing = unit(1.5, "lines")
    )

print(p)
dev.off()

#饼图
proportion_data <- epi_filtered@meta.data %>%
    group_by(sample_type, seurat_clusters) %>% summarise(count = n()) %>% mutate(proportion = count / sum(count))

npg_extended <- colorRampPalette(npg_pal)(13)
pdf("epi_sampletype_prop.pdf", width = 6000/300, height = 3000/300)  # 设置高分辨率和尺寸
ggplot(proportion_data, aes(x = "", y = proportion, fill = seurat_clusters)) +
  geom_bar(stat = "identity", width = 1) +         # 堆叠柱状图
  coord_polar(theta = "y") +                       # 转换为饼图
  scale_fill_manual(values = npg_extended) +       # 使用自定义颜色
  theme_void() +                                   # 空白背景
  labs(fill = "Cluster") +
  theme(
    legend.position = "right",                     # 图例放在右侧
    plot.title = element_blank(),                  # 移除标题
    # 分面标签设置
    strip.placement = "outside",                   # 标签放在绘图区域外
    strip.text = element_text(                     # 标签样式
      size = 40,                                   # 字体大小
      face = "bold",                               # 加粗
      margin = margin(b = 10)                      # 下方留白（避免与饼图重叠）
    ),
    legend.text = element_text(size = 36),         # 图例文本字体大小
    legend.title = element_text(size = 36)                 # 移除图例标题
  ) +
  facet_wrap(
    ~ sample_type,
    ncol = 3,
    strip.position = "bottom"                      # 标签放在下方
  )
dev.off()

#组织柱状堆叠图
cell_counts <- epi_filtered@meta.data %>%
  group_by(sample_type, seurat_clusters) %>%
  summarise(count = n()) %>%
  ungroup()

npg_extended <- colorRampPalette(npg_pal)(13)
sample_type_order <- c("Tumor", "Normal", "mLN", "nLN")
cell_counts$sample_type <- factor(cell_counts$sample_type, 
                                      levels = sample_type_order)
pdf("epi_sampletype_count.pdf", width = 6000/300, height = 3000/300)  # 设置高分辨率和尺寸


ggplot(cell_counts, aes(x = seurat_clusters, y = count, fill = sample_type)) +
  geom_bar(stat = "identity", position = "stack") +  # 堆叠柱状图
  labs(x = "Cluster", y = "Count", fill = "Sample Type") +  # 设置坐标轴和图例标题
  scale_fill_npg() +  # 使用ggsci的npg配色
  theme_classic() +  # 使用经典主题
  theme(
    axis.text.x = element_text(angle = 0, hjust =1, size = 28),  # 调整横轴标签角度
    axis.text.y = element_text(size = 28),
    axis.title.x = element_blank(),  # 横坐标标题字体大小
    axis.title.y = element_text(size = 40),  # 纵坐标标题字体大小
    legend.text = element_text(size = 36),  # 图例文本字体大小
    legend.title = element_text(size = 40),  # 图例标题字体大小
    legend.position = "right"  # 将图例放在顶部
  )
dev.off()

#病人来源柱状图
patient_distribution <- cell_counts %>%
  group_by(seurat_clusters) %>%
  summarise(
    patient_count = n_distinct(patients),
    patients = paste(sort(unique(patients)), collapse = ", ")
  ) %>%
  arrange(patient_count)

pdf("cluster_patient_distribution.pdf", width = 4000/300, height = 3000/300)
ggplot(patient_distribution, aes(x = reorder(seurat_clusters, -patient_count),  # 添加负号实现从大到小排序
                                y = patient_count, 
                                fill = patient_count)) +
  geom_col() +
  geom_text(aes(label = patient_count), vjust = -0.5, size = 6) +
  scale_fill_viridis_c(name = "Patient Count") +
  labs(x = "Cluster", y = "Number of Patients") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle =0, hjust = 1, size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    plot.title = element_text(size = 28, hjust = 0.5),
    legend.text = element_text(size = 16)
  )
dev.off()
#热图

library(tidyr)

# 创建病人-cluster矩阵
patient_cluster_matrix <- cell_counts %>%
  mutate(present = 1) %>%
  select(patients, seurat_clusters, present) %>%
  pivot_wider(names_from = seurat_clusters, values_from = present, values_fill = 0)

# 转换为矩阵用于热图
matrix_data <- as.matrix(patient_cluster_matrix[, -1])
rownames(matrix_data) <- patient_cluster_matrix$patients

pdf("patient_cluster_heatmap.pdf", width = 5000/300, height = 4000/300)
pheatmap::pheatmap(matrix_data,
                   color = colorRampPalette(c("white", "darkred"))(2),
                   cluster_rows = TRUE,
                   cluster_cols = TRUE,
                   show_rownames = TRUE,
                   show_colnames = TRUE,
                   fontsize = 12,
                   main = "Patient-Cluster Distribution")
dev.off()

#不同分期的堆叠图
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(13)
proportion_data <- epi_filtered@meta.data %>%
    group_by(period1, sample_type, seurat_clusters) %>% 
    summarise(count = n()) %>% 
    mutate(proportion = count / sum(count)) %>%
    # 转换罗马数字为阿拉伯数字
    mutate(period1 = case_when(
        period1 == "Ⅰ" ~ "1",
        period1 == "Ⅱ" ~ "2", 
        period1 == "Ⅲ" ~ "3",
        TRUE ~ as.character(period1)
    )) %>%
    # 创建组合变量
    mutate(period_sample = paste0(period1, "-", sample_type)) %>%
    # 设置组合变量的顺序
    mutate(period_sample = factor(period_sample, 
                                 levels = c("1-Normal", "1-Tumor", 
                                           "2-Normal", "2-Tumor", 
                                           "3-Normal", "3-Tumor", "mLN-mLN")))

pdf("epi_sampletype_period1.pdf", width = 8000/300, height = 3000/300) 

ggplot(proportion_data, aes(x = period_sample, y = proportion, fill = seurat_clusters)) +
  scale_fill_manual(values = npg_extended) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "", y = "Proportion", fill = "Cell Type") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 18, angle = 0, hjust = 1, vjust = 0.5),
    axis.line = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 28),
    legend.title = element_text(size = 28),
    axis.text.y = element_text(size = 18),
    axis.title.y = element_text(size = 20)
  )

dev.off()
saveRDS(epi_filtered, file = "epi_filtered.rds")

#轨迹
library(Seurat)
library(monocle)
library(ggplot2)
library(ggsci)
trace('project2MST', edit = T, where = asNamespace("monocle"))

expression_matrix <- LayerData(epi, assay = "RNA", layer = "data")
cell_metadata <- epi@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(expression_matrix))
rownames(gene_annotation) <- rownames(expression_matrix)

# 随机选择部分细胞
set.seed(123)
subset_cells <- sample(colnames(expression_matrix), size = 20000)  # 选择 20,000 个细胞
expression_matrix <- expression_matrix[, subset_cells]
cell_metadata <- cell_metadata[subset_cells, ]
cds <- newCellDataSet(expression_matrix,
                      phenoData = new("AnnotatedDataFrame", data = cell_metadata),
                      featureData = new("AnnotatedDataFrame", data = gene_annotation),
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds, min_expr = 0.1)
disp_table <- dispersionTable(cds)
ordering_genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
cds <- setOrderingFilter(cds, ordering_genes)

cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')
cds <- orderCells(cds)
saveRDS(cds, file = "epi_cds.rds")
cds <- orderCells(cds, root_state  = 4)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(13)
p <- plot_cell_trajectory(cds, color_by = "State") +
  facet_wrap("~State", nrow = 2) + scale_color_manual(values = npg_extended)
ggsave("epi_traj_state.png", plot = p)

p <- plot_cell_trajectory(cds, color_by = "seurat_clusters") +
  facet_wrap("~seurat_clusters", nrow = 3) + scale_color_manual(values = npg_extended)

ggsave("epi_traj_cluster.png", plot = p)


#基因表达
genes_to_plot <- c(
    "UBE2C", "UBE2T","BIRC5","PCLAF","HMGB3",
    "GLI2","ALCAM","NEAT1","EDN1",
    "LAMB3","KRT7","LMO7","CEACAM1",
    "SPRR1B","PRR1A","SPRR2A","SPRR2D","SPRR3",
    "IDO1","CXCL10","GBP5","CDK1"

)

pdf("dotplot.pdf", width = 8000/300, height = 3000/300)  # 设置高分辨率和尺寸
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