merged_seurat_obj <- readRDS("merged_anno_period.rds")
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggsci)
t_cells <- subset(merged_seurat_obj, subset = cell_type == "T cell")
t_cells <- NormalizeData(t_cells)
t_cells <- FindVariableFeatures(t_cells, nfeatures = 2000)
hvgs <- VariableFeatures(t_cells)
t_cells <- ScaleData(t_cells, features = hvgs)
t_cells <- RunPCA(t_cells, features = hvgs, npcs = 20)
library(harmony)
t_cells <- RunHarmony(t_cells, "sample_sources")
t_cells <- RunUMAP(t_cells, dims = 1:20, reduction = "harmony")
t_cells <- FindNeighbors(t_cells, dims = 1:20, reduction = "harmony")
t_cells <- FindClusters(t_cells, resolution = 0.3)
t_cells <- subset(t_cells, subset = seurat_clusters %in% c(0,1,2,3,4,5,6,7))
t_cells <- NormalizeData(t_cells)
t_cells <- FindVariableFeatures(t_cells, nfeatures = 2000)
hvgs <- VariableFeatures(t_cells)
t_cells <- ScaleData(t_cells, features = hvgs)
t_cells <- RunPCA(t_cells, features = hvgs, npcs = 20)
library(harmony)
t_cells <- RunHarmony(t_cells, "sample_sources")
t_cells <- RunUMAP(t_cells, dims = 1:20, reduction = "harmony")
t_cells <- FindNeighbors(t_cells, dims = 1:20, reduction = "harmony")
t_cells <- FindClusters(t_cells, resolution = 0.3)


# 获取图例的个数和名称长度
seurat_clusters <- as.character(unique(t_cells@meta.data$seurat_clusters))  # 转换为字符向量
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
npg_extended <- colorRampPalette(npg_pal)(9)
pdf("t_clusters.pdf", width = dynamic_width/300, height = base_height/300)
DimPlot(t_cells, reduction = "umap", label = TRUE, pt.size = 1, label.size = 8) +
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



#t_cell_markers <- FindAllMarkers(t_cells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "wilcox")
#t_significant_markers <- subset(t_cell_markers, p_val_adj < 0.05)
#write.csv(t_significant_markers, "t_all_marker.csv")
#t_significant_markers <- t_significant_markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
#write.csv(t_significant_markers, "t_top_marker.csv")

identity_mapping <- c(
  "0" = "GZMK+ CD8+Tem",
  "1" = "CCR7+ Tn",
  "2" = "FOXP3+ Treg", 
  "3" = "HAVCR2+ CD8+Tex",
  "4" = "Activited T",
  "5" = "FCER1G+ TNK",
  "6" = "TOX+ CD4+Tex",
  "7" = "FCGR3A+ TNK",
  "8" = "Stress Response T"
)



cell_type <- identity_mapping[t_cells@meta.data$seurat_clusters]
t_cells@meta.data$cell_type <- cell_type
saveRDS(t_cells,file="t_cells_anno.rds")
# 将 cell_type 列转换为因子，并指定顺序
#t_cells@meta.data$cell_type <- factor(t_cells@meta.data$cell_type, levels = identity_mapping)


npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(9)
# 获取唯一的细胞类型并转换为字符向量
cell_types <- as.character(unique(t_cells@meta.data$cell_type))
num_legend_items <- length(cell_types)  # 图例的个数
max_label_length <- max(nchar(cell_types))  # 图例名称的最大长度

# 动态计算图片尺寸
base_width <- 5000  # 基础宽度
base_height <- 5000  # 基础高度
legend_width_factor <- 100  # 每个图例项增加的宽度
label_length_factor <- 10  # 每个字符增加的宽度

# 计算动态宽度
dynamic_width <- base_width + (num_legend_items * legend_width_factor) + (max_label_length * label_length_factor)

# 保存图片
pdf("t_annotation.pdf", width = dynamic_width/300, height = base_height/300)
DimPlot(t_cells, reduction = "umap", label = TRUE, pt.size = 1, label.size = 8, group.by = "cell_type") +
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    ggtitle(NULL) +
    scale_color_manual(values = npg_extended) +
    coord_fixed(ratio = 1) +
    guides(color = guide_legend(title = NULL, override.aes = list(size = 5))) +
    theme(
        text = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size =40, color = "black"),
        axis.text.y = element_text(size = 40, color = "black"),
        axis.title.x = element_text(size = 56, face = "bold", color = "black"),
        axis.title.y = element_text(size = 56, face = "bold", color = "black", margin = margin(r = 20)),  # 增加右侧间距
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 40, face = "bold", color = "black"),
        legend.title = element_text(size = 40, face = "bold", color = "black"),
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
    ) +scale_x_continuous(limits = c(-10.5,10))
dev.off()

#比例柱状图
# 计算每个sample_type和period1组合中的cell_type比例
cell_proportions <- t_cells@meta.data %>%
  group_by(sample_type, period1, cell_type) %>%
  summarise(count = n()) %>%
  group_by(sample_type, period1) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()

# 转换罗马数字为阿拉伯数字
cell_proportions <- cell_proportions %>%
  mutate(period1 = case_when(
    period1 == "Ⅰ" ~ "1",
    period1 == "Ⅱ" ~ "2", 
    period1 == "Ⅲ" ~ "3",
    TRUE ~ as.character(period1)
  ))

# 创建组合标签并设置正确的顺序
cell_proportions <- cell_proportions %>%
  mutate(
    sample_period = case_when(
      sample_type == "Tumor" ~ paste0("Tumor_", period1),
      sample_type == "mLN" ~ "mLN",
      sample_type == "Normal" ~ paste0("Normal_", period1),
      sample_type == "nLN" ~ "nLN",
      TRUE ~ paste0(sample_type, "_", period1)
    )
  )

# 设置样本顺序
sample_period_order <- c(
  "Tumor_1", "Tumor_2", "Tumor_3",
  "mLN",
  "Normal_1", "Normal_2", "Normal_3", 
  "nLN"
)

cell_proportions$sample_period <- factor(cell_proportions$sample_period, 
                                        levels = sample_period_order)

# 设置颜色
npg_pal <- pal_npg()(10)
cell_type_colors <- colorRampPalette(npg_pal)(9)

# 绘制分面柱状图
p <- ggplot(cell_proportions, aes(x = cell_type, y = proportion, fill = cell_type)) +
  geom_col() +
  scale_fill_manual(values = cell_type_colors, name = "Cell Type") +
  labs(x = "", y = "Proportion") +
  theme_classic() +
  theme(
    # 横坐标设置 - 确保所有子图都有轴线
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.5),  # 确保横坐标轴线显示
    
    # 纵坐标设置
    axis.text.y = element_text(size = 9),
    axis.title.y = element_text(size = 11),
    
    # 图例位置 - 放在整个图形的顶部
    legend.position = "top",
    legend.direction = "horizontal",
    legend.justification = "center",
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10),
    legend.box = "horizontal",
    
    # 分面标题设置 - 移除方框，放在横坐标下方
    strip.background = element_blank(),  # 移除标题背景
    strip.text = element_text(size = 10, face = "bold", vjust = 1),  # 标题文字
    strip.placement = "outside",  # 标题放在绘图区域外
    
    # 其他设置
    panel.spacing = unit(0.8, "lines")
  ) +
  facet_wrap(~ sample_period, ncol = 4, scales = "free_x")  # 每行4张子图

# 显示图形
print(p)

# 保存图形
ggsave("t_cells_proportion_sorted.pdf", p, 
       width = 16, height = 10, dpi = 300)


#比例气泡图
# 计算每个sample_type和period1组合中的cell_type比例
cell_proportions <- t_cells@meta.data %>%
  group_by(sample_type, period1, cell_type) %>%
  summarise(count = n()) %>%
  group_by(sample_type, period1) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()

# 转换罗马数字为阿拉伯数字
cell_proportions <- cell_proportions %>%
  mutate(period1 = case_when(
    period1 == "Ⅰ" ~ "1",
    period1 == "Ⅱ" ~ "2", 
    period1 == "Ⅲ" ~ "3",
    TRUE ~ as.character(period1)
  ))

# 创建组合标签并设置正确的顺序
cell_proportions <- cell_proportions %>%
  mutate(
    sample_period = case_when(
      sample_type == "Tumor" ~ paste0("Tumor_", period1),
      sample_type == "mLN" ~ "mLN",
      sample_type == "Normal" ~ paste0("Normal_", period1),
      sample_type == "nLN" ~ "nLN",
      TRUE ~ paste0(sample_type, "_", period1)
    )
  )

# 设置样本顺序
sample_period_order <- c(
  "Tumor_1", "Tumor_2", "Tumor_3",
  "mLN",
  "Normal_1", "Normal_2", "Normal_3", 
  "nLN"
)

cell_proportions$sample_period <- factor(cell_proportions$sample_period, 
                                        levels = sample_period_order)

# 设置颜色
npg_pal <- pal_npg()(10)
cell_type_colors <- colorRampPalette(npg_pal)(length(unique(cell_proportions$cell_type)))

# 绘制气泡图
p <- ggplot(cell_proportions, aes(x = cell_type, y = sample_period, 
                                 size = proportion, color = cell_type)) +
  geom_point(alpha = 0.8) +  # 使用geom_point绘制气泡
  scale_size_continuous(name = "Proportion", 
                       range = c(1, 10),  # 调整气泡大小范围
                       limits = c(0, 1),
                       breaks = c(0.2, 0.4, 0.6, 0.8)) +  # 设置图例刻度
  scale_color_manual(values = cell_type_colors) +
  scale_y_discrete(limits = rev(levels(cell_proportions$sample_period))) +  # 反转y轴顺序
  labs(x = "Cell Type", y = "Sample Type & Period", 
       title = "T Cell Proportions Across Samples") +
  theme_classic() +  # 使用theme_classic，它自带坐标轴线
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
    axis.text.y = element_text(size = 11, face = "bold"),
    axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 12, face = "bold", margin = margin(r = 10)),
    legend.position = "right",  # 显示图例在右侧
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    # 移除网格线
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # 确保坐标轴框线显示
    axis.line = element_line(color = "black", linewidth = 0.5),
    panel.background = element_rect(fill = "white", color = "black", linewidth = 0.5),  # 添加面板边框
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  ) +
  guides(color = "none")  # 移除颜色图例，只保留大小图例

# 显示图形
print(p)

# 保存图形
ggsave("t_cells_bubble_plot.pdf", p, width = 16, height = 10, dpi = 300)















#不同组织中的分布
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(9)

t_cells$sample_type <- factor(t_cells$sample_type, 
                              levels = c("Normal", "Tumor", "mLN","nLN"))

pdf("t_cells_cluster_by_sampletype.pdf", width = 6000/300, height = 3000/300)

p <- DimPlot(t_cells, 
             reduction = "umap", 
             label = FALSE, 
             pt.size = 1, 
             group.by = "cell_type", 
             label.size = 8,
             split.by = "sample_type",
             ncol = 2) + 
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


#轨迹
library(Seurat)
library(monocle)
library(ggplot2)
library(ggsci)
trace('project2MST', edit = T, where = asNamespace("monocle"))

expression_matrix <- LayerData(t_cells, assay = "RNA", layer = "data")
cell_metadata <- t_cells@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(expression_matrix))
rownames(gene_annotation) <- rownames(expression_matrix)

# 随机选择部分细胞
set.seed(123)
subset_cells <- sample(colnames(expression_matrix), size = 10000)  # 选择 20,000 个细胞
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
saveRDS(cds, file = "t_cells_cds.rds")
cds <- orderCells(cds, root_state  = 4)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(9)
p <- plot_cell_trajectory(cds, color_by = "State") +
  facet_wrap("~State", nrow = 2) + scale_color_manual(values = npg_extended)
ggsave("t_cells_traj_state.png", plot = p)

p <- plot_cell_trajectory(cds, color_by = "cell_type") +
  facet_wrap("~cell_type", nrow = 3) + scale_color_manual(values = npg_extended)

ggsave("t_cells_traj_cell_type.pdf", plot = p,width=8,height=10)

p <- plot_cell_trajectory(cds, color_by = "Pseudotime") 

# 保存图片
ggsave("trajectory_plot_t_pseudo.png", plot = p, width = 16, height = 6)
