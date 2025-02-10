#plot according to sample types
p1 <- DimPlot(merged_seurat_obj, reduction = "umap", group.by = "sample_type", pt.size = 0.5) +
  # 设置横纵坐标标签
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  # 移除标题
  ggtitle(NULL) +
  # 移除图例标题（如果需要）
  guides(color = guide_legend(title = NULL))

print(p1)

#plot according to sample sources
p2 <- DimPlot(merged_seurat_obj, reduction = "umap", group.by = "sample_sources", pt.size = 0.5) +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  ggtitle(NULL) +
  guides(color = guide_legend(title = NULL))

print(p2)

#cluster
ElbowPlot(merged_seurat_obj)
merged_seurat_obj <- FindNeighbors(merged_seurat_obj, reduction = "harmony", dims = 1:20)
seq <- seq(0.3, 1.5, by = 0.1)
for (res in seq){
    merged_seurat_obj <- FindClusters(merged_seurat_obj, resolution = res)
}

library(clustree)
library(patchwork)
p3 <- clustree(merged_seurat_obj, prefix = 'RNA_snn_res.') + coord_flip()
p4 <- DimPlot(merged_seurat_obj, group.by = 'RNA_snn_res.0.5', label = T) +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  ggtitle(NULL) +
  guides(color = guide_legend(title = NULL))
#merged_seurat_obj <- FindClusters(merged_seurat_obj, resolution = 0.5)
#merged_seurat_obj <- RunUMAP(merged_seurat_obj, reduction = "harmony", dims = 1:20)
#DimPlot(merged_seurat_obj, reduction = "umap", label = TRUE)

#annotation
markers <- FindAllMarkers(object = merged_seurat_obj, 
                                  test.use = "roc", 
                                  only.pos = TRUE, 
                                  min.pct = 0.25, 
                                  thresh.use = 0.25)

significant_markers <- subset(markers, myAUC > 0.7)
write.csv(significant_markers,"marker.csv")

new.cluster.ids <- c("T cell", "T cell", "Neutrophils", "B cell", "T cell","NK", "Mono/Macro", "T cell", "mDC", "Plasma", "Neutrophils", "Mast cell", "Epithelial cell", "T cell", "B cell", "Mono/Macro", "pDC", "Fibroblast", "Proliferating cell", "NK")
names(new.cluster.ids) <- levels(merged_seurat_obj)
merged_seurat_obj <- RenameIdents(merged_seurat_obj, new.cluster.ids)
p5 <- DimPlot(merged_seurat_obj, reduction = "umap", label = TRUE) +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  ggtitle(NULL) +
  guides(color = guide_legend(title = NULL)) +
  theme(text = element_text(size = 8, face = "bold"))
merged_seurat_obj@meta.data$ident <- Idents(merged_seurat_obj)

#DEGs
genes_to_plot <- c(
    "CD3D", "CD3E", "CD8A", # T细胞（T）
    "CSF3R", "S100A8", "S100A9", # 中性粒细胞（neu）
    "CD19", "CD79A", "MS4A1", # B细胞（B）
    "GNLY", "KLRD1", "FGFBP2", # NK细胞（NK）
    "CD68", "C1QA", "C1QB", # 单核细胞（mono）
    "CD1E", "CD1C", "S100B", # 髓系树突细胞（mDC）
    "IGHG1", "IGHG3", "JCHAIN", # 浆细胞（plasma）
    "CPA3", "TPSAB1", "TPSB2", # 肥大细胞（mast cell）
    "KRT15", "KRT19", "KRT17", # 上皮细胞（Epi）
    "LILRA4", "CLEC4C", "CXCR3", # 浆细胞样树突细胞（pDC）
    "SERP2", "COL1A2", "COL3A1", # 成纤维细胞（fibro）
    "TOP2A", "TPX2", "MKI67" # 增殖细胞（proliferating）
)
DotPlot(merged_seurat_obj, features = genes_to_plot, group.by = "ident") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

FeaturePlot(merged_seurat_obj, features = genes_to_plot)

DoHeatmap(merged_seurat_obj, features = genes_to_plot, 
          group.by = "ident", group.bar.height = 0.01, size = 3, angle = 90, hjust = 0, label = TRUE) + scale_fill_gradient2(low = "blue", mid = "white", high = "red")

#cell proportion
identity_mapping <- c(
  "Sample1" = "S133T",
  "Sample2" = "S133N",
  "Sample3" = "S134T",
  "Sample4" = "S134N",
  "Sample5" = "S135T",
  "Sample6" = "S135N",
  "Sample7" = "S149T",
  "Sample8" = "S149N",
  "Sample9" = "S150T",
  "Sample10" = "S150N",
  "Sample11" = "S158T",
  "Sample12" = "S158N",
  "Sample13" = "S159T",
  "Sample14" = "S159N"
)
sample <- identity_mapping[merged_seurat_obj@meta.data$orig.ident]
merged_seurat_obj@meta.data$sample <- sample

proportion_data <- merged_seurat_obj@meta.data %>%
  group_by(sample, ident) %>% summarise(count = n()) %>% mutate(proportion = count / sum(count))

ggplot(proportion_data, aes(x = sample, y = proportion, fill = ident)) +
  geom_bar(stat = "identity", position = "stack") +  # 堆叠柱状图
  labs(x = "Sample Type", y = "Proportion", fill = "Cell Type") +  # 设置坐标轴和图例标题
  theme_classic() +  # 使用经典主题（独立坐标轴）
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # 调整横轴标签角度
    axis.line = element_line(color = "black"),  # 设置坐标轴颜色
    panel.grid.major = element_blank(),  # 移除主要网格线
    panel.grid.minor = element_blank()   # 移除次要网格线
  )

proportion_data2 <- merged_seurat_obj@meta.data %>%
  group_by(sample_type, ident) %>% summarise(count = n()) %>% mutate(proportion = count / sum(count))

ggplot(proportion_data2, aes(x = "", y = proportion, fill = ident)) +
  geom_bar(stat = "identity", width = 1) +  # 堆叠柱状图
  coord_polar(theta = "y") +  # 转换为极坐标，生成饼图
  labs(title = "Proportion of Cell Types", fill = "Cell Type") +  # 设置标题和图例标题
  theme_void() +  # 使用空白主题
  theme(
    legend.position = "right",  # 设置图例位置
    plot.title = element_text(hjust = 0.5)  # 设置标题居中
  ) +
  facet_wrap(~ sample_type, ncol = 2)

cell_counts <- merged_seurat_obj@meta.data %>%
  group_by(sample, ident) %>%
  summarise(count = n()) %>%
  ungroup()

cell_counts <- cell_counts %>%
  mutate(sample_type = ifelse(grepl("T$", sample), "Tumor", "Normal"))

sample_totals <- cell_counts %>%
  group_by(sample) %>%
  summarise(total = sum(count)) %>%
  ungroup()

cell_proportion <- cell_counts %>%
  left_join(sample_totals, by = "sample") %>%
  mutate(proportion = count / total)

ggplot(cell_proportion, aes(x = ident, y = proportion, fill = sample_type)) +
  geom_boxplot() +  # 绘制箱线图
  labs(x = "Cell Type", y = "Proportion", fill = "Sample Type") +  # 设置坐标轴和图例标题
  theme_classic() +  # 使用经典主题
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # 调整横轴标签角度
    axis.line = element_line(color = "black"),  # 设置坐标轴颜色
    panel.grid.major = element_blank(),  # 移除主要网格线
    panel.grid.minor = element_blank(),  # 移除次要网格线
    legend.position = "right"  # 设置图例位置在底部
  ) +
  scale_fill_manual(values = c("Tumor" = "red", "Normal" = "blue"))
