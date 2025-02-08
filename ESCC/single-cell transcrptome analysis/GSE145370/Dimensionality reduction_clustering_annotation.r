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
sample <- c("S133T", "S133N", "S134T", "S134N", "S135T", "S135N", "S149T", "S149N", "S150T", "S150N", "S158T", "S158N", "S159T", "S159N")
merged_seurat_obj@meta.data$sample <- sample[as.numeric(factor(merged_seurat_obj@meta.data$orig.ident))]

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

ggplot(proportion_data2, aes(x = "", y = count, fill = ident)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  facet_wrap(~ sample_type, nrow = 2) +  # 根据样本进行分面，可调整 nrow 参数
  theme_void() +
  labs(title = "不同样本的细胞组成饼图")
