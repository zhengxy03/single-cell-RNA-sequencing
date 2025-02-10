tumor <- subset(merged_seurat_obj, subset = sample_type == "Tumor")
identity_mapping <- c(
  "Sample1" = "ⅢB",
  "Sample3" = "ⅢA",
  "Sample5" = "ⅡA",
  "Sample7" = "ⅡA",
  "Sample9" = "ⅢB",
  "Sample11" = "ⅡB",
  "Sample13" = "ⅢB"
)
period <- identity_mapping[tumor@meta.data$orig.ident]
tumor@meta.data$period <- period

#proportion
cell_counts2 <- tumor@meta.data %>%
  group_by(sample, ident) %>%
  summarise(count = n()) %>%
  ungroup()

sample_totals <- cell_counts2 %>%
  group_by(sample) %>%
  summarise(total = sum(count)) %>%
  ungroup()

cell_proportion2 <- cell_counts2 %>%
  left_join(sample_totals, by = "sample") %>%
  mutate(proportion = count / total)


ggplot(cell_proportion2, aes(x = ident, y = proportion, fill = sample)) +
  geom_bar(stat = "identity", position = "dodge") +  # 绘制条形图
  labs(x = "Cell Type", y = "Proportion", fill = "Sample Type") +  # 设置坐标轴和图例标题
  theme_classic() +  # 使用经典主题
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # 调整横轴标签角度
    axis.line = element_line(color = "black"),  # 设置坐标轴颜色
    panel.grid.major = element_blank(),  # 移除主要网格线
    panel.grid.minor = element_blank(),  # 移除次要网格线
    legend.position = "right"  # 设置图例位置
  )


#t cells
t_cells <- subset(tumor, subset = ident == "T cell")
t_cells <- NormalizeData(t_cells)
t_cells <- FindVariableFeatures(t_cells, selection.method = "vst", nfeatures = 2000)
t_cells <- ScaleData(t_cells, features = rownames(t_cells))

t_cells <- RunPCA(t_cells, features = VariableFeatures(object = t_cells))
ElbowPlot(t_cells)
t_cells <- RunUMAP(t_cells, dims = 1:15)

t_cells <- FindNeighbors(t_cells, dims = 1:15)
seq <- seq(0.3, 1.5, by = 0.1)
for (res in seq){
    merged_seurat_obj <- FindClusters(merged_seurat_obj, resolution = res)
}
clustree(t_cells, prefix = 'RNA_snn_res.') + coord_flip()
DimPlot(t_cells, group.by = 'RNA_snn_res.0.5', label = T) +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  ggtitle(NULL) +
  guides(color = guide_legend(title = NULL))
#t_cells <- FindClusters(t_cells, resolution = 0.5)
DimPlot(t_cells, reduction = "umap", group.by = "seurat_clusters", label = TRUE, pt.size = 1.5) +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  ggtitle(NULL) +
  guides(color = guide_legend(title = NULL))

t_cell_markers <- FindAllMarkers(t_cells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "wilcox")
t_significant_markers <- subset(t_cell_markers, p_val_adj < 0.05)
write.csv(t_significant_markers, "t_all_marker.csv")

#t cells annotation
identity_mapping <- c("0" = "CD8+Tex", "1" = "CD8+Tem", "2" = "CD8+Tex", "3" = "CD8+Tcm", "4" = "CD4+Treg", "5" = "CD4+Treg", "6" = "CD8+Tex", "7" = "CD8+Trm", "8" = "CD4+Treg", "9" = "CD8+Teff", "10" = "TNK", "11" = "CD8+Teff", "12" = "CD4+Treg", "13" = "CD8+Tcm")

types <- identity_mapping[t_cells@meta.data$seurat_clusters]
t_cells@meta.data$types <- types
DimPlot(t_cells, reduction = "umap", group.by = "types", label = TRUE, pt.size = 1) +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  ggtitle(NULL) +
  guides(color = guide_legend(title = NULL))


#trajectory

library(slingshot)
library(ggsci)
#convert Seurat_obj to SingleCellExperiment_obj
tumor_sce <- as.SingleCellExperiment(tumor)

#Run slightshot
tumor_sce <- slingshot(tumor_sce, clusterLabels = "seurat_clusters", reducedDim = "UMAP")
trajectories <- slingCurves(tumor_sce)

#plot
umap_data <- as.data.frame(Embeddings(tumor, reduction = "umap"))
umap_data$period <- tumor$period
curve_data <- as.data.frame(trajectories[[1]]$s[trajectories[[1]]$ord, ])

ggplot(umap_data, aes(x = umap_1, y = umap_2, color = period)) +
  geom_point(size = 1, alpha = 0.8, shape = 16) +
  geom_path(data = curve_data, aes(x = umap_1, y = umap_2), 
            color = "black", size = 1.5, linetype = "solid") +
  scale_color_npg() +  # 使用 Nature Publishing Group 配色
  theme_classic() +
  ggtitle("Trajectory Analysis of Tumor Samples") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), 
    legend.position = "right",
    legend.title = element_blank()  # 隐藏图例标题
  )
