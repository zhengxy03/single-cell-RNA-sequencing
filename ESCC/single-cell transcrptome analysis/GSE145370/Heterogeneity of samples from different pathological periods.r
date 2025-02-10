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






t



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
