emt_genes <- c("VIM", "CDH1", "ZEB1", "SNAI1", "TWIST1", "FN1", "CDH2",
                "SNAI2", "ZEB2", "MMP2", "MMP9", "BCL9", "CTNNB1")

emt_genes <- emt_genes[emt_genes %in% rownames(epi)]

# 计算 EMT 评分
epi <- AddModuleScore(
  object = epi,
  features = list(emt_genes),
  name = "EMT_Score"
)
head(epi@meta.data)

npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(4)

epi@meta.data$EMT_Score1_scaled <- scale(epi@meta.data$EMT_Score1)
summary(epi@meta.data$EMT_Score1_scaled)

# 绘制 UMAP 分布图
p <- FeaturePlot(
  object = epi,          # Seurat 对象
  features = "EMT_Score1_scaled",            # EMT 评分列名
  label = TRUE,                       # 是否显示聚类标签
  repel = TRUE,                       # 避免标签重叠
  order = TRUE,                       # 按评分高低排序
  cols = c("#4DBBD5FF", "#E64B35FF"),
  pt.size = 0.5                       # 点的大小
) + 
  ggtitle("EMT Score UMAP Distribution")  # 添加标题

ggsave("emt-plot.png", plot = p)

# 绘制小提琴图
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(17)
p <- VlnPlot(
  object = epi,
  features = "EMT_Score1",
  group.by = "Epi_cluster",
  pt.size = 0,
  cols = npg_extended
)
ggsave("epi_emt-plot2.png", plot = p)

png("epi_emt_dotplot.png", width = 8000, height = 3000, res = 300)
DotPlot(epi, 
        features = emt_genes, 
        group.by = "Epi_cluster",
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