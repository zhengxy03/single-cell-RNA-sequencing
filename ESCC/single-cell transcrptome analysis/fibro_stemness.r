library(Seurat)
library(GSVA)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(ggsci)
fibroblasts <- readRDS("fibro_modify.rds")


stem_pathways <- read.gmt("stem_genesets.v2024.1.Hs.gmt")
stemness_genes <- unique(stem_pathways$gene)

fibroblasts_genes <- rownames(fibroblasts)
matched_genes <- stemness_genes[stemness_genes %in% fibroblasts_genes]
length(matched_genes)
stemness_genes <- matched_genes

fibroblasts <- AddModuleScore(fibroblasts, features = list(stemness_genes), name = "Stemness_Score")
head(fibroblasts@meta.data)
summary(fibroblasts@meta.data$Stemness_Score1)
#fibroblasts@meta.data$Stemness_Score_Z <- scale(fibroblasts@meta.data$Stemness_Score1)
#summary(fibroblasts@meta.data$Stemness_Score_Z)

npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(4)



p <- FeaturePlot(
  object = fibroblasts,
  features = "Stemness_Score1",
  label = TRUE,
  repel = TRUE,
  order = TRUE,
  pt.size = 1
) +
  scale_color_gradient2(
    low = "#4DBBD5FF",  # 低值颜色
    mid = "white",      # 中值颜色
    high = "#E64B35FF", # 高值颜色
    midpoint = 0,       # 中值点（对应Stemness_Score1的中值）
    limits = c(-0.6, 0.8)  # 固定图例范围
  )
ggsave("fibro_stem_umap.png", plot = p)

npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(15)
p <- VlnPlot(
  object = fibroblasts,
  features = "Stemness_Score_Z",
  group.by = "cell_type",
  pt.size = 0,
  cols = npg_extended
) +
  theme(legend.position = "none",
  axis.title.x = element_blank(), 
  axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("fibro_stem_vln.png", plot = p)


#stemness_genes <- c(
    "SOX2",    # 维持干性和多能性的关键转录因子
    "NANOG",   # 多能性标志基因

    "MYC",     # 调控细胞增殖和干性
    "KLF4",    # 多能性相关基因
    "CD44",    # 肿瘤干细胞标志物
    "ALDH1A1", # 醛脱氢酶，与干性相关
    "BMI1",    # 调控干细胞自我更新
    "NOTCH1",  # Notch信号通路，与干性维持相关
    "ESR1",     # 雌激素受体，与某些癌症的干性相关
    "KDM5B",
    "EPAS1",
    "CTNNB1",
    "HIF1A",
    "EZH2"
)


