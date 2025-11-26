# difference between different response status
use datasets GSE241934, GSE207422
## merge
```
library(Seurat)
seurat_obj1 <- readRDS("GSE241934.rds")
seurat_obj2 <- readRDS("GSE207422.rds")
seurat_obj3 <- readRDS("GSE146100.rds")
seurat_obj1$orig.ident <- "GSE241934"
seurat_obj2$orig.ident <- "GSE207422"
seurat_obj3$orig.ident <- "GSE146100"

merged_seurat_obj <- merge(seurat_obj1, y = c(seurat_obj2,seurat_obj3), add.cell.ids = c("GSE241934", "GSE207422","GSE146100"))
merged_seurat_obj <- JoinLayers(merged_seurat_obj)
#33460 features across 394115 samples within 1 assay

# 基因过滤
counts_matrix <- GetAssayData(merged_seurat_obj, slot = "counts", assay = "RNA")
expressed_cells_per_gene <- Matrix::rowSums(counts_matrix > 0)

min_cell_percentage <- 0.001
min_cell_count <- 10
keep_genes <- expressed_cells_per_gene >= min_cell_percentage * ncol(merged_seurat_obj) & 
              expressed_cells_per_gene >= min_cell_count
merged_seurat_obj <- merged_seurat_obj[keep_genes, ]
#20645 features across 394115 samples within 1 assay

#average genes per cell
genes_per_cell <- merged_seurat_obj$nFeature_RNA
average_genes <- mean(genes_per_cell)
#1822.567

columns_to_remove <- c(
  "RNA_snn_res.1", "RNA_snn_res.0.6", "RNA_snn_res.0.8",
  "Patient", "Pathology", "PD1.Antibody", "Chemotherapy", 
  "Pathologic.Response", "Residual.Tumor", "RECIST",
  "Resource","PD.L1.TPS","Pathological.Response.Rate",
  "major_cell_type","major.cell.type","seurat_clusters",
  "Cycles","cell.type","PD1","EGFR","Smoking_History"
)
merged_seurat_obj@meta.data <- merged_seurat_obj@meta.data[, 
  !colnames(merged_seurat_obj@meta.data) %in% columns_to_remove]
unique(colnames(merged_seurat_obj@meta.data))
print(table(merged_seurat_obj@meta.data$Response, useNA = "always"))

library(dplyr)
merged_seurat_obj@meta.data <- merged_seurat_obj@meta.data %>%
  mutate(Response = case_when(
    is.na(Response) & pathological_response == "pCR" ~ "YES",
    is.na(Response) & pathological_response %in% c("MPR", "non-MPR") ~ "NO",
    TRUE ~ Response  # 保持原有值不变
  ))
print(table(merged_seurat_obj@meta.data$Response, useNA = "always"))


library(harmony)
merged_seurat_obj <- NormalizeData(merged_seurat_obj)
merged_seurat_obj <- FindVariableFeatures(merged_seurat_obj, nfeatures = 2000)

hvgs <- VariableFeatures(merged_seurat_obj)
merged_seurat_obj <- ScaleData(merged_seurat_obj, features = hvgs)
merged_seurat_obj <- RunPCA(merged_seurat_obj, features = hvgs, npcs = 20)
merged_seurat_obj <- RunHarmony(merged_seurat_obj, "orig.ident")
saveRDS(merged_seurat_obj,file="response_pre.rds")
```
## clustering and annotation
```
png("elbowplot.png", width = 800, height = 600)
elbowplot <- ElbowPlot(merged_seurat_obj)
print(elbowplot)
dev.off()

library(ggplot2)
library(ggsci)
merged_seurat_obj <- RunUMAP(merged_seurat_obj, reduction = "harmony", dims = 1:20)
merged_seurat_obj <- FindNeighbors(merged_seurat_obj, reduction = "harmony", dims = 1:20)
merged_seurat_obj <- FindClusters(merged_seurat_obj, resolution = 0.2)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(19)
seurat_clusters <- as.character(unique(merged_seurat_obj@meta.data$seurat_clusters))  # 转换为字符向量
num_legend_items <- length(seurat_clusters)  # 图例的个数
max_label_length <- max(nchar(seurat_clusters))  # 图例名称的最大长度

base_width <- 3000  # 基础宽度
base_height <- 3000  # 基础高度
legend_width_factor <- 100  # 每个图例项增加的宽度
label_length_factor <- 10  # 每个字符增加的宽度
dynamic_width <- base_width + (num_legend_items * legend_width_factor) + (max_label_length * label_length_factor)
pdf("clusters.pdf", width = dynamic_width/300, height = base_height/300)  # 转换为英寸

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
        axis.title.y = element_text(size = 36, face = "bold", color = "black", margin = margin(r = 20)),
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

library(Seurat)
library(ggplot2)
library(ggsci)
library(dplyr)
merged_seurat_obj <- readRDS("response_umap.rds")
identity_mapping <- c(
    "0" = "T",
    "1" = "B",
    "2" = "Macrophage",
    "3" = "Fibroblast",
    "4" = "Alveolar_type_II_Epi",
    "5" = "Dedifferentiate_Epi",
    "6" = "Plasma",
    "7" = "Pericyte",
    "8" = "Neutrophil", 
    "9" = "Endothelial",
    "10" = "Monocyte/Macrophage",
    "11" = "Mast",
    "12" = "Proliferating",
    "13" = "Ciliated_1",
    "14" = "Basal_Epithelial",
    "15" = "Ciliated_2",
    "16" = "pDC",
    "17" = "Aerocyte",
    "18" = "Alveolar_type_II_Epi"
)

cell_type <- identity_mapping[merged_seurat_obj@meta.data$seurat_clusters]
merged_seurat_obj@meta.data$cell_type <- cell_type

npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(18)

cell_types <- unique(merged_seurat_obj@meta.data$cell_type)
num_legend_items <- length(cell_types)  # 图例的个数
max_label_length <- max(nchar(cell_types))  # 图例名称的最大长度

base_width <- 3000  # 基础宽度
base_height <- 3000  # 基础高度
legend_width_factor <- 100  # 每个图例项增加的宽度
label_length_factor <- 10  # 每个字符增加的宽度
dynamic_width <- base_width + (num_legend_items * legend_width_factor) + (max_label_length * label_length_factor)
pdf("annotation.pdf", width = dynamic_width/300, height = base_height/300)
DimPlot(merged_seurat_obj, reduction = "umap", label = TRUE, pt.size = 1, group.by = "cell_type", label.size = 4) +
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    ggtitle(NULL) +
    scale_color_manual(values = npg_extended) +
    coord_fixed(ratio = 1) +
    guides(color = guide_legend(title = NULL, override.aes = list(size = 5))) +
    theme(
        text = element_text(size = 8, face = "bold"),
        axis.text.x = element_text(size = 28, color = "black"),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title.x = element_text(size = 36, face = "bold", color = "black"),
        axis.title.y = element_text(size = 36, face = "bold", color = "black", margin = margin(r = 20)),  # 增加右侧间距
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 16, face = "bold", color = "black"),
        legend.title = element_text(size = 16, face = "bold", color = "black"),
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
saveRDS(merged_seurat_obj,file="response_anno.rds")
```
## celltype differnence between response and non-response
```
plot_data <- merged_seurat_obj
plot_data$Response_plot <- ifelse(
  plot_data$Response == "YES", "Good", "Poor"
)

npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(18)
pdf("annotation_by_response.pdf", width = 8000/300, height = 3000/300)

p <- DimPlot(plot_data, 
             reduction = "umap", 
             label = TRUE, 
             pt.size = 1, 
             group.by = "cell_type", 
             label.size = 4,
             split.by = "Response_plot",  # 使用临时列
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


#proportion
proportion_data <- plot_data@meta.data %>%
    group_by(Response_plot, cell_type) %>% summarise(count = n()) %>% mutate(proportion = count / sum(count))

npg_extended <- colorRampPalette(npg_pal)(18)
pdf("response_prop.pdf", width = 6000/300, height = 3000/300)  # 设置高分辨率和尺寸
ggplot(proportion_data, aes(x = "", y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", width = 1) +         # 堆叠柱状图
  coord_polar(theta = "y") +                       # 转换为饼图
  scale_fill_manual(values = npg_extended) +       # 使用自定义颜色
  theme_void() +                                   # 空白背景
  labs(fill = "Cell Type") +
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
    ~ Response_plot,
    ncol = 2,
    strip.position = "bottom"                      # 标签放在下方
  )
dev.off()

#boxplot
library(ggpubr)
library(dplyr)
library(ggsci)

cell_counts <- plot_data@meta.data %>%
  group_by(patients, cell_type, Response_plot) %>%
  summarise(count = n()) %>%
  ungroup()

# 计算每个样本的总细胞数
sample_totals <- cell_counts %>%
  group_by(patients) %>%
  summarise(total = sum(count)) %>%
  ungroup()

# 将样本总细胞数合并到 cell_counts 中
cell_counts <- cell_counts %>%
  left_join(sample_totals, by = "patients")

# 计算细胞类型比例
cell_proportions <- plot_data@meta.data %>%
  group_by(patients, cell_type, Response_plot) %>%
  summarise(count = n()) %>%  # 计算每种细胞类型的数量
  ungroup() %>%
  left_join(sample_totals, by = "patients") %>%  # 合并样本总细胞数
  mutate(proportion = count / total)  # 计算细胞类型比例

# 对每个细胞类型构建列联表并进行卡方检验
chisq_results <- cell_proportions %>%
  group_by(cell_type) %>%
  summarise(
    p_value = {
      # 提取正常样本和肿瘤样本的细胞数量和总细胞数
      normal_count = sum(count[Response_plot == "Good"])
      tumor_count = sum(count[Response_plot == "Poor"])
      normal_total = sum(total[Response_plot == "Good"])
      tumor_total = sum(total[Response_plot == "Poor"])
      
      # 构建 2x2 列联表
      cont_table <- matrix(
        c(normal_count, tumor_count, normal_total - normal_count, tumor_total - tumor_count),
        nrow = 2
      )
      
      # 进行卡方检验
      chisq.test(cont_table)$p.value
    }
  ) %>%
  ungroup()

# 添加显著性标记
chisq_results <- chisq_results %>%
  mutate(
    significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

# 绘制基于比例的箱线图并添加显著性标记
pdf("boxplot_proportion.pdf", width = 6000/300, height = 3000/300)
ggplot(cell_proportions, aes(x = cell_type, y = proportion, fill = Response_plot)) +
  geom_boxplot() +  # 绘制箱线图
  geom_text(
    data = chisq_results, 
    aes(x = cell_type, y = max(cell_proportions$proportion) * 1.05, label = significance),  # 调整 y 值
    size = 12, 
    vjust = 0.5,  # 调整 vjust 参数
    inherit.aes = FALSE  # 忽略父图层的 aes 映射
  ) +  # 添加显著性标记
  labs(x = "", y = "Cell Proportion", fill = "Response") +  # 设置坐标轴和图例标题
  theme_classic() +  # 使用经典主题
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 28),  # 调整横轴标签角度
    axis.text.y = element_text(size = 28),
    axis.title.y = element_text(size = 40),
    axis.line = element_line(color = "black"),  # 设置坐标轴颜色
    panel.grid.major = element_blank(),  # 移除主要网格线
    panel.grid.minor = element_blank(),  # 移除次要网格线
    legend.text = element_text(size = 36),         # 图例文本字体大小
    legend.title = element_text(size = 40),
    legend.position = "right"  # 设置图例位置在右侧
  ) +
  scale_fill_npg() +  # 使用 npg 配色
  scale_y_continuous(limits = c(0, 1))  # 设置 y 轴范围为 0 到 1
dev.off()
```
## LUAD and LUSC response difference
```
#LUAD
library(Seurat)
library(ggplot2)
library(ggsci)
library(dplyr)
# 创建整合的plot_data
plot_data <- merged_seurat_obj
plot_data$Response_plot <- ifelse(
  plot_data$Response == "YES", "Good", "Poor"
)
plot_data$cancer_response <- paste(plot_data$cancer_type, plot_data$Response_plot, sep = "_")

npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(18)
pdf("combined_annotation_by_cancer_response.pdf", width = 12000/300, height = 3000/300)

p <- DimPlot(plot_data, 
             reduction = "umap", 
             label = TRUE, 
             pt.size = 1, 
             group.by = "cell_type", 
             label.size = 4,
             split.by = "cancer_response",  # 同时按癌症类型和应答分组
             ncol = 4) +  # 改为4列：LUAD_Good, LUAD_Poor, LUSC_Good, LUSC_Poor
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

#stickplot
plot_data <- merged_seurat_obj
plot_data$Response_plot <- ifelse(
  plot_data$Response == "YES", "Good", "Poor"
)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(18)
proportion_data <- plot_data@meta.data %>%
    group_by(cancer_type, Response_plot, cell_type) %>%
    summarise(count = n()) %>%
    mutate(proportion = count / sum(count)) %>%
    # 创建组合变量
    mutate(cancer_sample = paste0(cancer_type, "-", Response_plot)) %>%
    # 设置组合变量的顺序
    mutate(cancer_sample = factor(cancer_sample,
                                 levels = c("LUAD-Good", "LUAD-Poor",
                                           "LUSC-Good", "LUSC-Poor")))

pdf("Response_cancer_type.pdf", width = 8000/300, height = 3000/300) 

ggplot(proportion_data, aes(x = cancer_sample, y = proportion, fill = cell_type)) +
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
```
## immune
```
immune <- subset(merged_seurat_obj,subset = cell_type %in% c("B","Macrophage","Mast","Monocyte/Macrophage","Neutrophil","pDC","Plasma","T"))
#20645 features across 318281 samples within 1 assay
saveRDS(immune,file="allcellsets_immune.rds")