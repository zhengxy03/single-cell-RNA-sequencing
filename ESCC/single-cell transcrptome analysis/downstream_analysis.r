#GSE203115
merged_seurat_obj <- NormalizeData(merged_seurat_obj)
merged_seurat_obj <- FindVariableFeatures(merged_seurat_obj, selection.method = "vst", nfeatures = 2000)
merged_seurat_obj <- ScaleData(merged_seurat_obj, features = rownames(merged_seurat_obj))

#PCA
merged_seurat_obj <- RunPCA(merged_seurat_obj, features = VariableFeatures(object = merged_seurat_obj))

#tSNE
merged_seurat_obj <- RunTSNE(merged_seurat_obj, dims = 1:10)

#UMAP
merged_seurat_obj <- RunUMAP(merged_seurat_obj, dims = 1:10)

#unsupervised clusting
merged_seurat_obj <- FindNeighbors(merged_seurat_obj, dims = 1:10)
merged_seurat_obj <- FindClusters(merged_seurat_obj, resolution = 0.2)

DimPlot(merged_seurat_obj, reduction = "umap", label = TRUE)

#differentiate non-responder and responder
DimPlot(merged_seurat_obj, 
        reduction = "umap", 
        group.by = "response_status",  
        cols = c("Non-Responder" = "red", "Responder" = "blue"),
        pt.size = 1)

#subsets of non-responder and responder
non_responder_obj <- subset(merged_seurat_obj, subset = response_status == "Non-Responder")
responder_obj <- subset(merged_seurat_obj, subset = response_status == "Responder")

num_clusters <- length(unique(non_responder_obj$seurat_clusters))
grey_colors <- rep("grey", num_clusters)
DimPlot(non_responder_obj, 
        reduction = "umap", 
        group.by = "seurat_clusters",  
        cols = grey_colors,  
        pt.size = 1.5,  
        label = TRUE) +  
  ggtitle("Non-Responder") +
  theme(legend.position = "none")

#DEGs
markers <- FindAllMarkers(merged_seurat_obj, 
                          only.pos = TRUE, 
                          min.pct = 0.1, 
                          logfc.threshold = 0.25, 
                          test.use = "wilcox")
#significant markers
significant_markers <- subset(markers, p_val_adj < 0.05)
write.csv(significant_markers,"marker.csv")

#cell annotation
new.cluster.ids <- c("T cell", "Fibroblast", "Myeloid", "Endothelia cell", "Epithelia cell","Epithelia cell", "B cell & mast cell", "SMC & pericyte cell", "plasma", "Double cell")
names(new.cluster.ids) <- levels(merged_seurat_obj)
merged_seurat_obj <- RenameIdents(merged_seurat_obj, new.cluster.ids)
DimPlot(merged_seurat_obj, reduction = "umap", label = TRUE, pt.size = 0.5)

#get T cells subset
t_cells <- subset(merged_seurat_obj, subset = seurat_clusters == "0")
t_cells <- NormalizeData(t_cells)
t_cells <- FindVariableFeatures(t_cells, selection.method = "vst", nfeatures = 2000)
t_cells <- ScaleData(t_cells, features = rownames(t_cells))

#PCA
t_cells <- RunPCA(t_cells, features = VariableFeatures(object = t_cells))

#UMAP
t_cells <- RunUMAP(t_cells, dims = 1:10)


#find T clusters
t_cells <- FindNeighbors(t_cells, dims = 1:10)
t_cells <- FindClusters(t_cells, resolution = 0.5)

DimPlot(t_cells, reduction = "umap", label = TRUE)
#annotation
t_cell_markers <- FindAllMarkers(t_cells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "wilcox")
t_significant_markers <- subset(t_cell_markers, p_val_adj < 0.05)
write.csv(t_significant_markers, "t_all_marker.csv")

#t cells annotation
new.cluster.ids <- c("CD8+Tem", "CD8+Tex", "CD4+Tcm", "CD4+Treg", "CD8+Tex", "CD8+Tex", "Cycling T", "Cycling T", "Cycling T", "NK")
names(new.cluster.ids) <- levels(t_cells)
t_cells <- RenameIdents(t_cells, new.cluster.ids)
DimPlot(t_cells, reduction = "umap", label = TRUE, pt.size = 0.5)

#VlnPlot
genes_to_plot <- c("GZMK", "CXCL13", "CTLA4", "FOXP3", "ICOS", "SELL", "HAVCR2", "PRDM1")
VlnPlot(t_cells, features = genes_to_plot, group.by = "seurat_clusters", pt.size = 0)

#bubble heat map
DotPlot(t_cells, features = genes_to_plot, group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Trajectory analysis
#Monocle2
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("monocle")
library(monocle)  #版本不兼容，看看以后会不会更新🌚

#Slightshot
BiocManager::install("slingshot")
library(slightshot)

CD8 <- subset(t_cells, subset = seurat_clusters %in% c(0, 1, 4, 5))
CD4 <- subset(t_cells, subset = seurat_clusters %in% c(2,3))
#convert Seurat_obj to SingleCellExperiment_obj
cd8_sce <- as.SingleCellExperiment(CD8)

#Run slightshot
cd8_sce <- slingshot(cd8_sce, clusterLabels = "seurat_clusters", reducedDim = "UMAP")
trajectories <- slingCurves(cd8_sce)

#plot
umap_data <- as.data.frame(Embeddings(CD8, reduction = "umap"))
umap_data$response_status <- CD8$response_status
curve_data <- as.data.frame(trajectories[[1]]$s[trajectories[[1]]$ord, ])

ggplot(umap_data, aes(x = umap_1, y = umap_2, color = response_status)) +
  geom_point(size = 2, alpha = 0.8, shape = 16) +
  geom_path(data = curve_data, aes(x = umap_1, y = umap_2), 
            color = "black", size = 1.5, linetype = "solid") +
  scale_color_npg() +  # 使用 Nature Publishing Group 配色
  theme_classic() +
  ggtitle("Trajectory Analysis of CD8 T Cells by Response Status") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), 
    legend.position = "right",
    legend.title = element_blank()  # 隐藏图例标题
  )




#GSE196756
#quality control
#calculate the proportion of mitochondria
merged_seurat_obj[["percent.mt"]] <- PercentageFeatureSet(merged_seurat_obj, pattern = "^MT-")
#filter low qyality cells
merged_seurat_obj <- subset(merged_seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 50)
print(paste("过滤后的细胞数量:", ncol(merged_seurat_obj)))
#33817

#further analysis
merged_seurat_obj <- NormalizeData(merged_seurat_obj)
merged_seurat_obj <- FindVariableFeatures(merged_seurat_obj, nfeatures = 2000)
merged_seurat_obj <- ScaleData(merged_seurat_obj)
merged_seurat_obj <- RunPCA(merged_seurat_obj, npcs = 50)

#correct the batch effects
library(harmony)
merged_seurat_obj <- RunHarmony(merged_seurat_obj, group.by.vars = "orig.ident", dims.use = 1:50)

#find clusters based on harmony space
merged_seurat_obj <- FindNeighbors(merged_seurat_obj, reduction = "harmony", dims = 1:50)
merged_seurat_obj <- FindClusters(merged_seurat_obj, resolution = 0.5)

#UMAP
merged_seurat_obj <- RunUMAP(merged_seurat_obj, reduction = "harmony", dims = 1:50)

#visualize based on clusters
DimPlot(merged_seurat_obj, reduction = "umap", label = TRUE, pt.size = 0.5)

#DEGs
markers <- FindAllMarkers(merged_seurat_obj, min.pct = 0.25, logfc.threshold = 0.25, only.pos = TRUE)
significant_markers <- subset(markers, p_val_adj < 0.05)
write.csv(significant_markers,"marker.csv")

#cell annotation
new.cluster.ids <- c("mast cells", "NK cells", "T cells", "B cells", "Fibroblast","Neutrophils", "Epithelial cells", "Myeloid cells", "T cells", "T cells", "Myeloid cells", "Epithelial cells", "Plasma", "Myeloid cells", "Fibroblast", "Endothelial cells", "Epithelial cells", "Fibroblast", "Fibroblast", "Mast cells", "Fibroblast", "Epithelial cells", "Epithelial cells", "Myeloid cells")
names(new.cluster.ids) <- levels(merged_seurat_obj)
merged_seurat_obj <- RenameIdents(merged_seurat_obj, new.cluster.ids)
DimPlot(merged_seurat_obj, reduction = "umap", label = TRUE, pt.size = 0.5)

#visualize based on sample types
DimPlot(merged_seurat_obj, reduction = "umap", group.by = "sample_type")

#heatmap of markers
gene_list <- c("IGHG1", "IGKC", "GNLY", "KLRD1", "KRT5", "KRT13", "TP63", "MS4A1", "CD79A", 
               "CD19", "BANK1", "CD2", "CD3D", "CD3E", "CD14", "CD68", "ENG", "PECAM1", 
               "VWF", "CPA3", "TPSAB1", "TPSB2", "COL1A1", "COL1A2", "COL14A1", "DCN", 
               "PDGFRB", "THY1", "S100A8", "FFAR2", "CSF3R", "FCGR3B", "CXCR2")
#add cell Idents into meta.data
merged_seurat_obj@meta.data$ident <- Idents(merged_seurat_obj)
#DoHeatmap
DoHeatmap(merged_seurat_obj, features = gene_list, 
          group.by = "ident", group.bar.height = 0.01, size = 3, angle = 90, hjust = 0, label = TRUE) + scale_fill_gradient2(low = "blue", mid = "white", high = "red")

#cell proportion
proportion_data <- merged_seurat_obj@meta.data %>%
  group_by(sample_type, ident) %>% summarise(count = n()) %>% mutate(proportion = count / sum(count))

ggplot(proportion_data, aes(x = sample_type, y = proportion, fill = ident)) +
  geom_bar(stat = "identity", position = "stack") +  # 堆叠柱状图
  labs(x = "Sample Type", y = "Proportion", fill = "Cell Type") +  # 设置坐标轴和图例标题
  theme_classic() +  # 使用经典主题（独立坐标轴）
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # 调整横轴标签角度
    axis.line = element_line(color = "black"),  # 设置坐标轴颜色
    panel.grid.major = element_blank(),  # 移除主要网格线
    panel.grid.minor = element_blank()   # 移除次要网格线
  )

#get Epi subset
Epi_cells <- subset(merged_seurat_obj, subset = ident == "Epithelial cells")
#3318
Epi_cells <- NormalizeData(Epi_cells)
Epi_cells <- FindVariableFeatures(Epi_cells, selection.method = "vst", nfeatures = 2000)
Epi_cells <- ScaleData(Epi_cells, features = rownames(Epi_cells))

#PCA
Epi_cells <- RunPCA(Epi_cells, features = VariableFeatures(object = Epi_cells))

#UMAP
Epi_cells <- RunUMAP(Epi_cells, dims = 1:10)

#find Epi clusters
Epi_cells <- FindNeighbors(Epi_cells, dims = 1:10)
Epi_cells <- FindClusters(Epi_cells, resolution = 0.25)

DimPlot(Epi_cells, reduction = "umap", group.by = "seurat_clusters", label = TRUE)

#subset of tumor and normal tissues
sample_types2 <- c("T", "N", "T", "N", "T", "N")
Epi_cells@meta.data$sample_type2 <- sample_types2[as.numeric(factor(Epi_cells@meta.data$orig.ident))]

markers <- FindMarkers(Epi_cells, 
                       ident.1 = "T", 
                       ident.2 = "N", 
                       group.by = "sample_type2",
                       test.use = "wilcox",  # 使用 Wilcoxon rank sum test
                       logfc.threshold = 0.25)



#volano figure
library(ggplot2)
library(dplyr)
library(ggrepel)

top_genes <- markers %>% arrange(p_val) %>% head(10)
markers$gene <- rownames(markers)
markers$neg_log10_pval <- -log10(markers$p_val)
markers$diffexpressed <- "No"
markers$diffexpressed[markers$avg_log2FC > 0.25 & markers$p_val < 0.05] <- "Up"
markers$diffexpressed[markers$avg_log2FC < -0.25 & markers$p_val < 0.05] <- "Down"

ggplot(markers, aes(x = avg_log2FC, y = neg_log10_pval, color = diffexpressed)) +
  geom_point(size = 1.5) + 
  scale_color_manual(values = c("Up" = "red", "Down" = "green", "No" = "gray")) + 
  labs(x = "Log2 Fold Change", y = "-Log10(P-value)", color = "Expression") + 
  theme_minimal() +
  theme(legend.position = "top") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") + 
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "black") +
  geom_text_repel(
    data = top_genes,  # 选择要标注的基因
    aes(label = gene),  # 标注基因名称
    color = "black",  # 文本颜色
    size = 3,  # 文本大小
    box.padding = 0.5,  # 文本与点之间的间距
    segment.color = "black",  # 短线颜色
    segment.size = 0.5,  # 短线粗细
    min.segment.length = 0,  # 短线最小长度
    max.overlaps = Inf  # 允许的最大重叠次数
  )

#enrichment analysis
#GSVA
BiocManager::install("GSVA")
library(GSVA)

#normalized expression matrix
expression_matrix <- as.matrix(Epi_cells@assays$RNA$data)

#import gene set
#https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp#H
gmt_file <- "../../h.all.v2024.1.Hs.symbols.gmt"
gene_sets <- GSEABase::getGmt(gmt_file)

#Run GSVA
gsvaPar <- ssgseaParam(exprData = expression_matrix, 
                       geneSets = gene_sets,
                       normalize = TRUE)
gsva_results <- GSVA::gsva(gsvaPar, verbose = FALSE)

#calculate enrichment score
group_info <- Epi_cells@meta.data$sample_type2
tumor_samples <- colnames(gsva_results)[group_info == "T"]
tumor_means <- rowMeans(gsva_results[, tumor_samples, drop = FALSE])

normal_samples <- colnames(gsva_results)[group_info == "N"]
normal_means <- rowMeans(gsva_results[, normal_samples, drop = FALSE])

enrichment_matrix <- cbind(T = tumor_means, N = normal_means)
rownames(enrichment_matrix) <- gsub("^HALLMARK_", "", rownames(enrichment_matrix))
#z-score normalized
enrichment_matrix_zscore <- t(scale(t(enrichment_matrix))) 

#heatmap
#pheatmap
install.packages("pheatmap")
library(pheatmap)

pheatmap(enrichment_matrix_zscore, 
         scale = "none",  # 不进行额外标准化（因为已经 z-score 标准化）
         cluster_rows = TRUE,  # 启用行聚类
         cluster_cols = TRUE,  # 启用列聚类
         treeheight_row = 0,  # 隐藏行聚类树
         treeheight_col = 0,  # 隐藏列聚类树
         color = colorRampPalette(c("blue", "white", "red"))(100),
         show_colnames = TRUE,  # 显示列名
         show_rownames = TRUE,  # 显示行名
         angle_col = 0)

#KM-plot
#http://kmplot.com/analysis/
#choose pan-cancer




#analysis fibroblast
Fibroblast_cells <- subset(merged_seurat_obj, subset = ident == "Fibroblast")
#3635
Fibroblast_cells <- NormalizeData(Fibroblast_cells)
Fibroblast_cells <- FindVariableFeatures(Fibroblast_cells, selection.method = "vst", nfeatures = 2000)
Fibroblast_cells <- ScaleData(Fibroblast_cells, features = rownames(Fibroblast_cells))

#PCA
Fibroblast_cells <- RunPCA(Fibroblast_cells, features = VariableFeatures(object = Fibroblast_cells))

#UMAP
Fibroblast_cells <- RunUMAP(Fibroblast_cells, dims = 1:10)

#find Epi clusters
Fibroblast_cells <- FindNeighbors(Fibroblast_cells, dims = 1:10)
Fibroblast_cells <- FindClusters(Fibroblast_cells, resolution = 0.25)

DimPlot(Fibroblast_cells, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
#visualize based on sample types
DimPlot(Fibroblast_cells, reduction = "umap", group.by = "sample_type")

#DEGs
markers <- FindAllMarkers(Fibroblast_cells, min.pct = 0.25, logfc.threshold = 0.25, only.pos = TRUE)
significant_markers <- subset(markers, p_val_adj < 0.05)
write.csv(significant_markers,"fibro_marker.csv")

#cell annotation
new.cluster.ids <- c("iCAFs", "iCAFs", "myCAFs", "iCAFs", "myCAFs", "iCAFs", "myCAFs", "iCAFs", "myCAFs", "myCAFs", "iCAFs", "iCAFs")
names(new.cluster.ids) <- levels(Fibroblast_cells)
Fibroblast_cells <- RenameIdents(Fibroblast_cells, new.cluster.ids)
Fibroblast_cells$Idents <- Idents(Fibroblast_cells)
DimPlot(Fibroblast_cells, reduction = "umap", group.by = "Idents", label = TRUE, pt.size = 0.5)

#featureplot
genes_to_plot <- c("CFD", "CXCL12", "ACTA2", "TAGLN")
FeaturePlot(Fibroblast_cells, features = genes_to_plot)



#T cells analysis
t_cells <- subset(merged_seurat_obj, subset = ident == "T cells")
#6427
t_cells <- NormalizeData(t_cells)
t_cells <- FindVariableFeatures(t_cells, selection.method = "vst", nfeatures = 2000)
t_cells <- ScaleData(t_cells, features = rownames(t_cells))

#PCA
t_cells <- RunPCA(t_cells, features = VariableFeatures(object = t_cells))

#UMAP
t_cells <- RunUMAP(t_cells, dims = 1:10)


#find T clusters
t_cells <- FindNeighbors(t_cells, dims = 1:10)
t_cells <- FindClusters(t_cells, resolution = 0.2)

DimPlot(t_cells, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
#annotation
t_cell_markers <- FindAllMarkers(t_cells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "wilcox")
t_significant_markers <- subset(t_cell_markers, p_val_adj < 0.05)
write.csv(t_significant_markers, "t_marker.csv")

#t cells annotation
new.cluster.ids <- c("CD4 CCR7", "CD4 TRBV28", "CD4 Treg", "CD4 IL6R", "CD4 IL17F", "CD4 HAVCR2", "CD8-1", "CD8-2", "CD8-2", "CD8-1")
names(new.cluster.ids) <- levels(t_cells)
t_cells <- RenameIdents(t_cells, new.cluster.ids)
t_cells$Idents <- Idents(t_cells)
DimPlot(t_cells, reduction = "umap", group.by = "Idents", label = TRUE, pt.size = 0.5)

#Pseudotime analysis
library(slingshot)

#Get UMAP coordinates
umap_coords <- Embeddings(t_cells, reduction = "umap")
cell_types <- Idents(t_cells)

#Trajectory analysis
slingshot_obj <- slingshot(umap_coords, clusterLabels = cell_types)

#Trajectory curves
curves <- slingCurves(slingshot_obj)
length(curves)
#2

#Pseudotime info
pseudotime <- slingPseudotime(slingshot_obj)
head(pseudotime)

#Visualize pseudotime
library(ggplot2)
library(viridis)
library(dplyr)
library(ggsci)
library(tidyr)

#Create a data frame for plotting
plot_data <- data.frame(
  UMAP1 = umap_coords[, 1], 
  UMAP2 = umap_coords[, 2], 
  Pseudotime = pseudotime[, 1], # 使用第一个轨迹的伪时间(改成2就是第二个轨迹)
  Cluster = cell_types
)

#calculate point boundary
x_min <- min(plot_data$UMAP1)
x_max <- max(plot_data$UMAP1)
y_min <- min(plot_data$UMAP2)
y_max <- max(plot_data$UMAP2)

#Plot with modified aesthetics
lineage_colors <- c("#CCCCCC", "#666666")  #浅灰和深灰

p1 <- ggplot(plot_data, aes(x = UMAP1, y = UMAP2, color = Pseudotime)) +
  geom_point(size = 1.5, alpha = 0.7, shape = 16) +  # 调整点的大小和透明度
  scale_color_gradientn(colors = c("#F7FBFF", "#6BAED6", "#08306B")) +  # 蓝色渐变
  theme_classic(base_size = 12) +  # 使用经典主题，设置基础字号
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # 添加实线边框
    axis.line = element_line(size = 0.8),  # 加粗坐标轴线
    axis.title = element_text(size = 14, face = "bold"),  # 坐标轴标题字体
    axis.text = element_text(size = 12, color = "black"),  # 坐标轴刻度字体
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  # 图标题字体
    legend.title = element_text(size = 12, face = "bold"),  # 图例标题字体
    legend.text = element_text(size = 10),  # 图例文字字体
    legend.position = "right"  # 图例位置
  ) +
  labs(title = "Pseudotime Trajectory of T Cells", x = "UMAP1", y = "UMAP2")

#Add trajectory curves with smaller arrows
for (i in seq_along(curves)) {
  curve_data <- as.data.frame(curves[[i]]$s[curves[[i]]$ord, ])
  
  # 裁剪轨迹线，确保不超出点的边界范围
  curve_data <- curve_data %>%
    filter(umap_1 >= x_min & umap_1 <= x_max & umap_2 >= y_min & umap_2 <= y_max)
  
  p1 <- p1 + geom_path(
    data = curve_data, 
    aes(x = umap_1, y = umap_2), 
    color = lineage_colors[i],  # 每个分支使用不同颜色
    size = 1.2,  # 轨迹线粗细
    arrow = arrow(type = "closed", length = unit(0.15, "inches"))  # 箭头大小
  )
}
print(p1)

#Calculate the center position of each cell type (for labeling)
cell_type_centers <- plot_data %>%
  group_by(Cluster) %>%
  summarise(
    UMAP1 = median(UMAP1),  # 使用中位数作为中心位置
    UMAP2 = median(UMAP2)
  )

p2 <- ggplot(plot_data, aes(x = UMAP1, y = UMAP2, color = cell_types)) +
  geom_point(size = 1.5, alpha = 0.7) +  # 调整点的大小和透明度
  scale_color_npg() +  # 使用 Nature 风格的配色方案
  theme_classic(base_size = 12) +  # 使用经典主题，设置基础字号
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # 添加实线边框
    axis.line = element_line(size = 0.8),  # 加粗坐标轴线
    axis.title = element_text(size = 14, face = "bold"),  # 坐标轴标题字体
    axis.text = element_text(size = 12, color = "black"),  # 坐标轴刻度字体
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  # 图标题字体
    legend.title = element_text(size = 12, face = "bold"),  # 图例标题字体
    legend.text = element_text(size = 10),  # 图例文字字体
    legend.position = "right"  # 图例位置
  ) +
  labs(title = "Cell Type Classification", x = "UMAP1", y = "UMAP2")

for (i in seq_along(curves)) {
  curve_data <- as.data.frame(curves[[i]]$s[curves[[i]]$ord, ])
  
  # 裁剪轨迹线，确保不超出点的边界范围
  curve_data <- curve_data %>%
    filter(umap_1 >= min(plot_data$UMAP1) & umap_1 <= max(plot_data$UMAP1) &
           umap_2 >= min(plot_data$UMAP2) & umap_2 <= max(plot_data$UMAP2))
  
  p2 <- p2 + geom_path(
    data = curve_data, 
    aes(x = umap_1, y = umap_2), 
    color = "black",  # 轨迹线颜色
    size = 1,  # 轨迹线粗细
    arrow = arrow(type = "closed", length = unit(0.1, "inches"))  # 箭头大小
  )
}

#add cell types annotation
p2 <- p2 + geom_text(
  data = cell_type_centers, 
  aes(x = UMAP1, y = UMAP2, label = Cluster), 
  color = "black",  # 标注颜色
  size = 2.5,  # 标注字体大小
  fontface = "bold",  # 标注字体加粗
  vjust = 1.5,  # 垂直调整位置
  hjust = 1.5   # 水平调整位置
)

print(p2)

#Plotting the pseudotime density of the first trajectory
p3 <- ggplot(plot_data, aes(x = Pseudotime, fill = Cluster)) +
  geom_density(alpha = 0.5, adjust = 4) +  # 绘制密度图，设置透明度, adjust：控制带宽的缩放因子。默认值为 1，增加该值可以平滑密度图
  scale_fill_npg() +  # 使用 Nature 风格的配色方案
  theme_classic(base_size = 12) +  # 使用经典主题，设置基础字号
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # 添加实线边框
    axis.line = element_line(size = 0.8),  # 加粗坐标轴线
    axis.title = element_text(size = 14, face = "bold"),  # 坐标轴标题字体
    axis.text = element_text(size = 12, color = "black"),  # 坐标轴刻度字体
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  # 图标题字体
    legend.title = element_text(size = 12, face = "bold"),  # 图例标题字体
    legend.text = element_text(size = 10),  # 图例文字字体
    legend.position = "right"  # 图例位置
  ) +
  labs(title = "Pseudotime Density Plot", x = "Pseudotime", y = "Density")
print(p3)

#split to 8 pics
p4 <- ggplot(plot_data, aes(x = Pseudotime, fill = Cluster)) +
  geom_density(alpha = 0.5, adjust = 4) +  # 填充颜色，设置透明度，调整带宽
  scale_fill_npg() +  # 使用 Nature 风格的配色方案
  facet_wrap(~ Cluster, ncol = 1) +  # 按细胞类型分面，一列显示
  theme_classic(base_size = 12) +  # 使用经典主题，设置基础字号
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # 添加实线边框
    axis.line = element_line(size = 0.8),  # 加粗坐标轴线
    axis.title = element_text(size = 14, face = "bold"),  # 坐标轴标题字体
    axis.text = element_text(size = 12, color = "black"),  # 坐标轴刻度字体
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  # 图标题字体
    legend.title = element_text(size = 12, face = "bold"),  # 图例标题字体
    legend.text = element_text(size = 10),  # 图例文字字体
    legend.position = "right",  # 图例位置
    axis.text.y = element_blank(),  # 去掉 y 轴刻度
    axis.ticks.y = element_blank(),  # 去掉 y 轴刻度线
    axis.title.y = element_blank(),  # 去掉 y 轴标题
    strip.background = element_blank(),  # 去掉分面标签背景
    strip.text = element_text(size = 12, face = "bold"),  # 分面标签字体
    panel.spacing = unit(0.1, "lines")  # 调整分面图之间的间距，使其更紧凑
  ) +
  labs(title = "Pseudotime Density Plot by Cell Type", x = "Pseudotime")

print(p4)

#enrichment analysis with pseudotime info
#extract cd4 T cells
cd4_cells <- subset(t_cells, idents = grep("^CD4", levels(Idents(t_cells)), value = TRUE))
cd4_cell_names <- colnames(cd4_cells)
#extract Lineage1 pseudotime
cd4_pseudotime <- pseudotime[cd4_cell_names, 1]
#extract hallmark genesets
hallmark_genesets <- read.gmt("../../h.all.v2024.1.Hs.symbols.gmt")
hallmark_pathway_genes <- split(hallmark_genesets$gene, hallmark_genesets$term)

#calculate hallmark pathway scores
hallmark_pathway_scores <- matrix(NA, nrow = length(hallmark_pathway_genes), ncol = ncol(cd4_cells))
rownames(hallmark_pathway_scores) <- names(hallmark_pathway_genes)
colnames(hallmark_pathway_scores) <- colnames(cd4_cells)

for (i in 1:length(hallmark_pathway_genes)) {
  genes <- hallmark_pathway_genes[[i]]
  genes <- genes[genes %in% rownames(cd4_cells)]  # 只保留表达矩阵中存在的基因
  if (length(genes) > 0) {
    hallmark_pathway_scores[i, ] <- colMeans(GetAssayData(cd4_cells, slot = "data")[genes, ])
  }
}

hallmark_pathway_scores <- hallmark_pathway_scores[rowSums(is.na(hallmark_pathway_scores)) == 0, ]

#ordered by pseudoime
ordered_cells <- order(cd4_pseudotime)
hallmark_pathway_scores_ordered <- hallmark_pathway_scores[, ordered_cells]

#z-score standardization
hallmark_pathway_scores_zscore <- t(scale(t(hallmark_pathway_scores_ordered)))

#calculate pathway mean scores & sort
hallmark_pathway_mean_scores <- rowMeans(hallmark_pathway_scores_zscore, na.rm = TRUE)

top_30_hallmark_pathways <- names(sort(hallmark_pathway_mean_scores, decreasing = TRUE))[1:30]
hallmark_pathway_scores_zscore_top30 <- hallmark_pathway_scores_zscore[top_30_hallmark_pathways, ]

#plot
library(pheatmap)

# pseudotime annotation
annotation_col <- data.frame(
  Pseudotime = cd4_pseudotime[ordered_cells]
)
rownames(annotation_col) <- colnames(hallmark_pathway_scores_zscore_top30)

pheatmap(
  hallmark_pathway_scores_zscore_top30,
  cluster_rows = TRUE,  # 对通路进行聚类
  cluster_cols = FALSE, # 不对细胞进行聚类（已按伪时间排序）
  show_colnames = FALSE, # 不显示细胞名称
  annotation_col = annotation_col,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  scale = "none",  # 已经 Z-score 标准化，无需再次标准化
  main = "Top 30 HALLMARK Pathway Score Dynamics Along Pseudotime"
)




# 将伪时间分为四个阶段
cd4_pseudotime_stage <- cut(
    cd4_pseudotime,
    breaks = quantile(cd4_pseudotime, probs = seq(0, 1, length.out = 5), na.rm = TRUE),  # 忽略缺失值
    include.lowest = TRUE,
    labels = c("Stage1", "Stage2", "Stage3", "Stage4")
)
cd4_cells$pseudotime_stage <- cd4_pseudotime_stage

# 设置 ident 为伪时间阶段
Idents(cd4_cells) <- cd4_cells$pseudotime_stage
