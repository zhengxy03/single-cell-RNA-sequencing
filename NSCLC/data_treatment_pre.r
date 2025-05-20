#GSE131907

library(Seurat)
library(dplyr)
library(ggplot2)

#data import
raw_counts <- readRDS("GSE131907_Lung_Cancer_raw_UMI_matrix.rds")
lung_matrix <- raw_counts[, grepl("_LUNG_", colnames(raw_counts))]

expression_matrix <- as.matrix(lung_matrix)
rownames(expression_matrix) <- rownames(lung_matrix)
lung_samples <- grep("_LUNG_", colnames(expression_matrix), value = TRUE)

sample_info <- data.frame(
  cell_id = lung_samples,
  tissue = "LUNG",
  sample_type = gsub(".*_(LUNG|LN|NS)_(.*)", "\\2", lung_samples),
  condition = gsub("([A-Z]+)([0-9]+)", "\\1", gsub(".*_(LUNG|LN|NS)_([A-Z]+[0-9]+)", "\\2", lung_samples)),
  sample_id = gsub("([A-Z]+)([0-9]+)", "\\2", gsub(".*_(LUNG|LN|NS)_([A-Z]+[0-9]+)", "\\2", lung_samples)), 
  stringsAsFactors = FALSE
)
head(sample_info)

seurat_obj <- CreateSeuratObject(
  counts = lung_matrix,
  meta.data = sample_info,
  project = "Lung_Cancer"
)
head(seurat_obj@meta.data)

#normalization and clustering
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 2000)
hvgs <- VariableFeatures(seurat_obj)

seurat_obj <- ScaleData(seurat_obj, features = hvgs)
seurat_obj <- RunPCA(seurat_obj, features = hvgs, npcs = 20)

seurat_obj <- FindNeighbors(seurat_obj, dims = 1:15)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)
DimPlot(seurat_obj, reduction = "umap", label = TRUE)

#cell type annotation
markers <- FindAllMarkers(seurat_obj, 
                          only.pos = TRUE, 
                          min.pct = 0.25, 
                          logfc.threshold = 0.25, 
                          test.use = "wilcox")

significant_markers <- subset(markers, p_val_adj < 0.05)
significant_markers <- significant_markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.csv(significant_markers,"marker_top.csv")

identity_mapping <- c(
    "0" = "Effector T cell",
    "1" = "NK cell",
    "2" = "Macrophage",
    "3" = "Naive T cell",
    "4" = "Myeloid cell",
    "5" = "B cell",
    "6" = "Alveolar type II cell",
    "7" = "Neutrophil",
    "8" = "Fibroblast",
    "9" = "Monocyte",
    "10" = "Mast cell",
    "11" = "Basal Epithelial cell",
    "12" = "Ciliated Epithelial cell_1",
    "13" = "Dendritic cell",
    "14" = "Endothelial cell",
    "15" = "Plasma cell",
    "16" = "Alveolar type I cell",
    "17" = "Proliferating cell_1",
    "18" = "Proliferating cell_2",
    "19" = "Ciliated Epithelial cell_2",
    "20" = "Mesothelial cell"
)
cell_type <- identity_mapping[seurat_obj@meta.data$seurat_clusters]
seurat_obj@meta.data$cell_type <- cell_type
cell_type_annotation <- DimPlot(seurat_obj, reduction = "umap", label = FALSE, group.by = "cell_type")
ggsave("celltype_annotation.png", plot = cell_type_annotation, width = 8, height = 6, dpi = 300)


#tumor
tumor <- subset(seurat_obj, subset = condition == "Tumor")
tumor_celltype <- DimPlot(tumor, reduction = "umap", label = FALSE, group.by = "cell_type")
ggsave("tumor_celltype.png", plot = tumor_celltype, width = 8, height = 4, dpi = 300)

#normal
normal <- subset(seurat_obj, subset = condition == "Normal")
normal_celltype <- DimPlot(normal, reduction = "umap", label = FALSE, group.by = "cell_type")
ggsave("normal_celltype.png", plot = normal_celltype, width = 8, height = 4, dpi = 300)

#featureplot
target_genes <- c("CDKN2A", "FZD10", "NOTCH1", "PDGFRA", "WNT7B")
tumor_feature_plot <- FeaturePlot(tumor, features = target_genes)
ggsave("tumor_featureplot.png", plot = tumor_feature_plot, width = 6 , height = 8, dpi = 300)
normal_feature_plot <- FeaturePlot(normal, features = target_genes)
ggsave("normal_featureplot.png", plot = normal_feature_plot, width = 6 , height = 8, dpi = 300)