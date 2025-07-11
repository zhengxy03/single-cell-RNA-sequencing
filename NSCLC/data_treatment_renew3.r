#GSE131907

library(Seurat)
library(dplyr)
library(ggplot2)
library(DESeq2)
library(edgeR)
library(limma)
library(tidyverse)
#data import
raw_counts <- readRDS("GSE131907_Lung_Cancer_raw_UMI_matrix.rds")

sample_info <- read.table(
  file = "GSE131907_Lung_Cancer_cell_annotation.txt",
  sep = "\t",
  header = TRUE,
  fill = TRUE,
  stringsAsFactors = FALSE
)

rownames(sample_info) <- sample_info$Index

seurat_obj <- CreateSeuratObject(
  counts = raw_counts,
  meta.data = sample_info,
  project = "Lung_Cancer"
)

str(seurat_obj)
head(seurat_obj@meta.data)



#lung
lung <- subset(seurat_obj, subset = Sample_Origin %in% c("nLung", "tLung"))
lung@meta.data$patient_id <- gsub("LUNG_([TN][0-9]+)", "P\\1", lung@meta.data$Sample)
lung@meta.data$patient_id <- gsub("T", "", lung@meta.data$patient_id)
lung@meta.data$patient_id <- gsub("N", "", lung@meta.data$patient_id)
lung_paired <- subset(lung, subset = patient_id %in% c("P06", "P08", "P09", "P18", "P19", "P20", "P28", "P30", "P31", "P34"))
# 过滤低表达基因（min.cells < 0.1%总细胞数）
min_cells_threshold <- ncol(lung_paired) * 0.001  # 0.1%的细胞数
lung_paired <- subset(
  lung_paired,
  features = rownames(lung_paired)[Matrix::rowSums(GetAssayData(lung_paired, assay = "RNA", layer = "counts") > 0) >= min_cells_threshold]
)

lung_paired@meta.data$Sample_Origin <- ifelse(
  lung_paired@meta.data$Sample_Origin == "nLung", "Normal",
  ifelse(lung_paired@meta.data$Sample_Origin == "tLung", "Tumor", lung_paired@meta.data$Sample_Origin)
)
lung_paired@meta.data$cell_type_pre <- lung_paired@meta.data$cell_type

lung_paired <- NormalizeData(
  lung_paired,
  normalization.method = "LogNormalize",
  scale.factor = 1e6,
  assay = "RNA"
)

lung_paired <- FindVariableFeatures(
  lung_paired,
  selection.method = "vst",
  nfeatures = 2000,
  mean.cutoff = c(0.0125, 3),
  dispersion.cutoff = c(0.5, Inf)
)
hvgs <- VariableFeatures(lung_paired)

lung_paired <- ScaleData(
  lung_paired,
  features = hvgs,
  do.scale = FALSE,
  do.center = TRUE,
  scale.max = 10
)

lung_paired <- RunPCA(
  lung_paired,
  features = hvgs,
  npcs = 50,
  verbose = TRUE
)

#tsne
n_significant_pcs <- 15

lung_paired <- RunTSNE(
  lung_paired,
  dims = 1:n_significant_pcs,
  perplexity = 30,
  verbose = TRUE
)

#umap
lung_paired <- RunUMAP(
  lung_paired,
  dims = 1:n_significant_pcs,
  verbose = TRUE
)

lung_paired <- FindNeighbors(lung_paired, dims = 1:15)
lung_paired <- FindClusters(lung_paired, resolution = 0.5)
DimPlot(lung_paired, reduction = "umap", label = TRUE)

umap_markers <- FindAllMarkers(lung_paired, 
                              only.pos = TRUE, 
                              min.pct = 0.25, 
                              logfc.threshold = 0.25, 
                              test.use = "wilcox")

umap_significant_markers <- subset(umap_markers, p_val_adj < 0.05)
umap_significant_markers <- umap_significant_markers %>% 
  group_by(cluster) %>% 
  top_n(n = 20, wt = avg_log2FC)
write.csv(umap_significant_markers, "umap_marker_top.csv")


identity_mapping <- c(
    "0" = "Naive T cell",
    "1" = "NK cell",
    "2" = "Alveolar Macrophage",
    "3" = "Effector T cell",
    "4" = "Transitional T cell",
    "5" = "Tumor-associated Macrophage cell",
    "6" = "Macrophage/Monocyte",
    "7" = "Alveolar type II cell",
    "8" = "Fibroblast",
    "9" = "Monocyte",
    "10" = "B cell",
    "11" = "Mast cell",
    "12" = "Ciliated Epithelial cell_1",
    "13" = "Conventional dendritic cell",
    "14" = "Endothelial cell",
    "15" = "Plasma cell",
    "16" = "Basal Epithelial cell_1",
    "17" = "Adenocarcinoma stem-like cell",
    "18" = "Ciliated Epithelial cell_2",
    "19" = "Alveolar type I cell",
    "20" = "Proliferating cell",
    "21" = "Basal Epithelial cell_2",
    "22" = "Plasmacytoid dendritic cell",
    "23" = "Ciliated Epithelial cell_3"
)
cell_type <- identity_mapping[lung_paired@meta.data$seurat_clusters]
lung_paired@meta.data$cell_type <- cell_type

lung_paired@meta.data$cell_type <- factor(
  lung_paired@meta.data$cell_type,
  levels = sort(unique(lung_paired@meta.data$cell_type))
)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(24)
cell_types <- sort(unique(lung_paired@meta.data$cell_type))
num_types <- length(cell_types)
color_map <- setNames(npg_extended[1:num_types], cell_types)


#tsne-annoation
tsne_coords <- Embeddings(lung_paired, reduction = "tsne")

tsne_data <- data.frame(
  tSNE_1 = tsne_coords[, 1],
  tSNE_2 = tsne_coords[, 2],
  cell_type = lung_paired$cell_type
)


p_tsne <- ggplot(tsne_data, aes(x = tSNE_1, y = tSNE_2, color = cell_type)) +
  geom_point(size = 0.001, shape =16) + 
  scale_color_manual(values = color_map) +
  labs(x="tSNE 1",y="tSNE 2") +
  theme(aspect.ratio = 1,
    plot.title = element_blank(),
    axis.text = element_text(colour = "#000000", face = "plain", size = 5.5),
    axis.title.y = element_text(colour = "#000000", face = "bold", size = 7),
    axis.title.x = element_text(colour = "#000000", face = "bold", size = 7),
    axis.ticks.length = unit(0.35, units = "mm"),
    axis.ticks = element_line(colour = "#000000", linewidth = 0.1167156),
    axis.line = element_line(colour = "#000000", linewidth = 0.1167156),
    legend.title = element_blank(),
    legend.margin = margin(t = 0,r = 0,b = 0,l = 0,unit = "cm"),
    legend.position = "right",
    legend.spacing = unit(-0.2, "cm"),
    legend.background = element_blank(),
    legend.key.spacing = unit(-0.2, "cm"),
    legend.key.size = unit(0.3, "cm"),
    legend.spacing.y = unit(-0.2, "cm"),
    legend.text = element_text(size = 5,margin = margin(r = 0.1, unit = "cm")),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank()
  )


pdf(
  file = "celltype_all_tsne.pdf",
  family = "Helvetica",
  width = 5.6,
  height = 3.2,
  useDingbats = FALSE
)
print(p_tsne)
dev.off()

#featureplot-tsne
tumor <- subset(lung_paired, subset = Sample_Origin == "Tumor")
normal <- subset(lung_paired, subset = Sample_Origin == "Normal")

target_genes <- c("CDKN2A", "NOTCH1", "PDGFRA", "WNT7B")
#tumor
pdf(file = "gene_feature_plots_tumor_tsne.pdf", 
    width = 3.2, 
    height = 3.2,
    family = "Helvetica",
    useDingbats = FALSE
)

tsne_coords <- Embeddings(tumor, reduction = "tsne") %>%
  as.data.frame() %>%
  setNames(c("tSNE_1", "tSNE_2")) %>%
  rownames_to_column("cell_id")

expr_matrix <- GetAssayData(tumor, slot = "data")

for (gene in target_genes) {
  gene_expr <- expr_matrix[gene, ] %>%
    as.data.frame() %>%
    rownames_to_column("cell_id") %>%
    setNames(c("cell_id", "expression"))

  plot_data <- left_join(tsne_coords, gene_expr, by = "cell_id") %>%
    mutate(group = ifelse(expression > 0, "expressed", "not_expressed")) %>%
    arrange(expression)
  
  p <- ggplot(plot_data, aes(x = tSNE_1, y = tSNE_2)) +
    geom_point(
      data = subset(plot_data, group == "not_expressed"),
      aes(color = expression),
      size = 0.001, 
      shape = 16
    ) +
    geom_point(
      data = subset(plot_data, group == "expressed"),
      aes(color = expression),
      size = 0.001, 
      shape = 16
    ) +
    scale_color_gradientn(
      colors = c("lightgrey", "blue")
    ) +
    ggtitle(gene) +
    labs(x = "tSNE 1", y = "tSNE 2") +
    theme(
      aspect.ratio = 1,
      plot.title = element_text(colour = "#000000", face = "bold", size = 7, hjust = 0.5),
      axis.text = element_text(colour = "#000000", face = "plain", size = 5.5),
      axis.title.y = element_text(colour = "#000000", face = "bold", size = 7),
      axis.title.x = element_text(colour = "#000000", face = "bold", size = 7),
      axis.ticks.length = unit(0.35, units = "mm"),
      axis.ticks = element_line(colour = "#000000", linewidth = 0.1167156),
      axis.line = element_line(colour = "#000000", linewidth = 0.1167156),
      legend.title = element_blank(),
      legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
      legend.position = "right",
      legend.spacing = unit(-0.2, "cm"),
      legend.background = element_blank(),
      legend.key.spacing = unit(-0.2, "cm"),
      legend.key.size = unit(0.3, "cm"),
      legend.spacing.y = unit(-0.2, "cm"),
      legend.text = element_text(size = 5, margin = margin(r = 0.1, unit = "cm")),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank()
    )
  print(p)
}
dev.off()
#normal
pdf(file = "gene_feature_plots_normal_tsne.pdf", 
    width = 3.2, 
    height = 3.2,
    family = "Helvetica",
    useDingbats = FALSE
)

tsne_coords <- Embeddings(normal, reduction = "tsne") %>%
  as.data.frame() %>%
  setNames(c("tSNE_1", "tSNE_2")) %>%
  rownames_to_column("cell_id")

expr_matrix <- GetAssayData(normal, slot = "data")

for (gene in target_genes) {
  gene_expr <- expr_matrix[gene, ] %>%
    as.data.frame() %>%
    rownames_to_column("cell_id") %>%
    setNames(c("cell_id", "expression"))

  plot_data <- left_join(tsne_coords, gene_expr, by = "cell_id") %>%
    mutate(group = ifelse(expression > 0, "expressed", "not_expressed")) %>%
    arrange(expression)
  
  p <- ggplot(plot_data, aes(x = tSNE_1, y = tSNE_2)) +
    geom_point(
      data = subset(plot_data, group == "not_expressed"),
      aes(color = expression),
      size = 0.001, 
      shape = 16
    ) +
    geom_point(
      data = subset(plot_data, group == "expressed"),
      aes(color = expression),
      size = 0.001, 
      shape = 16
    ) +
    scale_color_gradientn(
      colors = c("lightgrey", "blue")
    ) +
    ggtitle(gene) +
    labs(x = "tSNE 1", y = "tSNE 2") +
    theme(
      aspect.ratio = 1,
      plot.title = element_text(colour = "#000000", face = "bold", size = 7, hjust = 0.5),
      axis.text = element_text(colour = "#000000", face = "plain", size = 5.5),
      axis.title.y = element_text(colour = "#000000", face = "bold", size = 7),
      axis.title.x = element_text(colour = "#000000", face = "bold", size = 7),
      axis.ticks.length = unit(0.35, units = "mm"),
      axis.ticks = element_line(colour = "#000000", linewidth = 0.1167156),
      axis.line = element_line(colour = "#000000", linewidth = 0.1167156),
      legend.title = element_blank(),
      legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
      legend.position = "right",
      legend.spacing = unit(-0.2, "cm"),
      legend.background = element_blank(),
      legend.key.spacing = unit(-0.2, "cm"),
      legend.key.size = unit(0.3, "cm"),
      legend.spacing.y = unit(-0.2, "cm"),
      legend.text = element_text(size = 5, margin = margin(r = 0.1, unit = "cm")),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank()
    )
  print(p)
}
dev.off()



npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(24)
Cell_types <- sort(unique(lung_paired@meta.data$Cell_type))
num_types <- length(Cell_types)
color_map <- setNames(npg_pal[1:num_types], Cell_types)
#uamp-annotation
umap_coords <- Embeddings(lung_paired, reduction = "umap")

umap_data <- data.frame(
  UMAP_1 = umap_coords[, 1],
  UMAP_2 = umap_coords[, 2],
  cell_type = lung_paired$Cell_type
)

p_umap <- ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, color = cell_type)) +
  geom_point(size = 0.001, shape = 16) +
  scale_color_manual(values = color_map) +
  guides(color = guide_legend(override.aes = list(size = 2))) + 
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme(
    aspect.ratio = 1,
    plot.title = element_blank(),
    axis.text = element_text(colour = "#000000", face = "plain", size = 5.5),
    axis.title.y = element_text(colour = "#000000", face = "bold", size = 7),
    axis.title.x = element_text(colour = "#000000", face = "bold", size = 7),
    axis.ticks.length = unit(0.35, units = "mm"),
    axis.ticks = element_line(colour = "#000000", linewidth = 0.1167156),
    axis.line = element_line(colour = "#000000", linewidth = 0.1167156),
    legend.title = element_blank(),
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
    legend.position = "right",
    legend.spacing = unit(-0.2, "cm"),
    legend.background = element_blank(),
    legend.key.spacing = unit(-0.2, "cm"),
    legend.key.size = unit(0.3, "cm"),
    legend.spacing.y = unit(-0.2, "cm"),
    legend.text = element_text(size = 5, margin = margin(r = 0.1, unit = "cm")),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank()
  )

pdf(
  file = "Celltype_all_umap.pdf",
  family = "Helvetica",
  width = 5.6,
  height = 3.2,
  useDingbats = FALSE
)
print(p_umap)
dev.off()

#featureplot-umap
tumor <- subset(lung_paired, subset = Sample_Origin == "Tumor")
normal <- subset(lung_paired, subset = Sample_Origin == "Normal")

target_genes <- c("CDKN2A", "NOTCH1", "PDGFRA", "WNT7B")

#tumor
pdf(file = "gene_feature_plots_tumor_umap.pdf", 
    width = 3.2, 
    height = 3.2,
    family = "Helvetica",
    useDingbats = FALSE
)

umap_coords <- Embeddings(tumor, reduction = "umap") %>%
  as.data.frame() %>%
  setNames(c("UMAP_1", "UMAP_2")) %>% 
  rownames_to_column("cell_id")

expr_matrix <- GetAssayData(tumor, slot = "data")

for (gene in target_genes) {
  gene_expr <- expr_matrix[gene, ] %>%
    as.data.frame() %>%
    rownames_to_column("cell_id") %>%
    setNames(c("cell_id", "expression"))

  plot_data <- left_join(umap_coords, gene_expr, by = "cell_id") %>%
    mutate(group = ifelse(expression > 0, "expressed", "not_expressed")) %>%
    arrange(expression)
  
  p <- ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2)) + 
    geom_point(
      data = subset(plot_data, group == "not_expressed"),
      aes(color = expression),
      size = 0.001,
      shape = 16
    ) +
    geom_point(
      data = subset(plot_data, group == "expressed"),
      aes(color = expression),
      size = 0.001, 
      shape = 16
    ) +
    scale_color_gradientn(
      colors = c("lightgrey", "blue")
    ) +
    ggtitle(gene) +
    labs(x = "UMAP 1", y = "UMAP 2") +
    theme(
      aspect.ratio = 1,
      plot.title = element_text(colour = "#000000", face = "bold", size = 7, hjust = 0.5),
      axis.text = element_text(colour = "#000000", face = "plain", size = 5.5),
      axis.title.y = element_text(colour = "#000000", face = "bold", size = 7),
      axis.title.x = element_text(colour = "#000000", face = "bold", size = 7),
      axis.ticks.length = unit(0.35, units = "mm"),
      axis.ticks = element_line(colour = "#000000", linewidth = 0.1167156),
      axis.line = element_line(colour = "#000000", linewidth = 0.1167156),
      legend.title = element_blank(),
      legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
      legend.position = "right",
      legend.spacing = unit(-0.2, "cm"),
      legend.background = element_blank(),
      legend.key.spacing = unit(-0.2, "cm"),
      legend.key.size = unit(0.3, "cm"),
      legend.spacing.y = unit(-0.2, "cm"),
      legend.text = element_text(size = 5, margin = margin(r = 0.1, unit = "cm")),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank()
    )
  print(p)
}
dev.off()
#normal
pdf(file = "gene_feature_plots_normal_umap.pdf", 
    width = 3.2, 
    height = 3.2,
    family = "Helvetica",
    useDingbats = FALSE
)

umap_coords <- Embeddings(normal, reduction = "umap") %>%
  as.data.frame() %>%
  setNames(c("UMAP_1", "UMAP_2")) %>%
  rownames_to_column("cell_id")

expr_matrix <- GetAssayData(normal, slot = "data")

for (gene in target_genes) {
  gene_expr <- expr_matrix[gene, ] %>%
    as.data.frame() %>%
    rownames_to_column("cell_id") %>%
    setNames(c("cell_id", "expression"))

  plot_data <- left_join(umap_coords, gene_expr, by = "cell_id") %>%
    mutate(group = ifelse(expression > 0, "expressed", "not_expressed")) %>%
    arrange(expression)
  
  p <- ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(
      data = subset(plot_data, group == "not_expressed"),
      aes(color = expression),
      size = 0.001,
      shape = 16
    ) +
    geom_point(
      data = subset(plot_data, group == "expressed"),
      aes(color = expression),
      size = 0.001,
      shape = 16
    ) +
    scale_color_gradientn(
      colors = c("lightgrey", "blue")
    ) +
    ggtitle(gene) +
    labs(x = "UMAP 1", y = "UMAP 2") +
    theme(
      aspect.ratio = 1,
      plot.title = element_text(colour = "#000000", face = "bold", size = 7, hjust = 0.5),
      axis.text = element_text(colour = "#000000", face = "plain", size = 5.5),
      axis.title.y = element_text(colour = "#000000", face = "bold", size = 7),
      axis.title.x = element_text(colour = "#000000", face = "bold", size = 7),
      axis.ticks.length = unit(0.35, units = "mm"),
      axis.ticks = element_line(colour = "#000000", linewidth = 0.1167156),
      axis.line = element_line(colour = "#000000", linewidth = 0.1167156),
      legend.title = element_blank(),
      legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
      legend.position = "right",
      legend.spacing = unit(-0.2, "cm"),
      legend.background = element_blank(),
      legend.key.spacing = unit(-0.2, "cm"),
      legend.key.size = unit(0.3, "cm"),
      legend.spacing.y = unit(-0.2, "cm"),
      legend.text = element_text(size = 5, margin = margin(r = 0.1, unit = "cm")),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank()
    )
  print(p)
}
dev.off()



#distribution plot

lung_paired@meta.data$combined_group <- paste0(
  lung_paired@meta.data$patient_id, "_", 
  lung_paired@meta.data$Sample_Origin
)
lung_paired@meta.data$Sample_Origin <- factor(
  lung_paired@meta.data$Sample_Origin,
  levels = c("Tumor", "Normal"),
  ordered = TRUE
)
unique_patients <- sort(unique(lung_paired@meta.data$patient_id))
unique_sample_types <- levels(lung_paired@meta.data$Sample_Origin)
#tsne
tsne_coords <- Embeddings(lung_paired, reduction = "tsne")
tsne_data <- data.frame(
  tSNE_1 = tsne_coords[, 1],
  tSNE_2 = tsne_coords[, 2],
  cell_type = lung_paired$cell_type,
  patient_id = lung_paired$patient_id,
  Sample_Origin = lung_paired$Sample_Origin,
  combined_group = lung_paired$combined_group
)

for (sample in unique_sample_types) {
  pdf_file <- paste0("patient_celltype_tsne_", sample, ".pdf")
  
  pdf(file = pdf_file, 
      width = 5.6, 
      height = 3.2,
      family = "Helvetica",
      useDingbats = FALSE
  )
  
  for (patient in unique_patients) {
    patient_subset_data <- subset(
      tsne_data, 
      combined_group == paste0(patient, "_", sample)
    )
    
    if (nrow(patient_subset_data) == 0) next

    p <- ggplot(patient_subset_data, aes(x = tSNE_1, y = tSNE_2, color = cell_type)) +
      geom_point(size = 0.001, shape = 16) + 
      scale_color_manual(values = color_map) +
      labs(x = "tSNE 1", y = "tSNE 2") +
      ggtitle(paste(patient, sample, sep = " - ")) +
      guides(color=guide_legend(ncol=2)) +
      theme(
        aspect.ratio = 1,
        plot.title = element_text(
          colour = "#000000", 
          face = "bold", 
          size = 7, 
          hjust = 0.5
        ),
        axis.text = element_text(colour = "#000000", face = "plain", size = 5.5),
        axis.title.y = element_text(colour = "#000000", face = "bold", size = 7),
        axis.title.x = element_text(colour = "#000000", face = "bold", size = 7),
        axis.ticks.length = unit(0.35, units = "mm"),
        axis.ticks = element_line(colour = "#000000", linewidth = 0.1167156),
        axis.line = element_line(colour = "#000000", linewidth = 0.1167156),
        legend.title = element_blank(),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
        legend.position = "right",
        legend.spacing = unit(-0.2, "cm"),
        legend.background = element_blank(),
        legend.key.spacing = unit(-0.2, "cm"),
        legend.key.size = unit(0.3, "cm"),
        legend.spacing.y = unit(-0.2, "cm"),
        legend.text = element_text(size = 5, margin = margin(r = 0.1, unit = "cm")),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()
      )
    print(p)
  }
  dev.off()
}

#umap
umap_coords <- Embeddings(lung_paired, reduction = "umap")
umap_data <- data.frame(
  UMAP_1 = umap_coords[, 1],
  UMAP_2 = umap_coords[, 2],
  cell_type = lung_paired$cell_type,
  patient_id = lung_paired$patient_id,
  Sample_Origin = lung_paired$Sample_Origin,
  combined_group = lung_paired$combined_group
)

for (sample in unique_sample_types) {
  pdf_file <- paste0("patient_celltype_umap_", sample, ".pdf")
  
  pdf(file = pdf_file, 
      width = 5.6, 
      height = 3.2,
      family = "Helvetica",
      useDingbats = FALSE
  )
  
  for (patient in unique_patients) {
    patient_subset_data <- subset(
      umap_data, 
      combined_group == paste0(patient, "_", sample)
    )
    
    if (nrow(patient_subset_data) == 0) next

    p <- ggplot(patient_subset_data, aes(x = UMAP_1, y = UMAP_2, color = cell_type)) +
      geom_point(size = 0.001, shape = 16) + 
      scale_color_manual(values = color_map) +
      labs(x = "UMAP 1", y = "UMAP 2") +
      ggtitle(paste(patient, sample, sep = " - ")) +
      guides(color = guide_legend(ncol = 2)) +
      theme(
        aspect.ratio = 1,
        plot.title = element_text(
          colour = "#000000", 
          face = "bold", 
          size = 7, 
          hjust = 0.5
        ),
        axis.text = element_text(colour = "#000000", face = "plain", size = 5.5),
        axis.title.y = element_text(colour = "#000000", face = "bold", size = 7),
        axis.title.x = element_text(colour = "#000000", face = "bold", size = 7),
        axis.ticks.length = unit(0.35, units = "mm"),
        axis.ticks = element_line(colour = "#000000", linewidth = 0.1167156),
        axis.line = element_line(colour = "#000000", linewidth = 0.1167156),
        legend.title = element_blank(),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
        legend.position = "right",
        legend.spacing = unit(-0.2, "cm"),
        legend.background = element_blank(),
        legend.key.spacing = unit(-0.2, "cm"),
        legend.key.size = unit(0.3, "cm"),
        legend.spacing.y = unit(-0.2, "cm"),
        legend.text = element_text(size = 5, margin = margin(r = 0.1, unit = "cm")),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()
      )
    print(p)
  }
  dev.off()
}

#proportion
proportion_data <- lung_paired@meta.data %>%
  group_by(patient_id, Sample_Origin, cell_type) %>% 
  summarise(count = n()) %>% 
  mutate(proportion = count / sum(count)) %>%
  ungroup()

patients <- unique(proportion_data$patient_id)
fixed_cell_types <- sort(unique(proportion_data$cell_type))


pdf("proportion_all_patients_combined_CellType.pdf", width = 30, height = 10, family = "Helvetica")

for (patient in patients) {
  patient_data <- proportion_data %>% filter(patient_id == patient)
  sample_types <- unique(patient_data$Sample_Origin)
  pie_plots <- list()
  
  for (sample in sample_types) {
    sample_data <- patient_data %>% filter(Sample_Origin == sample)
    p <- ggplot(sample_data, aes(x = "", y = proportion, fill = cell_type)) +
      geom_bar(stat = "identity", width = 1) +
      coord_polar(theta = "y") +
      theme_void() +
      labs(title = sample, fill = "Cell Type") +
      scale_fill_manual(values = color_map) +
      theme(
        plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
        legend.position = "right",
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 22),
      ) + guides(fill = guide_legend(ncol = 1))
    pie_plots[[sample]] <- p
  }
  
  combined_pie <- wrap_plots(pie_plots, ncol = length(sample_types))
  
  bar_chart <- ggplot(patient_data, aes(x = Sample_Origin, y = proportion, fill = cell_type)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(x = "Sample Type", y = "Proportion", fill = "Cell Type") +
    theme_classic() +
    scale_fill_manual(values = color_map) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
      axis.title.x = element_text(size = 24),
      axis.title.y = element_text(size = 24),
      axis.text.y = element_text(size = 20),
      legend.position = "none"
    )
  
  patient_combined <- combined_pie | bar_chart + 
    plot_layout(widths = c(2, 1)) +
    plot_annotation(title = paste("Patient", patient, "- Cell Type Proportions"),
                    theme = theme(plot.title = element_text(size = 24, face = "bold", hjust = 0.5)))
  
  print(patient_combined)
}

dev.off()

proportion_data <- lung_paired@meta.data %>%
  group_by(Sample_Origin, cell_type) %>% 
  summarise(count = n()) %>% 
  mutate(proportion = count / sum(count)) %>%
  ungroup()

patients <- unique(proportion_data$patient_id)
fixed_cell_types <- sort(unique(proportion_data$cell_type))


pdf("proportion_combined_CellType.pdf", width = 30, height = 10, family = "Helvetica")

for (patient in patients) {
  patient_data <- proportion_data %>% filter(patient_id == patient)
  sample_types <- unique(patient_data$Sample_Origin)
  pie_plots <- list()
  
  for (sample in sample_types) {
    sample_data <- patient_data %>% filter(Sample_Origin == sample)
    p <- ggplot(sample_data, aes(x = "", y = proportion, fill = cell_type)) +
      geom_bar(stat = "identity", width = 1) +
      coord_polar(theta = "y") +
      theme_void() +
      labs(title = sample, fill = "Cell Type") +
      scale_fill_manual(values = color_map) +
      theme(
        plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
        legend.position = "right",
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 22),
      ) + guides(fill = guide_legend(ncol = 1))
    pie_plots[[sample]] <- p
  }
  
  combined_pie <- wrap_plots(pie_plots, ncol = length(sample_types))
  
  bar_chart <- ggplot(patient_data, aes(x = Sample_Origin, y = proportion, fill = cell_type)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(x = "Sample Type", y = "Proportion", fill = "Cell Type") +
    theme_classic() +
    scale_fill_manual(values = color_map) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
      axis.title.x = element_text(size = 24),
      axis.title.y = element_text(size = 24),
      axis.text.y = element_text(size = 20),
      legend.position = "none"
    )
  
  patient_combined <- combined_pie | bar_chart + 
    plot_layout(widths = c(2, 1)) +
    plot_annotation(title = paste("Patient", patient, "- Cell Type Proportions"),
                    theme = theme(plot.title = element_text(size = 24, face = "bold", hjust = 0.5)))
  
  print(patient_combined)
}

dev.off()

#pseudo bulk
meta <- lung_paired@meta.data %>%
  as_tibble(rownames = "cell") %>%
  mutate(Patient = patient_id, Condition = Sample_Origin) %>%
  select(cell, Patient, Condition)

counts_matrix <- GetAssayData(lung_paired, assay = "RNA", layer = "counts")

meta$sample_id <- paste0(meta$Patient, "_", meta$Condition)
mat_t <- t(counts_matrix)
pb_t <- rowsum(mat_t, group = meta$sample_id)
pb_counts <- t(pb_t)


coldata <- tibble(sample_id = colnames(pb_counts)) %>%
  separate(
    sample_id,
    into = c("Patient", "Condition"),
    sep = "_",
    remove = FALSE
  ) %>%
  mutate(
    Condition = factor(Condition, levels = c("Normal", "Tumor")),
    Patient = factor(Patient)
  ) %>%
  column_to_rownames("sample_id")

dds <- DESeqDataSetFromMatrix(
  countData = pb_counts,
  colData   = coldata,
  design    = ~ Patient + Condition
)
dds <- DESeq(dds)

res <- results(dds, contrast = c("Condition", "Tumor", "Normal"))
resultsNames(dds)
res <- lfcShrink(dds, coef = "Condition_Tumor_vs_Normal", type = "apeglm")

res_tbl <- as_tibble(res, rownames = "gene") %>% arrange(padj)

gene_list <- c("CDKN2A", "FZD10", "NOTCH1", "PDGFRA", "WNT7B")
res_focus <- res_tbl %>%
  filter(gene %in% gene_list) %>%
  select(gene, log2FoldChange, pvalue, padj)



#delete each cell type
cell_types <- unique(lung_paired@meta.data$cell_type)
all_results <- list()

for (ct in cell_types) {
  message("Processing without cell type: ", ct)
  

  lung_paired_filtered <- subset(lung_paired, subset = Cell_type != ct)
  
  meta <- lung_paired_filtered@meta.data %>%
    as_tibble(rownames = "cell") %>%
    mutate(Patient = patient_id, Condition = Sample_Origin) %>%
    select(cell, Patient, Condition)

  counts_matrix <- GetAssayData(lung_paired_filtered, assay = "RNA", layer = "counts")

  meta$sample_id <- paste0(meta$Patient, "_", meta$Condition)
  mat_t <- t(counts_matrix)
  pb_t <- rowsum(mat_t, group = meta$sample_id)
  pb_counts <- t(pb_t)

  coldata <- tibble(sample_id = colnames(pb_counts)) %>%
    separate(
      sample_id,
      into = c("Patient", "Condition"),
      sep = "_",
      remove = FALSE
    ) %>%
    mutate(
      Condition = factor(Condition, levels = c("Normal", "Tumor")),
      Patient = factor(Patient)
    ) %>%
    column_to_rownames("sample_id")

  dds <- DESeqDataSetFromMatrix(
    countData = pb_counts,
    colData   = coldata,
    design    = ~ Patient + Condition
  )
  dds <- DESeq(dds)

  res <- results(dds, contrast = c("Condition", "Tumor", "Normal"))
  res <- lfcShrink(dds, coef = "Condition_Tumor_vs_Normal", type = "apeglm")

  res_tbl <- as_tibble(res, rownames = "gene") %>% arrange(padj)
  

  res_focus <- res_tbl %>%
    filter(gene %in% gene_list) %>%
    select(gene, log2FoldChange, pvalue, padj) %>%
    mutate(Removed_Cell_Type = ct)
  

  all_results[[ct]] <- res_focus
}


combined_results <- bind_rows(all_results)
write.csv(combined_results, "differential_expression_by_cell_type_6.9.csv", row.names = FALSE)