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

seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 2000)
hvgs <- VariableFeatures(seurat_obj)

seurat_obj <- ScaleData(seurat_obj, features = hvgs)
seurat_obj <- RunPCA(seurat_obj, features = hvgs, npcs = 20)
seurat_obj <- RunTSNE(seurat_obj, dims = 1:10)

saveRDS(seurat_obj, file = "GSE131907_all.rds")

seurat_obj <- FindNeighbors(seurat_obj, dims = 1:15)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
p1 <- DimPlot(seurat_obj, reduction = "tsne", label = TRUE)
ggsave("clusters.png", plot = p1, width = 6, height = 4,dpi = 300)

markers <- FindAllMarkers(seurat_obj, 
                          only.pos = TRUE, 
                          min.pct = 0.25, 
                          logfc.threshold = 0.25, 
                          test.use = "wilcox")

significant_markers <- subset(markers, p_val_adj < 0.05)
significant_markers <- significant_markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.csv(significant_markers,"marker_top.csv")

identity_mapping <- c(
    "0" = "NK cell",
    "1" = "Naive T cell",
    "2" = "B cell",
    "3" = "Monocyte",
    "4" = "Macrophage",
    "5" = "Macrophage/Monocyte",
    "6" = "Alveolar type II cell_1",
    "7" = "Basal Epithelial cell_1",
    "8" = "Mesenchymal cell",
    "9" = "Squamous Epithelial cell",
    "10" = "Fibroblast",
    "11" = "Mast cell",
    "12" = "Plasma cell",
    "13" = "Nerve cell",
    "14" = "Basal Epithelial cell_2",
    "15" = "Proliferating cell",
    "16" = "Vascular Endothelial cell",
    "17" = "Airway Secretory cell",
    "18" = "Effector T cell",
    "19" = "Proliferating cell",
    "20" = "Club cell",
    "21" = "Circulating fetal cell",
    "22" = "Plasmacytoid Dendritic cell",
    "23" = "Brush cell",
    "24" = "Ciliated cell",
    "25" = "Oligodendrocyte",
    "26" = "Alveolar type I cell",
    "27" = "Tumor-associated Macrophage",
    "28" = "Dorsal Neuronal cell",
    "29" = "Mesothelial cells",
    "30" = "Alveolar type II cell_2",
    "31" = "Basal Epithelial cell_3"
)
cell_type <- identity_mapping[seurat_obj@meta.data$seurat_clusters]
seurat_obj@meta.data$cell_type <- cell_type
cell_type_annotation <- DimPlot(seurat_obj, reduction = "tsne", label = FALSE, group.by = "cell_type")
ggsave("celltype_annotation.png", plot = cell_type_annotation, width = 9, height = 4, dpi = 300)

#lung
lung <- subset(seurat_obj, subset = Sample_Origin %in% c("nLung", "tLung"))
lung@meta.data$patient_id <- gsub("LUNG_([TN][0-9]+)", "P\\1", lung@meta.data$Sample)
lung@meta.data$patient_id <- gsub("T", "", lung@meta.data$patient_id)
lung@meta.data$patient_id <- gsub("N", "", lung@meta.data$patient_id)
lung_paired <- subset(lung, subset = patient_id %in% c("P06", "P08", "P09", "P18", "P19", "P20", "P28", "P30", "P31", "P34"))

tumor <- subset(lung_paired, subset = Sample_Origin == "tLung")
normal <- subset(lung_paired, subset = Sample_Origin == "nLung")

patient_ids <- c("P06", "P08", "P09", "P18", "P19", "P20", "P28", "P30", "P31", "P34")
library(patchwork)


#distribution plot
lung_paired@meta.data$Sample_Origin <- lung_paired@meta.data$Condition
lung_paired@meta.data$Sample_Origin <- ifelse(
  lung_paired@meta.data$Sample_Origin == "nLung", "Normal",
  ifelse(lung_paired@meta.data$Sample_Origin == "tLung", "Tumor", lung_paired@meta.data$Sample_Origin)
)

lung_paired@meta.data$Sample_Origin <- factor(
  lung_paired@meta.data$Sample_Origin,
  levels = c("Tumor", "Normal"),
  ordered = TRUE
)


unique_patients <- sort(unique(lung_paired@meta.data$patient_id))
unique_sample_types <- levels(lung_paired@meta.data$Sample_Origin)

all_plot_list <- list()


for (patient in unique_patients) {

  for (sample in unique_sample_types) {
    patient_subset <- subset(lung_paired, subset = combined_group == paste0(patient, "_", sample))
 
    current_group <- paste(patient, sample, sep = "_")
    

    p <- DimPlot(
      patient_subset, 
      reduction = "tsne", 
      label = FALSE, 
      pt.size = 1, 
      group.by = "cell_type"
    ) +
    ggtitle(paste(patient, sample, sep = " - ")) +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 12, face = "bold")
    )
    
    all_plot_list[[current_group]] <- p
  }
}

combined_plot <- wrap_plots(all_plot_list, ncol = 2)
png("lung_paired_tsne_celltype_by_patient_sample.png", width = 6000, height = 12000, res = 300)
print(combined_plot)
dev.off()

#proportion
proportion_data <- lung_paired@meta.data %>%
  group_by(patient_id, Sample_Origin, cell_type) %>% 
  summarise(count = n()) %>% 
  mutate(proportion = count / sum(count)) %>%
  ungroup()

patients <- unique(proportion_data$patient_id)
all_plots <- list()
fixed_cell_types <- sort(unique(proportion_data$cell_type))

default_colors <- scales::hue_pal()(length(fixed_cell_types))
color_map <- setNames(default_colors, fixed_cell_types)

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
  
  all_plots[[patient]] <- patient_combined
}

final_combined <- wrap_plots(all_plots, ncol = 1) + 
  plot_annotation(title = "All Patients - Cell Type Proportions",
                  theme = theme(plot.title = element_text(size = 48, face = "bold", hjust = 0.5)))

png("all_patients_combined.png", width = 9000, height = 27000, res = 300)
print(final_combined)
dev.off()


#pseudobulk
pseudo_sce <- PseudobulkExpression(
  object = lung_paired, 
  group.by = c("Sample_Origin","patient_id", "cell_type"),
  assays = "RNA",                    
  method = "aggregate",
  layer = "counts"           
)
pseudo_sce_df <- as.data.frame(pseudo_sce)
write.csv(pseudo_sce_df, file = "pseudo_bulk.csv", row.names = TRUE)

#方法二
meta <- lung_paired@meta.data %>%
  as_tibble(rownames = "cell") %>%
  mutate(Patient = patient_id, Condition = Sample_Origin) %>%
  select(cell, Patient, Condition, cell_type)

counts_matrix <- as.matrix(GetAssayData(lung_paired, assay = "RNA", layer = "counts"))


meta$sample_id <- paste0(meta$Patient, "_", meta$Condition, "_", meta$cell_type)
mat_t <- t(counts_matrix)
pb_t <- rowsum(mat_t, group = meta$sample_id)
pb_counts <- t(pb_t)

coldata <- tibble(sample_id = colnames(pb_counts)) %>%
  separate(
    sample_id,
    into = c("Patient", "Condition", "cell_type"),
    sep = "_",
    remove = FALSE
  ) %>%
  mutate(
    Condition = factor(Condition, levels = c("Normal", "Tumor")),
    Patient = factor(Patient),
    cell_type = factor(cell_type)
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


#del each celltype
dds_original <- dds
cell_types <- unique(coldata$cell_type)
all_results <- list()
cell_types <- unique(meta$cell_type)
all_results <- list()

for (ct in cell_types) {
  message("Processing without cell type: ", ct)

  meta_filtered <- meta %>% filter(cell_type != ct)

  cells_to_keep <- meta_filtered$cell
  counts_matrix_filtered <- counts_matrix[, cells_to_keep]

  meta_filtered$sample_id <- paste0(meta_filtered$Patient, "_", meta_filtered$Condition, "_", meta_filtered$cell_type)
  mat_t_filtered <- t(counts_matrix_filtered)
  pb_t_filtered <- rowsum(mat_t_filtered, group = meta_filtered$sample_id)
  pb_counts_filtered <- t(pb_t_filtered)
 
  coldata_filtered <- tibble(sample_id = colnames(pb_counts_filtered)) %>%
    separate(
      sample_id,
      into = c("Patient", "Condition", "cell_type"),
      sep = "_",
      remove = FALSE
    ) %>%
    mutate(
      Condition = factor(Condition, levels = c("Normal", "Tumor")),
      Patient = factor(Patient),
      cell_type = factor(cell_type)
    ) %>%
    column_to_rownames("sample_id")
  

  dds_filtered <- DESeqDataSetFromMatrix(
    countData = pb_counts_filtered,
    colData   = coldata_filtered,
    design    = ~ Patient + Condition 
  )
  dds_filtered <- DESeq(dds_filtered)
  
  res_filtered <- results(dds_filtered, contrast = c("Condition", "Tumor", "Normal"))
  res_filtered <- lfcShrink(dds_filtered, coef = "Condition_Tumor_vs_Normal", type = "apeglm")

  res_tbl_filtered <- as_tibble(res_filtered, rownames = "gene") %>% arrange(padj)
  res_focus_filtered <- res_tbl_filtered %>%
    filter(gene %in% gene_list) %>%
    select(gene, log2FoldChange, pvalue, padj)
  

  res_focus_filtered$removed_cell_type <- ct
  

  all_results[[ct]] <- res_focus_filtered
}


combined_results <- bind_rows(all_results)
write.csv(combined_results, "differential_expression_by_cell_type.csv", row.names = FALSE)

#filter
counts_matrix <- GetAssayData(lung_paired, layer = "counts")
expr_cells <- which(counts_matrix["FZD10", ] > 0)
expr_ratio <- length(expr_cells) / ncol(lung_paired)

min_cell_percentage <- 0.001 

expressed_cells_per_gene <- Matrix::rowSums(counts_matrix > 0)
min_cells_required <- min_cell_percentage * ncol(lung_paired)
genes_to_keep <- expressed_cells_per_gene >= min_cells_required
lung_paired <- lung_paired[genes_to_keep, ]

