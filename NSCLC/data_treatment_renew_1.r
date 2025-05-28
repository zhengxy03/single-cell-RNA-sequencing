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


# 过滤低表达基因（min.cells < 0.1%总细胞数）
min_cells_threshold <- ncol(seurat_obj) * 0.001  # 0.1%的细胞数

seurat_obj <- subset(seurat_obj, features = rownames(seurat_obj)[Matrix::rowSums(raw_counts > 0) >= min_cells_threshold])




seurat_obj <- NormalizeData(
  seurat_obj,
  normalization.method = "LogNormalize",
  scale.factor = 1e6,
  assay = "RNA"
)


seurat_obj <- FindVariableFeatures(
  seurat_obj,
  selection.method = "vst",
  nfeatures = 2000,
  mean.cutoff = c(0.0125, 3),
  dispersion.cutoff = c(0.5, Inf)
)
hvgs <- VariableFeatures(seurat_obj)


seurat_obj <- ScaleData(
  seurat_obj,
  features = hvgs,
  do.scale = FALSE,
  do.center = TRUE,
  scale.max = 10
)


seurat_obj <- RunPCA(
  seurat_obj,
  features = hvgs,
  npcs = 50,
  verbose = TRUE
)


n_significant_pcs <- 15


seurat_obj <- RunTSNE(
  seurat_obj,
  dims = 1:n_significant_pcs,
  perplexity = 30,
  verbose = TRUE
)

npg_pal <- pal_npg()(10)
cell_types <- sort(unique(lung_paired@meta.data$Cell_type))
color_map <- setNames(rep(npg_pal, length.out = length(cell_types)), cell_types)


p <- DimPlot(
  lung_paired,
  reduction = "tsne",
  group.by = "Cell_type",
  label = FALSE,
  pt.size = 0.8,
  order = cell_types
) +
scale_color_manual(
  values = color_map,
  name = "Cell Type"
) +
theme(
  plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
  legend.text = element_text(size = 12),
  legend.title = element_text(size = 14, face = "bold"),
  legend.position = "right"
)


ggsave(
  filename = "lung_paired_tsne_CellType_npg.png",
  plot = p,
  width = 10,）
  height = 8,
  dpi = 300,
  bg = "white"
)

saveRDS(seurat_obj, file = "GSE131907_all_processed.rds")


seurat_obj <- FindNeighbors(
  seurat_obj,
  dims = 1:n_significant_pcs,
  k.param = 20
)

seurat_obj <- FindClusters(
  seurat_obj,
  resolution = 0.5,
  algorithm = 1  # Louvain算法
)

lung <- subset(seurat_obj, subset = Sample_Origin %in% c("nLung", "tLung"))
lung@meta.data$patient_id <- gsub("LUNG_([TN][0-9]+)", "P\\1", lung@meta.data$Sample)
lung@meta.data$patient_id <- gsub("T", "", lung@meta.data$patient_id)
lung@meta.data$patient_id <- gsub("N", "", lung@meta.data$patient_id)
lung_paired <- subset(lung, subset = patient_id %in% c("P06", "P08", "P09", "P18", "P19", "P20", "P28", "P30", "P31", "P34"))


#distribution plot

lung_paired@meta.data$Sample_Origin <- ifelse(
  lung_paired@meta.data$Sample_Origin == "nLung", "Normal",
  ifelse(lung_paired@meta.data$Sample_Origin == "tLung", "Tumor", lung_paired@meta.data$Sample_Origin)
)

lung_paired@meta.data$Sample_Origin <- factor(
  lung_paired@meta.data$Sample_Origin,
  levels = c("Tumor", "Normal"),
  ordered = TRUE
)
lung_paired@meta.data$combined_group <- paste0(
  lung_paired@meta.data$patient_id, "_", 
  lung_paired@meta.data$Sample_Origin
)
unique_patients <- sort(unique(lung_paired@meta.data$patient_id))
unique_sample_types <- levels(lung_paired@meta.data$Sample_Origin)

all_plot_list <- list()


cell_types <- sort(unique(lung_paired@meta.data$Cell_type))
npg_pal <- pal_npg()(10)

color_map <- setNames(rep(npg_pal, length.out = length(cell_types)), cell_types)

for (patient in unique_patients) {
    for (sample in unique_sample_types) {
        patient_subset <- subset(lung_paired, subset = combined_group == paste0(patient, "_", sample))
        
        current_group <- paste(patient, sample, sep = "_")

        p <- DimPlot(
            patient_subset, 
            reduction = "tsne", 
            label = FALSE, 
            pt.size = 1, 
            group.by = "Cell_type"
        ) +
            scale_color_manual(values = color_map) +
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

png("lung_paired_tsne_CellType_by_patient_sample2.png", width = 3000, height = 10000, res = 300)
print(combined_plot)
dev.off()


#proportion
proportion_data <- lung_paired@meta.data %>%
  group_by(patient_id, Sample_Origin, Cell_type) %>% 
  summarise(count = n()) %>% 
  mutate(proportion = count / sum(count)) %>%
  ungroup()

patients <- unique(proportion_data$patient_id)
all_plots <- list()
fixed_cell_types <- sort(unique(proportion_data$Cell_type))


npg_pal <- pal_npg()(10)

color_map <- setNames(rep(npg_pal, length.out = length(fixed_cell_types)), fixed_cell_types)

for (patient in patients) {
  patient_data <- proportion_data %>% filter(patient_id == patient)
  sample_types <- unique(patient_data$Sample_Origin)
  pie_plots <- list()
  
  for (sample in sample_types) {
    sample_data <- patient_data %>% filter(Sample_Origin == sample)
    p <- ggplot(sample_data, aes(x = "", y = proportion, fill = Cell_type)) +
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
  
  bar_chart <- ggplot(patient_data, aes(x = Sample_Origin, y = proportion, fill = Cell_type)) +
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

png("all_patients_combined_CellType.png", width = 9000, height = 27000, res = 300)
print(final_combined)
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
cell_types <- unique(lung_paired@meta.data$Cell_type)
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
write.csv(combined_results, "differential_expression_by_cell_type_5_26.csv", row.names = FALSE)