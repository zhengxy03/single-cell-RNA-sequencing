#GSE131907

library(Seurat)
library(dplyr)
library(ggplot2)

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
p1 <- DimPlot(seurat_obj, reduction = "tsne", label = TRUE)
ggsave("clusters.png", plot = p1, dpi = 300)
p2 <- DimPlot(seurat_obj, reduction = "tsne", label = FALSE, group.by = "Cell_subtype")
ggsave("annotation.png", plot = p2, dpi = 300)

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
generate_combined_plot <- function(data, condition_name) {
  plots <- list()
  
  for (patient in patient_ids) {
    pat <- subset(data, subset = patient_id == patient)

    p <- DimPlot(
      pat, 
      reduction = "tsne", 
      label = FALSE, 
      group.by = "cell_type",
    ) +
    ggtitle(patient) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 10),
      legend.position = "bottom"  
    )
    
    plots[[patient]] <- p
  }

  combined_plot <- wrap_plots(
    plots, 
    ncol = 2
  ) +
  plot_annotation(
    title = condition_name,
    theme = theme(
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      plot.margin = margin(10, 10, 10, 10)
    )
  )
  
  return(combined_plot)
}


tumor_combined <- generate_combined_plot(tumor, "Tumor")
normal_combined <- generate_combined_plot(normal, "Normal")


ggsave(
  "tumor_patients_combined.png",
  plot = tumor_combined,
  width = 24, height = 48, dpi = 300
)

ggsave(
  "normal_patients_combined.png",
  plot = normal_combined,
  width = 24, height = 48, dpi = 300
)


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