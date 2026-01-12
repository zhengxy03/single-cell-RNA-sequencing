library(Seurat)
seurat_obj1 <- readRDS("GSE243013_TN.rds")
seurat_obj1$orig.ident <- "GSE243013"
seurat_obj2 <- readRDS("GSE241934_TN_LUAD.rds")
seurat_obj2$orig.ident <- "GSE241934"
merged_seurat_obj <- merge(seurat_obj1, y = c(seurat_obj2), add.cell.ids = c("GSE243013","GSE241934"))
merged_seurat_obj <- JoinLayers(merged_seurat_obj)
#31852 features across 791397 samples within 1 assay
counts_matrix <- GetAssayData(merged_seurat_obj, slot = "counts", assay = "RNA")
expressed_cells_per_gene <- Matrix::rowSums(counts_matrix > 0)

min_cell_percentage <- 0.001
min_cell_count <- 10
keep_genes <- expressed_cells_per_gene >= min_cell_percentage * ncol(merged_seurat_obj) & 
              expressed_cells_per_gene >= min_cell_count
merged_seurat_obj <- merged_seurat_obj[keep_genes, ]
merged_seurat_obj
#16227 features across 292715 samples within 1 assay

merged_seurat_obj <- subset(merged_seurat_obj, 
                     subset = nFeature_RNA > 200 & 
                              nFeature_RNA < 5000)
merged_seurat_obj
#17076 features across 789425 samples within 1 assay
genes_per_cell <- merged_seurat_obj$nFeature_RNA
average_genes <- mean(genes_per_cell)
average_genes
#1573.534



table(merged_seurat_obj@meta.data$LN_group,merged_seurat_obj@meta.data$cancer_type)
#        LUAD   LUSC
#  LN-  75136  98539
#  LN+ 252320 363430

library(dplyr)

library(harmony)
merged_seurat_obj <- NormalizeData(merged_seurat_obj)
merged_seurat_obj <- FindVariableFeatures(merged_seurat_obj, nfeatures = 2000)

hvgs <- VariableFeatures(merged_seurat_obj)
merged_seurat_obj <- ScaleData(merged_seurat_obj, features = hvgs)
merged_seurat_obj <- RunPCA(merged_seurat_obj, features = hvgs, npcs = 20)
merged_seurat_obj <- RunHarmony(merged_seurat_obj, "orig.ident")

png("elbowplot.png", width = 800, height = 600)
elbowplot <- ElbowPlot(merged_seurat_obj)
print(elbowplot)
dev.off()

library(ggplot2)
library(ggsci)
merged_seurat_obj <- RunUMAP(merged_seurat_obj, reduction = "harmony", dims = 1:20)
merged_seurat_obj <- FindNeighbors(merged_seurat_obj, reduction = "harmony", dims = 1:20)
merged_seurat_obj <- FindClusters(merged_seurat_obj, resolution = 0.3)


markers <- FindAllMarkers(merged_seurat_obj, 
                          only.pos = TRUE, 
                          min.pct = 0.25, 
                          logfc.threshold = 0.25, 
                          test.use = "wilcox")
#significant markers
library(dplyr)
significant_markers <- subset(markers, p_val_adj < 0.05)
significant_markers <- significant_markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.csv(significant_markers,"marker_top_20.csv")
saveRDS(merged_seurat_obj,file="TN_umap.rds")


