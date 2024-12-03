setwd("//wsl.localhost/Ubuntu/home/zxy0303/project/scRNA")
library(Seurat)
counts_matrix1 <- read.table("scc13.txt", header = TRUE, sep = "\t")
seurat_obj1 <- CreateSeuratObject(counts = counts_matrix1, project = "scc13")
counts_matrix2 <- read.table("cal27.txt", header = TRUE, sep = "\t")
seurat_obj2 <- CreateSeuratObject(counts = counts_matrix2, project = "cal27")
seurat_obj1[["percent.mito"]] <- PercentageFeatureSet(seurat_obj1, pattern = "^hg19-MT-")
seurat_obj2[["percent.mito"]] <- PercentageFeatureSet(seurat_obj2, pattern = "^hg19-MT-")
seurat_obj1 <- subset(seurat_obj1, subset = nFeature_RNA > 200 & percent.mito < 10)
seurat_obj2 <- subset(seurat_obj2, subset = nFeature_RNA > 200 & percent.mito < 10)
merged_seurat_obj <- merge(seurat_obj1, y = c(seurat_obj2), add.cell.ids = c("scc13", "cal27"))

merged_seurat_obj[["RNA"]]$counts
LayerData(merged_seurat_obj, assay = "RNA", layer = "counts")
merged_seurat_obj <- JoinLayers(merged_seurat_obj)
dim(merged_seurat_obj[["RNA"]]$counts)

library(stringr)
phe = merged_seurat_obj@meta.data
table(phe$orig.ident)
View(phe)

merged_seurat_obj <- NormalizeData(merged_seurat_obj)
merged_seurat_obj <- FindVariableFeatures(merged_seurat_obj, selection.method = "vst", nfeatures = 2000, mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.5, Inf))
merged_seurat_obj <- ScaleData(merged_seurat_obj, vars.to.regress = c("nCount_RNA", "percent.mito"))
merged_seurat_obj <- RunPCA(merged_seurat_obj, features = VariableFeatures(object = merged_seurat_obj))
DimPlot(merged_seurat_obj, reduction = "pca")
ElbowPlot(merged_seurat_obj)
merged_seurat_obj <- RunUMAP(merged_seurat_obj, dims = 1:15)
DimPlot(merged_seurat_obj, reduction = "umap", label = TRUE)
merged_seurat_obj <- FindNeighbors(merged_seurat_obj, dims = 1:15)
merged_seurat_obj <- FindClusters(merged_seurat_obj, resolution = 0.15)
cluster_markers <- FindAllMarkers(merged_seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

cluster_markers %>% subset(p_val<0.05)
list_marker <- cluster_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
df_marker=data.frame(p_val = list_marker$p_val,
                     avg_log2FC = list_marker$avg_log2FC,
                     pct.1 = list_marker$pct.1,
                     pct.2 = list_marker$pct.2,
                     p_val_adj = list_marker$p_val_adj,
                     cluster = list_marker$cluster,
                     gene = list_marker$gene)
write.csv(df_marker,"marker.csv")