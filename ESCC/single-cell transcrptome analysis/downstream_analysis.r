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

#Monocle2
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("monocle")
library(monocle)

CD8 <- subset(t_cells, subset = seurat_clusters %in% c(0, 1, 4, 5))
CD4 <- subset(t_cells, subset = seurat_clusters %in% c(2,3))


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

markers$gene <- rownames(markers)
markers$neg_log10_pval <- -log10(markers$p_val)
markers$diffexpressed <- "No"
markers$diffexpressed[markers$avg_log2FC > 0.25 & markers$p_val < 0.05] <- "Up"
markers$diffexpressed[markers$avg_log2FC < -0.25 & markers$p_val < 0.05] <- "Down"