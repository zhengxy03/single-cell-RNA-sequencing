fibroblasts <- subset(merged_seurat_obj, subset = cell_type == "Fibroblast")
fibroblasts <- NormalizeData(fibroblasts)
fibroblasts <- FindVariableFeatures(fibroblasts, nfeatures = 2000)
hvgs <- VariableFeatures(fibroblasts)
fibroblasts <- ScaleData(fibroblasts, features = hvgs)
fibroblasts <- RunPCA(fibroblasts, features = hvgs, npcs = 20)

fibroblasts <- RunHarmony(fibroblasts, "sample_sources")
fibroblasts <- RunUMAP(fibroblasts, dims = 1:15, reduction = "harmony")
fibroblasts <- FindNeighbors(fibroblasts, dims = 1:15, reduction = "harmony")
fibroblasts <- FindClusters(fibroblasts, resolution = 0.3)

DimPlot(fibroblasts, reduction = "umap", label = TRUE) +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  ggtitle(NULL) +
  guides(color = guide_legend(title = NULL, override.aes = list(size = 5))) +
  theme(
    text = element_text(size = 8, face = "bold"),
    axis.text.x = element_text(size = 12), 
    axis.text.y = element_text(size = 12), 
    axis.title.x = element_text(size = 12), 
    axis.title.y = element_text(size = 12), 
    plot.title = element_text(size = 12),
    legend.text = element_text(size = 8), # 图例文本字体大小
    legend.title = element_text(size = 8) # 图例标题字体大小（如果有标题的话，这里原代码设置为NULL）
  )

fibroblast_markers <- FindAllMarkers(fibroblasts, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "wilcox")
fibroblast_significant_markers <- subset(fibroblast_markers, p_val_adj < 0.05)
#write.csv(fibroblast_significant_markers, "fibroblast_all_marker.csv")
fibroblast_significant_markers <- fibroblast_significant_markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
write.csv(fibroblast_significant_markers, "fibroblast_top_marker.csv")

neg_fibroblast_significant_markers <- fibroblast_significant_markers %>% 
  group_by(cluster) %>% 
  top_n(n = -50, wt = avg_log2FC)
write.csv(neg_fibroblast_significant_markers, "neg_fibroblast_top_marker.csv")

immune_cluster_ids <- c(9)
fibroblasts <- subset(fibroblasts, idents = setdiff(levels(Idents(fibroblasts)), immune_cluster_ids))

FeaturePlot(merged_seurat_obj, features = "CXCL5")

mapping <- c(
    "0" = "iCAFs",
    "1" = "CD248+vCAFs",
    "2" = ""

)