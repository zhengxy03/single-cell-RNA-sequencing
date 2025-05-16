#GSE117570
setwd("E:/project/nsclc/GSE117570_RAW")
library(Seurat)

#import raw data
counts_matrix1 <- read.csv("GSM3304007_P1_Tumor_processed_data.txt.gz")
counts_matrix2 <- read.csv("GSM3304008_P1_Normal_processed_data.txt.gz")
counts_matrix3 <- read.csv("GSM3304009_P2_Tumor_processed_data.txt.gz")
counts_matrix4 <- read.csv("GSM3304010_P2_Normal_processed_data.txt.gz")
counts_matrix5 <- read.csv("GSM3304011_P3_Tumor_processed_data.txt.gz")
counts_matrix6 <- read.csv("GSM3304012_P3_Normal_processed_data.txt.gz")
counts_matrix7 <- read.csv("GSM3304013_P4_Tumor_processed_data.txt.gz")
counts_matrix8 <- read.csv("GSM3304014_P4_Normal_processed_data.txt.gz")

seurat_obj1 <- CreateSeuratObject(counts = counts_matrix1)

samples <- data.frame(
  gsm_prefix = c(007, 008, 009, 010, 011, 012, 013, 014),
  patient = rep(1:4, each = 2),
  type = rep(c("Tumor", "Normal"), 4)
)

seurat_list <- list()

for (i in 1:nrow(samples)) {
  filename <- sprintf("GSM3304%03d_P%d_%s_processed_data.txt.gz",
                     samples$gsm_prefix[i],
                     samples$patient[i],
                     samples$type[i])

  counts <- read.table(filename, header = TRUE, row.names = 1)
  
  seurat_list[[i]] <- CreateSeuratObject(
    counts = counts,
    project = paste0("P", samples$patient[i], "_", samples$type[i])
  )
  
  assign(paste0("seurat_P", samples$patient[i], "_", samples$type[i]), seurat_list[[i]])
}


names(seurat_list) <- paste0("P", samples$patient, "_", samples$type)

#add metadata
for (i in 1:length(seurat_list)) {
  current_patient <- paste0("P", samples$patient[i])
  current_type <- samples$type[i]

  seurat_list[[i]]$patient <- current_patient
  seurat_list[[i]]$sample_type <- current_type
}

#merge all objs
merged_seurat_obj <- merge(x = seurat_list[[1]], 
                       y = seurat_list[2:length(seurat_list)],
                       add.cell.ids = names(seurat_list),
                       project = "CombinedAnalysis")
merged_seurat_obj <- JoinLayers(merged_seurat_obj)

merged_seurat_obj <- NormalizeData(merged_seurat_obj)
merged_seurat_obj <- FindVariableFeatures(merged_seurat_obj, nfeatures = 2000)

hvgs <- VariableFeatures(merged_seurat_obj)
merged_seurat_obj <- ScaleData(merged_seurat_obj, features = hvgs)
merged_seurat_obj <- RunPCA(merged_seurat_obj, features = hvgs, npcs = 20)



library(harmony)
merged_seurat_obj <- RunHarmony(merged_seurat_obj, "orig.ident")

ElbowPlot(merged_seurat_obj)
merged_seurat_obj <- FindNeighbors(merged_seurat_obj, reduction = "harmony", dims = 1:15)

merged_seurat_obj <- FindClusters(merged_seurat_obj, resolution = 0.5)
merged_seurat_obj <- RunUMAP(merged_seurat_obj, reduction = "harmony", dims = 1:20)

DimPlot(merged_seurat_obj, reduction = "umap", label = TRUE)

target_genes <- c("CDKN2A","FZD10","NOTCH1","PDGFRA","WNT7B")
FeaturePlot(merged_seurat_obj, features = target_genes)

identity_mapping <- c(
    "0" = "Effector T cell",
    "1" = "Alveolar type II cell",
    "2" = "Neutrophil",
    "3" = "Macrophage",
    "4" = "Dendritic cell",
    "5" = "Naive T cell",
    "6" = "NK cell",
    "7" = "Basal Epithelial cell",
    "8" = "Ciliated Epithelial cell",
    "9" = "B cell",
    "10" = "Fibroblast",
    "11" = "Airway secretory cell",
    "12" = "Ciliated Airway Epithelial cell",
    "13" = "Proliferating cell"
)
cell_type <- identity_mapping[seurat_obj@meta.data$seurat_clusters]
seurat_obj@meta.data$cell_type <- cell_type
DimPlot(seurat_obj, reduction = "umap", label = FALSE, group.by = "cell_type")







#GSE131907
raw_counts <- readRDS("GSE131907_Lung_Cancer_raw_UMI_matrix.rds")
lung_matrix <- raw_counts[, grepl("_LUNG_", colnames(raw_counts))]

merged_data <- readRDS("GSE131907_Lung_Cancer_normalized_log2TPM_matrix.rds")
str(merged_data)

expression_matrix <- as.matrix(merged_data)
rownames(expression_matrix) <- rownames(merged_data)
cell_ids <- colnames(merged_data)

lung_samples <- grep("_LUNG_", colnames(expression_matrix), value = TRUE)

lung_expr <- expression_matrix[, lung_samples]

sample_info <- data.frame(
  cell_id = lung_samples,
  tissue = "LUNG",

  sample_type = gsub(".*_(LUNG|LN|NS)_(.*)", "\\2", lung_samples),

  condition = gsub("([A-Z]+)([0-9]+)", "\\1", gsub(".*_(LUNG|LN|NS)_([A-Z]+[0-9]+)", "\\2", lung_samples)),

  sample_id = gsub("([A-Z]+)([0-9]+)", "\\2", gsub(".*_(LUNG|LN|NS)_([A-Z]+[0-9]+)", "\\2", lung_samples)),
  
  stringsAsFactors = FALSE
)
head(sample_info)

#sample_info <- read.table(
  file = "GSE131907_Lung_Cancer_cell_annotation.txt",
  sep = "\t",
  header = TRUE,
  fill = TRUE
)



#lung_samples <- sample_info[grepl("_LUNG_", sample_info$Index), ]

seurat_obj <- CreateSeuratObject(
  counts = lung_matrix,       # 使用原始计数矩阵（raw_counts）
  meta.data = sample_info,    # 添加元数据
  project = "Lung_Cancer"     # 项目名称
)
head(seurat_obj@meta.data)

seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 2000)

hvgs <- VariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = hvgs)
seurat_obj <- RunPCA(seurat_obj, features = hvgs, npcs = 20)

seurat_obj <- FindNeighbors(seurat_obj, dims = 1:15)

seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)

DimPlot(seurat_obj, reduction = "umap", label = TRUE)

target_genes <- c("CDKN2A","FZD10","NOTCH1","PDGFRA","WNT7B")
FeaturePlot(seurat_obj, features = target_genes)

markers <- FindAllMarkers(seurat_obj, 
                          only.pos = TRUE, 
                          min.pct = 0.25, 
                          logfc.threshold = 0.25, 
                          test.use = "wilcox")
library(dplyr)
significant_markers <- subset(markers, p_val_adj < 0.05)
significant_markers <- significant_markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.csv(significant_markers,"marker_top.csv")

identity_mapping <- c(
    "0" = "T cells",
    "1" = "T cells",
    "2" = "Myeloid cells",
    "3" = "T cells",
    "4" = "Myeloid cells",
    "5" = "B cells",
    "6" = "Epithelial cells",
    "7" = "Myeloid cells",
    "8" = "Fibroblasts",
    "9" = "Myeloid cells",
    "10" = "Mast cells",
    "11" = "Epithelial cells",
    "12" = "Epithelial cells",
    "13" = "Myeloid cells",
    "14" = "Endothelial cells",
    "15" = "Plasma",
    "16" = "Epithelial cells",
    "17" = "Proliferating cells",
    "18" = "Proliferating cells",
    "19" = "Undetermined",
    "20" = "Epithelial cells"
)

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
DimPlot(seurat_obj, reduction = "umap", label = FALSE, group.by = "cell_type")


levels(seurat_obj@meta.data$condition)


seurat_obj@meta.data$condition <- factor(
  seurat_obj@meta.data$condition,
  levels = c("Tumor", "Normal")  # 强制 T 在前
)

# 重新绘制 VlnPlot
VlnPlot(
  seurat_obj, 
  features = target_genes, 
  ncol = 2, 
  pt.size = 0, 
  group.by = "cell_type", 
  split.by = "condition"
) + 
  ylab(expression(log[2]("CPM/100 + 1"))) +
  theme(axis.title.y = element_text(size = 12)) + theme(legend.position = "right") 


target_genes <- c("CDKN2A")
plot_data <- FetchData(
  seurat_obj,
  vars = c(target_genes, "cell_type", "condition"),
  layer = "data"
)


plots <- lapply(target_genes, function(gene) {

  y_max <- max(plot_data[[gene]], na.rm = TRUE) * 1.05
  

  p <- ggplot(plot_data, aes(x = cell_type, y = .data[[gene]], fill = condition)) +
    geom_violin(scale = "width", trim = TRUE) +
    stat_compare_means(
      aes(group = condition),
      method = "wilcox.test",
      label = "p.format",
      label.y = y_max * 0.9
    ) +
    labs(y = expression(log[2]("CPM/100 + 1"))) +
    labs(title = gene) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.x = element_blank()
    )
  
  return(p)
})


wrap_plots(plots, ncol = 1)






plot_data <- FetchData(
    seurat_obj,
    vars = c(target_genes, "cell_type", "condition"),
    layer = "data"
)


cell_types <- unique(plot_data$cell_type)

conditions <- unique(plot_data$condition)




plots <- lapply(target_genes, function(gene) {
    p_values <- data.frame(
        cell_type = character(),
        p_value = numeric(),
        is_expressed = logical(),
        stringsAsFactors = FALSE
    )
    
    for(ct in cell_types) {
        ct_data <- plot_data[plot_data$cell_type == ct, ]
        
        # 提取两个条件下的基因表达值
        group1 <- ct_data[ct_data$condition == conditions[1], gene]
        group2 <- ct_data[ct_data$condition == conditions[2], gene]
        
        # 检测基因在该细胞类型中是否有表达（至少一组有非零值）
        is_expressed <- any(group1 > 0) || any(group2 > 0)
        
        # 如果有表达，则执行Wilcoxon检验
        if(is_expressed) {
            test_result <- wilcox.test(group1, group2)
            p_val <- round(test_result$p.value, 3)
        } else {
            p_val <- NA
        }
        
    
        p_values <- rbind(p_values, data.frame(
            cell_type = ct,
            p_value = p_val,
            is_expressed = is_expressed
        ))
    }

    y_max <- max(plot_data[[gene]], na.rm = TRUE) * 1.05
    
    p <- ggplot(plot_data, aes(x = cell_type, y = .data[[gene]], fill = condition)) +
        geom_violin(scale = "width", trim = TRUE) +
        labs(y = expression(log[2]("CPM/100 + 1"))) +
        labs(title = gene) +
        theme_classic() +
        theme(
            legend.position = "none",
            axis.text.x = element_blank(), 
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank()  
        )
    
    p <- p + geom_text(
        data = p_values[p_values$is_expressed, ],
        aes(x = cell_type, y = y_max * 0.9, label = paste("p =", p_value)),
        inherit.aes = FALSE,
        size = 3
    )
    
    return(p)
})


plots[[length(plots)]] <- plots[[length(plots)]] + 
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust =1, size = 14),
        axis.ticks.x = element_line()
    )

final_plot <- wrap_plots(plots, ncol = 1, guides = "collect") & 
    theme(legend.position = "right",
          legend.justification = "top",
          legend.box.margin = margin(0, 0, 0, 20)) & 
    labs(fill = "sample type") 

final_plot