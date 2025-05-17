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



seurat_obj <- CreateSeuratObject(
  counts = lung_matrix,
  meta.data = sample_info,
  project = "Lung_Cancer"
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
            p_val <- formatC(test_result$p.value, format = "e", digit = 2)
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



#boxplot
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
    
    # 检测基因在该细胞类型中是否有表达
    is_expressed <- any(group1 > 0) || any(group2 > 0)
    
    # 如果有表达，则执行Wilcoxon检验
    if(is_expressed) {
      test_result <- wilcox.test(group1, group2)
      p_val <- formatC(test_result$p.value, format = "e", digit = 2)
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
    geom_boxplot(

      width = 0.6,                 
      alpha = 0.7,                 
      position = position_dodge(0.7)
    ) +
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



#

target_gene <- "CDKN2A"
target_cell_type <- "B cell"


full_data <- FetchData(
    seurat_obj,
    vars = c(target_gene, "cell_type", "condition"),
    layer = "data"
) %>%
    filter(cell_type == target_cell_type)

plot_data <- full_data

conditions <- unique(full_data$condition)


if(nrow(plot_data) == 0) {
    stop(paste("No cells with non-zero", target_gene, "expression in", target_cell_type))
}


group1_full <- full_data[full_data$condition == conditions[1], target_gene]
group2_full <- full_data[full_data$condition == conditions[2], target_gene]


is_expressed_full <- any(group1_full > 0) || any(group2_full > 0)


if(is_expressed_full && all(table(full_data$condition) > 0)) {
    test_result <- wilcox.test(group1_full, group2_full)
    p_val <- formatC(test_result$p.value, format = "e", digit = 2)
    p_text <- paste("p =", p_val)
} else {
    p_text <- "Insufficient data for test"
}


y_max <- max(plot_data[[target_gene]], na.rm = TRUE) * 1.05


p <- ggplot(plot_data, aes(x = condition, y = .data[[target_gene]], fill = condition)) +
    geom_violin(scale = "width", trim = TRUE) +  # 显示完整小提琴图
    geom_boxplot(width = 0.2, fill = "white", outlier.size = 2) +  # 显示箱线图离群值
    labs(
        title = paste(target_gene, "expression in", target_cell_type),
        subtitle = "Non-zero expression only (p-value from full data)",  # 明确p值来源
        y = expression(log[2]("CPM/100 + 1")),
        x = "Condition"
    ) +
    theme_classic() +
    theme(
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray50"),
        axis.text.x = element_text(angle = 0, size = 12),
        axis.title = element_text(size = 14)
    ) +
    ylim(0, y_max) +  # 使用非零表达数据的y轴范围
    guides(fill = guide_legend(title = "Condition"))


p <- p + geom_text(
    data = data.frame(),
    aes(x = 1.5, y = y_max * 0.95, label = p_text),
    inherit.aes = FALSE,
    size = 5
)

p