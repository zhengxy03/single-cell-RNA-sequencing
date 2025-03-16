emt_genes <- c("VIM", "CDH1", "ZEB1", "SNAI1", "TWIST1", "FN1", "CDH2",
                "SNAI2", "ZEB2", "MMP2", "MMP9", "BCL9", "CTNNB1")

emt_genes <- emt_genes[emt_genes %in% rownames(merged_seurat_obj)]

#calculate emt score
merged_seurat_obj <- AddModuleScore(
  object = merged_seurat_obj,
  features = list(emt_genes),
  name = "EMT_Score"
)
head(seurat_obj@meta.data)



npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(4)

merged_seurat_obj@meta.data$EMT_Score1_scaled <- scale(merged_seurat_obj@meta.data$EMT_Score1)
summary(merged_seurat_obj@meta.data$EMT_Score1_scaled)

p <- FeaturePlot(
  object = merged_seurat_obj,          # Seurat 对象
  features = "EMT_Score1_scaled",            # EMT 评分列名
  label = TRUE,                       # 是否显示聚类标签
  repel = TRUE,                       # 避免标签重叠
  order = TRUE,                       # 按评分高低排序
  cols = c("#4DBBD5FF", "#E64B35FF"),
  pt.size = 0.5                       # 点的大小
) + 
  ggtitle("EMT Score UMAP Distribution")  # 添加标题

ggsave("emt-plot.png", plot = p)

npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(15)
p <- VlnPlot(
  object = merged_seurat_obj,
  features = "EMT_Score1",
  group.by = "cell_type",
  pt.size = 0,
  cols = npg_extended
)
ggsave("emt-plot2.png", plot = p)

#emt feature group
emt_cells <- colnames(merged_seurat_obj)[merged_seurat_obj$EMT_Score1 > 0.5]  # 调整阈值
emt_seurat_obj <- subset(merged_seurat_obj, cells = emt_cells)
saveRDS(emt_seurat_obj, file = "emt_seurat_obj.rds")

library(NMF)
#emt_expression_matrix <- as.matrix(GetAssayData(emt_seurat_obj, layer = "counts"))
#emt_expression_matrix <- emt_expression_matrix[rowSums(emt_expression_matrix) > 10, ]
#emt_expression_matrix <- emt_expression_matrix[rowSums(emt_expression_matrix > 0) >= 0.1 * ncol(emt_expression_matrix), ]
#nrow(emt_expression_matrix)

emt_seurat_obj <- FindVariableFeatures(emt_seurat_obj, nfeatures = 5000)
emt_expression_matrix <- GetAssayData(emt_seurat_obj, layer = "counts")[VariableFeatures(emt_seurat_obj), ]
summary(as.vector(emt_expression_matrix))
emt_expression_matrix <- as.matrix(emt_expression_matrix)
emt_expression_matrix <- log1p(emt_expression_matrix)


set.seed(123)
sampled_cells <- sample(colnames(emt_expression_matrix), 5000)
emt_expression_matrix <- emt_expression_matrix[, sampled_cells]

emt_expression_matrix <- emt_expression_matrix[rowSums(emt_expression_matrix) > 0, ]
emt_expression_matrix <- emt_expression_matrix[complete.cases(emt_expression_matrix), ]
emt_expression_matrix[is.na(emt_expression_matrix)] <- 0

set.seed(123)
library(parallel)
nmf_result <- nmf(
  emt_expression_matrix, 
  rank = 3,  # 降低 rank
  method = "lee",  # 使用更快的算法
  nrun = 5,  # 减少 nrun
  .pbackend = "mc",  # 启用并行计算
  .options = "v"     # 显示进度
)
#nmf_result <- nmf(emt_expression_matrix, rank = 3, method = "brunet", nrun = 5)
# 测试多个 k 值
#k_values <- 3:10
#nmf_eval <- nmf(emt_expression_matrix, rank = k_values, method = "brunet", nrun = 10)
#plot(nmf_eval)

summary(nmf_result)

#提取表达程序
gene_weights <- basis(nmf_result)
summary(as.vector(gene_weights))
print(head(rownames(gene_weights)))
cell_weights <- coef(nmf_result)
head(gene_weights)
head(cell_weights)
num_modules <- ncol(gene_weights)  

library(GSVA)
library(clusterProfiler)
library(org.Hs.eg.db)

hallmark_pathways <- read.gmt("h.all.v2024.1.Hs.symbols.gmt")
hallmark_list <- split(hallmark_pathways$gene, hallmark_pathways$term)
hallmark_pathways$term <- sub("^HALLMARK_", "", hallmark_pathways$term)
unique(hallmark_pathways$term)
head(hallmark_list)

#功能富集
top_genes_per_module <- apply(gene_weights, 2, function(x) {
  rownames(gene_weights)[order(x, decreasing = TRUE)][1:50]
})
print(top_genes_per_module)

enrichment_results <- lapply(top_genes_per_module, function(genes) {
  # 将基因符号转换为 Entrez ID
  entrez_ids <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID
  
  # 使用 GO 或 KEGG 进行富集分析
  enrich_result <- enrichGO(
    gene = entrez_ids,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",  # 生物过程
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05
  )
  
  return(enrich_result)
})

# 查看每个模块的富集结果
for (i in 1:length(enrichment_results)) {
  cat("Module", i, "Enrichment Results:\n")
  print(head(enrichment_results[[i]]))
  cat("\n")
}


module_similarity <- cor(gene_weights)
hc <- hclust(as.dist(1 - module_similarity), method = "complete")
plot(hc, main = "Module Hierarchical Clustering", xlab = "", sub = "")

emt_pathway_genes <- hallmark_pathways$gene[hallmark_pathways$term == "EPITHELIAL_MESENCHYMAL_TRANSITION"]
emt_core_genes <- lapply(top_genes_per_module, function(genes) {
  intersect(genes, emt_pathway_genes)
})


结果：
write.csv(gene_weights, file = "module_gene_weights.csv")
write.csv(cell_weights, file = "module_cell_weights.csv")

for (i in 1:length(enrichment_results)) {
  write.csv(enrichment_results[[i]], file = paste0("module_", i, "_enrichment_results.csv"))
}

#plot
library(pheatmap)
zscore_gene_weights <- t(scale(t(gene_weights)))

# 检查标准化后的值分布
summary(as.vector(zscore_gene_weights))

# 绘制热图
pheatmap(
  zscore_gene_weights,
  clustering_method = "complete",
  color = colorRampPalette(c("blue", "white", "red"))(100),
  main = "Z-score Normalized Gene Weights per Module"
)


ordered_genes <- order(apply(gene_weights, 1, which.max))  # 按最大权重模块排序
ordered_modules <- order(colnames(gene_weights))          # 按模块名称排序

# 重新排列 gene_weights
ordered_gene_weights <- gene_weights[ordered_genes, ordered_modules]





gene_names <- rownames(gene_weights)
top_genes_names <- lapply(top_genes_per_module, function(indices) {
  gene_names[indices]
})

#层次聚类
module_similarity <- cor(t(gene_weights))
hc <- hclust(as.dist(1 - module_similarity), method = "complete")


enrichment_results <- list()
for (module in 1:num_modules) {
  # 提取当前模块的前50个基因
  module_genes <- top_genes_per_module[, module]

  # 创建GSVA参数对象
  gsva_param <- gsvaParam(
    exprData = emt_expression_matrix,
    geneSets = hallmark_pathways,
    kcdf = "Gaussian"
  )

  # 使用GSVA进行富集分析
  gsva_result <- gsva(gsva_param)

  # 提取通路得分
  pathway_scores <- data.frame(colnames(gsva_result), t(gsva_result))
  colnames(pathway_scores) <- c("Pathway", "Score")

  # 按得分排序并保留前6个富集通路
  sorted_pathways <- pathway_scores[order(-pathway_scores$Score), ]
  top_pathways <- sorted_pathways[1:6, ]

  # 存储结果
  enrichment_results[[paste0("Module_", module)]] <- top_pathways
}

head(enrichment_results)

#EMT core genes
# 初始化一个列表来存储每个模块的EMT核心基因
emt_core_genes <- list()

for (module in 1:num_modules) {
  # 获取当前模块的富集结果
  module_result <- enrichment_results[[paste0("Module_", module)]]
  
  # 查找与EMT相关的通路
  emt_related <- grep("EPITHELIAL_MESENCHYMAL_TRANSITION", module_result$Pathway, ignore.case = TRUE)
  
  if (length(emt_related) > 0) {
    # 提取与EMT相关的通路
    emt_pathway <- module_result$Pathway[emt_related]
    
    # 提取当前模块的前50个基因
    module_genes <- top_genes_per_module[, module]
    
    # 将与EMT相关的基因存储为当前模块的EMT核心基因
    emt_core_genes[[paste0("Module_", module)]] <- module_genes
  }
}

# 查看EMT核心基因
print(emt_core_genes)


#plot
library(ggplot2)

# 绘制每个模块的前6个富集通路
for (module in 1:num_modules) {
  # 获取当前模块的富集结果
  module_result <- enrichment_results[[paste0("Module_", module)]]
  
  # 绘制条形图
  ggplot(module_result, aes(x = reorder(Pathway, Score), y = Score)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    xlab("Pathway") +
    ylab("Score") +
    ggtitle(paste0("Module ", module, " Enriched Pathways"))
}


















top_genes <- apply(gene_weights, 2, function(x) {
  rownames(gene_weights)[order(x, decreasing = TRUE)][1:50]
})

enrich_results <- lapply(1:ncol(top_genes), function(i) {
  genes <- top_genes[, i]
  enrichGO(genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH")
})

pathway_descriptions <- lapply(enrich_results, function(res) {
  if (is.null(res) || length(res) == 0) {
    return(character(0))
  }
  if (is.null(res$Description) || length(res$Description) == 0) {
    return(character(0))
  }
  res$Description
})

pathway_df <- data.frame(
  pathway = unlist(pathway_descriptions),
  stringsAsFactors = FALSE
)
write.csv(pathway_df, "enrich_results.csv", row.names = FALSE)

top_pathways <- lapply(enrich_results, function(res) {
  head(res$Description, 6)
})
emt_pathways <- grep("EPITHELIAL_MESENCHYMAL_TRANSITION", top_pathways, value = TRUE)
emt_core_genes <- unique(unlist(lapply(emt_pathways, function(pathway) {
  genes <- enrich_results[[which(top_pathways == pathway)]]$geneID
  strsplit(genes, "/")[[1]]
})))

