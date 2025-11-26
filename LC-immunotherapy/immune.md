# immune cells difference between response-different samples
## merge
```
library(Seurat)
seurat_obj1 <- readRDS("GSE243013.rds")
seurat_obj2 <- readRDS("allcellsets_immune.rds")
seurat_obj3 <- readRDS("GSE179994_pritumor.rds")

filter_seurat_with_your_criteria <- function(obj) {
  # 1. 先过滤低质量细胞（可选，但建议加，进一步减容）
  # （若你不想过滤细胞，可删除这步，直接保留所有细胞）
  obj <- subset(
    obj,
    subset = nFeature_RNA > 200 & nFeature_RNA < 6000
  )
  
  # 2. 按你的标准过滤基因（核心步骤）
  counts_matrix <- GetAssayData(obj, assay = "RNA", layer = "counts")  # 提取单个对象的counts
  expressed_cells_per_gene <- Matrix::rowSums(counts_matrix > 0)       # 计算该对象中每个基因的表达细胞数
  min_cell_percentage <- 0.001
  min_cell_count <- 10
  # 你的过滤条件（每个对象单独计算，因为ncol(obj)是单个对象的细胞数）
  keep_genes <- expressed_cells_per_gene >= min_cell_percentage * ncol(obj) & 
                expressed_cells_per_gene >= min_cell_count
  # 过滤该对象的基因
  obj_filt <- obj[keep_genes, ]
  
  # 3. 可选：用DietSeurat精简对象（删除冗余数据，进一步减容）
  obj_filt <- DietSeurat(
    obj_filt,
    assay = "RNA",
    counts = TRUE, data = TRUE  # 保留counts和data，后续合并用
  )
  
  return(obj_filt)
}

# 对3个原始对象分别执行过滤（用你的标准）
seurat_obj1_filt <- filter_seurat_with_your_criteria(seurat_obj1)
seurat_obj2_filt <- filter_seurat_with_your_criteria(seurat_obj2)
seurat_obj3_filt <- filter_seurat_with_your_criteria(seurat_obj3)



merged_seurat_obj <- merge(seurat_obj1, y = c(seurat_obj2,seurat_obj3))
merged_seurat_obj <- JoinLayers(merged_seurat_obj)
merged_seurat_obj

counts_matrix <- GetAssayData(merged_seurat_obj, slot = "counts", assay = "RNA")
expressed_cells_per_gene <- Matrix::rowSums(counts_matrix > 0)

min_cell_percentage <- 0.001
min_cell_count <- 10
keep_genes <- expressed_cells_per_gene >= min_cell_percentage * ncol(merged_seurat_obj) & 
              expressed_cells_per_gene >= min_cell_count
merged_seurat_obj <- merged_seurat_obj[keep_genes, ]
merged_seurat_obj

genes_per_cell <- merged_seurat_obj$nFeature_RNA
average_genes <- mean(genes_per_cell)
average_genes

library(harmony)
merged_seurat_obj <- NormalizeData(merged_seurat_obj)
merged_seurat_obj <- FindVariableFeatures(merged_seurat_obj, nfeatures = 2000)

hvgs <- VariableFeatures(merged_seurat_obj)
merged_seurat_obj <- ScaleData(merged_seurat_obj, features = hvgs)
merged_seurat_obj <- RunPCA(merged_seurat_obj, features = hvgs, npcs = 20)
merged_seurat_obj <- RunHarmony(merged_seurat_obj, "orig.ident")
saveRDS(merged_seurat_obj,file=immune_pre.rds")