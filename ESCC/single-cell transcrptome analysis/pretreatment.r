library(Seurat)

base_dir <- "/share/home/wangq/zxy/ESCC/treatment"
base_dir <- "E:/project/ESCC/treatment"
setwd("E:/project/ESCC/treatment")
dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)

sample_types <- c("tumor", "tumor", "tumor", "normal", 
                  "normal", "tumor", "tumor", "tumor", 
                  "tumor", "tumor", "normal", "tumor", 
                  "normal", "tumor", "tumor", "normal", 
                  "tumor", "normal", "tumor", "tumor", 
                  "normal", "tumor", "normal", "tumor", 
                  "normal", "normal", "tumor", "normal",
                  "tumor", "normal", "tumor", "tumor",
                  "normal")
sample_sources <- c("GSE160269", "GSE160269", "GSE160269", "GSE160269",
                    "GSE160269", "GSE160269", "GSE160269", "GSE160269",
                    "GSE160269", "GSE160269", "GSE160269", "GSE160269",
                    "GSE160269", "GSE160269", "PRJNA777911", "PRJNA777911", "PRJNA777911",
                    "PRJNA777911", "PRJNA777911", "PRJNA777911", "PRJNA777911",
                    "PRJNA777911", "PRJNA777911", "PRJNA777911", "PRJNA777911",
                    "PRJNA777911", "PRJNA777911", "PRJNA777911", "PRJNA777911",
                    "PRJNA777911", "PRJNA777911", "PRJNA777911", "PRJNA777911")

mapping <- c(
  "SRR15093493" = 1,
  "SRR15093496" = 2,
  "SRR15093498" = 3,
  "SRR15093530" = 4,
  "SRR15093533" = 5,
  "SRR15093535" = 6,
  "SRR15093537" = 7,
  "SRR15093559" = 1,
  "SRR15093563" = 2,
  "SRR15093565" = 3,
  "SRR15093585" = 4,
  "SRR15093589" = 6,
  "SRR15093591" = 5,
  "SRR15093592" = 7,
  "SRR16796861" = 8,
  "SRR16796862" = 9,
  "SRR16796863" = 10,
  "SRR16796865" = 11,
  "SRR16796868" = 12,
  "SRR16796871" = 13,
  "SRR16796872" = 14,
  "SRR16796874" = 15,
  "SRR16796875" = 16,
  "SRR16796876" = 17,
  "SRR16796877" = 18,
  "SRR16796879" = 19,
  "SRR16796880" = 20,
  "SRR16796881" = 21,
  "SRR16796884" = 22,
  "SRR16796885" = 23,
  "SRR16796886" = 24,
  "SRR16796890" = 25,
  "SRR16796891" = 26
)

seurat_objs <- list()

for (i in seq_along(dirs)) {
  dir <- dirs[i]
  print(paste("current_dir:", dir))
  counts_matrix <- Read10X(data.dir = dir)
  
  # 从目录名中提取样本编号
  sample_name <- basename(dir)
  sample_num <- mapping[sample_name]
  
  seurat_obj <- CreateSeuratObject(counts = counts_matrix, project = paste0("Sample", sample_num))
  
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA >= 400 & nFeature_RNA <= 7500 & 
                                     percent.mt < 10 & 
                                     nCount_RNA >= 500 & nCount_RNA <= 50000)
  
  num_cells <- ncol(seurat_obj)

  sample_type_vector <- rep(sample_types[i], num_cells)
  sample_sources_vector <- rep(sample_sources[i], num_cells)

  seurat_obj <- AddMetaData(seurat_obj, metadata = data.frame(
    sample_type = sample_type_vector,
    sample_sources = sample_sources_vector
  ))

  seurat_objs[[i]] <- seurat_obj
}

# 根据映射表生成正确的细胞ID
cell_ids <- sapply(dirs, function(dir) {
  sample_name <- basename(dir)
  sample_num <- mapping[sample_name]
  paste0("Sample", sample_num)
})

merged_seurat_obj <- merge(seurat_objs[[1]], y = seurat_objs[-1], add.cell.ids = paste0("Sample", 1:length(dirs)))
merged_seurat_obj <- JoinLayers(merged_seurat_obj)

# 基因过滤
counts_matrix <- GetAssayData(merged_seurat_obj, slot = "counts", assay = "RNA")
expressed_cells_per_gene <- Matrix::rowSums(counts_matrix > 0)

min_cell_percentage <- 0.001
min_cell_count <- 10
keep_genes <- expressed_cells_per_gene >= min_cell_percentage * ncol(merged_seurat_obj) & 
              expressed_cells_per_gene >= min_cell_count
merged_seurat_obj <- merged_seurat_obj[keep_genes, ]



merged_seurat_obj <- NormalizeData(merged_seurat_obj)
merged_seurat_obj <- FindVariableFeatures(merged_seurat_obj, nfeatures = 2000)

hvgs <- VariableFeatures(merged_seurat_obj)
merged_seurat_obj <- ScaleData(merged_seurat_obj, features = hvgs)
merged_seurat_obj <- RunPCA(merged_seurat_obj, features = hvgs, npcs = 20)

# 双细胞检测和去除
# 确定最佳的 pK 值
library(DoubletFinder)
sweep.res.list <- paramSweep(merged_seurat_obj, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

# 估计双细胞比例
#nExp_poi <- round(0.075 * ncol(merged_seurat_obj)) 

# 进行双细胞检测
merged_seurat_obj <- doubletFinder(merged_seurat_obj, PCs = 1:20, pN = 0.25, pK = as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)])), nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# 去除双细胞
#doublet_col <- colnames(merged_seurat_obj@meta.data)[grep("DF.classification_0.25_0.29_2544", colnames(merged_seurat_obj@meta.data))]
#merged_seurat_obj <- subset(merged_seurat_obj, subset = DF.classifications_0.25_0.29_2544 == "Singlet")

pK_value <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))

# 分批处理
batch_size <- 10000
num_cells <- ncol(merged_seurat_obj)
num_batches <- ceiling(num_cells / batch_size)

batch_results <- list()

for (i in 1:num_batches) {
  start_idx <- (i - 1) * batch_size + 1
  end_idx <- min(i * batch_size, num_cells)
  
  # 提取当前批次的细胞
  batch_cells <- subset(merged_seurat_obj, cells = colnames(merged_seurat_obj)[start_idx:end_idx])
  
  # 检测双细胞
  nExp_poi <- round(0.075 * ncol(batch_cells))  # 假设双细胞比例为 7.5%
  batch_cells <- doubletFinder(batch_cells, PCs = 1:20, pN = 0.25, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  
  # 提取双细胞检测结果
  doublet_col <- grep("DF.classifications", colnames(batch_cells@meta.data), value = TRUE)
  batch_results[[i]] <- batch_cells@meta.data[[doublet_col]]
  
  # 释放内存
  rm(batch_cells)
  gc()
}

# 合并所有批次的结果
#all_results <- do.call(rbind, batch_results)
# 合并所有批次的结果
all_results <- unlist(batch_results)

# 将双细胞检测结果添加到原始的 Seurat 对象中
merged_seurat_obj@meta.data$DF.classifications <- all_results

# 去除双细胞
merged_seurat_obj <- subset(merged_seurat_obj, subset = DF.classifications == "Singlet")


#average genes per cell
genes_per_cell <- merged_seurat_obj$nFeature_RNA
average_genes <- mean(genes_per_cell)


library(harmony)
merged_seurat_obj <- RunHarmony(merged_seurat_obj, "orig.ident")


samples <- unique(merged_seurat_obj@meta.data$orig.ident)
sample_numbers <- as.numeric(gsub("[^0-9]", "", samples))
samples_ordered <- samples[order(sample_numbers)]
non_numeric_samples <- samples[is.na(sample_numbers)]  # 提取没有数字部分的样本名称
samples_ordered <- c(samples_ordered, non_numeric_samples)

merged_seurat_obj@meta.data$orig.ident <- factor(
  merged_seurat_obj@meta.data$orig.ident,
  levels = samples_ordered  # 按数字顺序设置水平
)


# 保存 merged_seurat_obj 对象
saveRDS(merged_seurat_obj, file = "merged_seurat_obj.rds")