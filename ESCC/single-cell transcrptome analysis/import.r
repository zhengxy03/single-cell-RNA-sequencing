library(Seurat)
setwd("F:/matrix/GSE160269")

metadata <- read.csv("数据集导入.csv")

samples <- metadata$sample
seurat_list <- list()

for (sample in samples) {
  # 读取样本数据（假设每个样本的数据存储在当前目录下的文件夹中）
  data_dir <- sample  # 样本文件夹名与样本名相同
  data <- Read10X(data.dir = data_dir)
  
  # 创建 Seurat 对象
  seurat_obj <- CreateSeuratObject(counts = data, project = sample)
  
  # 计算线粒体基因的百分比
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  
  # 进行质量过滤
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA >= 400 & nFeature_RNA <= 7500 & 
                         percent.mt < 10 & 
                         nCount_RNA >= 500 & nCount_RNA <= 50000)
  
  # 将过滤后的 Seurat 对象添加到列表
  seurat_list[[sample]] <- seurat_obj
}

merged_seurat_obj <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = samples)

sample_metadata <- metadata[match(samples, metadata$sample), ]
sample_names <- merged_seurat_obj@meta.data$orig.ident
cell_metadata <- sample_metadata[match(sample_names, sample_metadata$sample), ]
merged_seurat_obj@meta.data <- cbind(merged_seurat_obj@meta.data, cell_metadata)
merged_seurat_obj <- JoinLayers(merged_seurat_obj)
head(merged_seurat_obj@meta.data)

merged_seurat_obj <- NormalizeData(merged_seurat_obj)
merged_seurat_obj <- FindVariableFeatures(merged_seurat_obj, nfeatures = 2000)

hvgs <- VariableFeatures(merged_seurat_obj)
merged_seurat_obj <- ScaleData(merged_seurat_obj, features = hvgs)
merged_seurat_obj <- RunPCA(merged_seurat_obj, features = hvgs, npcs = 20)

library(DoubletFinder)
sweep.res.list <- paramSweep(merged_seurat_obj, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)


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

all_results <- unlist(batch_results)

# 将双细胞检测结果添加到原始的 Seurat 对象中
merged_seurat_obj@meta.data$DF.classifications <- all_results

# 去除双细胞
merged_seurat_obj <- subset(merged_seurat_obj, subset = DF.classifications == "Singlet")



#remove period 0
merged_seurat_obj <- subset(merged_seurat_obj, subset = period1 != 0)

library(harmony)
merged_seurat_obj <- RunHarmony(merged_seurat_obj, "sample_sources")
