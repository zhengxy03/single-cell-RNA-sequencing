library(Seurat)

base_dir <- "/share/home/wangq/zxy/ESCC/treatment"
dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)

sample_types <- c("tumor", "tumor", "normal", "normal", 
                  "tumor", "normal", "normal", "tumor", 
                  "tumor", "tumor", "normal", "tumor", 
                  "normal", "tumor", "tumor", "normal", 
                  "tumor", "normal", "tumor", "tumor", 
                  "normal", "tumor", "normal", "tumor", 
                  "normal", "normal", "tumor", "normal",
                  "tumor", "normal", "tumor", "tumor",
                  "normal")
sample_sources <- c("GSE145370", "GSE145370", "GSE145370", "GSE145370",
                    "GSE160269", "GSE160269", "GSE160269", "GSE160269",
                    "GSE160269", "GSE160269", "GSE160269", "GSE160269",
                    "GSE160269", "GSE160269", "PRJNA777911", "PRJNA777911", "PRJNA777911",
                    "PRJNA777911", "PRJNA777911", "PRJNA777911", "PRJNA777911",
                    "PRJNA777911", "PRJNA777911", "PRJNA777911", "PRJNA777911",
                    "PRJNA777911", "PRJNA777911", "PRJNA777911", "PRJNA777911",
                    "PRJNA777911", "PRJNA777911", "PRJNA777911", "PRJNA777911")

mapping <- c(
  "SRR11094242" = 1,
  "SRR11094244" = 2,
  "SRR11094253" = 3,
  "SRR11094255" = 4,
  "SRR15093496" = 5,
  "SRR15093530" = 6,
  "SRR15093533" = 7,
  "SRR15093535" = 8,
  "SRR15093537" = 9,
  "SRR15093563" = 5,
  "SRR15093585" = 6,
  "SRR15093589" = 8,
  "SRR15093591" = 7,
  "SRR15093592" = 9,
  "SRR16796861" = 10,
  "SRR16796862" = 11,
  "SRR16796863" = 12,
  "SRR16796865" = 13,
  "SRR16796868" = 14,
  "SRR16796871" = 15,
  "SRR16796872" = 16,
  "SRR16796874" = 17,
  "SRR16796875" = 18,
  "SRR16796876" = 19,
  "SRR16796877" = 20,
  "SRR16796879" = 21,
  "SRR16796880" = 22,
  "SRR16796881" = 23,
  "SRR16796884" = 24,
  "SRR16796885" = 25,
  "SRR16796886" = 26,
  "SRR16796890" = 27,
  "SRR16796891" = 28
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
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 400 & percent.mt < 10)

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

#average genes per cell
genes_per_cell <- merged_seurat_obj$nFeature_RNA
average_genes <- mean(genes_per_cell)

merged_seurat_obj <- NormalizeData(merged_seurat_obj)
merged_seurat_obj <- FindVariableFeatures(merged_seurat_obj, nfeatures = 2000)

hvgs <- VariableFeatures(merged_seurat_obj)
merged_seurat_obj <- ScaleData(merged_seurat_obj, features = hvgs)
merged_seurat_obj <- RunPCA(merged_seurat_obj, features = hvgs, npcs = 20)

library(harmony)
merged_seurat_obj <- RunHarmony(merged_seurat_obj, "orig.ident")

# 保存 merged_seurat_obj 对象
saveRDS(merged_seurat_obj, file = "merged_seurat_obj.rds")