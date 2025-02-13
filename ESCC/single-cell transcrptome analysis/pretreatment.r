library(Seurat)

base_dir <- "/share/home/wangq/zxy/ESCC/treatment"
dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)

sample_types <- c("tumor", "normal", "normal", "tumor", 
                  "tumor", "normal", "tumor", "normal", 
                  "tumor", "tumor", "normal", "tumor", 
                  "normal", "tumor", "tumor", "normal", 
                  "tumor", "normal", "tumor", "normal", 
                  "normal", "tumor", "normal", "tumor", 
                  "normal", "tumor", "tumor", "normal")
sample_sources <- c("GSE160269", "GSE160269", "GSE160269", "GSE160269",
                    "GSE160269", "GSE160269", "GSE160269", "GSE160269",
                    "GSE160269", "PRJNA777911", "PRJNA777911", "PRJNA777911",
                    "PRJNA777911", "PRJNA777911", "PRJNA777911", "PRJNA777911",
                    "PRJNA777911", "PRJNA777911", "PRJNA777911", "PRJNA777911",
                    "PRJNA777911", "PRJNA777911", "PRJNA777911", "PRJNA777911",
                    "PRJNA777911", "PRJNA777911", "PRJNA777911", "PRJNA777911")

seurat_objs <- list()

for (i in seq_along(dirs)) {
  dir <- dirs[i]
  print(paste("current_dir:", dir))
  counts_matrix <- Read10X(data.dir = dir)
  seurat_obj <- CreateSeuratObject(counts = counts_matrix, project = paste0("Sample", i))
  
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