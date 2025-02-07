library(Seurat)
library(ggplot2)
library(dplyr)

base_dir <- "E:/project/ESCC/GSE145370/matrix/"
dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)

sample_types <- ifelse(1:length(dirs) %% 2 == 1, "Tumor", "Normal")
sample_sources <- c("S133", "S133", "S134", "S134", "S135", "S135", "S149", "S149", "S150", "S150", "S158", "S158", "S159", "S159")

seurat_objs <- list()
for (i in seq_along(dirs)) {
  dir <- dirs[i]
  print(paste("current_dir:", dir))
  counts_matrix <- Read10X(data.dir = dir)
  seurat_obj <- CreateSeuratObject(counts = counts_matrix, project = paste0("Sample", i))
  
  num_cells <- ncol(seurat_obj)

  sample_type_vector <- rep(sample_types[i], num_cells)
  sample_sources_vector <- rep(sample_sources[i], num_cells)

  seurat_obj <- AddMetaData(seurat_obj, metadata = data.frame(
    sample_type = sample_type_vector,
    sample_sources = sample_sources_vector
  ))

  seurat_objs[[i]] <- seurat_obj
}
print(seurat_objs)

merged_seurat_obj <- merge(seurat_objs[[1]], y = seurat_objs[-1], add.cell.ids = paste0("Sample", 1:length(dirs)))
merged_seurat_obj <- JoinLayers(merged_seurat_obj)

