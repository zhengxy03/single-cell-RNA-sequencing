library(Seurat)

base_dir <- "/share/home/wangq/zxy/ESCC/matrix/GSE160269"

setwd("/share/home/wangq/zxy/ESCC/matrix/GSE160269")
dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)


# 定义样本信息
period1 <- c("IIIC", "IIA", "IIIA", "IA", "IIIA", "IIA", "IB", "IB", "IA", "IB", 
             "IIIB", "IIB", "IIIB", "IIIA", "IIB", "IIA", "IIIA", "IIA", "IIA", 
             "IB", "IA", "IB", "IIA", "IIIB", "IIIA", "IIIA", "IIIC", "IIA", 
             "IIIB", "IIIB", "IB", "IIA", "IB", "IIA", "IIA", "IIA", "IIA", 
             "IIIC", "IA", "IA", "IIIA", "IIA", "IIB", "IIIB", "IIIA", "IIA", 
             "IIB", "IIIB", "IIIA", "IA", "IIIB", "IIIA", "IIA", "IIIA", "IIIB", 
             "IIIA", "IIIC", "IB", "IA", "IIIB", "IIA", "IIIB", "IIA", "IB", 
             "IIIC", "IIIA", "IIIA", "IA", "IIA", "IIA", "IB", "IB", "IA", 
             "IIIA", "IIIB", "IIB", "IIIA", "IIIB", "IIB", "IIA", "IIIA", "IIA", 
             "IIA", "IB", "IIIC", "IA", "IIIA", "IIA", "IIIB", "IIIA", "IIIC", 
             "IIA", "IA", "IA", "IIA", "IIB", "IIIB", "IIIA", "IIA", "IIB", 
             "IIIB", "IIIA", "IA", "IIIB", "IIA", "IIIB", "IIIA", "IIIA", "IIIB", 
             "IIIA", "IIIC", "IB", "IA", "IIIB", "IIA", "IIIB", "IB", "IIA", 
             "IA", "IA")
period2 <- c("Ⅲ", "Ⅱ", "Ⅲ", "Ⅰ", "Ⅲ", "Ⅱ", "Ⅰ", "Ⅰ", "Ⅰ", "Ⅰ", "Ⅲ", "Ⅱ", 
             "Ⅲ", "Ⅲ", "Ⅱ", "Ⅱ", "Ⅲ", "Ⅱ", "Ⅱ", "Ⅰ", "Ⅰ", "Ⅰ", "Ⅱ", "Ⅲ", 
             "Ⅲ", "Ⅲ", "Ⅲ", "Ⅱ", "Ⅲ", "Ⅲ", "Ⅰ", "Ⅱ", "Ⅰ", "Ⅱ", "Ⅱ", "Ⅱ", 
             "Ⅱ", "Ⅲ", "Ⅰ", "Ⅰ", "Ⅲ", "Ⅱ", "Ⅱ", "Ⅲ", "Ⅲ", "Ⅱ", "Ⅱ", "Ⅲ", 
             "Ⅲ", "Ⅰ", "Ⅲ", "Ⅲ", "Ⅱ", "Ⅲ", "Ⅲ", "Ⅲ", "Ⅲ", "Ⅰ", "Ⅰ", "Ⅲ", 
             "Ⅱ", "Ⅲ", "Ⅱ", "Ⅰ", "Ⅲ", "Ⅲ", "Ⅲ", "Ⅰ", "Ⅱ", "Ⅱ", "Ⅰ", "Ⅰ", 
             "Ⅰ", "Ⅲ", "Ⅲ", "Ⅱ", "Ⅲ", "Ⅲ", "Ⅱ", "Ⅱ", "Ⅲ", "Ⅱ", "Ⅱ", "Ⅰ", 
             "Ⅲ", "Ⅰ", "Ⅲ", "Ⅱ", "Ⅲ", "Ⅲ", "Ⅲ", "Ⅱ", "Ⅰ", "Ⅰ", "Ⅱ", "Ⅱ", 
             "Ⅲ", "Ⅲ", "Ⅱ", "Ⅱ", "Ⅲ", "Ⅲ", "Ⅲ", "Ⅲ", "Ⅰ", "Ⅰ", "Ⅲ", "Ⅱ", 
             "Ⅲ", "Ⅰ", "Ⅱ", "Ⅰ", "Ⅰ")
sample_type <- c("T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", 
                 "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", 
                 "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", 
                 "T", "T", "T", "T", "N", "T", "N", "N", "N", "T", "T", "T", 
                 "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", 
                 "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", 
                 "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", 
                 "T", "T", "T", "T", "T", "T", "T", "N", "N", "T", "N", "T", 
                 "T", "N", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", 
                 "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T")
orig_ident <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 
                19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 10, 31, 22, 32, 
                31, 33, 34, 35, 36, 37, 38, 39, 59, 33, 60, 61, 62, 40, 41, 42, 
                43, 44, 45, 46, 34, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 35, 
                57, 1, 5, 3, 4, 6, 2, 7, 8, 9, 39, 11, 12, 14, 13, 15, 16, 17, 
                18, 19, 20, 36, 21, 26, 23, 24, 25, 27, 28, 58, 58, 59, 60, 29, 
                62, 40, 41, 61, 42, 43, 44, 45, 47, 30, 48, 46, 49, 50, 51, 52, 
                53, 54, 55, 56, 57, 32, 37, 38)
patient <- c("P63", "P61", "P57", "P56", "P62", "P54", "P52", "P49", "P48", "P17", 
             "P44", "P42", "P47", "P40", "P39", "P38", "P37", "P36", "P32", "P31", 
             "P30", "P16", "P27", "P26", "P24", "P28", "P23", "P22", "P21", "P20", 
             "P17", "P15", "P16", "P19", "P15", "P12", "P11", "P10", "P8", "P2", 
             "P1", "P9", "P130", "P12", "P128", "P127", "P126", "P130", "P128", 
             "P127", "P126", "P107", "P104", "P94", "P11", "P91", "P89", "P87", 
             "P84", "P83", "P82", "P79", "P75", "P74", "P80", "P10", "P65", "P63", 
             "P62", "P57", "P56", "P54", "P61", "P52", "P49", "P48", "P9", "P44", 
             "P42", "P40", "P47", "P39", "P38", "P37", "P36", "P32", "P31", "P8", 
             "P30", "P28", "P27", "P26", "P24", "P23", "P22", "P76", "P76", "P130", 
             "P128", "P21", "P126", "P130", "P128", "P127", "P127", "P126", "P107", 
             "P104", "P91", "P20", "P89", "P94", "P87", "P84", "P83", "P82", "P79", 
             "P75", "P74", "P80", "P65", "P19", "P2", "P1")

# 定义 sample_sources 列
sample_sources <- rep("GSE160269", length(period1))

# 创建映射表
srr_names <- paste0("SRR15093488", seq(0, 123))
mapping <- setNames(seq_along(srr_names), srr_names)

# 假设你的数据目录
dirs <- list.files(path = ".", pattern = "^SRR15093", full.names = TRUE)

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

  period1_vector <- rep(period1[sample_num], num_cells)
  period2_vector <- rep(period2[sample_num], num_cells)
  sample_type_vector <- rep(sample_type[sample_num], num_cells)
  orig_ident_vector <- rep(orig_ident[sample_num], num_cells)
  patient_vector <- rep(patient[sample_num], num_cells)
  sample_sources_vector <- rep(sample_sources[sample_num], num_cells)

  seurat_obj <- AddMetaData(seurat_obj, metadata = data.frame(
    period1 = period1_vector,
    period2 = period2_vector,
    sample_type = sample_type_vector,
    orig_ident = orig_ident_vector,
    patient = patient_vector,
    sample_sources = sample_sources_vector
  ))

  seurat_objs[[i]] <- seurat_obj
}

# 根据映射表生成正确的细胞 ID
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

# 保存 merged_seurat_obj 对象
saveRDS(merged_seurat_obj, file = "merged_seurat_obj.rds")