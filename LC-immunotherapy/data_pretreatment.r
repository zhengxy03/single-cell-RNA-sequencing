#datasets pretreatment
#GSE303680
#

#GSE243013
#mv barcodes.csv.gz barcodes.tsv.gz
#mv features.csv.gz features.tsv.gz

library(Seurat)
library(Matrix)
library(readr)
dir <- "~/zxy/NSCLC/GSE243013"
matrix_data <- readMM(file.path(dir, "matrix.mtx.gz"))
barcodes <- read_tsv(file.path(dir, "barcodes.tsv.gz"), skip = 1, col_names = "barcode")$barcode
features <- read_tsv(file.path(dir, "features.tsv.gz"), skip = 1, col_names = "gene_symbol")$gene_symbol
print(paste("barcodes 数量:", length(barcodes)))
print(paste("features 数量:", length(features))) 
print(paste("matrix 维度:", dim(matrix_data)))

matrix_data_transposed <- t(matrix_data)

print("转置后维度:")
print(paste("转置Matrix dim:", paste(dim(matrix_data_transposed), collapse = " x ")))

rownames(matrix_data_transposed) <- features
colnames(matrix_data_transposed) <- barcodes

seurat_obj <- CreateSeuratObject(
  counts = matrix_data_transposed,
  project = "GSE243013"
)
library(readr)
library(Seurat)


metadata <- read_csv("metadata.csv.gz")

print("Metadata结构:")
print(dim(metadata))
print(colnames(metadata))
print(head(metadata))
common_cells <- intersect(rownames(metadata), colnames(seurat_obj))
print(paste("匹配的细胞数:", length(common_cells)))

rownames(metadata) <- metadata$cellID
metadata_to_add <- metadata[, -which(colnames(metadata) == "cellID")]
seurat_obj <- AddMetaData(seurat_obj, metadata = metadata_to_add)
unique_samples <- unique(seurat_obj$sampleID)

patient_mapping <- setNames(
  paste0("P", 1:length(unique_samples)),  # P1, P2, P3...
  unique_samples                           # 对应的sampleID
)
patients_vector <- patient_mapping[as.character(seurat_obj$sampleID)]
library(tibble)

# 创建包含细胞名和patients的tibble
patients_tibble <- tibble(
  cell_barcode = colnames(seurat_obj),
  patients = patients_vector
) %>%
  column_to_rownames("cell_barcode")

# 添加metadata
seurat_obj <- AddMetaData(
  seurat_obj,
  metadata = patients_tibble
)
#243
saveRDS(seurat_obj,file="GSE243013.rds") #only immune cells

filtered_samples <- read_csv("GSE243013_filtered.csv")
matching_samples <- intersect(seurat_obj$sampleID, filtered_samples$Tumor_Sample_Barcode)
print(paste("匹配的样本数:", length(matching_samples)))
#51
seurat_filtered <- subset(seurat_obj, subset = sampleID %in% matching_samples)
print(seurat_filtered)


library(dplyr)

stage_mapping <- setNames(filtered_samples$Stage, filtered_samples$Tumor_Sample_Barcode)
Stage <- stage_mapping[seurat_filtered@meta.data$sampleID]
seurat_filtered@meta.data$Stage <- Stage
saveRDS(seurat_filtered,file="GSE243013_TN.rds")
#31831 features across 256006 samples within 1 assay

#GSE241934
library(Seurat)
library(Matrix)
library(readr)
dir <- "~/zxy/NSCLC/GSE241934/RWC"
raw_data1 <- Read10X(data.dir = dir)
seurat_obj1 <- CreateSeuratObject(counts = raw_data1, project = "GSE241934_RWC")
metadata1 <- read_tsv("RWC/GSE241934_Real_Meta.txt.gz")
metadata_df <- as.data.frame(metadata1)
rownames(metadata_df) <- metadata_df$cellID
metadata_to_add <- metadata_df[, colnames(metadata_df) != "cellID", drop = FALSE]
common_cells <- intersect(rownames(metadata_to_add), colnames(seurat_obj1))
print(paste("匹配细胞数:", length(common_cells)))
seurat_obj1 <- AddMetaData(seurat_obj1, metadata = metadata_to_add)

dir <- "~/zxy/NSCLC/GSE241934/T"
raw_data2 <- Read10X(data.dir = dir)
seurat_obj2 <- CreateSeuratObject(counts = raw_data2, project = "GSE241934_T")
metadata2 <- read_tsv("T/GSE241934_IIT_Meta.txt.gz")
metadata_df <- as.data.frame(metadata2)
rownames(metadata_df) <- metadata_df$cellID
metadata_to_add <- metadata_df[, colnames(metadata_df) != "cellID", drop = FALSE]
common_cells <- intersect(rownames(metadata_to_add), colnames(seurat_obj2))
print(paste("匹配细胞数:", length(common_cells)))
seurat_obj2 <- AddMetaData(seurat_obj2, metadata = metadata_to_add)

merged_seurat_obj <- merge(seurat_obj1, y = c(seurat_obj2), add.cell.ids = c("RWC", "T"))
merged_seurat_obj <- JoinLayers(merged_seurat_obj)
LUAD <- subset(merged_seurat_obj, subset = Histology == "LUAD")
colnames(LUAD@meta.data)[colnames(LUAD@meta.data) == "Histology"] <- "cancer_type"
colnames(LUAD@meta.data)[colnames(LUAD@meta.data) == "Pathological.Response"] <- "pathological_response"

start_number <- 244  # 从P244开始
unique_samples <- unique(LUAD$sampleID)
#43
patient_mapping <- setNames(
  paste0("P", start_number:(start_number + length(unique_samples) - 1)),
  unique_samples
)
patients_vector <- patient_mapping[as.character(LUAD$sampleID)]
library(tibble)

# 创建包含细胞名和patients的tibble
patients_tibble <- tibble(
  cell_barcode = colnames(LUAD),
  patients = patients_vector
) %>%
  column_to_rownames("cell_barcode")

# 添加metadata
LUAD <- AddMetaData(
  LUAD,
  metadata = patients_tibble
)
#286

saveRDS(LUAD,file="GSE241934.rds")
#27693 features across 294914 samples within 1 assay

immune <- subset(LUAD, subset = major_cell_type %in% c("B","Mast","Myeloid","NK","T"))
#27693 features across 238087 samples within 1 assay
saveRDS(immune,file="GSE241934_immune.rds")

filtered_samples <- read_csv("GSE241934_filtered.csv")
matching_samples <- intersect(immune$sampleID, filtered_samples$SampleID)
print(paste("匹配的样本数:", length(matching_samples)))
#5
seurat_filtered <- subset(immune, subset = sampleID %in% matching_samples)
print(seurat_filtered)
#27693 features across 36709 samples within 1 assay
library(dplyr)

stage_mapping <- setNames(filtered_samples$Stage, filtered_samples$SampleID)
Stage <- stage_mapping[seurat_filtered@meta.data$sampleID]
seurat_filtered@meta.data$Stage <- Stage


saveRDS(seurat_filtered,file="GSE241934_TN.rds")





#GSE229353


#GSE207422
library(Seurat)
library(dplyr)
library(readr)

umi_matrix <- read.table("UMI_matrix.txt.gz", header = TRUE, row.names = 1, sep = "\t")
metadata <- read.csv("metadata_sorted.csv", header = TRUE, row.names = 1)

umi_prefixes <- sub("_\\d+$", "", colnames(umi_matrix))

# 检查匹配情况
print("UMI矩阵中的样本前缀:")
print(unique(umi_prefixes))
print("Metadata中的Sample:")
print(unique(metadata$Sample))

# 查找匹配的细胞
matching_cells <- colnames(umi_matrix)[umi_prefixes %in% metadata$Sample]

print(paste("找到", length(matching_cells), "个匹配的细胞"))

# 创建细胞到样本的映射
cell_to_sample <- setNames(umi_prefixes[umi_prefixes %in% metadata$Sample], 
                          matching_cells)

# 准备metadata
cell_metadata <- metadata[match(cell_to_sample, metadata$Sample), ]
rownames(cell_metadata) <- names(cell_to_sample)  # 使用完整的barcode作为行名

seurat_obj <- CreateSeuratObject(
  counts = umi_matrix[, matching_cells],
  meta.data = cell_metadata,
  project = "GSE207422"
)
filtered_samples <- read_csv("GSE207422_filtered.csv")
matching_samples <- intersect(seurat_obj$Sample, filtered_samples$Sample_ID)
print(paste("匹配的样本数:", length(matching_samples)))
#14
seurat_filtered <- subset(seurat_obj, subset = Sample %in% matching_samples)
print(seurat_filtered)
#24292 features across 87589 samples within 1 assay
saveRDS(seurat_filtered,file="GSE207422.rds")

seurat_filtered <- NormalizeData(seurat_filtered)
seurat_filtered <- FindVariableFeatures(seurat_filtered, nfeatures = 2000)

hvgs <- VariableFeatures(seurat_filtered)
seurat_filtered <- ScaleData(seurat_filtered, features = hvgs)
seurat_filtered <- RunPCA(seurat_filtered, features = hvgs, npcs = 20)
seurat_filtered <- RunUMAP(seurat_filtered, dims = 1:20)
seurat_filtered <- FindNeighbors(seurat_filtered, dims = 1:20)
seurat_filtered <- FindClusters(seurat_filtered, resolution = 0.5)
markers <- FindAllMarkers(seurat_filtered, 
                          only.pos = TRUE, 
                          min.pct = 0.1, 
                          logfc.threshold = 0.25, 
                          test.use = "wilcox")
#significant markers
library(dplyr)
significant_markers <- subset(markers, p_val_adj < 0.05)
significant_markers <- significant_markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.csv(significant_markers,"marker_top.csv")
identity_mapping <- c(
    "0" = "T cell",
    "1" = "DC",
    "2" = "Neutrophil",
    "3" = "B cell",
    "4" = "Macrophage",
    "5" = "Epithelial",
    "6" = "Plasma",
    "7" = "Proliferating cell",
    "8" = "Fibroblast",
    "9" = "Epithelial",
    "10" = "Epithelial",  
    "11" = "DC",  
    "12" = "Fibroblast",
    "13" = "Mast cell",
    "14" = "Ciliated cell",
    "15" = "Proliferating cell",    
    "16" = "Airway secretory cell",
    "17" = "Endothelial cell"
)
cell_type <- identity_mapping[seurat_filtered@meta.data$seurat_clusters]
seurat_filtered@meta.data$cell_type <- cell_type
cell_type_annotation <- DimPlot(seurat_filtered, reduction = "umap", label = FALSE, group.by = "cell_type")
library(ggplot2)
ggsave("celltype_annotation.png", plot = cell_type_annotation, width = 8, height = 6, dpi = 300)
immune <- subset(seurat_filtered, subset = cell_type %in% c("B cell","Mast cell","Macrophage","Neutrophil","T cell", "Plasma", "DC"))
#24292 features across 72172 samples within 1 assay
saveRDS(immune,file="GSE207422_immune.rds")

#GSE146100
library(Seurat)
library(dplyr)
library(readr)
library(data.table)

matrix <- fread("GSE146100.txt.gz", header = TRUE, sep = "\t")
dim(matrix)
#[1] 19799 11613
print(matrix[1:5, 1:5])

gene_names <- matrix[[1]]
counts_matrix <- as.matrix(matrix[, -1, with = FALSE])
rownames(counts_matrix) <- gene_names
seurat_obj <- CreateSeuratObject(
  counts = counts_matrix,
  project = "GSE146100"
)
#19799 features across 11612 samples within 1 assay

seurat_obj@meta.data$Response <- ifelse(seurat_obj@meta.data$orig.ident == "W2", "YES", "NO")
seurat_obj@meta.data$cancer_type <- "LUAD"
seurat_obj@meta.data$patients <- case_when(
  seurat_obj@meta.data$orig.ident == "W1" ~ "P301",
  seurat_obj@meta.data$orig.ident == "W2" ~ "P302", 
  seurat_obj@meta.data$orig.ident == "W3" ~ "P303"
)





#GSE179994
library(Seurat)
library(dplyr)
library(readr)

matrix <- readRDS("GSE179994.rds.gz")
seurat_obj <- CreateSeuratObject(
  counts = matrix, 
  project = "GSE179994"
)

metadata <- read_tsv("metadata.tsv.gz")
metadata_df <- as.data.frame(metadata)
rownames(metadata_df) <- metadata_df$cellid
metadata_to_add <- metadata_df[, colnames(metadata_df) != "cellid", drop = FALSE]

common_cells <- intersect(rownames(metadata_to_add), colnames(seurat_obj))
print(paste("匹配细胞数:", length(common_cells)))
seurat_obj <- AddMetaData(seurat_obj, metadata = metadata_to_add)
saveRDS(seurat_obj,file="GSE179994.rds")
#19790 features across 150849 samples within 1 assay
filtered_samples <- read_csv("GSE179994_filtered.csv")
matching_samples <- intersect(seurat_obj$patient, filtered_samples$Patient)
print(paste("匹配的样本数:", length(matching_samples)))
#8
seurat_filtered <- subset(seurat_obj, subset = patient %in% matching_samples)
print(seurat_filtered)
#19790 features across 55434 samples within 1 assay
library(dplyr)

cancer_type_mapping <- setNames(filtered_samples$cancer_type, filtered_samples$Patient)
cancer_type <- cancer_type_mapping[seurat_filtered@meta.data$patient]
seurat_filtered@meta.data$cancer_type <- cancer_type

Response_mapping <- setNames(filtered_samples$Response, filtered_samples$Patient)
Response <- Response_mapping[seurat_filtered@meta.data$patient]
seurat_filtered@meta.data$Response <- Response


patients_mapping <- setNames(filtered_samples$patients, filtered_samples$Patient)
patients <- patients_mapping[seurat_filtered@meta.data$patient]
seurat_filtered@meta.data$patients <- patients

saveRDS(seurat_filtered,file="GSE179994_pritumor.rds")


#GSE131933
library(data.table)
library(Seurat)
library(Matrix)

con <- gzfile("GSE131933_T1D7_gene_count.txt.gz", "r")
header_line <- readLines(con, n = 1)
close(con)

cell_names <- gsub("\"", "", header_line)  # 移除引号
cell_names <- unlist(strsplit(cell_names, " "))  # 按空格分割
cell_names <- cell_names[cell_names != ""]  # 移除空字符串

print(paste("细胞数量:", length(cell_names)))

count_data <- fread(
  "GSE131933_T1D7_gene_count.txt.gz", 
  header = FALSE,           # 没有标题行
  skip = 1,                 # 跳过第一行（细胞名）
  sep = " "                 # 空格分隔
)

print("数据维度:")
print(dim(count_data))
print("前5行前5列预览:")
print(count_data[1:5, 1:5])

# 提取基因名（第一列）
gene_names <- count_data$V1

# 创建表达矩阵（从第二列开始）
expression_matrix <- as.matrix(count_data[, -1])

# 设置行名（基因）和列名（细胞）
rownames(expression_matrix) <- gene_names
colnames(expression_matrix) <- cell_names  # 使用之前提取的cell_names

print("表达矩阵创建成功!")
print(paste("基因数:", nrow(expression_matrix)))
print(paste("细胞数:", ncol(expression_matrix)))
print("矩阵预览:")
print(expression_matrix[1:5, 1:5])

expression_sparse <- as(expression_matrix, "dgCMatrix")

# 创建Seurat对象
seurat_obj <- CreateSeuratObject(
  counts = expression_sparse,
  project = "GSE131933_T1D7"
)

