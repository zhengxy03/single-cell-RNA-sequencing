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
saveRDS(matrix_data,file="matrix.rds")