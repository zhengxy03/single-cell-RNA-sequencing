library(Seurat)
library(ggplot2)
library(dplyr)
#format 1:barcodes.tsv.gz(cell barcode),features.tsv.gz(gene ID),matrix.mtx.gz
#GSE203115
setwd("E:/project/ESCC/GSE203115/matrix")
dir1 = "E:/project/ESCC/GSE203115/matrix/Esc1"
dir2 = "E:/project/ESCC/GSE203115/matrix/Esc2"
dir3 = "E:/project/ESCC/GSE203115/matrix/Esc3"

counts_matrix1 <- Read10X(data.dir = dir1) #renamed and integrate three files
counts_matrix2 <- Read10X(data.dir = dir2)
counts_matrix3 <- Read10X(data.dir = dir3)

class(counts_matrix1) 
#get expressin matrix, one row for a gene and one column for a cell (in sparse matrix dgCMatrix format)

seurat_obj1 <- CreateSeuratObject(counts = counts_matrix1, project = "Sample1") #create Seurat object
seurat_obj2 <- CreateSeuratObject(counts = counts_matrix2, project = "Sample2")
seurat_obj3 <- CreateSeuratObject(counts = counts_matrix3, project = "Sample3")
merged_seurat_obj <- merge(seurat_obj1, y = c(seurat_obj2,seurat_obj3), add.cell.ids = c("Sample1", "Sample2", "Sample3"))
merged_seurat_obj <- JoinLayers(merged_seurat_obj)

merged_seurat_obj$response_status <- ifelse(
  merged_seurat_obj$orig.ident %in% c("Sample3"), "Non-Responder", "Responder"
)


#GSE196756
setwd("E:/project/ESCC/GSE196756/matrix/GSE196756_RAW")
dirs <- c(
  "E:/project/ESCC/GSE196756/matrix/GSE196756_RAW/GSM5900215_1T",
  "E:/project/ESCC/GSE196756/matrix/GSE196756_RAW/GSM5900216_1N",
  "E:/project/ESCC/GSE196756/matrix/GSE196756_RAW/GSM5900217_2T",
  "E:/project/ESCC/GSE196756/matrix/GSE196756_RAW/GSM5900218_2N",
  "E:/project/ESCC/GSE196756/matrix/GSE196756_RAW/GSM5900219_3T",
  "E:/project/ESCC/GSE196756/matrix/GSE196756_RAW/GSM5900220_3N"
)

seurat_objs <- list()

for (i in seq_along(dirs)) {
  dir <- dirs[i]
  print(paste("current_dir:", dir))
  counts_matrix <- Read10X(data.dir = dir)
  seurat_obj <- CreateSeuratObject(counts = counts_matrix, project = paste0("Sample", i))
  seurat_objs[[i]] <- seurat_obj
}
print(seurat_objs)


#format 2:counts_matrix.csv.gz
setwd()
counts_matrix2 <- read.csv()
seurat_obj2 <- CreateSeuratObject(counts = counts_matrix2)
