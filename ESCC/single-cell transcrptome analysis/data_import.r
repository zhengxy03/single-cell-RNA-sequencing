library(Seurat)
library(ggplot2)
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






#format 2:counts_matrix.csv.gz
setwd()
counts_matrix2 <- read.csv()
seurat_obj2 <- CreateSeuratObject(counts = counts_matrix2)
