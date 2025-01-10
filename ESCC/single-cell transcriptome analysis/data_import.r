library(Seurat)
#format 1:barcodes.tsv.gz(cell barcode),features.tsv.gz(gene ID),matrix.mtx.gz
#GSE203115
setwd("E:/project/ESCC/GSE203115/")
counts_matrix1 <- Read10X(".") #integrate three files
class(counts_matrix1) #get expressin matrix, one row for a gene and one column for a cell (in sparse matrix dgCMatrix format)
seurat_obj1 <- CreateSeuratObject(counts = counts_matrix1) #create Seurat object

#format 2:counts_matrix.csv.gz
setwd()
counts_matrix2 <- read.csv()
seurat_obj2 <- CreateSeuratObject(counts = counts_matrix2)
