# single-cell-RNA-sequencing
# 1 download data
GSE144240-GSE144236
## 1.1 
```
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE144nnn/GSE144236/suppl/GSE144236%5Fpatient%5Fmetadata%5Fnew.txt.gz
gzip -d 
if (!requireNamespace("Seurat", quietly = TRUE))
    install.packages("Seurat")
library(Seurat)

setwd("//wsl.localhost/Ubuntu/home/zxy0303/project/scRNA")
data <- read.table("patient_metadata_new.txt", header = TRUE, sep = "\t")

seurat_obj <- CreateSeuratObject(data)

seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & percent.mt < 10)

