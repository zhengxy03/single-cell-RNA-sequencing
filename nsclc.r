#GSE117570
setwd("E:/project/nsclc/GSE117570_RAW")
library(Seurat)

#import
counts_matrix1 <- read.csv("GSM3304007_P1_Tumor_processed_data.txt.gz")
counts_matrix2 <- read.csv("GSM3304008_P1_Normal_processed_data.txt.gz")
counts_matrix3 <- read.csv("GSM3304009_P2_Tumor_processed_data.txt.gz")
counts_matrix4 <- read.csv("GSM3304010_P2_Normal_processed_data.txt.gz")
counts_matrix5 <- read.csv("GSM3304011_P3_Tumor_processed_data.txt.gz")
counts_matrix6 <- read.csv("GSM3304012_P3_Normal_processed_data.txt.gz")
counts_matrix7 <- read.csv("GSM3304013_P4_Tumor_processed_data.txt.gz")
counts_matrix8 <- read.csv("GSM3304014_P4_Normal_processed_data.txt.gz")

seurat_obj1 <- CreateSeuratObject(counts = counts_matrix1)

samples <- data.frame(
  gsm_prefix = c(007, 008, 009, 010, 011, 012, 013, 014),
  patient = rep(1:4, each = 2),
  type = rep(c("Tumor", "Normal"), 4)
)

seurat_list <- list()

for (i in 1:nrow(samples)) {
  filename <- sprintf("GSM3304%03d_P%d_%s_processed_data.txt.gz",
                     samples$gsm_prefix[i],
                     samples$patient[i],
                     samples$type[i])

  counts <- read.table(filename, header = TRUE, row.names = 1)
  
  seurat_list[[i]] <- CreateSeuratObject(
    counts = counts,
    project = paste0("P", samples$patient[i], "_", samples$type[i])
  )
  
  assign(paste0("seurat_P", samples$patient[i], "_", samples$type[i]), seurat_list[[i]])
}


names(seurat_list) <- paste0("P", samples$patient, "_", samples$type)

#add metadata
for (i in 1:length(seurat_list)) {
  current_patient <- paste0("P", samples$patient[i])
  current_type <- samples$type[i]

  seurat_list[[i]]$patient <- current_patient
  seurat_list[[i]]$sample_type <- current_type
}

#merge all objs
