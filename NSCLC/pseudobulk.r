library(SingleCellExperiment)
library(glmGamPoi)
sce <- as.SingleCellExperiment(lung_paired)
pseudobulk_sce <- pseudobulk(
    sce,
    group_by = vars(patient_id, Sample_Origin, cell_type),
    n_cells = n(),
)

fit <- glm_gp(
  pseudobulk_sce,
  design = ~ Sample_Origin + patient_id * cell_type,
  size_factor = "ratio")
head(fit$Beta)

fit <- glm_gp(
  pseudobulk_sce,
  design = ~ Sample_Origin * patient_id * cell_type,
  size_factor = "ratio")
head(fit$Beta)


beta_matrix <- fit$Beta
write.csv(x = beta_matrix,file = "model_beta.csv", quote = FALSE, row.names = TRUE)


design_matrix <- model.matrix(~ Sample_Origin + patient_id * cell_type, data = colData(pseudobulk_sce))


zero_cols <- which(colSums(design_matrix) == 0)


if (length(zero_cols) > 0) {
  design_matrix <- design_matrix[, -zero_cols]
  warning(paste("Removed", length(zero_cols), "columns with all zeros."))
}


fit <- glm_gp(
  pseudobulk_sce,
  design = design_matrix,
  size_factor = "ratio"
)



qr_result <- qr(design_matrix)
keep_cols <- qr_result$pivot[1:qr_result$rank]  # 保留独立列
design_matrix_clean <- design_matrix[, keep_cols]


cat("矩阵秩：", qr(design_matrix_clean)$rank, "\n")
cat("矩阵列数：", ncol(design_matrix_clean), "\n")

fit <- glm_gp(
  pseudobulk_sce,
  design = design_matrix_clean,
  size_factor = "ratio"
)
