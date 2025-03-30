#expr_matrix <- GetAssayData(epi_T, layer = "counts")
# 获取表达矩阵
T_expr_matrix <- GetAssayData(T, layer = "counts")
epi_expr_matrix <- GetAssayData(epi, layer = "counts")

# 找到共同的基因
common_genes <- intersect(rownames(T_expr_matrix), rownames(epi_expr_matrix))

# 过滤表达矩阵，只保留共同基因
T_expr_matrix <- T_expr_matrix[common_genes, ]
epi_expr_matrix <- epi_expr_matrix[common_genes, ]

# 合并表达矩阵
expression_matrix <- cbind(T_expr_matrix, epi_expr_matrix)

epi_clusters <- epi$seurat_clusters
# 生成注释文件
annotations_file <- data.frame(
    cells = colnames(expression_matrix),  # 细胞名称
    cell_type = ifelse(colnames(expression_matrix) %in% colnames(T_expr_matrix), "T cells", 
                      paste0("Epi_Cluster_", epi_clusters[match(colnames(expression_matrix), colnames(epi_expr_matrix))]))  # 细胞类型
)

# 加载必要的库
library(infercnv)
library(AnnoProbe)
library(biomaRt)

# 获取基因名称
gene_names <- rownames(expression_matrix)

# 基因注释
gene_annotation <- annoGene(
    IDs = gene_names,  # 基因名称
    ID_type = "SYMBOL",  # 基因名称类型（如 SYMBOL 或 ENSEMBL）
    species = "human"  # 物种
)

# 生成基因位置参考文件
gene_location_ref <- data.frame(
    gene = gene_annotation$SYMBOL,  # 基因符号
    chr = gene_annotation$chr,  # 染色体名称
    start = gene_annotation$start,  # 基因起始位置
    stop = gene_annotation$end
)

# 保存基因位置参考文件
write.table(gene_location_ref, file = "gene_location_ref.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# 读取基因位置参考文件
gene_location_ref <- read.table("gene_location_ref.txt", header = FALSE, sep = "\t")

# 检查重复的基因名称
duplicated_genes <- gene_location_ref$V1[duplicated(gene_location_ref$V1)]
print(duplicated_genes)
# 去除重复的基因名称
gene_location_ref <- gene_location_ref[!duplicated(gene_location_ref$V1), ]

# 保存为文件
write.table(gene_location_ref, file = "gene_location_ref.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
groupFiles <- 'groupFiles.txt'
write.table(annotations_file, file = groupFiles, sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE)

# 3. 保存基因位置参考文件到文件
geneFile <- 'geneFile.txt'
write.table(gene_location_ref, file = geneFile, sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE)

# 4. 创建 infercnv 对象
infercnv_obj <- CreateInfercnvObject(
    raw_counts_matrix = expression_matrix,  # 表达矩阵文件路径
    annotations_file = groupFiles,  # 注释文件路径
    delim = "\t",  # 文件分隔符
    gene_order_file = geneFile,  # 基因位置参考文件路径
    ref_group_names = c("T cells")  # 以 T 细胞作为参考
)

infercnv_obj <- infercnv::run(
    infercnv_obj,  # infercnv 对象
    cutoff = 0.2,  # 设置 cutoff 值（默认值为 0.1）
    out_dir = "./infercnv",  # 输出目录
    cluster_by_groups = TRUE,  # 按组聚类
    denoise = TRUE,  # 去噪
    HMM = FALSE  # 使用 HMM 模型
)



#抽样
set.seed(123)

# 设置抽样数量（例如抽取 10,000 个细胞）
n_samples <- 10000
sampled_epi_cells <- sample(colnames(epi_expr_matrix), size = n_samples, replace = FALSE)
sampled_expression_matrix <- epi_expr_matrix[, sampled_epi_cells]
sampled_annotations_file <- data.frame(
    cells = sampled_epi_cells,
    cell_type = paste0("Epi_Cluster_", epi$seurat_clusters[match(sampled_epi_cells, colnames(epi))])
)
groupFiles <- "sampled_groupFiles.txt"
write.table(sampled_annotations_file, file = groupFiles, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

geneFile <- "geneFile.txt"
write.table(gene_location_ref, file = geneFile, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

T_ref_matrix <- T_expr_matrix[, ]

infercnv_obj <- CreateInfercnvObject(
    raw_counts_matrix = sampled_expression_matrix,
    annotations_file = groupFiles,
    delim = "\t",
    gene_order_file = geneFile,
    reference_group_matrix = T_ref_matrix  # 参考组仍为 T 细胞（需单独准备 T 细胞矩阵）
)

sampled_cells <- sample(colnames(expression_matrix), size = n_samples, replace = FALSE)

sampled_expression_matrix <- expression_matrix[, sampled_cells]
sampled_annotations_file <- annotations_file[annotations_file$cells %in% sampled_cells, ]

groupFiles <- "sampled_groupFiles.txt"
write.table(sampled_annotations_file, file = groupFiles, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

geneFile <- "geneFile.txt"
write.table(gene_location_ref, file = geneFile, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

infercnv_obj <- CreateInfercnvObject(
    raw_counts_matrix = sampled_expression_matrix,  # 直接使用表达矩阵
    annotations_file = groupFiles,  # 注释文件路径
    delim = "\t",  # 文件分隔符
    gene_order_file = geneFile,  # 基因位置参考文件路径
    ref_group_names = c("T cells")  # 以 T 细胞作为参考
)

infercnv_obj <- infercnv::run(
    infercnv_obj,
    cutoff = 0.2,  # 设置 cutoff 值
    out_dir = "infercnv_output",  # 输出目录
    cluster_by_groups = TRUE,  # 按组聚类
    denoise = TRUE,  # 去噪
    HMM = FALSE  # 是否使用 HMM 模型
)

#cnv score
cnv_scores <- infercnv_obj@expr.data
summary(as.vector(cnv_scores))
cnv_threshold <- quantile(cell_cnv_scores, probs = 0.90)
print(paste("CNV 得分 90% 分位数:", cnv_threshold))

malignant_cells <- names(cell_cnv_scores[cell_cnv_scores > cnv_threshold])
normal_cells <- names(cell_cnv_scores[cell_cnv_scores <= cnv_threshold])

print(paste("恶性细胞数量:", length(malignant_cells)))
print(paste("正常细胞数量:", length(normal_cells)))

malignant_epithelial_cells <- malignant_cells[malignant_cells %in% sampled_annotations_file$cells[sampled_annotations_file$cell_type == "Epithelial cells"]]
print(paste("恶性上皮细胞数量:", length(malignant_epithelial_cells)))
