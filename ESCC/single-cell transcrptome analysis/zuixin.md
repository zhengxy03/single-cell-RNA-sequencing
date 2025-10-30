# 导入 GSE160269.rds 和 PRJNA77911.rds
# 下载淋巴结相关数据
```bash
cd GSE203067
cat <<EOF > source.csv
SRR19216395,,HiSeq X Ten
SRR19216396,,HiSeq X Ten
SRR19216397,,HiSeq X Ten
SRR19216398,,HiSeq X Ten
SRR19216399,,HiSeq X Ten
SRR19216400,,HiSeq X Ten
SRR19216401,,HiSeq X Ten
EOF

anchr ena info | perl - -v source.csv > ena_info.yml
anchr ena prep | perl - ena_info.yml
mlr --itsv --ocsv cat ena_info.tsv > ena_info.csv

mlr --icsv --omd cat ena_info.csv

aria2c -x 4 -s 2 -c -i ena_info.ftp.txt

md5sum --check ena_info.md5.txt
```
# 上传淋巴结相关数据到超算并转换格式/重命名
```
bsub -q mpi -n 24 -J 97 "fastq-dump --gzip --split-files SRR19216397.lite.1"
```
# cell ranger处理
```
#!/bin/bash
# 定义样本列表
samples=(
  "SRR19216396"
)

# 定义公共参数
transcriptome="/share/home/wangq/zxy/ESCC/refgenome/refdata-gex-GRCh38-2024-A"
fastqs="/share/home/wangq/zxy/GSE203067"
nosecondary="--nosecondary"
create_bam="--create-bam=false"
localcores=24  # 设置使用的核心数

# 循环处理每个样本
for sample in "${samples[@]}"; do
    job_script="job_${sample}.sh"
    
    cat << EOF > $job_script
#!/bin/bash
#BSUB -q mpi
#BSUB -n $localcores  # 使用指定数量的核心
#BSUB -o ${sample}_%J.out
#BSUB -e ${sample}_%J.err
cellranger count --id=$sample \\
                 --transcriptome=$transcriptome \\
                 --fastqs=$fastqs \\
                 --sample=$sample \\
                 $nosecondary \\
                 $create_bam \\
                 --localcores=$localcores
EOF

    # 提交作业
    bsub < $job_script
done
```
# 下游分析
```
library(Seurat)
library(ggplot2)
library(dplyr)

setwd("/share/home/wangq/zxy/GSE203067/")

metadata <- read.csv("数据集导入.csv")

samples <- metadata$sample
seurat_list <- list()

for (sample in samples) {
  data_dir <- paste0("./", sample, "/outs/filtered_feature_bc_matrix")
  cat("样本:", sample, "\n")
  cat("目录路径:", data_dir, "\n")
  cat("目录存在:", dir.exists(data_dir), "\n")
  data <- Read10X(data.dir = data_dir)
  seurat_obj <- CreateSeuratObject(counts = data, project = sample)

  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA >= 400 & nFeature_RNA <= 7500 & 
                         percent.mt < 10 & 
                         nCount_RNA >= 500 & nCount_RNA <= 50000)

  seurat_list[[sample]] <- seurat_obj
}

merged_seurat_obj <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = samples)

sample_metadata <- metadata[match(samples, metadata$sample), ]
sample_names <- merged_seurat_obj@meta.data$orig.ident
cell_metadata <- sample_metadata[match(sample_names, sample_metadata$sample), ]
merged_seurat_obj@meta.data <- cbind(merged_seurat_obj@meta.data, cell_metadata)
merged_seurat_obj <- JoinLayers(merged_seurat_obj)
head(merged_seurat_obj@meta.data)
saveRDS(merged_seurat_obj, file = "GSE203067.rds")
```
# 加入patients元数据
```
merged_seurat_obj <- merge(GSE203067, 
                   y = c(GSE160269, PRJNA777911))
merged_seurat_obj <- JoinLayers(merged_seurat_obj)
patients_df <- read.csv("patients.csv")
print(head(patients_df))

merged_seurat_obj$patients <- patients_df$patients[match(merged_seurat_obj$sample, patients_df$sample)]
saveRDS(merged_seurat_obj, file = "merged_origin.rds")

```
# 去除双细胞、批次效应
```
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

library(DoubletFinder)
sweep.res.list <- paramSweep(merged_seurat_obj, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
merged_seurat_obj$DF.classifications <- NULL
pK_value <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))

# 分批处理
batch_size <- 10000
num_cells <- ncol(merged_seurat_obj)
num_batches <- ceiling(num_cells / batch_size)

batch_results <- list()

for (i in 1:num_batches) {
  start_idx <- (i - 1) * batch_size + 1
  end_idx <- min(i * batch_size, num_cells)
  
  # 提取当前批次的细胞
  batch_cells <- subset(merged_seurat_obj, cells = colnames(merged_seurat_obj)[start_idx:end_idx])
  
  # 检测双细胞
  nExp_poi <- round(0.075 * ncol(batch_cells))  # 假设双细胞比例为 7.5%
  batch_cells <- doubletFinder(batch_cells, PCs = 1:20, pN = 0.25, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  
  # 提取双细胞检测结果
  doublet_col <- grep("DF.classifications", colnames(batch_cells@meta.data), value = TRUE)
  batch_results[[i]] <- batch_cells@meta.data[[doublet_col]]
  
  # 释放内存
  rm(batch_cells)
  gc()
}

all_results <- unlist(batch_results)

# 将双细胞检测结果添加到原始的 Seurat 对象中
merged_seurat_obj@meta.data$DF.classifications <- all_results

# 去除双细胞
merged_seurat_obj <- subset(merged_seurat_obj, subset = DF.classifications == "Singlet")

#average genes per cell
genes_per_cell <- merged_seurat_obj$nFeature_RNA
average_genes <- mean(genes_per_cell)

library(harmony)
merged_seurat_obj <- RunHarmony(merged_seurat_obj, "orig.ident")


samples <- unique(merged_seurat_obj@meta.data$orig.ident)
sample_numbers <- as.numeric(gsub("[^0-9]", "", samples))
samples_ordered <- samples[order(sample_numbers)]
non_numeric_samples <- samples[is.na(sample_numbers)]  # 提取没有数字部分的样本名称
samples_ordered <- c(samples_ordered, non_numeric_samples)

merged_seurat_obj@meta.data$orig.ident <- factor(
  merged_seurat_obj@meta.data$orig.ident,
  levels = samples_ordered  # 按数字顺序设置水平
)


# 保存 merged_seurat_obj 对象
saveRDS(merged_seurat_obj, file = "merged_pre.rds")

merged_seurat_obj <- FindNeighbors(merged_seurat_obj, reduction = "harmony", dims = 1:20)

library(ggplot2)
library(ggsci)
merged_seurat_obj <- FindClusters(merged_seurat_obj, resolution = 0.2)
merged_seurat_obj <- RunUMAP(merged_seurat_obj, reduction = "harmony", dims = 1:20)
saveRDS(merged_seurat_obj, file = "merged_umap.rds")

npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(18)

# 获取图例的个数和名称长度
seurat_clusters <- as.character(unique(merged_seurat_obj@meta.data$seurat_clusters))  # 转换为字符向量
num_legend_items <- length(seurat_clusters)  # 图例的个数
max_label_length <- max(nchar(seurat_clusters))  # 图例名称的最大长度

# 动态计算图片尺寸
base_width <- 3000  # 基础宽度
base_height <- 3000  # 基础高度
legend_width_factor <- 100  # 每个图例项增加的宽度
label_length_factor <- 10  # 每个字符增加的宽度

# 计算动态宽度
dynamic_width <- base_width + (num_legend_items * legend_width_factor) + (max_label_length * label_length_factor)
pdf("clusters.pdf", width = dynamic_width/300, height = base_height/300)  # 转换为英寸

DimPlot(merged_seurat_obj, reduction = "umap", label = TRUE, pt.size = 1, label.size = 8) +
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    ggtitle(NULL) +
    scale_color_manual(values = npg_extended) +
    coord_fixed(ratio = 1) +
    guides(color = guide_legend(title = NULL, override.aes = list(size = 5))) +
    theme(
        text = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 28, color = "black"),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title.x = element_text(size = 36, face = "bold", color = "black"),
        axis.title.y = element_text(size = 36, face = "bold", color = "black", margin = margin(r = 20)),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 28, face = "bold", color = "black"),
        legend.title = element_text(size = 28, face = "bold", color = "black"),
        legend.position = "right",
        legend.box.margin = margin(0, 0, 0, 0),
        legend.key = element_blank(),
        legend.background = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5),
        aspect.ratio = 1,
        plot.margin = margin(10, 50, 10, 10)
    )

dev.off()

identity_mapping <- c(
    "0" = "T cell",
    "1" = "T cell",
    "2" = "B cell",
    "3" = "Epithelial cell",
    "4" = "Fibroblast",
    "5" = "Macrophage",
    "6" = "Endothelial cell",
    "7" = "Plasma",
    "8" = "Mast cell", 
    "9" = "Proliferating cell",
    "10" = "Pericyte",
    "11" = "Schwann",
    "12" = "Epithelial cell",
    "13" = "T cell"
)

cell_type <- identity_mapping[merged_seurat_obj@meta.data$seurat_clusters]
merged_seurat_obj@meta.data$cell_type <- cell_type
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(14)

# 获取图例的个数和名称长度
cell_types <- unique(merged_seurat_obj@meta.data$cell_type)
num_legend_items <- length(cell_types)  # 图例的个数
max_label_length <- max(nchar(cell_types))  # 图例名称的最大长度

# 动态计算图片尺寸
base_width <- 3000  # 基础宽度
base_height <- 3000  # 基础高度
legend_width_factor <- 100  # 每个图例项增加的宽度
label_length_factor <- 10  # 每个字符增加的宽度

# 计算动态宽度
dynamic_width <- base_width + (num_legend_items * legend_width_factor) + (max_label_length * label_length_factor)

# 导出图片
pdf("annotation.pdf", width = dynamic_width/300, height = base_height/300)
DimPlot(merged_seurat_obj, reduction = "umap", label = TRUE, pt.size = 1, group.by = "cell_type", label.size = 6) +
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    ggtitle(NULL) +
    scale_color_manual(values = npg_extended) +
    coord_fixed(ratio = 1) +
    guides(color = guide_legend(title = NULL, override.aes = list(size = 5))) +
    theme(
        text = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 28, color = "black"),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title.x = element_text(size = 36, face = "bold", color = "black"),
        axis.title.y = element_text(size = 36, face = "bold", color = "black", margin = margin(r = 20)),  # 增加右侧间距
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 28, face = "bold", color = "black"),
        legend.title = element_text(size = 28, face = "bold", color = "black"),
        legend.position = "right",
        legend.box.margin = margin(0, 0, 0, 0),
        legend.key = element_blank(),
        legend.background = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5),
        aspect.ratio = 1,
        plot.margin = margin(10, 50, 10, 10)
    )
dev.off()

#DEGs

genes_to_plot <- c(
    # 淋巴细胞系
    "CD3D", "CD3E", "NKG7",          # T细胞
    "BANK1", "CD79A", "MS4A1",        # B细胞

    "IGHG1", "IGHG3", "JCHAIN",      # 浆细胞
    
    # 髓系免疫细胞
    "C1QA", "AIF1", "LYZ",         # 巨噬细胞

    
    # 结构细胞
    "SFN", "KRT19", "KRT17",       # 上皮细胞
    "SFRP2", "COL1A2", "FGF7",     # 成纤维细胞
    "CLDN5", "PECAM1", "RAMP2",      # 内皮细胞
    "RGS5", "MCAM", "ACTA2",         # 周细胞
    
    # 特殊功能细胞
    "CPA3", "TPSAB1", "TPSB2",       # 肥大细胞
    "TOP2A", "STMN1", "MKI67",         # 增殖细胞
    "S100B", "SOX10", "SOX14"
)

pdf("dotplot.pdf", width = 8000/300, height = 3000/300)  # 设置高分辨率和尺寸
DotPlot(merged_seurat_obj, 
        features = genes_to_plot, 
        group.by = "cell_type",
        dot.scale = 15) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 28),  # 横坐标基因名字体大小
    axis.text.y = element_text(size = 28),  # 纵坐标细胞类型字体大小
    axis.title.x = element_blank(),  # 去除横坐标标题
    axis.title.y = element_text(size = 40),
    panel.grid.major = element_line(color = "grey70", linewidth = 0.5),  # 添加网格线
    panel.grid.minor = element_line(color = "grey80", linewidth = 0.3),  # 细网格线 
    legend.text = element_text(size = 28, face = "bold", color = "black"),
    legend.title = element_text(size = 28, face = "bold", color = "black"),
    legend.position = "right",
    legend.box.margin = margin(0, 0, 0, 0),
    legend.key = element_blank(),
    legend.background = element_blank()
  ) +
  scale_color_gsea()  # 使用ggsci的GSEA配色
dev.off()
```
# 不同组织中的细胞分布及比例
```
#plot celltype according to sample types

pdf("annotation_by_sampletype.pdf", width = 8000/300, height = 3000/300)

p <- DimPlot(merged_seurat_obj, 
             reduction = "umap", 
             label = TRUE, 
             pt.size = 1, 
             group.by = "cell_type", 
             label.size = 4,
             split.by = "sample_type",
             ncol = 4) +
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    ggtitle(NULL) +
    scale_color_manual(values = npg_extended) +
    coord_fixed(ratio = 1) +
    guides(color = guide_legend(title = NULL, override.aes = list(size = 5))) +
    theme(
        text = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 20, color = "black"),
        axis.text.y = element_text(size = 20, color = "black"),
        axis.title.x = element_text(size = 24, face = "bold", color = "black"),
        axis.title.y = element_text(size = 24, face = "bold", color = "black", margin = margin(r = 20)),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 20, face = "bold", color = "black"),
        legend.title = element_text(size = 20, face = "bold", color = "black"),
        legend.position = "right",
        legend.box.margin = margin(0, 0, 0, 0),
        legend.key = element_blank(),
        legend.background = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5),
        aspect.ratio = 1,
        plot.margin = margin(10, 50, 10, 10),
        strip.text = element_text(size = 16, face = "bold", 
                                 margin = margin(b = 15)),  # 增加分面标题下方的间距
        panel.spacing = unit(1.5, "lines")  # 增加分面面板之间的间距
    )

print(p)
dev.off()

#cell proportion based on sampletype
proportion_data <- merged_seurat_obj@meta.data %>%
    group_by(sample_type, cell_type) %>% summarise(count = n()) %>% mutate(proportion = count / sum(count))

npg_extended <- colorRampPalette(npg_pal)(14)
pdf("sampletype_prop.pdf", width = 6000/300, height = 3000/300)  # 设置高分辨率和尺寸
ggplot(proportion_data, aes(x = "", y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", width = 1) +         # 堆叠柱状图
  coord_polar(theta = "y") +                       # 转换为饼图
  scale_fill_manual(values = npg_extended) +       # 使用自定义颜色
  theme_void() +                                   # 空白背景
  labs(fill = "Cell Type") +
  theme(
    legend.position = "right",                     # 图例放在右侧
    plot.title = element_blank(),                  # 移除标题
    # 分面标签设置
    strip.placement = "outside",                   # 标签放在绘图区域外
    strip.text = element_text(                     # 标签样式
      size = 40,                                   # 字体大小
      face = "bold",                               # 加粗
      margin = margin(b = 10)                      # 下方留白（避免与饼图重叠）
    ),
    legend.text = element_text(size = 36),         # 图例文本字体大小
    legend.title = element_text(size = 36)                 # 移除图例标题
  ) +
  facet_wrap(
    ~ sample_type,
    ncol = 2,
    strip.position = "bottom"                      # 标签放在下方
  )
dev.off()

library(ggpubr)
library(dplyr)
library(ggsci)

cell_counts <- merged_seurat_obj@meta.data %>%
  group_by(patients, cell_type, sample_type) %>%
  summarise(count = n()) %>%
  ungroup()

# 计算每个样本的总细胞数
sample_totals <- cell_counts %>%
  group_by(patients) %>%
  summarise(total = sum(count)) %>%
  ungroup()

# 计算细胞类型比例
cell_proportions <- merged_seurat_obj@meta.data %>%
  group_by(patients, cell_type, sample_type) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  left_join(sample_totals, by = "patients") %>%
  mutate(proportion = count / total)

# === 修改：使用卡方检验比较四个样本类型的比例差异 ===
chisq_results <- cell_proportions %>%
  group_by(cell_type) %>%
  summarise(
    p_value = {
      # 构建每个细胞类型在四个样本类型中的观察值
      obs_counts <- c(
        sum(count[sample_type == "Tumor"]),
        sum(count[sample_type == "Normal"]), 
        sum(count[sample_type == "mLN"]),
        sum(count[sample_type == "nLN"])
      )
      
      # 计算四个样本类型的总细胞数（作为期望值的基础）
      total_counts <- c(
        sum(total[sample_type == "Tumor"]),
        sum(total[sample_type == "Normal"]),
        sum(total[sample_type == "mLN"]),
        sum(total[sample_type == "nLN"])
      )
      
      # 如果任何组的计数为0或总计数太小，返回NA
      if (sum(obs_counts) < 10 || any(obs_counts == 0)) {
        NA_real_
      } else {
        # 使用卡方检验比较观察比例与期望比例
        # 期望值基于总细胞数的分布
        expected_prop <- total_counts / sum(total_counts)
        chisq.test(obs_counts, p = expected_prop)$p.value
      }
    }
  ) %>%
  ungroup()

# 添加显著性标记
chisq_results <- chisq_results %>%
  mutate(
    significance = case_when(
      is.na(p_value) ~ "ns",
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**", 
      p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

# === 排序部分保持不变 ===
celltype_by_proportion <- cell_proportions %>%
  group_by(cell_type) %>%
  summarise(mean_prop = mean(proportion)) %>%
  arrange(desc(mean_prop)) %>%
  pull(cell_type)

sample_type_order <- c("Tumor", "Normal", "mLN", "nLN")

cell_proportions$cell_type <- factor(cell_proportions$cell_type, 
                                    levels = celltype_by_proportion)
cell_proportions$sample_type <- factor(cell_proportions$sample_type, 
                                      levels = sample_type_order)

if ("cell_type" %in% colnames(chisq_results)) {
  chisq_results$cell_type <- factor(chisq_results$cell_type, 
                                   levels = celltype_by_proportion)
}

# 绘制图形
pdf("箱图_proportion.pdf", width = 6000/300, height = 3000/300)
ggplot(cell_proportions, aes(x = cell_type, y = proportion, fill = sample_type)) +
  geom_boxplot() +
  geom_text(
    data = chisq_results, 
    aes(x = cell_type, y = max(cell_proportions$proportion, na.rm = TRUE) * 1.05, 
        label = significance),
    size = 12, 
    vjust = 0.5,
    inherit.aes = FALSE
  ) +
  labs(x = "", y = "Cell Proportion", fill = "Sample Type") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 28),
    axis.text.y = element_text(size = 28),
    axis.title.y = element_text(size = 40),
    axis.line = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 36),
    legend.title = element_text(size = 40),
    legend.position = "right"
  ) +
  scale_fill_npg() +
  scale_y_continuous(limits = c(0, 1))
dev.off()
```
# 分期中细胞比例变化
```
table(merged_seurat_obj@meta.data$period1)

# 将“ⅡⅠ”和“Ⅲ”统一合并为“Ⅲ”（或自定义名称如“三期”）
merged_seurat_obj@meta.data$period1 <- ifelse(
  merged_seurat_obj@meta.data$period1 %in% c("ⅡⅠ", "Ⅲ"),  # 包含两个三期
  "Ⅲ",  # 合并后统一为“Ⅲ”
  merged_seurat_obj@meta.data$period1  # 其他期别不变
)

# 检查合并结果
table(merged_seurat_obj@meta.data$period1)

merged_seurat_obj@meta.data$period1[merged_seurat_obj@meta.data$sample_type == "nLN"] <- "nLN"
merged_seurat_obj@meta.data$period1[merged_seurat_obj@meta.data$sample_type == "mLN"] <- "mLN"
saveRDS(merged_seurat_obj, file = "merged_anno_period.rds")

npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(14)
proportion_data <- merged_seurat_obj@meta.data %>%
    group_by(period1, sample_type, cell_type) %>% 
    summarise(count = n()) %>% 
    mutate(proportion = count / sum(count)) %>%
    # 转换罗马数字为阿拉伯数字
    mutate(period1 = case_when(
        period1 == "Ⅰ" ~ "1",
        period1 == "Ⅱ" ~ "2", 
        period1 == "Ⅲ" ~ "3",
        TRUE ~ as.character(period1)
    ))


# Tumor + mLN 图
tumor_data <- proportion_data %>% filter(sample_type == "Tumor")
pdf("period1_line_tumor.pdf", width = 4000/300, height = 3000/300)
ggplot(tumor_data, aes(x = period1, y = proportion, group = cell_type, color = cell_type)) +
    geom_line(linewidth = 1.5) +
    geom_point(size = 3) +
    scale_color_manual(values = npg_extended) +
    labs(x = "Period", y = "Cell Proportion", color = "Cell Type") +
    theme_classic() +
    theme(
        axis.text.x = element_text(size = 28),
        axis.text.y = element_text(size = 28),
        axis.title.x = element_text(size = 36),
        axis.title.y = element_text(size = 36),
        legend.text = element_text(size = 28),
        legend.title = element_text(size = 32)
    ) +
    scale_y_continuous(labels = scales::percent)
dev.off()

normal_data <- proportion_data %>% filter(sample_type == "Normal")
pdf("period1_line_Normal.pdf", width = 4000/300, height = 3000/300)
ggplot(normal_data, aes(x = period1, y = proportion, group = cell_type, color = cell_type)) +
    geom_line(linewidth = 1.5) +
    geom_point(size = 3) +
    scale_color_manual(values = npg_extended) +
    labs(x = "Period", y = "Cell Proportion", color = "Cell Type") +
    theme_classic() +
    theme(
        axis.text.x = element_text(size = 28),
        axis.text.y = element_text(size = 28),
        axis.title.x = element_text(size = 36),
        axis.title.y = element_text(size = 36),
        legend.text = element_text(size = 28),
        legend.title = element_text(size = 32)
    ) +
    scale_y_continuous(labels = scales::percent)
dev.off()

```



```

epi_clean <- subset(epi, subset = seurat_clusters %in% c(0,1,2,3))
epi_clean <- subset(epi_clean,subset = seurat_clusters %in% c(0,1,2,3,5,7,8,10))

epi_clean <- NormalizeData(epi_clean)
epi_clean <- FindVariableFeatures(epi_clean, nfeatures = 2000)
hvgs <- VariableFeatures(epi_clean)
epi_clean <- ScaleData(epi_clean, features = hvgs)
epi_clean <- RunPCA(epi_clean, features = hvgs, npcs = 20)
library(harmony)
epi_clean <- RunHarmony(epi_clean, "sample_sources")
epi_clean <- RunUMAP(epi_clean, dims = 1:20, reduction = "harmony")
epi_clean <- FindNeighbors(epi_clean, dims = 1:20, reduction = "harmony")
epi_clean <- FindClusters(epi_clean, resolution = 0.1)

seurat_clusters <- as.character(unique(epi_clean@meta.data$seurat_clusters))
# 计算图例的个数
num_legend_items <- length(seurat_clusters)
# 计算图例名称的最大长度
max_label_length <- max(nchar(seurat_clusters))

# 动态计算图片尺寸
base_width <- 3000  # 基础宽度
base_height <- 3000  # 基础高度
legend_width_factor <- 100  # 每个图例项增加的宽度
label_length_factor <- 10  # 每个字符增加的宽度

# 计算动态宽度
dynamic_width <- base_width + (num_legend_items * legend_width_factor) + (max_label_length * label_length_factor)

# 创建 NPG 颜色调色板
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(13)

# 打开 PNG 图形设备，设置图片尺寸和分辨率
pdf("epi_clusters.pdf", width = dynamic_width/300, height = base_height/300)

# 绘制 UMAP 图
DimPlot(epi_clean, reduction = "umap", label = TRUE, pt.size = 1, label.size = 8) +
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    ggtitle(NULL) +
    scale_color_manual(values = npg_extended) +
    coord_fixed(ratio = 1) +
    guides(color = guide_legend(title = NULL, override.aes = list(size = 5))) +
    theme(
        text = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 28, color = "black"),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title.x = element_text(size = 36, face = "bold", color = "black"),
        axis.title.y = element_text(size = 36, face = "bold", color = "black", margin = margin(r = 20)),  # 增加右侧间距
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 28, face = "bold", color = "black"),
        legend.title = element_text(size = 28, face = "bold", color = "black"),
        legend.position = "right",
        legend.box.margin = margin(0, 0, 0, 0),
        legend.key = element_blank(),
        legend.background = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5),
        aspect.ratio = 1,
        plot.margin = margin(10, 50, 10, 10)
    )

# 关闭图形设备
dev.off()

library(dplyr)

epi_markers <- FindAllMarkers(epi_clean, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "wilcox")
epi_significant_markers <- subset(epi_markers, p_val_adj < 0.05)
write.csv(epi_significant_markers, "epi_all_marker.csv")
epi_significant_markers <- epi_significant_markers %>% 
  group_by(cluster) %>% 
  top_n(n = 50, wt = avg_log2FC)
write.csv(epi_significant_markers, "epi_top_marker.csv")

identity_mapping <- c(
    "0" = "Proliferative epi",
    "1" = "Basal stem epi",
    "2" = "Differentiated epi",
    "3" = "Differentiated epi",
 
    "5" = "Inflamed epi", #应激响应/炎症状态的上皮细胞
    
    "7" = "Basal stem epi",    
    "8" = "Inflamed epi",

    "10" = "Inflamed epi"
)


identity_mapping <- c(
    "0" = "Proliferative epi",
    "1" = "Basal stem epi",
    "2" = "Differentiated epi",
    "3" = "Differentiated epi",
 
    "5" = "Inflamed epi", #应激响应/炎症状态的上皮细胞
    
    "7" = "Basal stem epi",    
    "8" = "Inflamed epi",

    "10" = "Inflamed epi"
)