# immune cells difference between LN metastasis and non-metastasis patients
## merge
```
library(Seurat)
seurat_obj1 <- readRDS("GSE243013_TN.rds")
seurat_obj1$orig.ident <- "GSE243013"
seurat_obj2 <- readRDS("GSE241934_TN.rds")
seurat_obj2$orig.ident <- "GSE241934"
merged_seurat_obj <- merge(seurat_obj1, y = c(seurat_obj2), add.cell.ids = c("GSE243013","GSE241934"))
merged_seurat_obj <- JoinLayers(merged_seurat_obj)
#31852 features across 292715 samples within 1 assay
counts_matrix <- GetAssayData(merged_seurat_obj, slot = "counts", assay = "RNA")
expressed_cells_per_gene <- Matrix::rowSums(counts_matrix > 0)

min_cell_percentage <- 0.001
min_cell_count <- 10
keep_genes <- expressed_cells_per_gene >= min_cell_percentage * ncol(merged_seurat_obj) & 
              expressed_cells_per_gene >= min_cell_count
merged_seurat_obj <- merged_seurat_obj[keep_genes, ]
#16227 features across 292715 samples within 1 assay

genes_per_cell <- merged_seurat_obj$nFeature_RNA
average_genes <- mean(genes_per_cell)
#1524

merged_seurat_obj@meta.data$Group <- ifelse(
  grepl("N0", merged_seurat_obj@meta.data$Stage), 
  "Group1", 
  "Group2"
)

table(merged_seurat_obj@meta.data$Stage, merged_seurat_obj@meta.data$Group)
#          Group1 Group2
#  T1bN1        0   6964
#  T1bN1M0      0   5644
#  T1bN2M0      0  29698
#  T1cN0M0   8154      0
#  T2aN0    13490      0
#  T2aN0M0  11251      0
#  T2bN0    12188      0
#  T2bN0M0  21596      0
#  T2N0M0   15483      0
#  T3N0     38859      0
#  T3N0M0   30978      0
#  T4N0     45155      0
#  T4N0M0   53255      0

library(dplyr)
if ("pathological_response" %in% colnames(merged_seurat_obj@meta.data)) {
    merged_seurat_obj@meta.data$Response <- ifelse(
      merged_seurat_obj@meta.data$pathological_response == "pCR", "YES", "NO"
    )
}

library(harmony)
merged_seurat_obj <- NormalizeData(merged_seurat_obj)
merged_seurat_obj <- FindVariableFeatures(merged_seurat_obj, nfeatures = 2000)

hvgs <- VariableFeatures(merged_seurat_obj)
merged_seurat_obj <- ScaleData(merged_seurat_obj, features = hvgs)
merged_seurat_obj <- RunPCA(merged_seurat_obj, features = hvgs, npcs = 20)
merged_seurat_obj <- RunHarmony(merged_seurat_obj, "orig.ident")

png("elbowplot.png", width = 800, height = 600)
elbowplot <- ElbowPlot(merged_seurat_obj)
print(elbowplot)
dev.off()

library(ggplot2)
library(ggsci)
merged_seurat_obj <- RunUMAP(merged_seurat_obj, reduction = "harmony", dims = 1:20)
merged_seurat_obj <- FindNeighbors(merged_seurat_obj, reduction = "harmony", dims = 1:20)
merged_seurat_obj <- FindClusters(merged_seurat_obj, resolution = 0.3)

npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(15)
seurat_clusters <- as.character(unique(merged_seurat_obj@meta.data$seurat_clusters))  # 转换为字符向量
num_legend_items <- length(seurat_clusters)  # 图例的个数
max_label_length <- max(nchar(seurat_clusters))  # 图例名称的最大长度

base_width <- 3000  # 基础宽度
base_height <- 3000  # 基础高度
legend_width_factor <- 100  # 每个图例项增加的宽度
label_length_factor <- 10  # 每个字符增加的宽度
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

markers <- FindAllMarkers(merged_seurat_obj, 
                          only.pos = TRUE, 
                          min.pct = 0.25, 
                          logfc.threshold = 0.25, 
                          test.use = "wilcox")
#significant markers
library(dplyr)
significant_markers <- subset(markers, p_val_adj < 0.05)
significant_markers <- significant_markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
write.csv(significant_markers,"marker_top_50.csv")
saveRDS(merged_seurat_obj,file="TN_umap.rds")


library(Seurat)
library(ggplot2)
library(ggsci)
library(dplyr)
merged_seurat_obj <- readRDS("TN_umap.rds")
identity_mapping <- c(
    "0" = "Tex",
    "1" = "Transitional_B",
    "2" = "Tn/Tm",
    "3" = "Treg_1",
    "4" = "NK",
    "5" = "Myeloid_1",
    "6" = "Myeloid_2",
    "7" = "Plasma",
    "8" = "Mast",
    "9" = "Germinal_Center_B ",
    "10" = "Proliferating",
    "11" = "DC",
    "12" = "Treg_2",
    "13" = "Neutrophil", 
    "14" = "pDC"
)

cell_type <- identity_mapping[merged_seurat_obj@meta.data$seurat_clusters]
merged_seurat_obj@meta.data$cell_type <- cell_type

npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(15)

cell_types <- unique(merged_seurat_obj@meta.data$cell_type)
num_legend_items <- length(cell_types)  # 图例的个数
max_label_length <- max(nchar(cell_types))  # 图例名称的最大长度

base_width <- 3000  # 基础宽度
base_height <- 3000  # 基础高度
legend_width_factor <- 100  # 每个图例项增加的宽度
label_length_factor <- 10  # 每个字符增加的宽度
dynamic_width <- base_width + (num_legend_items * legend_width_factor) + (max_label_length * label_length_factor)
pdf("annotation.pdf", width = dynamic_width/300, height = base_height/300)
DimPlot(merged_seurat_obj, reduction = "umap", label = TRUE, pt.size = 1, group.by = "cell_type", label.size = 4) +
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    ggtitle(NULL) +
    scale_color_manual(values = npg_extended) +
    coord_fixed(ratio = 1) +
    guides(color = guide_legend(title = NULL, override.aes = list(size = 5))) +
    theme(
        text = element_text(size = 8, face = "bold"),
        axis.text.x = element_text(size = 28, color = "black"),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title.x = element_text(size = 36, face = "bold", color = "black"),
        axis.title.y = element_text(size = 36, face = "bold", color = "black", margin = margin(r = 20)),  # 增加右侧间距
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 16, face = "bold", color = "black"),
        legend.title = element_text(size = 16, face = "bold", color = "black"),
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
saveRDS(merged_seurat_obj,file="TN_anno.rds")

```
## group information stastistic
```
plot_data <- merged_seurat_obj@meta.data %>%
  distinct(patients, cancer_type, Group, Response) %>%
  group_by(cancer_type, Group, Response) %>%
  summarise(patient_count = n(), .groups = 'drop')

plot_data <- plot_data %>%
  mutate(Group_Response = paste(Group, Response, sep = "-")) %>%
  mutate(Group_Response = factor(Group_Response,
                                levels = c("Group1-YES", "Group1-NO", 
                                         "Group2-YES", "Group2-NO")))

library(ggsci)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(10)

group_colors <- setNames(npg_extended, c("Group1-YES", "Group1-NO", "Group2-YES", "Group2-NO"))

print(plot_data)

pdf("cancer_type_group_response_patients_stacked.pdf", width = 4000/300, height = 3000/300)

ggplot(plot_data, aes(x = cancer_type, y = patient_count, fill = Group_Response)) +
  geom_bar(stat = "identity", position = "stack", width = 0.6) +
  scale_fill_manual(values = group_colors) +
  labs(
    x = "Cancer Type",
    y = "Number of Patients",
    fill = "Group-Response"
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 20, color = "black", face = "bold"),
    axis.text.y = element_text(size = 20, color = "black", face = "bold"),
    axis.title.x = element_text(size = 24, face = "bold", color = "black", margin = margin(t = 10)),
    axis.title.y = element_text(size = 24, face = "bold", color = "black", margin = margin(r = 20)),
    legend.text = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.position = "right",
    legend.key.size = unit(0.8, "cm"),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank()
  ) +
  # 添加数值标签
  geom_text(aes(label = patient_count), 
            position = position_stack(vjust = 0.5), 
            size = 6, 
            fontface = "bold",
            color = "white")

dev.off()
```
## LUAD and LUSC response difference
```

plot_data <- merged_seurat_obj

plot_data$cancer_response <- paste(plot_data$cancer_type, plot_data$Response, sep = "_")

plot_data$cancer_response <- factor(plot_data$cancer_response,
                                    levels = c("LUAD_YES", "LUAD_NO", "LUSC_YES", "LUSC_NO"))

npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(15)

pdf("combined_annotation_by_cancer_response.pdf", width = 12000/300, height = 3000/300)

p <- DimPlot(plot_data, 
             reduction = "umap", 
             label = TRUE, 
             pt.size = 1, 
             group.by = "cell_type", 
             label.size = 4,
             split.by = "cancer_response",  # 现在会按因子顺序分面
             ncol = 4) +  # 4列对应4个分组
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
                                 margin = margin(b = 15)),
        panel.spacing = unit(1.5, "lines")
    )

print(p)
dev.off()


plot_data <- merged_seurat_obj

npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(15)
proportion_data <- plot_data@meta.data %>%
    group_by(cancer_type, Response, cell_type) %>%
    summarise(count = n()) %>%
    mutate(proportion = count / sum(count)) %>%
    # 创建组合变量
    mutate(cancer_sample = paste0(cancer_type, "-", Response)) %>%
    # 设置组合变量的顺序
    mutate(cancer_sample = factor(cancer_sample,
                                 levels = c("LUAD-YES", "LUAD-NO",
                                           "LUSC-YES", "LUSC-NO")))

pdf("Response_cancer_type_pie.pdf", width = 8000/300, height = 3000/300) 

ggplot(proportion_data, aes(x = "", y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack", width = 1) +
  coord_polar(theta = "y") +  # 转换为饼图
  scale_fill_manual(values = npg_extended) +
  facet_wrap(~ cancer_sample, ncol = 4) +  # 4个饼图水平排列
  labs(x = "", y = "", fill = "Cell Type") +
  theme_void() +
  theme(
    legend.text = element_text(size = 28),
    legend.title = element_text(size = 28),
    strip.text = element_text(
      size = 20,
      face = "bold",
      margin = margin(b = 15)
    ),
    legend.position = "right"
  )

dev.off()


library(ggpubr)
library(dplyr)
library(ggsci)
library(tidyr)


current_cancer <- "LUAD"

cell_counts <- merged_seurat_obj@meta.data %>%
  filter(cancer_type == current_cancer) %>%
  group_by(patients, cell_type, Response) %>%
  summarise(count = n(), .groups = "drop")

sample_totals <- cell_counts %>%
  group_by(patients) %>%
  summarise(total = sum(count), .groups = "drop")

cell_proportions <- cell_counts %>%
  left_join(sample_totals, by = "patients") %>%
  mutate(proportion = count / total)

cell_proportions$Response <- factor(cell_proportions$Response, levels = c("YES", "NO"))

celltype_response_presence <- cell_proportions %>%
  group_by(cell_type, Response) %>%
  summarise(has_data = n() > 0, .groups = "drop") %>%
  pivot_wider(names_from = Response, values_from = has_data, values_fill = FALSE)

all_responses <- c("YES", "NO")

missing_in_NO <- celltype_response_presence %>%
  filter(YES == TRUE & NO == FALSE) %>%
  pull(cell_type)

missing_in_YES <- celltype_response_presence %>%
  filter(YES == FALSE & NO == TRUE) %>%
  pull(cell_type)

cat("在NO中缺失的细胞类型:", paste(missing_in_NO, collapse = ", "), "\n")
cat("在YES中缺失的细胞类型:", paste(missing_in_YES, collapse = ", "), "\n")


plot_data <- cell_proportions
if (length(missing_in_NO) > 0) {
  NO_patients <- unique(cell_proportions$patients[cell_proportions$Response == "NO"])
  missing_records_NO <- expand_grid(
    patients = NO_patients,
    cell_type = missing_in_NO,
    Response = "NO"
  ) %>%
  left_join(sample_totals, by = "patients") %>%
  mutate(count = 0, proportion = 0)
  
  plot_data <- plot_data %>% bind_rows(missing_records_NO)
}
if (length(missing_in_YES) > 0) {
  YES_patients <- unique(cell_proportions$patients[cell_proportions$Response == "YES"])
  missing_records_YES <- expand_grid(
    patients = YES_patients,
    cell_type = missing_in_YES,
    Response = "YES"
  ) %>%
  left_join(sample_totals, by = "patients") %>%
  mutate(count = 0, proportion = 0)
  
  plot_data <- plot_data %>% bind_rows(missing_records_YES)
}

all_celltypes <- unique(plot_data$cell_type)
plot_data$Response <- factor(plot_data$Response, levels = c("YES", "NO"))
full_combinations <- expand_grid(
  cell_type = all_celltypes,
  Response = factor(c("YES", "NO"), levels = c("YES", "NO"))
)

cell_response_counts <- cell_counts %>%
  group_by(cell_type, Response) %>%
  summarise(total_count = sum(count), .groups = "drop") %>%
  right_join(full_combinations, by = c("cell_type", "Response")) %>%
  mutate(total_count = replace(total_count, is.na(total_count), 0)) %>%
  arrange(cell_type, Response)

global_response_totals <- cell_counts %>%
  group_by(Response) %>%
  summarise(response_total = sum(count), .groups = "drop") %>%
  right_join(data.frame(Response = factor(c("YES", "NO"), levels = c("YES", "NO"))), by = "Response") %>%
  mutate(response_total = replace(response_total, is.na(response_total), 0)) %>%
  arrange(Response)

chisq_results <- cell_response_counts %>%
  group_by(cell_type) %>%
  summarise(
    p_value = {
      # 明确指定YES和NO的计数，不受顺序影响
      yes_count <- total_count[Response == "YES"]
      no_count <- total_count[Response == "NO"]
      yes_total <- global_response_totals$response_total[global_response_totals$Response == "YES"]
      no_total <- global_response_totals$response_total[global_response_totals$Response == "NO"]
      
      if (sum(c(yes_count, no_count)) < 5) {
        NA_real_
      } else {
        # 构建固定的2x2列联表
        # 第一行：YES组的该细胞类型计数 vs 其他细胞计数
        # 第二行：NO组的该细胞类型计数 vs 其他细胞计数
        cont_table <- matrix(c(yes_count, yes_total - yes_count,
                               no_count, no_total - no_count),
                             nrow = 2, byrow = FALSE)  # 固定按行填充
        
        if (any(cont_table < 0) || sum(cont_table) == 0) {
          NA_real_
        } else if (any(cont_table < 5)) {
          fisher.test(cont_table)$p.value
        } else {
          chisq.test(cont_table)$p.value
        }
      }
    },
    .groups = "drop"
  ) %>%
  mutate(
    significance = case_when(
      is.na(p_value) ~ "ns",
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**", 
      p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

# === 按首字母排序 ===
celltype_alphabetical <- sort(all_celltypes)
plot_data$cell_type <- factor(plot_data$cell_type, levels = celltype_alphabetical)
chisq_results$cell_type <- factor(chisq_results$cell_type, levels = celltype_alphabetical)

y_max <- max(plot_data$proportion, na.rm = TRUE) * 1.1

pdf("LUAD_response_boxplot.pdf", width = 6000/300, height = 3000/300)
ggplot(plot_data, aes(x = cell_type, y = proportion, fill = Response)) +
  geom_boxplot(width = 0.7) +
  geom_text(
    data = chisq_results, 
    aes(x = cell_type, y = y_max, label = significance),
    size = 12, 
    inherit.aes = FALSE
  ) +
  labs(x = "", y = "Cell Proportion", fill = "Response") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.y = element_text(size = 32),
    axis.line = element_line(color = "black"),
    panel.grid = element_blank(),
    legend.text = element_text(size = 36),
    legend.title = element_text(size = 40),
    legend.position = "right"
  ) +
  scale_fill_npg() +
  scale_y_continuous(limits = c(0, y_max))
dev.off()

```
# LUSC 2groups difference
```
LUSC <- subset(merged_seurat_obj, subset = cancer_type == "LUSC")
library(Seurat)
library(ggplot2)
library(ggsci)
library(dplyr)
plot_data <- LUSC

plot_data$group_response <- paste(plot_data$Group, plot_data$Response, sep = "_")


desired_order <- c("Group1_YES", "Group1_NO", 
                   "Group2_YES", "Group2_NO")

plot_data$group_response <- factor(plot_data$group_response,
                                   levels = desired_order)


print("存在的组合：")
print(table(plot_data$group_response))

npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(15)
pdf("LUSC_combined_annotation_by_group_response.pdf", width = 12000/300, height = 3000/300)

p <- DimPlot(plot_data, 
             reduction = "umap", 
             label = TRUE, 
             pt.size = 1, 
             group.by = "cell_type", 
             label.size = 4,
             split.by = "group_response",  # 现在会按因子顺序分面
             ncol = 4) +  # 固定4列：Group1_YES, Group1_NO, Group2_YES, Group2_NO
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
                                 margin = margin(b = 15)),
        panel.spacing = unit(1.5, "lines")
    )

print(p)
dev.off()


plot_data <- LUSC

npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(15)
proportion_data <- plot_data@meta.data %>%
    group_by(Group, Response, cell_type) %>%
    summarise(count = n()) %>%
    mutate(proportion = count / sum(count)) %>%
    # 创建组合变量
    mutate(group_sample = paste0(Group, "-", Response)) %>%
    # 设置组合变量的顺序
    mutate(group_sample = factor(group_sample,
                                 levels = c("Group1-YES", "Group1-NO",
                                           "Group2-YES", "Group2-NO")))

pdf("LUSC_Response_group_stacked.pdf", width = 8000/300, height = 3000/300) 

ggplot(proportion_data, aes(x = group_sample, y = proportion, fill = cell_type)) +
  scale_fill_manual(values = npg_extended) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "", y = "Proportion", fill = "Cell Type") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 18, angle = 0, hjust = 1, vjust = 0.5),
    axis.line = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 28),
    legend.title = element_text(size = 28),
    axis.text.y = element_text(size = 18),
    axis.title.y = element_text(size = 20)
  )

dev.off()


library(ggpubr)
library(dplyr)
library(ggsci)

cell_counts <- LUSC@meta.data %>%
  group_by(patients, cell_type, Group, Response) %>%
  summarise(count = n()) %>%
  ungroup()

# 计算每个样本的总细胞数
sample_totals <- cell_counts %>%
  group_by(patients) %>%
  summarise(total = sum(count)) %>%
  ungroup()

# 计算细胞类型比例
cell_proportions <- LUSC@meta.data %>%
  group_by(patients, cell_type, Group, Response) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  left_join(sample_totals, by = "patients") %>%
  mutate(proportion = count / total) %>%
  # 创建组合变量
  mutate(group_response = paste0(Group, "-", Response)) %>%
  # 设置组合变量的顺序
  mutate(group_response = factor(group_response,
                                 levels = c("Group1-YES", "Group1-NO",
                                           "Group2-YES", "Group2-NO")))

# === 修改卡方检验：比较四个癌症-应答组的比例差异 ===
chisq_results <- cell_proportions %>%
  group_by(cell_type) %>%
  summarise(
    p_value = {
      # 构建每个细胞类型在四个癌症-应答组中的观察值
      obs_counts <- c(
        sum(count[group_response == "Group1-YES"]),
        sum(count[group_response == "Group1-NO"]), 
        sum(count[group_response == "Group2-YES"]),
        sum(count[group_response == "Group2-NO"])
      )
      
      # 计算四个癌症-应答组的总细胞数（作为期望值的基础）
      total_counts <- c(
        sum(total[group_response == "Group1-YES"]),
        sum(total[group_response == "Group1-NO"]),
        sum(total[group_response == "Group2-YES"]),
        sum(total[group_response == "Group2-NO"])
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

# === 排序部分 ===
celltype_by_proportion <- cell_proportions %>%
  group_by(cell_type) %>%
  summarise(mean_prop = mean(proportion)) %>%
  arrange(desc(mean_prop)) %>%
  pull(cell_type)

# 设置因子顺序
cell_proportions$cell_type <- factor(cell_proportions$cell_type, 
                                    levels = celltype_by_proportion)

if ("cell_type" %in% colnames(chisq_results)) {
  chisq_results$cell_type <- factor(chisq_results$cell_type, 
                                   levels = celltype_by_proportion)
}

# 绘制图形
pdf("LUSC_box_proportion_group_response.pdf", width = 6000/300, height = 3000/300)
ggplot(cell_proportions, aes(x = cell_type, y = proportion, fill = group_response)) +
  geom_boxplot() +
  geom_text(
    data = chisq_results, 
    aes(x = cell_type, y = max(cell_proportions$proportion, na.rm = TRUE) * 1.05, 
        label = significance),
    size = 12, 
    vjust = 0.5,
    inherit.aes = FALSE
  ) +
  labs(x = "", y = "Cell Proportion", fill = "Group-Response") +
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
  scale_y_continuous(limits = c(0, 0.65))
dev.off()


#折线图
plot_data <- LUSC@meta.data %>%
  group_by(Group, Response, cell_type) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  # 计算每个Group-Response组合内的比例
  group_by(Group, Response) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()

# 绘制折线图 - 去掉分面标题方框
pdf("LUSC_group_response_lineplot.pdf", width = 8000/300, height = 4000/300)

ggplot(plot_data, aes(x = Response, y = proportion, group = cell_type, color = cell_type)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) +
  facet_wrap(~ Group, ncol = 2) +
  labs(
    x = "Response",
    y = "Cell Proportion",
    color = "Cell Type"
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 20, color = "black", face = "bold"),
    axis.text.y = element_text(size = 20, color = "black", face = "bold"),
    axis.title.x = element_text(size = 24, face = "bold", color = "black", margin = margin(t = 10)),
    axis.title.y = element_text(size = 24, face = "bold", color = "black", margin = margin(r = 20)),
    legend.text = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.position = "right",
    legend.key.size = unit(0.8, "cm"),
    # 修改分面标题样式 - 去掉方框
    strip.background = element_blank(),  # 去掉背景
    strip.text = element_text(size = 20, face = "bold", margin = margin(b = 10)),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.spacing = unit(1.5, "lines")
  ) +
  scale_color_manual(values = npg_extended) +
  scale_y_continuous(labels = scales::percent_format())

dev.off()
```

