```
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggsci)
merged_seurat_obj <- readRDS("TN_anno_new.rds")
myeloid <- subset(merged_seurat_obj,subset=cell_type %in% c("cDC1","Macrophage","pDC","Monocyte"))

myeloid <- NormalizeData(myeloid)
myeloid <- FindVariableFeatures(myeloid, nfeatures = 2000)
hvgs <- VariableFeatures(myeloid)
myeloid <- ScaleData(myeloid, features = hvgs)
myeloid <- RunPCA(myeloid, features = hvgs, npcs = 20)
#16227 features across 19200 samples within 1 assay
#only one dataset
#library(harmony)
#t_cells <- RunHarmony(t_cells, "orig.ident")
myeloid <- RunUMAP(myeloid, dims = 1:20)
myeloid <- FindNeighbors(myeloid, dims = 1:20)
myeloid <- FindClusters(myeloid, resolution = 0.3)
myeloid_markers <- FindAllMarkers(myeloid, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "wilcox")
myeloid_significant_markers <- subset(myeloid_markers, p_val_adj < 0.05)
#write.csv(myeloid_significant_markers, "myeloid_all_marker.csv")
myeloid_significant_markers <- myeloid_significant_markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.csv(myeloid_significant_markers, "myeloid_top_marker_20_filter.csv")
saveRDS(myeloid,file="myeloid.rds")

seurat_clusters <- as.character(unique(myeloid@meta.data$seurat_clusters))  # 转换为字符向量
num_legend_items <- length(seurat_clusters)  # 图例的个数
max_label_length <- max(nchar(seurat_clusters))  # 图例名称的最大长度

# 动态计算图片尺寸
base_width <- 3000  # 基础宽度
base_height <- 3000  # 基础高度
legend_width_factor <- 100  # 每个图例项增加的宽度
label_length_factor <- 10  # 每个字符增加的宽度

# 计算动态宽度
dynamic_width <- base_width + (num_legend_items * legend_width_factor) + (max_label_length * label_length_factor)

npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(16)
pdf("myeloid_clusters.pdf", width = dynamic_width/300, height = base_height/300)
DimPlot(myeloid_filter, reduction = "umap", label = TRUE, pt.size = 1, label.size = 8) +
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



library(ggplot2)
target_genes <- c("CD1C","CLEC9A","ZBTB46")
feature_plot <- FeaturePlot(myeloid, features = target_genes)
ggsave("myeloid_featureplot.png", plot = feature_plot, width = 6 , height = 8, dpi = 300)


myeloid_filter <- subset(myeloid,subset=seurat_clusters %in% c(0,1,2,4,5,6,8,9,10,11,12,13,14,15))
myeloid <- NormalizeData(myeloid)
myeloid <- FindVariableFeatures(myeloid, nfeatures = 2000)
hvgs <- VariableFeatures(myeloid)
myeloid <- ScaleData(myeloid, features = hvgs)
myeloid <- RunPCA(myeloid, features = hvgs, npcs = 20)
#16227 features across 16636 samples within 1 assa
#only one dataset
#library(harmony)
#t_cells <- RunHarmony(t_cells, "orig.ident")
myeloid <- RunUMAP(myeloid, dims = 1:20)
myeloid <- FindNeighbors(myeloid, dims = 1:20)
myeloid <- FindClusters(myeloid, resolution = 0.3)


identity_mapping <- c(
  "0" = "DC_cDC2_CD1C",
  "1" = "Inflam-TAM_VCAN",
  "2" = "Lipid-TAM_MARCO",
  "3" = "B",
  "4" = "Repair-TAM_SLC40A1",
  "5" = "IFN-TAM_CXCL9",
  "6" = "Inflam-TAM_CXCL8",
  "7" = "MT",
  "8" = "Prolif-TAM_MKI67",
  "9" = "DC_pDC_LILRA4",
  "10" = "DC_cDC1_CLEC9A",
  "11" = "Lipid-TAM_APOE",
  "12" = "Lipid-TAM_TREM2",  
  "13" = "DC_LC_CD207",
  "14" = "DC_mDC_LAMP3",
  "15" = "Lipid-TAM_SPP1"
  )

cell_type <- identity_mapping[myeloid@meta.data$seurat_clusters]
myeloid@meta.data$cell_type <- cell_type

identity_mapping <- c(
  "0" = "DC",
  "1" = "Macro",
  "2" = "Macro",
  "3" = "B",
  "4" = "Macro",
  "5" = "Macro",
  "6" = "Macro",
  "7" = "MT",
  "8" = "Macro",
  "9" = "DC",
  "10" = "DC", 
  "11" = "Macro",
  "12" = "Macro",
  "13" = "DC",
  "14" = "DC",
  "15" = "Macro"
)

major_cell <- identity_mapping[myeloid@meta.data$seurat_clusters]
myeloid@meta.data$major_cell <- major_cell

saveRDS(myeloid,file="myeloid_anno.rds")


```

# response
```
LUSC <- subset(myeloid,subset=cancer_type=="LUSC")
LUAD <- subset(myeloid, subset = cancer_type == "LUAD")
macro <- subset(LUSC,subset=major_cell=="Macro")
DC <- subset(LUSC,subset=major_cell=="DC")

macro <- subset(LUAD,subset=major_cell=="Macro")
DC <- subset(LUAD,subset=major_cell=="DC")
```


# 分布图(macro+DC)
```
library(Seurat)
library(ggplot2)
library(ggsci)
library(dplyr)

# 提取LUAD和LUSC数据
LUAD <- subset(myeloid, subset = cancer_type == "LUAD")
LUSC <- subset(myeloid, subset = cancer_type == "LUSC")

# 定义绘图函数
create_ln_group_umap <- function(seurat_obj, cancer_type, cell_type_name) {
  # 根据细胞类型筛选
  if (cell_type_name == "Macro") {
    plot_data <- subset(seurat_obj, subset = major_cell == "Macro")
  } else if (cell_type_name == "DC") {
    plot_data <- subset(seurat_obj, subset = major_cell == "DC")
  } else {
    stop("cell_type_name 必须是 'Macro' 或 'DC'")
  }
  
  # 检查数据
  if (ncol(plot_data) == 0) {
    cat(cancer_type, "中", cell_type_name, "没有数据，跳过\n")
    return(NULL)
  }
  
  # 设置LN_group因子顺序
  plot_data$LN_group <- factor(
    plot_data$LN_group,
    levels = c("LN+", "LN-")
  )
  
  # 获取颜色
  npg_pal <- pal_npg()(10)
  npg_extended <- colorRampPalette(npg_pal)(16)
  
  # 创建图形
  p <- DimPlot(plot_data, 
               reduction = "umap", 
               label = TRUE, 
               pt.size = 1, 
               group.by = "cell_type", 
               label.size = 4,
               split.by = "LN_group",  
               ncol = 2) + 
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    ggtitle(paste0(cancer_type, " ", cell_type_name)) +
    scale_color_manual(values = npg_extended) +
    coord_fixed(ratio = 1) +
    guides(color = guide_legend(title = NULL, override.aes = list(size = 5))) +
    theme(
      text = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(size = 20, color = "black"),
      axis.text.y = element_text(size = 20, color = "black"),
      axis.title.x = element_text(size = 24, face = "bold", color = "black"),
      axis.title.y = element_text(size = 24, face = "bold", color = "black", margin = margin(r = 20)),
      plot.title = element_text(size = 28, face = "bold", hjust = 0.5, margin = margin(b = 20)),
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
      plot.margin = margin(20, 50, 10, 10),
      strip.text = element_text(size = 22, face = "bold", 
                               margin = margin(b = 15)),
      panel.spacing = unit(1.5, "lines")
    )
  
  return(list(plot = p, cell_count = ncol(plot_data)))
}



# 定义要绘制的组合
combinations <- expand.grid(
  cancer_type = c("LUAD", "LUSC"),
  cell_type = c("Macro", "DC"),
  stringsAsFactors = FALSE
)

# 循环生成所有图形
for (i in 1:nrow(combinations)) {
  cancer_type <- combinations$cancer_type[i]
  cell_type <- combinations$cell_type[i]
  
  cat("\n处理:", cancer_type, "-", cell_type, "\n")
  
  # 获取相应的Seurat对象
  if (cancer_type == "LUAD") {
    seurat_obj <- LUAD
  } else {
    seurat_obj <- LUSC
  }
  
  # 创建图形
  result <- create_ln_group_umap(seurat_obj, cancer_type, cell_type)
  
  if (is.null(result)) {
    next
  }
  
  p <- result$plot
  cell_count <- result$cell_count
  
  cat("  细胞数量:", cell_count, "\n")
  
  # 1. 保存为PDF
  pdf_file <- paste0(cancer_type, "_", cell_type, "_annotation_by_group.pdf")
  pdf(pdf_file, width = 4000/300, height = 2000/300)
  print(p)
  dev.off()
  cat("  ✓ PDF: ", pdf_file, "\n")

  highres_tiff_file <- paste0(cancer_type, "_", cell_type, "_annotation_by_group.tiff")
  ggsave(highres_tiff_file,
         plot = p,
         width = 4000/300,
         height = 2000/300,
         dpi = 600,  # 更高分辨率
         device = "tiff",
         compression = "lzw",
         bg = "white")
  cat("  ✓ 高分辨率TIFF: ", highres_tiff_file, "\n")
}


```

# 箱图(LUSC)
```
library(ggpubr)
library(dplyr)
library(ggsci)
library(tidyr)
LUSC <- subset(myeloid,subset=cancer_type=="LUSC")
LUAD <- subset(myeloid,subset=cancer_type=="LUAD")
macro <- subset(LUSC,subset=major_cell=="Macro")
DC <- subset(LUSC,subset=major_cell=="DC")

# 统计每个样本-细胞类型-LN分组的细胞数
cell_counts <- macro@meta.data %>%
  group_by(patients, cell_type, LN_group) %>%  # Group改为LN_group
  summarise(count = n(), .groups = "drop")

# 计算每个样本的总细胞数
sample_totals <- cell_counts %>%
  group_by(patients) %>%
  summarise(total = sum(count), .groups = "drop")

# 计算细胞类型比例
cell_proportions <- cell_counts %>%
  left_join(sample_totals, by = "patients") %>%
  mutate(proportion = count / total)

# === 找出哪个细胞类型在LN-中缺失 ===
# 检查每个细胞类型在两个LN_group中的存在情况
celltype_group_presence <- cell_proportions %>%
  group_by(cell_type, LN_group) %>%  # Group改为LN_group
  summarise(has_data = n() > 0, .groups = "drop") %>%
  pivot_wider(names_from = LN_group, values_from = has_data, values_fill = FALSE)

# 找出在LN+中存在但在LN-中缺失的细胞类型
missing_in_ln_neg <- celltype_group_presence %>%
  filter(`LN+` == TRUE & `LN-` == FALSE) %>%
  pull(cell_type)

cat("在LN-中缺失的细胞类型:", paste(missing_in_ln_neg, collapse = ", "), "\n")

# === 只为缺失的细胞类型在LN-中添加0值 ===
if (length(missing_in_ln_neg) > 0) {
  # 获取LN-的所有患者
  ln_neg_patients <- unique(cell_proportions$patients[cell_proportions$LN_group == "LN-"])
  
  # 只为缺失的细胞类型创建0值记录
  missing_records <- expand_grid(
    patients = ln_neg_patients,
    cell_type = missing_in_ln_neg,
    LN_group = "LN-"
  ) %>%
  left_join(sample_totals, by = "patients") %>%
  mutate(count = 0, proportion = 0)
  
  # 合并到绘图数据
  plot_data <- cell_proportions %>% bind_rows(missing_records)
} else {
  plot_data <- cell_proportions
}

# === 关键修复：确保分组数量匹配 ===
# 获取所有细胞类型和LN_group
all_celltypes <- unique(plot_data$cell_type)
all_ln_groups <- c("LN+", "LN-")  # 明确指定两个分组
n_groups <- length(all_ln_groups)

# 构建完整的细胞类型-LN_group组合并填充0
full_combinations <- expand.grid(
  cell_type = all_celltypes,
  LN_group = all_ln_groups,
  stringsAsFactors = FALSE
)

# 按细胞类型和LN_group汇总计数
cell_group_counts <- cell_counts %>%
  group_by(cell_type, LN_group) %>%  # Group改为LN_group
  summarise(total_count = sum(count), .groups = "drop") %>%
  right_join(full_combinations, by = c("cell_type", "LN_group")) %>%  # Group改为LN_group
  mutate(total_count = replace(total_count, is.na(total_count), 0)) %>%
  arrange(cell_type, LN_group)  # 确保分组顺序一致

# 计算每个LN_group的总细胞数（全局）
global_group_totals <- cell_counts %>%
  group_by(LN_group) %>%  # Group改为LN_group
  summarise(group_total = sum(count), .groups = "drop") %>%
  right_join(data.frame(LN_group = all_ln_groups), by = "LN_group") %>%  # Group改为LN_group
  mutate(group_total = replace(group_total, is.na(group_total), 0)) %>%
  arrange(LN_group)  # 保持分组顺序一致

# === 修复卡方检验逻辑 ===
chisq_results <- cell_group_counts %>%
  group_by(cell_type) %>%
  summarise(
    p_value = {
      # 获取当前细胞类型在两个LN_group的计数
      obs <- c(
        total_count[LN_group == "LN+"],
        total_count[LN_group == "LN-"]
      )
      
      # 获取对应LN_group的总细胞数
      group_tot <- c(
        global_group_totals$group_total[global_group_totals$LN_group == "LN+"],
        global_group_totals$group_total[global_group_totals$LN_group == "LN-"]
      )
      
      # 检查数据有效性
      if (sum(obs) < 5 || length(obs) != 2) {
        NA_real_
      } else {
        # 构建2x2列联表（细胞类型计数 vs 其他细胞计数）
        other_counts <- group_tot - obs
        cont_table <- matrix(c(obs[1], other_counts[1],
                               obs[2], other_counts[2]),
                             nrow = 2, byrow = TRUE)
        
        # 检查列联表有效性
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
  # 添加显著性标记
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

# 设置LN_group的顺序
plot_data$LN_group <- factor(plot_data$LN_group, levels = c("LN+", "LN-"))

# === 绘制图形 ===
y_max <- max(plot_data$proportion, na.rm = TRUE) * 1.1

pdf("LUSC_macro_box_proportion_by_ln_group.pdf", width = 6000/300, height = 3000/300)
ggplot(plot_data, aes(x = cell_type, y = proportion, fill = LN_group)) +  # Group改为LN_group
  geom_boxplot(width = 0.7) +
  geom_text(
    data = chisq_results, 
    aes(x = cell_type, y = y_max, label = significance),
    size = 12, 
    inherit.aes = FALSE
  ) +
  labs(x = "", y = "Cell Proportion", fill = "LN Group") +  # Group改为LN Group
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

# 首先创建图形对象
p <- ggplot(plot_data, aes(x = cell_type, y = proportion, fill = LN_group)) +
  geom_boxplot(width = 0.7) +
  geom_text(
    data = chisq_results, 
    aes(x = cell_type, y = y_max, label = significance),
    size = 12, 
    inherit.aes = FALSE
  ) +
  labs(x = "", y = "Cell Proportion", fill = "LN Group") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.y = element_text(size = 32),
    axis.line = element_line(color = "black"),
    panel.grid = element_blank(),
    legend.text = element_text(size = 36),
    legend.title = element_text(size = 40),
    legend.position = "right",
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  scale_fill_npg() +
  scale_y_continuous(limits = c(0, y_max))

# 使用ggsave保存为TIFF
ggsave("LUSC_macro_box_proportion_by_ln_group.tiff",
       plot = p,
       width = 6000/300,  # 20英寸
       height = 3000/300, # 10英寸
       units = "in",
       dpi = 300,
       device = "tiff",
       compression = "lzw",
       bg = "white")

```
# 饼图
```
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(16)
library(dplyr)

LUAD <- subset(myeloid, subset = cancer_type == "LUAD")
LUSC <- subset(myeloid, subset = cancer_type == "LUSC")


# 先提取Macrophage相关的数据
macro_data <- LUSC@meta.data %>%
  filter(major_cell == "Macro") %>%  # 只保留Macrophage
  group_by(LN_group, major_cell, cell_type) %>%
  summarise(count = n()) %>%
  group_by(LN_group, major_cell) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup() %>%
  mutate(group_major = paste0(LN_group, "_", major_cell)) %>%
  mutate(group_major = factor(group_major,
                             levels = c("LN+_Macro", "LN-_Macro")))

# 画Macrophage的饼图
pdf("LUSC_Group_Macro_cell_type_pie.pdf", width = 4000/300, height = 3000/300) 

ggplot(macro_data, aes(x = "", y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack", width = 1) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = npg_extended) +
  facet_wrap(~ group_major, ncol = 2) +  # 改为2列
  labs(x = "", y = "", fill = "Cell Type", 
       title = "Macrophage Cell Type Composition") +
  theme_void() +
  theme(
    legend.text = element_text(size = 28),
    legend.title = element_text(size = 28),
    strip.text = element_text(
      size = 20,
      face = "bold",
      margin = margin(b = 15)
    ),
    legend.position = "right",
    plot.title = element_text(size = 32, face = "bold", hjust = 0.5, margin = margin(b = 20))
  )

dev.off()
#tiff
ggsave(
  filename = "LUSC_Group_Macro_cell_type_pie.tiff",
  plot = last_plot(),
  device = "tiff",
  width = 4000/300,  # 单位：英寸
  height = 3000/300,
  units = "in",
  dpi = 300,
  compression = "lzw",
  bg = "white"
)




#LUAD
macro_data <- LUAD@meta.data %>%
  filter(major_cell == "Macro") %>%  # 只保留Macrophage
  group_by(LN_group, major_cell, cell_type) %>%
  summarise(count = n()) %>%
  group_by(LN_group, major_cell) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup() %>%
  mutate(group_major = paste0(LN_group, "_", major_cell)) %>%
  mutate(group_major = factor(group_major,
                             levels = c("LN+_Macro", "LN-_Macro")))

# 画Macrophage的饼图
pdf("LUAD_Group_Macro_cell_type_pie.pdf", width = 4000/300, height = 3000/300) 

ggplot(macro_data, aes(x = "", y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack", width = 1) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = npg_extended) +
  facet_wrap(~ group_major, ncol = 2) +  # 改为2列
  labs(x = "", y = "", fill = "Cell Type", 
       title = "Macrophage Cell Type Composition") +
  theme_void() +
  theme(
    legend.text = element_text(size = 28),
    legend.title = element_text(size = 28),
    strip.text = element_text(
      size = 20,
      face = "bold",
      margin = margin(b = 15)
    ),
    legend.position = "right",
    plot.title = element_text(size = 32, face = "bold", hjust = 0.5, margin = margin(b = 20))
  )

dev.off()


ggsave(
  filename = "LUAD_Group_Macro_cell_type_pie.tiff",
  plot = last_plot(),
  device = "tiff",
  width = 4000/300,  # 单位：英寸
  height = 3000/300,
  units = "in",
  dpi = 300,
  compression = "lzw",
  bg = "white"
)


```

# macro二合一(LUSC)
```
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(patchwork)

# ===== 第一部分：准备柱状图数据 =====
prepare_bar_plot_data <- function(macro_obj) {
  cell_proportions_bar <- macro_obj@meta.data %>%
    group_by(LN_group, pathological_response, cell_type) %>%
    summarise(count = n(), .groups = 'drop') %>%
    group_by(LN_group, pathological_response) %>%
    mutate(proportion = count / sum(count)) %>%
    ungroup()
  
  # 设置LN组顺序：LN+在前，LN-在后
  cell_proportions_bar$LN_group <- factor(
    cell_proportions_bar$LN_group,
    levels = c("LN+", "LN-")
  )
  
  # 设置病理反应顺序
  cell_proportions_bar$pathological_response <- factor(
    cell_proportions_bar$pathological_response, 
    levels = c("pCR", "MPR", "non-MPR")
  )
  
  all_cell_types <- unique(macro_obj@meta.data$cell_type)
  cell_proportions_complete <- cell_proportions_bar %>%
    complete(LN_group, pathological_response, cell_type = all_cell_types, 
             fill = list(count = 0, proportion = 0))
  
  # 保持LN组和病理反应顺序
  cell_proportions_complete$LN_group <- factor(
    cell_proportions_complete$LN_group,
    levels = c("LN+", "LN-")
  )
  
  cell_proportions_complete$pathological_response <- factor(
    cell_proportions_complete$pathological_response, 
    levels = c("pCR", "MPR", "non-MPR")
  )
  
  celltype_alphabetical <- sort(all_cell_types)
  cell_proportions_complete$cell_type <- factor(
    cell_proportions_complete$cell_type, 
    levels = celltype_alphabetical
  )
  
  return(list(
    bar_data = cell_proportions_complete,
    celltypes = celltype_alphabetical
  ))
}

# ===== 第二部分：数据处理模块 =====

# 模块1: 提取和计算细胞计数
extract_cell_counts <- function(macro_obj, ln_group) {
  cell_counts <- macro_obj@meta.data %>%
    filter(LN_group == ln_group) %>%
    group_by(patients, cell_type, pathological_response) %>%
    summarise(count = n(), .groups = "drop")
  
  return(cell_counts)
}

# 模块2: 计算样本总数
calculate_sample_totals <- function(cell_counts) {
  sample_totals <- cell_counts %>%
    group_by(patients) %>%
    summarise(total = sum(count), .groups = "drop")
  
  return(sample_totals)
}

# 模块3: 计算细胞比例
calculate_cell_proportions <- function(cell_counts, sample_totals) {
  cell_proportions <- cell_counts %>%
    left_join(sample_totals, by = "patients") %>%
    mutate(proportion = count / total)
  
  # 设置病理反应顺序
  cell_proportions$pathological_response <- factor(
    cell_proportions$pathological_response, 
    levels = c("pCR", "MPR", "non-MPR")
  )
  
  return(cell_proportions)
}

# 模块4: 补充缺失数据
supplement_missing_data <- function(cell_proportions, sample_totals) {
  # 找出缺失的病理反应并添加0值
  celltype_response_presence <- cell_proportions %>%
    group_by(cell_type, pathological_response) %>%
    summarise(has_data = n() > 0, .groups = "drop") %>%
    pivot_wider(names_from = pathological_response, values_from = has_data, values_fill = FALSE)
  
  # 检查缺失的组合
  all_responses <- c("pCR", "MPR", "non-MPR")
  missing_combinations <- list()
  
  for (resp in all_responses) {
    missing_celltypes <- celltype_response_presence %>%
      filter(.data[[resp]] == FALSE) %>%
      pull(cell_type)
    
    if (length(missing_celltypes) > 0) {
      missing_combinations[[resp]] <- missing_celltypes
    }
  }
  
  plot_data <- cell_proportions
  
  # 为缺失的组合添加0值
  for (resp in names(missing_combinations)) {
    resp_patients <- unique(cell_proportions$patients[cell_proportions$pathological_response == resp])
    
    if (length(resp_patients) > 0) {
      missing_records <- expand_grid(
        patients = resp_patients,
        cell_type = missing_combinations[[resp]],
        pathological_response = resp
      ) %>%
      left_join(sample_totals, by = "patients") %>%
      mutate(count = 0, proportion = 0)
      
      plot_data <- plot_data %>% bind_rows(missing_records)
    }
  }
  
  return(plot_data)
}

# 模块5: 执行统计检验
perform_statistical_tests <- function(cell_counts, plot_data) {
  # 统计检验 - 三组比较
  all_celltypes <- unique(plot_data$cell_type)
  full_combinations <- expand_grid(
    cell_type = all_celltypes,
    pathological_response = factor(c("pCR", "MPR", "non-MPR"), 
                                  levels = c("pCR", "MPR", "non-MPR"))
  )
  
  cell_response_counts <- cell_counts %>%
    group_by(cell_type, pathological_response) %>%
    summarise(total_count = sum(count), .groups = "drop") %>%
    right_join(full_combinations, by = c("cell_type", "pathological_response")) %>%
    mutate(total_count = replace(total_count, is.na(total_count), 0))
  
  # 计算全局病理反应组总数
  global_response_totals <- cell_counts %>%
    group_by(pathological_response) %>%
    summarise(response_total = sum(count), .groups = "drop") %>%
    right_join(data.frame(pathological_response = factor(c("pCR", "MPR", "non-MPR"), 
                                                        levels = c("pCR", "MPR", "non-MPR"))), 
              by = "pathological_response") %>%
    mutate(response_total = replace(response_total, is.na(response_total), 0))
  
  # 卡方检验（三组比较）
  chisq_results <- cell_response_counts %>%
    group_by(cell_type) %>%
    summarise(
      p_value = calculate_p_value_for_celltype(cur_data(), global_response_totals),
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
  
  return(chisq_results)
}

# 模块6: 为单个细胞类型计算p值
calculate_p_value_for_celltype <- function(cell_data, global_response_totals) {
  # 获取三组的计数
  pcr_count <- cell_data$total_count[cell_data$pathological_response == "pCR"]
  mpr_count <- cell_data$total_count[cell_data$pathological_response == "MPR"]
  non_mpr_count <- cell_data$total_count[cell_data$pathological_response == "non-MPR"]
  
  pcr_total <- global_response_totals$response_total[global_response_totals$pathological_response == "pCR"]
  mpr_total <- global_response_totals$response_total[global_response_totals$pathological_response == "MPR"]
  non_mpr_total <- global_response_totals$response_total[global_response_totals$pathological_response == "non-MPR"]
  
  total_counts <- sum(c(pcr_count, mpr_count, non_mpr_count))
  
  if (total_counts < 15 || any(c(pcr_count, mpr_count, non_mpr_count) < 5)) {
    return(NA_real_)
  }
  
  # 构建3x2列联表
  other_pcr <- pcr_total - pcr_count
  other_mpr <- mpr_total - mpr_count
  other_non_mpr <- non_mpr_total - non_mpr_count
  
  cont_table <- matrix(c(pcr_count, other_pcr,
                        mpr_count, other_mpr,
                        non_mpr_count, other_non_mpr),
                      nrow = 3, ncol = 2, byrow = TRUE)
  
  if (any(cont_table < 0) || sum(cont_table) == 0) {
    return(NA_real_)
  } else if (any(cont_table < 5)) {
    return(fisher.test(cont_table)$p.value)
  } else {
    return(chisq.test(cont_table)$p.value)
  }
}

# 模块7: 处理单个LN组数据
process_ln_group_data <- function(macro_obj, ln_group) {
  # 提取细胞计数
  cell_counts <- extract_cell_counts(macro_obj, ln_group)
  
  # 计算样本总数
  sample_totals <- calculate_sample_totals(cell_counts)
  
  # 计算细胞比例
  cell_proportions <- calculate_cell_proportions(cell_counts, sample_totals)
  
  # 补充缺失数据
  plot_data <- supplement_missing_data(cell_proportions, sample_totals)
  
  # 执行统计检验
  chisq_results <- perform_statistical_tests(cell_counts, plot_data)
  
  return(list(
    plot_data = plot_data,
    chisq_results = chisq_results
  ))
}

# ===== 第三部分：主数据处理流程（确保LN+在前） =====
process_box_plot_data <- function(macro_obj, celltype_alphabetical) {
  # 【关键修改】按顺序处理：先LN+，后LN-
  ln_groups_ordered <- c("LN+", "LN-")
  
  # 初始化列表存储结果
  plot_data_list <- list()
  stats_list <- list()
  
  for (ln_group in ln_groups_ordered) {
    group_data <- process_ln_group_data(macro_obj, ln_group)
    
    plot_data_list[[ln_group]] <- group_data$plot_data %>% mutate(LN_group = ln_group)
    stats_list[[ln_group]] <- group_data$chisq_results %>% mutate(LN_group = ln_group)
  }
  
  # 合并数据，保持LN+在前
  plot_data_box_combined <- bind_rows(plot_data_list)
  chisq_results_combined <- bind_rows(stats_list)
  
  # 设置LN组因子顺序：LN+在前
  plot_data_box_combined$LN_group <- factor(
    plot_data_box_combined$LN_group,
    levels = c("LN+", "LN-")
  )
  
  chisq_results_combined$LN_group <- factor(
    chisq_results_combined$LN_group,
    levels = c("LN+", "LN-")
  )
  
  # 确保数据按病理反应顺序排序
  plot_data_box_combined <- plot_data_box_combined %>%
    arrange(LN_group, cell_type, pathological_response) %>%
    mutate(
      pathological_response = factor(pathological_response, 
                                    levels = c("pCR", "MPR", "non-MPR")),
      order_factor = interaction(cell_type, pathological_response, sep = "_")
    )
  
  # 设置细胞类型顺序
  plot_data_box_combined$cell_type <- factor(plot_data_box_combined$cell_type, 
                                            levels = celltype_alphabetical)
  chisq_results_combined$cell_type <- factor(chisq_results_combined$cell_type, 
                                            levels = celltype_alphabetical)
  
  # 计算y轴最大值
  y_max <- max(plot_data_box_combined$proportion, na.rm = TRUE) * 1.1
  
  return(list(
    box_data = plot_data_box_combined,
    stats_results = chisq_results_combined,
    y_max = y_max
  ))
}

# ===== 第四部分：绘图模块 =====

# 模块8: 创建柱状图
create_bar_plot <- function(bar_data, celltypes) {
  # 设置颜色
  response_colors <- c("pCR" = "#00A087FF", "MPR" = "#F39B7FFF", "non-MPR" = "#E64B35FF")
  
  p <- ggplot(bar_data, aes(x = cell_type, y = proportion, fill = pathological_response)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    scale_fill_manual(
      values = response_colors,
      breaks = c("pCR", "MPR", "non-MPR"),
      labels = c("pCR", "MPR", "non-MPR")
    ) +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(x = "", y = "Proportion", title = "Bar Plot") +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black", face = "bold"),
      axis.text.y = element_text(size = 12, color = "black", face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold", color = "black"),
      axis.line = element_line(color = "black"),
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(size = 14, face = "bold"),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
    ) +
    facet_grid(LN_group ~ ., scales = "free_x")
  
  return(p)
}

# 模块9: 创建箱图
create_box_plot <- function(box_data, stats_results, y_max, celltypes) {
  # 设置颜色
  response_colors <- c("pCR" = "#00A087FF", "MPR" = "#F39B7FFF", "non-MPR" = "#E64B35FF")
  
  p <- ggplot(box_data, aes(x = cell_type, y = proportion, fill = pathological_response)) +
    geom_boxplot(width = 0.7, aes(group = interaction(cell_type, pathological_response))) +
    # 为0值数据添加横线
    geom_segment(
      data = box_data %>% filter(proportion == 0),
      aes(x = as.numeric(cell_type) - 0.25, 
          xend = as.numeric(cell_type) + 0.25,
          y = 0, yend = 0),
      color = "black", linewidth = 1, inherit.aes = FALSE
    ) +
    # 添加显著性标记
    geom_text(
      data = stats_results, 
      aes(x = cell_type, y = y_max, label = significance),
      size = 6, 
      inherit.aes = FALSE
    ) +
    scale_fill_manual(
      values = response_colors,
      breaks = c("pCR", "MPR", "non-MPR"),
      labels = c("pCR", "MPR", "non-MPR")
    ) +
    scale_y_continuous(labels = scales::percent_format(), limits = c(0, y_max)) +
    labs(x = "Cell Type", y = "Cell Proportion", title = "Box Plot") +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black", face = "bold"),
      axis.text.y = element_text(size = 12, color = "black", face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold", color = "black"),
      axis.line = element_line(color = "black"),
      panel.grid = element_blank(),
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(size = 14, face = "bold"),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
    ) +
    facet_grid(LN_group ~ ., scales = "free_x")
  
  return(p)
}

# 模块10: 创建图例
create_legend <- function(bar_data) {
  response_colors <- c("pCR" = "#00A087FF", "MPR" = "#F39B7FFF", "non-MPR" = "#E64B35FF")
  
  legend_plot <- ggplot(bar_data, 
                        aes(x = cell_type, y = proportion, fill = pathological_response)) +
    geom_col(position = position_dodge(width = 0.8)) +
    scale_fill_manual(
      values = response_colors,
      breaks = c("pCR", "MPR", "non-MPR"),
      labels = c("pCR", "MPR", "non-MPR"),
      name = "Pathological Response"
    ) +
    theme(legend.position = "top",
          legend.text = element_text(size = 14, face = "bold"),
          legend.title = element_text(size = 16, face = "bold"))
  
  legend <- get_legend(legend_plot)
  return(legend)
}

# ===== 第五部分：主程序 =====

# 1. 准备柱状图数据
bar_data_result <- prepare_bar_plot_data(macro)
bar_data <- bar_data_result$bar_data
celltypes <- bar_data_result$celltypes

# 2. 处理箱图数据（确保LN+在前）
box_data_result <- process_box_plot_data(macro, celltypes)
box_data <- box_data_result$box_data
stats_results <- box_data_result$stats_results
y_max <- box_data_result$y_max

# 3. 创建图形
bar_plot <- create_bar_plot(bar_data, celltypes)
box_plot <- create_box_plot(box_data, stats_results, y_max, celltypes)
legend <- create_legend(bar_data)

# 4. 组合图形
combined_plot <- (bar_plot | box_plot) +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")

# 添加图例
final_plot <- wrap_plots(legend, combined_plot, ncol = 1, heights = c(0.05, 0.95))

# 5. 保存图形
pdf("LUSC_macro_combined_bar_box_plots.pdf", width = 6000/300, height = 3000/300)
print(final_plot)
dev.off()

# 同时保存TIFF格式
ggsave("LUSC_macro_combined_bar_box_plots.tiff",
       plot = final_plot,
       width = 6000/300,  # 20英寸 (6000/300 = 20)
       height = 3000/300, # 10英寸 (3000/300 = 10)
       units = "in",
       dpi = 300,
       device = "tiff",
       compression = "lzw",
       bg = "white")

cat("完成！已生成PDF和TIFF格式的图形。\n")
cat("图形顺序：LN+在上方，LN-在下方\n")
```
# 二合一（LUAD）
```
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(patchwork)

# ===== 第一部分：准备柱状图数据 =====
prepare_bar_plot_data <- function(macro_obj) {
  cell_proportions_bar <- macro_obj@meta.data %>%
    group_by(LN_group, pathological_response, cell_type) %>%
    summarise(count = n(), .groups = 'drop') %>%
    group_by(LN_group, pathological_response) %>%
    mutate(proportion = count / sum(count)) %>%
    ungroup()
  
  # 设置LN组顺序：LN+在前，LN-在后
  cell_proportions_bar$LN_group <- factor(
    cell_proportions_bar$LN_group,
    levels = c("LN+", "LN-")
  )
  
  # 设置病理反应顺序
  cell_proportions_bar$pathological_response <- factor(
    cell_proportions_bar$pathological_response, 
    levels = c("pCR", "MPR", "non-MPR")
  )
  
  all_cell_types <- unique(macro_obj@meta.data$cell_type)
  
  # 补全数据，包括缺失的病理反应
  cell_proportions_complete <- cell_proportions_bar %>%
    complete(LN_group, pathological_response, cell_type = all_cell_types, 
             fill = list(count = 0, proportion = 0))
  
  # 保持LN组和病理反应顺序
  cell_proportions_complete$LN_group <- factor(
    cell_proportions_complete$LN_group,
    levels = c("LN+", "LN-")
  )
  
  cell_proportions_complete$pathological_response <- factor(
    cell_proportions_complete$pathological_response, 
    levels = c("pCR", "MPR", "non-MPR")
  )
  
  celltype_alphabetical <- sort(all_cell_types)
  cell_proportions_complete$cell_type <- factor(
    cell_proportions_complete$cell_type, 
    levels = celltype_alphabetical
  )
  
  return(list(
    bar_data = cell_proportions_complete,
    celltypes = celltype_alphabetical
  ))
}

# ===== 第二部分：数据处理模块 =====

# 模块1: 提取和计算细胞计数
extract_cell_counts <- function(macro_obj, ln_group) {
  cell_counts <- macro_obj@meta.data %>%
    filter(LN_group == ln_group) %>%
    group_by(patients, cell_type, pathological_response) %>%
    summarise(count = n(), .groups = "drop")
  
  return(cell_counts)
}

# 模块2: 计算样本总数
calculate_sample_totals <- function(cell_counts) {
  sample_totals <- cell_counts %>%
    group_by(patients) %>%
    summarise(total = sum(count), .groups = "drop")
  
  return(sample_totals)
}

# 模块3: 计算细胞比例 - 不添加缺失的pCR数据
calculate_cell_proportions <- function(cell_counts, sample_totals) {
  cell_proportions <- cell_counts %>%
    left_join(sample_totals, by = "patients") %>%
    mutate(proportion = count / total)
  
  # 设置病理反应顺序
  cell_proportions$pathological_response <- factor(
    cell_proportions$pathological_response, 
    levels = c("pCR", "MPR", "non-MPR")
  )
  
  return(cell_proportions)
}

# 模块4-1: 检查数据是否足够进行统计检验
check_data_sufficiency <- function(ct_data, cell_type) {
  # 找出当前细胞类型实际存在的病理反应组
  ct_responses <- unique(ct_data$pathological_response)
  ct_responses <- as.character(ct_responses[!is.na(ct_responses)])
  
  # 计算每个组有数据的样本数
  group_counts <- sapply(ct_responses, function(resp) {
    sum(ct_data$proportion[ct_data$pathological_response == resp] > 0, na.rm = TRUE)
  })
  
  # 计算总的有数据的样本数
  total_valid_samples <- sum(ct_data$proportion > 0, na.rm = TRUE)
  
  list(
    valid_groups = ct_responses,
    group_counts = group_counts,
    total_samples = total_valid_samples,
    is_sufficient = length(ct_responses) >= 2 && total_valid_samples >= 5
  )
}

# 模块4-2: 准备列联表数据
prepare_contingency_table <- function(cell_counts, cell_type, valid_groups) {
  # 获取当前细胞类型的计数
  ct_counts <- cell_counts %>%
    filter(cell_type == !!cell_type) %>%
    group_by(pathological_response) %>%
    summarise(total_count = sum(count), .groups = "drop")
  
  # 获取该LN组的总细胞数
  group_totals <- cell_counts %>%
    group_by(pathological_response) %>%
    summarise(group_total = sum(count), .groups = "drop")
  
  # 合并数据，确保所有反应类型都有值
  all_responses <- c("pCR", "MPR", "non-MPR")
  ct_counts_complete <- data.frame(pathological_response = all_responses) %>%
    left_join(ct_counts, by = "pathological_response") %>%
    left_join(group_totals, by = "pathological_response") %>%
    mutate(
      total_count = ifelse(is.na(total_count), 0, total_count),
      group_total = ifelse(is.na(group_total), 0, group_total),
      other_count = group_total - total_count
    )
  
  # 筛选实际存在的组别
  existing_data <- ct_counts_complete %>%
    filter(pathological_response %in% valid_groups)
  
  # 构建列联表
  cont_table <- matrix(
    c(existing_data$total_count, existing_data$other_count),
    nrow = length(valid_groups),
    ncol = 2,
    byrow = FALSE
  )
  
  return(cont_table)
}

# 模块4-3: 执行统计检验
perform_statistical_test <- function(cont_table, valid_groups, cell_type, ln_group) {
  p_value <- NA_real_
  
  # 检查数据有效性
  if (sum(cont_table) > 0 && all(cont_table >= 0)) {
    # 根据组数选择检验方法
    if (length(valid_groups) == 2) {
      # 2×2列联表
      if (any(cont_table < 5)) {
        # 使用Fisher精确检验
        tryCatch({
          p_value <- fisher.test(cont_table)$p.value
          cat(ln_group, "-", cell_type, ": 使用Fisher精确检验\n")
        }, error = function(e) {
          cat(ln_group, "-", cell_type, ": Fisher检验出错:", e$message, "\n")
        })
      } else {
        # 使用卡方检验
        tryCatch({
          p_value <- chisq.test(cont_table)$p.value
          cat(ln_group, "-", cell_type, ": 使用卡方检验\n")
        }, error = function(e) {
          cat(ln_group, "-", cell_type, ": 卡方检验出错:", e$message, "\n")
        })
      }
    } else if (length(valid_groups) == 3) {
      # 3×2列联表
      if (any(cont_table < 5)) {
        # 使用Fisher精确检验
        tryCatch({
          p_value <- fisher.test(cont_table)$p.value
          cat(ln_group, "-", cell_type, ": 使用Fisher精确检验\n")
        }, error = function(e) {
          cat(ln_group, "-", cell_type, ": Fisher检验出错:", e$message, "\n")
        })
      } else {
        # 使用卡方检验
        tryCatch({
          p_value <- chisq.test(cont_table)$p.value
          cat(ln_group, "-", cell_type, ": 使用卡方检验\n")
        }, error = function(e) {
          cat(ln_group, "-", cell_type, ": 卡方检验出错:", e$message, "\n")
        })
      }
    }
  }
  
  return(p_value)
}

# 模块4-4: 主统计检验函数
perform_statistical_tests_main <- function(cell_counts, plot_data, ln_group) {
  cat("\n=== 处理", ln_group, "组 ===\n")
  
  # 对每个细胞类型执行统计检验
  all_celltypes <- unique(plot_data$cell_type)
  
  chisq_results <- data.frame()
  
  for (ct in all_celltypes) {
    # 获取当前细胞类型的数据
    ct_data <- plot_data %>% filter(cell_type == ct)
    
    # 检查数据是否足够
    data_check <- check_data_sufficiency(ct_data, ct)
    
    p_value <- NA_real_
    
    if (data_check$is_sufficient) {
      cat(ln_group, "-", ct, ": ", length(data_check$valid_groups), 
          "个有数据组 (", paste(data_check$valid_groups, collapse = ","), 
          "), ", data_check$total_samples, "个有数据样本\n")
      
      # 准备列联表
      cont_table <- prepare_contingency_table(cell_counts, ct, data_check$valid_groups)
      
      # 执行统计检验
      p_value <- perform_statistical_test(cont_table, data_check$valid_groups, ct, ln_group)
    } else {
      cat(ln_group, "-", ct, ": 数据不足，跳过统计检验（",
          length(data_check$valid_groups), "个有数据组, ",
          data_check$total_samples, "个有数据样本）\n")
    }
    
    # 添加结果
    result_row <- data.frame(
      cell_type = ct,
      p_value = p_value,
      significance = ifelse(is.na(p_value) || p_value >= 0.05, "ns",
                           ifelse(p_value < 0.001, "***",
                                  ifelse(p_value < 0.01, "**", "*")))
    )
    
    chisq_results <- bind_rows(chisq_results, result_row)
  }
  
  return(chisq_results)
}

# 模块5: 处理单个LN组数据
process_ln_group_data <- function(macro_obj, ln_group) {
  # 提取细胞计数
  cell_counts <- extract_cell_counts(macro_obj, ln_group)
  
  if (nrow(cell_counts) == 0) {
    cat(ln_group, "组没有数据\n")
    return(list(
      plot_data = data.frame(),
      chisq_results = data.frame()
    ))
  }
  
  # 计算样本总数
  sample_totals <- calculate_sample_totals(cell_counts)
  
  # 计算细胞比例
  cell_proportions <- calculate_cell_proportions(cell_counts, sample_totals)
  
  # 执行统计检验
  chisq_results <- perform_statistical_tests_main(cell_counts, cell_proportions, ln_group)
  
  return(list(
    plot_data = cell_proportions,
    chisq_results = chisq_results
  ))
}

# ===== 第三部分：主数据处理流程 =====
process_box_plot_data <- function(macro_obj, celltype_alphabetical) {
  # 按顺序处理：先LN+，后LN-
  ln_groups_ordered <- c("LN+", "LN-")
  
  # 初始化列表存储结果
  plot_data_list <- list()
  stats_list <- list()
  
  for (ln_group in ln_groups_ordered) {
    group_data <- process_ln_group_data(macro_obj, ln_group)
    
    if (nrow(group_data$plot_data) > 0) {
      plot_data_list[[ln_group]] <- group_data$plot_data %>% mutate(LN_group = ln_group)
      stats_list[[ln_group]] <- group_data$chisq_results %>% mutate(LN_group = ln_group)
    }
  }
  
  # 合并数据
  if (length(plot_data_list) > 0) {
    plot_data_box_combined <- bind_rows(plot_data_list)
    chisq_results_combined <- bind_rows(stats_list)
    
    # 设置LN组因子顺序：LN+在前
    plot_data_box_combined$LN_group <- factor(
      plot_data_box_combined$LN_group,
      levels = c("LN+", "LN-")
    )
    
    chisq_results_combined$LN_group <- factor(
      chisq_results_combined$LN_group,
      levels = c("LN+", "LN-")
    )
    
    # 设置细胞类型顺序
    plot_data_box_combined$cell_type <- factor(plot_data_box_combined$cell_type, 
                                              levels = celltype_alphabetical)
    chisq_results_combined$cell_type <- factor(chisq_results_combined$cell_type, 
                                              levels = celltype_alphabetical)
    
    # 设置病理反应顺序
    plot_data_box_combined$pathological_response <- factor(
      plot_data_box_combined$pathological_response,
      levels = c("pCR", "MPR", "non-MPR")
    )
    
    return(list(
      box_data = plot_data_box_combined,
      stats_results = chisq_results_combined
    ))
  } else {
    return(NULL)
  }
}

# ===== 第四部分：绘图模块 =====

# 模块6: 创建柱状图
create_bar_plot <- function(bar_data) {
  # 设置颜色
  response_colors <- c("pCR" = "#00A087FF", "MPR" = "#F39B7FFF", "non-MPR" = "#E64B35FF")
  
  p <- ggplot(bar_data, aes(x = cell_type, y = proportion, fill = pathological_response)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    scale_fill_manual(
      values = response_colors,
      breaks = c("pCR", "MPR", "non-MPR"),
      labels = c("pCR", "MPR", "non-MPR")
    ) +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(x = "", y = "Proportion", title = "Bar Plot") +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black", face = "bold"),
      axis.text.y = element_text(size = 12, color = "black", face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold", color = "black"),
      axis.line = element_line(color = "black"),
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(size = 14, face = "bold"),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
    ) +
    facet_grid(LN_group ~ ., scales = "free_x")
  
  return(p)
}

# 模块7: 创建箱图 - 统一y轴最大值，标记放在最上面
create_box_plot <- function(box_data, stats_results) {
  if (is.null(box_data) || nrow(box_data) == 0) {
    return(ggplot() + theme_void() + labs(title = "No box plot data available"))
  }
  
  # 设置颜色
  response_colors <- c("pCR" = "#00A087FF", "MPR" = "#F39B7FFF", "non-MPR" = "#E64B35FF")
  
  # 计算整个数据集的最大值
  global_y_max <- 1.0
  
  # 为显著性标记设置固定的y位置（在图形最上方）
  stats_results$y_position <- global_y_max
  
  p <- ggplot(box_data, aes(x = cell_type, y = proportion, fill = pathological_response)) +
    geom_boxplot(width = 0.7, aes(group = interaction(cell_type, pathological_response))) +
    # 添加显著性标记 - 统一放在最上方
    geom_text(
      data = stats_results, 
      aes(x = cell_type, y = y_position, label = significance),
      size = 6, 
      inherit.aes = FALSE
    ) +
    scale_fill_manual(
      values = response_colors,
      breaks = c("pCR", "MPR", "non-MPR"),
      labels = c("pCR", "MPR", "non-MPR")
    ) +
    scale_y_continuous(
      labels = scales::percent_format(),
      limits = c(0, 1.0),
      expand = expansion(mult = c(0, 0.1))
    ) +
    labs(x = "Cell Type", y = "Cell Proportion", title = "Box Plot") +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black", face = "bold"),
      axis.text.y = element_text(size = 12, color = "black", face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold", color = "black"),
      axis.line = element_line(color = "black"),
      panel.grid = element_blank(),
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(size = 14, face = "bold"),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
    ) +
    facet_grid(LN_group ~ ., scales = "free_x")
  
  return(p)
}

# 模块8: 创建图例
create_legend <- function(bar_data) {
  response_colors <- c("pCR" = "#00A087FF", "MPR" = "#F39B7FFF", "non-MPR" = "#E64B35FF")
  
  legend_plot <- ggplot(bar_data, 
                        aes(x = cell_type, y = proportion, fill = pathological_response)) +
    geom_col(position = position_dodge(width = 0.8)) +
    scale_fill_manual(
      values = response_colors,
      breaks = c("pCR", "MPR", "non-MPR"),
      labels = c("pCR", "MPR", "non-MPR"),
      name = "Pathological Response"
    ) +
    theme(legend.position = "top",
          legend.text = element_text(size = 14, face = "bold"),
          legend.title = element_text(size = 16, face = "bold"))
  
  legend <- get_legend(legend_plot)
  return(legend)
}

# ===== 第五部分：主程序 =====

# 1. 准备柱状图数据
cat("处理LUAD数据集...\n")
macro <- subset(LUAD, subset = major_cell == "Macro")
bar_data_result <- prepare_bar_plot_data(macro)
bar_data <- bar_data_result$bar_data
celltypes <- bar_data_result$celltypes

# 2. 处理箱图数据
box_data_result <- process_box_plot_data(macro, celltypes)

# 3. 创建图形
bar_plot <- create_bar_plot(bar_data)

if (!is.null(box_data_result)) {
  box_plot <- create_box_plot(box_data_result$box_data, box_data_result$stats_results)
} else {
  box_plot <- ggplot() + 
    theme_void() + 
    labs(title = "No box plot data available")
}

legend <- create_legend(bar_data)

# 4. 组合图形
if (!is.null(box_data_result) && nrow(box_data_result$box_data) > 0) {
  combined_plot <- (bar_plot | box_plot) +
    plot_layout(guides = "collect") &
    theme(legend.position = "top")
  
  # 添加图例
  final_plot <- wrap_plots(legend, combined_plot, ncol = 1, heights = c(0.05, 0.95))
} else {
  # 只有柱状图
  final_plot <- wrap_plots(legend, bar_plot, ncol = 1, heights = c(0.05, 0.95))
}

# 5. 保存图形
pdf("LUAD_macro_combined_bar_box_plots.pdf", width = 6000/300, height = 3000/300)
print(final_plot)
dev.off()

# 同时保存TIFF格式
ggsave("LUAD_macro_combined_bar_box_plots.tiff",
       plot = final_plot,
       width = 6000/300,
       height = 3000/300,
       units = "in",
       dpi = 300,
       device = "tiff",
       compression = "lzw",
       bg = "white")
```