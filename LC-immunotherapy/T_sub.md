# T sub_cell_type
```
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggsci)
merged_seurat_obj <- readRDS("TN_anno.rds")
LUSC <- subset(merged_seurat_obj,subset=cancer_type=="LUSC")

t_cells <- subset(LUSC,subset=cell_type %in% c("Treg_1","Tex","Treg_2","Tn/Tm"))
t_cells <- NormalizeData(t_cells)
t_cells <- FindVariableFeatures(t_cells, nfeatures = 2000)
hvgs <- VariableFeatures(t_cells)
t_cells <- ScaleData(t_cells, features = hvgs)
t_cells <- RunPCA(t_cells, features = hvgs, npcs = 20)
#16227 features across 95099 samples within 1 assay
#only one dataset
#library(harmony)
#t_cells <- RunHarmony(t_cells, "orig.ident")
t_cells <- RunUMAP(t_cells, dims = 1:20)
t_cells <- FindNeighbors(t_cells, dims = 1:20)
t_cells <- FindClusters(t_cells, resolution = 0.3)

seurat_clusters <- as.character(unique(t_cells@meta.data$seurat_clusters))  # 转换为字符向量
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
npg_extended <- colorRampPalette(npg_pal)(15)
pdf("t_clusters.pdf", width = dynamic_width/300, height = base_height/300)
DimPlot(t_cells, reduction = "umap", label = TRUE, pt.size = 1, label.size = 8) +
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

t_cell_markers <- FindAllMarkers(t_cells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "wilcox")
t_significant_markers <- subset(t_cell_markers, p_val_adj < 0.05)
#write.csv(t_significant_markers, "t_all_marker.csv")
t_significant_markers <- t_significant_markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.csv(t_significant_markers, "t_top_marker_20.csv")

identity_mapping <- c(
  "0" = "CD4T_Tn/Tm_TCF7",
  "1" = "T_Tn/Tm_ANXA1",
  "2" = "CD8T_Trm_ITGAE",
  "3" = "CD8T_Tem_GZMK",
  "4" = "CD8T_Tem/Trm_IFNG", 
  "5" = "CD4T_Treg_FOXP3",
  "6" = "CD4T_Tn/Tm_BRAF",
  "7" = "CD4T_Tfh_CXCL13",
  "8" = "T_Tn/Tm_RELB",
  "9" = "CD8T_Tem_NR4A3",
  "10" = "CD4T_Tn/Tm_CCR7",
  "11" = "T_Tex/Tfh_LAG3",
  "12" = "CD8T_Tem_GZMK",
  "13"= "ILC3_KIT",
  "14" = "CD8T_Tex_HAVCR2"
)

cell_type <- identity_mapping[t_cells@meta.data$seurat_clusters]
t_cells@meta.data$cell_type <- cell_type
saveRDS(t_cells,file="t_cells_anno.rds")

library(ggplot2)
target_genes <- c("CD3","CD4", "CD8A", "CD8B")
feature_plot <- FeaturePlot(t_cells, features = target_genes)
ggsave("t_cells_featureplot.png", plot = feature_plot, width = 6 , height = 8, dpi = 300)


library(ggsci)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(15)
# 获取唯一的细胞类型并转换为字符向量
cell_types <- as.character(unique(t_cells@meta.data$cell_type))

num_legend_items <- length(cell_types)  # 图例的个数
max_label_length <- max(nchar(cell_types))  # 图例名称的最大长度

# 动态计算图片尺寸
base_width <- 5000  # 基础宽度
base_height <- 5000  # 基础高度
legend_width_factor <- 100  # 每个图例项增加的宽度
label_length_factor <- 10  # 每个字符增加的宽度

dynamic_width <- base_width + (num_legend_items * legend_width_factor) + (max_label_length * label_length_factor)

# 保存图片
pdf("t_annotation.pdf", width = dynamic_width/300, height = base_height/300)
DimPlot(t_cells, reduction = "umap", label = TRUE, pt.size = 1, label.size = 6, group.by = "cell_type") +
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    ggtitle(NULL) +
    scale_color_manual(values = npg_extended) +
    coord_fixed(ratio = 1) +
    guides(color = guide_legend(title = NULL, override.aes = list(size = 5))) +
    theme(
        text = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size =40, color = "black"),
        axis.text.y = element_text(size = 40, color = "black"),
        axis.title.x = element_text(size = 56, face = "bold", color = "black"),
        axis.title.y = element_text(size = 56, face = "bold", color = "black", margin = margin(r = 20)),  # 增加右侧间距
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 24, face = "bold", color = "black"),
        legend.title = element_text(size = 24, face = "bold", color = "black"),
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
    ) +scale_x_continuous(limits = c(-10.5,10))
dev.off()
```
# group difference
```
#分布图
library(Seurat)
library(ggplot2)
library(ggsci)
library(dplyr)
plot_data <- t_cells

plot_data$group_response <- paste(plot_data$Group, plot_data$Response, sep = "_")
plot_data$group_response <- factor(plot_data$group_response,
                                  levels = c("Group1_YES", "Group1_NO",
                                            "Group2_YES", "Group2_NO"))
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(15)
pdf("t_combined_annotation_by_group_response.pdf", width = 12000/300, height = 3000/300)

p <- DimPlot(plot_data, 
             reduction = "umap", 
             label = TRUE, 
             pt.size = 1, 
             group.by = "cell_type", 
             label.size = 4,
             split.by = "group_response",  
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
                                 margin = margin(b = 15)),
        panel.spacing = unit(1.5, "lines")
    )

print(p)
dev.off()

#折线图
plot_data <- t_cells@meta.data %>%
  group_by(Group, Response, cell_type) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  # 计算每个Group-Response组合内的比例
  group_by(Group, Response) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()

# 设置Response因子的顺序 - YES在前，NO在后
plot_data$Response <- factor(plot_data$Response, levels = c("YES", "NO"))

# 绘制折线图 - 去掉分面标题方框
pdf("t_group_response_lineplot.pdf", width = 2000/300, height = 4000/300)

ggplot(plot_data, aes(x = Response, y = proportion, group = cell_type, color = cell_type)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) +
  facet_wrap(~ Group, ncol = 1, scales = "free_y") +  # 添加 scales = "free_y"
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

#柱状图
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggsci)

cell_proportions <- t_cells@meta.data %>%
  group_by(Group, cell_type) %>%
  summarise(count = n()) %>%
  group_by(Group) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()

npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(15)

all_cell_types <- unique(t_cells@meta.data$cell_type)
cell_proportions_complete <- cell_proportions %>%
  complete(Group, cell_type = all_cell_types, fill = list(count = 0, proportion = 0))

# 绘制分面柱状图 - 图例多行显示
pdf("t_cells_group_proportion.pdf", width = 3000/300, height = 3000/300)  # 适当增加高度

p <- ggplot(cell_proportions_complete, aes(x = cell_type, y = proportion, fill = cell_type)) +
  geom_col() +
  scale_fill_manual(values = npg_extended, name = "Cell Type") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "", y = "Proportion") +
  theme_classic() +
  theme(
    # 横坐标设置
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    
    # 纵坐标设置
    axis.text.y = element_text(size = 18, color = "black", face = "bold"),
    axis.title.y = element_text(size = 20, face = "bold", color = "black", margin = margin(r = 15)),
    
    # 图例位置 - 多行显示
    legend.position = "top",
    legend.direction = "horizontal",  # 保持水平方向
    legend.justification = "center",
    legend.text = element_text(size =22, face = "bold"),
    legend.title = element_blank(),
    legend.box = "vertical",  # 改为垂直排列多行
    legend.key.size = unit(0.8, "cm"),
    legend.spacing.x = unit(0.5, "cm"),
    legend.box.margin = margin(10, 10, 10, 10),
    legend.box.just = "center",
    
    # 分面标题设置
    strip.background = element_blank(),
    strip.text = element_text(size = 18, face = "bold", vjust = 1),
    
    # 其他设置
    panel.spacing = unit(1, "lines"),
    plot.margin = margin(100, 20, 20, 20)  # 增大顶部边距给多行图例更多空间
  ) +
  guides(fill = guide_legend(nrow = 5, byrow = TRUE)) +  # 关键：设置图例显示为3行
  facet_wrap(~ Group, ncol = 1, scales = "free_x")

print(p)
dev.off()

#图例在右
pdf("t_cells_group_proportion.pdf", width = 5000/300, height = 3000/300)

p <- ggplot(cell_proportions_complete, aes(x = cell_type, y = proportion, fill = cell_type)) +
    geom_col() +
    scale_fill_manual(values = npg_extended, name = "Cell Type") +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(x = "", y = "Proportion") +
    theme_classic() +
    theme(
        # 横坐标设置
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        
        # 纵坐标设置
        axis.text.y = element_text(size = 18, color = "black", face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold", color = "black", margin = margin(r = 15)),
        
        # 图例位置 - 右侧两列显示
        legend.position = "right",  # 改为右侧
        legend.direction = "vertical",  # 改为垂直方向
        legend.justification = "top",
        legend.text = element_text(size = 16, face = "bold"),  # 适当减小文字大小
        legend.title = element_blank(),
        legend.box = "vertical",
        legend.key.size = unit(0.6, "cm"),  # 适当减小图例键大小
        legend.spacing.y = unit(0.3, "cm"),  # 垂直间距
        legend.box.margin = margin(10, 10, 10, 10),
        
        # 分面标题设置
        strip.background = element_blank(),
        strip.text = element_text(size = 18, face = "bold", vjust = 1),
        
        # 其他设置
        panel.spacing = unit(1, "lines"),
        plot.margin = margin(20, 20, 20, 20)  # 调整边距
    ) +
    guides(fill = guide_legend(ncol = 1, byrow = TRUE)) +  # 关键：设置图例显示为2列
    facet_wrap(~ Group, ncol = 1, scales = "free_x")

print(p)
dev.off()

#箱图
library(ggpubr)
library(dplyr)
library(ggsci)
library(tidyr)

cell_counts <- t_cells@meta.data %>%
  group_by(patients, cell_type, Group) %>%
  summarise(count = n(), .groups = "drop")

# 计算每个样本的总细胞数
sample_totals <- cell_counts %>%
  group_by(patients) %>%
  summarise(total = sum(count), .groups = "drop")

# 计算细胞类型比例
cell_proportions <- cell_counts %>%
  left_join(sample_totals, by = "patients") %>%
  mutate(proportion = count / total)

# === 找出哪个细胞类型在Group2中缺失 ===
# 检查每个细胞类型在两个Group中的存在情况
celltype_group_presence <- cell_proportions %>%
  group_by(cell_type, Group) %>%
  summarise(has_data = n() > 0, .groups = "drop") %>%
  pivot_wider(names_from = Group, values_from = has_data, values_fill = FALSE)

# 找出在Group1中存在但在Group2中缺失的细胞类型
missing_in_group2 <- celltype_group_presence %>%
  filter(Group1 == TRUE & Group2 == FALSE) %>%
  pull(cell_type)

cat("在Group2中缺失的细胞类型:", paste(missing_in_group2, collapse = ", "), "\n")

# === 只为缺失的细胞类型在Group2中添加0值 ===
if (length(missing_in_group2) > 0) {
  # 获取Group2的所有患者
  group2_patients <- unique(cell_proportions$patients[cell_proportions$Group == "Group2"])
  
  # 只为缺失的细胞类型创建0值记录
  missing_records <- expand_grid(
    patients = group2_patients,
    cell_type = missing_in_group2,
    Group = "Group2"
  ) %>%
  left_join(sample_totals, by = "patients") %>%
  mutate(count = 0, proportion = 0)
  
  # 合并到绘图数据
  plot_data <- cell_proportions %>% bind_rows(missing_records)
} else {
  plot_data <- cell_proportions
}

# === 关键修复：确保分组数量匹配 ===
# 获取所有细胞类型和分组
all_celltypes <- unique(plot_data$cell_type)
all_groups <- unique(plot_data$Group)
n_groups <- length(all_groups)

# 构建完整的细胞类型-分组组合并填充0
full_combinations <- expand.grid(
  cell_type = all_celltypes,
  Group = all_groups,
  stringsAsFactors = FALSE
)

# 按细胞类型和分组汇总计数
cell_group_counts <- cell_counts %>%
  group_by(cell_type, Group) %>%
  summarise(total_count = sum(count), .groups = "drop") %>%
  right_join(full_combinations, by = c("cell_type", "Group")) %>%
  mutate(total_count = replace(total_count, is.na(total_count), 0)) %>%
  arrange(cell_type, Group)  # 确保分组顺序一致

# 计算每个分组的总细胞数（全局）
global_group_totals <- cell_counts %>%
  group_by(Group) %>%
  summarise(group_total = sum(count), .groups = "drop") %>%
  right_join(data.frame(Group = all_groups), by = "Group") %>%
  mutate(group_total = replace(group_total, is.na(group_total), 0)) %>%
  arrange(Group)  # 保持分组顺序一致

# === 修复卡方检验逻辑 ===
chisq_results <- cell_group_counts %>%
  group_by(cell_type) %>%
  summarise(
    p_value = {
      # 获取当前细胞类型在所有分组的计数（顺序与global_group_totals一致）
      obs <- total_count[match(all_groups, Group)]
      
      # 获取对应分组的总细胞数
      group_tot <- global_group_totals$group_total
      
      # 检查数据有效性
      if (sum(obs) < 5 || length(obs) != length(group_tot)) {
        NA_real_
      } else {
        # 计算期望比例（基于全局分组总细胞数）
        expected_prop <- group_tot / sum(group_tot)
        
        # 处理两组情况（Fisher精确检验优先）
        if (length(obs) == 2) {
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
        } else {
          # 多组情况使用卡方检验
          if (any(obs * expected_prop < 5)) {
            NA_real_  # 期望频数不足时返回NA
          } else {
            chisq.test(obs, p = expected_prop)$p.value
          }
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

# === 绘制图形 ===
y_max <- max(plot_data$proportion, na.rm = TRUE) * 1.1

pdf("t_cells_箱图_proportion_by_group.pdf", width = 6000/300, height = 3000/300)
ggplot(plot_data, aes(x = cell_type, y = proportion, fill = Group)) +
  geom_boxplot(width = 0.7) +
  geom_text(
    data = chisq_results, 
    aes(x = cell_type, y = y_max, label = significance),
    size = 12, 
    inherit.aes = FALSE
  ) +
  labs(x = "", y = "Cell Proportion", fill = "Group") +
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


#箱图和柱状图二合一
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(patchwork)

# 定义细胞类型分类
identity_mapping <- c(
  "0" = "CD4T_Tn/Tm_TCF7",
  "1" = "T_Tn/Tm_ANXA1",
  "2" = "CD8T_Trm_ITGAE", 
  "3" = "CD8T_Tem_GZMK",
  "4" = "CD8T_Tem/Trm_IFNG",
  "5" = "CD4T_Treg_FOXP3",
  "6" = "CD4T_Tn/Tm_BRAF",
  "7" = "CD4T_Tfh_CXCL13",
  "8" = "T_Tn/Tm_RELB",
  "9" = "CD8T_Tem_NR4A3",
  "10" = "CD4T_Tn/Tm_CCR7",
  "11" = "T_Tex/Tfh_LAG3",
  "12" = "CD8T_Tem_GZMK",
  "13" = "ILC3_KIT",
  "14" = "CD8T_Tex_HAVCR2"
)

# ===== 第一部分：准备柱状图数据 =====
cell_proportions_bar <- t_cells@meta.data %>%
  group_by(Group, cell_type) %>%
  summarise(count = n()) %>%
  group_by(Group) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()

# 分类细胞类型
cell_proportions_bar <- cell_proportions_bar %>%
  mutate(
    cell_category = case_when(
      grepl("^CD4T", cell_type) ~ "CD4 T Cells",
      grepl("^CD8T", cell_type) ~ "CD8 T Cells", 
      TRUE ~ "Other T Cells"
    ),
    cell_category = factor(cell_category, levels = c("CD4 T Cells", "CD8 T Cells", "Other T Cells"))
  )

# 转换为宽格式
cell_wide <- cell_proportions_bar %>%
  select(Group, cell_type, proportion, cell_category) %>%
  pivot_wider(names_from = Group, values_from = proportion, values_fill = 0)

# 转换回长格式用于绘图
plot_data_bar <- cell_wide %>%
  pivot_longer(cols = c(Group1, Group2), 
               names_to = "Group", 
               values_to = "proportion") %>%
  mutate(Group = factor(Group, levels = c("Group1", "Group2")))

# ===== 第二部分：准备箱图数据 =====
cell_counts <- t_cells@meta.data %>%
  group_by(patients, cell_type, Group) %>%
  summarise(count = n(), .groups = "drop")

# 计算每个样本的总细胞数
sample_totals <- cell_counts %>%
  group_by(patients) %>%
  summarise(total = sum(count), .groups = "drop")

# 计算细胞类型比例
cell_proportions_box <- cell_counts %>%
  left_join(sample_totals, by = "patients") %>%
  mutate(proportion = count / total)

# 分类细胞类型
cell_proportions_box <- cell_proportions_box %>%
  mutate(
    cell_category = case_when(
      grepl("^CD4T", cell_type) ~ "CD4 T Cells",
      grepl("^CD8T", cell_type) ~ "CD8 T Cells", 
      TRUE ~ "Other T Cells"
    ),
    cell_category = factor(cell_category, levels = c("CD4 T Cells", "CD8 T Cells", "Other T Cells"))
  )

# 找出在Group2中缺失的细胞类型并添加0值
celltype_group_presence <- cell_proportions_box %>%
  group_by(cell_type, Group, cell_category) %>%
  summarise(has_data = n() > 0, .groups = "drop") %>%
  pivot_wider(names_from = Group, values_from = has_data, values_fill = FALSE)

missing_in_group2 <- celltype_group_presence %>%
  filter(Group1 == TRUE & Group2 == FALSE) %>%
  pull(cell_type)

if (length(missing_in_group2) > 0) {
  group2_patients <- unique(cell_proportions_box$patients[cell_proportions_box$Group == "Group2"])
  
  missing_records <- expand_grid(
    patients = group2_patients,
    cell_type = missing_in_group2,
    Group = "Group2"
  ) %>%
  left_join(sample_totals, by = "patients") %>%
  mutate(count = 0, proportion = 0)
  
  missing_records <- missing_records %>%
    mutate(
      cell_category = case_when(
        grepl("^CD4T", cell_type) ~ "CD4 T Cells",
        grepl("^CD8T", cell_type) ~ "CD8 T Cells", 
        TRUE ~ "Other T Cells"
      )
    )
  
  plot_data_box <- cell_proportions_box %>% bind_rows(missing_records)
} else {
  plot_data_box <- cell_proportions_box
}

# 统计检验
all_celltypes <- unique(plot_data_box$cell_type)
all_groups <- unique(plot_data_box$Group)

full_combinations <- expand.grid(
  cell_type = all_celltypes,
  Group = all_groups,
  stringsAsFactors = FALSE
)

cell_group_counts <- cell_counts %>%
  group_by(cell_type, Group) %>%
  summarise(total_count = sum(count), .groups = "drop") %>%
  right_join(full_combinations, by = c("cell_type", "Group")) %>%
  mutate(total_count = replace(total_count, is.na(total_count), 0)) %>%
  arrange(cell_type, Group)

global_group_totals <- cell_counts %>%
  group_by(Group) %>%
  summarise(group_total = sum(count), .groups = "drop") %>%
  right_join(data.frame(Group = all_groups), by = "Group") %>%
  mutate(group_total = replace(group_total, is.na(group_total), 0)) %>%
  arrange(Group)

chisq_results <- cell_group_counts %>%
  group_by(cell_type) %>%
  summarise(
    p_value = {
      obs <- total_count[match(all_groups, Group)]
      group_tot <- global_group_totals$group_total
      
      if (sum(obs) < 5 || length(obs) != length(group_tot)) {
        NA_real_
      } else {
        expected_prop <- group_tot / sum(group_tot)
        
        if (length(obs) == 2) {
          other_counts <- group_tot - obs
          cont_table <- matrix(c(obs[1], other_counts[1],
                                 obs[2], other_counts[2]),
                               nrow = 2, byrow = TRUE)
          
          if (any(cont_table < 0) || sum(cont_table) == 0) {
            NA_real_
          } else if (any(cont_table < 5)) {
            fisher.test(cont_table)$p.value
          } else {
            chisq.test(cont_table)$p.value
          }
        } else {
          if (any(obs * expected_prop < 5)) {
            NA_real_
          } else {
            chisq.test(obs, p = expected_prop)$p.value
          }
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

chisq_results <- chisq_results %>%
  mutate(
    cell_category = case_when(
      grepl("^CD4T", cell_type) ~ "CD4 T Cells",
      grepl("^CD8T", cell_type) ~ "CD8 T Cells", 
      TRUE ~ "Other T Cells"
    ),
    cell_category = factor(cell_category, levels = c("CD4 T Cells", "CD8 T Cells", "Other T Cells"))
  )

# 按首字母排序
celltype_alphabetical <- sort(all_celltypes)
plot_data_bar$cell_type <- factor(plot_data_bar$cell_type, levels = celltype_alphabetical)
plot_data_box$cell_type <- factor(plot_data_box$cell_type, levels = celltype_alphabetical)
chisq_results$cell_type <- factor(chisq_results$cell_type, levels = celltype_alphabetical)

# ===== 第三部分：创建每个类别的图形 =====
create_plots_for_category <- function(category) {
  
  # 筛选当前类别的数据
  bar_data <- plot_data_bar %>% filter(cell_category == category)
  box_data <- plot_data_box %>% filter(cell_category == category)
  chisq_data <- chisq_results %>% filter(cell_category == category)
  
  y_max_box <- max(box_data$proportion, na.rm = TRUE) * 1.1
  
  # 柱状图
  p_bar <- ggplot(bar_data, aes(x = cell_type, y = proportion, fill = Group)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    scale_fill_manual(values = c("Group1" = "#E64B35FF", "Group2" = "#4DBBD5FF")) +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(x = "", y = "Proportion", title = paste(category, "- Bar Plot")) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 28, color = "black", face = "bold"),
      axis.text.y = element_text(size =28, color = "black", face = "bold"),
      axis.title.y = element_text(size = 28, face = "bold", color = "black"),
      legend.position = "none",
      plot.title = element_text(size =28, face = "bold", hjust = 0.5)
    )
  
  # 箱图
  p_box <- ggplot(box_data, aes(x = cell_type, y = proportion, fill = Group)) +
    geom_boxplot(width = 0.7) +
    geom_text(
      data = chisq_data, 
      aes(x = cell_type, y = y_max_box, label = significance),
      size = 4, 
      inherit.aes = FALSE
    ) +
    scale_fill_manual(values = c("Group1" = "#E64B35FF", "Group2" = "#4DBBD5FF")) +
    scale_y_continuous(labels = scales::percent_format(), limits = c(0, y_max_box)) +
    labs(x = "", y = "Cell Proportion", title = paste(category, "- Box Plot")) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 28, face = "bold"),
      axis.text.y = element_text(size = 28, face = "bold"),
      axis.title.y = element_text(size = 28, face = "bold"),
      legend.position = "none",
      plot.title = element_text(size =28, face = "bold", hjust = 0.5)
    )
  
  return(list(bar = p_bar, box = p_box))
}

# 为每个类别创建图形
cd4_plots <- create_plots_for_category("CD4 T Cells")
cd8_plots <- create_plots_for_category("CD8 T Cells")
other_plots <- create_plots_for_category("Other T Cells")

# ===== 第四部分：组合图形 =====
# 创建图例
legend_plot <- ggplot(plot_data_bar, aes(x = cell_type, y = proportion, fill = Group)) +
  geom_col() +
  scale_fill_manual(values = c("Group1" = "#E64B35FF", "Group2" = "#4DBBD5FF")) +
  theme(legend.position = "top",
        legend.text = element_text(size = 28, face = "bold"),
        legend.title = element_text(size = 28, face = "bold"))

legend <- get_legend(legend_plot)

# 组合图形：三行两列
combined_plot <- (cd4_plots$bar | cd4_plots$box) /
                 (cd8_plots$bar | cd8_plots$box) /
                 (other_plots$bar | other_plots$box)

# 添加图例和标题
final_plot <- wrap_plots(legend, combined_plot, ncol = 1, heights = c(0.05, 0.95))

# 保存图形
pdf("t_cells_combined_bar_box.pdf", width = 8000/300, height = 8000/300)
print(final_plot)
dev.off()

#气泡图
# 准备数据
cell_proportions <- t_cells@meta.data %>%
  group_by(Group, Response, cell_type) %>%
  summarise(count = n()) %>%
  group_by(Group, Response) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()

# 设置Response的顺序 - YES在前，NO在后
cell_proportions$Response <- factor(cell_proportions$Response, levels = c("YES", "NO"))

# 创建组合标签并设置正确的顺序
cell_proportions <- cell_proportions %>%
  mutate(
    group_response = paste0(Group, "_", Response)
  )

# 设置分组顺序
group_response_order <- c(
  "Group1_YES", "Group1_NO",
  "Group2_YES", "Group2_NO"
)

cell_proportions$group_response <- factor(cell_proportions$group_response, 
                                         levels = group_response_order)

# 设置颜色
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(15)

# 绘制气泡图
pdf("t_cells_bubble_plot.pdf", width =7000/300, height = 2000/300)

p <- ggplot(cell_proportions, aes(x = cell_type, y = group_response, 
                                 size = proportion, color = cell_type)) +
  geom_point(alpha = 0.8) +
  scale_size_continuous(name = "Proportion", 
                       range = c(1, 10),  # 气泡大小范围
                       limits = c(0, 1),  # 固定上限为20%
                       breaks = c(0.05, 0.1, 0.2, 0.4,0.8),  # 固定显示5%, 10%, 15%, 20%
                       labels = c("5%", "10%", "20%", "40%","80%")) +  # 自定义标签
  scale_color_manual(values = npg_extended) +
  scale_y_discrete(limits = rev(levels(cell_proportions$group_response))) +
  labs(x = "Cell Type", y = "Group & Response", 
       title = "T Cells Proportions Across Groups and Responses") +
  theme_classic() +
  theme(
    text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20, face = "bold", color = "black"),
    axis.text.y = element_text(size = 20, face = "bold", color = "black"),
    axis.title.x = element_text(size = 24, face = "bold", color = "black", margin = margin(t = 10)),
    axis.title.y = element_text(size = 24, face = "bold", color = "black", margin = margin(r = 20)),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16, face = "bold"),
    legend.position = "right",
    legend.key.size = unit(0.8, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    panel.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    plot.title = element_text(size = 28, face = "bold", hjust = 0.5)
  ) +
  guides(color = "none")

print(p)
dev.off()

```
# response
```
#箱图
#箱图 - Group1的Response
library(ggpubr)
library(dplyr)
library(ggsci)
library(tidyr)

# 只处理Group1的数据
current_group <- "Group1"

cat("正在处理:", current_group, "\n")

# 筛选Group1的数据
cell_counts <- t_cells@meta.data %>%
  filter(Group == current_group) %>%
  group_by(patients, cell_type, Response) %>%
  summarise(count = n(), .groups = "drop")

# 计算每个样本的总细胞数（在Group1内）
sample_totals <- cell_counts %>%
  group_by(patients) %>%
  summarise(total = sum(count), .groups = "drop")

# 计算细胞类型比例
cell_proportions <- cell_counts %>%
  left_join(sample_totals, by = "patients") %>%
  mutate(proportion = count / total)

# === 设置Response的顺序为YES, NO ===
cell_proportions$Response <- factor(cell_proportions$Response, levels = c("YES", "NO"))

# === 找出哪个细胞类型在某个Response中缺失 ===
celltype_response_presence <- cell_proportions %>%
  group_by(cell_type, Response) %>%
  summarise(has_data = n() > 0, .groups = "drop") %>%
  pivot_wider(names_from = Response, values_from = has_data, values_fill = FALSE)

# 获取所有Response类型（按照YES, NO顺序）
all_responses <- c("YES", "NO")

# 找出在某个Response中存在但在其他Response中缺失的细胞类型
missing_in_NO <- celltype_response_presence %>%
  filter(YES == TRUE & NO == FALSE) %>%
  pull(cell_type)

missing_in_YES <- celltype_response_presence %>%
  filter(YES == FALSE & NO == TRUE) %>%
  pull(cell_type)

cat("在NO中缺失的细胞类型:", paste(missing_in_NO, collapse = ", "), "\n")
cat("在YES中缺失的细胞类型:", paste(missing_in_YES, collapse = ", "), "\n")

# === 为缺失的细胞类型添加0值 ===
plot_data <- cell_proportions

# 为NO缺失的细胞类型在NO中添加0值
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

# 为YES缺失的细胞类型在YES中添加0值
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

# === 关键修复：确保分组数量匹配 ===
all_celltypes <- unique(plot_data$cell_type)

# 确保Response顺序为YES, NO
plot_data$Response <- factor(plot_data$Response, levels = c("YES", "NO"))

# 构建完整的细胞类型-Response组合并填充0
full_combinations <- expand_grid(
  cell_type = all_celltypes,
  Response = factor(c("YES", "NO"), levels = c("YES", "NO"))
)

# 按细胞类型和Response汇总计数
cell_response_counts <- cell_counts %>%
  group_by(cell_type, Response) %>%
  summarise(total_count = sum(count), .groups = "drop") %>%
  right_join(full_combinations, by = c("cell_type", "Response")) %>%
  mutate(total_count = replace(total_count, is.na(total_count), 0)) %>%
  arrange(cell_type, Response)

# 计算每个Response的总细胞数（在Group1内）
global_response_totals <- cell_counts %>%
  group_by(Response) %>%
  summarise(response_total = sum(count), .groups = "drop") %>%
  right_join(data.frame(Response = factor(c("YES", "NO"), levels = c("YES", "NO"))), by = "Response") %>%
  mutate(response_total = replace(response_total, is.na(response_total), 0)) %>%
  arrange(Response)

# === 修复卡方检验逻辑 - 确保顺序不影响结果 ===
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

# === 绘制图形 ===
y_max <- max(plot_data$proportion, na.rm = TRUE) * 1.1

pdf("t_cells_Group1_response_箱图.pdf", width = 6000/300, height = 3000/300)
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

cat("已完成: Group1\n")

#箱图 - Group2的Response
library(ggpubr)
library(dplyr)
library(ggsci)
library(tidyr)

# 只处理Group2的数据
current_group <- "Group2"

cat("正在处理:", current_group, "\n")

# 筛选Group2的数据
cell_counts <- t_cells@meta.data %>%
  filter(Group == current_group) %>%
  group_by(patients, cell_type, Response) %>%
  summarise(count = n(), .groups = "drop")

# 计算每个样本的总细胞数（在Group2内）
sample_totals <- cell_counts %>%
  group_by(patients) %>%
  summarise(total = sum(count), .groups = "drop")

# 计算细胞类型比例
cell_proportions <- cell_counts %>%
  left_join(sample_totals, by = "patients") %>%
  mutate(proportion = count / total)

# === 设置Response的顺序为YES, NO ===
cell_proportions$Response <- factor(cell_proportions$Response, levels = c("YES", "NO"))

# === 找出哪个细胞类型在某个Response中缺失 ===
celltype_response_presence <- cell_proportions %>%
  group_by(cell_type, Response) %>%
  summarise(has_data = n() > 0, .groups = "drop") %>%
  pivot_wider(names_from = Response, values_from = has_data, values_fill = FALSE)

# 获取所有Response类型（按照YES, NO顺序）
all_responses <- c("YES", "NO")

# 找出在某个Response中存在但在其他Response中缺失的细胞类型
missing_in_NO <- celltype_response_presence %>%
  filter(YES == TRUE & NO == FALSE) %>%
  pull(cell_type)

missing_in_YES <- celltype_response_presence %>%
  filter(YES == FALSE & NO == TRUE) %>%
  pull(cell_type)

cat("在NO中缺失的细胞类型:", paste(missing_in_NO, collapse = ", "), "\n")
cat("在YES中缺失的细胞类型:", paste(missing_in_YES, collapse = ", "), "\n")

# === 为缺失的细胞类型添加0值 ===
plot_data <- cell_proportions

# 为NO缺失的细胞类型在NO中添加0值
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

# 为YES缺失的细胞类型在YES中添加0值
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

# === 关键修复：确保分组数量匹配 ===
all_celltypes <- unique(plot_data$cell_type)

# 确保Response顺序为YES, NO
plot_data$Response <- factor(plot_data$Response, levels = c("YES", "NO"))

# 构建完整的细胞类型-Response组合并填充0
full_combinations <- expand_grid(
  cell_type = all_celltypes,
  Response = factor(c("YES", "NO"), levels = c("YES", "NO"))
)

# 按细胞类型和Response汇总计数
cell_response_counts <- cell_counts %>%
  group_by(cell_type, Response) %>%
  summarise(total_count = sum(count), .groups = "drop") %>%
  right_join(full_combinations, by = c("cell_type", "Response")) %>%
  mutate(total_count = replace(total_count, is.na(total_count), 0)) %>%
  arrange(cell_type, Response)

# 计算每个Response的总细胞数（在Group2内）
global_response_totals <- cell_counts %>%
  group_by(Response) %>%
  summarise(response_total = sum(count), .groups = "drop") %>%
  right_join(data.frame(Response = factor(c("YES", "NO"), levels = c("YES", "NO"))), by = "Response") %>%
  mutate(response_total = replace(response_total, is.na(response_total), 0)) %>%
  arrange(Response)

# === 修复卡方检验逻辑 - 确保顺序不影响结果 ===
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

# === 绘制图形 ===
y_max <- max(plot_data$proportion, na.rm = TRUE) * 1.1

pdf("t_cells_Group2_response_箱图.pdf", width = 6000/300, height = 3000/300)
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

cat("已完成: Group2\n")

#桑葚图
# 准备包含细胞类型的数据
library(ggalluvial)
sankey_data_celltype <- t_cells@meta.data %>%
  group_by(Group, cell_type) %>%  # 移除Response
  summarise(cell_count = n(), .groups = 'drop') %>%
  group_by(Group) %>%  # 只按Group分组
  mutate(total = sum(cell_count),
         proportion = cell_count / total) %>%
  ungroup()

# 两层桑葚图
pdf("t_cells_sankey_2layer.pdf", width = 6000/300, height = 4000/300)

ggplot(sankey_data_celltype,
       aes(axis1 = Group, axis2 = cell_type,  # 只有两层
           y = cell_count)) +
  geom_alluvium(aes(fill = cell_type), width = 1/8) +
  geom_stratum(width = 1/8, fill = "grey", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), 
            size = 8, fontface = "bold") +  # 增大文字大小
  scale_x_discrete(limits = c("Group", "Cell Type"),  # 只有两个层级
                   expand = c(0.1, 0.1)) +
  scale_fill_manual(values = npg_extended) +
  labs(title = "Cell Distribution: Group → Cell Type",  # 修改标题
       y = "Number of Cells") +
  theme_void() +
  theme(
    plot.title = element_text(size = 32, face = "bold", hjust = 0.5, margin = margin(b = 20)),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18),
    legend.position = "right",
    plot.margin = margin(20, 20, 20, 20)
  )

dev.off()


library(ggalluvial)
library(dplyr)
library(ggsci)

# 准备包含patients的数据
sankey_data_patients <- t_cells@meta.data %>%
  group_by(patients, Group, Response, cell_type) %>%
  summarise(cell_count = n(), .groups = 'drop') %>%
  group_by(patients, Group, Response) %>%
  mutate(total = sum(cell_count),
         proportion = cell_count / total) %>%
  ungroup()

# 四层桑葚图
pdf("t_cells_sankey_4layer_patients.pdf", width = 10000/300, height = 6000/300)

ggplot(sankey_data_patients,
       aes(axis1 = patients, axis2 = Group, axis3 = Response, axis4 = cell_type, 
           y = cell_count)) +
  geom_alluvium(aes(fill = cell_type), width = 1/12) +
  geom_stratum(width = 1/12, fill = "grey", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), 
            size = 3, fontface = "bold") +
  scale_x_discrete(limits = c("Patients", "Group", "Response", "Cell Type"), 
                   expand = c(0.02, 0.02)) +
  scale_fill_manual(values = npg_extended) +
  labs(title = "Cell Distribution: Patients → Group → Response → Cell Type",
       y = "Number of Cells") +
  theme_void() +
  theme(
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.position = "right"
  )

dev.off()
```
# CD4/CD8
```
#箱图
library(ggpubr)
library(dplyr)
library(ggsci)
library(tidyr)
library(patchwork)

# 定义细胞类型分类函数
classify_cell_type <- function(cell_type) {
  case_when(
    grepl("^CD4T", cell_type) ~ "CD4 T Cells",
    grepl("^CD8T", cell_type) ~ "CD8 T Cells", 
    TRUE ~ "Other T Cells"
  )
}

# 处理两组数据的函数
process_group_data <- function(current_group) {
  cat("正在处理:", current_group, "\n")
  
  # 筛选数据
  cell_counts <- t_cells@meta.data %>%
    filter(Group == current_group) %>%
    group_by(patients, cell_type, Response) %>%
    summarise(count = n(), .groups = "drop")
  
  # 计算每个样本的总细胞数
  sample_totals <- cell_counts %>%
    group_by(patients) %>%
    summarise(total = sum(count), .groups = "drop")
  
  # 计算细胞类型比例
  cell_proportions <- cell_counts %>%
    left_join(sample_totals, by = "patients") %>%
    mutate(proportion = count / total)
  
  # === 关键修改：确保Response按YES/NO顺序 ===
  cell_proportions$Response <- factor(cell_proportions$Response, levels = c("YES", "NO"))
  
  # 分类细胞类型
  cell_proportions <- cell_proportions %>%
    mutate(cell_category = classify_cell_type(cell_type))
  
  # 找出缺失的细胞类型
  celltype_response_presence <- cell_proportions %>%
    group_by(cell_type, Response) %>%
    summarise(has_data = n() > 0, .groups = "drop") %>%
    pivot_wider(names_from = Response, values_from = has_data, values_fill = FALSE)
  
  missing_in_NO <- celltype_response_presence %>%
    filter(YES == TRUE & NO == FALSE) %>%
    pull(cell_type)
  
  missing_in_YES <- celltype_response_presence %>%
    filter(YES == FALSE & NO == TRUE) %>%
    pull(cell_type)
  
  plot_data <- cell_proportions
  
  # 为缺失的细胞类型添加0值（用于创建横线占位）
  if (length(missing_in_NO) > 0) {
    NO_patients <- unique(cell_proportions$patients[cell_proportions$Response == "NO"])
    missing_records_NO <- expand_grid(
      patients = NO_patients,
      cell_type = missing_in_NO,
      Response = "NO"
    ) %>%
    left_join(sample_totals, by = "patients") %>%
    mutate(count = 0, proportion = 0, cell_category = classify_cell_type(cell_type))
    
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
    mutate(count = 0, proportion = 0, cell_category = classify_cell_type(cell_type))
    
    plot_data <- plot_data %>% bind_rows(missing_records_YES)
  }
  
  # 统计检验
  all_celltypes <- unique(plot_data$cell_type)
  
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
        yes_count <- total_count[Response == "YES"]
        no_count <- total_count[Response == "NO"]
        yes_total <- global_response_totals$response_total[global_response_totals$Response == "YES"]
        no_total <- global_response_totals$response_total[global_response_totals$Response == "NO"]
        
        if (sum(c(yes_count, no_count)) < 5) {
          NA_real_
        } else {
          cont_table <- matrix(c(yes_count, yes_total - yes_count,
                                 no_count, no_total - no_count),
                               nrow = 2, byrow = FALSE)
          
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
  
  # 添加细胞类别信息
  chisq_results <- chisq_results %>%
    mutate(cell_category = classify_cell_type(cell_type))
  
  return(list(plot_data = plot_data, chisq_results = chisq_results))
}

# 处理两组数据
group1_data <- process_group_data("Group1")
group2_data <- process_group_data("Group2")

# 合并数据并添加组别信息
plot_data_combined <- bind_rows(
  group1_data$plot_data %>% mutate(Group = "Group1"),
  group2_data$plot_data %>% mutate(Group = "Group2")
)

chisq_results_combined <- bind_rows(
  group1_data$chisq_results %>% mutate(Group = "Group1"),
  group2_data$chisq_results %>% mutate(Group = "Group2")
)

# 获取所有细胞类型并按类别分组
cell_categories <- c("CD4 T Cells", "CD8 T Cells", "Other T Cells")

# 为每个类别获取完整的细胞类型列表（合并两个Group的所有细胞类型）
cd4_cells_all <- sort(unique(plot_data_combined$cell_type[plot_data_combined$cell_category == "CD4 T Cells"]))
cd8_cells_all <- sort(unique(plot_data_combined$cell_type[plot_data_combined$cell_category == "CD8 T Cells"]))
other_cells_all <- sort(unique(plot_data_combined$cell_type[plot_data_combined$cell_category == "Other T Cells"]))

# 为每个Group和细胞类别创建独立的图形
create_group_category_plot <- function(group_name, category, cell_types, show_x_axis = FALSE) {
  
  # 筛选当前Group和类别的数据
  group_category_data <- plot_data_combined %>% 
    filter(Group == group_name, cell_category == category)
  
  group_category_chisq <- chisq_results_combined %>% 
    filter(Group == group_name, cell_category == category)
  
  # === 关键修改：确保Response因子水平 ===
  group_category_data$Response <- factor(group_category_data$Response, levels = c("YES", "NO"))
  
  # 设置细胞类型的因子水平
  group_category_data$cell_type <- factor(group_category_data$cell_type, levels = cell_types)
  group_category_chisq$cell_type <- factor(group_category_chisq$cell_type, levels = cell_types)
  
  # 计算y轴最大值
  y_max <- ifelse(nrow(group_category_data) > 0, 
                  max(group_category_data$proportion, na.rm = TRUE) * 1.1, 0.1)
  
  # 为proportion=0的数据添加横线
  zero_data <- group_category_data %>% filter(proportion == 0)
  
  # 创建当前Group和类别的图形
  p <- ggplot(group_category_data, aes(x = cell_type, y = proportion, fill = Response)) +
    geom_boxplot(width = 0.7) +
    # 添加0值横线
    geom_segment(
      data = zero_data,
      aes(x = as.numeric(cell_type) - 0.35, xend = as.numeric(cell_type) + 0.35,
          y = 0, yend = 0),
      color = "black", linewidth = 1, inherit.aes = FALSE
    ) +
    geom_text(
      data = group_category_chisq, 
      aes(x = cell_type, y = y_max, label = significance),
      size = 6, 
      inherit.aes = FALSE
    ) +
    labs(x = "", y = "Cell Proportion", fill = "Response") +
    theme_classic() +
    theme(
      axis.text.x = if (show_x_axis) {
        element_text(angle = 45, hjust = 1, size = 12, face = "bold")
      } else {
        element_blank()
      },
      axis.text.y = element_text(size = 12, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold"),
      axis.line = element_line(color = "black"),
      panel.grid = element_blank(),
      legend.text = element_text(size = 12, face = "bold"),
      legend.title = element_text(size = 14, face = "bold"),
      legend.position = "none",
      axis.ticks.x = if (show_x_axis) element_line() else element_blank(),
      plot.margin = margin(5, 5, 5, 5)
    ) +
    # === 关键修改：使用breaks确保图例顺序 ===
    scale_fill_manual(
      values = c("YES" = "#E64B35FF", "NO" = "#4DBBD5FF"),
      breaks = c("YES", "NO")  # 明确指定图例顺序
    ) +
    scale_y_continuous(labels = scales::percent_format()) +
    scale_x_discrete(drop = FALSE)
  
  return(p)
}

# 创建所有子图
# Group1: CD4, CD8, Other (都不显示横坐标)
g1_cd4 <- create_group_category_plot("Group1", "CD4 T Cells", cd4_cells_all, FALSE)
g1_cd8 <- create_group_category_plot("Group1", "CD8 T Cells", cd8_cells_all, FALSE)
g1_other <- create_group_category_plot("Group1", "Other T Cells", other_cells_all, FALSE)

# Group2: CD4, CD8, Other (都显示横坐标)
g2_cd4 <- create_group_category_plot("Group2", "CD4 T Cells", cd4_cells_all, TRUE)
g2_cd8 <- create_group_category_plot("Group2", "CD8 T Cells", cd8_cells_all, TRUE)
g2_other <- create_group_category_plot("Group2", "Other T Cells", other_cells_all, TRUE)

# 创建图例 - 确保YES/NO顺序
legend_plot <- ggplot(plot_data_combined, aes(x = cell_type, y = proportion, fill = Response)) +
  geom_boxplot() +
  # 明确设置图例顺序
  scale_fill_manual(
    values = c("YES" = "#E64B35FF", "NO" = "#4DBBD5FF"),
    breaks = c("YES", "NO")  # YES在前，NO在后
  ) +
  theme(legend.position = "top",
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"))

legend <- get_legend(legend_plot)

# 添加类别标签
cd4_label <- ggplot() + 
  annotate("text", x = 0.5, y = 0.5, label = "CD4 T Cells", size = 5, fontface = "bold") + 
  theme_void()

cd8_label <- ggplot() + 
  annotate("text", x = 0.5, y = 0.5, label = "CD8 T Cells", size = 5, fontface = "bold") + 
  theme_void()

other_label <- ggplot() + 
  annotate("text", x = 0.5, y = 0.5, label = "Other T Cells", size = 5, fontface = "bold") + 
  theme_void()

# 组合图形
combined_plot <- (cd4_label & theme(plot.margin = margin(0,0,0,0))) / 
                 (g1_cd4 | g1_cd8 | g1_other) / 
                 (g2_cd4 | g2_cd8 | g2_other) +
  plot_layout(heights = c(0.05, 0.475, 0.475))

# 添加图例和Group标签
final_plot <- wrap_plots(
  legend,
  ggplot() + 
    annotate("text", x = 0.2, y = 0.5, label = "Group1", size = 6, fontface = "bold", hjust = 0) + 
    theme_void(),
  combined_plot,
  ggplot() + 
    annotate("text", x = 0.2, y = 0.5, label = "Group2", size = 6, fontface = "bold", hjust = 0) + 
    theme_void(),
  ncol = 1, 
  heights = c(0.05, 0.05, 0.8, 0.05)
)

# 保存图形
pdf("t_cells_CD48_combined_response_boxplots.pdf", width = 4000/300, height = 4000/300)
print(final_plot)
dev.off()

cat("已完成合并图形\n")


```