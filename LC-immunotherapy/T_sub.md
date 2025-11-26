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
pdf("t_group_response_lineplot.pdf", width = 8000/300, height = 4000/300)

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

#柱状图
library(tidyr)
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

# 绘制分面柱状图
pdf("t_cells_group_proportion.pdf", width = 6000/300, height = 4000/300)

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
    axis.text.y = element_text(size = 20, color = "black", face = "bold"),
    axis.title.y = element_text(size = 24, face = "bold", color = "black", margin = margin(r = 20)),
    
    # 图例位置
    legend.position = "top",
    legend.direction = "horizontal",
    legend.justification = "center",
    legend.text = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.box = "horizontal",
    legend.key.size = unit(0.8, "cm"),
    
    # 分面标题设置
    strip.background = element_blank(),
    strip.text = element_text(size = 20, face = "bold", vjust = 1),
    strip.placement = "outside",
    
    # 其他设置
    panel.spacing = unit(1.5, "lines")
  ) +
  facet_wrap(~ Group, ncol = 1, scales = "free_x")

print(p)
dev.off()

#箱图
library(ggpubr)
library(dplyr)
library(ggsci)

cell_counts <- t_cells@meta.data %>%
  group_by(patients, cell_type, Group) %>%
  summarise(count = n()) %>%
  ungroup()

# 计算每个样本的总细胞数
sample_totals <- cell_counts %>%
  group_by(patients) %>%
  summarise(total = sum(count)) %>%
  ungroup()

# 计算细胞类型比例
cell_proportions <- t_cells@meta.data %>%
  group_by(patients, cell_type, Group) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  left_join(sample_totals, by = "patients") %>%
  mutate(proportion = count / total)

# === 修改：使用卡方检验比较两个Group的比例差异 ===
chisq_results <- cell_proportions %>%
  group_by(cell_type) %>%
  summarise(
    p_value = {
      # 构建每个细胞类型在两个Group中的观察值
      obs_counts <- c(
        sum(count[Group == "Group1"]),
        sum(count[Group == "Group2"])
      )
      
      # 计算两个Group的总细胞数（作为期望值的基础）
      total_counts <- c(
        sum(total[Group == "Group1"]),
        sum(total[Group == "Group2"])
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

# === 修改：按首字母顺序排序 ===
celltype_alphabetical <- cell_proportions %>%
  distinct(cell_type) %>%
  arrange(cell_type) %>%
  pull(cell_type)

cell_proportions$cell_type <- factor(cell_proportions$cell_type, 
                                    levels = celltype_alphabetical)

if ("cell_type" %in% colnames(chisq_results)) {
  chisq_results$cell_type <- factor(chisq_results$cell_type, 
                                   levels = celltype_alphabetical)
}

# 绘制图形
pdf("t_cells_箱图_proportion_by_group.pdf", width = 6000/300, height = 3000/300)
ggplot(cell_proportions, aes(x = cell_type, y = proportion, fill = Group)) +
  geom_boxplot() +
  geom_text(
    data = chisq_results, 
    aes(x = cell_type, y = max(cell_proportions$proportion, na.rm = TRUE) * 1.05, 
        label = significance),
    size = 12, 
    vjust = 0.5,
    inherit.aes = FALSE
  ) +
  labs(x = "", y = "Cell Proportion", fill = "Group") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.y = element_text(size = 32),
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
pdf("t_cells_bubble_plot.pdf", width = 8000/300, height = 4000/300)

p <- ggplot(cell_proportions, aes(x = cell_type, y = group_response, 
                                 size = proportion, color = cell_type)) +
  geom_point(alpha = 0.8) +  # 使用geom_point绘制气泡
  scale_size_continuous(name = "Proportion", 
                       range = c(1, 10),  # 调整气泡大小范围
                       limits = c(0, 1),
                       breaks = c(0.2, 0.4, 0.6, 0.8),
                       labels = scales::percent_format()) +  # 设置图例刻度
  scale_color_manual(values = npg_extended) +
  scale_y_discrete(limits = rev(levels(cell_proportions$group_response))) +  # 反转y轴顺序
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
    # 移除网格线
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # 确保坐标轴框线显示
    axis.line = element_line(color = "black", linewidth = 0.5),
    panel.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    plot.title = element_text(size = 28, face = "bold", hjust = 0.5)
  ) +
  guides(color = "none")  # 移除颜色图例，只保留大小图例

print(p)
dev.off()
```
# response
```
#箱图
library(ggpubr)
library(dplyr)
library(ggsci)

# 准备数据
cell_counts <- t_cells@meta.data %>%
  group_by(patients, cell_type, Group, Response) %>%
  summarise(count = n()) %>%
  ungroup()

# 计算每个样本的总细胞数
sample_totals <- cell_counts %>%
  group_by(patients) %>%
  summarise(total = sum(count)) %>%
  ungroup()

# 计算细胞类型比例
cell_proportions <- t_cells@meta.data %>%
  group_by(patients, cell_type, Group, Response) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  left_join(sample_totals, by = "patients") %>%
  mutate(proportion = count / total)

# 按首字母顺序排序细胞类型
celltype_alphabetical <- cell_proportions %>%
  distinct(cell_type) %>%
  arrange(cell_type) %>%
  pull(cell_type)

cell_proportions$cell_type <- factor(cell_proportions$cell_type, 
                                    levels = celltype_alphabetical)

# 设置Response顺序 - YES在前，NO在后
cell_proportions$Response <- factor(cell_proportions$Response, levels = c("YES", "NO"))

# === Group1的箱线图 ===
group1_data <- cell_proportions %>% filter(Group == "Group1")

# Group1的统计检验
chisq_results_group1 <- group1_data %>%
  group_by(cell_type) %>%
  summarise(
    p_value = {
      obs_counts <- c(
        sum(count[Response == "YES"]),
        sum(count[Response == "NO"])
      )
      
      total_counts <- c(
        sum(total[Response == "YES"]),
        sum(total[Response == "NO"])
      )
      
      if (sum(obs_counts) < 10 || any(obs_counts == 0)) {
        NA_real_
      } else {
        expected_prop <- total_counts / sum(total_counts)
        chisq.test(obs_counts, p = expected_prop)$p.value
      }
    }
  ) %>%
  ungroup() %>%
  mutate(
    significance = case_when(
      is.na(p_value) ~ "ns",
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**", 
      p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

chisq_results_group1$cell_type <- factor(chisq_results_group1$cell_type, 
                                        levels = celltype_alphabetical)

# 绘制Group1的箱线图
pdf("t_cells_Group1_Response_boxplot.pdf", width = 6000/300, height = 3000/300)
ggplot(group1_data, aes(x = cell_type, y = proportion, fill = Response)) +
  geom_boxplot() +
  geom_text(
    data = chisq_results_group1, 
    aes(x = cell_type, y = max(group1_data$proportion, na.rm = TRUE) * 1.05, 
        label = significance),
    size = 12, 
    vjust = 0.5,
    inherit.aes = FALSE
  ) +
  labs(x = "", y = "Cell Proportion", fill = "Response", 
       title = "Group1: Response Comparison") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.y = element_text(size = 32),
    axis.line = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 36),
    legend.title = element_text(size = 40),
    legend.position = "right",
    plot.title = element_text(size = 32, face = "bold", hjust = 0.5)
  ) +
  scale_fill_npg() +
  scale_y_continuous(limits = c(0, 1))
dev.off()

# === Group2的箱线图 ===
group2_data <- cell_proportions %>% filter(Group == "Group2")

# Group2的统计检验
chisq_results_group2 <- group2_data %>%
  group_by(cell_type) %>%
  summarise(
    p_value = {
      obs_counts <- c(
        sum(count[Response == "YES"]),
        sum(count[Response == "NO"])
      )
      
      total_counts <- c(
        sum(total[Response == "YES"]),
        sum(total[Response == "NO"])
      )
      
      if (sum(obs_counts) < 10 || any(obs_counts == 0)) {
        NA_real_
      } else {
        expected_prop <- total_counts / sum(total_counts)
        chisq.test(obs_counts, p = expected_prop)$p.value
      }
    }
  ) %>%
  ungroup() %>%
  mutate(
    significance = case_when(
      is.na(p_value) ~ "ns",
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**", 
      p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

chisq_results_group2$cell_type <- factor(chisq_results_group2$cell_type, 
                                        levels = celltype_alphabetical)

# 绘制Group2的箱线图
pdf("t_cells_Group2_Response_boxplot.pdf", width = 6000/300, height = 3000/300)
ggplot(group2_data, aes(x = cell_type, y = proportion, fill = Response)) +
  geom_boxplot() +
  geom_text(
    data = chisq_results_group2, 
    aes(x = cell_type, y = max(group2_data$proportion, na.rm = TRUE) * 1.05, 
        label = significance),
    size = 12, 
    vjust = 0.5,
    inherit.aes = FALSE
  ) +
  labs(x = "", y = "Cell Proportion", fill = "Response",
       title = "Group2: Response Comparison") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.y = element_text(size = 32),
    axis.line = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 36),
    legend.title = element_text(size = 40),
    legend.position = "right",
    plot.title = element_text(size = 32, face = "bold", hjust = 0.5)
  ) +
  scale_fill_npg() +
  scale_y_continuous(limits = c(0, 1))
dev.off()


#桑葚图
# 准备包含细胞类型的数据
sankey_data_celltype <- t_cells@meta.data %>%
  group_by(Group, Response, cell_type) %>%
  summarise(cell_count = n(), .groups = 'drop') %>%
  group_by(Group, Response) %>%
  mutate(total = sum(cell_count),
         proportion = cell_count / total) %>%
  ungroup()

# 三层桑葚图
pdf("t_cells_sankey_3layer.pdf", width = 8000/300, height = 5000/300)

ggplot(sankey_data_celltype,
       aes(axis1 = Group, axis2 = Response, axis3 = cell_type, 
           y = cell_count)) +
  geom_alluvium(aes(fill = cell_type), width = 1/8) +
  geom_stratum(width = 1/8, fill = "grey", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), 
            size = 6, fontface = "bold") +
  scale_x_discrete(limits = c("Group", "Response", "Cell Type"), 
                   expand = c(0.05, 0.05)) +
  scale_fill_manual(values = npg_extended) +
  labs(title = "Cell Distribution: Group → Response → Cell Type",
       y = "Number of Cells") +
  theme_void() +
  theme(
    plot.title = element_text(size = 28, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18),
    legend.position = "right"
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