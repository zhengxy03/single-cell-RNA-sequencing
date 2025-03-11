#sample_type
# 获取唯一的 sample_type 值，并按字母顺序排序
unique_sample_types <- sort(unique(epi@meta.data$sample_type))  # 按字母顺序排序

# 如果需要手动指定顺序，可以这样做：
# unique_sample_types <- c("Tumor", "Normal", "Blood")

# 创建一个空列表，用于存储每个 sample_type 的图
plot_list <- list()

# 遍历每个 sample_type，生成单独的图
for (sample in unique_sample_types) {
    # 创建颜色映射：当前 sample_type 为 npg 红色，其他为灰色
    color_mapping <- setNames(
        ifelse(unique_sample_types == sample, pal_npg()(1), "gray"),  # npg 红色
        unique_sample_types
    )
    
    # 绘制 DimPlot
    p <- DimPlot(epi, reduction = "umap", label = FALSE, pt.size = 1, group.by = "sample_type") +
        scale_color_manual(values = color_mapping) +  # 使用自定义颜色映射
        ggtitle(sample) +  # 设置标题为当前 sample_type
        theme(
            legend.position = "none",  # 隐藏图例
            plot.title = element_text(size = 24, face = "bold", hjust = 0.5),  # 标题样式
            axis.title.x = element_text(size = 20, face = "bold", color = "black"),  # X 轴标题
            axis.title.y = element_text(size = 20, face = "bold", color = "black"),  # Y 轴标题
            axis.text.x = element_text(size = 16, color = "black"),  # X 轴刻度
            axis.text.y = element_text(size = 16, color = "black")   # Y 轴刻度
        )
    
    # 将图添加到列表中
    plot_list[[sample]] <- p
}

# 使用 patchwork 将图形排列在一起（每行两个）
combined_plot <- wrap_plots(plot_list, ncol = 2)  

# 保存为 PNG 文件
png("epi_sample_type_umap.png", width = 6000, height = 3000, res = 300)
print(combined_plot)
dev.off()
