
target_gene <- "CDKN2A"

full_data <- FetchData(
    seurat_obj,
    vars = c(target_gene, "cell_type", "condition"),
    layer = "data"
)

cell_types <- unique(full_data$cell_type)
conditions <- unique(full_data$condition)

cell_plots <- lapply(cell_types, function(ct) {
    ct_full_data <- full_data %>%
        filter(cell_type == ct)

    if(nrow(ct_full_data) == 0) {
        p <- ggplot() + 
            geom_blank() +
            labs(title = ct) +
            theme_void() +
            theme(plot.title = element_text(hjust = 0.5, size = 12))
        return(p)
    }
    
    group1_full <- ct_full_data[ct_full_data$condition == conditions[1], target_gene]
    group2_full <- ct_full_data[ct_full_data$condition == conditions[2], target_gene]

    is_expressed_full <- any(group1_full > 0) || any(group2_full > 0)

    if(is_expressed_full && all(table(ct_full_data$condition) > 0)) {
        test_result <- wilcox.test(group1_full, group2_full)
        p_val <- formatC(test_result$p.value, format = "e", digit = 2)
        p_text <- paste("p =", p_val)
    } else {
        p_text <- "Insufficient data"
    }

    y_max <- max(ct_full_data[[target_gene]], na.rm = TRUE) * 1.05

    p <- ggplot(ct_full_data, aes(x = condition, y = .data[[target_gene]], fill = condition)) +
        geom_boxplot(width = 0.6, outlier.size = 2, alpha = 0.7) +
        labs(title = ct) +
        theme_classic() +
        theme(
            legend.position = "none",
            plot.title = element_text(hjust = 0.5, size = 12),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 10)
        ) +
        ylim(0, y_max)
    

    p <- p + geom_text(
        data = data.frame(),
        aes(x = 1.5, y = y_max * 0.95, label = p_text),
        inherit.aes = FALSE,
        size = 3
    )

    if(ct == cell_types[1]) {
        p <- p + labs(y = expression(log[2]("CPM/100 + 1")))
    } else {
        p <- p + labs(y = "")
    }
    
    return(p)
})

combined_plot <- wrap_plots(cell_plots, ncol = 3) +  # 每行显示3个子图
    plot_annotation(
        title = paste(target_gene, "expression by cell type"),
        theme = theme(
            plot.title = element_text(hjust = 0, size = 16, face = "bold")  # 大标题显示在左上角
        )
    )

combined_plot <- combined_plot + plot_annotation(
    theme = theme(legend.position = "bottom")
) & guides(fill = guide_legend(title = "Condition"))

print(combined_plot)




#隐藏离群值

target_gene <- "CDKN2A"

full_data <- FetchData(
    seurat_obj,
    vars = c(target_gene, "cell_type", "condition"),
    layer = "data"
)

cell_types <- unique(full_data$cell_type)
conditions <- unique(full_data$condition)

# 定义函数识别离群点（基于Tukey方法）
identify_outliers <- function(x) {
  q1 <- quantile(x, 0.25, na.rm = TRUE)
  q3 <- quantile(x, 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  lower <- q1 - 1.5 * iqr
  upper <- q3 + 1.5 * iqr
  return(x < lower | x > upper)
}


cell_plots <- lapply(cell_types, function(ct) {

    ct_full_data <- full_data %>%
        filter(cell_type == ct)

    if(nrow(ct_full_data) == 0) {
        p <- ggplot() + 
            geom_blank() +
            labs(title = ct) +
            theme_void() +
            theme(plot.title = element_text(hjust = 0.5, size = 12))
        return(p)
    }

    group1_full <- ct_full_data[ct_full_data$condition == conditions[1], target_gene]
    group2_full <- ct_full_data[ct_full_data$condition == conditions[2], target_gene]
    
    is_expressed_full <- any(group1_full > 0) || any(group2_full > 0)
    
    if(is_expressed_full && all(table(ct_full_data$condition) > 0)) {
        test_result <- wilcox.test(group1_full, group2_full)
        p_val <- formatC(test_result$p.value, format = "e", digit = 2)
        p_text <- paste("p =", p_val)
    } else {
        p_text <- "Insufficient data"
    }
    
    filtered_data <- ct_full_data %>%
        group_by(condition) %>%
        mutate(is_outlier = identify_outliers(.data[[target_gene]])) %>%
        ungroup()
    

    y_max <- max(filtered_data[!filtered_data$is_outlier,][[target_gene]], na.rm = TRUE) * 1.05
    

    p <- ggplot(ct_full_data, aes(x = condition, y = .data[[target_gene]], fill = condition)) +
        geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.7) + 
        labs(title = ct) + 
        theme_classic() +
        theme(
            legend.position = "none", 
            plot.title = element_text(hjust = 0.5, size = 12),
            axis.text.x = element_blank(), 
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 10)
        ) +
        ylim(0, y_max)
    

    p <- p + geom_text(
        data = data.frame(),
        aes(x = 1.5, y = y_max * 0.95, label = p_text),
        inherit.aes = FALSE,
        size = 3
    )
    

    if(ct == cell_types[1]) {
        p <- p + labs(y = expression(log[2]("CPM/100 + 1")))
    } else {
        p <- p + labs(y = "")
    }
    
    return(p)
})


combined_plot <- wrap_plots(cell_plots, ncol = 3) +  # 每行显示3个子图
    plot_annotation(
        title = paste(target_gene, "expression by cell type"),
        theme = theme(
            plot.title = element_text(hjust = 0, size = 16, face = "bold")  # 大标题显示在左上角
        )
    )

combined_plot <- combined_plot + plot_annotation(
    theme = theme(legend.position = "bottom")
) & guides(fill = guide_legend(title = "Condition"))

print(combined_plot)


#barplot

target_gene <- "CDKN2A"
target_cell_type <- "Effector T cell"

full_data <- FetchData(
    seurat_obj,
    vars = c(target_gene, "cell_type", "condition"),
    layer = "data"
) %>%
    filter(cell_type == target_cell_type)

conditions <- unique(full_data$condition)
group1 <- full_data[full_data$condition == conditions[1], target_gene]
group2 <- full_data[full_data$condition == conditions[2], target_gene]


is_expressed <- any(group1 > 0) || any(group2 > 0)


if (is_expressed && all(table(full_data$condition) > 0)) {
    test_result <- wilcox.test(group1, group2)
    p_val <- formatC(test_result$p.value, format = "e", digit = 2)
    p_text <- paste("p =", p_val)
} else {
    p_text <- "Insufficient data"
}


summary_data <- full_data %>%
    group_by(condition) %>%
    summarise(
        mean_expr = mean(.data[[target_gene]], na.rm = TRUE),
        sd_expr = sd(.data[[target_gene]], na.rm = TRUE)
    )


p <- ggplot(summary_data, aes(x = condition, y = mean_expr, fill = condition)) +
    geom_bar(stat = "identity", width = 0.7, alpha = 0.7) +

    geom_errorbar(aes(ymin = mean_expr - sd_expr, ymax = mean_expr + sd_expr),
                  width = 0.3, color = "black", size = 0.8) +

    labs(
        title = paste(target_gene, "in", target_cell_type),
        subtitle = p_text,
        y = expression(log[2]("CPM/100 + 1"))
    ) +
    theme_classic() +
    theme(
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, size = 12),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()
    )


print(p)



target_gene <- "CDKN2A"

full_data <- FetchData(
    seurat_obj,
    vars = c(target_gene, "cell_type", "condition"),
    layer = "data"
)


cell_types <- unique(full_data$cell_type)

cell_plots <- list()
for (ct in cell_types) {
    ct_full_data <- full_data %>%
        filter(cell_type == ct)

    if (nrow(ct_full_data) == 0) {
        warning(paste("没有找到", ct, "的数据"))
        next
    }
    
    conditions <- unique(ct_full_data$condition)
    group1 <- ct_full_data[ct_full_data$condition == conditions[1], target_gene]
    group2 <- ct_full_data[ct_full_data$condition == conditions[2], target_gene]
    

    is_expressed <- any(group1 > 0) || any(group2 > 0)
    
    if (is_expressed && all(table(ct_full_data$condition) > 0)) {
        test_result <- wilcox.test(group1, group2)
        p_val <- formatC(test_result$p.value, format = "e", digit = 2)
        p_text <- paste("p =", p_val)
    } else {
        p_text <- "Insufficient data"
    }
    
    summary_data <- ct_full_data %>%
        group_by(condition) %>%
        summarise(
            mean_expr = mean(.data[[target_gene]], na.rm = TRUE),
            sd_expr = sd(.data[[target_gene]], na.rm = TRUE)
        )
    
    p <- ggplot(summary_data, aes(x = condition, y = mean_expr, fill = condition)) +
        geom_bar(stat = "identity", width = 0.7, alpha = 0.7) +
        geom_segment(
            aes(
                x = as.numeric(factor(condition)),  # 转换为数值位置
                xend = as.numeric(factor(condition)),
                y = mean_expr,
                yend = mean_expr + sd_expr
                ),
        color = "black",
        size = 0.8
        ) +
  
  # 添加顶部短横线
        geom_segment(
                    aes(
                        x = as.numeric(factor(condition)) - 0.15,
                        xend = as.numeric(factor(condition)) + 0.15,
                        y = mean_expr + sd_expr,
                        yend = mean_expr + sd_expr
                        ),
                    color = "black",
                    size = 0.8
                    ) +
        
        labs(
            title = ct,
            subtitle = p_text,
            y = expression(log[2]("CPM/100 + 1"))
        ) +
        theme_classic() +
        theme(
            legend.position = "none",
            plot.title = element_text(hjust = 0.5, size = 12),

            axis.title.x = element_blank()
        )

    cell_plots[[ct]] <- p
}



combined_plot <- wrap_plots(cell_plots, ncol = 7) +
    plot_annotation(title = paste(target_gene, "expression by cell type"),
                    theme = theme(plot.title = element_text(hjust = 0, size = 16, face = "bold")))
print(combined_plot)


#显著性分类
target_gene <- "CDKN2A"

target_gene <- "WNT7B"

full_data <- FetchData(
    seurat_obj,
    vars = c(target_gene, "cell_type", "condition"),
    layer = "data"
)

cell_types <- unique(full_data$cell_type)


ns_plots <- list()       # 不显著
up_plots <- list()       # 显著升高
down_plots <- list()     # 显著降低


p_threshold <- 0.05

for (ct in cell_types) {
    ct_full_data <- full_data %>%
        filter(cell_type == ct)
    
    if (nrow(ct_full_data) == 0) {
        warning(paste("没有找到", ct, "的数据"))
        next
    }
    
    conditions <- unique(ct_full_data$condition)
    
    
    if (length(conditions) != 2) {
        warning(paste(ct, "的条件数量不为2，跳过"))
        next
    }
    
    
    if ("Tumor" %in% conditions && "Normal" %in% conditions) {
        tumor_group <- ct_full_data[ct_full_data$condition == "Tumor", target_gene]
        normal_group <- ct_full_data[ct_full_data$condition == "Normal", target_gene]
        
        # 检查每组至少有2个样本（Wilcoxon检验的最低要求）
        if (length(tumor_group) < 2 || length(normal_group) < 2) {
            warning(paste(ct, "的样本量不足，无法进行Wilcoxon检验"))
            p_text <- "Insufficient samples"
            group <- "ns"
        } else {
            is_expressed <- any(tumor_group > 0) || any(normal_group > 0)
            
            if (is_expressed) {
                # 执行Wilcoxon检验并处理可能的错误
                tryCatch({
                    test_result <- wilcox.test(tumor_group, normal_group)
                    p_val <- test_result$p.value
                    
                    # 检查p值是否为NA
                    if (is.na(p_val)) {
                        warning(paste(ct, "的Wilcoxon检验无法计算p值"))
                        p_text <- "p-value unavailable"
                        group <- "ns"
                    } else {
                        p_text <- paste("p =", formatC(p_val, format = "e", digits = 2))
                        
                        # 计算均值差异判断表达趋势
                        mean_tumor <- mean(tumor_group, na.rm = TRUE)
                        mean_normal <- mean(normal_group, na.rm = TRUE)
                        fold_change <- mean_tumor / mean_normal
                        
                        # 确定分组
                        if (p_val >= p_threshold) {
                            group <- "ns"  # 不显著
                        } else if (fold_change > 1) {
                            group <- "up"  # 显著升高
                        } else {
                            group <- "down"  # 显著降低
                        }
                    }
                }, error = function(e) {
                    warning(paste(ct, "的Wilcoxon检验出错:", e$message))
                    p_text <- "Test error"
                    group <- "ns"
                })
            } else {
                p_text <- "No expression"
                group <- "ns"
            }
        }
    } else {
        warning(paste("未找到Tumor或Normal条件，跳过", ct))
        next
    }
    
    summary_data <- ct_full_data %>%
        group_by(condition) %>%
        summarise(
            mean_expr = mean(.data[[target_gene]], na.rm = TRUE),
            sd_expr = sd(.data[[target_gene]], na.rm = TRUE)
        )
    
    p <- ggplot(summary_data, aes(x = condition, y = mean_expr, fill = condition)) +
        geom_bar(stat = "identity", width = 0.7, alpha = 0.7) +
        geom_segment(
            aes(
                x = as.numeric(factor(condition)),  # 转换为数值位置
                xend = as.numeric(factor(condition)),
                y = mean_expr,
                yend = mean_expr + sd_expr
            ),
            color = "black",
            size = 0.8
        ) +
        
        # 添加顶部短横线
        geom_segment(
            aes(
                x = as.numeric(factor(condition)) - 0.15,
                xend = as.numeric(factor(condition)) + 0.15,
                y = mean_expr + sd_expr,
                yend = mean_expr + sd_expr
            ),
            color = "black",
            size = 0.8
        ) +
        labs(
            title = ct,
            subtitle = p_text,
            y = expression(log[2]("CPM/100 + 1"))
        ) +
        scale_x_discrete(limits = c("Tumor", "Normal")) +
        theme_classic() +
        theme(
            legend.position = "none",
            plot.title = element_text(hjust = 0.5, size = 12),
            axis.title.x = element_blank()
        )
    
    
    if (group == "ns") {
        ns_plots[[ct]] <- p
    } else if (group == "up") {
        up_plots[[ct]] <- p
    } else if (group == "down") {
        down_plots[[ct]] <- p
    }
}


if (length(ns_plots) > 0) {
    ns_combined <- wrap_plots(ns_plots, ncol = 7) +
        plot_annotation(title = paste(target_gene, "Expression: Not Significant"),
                        theme = theme(plot.title = element_text(hjust = 0, size = 16, face = "bold")))
    print(ns_combined)
}

if (length(up_plots) > 0) {
    up_combined <- wrap_plots(up_plots, ncol = 7) +
        plot_annotation(title = paste(target_gene, "Expression: Significantly Upregulated in Tumor"),
                        theme = theme(plot.title = element_text(hjust = 0, size = 16, face = "bold")))
    print(up_combined)
}

if (length(down_plots) > 0) {
    down_combined <- wrap_plots(down_plots, ncol = 7) +
        plot_annotation(title = paste(target_gene, "Expression: Significantly Downregulated in Tumor"),
                        theme = theme(plot.title = element_text(hjust = 0, size = 16, face = "bold")))
    print(down_combined)
}


calc_grid <- function(n, ncol = 7) {
    nrow <- ceiling(n / ncol)
    return(list(nrow = nrow, ncol = ncol))
}


create_panel <- function(plot_list, title, ncol = 7) {
    if (length(plot_list) == 0) return(NULL)
    
    grid <- calc_grid(length(plot_list), ncol)
    
    
    title_plot <- ggplot() + 
        annotate("text", x = 0.5, y = 0.5, label = title, 
                 size = 5, fontface = "bold", hjust = 0.5) +
        theme_void() +
        theme(plot.margin = margin(5, 0, 5, 5))
    
    # 组合标题和子图
    wrap_plots(
        title_plot,
        wrap_plots(plot_list, ncol = ncol),
        ncol = 1,
        heights = c(0.5, grid$nrow)
    )
}

# 创建所有面板并组合
create_combined_plot <- function(up_plots, down_plots, ns_plots, target_gene, ncol = 7) {
    # 生成带标题的面板
    panel_up <- create_panel(up_plots, "UPREGULATED IN TUMOR", ncol)
    panel_down <- create_panel(down_plots, "DOWNREGULATED IN TUMOR", ncol)
    panel_ns <- create_panel(ns_plots, "NOT SIGNIFICANT", ncol)
    
    # 计算总高度
    grid_up <- calc_grid(length(up_plots), ncol)
    grid_down <- calc_grid(length(down_plots), ncol)
    grid_ns <- calc_grid(length(ns_plots), ncol)
    
    heights <- c(
        ifelse(!is.null(panel_up), 0.5 + grid_up$nrow, 0),
        ifelse(!is.null(panel_down), 0.5 + grid_down$nrow, 0),
        ifelse(!is.null(panel_ns), 0.5 + grid_ns$nrow, 0)
    )
    
    # 垂直拼接所有面板
    combined <- wrap_plots(
        list(panel_up, panel_down, panel_ns),
        ncol = 1,
        heights = heights
    ) + 
        plot_annotation(
            title = paste("DIFFERENTIAL EXPRESSION OF", target_gene),
            theme = theme(
                plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(b = 10))
            )
        )
    
    # 计算合适的绘图设备尺寸
    base_width <- 14
    base_height <- 2 * (sum(heights) - 0.5)  # 减去标题的0.5
    
    return(list(
        plot = combined,
        width = base_width,
        height = base_height
    ))
}

# 使用函数创建并保存图片
result <- create_combined_plot(up_plots, down_plots, ns_plots, target_gene)

if (!is.null(result$plot)) {
    print(result$plot)
    ggsave(paste0(target_gene, "_combined_expression.png"), 
           plot = result$plot, 
           width = result$width, 
           height = result$height, 
           dpi = 300)
}


