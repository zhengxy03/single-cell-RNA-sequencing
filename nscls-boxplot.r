
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
    
    ct_plot_data <- ct_full_data %>%
        filter(.data[[target_gene]] > 0)
    
    if(nrow(ct_plot_data) == 0) {
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
    
    y_max <- max(ct_plot_data[[target_gene]], na.rm = TRUE) * 1.05
    

    p <- ggplot(ct_plot_data, aes(x = condition, y = .data[[target_gene]], fill = condition)) +
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


combined_plot <- wrap_plots(cell_plots, ncol = 3) + 
    plot_annotation(
        title = paste(target_gene, "expression by cell type"),
        theme = theme(
            plot.title = element_text(hjust = 0, size = 16, face = "bold")  
        )
    )


combined_plot <- combined_plot + plot_annotation(
    theme = theme(legend.position = "bottom")
) & guides(fill = guide_legend(title = "Condition"))


print(combined_plot)    