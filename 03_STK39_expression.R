# ============================================================
# Part 3. 
# ============================================================

setwd('/work/home/20260401_STK39')

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(ggpubr)

# -------------------- 1.  --------------------
group1_samples <- c(
  "GSM8037607", "GSM8037621", "GSM8037610",
  "GSM8037622", "GSM8037623", "GSM8037608",
  "GSM8037611", "GSM8037624", "GSM8037625",
  "GSM8037619"
)

group2_samples <- c(
  "GSM8037617", "GSM8037618",
  "GSM8037630", "GSM8037628"
)

group1_name <- "Young"
group2_name <- "Old"

# -------------------- 2.  --------------------
integrated <- readRDS(
  "TESTIS_FINAL_ANNOTATED.rds"
)

if (!"STK39" %in% rownames(integrated)) {
  stop("未检测到")
}

# -------------------- 3.  --------------------
if ("RNA" %in% names(integrated@assays)) {
  DefaultAssay(integrated) <- "RNA"
} else if ("SCT" %in% names(integrated@assays)) {
  DefaultAssay(integrated) <- "SCT"
} else {
  DefaultAssay(integrated) <- names(integrated@assays)[1]
}

stk39_expr <- FetchData(
  integrated,
  vars = "STK39"
)

integrated$STK39_expr <- stk39_expr$STK39

expr_data <- data.frame(
  cell_type = integrated$cell_type,
  sample = integrated$orig.ident,
  STK39 = integrated$STK39_expr,
  stringsAsFactors = FALSE
)

expr_data <- expr_data[
  !grepl("Unknown", expr_data$cell_type),
]

# -------------------- 4.  --------------------
expr_data$group <- ifelse(
  expr_data$sample %in% group1_samples,
  group1_name,
  ifelse(
    expr_data$sample %in% group2_samples,
    group2_name,
    "Other"
  )
)

expr_data <- expr_data[
  expr_data$group %in% c(
    group1_name,
    group2_name
  ),
]

expr_data$group <- factor(
  expr_data$group,
  levels = c(
    group1_name,
    group2_name
  )
)

integrated$group <- ifelse(
  integrated$orig.ident %in% group1_samples,
  group1_name,
  ifelse(
    integrated$orig.ident %in% group2_samples,
    group2_name,
    "Other"
  )
)

integrated$group <- factor(
  integrated$group,
  levels = c(
    group1_name,
    group2_name,
    "Other"
  )
)

# -------------------- 5.  --------------------
cell_types <- unique(
  expr_data$cell_type
)

group_colors <- c(
  "#1f77b4",
  "#ff7f0e"
)

stat_results_list <- list()
violin_plot_list <- list()

for (ct in cell_types) {
  ct_data <- expr_data[
    expr_data$cell_type == ct,
  ]

  if (nrow(ct_data) == 0) {
    next
  }

  ct_data$group <- factor(
    ct_data$group,
    levels = c(
      group1_name,
      group2_name
    )
  )

  if (length(unique(ct_data$group)) != 2) {
    next
  }

  group1_data <- ct_data$STK39[
    ct_data$group == group1_name
  ]
  group2_data <- ct_data$STK39[
    ct_data$group == group2_name
  ]

  wilcox_test <- wilcox.test(
    group1_data,
    group2_data,
    exact = FALSE
  )

  mean_diff <- mean(
    group1_data,
    na.rm = TRUE
  ) - mean(
    group2_data,
    na.rm = TRUE
  )

  stat_results_list[[ct]] <- data.frame(
    cell_type = ct,
    test = "Wilcoxon rank-sum test",
    statistic = unname(
      wilcox_test$statistic
    ),
    p_value = wilcox_test$p.value,
    mean_group1 = mean(
      group1_data,
      na.rm = TRUE
    ),
    mean_group2 = mean(
      group2_data,
      na.rm = TRUE
    ),
    mean_difference = mean_diff,
    n_cells_group1 = length(group1_data),
    n_cells_group2 = length(group2_data),
    n_total = nrow(ct_data)
  )

  group_means <- ct_data %>%
    group_by(group) %>%
    summarise(
      mean_expr = mean(
        STK39,
        na.rm = TRUE
      ),
      .groups = "drop"
    )

  y_range <- range(
    ct_data$STK39,
    na.rm = TRUE
  )

  y_span <- y_range[2] - y_range[1]
  if (!is.finite(y_span) || y_span == 0) {
    y_span <- max(
      abs(y_range),
      1
    ) * 0.15
  }

  y_padding <- y_span * 0.15
  y_lims <- c(
    y_range[1] - y_padding,
    y_range[2] + y_padding
  )

  p <- ggplot(
    ct_data,
    aes(
      x = group,
      y = STK39,
      fill = group
    )
  ) +
    geom_violin(
      scale = "width",
      alpha = 0.7,
      trim = TRUE
    ) +
    geom_point(
      data = group_means,
      aes(
        x = group,
        y = mean_expr,
        fill = group
      ),
      shape = 21,
      size = 4,
      color = "black",
      stroke = 1
    ) +
    stat_compare_means(
      method = "wilcox.test",
      comparisons = list(
        c(
          group1_name,
          group2_name
        )
      ),
      label = "p.format",
      label.x.npc = "center",
      label.y = max(y_lims) * 0.95,
      size = 4,
      vjust = 1
    ) +
    stat_compare_means(
      method = "wilcox.test",
      comparisons = list(
        c(
          group1_name,
          group2_name
        )
      ),
      label = "p.signif",
      hide.ns = TRUE,
      tip.length = 0.01,
      step.increase = 0.08,
      size = 6,
      vjust = 0.5
    ) +
    scale_fill_manual(
      values = group_colors
    ) +
    scale_y_continuous(
      limits = y_lims,
      expand = expansion(
        mult = c(0.05, 0.15)
      )
    ) +
    labs(
      title = paste(
        "STK39 Expression in",
        ct
      ),
      subtitle = paste0(
        "Young: ",
        sum(ct_data$group == group1_name),
        " cells | Old: ",
        sum(ct_data$group == group2_name),
        " cells\nPoints show group means"
      ),
      x = "Group",
      y = "STK39 Expression",
      fill = "Group"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        size = 14,
        face = "bold"
      ),
      plot.subtitle = element_text(
        hjust = 0.5,
        size = 11
      ),
      axis.text.x = element_text(
        size = 12
      ),
      axis.text.y = element_text(
        size = 12
      ),
      axis.title = element_text(
        size = 13
      ),
      legend.position = "right"
    )

  violin_plot_list[[ct]] <- p

  safe_ct_name <- gsub(
    "[^a-zA-Z0-9]",
    "_",
    ct
  )

  ggsave(
    paste0(
      "STK39_EXPRESSION_VIOLIN_",
      safe_ct_name,
      "_TWO_GROUPS.pdf"
    ),
    p,
    width = 8,
    height = 8
  )

  ggsave(
    paste0(
      "STK39_EXPRESSION_VIOLIN_",
      safe_ct_name,
      "_TWO_GROUPS.png"
    ),
    p,
    width = 8,
    height = 8,
    dpi = 150
  )
}

# -------------------- 6.  --------------------
if (length(stat_results_list) > 0) {
  stat_results_df <- do.call(
    rbind,
    stat_results_list
  )

  stat_results_df$significance <- ifelse(
    stat_results_df$p_value < 0.001,
    "***",
    ifelse(
      stat_results_df$p_value < 0.01,
      "**",
      ifelse(
        stat_results_df$p_value < 0.05,
        "*",
        "ns"
      )
    )
  )

  write.csv(
    stat_results_df,
    "STK39_CELLTYPE_GROUP_COMPARISONS_WITH_SIGNIFICANCE.csv",
    row.names = FALSE
  )
}

# -------------------- 7.  --------------------
if (length(violin_plot_list) > 0) {
  n_plots <- length(
    violin_plot_list
  )
  n_cols <- 2
  n_rows <- ceiling(
    n_plots / n_cols
  )

  violin_for_combined <- lapply(
    violin_plot_list,
    function(p) {
      p + theme(
        legend.position = "none"
      )
    }
  )

  combined_plot <- wrap_plots(
    violin_for_combined,
    ncol = n_cols,
    nrow = n_rows
  ) +
    plot_annotation(
      title = paste(
        "STK39 Expression:",
        group1_name,
        "vs",
        group2_name
      ),
      subtitle = paste(
        "Wilcoxon rank-sum test |",
        n_plots,
        "cell types"
      ),
      tag_levels = "A"
    )

  ggsave(
    "STK39_EXPRESSION_VIOLIN_COMBINED_TWO_GROUPS_2COL.pdf",
    combined_plot,
    width = 14,
    height = n_rows * 7
  )

  ggsave(
    "STK39_EXPRESSION_VIOLIN_COMBINED_TWO_GROUPS_2COL.png",
    combined_plot,
    width = 14,
    height = n_rows * 7,
    dpi = 150
  )
}

# -------------------- 8.  --------------------
p_stk39_umap <- FeaturePlot(
  integrated,
  features = "STK39",
  reduction = "umap",
  min.cutoff = "q10",
  max.cutoff = "q90",
  pt.size = 0.5,
  order = TRUE,
  cols = c(
    "lightgrey",
    "red"
  )
) +
  ggtitle(
    "STK39 Expression on UMAP"
  ) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      size = 16,
      face = "bold"
    )
  )

ggsave(
  "STK39_UMAP_EXPRESSION.pdf",
  p_stk39_umap,
  width = 10,
  height = 8
)

ggsave(
  "STK39_UMAP_EXPRESSION.png",
  p_stk39_umap,
  width = 10,
  height = 8,
  dpi = 300
)

# -------------------- 9.  --------------------
two_group_subset <- subset(
  integrated,
  subset =
    orig.ident %in%
    c(
      group1_samples,
      group2_samples
    )
)

p_stk39_umap_two_groups <- FeaturePlot(
  two_group_subset,
  features = "STK39",
  reduction = "umap",
  min.cutoff = "q10",
  max.cutoff = "q90",
  pt.size = 0.5,
  order = TRUE,
  cols = c(
    "lightgrey",
    "red"
  )
) +
  ggtitle(
    "STK39 Expression on UMAP: Young + Old Samples"
  ) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      size = 16,
      face = "bold"
    )
  )

ggsave(
  "STK39_UMAP_EXPRESSION_YOUNG_OLD.pdf",
  p_stk39_umap_two_groups,
  width = 10,
  height = 8
)

ggsave(
  "STK39_UMAP_EXPRESSION_YOUNG_OLD.png",
  p_stk39_umap_two_groups,
  width = 10,
  height = 8,
  dpi = 300
)

# -------------------- 10.--------------------
saveRDS(
  integrated,
  "TESTIS_FINAL_ANNOTATED_WITH_STK39_GROUP.rds"
)
