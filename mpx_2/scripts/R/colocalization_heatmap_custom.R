
ColocalizationHeatmap_custom <- function (
    data, 
    marker1_col = "marker_1", 
    marker2_col = "marker_2", 
    value_col = "estimate", 
    size_col = "p_adj", 
    size_col_transform = NULL, 
    size_range = c(0.1, 3), 
    colors = c("#053061", "#2166AC","#4393C3", "#92C5DE", "#D1E5F0", 
               "#F7F7F7", 
               "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"), 
    cluster_rows = TRUE, 
    cluster_cols = TRUE, 
    clustering_distance_rows = "euclidean", 
    clustering_distance_cols = "euclidean", 
    clustering_method = "complete", 
    type = c("tiles", "dots"), 
    return_plot_data = FALSE, 
    symmetrise = TRUE, 
    marker_order = NULL, # NEW
    legend_range = NULL, 
    legend_title = "", ...) {
  
  if (any(size_range[1] < 0)) {
    cli::cli_abort(c(i = "{.arg size_range} must have non-negative values",
                     x = "Found {.val {sum(size_range < 0)}} negative values"))
  }
  
  if (is.grouped_df(data)) {
    cli::cli_warn(c(i = "The input data is grouped and will be ungrouped."))
    data <- data %>% ungroup()
  }
  
  cols_keep <- c(marker1_col, marker2_col, value_col)
  numeric_cols <- value_col
  
  if (type == "dots") {
    # assert_col_in_data(size_col, data)
    cols_keep <- c(cols_keep, size_col)
    numeric_cols <- c(numeric_cols, size_col)
  }
  
  plot_data <- data %>% select(all_of(cols_keep))
  
  dup_rows <- duplicated(select(plot_data, all_of(c(marker1_col, marker2_col))))
  
  if (sum(dup_rows) > 0) {
    dup_pairs <- select(plot_data, all_of(c(marker1_col, 
                                            marker2_col)
    )
    ) %>% 
      apply(1, paste, collapse = "/")
    cli::cli_abort(c(i = "Each row must represent a unique {.var {marker1_col}}/{.var {marker2_col}} pair", 
                     x = "Found duplicated pairs: {.val {dup_pairs}}"))
  }
  if (symmetrise) {
    plot_data <- plot_data %>% bind_rows(plot_data %>% rename(`:=`(!!sym(marker1_col), 
                                                                   !!sym(marker2_col)), `:=`(!!sym(marker2_col), !!sym(marker1_col)))) %>% 
      distinct()
  }
  if (is.null(legend_range)) {
    legend_range <- max(abs(plot_data %>% pull(all_of(value_col)))) * 
      c(-1, 1)
  }
  
  plot_data <- plot_data %>% 
    mutate(across(all_of(value_col), 
                  ~case_when(.x < legend_range[1] ~ legend_range[1], 
                             .x > legend_range[2] ~ legend_range[2], 
                             TRUE ~ .x)))
  
  if (!is.null(size_col_transform)) {
    plot_data <- plot_data %>% mutate(`:=`(!!sym(size_col), 
                                           size_col_transform(!!sym(size_col))))
    size_col_label <- paste0(size_col, "_transformed")
  }
  else {
    size_col_label <- size_col
  }
  if (type == "tiles" || (cluster_rows || cluster_cols)) {
    plot_data_wide <- plot_data %>% {
      if (!is.null(size_col) && (size_col %in% names(.))) {
        select(., -all_of(size_col))
      }
      else {
        .
      }
    } %>% 
      pivot_wider(names_from = all_of(marker2_col), values_from = all_of(value_col), 
                  values_fill = 0) %>% column_to_rownames(marker1_col) %>% 
      as.matrix()
    if (cluster_rows && (type == "dots")) {
      rows_clust <- dist(plot_data_wide, method = clustering_distance_rows) %>% 
        hclust(method = clustering_method)
      plot_data <- plot_data %>% mutate(`:=`(!!sym(marker1_col), 
                                             factor(!!sym(marker1_col), levels = with(rows_clust, 
                                                                                      labels[order]))))
    }
    if (cluster_cols && (type == "dots")) {
      if (symmetrise && cluster_rows) {
        cols_clust <- rows_clust
      }
      else {
        cols_clust <- dist(t(plot_data_wide), method = clustering_distance_cols) %>% 
          hclust(method = clustering_method)
      }
      plot_data <- plot_data %>% mutate(`:=`(!!sym(marker2_col), 
                                             factor(!!sym(marker2_col), levels = with(cols_clust, 
                                                                                      labels[order]))))
    }
  }
  
  #' New: Add functionality to order factors by provided list of markers - Quick version (not yet added all possible tests)
  if (type == "dots" && (!cluster_rows || !cluster_cols) && !is.null(marker_order)) {
    plot_data_wide <- plot_data %>% {
      if (!is.null(size_col) && (size_col %in% names(.))) {
        select(., -all_of(size_col))
      }
      else {
        .
      }
    } %>% 
      pivot_wider(names_from = all_of(marker2_col), values_from = all_of(value_col), 
                  values_fill = 0) %>% column_to_rownames(marker1_col) %>% 
      as.matrix()
    
    plot_data <- plot_data %>% 
      mutate(`:=`(!!sym(marker1_col), 
                  factor(!!sym(marker1_col), levels = marker_order))
      ) |> 
      mutate(`:=`(!!sym(marker2_col), 
                  factor(!!sym(marker2_col), levels = marker_order)))
  }
  
  if (return_plot_data) {
    if (type == "dots") {
      return(plot_data)
    }
    else {
      return(plot_data_wide)
    }
  }
  
  
  if (type == "dots") {
    p <- ggplot(plot_data, aes(!!sym(marker1_col), !!sym(marker2_col), 
                               fill = !!sym(value_col), size = !!sym(size_col))) + 
      geom_point(shape = 21) + scale_size(range = size_range) + 
      scale_y_discrete(limits = rev) + scale_x_discrete(position = "top") + 
      scale_fill_gradientn(colors = colors, limits = legend_range) + 
      theme_bw() + 
      theme(
        axis.text.x = element_text(angle = 45, hjust = 0),
        axis.text = element_text(color = "black"),
        axis.title = element_blank(),
        axis.ticks = element_line(linewidth=0.25, color = "black"),
        axis.line = element_blank(),
        panel.grid = element_line(linewidth=0.15, color = "black"),
        panel.background = element_rect(linewidth=0.5, color = "black", fill = NULL),
        panel.border = element_blank(),
        legend.key = element_blank()
      ) + 
      coord_fixed() + 
      labs(size = size_col_label)
  }
  else {
    p <- plot_data_wide %>% ComplexHeatmap::pheatmap(breaks = seq(legend_range[1], 
                                                                  legend_range[2], length.out = 101), color = colorRampPalette(colors)(100), 
                                                     heatmap_legend_param = list(title = legend_title), 
                                                     cluster_rows = cluster_rows, cluster_cols = cluster_cols, 
                                                     clustering_distance_rows = clustering_distance_rows, 
                                                     clustering_distance_cols = clustering_distance_cols, 
                                                     clustering_method = clustering_method, ...)
  }
  return(p)
}
