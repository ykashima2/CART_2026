#' Validate inputs for DifferentialBubblePlot
#'
#' @param data Data frame to validate
#' @param effect_col Character. Column name containing effect size values
#' @param type Character. Type of analysis
#' @param alpha_cutoff Numeric. Significance cutoff
#' @param effect_cutoff Numeric. Effect size cutoff
#' @param color_palette Character vector. Color gradient
#' @param color_range Numeric vector. Color scale range
#' @param group_var Character. Grouping variable
#' @param facet_var Character. Faceting variable
#' @param size_range Numeric vector. Point size range
#' @param clustering_method Character. Clustering method
#' @param has_pvalues Logical. Whether data has p-values
#' @param marker_order Character vector. Marker order
#' @param show_dendrogram Logical. Show dendrogram
#' @param base_size Numeric. Base font size
#' @param pvalue_col Character. P-value column name
#' @param size_var Character. Custom size variable column name
#'
#' @return Nothing (invisible). Throws errors if validation fails
#' @keywords internal
#' @noRd
.validate_bubble_plot_input <- function(
  data,
  effect_col,
  type,
  alpha_cutoff,
  effect_cutoff,
  color_palette,
  color_range,
  group_var,
  facet_var,
  size_range,
  clustering_method,
  has_pvalues,
  marker_order,
  show_dendrogram,
  base_size,
  pvalue_col,
  size_var
) {
  if (!is.null(type) && type == "proximity") {
    required_cols <- c(effect_col)
  } else {
    required_cols <- c("marker", effect_col, group_var)
    if (has_pvalues) {
      p_col <- if (!is.null(pvalue_col)) pvalue_col else "p_adj"
      required_cols <- c(required_cols, p_col)
    }
  }

  required_cols <- required_cols[!is.null(required_cols)]

  if (length(required_cols) > 0) {
    sapply(
      required_cols,
      pixelatorR:::assert_col_in_data,
      data = data
    )
  }

  pixelatorR:::assert_class(data, "data.frame")
  pixelatorR:::assert_col_class(
    effect_col,
    data,
    "numeric",
    allow_null = TRUE
  )
  pixelatorR:::assert_single_value(effect_col, "string", allow_null = TRUE)
  pixelatorR:::assert_single_value(type, "string", allow_null = TRUE)
  pixelatorR:::assert_single_value(alpha_cutoff, "numeric")
  pixelatorR:::assert_vector(color_palette, "character")
  pixelatorR:::assert_within_limits(length(color_palette), c(2, Inf))
  pixelatorR:::assert_vector(color_range, "numeric", n = 2, allow_null = TRUE)
  # Only validate group_var for non-proximity analysis
  if (is.null(type) || type != "proximity") {
    pixelatorR:::assert_col_in_data(group_var, data, allow_null = TRUE)
    pixelatorR:::assert_col_class(group_var, data, "character")
  }
  pixelatorR:::assert_col_in_data(facet_var, data, allow_null = TRUE)
  pixelatorR:::assert_single_value(facet_var, "string", allow_null = TRUE)

  if (!is.null(facet_var)) {
    pixelatorR:::assert_col_class(facet_var, data, "character")
  }

  # Validate dendrogram parameters
  pixelatorR:::assert_single_value(show_dendrogram, "bool")
  pixelatorR:::assert_single_value(base_size, "numeric")

  # Validate pvalue_col parameter
  pixelatorR:::assert_single_value(pvalue_col, "string", allow_null = TRUE)
  if (!is.null(pvalue_col) && has_pvalues) {
    pixelatorR:::assert_col_in_data(pvalue_col, data)
    pixelatorR:::assert_col_class(pvalue_col, data, "numeric")
  }

  # Validate size_var parameter
  pixelatorR:::assert_single_value(size_var, "string", allow_null = TRUE)
  if (!is.null(size_var)) {
    pixelatorR:::assert_col_in_data(size_var, data)
    pixelatorR:::assert_col_class(size_var, data, "numeric")
  }

  if (!is.null(type)) {
    type <- match.arg(
      type,
      choices = c("abundance", "polarity", "proximity")
    )
  }

  invisible()
}

#' Select appropriate effect column based on analysis type
#'
#' @param type Character or NULL. Type of analysis
#' @param effect_col Character or NULL. User-specified effect column
#'
#' @return Character. The selected effect column name
#' @keywords internal
#' @noRd
.select_effect_col <- function(data, type, effect_col) {
  if (type == "proximity" && !is.null(type) && is.null(effect_col)) {
    proximity_effect_cols <- c("difference", "estimate", "diff_median")
    effect_col_in_data <- proximity_effect_cols %in% colnames(data)
    effect_col <- proximity_effect_cols[effect_col_in_data]
  } else if (type != "proximity" && !is.null(type)) {
    effect_col <- switch(
      type,
      "abundance" = "difference",
      "polarity" = "estimate",
      effect_col
    )
  } else if (is.null(effect_col)) {
    cli::cli_abort(
      "{.arg type} is not provided, {.arg effect_col} must be specified."
    )
  }

  return(effect_col)
}

.check_bubble_plot_packages <- function(show_dendrogram) {
  if (!show_dendrogram) {
    return(TRUE)
  }

  # Check if required packages are available
  expect_ggdendro()

  return(TRUE)
}

#' Create proximity bubble plot
#'
#' @param data Data frame with marker_1, marker_2 columns
#' @param effect_col Character. Effect column name
#' @param alpha_cutoff Numeric. Significance cutoff
#' @param effect_cutoff Numeric. Effect cutoff
#' @param color_palette Character vector. Color palette
#' @param color_range Numeric vector. Color range
#' @param size_range Numeric vector. Size range
#' @param clustering_method Character. Clustering method
#' @param has_pvalues Logical. Whether data has p-values
#' @param marker_order Character vector. Marker order
#' @param show_dendrogram Logical. Show dendrogram
#' @param base_size Numeric. Base font size
#' @param pvalue_col Character. P-value column name
#'
#' @return ggplot2 object or combined plot with dendrogram
#' @keywords internal
#' @noRd
.create_bubble_plot_proximity <- function(
  data,
  effect_col,
  alpha_cutoff,
  effect_cutoff,
  color_palette,
  color_range,
  size_range,
  clustering_method,
  has_pvalues,
  marker_order,
  show_dendrogram,
  base_size,
  pvalue_col,
  size_var
) {
  # Validate required columns for proximity analysis
  required_cols <- c("marker_1", "marker_2")
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    cli::cli_abort(
      c(
        "x" = "For {.arg type} = 'proximity', data must contain columns: {.val {required_cols}}",
        "i" = "Missing columns: {.val {missing_cols}}"
      )
    )
  }

  # Set default effect_col for proximity if not specified
  if (is.null(effect_col)) {
    effect_col <- "estimate"
  }

  # Use provided pvalue_col or auto-detect p-value column if has_pvalues is TRUE
  p_value_col <- NULL
  if (has_pvalues) {
    if (!is.null(pvalue_col)) {
      # Use user-specified column
      p_value_col <- pvalue_col
    } else {
      # Auto-detect p-value column
      if ("p_adj" %in% names(data)) {
        p_value_col <- "p_adj"
      } else if ("join_count_p" %in% names(data)) {
        p_value_col <- "join_count_p"
      } else {
        cli::cli_warn(
          "No p-value column found. Setting has_pvalues = FALSE."
        )
        has_pvalues <- FALSE
      }
    }
  }

  # Filter data by significance and effect size
  if (has_pvalues) {
    significant_data <- data |>
      filter(
        !!sym(p_value_col) <= alpha_cutoff,
        abs(!!sym(effect_col)) >= effect_cutoff
      )
  } else {
    significant_data <- data |>
      filter(
        abs(!!sym(effect_col)) >= effect_cutoff
      )
  }

  # Create symmetric data
  if (has_pvalues) {
    symmetric_data <- significant_data |>
      bind_rows(
        tibble(
          marker_1 = significant_data$marker_2,
          marker_2 = significant_data$marker_1,
          !!effect_col := significant_data[[effect_col]],
          !!p_value_col := significant_data[[p_value_col]]
        )
      ) |>
      mutate(
        !!p_value_col := ifelse(
          !!sym(p_value_col) == 0,
          .Machine$double.xmin,
          !!sym(p_value_col)
        )
      )
  } else {
    symmetric_data <- significant_data |>
      bind_rows(
        tibble(
          marker_1 = significant_data$marker_2,
          marker_2 = significant_data$marker_1,
          !!effect_col := significant_data[[effect_col]]
        )
      )
  }

  # Cap effect values if color_range is provided
  if (!is.null(color_range)) {
    symmetric_data <- symmetric_data |>
      mutate(
        !!effect_col := pmin(
          pmax(!!sym(effect_col), color_range[1]),
          color_range[2]
        )
      )
  }

  # Use provided marker order or calculate new one using hierarchical clustering
  marker_hc <- NULL # Store hierarchical clustering for dendrogram
  if (is.null(marker_order)) {
    unique_markers <- sort(
      unique(c(symmetric_data$marker_1, symmetric_data$marker_2))
    )
    marker_matrix <- matrix(
      NA,
      nrow = length(unique_markers),
      ncol = length(unique_markers),
      dimnames = list(unique_markers, unique_markers)
    )

    # Fill matrix with effect values
    for (i in seq_len(nrow(symmetric_data))) {
      marker_matrix[
        symmetric_data$marker_1[i],
        symmetric_data$marker_2[i]
      ] <- symmetric_data[[effect_col]][i]
    }
    marker_matrix[is.na(marker_matrix)] <- 0

    # Perform hierarchical clustering
    marker_hc <- hclust(
      dist(marker_matrix),
      method = clustering_method
    )
    marker_order <- marker_hc$labels[marker_hc$order]
  } else {
    # Validate that provided markers exist in the data
    all_markers_in_data <- sort(unique(c(
      symmetric_data$marker_1,
      symmetric_data$marker_2
    )))
    missing_markers <- setdiff(marker_order, all_markers_in_data)
    if (length(missing_markers) > 0) {
      cli::cli_warn(
        "The following markers in marker_order are not found in the data: {.val {missing_markers}}"
      )
    }
    # If custom order is provided, don't show dendrograms (no clustering performed)
    if (show_dendrogram) {
      cli::cli_warn(
        "Dendrogram cannot be shown when custom marker_order is provided. Set show_dendrogram = FALSE."
      )
      show_dendrogram <- FALSE
    }
  }

  # Create all combinations for complete matrix
  all_combinations <- expand.grid(
    marker_1 = marker_order,
    marker_2 = marker_order,
    stringsAsFactors = FALSE
  )

  all_combinations$key <- paste(
    all_combinations$marker_1,
    all_combinations$marker_2,
    sep = "_"
  )
  symmetric_data$key <- paste(
    symmetric_data$marker_1,
    symmetric_data$marker_2,
    sep = "_"
  )

  plot_data <- left_join(all_combinations, symmetric_data, by = "key")

  plot_data <- plot_data |>
    mutate(
      marker_1 = marker_1.x,
      marker_2 = marker_2.x
    ) |>
    select(-marker_1.x, -marker_2.x, -marker_1.y, -marker_2.y, -key)

  plot_data <- plot_data |>
    mutate(
      !!effect_col := ifelse(
        is.na(!!sym(effect_col)),
        0,
        !!sym(effect_col)
      )
    )

  # Create proximity plot
  p <- ggplot(
    plot_data,
    aes(
      x = as.numeric(factor(marker_2, levels = marker_order)),
      y = as.numeric(factor(marker_1, levels = marker_order)),
      fill = !!sym(effect_col)
    )
  )

  # Add stroke column to plot_data for conditional styling
  plot_data <- plot_data |>
    mutate(
      stroke_width = ifelse(marker_1 == marker_2, 0.8, 0.25)
    )

  if (has_pvalues) {
    if (!is.null(size_var)) {
      # Use custom size variable without transformation
      # Check if values are very small (possible p-values) and warn
      size_values <- plot_data[[size_var]]
      if (any(size_values > 0 & size_values < 0.1, na.rm = TRUE)) {
        cli::cli_warn(
          c(
            "!" = "Values in {.arg size_var} '{size_var}' appear to be very small (< 0.1).",
            "i" = paste(
              "If these are p-values, consider applying -log10 transformation",
              "before plotting for better visualization."
            )
          )
        )
      }

      p <- p +
        geom_point(
          data = plot_data,
          aes(
            size = !!sym(size_var),
            stroke = stroke_width
          ),
          shape = 21,
          alpha = 0.8,
          na.rm = FALSE
        ) +
        scale_size_continuous(
          name = size_var,
          range = size_range
        ) +
        guides(stroke = "none")
    } else {
      # Use default behavior: -log10 transformation of p-values
      p <- p +
        geom_point(
          data = plot_data,
          aes(
            size = -log10(!!sym(p_value_col)),
            stroke = stroke_width
          ),
          shape = 21,
          alpha = 0.8,
          na.rm = FALSE
        ) +
        scale_size_continuous(
          name = "-log10(p.adj)",
          range = size_range,
          limits = c(
            1.2,
            max(-log10(symmetric_data[[p_value_col]]))
          )
        ) +
        guides(stroke = "none")
    }
  } else {
    if (!is.null(size_var)) {
      # Use custom size variable when p-values are not available
      # Check if values are very small and warn
      size_values <- plot_data[[size_var]]
      if (any(size_values > 0 & size_values < 0.1, na.rm = TRUE)) {
        cli::cli_warn(
          c(
            "!" = "Values in {.arg size_var} '{size_var}' appear to be very small (< 0.1).",
            "i" = paste(
              "If these are p-values, consider applying -log10 transformation",
              "before plotting for better visualization."
            )
          )
        )
      }

      p <- p +
        geom_point(
          data = plot_data,
          aes(
            size = !!sym(size_var),
            stroke = stroke_width
          ),
          shape = 21,
          alpha = 0.8,
          na.rm = FALSE
        ) +
        scale_size_continuous(
          name = size_var,
          range = size_range
        ) +
        guides(stroke = "none")
    } else {
      # Use fixed size when p-values are not available
      p <- p +
        geom_point(
          data = plot_data,
          aes(stroke = stroke_width),
          size = mean(size_range),
          shape = 21,
          alpha = 0.8,
          na.rm = FALSE
        ) +
        guides(stroke = "none")
    }
  }

  # Set color limits
  if (!is.null(color_range)) {
    color_limits <- color_range
  } else {
    value_max <- max(
      abs(symmetric_data[[effect_col]]),
      na.rm = TRUE
    )
    color_limits <- c(-value_max, value_max)
  }

  color_ramp_palette <- colorRampPalette(color_palette)

  p <- p +
    scale_fill_gradientn(
      name = effect_col,
      colours = color_ramp_palette(100),
      limits = color_limits,
      na.value = "grey80",
      oob = scales::squish
    ) +
    scale_x_continuous(
      breaks = seq_along(marker_order),
      labels = marker_order,
      limits = c(0.5, length(marker_order) + 0.5),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      breaks = seq_along(marker_order),
      labels = marker_order,
      limits = c(0.5, length(marker_order) + 0.5),
      expand = c(0, 0)
    ) +
    theme_bw(base_size = base_size) +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1
      ),
      axis.text.y = element_text(angle = 0, hjust = 1),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    ) +
    labs(
      x = NULL,
      y = NULL
    )

  # Add dendrograms if requested for proximity analysis
  if (show_dendrogram && !is.null(marker_hc)) {
    # Check if required packages are available
    if (!.check_bubble_plot_packages(show_dendrogram)) {
      return(p)
    }

    # Create top dendrogram plot
    dend_data <- ggdendro::dendro_data(
      marker_hc,
      type = "rectangle"
    )
    top_dend <- ggplot(ggdendro::segment(dend_data)) +
      geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
      scale_x_continuous(
        breaks = seq_along(marker_order),
        labels = marker_order,
        limits = c(0.5, length(marker_order) + 0.5),
        expand = c(0, 0)
      ) +
      theme_void() +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(5, 5, 0, 5) # top, right, bottom, left
      )

    combined_plot <- top_dend /
      p +
      patchwork::plot_layout(heights = c(1, 4))
    return(combined_plot)
  }

  return(p)
}

#' Create abundance/polarity bubble plot
#'
#' @param data Data frame with marker, target columns
#' @param effect_col Character. Effect column name
#' @param alpha_cutoff Numeric. Significance cutoff
#' @param effect_cutoff Numeric. Effect cutoff
#' @param color_palette Character vector. Color palette
#' @param color_range Numeric vector. Color range
#' @param group_var Character. Grouping variable
#' @param facet_var Character. Faceting variable
#' @param size_range Numeric vector. Size range
#' @param clustering_method Character. Clustering method
#' @param has_pvalues Logical. Whether data has p-values
#' @param marker_order Character vector. Marker order
#' @param show_dendrogram Logical. Show dendrogram
#' @param base_size Numeric. Base font size
#' @param pvalue_col Character. P-value column name
#'
#' @return ggplot2 object or combined plot with dendrogram
#' @keywords internal
#' @noRd
.create_bubble_plot_abundance_polarity <- function(
  data,
  effect_col,
  alpha_cutoff,
  effect_cutoff,
  color_palette,
  color_range,
  group_var,
  facet_var,
  size_range,
  clustering_method,
  has_pvalues,
  marker_order,
  show_dendrogram,
  base_size,
  pvalue_col,
  size_var
) {
  # Check for duplicate marker-target combinations when group_var is not provided or equals "target"
  if (group_var == "target") {
    duplicate_check <- data |>
      group_by(marker, target) |>
      summarise(count = n(), .groups = "drop") |>
      filter(count > 1)

    if (nrow(duplicate_check) > 0) {
      cli::cli_abort(
        c(
          "x" = "Duplicate marker-target combinations detected.",
          "i" = "Use {.arg group_var} to specify another column for grouping."
        )
      )
    }
  }

  # Determine which p-value column to use
  p_value_col <- if (!is.null(pvalue_col)) pvalue_col else "p_adj"

  # Get significant markers
  if (has_pvalues) {
    significant_markers <- data |>
      filter(
        !!sym(p_value_col) <= alpha_cutoff,
        abs(!!sym(effect_col)) >= effect_cutoff
      ) |>
      pull(marker) |>
      unique()
  } else {
    significant_markers <- data |>
      filter(
        abs(!!sym(effect_col)) >= effect_cutoff
      ) |>
      pull(marker) |>
      unique()
  }

  if (length(significant_markers) == 0) {
    cli::cli_abort(
      "No markers found that satisfy the {.arg alpha_cutoff} or {.arg effect_cutoff} thresholds."
    )
  }

  select_cols <- c(group_var, "marker", effect_col)
  if (has_pvalues) {
    select_cols <- c(select_cols, p_value_col)
  }
  # Add facet_var if specified
  if (!is.null(facet_var)) {
    select_cols <- c(select_cols, facet_var)
  }

  # Add target column if it exists and is different from group_var (for potential auto-faceting)
  if ("target" %in% colnames(data) && group_var != "target") {
    select_cols <- c(select_cols, "target")
  }

  # Add size_var column if specified
  if (!is.null(size_var)) {
    select_cols <- c(select_cols, size_var)
  }

  data_matrix <- data |>
    filter(marker %in% significant_markers) |>
    select(all_of(select_cols))

  # Check if we should automatically facet by target (only if target column exists and is different from group_var)
  if (
    is.null(facet_var) &&
      group_var != "target" &&
      "target" %in% colnames(data_matrix)
  ) {
    unique_targets <- unique(data_matrix$target)
    if (length(unique_targets) > 1) {
      facet_var <- "target"
      cli::cli_inform(
        c(
          "i" = "Automatically faceting due to multiple unique {.arg target} values."
        )
      )
    }
  }

  # Create a matrix for clustering
  clustering_matrix <- data_matrix |>
    group_by(!!sym(group_var), marker) |>
    summarise(
      avg_value = mean(!!sym(effect_col), na.rm = TRUE),
      .groups = "drop"
    ) |>
    pivot_wider(
      names_from = marker,
      values_from = avg_value,
      values_fill = 0
    )

  # Get the clustering order
  marker_cols <- clustering_matrix |>
    select(-!!sym(group_var)) |>
    as.matrix()

  target_cols <- t(marker_cols)

  # Use provided marker order or cluster the plot data
  marker_clusters <- NULL
  target_clusters <- NULL

  if (is.null(marker_order)) {
    # Use hierarchical clustering for marker order
    if (nrow(marker_cols) >= 2 && ncol(marker_cols) >= 2) {
      marker_clusters <- hclust(
        dist(marker_cols),
        method = clustering_method
      )
      target_clusters <- hclust(
        dist(target_cols),
        method = clustering_method
      )

      marker_order <- rownames(target_cols)[target_clusters$order]
      target_order <- clustering_matrix[[group_var]][marker_clusters$order]
    } else {
      marker_order <- sort(colnames(marker_cols))
      target_order <- sort(clustering_matrix[[group_var]])
    }
  } else {
    # Use provided marker order and cluster targets only
    # Validate that provided markers exist in the data
    available_markers <- colnames(marker_cols)
    missing_markers <- setdiff(marker_order, available_markers)
    if (length(missing_markers) > 0) {
      cli::cli_warn(
        "The following markers in marker_order are not found in the data: {.val {missing_markers}}"
      )
      marker_order <- marker_order[marker_order %in% available_markers]
    }

    # Still cluster targets for optimal grouping
    if (nrow(marker_cols) >= 2) {
      marker_clusters <- hclust(
        dist(marker_cols),
        method = clustering_method
      )
      target_order <- clustering_matrix[[group_var]][marker_clusters$order]
    } else {
      target_order <- sort(clustering_matrix[[group_var]])
    }
  }

  plot_data <- data_matrix |>
    mutate(
      !!sym(group_var) := factor(!!sym(group_var), levels = target_order),
      marker = factor(marker, levels = marker_order)
    )

  if (has_pvalues) {
    plot_data_filtered <- plot_data |>
      filter(
        !!sym(p_value_col) <= alpha_cutoff,
        abs(!!sym(effect_col)) >= effect_cutoff
      )

    # Determine size variable for plotting
    if (!is.null(size_var)) {
      # Use custom size variable without transformation
      # Check if values are very small (possible p-values) and warn
      size_values <- plot_data_filtered[[size_var]]
      if (any(size_values > 0 & size_values < 0.1, na.rm = TRUE)) {
        cli::cli_warn(
          c(
            "!" = "Values in {.arg size_var} '{size_var}' appear to be very small (< 0.1).",
            "i" = paste(
              "If these are p-values, consider applying -log10 transformation",
              "before plotting for better visualization."
            )
          )
        )
      }

      plot <- ggplot(
        plot_data_filtered,
        aes(
          x = marker,
          y = !!sym(group_var),
          size = !!sym(size_var)
        )
      )
    } else {
      # Use default behavior: -log10 transformation of p-values
      plot_data_filtered <- plot_data_filtered |>
        mutate(
          log10_p = ifelse(
            !!sym(p_value_col) == 0,
            324,
            -log10(!!sym(p_value_col))
          )
        )

      plot <- ggplot(
        plot_data_filtered,
        aes(
          x = marker,
          y = !!sym(group_var),
          size = log10_p
        )
      )
    }
  } else {
    plot_data_filtered <- plot_data |>
      filter(
        abs(!!sym(effect_col)) >= effect_cutoff
      )

    if (!is.null(size_var)) {
      # Use custom size variable when p-values are not available
      # Check if values are very small and warn
      size_values <- plot_data_filtered[[size_var]]
      if (any(size_values > 0 & size_values < 0.1, na.rm = TRUE)) {
        cli::cli_warn(
          c(
            "!" = "Values in {.arg size_var} '{size_var}' appear to be very small (< 0.1).",
            "i" = paste(
              "If these are p-values, consider applying -log10 transformation",
              "before plotting for better visualization."
            )
          )
        )
      }

      plot <- ggplot(
        plot_data_filtered,
        aes(
          x = marker,
          y = !!sym(group_var),
          size = !!sym(size_var)
        )
      )
    } else {
      plot <- ggplot(
        plot_data_filtered,
        aes(
          x = marker,
          y = !!sym(group_var)
        )
      )
    }
  }

  if (has_pvalues) {
    plot <- plot +
      geom_point(
        aes(fill = !!sym(effect_col)),
        shape = 21,
        alpha = 0.8
      )
  } else {
    if (!is.null(size_var)) {
      # Use size_var for sizing when has_pvalues is FALSE
      plot <- plot +
        geom_point(
          aes(fill = !!sym(effect_col)),
          shape = 21,
          alpha = 0.8
        )
    } else {
      # Use fixed size when no size variable is provided
      plot <- plot +
        geom_point(
          aes(fill = !!sym(effect_col)),
          shape = 21,
          alpha = 0.8,
          size = mean(size_range)
        )
    }
  }

  if (!is.null(color_range)) {
    if (length(color_range) != 2) {
      cli::cli_abort(
        "color_range must be NULL or a numeric vector of length 2"
      )
    }
    color_limits <- color_range
  } else {
    value_max <- max(abs(plot_data_filtered[[effect_col]]), na.rm = TRUE)
    color_limits <- c(-value_max, value_max)
  }

  color_ramp_palette <- colorRampPalette(color_palette)

  plot <- plot +
    scale_fill_gradientn(
      name = effect_col,
      colours = color_ramp_palette(100),
      limits = color_limits,
      oob = scales::squish
    )

  if (has_pvalues || !is.null(size_var)) {
    # Determine the legend name for size scale
    size_legend_name <- if (!is.null(size_var)) {
      size_var
    } else {
      "-log10(p.adj)"
    }

    plot <- plot +
      scale_size_continuous(
        name = size_legend_name,
        range = size_range
      )
  }

  plot <- plot +
    theme_bw(base_size = base_size) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      legend.position = "right",
      plot.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      panel.grid.major = element_line(
        color = "gray80",
        linetype = "dotted",
        linewidth = 0.5
      ),
      panel.grid.minor = element_line(
        color = "gray80",
        linetype = "dotted",
        linewidth = 0.5
      )
    ) +
    labs(
      x = NULL,
      y = NULL
    )

  # Add faceting if facet_var is specified
  if (!is.null(facet_var)) {
    plot <- plot + facet_wrap(vars(!!sym(facet_var)), scales = "fixed")
  }

  # Add top dendrogram if requested for abundance/polarity analysis.
  should_show_top_dendrogram_ab_pol <- show_dendrogram &&
    !is.null(target_clusters)

  if (should_show_top_dendrogram_ab_pol) {
    # Check if required packages are available
    if (!.check_bubble_plot_packages(show_dendrogram)) {
      return(plot)
    }

    # Adjust main plot x-scale for alignment with top dendrogram.
    # The 'marker_order_for_plot_axis' variable (which is the internal 'marker_order' variable,
    # potentially updated by target_clusters) should be used by the main plot's x-axis.
    # This ensures the main plot's x-axis matches the dendrogram's.
    plot <- plot + scale_x_discrete(expand = expansion(mult = 0, add = 0.5))

    # Top dendrogram for markers (horizontal)
    top_dend_data <- ggdendro::dendro_data(
      target_clusters,
      type = "rectangle"
    )
    # Use the order and labels directly from target_clusters for the dendrogram axis
    dend_marker_order <- target_clusters$labels[target_clusters$order]
    top_dend <- ggplot(ggdendro::segment(top_dend_data)) +
      geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
      scale_x_continuous(
        breaks = seq_along(dend_marker_order),
        labels = dend_marker_order,
        limits = c(0.5, length(dend_marker_order) + 0.5),
        expand = c(0, 0)
      ) +
      theme_void() +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(5, 5, 0, 5) # top, right, bottom, left
      )

    # Combine plots using patchwork
    combined_plot <- top_dend /
      plot +
      patchwork::plot_layout(heights = c(1, 4))

    return(combined_plot)
  }

  return(plot)
}

#' Create a Differential Bubble Plot
#'
#' @description
#' Creates a bubble plot to visualize differential abundance or polarity data.
#' The plot shows markers (x-axis) versus targets or groups (y-axis) with point size
#' representing significance and color representing the magnitude and direction of change.
#'
#' @param data Data frame. Output of `RunDAA` or `RunDPA`. For proximity analysis, data should contain
#'   `marker_1` and `marker_2` columns instead of a single `marker` column.
#' @param effect_col Character. Column name in `data` containing the effect size values (e.g., log fold change).
#'   Default is NULL.
#' @param type Character or NULL. Type of analysis to visualize. Can be NULL, "abundance", "polarity", or "proximity".
#'   If NULL, `effect_col` must be specified. If "abundance", `effect_col` is automatically set to "difference".
#'   If "polarity", `effect_col` is automatically set to "estimate". If "proximity", creates a symmetric
#'   colocalization heatmap and expects data with `marker_1` and `marker_2` columns. Default is NULL.
#' @param alpha_cutoff Numeric. Significance cutoff for adjusted p-values. Default is 0.05.
#' @param effect_cutoff Numeric. Minimum absolute effect size to consider. Default is 0.1.
#' @param color_palette Character vector. Three colors for the gradient: negative, neutral, positive.
#' Must be a valid argument to `col2rgb()`.
#'   Default is c("#77BDDD", "white", "#FB2340").
#' @param color_range Numeric vector of length 2 or NULL. Range for color scale. If NULL, a symmetric
#'   range based on the data will be used. Default is NULL.
#' @param group_var Character. Column name in `data` to use for grouping on y-axis.
#'   Default is "target". Ignored for proximity analysis.
#' @param facet_var Character or NULL. Optional column name in `data` to use for faceting. Default is NULL.
#'   If NULL and there are multiple targets when group_var is not "target", will automatically facet by target.
#' @param size_range Numeric vector of length 2. Range for point sizes. Default is c(2, 10).
#' @param size_var Character or NULL. Column name in `data` to use for bubble sizing.
#'   If NULL, uses the default behavior: -log10 transformation of p-values when available.
#'   If specified, values are used directly without transformation. A warning is issued
#'   if values appear to be very small (< 0.1), suggesting they may be p-values that
#'   should be -log10 transformed. Default is NULL.
#' @param clustering_method Character. The agglomeration method for hierarchical clustering.
#'   Options: "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", or "centroid".
#'   Default is "ward.D2".
#' @param has_pvalues Logical. Whether the data contains p-values for significance testing.
#'   If TRUE, point sizes will be based on -log10(p-values). If FALSE, uses fixed point sizes.
#'   Default is TRUE.
#' @param marker_order Character vector or NULL. Optional character vector specifying the marker order.
#'   If NULL, markers will be ordered using hierarchical clustering. For proximity analysis, this
#'   controls the order of both x and y axes. Default is NULL.
#' @param show_dendrogram Logical. Whether to show hierarchical clustering dendrogram above the plot.
#'   Default is FALSE.
#' @param base_size Numeric. Base font size for the plot theme. Default is 12.
#' @param pvalue_col Character or NULL. Column name in `data` containing p-values.
#'   If NULL, defaults to \"p_adj\" for standard adjusted p-values. For proximity analysis, also accepts
#'   \"join_count_p\" if \"p_adj\" is not available. Default is NULL.
#'
#' @details
#' The function performs hierarchical clustering on both markers and grouping variables to organize
#' the plot. It filters data based on the provided significance and effect cutoffs.
#' Point sizes represent -log10(adjusted p-value) by default, or custom values when size_var is specified.
#' Colors represent the effect size.
#'
#' For proximity analysis, the function creates a symmetric heatmap where marker pairs are plotted
#' as points with colors representing colocalization scores and sizes representing significance.
#' The data is automatically symmetrized and organized using hierarchical clustering.
#'
#' The plot is returned as a ggplot2 object, allowing further customization with
#' additional ggplot2 parameters.
#'
#' @return A ggplot2 object. For proximity analysis, returns a symmetric heatmap with hierarchical clustering.
#'
#' @examples
#' \dontrun{
#' library(pixelatorR)
#' pxl_file <- system.file("extdata/five_cells",
#'   "five_cells.pxl",
#'   package = "pixelatorR"
#' )
#'
#' se <- ReadMPX_Seurat(pxl_file)
#' se <- merge(se, rep(list(se), 9), add.cell.ids = LETTERS[1:10])
#' se$sample <- c("T", "R", "C", "C", "C") %>% rep(times = 10)
#' se <- Seurat::NormalizeData(se, normalization.method = "CLR", margin = 2) |>
#'   JoinLayers()
#'
#' daa_markers <- RunDAA(se,
#'   contrast_column = "sample",
#'   targets = "T", reference = "C"
#' )
#'
#'
#' DifferentialBubblePlot(
#'   daa_markers,
#'   type = "abundance",
#'   alpha_cutoff = 0.05,
#'   effect_cutoff = 0.5
#' )
#' }
#'
#' @export
#'
DifferentialBubblePlot <- function(
  data,
  effect_col = NULL,
  type = NULL,
  alpha_cutoff = 0.05,
  effect_cutoff = 0.1,
  color_palette = c("#77BDDD", "white", "#FB2340"),
  color_range = NULL,
  group_var = "target",
  facet_var = NULL,
  size_range = c(2, 10),
  clustering_method = "ward.D2",
  has_pvalues = TRUE,
  marker_order = NULL,
  show_dendrogram = FALSE,
  base_size = 12,
  pvalue_col = NULL,
  size_var = NULL
) {
  # Validate all inputs
  .validate_bubble_plot_input(
    data = data,
    effect_col = effect_col,
    type = type,
    alpha_cutoff = alpha_cutoff,
    effect_cutoff = effect_cutoff,
    color_palette = color_palette,
    color_range = color_range,
    group_var = group_var,
    facet_var = facet_var,
    size_range = size_range,
    clustering_method = clustering_method,
    has_pvalues = has_pvalues,
    marker_order = marker_order,
    show_dendrogram = show_dendrogram,
    base_size = base_size,
    pvalue_col = pvalue_col,
    size_var = size_var
  )

  # Select the appropriate effect column
  effect_col <- .select_effect_col(
    data = data,
    type = type,
    effect_col = effect_col
  )

  # Handle proximity type
  if (!is.null(type) && type == "proximity") {
    return(.create_bubble_plot_proximity(
      data = data,
      effect_col = effect_col,
      alpha_cutoff = alpha_cutoff,
      effect_cutoff = effect_cutoff,
      color_palette = color_palette,
      color_range = color_range,
      size_range = size_range,
      clustering_method = clustering_method,
      has_pvalues = has_pvalues,
      marker_order = marker_order,
      show_dendrogram = show_dendrogram,
      base_size = base_size,
      pvalue_col = pvalue_col,
      size_var = size_var
    ))
  }

  # Handle abundance/polarity types
  return(.create_bubble_plot_abundance_polarity(
    data = data,
    effect_col = effect_col,
    alpha_cutoff = alpha_cutoff,
    effect_cutoff = effect_cutoff,
    color_palette = color_palette,
    color_range = color_range,
    group_var = group_var,
    facet_var = facet_var,
    size_range = size_range,
    clustering_method = clustering_method,
    has_pvalues = has_pvalues,
    marker_order = marker_order,
    show_dendrogram = show_dendrogram,
    base_size = base_size,
    pvalue_col = pvalue_col,
    size_var = size_var
  ))
}
