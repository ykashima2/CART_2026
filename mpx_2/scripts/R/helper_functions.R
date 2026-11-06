library(ggplot2)
library(dplyr)
library(patchwork)
library(MASS)

retro_pal_2 <- colorRampPalette(c("#77BDDD", "white", "#FB2340"))


# Utility functions ===================================================================================================

# preprocess data ----------------------------------------------------------------------------------------------------
process_data <- function(data) {
  # Preprocess
  data <- data |>
    Seurat::NormalizeData(
      normalization.method = "CLR",
      margin = 2
    ) |>
    Seurat::FindVariableFeatures(nFeatures = 50) |>
    Seurat::ScaleData(
      do.scale = TRUE,
      do.center = TRUE
    ) |>
    Seurat::RunPCA() |>
    RunUMAP(dims = 1:20, reduction = "pca")

  # Harmony and clustering
  data <- RunHarmony(data, group.by.vars = "sampleID", plot_convergence = FALSE)

  data <- data |>
    RunUMAP(reduction = "harmony", dims = 1:20) |>
    FindNeighbors(reduction = "harmony", dims = 1:20) |>
    FindClusters(random.seed = 1, resolution = 1)

  # data <- data |>
  #   subset(
  #     is_potential_doublet == FALSE
  #   )
  return(data)
}


# filter markers based on abundance of isotype controls ----------------------------------------------------------------
#' Filter markers based on isotype control abundance
#'
#' @description
#' Filters markers by comparing their median abundance to isotype control abundance.
#' Removes markers with median abundance below a threshold relative to isotype controls.
#' Also removes isotype controls and ACTB markers.
#'
#' @param data A Seurat object containing marker expression data
#' @param group_vars Named list of metadata variables and values to filter cells by (optional)
#' @param ratio Numeric threshold for marker/isotype ratio filtering (default: 1.5)
#'
#' @return Character vector of marker names that pass filtering criteria, or NULL if error
#'
#' @details
#' The function:
#' 1. If group_vars provided, checks if specified metadata variables exist
#' 2. If group_vars provided, filters cells based on metadata criteria
#' 3. Calculates median expression for each marker
#' 4. Filters markers with median > ratio * isotype_median
#' 5. Removes isotype controls and ACTB
#'
#' @examples
#' \dontrun{
#' # With grouping
#' markers <- filter_markers(seurat_obj,
#'   group_vars = c(celltype = "T cells"),
#'   ratio = 2
#' )
#'
#' # Without grouping
#' markers <- filter_markers(seurat_obj, ratio = 2)
#' }
filter_markers <- function(data, group_vars = NULL, ratio = 1.5) {
  # If group_vars provided, check if all variables exist in metadata
  if (!is.null(group_vars)) {
    missing_vars <- setdiff(names(group_vars), colnames(data@meta.data))
    if (length(missing_vars) > 0) {
      cli::cli_warn(paste(
        "Missing variables in metadata:",
        paste(missing_vars, collapse = ", ")
      ))
      return(NULL)
    }

    # Create filtering expression
    filter_expr <- paste(
      paste0(names(group_vars), " == '", group_vars, "'"),
      collapse = " & "
    )
  }

  tryCatch(
    {
      # Initialize cells_to_keep
      if (!is.null(group_vars)) {
        # Use FetchData to get the cells meeting our criteria
        cells_to_keep <- rownames(
          FetchData(data, vars = names(group_vars)) %>%
            filter(eval(parse(text = filter_expr)))
        )

        if (length(cells_to_keep) == 0) {
          cli_warn("No cells match the filtering criteria")
          return(NULL)
        }

        # Subset the data to keep only filtered cells
        data_subset <- data[, cells_to_keep]
      } else {
        # Use all cells if no grouping provided
        data_subset <- data
      }

      # Calculate medians and filter markers
      LayerData(
        data_subset,
        layer = "counts"
      ) %>%
        apply(MARGIN = 1, median) %>%
        enframe("marker", "median") %>%
        mutate(isotype_median = median(median[grepl("mIgG", marker)])) %>%
        filter(median > ratio * isotype_median, !grepl("mIgG|ACTB", marker)) %>%
        pull(marker)
    },
    error = function(e) {
      error_msg <- if (!is.null(group_vars)) {
        paste(
          "Error in filtering with conditions:",
          paste(names(group_vars), group_vars, sep = "=", collapse = ", "),
          ":",
          e$message
        )
      } else {
        paste("Error in filtering:", e$message)
      }
      cli::cli_warn(error_msg)
      return(NULL)
    }
  )
}

# bubble plot (abundance and polarity) ---------------------------------
#' Generate Differential Bubble Plot
#'
#' @description
#' Creates a bubble plot visualizing differential analysis results with
#' hierarchically clustered markers and targets. The plot shows differences between
#' conditions with bubble size representing significance and color representing
#' effect size.
#'
#' @param daa_results A data frame containing differential analysis results
#'   with columns: marker, target, p_adj
#' @param value_col Name of column containing values to plot (default: "difference")
#' @param p_adj_threshold Numeric threshold for adjusted p-value filtering (default: 0.05)
#' @param effect_threshold Numeric threshold for minimum absolute effect size (default: 0.1)
#' @param color_palette Function that generates color palette (default: retro_pal_2)
#' @param base_size Base font size for theme (default: 14)
#' @param title Plot title (default: "DA: treatment vs control")
#' @param grouping_var Optional grouping variable for faceting (default: NULL)
#' @param size_range Optional numeric vector of length 2 specifying the min and max values for point sizes (default: c(3, 12))
#'
#' @return A ggplot object representing the bubble plot
#'
#' @details
#' The function:
#' 1. Filters significant markers based on p-value and difference thresholds
#' 2. Performs hierarchical clustering on both markers and targets
#' 3. Creates a bubble plot with clustered ordering
#' 4. Customizes the plot appearance with specified themes and scales
#'
#' @examples
#' \dontrun{
#' plot <- differential_bubble_plot(
#'   data = my_results,
#'   value_col = "difference",
#'   p_adj_threshold = 0.05,
#'   difference_threshold = 0.1
#' )
#' }
differential_bubble_plot <- function(
  data,
  value_col = "difference",
  p_adj_threshold = 0.05,
  effect_threshold = 0.1,
  color_palette = retro_pal_2,
  base_size = 14,
  title = "DA: treatment vs control",
  grouping_var = NULL,
  size_range = c(3, 12)
) {
  # [Previous validation code remains the same...]
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }

  if (is.null(value_col) || !value_col %in% names(data)) {
    stop("'value_col' must be a valid column name in the data")
  }

  if (!is.null(grouping_var) && !grouping_var %in% names(data)) {
    stop("'grouping_var' must be NULL or a valid column name in the data")
  }

  required_cols <- c("marker", "target", "p_adj")
  missing_cols <- required_cols[!required_cols %in% names(data)]
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # Get significant markers
  significant_markers <- data |>
    filter(
      p_adj < p_adj_threshold,
      !grepl("mIgG|ACTB", marker),
      abs(!!sym(value_col)) > effect_threshold
    ) |>
    pull(marker) |>
    unique()

  if (length(significant_markers) == 0) {
    stop("No significant markers found with current thresholds")
  }

  # Create marker matrix for clustering
  if (!is.null(grouping_var)) {
    data_matrix <- data |>
      filter(marker %in% significant_markers) |>
      select(target, !!sym(grouping_var), marker, !!sym(value_col)) |>
      pivot_wider(
        names_from = marker,
        values_from = !!sym(value_col),
        values_fill = 0
      )
  } else {
    data_matrix <- data |>
      filter(marker %in% significant_markers) |>
      select(target, marker, !!sym(value_col)) |>
      pivot_wider(
        names_from = marker,
        values_from = !!sym(value_col),
        values_fill = 0
      )
  }

  # Cluster markers
  marker_cols <- if (!is.null(grouping_var)) {
    select(data_matrix, -target, -!!sym(grouping_var))
  } else {
    select(data_matrix, -target)
  }

  col_order <- marker_cols |>
    t() |>
    dist() |>
    hclust(method = "ward.D2") |>
    with(labels[order])

  # Cluster targets/groups
  row_order <- if (!is.null(grouping_var)) {
    data_matrix |>
      column_to_rownames(var = grouping_var) |>
      select(-target) |>
      dist() |>
      hclust(method = "ward.D2") |>
      with(labels[order])
  } else {
    data_matrix$target
  }

  # Generate plot data
  plot_data <- data |>
    filter(
      marker %in% significant_markers
    ) |>
    mutate(
      marker = factor(marker, levels = col_order),
      # Handle p-values of 0 by setting them to the smallest possible p-value
      log10_p = ifelse(p_adj == 0, 324, -log10(p_adj)), # 324 is -log10 of smallest possible double
      # Set size and fill to NA for non-significant results or those below effect threshold
      log10_p = ifelse(
        p_adj > p_adj_threshold | abs(!!sym(value_col)) <= effect_threshold,
        NA,
        log10_p
      ),
      plot_value = ifelse(
        p_adj > p_adj_threshold | abs(!!sym(value_col)) <= effect_threshold,
        NA,
        !!sym(value_col)
      )
    )

  if (!is.null(grouping_var)) {
    plot_data <- plot_data |>
      mutate(
        !!sym(grouping_var) := factor(!!sym(grouping_var), levels = row_order)
      )
  } else {
    plot_data <- plot_data |>
      mutate(target = factor(target, levels = row_order))
  }

  # Create plot - Fixed aesthetic mapping
  plot_mapping <- if (!is.null(grouping_var)) {
    aes(
      x = marker,
      y = !!sym(grouping_var),
      fill = plot_value,
      size = log10_p
    )
  } else {
    aes(
      x = marker,
      y = target,
      fill = plot_value,
      size = log10_p
    )
  }

  ggplot(plot_data, plot_mapping) +
    geom_point(shape = 21, alpha = 0.8) +
    theme_bw(base_size = base_size) +
    scale_fill_gradientn(
      colours = color_palette(100),
      limits = c(
        -max(abs(as.matrix(marker_cols))),
        max(abs(as.matrix(marker_cols)))
      )
    ) +
    scale_size_continuous(
      name = "-log10(p.adj)",
      range = size_range
    ) +
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
      fill = value_col,
      size = "p.adj",
      title = title,
      x = NULL,
      y = NULL
    )
}


# colocalization bubble plot ---------------------------------
#' Create a Colocalization Heatmap with Hierarchical Clustering
#'
#' This function creates a symmetrical colocalization heatmap with hierarchical clustering
#' for marker relationships. The plot includes adjustable point sizes based on p-values
#' and a color gradient for correlation estimates.
#'
#' @param data A data frame containing columns:
#'   - marker_1: First marker name
#'   - marker_2: Second marker name
#'   - estimate: Correlation or association estimate
#'   - p_adj: Adjusted p-value
#' @param color_palette A vector of 3 colors for low, mid, and high values (default: retro_pal_2(3))
#' @param base_size Base font size for the plot (default: 20)
#' @param p_adj_threshold Significance threshold for adjusted p-values (default: 0.05)
#' @param effect_threshold Minimum absolute effect size threshold (default: 0.1)
#' @param marker_order Optional character vector specifying the marker order
#' @param color_limits Optional numeric vector of length 2 specifying the min and max values for color scale (default: NULL)
#' @param size_range Optional numeric vector of length 2 specifying the min and max values for point sizes (default: c(1, 6))
#'
#' @return A ggplot object containing the correlation heatmap
#'
#' @examples
#' data <- data.frame(
#'   marker_1 = c("A", "B"),
#'   marker_2 = c("C", "D"),
#'   estimate = c(0.5, -0.3),
#'   p_adj = c(0.01, 0.05)
#' )
#' colocalization_bubble_plot(data)
colocalization_bubble_plot <- function(
  data,
  color_palette = colorRampPalette(c("#77BDDD", "white", "#FB2340")),
  base_size = 20,
  p_column = NULL,
  p_adj_threshold = 0.05,
  effect_column = "estimate",
  effect_threshold = 0.1,
  marker_order = NULL,
  color_range = NULL,
  size_range = c(1, 6)
) {
  # Check if p_column is provided
  if (is.null(p_column)) {
    stop("p_column must be specified to calculate bubble sizes using -log10(p)")
  }

  # Filter data by significance and effect size
  significant_data <- data |>
    filter(
      !!sym(p_column) < p_adj_threshold,
      abs(!!sym(effect_column)) > effect_threshold
    )

  # Create symmetric data
  # Create new data frame for binding
  new_data <- data.frame(
    marker_1 = significant_data$marker_2,
    marker_2 = significant_data$marker_1
  )
  new_data[[effect_column]] <- significant_data[[effect_column]]
  new_data[[p_column]] <- significant_data[[p_column]]

  symmetric_data <- bind_rows(significant_data, new_data)

  # Cap values to color range if specified
  if (!is.null(color_range)) {
    if (length(color_range) != 2) {
      stop(
        "color_range must be a vector of length 2 specifying lower and upper bounds"
      )
    }
    symmetric_data[[effect_column]] <- pmin(
      pmax(
        symmetric_data[[effect_column]],
        color_range[1]
      ),
      color_range[2]
    )
  }

  # Use provided marker order or calculate new one
  if (is.null(marker_order)) {
    unique_markers <- sort(unique(c(
      symmetric_data$marker_1,
      symmetric_data$marker_2
    )))
    marker_matrix <- matrix(
      NA,
      nrow = length(unique_markers),
      ncol = length(unique_markers),
      dimnames = list(unique_markers, unique_markers)
    )

    # Fill matrix with estimates
    for (i in 1:nrow(symmetric_data)) {
      marker_matrix[symmetric_data$marker_1[i], symmetric_data$marker_2[i]] <-
        symmetric_data[[effect_column]][i]
      # Add symmetric value
      marker_matrix[symmetric_data$marker_2[i], symmetric_data$marker_1[i]] <-
        symmetric_data[[effect_column]][i]
    }

    # Set diagonal to 0 and replace NAs
    diag(marker_matrix) <- 0
    marker_matrix[is.na(marker_matrix)] <- 0

    # Perform hierarchical clustering
    d <- dist(marker_matrix, method = "euclidean")
    hc <- hclust(d, method = "ward.D2")
    marker_order <- hc$labels[hc$order]
  }

  # Convert markers to factors with the clustered order
  symmetric_data$marker_1 <- factor(
    symmetric_data$marker_1,
    levels = marker_order
  )
  symmetric_data$marker_2 <- factor(
    symmetric_data$marker_2,
    levels = marker_order
  )

  # Calculate size based on -log10(p-value)
  symmetric_data$size <- -log10(symmetric_data[[p_column]])

  # Create the plot
  p <- ggplot(
    symmetric_data,
    aes(
      x = marker_1,
      y = marker_2,
      size = size,
      fill = .data[[effect_column]]
    )
  ) +
    geom_point(shape = 21, stroke = 0.25) +
    scale_size(range = size_range, name = "-log10(p-value)") +
    scale_fill_gradientn(
      colors = color_palette(100),
      limits = color_range
    ) +
    theme_minimal(base_size = base_size) +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1
      ),
      panel.grid.major = element_line(color = "gray90", linetype = "dotted"),
      panel.grid.minor = element_blank()
    ) +
    labs(
      x = NULL,
      y = NULL,
      fill = "Estimate"
    )

  return(p)
}
# function to plot 3D layouts side by side ---------------------------------
multi_panel_plot <- function(
  pop,
  components_to_plot,
  markers,
  n = 1,
  node_size = 4,
  rotation = 0,
  ncol = 2
) {
  # Input validation
  if (rotation < 0 || rotation > 360) {
    stop("Rotation must be between 0 and 360 degrees")
  }

  if (length(markers) == 0) {
    stop("At least one marker must be provided")
  }

  if (ncol < 1) {
    stop("ncol must be at least 1")
  }

  # Get cell graph data
  cg <- CellGraphs(pop)[[components_to_plot[n]]]
  counts <- CellGraphData(cg, "counts")
  g <- CellGraphData(cg, "cellgraph")
  layout <- CellGraphData(cg, "layout")$wpmds_3d

  # Create data frames for plotting with rotated coordinates
  if (rotation != 0) {
    # Convert rotation to radians
    angle_rad <- rotation * pi / 180

    # Rotation matrix around z-axis
    rot_matrix <- matrix(
      c(
        cos(angle_rad),
        -sin(angle_rad),
        0,
        sin(angle_rad),
        cos(angle_rad),
        0,
        0,
        0,
        1
      ),
      nrow = 3,
      byrow = TRUE
    )

    # Apply rotation to all three dimensions
    rotated_coords <- as.matrix(layout) %*% rot_matrix

    # Update layout coordinates
    layout <- rotated_coords
  }

  # Create custom color scale
  custom_colors <- c("#DFE4E0", "#D1C7D0", "#981E5C")

  # Create a list to store all plots
  plot_list <- list()

  # Create a plot for each marker
  for (marker in markers) {
    if (marker != "") {
      # Skip empty marker names
      # Calculate Gi statistics for the marker
      gi_res <- local_G(g, counts, k = 4, return_p_vals = TRUE)
      vals <- gi_res$gi_mat[, marker]
      vals[vals > quantile(vals, 0.99)] <- quantile(vals, 0.99)
      max_abs_val <- max(abs(vals))

      # Create data frame for plotting
      plot_data <- data.frame(
        x = layout[, 1],
        y = layout[, 2],
        value = vals
      )

      # Create plot
      p <- ggplot(plot_data, aes(x = x, y = y, color = value)) +
        geom_point(size = node_size) +
        scale_color_gradientn(
          colors = custom_colors,
          limits = c(-max_abs_val, max_abs_val)
        ) +
        theme_minimal() +
        theme(
          axis.text = element_blank(),
          axis.title = element_blank(),
          panel.grid = element_blank(),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 14)
        ) +
        ggtitle(marker) +
        coord_fixed()

      plot_list[[marker]] <- p
    }
  }

  # Combine plots using patchwork with specified number of columns
  if (length(plot_list) > 0) {
    combined_plot <- wrap_plots(plot_list, ncol = ncol)
    return(combined_plot)
  } else {
    stop("No valid markers provided")
  }
}


# function to load a list of packages ----------------------------------------
load_packages <- function(packages) {
  for (package in packages) {
    library(package, character.only = TRUE)
  }
}

# read markdown tables in a df
read_markdown <- function(file, trim_ws = TRUE, ...) {
  if (length(file) > 1) {
    lines <- file
  } else {
    lines <- readr::read_lines(file)
  }
  lines <- lines[!grepl("^[\\:\\s\\+\\-\\=\\_\\|]*$", lines, perl = TRUE)]
  lines <- gsub("(^\\s*?\\|)|(\\|\\s*?$)", "", lines)
  readr::read_delim(
    paste(lines, collapse = "\n"),
    delim = "|",
    trim_ws = trim_ws,
    ...
  )
}


# density scatter plot ------------------------------------------------------
.get2Ddensity <- function(x, y, n = 500, ...) {
  dens_grid <- MASS::kde2d(x, y, n = n, ...)
  dens <- dens_grid$z[cbind(
    findInterval(x, dens_grid$x),
    findInterval(y, dens_grid$y)
  )]
  return(dens)
}

DensityScatterPlotDF <- function(
  data,
  x,
  y,
  facet_vars = NULL,
  plot_gate = NULL,
  scale_density = TRUE,
  margin_density = FALSE,
  coord_fixed = TRUE,
  pt_size = 1,
  alpha = 1,
  grid_n = 500,
  colors = NULL,
  y_limits = NULL,
  ...
) {
  # Validate input parameters
  stopifnot(
    "'data' must be a data.frame" = is.data.frame(data),
    "'x' must be a character of length 1" = is.character(x) && length(x) == 1,
    "'y' must be a character of length 1" = is.character(y) && length(y) == 1,
    "'x' must be available in the data" = x %in% colnames(data),
    "'y' must be available in the data" = y %in% colnames(data),
    "'facet_vars' must either be NULL or be a character vector with 1 or 2 elements" = is.null(
      facet_vars
    ) |
      length(facet_vars) %in% 1:2,
    "Variables in 'facet_vars' must be available in the data" = is.null(
      facet_vars
    ) |
      all(facet_vars %in% colnames(data))
  )

  # Prepare plot data
  plot_data <- data %>%
    rename(marker1 = !!x, marker2 = !!y) %>%
    group_by_at(facet_vars) %>%
    mutate(dens = .get2Ddensity(marker1, marker2, n = grid_n, ...)) %>%
    ungroup()

  if (isTRUE(scale_density)) {
    plot_data <- plot_data %>%
      group_by_at(facet_vars) %>%
      mutate(dens = dens / max(dens)) %>%
      ungroup()
  }

  # Set plot theme
  if (is.null(facet_vars)) {
    plot_theme <- theme_bw() +
      theme(
        panel.grid = element_blank(),
        panel.spacing = unit(0, "lines"),
        plot.margin = unit(c(0, 0, 0, 0), "lines")
      )
  } else {
    plot_theme <- theme_bw() +
      theme(panel.grid = element_blank())
  }

  # Determine plot range
  x_range <- range(plot_data$marker1)
  y_range <- if (is.null(y_limits)) range(plot_data$marker2) else y_limits

  if (!is.null(plot_gate)) {
    x_range <- range(c(x_range, plot_gate$xmin, plot_gate$xmax))
    if (is.null(y_limits)) {
      y_range <- range(c(y_range, plot_gate$ymin, plot_gate$ymax))
    }
  }

  # Make plot
  plot <- plot_data %>%
    ggplot(aes(marker1, marker2, color = dens)) +
    geom_hline(yintercept = 0, color = "gray") +
    geom_vline(xintercept = 0, color = "gray") +
    geom_point(
      size = pt_size,
      alpha = alpha,
      show.legend = FALSE
    ) +
    scale_x_continuous(limits = x_range) +
    scale_y_continuous(limits = y_range) +
    plot_theme +
    ggplot2::labs(x = x, y = y)

  if (isTRUE(coord_fixed) && isFALSE(margin_density)) {
    plot <- plot + coord_fixed()
  }

  if (!is.null(colors)) {
    plot <- plot + scale_color_gradientn(colors = colors)
  } else {
    plot <- plot + scale_color_viridis_c(option = "H")
  }

  # Facet
  if (!is.null(facet_vars)) {
    if (length(facet_vars) == 1) {
      plot <- plot + facet_grid(rows = vars(!!!syms(facet_vars)))
    }
    if (length(facet_vars) == 2) {
      plot <- plot +
        facet_grid(
          rows = vars(!!!syms(facet_vars[1])),
          cols = vars(!!!syms(facet_vars[2]))
        )
    }
  }

  # Add gates
  if (!is.null(plot_gate)) {
    if (!is.null(facet_vars)) {
      join_vars <- intersect(facet_vars, colnames(plot_gate))
      if (length(join_vars) > 0) {
        gate_label <- plot_gate %>%
          left_join(plot_data, by = join_vars)
      } else {
        gate_label <- plot_gate %>%
          cross_join(plot_data)
      }
    } else {
      gate_label <- plot_data %>%
        cross_join(plot_gate)
    }

    gate_label <- gate_label %>%
      mutate(
        in_gate = marker1 > xmin &
          marker1 < xmax &
          marker2 > ymin &
          marker2 < ymax
      ) %>%
      group_by_at(c(facet_vars, "in_gate", "xmin", "xmax", "ymin", "ymax")) %>%
      summarise(
        n = n(),
        .groups = "drop"
      ) %>%
      group_by_at(c(facet_vars, "xmin", "xmax", "ymin", "ymax")) %>%
      mutate(p = 100 * n / sum(n)) %>%
      ungroup() %>%
      filter(in_gate) %>%
      mutate(
        label = paste(round(p, 1), "%"),
        x = (xmax + xmin) / 2,
        y = ymax
      )

    plot <- plot +
      geom_rect(
        data = plot_gate,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        inherit.aes = FALSE,
        fill = "transparent",
        color = "black",
        linetype = "dashed"
      ) +
      geom_text(
        data = gate_label,
        aes(x = x, y = y, label = label),
        inherit.aes = FALSE,
        color = "#FF4E88",
        vjust = 1,
        size = 6
      )
  }

  # Add marginal density
  if (isTRUE(margin_density)) {
    density_plot_x <- plot_data %>%
      ggplot(aes(marker1)) +
      geom_density(
        fill = "lightgray",
        color = NA
      ) +
      scale_x_continuous(limits = x_range) +
      scale_y_continuous(expand = expansion(0)) +
      theme_void() +
      theme(
        panel.spacing = unit(0, "lines"),
        plot.margin = unit(c(0, 0, 0, 0), "lines")
      )

    density_plot_y <- plot_data %>%
      ggplot(aes(marker2)) +
      geom_density(
        fill = "lightgray",
        color = NA
      ) +
      scale_x_continuous(limits = y_range) +
      scale_y_continuous(expand = expansion(0)) +
      theme_void() +
      theme(
        panel.spacing = unit(0, "lines"),
        plot.margin = unit(c(0, 0, 0, 0), "lines")
      ) +
      coord_flip()

    plot <- density_plot_x +
      plot_spacer() +
      plot +
      density_plot_y +
      plot_layout(
        widths = c(5, 1),
        heights = c(1, 5)
      )
  }

  return(plot)
}


# polarization plot ---------------------------------------------------------
generate_polarization_plot <- function(
  type_of_cells,
  targets,
  reference,
  contrast_column = "culture",
  pol_scores_filtered,
  mpx_data
) {
  pol_test <- RunDPA(
    pol_scores_filtered %>%
      filter(manual_annotation %in% type_of_cells),
    target = targets,
    reference = reference,
    contrast_column = contrast_column
  )

  median_clr_values <- FetchData(
    mpx_data,
    vars = c("sample", "manual_annotation", "culture", rownames(mpx_data)),
    layer = "data"
  ) %>%
    pivot_longer(
      cols = -c(sample, manual_annotation, culture),
      names_to = "marker",
      values_to = "abundance"
    ) %>%
    filter(culture %in% c(reference, targets)) %>%
    group_by(sample, culture, marker) %>%
    summarise(median = median(abundance))

  pol_test <- pol_test %>%
    left_join(
      median_clr_values,
      by = c("marker" = "marker", "target" = "culture")
    )

  pol_plot <- ggplot(
    pol_test,
    aes(estimate, -log10(p_adj), color = median, size = median, label = marker)
  ) +
    geom_point() +
    scale_color_gradient2(
      low = "#0066FFFF",
      mid = "#74D055FF",
      high = "#FF0000FF",
      midpoint = 3,
      name = "Median Abundance (CLR)"
    ) +
    geom_vline(
      xintercept = c(-1.96, 1.96),
      color = "grey",
      linetype = "dashed"
    ) +
    geom_hline(yintercept = -log10(0.05), color = "grey", linetype = "dashed") + # Significance line
    geom_vline(xintercept = 0, color = "black", linetype = "solid") +
    geom_text_repel(max.overlaps = 5, size = 6) + # Increased font size for text labels
    theme_bw() +
    facet_wrap(~target, scales = "free") +
    # ggplot2::labs(title = paste("Polarization Analysis in", type_of_cells, "Cells (TIL: 10049)"), x = "Effect Size versus Not Co-cultured (Estimate)", y = "-log10(Adjusted P-value)") +
    theme(
      axis.text = element_text(size = 18),
      axis.title = element_text(size = 20),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 16),
      strip.text = element_text(size = 18, face = "bold"),
      strip.background = element_rect(fill = "white")
    ) + # Increased font sizes for axis text and titles
    theme(legend.position = "bottom") +
    scale_size_continuous(name = "Median Abundance (CLR)")

  return(list(pol_test = pol_test, pol_plot = pol_plot))
}

type_of_cells <- "CD8T"
targets <- c("10049til_10049cc", "10049til_10049cc_sgB2M", "10049til_10378cc")
reference <- "10049til_no_cc"
contrast_column <- "culture"

# pol_plot_cd8_49 <- generate_polarization_plot(type_of_cells, targets, reference, contrast_column, pol_scores_filtered, mpx_data)

# Colocalization -----------------------------------------------------------

generate_colocalization_plot <- function(
  type_of_cells,
  targets,
  reference,
  contrast_column = "culture",
  col_scores_filtered,
  mpx_data
) {
  col_test <- RunDCA(
    col_scores_filtered %>%
      filter(manual_annotation %in% type_of_cells),
    target = targets,
    reference = reference,
    contrast_column = contrast_column,
    cl = 8
  )

  median_clr_values <- FetchData(
    mpx_data,
    vars = c("sample", "manual_annotation", "culture", rownames(mpx_data)),
    layer = "data"
  ) %>%
    pivot_longer(
      cols = -c(sample, manual_annotation, culture),
      names_to = "marker",
      values_to = "abundance"
    ) %>%
    filter(culture %in% c(reference, targets)) %>%
    group_by(sample, culture, marker) %>%
    summarise(median = median(abundance))

  col_test <- col_test %>%
    left_join(
      median_clr_values,
      by = c("marker_1" = "marker", "target" = "culture")
    ) %>%
    mutate(marker_1_abund = median) %>%
    left_join(
      median_clr_values,
      by = c("marker_2" = "marker", "target" = "culture")
    ) %>%
    mutate(marker_2_abund = median.y) %>%
    mutate(abund_avg = (marker_1_abund + marker_2_abund) / 2) %>%
    dplyr::select(
      estimate,
      target,
      reference,
      p_adj,
      marker_1,
      marker_2,
      abund_avg
    )

  col_plot <- ggplot(
    col_test,
    aes(
      estimate,
      -log10(p_adj),
      color = abund_avg,
      size = abund_avg,
      label = paste(marker_1, marker_2, sep = "/")
    )
  ) +
    geom_point() +
    scale_color_gradient2(
      low = "#0066FFFF",
      mid = "#74D055FF",
      high = "#FF0000FF",
      midpoint = 3,
      name = "Median Abundance (CLR)"
    ) +
    geom_vline(
      xintercept = c(-1.96, 1.96),
      color = "grey",
      linetype = "dashed"
    ) +
    geom_hline(yintercept = -log10(0.05), color = "grey", linetype = "dashed") + # Significance line
    geom_vline(xintercept = 0, color = "black", linetype = "solid") +
    geom_text_repel(max.overlaps = 5, size = 6) + # Increased font size for text labels
    theme_bw() +
    facet_wrap(~target, scales = "free") +
    theme(
      axis.text = element_text(size = 18),
      axis.title = element_text(size = 20),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 16),
      strip.text = element_text(size = 18, face = "bold"),
      strip.background = element_rect(fill = "white")
    ) + # Increased font sizes for axis text and titles
    theme(legend.position = "bottom") +
    scale_size_continuous(name = "Median Abundance (CLR)")

  return(list(col_test = col_test, col_plot = col_plot))
}

type_of_cells <- "CD8T"
targets <- c("10049til_10049cc", "10049til_10049cc_sgB2M", "10049til_10378cc")
reference <- "10049til_no_cc"
contrast_column <- "culture"

# pol_plot_cd8_49 <- generate_polarization_plot(type_of_cells, targets, reference, contrast_column, pol_scores_filtered, mpx_data)

# Col scores markers --------------------------------------------------------

get_col_scores <- function(col_scores_filtered, culture_type, cell_type) {
  markers_1 <- col_scores_filtered %>%
    filter(culture %in% culture_type & manual_annotation %in% cell_type) %>%
    ungroup() %>%
    dplyr::select(marker_1) %>%
    unique() %>%
    pull(marker_1)
  markers_2 <- col_scores_filtered %>%
    filter(culture %in% culture_type & manual_annotation %in% cell_type) %>%
    ungroup() %>%
    dplyr::select(marker_2) %>%
    unique() %>%
    pull(marker_2)
  all_markers <- unique(c(markers_1, markers_2))

  return(all_markers)
}

# Col score converter
col_mat_gen <- function(data, all_markers) {
  data %>%
    dplyr::select(marker_1, marker_2, pearson_z) %>%
    group_by(marker_1, marker_2) %>%
    summarise(median = median(pearson_z)) %>%
    mutate(median = ifelse(median > 15, 15, median)) %>%
    mutate(median = ifelse(median < -15, -15, median)) %>%
    {
      bind_rows(
        .,
        rename(., marker_2 = marker_1, marker_1 = marker_2),
        tibble(marker_1 = all_markers, marker_2 = all_markers, median = NA)
      )
    } %>%
    pivot_wider(
      names_from = marker_2,
      values_from = median,
      values_fill = list(median = NA)
    ) %>%
    tibble::column_to_rownames("marker_1") %>%
    as.matrix()
}


# Differential Abundance ---------------------------------------------------

# Function to generate differential abundance volcanos
generate_da_volcano <- function(
  mpx_data,
  target_list,
  reference,
  target_col,
  cell_type,
  output_path
) {
  # Filter data for specified targets, reference, and cell type
  tmp <- mpx_data %>%
    subset(subset = culture %in% c(target_list, reference)) %>%
    subset(subset = manual_annotation %in% cell_type)

  # Initialize list to store differential abundance results
  da_proteins <- list()

  # Loop through each target to find markers
  for (target in target_list) {
    # Set identity for the current target
    tmp <- SetIdent(tmp, value = target_col)

    # Find markers for the current target
    da_proteins[[target]] <- tmp %>%
      subset(subset = culture == c(target, reference)) %>%
      FindMarkers(
        test.use = "wilcox",
        ident.1 = target,
        ident.2 = reference,
        mean.fxn = rowMeans,
        fc.name = "difference",
        min.pct = 0,
        logfc.threshold = 0
      ) %>%
      rownames_to_column("marker")
  }

  # Bind rows to combine results from all targets
  da_proteins <- bind_rows(da_proteins, .id = "target")

  # Fetch and prepare CLR data for median abundance
  clr_data <- FetchData(
    mpx_data,
    vars = c("sample", "manual_annotation", "culture", rownames(mpx_data)),
    layer = "data"
  ) %>%
    pivot_longer(
      cols = -c(sample, manual_annotation, culture),
      names_to = "marker",
      values_to = "abundance"
    ) %>%
    filter(culture %in% c(reference, targets)) %>%
    group_by(sample, culture, marker) %>%
    summarise(median = median(abundance))

  # Join CLR data with differential abundance results
  da_proteins <- da_proteins %>%
    left_join(clr_data, by = c("marker" = "marker", "target" = "culture")) %>%
    filter(!marker %in% c("mIgG1", "mIgG2a", "mIgG2b"))

  # Plotting the volcano plot
  da_volcano_plot <- ggplot(
    da_proteins,
    aes(
      x = difference,
      y = -log10(p_val_adj),
      color = median,
      size = median,
      label = marker
    )
  ) +
    geom_point() +
    scale_color_gradient2(
      low = "#0066FFFF",
      mid = "#74D055FF",
      high = "#FF0000FF",
      midpoint = 3,
      name = "Median Abundance (CLR)"
    ) +
    scale_size(range = c(1, 4)) + # Size scale for points
    theme_bw() +
    geom_hline(yintercept = -log10(0.05), color = "grey", linetype = "solid") + # Significance line
    geom_vline(xintercept = 0, color = "grey", linetype = "solid") +
    geom_text_repel(max.overlaps = 4, size = 5) +
    geom_vline(xintercept = c(-0.5, 0.5), color = "grey", linetype = "dashed") +
    labs(
      x = "Difference",
      y = "-log10(Adjusted P-value)",
      color = "Median Abundance",
      size = "Median Abundance"
    ) +
    facet_wrap(~target, scales = "free") +
    theme(
      axis.text = element_text(size = 18),
      axis.title = element_text(size = 20),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 16),
      strip.text = element_text(size = 18, face = "bold"),
      strip.background = element_rect(fill = "white")
    ) + # Increased font sizes for axis text and titles
    theme(legend.position = "bottom") +
    scale_size_continuous(name = "Median Abundance (CLR)")

  # Save the plot
  ggsave(output_path, da_volcano_plot, width = 16, height = 8)
}

# generate_da_volcano(mpx_data, target_list, reference, target_col, cell_type, "PCA_0012/results/Differential_abundance/da_volcano_cd8_t_49.png")

#' Create a Lollipop Plot for Differential Abundance Analysis
#'
#' @description
#' Creates a lollipop plot to visualize differential abundance analysis results. The plot shows markers
#' on the y-axis and their corresponding values on the x-axis, with points connected to the baseline
#' by lines. Different conditions are represented by different colors, and significance is shown
#' through point and line transparency.
#'
#' @param data A data frame containing the analysis results
#' @param condition_col Character. Name of the column containing condition labels
#' @param value_col Character. Name of the column containing the values to plot
#' @param top_n Integer. Number of top markers to display (default: 20)
#' @param consistency Character. One of "all", "consistent", or "inconsistent". Filters markers based on their consistency across conditions
#' @param threshold Numeric. Threshold for considering an effect size as meaningful (default: 0.1)
#' @param p_adj_threshold Numeric. Adjusted p-value threshold for significance (default: 0.05)
#' @param point_size Numeric. Size of the points in the plot (default: 2)
#' @param line_size Numeric. Thickness of the lines in the plot (default: 1)
#' @param show_zero_line Logical. Whether to show a dashed line at y=0 (default: TRUE)
#' @param title Character. Main title for the plot (default: "Differential abundance analysis")
#' @param y_lab Character. Label for the y-axis (default: NULL)
#' @param color_palette Named vector. Colors for different conditions (default: specific colors for ibrutinib, pirtobrutinib, and venetoclax)
#'
#' @details
#' The function requires the following columns in the input data frame:
#' - marker: Marker identifiers
#' - p_adj: Adjusted p-values
#' - The columns specified in condition_col and value_col
#'
#' Markers are considered consistent when they:
#' - Are significant across all conditions (p_adj < p_adj_threshold)
#' - Show effects in the same direction (all positive or all negative)
#' - Have effect sizes above the specified threshold
#'
#' @return A ggplot object representing the lollipop plot
#'
#' @examples
#' \dontrun{
#' lollipop_plot(
#'   data = results_df,
#'   condition_col = "treatment",
#'   value_col = "log2FC",
#'   top_n = 15,
#'   consistency = "consistent",
#'   threshold = 0.1
#' )
#' }
#'
lollipop_plot <- function(
  data,
  condition_col = "treatment",
  value_col = "difference",
  top_n = 20,
  consistency = c("all", "consistent", "inconsistent"),
  threshold = 0.1,
  p_adj_threshold = 0.05,
  point_size = 2,
  line_size = 1,
  show_zero_line = TRUE,
  title = "Differential abundance analysis",
  y_lab = NULL,
  color_palette = c(
    "ibrutinib" = "#E69F00",
    "pirtobrutinib" = "#56B4E9",
    "venetoclax" = "#009E73"
  )
) {
  # Input validation
  consistency <- match.arg(consistency)
  required_cols <- c("marker", "p_adj", condition_col, value_col)
  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # Calculate consistency for each marker with the threshold condition
  marker_consistency <- data %>%
    group_by(marker) %>%
    summarize(
      n_significant = sum(p_adj < p_adj_threshold),
      n_conditions = n(),
      has_different_signs = min(.data[[value_col]]) * max(.data[[value_col]]) <
        0,
      has_small_effect = any(abs(.data[[value_col]]) < threshold),
      consistent = if_else(
        n_significant == n_conditions &&
          !has_different_signs &&
          !has_small_effect,
        "consistent",
        "inconsistent"
      ),
      max_abs_difference = max(abs(.data[[value_col]]))
    ) %>%
    select(
      -has_different_signs,
      -has_small_effect,
      -n_significant,
      -n_conditions
    )

  if (consistency != "all") {
    marker_consistency <- marker_consistency %>%
      filter(consistent == consistency)
  }

  if (nrow(marker_consistency) == 0) {
    stop("No markers found matching the specified criteria")
  }

  top_markers <- marker_consistency %>%
    top_n(top_n, max_abs_difference) %>%
    arrange(desc(max_abs_difference)) %>%
    pull(marker)

  # Add significance information to plot data
  plot_data <- data %>%
    filter(marker %in% top_markers) %>%
    mutate(
      marker = factor(marker, levels = rev(top_markers)),
      condition = factor(.data[[condition_col]]),
      is_significant = p_adj < p_adj_threshold,
      # Create alpha values based on significance
      point_alpha = if_else(is_significant, 1, 0.2)
    ) %>%
    filter(!grepl("^mIgG|ACTB", marker))

  # Create subtitle based on consistency selection
  subtitle <- switch(
    consistency,
    "all" = "Top markers (abs. difference)",
    "consistent" = sprintf(
      "Top consistent markers (abs. difference ≥ %.2f)",
      threshold
    ),
    "inconsistent" = sprintf(
      "Top inconsistent markers (abs. difference < %.2f or mixed signs)",
      threshold
    )
  )

  # Build plot
  p <- ggplot(
    plot_data,
    aes(x = marker, y = .data[[value_col]], color = condition)
  )

  if (show_zero_line) {
    p <- p + geom_hline(yintercept = 0, linetype = "dashed", color = "gray50")
  }

  p <- p +
    # Add alpha aesthetic to points and lines
    geom_point(
      aes(alpha = point_alpha),
      position = position_dodge(width = 0.5),
      size = point_size
    ) +
    geom_linerange(
      aes(ymin = 0, ymax = .data[[value_col]], alpha = point_alpha),
      position = position_dodge(width = 0.5),
      size = line_size
    ) +
    # Remove alpha from legend
    guides(alpha = "none") +
    coord_flip() +
    scale_color_manual(
      values = color_palette,
      name = condition_col
    ) +
    # Define specific alpha values
    scale_alpha_identity() +
    labs(
      title = title,
      subtitle = subtitle,
      x = "",
      y = y_lab
    ) +
    theme_minimal(base_size = 16) +
    guides(color = guide_legend(reverse = FALSE)) +
    theme(
      legend.position = "bottom",
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
      axis.title = element_text(size = 16)
    )

  return(p)
}

polarity_abundance_plot <- function(pop, dpa_results, daa_results) {
  # Join DPA and DAA data
  dpa_join_data <- dpa_results |>
    select(target, marker, estimate, p_adj) |>
    filter(target == pop) |>
    rename(estimate_dpa = estimate, p_adj_dpa = p_adj)

  daa_join_data <- daa_results |>
    select(target, marker, difference, p_adj) |>
    filter(target == pop) |>
    rename(estimate_daa = difference, p_adj_daa = p_adj)

  plot_data <- dpa_join_data |>
    left_join(daa_join_data, by = "marker")

  abund_max <- max(abs(plot_data$estimate_daa)) * 1.1
  pol_max <- max(abs(plot_data$estimate_dpa)) * 1.1

  plot_data |>
    filter(!grepl("mIgG|ACTB", marker)) |>
    mutate(
      significance = case_when(
        p_adj_dpa < 0.05 & p_adj_daa < 0.05 ~ "Both significant",
        p_adj_dpa < 0.05 ~ "Polarity significant",
        p_adj_daa < 0.05 ~ "Abundance significant",
        TRUE ~ "Not significant"
      ),
      significance = factor(
        significance,
        levels = c(
          "Both significant",
          "Polarity significant",
          "Abundance significant",
          "Not significant"
        )
      ),
      avg_estimate = abs((estimate_dpa + estimate_daa) / 2),
      # Scale text size based on average estimate
      text_size = scales::rescale(avg_estimate, to = c(2, 8))
    ) |>
    ggplot(aes(x = estimate_dpa, y = estimate_daa, color = significance)) +
    geom_point(aes(size = avg_estimate * 0.5)) +
    ggpubr::theme_pubr(base_size = 20) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    ggrepel::geom_text_repel(
      data = . %>% filter(significance != "Not significant"),
      aes(label = marker, size = text_size),
      nudge_y = 0.05,
      max.overlaps = 15
    ) +
    scale_color_manual(
      values = c(
        "Both significant" = "#E41A1C",
        "Polarity significant" = "#377EB8",
        "Abundance significant" = "#4DAF4A",
        "Not significant" = "grey70"
      )
    ) +
    guides(
      size = "none",
      color = guide_legend(override.aes = list(size = 5))
    ) + # Remove size from legend
    lims(x = c(-pol_max, pol_max), y = c(-abund_max, abund_max)) +
    labs(
      x = "Polarity difference",
      y = "Abundance difference",
      color = ""
    )
}

# Save output figures
export_plot <- function(
  filename,
  plot = ggplot2::last_plot(),
  width = 10,
  height = 10,
  create_dir = TRUE,
  file_formats = getOption("export_plot.file_formats", c("png", "pdf")),
  overwrite = getOption("export_plot.overwrite", TRUE)
) {
  # assert_single_value(filename, type = "string")
  # assert_class(plot, "ggplot")
  # assert_single_value(width, type = "numeric")
  # assert_single_value(height, type = "numeric")
  # assert_single_value(create_dir, type = "bool")
  # assert_vector(file_formats, type = "character", n = 1)
  # assert_single_value(overwrite, type = "bool")
  ext <- fs::path_ext(filename)
  if (nchar(ext) > 0) {
    file_formats <- ext
    filename <- fs::path_ext_remove(filename)
  }
  for (format in file_formats) {
    out_file <- paste0(filename, ".", format)
    if (file.exists(out_file) && !overwrite) {
      rlang::warn(glue::glue("Skipping existing file: {out_file}"))
      next
    }
    ggsave(
      filename = out_file,
      plot = plot,
      width = width,
      height = height,
      create.dir = create_dir,
      limitsize = FALSE
    )
  }
  invisible(NULL)
}
