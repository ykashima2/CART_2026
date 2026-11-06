lollipop_plot <- function(data,
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
                          color_palette = NULL) {
    # Input validation
    consistency <- match.arg(consistency)
    required_cols <- c("marker", condition_col, value_col)
    missing_cols <- setdiff(required_cols, colnames(data))
    if (length(missing_cols) > 0) {
        stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
    }

    # Check if p_adj column exists
    has_p_values <- "p_adj" %in% colnames(data)

    # Calculate consistency for each marker with the threshold condition
    if (has_p_values) {
        marker_consistency <- data |>
            group_by(marker) |>
            summarize(
                n_significant = sum(p_adj < p_adj_threshold),
                n_conditions = n(),
                has_different_signs = min(.data[[value_col]]) * max(.data[[value_col]]) < 0,
                has_small_effect = any(abs(.data[[value_col]]) < threshold),
                consistent = if_else(
                    n_significant == n_conditions &&
                        !has_different_signs &&
                        !has_small_effect,
                    "consistent",
                    "inconsistent"
                ),
                max_abs_difference = max(abs(.data[[value_col]]))
            ) |>
            select(-has_different_signs, -has_small_effect, -n_significant, -n_conditions)
    } else {
        marker_consistency <- data |>
            group_by(marker) |>
            summarize(
                has_different_signs = min(.data[[value_col]]) * max(.data[[value_col]]) < 0,
                has_small_effect = any(abs(.data[[value_col]]) < threshold),
                consistent = if_else(
                    !has_different_signs &&
                        !has_small_effect,
                    "consistent",
                    "inconsistent"
                ),
                max_abs_difference = max(abs(.data[[value_col]]))
            ) |>
            select(-has_different_signs, -has_small_effect)
    }

    if (consistency != "all") {
        marker_consistency <- marker_consistency |>
            filter(consistent == consistency)
    }

    if (nrow(marker_consistency) == 0) {
        stop("No markers found matching the specified criteria")
    }

    top_markers <- marker_consistency |>
        top_n(top_n, max_abs_difference) |>
        arrange(desc(max_abs_difference)) |>
        pull(marker)

    # Add significance information to plot data
    plot_data <- data |>
        filter(marker %in% top_markers) |>
        mutate(
            marker = factor(marker, levels = rev(top_markers)),
            condition = factor(.data[[condition_col]]),
            point_alpha = if (has_p_values) {
                if_else(p_adj < p_adj_threshold, 1, 0.2)
            } else {
                1 # If no p-values, all points are fully opaque
            }
        ) |>
        filter(!grepl("^mIgG|ACTB", marker))

    # Create subtitle based on consistency selection
    subtitle <- switch(consistency,
        "all" = "Top markers (abs. difference)",
        "consistent" = sprintf("Top consistent markers (abs. difference ≥ %.2f)", threshold),
        "inconsistent" = sprintf("Top inconsistent markers (abs. difference < %.2f or mixed signs)", threshold)
    )

    # Build plot
    p <- ggplot(plot_data, aes(x = marker, y = .data[[value_col]], color = condition))

    if (show_zero_line) {
        p <- p + geom_hline(yintercept = 0, linetype = "dashed", color = "gray50")
    }

    p <- p +
        geom_point(aes(alpha = point_alpha),
            position = position_dodge(width = 0.5),
            size = point_size
        ) +
        geom_linerange(aes(ymin = 0, ymax = .data[[value_col]], alpha = point_alpha),
            position = position_dodge(width = 0.5),
            size = line_size
        ) +
        guides(alpha = "none") +
        coord_flip()

    # Add color scale based on whether color_palette is provided
    p <- p + if (is.null(color_palette)) {
        scale_color_brewer(palette = "Set1", name = condition_col)
    } else {
        # Verify that all conditions have colors
        conditions <- unique(plot_data$condition)
        missing_conditions <- setdiff(conditions, names(color_palette))
        if (length(missing_conditions) > 0) {
            stop("Missing colors for conditions: ", paste(missing_conditions, collapse = ", "))
        }
        scale_color_manual(values = color_palette, name = condition_col)
    }

    p <- p +
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
