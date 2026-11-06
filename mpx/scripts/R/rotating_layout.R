# library(magrittr)

rotating_layout_old <- function(
    data,
    file,
    cell_col = NULL,
    marker_col = NULL,
    pt_opacity = 0.7,
    pt_size = 1,
    colors = RColorBrewer::brewer.pal(9, "Blues"),
    max_degree = 360,
    center_zero = FALSE,
    frames = 500,
    pad = 0.1,
    show_first_frame = TRUE,
    width = 540,
    height = 500,
    res = 150,
    delay = 1 / 20,
    ggplot_theme = NULL,
    title = "",
    bg = "white",
    label_grid_axes = TRUE,
    margin_widths = c(0.1, 0.1),
    use_facet_grid = FALSE,
    flip = FALSE,
    graphics_use = c("ggplot2", "base"),
    boomerang = FALSE,
    cl = NULL,
    keep_frames = FALSE
) {
    if (fs::path_ext(file) == "gif") {
        rlang::check_installed("gifski")
    } else {
        rlang::check_installed("av")
    }

    # ragg is slightly faster than png, but needs to be installed
    if (!requireNamespace("ragg")) {
        warn(glue("For faster rendering, install the 'ragg' package."))
        dev_png <- png
    } else {
        dev_png <- ragg::agg_png
    }

    graphics_use <- match.arg(graphics_use, choices = c("ggplot2", "base"))
    if (graphics_use == "base") {
        cli_alert_warning(
            glue(
                "The following parameters are not supported with graphics_use = 'base':\n",
                "ggplot_theme, label_grid_axes, title\n",
                "Titles, text annotations and color bar will be missing from the plot."
            )
        )
    }

    .validate_render_rotating_layout_input_params(
        data,
        marker_col,
        cell_col,
        pt_opacity,
        pt_size,
        colors,
        max_degree,
        center_zero,
        frames,
        pad,
        show_first_frame,
        width,
        height,
        res,
        delay,
        ggplot_theme,
        title,
        bg,
        label_grid_axes,
        margin_widths,
        use_facet_grid,
        flip,
        boomerang
    )

    # Set variables if NULL
    cell_col <- cell_col %||% "cell"
    marker_col <- marker_col %||% "marker"
    # Add missing columns to data required for downstream processing
    if (!cell_col %in% names(data)) {
        data <- data %>%
            mutate(!!sym(cell_col) := cell_col)
    }
    if (!marker_col %in% names(data)) {
        data <- data %>%
            mutate(!!sym(marker_col) := marker_col)
    }

    # Create a vector of angles for each frame
    angles <- seq(0, max_degree, length.out = frames)

    # Create chunks to process
    indices_split <- seq_along(angles) %>% split(ceiling(seq_along(angles) / 3))

    # Create temporary directory where the frames will be saved
    tmp_dir <- fs::file_temp() %>% stringr::str_replace("file", "dir")
    dir.create(tmp_dir)

    # Get cell axes ranges
    # First, we group the tibble by the cell column, then we
    # split the tibble into a list of tibbles, one for each cell.
    # We then rescale the x, y, and z coordinates of each cell
    # such that the maximum absolute value of the coordinates is 1.
    # This is to ensure that all cells will fit to the same plot area.
    xyz_grouped <- data %>%
        select(x, y, z, node_val, all_of(c(marker_col, cell_col))) %>%
        group_by(across(all_of(cell_col)))
    group_d <- xyz_grouped %>% group_keys()
    xyz_list <- xyz_grouped %>%
        group_split() %>%
        `names<-`(
            if (nrow(group_d) == 1) "cell" else as.character(group_d[[1]])
        )
    xyz_list <- lapply(xyz_list, function(xyz) {
        xyz_scaled <- xyz %>%
            select(x, y, z) %>%
            scale_layout(percentile = 1)
        xyz[, c("x", "y", "z")] <- xyz_scaled
        return(xyz)
    })

    # Set the axis limits
    padded_max_lim <- 1 + 2 * pad
    xyz_limits <- c(-padded_max_lim, padded_max_lim)

    if (inherits(data$node_val, "numeric")) {
        # Get marker value limits
        # First, we group the tibble by the marker column, then we
        # split the tibble into a list of tibbles, one for each marker.
        # We then summarize the minimum and maximum values of the node_val
        # column for each marker. This is to make sure that each marker
        # has the same fixed color scale across different cells.
        marker_limits <- data %>%
            group_by(!!sym(marker_col)) %>%
            summarize(
                min = min(node_val, na.rm = TRUE),
                max = max(node_val, na.rm = TRUE),
                .groups = "keep"
            )
        marker_limits <- marker_limits %>%
            group_split() %>%
            `names<-`(marker_limits[[1]]) %>%
            lapply(function(x) {
                x %>%
                    select(min, max) %>%
                    as.numeric()
            })
    } else {
        marker_limits <- NULL
    }

    # Create a nested list of tibbles, one for each cell, and one for each marker
    xyz_list_nested <- lapply(xyz_list, function(xyz) {
        xyz_list_cell <- xyz %>% group_by(across(all_of(marker_col)))
        group_d <- xyz_list_cell %>% group_keys()
        xyz_list_cell <- xyz_list_cell %>%
            group_split() %>%
            `names<-`(
                if (nrow(group_d) == 1) "marker" else as.character(group_d[[1]])
            )
        return(xyz_list_cell)
    })

    # Show the first frame if requested
    if (show_first_frame && interactive()) {
        .render_first_frame(
            xyz_list_nested,
            marker_limits,
            pt_opacity,
            pt_size,
            marker_col,
            cell_col,
            colors,
            center_zero,
            xyz_limits,
            tmp_dir,
            width,
            height,
            res,
            ggplot_theme,
            title,
            bg,
            use_facet_grid,
            flip,
            label_grid_axes,
            margin_widths,
            graphics_use,
            dev_png
        )
    }

    cli::cli_h3("Rendering frames...")
    # Loop over the indices chunks and create the frames
    dump <- pbapply::pblapply(
        indices_split,
        function(indices) {
            for (i in indices) {
                # Apply rotation and arrange values by the y axis
                # such that the nodes are stacked from back to front
                xyz_list_cur <- lapply(xyz_list_nested, function(xyz_cell) {
                    xyz_cell <- lapply(xyz_cell, function(xyz) {
                        xyz_rotated <- .rotate_around_z(
                            xyz %>% select(x, y, z),
                            angles[i]
                        )
                        xyz[, c("x", "y", "z")] <- xyz_rotated
                        xyz <- xyz %>%
                            arrange(pick(all_of("y")))
                        return(xyz)
                    })
                    return(xyz_cell)
                })

                # Current frame PNG file
                file <- file.path(
                    tmp_dir,
                    paste0(
                        "plot_",
                        stringr::str_pad(i, side = "left", width = 4, "0"),
                        ".png"
                    )
                )

                # Create and export plot using ggplot2
                if (graphics_use == "ggplot2") {
                    p <- .scatter_3d_plot(
                        xyz_list_cur,
                        cell_col,
                        marker_col,
                        marker_limits,
                        colors,
                        center_zero,
                        xyz_limits,
                        pt_opacity,
                        pt_size,
                        use_facet_grid,
                        flip,
                        label_grid_axes,
                        margin_widths
                    )
                    # Add ggplot theme
                    if (!is.null(ggplot_theme)) {
                        p <- p & ggplot_theme
                    }
                    if (title != "") {
                        p <- p + plot_annotation(title = title)
                    }
                    ggsave(
                        file,
                        p,
                        width = width,
                        height = height,
                        dpi = res,
                        device = dev_png,
                        units = "px",
                        bg = bg
                    )
                }

                # Create and export plot using base R graphics
                if (graphics_use == "base") {
                    dev_png(
                        file,
                        width = width,
                        height = height,
                        res = res,
                        units = "px",
                        bg = bg
                    )
                    .scatter_3d_plot_base_R(
                        xyz_list_cur,
                        marker_limits,
                        colors,
                        center_zero,
                        xyz_limits,
                        pt_opacity,
                        pt_size,
                        bg,
                        flip
                    )
                    dev.off()
                }
            }
        },
        cl = cl
    )

    # Create video
    cli::cli_h3("Creating video...")
    png_files <- fs::dir_ls(tmp_dir, glob = "*.png")
    if (boomerang) {
        png_files <- c(png_files, rev(png_files))
    }

    if (fs::path_ext(file) == "gif") {
        # Use gifski if the file extension is .gif
        out_path <- gifski::gifski(
            png_files,
            file,
            width = width,
            height = height,
            delay = delay
        )
    } else {
        # Use av if the file extension is not .gif
        try_message <- try(
            {
                av::av_encode_video(png_files, file, framerate = 1 / delay)
            },
            silent = TRUE
        )
        # If the video creation fails, the most likely reason is that
        # the format is not supported. In this case, we provide a warning
        # and return the path to the directory containing the frames.
        if (inherits(try_message, "try-error")) {
            cli::cli_alert_danger(
                glue(
                    "{try_message}\n",
                    "Failed to generate the video. ",
                    "Make sure to use an appropriate file extension. \n",
                    "The frames can be found in:\n {tmp_dir}"
                )
            )
            return(invisible(NULL))
        }
    }

    # Remove temporary directory to clean up
    if (!keep_frames) {
        try_message <- try(
            {
                fs::dir_delete(tmp_dir)
            },
            silent = TRUE
        )
        if (inherits(try_message, "try-error")) {
            warn(glue("Failed to remove temporary directory:\n {tmp_dir}"))
        }
    } else {
        cli::cli_alert_info(glue("The frames can be found in:\n {tmp_dir}"))
    }

    cli::cli_alert_success(glue("The video has been saved to:\n {file}"))
}


.scatter_3d_plot <- function(
    xyz_list_nested,
    cell_col,
    marker_col,
    marker_limits,
    colors,
    center_zero,
    xyz_limits,
    pt_opacity,
    pt_size,
    use_facet_grid = FALSE,
    flip = FALSE,
    label_grid_axes = TRUE,
    margin_widths = c(0.1, 0.1)
) {
    if (use_facet_grid) {
        xyz_collapsed <- lapply(xyz_list_nested, function(xyz_cell) {
            xyz_cell %>%
                unname() %>%
                bind_rows()
        }) %>%
            unname() %>%
            bind_rows()
        p <- ggplot(xyz_collapsed, aes(x, z, size = y, color = node_val)) +
            geom_point(alpha = pt_opacity) +
            scale_size(range = c(0, pt_size), limits = xyz_limits) +
            coord_fixed() +
            theme_void() +
            guides(size = "none") +
            scale_x_continuous(limits = xyz_limits, expand = expansion()) +
            scale_y_continuous(limits = xyz_limits, expand = expansion()) +
            theme(
                text = element_text(color = "white"),
                strip.text = element_text(color = "white"),
                legend.text = element_text(color = "white"),
                legend.title = element_text(color = "white")
            ) +
            {
                if (inherits(xyz_collapsed$node_val, "numeric")) {
                    max_abs_val <- max(abs(xyz_collapsed$node_val))
                    lims <- if (center_zero) {
                        c(-max_abs_val, max_abs_val)
                    } else {
                        c(0, max_abs_val)
                    }
                    scale_color_gradientn(
                        colours = colors,
                        limits = lims,
                        oob = scales::squish
                    )
                } else {
                    scale_color_manual(values = colors)
                }
            } +
            {
                if (flip) {
                    facet_grid(reformulate(marker_col, cell_col), switch = "y")
                } else {
                    facet_grid(reformulate(cell_col, marker_col), switch = "y")
                }
            }
    }

    if (!use_facet_grid) {
        # Create a list of ggplot objects per cell
        cell_plots <- lapply(names(xyz_list_nested), function(cell_id) {
            xyz_list_cell <- xyz_list_nested[[cell_id]]
            # Create a list of ggplot objects per marker
            marker_plots <- lapply(names(xyz_list_cell), function(marker_id) {
                p <- ggplot(
                    xyz_list_cell[[marker_id]] %>%
                        mutate(row_text = marker_id),
                    aes(x, z, size = y, color = node_val)
                ) +
                    geom_point(alpha = pt_opacity) +
                    scale_size(range = c(0, pt_size), limits = xyz_limits) +
                    coord_fixed() +
                    theme_void() +
                    guides(size = "none") +
                    scale_x_continuous(
                        limits = xyz_limits,
                        expand = expansion()
                    ) +
                    scale_y_continuous(
                        limits = xyz_limits,
                        expand = expansion()
                    ) +
                    theme(
                        text = element_text(color = "white"),
                        strip.text = element_text(color = "white"),
                        legend.text = element_text(color = "white"),
                        legend.title = element_text(color = "white")
                    ) +
                    {
                        if (
                            inherits(
                                xyz_list_cell[[marker_id]]$node_val,
                                "numeric"
                            )
                        ) {
                            max_abs_val <- max(abs(
                                xyz_list_cell[[marker_id]]$node_val
                            ))
                            lims <- if (center_zero) {
                                c(-max_abs_val, max_abs_val)
                            } else {
                                c(0, max_abs_val)
                            }
                            scale_color_gradientn(
                                colours = colors,
                                limits = lims,
                                oob = scales::squish,
                                na.value = colors[1]
                            )
                        } else {
                            scale_color_manual(values = colors)
                        }
                    }
                return(p)
            })

            if (flip) {
                p_markers <- wrap_plots(marker_plots, nrow = 1)
            } else {
                p_markers <- wrap_plots(marker_plots, ncol = 1)
            }

            return(p_markers)
        })

        # Wrap marker plots along rows or columns
        if (flip) {
            p <- wrap_plots(cell_plots, ncol = 1) &
                theme(legend.position = "right")
        } else {
            p <- wrap_plots(cell_plots, nrow = 1)
        }
        # Collect guides such that the color legend only appear on
        # the right side of the plot
        p <- p + plot_layout(guides = "collect")

        # row and col labels
        if (
            label_grid_axes &&
                ((length(xyz_list_nested) > 1) ||
                    (length(xyz_list_nested[[1]]) > 1))
        ) {
            # Create lists of row and column labels (grobs)
            # that can be added to the final patchwork
            marker_labels <- lapply(names(xyz_list_nested[[1]]), function(lbl) {
                wrap_elements(
                    grid::textGrob(
                        label = lbl,
                        gp = grid::gpar(col = "white")
                    )
                ) +
                    theme(
                        plot.background = element_rect(
                            fill = "transparent",
                            colour = NA
                        )
                    )
            })
            cell_labels <- lapply(names(xyz_list_nested), function(lbl) {
                wrap_elements(
                    grid::textGrob(
                        label = lbl,
                        gp = grid::gpar(col = "white")
                    )
                ) +
                    theme(
                        plot.background = element_rect(
                            fill = "transparent",
                            colour = NA
                        )
                    )
            })

            if (!flip) {
                col_labels <- cell_labels
                row_labels <- marker_labels
            } else {
                col_labels <- marker_labels
                row_labels <- cell_labels
            }

            col_width <- (1 - margin_widths[2]) / length(col_labels)
            row_width <- (1 - margin_widths[1]) / length(row_labels)

            # Create a patchwork with the column labels, leaving an empty
            # space for the upper left corner of the grid
            col_labels <- c(list(plot_spacer()), col_labels) %>%
                wrap_plots(
                    nrow = 1,
                    widths = c(
                        margin_widths[2],
                        rep(col_width, length(col_labels))
                    )
                )

            row_labels <- wrap_plots(row_labels, ncol = 1)

            # Add row labels to the patchwork on the left side
            p <- wrap_plots(
                row_labels,
                p,
                widths = c(margin_widths[1], 1 - margin_widths[1])
            )
            # Add column labels to the patchwork on the top
            p <- wrap_plots(
                col_labels,
                p,
                ncol = 1,
                heights = c(margin_widths[2], 1 - margin_widths[2])
            )
        }
    }

    # Set background color to transparent so that the
    # background color for the frames can be controlled with
    # the png device
    p <- p &
        theme(
            plot.background = element_rect(fill = "transparent", color = NA),
            text = element_text(color = "white"),
            strip.text = element_text(color = "white"),
            legend.text = element_text(color = "white"),
            legend.title = element_text(color = "white")
        )

    return(p)
}


.scatter_3d_plot_base_R <- function(
    xyz_list_nested,
    marker_limits,
    colors,
    center_zero,
    xyz_limits,
    pt_opacity,
    pt_size,
    bg,
    flip = FALSE
) {
    # Save the current par settings
    # and reset when the function exits
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))

    if (flip) {
        ar <- c(length(xyz_list_nested[[1]]), length(xyz_list_nested))
        lyt <- matrix(1:(ar[1] * ar[2]), nrow = ar[2], byrow = FALSE)
    } else {
        ar <- c(length(xyz_list_nested), length(xyz_list_nested[[1]]))
        lyt <- matrix(1:(ar[1] * ar[2]), nrow = ar[1], byrow = TRUE)
    }

    # Create a grid for plotting
    par(
        mar = c(0, 0, 0, 0), # Remove margins
        pty = "s", # square aspect ratio <-> coord_fixed()
        bg = ifelse(bg == "transparent", NA, bg)
    )

    layout(mat = lyt)

    for (marker_id in names(xyz_list_nested[[1]])) {
        for (cell_id in names(xyz_list_nested)) {
            df <- xyz_list_nested[[cell_id]][[marker_id]]

            max_radius <- max(sqrt(rowSums((df %>% select(x, y, z))^2)))

            z_norm <- scales::rescale(
                df$z,
                from = c(-max_radius, max_radius),
                to = c(3, 1)
            )
            apparent_sizes <- sqrt(1 / (z_norm^2))

            # Define node colors based on the type of node_val
            if (inherits(df$node_val, "numeric")) {
                if (center_zero) {
                    max_abs_node_val <- max(abs(df$node_val))
                    cols <- scales::col_numeric(
                        domain = c(-max_abs_node_val, max_abs_node_val),
                        palette = colors
                    )(df$node_val)
                } else {
                    cols <- scales::col_numeric(
                        domain = c(0, max(abs(df$node_val))),
                        palette = colors
                    )(df$node_val)
                }
            } else {
                if (inherits(df$node_val, "character")) {
                    df$node_val <- as.factor(df$node_val)
                }
                cols <- scales::col_factor(
                    domain = levels(df$node_val),
                    palette = colors
                )(df$node_val)
            }

            # Add transparency to the colors
            cols <- scales::alpha(cols, alpha = pt_opacity)

            # Scale the apparent size of the points based on the y-axis values
            # plot() uses the cex parameter to scale the size of the points
            # which is based on the radius of the circle.
            # First, we need to normalize the y-axis values to a range of 3 to 1,
            # meaning that the point closest to the camera is 1 distance away and the
            # point furthest away is 3 distances away.
            y_norm <- scales::rescale(
                df$y,
                from = c(-max_radius, max_radius),
                to = c(3, 1)
            )
            # Then we use the inverse square law to adjust the area of the points
            # and take the square root to get the radii. By multiplying the radii
            # by pt_size, we can adjust the size of the points while retaining the
            # perceived size differences.
            apparent_sizes <- sqrt(1 / (y_norm^2)) * pt_size

            # TODO: Add column and ro labels
            plot(
                df$x,
                df$z,
                cex = apparent_sizes,
                pch = 16,
                col = cols,
                xlab = "",
                ylab = "",
                axes = FALSE,
                xlim = xyz_limits,
                ylim = xyz_limits,
            )
        }
    }
}


.create_rotation_matrix <- function(angle_degrees) {
    angle_radians <- angle_degrees * pi / 180
    matrix(
        c(
            cos(angle_radians),
            -sin(angle_radians),
            0,
            sin(angle_radians),
            cos(angle_radians),
            0,
            0,
            0,
            1
        ),
        nrow = 3,
        byrow = TRUE
    )
}


.rotate_around_z <- function(coords, angle_degrees) {
    rotation_matrix <- .create_rotation_matrix(angle_degrees)
    as_tibble(
        t(rotation_matrix %*% t(coords %>% as.matrix())),
        .name_repair = ~ c("x", "y", "z")
    )
}


scale_layout <- function(
    layout,
    percentile = 0.95
) {
    stopifnot(
        "'layout' must be a data.frame like object" = inherits(
            layout,
            c("tbl_df", "data.frame")
        ),
        "'percentile' must be a numeric value between 0 and 1" = is.numeric(
            percentile
        ) &&
            between(percentile, 0, 1)
    )

    layout <- layout %>%
        select(any_of(c("x", "y", "z")))

    if (!ncol(layout) %in% c(2, 3)) {
        abort("The layout must have 2 or 3 columns.")
    }

    radii <- sqrt(rowSums(layout^2))
    perc <- quantile(radii, percentile)
    layout_scaled <- layout %>%
        mutate(across(everything(), ~ .x / perc))

    return(layout_scaled)
}


.render_first_frame <- function(
    xyz_list_nested,
    marker_limits,
    pt_opacity,
    pt_size,
    marker_col,
    cell_col,
    colors,
    center_zero,
    xyz_limits,
    tmp_dir,
    width,
    height,
    res,
    ggplot_theme,
    title,
    bg,
    use_facet_grid,
    flip,
    label_grid_axes,
    margin_widths,
    graphics_use,
    dev_png
) {
    rlang::check_installed("png")
    file <- file.path(tmp_dir, "plot_0000.png")

    # Use the more flexible but slower ggplot2 graphics
    if (graphics_use == "ggplot2") {
        p <- .scatter_3d_plot(
            xyz_list_nested,
            cell_col,
            marker_col,
            marker_limits,
            colors,
            center_zero,
            xyz_limits,
            pt_opacity,
            pt_size,
            use_facet_grid,
            flip,
            label_grid_axes,
            margin_widths
        )

        # Add theme
        if (!is.null(ggplot_theme)) {
            p <- p & ggplot_theme
        }
        if (title != "") {
            p <- p + plot_annotation(title = title)
        }

        # Export first frame
        ggsave(
            file,
            p,
            width = width,
            height = height,
            dpi = res,
            device = dev_png,
            units = "px",
            bg = bg
        )
    }

    # Use the faster but less flexible base graphics
    if (graphics_use == "base") {
        dev_png(
            file,
            width = width,
            height = height,
            res = res,
            units = "px",
            bg = bg
        )
        .scatter_3d_plot_base_R(
            xyz_list_nested,
            marker_limits,
            colors,
            center_zero,
            xyz_limits,
            pt_opacity,
            pt_size,
            bg,
            flip
        )
        dev.off()
    }

    # Display the first frame in the viewer
    p <- png::readPNG(file) %>% as.raster()
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))
    par(mar = c(0, 0, 0, 0))
    plot(p)
    cli_alert(glue(
        "The first frame has been rendered in ",
        "the viewer. Do you want to continue?"
    ))

    # Ask user if they want to continue
    if (utils::menu(c("Yes", "No")) != 1) {
        abort("Aborted by user", call = caller_env())
    }

    dev.off()
    fs::file_delete(file)
}


.areColors <- function(colors) {
    sapply(colors, function(X) {
        tryCatch(is.matrix(col2rgb(X)), error = function(e) FALSE)
    })
}


.validate_render_rotating_layout_input_params <- function(
    data,
    marker_col,
    cell_col,
    pt_opacity,
    pt_size,
    colors,
    max_degree,
    center_zero,
    frames,
    pad,
    show_first_frame,
    width,
    height,
    res,
    delay,
    ggplot_theme,
    title,
    bg,
    label_grid_axes,
    margin_widths,
    use_facet_grid,
    flip,
    boomerang
) {
    stopifnot(
        "'data' must be a 'tbl_df' object" = inherits(data, "tbl_df"),
        "'pt_opacity' must be a numeric value between 0 and 1" = is.numeric(
            pt_opacity
        ) &&
            between(pt_opacity, 0, 1),
        "'pt_size' must be a positive numeric value between 0 and 5" = is.numeric(
            pt_size
        ) &&
            between(pt_size, 0, 5),
        "'colors' must be a vector of color names" = inherits(
            colors,
            "character"
        ),
        "'colors' can only contain valid color strings" = all(.areColors(
            colors
        )),
        "'max_degree' must be a numeric value between 90 and 360" = is.numeric(
            max_degree
        ) &&
            between(max_degree, 90, 360),
        "'center_zero' must be either TRUE or FALSE" = is.logical(center_zero),
        "'frames' must be a positive numeric value" = is.numeric(frames) &&
            frames > 0,
        "'pad' must be a numeric value between 0 and 1" = is.numeric(pad) &&
            between(pad, 0, 1),
        "'show_first_frame' must be either TRUE or FALSE" = is.logical(
            show_first_frame
        ),
        "'width' must be a positive numeric value" = is.numeric(width) &&
            width > 0,
        "'height' must be a positive numeric value" = is.numeric(height) &&
            height > 0,
        "'res' must be a positive numeric value" = is.numeric(res) && res > 0,
        "'delay' must be a value between 0.01 and 1" = is.numeric(delay) &&
            between(delay, 0.01, 1),
        "'title' must be a character string" = is.character(title),
        "'bg' must be a valid color" = length(bg) == 1 && .areColors(bg),
        "'label_grid_axes' must be either TRUE or FALSE" = is.logical(
            label_grid_axes
        ),
        "'margin_widths' must be a numeric vector of length 2 with values between 0 and 1" = is.numeric(
            margin_widths
        ) &&
            length(margin_widths) == 2 &&
            all(between(margin_widths, 0, 1)),
        "'use_facet_grid' must be either TRUE or FALSE" = is.logical(
            use_facet_grid
        ),
        "'flip' must be either TRUE or FALSE" = is.logical(flip),
        "'boomerang' must be either TRUE or FALSE" = is.logical(boomerang)
    )

    missing_cols <- setdiff(c("x", "y", "x", "node_val"), names(data))
    if (length(missing_cols) > 0) {
        abort(glue("Missing column(s): {paste(missing_cols, collapse = ', ')}"))
    }

    # Check marker_col and cell_col
    if (!is.null(marker_col)) {
        stopifnot(
            "'marker_col' must be a character string" = is.character(
                marker_col
            ),
            "'marker_col' must be a column name in 'data'" = marker_col %in%
                names(data)
        )
    }
    if (!is.null(cell_col)) {
        stopifnot(
            "'cell_col' must be a character string" = is.character(cell_col),
            "'cell_col' must be a column name in 'data'" = cell_col %in%
                names(data)
        )
    }
    if (!is.null(marker_col) && !is.null(cell_col)) {
        stopifnot(
            "'marker_col' and 'cell_col' cannot be the same" = marker_col !=
                cell_col
        )
    }

    # Check column types
    col_classes <- sapply(data[, c("x", "y", "z")], class)
    if (any(col_classes != "numeric")) {
        abort("Columns 'x', 'y', and 'z' must be numeric")
    }
    if (inherits(data$node_val, c("character", "factor"))) {
        if (length(colors) < length(unique(data$node_val))) {
            warn(
                "Not enough colors for node_val groups. Expanding colors vector."
            )
            colors <- colorRampPalette(colors)(length(unique(data$node_val)))
        }
        if (inherits(data$node_val, "character")) {
            data <- data %>%
                mutate(node_val = factor(node_val, levels = unique(node_val)))
        }
    }
    if (inherits(data$node_val, "numeric")) {
        if (any(is.na(data$node_val))) {
            abort("Column 'node_val' cannot contain missing values")
        }
    }
    if (!is.null(ggplot_theme)) {
        if (!inherits(ggplot_theme, "theme")) {
            abort("Argument 'ggplot_theme' must be a 'theme' object")
        }
    }
}


## Static version =====================================================
render_static_layout <- function(
    data,
    file = NULL,
    pt_opacity = 0.7,
    pt_size = 1,
    colors = c("grey30", "#3c4c5f", "#f5c96a"),
    angle = 0,
    center_zero = FALSE,
    pad = 0.1,
    width = 540,
    height = 500,
    res = 150,
    ggplot_theme = NULL,
    title = "",
    bg = "white"
) {
    # Validate input data structure
    required_cols <- c("x", "y", "z", "marker", "node_val", "component")
    if (!all(required_cols %in% names(data))) {
        missing <- setdiff(required_cols, names(data))
        stop(sprintf(
            "Missing required columns: %s",
            paste(missing, collapse = ", ")
        ))
    }

    # Group data by component and marker
    xyz_grouped <- data |>
        group_by(component, marker)

    # Scale coordinates for each component
    xyz_list <- xyz_grouped |>
        group_split() |>
        lapply(function(xyz) {
            xyz_scaled <- xyz |>
                select(x, y, z) |>
                scale_layout(percentile = 1)
            xyz[, c("x", "y", "z")] <- xyz_scaled
            return(xyz)
        })

    # Set axis limits with padding
    padded_max_lim <- 1 + 2 * pad
    xyz_limits <- c(-padded_max_lim, padded_max_lim)

    # Apply rotation to coordinates
    xyz_rotated <- lapply(xyz_list, function(xyz) {
        rotated <- .rotate_around_z(xyz |> select(x, y, z), angle)
        xyz[, c("x", "y", "z")] <- rotated
        xyz |> arrange(y) # Order by depth
    })

    # Combine all components into one plot
    plot_data <- bind_rows(xyz_rotated)

    # Create plot
    p <- ggplot(plot_data, aes(x = x, y = z, size = -y, color = node_val)) +
        geom_point(alpha = pt_opacity) +
        scale_size(range = c(0, pt_size), guide = "none") +
        scale_color_gradientn(
            colors = colors,
            limits = if (center_zero) {
                max_abs <- max(abs(data$node_val))
                c(-max_abs, max_abs)
            } else {
                NULL
            },
            oob = scales::squish
        ) +
        coord_fixed() +
        facet_wrap(~marker, nrow = 1) +
        theme_void() +
        theme(
            strip.text = element_text(face = "bold"),
            legend.position = "right"
        )

    # Add custom theme if provided
    if (!is.null(ggplot_theme)) {
        p <- p + ggplot_theme
    }

    # Add title if provided
    if (title != "") {
        p <- p + ggtitle(title)
    }

    # Save plot if file path provided
    if (!is.null(file)) {
        ggsave(
            file,
            p,
            width = width,
            height = height,
            dpi = res,
            units = "px",
            bg = bg
        )
    }

    return(p)
}

# Add y-axis rotation matrix function
.create_y_rotation_matrix <- function(angle_degrees) {
    angle_rad <- angle_degrees * pi / 180
    matrix(
        c(
            cos(angle_rad),
            0,
            sin(angle_rad),
            0,
            1,
            0,
            -sin(angle_rad),
            0,
            cos(angle_rad)
        ),
        nrow = 3,
        byrow = TRUE
    )
}

# Add y-axis rotation function
.rotate_around_y <- function(coords, angle_degrees) {
    rotation_matrix <- .create_y_rotation_matrix(angle_degrees)
    as_tibble(
        t(rotation_matrix %*% t(coords |> as.matrix())),
        .name_repair = ~ c("x", "y", "z")
    )
}

.trim_quantiles <- function(
    x,
    bottom_q = 0,
    top_q = 0.99
) {
    # Handle matrix input
    if (is.matrix(x)) {
        # Apply trimming to each column separately
        for (i in 1:ncol(x)) {
            col_data <- x[, i]

            # Calculate quantiles for this column
            low_thr <- quantile(col_data, probs = bottom_q)
            high_thr <- quantile(col_data, probs = top_q)

            # Trim data for this column
            x[col_data < low_thr, i] <- low_thr
            x[col_data > high_thr, i] <- high_thr
        }
    } else {
        # Original behavior for numeric vectors
        low_thr <- quantile(x, probs = bottom_q)
        high_thr <- quantile(x, probs = top_q)

        x[x < low_thr] <- low_thr
        x[x > high_thr] <- high_thr
    }

    return(x)
}
