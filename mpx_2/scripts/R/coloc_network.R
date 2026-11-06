#' Colocalization network plots
#'
#' Helpers to build the network figures of marker colocalization (mean pearson_z)
#' for a set of markers in one cell type. `PlotColocNetwork()` runs the full
#' workflow and exports two figures:
#'   - `{file_stem}`: all samples as facets of one ggraph network
#'   - `{file_stem}_styled`: one row per donor (t0 -> t1) with a shared legend

default_network_gradient <- c(
  "#1F395F",
  "#44628E",
  "#718BB2",
  "#A8B9D1",
  "#DBE1EA",
  "#F1EEE9",
  "#F0D7E0",
  "#E3A6B8",
  "#CB6E8B",
  "#A23F5E",
  "#781534"
)

#' Summarize colocalization scores per sample for a marker set
#'
#' @param data Seurat object with PNA/MPX data.
#' @param markers Markers to include (pairs between these markers are kept).
#' @param cell_type Cell type (value of the `celltype` metadata column).
#' @param samples Values of `group_col` (e.g. "donorD_t0") to include. NULL = all.
#' @param group_col Metadata column defining the samples.
#' @param scores Optional precomputed `ColocalizationScores()` result with
#'   columns marker_1, marker_2, pearson_z, celltype and `group_col`. Avoids
#'   recomputing the scores for every figure.
#' @param min_cells Minimum number of cells of `cell_type` in a sample for the
#'   sample to be included in the plots. Samples below it stay in the summary
#'   with `included = FALSE`.
#'
#' @return A tibble with the per-sample summary from `SummarizeProximityScores()`
#'   (mean_join_count_z, pct_detected, ...), the number of cells of `cell_type`
#'   per sample (`n_cells`) and the inclusion status (`included`, n_cells >=
#'   min_cells). When `group_col` is "sample_alias" (format donor_timepoint),
#'   `donor` and `timepoint` columns are added.
SummarizeColocNetworkData <- function(
  data,
  markers,
  cell_type,
  samples = NULL,
  group_col = "sample_alias",
  scores = NULL,
  min_cells = 1
) {
  if (is.null(scores)) {
    # Work on a copy so the caller's object is not modified
    data$component <- NULL
    scores <- ColocalizationScores(
      data,
      meta_data_columns = c(group_col, "celltype")
    )
  }

  proximity_data <- scores |>
    dplyr::filter(
      celltype == cell_type,
      marker_1 %in% markers,
      marker_2 %in% markers,
      marker_1 != marker_2
    )
  if (!is.null(samples)) {
    proximity_data <- proximity_data |>
      dplyr::filter(.data[[group_col]] %in% samples)
  }
  proximity_data <- dplyr::collect(proximity_data)
  if (nrow(proximity_data) == 0) {
    stop(
      "No colocalization scores left after filtering on cell type, markers and samples."
    )
  }

  summary_data <- SummarizeProximityScores(
    proximity_data |> dplyr::rename(join_count_z = pearson_z),
    group_vars = group_col,
    proximity_metric = "join_count_z"
  )

  # Cells of this cell type per sample, and inclusion status
  cell_counts <- data[[]] |>
    dplyr::filter(celltype == cell_type) |>
    dplyr::count(.data[[group_col]], name = "n_cells")
  # SummarizeProximityScores() returns its own n_cells (cells with scores for
  # these markers); replace it with the cell count from the metadata
  summary_data <- summary_data |>
    dplyr::select(-dplyr::any_of("n_cells")) |>
    dplyr::left_join(cell_counts, by = group_col) |>
    dplyr::mutate(
      n_cells = dplyr::coalesce(n_cells, 0L),
      included = n_cells >= min_cells
    )
  excluded <- summary_data |>
    dplyr::filter(!included) |>
    dplyr::distinct(.data[[group_col]], n_cells)
  if (nrow(excluded) > 0) {
    message(
      "Samples below min_cells (",
      min_cells,
      ") flagged as excluded: ",
      paste0(excluded[[group_col]], " (", excluded$n_cells, ")", collapse = ", ")
    )
  }

  if (group_col == "sample_alias") {
    summary_data <- summary_data |>
      tidyr::separate(
        sample_alias,
        into = c("donor", "timepoint"),
        sep = "_",
        remove = FALSE
      )
  }
  summary_data
}

#' Keep the samples flagged as included in a network summary
#'
#' @param min_cells Optional extra cutoff applied on `n_cells`, on top of the
#'   `included` column (NULL = use `included` only).
#' @noRd
keep_included_samples <- function(summary_data, min_cells = NULL) {
  keep <- summary_data$included
  if (!is.null(min_cells)) {
    keep <- keep & summary_data$n_cells >= min_cells
  }
  out <- summary_data[keep, ]
  if (nrow(out) == 0) {
    stop("No samples pass the minimum cell count.")
  }
  out
}

#' Network overview: all included samples as facets
#'
#' @param summary_data Output of `SummarizeColocNetworkData()`.
#' @inheritParams PlotColocNetwork
PlotColocNetworkOverview <- function(
  summary_data,
  cell_type,
  group_col = "sample_alias",
  min_pct_detected = 0.5,
  min_abs_z = 0,
  only_positive = FALSE,
  gradient_palette = default_network_gradient,
  min_cells = NULL
) {
  excluded_label <- summary_data |>
    dplyr::anti_join(
      keep_included_samples(summary_data, min_cells),
      by = group_col
    ) |>
    dplyr::pull(dplyr::all_of(group_col)) |>
    unique() |>
    sort() |>
    paste(collapse = ", ")
  summary_data <- keep_included_samples(summary_data, min_cells)

  df_net <- summary_data |>
    dplyr::filter(
      if (only_positive) {
        mean_join_count_z > min_abs_z
      } else {
        abs(mean_join_count_z) > min_abs_z
      }
    ) |>
    dplyr::filter(pct_detected > min_pct_detected) |>
    dplyr::mutate(from = marker_1, to = marker_2) |>
    dplyr::select(from, to, mean_join_count_z, dplyr::all_of(group_col))
  if (nrow(df_net) == 0) {
    stop("No edges pass the min_abs_z / min_pct_detected filters.")
  }

  max_z <- max(abs(df_net$mean_join_count_z), na.rm = TRUE)
  samples_label <- paste(sort(unique(df_net[[group_col]])), collapse = ", ")

  tidygraph::as_tbl_graph(df_net) |>
    ggraph::ggraph(layout = "kk") +
    ggraph::geom_edge_link(aes(
      color = mean_join_count_z,
      width = abs(mean_join_count_z) * 0.6
    )) +
    ggraph::geom_node_label(aes(label = name), size = 2.5) +
    ggraph::scale_edge_color_gradientn(
      colours = gradient_palette,
      limits = c(-max_z, max_z)
    ) +
    scale_x_continuous(expand = c(0.2, 0.2)) +
    scale_y_continuous(expand = c(0.2, 0.2)) +
    ggraph::facet_edges(stats::as.formula(paste0("~", group_col)), nrow = 2) +
    theme(
      aspect.ratio = 1,
      panel.background = element_rect(fill = "white", color = "black"),
      panel.spacing = unit(1, "lines")
    ) +
    labs(
      edge_color = "Mean\nPearson Z",
      title = glue::glue("Colocalization Score Network {cell_type}"),
      subtitle = glue::glue(
        "Metric threshold{if (only_positive) ' (only positive)' else ' (absolute value)'}: {min_abs_z}, Percentage of cells above threshold: {min_pct_detected}"
      ),
      caption = glue::glue(
        "Mean pearson Z per sample; edges shown where the pair is detected in > {min_pct_detected * 100}% of cells.
        Samples: {samples_label}{if (nzchar(excluded_label)) paste0('\nExcluded (too few cells): ', excluded_label) else ''}"
      )
    ) +
    ggraph::scale_edge_width(guide = "none")
}

#' Styled network figure: one row per donor (T0 -> T1) with a shared legend
#'
#' Uses all pairs of the summary (no min_pct_detected / min_abs_z filter), as in
#' the original figures.
#'
#' @param summary_data Output of `SummarizeColocNetworkData()` with
#'   group_col = "sample_alias" (needs donor and timepoint columns).
#' @inheritParams PlotColocNetwork
PlotColocNetworkStyled <- function(
  summary_data,
  value_limit = 5,
  layout_seed = 1,
  legend_sample = NULL,
  plot_margin = ggplot2::margin(4, 4, 4, 4),
  min_cells = NULL
) {
  summary_data <- keep_included_samples(summary_data, min_cells)
  edge_data <- summary_data |>
    dplyr::transmute(
      sample_alias,
      donor = stringr::str_remove(donor, "^donor"), # "donorD" -> "D"
      timepoint,
      from = marker_1,
      to = marker_2,
      mean_join_count_z
    )

  # Fixed node layout shared by all panels, normalised to 0-1
  g <- igraph::graph_from_data_frame(
    edge_data |> dplyr::distinct(from, to),
    directed = FALSE
  )
  set.seed(layout_seed)
  layout <- igraph::layout_with_fr(g)
  layout_norm <- apply(layout, 2, function(v) (v - min(v)) / (max(v) - min(v)))
  node_coords <- tibble::tibble(
    node = igraph::V(g)$name,
    x = layout_norm[, 1],
    y = layout_norm[, 2]
  )

  plot_network <- function(df) {
    df <- df |>
      dplyr::left_join(node_coords, by = c("from" = "node")) |>
      dplyr::left_join(
        node_coords,
        by = c("to" = "node"),
        suffix = c("", "_end")
      ) |>
      dplyr::rename(xend = x_end, yend = y_end)

    ggplot() +
      geom_segment(
        data = df,
        aes(x = x, y = y, xend = xend, yend = yend, color = mean_join_count_z),
        lineend = "round",
        linewidth = 5
      ) +
      geom_label(
        data = node_coords,
        aes(x = x, y = y, label = node),
        size = 6.2,
        fontface = "bold",
        label.padding = unit(0.35, "lines"),
        label.r = unit(0.1, "lines"),
        fill = "white"
      ) +
      scale_color_distiller(
        palette = "RdBu",
        direction = -1,
        limits = c(-value_limit, value_limit),
        oob = scales::squish
      ) +
      xlim(-0.25, 1.25) +
      ylim(-0.25, 1.25) +
      coord_fixed(clip = "off") +
      theme_void() +
      theme(legend.position = "none", plot.margin = plot_margin)
  }

  timepoint_label <- c(t0 = "Pre-infusion (T0)", t1 = "Post-infusion (T1)")

  make_donor_row <- function(donor_id) {
    d <- edge_data |>
      dplyr::filter(donor == donor_id) |>
      dplyr::mutate(
        mean_join_count_z = dplyr::case_when(
          mean_join_count_z > value_limit ~ value_limit,
          mean_join_count_z < -value_limit ~ -value_limit,
          TRUE ~ mean_join_count_z
        )
      )

    tps <- intersect(c("t0", "t1"), unique(d$timepoint))
    panels <- lapply(tps, function(tp) {
      plot_network(d |> dplyr::filter(timepoint == tp)) +
        labs(subtitle = timepoint_label[[tp]]) +
        theme(
          plot.subtitle = element_text(
            hjust = 0.5,
            size = 19,
            face = "plain",
            margin = margin(b = 6)
          )
        )
    })

    if (length(panels) == 2) {
      arrow_panel <- ggplot() +
        geom_segment(
          aes(x = 0, xend = 1, y = 0, yend = 0),
          linewidth = 2.2,
          color = "black",
          arrow = arrow(length = unit(0.28, "inches"), type = "closed")
        ) +
        xlim(0, 1) +
        ylim(-0.3, 0.3) +
        theme_void()
      row <- panels[[1]] +
        arrow_panel +
        panels[[2]] +
        patchwork::plot_layout(widths = c(1, 0.35, 1))
    } else {
      row <- panels[[1]]
    }

    row_title <- ggplot() +
      annotate(
        "text",
        x = 0,
        y = 0,
        label = paste("Donor", donor_id),
        size = 9,
        fontface = "bold"
      ) +
      theme_void()

    row_title / row + patchwork::plot_layout(heights = c(0.18, 1))
  }

  donor_ids <- sort(unique(edge_data$donor))
  donor_rows <- lapply(donor_ids, make_donor_row)
  final_plot <- patchwork::wrap_plots(donor_rows, ncol = 1) +
    patchwork::plot_layout(heights = rep(1, length(donor_ids)))

  # Legend built separately and appended at the bottom (horizontal)
  if (is.null(legend_sample)) {
    legend_sample <- edge_data$sample_alias[[1]]
  } else if (!legend_sample %in% edge_data$sample_alias) {
    stop("legend_sample '", legend_sample, "' is not an included sample.")
  }
  legend_plot <- plot_network(
    edge_data |> dplyr::filter(sample_alias == legend_sample)
  ) +
    scale_color_distiller(
      name = "Mean Pearson Z     ",
      palette = "RdBu",
      direction = -1,
      limits = c(-value_limit, value_limit),
      oob = scales::squish
    ) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 17),
      legend.text = element_text(size = 14)
    )
  legend_only <- cowplot::get_legend(legend_plot)

  patchwork::wrap_elements(final_plot) /
    patchwork::wrap_elements(legend_only) +
    patchwork::plot_layout(heights = c(1, 0.05))
}

#' Colocalization network figures for a marker set in one cell type
#'
#' @param data Seurat object (e.g. data_t0t1_b1b2).
#' @param markers Markers to include in the network.
#' @param cell_type Cell type to plot.
#' @param samples Values of `group_col` (sample_alias, e.g. "donorD_t0") to
#'   include. NULL = all samples.
#' @param group_col Metadata column defining samples. The styled figure needs
#'   "sample_alias" (donor_timepoint) and is skipped otherwise.
#' @param min_pct_detected Edges are kept (overview plot) if the pair is
#'   detected in more than this fraction of cells.
#' @param min_abs_z Edges are kept (overview plot) if abs(mean pearson_z) is
#'   above this value.
#' @param only_positive If TRUE, keep only edges with mean pearson_z > min_abs_z.
#' @param value_limit Colour scale cap (pearson_z) in the styled figure.
#' @param layout_seed Seed for the node layout of the styled figure.
#' @param plot_margin Margin around each network panel of the styled figure,
#'   a `ggplot2::margin()` (default `margin(4, 4, 4, 4)`).
#' @param scores Optional precomputed `ColocalizationScores()` result.
#' @param legend_sample sample_alias used to draw the legend of the styled
#'   figure. Defaults to the first included sample.
#' @param min_cells Minimum number of cells of `cell_type` per sample. Samples
#'   below it are left out of both figures but stay in the returned `summary`
#'   with `included = FALSE` (columns `n_cells`, `included`).
#' @param gradient_palette Colours for the overview edge scale.
#' @param out_dir,file_stem Output location; NULL out_dir = do not export.
#' @param overview_size,styled_size c(width, height) of the exported figures.
#'
#' @return A list with `summary`, `overview` and `styled` (NULL if skipped).
PlotColocNetwork <- function(
  data,
  markers,
  cell_type,
  samples = NULL,
  group_col = "sample_alias",
  min_pct_detected = 0.5,
  min_abs_z = 0,
  only_positive = FALSE,
  value_limit = 5,
  layout_seed = 1,
  scores = NULL,
  legend_sample = NULL,
  plot_margin = ggplot2::margin(4, 4, 4, 4),
  min_cells = 1,
  gradient_palette = default_network_gradient,
  out_dir = NULL,
  file_stem = NULL,
  overview_size = c(12, 6),
  styled_size = c(10, 35)
) {
  summary_data <- SummarizeColocNetworkData(
    data,
    markers = markers,
    cell_type = cell_type,
    samples = samples,
    group_col = group_col,
    scores = scores,
    min_cells = min_cells
  )

  overview <- PlotColocNetworkOverview(
    summary_data,
    cell_type = cell_type,
    group_col = group_col,
    min_pct_detected = min_pct_detected,
    min_abs_z = min_abs_z,
    only_positive = only_positive,
    gradient_palette = gradient_palette,
    min_cells = min_cells
  )

  styled <- NULL
  if (group_col == "sample_alias") {
    styled <- PlotColocNetworkStyled(
      summary_data,
      value_limit = value_limit,
      layout_seed = layout_seed,
      legend_sample = legend_sample,
      plot_margin = plot_margin,
      min_cells = min_cells
    )
  } else {
    message("group_col is not 'sample_alias': skipping the styled figure")
  }

  if (!is.null(out_dir) && !is.null(file_stem)) {
    export_plot(
      filename = file.path(out_dir, file_stem),
      plot = overview,
      width = overview_size[1],
      height = overview_size[2]
    )
    if (!is.null(styled)) {
      export_plot(
        filename = file.path(out_dir, paste0(file_stem, "_styled")),
        plot = styled,
        width = styled_size[1],
        height = styled_size[2]
      )
    }
  }

  return(list(summary = summary_data, overview = overview, styled = styled))
}
