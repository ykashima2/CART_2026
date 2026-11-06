#' Filter markers based on isotype control abundance or sample prevalence
#'
#' @description
#' Filters markers by comparing their median (or mean) abundance to isotype control abundance (default behavior).
#' Alternatively, filters markers based on their prevalence across samples using a fixed threshold.
#' Removes isotype controls and markers matching `exclude_regex`.
#'
#' @param data A Seurat object containing marker expression data
#' @param group_vars Character string. Metadata column to group by (e.g. "seurat_clusters"). If provided, the function returns a list of results for each group. If NULL (default), operates on the full data.
#' @param ratio Numeric threshold applied to the isotype statistic (default: 1.5). Used if `threshold` is NULL.
#' @param threshold Numeric threshold value. If provided, overrides the isotype-based threshold calculation.
#' @param min_percent Minimum fraction of cells whose expression must exceed the threshold for a marker to be retained (default: 0.1). Used in global filtering mode.
#' @param min_samples Integer. Minimum number of samples where the marker's statistic exceeds the threshold. If provided, enables sample-based filtering mode.
#' @param sample_col Character string. Column name in metadata defining sample groups. Required if `min_samples` is used.
#' @param statistic Character string. Statistic to calculate: "median" (default) or "mean". Affects global statistics, per-sample statistics, and isotype reference calculations.
#' @param layer Character string specifying which assay slot/layer to pull from, e.g. "counts" or "data" (default: "counts")
#' @param isotype_regex Regular expression identifying isotype control rows (default: "mIgG")
#' @param exclude_regex Regular expression of markers to discard from the result (default: "mIgG|ACTB")
#' @param return_stats Logical; if TRUE returns a tibble with per-marker statistics instead of a character vector (default: FALSE)
#'
#' @return If `group_vars` is NULL: Character vector of marker names (or tibble if `return_stats=TRUE`).
#'         If `group_vars` is provided: Named list where each element is the result for that group.
#'
#' @details
#' The function operates in two modes:
#'
#' **1. Global Filtering Mode (Default)**
#' - Calculates global statistic (mean or median) for each marker.
#' - If `threshold` is NULL, calculates isotype statistic (mean or median) and sets `threshold = ratio * isotype_stat`.
#' - Retains markers where `statistic > threshold` AND `percent_above >= min_percent`.
#'
#' **2. Sample-Based Filtering Mode (if `min_samples` is provided)**
#' - Requires `sample_col` to be present in metadata.
#' - If `threshold` is NULL, calculates it from isotypes as above.
#' - For each marker, calculates the `statistic` (mean or median) within each sample.
#' - Counts how many samples have `statistic > threshold`.
#' - Retains markers where `n_samples_pass >= min_samples`.
#'
#' In both modes, markers matching `exclude_regex` are removed.
#'
#' @examples
#' \dontrun{
#' # Filter global dataset
#' markers <- FilterMarkers(seurat_obj, ratio = 2)
#'
#' # Filter per cluster
#' markers_list <- FilterMarkers(
#'   seurat_obj,
#'   group_vars = "seurat_clusters",
#'   ratio = 2
#' )
#'
#' # Use the normalised "data" layer and get full statistics
#' stats <- FilterMarkers(
#'   seurat_obj,
#'   ratio = 2,
#'   layer = "data",
#'   return_stats = TRUE
#' )
#' }
FilterMarkers <- function(
  data,
  group_vars = NULL,
  ratio = 1.5,
  threshold = NULL,
  min_percent = 0.1,
  min_samples = NULL,
  sample_col = NULL,
  statistic = "median",
  layer = "counts",
  isotype_regex = "mIgG",
  exclude_regex = "mIgG|ACTB",
  return_stats = FALSE
) {
  # ---------------------------------------------------------------------------
  # 0. Validate arguments -----------------------------------------------------
  # ---------------------------------------------------------------------------
  statistic <- match.arg(statistic, c("median", "mean"))

  if (!is.null(min_samples)) {
    if (is.null(sample_col)) {
      cli::cli_abort(
        "If 'min_samples' is specified, 'sample_col' must be provided."
      )
    }
    if (!sample_col %in% colnames(data@meta.data)) {
      cli::cli_abort(sprintf(
        "Sample column '%s' not found in metadata.",
        sample_col
      ))
    }
  }

  # ---------------------------------------------------------------------------
  # 1. Handle Grouping --------------------------------------------------------
  # ---------------------------------------------------------------------------
  if (!is.null(group_vars)) {
    # Check if group_vars is a valid column
    if (!group_vars %in% colnames(data@meta.data)) {
      cli::cli_abort(sprintf(
        "Grouping variable '%s' not found in metadata.",
        group_vars
      ))
    }

    # Iterate over groups
    unique_groups <- unique(data@meta.data[[group_vars]])
    # Sort for consistent order if possible
    if (!is.null(levels(unique_groups))) {
      unique_groups <- levels(unique_groups)[
        levels(unique_groups) %in% unique_groups
      ]
    } else {
      unique_groups <- sort(unique_groups)
    }

    results <- list()

    cli::cli_progress_bar(
      name = sprintf("Filtering markers by %s", group_vars),
      format = "{cli::pb_name} {cli::pb_bar} {cli::pb_percent} | {cli::pb_eta_str}",
      total = length(unique_groups)
    )

    for (grp in unique_groups) {
      cli::cli_progress_update()
      # Subset data for the group
      cells_in_group <- rownames(data@meta.data[
        data@meta.data[[group_vars]] == grp,
        ,
        drop = FALSE
      ])
      data_subset <- data[, cells_in_group]

      # Run core logic
      res <- tryCatch(
        {
          .compute_marker_stats(
            data_subset,
            ratio = ratio,
            threshold = threshold,
            min_percent = min_percent,
            min_samples = min_samples,
            sample_col = sample_col,
            statistic = statistic,
            layer = layer,
            isotype_regex = isotype_regex,
            exclude_regex = exclude_regex,
            return_stats = return_stats
          )
        },
        error = function(e) {
          cli::cli_warn(sprintf("Error in group '%s': %s", grp, e$message))
          NULL
        }
      )

      results[[as.character(grp)]] <- res
    }

    return(results)
  } else {
    # -------------------------------------------------------------------------
    # 2. Run on Full Data -----------------------------------------------------
    # -------------------------------------------------------------------------
    return(.compute_marker_stats(
      data,
      ratio = ratio,
      threshold = threshold,
      min_percent = min_percent,
      min_samples = min_samples,
      sample_col = sample_col,
      statistic = statistic,
      layer = layer,
      isotype_regex = isotype_regex,
      exclude_regex = exclude_regex,
      return_stats = return_stats
    ))
  }
}

#' Core logic for marker filtering
#' @keywords internal
.compute_marker_stats <- function(
  data,
  ratio,
  threshold,
  min_percent,
  min_samples,
  sample_col,
  statistic,
  layer,
  isotype_regex,
  exclude_regex,
  return_stats
) {
  # -----------------------------------------------------------------------
  # 3. Retrieve count matrix ---------------------------------------------
  # -----------------------------------------------------------------------
  counts_matrix <- LayerData(data, layer = layer)
  counts_mat <- as.matrix(counts_matrix)

  # -----------------------------------------------------------------------
  # 4. Identify isotype rows and compute/set threshold --------------------
  # -----------------------------------------------------------------------

  isotype_stat <- NA
  isotype_idx <- grepl(isotype_regex, rownames(counts_mat))

  if (any(isotype_idx)) {
    if (statistic == "median") {
      isotype_stat <- stats::median(counts_mat[isotype_idx, ])
    } else {
      isotype_stat <- mean(counts_mat[isotype_idx, ])
    }
  }

  # Set threshold if not manually provided
  if (is.null(threshold)) {
    if (!any(isotype_idx)) {
      cli::cli_abort(c(
        "x" = sprintf(
          "No isotype control rows found using regex '%s'.",
          isotype_regex
        ),
        "i" = "Cannot calculate automatic threshold. Please provide a manual 'threshold' or check your regex."
      ))
    }

    threshold <- ratio * isotype_stat
  }

  # -----------------------------------------------------------------------
  # 5. Compute Statistics -------------------------------------------------
  # -----------------------------------------------------------------------

  # Global statistics
  if (statistic == "median") {
    global_stats <- matrixStats::rowMedians(counts_mat)
  } else {
    global_stats <- rowMeans(counts_mat)
  }

  percent_above <- rowMeans(counts_mat > threshold)

  # Sample-based statistics (if used)
  n_samples_pass <- NULL

  if (!is.null(min_samples)) {
    sample_ids <- data@meta.data[[sample_col]]
    unique_samples <- unique(sample_ids)

    n_samples_pass <- rep(0, nrow(counts_mat))

    for (s in unique_samples) {
      cells_in_sample <- which(sample_ids == s)

      if (length(cells_in_sample) == 0) {
        next
      }

      mat_s <- counts_mat[, cells_in_sample, drop = FALSE]

      if (statistic == "median") {
        stats_s <- matrixStats::rowMedians(mat_s)
      } else {
        stats_s <- rowMeans(mat_s)
      }

      # Increment count if statistic exceeds threshold
      n_samples_pass <- n_samples_pass + (stats_s > threshold)
    }
  }

  # -----------------------------------------------------------------------
  # 6. Assemble and Filter Results ----------------------------------------
  # -----------------------------------------------------------------------

  markers_stats <- tibble::tibble(
    marker = rownames(counts_matrix)
  )

  # Add dynamic statistic columns
  markers_stats[[statistic]] <- global_stats
  markers_stats[[paste0("isotype_", statistic)]] <- isotype_stat

  # Add remaining columns
  markers_stats$percent_above <- percent_above
  markers_stats$threshold_used <- threshold

  # Add sample stats if available
  if (!is.null(n_samples_pass)) {
    markers_stats$n_samples_pass <- n_samples_pass
  }

  # -----------------------------------------------------------------------
  # 7. Apply Filters -----------------------------------------------------
  # -----------------------------------------------------------------------
  if (!is.null(min_samples)) {
    # Sample-based filtering mode
    markers_stats <- markers_stats %>%
      dplyr::filter(
        n_samples_pass >= min_samples,
        !grepl(exclude_regex, marker)
      )
  } else {
    # Global filtering mode (default)
    markers_stats <- markers_stats %>%
      dplyr::filter(
        .data[[statistic]] > threshold,
        percent_above >= min_percent,
        !grepl(exclude_regex, marker)
      )
  }

  if (return_stats) {
    return(markers_stats)
  } else {
    return(markers_stats$marker)
  }
}
