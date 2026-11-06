# Helpers for selecting markers by abundance (used by module 05 and reusable
# in other modules)

# Mean CLR per sample x cell type x marker, for samples with enough cells
compute_sample_clr_means <- function(
  object,
  cell_types = NULL,
  group_col = "cell_annotation",
  sample_col = "sample_alias",
  layer = "data",
  min_cells_per_sample = 15L
) {
  meta <- object@meta.data |>
    tibble::rownames_to_column("cell_id") |>
    dplyr::mutate(
      sample = as.character(.data[[sample_col]]),
      group = as.character(.data[[group_col]])
    )
  if (!is.null(cell_types)) {
    meta <- dplyr::filter(meta, group %in% cell_types)
  }

  clr_mat <- SeuratObject::LayerData(object, layer = layer)
  group_id <- paste(meta$sample, meta$group, sep = "||")
  sums <- rowsum(as.matrix(Matrix::t(clr_mat[, meta$cell_id])), group_id)
  means <- sums / as.vector(table(group_id)[rownames(sums)])

  tibble::as_tibble(means, rownames = "id") |>
    tidyr::pivot_longer(-id, names_to = "marker", values_to = "mean_clr") |>
    tidyr::separate(id, into = c("sample", "group"), sep = "\\|\\|") |>
    dplyr::left_join(
      dplyr::count(meta, sample, group, name = "n_cells"),
      by = c("sample", "group")
    ) |>
    dplyr::filter(n_cells >= min_cells_per_sample)
}

# Summary per cell type and marker: lowest sample mean CLR and pass/fail
summarise_marker_clr <- function(sample_means, clr_cutoff = 0.5) {
  sample_means |>
    dplyr::summarise(
      min_sample_mean_clr = min(mean_clr),
      n_samples = dplyr::n(),
      .by = c(group, marker)
    ) |>
    dplyr::mutate(pass = min_sample_mean_clr > clr_cutoff)
}

# Filter a named list of marker vectors (one per cell type): keep markers whose
# mean CLR is > clr_cutoff in every sample with >= min_cells_per_sample cells.
# `always_keep` markers bypass the cutoff. Returns the filtered list, with the
# per-marker summary attached as attribute "summary".
select_markers_by_clr <- function(
  object,
  markers_to_test,
  clr_cutoff = 0.5,
  min_cells_per_sample = 15L,
  always_keep = character(0),
  ...
) {
  summary_tbl <- compute_sample_clr_means(
    object,
    cell_types = names(markers_to_test),
    min_cells_per_sample = min_cells_per_sample,
    ...
  ) |>
    summarise_marker_clr(clr_cutoff)

  no_samples <- setdiff(names(markers_to_test), unique(summary_tbl$group))
  if (length(no_samples) > 0) {
    warning(
      "No sample has >= ", min_cells_per_sample, " cells for: ",
      paste(no_samples, collapse = ", "),
      ". Only `always_keep` markers are retained for these cell types."
    )
  }

  selected <- lapply(names(markers_to_test), function(ct) {
    passing <- summary_tbl |>
      dplyr::filter(group == ct, pass) |>
      dplyr::pull(marker)
    sort(intersect(markers_to_test[[ct]], union(passing, always_keep)))
  }) |>
    setNames(names(markers_to_test))

  attr(selected, "summary") <- summary_tbl
  selected
}
