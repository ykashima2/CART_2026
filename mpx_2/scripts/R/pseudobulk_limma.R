#' Pseudobulk differential abundance with limma
#'
#' Sample-level (pseudobulk) differential abundance testing for MPX data,
#' ported from the pseudobulk limma path of the PNA proxiome analysis template
#' (pna_b1b2/modules/05_statistical_testing.qmd, sections 5.4-5.5) and
#' simplified to a single pairwise contrast.
#'
#' Per cell type, CLR values are averaged per sample (samples with fewer than
#' `min_cells` cells of that type are dropped), and limma is fitted across the
#' samples. This treats samples, not cells, as the replicates, so it does not
#' overstate significance the way cell-level tests (e.g. `RunDAA()`) can.
#'
#' All functions use explicit namespaces (dplyr::, tidyr::, ...) because this
#' file is sourced before the analysis script loads its packages.

#' Mean CLR per sample and cell type
#'
#' @param object Seurat object with CLR-normalised MPX data.
#' @param group_col Metadata column defining the cell types to test separately.
#' @param sample_col Metadata column defining samples (the replicates).
#' @param contrast_column Metadata column holding the conditions to compare.
#' @param covariate_cols Metadata columns to carry along as sample covariates.
#'   They must be constant within a sample.
#' @param features Markers to include. NULL (default) uses all markers.
#' @param assay,layer Assay and layer to fetch values from.
#' @param min_cells Minimum cells per sample and cell type for the sample to be
#'   included for that cell type.
#'
#' @return A list with `sample_means` (long table: group, sample, condition,
#'   covariates, n_cells, marker, mean_clr) and `excluded` (sample x cell type
#'   combinations below `min_cells`, with their cell counts).
pseudobulk_mean_clr <- function(
  object,
  group_col = "celltype",
  sample_col = "sample_alias",
  contrast_column = "timepoint",
  covariate_cols = NULL,
  features = NULL,
  assay = "mpxCells",
  layer = "data",
  min_cells = 20
) {
  meta_cols <- unique(c(group_col, sample_col, contrast_column, covariate_cols))
  missing_cols <- setdiff(meta_cols, colnames(object[[]]))
  if (length(missing_cols) > 0) {
    stop(
      "Missing metadata column(s): ",
      paste(missing_cols, collapse = ", ")
    )
  }
  if (is.null(features)) {
    features <- rownames(object[[assay]])
  }
  features <- setdiff(as.character(features), meta_cols)

  meta <- object[[]][, meta_cols, drop = FALSE] |>
    tibble::as_tibble(rownames = "cell_id") |>
    dplyr::mutate(dplyr::across(dplyr::all_of(meta_cols), as.character))

  # Covariates must be sample-level: one value per sample
  if (length(covariate_cols) > 0) {
    n_levels <- meta |>
      dplyr::summarise(
        dplyr::across(dplyr::all_of(covariate_cols), dplyr::n_distinct),
        .by = dplyr::all_of(sample_col)
      )
    varying <- covariate_cols[vapply(
      covariate_cols,
      \(col) any(n_levels[[col]] > 1),
      logical(1)
    )]
    if (length(varying) > 0) {
      stop(
        "Covariate(s) vary within samples: ",
        paste(varying, collapse = ", ")
      )
    }
  }

  cell_counts <- meta |>
    dplyr::count(
      dplyr::across(dplyr::all_of(c(
        group_col,
        meta_cols[meta_cols != group_col]
      ))),
      name = "n_cells"
    )
  excluded <- cell_counts |> dplyr::filter(.data$n_cells < min_cells)
  kept <- cell_counts |> dplyr::filter(.data$n_cells >= min_cells)

  values <- SeuratObject::FetchData(
    object,
    vars = features,
    assay = assay,
    layer = layer
  ) |>
    tibble::as_tibble(rownames = "cell_id")

  sample_means <- meta |>
    dplyr::semi_join(kept, by = c(group_col, sample_col)) |>
    dplyr::inner_join(values, by = "cell_id") |>
    tidyr::pivot_longer(
      dplyr::all_of(features),
      names_to = "marker",
      values_to = "clr"
    ) |>
    dplyr::summarise(
      mean_clr = mean(.data$clr, na.rm = TRUE),
      .by = dplyr::all_of(c(meta_cols, "marker"))
    ) |>
    dplyr::left_join(kept, by = meta_cols)

  list(sample_means = sample_means, excluded = excluded)
}


#' Fit limma for one pairwise contrast on a features x samples matrix
#'
#' Port of `pseudobulk_limma_fit()` from the PNA module 05, for a single
#' target-vs-reference contrast.
#'
#' @param expr_mat Features x samples matrix of sample-level means.
#' @param sample_info Data frame with one row per column of `expr_mat`, in the
#'   same order, holding `sample_col`, `contrast_column` and any covariates.
#' @param covariate_cols Fixed-effect covariates added to the design.
#' @param block_col Optional blocking column (e.g. donor) for repeated measures,
#'   handled with `limma::duplicateCorrelation` instead of as fixed effects.
#' @param label Label used in messages (e.g. the cell type).
#'
#' @return A tibble with feature, difference (logFC), p, target, reference.
pseudobulk_limma_fit <- function(
  expr_mat,
  sample_info,
  sample_col = "sample_alias",
  contrast_column = "timepoint",
  reference = "t0",
  target = "t1",
  covariate_cols = NULL,
  block_col = NULL,
  label = NA_character_
) {
  sample_info <- as.data.frame(sample_info, stringsAsFactors = FALSE)
  stopifnot(identical(
    colnames(expr_mat),
    as.character(sample_info[[sample_col]])
  ))

  condition <- factor(
    as.character(sample_info[[contrast_column]]),
    levels = c(reference, target)
  )
  if (anyNA(condition)) {
    stop(
      "sample_info contains conditions other than ",
      reference,
      " and ",
      target
    )
  }

  # One coefficient per condition (no intercept) so the contrast is explicit
  design <- stats::model.matrix(~ 0 + condition)
  colnames(design) <- levels(condition)
  if (length(covariate_cols) > 0) {
    covar_formula <- stats::as.formula(paste0(
      "~ ",
      paste(covariate_cols, collapse = " + ")
    ))
    # Drop the covariate intercept so the design is not over-parameterised
    covar_design <- stats::model.matrix(covar_formula, data = sample_info)[,
      -1,
      drop = FALSE
    ]
    design <- cbind(design, covar_design)
  }

  design_qr <- qr(design)
  if (design_qr$rank < ncol(design)) {
    aliased <- colnames(design)[
      design_qr$pivot[(design_qr$rank + 1L):ncol(design)]
    ]
    warning(
      "Pseudobulk design for ",
      label,
      " is rank-deficient (rank ",
      design_qr$rank,
      " of ",
      ncol(design),
      " columns). Aliased: ",
      paste(aliased, collapse = ", "),
      ". The contrast may be returned as NA.",
      call. = FALSE
    )
  }

  # Numeric contrast (target - reference). makeContrasts() parses names as R
  # expressions and fails on hyphens or leading digits.
  contrast_matrix <- matrix(
    0,
    nrow = ncol(design),
    ncol = 1,
    dimnames = list(colnames(design), paste0(target, "-", reference))
  )
  contrast_matrix[target, 1] <- 1
  contrast_matrix[reference, 1] <- -1

  # Optional repeated-measures blocking (e.g. donor), as in the PNA module 05
  block <- NULL
  if (!is.null(block_col) && nzchar(block_col)) {
    candidate <- as.character(sample_info[[block_col]])
    if (anyNA(candidate)) {
      warning(
        "Blocking column ",
        block_col,
        " has missing values; fitting without blocking.",
        call. = FALSE
      )
    } else if (dplyr::n_distinct(candidate) < nrow(sample_info)) {
      block <- candidate
    } else {
      message(
        label,
        ": no repeated ",
        block_col,
        " values; fitting without blocking."
      )
    }
  }

  if (is.null(block)) {
    fit <- limma::lmFit(expr_mat, design)
  } else {
    dupcor <- try(
      limma::duplicateCorrelation(expr_mat, design, block = block),
      silent = TRUE
    )
    consensus <- if (inherits(dupcor, "try-error")) {
      NA_real_
    } else {
      dupcor$consensus.correlation
    }
    if (!is.finite(consensus)) {
      warning(
        "duplicateCorrelation failed for ",
        label,
        "; fitting without blocking.",
        call. = FALSE
      )
      fit <- limma::lmFit(expr_mat, design)
    } else {
      message(label, ": duplicateCorrelation consensus ", signif(consensus, 3))
      fit <- limma::lmFit(
        expr_mat,
        design,
        block = block,
        correlation = consensus
      )
    }
  }

  # robust = TRUE down-weights features with outlying variances
  fit2 <- limma::eBayes(
    limma::contrasts.fit(fit, contrast_matrix),
    robust = TRUE
  )

  limma::topTable(fit2, coef = 1, number = Inf, sort.by = "none") |>
    tibble::as_tibble(rownames = "feature") |>
    dplyr::transmute(
      feature = .data$feature,
      difference = .data$logFC,
      p = .data$P.Value,
      target = target,
      reference = reference
    )
}


#' Run pseudobulk limma for each cell type (internal)
#'
#' Shared by RunPseudobulkDAA() and RunPseudobulkDCA(). For each cell type:
#' checks the samples per group, drops single-level covariates, builds the
#' features x samples matrix and fits limma. FDR is applied across all
#' features and cell types.
#'
#' @param sample_means Long table with one row per cell type, sample and
#'   feature: `group_col`, `sample_col`, `contrast_column`, covariates,
#'   `n_cells`, `feature_col` and `value_col`.
#' @param values_fill Value for features missing in a sample (NULL keeps NA).
#'   Colocalization uses 0, as in the PNA analysis.
#'
#' @return A list with `results` (one row per cell type and feature) and
#'   `skipped`.
run_pseudobulk_limma_groups <- function(
  sample_means,
  feature_col,
  value_col,
  contrast_column,
  reference,
  target,
  group_col,
  sample_col,
  covariate_cols = NULL,
  block_col = NULL,
  values_fill = NULL,
  min_samples_per_group = 2
) {
  skipped <- list()
  results <- list()

  for (cell_type in sort(unique(sample_means[[group_col]]))) {
    ct_means <- sample_means |> dplyr::filter(.data[[group_col]] == cell_type)
    sample_info <- ct_means |>
      dplyr::distinct(dplyr::across(
        dplyr::all_of(c(
          sample_col,
          contrast_column,
          covariate_cols,
          block_col,
          "n_cells"
        ))
      )) |>
      dplyr::arrange(.data[[sample_col]])

    n_ref <- sum(sample_info[[contrast_column]] == reference)
    n_tgt <- sum(sample_info[[contrast_column]] == target)
    if (n_ref < min_samples_per_group || n_tgt < min_samples_per_group) {
      skipped[[cell_type]] <- tibble::tibble(
        celltype = cell_type,
        n_samples_reference = n_ref,
        n_samples_target = n_tgt,
        reason = paste0(
          "fewer than ",
          min_samples_per_group,
          " samples in a group"
        )
      )
      next
    }

    # Drop covariates with a single level in this cell type's samples: they
    # carry no information and would make the design rank-deficient
    covs_used <- covariate_cols[vapply(
      covariate_cols,
      \(col) dplyr::n_distinct(sample_info[[col]]) > 1,
      logical(1)
    )]
    covs_dropped <- setdiff(covariate_cols, covs_used)
    if (length(covs_dropped) > 0) {
      message(
        cell_type,
        ": covariate(s) with a single level dropped: ",
        paste(covs_dropped, collapse = ", ")
      )
    }

    n_terms <- 2 +
      sum(vapply(
        covs_used,
        \(col) dplyr::n_distinct(sample_info[[col]]) - 1,
        numeric(1)
      ))
    if (nrow(sample_info) <= n_terms) {
      skipped[[cell_type]] <- tibble::tibble(
        celltype = cell_type,
        n_samples_reference = n_ref,
        n_samples_target = n_tgt,
        reason = paste0("not more samples than model terms (", n_terms, ")")
      )
      next
    }

    expr_mat <- ct_means |>
      dplyr::select(dplyr::all_of(c(feature_col, sample_col, value_col))) |>
      tidyr::pivot_wider(
        names_from = dplyr::all_of(sample_col),
        values_from = dplyr::all_of(value_col),
        values_fill = values_fill
      ) |>
      tibble::column_to_rownames(feature_col) |>
      as.matrix()
    expr_mat <- expr_mat[,
      as.character(sample_info[[sample_col]]),
      drop = FALSE
    ]

    res <- pseudobulk_limma_fit(
      expr_mat,
      sample_info,
      sample_col = sample_col,
      contrast_column = contrast_column,
      reference = reference,
      target = target,
      covariate_cols = covs_used,
      block_col = block_col,
      label = cell_type
    )
    results[[cell_type]] <- res |>
      dplyr::rename(!!feature_col := "feature") |>
      dplyr::mutate(
        celltype = cell_type,
        n_samples_target = n_tgt,
        n_samples_reference = n_ref,
        covariates = if (length(covs_used) > 0) {
          paste(covs_used, collapse = ", ")
        } else {
          "none"
        },
        method = "limma_pseudobulk"
      )
  }

  results <- dplyr::bind_rows(results)
  if (nrow(results) > 0) {
    # One FDR family per contrast: all features across all cell types
    results <- results |>
      dplyr::mutate(p_adj = stats::p.adjust(.data$p, method = "fdr")) |>
      dplyr::relocate(
        dplyr::all_of(feature_col),
        "celltype",
        "difference",
        "p",
        "p_adj",
        "target",
        "reference"
      )
  }

  list(results = results, skipped = dplyr::bind_rows(skipped))
}


#' Pseudobulk differential abundance per cell type
#'
#' For each cell type: mean CLR per sample (samples with at least `min_cells`
#' cells), then a limma test of `target` vs `reference`. Unpaired by default;
#' pass `block_col` (e.g. "donor") for a paired / repeated-measures design.
#'
#' @param object Seurat object with CLR-normalised MPX data.
#' @param contrast_column,reference,target Condition column and the two levels
#'   to compare. `difference` is target - reference.
#' @param group_col Metadata column with the cell types to test separately.
#' @param sample_col Metadata column defining samples (the replicates).
#' @param covariate_cols Fixed-effect covariates (e.g. "batch"). A covariate
#'   with a single level among a cell type's samples is dropped for that cell
#'   type, with a message.
#' @param block_col Optional blocking column for repeated measures.
#' @param features Markers to test. NULL (default) tests all markers.
#' @param min_cells Minimum cells per sample and cell type (default 30).
#' @param min_samples_per_group Minimum samples in each of target and reference
#'   for a cell type to be tested (default 2).
#' @param assay,layer Assay and layer holding the CLR values.
#'
#' @return A list with
#'   - `results`: marker, celltype (the `group_col` value), difference, p,
#'     p_adj, target, reference, n_samples_target, n_samples_reference,
#'     covariates, method. p_adj is FDR across all markers and cell types.
#'   - `sample_means`: the sample-level means used as limma input.
#'   - `excluded_samples`: sample x cell type combinations below `min_cells`.
#'   - `skipped`: cell types that were not tested, and why.
RunPseudobulkDAA <- function(
  object,
  contrast_column = "timepoint",
  reference = "t0",
  target = "t1",
  group_col = "celltype",
  sample_col = "sample_alias",
  covariate_cols = NULL,
  block_col = NULL,
  features = NULL,
  min_cells = 30,
  min_samples_per_group = 2,
  assay = "mpxCells",
  layer = "data"
) {
  pb <- pseudobulk_mean_clr(
    object,
    group_col = group_col,
    sample_col = sample_col,
    contrast_column = contrast_column,
    covariate_cols = unique(c(covariate_cols, block_col)),
    features = features,
    assay = assay,
    layer = layer,
    min_cells = min_cells
  )
  sample_means <- pb$sample_means |>
    dplyr::filter(.data[[contrast_column]] %in% c(reference, target))

  fit <- run_pseudobulk_limma_groups(
    sample_means,
    feature_col = "marker",
    value_col = "mean_clr",
    contrast_column = contrast_column,
    reference = reference,
    target = target,
    group_col = group_col,
    sample_col = sample_col,
    covariate_cols = covariate_cols,
    block_col = block_col,
    min_samples_per_group = min_samples_per_group
  )

  list(
    results = fit$results,
    sample_means = sample_means,
    excluded_samples = pb$excluded,
    skipped = fit$skipped
  )
}


#' Mean colocalization score (pearson_z) per sample and marker pair
#'
#' For one cell type, computes MPX colocalization scores and averages pearson_z
#' per sample and marker pair. Cells in which a pair is not observed count as 0
#' (the `include_missing_obs = TRUE` definition of SummarizeProximityScores()),
#' with the denominator taken from the actual number of cells per sample, so
#' filtering to selected markers cannot change it.
#'
#' @param object Seurat object subset to one cell type.
#' @param markers Markers to keep; pairs with both markers in this set are used.
#' @param samples Samples to keep (values of `sample_col`).
#' @param sample_col,meta_cols Sample column and further metadata columns
#'   (condition, covariates) carried into the output.
#'
#' @return Long table: sample and metadata columns, n_cells, marker_1,
#'   marker_2, pair ("marker_1__marker_2", markers in alphabetical order),
#'   mean_pearson_z, pct_detected.
pseudobulk_mean_coloc <- function(
  object,
  markers,
  samples,
  sample_col = "sample_alias",
  meta_cols = NULL
) {
  meta_cols <- unique(c(sample_col, meta_cols))
  object <- subset(
    object,
    cells = colnames(object)[object[[]][[sample_col]] %in% samples]
  )

  # Actual cells per sample: the denominator for the mean
  n_cells <- object[[]] |>
    tibble::as_tibble() |>
    dplyr::mutate(dplyr::across(dplyr::all_of(meta_cols), as.character)) |>
    dplyr::count(dplyr::across(dplyr::all_of(meta_cols)), name = "n_cells")

  # ColocalizationScores() errors if the object has a "component" meta column
  object$component <- NULL

  pixelatorR::ColocalizationScores(object, meta_data_columns = meta_cols) |>
    dplyr::filter(
      .data$marker_1 %in% markers,
      .data$marker_2 %in% markers,
      .data$marker_1 != .data$marker_2
    ) |>
    dplyr::mutate(
      dplyr::across(dplyr::all_of(meta_cols), as.character),
      # One fixed order per pair, so a pair cannot appear twice
      m1 = pmin(.data$marker_1, .data$marker_2),
      m2 = pmax(.data$marker_1, .data$marker_2)
    ) |>
    dplyr::summarise(
      sum_pearson_z = sum(.data$pearson_z, na.rm = TRUE),
      n_detected = dplyr::n(),
      .by = dplyr::all_of(c(meta_cols, "m1", "m2"))
    ) |>
    dplyr::rename(marker_1 = "m1", marker_2 = "m2") |>
    dplyr::left_join(n_cells, by = meta_cols) |>
    dplyr::mutate(
      pair = paste(.data$marker_1, .data$marker_2, sep = "__"),
      mean_pearson_z = .data$sum_pearson_z / .data$n_cells,
      pct_detected = .data$n_detected / .data$n_cells
    ) |>
    dplyr::select(-"sum_pearson_z", -"n_detected")
}


#' Pseudobulk differential colocalization per cell type
#'
#' For each cell type: samples with at least `min_cells` cells; markers with a
#' mean CLR at or above `min_mean_clr` in every one of those samples; mean pearson_z
#' per sample for all pairs of those markers; then a limma test of `target` vs
#' `reference`, as in RunPseudobulkDAA(). A pair missing in a sample enters the
#' test as 0, as in the PNA analysis.
#'
#' @inheritParams RunPseudobulkDAA
#' @param min_mean_clr Abundance filter: a marker is included for a cell type
#'   only if its mean CLR is at or above this value in all tested samples (default 1).
#'
#' @return A list with
#'   - `results`: pair, marker_1, marker_2, celltype, difference (mean
#'     pearson_z, target - reference), p, p_adj (FDR across all pairs and cell
#'     types), target, reference, n_samples_target, n_samples_reference,
#'     covariates, method.
#'   - `sample_means`: per-sample mean pearson_z used as limma input.
#'   - `excluded_samples`: sample x cell type combinations below `min_cells`.
#'   - `markers_tested`: per cell type, the markers passing the CLR filter.
#'   - `skipped`: cell types that were not tested, and why.
RunPseudobulkDCA <- function(
  object,
  contrast_column = "timepoint",
  reference = "t0",
  target = "t1",
  group_col = "celltype",
  sample_col = "sample_alias",
  covariate_cols = NULL,
  block_col = NULL,
  features = NULL,
  min_cells = 30,
  min_mean_clr = 1,
  min_samples_per_group = 2,
  assay = "mpxCells",
  layer = "data"
) {
  meta_cols <- unique(c(contrast_column, covariate_cols, block_col))

  # Samples passing the cell filter, and their mean CLR per marker
  pb <- pseudobulk_mean_clr(
    object,
    group_col = group_col,
    sample_col = sample_col,
    contrast_column = contrast_column,
    covariate_cols = unique(c(covariate_cols, block_col)),
    features = features,
    assay = assay,
    layer = layer,
    min_cells = min_cells
  )
  clr_means <- pb$sample_means |>
    dplyr::filter(.data[[contrast_column]] %in% c(reference, target))

  # Abundance filter: mean CLR above the threshold in every tested sample
  markers_tested <- clr_means |>
    dplyr::summarise(
      n_samples = dplyr::n_distinct(.data[[sample_col]]),
      min_sample_mean_clr = min(.data$mean_clr),
      .by = dplyr::all_of(c(group_col, "marker"))
    ) |>
    dplyr::filter(.data$min_sample_mean_clr >= min_mean_clr)

  coloc_means <- list()
  early_skipped <- list()
  for (cell_type in sort(unique(clr_means[[group_col]]))) {
    ct_markers <- markers_tested |>
      dplyr::filter(.data[[group_col]] == cell_type) |>
      dplyr::pull("marker")
    ct_samples <- clr_means |>
      dplyr::filter(.data[[group_col]] == cell_type) |>
      dplyr::pull(dplyr::all_of(sample_col)) |>
      unique()
    ct_n <- clr_means |>
      dplyr::filter(.data[[group_col]] == cell_type) |>
      dplyr::distinct(dplyr::across(dplyr::all_of(c(
        sample_col,
        contrast_column
      ))))
    n_ref <- sum(ct_n[[contrast_column]] == reference)
    n_tgt <- sum(ct_n[[contrast_column]] == target)

    # Skip before computing colocalization scores when the cell type cannot be
    # tested anyway (same rules as run_pseudobulk_limma_groups())
    skip_reason <- if (
      n_ref < min_samples_per_group || n_tgt < min_samples_per_group
    ) {
      paste0("fewer than ", min_samples_per_group, " samples in a group")
    } else if (length(ct_markers) < 2) {
      paste0(
        "fewer than 2 markers with mean CLR > ",
        min_mean_clr,
        " in all samples"
      )
    } else {
      NULL
    }
    if (!is.null(skip_reason)) {
      message(cell_type, ": not tested (", skip_reason, ")")
      early_skipped[[cell_type]] <- tibble::tibble(
        celltype = cell_type,
        n_samples_reference = n_ref,
        n_samples_target = n_tgt,
        reason = skip_reason
      )
      next
    }
    message(
      cell_type,
      ": ",
      length(ct_markers),
      " markers, ",
      choose(length(ct_markers), 2),
      " pairs, ",
      length(ct_samples),
      " samples"
    )
    ct_object <- subset(
      object,
      cells = colnames(object)[object[[]][[group_col]] == cell_type]
    )
    coloc_means[[cell_type]] <- pseudobulk_mean_coloc(
      ct_object,
      markers = ct_markers,
      samples = ct_samples,
      sample_col = sample_col,
      meta_cols = meta_cols
    ) |>
      dplyr::mutate(!!group_col := cell_type)
    rm(ct_object)
    gc()
  }
  sample_means <- dplyr::bind_rows(coloc_means)

  fit <- if (nrow(sample_means) > 0) {
    run_pseudobulk_limma_groups(
      sample_means,
      feature_col = "pair",
      value_col = "mean_pearson_z",
      contrast_column = contrast_column,
      reference = reference,
      target = target,
      group_col = group_col,
      sample_col = sample_col,
      covariate_cols = covariate_cols,
      block_col = block_col,
      values_fill = 0,
      min_samples_per_group = min_samples_per_group
    )
  } else {
    list(results = tibble::tibble(), skipped = tibble::tibble())
  }

  results <- fit$results
  if (nrow(results) > 0) {
    results <- results |>
      dplyr::left_join(
        sample_means |>
          dplyr::distinct(.data$pair, .data$marker_1, .data$marker_2),
        by = "pair"
      ) |>
      dplyr::relocate("pair", "marker_1", "marker_2")
  }

  list(
    results = results,
    sample_means = sample_means,
    excluded_samples = pb$excluded,
    markers_tested = markers_tested,
    skipped = dplyr::bind_rows(dplyr::bind_rows(early_skipped), fit$skipped)
  )
}
