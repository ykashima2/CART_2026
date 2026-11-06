# Helpers for the batched platelet clean-up (module 02b_platelet_cleanup.qmd)
#
# Strategy (tested on 100 cells in module 02):
#   1. segment every cell graph into T cell / platelet / interface / other nodes
#      with protein weights `w` (pixelatorR::segment_cell)
#   2. compute the platelet fraction of every cell and drop cells above a cutoff
#   3. keep only the T cell compartment of the remaining cells
#   4. recompute counts (partition_counts) and proximity scores (analytical)
#
# The work is done in batches (cells of one sample at a time) and every batch is
# written to disk, so that memory use stays bounded and a stopped run can be
# resumed. `keep` and `contam` are the compartment names given by the columns of
# `w` ("tcell" and "platelet" in this project).

batch_file_suffixes <- c(
  counts = "counts.rds",
  proximity = "proximity.rds",
  qc = "qc.rds",
  cellgraphs = "cellgraphs.rds"
)

# Write an object under a temporary name and rename it when complete, so that a
# file that exists is always a complete file
save_atomic <- function(object, path) {
  tmp <- paste0(path, ".tmp")
  saveRDS(object, tmp)
  file.rename(tmp, path)
  invisible(path)
}

batch_paths <- function(out_dir, batch_id) {
  stats::setNames(
    file.path(out_dir, "batches", paste0(batch_id, "_", batch_file_suffixes)),
    names(batch_file_suffixes)
  )
}

batch_is_done <- function(out_dir, batch_id) {
  all(file.exists(batch_paths(out_dir, batch_id)))
}

#' Split cells into batches of at most `batch_size` cells from a single sample
#'
#' @param meta Cell metadata (`pg_data[[]]`) with cell names as rownames.
#' @param cells Cells to process.
#' @param batch_size Maximum number of cells per batch.
#' @param sample_col Metadata column with the sample (one PXL file per sample).
#' @return A tibble with batch_id, sample, n_cells and a list-column `cells`.
make_batch_manifest <- function(
  meta,
  cells,
  batch_size,
  sample_col = "sample_alias"
) {
  tibble::tibble(
    cell = cells,
    sample = as.character(meta[cells, sample_col])
  ) |>
    dplyr::arrange(sample, cell) |>
    dplyr::mutate(
      batch_in_sample = ceiling(dplyr::row_number() / batch_size),
      .by = sample
    ) |>
    dplyr::summarise(
      cells = list(cell),
      n_cells = dplyr::n(),
      .by = c(sample, batch_in_sample)
    ) |>
    dplyr::mutate(batch_id = sprintf("batch_%04d", dplyr::row_number())) |>
    dplyr::select(batch_id, sample, batch_in_sample, n_cells, cells)
}

#' Resident memory of the current R process in bytes (macOS and Linux)
process_rss_bytes <- function() {
  rss_kb <- suppressWarnings(as.numeric(system2(
    "ps",
    c("-o", "rss=", "-p", Sys.getpid()),
    stdout = TRUE
  )))
  rss_kb * 1024
}

#' Load the cell graphs of pilot cells and measure the memory they use
#'
#' The graphs keep part of their data in C memory (igraph), which
#' `object.size()` does not see, so the resident memory of the R process is
#' measured before and after loading.
#'
#' @param pg_data Seurat object.
#' @param cells Pilot cells.
#' @return Bytes of memory per cell.
measure_graph_memory <- function(pg_data, cells) {
  invisible(gc(verbose = FALSE))
  before <- process_rss_bytes()
  se <- subset(pg_data, cells = cells)
  se <- pixelatorR::LoadCellGraphs(se, verbose = FALSE)
  cg_list <- pixelatorR::CellGraphs(se)
  after <- process_rss_bytes()
  rm(se, cg_list)
  invisible(gc(verbose = FALSE))
  max(after - before, 0) / length(cells)
}

#' Batch size from the memory use per cell
#'
#' @param bytes_per_cell Memory of the loaded cell graphs per cell, from
#'   `measure_graph_memory()`.
#' @param memory_budget_gb Memory that one batch may use.
#' @param safety Peak memory relative to the loaded graphs (the segmented and
#'   filtered copies exist at the same time).
calibrate_batch_size <- function(
  bytes_per_cell,
  memory_budget_gb = 8,
  safety = 3,
  min_size = 50,
  max_size = 2000
) {
  size <- floor(memory_budget_gb * 1e9 / (bytes_per_cell * safety))
  list(
    bytes_per_cell = bytes_per_cell,
    batch_size = max(min_size, min(max_size, size))
  )
}

#' Compartment statistics of one segmented cell graph (one row per cell)
#'
#' @param cg A segmented CellGraph (node attribute `compartment`).
#' @param component Cell name.
#' @param n_umi Number of UMIs of the cell before clean-up.
cell_compartment_qc <- function(
  cg,
  component,
  n_umi,
  keep = "tcell",
  contam = "platelet"
) {
  g <- pixelatorR::CellGraphData(cg, "cellgraph")
  compartment <- igraph::vertex_attr(g, "compartment")
  levels_all <- unique(c(keep, contam, "interface", "other"))
  n_by_comp <- table(factor(compartment, levels = levels_all))

  contam_idx <- which(compartment == contam)
  patch_sizes <- if (length(contam_idx) > 0) {
    igraph::components(igraph::induced_subgraph(g, contam_idx))$csize
  } else {
    integer(0)
  }

  tibble::tibble(
    component = component,
    n_nodes = length(compartment),
    n_umi_before = n_umi,
    n_nodes_keep = unname(n_by_comp[keep]),
    n_nodes_contam = unname(n_by_comp[contam]),
    n_nodes_interface = unname(n_by_comp["interface"]),
    n_nodes_other = unname(n_by_comp["other"]),
    n_patches = length(patch_sizes),
    largest_patch = if (length(patch_sizes) > 0) max(patch_sizes) else 0L,
    # Fraction of the cell covered by contaminating patches: the test in module
    # 02 divided the patch sizes by n_umi; the node-based version is kept too
    contam_fraction_umi = sum(patch_sizes) / n_umi,
    contam_fraction_nodes = sum(patch_sizes) / length(compartment)
  )
}

#' Process one batch: segment, QC, filter, recompute counts and proximity scores
#'
#' Results are written to `out_dir/batches/<batch_id>_*.rds`.
#'
#' @param batch One row of the manifest.
#' @param pg_data Seurat object that is the source of the cells and cell graphs.
#' @param w Protein weight matrix from `cc_protein_weights()`.
#' @param max_contam_fraction Cells with a contaminating fraction at or above
#'   this value are removed.
#' @param fraction_basis "umi" or "nodes": denominator of the fraction.
process_batch <- function(
  batch,
  pg_data,
  w,
  out_dir,
  keep = "tcell",
  contam = "platelet",
  max_contam_fraction = 0.2,
  fraction_basis = c("umi", "nodes"),
  detect_interface = TRUE,
  min_comp_size = 50,
  verbose = TRUE,
  step = ""
) {
  fraction_basis <- match.arg(fraction_basis)
  batch_id <- batch$batch_id
  cells <- batch$cells[[1]]
  t0 <- Sys.time()

  # Progress message prefixed with the position in the run, e.g. "[3/32] "
  say <- function(msg) {
    if (verbose) cli::cli_alert_info("{step}{batch_id}: {msg}")
  }

  # Subsetting an object to cells it does not contain fails with an obscure
  # error when the object holds graphs, so check the cells first
  missing_cells <- setdiff(cells, colnames(pg_data))
  if (length(missing_cells) > 0) {
    stop(
      batch_id,
      ": ",
      length(missing_cells),
      " of ",
      length(cells),
      " cells of this batch are not in the object passed as `pg_data`. ",
      "Subset the manifest to the same samples as the object."
    )
  }

  # 1. Load the cell graphs of the batch
  say(paste0(
    "loading cell graphs of ",
    length(cells),
    " cells (",
    batch$sample,
    ")..."
  ))
  se <- subset(pg_data, cells = cells)
  se <- pixelatorR::LoadCellGraphs(se, verbose = FALSE)
  cg_list <- pixelatorR::CellGraphs(se)
  se <- pixelatorR::RemoveCellGraphs(se)
  cg_list <- cg_list[cells] # redundant?
  n_umi <- stats::setNames(se$n_umi, colnames(se))[cells]

  # 2. Segment every graph; a failing graph is recorded and removed, not fatal
  say(paste0("segmenting ", length(cells), " cells..."))
  seg <- lapply(cells, function(nm) {
    tryCatch(
      pixelatorR::segment_cell(
        cg_list[[nm]],
        w,
        detect_interface = detect_interface,
        keep_largest_comp = FALSE,
        min_comp_size = min_comp_size,
        verbose = FALSE
      ),
      error = function(e) {
        structure(list(), class = "segment_error", msg = conditionMessage(e))
      }
    )
  }) |>
    stats::setNames(cells)
  rm(cg_list)
  failed <- vapply(seg, inherits, logical(1), what = "segment_error")

  # 3. One QC row per cell (cells that failed get status "segment_error")
  qc <- dplyr::bind_rows(lapply(cells, function(nm) {
    if (failed[[nm]]) {
      return(tibble::tibble(
        component = nm,
        n_umi_before = n_umi[[nm]],
        status = "segment_error"
      ))
    }
    cell_compartment_qc(seg[[nm]], nm, n_umi[[nm]], keep, contam) |>
      dplyr::mutate(status = "ok")
  }))
  stopifnot(setequal(qc$component, cells)) # no cell may go missing
  qc$contam_fraction <- if (fraction_basis == "umi") {
    qc$contam_fraction_umi
  } else {
    qc$contam_fraction_nodes
  }
  qc$batch_id <- batch_id
  qc$sample <- batch$sample
  qc$keep_cell <- qc$status == "ok" &
    !is.na(qc$contam_fraction) &
    qc$contam_fraction < max_contam_fraction &
    qc$n_nodes_keep > 0
  qc <- qc[match(cells, qc$component), ]

  keep_cells <- qc$component[qc$keep_cell]

  # 4. Counts of the kept compartment (before the graphs are reduced)
  # partition_counts() only returns the markers present in a graph, so every
  # cell is aligned to the full marker panel (missing markers = 0)
  features <- rownames(se)
  counts_list <- lapply(keep_cells, function(nm) {
    pc <- pixelatorR::partition_counts(
      seg[[nm]],
      partition_column = "compartment"
    )
    v <- stats::setNames(numeric(length(features)), features)
    present <- intersect(colnames(pc), features)
    v[present] <- as.numeric(as.matrix(pc[keep, present, drop = FALSE]))
    out <- Matrix::Matrix(
      v,
      ncol = 1,
      sparse = TRUE,
      dimnames = list(features, nm)
    )
    out
  })
  counts <- if (length(counts_list) > 0) {
    do.call(cbind, counts_list)
  } else {
    NULL
  }

  # 5. Reduce the kept graphs to the kept compartment
  cg_keep <- lapply(keep_cells, function(nm) {
    cg <- seg[[nm]]
    g <- pixelatorR::CellGraphData(cg, "cellgraph")
    nodes_keep <- igraph::V(g)$name[
      igraph::vertex_attr(g, "compartment") == keep
    ]
    subset(cg, nodes = nodes_keep)
  }) |>
    stats::setNames(keep_cells)
  rm(seg)

  # 6. Proximity scores for the filtered graphs
  prox <- NULL
  if (length(keep_cells) > 0) {
    say(paste0(
      "recomputing proximity scores for ",
      length(keep_cells),
      " kept cells (",
      length(cells) - length(keep_cells),
      " removed)..."
    ))
    assay <- pixelatorR::CreatePNAAssay5(
      counts = as.matrix(counts),
      cellgraphs = cg_keep
    )
    prox <- pixelatorR::ComputeProximityScores(assay, mode = "analytical")
    rm(assay)
    qc$n_umi_after <- NA_real_
    qc$n_umi_after[match(keep_cells, qc$component)] <- Matrix::colSums(counts)
  }

  # 7. Write the batch
  paths <- batch_paths(out_dir, batch_id)
  dir.create(dirname(paths[[1]]), recursive = TRUE, showWarnings = FALSE)
  save_atomic(counts, paths[["counts"]])
  save_atomic(prox, paths[["proximity"]])
  save_atomic(cg_keep, paths[["cellgraphs"]])
  save_atomic(qc, paths[["qc"]]) # written last: a finished batch has all files

  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  rss_gb <- round(process_rss_bytes() / 1e9, 1)
  if (verbose) {
    cli::cli_alert_success(
      "{step}{batch_id} ({batch$sample}): {length(keep_cells)}/{length(cells)} cells kept, {round(elapsed, 1)} min, memory {rss_gb} GB"
    )
  }
  rm(cg_keep, counts, prox, se)
  invisible(gc(verbose = FALSE))

  tibble::tibble(
    batch_id = batch_id,
    sample = batch$sample,
    n_in = length(cells),
    n_kept = length(keep_cells),
    n_failed = sum(failed),
    minutes = elapsed,
    rss_gb = rss_gb,
    status = "ok",
    error = NA_character_,
    finished = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  )
}

#' Format minutes as e.g. "45 s", "12.3 min" or "2.5 h"
format_minutes <- function(minutes) {
  if (is.na(minutes)) {
    "unknown"
  } else if (minutes < 1) {
    paste0(round(minutes * 60), " s")
  } else if (minutes < 90) {
    paste0(round(minutes, 1), " min")
  } else {
    paste0(round(minutes / 60, 1), " h")
  }
}

#' Run all batches of a manifest, skipping batches that are already finished
#'
#' Progress, an estimate of the remaining time and a final summary are printed.
#' A failing batch is logged and the run continues with the next batch (it is
#' retried the next time, since only finished batches are skipped); set
#' `stop_on_error = TRUE` to stop at the first failure.
#'
#' @param n_max Maximum number of batches to run in this call.
#' @param stop_on_error Stop at the first failing batch.
#' @param ... Passed on to `process_batch()`.
#' @return The per-batch results of this call (invisibly), also appended to
#'   `out_dir/batch_log.csv`.
run_batches <- function(
  manifest,
  pg_data,
  w,
  out_dir,
  ...,
  n_max = Inf,
  stop_on_error = FALSE
) {
  dir.create(
    file.path(out_dir, "batches"),
    recursive = TRUE,
    showWarnings = FALSE
  )
  is_done <- vapply(
    manifest$batch_id,
    batch_is_done,
    logical(1),
    out_dir = out_dir
  )
  todo <- utils::head(manifest[!is_done, ], n_max)
  n_todo <- nrow(todo)

  cli::cli_alert_success(
    "Platelet clean-up: {nrow(manifest)} batches ({format(sum(manifest$n_cells), big.mark = ',')} cells, {length(unique(manifest$sample))} samples); {sum(is_done)} already done, {n_todo} to run"
  )
  if (n_todo == 0) {
    cli::cli_alert_success("Nothing to do: all requested batches are finished.")
    return(invisible(NULL))
  }

  log_path <- file.path(out_dir, "batch_log.csv")
  t_start <- Sys.time()
  cells_done <- 0
  results <- vector("list", n_todo)

  for (i in seq_len(n_todo)) {
    batch <- todo[i, ]
    prefix <- paste0("[", i, "/", n_todo, "] ")

    res <- tryCatch(
      process_batch(batch, pg_data, w, out_dir, step = prefix, ...),
      error = function(e) {
        if (stop_on_error) {
          stop(e)
        }
        cli::cli_alert_danger(
          "{prefix}{batch$batch_id} ({batch$sample}) failed: {conditionMessage(e)}"
        )
        # Same columns as a successful batch, so that the log stays consistent
        tibble::tibble(
          batch_id = batch$batch_id,
          sample = batch$sample,
          n_in = batch$n_cells,
          n_kept = NA_integer_,
          n_failed = NA_integer_,
          minutes = NA_real_,
          rss_gb = NA_real_,
          status = "error",
          error = conditionMessage(e),
          finished = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
        )
      }
    )
    results[[i]] <- res

    # Append to the log (rewritten as a whole so the columns always match)
    # (all columns as text, so that types never conflict between runs)
    old_log <- if (file.exists(log_path)) {
      readr::read_csv(
        log_path,
        col_types = readr::cols(.default = readr::col_character())
      )
    }
    new_row <- dplyr::mutate(
      res,
      dplyr::across(dplyr::everything(), as.character)
    )
    readr::write_csv(dplyr::bind_rows(old_log, new_row), log_path)

    # Progress and estimate of the remaining time from cells per minute so far
    if (res$status == "ok") {
      cells_done <- cells_done + batch$n_cells
    }
    elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
    cells_left <- sum(todo$n_cells[seq_len(n_todo) > i])
    if (i < n_todo && cells_done > 0) {
      eta <- cells_left / (cells_done / elapsed)
      cli::cli_text(
        "   Elapsed {format_minutes(elapsed)}, about {format_minutes(eta)} left ({round(cells_done / elapsed)} cells/min)"
      )
    }
  }

  # Summary
  results <- dplyr::bind_rows(results)
  n_ok <- sum(results$status == "ok")
  n_err <- sum(results$status == "error")
  total <- format_minutes(as.numeric(difftime(
    Sys.time(),
    t_start,
    units = "mins"
  )))
  cells_in <- sum(results$n_in[results$status == "ok"])
  cells_kept <- sum(results$n_kept, na.rm = TRUE)

  if (n_err == 0) {
    cli::cli_alert_success(
      "Finished: {n_ok}/{n_todo} batches succeeded in {total}. {format(cells_kept, big.mark = ',')} of {format(cells_in, big.mark = ',')} cells kept."
    )
  } else {
    cli::cli_alert_danger(
      "Finished with errors: {n_ok}/{n_todo} batches succeeded, {n_err} failed in {total}. Failed: {paste(results$batch_id[results$status == 'error'], collapse = ', ')}. Rerun to retry them (finished batches are skipped)."
    )
  }

  n_left <- sum(
    !vapply(
      manifest$batch_id,
      batch_is_done,
      logical(1),
      out_dir = out_dir
    )
  )
  if (n_left == 0) {
    cli::cli_alert_success(
      "All {nrow(manifest)} batches are finished: run the assembly chunk."
    )
  } else {
    cli::cli_alert_info(
      "{n_left} of {nrow(manifest)} batches are not finished yet."
    )
  }

  invisible(results)
}

#' Combine the batch results into one cleaned Seurat object (without cell graphs)
#'
#' @param meta Cell metadata of the source object (`pg_data[[]]`).
#' @param assay_name Name of the assay in the new object (use the assay name of
#'   the source object, e.g. `DefaultAssay(pg_data)`).
#' @return A list with `object` (Seurat), `qc` and `n_batches`.
assemble_cleaned_object <- function(
  manifest,
  out_dir,
  meta,
  assay_name = "PNA"
) {
  done <- vapply(
    manifest$batch_id,
    batch_is_done,
    logical(1),
    out_dir = out_dir
  )
  if (!all(done)) {
    stop(
      sum(!done),
      " batches are not finished: ",
      paste(manifest$batch_id[!done], collapse = ", ")
    )
  }
  read_all <- function(type) {
    lapply(manifest$batch_id, function(id) {
      readRDS(batch_paths(out_dir, id)[[type]])
    })
  }

  counts <- Filter(Negate(is.null), read_all("counts"))
  counts <- do.call(cbind, counts)
  prox <- dplyr::bind_rows(read_all("proximity"))
  qc <- dplyr::bind_rows(read_all("qc"))

  assay <- pixelatorR::CreatePNAAssay5(
    counts = as.matrix(counts),
    cellgraphs = stats::setNames(
      rep(list(NULL), ncol(counts)),
      colnames(counts)
    )
  )
  pixelatorR::ProximityScores(assay) <- prox
  obj <- SeuratObject::CreateSeuratObject(assay, assay = assay_name)

  meta_add <- meta[
    colnames(obj),
    setdiff(colnames(meta), colnames(obj[[]])),
    drop = FALSE
  ]
  qc_add <- as.data.frame(qc[match(colnames(obj), qc$component), ])
  rownames(qc_add) <- colnames(obj)
  qc_add <- qc_add[,
    setdiff(
      colnames(qc_add),
      c(colnames(obj[[]]), colnames(meta_add), "component", "sample")
    ),
    drop = FALSE
  ]
  obj <- SeuratObject::AddMetaData(obj, cbind(meta_add, qc_add))

  list(object = obj, qc = qc, n_batches = nrow(manifest))
}
