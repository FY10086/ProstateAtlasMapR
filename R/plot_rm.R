#' Histogram of RM prediction scores
#'
#' @param query Seurat object after [map_query()].
#' @param annotation_level `"H1"`, `"H2"`, or `"H3"`.
#' @param breaks Number of histogram bins. Default `50`.
#'
#' @return A ggplot object.
#' @export
plot_rm_score_hist <- function(query,
                               annotation_level = c("H1", "H2", "H3"),
                               breaks = 50) {
  .check_ggplot()
  annotation_level <- .match_level(annotation_level)
  score_col <- .rm_score_col(annotation_level)
  if (!score_col %in% colnames(query[[]])) {
    cli::cli_abort("Column {.field {score_col}} not found. Run {.fn map_query} first.")
  }
  scores <- query[[score_col, drop = TRUE]]
  ggplot2::ggplot(data.frame(score = scores), ggplot2::aes(score)) +
    ggplot2::geom_histogram(bins = breaks, fill = "grey80", colour = "white") +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::labs(
      title = paste0(annotation_level, " prediction score"),
      x = score_col,
      y = "Number of cells"
    )
}

#' Feature plot of RM prediction scores on the projected UMAP
#'
#' @param query Seurat object after [map_query()].
#' @param annotation_level `"H1"`, `"H2"`, or `"H3"`.
#' @param reduction Dimensional reduction to plot. Default `"umap.rpca"`.
#'
#' @return A ggplot object (Seurat FeaturePlot).
#' @export
plot_rm_score <- function(query,
                          annotation_level = c("H1", "H2", "H3"),
                          reduction = "umap.rpca") {
  .check_ggplot()
  annotation_level <- .match_level(annotation_level)
  score_col <- .rm_score_col(annotation_level)
  if (!score_col %in% colnames(query[[]])) {
    cli::cli_abort("Column {.field {score_col}} not found. Run {.fn map_query} first.")
  }
  if (!reduction %in% Seurat::Reductions(query)) {
    cli::cli_abort("Reduction {.val {reduction}} not found in query.")
  }
  Seurat::FeaturePlot(
    query,
    features = score_col,
    reduction = reduction,
    label = FALSE
  )
}

#' DimPlot of predicted labels with low-score cells as Unknown
#'
#' Cells with prediction score below `score_threshold` are labelled
#' `"Unknown"` (grey). Set `score_threshold = NULL` to skip thresholding.
#'
#' **H1 / H2:** compact Seurat `DimPlot` with legend on the right
#' (`repr.plot` 12 x 10).
#' **H3:** same style as `09-RM.ipynb` — UMAP and legend are split with
#' `cowplot::plot_grid()` (`repr.plot` 22 x 10) so long subtype names
#' do not crush the map.
#'
#' @param query Seurat object after [map_query()].
#' @param annotation_level `"H1"`, `"H2"`, or `"H3"`.
#' @param score_threshold Score cutoff for Unknown. Default `0.6`.
#'   Use `NULL` to keep all predicted labels.
#' @param reduction Dimensional reduction to plot. Default `"umap.rpca"`.
#'
#' @return A ggplot object (H1/H2) or a cowplot / ggplot grob (H3).
#' @export
plot_rm_prediction <- function(query,
                               annotation_level = c("H1", "H2", "H3"),
                               score_threshold = 0.6,
                               reduction = "umap.rpca") {
  .check_ggplot()
  annotation_level <- .match_level(annotation_level)
  lab_col <- .rm_label_col(annotation_level)
  score_col <- .rm_score_col(annotation_level)

  if (!lab_col %in% colnames(query[[]])) {
    cli::cli_abort("Column {.field {lab_col}} not found. Run {.fn map_query} first.")
  }
  if (!reduction %in% Seurat::Reductions(query)) {
    cli::cli_abort("Reduction {.val {reduction}} not found in query.")
  }

  labs <- as.character(query[[lab_col, drop = TRUE]])
  if (!is.null(score_threshold)) {
    if (!score_col %in% colnames(query[[]])) {
      cli::cli_abort("Column {.field {score_col}} not found (needed for score_threshold).")
    }
    scores <- query[[score_col, drop = TRUE]]
    labs[!is.na(scores) & scores < score_threshold] <- "Unknown"
  }
  labs[is.na(labs)] <- "Unknown"

  types_main <- sort(setdiff(unique(labs), "Unknown"))
  cols <- .celltype_cols(c(types_main, "Unknown"))
  cols <- cols[intersect(names(cols), c(types_main, "Unknown"))]
  plot_col <- paste0(".patlas_", lab_col, "_plot")
  ## local copy so we do not leave a permanent meta column on the caller's object
  query_plot <- query
  query_plot[[plot_col]] <- factor(labs, levels = c(types_main, "Unknown"))

  title <- if (!is.null(score_threshold)) {
    paste0(lab_col, "  (score < ", score_threshold, " → Unknown)")
  } else {
    lab_col
  }

  ## Jupyter / IRkernel display size (also used as a layout hint)
  .set_rm_repr_size(annotation_level)

  if (identical(annotation_level, "H3")) {
    return(.plot_rm_prediction_h3(query_plot, plot_col, cols, reduction, title))
  }

  Seurat::DimPlot(
    query_plot,
    group.by = plot_col,
    reduction = reduction,
    label = FALSE
  ) +
    ggplot2::scale_color_manual(values = cols, na.value = "grey80", drop = FALSE) +
    ggplot2::ggtitle(title) +
    ggplot2::theme(aspect.ratio = 1)
}

#' Set IRkernel / Jupyter figure size for RM prediction plots
#' @noRd
.set_rm_repr_size <- function(annotation_level) {
  if (identical(annotation_level, "H3")) {
    options(repr.plot.width = 22, repr.plot.height = 10)
  } else {
    ## H1 / H2
    options(repr.plot.width = 12, repr.plot.height = 10)
  }
  invisible(TRUE)
}

#' H3 layout: DimPlot | legend via cowplot (09-RM.ipynb)
#' @noRd
.plot_rm_prediction_h3 <- function(query, group_var, cols, reduction, title) {
  if (!requireNamespace("cowplot", quietly = TRUE)) {
    cli::cli_abort(
      "Package {.pkg cowplot} is required for H3 {.fn plot_rm_prediction}. Install it with {.code install.packages(\"cowplot\")}."
    )
  }

  n_main <- sum(names(cols) != "Unknown")
  legend_ncol <- if (n_main + 1L <= 30L) 1L else 2L
  legend_text_size <- if (legend_ncol == 1L) 13 else 15
  legend_pt_size <- if (legend_ncol == 1L) 4 else 4.5
  ## tuned for repr.plot.width = 22
  rel_w <- if (legend_ncol == 1L) c(1, 0.85) else c(1, 0.95)

  p_main <- Seurat::DimPlot(
    query,
    group.by = group_var,
    reduction = reduction,
    label = FALSE
  ) +
    ggplot2::scale_color_manual(values = cols, drop = FALSE) +
    ggplot2::ggtitle(title) +
    Seurat::NoLegend()

  ## Build legend from a minimal ggplot — Seurat DimPlot themes break
  ## cowplot::get_legend() text rendering in ggplot2 >= 3.5.
  lv <- names(cols)
  d_leg <- data.frame(
    grp = factor(lv, levels = lv),
    stringsAsFactors = FALSE
  )
  p_leg_src <- ggplot2::ggplot(d_leg, ggplot2::aes(x = 1, y = 1, colour = grp)) +
    ggplot2::geom_point(size = legend_pt_size) +
    ggplot2::scale_colour_manual(values = cols, drop = FALSE, name = NULL) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        ncol = legend_ncol,
        override.aes = list(size = legend_pt_size)
      )
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = "right",
      legend.text = ggplot2::element_text(size = legend_text_size, colour = "black"),
      legend.key = ggplot2::element_rect(fill = NA, colour = NA),
      legend.background = ggplot2::element_rect(fill = "white", colour = NA),
      legend.margin = ggplot2::margin(0, 4, 0, 4)
    )

  p_legend <- cowplot::get_legend(p_leg_src)

  cowplot::plot_grid(
    p_main,
    p_legend,
    ncol = 2,
    rel_widths = rel_w,
    align = "h",
    axis = "tb"
  )
}
