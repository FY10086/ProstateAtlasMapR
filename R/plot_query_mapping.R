#' Visualize a mapped query alongside the reference atlas
#'
#' Draws the reference cells (colored by the chosen annotation level) and
#' the query cells (colored by their predicted labels) on the same UMAP
#' coordinate system produced by [map_query()], for a quick visual sanity
#' check of the mapping result.
#'
#' @param reference A `patlas_reference` object, or a path to one.
#' @param query A query `Seurat` object that has already been processed by
#'   [map_query()] (i.e. it has a reduction named `reference$umap_reduction`
#'   and a `predicted.<annotation column>` metadata column).
#' @param annotation_level Which annotation level was used for mapping. One
#'   of `"H1"`, `"H2"`, `"H3"`.
#' @param pt.size.reference,pt.size.query Point sizes for the reference and
#'   query cells, passed to [Seurat::DimPlot()].
#' @param ... Additional arguments passed to [Seurat::DimPlot()].
#'
#' @return A `patchwork`/`ggplot` object with two side-by-side panels:
#'   reference cells and query cells on the same UMAP coordinates.
#' @export
plot_query_mapping <- function(reference,
                                query,
                                annotation_level = c("H1", "H2", "H3"),
                                pt.size.reference = 0.1,
                                pt.size.query = 0.6,
                                ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE) || !requireNamespace("patchwork", quietly = TRUE)) {
    cli::cli_abort("Packages {.pkg ggplot2} and {.pkg patchwork} are required for {.fn plot_query_mapping}.")
  }
  reference <- .load_reference(reference)
  annotation_level <- .match_level(annotation_level)
  .check_level_available(annotation_level, reference)
  annotation_col <- .annotation_col(annotation_level)
  predicted_col <- paste0("predicted.", annotation_col)

  if (!predicted_col %in% colnames(query@meta.data)) {
    cli::cli_abort("Column {.field {predicted_col}} not found in query; did you run {.fn map_query} with {.code annotation_level = \"{annotation_level}\"} first?")
  }
  if (!reference$umap_reduction %in% Seurat::Reductions(query)) {
    cli::cli_abort("Reduction {.field {reference$umap_reduction}} not found in query; did you run {.fn map_query} first?")
  }

  p1 <- Seurat::DimPlot(
    reference$seurat_obj,
    reduction = reference$umap_reduction,
    group.by = annotation_col,
    pt.size = pt.size.reference,
    ...
  ) + ggplot2::ggtitle(paste0("Reference (", annotation_col, ")"))

  p2 <- Seurat::DimPlot(
    query,
    reduction = reference$umap_reduction,
    group.by = predicted_col,
    pt.size = pt.size.query,
    ...
  ) + ggplot2::ggtitle(paste0("Query (", predicted_col, ")"))

  patchwork::wrap_plots(p1, p2)
}
