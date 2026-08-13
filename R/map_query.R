#' Map a query Seurat object onto the prostate atlas reference
#'
#' Runs Seurat's standard anchor-based reference-mapping workflow
#' ([Seurat::FindTransferAnchors()] + [Seurat::MapQuery()]) to transfer cell
#' type labels from the reference atlas onto a new (query) single-cell
#' dataset, and to project the query cells into the reference's published
#' UMAP coordinate system.
#'
#' @param reference A `patlas_reference` object (from [build_reference()] or
#'   [get_reference()]), or a path to a `.qs` file containing one.
#' @param query A `Seurat` object to map onto the reference.
#' @param annotation_level Which annotation level to transfer: one of
#'   `"H1"`, `"H2"`, `"H3"`. Must be one of `reference$annotation_levels`.
#' @param query_assay Assay in `query` holding gene expression. If `NULL`
#'   (default), uses `"RNA"` when present, otherwise [Seurat::DefaultAssay()].
#'   For Seurat sketch objects this is often `"sketch"`.
#' @param dims Number of reference dimensions to use for anchoring. Defaults
#'   to `reference$npcs`.
#' @param normalization.method Passed to [Seurat::FindTransferAnchors()];
#'   `query` is normalized automatically with this method if it has not
#'   been already.
#' @param k.anchor,k.score Passed to [Seurat::FindTransferAnchors()].
#' @param k.weight Passed to [Seurat::MapQuery()] / `TransferData()`.
#' @param verbose Print progress messages.
#'
#' @return The `query` Seurat object with two new columns added to
#'   `meta.data` (`predicted.celltype_manual_H{level}` and
#'   `predicted.celltype_manual_H{level}.score`), and a new dimensional
#'   reduction -- named the same as `reference$umap_reduction` (e.g.
#'   `"umap.rpca"`) -- holding the query cells projected into the
#'   reference's published UMAP coordinate system.
#' @export
map_query <- function(reference,
                       query,
                       annotation_level = c("H1", "H2", "H3"),
                       query_assay = NULL,
                       dims = NULL,
                       normalization.method = "LogNormalize",
                       k.anchor = 5,
                       k.score = 30,
                       k.weight = 100,
                       verbose = TRUE) {
  reference <- .load_reference(reference)
  annotation_level <- .match_level(annotation_level)
  .check_level_available(annotation_level, reference)
  annotation_col <- .annotation_col(annotation_level)

  if (!inherits(query, "Seurat")) {
    cli::cli_abort("{.arg query} must be a Seurat object.")
  }
  if (is.null(query_assay)) {
    query_assay <- if ("RNA" %in% Seurat::Assays(query)) {
      "RNA"
    } else {
      Seurat::DefaultAssay(query)
    }
    if (verbose) {
      cli::cli_alert_info("Using query assay {.val {query_assay}}.")
    }
  }
  if (!query_assay %in% Seurat::Assays(query)) {
    cli::cli_abort("Assay {.val {query_assay}} not found in query. Available: {.val {Seurat::Assays(query)}}")
  }

  ref_seu <- reference$seurat_obj
  if (is.null(dims)) dims <- reference$npcs

  keep_cells <- colnames(ref_seu)[!is.na(ref_seu@meta.data[[annotation_col]])]
  n_dropped <- ncol(ref_seu) - length(keep_cells)
  if (n_dropped > 0) {
    if (verbose) {
      cli::cli_alert_warning("Excluding {n_dropped} reference cell(s) with NA in {.field {annotation_col}}.")
    }
    ref_seu <- subset(ref_seu, cells = keep_cells)
  }

  Seurat::DefaultAssay(query) <- query_assay
  ## Seurat v5 split layers (data.GSE*, counts.GSE*, ...) must be joined,
  ## otherwise LayerData()/FindTransferAnchors may use only the first layer.
  query <- .ensure_joined(query, assay = query_assay)

  data_layer <- tryCatch(
    SeuratObject::LayerData(query, assay = query_assay, layer = "data"),
    error = function(e) NULL
  )
  if (is.null(data_layer) || length(data_layer) == 0) {
    if (verbose) cli::cli_alert_info("Normalizing query data ({.val {normalization.method}})...")
    query <- Seurat::NormalizeData(query, normalization.method = normalization.method, verbose = FALSE)
  }

  if (verbose) {
    cli::cli_alert_info("Finding transfer anchors (reference.reduction = {.val {reference$reduction}}, dims = 1:{dims})...")
  }
  anchors <- Seurat::FindTransferAnchors(
    reference = ref_seu,
    query = query,
    reference.assay = SeuratObject::DefaultAssay(ref_seu),
    query.assay = query_assay,
    reference.reduction = reference$reduction,
    dims = seq_len(dims),
    normalization.method = normalization.method,
    k.anchor = k.anchor,
    k.score = k.score,
    verbose = verbose
  )

  refdata <- stats::setNames(list(annotation_col), annotation_col)

  if (verbose) {
    cli::cli_alert_info("Transferring {.field {annotation_col}} labels and projecting onto {.val {reference$umap_reduction}}...")
  }
  ## By default Seurat::ProjectUMAP() names the query's projected UMAP
  ## "ref.umap"; we force it to reuse the reference's own reduction name/key
  ## (e.g. "umap.rpca") so reference and query embeddings can be combined
  ## directly for plotting.
  umap_key <- Seurat::Key(ref_seu[[reference$umap_reduction]])
  query <- Seurat::MapQuery(
    anchorset = anchors,
    query = query,
    reference = ref_seu,
    refdata = refdata,
    reference.reduction = reference$reduction,
    reduction.model = reference$umap_reduction,
    transferdata.args = list(k.weight = k.weight),
    projectumap.args = list(
      reduction.name = reference$umap_reduction,
      reduction.key = umap_key
    ),
    verbose = verbose
  )

  query
}
