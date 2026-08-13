#' Build a lightweight, distributable reference from the prostate atlas
#'
#' Loads the full prostate single-cell atlas (typically `QC_fin.qs`),
#' attaches a projection-capable model to the existing UMAP reduction
#' (`umap.rpca` by default) **without altering its coordinates**, strips the
#' object down to the minimum needed for Seurat's anchor-based reference
#' mapping workflow, and optionally caches the lean result to disk. This
#' lean object -- not the original multi-gigabyte atlas -- is the artifact
#' meant to be shared/distributed (see [get_reference()]).
#'
#' @details
#' The atlas's `umap.rpca` reflects a round of post-hoc, low-quality-cell
#' filtering performed *after* the UMAP embedding was computed, so its
#' coordinates cannot be reproduced by simply re-running
#' [Seurat::RunUMAP()] on the current cell set. Instead, this function
#' attaches a `uwot` transform model built with `init = <existing
#' coordinates>` and `n_epochs = 0`, which returns the supplied coordinates
#' completely unchanged while still building the nearest-neighbor structure
#' required to embed brand-new query cells later on.
#'
#' @param qc_path Path to the source Seurat object (`.qs` file), typically
#'   the full prostate atlas `QC_fin.qs`.
#' @param reduction Name of the dimensional reduction with feature loadings
#'   used to anchor/project new query cells. Default `"rpca"`.
#' @param umap_reduction Name of the UMAP reduction to attach a projection
#'   model to; its existing coordinates are preserved exactly. Default
#'   `"umap.rpca"`.
#' @param npcs Number of dimensions of `reduction` to use. Default `50`.
#' @param n.neighbors,metric Only used to build the neighbor structure that
#'   supports projecting *new* query cells; they do not affect the existing
#'   reference cells' coordinates. Defaults match [Seurat::RunUMAP()].
#' @param annotation_levels Which annotation level(s) (`celltype_manual_H1`
#'   / `H2` / `H3`) to retain in the reference. All three are kept by
#'   default so a single reference file supports every level; the level to
#'   use for a given mapping call is chosen later in [map_query()].
#' @param keep_meta_cols Additional `meta.data` columns to retain besides
#'   the annotation columns (e.g. `"sample.ID"`). `NULL` by default.
#' @param assay Assay holding gene expression. Default `"RNA"`.
#' @param save_path If not `NULL`, save the resulting reference with
#'   `qs::qsave(..., preset = "high")` to this path.
#' @param verbose Print progress messages.
#'
#' @return An object of class `patlas_reference`.
#' @export
build_reference <- function(qc_path,
                             reduction = "rpca",
                             umap_reduction = "umap.rpca",
                             npcs = 50,
                             n.neighbors = 30,
                             metric = "cosine",
                             annotation_levels = c("H1", "H2", "H3"),
                             keep_meta_cols = NULL,
                             assay = "RNA",
                             save_path = NULL,
                             verbose = TRUE) {
  if (!file.exists(qc_path)) {
    cli::cli_abort("File not found: {.path {qc_path}}")
  }
  annotation_levels <- match.arg(annotation_levels, choices = .pa_levels, several.ok = TRUE)
  annotation_cols <- .annotation_col(annotation_levels)

  if (verbose) {
    cli::cli_alert_info("Reading {.path {qc_path}} (this can take a while for large atlases)...")
  }
  seu <- qs::qread(qc_path)
  if (!inherits(seu, "Seurat")) {
    cli::cli_abort("Object read from {.path {qc_path}} is not a Seurat object.")
  }

  missing_cols <- setdiff(annotation_cols, colnames(seu@meta.data))
  if (length(missing_cols) > 0) {
    cli::cli_abort("Missing annotation column(s) in meta.data: {.val {missing_cols}}")
  }
  for (col in annotation_cols) {
    n_na <- sum(is.na(seu@meta.data[[col]]))
    if (verbose && n_na > 0) {
      cli::cli_alert_warning("{n_na} cell(s) have NA in {.field {col}}; they will be dropped when this level is selected in {.fn map_query}.")
    }
  }

  if (verbose) cli::cli_alert_info("Ensuring assay {.val {assay}} layers are joined...")
  seu <- .ensure_joined(seu, assay = assay)

  if (!reduction %in% Seurat::Reductions(seu)) {
    cli::cli_abort("Reduction {.val {reduction}} not found. Available: {.val {Seurat::Reductions(seu)}}")
  }
  loadings <- Seurat::Loadings(seu[[reduction]])
  if (length(loadings) == 0 || nrow(loadings) == 0) {
    cli::cli_abort(c(
      "Reduction {.val {reduction}} has no feature loadings.",
      "i" = "{.fn Seurat::FindTransferAnchors} needs a reduction with loadings to project new query cells."
    ))
  }
  if (!umap_reduction %in% Seurat::Reductions(seu)) {
    cli::cli_abort("UMAP reduction {.val {umap_reduction}} not found. Available: {.val {Seurat::Reductions(seu)}}")
  }

  if (verbose) {
    cli::cli_alert_info("Attaching a projection model to {.val {umap_reduction}} (coordinates will not change)...")
  }
  seu <- .attach_umap_model(
    seu,
    reduction = reduction,
    umap_reduction = umap_reduction,
    npcs = npcs,
    n.neighbors = n.neighbors,
    metric = metric,
    verbose = verbose
  )
  saved_model <- Seurat::Misc(seu[[umap_reduction]], slot = "model")

  ref_features <- intersect(rownames(loadings), rownames(seu[[assay]]))
  if (verbose) {
    cli::cli_alert_info("Trimming reference to {length(ref_features)} anchor feature(s); dropping scale.data / graphs / unused reductions...")
  }
  seu_diet <- Seurat::DietSeurat(
    seu,
    layers = c("counts", "data"),
    features = ref_features,
    assays = assay,
    dimreducs = c(reduction, umap_reduction),
    graphs = NULL
  )
  ## DietSeurat() is not guaranteed to preserve Misc() on kept reductions;
  ## re-attach the projection model explicitly to be safe (Seurat warns
  ## about overwriting existing misc data here, which is expected/harmless).
  suppressWarnings(
    Seurat::Misc(seu_diet[[umap_reduction]], slot = "model") <- saved_model
  )

  meta_keep <- intersect(unique(c(annotation_cols, keep_meta_cols)), colnames(seu_diet@meta.data))
  seu_diet@meta.data <- seu_diet@meta.data[, meta_keep, drop = FALSE]

  ref <- new_patlas_reference(
    seurat_obj = seu_diet,
    annotation_levels = annotation_levels,
    reduction = reduction,
    umap_reduction = umap_reduction,
    npcs = npcs,
    source_path = qc_path
  )

  if (!is.null(save_path)) {
    if (verbose) cli::cli_alert_info("Saving lean reference to {.path {save_path}}...")
    qs::qsave(ref, save_path, preset = "high")
  }

  if (verbose) print(ref)

  ref
}
