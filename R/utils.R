#' Available annotation levels in the atlas
#' @noRd
.pa_levels <- c("H1", "H2", "H3")

#' Map an annotation level (H1/H2/H3) to its metadata column name
#' @noRd
.annotation_col <- function(level) {
  paste0("celltype_manual_", level)
}

#' Validate and normalize a single annotation level argument
#' @noRd
.match_level <- function(annotation_level) {
  annotation_level <- match.arg(annotation_level, choices = .pa_levels)
  annotation_level
}

#' Validate that a level is available in a given reference
#' @noRd
.check_level_available <- function(annotation_level, reference) {
  if (!annotation_level %in% reference$annotation_levels) {
    cli::cli_abort(c(
      "Annotation level {.val {annotation_level}} is not available in this reference.",
      "i" = "Available level(s): {.val {reference$annotation_levels}}"
    ))
  }
  invisible(TRUE)
}

#' Ensure an assay's layers are joined (Seurat v5 split-layer objects)
#' @noRd
.ensure_joined <- function(seu, assay = "RNA") {
  layers <- tryCatch(SeuratObject::Layers(seu[[assay]]), error = function(e) NULL)
  needs_join <- !is.null(layers) && any(grepl("^(counts|data)\\.", layers))
  if (needs_join) {
    seu <- SeuratObject::JoinLayers(seu, assay = assay)
  }
  seu
}

#' Attach a projection-capable UMAP model to an existing reduction without
#' moving any of the existing cells' coordinates.
#'
#' Uses `uwot::umap(..., init = <existing coordinates>, n_epochs = 0)`: with
#' zero optimization epochs, uwot returns the supplied `init` unchanged as
#' the final embedding, while still building the nearest-neighbor structure
#' needed to project brand-new points into that same space later on.
#'
#' @noRd
.attach_umap_model <- function(seu,
                                reduction,
                                umap_reduction,
                                npcs,
                                n.neighbors = 30,
                                metric = "cosine",
                                verbose = TRUE) {
  src_emb <- Seurat::Embeddings(seu[[reduction]])
  if (npcs > ncol(src_emb)) {
    cli::cli_abort("npcs ({npcs}) exceeds the number of dimensions in reduction {.val {reduction}} ({ncol(src_emb)}).")
  }
  src_emb <- src_emb[, seq_len(npcs), drop = FALSE]
  target_emb <- Seurat::Embeddings(seu[[umap_reduction]])

  common_cells <- intersect(rownames(src_emb), rownames(target_emb))
  if (length(common_cells) != ncol(seu)) {
    cli::cli_abort("Cells in {.val {reduction}} and {.val {umap_reduction}} do not match; cannot attach a projection model.")
  }
  src_emb <- src_emb[common_cells, , drop = FALSE]
  target_emb <- target_emb[common_cells, , drop = FALSE]

  model <- uwot::umap(
    X = src_emb,
    init = target_emb,
    n_epochs = 0,
    n_neighbors = n.neighbors,
    metric = metric,
    ret_model = TRUE,
    verbose = FALSE
  )

  diff <- max(abs(as.matrix(model$embedding) - as.matrix(target_emb)))
  if (diff > 1e-6) {
    cli::cli_warn(c(
      "Reconstructed UMAP embedding differs from the original {.val {umap_reduction}} coordinates.",
      "i" = "Max absolute difference: {diff} (expected ~0). Proceeding anyway."
    ))
  } else if (verbose) {
    cli::cli_alert_success("{.val {umap_reduction}} coordinates verified unchanged (max abs diff = {format(diff, scientific = TRUE)}).")
  }

  new_reduc <- Seurat::CreateDimReducObject(
    embeddings = as.matrix(target_emb),
    key = Seurat::Key(seu[[umap_reduction]]),
    assay = Seurat::DefaultAssay(seu[[umap_reduction]])
  )
  Seurat::Misc(new_reduc, slot = "model") <- model
  seu[[umap_reduction]] <- new_reduc
  seu
}

#' Construct a `patlas_reference` object
#' @noRd
new_patlas_reference <- function(seurat_obj,
                                  annotation_levels,
                                  reduction,
                                  umap_reduction,
                                  npcs,
                                  source_path) {
  structure(
    list(
      seurat_obj = seurat_obj,
      annotation_levels = annotation_levels,
      reduction = reduction,
      umap_reduction = umap_reduction,
      npcs = npcs,
      source_path = source_path,
      built_at = Sys.time(),
      pkg_version = tryCatch(
        as.character(utils::packageVersion("ProstateAtlasMapR")),
        error = function(e) NA_character_
      )
    ),
    class = "patlas_reference"
  )
}

#' Load a reference from a `patlas_reference` object or a path to one
#' @noRd
.load_reference <- function(reference) {
  if (is.character(reference)) {
    if (!file.exists(reference)) {
      cli::cli_abort("Reference file not found: {.path {reference}}")
    }
    reference <- qs::qread(reference)
  }
  if (!inherits(reference, "patlas_reference")) {
    cli::cli_abort("{.arg reference} must be a {.cls patlas_reference} object or a path to one (see {.fn build_reference}/{.fn get_reference}).")
  }
  reference
}

#' @export
print.patlas_reference <- function(x, ...) {
  cli::cli_h3("<patlas_reference>")
  cli::cli_ul(c(
    "Cells: {.val {ncol(x$seurat_obj)}}",
    "Features: {.val {nrow(x$seurat_obj)}}",
    "Annotation level(s): {.val {x$annotation_levels}}",
    "Anchor reduction: {.val {x$reduction}} ({.val {x$npcs}} dims)",
    "UMAP reduction: {.val {x$umap_reduction}}",
    "Built at: {.val {format(x$built_at)}}",
    "Source: {.val {x$source_path}}",
    "Package version: {.val {x$pkg_version}}"
  ))
  invisible(x)
}
