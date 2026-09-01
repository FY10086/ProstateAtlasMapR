#' Run RCTD spatial deconvolution against the prostate atlas
#'
#' Wraps [spacexr::create.RCTD()] + [spacexr::run.RCTD()] using the same
#' prostate atlas reference as [map_query()]. Choose which annotation level
#' (`H1` / `H2` / `H3`) to use as the RCTD cell-type labels.
#'
#' @details
#' The lean `patlas_reference` (from [build_reference()] / [get_reference()])
#' keeps the ~2,000 anchor genes used for Seurat mapping. That gene set is
#' usually enough for RCTD, but if you need a fuller transcriptome you can
#' pass a Seurat object (or `.qs` path) that still has `counts` plus the
#' `celltype_manual_H*` columns.
#'
#' Spatial coordinates resolution order:
#' 1. `coords` as a reduction name (e.g. `"spatial"`) → first two embedding
#'    dimensions;
#' 2. `coords` as a `data.frame`/`matrix` with columns `x`/`y`;
#' 3. otherwise [Seurat::GetTissueCoordinates()] (optionally via `image`).
#' Set `use_fake_coords = TRUE` only when you want proportions without a real
#' spatial layout (coords are still required by spacexr internally).
#'
#' @param reference A `patlas_reference`, a path to one, a Seurat object with
#'   atlas annotations, or a path to such a Seurat `.qs` file.
#' @param spatial A Seurat object holding spatial transcriptomics data
#'   (Visium / other spot-based assays with tissue coordinates).
#' @param annotation_level Which atlas level to use as RCTD cell types:
#'   `"H1"`, `"H2"`, or `"H3"`.
#' @param assay Assay in `spatial` with raw counts. If `NULL`, uses `"Spatial"`
#'   when present, otherwise [Seurat::DefaultAssay()].
#' @param image Name of the Seurat image to pull coordinates from; passed to
#'   [Seurat::GetTissueCoordinates()]. `NULL` uses the default image. Ignored
#'   when `coords` is supplied.
#' @param coords Spot coordinates. Preferred: a **reduction name** in
#'   `spatial` (e.g. `"spatial"`), using the first two dimensions as `x`/`y`.
#'   Alternatively a `data.frame`/`matrix` with columns `x` and `y` and
#'   rownames matching `colnames(spatial)`. If `NULL`, falls back to
#'   [Seurat::GetTissueCoordinates()].
#' @param use_fake_coords If `TRUE`, ignore real coordinates and use a
#'   placeholder grid (see [spacexr::SpatialRNA()]). Default `FALSE`.
#' @param doublet_mode RCTD mode passed to [spacexr::run.RCTD()]:
#'   `"doublet"` (default; Visium-like), `"full"`, or `"multi"`.
#' @param max_cores Parallel cores for [spacexr::create.RCTD()]. Default `4`.
#' @param n_max_cells Max cells kept per cell type when building the spacexr
#'   `Reference` (downsampling). Default `10000`.
#' @param min_UMI Minimum UMI for reference cells. Default `100`.
#' @param CELL_MIN_INSTANCE Minimum cells required per cell type in
#'   [spacexr::create.RCTD()]. Default `25`. Lower for rare H3 types if needed.
#' @param gene_cutoff,fc_cutoff,gene_cutoff_reg,fc_cutoff_reg,UMI_min,UMI_max,counts_MIN,UMI_min_sigma,CONFIDENCE_THRESHOLD,DOUBLET_THRESHOLD,MAX_MULTI_TYPES
#'   Passed through to [spacexr::create.RCTD()].
#' @param class_level Optional coarser annotation level (`"H1"` / `"H2"` /
#'   `"H3"`) used to build spacexr's `class_df` so RCTD can report
#'   class-level confidence. Must differ from `annotation_level` and be
#'   present in `reference`. `NULL` (default) skips `class_df`.
#' @param assay_name Name of the new assay written onto `spatial` that holds
#'   cell-type weights (features = cell types, cells = spots). Default
#'   `"RCTD"`. Set to `NULL` to skip writing an assay.
#' @param store_rctd If `TRUE` (default), store the full spacexr `RCTD`
#'   object in `spatial@misc$rctd`.
#' @param verbose Print progress messages.
#'
#' @return The `spatial` Seurat object with:
#'   * meta.data columns `rctd_spot_class`, `rctd_first_type`,
#'     `rctd_second_type` (doublet mode; `NA` when unavailable);
#'   * `@misc$rctd_weights`: spots x cell-types matrix (original type names);
#'   * an assay of weights named `assay_name` (when not `NULL`; `_` in type
#'     names is replaced by `-` to satisfy Seurat feature-name rules);
#'   * optionally the raw `RCTD` object in `@misc$rctd`.
#'
#' @export
run_rctd <- function(reference,
                     spatial,
                     annotation_level = c("H1", "H2", "H3"),
                     assay = NULL,
                     image = NULL,
                     coords = NULL,
                     use_fake_coords = FALSE,
                     doublet_mode = c("doublet", "full", "multi"),
                     max_cores = 4,
                     n_max_cells = 10000,
                     min_UMI = 100,
                     CELL_MIN_INSTANCE = 25,
                     gene_cutoff = 0.000125,
                     fc_cutoff = 0.5,
                     gene_cutoff_reg = 2e-04,
                     fc_cutoff_reg = 0.75,
                     UMI_min = 100,
                     UMI_max = 2e+07,
                     counts_MIN = 10,
                     UMI_min_sigma = 300,
                     CONFIDENCE_THRESHOLD = 5,
                     DOUBLET_THRESHOLD = 20,
                     MAX_MULTI_TYPES = 4,
                     class_level = NULL,
                     assay_name = "RCTD",
                     store_rctd = TRUE,
                     verbose = TRUE) {
  .check_spacexr()
  annotation_level <- .match_level(annotation_level)
  doublet_mode <- match.arg(doublet_mode)

  if (!inherits(spatial, "Seurat")) {
    cli::cli_abort("{.arg spatial} must be a Seurat object.")
  }

  ref_seu <- .as_reference_seurat(reference, annotation_level)
  annotation_col <- .annotation_col(annotation_level)

  if (is.null(assay)) {
    assay <- if ("Spatial" %in% Seurat::Assays(spatial)) {
      "Spatial"
    } else {
      Seurat::DefaultAssay(spatial)
    }
    if (verbose) cli::cli_alert_info("Using spatial assay {.val {assay}}.")
  }
  if (!assay %in% Seurat::Assays(spatial)) {
    cli::cli_abort("Assay {.val {assay}} not found in spatial. Available: {.val {Seurat::Assays(spatial)}}")
  }

  spatial <- .ensure_joined(spatial, assay = assay)
  ref_seu <- .ensure_joined(ref_seu, assay = SeuratObject::DefaultAssay(ref_seu))

  class_df <- .build_class_df(ref_seu, annotation_level, class_level)

  if (verbose) {
    cli::cli_alert_info(
      "Building spacexr Reference from {.field {annotation_col}} (n_max_cells = {n_max_cells})..."
    )
  }
  rctd_ref <- .seurat_to_rctd_reference(
    ref_seu,
    annotation_col = annotation_col,
    n_max_cells = n_max_cells,
    min_UMI = min_UMI
  )

  if (verbose) cli::cli_alert_info("Building spacexr SpatialRNA...")
  spatial_rna <- .seurat_to_spatial_rna(
    spatial,
    assay = assay,
    image = image,
    coords = coords,
    use_fake_coords = use_fake_coords
  )

  n_genes_ref <- nrow(rctd_ref@counts)
  if (verbose && n_genes_ref <= 3000) {
    cli::cli_alert_warning(
      "Reference has only {n_genes_ref} gene(s). Lean atlas references keep ~2k anchor genes; pass a fuller Seurat object if you want a larger RCTD gene universe."
    )
  }

  if (verbose) {
    cli::cli_alert_info("create.RCTD() (max_cores = {max_cores})...")
  }
  rctd <- spacexr::create.RCTD(
    spatialRNA = spatial_rna,
    reference = rctd_ref,
    max_cores = max_cores,
    gene_cutoff = gene_cutoff,
    fc_cutoff = fc_cutoff,
    gene_cutoff_reg = gene_cutoff_reg,
    fc_cutoff_reg = fc_cutoff_reg,
    UMI_min = UMI_min,
    UMI_max = UMI_max,
    counts_MIN = counts_MIN,
    UMI_min_sigma = UMI_min_sigma,
    class_df = class_df,
    CELL_MIN_INSTANCE = CELL_MIN_INSTANCE,
    MAX_MULTI_TYPES = MAX_MULTI_TYPES,
    CONFIDENCE_THRESHOLD = CONFIDENCE_THRESHOLD,
    DOUBLET_THRESHOLD = DOUBLET_THRESHOLD
  )

  if (verbose) {
    cli::cli_alert_info("run.RCTD(doublet_mode = {.val {doublet_mode}})...")
  }
  rctd <- spacexr::run.RCTD(rctd, doublet_mode = doublet_mode)

  spatial <- .attach_rctd_results(
    spatial,
    rctd,
    doublet_mode = doublet_mode,
    assay_name = assay_name,
    store_rctd = store_rctd,
    annotation_level = annotation_level,
    verbose = verbose
  )

  if (verbose) cli::cli_alert_success("RCTD finished.")
  spatial
}

#' Require spacexr without making it a hard Import
#' @noRd
.check_spacexr <- function() {
  if (!requireNamespace("spacexr", quietly = TRUE)) {
    cli::cli_abort(c(
      "Package {.pkg spacexr} is required for {.fn run_rctd}.",
      "i" = "Install with {.code install.packages(\"spacexr\")} or from https://github.com/dmcable/spacexr"
    ))
  }
  invisible(TRUE)
}

#' Coerce reference argument to a Seurat object with the requested annotation
#' @noRd
.as_reference_seurat <- function(reference, annotation_level) {
  if (inherits(reference, "patlas_reference")) {
    .check_level_available(annotation_level, reference)
    return(reference$seurat_obj)
  }
  if (is.character(reference)) {
    if (!file.exists(reference)) {
      cli::cli_abort("Reference file not found: {.path {reference}}")
    }
    reference <- qs::qread(reference)
    if (inherits(reference, "patlas_reference")) {
      .check_level_available(annotation_level, reference)
      return(reference$seurat_obj)
    }
  }
  if (inherits(reference, "Seurat")) {
    col <- .annotation_col(annotation_level)
    if (!col %in% colnames(reference@meta.data)) {
      cli::cli_abort(c(
        "Seurat reference is missing {.field {col}}.",
        "i" = "Expected atlas-style meta.data columns celltype_manual_H1/H2/H3."
      ))
    }
    return(reference)
  }
  cli::cli_abort(c(
    "{.arg reference} must be a {.cls patlas_reference}, a Seurat object, or a path to either.",
    "i" = "See {.fn get_reference} / {.fn build_reference}."
  ))
}

#' Build spacexr::Reference from a Seurat object
#' @noRd
.seurat_to_rctd_reference <- function(seu,
                                      annotation_col,
                                      n_max_cells = 10000,
                                      min_UMI = 100) {
  assay <- SeuratObject::DefaultAssay(seu)
  counts <- tryCatch(
    SeuratObject::LayerData(seu, assay = assay, layer = "counts"),
    error = function(e) {
      cli::cli_abort(c(
        "Could not extract {.field counts} from reference assay {.val {assay}}.",
        "x" = conditionMessage(e),
        "i" = "RCTD requires raw integer counts."
      ))
    }
  )
  counts <- methods::as(counts, "dgCMatrix")

  cell_types <- seu[[annotation_col, drop = TRUE]]
  names(cell_types) <- colnames(seu)
  keep <- !is.na(cell_types) & cell_types != ""
  n_drop <- sum(!keep)
  if (n_drop > 0) {
    cli::cli_alert_warning(
      "Excluding {n_drop} reference cell(s) with NA/empty {.field {annotation_col}}."
    )
  }
  cell_types <- cell_types[keep]
  counts <- counts[, names(cell_types), drop = FALSE]

  ## spacexr forbids '/' in cell-type names
  if (any(grepl("/", as.character(cell_types), fixed = TRUE))) {
    cli::cli_abort(c(
      "Cell type names in {.field {annotation_col}} contain '/' which spacexr forbids.",
      "i" = "Rename those labels before running {.fn run_rctd}."
    ))
  }

  cell_types <- factor(cell_types)
  spacexr::Reference(
    counts = counts,
    cell_types = cell_types,
    n_max_cells = n_max_cells,
    min_UMI = min_UMI
  )
}

#' Build spacexr::SpatialRNA from a Seurat spatial object
#' @noRd
.seurat_to_spatial_rna <- function(spatial,
                                   assay,
                                   image = NULL,
                                   coords = NULL,
                                   use_fake_coords = FALSE) {
  counts <- tryCatch(
    SeuratObject::LayerData(spatial, assay = assay, layer = "counts"),
    error = function(e) {
      cli::cli_abort(c(
        "Could not extract {.field counts} from spatial assay {.val {assay}}.",
        "x" = conditionMessage(e)
      ))
    }
  )
  counts <- methods::as(counts, "dgCMatrix")

  if (isTRUE(use_fake_coords)) {
    return(spacexr::SpatialRNA(coords = NULL, counts = counts, use_fake_coords = TRUE))
  }

  coords <- .get_spatial_coords(spatial, image = image, coords = coords)
  common <- intersect(rownames(coords), colnames(counts))
  if (length(common) == 0) {
    cli::cli_abort("No overlapping barcodes between spatial coordinates and counts.")
  }
  if (length(common) < ncol(counts)) {
    cli::cli_alert_warning(
      "Dropping {ncol(counts) - length(common)} spot(s) without coordinates."
    )
  }
  counts <- counts[, common, drop = FALSE]
  coords <- coords[common, , drop = FALSE]
  spacexr::SpatialRNA(coords = coords, counts = counts)
}

#' Normalize Seurat tissue coordinates to a coords data.frame with x/y
#'
#' `coords` may be:
#' - `NULL`: [Seurat::GetTissueCoordinates()]
#' - character length 1: a DimReduc name; first two dims → x/y
#' - data.frame/matrix with columns x and y
#'
#' @noRd
.get_spatial_coords <- function(spatial, image = NULL, coords = NULL) {
  if (!is.null(coords)) {
    if (is.character(coords)) {
      if (length(coords) != 1L || is.na(coords) || !nzchar(coords)) {
        cli::cli_abort("{.arg coords} reduction name must be a single non-empty string.")
      }
      if (!coords %in% Seurat::Reductions(spatial)) {
        cli::cli_abort(c(
          "Reduction {.val {coords}} not found in {.arg spatial}.",
          "i" = "Available: {.val {Seurat::Reductions(spatial)}}"
        ))
      }
      emb <- Seurat::Embeddings(spatial, reduction = coords)
      if (ncol(emb) < 2L) {
        cli::cli_abort(
          "Reduction {.val {coords}} has {ncol(emb)} dimension(s); need at least 2 for x/y."
        )
      }
      out <- data.frame(
        x = emb[, 1L],
        y = emb[, 2L],
        row.names = rownames(emb)
      )
      return(out)
    }

    if (!is.data.frame(coords) && !is.matrix(coords)) {
      cli::cli_abort(c(
        "{.arg coords} must be a reduction name, or a data.frame/matrix with columns x and y.",
        "i" = "Example: {.code coords = \"spatial\"}"
      ))
    }
    coords <- as.data.frame(coords)
    if (!all(c("x", "y") %in% colnames(coords))) {
      ## allow 2-column embeddings without renaming (e.g. coords_1 / coords_2)
      if (ncol(coords) >= 2L && all(vapply(coords[, 1:2, drop = FALSE], is.numeric, logical(1)))) {
        coords <- data.frame(
          x = coords[[1L]],
          y = coords[[2L]],
          row.names = rownames(coords)
        )
      } else {
        cli::cli_abort("{.arg coords} must have columns {.field x} and {.field y} (or two leading numeric columns).")
      }
    }
    if (is.null(rownames(coords)) || anyDuplicated(rownames(coords))) {
      cli::cli_abort("{.arg coords} must have unique rownames matching spot barcodes.")
    }
    return(coords[, c("x", "y"), drop = FALSE])
  }

  tc <- tryCatch(
    {
      if (is.null(image)) {
        Seurat::GetTissueCoordinates(spatial)
      } else {
        Seurat::GetTissueCoordinates(spatial, image = image)
      }
    },
    error = function(e) {
      cli::cli_abort(c(
        "Could not obtain tissue coordinates via {.fn Seurat::GetTissueCoordinates}.",
        "x" = conditionMessage(e),
        "i" = "Pass {.arg coords} as a reduction name (e.g. {.code \"spatial\"}), a data.frame with x/y, or set {.arg use_fake_coords = TRUE}."
      ))
    }
  )
  tc <- as.data.frame(tc)

  if ("cell" %in% colnames(tc) && (is.null(rownames(tc)) || !all(tc$cell %in% colnames(spatial)))) {
    rownames(tc) <- as.character(tc$cell)
  }
  if (is.null(rownames(tc)) || !any(rownames(tc) %in% colnames(spatial))) {
    cli::cli_abort(c(
      "Tissue coordinates lack recognizable barcode rownames.",
      "i" = "Pass {.arg coords} as a reduction name or a data.frame with rownames = spot barcodes."
    ))
  }

  if (all(c("x", "y") %in% colnames(tc))) {
    out <- tc[, c("x", "y"), drop = FALSE]
  } else if (all(c("imagecol", "imagerow") %in% colnames(tc))) {
    out <- data.frame(x = tc$imagecol, y = tc$imagerow, row.names = rownames(tc))
  } else if (ncol(tc) >= 2) {
    ## Seurat v5 Visium often returns first two numeric columns as spatial axes
    num_cols <- names(tc)[vapply(tc, is.numeric, logical(1))]
    if (length(num_cols) < 2) {
      cli::cli_abort("Could not find two numeric coordinate columns in GetTissueCoordinates() output.")
    }
    out <- data.frame(
      x = tc[[num_cols[1]]],
      y = tc[[num_cols[2]]],
      row.names = rownames(tc)
    )
  } else {
    cli::cli_abort("Unexpected GetTissueCoordinates() format; pass {.arg coords} explicitly.")
  }
  out
}

#' Optional class_df for spacexr from a coarser annotation level
#' @noRd
.build_class_df <- function(ref_seu, annotation_level, class_level) {
  if (is.null(class_level)) return(NULL)
  class_level <- .match_level(class_level)
  if (identical(class_level, annotation_level)) {
    cli::cli_abort("{.arg class_level} must differ from {.arg annotation_level}.")
  }
  fine_col <- .annotation_col(annotation_level)
  coarse_col <- .annotation_col(class_level)
  if (!coarse_col %in% colnames(ref_seu@meta.data)) {
    cli::cli_abort("class_level column {.field {coarse_col}} not found in reference meta.data.")
  }
  fine <- as.character(ref_seu[[fine_col, drop = TRUE]])
  coarse <- as.character(ref_seu[[coarse_col, drop = TRUE]])
  keep <- !is.na(fine) & !is.na(coarse) & fine != "" & coarse != ""
  mapping <- unique(data.frame(
    cell_type = fine[keep],
    class = coarse[keep],
    stringsAsFactors = FALSE
  ))
  ## if a fine type maps to multiple coarse classes, keep the majority vote
  if (anyDuplicated(mapping$cell_type)) {
    tab <- as.data.frame(table(
      cell_type = fine[keep],
      class = coarse[keep]
    ), stringsAsFactors = FALSE)
    tab <- tab[tab$Freq > 0, , drop = FALSE]
    tab <- tab[order(tab$cell_type, -tab$Freq), , drop = FALSE]
    mapping <- tab[!duplicated(tab$cell_type), c("cell_type", "class"), drop = FALSE]
    cli::cli_alert_warning(
      "Some {.field {fine_col}} labels map to multiple {.field {coarse_col}} classes; using majority vote for class_df."
    )
  }
  rownames(mapping) <- mapping$cell_type
  mapping
}

#' Write RCTD results back onto the spatial Seurat object
#' @noRd
.attach_rctd_results <- function(spatial,
                                 rctd,
                                 doublet_mode,
                                 assay_name,
                                 store_rctd,
                                 annotation_level,
                                 verbose = TRUE) {
  meta <- .rctd_results_meta(rctd, doublet_mode)
  common <- intersect(rownames(meta), colnames(spatial))
  ## Seurat v5: in-place meta.data[i] <- does not stick; assign full vectors.
  md <- spatial[[]]
  for (col in colnames(meta)) {
    vals <- rep(NA_character_, nrow(md))
    names(vals) <- rownames(md)
    if (length(common) > 0) {
      vals[common] <- as.character(meta[common, col])
    }
    md[[col]] <- unname(vals[rownames(md)])
    ## also store level-prefixed columns (rctd_H2_first_type, ...)
    pref <- sub("^rctd_", paste0("rctd_", annotation_level, "_"), col)
    if (!identical(pref, col)) {
      md[[pref]] <- unname(vals[rownames(md)])
    }
  }
  spatial[[]] <- md

  weights <- .rctd_weight_matrix(rctd, doublet_mode)
  if (!is.null(weights)) {
    ## authoritative matrix (spots x cell types) keeps original type names
    spatial@misc$rctd_weights <- weights
    ## level-prefixed copy so multi-level runs can be plotted selectively
    spatial@misc[[paste0("rctd_weights_", annotation_level)]] <- weights
  }
  if (!is.null(assay_name) && !is.null(weights)) {
    ## Seurat forbids '_' in feature names; sanitize for the assay only.
    w <- weights
    orig_types <- colnames(w)
    colnames(w) <- gsub("_", "-", orig_types, fixed = TRUE)
    if (any(orig_types != colnames(w))) {
      spatial@misc$rctd_weight_feature_map <- data.frame(
        assay_feature = colnames(w),
        cell_type = orig_types,
        stringsAsFactors = FALSE
      )
      if (verbose) {
        cli::cli_alert_warning(
          "Cell-type names contained '_'; assay {.val {assay_name}} features use '-' instead. Original names are in {.code spatial@misc$rctd_weights}."
        )
      }
    }
    w_assay <- Matrix::t(methods::as(w, "dgCMatrix"))
    spatial[[assay_name]] <- SeuratObject::CreateAssayObject(counts = w_assay)
    if (verbose) {
      cli::cli_alert_info(
        "Wrote weight assay {.val {assay_name}} ({nrow(w_assay)} cell type(s) x {ncol(w_assay)} spot(s))."
      )
    }
  }

  spatial@misc$rctd_annotation_level <- annotation_level
  spatial@misc$rctd_doublet_mode <- doublet_mode
  if (isTRUE(store_rctd)) {
    spatial@misc$rctd <- rctd
    spatial@misc[[paste0("rctd_", annotation_level)]] <- rctd
  }
  spatial
}

#' Extract doublet-mode style meta columns from an RCTD object
#' @noRd
.rctd_results_meta <- function(rctd, doublet_mode) {
  if (doublet_mode == "doublet" && !is.null(rctd@results$results_df)) {
    df <- as.data.frame(rctd@results$results_df)
    out <- data.frame(
      rctd_spot_class = as.character(df$spot_class),
      rctd_first_type = as.character(df$first_type),
      rctd_second_type = as.character(df$second_type),
      row.names = rownames(df),
      stringsAsFactors = FALSE
    )
    return(out)
  }

  ## full / multi: best-effort first type = max weight
  w <- .rctd_weight_matrix(rctd, doublet_mode)
  if (is.null(w) || nrow(w) == 0) {
    return(data.frame(
      rctd_spot_class = character(),
      rctd_first_type = character(),
      rctd_second_type = character()
    ))
  }
  first <- colnames(w)[max.col(as.matrix(w), ties.method = "first")]
  data.frame(
    rctd_spot_class = if (doublet_mode == "full") "full" else "multi",
    rctd_first_type = first,
    rctd_second_type = NA_character_,
    row.names = rownames(w),
    stringsAsFactors = FALSE
  )
}

#' Spot x cell-type weight matrix (rows = spots)
#' @noRd
.rctd_weight_matrix <- function(rctd, doublet_mode) {
  if (doublet_mode == "full" && !is.null(rctd@results$weights)) {
    w <- as.matrix(rctd@results$weights)
    rs <- rowSums(w)
    rs[rs == 0] <- 1
    return(w / rs)
  }

  if (doublet_mode == "doublet" &&
      !is.null(rctd@results$results_df) &&
      !is.null(rctd@results$weights_doublet)) {
    df <- rctd@results$results_df
    wd <- as.matrix(rctd@results$weights_doublet)
    ctypes <- rctd@cell_type_info$info[[2]]
    if (is.null(ctypes)) ctypes <- levels(df$first_type)
    w <- matrix(0, nrow = nrow(df), ncol = length(ctypes))
    rownames(w) <- rownames(df)
    colnames(w) <- ctypes
    for (i in seq_len(nrow(df))) {
      sc <- as.character(df$spot_class[i])
      if (sc == "reject") next
      t1 <- as.character(df$first_type[i])
      if (!is.na(t1) && t1 %in% colnames(w)) {
        w[i, t1] <- wd[i, 1]
      }
      if (sc %in% c("doublet_certain", "doublet_uncertain")) {
        t2 <- as.character(df$second_type[i])
        if (!is.na(t2) && t2 %in% colnames(w)) {
          w[i, t2] <- wd[i, 2]
        }
      }
    }
    rs <- rowSums(w)
    rs[rs == 0] <- 1
    return(w / rs)
  }

  if (doublet_mode == "multi" && is.list(rctd@results) && length(rctd@results) > 0) {
    ## multi mode: list per barcode with sub_weights / cell_type_list
    res <- rctd@results
    if (!is.null(res$results_df) || !is.null(res$weights)) {
      ## unexpected shape; fall through
    } else {
      barcodes <- names(res)
      ctypes <- sort(unique(unlist(lapply(res, function(x) x$cell_type_list), use.names = FALSE)))
      if (length(barcodes) > 0 && length(ctypes) > 0) {
        w <- matrix(0, nrow = length(barcodes), ncol = length(ctypes))
        rownames(w) <- barcodes
        colnames(w) <- ctypes
        for (bc in barcodes) {
          ct <- res[[bc]]$cell_type_list
          sw <- as.numeric(res[[bc]]$sub_weights)
          if (length(ct) && length(sw) == length(ct)) {
            w[bc, ct] <- sw
          }
        }
        rs <- rowSums(w)
        rs[rs == 0] <- 1
        return(w / rs)
      }
    }
  }

  NULL
}
