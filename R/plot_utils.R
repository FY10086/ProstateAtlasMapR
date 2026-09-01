#' Require plotting packages
#' @noRd
.check_ggplot <- function() {
  missing <- character()
  for (pkg in c("ggplot2", "patchwork", "ggsci", "scales")) {
    if (!requireNamespace(pkg, quietly = TRUE)) missing <- c(missing, pkg)
  }
  if (length(missing)) {
    cli::cli_abort(
      "Install required plotting package(s): {.pkg {missing}}"
    )
  }
  invisible(TRUE)
}

#' Annotation / score column names for RM results
#' @noRd
.rm_label_col <- function(annotation_level) {
  paste0("predicted.", .annotation_col(annotation_level))
}

#' @noRd
.rm_score_col <- function(annotation_level) {
  paste0(.rm_label_col(annotation_level), ".score")
}

#' Cell-type palette: D3 if <=20, else IGV ramp; Unknown/Other reserved greys
#' @noRd
.celltype_cols <- function(types) {
  types <- as.character(types)
  types <- types[!is.na(types) & nzchar(types)]
  types_main <- setdiff(unique(types), c("Unknown", "Other"))
  n <- length(types_main)
  if (n == 0) {
    return(c(Unknown = "grey72", Other = "grey58"))
  }
  if (n <= 20) {
    cols <- stats::setNames(ggsci::pal_d3("category20")(n), types_main)
  } else {
    base <- ggsci::pal_igv("default")(51)
    cols <- stats::setNames(grDevices::colorRampPalette(base)(n), types_main)
  }
  c(cols, Unknown = "grey72", Other = "grey58")
}

#' Shorten long cell-type labels for legends (aligned with 10-spatial.ipynb)
#'
#' 1) If >5 words: keep first + "..." + last 3 words.
#' 2) If still longer than `max_chars` (typical H3 marker+description):
#'    fall back to the marker prefix before the first space.
#' @noRd
.short_labels <- function(types, max_chars = 28) {
  vapply(types, function(x) {
    w <- strsplit(x, " ", fixed = TRUE)[[1]]
    lab <- if (length(w) <= 5) {
      x
    } else {
      paste(c(w[1], "...", utils::tail(w, 3)), collapse = " ")
    }
    if (nchar(lab) > max_chars) {
      prefix <- sub(" .*$", "", x)
      lab <- if (nzchar(prefix) && nchar(prefix) < nchar(lab)) prefix else lab
    }
    if (nchar(lab) > max_chars) {
      lab <- paste0(substr(lab, 1L, max_chars - 1L), "\u2026")
    }
    lab
  }, character(1), USE.NAMES = TRUE)
}

#' Extract RCTD weights + coords from a Seurat object after run_rctd()
#' @noRd
.rctd_extract <- function(spatial, reduction = "physical", level = NULL) {
  if (!inherits(spatial, "Seurat")) {
    cli::cli_abort("{.arg spatial} must be a Seurat object.")
  }
  if (is.null(level)) {
    level <- spatial@misc$rctd_annotation_level
    if (is.null(level)) level <- "H2"
  } else {
    level <- .match_level(level)
  }

  stored <- spatial@misc$rctd_annotation_level
  W_raw <- spatial@misc[[paste0("rctd_weights_", level)]]

  ## unprefixed matrix from the matching / only run
  if (is.null(W_raw) && !is.null(spatial@misc$rctd_weights)) {
    pref_keys <- grep("^rctd_weights_", names(spatial@misc), value = TRUE)
    if (identical(stored, level) || is.null(stored) || length(pref_keys) == 0L) {
      W_raw <- spatial@misc$rctd_weights
    }
  }

  if (is.null(W_raw)) {
    assay_nm <- paste0("RCTD_", level)
    if (assay_nm %in% Seurat::Assays(spatial)) {
      W_raw <- t(as.matrix(Seurat::GetAssayData(spatial, assay = assay_nm, layer = "counts")))
    } else if ("RCTD" %in% Seurat::Assays(spatial) &&
               (identical(stored, level) || is.null(stored))) {
      W_raw <- t(as.matrix(Seurat::GetAssayData(spatial, assay = "RCTD", layer = "counts")))
    }
  }

  if (is.null(W_raw)) {
    avail <- grep("^rctd_weights", names(spatial@misc), value = TRUE)
    assays <- grep("^RCTD", Seurat::Assays(spatial), value = TRUE)
    cli::cli_abort(c(
      "No RCTD weights found for annotation level {.val {level}}.",
      "i" = "Expected {.code @misc$rctd_weights_{level}} or unprefixed {.code rctd_weights} from that run.",
      if (!is.null(stored)) {
        c("i" = "Last {.fn run_rctd} level on this object: {.val {stored}}.")
      } else {
        NULL
      },
      if (length(avail)) {
        c("i" = "Weight keys in {.code @misc}: {.val {avail}}.")
      } else {
        NULL
      },
      if (length(assays)) {
        c("i" = "RCTD assays: {.val {assays}}.")
      } else {
        NULL
      },
      "i" = "Re-run {.code run_rctd(..., annotation_level = \"{level}\")} or pass the matching level."
    ))
  }

  if (!reduction %in% Seurat::Reductions(spatial)) {
    cli::cli_abort(c(
      "Reduction {.val {reduction}} not found.",
      "i" = "Available: {.val {Seurat::Reductions(spatial)}}"
    ))
  }

  W <- W_raw[rowSums(W_raw, na.rm = TRUE) > 0, , drop = FALSE]
  emb <- as.data.frame(Seurat::Embeddings(spatial, reduction = reduction))[, 1:2, drop = FALSE]
  colnames(emb) <- c("x", "y")
  emb$barcode <- rownames(emb)

  md <- spatial[[]]
  col_first <- paste0("rctd_", level, "_first_type")
  col_class <- paste0("rctd_", level, "_spot_class")
  if (!col_first %in% colnames(md)) col_first <- "rctd_first_type"
  if (!col_class %in% colnames(md)) col_class <- "rctd_spot_class"

  emb$first_type <- if (col_first %in% colnames(md)) md[emb$barcode, col_first] else NA_character_
  emb$spot_class <- if (col_class %in% colnames(md)) md[emb$barcode, col_class] else NA_character_
  emb$spot_class[is.na(emb$spot_class)] <- "filtered"

  Wa <- matrix(
    NA_real_,
    nrow = nrow(emb),
    ncol = ncol(W_raw),
    dimnames = list(emb$barcode, colnames(W_raw))
  )
  common <- intersect(rownames(W), emb$barcode)
  Wa[common, ] <- W[common, , drop = FALSE]

  mean_w <- sort(colMeans(W), decreasing = TRUE)
  n_first <- table(
    factor(
      emb$first_type[!emb$spot_class %in% c("reject", "filtered") & !is.na(emb$first_type)],
      levels = names(mean_w)
    )
  )

  list(
    W = W,
    Wa = Wa,
    xy = emb,
    mean_w = mean_w,
    n_first = n_first,
    cols_ct = .celltype_cols(names(mean_w)),
    lab_ct = .short_labels(colnames(W_raw)),
    level = level,
    mode = spatial@misc$rctd_doublet_mode %||% "unknown"
  )
}

#' Jupyter / IRkernel figure size for RCTD discrete maps (align with RM)
#'
#' H3 spatial maps set a provisional size; [.cowplot_side_legend()] then
#' retunes `repr.plot.width` / height to fit a single-column full-name legend.
#' @noRd
.set_rctd_repr_size <- function(annotation_level, kind = c("spatial", "pie")) {
  kind <- match.arg(kind)
  if (identical(annotation_level, "H3")) {
    ## provisional; .cowplot_side_legend() expands for long single-column legends
    ht <- if (identical(kind, "pie")) 14 else 13
    options(repr.plot.width = 22, repr.plot.height = ht)
  } else {
    if (identical(kind, "pie")) {
      options(repr.plot.width = 13, repr.plot.height = 11)
    } else {
      options(repr.plot.width = 12, repr.plot.height = 10)
    }
  }
  invisible(TRUE)
}

#' Wrap long legend labels (kept for callers that want it); H3 side legends
#' prefer unwrapped names and a wider canvas instead.
#' @noRd
.wrap_legend_labels <- function(labels, width = 32L) {
  vapply(as.character(labels), function(s) {
    if (!nzchar(s) || nchar(s) <= width) {
      return(s)
    }
    paste(strwrap(s, width = width), collapse = "\n")
  }, character(1), USE.NAMES = FALSE)
}

#' Cowplot side legend for long discrete *fill* labels (pie charts)
#' @noRd
.cowplot_fill_legend <- function(p_main, fill_vals, legend_labels, title = NULL) {
  .cowplot_side_legend(
    p_main = p_main,
    values = fill_vals,
    legend_labels = legend_labels,
    title = title,
    aesthetic = "fill"
  )
}

#' Cowplot side legend for long discrete colour labels (same idea as RM H3)
#' @noRd
.cowplot_map_legend <- function(p_main, cols, legend_labels = NULL, title = NULL) {
  .cowplot_side_legend(
    p_main = p_main,
    values = cols,
    legend_labels = legend_labels,
    title = title,
    aesthetic = "colour"
  )
}

#' Shared map | legend layout; single-column full names, large canvas
#'
#' Always one legend column (no 2-col split). Uses `theme(aspect.ratio)` so the
#' tissue fills the map column; sizes `repr.plot.*` to map + measured legend.
#' @noRd
.cowplot_side_legend <- function(p_main,
                                 values,
                                 legend_labels = NULL,
                                 title = NULL,
                                 aesthetic = c("colour", "fill")) {
  if (!requireNamespace("cowplot", quietly = TRUE)) {
    cli::cli_abort(
      "Package {.pkg cowplot} is required for H3 layouts. Install with {.code install.packages(\"cowplot\")}."
    )
  }
  aesthetic <- match.arg(aesthetic)
  lv <- names(values)
  if (is.null(legend_labels)) legend_labels <- lv
  ## Full names on one line — wider figure instead of wrapping / 2-col split
  legend_labels <- as.character(legend_labels)
  n_item <- length(lv)
  legend_ncol <- 1L
  legend_text_size <- if (n_item <= 18L) {
    12
  } else if (n_item <= 28L) {
    10.5
  } else if (n_item <= 40L) {
    9.5
  } else {
    8.5
  }
  legend_pt_size <- if (n_item <= 28L) 3.8 else 3.2

  ## Swap coord_fixed → aspect.ratio so the tissue fills the map column.
  asp <- .spatial_aspect_ratio(p_main)
  p_main <- .drop_coord_fixed(p_main) +
    ggplot2::theme(
      legend.position = "none",
      aspect.ratio = asp,
      plot.margin = ggplot2::margin(6, 6, 6, 6)
    )
  if (!is.null(title)) {
    p_main <- p_main + ggplot2::ggtitle(title)
  }

  d_leg <- data.frame(grp = factor(lv, levels = lv), stringsAsFactors = FALSE)
  if (identical(aesthetic, "fill")) {
    p_leg_src <- ggplot2::ggplot(d_leg, ggplot2::aes(x = 1, y = 1, fill = grp)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_manual(
        values = values,
        labels = legend_labels,
        drop = FALSE,
        name = NULL
      ) +
      ggplot2::guides(fill = ggplot2::guide_legend(ncol = legend_ncol))
  } else {
    p_leg_src <- ggplot2::ggplot(d_leg, ggplot2::aes(x = 1, y = 1, colour = grp)) +
      ggplot2::geom_point(size = legend_pt_size) +
      ggplot2::scale_colour_manual(
        values = values,
        labels = legend_labels,
        drop = FALSE,
        name = NULL
      ) +
      ggplot2::guides(
        colour = ggplot2::guide_legend(
          ncol = legend_ncol,
          override.aes = list(size = legend_pt_size)
        )
      )
  }
  p_leg_src <- p_leg_src +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = "right",
      legend.text = ggplot2::element_text(
        size = legend_text_size,
        colour = "black",
        lineheight = 1
      ),
      legend.key = ggplot2::element_rect(fill = NA, colour = NA),
      legend.key.spacing.y = grid::unit(2, "pt"),
      legend.background = ggplot2::element_rect(fill = "white", colour = NA),
      legend.margin = ggplot2::margin(0, 4, 0, 4),
      legend.box.margin = ggplot2::margin(0, 0, 0, 0),
      plot.margin = ggplot2::margin(0, 0, 0, 0)
    )

  leg <- cowplot::get_legend(p_leg_src)
  leg_cm <- tryCatch(
    grid::convertWidth(grid::grobWidth(leg), "cm", valueOnly = TRUE),
    error = function(e) 9
  )
  if (!is.finite(leg_cm) || leg_cm < 1) leg_cm <- 9
  ## Allow wide single-column H3 names (previously capped too tight → 2-col)
  leg_cm <- min(max(leg_cm + 0.4, 5), 16)
  leg_in <- leg_cm / 2.54

  ## Tall enough for one column of long labels; large map like RM H3
  fig_h <- max(13, min(24, 0.38 * n_item + 3.5))
  map_in <- max(11, fig_h * 0.85)
  options(
    repr.plot.width = map_in + leg_in + 0.25,
    repr.plot.height = fig_h
  )

  cowplot::plot_grid(
    p_main,
    leg,
    nrow = 1,
    rel_widths = c(map_in, leg_in),
    align = "h",
    axis = "tb"
  )
}

#' Panel height/width from a spatial ggplot's x/y ranges
#' @noRd
.spatial_aspect_ratio <- function(p) {
  grab_xy <- function(df) {
    if (is.null(df) || !is.data.frame(df)) {
      return(NULL)
    }
    xn <- intersect(c("x", "X", "imagecol"), names(df))
    yn <- intersect(c("y", "Y", "imagerow"), names(df))
    if (!length(xn) || !length(yn)) {
      return(NULL)
    }
    xr <- range(df[[xn[[1]]]], na.rm = TRUE)
    yr <- range(df[[yn[[1]]]], na.rm = TRUE)
    c(dx = diff(xr), dy = diff(yr))
  }

  spans <- grab_xy(p$data)
  if (is.null(spans) && length(p$layers)) {
    for (ly in p$layers) {
      spans <- grab_xy(ly$data)
      if (!is.null(spans)) break
    }
  }
  if (is.null(spans)) {
    built <- tryCatch(ggplot2::ggplot_build(p), error = function(e) NULL)
    if (!is.null(built) && length(built$layout$panel_params)) {
      pp <- built$layout$panel_params[[1]]
      xr <- pp$x.range %||% pp$x$range
      yr <- pp$y.range %||% pp$y$range
      if (!is.null(xr) && !is.null(yr)) {
        spans <- c(dx = diff(xr), dy = diff(yr))
      }
    }
  }
  if (is.null(spans)) {
    return(1)
  }
  dx <- spans[["dx"]]
  dy <- spans[["dy"]]
  if (!is.finite(dx) || !is.finite(dy) || dx <= 0 || dy <= 0) {
    return(1)
  }
  dy / dx
}

#' Remove CoordFixed so theme(aspect.ratio) can fill the map column
#' @noRd
.drop_coord_fixed <- function(p) {
  if (!inherits(p, "ggplot")) {
    return(p)
  }
  if (inherits(p$coordinates, "CoordFixed")) {
    ## Keep any existing limits from scales; plain Cartesian is enough.
    p$coordinates <- ggplot2::coord_cartesian(default = TRUE)
  }
  p
}

#' Default prostate marker genes for weight-vs-marker plots (H2-oriented)
#' @noRd
.default_rctd_markers <- function() {
  c(
    "Luminal cell" = "KLK3",
    "Basal cell" = "KRT14",
    "Club cell" = "PIGR",
    "Hillock cell" = "KRT13",
    "Neuroendocrine cell" = "ASCL1",
    "Endothelial cell" = "VWF",
    "Fibroblast" = "DCN",
    "Smooth muscle cell" = "ACTA2",
    "Pericyte" = "RGS5",
    "CD4 T cell" = "IL7R",
    "CD8 T cell" = "CD8A",
    "Natural killer cell" = "NKG7",
    "B cell" = "MS4A1",
    "Macrophage" = "C1QA",
    "Monocyte" = "FCN1",
    "Dendritic cell" = "FCER1A",
    "Mast cell" = "TPSAB1",
    "Schwann cell" = "S100B",
    "Gamma delta T cell" = "TRDC"
  )
}

#' Shared spatial theme (align with 10-spatial.ipynb `th_sp`)
#' @noRd
.theme_spatial <- function() {
  ggplot2::theme_void(base_size = 13) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, colour = "grey35", size = ggplot2::rel(0.8)),
      legend.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(0.85)),
      legend.position = "right",
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.margin = ggplot2::margin(6, 8, 6, 6)
    )
}

#' Auto point size from n spots
#' @noRd
.auto_pt_size <- function(n) {
  max(0.5, min(3.2, 150 / sqrt(max(n, 1))))
}

#' Median nearest-neighbour spacing (pie radius / Visium pitch)
#' @noRd
.spot_pitch <- function(xy) {
  xy <- as.data.frame(xy)[, 1:2, drop = FALSE]
  colnames(xy) <- c("x", "y")
  xy <- xy[stats::complete.cases(xy), , drop = FALSE]
  n <- nrow(xy)
  if (n < 2) {
    return(1)
  }
  if (requireNamespace("RANN", quietly = TRUE)) {
    d2 <- RANN::nn2(xy, k = 2)$nn.dists[, 2]
    p <- stats::median(d2, na.rm = TRUE)
    if (is.finite(p) && p > 0) {
      return(p)
    }
  }
  ## fallback without RANN: subsample pairwise minima
  idx <- if (n > 800) sample.int(n, 800) else seq_len(n)
  d <- as.matrix(stats::dist(xy[idx, , drop = FALSE]))
  diag(d) <- Inf
  p <- stats::median(apply(d, 1, min), na.rm = TRUE)
  if (!is.finite(p) || p <= 0) {
    dx <- diff(range(xy$x))
    dy <- diff(range(xy$y))
    p <- max(dx, dy) / sqrt(n) * 0.5
  }
  p
}

#' Finish spatial ggplot: optional y-flip + coord_fixed (true tissue aspect)
#'
#' Prefer `coord_fixed()` over `theme(aspect.ratio = 1)` so non-square Visium
#' sections are not stretched, matching 10-spatial.ipynb.
#' @noRd
.finish_spatial <- function(p, reverse_y = TRUE) {
  if (isTRUE(reverse_y)) {
    p <- p + ggplot2::scale_y_reverse()
  }
  p + ggplot2::coord_fixed()
}

#' Pick an expression assay (skip RCTD weight assays)
#' @noRd
.expr_assay <- function(spatial, assay = NULL) {
  assays <- Seurat::Assays(spatial)
  if (!is.null(assay) && assay %in% assays && !grepl("^RCTD", assay)) {
    return(assay)
  }
  for (a in c("RNA", "Spatial", assays)) {
    if (a %in% assays && !grepl("^RCTD", a)) {
      return(a)
    }
  }
  cli::cli_abort("No expression assay found (RNA/Spatial). Cannot plot markers.")
}

#' Ensure normalised data layer exists; return spatial + assay
#' @noRd
.ensure_normalized <- function(spatial, assay = NULL) {
  assay <- .expr_assay(spatial, assay)
  spatial <- .ensure_joined(spatial, assay = assay)
  layers <- tryCatch(
    SeuratObject::Layers(spatial[[assay]]),
    error = function(e) character()
  )
  has_data <- any(layers == "data" | grepl("^data(\\.|$)", layers))
  if (!has_data) {
    has_counts <- any(layers == "counts" | grepl("^counts(\\.|$)", layers))
    if (!has_counts) {
      cli::cli_abort(c(
        "Assay {.val {assay}} has no {.code counts}/{.code data} layer to normalise.",
        "i" = "Available layers: {.val {layers}}"
      ))
    }
    cli::cli_inform("Normalising assay {.val {assay}} (no {.code data} layer found).")
    spatial <- Seurat::NormalizeData(spatial, assay = assay, verbose = FALSE)
  }
  list(spatial = spatial, assay = assay)
}

#' Fetch one gene's log-normalised expression aligned to barcodes
#' @noRd
.fetch_expr <- function(spatial, gene, barcodes, assay) {
  mat <- tryCatch(
    SeuratObject::LayerData(spatial, assay = assay, layer = "data"),
    error = function(e) NULL
  )
  if (is.null(mat) || !gene %in% rownames(mat)) {
    return(rep(NA_real_, length(barcodes)))
  }
  as.numeric(mat[gene, barcodes, drop = TRUE])
}

#' Resolve marker gene for a cell type (H2 name / H3 prefix genes)
#' @noRd
.resolve_marker <- function(celltype, markers, available_genes = NULL) {
  if (is.null(celltype) || !nzchar(celltype) || is.na(celltype)) {
    return(NA_character_)
  }
  ## 1) exact key
  if (celltype %in% names(markers)) {
    return(unname(markers[[celltype]]))
  }
  ## 2) H3 style "GENE1+GENE2+ description" → first gene present in assay
  prefix <- sub(" .*$", "", celltype)
  genes <- strsplit(gsub("[-+]+$", "", prefix), "[-+]")[[1]]
  genes <- genes[nzchar(genes)]
  if (!is.null(available_genes) && length(genes)) {
    for (g in genes) {
      if (g %in% available_genes) {
        return(g)
      }
    }
  }
  ## 3) case-insensitive: H2 label is a substring of H3 name
  ct_l <- tolower(celltype)
  hit <- names(markers)[
    vapply(
      names(markers),
      function(nm) grepl(tolower(nm), ct_l, fixed = TRUE),
      logical(1)
    )
  ]
  if (length(hit)) {
    g <- unname(markers[[hit[[1]]]])
    if (is.null(available_genes) || g %in% available_genes) {
      return(g)
    }
  }
  if (length(genes)) {
    return(genes[[1]])
  }
  NA_character_
}

`%||%` <- function(x, y) if (is.null(x)) y else x
