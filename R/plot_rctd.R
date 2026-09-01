#' Spatial map of RCTD dominant cell type (first_type)
#'
#' Reject / missing calls are shown as `"Unknown"`. Rare types below
#' `min_pct` of spots can be merged into `"Other"`.
#'
#' **H1 / H2:** original compact map with legend on the right
#' (`repr.plot` 12 x 10).
#' **H3:** map and legend split with measured legend width (`repr.plot`
#' height-matched) so long subtype names do not crush the tissue or leave a
#' large gap beside `coord_fixed` spatial panels.
#'
#' @param spatial Seurat object after [run_rctd()].
#' @param annotation_level Atlas level `"H1"`, `"H2"`, or `"H3"`. Default
#'   `NULL` uses the level stored by the last [run_rctd()] call (else `"H2"`).
#' @param reduction Spatial coordinates reduction. Default `"physical"`.
#' @param min_pct Merge types below this percent of spots into Other.
#'   Default `0` (no merge).
#' @param reverse_y Flip the y axis (common for Visium). Default `TRUE`.
#'
#' @return A ggplot object (H1/H2) or cowplot grob (H3).
#' @export
plot_rctd_first_type <- function(spatial,
                                 annotation_level = NULL,
                                 reduction = "physical",
                                 min_pct = 0,
                                 reverse_y = TRUE) {
  .check_ggplot()
  dat <- .rctd_extract(spatial, reduction = reduction, level = annotation_level)
  xy <- dat$xy
  level <- dat$level
  .set_rctd_repr_size(level)

  ft <- as.character(xy$first_type)
  ft[xy$spot_class %in% c("reject", "filtered") | is.na(ft)] <- "Unknown"
  tab <- sort(table(ft[ft != "Unknown"]), decreasing = TRUE)
  big <- names(tab)[as.numeric(tab) / length(ft) * 100 >= min_pct]
  ft[!ft %in% c(big, "Unknown")] <- "Other"
  lv <- c(big, intersect(c("Other", "Unknown"), unique(ft)))

  d <- xy
  d$grp <- factor(ft, levels = lv)
  cols <- dat$cols_ct
  cols <- cols[intersect(names(cols), lv)]
  if ("Other" %in% lv) cols["Other"] <- "grey58"
  if ("Unknown" %in% lv) cols["Unknown"] <- "grey72"

  counts <- as.numeric(table(d$grp))
  names(counts) <- levels(d$grp)
  title <- paste0("celltype - ", level)

  ## H3: cowplot map | legend (full names), same idea as plot_rm_prediction
  if (identical(level, "H3")) {
    p_main <- ggplot2::ggplot(d, ggplot2::aes(x = x, y = y, colour = grp)) +
      ggplot2::geom_point(size = .auto_pt_size(nrow(d)), shape = 16) +
      ggplot2::scale_colour_manual(values = cols, name = NULL, drop = FALSE) +
      ggplot2::labs(title = title) +
      .theme_spatial()
    p_main <- .finish_spatial(p_main, reverse_y)
    return(.cowplot_map_legend(p_main, cols, legend_labels = lv, title = title))
  }

  ## H1 / H2: original right-side legend with shortened labels + counts
  lab <- c(dat$lab_ct, Other = "Other", Unknown = "Unknown")
  legend_lab <- sprintf("%s  (%s)", lab[lv], scales::comma(counts[lv]))
  p <- ggplot2::ggplot(d, ggplot2::aes(x = x, y = y, colour = grp)) +
    ggplot2::geom_point(size = .auto_pt_size(nrow(d)), shape = 16) +
    ggplot2::scale_colour_manual(
      values = cols,
      name = NULL,
      labels = legend_lab,
      drop = FALSE
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        ncol = max(1L, ceiling(length(lv) / 20)),
        override.aes = list(size = 3.2)
      )
    ) +
    ggplot2::labs(title = title) +
    .theme_spatial() +
    ggplot2::theme(
      legend.position = "right",
      legend.text = ggplot2::element_text(size = 9)
    )

  .finish_spatial(p, reverse_y)
}

#' RCTD composition summary: mean weight + dominant-type counts
#'
#' Side-by-side lollipop (mean weight %) and bar (n spots as first_type).
#'
#' @param spatial Seurat object after [run_rctd()].
#' @param annotation_level Atlas level `"H1"`, `"H2"`, or `"H3"`. Default
#'   `NULL` uses the level stored by the last [run_rctd()] call (else `"H2"`).
#' @param n_types Max cell types to show (by mean weight). Default `25`.
#'
#' @return A patchwork object.
#' @export
plot_rctd_composition <- function(spatial, annotation_level = NULL, n_types = 25) {
  .check_ggplot()
  dat <- .rctd_extract(spatial, level = annotation_level)
  ## composition is wide; keep a readable notebook size
  if (identical(dat$level, "H3")) {
    options(repr.plot.width = 16, repr.plot.height = 8)
  } else {
    options(repr.plot.width = 13, repr.plot.height = 7)
  }
  mean_w <- dat$mean_w
  kept <- utils::head(names(mean_w)[mean_w > 1e-4], n_types)
  if (!length(kept)) {
    cli::cli_abort("No cell types with mean weight > 1e-4 to plot.")
  }
  d <- data.frame(
    celltype = factor(kept, levels = rev(kept)),
    mean_w = as.numeric(mean_w[kept]) * 100,
    n_first = as.numeric(dat$n_first[kept]),
    stringsAsFactors = FALSE
  )

  th <- ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = 12),
      plot.margin = ggplot2::margin(6, 10, 6, 6)
    )

  p_lolli <- ggplot2::ggplot(d, ggplot2::aes(x = mean_w, y = celltype, colour = celltype)) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, xend = mean_w, yend = celltype),
      linewidth = 0.9
    ) +
    ggplot2::geom_point(size = 3.2) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.1f", mean_w)),
      hjust = 0,
      nudge_x = max(d$mean_w) * 0.02,
      size = 3.1,
      colour = "grey25"
    ) +
    ggplot2::scale_colour_manual(values = dat$cols_ct, guide = "none") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.13))) +
    ggplot2::scale_y_discrete(labels = dat$lab_ct) +
    ggplot2::labs(x = "Mean weight across spots (%)", y = NULL, title = "Cell-type composition") +
    th +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())

  p_first <- ggplot2::ggplot(d, ggplot2::aes(x = n_first, y = celltype, fill = celltype)) +
    ggplot2::geom_col(width = 0.68) +
    ggplot2::geom_text(
      ggplot2::aes(label = scales::comma(n_first)),
      hjust = 0,
      nudge_x = max(d$n_first, na.rm = TRUE) * 0.012,
      size = 3.1,
      colour = "grey25"
    ) +
    ggplot2::scale_fill_manual(values = dat$cols_ct, guide = "none") +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.15)),
      labels = scales::comma
    ) +
    ggplot2::labs(x = "Dominant spots (first_type)", y = NULL, title = "Dominant cell type per spot") +
    th +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank()
    )

  ## left panel needs more width for long cell-type names (notebook: p_lolli | p_first)
  p_lolli + p_first + patchwork::plot_layout(widths = c(1.35, 1))
}

#' Faceted spatial maps of top RCTD cell-type weights
#'
#' @param spatial Seurat object after [run_rctd()].
#' @param annotation_level Atlas level `"H1"`, `"H2"`, or `"H3"`. Default
#'   `NULL` uses the level stored by the last [run_rctd()] call (else `"H2"`).
#' @param reduction Spatial coordinates reduction. Default `"physical"`.
#' @param n_types Number of top types (by mean weight). Default `11`.
#' @param ncol Facet columns. Default `3`.
#' @param reverse_y Flip the y axis. Default `TRUE`.
#'
#' @return A ggplot object.
#' @export
plot_rctd_weights <- function(spatial,
                              annotation_level = NULL,
                              reduction = "physical",
                              n_types = 11,
                              ncol = 3,
                              reverse_y = TRUE) {
  .check_ggplot()
  dat <- .rctd_extract(spatial, reduction = reduction, level = annotation_level)
  cts <- utils::head(names(dat$mean_w)[dat$mean_w > 1e-4], n_types)
  if (!length(cts)) {
    cli::cli_abort("No cell types with mean weight > 1e-4 to plot.")
  }
  ## notebook: width=13, height=4.3 * n_rows + 0.8
  options(
    repr.plot.width = 13,
    repr.plot.height = 4.3 * ceiling(length(cts) / ncol) + 0.8
  )
  xy <- dat$xy
  Wa <- dat$Wa

  d <- data.frame(
    x = rep(xy$x, length(cts)),
    y = rep(xy$y, length(cts)),
    ct = factor(rep(cts, each = nrow(xy)), levels = cts),
    val = as.numeric(Wa[, cts, drop = FALSE]),
    stringsAsFactors = FALSE
  )
  d <- d[order(d$val, na.last = FALSE), , drop = FALSE]
  w_hi <- as.numeric(stats::quantile(d$val, 0.99, na.rm = TRUE))
  if (!is.finite(w_hi) || w_hi <= 0) w_hi <- 1

  ## H3 facets: keep readable shortened strip labels
  facet_lab <- stats::setNames(
    sprintf("%s  (mean %.1f%%)", dat$lab_ct[cts], dat$mean_w[cts] * 100),
    cts
  )
  pal_weight <- c("grey93", "#FFF3B0", "#FDC966", "#F98E52", "#E03B3B", "#8C0F1A")

  p <- ggplot2::ggplot(d, ggplot2::aes(x = x, y = y, colour = val)) +
    ggplot2::geom_point(size = .auto_pt_size(nrow(xy)) * 0.8, shape = 16) +
    ggplot2::facet_wrap(
      ~ct,
      ncol = ncol,
      labeller = ggplot2::labeller(ct = facet_lab)
    ) +
    ggplot2::scale_colour_gradientn(
      colours = pal_weight,
      limits = c(0, w_hi),
      oob = scales::squish,
      na.value = "grey92",
      name = "Weight",
      trans = "sqrt",
      breaks = c(0, 0.05, 0.15, 0.3, 0.6),
      guide = ggplot2::guide_colourbar(barheight = grid::unit(3.2, "cm"))
    ) +
    ggplot2::labs(title = paste0("RCTD weights - top ", length(cts), " types (", dat$level, ")")) +
    .theme_spatial() +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold", size = 10, margin = ggplot2::margin(2, 2, 6, 2)),
      panel.spacing = grid::unit(1.1, "lines")
    )

  .finish_spatial(p, reverse_y)
}

#' RCTD weight vs marker gene expression (side by side)
#'
#' @param spatial Seurat object after [run_rctd()] (must still have expression assay).
#' @param annotation_level Atlas level `"H1"`, `"H2"`, or `"H3"`. Default
#'   `NULL` uses the level stored by the last [run_rctd()] call (else `"H2"`).
#' @param reduction Spatial coordinates reduction. Default `"physical"`.
#' @param n_types Number of top types. Default `6`.
#' @param markers Named character vector: `celltype = gene`. `NULL` uses
#'   built-in prostate markers.
#' @param assay Expression assay. `NULL` = DefaultAssay.
#' @param reverse_y Flip the y axis. Default `TRUE`.
#'
#' @return A patchwork object.
#' @export
plot_rctd_weight_vs_marker <- function(spatial,
                                       annotation_level = NULL,
                                       reduction = "physical",
                                       n_types = 6,
                                       markers = NULL,
                                       assay = NULL,
                                       reverse_y = TRUE) {
  .check_ggplot()
  dat <- .rctd_extract(spatial, reduction = reduction, level = annotation_level)
  cts <- utils::head(names(dat$mean_w)[dat$mean_w > 1e-4], n_types)
  if (!length(cts)) {
    cli::cli_abort("No cell types with mean weight > 1e-4 to plot.")
  }
  if (is.null(markers)) markers <- .default_rctd_markers()

  ## Prefer RNA/Spatial; auto-NormalizeData if no "data" layer (notebook does the same)
  normed <- .ensure_normalized(spatial, assay = assay)
  spatial <- normed$spatial
  assay <- normed$assay

  ## notebook: width=13, height=6.2 * n_types + 0.6 (two spatial panels)
  options(
    repr.plot.width = 13,
    repr.plot.height = 6.2 * length(cts) + 0.6
  )

  .marker_placeholder <- function(msg = "No marker gene\navailable in assay") {
    ggplot2::ggplot() +
      ggplot2::annotate(
        "text",
        x = 0.5,
        y = 0.5,
        label = msg,
        size = 4,
        colour = "grey40"
      ) +
      ggplot2::theme_void() +
      ggplot2::labs(title = "Marker expression")
  }

  pal_weight <- c("grey93", "#FFF3B0", "#FDC966", "#F98E52", "#E03B3B", "#8C0F1A")
  pal_expr <- c("grey93", "#DEEBF7", "#9ECAE1", "#4292C6", "#2171B5", "#08306B")
  xy <- dat$xy
  Wa <- dat$Wa
  pt <- .auto_pt_size(nrow(xy)) * 0.8

  plist <- lapply(cts, function(celltype) {
    w_vec <- Wa[, celltype]
    d_one <- data.frame(xy[, c("x", "y")], val = w_vec)
    d_one <- d_one[order(d_one$val, na.last = FALSE), , drop = FALSE]
    w_lim <- as.numeric(stats::quantile(w_vec, 0.99, na.rm = TRUE))
    if (!is.finite(w_lim) || w_lim <= 0) w_lim <- 1

    ct_lab <- dat$lab_ct[[celltype]] %||% celltype
    p_w <- ggplot2::ggplot(d_one, ggplot2::aes(x = x, y = y, colour = val)) +
      ggplot2::geom_point(size = pt, shape = 16) +
      ggplot2::scale_colour_gradientn(
        colours = pal_weight,
        na.value = "grey92",
        name = "Weight",
        limits = c(0, w_lim),
        oob = scales::squish,
        guide = ggplot2::guide_colourbar(barheight = grid::unit(2.4, "cm"))
      ) +
      ggplot2::labs(
        title = paste0(ct_lab, " - RCTD weight"),
        subtitle = sprintf(
          "mean %.1f%% | first_type in %s spots",
          dat$mean_w[[celltype]] * 100,
          scales::comma(as.numeric(dat$n_first[[celltype]]))
        )
      ) +
      .theme_spatial()
    p_w <- .finish_spatial(p_w, reverse_y)

    marker <- .resolve_marker(
      celltype,
      markers,
      available_genes = rownames(spatial[[assay]])
    )
    if (is.na(marker) || !nzchar(marker) || !marker %in% rownames(spatial[[assay]])) {
      ## keep 2-column layout even without a marker
      return(p_w + .marker_placeholder() + patchwork::plot_layout(widths = c(1, 1)))
    }

    expr <- .fetch_expr(spatial, marker, xy$barcode, assay)
    if (all(is.na(expr))) {
      return(
        p_w +
          .marker_placeholder(paste0(marker, "\n(no expression values)")) +
          patchwork::plot_layout(widths = c(1, 1))
      )
    }
    r_sp <- suppressWarnings(stats::cor(w_vec, expr, method = "spearman", use = "complete.obs"))
    d_mk <- data.frame(xy[, c("x", "y")], expr = expr)
    d_mk <- d_mk[order(d_mk$expr), , drop = FALSE]
    e_lim <- as.numeric(stats::quantile(expr, c(0.05, 0.99), na.rm = TRUE))
    if (any(!is.finite(e_lim))) e_lim <- range(expr, na.rm = TRUE)

    p_mk <- ggplot2::ggplot(d_mk, ggplot2::aes(x = x, y = y, colour = expr)) +
      ggplot2::geom_point(size = pt, shape = 16) +
      ggplot2::scale_colour_gradientn(
        colours = pal_expr,
        na.value = "grey92",
        name = "Expr",
        limits = e_lim,
        oob = scales::squish,
        guide = ggplot2::guide_colourbar(barheight = grid::unit(2.4, "cm"))
      ) +
      ggplot2::labs(
        title = paste0(marker, " - log-normalised expression"),
        subtitle = sprintf("Spearman r with %s weight = %.2f", ct_lab, r_sp)
      ) +
      .theme_spatial()
    p_mk <- .finish_spatial(p_mk, reverse_y)

    p_w + p_mk + patchwork::plot_layout(widths = c(1, 1))
  })

  patchwork::wrap_plots(plist, ncol = 1) +
    patchwork::plot_annotation(
      title = paste0("RCTD weight vs marker - top ", length(cts), " types (", dat$level, ")"),
      theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 14))
    )
}

#' Per-spot pie charts of RCTD cell-type composition
#'
#' @param spatial Seurat object after [run_rctd()].
#' @param annotation_level Atlas level `"H1"`, `"H2"`, or `"H3"`. Default
#'   `NULL` uses the level stored by the last [run_rctd()] call (else `"H2"`).
#' @param reduction Spatial coordinates reduction. Default `"physical"`.
#' @param top_n Number of types to show; remaining mass merged as Other.
#'   Default `12`.
#' @param reverse_y Flip the y axis. Default `TRUE`.
#'
#' @return A ggplot object.
#' @export
plot_rctd_pie <- function(spatial,
                          annotation_level = NULL,
                          reduction = "physical",
                          top_n = 12,
                          reverse_y = TRUE) {
  .check_ggplot()
  if (!requireNamespace("scatterpie", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg scatterpie} is required for {.fn plot_rctd_pie}.")
  }
  dat <- .rctd_extract(spatial, reduction = reduction, level = annotation_level)
  .set_rctd_repr_size(dat$level, kind = "pie")
  top_types <- utils::head(names(dat$mean_w), top_n)
  W <- dat$W
  xy <- dat$xy

  mat <- W
  other_cols <- setdiff(colnames(mat), top_types)
  comp <- cbind(
    mat[, top_types, drop = FALSE],
    Other = if (length(other_cols)) rowSums(mat[, other_cols, drop = FALSE]) else 0
  )
  rs <- rowSums(comp)
  rs[rs == 0] <- 1
  comp <- comp / rs

  pie_cols <- colnames(comp)
  safe <- make.names(pie_cols, unique = TRUE)
  d_pie <- cbind(xy[rownames(W), c("x", "y"), drop = FALSE], stats::setNames(as.data.frame(comp), safe))
  ## notebook: pitch <- median(nn distance); r <- pitch * 0.48
  pitch <- .spot_pitch(d_pie[, c("x", "y")])
  d_pie$r <- pitch * 0.48

  fill_vals <- c(unname(dat$cols_ct[top_types]), "grey80")
  names(fill_vals) <- safe
  ## H3: full names; H1/H2: shortened
  if (identical(dat$level, "H3")) {
    fill_labs <- unname(c(top_types, "Other"))
  } else {
    fill_labs <- unname(c(dat$lab_ct[top_types], "Other"))
  }
  title <- paste0("Per-spot cell-type composition (", dat$level, ")")

  p <- ggplot2::ggplot() +
    scatterpie::geom_scatterpie(
      data = d_pie,
      ggplot2::aes(x = x, y = y, r = r),
      cols = safe,
      colour = NA
    ) +
    ggplot2::scale_fill_manual(values = fill_vals, labels = fill_labs, name = NULL) +
    ggplot2::labs(title = title) +
    .theme_spatial() +
    ggplot2::theme(legend.key.size = grid::unit(0.45, "cm"))

  p <- .finish_spatial(p, reverse_y)

  if (identical(dat$level, "H3")) {
    return(.cowplot_fill_legend(p, fill_vals, fill_labs, title = title))
  }
  p
}

#' RCTD QC maps: dominant weight and Shannon effective type number
#'
#' @param spatial Seurat object after [run_rctd()].
#' @param annotation_level Atlas level `"H1"`, `"H2"`, or `"H3"`. Default
#'   `NULL` uses the level stored by the last [run_rctd()] call (else `"H2"`).
#' @param reduction Spatial coordinates reduction. Default `"physical"`.
#' @param reverse_y Flip the y axis. Default `TRUE`.
#'
#' @return A patchwork object.
#' @export
plot_rctd_qc <- function(spatial,
                         annotation_level = NULL,
                         reduction = "physical",
                         reverse_y = TRUE) {
  .check_ggplot()
  dat <- .rctd_extract(spatial, reduction = reduction, level = annotation_level)
  ## notebook: 11 x 5.5
  options(repr.plot.width = 11, repr.plot.height = 5.5)
  xy <- dat$xy
  W <- dat$W

  dom_all <- stats::setNames(rep(NA_real_, nrow(xy)), xy$barcode)
  eff_all <- dom_all
  dom_all[rownames(W)] <- apply(W, 1, max)
  P <- W
  P[P <= 0] <- NA
  eff_all[rownames(W)] <- exp(-rowSums(P * log(P), na.rm = TRUE))

  d <- data.frame(
    x = xy$x,
    y = xy$y,
    dom = as.numeric(dom_all),
    eff = as.numeric(eff_all)
  )
  pal_weight <- c("grey93", "#FFF3B0", "#FDC966", "#F98E52", "#E03B3B", "#8C0F1A")
  pt <- .auto_pt_size(nrow(d))

  p_dom <- ggplot2::ggplot(
    d[order(d$dom, na.last = FALSE), , drop = FALSE],
    ggplot2::aes(x = x, y = y, colour = dom)
  ) +
    ggplot2::geom_point(size = pt, shape = 16) +
    ggplot2::scale_colour_gradientn(
      colours = pal_weight,
      na.value = "grey92",
      name = "Max\nweight",
      guide = ggplot2::guide_colourbar(barheight = grid::unit(2.6, "cm"))
    ) +
    ggplot2::labs(title = "Dominant cell-type fraction") +
    .theme_spatial()

  p_eff <- ggplot2::ggplot(
    d[order(d$eff, na.last = FALSE), , drop = FALSE],
    ggplot2::aes(x = x, y = y, colour = eff)
  ) +
    ggplot2::geom_point(size = pt, shape = 16) +
    ggplot2::scale_colour_viridis_c(
      option = "mako",
      direction = -1,
      na.value = "grey92",
      name = "Eff.\ntypes",
      guide = ggplot2::guide_colourbar(barheight = grid::unit(2.6, "cm"))
    ) +
    ggplot2::labs(title = "Shannon entropy") +
    .theme_spatial()

  p_dom <- .finish_spatial(p_dom, reverse_y)
  p_eff <- .finish_spatial(p_eff, reverse_y)

  ## notebook: p_dom | p_eff
  p_dom + p_eff + patchwork::plot_layout(widths = c(1, 1))
}
