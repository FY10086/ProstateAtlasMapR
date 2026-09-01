skip_if_not_installed("spacexr")

test_that("run_rctd writes doublet labels and a weight assay", {
  ref_seu <- make_synthetic_rctd_reference()
  spat <- make_synthetic_spatial(rownames(ref_seu), n_spots = 24)

  out <- run_rctd(
    reference = ref_seu,
    spatial = spat$seurat,
    annotation_level = "H1",
    coords = "spatial",
    doublet_mode = "doublet",
    max_cores = 1,
    n_max_cells = 50,
    min_UMI = 1,
    CELL_MIN_INSTANCE = 10,
    UMI_min = 1,
    UMI_min_sigma = 1,
    counts_MIN = 1,
    gene_cutoff = 0,
    fc_cutoff = 0.1,
    gene_cutoff_reg = 0,
    fc_cutoff_reg = 0.1,
    assay_name = "RCTD",
    store_rctd = TRUE,
    verbose = FALSE
  )

  expect_true(inherits(out, "Seurat"))
  expect_true("rctd_spot_class" %in% colnames(out[[]]))
  expect_true("rctd_first_type" %in% colnames(out[[]]))
  expect_true("RCTD" %in% Seurat::Assays(out))
  expect_true(!is.null(out@misc$rctd))
  expect_true(!is.null(out@misc$rctd_weights))
  expect_true(any(!is.na(out$rctd_first_type)))
  expect_true(all(out$rctd_spot_class %in% c(
    "singlet", "doublet_certain", "doublet_uncertain", "reject"
  )))
})

test_that("run_rctd accepts patlas_reference and class_level", {
  atlas <- make_synthetic_atlas(n_genes = 180, n_cells_per_batch = 80, npcs = 8)
  ref <- build_reference(
    qc_path = {
      tmp <- tempfile(fileext = ".qs")
      qs::qsave(atlas, tmp)
      tmp
    },
    npcs = 8,
    n.neighbors = 15,
    verbose = FALSE
  )
  spat <- make_synthetic_spatial(rownames(ref$seurat_obj), n_spots = 20)

  out <- run_rctd(
    reference = ref,
    spatial = spat$seurat,
    annotation_level = "H2",
    class_level = "H1",
    coords = spat$coords,
    doublet_mode = "full",
    max_cores = 1,
    n_max_cells = 40,
    min_UMI = 1,
    CELL_MIN_INSTANCE = 5,
    UMI_min = 1,
    UMI_min_sigma = 1,
    counts_MIN = 1,
    gene_cutoff = 0,
    fc_cutoff = 0.1,
    gene_cutoff_reg = 0,
    fc_cutoff_reg = 0.1,
    store_rctd = FALSE,
    verbose = FALSE
  )

  expect_equal(out@misc$rctd_annotation_level, "H2")
  expect_equal(out@misc$rctd_doublet_mode, "full")
  expect_true("RCTD" %in% Seurat::Assays(out))
  expect_null(out@misc$rctd)
})
