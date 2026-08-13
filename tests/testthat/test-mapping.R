test_that("build_reference() attaches a UMAP model without moving existing coordinates", {
  skip_if_not_installed("uwot")

  atlas <- make_synthetic_atlas(seed = 1, npcs = 8)
  tmp_qs <- withr::local_tempfile(fileext = ".qs")
  qs::qsave(atlas, tmp_qs)

  orig_umap <- Seurat::Embeddings(atlas[["umap.rpca"]])

  ref <- suppressMessages(build_reference(
    qc_path = tmp_qs,
    npcs = 8,
    save_path = NULL,
    verbose = FALSE
  ))

  expect_s3_class(ref, "patlas_reference")
  expect_setequal(ref$annotation_levels, c("H1", "H2", "H3"))
  expect_true(ref$umap_reduction %in% Seurat::Reductions(ref$seurat_obj))

  new_umap <- Seurat::Embeddings(ref$seurat_obj[["umap.rpca"]])
  common <- intersect(rownames(orig_umap), rownames(new_umap))
  expect_true(length(common) > 0)
  expect_equal(
    unname(as.matrix(new_umap[common, , drop = FALSE])),
    unname(as.matrix(orig_umap[common, , drop = FALSE])),
    tolerance = 1e-8
  )

  ## the projection model must have been attached
  model <- Seurat::Misc(ref$seurat_obj[["umap.rpca"]], slot = "model")
  expect_true(length(model) > 0)

  ## metadata should only contain the retained annotation columns
  expect_setequal(colnames(ref$seurat_obj@meta.data), c("celltype_manual_H1", "celltype_manual_H2", "celltype_manual_H3"))
})

test_that("map_query() transfers labels and projects onto the reference UMAP", {
  skip_if_not_installed("uwot")

  atlas <- make_synthetic_atlas(seed = 2, npcs = 8)
  tmp_qs <- withr::local_tempfile(fileext = ".qs")
  qs::qsave(atlas, tmp_qs)

  ref <- suppressMessages(build_reference(qc_path = tmp_qs, npcs = 8, verbose = FALSE))
  before_umap <- Seurat::Embeddings(ref$seurat_obj[["umap.rpca"]])

  query <- make_synthetic_query(gene_names = rownames(ref$seurat_obj), n_cells = 40, seed = 10)

  for (lvl in c("H1", "H2")) {
    mapped <- suppressMessages(suppressWarnings(map_query(
      reference = ref,
      query = query,
      annotation_level = lvl,
      k.anchor = 3,
      k.score = 10,
      k.weight = 15,
      verbose = FALSE
    )))

    pred_col <- paste0("predicted.celltype_manual_", lvl)
    expect_true(pred_col %in% colnames(mapped@meta.data))
    expect_true(paste0(pred_col, ".score") %in% colnames(mapped@meta.data))
    expect_false(any(is.na(mapped@meta.data[[pred_col]])))

    expect_true("umap.rpca" %in% Seurat::Reductions(mapped))
    expect_equal(ncol(Seurat::Embeddings(mapped[["umap.rpca"]])), 2)
    expect_equal(nrow(Seurat::Embeddings(mapped[["umap.rpca"]])), ncol(query))
  }

  ## the reference's own umap.rpca coordinates must be untouched after mapping
  after_umap <- Seurat::Embeddings(ref$seurat_obj[["umap.rpca"]])
  expect_equal(unname(as.matrix(after_umap)), unname(as.matrix(before_umap)), tolerance = 1e-8)
})

test_that("map_query() drops NA-annotated reference cells and warns", {
  skip_if_not_installed("uwot")

  atlas <- make_synthetic_atlas(seed = 3, npcs = 8)
  tmp_qs <- withr::local_tempfile(fileext = ".qs")
  qs::qsave(atlas, tmp_qs)
  ref <- suppressMessages(build_reference(qc_path = tmp_qs, npcs = 8, verbose = FALSE))

  expect_true(anyNA(ref$seurat_obj$celltype_manual_H3))

  query <- make_synthetic_query(gene_names = rownames(ref$seurat_obj), n_cells = 20, seed = 11)
  msgs <- character()
  mapped <- withCallingHandlers(
    suppressWarnings(map_query(
      ref, query,
      annotation_level = "H3",
      k.anchor = 3, k.score = 10, k.weight = 10,
      verbose = TRUE
    )),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  expect_true(any(grepl("NA", msgs, fixed = TRUE)))
  expect_true("predicted.celltype_manual_H3" %in% colnames(mapped@meta.data))
})

test_that("map_query() errors on an annotation level absent from the reference", {
  skip_if_not_installed("uwot")

  atlas <- make_synthetic_atlas(seed = 4, npcs = 8)
  tmp_qs <- withr::local_tempfile(fileext = ".qs")
  qs::qsave(atlas, tmp_qs)
  ref <- suppressMessages(build_reference(qc_path = tmp_qs, npcs = 8, annotation_levels = "H1", verbose = FALSE))

  query <- make_synthetic_query(gene_names = rownames(ref$seurat_obj), n_cells = 10, seed = 12)
  expect_error(map_query(ref, query, annotation_level = "H2", verbose = FALSE), "not available")
})

test_that("print.patlas_reference() runs without error", {
  skip_if_not_installed("uwot")
  atlas <- make_synthetic_atlas(seed = 5, npcs = 8)
  tmp_qs <- withr::local_tempfile(fileext = ".qs")
  qs::qsave(atlas, tmp_qs)
  ref <- suppressMessages(build_reference(qc_path = tmp_qs, npcs = 8, verbose = FALSE))
  expect_message(print(ref), "patlas_reference")
})
