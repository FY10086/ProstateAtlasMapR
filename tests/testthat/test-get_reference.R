test_that("get_reference() loads a local path directly without downloading", {
  ref_obj <- new_patlas_reference(
    seurat_obj = SeuratObject::CreateSeuratObject(counts = matrix(1:4, 2, 2, dimnames = list(c("g1", "g2"), c("c1", "c2")))),
    annotation_levels = "H1",
    reduction = "rpca",
    umap_reduction = "umap.rpca",
    npcs = 8,
    source_path = "dummy"
  )
  tmp <- withr::local_tempfile(fileext = ".qs")
  qs::qsave(ref_obj, tmp)

  loaded <- get_reference(path = tmp, verbose = FALSE)
  expect_s3_class(loaded, "patlas_reference")
  expect_equal(loaded$reduction, "rpca")
})

test_that("get_reference() downloads, verifies checksum, and caches", {
  skip_if_not_installed("openssl")

  ref_obj <- new_patlas_reference(
    seurat_obj = SeuratObject::CreateSeuratObject(counts = matrix(1:4, 2, 2, dimnames = list(c("g1", "g2"), c("c1", "c2")))),
    annotation_levels = "H1",
    reduction = "rpca",
    umap_reduction = "umap.rpca",
    npcs = 8,
    source_path = "dummy"
  )
  src_dir <- withr::local_tempdir()
  src_file <- file.path(src_dir, "ref_v1.qs")
  qs::qsave(ref_obj, src_file)
  url <- paste0("file://", src_file)
  con <- file(src_file, "rb")
  sha256 <- as.character(openssl::sha256(con))
  close(con)

  cache_dir <- withr::local_tempdir()

  loaded <- get_reference(url = url, sha256 = sha256, cache_dir = cache_dir, verbose = FALSE)
  expect_s3_class(loaded, "patlas_reference")
  cached_file <- file.path(cache_dir, "ref_v1.qs")
  expect_true(file.exists(cached_file))

  mtime_before <- file.info(cached_file)$mtime
  Sys.sleep(1)
  loaded2 <- get_reference(url = url, sha256 = sha256, cache_dir = cache_dir, verbose = FALSE)
  mtime_after <- file.info(cached_file)$mtime
  expect_equal(mtime_before, mtime_after) # should reuse cache, not re-download
})

test_that("get_reference() aborts on checksum mismatch", {
  skip_if_not_installed("openssl")

  ref_obj <- new_patlas_reference(
    seurat_obj = SeuratObject::CreateSeuratObject(counts = matrix(1:4, 2, 2, dimnames = list(c("g1", "g2"), c("c1", "c2")))),
    annotation_levels = "H1",
    reduction = "rpca",
    umap_reduction = "umap.rpca",
    npcs = 8,
    source_path = "dummy"
  )
  src_dir <- withr::local_tempdir()
  src_file <- file.path(src_dir, "ref_bad.qs")
  qs::qsave(ref_obj, src_file)
  url <- paste0("file://", src_file)

  cache_dir <- withr::local_tempdir()
  expect_error(
    get_reference(url = url, sha256 = strrep("0", 64), cache_dir = cache_dir, verbose = FALSE),
    "checksum"
  )
})

test_that("get_reference() falls back to baked-in Release URL when options are NULL", {
  withr::local_options(list(
    ProstateAtlasMapR.reference_url = NULL,
    ProstateAtlasMapR.reference_sha256 = NULL
  ))
  ## Should not abort on missing URL; resolves to GitHub Release default.
  ## Do not download: intercept by passing a tiny local file via path instead
  ## and only check that url resolution would succeed through the helpers.
  expect_identical(
    ProstateAtlasMapR:::.default_reference_url(),
    "https://github.com/FY10086/ProstateAtlasMapR/releases/download/v0.1.0/prostate_atlas_reference.qs"
  )
  expect_true(nzchar(ProstateAtlasMapR:::.default_reference_sha256()))
})
