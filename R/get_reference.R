#' Get the prebuilt prostate atlas reference (download + cache, or load local)
#'
#' The lean reference produced by [build_reference()] is distributed
#' separately from this package's source code (e.g. via Zenodo or GitHub
#' Releases), since it is still too large to ship inside an R package. This
#' function downloads it once, verifies its integrity, caches it locally,
#' and returns the loaded `patlas_reference` object -- or simply loads a
#' local copy if you already have one.
#'
#' @param path Path to an existing local reference file. If supplied, the
#'   file is read directly and no download is attempted.
#' @param url URL to download the reference from. Defaults to
#'   `getOption("ProstateAtlasMapR.reference_url")`.
#' @param sha256 Expected SHA-256 checksum (hex string) used to verify the
#'   downloaded file's integrity. Defaults to
#'   `getOption("ProstateAtlasMapR.reference_sha256")`. If `NULL`, the
#'   checksum step is skipped (with a warning).
#' @param cache_dir Directory used to cache the downloaded file. Defaults to
#'   a per-user cache directory via [tools::R_user_dir()].
#' @param force Force re-download even if a valid cached copy already
#'   exists. Default `FALSE`.
#' @param verbose Print progress messages.
#'
#' @return An object of class `patlas_reference`.
#' @export
get_reference <- function(path = NULL,
                           url = getOption("ProstateAtlasMapR.reference_url"),
                           sha256 = getOption("ProstateAtlasMapR.reference_sha256"),
                           cache_dir = tools::R_user_dir("ProstateAtlasMapR", "cache"),
                           force = FALSE,
                           verbose = TRUE) {
  if (!is.null(path)) {
    if (!file.exists(path)) {
      cli::cli_abort("File not found: {.path {path}}")
    }
    return(.load_reference(path))
  }

  if (is.null(url)) {
    cli::cli_abort(c(
      "No download URL available.",
      "i" = "Either pass {.arg path} to a local reference file, or set {.code options(ProstateAtlasMapR.reference_url = <url>)} (and optionally {.code ProstateAtlasMapR.reference_sha256})."
    ))
  }

  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  }
  dest <- file.path(cache_dir, .cache_filename(url))

  use_cache <- !force && file.exists(dest) && .checksum_ok(dest, sha256)
  if (use_cache) {
    if (verbose) cli::cli_alert_info("Using cached reference at {.path {dest}}.")
  } else {
    if (verbose) cli::cli_alert_info("Downloading reference from {.url {url}} ...")
    tmp <- paste0(dest, ".part")
    curl::curl_download(url, tmp, quiet = !verbose)
    if (!.checksum_ok(tmp, sha256)) {
      file.remove(tmp)
      cli::cli_abort("Downloaded file failed checksum verification; aborting.")
    }
    file.rename(tmp, dest)
  }

  .load_reference(dest)
}

#' @rdname get_reference
#' @export
download_reference <- get_reference

#' Derive a cache filename from a URL
#' @noRd
.cache_filename <- function(url) {
  bn <- basename(sub("[?#].*$", "", url))
  if (!nzchar(bn)) bn <- "reference.qs"
  if (!grepl("\\.qs$", bn)) bn <- paste0(bn, ".qs")
  bn
}

#' Verify a file's SHA-256 checksum, if one was supplied
#' @noRd
.checksum_ok <- function(file, sha256) {
  if (!file.exists(file)) return(FALSE)
  if (is.null(sha256)) return(TRUE)
  if (!requireNamespace("openssl", quietly = TRUE)) {
    cli::cli_warn("Package {.pkg openssl} is not installed; skipping checksum verification.")
    return(TRUE)
  }
  con <- file(file, "rb")
  on.exit(close(con), add = TRUE)
  actual <- as.character(openssl::sha256(con))
  identical(tolower(actual), tolower(sha256))
}
