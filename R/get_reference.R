#' Default hosted lean-reference URL (GitHub Release v0.1.0)
#' @noRd
.default_reference_url <- function() {
  "https://github.com/FY10086/ProstateAtlasMapR/releases/download/v0.1.0/prostate_atlas_reference.qs"
}

#' Default SHA-256 for [.default_reference_url()]
#' @noRd
.default_reference_sha256 <- function() {
  "c83d94b5326aad495cb854f98fbb31ca10a0b1143eaf45bab055f368a8071d5f"
}

#' Get the prebuilt prostate atlas reference (download + cache, or load local)
#'
#' The lean reference produced by [build_reference()] is hosted as a
#' [GitHub Release](https://github.com/FY10086/ProstateAtlasMapR/releases)
#' asset (`prostate_atlas_reference.qs`, ~0.9 GB), not inside the R package
#' source. By default this function downloads that file once, verifies its
#' SHA-256 checksum, caches it under [tools::R_user_dir()], and returns the
#' loaded `patlas_reference` object. Pass `path` to use a local copy instead.
#'
#' @param path Path to an existing local reference file. If supplied, the
#'   file is read directly and no download is attempted.
#' @param url URL to download the reference from. Defaults to
#'   `getOption("ProstateAtlasMapR.reference_url")`, which is set on package
#'   load to the v0.1.0 GitHub Release asset.
#' @param sha256 Expected SHA-256 checksum (hex string). Defaults to
#'   `getOption("ProstateAtlasMapR.reference_sha256")`. Verification requires
#'   the **openssl** package. If `NULL`, checksum verification is skipped
#'   (with a warning).
#' @param cache_dir Directory used to cache the downloaded file. Defaults to
#'   a per-user cache directory via [tools::R_user_dir()].
#' @param force Force re-download even if a valid cached copy already
#'   exists. Default `FALSE`.
#' @param verbose Print progress messages.
#'
#' @return An object of class `patlas_reference`.
#'
#' @examples
#' \dontrun{
#' ## Download once (~0.9 GB), then reuse the cache:
#' ref <- get_reference()
#'
#' ## Or point at a local copy:
#' ref <- get_reference(path = "prostate_atlas_reference.qs")
#' }
#'
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
    if (verbose) {
      cli::cli_alert_info("Loading local reference from {.path {path}}.")
    }
    return(.load_reference(path))
  }

  if (is.null(url) || !nzchar(url)) {
    cli::cli_abort(c(
      "No download URL available.",
      "i" = "Pass {.arg path} to a local {.code .qs} file, or set {.code options(ProstateAtlasMapR.reference_url = <url>)}.",
      "i" = "Default release: {.url https://github.com/FY10086/ProstateAtlasMapR/releases/tag/v0.1.0}"
    ))
  }

  if (is.null(sha256) || !nzchar(sha256)) {
    cli::cli_warn(c(
      "No SHA-256 provided; downloaded file will not be integrity-checked.",
      "i" = "Set {.code options(ProstateAtlasMapR.reference_sha256 = <hex>)} or pass {.arg sha256}."
    ))
    sha256 <- NULL
  }

  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  }
  dest <- file.path(cache_dir, .cache_filename(url))

  use_cache <- isFALSE(force) && file.exists(dest) && .checksum_ok(dest, sha256)
  if (use_cache) {
    if (verbose) {
      cli::cli_alert_success("Using cached reference at {.path {dest}}.")
    }
  } else {
    if (file.exists(dest) && isTRUE(force) && verbose) {
      cli::cli_alert_info("Re-downloading reference ({.arg force = TRUE}) ...")
    } else if (file.exists(dest) && verbose) {
      cli::cli_alert_warning("Cached file failed checksum; re-downloading ...")
    } else if (verbose) {
      cli::cli_alert_info(
        "Downloading reference (~0.9 GB) from {.url {url}} ..."
      )
    }
    .download_reference_file(url, dest, sha256 = sha256, verbose = verbose)
    if (verbose) {
      sz <- file.info(dest)$size
      cli::cli_alert_success(
        "Saved to {.path {dest}} ({(.format_bytes(sz))})."
      )
    }
  }

  .load_reference(dest)
}

#' Download + verify into dest (atomic .part then rename)
#' @noRd
.download_reference_file <- function(url, dest, sha256 = NULL, verbose = TRUE) {
  tmp <- paste0(dest, ".part")
  if (file.exists(tmp)) {
    unlink(tmp)
  }
  on.exit({
    if (file.exists(tmp)) unlink(tmp)
  }, add = TRUE)

  if (grepl("^file://", url, ignore.case = TRUE)) {
    ## Local file URLs (used in tests): copy instead of curl
    src <- sub("^file://", "", url, ignore.case = TRUE)
    src <- sub("^localhost", "", src)
    if (!file.exists(src)) {
      cli::cli_abort("Local file URL not found: {.path {src}}")
    }
    ok <- file.copy(src, tmp, overwrite = TRUE)
    if (!isTRUE(ok)) {
      cli::cli_abort("Failed to copy local reference from {.path {src}}.")
    }
  } else {
    h <- curl::new_handle()
    ## Large asset: no overall timeout; follow GitHub Release redirects
    curl::handle_setopt(
      h,
      followlocation = TRUE,
      connecttimeout = 60,
      timeout = 0,
      noprogress = !verbose
    )

    tryCatch(
      curl::curl_download(url, tmp, quiet = !verbose, handle = h),
      error = function(e) {
        cli::cli_abort(c(
          "Download failed: {conditionMessage(e)}",
          "i" = "URL: {.url {url}}",
          "i" = "Or download manually and pass {.code get_reference(path = ...)}."
        ))
      }
    )
  }

  if (!file.exists(tmp) || isTRUE(file.info(tmp)$size < 1)) {
    cli::cli_abort("Download produced an empty file; aborting.")
  }

  if (!.checksum_ok(tmp, sha256)) {
    unlink(tmp)
    cli::cli_abort(c(
      "Downloaded file failed checksum (SHA-256) verification; aborting.",
      "i" = "Expected: {.val {sha256}}",
      "i" = "Re-try with {.code force = TRUE}, or re-download from the Release page."
    ))
  }

  if (file.exists(dest)) {
    unlink(dest)
  }
  if (!file.rename(tmp, dest)) {
    ## cross-device rename fallback
    file.copy(tmp, dest, overwrite = TRUE)
    unlink(tmp)
  }
  invisible(dest)
}

#' Derive a cache filename from a URL
#' @noRd
.cache_filename <- function(url) {
  bn <- basename(sub("[?#].*$", "", url))
  if (!nzchar(bn)) bn <- "prostate_atlas_reference.qs"
  if (!grepl("\\.qs$", bn, ignore.case = TRUE)) {
    bn <- paste0(bn, ".qs")
  }
  bn
}

#' Human-readable byte size
#' @noRd
.format_bytes <- function(n) {
  if (is.na(n) || n < 0) {
    return("unknown size")
  }
  units <- c("B", "KB", "MB", "GB", "TB")
  i <- 1L
  x <- as.numeric(n)
  while (x >= 1024 && i < length(units)) {
    x <- x / 1024
    i <- i + 1L
  }
  sprintf("%.1f %s", x, units[[i]])
}

#' Verify a file's SHA-256 checksum, if one was supplied
#' @noRd
.checksum_ok <- function(file, sha256) {
  if (!file.exists(file)) {
    return(FALSE)
  }
  if (is.null(sha256) || !nzchar(sha256)) {
    return(TRUE)
  }
  if (!requireNamespace("openssl", quietly = TRUE)) {
    cli::cli_abort(c(
      "Package {.pkg openssl} is required to verify the reference checksum.",
      "i" = "Install with {.code install.packages(\"openssl\")}."
    ))
  }
  con <- file(file, "rb")
  on.exit(close(con), add = TRUE)
  ## openssl returns a 'hash' character; strip class on both sides
  actual <- tolower(as.vector(as.character(openssl::sha256(con))))
  expected <- tolower(as.vector(as.character(sha256)))
  identical(actual, expected)
}
