.onLoad <- function(libname, pkgname) {
  ## Always fill missing/NULL options with the GitHub Release defaults.
  ## (Old sessions may have left these names unset or NULL.)
  if (is.null(getOption("ProstateAtlasMapR.reference_url")) ||
      !nzchar(getOption("ProstateAtlasMapR.reference_url"))) {
    options(ProstateAtlasMapR.reference_url = .default_reference_url())
  }
  if (is.null(getOption("ProstateAtlasMapR.reference_sha256")) ||
      !nzchar(getOption("ProstateAtlasMapR.reference_sha256"))) {
    options(ProstateAtlasMapR.reference_sha256 = .default_reference_sha256())
  }
  invisible()
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "ProstateAtlasMapR: get_reference() downloads the lean atlas (~0.9 GB) ",
    "from GitHub Releases on first use (cached afterwards). ",
    "Or pass get_reference(path = \".../prostate_atlas_reference.qs\")."
  )
}
