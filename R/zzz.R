.onLoad <- function(libname, pkgname) {
  ## Default: GitHub Release v0.1.0 lean reference (~0.9 GB).
  ## Users may override with options(...) before or after loading.
  op <- options()
  op.pkg <- list(
    ProstateAtlasMapR.reference_url = .default_reference_url(),
    ProstateAtlasMapR.reference_sha256 = .default_reference_sha256()
  )
  toset <- !(names(op.pkg) %in% names(op))
  if (any(toset)) {
    options(op.pkg[toset])
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
