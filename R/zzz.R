.onLoad <- function(libname, pkgname) {
  ## Prebuilt lean reference is hosted separately (GitHub Releases / Zenodo).
  ## Fill these in after you upload prostate_atlas_reference.qs, e.g.:
  ##   options(
  ##     ProstateAtlasMapR.reference_url = "https://github.com/<org>/ProstateAtlasMapR/releases/download/v0.1.0/prostate_atlas_reference.qs",
  ##     ProstateAtlasMapR.reference_sha256 = "<sha256 hex>"
  ##   )
  op <- options()
  op.pkg <- list(
    ProstateAtlasMapR.reference_url = NULL,
    ProstateAtlasMapR.reference_sha256 = NULL
  )
  toset <- !(names(op.pkg) %in% names(op))
  if (any(toset)) options(op.pkg[toset])
  invisible()
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "ProstateAtlasMapR: use get_reference(path = ...) for a local lean reference, ",
    "or set options(ProstateAtlasMapR.reference_url = <url>) after publishing the data release."
  )
}
