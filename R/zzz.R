#' @importFrom rlang .data
NULL

.onLoad <- function(libname, pkgname) {
  prolfqua::register_facade(
    "lm_missing",
    class = "ContrastsLMMissingFacade",
    needs = "same",
    package = "prolfquatakeoutmissing"
  )
  invisible()
}
