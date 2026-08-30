#' Summarize an FPGS Fit
#'
#' Displays treatment effect estimates obtained from the available
#' FPGS estimators.
#'
#' @param object An object of class \code{"fpgs"} returned by \code{fpgs()}.
#' @param ... Additional arguments passed to the summary method.
#'
#' @return Invisibly returns the results from \code{fpgs_all()}.
#'
#' @export
summary.fpgs <- function(object, ...) {

  res <- fpgs_ate(object)

  cat("\nFPGS Estimates\n")
  cat("----------------------------------\n")

  print(res)

  invisible(res)
}
