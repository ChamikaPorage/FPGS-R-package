#' Summarize an FPGS Object
#'
#' Computes and displays treatment effect estimates obtained using the
#' available FPGS-based estimators.
#'
#' @param object An object of class \code{"fpgs"} returned by \code{fpgs()}.
#' @param ... Additional arguments passed to \code{fpgs_all()}.
#'
#' @return A data frame containing the FPGS-based treatment effect estimates,
#'   returned invisibly.
#'
#' @export
summary.fpgs <- function(object, ...) {
  
  res <- fpgs_all(object, ...)
  
  cat("\nFPGS Treatment Effect Estimates\n")
  cat("----------------------------------\n")
  
  print(res)
  
  invisible(res)
}
