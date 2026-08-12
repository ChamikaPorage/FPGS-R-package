#' Print an FPGS Object
#'
#' Prints a concise summary of an estimated Full Prognostic Score object.
#'
#' @param x An object of class \code{"fpgs"} returned by \code{fpgs()}.
#' @param ... Additional arguments passed to the print method.
#'
#' @return The input object \code{x}, returned invisibly.
#'
#' @export
print.fpgs <- function(x, ...) {
  
  cat("\nFull Prognostic Score Model\n")
  cat("----------------------------------\n")
  
  cat("Outcome:", x$outcome, "\n")
  cat("Treatment:", x$treatment, "\n")
  
  if (!is.null(x$outcome_type)) {
    cat("Outcome type:", x$outcome_type, "\n")
  }
  
  if (!is.null(x$method)) {
    cat("Method:", x$method, "\n")
  }
  
  if (!is.null(x$folds)) {
    cat("Cross-fitting folds:", x$folds, "\n")
  }
  
  cat("Number of covariates:", length(x$covariates), "\n")
  
  invisible(x)
}
