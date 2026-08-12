#' @export
print.fpgs <- function(x, ...) {
  
  cat("\nFull Prognostic Score Model\n")
  cat("----------------------------------\n")
  
  cat("Outcome:", x$outcome, "\n")
  cat("Treatment:", x$treatment, "\n")
  
  if (!is.null(x$outcome_type))
    cat("Outcome type:", x$outcome_type, "\n")
  
  if (!is.null(x$method))
    cat("Method:", x$method, "\n")
  
  cat("Number of covariates:",
      length(x$covariates), "\n")
  
  invisible(x)
}