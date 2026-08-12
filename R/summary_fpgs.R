#' @export
summary.fpgs <- function(object, ...) {
  
  res <- fpgs_all(object)
  
  cat("\nFPGS Estimates\n")
  cat("----------------------------------\n")
  
  print(res)
  
  invisible(res)
}