#' FPGS Regression Imputation Estimator
#'
#' @param fit An object returned by fpgs()
#'
#' @return An object of class "fpgs_ri"
#'
#' @export
fpgs_ri <- function(fit) {
  
  if (!inherits(fit, "fpgs")) {
    stop("fit must be an object returned by fpgs()")
  }
  
  tau_i <- fit$mu1_hat - fit$mu0_hat
  
  ate <- mean(tau_i)
  
  out <- list(
    estimate = ate,
    mu1_hat = fit$mu1_hat,
    mu0_hat = fit$mu0_hat,
    individual_effects = tau_i,
    method = "FPGS regression imputation"
  )
  
  class(out) <- "fpgs_ri"
  out
}