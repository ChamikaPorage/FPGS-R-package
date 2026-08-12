#' FPGS Regression Imputation Estimator
#'
#' Estimates the average treatment effect (ATE) using regression imputation
#' based on the estimated FPGS.
#'
#' @param fit An object of class \code{"fpgs"} returned by \code{fpgs()}.
#'
#' @return An object of class \code{"fpgs_ri"} containing the ATE estimate,
#'   estimated conditional outcome means, and individual treatment effects.
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
