#' FPGS Weighting Estimator
#'
#' Estimates ATE using propensity scores estimated from FPGS.
#'
#' @param fit An object returned by fpgs()
#' #' @param type Weighting estimator: "ipw", "normalized", or "aipw"
#'
#' @export
fpgs_weight <- function(fit, type = c("ipw", "normalized", "aipw")) {
  
  type <- match.arg(type)
  
  if (!inherits(fit, "fpgs")) {
    stop("fit must be an object returned by fpgs()")
  }
  
  dat <- fit$data
  Y <- dat[[fit$outcome]]
  Tr <- dat[[fit$treatment]]
  
  fpgs_dat <- data.frame(
    Tr = Tr,
    mu0_hat = fit$mu0_hat,
    mu1_hat = fit$mu1_hat
  )
  
  ps_fit <- glm(
    Tr ~ mu0_hat + mu1_hat,
    data = fpgs_dat,
    family = binomial()
  )
  
  e_hat <- predict(ps_fit, type = "response")
  
  e_hat <- pmin(pmax(e_hat, 1e-6), 1 - 1e-6)
  
  if (type == "ipw") {
    mu1 <- mean(Tr * Y / e_hat)
    mu0 <- mean((1 - Tr) * Y / (1 - e_hat))
  }
  
  if (type == "normalized") {
    mu1 <- sum(Tr * Y / e_hat) / sum(Tr / e_hat)
    mu0 <- sum((1 - Tr) * Y / (1 - e_hat)) / sum((1 - Tr) / (1 - e_hat))
  }
  
  if (type == "aipw") {
    m1 <- fit$mu1_hat
    m0 <- fit$mu0_hat
    
    mu1 <- mean(m1 + Tr * (Y - m1) / e_hat)
    mu0 <- mean(m0 + (1 - Tr) * (Y - m0) / (1 - e_hat))
  }
  
  ate <- mu1 - mu0
  
  
  out <- list(
    estimate = ate,
    mu1 = mu1,
    mu0 = mu0,
    propensity_score = e_hat,
    ps_model = ps_fit,
    type = type,
    method = paste("FPGS", type, "weighting")
  )
  
  class(out) <- "fpgs_weight"
  
  out
}