#' Estimate Binary FPGS Using Parametric Regression
#'
#' @param data Data frame containing the observed data.
#' @param outcome Binary outcome variable name.
#' @param treatment Treatment variable name.
#' @param covariates Covariate names. If NULL, all variables other than
#'   the outcome and treatment are used.
#'
#' @return An object of class \code{"fpgs"}.
#'
#' @keywords internal

fpgs_binary <- function(data,
                        outcome,
                        treatment,
                        covariates = NULL) {
  
  dat <- as.data.frame(data)
  
  if (is.null(covariates)) {
    covariates <- setdiff(names(dat), c(outcome, treatment))
  }
  
  tr <- dat[[treatment]]
  
  form <- as.formula(
    paste(outcome, "~", paste(covariates, collapse = " + "))
  )
  
  fit1 <- glm(
    form,
    data = dat[tr == 1, , drop = FALSE],
    family = binomial()
  )
  
  fit0 <- glm(
    form,
    data = dat[tr == 0, , drop = FALSE],
    family = binomial()
  )
  
  mu1_hat <- predict(fit1, newdata = dat, type = "response")
  mu0_hat <- predict(fit0, newdata = dat, type = "response")
  
  out <- list(
    data = dat,
    outcome = outcome,
    treatment = treatment,
    covariates = covariates,
    mu0_hat = mu0_hat,
    mu1_hat = mu1_hat,
    fpgs = data.frame(
      mu0_hat = mu0_hat,
      mu1_hat = mu1_hat
    ),
    outcome_type = "binary",
    method = "parametric"
  )
  
  class(out) <- "fpgs"
  
  out
}
