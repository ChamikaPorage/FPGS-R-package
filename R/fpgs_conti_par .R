#' Estimate Full Prognostic Scores
#'
#' Estimates the prognostic scores
#' E[Y(0)|X] and E[Y(1)|X]
#'
#' @param data Data frame
#' @param outcome Outcome variable name
#' @param treatment Treatment variable name
#' @param covariates Covariate names
#'
#' @return An object of class "fpgs"
#'
#' @export

fpgs_continuous_parametric <- function(data,
                                       outcome,
                                       treatment,
                                       covariates = NULL) {
  
  data <- as.data.frame(data)
  
  if (!outcome %in% names(data)) {
    stop("outcome variable not found in data")
  }
  
  if (!treatment %in% names(data)) {
    stop("treatment variable not found in data")
  }
  
  if (is.null(covariates)) {
    covariates <- setdiff(names(data), c(outcome, treatment))
  }
  
  if (length(covariates) == 0) {
    stop("No covariates found.")
  }
  
  tr <- data[[treatment]]
  
  if (!all(tr %in% c(0, 1))) {
    stop("treatment must be coded as 0/1")
  }
  
  form <- as.formula(
    paste(outcome, "~", paste(covariates, collapse = " + "))
  )
  
  fit1 <- lm(form, data = data[tr == 1, , drop = FALSE])
  fit0 <- lm(form, data = data[tr == 0, , drop = FALSE])
  
  mu1_hat <- predict(fit1, newdata = data)
  mu0_hat <- predict(fit0, newdata = data)
  
  out <- list(
    data = data,
    outcome = outcome,
    treatment = treatment,
    covariates = covariates,
    mu0_hat = mu0_hat,
    mu1_hat = mu1_hat,
    fpgs = data.frame(mu0_hat = mu0_hat, mu1_hat = mu1_hat),
    outcome_type = "continuous",
    method = "parametric"
  )
  
  class(out) <- "fpgs"
  out
}