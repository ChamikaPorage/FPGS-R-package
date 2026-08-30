#' FPGS Weighting Estimator
#'
#' Estimates the average treatment effect (ATE) using weighting methods based
#' on propensity scores estimated from the two-dimensional Full Prognostic
#' Score (FPGS).
#'
#' @param fit An object of class \code{"fpgs"} returned by \code{fpgs()}.
#' @param type Character string specifying the weighting estimator:
#'   \code{"ipw"}, \code{"nipw"}, or \code{"aipw"}.
#'
#' @return An object of class \code{"fpgs_weight"} containing the ATE estimate,
#'   treatment-specific means, estimated propensity scores, and fitted models.
#'
#' @details
#' The propensity score is estimated using the two estimated FPGS components.
#' For AIPW, additional treatment-specific outcome models are fitted using
#' the estimated FPGS components as predictors, and the resulting fitted
#' outcome means are used in the augmented estimator.
#'
#' @export
fpgs_weight <- function(fit, type = c("ipw", "nipw", "aipw")) {

  type <- match.arg(type)

  if (!inherits(fit, "fpgs")) {
    stop("fit must be an object returned by fpgs()")
  }

  data <- fit$data
  Y <- data[[fit$outcome]]
  Tr <- data[[fit$treatment]]

  fpgs_data <- data.frame(
    Y = Y,
    Tr = Tr,
    mu0_hat = fit$mu0_hat,
    mu1_hat = fit$mu1_hat
  )

  # FPGS-based propensity score model
  ps_fit <- stats::glm(
    Tr ~ mu0_hat + mu1_hat,
    data = fpgs_data,
    family = stats::binomial()
  )

  e_hat <- stats::predict(
    ps_fit,
    type = "response"
  )

  e_hat <- pmin(
    pmax(e_hat, 1e-6),
    1 - 1e-6
  )

  if (type == "ipw") {

    mu1 <- mean(
      Tr * Y / e_hat
    )

    mu0 <- mean(
      (1 - Tr) * Y / (1 - e_hat)
    )
  }

  if (type == "nipw") {

    mu1 <- sum(Tr * Y / e_hat) /
      sum(Tr / e_hat)

    mu0 <- sum((1 - Tr) * Y / (1 - e_hat)) /
      sum((1 - Tr) / (1 - e_hat))
  }

  outcome_model1 <- NULL
  outcome_model0 <- NULL
  outcome_mean1 <- NULL
  outcome_mean0 <- NULL

  if (type == "aipw") {

    # Fit treatment-specific outcome models using the estimated FPGS
    if (fit$outcome_type == "binary") {

      outcome_model1 <- stats::glm(
        Y ~ mu0_hat + mu1_hat,
        data = fpgs_data,
        subset = Tr == 1,
        family = stats::binomial()
      )

      outcome_model0 <- stats::glm(
        Y ~ mu0_hat + mu1_hat,
        data = fpgs_data,
        subset = Tr == 0,
        family = stats::binomial()
      )

    } else {

      outcome_model1 <- stats::lm(
        Y ~ mu0_hat + mu1_hat,
        data = fpgs_data,
        subset = Tr == 1
      )

      outcome_model0 <- stats::lm(
        Y ~ mu0_hat + mu1_hat,
        data = fpgs_data,
        subset = Tr == 0
      )
    }

    outcome_mean1 <- stats::predict(
      outcome_model1,
      newdata = fpgs_data,
      type = "response"
    )

    outcome_mean0 <- stats::predict(
      outcome_model0,
      newdata = fpgs_data,
      type = "response"
    )

    mu1 <- mean(
      outcome_mean1 +
        Tr * (Y - outcome_mean1) / e_hat
    )

    mu0 <- mean(
      outcome_mean0 +
        (1 - Tr) * (Y - outcome_mean0) / (1 - e_hat)
    )
  }

  ate <- mu1 - mu0

  out <- list(
    estimate = ate,
    mu1 = mu1,
    mu0 = mu0,
    propensity_score = e_hat,
    ps_model = ps_fit,
    outcome_model1 = outcome_model1,
    outcome_model0 = outcome_model0,
    outcome_mean1 = outcome_mean1,
    outcome_mean0 = outcome_mean0,
    type = type,
    method = paste("FPGS", type, "weighting")
  )

  class(out) <- "fpgs_weight"

  out
}
