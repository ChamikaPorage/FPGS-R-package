#' FPGS Matching Estimator
#'
#' Estimates the average treatment effect (ATE) by matching on the
#' two-dimensional estimated Full Prognostic Score.
#'
#' @param fit An object of class \code{"fpgs"} returned by \code{fpgs()}.
#' @param M Number of matches used for each unit. Default is 1.
#'
#' @return An object of class \code{"fpgs_match"} containing the ATE estimate,
#'   estimated variance and standard error, number of matches, and the
#'   matching object returned by \code{Matching::Match()}.
#'
#' @details
#' Matching is performed directly on the two-dimensional estimated FPGS using
#' \code{Matching::Match()}. The analytic standard error is obtained directly
#' from \code{Matching::Match()}, and the variance is calculated as the squared
#' standard error.
#'
#' For FPGS estimated using random forests, bootstrap inference may be used
#' to account for uncertainty in estimation of the FPGS.
#'
#' @export
fpgs_match <- function(fit, M = 1) {

  if (!inherits(fit, "fpgs")) {
    stop("fit must be an object returned by fpgs()")
  }

  if (!requireNamespace("Matching", quietly = TRUE)) {
    stop(
      "Package 'Matching' is required. ",
      "Install it with install.packages('Matching')."
    )
  }

  dat <- fit$data
  Y <- dat[[fit$outcome]]
  Tr <- dat[[fit$treatment]]

  X <- as.matrix(fit$fpgs)
  storage.mode(X) <- "double"

  match_obj <- Matching::Match(
    Y = Y,
    Tr = Tr,
    X = X,
    estimand = "ATE",
    M = M
  )

  estimate <- unname(match_obj$est)
  se <- unname(match_obj$se.standard)
  variance <- se^2

  out <- list(
    estimate = estimate,
    variance = variance,
    se = se,
    M = M,
    match_object = match_obj,
    method = "FPGS matching",
    estimand = "ATE"
  )

  class(out) <- "fpgs_match"

  out
}
