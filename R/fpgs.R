#' Estimate the Full Prognostic Score
#'
#' Estimates the two-dimensional Full Prognostic Score (FPGS) using
#' parametric or random forest outcome models for continuous or binary
#' outcomes.
#'
#' @param data A data frame containing the observed data.
#' @param outcome Character string giving the name of the outcome variable.
#' @param treatment Character string giving the name of the treatment variable.
#' @param covariates Character vector giving the names of the pre treatment variables.
#' @param outcome_type Type of outcome: \code{"continuous"} or \code{"binary"}.
#' @param method Estimation method: \code{"parametric"} or \code{"rf"}.
#'   where \code{"rf"} denotes random forest.
#' @param folds Number of folds used for cross-fitting when
#'   \code{method = "rf"}. Default is 5.
#' @param ... Additional arguments passed to the random forest fitting
#'   functions when \code{method = "rf"}, such as \code{num.trees},
#'   \code{mtry}, and \code{min.node.size}.
#'
#' @return An object of class \code{"fpgs"} containing the estimated
#'   FPGS and information required by the downstream treatment effect
#'   estimators.
#'
#' @export
fpgs <- function(data,
                 outcome,
                 treatment,
                 covariates = NULL,
                 outcome_type = c("continuous", "binary"),
                 method = c("parametric", "rf"),
                 folds = 5,
                 ...) {

  outcome_type <- match.arg(outcome_type)
  method <- match.arg(method)

  if (outcome_type == "continuous" && method == "parametric") {

    fit <- fpgs_continuous_parametric(
      data,
      outcome,
      treatment,
      covariates
    )
  }

  if (outcome_type == "continuous" && method == "rf") {

    fit <- fpgs_continuous_rf(
      data,
      outcome,
      treatment,
      covariates,
      folds = folds,
      ...
    )
  }

  if (outcome_type == "binary" && method == "parametric") {

    fit <- fpgs_binary_parametric(
      data,
      outcome,
      treatment,
      covariates
    )
  }

  if (outcome_type == "binary" && method == "rf") {

    fit <- fpgs_binary_rf(
      data,
      outcome,
      treatment,
      covariates,
      folds = folds,
      ...
    )
  }

  # Store information needed by downstream estimators
  fit$outcome_type <- outcome_type
  fit$method <- method

  fit
}
