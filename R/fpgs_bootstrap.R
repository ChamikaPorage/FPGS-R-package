#' Bootstrap Inference for FPGS Estimators
#'
#' Computes bootstrap standard errors and percentile confidence intervals
#' for FPGS treatment effect estimators.
#'
#' @param fit An object of class \code{"fpgs"} returned by \code{fpgs()}.
#' @param estimator Character string specifying the estimator. One of
#'   \code{"ri"}, \code{"matching"}, \code{"ipw"}, \code{"nipw"},
#'   \code{"aipw"}, or \code{"stratification"}.
#' @param B Number of bootstrap replications. Default is 500.
#' @param conf.level Confidence level for the percentile bootstrap interval.
#'   Default is 0.95.
#' @param seed Optional random seed.
#' @param verbose Logical. If \code{TRUE}, prints bootstrap progress.
#' @param ... Additional arguments passed to the treatment effect estimator,
#'   such as \code{M} for \code{fpgs_match()} or \code{n_bins} and
#'   \code{min_cell} for \code{fpgs_stratify()}.
#'
#' @return An object of class \code{"fpgs_bootstrap"} containing the point
#'   estimate, bootstrap variance and standard error, percentile confidence
#'   interval, bootstrap estimates, and number of valid bootstrap replications.
#'
#' @details
#' In each bootstrap replication, observations are sampled with replacement.
#' The FPGS is then re-estimated using the same outcome type and estimation
#' method as in the original fit, after which the selected treatment effect
#' estimator is recomputed.
#'
#' Bootstrap standard errors are calculated as the standard deviation of the
#' valid bootstrap estimates. Percentile confidence intervals are obtained
#' from the corresponding empirical bootstrap quantiles.
#'
#' @export
fpgs_bootstrap <- function(fit,
                           estimator = c(
                             "ri",
                             "matching",
                             "ipw",
                             "nipw",
                             "aipw",
                             "stratification"
                           ),
                           B = 500,
                           conf.level = 0.95,
                           seed = NULL,
                           verbose = TRUE,
                           ...) {

  estimator <- match.arg(estimator)

  if (!inherits(fit, "fpgs")) {
    stop("fit must be an object returned by fpgs()", call. = FALSE)
  }

  if (!is.numeric(B) || length(B) != 1 || B < 2) {
    stop("B must be a number greater than or equal to 2", call. = FALSE)
  }

  if (!is.numeric(conf.level) || length(conf.level) != 1 ||
      conf.level <= 0 || conf.level >= 1) {
    stop("conf.level must be a number between 0 and 1", call. = FALSE)
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  data <- as.data.frame(fit$data)
  n <- nrow(data)

  estimate <- .fpgs_apply_estimator(fit, estimator, ...)

  bootstrap_estimates <- rep(NA_real_, B)

  for (b in seq_len(B)) {

    if (verbose && (b == 1 || b %% 25 == 0 || b == B)) {
      message("Bootstrap replication ", b, " of ", B)
    }

    bootstrap_index <- sample(seq_len(n), size = n, replace = TRUE)
    bootstrap_data <- data[bootstrap_index, , drop = FALSE]

    bootstrap_estimates[b] <- tryCatch(
      {
        bootstrap_fit <- .fpgs_refit(fit, bootstrap_data)

        .fpgs_apply_estimator(
          bootstrap_fit,
          estimator,
          ...
        )
      },
      error = function(e) {
        NA_real_
      }
    )
  }

  valid_estimates <- bootstrap_estimates[
    is.finite(bootstrap_estimates)
  ]

  B_valid <- length(valid_estimates)

  if (B_valid < 2) {
    stop(
      "Fewer than two valid bootstrap estimates were obtained.",
      call. = FALSE
    )
  }

  if (B_valid < 0.9 * B) {
    warning(
      "More than 10% of bootstrap replications failed. ",
      "Bootstrap inference may be unreliable.",
      call. = FALSE
    )
  }

  se <- stats::sd(valid_estimates)
  variance <- se^2

  alpha <- 1 - conf.level

  conf_int <- stats::quantile(
    valid_estimates,
    probs = c(alpha / 2, 1 - alpha / 2),
    na.rm = TRUE,
    names = FALSE
  )

  out <- list(
    estimate = estimate,
    variance = variance,
    se = se,
    conf.int = conf_int,
    conf.level = conf.level,
    boot_estimates = valid_estimates,
    boot_estimates_all = bootstrap_estimates,
    B = as.integer(B),
    B_valid = B_valid,
    estimator = estimator,
    outcome_type = fit$outcome_type,
    fpgs_method = fit$method,
    method = paste("Bootstrap inference for FPGS", estimator)
  )

  class(out) <- "fpgs_bootstrap"

  out
}


# Internal helper for refitting the FPGS in a bootstrap sample
.fpgs_refit <- function(fit, data) {

  if (fit$outcome_type == "continuous" &&
      fit$method == "parametric") {

    return(
      fpgs_continuous_parametric(
        data = data,
        outcome = fit$outcome,
        treatment = fit$treatment,
        covariates = fit$covariates
      )
    )
  }

  if (fit$outcome_type == "binary" &&
      fit$method == "parametric") {

    return(
      fpgs_binary_parametric(
        data = data,
        outcome = fit$outcome,
        treatment = fit$treatment,
        covariates = fit$covariates
      )
    )
  }

  if (fit$outcome_type == "continuous" &&
      fit$method == "rf") {

    return(
      fpgs_continuous_rf(
        data = data,
        outcome = fit$outcome,
        treatment = fit$treatment,
        covariates = fit$covariates,
        folds = fit$folds,
        mtry = fit$mtry
      )
    )
  }

  if (fit$outcome_type == "binary" &&
      fit$method == "rf") {

    return(
      fpgs_binary_rf(
        data = data,
        outcome = fit$outcome,
        treatment = fit$treatment,
        covariates = fit$covariates,
        folds = fit$folds,
        num.trees = fit$num.trees,
        mtry = fit$mtry,
        min.node.size = fit$min.node.size
      )
    )
  }

  stop("Unsupported FPGS object.", call. = FALSE)
}


# Internal helper for applying one treatment effect estimator
.fpgs_apply_estimator <- function(fit, estimator, ...) {

  if (estimator == "ri") {
    return(fpgs_ri(fit)$estimate)
  }

  if (estimator == "matching") {
    return(fpgs_match(fit, ...)$estimate)
  }

  if (estimator %in% c("ipw", "nipw", "aipw")) {
    return(fpgs_weight(fit, type = estimator)$estimate)
  }

  if (estimator == "stratification") {
    return(fpgs_stratify(fit, ...)$estimate)
  }

  stop("Unknown estimator.", call. = FALSE)
}


#' Print Method for FPGS Bootstrap Inference
#'
#' @param x An object of class \code{"fpgs_bootstrap"}.
#' @param ... Additional arguments, currently ignored.
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
print.fpgs_bootstrap <- function(x, ...) {

  cat("\n", x$method, "\n", sep = "")
  cat("----------------------------------\n")

  cat("Estimator:", x$estimator, "\n")
  cat("FPGS method:", x$fpgs_method, "\n")
  cat("Outcome type:", x$outcome_type, "\n")
  cat("Estimate:", x$estimate, "\n")
  cat("Bootstrap SE:", x$se, "\n")
  cat(
    paste0(100 * x$conf.level, "% CI: ["),
    x$conf.int[1], ", ", x$conf.int[2], "]\n",
    sep = ""
  )
  cat(
    "Bootstrap replications:",
    x$B_valid, "valid out of", x$B, "\n"
  )

  invisible(x)
}


#' Summary Method for FPGS Bootstrap Inference
#'
#' @param object An object of class \code{"fpgs_bootstrap"}.
#' @param ... Additional arguments, currently ignored.
#'
#' @return A data frame containing the bootstrap inference results.
#'
#' @export
summary.fpgs_bootstrap <- function(object, ...) {

  data.frame(
    estimator = object$estimator,
    fpgs_method = object$fpgs_method,
    outcome_type = object$outcome_type,
    estimate = object$estimate,
    variance = object$variance,
    se = object$se,
    ci_lower = object$conf.int[1],
    ci_upper = object$conf.int[2],
    conf_level = object$conf.level,
    B = object$B,
    B_valid = object$B_valid
  )
}


#' Confidence Interval for FPGS Bootstrap Inference
#'
#' @param object An object of class \code{"fpgs_bootstrap"}.
#' @param parm Currently ignored.
#' @param level Confidence level. Defaults to the level stored in
#'   \code{object}.
#' @param ... Additional arguments, currently ignored.
#'
#' @return A numeric vector containing the lower and upper confidence limits.
#'
#' @export
confint.fpgs_bootstrap <- function(object,
                                   parm = NULL,
                                   level = object$conf.level,
                                   ...) {

  if (!is.numeric(level) || length(level) != 1 ||
      level <= 0 || level >= 1) {
    stop("level must be a number between 0 and 1", call. = FALSE)
  }

  alpha <- 1 - level

  conf_int <- stats::quantile(
    object$boot_estimates,
    probs = c(alpha / 2, 1 - alpha / 2),
    na.rm = TRUE,
    names = FALSE
  )

  names(conf_int) <- c("lower", "upper")

  conf_int
}


#' Extract FPGS Bootstrap Estimate
#'
#' @param object An object of class \code{"fpgs_bootstrap"}.
#' @param ... Additional arguments, currently ignored.
#'
#' @return The estimated average treatment effect.
#'
#' @export
coef.fpgs_bootstrap <- function(object, ...) {

  object$estimate
}


#' Variance-Covariance Matrix for FPGS Bootstrap Estimate
#'
#' @param object An object of class \code{"fpgs_bootstrap"}.
#' @param ... Additional arguments, currently ignored.
#'
#' @return A one-by-one matrix containing the bootstrap variance estimate.
#'
#' @export
vcov.fpgs_bootstrap <- function(object, ...) {

  variance_matrix <- matrix(
    object$variance,
    nrow = 1,
    ncol = 1
  )

  rownames(variance_matrix) <- "ATE"
  colnames(variance_matrix) <- "ATE"

  variance_matrix
}
