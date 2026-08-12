#' Bootstrap inference for FPGS estimators
#'
#' Computes bootstrap standard errors and confidence intervals for FPGS
#' estimators by resampling observations, refitting the FPGS model, and
#' recomputing the selected treatment effect estimator in each bootstrap sample.
#'
#' @param fit An object returned by [fpgs()].
#' @param estimator Character string specifying the estimator. One of
#'   "ri", "matching", "stratification", "ipw", "normalized", or "aipw".
#' @param B Number of bootstrap replications. Default is 500.
#' @param conf.level Confidence level for the percentile bootstrap interval.
#'   Default is 0.95.
#' @param seed Optional random seed.
#' @param verbose Logical. If TRUE, prints progress.
#' @param ... Additional arguments passed to estimator functions, for example
#'   M for [fpgs_match()] or n_strata for [fpgs_stratify()].
#'
#' @return An object of class "fpgs_bootstrap".
#'
#' @export
fpgs_bootstrap <- function(fit,
                           estimator = c("ri", "matching",
                                         "ipw", "normalized", "aipw", "stratification"),
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

  dat <- as.data.frame(fit$data)
  n <- nrow(dat)

  theta_hat <- .fpgs_apply_estimator(fit, estimator, ...)

  boot_est <- rep(NA_real_, B)

  for (b in seq_len(B)) {

    if (verbose && (b == 1 || b %% 25 == 0 || b == B)) {
      message("Bootstrap replication ", b, " of ", B)
    }

    idx <- sample(seq_len(n), size = n, replace = TRUE)
    dat_b <- dat[idx, , drop = FALSE]

    boot_est[b] <- tryCatch({
      fit_b <- .fpgs_refit(fit, dat_b)
      .fpgs_apply_estimator(fit_b, estimator, ...)
    }, error = function(e) {
      NA_real_
    })
  }

  boot_est_valid <- boot_est[is.finite(boot_est)]

  if (length(boot_est_valid) < 2) {
    stop("Fewer than two valid bootstrap estimates were obtained.", call. = FALSE)
  }

  alpha <- 1 - conf.level

  se <- stats::sd(boot_est_valid)
  ci <- stats::quantile(
    boot_est_valid,
    probs = c(alpha / 2, 1 - alpha / 2),
    na.rm = TRUE,
    names = FALSE
  )

  out <- list(
    estimate = theta_hat,
    se = se,
    conf.int = ci,
    conf.level = conf.level,
    boot_estimates = boot_est_valid,
    boot_estimates_all = boot_est,
    B = as.integer(B),
    B_valid = length(boot_est_valid),
    estimator = estimator,
    outcome_type = fit$outcome_type,
    fpgs_method = fit$method,
    method = paste("Bootstrap inference for FPGS", estimator)
  )

  class(out) <- "fpgs_bootstrap"
  out
}


# Internal helper: refit the same FPGS model on a bootstrap sample.
.fpgs_refit <- function(fit, data) {

  if (fit$outcome_type == "continuous" && fit$method == "parametric") {
    return(
      fpgs_continuous_parametric(
        data = data,
        outcome = fit$outcome,
        treatment = fit$treatment,
        covariates = fit$covariates
      )
    )
  }

  if (fit$outcome_type == "binary" && fit$method == "parametric") {
    return(
      fpgs_binary(
        data = data,
        outcome = fit$outcome,
        treatment = fit$treatment,
        covariates = fit$covariates
      )
    )
  }

  if (fit$outcome_type == "continuous" && fit$method == "rf") {
    return(
      fpgs_rf(
        data = data,
        outcome = fit$outcome,
        treatment = fit$treatment,
        covariates = fit$covariates,
        folds = fit$folds,
        mtry = fit$mtry
      )
    )
  }

  if (fit$outcome_type == "binary" && fit$method == "rf") {
    return(
      fpgs_binary_rf(
        data = data,
        outcome = fit$outcome,
        treatment = fit$treatment,
        covariates = fit$covariates,
        folds = fit$folds,
        mtry = fit$mtry
      )
    )
  }

  stop("Unsupported FPGS object.", call. = FALSE)
}


# Internal helper: apply one ATE estimator to an FPGS object.
.fpgs_apply_estimator <- function(fit, estimator, ...) {

  if (estimator == "ri") {
    return(fpgs_ri(fit)$estimate)
  }

  if (estimator == "matching") {
    return(fpgs_match(fit, ...)$estimate)
  }

  if (estimator %in% c("ipw", "normalized", "aipw")) {
    return(fpgs_weight(fit, type = estimator)$estimate)
  }

  if (estimator == "stratification") {
    return(fpgs_stratify(fit, ...)$estimate)
  }
  stop("Unknown estimator.", call. = FALSE)
}


#' Print method for FPGS bootstrap inference
#'
#' @param x An object of class "fpgs_bootstrap".
#' @param ... Additional arguments, currently ignored.
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
  cat("Bootstrap replications:", x$B_valid, "valid out of", x$B, "\n")

  invisible(x)
}


#' Summary method for FPGS bootstrap inference
#'
#' @param object An object of class "fpgs_bootstrap".
#' @param ... Additional arguments, currently ignored.
#'
#' @export
summary.fpgs_bootstrap <- function(object, ...) {

  out <- data.frame(
    estimator = object$estimator,
    fpgs_method = object$fpgs_method,
    outcome_type = object$outcome_type,
    estimate = object$estimate,
    se = object$se,
    ci_lower = object$conf.int[1],
    ci_upper = object$conf.int[2],
    conf_level = object$conf.level,
    B = object$B,
    B_valid = object$B_valid
  )

  out
}


#' Confidence interval for FPGS bootstrap inference
#'
#' @param object An object of class "fpgs_bootstrap".
#' @param parm Currently ignored.
#' @param level Confidence level. Defaults to the level stored in object.
#' @param ... Additional arguments, currently ignored.
#'
#' @export
confint.fpgs_bootstrap <- function(object,
                                   parm = NULL,
                                   level = object$conf.level,
                                   ...) {

  if (!is.numeric(level) || length(level) != 1 || level <= 0 || level >= 1) {
    stop("level must be a number between 0 and 1", call. = FALSE)
  }

  alpha <- 1 - level

  ci <- stats::quantile(
    object$boot_estimates,
    probs = c(alpha / 2, 1 - alpha / 2),
    na.rm = TRUE,
    names = FALSE
  )

  names(ci) <- c("lower", "upper")
  ci
}

#' Extract FPGS bootstrap estimate
#'
#' @param object An object of class "fpgs_bootstrap".
#' @param ... Additional arguments, currently ignored.
#'
#' @export
coef.fpgs_bootstrap <- function(object, ...) {
  object$estimate
}


#' Variance-covariance matrix for FPGS bootstrap estimate
#'
#' @param object An object of class "fpgs_bootstrap".
#' @param ... Additional arguments, currently ignored.
#'
#' @export
vcov.fpgs_bootstrap <- function(object, ...) {
  out <- matrix(object$se^2, nrow = 1, ncol = 1)
  rownames(out) <- colnames(out) <- "ATE"
  out
}