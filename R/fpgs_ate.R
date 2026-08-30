#' Estimate the ATE Using FPGS Methods
#'
#' Runs regression imputation (RI), matching, inverse probability weighting
#' (IPW), normalized inverse probability weighting (NIPW), augmented inverse
#' probability weighting (AIPW), and stratification estimators using an
#' estimated FPGS.
#'
#' @param fit An object of class \code{"fpgs"} returned by \code{fpgs()}.
#' @param inference Character string specifying the type of inference.
#'   Options are \code{"none"}, \code{"analytic"}, or \code{"bootstrap"}.
#'   Default is \code{"none"}.
#' @param B Number of bootstrap replications when
#'   \code{inference = "bootstrap"}. Default is 200.
#' @param conf.level Confidence level for confidence intervals.
#'   Default is 0.95.
#' @param M Number of matches used by \code{fpgs_match()}. Default is 1.
#' @param n_bins Number of quantile-based bins used for each FPGS component
#'   in \code{fpgs_stratify()}. Default is 4.
#' @param min_cell Minimum recommended number of treated and untreated
#'   observations within each retained stratum. Default is 3.
#'
#' @return A data frame containing treatment effect estimates. When analytic
#'   or bootstrap inference is requested, estimated variances, standard errors,
#'   and confidence intervals are also returned.
#'
#' @details
#' If \code{inference = "none"}, only point estimates are returned.
#'
#' If \code{inference = "analytic"}, analytic inference is used. This option
#' is available only for parametrically estimated FPGS models. Regression
#' imputation and weighting estimators use sandwich variance estimation,
#' matching uses the standard error returned by \code{Matching::Match()},
#' and stratification uses an analytic within-stratum variance estimator.
#'
#' If \code{inference = "bootstrap"}, bootstrap standard errors and percentile
#' confidence intervals are calculated by resampling observations, refitting
#' the FPGS, and recomputing each treatment effect estimator. The default
#' number of bootstrap replications is 200 and can be changed using
#' the \code{B} argument.
#'
#' @export
fpgs_ate <- function(fit,
                     inference = c("none", "analytic", "bootstrap"),
                     B = 200,
                     conf.level = 0.95,
                     M = 1,
                     n_bins = 4,
                     min_cell = 3) {

  if (!inherits(fit, "fpgs")) {
    stop("fit must be an object returned by fpgs()")
  }

  inference <- match.arg(inference)

  if (inference == "analytic" && fit$method != "parametric") {
    stop(
      "Analytic inference is only available for parametric FPGS fits. ",
      "Use inference = 'bootstrap' for random forest fits."
    )
  }

  if (!is.numeric(conf.level) || length(conf.level) != 1 ||
      conf.level <= 0 || conf.level >= 1) {
    stop("conf.level must be a number between 0 and 1")
  }

  # Run the six treatment effect estimators
  ri <- fpgs_ri(fit)

  matching <- fpgs_match(fit, M = M)

  ipw <- fpgs_weight(fit, type = "ipw")

  nipw <- fpgs_weight(fit, type = "nipw")

  aipw <- fpgs_weight(fit, type = "aipw")

  stratification <- fpgs_stratify(fit, n_bins = n_bins, min_cell = min_cell)

  estimator_names <- c("RI", "Matching", "IPW", "NIPW", "AIPW", "Stratification")

  estimates <- c(
    ri$estimate,
    matching$estimate,
    ipw$estimate,
    nipw$estimate,
    aipw$estimate,
    stratification$estimate
  )

  # Point estimates only
  if (inference == "none") {

    return(
      data.frame(
        estimator = estimator_names,
        estimate = estimates
      )
    )
  }

  # Analytic inference
  if (inference == "analytic") {

    ri_var <- ri_variance(fit)

    ipw_var <- weight_variance(fit, type = "ipw")

    nipw_var <- weight_variance(fit, type = "nipw")

    aipw_var <- weight_variance(fit, type = "aipw")

    variances <- c(
      ri_var$variance,
      matching$variance,
      ipw_var$variance,
      nipw_var$variance,
      aipw_var$variance,
      stratification$variance
    )

    standard_errors <- c(
      ri_var$se,
      matching$se,
      ipw_var$se,
      nipw_var$se,
      aipw_var$se,
      stratification$se
    )

    z_value <- stats::qnorm(
      1 - (1 - conf.level) / 2
    )

    ci_lower <- estimates - z_value * standard_errors
    ci_upper <- estimates + z_value * standard_errors

    return(
      data.frame(
        estimator = estimator_names,
        estimate = estimates,
        variance = variances,
        se = standard_errors,
        ci_lower = ci_lower,
        ci_upper = ci_upper
      )
    )
  }

  # Bootstrap inference
  bootstrap_results <- list(

    RI = fpgs_bootstrap(
      fit,
      estimator = "ri",
      B = B,
      conf.level = conf.level
    ),

    Matching = fpgs_bootstrap(
      fit,
      estimator = "matching",
      B = B,
      conf.level = conf.level,
      M = M
    ),

    IPW = fpgs_bootstrap(
      fit,
      estimator = "ipw",
      B = B,
      conf.level = conf.level
    ),

    NIPW = fpgs_bootstrap(
      fit,
      estimator = "nipw",
      B = B,
      conf.level = conf.level
    ),

    AIPW = fpgs_bootstrap(
      fit,
      estimator = "aipw",
      B = B,
      conf.level = conf.level
    ),

    Stratification = fpgs_bootstrap(
      fit,
      estimator = "stratification",
      B = B,
      conf.level = conf.level,
      n_bins = n_bins,
      min_cell = min_cell
    )
  )

  do.call(
    rbind,
    lapply(names(bootstrap_results), function(name) {

      result <- bootstrap_results[[name]]

      data.frame(
        estimator = name,
        estimate = result$estimate,
        variance = result$variance,
        se = result$se,
        ci_lower = result$conf.int[1],
        ci_upper = result$conf.int[2]
      )
    })
  )
}
