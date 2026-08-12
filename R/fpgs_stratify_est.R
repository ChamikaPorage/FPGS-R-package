#' FPGS Stratification Estimator
#'
#' Estimates the average treatment effect (ATE) by stratifying on the
#' two-dimensional estimated FPGS.
#'
#' @param fit An object of class \code{"fpgs"} returned by \code{fpgs()}.
#' @param n_bins Number of quantile-based bins used for each FPGS component.
#'   Default is 5.
#' @param min_cell Minimum recommended number of treated and untreated
#'   observations within each stratum. Default is 5.
#'
#' @return An object of class \code{"fpgs_stratify"} containing the ATE
#'   estimate, stratum assignments, and stratum-specific treatment effects.
#'
#' @export
#'
fpgs_stratify <- function(fit,
                          n_bins = 5,
                          min_cell = 5) {

  if (!inherits(fit, "fpgs")) {
    stop("fit must be an object returned by fpgs().")
  }

  dat <- fit$data
  Y <- dat[[fit$outcome]]
  Tr <- dat[[fit$treatment]]
  N <- nrow(dat)

  # Quantile cut points for each FPGS component
  q0 <- unique(
    quantile(
      fit$mu0_hat,
      probs = seq(0, 1, length.out = n_bins + 1),
      na.rm = TRUE
    )
  )

  q1 <- unique(
    quantile(
      fit$mu1_hat,
      probs = seq(0, 1, length.out = n_bins + 1),
      na.rm = TRUE
    )
  )

  # Create bins for mu0_hat and mu1_hat
  bin0 <- cut(
    fit$mu0_hat,
    breaks = q0,
    include.lowest = TRUE,
    labels = FALSE
  )

  bin1 <- cut(
    fit$mu1_hat,
    breaks = q1,
    include.lowest = TRUE,
    labels = FALSE
  )

  # Cross-classify the two FPGS components
  strata <- interaction(bin0, bin1, drop = TRUE)

  tmp <- data.frame(
    Y = Y,
    Tr = Tr,
    strata = strata
  )

  # Estimate treatment effect within each stratum
  stratum_effects <- lapply(
    split(tmp, tmp$strata),
    function(d) {

      n <- nrow(d)
      n1 <- sum(d$Tr == 1)
      n0 <- sum(d$Tr == 0)

      # A treatment effect cannot be estimated if one group is absent
      if (n1 == 0 || n0 == 0) {
        return(
          data.frame(
            n = n,
            n_treated = n1,
            n_control = n0,
            mean_treated = NA_real_,
            mean_control = NA_real_,
            te = NA_real_,
            weight = n / N,
            weighted_te = NA_real_
          )
        )
      }

      mean_treated <- mean(d$Y[d$Tr == 1])
      mean_control <- mean(d$Y[d$Tr == 0])

      te <- mean_treated - mean_control

      data.frame(
        n = n,
        n_treated = n1,
        n_control = n0,
        mean_treated = mean_treated,
        mean_control = mean_control,
        te = te,
        weight = n / N,
        weighted_te = (n / N) * te
      )
    }
  )

  stratum_effects <- do.call(
    rbind,
    stratum_effects
  )

  # Check whether every stratum contains both treatment groups
  invalid <- is.na(stratum_effects$te)

  if (any(invalid)) {
    stop(
      "Some strata contain only treated or only untreated individuals. ",
      "Reduce n_bins or inspect treatment overlap."
    )
  }

  # Warn about very small cells
  small_cells <- with(
    stratum_effects,
    n_treated < min_cell | n_control < min_cell
  )

  if (any(small_cells)) {
    warning(
      "Some strata contain fewer than ",
      min_cell,
      " treated or untreated observations. ",
      "The stratification estimate may be unstable."
    )
  }

  # Weighted average of stratum-specific treatment effects
  ate <- sum(
    stratum_effects$weighted_te
  )

  out <- list(
    estimate = ate,
    n_bins = n_bins,
    min_cell = min_cell,
    strata = strata,
    stratum_effects = stratum_effects,
    method = "FPGS stratification"
  )

  class(out) <- "fpgs_stratify"

  out
}
