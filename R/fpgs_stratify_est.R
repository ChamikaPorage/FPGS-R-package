fpgs_stratify <- function(fit, n_strata = 5) {

  if (!inherits(fit, "fpgs")) {
    stop("fit must be an object returned by fpgs()")
  }

  dat <- fit$data
  Y <- dat[[fit$outcome]]
  Tr <- dat[[fit$treatment]]
  N <- nrow(dat)

  q0 <- unique(quantile(
    fit$mu0_hat,
    probs = seq(0, 1, length.out = n_strata + 1),
    na.rm = TRUE
  ))

  q1 <- unique(quantile(
    fit$mu1_hat,
    probs = seq(0, 1, length.out = n_strata + 1),
    na.rm = TRUE
  ))

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

  strata <- interaction(bin0, bin1, drop = TRUE)

  tmp <- data.frame(
    Y = Y,
    Tr = Tr,
    strata = strata
  )

  stratum_effects <- lapply(split(tmp, tmp$strata), function(d) {

    n <- nrow(d)

    if (sum(d$Tr == 1) == 0 || sum(d$Tr == 0) == 0) {
      return(data.frame(
        n = n,
        mean_treated = NA_real_,
        mean_control = NA_real_,
        te = NA_real_,
        weight = n / N,
        weighted_te = NA_real_
      ))
    }

    mean_treated <- mean(d$Y[d$Tr == 1])
    mean_control <- mean(d$Y[d$Tr == 0])
    te <- mean_treated - mean_control

    data.frame(
      n = n,
      mean_treated = mean_treated,
      mean_control = mean_control,
      te = te,
      weight = n / N,
      weighted_te = (n / N) * te
    )
  })

  stratum_effects <- do.call(rbind, stratum_effects)

  ate <- sum(stratum_effects$weighted_te, na.rm = TRUE)

  out <- list(
    estimate = ate,
    n_strata = n_strata,
    strata = strata,
    stratum_effects = stratum_effects,
    method = "FPGS stratification"
  )

  class(out) <- "fpgs_stratify"

  out
}
