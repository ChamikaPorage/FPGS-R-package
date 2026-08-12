#' Run All FPGS Estimators
#'
#' @param fit An object returned by fpgs()
#'
#' @return Data frame containing ATE estimates
#'
#' @export

fpgs_all <- function(fit,
                     inference = FALSE,
                     B = 500,
                     conf.level = 0.95,
                     ...) {

  if (!inference) {


    ri <- fpgs_ri(fit)
    mat <- fpgs_match(fit)
    ipw <- fpgs_weight(fit, type = "ipw")
    nipw <- fpgs_weight(fit, type = "normalized")
    aipw <- fpgs_weight(fit, type = "aipw")
    strat <- fpgs_stratify(fit)

    return(
      data.frame(
        estimator = c(
          "RI",
          "Matching",
          "IPW",
          "NIPW",
          "AIPW",
          "Stratification"
        ),
        estimate = c(
          ri$estimate,
          mat$estimate,
          ipw$estimate,
          nipw$estimate,
          aipw$estimate,
          strat$estimate
        )
      )
    )


  }

  res_list <- list(


    RI = fpgs_bootstrap(
      fit,
      estimator = "ri",
      B = B,
      conf.level = conf.level,
      ...
    ),

    Matching = fpgs_bootstrap(
      fit,
      estimator = "matching",
      B = B,
      conf.level = conf.level,
      ...
    ),

    IPW = fpgs_bootstrap(
      fit,
      estimator = "ipw",
      B = B,
      conf.level = conf.level
    ),

    NIPW = fpgs_bootstrap(
      fit,
      estimator = "normalized",
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
      ...
    )


  )

  do.call(
    rbind,
    lapply(names(res_list), function(nm) {


      x <- res_list[[nm]]

      data.frame(
        estimator = nm,
        estimate = x$estimate,
        se = x$se,
        ci_lower = x$conf.int[1],
        ci_upper = x$conf.int[2]
      )
    })


  )
}
