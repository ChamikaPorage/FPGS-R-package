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
    return(
      fpgs_continuous_parametric(
        data,
        outcome,
        treatment,
        covariates
      )
    )
  }
  
  if (outcome_type == "continuous" && method == "rf") {
    return(
      fpgs_rf(
        data,
        outcome,
        treatment,
        covariates,
        folds = folds,
        ...
      )
    )
  }
  
  if (outcome_type == "binary" && method == "parametric") {
    return(
      fpgs_binary(
        data,
        outcome,
        treatment,
        covariates
      )
    )
  }
  
  if (outcome_type == "binary" && method == "rf") {
    return(
      fpgs_binary_rf(
        data,
        outcome,
        treatment,
        covariates,
        folds = folds,
        ...
      )
    )
  }
}