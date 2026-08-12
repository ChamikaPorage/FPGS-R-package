#' Estimate Binary FPGS using Random Forest with Cross-Fitting
#'
#' @param data Data frame
#' @param outcome Binary outcome variable name
#' @param treatment Treatment variable name
#' @param covariates Covariate names
#' @param folds Number of folds
#'
#' @export
fpgs_binary_rf <- function(data,
                           outcome,
                           treatment,
                           covariates = NULL,
                           folds = 5,
                           num.trees = 300,
                           mtry = NULL,
                           min.node.size = 5) {
  
  if (!requireNamespace("ranger", quietly = TRUE)) {
    stop("Package 'ranger' is required.")
  }
  
  if (!requireNamespace("caret", quietly = TRUE)) {
    stop("Package 'caret' is required.")
  }
  
  dat <- as.data.frame(data)
  
  if (is.null(covariates)) {
    covariates <- setdiff(names(dat), c(outcome, treatment))
  }
  
  if (is.null(mtry)) {
    mtry <- max(1, floor(sqrt(length(covariates))))
  }
  
  Y <- dat[[outcome]]
  Tr <- dat[[treatment]]
  n <- nrow(dat)
  
  dat[[outcome]] <- factor(dat[[outcome]], levels = c(0, 1))
  
  fold_id <- caret::createFolds(Y, k = folds, list = TRUE)
  
  mu0_hat <- rep(NA_real_, n)
  mu1_hat <- rep(NA_real_, n)
  
  model_data <- dat[, c(outcome, treatment, covariates), drop = FALSE]
  
  form <- as.formula(
    paste(outcome, "~", paste(covariates, collapse = " + "))
  )
  
  for (k in seq_len(folds)) {
    
    train_idx <- unlist(fold_id[-k])
    valid_idx <- fold_id[[k]]
    
    train_data <- model_data[train_idx, , drop = FALSE]
    valid_data <- model_data[valid_idx, , drop = FALSE]
    
    train_data_0 <- train_data[train_data[[treatment]] == 0, , drop = FALSE]
    train_data_1 <- train_data[train_data[[treatment]] == 1, , drop = FALSE]
    
    rf0 <- ranger::ranger(
      formula = form,
      data = train_data_0,
      probability = TRUE,
      num.trees = num.trees,
      mtry = mtry,
      min.node.size = min.node.size
    )
    
    rf1 <- ranger::ranger(
      formula = form,
      data = train_data_1,
      probability = TRUE,
      num.trees = num.trees,
      mtry = mtry,
      min.node.size = min.node.size
    )
    
    pred0 <- predict(rf0, data = valid_data)$predictions
    pred1 <- predict(rf1, data = valid_data)$predictions
    
    mu0_hat[valid_idx] <- pred0[, "1"]
    mu1_hat[valid_idx] <- pred1[, "1"]
  }
  
  out <- list(
    data = data,
    outcome = outcome,
    treatment = treatment,
    covariates = covariates,
    mu0_hat = mu0_hat,
    mu1_hat = mu1_hat,
    fpgs = data.frame(
      mu0_hat = mu0_hat,
      mu1_hat = mu1_hat
    ),
    outcome_type = "binary",
    method = "rf",
    learner = "random forest",
    crossfit = TRUE,
    folds = folds,
    mtry = mtry
  )
  
  class(out) <- "fpgs"
  
  out
}