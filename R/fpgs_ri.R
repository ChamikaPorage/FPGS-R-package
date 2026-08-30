#' FPGS Regression Imputation Estimator
#'
#' Estimates the average treatment effect using regression imputation
#' based on the estimated Full Prognostic Score.
#'
#' For parametric FPGS estimation, linear regression is used in the second
#' stage. For random-forest FPGS estimation, random forest with cross-fitting
#' is used in the second stage.
#'
#' @param fit An object of class \code{"fpgs"} returned by \code{fpgs()}.
#'
#' @return An object of class \code{"fpgs_ri"} containing the ATE estimate
#'   and predicted potential outcomes.
#'
#' @export
fpgs_ri <- function(fit) {
  
  if (!inherits(fit, "fpgs")) {
    stop("fit must be an object returned by fpgs().")
  }
  
  dat <- fit$data
  Y <- dat[[fit$outcome]]
  Tr <- dat[[fit$treatment]]
  
  fpgs_dat <- data.frame(
    Y = Y,
    Tr = Tr,
    mu0_hat = fit$mu0_hat,
    mu1_hat = fit$mu1_hat
  )
  
  # Parametric second-stage regression
  if (fit$method == "parametric") {
    
    model0 <- lm(
      Y ~ mu0_hat + mu1_hat,
      data = fpgs_dat[Tr == 0, , drop = FALSE]
    )
    
    model1 <- lm(
      Y ~ mu0_hat + mu1_hat,
      data = fpgs_dat[Tr == 1, , drop = FALSE]
    )
    
    pred0 <- predict(model0, newdata = fpgs_dat)
    pred1 <- predict(model1, newdata = fpgs_dat)
    
    ate <- mean(pred1 - pred0)
    
    out <- list(
      estimate = ate,
      pred0 = pred0,
      pred1 = pred1,
      model0 = model0,
      model1 = model1,
      second_stage = "linear regression",
      method = "FPGS regression imputation"
    )
  }
  
  # Random-forest second-stage regression
  if (fit$method == "rf") {
    
    if (!requireNamespace("ranger", quietly = TRUE)) {
      stop("Package 'ranger' is required.")
    }
    
    if (!requireNamespace("caret", quietly = TRUE)) {
      stop("Package 'caret' is required.")
    }
    
    folds <- fit$folds
    n <- nrow(fpgs_dat)
    
    fold_id <- caret::createFolds(
      Y,
      k = folds,
      list = TRUE
    )
    
    pred0 <- rep(NA_real_, n)
    pred1 <- rep(NA_real_, n)
    
    for (k in seq_len(folds)) {
      
      train_index <- unlist(fold_id[-k])
      test_index <- fold_id[[k]]
      
      train_data <- fpgs_dat[train_index, , drop = FALSE]
      test_data <- fpgs_dat[test_index, , drop = FALSE]
      
      train0 <- train_data[
        train_data$Tr == 0,
        ,
        drop = FALSE
      ]
      
      train1 <- train_data[
        train_data$Tr == 1,
        ,
        drop = FALSE
      ]
      
      model0 <- ranger::ranger(
        Y ~ mu0_hat + mu1_hat,
        data = train0,
        num.trees = 300,
        mtry = 2
      )
      
      model1 <- ranger::ranger(
        Y ~ mu0_hat + mu1_hat,
        data = train1,
        num.trees = 300,
        mtry = 2
      )
      
      pred0[test_index] <- predict(
        model0,
        data = test_data
      )$predictions
      
      pred1[test_index] <- predict(
        model1,
        data = test_data
      )$predictions
    }
    
    ate <- mean(pred1 - pred0)
    
    out <- list(
      estimate = ate,
      pred0 = pred0,
      pred1 = pred1,
      folds = folds,
      second_stage = "random forest",
      method = "FPGS regression imputation"
    )
  }
  
  class(out) <- "fpgs_ri"
  
  out
}
