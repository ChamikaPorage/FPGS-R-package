gen_data_bin <- function(N,
                         seed = NULL,
                         alpha0 = -1.5,
                         beta0 = -2,
                         gamma = log(1.5),
                         delta = log(2)) {
  if (!is.null(seed)) set.seed(seed)
  
  Alpha_vec <- c(
    log(1.5), log(2), log(2.5),
    log(1.5), log(2), log(2.5),
    0, 0, 0
  )
  
  Beta_vec <- c(
    log(1.5), log(2), log(2.5),
    0, 0, 0,
    log(1.5), log(2), log(2.5)
  )
  
  X <- matrix(rnorm(N * 9), ncol = 9)
  colnames(X) <- paste0("x", 1:9)
  
  ps_true <- plogis(alpha0 + drop(X %*% Alpha_vec))
  Tr <- rbinom(N, 1, ps_true)
  
  x2 <- X[, "x2"]
  lp0 <- beta0 + drop(X %*% Beta_vec)
  pY0 <- plogis(lp0)
  pY1 <- plogis(lp0 + gamma + delta * x2)
  
  Y <- rbinom(N, 1, ifelse(Tr == 1, pY1, pY0))
  
  data.frame(
    Y = Y,
    Tr = Tr,
    X,
    check.names = FALSE
  )
}
