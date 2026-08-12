###Data generation process 1####

gen_data_1 <- function(N, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  X1 <- rnorm(N)
  X2 <- rnorm(N)
  X3 <- rnorm(N)
  X4 <- rnorm(N)
  X5 <- rbinom(N, 1, 0.5)
  X6 <- rnorm(N)
  
  E1 <- rbinom(N, 1, 0.5)
  E2 <- rbinom(N, 1, 0.5)
  E3 <- rbinom(N, 1, 0.5)
  
  p_t <- plogis(0.1 * X1 - 0.1 * X2 + 1.1 * X3 - 1.1 * X4 + 0.4 * X5)
  Tr <- rbinom(N, 1, p_t)
  
  mu <- -3.85 + 0.5 * X1 - 2 * X2 - 0.5 * X3 + 2 * X4 +
    X6 - E1 - 2 * E3 + 5 * Tr + Tr * E1 + 4 * Tr * E2 - 4 * Tr * E3
  
  Y <- rnorm(N, mu, 1)
  
  data.frame(Y, Tr, X1, X2, X3, X4, X5, X6, E1, E2, E3)
}
