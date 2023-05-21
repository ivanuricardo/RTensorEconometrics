test_that("C1_regression", {
  library(TensorEconometrics)
  set.seed(20230521)
  Y <- rand_tensor(c(161, 32, 3))
  X <- rand_tensor(c(161, 32, 3))
  obs_dim_X <- 1
  obs_dim_Y <- 1
  R <- 4
  
  cp_regression(Y, X, R, obs_dim_X, obs_dim_Y)
  
  
  
})