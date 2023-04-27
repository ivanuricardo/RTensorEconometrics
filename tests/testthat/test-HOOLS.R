test_that("OLS_match", {
  set.seed(111111) # for reproducibility
  var_ols <- function(A){
  dimensions_A <- dim(A)
  original_A <- A[-1,]
  lagged_A <- A[1:(dimensions_A[1]-1),]
  ols_hat <- solve(t(lagged_A)%*%lagged_A) %*% (t(lagged_A)%*%original_A)
  return(t(ols_hat))
  }
  library(rTensor)
  
  # Define the VAR/MAR coefficients
  vec_parameters <- matrix(c(0.2, 0.1, -0.2, -0.4, 0.1, 0.2, 0.5, -0.1,
                             0.4, 0.2, 0.6, -0.3, -0.1, 0.3, -0.11, 0.2,
                             0.1, 0.2, 0.5, 0.15, -0.5, -0.1, -0.1, 0.2,
                             0.23, -0.13, -0.1, 0.4, -0.3, 0.1, -0.3,
                             -0.12, -0.3, 0.4, -0.4, 0.4), nrow = 6, ncol = 6)
  
  sim_mat <- matrix_ols_sim(100, vec_parameters, 3, 2)
  
  sim_matY <- sim_mat$matrix_data[2:100,,]  # Lose one observation for lag
  sim_matX <- sim_mat$matrix_data[1:99,,]  # lagged data
  
  
  result <- HOOLS(as.tensor(sim_matY), as.tensor(sim_matX), obs_dim_Y = 1, obs_dim_X = 1)
  expected <- var_ols(sim_mat$vector_data)
  
  expect_equal(t(matrix(result@data, nrow = 6)), expected)
  
  # Now check whether alternate obs orderings also work
  alt_Y <- aperm(sim_mat$matrix_data[2:100,,], c(2,3,1))
  alt_X <- aperm(sim_mat$matrix_data[1:99,,], c(2,3,1))
  
  result2 <- HOOLS(as.tensor(alt_Y), as.tensor(alt_X), obs_dim_Y = 3, obs_dim_X = 3)
  expect_equal(t(matrix(result2@data, nrow = 6)), expected)
})

test_that("Higher_Order_tensors", {
  set.seed(12345678) # for reproducibility
  library(bigtime)
  library(rTensor)
  var_ols <- function(A){
    dimensions_A <- dim(A)
    original_A <- A[-1,]
    lagged_A <- A[1:(dimensions_A[1]-1),]
    ols_hat <- solve(t(lagged_A)%*%lagged_A) %*% (t(lagged_A)%*%original_A)
    return(t(ols_hat))
  }
  
  # Define the VAR/MAR coefficients
  sim_data <- simVAR(periods = 100, k = 24, p = 1, sparsity_pattern = "none", seed = 12345678)
 
  sim_array <- array(sim_data$Y, dim = c(100, 3, 2, 4))
  sim_tensY <- sim_array[2:100,,,]
  sim_tensX <- sim_array[1:99,,,]
 
  result <- HOOLS(as.tensor(sim_tensY), as.tensor(sim_tensX), obs_dim_Y = 1, obs_dim_X = 1)
  expected <- var_ols(sim_data$Y)
  expect_equal(t(matrix(result@data, nrow = 24, dimnames = dimnames(expected))), expected)
})

