test_that("OLS_match", {
  set.seed(111111) # for reproducibility
  var_ols <- function(A){
    dimensions_A <- dim(A)
    original_A <- A[,-1]
    lagged_A <- A[,1:(dimensions_A[2]-1)]
    ols_hat <- (original_A%*%t(lagged_A)) %*% solve(lagged_A%*%t(lagged_A))
    return(ols_hat)
  }
  
  # Define the VAR/MAR coefficients
  vec_parameters <- matrix(c(0.2, 0.1, -0.2, -0.4, 0.1, 0.2, 0.5, -0.1,
                             0.4, 0.2, 0.6, -0.3, -0.1, 0.3, -0.11, 0.2,
                             0.1, 0.2, 0.5, 0.15, -0.5, -0.1, -0.1, 0.2,
                             0.23, -0.13, -0.1, 0.4, -0.3, 0.1, -0.3,
                             -0.12, -0.3, 0.4, -0.4, 0.4), nrow = 6, ncol = 6)
  
  sim_mat <- matrix_ols_sim(100, vec_parameters, 3, 2)
  
  sim_matY <- sim_mat$matrix_data[,,2:100]  # Lose one observation for lag
  sim_matX <- sim_mat$matrix_data[,,1:99]  # lagged data
  
  
  result <- HOOLS(as.tensor(sim_matY), as.tensor(sim_matX))
  expected <- var_ols(sim_mat$vector_data)
  
  expect_equal(unfold(result, c(1,2), c(3,4))@data, expected)
  
  # Now check whether alternate obs orderings also work
  alt_Y <- aperm(sim_mat$matrix_data[,,2:100], c(2,1,3))
  alt_X <- aperm(sim_mat$matrix_data[,,1:99], c(2,1,3))
  
  result2 <- HOOLS(as.tensor(alt_Y), as.tensor(alt_X))
  expect_equal(matrix(result2@data, nrow = 6), expected)
})
