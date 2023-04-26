test_that("MatrixOLS-Sim", {
  parameter_matrix <- matrix(runif(36, min = -1, max = 1), nrow = 6)
  obs <- 1000
  
  result <- matrix_ols_sim(obs = obs, parameter_matrix = parameter_matrix, num_rows = 3, num_cols = 2)
  dim_results <- dim(result$matrix_data)
  
  expect_equal(dim_results, c(1000, 3, 2))
})