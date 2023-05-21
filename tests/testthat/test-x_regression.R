test_that("x_regression_output", {
  set.seed(20230521)
  
  least_squares <- function(X, y) {
    return(solve(crossprod(X)) %*% crossprod(X, y))
  }
  
  # num of matrices in the list is num of dimensions in B
  # num of columns in each matrix is R = 4
  # num of rows is the dimension value for X
  Y <- rand_tensor(c(100, 7, 3))
  X <- rand_tensor(c(100, 7, 3))
  R <- 6
  
  init_list <- list(
    matrix(rnorm(X@modes[2] * R), nrow = X@modes[2]),
    matrix(rnorm(X@modes[3] * R), nrow = X@modes[3]),
    matrix(rnorm(Y@modes[2] * R), nrow = Y@modes[2]),
    matrix(rnorm(Y@modes[3] * R), nrow = Y@modes[3])
  )
  
  # Find the result for the regression
  x_regression(init_list, Y, X, R, 1)
  
  # Should be the same as if I ran the regression 
})

test_that("x_regression_cp_list", {
  
})