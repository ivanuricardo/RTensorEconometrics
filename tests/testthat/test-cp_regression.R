test_that("cp_regression_identification", {
  skip(message = "Skipping, remove for official test")
  # Initialize data tensors and R
  set.seed(20230527)
  X <- rand_tensor(c(100, 4, 3))
  Y <- rand_tensor(c(100, 4, 3))
  R <- 6
  obs_dim_X <- 1
  obs_dim_Y <- 1
  
  # Estimate CP regression
  cp_est <- cp_regression(Y, X, R, 1, 1)
  
  # Ensure length and ordering identification is met
  # First, lengths of vectors are one in line with Kolda and Bader 2008
  norm_vec <- c()
  for (r in 1:R) {
    factor_columns <- lapply(cp_est$factor_mat, function(n) n[,r])
    factor_norms <- sapply(factor_columns, function(x) norm(x, type = "2"))
    norm_vec <- append(norm_vec, factor_norms)
  }
  expect_equal(norm_vec, rep(1, 24))
  
  # Next, ensure ordering is in line Kolda and Bader 2008
  
  expect_true(all(diff(cp_est$lambdas) <= 0))
})

test_that("reconstruct_cp same as output B", {
  set.seed(20230614)
  X <- rnorm_tnsr(c(100,2,3))
  Y <- rnorm_tnsr(c(100,2,3))
  R <- 3
  obs_dim_X <- 1
  obs_dim_Y <- 1
  
  cp_est <- cp_regression(Y,X,R,obs_dim_X, obs_dim_Y)
  reconstructed <- reconstruct_cp(cp_est$factor_mat[[1]], cp_est$factor_mat[[2]],
                                  cp_est$factor_mat[[3]], cp_est$factor_mat[[4]],
                                  r=R, lambda = cp_est$lambdas)
  expect_equal(as.tensor(reconstructed), cp_est$B)
})