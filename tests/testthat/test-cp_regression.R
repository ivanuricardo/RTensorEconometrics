test_that("cp_regression_identification", {
  # Initialize data tensors and R
  set.seed(20230527)
  X <- rand_tensor(c(100, 4, 3))
  Y <- rand_tensor(c(100, 4, 3))
  R <- 6
  
  # Estimate CP regression
  cp_est <- cp_regression(Y, X, R, 1, 1)
  
  # Ensure length and ordering identification is met
  # First, lengths of vectors are one in line with Lock 2018
  norm_vec <- c()
  for (r in 1:R) {
    factor_columns <- lapply(cp_est$factor_mat, function(n) n[,r])
    factor_norms <- sapply(factor_columns, function(x) norm(x, type = "2"))
    norm_vec <- append(norm_vec, factor_norms)
  }
  expect_equal(norm_vec, rep(1, 24))
  
  # Next, ensure ordering is in line with Lock 2018
  u1 <- norm_vec[1:6]
  u2 <- norm_vec[7:12]
  v1 <- norm_vec[13:18]
  v2 <- norm_vec[19:24]
  
  expect_true(all(diff(u1) <= 0))
})