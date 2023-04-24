# Preset Parameter tensor and matrix for OLS

test_that("OLS_match", {
  vec_parameters <- matrix(c(0.2, 0.1, -0.2, -0.4, 0.1, 0.2, 0.5, -0.1,
                           0.4, 0.2, 0.6, -0.3, -0.1, 0.3, -0.11, 0.2,
                           0.1, 0.2, 0.5, 0.15, -0.5, -0.1, -0.1, 0.2,
                           0.23, -0.13, -0.1, 0.4, -0.3, 0.1, -0.3,
                           -0.12, -0.3, 0.4, -0.4, 0.4), nrow = 6, ncol = 6)
  mat_parameters <- array(vec_parameters, dim = c(3,2,3,2))
   
  
})