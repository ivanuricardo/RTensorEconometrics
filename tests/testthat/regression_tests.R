source("R/tensor_regressions.R")
library(testthat)

test_that("Tensor inverse returns identity", {
  expect_equal(tensor_inverse())
})

####### Testing Data
# I generate tensor data from a TAR model with the following matrix of parameters
library(tsDyn)
set.seed(20230412)
mat_par <- matrix(c(0.2, 0.1, -0.2, -0.4, 0.1, 0.2, 0.5, -0.1, 0.4, 0.2, 0.6,
                    -0.3, -0.1, 0.3, -0.11, 0.2, 0.1, 0.2, 0.5, 0.15, -0.5, -0.1,
                    -0.1, 0.2, 0.23, -0.13, -0.1, 0.4, -0.3, 0.1, -0.3, -0.12, -0.3,
                    0.4, -0.4, 0.4), nrow = 6, ncol = 6)
arr_par <- array(mat_par, dim = c(3,2,3,2))
sim_data <- VAR.sim(mat_par, n = 1000, include = 'none')
mat_data <- array(t(sim_data), dim = c(3,2,1000))

# Suppose I take the lagged version of mat_data as predictor and take the response
# as mat_data with the first observation removed.

response_data <- mat_data[,,-1]
predictor_data <- mat_data[,,1:999]

# Apply least squares
mat_parameters <- HOOLS(response_data, predictor_data)
mat_parameters
arr_par
