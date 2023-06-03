library(TensorEconometrics)

set.seed(20230603)
X <- rand_tensor(c(100, 10, 8))
Y <- rand_tensor(c(100, 10, 8))
R <- c(5,4,3,2)
obs_dim_X <- 1
obs_dim_Y <- 1
convThresh <- 1e-04
max_iter <- 400

# Tucker Regression

xt_regression <- function () {
  
}

yt_regression <- function () {
  
}

tuck_conv <- function() {
  
}

tucker_regression <- function(Y, X, R, obs_dim_X, obs_dim_Y, convThresh = 1e-05, 
                              max_iter = 400, seed = 0) {
  if (seed > 0) set.seed(seed)
  
  # Generate initial random tensor and random CP decomposition
  init_B <- rand_tensor(c(X@modes[-obs_dim_X], Y@modes[-obs_dim_Y]))
  init_list <- list(
    matrix(rnorm(X@modes[2] * R), nrow = X@modes[2]),
    matrix(rnorm(X@modes[3] * R), nrow = X@modes[3]),
    matrix(rnorm(Y@modes[2] * R), nrow = Y@modes[2]),
    matrix(rnorm(Y@modes[3] * R), nrow = Y@modes[3])
  )
  
}