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
      matrix(rnorm(X@modes[2] * R[1]), nrow = X@modes[2]),
      matrix(rnorm(X@modes[3] * R[2]), nrow = X@modes[3]),
      matrix(rnorm(Y@modes[2] * R[3]), nrow = Y@modes[2]),
      matrix(rnorm(Y@modes[3] * R[4]), nrow = Y@modes[3])
  )
  
  converged <- FALSE
  num_iter <- 0
  while (num_iter < max_iter) {
    num_iter <- num_iter + 1
    
    for (dim in 1:init_B@num_modes) {
      if (dim < (init_B@num_modes/2 + 1)) {
        x_reg_list <- xt_regression(init_list = init_list, Y = Y, X = X,
                                   R = R, idx = dim)
        pre_init_list <- init_list
        init_list[[dim]] <- x_reg_list$B1
        lambdas <- x_reg_list$x_lambdas
        if (convergence_func(pre_init_list, init_list, dim, convThresh)) {
          converged <- TRUE
          break
        }
      } else {
        y_reg_list <- yt_regression(init_list = init_list, Y = Y, X = X,
                                         R = R, Ddims= Ddims, idx = dim)
        pre_init_list <- init_list
        init_list[[dim]] <- y_reg_list$B3
        lambdas <- y_reg_list$y_lambdas
        if (convergence_func(pre_init_list, init_list, dim, convThresh)) {
          converged <- TRUE
          break
        }
      }
    }
  }
}