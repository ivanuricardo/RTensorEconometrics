# Tensor Regressions

#' Simulate data from a VAR process and estimate the OLS coefficients
#'
#' This function simulates data from a VAR process with coefficients given by a parameter matrix and then
#' estimates the OLS coefficients using matrix algebra.
#' 
#' @param obs An integer specifying the number of observations to simulate.
#' @param parameter_matrix A num_rows * num_cols by num_rows * num_cols parameter matrix for the VAR process.
#' @param num_rows An integer specifying the number of rows in the parameter matrix.
#' @param num_cols An integer specifying the number of columns in the parameter matrix.
#'
#' @return A list containing the simulated data as a 3-mode tensor with dimensions obs x num_rows x num_cols and
#' as a matrix with dimensions obs x (num_rows * num_cols).
#' 
#' @details This function simulates data from a VAR process of the form X_t = A_1 X_{t-1} + ... + A_p X_{t-p} + u_t,
#' where X_t is a num_rows x num_cols matrix of observed variables, A_1, ..., A_p are num_rows x num_cols parameter
#' matrices, p is the order of the VAR process, and u_t is a num_rows x num_cols matrix of errors with independent,
#' identically distributed normal entries. The parameter matrix is given by vec(A_1), vec(A_2), ..., vec(A_p),
#' where vec(A_i) denotes the vectorization of the matrix A_i. The OLS coefficients for the VAR process are then
#' estimated using matrix algebra.
#' 
#' @examples
#' # Simulate data from a VAR process with 2 variables and order 1
#' parameter_matrix <- matrix(c(0.5, 0.3, 0.2, 0.4), nrow = 2)
#' matrix_ols_sim(obs = 100, parameter_matrix = parameter_matrix, num_rows = 2, num_cols = 1)
#'
#' @importFrom stats rnorm
#' @export
#' @rdname matrix_ols_sim
matrix_ols_sim <- function(obs, parameter_matrix, num_rows, num_cols) {
  mat_parameters <- aperm(array(parameter_matrix, dim = c(num_rows,num_cols,num_rows,num_cols)), c(1,3,2,4))
  row_col_prod <- num_rows*num_cols
  
  # Generate random normal variables for the initial values
  init_vals <- matrix(rnorm(row_col_prod), nrow = 1)
  
  # Simulate the VAR process
  sim_data <- matrix(0, nrow = obs, ncol = row_col_prod)
  for (i in 2:obs) {
    sim_data[i,] <- parameter_matrix %*% sim_data[i-1,] + rnorm(6)
  }
  sim_data[1, ] <- init_vals
  sim_mat <- array(sim_data, dim = c(obs,3,2))
  return(list(matrix_data = sim_mat, vector_data = sim_data))
}

#' Compute the Higher-Order Ordinary Least Squares (HOOLS) estimator
#'
#' The HOOLS estimator is a generalization of the ordinary least squares (OLS) estimator to higher-order tensors.
#' 
#' @param Y A tensor of response variables.
#' @param X A tensor of predictor variables.
#' @param obs_dim_Y An integer indicating the dimension of \code{Y} along which observations are stacked.
#' @param obs_dim_X An integer indicating the dimension of \code{X} along which observations are stacked.
#'
#' @return A tensor of HOOLS estimators.
#' 
#' @details This function computes the HOOLS estimator for a tensor regression model. The HOOLS estimator is defined
#' as the solution to the equation X = X_1 B_1 + X_2 B_2 + ... + X_r B_r, where X is a tensor of predictor variables,
#' X_1, X_2, ..., X_r are subtensors of X obtained by fixing the values of certain indices, and B_1, B_2, ..., B_r
#' are tensors of HOOLS estimators.
#' 
#' @examples
#' # Generate random tensor data
#' library(rTensor)
#' Y <- rand_tensor(c(100, 3,4))
#' X <- rand_tensor(c(100, 3,2))
#' 
#' # Compute HOOLS estimator with observations along the first dimension
#' HOOLS(Y, X, obs_dim_Y = 1, obs_dim_X = 1)
#'
#' @importFrom ttt tensor_inverse
#' @export
#' @rdname HOOLS
HOOLS <- function(Y, X, obs_dim_Y = length(dim(Y)), obs_dim_X = length(dim(X))) {
  # Check input types
  if (!inherits(Y, "Tensor") || !inherits(X, "Tensor")) {
    stop("Y and X must be tensors")
  }
  
  # Compute numerator and denominator tensors
  numerator <- ttt(X, Y, obs_dim_Y, obs_dim_X)
  denominator <- ttt(X, X, obs_dim_X, obs_dim_X)
  
  # Invert the denominator tensor
  inverted_den <- tensor_inverse(denominator)
  
  # Compute HOOLS estimator
  number_dims <- length(dim(numerator))
  ols_hat <- ttt(inverted_den, numerator, (number_dims/2 + 1):number_dims, 1:(number_dims/2))
  
  return(ols_hat)
}

#' Regression associated with X
#'
#' This function computes the regression for the X portion of the CP 
#' regression. This is part one of a two part scheme for the parameters
#' of the CP regression. 
#'
#' @param init_list The initial list of CP decomposition factors.
#' @param Y The target tensor.
#' @param X The input tensor.
#' @param R The rank of the CP decomposition.
#' @param idx The index associated with the X regression
#'
#' @return The factor matrix associated with the X regression
#' 
#' @seealso
#' \code{\link{y_regression}}, \code{\link{conv_cond}}
#' \code{\link{cp_regression}}
#'
#' @export
x_regression <- function(init_list, Y, X, R, idx) {
  X_reg <- matrix(nrow = prod(Y@modes), ncol = 0)
  
  for (r in 1:R) {
    reduced_CP_list <- init_list[-idx]
    factor_column <- lapply(reduced_CP_list, function(m) m[, r])
    outer_product <- Reduce(function(x, y) x %o% y, factor_column)
    Cr <- ttt(X, as.tensor(outer_product), alongA = (Y@num_modes+1-idx),
              alongB = 1)
    unfolded_Cr <- unfold(Cr, c(1, 3, 4), 2)
    X_reg <- cbind(X_reg, unfolded_Cr@data)
  }
  
  vec_B1 <- solve(crossprod(X_reg)) %*% (crossprod(X_reg, vec(Y)))
  B1 <- matrix(vec_B1, ncol = R)
  
  # Extract lambdas according to euclidean normalization
  x_lambdas <- apply(B1, 2, function(x) sqrt(sum(x^2)))
  B1_norm <- sweep(B1, 2, x_lambdas, "/")
  
  # Sort matrix in decreasing order
  sorted_idx <- order(x_lambdas, decreasing = TRUE)
  B1_permuted <- B1_norm[, sorted_idx]
  
  return(list(B1 = B1_permuted, x_lambdas = x_lambdas[sorted_idx]))
}

#' Regression associated with Y
#'
#' This function computes the regression for the Y portion of the CP 
#' regression. This is part two of a two part scheme for the parameters
#' of the CP regression. 
#'
#' @param init_list The initial list of CP decomposition factors.
#' @param Y The response tensor
#' @param X The predictor tensor
#' @param R The CP rank
#' @param Ddims Vector with dimensions of D
#' @param idx Index for Y tensor associated with CP regression
#'
#' @return The D1 regression factor matrix.
#'
#' @seealso
#' \code{\link{x_regression}}, \code{\link{conv_cond}}
#' \code{\link{cp_regression}}
#'
#' @export
y_regression <- function(init_list, Y, X, R, Ddims, idx) {
  D1 <- matrix(nrow = Ddims[(idx-2)], ncol = 0)
  y_idx <- idx - 1
  norm_cols <- function(mat) norm(mat, type = "2")
  for (r in 1:R) {
    reduced_CP_list <- init_list[-idx]
    factor_column <- lapply(reduced_CP_list, function(n) n[, r])
    outer_product <- Reduce(function(x, y) x %o% y, factor_column)
    Dr <- ttt(X, as.tensor(outer_product), alongA = 2:3, alongB = 1:2)
    D1 <- cbind(D1, vec(Dr))
  }
  
  # Unfold Y along the dimensions not considered
  Y_unfolded <- unfold(Y, row_idx = setdiff(1:3, y_idx),
                       col_idx = y_idx)@data
  
  # Compute the estimate for B^(3) through OLS
  B3 <- t(solve(crossprod(D1)) %*% t(D1) %*% Y_unfolded)
  
  # Extract lambdas to normalize via Euclidean norm
  y_lambdas <- apply(B3, 2, norm_cols)
  B3_norm <- sweep(B3, 2, y_lambdas, "/")
  
  # Permute the matrix via norms
  sorted_idx <- order(y_lambdas, decreasing = TRUE)
  B3_permuted <- B3_norm[, sorted_idx]
  
  return(list(B3 = B3_permuted, y_lambdas = y_lambdas[sorted_idx]))
}

convergence_func <- function(pre_init_list, init_list, idx, convThresh) {
  fnorm_conv <- norm(pre_init_list[[idx]] - init_list[[idx]])
  if(fnorm_conv < convThresh) return(TRUE)
  return(FALSE)
}

#' Perform CP regression
#'
#' This function performs CP regression to estimate the tensor B given
#' the input tensor Y and covariate tensor X.
#'
#' @param Y The input tensor.
#' @param X The covariate tensor.
#' @param R The rank of the CP decomposition.
#' @param obs_dim_X The observed dimensions of the covariate tensor X.
#' @param obs_dim_Y The observed dimensions of the input tensor Y.
#' @param convThresh The convergence threshold for the CP regression algorithm.
#' @param max_iter The maximum number of iterations.
#' @param seed The random seed for reproducibility.
#'
#' @return A list containing the estimated tensor B, factor matrices, 
#' convergence status, and the number of iterations.
#'
#' @examples
#' Y <- # input tensor
#' X <- # covariate tensor
#' R <- 3
#' obs_dim_X <- # observed dimensions of X
#' obs_dim_Y <- # observed dimensions of Y
#' cp_regression(Y, X, R, obs_dim_X, obs_dim_Y, convThresh = 1e-05,
#'               max_iter = 500, seed = 0)
#'
#' @export
#' @seealso
#' \code{\link{y_regression}}, \code{\link{x_regression}} 
#' \code{\link{conv_cond}}
#' 
#' @export
cp_regression <- function(Y, X, R, obs_dim_X, obs_dim_Y, convThresh = 1e-05,
                          max_iter = 500, seed = 0) {
  if (seed > 0) set.seed(seed)
  
  # Generate initial random tensor and random CP decomposition
  init_B <- rand_tensor(c(X@modes[-obs_dim_X], Y@modes[-obs_dim_Y]))
  init_list <- list(
    matrix(rnorm(X@modes[2] * R), nrow = X@modes[2]),
    matrix(rnorm(X@modes[3] * R), nrow = X@modes[3]),
    matrix(rnorm(Y@modes[2] * R), nrow = Y@modes[2]),
    matrix(rnorm(Y@modes[3] * R), nrow = Y@modes[3])
  )
  
  converged <- FALSE
  num_iter <- 0
  while (num_iter < max_iter) {
    num_iter <- num_iter + 1
    Ddims <- c(Y@modes[1] * Y@modes[3], Y@modes[1] * Y@modes[2])
    
    for (idx in 1:init_B@num_modes) {
      if (idx < (init_B@num_modes/2 + 1)) {
        x_reg_list <- x_regression(init_list = init_list, Y = Y, X = X,
                                   R = R, idx = idx)
        pre_init_list <- init_list
        init_list[[idx]] <- x_reg_list$B1
        lambdas <- x_reg_list$x_lambdas
        if (convergence_func(pre_init_list, init_list, idx, convThresh)) {
          converged <- TRUE
          break
        }
      } else {
        y_reg_list <- y_regression(init_list = init_list, Y = Y, X = X,
                                         R = R, Ddims= Ddims, idx = idx)
        pre_init_list <- init_list
        init_list[[idx]] <- y_reg_list$B3
        lambdas <- y_reg_list$y_lambdas
        if (convergence_func(pre_init_list, init_list, idx, convThresh)) {
          converged <- TRUE
          break
        }
      }
    }
    if (converged) {
      break  # Exit the loop if converged
    }
  }
  
  return(list(B = init_B, factor_mat = init_list, converged = converged,
              num_iter = num_iter, lambdas = lambdas))
}

# Tucker Regression
xt_regression <- function (X, Y, init_list, idx) {
  omitted_tensor <- ttm(ttm(ttm(init_list[[1]], init_list[[(4-idx)]], (3-idx)), 
                            init_list[[4]], 3), init_list[[5]], 4)
  H <- ttt(X, omitted_tensor, alongA = (4-idx), alongB = (3-idx))
  unfolded_H <- unfold(H, row_idx = c(2,3), col_idx = c(1,4,5))@data
  vec_U <- solve(unfolded_H %*% t(unfolded_H)) %*% unfolded_H %*% vec(Y)
  return(matrix(vec_U, ncol = omitted_tensor@modes[idx]))
}

yt_regression <- function (X, Y, init_list, idx) {
  # This alternates between solving for U3 then solving for U4
  # In order for this to work, we then first omit U3 (include U4) then omit
  # U4 (include U3). Additionally, first multiply over 4 then 3
  # NOTE: init_list includes the core tensor, so shift idx by one
  omitted_tensor <- ttm(ttm(ttm(init_list[[1]], init_list[[2]], 1),
                            init_list[[3]], 2), init_list[[(8-idx)]], (7-idx))
  O <- ttt(X, omitted_tensor, alongA = c(2,3), alongB = c(1,2))
  unfold_dim <- idx-1
  unfolded_O <- unfold(O, row_idx = unfold_dim,
                       col_idx = setdiff(1:3, unfold_dim))@data
  V <- solve(unfolded_O %*% t(unfolded_O)) %*% unfolded_O %*%
    t(unfold(Y, row_idx = unfold_dim, col_idx = setdiff(1:3, unfold_dim))@data)
  return(t(V))
}

core_regression <- function(X, Y, init_list) {
  Xstar <- ttm(ttm(X, t(init_list[[2]]), 2), t(init_list[[3]]), 3)
  Ystar <- ttm(ttm(Y, MASS::ginv(init_list[[4]]), 2),
               MASS::ginv(init_list[[5]]), 3)
  unfolded_Xstar <- unfold(Xstar, row_idx = 1, col_idx = c(2,3))@data
  unfolded_Ystar <- unfold(Ystar, row_idx = 1, col_idx = c(2,3))@data
  G <- solve(crossprod(unfolded_Xstar)) %*% 
    crossprod(unfolded_Xstar, unfolded_Ystar)
  return(G)
}

tuck_conv <- function(pre_init_list, init_list, idx, convThresh) {
  fnorm_conv <- norm(pre_init_list[[idx+1]] - init_list[[idx+1]])
  if(fnorm_conv < convThresh) return(TRUE)
  return(FALSE)
}

#' Perform Tucker regression
#'
#' This function performs Tucker regression on a given tensor, decomposing it
#' into a low-rank tensor and factor matrices.
#'
#' @param Y The target tensor.
#' @param X The predictor tensor.
#' @param R A vector specifying the desired rank of each mode in the Tucker
#'   decomposition.
#' @param obs_dim_X A vector specifying the observed dimensions of the predictor
#'   tensor.
#' @param obs_dim_Y A vector specifying the observed dimensions of the target
#'   tensor.
#' @param convThresh The convergence threshold for the Tucker decomposition.
#' @param max_iter The maximum number of iterations for convergence.
#' @param seed An optional seed for reproducible results.
#'
#' @return A list containing the decomposed components and the rebuilt tensor.
#'
#' @export
tucker_regression <- function(Y, X, R, obs_dim_X, obs_dim_Y, convThresh = 1e-05, 
                              max_iter = 400, seed = 0) {
  if (seed > 0) set.seed(seed)
  
  # Generate initial random tensor and random CP decomposition
  init_B <- rand_tensor(c(X@modes[-obs_dim_X], Y@modes[-obs_dim_Y]))
  # This is a list of a core tensor plus factor matrices
  init_list <- list(
      rand_tensor(R),
      matrix(rnorm(X@modes[2] * R[1]), nrow = X@modes[2]),
      matrix(rnorm(X@modes[3] * R[2]), nrow = X@modes[3]),
      matrix(rnorm(Y@modes[2] * R[3]), nrow = Y@modes[2]),
      matrix(rnorm(Y@modes[3] * R[4]), nrow = Y@modes[3])
  )
  
  converged <- FALSE
  num_iter <- 0
  while (num_iter < max_iter) {
    num_iter <- num_iter + 1
    
    for (idx in 1:init_B@num_modes) {
      if (idx < (init_B@num_modes/2 + 1)) {
        x_results <- xt_regression(init_list = init_list, Y = Y, X = X,
                                   idx = idx)
        pre_init_list <- init_list
        init_list[[idx+1]] <- x_results
        if (tuck_conv(pre_init_list, init_list, idx, convThresh)) {
          converged <- TRUE
          break
        }
      } else {
        y_results <- yt_regression(init_list = init_list, Y = Y, X = X,
                                         idx = idx)
        pre_init_list <- init_list
        init_list[[idx+1]] <- y_results
        if (tuck_conv(pre_init_list, init_list, idx, convThresh)) {
          converged <- TRUE
          break
        }
      }
    }
    flattened_core_update <- core_regression(X = X, Y = Y, init_list = init_list)
    pre_init_list <- init_list
    core_update <- as.tensor(array(flattened_core_update, dim = R))
    init_list[[1]] <- core_update
    
    fnorm_core <- fnorm(pre_init_list[[1]] - init_list[[1]])
    if (fnorm_core < convThresh) converged <- TRUE
    
    if (converged) break  # Exit the loop if converged
  }
  
  return(list(components=init_list, rebuilt_tnsr=tucker_rebuild(init_list)))
}
